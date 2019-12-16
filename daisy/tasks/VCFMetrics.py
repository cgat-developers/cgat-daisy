import os
import re
import json
import pandas
import collections

try:
    import pysam
except ImportError:
    pass

from .Runner import compare_sets, dict_to_labelvalue, is_true, is_set, \
    resolve_argument, get_associated_file

from .MetricRunner import MetricRunner
from daisy.toolkit import update_namedtuple, redirect2mounts
from daisy.tasks import redirect_defaults2mountpoint
from daisy.tasks.Utils import match_sequence_dictionaries

import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools


def build_reference_fasta_map(reference_fasta_map_param):
    """build a dictionary mapping assembly keys to file paths with reference sequence.

    This method takes a hierachical or no-hierarchical dictionary and
    builds reduces it to a non-hierarchical dictionary.
    """
    if reference_fasta_map_param is None:
        return reference_fasta_map_param

    if not isinstance(reference_fasta_map_param, collections.Mapping):
        raise ValueError("expected a dictionary")

    if IOTools.is_nested(reference_fasta_map_param):
        return dict([(x, y["path"]) for x, y in reference_fasta_map_param.items()])
    else:
        return reference_fasta_map_param


def get_reference_for_vcf(infile, fastafiles):
    """deduce reference sequence used within VCF files.

    Note that the VCF file must be indexed with tabix.

    The method works by attempting to read contig names and lengths
    from the VCF file (label: ``contig``). It will then compare the
    sequence dictionary in the vcf file with a list of fastafiles. The
    comparison will stop at the first match that is found.

    The method will also attempt to find entries for the reference
    genome in the VCF (label: ``reference`` or ``assembly``).

    If no lengths are defined, the method will use chromosome names
    only by attempting to fetch variants according to the chromosome
    lists in the VCF files.

    :param infile: :term:`VCF` formatted file
    :param fastafiles: list of :term:`fasta` formatted files. The
        fasta files need to indexed with samtools faidx.
    :param build_map: dictionary mapping builds to references

    :return: a named tuple with fields `filename`, `diffs`,
        `assembly`, `samples`. The first is the filename if
        found, otherwise None.  If not found, diffs is a list of all
        input files with a list missing contigs or length mismatches.
        If given, reference contains the reference cleaned from the
        VCF.

    """

    assembly = None
    fastafn = None

    # Temporary fix: see issue SYS-517
    if not os.path.exists(infile):
        E.warn("could not find file {}".format(bamfile))

    with pysam.VariantFile(infile) as inf:
        sequence_records = []
        for record in inf.header.records:
            if record.key == "contig":
                sequence_records.append((record["ID"],
                                         int(record.get("length", 0))))
            else:
                s = str(record)
                if s.startswith("##reference="):
                    assembly = re.match("##reference=(.*)", s).groups()[0]
                if s.startswith("##assembly="):
                    assembly = re.match("##assembly=(.*)", s).groups()[0]
                samples = inf.header.samples

        vcfdict = dict(sequence_records)

    fastafn, diffs = match_sequence_dictionaries(
        vcfdict, fastafiles)

    if fastafn is None:
        E.warn(
            "no reference match could be found for {}".format(infile))

    result_type = collections.namedtuple(
        "VCFInfo",
        ["reference_fasta", "diffs", "assembly", "samples"])

    return result_type(fastafn, diffs, assembly, samples)


class Preprocessor():

    def build_statement_with_preprocessing(self,
                                           infile, outfile, params,
                                           main_statement):

        tmpfile, pre_statement, post_statement = self.pre_process(infile,
                                                                  outfile,
                                                                  params)
        if main_statement:
            main_statement = main_statement.replace(infile, tmpfile)
        else:
            # if no main statement is given or it is empty, for
            # example if Preprocessor is part of an import task,
            # simply use last result
            pre_statement = pre_statement.replace(tmpfile, outfile)
            main_statement = ""

        stmts = [pre_statement, main_statement, post_statement]
        return " && ".join([x for x in stmts if x])

    def run_with_preprocessing(self, infile, outfile, params, main_statement,
                               **kwargs):
        return P.run(
            self.build_statement_with_preprocessing(
                infile, outfile, params, main_statement),
            **kwargs)


class VCFPreprocessor(Preprocessor):
    """Preprocessor for VCF metrics.

    This class can apply a set of operations on a VCF
    file before it is submitted to a tool or metric.

    add_chr
       add "chr" prefix to chromosome names
    apply_substitution
       apply a perl substitution, such as 's/TC=/DP=/`.
    copy_vcf
       copy vcf file to temporary dir before processing.
    filter_exclude
       apply bcftools filter command with an exclude option
    filter_samples
       apply bcftools filter command and select a set of samples
    filter_vcf
       apply bcftools filter command with an include option
    remove_ambiguity_codes
       remove variants that contain variants with ambiguity codes, i.e.,
       those that are not [ACGT].
    remove_chr
       remove a "chr" prefix from chromosome names.
    remove_multi
       remove loci with overlapping variants
    remove_samples
       remove a set of samples from VCF
    restrict_bed
       restrict VCF file to regions in a bed file
    restrict_region
       restrict VCF file to region string
    restrict_canonical
       restrict VCF file to canonical chromosomes. A canonical chromosome is
       a single number, X, Y, M or MT.
    unset_filter
       set all values in the FILTER column to PASS
    vcfallelicprimitives
       apply vcflib's 'vcfallelicprimitives' to VCF file.
    vcf2vcf
       apply daisy's vcf2vcf tool to VCF file.

    Best used as a Mixin class to a VCF tool/metric.
    """

    add_chr = False
    apply_substitution = None
    copy_vcf = False
    filter_exclude = None
    filter_samples = None
    filter_vcf = None
    reformat = False
    remove_ambiguity_codes = False
    remove_chr = False
    remove_multi = False
    rename_samples = False
    map_unknown_genotypes_to_reference = False
    restrict_bed = None
    restrict_canonical = False
    restrict_region = None
    unset_filter = False
    vcfallelicprimitives = False
    vcf2vcf = None

    path_bcftools = "bcftools"
    path_vcfallelicprimities = "vcfallelicprimitives"
    path_vcfcreatemulti = "vcfcreatemulti"

    def pre_process(self, infile, outfile, params):

        statements = []

        if infile.endswith(".vcf.gz"):
            suffix = "vcf.gz"
            compress_type = "z"
        elif infile.endswith(".bcf"):
            suffix = "bcf"
            compress_type = "b"
        else:
            raise ValueError("unknown suffix for file {}".format(infile))

        # indexing creates .tbi files (for rtg vcfeval)
        infile = IOTools.snip(infile, ".{}".format(suffix))
        if params.copy_vcf:
            statements.append(
                "cp @IN@.{suffix} @OUT@.{suffix}; "
                "cp @IN@.{suffix}.* @OUT@.{suffix}.*".format(**locals()))

        if params.apply_substitution:
            statements.append(
                "{params.path_bcftools} view @IN@.{suffix} "
                "| perl -p -e '{params.apply_substitution}' "
                "2> {outfile}_apply_substitution.log "
                "| {params.path_bcftools} view --output-type {compress_type} "
                "> @OUT@.{suffix}".format(**locals()))

        if params.filter_samples:
            statements.append(
                "{params.path_bcftools} convert "
                "--samples {params.filter_samples} "
                "@IN@.{suffix} "
                "-O z "
                "2> {outfile}_filter_samples.log "
                "> @OUT@.{suffix}".format(**locals()))

        if params.rename_samples:
            with IOTools.open_file(outfile + ".sample_map.tsv", "w") as outf:
                maps = params.rename_samples.split(",")
                for m in maps:
                    if "=" in m:
                        old, new = [x.strip() for x in m.split("=")]
                        outf.write("\t".join((old, new)) + "\n")
                    else:
                        outf.write(m + "\n")

            # reheader outputs compressed vcf
            statements.append(
                "{params.path_bcftools} reheader "
                "--samples {outfile}.sample_map.tsv "
                "@IN@.{suffix} "
                "2> {outfile}_rename_samples.log "
                "> @OUT@.{suffix}".format(**locals()))

        if params.remove_multi:
            statements.append(
                "{params.path_bcftools} view @IN@.{suffix} "
                "| {params.path_vcfcreatemulti} "
                "2> {outfile}_remove_multi.log "
                "| grep -v \"combined=\" "
                "| bgzip "
                "> @OUT@.{suffix}".format(**locals()))

        if params.map_unknown_genotypes_to_reference:
            statements.append(
                "{params.path_bcftools} view @IN@.{suffix} "
                "| perl -p -e \"s/[.][|]1/0|1/g; s/1[|][.]/1|0/g\" "
                "2> {outfile}_map_unknown_genotypes.log "
                "| bgzip "
                "> @OUT@.{suffix}".format(**locals()))

        if params.filter_vcf:
            filter_statement = re.sub("%", "%%", params.filter_vcf)
            statements.append(
                "{params.path_bcftools} filter "
                "--include \"{filter_statement}\" "
                "@IN@.{suffix} "
                "-O z "
                "2> {outfile}_filter.log "
                "> @OUT@.{suffix}".format(**locals()))

        if is_true(params.unset_filter):
            statements.append(
                "{params.path_bcftools} view @IN@.{suffix} "
                "| awk -v OFS=\"\\t\" '!/^#/ {{$7=\"PASS\"}} {{print}}' "
                "| bgzip "
                "2> {outfile}_unset_filter.log "
                "> @OUT@.{suffix}".format(**locals()))

        if params.filter_exclude:
            filter_statement = re.sub("%", "%%", params.filter_exclude)
            statements.append(
                "{params.path_bcftools} filter "
                "--exclude \"{filter_statement}\" "
                "@IN@.{suffix} "
                "-O z "
                "2> {outfile}_filter_exclude.log "
                "> @OUT@.{suffix}".format(**locals()))

        if is_true(params.add_chr):
            statements.append(
                "{params.path_bcftools} view @IN@.{suffix} "
                "| perl -pe '$_ = \"chr\" . $_ unless (/^#/); s/<ID=/<ID=chr/ "
                "if (/^##contig/)'"
                "| {params.path_bcftools} view --output-type {compress_type} "
                "2> {outfile}_add_chr.log "
                "> @OUT@.{suffix}".format(**locals()))

        if is_true(params.remove_chr):
            # also substitute chrM to MT.
            statements.append(
                "{params.path_bcftools} view @IN@.{suffix} "
                "| perl -pe "
                "'s/^chrM/chrMT/ unless /^#/; "
                "s/chrM/chrMT/ if (/^##contig/); "
                "s/^chr// unless (/^#/); "
                "s/<ID=chr/<ID=/ if (/^##contig/)'"
                "| {params.path_bcftools} view --output-type {compress_type} "
                "2> {outfile}_remove_chr.log "
                "> @OUT@.{suffix}".format(**locals()))

        if is_true(params.vcf2vcf):
            statements.append(
                "daisy vcf2vcf "
                "{params.vcf2vcf} "
                "--log {outfile}_vcf2vcf.log "
                "@IN@.{suffix} "
                "| vcf-sort | {params.path_bcftools} view --output-type {compress_type} "
                ">& @OUT@.{suffix}".format(**locals()))

        if is_true(params.restrict_region):
            statements.append(
                "if [[ ! -e @IN@.{suffix}.csi || ! -e @IN@.{suffix}.tbi ]]; "
                "    then {params.path_bcftools} index --tbi -f @IN@.{suffix}; fi; "
                "{params.path_bcftools} view "
                "--regions {params.restrict_region} "
                "@IN@.{suffix} "
                "-O v "
                "2> {outfile}_filter_restrict_region.log "
                "| {params.path_bcftools} view --output-type {compress_type} "
                "> @OUT@.{suffix}".format(**locals()))

        if is_true(params.restrict_canonical):
            statements.append(
                "{params.path_bcftools} view @IN@.{suffix} "
                "| awk '/^#/ {{print; next}} "
                "/^([0-9]+)\\t/ || /^(chr[0-9]+)\\t/ {{print; next}} "
                "/^([XYM])\\t/ || /^(chr[XYM])\\t/ {{print; next}} "
                "$1 == \"MT\" || $1 == \"chrMT\" {{print; next}} ' "
                "2> {outfile}_filter_restrict_canonical.log "
                "| {params.path_bcftools} view --output-type {compress_type} "
                "> @OUT@.{suffix}".format(**locals()))

        if params.restrict_bed:
            statements.append(
                "if [[ ! -e @IN@.{suffix}.csi || ! -e @IN@.{suffix}.tbi ]]; "
                "    then {params.path_bcftools} index --tbi -f @IN@.{suffix}; fi; "
                "{params.path_bcftools} filter "
                "--regions-file {params.restrict_bed} "
                "@IN@.{suffix} "
                "-O v "
                "2> {outfile}_filter_restrict_bed.log "
                "| awk 'BEGIN {{cmd=sprintf(\"sort -k1,1 -k2,2n \"); }} "
                "/^#/ {{ print; next; }} {{ print | cmd; next }} "
                "END {{ close(cmd); }} '"
                "| {params.path_bcftools} view --output-type {compress_type} "
                "> @OUT@.{suffix}".format(**locals()))

        if is_true(params.vcfallelicprimitives):
            # do not use vcfsort from vcflib, it does not work
            # from stream
            statements.append(
                "{params.path_bcftools} view @IN@.{suffix} "
                "| vcfallelicprimitives --keep-geno --keep-info "
                "| vcfstreamsort "
                "| {params.path_bcftools} view --output-type {compress_type} "
                "> @OUT@.{suffix}")

        if is_true(params.remove_ambiguity_codes):
            statements.append(
                "{params.path_bcftools} filter -e \"REF ~ '[^ACGT]'\" "
                "@IN@.{suffix} "
                "-O z "
                "> @OUT@.{suffix}".format(**locals()))

        if not statements:
            return infile + "." + suffix, "", ""

        filename, build_statement, cleanup_statement = P.join_statements(
            statements, infile)
        filename += ".{suffix}".format(**locals())
        build_statement += (
            "; {params.path_bcftools} index --tbi -f {filename} >& {outfile}.index.log"
            .format(**locals()))

        return filename, build_statement, cleanup_statement


class MetricRunnerVCF(MetricRunner, VCFPreprocessor):
    pass


class MetricRunnerVCFTools(MetricRunnerVCF):
    name = "vcftools"
    path = "vcftools"
    method = None

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stdout=True).strip()
        return re.search(r"VCFtools (\S+)", help_string).groups()[0]

    def run(self, infile, outfile, params):

        assert self.method is not None
        retval = self.run_with_preprocessing(
            infile, outfile, params,
            "{params.path} "
            "--gzvcf {infile} "
            "{params.options} "
            "--{params.method} "
            "--stdout "
            "> {outfile}"
            .format(**locals()))

        return retval


class run_metric_vcftools_tstv_summary(MetricRunnerVCFTools):
    """run vcftools TsTv-summary on a :term:`vcf` file.

    This metric permits VCF preprocessing before running
    the analysis.

    *Columns*

    model
        Substitutions, such as AG, AC, AT, GC, GT, Ts, Tv
    count
        Number of substitutions.
    """

    name = "vcftools_tstv_summary"
    method = "TsTv-summary"


class run_metric_vcftools_tstv_by_count(MetricRunnerVCFTools):
    """run vcftools TsTv-by-count on a :term:`vcf` file.

    This metric permits VCF preprocessing before running
    the analysis.

    *Columns*

    alt_allele_count
        Number of alternative alleles.
    n_ts
        Number of transitions at sites with alt_allele_count.
    n_tv
        Number of transversions at sites with alt_allele_count.
    ts/tv
        Transition/transversion ratio at sites with alt_allele_count.
    """

    name = "vcftools_tstv_by_count"
    method = "TsTv-by-count"


class MetricRunnerBCFTools(MetricRunnerVCF):
    path = "bcftools"

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search(r"Version: (\S+)", help_string).groups()[0]

    def split_output(self, outfile):

        def _split_output(lines):
            is_comment = True
            body, header = [], []
            for line in lines:
                if line.startswith("#"):
                    if not is_comment:
                        yield header, body
                        header, body = [], []
                    header.append(line)
                    is_comment = True
                else:
                    fields = line[:-1].split("\t")
                    # sanitize summary fields
                    fields[2] = re.sub(":", "", re.sub("[ -]", "_", fields[2]))
                    body.append(fields)
                    is_comment = False
            yield header, body

        # split into separate files for upload
        with open(outfile) as inf:
            for header, body in _split_output(inf):
                header = re.sub(r"\[\d+\]", "", header[-1][2:])
                fields = header.split("\t")
                identifier = fields.pop(0)
                try:
                    tablename = self.map_section_to_table[identifier]
                except KeyError:
                    continue

                output_file = self.map_table_to_file(tablename, outfile)

                with open(output_file, "w") as outf:
                    outf.write("\t".join(fields))
                    # remove first column, which contains the identifier
                    outf.write("\n".join(
                        ["\t".join(x[1:]) for x in body]) + "\n")


class run_metric_bcftools_stats(MetricRunnerBCFTools):
    """run bcftools stats on a :term:`vcf` file.

    :param infile: Input file in :term:`vcf` format
    :param outfile: Output file in :term:`tsv` format

    This metric permits VCF preprocessing before running
    the analysis.

    This metric outputs several tables.

    allele_frequency
       id
       allele frequency
       number of SNPs
       number of transitions
       number of transversions number of indels
       repeat-consistent
       repeat-inconsistent
       not applicable

    depth_distribution
       id
       bin
       number of genotypes
       fraction of genotypes _percent_
       number of sites fraction of sites _percent_

    indel_context_length
       id
       length of repeat element
       repeat-consistent deletions
       repeat-inconsistent deletions
       consistent insertions
       inconsistent insertions c/_c+i_ ratio

    indel_context_summary
       id
       repeat-consistent
       repeat-inconsistent
       not applicable
       c/_c+i_ ratio

    indel_distribution
       id
       length _deletions negative_
       count

    quality
       id
       Quality number of SNPs
       number of transitions _1st ALT_ number of transversions _1st ALT_
       number of indels

    singleton_stats
       id
       allele count
       number of SNPs
       number of transitions
       number of transversions number of indels
       repeat-consistent
       repeat-inconsistent
       not applicable

    substitution_types
       id
       type
       count

    summary_numbers
       id
       key
       value

    indel_frameshifts
       id
       in_frame
       out_frame
       not_applicable
       out__in_out__ratio
       in_frame_1st_ALT_
       out_frame_1st_ALT
       not_applicable__1st_ALT_
       out__in_out__ratio__1st ALT_

    per_sample_counts
       id
       sample
       nRefHom
       nNonRefHom
       nHets
       nTransitions
       nTransversions
       nIndels
       average_depth
       nSingletons

    per_sample_indels
       id
       sample
       in_frame
       out_frame
       not_applicable
       out__in_out__ratio
       nHets
       nAA

    tstv
       id
       ts
       tv
       ts/tv
       ts _1st ALT_
       tv _1st ALT_
       ts/tv _1st ALT_
    """

    name = "bcftools_stats"
    method = "stats"

    map_section_to_table = {
        'SN': "bcftools_stats_summary_numbers",
        'TSTV': "bcftools_stats_tstv",
        'FS': "indel_frameshifts",
        'ICS': "bcftools_stats_indel_context_summary",
        'ICL': "bcftools_stats_indel_context_length",
        'SiS': "bcftools_stats_singleton_stats",
        'AF': "bcftools_stats_allele_frequency",
        'QUAL': "bcftools_stats_quality",
        'IDD': "bcftools_stats_indel_distribution",
        'ST': "bcftools_stats_substitution_types",
        'DP': "bcftools_stats_depth_distribution",
        'PSC': "bcftools_stats_per_sample_counts",
        'PSI': "bcftools_stats_per_sample_indels",
        'HWE': "bcftools_stats_hw_equilibrium",
    }

    bcftools_restrict_to_region = None

    # if given, used to set --fasta-ref
    reference_fasta = None
    # if given, used to set --exons
    exons = None
    # if given, used to set --regions-file
    regions_file = None
    # if given, used to set --samples option. If set to "auto"
    # all samples will be analyzed - if present.
    samples = None

    def __call__(self, *args, **kwargs):
        self.tablenames = list(self.map_section_to_table.values())
        MetricRunner.__call__(self, *args, **kwargs)

    def run(self, infile, outfile, params):

        if IOTools.is_empty(infile):
            E.warn("file {} is empty - no processing".format(infile))
            IOTools.touch_file(outfile)
            return

        infile_tbi = infile + ".tbi"
        if not os.path.exists(infile_tbi):
            E.warn("VCF file {} is not tabixed - no processing".format(infile))
            IOTools.touch_file(outfile)
            return

        options = params.options
        retvals = []
        reference_label = None

        if params.reference_fasta == "auto" or params.samples == "auto":
            fasta = resolve_argument(list(params.reference_fasta_map.values()),
                                     ",").split(",")
            map_path2name = dict([(x[1], x[0]) for x in list(params.reference_fasta_map.items())])
            vcf_info = get_reference_for_vcf(infile, fasta)
            if not vcf_info:
                E.warn("attempted to detect reference fasta, but unable to do so. "
                       "diffs: {}".format(vcf_info.diffs))

        if params.reference_fasta:
            if params.reference_fasta == "auto":
                if vcf_info:
                    options += " --fasta-ref {}".format(vcf_info.reference_fasta)
                    reference_label = map_path2name[vcf_info.reference_fasta]
            else:
                options += " --fasta-ref {}".format(params.reference_fasta)
                reference_label = map_path2name.get(params.reference_fasta, None)

        if params.regions_file:
            regions_file = get_associated_file(params,
                                               reference_label,
                                               "regions_file")
            options += " --regions-file {}".format(regions_file)

        if params.exons:
            exons = get_associated_file(params,
                                        reference_label,
                                        "exons")
            options += " --exons {}".format(exons)

        if params.samples:
            if params.samples == "auto":
                if len(vcf_info.samples) > 0:
                    options += " --samples -"
            else:
                options += "--samples {}".format(params.samples)

        if is_set(params.bcftools_restrict_to_region):
            region = params.bcftools_restrict_to_region.strip()
            # This convoluted processing is necessary as giving both
            # -R and -r at the same time to bcftools does not work and
            # I want to use the same syntax with or without a
            # region-file.  VCFMetrics.restrict_region implements a
            # generic solution by subsetting the VCF, but for large
            # data sets this is a waste of effort.

            # modify options string if a regions-file is given and it
            # should be restricted to a particular genomic
            # region. Note that the region-file needs to be indexed by
            # tabix.
            if "--regions-file" in options or "-r " in options:
                try:
                    orig_bedfile = re.search(r"(-r |--regions-file=*)\s*(\S+)", options).groups()[1]
                except AttributeError:
                    raise ValueError("could not extract regions file from '{}'".format(options))
                new_bedfile = outfile + ".regions.gz"
                retvals.extend(P.run("tabix {} {} > {}".format(
                    orig_bedfile,
                    region,
                    new_bedfile)))
                options = re.sub(orig_bedfile, new_bedfile, options)
            else:
                # otherwise, simply add to option string.
                options += " --regions={}".format(region)

        retvals.extend(self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{params.path} stats "
            "{options} "
            "{infile} "
            "2> {outfile}.log "
            "> {outfile}".format(**locals())))

        self.split_output(outfile)

        return retvals


class run_metric_bcftools_guess_ploidy(MetricRunnerBCFTools):
    """run bcftools guess-ploidy on a :term:`vcf` file.

    :param infile: Input file in :term:`vcf` format
    :param outfile: Output file in :term:`tsv` format

    """

    name = "bcftools_vcf2sex"
    plugin_options = ""

    def run(self, infile, outfile, params):

        with IOTools.open_file(outfile, "w") as outf:
            outf.write("sample\tsex\n")

        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{params.path} plugin "
            "guess-ploidy "
            "{params.options} "
            "{infile} "
            "-- "
            "{params.plugin_options} "
            "2> {outfile}.log "
            ">> {outfile}".format(**locals()))
        return retval


class run_metric_bcftools_vcf2sex(run_metric_bcftools_guess_ploidy):
    """run bcftools vcf2sex on a :term:`vcf` file (deprecated).

    :param infile: Input file in :term:`vcf` format
    :param outfile: Output file in :term:`tsv` format

    Note that vcf2sex has been deprecated as of bcftools 1.4 in favour
    of guess-ploidy. This metric is present for backwards
    compatibility, use bcftools_guess_ploidy instead.

    """

    name = "bcftools_vcf2sex"


class run_metric_bcftools_gtcheck(MetricRunnerBCFTools):
    """run bcftools gtcheck on a :term:`vcf` file.

    :param infile: Input file in :term:`vcf` format
    :param outfile: Output file in :term:`tsv` format

    """

    name = "bcftools_gtcheck"

    map_section_to_table = {
        'SM': "bcftools_gtcheck_sample",
        'CN': "bcftools_gtcheck_concordance",
    }

    def __call__(self, *args, **kwargs):
        self.tablenames = list(self.map_section_to_table.values())
        MetricRunner.__call__(self, *args, **kwargs)

    def run(self, infile, outfile, params):

        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{params.path} gtcheck "
            "{params.options} "
            "{infile} "
            "2> {outfile}.log "
            "> {outfile}.tmp".format(**locals()))

        P.run("mv {}.tmp {}".format(outfile, outfile))

        self.split_output(outfile)

        return retval


class run_metric_bcftools_roh(MetricRunnerBCFTools):
    """run bcftools roh on a :term:`vcf` file.

    :param infile: Input file in :term:`vcf` format
    :param outfile: Output file in :term:`tsv` format
    """

    name = "bcftools_roh"
    reference_fasta = None

    map_section_to_table = {
        'RG': "bcftools_roh_regions",
        'ST': "bcftools_roh_states",
    }

    def __call__(self, *args, **kwargs):
        self.tablenames = list(self.map_section_to_table.values())
        MetricRunner.__call__(self, *args, **kwargs)

    def split_output(self, outfile):
        # roh output is interspersed

        outfiles = {}

        with open(outfile) as inf:
            is_comment = True
            for line in inf:
                if line.startswith("#"):
                    # comment section ends at first empty comment line
                    if line[1:-1].strip() == "":
                        is_comment = False
                        continue
                    if is_comment:
                        continue
                    header = re.sub(r"\[\d+\]", "", line[2:-1]).split("\t")
                    identifier = header.pop(0)
                    try:
                        tablename = self.map_section_to_table[identifier]
                    except KeyError:
                        continue
                    outfiles[identifier] = IOTools.open_file(
                        self.map_table_to_file(tablename, outfile), "w")
                    outfiles[identifier].write("\t".join(header) + "\n")
                else:
                    data = line[:-1].split("\t")
                    identifier = data.pop(0)
                    outfiles[identifier].write("\t".join(data) + "\n")

    def run(self, infile, outfile, params):

        if params.reference_fasta is None:
            raise ValueError("{} requires reference_fasta".format(self.name))

        with pysam.VariantFile(infile) as vcf:
            samples = [x for x in vcf.header.samples if x != "somatic"]

        samples = ",".join(samples)
        retvals = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{params.path} roh "
            "--samples {samples} "
            "{params.options} "
            "{infile} "
            "2> {outfile}.err "
            "> {outfile}".format(**locals()))

        self.split_output(outfile)

        return retvals


def restrict_bed(outfile, bedfile, vcffile,
                 remove_chr=False, add_chr=False):
    """restrict analysis to the chromosomes in the VCF.

    This is important in order to get accurated false-negative counts
    in a per-chromosome analysis.

    Note that this only works for a per-chromosome analysis - if the
    VCF is only a part of a chromosome the counts will still be an
    over-estimate.

    """
    # This method is a patch, should be better placed in a mixin.
    if IOTools.is_empty(vcffile):
        IOTools.touch(outfile)
        return outfile

    if not os.path.exists(bedfile):
        raise OSError("file {} does not exists".format(bedfile))

    tbx = pysam.TabixFile(vcffile)

    contigs = tbx.contigs
    tbx.close()

    if remove_chr:
        contigs = [re.sub("chr", "", x) for x in contigs]
    if add_chr:
        contigs = ["chr" + x for x in contigs]

    if len(contigs) == 0:
        raise ValueError("no contigs to filter for")

    contig_filter = " || ".join(
        ['$1 == "{}"'.format(x) for x in contigs])

    P.run("zcat {bedfile} "
          "| awk '{contig_filter}' "
          "| bgzip > {outfile}.tmp; ".format(**locals()))
    P.run("mv {}.tmp {}".format(outfile, outfile))
    P.run("tabix -p bed {outfile} ".format(**locals()))

    return outfile


def create_genome_bed(outfile, vcffile, fastafile,
                      remove_chr=False, add_chr=False,
                      restrict_canonical=True):
    """

    restrict to canonical chromosomes.
    """

    if not IOTools.is_empty(vcffile):
        tbx = pysam.TabixFile(vcffile)
        vcf_contigs = tbx.contigs
        if remove_chr:
            vcf_contigs = [re.sub("chr", "", x) for x in vcf_contigs]
        elif add_chr:
            vcf_contigs = ["chr" + x for x in vcf_contigs]
        vcf_contigs = set(vcf_contigs)
        tbx.close()
    else:
        vcf_contigs = None

    if restrict_canonical:
        canonical = set([str(x) for x in range(0, 22)] +
                        ['X', 'Y', 'M', 'MT', 'Mt'])
        if add_chr:
            canonical = set(["chr" + x for x in canonical])
    else:
        canonical = None

    fn = IOTools.snip(outfile, ".gz")
    with pysam.FastaFile(fastafile) as inf, IOTools.open_file(fn, "w") as outf:
        for contig, size in zip(inf.references,
                                inf.lengths):
            if vcf_contigs and contig not in vcf_contigs:
                continue
            if canonical and contig not in canonical:
                continue
            outf.write("{}\t{}\t{}\n".format(
                contig, 0, size))

    pysam.tabix_index(fn, preset="bed", force=True)
    return outfile


class run_metric_bcftools_stats_compare_vcf_to_truth_vcf(MetricRunnerBCFTools):
    """run bcftools stats on a :term:`vcf` file
    against a reference :term:`vcf` file.

    This metric permits VCF preprocessing before running
    the analysis.
    """

    name = "bcftools_stats_compare_vcf_to_truth_vcf"

    map_section_to_table = {
        'SN': "bcftools_stats_summary_numbers",
        'TSTV': "bcftools_stats_tstv",
        'FS': "indel_frameshifts",
        'ICS': "bcftools_stats_indel_context_summary",
        'ICL': "bcftools_stats_indel_context_length",
        'SiS': "bcftools_stats_singleton_stats",
        'AF': "bcftools_stats_allele_frequency",
        'QUAL': "bcftools_stats_quality",
        'IDD': "bcftools_stats_indel_distribution",
        'ST': "bcftools_stats_substitution_types",
        'DP': "bcftools_stats_depth_distribution",
        'PSC': "bcftools_stats_per_sample_counts",
        'PSI': "bcftools_stats_per_sample_indels",
        'HWE': "bcftools_stats_hw_equilibrium",
        'NRDs': "bcftools_stats_nrd_snps",
        'NRDi': "bcftools_stats_nrd_indels",
    }

    bcftools_restrict_to_region = None
    reference_fasta = None

    def run(self, infile, outfile, params):

        if params.reference_fasta is None:
            raise ValueError("missing input parameter 'reference_fasta'")

        if IOTools.is_empty(infile):
            E.warn("file {} is empty - no processing".format(infile))
            IOTools.touch_file(outfile)
            return

        options = params.options
        retvals = []
        if is_set(params.bcftools_restrict_to_region):
            region = params.bcftools_restrict_to_region.strip()
            # This convoluted processing is necessary as giving both
            # -R and -r at the same time to bcftools does not work and
            # I want to use the same syntax with or without a
            # region-file.  VCFMetrics.restrict_region implements a
            # generic solutions by subsetting the VCF, but for large
            # data sets this is a waste of effort.

            # modify options string if a regions-file is given and it
            # should be restricted to a particular genomic
            # region. Note that the region-file needs to be indexed by
            # tabix.
            if "--regions-file" in options or "-r " in options:
                try:
                    orig_bedfile = re.search(r"(-r |--regions-file=*)\s*(\S+)", options).groups()[1]
                except AttributeError:
                    raise ValueError("could not extract regions file from '{}'".format(options))
                new_bedfile = outfile + ".regions.gz"
                retvals.extend(P.run("tabix {} {} > {}".format(
                    orig_bedfile,
                    region,
                    new_bedfile)))
                options = re.sub(orig_bedfile, new_bedfile, options)
            else:
                # otherwise, simply add to option string.
                options += " --regions={}".format(region)

        retvals.extend(self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{params.path} stats "
            "-s - "
            "--fasta-ref {params.reference_fasta} "
            "{options} "
            "{infile} "
            "{params.reference_vcf} "
            "2> {outfile}.log "
            "> {outfile}.tmp".format(**locals())))
        P.run("mv {}.tmp {}".format(outfile, outfile))

#        self.split_output(outfile)

        return retvals


class run_metric_useq_vcfcomparator(MetricRunnerVCF):
    path = "/data/install/Free/useq-8.9.6/Apps/VCFComparator"
    name = "useq_vcfcomparator"

    def get_version(self):
        help_string = E.run("java -jar {self.path} -h".format(**locals()),
                            return_stdout=True).strip()
        return re.search(r"USeq_(\S+)", help_string).groups()[0]

    def run(self, infile, outfile, params):

        outfile_regions = outfile + ".bed.gz"
        restrict_bed(outfile_regions,
                     params.callable_bed,
                     infile,
                     remove_chr=params.remove_chr,
                     add_chr=params.add_chr)

        return self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "java -jar -Xmx4g {params.path} "
            "-a {params.reference_vcf} "
            "-b {outfile_regions} "
            "-c {infile} "
            "-d {outfile_regions} "
            "-p {outfile}.out "
            "{params.options} "
            ">& {outfile}.log".format(**locals()))


class run_metric_bcftools_query(MetricRunnerBCFTools):
    """output counts of filter codes from a :term:`vcf` file.

    :param infile: Input file in :term:`vcf` format
    :param outfile: Output file in :term:`tsv` format

    force_unique : string
        If set, enforce uniqueness of values. Possible values are
        `first`: take first value, `last`: take last value,
        `ignore`: ignore any rows with multiple values.

    add_variant_type: bool
        If set, add the variant type. The variant type is deduced
        by looking at ALT versus REF and the following decision rules:
        Multiple alleles: MUL
        Two single letters: SNP
        Two strings of different length: INS/DEL
        Two strings of the same length: MNP

    This metric permits VCF preprocessing before running
    the analysis.
    """
    name = "bcftools_query"
    fields = "FILTER"

    force_unique = None

    add_variant_type = False

    def run(self, infile, outfile, params):

        fields = self.fields.split(",")
        headers = fields[:]

        variant_type_statement = ""

        if params.add_variant_type:
            index_ref = len(fields) + 1
            index_alt = index_ref + 1
            index_type = index_alt + 1
            fields.extend(["REF", "ALT"])
            headers.extend(["REF", "ALT", "TYPE"])
            variant_type_statement = (
                "| awk '{{t=\"?\"}} "
                " length(${index_alt})==length(${index_ref}) {{t=\"MNP\"}} "
                " length(${index_alt})> length(${index_ref}) {{t=\"INS\"}} "
                " length(${index_alt})< length(${index_ref}) {{t=\"DEL\"}} "
                "(length(${index_alt})==1) && "
                "(length(${index_ref})==1)                   {{t=\"SNP\"}} "
                "${index_alt} ~ /,/                          {{t=\"MUL\"}} "
                "{{print $0 \"\\t\" t}} '".format(**locals())
            )

        filter_option = '-f "[%%{}]\\n"'.format(
            "\\t%%".join([x.strip() for x in fields]))

        with IOTools.open_file(outfile, "w") as outf:
            outf.write("\t".join(headers) + "\n")

        filter_statement = ""
        if params.force_unique:
            if params.force_unique == "ignore":
                filter_statement = '| grep -v ","'
            elif params.force_unique == "first":
                filter_statement = r'| perl -p -e "s/,\S+//g"'
            elif params.force_unique == "last":
                filter_statement = r'| perl -p -e "s/\S+,//g"'
            else:
                raise ValueError("unknown make_unique value {}".format(
                    self.make_unique))

        options = params.options
        if "%" in options and "%%" not in options:
            options = re.sub("%", "%%", options)

        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{params.path} query "
            "{options} "
            "{filter_option} "
            "{infile} "
            "2> {outfile}.log "
            "{variant_type_statement} "
            "{filter_statement} "
            ">> {outfile}.tmp".format(**locals()))
        P.run("mv {}.tmp {}".format(outfile, outfile))

        return retval


class run_metric_bcftools_query_sample(MetricRunnerBCFTools):
    """run bcftools query on a :term:`vcf` file.

    Output a table with particular metrics per sample. The
    resulting table is melted and has the following columns:

    *Columns*
    sample
        sample name
    format
        the metric name, corresponding to a format id within the
        VCF.
    value
        metric value

    :param infile: Input file in :term:`vcf` format
    :param outfile: Output file in compressed :term:`tsv`
         format

    """

    name = "bcftools_query_sample"

    header = None
    known_columns = ["CHROM", "POS", "REF", "ALT", "FILTER"]

    melt = True

    def run(self, infile, outfile, params):

        if params.header:
            with IOTools.open_file(outfile, "w") as outf:
                outf.write(params.header)

        # quote percent chars if not already quoted
        options = params.options
        if "%" in params.options and "%%" not in params.options:
            options = re.sub("%", "%%", options)

        if params.melt:
            statement = (
                "{params.path} query -H "
                "{options} "
                "{infile} "
                "2> {outfile}.log "
                r"| perl -p -e 's/\[\d+\]//g if (/^#/); s/^# *//' "
                "> {outfile}.tmp".format(**locals()))
        else:
            statement = (
                "{params.path} query -H "
                "{options} "
                "{infile} "
                "2> {outfile}.log "
                r"| perl -p -e 's/\[\d+\]//g if (/^#/); s/:\S+//g; s/^# *//' "
                "> {outfile}".format(**locals()))

        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            statement)

        if params.melt:
            table = pandas.read_csv(outfile + ".tmp", sep="\t")
            id_vars = [x for x in self.known_columns if x in table.columns]
            melted = pandas.melt(table, id_vars=id_vars)
            split = pandas.DataFrame.from_records(
                list(melted["variable"].str.split(":")),
                columns=["sample", "format"])
            melted.drop("variable", inplace=True, axis=1)

            df = pandas.concat([melted, split], axis=1)
            df = df[id_vars + ["sample", "format", "value"]]
            df = df[df["sample"] != "Unnamed"]

            df.to_csv(outfile, sep="\t", index=False)
            os.unlink(outfile + ".tmp")

        return retval


class run_metric_daisy_vcf_stats(MetricRunnerVCF):
    """run vcfstats on a :term:`vcf` file.

    :param infile: Input file in :term:`vcf` format
    :param outfile: Multiple output files in compressed :term:`tsv`
         format

    """

    name = "daisy_vcf_stats"
    path = "daisy"

    tablenames = ["daisy_vcf_stats_mutation_profile",
                  "daisy_vcf_stats_mutation_profile_signatures",
                  "daisy_vcf_stats_kinship",
                  "daisy_vcf_stats_format_per_sample",
                  "daisy_vcf_stats_gc_context",
                  "daisy_vcf_stats_gc_dp_prof",
                  "daisy_vcf_stats_gc_dp_profile_stats",
                  "daisy_vcf_stats_format_unset_samples",
                  "daisy_vcf_stats_format_unset_sites"]

    upload_melted = {
        "daisy_vcf_stats_format_unset_samples":
        {"id_vars": ["FORMAT"],
         "var_name": "sample",
         "value_name": "unset"},
        "daisy_vcf_stats_gc_context":
        {"id_vars": ["percent_gc"],
         "var_name": "sample",
         "value_name": "counts"},
    }

    upload_as_matrix = {
        "daisy_vcf_stats_format_per_sample": "FORMAT",
    }

    reference_fasta = None

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        if params.reference_fasta is None:
            raise ValueError("please specify a reference fasta file")

        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{params.path} vcf-stats "
            "--output-filename-pattern={outfile}.daisy_vcf_stats_%%s.tsv "
            "--force-output "
            "--input-fasta={params.reference_fasta} "
            "{params.options} "
            "{infile} "
            ">& {outfile} ".format(**locals()))

        return retval


class run_metric_daisy_vcf_compare_phase(MetricRunnerVCF):
    """run vcf-compare-phase on a :term:`vcf` file against
    a standard of truth.

    :param infile: Input file in :term:`vcf` format
    :param outfile: Multiple output files in compressed :term:`tsv`
         format

    """

    name = "daisy_vcf_compare_phase"
    path = "daisy"

    reference_vcf = None
    reference_fasta = None

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{params.path} vcf-compare-phase "
            "--input-vcf {infile} "
            "--truth-vcf {params.reference_vcf} "
            "--input-fasta {params.reference_fasta} "
            "--log {outfile}.log "
            "{params.options} "
            "{infile} "
            "> {outfile}.tmp ".format(**locals()))
        P.run("mv {}.tmp {}".format(outfile, outfile))

        return retval


class run_metric_whatshap_compare_phase(MetricRunnerVCF):
    """run whatshap comparison for phasing for a phased file against a truth phased vcf

    """

    name = "whatshap_compare_phase"
    path = "whatshap"

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stdout=True).strip()
        return help_string

    def run(self, infile, outfile, params):

        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{params.path} compare "
            "{infile} "
            "{params.reference_vcf} "
            "2> {outfile}.tmp".format(**locals()))
        P.run("mv {}.tmp {}".format(outfile, outfile))

        return retval


class run_metric_daisy_vcfqc_report(MetricRunnerVCF):
    """run vcfqc-report on a :term:`vcf` file.
    """
    name = "daisy_vcfqc_report"
    path = "daisy"

    pedigree_file = None

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        if params.pedigree_file:
            with IOTools.open_file(params.pedigree_file) as inf:
                pedigrees = [x.split("\t") for x in inf]

            for family_id, individual_id, paternal_id, maternal_id, sex, phenotype in pedigrees:
                if family_id in outfile:
                    pedigree_option = "--pedigree={individual_id},{maternal_id},{paternal_id} "
                    "--sex={sex} ".format(**locals())
                    break
            else:
                raise ValueError(
                    "could not find pedigree information for {}".format(outfile))

        else:
            pedigree_option = ''

        outdir = os.path.dirname(outfile)

        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{params.path} vcfqc-report "
            "--output-dir={outdir} "
            "{params.options} "
            "{pedigree_option} "
            "{infile} "
            ">& {outfile} ".format(**locals()),
            job_memory="32G"
        )


class run_metric_variant_effect_predictor(MetricRunnerVCF):
    """run ENSEMBL's variant effect predictor.

    :param infile: Input file in :term:`vcf` format
    :param outfile: Multiple output files in compressed :term:`tsv`
         format

    """
    name = "variant_effect_predictor"
    path = "/data/install/Free/ensembl-tools-release-80/scripts/variant_effect_predictor/variant_effect_predictor.pl"

    dir_ensembl = "/data/library/Ensembl"
    dir_plugin = "/data/install/Free/ensembl-tools-release-80/plugins/"
    reference_fasta = None
    assembly = None

    def get_version(self):
        help_string = E.run("perl {self.path}".format(**locals()),
                            return_stdout=True).strip()
        return re.search(r"version (\S+)", help_string).groups()[0]

    def run(self, infile, outfile, params):

        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "perl {params.path} "
            "--offline "
            "--cache "
            "--dir {params.dir_ensembl} "
            "--assembly {params.assembly} "
            "--fasta {params.reference_fasta} "
            "--dir_plugins {params.dir_plugin} "
            "--force_overwrite "
            "--format vcf "
            "-i {infile} "
            "--vcf "
            "-o STDOUT "
            "2> {outfile}.log "
            "| bgzip "
            "> {outfile}.tmp ".format(**locals()))
        P.run("mv {}.tmp {}".format(outfile, outfile))

        return retval


class run_metric_vcf_identity(MetricRunnerVCFTools):
    """
    Create a copy of output file, applying any preprocessing
    options available to VCF files.
    """

    name = "vcf_identity"

    def run(self, infile, outfile, params):

        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "cp {infile} {outfile}".format(**locals()))

        return retval


class run_metric_daisy_vcf2reference(MetricRunnerVCF):
    """get genome build and GENOMICS reference for a VCF file.

    .. note::

        This tool runs locally.

    :param infile: Input file in :term:`VCF` format
    :param outfile: Output file in :term:`tsv` format

    """

    name = "daisy_vcf2reference"
    path = "daisy"

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        if params.reference_fasta_map is None:
            raise ValueError("vcf2reference requires a reference sequence map")

        infile_tbi = infile + ".tbi"
        if not os.path.exists(infile_tbi):
            E.warn("VCF file {} is not tabixed - no processing".format(infile))
            IOTools.touch_file(outfile)
            return

        fasta = resolve_argument(list(
            params.reference_fasta_map.values()), ",").split(",")

        vcf_info = get_reference_for_vcf(infile, fasta)
        if vcf_info.reference_fasta is None:
            if vcf_info.diffs is None:
                reference_fasta = "corrupted"
            else:
                reference_fasta = "unknown"
                E.debug("differences: {}".format(str(vcf_info.diffs)))
            path = ""
        else:
            reference_fasta = vcf_info.reference_fasta
            map_path2name = dict([(x[1], x[0]) for x in list(params.reference_fasta_map.items())])
            path = map_path2name.get(reference_fasta,
                                     os.path.basename(reference_fasta))

        # todo: add build number?
        with IOTools.open_file(outfile, "w") as outf:
            outf.write("filename\treference\tpath\tassembly\tnsamples\n")
            outf.write("\t".join(map(str,
                                     [infile,
                                      reference_fasta,
                                      path,
                                      vcf_info.assembly if vcf_info.assembly else "",
                                      len(vcf_info.samples)])) + "\n")

        return None


class run_metric_daisy_vcf2tsv(MetricRunnerVCF):
    """get genome build and GENOMICS reference for a VCF file.

    .. note::

        This tool runs locally.

    :param infile: Input file in :term:`VCF` format
    :param outfile: Output file in :term:`tsv` format

    """

    name = "daisy_vcf2tsv"
    path = "daisy"

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        return self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "daisy vcf2tsv "
            "{params.options} "
            "--log {outfile}.log "
            "{infile} "
            "> {outfile}".format(**locals()))


class MetricRunnerVCFBedtools(MetricRunnerVCF):
    path = "bedtools"

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stdout=True).strip()
        return re.search(r"bedtools (\S+)", help_string).groups()[0]


class run_metric_vcf_bedtools_jaccard(MetricRunnerVCFBedtools):
    """run bedtools jaccard on a :term:`vcf` file comparing
    against a rference :term:`vcf` file.

    :param infile: Input file in :term:`vcf` format
    :param outfile: Output file in :term:`tsv` format

    *Columns*

    intersection
        number of nucleotides in the intersection of both interval sets
    union-intersection
        number of nucleotides in the union of both interval sets minus the
        intersection
    jaccard
        jaccard coefficiont, the ratio of intersection and union-intersection.
    n_intersections
        number of intersections
    """

    name = "vcf_bedtools_jaccard"
    reference_vcf = None

    def run(self, infile, outfile, params):

        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{params.path} jaccard "
            "-a {infile} "
            "-b {params.reference_vcf} "
            "{params.options} "
            "2> {outfile}.log "
            "> {outfile}"
            .format(**locals()))

        return retval
