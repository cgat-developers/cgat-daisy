import re
import os
import pysam

from .MetricRunner import MetricRunner
from .Runner import resolve_argument, is_true, get_associated_file
import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools
from .VCFMetrics import Preprocessor, build_reference_fasta_map
from daisy.tasks.Utils import match_sequence_dictionaries


def get_reference_for_bam(bamfile, fastafiles):
    """deduce reference sequence used within BAM files.

    This method compares the sequence dictionary in the bamfile with a
    list of fastafiles. The comparison will stop at the first match
    that is found.

    :param bamfile: :term:`BAM` formatted file
    :param fastafiles: list of :term:`fasta` formatted files. The
        fasta files need to indexed with samtools faidx.
    :return: a tuple (filename, diffs). The first is the filename if
        found, otherwise None.  If not found, diffs is a list of all
        input files with a list missing contigs or length mismatches.

    """

    diffs = []

    # Temporary fix: see issue SYS-517
    if not os.path.exists(bamfile):
        E.warn("could not find file {}".format(bamfile))

    try:
        with pysam.AlignmentFile(bamfile, check_sq=False) as inf:
            sequence_dict = dict(list(zip(inf.references, inf.lengths)))
    except IOError as ex:
        E.warn("could not open bamfile {}: {}".format(bamfile, ex))
        return None, None

    fastafn, diffs = match_sequence_dictionaries(
        sequence_dict, fastafiles)

    return fastafn, diffs


class BAMPreprocessor(Preprocessor):
    """Preprocessor for BAM metrics.

    This class can apply a set of operations on a BAM
    file before it is submitted to a tool or metric.

    remove_chr
       remove a "chr" prefix from chromosome names.

    """

    copy_bam = False
    split_bam = None
    region = None
    bam2bam = None
    remove_chr = None
    shift_quality = None

    def pre_process(self, infile, outfile, params):

        statements = []
        infile = IOTools.snip(infile, ".bam")
        tmpdir = P.get_parameters_as_namedtuple().tmpdir
        outprefix = os.path.basename(os.path.dirname(outfile))

        if params.copy_bam:
            statements.append(
                "cp @IN@.bam @OUT@.bam; "
                "cp @IN@.bam.bai @OUT@.bam.bai")

        if params.split_bam:
            statements.append(
                "daisy bam2bam-split-reads "
                "-i @IN@.bam "
                "-o - "
                "{params.split_bam} "
                "--log={outfile}_split_bam.log "
                "2> {outfile}_split_bam.err "
                "> @OUT@.bam; ".format(**locals()))

        if params.bam2bam:
            statements.append(
                "daisy bam2bam "
                "--stdin=@IN@.bam "
                "{params.bam2bam} "
                "--log={outfile}_bam2bam.log "
                "2> {outfile}_bam2bam.err "
                "> @OUT@.bam; ".format(**locals()))

        if params.region:
            statements.append(
                "samtools view -b @IN@.bam {} > @OUT@.bam".format(params.region))

        if params.shift_quality:
            statements.append(
                "samtools view -h @IN@.bam "
                "| perl -lane "
                "'if(/^@/) {{print; next;}} "
                "@qual=split(//, $F[10]); "
                "$_=chr(ord($_)+{}) for (@qual); "
                "$F[10]=join(\"\",@qual); "
                "print join(\"\\t\", @F)' "
                "| samtools view -bS > @OUT@.bam".format(params.shift_quality))

        if is_true(params.remove_chr):
            # also substitute chrM to MT.
            statements.append(
                "samtools view -h @IN@.bam "
                "| awk -v OFS='\\t' '"
                "$1 == \"@SQ\" "
                "{{ gsub(\"chrM\", \"chrMT\", $2); "
                "   gsub(\"chr\", \"\", $2); print; next }} "
                "{{ gsub(\"chrM\", \"chrMT\", $3); "
                "   gsub(\"chr\", \"\", $3); print; next}} '"
                "| samtools view -bS - "
                "2> {outfile}_remove_chr.log "
                "> @OUT@.bam; ".format(**locals()))

        if not statements:
            return infile + ".bam", "", ""

        filename, build_statement, cleanup_statement = P.join_statements(
            statements, infile)
        filename += ".bam"
        build_statement += (
            "; samtools index {filename} >& {outfile}.index.log"
            .format(**locals()))

        return filename, build_statement, cleanup_statement


class MetricRunnerSamtools(MetricRunner):
    path = "samtools"

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]


class run_metric_samtools_idxstats(MetricRunnerSamtools):
    """run samtools idxstats on a :term:`bam` file.

    :param infile: Input file in :term:`bam` format
    :param outfile: Output file in :term:`tsv` format

    *Columns*

    chromosome
        chromosome name
    size
        chromosome length
    mapped
        number of reads mapped to chromosome
    unmapped
        number of reads unmapped, but assigned to chromosome
    """

    name = "samtools_idxstats"

    def run(self, infile, outfile, params):

        with open(outfile, "w") as outf:
            outf.write("chromosome\tsize\tmapped\tunmapped\n")

        try:
            retval = P.run("{params.path} idxstats "
                           "{infile} "
                           "2> {outfile}.log "
                           ">> {outfile}; "
                           .format(**locals()))
        except OSError as e:
            E.warn("input  file {} gave the following errors: {}".format(
                infile, str(e)))
            retval = None

        return retval


class run_metric_samtools_header(MetricRunnerSamtools):
    """get information in the header of a bam-file.

    :param infile: Input file in :term:`bam` format
    :param outfile: Output file in :term:`tsv` format

    *Columns*

    header_tag
        :term:`bam` header tag such as `SQ`
    tag
        :term:`bam` tag, such as `LN`
    lineno:
        line number in header of tag. Can be used to group tags
        that belong together.
    value
       contents of tag.

    """

    name = "samtools_header"

    def run(self, infile, outfile, params):

        try:
            retval = P.run("{params.path} view -H  "
                           "{infile} "
                           "2> {outfile}.log "
                           "> {outfile}.tmp; "
                           .format(**locals()))
        except OSError as e:
            E.warn("input file {} gave the following errors: {}".format(
                infile, str(e)))

        with open(outfile, "w") as outf, open(outfile + ".tmp") as inf:
            outf.write("header_tag\ttag\tlineno\tvalue\n")
            for lineno, line in enumerate(inf):
                fields = line[1:-1].split("\t")
                header_tag = fields[0]
                if header_tag == "CO":
                    # Do not split comment lines
                    outf.write("\t".join(
                        (header_tag, "", str(lineno),
                         "\t".join(fields[1:]))) + "\n")
                else:
                    for field in fields[1:]:
                        sub_tag, content = field.split(":", 1)
                        outf.write("\t".join(
                            (header_tag, sub_tag,
                             str(lineno), content)) + "\n")

        os.unlink(outfile + ".tmp")
        return retval


class run_metric_samtools_flagstat(MetricRunnerSamtools):
    """run samtools flagstat on a bam-file.

    :param infile: Input file in :term:`bam` format
    :param outfile: Output file in :term:`tsv` format

    *Columns*

    counts
        Number of reads in a category (passing QC)
    counts_fail
        Number reads failing QC in a category
    category
        Category of counts, such as 'mapped', 'properly paired'
    """

    name = "samtools_flagstat"

    def run(self, infile, outfile, params):

        with open(outfile, "w") as outf:
            outf.write("counts\tcounts_fail\tcategory\n")

        try:
            retval = P.run("{params.path} flagstat "
                           "{infile} "
                           "2> {outfile}.log "
                           "| perl -p -e 's/ \+ /\\t/; s/ /\\t/; s/\\(.*//' "
                           ">> {outfile}; "
                           .format(**locals()))
        except OSError as e:
            E.warn("input  file {} gave the following errors: {}".format(
                infile, str(e)))

        return retval


class run_metric_samtools_stats(MetricRunnerSamtools):
    """run samtools flagstat on a bam-file.

    :param infile: Input file in :term:`bam` format
    :param outfile: Output file in :term:`tsv` format

    This metric outputs several tables:

    ACGT_PER_CYCLE
        cycle
            Machine cycle
        percent_A
            Percent bases A
        percent_C
            Percent bases C
        percent_G
            Percent bases G
        percent_T
            Percent bases T

    GC_depth
        percent_GC
            Percent GC
        unique_percentile
        depth_10th_percentile
        depth_25th_percentile
        depth_50th_percentile
        depth_75th_percentile
        depth_90th_percentile

    chksums
        read_names
            Checksum of read names
        sequences
            Checksum of reference sequences
        qualities
            Checksum of read qualities

    coverage_distribution
        range
            .
        coverage
            .
        counts
            .

    first_fragment_qualities
        cycle
            .
        quality
            .

    indel_distribution
        length
            size of indel
        insertions
            number of insertions of this size
        deletions
            number of deletions of this size

    indels_per_cycle
        cyle
            .
        insertions_fwd
            .
        insertions_rev
            .
        deletions_fwd
            .
        deletions_rev
            .

    insert_sizes
        pairs_total
            .
        inward_pairs
            .
        outward_pairs
            .
        other_pairs
            .

    last_fragment_GC_content
        percent_GC
            .
        count
            .

    last_fragment_qualities
        cycle
            .
        quality
            .

    read_lengths
        read_length
            .
        count
            .

    summary_numbers
        category
            category, such as "Raw total sequences"
        count
            number of counts in a category

    Some notes on summaries:

    mismatches
       based on NM tag. This includes indels.
    bases mapped
       based on sequence length as recorded in BAM file. Only mapped
       reads are included and only the length of the primary alignment
       is used to avoid double counting.
    bases mapped (cigar)
       mapped bases determined by CIGAR string. Include both MATCH
       and INSERT states.
    total length
       sum of query lengths (including any sequence that has been
       clipped) including mapped and unmapped reads. Only primary
       alignment is used to avoid double counting.

    Samtools ignores secondary alignments for most statistics.

    """

    name = "samtools_stats"

    _map_section_to_table = {
        'CHK': ("samtools_stats_chksums",
                ("read_names", "sequences", "qualities")),
        'SN': ("samtools_stats_summary_numbers",
               ("category", "count")),
        'FFQ': ("samtools_stats_first_fragment_qualities",
                ("cycle", "VAR_quality")),
        'LFQ': ("samtools_stats_last_fragment_qualities",
                ("cycle", "VAR_quality")),
        'GCF': ("samtools_stats_first_fragment_GC_content",
                ("percent_GC", "count")),
        'GCF': ("samtools_stats_last_fragment_GC_content",
                ("percent_GC", "count")),
        'GCC': ("samtools_stats_ACGT_PER_CYCLE",
                ("cycle", "percent_A", "percent_C", "percent_G", "percent_T")),
        'IS': ("samtools_stats_insert_sizes",
               ("pairs_total", "inward_pairs",
                "outward_pairs", "other_pairs")),
        'RL': ("samtools_stats_read_lengths",
               ("read_length", "count")),
        'ID': ("samtools_stats_indel_distribution",
               ("length", "insertions", "deletions")),
        'IC': ("samtools_stats_indels_per_cycle",
               ("cyle", "insertions_fwd", "insertions_rev",
                "deletions_fwd", "deletions_rev")),
        'COV': ("samtools_stats_coverage_distribution",
                ("range", "coverage", "counts")),
        'GCD': ("samtools_stats_GC_depth",
                ("percent_GC",
                 "unique_percentile",
                 "depth_10th_percentile",
                 "depth_25th_percentile",
                 "depth_50th_percentile",
                 "depth_75th_percentile",
                 "depth_90th_percentile"))
    }

    reference_fasta = None
    reference_fasta_map = None
    target_regions = None
    target_regions_map = None

    ignore_missing_reference = False

    def __call__(self, *args, **kwargs):
        self.tablenames = [x[0] for x in list(self._map_section_to_table.values())]
        MetricRunner.__call__(self, *args, **kwargs)

    def run(self, infile, outfile, params):

        options = []
        reference_fasta = params.reference_fasta
        reference_fasta_map = build_reference_fasta_map(params.reference_fasta_map)
        reference_label = None
        use_target_regions = True
        if params.reference_fasta:
            map_path2name = dict([(x[1], x[0]) for x in list(reference_fasta_map.items())])
            if params.reference_fasta == "auto":

                fasta = resolve_argument(
                    list(reference_fasta_map.values()), ",").split(",")

                reference_fasta, diffs = get_reference_for_bam(
                    infile,
                    fastafiles=fasta)

                if reference_fasta:
                    options.append("--ref-seq {}".format(reference_fasta))
                    reference_label = map_path2name[reference_fasta]
                elif diffs:
                    E.warn(
                        "attempted to detect reference fasta, but unable to do so. "
                        "diffs: {}".format(diffs))
                else:
                    E.warn("sequence dict is empty, BAM likely to be empty. "
                           "target_regions will be ignored")
                    use_target_regions = False
            else:
                options.append("--ref-seq {}".format(params.reference_fasta))
                reference_label = map_path2name.get(params.reference_fasta, None)

        if params.target_regions and use_target_regions:
            target_regions = get_associated_file(params,
                                                 reference_label,
                                                 "target_regions")
            # convert to 1-based coordinates and decompress
            if target_regions.endswith(".bed.gz"):
                target_regions = (
                    "<(zcat {} "
                    "| awk '{{printf(\"%%s\\t%%i\\t%%i\\n\", $1, $2+1, $3)}}')".format(
                        target_regions))
            options.append("--target-regions {}".format(target_regions))

        options = " ".join(options)
        if not os.path.exists(outfile + ".tmp"):
            try:
                retval = P.run("{params.path} stats "
                               "{self.options} "
                               "{options} "
                               "{infile} "
                               "2> {outfile}.log "
                               "> {outfile}.tmp; "
                               .format(**locals()), job_memory="16G")
            except OSError as e:
                E.warn("input file {} gave the following errors: {}".format(
                    infile, str(e)))
                return None
        else:
            retval = None

        def split_output(lines):
            is_comment = True
            section, body = None, []
            for line in lines:
                if line.startswith("#"):
                    if body:
                        yield section, body
                    body = []
                    is_comment = True
                else:
                    # the following preserves new-line
                    line = re.sub("\t#.*", "", line)
                    fields = line[:-1].split("\t")
                    section = fields[0]
                    body.append(fields[1:])
                    is_comment = False

            if body:
                yield section, body

        # split into separate files for upload
        with IOTools.open_file(outfile + ".tmp") as inf:
            for section, body in split_output(inf):
                try:
                    tablename, columns = self._map_section_to_table[section]
                except KeyError:
                    continue

                output_file = self.map_table_to_file(tablename, outfile)
                with IOTools.open_file(output_file, "w") as outf:

                    if len(columns) > 1 and columns[1].startswith("VAR_"):
                        outf.write("{}\t{}\n".format(
                            columns[0], columns[1][4:]))
                        for data in body:
                            outf.write("{}\t{}\n".format(
                                data[0], ",".join(data)))
                    else:
                        outf.write("\t".join(columns) + "\n")
                        # remove first column, which contains the identifier
                        outf.write("\n".join(
                            ["\t".join(x) for x in body]) + "\n")

        os.rename(outfile + ".tmp", outfile)

        return retval


class run_metric_bam_fastqc(MetricRunner):
    name = "bam_fastqc"
    path = "fastqc"

    table_suffixes = [
        "basic_statistics",
        "per_base_sequence_quality",
        "per_sequence_quality_scores",
        "per_base_sequence_content",
        "per_sequence_gc_content",
        "per_base_n_content",
        "per_tile_sequence_quality",
        "sequence_length_distribution",
        "sequence_duplication_levels",
        "overrepresented_sequences",
        "adapter_content",
        "kmer_content",
        "summary"]

    blob_globs = [os.path.join("result_fastqc",
                               "Images",
                               "*.png")]

    def __init__(self, *args, **kwargs):
        self.tablenames = ["{}_{}".format(self.name, x)
                           for x in self.table_suffixes]
        MetricRunner.__init__(self, *args, **kwargs)

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stdout=True).strip()
        return re.search(r"FastQC (\S+)", help_string).groups()[0]

    def run(self, infile, outfile, params):
        # TODO: bam_fastqc_sequence_length_distribution.tsv may
        # contain ranges such as '30-31'. Convert to beginning of
        # range like in this perl command:
        #
        # perl -p -i -e "s/\-\d+//"
        # *.dir/bam_fastqc.dir/bam_fastqc.tsv.bam_fastqc_sequence_length_distribution.tsv

        if infile.endswith(".gz"):
            prefix = IOTools.snip(os.path.basename(infile[:-3]))
        else:
            prefix = IOTools.snip(os.path.basename(infile))

        outdir = os.path.dirname(outfile)

        datafile = os.path.join(
            outdir,
            "{}_fastqc".format(prefix),
            "fastqc_data.txt")

        if not os.path.exists(datafile):
            if not os.path.exists(outdir):
                os.makedirs(outdir)

            retval = P.run(
                "{params.path} "
                "{params.options} "
                "--extract "
                "--outdir {outdir} "
                "{infile} "
                ">& {outfile} ".format(**locals()), **params._asdict())
        else:
            IOTools.touch_file(outfile)
            retval = None

        def _split_output(lines):
            body, header, section, status = [], None, None, None
            for line in lines:
                if line.startswith("##FastQC"):
                    continue
                elif line.startswith("#"):
                    header, body = line[1:-1].split("\t"), []
                elif line.startswith(">>END_MODULE"):
                    yield section, header, body, status
                    body, header, section, status = [], None, None, None
                elif line.startswith(">>"):
                    section, status = line[2:-1].split("\t")
                else:
                    fields = line[:-1].split("\t")
                    body.append(fields)

        # split into separate files for upload
        summary_data = []
        with IOTools.open_file(datafile) as inf:
            for section, header, body, status in _split_output(inf):
                if len(body) == 0:
                    continue
                summary_data.append((section, status))
                tablename = "{}_".format(self.name) + re.sub(" ", "_", section).lower()
                if tablename not in self.tablenames:
                    raise ValueError(
                        "unknown tablename {}, expected one of {}".
                        format(tablename, self.tablenames))
                output_file = ".".join((outfile, tablename, "tsv"))
                with open(output_file, "w") as outf:
                    outf.write("\t".join([x.lower() for x in header]) + "\n")
                    # remove first column, which contains the identifier
                    outf.write("\n".join(
                        ["\t".join(x) for x in body]) + "\n")

        output_file = ".".join((outfile,
                                "{}_summary".format(self.name), "tsv"))
        with IOTools.open_file(output_file, "w") as outf:
            outf.write("section\tstatus\n")
            for section, status in summary_data:
                outf.write("{}\t{}\n".format(section, status))

        return retval


class MetricRunnerPicard(MetricRunner):
    path = "/data/install/Free/picard-tools-1.140/picard.jar"

    def get_version(self):
        return E.run(
            "java -jar {self.path} SortVcf --version".format(**locals()),
            return_stderr=True).strip()


class run_metric_samtools_mpileup_histogram(MetricRunnerSamtools):
    """run samtools mpileup on a :term:`bam` file and compute
    histogram with coverage statistics.

    To include all reads in the pileup, use the options
    `-Q 0 -B -A`.

    :param infile: Input file in :term:`bam` format
    :param outfile: Output file in :term:`tsv` format

    *Columns*

    total_counts
        total number of reads aligning at a position
    reference_counts
        number of reference reads aligning to position
    position_counts
        number of positions that have total_counts/reference_counts
    """

    name = "samtools_mpileup_histogram"
    reference_fasta = None

    def run(self, infile, outfile, params):

        if params.reference_fasta is None:
            raise ValueError("samtools_mpileup_histogram requires a reference_fasta option")
        # samtools pileup format:
        # . and , are matches on forward/reverse strand
        # a ^ followed by a quality char is a read start. To avoid
        # double counting, remove all ^. and ^, dimers.

        retval = P.run(
            "{params.path} mpileup "
            "-f {params.reference_fasta} "
            "{params.options} "
            "{infile} "
            "2> {outfile}.pileup.log "
            "| awk 'BEGIN {{ "
            "   printf(\"total_counts\\treference_counts\\tposition_counts\\n\") "
            "}} "
            "{{ split($5, c, \"\"); "
            "  total_counts=$4; "
            "  for (i=1; i < length($5); i++) {{"
            "    counts[c[i]] += 1; counts[c[i]c[i+1]] += 1 "
            "  }}; "
            "  counts[c[length($5)]] += 1; "
            "  r=counts[\".\"]+counts[\",\"]-counts[\"^.\"]-counts[\"^,\"]; "
            "  if (total_counts < r) {{ "
            "     printf(\"invalid counts %%i < %%i (%%i, %%i, %%i, %%i)\\n\", "
            "     total_counts, r, "
            "     counts[\".\"], counts[\",\"], counts[\"^.\"], counts[\"^,\"]) "
            "     > \"{outfile}.err\"; delete counts; next; "
            "   }} "
            "   hist[total_counts][r] += 1; "
            "   delete counts; "
            "   next;"
            "}} "
            "END {{ "
            "  for (i in hist) for (j in hist[i]) "
            "    printf(\"%%i\\t%%i\\t%%i\\n\", i, j, hist[i][j]); "
            "}} '"
            "2> {outfile}.awk.log "
            "> {outfile}; "
            .format(**locals()))

        return retval


class run_metric_picard_collectmultiplemetrics(MetricRunnerPicard):
    name = "picard_collectmultiplemetrics"

    tablenames = [
        "picard_alignment_summary_metrics",
        "picard_insert_size_metrics",
        "picard_insert_size_histogram",
        "picard_base_distribution_by_cycle_metrics",
        "picard_quality_by_cycle_histogram",
        "picard_quality_distribution_histogram"
    ]

    def __call__(self, *args, **kwargs):
        MetricRunner.__call__(self, *args, **kwargs)

    def run(self, infile, outfile, params):

        if "reference_fasta" in params._fields:
            reference_fasta = "REFERENCE_SEQUENCE={}".format(
                params.reference_fasta)
        else:
            reference_fasta = ""

        # command can fail when no output is produced, but still produce output
        # 12G is required for java overhead
        retval = P.run(
            "java -Xmx8000m -jar {params.path} "
            "CollectMultipleMetrics "
            "{reference_fasta} "
            "INPUT={infile} "
            "TMP_DIR=%(tmpdir)s "
            "{params.options} "
            "OUTPUT={outfile} "
            ">& {outfile} ".format(**locals()),
            job_memory="12G",
            ignore_errors=True)

        def get_section(section, data):
            pattern = "## {}".format(section)
            keep = False
            result = []
            for line in data:
                if line.startswith("##"):
                    if line.startswith(pattern):
                        keep = True
                    else:
                        keep = False
                if keep:
                    result.append(line)
            return result

        for tablename in self.tablenames:
            filename = re.sub("histogram", "metrics", tablename)
            raw = filename[len("picard_"):]
            src = outfile + "." + raw
            dest = outfile + "." + tablename + ".tsv"

            if not os.path.exists(src):
                E.warn("no file {}, ignored".format(src))
                continue

            with IOTools.open_file(src) as inf:
                data = inf.readlines()

            if tablename.endswith("metrics"):
                data = get_section("METRICS", data)
            elif tablename.endswith("histogram"):
                data = get_section("HISTOGRAM", data)

            with IOTools.open_file(dest, "w") as outf:
                outf.write("".join(data))

        return retval


class run_metric_daisy_bam_compare_alignments(MetricRunner):
    """run daisy bam_compare_alignments on a :term:`BAM` file.

    :param infile: Input file in :term:`BAM` format
    :param outfile: Multiple output files in compressed :term:`tsv`
         format

    """

    name = "daisy_bam_compare_alignments"
    path = "daisy"

    tablenames = ["daisy_bam_compare_alignments_mapped",
                  "daisy_bam_compare_alignments_overlap",
                  "daisy_bam_compare_alignments_summary"]

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        if not os.path.exists(params.reference_bam):
            raise OSError(
                "reference bam file {} does not exist".format(
                    params.reference_bam))

        tmpdir = P.get_temp_dir(clear=True)

        statement = (
            "mkdir {tmpdir}; "
            "samtools sort -n {infile} > {tmpdir}/comp.bam; "
            "samtools sort -n {params.reference_bam} > {tmpdir}/ref.bam; "
            "{params.path} bam-compare-alignments "
            "--output-filename-pattern={outfile}.daisy_bam_compare_alignments_%%s.tsv "
            "{params.options} "
            "--input-bam={tmpdir}/comp.bam "
            "--reference-bam={tmpdir}/ref.bam "
            ">& {outfile}; "
            "rm -rf {tmpdir}; "
            .format(**locals()))

        retval = P.run(statement)

        return retval


class run_metric_daisy_bam2stats(MetricRunner):
    """run daisy bam2stats tools on a :term:`BAM` file.

    :param infile: Input file in :term:`BAM` format
    :param outfile: Multiple output files in :term:`tsv`
         format

    """

    name = "daisy_bam2stats"
    path = "daisy"

    tablenames = ["daisy_bam2stats_counts",
                  # "daisy_bam2stats_details",  # ignore details table, potentially very large
                  "daisy_bam2stats_histogram",
                  "daisy_bam2stats_summaries",
                  "daisy_bam2stats_nm",
                  "daisy_bam2stats_mapq"]

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        statement = (
            "{params.path} bam2stats "
            "--output-filename-pattern={outfile}.daisy_bam2stats_%%s.tsv "
            "--force-output "
            "{params.options} "
            "{infile} "
            "--log {outfile} "
            "> {outfile}.daisy_bam2stats_counts.tsv"
            .format(**locals()))

        retval = P.run(statement)

        return retval


class run_metric_daisy_bam2reference(MetricRunner):
    """deduce reference sequence for a BAM file.

    .. note::

        This tool runs locally.

    :param infile: Input file in :term:`BAM` format
    :param outfile: Multiple output files in :term:`tsv`
         format

    """

    name = "daisy_bam2reference"
    path = "daisy"

    reference_fasta_map = None

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        if params.reference_fasta_map is None:
            raise ValueError("bam2reference requires a reference sequence map")

        reference_fasta_map = build_reference_fasta_map(params.reference_fasta_map)

        fasta = resolve_argument(list(reference_fasta_map.values()), ",").split(",")
        retval, diff = get_reference_for_bam(infile, fasta)
        if retval is None:
            if diff is None:
                retval = "corrupted"
            else:
                retval = "unknown"
                E.debug("differences: {}".format(str(diff)))
            path = ""
        else:
            map_path2name = dict([(x[1], x[0]) for x in list(reference_fasta_map.items())])
            path = map_path2name.get(retval, os.path.basename(retval))

        with IOTools.open_file(outfile, "w") as outf:
            outf.write("filename\treference\tpath\n")
            outf.write("\t".join((infile, retval, path)) + "\n")

        return None


class run_metric_daisy_bam2window_stats(MetricRunner):
    """run daisy bam2window-stats tools on a :term:`BAM` file.

    :param infile: Input file in :term:`BAM` format
    :param outfile: Multiple output files in :term:`tsv`
         format

    """

    name = "daisy_bam2window_stats"
    path = "daisy"

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        statement = (
            "{params.path} bam2window-stats "
            "{params.options} "
            "{infile} "
            "--log {outfile}.log "
            "2> {outfile}.err "
            "> {outfile}"
            .format(**locals()))

        retval = P.run(statement)

        return retval


class run_metric_bam_bedtools_coverage(MetricRunner):
    """run bedtools coverage on a :term:`BAM` file.

    :param infile: Input file in :term:`BAM` format
    :param outfile: Multiple output files in :term:`tsv`
         format

    """

    name = "bam_bedtools_coverage"
    path = "bedtools"

    reference_vcf = None

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):
        if params.reference_vcf is None:
            raise ValueError("tool requires reference_vcf to be set")

        statement = (
            "{params.path} coverage "
            "-a {infile} "
            "-b {params.reference_vcf} "
            "{params.options} "
            "> {outfile}"
            .format(**locals()))

        retval = P.run(statement)

        return retval


class run_metric_samtools_depth(MetricRunnerSamtools):
    """run samtools depth on a :term:`BAM` file.

    :param infile: Input file in :term:`BAM` format
    :param outfile: Multiple output files in :term:`tsv`
         format

    """

    name = "samtools_depth"
    path = "samtools"

    def run(self, infile, outfile, params):
        statement = (
            "{params.path} depth "
            "-a "
            "{params.options} "
            "{infile} "
            "> {outfile}"
            .format(**locals()))

        return P.run(statement)


class run_metric_daisy_bam2depth(MetricRunner):
    """run daisy bam2depth tool on a :term:`BAM` file.

    :param infile: Input file in :term:`BAM` format
    :param outfile: A tsv separated file with depth histogram.
    """

    name = "daisy_bam2depth"
    path = "daisy"

    reference_fasta = None

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        if params.reference_fasta is None:
            raise ValueError("daisy_bam2depth requires a reference_fasta option")

        statement = (
            "{params.path} bam2depth "
            "--input-filename-fasta={params.reference_fasta} "
            "{params.options} "
            "{infile} "
            "--log {outfile}.log "
            "2> {outfile}.err "
            "> {outfile}"
            .format(**locals()))

        return P.run(statement)
