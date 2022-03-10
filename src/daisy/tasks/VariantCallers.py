import os
import re
import pysam
import shutil

from daisy.toolkit import parse_region_string
from .Runner import resolve_argument
from .ToolRunner import ToolRunner
import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools
from .BAMMetrics import get_reference_for_bam
from distutils.version import LooseVersion


def get_reference(params):
    bams = resolve_argument(params.bam, ",").split(",")
    fastas = [params.reference_fasta]
    try:
        fastas.extend(params.alt_reference_fasta.split(","))
    except AttributeError:
        pass

    # check only first bam file, assume that they are consistent,
    # but this could be a separate check.
    fasta, diff = get_reference_for_bam(bams[0], fastas)
    if diff:
        raise ValueError(
            "could not match BAM file to references: "
            "{}".format(",".join(
                ["{}:missing_contig={},length_mismatch={}".format(*x)
                 for x in diff])))
    return fasta


class VariantCaller(ToolRunner):
    expected = ["reference_fasta", "bam"]
    output = "result.vcf.gz"

    def run_parallel_by_region(self,
                               statement,
                               bam,
                               reference_fasta,
                               outfile,
                               params,
                               region_option,
                               **kwargs):
        """
        region_option : string
            Needs to contain interpolatable variables contig, start and end. This
            should be long-form option.

        """
        # NOTE: at the moment this method assumes that the variant caller
        # accepts pythonic coordinates
        retvals = []
        statements = []
        files_to_merge = []
        jobsfile = outfile + ".jobs"

        # check if region limited in statement=
        if "{contig}" not in region_option or "{start}" not in region_option or "{end}" not in region_option:
            raise ValueError(
                "region_option requires {{contig}}, {{start}} and {{end}} fields, got {}".format(
                    region_option))

        region_prefix = region_option[:region_option.index("{")]
        if re.search(region_prefix, statement):
            region = re.search("{}[= ]*(\S+)".format(region_prefix),
                               params.options).groups()[0]
            filter_contig, filter_start, filter_end = parse_region_string(
                region)
        else:
            filter_contig, filter_start, filter_end = None, None, None

        statement = re.sub("{}[= ]*\S+".format(region_prefix), "",
                           statement)
        statement += " {}".format(region_option)

        statements = []

        with pysam.FastaFile(reference_fasta) as fastaf:
            for contig, length in zip(fastaf.references,
                                      fastaf.lengths):
                if filter_contig and contig != filter_contig:
                    continue
                begin_range = filter_start if filter_start else 0
                end_range = filter_end if filter_end else length

                for start in range(begin_range, end_range,
                                   params.chunk_size):
                    fn = os.path.join(
                        outfile + ".chunk_{}_{:08}.vcf.gz".format(
                            contig, start))
                    files_to_merge.append(fn)
                    if os.path.exists(fn):
                        continue
                    end = min(start + params.chunk_size, length)
                    statements.append(
                        statement.format(contig=contig, start=start, end=end) +
                        " 2> {fn}.log "
                        "| bgzip "
                        "> {fn}\n".format(**locals()))

        retvals = P.run(statements, job_array=True, **kwargs)

        fn = " ".join(files_to_merge)
        statement = (
            "zcat {fn} "
            "| vcffirstheader "
            "2> {outfile}.vcffirstheader.log "
            "| vcfstreamsort -w 1000 "
            "2> {outfile}.vcfstreamsort.log "
            "| vcfuniq "
            "2> {outfile}.vcfuniq.log "
            "| bgzip "
            "2> {outfile}.bgzip.log "
            "> {outfile}; "
            "tabix -p vcf {outfile} "
            "2> {outfile}.tabix.log; "
            "rm -f {fn} "
            "".format(**locals()))

        retvals.extend(P.run(statement))
        return retvals


class run_tool_freebayes(VariantCaller):
    name = "freebayes"
    path = "freebayes"

    parallel = False
    # chunk size in bp for parallelization.
    chunk_size = 1000000

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stdout=True).strip()
        return re.search("version:\s+(\S+)", help_string).groups()[0]

    def run(self, outfile, params):

        bam = resolve_argument(params.bam, sep=" ")
        reference_fasta = get_reference(params)

        if params.parallel:
            statements = []
            files_to_merge = []
            jobsfile = outfile + ".jobs"

            if re.search("--region", params.options):
                region = re.search("--region[= ]*(\S+)",
                                   params.options).groups()[0]
                filter_contig, filter_start, filter_end = parse_region_string(
                    region)
            else:
                filter_contig, filter_start, filter_end = None, None, None

            plain_options = re.sub("--region[= ]\S+", "", params.options)

            statements = []

            with pysam.FastaFile(reference_fasta) as fastaf:
                for contig, length in zip(fastaf.references,
                                          fastaf.lengths):
                    if filter_contig and contig != filter_contig:
                        continue
                    begin_range = filter_start if filter_start else 0
                    end_range = filter_end if filter_end else length

                    for start in range(begin_range, end_range, params.chunk_size):
                        fn = os.path.join(
                            outfile + ".chunk_{}_{:08}.vcf.gz".format(contig, start))
                        files_to_merge.append(fn)
                        if os.path.exists(fn):
                            continue
                        end = min(start + params.chunk_size, length)
                        statements.append(
                            "{params.path} "
                            "--fasta-reference {reference_fasta} "
                            "--region {contig}:{start}-{end} "
                            "{plain_options} "
                            "{bam} "
                            "2> {fn}.log "
                            "| bgzip "
                            "> {fn}\n".format(**locals()))

            retvals = P.run(statements, job_array=True)

            fn = " ".join(files_to_merge)
            statement = (
                "zcat {fn} "
                "| vcffirstheader "
                "2> {outfile}.vcffirstheader.log "
                "| vcfstreamsort -w 1000 "
                "2> {outfile}.vcfstreamsort.log "
                "| vcfuniq "
                "2> {outfile}.vcfuniq.log "
                "| bgzip "
                "2> {outfile}.bgzip.log "
                "> {outfile}; "
                "tabix -p vcf {outfile} "
                "2> {outfile}.tabix.log; "
                "rm -f {fn} "
                "".format(**locals()))

            retvals.extend(P.run(statement))

        else:
            # limit number of jobs to node to limit I/O
            job_threads = 2

            retvals = P.run("{params.path} "
                            "--fasta-reference {reference_fasta} "
                            "{params.options} "
                            "{bam} "
                            "2> {outfile}.log "
                            "| bgzip "
                            "> {outfile}; "
                            "tabix -p vcf {outfile}"
                            .format(**locals()), **params._asdict())

        if "set_filter_exclude" in params._fields:
            with IOTools.open_file(outfile + ".header.vcf", "w") as outf:
                outf.write(
                    "##FILTER=<ID=HARD,Description=\"Variant fails hard filters: {}\"> "
                    .format(params.set_filter_exclude))
            job_threads = 1
            # note to include in first step as these will be set to value "HARD"
            retvals.extend(
                P.run("bcftools query "
                      "--include \"{params.set_filter_exclude}\" "
                      "-f \"%%CHROM\\t%%POS\\tHARD\\n\" "
                      "{outfile}.save.vcf.gz "
                      "| bgzip > {outfile}.tab.gz; "
                      "tabix -s 1 -b 2 -e 2 {outfile}.tab.gz; "
                      "bcftools annotate "
                      "-a {outfile}.tab.gz "
                      "-c CHROM,POS,FILTER "
                      "--header-lines {outfile}.header.vcf "
                      "{outfile}.save.vcf.gz "
                      "| bgzip > {outfile}.new.vcf.gz; "
                      "mv {outfile}.new.vcf.gz {outfile}; "
                      "tabix -f -p vcf {outfile} ".format(**locals())))

        return retvals


class run_tool_samtools(VariantCaller):

    name = "samtools"
    path = "samtools"

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]

    def run(self, outfile, params):

        bam = resolve_argument(params.bam, sep=" ")
        reference_fasta = get_reference(params)

        return P.run("{params.path} mpileup "
                     "-v "
                     "-f {reference_fasta} "
                     "{params.options} "
                     "{bam} "
                     "2> {outfile}.log "
                     "> {outfile}; "
                     "tabix -p vcf {outfile} "
                     .format(**locals()))


class run_tool_bcftools(VariantCaller):

    name = "bcftools"
    path = "bcftools"
    path_samtools = "samtools"
    samtools_options = ""
    options = "--multiallelic-caller"

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]

    def run(self, outfile, params):

        bam = resolve_argument(params.bam, sep=" ")
        reference_fasta = get_reference(params)

        # warning: requires -m or -c in the options
        if "--multiallelic-caller" not in params.options and \
           "-m" not in params.options and \
           "-c" not in params.options and \
           "--consensus-caller" not in params.options:
            E.warn("bcftools call requires -m or -c, got {}".format(
                params.options))

        # limit number of jobs to node to limit I/O
        job_threads = 4

        return P.run("{params.path_samtools} mpileup "
                     "-ug "
                     "-f {reference_fasta} "
                     "{params.samtools_options} "
                     "{bam} "
                     "2> {outfile}.pileup.log "
                     "| {params.path} call "
                     "--variants-only "
                     "--output-type z "
                     "{params.options} "
                     "2> {outfile}.call.log "
                     "> {outfile}; "
                     "tabix -p vcf {outfile} "
                     .format(**locals()))


class run_tool_octopus(VariantCaller):

    name = "octopus"
    path = "octopus"

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stdout=True).strip()
        return re.search("octopus (.+)$",
                         help_string).groups()[0]

    def run(self, outfile, params):

        if "--threads" in params.options:
            job_threads = int(re.search("--threads[= ]*(\d+)",
                                        params.options).groups()[0])

        bam = resolve_argument(params.bam, sep=" ")
        reference_fasta = get_reference(params)
        tmpdir = P.get_temp_filename(clear=True)

        # Octopus seems to have issues if ulimit's are set.
        return P.run("{params.path} "
                     "--working-directory $TMPDIR "
                     "--reads {bam} "
                     "--reference {reference_fasta} "
                     "--output {outfile} "
                     "{params.options} "
                     ">& {outfile}.log "
                     .format(**locals()),
                     job_memory="unlimited")


class run_tool_platypus(VariantCaller):

    name = "platypus"
    path = "Platypus.py"

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stdout=True).strip()
        return "undefined"

    def run(self, outfile, params):

        if "--nCPU" in params.options:
            job_threads = int(re.search("--nCPU\s*(\d+)",
                                        params.options).groups()[0])

        bam = resolve_argument(params.bam)
        reference_fasta = get_reference(params)

        tmpfile = P.get_temp_filename(clear=True)

        return P.run("{params.path} callVariants "
                     "--bamFiles {bam} "
                     "--refFile {reference_fasta} "
                     "--output {tmpfile} "
                     "{params.options} "
                     ">& {outfile}.log; "
                     "bgzip {tmpfile}; "
                     "tabix -p vcf {tmpfile}.gz; "
                     "mv {tmpfile}.gz  {outfile}; "
                     "mv {tmpfile}.gz.tbi {outfile}.tbi; "
                     .format(**locals()))


class run_tool_breakdancer(VariantCaller):

    name = "breakdancer"
    path = "breakdancer-max"
    path_cfg = ("/data/install/Free/breakdancer-1.3.6/lib/"
                "breakdancer-max1.4.5-unstable-60-3876c5f/bam2cfg.pl")

    Cfg_options = ""

    def get_version(self):
        help_string = E.run("{self.path} -h".format(**locals()),
                            return_stderr=True).strip()
        if help_string and "not found" not in help_string:
            return re.search("breakdancer-max version (.+)$",
                             help_string).groups()[0]
        else:
            raise ValueError("breakdancer not found at/as {}: {}".format(
                self.path, help_string))

    def run(self, outfile, params):

        bam = resolve_argument(params.bam, sep=" ")

        retvals = []
        retvals.append(P.run("{params.path_cfg} {params.cfg_options} "
                             "{bam} "
                             "2> {outfile}.cfg.log "
                             "> {outfile}.cfg ".format(**locals())))

        retvals.append(P.run("{params.path} "
                             "-d {outfile} "
                             "-g {outfile}.bed "
                             "{params.options} "
                             "{outfile}.cgf "
                             "2> {outfile}.log "
                             "> {outfile} ".format(**locals())))
        return retvals


class run_tool_delly(VariantCaller):

    name = "delly"
    path = "delly"
    path_vcf_concat = "vcf-concat"
    path_vcf_sort = "vcf-sort"
    variant_types = "DEL,INS,DUP"

    def get_version(self):
        help_string = E.run("{self.path} ".format(**locals()),
                            return_stdout=True,
                            on_error="ignore").strip()
        if help_string:
            return re.search("Delly \(Version: (\S+)\)",
                             help_string).groups()[0]
        else:
            raise ValueError("delly not found at/as {}".format(self.path))

    def run(self, outfile, params):

        bam = resolve_argument(params.bam, sep=" ")
        reference_fasta = get_reference(params)

        stments, retvals = [], []
        variant_types = [x.strip() for x in params.variant_types.split(",")]

        for variant_type in variant_types:
            stments.append(
                "{params.path} "
                "--type {variant_type} "
                "--genome {reference_fasta} "
                "--outfile {outfile}.{variant_type}.vcf "
                "{params.options} "
                "{bam} "
                ">& {outfile}.{variant_type}.log; "
                "bgzip -f {outfile}.{variant_type}.vcf; "
                "tabix -f -p vcf {outfile}.{variant_type}.vcf.gz".format(
                    **locals()))

        retvals.extend(P.run(stments))

        vcf_files = " ".join([outfile + "." + x + ".vcf.gz"
                              for x in variant_types])
        retvals.append(P.run(
            "{params.path_vcf_concat} "
            "{vcf_files} "
            "| {params.path_vcf_sort} "
            "| bgzip "
            "> {outfile}; "
            "tabix -fp vcf {outfile}".format(**locals())))

        return retvals


class run_tool_lumpy(VariantCaller):
    name = "lumpy"
    path = "lumpyexpress"

    def get_version(self):
        help_string = E.run("{self.path} ".format(**locals()),
                            return_stdout=True,
                            on_error="ignore").strip()
        # lumpy express without arguments ends in error
        if help_string:
            raise NotImplementedError()
            return re.search(r"lumpy \(Version: (\S+)\)",
                             help_string).groups()[0]
        else:
            raise ValueError("lumpy not found at/as {}".format(self.path))

    def run(self, outfile, params):

        bam = resolve_argument(params.bam, sep=",")

        # "-T {outfile}.tmpdir -k "

        outfile = IOTools.snip(outfile, ".gz")
        # note that lumpy removes the temporary directory
        # after running, thus make sure it is unique and exists
        return P.run("{params.path} "
                     "-B {bam} "
                     "-o {outfile} "
                     "-T %(tmpdir)s_{self.__name__} "
                     "-v "
                     "{params.options} "
                     ">& {outfile}.log; "
                     "vcf-sort {outfile} "
                     "| bgzip > {outfile}.gz; "
                     "tabix -p vcf {outfile}.gz"
                     .format(**locals()))


class run_tool_delly_cnv(run_tool_delly):
    """run delly in CNV mode.

    Use delly to predict structural variants, but
    output these as a BED-file.
    """
    name = "delly_cnv"
    path = "delly"
    path_vcf_concat = "vcf-concat"
    path_vcf_sort = "vcf-sort"
    path_bcftools = "bcftools"
    variant_types = "DEL,INS,DUP"
    bcftools_options = ""

    output = "result.bed.gz"

    def run(self, outfile, params):

        retvals = []
        prefix = IOTools.snip(outfile, ".bed.gz")
        vcffile = prefix + ".vcf.gz"
        if not os.path.exists(vcffile):
            retvals.extend(run_tool_delly.run(self, vcffile, params))

        statements = []

        statements.append(
            "{self.path_bcftools} query "
            "{params.bcftools_options} "
            "-f \"%%CHROM\\t%%POS\\t%%END\\t%%SVTYPE\\n\" "
            "{vcffile} "
            "| awk -v OFS='\\t' '$3 != \".\" {{ switch ($4) {{"
            "case \"DEL\": $5=0; break; "
            "case \"DUP\": $5=3; break; "
            "case \"INS\": next; break; "
            "}}; print }}' "
            "| bgzip "
            "> {outfile}".format(**locals()))
        statements.append(
            "tabix -f -p bed {outfile}".format(**locals()))

        statement = "; ".join(statements)
        retvals.append(P.run(statement))

        return retvals


class run_tool_daisy_ont_phase_variants(VariantCaller):

    name = "daisy_ont_phase_variants"
    path = "daisy"

    expected = ["bam", "vcf"]

    def get_version(self):
        return "builtin"

    def run(self, outfile, params):

        bam = resolve_argument(params.bam, sep=",")

        retval = P.run(
            "{params.path} ont-phase-variants "
            "{params.options} "
            "--input-bam={bam} "
            "--input-vcf={params.vcf} "
            "--log={outfile}.log "
            "| bgzip "
            "> {outfile}; "
            "tabix -p vcf {outfile}".format(**locals()))

        return retval


class run_tool_manta_strelka(VariantCaller):
    name = "manta_strelka"
    path_manta = "configManta.py"
    path_strelka = "configureStrelkaGermlineWorkflow.py"

    manta_options = ""
    strelka_options = ""

    manta_threads = 4
    strelka_threads = 4

    manta_memory_gb = 8
    strelka_memory_gb = 8

    is_exome = False

    def get_version(self):
        help_string1 = E.run("{self.path_manta} --version".format(**locals()),
                             return_stdout=True).strip()

        help_string2 = E.run("{self.path_strelka} --version".format(**locals()),
                             return_stdout=True).strip()
        return "-".join([help_string1, help_string2])

    def run(self, outfile, params):

        bams = resolve_argument(params.bam, " ").split(" ")
        reference_fasta = get_reference(params)
        outdir = os.path.dirname(os.path.abspath(outfile))

        if "--numberOfJobs" in params.options:
            job_threads = int(re.search("--numberOfJobs\s*(\d+)",
                                        params.options).groups()[0])

        if "Somatic" in params.path_strelka:
            if len(bams) != 2:
                raise ValueError(
                    "somatic workflow expects 2 BAM files, but got {}: {}".format(
                        len(bams), bams))
            bamstring = "--normalBam {} --tumorBam {}".format(*bams)
            vcf_output = os.path.join(
                outdir, "strelka.dir", "results", "variants",
                "somatic.snvs.vcf.gz")
        elif "Germline" in params.path_strelka:
            bamstring = " ".join(["--bam={}".format(x) for x in bams])
            vcf_output = os.path.join(
                outdir, "strelka.dir", "results", "variants",
                "variants.vcf.gz")
        else:
            raise NotImplementedError(
                "strelka workflow {} not supported".format(
                    params.path_strelka))

        indel_candidates_output = os.path.join(outdir, "manta.dir",
                                               "results", "variants",
                                               "candidateSmallIndels.vcf.gz")
        retvals = []
        if not os.path.exists(indel_candidates_output):
            statement = (
                "{params.path_manta} "
                "{bamstring} "
                "--referenceFasta={reference_fasta} "
                "{params.manta_options} "
                "--runDir={outdir}/manta.dir "
                ">& {outdir}/manta-config.log".format(**locals()))

            retvals.extend(P.run(statement, to_cluster=False))

            statement = (
                "{outdir}/manta.dir/runWorkflow.py "
                "--mode=local "
                "--jobs={params.manta_threads} "
                "--memGb={params.manta_memory_gb} "
                ">& {outdir}/manta-run.log".format(**locals()))

            retvals.extend(P.run(
                statement,
                job_threads=params.manta_threads,
                job_memory="{}G".format(params.manta_memory_gb / params.manta_threads)))

        if not os.path.exists(vcf_output):
            statement = (
                "{params.path_strelka} "
                "{bamstring} "
                "--referenceFasta={reference_fasta} "
                "--indelCandidates={indel_candidates_output} "
                "{params.strelka_options} "
                "--runDir={outdir}/strelka.dir "
                ">& {outdir}/strelka-config.log".format(**locals()))

            retvals.extend(P.run(statement, to_cluster=False))

            statement = (
                "{outdir}/strelka.dir/runWorkflow.py "
                "--mode=local "
                "--jobs={params.strelka_threads} "
                "--memGb={params.strelka_memory_gb} "
                ">& {outdir}/strelka-run.log".format(**locals()))

            retvals.extend(P.run(
                statement,
                job_threads=params.strelka_threads,
                job_memory="{}G".format(params.strelka_memory_gb / params.strelka_threads)))

        if "Somatic" in params.path_strelka:
            vcf_indel_output = os.path.join(os.path.dirname(vcf_output),
                                            "somatic.indels.vcf.gz")
            statement = (
                "bcftools concat --allow-overlaps "
                "{vcf_output} {vcf_indel_output} "
                "2> {outfile}.concat.log "
                "| daisy vcf2vcf "
                "--method=add-strelka-genotype "
                "--log={outfile}.add_genotype.log "
                "| bgzip > {outfile}; "
                "tabix -p vcf {outfile}".format(**locals()))
            retvals.extend(P.run(statement))
        else:
            shutil.move(vcf_output, outfile)
            os.symlink(os.path.abspath(outfile), vcf_output)
            shutil.move(vcf_output + ".tbi", outfile + ".tbi")
            os.symlink(os.path.abspath(outfile) + ".tbi", vcf_output + ".tbi")

        return retvals


class run_tool_whatshap_ont_phase_variants(VariantCaller):

    name = "whatshap_ont_phase_variants"
    path = "xxx"

    expected = ["bam", "vcf"]

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stdout=True).strip()
        return help_string

    def run(self, outfile, params):

        bam = resolve_argument(params.bam, sep=",")
        outfile_no_gz = outfile.split(".gz")[0]
        retval = P.run(
            "{params.path} phase "
            "{params.options} "
            "--output {outfile} "
            "--output-read-list {outfile}.read-list "
            "{params.vcf} "
            "{bam} "
            " 2> {outfile}.log && "
            " gzip -d {outfile} && bgzip {outfile_no_gz} && "
            "tabix -p vcf {outfile}".format(**locals()))

        return retval


class run_tool_daisy_bam_pileup2tsv(VariantCaller):

    name = "daisy_bam_pileup2tsv"
    path = "daisy"

    expected = ["bam", "reference_fasta"]

    def get_version(self):
        return "builtin"

    def run(self, outfile, params):

        bam = resolve_argument(params.bam, sep=",")

        retval = P.run(
            "{params.path} bam-pileup2tsv "
            "--reference-fasta={params.reference_fasta} "
            "--log={outfile}.log "
            "--method=depth-vcf "
            "{params.options} "
            "{bam} "
            "| bgzip "
            "> {outfile}; "
            "tabix -p vcf {outfile}".format(**locals()))

        return retval
