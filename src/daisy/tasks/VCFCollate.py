import re
import os
import pysam
import shutil
import pandas

from .CollateRunner import CollateRunner
import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools
import daisy.toolkit as tk


class run_collate_bcftools_merge(CollateRunner):
    path = "bcftools"
    name = "bcftools_merge"

    set_missing_genotype_to_reference = False

    restrict_to_all = False

    remove_fields = None

    block_size = 1000

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]

    def run(self, infiles, outfile, params):

        tmpdir = P.get_temp_dir(clear=True)

        statements = ["mkdir {}".format(tmpdir)]

        if params.remove_fields:
            cleanup_statement = (
                "| {params.path} annotate "
                "-x {params.remove_fields} "
                "2> {outfile}_annotate.log ".format(**locals()))
        else:
            cleanup_statement = ""

        # the current pattern is properly overly specific and
        # substitutes ./. with 0/0
        if params.set_missing_genotype_to_reference:
            set_genotype = "| perl -p -e 's/\.\/\./0\/0/g'"
        else:
            set_genotype = ""

        with IOTools.open_file(outfile + ".filelist_blocks", "w") as blockf:

            for start in range(0, len(infiles), self.block_size):

                fn = outfile + ".filelist_{}".format(start)
                fn_vcf = os.path.join(tmpdir,
                                      "block_{}.vcf.gz".format(start))
                with IOTools.open_file(fn, "w") as outf:
                    end = start + self.block_size
                    outf.write(
                        "\n".join(infiles[start:end]) + "\n")

                statements.append(
                    "{params.path} merge "
                    "{params.options} "
                    "-O v "
                    "--file-list {outfile}.filelist_{start} "
                    "2> {outfile}_merge_{start}.log "
                    "{cleanup_statement} "
                    "{set_genotype} "
                    "| bgzip "
                    "> {fn_vcf}; "
                    "tabix -p vcf {fn_vcf}".format(**locals()))

                blockf.write(fn_vcf + "\n")

        if params.restrict_to_all:
            filter_statement = (
                "| {params.path} filter "
                "--include \"FORMAT/GT != '.'\" "
                "-O v "
                "2> {outfile}_filter.log ".format(**locals()))
        else:
            filter_statement = ""

        statements.append(
            "{params.path} merge "
            "{params.options} "
            "-O v "
            "--file-list {outfile}.filelist_blocks "
            "2> {outfile}_merge.log "
            "{filter_statement} "
            "| bgzip "
            "> {outfile}; "
            "tabix -p vcf {outfile} ".format(**locals()))

        statements.append("rm -rf {}".format(tmpdir))

        statement = "; ".join(statements)

        retvals = P.run(statement,
                        **params._asdict())

        return retvals


class run_collate_vcftools_concat(CollateRunner):
    path = "vcf-concat"
    name = "vcftools_concat"
    path_vcf_sort = "vcf-sort"

    def get_version(self):
        help_string = E.run("{self.path} -h".format(**locals()),
                            return_stderr=True).strip()
        if "vcf-concat [OPTIONS]" in help_string:
            return "unknown"
        else:
            raise ValueError("vcf-concat not found")

    def run(self, infiles, outfile, params):

        files = " ".join(infiles)

        # merge files and remove all variants that have not been called
        # in any of the samples.
        return P.run(
            "{params.path} "
            "{params.options} "
            "{files} "
            "| {params.path_vcf_sort} "
            "| bgzip "
            "> {outfile}; "
            "tabix -p vcf {outfile}; ".format(**locals()))


class run_collate_bcftools_intersect(CollateRunner):
    """run bcftools isec outputting the shared variants

    :param infiles: Input files in :term:`vcf` format
    :param outfile: Output file in :term:`vcf` format

    This tool runs bcftools isec  to output
    the variants common to each of the two files.
    """

    path = "bcftools"
    name = "bcftools_intersect"

    expected = ["vcf"]
    output = ["test_unique.dir/result.vcf.gz",
              "comp_unique.dir/result.vcf.gz",
              "test_shared.dir/result.vcf.gz",
              "comp_shared.dir/result.vcf.gz"]

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]

    def run(self, infiles, outfiles, params):

        vcfs = infiles
        if len(vcfs) != 2:
            raise ValueError("expected 2 VCF files, received {}".format(vcfs))
        vcf = " ".join(infiles)

        if isinstance(outfiles, str):
            # files not known to ruffus, so expect a glob expression such as
            # \2.dir/*.dir/*.vcf.gz
            outdir = os.path.dirname(os.path.dirname(outfiles))
        else:
            outdir = os.path.commonprefix(outfiles)

        outfile = os.path.join(outdir, "result.log")

        retval = P.run(
            "{params.path} isec "
            "{params.options} "
            "--output-type z "
            "--prefix {outdir} "
            "{vcf} "
            "&> {outfile} "
            .format(**locals()))

        f = ["000{}.vcf.gz".format(x) for x in range(4)]
        self.distribute_results(outdir, list(zip(f, self.output)))

        f = ["000{}.vcf.gz.tbi".format(x) for x in range(4)]
        ff = [x + ".tbi" for x in self.output]
        self.distribute_results(outdir, list(zip(f, ff)))

        return retval


class run_collate_bcftools_concat(CollateRunner):
    path = "bcftools"
    name = "bcftools_concat"
    path_tabix = "tabix"

    def get_version(self):
        return "unknown"

    def run(self, infiles, outfile, params):

        files = " ".join(tk.sort_by_chromosome(infiles))

        # merge files and remove all variants that have not been called
        # in any of the samples.
        return P.run(
            "{params.path} concat "
            "{params.options} "
            "{files} "
            "-Oz -o {outfile}; "
            "{path_tabix} -p vcf {outfile}; ".format(**locals()))


class run_collate_vcflib_mergesites(CollateRunner):
    path = "bcftools"
    name = "vcflib_mergesites"
    job_threads = 4

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]

    def run(self, infiles, outfile, params):

        files = " ".join(infiles)

        job_threads = params.job_threads

        # todo:
        # 1. add header.
        # 2. do batch+merge sort in order to avoid hitting temporary space limits.
        # 3. remove unnecessary info fields while sorting, add them later.

        tmpdir = P.get_temp_filename()
        retval = P.run(
            "mkdir {tmpdir}; "
            "bcftools view -h {infiles[0]} "
            "| cut -f 1-10 "
            "| bgzip > {outfile}; "
            "zcat {files} "
            "| awk -v OFS='\\t' "
            "'!/^#/ && $5 != \"<NON_REF>\" "
            "{{$8=\".\";$9=\".\";$6=\".\";$7=\"GT\";$10=\".\"; print}}' "
            "2> {outfile}.filter.log "
            "| sort -k1,1V -k2,2n "
            "--parallel {job_threads} "
            "-T {tmpdir} "
            "2> {outfile}.sort.log "
            "| uniq "
            "| bgzip "
            ">> {outfile}; "
            "tabix -p vcf {outfile}; "
            "rm -rf {tmpdir} ".format(**locals()))


class run_collate_illumina_agg(CollateRunner):
    path = "agg"
    name = "illumina_agg"

    block_size = 500
    regex_filename = "(\S+).dir/result.vcf.gz"

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version:\s*(\S+)", help_string).groups()[0]

    def run(self, infiles, outfile, params):

        phase1_statements = []
        outdir = os.path.dirname(outfile)
        retvals = []
        vcf_files = []
        for fn in infiles:
            prefix = re.search(params.regex_filename, fn).groups()[0]

            out_fn = os.path.join(outdir, "file_{}".format(prefix))
            vcf_files.append(out_fn)

            if os.path.exists(out_fn + ".bcf"):
                continue

            phase1_statements.append(
                "{self.path} ingest1 "
                "--output {out_fn} "
                "--fasta-ref {params.reference_fasta} "
                "{fn} "
                ">& {out_fn}.log; ".format(**locals()))

        phase2_statements = []
        block_files = []
        for start in range(0, len(vcf_files), self.block_size):

            out_fn = os.path.join(outdir,
                                  "block_{}".format(start))

            block_files.append(out_fn)

            if os.path.exists(out_fn + ".bcf"):
                continue

            end = start + self.block_size
            files = " ".join(["{}.bcf".format(x) for x in vcf_files[start:end]])

            phase2_statements.append(
                "{self.path} ingest2 "
                "--output {out_fn} "
                "{files} "
                ">& {out_fn}.log; ".format(**locals()))

        if phase2_statements:
            if phase1_statements:
                retvals.extend(P.run(phase1_statements, job_memory="4G"))
            else:
                E.warn("all files complete for phase 1")
            retvals.extend(P.run(phase2_statements, job_memory="4G"))
        else:
            E.warn("all files complete for phase 2")

        with pysam.VariantFile(block_files[0] + ".bcf") as bcf_file:
            contigs = list(bcf_file.header.contigs)

        files = " ".join(["{}.bcf".format(x) for x in block_files])
        phase3_statements = []
        chromosome_files = []
        for contig in contigs:

            out_fn = os.path.join(outdir, "chr_{}.bcf".format(contig))
            chromosome_files.append(out_fn)
            if os.path.exists(out_fn):
                continue

            phase3_statements.append(
                "{self.path} genotype "
                "--thread 4 "
                "--output-file {out_fn} "
                "--output-type b "
                "-r {contig} "
                "{files} "
                ">& {out_fn}.log; "
                "bcftools index {out_fn}".format(**locals()))

        retvals.extend(P.run(phase3_statements, job_memory="4G", job_threads=4))

        if phase3_statements or not os.path.exists(outfile):
            files = " ".join(chromosome_files)
            retvals.extend(P.run(
                "bcftools concat "
                "-o {outfile} "
                "-O z "
                "{files} "
                ">& {outfile}_concat.log; "
                "tabix -p vcf {outfile}".format(**locals())))

        return retvals


class run_collate_gatk_combine_gvcfs(CollateRunner):

    name = "gatk_combine_gvcfs"

    path = ("/data/install/Licensed/gatk/"
            "GenomeAnalysisTK-2015.1.1-0-g002d5df/GenomeAnalysisTK.jar")

    block_size = 200

    def get_version(self):
        return E.run("java -jar {self.path} --version".format(**locals()),
                     return_stdout=True).strip()

    def run(self, infiles, outfile, params):

        statements = []
        outdir = os.path.dirname(outfile)
        temp_files = []
        for start in range(0, len(infiles), self.block_size):

            fn_vcf = os.path.join(outdir,
                                  "block_{}.vcf.gz".format(start))
            temp_files.append(fn_vcf)

            if os.path.exists(fn_vcf):
                continue

            end = start + self.block_size
            files = " ".join(["--variant {}".format(x)
                              for x in infiles[start:end]])

            statements.append(
                "java "
                "-Djava.io.tmpdir=%(tmpdir)s "
                "-jar {params.path} "
                "-T CombineGVCFs "
                "-R {params.reference_fasta} "
                "{params.options} "
                "{files} "
                "--out {fn_vcf} "
                "--log_to_file {fn_vcf}.log "
                ">& {fn_vcf}.err; ".format(**locals()))

        retvals = P.run(statements, job_memory="28G")
        files = " ".join(["--variant {}".format(x) for x in temp_files])

        statement = (
            "java "
            "-Djava.io.tmpdir=%(tmpdir)s "
            "-jar {params.path} "
            "-T GenotypeGVCFs "
            "-R {params.reference_fasta} "
            "{params.options} "
            "{files} "
            "--out {outfile} "
            "--log_to_file {outfile}.log "
            ">& {outfile}.err; ".format(**locals()))

        retvals.append(P.run(statement, job_memory="28G"))

        return retvals


class run_collate_somatic_variant_detector(CollateRunner):
    """run somatic variant detector.

    In the VCF, the order of samples should be normal, tumour.

    :param infile: Input file in :term:`vcf` format
    :param outfile: Output file in :term:`vcf` format

    """
    name = "somatic_variant_detector"
    path = "somatic_variant_detector"   # weInterpret/ie/tools/somatic_variant_detector

    def get_version(self):
        return "weInterpret 0.8.1"

    def run(self, infiles, outfile, params):

        if isinstance(infiles, list) or isinstance(infiles, tuple):
            if len(infiles) > 1:
                raise NotImplementedError(
                    "collated somatic variant detection of multiple VCF files not implemented")
            infile = infiles[0]
        else:
            infile = infiles

        with pysam.VariantFile(infile) as inf:
            samples = list(inf.header.samples)
            if len(samples) != 2:
                raise ValueError("expected only two samples in VCF, got {}: {}".format(
                    len(samples), ",".join(samples)))
            normal_sample_id, tumour_sample_id = samples

        statement = (
            "{params.path} "
            "{params.options} "
            "-i {infile} "
            "-o {outfile} "
            "-n {normal_sample_id} "
            "-t {tumour_sample_id} "
            "2> {outfile}.log ".format(**locals()))

        return P.run(statement)


class run_collate_daisy_plot_variant_stats(CollateRunner):
    """plot variant statistics.

    :param infile: Input file in :term:`VCF` format
    :param outfile: Output file in :term:`tsv` format

    """

    name = "daisy_plot_variant_stats"
    path = "daisy"

    add_glob = ""

    blob_globs = [os.path.join("*.png")]

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stdout=True).strip()
        return "builtin"

    def run(self, infiles, outfile, params):

        infiles = " ".join([x + params.add_glob for x in infiles])

        statement = (
            "daisy plot-variant-stats "
            "{params.options} "
            "--output-filename-pattern={outfile}.%%s.png "
            "{infiles} "
            "> {outfile}".format(**locals()))
        return P.run(statement)
