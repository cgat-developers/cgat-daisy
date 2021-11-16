import re
import os

from .ToolRunner import ToolRunner
from .CollateRunner import CollateRunner
from .MetricRunner import MetricRunner

import cgatcore.pipeline as P
import cgatcore.iotools as IOTools
import cgatcore.experiment as E

from .VariantCallers import VariantCaller, get_reference
from .Runner import resolve_argument


class run_tool_ont_nanonet(ToolRunner):

    name = "ont_nanonet"
    expected = ["tgz"]
    output = "result.fastq.gz"
    path = "daisy"

    total_size = 100
    batch_size = 10
    num_threads = 5
    # Note that without specifying a chemistry, nanonet fails
    options = "--chemistry r9.4 --no-event_detect"

    def get_version(self):
        return "builtin"

    def run(self, outfile, params):

        statements = []

        batch_outfiles = []
        for batch_id, batch_start in enumerate(range(
                0, params.total_size, params.batch_size)):

            batch_outfile = "{outfile}_{batch_id}.fastq.gz".format(
                **locals())
            batch_outfiles.append(batch_outfile)

            if os.path.exists(batch_outfile):
                continue

            statements.append(
                "daisy ont-nanonet "
                "--log={batch_outfile}.log "
                "--batch-number={batch_id} "
                "--batch-size={params.batch_size} "
                "--num-threads={params.num_threads} "
                "{params.options} "
                "{params.tgz} "
                "2> {batch_outfile}.err "
                "| gzip "
                "> {batch_outfile}".format(
                    **locals()))

        retvals = P.run(statements,
                        job_array=True,
                        job_threads=params.num_threads)

        batch_outfiles = " ".join(batch_outfiles)
        statement = ("cat {batch_outfiles} > {outfile}; ".format(**locals()))

        #              "rm -f {batch_outfiles}".format(
        #            **locals()))
        retvals.append(P.run(statement))
        return retvals


class run_collate_ont_classify(CollateRunner):
    path = "daisy"
    name = "ont_daisy_classify"

    min_average_quality = 14
    min_length = 1000
    min_size_bytes = 0
    newer_than = None

    def get_version(self):
        return "builtin"

    def run(self, infiles, outfile, params):

        if not outfile.endswith("-pass.fastq.gz"):
            raise ValueError("outfile must end in -pass.fastq.gz, got {}".format(
                outfile))

        if params.min_size_bytes:
            before = len(infiles)
            infiles = [x for x in infiles if os.path.getsize(x) >= params.min_size_bytes]
            E.debug("removing small files: after={}, before={}, removed={}".format(
                    len(infiles), before, before - len(infiles)))

        if params.newer_than:
            before = len(infiles)
            cutoff = os.path.getmtime(params.newer_than)
            infiles = [x for x in infiles if os.path.getmtime(x) > cutoff]
            E.debug("removing old files: after={}, before={}, removed={}".format(
                    len(infiles), before, before - len(infiles)))

        if len(infiles) == 0:
            E.warn("no files left after filtering, creating empty file")
            IOTools.touch_file(outfile)
            return

        infiles = " ".join(infiles)

        outfile_fail = IOTools.snip(outfile, "-pass.fastq.gz") + "-fail.fastq.gz"

        statement = (
            "zcat {infiles} "
            "| daisy fastq2fastq "
            "--method=filter-ONT "
            "--min-average-quality={params.min_average_quality} "
            "--log={outfile}.log "
            "--min-length={params.min_length} "
            "--output-removed-fastq={outfile_fail} "
            "- "
            "| gzip "
            "> {outfile}".format(**locals()))
        return P.run(statement)


class run_metric_ont_classify(MetricRunner):
    path = "daisy"
    name = "ont_daisy_classify"

    min_average_quality = 14
    min_length = 1000

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        outfile_pass = IOTools.snip(outfile, ".tsv") + "-pass.fastq.gz"
        outfile_fail = IOTools.snip(outfile, ".tsv") + "-fail.fastq.gz"

        statement = (
            "zcat {infile} "
            "| daisy fastq2fastq "
            "--method=filter-ONT "
            "--min-average-quality={params.min_average_quality} "
            "--log={outfile}.log "
            "--min-length={params.min_length} "
            "--output-removed-fastq={outfile_fail} "
            "--output-stats-tsv={outfile} "
            "- "
            "| gzip "
            "> {outfile_pass} "
            "".format(**locals()))
        return P.run(statement)


class run_metric_ont_variant_depth_ratio(MetricRunner):
    """use freebayes to return a table with ref and alternate allele
    counts at positions given in a reference VCF file.

    If sample_size is given, a sample of homozygous reference alleles
    is added.
    """

    path = "samtools"

    path_freebayes = "freebayes"
    path_bcftools = "bcftools"

    reference_fasta = None
    reference_vcf = None
    ref_sample_size = None

    name = "ont_variant_depth_ratio"

    options_freebayes = ("--no-indels  --no-mnps --no-complex "
                         "--haplotype-length 0 --pooled-continuous "
                         "--min-alternate-fraction 0")

    options_bcftools = ""

    def get_version(self):
        help_string = E.run("{self.path_freebayes} --version".format(**locals()),
                            return_stdout=True).strip()
        return re.search("version:\s+(\S+)", help_string).groups()[0]

    def run(self, infile, outfile, params):

        if params.reference_fasta is None:
            raise ValueError("ont_variant_depth_ratio requires reference_fasta to be set")

        if params.reference_vcf is None:
            raise ValueError("ont_variant_depth_ratio requires reference_vcf to be set")

        statement = []
        if params.ref_sample_size is not None:
            reference_vcf = outfile + ".ref_sample.vcf.gz"
            statement.append(
                "daisy fasta2vcf "
                "--log={outfile}.fasta2vcf.log "
                "--sample-size={params.ref_sample_size} {params.reference_fasta} "
                "| bgzip "
                "> {outfile}.fasta2vcf.vcf.gz; "
                "tabix -p vcf {outfile}.fasta2vcf.vcf.gz; "
                "bcftools concat --allow-overlap "
                "{params.reference_vcf} "
                "{outfile}.fasta2vcf.vcf.gz "
                "| bgzip "
                "> {reference_vcf}; "
                "tabix -p vcf {reference_vcf} "
                .format(**locals()))
        else:
            reference_vcf = params.reference_vcf

        statement.append(
            "{params.path_freebayes} "
            "-f {params.reference_fasta} "
            "--variant-input {reference_vcf} "
            "--only-use-input-alleles "
            "{params.options_freebayes} "
            "{infile} "
            "| bgzip "
            "> {outfile}.genotyped.vcf.gz; ".format(**locals()))

        # "tabix -p vcf {outfile}.genotyped.vcf.gz; "
        # "{params.path_bcftools} view {params.options_bcftools} "
        # "{reference_vcf} "
        # "| bgzip > {outfile}.ref.vcf.gz; "
        # "tabix -p vcf {outfile}.ref.vcf.gz; "
        # "{params.path_bcftools} query -f \"%%CHROM\\t%%POS\\t[%%GT]\\t[%%DPR]\\n\" "
        # "{outfile}.genotyped.vcf.gz > {outfile}.genotyped.tsv; "
        # "{params.path_bcftools} query -f \"%%CHROM\\t%%POS\\t[%%GT]\\n\" "
        # "{outfile}.ref.tsv; "
        # "join -1 2 -2 2 {outfile}.ref.tsv {outfile}.genotyped.tsv "
        # "| perl -p -e \"s/[, ]/\\t/g\" "
        # "| cut -f 1,3,5,6,7 "
        # "| grep -v '\.' "
        # "> {outfile}".format(**locals()))

        statement = ";".join(statement)
        return P.run(statement)


class run_tool_nanopolish_variants(VariantCaller):
    name = "ont_nanopolish_variants"
    path = "nanopolish"

    fast5_path = "/data/gru/ont/andreas/fast5"

    chunk_size = 1000000

    def get_version(self):
        help_string = E.run("{self.path} variants --version".format(**locals()),
                            return_stdout=True).strip()
        return re.search("Version (\S+)", help_string).groups()[0]

    def run(self, outfile, params):

        bam = resolve_argument(params.bam)
        reference_fasta = get_reference(params)

        filtered_bam = outfile + ".filtered.bam"
        filtered_fasta = outfile + ".filtered.fasta"
        retvals = []

        if not os.path.exists(filtered_bam):
            statement = (
                "samtools view -b -F 256 {bam} > {filtered_bam}; "
                "samtools index {filtered_bam}".format(**locals()))

            retvals.extend(P.run(statement))

        if not os.path.exists(filtered_fasta):
            statement = (
                "samtools view {filtered_bam} "
                "| cut -f 1,10 "
                "| sort "
                "| uniq "
                "| sed -r 's:(\\S+)\\t:>\\1 {params.fast5_path}/\\1.fast5\\n:' "
                "> {filtered_fasta}".format(**locals()))
            retvals.extend(P.run(statement))

        statement = (
            "{params.path} "
            "variants "
            "--reads={filtered_fasta} "
            "--genome={reference_fasta} "
            "--bam={filtered_bam} "
            "{params.options} ".format(**locals()))

        self.run_parallel_by_region(
            statement,
            bam,
            reference_fasta,
            outfile,
            params,
            region_option="--window={contig}:{start}-{end}",
            job_memory="16G")

        return retvals


class run_tools_ont_large_variant_caller(ToolRunner):
    """run

    Note that this script requires samtools to be on the PATH.
    """
    name = "ont_large_variant_caller"
    path = "large-variant-caller"
    expected = ["bam"]
    output = "result.bed.gz"

    max_segment_length = 500000

    def get_version(self):
        return "builtin"

    def run(self, outfile, params):

        bam = resolve_argument(params.bam)

        statement = (
            "daisy ont-large-variant-caller "
            "--log={outfile}.log "
            "--bamfile={bam} "
            "{params.options} "
            "| uniq "
            "| bgzip "
            "> {outfile}.all.bed.gz; "
            "zcat {outfile}.all.bed.gz "
            "| awk '$3 - $2 < {params.max_segment_length}' "
            "| bgzip "
            "> {outfile}; "
            .format(**locals()))

        return P.run(statement)
