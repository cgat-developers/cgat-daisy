import pysam
import os

from .ToolRunner import ToolRunner
import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools
from .Runner import resolve_argument


class run_tool_pbsim(ToolRunner):
    """run the pbsim simulator to simulate sequenes.

    The simulator expects two genomic sequences corresponding
    to two haplotypes.

    This simulator samples read lengths and quality scores
    from a fastq file and outputs a file called "result.bam"
    with the aligned reads and a file "result.fastq.gz" with
    the read sequences.
    """

    path = "pbsim"
    name = "pbsim_fastq"

    expected = ["reference_fasta", "fastq"]
    output = "result.bam"
    set_quality_score = False

    use_sample_method = True

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        if "USAGE: pbsim" in help_string:
            return "unknown"
        else:
            raise ValueError("pbsim not found")

    def run(self, outfile, params):

        reference_fasta = resolve_argument(params.reference_fasta, ",").split(",")
        if len(reference_fasta) == 2:
            reference1, reference2 = reference_fasta
        else:
            raise NotImplementedError()

        outfile = os.path.abspath(outfile)
        outfile_fastq = os.path.join(os.path.dirname(outfile), "result.fastq.gz")

        # build BAM header
        with pysam.FastaFile(reference1) as inf:
            with IOTools.open_file(outfile + ".header.sam", "w") as outf:
                outf.write("@HD\tVN:1.3\tSO:unsorted\n")
                for contig, length in zip(inf.references, inf.lengths):
                    outf.write("@SQ\tSN:{}\tLN:{}\n".format(contig, length))

        # not enough space on tmp
        # tmpdir = P.get_temp_filename(clear=True)
        tmpdir = os.path.join(os.path.dirname(outfile), "tmp")
        statements = []
        statements.append("mkdir -p {tmpdir}")

        if params.use_sample_method:
            fastq_filename = os.path.join(tmpdir, "tmp.fastq")
            if not os.path.exists(fastq_filename):
                if params.set_quality_score and params.set_quality_score.strip():

                    statements.append(
                        "daisy fastq2fastq "
                        "--quality-offset={params.set_quality_score} "
                        "--log={outfile}.fastq.log "
                        "{params.fastq} "
                        "> {fastq_filename}".format(
                            **locals()))
                else:
                    statements.append(
                        "zcat {params.fastq} > {tmpdir}/tmp.fastq")

            statements.append("cd {tmpdir}")

            statements.append(
                "{params.path} "
                "{params.options} "
                "--sample-fastq={tmpdir}/tmp.fastq "
                "{reference1} "
                "--prefix=H1 "
                ">& {outfile}.pbsim1.log")

            statements.append(
                "{params.path} "
                "{params.options} "
                "--sample-fastq={tmpdir}/tmp.fastq "
                "{reference2} "
                "--prefix=H2 "
                ">& {outfile}.pbsim2.log")
        else:
            statements.append("cd {tmpdir}")

            statements.append(
                "{params.path} "
                "{params.options} "
                "--prefix=H1 "
                "{reference1} "
                ">& {outfile}.pbsim1.log")
            statements.append(
                "{params.path} "
                "{params.options} "
                "--prefix=H2 "
                "{reference2} "
                ">& {outfile}.pbsim2.log")

        statements.append(
            "daisy fastq2fastq "
            "--input-fastq-file={tmpdir}/H1_0001.fastq "
            "--output-removed-tsv={outfile}.removed1 "
            "--set-prefix=H1 "
            "--method=filter-N "
            "--log={outfile}.log "
            "| gzip "
            "> {outfile_fastq}")

        statements.append(
            "cat {tmpdir}/H1_0001.maf "
            "| daisy maf2maf "
            "--input-filter-tsv={outfile}.removed1 "
            "--log={outfile}.maf.log "
            "--set-prefix=H1 "
            "> {tmpdir}/tmp.maf")

        statements.append(
            "daisy fastq2fastq "
            "--input-fastq-file={tmpdir}/H2_0001.fastq "
            "--output-removed-tsv={outfile}.removed2 "
            "--method=filter-N "
            "--set-prefix=H2 "
            "--log={outfile}.log "
            "| gzip "
            ">> {outfile_fastq}")

        statements.append(
            "cat {tmpdir}/H2_0001.maf "
            "| daisy maf2maf "
            "--input-filter-tsv={outfile}.removed1 "
            "--log={outfile}.maf.log "
            "--set-prefix=H2 "
            ">> {tmpdir}/tmp.maf")

        # generalize for chromosomes
        statements.append(
            "maf-convert sam {tmpdir}/tmp.maf "
            "| grep -v '^@' "
            "| perl -p -e \"s/ref/22/\" "
            ">> {tmpdir}/tmp.sam")

        statements.append(
            "cat {outfile}.header.sam {tmpdir}/tmp.sam "
            "| samtools view -bS "
            "| samtools sort -T {tmpdir}/ -O bam - > {outfile}")

        statements.append(
            "samtools index {outfile}")

        statements.append(
            "rm -rf {tmpdir}")

        statement = "; ".join(statements).format(**locals())

        return P.run(statement)
