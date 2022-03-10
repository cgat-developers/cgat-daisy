import os
import re

from .Runner import resolve_argument
from .MetricRunner import MetricRunner
import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools


class run_metric_last_dotplot(MetricRunner):
    name = "last_dotplot"

    path_lastal = "lastal"
    path_lastsplit = "last-split"
    path_lastdotplot = "last-dotplot"
    path_mafsort = "maf-sort"
    lastal_options = "-m50 -E0.05"
    lastsplit_options = "-m1"

    min_contig_length = 10000

    blob_globs = [os.path.join("*.png")]

    reference_database = None

    def get_version(self):
        help_string = E.run("{self.path_lastal} -V".format(**locals()),
                            return_stdout=True).strip()
        return re.search("lastal (.+)", help_string).groups()[0]

    def run(self, infile, outfile, params):

        if params.reference_database is None:
            raise ValueError("please provide a reference database")

        statement = (
            "{params.path_lastal} {params.lastal_options} "
            "{params.reference_database} {infile} "
            "| {params.path_lastsplit} {params.lastsplit_options} "
            "| {params.path_mafsort} "
            "| gzip "
            "> {outfile}.maf.gz; "
            "{params.path_lastdotplot} "
            "<(zcat {outfile}.maf.gz "
            "| daisy maf2maf --log={outfile}.filter.log --min-length={params.min_contig_length} ) "
            "{outfile}.png "
            .format(**locals()))

        retval = P.run(statement, job_memory="15G")
        IOTools.touch_file(outfile)
        return retval


class run_metric_quast(MetricRunner):
    name = "quast"
    path = "quast.py"

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stderr=True).strip()
        if help_string and "not found" not in help_string:
            return re.search("QUAST (.+)", help_string).groups()[0]
        else:
            raise ValueError("QUAST not found at/as {}: {}".format(
                self.path, help_string))

    def run(self, infile, outfile, params):

        if "--threads" in params.options:
            job_threads = int(re.search("--threads\s*[=]*(\d+)",
                                        params.options).groups()[0])

        outdir = os.path.dirname(outfile)

        infiles = resolve_argument(infile, sep=" ")

        statement = (
            "{params.path} "
            "{params.options} "
            "-R {params.reference_fasta} "
            "--output-dir {outdir} "
            "{infiles} "
            ">& {outfile}.log2; "
            "perl -p -e 's/# /n/g; s/\(//g; s/\)//g; s/ /_/g' "
            "< {outdir}/transposed_report.tsv "
            "> {outfile}".format(**locals()))

        retval = P.run(statement, job_memory="16G")
        return retval


class run_metric_mummer_dotplot(MetricRunner):
    name = "mummer_dotplot"

    path_nucmer = "nucmer"
    path_dnadiff = "dnadiff"
    path_mummerplot = "mummerplot"

    blob_globs = [os.path.join("*.png")]

    reference_fasta = None

    def get_version(self):
        help_string = E.run("{self.path_nucmer} -v".format(**locals()),
                            return_stderr=True).strip()
        if help_string and "not found" not in help_string:
            return re.search("version (.+)", help_string).groups()[0]
        else:
            raise ValueError("mummer not found at/as {}: {}".format(
                self.path_nucmer, help_string))

    def run(self, infile, outfile, params):

        if params.reference_fasta is None:
            raise ValueError("please provide a reference database")

        statement = (
            "{params.path_nucmer} -p {outfile} {params.reference_fasta} {infile} >& {outfile}.nucmer; "
            "{params.path_dnadiff} -p {outfile} -d {outfile}.delta >& {outfile}.dnadiff; "
            "{params.path_mummerplot} --large --fat --png {outfile}.1delta >& {outfile}.mummerplot"
            .format(**locals()))

        retval = P.run(statement)
        IOTools.touch_file(outfile)
        return retval
