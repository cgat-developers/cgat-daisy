import re
import pysam

from .MetricRunner import MetricRunner

import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools


class run_metric_daisy_fasta2stats(MetricRunner):
    """run daisy bam2stats tools on a :term:`BAM` file.

    :param infile: Input file in :term:`FASTA` format
    :param outfile: Multiple output files in :term:`tsv`
         format
    """

    name = "daisy_fasta2stats"
    path = "daisy"

    tablenames = ["daisy_fasta2stats_summary",
                  "daisy_fasta2stats_sequences"]

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        statement = (
            "{params.path} fasta2stats "
            "--output-filename-sequences={outfile}.daisy_fasta2stats_sequences.tsv "
            "--log {outfile} "
            "{infile} "
            "> {outfile}.daisy_fasta2stats_summary.tsv "
            .format(**locals()))

        return P.run(statement)
