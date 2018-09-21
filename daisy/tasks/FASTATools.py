from .ToolRunner import ToolRunner

import cgatcore.pipeline as P
import cgatcore.experiment as E


class run_tool_daisy_fasta2fasta(ToolRunner):
    """run daisy fasta2fasta tools on a :term:`FASTA` file.

    :param infile: Input file in :term:`FASTA` format
    :param outfile: Output file in :term:`FASTA` format
    """

    name = "daisy_fasta2fasta"
    path = "daisy"

    expected = ["fasta"]
    output = "result.fasta"

    methods = []

    def get_version(self):
        return "builtin"

    def run(self, outfile, params):

        if not params.methods:
            raise ValueError("fasta2fasta requires one or more methods")

        if isinstance(params.methods, list):
            methods = ",".join(params.methods)
        else:
            methods = params.methods

        statement = (
            "{params.path} fasta2fasta "
            "--log={outfile}.log "
            "--method={methods} "
            "{params.options} "
            "{params.fasta} "
            "> {outfile}; "
            "samtools faidx {outfile} "
            .format(**locals()))

        return P.run(statement)
