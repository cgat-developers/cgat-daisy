import re

from .Runner import resolve_argument
from .ToolRunner import ToolRunner
import cgatcore.pipeline as P


class run_tool_LDSC(ToolRunner):
    name = "ldsc"
    output = "log"
    path = "/home/steven/ldsc/ldsc.py"
    expected = ["input_sumstats"]

    def get_version(self):
        return "1.0.0"

    def run(self, outfile, params):
        options = params.options
        outdir = re.sub(r".*outdir ", "", options)
        options = re.sub(r"--outdir (.*)", "", options)

        other_sumstats = re.sub(r".*--rg ", "", options)
        if other_sumstats == options:
            other_sumstats = ""
        options = re.sub(r"--rg .*sumstats.gz", "--rg", options)

        outputfile = outfile
        outputfile = re.sub(r"(.*)/", r"\1/" + outdir, outputfile)

        retval = P.run(". ../env/bin/activate; "
                       "{params.path} "
                       "{options} "
                       "{params.input_sumstats}"
                       "{other_sumstats} "
                       "--out {outputfile}"
                       .format(**locals()))

        return retval
