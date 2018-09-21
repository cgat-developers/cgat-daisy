import os
import re

from .ToolRunner import ToolRunner
import cgatcore.pipeline as P
import cgatcore.experiment as E


class run_tool_daisy_finemap(ToolRunner):
    name = "daisy_finemap"
    path = "finemap"

    expected = ["zfile"]
    output = "result.tsv"

    def get_version(self):
        help_string = E.run("{self.path} --help".format(**locals()),
                            return_stdout=True).strip()
        return re.search("Welcome to FINEMAP (\S+)", help_string).groups()[0]

    def run(self, outfile, params):

        prefix = os.path.splitext(os.path.abspath(params.zfile))[0]

        statement = (
            "daisy finemap "
            "--finemap-path={params.path} "
            "--input-prefix={prefix} "
            "--passthrough='{params.options}' "
            "--log={outfile}.log "
            "> {outfile} ".format(**locals()))

        return P.run(statement)
