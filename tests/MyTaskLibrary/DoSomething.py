import cgatcore.pipeline as P
from daisy.tasks.ToolRunner import ToolRunner


class run_tool_my_tasklibrary_cat(ToolRunner):
    name = "my_tasklibrary_cat"
    path = "cat"
    expected = ["data"]

    def get_version(self):
        return "builtin"

    def run(self, outfile, params):
        return P.run("{params.path} {params.data} > {outfile}".format(
                     **locals()))
