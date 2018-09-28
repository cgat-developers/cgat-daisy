import cgatcore.pipeline as P

from daisy.tasks.MetricRunner import MetricRunner


class run_metric_my_tasklibrary_count(MetricRunner):
    name = "my_tasklibrary_count"
    path = "wc"

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):
        return P.run(
            "{params.path} -l < {infile} > {outfile}".format(
                **locals()))
