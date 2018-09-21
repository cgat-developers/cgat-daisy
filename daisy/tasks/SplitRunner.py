import os

from daisy.tasks.Runner import Runner, collect_file_meta_information
from daisy.toolkit import as_namedtuple
import cgatcore.pipeline as P
import cgatcore.experiment as E


class SplitRunner(Runner):

    action = "split"
    tablenames = None
    _input_alias = None
    _input_regex = None

    def __init__(self, *args, **kwargs):
        Runner.__init__(self, *args, **kwargs)

    def __call__(self, infiles, outfile, only_info=False):

        # NOTE: extras not implemented in ruffus 2.6.3, thus
        # use parameter:
        only_info = "only_info" in P.PARAMS

        self.input_files = collect_file_meta_information({"in": infiles})
        self.input_alias = self.build_alias(str(infiles),
                                            regex=self._input_regex,
                                            alias=self._input_alias)

        is_empty_outfile = outfile == []
        # patch for missing outfiles when number of outfiles is not
        # known (is it a ruffus thing?)
        if is_empty_outfile:
            assert isinstance(infiles, str)
            outdir = "{}/{}.dir".format(
                os.path.dirname(infiles),
                self.name)
            outfile = os.path.join(outdir, self.output)
            # dummy.info gets removed and then added
            self.save_meta(os.path.join(outdir, "dummy.info"))
        else:
            self.save_meta(outfile)

        if only_info:
            E.warn(
                "only_info - meta information in {} has been updated".format(
                    self.build_meta_filename(outfile, "benchmark.info")))
            return

        params = self.build_params()
        benchmark = self.run(infiles, outfile, as_namedtuple(params))

        if not is_empty_outfile:
            self.save_benchmark(
                outfile,
                benchmark)

    def run(self, outfile):
        raise NotImplementedError(
            "run method needs to be implemented in derived class")
