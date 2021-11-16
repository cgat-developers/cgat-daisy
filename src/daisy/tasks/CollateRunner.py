from daisy.tasks.Runner import Runner, collect_file_meta_information
from daisy.toolkit import as_namedtuple
import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools
import os
import glob
import re


class CollateRunner(Runner):

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

        self.set_replication_id_from_regex(infiles)

        if isinstance(outfile, list):
            outdir = [os.path.dirname(x) for x in outfile]
            basefile = os.path.commonprefix(outfile)
        else:
            outdir = os.path.dirname(outfile)
            basefile = outfile

        kwargs = {'output_file': outfile,
                  'input_files': infiles,
                  'outdir': outdir}

        self.save_meta(outfile, **kwargs)

        if only_info:
            E.warn(
                "only_info - meta information in {} has been updated".format(
                    os.path.join(os.path.dirname(basefile), "benchmark.info")))
            return

        params = self.build_params(output_file=outfile)
        benchmark = self.run(infiles, outfile, as_namedtuple(params))

        self.save_benchmark(
            basefile,
            benchmark)

    def run(self, outfile):
        raise NotImplementedError(
            "run method needs to be implemented in derived class")


class run_collate_link_output(CollateRunner):

    path = ""
    name = "link_output"

    suffixes = [".tbi", ".bai"]

    regex = "([^/]+).dir/[^/]+\.(\S+)"
    pattern_out = "\\1.\\2"

    def get_version(self):
        return "builtin"

    def run(self, infiles, outfile, params):

        def _link(infile, outfile):
            if os.path.exists(os.path.abspath(outfile)):
                return

            dirname = os.path.dirname(outfile)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
            os.symlink(infile, os.path.abspath(outfile))

        rx = re.compile(params.regex)

        outfiles = []
        for infile in infiles:

            outpath = os.path.join(
                os.path.dirname(outfile),
                rx.search(infile).expand(params.pattern_out))

            for suffix in self.suffixes:
                for fn in glob.glob(infile + suffix):
                    _link(fn, outpath + suffix)
            _link(os.path.abspath(infile), outpath)
            outfiles.append(outpath)

        with IOTools.open_file(outfile, "w") as outf:
            outf.write("\n".join(outfiles) + "\n")


class run_collate_cat(CollateRunner):

    path = ""
    name = "cat"

    regex = "([^/]+).dir/[^/]+\.(\S+)"
    pattern_out = "\\1.\\2"

    def get_version(self):
        return "builtin"

    def run(self, infiles, outfile, params):

        infiles = " ".join(infiles)
        statement = (
            "cat {infiles} > {outfile}".format(**locals()))
        return P.run(statement)
