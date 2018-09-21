import os
import re

from daisy.tasks.Runner import Runner
from daisy.toolkit import as_namedtuple
import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools


class MetricRunner(Runner):

    action = "metric"

    # a list of tables that this metric will produce.
    # If not set by child, the tablename will be derived
    # from the name of the metric.
    tablenames = None

    # a dictionary with index specifications. Each dictionary key is
    # an index name and each specification is a tuple of the table
    # name and the fields that the index should be created
    # on. Multiple fields can be provided as a "," separated string,
    # such as "variant_id,pos"
    tableindices = None

    def __init__(self, *args, **kwargs):
        Runner.__init__(self, *args, **kwargs)
        # initialize table names - shared across class
        if self.tablenames is None:
            self.tablenames = [self.name]

        self.ignore = kwargs.get("ignore", None)

        # make sure that tables are lower case, otherwise issue
        # with mixed-case table names in postgres
        self.tablenames = [x.lower() for x in self.tablenames]

    def __call__(self, infiles, outfile, only_info=False):

        # NOTE: extras not implemented in ruffus 2.6.3, thus
        # use parameter:
        only_info = "only_info" in P.PARAMS

        # ensure output directory exists.
        # This should be done on the pipeline level, but
        # ruffus currently seems not to allow this.
        outdir = os.path.dirname(outfile)
        if outdir and not os.path.exists(outdir):
            os.makedirs(outdir)

        output_files = [self.map_table_to_file(x, outfile)
                        for x in self.tablenames]

        kwargs = {'output_files': output_files,
                  'input_files': infiles,
                  'outdir': outdir}

        if self._runtime_regex:
            kwargs["alias"] = self.build_alias(str(infiles),
                                               regex=self._runtime_regex,
                                               alias=self._runtime_alias)

        self.save_meta(outfile, **kwargs)

        if self.ignore:
            found = False
            for i in self.ignore:
                if i in outdir:
                    found = True
                    break

            if found:
                E.warn("skipping task {} at runtime, an empty file is created".format(
                    outfile))
                IOTools.touch_file(outfile)
                return

        # if self.runtime_filter:
        # TODO: create empty outfile if regex matches
        #    pass

        if only_info:
            E.warn(
                "only_info - meta information in {} has been updated".format(
                    IOTools.snip(outfile) + ".info"))
            return

        # AH: duplicated from above?
        params = self.build_params(output_files=output_files)

        on_error_options = ["raise", "ignore"]
        on_error = params.get("on_error", "raise")
        if on_error not in on_error_options:
            raise ValueError(
                "unknown option to 'on_error': '{}' "
                "should be one of '{}'".format(
                    on_error, ",".join(on_error_options)))

        if self.ignore_task(infiles, outfile, params):
            return

        # deal with placeholder files created by identity that are
        # located on a remote mount point
        def map_to_mount(fn):
            if os.path.exists(fn + ".mnt"):
                if not P.PARAMS["mount_point"]:
                    raise ValueError(
                        "encountered mounted file {}, but no mount point present".format(fn))
                with open(fn + ".mnt") as inf:
                    mount_path = inf.read()
                return os.path.join(P.PARAMS["mount_point"], mount_path)
            else:
                return fn

        # replace infiles with mount locations if necessary
        if isinstance(infiles, list):
            infiles = [map_to_mount(x) for x in infiles]
        else:
            infiles = map_to_mount(infiles)

        try:
            benchmark = self.run(infiles, outfile, as_namedtuple(params))
        except Exception as ex:
            on_error = params.get("on_error", "raise")
            if on_error == "raise":
                raise
            elif on_error == "ignore":
                E.warn("error occured during execution of {} but will be ignored:\n{}".format(
                    self.__name__, ex))
                E.warn("an empty output file {} will be created.".format(
                    outfile))
                IOTools.touch_file(outfile)
                benchmark = None

        if benchmark:
            self.save_benchmark(outfile,
                                benchmark)

    def map_table_to_file(self, tablename, outfile):
        if len(self.tablenames) == 1 and self.tablenames[0] == self.name.lower():
            return outfile
        return outfile + "." + tablename + ".tsv"

    def run(self, outfile):
        raise NotImplementedError(
            "run method needs to be implemented in derived class")


class MetricAgainstReference(MetricRunner):

    def run(self, infile, outfile, params):

        return P.run("{params.path} "
                     "{params.options} "
                     "{infile} "
                     "{params.reference_data} > {outfile}"
                     .format(**locals()))


class run_metric_lengths(MetricAgainstReference):
    name = "lengths"
    path = "daisy compare-lengths"

    def get_version(self):
        return "builtin"


class run_metric_chars(MetricAgainstReference):
    name = "chars"
    path = "daisy compare-chars"

    def get_version(self):
        return "builtin"


class run_metric_filestat(MetricRunner):
    name = "filestat"
    path = "stat"

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stdout=True).strip()
        return re.search("stat \(GNU coreutils\) (\S+)", help_string).groups()[0]

    def run(self, infile, outfile, params):
        return P.run("{params.path} "
                     "--printf='filename\\tsize\\tepoch_modified"
                     "\\tmodified"
                     "\\n"
                     "%%n\\t%%s\\t%%Y\\t%%y\\n' "
                     "{infile} > {outfile}".format(**locals()))


class run_metric_daisy_table2stats(MetricRunner):
    name = "daisy_table2stats"
    path = "daisy table2stats"

    # If glob is given, apply to multiple tables that a tool
    # might have output. The glob expression is appended to the
    # directory of the input filename.
    glob = None

    def get_version(self):
        return "builtin"

    def run(self, infile, outfile, params):

        if params.glob is not None:
            statement = ("{params.path} "
                         "{params.options} "
                         "--input-filename {glob_expression} "
                         "--log {outfile}.log "
                         "2> {outfile}.err "
                         "> {outfile} ".format(
                             glob_expression=os.path.join(
                                 os.path.abspath(os.path.dirname(infile)), params.glob),
                             **locals()))
        else:
            statement = ("{params.path} "
                         "{params.options} "
                         "--input-filename {infile} "
                         "--log {outfile}.log "
                         "2> {outfile}.err "
                         "> {outfile} ".format(**locals()))

        return P.run(statement)
