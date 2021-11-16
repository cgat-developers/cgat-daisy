import inspect
import json
import os
import re
import sre_constants
import collections
import shutil

import cgatcore.experiment as E
import cgatcore.pipeline as P
import cgatcore.iotools as IOTools
from daisy.toolkit import hash
from daisy.tasks import get_default_params


def resolve_argument(argument, sep=","):
    """if argument is a container type (dict, list, tuple)
    resolve its contents to comma-separated list.
    """
    if isinstance(argument, dict):
        if len(argument) != 1:
            raise ValueError(
                "expected a single entry dictionary, got '{}'"
                .format(argument))
        return sep.join(x[2] for x in IOTools.nested_iter(argument))
    elif isinstance(argument, list) or isinstance(argument, tuple):
        return sep.join(argument)
    # special treatment for output from run_collate_link_output
    elif "filelist" in argument:
        f = [x.strip() for x in IOTools.open_file(argument).readlines() if not x.startswith("#")]
        return sep.join([x for x in f if x])

    return argument


def dict_to_labelvalue(d):
    """return label, value of a single entry dictionary.

    If the dictionary contains more than one key/value pair,
    return the complete dictionary.
    """
    if len(d) == 1:
        return list(d.keys())[0], list(d.values())[0]
    else:
        return "all", d


def compare_sets(a, b):
    """return missing and additional elements in two containers.
    """
    expected = set(a)
    received = set(b)
    return expected.difference(received), received.difference(expected)


def is_true(param):
    """return True if a flag has been set."""
    return param and param not in ("0", "F")


def is_set(param):
    """return True if a value has been set."""
    if param is None or param.strip() == "":
        return False
    return True


def always_succeed(f):
    def inner(self, outfile, *args, **kwargs):
        try:
            f()
        except Exception as e:
            E.warn("received exception {} - touching {}".format(str(e),
                                                                outfile))
        IOTools.touch_file(outfile)
    return inner


def collect_file_meta_information(file_dict, nostats=False):
    """collect meta information on files

    Arg:
       file_dict(dict) : nested dictionary

    Returns:
       info(list)
    """
    results = []

    for d, key, filenames in IOTools.nested_iter(file_dict):
        if filenames is None:
            continue

        if isinstance(filenames, str):
            filenames = filenames.split(",")

        filenames = [x.strip() for x in filenames]
        for filename in filenames:
            abspath = os.path.realpath(filename)
            if nostats:
                st_size, st_mtime, st_ctime = 0, 0, 0
            else:
                if not os.path.exists(abspath):
                    raise OSError("file {} does not exist".format(filename))
                s = os.stat(filename)
                st_size, st_mtime, st_ctime = s.st_size, s.st_mtime, s.st_ctime

            results.append(
                collections.OrderedDict(list(zip(
                    ("path",
                     "abspath",
                     "size",
                     "modification_time",
                     "creation_time"),
                    (filename,
                     abspath,
                     st_size,
                     st_mtime,
                     st_ctime)))))

    return results


def get_associated_file(params, reference_label, option_field):
    """
    :param params: parameters for task.
    :param reference_label: lookup to perform for in option_field.
    :param option_field: the name of the option as defined in
        the parameter class.
    """
    option_value = getattr(params, option_field)
    option_map = option_field + "_map"
    option_defaults = "default_" + option_field

    if option_value.startswith("auto"):
        if reference_label:
            values_dict = getattr(params, option_map).get(
                option_value, None)

            if values_dict is None:
                raise KeyError(
                    "can not find '{}' section '{}', "
                    "expected one of {}".format(
                        option_field,
                        option_value,
                        str(sorted(getattr(params, option_defaults).keys()))))

            if reference_label not in values_dict:
                raise KeyError(
                    "can not find regions for reference '{}', "
                    "expected one of {}".format(
                        reference_label,
                        str(sorted(values_dict.keys()))))
            return values_dict[reference_label]
        else:
            raise ValueError(
                "auto detection of regions requires auto-detection "
                "of reference sequence (set 'reference_fasta: auto') "
                "and successful "
                "detection of the reference sequence")
    else:
        # simple pass-through
        return option_value


# A version cache to prevent multiple checks form the same task
VERSION_CACHE = {}


# Be careful with state with these functors, the same object
# will be called multiple times in a pipeline.
class Runner(object):

    # The name identifies a section in the configuration file
    name = "Runner"

    # Tool command line options
    options = ""

    def __init__(self, **kwargs):

        # list of files that this tool depends upon, for example
        # configuration files.
        self.dependencies = []

        self._regex = None
        self._alias = None
        self._name = None
        self._explicit_alias = []
        self._runtime_alias = None
        self._runtime_regex = None
        self._replicate_regex = None
        self._ignore = None
        self._replication_id = None

        # add configuration values
        for key, value in list(kwargs.items()):
            # skip empty values as they cause problems
            # if they are later used in string interpolation
            # and a None appears instead of ""
            if value is None:
                continue
            if key in ("regex", "alias", "name",
                       "runtime_regex", "runtime_alias",
                       "ignore", "replicate_regex"):
                # add private functions here. Prevents for example
                # regex matching itself.
                setattr(self, "_" + key, value)
                del kwargs[key]
            else:
                if isinstance(value, str) and "alias=" in value:
                    x = re.search("alias=([^;]+);", value)
                    self._explicit_alias.append(x.groups()[0])
                    value = " ".join((value[:x.start()].strip(),
                                      value[x.end():].strip())).strip()
                setattr(self, key, value)

        self.option_hash = hash(str(sorted(kwargs.items())))
        self._option_dict = kwargs
        self.setup()

        # create a unique name by adding the parameterization
        if self._name:
            self.__name__ = self._name
        else:
            self.__name__ = self.name + "_" + self.option_hash

    def set_mountpoint(self, mountpoint):
        self.mountpoint = mountpoint

    def build_alias(self, option_string, regex=None, alias=None):
        name = ""

        # remove quotes and other characters to simplify writing
        # regex expressions
        option_string = re.sub("[\"':{}]", " ", str(option_string))

        if regex:
            name = re.search(regex, option_string)

            if name:
                if alias:
                    name = name.expand(alias)
                else:
                    name = "_".join(name.groups())
            else:
                raise ValueError(
                    "regex '{}' did not match options '{}'".format(
                        regex, option_string))

        return name

    def get_plain_name(self, optional_data=None):

        parts = []
        if self._explicit_alias:
            parts.extend(self._explicit_alias)

        if self._regex:
            option_string = [str(x) for x in list(self._option_dict.values())
                             if x is not None]
            if optional_data:
                option_string.append(str(optional_data))
            option_string = " ".join(option_string)
            parts.append(self.build_alias(option_string,
                                          alias=self._alias,
                                          regex=self._regex))
        return "_".join(parts)

    def set_replication_id(self, replication_id):
        self._replication_id = replication_id

    def set_replication_id_from_regex(self, infiles):
        if self._replicate_regex:
            try:
                replication_ids = [re.search(self._replicate_regex, x).groups()[0] for x in infiles]
            except AttributeError:
                raise ValueError("replication id could not be extracted from {} with regex {}".format(
                    infiles, self._replicate_regex))
            replication_ids = list(set(replication_ids))
            if len(replication_ids) > 1:
                raise ValueError("received multiple replication ids from {} with regex {}".format(
                    infiles, self._replicate_regex))
            try:
                self._replication_id = int(replication_ids[0])
            except ValueError:
                raise ValueError("non-numerical replication id: {}".format(replication_ids[0]))
        else:
            self._replication_id = 1

    def register_input(self,
                       input_files=None,
                       alias=None,
                       regex=None,
                       make_unique=False,
                       is_test=False,
                       replicate_alias=None):
        """notify task of input files.

        Given a dictionary mapping slots to input files check that
        each input files exists, collect meta information and build
        input aliases.

        """

        if self._name:
            parts = [self._name]
        else:
            parts = [self.name]

        if self._replication_id:
            if replicate_alias:
                if "\\1" not in replicate_alias:
                    raise ValueError("replicate_alias must contain a '\\1', but is '{}'".format(
                        replicate_alias))
                parts.append(re.sub("\\\\1", str(self._replication_id), replicate_alias))
            else:
                parts.append(str(self._replication_id))

        # build a unique name with human readable components
        plain_name = self.get_plain_name()
        if plain_name:
            parts.append(plain_name)

        if make_unique:
            parts.append(self.option_hash)

        self.input_alias = None

        if input_files is not None:

            if "name" in input_files:
                # use explicitely given name
                self.input_alias = input_files["name"]

                input_files = dict([(x, y) for x, y in list(input_files.items())
                                    if x != "name"])

            self.input_labels = list(input_files.keys())
            if is_test:
                self.input_files = collect_file_meta_information(
                    input_files, nostats=True)
            else:
                self.input_files = collect_file_meta_information(input_files)

            if self.input_alias is None:
                # derive name from input files and regular expressions
                # that have been given.

                # build input alias from input files dictionary using regex
                try:
                    self.input_alias = self.build_alias(str(input_files),
                                                        regex=regex,
                                                        alias=alias)
                except ValueError:
                    if is_test:
                        self.input_alias = ""
                    else:
                        raise

                # fall-back on group labels
                if self.input_alias == "":
                    group_labels = []
                    for k, v in list(input_files.items()):
                        if isinstance(v, dict):
                            group_labels.extend(list(v.keys()))
                    self.input_alias = "".join(group_labels)

            if self.input_alias:
                parts.append(self.input_alias)

            if make_unique:
                parts.append(hash(str(sorted(
                    [x["abspath"] for x in self.input_files]))))

        self.__name__ = "_".join(parts)

    def validate_input(self):

        missing, extra = compare_sets(
            set(self.expected),
            set(self.input_labels))

        if missing:
            raise ValueError("missing input for {}: {}".format(
                self.name, list(missing)))

    def instantiate_input(self, input_files):

        self.validate_input()
        for l, v in zip(self.input_labels, input_files):
            setattr(self, l, v)

    def build_params(self, **kwargs):
        """build parameter dictionary.

        This method applies the following tasks to the parameter
        dictionary:

        1. add default parameters from global config file

        2. interpolate task specific parameters.

        3. for any <param>, check if "<param>_regex exists. If
           it does, apply regular expression on output_files and
           interpolate with the value of <param>. <param>_regex needs to
           contain one or more grouping statements and <param> needs
           to contain the appropriate number of backreferences (\1,
           \2).

        """

        # insert defaults from function body
        d = inspect.getmembers(self, lambda a: not(inspect.isroutine(a)))
        d = dict([x for x in d if not x[0].startswith("_")])

        # collect version
        key = self.__class__.__name__
        version = VERSION_CACHE.get(key, self.get_version())
        d["version"] = version
        VERSION_CACHE[key] = version

        default_params = get_default_params()
        if self.name in default_params:
            d.update(default_params[self.name])

        # update with specifically given paramters.
        d.update(kwargs)

        if "task_specific" in d and ("output_files" in d or
                                     "output_file" in d):
            files = d.get("output_files", []) + [d.get("output_file", None)]
            for f in [x for x in files if x is not None]:
                # sort so that more general (shorter) keys come first
                for key in sorted(d["task_specific"].keys()):
                    if re.search(key, f):
                        for option, value in list(d["task_specific"][key].items()):
                            d[option] = value

        if "output_file" in d or "output_files" in d:
            if "output_files" in d:
                to_match = d["output_files"]
            else:
                to_match = [d["output_file"]]
            for key, value in list(d.items()):
                if key + "_regex" in d:
                    try:
                        rx = re.compile(d[key + "_regex"])
                    except sre_constants.error as e:
                        raise ValueError(
                            "error in creating regular expression for {} '{}': {}".format(
                                key, d[key + "_regex"], str(e)))
                    for o in to_match:
                        # use the first output file that matches
                        r = rx.search(o)
                        if r:
                            d[key] = r.expand(value)
                            break
                    else:
                        raise ValueError(
                            "regular expresssion '{}' did not match any parameters in {}".format(
                                d[key + "_regex"],
                                ",".join(["'{}'".format(x) for x in to_match])))
        return d

    def build_meta_filename(self, outfile, name):

        if isinstance(outfile, list):
            # multiple output files end up in different subdirs, take
            # the root of these
            dirname = os.path.commonprefix(outfile)
        else:
            dirname = os.path.dirname(outfile)

        return os.path.join(dirname, name)

    def save_meta(self, outfile, **kwargs):

        d = self.build_params(**kwargs)

        if "name" not in d:
            raise ValueError("tool/metric without name: {}".format(d))

        if "version" not in d:
            raise ValueError("tool {} has no version: {}".format(
                d["name"], d))

        if "alias" not in d:
            # submit d in order to look for input files.
            d["alias"] = self.get_plain_name(d)

        if self._replication_id:
            d["replication_id"] = self._replication_id
        else:
            d["replication_id"] = 1

        filename = self.build_meta_filename(outfile, "benchmark.info")
        if not os.path.exists(os.path.dirname(filename)):
            os.makedirs(os.path.dirname(filename))

        with open(filename, "w") as outf:
            outf.write(json.dumps(d, sort_keys=True))

    def save_benchmark(self, outfile, benchmark):

        if not isinstance(benchmark, list):
            benchmark = [benchmark]

        # flatten if nested list and remove None
        benchmark = [x for x in IOTools.flatten(benchmark, ltypes=(list,))
                     if x is not None]

        filename = self.build_meta_filename(outfile, "benchmark.bench")

        if not benchmark:
            E.warn("could not save benchmark info to {}".format(filename))
            return

        try:
            header = benchmark[0]._fields
        except AttributeError as ex:
            E.warn("could not save benchmark timings for {}:"
                   " {} from {}".format(outfile,
                                        str(ex),
                                        str(benchmark[0])))
            return
        
        header_line = "\t".join(header) + "\n"
        # append to existing file
        if os.path.exists(filename):
            open_mode = "a"
            first_line = IOTools.get_first_line(filename)
            if first_line != header_line:
                raise ValueError(
                    "header mismatch when appending to existing file {}: got '{}', expected '{}'".format(
                        filename, first_line, "\t".join(header_line)))
        else:
            open_mode = "w"

        with open(filename, open_mode) as outf:
            if open_mode == "w":
                outf.write(header_line)
            for b in benchmark:
                outf.write("\t".join(map(str, b)) + "\n")

    def ignore_task(self, infiles, outfiles, params):
        """return True if task should be ignored.

        This method will also create the output file(s).
        """
        if self._ignore:
            m = str(outfiles)
            for ignore in IOTools.val2list(self._ignore):
                if ignore in m:
                    E.warn(
                        "task {} will be ignored".format(self.__name__))
                    for f in IOTools.val2list(outfiles):
                        E.info("creating empty file {}".format(f))
                        IOTools.touch_file(f)
                    return True
        return False

    def distribute_results(self, workdir, pairs, statement=None):
        """distribute results from a task into separate output directories.

        Arguments
        ---------
        workdir : string
            working directory
        pairs : list
            tuples of input/output filenames
        statement : string
            optional statement to be executed to transform input
            to output. If not given, the files are simply moved.

        """
        statements = []
        for infile, outfile in pairs:
            infile = os.path.join(workdir, infile)
            outfile = os.path.join(workdir, outfile)
            if not os.path.exists(infile):
                raise ValueError(
                    "expected file {} does not exist".format(infile))
            if not os.path.exists(os.path.dirname(outfile)):
                os.makedirs(os.path.dirname(outfile))
            if statement is None:
                shutil.move(infile, outfile)
            else:
                statements.append(statement.format(**locals()))
        if statements:
            return P.run(statements)

    def setup(self):
        pass


class EmptyRunner(Runner):

    def __call__(self, *args, **kwargs):
        pass
