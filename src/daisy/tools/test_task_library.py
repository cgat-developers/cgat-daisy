"""test-task-library
==========================

Perform a check of tools and metrics implemented in the task library.

To execute all tests, type::

    daisy check-task-library

.. note::

    This tool requires a collection of reference data. See the
    :file:`check_task_library.yml` that sits alongside this script
    for configuration.

"""

import sys
import os
import unittest
import inspect
import yaml
import shutil
import tempfile
import re
import cgatcore.experiment as E
import cgatcore.pipeline as P
import cgatcore.iotools as IOTools

from daisy.toolkit import arvados_enabled

# import tasks to apply in this pipeline
from daisy.tasks import map_tool_to_runner, \
    map_metric_to_runner, \
    map_collate_to_runner, \
    map_split_to_runner, \
    get_task_version


class TestTool(unittest.TestCase):

    def _test_task_has_expected(self, taskf):
        self.assertTrue(len(taskf.expected) > 0)

    def _test_task_has_output(self, taskf):
        self.assertTrue(hasattr(taskf, "output"))
        self.assertTrue(type(taskf.output) in (str, list))

    def _test_task_has_get_version(self, taskf):
        """fail if ValueError is raised with message
        'no version defined for xyz'"""
        version = get_task_version(taskf)

    def _test_task_will_run(self, taskf):

        # check if task is installed
        version = get_task_version(taskf)
        if version is None:
            self.skipTest("tools for task {} not available".format(
                taskf.name))
            return

        # define input/output files
        tool_config = self.test_config["tool"]
        input_files = {}
        for expected in taskf.expected:
            if taskf.name in tool_config:
                # task specific input files
                p = tool_config[taskf.name].get(
                    expected,
                    tool_config.get(expected, None))
            else:
                # generic input files
                p = tool_config.get(expected, None)
            if p is None:
                self.skipTest("data for input slot {} not provided for {}".format(
                    expected, taskf.name))
                return

            input_files[expected] = p

        tmpdir = tempfile.mkdtemp(dir=".",
                                  prefix="tmp_{}".format(taskf.name))
        if isinstance(taskf.output, str):
            outfile = os.path.join(tmpdir, taskf.output)
        else:
            outfile = [os.path.join(tmpdir, x) for x in taskf.output]

        # instantiate task
        t = taskf()
        # set custom options for test
        if taskf.name in tool_config:
            for key, value in tool_config[taskf.name].items():
                setattr(t, key, value)

        # run task
        t.register_input(input_files)
        t(input_files.values(), outfile)

        # check if tool produced non-zero output
        if isinstance(outfile, list):
            for x in outfile:
                self.assertTrue(os.path.exists(x))
                self.assertGreater(os.path.getsize(x), 0)
        else:
            self.assertTrue(os.path.exists(outfile))
            self.assertGreater(os.path.getsize(outfile), 0)

        # cleanup
        try:
            shutil.rmtree(tmpdir)
        except OSError as ex:
            E.warn("could not remove {}: {}".format(tmpdir, ex))


class TestMetric(unittest.TestCase):

    def _test_task_has_get_version(self, taskf):
        """fail if ValueError is raised with message
        'no version defined for xyz'"""
        version = get_task_version(taskf)

    def _test_task_will_run(self, taskf):
        # check if task is installed
        version = get_task_version(taskf)
        if version is None:
            self.skipTest("tools for task {} not available".format(
                taskf.name))
            return

        # define input/output files
        metric_config = self.test_config["metric"]
        infiles = None
        for key, values in metric_config.items():
            if key == taskf.name:
                infiles = values.get("files", None)
            elif "patterns" in values:
                for pattern in values["patterns"]:
                    if re.search(pattern, taskf.name):
                        infiles = values.get("files", None)
                        break
            if infiles:
                break

        if infiles is None:
            self.skipTest("no input files specified for {}".format(taskf.name))
            return

        tmpdir = tempfile.mkdtemp(dir=".",
                                  prefix="tmp_{}_".format(taskf.name))
        outfile = os.path.join(tmpdir, "output.tsv")

        # instantiate task
        task = taskf()
        # set custom options for test
        if taskf.name in metric_config:
            for key, value in metric_config[taskf.name].items():
                setattr(task, key, value)

        # run task
        task(infiles, outfile)

        # check if tool produced non-zero output
        self.assertTrue(os.path.exists(outfile))
        self.assertGreater(os.path.getsize(outfile), 0)

        # cleanup
        try:
            shutil.rmtree(tmpdir)
        except OSError as ex:
            E.warn("could not remove {}: {}".format(tmpdir, ex))


def add_tests(name, taskf, testclass):

    templates = [x for x in dir(testclass) if x.startswith("_test")]
    for template in templates:
        test_name = template[1:]
        full_name = test_name + "_" + name
        # create a function calling the template function
        exec("def {}_f(self): return self.{}(taskf)".format(
            full_name, template), locals())
        # hook function into TestTool
        setattr(testclass, full_name, locals()[full_name + "_f"])


def clear_tests(testclass):
    for name, method in inspect.getmembers(testclass):
        if name.startswith("test_task_"):
            delattr(testclass, name)


def main(argv=None):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-r", "--restrict-regex", dest="restrict_regex", action="append",
        help="pattern to restrict tests to certain tools/metrics. "
        "Can be specified multiple times [%default]")

    parser.add_option(
        "--data-directory", dest="data_directory",
        help="directory with sample data sets. This will override the default "
        "datadir in the configuration file and the environment variable "
        "DAISY_TEST_DATADIR [%default]")

    parser.add_option(
        "--library-directory", dest="library_directory",
        action="append",
        help="directory tasks functions. Will be added to the built-in "
        "and the one specified in DAISY_TASKLIBRARY environment variable "
        "[%default]")

    parser.add_option(
        "--always-mount", dest="always_mount",
        action="store_true",
        help="force mounting of arvados keep [%default]")

    parser.add_option(
        "--keep-failed-temp", dest="keep_failed_temp",
        action="store_true",
        help="keep temporary files of failed tests [%default]")

    parser.set_defaults(
        restrict_regex=[],
        always_mount=False,
        data_directory=None,
        keep_failed_temp=False,
        library_directories=[],
    )

    (options, args) = E.start(parser,
                              argv=argv,
                              add_output_options=True)

    P.get_parameters()

    # load the built-in tests
    filenames = [os.path.join(os.path.dirname(os.path.dirname(__file__)),
                              "tasks", "test_task_library.yml")]
    if "DAISY_TASKLIBRARY" in os.environ:
        filenames.append(os.path.join(os.environ["DAISY_TASKLIBRARY"],
                                      "test_task_library.yml"))
    filenames.extend(options.library_directories)

    master_config = None
    for fn in filenames:
        if not os.path.exists(fn):
            E.warn("file {} does not exist".format(fn))
            continue
        with IOTools.open_file(fn) as inf:
            raw_txt = inf.read()
            test_config = yaml.load(raw_txt, Loader=yaml.FullLoader)
            if test_config is None:
                E.warn("file {} is empty".format(fn))
                continue

            data_directory = os.environ.get(
                "DAISY_TEST_DATADIR",
                test_config.get("data_directory"))

            if options.data_directory:
                data_directory = options.data_directory

            # reload config with placeholders replaced
            test_config = yaml.load(
                re.sub("DATADIR", data_directory, raw_txt),
                Loader=yaml.FullLoader)
            if master_config is None:
                master_config = test_config
            else:
                # add additional tool/test metrics
                master_config["tool"].update(test_config.get("tool", {}))
                master_config["metric"].update(test_config.get("metric", {}))

    for test_section, testclass, map_name_to_runner in [
            ("tool", TestTool, map_tool_to_runner),
            ("metric", TestMetric, map_metric_to_runner)]:

        ignore = master_config[test_section].get("ignore", [])
        # propagate config variables
        testclass.test_config = master_config

        for task, taskf in sorted(map_name_to_runner.items()):
            found = False
            for to_ignore in ignore:
                if re.match(to_ignore, task):
                    found = True
            if found:
                continue
            if options.restrict_regex:
                take = False
                for x in options.restrict_regex:
                    if re.search(x, task):
                        take = True
                if not take:
                    continue
            add_tests(task, taskf, testclass)

    failed = False
    with arvados_enabled(always_mount=options.always_mount):
        for testclass in [TestTool, TestMetric]:
            suite = unittest.TestLoader().loadTestsFromTestCase(
                testclass)
            result = unittest.TextTestRunner(verbosity=2).run(suite)
            failed |= not result.wasSuccessful()

            # remove all tests in test class - necessary if function is
            # called repeatedly
            clear_tests(testclass)

    E.stop()
    return failed


if __name__ == "__main__":
    sys.exit(main())
