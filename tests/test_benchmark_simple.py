import unittest
import os
import itertools
import cgatcore.experiment as E
from TestUtils import BaseTest


class TestBenchmarkSimpleBasicWorkflow(BaseTest):

    filename = "benchmark_basic_workflow.yml"

    expected_tools = ["modify", "revert"]
    expected_metrics = ["chars", "lengths"]
    aliases = {}

    def setUp(self):
        BaseTest.setUp(self)

        with open(os.path.join(self.work_dir, "data1.txt"), "w") as outf:
            outf.write("Hello world")

        with open(os.path.join(self.work_dir, "data2.txt"), "w") as outf:
            outf.write("Goodbye world")

    def check_files_present(self, combinations, aliases=None):

        def _check(fn):
            afn = os.path.join(self.work_dir, fn)
            self.assertTrue(
                os.path.exists(afn),
                "expected output files {} is missing".format(afn))

        for tool, metric in combinations:
            tool_aliases = aliases.get(tool, tool)
            if isinstance(tool_aliases, str):
                tool_aliases = [tool_aliases]

            metric_aliases = aliases.get(metric, metric)
            if isinstance(metric_aliases, str):
                metric_aliases = [metric_aliases]

            for tool_alias, metric_alias in itertools.product(tool_aliases,
                                                              metric_aliases):

                tooldir = "{}.dir".format(tool_alias)
                metricdir = "{}.dir".format(metric_alias)

                _check(os.path.join(os.path.join(tooldir, "benchmark.info")))
                _check(os.path.join(os.path.join(tooldir, "benchmark.bench")))
                _check(
                    os.path.join(os.path.join(tooldir, metricdir, "benchmark.info")))
                _check(
                    os.path.join(os.path.join(tooldir, metricdir, "benchmark.bench")))
                _check(
                    os.path.join(os.path.join(tooldir, metricdir, "{}.tsv".format(metric))))
                _check(
                    os.path.join(os.path.join(tooldir, metricdir, "{}.tsv".format(metric))))

    def test_workflow_make(self):

        logfile = os.path.join(self.work_dir, "output.log")
        statement = (
            "daisy "
            "benchmark-simple "
            "--debug "
            "--log={} "
            "--config-file={} "
            "--work-dir={} "
            "--local "
            "make all".format(
                logfile,
                os.path.join(os.path.dirname(__file__), self.filename),
                self.work_dir))

        p = E.run(statement, return_popen=True)
        stdout, stderr = p.communicate()

        if os.path.exists(logfile):
            with open(logfile) as inf:
                logdata = inf.read()
        else:
            logdata = "no log data available at {}".format(logfile)

        self.assertEqual(
            p.returncode, 0, msg="stderr = {}, log={}".format(stderr, logdata))

        self.check_files_present(
            itertools.product(self.expected_tools, self.expected_metrics), aliases=self.aliases)


class TestBenchmarkSimpleBasicWorkflowShow(TestBenchmarkSimpleBasicWorkflow):

    def test_workflow_show(self):

        statement = (
            "daisy benchmark-simple "
            "--debug "
            "--log={} "
            "--config-file={} "
            "--work-dir={} "
            "--local "
            "show all".format(
                os.path.join(self.work_dir, "output.log"),
                os.path.join(os.path.dirname(__file__), self.filename),
                self.work_dir))

        p = E.run(statement, return_popen=True)
        stdout, stderr = p.communicate()

        self.assertEqual(
            p.returncode, 0, msg="stderr = {}".format(stderr))

    def test_workflow_make(self):
        return


class TestBenchmarkToolEmptyOption(TestBenchmarkSimpleBasicWorkflow):

    filename = "benchmark_tool_empty_option.yml"


class TestBenchmarkToolNoOption(TestBenchmarkSimpleBasicWorkflow):

    filename = "benchmark_tool_no_option.yml"


class TestBenchmarkExplicitAlias(TestBenchmarkSimpleBasicWorkflow):

    filename = "benchmark_explicit_alias.yml"

    aliases = {"modify": "modify_with_alias",
               "chars": "chars_with_alias"}


class TestBenchmarkRegexAlias(TestBenchmarkSimpleBasicWorkflow):

    filename = "benchmark_regex_alias.yml"

    aliases = {"modify": ["modify_with_alias_1", "modify_with_alias_2"],
               "chars": ["chars_with_alias_A", "chars_with_alias_B"]}


# class TestBenchmarkExplicitAliasInput(TestBenchmarkSimpleBasicWorkflow):

#     filename = "benchmark_explicit_alias_input.yml"

#     aliases = {"modify": ["modify_with_input_alias1", "modify_with_input_alias2"],
#                "revert": ["revert_with_input_alias1", "revert_with_input_alias2"]}


class TestBenchmarkRegexAliasInput(TestBenchmarkSimpleBasicWorkflow):

    filename = "benchmark_regex_alias_input.yml"

    aliases = {
        "modify": ["modify_with_input_alias1", "modify_with_input_alias2"],
        "revert": ["revert_with_input_alias1", "revert_with_input_alias2"]}


class TestBenchmarkReplication(TestBenchmarkSimpleBasicWorkflow):

    filename = "benchmark_replication.yml"

    aliases = {
        "modify": ["modify_1", "modify_2"],
        "revert": ["revert_1", "revert_2"]}


if __name__ == "__main__":
    unittest.main()
