"""Test cases for the Pipeline.Execution module."""

import unittest
import os

import cgatcore.experiment as E
import cgatcore.iotools as iotools
from TestUtils import BaseTest


class TestTaskLibrary(BaseTest):

    def setUp(self):
        BaseTest.setUp(self)
        self.library = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "MyTaskLibrary")

    def test_task_library_contains_default_only(self):
        stdout = E.run("daisy run-task --list-tasks -v 0",
                       return_stdout=True)
        tasks = stdout.splitlines()
        self.assertGreater(len(tasks), 0)
        self.assertFalse("tool_my_tasklibrary_cat" in tasks)
        self.assertFalse("metric_my_tasklibrary_count" in tasks)

    def test_task_library_contains_user_tasks(self):
        stdout = E.run(
            "export DAISY_TASKLIBRARY={} && "
            "daisy run-task --list-tasks -v 0".format(self.library),
            return_stdout=True)
        tasks = stdout.splitlines()

        self.assertGreater(len(tasks), 0)
        self.assertTrue("tool_my_tasklibrary_cat" in tasks)
        self.assertTrue("metric_my_tasklibrary_count" in tasks)

    def test_task_library_default_tasks_can_be_executed(self):
        infile = os.path.abspath(__file__)
        outfile = os.path.join(self.work_dir, "out")
        env = os.environ.copy()
        env["DAISY_TASKLIBRARY"] = self.library
        statement = (
            "daisy run-task "
            "--local "
            "--task=metric_filestat "
            "--input-file={} "
            "--output-file={}".format(infile, outfile))

        p = E.run(statement, env=env, return_popen=True)
        stdout, stderr = p.communicate()
        self.assertEqual(p.returncode, 0,
                         msg="stderr = {}".format(stderr))

    def test_task_library_user_tool_can_be_executed(self):
        infile = os.path.abspath(__file__)
        outfile = os.path.join(self.work_dir, "out")

        env = os.environ.copy()
        env["DAISY_TASKLIBRARY"] = self.library

        statement = (
            "daisy run-task "
            "--local "
            "--task=tool_my_tasklibrary_cat "
            "--input-file={} "
            "--input-slot=data "
            "--output-file={}".format(infile, outfile))

        p = E.run(statement, env=env, return_popen=True)
        stdout, stderr = p.communicate()
        self.assertEqual(p.returncode, 0,
                         msg="stderr = {}".format(stderr))

        with iotools.open_file(infile) as inf1, iotools.open_file(outfile) as inf2:
            data1 = inf1.readlines()
            data2 = inf2.readlines()
        self.assertEqual(data1, data2)

    def test_task_library_user_metric_can_be_executed(self):
        infile = os.path.abspath(__file__)
        outfile = os.path.join(self.work_dir, "out")
        env = os.environ.copy()
        env["DAISY_TASKLIBRARY"] = self.library
        statement = (
            "daisy run-task "
            "--local "
            "--task=metric_my_tasklibrary_count "
            "--input-file={} "
            "--input-slot=data "
            "--output-file={}".format(infile, outfile))

        p = E.run(statement, env=env, return_popen=True)
        stdout, stderr = p.communicate()
        self.assertEqual(p.returncode, 0,
                         msg="stderr = {}".format(stderr))

        with iotools.open_file(infile) as inf1, iotools.open_file(outfile) as inf2:
            data1 = inf1.readlines()
            data2 = inf2.readlines()
        self.assertEqual(len(data1), int(data2[0]))


if __name__ == "__main__":
    unittest.main()
