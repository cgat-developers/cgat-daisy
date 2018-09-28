import os
import pandas
import getpass
import re
import unittest
import tempfile
import shutil
from io import StringIO
import paramiko
from paramiko.ssh_exception import SSHException
from contextlib import contextmanager
from pandas.util.testing import assert_frame_equal
import cgatcore.pipeline as P


class BaseTest(unittest.TestCase):
    def setUp(self):
        self.work_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.work_dir)


class TableTest(BaseTest):

    def check_dataframes_are_equal(self, observed, expected, label="",
                                   sort_column=None):

        if sort_column:
            observed = observed.sort_values(sort_column)
            expected = expected.sort_values(sort_column)

        self.assertEqual(sorted(list(observed.columns)),
                         sorted(list(expected.columns)))

        is_equal = observed.equals(expected)
        unequal_columns = []
        if not is_equal:
            for column in observed.columns:
                if not observed[column].equals(expected[column]):
                    unequal_columns.append(column)

        self.assertTrue(
            is_equal,
            "output unequal for {} (mismatching columns: {})\n"
            "observed:\n{}\nexpected:\n{}".format(
                label,
                unequal_columns,
                observed,
                expected))


class MultipleTablesTest(TableTest):

    def check_files_present(self, output_pattern, expected):

        for x in expected:
            fn = output_pattern % x
            self.assertTrue(
                os.path.exists(fn),
                "expected output file {} is missing".format(fn)
            )

    def check_tables_identical(self, output_pattern, comp_output,
                               sort=False):

        for x, comp_data in comp_output.items():
            data = pandas.read_csv(output_pattern % x, sep="\t")
            comp = pandas.read_csv(StringIO(re.sub(" +", ",", comp_data)))
            if sort:
                data = data.sort_values(list(data.columns)).reset_index(drop=True)
                comp = comp.sort_values(list(comp.columns)).reset_index(drop=True)
            try:
                assert_frame_equal(data, comp)
            except AssertionError as ex:
                raise AssertionError(
                    "{}\nobserved in {}\n{}\nexpected=\n{}\n".format(
                        str(ex),
                        output_pattern % x,
                        data,
                        comp))


@contextmanager
def run_on_cluster(to_cluster):
    if to_cluster:
        P.start_session()
        try:
            yield
        finally:
            P.close_session()
    else:
        yield


def remote_file_exists(filename, hostname=None, expect=False):
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        ssh.connect(hostname, username=getpass.getuser())
    except SSHException as ex:
        # disable test on VM, key issues.
        return expect

    stdin, stdout, ssh_stderr = ssh.exec_command("ls -d {}".format(filename))
    out = stdout.read().decode("utf-8")
    return out.strip() == filename
