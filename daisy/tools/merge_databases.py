import sys
import re
import os
import glob
import collections
import copy
import itertools
import sqlite3

import cgatcore.experiment as E


def apply_run_filter(run_id, method,
                     min_id, max_id):
    if method is None:
        return False
    elif method == "first":
        return run_id != min_id
    elif method == "last":
        return run_id != max_id


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-d", "--database", dest="databases", action="append",
        help="only show what will be done, don't do it [%default]")

    parser.add_option(
        "-f", "--filter", dest="filter_method", type="choice",
        choices=("first", "last"),
        help="only input a selection of results")

    parser.add_option(
        "-t", "--target", dest="target_database", type="string",
        help="the target database [%default]")

    parser.set_defaults(
        filter_method=None,
        databases=[])

    (options, args) = E.start(parser, argv)

    run_id_offset = 0
    instance_id_offset = 0

    for database in options.databases:
        source_db = sqlite3.connect(database)
        is_instance = False
        is_run = False

        cc = source_db.cursor()
        min_run_id = cc.execute("SELECT MIN (id) FROM run").fetchall()[0][0]
        max_run_id = cc.execute("SELECT MAX (id) FROM run").fetchall()[0][0]
        max_instance_id = cc.execute("SELECT MAX (id) FROM instance").fetchall()[0][0]

        E.info("{}: min_run_id={}, max_run_id={}, max_instance_id={}".format(
            database, min_run_id, max_run_id, max_instance_id))

        for line in source_db.iterdump():

            if line.startswith("CREATE TABLE"):
                try:
                    tablename = re.search("CREATE TABLE \"(\S+)\"", line).groups()[0]
                except AttributeError:
                    tablename = re.search("CREATE TABLE (\S+)", line).groups()[0]

                is_instance = False
                is_run = False
                if tablename == "run":
                    offset = run_id_offset
                    pos = "first"
                    is_run = True
                elif tablename == "tags":
                    offset = run_id_offset
                    pos = "first"
                elif tablename == "instance":
                    is_instance = True
                elif tablename == "tool_timings":
                    offset = instance_id_offset
                    pos = "last"
                elif tablename == "metric_timings":
                    offset = instance_id_offset
                    pos = "last"
                else:
                    # metric table
                    offset = instance_id_offset
                    pos = "last"

            elif line.startswith("INSERT INTO"):

                if is_instance:
                    i, n = re.search("VALUES\((\d+),(\d+),", line).groups()
                    if apply_run_filter(n, options.filter_method, min_run_id, max_run_id):
                        line = None
                    else:
                        line = re.sub("VALUES\({},{},".format(i, n),
                                      "VALUES({},{},".format(
                                          int(i) + instance_id_offset,
                                          int(n) + run_id_offset),
                                      line)
                else:
                    if pos == "last":
                        n = re.search(",(\d+)\)", line).groups()[0]
                        line = re.sub(",{}\)".format(n),
                                      ",{})".format(int(n) + offset),
                                      line)
                    elif pos == "first":
                        n = re.search("VALUES\((\d+),", line).groups()[0]
                        line = re.sub("VALUES\({},".format(n),
                                      "VALUES({},".format(int(n) + offset),
                                      line)
                        if is_run:
                            if apply_run_filter(n, options.filter_method, min_run_id, max_run_id):
                                line = None

            if line is not None:
                print(line)

        cc = source_db.cursor()
        run_id_offset += max_run_id
        instance_id_offset += max_instance_id

        E.info("{}: updated offsets to run_id={}, instance_id={}".format(
            database, run_id_offset, instance_id_offset))

    E.stop()

if __name__ == "__main__":
    sys.exit(main())
