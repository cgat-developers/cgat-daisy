"""Rename files in a benchmark to fit new configuration
=======================================================

When a benchark configuration changes, the data layout might change
even though the same computation will be performed.

This tool aims to rename directories to fit the new layout.

Use with caution. To apply, simply execute in a project directory::

   daisy benchmark-rename-files

"""

import sys
import os
import glob
import collections
import copy
import itertools

import cgatcore.pipeline as P

import cgatcore.experiment as E
import cgatcore.iotools as IOTools

import daisy.workflow as workflow
import daisy.toolkit as toolkit

# import tasks to apply in this pipeline
from daisy.tasks import map_tool_to_runner, \
    map_metric_to_runner, \
    map_collate_to_runner


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-n", "--dry-run", dest="dry_run", action="store_true",
        help="only show what will be done, don't do it [%default]")

    parser.add_option(
        "-l", "--link", dest="link", action="store_true",
        help="link instead of rename [%default]")

    parser.set_defaults(
        dry_run=False,
        link=False)

    (options, args) = E.start(parser, argv)

    config = P.get_parameters("benchmark.yml")

    old_data, new_data = [], []

    for old_info in glob.glob("*.dir/tool.info"):
        old_dir, old_file = os.path.split(old_info)
        old_info = toolkit.read_data(old_info)
        old_data.append((old_dir, old_info))

    tool_functions = workflow.build_tool_functions(map_tool_to_runner, config)

    config_files = workflow.expand_globs(config["input"])
    input_combos = workflow.build_combinations(config_files)

    map_property_to_dir = collections.defaultdict(list)

    for toolf, input_files in itertools.product(tool_functions, input_combos):

        # create a copy of the task function and give it its unique name
        # by mangling it with the input_files
        taskf = copy.copy(toolf)
        taskf.register_input(input_files)
        result_dir = os.path.basename(os.path.join(taskf.__name__ + ".dir"))
        new_data.append((result_dir, taskf))
        for a, x, y in IOTools.nested_iter(taskf.input_files):
            map_property_to_dir[(x, y)].append(result_dir)
        map_property_to_dir[("name", taskf.name)].append(result_dir)
        for x, y in list(taskf._option_dict.items()):
            map_property_to_dir[(x, y)].append(result_dir)

    # match by input_files
    options.stdout.write("\t".join(("old", "new", "matching")) + "\n")

    for old_dir, old_info in old_data:
        targets = []
        for a, x, y in IOTools.nested_iter(old_info["input_files"]):
            if (x, y) in map_property_to_dir:
                targets.extend(map_property_to_dir[(x, y)])
        for x, y in list(old_info.items()):
            try:
                targets.extend(map_property_to_dir[(x, y)])
            except TypeError:
                pass

        counts = collections.Counter(targets)
        max_count = max(counts.values())
        max_count_items = [x for x, y in list(counts.items()) if y == max_count]

        if len(max_count_items) > 1:
            E.warn("multiple matches for {}, ignored".format(old_dir))
            continue

        new_dir = max_count_items[0]

        options.stdout.write("\t".join(map(str, (
            old_dir, new_dir, max_count))) + "\n")

        if os.path.exists(new_dir):
            raise ValueError("directory {} already exists".format(new_dir))

        if options.dry_run:
            continue

        if options.link:
            os.symlink(old_dir, new_dir)
        else:
            os.rename(old_dir, new_dir)

    E.stop()

if __name__ == "__main__":
    sys.exit(main())
