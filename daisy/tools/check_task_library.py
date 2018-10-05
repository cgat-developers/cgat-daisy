"""check-task-library
==========================

Perform a quick check of tools and metrics in the task library.

The script outputs a tab-separated table with three columns:

.. csv-table::
   :header: "column", "description"

   "category", "task category (tool|metric|...)"
   "name", "task name"
   "version", "version string. Set to 'unavailable' if tools in to found"

"""

import sys
import cgatcore.experiment as E

# import tasks to apply in this pipeline
from daisy.tasks import map_tool_to_runner, \
    map_metric_to_runner, \
    map_collate_to_runner, \
    map_split_to_runner


def main(argv=None):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    (options, args) = E.start(parser,
                              argv=argv,
                              add_output_options=True)

    total_counter = E.Counter()
    table = []

    for section, map_task2runner in [
            ("tool", map_tool_to_runner),
            ("metric", map_metric_to_runner),
            ("split", map_split_to_runner),
            ("collate", map_collate_to_runner)]:
        E.debug("processing section: {}".format(section))
        counter = E.Counter()

        for task, taskf in sorted(map_task2runner.items()):
            counter.ntasks += 1
            comments = []
            try:
                version = taskf().get_version()
                counter.version_ok += 1
            except Exception:
                version = ""
                comments.append("unavailable")
                counter.version_fail += 1

            comments = "; ".join(comments)
            table.append(
                (section, task, version,
                 comments))

        E.info("{}: {}".format(section, counter))
        total_counter += counter

    options.stdout.write("section\ttask\tversion\tcomments\n")
    for row in table:
        options.stdout.write("\t".join(map(
            str, row)) + "\n")

    E.info("{}: {}".format("total", counter))
    E.stop()


if __name__ == "__main__":
    sys.exit(main())
