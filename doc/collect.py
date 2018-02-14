"""Collect daisy tools into docs
================================

To use, type:

   python collect.py

"""

import sys
import glob
import os
import re
import daisy.Experiment as E
import daisy.IOTools as IOTools

TEMPLATE_TOOL = """
========================================================
{tool_name}
========================================================

.. automodule:: daisy.tools.{tool_module}

Command line usage
==================

.. program-output:: daisy {tool_name} -?
"""


def main(argv=None):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    (options, args) = E.start(parser,
                              argv=argv, add_output_options=True)

    tools = glob.glob(
        os.path.join(os.path.dirname(__file__),
                     "..", "src", "daisy", "tools", "*.py"))

    counter = E.Counter()
    for tool in tools:
        counter.found += 1
        tool_module = re.sub(".py", "", os.path.basename(tool))
        tool_name = re.sub("_", "-", tool_module)
        if tool_name in ("__init__", "cli"):
            c.ignored += 1
            continue

        dest = os.path.join("tools", "{}.rst".format(tool_name))

        if os.path.exists(dest) and not options.output_force:
            counter.skipped += 1
            continue

        with IOTools.openFile(dest, "w") as outf:
            outf.write(TEMPLATE_TOOL.format(
                **locals()))

        counter.new += 1

    E.info(counter)
    E.stop()


if __name__ == "__main__":
    sys.exit(main())
