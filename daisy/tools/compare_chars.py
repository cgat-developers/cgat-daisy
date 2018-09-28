#!/usr/bin/env python

"""Prototype of a metric collector"""

import sys
from collections import Counter
import cgatcore.experiment as E
import cgatcore.iotools as IOTools


def main(argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-o", "--option", dest="option", type="string")

    (options, args) = E.start(parser, argv)

    with IOTools.open_file(args[0]) as inf:
        data = "".join(inf.readlines()).strip()
    with IOTools.open_file(args[1]) as inf:
        reference = "".join(inf.readlines()).strip()

    data_counts = Counter(data)
    ref_counts = Counter(reference)

    keys = set(list(data_counts.keys()) + list(ref_counts.keys()))

    options.stdout.write("key\tinput\treference\n")
    for key in sorted(keys):
        options.stdout.write(
            "\t".join((key, str(data_counts[key]), str(ref_counts[key]))) + "\n")

    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
