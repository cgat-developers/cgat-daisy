#!/usr/bin/env python

"""Prototype of a metric collector"""

import sys
import cgatcore.experiment as E
import cgatcore.iotools as IOTools


def main(argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    (options, args) = E.start(parser, argv)

    with IOTools.open_file(args[0]) as inf:
        data = "".join(inf.readlines()).strip()
    with IOTools.open_file(args[1]) as inf:
        reference = "".join(inf.readlines()).strip()

    options.stdout.write("key\tinput\treference\n")
    options.stdout.write("length\t{}\t{}\n".format(len(data), len(reference)))

    E.stop()

if __name__ == "__main__":
    sys.exit(main())
