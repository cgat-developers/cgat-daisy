#!/usr/bin/env python

"""Tool prototype"""

import sys
import cgatcore.experiment as E


def main(argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    (options, args) = E.start(parser, argv)

    data = "".join(open(args[0]).readlines())

    print(data[::-1])

if __name__ == "__main__":
    sys.exit(main())
