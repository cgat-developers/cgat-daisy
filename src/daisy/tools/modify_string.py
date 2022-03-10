"""Prototype of a tool"""

import sys
import re
import cgatcore.experiment as E


def main(argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-o", "--option", dest="option", type="string")

    (options, args) = E.start(parser, argv)

    data = "".join(open(args[0]).readlines())

    print(re.sub("o", "a", data))

    E.stop()

if __name__ == "__main__":
    sys.exit(main())
