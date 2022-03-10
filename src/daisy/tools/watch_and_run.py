"""continuously run daisy pipeline."""

import os
import sys
import time
import glob
import yaml
import cgatcore.experiment as E
import cgatcore.iotools as IOTools


def main(argv=sys.argv):

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-d", "--data-directory", dest="data_directory", type="string",
        help="directory in which to create links. "
        "[%default]")

    parser.add_option(
        "-f", "--force-first", dest="force_first", action="store_true",
        help="force running of pipeline in first instance [%default]")

    parser.add_option(
        "-k", "--keep-level", dest="keep_level", type="int",
        help="level to keep for the directories that are being watched [%default]")

    parser.add_option(
        "-c", "--command", dest="command", type="string",
        help="command to run when new data appears [%default]")

    parser.set_defaults(
        data_directory="data",
        last_update=3600,
        input_fastq_file=None,
        keep_level=3,
        force_first=False,
        sleep=60,
        command="daisy run -v 5 -p 100 make all",
    )

    (options, args) = E.start(parser, argv)

    if not os.path.exists("benchmark.yml"):
        raise ValueError("config file {} does not exist".format(
            "benchmark.yml"))

    with IOTools.open_file("benchmark.yml") as inf:
        config = yaml.load(inf, Loader=yaml.FullLoader)

    if "watch" not in config:
        raise ValueError("config file needs to contain a 'watch' section")

    if isinstance(config["watch"], list):
        watchlist = config["watch"]
    else:
        watchlist = [config["watch"]]

    E.info("watching with {} glob expressions".format(len(watchlist)))

    while 1:

        current_time = time.time()

        c = E.Counter()

        iteration = 1

        for glob_expr in watchlist:
            filenames = glob.glob(glob_expr)
            E.debug("found {} files for {}".format(len(filenames),
                                                   glob_expr))

            for fn in filenames:
                c.found += 1
                parts = os.path.abspath(fn).split(os.sep)

                dest_fn = os.path.abspath(
                    os.path.join(
                        options.data_directory,
                        os.sep.join(parts[-options.keep_level:])))

                dirname = os.path.dirname(dest_fn)
                if not os.path.exists(dirname):
                    E.info("creating new directory {}".format(dirname))
                    os.makedirs(dirname)

                if not os.path.exists(dest_fn):
                    modification_time = os.path.getmtime(fn)
                    timedelta = current_time - modification_time
                    if timedelta > options.last_update:
                        E.info("new file detected, creating link: {}".format(dest_fn))
                        c.new_file_create += 1
                        os.symlink(os.path.abspath(fn), dest_fn)
                    else:
                        E.info(
                            "new file detected, but too recent ({}s): {}".format(
                                timedelta,
                                dest_fn))
                        c.new_file_wait += 1
                else:
                    c.existing += 1

        E.info("iteration {}: {}".format(iteration, str(c)))

        if iteration == 1 and options.force_first:
            E.run(options.command)
        elif c.new_file_create == 0:
            E.info("found no new files, waiting for {} seconds".format(options.sleep))
            time.sleep(options.sleep)
        else:
            E.run(options.command)

        iteration += 1

if __name__ == "__main__":
    sys.exit(main(sys.argv))
