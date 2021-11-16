"""Run a task in the task library
=================================

This command provides command line access to task from the task
library.

To run a task from the task library, specify the task, input
and output files. For example::

   benchmark run-task --task=metric_bcftools_stats
     --input-file=x.vcf.gz
     --output-file=test.tsv

This command will run the metric bcftools_stats on the input file
x.vcf.gz and output test.tsv as well as some additional files.

To supply command line options to the tool, use `--` followed by the options::

   benchmark run-task --task=metric_bcftools_stats
     --input-file=x.vcf.gz
     --output-file=test.tsv
     -- --options='--fasta-ref=/data/library/reference/hs37d5/hs37d5.fa'

Running tools requires setting input slots for the input files::

   benchmark run-task --task=tool_weCall \
      --input-file=/data/library/SmallVariantCaller/bam/37/WES/final_bam/NA12878.bam
      --input-slot=bam
      --input-file=/data/library/reference/hs37d5/hs37d5.fa
      --input-slot=reference_fasta
      --output-file=test.vcf.gz
      -- --path=/local/scratch/test_mapping/weCall-0.8.1-Linux-x86_64-Jul-29-10/bin/weCall
         --options=--regions=20

The names of the slots depend on the type of tool invoked.

"""

import sys
import tempfile
import os
import shutil
import signal
import cgatcore.experiment as E
import cgatcore.pipeline as P
from daisy.toolkit import redirect2mounts

from daisy.tasks import map_tool_to_runner, \
    map_metric_to_runner, \
    map_collate_to_runner, \
    map_split_to_runner, \
    redirect_defaults2mountpoint


def cleanup(signal=None, frame=None):

    mountpoint = P.PARAMS.get("mount_point", None)
    if mountpoint:
        E.run("fusermount -u {}".format(mountpoint))
        E.debug("finally: unmounted arvados at {}".format(mountpoint))
        shutil.rmtree(mountpoint)
        P.PARAMS["mount_point"] = None


def parse_args(args):
    name = None
    for arg in args:
        if arg.startswith("--"):
            if "=" in arg:
                name, value = arg[2:].split("=", 1)
                yield name, value
                name = None
            else:
                name = arg[2:]
        else:
            yield (name, arg)
            name = None
    if name is not None:
        yield (name, None)


def main(argv=sys.argv):

    TASKS = {}
    for label, collection in [("tool", map_tool_to_runner),
                              ("metric", map_metric_to_runner),
                              ("collate", map_collate_to_runner),
                              ("split", map_split_to_runner)]:
        for key, f in list(collection.items()):
            k = "{}_{}".format(label, key)
            if k in TASKS:
                raise ValueError("duplicate keys in TASK: {} {} {}"
                                 .format(k, TASKS[k], f))
            TASKS[k] = f

    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option(
        "-i", "--input-file", dest="input_files", type="string", action="append",
        help="input file. Can be used more than once [%default]")

    parser.add_option(
        "-s", "--input-slot", dest="input_slots", type="string", action="append",
        help="input slot. Must be used as often as input_files for tools [%default]")

    parser.add_option(
        "-o", "--output-file", dest="output_files", type="string", action="append",
        help="output file. Can be used more than once [%default]")

    parser.add_option(
        "-n", "--dry-run", dest="dry_run", action="store_true",
        help="show statement to be executed, do not execute [%default]")

    parser.add_option(
        "--engine", dest="engine", type="choice",
        choices=("local", "arvados"),
        help="engine to use [%default]")

    parser.add_option(
        "-t", "--task", dest="task", type="choice",
        choices=sorted(TASKS.keys()),
        help="task to run [%default]")

    parser.add_option(
        "-l", "--list-tasks", dest="list_tasks", action="store_true",
        help="list all available tasks and exit [%default]")

    parser.add_option(
        "--always-mount", dest="always_mount",
        action="store_true",
        help="force mounting of arvados keep [%default]")

    parser.set_defaults(
        input_files=[],
        input_slots=[],
        output_files=[],
        engine="local",
        dry_run=False,
        task=None,
        always_mount=False,
    )

    (options, args) = E.start(parser, argv, add_cluster_options=True)

    if options.list_tasks:
        options.stdout.write("available_tasks\n{}\n".format(
            "\n".join(sorted(TASKS.keys()))))
        E.stop()
        return

    if len(options.input_files) == 0:
        raise ValueError("no input files specified, use --input-file")

    if len(options.output_files) == 0:
        raise ValueError("no output files specified, use --output-file")

    if options.task is None:
        raise ValueError("please specify a task to run (--task)")

    P.get_parameters()

    if options.engine == "arvados":

        raise ValueError("arvados support disabled")
        # crunch_json = Arvados.build_crunch_script(argv)
        crunch_json = None
        retval = E.run('arv-crunch-job --job="$(cat {})"'.format(crunch_json))

        if retval != 0:
            raise ValueError("error while executing")

        os.unlink(crunch_json)
        E.stop()
        return retval

    # Start SGE session
    if not options.without_cluster:
        P.start_session()

    params = dict(parse_args(args))

    signal.signal(signal.SIGINT, cleanup)

    # redirect all mount points in parameters and input files.
    mountpoint = redirect2mounts(
        [params, options.input_files],
        always_mount=options.always_mount)
    mountpoint = redirect_defaults2mountpoint(mountpoint)
    # register mountpoint with pipeline
    P.PARAMS["mount_point"] = mountpoint
    P.PARAMS["dryrun"] = options.dry_run

    try:
        # instantiate task runner
        runner = TASKS[options.task](**params)

        if len(options.output_files) == 0:
            tmpfile = tempfile.NamedTemporaryFile(delete=False)
            os.unlink(tmpfile.name)
            options.output_files.append(tmpfile.name)

        if options.task.startswith("tool"):
            if len(options.input_slots) != len(options.input_files):
                raise ValueError(
                    "for tools, provide the same number as input slots as there"
                    "are input files (--input-slots)")

            input_files = dict(zip(options.input_slots,
                                   options.input_files))

            runner.register_input(input_files)
            # check if expected is in params
            runner(list(input_files.values()),
                   options.output_files[0])
        elif options.task.startswith("metric"):
            runner(options.input_files[0],
                   options.output_files[0])
        elif options.task.startswith("collate"):
            runner(
                options.input_files,
                options.output_files[0])
        elif options.task.startswith("split"):
            runner(
                options.input_files[0],
                options.output_files)

        # stop SGE session
        P.close_session()

    finally:
        cleanup()

    E.stop()

if __name__ == "__main__":
    sys.exit(main())
