"""Workflow construction
========================

This module contains methods for constructing workflows.
The workflow manager currently in use is :ref:`ruffus`.

API
---
"""


import os
import ruffus
import itertools
import copy
import re
import glob
import collections
import pandas as pd

import cgatcore.iotools as IOTools
import cgatcore.experiment as E
import cgatcore.pipeline as P
from daisy.storage import upload_result, export_result
from daisy.combinations import build_combinations
from daisy.tasks.Runner import EmptyRunner
from daisy.tasks.ToolRunner import run_tool_identity


def build_output(taskf, result_dir):

    multiple_outputs = isinstance(taskf.output, list)
    if not multiple_outputs:
        flexible_outputs = re.search("\{\S*\}", taskf.output) is not None
    else:
        flexible_outputs = False

    if multiple_outputs and flexible_outputs:
        raise NotImplementedError(
            "mixing multiple and variable outputs not supported")

    if multiple_outputs:
        output = [os.path.join(result_dir, x) for x in taskf.output]
        suffixes = [re.sub("^[^.]*\.", "",
                           os.path.basename(x)) for x in output]
    elif flexible_outputs:
        output = os.path.join(result_dir, taskf.output)
        output = re.sub("\{\S*\}", "*", output)
        suffixes = [re.sub("^[^.]*\.", "", os.path.basename(output))]
    else:
        output = os.path.join(result_dir, taskf.output)
        suffixes = [re.sub("^[^.]*\.", "", os.path.basename(output))]

    suffix = None
    for s in suffixes:
        if suffix is None:
            suffix = s
        elif s != suffix:
            raise ValueError(
                "tool {} produces output files of different type, "
                "got {}, expected {}".format(taskf.name, s, suffix))

    return output, multiple_outputs, flexible_outputs, suffix


def group_files(config_files, group_regex, group_alias="\\1"):
    """group input files by regular expression"""

    rx = re.compile(group_regex)

    for key, files in list(config_files.items()):

        if isinstance(files, list):
            groups = collections.defaultdict(list)
            unmatched = []
            for fn in sorted(files):
                r = rx.search(fn)
                if r is None:
                    unmatched.append(fn)
                    continue
                group_name = r.expand(group_alias)
                groups[group_name].append(fn)

            if len(unmatched) == len(files):
                pass
            elif len(unmatched) == 0:
                config_files[key] = [{x: y} for x, y in list(groups.items())]
            else:
                raise ValueError(
                    "input files not matching regular expression {}: {}"
                    .format(group_regex, str(unmatched)))

    return config_files


def expand_globs(config, is_test=False):
    """detect and expand glob expressions in the input section.

    A glob expression is any filename that contains a '*'. Multiple
    glob expressions can be combined on the same line by a ','.

    A "find" expression is detected starting with 'find'. These
    expressions will be evaluated in a shell and the results insterted
    into the dictionary.

    If a filename starts with "file=", the contents of the file
    following the "=" are read and inserted. Multiple files can be
    separated by a ','.

    If a glob or find expression is evaluated to nothing, an exception
    is raised unless ``is_test`` is set. In that case, two files will be
    returned called "test1" and "test2".
    """

    for d, key, value in IOTools.nested_iter(config):
        if isinstance(value, str) and not (isinstance(key, str) and key.endswith("_regex")):
            if value.startswith("find"):
                try:
                    data = E.run(value, return_stdout=True)
                except Exception as e:
                    data = e.output
                d[key] = [x for x in data.split("\n") if x]
            elif "*" in value:
                if "," in value:
                    v = [glob.glob(x.strip()) for x in value.split(",")]
                    v = [item for sublist in v for item in sublist]
                else:
                    v = glob.glob(value)
                d[key] = v
            elif value.startswith("file="):
                filenames = [x.strip() for x in value.split("=")[1].split(",")]
                paths = []
                for fn in filenames:
                    with IOTools.open_file(fn) as inf:
                        paths.extend([x.strip() for x in inf if x.strip()])
                d[key] = paths
            if len(d[key]) == 0:
                if not is_test:
                    raise ValueError(
                        "expression '{}' expanded to nothing".format(value))
                else:
                    # insert some random files for testing purposes:
                    if "*" in value:
                        # replace glob expressions
                        value = re.sub(",.*", "", value)
                        d[key] = [re.sub("[*]", "test1", value),
                                  re.sub("[*]", "test2", value)]
                    else:
                        if "bam" in value:
                            d[key] = ["test1.bam", "test2.bam"]
                        elif "vcf" in value:
                            d[key] = ["test1.vcf.gz", "test2.vcf.gz"]
                        else:
                            d[key] = ["test1.txt", "test2.txt"]
    return config


def check_unique(tool_functions,
                 input_combos=None,
                 input_regex=None,
                 input_alias=None,
                 is_test=False):
    # compute a list of task names
    names = []
    if input_combos:
        for toolf, input_files in itertools.product(tool_functions,
                                                    input_combos):
            taskf = copy.copy(toolf)
            taskf.register_input(input_files,
                                 regex=input_regex,
                                 alias=input_alias,
                                 is_test=is_test)
            names.append(taskf.__name__)
    else:
        for toolf in tool_functions:
            taskf = copy.copy(toolf)
            taskf.register_input(regex=input_regex,
                                 alias=input_alias,
                                 is_test=is_test)
            names.append(taskf.__name__)

    counts = collections.Counter(names)
    for name, count in list(counts.items()):
        if count > 1:
            make_unique = True
            P.get_logger().debug(
                "adding hash identifier because of duplicate name: {}={}".format(name, count))
            break
    else:
        make_unique = False

    return make_unique


def expand_generators(config):
    """expand generator expressions in option lists.

    A generator expression are valid python syntax and
    has the following syntax::

      options: generate=["--chrom={}".format(x) for x in [1,2,3,4,5]]

    """

    to_delete = []
    for d, key, value in IOTools.nested_iter(config):
        if isinstance(value, str):
            if value.startswith("generate="):
                expression = re.sub("^generate=\s*", "", value)
                if expression.startswith("'") and expression.startswith("'"):
                    expression = expression[1:-1]
                try:
                    argument_list = eval(expression)
                except SyntaxError as ex:
                    raise ValueError(
                        "error occured while evaluating generator "
                        "expression {}: {}".format(expression, ex))
                if isinstance(d, list):
                    d.extend(argument_list)
                    to_delete.append((d, key))
                else:
                    d[key] = argument_list

    for d, key in to_delete[::-1]:
        del d[key]

    return config


def build_tool_functions(map_tool_to_runner, config):

    try:
        tools = config["setup"]["tools"]
    except (KeyError, TypeError) as ex:
        raise KeyError("configuration file requires a 'setup:tools' section: {} '{}'".format(
            ex, config))

    # create all vs all combinations for tools, tool options
    # and input files
    tool_functions = []

    for tool in tools:
        try:
            toolc = map_tool_to_runner[tool]
        except KeyError:
            raise KeyError("unknown tool '{}', available tools are {}".format(
                tool, ", ".join(sorted(map_tool_to_runner.keys()))))
        if toolc.name in config:
            conf = config[toolc.name]
        else:
            conf = {}
        conf = expand_generators(conf)
        configurations = build_combinations(conf)
        for configuration in configurations:
            tool_functions.append(toolc(**configuration))
    return tool_functions


def add_tools_to_pipeline(pipeline,
                          map_tool_to_runner,
                          config=None,
                          input_files=None,
                          **kwargs):
    """add tools to a workflow pipeline.

    This function adds for each input and tool combination
    a task to the workflow.

    The configuration dictionary should contain the following
    sections:

    input:
       Configuration of input files. Key/value pairs and possibly
       hierarchical.

       The following keys are optional:
          regex
          alias
          group_regex
          group_alias

    tool:
       A list of tools to apply.

    A typical configuration dictionary might look like this::

        {"input": {"bam": "*.bam"}, "tool": ["bwa_mem", "isaac"]}

    Arguments
    ---------
    pipeline : object
        The ruffus pipeline that tasks will be added to.
    map_tool_to_runner: dict
        Dictionary mapping tools to functions in the
        :ref:`tasks`.
    config: dict
        Configuration dictionary.
    input_files: list
        List of (optional) input files.
    """
    tool_functions = build_tool_functions(map_tool_to_runner, config)

    if "input" not in config:
        raise KeyError("configuration file requires an 'input' section")

    if config["input"] is None:
        raise ValueError("input section is empty")

    input_regex = config["input"].pop("regex", None)
    input_alias = config["input"].pop("alias", None)
    replicate_alias = config["input"].pop("replicate_alias", None)
    input_group_regex = config["input"].pop("group_regex", None)
    input_group_alias = config["input"].pop("group_alias", "\\1")

    ignore = config["setup"].get("ignore", [])
    ignore.extend(config["input"].get("ignore", []))

    do_replication = config["setup"].pop("replication", None)
    if do_replication:
        replications = int(do_replication)
        P.get_logger().info("running experiment with {} replications".format(replications))
    else:
        replications = 1

    # update selected fields for testing purposes
    is_test = "is_test" in config
    if "test" in config["input"]:
        config["input"].update(config["input"]["test"])
        del config["input"]["test"]

    # build input/tool combinations, optionally grouping them
    config_files = expand_globs(config["input"], is_test=is_test)
    if input_group_regex:
        config_files = group_files(config_files,
                                   input_group_regex,
                                   input_group_alias)

    input_combos = build_combinations(config_files)
    tool_runners = []

    make_unique = check_unique(tool_functions,
                               input_combos=input_combos,
                               input_regex=input_regex,
                               input_alias=input_alias,
                               is_test=is_test)

    suffix = None

    for toolf, input_files in itertools.product(tool_functions, input_combos):

        for replication_idx in range(replications):
            # create a copy of the task function and give it its unique name
            # by mangling it with the input_files
            taskf = copy.copy(toolf)

            if do_replication:
                taskf.set_replication_id(replication_idx + 1)

            taskf.register_input(input_files,
                                 regex=input_regex,
                                 alias=input_alias,
                                 make_unique=make_unique,
                                 is_test=is_test,
                                 replicate_alias=replicate_alias)

            if "name" in input_files:
                # create copy of input_files without name, do
                # not modify original as different tools require
                # the 'name'
                input_files = dict([(x, y) for x, y in list(input_files.items())
                                    if x != "name"])

            result_dir = os.path.join(taskf.__name__ + ".dir")

            found = False

            for i in IOTools.val2list(ignore):
                if i in result_dir:
                    P.get_logger().warn(
                        "the following task will be ignored: "
                        "{} matching {}".format(
                            result_dir, i))
                    found = True
            if found:
                continue

            output, multiple_outputs, flexible_outputs, _suffix = \
                build_output(taskf, result_dir)
            if suffix is None:
                suffix = _suffix
            elif suffix != _suffix:
                raise ValueError(
                    "tools produce output files of different type, "
                    "got {}, expected {}".format(_suffix, suffix))

            tool_task = pipeline.merge(
                task_func=taskf,
                input=list(input_files.values()),
                output=output,
                **kwargs).mkdir(result_dir)

            # if there are multilpe output files, split the task so that
            # each output file will be processed separately further down the
            # pipeline.
            if multiple_outputs:
                f = EmptyRunner()
                f.__name__ = taskf.__name__ + "_split"
                tool_task = pipeline.split(
                    task_func=f,
                    input=tool_task,
                    output=output)

            tool_runners.append(tool_task)

    # convenience target
    f = EmptyRunner()
    f.__name__ = "tools"
    pipeline.merge(task_func=f,
                   input=tool_runners,
                   output=None)

    return suffix, tool_runners


def add_external_data_to_pipeline(pipeline,
                                  config=None,
                                  **kwargs):

    is_test = "is_test" in config

    external_config = config["setup"]["external"]

    input_regex = external_config["input"].pop("regex", None)
    input_alias = external_config["input"].pop("alias", None)
    input_group_regex = external_config["input"].pop("group_regex", None)
    input_group_alias = external_config["input"].pop("group_alias", "\\1")

    config_files = expand_globs(external_config["input"], is_test=is_test)

    if input_group_regex:
        config_files = group_files(config_files,
                                   input_group_regex,
                                   input_group_alias)

    input_combos = build_combinations(config_files)
    tool_runners = []

    suffix = None

    toolf = run_tool_identity(**external_config)

    for input_files in input_combos:

        taskf = copy.copy(toolf)

        taskf.register_input(input_files,
                             regex=input_regex,
                             alias=input_alias,
                             is_test=is_test)

        result_dir = os.path.join(taskf.__name__ + ".dir")

        output, multiple_outputs, flexible_outputs, _suffix = \
            build_output(taskf, result_dir)
        if suffix is None:
            suffix = _suffix
        elif suffix != _suffix:
            raise ValueError(
                "tools produce output files of different type, "
                "got {}, expected {}".format(_suffix, suffix))

        tool_task = pipeline.merge(
            task_func=taskf,
            input=list(input_files.values()),
            output=output,
            **kwargs).mkdir(result_dir)

        tool_runners.append(tool_task)

    return tool_runners


def add_collations_to_pipeline(pipeline,
                               map_tool_to_runner,
                               collations,
                               tasks=None,
                               config=None,
                               **kwargs):

    runners = []

    ignore = config["setup"].get("ignore", [])
    ignore.extend(config["input"].get("ignore", []))

    for coll in collations:

        if coll not in config:
            raise KeyError(
                "configuration file requires a section for '{}'".format(coll))

        coll_info = config[coll]

        for keyword in ("runner", "regex_in", "pattern_out"):
            if keyword not in coll_info:
                raise ValueError("section {} is missing required keyword '{}'".format(
                    coll, keyword))

        runner_options = config.get(coll_info["runner"], {})
        runner_name = runner_options.get("name", coll_info["runner"]).strip()

        colcc = map_tool_to_runner[runner_name]
        taskf = colcc(**runner_options)

        # automatically set alias through regex (required field)
        taskf._input_regex = coll_info.get("regex", None)
        taskf._input_alias = coll_info.get("alias", None)
        taskf._replicate_regex = coll_info.get("regex_replicate", None)
        taskf.__name__ = coll

        if tasks is not None:
            input_tasks = tasks
        elif "glob" in coll_info:
            input_tasks = coll_info["glob"]
        else:
            raise ValueError("need either tasks or glob expression "
                             "for collation")

        filter_regex = ruffus.regex(coll_info["regex_in"])

        filter_regex = ruffus.regex(coll_info["regex_in"])
        result_dir = os.path.join(coll + ".dir")

        output_pattern = coll_info["pattern_out"]
        output_prefix = r"{}/{}".format(result_dir, output_pattern)
        output_dir = os.path.dirname(output_prefix)

        if hasattr(taskf, "output"):
            output, multiple_outputs, flexible_outputs, _suffix = \
                build_output(taskf, output_dir)
        else:
            multiple_outputs = False
            output = output_prefix

        found = False
        for i in IOTools.val2list(ignore):
            if i in result_dir:
                P.get_logger().warn(
                    "the following task will be ignored: "
                    "{} matching {}".format(
                        result_dir, i))
                found = True
        if found:
            continue

        metric_task = pipeline.collate(
            task_func=taskf,
            input=input_tasks,
            filter=filter_regex,
            output=output,
            **kwargs).mkdir(
                input_tasks,
                filter_regex,
                output_dir)

        if multiple_outputs:
            f = EmptyRunner()
            f.__name__ = taskf.__name__ + "_passthrough"
            output = [re.sub(r"\\\d+", "*", x) for x in output]
            metric_task = pipeline.split(
                task_func=f,
                input=metric_task,
                output=output)

        runners.append(metric_task)

    return runners


def add_splits_to_pipeline(pipeline,
                           map_tool_to_runner,
                           tool_runners,
                           splits,
                           tasks=None,
                           config=None,
                           **kwargs):

    split_functions = []
    for split in splits:
        if split not in config:
            raise KeyError(
                "configuration file requires a section for '{}'".format(split))

        split_info = config[split]
        splitcc = map_tool_to_runner[split_info["runner"].strip()]

        if "runner" in split_info:
            runner_conf = config.get(split_info["runner"], {})
        else:
            runner_conf = {}

        configurations = build_combinations(runner_conf)
        for configuration in configurations:
            split_functions.append((splitcc(**configuration), split_info))

    suffix = None
    runners = []

    for taskf, split_info in split_functions:

        # automatically set alias through regex (required field)
        taskf._input_regex = split_info.get("regex", None)
        taskf._input_alias = split_info.get("alias", None)

        taskf.register_input()
        unique_name = taskf.__name__
        result_dir = os.path.join(unique_name + ".dir")
        output, multiple_outputs, flexible_outputs, _suffix = \
            build_output(taskf, result_dir)
        if suffix is None:
            suffix = _suffix
        elif suffix != _suffix:
            raise ValueError(
                "tools produce output files of different type, "
                "got {}, expected {}".format(_suffix, suffix))

        filter_regex = ruffus.regex("(.*)/(.*).{}".format(suffix))

        if multiple_outputs:
            output = [r"\1/{}".format(x) for x in output]
            outdirs = set([os.path.dirname(x) for x in output])
        elif flexible_outputs:
            output = r"\1/{}".format(output)
            outdirs = os.path.dirname(os.path.dirname(output))

        if multiple_outputs and len(outdirs) == 1:
            raise ValueError("multiple outputs going into the same directory")

        split_task = pipeline.subdivide(
            task_func=taskf,
            input=tool_runners,
            filter=filter_regex,
            output=output,
            **kwargs).mkdir(tool_runners, filter_regex, outdirs)

        runners.append(split_task)

    return runners


def add_metrics_to_pipeline(pipeline,
                            metrics,
                            map_metric_to_runner,
                            tool_runners,
                            suffix="tsv",
                            prefix=None,
                            config=None,
                            **kwargs):

    single_input_metric_functions = []

    for metric in metrics:
        metricc = map_metric_to_runner[metric.strip()]
        if metricc.name in config:
            conf = config[metricc.name]
        else:
            conf = {}

        conf = expand_generators(conf)
        configurations = build_combinations(conf)
        for configuration in configurations:
            single_input_metric_functions.append(metricc(**configuration))

    make_unique = check_unique(single_input_metric_functions,
                               input_combos=None,
                               input_regex=None,
                               input_alias=None)

    metric_runners = []
    for taskf in single_input_metric_functions:

        ignore = config.get(taskf.name, {}).get("ignore", [])
        taskf.register_input(make_unique=make_unique)
        unique_name = taskf.__name__

        # make task name unique by adding 'prefix' as this method might
        # be called multiple times for straight, collated and split tasks
        if prefix:
            taskf.__name__ = prefix + taskf.__name__

        filter_regex = ruffus.regex("(.*)/(.*).{}".format(suffix))
        result_dir = os.path.join(unique_name + ".dir")
        output = r"\1/{}/{}.tsv".format(result_dir, taskf.name)

        found = False
        # Note that ignore will only work on the static parts of a task
        # as result_dir contains a pattern that will be filled in at runtime,
        # e.g. \1/echidna_test.dir/echidna_test.tsv.
        for i in ignore:
            if i in result_dir:
                P.get_logger().warn(
                    "the following task will be ignored: "
                    "{} matching {}".format(
                        result_dir, i))
                found = True

        if found:
            continue

        metric_task = pipeline.transform(
            task_func=taskf,
            input=tool_runners,
            filter=filter_regex,
            output=output,
            **kwargs)

        metric_runners.append(metric_task)

    f = EmptyRunner()
    if prefix:
        f.__name__ = prefix + "metrics"
    else:
        f.__name__ = "metrics"
    pipeline.merge(task_func=f,
                   input=metric_runners,
                   output=None)

    return metric_runners


def add_upload_to_pipeline(pipeline, metric_runners, CONFIG):

    if "title" not in CONFIG:
        raise ValueError("refusing to upload data without a 'title'")

    if "description" not in CONFIG:
        raise ValueError("refusing to upload data without a 'description'")

    if "tags" not in CONFIG:
        raise ValueError("refusing to upload data without 'tags'")

    upload_task = pipeline.merge(
        task_func=upload_result,
        input=metric_runners,
        output="results.commit",
        extras=[CONFIG]).jobs_limit(1, "DB")

    pipeline.merge(
        task_func=EmptyRunner(name="upload"),
        input=upload_task,
        output="upload")


def add_export_to_pipeline(pipeline, tool_runners, suffix,
                           config, **kwargs):

    conf = config.get("export", {})
    prefix = conf.get("prefix", "").strip()

    result_dir = "export.dir"
    filter_regex = ruffus.regex("(.*).dir/(.*).{}".format(suffix))
    output = r"{}/{}\1.{}".format(result_dir, prefix, suffix)

    export_result.__name__ = "export"

    pipeline.transform(
        task_func=export_result,
        input=tool_runners,
        filter=filter_regex,
        output=output,
        **kwargs).mkdir(result_dir)


def add_all_task_to_pipeline(pipeline, metric_runners):
    pipeline.merge(
        task_func=EmptyRunner(name="all"),
        input=metric_runners,
        output="all")


if __name__ == "__main__":
    import doctest
    doctest.testmod()
