"""A basic benchmarking workflow
=================================

Implement the basic benchmarking workflow:

1. Apply a set of tools to a collection of input files.

2. Apply a set of metrics on each of the outputs from the
   tools

Restrictions
------------

All tools need to work on the same type of input data (e.g. :term:`BAM`) and
produce the same type of output (e.g. :term:`VCF`).

"""
import sys
import ruffus
import collections
import cgatcore.pipeline as P
import cgatcore.experiment as E
from daisy.toolkit import arvados_enabled

from daisy.workflow import add_tools_to_pipeline, \
    add_metrics_to_pipeline, \
    add_external_data_to_pipeline, \
    add_collations_to_pipeline, \
    add_splits_to_pipeline, \
    add_upload_to_pipeline, \
    add_export_to_pipeline, \
    add_all_task_to_pipeline

# import tasks to apply in this pipeline
from daisy.tasks import map_tool_to_runner, \
    map_metric_to_runner, \
    map_collate_to_runner, \
    map_split_to_runner, \
    redirect_defaults2mountpoint


def main(argv):

    options = P.initialize(argv, config_file="benchmark.yml")

    # compatibility with cgatcore < 0.6.3
    if isinstance(options, tuple):
        options = options[0]

    # not sure what this does
    # if not options.config_file:
    #     P.get_parameters(options.config_file)
    # else:
    #     sys.exit(P.main(options, args))

    params = P.get_params()
    
    with arvados_enabled(always_mount=options.always_mount):
        mountpoint = params.get("mount_point", None)
        if mountpoint:
            redirect_defaults2mountpoint(mountpoint)

        # A selection of command line arguments are added to PARAMS
        # as 'extras' not implemented in ruffus 2.6.3
        kwargs = collections.defaultdict(dict)
        if options.only_info:
            kwargs["extras"].update({'only_info': True})
            P.PARAMS["only_info"] = True
        if options.is_test:
            kwargs["extras"].update({'is_test': True})
            P.PARAMS["is_test"] = True

        E.debug("construction of workflow started")
        pipeline = ruffus.Pipeline('benchmark')
        # Tool execution
        suffix, tool_runners = add_tools_to_pipeline(
            pipeline,
            map_tool_to_runner,
            config=P.PARAMS,
            **kwargs)

        E.debug("added {} tools to workflow".format(len(tool_runners)))
        # Optionally, add externally computed files as
        # pseudo-tools:
        if "external" in P.PARAMS["setup"]:
            external_runners = add_external_data_to_pipeline(
                pipeline,
                config=P.PARAMS,
                **kwargs)
            tool_runners.extend(external_runners)

        # Optionally, combine tool runs into aggregate
        # outputs. The type of the output is preserved
        # (VCF -> VCF, etc.)
        # For example, call individual members in a trio
        # and then build a combined VCF to analyse mendelian
        # inconsistencies.
        if "collate" in P.PARAMS["setup"]:
            collate_runners = add_collations_to_pipeline(
                pipeline,
                map_collate_to_runner,
                P.PARAMS["setup"]["collate"],
                tasks=tool_runners,
                config=P.PARAMS)
            if P.PARAMS["setup"].get("only_collate", False):
                tool_runners = []
            if P.PARAMS["setup"].get("no_collate_metrics", False):
                collate_runners = []
            E.debug("added {} collators to workflow".format(len(collate_runners)))
        else:
            collate_runners = []

        # Optionally, split up the output before applying
        # additional analyses. The type of the output is preserved
        # (VCF -> VCF, etc).
        # For example, identify false positives, false negatives
        # and true positives and collect metrics individually.
        if "split" in P.PARAMS["setup"]:
            split_runners = add_splits_to_pipeline(
                pipeline,
                map_split_to_runner,
                tool_runners,
                P.PARAMS["setup"]["split"],
                tasks=tool_runners,
                config=P.PARAMS)
            if P.PARAMS["setup"].get("only_split", False):
                tool_runners = []
            E.debug("added {} splitters to workflow".format(len(split_runners)))
        else:
            split_runners = []

        metric_runners = []
        for prefix, r in zip(["tool", "collate", "split"],
                             [tool_runners, collate_runners, split_runners]):
            if not r:
                continue

            metrics = None

            if prefix == "collate" and "collate_metrics" in P.PARAMS["setup"]:
                metrics = P.PARAMS["setup"]["collate_metrics"]
            elif prefix == "split" and "split_metrics" in P.PARAMS["setup"]:
                metrics = P.PARAMS["setup"]["split_metrics"]
            elif "metrics" in P.PARAMS["setup"]:
                metrics = P.PARAMS["setup"]["metrics"]
            else:
                raise KeyError("configuration file requires a 'setup:metrics' section")

            # Metric execution
            mm = add_metrics_to_pipeline(
                pipeline,
                metrics,
                map_metric_to_runner,
                r,
                suffix=suffix,
                prefix=prefix + "_",
                config=P.PARAMS,
                **kwargs)

            if len(mm) == 0:
                raise ValueError("workflow construction error: "
                                 "no metric tasks result for metrics {}".format(
                                     metrics))

            metric_runners.extend(mm)
            E.debug("added {} {}_metrics to workflow".format(len(mm), prefix))

        # add plot task
        if "aggregate" in P.PARAMS["setup"]:
            aggregate_metrics = add_collations_to_pipeline(
                pipeline,
                map_collate_to_runner,
                P.PARAMS["setup"]["aggregate"],
                metric_runners,
                config=P.PARAMS)

            E.debug("added metric aggregation to workflow")
        else:
            aggregate_metrics = []

        add_upload_to_pipeline(pipeline,
                               metric_runners + aggregate_metrics,
                               P.PARAMS)
        E.debug("added upload to workflow".format(prefix))

        # add export task
        export = P.PARAMS["setup"].get(
            "export",
            ["tools", "collate", "split"])
        map_export2runner = {
            "collate": collate_runners,
            "tools": tool_runners,
            "split": split_runners}

        export_runners = []
        for e in export:
            try:
                export_runners.extend(map_export2runner[e])
            except KeyError:
                raise KeyError("unknown export section: {}".format(e))

        add_export_to_pipeline(
            pipeline,
            export_runners,
            suffix=suffix,
            config=P.PARAMS)

        E.debug("added export to workflow")

        add_all_task_to_pipeline(pipeline, metric_runners + aggregate_metrics)

        # Collate output files to facilitate analysis
        if "collation" in P.PARAMS:
            collators = add_collations_to_pipeline(
                pipeline,
                map_collate_to_runner,
                P.PARAMS["collation"],
                config=P.PARAMS)

        E.debug("construction of workflow completed")

        E.debug("starting workflow")
        P.run_workflow(options, pipeline=pipeline)

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv[:]))
