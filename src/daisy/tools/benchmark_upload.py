"""Upload results from a benchmark run
======================================

This script uploads the results of a benchmark run without running the
full workflow. This might be useful if the workflow is complete and
the data need to be uploaded, but the original data have since gone or
been updated.

note::

   Preferably use ``daisy run upload`` to upload data. The implementation
   here is for cases where patching is required.

"""
import sys
import os
import cgatcore.experiment as E
import cgatcore.pipeline as P
import cgatcore.iotools as IOTools
from daisy.storage import upload_result


def main(argv):

    def _add_input(parser):
        parser.add_option("--data-dir", default=".")
        parser.add_option("--force", default=False, action="store_true")
        parser.add_option("--min-depth", default=0, type="int")
        parser.add_option("--follow-links", default=False, action="store_true")
        parser.add_option("--limit-metrics", default=0, type="int")
        parser.add_option("--output-filename-metrics")
        parser.add_option("--input-filename-metrics")

    P.initialize(argv, callback=_add_input)
    options = E.get_args()

    if options.config_file:
        PARAMS = P.get_parameters(options.config_file)
    else:
        sys.exit(P.main(options))

    if os.path.exists("results.commit"):
        if not options.force:
            raise ValueError(
                "a results.commit file already exists. Please remove "
                "before uploading.")

    data_dir = os.path.abspath(options.data_dir)
    if options.input_filename_metrics:
        with IOTools.open_file(options.input_filename_metrics) as inf:
            infiles = [x.strip() for x in inf.readlines() if x.strip()]
        if options.limit_metrics:
            infiles = infiles[:options.limit_metrics]
    else:
        E.info(f"collecting files to upload starting in {data_dir}")
        infiles = []
        for root, dirs, files in os.walk(data_dir, followlinks=options.follow_links):
            E.debug(f"working on {root}: dirs={len(dirs)}, files={len(files)}")
            # ignore first level (tools) (needs better check)
            depth = root[len(data_dir):].count(os.sep)
            if "benchmark.info" in files:
                if depth <= options.min_depth:
                    E.info(f"skipping - depth not high enough: {depth}")
                else:
                    infiles.append(os.path.join(root, "benchmark.info"))

            if options.limit_metrics and len(infiles) > options.limit_metrics:
                E.info(f"stopping collection as {len(infiles)} reached")
                break

    E.info("found a potential {} benchmark.info files to upload".format(len(infiles)))
    if options.output_filename_metrics:
        with IOTools.open_file(options.output_filename_metrics, "w") as outf:
            outf.write("\n".join(infiles) + "\n")

    # find all files of interest
    oldwd = os.getcwd()
    os.chdir(data_dir)
    upload_result(infiles, "results.commit", PARAMS)
    os.chdir(oldwd)

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv[:]))
