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

from daisy.storage import upload_result


def main(argv):

    options, args = P.parse_commandline(argv)

    if options.config_file:
        PARAMS = P.get_parameters(options.config_file)
    else:
        sys.exit(P.main(options, args))

    if os.path.exists("results.commit"):
        raise ValueError(
            "a results.commit file already exists. Please remove "
            "before uploading.")

    E.info("collecting files to upload")
    infiles = []
    for root, dirs, files in os.walk("."):
        # ignore first level (tools)
        if root.count(os.sep) == 1:
            continue
        if "benchmark.info" in files:
            infiles.append(os.path.join(root, "benchmark.info"))

    E.info("found a potential {} benchmark.info files to upload".format(len(infiles)))

    # find all files of interest
    upload_result(infiles, "results.commit", PARAMS)

    E.stop()


if __name__ == "__main__":
    sys.exit(main(sys.argv[:]))
