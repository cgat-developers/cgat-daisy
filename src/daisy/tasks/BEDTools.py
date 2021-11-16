import re

from .CollateRunner import CollateRunner
from .ToolRunner import ToolRunner
from .Runner import resolve_argument
import cgatcore.pipeline as P
import cgatcore.experiment as E


class run_collate_bed_cat(CollateRunner):
    """concatenate multiple bed files.

    :param infiles: Input files in :term:`bed` format
    :param outfile: Output file in :term:`bed` format

    The output file will be coordinate sorted and indexed.
    """

    name = "bed_cat"

    version = "unknown"

    def run(self, infiles, outfile, params):

        files = " ".join(infiles)

        return P.run("zcat {files} "
                     "| sort -k 1,1 -k2,2n "
                     "| bgzip > {outfile}; "
                     "tabix -p bed {outfile} ".format(**locals()))


class CollateRunnerBedtools(CollateRunner):
    path = "bedtools"

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stdout=True).strip()
        return re.search("bedtools (\S+)", help_string).groups()[0]


class run_collate_bedtools_merge(CollateRunnerBedtools):
    """run bedtools merge command on multiple input files.

    :param infiles: Input files in :term:`bed` format
    :param outfile: Output file in :term:`bed` format

    The input files will be sorted by position before calling bedtools
    merge.

    """

    name = "bedtools_merge"

    def run(self, infiles, outfile, params):

        files = " ".join(infiles)

        return P.run("zcat {files} "
                     "| sort -k 1,1 -k2,2n "
                     "| {params.path} merge {params.options} "
                     "| bgzip > {outfile}; "
                     "tabix -p bed {outfile} ".format(**locals()))


class ToolRunnerBedtools(ToolRunner):
    path = "bedtools"

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stdout=True).strip()
        return re.search("bedtools (\S+)", help_string).groups()[0]


class run_tool_gtf2cds_bed(ToolRunnerBedtools):
    """extract CDS regions from a GTF file.

    The BED file output contains only the region coordinates
    and gene ids, etc.
    """

    name = "gtf2cds_bed"

    expected = ["gtf"]
    output = "result.bed.gz"

    merge_regions = False
    merge_options = ""

    def run(self, outfile, params):

        if params.merge_regions:
            merge_command = "| bedtools merge {params.merge_options} -i -".format(
                **locals())
        else:
            merge_command = ""

        gtf = resolve_argument(params.gtf, " ")
        statement = (
            "zless {gtf} "
            "| awk '$3 == \"CDS\" "
            "{{printf(\"%%s\\t%%i\\t%%i\\n\", "
            "$1, $4-1, $5) }}' "
            "| sort -k1,1 -k2,2n "
            "{merge_command} "
            "| bgzip "
            "> {outfile}; "
            "tabix -p bed {outfile}".format(**locals()))

        return P.run(statement)
