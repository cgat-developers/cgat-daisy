import re
import pysam

from .MetricRunner import MetricRunner
from .SplitRunner import SplitRunner

import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools


def standardise_bed_files(
        outfile1, outfile2,
        bedfile1, bedfile2):
    """create two bed files containing only chromosomes that are
    shared between both sets.

    Ensure files are sorted in the same order and remove
    any "chr" prefixes.

    """
    tbx = pysam.TabixFile(bedfile1)
    contigs_a = set([re.sub("^chr", "", x) for x in tbx.contigs])
    tbx.close()

    tbx = pysam.TabixFile(bedfile2)
    contigs_b = set([re.sub("^chr", "", x) for x in tbx.contigs])
    tbx.close()

    common = contigs_a.intersection(contigs_b)

    contig_filter = " || ".join(
        ['$1 == "{}"'.format(x) for x in common])

    stmnt = ("zcat {bedfile1} "
             "| perl -p -e 's/^chr//' "
             "| awk '{contig_filter}' "
             "| sort -k1,1 -k2,2n "
             "| bgzip > {outfile1}; "
             "tabix -p bed {outfile1}; "
             "zcat {bedfile2} "
             "| perl -p -e 's/^chr//' "
             "| awk '{contig_filter}' "
             "| sort -k1,1 -k2,2n "
             "| bgzip > {outfile2}; "
             "tabix -p bed {outfile2}").format(**locals())

    return stmnt


class MetricRunnerBedtools(MetricRunner):
    path = "bedtools"

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stdout=True).strip()
        return re.search("bedtools (\S+)", help_string).groups()[0]


class run_metric_bedtools_jaccard(MetricRunnerBedtools):
    """run bedtools jaccard on a :term:`bed` file comparing
    against a rerference :term:`bed` file.

    :param infile: Input file in :term:`bed` format
    :param outfile: Output file in :term:`tsv` format

    *Columns*

    intersection
        number of nucleotides in the intersection of both interval sets
    union-intersection
        number of nucleotides in the union of both interval sets minus the
        intersection
    jaccard
        jaccard coefficiont, the ratio of intersection and union-intersection.
    n_intersections
        number of intersections
    """

    name = "bedtools_jaccard"
    reference_bed = None

    def run(self, infile, outfile, params):

        if params.reference_bed is None:
            raise ValueError("{} requires reference_bed to be set".format(
                self.name))

        # jaccard requires a consistent sort order, so sort both
        # bed files:
        tmpf = P.get_temp_filename(clear=True)

        tmpf1, tmpf2 = tmpf + "_a.bed.gz", tmpf + "_b.bed.gz"
        stmnt = standardise_bed_files(tmpf1, tmpf2, infile,
                                      params.reference_bed)

        retval = P.run(
            "{stmnt}; "
            "{params.path} jaccard "
            "-a {tmpf1} -b {tmpf2} "
            "{params.options} "
            "2> {outfile}.log "
            ">> {outfile}; "
            "rm -f {tmpf1} {tmpf2}"
            .format(**locals()))

        return retval


class run_metric_bedtools_summarise(MetricRunnerBedtools):
    """run bedtools groupby on a :term:`bed` file.

    :param infile: Input file in :term:`bed` format
    :param outfile: Output file in :term:`tsv` format

    *Columns*

    """

    name = "bedtools_summarise"

    def run(self, infile, outfile, params):

        with IOTools.open_file(outfile, "w") as outf:
            outf.write(
                "contig\tcount\tsum\tmin\tmax\tmean\t"
                "median\tstddev\tcollapse\n")

        retval = P.run(
            "zcat {infile} "
            "| awk '{{printf(\"%%s\\t%%i\\n\", $1, $3-$2); "
            " printf(\"total\\t%%i\\n\", $3-$2)}}' "
            "| sort -k1,1 "
            "| {params.path} groupby "
            "-g 1 "
            "-c 2 "
            "-o count,sum,min,max,mean,median,stddev,collapse "
            "{params.options} "
            "2> {outfile}.log "
            ">> {outfile}; "
            .format(**locals()))

        return retval


class run_metric_bedtools_intersection(MetricRunnerBedtools):
    """compute intersection between two bed files outputting
    shared and unique intervals.

    *Columns*

    category
       Category of counts (shared|unique_test|unique_truth)
    count
       Number of intervals

    """

    name = "bedtools_intersect"

    sections = ["shared", "unique_truth", "unique_test"]

    reference_bed = None

    def run(self, infile, outfile, params):

        if params.reference_bed is None:
            raise ValueError("{} requires reference_bed to be set".format(
                self.name))

        # requires a consistent sort order, so sort both files.
        # It also requires the chromosome content to be identical,
        # so restrict output to common sets.
        tmpf = P.get_temp_filename(clear=True)

        tmpf_test, tmpf_truth = tmpf + "_a.bed.gz", tmpf + "_b.bed.gz"
        stmnt = standardise_bed_files(tmpf_test,
                                      tmpf_truth,
                                      infile,
                                      params.reference_bed)

        statements = [stmnt]
        statements.append(
            "{params.path} intersect "
            "-a {tmpf_test} "
            "-b {tmpf_truth} "
            "-wa "
            "| bgzip "
            "> {outfile}.shared.bed.gz"
        )
        statements.append(
            "{params.path} intersect "
            "-a {tmpf_test} "
            "-b {tmpf_truth} "
            "-wa -v"
            "| bgzip "
            "> {outfile}.unique_test.bed.gz"
        )
        statements.append(
            "{params.path} intersect "
            "-b {tmpf_test} "
            "-a {tmpf_truth} "
            "-wa -v"
            "| bgzip "
            "> {outfile}.unique_truth.bed.gz"
        )
        statements.append(
            "rm -f {tmpf_test} {tmpf_truth}")

        for section in self.sections:
            statements.append(
                "tabix -p bed {outfile}.{section}.bed.gz".format(**locals()))

        statement = "; ".join(statements)
        retval = P.run(statement.format(**locals()))

        # these are small files, so doing it here. Implement tabix.count()
        # method
        counts = dict()
        for section in self.sections:
            # with pysam.Tabixfile(outfile + "." + section + ".bed.gz") as inf:
            inf = pysam.Tabixfile(outfile + "." + section + ".bed.gz")
            counts[section] = len(list(inf.fetch()))
            inf.close()

        with IOTools.open_file(outfile, "w") as outf:
            outf.write("section\tcounts\n")
            outf.write("\n".join(["\t".join(
                map(str, x)) for x in list(counts.items())]) + "\n")

        return retval


class run_metric_bedtools_intersection_and_annotate(
        run_metric_bedtools_intersection):
    name = "bedtools_intersect_and_annotate"
    gat_path = "gat-run.py"

    tablenames = ["bedtools_intersect_and_annotate_counts",
                  "bedtools_intersect_and_annotate_enrichment"]

    annotations_bed = None
    workspace_bed = None

    def get_version(self):
        version_bedtools = run_metric_bedtools_intersection.get_version(self)
        help_string = E.run(
            "{self.gat_path} --version 2> /dev/null".format(**locals()),
            return_stdout=True).strip()
        return "{} {}".format(
            version_bedtools,
            re.search(r"gat-run.py version: (\S+):", help_string).groups()[0])

    def run(self, infile, outfile, params):

        if params.annotations_bed is None:
            raise ValueError("{} requires annotations_bed to be set".format(
                self.name))

        if params.workspace_bed is None:
            raise ValueError("{} requires workspace_bed to be set".format(
                self.name))

        retval = run_metric_bedtools_intersection.run(
            self, infile, outfile, params)
        retvals = [retval]

        statements = [
            "mv {outfile} {outfile}.bedtools_intersect_and_annotate_counts.tsv"
            .format(**locals())]
        bed_files = []
        for section in self.sections:
            tmpf = P.get_temp_filename(clear=True) + "-" + section + ".gz"
            statements.append(
                "zcat {outfile}.{section}.bed.gz "
                "| awk -v OFS='\\t' '{{ $4 = \"{section}\"; print }}' "
                "| bgzip > {tmpf}".format(**locals()))

            bed_files.append(tmpf)

        segment_files = " ".join(
            ["--segment-bed-file={}".format(x) for x in bed_files])

        statements.append(
            "{params.gat_path} "
            "{segment_files} "
            "--with-segment-tracks "
            "--annotation-bed-file={params.annotations_bed} "
            "--workspace-bed-file={params.workspace_bed} "
            "--log={outfile} "
            "{params.options} "
            "> {outfile}.bedtools_intersect_and_annotate_enrichment.tsv"
            .format(**locals()))

        for f in bed_files:
            statements.append("rm -f {}".format(f))

        statement = "; ".join(statements)
        retvals.append(P.run(statement))

        return retvals


class run_metric_gat_enrichment(MetricRunnerBedtools):
    """compute overlap and significance of overlap between
    two sets of intervals.

    *Columns*

    """
    path = "gat-run.py"
    name = "gat_enrichment"

    def get_version(self):
        help_string = E.run("{self.path} --version 2> /dev/null".format(
            **locals()), return_stdout=True).strip()
        return re.search("gat-run.py version: (\S+):", help_string).groups()[0]

    def run(self, infile, outfile, params):

        tmpf = P.get_temp_filename(clear=True)

        tmpf_test, tmpf_truth = tmpf + "_a.bed.gz", tmpf + "_b.bed.gz"
        stmnt = standardise_bed_files(tmpf_test,
                                      tmpf_truth,
                                      infile,
                                      params.annotations_bed)
        statements = [stmnt]
        statements.append(
            "{params.path} "
            "--segment-bed-file={tmpf_test} "
            "--ignore-segment-tracks "
            "--annotation-bed-file={tmpf_truth} "
            "--workspace-bed-file={params.workspace_bed} "
            "--log={outfile}.log "
            "{params.options} "
            "> {outfile}"
        )

        statement = "; ".join(statements)
        return P.run(statement.format(**locals()))


class SplitRunnerBedtools(SplitRunner):
    path = "bedtools"

    def get_version(self):
        help_string = E.run("{self.path} --version".format(**locals()),
                            return_stdout=True).strip()
        return re.search("bedtools (\S+)", help_string).groups()[0]


class run_split_bedtools_intersection(SplitRunnerBedtools):
    """compute intersection between two bed files outputting
    shared and unique intervals.

    """

    name = "bedtools_intersect"

    output = ["shared.dir/result.bed.gz",
              "unique_test.dir/result.bed.gz",
              "unique_truth.dir/result.bed.gz"]

    def run(self, infile, outfiles, params):

        # requires a consistent sort order, so sort both files.
        # It also requires the chromosome content to be identical,
        # so restrict output to common sets.
        tmpf = P.get_temp_filename(clear=True)

        outfile_shared, outfile_test, outfile_truth = outfiles

        tmpf_test, tmpf_truth = tmpf + "_a.bed.gz", tmpf + "_b.bed.gz"
        stmnt = standardise_bed_files(tmpf_test,
                                      tmpf_truth,
                                      infile,
                                      params.reference_bed)

        statements = [stmnt]
        statements.append(
            "{params.path} intersect "
            "-a {tmpf_test} "
            "-b {tmpf_truth} "
            "-wa "
            "| bgzip "
            "> {outfile_shared} "
        )
        statements.append(
            "{params.path} intersect "
            "-a {tmpf_test} "
            "-b {tmpf_truth} "
            "-wa -v"
            "| bgzip "
            "> {outfile_test}"
        )
        statements.append(
            "{params.path} intersect "
            "-b {tmpf_test} "
            "-a {tmpf_truth} "
            "-wa -v"
            "| bgzip "
            "> {outfile_truth}"
        )
        statements.append(
            "rm -f {tmpf_test} {tmpf_truth}")

        for f in outfiles:
            statements.append(
                "tabix -f -p bed {}".format(f))

        statement = "; ".join(statements)
        retval = P.run(statement.format(**locals()))

        return retval
