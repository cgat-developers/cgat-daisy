import re
import os
import pysam
import shutil
import pandas

from .Runner import resolve_argument, is_true
from .ToolRunner import ToolRunner
from .VCFMetrics import VCFPreprocessor, restrict_bed, create_genome_bed
import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools


class ToolRunnerVCF(ToolRunner, VCFPreprocessor):
    pass


class run_tool_bcftools_intersect(ToolRunnerVCF):
    """run bcftools isec --complement

    :param infiles: Input files in :term:`vcf` format
    :param outfile: Output file in :term:`vcf` format

    This tool runs bcftools isec --complement to output
    the variants unique to each of the two files.
    """

    path = "bcftools"
    name = "bcftools_intersect"

    expected = ["vcf"]
    output = ["test_unique.dir/result.vcf.gz",
              "comp_unique.dir/result.vcf.gz",
              "test_shared.dir/result.vcf.gz",
              "comp_shared.dir/result.vcf.gz"]

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]

    def run(self, outfiles, params):

        vcf = resolve_argument(params.vcf, sep=" ")
        vcfs = vcf.split(" ")
        if len(vcfs) != 2:
            raise ValueError("expected 2 VCF files, received {}".format(vcfs))

        outdir = os.path.commonprefix(outfiles)
        outfile = os.path.join(outdir, "result.log")

        retval = P.run(
            "{params.path} isec "
            "{params.options} "
            "--output-type z "
            "--prefix {outdir} "
            "{vcf} "
            "&> {outfile} "
            .format(**locals()))

        f = ["000{}.vcf.gz".format(x) for x in range(4)]
        self.distribute_results(outdir, list(zip(f, self.output)))

        f = ["000{}.vcf.gz.tbi".format(x) for x in range(4)]
        ff = [x + ".tbi" for x in self.output]
        self.distribute_results(outdir, list(zip(f, ff)))

        return retval


class run_tool_bcftools_intersect_lookup(ToolRunnerVCF):
    """run bcftools isec --complement on two filtered
    files reporting the results in unfiltered output.

    :param infiles: Input files in :term:`vcf` format
    :param outfile: Output file in :term:`vcf` format

    This tool expects three VCF:
    1. A raw, unfiltered VCF file
    2. A filtered version of 1.
    3. A reference set of variants.

    This tool runs bcftools isec --complement on 2 and on 3
    and then outputs the variants at the corresponding positions
    from 1.
    """

    path = "bcftools"
    name = "bcftools_intersect_lookup"
    path_bedtools = "bedtools"

    expected = ["vcf"]
    output = ["test_unique.dir/result.vcf.gz",
              "comp_unique.dir/result.vcf.gz",
              "test_shared.dir/result.vcf.gz"]

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]

    def run(self, outfiles, params):

        vcf = resolve_argument(params.vcf, sep=" ")
        vcfs = vcf.split(" ")
        if len(vcfs) != 3:
            raise ValueError("expected 3 VCF files, received {}".format(vcfs))

        lookup_vcf = vcfs.pop(0)
        vcf = " ".join(vcfs)

        outdir = os.path.commonprefix(outfiles)
        outfile = os.path.join(outdir, "result.log")

        retvals = []
        retvals.append(P.run(
            "{params.path} isec "
            "{params.options} "
            "--output-type z "
            "--prefix {outdir} "
            "{vcf} "
            "&> {outfile} "
            .format(**locals())))

        retvals.extend(self.distribute_results(
            outdir,
            list(zip(["0000.vcf.gz", "0001.vcf.gz", "0002.vcf.gz"], self.output)),
            statement="{params.path_bedtools} intersect "
            "-a {lookup_vcf} "
            "-b {{infile}} "
            "-wa "
            "-header "
            "| bgzip > {{outfile}}".format(**locals())))

        return retvals


class run_tool_annotate_with_truth(ToolRunnerVCF):
    """run bcftools isec --complement on two filtered
    files reporting the results in unfiltered output.

    :param infiles: Input files in :term:`vcf` format
    :param outfile: Output file in :term:`vcf` format

    This tool expects three VCF:
    1. A raw, unfiltered VCF file (test)
    2. A set of three files from the test file : FP/FN/TP
    3. A set of three files from a reference file: FP/FN/TP

    This tool will produce a single VCF with an additional
    tag in the INFO field called AS (for assessment). Its
    values are
    * TPS: Shared true positive
    * FNU: Unique false negative
    * FNS: Shared false negative
    * FNO: False negative unique to other set
    * FPU: Unique false positive
    * FPS: Shared false positive
    * FPU: False positive unique to other set

    Variants with missing quality code are typically true
    negatives and appear when calling in a trio and the
    particular variant is not present in the sample under
    consideration.

    """

    path = "bcftools"
    name = "annotate_with_truth"
    path_bedtools = "bedtools"

    expected = ["vcf"]
    output = "result.vcf.gz"

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]

    def run(self, outfile, params):

        try:
            vcf_target = params.vcf["target"]
            test_fp = params.vcf["test"]["fp"]
            test_fn = params.vcf["test"]["fn"]
            test_tp = params.vcf["test"]["tp"]
            comp_fp = params.vcf["compare"]["fp"]
            comp_fn = params.vcf["compare"]["fn"]
            comp_tp = params.vcf["compare"]["tp"]
        except KeyError as msg:
            raise ValueError(
                "missing input data: {}".format(msg))

        tmpdir = P.get_temp_filename(clear=True)

        outdir = os.path.dirname(outfile)
        bedfile = os.path.join(tmpdir, "annotations.bed.gz")
        bedfile_sorted = os.path.join(outdir, "annotations.bed.gz")

        header = os.path.join(outdir, "header.txt")
        with open(header, "w") as outf:
            outf.write(
                '##INFO=<ID=AS,Number=.,Type=String,'
                'Description="Assessment code. Combination of FP/FN/TP and '
                'U for unique, O for other and S for shared.">')

        statements = ["mkdir {tmpdir}".format(**locals())]
        toprocess = []
        for a, b, label in zip(
                (test_fp, test_fn, test_tp),
                (comp_fp, comp_fn, comp_tp),
                ("FP", "FN", "TP")):
            statements.append(
                "{params.path} isec "
                "--output-type z "
                "--prefix {tmpdir}/{label} "
                "{a} {b}"
                "&> {outfile}.isec_{label}.log "
                .format(**locals()))
            toprocess.append(
                (os.path.join(tmpdir, label, "0000.vcf.gz"),
                 label + "U"))
            toprocess.append(
                (os.path.join(tmpdir, label, "0001.vcf.gz"),
                 label + "O"))
            toprocess.append(
                (os.path.join(tmpdir, label, "0002.vcf.gz"),
                 label + "S"))

        # TPO = FNU
        # TPU = FNO
        toprocess = [x for x in toprocess if x[1] not in ["TPO", "TPU"]]
        # files to keep, these are variants that will be not in the vcf
        # file that is being annotated.

        keep = ["FNS", "FNU", "FNO"]
        for f, label in toprocess:
            statements.append(
                "zcat {f} "
                "| awk '!/^#/ "
                "{{printf(\"%%s\\t%%i\\t%%i\\t{label}\\n\", $1, $2-1, $2) }}'"
                "| bgzip "
                ">> {bedfile} "
                .format(**locals()))
            if label in keep:
                statements.append(
                    "cp {f} {outfile}.{label}.vcf.gz".format(**locals()))

        statements.append(
            "zcat {bedfile} "
            "| sort -k1,1 -k2,2n "
            "| bedtools merge -i stdin -c 4 -o distinct -delim ',' "
            "2> {bedfile_sorted}.log "
            "| bgzip "
            "> {bedfile_sorted}".format(**locals()))
        statements.append(
            "tabix -p bed {bedfile_sorted}".format(**locals()))

        statements.append(
            "bcftools annotate "
            "--annotations={bedfile_sorted} "
            "--columns=CHROM,FROM,TO,AS "
            "--header-lines {header} "
            "--output-type z "
            "{vcf_target} "
            "2> {outfile}.log "
            "> {outfile}; "
            "tabix -p vcf {outfile} ".format(**locals()))

        statements.append("rm -rf {tmpdir}".format(**locals()))

        statement = "; ".join(statements)

        return self.run_with_preprocessing(
            vcf_target, outfile, params,
            statement)


class run_tool_import_vcf(ToolRunnerVCF):

    name = "import_vcf"

    expected = ["vcf"]
    output = "result.vcf.gz"

    def get_version(self):
        return "builtin"

    def run(self, outfile, params):

        return self.run_with_preprocessing(
            params.vcf,
            outfile,
            params,
            "cp {params.vcf} {outfile}; "
            "tabix -p vcf {outfile} >& {outfile}.tabix.log".format(**locals()))


class run_tool_bcftools_subset_and_intersect(ToolRunnerVCF):
    """run bcftools view -S and then bcftools isec

    :param infiles: Input files in :term:`vcf` format
    :param outfile: Output file in :term:`vcf` format

    This tool expects three VCFs:
    1. A raw, unfiltered VCF file ("primary_vcf")
    2. A reference VCF (with samples, "filter_vcf")

    The first VCF is subset down to the samples and sites in
    the second VCF
    """

    path = "bcftools"
    name = "bcftools_subset_and_intersect"

    expected = ["primary_vcf", "filter_vcf"]
    output = "result.vcf.gz"

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]

    def run(self, outfile, params):

        outfile = os.path.abspath(outfile)
        if params.primary_vcf is None:
            raise ValueError("expected primary_vcf, received {}".format(params.primary_vcf))
        if params.filter_vcf is None:
            raise ValueError("expected filter_vcf, received {}".format(params.filter_vcf))

        primary_vcf = os.path.abspath(params.primary_vcf)
        filter_vcf = os.path.abspath(params.filter_vcf)
        outdir = os.path.dirname(outfile)

        retval = P.run(
            "( "
            "cd {outdir} && "
            "{params.path} query -l {filter_vcf} > subset_samples "
            "&& {params.path} view {params.options} --force-samples "
            "-S subset_samples {primary_vcf} -Ob -o test.subset_samples.bcf "
            "&& {params.path} index test.subset_samples.bcf "
            "&& {params.path} isec {params.options} -n=2 "
            "--prefix isec "
            "test.subset_samples.bcf {filter_vcf} "
            "--output-type z "
            "&& mv -f isec/0000.vcf.gz {outfile} "
            "&& tabix {outfile} "
            ") &> {outfile}.log "
            .format(**locals()))

        return retval


class run_tool_bcftools_filter(ToolRunnerVCF):
    """run bcftools filter

    :param infiles: Input files in :term:`vcf` format
    :param outfile: Output file in :term:`vcf` format

    This tool runs bcftools filter {options} to filter input files.
    """

    path = "bcftools"
    name = "bcftools_filter"

    expected = ["vcf"]
    output = "result.vcf.gz"

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]

    def run(self, outfile, params):

        if params.vcf is None:
            raise ValueError("expected vcf, received {}".format(params.vcf))

        vcf = resolve_argument(params.vcf, sep=" ")
        outfile = os.path.abspath(outfile)
        outdir = os.path.dirname(outfile)

        retval = P.run(
            "{params.path} filter {params.options} {vcf} -Oz -o {outfile} "
            "&> {outfile}.log"
            .format(**locals()))

        return retval


class run_tool_bcftools_merge(ToolRunnerVCF):
    """run bcftools filter

    :param infiles: Input files in :term:`vcf` format
    :param outfile: Output file in :term:`vcf` format

    This tool runs `bcftools merge {options}` to merge input files.
    """

    path = "bcftools"
    name = "bcftools_merge"

    expected = ["vcfs"]
    output = "result.vcf.gz"

    def get_version(self):
        help_string = E.run("{self.path}".format(**locals()),
                            return_stderr=True).strip()
        return re.search("Version: (\S+)", help_string).groups()[0]

    def run(self, outfile, params):

        if params.vcfs is None:
            raise ValueError("expected vcfs, received {}".format(params.vcf))

        vcfs = resolve_argument(params.vcfs, sep=" ")
        outfile = os.path.abspath(outfile)

        retval = P.run(
            "{params.path} merge {params.options} {vcfs} -Oz -o {outfile} "
            "&> {outfile}.log"
            .format(**locals()))

        return retval
