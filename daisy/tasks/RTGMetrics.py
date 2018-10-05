import re
import os
import shutil
import pandas
from distutils.version import LooseVersion
import pysam
import cgatcore.experiment as E
import cgatcore.iotools as IOTools
from daisy.toolkit import update_namedtuple
from daisy.tasks.VCFMetrics import MetricRunnerVCF, restrict_bed


class MetricRunnerVCFRTG(MetricRunnerVCF):
    path = "rtg"

    def get_version(self):
        help_string = E.run("{self.path} version 2> /dev/null".format(**locals()),
                            return_stdout=True,
                            on_error="ignore").strip()
        if help_string and "not found" not in help_string:
            return re.search(r"Product: RTG Tools (\S+)", help_string).groups()[0]
        else:
            raise ValueError("rtg not found at/as {}: {}".format(
                self.path, help_string))


class run_metric_rtg_vcfeval(MetricRunnerVCFRTG):
    """run vcfeval tools by RTG on a :term:`vcf` file.

    :param infile: Input file in :term:`vcf` format
    :param outfile: Multiple output files in compressed :term:`tsv`
         format

    reference_sdf : string
        Filename with reference sequence in SDF format.

    reference_vcf: string
        Filename with reference variant calls in VCF format.

    callable_bed: string
        Filename with callable regions in BED format.

    """

    name = "rtg_vcfeval"

    callable_bed = None
    reference_vcf = None
    reference_sdf = None

    def run(self, infile, outfile, params):

        if params.reference_vcf is None:
            raise ValueError("missing input parameter 'reference_vcf'")
        if params.reference_sdf is None:
            raise ValueError("missing input parameter 'reference_sdf'")
        if params.callable_bed is None:
            raise ValueError("missing input parameter 'callable_bed'")

        outfile_regions = outfile + ".bed.gz"
        restrict_bed(outfile_regions,
                     params.callable_bed,
                     infile,
                     remove_chr=params.remove_chr,
                     add_chr=params.add_chr)

        outputdir = os.path.join(os.path.dirname(outfile),
                                 "vcfeval.dir")
        if os.path.exists(outputdir):
            shutil.rmtree(outputdir)

        if LooseVersion(self.get_version()) < LooseVersion("3.7"):
            bed_options = "--bed-regions={}".format(params.callable_bed)
            output_columns = ["threshold",
                              "true_positive_count",
                              "false_positive_count",
                              "false_negative_count",
                              "false_discovery_rate",
                              "false_negative_rate",
                              "f_measure"]
        else:
            bed_options = "--evaluation-regions={}".format(params.callable_bed)
            output_columns = ["threshold",
                              "true_positive_baseline",
                              "true_positive_count",
                              "false_positive_count",
                              "false_negative_count",
                              "false_discovery_rate",
                              "false_negative_rate",
                              "f_measure"]

        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{params.path} vcfeval "
            "--calls={infile} "
            "--baseline={params.reference_vcf} "
            "--template={params.reference_sdf} "
            "{bed_options} "
            "--output={outputdir} "
            "{params.options} "
            ">& {outfile}.log ".format(**locals()),
            job_memory="unlimited")

        with IOTools.open_file(os.path.join(outputdir, "summary.txt")) as inf:
            with IOTools.open_file(outfile, "w") as outf:
                table = []
                for line in inf:
                    if line.startswith("-"):
                        continue
                    line = re.sub("^ +", "", line)
                    line = re.sub(" +", "\t", line)
                    fields = line[:-1].split("\t")
                    table.append(fields)

                df = pandas.DataFrame(table[1:], columns=table[0])
                df.columns = output_columns
                # convert precision and sensitivity
                df["false_discovery_rate"] = 1.0 - df["false_discovery_rate"].astype(float)
                df["false_negative_rate"] = 1.0 - df["false_negative_rate"].astype(float)
                df.to_csv(outf, sep="\t", index=False)

        return retval


class run_metric_rtg_vcfdiff(MetricRunnerVCFRTG):
    """run vcfdiff tools by RTG on a :term:`vcf` file.

    :param infile: Input file in :term:`vcf` format
    :param outfile: Multiple output files in compressed :term:`tsv`
         format

    reference_sdf : string
        Filename with reference sequence in SDF format.

    reference_vcf: string
        Filename with reference variant calls in VCF format.

    callable_bed: string
        Filename with callable regions in BED format.

    reference_fasta: string
        Fasta file with genomic reference. Used if no
        callable_bed is given to use the full genome
        in the comparison.

    """

    name = "rtg_vcfdiff"

    callable_bed = None
    reference_vcf = None
    reference_sdf = None
    reference_fasta = None

    def run(self, infile, outfile, params):

        if params.reference_vcf is None:
            raise ValueError("missing input parameter 'reference_vcf'")
        if params.reference_sdf is None:
            raise ValueError("missing input parameter 'reference_sdf'")
        if params.callable_bed is None and params.reference_fasta is None:
            raise ValueError(
                "missing input parameter: either 'callable_bed' or "
                "'reference_fasta' is needed")

        outfile_regions = outfile + ".bed.gz"

        if "callable_bed" in params._fields:
            restrict_bed(outfile_regions,
                         params.callable_bed,
                         infile,
                         remove_chr=params.remove_chr,
                         add_chr=params.add_chr)
        else:
            create_genome_bed(outfile_regions,
                              infile,
                              params.reference_fasta,
                              remove_chr=params.remove_chr,
                              add_chr=params.add_chr)

        with pysam.VariantFile(params.reference_vcf.strip()) as inf:
            try:
                # in some pathological VCF (multiple headers), sample
                # names are not properly read in pysam.
                sample_name = list(inf.header.samples)[0]
            except IndexError:
                sample_name = "TOCOMPARE"

            params = update_namedtuple(params, rename_samples=sample_name)

        outfile_reference = outfile + ".ref.vcf.gz"
        preprocess_reference = self.build_statement_with_preprocessing(
            params.reference_vcf,
            outfile_reference,
            params,
            "mv {params.reference_vcf} {outfile_reference}; "
            "tabix -f -p vcf {outfile_reference}".format(
                **locals()))

        outputdir = os.path.join(os.path.dirname(outfile),
                                 "vcfeval.dir")

        if os.path.exists(outputdir):
            shutil.rmtree(outputdir)

        # The java VM does not work with the ulimit -v and ulimit -h
        # options.
        retval = self.run_with_preprocessing(
            infile,
            outfile,
            params,
            "{preprocess_reference}; "
            "{params.path} vcfeval "
            "--calls={infile} "
            "--baseline={outfile_reference} "
            "--template={params.reference_sdf} "
            "--bed-regions={outfile_regions} "
            "--output={outputdir} "
            "--sample={sample_name} "
            ">& {outfile}.log; "
            "rm -f {outfile_reference} {outfile_reference}.tbi ".format(**locals()),
            job_memory="unlimited",
        )

        with IOTools.open_file(os.path.join(outputdir, "summary.txt")) as inf:
            with IOTools.open_file(outfile, "w") as outf:
                table = []
                for line in inf:
                    if line.startswith("-"):
                        continue
                    line = re.sub("^ +", "", line)
                    line = re.sub(" +", "\t", line)
                    fields = line[:-1].split("\t")
                    table.append(fields)

                df = pandas.DataFrame(table[1:], columns=table[0])
                df.columns = ["threshold",
                              "true_positive_count",
                              "false_positive_count",
                              "false_negative_count",
                              "false_discovery_rate",
                              "false_negative_rate",
                              "f_measure"]
                # convert precision and sensitivity
                df["false_discovery_rate"] = 1.0 - df["false_discovery_rate"].astype(float)
                df["false_negative_rate"] = 1.0 - df["false_negative_rate"].astype(float)
                df.to_csv(outf, sep="\t", index=False)

        return retval
