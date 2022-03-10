import re
import os

from .ToolRunner import ToolRunner
from .Runner import resolve_argument
import cgatcore.pipeline as P
import cgatcore.experiment as E
import cgatcore.iotools as IOTools
from .VCFMetrics import VCFPreprocessor
from .VariantCallers import get_reference


class GATKRunner(ToolRunner):
    path = "/data/install/Licensed/gatk/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar"

    def get_version(self):
        return E.run("java -jar {self.path} --version".format(**locals()),
                     return_stdout=True).strip()

    def run_statements(self, stmnts, **kwargs):

        stmnts = [x for x in stmnts if x]

        filename, main_statement, post_statement = P.join_statements(
            stmnts, infile=None)

        stmnt = " ; ".join([x for x in [main_statement, post_statement] if x])

        job_threads1, job_threads2 = 1, 1
        if "--num_threads" in stmnt:
            try:
                job_threads1 = max(
                    [int(x) for x in re.search("--num_threads\s*(\d+)",
                                               stmnt).groups()])
            except AttributeError:
                pass

        if "--num_cpu_threads_per_data" in stmnt:
            try:
                job_threads2 = max(
                    [int(x) for x in
                     re.search("--num_cpu_threads_per_data_thread\s*(\d+)",
                               stmnt).groups()])
            except AttributeError:
                pass

        job_threads = max(job_threads1, job_threads2)
        return P.run(stmnt, **kwargs)

    def build_calibration_workflow(self, outfile, prefix, vcfinput, params):

        try:
            variant_mode = params.variant_mode
        except AttributeError:
            variant_mode = "all"

        stmnts = []
        do_snps = variant_mode == "all" or variant_mode == "snps"
        do_indels = variant_mode == "all" or variant_mode == "indels"
        reference_fasta = get_reference(params)

        if do_snps:
            stmnts.append(
                "java "
                "-Djava.io.tmpdir=%(tmpdir)s "
                "-jar {self.path} "
                "--analysis_type VariantRecalibrator "
                "--input {vcfinput} "
                "--reference_sequence {reference_fasta} "
                "--recal_file {prefix}.snp.recal "
                "--tranches_file {prefix}.snp.tranches "
                "--mode SNP "
                "--logging_level INFO "
                "--log_to_file {prefix}.VariantRecalibratorSNP.log "
                "{params.variantrecalibrator_snp} "
                ">& {prefix}.VariantRecalibratorSNP.err".format(**locals()))

        if do_indels:
            stmnts.append(
                "java "
                "-Djava.io.tmpdir=%(tmpdir)s "
                "-jar {self.path} "
                "--analysis_type VariantRecalibrator "
                "--input {vcfinput} "
                "--reference_sequence {reference_fasta} "
                "--recal_file {prefix}.indel.recal "
                "--tranches_file {prefix}.indel.tranches "
                "--mode INDEL "
                "--logging_level INFO "
                "--log_to_file {prefix}.VariantRecalibratorIndel.log "
                "{params.variantrecalibrator_indel} "
                ">& {prefix}.VariantRecalibratorIndel.err".format(**locals()))

        if do_snps:
            if do_indels:
                output = prefix + ".snp_recal.vcf.gz"
            else:
                output = outfile

            stmnts.append(
                "java "
                "-Djava.io.tmpdir=%(tmpdir)s "
                "-jar {self.path} "
                "--analysis_type ApplyRecalibration "
                "--input {vcfinput} "
                "--reference_sequence {reference_fasta} "
                "--recal_file {prefix}.snp.recal "
                "--tranches_file {prefix}.snp.tranches "
                "--mode SNP "
                "--out {output} "
                "--log_to_file {prefix}.ApplyRecalibrationSNP.log "
                "{params.applyrecalibration_snp} "
                ">& {prefix}.ApplyRecalibrationSNP.err".format(**locals()))

            if do_indels:
                vcfinput = output

        if do_indels:
            stmnts.append(
                "java "
                "-Djava.io.tmpdir=%(tmpdir)s "
                "-jar {self.path} "
                "--analysis_type ApplyRecalibration "
                "--input {vcfinput} "
                "--reference_sequence {reference_fasta} "
                "--recal_file {prefix}.indel.recal "
                "--tranches_file {prefix}.indel.tranches "
                "--mode INDEL "
                "--out {outfile} "
                "--log_to_file {prefix}.ApplyRecalibrationIndel.log "
                "{params.applyrecalibration_indel} "
                ">& {prefix}.ApplyRecalibrationIndel.err".format(**locals()))

        return stmnts


class run_tool_bam_GATKPreprocess(GATKRunner):

    expected = ["reference_fasta", "bam"]
    output = "result.bam"

    name = "gatk_preprocess"

    print_reads = ""
    indelrealigner = ""
    realignertargetcreator = ""

    def run(self, outfile, params):

        bam = resolve_argument(params.bam)
        reference_fasta = get_reference(params)
        stmnts = []
        prefix = IOTools.snip(outfile, ".bam")

        stmnts.append(
            "java "
            "-Djava.io.tmpdir=%(tmpdir)s "
            "-jar {self.path} "
            "--analysis_type RealignerTargetCreator "
            "--input_file {bam} "
            "--reference_sequence {reference_fasta} "
            "--logging_level INFO "
            "--log_to_file {outfile}.RealignerTargetCreator.log "
            "{params.realignertargetcreator} "
            "--out {outfile}.realign.intervals "
            ">& {outfile}.RealignerTargetCreator.err".format(**locals()))

        stmnts.append(
            "java "
            "-Djava.io.tmpdir=%(tmpdir)s "
            "-jar {self.path} "
            "--analysis_type IndelRealigner "
            "--input_file {bam} "
            "--reference_sequence {reference_fasta} "
            "--targetIntervals {outfile}.realign.intervals "
            "--logging_level INFO "
            "--log_to_file {outfile}.IndelRealigner.log "
            "{params.indelrealigner} "
            "--out @OUT@.bam "
            ">& {outfile}.IndelRealigner.err".format(**locals()))

        stmnts.append(
            "java "
            "-Djava.io.tmpdir=%(tmpdir)s "
            "-jar {self.path} "
            "--analysis_type BaseRecalibrator "
            "--input_file @IN@.bam "
            "--reference_sequence {reference_fasta} "
            "--logging_level INFO "
            "{params.baserecalibrator} "
            "--log_to_file {outfile}.BaseRecalibrator.log "
            "--out {outfile}.recal_data.table "
            ">& {outfile}.BaseRecalibrator.err".format(**locals()))

        stmnts.append(
            "java "
            "-Djava.io.tmpdir=%(tmpdir)s "
            "-jar {self.path} "
            "--analysis_type PrintReads "
            "--input_file @IN-1@.bam "
            "--reference_sequence {reference_fasta} "
            "--BQSR {outfile}.recal_data.table "
            "--logging_level INFO "
            "--log_to_file {outfile}.PrintReads.log "
            "--out {outfile} "
            ">& {outfile}.PrintReads.err".format(**locals()))

        stmnts.append(
            "mv {prefix}.bai {outfile}.bam.bai")

        return self.run_statements(stmnts, job_memory="3G")


class run_tool_GATKHaplotypeCaller(GATKRunner):
    """one-step calling workflow using the GATK HaplotypeCaller

    Variants are called with the GATK HaplotypeCaller. SNPs
    and indel recalibration metrics are computed individually
    and then applied in succession.

    Note that the chromosome names in the GATK resources are 1, 2, 3,
    thus make sure the BAM input files follow the same naming.

    The full workflow is:

    HaplotypeCaller
    VariantRecalibrator --mode=SNP
    VariantRecalibrator --mode=INDEL
    ApplyRelibration --mode=SNP
    ApplyRelibration --mode=INDEL

    """
    name = "gatk_haplotype_caller"

    expected = ["reference_fasta", "bam"]
    output = "result.vcf.gz"

    haplotypecaller = ""
    variantrecalibrator_snp = ""
    variantrecalibrator_indel = ""
    applyrecalibration_snp = ""
    applyrecalibration_indel = ""

    def run(self, outfile, params):

        bam = resolve_argument(params.bam)
        reference_fasta = get_reference(params)

        stmnts = []

        prefix = IOTools.snip(outfile, ".vcf.gz")
        vcf_output = prefix + ".raw.vcf.gz"

        if not os.path.exists(vcf_output):
            stmnts.append(
                "java "
                "-Djava.io.tmpdir=%(tmpdir)s "
                "-jar {self.path} "
                "--analysis_type HaplotypeCaller "
                "--input_file {bam} "
                "--reference_sequence {reference_fasta} "
                "--logging_level INFO "
                "--log_to_file {outfile}.HaplotypeCaller.log "
                "{params.haplotypecaller} "
                "--out {vcf_output} "
                ">& {prefix}.HaplotypeCaller.err".format(**locals()))
        else:
            E.warn("output file {vcf_output} already exists - "
                   "it will not be recomputed".format(**locals()))

        stmnts.extend(self.build_calibration_workflow(
            outfile, prefix, vcf_output, params))

        return self.run_statements(stmnts, job_memory="5G")


class run_tool_GATKHaplotypeJointCaller(GATKRunner):
    """one-step calling workflow using the GATK HaplotypeCaller

    Variants are called with the GATK HaplotypeCaller. SNPs
    and indel recalibration metrics are computed individually
    and then applied in succession.

    Note that the chromosome names in the GATK resources are 1, 2, 3,
    thus make sure the BAM input files follow the same naming.

    The full workflow is:

    HaplotypeCaller
    GenotypeGVCFs
    VariantRecalibrator --mode=SNP
    VariantRecalibrator --mode=INDEL
    ApplyRelibration --mode=SNP
    ApplyRelibration --mode=INDEL

    """
    name = "gatk_haplotype_joint_caller"

    expected = ["reference_fasta", "bam"]
    output = "result.vcf.gz"

    haplotypecaller = ""
    genotypegvcfs = ""
    variantrecalibrator_snp = ""
    variantrecalibrator_indel = ""
    applyrecalibration_snp = ""
    applyrecalibration_indel = ""

    def run(self, outfile, params):

        prefix = IOTools.snip(outfile, ".vcf.gz")
        bams = resolve_argument(params.bam, ",")
        reference_fasta = get_reference(params)

        statements, gvcfs = [], []
        # TODO: sort out multi-threading
        for idx, bam in enumerate(bams.split(",")):
            output = prefix + "." + str(idx) + ".g.vcf"
            gvcfs.append(output)

            if os.path.exists(output):
                E.info("{} already exists - skipped".format(output))
                continue

            statements.append(
                "java "
                "-Djava.io.tmpdir=%(tmpdir)s "
                "-jar {self.path} "
                "--analysis_type HaplotypeCaller "
                "--input_file {bam} "
                "--reference_sequence {reference_fasta} "
                "--emitRefConfidence GVCF "
                "--logging_level INFO "
                "--log_to_file {prefix}.HaplotypeCaller.{idx}.log "
                "{params.haplotypecaller} "
                "--out {output} "
                ">& {prefix}.HaplotypeCaller.{idx}.err".format(**locals()))

        if statements:
            self.run_statements(statements, job_memory="4G")

        stmnts = []
        gvcfs = " ".join(["--variant {}".format(x) for x in gvcfs])
        vcf_output = prefix + ".raw.vcf.gz"
        stmnts.append(
            "java "
            "-Djava.io.tmpdir=%(tmpdir)s "
            "-jar {self.path} "
            "--analysis_type GenotypeGVCFs "
            "--reference_sequence {reference_fasta} "
            "{gvcfs} "
            "--logging_level INFO "
            "--log_to_file {prefix}.GenotypeGVCFs.log "
            "{params.genotypegvcfs} "
            "--out {vcf_output} "
            ">& {prefix}.GenotypeGVCFs".format(**locals()))

        stmnts.extend(self.build_calibration_workflow(
            outfile, prefix, vcf_output, params))

        return self.run_statements(stmnts, job_memory="4G")


class run_tool_GATKVQSRFilter(GATKRunner, VCFPreprocessor):
    """apply VSQR filter to a VCF file

    Before running the VQSR filter, the VCF fill will be
    annotated with the GATK/VariantAnnotator tool.

    Note that GATK VariantAnnotator does not have the ability
    to annotate indels with the required metrics for VQSR, see here:

    http://gatkforums.broadinstitute.org/discussion/5177/sor-annotation-for-indels#

    The recommendation is to genotype with the GATK haplotype caller.

    Any values in the FILTER column will be set to PASS before running
    it through VQSR.
    """
    name = "gatk_vqsr"

    expected = ["vcf", "bam", "reference_fasta"]
    output = "result.vcf.gz"

    variantrecalibrator_snp = ""
    variantrecalibrator_indel = ""
    applyrecalibration_snp = ""
    applyrecalibration_indel = ""

    # valid values are: snps, indels, all
    variant_mode = "snps"

    def run(self, outfile, params):

        prefix = IOTools.snip(outfile, ".vcf.gz")

        bam = resolve_argument(params.bam, sep=",")
        reference_fasta = get_reference(params)

        bam = " ".join(["--input_file {}".format(x) for x in bam.split(",")])
        stmnts = []
        if not os.path.exists(prefix + ".annotated.vcf.gz"):
            tmpfile, pre_statement, post_statement = self.pre_process(
                params.vcf, outfile, params)

            stmnts.append(pre_statement)
            stmnts.append(
                "java "
                "-Djava.io.tmpdir=%(tmpdir)s "
                "-jar {self.path} "
                "--analysis_type VariantAnnotator "
                "--variant {tmpfile} "
                "{bam} "
                "--reference_sequence {reference_fasta} "
                "--logging_level INFO "
                "--log_to_file {prefix}.VariantAnnotator.log "
                "--annotation FisherStrand "
                "--annotation StrandOddsRatio "
                "--annotation ReadPosRankSumTest "
                "--annotation RMSMappingQuality "
                "--annotation MappingQualityRankSumTest "
                "{params.options} "
                "--out {prefix}.annotated.vcf.gz "
                ">& {prefix}.VariantAnnotator.err".format(**locals()))

            stmnts.extend(self.build_calibration_workflow(
                outfile, prefix, prefix + ".annotated.vcf.gz", params))

            stmnts.append(post_statement)
        else:
            E.warn("using pre-existing file {} with annotated variants".format(
                prefix + ".annotated.vcf.gz"))

            stmnts.extend(self.build_calibration_workflow(
                outfile, prefix, prefix + ".annotated.vcf.gz", params))

        return self.run_statements(stmnts, job_memory="3G")
