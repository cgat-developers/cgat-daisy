from .BAMMetrics import run_metric_bam_fastqc


class run_metric_fastq_fastqc(run_metric_bam_fastqc):
    name = "fastq_fastqc"
