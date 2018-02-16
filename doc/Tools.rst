===========
Daisy tools
===========

The daisy framework includes a variety of tools and utility
scripts. These tools are available through the :file:`daisy`
command line interface. To list all available tools, type::

   daisy

To run a particular tools, type::

   daisy <toolname>

for example, to get help on the `bam2stats` tool, type::

   daisy bam2stats --help

Both daisy-chains and simple command line utilities are
available through this interface.

Workflow tools
==================

.. toctree::
   :maxdepth: 1

   tools/benchmark-rename-files.rst
   tools/benchmark-simple.rst
   tools/benchmark-upload.rst
   tools/check-task-library.rst
   tools/merge-databases.rst
   tools/run-etm.rst
   tools/run-task.rst
   tools/run.rst
   tools/test-task-library.rst
   tools/profile-chain.rst
   tools/watch-and-run.rst

Genomics tools
==============

Tools in this section manipulate and compute metrics from genomic data
sets.

.. toctree::
   :maxdepth: 1
	      
   tools/bam-compare-alignments.rst
   tools/bam-pileup2tsv.rst
   tools/bam2bam-split-reads.rst
   tools/bam2bam.rst
   tools/bam2depth.rst
   tools/bam2stats.rst
   tools/bam2window-stats.rst
   tools/bed-vs-bed.rst
   tools/bed2plot.rst
   tools/fasta2fasta.rst
   tools/fasta2stats.rst
   tools/fasta2vcf.rst
   tools/fastq2fasta.rst
   tools/fastq2fastq.rst
   tools/fastq2tsv.rst
   tools/maf2maf.rst
   tools/table2stats.rst
   tools/vcf-assessment.rst
   tools/vcf-compare-phase.rst
   tools/vcf-stats.rst
   tools/vcf-vs-vcf.rst
   tools/vcf2tsv.rst
   tools/vcf2vcf.rst
   tools/vcfqc-report.rst
   tools/version.rst

