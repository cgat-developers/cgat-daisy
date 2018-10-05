===================
Daisy Documentation
===================

.. image:: daisies.png

Daisy is a framework to perform computational experiments efficiently,
reproducibly, and at scale. An experiment is defined by an
experimental design in :term:`yaml` format that describes one or more
tools to be run on one or more data sets and collecting on or more
metrics from the results.

At its simplest, an experimental design would look like this::

   setup:
	tools:
      	- freebayes
      	- weCall
   	metrics:
      	- bcftools_stats
   input:
      bam: *.bam
      reference_fasta: hs37d5.fa

This would apply the variant callers `freebayes` and `weCall` to all
files ending in `.bam` in the current directory and running the tool
`bctools stats` on the results. This is a typical benchmarking
scenario comparing to variant callers.

The same experimental design can be used for parameter optimization,
for example::

   setup:
	tools:
      	- freebayes
   	metrics:
      	- bcftools_stats
   input:
      bam: *.bam
      reference_fasta: hs37d5.fa
   freebayes:
      options:
	- --haplotype-length 50
	- --haplotype-length 100

This design would run freebayes with two different options.

While originally developed for benchmarking, we have found the same
framework useful in processing large data sets. For example, to map a
multiple read data sets from a high-throughput sequencing experiment,
you could say::

   setup:
	tools:
      	- bwa_mem
   	metrics:
      	- samtools_stats
	
   input:
      fastq: *.fastq.gz
      reference_fasta: hs37d5.fa
      group_regex: /([^/]+)-R.*.fastq.gz
      group_alias: \1

The regular expression ensures that the two components of a sample
(named `<sample>-R1.fastq.gz` and `<sample>-R2.fastq.gz` are grouped
and supplied together to the mapping tool `bwa mem`.  The `samtools
stats` command is run to provide QC metrics such as the proportion of
reads mapped to the reference genome.

Contents
========

.. toctree::
   :maxdepth: 2

   Installation.rst
   Usage.rst
   Reference.rst
   Development.rst
   Tutorials.rst
   Glossary.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

