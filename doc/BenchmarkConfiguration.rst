.. _configuration:

=============
Configuration
=============

Benchmark reads its configuration from configuration files in
:term:`yaml` format. The following configuration files are read:

1. A default configuration file within the package which sets 
   default values for a variety of options.

2. A global configuration file in the user's home directory called
   :file:`~/.genomics_benchmark.yml`. This is a good place to put
   site-wide options for example to control which cluster queue to
   use.

3. A specific configuration file in the current directory called
   :file:`benchmark.yml` which describes the workflow.

A configuration setting in a later file overrides an earlier one.

The benchmark configuration file
================================

A benchmark configuration has the following mandatory sections:

.. code-block:: yaml

  title : >-
    A title describing the benchmark set-up

  description: >-
    A verbose, high-level description of the benchmark - its objective
    and rationale, the data sets, methods and metrics.

  tags:
    - A tag describing the benchmark
    - Another tag describing the benchmark
    - ...

  setup:
    # list of tools
    tools:
      - tool1
      - tool2

    # list of metrics
    metrics:
      - metric1
      - metric2

  # input data sets
  input:

    reference: ref.fa

    bam:
      - A.bam
      - B.bam

  # set options for tool1
  tool1:
     options: --regions=chr1

The labels in the ``input`` section need to correspond to the slots
defined by the :attr:`expected` attribute in the :class:`.ToolRunner`
that have been listed in the setup section. For example, variant
callers are derived from a class :class:`.VariantCaller` that defines
the common data types for all callers::

    class VariantCaller(ToolRunner):
        expected = ["reference", "bam"]

Based on the benchmark configuration file, a workflow is created that
will instantiate all possible combinations of tools, metrics and
data sets and execute them. The rules for building combinations are as
follows:

1. Each tool is combined with each input data set.

2. If an input data slot contains a list of values, each tool will
   be run on each data set.

3. Each combination/tool will be run against each metric.

In the example above, the following will be executed::

  tool1 X metric1 X ref.fa + A.bam
  tool1 X metric1 X ref.fa + B.bam
  tool1 X metric2 X ref.fa + A.bam
  tool1 X metric2 X ref.fa + B.bam
  tool2 X metric1 X ref.fa + A.bam
  tool2 X metric1 X ref.fa + B.bam
  tool2 X metric2 X ref.fa + A.bam
  tool2 X metric2 X ref.fa + B.bam

It is possible to group input sets. Variant callers typically accept
several .bam files for joint calling. To implement this, group bam
files in an additional level:

.. code-block:: yaml

    bam:
      - pedigree1
        - A.bam
        - B.bam
      - pedigree2
        - C.bam
        - D.bam

This will result in the following combinations::

  tool1 X metric1 X ref.fa + (A.bam + B.bam)
  tool1 X metric1 X ref.fa + (C.bam + D.bam)
  tool1 X metric2 X ref.fa + (A.bam + B.bam)
  tool1 X metric2 X ref.fa + (C.bam + D.bam)
  tool2 X metric1 X ref.fa + (A.bam + B.bam)
  tool2 X metric1 X ref.fa + (C.bam + D.bam)
  tool2 X metric2 X ref.fa + (A.bam + B.bam)
  tool2 X metric2 X ref.fa + (C.bam + D.bam)

For this mechanism to work, the :term:`tool` needs to be aware
that it might receive a single or multiple files. The method
:func:`.resolve_argument` helps here. In the example below, the
tool expects a `,` separated list of input files::

    def run(self, outfile, params):
        bam = resolve_argument(params.bam, sep=",")
        retval = P.run("{params.path} "
                       "--inputs {bam} "
		       "> {outfile} ")

Tool/metric configuration
=========================

Tools and metrics can receive optional (or required) configuration
arguments in their own sections. The configuration options are grouped
into sections within the configuration file named according to the
metric or tool:

.. code-block:: yaml

   tool1:
      options: --region=chr1
      
   metric1:
      reference: ref.fa

This will provide the option ``--verbose`` when running `tool1` and the
parameter ``reference`` to `metric1`. Note that the tool and metric
runner need to be aware of these options. See more about writing
tools and metrics in :ref:`tasklibrary`.

Multiple versions can be specified to provide an additional level of
combinations. For example:

.. code-block:: yaml

   tool1:
      options:
        - --region=chr1
        - --region=chr2

   metric1:
      reference:
        - ref.fa
	- other_ref.fa

will run `tool1` with options ``--region1`` and ``--region2`` and
`metric1` with two different reference data sets. Shared options can
be specified using the ``prefix`` special command.

.. code-block:: yaml

   tool1:
      options:
        - prefix=--verbose
        - --region=chr1
        - --region=chr2

By default, tools and metrics are expected to reside in the user's
:envvar:`PATH` variable. To run a particular version of a tool, use
the `path` configuration value:

.. code-block:: yaml
	
   weCall:
      path: /path/to/weCall/bin/weCall

Note that this can also be multiplexed. To run several versions of
a tool in a benchmark, type:

.. code-block:: yaml

   weCall:
      path:
         - /path/to/weCall-old/bin/weCall
         - /path/to/weCall-new/bin/weCall

Note that this assumes that the executables are entirely
self-contained and automatically pick up references relative to their
location.

Automatic file expansion
========================

To help with the combinatorics, the benchmark file is
aware of glob and find expressions. For example:

.. code-block:: yaml

   input:
      file: find /data/library -name "*.bam"

Will execute the unix ``find`` command and enter all files that
have been found into the daisy.

Filenames containing a `*` are interpreted as glob expressions:

.. code-block:: yaml

    input:
       file: /data/library/1000Genomes/LowCovChr20BAMs/CEU_chr20/NA127*.bam

.. _collation:

Collation
=========

Occasionally, tools need to be run individually, but metrics are
computed on an aggregation of the tool output. For example, you might
want to call variants across a population, but then compute allele
frequencies on the aggregate VCF. For such a workflow, define a
:ref:`collate` task::

  setup:

    tools:
      - weCall
    collate:
      - mergegvcf_agg
    metrics:
      - bcftools_stats

  input:
    reference_fasta: /data/library/reference/hs37d5/hs37d5.fa
    bam: sample*.bam
    regex: (\S+).bam

  mergegvcf_agg:
    regex_in: (\S+).dir/result.vcf.gz
    pattern_out: result.vcf.gz
    runner: illumina_agg

  illumina_agg:
    reference_fasta: /data/library/reference/hs37d5/hs37d5.fa

The workflow above will run weCall on all bam files matching the glob
expression. The output will then be submitted to a collate task called
``mergegvcf_agg``. The task describes how input files should be
grouped (``regex_in`` and ``pattern_out``) and which tool should be
used for merging (``runner``). The tool (``illumina_agg``) is then
configured in a separate section.

.. _splitting:

Splitting
=========

The output of tools may be split in order to compute metrics on parts
of the output separately. For example, the following will split the
output by chromosome and then apply all metrics on both the original
output and all the split files::

  setup:
    ...

    split:
      - split_by_chrom

  split_by_chrom:
    runner: vcf_by_chromosome

.. _exporting:

Exporting
=========

The benchmark system can export tool data for further use. To export,
simply type::

    benchmark run make export

This will move all output files into a directory called
:file:`export.dir` and place symbolic links into the pipeline
directives to preserve workflow state.

Files in the directory :file:`export.dir` will be renamed to label
them according to the experiment. For example,
:file:`weCall_NA12878.dir/result.vcf.gz` will become
:file:`export.dir/weCall_NA12878.vcf.gz`.

The :term:`export` target is a convenience function to collect all the
tool data computed in an experiment if the tool data is of further
interest, for example for additional processing in other benchmarks.

Global configuration
====================

Below is a configuration values for interfacing Benchmark with the
system.

Cluster
-------

Options to configure behaviour for running jobs on the cluster are
in the section ``cluster``. The default values are::

   cluster:
      queue: main.q
      priority: -1
      num_jobs: 100
      memory_resource: h_vmem
      memory_default: 4G
      options: ""
      parallel_environment: smp

Note that some cluster options can be overridden at the command
line. For example, ``--cluster-queue=slow.q`` will send jobs to
``slow.q``.  The options ``--local`` will run jobs without the
queuing system.

Database
--------

Database access is implemented through setting a database URL.
The default is:

.. code-block:: yaml

   database:
      url: postgres://localhost:5432/Benchmark

With postgresql_ it is possible to use schema to organise metric
tables. To use a schema, use:

.. code-block:: yaml
		
   database:
     url: postgresql://andreas@trafalgar.camdc.genomicsplc.com/benchmark
     schema: cnv

