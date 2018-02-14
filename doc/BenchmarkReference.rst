.. _reference:

=======================
Reference
=======================

The :ref:`reference` lists all the :term:`sections` and
:term:`directives` that may appear in a benchmark configuration file.

.. _title:

title
=====

String (required)

The :ref:`title` section contains the title. For example::

  title : >-
    Evaluate variant caller performance

.. _description:

description
===========

String (required)

The :ref:`description` section contains a description of the
experiment. For example::

   description: >-
      This experiment calls variants on the Platinum data set using a
      variety of variant callers. Performance is evaluated through
      comparison to the GiaB Standard of Truth.

.. _tags:

tags
====

List of strings (required)

The :ref:`tags` section contains a list of strings that can be used
to categorize the experiment. For example::

   tags:
     - Variant calling
     - Genome In a Bottle
     - Platinum data

.. _setup:

setup
=====

Collection of directives (required)

.. glossary::

   tools
      (required)
      A list of tools that should be executed on the input data, for
      example::

	tools:
	   - weCall
	   - platypus

      Each tool can be parameterized in a separate :ref:`tool_task`
      section.

   metrics 
      (required)
      A list of metrics that should be executed on the output of
      the tools, for example::

	metrics:
	   - bcftools_stats
           - rtg_vcfeval

      Each tool can be parameterized in a separate :ref:`metric_task`
      section.

   collate
      (optional)
      A list of collation tasks that should be executed on the tool
      output. See :ref:`collation` for an example and :ref:`collate-task`.

   split
      (optional)
      A list of split tasks that should be executed on the tool output.
      See :ref:`splitting` for an example.

   ignore 
      (optional)
      A list of input/tool/metric names that should be ignored.

   only_collate
      (optional)

      Flag. Do not compute any metrics on the output of tool tasks, but only
      on the output after collation.
      
   collate_metrics
      (optional)

      A list of metrics that will only be applied to the output of collation tasks.

   split_metrics
      (optional)
      
      A list of metrics that will only be applied to the output of split tasks.

   external
      (optional)
   
      Add external data to the daisy. Metrics are computed on the
      external data alongside the output by the tools run by the
      benchmark system. This section must contain an :ref:`input`
      section and an `output` directive::

        setup:
	  tools:
	    - bwa_mem
	  metrics:
	    - samtools_stats
	  external:
	    input:
              file: "/other_data/2016_sample*.bam"
              regex: 2016_(sample.*).bam
	    output: result.bam
	    add_glob: .bai

       The functionaly of the external section is similar to the
       :ref:`run_tool_identity` and in fact uses the same
       implementation and accepts the same options.  The difference is
       that the :term:`external` directive permits computing metrics
       on data created by the pipeline alongside external data, while
       the :ref:`run_tool_identity` assumes that only metrics will be
       applied.

   export
       (optional)

       A list of tasks which output should be exported. By default,
       the output of the :term:`tools`, :term:`collate` and the :term:`split`
       tasks will be exported. To only export the output of the :term:`tools`
       section, use::

	  setup:
	    export:
	      - tools
       
       Additionally, there are directives for determining the name of
       output files when exporting tool output data. Currently defined
       is:

       .. glossary::
          prefix
	     (optional) add a prefix to exported output files.

       For more information, see :ref:`exporting`.

.. _input:

input
=====

Collection of directives (required)

The input section contains the filenames for the input data. Files are
labeled according to slots defined by a tool, for example::

   input:
      reference_fasta: hg19.fa
      bam:
         - individual1.bam
	 - individual2.bam

Filenames can contain keywords such as :term:`glob` or :term:`find` to
refer to a group of files. Each file encountered by a :term:`glob` or
:term:`find` expression will be added as an item to the list of files
to be processed and thus each file will instantiate a separate task.
Multiple glob statements can be separated by a comma (``,``)::

   input:
      reference_fasta: hg19.fa
      bam: individual*.bam, sample*.bam

In addition, the :ref:`input` section may contain additional
directives.

.. glossary::

   regex 
      (optional)

      A string containg a regular expression to extract a name from a
      filename. The regular expression should contain at least one
      ``()``-group. For example::

        input:
           reference_fasta: hg19.fa
           bam:
             - family1_individual1.bam
	     - family2_individual2.bam
	   regex: (\S+)_(\S+).bam
	   alias: \2
	
      will set the aliases ``indivial1`` and ``individual2`` ignoring
      the family.

   alias 
     (optional)

      A pattern that can be used to build a name from the regular
      expression given by :term:`regex`. The default is to concatenate
      all groups in the regular expression separated by an underscore
      (``_``). See :term:`regex`.

   groupby 
      (optional)

      Either ``option`` or ``label``. This option determines how input
      files should be grouped. The default is ``option``, so that
      files will be grouped across labels. For example::
      
         input:
           reference_fasta: hg19.fa
           bam:
	     pair1:
               - individual1.bam
	       - individual2.bam
	     pair2:
               - individual1.bam
	       - individual2.bam
      
      will result in the following pairs::

         {"reference_fasta": "hg19.fa", "pair1": {"bam": ("individual1.bam", "individual2.bam")}}
         {"reference_fasta": "hg19.fa", "pair2": {"bam": ("individual3.bam", "individual4.bam")}}

      while ::

         input:
	   pair1:
             reference_fasta: hg19.fa
             bam:
               - individual1.bam
	       - individual2.bam
	   pair2:
             reference_fasta: hg38.fa
             bam:
               - individual3.bam
	       - individual4.bam
	   groupby: label

      will result in::

         {"pair1": {"reference_fasta": "hg19.fa", "bam": ("individual1.bam", "individual2.bam")}}
         {"pair2": {"reference_fasta": "hg38.fa", "bam": ("individual3.bam", "individual4.bam")}}

   group_regex 
     (optional)

     A regular expression used to group input files. For example, if you are interested in calling
     variants inside families and the files are named ``family-sample.bam``, use::

       input:
          bam: *.bam
          group_regex: (\S+)-(\S+).bam
	  group_alias: \1

   group_alias 
     (optional)

     String used to build an alias for a group. See :term:`group_regex`.

   ignore
      (optional)

      Ignore a particular tool or metric. This directives accepts a list
      of patterns::

        ignore:
   	  - gatk_haplotype_caller_WES_NA12891_remapped_dedup

.. _tool_task:

tool-task
=========

A tool task paramaterizes a tool further. For example::

   setup:
     tools:
       - weCall

   weCall:
     options: =regions=1

will run weCall only on chromosome 1.

The benchmark system allows the user the specify multiple alternative
configurations of a tool. Thus, if given a list of configurations, all
of these will be run alongside each other. For example, the following
will run the tool ``weCall`` twice, once on chromosome 1 and once on
chromosome 2::

   weCall:
     options:
        - =regions=1
        - =regions=2

The system requires unique names for each task. By default, these will
be created through hashing the options. To define names explicitely to
facilitate further analysis, use the :term:`alias=` directive. Instead
of setting them explicitely, aliases can be derived automatically from
option names using the :term:`regex` and :term:`alias` directives.

Option strings can be created programmatically with the :term:`generate=`
directive. The full list of directives is below:

.. glossary::

   prefix=

      (optional)
      Shared list of values for a particular option.

      Default options that are common to all tasks can be specified with the 
      :term:`prefix=` directive::

	 weCall:
	   options:
	      - prefix==jobThreads=10
	      - =regions=1
	      - =regions=2
   alias=
      (optional)
      Set an explicit alias for an option::

	weCall:
	  options:
	     - prefix==jobThreads=10
	     - alias=chr1; =regions=1
	     - alias=chr2; =regions=2

   regex
      (optional)
      Regular expression to derive a name using the options submitted to
      the task::

	weCall:
	  options:
	     - prefix==jobThreads=10
	     - =regions=1
	     - =regions=2
	  regex: =regions=(\S+)
	  alias: chr\1

   alias
      (optional)
      String used to build a name from the parts extracted by a regular
      expression (see :term:`regex`).

   generate=
      (optional)
      Generator expression to create a list of options automatically::

        weCall:
          options:
            - prefix==jobThreads=10
            - generate=["alias=chr{}; {}".format(x, x) for x in [1, 2]]

   ignore
      (optional)
      Ignore a particular tool or metric. This directives accepts a list
      of patterns. Any task matching that contains any of the strings in the
      list will be ignored::

        ignore:
   	  - gatk_haplotype_caller_WES_NA12891_remapped_dedup

   task_specific:
      (optional)
      Apply task specific options to a particular command. This directive
      accepts a collection of patterns and appropriate parameters. For example,
      to apply additional filters to metrics compute on freebayes output, use::

        task_specific:
	  freebayes.*: 
             filter_exclude: "FORMAT/GT == '.' || DP < 5 || QUAL < 20 || N_ALT >= 2"
	
.. _metric_task:

metric-task
===========

A tool task paramaterizes a metric further. This section is
identical to :ref:`tool_task`.


.. _collate_task:

collate-task
============

A collate task describes how output data should be grouped. See
:ref:`collation` for an example.

.. glossary::

   regex_in
      (required) regular expression that determines how files should be grouped.

   pattern_out
      (required) output pattern. If all files should be merged, this will simply
      be the filename used by the preceeding tools, for example::

      regex_in: (\S+).dir/result.vcf.gz
      pattern_out: result.vcf.gz

   runner
      (required) name of the tool to be run for combining multiple files
      into one.

.. _split:

split-task
==========

A split task describes how output data should be split before
computing metrics. See :ref:`splitting` for an example.

.. glossary::

   runner
      (required) name of the tool to be run for combining multiple files
      into one.
   
.. _database:

database
========

Collection of directives (optional)

This section contains directives with database connection details.

.. glossary::

   url (optional)
      Database URL. See `here <http://docs.sqlalchemy.org/en/latest/core/engines.html#database-urls>`_ for
      a list of accepted formats. The system is currently tested with sqlite and postgres.

   schema (optional)
      Database schema to use for data tables. If not given or the database
      does not support schemas, the data tables will sit alongside the meta
      tables in the database.

.. _cluster:

cluster
=======

Collection of directives (optional)

A collection of options to specify cluster parameters. Typically,
parameters are set with either defaults hardcoded or in a
user-specific configuration file. If there are experiment specific
options, they can also be specified in the :file:`benchmark.yml` file.

.. glossary::

   queue
     (optional) The cluster queue.
     
   priority
     (optional) The job priority. This should be a negative number.

   num_jobs
     (optional) Number of jobs to submit in parallel to the queueing system.

   memory_resource
      (optional) Name of the memory resource

   memory_default
      (optional) Default amount of memory to allocate

   parallel_environment
      (optional) Name of the parallel environment to use for multi-threaded
      applications.

   options
      (optional) Generic options to use for job submissions.
