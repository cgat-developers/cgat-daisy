==================
Using the platform
==================

This section provides a brief tutorial on how to run the Benchmark
suite. 

Creating a benchmark
====================

To set up a benchmark, create an empty directory and create a
configuration file called :file:`benchmark.yml`, for example:

.. code-block:: yaml

   setup:
     suffix: vcf.gz

     tools:
       - weCall
       - samtools
       - platypus
     metrics:
       - vcftools_tstv_summary
       - vcftools_tstv_by_count

   input:
     bam: /data/library/PlatinumGenomes/ERP001960/NA12878_S1.bam
     reference: /data/library/reference/37hg19_chr/37hg19_chr.fa

   weCall:
     options: --regions chr20:0-1000000

   samtools:
     options: -r chr20:0-1000000

   platypus:
     options: --regions chr20:0-1000000

The `setup` section lists the tools (`weCall`, `samtools`, `platypus`)
we want to run on our input data and the metrics
(`vcftools_tstv_summary`, `vcftools_tstv_by_count`) we want to collect
from the output. The field `suffix` sets the file suffix to be used
for the tool output.

The input section lists the input data sets we want the tools to apply
to. Here, there are two input files, a :term:`bam` formatted file with
short read data and a genomic reference sequence in :term:`fasta`
format.

The additional sections provide parameterizations for tools or metrics. In
the example above, we limit the tools to a particular genomic region.

Running a benchmark
===================

To run a benchmark, type::

    python <SRCDIR>/bin/daisy run -v 5 -p 10 make

SRCDIR is the location of the repository. Type::

    python <SRCDIR>/bin/daisy run --help

to get a list of command line options to control execution behaviour.

:file:`benchmark.py` implements a simple 3 step workflow. The script
will run the callers in parallel on the cluster (if available) and
once finished, apply the metric tools. The data is then uploaded to a
database. The directory will look like this::

    daisy.yml
    pipeline.log
    weCall_e8a713_57e5b3.dir
    weCall_e8a713_57e5b3.dir/tool.info
    weCall_e8a713_57e5b3.dir/result.vcf.log2
    weCall_e8a713_57e5b3.dir/result.vcf.log
    weCall_e8a713_57e5b3.dir/result.vcf.gz
    weCall_e8a713_57e5b3.dir/result.vcf.gz.tbi
    weCall_e8a713_57e5b3.dir/tool.bench
    weCall_e8a713_57e5b3.dir/vcftools_tstv_by_count_d75171.dir
    weCall_e8a713_57e5b3.dir/vcftools_tstv_by_count_d75171.dir/vcftools_tstv_by_count.info
    weCall_e8a713_57e5b3.dir/vcftools_tstv_by_count_d75171.dir/vcftools_tstv_by_count.tsv
    weCall_e8a713_57e5b3.dir/vcftools_tstv_by_count_d75171.dir/vcftools_tstv_by_count.bench
    weCall_e8a713_57e5b3.dir/vcftools_tstv_summary_d75171.dir
    weCall_e8a713_57e5b3.dir/vcftools_tstv_summary_d75171.dir/vcftools_tstv_summary.info
    weCall_e8a713_57e5b3.dir/vcftools_tstv_summary_d75171.dir/vcftools_tstv_summary.tsv
    weCall_e8a713_57e5b3.dir/vcftools_tstv_summary_d75171.dir/vcftools_tstv_summary.bench
    samtools_a17daa_57e5b3.dir
    samtools_a17daa_57e5b3.dir/tool.info
    samtools_a17daa_57e5b3.dir/result.vcf.gz.log
    samtools_a17daa_57e5b3.dir/result.vcf.gz
    samtools_a17daa_57e5b3.dir/result.vcf.gz.tbi
    samtools_a17daa_57e5b3.dir/tool.bench
    samtools_a17daa_57e5b3.dir/vcftools_tstv_by_count_d75171.dir
    samtools_a17daa_57e5b3.dir/vcftools_tstv_by_count_d75171.dir/vcftools_tstv_by_count.info
    samtools_a17daa_57e5b3.dir/vcftools_tstv_by_count_d75171.dir/vcftools_tstv_by_count.tsv
    samtools_a17daa_57e5b3.dir/vcftools_tstv_by_count_d75171.dir/vcftools_tstv_by_count.bench
    samtools_a17daa_57e5b3.dir/vcftools_tstv_summary_d75171.dir
    samtools_a17daa_57e5b3.dir/vcftools_tstv_summary_d75171.dir/vcftools_tstv_summary.info
    samtools_a17daa_57e5b3.dir/vcftools_tstv_summary_d75171.dir/vcftools_tstv_summary.tsv
    samtools_a17daa_57e5b3.dir/vcftools_tstv_summary_d75171.dir/vcftools_tstv_summary.bench
    platypus_e8a713_57e5b3.dir
    platypus_e8a713_57e5b3.dir/tool.info
    platypus_e8a713_57e5b3.dir/result.vcf.log
    platypus_e8a713_57e5b3.dir/result.vcf.gz
    platypus_e8a713_57e5b3.dir/result.vcf.gz.tbi
    platypus_e8a713_57e5b3.dir/tool.bench
    platypus_e8a713_57e5b3.dir/vcftools_tstv_by_count_d75171.dir
    platypus_e8a713_57e5b3.dir/vcftools_tstv_by_count_d75171.dir/vcftools_tstv_by_count.info
    platypus_e8a713_57e5b3.dir/vcftools_tstv_by_count_d75171.dir/vcftools_tstv_by_count.tsv
    platypus_e8a713_57e5b3.dir/vcftools_tstv_by_count_d75171.dir/vcftools_tstv_by_count.bench
    platypus_e8a713_57e5b3.dir/vcftools_tstv_summary_d75171.dir
    platypus_e8a713_57e5b3.dir/vcftools_tstv_summary_d75171.dir/vcftools_tstv_summary.info
    platypus_e8a713_57e5b3.dir/vcftools_tstv_summary_d75171.dir/vcftools_tstv_summary.tsv
    platypus_e8a713_57e5b3.dir/vcftools_tstv_summary_d75171.dir/vcftools_tstv_summary.bench
    shell.log
    results.commit

There are three output directories starting with the tool-name. These
contain a file :file:`result.vcf.gz` with the tool output. Each of
these directories in turn contains further subdirectories with the
output of the metrics that have been run onto the tool outputs.

The log file :file:`pipeline.log` contains the commands executed and
will report any errors encountered. If all has been succesful, it
should contain a line such as::

    ## 2015-11-04 17:01:09,436 INFO uploaded results under run_id 112

meaning that our data have been uploaded to the database and are
accessible under run_id 112.

Adding a new tool
=================

Let us add the variant caller freebayes_ to our daisy. We have
installed it and it is on our path. The command ``freebayes -h``
provides us with its commandline options. Its basic usage is::

    freebayes -f input.fa [OPTIONS] input.bam >output.vcf                                                                              

To enable freebayes_, add the following code to a python module in
Benchmark.TaskLibrary:

.. code-block:: python

   from .Runner import resolve_argument
   from .ToolRunner import ToolRunner
   import Benchmark.Experiment as E
   import Benchmark.Pipeline as P

   class run_tool_freebayes(ToolRunner):
       name = "freebayes"
       path = "freebayes"
       expected = ["bam",  "reference"]

       @property
       def version(self):
	   help_string = E.run("{self.path} --version".format(**locals()),
			       return_stdout=True).strip()
	   return re.search("version:\s+(\S+)", help_string).groups()[0]

       def run(self, outfile, params):

	   bam = resolve_argument(params.bam, sep=" ")

	   return P.run("{params.path} "
			"--fasta-reference {params.reference} "
			"{params.options} "
			"{bam} "
			"| bgzip "
			"> {outfile}; "
			"tabix -p vcf {outfile}"
			.format(**locals()))

The first lines import functions and classes within the benchmark suite.

.. code-block:: python

   from .Runner import resolve_argument
   from .ToolRunner import ToolRunner
   import Benchmark.Experiment as E
   import Benchmark.Pipeline as P

The next section defines our task. 

.. code-block:: python

   class run_tool_freebayes(ToolRunner):
       name = "freebayes"
       expected = ["bam",  "reference"]
       path = "freebayes"

The task`s name ``run_tool_freebayes`` makes sure that our task is
automatically identified as a tool within our Task Library. The
attribute :attr:`name` links this task with a name in the
configuration file. 

The parameter :attr:`expected` lists the input data that our tool
expects. The section ``input`` in the :file:`benchmark.yml` file
needs to provide these. Finally, :attr:`path` identifies the name 
of the executable.

The next section implements a command line call to obtain the
version of the tool. Every task function should provide this.

.. code-block:: python

       @property
       def version(self):
	   help_string = E.run("{self.path} --version".format(**locals()),
			       return_stdout=True).strip()
	   return re.search("version:\s+(\S+)", help_string).groups()[0]

Finally, the tool will be exectuted in the :meth:`run()` method:

.. code-block:: python

       def run(self, outfile, params):

	   bam = resolve_argument(params.bam, sep=" ")

	   return P.run("{params.path} "
			"--fasta-reference {params.reference} "
			"{params.options} "
			"{bam} "
			"| bgzip "
			"> {outfile}; "
			"tabix -p vcf {outfile}"
			.format(**locals()))

Basically, a command line statement is built from arguments
supplied to the task, the output file and a class representing
the options supplied to this method. Note the reference to
``params.reference`` and ``params.path`` to access these.
The command line statement is sent to the :meth:`Pipeline.run()`
method to execute it either on the cluster or locally, depending on
the input sections.

Just by adding this section of the code to our daisy.yml file
we now have the tool freebayes_ available and we can add it to
our comparison:

.. code-block:: yaml

   ...
   setup:
     suffix: vcf.gz

     tools:
       - weCall
       - samtools
       - platypus
       - freebayes
    ...

This is all that is required. Note that when re-running the pipeline::

    python <SRCDIR>/bin/daisy run -v 5 -p 10 make

only freebayes_ will be executed as the system detects that the
previous files are up-to-date and need not be recomputed. Note that
we can also supply options to freebayes:

.. code-block:: yaml

   ...

   freebayes:
       options: --report-genotype-likelihood-max
   ...

Closing remarks
===============

The benchmark system has been designed to be easy to use while
at the same time providing maximum flexibility. Thus, quite a few
things are happening behind the scenes. In particular, look out for
the following features:

1. Collating output files for easier analysis.
2. Task specific parameterization

.. _freebayes: https://github.com/ekg/freebayes
