.. _tasklibrary:

================
Task Library
================

Tasks are python functors and are defined in the TaskLibrary
sub-package within the Benchmark distribution. Functors follow
the naming convention

``run_tool_xyz``
    for a task running a tool called xyz

``run_metric_xyz``
   for a metric running a tool called xyz

Tools and metrics are auto-discovered by the package, hence the
requirement of a naming convention.

Tools and metrics are organised in a class hierarchy that is rooted at
the class :class:`.Runner`. Immediately derived from these are
:class:`.ToolRunner` and :class:`.MetricRunner`. User defined tasks
should be derived from the latter these.

The base classes take care of interfacing with the benchmark system
such as collecting dependencies and runtime statistics.  A tool- or
metric-runner should implement a run() method. The principal
difference between a tool and metric is the call interface for the
run() method. 

A :term:`tool` expects two arguments outfile and params, which
are the output file name and a parameter dictionary. The input file names
are part of the parameter dictionary and can be referred to by
name. The expected names are listed in a class attribute, for
example::

    class run_tool_listdir(outfile, params):
        expected = ["directory"]
	path = "ls"

	def run(self, outfile, params):
	    return P.run("{params.path} -1 {params.directory} "
	    "> {outfile}".format(**locals()))

A :term:`metric` expects three arguments, outfile and params as
before, but also an input file name, for example::

   class run_metric_count(infile, outfile, params):
       path = "wc"
     
       def run(self, infile, outfile, params):
           return P.run("{params.path} -l {infile} "
	   "> {outfile}".format(**locals()))

The distinction exists because the input for tools are files outside
the workflow and static, while the input for a metric are files that are created
within the workflow and are thus dynamically generated.

Both :term:`metric` and :term:`tool` should return the value of the
:meth:`.run` call, which contains runtime information such as the time
to execute the command and the host name.

Parameterisation of tasks
================================

A task can be parameterised by defining member variables
in the task runner::

   class run_metric_count(infile, outfile, params):
       path = "wc"
       counting_method = "lines"

       def run(self, infile, outfile, params):
       
           if params.counting_method == "lines":
	       opt = "-l"
	   elif params.counting_method == "chars":
	       opt = "-c"
	   else:
	       raise ValueError("unknown counting method")
	       
           return P.run("{params.path} {opt} {infile} "
	   "> {outfile}".format(**locals()))

Note how the method refers to the options through the ``params``
variable. This is necessary as the class attributes are shared between
all tasks, while the params contain the options specific to a
task. Using ``params`` permits setting task specific parameters for
a tool in the configuration file.

A parameter called ``options`` is always defined and can be used
directly. This is for generic command-line options. The example above
could thus simply be::

   class run_metric_count(infile, outfile, params):
       path = "wc"
       counting_method = "lines"

       def run(self, infile, outfile, params):
       
           return P.run("{params.path} {params.options} {infile} "
	   "> {outfile}".format(**locals()))

It is generally preferred to keep the number of explicit options to a
minimum and rather let the user set command line options for
the specific tool. 

A use case where it additional options are useful is when the command
line statement will take different shape depending on options. For
example, a metric might offer pre-filters for computing the metric. In
this case it is important to remember that tasks are functors and thus
should not contain state information that is specific to a task.

Additional topics to be written:

* specifying multi-threaded execution
* pre-processing
* database interface

Data upload
===========

The function of a metric task is to run a metric on some data and
output the data in well-formed, tab-separated tables. The benchmark
system then takes care of uploading the data into the database. For
general notes about the database and how data are organized, see
:ref:`data_organization`.

Standard case
-------------

The simplest metric outputs a single table into ``outfile``. In the
example above::

   class run_metric_count(infile, outfile, params):
       path = "wc"
       counting_method = "lines"

       def run(self, infile, outfile, params):
       
           return P.run("{params.path} {params.options} {infile} "
	   "> {outfile}".format(**locals()))

``outfile`` will be created in the approriate location and could look
like this::

   >cat tool.dir/count.dir/count.tsv
   word   count
   hello  1
   world  1


This table on file will automatically be uploaded into a database
table called ``count``::

   > SELECT * FROM count;
   word   count    instance_id
   hello  1        5
   world  1        5

Note how the column :term:`instance_id` has been added to link the
results to a specific :term:`instance`.

Multiple outfiles
-----------------

Frequently, a metric might output multiple tables. To register these
into the system, define them as separate tablenames::

   class run_metric_count(infile, outfile, params):
       path = "wc"
       counting_method = "lines"

       tablenames = ["count_lines", "count_words"]

       def run(self, infile, outfile, params):
       
           retvals = P.run("{params.path} --count-words {params.options} {infile} "
	   "> {outfile}.count_words.tsv".format(**locals()))

           retvals.extend(P.run("{params.path} --count-lines {params.options} {infile} "
	   "> {outfile}.count_lines.tsv".format(**locals())))

	   return retvals
	   
Note how the names of the table (e.g. ``count_lines``) imply the name
of the ouput file
(:file:`tool.dir/count.dir/count.tsv.count_lines.tsv`).

.. note::

   If the naming convention is not followed, existing tables will not
   be picked up.

A missing output file is ignored. This accomodates metrics that create
multiple output files of which some are optional.

Transforming data before upload
-------------------------------

Before uploading, the Benchmark system can apply some basic table transformations.
These are registered in the task function. 

.. glossary::

   upload_melted

      melt the table before uploading. This is important for metrics
      that output a variable number of columns, for example, if a tool
      outputs one column per sample the database. The argument is a
      dictionary mapping the table name to parameters in the 
      `pandas melt <http://pandas.pydata.org/pandas-docs/stable/generated/pandas.melt.html>`_
      function. For example::
      
         upload_melted = {
            "benchmark_vcf_stats_format_unset_samples":
            {"id_vars": ["FORMAT"],
             "var_name": "sample",
             "value_name": "unset"},
            }
	    
   upload_transpose
       
      transpose a table before uploading. The argument is a list of tables
      to transpose::

         upload_transpose = ["count_lines"]

   upload_separate

       upload each data point into a separate table. By default and design, 
       data from the same metric are stored in the same table. Using this option,
       each run will create a separate table which is the name of the metric
       suffixed with the :term:`instance_id`. Use sparingly, as a large number
       of tables will clutter the database. The argument to this option is 
       a list of tables to uload separately::

          upload_separate = ["count_lines"]
      
   upload_normalize

       normalize a column in a table. This is useful if the table contains
       one or more columns with categorical data. To save disk space, all the
       categorical data will be replaced by integer levels and the mapping
       between level and values will be stored in a separate table. The argument
       to this option is a dictionary of tables and a list of columns that
       should be stored as factors::
       
         upload_normalize = {"count_words": ["word"]}
	 
       In our example, this will create an additional table
       ``count_factors``.  Note that factor values are only consistent
       within an instance, but not across instances or
       experiments. Multiple columns will be normalized together as a
       combination of values.

   upload_as_matrix

       upload a table as a matrix. The argument is a dictionary of
       tables. Optionally, a group key can be specified::

          upload_as_matrix = {
            "benchmark_vcf_stats_format_per_sample": "FORMAT",
          }

       This transformation assumes that the resulting data is uniform,
       i.e. that all columns have the same type. Uploading as a matrix
       is advisable for data that have data dependent labels in both
       rows and columns. Melting such data will result in a massively
       inflated data size that needs to be stored. The resulting table
       will contain the fields ``instance_id``, ``rows``, ``columns``,
       ``data`` and ``dtype`` (data type), where rows and columns are
       ``,`` separated lists of values.

Contents of the task library
================================
The inheritance diagram of the TaskLibrary is below:

.. inheritance-diagram::
   daisy.tasks.Runner
   daisy.tasks.ToolRunner
   daisy.tasks.MetricRunner
   daisy.tasks.BAMTools
   daisy.tasks.BAMMetrics
   daisy.tasks.VariantCallers
   daisy.tasks.VCFTools
   daisy.tasks.VCFMetrics

Task library methods
------------------------------------------------

This section lists modules containing tool runners and metric runners.

.. automodule:: daisy.tasks.BAMTools
   :undoc-members:
   :members:

.. automodule:: daisy.tasks.BAMMetrics
   :undoc-members:
   :members:

.. automodule:: daisy.tasks.VCFMetrics
   :undoc-members:
   :members:

.. automodule:: daisy.tasks.VariantCallers
   :undoc-members:
   :members:

The task library engine
===========================

This section lists modules that are part of the task library engine

.. automodule:: daisy.tasks.Runner
   :undoc-members:
   :members:

.. automodule:: daisy.tasks.MetricRunner
   :undoc-members:
   :members:

.. automodule:: daisy.tasks.ToolRunner
   :undoc-members:
   :members:
