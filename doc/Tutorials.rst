=========
Tutorials
=========

Creating a task
===============

This tutorial shows how to add a new task into the :ref:`TaskLibrary`
so that it can be executed as part of a workflow vi `:doc:`tools/run`
or as a stand-alone task via :doc:`tools/run-task`

To put it quite simply, a task is a wrapper class around a command
line statement. There are two types of tasks within the task library
that have slightly different signatures. The simplest is a
:term:`metric` that takes one input file and outputs one output file
in :term:`tsv` format. To add a task, simply create a new file called
:file:`MyTask.py` inside the :file:`TaskLibrary` directory within the
repository::

    import daisy.Pipeline as P
    from .MetricRunner import MetricRunner

    class run_metric_filestat(MetricRunner):
	name = "filestat"
	path = "stat"

	def run(self, infile, outfile, params):
	    return P.run("{params.path} --printf='filename\\tsize\\n%%n\\t%%s\\n' "
			 "{infile} > {outfile}".format(**locals()))

The metric will be auto-discovered.




