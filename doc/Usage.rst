============
Usage
============

The benchmarking system is running tools and metrics and multiple input
data and organising the data into a database. Each benchmark is
defined by an experimental design that is described in a yaml
configuration file.

1. Running a collection of tools on the same data set and comparing
   the output.

2. Running one tool on different data sets to gauge its performance.

3. Running one tool with different options on the same data set to
   optimise its performance.

The benchmarking platform provides a system to run these benchmarking
workflows in a systematic fashion. Its hallmarks are:

1. Benchmarks are set-up by human readable configuration files.

2. Tools and metrics are organised in a :ref:`Task Library` to provide
   standardised and reusable components.

3. Tasks can be executed on a cluster to run computations at
   scale. The workflows are controlled by a command-line interface.

4. Results of benchmarks are uploaded to a database to allow data
   mining and tracking of benchmark results over time.


.. toctree::
   :maxdepth: 2
  
   BenchmarkUsage.rst
   BenchmarkConfiguration.rst
   BenchmarkDatabase.rst
   BenchmarkReference.rst




