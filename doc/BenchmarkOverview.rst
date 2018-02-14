========
Overview
========

The purpose of the Benchmark library is to provide a platform
for running and managing benchmarks. Typical scenarios are:

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

3. Tasks are executed on the cluster for improved performance. The
   workflows are controlled by a command-line interface.

4. Results of benchmarks are uploaded to a database to allow data
   mining and tracking of benchmark results over time.
   


