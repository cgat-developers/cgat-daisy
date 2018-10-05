=======
Testing
=======

The benchmark system is covered by acceptance tests. These can be
found in the directory :file:`test`. To execute all acceptance tests,
type::

   nosetests test

Acceptance tests
================

To run the acceptance tests, type::

   nosetests test/acceptance

Regression tests
================

The testing setup will check all workflows stored in the benchmark
library (:file:`daisy-chains`) for compatibility with the current
benchmark system. It does this by checking if the "show all" command
will return without failure. 

In most cases, a configuration file should be able to pass the
test-suite without issues. However, occasionally glob and regex
expressions might cause errors as they assume files or directories
that are not present during the tests. To circumvent these issues, add
a :term:`test` directive to the configuration file specifying explicit
filenames::

  fastq: ../data/HG*.fastq.gz

  test:
    fastq:
      ../data/HG001-R1.fastq.gz
      ../data/HG001-R2.fastq.gz
      ../data/HG002-R1.fastq.gz
      ../data/HG002-R2.fastq.gz

  group_regex: /([^/]+)-R.*.fastq.gz

The :term:`test` directive is evaluated when ``--is-test`` is specified
at the command line and replaces any fields of the same name in the
:term:`input` section. Here, the glob expression is replaced by four
dummy filenames that ensure that the :term:`group_regex` directive
will match. 

The testing system will also apply regression tests to make sure that
the workflows have not changed. It does this by comparing the current
output to previous workflow output. Previous workflow outputs are
stored in the directory :file:`test/regression/daisy-chains`,
which mirrors in its layout the directory :file:`daisy-chains`.

If a workflow changes, simply delete the existing file in
:file:`test/regression/daisy-chains` and rerun the tests, a new
copy reflecting the latest workflow with be created automatically.
The directory `test/regression/daisy-chains` is under version
control and any changes need to commited. To re-create all regression
test data, type::

   rm -rf test/regression/daisy-chains
   nosetests test/regression
