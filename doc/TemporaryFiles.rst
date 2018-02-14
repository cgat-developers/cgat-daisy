============================
Working with temporary files
============================

Using $TMPDIR
=============

Daisy will set the TMPDIR environment variable in its job scripts
to a location that is unique for each job. The location will be
created before a job and removed after a job exists. For example::

    def run(self, outfile, params):
	statement = (
	    "dosth > $TMPDIR/myfile; "
	    "wc -l $TMPDIR/myfile > {outfile}"
	    .format(**locals()))
	return P.run(statement)

Note that :file:``myfile`` will be cleaned up automatically after
the job completes.

Using daisy API
===============

The daisy toolkit provides several API functions to work with
temporary files. The functions :func:`getTempFilename`,
:func:`getTempFile` and :func:`getTempDir` provide these. These
functions are aware of the temporary storage locations either
specified in configuration files or on the command line and
distinguish between the ``private`` locations that are visible only
within a particular compute node, and ``shared`` locations that are
visible between compute nodes and typically on a network mounted
location.

Note that all these functions will not clean up the temporary files
and directories have been created. For example::

    def run(self, outfile, params):
        tmpdir = P.getTempFilename(clear=True)
	statement = (
	    "mkdir {tmpdir}; "
	    "dosth > {tmpdir}/myfile; "
	    "rm -rf {tmpdir}").format(**locals()))
	return P.run(statement)

Note how the user code both creates and removes the temporary
directory. If there is an error in ``dosth``, temporary files will
remain. This can be desired behaviour, for example in long-running
jobs where part of the computation can be salvaged if aborted.
