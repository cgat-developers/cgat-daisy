===============
Installation
===============

Code installation
=================

To install the benchmarking suite from the repository, clone the
latest version from github::

    git clone https://github.com/cgat-developers/cgat-daisy
    cd cgat-daisy
    python setup.py install

Alternatively, daisy release packages are available via pypi::

   pip install cgat-daisy

and conda::

   conda install -c bioconda cgat-daisy


Configuration
=============

Daisy connects to the cluster. To set global configuration values, it
looks for the file :file:`daisy.yml` in the user's home directory. If
it exists it will be used to set default values. A typical
configuration file might look like this::

  cluster:
    queue: all.q
    options: ""
    priority: -10
    memory_resource: h_vmem
    memory_default: 4G
    parallel_environment: shmem

For more information, see :term:`Parameterisation`.
    
External dependencies
======================

System packages
--------------------

Daisy requires certain system packages installed. Below is
an (incomplete) list:

* libsqlite3-dev
* libfreetype6
* libfreetype6-dev
* libpq-dev
* gcc
* libjpeg-dev
* libx11-dev

.. clustersetup:

Using the Cluster
------------------------------------

In order to use the cluster, libdrmaa.so needs to be installed and in
:envvar:`LD_LIBRARY_PATH`.  For SGE, the environment variable
:envvar:`SGE_ROOT` needs to set, for example::

   export SGE_ROOT=/var/lib/gridengine

Tools and Metrics
------------------------------------------

The :ref:`tasklibrary` within the repository contains :term:`Runners`
to execute tools and metrics, but these are only wrappers. The actual
tools and metrics need to be installed separately.

By default, tools and metrics are assumed to be in your :envvar:`PATH`
and accessible from both the submit and execution hosts. The system
also assumes that files on the submit and execution hosts are
accessible from the same paths.

Database setup
==============

Daisy can work with several database engines. For sqlite_, no
additional setup is required. For postgresql_, Daisy needs to
connect to an existing database. The database needs to be set up
prior to running a daisy and appropriate user permissions need
to have been set up. 

With postgresql it is possible to use schema to organise metric
tables. Schemas also need to be created beforehand, they will
not be created by Daisy.
