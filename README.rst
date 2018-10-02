.. image:: https://travis-ci.org/cgat-developers/cgat-daisy.svg?branch=master
    :target: https://travis-ci.org/cgat-developers/cgat-daisy

================
Daisy/CGAT Bench
================

Daisy is a system to design and execute benchmarking tasks.

After installation, use the ``daisy`` command to see how to use them.

For questions, please open a discussion on the GitHub 
`issue page <https://github.com/cgat-developers/cgat-daisy/issues>`_.

Installation
============

Install using Conda
-------------------

The preferred method to install Daisy is using the installation script, which uses
`Conda <https://conda.io>`_.

Here are the steps::

        # download installation script:
        curl -O https://raw.githubusercontent.com/cgat-developers/cgat-daisy/master/install-CGAT-tools.sh

        # see help:
        bash install-CGAT-tools.sh

        # install the development version (recommended, no production version yet):
        bash install-CGAT-tools.sh --devel [--location </full/path/to/folder/without/trailing/slash>]

        # the code is downloaded in zip format by default. If you want to get a git clone, use:
        --git # for an HTTPS clone
        --git-ssh # for a SSH clone (you need to be a cgat-developer contributor on GitHub to do this)

        # enable the conda environment as requested by the installation script
        # NB: you probably want to automate this by adding the instructions below to your .bashrc
        source </full/path/to/folder/without/trailing/slash>/conda-install/etc/profile.d/conda.sh
        conda activate base
        conda activate cgat-a

        # finally, please run the cgatflow command-line tool to check the installation:
        cgat --help

The installation script will put everything under the specified location. The aim is to provide a portable
installation that does not interfere with the existing software. As a result, you will get a conda environment
working with CGAT Apps which can be enabled on demand according to your needs.

Usage
=====

Run the ``daisy --help`` for the tools that are available.

Documentation
=============

Full documentation is on RTD (WIP)


.. _cgat-core: https://github.com/cgat-developers/cgat-core
