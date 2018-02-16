===================
Daisy Documentation
===================

.. image:: daisies.png

Daisy is a framework and toolkit to perform computational analyses
efficiently and reproducibly. It's aim is to build:

* Abstract computational envirorment to achieve portability,
* Provide logging and runtime metrics,
* Provide configuration and parameterization,
* Provide Consistent CLI for tools and workflows,
* Provide templates for tools and pipelines

.. image:: daisy_components.png

Daisy is based on a foundation of powerful libraries for data
handling, data management, job execution, logging. On top of this
foundation sits the the daisy engine that uses these libraries to
provide a consistent environment for users to build and execute
workflows and execute tasks.

Daisy applications are built using the daisy engine. The applications
are stand-alone tools (:ref:`Toolkit`) with a specific functionality
or might be full workflows (:ref:`Daisy chains`) running daisy tools
or 3rd-part applications

This documentation provides usage instruction for all of Daisy's
components.

Contents:

.. toctree::
   :maxdepth: 2

   Installation.rst
   DaisyChains.rst
   Tools.rst
   Engine.rst
   Development.rst
   Tutorials.rst
   Glossary.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

