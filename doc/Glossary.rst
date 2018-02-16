=========
Glossary
=========

.. glossary::

   tool
       A task running a tool. Tools typically perform data
       transformation, for example a variant caller derives variant
       calls from reads aligned in a :term:`bam` file and outputs
       them in a :term:`vcf` file.

   metric
      A task computing a metric on a file. The output is one or
      more tab- separated files in :term:`tsv` format.

   bam
      File format for storing aligned short-read data.

   vcf
     File format for storing variant calls. 

   fasta
     File format for storing sequence data.

   tsv
     Tab separated file format. A :term:`tsv` contains a header and
     rows. Comments start with `#` at the beginning of a line.

   runners
   runner
      Python function responsible for running a :term:`metric` or :term:`tool`.
