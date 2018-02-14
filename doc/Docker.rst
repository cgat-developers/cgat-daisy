======
Docker
======

Building the daisy docker image
===============================

To build a container, run::

   cd daisy
   scripts/make-docker-container.sh daisy

Using the daisy docker container
================================

Running a tool with data from stdin::

  echo -e "abc\n2\n" |
  docker run -i daisy-docker-0.2.0-0:latest \
     "daisy table2stats"

Note that the daisy command needs to be put into quotes.
     
Running a tool with data on file system::

  docker run -v /local/scratch/project/ONT/na12878-chr22:/data \
  daisy-docker-0.2.0-0:latest \
  "daisy bam2stats /data/NA12878-Xten_remapped-chr22.bam"

Getting an interactive shell within the daisy container::

  docker run -it daisy-docker-0.2.0-0:latest "source activate daisy && /bin/bash"
