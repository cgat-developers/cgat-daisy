==========
Deployment
==========

This section describes the steps to build the daisy deployables and
how to test and install them.

Building locally
================

Stage 1 - Building source and deployment tar balls
--------------------------------------------------

Stage 1 builds a source tar-ball and deployment scripts tar ball.

To deploy the benchmark framework, execute the following commands::

    export TEST_PROCESSES=4 && tree && ./daisy/daisy/scripts/test-and-package.sh

This will create three tar-balls::

    daisy-0.2.0-0-tools.tar.gz
    daisy-0.2.0-0-deployment-scripts.tar.gz
    daisy-0.2.0-0-standalone.tar.gz

The ``standalone`` tar-ball contains the source of the daisy system
including all dependencies. This tar-ball can be installed without
access to artifactory.

The ``deployment-scripts`` tar-ball contains the scripts to deploy the
build on OpenStack.

Installing a release tar-ball
-----------------------------

Installation through artifactory
++++++++++++++++++++++++++++++++

To install a release tar-ball, untar the tools tar-ball and point
it towards an articactory or conda channel. The channel can be on
artifactory or it can be local repository. For example, to test the
deployable locally after the deployables have been built:

   cd ..
   mkdir testdir
   tar xzf ../daisy/daisy-0.2.0-0-tools.tar.gz
   daisy/install ../daisy/conda-bld/linux-64

   source daisy/activate

Installation of stand-alone deployable
++++++++++++++++++++++++++++++++++++++


Stage 3 - Running accepance tests manually
------------------------------------------

To run the acceptance tests manually, you will need the tests tar-ball. Untar, configure and install::

   tar -xvzf daisy/daisy-0.0.0-tests.tar.gz 
   cd deployable-test/daisy-tests
   ./configure-main && ./install-main

Stage 3 - Running accepance tests
---------------------------------

To run the acceptance tests manually, you will need the tests tar-ball. Untar, configure and install::

   tar -xvzf daisy/daisy-0.0.0-tests.tar.gz 
   cd deployable-test/daisy-tests
   ./configure-main && ./install-main
   scripts/acceptance-test.sh

Stage 3 - Running regression tests
-----------------------------------

To run the regression tests manually, you will need the tests tar-ball
and the library tar-ball. Untar, configure and install::

   tar -xvzf daisy/daisy-0.0.0-tests.tar.gz 
   tar -xvzf daisy/daisy-1-960a57f-library.tar.gz
   ln -s ../daisy-library deployable-test/
   cd deployable-test/daisy-tests
   ./configure-main && ./install-main
   scripts/regression-test.sh

OpenStack
=========

Stage 2 - Building a release on OpenStack
-----------------------------------------  

To build the release on open-stack, untar the deployment-scripts tar ball
and call the script :file:`start-release` giving it the target release version
and the source tar ball as arguments::

   tar -xvzf daisy-...-deployable.tar.gz --strip-components=1
   export OS_USERNAME=andreas; export OS_PASSWORD="password"; deployment/start-release.sh 1.0 daisy-1-ab853f49...-source.tar.gz

Stage 3 - Running tests
-----------------------

To run the acceptance and regression tests, use::

   tar -xvzf daisy-...-deployable.tar.gz --strip-components=1
   export OS_USERNAME=andreas; export OS_PASSWORD="password"; deployment/start-acceptance.sh 1.0 1.0 daisy-1-ab853f49...-tests.tar.gz

For regression testing, you will need to provide both the tests deployable as well as
the library deployable::

   export OS_USERNAME=andreas; export OS_PASSWORD="password"; deployment/start-regression.sh 1.0 1.0 daisy-1-ab853f49...-tests.tar.gz daisy-1-ab853f49-library.tar.gz

