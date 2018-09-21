import sysconfig
import sys
import os
import subprocess
import re
import setuptools
from setuptools import setup, find_packages

########################################################################
#######################################################################
# Check for dependencies
#
# Is there a way to do this more elegantly?
#     1. Run "pip install numpy"
#     2. Wrap inside functions (works for numpy/pysam, but not cython)
try:
    import numpy
except ImportError:
    raise ImportError(
        "the CGAT code collection requires numpy to be installed "
        "before running setup.py (pip install numpy)")

try:
    import pysam
except ImportError:
    raise ImportError(
        "the CGAT code collection requires pysam to "
        "be installed before running setup.py (pip install pysam)")

from distutils.version import LooseVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    print("Version detected:", LooseVersion(setuptools.__version__))
    raise ImportError(
        "the CGAT code collection requires setuptools 1.1 higher")

########################################################################
########################################################################
IS_OSX = sys.platform == 'darwin'

########################################################################
########################################################################
# collect version
sys.path.insert(0, "daisy")
from version import __version__ as daisy_version

daisy_packages = find_packages()
daisy_package_dirs = {'daisy': 'daisy'}

##########################################################
##########################################################
# Classifiers
classifiers = """
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved
Programming Language :: Python
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

setup(
    # package information
    name='daisy',
    version=daisy_version,
    description='daisy : data',
    author='Andreas Heger',
    author_email='andreas.heger@gmail.com',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description='daisy : data ',
    classifiers=[_f for _f in classifiers.split("\n") if _f],
    url="http://www.cgat.org/cgat/Tools/",
    # package contents
    packages=daisy_packages,
    package_dir=daisy_package_dirs,
    include_package_data=True,
    entry_points={
        'console_scripts': ['daisy = daisy.tools.cli:main']
    },
    # other options
    zip_safe=False,
    test_suite="tests",
)
