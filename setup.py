import sys
import setuptools
from setuptools import setup, find_packages

from distutils.version import LooseVersion
if LooseVersion(setuptools.__version__) < LooseVersion('1.1'):
    print("Version detected:", LooseVersion(setuptools.__version__))
    raise ImportError(
        "the CGAT code collection requires setuptools 1.1 higher")

# collect version
sys.path.insert(0, "daisy")
from version import __version__ as daisy_version

daisy_packages = find_packages()
daisy_package_dirs = {'daisy': 'daisy'}

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
    name='cgat-daisy',
    version=daisy_version,
    description='daisy : data i-something system',
    author='Andreas Heger',
    author_email='andreas.heger@gmail.com',
    license="MIT",
    platforms=["any"],
    keywords="computational genomics",
    long_description='daisy : data ',
    classifiers=[_f for _f in classifiers.split("\n") if _f],
    url="http://www.cgat.org/cgat/Tools/",
    packages=daisy_packages,
    package_dir=daisy_package_dirs,
    include_package_data=True,
    entry_points={
        'console_scripts': ['daisy = daisy.tools.cli:main'],
    },
    zip_safe=False,
    test_suite="tests",
)
