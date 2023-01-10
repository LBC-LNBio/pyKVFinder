#!/usr/bin/env python

# System imports
import sys
import pathlib
from setuptools import setup, Extension, dist

# pyKVFinder information
from pyKVFinder import __name__, __version__

# Prepare reqs from requirements.txt
with open("requirements.txt") as f:
    reqs = f.read().splitlines()

# Get the long description from the README file
long_description = (pathlib.Path(__file__).parent.resolve() / "README.rst").read_text(
    encoding="utf-8"
)

# Third-party modules - we depend on numpy for everything
np_req = [req for req in reqs if req.find("numpy") != -1]
dist.Distribution().fetch_build_eggs(np_req)
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# Extension modules
_pyKVFinder = Extension(
    name="_pyKVFinder",
    sources=["C/pyKVFinder.i", "C/pyKVFinder.c"],
    include_dirs=[numpy_include, "C"],
    extra_compile_args=["-fopenmp", "-Ofast", "-lm"],
    extra_link_args=["-lgomp"],
)

# Setup
setup(
    name=__name__,
    version=__version__,
    description="Python package to detect and characterize cavities in biomolecular structures.",
    # This is an optional longer description of your project that represents
    # the body of text which users will see when they visit PyPI.
    long_description=long_description,
    long_description_content_type="text/x-rst",
    # This field corresponds to the "Home-Page" metadata field:
    url="https://github.com/LBC-LNBio/pyKVFinder",
    # Authors information
    author="João Victor da Silva Guerra",  # Optional
    author_email="jvsguerra@gmail.com",
    # Classifiers help users find your project by categorizing it.
    classifiers=[  # Optional
        # How mature is this project? Common values are
        #   3 - Alpha
        #   4 - Beta
        #   5 - Production/Stable
        "Development Status :: 4 - Beta",
        # Indicate who your project is intended for
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        # Pick your license as you wish
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        # Specify the Python versions you support here. In particular, ensure
        # that you indicate you support Python 3. These classifiers are *not*
        # checked by 'pip install'. See instead 'python_requires' below.
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3 :: Only",
    ],
    # Keywords
    keywords="structural biology, proteins, cavity detection, cavity characterization",
    # Extension C modules
    ext_modules=[_pyKVFinder],
    # Python package configuration
    packages=["pyKVFinder"],
    # Python versions support
    python_requires=">=3.8, <4",
    # This field lists other packages that your project depends on to run.
    install_requires=reqs,
    # If there are data files included in your packages that need to be
    # installed, specify them here.
    include_package_data=True,
    # Command-line interface
    entry_points={"console_scripts": ["pyKVFinder=pyKVFinder.main:cli"]},
    # List additional URLs that are relevant to your project as a dict.
    project_urls={  # Optional
        "Source": "https://github.com/LBC-LNBio/pyKVFinder/",
        "Issues": "https://github.com/LBC-LNBio/pyKVFinder/issues",
    },
)
