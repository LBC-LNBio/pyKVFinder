# System imports
from setuptools import setup, Extension, dist

# pyKVFinder information
from pyKVFinder import __name__, __version__

# Prepare reqs from requirements.txt
with open('requirements.txt') as f:
    reqs = f.read().splitlines()

# Third-party modules - we depend on numpy for everything
np_req = [req for req in reqs if req.find('numpy') != -1]
dist.Distribution().fetch_build_eggs(np_req)
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# Extension modules
_grid = Extension(
    name="_grid",
    sources=["src/grid.i", "src/grid.c"],
    include_dirs=[numpy_include, 'src'],
    extra_compile_args=['-fopenmp', '-Ofast', '-lm'],
    extra_link_args=['-lgomp'],
)

# Setup
setup(
    name=__name__,
    version=__version__,
    ext_modules=[_grid],
    include_package_data=True,
    install_requires=reqs,
    packages=['pyKVFinder'],
    entry_points={'console_scripts': ['pyKVFinder=pyKVFinder.main:cli']},
)
