# System imports
from distutils.core import setup, Extension

# parKVFinder information
name = "parKVFinder"
version = "1.1"

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# Extension modules
_grid = Extension(
    name="_gridprocessing",
    sources=["src/grid.i", "src/grid.c"],
    include_dirs=[numpy_include],
    extra_compile_args=['-fopenmp', '-lm'],
    extra_link_args=['-lgomp'],
    # swig_opts=['-c']
)

# Prepare reqs from requirements.txt
with open('requirements.txt') as f:
    reqs = f.read().splitlines()

# Setup
setup(
    name=name, 
    version=version,
    ext_modules=[_grid],
    include_package_data=True,
    install_requires=reqs,
    packages=['pyKVFinder']
    # include_dirs=['src'],
    # extra_compile_args=['-fopenmp'],
    # extra_link_args=['-lgomp'],
    # swig_opts=['-threads']
)

# setup(
#     name=name, 
#     version=version,
#     # distutils detects .i files and compiles them automatically
#     ext_modules=[
#         Extension(
#             name='_parKVFinder', # SWIG requires _ as a prefix for the module name
#             sources=[
#                 "parKVFinder.i", 
#                 "src/parKVFinder.c",
#                 "src/matrixprocessing.c",
#                 "src/pdbprocessing.c",
#                 "src/resultsprocessing.c",
#                 "src/tomlprocessing.c",
#                 "src/dictionaryprocessing.c",
#                 "src/argparser.c"
#                 ],
#             include_dirs=['src'],
#             extra_compile_args=["-fopenmp"],
#             extra_link_args=['-lgomp'],
#             swig_opts=['-threads']
#         )
#     ]
#     # entry_points={
#     #     'console_scripts': [
#     #         'parKVFinder=parKVFinder.parKVFinder:run',
#     #     ],
#     # },
# )