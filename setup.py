import sys
from setuptools import Extension, setup


class get_numpy_include(object):
    def __str__(self):
        import numpy

        return numpy.get_include()


setup(
    ext_modules=[
        Extension(
            name="_pyKVFinder",
            sources=["C/pyKVFinder.i", "C/pyKVFinder.c"],
            include_dirs=[get_numpy_include(), "C"],
            extra_compile_args=["-fopenmp", "-Ofast", "-lm"],
            extra_link_args=["-lgomp", "-static"]
            if sys.platform != "linux"
            else ["-lgomp"],
        ),
    ]
)
