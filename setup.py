import sys
from setuptools import Extension, setup
import numpy


class GetNumpyInclude:
    def __str__(self):
        return numpy.get_include()


def get_extra_link_args():
    if sys.platform == "darwin":
        return ["-L/opt/homebrew/opt/libomp/lib", "-lomp"]
    elif sys.platform != "linux":
        return ["-lgomp", "-static"]
    return ["-lgomp"]


def get_extra_compile_args():
    if sys.platform == "darwin":
        return ["-fopenmp", "-O3", "-march=native", "-ffast-math", "-DNDEBUG", "-Xpreprocessor"]
    return ["-fopenmp", "-Ofast", "-lm"]


setup(
    ext_modules=[
        Extension(
            name="_pyKVFinder",
            sources=["C/pyKVFinder.i", "C/pyKVFinder.c"],
            include_dirs=[str(GetNumpyInclude()), "C"],
            extra_compile_args=get_extra_compile_args(),
            extra_link_args=get_extra_link_args(),
        ),
    ]
)
