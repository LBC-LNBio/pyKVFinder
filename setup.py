import os
import sys
from setuptools import Extension, setup
import numpy


class GetNumpyInclude:
    def __str__(self):
        return numpy.get_include()


def get_extra_link_args():
    if sys.platform == "darwin":
        arch = os.uname().machine
        if arch == "arm64":
            extra_link_args = ["-L/opt/homebrew/opt/libomp/lib", "-lomp"]
        else:
            extra_link_args = ["-L/usr/local/opt/libomp/lib", "-lomp"]
    elif sys.platform == "linux":
        extra_link_args = ["-lgomp"]
    elif sys.platform == "win32":
        extra_link_args = ["/openmp", "/O2"]
    else:
        extra_link_args = []
    return extra_link_args


def get_extra_compile_args():
    if sys.platform == "darwin":
        extra_compile_args = ["-Xpreprocessor", "-fopenmp=libomp", "-O3", "-ffast-math"]
    elif sys.platform == "linux":
        extra_compile_args = ["-fopenmp", "-Ofast"]
    elif sys.platform == "win32":
        extra_compile_args = []
    else:
        extra_compile_args = []
    return extra_compile_args


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
