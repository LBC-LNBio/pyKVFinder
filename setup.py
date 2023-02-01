from setuptools import Extension, setup

setup(
    ext_modules=[
        Extension(
            name="_pyKVFinder",
            sources=["C/pyKVFinder.i", "C/pyKVFinder.c"],
            include_dirs=["C"],
            extra_compile_args=["-fopenmp", "-Ofast", "-lm"],
            extra_link_args=["-lgomp"],
        ),
    ]
)
