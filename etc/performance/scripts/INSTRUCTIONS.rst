HOW TO RUN BENCHMARKING ?
=========================

Requirements
------------
* Python 3
    * ``biopython``: ``pip install biopython``
    * ``pandas``: ``pip install pandas``
* R
    * ``ggplot2``: install.packages('ggplot2')
* parKVFinder
* pyKVFinder

Steps
-----

1. Download and compile parKVFinder software

``./benchmarking/compile_parKVFinder.sh``

:Note: 
  * Adjust list of threads to be applied in parKVFider in nthreads array ``nthreads``;
  * Set path of ``performance`` directory in ``CWD``;
  * Set path to install parKVFinder repository in ``PARKVFINDER_INSTALLATION``.

2. Run benchmarking script: Execute parKVFinder and pyKVFinder with different number of threads and compile results (time and speedup).

``./runbenchmarking.sh``

:Note:
  * Adjust list of threads to be applied in parKVFider in nthreads array ``nthreads``;
  * Set path of ``performance`` directory in ``CWD``;
  * Set path to install parKVFinder repository in ``PARKVFINDER_INSTALLATION``.

Glossary
--------

* Speedup 1 
* Speedup 2
* Speedup 3
