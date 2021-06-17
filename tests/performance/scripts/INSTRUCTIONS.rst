HOW TO RUN BENCHMARKING ?
=========================

Requirements
------------
* Python 3
    * ``biopython``: ``pip install biopython``
    * ``pandas``: ``pip install pandas``
* R
    * ``ggplot2``: install.packages('ggplot2')
    * ``tidyr``: install.packages('tidyr')
    * ``ggpmisc``: install.packages('ggpmisc')
    * ``RColorBrewer``: install.packages('RColorBrewer')
* parKVFinder
* pyKVFinder

Steps
-----

1. Download and compile parKVFinder software

``./compile_parKVFinder.sh``

:Note:
  * Define path of ``pyKVFinder/tests/performance`` directory in ``OUTPUT_DIRECTORY``;
  * Define path to install parKVFinder repository in ``PARKVFINDER_LOCATION``;
  * Check and tune list of OpenMP threads to be applied in parKVFider in nthreads array ``nthreads``.

2. Run benchmarking script: Execute parKVFinder and pyKVFinder with different number of threads and compile results (time and speedup).

``./runbenchmarking.sh OUTPUT_DIRECTORY``

:Example: 
  ``./runbenchmarking.sh ../``: Running runbenchmarking.sh script from ``scripts`` directory.

:Note:
  * Define path of ``pyKVFinder/tests/performance`` directory in ``OUTPUT_DIRECTORY``;
  * Check and tune list of threads to be applied in parKVFider in nthreads array ``nthreads``.

Output
------

The benchmarking generates: 

* `raw` directory: ``pyKVFinder/tests/performance/raw``
  Contains the time to execute parKVFinder and pyKVFinder (default, depth, hydropathy) with different number of threads (i.e. 1, 2, 4, 8, 12, 16, 20 and 24 OpenMP threads) for kv1000 proteins with 3 replicates each. The files are organized in two directories: ``parKVFinder``, ``pyKVFinder``, ``pyKVFinder-depth``, ``pyKVFinder-hydropathy``. Each one contain CSV files with the pdb, average time and standard deviation for different number of threads (i.e. 01.csv, 02.csv, 04.csv, 08.csv, 12.csv, 16.csv, 20.csv and 24.csv).

* ``data`` directory: ``pyKVFinder/tests/performance/data``
  Contains processed time data and plots of the benchmarking.

  * ``time.csv``: file contains the processed data from ``raw`` directory.

    :header:
      ``,pdb,natoms,parKVFinder_01,parKVFinder_02,parKVFinder_04,parKVFinder_08,parKVFinder_12,parKVFinder_16,parKVFinder_20,parKVFinder_24,pyKVFinder_01,pyKVFinder_02,pyKVFinder_04,pyKVFinder_08,pyKVFinder_12,pyKVFinder_16,pyKVFinder_20,pyKVFinder_24,pyKVFinder-depth_01,pyKVFinder-depth_02,pyKVFinder-depth_04,pyKVFinder-depth_08,pyKVFinder-depth_12,pyKVFinder-depth_16,pyKVFinder-depth_20,pyKVFinder-depth_24,pyKVFinder-hydropathy_01,pyKVFinder-hydropathy_02,pyKVFinder-hydropathy_04,pyKVFinder-hydropathy_08,pyKVFinder-hydropathy_12,pyKVFinder-hydropathy_16,pyKVFinder-hydropathy_20,pyKVFinder-hydropathy_24``

  * ``plots`` directory: contains time, speedup and effiency plots of the benchmarking. 

    * ``comparison``: contains time, speedup and effiency plots of pyKVFinder (default, depth, hydropathy) software against parKVFinder.
      
      :metrics:
        * ``speedup``: ``parKVFinder`` time with ``x`` thread divided by ``pyKVFinder`` time with ``x`` threads
        * ``effiency``: ``speedup```of ``x`` threads divided by ``x`` threads

    * ``pyKVFinder`` directory: contains time, speedup and effiency plots of pyKVFinder (default) software against itself.
      
      :metrics:
        * ``speedup``: ``pyKVFinder`` time with 1 thread divided by ``pyKVFinder`` time with ``x`` threads
        * ``effiency``: ``speedup```of ``x`` threads divided by ``x`` threads

    * ``pyKVFinder-depth`` directory: contains time, speedup and effiency plots of pyKVFinder (depth) software against pyKVFinder (default).
      
      :metrics:
        * ``speedup``: ``pyKVFinder --depth`` time with ``x`` thread divided by ``pyKVFinder`` time with ``x`` threads
        * ``effiency``: ``speedup```of ``x`` threads divided by ``x`` threads

    * ``pyKVFinder-hydropathy`` directory: contains time, speedup and effiency plots of pyKVFinder (hydropathy) software against pyKVFinder (default).
      
      :metrics:
        * ``speedup``: ``pyKVFinder --hydropathy`` time with ``x`` thread divided by ``pyKVFinder`` time with  threads
        * ``effiency``: ``speedup```of ``x`` threads divided by ``x`` threads

* ``kv1000`` directory: ``pyKVFinder/tests/performance/kv1000``
  Contains 1000 proteins of the kv1000 dataset (https://github.com/jvsguerra/kv1000.git)
