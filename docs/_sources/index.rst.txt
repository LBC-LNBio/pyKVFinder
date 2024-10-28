######################################
Welcome to pyKVFinder's documentation!
######################################

Welcome to the **Python-C parallel KVFinder (pyKVFinder)** documentation, this page was built to help you get started with our cavity detection and characterization package.

**pyKVFinder** is an open-source Python package designed for detecting and characterizing cavities in biomolecular structures. It employs a geometric grid-and-sphere-based method with a dual-probe system to detect cavities and provides detailed information about their spatial, depth, constitutional, and hydropathy characteristics. In addition to cavity detection, **pyKVFinder** also estimates the molecular volume, using van der Waals (vdW) surface, solvent excluded surface (SES), and solvent accessible surface (SAS) to represent the molecular surface. 

.. Headings:
..  # with overline, for parts
..  * with overline, for chapters
..  = for sections
..  - for subsections
..  ^ for subsubsections
..  " for paragraphs

.. toctree::
   :maxdepth: 1
   :caption: Python Package

   Installation <package/installation/index>
   Tutorial <package/tutorial/index>
   API Reference <package/api_reference/index>
   Examples <package/examples/index>
   GitHub repository <https://github.com/LBC-LNBio/pyKVFinder>

.. toctree::
   :maxdepth: 1
   :caption: Plugins

   Command-line interface <plugins/cli/index>
   PyMOL pyKVFinder Tools <plugins/pymol/index>
   ChimeraX <plugins/chimerax/index>

.. toctree::
   :maxdepth: 1
   :caption: About

   about/index
