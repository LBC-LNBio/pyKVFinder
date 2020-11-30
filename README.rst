**********
pyKVFinder
**********

.. image:: https://img.shields.io/pypi/v/pyKVFinder
    :target: https://pypi.org/project/pyKVFinder/

.. image:: https://img.shields.io/pypi/pyversions/pyKVFinder
    :target: https://pypi.org/project/pyKVFinder/


A python package for detecting and characterizing biomolecular cavities.

See also:

* `pyKVFinder GitHub repository <https://github.com/LBC-LNBio/pyKVFinder/>`_
* `pyKVFinder wiki <https://github.com/LBC-LNBio/pyKVFinder/wiki>`_

Installation
============

To install the latest release on `PyPI <https://pypi.org/project/pyKVFinder>`_, 
run:

::

  pip install pyKVFinder

Or to install the latest developmental version, run:

::

  git clone https://github.com/LBC-LNBio/pyKVFinder.git
  pip install -e pyKVFinder


API Reference
=============

``pyKVFinder.read_vdw(fn = "vdw.dat")``
  Read van der Waals radii from .dat file

  :Args:
    * ``fn`` : *str, default 'vdw.dat'*
        A path to a van der Waals radii file

  :Returns:
    ``vdw`` : *dict*
        A dictionary containing radii values (vdw[resname][atom])

``pyKVFinder.read_pdb(fn, vdw)``
  Read PDB file into numpy arrays

  :Args:
    * ``fn`` : *str*
        A path to a PDB file
    * ``vdw`` : *dict*
        A dictionary with radii values (vdw[resname][atom])

  :Returns:
    ``(pdb, xyzr)`` : *tuple*
        A tuple with two elements: a numpy array with resnum, chain, resname and atom name, and a numpy array with xyz coordinates and radii values

``pyKVFinder.get_vertices(xyzr, probe_out = 4.0, step = 0.6)``
    Get 3D grid vertices

  :Args:
    * ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz coordinates and radii values
    * ``probe_out`` : *float, default 4.0*
        Probe Out size (A)
    * ``step`` : *float, default 0.6*
        Grid spacing (A)

  :Returns:
    ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)

``pyKVFinder.get_vertices_from_file(fn, pdb, xyzr, step = 0.6, probe_in = 1.4, probe_out = 4.0, nthreads)``
    Gets 3D grid vertices from box configuration file and selects atoms inside custom 3D grid

  :Args:
    * ``fn``: *str*
        A path to a box configuration file (TOML-formatted)
    * ``pdb`` : *numpy.ndarray*
        A numpy array with resnum, chain, resname and atom name
    * ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz coordinates and radii values
    * ``step`` : *float, default 0.6*
        Grid spacing (A)
    * ``probe_in``: *float, default 1.4*
        Probe In size (A)
    * ``probe_out`` : *float, default 4.0*
        Probe Out size (A)
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads

  :Returns:
    ``(vertices, pdb, xyzr, sincos, nx, ny, nz)`` : *tuple*
        A tuple with five elements: a numpy array with xyz vertices coordinates (vertices), a numpy array with resnum, chain, resname and atom name (pdb), a numpy array with xyz coordinates and radii values (xyzr), a numpy array with sine and cossine of 3D grid angles (sincos), x grid units (nx), y grid units (ny) and z grid units (nz)

``pyKVFinder.get_dimensions(vertices, step = 0.6)``
  Gets dimensions of 3D grid from vertices

  :Args:
    * ``vertices`` (numpy.ndarray): A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)
    * ``step`` (float): Grid spacing (A)

  :Returns:
    * A tuple with three elements: x grid units (nx), y grid units (ny) and z grid units (nz)

``pyKVFinder.get_sincos(vertices)``
  Gets sine and cossine of 3D grid angles (a, b)

  :Args:
    * ``vertices`` (numpy.ndarray): A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)

  :Returns:
    A numpy array with sine and cossine of 3D grid angles (a, b)

``pyKVFinder.detect(nx, ny, nz, xyzr, vetices, sincos, step = 0.6, probe_in = 1.4, probe_out = 4.0, removal_distance = 2.4, volume_cutoff = 5.0, lxyzr = None, ligand_cutoff = 5.0, box_adjustment = False, surface = 'SES', nthreads, verbose = False)``
  Detects biomolecular cavities

  :Args:
    * nx (int): x 3D grid units
    * ny (int): y 3D grid units
    * nz (int): z 3D grid units
    * ``xyzr`` (numpy.ndarray): A numpy array with xyz coordinates and radii values
    * ``vertices`` (numpy.ndarray): A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)
    * ``sincos`` (numpy.ndarray): 
        A numpy array with sine and cossine of 3D grid angles (a, b)
    * ``step`` (float): Grid spacing (A)
    * ``probe_in`` (float): Probe In size (A)
    * ``probe_out`` (float): Probe Out size (A)
    * ``removal_distance`` ()
    * ``nthreads`` (int): Number of threads



Command Line Interface
======================






..   :Note:
..     Box Configuration File Template:
..       [box]

..       p1 = [x1, y1, z1]
      
..       p2 = [x2, y2, z2]
      
..       p3 = [x3, y3, z3]
      
..       p4 = [x4, y4, z4]


parKVFinder in Python v3
========================

This is pyKVFinder Python3 library.

Comments
--------

The surface representations are: SES (solvent excluded surface) and SAS (solvent accessible surface). Also, vdW representation could also be achieved by using Probe In of 0.0 Angstrom; however, not useful for cavity detection.