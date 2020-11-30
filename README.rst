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

``pyKVFinder.pyKVFinder(pdb, ligand, dictionary, box, step = 0.6, probe_in = 1.4, probe_out = 4.0, removal_distance = 2.4, volume_cutoff = 5.0, ligand_cutoff = 5.0, surface = 'SES', ignore_backbone = False, nthreads = os.cpu_count() - 1, verbose = False)``
  Detects and characterizes cavities (volume, area and interface residues)

  :Args:
    * ``pdb`` : *str*
        A path to input PDB file
    * ``ligand`` : *str*
        A path to ligand PDB file
    * ``dictionary`` : *str, default 'vdw.dat'*
        A path to a van der Waals radii file
    * ``fn``: *str*
        A path to a box configuration file (TOML-formatted)
    * ``step`` : *float, default 0.6*
        Grid spacing (A)
    * ``probe_in`` : *float, default 1.4*
        Probe In size (A)
    * ``probe_out`` : *float, default 4.0*
        Probe Out size (A)
    * ``removal_distance`` : *float, default 2.4*
        Length to be removed from the cavity-bulk frontier (A)
    * ``volume_cutoff`` : *float, default 5.0*
        Cavities volume filter (A3)        
    * ``ligand_cutoff`` : *float, default 5.0*
        Radius value to limit a space around a ligand (A)
    * ``surface`` : *str, default 'SES'*
        SES (Solvent Excluded Surface) or SAS (Solvent Accessible Surface)
    * ``ignore_backbone`` :  *bool, default False*
        Whether to ignore backbone atoms (C, CA, N, O) when defining interface residues
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads
    * ``verbose`` : *bool, default False*
        Print extra information to standard output

  :Returns:
    ``pyKVFinderResults`` : *object*
        A class that contains cavities and surface points 3D grids, volume, area and interface residues per cavity, 3D grid vertices, grid spacing and number of cavities

``pyKVFinder.read_vdw(fn = "vdw.dat")``
  Reads van der Waals radii from .dat file.

  :Args:
    * ``fn`` : *str, default 'vdw.dat'*
        A path to a van der Waals radii file

  :Returns:
    ``vdw`` : *dict*
        A dictionary containing radii values (vdw[resname][atom])

``pyKVFinder.read_pdb(fn, vdw)``
  Reads PDB file into numpy arrays.

  :Args:
    * ``fn`` : *str*
        A path to a PDB file
    * ``vdw`` : *dict*
        A dictionary with radii values (vdw[resname][atom])

  :Returns:
    ``(pdb, xyzr)`` : *tuple*
        A tuple with two elements: a numpy array with resnum, chain, resname and atom name, and a numpy array with xyz coordinates and radii values

``pyKVFinder.get_vertices(xyzr, probe_out = 4.0, step = 0.6)``
  Gets 3D grid vertices.

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

``pyKVFinder.get_vertices_from_file(fn, pdb, xyzr, step = 0.6, probe_in = 1.4, probe_out = 4.0, nthreads = os.cpu_count() - 1)``
  Gets 3D grid vertices from box configuration file, selects atoms inside custom 3D grid, define sine and cosine of 3D grid angles and define xyz grid units.

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
  Gets dimensions of 3D grid from vertices.

  :Args:
    * ``vertices`` : * numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)
    * ``step`` : *float*
        Grid spacing (A)

  :Returns:
    ``(nx, ny, nz)`` : *tuple*
        A tuple with three elements: x grid units (nx), y grid units (ny) and z grid units (nz)

``pyKVFinder.get_sincos(vertices)``
  Gets sine and cossine of 3D grid angles (a, b).

  :Args:
    * ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)

  :Returns:
    ``sincos`` : *numpy.ndarray*
        A numpy array with sine and cossine of 3D grid angles (a, b)

``pyKVFinder.detect(nx, ny, nz, xyzr, vetices, sincos, step = 0.6, probe_in = 1.4, probe_out = 4.0, removal_distance = 2.4, volume_cutoff = 5.0, lxyzr = None, ligand_cutoff = 5.0, box_adjustment = False, surface = 'SES', nthreads = os.cpu_count() - 1, verbose = False)``
  Detects biomolecular cavities.

  :Args:
    * ``nx`` : *int*
        x 3D grid units
    * ``nx`` : *int*
        y 3D grid units
    * ``nx`` : *int*
        z 3D grid units
    * ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz coordinates and radii values
    * ``vertices``: *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)
    * ``sincos``: *numpy.ndarray*
        A numpy array with sine and cossine of 3D grid angles (a, b)
    * ``step`` : *float, default 0.6*
        Grid spacing (A)
    * ``probe_in`` : *float, default 1.4*
        Probe In size (A)
    * ``probe_out`` : *float, default 4.0*
        Probe Out size (A)
    * ``removal_distance`` : *float, default 2.4*
        Length to be removed from the cavity-bulk frontier (A)
    * ``volume_cutoff`` : *float, default 5.0*
        Cavities volume filter (A3)
    * ``lxyzr`` : *numpy.ndarray*
        A numpy array with xyz coordinates and radii values of ligand atoms
    * ``ligand_cutoff`` : *float, default 5.0*
        Radius value to limit a space around a ligand (A)
    * ``box_adjustment`` :  *bool, default False*
        Whether a custom 3D grid is applied
    * ``surface`` : *str, default 'SES'*
        SES (Solvent Excluded Surface) or SAS (Solvent Accessible Surface)
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads
    * ``verbose`` : *bool, default False*
        Print extra information to standard output

  :Returns:
    ``(ncav, cavities)`` : *tuple*
        A tuple with two elements: number of cavities detected (ncav) and a numpy array with cavities (cavity points >= 2; cavities[nx][ny][nz])

``pyKVFinder.spatial(cavities, ncav, step = 0.6, nthreads = os.cpu_count() - 1, verbose = False)``
  Spatial characterization (volume and area) of the detected cavities.

  :Args:
    * ``cavities`` : *numpy.ndarray*
        A numpy array with cavities (cavity points >= 2; cavities[nx][ny][nz])
    * ``ncav`` : *int*
        Number of cavities in ``cavities`` numpy array
    * ``step`` : *float, default 0.6*
        Grid spacing (A)
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads
    * ``verbose`` : *bool, default False*
        Print extra information to standard output

  :Returns:
    ``(surface, volume, area)`` : *tuple*
        A tuple with three elements:  numpy array with surface points of cavities (surface points >= 2; surface[nx][ny][nz]), a dictionary with volume of each detected cavity and a dictionary with area of each detected cavity

``pyKVFinder.constitutional(cavities, pdb, xyzr, vertices, sincos, ncav, step = 0.6, probe_in = 1.4, ignore_backbone = False, nthreads = os.cpu_count() - 1, verbose = False)``
  Constitutional characterization (interface residues) of the detected cavities

  :Args:
    * ``cavities`` : *numpy.ndarray*
        A numpy array with cavities (cavity points >= 2; cavities[nx][ny][nz])
    * ``pdb`` : *numpy.ndarray*
        A numpy array with resnum, chain, resname and atom name
    * ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz coordinates and radii values
    * ``vertices``: *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)
    * ``sincos``: *numpy.ndarray*
        A numpy array with sine and cossine of 3D grid angles (a, b)
    * ``ncav`` : *int*
        Number of cavities in ``cavities`` numpy array
    * ``step`` : *float, default 0.6*
        Grid spacing (A)
    * ``probe_in`` : *float, default 1.4*
        Probe In size (A)
    * ``ignore_backbone`` :  *bool, default False*
        Whether to ignore backbone atoms (C, CA, N, O) when defining interface residues
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads
    * ``verbose`` : *bool, default False*
        Print extra information to standard output

  :Returns:
    ``residues`` : *dict*
        A dictionary with cavity name/list of interface residues pairs

``pyKVFinder.export(fn, cavities, surface, vertices, sincos, ncav, step = 0.6, nthreads = os.cpu_count() - 1)``
  Exports cavities to PDB-formatted file.

  :Args:
    * ``fn``: *str*
        A path to PDB file for writing cavities
    * ``cavities`` : *numpy.ndarray*
        A numpy array with cavities (cavity points >= 2; cavities[nx][ny][nz])
    * ``surface`` : *numpy.ndarray*
        A numpy array with surface points of cavities (cavity points >= 2; cavities[nx][ny][nz])
    * ``vertices``: *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)
    * ``sincos``: *numpy.ndarray*
        A numpy array with sine and cossine of 3D grid angles (a, b)
    * ``ncav`` : *int*
        Number of cavities in ``cavities`` and ``surface`` numpy arrays
    * ``step`` : *float, default 0.6*
        Grid spacing (A)
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads
  
  :Returns:
    A file with PDB-formatted data corresponding to cavity points

``pyKVFinder.write_results(fn, pdb, ligand, output, volume = None, area = None, residues = None, step = 0.6)``
  Writes file paths and cavity characterization to TOML-formatted file

  :Args:
    * ``fn``: *str*
        A path to TOML-formatted file for writing file paths and cavity characterization (volume, area and interface residues) per cavity detected
    * ``pdb`` : *str*
        A path to input PDB file
    * ``ligand`` : *str*
        A path to ligand PDB file
    * ``output`` :  *str*
        A path to cavity PDB file
    * ``volume`` : *dict*
        A dictionary with volume of each detected cavity
    * ``area`` : *dict*
        A dictionary with area of each detected cavity
    * ``residues`` : *dict*
        A dictionary with interface residues of each detected cavity
    * ``step`` : *float, default 0.6*
        Grid spacing (A)

  :Returns:
    A file with TOML-formatted data corresponding to file paths and cavity characterization per detected cavity

``pyKVFinder.pyKVFinderResults(cavities, surface, volume, area, residues, _vertices, _step, _ncav)``
  A class containing pyKVFinder detection and characterization results.

  :Attributes:
    * ``cavities`` : *numpy.ndarray*
        A numpy array with cavities (cavity points >= 2; cavities[nx][ny][nz])
    * ``surface`` : *numpy.ndarray*
        A numpy array with surface points of cavities (cavity points >= 2; cavities[nx][ny][nz])
    * ``volume`` : *dict*
        A dictionary with volume of each detected cavity
    * ``area`` : *dict*
        A dictionary with area of each detected cavity
    * ``residues`` : *dict*
        A dictionary with interface residues of each detected cavity
    * ``_vertices``: *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)
    * ``_step`` : *float, default 0.6*
        Grid spacing (A)
    * ``_ncav`` : *int*
        Number of cavities in ``cavities`` and ``surface`` numpy arrays

  :Methods:
    * ``export(fn = 'cavity.pdb', nthreads = os.cpu_count() - 1)``
        Exports cavities to PDB-formatted file
    * ``write(fn = 'results.toml')``
        Writes file paths and cavity characterization to TOML-formatted file
    * export_all(fn = 'results.toml', output = 'cavity.pdb', nthreads = os.cpu_count() - 1)
        Exports cavities to PDB-formatted file and writes results to TOML-formatted file

``pyKVFinder.pyKVFinderResults.export(fn = 'cavity.pdb', nthreads = os.cpu_count() - 1)``
  Exports cavities to PDB-formatted file

  :Args:
    * ``fn`` : *str, default 'cavity.pdb'*
        A path to PDB file for writing cavities
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads

  :Returns:
    A file with TOML-formatted data corresponding to file paths and cavity characterization per detected cavity

``pyKVFinder.pyKVFinderResults.write(fn = 'results.toml)``
  Writes file paths and cavity characterization to TOML-formatted file

  :Args:
    * ``fn`` : *str, default 'results.toml'*
        A path to TOML-formatted file for writing file paths and cavity characterization (volume, area and interface residues) per cavity detected

  :Returns:
    A file with TOML-formatted data corresponding to file paths and cavity characterization per detected cavity

``pyKVFinder.pyKVFinderResults.export_all(fn = 'results.toml', output = 'cavity.pdb', nthreads = os.cpu_count() - 1)``
  Exports cavities to PDB-formatted file and writes results to TOML-formatted file

  :Args:
    * ``fn`` : *str, default 'results.toml'*
        A path to TOML-formatted file for writing file paths and cavity characterization (volume, area and interface residues) per cavity detected
    * ``output`` : *str, default 'cavity.pdb'*
        A path to PDB file for writing cavities
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads

  :Returns:
    A file with PDB-formatted data corresponding to cavity points and a file with TOML-formatted data corresponding to file paths and cavity characterization per detected cavity

Van der Waals Radii File Template
=================================

The van der Waals radii file define the radius values for each residue and when not defined, it uses a generic value based on the atom type. However, the user can define its own file with a mandatory format and pass it to pyKVFinder. The format is shown above:

.. code-block::

  >RES
  C       1.66
  CA      2.00
  N       1.97
  O       1.69
  H       0.91

  :Note:
    The residue name should be in the standard PDB format and each radius value is separated by two tab characters of the atom name.

Box Configuration File Template
================================

There are two methods for defining a custom 3D grid in pyKVFinder.

The first directly defines four vertices of the 3D grid (origin, X-axis, Y-axis and Z-axis), the template is shown above:

.. code-block:: TOML

  [box]
  p1 = [x1, y1, z1]
  p2 = [x2, y2, z2]
  p3 = [x3, y3, z3]
  p4 = [x4, y4, z4]


The second defines a list of residues and a padding, the template is shown above:

.. code-block:: TOML

  [box]
  residues = [ ["resname", "chain",], ["resname", "chain",], ]
  padding =  3.5


Command Line Interface
======================

pyKVFinder Command Line Interface (CLI) aims to direct interaction between pyKVFinder and users.

.. code-block:: bash

  $ pyKVFinder
  Usage: pyKVFinder [-h] [-v] [--version] [-b <str>] [-O <path>] [--nthreads <int>] [-d <file>] [-s <float>] [-i <float>] [-o <float>] [-V <float>] [-R <float>] [-S <enum>]
                   [--ignore_backbone] [-B <.toml>] [-L <.pdb>] [--ligand_cutoff <float>]
                   <.pdb> 


The available options (short and long options) are shown above:

  * ``-h`` or ``--help``: Show help message.
    .. code-block:: bash

      pyKVFinder -h
      pyKVFinder --help

  * ``-v`` or ``--version``: Display parKVFinder version.
    .. code-block:: bash

      pyKVFinder -v
      pyKVFinder --version


