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

``pyKVFinder.pyKVFinder(pdb, ligand, dictionary, box, step = 0.6, probe_in = 1.4, probe_out = 4.0, removal_distance = 2.4, volume_cutoff = 5.0, ligand_cutoff = 5.0, include_depth = False, include_hydropathy = False, hydrophobicity_scale = 'EisenbergWeiss', surface = 'SES', ignore_backbone = False, nthreads = os.cpu_count() - 1, verbose = False)``
  Detects and characterizes cavities (volume, area, depth [optional], hydropathy [optional] and interface residues).

  :Args:
    * ``pdb`` : *str*
        A path to input PDB file
    * ``ligand`` : *str*
        A path to ligand PDB file
    * ``dictionary`` : *str, default 'vdw.dat'*
        A path to a van der Waals radii file
    * ``fn`` : *str*
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
    * ``include_depth`` : *bool, default False*
        Whether to characterize the depth of the detected cavities
    * ``include_hydropathy`` : *bool, default False*
        Whether to characterize the hydropathy of the detected cavities
    * ``hydrophobicity_scale`` : *str, default EisenbergWeiss*
        Name of a native hydrophobicity scale (EisenbergWeiss, HessaHeijne, KyteDoolitte, MoonFleming, WimleyWhite, ZhaoLondon) or a path to a TOML-formatted file with a custom hydrophobicity scale.
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
        A class that contains cavities 3D grid, surface points 3D grid, 3D grid of cavity points depth, 3D grid of surface points mapped with a hydrophobicity scale, volume, area, maximum depth and average depth, average hydropathy, and interface residues and their frequencies (residues and classes of residues) per cavity, 3D grid vertices, grid spacing and number of cavities

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
    ``(resinfo, xyzr)`` : *tuple*
        A tuple with two elements: a numpy array with residue number, chain, residue name and atom name, and a numpy array with xyz coordinates and radii values

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

``pyKVFinder.get_vertices_from_file(fn, resinfo, xyzr, step = 0.6, probe_in = 1.4, probe_out = 4.0, nthreads = os.cpu_count() - 1)``
  Gets 3D grid vertices from box configuration file, selects atoms inside custom 3D grid, define sine and cosine of 3D grid angles and define xyz grid units.

  :Args:
    * ``fn`` : *str*
        A path to a box configuration file (TOML-formatted)
    * ``resinfo`` : *numpy.ndarray*
        A numpy array with residue number, chain, residue name and atom name
    * ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz coordinates and radii values
    * ``step`` : *float, default 0.6*
        Grid spacing (A)
    * ``probe_in`` : *float, default 1.4*
        Probe In size (A)
    * ``probe_out`` : *float, default 4.0*
        Probe Out size (A)
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads

  :Returns:
    ``(vertices, resinfo, xyzr, sincos, nx, ny, nz)`` : *tuple*
        A tuple with five elements: a numpy array with xyz vertices coordinates (vertices), a numpy array with residue number, chain, residue name and atom name (resinfo), a numpy array with xyz coordinates and radii values (xyzr), a numpy array with sine and cossine of 3D grid angles (sincos), x grid units (nx), y grid units (ny) and z grid units (nz)

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
  Detects biomolecular cavities. Cavity points that belongs to the same cavity are assigned with an integer in the grid. Biomolecule points = 0, Unsigned cavity points = 1, Assigned cavity points >= 2 and Bulk points = -1.

  :Args:
    * ``nx`` : *int*
        x 3D grid units
    * ``nx`` : *int*
        y 3D grid units
    * ``nx`` : *int*
        z 3D grid units
    * ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz coordinates and radii values
    * ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)
    * ``sincos`` : *numpy.ndarray*
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

  :Note:
    * `Biomolecule point`: `0`
    * `Unassigned cavity point`: `1`
    * `Cavity point`: `>=2`
    * `Bulk point`: `-1`

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
        A tuple with three elements: numpy array with surface points of cavities (surface points >= 2; surface[nx][ny][nz]), a dictionary with volume of each detected cavity and a dictionary with area of each detected cavity

``pyKVFinder.depth(cavities, ncav, step = 0.6, nthreads = os.cpu_count() - 1, verbose = False)``
  Characterization of the depth of the detected cavities, including depth per cavity point and maximum and average depths of detected cavities.

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
    ``(depths, max_depth, avg_depth)`` : *tuple*
        A tuple with three elements: numpy array with depth of cavity points (depth[nx][ny][nz]), a dictionary with maximum depth of each detected cavity and a dictionary with average depth of each detected cavity

``pyKVFinder.constitutional(cavities, resinfo, xyzr, vertices, sincos, ncav, step = 0.6, probe_in = 1.4, ignore_backbone = False, nthreads = os.cpu_count() - 1, verbose = False)``
  Constitutional characterization (interface residues) of the detected cavities.

  :Args:
    * ``cavities`` : *numpy.ndarray*
        A numpy array with cavities (cavity points >= 2; cavities[nx][ny][nz])
    * ``resinfo`` : *numpy.ndarray*
        A numpy array with residue number, chain, residue name and atom name
    * ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz coordinates and radii values
    * ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)
    * ``sincos`` : *numpy.ndarray*
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

``pyKVFinder.hydropathy(surface, resinfo, xyzr, vertices, sincos, ncav, step = 0.6, probe_in = 1.4, hydrophobicity_scale = 'EisenbergWeiss', ignore_backbone = False, nthreads = os.cpu_count() - 1, verbose = False)``
  Hydropathy characterization of the detected cavities. Map a hydrophobicity scale per surface point and calculate average hydropathy of detected cavities.

  :Args:
    * ``surface`` : *numpy.ndarray*
        A numpy array with surface points of cavities (surface points >= 2; surface[nx][ny][nz])
    * ``resinfo`` : *numpy.ndarray*
        A numpy array with residue number, chain, residue name and atom name
    * ``xyzr`` : *numpy.ndarray*
        A numpy array with xyz coordinates and radii values
    * ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)
    * ``sincos`` : *numpy.ndarray*
        A numpy array with sine and cossine of 3D grid angles (a, b)
    * ``ncav`` : *int*
        Number of cavities in ``cavities`` numpy array
    * ``step`` : *float, default 0.6*
        Grid spacing (A)
    * ``probe_in`` : *float, default 1.4*
        Probe In size (A)
    * ``hydrophobicity_scale`` : *str, default 'EisenbergWeiss'*
        Name of a native hydrophobicity scale (EisenbergWeiss, HessaHeijne, KyteDoolitte, MoonFleming, WimleyWhite, ZhaoLondon) or a path to a TOML-formatted file with a custom hydrophobicity scale file.
    * ``ignore_backbone`` :  *bool, default False*
        Whether to ignore backbone atoms (C, CA, N, O) when defining interface residues
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads
    * ``verbose`` : *bool, default False*
        Print extra information to standard output

  :Returns:
    ``(scales, avg_hydropathy)`` : *tuple*
        A tuple with two elements: numpy array with hydrophobicity scale values mapped at surface points (scales[nx][ny][nz]) and a dictionary with average hydropathy of each detected cavity and range of the hydrophobicity scale mapped

``pyKVFinder.export(fn, cavities, surface, vertices, sincos, ncav, step = 0.6, B = None, output_hydropathy = 'hydropathy.pdb', scales = None, nthreads = os.cpu_count() - 1, append = False)``
  Exports cavities to PDB-formatted file with variable (B; optional) as B-factor, and hydropathy to PDB-formatted file as B-factor at surface points.

  :Args:
    * ``fn`` : *str*
        A path to PDB file for writing cavities
    * ``cavities`` : *numpy.ndarray*
        A numpy array with cavities (cavity points >= 2; cavities[nx][ny][nz])
    * ``surface`` : *numpy.ndarray*
        A numpy array with surface points of cavities (surface points >= 2; surface[nx][ny][nz])
    * ``vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)
    * ``sincos`` : *numpy.ndarray*
        A numpy array with sine and cossine of 3D grid angles (a, b)
    * ``ncav`` : *int*
        Number of cavities in ``cavities`` and ``surface`` numpy arrays
    * ``step`` : *float, default 0.6*
        Grid spacing (A)
    * ``B`` : *numpy.ndarray*
        B-factor for cavity points (B[nx][ny][nz])
    * ``output_hydropathy`` :  *str, default 'hydropathy.pdb'*
        A path to hydropathy PDB file (surface points mapped with a hydrophobicity scale)
    * ``scales``: *numpy.ndarray, default None*
        Hydrophobicity scale values mapped at surface points (scales[nx][ny][nz])
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads
    * ``append`` : *bool, default False*
        Append cavities to PDB file
  
  :Returns:
    A file with PDB-formatted data corresponding to cavity points (H), surface points (HA) and a target variable (B) as B-factor, and (optional) a file with PDB-formatted data corresponding to hydropathy mapped as B-factor at surface points (HA).

``pyKVFinder.calculate_frequencies(residues)``
  Calculate frequencies of residues and class of residues (R1, R2, R3, R4 and R5) for detected cavities.

  :Args:
    * ``residues`` : *dict*
        A dictionary with interface residues of each detected cavity

  :Returns:
    * ``frequencies`` : *dict*
        A dictionary with frequencies of interface residues and classes of residues of each detected cavity

  :Note:
    * ``R1`` : Alipathic apolar
        Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine
    * ``R2`` : Aromatic
        Phenylalanine, Tryptophan, Tyrosine
    * ``R3`` : Polar Uncharged
        Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine
    * ``R4`` : Negatively charged
        Aspartate, Glutamate
    * ``R5`` : Positively charged
        Arginine, Histidine, Lysine
    * ``RX`` : Non-standard
        Non-standard residues

``pyKVFinder.plot_frequencies(residues, fn = 'histograms.pdf)``
  Plot histograms of calculated frequencies (residues and classes of residues) for each detected cavity in a target PDF file.

  :Args:
    * ``frequencies`` : *dict*
        A dictionary with frequencies of interface residues and classes of residues of each detected cavity
    * ``fn`` : *str, default 'histograms.pdf'* 
        A path to PDF file for plotting histograms of frequencies.

  :Returns:
    A PDF file with histograms of calculated frequencies (residues and classes of residues) of each detected cavity.

  :Note:
    * ``R1`` : Alipathic apolar
        Alanine, Glycine, Isoleucine, Leucine, Methionine, Valine
    * ``R2`` : Aromatic
        Phenylalanine, Tryptophan, Tyrosine
    * ``R3`` : Polar Uncharged
        Asparagine, Cysteine, Glutamine, Proline, Serine, Threonine
    * ``R4`` : Negatively charged
        Aspartate, Glutamate
    * ``R5`` : Positively charged
        Arginine, Histidine, Lysine
    * ``RX`` : Non-standard
        Non-standard residues

``pyKVFinder.write_results(fn, pdb, ligand, output, output_hydropathy = None, volume = None, area = None, max_depth = None, avg_depth = None, avg_hydropathy = None, residues = None, step = 0.6)``
  Writes file paths and cavity characterization to TOML-formatted file.

  :Args:
    * ``fn`` : *str*
        A path to TOML-formatted file for writing file paths and cavity characterization (volume, area, depth [optional] and interface residues) per cavity detected
    * ``pdb`` : *str*
        A path to input PDB file
    * ``ligand`` : *str*
        A path to ligand PDB file
    * ``output`` :  *str*
        A path to cavity PDB file
    * ``output_hydropathy`` :  *str*
        A path to hydropathy PDB file (surface points mapped with a hydrophobicity scale)
    * ``volume`` : *dict*
        A dictionary with volume of each detected cavity
    * ``area`` : *dict*
        A dictionary with area of each detected cavity
    * ``max_depth`` : *dict*
        A dictionary with maximum depth of each detected cavity
    * ``avg_depth`` : *dict*
        A dictionary with average depth of each detected cavity
    * ``avg_hydropathy`` : *dict*
        A dictionary with average hydropathy of each detected cavity and range of the hydrophobicity scale mapped
    * ``residues`` : *dict*
        A dictionary with interface residues of each detected cavity
    * ``step`` : *float, default 0.6*
        Grid spacing (A)

  :Returns:
    A file with TOML-formatted data corresponding to file paths and cavity characterization per detected cavity

``pyKVFinder.pyKVFinderResults(cavities, surface, depths, scales, volume, area, max_depth, avg_depth, avg_hydropathy, residues, _vertices, _step, _ncav, _pdb = None, _ligand = None)``
  A class containing pyKVFinder detection and characterization results.

  :Attributes:
    * ``cavities`` : *numpy.ndarray*
        A numpy array with cavities (cavity points >= 2; cavities[nx][ny][nz])
    * ``surface`` : *numpy.ndarray*
        A numpy array with surface points of cavities (surface points >= 2; surface[nx][ny][nz])
    * ``depths`` : *numpy.ndarray*
        A numpy array with depth of cavity points (depth[nx][ny][nz])
    * ``scales``: *numpy.ndarray*
        Hydrophobicity scale values mapped at surface points (scales[nx][ny][nz])
    * ``volume`` : *dict*
        A dictionary with volume of each detected cavity
    * ``area`` : *dict*
        A dictionary with area of each detected cavity
    * ``max_depth`` : *dict*
        A dictionary with maximum depth of each detected cavity
    * ``avg_depth`` : *dict*
        A dictionary with average depth of each detected cavity
    * ``avg_hydropathy`` : *dict*
        A dictionary with average hydropathy of each detected cavity and range of the hydrophobicity scale mapped
    * ``residues`` : *dict*
        A dictionary with interface residues of each detected cavity
    * ``frequency`` : *dict*
        A dictionary with frequency of residues and class of residues of each detected cavity
    * ``_vertices`` : *numpy.ndarray*
        A numpy array with xyz vertices coordinates (origin, Xmax, Ymax, Zmax)
    * ``_step`` : *float, default 0.6*
        Grid spacing (A)
    * ``_ncav`` : *int*
        Number of cavities in ``cavities`` and ``surface`` numpy arrays
    * ``_pdb`` : *str*
        A path to input PDB file
    * ``_ligand`` : *str*
        A path to ligand PDB file

  :Methods:
    * ``export(output = 'cavity.pdb', output_hydropathy = 'hydropathy.pdb', nthreads = os.cpu_count() - 1)``
        Exports cavities to PDB-formatted file with variable (B; optional) as B-factor, and hydropathy to PDB-formatted file as B-factor at surface points (scales; optional)
    * ``write(fn = 'results.toml', output = None, output_hydropathy = None)``
        Writes TOML-formatted results file
    * ``plot_frequencies(pdf = 'histogram.pdf')``
        Plot histograms of frequencies in PDF file
    * ``export_all(fn = 'results.toml', output = 'cavity.pdb', output_hydropathy = 'hydropathy.pdb', include_frequencies_pdf = False, pdf = 'histogtrams.pdf', nthreads = os.cpu_count() - 1)``
        Exports cavities to PDB-formatted file with variable (B; optional) as B-factor, hydropathy to PDB-formatted file as B-factor at surface points (scales; optional), and writes TOML-formatted results file..Also includes a flag to plot histograms of frequencies (residues and classes of residues)

``pyKVFinder.pyKVFinderResults.export(output = 'cavity.pdb', output_hydropathy = 'hydropathy.pdb', nthreads = os.cpu_count() - 1)``
  Exports cavities to PDB-formatted file with variable (B; optional) as B-factor, and hydropathy to PDB-formatted file as B-factor at surface points (scales; optional).

  :Args:
    * ``output`` : *str, default 'cavity.pdb'*
        A path to PDB file for writing cavities
    * ``output_hydropathy`` : *str, default 'hydropathy.pdb'*
        A path to PDB file for writing hydropathy at surface points
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads

  :Returns:
    A file with PDB-formatted data corresponding to cavity points (H), surface points (HA) and a target variable (B) as B-factor, and (optional) a file with PDB-formatted data corresponding to hydropathy mapped as B-factor at surface points (HA).

``pyKVFinder.pyKVFinderResults.write(fn = 'results.toml, output = None, output_hydropathy = None)``
  Writes file paths and cavity characterization to TOML-formatted file.

  :Args:
    * ``fn`` : *str, default 'results.toml'*
        A path to TOML-formatted file for writing file paths and cavity characterization (volume, area, depth and interface residues) per cavity detected
    * ``output`` : *str, default None*
        A path to a cavity PDB file
    * ``output_hydropathy`` : *str, default None*
        A path to PDB file for writing hydropathy at surface points

  :Returns:
    A file with TOML-formatted data corresponding to file paths and cavity characterization per detected cavity

``pyKVFinder.plot_frequencies(pdf = 'histogram.pdf')``
  Plot histograms of frequencies in PDF file

  :Args:
    * ``pdf`` : *str, default 'histograms.pdf'*
        A path to a PDF file
  
  :Returns:
    A PDF file with histograms of calculated frequencies (residues and classes of residues) of each detected cavity.

``pyKVFinder.pyKVFinderResults.export_all(fn = 'results.toml', output = 'cavity.pdb', output_hydropathy = 'hydropathy.pdb', include_frequencies_pdf = False, pdf = 'histogtrams.pdf', nthreads = os.cpu_count() - 1)``
  Exports cavities to PDB-formatted file and writes results to TOML-formatted file.

  :Args:
    * ``fn`` : *str, default 'results.toml'*
        A path to TOML-formatted file for writing file paths and cavity characterization (volume, area and interface residues) per cavity detected
    * ``output`` : *str, default 'cavity.pdb'*
        A path to PDB file for writing cavities
    * ``output_hydropathy`` : *str, default 'hydropathy.pdb'*
        A path to PDB file for writing hydropathy at surface points
    * ``include_frequencies_pdf`` : *bool, default False*
        Whether to plot frequencies (residues and classes of residues) to PDF file
    * ``pdf`` : *str, default 'histograms.pdf'*
        A path to a PDF file
    * ``nthreads`` : *int, default 'number of cpus - 1'*
        Number of threads

  :Returns:
    A file with PDB-formatted data corresponding to cavity points (H), surface points (HA) and a target variable (B) as B-factor, (optional) a file with PDB-formatted data corresponding to hydropathy mapped as B-factor at surface points (HA), a file with TOML-formatted data corresponding to file paths and cavity characterization per detected cavity and (optinal) a PDF file with histograms of calculated frequencies (residues and classes of residues) of each detected cavity.

Van der Waals Radii File Template
=================================

The van der Waals radii file define the radius values for each residue and when not defined, it uses a generic value based on the atom type. However, the user can define its own file with a mandatory format and pass it to pyKVFinder. The format is shown below:

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
===============================

There are two methods for defining a custom 3D grid in pyKVFinder.

The first directly defines four vertices of the 3D grid (origin, X-axis, Y-axis and Z-axis), an example is shown below:

.. code-block:: TOML

  [box]
  # px = [x, y, z]
  p1 = [0.0, 0.0, 0.0]
  p2 = [1.0, 0.0, 0.0]
  p3 = [0.0, 1.0, 0.0]
  p4 = [0.0, 0.0, 1.0]


The second defines a list of residues and a padding, the template is shown below:

.. code-block:: TOML

  [box]
  residues = [ ["resname", "chain",], ["resname", "chain",], ]
  padding =  3.5


Hydrophobicity Scale File Template
==================================

The hydrophobicity scale file define the scale values for each residue and when not defined, it assigns 0.0 to missing residues. There are five native hydrophobicity scales: EisenbergWeiss, HessaHeijne, KyteDoolitte, MoonFleming, WimleyWhite and ZhaoLondon. However, the user can define its own file with a mandatory format and pass it to pyKVFinder. The format is shown below:

.. code-block:: TOML

    [EisenbergWeiss]
    ALA = -0.64
    ARG = 2.6
    ASN = 0.8
    ASP = 0.92
    CYS = -0.3
    GLN = 0.87
    GLU = 0.76
    GLY = -0.49
    HIS = 0.41
    ILE = -1.42
    LEU = -1.09
    LYS = 1.54
    MET = -0.66
    PHE = -1.22
    PRO = -0.12
    SER = 0.18
    THR = 0.05
    TRP = -0.83
    TYR = -0.27
    VAL = -1.11


Command Line Interface
======================

pyKVFinder's Command Line Interface (CLI) aims to direct the interaction between pyKVFinder and users.

.. code-block:: bash

  $ pyKVFinder
    Usage: pyKVFinder [-h] [-v] [--version] [-b <str>] [-O <str>] [--nthreads <int>] [-d <str>] [-s <float>] [-i <float>] [-o <float>] [-V <float>] [-R <float>] [-S <str>] [--ignore_backbone]
                      [-D] [--plot_frequencies] [-B <.toml>] [-L <.pdb>] [--ligand_cutoff <float>]
                      <.pdb>

The positional arguments are:

* ``<.pdb>``: A path to a target PDB file.
  
  .. code-block:: bash
    
    $ pyKVFinder <.pdb>

The optional arguments are:

* ``-h`` or ``--help``: Show help message.
  
  .. code-block:: bash
    
    $ pyKVFinder -h
    $ pyKVFinder --help

* ``--version``: Display parKVFinder version.
  
  .. code-block:: bash

    $ pyKVFinder --version

* ``-v`` or ``--verbose``: Print extra information to standard output.

  .. code-block:: bash

    $ pyKVFinder <.pdb> --verbose

:Default: ``None``

* ``-b <str>`` or ``--base_name <str>``: A prefix for output files.

  .. code-block:: bash

    $ pyKVFinder <.pdb> -b <str>
    $ pyKVFinder <.pdb> --base_name <str>

:Default: Prefix of target PDB file (<.pdb>)

* ``-O <str>`` or ``--output_directory <str>``: A path to a directory for output files.

  .. code-block:: bash

    $ pyKVFinder <.pdb> -O <str>
    $ pyKVFinder <.pdb> --output_directory <str>

:Default: Current working directory

* ``--nthreads <int>``: Number of threads to apply in parallel routines.

  .. code-block:: bash

    $ pyKVFinder <.pdb> --nthreads <int>

:Default: ``os.cpu_count() - 1``

The arguments for adjusting biomolecular detection are:

* ``-d <str>`` or ``--dictionary <str>``: A path to a van der Waals radii file (see template).

  .. code-block:: bash

    $ pyKVFinder <.pdb> -d <str>
    $ pyKVFinder <.pdb> --dictionary <str>

:Default: ``vdw.dat``

* ``-s <float>`` or ``--step <float>``: Grid spacing (A).

  .. code-block:: bash

    $ pyKVFinder <.pdb> -s <float>
    $ pyKVFinder <.pdb> --step <float>

:Default: ``0.6``

* ``-i <float>`` or ``--probe_in <float>``: Probe In size (A).

  .. code-block:: bash

    $ pyKVFinder <.pdb> -i <float>
    $ pyKVFinder <.pdb> --probe_in <float>

:Default: ``1.4``

* ``-o <float>`` or ``--probe_out <float>``: Probe Out size (A).

  .. code-block:: bash

    $ pyKVFinder <.pdb> -o <float>
    $ pyKVFinder <.pdb> --probe_out <float>

:Default: ``4.0``

* ``-V <float>`` or ``--volume_cutoff <float>``: Cavities volume filter (A3).

  .. code-block:: bash

    $ pyKVFinder <.pdb> -V <float>
    $ pyKVFinder <.pdb> --volume_cutoff <float>

:Default: ``5.0``

* ``-R <float>`` or ``--removal_distance <float>``: Length to be removed from the cavity-bulk frontier (A).

  .. code-block:: bash

    $ pyKVFinder <.pdb> -R <float>
    $ pyKVFinder <.pdb> --removal_distance <float>

:Default: ``2.4``

* ``-S <str>`` or ``--surface <str>``: A surface representation. Options are: ``SES`` and ``SAS``. SES specifies solvent excluded surface and SAS specifies solvent accessible surface.

  .. code-block:: bash

    $ pyKVFinder <.pdb> -S <str>
    $ pyKVFinder <.pdb> --surface <str>

:Default: ``SES``

* ``--ignore_backbone``: Ignore backbone contacts to cavity when defining interface residues.

  .. code-block:: bash

    $ pyKVFinder <.pdb> --ignore_backbone

:Default: ``None``

The parameters for additional characterization are:

* ``--D or --depth``: Characterize the depth of the detected cavities. This mode includes depth of each cavity point as the B-factor in the cavity PDB file and maximum and average depth of the detected cavities in the results file.

  .. code-block:: bash

    $ pyKVFinder <.pdb> -D
    $ pyKVFinder <.pdb> --depth

:Default: ``None``

* ``--plot_frequencies``: Plot histograms of calculated frequencies (residues and classes of residues) of the detected cavities in a PDF file. The classes of residues are aliphatic apolar (R1), aromatic (R2), polar uncharged (R3), negatively charged (R4), positively charged (R5) and non-standard (RX) residues.

  .. code-block:: bash

    $ pyKVFinder <.pdb> --plot_frequencies

:Default: ``None``

* ``--hydropathy [{EisenbergWeiss, HessaHeijne, KyteDoolitte, MoonFleming, WimleyWhite, ZhaoLondon, <.toml>}]``: Characterize the hydropathy of the detected cavities. This mode maps a target hydrophobicity scale as B-factor at surface points of the detected cavities. Also, it calculates the average hydropathy of each detected cavity. The constant hydrophobicity scale is EisenbergWeiss.

  .. code-block:: bash

    $ pyKVFinder <.pdb> --hydropathy

In addition, the user can define one of the native hydrophobicity scale. The native hydrophobicity scales are: EisenbergWeiss, HessaHeijne, KyteDoolitte, MoonFleming, WimleyWhite and ZhaoLondon.

  .. code-block:: bash

    $ pyKVFinder <.pdb> --hydropathy EisenbergWeiss
    $ pyKVFinder <.pdb> --hydropathy HessaHeijne
    $ pyKVFinder <.pdb> --hydropathy KyteDoolitte
    $ pyKVFinder <.pdb> --hydropathy MoonFleming
    $ pyKVFinder <.pdb> --hydropathy WimleyWhite
    $ pyKVFinder <.pdb> --hydropathy ZhaoLondon

Further, the user can also define a custom hydrophobicity scale file via a TOML-formatted file (see template).

  .. code-block:: bash

    $ pyKVFinder <.pdb> --hydropathy <.toml>

:Default: ``None``
:Constant: ``EisenbergWeiss``

The box adjustment argument is:

* ``-B <.toml>`` or ``--box <.toml>``: A path to TOML-formatted file with box parameters (see template). Adjust the 3D grid based on a list of residues (["resnum", "chain"]) and a padding or a set of four vertices (p1: origin, p2: X-axis max, p3: Y-axis max, p4: Z-axis max) with xyz coordinates ([x, y, z]).

  .. code-block:: bash

    $ pyKVFinder <.pdb> -B <.toml>
    $ pyKVFinder <.pdb> --box <.toml>

:Default: ``None``

The ligand adjustment arguments are:

* ``-L <.pdb>`` or ``--ligand <.pdb>``: A path to a ligand PDB file to limit the cavities within a radius (ligand_cutoff) around it.

  .. code-block:: bash

    $ pyKVFinder <.pdb> -L <.pdb>
    $ pyKVFinder <.pdb> --ligand <.pdb>

:Default: ``None``

* ``--ligand_cutoff <float>``: A radius value to limit a space around the defined ligand.

  .. code-block:: bash

    $ pyKVFinder <.pdb> -L <.pdb> --ligand_cutoff <float>

:Default: ``5.0``

Licensing
=========

This project is released under the terms of the GNU General Public License. View
*LICENSE.txt* for more information.
