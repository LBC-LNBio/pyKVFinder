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


Usage examples
==============

For an complete example of usage, we will display usage examples with a pdb from [here](https://github.com/LBC-LNBio/pyKVFinder/blob/master/tests/data/1FMO.pdb).

First of all, import pyKVFinder module on python:

.. code-block:: python

  >>> import pyKVFinder

Standard pipeline
-----------------

The standard pipeline for cavity detection with spatial and constitutional characterization (volume, area and interface residues) can be run at once with one command:

.. code-block:: python

  >>> results = pyKVFinder.pyKVFinder('1FMO.pdb')
  >>> results
  <pyKVFinderResults object>

Inside the `pyKVFinderResults object`, cavity and surface points, volume, area, and interface residues and their frequencies are stored as attributes. Below, we show how to access them:

.. code-block:: python
  
  >>> results.cavities
  array([[[-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1],
          ...,
          [-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1]],

        ...,

        [[-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1],
          ...,
          [-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)
  >>> results.surface
  array([[[-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1],
          ...,
          [-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1]],

        ...,

        [[-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1],
          ...,
          [-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1],
          [-1, -1, -1, ..., -1, -1, -1]]], dtype=int32)
  >>> results.volume
  {'KAA': 137.16, 'KAB': 47.52, 'KAC': 66.96, 'KAD': 8.21, 'KAE': 43.63, 'KAF': 12.53, 'KAG': 6.26, 'KAH': 520.13, 'KAI': 12.31, 'KAJ': 26.57, 'KAK': 12.31, 'KAL': 33.91, 'KAM': 23.11, 'KAN': 102.82, 'KAO': 6.05, 'KAP': 15.55, 'KAQ': 7.99, 'KAR': 7.78}
  >>> results.area
  {'KAA': 120.52, 'KAB': 58.76, 'KAC': 72.06, 'KAD': 17.62, 'KAE': 56.44, 'KAF': 22.53, 'KAG': 15.38, 'KAH': 489.25, 'KAI': 29.87, 'KAJ': 44.85, 'KAK': 30.58, 'KAL': 43.59, 'KAM': 45.25, 'KAN': 129.77, 'KAO': 11.57, 'KAP': 24.8, 'KAQ': 12.59, 'KAR': 15.97}
  >>> results.residues
  {'KAA': [['14', 'E', 'SER'], ['15', 'E', 'VAL'], ['18', 'E', 'PHE'], ['19', 'E', 'LEU'], ['100', 'E', 'PHE'], ['152', 'E', 'LEU'], ['155', 'E', 'GLU'], ['156', 'E', 'TYR'], ['292', 'E', 'LYS'], ['302', 'E', 'TRP'], ['303', 'E', 'ILE'], ['306', 'E', 'TYR']], 'KAB': [['18', 'E', 'PHE'], ['22', 'E', 'ALA'], ['25', 'E', 'ASP'], ['26', 'E', 'PHE'], ['29', 'E', 'LYS'], ['97', 'E', 'ALA'], ['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['156', 'E', 'TYR']], 'KAC': [['141', 'E', 'PRO'], ['142', 'E', 'HIS'], ['144', 'E', 'ARG'], ['145', 'E', 'PHE'], ['148', 'E', 'ALA'], ['299', 'E', 'THR'], ['300', 'E', 'THR'], ['305', 'E', 'ILE'], ['310', 'E', 'VAL'], ['311', 'E', 'GLU'], ['313', 'E', 'PRO']], 'KAD': [['122', 'E', 'TYR'], ['124', 'E', 'ALA'], ['176', 'E', 'GLN'], ['318', 'E', 'PHE'], ['320', 'E', 'GLY'], ['321', 'E', 'PRO'], ['322', 'E', 'GLY'], ['323', 'E', 'ASP']], 'KAE': [['95', 'E', 'LEU'], ['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['100', 'E', 'PHE'], ['103', 'E', 'LEU'], ['104', 'E', 'VAL'], ['105', 'E', 'LYS'], ['106', 'E', 'LEU']], 'KAF': [['123', 'E', 'VAL'], ['124', 'E', 'ALA'], ['175', 'E', 'ASP'], ['176', 'E', 'GLN'], ['181', 'E', 'GLN']], 'KAG': [['34', 'E', 'SER'], ['37', 'E', 'THR'], ['96', 'E', 'GLN'], ['106', 'E', 'LEU'], ['107', 'E', 'GLU'], ['108', 'E', 'PHE'], ['109', 'E', 'SER']], 'KAH': [['49', 'E', 'LEU'], ['50', 'E', 'GLY'], ['51', 'E', 'THR'], ['52', 'E', 'GLY'], ['53', 'E', 'SER'], ['54', 'E', 'PHE'], ['55', 'E', 'GLY'], ['56', 'E', 'ARG'], ['57', 'E', 'VAL'], ['70', 'E', 'ALA'], ['72', 'E', 'LYS'], ['74', 'E', 'LEU'], ['84', 'E', 'GLN'], ['87', 'E', 'HIS'], ['88', 'E', 'THR'], ['91', 'E', 'GLU'], ['104', 'E', 'VAL'], ['120', 'E', 'MET'], ['121', 'E', 'GLU'], ['122', 'E', 'TYR'], ['123', 'E', 'VAL'], ['127', 'E', 'GLU'], ['166', 'E', 'ASP'], ['168', 'E', 'LYS'], ['170', 'E', 'GLU'], ['171', 'E', 'ASN'], ['173', 'E', 'LEU'], ['183', 'E', 'THR'], ['184', 'E', 'ASP'], ['186', 'E', 'GLY'], ['187', 'E', 'PHE'], ['201', 'E', 'THR'], ['327', 'E', 'PHE']], 'KAI': [['131', 'E', 'HIS'], ['138', 'E', 'PHE'], ['142', 'E', 'HIS'], ['146', 'E', 'TYR'], ['174', 'E', 'ILE'], ['314', 'E', 'PHE']], 'KAJ': [['33', 'E', 'PRO'], ['89', 'E', 'LEU'], ['92', 'E', 'LYS'], ['93', 'E', 'ARG'], ['96', 'E', 'GLN'], ['349', 'E', 'GLU'], ['350', 'E', 'PHE']], 'KAK': [['157', 'E', 'LEU'], ['162', 'E', 'LEU'], ['163', 'E', 'ILE'], ['164', 'E', 'TYR'], ['185', 'E', 'PHE'], ['188', 'E', 'ALA']], 'KAL': [['49', 'E', 'LEU'], ['127', 'E', 'GLU'], ['129', 'E', 'PHE'], ['130', 'E', 'SER'], ['326', 'E', 'ASN'], ['327', 'E', 'PHE'], ['328', 'E', 'ASP'], ['330', 'E', 'TYR']], 'KAM': [['51', 'E', 'THR'], ['55', 'E', 'GLY'], ['56', 'E', 'ARG'], ['73', 'E', 'ILE'], ['74', 'E', 'LEU'], ['75', 'E', 'ASP'], ['115', 'E', 'ASN'], ['335', 'E', 'ILE'], ['336', 'E', 'ARG']], 'KAN': [['165', 'E', 'ARG'], ['166', 'E', 'ASP'], ['167', 'E', 'LEU'], ['199', 'E', 'CYS'], ['200', 'E', 'GLY'], ['201', 'E', 'THR'], ['204', 'E', 'TYR'], ['205', 'E', 'LEU'], ['206', 'E', 'ALA'], ['209', 'E', 'ILE'], ['219', 'E', 'VAL'], ['220', 'E', 'ASP'], ['223', 'E', 'ALA']], 'KAO': [['48', 'E', 'THR'], ['51', 'E', 'THR'], ['56', 'E', 'ARG'], ['330', 'E', 'TYR'], ['331', 'E', 'GLU']], 'KAP': [['222', 'E', 'TRP'], ['238', 'E', 'PHE'], ['253', 'E', 'GLY'], ['254', 'E', 'LYS'], ['255', 'E', 'VAL'], ['273', 'E', 'LEU']], 'KAQ': [['207', 'E', 'PRO'], ['208', 'E', 'GLU'], ['211', 'E', 'LEU'], ['213', 'E', 'LYS'], ['275', 'E', 'VAL'], ['277', 'E', 'LEU']], 'KAR': [['237', 'E', 'PRO'], ['238', 'E', 'PHE'], ['249', 'E', 'LYS'], ['254', 'E', 'LYS'], ['255', 'E', 'VAL'], ['256', 'E', 'ARG']]}
  >>> results.frequencies
  {'KAA': {'RESIDUES': {'GLU': 1, 'ILE': 1, 'LEU': 1, 'LYS': 1, 'PHE': 1, 'SER': 1, 'TRP': 1, 'TYR': 1, 'VAL': 1}, 'CLASS': {'R1': 3, 'R2': 3, 'R3': 1, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAB': {'RESIDUES': {'ALA': 1, 'ASN': 1, 'ASP': 1, 'LYS': 1, 'PHE': 1, 'TYR': 1, 'VAL': 1}, 'CLASS': {'R1': 2, 'R2': 2, 'R3': 1, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAC': {'RESIDUES': {'ALA': 1, 'ARG': 1, 'GLU': 1, 'HIS': 1, 'ILE': 1, 'PHE': 1, 'PRO': 1, 'THR': 1, 'VAL': 1}, 'CLASS': {'R1': 4, 'R2': 1, 'R3': 1, 'R4': 1, 'R5': 2, 'RX': 0}}, 'KAD': {'RESIDUES': {'ALA': 1, 'ASP': 1, 'GLN': 1, 'GLY': 1, 'PHE': 1, 'PRO': 1, 'TYR': 1}, 'CLASS': {'R1': 3, 'R2': 2, 'R3': 1, 'R4': 1, 'R5': 0, 'RX': 0}}, 'KAE': {'RESIDUES': {'ASN': 1, 'LEU': 1, 'LYS': 1, 'PHE': 1, 'VAL': 1}, 'CLASS': {'R1': 2, 'R2': 1, 'R3': 1, 'R4': 0, 'R5': 1, 'RX': 0}}, 'KAF': {'RESIDUES': {'ALA': 1, 'ASP': 1, 'GLN': 1, 'VAL': 1}, 'CLASS': {'R1': 2, 'R2': 0, 'R3': 1, 'R4': 1, 'R5': 0, 'RX': 0}}, 'KAG': {'RESIDUES': {'GLN': 1, 'GLU': 1, 'LEU': 1, 'PHE': 1, 'SER': 1, 'THR': 1}, 'CLASS': {'R1': 1, 'R2': 1, 'R3': 3, 'R4': 1, 'R5': 0, 'RX': 0}}, 'KAH': {'RESIDUES': {'ALA': 1, 'ARG': 1, 'ASN': 1, 'ASP': 1, 'GLN': 1, 'GLU': 1, 'GLY': 1, 'HIS': 1, 'LEU': 1, 'LYS': 1, 'MET': 1, 'PHE': 1, 'SER': 1, 'THR': 1, 'TYR': 1, 'VAL': 1}, 'CLASS': {'R1': 4, 'R2': 2, 'R3': 5, 'R4': 2, 'R5': 3, 'RX': 0}}, 'KAI': {'RESIDUES': {'HIS': 1, 'ILE': 1, 'PHE': 1, 'TYR': 1}, 'CLASS': {'R1': 1, 'R2': 2, 'R3': 0, 'R4': 0, 'R5': 1, 'RX': 0}}, 'KAJ': {'RESIDUES': {'ARG': 1, 'GLN': 1, 'GLU': 1, 'LEU': 1, 'LYS': 1, 'PHE': 1, 'PRO': 1}, 'CLASS': {'R1': 2, 'R2': 1, 'R3': 1, 'R4': 1, 'R5': 2, 'RX': 0}}, 'KAK': {'RESIDUES': {'ALA': 1, 'ILE': 1, 'LEU': 1, 'PHE': 1, 'TYR': 1}, 'CLASS': {'R1': 3, 'R2': 2, 'R3': 0, 'R4': 0, 'R5': 0, 'RX': 0}}, 'KAL': {'RESIDUES': {'ASN': 1, 'ASP': 1, 'GLU': 1, 'LEU': 1, 'PHE': 1, 'SER': 1, 'TYR': 1}, 'CLASS': {'R1': 1, 'R2': 2, 'R3': 2, 'R4': 2, 'R5': 0, 'RX': 0}}, 'KAM': {'RESIDUES': {'ARG': 1, 'ASN': 1, 'ASP': 1, 'GLY': 1, 'ILE': 1, 'LEU': 1, 'THR': 1}, 'CLASS': {'R1': 3, 'R2': 0, 'R3': 2, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAN': {'RESIDUES': {'ALA': 1, 'ARG': 1, 'ASP': 1, 'CYS': 1, 'GLY': 1, 'ILE': 1, 'LEU': 1, 'THR': 1, 'TYR': 1, 'VAL': 1}, 'CLASS': {'R1': 5, 'R2': 1, 'R3': 2, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAO': {'RESIDUES': {'ARG': 1, 'GLU': 1, 'THR': 1, 'TYR': 1}, 'CLASS': {'R1': 0, 'R2': 1, 'R3': 1, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAP': {'RESIDUES': {'GLY': 1, 'LEU': 1, 'LYS': 1, 'PHE': 1, 'TRP': 1, 'VAL': 1}, 'CLASS': {'R1': 3, 'R2': 2, 'R3': 0, 'R4': 0, 'R5': 1, 'RX': 0}}, 'KAQ': {'RESIDUES': {'GLU': 1, 'LEU': 1, 'LYS': 1, 'PRO': 1, 'VAL': 1}, 'CLASS': {'R1': 3, 'R2': 0, 'R3': 0, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAR': {'RESIDUES': {'ARG': 1, 'LYS': 1, 'PHE': 1, 'PRO': 1, 'VAL': 1}, 'CLASS': {'R1': 2, 'R2': 1, 'R3': 0, 'R4': 0, 'R5': 2, 'RX': 0}}}

With these attributes, we can write the detected cavities and the characterization to files. Further, we can set a flag to plot the histograms of the frequencies in a PDF file. Below, we illustrate the usage:

.. code-block:: python

  >>> results.export_all(fn='results.toml', output='cavity.pdb', include_frequencies_pdf=True, pdf='histograms.pdf')

We note that the above implementation of the pyKVFinder function uses default parameter specifications, and that therefore parameters can be adjusted to users’ needs.

For more information on `pyKVFinder function` and `pyKVFinderResults object` attributes and methods, please refer to the API Reference.

Full pipeline
-------------

The full pipeline for cavity dectection and full characterization (volume, area, depth, hydropathy and interface residues) can be run with one command by setting some parameters of ``pyKVFinder`` function:

.. code-block:: python

  >>> results = pyKVFinder.pyKVFinder('tests/data/1FMO.pdb', include_depth=True, include_hydropathy=True,  hydrophobicity_scale='EisenbergWeiss')

Inside the `pyKVFinderResults object`, in addition to cavity and surface points, volume, area, and interface residues and their frequencies showed above, depth and hydropathy points, average depth, maximum depth and average hydropathy are also stored as attributes. Below, we show how to access them:

.. code-block:: python

  >>> results.depths
  array([[[0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.],
          ...,
          [0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.]],

        ...,

        [[0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.],
          ...,
          [0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.]]])
  >>> results.scales
  array([[[0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.],
          ...,
          [0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.]],

        ...,

        [[0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.],
          ...,
          [0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.],
          [0., 0., 0., ..., 0., 0., 0.]]])
  >>> results.avg_depth
  {'KAA': 1.28, 'KAB': 0.86, 'KAC': 0.67, 'KAD': 0.29, 'KAE': 0.98, 'KAF': 0.24, 'KAG': 0.1, 'KAH': 3.75, 'KAI': 1.5, 'KAJ': 0.96, 'KAK': 0.0, 'KAL': 1.0, 'KAM': 0.24, 'KAN': 0.0, 'KAO': 0.29, 'KAP': 0.7, 'KAQ': 0.22, 'KAR': 0.12}
  >>> results.max_depth
  {'KAA': 3.79, 'KAB': 2.68, 'KAC': 2.62, 'KAD': 0.85, 'KAE': 3.0, 'KAF': 0.85, 'KAG': 0.6, 'KAH': 10.73, 'KAI': 2.55, 'KAJ': 2.24, 'KAK': 0.0, 'KAL': 3.0, 'KAM': 1.2, 'KAN': 0.0, 'KAO': 1.04, 'KAP': 2.08, 'KAQ': 0.85, 'KAR': 0.6}
  >>> results.avg_hydropathy
  {'KAA': -0.73, 'KAB': -0.06, 'KAC': -0.07, 'KAD': -0.62, 'KAE': -0.81, 'KAF': -0.14, 'KAG': -0.33, 'KAH': -0.16, 'KAI': -0.4, 'KAJ': 0.62, 'KAK': -0.99, 'KAL': 0.35, 'KAM': -0.33, 'KAN': 0.18, 'KAO': 0.88, 'KAP': -0.96, 'KAQ': 0.48, 'KAR': 0.24, 'EisenbergWeiss': [-1.42, 2.6]}

With these attributes, we can write the detected cavities with depth annotated on B-factor, the surface cavity points with hydrophobicity scale annotated on B-factor, and the characterization to files. Below, we illustrate the usage:

.. code-block:: python

  >>> results.export_all(fn='results.toml', output='cavity.pdb', output_hydropathy = 'hydropathy.pdb', include_frequencies_pdf=False)

We note that the above implementation of the pyKVFinder function uses default parameter specifications, except for include_depth and include_hydropathy parameters, and that therefore parameters can be adjusted to users’ needs.

For more information on `pyKVFinder function` and `pyKVFinderResults object` attributes and methods, please refer to the API Reference.

Separated steps
---------------

If you prefer, instead of running ``pyKVFinder`` function, you can apply the cavity detection and characterization in a step-by-step fashion. Below we describe each function in detail.

1. Loading van der Waals radii dictionary

The van der Waals radii file define the radius values for each residue and when not defined, it uses a generic value based on the atom type. ``pyKVFinder.read_vdw`` takes a `.dat` file and returns a dictionary contaning radii values for each atom of each residue.

.. code-block:: python

  >>> vdw = read_vdw()
  >>> vdw
  {'ALA': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB1': 1.487, '1HB': 1.487, 'HB2': 1.487, '2HB': 1.487, 'HB3': 1.487, '3HB': 1.487, 'C': 1.908, 'O': 1.6612}, 'ARG': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'HD2': 1.387, '1HD': 1.387, '2HD': 1.387, 'HD3': 1.387, 'HD1': 1.387, 'NE': 1.75, 'HE': 0.6, 'CZ': 1.908, 'NH1': 1.75, 'HH11': 0.6, '1HH1': 0.6, 'HH12': 0.6, '2HH1': 0.6, 'NH2': 1.75, 'HH21': 0.6, '2HH2': 0.6, 'HH22': 0.6, '1HH2': 0.6, 'C': 1.908, 'O': 1.6612}, 'ASH': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'OD1': 1.6612, 'OD2': 1.721, 'HD2': 0.0001, 'C': 1.908, 'O': 1.6612}, 'ASN': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'OD1': 1.6612, 'ND2': 1.824, 'HD21': 0.6, '1HD2': 0.6, 'HD22': 0.6, '2HD2': 0.6, 'C': 1.908, 'O': 1.6612}, 'ASP': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'OD1': 1.6612, 'OD2': 1.6612, 'C': 1.908, 'O': 1.6612}, 'CYM': {'N': 1.824, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB3': 1.387, 'HB2': 1.387, 'SG': 2.0, 'C': 1.908, 'O': 1.6612}, 'CYS': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, '2HB': 1.387, '1HB': 1.387, 'HB3': 1.387, 'HB1': 1.387, 'SG': 2.0, 'HG': 0.6, 'C': 1.908, 'O': 1.6612}, 'CYX': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, 'HB3': 1.387, 'SG': 2.0, 'C': 1.908, 'O': 1.6612}, 'GLH': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'HG2': 1.487, 'HG3': 1.487, 'CD': 1.908, 'OE1': 1.6612, 'OE2': 1.721, 'HE2': 0.0001, 'C': 1.908, 'O': 1.6612}, 'GLN': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'OE1': 1.6612, 'NE2': 1.824, 'HE21': 0.6, '1HE2': 0.6, 'HE22': 0.6, '2HE2': 0.6, 'C': 1.908, 'O': 1.6612}, 'GLU': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'OE1': 1.6612, 'OE2': 1.6612, 'C': 1.908, 'O': 1.6612}, 'GLY': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA2': 1.387, 'HA1': 1.387, '1HA': 1.387, '2HA': 1.387, 'HA3': 1.387, 'C': 1.908, 'O': 1.6612}, 'HID': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'HD1': 0.6, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'HIE': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'HE2': 0.6, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'HIP': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'HD1': 0.6, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'HE2': 0.6, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'ILE': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.487, 'CG2': 1.908, 'HG21': 1.487, '1HG2': 1.487, 'HG22': 1.487, '2HG2': 1.487, 'HG23': 1.487, '3HG2': 1.487, 'CG1': 1.908, 'HG12': 1.487, '2HG1': 1.487, 'HG13': 1.487, 'HG11': 1.487, '1HG1': 1.487, 'CD1': 1.908, 'HD11': 1.487, '1HD1': 1.487, 'HD12': 1.487, '2HD1': 1.487, 'HD13': 1.487, '3HD1': 1.487, 'C': 1.908, 'O': 1.6612}, 'LEU': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG': 1.487, 'CD1': 1.908, 'HD11': 1.487, '1HD1': 1.487, 'HD12': 1.487, '2HD1': 1.487, 'HD13': 1.487, '3HD1': 1.487, 'CD2': 1.908, 'HD21': 1.487, '1HD2': 1.487, 'HD22': 1.487, '2HD2': 1.487, 'HD23': 1.487, '3HD2': 1.487, 'C': 1.908, 'O': 1.6612}, 'LYN': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'HG2': 1.487, 'HG3': 1.487, 'CD': 1.908, 'HD2': 1.487, 'HD3': 1.487, 'CE': 1.908, 'HE2': 1.1, 'HE3': 1.1, 'NZ': 1.824, 'HZ2': 0.6, 'HZ3': 0.6, 'C': 1.908, 'O': 1.6612}, 'LYS': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CD': 1.908, 'HD2': 1.487, '1HD': 1.487, '2HD': 1.487, 'HD3': 1.487, 'HD1': 1.487, 'CE': 1.908, 'HE2': 1.1, '2HE': 1.1, 'HE3': 1.1, '1HE': 1.1, 'HE1': 1.1, 'NZ': 1.824, 'HZ1': 0.6, '1HZ': 0.6, 'HZ2': 0.6, '2HZ': 0.6, 'HZ3': 0.6, '3HZ': 0.6, 'C': 1.908, 'O': 1.6612}, 'MET': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'HG2': 1.387, '2HG': 1.387, 'HG3': 1.387, 'HG1': 1.387, '1HG': 1.387, 'SD': 2.0, 'CE': 1.908, 'HE1': 1.387, '1HE': 1.387, 'HE2': 1.387, '2HE': 1.387, 'HE3': 1.387, '3HE': 1.387, 'C': 1.908, 'O': 1.6612}, 'PHE': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'HZ': 1.459, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'C': 1.908, 'O': 1.6612}, 'PRO': {'N': 1.824, 'CD': 1.908, 'HD2': 1.387, '1HD': 1.387, '2HD': 1.387, 'HD3': 1.387, 'HD1': 1.387, 'CG': 1.908, 'HG2': 1.487, '2HG': 1.487, 'HG3': 1.487, 'HG1': 1.487, '1HG': 1.487, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CA': 1.908, 'HA': 1.387, 'C': 1.908, 'O': 1.6612}, 'SER': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, '2HB': 1.387, '1HB': 1.387, 'HB3': 1.387, 'HB1': 1.387, 'OG': 1.721, 'HG': 0.0001, 'C': 1.908, 'O': 1.6612}, 'THR': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.387, 'CG2': 1.908, 'HG21': 1.487, '1HG2': 1.487, 'HG22': 1.487, '2HG2': 1.487, 'HG23': 1.487, '3HG2': 1.487, 'OG1': 1.721, 'HG1': 0.0001, 'C': 1.908, 'O': 1.6612}, 'TRP': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.85, 'CD1': 2.0, 'HD1': 1.409, 'NE1': 1.75, 'HE1': 0.6, 'CE2': 1.85, 'CZ2': 1.908, 'HZ2': 1.459, 'CH2': 1.908, 'HH2': 1.459, 'CZ3': 1.908, 'HZ3': 1.459, 'CE3': 1.908, 'HE3': 1.459, 'CD2': 1.85, 'C': 1.908, 'O': 1.6612}, 'TYR': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'OH': 1.721, 'HH': 0.0001, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'C': 1.908, 'O': 1.6612}, 'VAL': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.487, 'CG1': 1.908, 'CG2': 1.908, 'HG11': 1.487, '1HG2': 1.487, '1HG1': 1.487, 'HG21': 1.487, 'HG12': 1.487, '2HG1': 1.487, 'HG22': 1.487, '2HG2': 1.487, 'HG13': 1.487, '3HG2': 1.487, '3HG1': 1.487, 'HG23': 1.487, 'C': 1.908, 'O': 1.6612}, 'HIS': {'N': 1.824, 'H': 0.6, 'HN': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, '2HB': 1.487, '1HB': 1.487, 'HB3': 1.487, 'HB1': 1.487, 'CG': 1.85, 'ND1': 1.75, 'HD1': 0.6, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'CD2': 2.0, 'HD2': 1.409, 'C': 1.908, 'O': 1.6612}, 'PTR': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'OH': 1.6837, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'SEP': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, 'HB3': 1.387, '1HB': 1.387, '2HB': 1.387, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'TPO': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.387, 'CG2': 1.908, 'HG21': 1.487, 'HG22': 1.487, 'HG23': 1.487, '1HG2': 1.487, '2HG2': 1.487, '3HG2': 1.487, 'OG1': 1.6837, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'H2D': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.85, 'ND1': 1.75, 'CE1': 1.85, 'HE1': 1.359, 'NE2': 1.75, 'HE2': 0.6, 'CD2': 2.0, 'HD2': 1.409, 'P': 2.1, 'O1P': 1.85, 'O2P': 1.85, 'O3P': 1.85, 'C': 1.908, 'O': 1.6612}, 'Y1P': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.487, 'HB3': 1.487, 'CG': 1.908, 'CD1': 1.908, 'HD1': 1.459, 'CE1': 1.908, 'HE1': 1.459, 'CZ': 1.908, 'CE2': 1.908, 'HE2': 1.459, 'CD2': 1.908, 'HD2': 1.459, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.721, 'O2P': 1.6612, 'O3P': 1.6612, 'H1P': 0.0001, 'C': 1.908, 'O': 1.6612}, 'T1P': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB': 1.387, 'CG2': 1.908, 'HG21': 1.487, 'HG22': 1.487, 'HG23': 1.487, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.721, 'O2P': 1.6612, 'O3P': 1.6612, 'H1P': 0.0001, 'C': 1.908, 'O': 1.6612}, 'S1P': {'N': 1.824, 'H': 0.6, 'CA': 1.908, 'HA': 1.387, 'CB': 1.908, 'HB2': 1.387, 'HB3': 1.387, 'OG': 1.6837, 'P': 2.1, 'O1P': 1.721, 'O2P': 1.6612, 'O3P': 1.6612, 'H1P': 0.0001, 'C': 1.908, 'O': 1.6612}, 'GEN': {'AC': 2.0, 'AG': 1.72, 'AL': 2.0, 'AM': 2.0, 'AR': 1.88, 'AS': 1.85, 'AT': 2.0, 'AU': 1.66, 'B': 2.0, 'BA': 2.0, 'BE': 2.0, 'BH': 2.0, 'BI': 2.0, 'BK': 2.0, 'BR': 1.85, 'C': 1.66, 'CA': 2.0, 'CD': 1.58, 'CE': 2.0, 'CF': 2.0, 'CL': 1.75, 'CM': 2.0, 'CO': 2.0, 'CR': 2.0, 'CS': 2.0, 'CU': 1.4, 'DB': 2.0, 'DS': 2.0, 'DY': 2.0, 'ER': 2.0, 'ES': 2.0, 'EU': 2.0, 'F': 1.47, 'FE': 2.0, 'FM': 2.0, 'FR': 2.0, 'GA': 1.87, 'GD': 2.0, 'GE': 2.0, 'H': 0.91, 'HE': 1.4, 'HF': 2.0, 'HG': 1.55, 'HO': 2.0, 'HS': 2.0, 'I': 1.98, 'IN': 1.93, 'IR': 2.0, 'K': 2.75, 'KR': 2.02, 'LA': 2.0, 'LI': 1.82, 'LR': 2.0, 'LU': 2.0, 'MD': 2.0, 'MG': 1.73, 'MN': 2.0, 'MO': 2.0, 'MT': 2.0, 'N': 1.97, 'NA': 2.27, 'NB': 2.0, 'ND': 2.0, 'NE': 1.54, 'NI': 1.63, 'NO': 2.0, 'NP': 2.0, 'O': 1.69, 'OS': 2.0, 'P': 2.1, 'PA': 2.0, 'PB': 2.02, 'PD': 1.63, 'PM': 2.0, 'PO': 2.0, 'PR': 2.0, 'PT': 1.72, 'PU': 2.0, 'RA': 2.0, 'RB': 2.0, 'RE': 2.0, 'RF': 2.0, 'RH': 2.0, 'RN': 2.0, 'RU': 2.0, 'S': 2.09, 'SB': 2.0, 'SC': 2.0, 'SE': 1.9, 'SG': 2.0, 'SI': 2.1, 'SM': 2.0, 'SN': 2.17, 'SR': 2.0, 'TA': 2.0, 'TB': 2.0, 'TC': 2.0, 'TE': 2.06, 'TH': 2.0, 'TI': 2.0, 'TL': 1.96, 'TM': 2.0, 'U': 1.86, 'V': 2.0, 'W': 2.0, 'XE': 2.16, 'Y': 2.0, 'YB': 2.0, 'ZN': 1.39, 'ZR': 2.0}}

.. note::

  The function takes the `builtin dicitonary <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/vdw.dat>`_ when a *.dat* file is not specified. Otherwise, user must specify a *.dat* file following `Van der Waals Radii Template <https://github.com/LBC-LNBio/pyKVFinder#van-der-waals-radii-file-template>`_.


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
        Exports cavities to PDB-formatted file with variable (B; optional) as B-factor, hydropathy to PDB-formatted file as B-factor at surface points (scales; optional), and writes TOML-formatted results file.Also includes a flag to plot histograms of frequencies (residues and classes of residues)

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

  :Default: ``False``

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
