pyKVFinder.write_results
========================

Writes file paths and cavity characterization to TOML-formatted file.

.. code-block:: python

    pyKVFinder.write_results(fn, pdb, ligand, output, output_hydropathy = None, volume = None, area = None, max_depth = None, avg_depth = None, avg_hydropathy = None, residues = None, frequencies = None, step = 0.6)  

:Args:

    ``fn`` : *str*
        A path to TOML-formatted file for writing file paths and cavity characterization (volume, area, depth [optional] and interface residues) per cavity detected
    ``pdb`` : *str*
        A path to input PDB file
    ``ligand`` : *str*
        A path to ligand PDB file
    ``output`` :  *str*
        A path to cavity PDB file
    ``output_hydropathy`` :  *str*
        A path to hydropathy PDB file (surface points mapped with a hydrophobicity scale)
    ``volume`` : *dict*
        A dictionary with volume of each detected cavity
    ``area`` : *dict*
        A dictionary with area of each detected cavity
    ``max_depth`` : *dict*
        A dictionary with maximum depth of each detected cavity
    ``avg_depth`` : *dict*
        A dictionary with average depth of each detected cavity
    ``avg_hydropathy`` : *dict*
        A dictionary with average hydropathy of each detected cavity and range of the hydrophobicity scale mapped
    ``residues`` : *dict*
        A dictionary with interface residues of each detected cavity
    ``frequencies`` : *dict*
        A dictionary with frequencies of interface residues and classes of residues of each detected cavity
    ``step`` : *float, default 0.6*
        Grid spacing (A)

:Returns:

    File with TOML-formatted data corresponding to file paths and cavity characterization per detected cavity.

.. note::

    The cavity nomenclature is based on the integer label. The cavity marked with 2, the first integer corresponding to a cavity, is KAA, the cavity marked with 3 is KAB, the cavity marked with 4 is KAC and so on. 

.. seealso::

    * `pyKVFinder.detect <detect.html>`_
    * `pyKVFinder.spatial <spatial.html>`_
    * `pyKVFinder.depth <depth.html>`_
    * `pyKVFinder.constitutional <constitutional.html>`_
    * `pyKVFinder.hydropathy <hydropathy.html>`_
    * `pyKVFinder.export <export.html>`_

.. raw:: html

  <h4><u>Example</u></h4>

With the cavity and surface points identified and depth and hydrophobicity scale mapped in the 3D grid, we can:

* Write cavity detection:

.. code-block:: python

  >>> from pyKVFinder import write_results
  >>> import os
  >>> pdb = os.path.join(os.path.dirname(pyKVFinder.__file__), 'data', 'tests', '1FMO.pdb')
  >>> write_results('results.toml', pdb=pdb, ligand=None, output='cavities.pdb')

* Write spatial characterization

.. code-block:: python

    >>> volume
    {'KAA': 137.16, 'KAB': 47.52, 'KAC': 66.96, 'KAD': 8.21, 'KAE': 43.63, 'KAF': 12.53, 'KAG': 6.26, 'KAH': 520.13, 'KAI': 12.31, 'KAJ': 26.57, 'KAK': 12.31, 'KAL': 33.91, 'KAM': 23.11, 'KAN': 102.82, 'KAO': 6.05, 'KAP': 15.55, 'KAQ': 7.99, 'KAR': 7.78}
    >>> area
    {'KAA': 120.52, 'KAB': 58.76, 'KAC': 72.06, 'KAD': 17.62, 'KAE': 56.44, 'KAF': 22.53, 'KAG': 15.38, 'KAH': 489.25, 'KAI': 29.87, 'KAJ': 44.85, 'KAK': 30.58, 'KAL': 43.59, 'KAM': 45.25, 'KAN': 129.77, 'KAO': 11.57, 'KAP': 24.8, 'KAQ': 12.59, 'KAR': 15.97}
    >>> write_results('results.toml', pdb=pdb, ligand=None, output=None, volume=volume, area=area)

* Write constitutional characterization

.. code-block:: python

    >>> residues
    {'KAA': [['14', 'E', 'SER'], ['15', 'E', 'VAL'], ['18', 'E', 'PHE'], ['19', 'E', 'LEU'], ['100', 'E', 'PHE'], ['152', 'E', 'LEU'], ['155', 'E', 'GLU'], ['156', 'E', 'TYR'], ['292', 'E', 'LYS'], ['302', 'E', 'TRP'], ['303', 'E', 'ILE'], ['306', 'E', 'TYR']], 'KAB': [['18', 'E', 'PHE'], ['22', 'E', 'ALA'], ['25', 'E', 'ASP'], ['26', 'E', 'PHE'], ['29', 'E', 'LYS'], ['97', 'E', 'ALA'], ['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['156', 'E', 'TYR']], 'KAC': [['141', 'E', 'PRO'], ['142', 'E', 'HIS'], ['144', 'E', 'ARG'], ['145', 'E', 'PHE'], ['148', 'E', 'ALA'], ['299', 'E', 'THR'], ['300', 'E', 'THR'], ['305', 'E', 'ILE'], ['310', 'E', 'VAL'], ['311', 'E', 'GLU'], ['313', 'E', 'PRO']], 'KAD': [['122', 'E', 'TYR'], ['124', 'E', 'ALA'], ['176', 'E', 'GLN'], ['318', 'E', 'PHE'], ['320', 'E', 'GLY'], ['321', 'E', 'PRO'], ['322', 'E', 'GLY'], ['323', 'E', 'ASP']], 'KAE': [['95', 'E', 'LEU'], ['98', 'E', 'VAL'], ['99', 'E', 'ASN'], ['100', 'E', 'PHE'], ['103', 'E', 'LEU'], ['104', 'E', 'VAL'], ['105', 'E', 'LYS'], ['106', 'E', 'LEU']], 'KAF': [['123', 'E', 'VAL'], ['124', 'E', 'ALA'], ['175', 'E', 'ASP'], ['176', 'E', 'GLN'], ['181', 'E', 'GLN']], 'KAG': [['34', 'E', 'SER'], ['37', 'E', 'THR'], ['96', 'E', 'GLN'], ['106', 'E', 'LEU'], ['107', 'E', 'GLU'], ['108', 'E', 'PHE'], ['109', 'E', 'SER']], 'KAH': [['49', 'E', 'LEU'], ['50', 'E', 'GLY'], ['51', 'E', 'THR'], ['52', 'E', 'GLY'], ['53', 'E', 'SER'], ['54', 'E', 'PHE'], ['55', 'E', 'GLY'], ['56', 'E', 'ARG'], ['57', 'E', 'VAL'], ['70', 'E', 'ALA'], ['72', 'E', 'LYS'], ['74', 'E', 'LEU'], ['84', 'E', 'GLN'], ['87', 'E', 'HIS'], ['88', 'E', 'THR'], ['91', 'E', 'GLU'], ['104', 'E', 'VAL'], ['120', 'E', 'MET'], ['121', 'E', 'GLU'], ['122', 'E', 'TYR'], ['123', 'E', 'VAL'], ['127', 'E', 'GLU'], ['166', 'E', 'ASP'], ['168', 'E', 'LYS'], ['170', 'E', 'GLU'], ['171', 'E', 'ASN'], ['173', 'E', 'LEU'], ['183', 'E', 'THR'], ['184', 'E', 'ASP'], ['186', 'E', 'GLY'], ['187', 'E', 'PHE'], ['201', 'E', 'THR'], ['327', 'E', 'PHE']], 'KAI': [['131', 'E', 'HIS'], ['138', 'E', 'PHE'], ['142', 'E', 'HIS'], ['146', 'E', 'TYR'], ['174', 'E', 'ILE'], ['314', 'E', 'PHE']], 'KAJ': [['33', 'E', 'PRO'], ['89', 'E', 'LEU'], ['92', 'E', 'LYS'], ['93', 'E', 'ARG'], ['96', 'E', 'GLN'], ['349', 'E', 'GLU'], ['350', 'E', 'PHE']], 'KAK': [['157', 'E', 'LEU'], ['162', 'E', 'LEU'], ['163', 'E', 'ILE'], ['164', 'E', 'TYR'], ['185', 'E', 'PHE'], ['188', 'E', 'ALA']], 'KAL': [['49', 'E', 'LEU'], ['127', 'E', 'GLU'], ['129', 'E', 'PHE'], ['130', 'E', 'SER'], ['326', 'E', 'ASN'], ['327', 'E', 'PHE'], ['328', 'E', 'ASP'], ['330', 'E', 'TYR']], 'KAM': [['51', 'E', 'THR'], ['55', 'E', 'GLY'], ['56', 'E', 'ARG'], ['73', 'E', 'ILE'], ['74', 'E', 'LEU'], ['75', 'E', 'ASP'], ['115', 'E', 'ASN'], ['335', 'E', 'ILE'], ['336', 'E', 'ARG']], 'KAN': [['165', 'E', 'ARG'], ['166', 'E', 'ASP'], ['167', 'E', 'LEU'], ['199', 'E', 'CYS'], ['200', 'E', 'GLY'], ['201', 'E', 'THR'], ['204', 'E', 'TYR'], ['205', 'E', 'LEU'], ['206', 'E', 'ALA'], ['209', 'E', 'ILE'], ['219', 'E', 'VAL'], ['220', 'E', 'ASP'], ['223', 'E', 'ALA']], 'KAO': [['48', 'E', 'THR'], ['51', 'E', 'THR'], ['56', 'E', 'ARG'], ['330', 'E', 'TYR'], ['331', 'E', 'GLU']], 'KAP': [['222', 'E', 'TRP'], ['238', 'E', 'PHE'], ['253', 'E', 'GLY'], ['254', 'E', 'LYS'], ['255', 'E', 'VAL'], ['273', 'E', 'LEU']], 'KAQ': [['207', 'E', 'PRO'], ['208', 'E', 'GLU'], ['211', 'E', 'LEU'], ['213', 'E', 'LYS'], ['275', 'E', 'VAL'], ['277', 'E', 'LEU']], 'KAR': [['237', 'E', 'PRO'], ['238', 'E', 'PHE'], ['249', 'E', 'LYS'], ['254', 'E', 'LYS'], ['255', 'E', 'VAL'], ['256', 'E', 'ARG']]}
    >>> frequencies
    {'KAA': {'RESIDUES': {'GLU': 1, 'ILE': 1, 'LEU': 2, 'LYS': 1, 'PHE': 2, 'SER': 1, 'TRP': 1, 'TYR': 2, 'VAL': 1}, 'CLASS': {'R1': 4, 'R2': 5, 'R3': 1, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAB': {'RESIDUES': {'ALA': 2, 'ASN': 1, 'ASP': 1, 'LYS': 1, 'PHE': 2, 'TYR': 1, 'VAL': 1}, 'CLASS': {'R1': 3, 'R2': 3, 'R3': 1, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAC': {'RESIDUES': {'ALA': 1, 'ARG': 1, 'GLU': 1, 'HIS': 1, 'ILE': 1, 'PHE': 1, 'PRO': 2, 'THR': 2, 'VAL': 1}, 'CLASS': {'R1': 5, 'R2': 1, 'R3': 2, 'R4': 1, 'R5': 2, 'RX': 0}}, 'KAD': {'RESIDUES': {'ALA': 1, 'ASP': 1, 'GLN': 1, 'GLY': 2, 'PHE': 1, 'PRO': 1, 'TYR': 1}, 'CLASS': {'R1': 4, 'R2': 2, 'R3': 1, 'R4': 1, 'R5': 0, 'RX': 0}}, 'KAE': {'RESIDUES': {'ASN': 1, 'LEU': 3, 'LYS': 1, 'PHE': 1, 'VAL': 2}, 'CLASS': {'R1': 5, 'R2': 1, 'R3': 1, 'R4': 0, 'R5': 1, 'RX': 0}}, 'KAF': {'RESIDUES': {'ALA': 1, 'ASP': 1, 'GLN': 2, 'VAL': 1}, 'CLASS': {'R1': 2, 'R2': 0, 'R3': 2, 'R4': 1, 'R5': 0, 'RX': 0}}, 'KAG': {'RESIDUES': {'GLN': 1, 'GLU': 1, 'LEU': 1, 'PHE': 1, 'SER': 2, 'THR': 1}, 'CLASS': {'R1': 1, 'R2': 1, 'R3': 4, 'R4': 1, 'R5': 0, 'RX': 0}}, 'KAH': {'RESIDUES': {'ALA': 1, 'ARG': 1, 'ASN': 1, 'ASP': 2, 'GLN': 1, 'GLU': 4, 'GLY': 4, 'HIS': 1, 'LEU': 3, 'LYS': 2, 'MET': 1, 'PHE': 3, 'SER': 1, 'THR': 4, 'TYR': 1, 'VAL': 3}, 'CLASS': {'R1': 11, 'R2': 4, 'R3': 8, 'R4': 6, 'R5': 4, 'RX': 0}}, 'KAI': {'RESIDUES': {'HIS': 2, 'ILE': 1, 'PHE': 2, 'TYR': 1}, 'CLASS': {'R1': 1, 'R2': 3, 'R3': 0, 'R4': 0, 'R5': 2, 'RX': 0}}, 'KAJ': {'RESIDUES': {'ARG': 1, 'GLN': 1, 'GLU': 1, 'LEU': 1, 'LYS': 1, 'PHE': 1, 'PRO': 1}, 'CLASS': {'R1': 2, 'R2': 1, 'R3': 1, 'R4': 1, 'R5': 2, 'RX': 0}}, 'KAK': {'RESIDUES': {'ALA': 1, 'ILE': 1, 'LEU': 2, 'PHE': 1, 'TYR': 1}, 'CLASS': {'R1': 4, 'R2': 2, 'R3': 0, 'R4': 0, 'R5': 0, 'RX': 0}}, 'KAL': {'RESIDUES': {'ASN': 1, 'ASP': 1, 'GLU': 1, 'LEU': 1, 'PHE': 2, 'SER': 1, 'TYR': 1}, 'CLASS': {'R1': 1, 'R2': 3, 'R3': 2, 'R4': 2, 'R5': 0, 'RX': 0}}, 'KAM': {'RESIDUES': {'ARG': 2, 'ASN': 1, 'ASP': 1, 'GLY': 1, 'ILE': 2, 'LEU': 1, 'THR': 1}, 'CLASS': {'R1': 4, 'R2': 0, 'R3': 2, 'R4': 1, 'R5': 2, 'RX': 0}}, 'KAN': {'RESIDUES': {'ALA': 2, 'ARG': 1, 'ASP': 2, 'CYS': 1, 'GLY': 1, 'ILE': 1, 'LEU': 2, 'THR': 1, 'TYR': 1, 'VAL': 1}, 'CLASS': {'R1': 7, 'R2': 1, 'R3': 2, 'R4': 2, 'R5': 1, 'RX': 0}}, 'KAO': {'RESIDUES': {'ARG': 1, 'GLU': 1, 'THR': 2, 'TYR': 1}, 'CLASS': {'R1': 0, 'R2': 1, 'R3': 2, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAP': {'RESIDUES': {'GLY': 1, 'LEU': 1, 'LYS': 1, 'PHE': 1, 'TRP': 1, 'VAL': 1}, 'CLASS': {'R1': 3, 'R2': 2, 'R3': 0, 'R4': 0, 'R5': 1, 'RX': 0}}, 'KAQ': {'RESIDUES': {'GLU': 1, 'LEU': 2, 'LYS': 1, 'PRO': 1, 'VAL': 1}, 'CLASS': {'R1': 4, 'R2': 0, 'R3': 0, 'R4': 1, 'R5': 1, 'RX': 0}}, 'KAR': {'RESIDUES': {'ARG': 1, 'LYS': 2, 'PHE': 1, 'PRO': 1, 'VAL': 1}, 'CLASS': {'R1': 2, 'R2': 1, 'R3': 0, 'R4': 0, 'R5': 3, 'RX': 0}}}
    >>> write_results('results.toml', pdb=pdb, ligand=None, output=output_cavity, volume=volume, area=area, residues=residues, frequencies=frequencies)

* Write depth characterization

.. code-block:: python

    >>> max_depth
    {'KAA': 3.79, 'KAB': 2.68, 'KAC': 2.62, 'KAD': 0.85, 'KAE': 3.0, 'KAF': 0.85, 'KAG': 0.6, 'KAH': 10.73, 'KAI': 2.55, 'KAJ': 2.24, 'KAK': 0.0, 'KAL': 3.0, 'KAM': 1.2, 'KAN': 0.0, 'KAO': 1.04, 'KAP': 2.08, 'KAQ': 0.85, 'KAR': 0.6}
    >>> avg_depth
    {'KAA': 1.28, 'KAB': 0.86, 'KAC': 0.67, 'KAD': 0.29, 'KAE': 0.98, 'KAF': 0.24, 'KAG': 0.1, 'KAH': 3.75, 'KAI': 1.5, 'KAJ': 0.96, 'KAK': 0.0, 'KAL': 1.0, 'KAM': 0.24, 'KAN': 0.0, 'KAO': 0.29, 'KAP': 0.7, 'KAQ': 0.22, 'KAR': 0.12}
    >>> write_results('results.toml', pdb=pdb, ligand=None, output=None, max_depth=max_depth, avg_depth=avg_depth)

* Write hydropathy characterization

.. code-block:: python

  >>> avg_hydropathy
  {'KAA': -0.73, 'KAB': -0.06, 'KAC': -0.07, 'KAD': -0.62, 'KAE': -0.81, 'KAF': -0.14, 'KAG': -0.33, 'KAH': -0.16, 'KAI': -0.4, 'KAJ': 0.62, 'KAK': -0.99, 'KAL': 0.35, 'KAM': -0.33, 'KAN': 0.18, 'KAO': 0.88, 'KAP': -0.96, 'KAQ': 0.48, 'KAR': 0.24, 'EisenbergWeiss': [-1.42, 2.6]}
  >>> write_results('results.toml', pdb=pdb, ligand=None, output=None, output_hydropathy='hydropathy.pdb', avg_hydropathy=avg_hydropathy)

* Write all

.. code-block:: python

  >>> write_results('results.toml', pdb=pdb, ligand=None, output='cavities.pdb', output_hydropathy='hydropathy.pdb', volume=volume, area=area, max_depth=max_depth, avg_depth=avg_depth, avg_hydropathy=avg_hydropathy, residues=residues, frequencies=frequencies)
