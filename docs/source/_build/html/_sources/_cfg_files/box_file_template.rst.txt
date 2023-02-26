Box configuration file template
###############################

There are three methods for defining a custom 3D grid in pyKVFinder.

The first directly defines four vertices of the 3D grid (origin, X-axis, Y-axis and Z-axis), an example is shown below:

.. code-block:: TOML

  [box]
  # px = [x, y, z]
  p1 = [0.0, 0.0, 0.0]
  p2 = [1.0, 0.0, 0.0]
  p3 = [0.0, 1.0, 0.0]
  p4 = [0.0, 0.0, 1.0]

Example: `custom-box.toml <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/tests/custom-box.toml>`_.

The second defines a list of residues and a padding, the template is shown below:

.. code-block:: TOML

  [box]
  residues = [ ["resnum", "chain", "resname",], ["resnum", "chain", "resname",], ]
  padding =  3.5

Example: `residues-box.toml <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/tests/residues-box.toml>`_.

The third uses `parKVFinder <https://github.com/LBC-LNBio/parKVFinder>`_'s TOML-formatted parameters file created by its PyMOL plugin.

.. code-block:: python

	[SETTINGS.visiblebox.p1]
	x = 0.0
	y = 0.0
	z = 0.0

	[SETTINGS.visiblebox.p2]
	x = 1.0
	y = 0.0
	z = 0.0

	[SETTINGS.visiblebox.p3]
	x = 0.0
	y = 1.0
	z = 0.0

	[SETTINGS.visiblebox.p4]
	x = 0.0
	y = 0.0
	z = 1.0
