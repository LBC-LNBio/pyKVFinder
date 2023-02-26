Van der Waals radii file template
#################################

The van der Waals radii file define the radius values for each atom by residue and when not defined, it uses a generic value based on the atom type. The package contains a built-in van der Waals radii file: `vdw.dat <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/vdw.dat>`_. However, the user can define its own file with a mandatory format and pass it to pyKVFinder. The format is shown below:

.. code-block:: none

    >RES
    C       1.66
    CA      2.00
    N       1.97
    O       1.69
    H       0.91

.. warning::
  
    The residue name should be in the standard PDB format and each radius value is separated by two tab characters of the atom name.
