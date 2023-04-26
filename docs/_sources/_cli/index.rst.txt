.. raw:: html

    <h1>Command-line interface</h1>

In addition to pyKVFinder package, a command-line interface (CLI) is available to ease cavity detection and characterization with the full set of customizable parameters.

.. code-block:: bash

  $ pyKVFinder
  Usage: pyKVFinder [-h] [-v] [--version] [-b <str>] [-O <str>]
                    [--nthreads <int>] [-d <str>] [-s <float>] [-i <float>]
                    [-o <float>] [-V <float>] [-R <float>] [-S <str>]
                    [--ignore_backbone] [-D] [--plot_frequencies]
                    [--hydropathy [{EisenbergWeiss, HessaHeijne, KyteDoolittle, 
                    MoonFleming, RadzickaWolfenden, WimleyWhite, ZhaoLondon, <.toml>}]]
                    [-B <.toml>] [-L (<.pdb> | <.xyz>)] [--ligand_cutoff <float>]
                    (<.pdb> | <.xyz>)

Positional arguments
--------------------

The positional arguments are:

* ``(<.pdb> | <.xyz>)``: A path to a target PDB or XYZ file.
  
  .. code-block:: bash
    
    $ pyKVFinder <.pdb>

Optional arguments
------------------

The optional arguments are:

* ``-h`` or ``--help``: Show help message.
  
  .. code-block:: bash
    
    $ pyKVFinder -h
    $ pyKVFinder --help

* ``--version``: Display pyKVFinder version.
  
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

Detection and characterization
******************************

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

Extra characterization
**********************

The parameters for additional characterization are:

* ``--D or --depth``: Characterize the depth of the detected cavities. This mode includes depth of each cavity point as the B-factor in the cavity PDB file and maximum and average depth of the detected cavities in the results file.

  .. code-block:: bash

    $ pyKVFinder <.pdb> -D
    $ pyKVFinder <.pdb> --depth

  :Default: ``None``

* ``--plot_frequencies``: Plot bar charts of calculated frequencies (residues and classes of residues) of the detected cavities in a PDF file. The classes of residues are aliphatic apolar (R1), aromatic (R2), polar uncharged (R3), negatively charged (R4), positively charged (R5) and non-standard (RX) residues.

  .. code-block:: bash

    $ pyKVFinder <.pdb> --plot_frequencies

  :Default: ``None``

* ``--hydropathy [{EisenbergWeiss, HessaHeijne, KyteDoolittle, MoonFleming, WimleyWhite, ZhaoLondon, <.toml>}]``: Characterize the hydropathy of the detected cavities. This mode maps a target hydrophobicity scale as B-factor at surface points of the detected cavities. Also, it calculates the average hydropathy of each detected cavity. The constant hydrophobicity scale is EisenbergWeiss.

  .. code-block:: bash

    $ pyKVFinder <.pdb> --hydropathy

  In addition, the user can define one of the built-in hydrophobicity scale. The built-in hydrophobicity scales are: EisenbergWeiss, HessaHeijne, KyteDoolittle, MoonFleming, RadzickaWolfenden, WimleyWhite and ZhaoLondon.

  .. code-block:: bash

    $ pyKVFinder <.pdb> --hydropathy EisenbergWeiss
    $ pyKVFinder <.pdb> --hydropathy HessaHeijne
    $ pyKVFinder <.pdb> --hydropathy KyteDoolittle
    $ pyKVFinder <.pdb> --hydropathy MoonFleming
    $ pyKVFinder <.pdb> --hydropathy RadzickaWolfenden
    $ pyKVFinder <.pdb> --hydropathy WimleyWhite
    $ pyKVFinder <.pdb> --hydropathy ZhaoLondon

  Further, the user can also define a custom hydrophobicity scale file via a TOML-formatted file (see template).

  .. code-block:: bash

    $ pyKVFinder <.pdb> --hydropathy <.toml>

  :Default: ``None``
  :Constant: ``EisenbergWeiss``

Box adjusment
*************

The box adjustment argument is:

* ``-B <.toml>`` or ``--box <.toml>``: A path to TOML-formatted file with box parameters (see template). Adjust the 3D grid based on a list of residues (["resnum", "chain"]) and a padding or a set of four vertices (p1: origin, p2: X-axis max, p3: Y-axis max, p4: Z-axis max) with xyz coordinates ([x, y, z]).

  .. code-block:: bash

    $ pyKVFinder <.pdb> -B <.toml>
    $ pyKVFinder <.pdb> --box <.toml>

  :Default: ``None``

Ligand adjustment
*****************

The ligand adjustment arguments are:

* ``-L (<.pdb> | <.xyz>)`` or ``--ligand (<.pdb> | <.xyz>)``: A path to a ligand PDB or XYZ file to limit the cavities within a radius (ligand_cutoff) around it.

  .. code-block:: bash

    $ pyKVFinder <.pdb> -L <.pdb>
    $ pyKVFinder <.pdb> --ligand <.pdb>

  :Default: ``None``

* ``--ligand_cutoff <float>``: A radius value to limit a space around the defined ligand.

  .. code-block:: bash

    $ pyKVFinder <.pdb> -L <.pdb> --ligand_cutoff <float>

  :Default: ``5.0``
