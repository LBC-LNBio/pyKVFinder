Hydrophobicity scale file template
##################################

The hydrophobicity scale file defines the name of the scale and the hydrophobicity value for each residue and when not defined, it assigns zero to the missing residues. The package contains six built-in hydrophobicity scales: `Eisenberg & Weiss <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/EisenbergWeiss.toml>`_ [1], `Hessa & Heijne <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/HessaHeijne.toml>`_ [2], `Kyte & Doolittle <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/KyteDoolittle.toml>`_ [3], `Moon & Fleming <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/MoonFleming.toml>`_ [4], `Wimley & White <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/WimleyWhite.toml>`_ [5] and `Zhao & London <https://github.com/LBC-LNBio/pyKVFinder/blob/master/pyKVFinder/data/ZhaoLondon.toml>`_ [6]. However, the user can define its own file with a mandatory format and pass it to pyKVFinder. The format is shown below:

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

.. raw:: html

    <h4><u>References</u></h4>

1. Eisenberg D, Weiss RM, Terwilliger TC. The hydrophobic moment detects periodicity in protein hydrophobicity. Proceedings of the National Academy of Sciences. 1984;81. 

2. Hessa T, Kim H, Bihlmaier K, Lundin C, Boekel J, Andersson H, et al. Recognition of transmembrane helices by the endoplasmic reticulum translocon. Nature. 2005;433. 

3. Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. Journal of Molecular Biology. 1982;157. 

4. Moon CP, Fleming KG. Side-chain hydrophobicity scale derived from transmembrane protein folding into lipid bilayers. Proceedings of the National Academy of Sciences. 2011;108. 

5. Wimley WC, White SH. Experimentally determined hydrophobicity scale for proteins at membrane interfaces. Nature Structural & Molecular Biology. 1996;3. 

6. Zhao G, London E. An amino acid “transmembrane tendency” scale that approaches the theoretical limit to accuracy for prediction of transmembrane helices: Relationship to biological hydrophobicity. Protein Science. 2006;15. 
