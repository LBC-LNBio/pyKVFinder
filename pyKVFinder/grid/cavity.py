import numpy


def _select_cavities(cavities: numpy.ndarray, selection: list[int]) -> numpy.ndarray:
    """Select cavities in the 3D grid by cavity labels.

    Parameters
    ----------
    cavities : numpy.ndarray
        Cavity points in the 3D grid (cavities[nx][ny][nz]).
        Cavities array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule points;

            * 1: empty space points;

            * >=2: cavity points.
    selection : list[int]
        A list of integer labels of each cavity to be selected.

    Returns
    -------
    selected : numpy.ndarray
        Selected cavity points in the 3D grid.
        Selected cavities array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule points;

            * 1: empty space points or unselected cavity points;

            * >=2: selected cavity points.
    """
    # Copy cavities
    selected = numpy.copy(cavities)

    # When outside selection, change cavities tags to 1
    for cav in range(2, cavities.max() + 1):
        if cav not in selection:
            selected[cavities == cav] = 1

    return selected


def _get_cavity_name(index: int) -> str:
    """Get cavity name, eg KAA, KAB, and so on, based on the index.

    Naming convention:
    - 0-675   -> KAA ... KZZ
    - 676-1351 -> Kaa ... Kzz

    Parameters
    ----------
    index : int
        Index in the dictionary.

    Returns
    -------
    cavity_name : str
        Cavity name
    """
    # Get block and offset
    block = index // (26 * 26)
    offset = index % (26 * 26)

    # Choose ASCII base: uppercase or lowercase
    base = 65 if block == 0 else 97  # 'A' or 'a'

    # Get first and second characters
    first = chr(base + (offset // 26))
    second = chr(base + (offset % 26))

    # Build cavity name
    cavity_name = f"K{first}{second}"

    return cavity_name


def _get_cavity_label(cavity_name: str) -> int:
    """Get cavity label, eg 2, 3, and so on, based on the cavity name.

    Naming convention:
    - KAA-KZZ -> 2-677
    - Kaa-Kzz -> 678-1353

    Parameters
    ----------
    cavity_name : str
        Cavity name

    Returns
    -------
    cavity_label : int
        Integer label of each cavity.

    Raises
    ------
    ValueError
        Invalid cavity name: `cavity_name`.
    """
    # Check cavity name
    if cavity_name[0] != "K":
        raise ValueError(f"Invalid cavity name: {cavity_name}.")

    # Get cavity label
    c1, c2 = cavity_name[1], cavity_name[2]

    # Uppercase block: KAA–KZZ
    if c1.isupper() and c2.isupper():
        base = 65  # 'A'
        block_offset = 0

    # Lowercase block: Kaa–Kzz
    elif c1.islower() and c2.islower():
        base = 97  # 'a'
        block_offset = 26 * 26

    # +2 preserves the original labeling convention
    cavity_label = (
        block_offset
        + (ord(c1) - base) * 26
        + (ord(c2) - base)
    ) + 2

    return cavity_label
