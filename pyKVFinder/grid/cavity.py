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

    Parameters
    ----------
    index : int
        Index in the dictionary.

    Returns
    -------
    cavity_name : str
        Cavity name
    """
    cavity_name = f"K{chr(65 + int(index / 26) % 26)}{chr(65 + (index % 26))}"
    return cavity_name


def _get_cavity_label(cavity_name: str) -> int:
    """Get cavity label, eg 2, 3, and so on, based on the cavity name.

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
    cavity_label = (ord(cavity_name[1]) - 65) * 26 + (ord(cavity_name[2]) - 65) + 2

    return cavity_label
