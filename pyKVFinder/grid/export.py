import os
import pathlib
import tempfile

import numpy

from .cavity import _get_cavity_label, _select_cavities
from .geometry import _get_sincos
from .openings import _get_opening_label


def export(
    fn: str | pathlib.Path | None,
    cavities: numpy.ndarray,
    surface: numpy.ndarray | None,
    vertices: numpy.ndarray | list[list[float]],
    step: float | int = 0.6,
    B: numpy.ndarray | None = None,
    Q: numpy.ndarray | None = None,
    selection: list[int] | list[str] | None = None,
    nthreads: int | None = None,
    append: bool = False,
    model: int = 0,
) -> str | None:
    """Export cavitiy (H) and surface (HA) points to PDB-formatted file with
    a variable (B; optional) in B-factor column, and hydropathy to
    PDB-formatted file in B-factor column at surface points (HA).

    Parameters
    ----------
    fn : str | pathlib.Path | None
        A path to PDB file for writing cavities. If None, return a raw string with the PDB-formatted file.
    cavities : numpy.ndarray
        Cavity points in the 3D grid (cavities[nx][ny][nz]).
        Cavities array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule points;

            * 1: empty space points;

            * >=2: cavity points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.
    surface : numpy.ndarray, optional
        Surface points in the 3D grid (surface[nx][ny][nz]). If None, surface
        is a numpy.zeros array with same shape of cavities.
        Surface array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: biomolecule or empty space points;

            * >=2: surface points.

        The empty space points are regions that do not meet the chosen
        volume cutoff to be considered a cavity.
    vertices : numpy.ndarray | list[list[float]]
        A numpy.ndarray or a list with xyz vertices coordinates (origin,
        X-axis, Y-axis, Z-axis).
    step : float | int, optional
        Grid spacing (A), by default 0.6.
    B : numpy.ndarray | None, optional
        A numpy.ndarray with values to be mapped on B-factor column in cavity
        points (B[nx][ny][nz]), by default None.
    Q : numpy.ndarray | None, optional
        A numpy.ndarray with hydrophobicity scale values to be mapped on
        B-factor column in surface points (Q[nx][ny][nz]), by default
        None.
    selection : list[int] | list[str] | None, optional
        A list of integer labels or a list of cavity names to be selected, by default None.
    nthreads : int | None, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    append : bool, optional
        Whether to append cavities to the PDB file, by default False.
    model : int, optional
        Model number, by default 0.

    Returns
    -------
    str | None
        A raw string with the PDB-formatted file.

    Raises
    ------
    TypeError
        "`fn` must be a string, a pathlib.Path or None. If None, return a raw string with the PDB-formatted file."
    TypeError
        `cavities` must be a numpy.ndarray.
    ValueError
        `cavities` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `surface` must be a numpy.ndarray.
    ValueError
        `surface` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `vertices` must be a list or a numpy.ndarray.
    ValueError
        `vertices` has incorrect shape. It must be (4, 3).
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.
    TypeError
        `B` must be a numpy.ndarray.
    ValueError
        `B` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `Q` must be a numpy.ndarray.
    ValueError
        `Q` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `selection` must be a list of strings (cavity names) or integers (cavity labels).
    ValueError
        Invalid `selection`: {selection}.
    TypeError
        `nthreads` must be a positive integer.
    ValueError
        `nthreads` must be a positive integer.
    TypeError
        `append` must be a boolean.
    TypeError
        `model` must be a integer.
    RuntimeError
        User must define `surface` when not defining `cavities`.

    Note
    ----
    The cavity nomenclature is based on the integer label. The cavity marked
    with 2, the first integer corresponding to a cavity, is KAA, the cavity
    marked with 3 is KAB, the cavity marked with 4 is KAC and so on. After KZZ,
    the next cavity is Kaa, Kab, and so on. The naming convention supports up to
    1352 cavities.

    See Also
    --------
    detect
    spatial
    depth
    hydropathy
    write_results

    Example
    -------
    With the cavity and surface points identified and depth and hydrophobicity scale mapped in the 3D grid, we can:

    * Export cavity points

    >>> from pyKVFinder import export
    >>> export('cavity_wo_surface.pdb', cavities, None, vertices)

    * Export cavity and surface points

    >>> export('cavities.pdb', cavities, surface, vertices)

    * Export cavity and surface points with depth mapped on them

    >>> export('cavities_with_depth.pdb', cavities, surface, vertices, B=depths)

    * Export surface points with hydrophobicity_scale mapped on them

    >>> export('cavities_with_hydropathy.pdb', cavities, surface, vertices, Q=scales)

    * Export all

    >>> export('cavities.pdb', cavities, surface, vertices, B=depths, Q=scales)

    * Export to a raw string

    >>> string = export(None, cavities, surface, vertices, B=depths, Q=scales)
    >>> print(string)
    ATOM     1  H   KAA     0       0.000   0.000   0.000  1.00  0.00
    ...
    ATOM  1000  H   KAA     0       0.000   0.000   0.000  1.00  0.00
    """
    from pyKVFinder._pyKVFinder import _export

    # Check arguments
    if not isinstance(fn, (str, pathlib.Path, type(None))):
        raise TypeError(
            "`fn` must be a string, a pathlib.Path or None. If None, return a raw string with the PDB-formatted file."
        )
    if isinstance(fn, (str, pathlib.Path)):
        os.makedirs(os.path.abspath(os.path.dirname(fn)), exist_ok=True)
    if type(cavities) not in [numpy.ndarray]:
        raise TypeError("`cavities` must be a numpy.ndarray.")
    elif len(cavities.shape) != 3:
        raise ValueError("`cavities` has the incorrect shape. It must be (nx, ny, nz).")
    if surface is not None:
        if type(surface) not in [numpy.ndarray]:
            raise TypeError("`surface` must be a numpy.ndarray.")
        elif len(surface.shape) != 3:
            raise ValueError(
                "`surface` has the incorrect shape. It must be (nx, ny, nz)."
            )
    if type(vertices) not in [numpy.ndarray, list]:
        raise TypeError("`vertices` must be a list or a numpy.ndarray.")
    elif numpy.asarray(vertices).shape != (4, 3):
        raise ValueError("`vertices` has incorrect shape. It must be (4, 3).")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")
    if B is not None:
        if type(B) not in [numpy.ndarray]:
            raise TypeError("`B` must be a numpy.ndarray.")
        elif len(B.shape) != 3:
            raise ValueError("`B` has the incorrect shape. It must be (nx, ny, nz).")
    if Q is not None:
        if type(Q) not in [numpy.ndarray]:
            raise TypeError("`Q` must be a numpy.ndarray.")
        elif len(Q.shape) != 3:
            raise ValueError("`Q` has the incorrect shape. It must be (nx, ny, nz).")
    if selection is not None:
        # Check selection types
        if all(isinstance(x, int) for x in selection):
            pass
        elif all(isinstance(x, str) for x in selection):
            selection = [_get_cavity_label(sele) for sele in selection]
        else:
            raise TypeError(
                "`selection` must be a list of strings (cavity names) or integers (cavity labels)."
            )
        # Check if selection includes valid cavity labels
        if any(x < 2 for x in selection):
            raise ValueError(f"Invalid `selection`: {selection}.")
    if nthreads is None:
        nthreads = os.cpu_count() - 1
    else:
        if type(nthreads) not in [int]:
            raise TypeError("`nthreads` must be a positive integer.")
        elif nthreads <= 0:
            raise ValueError("`nthreads` must be a positive integer.")
    if type(append) not in [bool]:
        raise TypeError("`append` must be a boolean.")
    if type(model) not in [int]:
        raise TypeError("`model` must be a integer.")

    # Convert types
    if isinstance(vertices, list):
        vertices = numpy.asarray(vertices)
    if isinstance(step, int):
        step = float(step)

    # Convert numpy.ndarray data types
    vertices = vertices.astype("float64") if vertices.dtype != "float64" else vertices
    if B is not None:
        B = B.astype("float64") if B.dtype != "float64" else B
    if Q is not None:
        Q = Q.astype("float64") if Q.dtype != "float64" else Q

    # Get sincos: sine and cosine of the grid rotation angles (sina, cosa, sinb, cosb)
    sincos = _get_sincos(vertices)

    # Unpack vertices
    P1, _, _, _ = vertices

    # If surface is None, create a zeros grid
    if surface is None:
        surface = numpy.zeros(cavities.shape, dtype="int32")

    # If B is None, create a zeros grid
    if B is None:
        B = numpy.zeros(cavities.shape, dtype="float64")

    # If Q is None, create an ones grid
    if Q is None:
        Q = numpy.ones(cavities.shape, dtype="float64")

    # Select cavities
    if selection is not None:
        cavities = _select_cavities(cavities, selection)
        surface = _select_cavities(surface, selection)

    # Get number of cavities
    ncav = int(cavities.max() - 1)

    # Export cavities
    if isinstance(fn, type(None)):
        with tempfile.NamedTemporaryFile(delete=False) as temp:
            _export(
                temp.name,
                cavities,
                surface,
                B,
                Q,
                P1,
                sincos,
                step,
                ncav,
                nthreads,
                append,
                model,
            )
            # move the file pointer to the beginning of the file
            temp.seek(0)
            # get string from temp file
            string = temp.read().decode("utf-8")
            # close file
            temp.close()
            # delete file
            if os.path.exists(temp.name):
                os.remove(temp.name)

            return string
    else:
        _export(
            fn, cavities, surface, B, Q, P1, sincos, step, ncav, nthreads, append, model
        )


def export_openings(
    fn: str | pathlib.Path | None,
    openings: numpy.ndarray,
    vertices: numpy.ndarray | list[list[float]],
    step: float | int = 0.6,
    selection: list[int] | list[str] | None = None,
    nthreads: int | None = None,
    append: bool = False,
    model: int = 0,
) -> None:
    """Export opening points (H) to a PDB-formatted file.

    Parameters
    ----------
    fn : str | pathlib.Path | None
        A path to PDB file for writing openings.
    openings : numpy.ndarray
        Openings points in the 3D grid (openings[nx][ny][nz]).
        Openings array has integer labels in each position, that are:

            * -1: bulk points;

            * 0: cavity or biomolecule points;

            * 1: empty space points;

            * >=2: Opening points.

        The empty space points are regions that do not meet the chosen
        openings cutoff to be considered an opening.
    vertices : numpy.ndarray | list[list[float]]
        A numpy.ndarray or a list with xyz vertices coordinates (origin,
        X-axis, Y-axis, Z-axis).
    step : float | int, optional
        Grid spacing (A), by default 0.6.
    selection : list[int] | list[str] | None, optional
        A list of integer labels or a list of opening names to be selected, by default None.
    nthreads : int | None, optional
        Number of threads, by default None. If None, the number of threads is
        `os.cpu_count() - 1`.
    append : bool, optional
        Whether to append openings to the PDB file, by default False.
    model : int, optional
        Model number, by default 0.

    Raises
    ------
    TypeError
        `openings` must be a numpy.ndarray.
    ValueError
        `openings` has the incorrect shape. It must be (nx, ny, nz).
    TypeError
        `vertices` must be a list or a numpy.ndarray.
    ValueError
        `vertices` has incorrect shape. It must be (4, 3).
    TypeError
        `step` must be a positive real number.
    ValueError
        `step` must be a positive real number.
    TypeError
        `selection` must be a list of strings (opening names) or integers (opening labels).
    ValueError
        Invalid `selection`: {selection}.
    TypeError
        `nthreads` must be a positive integer.
    ValueError
        `nthreads` must be a positive integer.
    TypeError
        `append` must be a boolean.
    TypeError
        `model` must be a integer.
    TypeError
        `fn` must be a string or pathlib.Path.

    Note
    ----
    The opening nomenclature is based on the integer label. The opening marked
    with 2, the first integer corresponding to a opening, is OAA, the opening
    marked with 3 is OAB, the opening marked with 4 is OAC and so on.

    See Also
    --------
    export
    detect
    depths
    openings

    Example
    -------
    With the opening points identified with ``openings``, we can export them to a PDB-formatted file:

    >>> from pyKVFinder import export_openings
    >>> export_openings('openings.pdb', openings, vertices)
    """
    from pyKVFinder._pyKVFinder import _export_openings

    # Check arguments
    if type(openings) not in [numpy.ndarray]:
        raise TypeError("`openings` must be a numpy.ndarray.")
    elif len(openings.shape) != 3:
        raise ValueError("`openings` has the incorrect shape. It must be (nx, ny, nz).")
    if type(vertices) not in [numpy.ndarray, list]:
        raise TypeError("`vertices` must be a list or a numpy.ndarray.")
    elif numpy.asarray(vertices).shape != (4, 3):
        raise ValueError("`vertices` has incorrect shape. It must be (4, 3).")
    if type(step) not in [float, int]:
        raise TypeError("`step` must be a positive real number.")
    elif step <= 0.0:
        raise ValueError("`step` must be a positive real number.")
    if selection is not None:
        # Check selection types
        if all(isinstance(x, int) for x in selection):
            pass
        elif all(isinstance(x, str) for x in selection):
            selection = [_get_opening_label(sele) for sele in selection]
        else:
            raise TypeError(
                "`selection` must be a list of strings (cavity names) or integers (cavity labels)."
            )
        # Check if selection includes valid cavity labels
        if any(x < 2 for x in selection):
            raise ValueError(f"Invalid `selection`: {selection}.")
    if nthreads is None:
        nthreads = os.cpu_count() - 1
    else:
        if type(nthreads) not in [int]:
            raise TypeError("`nthreads` must be a positive integer.")
        elif nthreads <= 0:
            raise ValueError("`nthreads` must be a positive integer.")
    if type(append) not in [bool]:
        raise TypeError("`append` must be a boolean.")
    if type(model) not in [int]:
        raise TypeError("`model` must be a integer.")

    # Convert types
    if isinstance(vertices, list):
        vertices = numpy.asarray(vertices)
    if isinstance(step, int):
        step = float(step)

    # Convert numpy.ndarray data types
    vertices = vertices.astype("float64") if vertices.dtype != "float64" else vertices

    # Get sincos: sine and cossine of the grid rotation angles (sina, cosa, sinb, cosb)
    sincos = _get_sincos(vertices)

    # Create base directories of results
    if fn is not None:
        if type(fn) not in [str, pathlib.Path]:
            raise TypeError("`fn` must be a string or a pathlib.Path.")
        os.makedirs(os.path.abspath(os.path.dirname(fn)), exist_ok=True)

    # Unpack vertices
    P1, _, _, _ = vertices

    # Select cavities
    if selection is not None:
        openings = _select_cavities(openings, selection)

    # Get number of openings
    nopenings = int(openings.max() - 1)

    # Export openings
    _export_openings(fn, openings, P1, sincos, step, nopenings, nthreads, append, model)
