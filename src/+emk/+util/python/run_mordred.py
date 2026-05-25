"""
run_mordred.py

Mordred molecular descriptor calculation helper for EasyMolKit.
Called from MATLAB via:
  m = py.importlib.import_module("run_mordred");
  m.mordred_list_names()
  m.mordred_calculate(mol, names)
  m.mordred_batch(mols, names)

Reduces IPC round-trips by computing all descriptors inside a single
Python call rather than looping from MATLAB.

ADR-002 rev.3: IPC minimisation principle for batch APIs.
"""

# ---------------------------------------------------------------------------
# Module-level cache: reuse Calculator across calls within the same session
# ---------------------------------------------------------------------------
_CALC_ALL    = None   # Calculator(all 2D descriptors, ignore_3D=True)
_ALL_NAMES   = None   # sorted list of all descriptor name strings


def _get_calc_all():
    """Return (and cache) a Calculator for all 2D descriptors."""
    global _CALC_ALL, _ALL_NAMES
    if _CALC_ALL is None:
        from mordred import Calculator, descriptors as _all_descs
        _CALC_ALL  = Calculator(_all_descs, ignore_3D=True)
        _ALL_NAMES = sorted(str(d) for d in _CALC_ALL.descriptors)
    return _CALC_ALL


def _build_filtered_calc(names_list):
    """Return a Calculator containing only the named descriptors."""
    calc_all = _get_calc_all()
    # Build name->descriptor mapping once
    desc_map = {str(d): d for d in calc_all.descriptors}
    selected = [desc_map[n] for n in names_list if n in desc_map]
    from mordred import Calculator
    return Calculator(selected, ignore_3D=True)


def _result_to_dict(result, calc):
    """Convert a Mordred Result to a plain dict of {name: float|nan}."""
    _NAN = float("nan")
    out = {}
    for d in calc.descriptors:
        key = str(d)
        val = result[d]
        try:
            out[key] = float(val)
        except (TypeError, ValueError):
            out[key] = _NAN
    return out


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def mordred_list_names():
    """Return sorted list of all 2D Mordred descriptor name strings.

    Returns
    -------
    list of str
        Names can be passed to mordred_calculate / mordred_batch as the
        ``names`` argument to request a subset of descriptors.
    """
    _get_calc_all()
    return list(_ALL_NAMES)


def mordred_calculate(mol, names):
    """Compute Mordred descriptors for a single molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        RDKit Mol object.  None yields an empty dict.
    names : iterable of str or empty list
        Descriptor names to compute.  Empty list means all 2D descriptors.

    Returns
    -------
    dict mapping descriptor name (str) -> float.
        Failed / non-numeric values are represented as float('nan').
    """
    if mol is None:
        return {}

    names_list = list(names)
    if len(names_list) == 0:
        calc = _get_calc_all()
    else:
        calc = _build_filtered_calc(names_list)

    try:
        result = calc(mol)
    except Exception:
        return {str(d): float("nan") for d in calc.descriptors}

    return _result_to_dict(result, calc)


def mordred_batch(mols, names):
    """Compute Mordred descriptors for multiple molecules.

    Parameters
    ----------
    mols : iterable of rdkit.Chem.rdchem.Mol or None
        Molecule objects to process.  None entries yield NaN rows.
    names : iterable of str or empty list
        Descriptor names to compute.  Empty list means all 2D descriptors.

    Returns
    -------
    list of dict
        Each element is {descriptor_name: float}.  len == len(mols).
        Failed descriptors are float('nan').
    """
    mols_list  = list(mols)
    names_list = list(names)

    if len(names_list) == 0:
        calc = _get_calc_all()
    else:
        calc = _build_filtered_calc(names_list)

    results = []
    for mol in mols_list:
        if mol is None:
            results.append({str(d): float("nan") for d in calc.descriptors})
            continue
        try:
            result = calc(mol)
            results.append(_result_to_dict(result, calc))
        except Exception:
            results.append({str(d): float("nan") for d in calc.descriptors})

    return results


def mordred_calculate_list(mol, names):
    """Compute Mordred descriptors and return as (names_list, values_list).

    Designed for efficient MATLAB integration: a single Python call returns
    both column names and values, eliminating N IPC round-trips per molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol or None
        RDKit Mol object.  None yields empty lists.
    names : iterable of str or empty list
        Descriptor names to compute.  Empty list means all 2D descriptors.

    Returns
    -------
    tuple of (list of str, list of float)
        names_list  - descriptor names (same order as values_list)
        values_list - corresponding double values (NaN for failures)
    """
    names_list = list(names)
    if len(names_list) == 0:
        _get_calc_all()
        names_list = list(_ALL_NAMES)

    if mol is None:
        return (names_list, [float("nan")] * len(names_list))

    d = mordred_calculate(mol, names_list)
    values = [d.get(n, float("nan")) for n in names_list]
    return (names_list, values)


def mordred_batch_matrix(mols, names):
    """Compute Mordred descriptors for multiple molecules as a matrix.

    Designed for efficient MATLAB integration: a single Python call returns
    the full N x M matrix, eliminating N x M IPC round-trips.

    Parameters
    ----------
    mols : iterable of rdkit.Chem.rdchem.Mol or None
        Molecule objects.  None entries yield NaN rows.
    names : iterable of str or empty list
        Descriptor names to compute.  Empty list means all 2D descriptors.

    Returns
    -------
    tuple of (list of str, list of list of float)
        names_list - descriptor names (length M)
        matrix     - list of N rows, each a list of M floats
    """
    names_list = list(names)
    if len(names_list) == 0:
        _get_calc_all()
        names_list = list(_ALL_NAMES)

    dicts = mordred_batch(mols, names_list)
    matrix = [[d.get(n, float("nan")) for n in names_list] for d in dicts]
    return (names_list, matrix)
