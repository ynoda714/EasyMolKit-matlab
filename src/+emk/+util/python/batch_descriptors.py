"""
batch_descriptors.py

Batch descriptor calculation helper for EasyMolKit batchCalculate.
Called from MATLAB via py.importlib.import_module("batch_descriptors").

Reduces IPC round-trips from N (one per molecule) to 1 by computing
all descriptors for all molecules inside a single Python call.

ADR-002 rev.3: IPC minimisation principle for batch APIs.
"""


def batch_calculate(mols, descriptor_names):
    """Compute descriptors for a list of Mol objects in a single Python call.

    Parameters
    ----------
    mols : iterable of rdkit.Chem.Mol
        Molecule objects to process.  None entries yield NaN rows.
    descriptor_names : iterable of str
        Names of descriptors to compute.  Must match the supported set
        defined in echem.descriptor.calculate.

    Returns
    -------
    list of list of float
        Outer length = len(mols).  Inner length = len(descriptor_names).
        Failed or missing descriptor values are represented as float('nan').
    """
    from rdkit.Chem import Descriptors, rdMolDescriptors

    _NAN = float("nan")

    # Convert iterables to lists once (MATLAB may pass py.list objects)
    mols_list  = list(mols)
    names_list = list(descriptor_names)

    # Pre-resolve descriptor functions once (not per molecule)
    desc_funcs = [_resolve_descriptor(n, Descriptors, rdMolDescriptors)
                  for n in names_list]

    results = []
    for mol in mols_list:
        if mol is None:
            results.append([_NAN] * len(names_list))
            continue
        row = []
        for fn in desc_funcs:
            if fn is None:
                row.append(_NAN)
            else:
                try:
                    val = fn(mol)
                    row.append(float(val) if val is not None else _NAN)
                except Exception:
                    row.append(_NAN)
        results.append(row)

    return results


def _resolve_descriptor(name, Descriptors, rdMolDescriptors):
    """Return the callable for a descriptor name, mirroring computeOne_() in
    echem.descriptor.calculate.  Returns None if the name is unrecognised.
    """
    _direct_map = {
        "MolWt":              lambda mol: Descriptors.MolWt(mol),
        "ExactMolWt":         lambda mol: Descriptors.ExactMolWt(mol),
        "LogP":               lambda mol: Descriptors.MolLogP(mol),
        "TPSA":               lambda mol: Descriptors.TPSA(mol),
        "NumHAcceptors":      lambda mol: rdMolDescriptors.CalcNumHBA(mol),
        "NumHDonors":         lambda mol: rdMolDescriptors.CalcNumHBD(mol),
        "NumRotatableBonds":  lambda mol: rdMolDescriptors.CalcNumRotatableBonds(mol),
        "RingCount":          lambda mol: rdMolDescriptors.CalcNumRings(mol),
        "FractionCSP3":       lambda mol: Descriptors.FractionCSP3(mol),
        "HeavyAtomCount":     lambda mol: mol.GetNumHeavyAtoms(),
    }
    return _direct_map.get(name, None)
