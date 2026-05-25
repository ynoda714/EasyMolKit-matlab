function smiles = toSmiles(mol)
% toSmiles  Convert an RDKit Mol object to a canonical SMILES string.
%
%   smiles = emk.mol.toSmiles(mol)
%
%   Converts a Python RDKit Mol reference to its canonical SMILES string
%   using RDKit's MolToSmiles().  The canonical form is deterministic for a
%   given RDKit version: the same molecule always yields the same string.
%   Note that canonical SMILES may differ from the original input SMILES
%   (e.g., "OCC" is normalised to "CCO").
%
%   Arguments:
%     mol     - py.rdkit.Chem.rdchem.Mol  RDKit molecule object (Python ref)
%
%   Returns:
%     smiles  - string  Canonical SMILES string
%
%   Errors:
%     emk:mol:toSmiles:invalidInput  - mol is not a py.rdkit.Chem.rdchem.Mol
%     emk:mol:toSmiles:rdkitError    - unexpected Python exception
%
%   Example:
%     mol    = emk.mol.fromSmiles("OCC");
%     smiles = emk.mol.toSmiles(mol);   % returns "CCO" (canonical)
%
%   See also: emk.mol.fromSmiles, emk.mol.isValid

    % --- Input validation ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:mol:toSmiles:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    logDebug("toSmiles: converting mol with %d atoms", double(mol.GetNumAtoms()));

    % --- Call RDKit ---
    % Use emk.util.rdkitModule() to get the Chem module via importlib.
    % Direct py.rdkit.* access fails (TypeError); importlib avoids this.
    mods = emk.util.rdkitModule();
    try
        pySmiles = mods.Chem.MolToSmiles(mol);
    catch ME
        error("emk:mol:toSmiles:rdkitError", ...
            "RDKit MolToSmiles raised an exception: %s", ME.message);
    end

    smiles = string(pySmiles);

    logDebug("toSmiles: result = '%s'", smiles);
end
