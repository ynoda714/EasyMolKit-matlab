function mol = fromSmiles(smiles)
% fromSmiles  Create an RDKit Mol object from a SMILES string.
%
%   mol = emk.mol.fromSmiles(smiles)
%
%   Parses a SMILES string using RDKit and returns a Python Mol reference.
%   The returned object is held as a py.rdkit.Chem.rdchem.Mol reference
%   (ADR-002: Python object reference retention).  Downstream callers pass
%   this reference directly to emk.descriptor.*, emk.fingerprint.*, etc.
%   without converting to a MATLAB native type.
%
%   Arguments:
%     smiles  - (string | char) SMILES string, e.g. "CCO"
%
%   Returns:
%     mol     - py.rdkit.Chem.rdchem.Mol  RDKit molecule object (Python ref)
%
%   Errors:
%     emk:mol:fromSmiles:invalidInput   - smiles is not a string/char scalar
%     emk:mol:fromSmiles:invalidSmiles  - RDKit returned None (parse failed)
%     emk:mol:fromSmiles:rdkitError     - unexpected Python exception
%
%   Example:
%     mol = emk.mol.fromSmiles("CCO");         % ethanol
%     mol = emk.mol.fromSmiles("c1ccccc1");    % benzene
%
%   See also: emk.mol.toSmiles, emk.mol.isValid, emk.descriptor.molWeight

    % --- Input validation ---
    if ~(ischar(smiles) || isStringScalar(smiles))
        error("emk:mol:fromSmiles:invalidInput", ...
            "smiles must be a string scalar, got: %s", class(smiles));
    end

    smiles = string(smiles);

    % Guard: reject empty / whitespace-only strings BEFORE calling RDKit.
    % RDKit's MolFromSmiles("") behaviour is version-dependent: some versions
    % return None (=> invalidSmiles), others return an empty Mol object with
    % 0 atoms (=> silent success that misleads callers).  An explicit check
    % ensures consistent semantics across all RDKit versions.
    if strlength(strtrim(smiles)) == 0
        error("emk:mol:fromSmiles:invalidSmiles", ...
            "Invalid SMILES: empty or whitespace-only string is not a valid molecule.");
    end

    logDebug("fromSmiles: input = '%s'", smiles);

    % --- Call RDKit ---
    % Use emk.util.rdkitModule() to get the Chem module via importlib.
    % Direct py.rdkit.Chem.* access fails in MATLAB (TypeError: rdkit object
    % is not callable); importlib resolves the module reference correctly.
    mods = emk.util.rdkitModule();
    try
        mol = mods.Chem.MolFromSmiles(smiles);
    catch ME
        error("emk:mol:fromSmiles:rdkitError", ...
            "RDKit MolFromSmiles raised an exception for '%s': %s", ...
            smiles, ME.message);
    end

    % --- None check: RDKit returns None for an unparseable SMILES ---
    % In MATLAB's Python interface, a Python None value is represented as
    % py.NoneType.  isa() is the idiomatic check.
    if isa(mol, "py.NoneType")
        error("emk:mol:fromSmiles:invalidSmiles", ...
            "Invalid SMILES: RDKit could not parse '%s'.", smiles);
    end

    logDebug("fromSmiles: success, numAtoms = %d", double(mol.GetNumAtoms()));
end
