function fp = maccs(mol)
% maccs  Generate a MACCS keys fingerprint for a molecule.
%
%   fp = emk.fingerprint.maccs(mol)
%
%   Computes a MACCS (Molecular ACCess System) 166-key fingerprint using
%   RDKit's GenMACCSKeys.  The returned bit vector has 167 bits (index 0 is
%   unused; bits 1-166 correspond to the 166 MACCS structural keys).
%
%   Unlike emk.fingerprint.morgan, MACCS keys are predefined structural
%   patterns (SMARTS-based) rather than circular hash encodings.  This makes
%   them interpretable: each bit position has a documented chemical meaning.
%   MACCS keys are most useful for scaffold-level similarity and
%   pharmacophore-feature comparisons.
%
%   The fingerprint is returned as a Python ExplicitBitVect reference
%   (ADR-002: Python object reference retention).
%   Pass to emk.similarity.tanimoto() for similarity scoring, or convert
%   to a MATLAB logical array via emk.fingerprint.toArray().
%
%   Arguments:
%     mol  - py.rdkit.Chem.rdchem.Mol  RDKit molecule object (Python ref)
%
%   Returns:
%     fp   - py.rdkit.DataStructs.cDataStructs.ExplicitBitVect
%            167-bit MACCS keys fingerprint (Python reference, ADR-002).
%            Bit 0 is always 0 (unused). Bits 1-166 are MACCS structural keys.
%            Use emk.fingerprint.toArray(fp) to convert to logical(1, 167).
%
%   Errors:
%     emk:fingerprint:maccs:invalidInput  - mol is not a Mol object
%     emk:fingerprint:maccs:rdkitError    - unexpected Python exception
%
%   Example:
%     mol = emk.mol.fromSmiles("c1ccccc1");   % benzene
%     fp  = emk.fingerprint.maccs(mol);       % 167-bit MACCS keys FP
%     bits = emk.fingerprint.toArray(fp);     % logical(1, 167)
%
%   See also: emk.fingerprint.morgan, emk.fingerprint.toArray,
%             emk.similarity.tanimoto, emk.mol.fromSmiles

    % --- Input validation: mol type (before any RDKit call) ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:fingerprint:maccs:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    logDebug("maccs: mol has %d heavy atoms", double(mol.GetNumHeavyAtoms()));

    % --- Call RDKit ---
    % GenMACCSKeys returns a 167-bit ExplicitBitVect.
    % Use emk.util.rdkitModule() to get the module via importlib.
    % Direct py.rdkit.* access fails (TypeError); importlib avoids this.
    mods = emk.util.rdkitModule();
    try
        fp = mods.MACCSkeys.GenMACCSKeys(mol);
    catch ME
        error("emk:fingerprint:maccs:rdkitError", ...
            "RDKit GenMACCSKeys raised an exception: %s", ...
            ME.message);
    end

    logDebug("maccs: generated fp with %d bits, %d bits ON", ...
        double(fp.GetNumBits()), double(fp.GetNumOnBits()));
end
