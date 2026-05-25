function fp = morgan(mol, options)
% morgan  Generate a Morgan (ECFP) circular fingerprint for a molecule.
%
%   fp = emk.fingerprint.morgan(mol)
%   fp = emk.fingerprint.morgan(mol, Radius=2, NBits=2048)
%
%   Computes a folded Morgan fingerprint bit vector using RDKit's
%   GetMorganFingerprintAsBitVect.  The fingerprint is returned as a Python
%   ExplicitBitVect reference (ADR-002: Python object reference retention).
%   Pass this object to emk.similarity.tanimoto() for similarity scoring,
%   or convert to a MATLAB logical array via emk.fingerprint.toArray().
%
%   Radius=2 (default) corresponds to ECFP4, which is the most widely used
%   setting in drug-discovery applications.  Radius=3 corresponds to ECFP6.
%
%   Arguments:
%     mol    - py.rdkit.Chem.rdchem.Mol  RDKit molecule object (Python ref)
%     Radius - int  (default: 2)   Morgan radius.  2 => ECFP4, 3 => ECFP6.
%     NBits  - int  (default: 2048) Bit vector length.
%
%   Returns:
%     fp     - py.rdkit.DataStructs.cDataStructs.ExplicitBitVect
%              Folded bit-vector fingerprint (Python reference, ADR-002).
%              Use emk.fingerprint.toArray(fp) to convert to logical array.
%
%   Errors:
%     emk:fingerprint:morgan:invalidInput  - mol is not a Mol object
%     emk:fingerprint:morgan:rdkitError    - unexpected Python exception
%
%   Example:
%     mol = emk.mol.fromSmiles("CCO");
%     fp  = emk.fingerprint.morgan(mol);                        % ECFP4/2048
%     fp  = emk.fingerprint.morgan(mol, Radius=3, NBits=1024);  % ECFP6/1024
%
%   See also: emk.fingerprint.toArray, emk.similarity.tanimoto,
%             emk.mol.fromSmiles

    arguments
        mol
        options.Radius (1,1) {mustBeInteger, mustBeNonnegative} = 2
        options.NBits  (1,1) {mustBeInteger, mustBePositive}    = 2048
    end

    % --- Input validation: mol type (before any RDKit call) ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:fingerprint:morgan:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    radius = int32(options.Radius);
    nbits  = int32(options.NBits);

    logDebug("morgan: radius=%d, nBits=%d, mol has %d heavy atoms", ...
        double(radius), double(nbits), double(mol.GetNumHeavyAtoms()));

    % --- Call RDKit ---
    % Use the rdFingerprintGenerator API (recommended since RDKit 2022+).
    % Use emk.util.rdkitModule() to get the module via importlib.
    % Direct py.rdkit.* access fails (TypeError); importlib avoids this.
    mods = emk.util.rdkitModule();
    try
        gen = mods.rdFpGen.GetMorganGenerator( ...
            pyargs("radius", radius, "fpSize", nbits));
        fp = gen.GetFingerprint(mol);
    catch ME
        error("emk:fingerprint:morgan:rdkitError", ...
            "RDKit MorganGenerator raised an exception: %s", ...
            ME.message);
    end

    logDebug("morgan: generated fp with %d bits, %d bits ON", ...
        double(fp.GetNumBits()), double(fp.GetNumOnBits()));
end
