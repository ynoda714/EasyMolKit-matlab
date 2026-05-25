function score = tanimoto(fp1, fp2)
% tanimoto  Compute the Tanimoto similarity between two fingerprints.
%
%   score = emk.similarity.tanimoto(fp1, fp2)
%
%   Computes the Tanimoto coefficient (Jaccard index for binary vectors)
%   between two Python fingerprint bit vectors using RDKit's
%   DataStructs.TanimotoSimilarity.  Both fingerprints must have the same
%   bit length; mixing Morgan (2048-bit) and MACCS (167-bit) fingerprints
%   will raise an error from RDKit.
%
%   Mathematical definition:
%
%     T(A, B) = |A AND B| / |A OR B| = c / (a + b - c)
%
%   where a = number of ON bits in A, b = in B, c = ON bits in both.
%   Returns 0.0 when both vectors are all-zero (by convention).
%   Value range: [0, 1].  1.0 = identical, 0.0 = no common ON bits.
%
%   Arguments:
%     fp1   - py.rdkit.DataStructs.cDataStructs.ExplicitBitVect
%             First fingerprint (Python reference from emk.fingerprint.*)
%     fp2   - py.rdkit.DataStructs.cDataStructs.ExplicitBitVect
%             Second fingerprint (Python reference from emk.fingerprint.*)
%
%   Returns:
%     score - double scalar in [0, 1].  Tanimoto similarity.
%
%   Errors:
%     emk:similarity:tanimoto:invalidInput  - fp1 or fp2 is not Python object
%     emk:similarity:tanimoto:rdkitError    - size mismatch or Python exception
%
%   Example:
%     mol1  = emk.mol.fromSmiles("CCO");
%     mol2  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
%     fp1   = emk.fingerprint.morgan(mol1);
%     fp2   = emk.fingerprint.morgan(mol2);
%     score = emk.similarity.tanimoto(fp1, fp2);  % double in [0, 1]
%
%   See also: emk.fingerprint.morgan, emk.fingerprint.maccs,
%             emk.fingerprint.toArray, emk.mol.fromSmiles

    % --- Input validation: must be Python objects (before any RDKit call) ---
    % Non-Python objects (double, string, char, logical, etc.) are rejected
    % here so that the error fires without requiring RDKit.
    if ~startsWith(class(fp1), "py.")
        error("emk:similarity:tanimoto:invalidInput", ...
            "fp1 must be a Python fingerprint object (ExplicitBitVect), got: %s", ...
            class(fp1));
    end
    if ~startsWith(class(fp2), "py.")
        error("emk:similarity:tanimoto:invalidInput", ...
            "fp2 must be a Python fingerprint object (ExplicitBitVect), got: %s", ...
            class(fp2));
    end

    logDebug("tanimoto: fp1 class=%s, fp2 class=%s", class(fp1), class(fp2));

    % --- Call RDKit ---
    % Use emk.util.rdkitModule() to get DataStructs via importlib.
    % Direct py.rdkit.* access fails (TypeError); importlib avoids this.
    mods = emk.util.rdkitModule();
    try
        score = double(mods.DataStructs.TanimotoSimilarity(fp1, fp2));
    catch ME
        error("emk:similarity:tanimoto:rdkitError", ...
            "RDKit TanimotoSimilarity raised an exception: %s", ...
            ME.message);
    end

    logDebug("tanimoto: score=%.4f", score);
end
