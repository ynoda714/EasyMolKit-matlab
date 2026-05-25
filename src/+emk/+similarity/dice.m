function score = dice(fp1, fp2)
% dice  Compute the Dice similarity between two fingerprints.
%
%   score = emk.similarity.dice(fp1, fp2)
%
%   Computes the Dice coefficient (Sorensen-Dice index for binary vectors)
%   between two Python fingerprint bit vectors using RDKit's
%   DataStructs.DiceSimilarity.  Both fingerprints must have the same
%   bit length; mixing Morgan (2048-bit) and MACCS (167-bit) fingerprints
%   will raise an error from RDKit.
%
%   Mathematical definition:
%
%     D(A, B) = 2|A AND B| / (|A| + |B|) = 2c / (a + b)
%
%   where a = number of ON bits in A, b = in B, c = ON bits in both.
%   Returns 0.0 when both vectors are all-zero (by convention).
%   Value range: [0, 1].  1.0 = identical, 0.0 = no common ON bits.
%
%   Note: Dice >= Tanimoto always holds for binary vectors.
%
%   Arguments:
%     fp1   - py.rdkit.DataStructs.cDataStructs.ExplicitBitVect
%             First fingerprint (Python reference from emk.fingerprint.*)
%     fp2   - py.rdkit.DataStructs.cDataStructs.ExplicitBitVect
%             Second fingerprint (Python reference from emk.fingerprint.*)
%
%   Returns:
%     score - double scalar in [0, 1].  Dice similarity.
%
%   Errors:
%     emk:similarity:dice:invalidInput  - fp1 or fp2 is not Python object
%     emk:similarity:dice:rdkitError    - size mismatch or Python exception
%
%   Example:
%     mol1  = emk.mol.fromSmiles("CCO");
%     mol2  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
%     fp1   = emk.fingerprint.morgan(mol1);
%     fp2   = emk.fingerprint.morgan(mol2);
%     score = emk.similarity.dice(fp1, fp2);  % double in [0, 1]
%
%   See also: emk.similarity.tanimoto, emk.fingerprint.morgan,
%             emk.fingerprint.maccs, emk.fingerprint.toArray

    % --- Input validation: must be Python objects (before any RDKit call) ---
    if ~startsWith(class(fp1), "py.")
        error("emk:similarity:dice:invalidInput", ...
            "fp1 must be a Python fingerprint object (ExplicitBitVect), got: %s", ...
            class(fp1));
    end
    if ~startsWith(class(fp2), "py.")
        error("emk:similarity:dice:invalidInput", ...
            "fp2 must be a Python fingerprint object (ExplicitBitVect), got: %s", ...
            class(fp2));
    end

    logDebug("dice: fp1 class=%s, fp2 class=%s", class(fp1), class(fp2));

    % --- Call RDKit ---
    % Use emk.util.rdkitModule() to get DataStructs via importlib.
    % Direct py.rdkit.* access fails (TypeError); importlib avoids this.
    mods = emk.util.rdkitModule();
    try
        score = double(mods.DataStructs.DiceSimilarity(fp1, fp2));
    catch ME
        error("emk:similarity:dice:rdkitError", ...
            "RDKit DiceSimilarity raised an exception: %s", ...
            ME.message);
    end

    logDebug("dice: score=%.4f", score);
end
