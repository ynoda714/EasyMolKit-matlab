function tf = isValid(smiles)
% isValid  Check whether a SMILES string represents a valid molecule.
%
%   tf = emk.mol.isValid(smiles)
%
%   Returns true if smiles can be parsed by RDKit into a valid Mol object,
%   false otherwise.  Unlike fromSmiles, this function never throws for
%   parse failures or unsupported input types; it returns false instead.
%   Unexpected Python exceptions (environment errors) are wrapped and
%   re-thrown as rdkitError to surface genuine RDKit installation issues.
%
%   Arguments:
%     smiles  - (string | char) SMILES string to validate.
%               Non-string/char inputs return false without error.
%
%   Returns:
%     tf      - logical scalar. true if parseable, false otherwise.
%
%   Errors:
%     emk:mol:isValid:rdkitError  - unexpected Python exception
%                                     (indicates an RDKit environment issue,
%                                     not a parse failure)
%
%   Example:
%     tf = emk.mol.isValid("CCO");              % true
%     tf = emk.mol.isValid("XYZ");              % false
%     tf = emk.mol.isValid("");                 % false
%     tf = emk.mol.isValid(42);                 % false (no throw)
%
%   See also: emk.mol.fromSmiles, emk.mol.toSmiles

    % --- Non-string inputs are never valid SMILES; return false, no throw ---
    if ~(ischar(smiles) || isStringScalar(smiles))
        logDebug("isValid: non-string input (%s) -> false", class(smiles));
        tf = false;
        return;
    end

    logDebug("isValid: checking '%s'", string(smiles));

    % --- Delegate to fromSmiles; treat parse errors as false ---
    %
    % fromSmiles throws:
    %   invalidInput  - cannot happen here; type check above already passed
    %   invalidSmiles - empty string, whitespace, or RDKit returns None
    %                   -> return false (expected validation failure)
    %   rdkitError    - unexpected Python exception
    %                   -> re-throw so the caller knows RDKit is broken
    try
        emk.mol.fromSmiles(smiles);
        tf = true;
    catch ME
        if strcmp(ME.identifier, "emk:mol:fromSmiles:rdkitError")
            error("emk:mol:isValid:rdkitError", ...
                "RDKit raised an exception while validating '%s': %s", ...
                string(smiles), ME.message);
        end
        tf = false;
    end

    logDebug("isValid: '%s' -> %d", string(smiles), tf);
end
