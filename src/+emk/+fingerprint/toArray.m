function bits = toArray(fp)
% toArray  Convert a Python fingerprint bit vector to a MATLAB logical array.
%
%   bits = emk.fingerprint.toArray(fp)
%
%   Converts a Python ExplicitBitVect (returned by emk.fingerprint.morgan,
%   emk.fingerprint.maccs, etc.) to a MATLAB logical row vector of length N,
%   where N = total bit count.  Element bits(i) is true when Python bit (i-1)
%   is set in the fingerprint (Python 0-indexed -> MATLAB 1-indexed).
%
%   Conversion strategy: a single ToBitString() IPC call fetches all bits
%   as a '0'/'1' character string, then char comparison converts to logical.
%   This minimizes IPC round-trips (ADR-002 IPC minimization principle).
%
%   Arguments:
%     fp    - Python fingerprint object (ExplicitBitVect)
%             Returned by emk.fingerprint.morgan() or emk.fingerprint.maccs()
%
%   Returns:
%     bits  - logical(1, N)  Row vector of length N = total bit count.
%             bits(i) == true means bit index (i-1) is ON in the fingerprint.
%
%   Errors:
%     emk:fingerprint:toArray:invalidInput  - fp is not a Python fingerprint
%     emk:fingerprint:toArray:rdkitError    - unexpected Python exception
%
%   Example:
%     mol  = emk.mol.fromSmiles("CCO");
%     fp   = emk.fingerprint.morgan(mol);
%     bits = emk.fingerprint.toArray(fp);  % logical(1, 2048)
%     fprintf("%d bits ON out of %d\n", sum(bits), numel(bits));
%     onIdx = find(bits);  % 1-indexed positions of ON bits
%
%   See also: emk.fingerprint.morgan, emk.fingerprint.maccs,
%             emk.similarity.tanimoto

    % --- Input validation: must be a Python object ---
    % Non-Python objects (double, string, char, logical, etc.) are rejected
    % before any Python call so that the error fires with no RDKit required.
    if ~startsWith(class(fp), "py.")
        error("emk:fingerprint:toArray:invalidInput", ...
            "fp must be a Python fingerprint object (ExplicitBitVect), got: %s", ...
            class(fp));
    end

    logDebug("toArray: converting fingerprint of class '%s'", class(fp));

    % --- Convert via ToBitString ---
    % ToBitString() returns a Python string of '0' and '1' characters,
    % e.g. "00001000101..." of length GetNumBits().
    % One IPC call fetches all bits at once (ADR-002 IPC minimization).
    %
    % If fp is a Python object but lacks ToBitString() (e.g. a Mol object
    % passed by mistake), the catch block converts the Python AttributeError
    % into a clear invalidInput error.
    try
        bitStr = char(string(fp.ToBitString()));
    catch ME
        error("emk:fingerprint:toArray:invalidInput", ...
            "fp does not support ToBitString(). " + ...
            "Expected ExplicitBitVect from emk.fingerprint.*, got: %s. " + ...
            "Detail: %s", class(fp), ME.message);
    end

    % Convert char array ('0'/'1') to logical row vector.
    % char comparison is O(N) and requires no additional Python calls.
    bits = bitStr == '1';

    logDebug("toArray: %d/%d bits ON", sum(bits), numel(bits));
end
