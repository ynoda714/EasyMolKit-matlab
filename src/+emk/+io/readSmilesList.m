function mols = readSmilesList(filePath)
% readSmilesList  Read a SMILES list file and return Mol objects.
%
%   mols = emk.io.readSmilesList(filePath)
%
%   Reads a plain-text SMILES file (one SMILES per line) and converts each
%   SMILES string to a py.rdkit.Chem.rdchem.Mol object via emk.mol.fromSmiles.
%
%   File format rules:
%     - Blank lines are skipped.
%     - Lines starting with '#' are treated as comments and skipped.
%     - Tab- or space-separated format is supported; only the first token
%       (the SMILES string) is used -- the optional name column is ignored.
%     - SMILES strings that fail to parse are skipped and logged as warnings.
%
%   Arguments:
%     filePath  - (string | char) Path to the SMILES list file.
%                 Absolute or relative path.  Must exist.
%
%   Returns:
%     mols      - (1 x N cell) Cell array of py.rdkit.Chem.rdchem.Mol objects
%                 (Python references).  N = number of successfully parsed mols.
%                 Empty cell {} if no valid molecules are found.
%
%   Errors:
%     emk:io:readSmilesList:invalidInput    - filePath is not a string/char scalar
%     emk:io:readSmilesList:fileNotFound    - file does not exist
%     emk:io:readSmilesList:allLinesFailed  - all non-comment SMILES lines
%                                               failed to parse
%
%   Example:
%     mols = emk.io.readSmilesList("data/list/sample.txt");
%     mw   = emk.descriptor.molWeight(mols{1});
%
%   See also: emk.io.readSdf, emk.mol.fromSmiles, emk.descriptor.calculate

    % --- Input validation ---
    if ~(ischar(filePath) || isStringScalar(filePath))
        error("emk:io:readSmilesList:invalidInput", ...
            "filePath must be a string scalar or char, got: %s", class(filePath));
    end

    filePath = string(filePath);

    if ~isfile(filePath)
        error("emk:io:readSmilesList:fileNotFound", ...
            "SMILES list file not found: %s", filePath);
    end

    logDebug("readSmilesList: reading '%s'", filePath);

    % --- Read all lines ---
    lines = readlines(filePath);

    % --- Parse each line and convert to Mol ---
    mols    = {};
    nSkip   = 0;
    nActive = 0;   % non-blank, non-comment lines attempted
    for i = 1:numel(lines)
        raw = strtrim(lines(i));

        % Skip blank lines and comment lines
        if strlength(raw) == 0 || startsWith(raw, "#")
            continue;
        end
        nActive = nActive + 1;

        % Extract first whitespace-delimited token as the SMILES string
        tokens = strsplit(raw, {' ', sprintf('\t')}, "CollapseDelimiters", true);
        smi    = tokens(1);

        % Convert SMILES to Mol; skip and warn on failure
        try
            mol = emk.mol.fromSmiles(smi);
            mols{end+1} = mol; %#ok<AGROW>
        catch ME
            nSkip = nSkip + 1;
            logWarn("readSmilesList: line %d, SMILES '%s' skipped (%s)", ...
                i, smi, ME.message);
        end
    end

    % --- M3-1: Error if every attempted line failed ---
    if nActive > 0 && isempty(mols)
        error("emk:io:readSmilesList:allLinesFailed", ...
            "All %d SMILES line(s) in '%s' failed to parse. " + ...
            "Verify the file format or check for invalid SMILES strings.", ...
            nActive, filePath);
    end

    logInfo("readSmilesList: loaded %d molecules from '%s' (%d skipped)", ...
        numel(mols), filePath, nSkip);
end
