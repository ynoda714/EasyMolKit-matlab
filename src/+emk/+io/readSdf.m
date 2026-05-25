function mols = readSdf(filePath)
% readSdf  Read an SDF file and return a cell array of RDKit Mol objects.
%
%   mols = emk.io.readSdf(filePath)
%
%   Reads molecules from an SD file (Structure-Data File) using RDKit's
%   SDMolSupplier.  Molecules that RDKit cannot parse are silently skipped
%   and a warning is logged for each skipped entry.
%
%   Arguments:
%     filePath  - (string | char) Path to the SDF file.
%                 Absolute or relative path.  Must exist.
%
%   Returns:
%     mols      - (1 x N cell) Cell array of py.rdkit.Chem.rdchem.Mol objects
%                 (Python references).  N = number of successfully parsed mols.
%                 Empty cell {} if no valid molecules are found.
%
%   Errors:
%     emk:io:readSdf:invalidInput  - filePath is not a string/char scalar
%     emk:io:readSdf:fileNotFound  - file does not exist
%     emk:io:readSdf:rdkitError    - unexpected Python exception
%
%   Example:
%     mols = emk.io.readSdf("data/sample.sdf");
%     mw   = emk.descriptor.molWeight(mols{1});
%
%   See also: emk.io.writeSdf, emk.mol.fromSmiles, emk.descriptor.calculate

    % --- Input validation ---
    if ~(ischar(filePath) || isStringScalar(filePath))
        error("emk:io:readSdf:invalidInput", ...
            "filePath must be a string scalar, got: %s", class(filePath));
    end

    filePath = string(filePath);

    if ~isfile(filePath)
        error("emk:io:readSdf:fileNotFound", ...
            "SDF file not found: %s", filePath);
    end

    logDebug("readSdf: reading '%s'", filePath);

    % --- Read SDF via MolFromMolBlock (OutOfProcess-safe) ---
    % SDMolSupplier returns a Boost.Python iterable object. MATLAB's
    % OutOfProcess IPC cannot deserialize it (pickle fails on Boost.Python
    % C++ extension types), causing "class 'py.Boost.Python.class' not found"
    % Instead, read the SDF as plain text in MATLAB,
    % split by the '$$$$' record separator, and call MolFromMolBlock() for
    % each block -- which returns a single Mol reference (IPC-safe, same
    % pattern as fromSmiles).
    mods = emk.util.rdkitModule();

    try
        sdfText   = fileread(char(filePath));
        % Split on the full record terminator "$$$$\n" so that each block
        % retains its leading newline (= empty molecule-name line in V2000
        % format). Using bare "$$$$" as the delimiter causes strtrim to strip
        % that leading newline, which shifts V2000 header rows and makes
        % MolFromMolBlock return None for every molecule.
        rawBlocks = strsplit(sdfText, sprintf('$$$$\n'));
    catch ME
        error("emk:io:readSdf:rdkitError", ...
            "Failed to read SDF file '%s': %s", filePath, ME.message);
    end

    % --- Parse each mol block ---
    mols  = {};
    nSkip = 0;
    for i = 1:numel(rawBlocks)
        blk = rawBlocks{i};
        if strlength(strtrim(blk)) == 0
            continue;   % trailing empty element after last $$$$
        end
        try
            mol = mods.Chem.MolFromMolBlock(char(blk));
        catch ME
            nSkip = nSkip + 1;
            logWarn("readSdf: block #%d MolFromMolBlock error: %s (skipping)", i, ME.message);
            continue;
        end
        if isa(mol, "py.NoneType")
            nSkip = nSkip + 1;
            logWarn("readSdf: molecule #%d could not be parsed (None), skipping", i);
        else
            mols{end+1} = mol; %#ok<AGROW>
        end
    end

    logInfo("readSdf: loaded %d molecules from '%s' (%d skipped)", ...
        numel(mols), filePath, nSkip);
end
