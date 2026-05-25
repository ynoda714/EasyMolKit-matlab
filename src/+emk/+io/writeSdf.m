function writeSdf(mols, filePath)
% writeSdf  Write a cell array of RDKit Mol objects to an SDF file.
%
%   emk.io.writeSdf(mols, filePath)
%
%   Writes one or more RDKit Mol objects to an SD file (Structure-Data File)
%   using RDKit's SDWriter.  If the parent directory does not exist, an
%   error is raised.  Existing files are overwritten.
%
%   Arguments:
%     mols      - (cell array) Cell array of py.rdkit.Chem.rdchem.Mol objects.
%                 Must not be empty.  Each element must be a Python Mol object.
%     filePath  - (string | char) Path where the SDF file will be written.
%                 The parent directory must already exist.
%
%   Errors:
%     emk:io:writeSdf:invalidInput    - mols or filePath has wrong type/shape
%     emk:io:writeSdf:invalidMol      - an element in mols is not a Python Mol
%     emk:io:writeSdf:dirNotFound     - parent directory does not exist
%     emk:io:writeSdf:rdkitError      - unexpected Python exception
%
%   Example:
%     mol1 = emk.mol.fromSmiles("CCO");
%     mol2 = emk.mol.fromSmiles("c1ccccc1");
%     runDir = makeRunDir();
%     emk.io.writeSdf({mol1, mol2}, fullfile(runDir, "output.sdf"));
%
%   See also: emk.io.readSdf, emk.mol.fromSmiles

    % --- Input validation: filePath ---
    if ~(ischar(filePath) || isStringScalar(filePath))
        error("emk:io:writeSdf:invalidInput", ...
            "filePath must be a string scalar, got: %s", class(filePath));
    end

    filePath = string(filePath);

    % --- Input validation: mols ---
    if ~iscell(mols)
        error("emk:io:writeSdf:invalidInput", ...
            "mols must be a cell array, got: %s", class(mols));
    end
    if isempty(mols)
        error("emk:io:writeSdf:invalidInput", ...
            "mols must not be empty");
    end

    % --- Validate each Mol element before opening the writer ---
    for i = 1:numel(mols)
        if ~startsWith(class(mols{i}), "py.")
            error("emk:io:writeSdf:invalidMol", ...
                "mols{%d} is not a Python Mol object, got: %s", i, class(mols{i}));
        end
    end

    % --- Validate parent directory ---
    parentDir = fileparts(filePath);
    if strlength(parentDir) > 0 && ~isfolder(parentDir)
        error("emk:io:writeSdf:dirNotFound", ...
            "Parent directory does not exist: %s", parentDir);
    end

    logDebug("writeSdf: writing %d molecules to '%s'", numel(mols), filePath);

    % --- Write SDF via MolToMolBlock (OutOfProcess-safe) ---
    % SDWriter.write(mol) cannot pass Boost.Python C++ extension objects
    % across MATLAB's OutOfProcess IPC boundary (pickle cannot serialize
    % rdchem.Mol; MATLAB raises "class 'py.Boost.Python.class' not found").
    % MolToMolBlock(mol) returns a plain Python string, which safely crosses
    % the IPC boundary. We assemble the SDF file (molblock + "$$$$\n") in
    % MATLAB using native file I/O.
    mods = emk.util.rdkitModule();
    fid  = fopen(char(filePath), 'w');
    if fid == -1
        error("emk:io:writeSdf:rdkitError", ...
            "Cannot open file for writing: %s", filePath);
    end
    try
        for i = 1:numel(mols)
            mb = char(mods.Chem.MolToMolBlock(mols{i}));
            fwrite(fid, [mb, sprintf('$$$$\n')], 'char');
        end
        fclose(fid);
    catch ME
        fclose(fid);
        error("emk:io:writeSdf:rdkitError", ...
            "RDKit MolToMolBlock raised an exception: %s", ME.message);
    end

    logInfo("writeSdf: wrote %d molecules to '%s'", numel(mols), filePath);
end
