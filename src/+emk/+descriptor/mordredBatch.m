function tbl = mordredBatch(mols, descriptorNames)
% mordredBatch  Compute Mordred descriptors for multiple molecules.
%
%   tbl = emk.descriptor.mordredBatch(mols)
%   tbl = emk.descriptor.mordredBatch(mols, ["MW","ALogP","nRot"])
%
%   Computes 2D molecular descriptors using the Mordred library for each
%   molecule in a cell array and returns a MATLAB table.
%   Invalid or null molecules yield NaN rows (a warning is logged).
%
%   Requires: emk.setup.installExtra("mordred") before first use.
%   Python/RDKit must also be initialised (emk.setup.initPython).
%
%   Performance: Uses a Python batch helper (run_mordred.py) to compute all
%   descriptors for all molecules in a single IPC round-trip (ADR-002 rev.3).
%
%   Arguments:
%     mols            - (1 x N cell) Cell array of py.rdkit.Chem.rdchem.Mol.
%                       Created by emk.mol.fromSmiles or emk.io.readSdf.
%     descriptorNames - (string array, optional) Names of descriptors to
%                       compute.  Default: all ~1800 2D Mordred descriptors.
%                       Call emk.descriptor.mordredNames() for the full list.
%
%   Returns:
%     tbl  - (table) N-row table.  Column names = descriptor names.
%            Non-numeric or failed descriptors appear as NaN.
%
%   Errors:
%     emk:descriptor:mordredBatch:invalidInput    - mols is not a cell array
%     emk:descriptor:mordredBatch:libraryNotFound - mordred not installed
%     emk:descriptor:mordredBatch:allMolsFailed   - all molecules failed
%     emk:descriptor:mordredBatch:pythonError     - unexpected Python exception
%
%   Example:
%     smilesList = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
%     mols = cellfun(@(s) emk.mol.fromSmiles(s), smilesList, ...
%                    "UniformOutput", false);
%     tbl  = emk.descriptor.mordredBatch(mols, ["MW","ALogP","nRot"]);
%     disp(tbl);
%
%   See also: emk.descriptor.mordred, emk.descriptor.mordredNames,
%             emk.descriptor.batchCalculate, emk.setup.installExtra

    % --- Input validation: mols (no Python required) ---
    if ~iscell(mols)
        error("emk:descriptor:mordredBatch:invalidInput", ...
            "mols must be a cell array, got: %s", class(mols));
    end

    % --- Optional: descriptorNames ---
    if nargin < 2 || isempty(descriptorNames)
        descriptorNames = string.empty;   % Python helper returns all 2D descriptors
    else
        if ischar(descriptorNames)
            descriptorNames = string(descriptorNames);
        end
        descriptorNames = reshape(string(descriptorNames), 1, []);
    end

    nMols = numel(mols);
    logDebug("mordredBatch: %d mol(s)", nMols);

    % --- Load Python helper ---
    mod = loadMordredHelper_();

    % --- Type-check mols in MATLAB (avoids passing invalid objects to Python) ---
    validIdx = false(1, nMols);
    for i = 1:nMols
        if isa(mols{i}, "py.rdkit.Chem.rdchem.Mol")
            validIdx(i) = true;
        else
            logWarn("mordredBatch: mol %d is not a valid Mol object (%s) -- filled with NaN", ...
                i, class(mols{i}));
        end
    end

    nValid    = sum(validIdx);
    nTypeSkip = sum(~validIdx);

    % --- Build py.list of valid mols (None placeholder for invalid) ---
    pyMols = py.list();
    for i = 1:nMols
        if validIdx(i)
            pyMols.append(mols{i});
        else
            pyMols.append(py.None);
        end
    end

    pyNames = py.list();
    for k = 1:numel(descriptorNames)
        pyNames.append(char(descriptorNames(k)));
    end

    % --- Single IPC call: mordred_batch_matrix returns (names, matrix) ---
    try
        pyTuple   = mod.mordred_batch_matrix(pyMols, pyNames);
        tupleCell = cell(pyTuple);
        colNames  = string(cell(tupleCell{1}));
        pyMatrix  = tupleCell{2};
    catch ME
        error("emk:descriptor:mordredBatch:pythonError", ...
            "Mordred batch helper raised an exception: %s", ME.message);
    end

    nDesc     = numel(colNames);
    data      = NaN(nMols, nDesc);
    pyRowCell = cell(pyMatrix);

    nPySkip = 0;
    for i = 1:nMols
        row = double(py.array.array("d", pyRowCell{i}));
        data(i, :) = row;
        if validIdx(i) && all(isnan(row))
            nPySkip = nPySkip + 1;
            logWarn("mordredBatch: descriptor computation failed for mol %d -- filled with NaN", i);
        end
    end

    nSkip = nTypeSkip + nPySkip;

    % --- M3-1 equivalent: error when ALL molecules produced NaN ---
    if nMols > 0 && nSkip == nMols
        error("emk:descriptor:mordredBatch:allMolsFailed", ...
            "All %d molecule(s) failed Mordred descriptor calculation. " + ...
            "Ensure mols contains valid py.rdkit.Chem.rdchem.Mol objects and " + ...
            "Mordred is installed (emk.setup.installExtra(""mordred"")).", nMols);
    end

    if nSkip > 0
        logWarn("mordredBatch: %d / %d molecule(s) produced NaN rows", nSkip, nMols);
    end

    % --- Build output table ---
    tbl = array2table(data, "VariableNames", cellstr(colNames));

    logDebug("mordredBatch: done (%d rows, %d descriptor cols).", nMols, nDesc);
end

% =========================================================================
% Private helpers
% =========================================================================

function mod = loadMordredHelper_()
% Import run_mordred Python module and verify mordred is available.
% Cached after the first successful import (persistent variable).

    persistent cachedMod;
    if ~isempty(cachedMod)
        mod = cachedMod;
        return;
    end

    % Verify mordred is installed before attempting import
    try
        py.importlib.import_module("mordred");
    catch ME
        if contains(ME.message, "ModuleNotFoundError") || ...
                contains(ME.message, "No module named")
            error("emk:descriptor:mordredBatch:libraryNotFound", ...
                "Mordred is not installed in the active Python environment.\n" + ...
                "Install it with: emk.setup.installExtra(""mordred"")");
        end
        rethrow(ME);
    end

    % Load run_mordred helper (adds src/+emk/+util/python/ to sys.path)
    try
        thisFile  = mfilename("fullpath");
        % mordredBatch.m resides at src/+emk/+descriptor/
        % Helper dir is at         src/+emk/+util/python/
        helperDir = fullfile(fileparts(fileparts(thisFile)), "+util", "python");
        py.sys.path().insert(int32(0), helperDir);
        mod       = py.importlib.import_module("run_mordred");
        cachedMod = mod;
        logDebug("mordredBatch: helper loaded from '%s'", helperDir);
    catch ME
        error("emk:descriptor:mordredBatch:pythonError", ...
            "Failed to load run_mordred Python helper: %s", ME.message);
    end
end
