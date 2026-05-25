function tbl = batchCalculate(mols, descriptorNames)
% batchCalculate  Compute descriptors for multiple molecules and return a table.
%
%   tbl = emk.descriptor.batchCalculate(mols)
%   tbl = emk.descriptor.batchCalculate(mols, ["MolWt","LogP","TPSA"])
%
%   Computes one or more RDKit-based physicochemical descriptors for each
%   molecule in a cell array and returns the results as a MATLAB table.
%   Invalid or null molecules in the cell array are skipped; the
%   corresponding row is filled with NaN values and a warning is logged.
%
%   Performance: Uses a Python batch helper (batch_descriptors.py) to
%   compute all descriptors in a single IPC round-trip per ADR-002 rev.3.
%   Falls back to a per-molecule loop if the helper is unavailable.
%
%   Arguments:
%     mols            - (1 x N cell) Cell array of py.rdkit.Chem.rdchem.Mol
%                       objects, e.g., from emk.io.readSdf or a loop over
%                       emk.mol.fromSmiles.
%     descriptorNames - (string array, optional) Names of descriptors to
%                       compute.  Defaults to all 10 supported descriptors.
%
%   Returns:
%     tbl  - (table) N-row table.  Column names = descriptor names.
%            Each cell contains a double scalar.  Rows for invalid
%            molecules contain NaN in all descriptor columns.
%
%   Supported descriptor names: same as emk.descriptor.calculate.
%
%   Errors:
%     emk:descriptor:batchCalculate:invalidInput    - mols is not a cell array
%     emk:descriptor:batchCalculate:unknownDescriptor - unrecognized name(s)
%     emk:descriptor:batchCalculate:allMolsFailed   - all molecules failed
%
%   Example:
%     smilesList = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
%     mols = cellfun(@(s) emk.mol.fromSmiles(s), smilesList, ...
%                    "UniformOutput", false);
%     tbl  = emk.descriptor.batchCalculate(mols);
%     tbl  = emk.descriptor.batchCalculate(mols, ["MolWt","LogP","TPSA"]);
%
%   See also: emk.descriptor.calculate, emk.io.readSdf,
%             emk.io.readSmilesList

    SUPPORTED_NAMES = supportedNames_();

    % --- Input validation: mols ---
    if ~iscell(mols)
        error("emk:descriptor:batchCalculate:invalidInput", ...
            "mols must be a cell array, got: %s", class(mols));
    end

    % --- Input validation: descriptorNames ---
    if nargin < 2 || isempty(descriptorNames)
        descriptorNames = SUPPORTED_NAMES;
    else
        if ischar(descriptorNames)
            descriptorNames = string(descriptorNames);
        end
        descriptorNames = reshape(string(descriptorNames), 1, []);
        unknown = descriptorNames(~ismember(descriptorNames, SUPPORTED_NAMES));
        if ~isempty(unknown)
            error("emk:descriptor:batchCalculate:unknownDescriptor", ...
                "Unknown descriptor(s): [%s]. Supported: [%s]", ...
                strjoin(unknown, ", "), strjoin(SUPPORTED_NAMES, ", "));
        end
    end

    nMols = numel(mols);
    nDesc = numel(descriptorNames);

    logDebug("batchCalculate: %d mol(s), %d descriptor(s)", nMols, nDesc);

    % --- Preallocate output matrix (NaN = placeholder for invalid mols) ---
    data  = NaN(nMols, nDesc);
    nSkip = 0;

    % --- Compute descriptors ---
    if nMols > 0
        [data, nSkip] = computeBatch_(mols, descriptorNames, nMols, nDesc, data);
    end

    % --- M3-1: Error if every molecule failed ---
    if nMols > 0 && nSkip == nMols
        error("emk:descriptor:batchCalculate:allMolsFailed", ...
            "All %d molecule(s) failed descriptor calculation. " + ...
            "Ensure mols contains valid py.rdkit.Chem.rdchem.Mol objects.", nMols);
    end

    if nSkip > 0
        logWarn("batchCalculate: %d / %d molecule(s) produced NaN rows", nSkip, nMols);
    end

    % --- Build output table ---
    tbl = array2table(data, "VariableNames", cellstr(descriptorNames));

    logDebug("batchCalculate: done (%d rows, %d cols)", nMols, nDesc);
end

% =========================================================================
% Private helpers
% =========================================================================

function [data, nSkip] = computeBatch_(mols, descriptorNames, nMols, nDesc, data)
% Attempt Python batch computation (single IPC call per ADR-002 rev.3).
% Falls back to per-molecule loop when the batch helper is unavailable.

    % -- Load Python batch helper (M3-2: IPC optimisation) --
    helperMod = loadBatchHelper_();
    if isempty(helperMod)
        [data, nSkip] = computeLoop_(mols, descriptorNames, nMols, nDesc, data);
        return;
    end

    % -- Type-check in MATLAB (no Python call required) --
    validIdx = false(1, nMols);
    for i = 1:nMols
        if isa(mols{i}, "py.rdkit.Chem.rdchem.Mol")
            validIdx(i) = true;
        else
            logWarn("batchCalculate: mol %d is not a valid Mol object (%s) -- filled with NaN", ...
                i, class(mols{i}));
        end
    end

    nTypeSkip = sum(~validIdx);
    nValid    = sum(validIdx);

    if nValid == 0
        nSkip = nMols;
        return;
    end

    % -- Build py.list of valid mols --
    validMols = mols(validIdx);
    pyMols    = py.list();
    for k = 1:nValid
        pyMols.append(validMols{k});
    end

    pyNames = py.list();
    for k = 1:nDesc
        pyNames.append(char(descriptorNames(k)));
    end

    % -- Single IPC call to Python batch helper --
    try
        pyResults  = helperMod.batch_calculate(pyMols, pyNames);
        pyRowsCell = cell(pyResults);
        vi     = 1;
        nPySkip = 0;
        for i = 1:nMols
            if validIdx(i)
                rowVec    = double(py.array.array("d", pyRowsCell{vi}));
                data(i,:) = rowVec;
                if all(isnan(rowVec))
                    nPySkip = nPySkip + 1;
                    logWarn("batchCalculate: descriptor calculation failed for mol %d -- filled with NaN", i);
                end
                vi = vi + 1;
            end
        end
        nSkip = nTypeSkip + nPySkip;

    catch ME
        logWarn("batchCalculate: Python batch helper raised an exception (%s) -- using loop fallback", ...
            ME.message);
        [data, nSkip] = computeLoop_(mols, descriptorNames, nMols, nDesc, data);
    end
end

% -------------------------------------------------------------------------
function [data, nSkip] = computeLoop_(mols, descriptorNames, nMols, nDesc, data)
% Per-molecule fallback loop (N IPC round-trips).

    nSkip = 0;
    for i = 1:nMols
        mol = mols{i};
        if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
            logWarn("batchCalculate: mol %d is not a valid Mol object (%s) -- filled with NaN", ...
                i, class(mol));
            nSkip = nSkip + 1;
            continue;
        end
        try
            desc = emk.descriptor.calculate(mol, descriptorNames);
            for j = 1:nDesc
                data(i, j) = desc.(descriptorNames(j));
            end
        catch ME
            logWarn("batchCalculate: descriptor calculation failed for mol %d -- filled with NaN: %s", ...
                i, ME.message);
            nSkip = nSkip + 1;
        end
    end
end

% -------------------------------------------------------------------------
function helperMod = loadBatchHelper_()
% Import batch_descriptors Python module (cached after first call).
% Adds src/+emk/+util/python/ to sys.path on first call.
% ADR-002 rev.3: IPC minimisation for batch APIs.

    persistent cachedMod;
    if ~isempty(cachedMod)
        helperMod = cachedMod;
        return;
    end

    helperMod = [];
    try
        thisFile  = mfilename("fullpath");
        % batchCalculate.m resides at src/+emk/+descriptor/
        % Helper dir is at      src/+emk/+util/python/
        helperDir = fullfile(fileparts(fileparts(thisFile)), "+util", "python");
        py.sys.path().insert(int32(0), helperDir);
        mod       = py.importlib.import_module("batch_descriptors");
        cachedMod = mod;
        helperMod = mod;
        logDebug("batchCalculate: batch helper loaded from '%s'", helperDir);
    catch ME
        logWarn("batchCalculate: batch helper unavailable (%s) -- using loop fallback", ME.message);
    end
end

% -------------------------------------------------------------------------
function names = supportedNames_()
    names = ["MolWt", "ExactMolWt", "LogP", "TPSA", ...
             "NumHAcceptors", "NumHDonors", "NumRotatableBonds", ...
             "RingCount", "FractionCSP3", "HeavyAtomCount"];
end
