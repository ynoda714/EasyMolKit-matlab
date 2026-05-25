function desc = mordred(mol, descriptorNames)
% mordred  Compute Mordred molecular descriptors for a single molecule.
%
%   desc = emk.descriptor.mordred(mol)
%   desc = emk.descriptor.mordred(mol, ["MW","ALogP","nRot"])
%
%   Computes 2D molecular descriptors using the Mordred library (mordredcommunity
%   fork).  When descriptorNames is omitted all available 2D descriptors
%   (~1800 values) are computed.  Failed descriptors are represented as NaN.
%
%   Requires: emk.setup.installExtra("mordred") before first use.
%   Python/RDKit must also be initialised (emk.setup.initPython).
%
%   Arguments:
%     mol             - py.rdkit.Chem.rdchem.Mol  RDKit molecule object.
%     descriptorNames - string array (optional)  Names of descriptors to
%                       compute.  Call emk.descriptor.mordredNames() to
%                       obtain the full list.  Default: all 2D descriptors.
%
%   Returns:
%     desc  - struct  Field names = descriptor names, values = double scalar.
%             Non-numeric or failed descriptors appear as NaN.
%
%   Errors:
%     emk:descriptor:mordred:invalidInput    - mol is not a Mol object
%     emk:descriptor:mordred:libraryNotFound - mordred not installed
%     emk:descriptor:mordred:pythonError     - unexpected Python exception
%
%   Performance note:
%     The first call builds and caches a Mordred Calculator inside Python.
%     Subsequent calls are faster.  Use mordredBatch() for multiple molecules.
%
%   Example:
%     mol  = emk.mol.fromSmiles("CCO");
%     desc = emk.descriptor.mordred(mol);
%     desc = emk.descriptor.mordred(mol, ["MW", "ALogP", "nRot"]);
%     fprintf("MW = %.3f  ALogP = %.3f\n", desc.MW, desc.ALogP);
%
%   See also: emk.descriptor.mordredBatch, emk.descriptor.mordredNames,
%             emk.descriptor.calculate, emk.setup.installExtra

    % --- Input validation: mol type (no Python required) ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:descriptor:mordred:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    % --- Optional: descriptorNames ---
    if nargin < 2 || isempty(descriptorNames)
        descriptorNames = string.empty;   % Python helper interprets [] as "all"
    else
        if ischar(descriptorNames)
            descriptorNames = string(descriptorNames);
        end
        descriptorNames = reshape(string(descriptorNames), 1, []);
    end

    % --- Load Python helper ---
    mod = loadMordredHelper_();

    % --- Call Python: single IPC round-trip via mordred_calculate_list ---
    logDebug("mordred: computing descriptors for mol (%d heavy atoms)", ...
        double(mol.GetNumHeavyAtoms()));

    pyNames = py.list();
    for k = 1:numel(descriptorNames)
        pyNames.append(char(descriptorNames(k)));
    end

    try
        pyTuple    = mod.mordred_calculate_list(mol, pyNames);
        tupleCell  = cell(pyTuple);
        colNames   = string(cell(tupleCell{1}));   % descriptor names
        values     = double(py.array.array("d", tupleCell{2}));
    catch ME
        error("emk:descriptor:mordred:pythonError", ...
            "Mordred raised an exception: %s", ME.message);
    end

    % --- Build output struct ---
    nDesc = numel(colNames);
    desc  = struct();
    for k = 1:nDesc
        desc.(colNames(k)) = values(k);
    end

    logDebug("mordred: computed %d descriptor(s).", nDesc);
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
            error("emk:descriptor:mordred:libraryNotFound", ...
                "Mordred is not installed in the active Python environment.\n" + ...
                "Install it with: emk.setup.installExtra(""mordred"")");
        end
        rethrow(ME);
    end

    % Load run_mordred helper (adds src/+emk/+util/python/ to sys.path)
    try
        thisFile  = mfilename("fullpath");
        % mordred.m resides at src/+emk/+descriptor/
        % Helper dir is at      src/+emk/+util/python/
        helperDir = fullfile(fileparts(fileparts(thisFile)), "+util", "python");
        py.sys.path().insert(int32(0), helperDir);
        mod       = py.importlib.import_module("run_mordred");
        cachedMod = mod;
        logDebug("mordred: helper loaded from '%s'", helperDir);
    catch ME
        error("emk:descriptor:mordred:pythonError", ...
            "Failed to load run_mordred Python helper: %s", ME.message);
    end
end
