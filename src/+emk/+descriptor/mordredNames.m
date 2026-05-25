function names = mordredNames()
% mordredNames  Return all available 2D Mordred descriptor names.
%
%   names = emk.descriptor.mordredNames()
%
%   Returns a sorted string array of all descriptor names that can be
%   passed to emk.descriptor.mordred() or emk.descriptor.mordredBatch().
%
%   Requires: emk.setup.installExtra("mordred") before first use.
%
%   Returns:
%     names  - (1 x N string) Sorted descriptor name strings.
%              Pass any subset to mordred() or mordredBatch() as the
%              second argument to compute only selected descriptors.
%
%   Errors:
%     emk:descriptor:mordredNames:libraryNotFound - mordred not installed
%     emk:descriptor:mordredNames:pythonError     - unexpected Python exception
%
%   Example:
%     allNames = emk.descriptor.mordredNames();
%     fprintf("Available Mordred descriptors: %d\n", numel(allNames));
%     % Select weight- and lipophilicity-related descriptors:
%     selected = allNames(startsWith(allNames, "MW") | ...
%                         startsWith(allNames, "ALogP"));
%     mol  = emk.mol.fromSmiles("CCO");
%     desc = emk.descriptor.mordred(mol, selected);
%
%   See also: emk.descriptor.mordred, emk.descriptor.mordredBatch,
%             emk.setup.installExtra

    % Verify mordred availability
    try
        py.importlib.import_module("mordred");
    catch ME
        if contains(ME.message, "ModuleNotFoundError") || ...
                contains(ME.message, "No module named")
            error("emk:descriptor:mordredNames:libraryNotFound", ...
                "Mordred is not installed in the active Python environment.\n" + ...
                "Install it with: emk.setup.installExtra(""mordred"")");
        end
        rethrow(ME);
    end

    % Load run_mordred helper
    try
        thisFile  = mfilename("fullpath");
        helperDir = fullfile(fileparts(fileparts(thisFile)), "+util", "python");
        py.sys.path().insert(int32(0), helperDir);
        mod   = py.importlib.import_module("run_mordred");
        names = string(cell(mod.mordred_list_names()));
    catch ME
        error("emk:descriptor:mordredNames:pythonError", ...
            "Failed to retrieve Mordred descriptor names: %s", ME.message);
    end
end
