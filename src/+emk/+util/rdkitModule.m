function mods = rdkitModule()
% rdkitModule  Return cached references to rdkit submodules.
%
%   mods = emk.util.rdkitModule()
%
%   Uses py.importlib.import_module() to reliably load rdkit submodules.
%   Direct py.rdkit.* access fails in MATLAB with TypeError because MATLAB
%   attempts to call the module object as a callable; importlib avoids this.
%   Module references are cached in a persistent variable so that the Python
%   IPC round-trip occurs only once per MATLAB session.
%
%   Returns a struct with fields:
%     mods.Chem             - rdkit.Chem module
%     mods.DataStructs      - rdkit.DataStructs module
%     mods.Descriptors      - rdkit.Chem.Descriptors module
%     mods.rdMolDescriptors - rdkit.Chem.rdMolDescriptors module
%     mods.rdFpGen          - rdkit.Chem.rdFingerprintGenerator module
%     mods.MACCSkeys        - rdkit.Chem.MACCSkeys module
%
%   Errors:
%     emk:util:rdkitModule:importFailed - py.importlib.import_module failed
%
%   Notes:
%     - The cache is invalidated when the function is cleared (clear all /
%       clear functions).  If pyenv is restarted, call
%       clear emk.util.rdkitModule to force re-import.
%     - scaffold.m and draw2d.m use their own local importlib calls for
%       modules not listed here (MurckoScaffold, AllChem, Draw).

    persistent cache;

    if ~isempty(cache)
        mods = cache;
        return;
    end

    logDebug("rdkitModule: importing rdkit submodules via importlib");

    try
        s.Chem             = py.importlib.import_module("rdkit.Chem");
        s.DataStructs      = py.importlib.import_module("rdkit.DataStructs");
        s.Descriptors      = py.importlib.import_module("rdkit.Chem.Descriptors");
        s.rdMolDescriptors = py.importlib.import_module("rdkit.Chem.rdMolDescriptors");
        s.rdFpGen          = py.importlib.import_module("rdkit.Chem.rdFingerprintGenerator");
        s.MACCSkeys        = py.importlib.import_module("rdkit.Chem.MACCSkeys");
    catch ME
        error("emk:util:rdkitModule:importFailed", ...
            "Failed to import rdkit submodules: %s\n" + ...
            "Ensure RDKit is installed (run emk.setup.install or emk.setup.installOnline).", ...
            ME.message);
    end

    cache = s;
    mods  = s;

    logDebug("rdkitModule: all submodules cached successfully");
end
