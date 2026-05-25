function result = verify()
% verify  Validate the Python environment and RDKit availability.
%
%   result = emk.setup.verify()
%
%   Checks that:
%     1. Python is configured in pyenv (or can be initialised via initPython)
%     2. rdkit.Chem can be imported from the configured Python environment
%
%   Returns:
%     result.python  (logical) - true if Python is configured and responsive
%     result.rdkit   (logical) - true if rdkit.Chem can be imported
%     result.version (string)  - Python version string, "" if unavailable
%
%   This function is non-throwing.  On failure it returns a struct with false
%   fields and logs a warning rather than raising an error.
%   Called automatically at the end of emk.setup.install().
%   Can also be called standalone to check the current environment status.
%
%   Error IDs: none (non-throwing by design)
%   See also: emk.setup.install, emk.setup.initPython

    result.python  = false;
    result.rdkit   = false;
    result.version = "";

    % Ensure Python is initialised
    pe = pyenv();
    if strcmp(string(pe.Status), "NotLoaded")
        try
            emk.setup.initPython();
        catch ME
            logWarn("verify: initPython failed: %s", ME.message);
            return;
        end
        pe = pyenv();
    end

    result.python  = true;
    result.version = string(pe.Version);
    logInfo("verify: Python %s OK (mode: %s)", result.version, string(pe.ExecutionMode));

    % Test RDKit import.
    % py.importlib.import_module starts the OutOfProcess Python worker if not
    % already running.  A missing or broken RDKit package raises a MATLAB
    % exception wrapping the Python ImportError; we catch it and log a warning
    % rather than throwing so that callers can decide how to handle the failure.
    try
        py.importlib.import_module("rdkit.Chem");
        result.rdkit = true;
        logInfo("verify: RDKit import OK.");
    catch ME
        logWarn("verify: RDKit import failed: %s", ME.message);
    end
end
