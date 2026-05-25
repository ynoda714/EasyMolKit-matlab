function useExternal(pythonPath)
% useExternal  Connect EasyMolKit to an external Python environment (Track 2).
%
%   emk.setup.useExternal(pythonPath)
%
%   Configures pyenv to use the specified external CPython executable instead
%   of the repo-local Embedded Python.  Intended for Track 2 libraries (Open
%   Babel, MDAnalysis, PyMOL-OSS) that require a full CPython installation
%   and cannot be added to the Embedded Python environment.
%
%   If Python is already loaded in this MATLAB session, the call is silently
%   ignored with a warning (pyenv cannot be reconfigured once loaded).
%
%   After calling useExternal(), use emk.setup.validate() to confirm the
%   desired libraries are present in the external environment.
%
%   Arguments:
%     pythonPath  (string | char) - absolute path to the Python executable
%                                   (e.g., "C:\envs\myenv\python.exe")
%
%   Error IDs:
%     emk:setup:useExternal:invalidInput  - pythonPath is not a string/char
%                                           or is empty
%     emk:setup:useExternal:fileNotFound  - executable not found at path
%     emk:setup:useExternal:pyenvFailed   - pyenv() call threw an error
%
%   See also: emk.setup.validate, emk.setup.recipe, emk.setup.initPython

    % --- Input validation (manual; handles string and char) ---
    if ~(ischar(pythonPath) || isStringScalar(pythonPath))
        error("emk:setup:useExternal:invalidInput", ...
            "pythonPath must be a string or char scalar, got %s.", class(pythonPath));
    end

    pythonPath = string(pythonPath);

    if strlength(strtrim(pythonPath)) == 0
        error("emk:setup:useExternal:invalidInput", ...
            "pythonPath must be a non-empty path string.");
    end

    % --- Idempotency guard ---
    % pyenv(Version=...) cannot be changed once Python is Loaded or Terminating.
    pe = pyenv();
    if ~strcmp(string(pe.Status), "NotLoaded")
        logWarn("useExternal: Python %s already active (Status: %s). " + ...
            "Cannot change pyenv after Python is loaded -- ignoring.", ...
            string(pe.Version), string(pe.Status));
        return;
    end

    % --- Verify executable exists ---
    if ~isfile(pythonPath)
        error("emk:setup:useExternal:fileNotFound", ...
            "Python executable not found: %s\n" + ...
            "Ensure the path points to a valid Python installation.", pythonPath);
    end

    % --- Configure pyenv ---
    cfg      = loadConfig();
    execMode = cfg.python.execution_mode;

    logInfo("useExternal: Connecting to external Python: %s", pythonPath);

    try
        pe = pyenv(Version=pythonPath, ExecutionMode=execMode);
        logInfo("useExternal: Python %s configured (mode: %s).", ...
            string(pe.Version), string(pe.ExecutionMode));
    catch ME
        error("emk:setup:useExternal:pyenvFailed", ...
            "pyenv configuration failed: %s", ME.message);
    end
end
