function initPython()
% initPython  Configure pyenv for EasyMolKit (OutOfProcess mode).
%
%   emk.setup.initPython()
%
%   Platform-aware Python environment initialization:
%     Desktop : uses repo-local python_env/python.exe (Embedded Python 3.10)
%     Online  : uses system Python (pre-installed in MATLAB Online)
%
%   pyenv(Version, ExecutionMode) can only be changed before Python is first
%   used in a MATLAB session (Status == "NotLoaded").  A second call when
%   Python is already active is silently ignored to allow safe repeated calls
%   from scripts.
%
%   Error IDs:
%     emk:setup:initPython:notInstalled  - python_env/ not found (Desktop)
%     emk:setup:initPython:pyenvFailed   - pyenv() call threw an error

    % --- Double-call guard: must be FIRST to guarantee idempotency ---
    % pyenv() with no args is always safe (does not start Python).
    % Any Status other than "NotLoaded" means Python is already configured
    % for this MATLAB session; further calls to pyenv(Version=...) would fail.
    pe = pyenv();
    if ~strcmp(string(pe.Status), "NotLoaded")
        logInfo("initPython: Python %s already active (Status: %s) -- skipping.", ...
            string(pe.Version), string(pe.Status));
        return;
    end

    cfg    = loadConfig();
    online = emk.util.isOnline();

    % --- Determine target Python executable ---
    if online
        pyPath = "";   % MATLAB Online: use system Python, no path needed
        logInfo("initPython: MATLAB Online detected -- using system Python");
    else
        % Desktop: resolve project root from this file's location (ADR-005).
        % src/+emk/+setup/initPython.m -> 4 fileparts -> project root
        projectRoot = fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))));

        % Track 2 override: if settings.json defines python.external_path,
        % delegate to useExternal() and return (ADR-007).
        if isfield(cfg, "python") && isfield(cfg.python, "external_path")
            extPath = string(cfg.python.external_path);
            if strlength(strtrim(extPath)) > 0
                % Resolve relative path against project root.
                if ~isAbsolutePath_(extPath)
                    extPath = fullfile(projectRoot, char(extPath));
                end
                if ~isfile(char(extPath))
                    logWarn("initPython: Track 2 external_path not found: %s", extPath);
                    logWarn("             Falling back to embedded Python.");
                else
                    logInfo("initPython: Track 2 external Python detected: %s", extPath);
                    emk.setup.useExternal(extPath);
                    return;
                end
            end
        end

        % Default: Embedded Python (Track 1).
        pyPath = fullfile(projectRoot, cfg.python.embedded_dir, "python.exe");
        if ~isfile(pyPath)
            error("emk:setup:initPython:notInstalled", ...
                "Embedded Python not found at: %s\n" + ...
                "Run emk.setup.install() first.", pyPath);
        end
        logInfo("initPython: Desktop mode -- embedded Python: %s", pyPath);
    end

    % --- Configure pyenv ---
    execMode = cfg.python.execution_mode;
    try
        % NOTE: pyPath is a string scalar, not a char array.
        % isempty("") returns false in MATLAB (1x1 string scalar is not empty).
        % Use strlength to correctly detect "no path specified" (Online case).
        if strlength(pyPath) == 0
            pe = pyenv(ExecutionMode=execMode);
        else
            pe = pyenv(Version=pyPath, ExecutionMode=execMode);
        end
        logInfo("initPython: Python %s configured (mode: %s)", ...
            string(pe.Version), string(pe.ExecutionMode));
    catch ME
        error("emk:setup:initPython:pyenvFailed", ...
            "pyenv configuration failed: %s", ME.message);
    end
end

% =========================================================================
% Private helpers
% =========================================================================

function tf = isAbsolutePath_(p)
% Return true when p is an absolute Windows or UNC path.
% Covers: "C:\...", "\\server\...", "/..." (forward-slash absolute).
    p = char(p);
    tf = (numel(p) >= 2 && p(2) == ':') || ...  % Drive-letter: C:\
         (numel(p) >= 2 && p(1) == '\' && p(2) == '\') || ...  % UNC: \\
         (numel(p) >= 1 && p(1) == '/');         % Unix-style /
end
