function installTrack2(name, options)
% installTrack2  Install a Track 2 library into a dedicated venv (Desktop only).
%
%   emk.setup.installTrack2(name)
%   emk.setup.installTrack2(name, BasePython="C:\Python310\python.exe")
%
%   Creates a repo-local virtual environment at python_env_t2/<name>/,
%   installs the requested GPL library into it via pip, writes the venv
%   python.exe path to config/settings.json (python.external_path), and
%   connects EasyMolKit to the new environment via useExternal().
%
%   On subsequent MATLAB sessions, emk.setup.initPython() automatically
%   detects python.external_path and calls useExternal() without any extra
%   user action.
%
%   Supported library names:
%     "mdanalysis"  - MD trajectory analysis (GPLv2+, ~50 MB)
%     "pymol"       - 3D molecular visualization (PSF/BSD, ~200 MB)
%
%   Open Babel is NOT supported here because it requires a Windows MSI
%   installer.  Use emk.setup.recipe("openbabel") for manual instructions.
%
%   Arguments:
%     name          (1,1) string - library name ("mdanalysis" or "pymol")
%     BasePython    (1,1) string - absolute path to a CPython 3.10+ executable
%                                  used to create the venv.  When omitted,
%                                  auto-detected from PATH ("python" or "py").
%     Config        struct       - pre-loaded config from loadConfig().
%
%   Error IDs:
%     emk:setup:installTrack2:notDesktop         - called on MATLAB Online
%     emk:setup:installTrack2:unknownLibrary     - unsupported library name
%     emk:setup:installTrack2:basePythonNotFound - base Python not found
%     emk:setup:installTrack2:venvFailed         - venv creation failed
%     emk:setup:installTrack2:installFailed      - pip install failed
%     emk:setup:installTrack2:importVerifyFailed - post-install import check failed
%     emk:setup:installTrack2:settingsWriteFailed - settings.json update failed
%
%   See also: emk.setup.recipe, emk.setup.useExternal, emk.setup.validate

    arguments
        name               (1,1) string
        options.BasePython (1,1) string = ""
        options.Config     struct        = struct()
    end

    % STEP 1: Validate library name (before platform check for precise error).
    info = resolveLibInfo_(name);

    % STEP 2: Desktop-only guard.
    if emk.util.isOnline()
        error("emk:setup:installTrack2:notDesktop", ...
            "installTrack2() is only available on MATLAB Desktop.\n" + ...
            "On MATLAB Online, install packages manually:\n" + ...
            "  !pip install --user %s", info.pypiName);
    end

    % STEP 3: Load config.
    if isempty(fieldnames(options.Config))
        cfg = loadConfig();
    else
        cfg = options.Config;
    end

    % STEP 4: Resolve base Python executable.
    basePy = resolveBasePython_(options.BasePython);

    % STEP 5: Build venv path (absolute, mfilename-based per ADR-005/ADR-007).
    projectRoot = fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))));
    venvDir     = fullfile(projectRoot, "python_env_t2", char(name));
    venvPy      = fullfile(venvDir, "Scripts", "python.exe");
    pipExe      = fullfile(venvDir, "Scripts", "pip.exe");

    % STEP 6: Show GPL warning.
    logWarn("installTrack2: %s is licensed under %s.", name, info.license);
    logWarn("installTrack2: EasyMolKit (MIT) is not affected, but your scripts");
    logWarn("installTrack2: that import %s may be subject to %s terms.", ...
        name, info.license);

    % STEP 7: Create venv (idempotent -- skip if already exists).
    if isfile(venvPy)
        logInfo("installTrack2: venv already exists at %s -- skipping creation.", venvDir);
    else
        logInfo("installTrack2: Creating venv at %s ...", venvDir);
        cmd = ['"', char(basePy), '" -m venv "', char(venvDir), '"'];
        [status, output] = system(cmd);
        if status ~= 0
            error("emk:setup:installTrack2:venvFailed", ...
                "venv creation failed (exit code %d):\n%s\n" + ...
                "Base Python: %s\nTarget: %s", ...
                status, output, basePy, venvDir);
        end
        logInfo("installTrack2: venv created.");
    end

    % STEP 8: pip install each package.
    proxyArg = buildPipProxyArg_(cfg);
    if info.useGohlkeWheel
        % pymol: use Gohlke prebuilt wheels which bundle all native DLLs.
        % pymol-open-source on PyPI causes DLL load failures on Windows.
        wheelUrl = buildGohlkeWheelUrl_(venvPy);
        logInfo("installTrack2: Installing pymol from Gohlke prebuilt wheel ...");
        logInfo("installTrack2:   %s", wheelUrl);
        cmd = ['"', char(pipExe), '" install "', char(wheelUrl), '"' ...
            ' --no-warn-script-location', char(proxyArg)];
        [status, output] = system(cmd);
        if status ~= 0
            error("emk:setup:installTrack2:installFailed", ...
                "pip install pymol (Gohlke wheel) failed (exit code %d):\n%s", ...
                status, output);
        end
        logInfo("installTrack2: pymol installed.");
    else
        for i = 1:numel(info.pipPackages)
            pkg = info.pipPackages(i);
            logInfo("installTrack2: Installing %s ...", pkg);
            cmd = ['"', char(pipExe), '" install ', char(pkg), ...
                ' --no-warn-script-location', char(proxyArg)];
            [status, output] = system(cmd);
            if status ~= 0
                error("emk:setup:installTrack2:installFailed", ...
                    "pip install %s failed (exit code %d):\n%s", pkg, status, output);
            end
            logInfo("installTrack2: %s installed.", pkg);
        end
    end

    % STEP 9: Verify import.
    verifyCmd = ['"', char(venvPy), '" -c "import ', char(info.importName), '" 2>&1'];
    [vstatus, vout] = system(verifyCmd);
    if vstatus ~= 0
        error("emk:setup:installTrack2:importVerifyFailed", ...
            "Import verification failed for '%s' (import %s):\n%s", ...
            name, info.importName, vout);
    end
    logInfo("installTrack2: import %s OK.", info.importName);

    % STEP 10: Write external_path to config/settings.json.
    settingsPath = fullfile(projectRoot, "config", "settings.json");
    writeExternalPath_(settingsPath, venvPy);

    % STEP 11: Connect to the new venv.
    emk.setup.useExternal(venvPy);

    logInfo("installTrack2: %s ready.", name);
    logInfo("installTrack2: On next session, emk.setup.initPython() auto-detects.");
end

% =========================================================================
% Private helpers
% =========================================================================

function info = resolveLibInfo_(name)
% Map library name to pip packages, import name, and license string.
    switch name
        case "mdanalysis"
            info = struct( ...
                pipPackages    = ["MDAnalysis"], ...
                pypiName       = "MDAnalysis", ...
                importName     = "MDAnalysis", ...
                license        = "GPLv2+", ...
                useGohlkeWheel = false);
        case "pymol"
            % pymol-open-source from PyPI has DLL load failures on Windows
            % because it does not bundle native DLLs (glew, freeglut, etc.).
            % Use Christoph Gohlke's prebuilt wheels instead -- they include
            % all required DLLs and work on Windows without extra system packages.
            % Reference: https://github.com/cgohlke/pymol-open-source-wheels
            info = struct( ...
                pipPackages    = string.empty(1,0), ...
                pypiName       = "pymol-open-source", ...
                importName     = "pymol", ...
                license        = "PSF/BSD hybrid", ...
                useGohlkeWheel = true);
        otherwise
            error("emk:setup:installTrack2:unknownLibrary", ...
                "Unknown Track 2 library: '%s'.\n" + ...
                "Supported names: mdanalysis, pymol.\n" + ...
                "For Open Babel (MSI required), see emk.setup.recipe(""openbabel"").", ...
                name);
    end
end

% -------------------------------------------------------------------------
function basePy = resolveBasePython_(hint)
% Return the absolute path to a usable base CPython executable.
% Priority: explicit hint > "py" launcher > "python" in PATH.
    if strlength(strtrim(hint)) > 0
        basePy = hint;
        if ~isfile(basePy)
            error("emk:setup:installTrack2:basePythonNotFound", ...
                "Specified BasePython not found: %s\n" + ...
                "Provide the absolute path to a CPython 3.10+ executable.", basePy);
        end
        return;
    end

    % Try Windows Python Launcher ("py") first (most reliable on Windows).
    [st, out] = system("py -3 --version 2>&1");
    if st == 0
        [st2, pyPath] = system("py -3 -c ""import sys; print(sys.executable)"" 2>&1");
        if st2 == 0
            basePy = string(strtrim(pyPath));
            return;
        end
    end

    % Fall back to "python" in PATH.
    [st, out] = system("python --version 2>&1");
    if st == 0
        [st2, pyPath] = system("python -c ""import sys; print(sys.executable)"" 2>&1");
        if st2 == 0
            basePy = string(strtrim(pyPath));
            return;
        end
    end

    error("emk:setup:installTrack2:basePythonNotFound", ...
        "Could not find a CPython executable in PATH.\n" + ...
        "Install Python 3.10+ from https://www.python.org/downloads/ and\n" + ...
        "pass the path explicitly:\n" + ...
        "  emk.setup.installTrack2(""%s"", BasePython=""C:\\Python310\\python.exe"")", ...
        out);
end

% -------------------------------------------------------------------------
function url = buildGohlkeWheelUrl_(venvPy)
% Return the Gohlke prebuilt wheel URL for PyMOL matching the venv Python version.
% Gohlke wheels bundle all native DLLs (glew, freeglut, freetype, etc.) and avoid
% the DLL load failure that occurs with pymol-open-source from PyPI on Windows.
    [~, verOut] = system(['"', char(venvPy), '" --version 2>&1']);
    tok = regexp(strtrim(verOut), 'Python (\d+)\.(\d+)', 'tokens', 'once');
    if isempty(tok)
        error('emk:setup:installTrack2:gohlkeWheelNotFound', ...
            'Could not detect Python version for Gohlke wheel selection.\nOutput: %s', ...
            strtrim(verOut));
    end
    major = str2double(tok{1});
    minor = str2double(tok{2});
    base  = 'https://github.com/cgohlke/pymol-open-source-wheels/releases/download/';
    switch sprintf('%d.%d', major, minor)
        case '3.10'
            url = string(base) + 'v2025.2.2/pymol-3.1.0-cp310-cp310-win_amd64.whl';
        case '3.11'
            url = string(base) + 'v2026.4.5/pymol-3.2.0a0-cp311-cp311-win_amd64.whl';
        case '3.12'
            url = string(base) + 'v2026.4.5/pymol-3.2.0a0-cp312-cp312-win_amd64.whl';
        case '3.13'
            url = string(base) + 'v2026.4.5/pymol-3.2.0a0-cp313-cp313-win_amd64.whl';
        otherwise
            error('emk:setup:installTrack2:gohlkeWheelNotFound', ...
                'No Gohlke PyMOL wheel for Python %d.%d.\nSupported: 3.10, 3.11, 3.12, 3.13.\n' + ...
                'See: https://github.com/cgohlke/pymol-open-source-wheels', ...
                major, minor);
    end
end

% -------------------------------------------------------------------------
function arg = buildPipProxyArg_(cfg)
% Return ' --proxy "<url>"' for pip, or "" if no proxy is configured.
    if isfield(cfg, "python") && isfield(cfg.python, "proxy") && ...
            ~isempty(cfg.python.proxy) && cfg.python.proxy ~= ""
        arg = sprintf(' --proxy "%s"', cfg.python.proxy);
    else
        arg = "";
    end
end

% -------------------------------------------------------------------------
function writeExternalPath_(settingsPath, venvPy)
% Persist venvPy as python.external_path in config/settings.json.
% Creates the file if absent; updates the field if present.
    try
        if isfile(settingsPath)
            raw = jsondecode(fileread(settingsPath));
        else
            raw = struct();
        end

        if ~isfield(raw, "python")
            raw.python = struct();
        end
        raw.python.external_path = char(venvPy);

        fid = fopen(settingsPath, "w");
        if fid == -1
            error("emk:setup:installTrack2:settingsWriteFailed", ...
                "Cannot open %s for writing.", settingsPath);
        end
        fprintf(fid, "%s\n", jsonencode(raw, PrettyPrint=true));
        fclose(fid);
        logInfo("installTrack2: external_path written to %s.", settingsPath);
    catch ME
        if strcmp(ME.identifier, "emk:setup:installTrack2:settingsWriteFailed")
            rethrow(ME);
        end
        error("emk:setup:installTrack2:settingsWriteFailed", ...
            "Failed to write external_path to settings.json: %s", ME.message);
    end
end
