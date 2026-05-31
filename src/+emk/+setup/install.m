function install(options)
% install  Deploy Embedded Python 3.10 + RDKit to python_env/ (Desktop only).
%
%   emk.setup.install()
%   emk.setup.install(Config=cfg)
%   emk.setup.install(PythonVersion="3.10")
%
%   Installs a repo-local, fully isolated Python environment.
%   Does NOT modify system PATH, registry, or any existing Python install.
%   Safe to call multiple times (idempotent).
%   See ADR-001 rev.3 for rationale and design decisions.
%
%   When Config is supplied (recommended), install() uses that struct instead
%   of calling emkLoadConfig() internally.  This lets callers set options inline
%   in main*.m without requiring a settings.json file.
%   cfg.optionalLibraries.<name>=true entries are installed after RDKit.
%
%   Processing steps:
%     1. Online guard        - Desktop only; MATLAB Online -> error
%     2. Path length check   - warn >200 chars, error >240 chars (ADR-001 rev.3)
%     3. Idempotency check   - skip if python_env/ already has correct install
%     4. Download Python zip - python-3.10.x-embed-amd64.zip from python.org
%     5. Extract             - unzip into python_env/
%     6. Enable site-pkgs    - edit python31x._pth: uncomment "import site"
%     7. Bootstrap pip       - download get-pip.py, run with python.exe, clean up
%     8. Install RDKit       - pip install rdkit-pypi==<cfg.rdkit.version>
%     9. Version check       - pip show rdkit-pypi; warn on mismatch
%    10. Configure pyenv     - call emk.setup.initPython()
%    11. Smoke verify        - call emk.setup.verify()
%    12. Optional libraries  - installExtra() for cfg.optionalLibraries.<name>=true
%
%   Arguments:
%     Config        struct  - pre-loaded config struct from emkLoadConfig().
%                             When omitted, emkLoadConfig() is called internally.
%     PythonVersion string  - Python version to download (default: "3.10").
%
%   File artifacts (all under python_env/):
%     python.exe                       - Python executable
%     python31x._pth                   - edited to enable site-packages
%     Lib/site-packages/rdkit/         - RDKit package files
%
%   Error IDs:
%     emk:setup:install:notDesktop         - called on MATLAB Online
%     emk:setup:install:pathTooLong        - install path > 240 chars
%     emk:setup:install:downloadFailed     - network download failed
%     emk:setup:install:extractFailed      - zip extraction or _pth edit failed
%     emk:setup:install:pipBootstrapFailed - get-pip.py execution failed
%     emk:setup:install:rdkitInstallFailed - pip install rdkit-pypi failed
%
%   See also: emk.setup.installOnline, emk.setup.initPython, emk.setup.verify

    arguments
        options.Config        struct = struct()
        options.PythonVersion (1,1) string = "3.10"
    end

    logInfo("install: Starting Embedded Python + RDKit setup...");

    % STEP 1: Guard -- Desktop only
    if emk.util.isOnline()
        error("emk:setup:install:notDesktop", ...
            "install() is for Desktop only. " + ...
            "Use emk.setup.installOnline() on MATLAB Online.");
    end

    % Resolve config: use supplied Config or load from file/defaults
    if isempty(fieldnames(options.Config))
        cfg = emkLoadConfig();
    else
        cfg = options.Config;
    end

    % STEP 2: Path length check (ADR-001 rev.3, Risk #10)
    checkPathLength_(cfg);

    % STEP 3: Idempotency check
    alreadyInstalled = isInstalled_(cfg);
    if alreadyInstalled
        logInfo("install: python_env/ already installed at correct version -- skipping download.");
        emk.setup.initPython();
        emk.setup.verify();
        % Fall through to STEP 12 to install any newly enabled optional libraries.
    end

    if ~alreadyInstalled
        % STEP 4: Download Python embeddable zip
        zipPath = downloadPython_(options.PythonVersion, cfg);

        % STEP 5: Extract to python_env/
        extractPython_(zipPath, cfg);

        % STEP 6: Edit _pth file to enable site-packages (required for pip)
        enableImportSite_(cfg);

        % STEP 7: Bootstrap pip via get-pip.py
        bootstrapPip_(cfg);

        % STEP 8: Install RDKit at pinned version
        installRdkit_(cfg);

        % STEP 9: Verify installed version matches config
        verifyRdkitVersion_(cfg);

        % STEP 10: Configure pyenv
        emk.setup.initPython();

        % STEP 11: Final smoke test
        result = emk.setup.verify();
        if ~result.python || ~result.rdkit
            logWarn("install: Verify failed after setup (python=%d, rdkit=%d).", ...
                result.python, result.rdkit);
        else
            logInfo("install: Embedded Python + RDKit setup complete.");
        end
    end

    % STEP 12: Expand use-case groups into individual library flags, then install
    cfg = emk.setup.expandUseCases(cfg);
    if isfield(cfg, "optionalLibraries")
        libs = fieldnames(cfg.optionalLibraries);
        for k = 1:numel(libs)
            libName = libs{k};
            if islogical(cfg.optionalLibraries.(libName)) && ...
                    cfg.optionalLibraries.(libName)
                logInfo("install: Installing optional library '%s' ...", libName);
                try
                    emk.setup.installExtra(libName, Config=cfg);
                catch ME
                    logWarn("install: Optional library '%s' failed: %s", ...
                        libName, ME.message);
                end
            end
        end
    end
end

% =========================================================================
% Private helpers
% =========================================================================

function checkPathLength_(cfg)
% Check Windows MAX_PATH limit for the deepest install path (ADR-001 rev.3).
% Warn if > 200 chars; hard error if > 240 chars.
    targetPath = fullfile(pwd, cfg.python.embedded_dir, ...
        "Lib", "site-packages", "rdkit");
    n = strlength(targetPath);
    if n > 240
        error("emk:setup:install:pathTooLong", ...
            "Install path length (%d chars) exceeds 240-char limit (MAX_PATH risk).\n" + ...
            "Move the project to a shorter directory path.\n  Path: %s", n, targetPath);
    elseif n > 200
        logWarn("install: Path length %d > 200 chars  EWindows MAX_PATH issues possible.", n);
    end
end

% -------------------------------------------------------------------------
function tf = isInstalled_(cfg)
% Return true if python.exe exists AND rdkit is installed at the
% exact version specified in cfg.rdkit.version (checked via pip show).
    pyExe = fullfile(pwd, cfg.python.embedded_dir, "python.exe");
    if ~isfile(pyExe)
        tf = false;
        return;
    end

    cmd = sprintf('"%s" -m pip show rdkit 2>nul', pyExe);
    [status, output] = system(cmd);
    if status ~= 0
        tf = false;
        return;
    end

    tok = regexp(output, 'Version:\s*([\d.]+)', 'tokens', 'once');
    if isempty(tok)
        tf = false;
        return;
    end

    installed = tok{1};
    tf = strcmp(normalizeVer_(installed), normalizeVer_(cfg.rdkit.version));
    if ~tf
        logInfo("install: Found rdkit %s but expected %s  -- will reinstall.", ...
            installed, cfg.rdkit.version);
    end
end

% -------------------------------------------------------------------------
function zipPath = downloadPython_(version, cfg)
% Download Python embeddable zip from python.org to tempdir.
% Returns the local path to the downloaded zip.
%
% Python 3.10 patch version is pinned for reproducibility (ADR-001 rev.3).
% Update patchVersion when a new 3.10.x is needed.
    patchVersion = resolvePatchVersion_(version);
    zipName = sprintf("python-%s-embed-amd64.zip", patchVersion);
    baseUrl = sprintf("https://www.python.org/ftp/python/%s/%s", patchVersion, zipName);
    zipPath = fullfile(tempdir, zipName);

    logInfo("install: Downloading Python %s ...", patchVersion);

    try
        websave(zipPath, baseUrl, buildWebOptions_(cfg));
    catch ME
        error("emk:setup:install:downloadFailed", ...
            "Failed to download Python %s: %s\n  URL: %s", patchVersion, ME.message, baseUrl);
    end

    logInfo("install: Python zip downloaded to %s", zipPath);
end

% -------------------------------------------------------------------------
function patchVer = resolvePatchVersion_(majorMinor)
% Map major.minor string to the pinned patch version.
% Only "3.10" is supported (ADR-001 rev.3: Python 3.10 fixed).
    switch majorMinor
        case "3.10"
            patchVer = "3.10.11";
        otherwise
            error("emk:setup:install:downloadFailed", ...
                "Unsupported Python version: '%s'. Only '3.10' is supported.", majorMinor);
    end
end

% -------------------------------------------------------------------------
function extractPython_(zipPath, cfg)
% Extract Python embeddable zip into python_env/.
    destDir = fullfile(pwd, cfg.python.embedded_dir);

    logInfo("install: Extracting Python to %s ...", destDir);
    try
        unzip(zipPath, destDir);
    catch ME
        error("emk:setup:install:extractFailed", ...
            "Failed to extract Python zip: %s", ME.message);
    end

    if isfile(zipPath)
        delete(zipPath);
    end

    logInfo("install: Extraction complete.");
end

% -------------------------------------------------------------------------
function enableImportSite_(cfg)
% Edit python31x._pth to uncomment "#import site".
% Embedded Python disables site-packages by default; this one edit enables
% pip and installed packages (including rdkit).
    embDir   = fullfile(pwd, cfg.python.embedded_dir);
    pthFiles = dir(fullfile(embDir, "python3*._pth"));

    if isempty(pthFiles)
        error("emk:setup:install:extractFailed", ...
            "No python3xx._pth file found in %s after extraction.", embDir);
    end

    pthPath = fullfile(embDir, pthFiles(1).name);
    raw     = fileread(pthPath);

    % Uncomment "#import site" -> "import site"
    updated = regexprep(raw, "#import site", "import site");

    if strcmp(raw, updated)
        logWarn("install: '#import site' not found in %s  Emay already be enabled.", ...
            pthFiles(1).name);
    end

    fid = fopen(pthPath, "w");
    if fid == -1
        error("emk:setup:install:extractFailed", ...
            "Cannot open %s for writing.", pthPath);
    end
    fwrite(fid, updated, "char");
    fclose(fid);

    logInfo("install: Site-packages enabled in %s", pthFiles(1).name);
end

% -------------------------------------------------------------------------
function bootstrapPip_(cfg)
% Download get-pip.py, run it with embedded python.exe, then clean up.
    getPipUrl  = "https://bootstrap.pypa.io/get-pip.py";
    getPipPath = fullfile(tempdir, "get-pip.py");
    pyExe      = fullfile(pwd, cfg.python.embedded_dir, "python.exe");

    logInfo("install: Downloading get-pip.py ...");
    try
        websave(getPipPath, getPipUrl, buildWebOptions_(cfg));
    catch ME
        error("emk:setup:install:downloadFailed", ...
            "Failed to download get-pip.py: %s", ME.message);
    end

    logInfo("install: Bootstrapping pip ...");
    cmd = ['"', char(pyExe), '" "', char(getPipPath), '"', char(buildPipProxyArg_(cfg))];
    [status, output] = system(cmd);

    if isfile(getPipPath)
        delete(getPipPath);
    end

    if status ~= 0
        if contains(output, "403") || contains(output, "Forbidden")
            error("emk:setup:install:pipBootstrapFailed", ...
                ["pip bootstrap failed: PyPI access blocked (HTTP 403).\n" ...
                 "Your organization network policy may be blocking files.pythonhosted.org.\n" ...
                 "Alternative: use MATLAB Online and run emk.setup.installOnline()\n\n" ...
                 "Details:\n%s"], output);
        end
        error("emk:setup:install:pipBootstrapFailed", ...
            "pip bootstrap failed (exit code %d):\n%s", status, output);
    end

    logInfo("install: pip bootstrap complete.");
end

% -------------------------------------------------------------------------
function installRdkit_(cfg)
% Install rdkit at the version pinned in cfg.rdkit.version.
    pyExe  = fullfile(pwd, cfg.python.embedded_dir, "python.exe");
    rdkVer = cfg.rdkit.version;

    logInfo("install: Installing rdkit==%s (this may take a few minutes) ...", rdkVer);

    cmd = ['"', char(pyExe), '" -m pip install rdkit==', ...
        char(rdkVer), ' --no-warn-script-location', char(buildPipProxyArg_(cfg))];
    [status, output] = system(cmd);

    if status ~= 0
        if contains(output, "403") || contains(output, "Forbidden")
            error("emk:setup:install:rdkitInstallFailed", ...
                ["pip install rdkit failed: PyPI access blocked (HTTP 403).\n" ...
                 "Your organization network policy may be blocking files.pythonhosted.org.\n" ...
                 "Alternative: use MATLAB Online and run emk.setup.installOnline()\n\n" ...
                 "Details:\n%s"], output);
        end
        error("emk:setup:install:rdkitInstallFailed", ...
            "pip install rdkit==%s failed (exit code %d):\n%s", rdkVer, status, output);
    end

    logInfo("install: rdkit==%s installed.", rdkVer);
end

% -------------------------------------------------------------------------
function verifyRdkitVersion_(cfg)
% Run "pip show rdkit" and warn if installed version != cfg.rdkit.version.
    pyExe    = fullfile(pwd, cfg.python.embedded_dir, "python.exe");
    expected = cfg.rdkit.version;

    cmd = sprintf('"%s" -m pip show rdkit', pyExe);
    [~, output] = system(cmd);

    tok = regexp(output, 'Version:\s*([\d.]+)', 'tokens', 'once');
    if isempty(tok)
        logWarn("install: Could not parse rdkit version from pip show output.");
        return;
    end

    installed = tok{1};
    if ~strcmp(normalizeVer_(installed), normalizeVer_(expected))
        logWarn("install: Installed rdkit '%s' differs from expected '%s'.", ...
            installed, expected);
    else
        logInfo("install: rdkit version verified: %s", installed);
    end
end

% -------------------------------------------------------------------------
function v = normalizeVer_(raw)
% Normalize a PEP-440 version string by removing leading zeros from each
% numeric segment so that "2024.03.6" and "2024.3.6" compare as equal.
    parts = strsplit(strtrim(string(raw)), '.');
    normed = arrayfun(@(p) string(str2double(p)), parts);
    v = strjoin(normed, '.');
end

% -------------------------------------------------------------------------
function wopts = buildWebOptions_(cfg)
% Build weboptions with optional proxy and a 60-second timeout.
    if ~isempty(cfg.python.proxy) && cfg.python.proxy ~= ""
        wopts = weboptions("ProxyServer", cfg.python.proxy, "Timeout", 60);
    else
        wopts = weboptions("Timeout", 60);
    end
end

% -------------------------------------------------------------------------
function arg = buildPipProxyArg_(cfg)
% Return ' --proxy "<url>"' string for pip, or "" if no proxy is configured.
    if ~isempty(cfg.python.proxy) && cfg.python.proxy ~= ""
        arg = sprintf(' --proxy "%s"', cfg.python.proxy);
    else
        arg = "";
    end
end
