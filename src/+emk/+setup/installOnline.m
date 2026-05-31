function installOnline(options)
% installOnline  Install RDKit on MATLAB Online via pip (Online only).
%
%   emk.setup.installOnline()
%   emk.setup.installOnline(Config=cfg)
%
%   Bootstraps pip if necessary and installs rdkit-pypi at the pinned
%   version from config/settings.json.  Safe to call every session:
%   if the correct version is already installed the pip step is skipped.
%   See ADR-001 rev.3 for rationale.
%
%   When Config is supplied (recommended via main_rdkit.m), installOnline
%   also installs optional libraries listed in cfg.optionalLibraries.
%
%   Processing steps:
%     1. Online guard     - Online only; Desktop -> error
%     2. Version check    - skip pip install if rdkit-pypi version matches config
%     3. Bootstrap pip    - websave get-pip.py; run python get-pip.py --user
%     4. Install RDKit    - pip install rdkit-pypi==<version> --user
%     5. Configure pyenv  - call emk.setup.initPython()
%     6. Insert sys.path  - add site-packages to py.sys.path (required for --user)
%     7. Smoke verify     - call emk.setup.verify()
%     8. Optional libs    - pip install --user for cfg.optionalLibraries.<name>=true
%
%   File artifacts: none (all files installed to ~/.local/ on the remote host)
%
%   Arguments:
%     Config  struct - pre-loaded config struct from emkLoadConfig().
%                      When omitted, emkLoadConfig() is called internally.
%
%   Error IDs:
%     emk:setup:installOnline:notOnline          - called on Desktop
%     emk:setup:installOnline:pipBootstrapFailed - get-pip.py execution failed
%     emk:setup:installOnline:rdkitInstallFailed - pip install rdkit-pypi failed
%
%   See also: emk.setup.install, emk.setup.initPython, emk.setup.verify

    arguments
        options.Config struct = struct()
    end

    logInfo("installOnline: Starting RDKit setup for MATLAB Online...");

    % STEP 1: Guard -- Online only
    if ~emk.util.isOnline()
        error("emk:setup:installOnline:notOnline", ...
            "installOnline() is for MATLAB Online only. " + ...
            "Use emk.setup.install() on Desktop.");
    end

    % Use supplied Config or load from file/defaults
    if isempty(fieldnames(options.Config))
        cfg = emkLoadConfig();
    else
        cfg = options.Config;
    end

    % STEP 2: Version check -- skip pip install if already at correct version
    if isRdkitInstalled_(cfg)
        logInfo("installOnline: rdkit-pypi %s already installed -- skipping pip.", ...
            cfg.rdkit.version);
    else
        % STEP 3: Bootstrap pip via get-pip.py
        bootstrapPip_(cfg);

        % STEP 4: Install RDKit at pinned version
        installRdkit_(cfg);
    end

    % STEP 5: Configure pyenv (must precede py.sys.path access)
    emk.setup.initPython();

    % STEP 6: Insert site-packages into py.sys.path
    % Required: --user installed packages are not visible to MATLAB's pyenv
    % without explicit path insertion (MLChem.m prior art; ADR-001 rev.3).
    insertSysPath_(cfg);

    % STEP 7: Final smoke test
    result = emk.setup.verify();
    if ~result.python || ~result.rdkit
        logWarn("installOnline: Verify failed (python=%d, rdkit=%d).", ...
            result.python, result.rdkit);
    else
        logInfo("installOnline: RDKit setup for MATLAB Online complete.");
    end

    % STEP 8: Expand use-case groups into individual library flags, then install
    cfg = emk.setup.expandUseCases(cfg);
    if isfield(cfg, "optionalLibraries")
        libs = fieldnames(cfg.optionalLibraries);
        for k = 1:numel(libs)
            libName = libs{k};
            if islogical(cfg.optionalLibraries.(libName)) && ...
                    cfg.optionalLibraries.(libName)
                logInfo("installOnline: Installing optional library '%s' ...", libName);
                try
                    installOptionalLib_(libName, cfg);
                catch ME
                    logWarn("installOnline: Optional library '%s' failed: %s", ...
                        libName, ME.message);
                end
            end
        end
    end
end

% =========================================================================
% Private helpers
% =========================================================================

function tf = isRdkitInstalled_(cfg)
% Return true if pip is available and rdkit is installed at the
% exact version specified in cfg.rdkit.version (checked via pip show).
    pipPath = getPipPath_();
    if ~isfile(pipPath)
        tf = false;
        return;
    end

    cmd = sprintf('"%s" show rdkit 2>/dev/null', pipPath);
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
    % Normalize both versions before comparison: pip (PEP 440) drops leading
    % zeros (e.g. "2024.03.6" -> "2024.3.6"), so a plain strcmp would always
    % mismatch and trigger an unnecessary reinstall.
    tf = strcmp(normalizeVer_(installed), normalizeVer_(cfg.rdkit.version));
    if ~tf
        logInfo("installOnline: Found rdkit %s but expected %s -- will reinstall.", ...
            installed, cfg.rdkit.version);
    end
end

% -------------------------------------------------------------------------
function bootstrapPip_(cfg)
% Download get-pip.py and bootstrap pip with --user install.
    getPipUrl  = "https://bootstrap.pypa.io/get-pip.py";
    getPipPath = fullfile(tempdir, "get-pip.py");

    logInfo("installOnline: Downloading get-pip.py ...");
    try
        websave(getPipPath, getPipUrl, buildWebOptions_(cfg));
    catch ME
        error("emk:setup:installOnline:pipBootstrapFailed", ...
            "Failed to download get-pip.py: %s", ME.message);
    end

    logInfo("installOnline: Bootstrapping pip ...");
    cmd = sprintf('python "%s" --user 2>&1', getPipPath);
    [status, output] = system(cmd);

    if isfile(getPipPath)
        delete(getPipPath);
    end

    if status ~= 0
        error("emk:setup:installOnline:pipBootstrapFailed", ...
            "get-pip.py bootstrap failed (exit code %d):\n%s", status, output);
    end

    logInfo("installOnline: pip bootstrap complete.");
end

% -------------------------------------------------------------------------
function installRdkit_(cfg)
% Run pip install rdkit at the pinned version (ADR-001 rev.3).
    pipPath  = getPipPath_();
    rdkitSpec = sprintf("rdkit==%s", cfg.rdkit.version);

    logInfo("installOnline: Installing %s (this may take a few minutes) ...", rdkitSpec);

    cmd = sprintf('"%s" install %s --user%s 2>&1', ...
        pipPath, rdkitSpec, buildPipProxyArg_(cfg));
    [status, output] = system(cmd);

    if status ~= 0
        error("emk:setup:installOnline:rdkitInstallFailed", ...
            "pip install %s failed (exit code %d):\n%s", rdkitSpec, status, output);
    end

    logInfo("installOnline: %s installed.", rdkitSpec);
end

% -------------------------------------------------------------------------
function insertSysPath_(cfg)
% Insert the --user site-packages directory into py.sys.path at index 0.
% Required because MATLAB's pyenv does not automatically pick up packages
% installed with "pip install --user" (ADR-001 rev.3, MLChem.m prior art).
    homeDir = getenv("HOME");
    if isempty(homeDir)
        homeDir = "/home/matlab";   % MATLAB Online default fallback
        logWarn("installOnline: HOME env var not set; using fallback: %s", homeDir);
    end

    pyVer        = cfg.python.version;   % e.g., "3.10"
    sitePackages = fullfile(homeDir, ".local", "lib", ...
        sprintf("python%s", pyVer), "site-packages");

    logInfo("installOnline: Inserting site-packages path: %s", sitePackages);
    try
        py.sys.path().insert(int32(0), sitePackages);
    catch ME
        logWarn("installOnline: py.sys.path insert failed: %s", ME.message);
    end
end

% -------------------------------------------------------------------------
function p = getPipPath_()
% Return the expected pip executable path for --user installs on MATLAB Online.
    homeDir = getenv("HOME");
    if isempty(homeDir)
        homeDir = "/home/matlab";
    end
    p = fullfile(homeDir, ".local", "bin", "pip");
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
function v = normalizeVer_(raw)
% Normalize a PEP-440 version string by removing leading zeros from each
% numeric segment so that "2024.03.6" and "2024.3.6" compare as equal.
    parts = strsplit(strtrim(string(raw)), '.');
    normed = arrayfun(@(p) string(str2double(p)), parts);
    v = strjoin(normed, '.');
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

% -------------------------------------------------------------------------
function installOptionalLib_(name, cfg)
% Install a Track 1 optional library via pip --user on MATLAB Online (Linux).
% Uses pip from ~/.local/bin/pip (bootstrapped by bootstrapPip_).
    pipName = resolveOnlinePipName_(name);
    if strlength(pipName) == 0
        logWarn("installOnline: Unknown optional library '%s' -- skipping.", name);
        return;
    end

    pipPath   = getPipPath_();
    proxyArg  = buildPipProxyArg_(cfg);
    cmd = sprintf('"%s" install %s --user%s 2>&1', pipPath, pipName, proxyArg);
    [status, output] = system(cmd);

    if status ~= 0
        error("emk:setup:installOnline:optionalInstallFailed", ...
            "pip install %s failed (exit code %d):\n%s", pipName, status, output);
    end

    logInfo("installOnline: %s installed.", name);
end

% -------------------------------------------------------------------------
function pipName = resolveOnlinePipName_(name)
% Map user-facing library name to pip package specifier for MATLAB Online (Linux).
% Returns empty string for unrecognised names (caller skips gracefully).
    switch name
        case "pubchempy",       pipName = "pubchempy";
        case "mordred",         pipName = "mordredcommunity";
        case "biopython",       pipName = "biopython";
        case "meeko",           pipName = "meeko gemmi";
        case "vina",            pipName = "vina";
        case "pdbfixer",        pipName = "pdbfixer";
        case "scipy",           pipName = "scipy";
        case "torch",           pipName = "torch --index-url https://download.pytorch.org/whl/cpu";
        case "torch_geometric", pipName = "torch_geometric";
        case "transformers",    pipName = "transformers";
        case "datasets",        pipName = "datasets";
        otherwise,              pipName = "";
    end
end
