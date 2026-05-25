function installExtra(name, options)
% installExtra  Install a Track 1 optional library into Embedded Python.
%
%   emk.setup.installExtra(name)
%   emk.setup.installExtra(name, Config=cfg)
%
%   Installs a Pure Python library into the repo-local Embedded Python
%   environment (python_env/) using pip.  Desktop only; MATLAB Online is
%   not supported because Embedded Python is not available there.
%
%   The pinned version is read from cfg.extraLibraries.<name> (when Config
%   is supplied) or from config/settings.json > extraLibraries.<name>.
%   If the field is absent, the latest available version is installed.
%
%   Supported library names:
%     "pubchempy"       - PubChem REST API client (MIT, ~1 MB)
%     "mordred"         - 1800+ molecular descriptors (mordredcommunity, BSD-3, ~5 MB)
%     "biopython"       - PDB / sequence analysis (Biopython License, ~50 MB)
%     "meeko"           - Ligand PDBQT preparation for AutoDock Vina (LGPL-2.1, ~1 MB)
%     "vina"            - AutoDock Vina docking engine Python API (Apache-2.0, ~30 MB)
%     "pdbfixer"        - PDB structure preparation; auto-installs openmm (MIT, ~70 MB)
%     "torch"           - PyTorch CPU-only (BSD-3, ~800 MB)  -- required for R08, R09
%     "torch_geometric" - PyG graph neural networks (MIT)    -- requires torch
%     "transformers"    - HuggingFace Transformers (Apache-2.0) -- for R09
%     "datasets"        - HuggingFace Datasets (Apache-2.0, ~50 MB) -- pairs with transformers
%
%   Arguments:
%     name    (1,1) string - library name (see above)
%     Config  struct       - pre-loaded config struct from loadConfig().
%                            When omitted, loadConfig() is called internally.
%
%   Error IDs:
%     emk:setup:installExtra:notDesktop         - called on MATLAB Online
%     emk:setup:installExtra:unknownLibrary     - unsupported name
%     emk:setup:installExtra:installFailed      - pip install failed
%     emk:setup:installExtra:importVerifyFailed - import check failed
%
%   See also: emk.setup.install, emk.setup.validate, emk.setup.recipe

    arguments
        name            (1,1) string
        options.Config  struct = struct()
    end

    % STEP 1: Validate library name (platform-independent; unknownLibrary before
    % any platform check so the caller gets the most precise error first).
    info = resolveLibInfo_(name);

    % STEP 2: Desktop-only guard (valid library names only reach here)
    if emk.util.isOnline()
        error("emk:setup:installExtra:notDesktop", ...
            "installExtra() requires Embedded Python (Desktop only).\n" + ...
            "On MATLAB Online, install packages manually:\n" + ...
            "  !~/.local/bin/pip install --user <package>");
    end

    % STEP 3: Locate Embedded Python executable
    if isempty(fieldnames(options.Config))
        cfg = loadConfig();   % no Config supplied -- load from file/defaults
    else
        cfg = options.Config; % use caller-supplied config (avoids re-load)
    end
    projectRoot = fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))));
    pyExe       = fullfile(projectRoot, cfg.python.embedded_dir, "python.exe");

    if ~isfile(pyExe)
        error("emk:setup:installExtra:installFailed", ...
            "Embedded Python not found at: %s\n" + ...
            "Run emk.setup.install() first.", pyExe);
    end

    % STEP 4: Build pip install specifier (with optional version pin)
    ver = resolveVersion_(cfg, name);
    if strlength(ver) > 0
        pipSpec = info.pipName + "==" + ver;
    else
        pipSpec = info.pipName;
    end

    % vina has no Windows wheels on PyPI -- redirect immediately rather than
    % attempting a futile source build that requires Boost C++.
    % Return normally (no error) so callers can continue installing other libraries.
    if name == "vina" && ispc()
        logWarn("installExtra: vina has no pre-built Windows wheels on PyPI.");
        logWarn("installExtra: Use MATLAB Online for docking workflows.");
        emk.setup.recipe("docking");
        return;
    end

    % pdbfixer requires openmm, whose _openmm.pyd is blocked by Windows Smart App
    % Control (SAC).  This is a Windows security policy that cannot be bypassed in
    % code.  Return normally so callers can continue with other libraries.
    if name == "pdbfixer" && ispc()
        logWarn("installExtra: pdbfixer requires openmm, which is blocked by Windows Smart App Control.");
        logWarn("installExtra: Use MATLAB Online for docking workflows.");
        emk.setup.recipe("docking");
        return;
    end

    if name == "torch"
        logWarn("installExtra: torch CPU-only is a large download (~800 MB). This may take several minutes.");
    elseif name == "pdbfixer"
        logWarn("installExtra: pdbfixer automatically installs openmm (~70 MB). This may take a few minutes.");
    elseif name == "meeko"
        % LGPL-2.1: dynamic import from EasyMolKit (MIT) is safe.
        % See docs/compliance.md CL-7 for full usage rationale.
        logWarn("installExtra: meeko is licensed under LGPL-2.1. " + ...
            "Dynamic import from EasyMolKit (MIT) is safe. " + ...
            "See docs/compliance.md CL-7 for details.");
    end
    logInfo("installExtra: Installing %s (%s) ...", name, pipSpec);

    % STEP 5: pip install (append pipSuffix for packages that need extra pip options,
    %   e.g. torch requires --index-url to select the CPU-only wheel)
    proxyArg = buildPipProxyArg_(cfg);

    % torch_geometric requires dynamic URL construction based on the installed
    % torch version and installs companion packages (torch_scatter, torch_sparse,
    % torch_cluster) in a single pip call.
    if name == "torch_geometric"
        installTorchGeometric_(pyExe, proxyArg);
        return;
    end

    % meeko requires gemmi (Mozilla Public License 2.0) which is not declared
    % as a pip dependency in some meeko versions.  Pre-install it explicitly.
    if name == "meeko"
        logInfo("installExtra: Pre-installing gemmi (required by meeko) ...");
        cmdGemmi = ['"', char(pyExe), '" -m pip install gemmi' ...
            ' --no-warn-script-location', char(proxyArg)];
        [gStatus, gOut] = system(cmdGemmi);
        if gStatus ~= 0
            error("emk:setup:installExtra:installFailed", ...
                "pip install gemmi (meeko pre-dep) failed (exit code %d):\n%s", ...
                gStatus, gOut);
        end
        logInfo("installExtra: gemmi installed successfully.");
    end

    pipSuffixStr = "";
    if strlength(info.pipSuffix) > 0
        pipSuffixStr = " " + info.pipSuffix;
    end
    cmd = ['"', char(pyExe), '" -m pip install ', char(pipSpec), ...
        ' --no-warn-script-location', char(proxyArg), char(pipSuffixStr)];
    [status, output] = system(cmd);

    if status ~= 0
        vinaHint = "";
        if name == "vina" && ispc()
            vinaHint = newline() + newline() + ...
                "HINT: AutoDock Vina Python bindings (vina) have no pre-built Windows " + newline() + ...
                "wheels on PyPI (Linux/macOS only). Options:" + newline() + ...
                "  1. Use the standalone Windows binary from GitHub Releases:" + newline() + ...
                "       https://github.com/ccsb-scripps/AutoDock-Vina/releases" + newline() + ...
                "  2. Install via conda-forge (requires Anaconda/Miniconda):" + newline() + ...
                "       conda install -c conda-forge vina" + newline() + ...
                "  3. Run docking on Linux/macOS or WSL2 where the wheel is available.";
        end
        error("emk:setup:installExtra:installFailed", ...
            "pip install %s failed (exit code %d):\n%s%s", pipSpec, status, output, vinaHint);
    end

    logInfo("installExtra: %s installed successfully.", name);

    % STEP 6: Unblock downloaded DLLs (Zone.Identifier / Mark of the Web)
    % pip-downloaded .pyd and .dll files may carry a Zone.Identifier alternate
    % data stream (set by browsers / some download managers).  Unblock-File
    % removes the zone marker without altering the file content.
    % Note: Smart App Control (SAC) reputation-based blocking is a separate
    % mechanism and requires disabling SAC in Windows Settings.
    unblockPydFiles_(pyExe);

    % Write/refresh sitecustomize.py so every Python session (including the
    % MATLAB pyenv process) auto-registers DLL search dirs on startup.
    % This is the permanent fix for the Python 3.8+ DLL search path issue.
    writeSiteCustomize_(pyExe);

    % STEP 7: Verify import
    % On Windows, pre-register site-packages subdirectories as DLL search paths
    % via a temp .py file (avoids cmd.exe quoting issues with -c "...").
    [vstatus, vout] = runImportCheck_(pyExe, info.importName);

    if vstatus ~= 0
        hint = "";
        if ispc() && contains(vout, char([12450 12503 12522 12465 12540 12471 12519 12531 21046 24481 12509 12522 12471 12540]))
            % char([...]) = Japanese: 'Application Control Policy' (Windows JA-JP SAC error)
            hint = newline() + newline() + ...
                "HINT: Smart App Control (SAC) is blocking an unsigned .pyd file." + newline() + ...
                "To fix: Settings > Windows Security > App & Browser Control" + newline() + ...
                "  > Smart App Control > Off  (permanent; cannot revert to Evaluation)." + newline() + ...
                "After turning SAC off, re-run emk.setup.installExtra(""" + name + """) " + ...
                "-- no reinstallation needed.";
        elseif ispc() && (contains(vout, "DLL load failed") || contains(vout, "_openmm"))
            hint = newline() + newline() + ...
                "HINT: A native DLL dependency could not be loaded." + newline() + ...
                "If the library uses OpenMM (e.g. pdbfixer), try:" + newline() + ...
                "  1. Install Microsoft Visual C++ Redistributable (2019 or later):" + newline() + ...
                "       https://aka.ms/vs/17/release/vc_redist.x64.exe" + newline() + ...
                "  2. Restart MATLAB and re-run emk.setup.installExtra(""" + name + """).";
        end
        error("emk:setup:installExtra:importVerifyFailed", ...
            "Import verification failed for '%s' (import %s):\n%s%s", ...
            name, info.importName, vout, hint);
    end

    logInfo("installExtra: import %s OK.", info.importName);
end

% =========================================================================
% Private helpers
% =========================================================================

function info = resolveLibInfo_(name)
% Map user-facing library name to pip package name, Python import name,
% and optional pip suffix arguments (e.g. --index-url for CPU-only torch).
    switch name
        case "pubchempy"
            info = struct(pipName="pubchempy",        importName="pubchempy",       pipSuffix="");
        case "mordred"
            info = struct(pipName="mordredcommunity",  importName="mordred",         pipSuffix="");
        case "biopython"
            info = struct(pipName="biopython",         importName="Bio",             pipSuffix="");
        case "torch"
            % CPU-only wheel -- avoids CUDA download (~800 MB vs ~2 GB with CUDA)
            info = struct(pipName="torch",             importName="torch",           ...
                pipSuffix="--index-url https://download.pytorch.org/whl/cpu");
        case "torch_geometric"
            % Requires torch to be installed first (handled by field order in cfg)
            info = struct(pipName="torch_geometric",   importName="torch_geometric", pipSuffix="");
        case "transformers"
            info = struct(pipName="transformers",       importName="transformers",    pipSuffix="");
        case "datasets"
            info = struct(pipName="datasets",            importName="datasets",        pipSuffix="");
        case "meeko"
            % LGPL-2.1: Python dynamic import does not trigger copyleft.
            % See docs/compliance.md CL-7 for usage rationale.
            info = struct(pipName="meeko",              importName="meeko",           pipSuffix="");
        case "gemmi"
            % gemmi is a pre-requisite of meeko (not declared as a pip dependency
            % in some versions of meeko).  Can also be installed standalone.
            info = struct(pipName="gemmi",              importName="gemmi",           pipSuffix="");
        case "vina"
            % vina has no Windows wheels on PyPI (Linux/macOS only).
            % On Windows, installExtra redirects to recipe("vina") before reaching here.
            % --only-binary :all: prevents a futile source build on Linux if wheel is absent.
            info = struct(pipName="vina",               importName="vina",            pipSuffix="--only-binary :all:");
        case "pdbfixer"
            % pdbfixer automatically installs openmm >= 8.2 as a dependency.
            % openmm: MIT + LGPL; dynamic import is safe for MIT projects.
            info = struct(pipName="pdbfixer",           importName="pdbfixer",        pipSuffix="");
        case "scipy"
            info = struct(pipName="scipy",               importName="scipy",           pipSuffix="");
        case "prody"
            % prody is required by meeko mk_prepare_receptor (--read_with_prody).
            info = struct(pipName="prody",               importName="prody",           pipSuffix="");
        otherwise
            error("emk:setup:installExtra:unknownLibrary", ...
                "Unknown library: '%s'.\n" + ...
                "Supported names: pubchempy, mordred, biopython, " + ...
                "scipy, meeko, gemmi, prody, vina, pdbfixer, " + ...
                "torch, torch_geometric, transformers, datasets.\n" + ...
                "For Track 2 / manual libraries (vina on Windows, openbabel, mdanalysis, pymol),\n" + ...
                "see emk.setup.recipe(name) for installation instructions.", name);
    end
end

% -------------------------------------------------------------------------
function unblockPydFiles_(pyExe)
% Remove Zone.Identifier (Mark of the Web) from .pyd and .dll files inside
% python_env/Lib/site-packages/ using PowerShell Unblock-File.
%
% NOTE: This only clears the Zone.Identifier NTFS alternate data stream
% (set when files are downloaded via browser or some download managers).
% It does NOT help with Smart App Control (SAC) reputation-based blocking,
% which evaluates Authenticode signatures and Microsoft cloud intelligence.
% If SAC is the cause, the user must disable it in Windows Settings.
%
% Implementation note:
%   Inline -Command "..." strings are passed through cmd.exe before PowerShell,
%   causing double-quote / brace quoting conflicts.  We write a .ps1 temp file
%   and execute with -File to avoid all inline quoting issues entirely.
    if ~ispc()
        return;
    end
    sitePackages = fullfile(fileparts(char(pyExe)), "Lib", "site-packages");
    if ~isfolder(sitePackages)
        logWarn("installExtra: site-packages not found at: %s -- skipping unblock.", sitePackages);
        return;
    end

    % Build PowerShell script as a MATLAB char array.
    % Single quotes inside char literals: '' -> '
    nl = newline();
    script = [ ...
        '$n = 0', nl, ...
        '$dir = ''', char(sitePackages), '''', nl, ...
        'Get-ChildItem -Recurse -LiteralPath $dir -ErrorAction SilentlyContinue |', nl, ...
        '    Where-Object { $_.Extension -eq ''.pyd'' -or $_.Extension -eq ''.dll'' } |', nl, ...
        '    ForEach-Object { Unblock-File -LiteralPath $_.FullName -ErrorAction SilentlyContinue; $n++ }', nl, ...
        'Write-Output ("Unblocked " + $n + " files")', nl];

    tmpScript = [tempname(), '.ps1'];
    fid = fopen(tmpScript, 'w');
    fwrite(fid, script);
    fclose(fid);

    cmd = sprintf('powershell -NoProfile -ExecutionPolicy Bypass -File "%s" 2>&1', tmpScript);
    [exitCode, output] = system(cmd);
    delete(tmpScript);

    if exitCode ~= 0
        logWarn("installExtra: Unblock-File failed (exit %d): %s", exitCode, strtrim(output));
    else
        logInfo("installExtra: %s", strtrim(output));
    end
end

% -------------------------------------------------------------------------
function writeSiteCustomize_(pyExe)
% Write (or refresh) sitecustomize.py in site-packages so that every
% Python session started from this Embedded Python installation
% automatically registers DLL search directories on startup.
%
% Python 3.8+ no longer searches PATH for DLLs loaded by .pyd extension
% modules.  Packages like openmm bundle DLLs in subdirectories
% (e.g. OpenMM.libs/lib/).  Without explicit os.add_dll_directory() calls,
% these DLLs cannot be found.
%
% sitecustomize.py is loaded automatically by the site module during
% Python startup, before user code runs.  Placing it in site-packages
% ensures it is found even for Embedded Python.
    if ~ispc()
        return;
    end
    sitePackages = fullfile(fileparts(char(pyExe)), "Lib", "site-packages");
    if ~isfolder(sitePackages)
        return;
    end

    scPath = fullfile(sitePackages, "sitecustomize.py");
    nl = newline();
    code = [...
        '# EasyMolKit auto-generated -- do not edit manually.', nl, ...
        '# Registers DLL search directories for Windows Python 3.8+ compatibility.', nl, ...
        '# Required by native extension modules that bundle DLLs in subdirectories', nl, ...
        '# of site-packages (e.g. OpenMM.libs/lib/, rdkit.libs/, etc.).', nl, ...
        'import os', nl, nl, ...
        'def _emk_register_dll_dirs():', nl, ...
        '    # Use __file__ (path of this sitecustomize.py) to locate the correct', nl, ...
        '    # site-packages directory.  Avoids picking up the wrong site-packages', nl, ...
        '    # when multiple Python installs coexist (e.g. roaming user install).', nl, ...
        '    sp = os.path.dirname(os.path.abspath(__file__))', nl, ...
        '    for root, dirs, files in os.walk(sp):', nl, ...
        '        if any(f.lower().endswith(''.dll'') for f in files):', nl, ...
        '            try:', nl, ...
        '                os.add_dll_directory(root)', nl, ...
        '            except Exception:', nl, ...
        '                pass', nl, nl, ...
        'if hasattr(os, ''add_dll_directory''):', nl, ...
        '    _emk_register_dll_dirs()', nl];

    fid = fopen(scPath, 'w');
    if fid == -1
        logWarn("installExtra: Could not write sitecustomize.py to: %s", scPath);
        return;
    end
    fwrite(fid, code);
    fclose(fid);
    logInfo("installExtra: sitecustomize.py written -- DLL auto-registration enabled for all Python sessions.");
end

% -------------------------------------------------------------------------
function [vstatus, vout] = runImportCheck_(pyExe, importName)
% Verify that a Python library can be imported.
%
% On Windows (Python 3.8+), .pyd extension modules no longer search PATH
% for their DLL dependencies.  Packages like openmm bundle their DLLs in
% subdirectories (e.g. openmm/, openmm.libs/, openmm/.libs/) that must be
% registered explicitly via os.add_dll_directory().
%
% We walk ALL of site-packages recursively and register every directory
% that contains at least one .dll file.  This is more thorough than adding
% only immediate subdirectories and covers all delvewheel / auditwheel
% bundling layouts.
%
% Using a temp .py file (not -c "...") avoids all cmd.exe quoting issues.
    pyStr = char(pyExe);
    impStr = char(importName);
    if ispc()
        nl = newline();
        code = [...
            'import os, sysconfig', nl, ...
            '# Use sysconfig to locate the current interpreter''s site-packages.', nl, ...
            '# Avoids picking up a wrong site-packages (e.g. roaming user install)', nl, ...
            '# when multiple Python installs coexist and sys.path has both.', nl, ...
            'sp = sysconfig.get_path(''purelib'')', nl, ...
            'if sp and os.path.isdir(sp):', nl, ...
            '    for root, dirs, files in os.walk(sp):', nl, ...
            '        if any(f.lower().endswith(''.dll'') for f in files):', nl, ...
            '            try:', nl, ...
            '                os.add_dll_directory(root)', nl, ...
            '            except Exception:', nl, ...
            '                pass', nl, ...
            'import ', impStr, nl];
        tmpPy = [tempname(), '.py'];
        fid = fopen(tmpPy, 'w');
        fwrite(fid, code);
        fclose(fid);
        [vstatus, vout] = system(['"', pyStr, '" "', tmpPy, '" 2>&1']);
        delete(tmpPy);
    else
        [vstatus, vout] = system(['"', pyStr, '" -c "import ', impStr, '" 2>&1']);
    end
end

% -------------------------------------------------------------------------
function ver = resolveVersion_(cfg, name)
% Return the pinned version string from config, or "" if not specified.
    ver = "";
    if isfield(cfg, "extraLibraries") && isfield(cfg.extraLibraries, char(name))
        ver = string(cfg.extraLibraries.(char(name)));
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
function installTorchGeometric_(pyExe, proxyArg)
% Detect the installed torch version, build the PyG find-links URL, and
% install torch_geometric + companion packages (torch_scatter, torch_sparse,
% torch_cluster) in a single pip call.
%
% Requires torch to be installed first (emk.setup.installExtra("torch")).
% The PyG find-links URL pattern is:
%   https://data.pyg.org/whl/torch-<X.Y.Z>+cpu.html
% Reference: https://data.pyg.org/whl/

    % Detect torch version
    torchVerCmd = ['"', char(pyExe), '" -c "import torch; print(torch.__version__)" 2>&1'];
    [st, torchVerRaw] = system(torchVerCmd);
    if st ~= 0
        error("emk:setup:installExtra:installFailed", ...
            "torch is not installed in the Embedded Python environment.\n" + ...
            "Install torch first:\n" + ...
            "  emk.setup.installExtra(""torch"")");
    end

    % Strip build suffix, e.g. "2.5.0+cpu" -> "2.5.0"
    torchVer = strtrim(torchVerRaw);
    torchVer = regexprep(torchVer, '\+.*$', '');
    pygUrl   = "https://data.pyg.org/whl/torch-" + torchVer + "+cpu.html";
    logInfo("installExtra: Detected torch %s, PyG wheel URL: %s", torchVer, pygUrl);
    logWarn("installExtra: torch_geometric + extras may take several minutes (~300 MB).");

    % Install torch_geometric, torch_scatter, torch_sparse, torch_cluster
    pkgs = "torch_geometric torch_scatter torch_sparse torch_cluster";
    cmd  = ['"', char(pyExe), '" -m pip install ', char(pkgs), ...
        ' --no-warn-script-location', char(proxyArg), ' -f ', char(pygUrl)];
    logInfo("installExtra: Installing torch_geometric + scatter/sparse/cluster ...");
    [status, output] = system(cmd);
    if status ~= 0
        error("emk:setup:installExtra:installFailed", ...
            "pip install torch_geometric failed (exit code %d):\n%s", status, output);
    end
    logInfo("installExtra: torch_geometric + scatter/sparse/cluster installed successfully.");

    % Unblock downloaded DLLs before import
    unblockPydFiles_(pyExe);

    % Write/refresh sitecustomize.py for DLL auto-registration.
    writeSiteCustomize_(pyExe);

    % Verify import
    verifyCmd = ['"', char(pyExe), '" -c "import torch_geometric" 2>&1'];
    [vstatus, vout] = system(verifyCmd);
    if vstatus ~= 0
        hint = "";
        if ispc() && contains(vout, char([12450 12503 12522 12465 12540 12471 12519 12531 21046 24481 12509 12522 12471 12540]))
            % char([...]) = Japanese: 'Application Control Policy' (Windows JA-JP SAC error)
            hint = newline() + newline() + ...
                "HINT: Smart App Control (SAC) is blocking an unsigned .pyd file." + newline() + ...
                "To fix: Settings > Windows Security > App & Browser Control" + newline() + ...
                "  > Smart App Control > Off  (permanent; cannot revert to Evaluation)." + newline() + ...
                "After turning SAC off, re-run emk.setup.installExtra(""torch_geometric"") " + ...
                "-- no reinstallation needed.";
        end
        error("emk:setup:installExtra:importVerifyFailed", ...
            "Import verification failed for torch_geometric:\n%s%s", vout, hint);
    end
    logInfo("installExtra: import torch_geometric OK.");
end
