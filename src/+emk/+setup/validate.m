function T = validate(varargin)
% validate  Diagnose installed library versions in the active Python environment.
%
%   T = emk.setup.validate()
%   T = emk.setup.validate(Libraries=["rdkit","pubchempy"])
%
%   Checks whether each specified library is installed and returns its
%   version using "python -m pip show".  Non-throwing: unavailable libraries
%   appear as Installed=false rather than causing an error.
%
%   Arguments (name-value):
%     Libraries  (string array) - library names to check.
%                Default: all known libraries (rdkit, pubchempy, mordred,
%                biopython, scipy, meeko, vina, pdbfixer,
%                openbabel, mdanalysis, pymol).
%
%   Returns:
%     T  table with columns:
%        Library   (string)  - friendly library name
%        Installed (logical) - true when pip reports the package
%        Version   (string)  - version string, or "" if not installed
%        Track     (string)  - "1" (Embedded / pip) or "2" (External Python)
%
%   Error IDs:
%     emk:setup:validate:invalidInput - Libraries is not a string array
%
%   See also: emk.setup.installExtra, emk.setup.useExternal, emk.setup.recipe

    % --- Parse name-value argument manually to support our error ID ---
    p = inputParser();
    p.addParameter("Libraries", string.empty, @(x) true);   % validate below
    p.parse(varargin{:});
    libs = p.Results.Libraries;

    % Type validation (throws our error ID, not MATLAB's generic one)
    if ~isstring(libs)
        error("emk:setup:validate:invalidInput", ...
            "Libraries must be a string array, got %s.", class(libs));
    end

    % Default: all known libraries
    if isempty(libs)
        libs = defaultLibraryList_();
    end

    % Flatten to column vector for uniform iteration
    libs = libs(:);

    % --- Locate Python executables (non-throwing) ---
    % Track 1: embedded/active pyenv executable
    % Track 2: external_path from settings.json (may differ from active pyenv)
    pyExe1 = resolvePythonExe_();
    pyExe2 = resolveExternalPyExe_();  % "" when no external_path configured

    % --- Build result arrays ---
    nLibs     = numel(libs);
    Library   = strings(nLibs, 1);
    Installed = false(nLibs, 1);
    Version   = strings(nLibs, 1);
    Track     = strings(nLibs, 1);

    for i = 1:nLibs
        libName    = libs(i);
        meta       = resolveLibMeta_(libName);
        Library(i) = libName;
        Track(i)   = meta.track;

        % Choose the right Python for each track.
        % Track 2 libraries live in an external venv and must be checked
        % using external_path, not the embedded Python (pyenv.Executable).
        if strcmp(meta.track, "2")
            pyExe = pyExe2;
        else
            pyExe = pyExe1;
        end

        if strlength(pyExe) == 0
            % Python executable could not be determined; skip pip check
            continue;
        end

        [ok, ver]    = pipShow_(pyExe, meta.pipShowName);
        Installed(i) = ok;
        Version(i)   = ver;
    end

    T = table(Library, Installed, Version, Track);

    logInfo("validate: Checked %d %s.", nLibs, ...
        choose_(nLibs == 1, "library", "libraries"));
end

% =========================================================================
% Private helpers
% =========================================================================

function pyExe = resolveExternalPyExe_()
% Return the external Python executable path from settings.json (non-throwing).
% Used for Track 2 library checks.  Returns "" when not configured.
    pyExe = "";
    try
        projectRoot = fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))));
        settingsPath = fullfile(projectRoot, "config", "settings.json");
        if ~isfile(settingsPath)
            return;
        end
        raw = jsondecode(fileread(settingsPath));
        if isfield(raw, "python") && isfield(raw.python, "external_path")
            candidate = string(raw.python.external_path);
            if strlength(candidate) > 0 && isfile(candidate)
                pyExe = candidate;
            end
        end
    catch
        % Ignore all errors; return ""
    end
end

% -------------------------------------------------------------------------
function pyExe = resolvePythonExe_()
% Return the active Python executable path (non-throwing).
% Priority: pyenv.Executable (if loaded) > initPython + pyenv > embedded fallback.
    pyExe = "";

    try
        pe = pyenv();
        if strcmp(string(pe.Status), "NotLoaded")
            try
                emk.setup.initPython();
                pe = pyenv();
            catch
                % initPython may fail if embedded Python is missing
            end
        end
        exePath = string(pe.Executable);
        if strlength(exePath) > 0 && isfile(exePath)
            pyExe = exePath;
            return;
        end
    catch
        % Ignore all errors; fall through to embedded fallback
    end

    % Last resort: try embedded Python directly
    try
        projectRoot = fileparts(fileparts(fileparts(fileparts(mfilename("fullpath")))));
        cfg = emkLoadConfig();
        candidate = fullfile(projectRoot, cfg.python.embedded_dir, "python.exe");
        if isfile(candidate)
            pyExe = candidate;
        end
    catch
        % Give up
    end
end

% -------------------------------------------------------------------------
function [ok, ver] = pipShow_(pyExe, pkgName)
% Run "python -m pip show <pkgName>" and parse the Version field.
    cmd = sprintf('"%s" -m pip show %s 2>nul', pyExe, pkgName);
    [status, output] = system(cmd);

    if status ~= 0
        ok  = false;
        ver = "";
        return;
    end

    tok = regexp(output, 'Version:\s*(\S+)', 'tokens', 'once');
    if isempty(tok)
        ok  = false;
        ver = "";
    else
        ok  = true;
        ver = string(tok{1});
    end
end

% -------------------------------------------------------------------------
function libs = defaultLibraryList_()
    libs = ["rdkit"; "pubchempy"; "mordred"; "biopython"; ...
            "scipy"; "meeko"; "vina"; "pdbfixer"; ...
            "openbabel"; "mdanalysis"; "pymol"];
end

% -------------------------------------------------------------------------
function meta = resolveLibMeta_(name)
% Return pip-show package name and track classification for a known library.
    switch name
        case "rdkit"
            meta = struct(pipShowName="rdkit",           track="1");
        case "pubchempy"
            meta = struct(pipShowName="pubchempy",        track="1");
        case "mordred"
            meta = struct(pipShowName="mordredcommunity", track="1");
        case "biopython"
            meta = struct(pipShowName="biopython",        track="1");
        case "scipy"
            meta = struct(pipShowName="scipy",             track="1");
        case "meeko"
            meta = struct(pipShowName="meeko",             track="1");
        case "vina"
            meta = struct(pipShowName="vina",              track="1");
        case "pdbfixer"
            meta = struct(pipShowName="pdbfixer",          track="1");
        case "openbabel"
            meta = struct(pipShowName="openbabel",        track="2");
        case "mdanalysis"
            meta = struct(pipShowName="MDAnalysis",       track="2");
        case "pymol"
            meta = struct(pipShowName="pymol",            track="2");
        otherwise
            % Unknown library: use name as-is with unknown track
            meta = struct(pipShowName=char(name),         track="?");
    end
end

% -------------------------------------------------------------------------
function s = choose_(cond, a, b)
% Return a if cond is true, b otherwise (ternary helper).
    if cond
        s = a;
    else
        s = b;
    end
end
