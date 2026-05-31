function cfg = emkLoadConfig()
% emkLoadConfig  Load EasyMolKit configuration with environment variable overrides.
%
%   cfg = emkLoadConfig()
%
%   Priority (highest to lowest):
%     1. Environment variables: EMK_<SECTION>_<KEY>  (e.g., EMK_PYTHON_VERSION)
%     2. config/settings.json  (if present)
%     3. Built-in defaults
%
%   Returns:
%     cfg (struct) - Configuration struct with fields:
%       .python.version         - Python version string (default: "3.10")
%       .python.embedded_dir    - Embedded Python directory (default: "python_env")
%       .python.execution_mode  - pyenv execution mode (default: "OutOfProcess")
%       .python.proxy           - Proxy URL for pip (default: "")
%       .rdkit.version          - RDKit version to install (default: "2024.03.6")
%       .rdkit.auto_init        - Auto-init Python on setup (default: true)
%       .runtime.eval_mode      - "prod" or "ref" (default: "prod")
%       .runtime.text_only_mode - Suppress figure output (default: true)
%       .output.root_dir        - Intermediate output dir (default: "result/intermediate")
%       .run.root_dir           - Run output dir (default: "result/runs")
%       .run.publish_latest     - Publish latest symlink (default: true)

    cfg = buildDefaults();
    cfg = applyJsonFile(cfg);
    cfg = applyEnvVars(cfg);
end

% =========================================================================
function cfg = buildDefaults()
    cfg.python.version        = "3.10";
    cfg.python.embedded_dir   = "python_env";
    cfg.python.execution_mode = "OutOfProcess";
    cfg.python.proxy          = "";
    cfg.python.external_path  = "";   % Track 2: path to external venv python.exe (ADR-007)

    cfg.rdkit.version         = "2024.03.6";
    cfg.rdkit.auto_init       = true;

    cfg.runtime.eval_mode       = "prod";
    cfg.runtime.text_only_mode  = true;

    cfg.output.root_dir = "result/intermediate";

    cfg.run.root_dir        = "result/runs";
    cfg.run.publish_latest  = true;

    cfg.extraLibraries.pubchempy = "1.0.4";
    cfg.extraLibraries.mordred   = "2.0.7";
    cfg.extraLibraries.biopython = "1.84";

    % Use-case group toggles (default: all off).
    % Each group expands into the corresponding cfg.optionalLibraries.* flags
    % via emk.setup.expandUseCases().  Set to true in main*.m or settings.json.
    %   qsar    -> pubchempy + mordred
    %   bio     -> biopython
    %   ml      -> torch (CPU-only) + torch_geometric
    %   nlp     -> transformers + datasets
    %   docking -> scipy + meeko + vina + pdbfixer  (Online only)
    cfg.useCase.qsar    = false;
    cfg.useCase.bio     = false;
    cfg.useCase.ml      = false;
    cfg.useCase.nlp     = false;
    cfg.useCase.docking = false;

    % Optional library installation toggles (default: all off).
    % Prefer cfg.useCase.* above for group installs; set individual flags here
    % only when you need a single library without its group companions.
    cfg.optionalLibraries.pubchempy       = false;
    cfg.optionalLibraries.mordred         = false;
    cfg.optionalLibraries.biopython       = false;
    cfg.optionalLibraries.torch           = false;
    cfg.optionalLibraries.torch_geometric = false;
    cfg.optionalLibraries.transformers    = false;
    cfg.optionalLibraries.datasets        = false;
    cfg.optionalLibraries.scipy           = false;
    cfg.optionalLibraries.meeko           = false;
    cfg.optionalLibraries.vina            = false;
    cfg.optionalLibraries.pdbfixer        = false;
end

% =========================================================================
function cfg = applyJsonFile(cfg)
    jsonPath = fullfile("config", "settings.json");
    if ~isfile(jsonPath)
        logDebug("emkLoadConfig: settings.json not found, using defaults");
        return;
    end
    try
        raw = jsondecode(fileread(jsonPath));
        cfg = mergeStruct(cfg, raw);
        logDebug("emkLoadConfig: loaded %s", jsonPath);
    catch ME
        logWarn("emkLoadConfig: failed to parse settings.json (%s), using defaults", ME.message);
    end
end

% =========================================================================
function cfg = applyEnvVars(cfg)
% Apply environment variable overrides.
% Naming convention: EMK_<SECTION>_<KEY> where SECTION and KEY are UPPER_CASE.
% Example: EMK_PYTHON_VERSION overrides cfg.python.version

    sections = fieldnames(cfg);
    for si = 1:numel(sections)
        sec = sections{si};
        if ~isstruct(cfg.(sec))
            continue;
        end
        keys = fieldnames(cfg.(sec));
        for ki = 1:numel(keys)
            key = keys{ki};
            envName = upper("EMK_" + sec + "_" + key);
            val = getenv(envName);
            if isempty(val)
                continue;
            end
            orig = cfg.(sec).(key);
            if islogical(orig)
                cfg.(sec).(key) = strcmp(lower(val), "true") || strcmp(val, "1");
            elseif isnumeric(orig)
                parsed = str2double(val);
                if ~isnan(parsed)
                    cfg.(sec).(key) = parsed;
                end
            else
                cfg.(sec).(key) = string(val);
            end
            logDebug("emkLoadConfig: %s overridden by env var %s", sec + "." + key, envName);
        end
    end
end

% =========================================================================
function out = mergeStruct(base, override)
% mergeStruct  Merge override fields into base; only update fields present in base.
%   Fields absent from base (e.g., JSON comment keys) are ignored.

    out = base;
    fields = fieldnames(override);
    for k = 1:numel(fields)
        f = fields{k};
        if ~isfield(base, f)
            continue; % skip unknown / comment fields
        end
        if isstruct(base.(f)) && isstruct(override.(f))
            out.(f) = mergeStruct(base.(f), override.(f));
        else
            val = override.(f);
            % Coerce char to string when the default type is string
            if isStringScalar(base.(f)) && ischar(val)
                val = string(val);
            end
            out.(f) = val;
        end
    end
end
