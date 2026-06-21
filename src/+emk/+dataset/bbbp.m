function tbl = bbbp(varargin)
% bbbp  Load the BBBP blood-brain barrier permeability dataset.
%
%   tbl = emk.dataset.bbbp()
%   tbl = emk.dataset.bbbp(CacheDir=path)
%   tbl = emk.dataset.bbbp(ForceDownload=true)
%
%   Downloads the BBBP dataset on first call and caches it locally.
%   Returns a MATLAB table of 2039 drugs labelled for CNS penetrance.
%
%   Source: Martins et al. (2012); DeepChem/MoleculeNet distribution.
%   License: See THIRD_PARTY_NOTICES.md.
%
%   Arguments:
%     CacheDir       - string (optional)  Default: <project_root>/data/benchmark/
%     ForceDownload  - logical (optional, default false)
%
%   Returns:
%     tbl - table(2039x3) with columns:
%       SMILES  (string)   Canonical SMILES
%       Name    (string)   Compound name
%       BBB     (logical)  true = BBB-permeable (CNS+), false = non-permeable
%
%   Errors:
%     emk:dataset:bbbp:invalidInput   - invalid option
%     emk:dataset:bbbp:downloadFailed - network error
%     emk:dataset:bbbp:parseFailed    - CSV parsing error
%
%   Example:
%     tbl = emk.dataset.bbbp();
%     fprintf("BBB+ : %d / %d\n", sum(tbl.BBB), height(tbl));
%     % Apply Lipinski filter to BBB+ subset:
%     pos = tbl(tbl.BBB, :);
%
%   References:
%     Martins, I.F. et al. (2012). A Bayesian Approach to in Silico
%       Blood-Brain Barrier Penetration Modeling. J. Chem. Inf. Model.
%       52(6):1686-1697. DOI: 10.1021/ci300124c
%     Wu, Z. et al. (2018). MoleculeNet. Chem. Sci. 9:513-530.
%       DOI: 10.1039/C7SC02664A
%
%   See also: emk.dataset.esol, emk.dataset.freesolv, emk.dataset.tox21

    % --- Parse optional name-value arguments ---
    cacheDir      = "";
    forceDownload = false;
    if nargin > 0
        if mod(numel(varargin), 2) ~= 0
            error("emk:dataset:bbbp:invalidInput", ...
                "Optional arguments must be name-value pairs.");
        end
        for k = 1:2:numel(varargin)
            argName = string(varargin{k});
            switch lower(argName)
                case "cachedir"
                    cacheDir = string(varargin{k+1});
                case "forcedownload"
                    forceDownload = logical(varargin{k+1});
                otherwise
                    error("emk:dataset:bbbp:invalidInput", ...
                        "Unknown option: '%s'", varargin{k});
            end
        end
    end

    % --- Resolve cache directory ---
    if strlength(cacheDir) == 0
        cacheDir = defaultCacheDir_();
    end
    if ~isfolder(cacheDir)
        mkdir(cacheDir);
    end

    csvFile = fullfile(cacheDir, "bbbp.csv");

    % --- Download if needed ---
    if ~isfile(csvFile) || forceDownload
        url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv";
        logInfo("dataset.bbbp: downloading BBBP dataset from DeepChem...");
        try
            websave(csvFile, url);
            logInfo("dataset.bbbp: saved to %s", csvFile);
        catch ME
            error("emk:dataset:bbbp:downloadFailed", ...
                "Failed to download BBBP dataset from %s: %s", url, ME.message);
        end
    else
        logDebug("dataset.bbbp: loading from cache %s", csvFile);
    end

    % --- Parse CSV ---
    try
        raw = readtable(csvFile, "TextType", "string", ...
                        "VariableNamingRule", "preserve");
    catch ME
        error("emk:dataset:bbbp:parseFailed", ...
            "Failed to parse %s: %s", csvFile, ME.message);
    end

    % --- Standardize columns ---
    % DeepChem BBBP columns: "num", "name", "p_np", "smiles"
    %   p_np: 1 = permeable (BBB+), 0 = non-permeable (BBB-)
    tbl = table();
    tbl.SMILES = raw.("smiles");
    tbl.Name   = raw.("name");
    tbl.BBB    = logical(raw.("p_np"));

    logInfo("dataset.bbbp: loaded %d molecules (%d BBB+, %d BBB-)", ...
        height(tbl), sum(tbl.BBB), sum(~tbl.BBB));
end

% -------------------------------------------------------------------------
function d = defaultCacheDir_()
    logPath = which("logInfo");
    if isempty(logPath)
        d = fullfile(pwd(), "data", "benchmark");
        return;
    end
    [srcUtil, ~] = fileparts(logPath);
    [src,     ~] = fileparts(srcUtil);
    [root,    ~] = fileparts(src);
    d = fullfile(root, "data", "benchmark");
end
