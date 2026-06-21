function tbl = esol(varargin)
% esol  Load the ESOL (Delaney) aqueous solubility dataset.
%
%   tbl = emk.dataset.esol()
%   tbl = emk.dataset.esol(CacheDir=path)
%   tbl = emk.dataset.esol(ForceDownload=true)
%
%   Downloads the ESOL dataset on first call and caches it locally.
%   Returns a MATLAB table of 1128 small organic molecules with
%   measured aqueous solubility (logS) and predicted values.
%
%   Source: Delaney (2004), ESOL dataset hosted by DeepChem/MoleculeNet.
%   License: Public domain (widely used benchmark; no explicit license
%             for the Delaney data itself; see THIRD_PARTY_NOTICES.md).
%
%   Arguments:
%     CacheDir       - string (optional)  Directory to store the cached CSV.
%                      Default: <project_root>/data/benchmark/
%     ForceDownload  - logical (optional, default false)  Re-download even
%                      if a cached file already exists.
%
%   Returns:
%     tbl - table(1128x5) with columns:
%       SMILES        (string)  Canonical SMILES
%       Name          (string)  Compound name / identifier
%       logS          (double)  Measured log10 aqueous solubility [mol/L]
%       logS_Delaney  (double)  Delaney ESOL predicted logS [mol/L]
%       MolWt         (double)  Molecular weight from dataset [g/mol]
%
%   Errors:
%     emk:dataset:esol:invalidInput    - invalid option argument
%     emk:dataset:esol:downloadFailed  - network error during download
%     emk:dataset:esol:parseFailed     - CSV parsing error
%
%   Example:
%     tbl = emk.dataset.esol();
%     % Filter to moderate solubility window
%     mod = tbl(tbl.logS >= -4 & tbl.logS <= -2, :);
%     summary(mod)
%
%   References:
%     Delaney, J.S. (2004). ESOL: Estimating Aqueous Solubility Directly
%       from Molecular Structure. J. Chem. Inf. Comput. Sci. 44(3):1000-1005.
%       DOI: 10.1021/ci034243x
%     Wu, Z. et al. (2018). MoleculeNet: A Benchmark for Molecular Machine
%       Learning. Chem. Sci. 9:513-530. DOI: 10.1039/C7SC02664A
%
%   See also: emk.dataset.freesolv, emk.dataset.bbbp, emk.dataset.tox21

    % --- Parse optional name-value arguments ---
    cacheDir      = "";
    forceDownload = false;
    if nargin > 0
        if mod(numel(varargin), 2) ~= 0
            error("emk:dataset:esol:invalidInput", ...
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
                    error("emk:dataset:esol:invalidInput", ...
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

    csvFile = fullfile(cacheDir, "esol.csv");

    % --- Download if needed ---
    if ~isfile(csvFile) || forceDownload
        url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv";
        logInfo("dataset.esol: downloading ESOL dataset from DeepChem...");
        try
            websave(csvFile, url);
            logInfo("dataset.esol: saved to %s", csvFile);
        catch ME
            error("emk:dataset:esol:downloadFailed", ...
                "Failed to download ESOL dataset from %s: %s", url, ME.message);
        end
    else
        logDebug("dataset.esol: loading from cache %s", csvFile);
    end

    % --- Parse CSV ---
    try
        raw = readtable(csvFile, "TextType", "string", ...
                        "VariableNamingRule", "preserve");
    catch ME
        error("emk:dataset:esol:parseFailed", ...
            "Failed to parse %s: %s", csvFile, ME.message);
    end

    % --- Standardize column names ---
    % Original DeepChem ESOL columns:
    %   "Compound ID", "ESOL predicted log solubility in mols per litre",
    %   "Minimum Degree", "Molecular Weight", "Number of H-Bond Donors",
    %   "Number of Rings", "Number of Rotatable Bonds", "Polar Surface Area",
    %   "measured log solubility in mols per litre", "smiles"
    tbl = table();
    tbl.SMILES       = raw.("smiles");
    tbl.Name         = raw.("Compound ID");
    tbl.logS         = raw.("measured log solubility in mols per litre");
    tbl.logS_Delaney = raw.("ESOL predicted log solubility in mols per litre");
    tbl.MolWt        = raw.("Molecular Weight");

    logInfo("dataset.esol: loaded %d molecules", height(tbl));
end

% -------------------------------------------------------------------------
function d = defaultCacheDir_()
% Resolve <project_root>/data/benchmark/ relative to logInfo.m location.
    logPath = which("logInfo");
    if isempty(logPath)
        d = fullfile(pwd(), "data", "benchmark");
        return;
    end
    % logInfo is at <root>/src/util/logInfo.m  -> 3 fileparts to root
    [srcUtil, ~] = fileparts(logPath);
    [src,     ~] = fileparts(srcUtil);
    [root,    ~] = fileparts(src);
    d = fullfile(root, "data", "benchmark");
end
