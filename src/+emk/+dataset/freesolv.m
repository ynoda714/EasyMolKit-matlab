function tbl = freesolv(varargin)
% freesolv  Load the FreeSolv free energy of hydration dataset.
%
%   tbl = emk.dataset.freesolv()
%   tbl = emk.dataset.freesolv(CacheDir=path)
%   tbl = emk.dataset.freesolv(ForceDownload=true)
%
%   Downloads the FreeSolv dataset on first call and caches it locally.
%   Returns a MATLAB table of 642 small molecules with experimental and
%   computed free energies of hydration.
%
%   Source: Mobley & Guthrie (2014); DeepChem/MoleculeNet distribution.
%   License: CC0 (public domain).  See THIRD_PARTY_NOTICES.md.
%
%   Arguments:
%     CacheDir       - string (optional)  Local cache directory.
%                      Default: <project_root>/data/benchmark/
%     ForceDownload  - logical (optional, default false)
%
%   Returns:
%     tbl - table(642x5) with columns:
%       SMILES        (string)  Canonical SMILES
%       Name          (string)  IUPAC name
%       DeltaG_exp    (double)  Experimental hydration free energy [kcal/mol]
%       DeltaG_calc   (double)  Calculated (GAFF+TIP3P) free energy [kcal/mol]
%       DeltaG_exp_sem (double) Standard error of the mean [kcal/mol]
%
%   Errors:
%     emk:dataset:freesolv:invalidInput   - invalid option argument
%     emk:dataset:freesolv:downloadFailed - network error
%     emk:dataset:freesolv:parseFailed    - CSV parsing error
%
%   Example:
%     tbl = emk.dataset.freesolv();
%     % Scatter: calculated vs experimental
%     scatter(tbl.DeltaG_exp, tbl.DeltaG_calc, "filled");
%     xlabel("Experimental \DeltaG [kcal/mol]");
%     ylabel("Calculated \DeltaG [kcal/mol]");
%
%   References:
%     Mobley, D.L. & Guthrie, J.P. (2014). FreeSolv: a database of
%       experimental and calculated hydration free energies, with input
%       files. J. Comput.-Aided Mol. Des. 28(7):711-720.
%       DOI: 10.1007/s10822-014-9747-x
%     Wu, Z. et al. (2018). MoleculeNet. Chem. Sci. 9:513-530.
%       DOI: 10.1039/C7SC02664A
%
%   See also: emk.dataset.esol, emk.dataset.bbbp, emk.dataset.tox21

    % --- Parse optional name-value arguments ---
    cacheDir      = "";
    forceDownload = false;
    if nargin > 0
        if mod(numel(varargin), 2) ~= 0
            error("emk:dataset:freesolv:invalidInput", ...
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
                    error("emk:dataset:freesolv:invalidInput", ...
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

    csvFile = fullfile(cacheDir, "freesolv.csv");

    % --- Download if needed ---
    if ~isfile(csvFile) || forceDownload
        url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/SAMPL.csv";
        logInfo("dataset.freesolv: downloading FreeSolv dataset from DeepChem...");
        try
            websave(csvFile, url);
            logInfo("dataset.freesolv: saved to %s", csvFile);
        catch ME
            error("emk:dataset:freesolv:downloadFailed", ...
                "Failed to download FreeSolv dataset from %s: %s", url, ME.message);
        end
    else
        logDebug("dataset.freesolv: loading from cache %s", csvFile);
    end

    % --- Parse CSV ---
    try
        raw = readtable(csvFile, "TextType", "string", ...
                        "VariableNamingRule", "preserve");
    catch ME
        error("emk:dataset:freesolv:parseFailed", ...
            "Failed to parse %s: %s", csvFile, ME.message);
    end

    % --- Standardize columns ---
    % DeepChem SAMPL (FreeSolv) columns:
    %   "iupac", "smiles", "expt", "calc", "expt_unc" (may vary by version)
    % We handle both possible column name variants.
    varNames = string(raw.Properties.VariableNames);
    tbl = table();

    tbl.SMILES = getCol_(raw, varNames, ["smiles", "SMILES"]);
    tbl.Name   = getCol_(raw, varNames, ["iupac", "IUPAC", "name", "Name"]);

    expCol  = getColName_(varNames, ["expt", "experimental", "exp"]);
    calcCol = getColName_(varNames, ["calc", "calculated", "GAFF"]);
    semCol  = getColName_(varNames, ["expt_unc", "exp_unc", "uncertainty"]);

    tbl.DeltaG_exp  = raw.(expCol);
    tbl.DeltaG_calc = raw.(calcCol);
    if ~isempty(semCol)
        tbl.DeltaG_exp_sem = raw.(semCol);
    else
        tbl.DeltaG_exp_sem = NaN(height(raw), 1);
    end

    logInfo("dataset.freesolv: loaded %d molecules", height(tbl));
end

% -------------------------------------------------------------------------
function col = getCol_(raw, varNames, candidates)
    for c = candidates
        idx = find(lower(varNames) == lower(c), 1);
        if ~isempty(idx)
            col = raw.(varNames(idx));
            return;
        end
    end
    col = repmat("", height(raw), 1);
end

function name = getColName_(varNames, candidates)
    for c = candidates
        idx = find(lower(varNames) == lower(c), 1);
        if ~isempty(idx)
            name = varNames(idx);
            return;
        end
    end
    name = "";
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
