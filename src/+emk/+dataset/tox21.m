function tbl = tox21(varargin)
% tox21  Load the Tox21 multi-task toxicity dataset.
%
%   tbl = emk.dataset.tox21()
%   tbl = emk.dataset.tox21(CacheDir=path)
%   tbl = emk.dataset.tox21(ForceDownload=true)
%
%   Downloads the Tox21 dataset on first call and caches it locally.
%   Returns a MATLAB table of ~7831 compounds labelled across 12 in-vitro
%   toxicology assays.  Many labels are missing (NaN) due to the
%   multi-lab nature of the dataset.
%
%   Source: Tox21 Data Challenge 2014; MoleculeNet distribution.
%   License: Public domain (US federal government data).
%            See THIRD_PARTY_NOTICES.md.
%
%   Arguments:
%     CacheDir       - string (optional)  Default: <project_root>/data/benchmark/
%     ForceDownload  - logical (optional, default false)
%
%   Returns:
%     tbl - table(N x 14) with columns:
%       SMILES   (string)  Canonical SMILES
%       MolID    (string)  Molecule identifier (mol_id)
%       NR_AR    (double)  Androgen receptor (1=active, 0=inactive, NaN=missing)
%       NR_AR_LBD (double) Androgen receptor LBD
%       NR_AhR   (double)  AhR agonist
%       NR_Aromatase (double) Aromatase inhibitor
%       NR_ER    (double)  Estrogen receptor alpha
%       NR_ER_LBD (double) Estrogen receptor alpha LBD
%       NR_PPAR_gamma (double) PPAR-gamma agonist
%       SR_ARE   (double)  Antioxidant response element
%       SR_ATAD5 (double)  Genotoxicity (ATAD5)
%       SR_HSE   (double)  Heat shock response
%       SR_MMP   (double)  Mitochondrial membrane potential
%       SR_p53   (double)  p53 signaling
%
%   Errors:
%     emk:dataset:tox21:invalidInput   - invalid option
%     emk:dataset:tox21:downloadFailed - network error
%     emk:dataset:tox21:parseFailed    - CSV parsing error
%
%   Example:
%     tbl = emk.dataset.tox21();
%     % Count labelled entries per endpoint
%     endpoints = tbl.Properties.VariableNames(3:end);
%     for e = endpoints
%       nLabelled = sum(~isnan(tbl.(e{1})));
%       fprintf("%s: %d labelled\n", e{1}, nLabelled);
%     end
%
%   References:
%     Tox21 Data Challenge. http://tripod.nih.gov/tox21/challenge/ (2014)
%     Huang, R. et al. (2016). Tox21Challenge to Build Predictive Models
%       of Nuclear Receptor and Stress Response Pathways as Mediated by
%       Exposure to Environmental Chemicals and Drugs. Front. Environ.
%       Sci. 3:85. DOI: 10.3389/fenvs.2015.00085
%     Wu, Z. et al. (2018). MoleculeNet. Chem. Sci. 9:513-530.
%       DOI: 10.1039/C7SC02664A
%
%   See also: emk.dataset.esol, emk.dataset.freesolv, emk.dataset.bbbp

    % --- Parse optional name-value arguments ---
    cacheDir      = "";
    forceDownload = false;
    if nargin > 0
        if mod(numel(varargin), 2) ~= 0
            error("emk:dataset:tox21:invalidInput", ...
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
                    error("emk:dataset:tox21:invalidInput", ...
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

    csvFile = fullfile(cacheDir, "tox21.csv");

    % --- Download if needed ---
    if ~isfile(csvFile) || forceDownload
        url = "https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/tox21.csv.gz";
        gzFile = fullfile(cacheDir, "tox21.csv.gz");
        logInfo("dataset.tox21: downloading Tox21 dataset from DeepChem...");
        try
            websave(gzFile, url);
            % Decompress .gz -> .csv (gunzip removes the .gz suffix)
            gunzip(gzFile, cacheDir);
            logInfo("dataset.tox21: saved to %s", csvFile);
        catch ME
            error("emk:dataset:tox21:downloadFailed", ...
                "Failed to download Tox21 dataset from %s: %s", url, ME.message);
        end
    else
        logDebug("dataset.tox21: loading from cache %s", csvFile);
    end

    % --- Parse CSV ---
    try
        raw = readtable(csvFile, "TextType", "string", ...
                        "VariableNamingRule", "modify");
    catch ME
        error("emk:dataset:tox21:parseFailed", ...
            "Failed to parse %s: %s", csvFile, ME.message);
    end

    % --- Build standardized table ---
    % Tox21 columns (MoleculeNet): smiles, mol_id, NR-AR, NR-AR-LBD,
    %   NR-AhR, NR-Aromatase, NR-ER, NR-ER-LBD, NR-PPAR-gamma,
    %   SR-ARE, SR-ATAD5, SR-HSE, SR-MMP, SR-p53
    tbl = table();

    % SMILES and ID columns
    varNames = string(raw.Properties.VariableNames);
    smilesCol = findColName_(varNames, ["smiles", "SMILES"]);
    idCol     = findColName_(varNames, ["mol_id", "molID", "ID"]);

    tbl.SMILES = raw.(smilesCol);
    tbl.MolID  = raw.(idCol);

    % Toxicity endpoint columns: map hyphenated names to valid MATLAB names
    endpointMap = {
        "NR-AR",         "NR_AR";
        "NR-AR-LBD",     "NR_AR_LBD";
        "NR-AhR",        "NR_AhR";
        "NR-Aromatase",  "NR_Aromatase";
        "NR-ER",         "NR_ER";
        "NR-ER-LBD",     "NR_ER_LBD";
        "NR-PPAR-gamma", "NR_PPAR_gamma";
        "SR-ARE",        "SR_ARE";
        "SR-ATAD5",      "SR_ATAD5";
        "SR-HSE",        "SR_HSE";
        "SR-MMP",        "SR_MMP";
        "SR-p53",        "SR_p53"
    };

    for i = 1:size(endpointMap, 1)
        srcName = endpointMap{i, 1};
        dstName = endpointMap{i, 2};
        % readtable with "modify" replaces '-' with '_', so column names may differ
        % Try both hyphenated and modified versions
        srcCandidates = [srcName, strrep(srcName, "-", "_")];
        colFound = "";
        for c = srcCandidates
            if ismember(c, varNames)
                colFound = c;
                break;
            end
        end
        if strlength(colFound) > 0
            rawCol = raw.(colFound);
            % Convert to double; empty strings -> NaN
            if ischar(rawCol) || iscell(rawCol) || isstring(rawCol)
                nums = double(rawCol);
                nums(rawCol == "" | rawCol == "nan") = NaN;
                tbl.(dstName) = nums;
            else
                tbl.(dstName) = double(rawCol);
            end
        else
            tbl.(dstName) = NaN(height(raw), 1);
            logWarn("dataset.tox21: column '%s' not found in CSV; filled with NaN", srcName);
        end
    end

    logInfo("dataset.tox21: loaded %d molecules, 12 toxicity endpoints", height(tbl));
end

% -------------------------------------------------------------------------
function name = findColName_(varNames, candidates)
    for c = candidates
        if ismember(c, varNames)
            name = c;
            return;
        end
    end
    % Fallback: return first candidate (will cause an error on access if missing)
    name = candidates(1);
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
