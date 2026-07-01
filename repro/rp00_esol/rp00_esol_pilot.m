% rp00_esol_pilot.m  RP00: ESOL Aqueous Solubility -- Calibration Run
%
% Reproduces the linear regression model from Delaney (2004) and
% calibrates RF03 thresholds for RP01+ verification runs.
%
% ROLE: Calibration run (not a verification run).
%   Thresholds derived here (RMSE_CV<=1.20, R2_CV>=0.75) are applied as
%   binding criteria from RP01 onward. RP00 itself is exempt from
%   pass/fail evaluation. See README.md Verification section.
%
%   Paper:   Delaney, J.S. (2004). ESOL: Estimating Aqueous Solubility
%            Directly from Molecular Structure.
%            J. Chem. Inf. Comput. Sci. 44(3):1000-1005.
%            DOI: 10.1021/ci034243x
%
%   Model:   logS ~ LogP + MolWt + NumRotatableBonds + AromaticProportion
%            Linear regression, 5-fold cross-validation.
%
%   Deliverables (M-REPRO-PILOT / updated in M-REPRO-FOUND):
%     RF01 -- Standard reproduction templates: see README.md, repro/Template.md, and repro/Template.jp.md
%     RF02 -- Version lock: emk.setup.snapshot() + emk.setup.lockfile().
%             Captured: MATLAB, Python, RDKit, toolbox versions, dataset
%             SHA-256, git commit.
%     RF03 -- Calibrated thresholds: RMSE_CV <= 1.20, R^2_CV >= 0.75
%
%   Prerequisites:
%     emk.setup.install() has been run once.
%     Statistics and Machine Learning Toolbox (fitlm, cvpartition).
%
%   Run: Ctrl+Enter in MATLAB with this file open (project root as CWD).

%% Section 0: Setup & Environment
% Bootstrap: navigate to project root regardless of launch directory.
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("RP00", "Section 0: Setup & Environment", "Repro Pilot");
emk.setup.initPython();

emk.mol.fromSmiles("C");  % warm-up: initializes Python/RDKit IPC process

% Capture environment info for RF02 version lock
snap = emk.setup.snapshot();

% Seed MATLAB's RNG for reproducibility. This seed covers warm-up only;
% Section 4 re-seeds immediately before cvpartition to guard against any
% MATLAB random state consumed between here and there (e.g., by emk.* internals
% in future versions -- current pipeline has no such consumption).
% Python/RDKit random state is independent and not seeded here; any future
% Python-side randomness must be seeded explicitly and documented.
rng(42, "twister");

logInfo("RP00 setup complete.");

%% Section 1: Load ESOL Dataset
logSection("RP00", "Section 1: Load ESOL Dataset", "Repro Pilot");

tbl    = emk.dataset.esol();
nTotal = height(tbl);
logInfo("Loaded %d molecules (ESOL / Delaney 2004)", nTotal);
logInfo("logS range: [%.2f, %.2f] log mol/L", min(tbl.logS), max(tbl.logS));

% RF02 dataset integrity: record SHA-256 of cached CSV file.
csvFile     = fullfile(pwd(), "data", "benchmark", "esol.csv");
datasetHash = computeFileHash_(csvFile);
logInfo("Dataset SHA-256: %s", datasetHash);

%% Section 2: Parse SMILES & Compute Descriptors
logSection("RP00", "Section 2: Parse SMILES & Compute Descriptors", "Repro Pilot");

% Parse all SMILES to RDKit mol objects
mols      = cell(1, nTotal);
validMask = false(1, nTotal);

for i = 1:nTotal
    try
        mols{i}      = emk.mol.fromSmiles(tbl.SMILES(i));
        validMask(i) = true;
    catch
        logWarn("SMILES parse failed [%d]: %s", i, tbl.SMILES(i));
    end
end

nValid    = sum(validMask);
if nValid == 0
    error("emk:rp00:parseSmiles:noValidMols", "No valid molecules parsed.");
end
validMols = mols(validMask);
yAll      = tbl.logS(validMask);
logInfo("Parsed %d / %d SMILES successfully", nValid, nTotal);

% Batch-compute Delaney (2004) base descriptors via RDKit
% LogP  = RDKit Crippen-Wildman MolLogP (proxy for Delaney's clogP)
% MolWt = average molecular weight (IUPAC atomic weights, implicit H)
% NumRotatableBonds = strict SMARTS (excludes terminal single bonds)
% HeavyAtomCount    = denominator for AromaticProportion
descNames = ["LogP", "MolWt", "NumRotatableBonds", "HeavyAtomCount"];
logInfo("Computing %d descriptors for %d molecules ...", numel(descNames), nValid);
descTbl = emk.descriptor.batchCalculate(validMols, descNames);

% AromaticProportion = NumAromaticAtoms / HeavyAtomCount  [Delaney 2004, Eq. 1]
descTbl.AromaticProportion = calcAromaticProportion_(validMols, descTbl.HeavyAtomCount);
logInfo("Descriptor computation complete.");
logInfo("  LogP range:    [%.2f, %.2f]", min(descTbl.LogP), max(descTbl.LogP));
logInfo("  AroProp range: [%.3f, %.3f]", ...
    min(descTbl.AromaticProportion), max(descTbl.AromaticProportion));

% RotBonds sensitivity: strict (default) vs loose (strict=False) definition.
% Delaney (2004) likely used a non-strict definition that includes terminal
% single bonds (C-OH, C-NH2), which strict SMARTS excludes.
% Non-significant p=0.85 and 25x coefficient underestimate (+0.0026 vs +0.066)
% may reflect this definition mismatch.
nRotStrict = descTbl.NumRotatableBonds;
nRotLoose  = calcRotBondsLoose_(validMols);
if ~any(isnan(nRotLoose))
    nDiff = sum(nRotStrict ~= nRotLoose);
    logInfo("RotBonds strict vs loose: %d/%d mols differ, mean_delta=+%.2f bonds", ...
        nDiff, nValid, mean(nRotLoose - nRotStrict));
    logInfo("  Strict range: [%.0f, %.0f];  Loose range: [%.0f, %.0f]", ...
        min(nRotStrict), max(nRotStrict), min(nRotLoose), max(nRotLoose));
end

% Guard: exclude molecules with NaN descriptors before regression.
% Without this, yPredCV stays at 0 for NaN-descriptor test-fold indices,
% corrupting pooled OOF RMSE and R^2. Track excluded indices for Section 6.
descNanMask = any(ismissing(descTbl), 2);
if any(descNanMask)
    nDescNaN     = sum(descNanMask);
    validIdxAll  = find(validMask);
    descNanInTbl = validIdxAll(descNanMask);  % original indices into tbl
    validMask(descNanInTbl) = false;           % mark for excluded.csv in Section 6
    descTbl   = descTbl(~descNanMask, :);
    validMols = validMols(~descNanMask);
    yAll      = yAll(~descNanMask);
    nValid    = numel(yAll);
    logWarn("Excluded %d molecules with NaN descriptors; n_regression=%d.", nDescNaN, nValid);
else
    descNanInTbl = zeros(0, 1);
end
% finalValidMask: valid after both filtering stages (SMILES parse + descriptor NaN). Used in Section 6.
finalValidMask = validMask;

%% Section 3: Linear Regression Model (Delaney 2004)
logSection("RP00", "Section 3: Linear Regression Model", "Repro Pilot");

% Feature matrix: [LogP | MolWt | RotBonds | AromaticProportion]
% RotBonds (strict) is retained to match Delaney (2004) four-descriptor
% structure, even though it is non-significant here (p=0.85).
% See README.md Discussion for RotBonds definition analysis.
X = [descTbl.LogP, descTbl.MolWt, descTbl.NumRotatableBonds, descTbl.AromaticProportion];
y = yAll;

mdlFull   = fitlm(X, y, "VarNames", {'LogP', 'MolWt', 'RotBonds', 'AroProp', 'logS'});
rmseTrain = sqrt(mdlFull.MSE);
r2Train   = mdlFull.Rsquared.Ordinary;
logInfo("Full model (all %d mols): RMSE=%.4f, R^2=%.4f", nValid, rmseTrain, r2Train);
disp(mdlFull);

%% Section 4: 5-Fold Cross-Validation
logSection("RP00", "Section 4: 5-Fold Cross-Validation", "Repro Pilot");

rng(42, "twister");  % Re-seed before fold assignment; guards against MATLAB RNG state consumed between Section 0 and here (emk.* internals may consume random state in future versions).
cv      = cvpartition(nValid, "KFold", 5);
rmses   = zeros(5, 1);
yPredCV = zeros(nValid, 1);  % accumulate out-of-fold predictions for pooled R^2

for fold = 1:5
    trainIdx = training(cv, fold);
    testIdx  = test(cv, fold);

    mdlFold = fitlm(X(trainIdx, :), y(trainIdx));
    yPred   = predict(mdlFold, X(testIdx, :));

    yPredCV(testIdx) = yPred;

    res         = y(testIdx) - yPred;
    rmses(fold) = sqrt(mean(res .^ 2));

    logInfo("  Fold %d: RMSE=%.4f", fold, rmses(fold));
end

% Fold-average RMSE (reference only): retained for historical comparison. RF03 criterion uses pooled OOF (rmseCV).
rmseCVFoldAvg = mean(rmses);
rmseCVStd     = std(rmses);

% Pooled OOF RMSE: primary metric, consistent with the pooled OOF R^2 method below.
rmseCV = sqrt(mean((y - yPredCV) .^ 2));

% Pooled OOF R^2: aggregate all out-of-fold predictions, then compute R^2
% once against the global mean of y. Avoids fold-size/variance bias
% that arises from per-fold baselines and simple averaging.
resAll = y - yPredCV;
r2CV   = 1 - sum(resAll .^ 2) / sum((y - mean(y)) .^ 2);

logInfo("5-fold CV: RMSE=%.4f (pooled OOF)", rmseCV);
logInfo("5-fold CV: RMSE=%.4f +/- %.4f (fold-avg; reference only -- not the RF03 criterion)", rmseCVFoldAvg, rmseCVStd);
logInfo("5-fold CV: R^2=%.4f (pooled OOF vs global mean)", r2CV);

%% Section 5: Verification (RF03 -- Calibration Only)
logSection("RP00", "Section 5: Verification (Calibration)", "Repro Pilot");

% RF03 thresholds loaded from lock_template.json (single source of truth).
% Tolerance rationale:
%   RMSE <= 1.20: Delaney (2004) ref 0.996 + tolerance 0.20 for RDKit LogP
%   R^2  >= 0.75: RDKit Crippen-Wildman LogP achieves ~0.76 (vs paper ~0.84)
%   R^2 gap is due to LogP implementation: Crippen-Wildman vs Kowwin/ALogPS.
templatePath = fullfile("repro", "rp00_esol", "lock_template.json");
if isfile(templatePath)
    tmpl     = jsondecode(fileread(templatePath));
    rf03crit = tmpl.rf03_criteria;
else
    rf03crit.rmse_cv          = struct("upper", 1.20);
    rf03crit.r2_cv            = struct("lower", 0.75);
    rf03crit.delaney_rmse_ref = 0.996;  % fallback -- must match lock_template.json rf03_criteria.delaney_rmse_ref
    logWarn("lock_template.json not found; using hardcoded RF03 thresholds.");
end
DELANEY_RMSE_REF = rf03crit.delaney_rmse_ref;  % Delaney (2004) Table 2, training RMSE (single source: lock_template.json)

rf03met.rmse_cv  = rmseCV;
rf03met.r2_cv    = r2CV;

rf03result = emk.repro.verify(rf03met, rf03crit);
disp(rf03result.report);

logInfo("Delaney ref RMSE: %.3f (training set, Table 2)", DELANEY_RMSE_REF);

% RP00 is the calibration run. The thresholds above were derived from
% these results. Reporting metric status is informational only.
% RP01 onward uses these thresholds as binding pass/fail criteria.
logInfo("==> RP00 CALIBRATION RUN: thresholds established (not a verification verdict).");
logInfo("    RMSE_CV=%.4f (thr<=%.2f), R2_CV=%.4f (thr>=%.2f)", ...
    rmseCV, rf03crit.rmse_cv.upper, r2CV, rf03crit.r2_cv.lower);

%% Section 6: Save Results
logSection("RP00", "Section 6: Save Results", "Repro Pilot");

runDir = makeRunDir("Prefix", "rp00_esol_pilot");

% Predictions CSV
yPredFull = predict(mdlFull, X);
outTbl = tbl(finalValidMask, ["SMILES", "Name", "logS"]);
outTbl.logS_pred = yPredFull;
outTbl.residual  = y - yPredFull;
outTbl = [outTbl, descTbl];
writetable(outTbl, fullfile(runDir, "predictions.csv"));
logInfo("Predictions saved: %s/predictions.csv", runDir);

% Save excluded molecules with per-molecule reason (SMILES parse failure or descriptor NaN).
if any(~finalValidMask)
    excludedTbl = tbl(~finalValidMask, ["SMILES", "Name", "logS"]);
    nExcluded   = height(excludedTbl);
    reasons     = repmat("SMILES_parse_failed", nExcluded, 1);
    if ~isempty(descNanInTbl)
        excOrigIdx = find(~finalValidMask);
        reasons(ismember(excOrigIdx, descNanInTbl)) = "descriptor_NaN";
    end
    excludedTbl.reason = reasons;
    writetable(excludedTbl, fullfile(runDir, "excluded.csv"));
    logInfo("Excluded molecules: %s/excluded.csv (%d entries)", runDir, nExcluded);
end

% Scatter plot: measured vs. predicted
fig = figure("Name", "RP00 ESOL: Measured vs Predicted");
scatter(y, yPredFull, 25, "filled", "MarkerFaceAlpha", 0.5);
hold on;
lo = floor(min([y; yPredFull]));
hi = ceil(max([y; yPredFull]));
plot([lo, hi], [lo, hi], "k--", "LineWidth", 1.2);
xlabel("Measured logS (log mol/L)");
ylabel("Predicted logS (log mol/L)");
title(sprintf("RP00 ESOL Full-dataset fit (not CV)  RMSE=%.3f  R^2=%.3f", rmseTrain, r2Train));
saveas(fig, fullfile(runDir, "predicted_vs_actual.png"));
close(fig);
logInfo("Figure saved: %s/predicted_vs_actual.png", runDir);

% Metrics JSON
metrics = struct( ...
    "rmse_train",           rmseTrain, ...
    "r2_train",             r2Train, ...
    "rmse_cv",              rmseCV, ...
    "rmse_cv_fold_avg",     rmseCVFoldAvg, ...
    "rmse_cv_std",          rmseCVStd, ...
    "r2_cv_pooled",         r2CV, ...
    "n_molecules",          nValid, ...
    "delaney_rmse_ref",     DELANEY_RMSE_REF, ...
    "rf03_rmse_tolerance",  rf03crit.rmse_cv.upper - DELANEY_RMSE_REF, ...
    "rf03_criteria",        rf03crit, ...
    "rf03_thresholds_met",  rf03result.pass, ...
    "is_calibration_run",   true);
metricsPath = fullfile(runDir, "metrics.json");
fid = fopen(metricsPath, "w");
if fid == -1
    error("emk:rp00:saveMetrics:fopenFailed", "Cannot open for writing: %s", metricsPath);
end
fprintf(fid, "%s\n", jsonencode(metrics, "PrettyPrint", true));
fclose(fid);
logInfo("Metrics saved: %s/metrics.json", runDir);

% Lock snapshot (RF02)
snap.run_date           = char(datetime("now", "Format", "yyyy-MM-dd'T'HH:mm:ss"));
snap.run_dir            = runDir;
snap.dataset_sha256     = datasetHash;
snap.is_calibration_run = true;
snap.rf03_thresholds    = rf03crit;
emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));

logInfo("RP00 ESOL calibration run complete. Output: %s", runDir);

% ==========================================================================
% Local helper functions
% ==========================================================================

function prop = calcAromaticProportion_(mols, heavyAtomCounts)
% Compute Delaney (2004) AromaticProportion = NumAromaticAtoms / HeavyAtomCount.
% Uses a single pyrun batch call to count aromatic atoms across all molecules.
    n = numel(mols);
    prop = zeros(n, 1);

    pyMolList = py.list();
    for i = 1:n
        pyMolList.append(mols{i});
    end

    try
        pyAros = pyrun( ...
            "aros = [float(sum(a.GetIsAromatic() for a in mol.GetAtoms())) for mol in mols]", ...
            "aros", mols=pyMolList);
        % Cast to py.list for robustness across Python/RDKit versions.
        if ~isa(pyAros, 'py.list')
            pyAros = py.list(pyAros);
        end
        aros = double(py.array.array("d", pyAros));
        if numel(aros) ~= n
            error("emk:rp00:calcAroProp:sizeMismatch", ...
                "pyrun returned %d values, expected %d", numel(aros), n);
        end
    catch ME
        error("emk:rp00:calcAroProp:pyrunFailed", ...
            "AromaticProportion batch failed: %s", ME.message);
    end

    for i = 1:n
        nHeavy = heavyAtomCounts(i);
        if isnan(nHeavy) || nHeavy <= 0
            warning("emk:rp00:calcAroProp:zeroHeavyAtoms", ...
                "Molecule %d has HeavyAtomCount=%g; AroProp set to NaN.", i, nHeavy);
            prop(i) = NaN;
        else
            prop(i) = aros(i) / nHeavy;
        end
    end
end

function nRot = calcRotBondsLoose_(mols)
% Compute NumRotatableBonds with strict=False (non-strict / loose definition).
% Loose mode includes terminal single bonds (C-OH, C-NH2) that strict excludes.
% Used for sensitivity analysis of the RotBonds definition vs. Delaney (2004).
    n = numel(mols);
    pyMolList = py.list();
    for i = 1:n
        pyMolList.append(mols{i});
    end
    try
        pyRots = pyrun( ...
            "from rdkit.Chem import rdMolDescriptors; rots = [float(rdMolDescriptors.CalcNumRotatableBonds(m, False)) for m in mols]", ...
            "rots", mols=pyMolList);
        if ~isa(pyRots, 'py.list')
            pyRots = py.list(pyRots);
        end
        nRot = double(py.array.array("d", pyRots));
    catch ME
        nRot = nan(n, 1);
        logWarn("RotBonds loose-mode batch failed (%s); sensitivity check skipped.", ME.message);
    end
end

function h = computeFileHash_(fpath)
% Compute SHA-256 hex digest of a file for RF02 dataset integrity verification.
% Supports Windows (PowerShell Get-FileHash), Linux/MATLAB Online (sha256sum),
% and macOS (shasum -a 256). Returns "unavailable" with a warning on failure.
    fpathChar = char(fpath);
    if ispc()
        fpathSafe = strrep(fpathChar, "'", "''");
        % -LiteralPath prevents glob/wildcard expansion; single-quote escaping handles embedded quotes.
        % Backtick/dollar injection risk is negligible because Get-FileHash does not evaluate the path as code.
        cmd = sprintf( ...
            'powershell -Command "(Get-FileHash -LiteralPath ''%s'' -Algorithm SHA256).Hash"', ...
            fpathSafe);
        [status, out] = system(cmd);
        if status == 0
            h = string(upper(strtrim(out)));
            return;
        end
    else
        % Linux (MATLAB Online): sha256sum; macOS fallback: shasum -a 256.
        % Single-quote the path to avoid backslash-escape expansion issues in MATLAB sprintf.
        % Paths with embedded single quotes are escaped by doubling (same as PowerShell convention).
        fpathSafe = strrep(fpathChar, "'", "''");
        [status, out] = system(sprintf("sha256sum '%s'", fpathSafe));
        if status ~= 0
            [status, out] = system(sprintf("shasum -a 256 '%s'", fpathSafe));
        end
        if status == 0
            parts = strsplit(strtrim(out));
            if numel(parts) >= 1
                h = string(upper(parts{1}));
                return;
            end
        end
    end
    warning("emk:rp00:computeFileHash:failed", ...
        "SHA-256 hash unavailable for '%s' (platform: %s).", fpathChar, computer());
    h = "unavailable";
end
