% rp01_esol.m  RP01: ESOL QSAR Regression + L05 Extended
%
% Reproduces the Delaney (2004) ESOL QSAR model (4-feature linear
% regression: LogP, MolWt, RotBonds, AroProp) and extends with 5
% additional descriptors (L05: TPSA, HBD, HBA, FracCSP3, QED).
%
% Includes M-REPRO-AUDIT diagnostics:
%   1. VIF for Model A (4-feat) and Model B (9-feat) — multicollinearity.
%   2. Paired t-test: Model A vs B per-fold RMSE and R^2.
%   3. TPSA sign verification: partial coefficient vs marginal correlation.
%      (Actual values vary by run -- see execution log for r and coef.)
%
%   Run: Ctrl+Enter in MATLAB with project root as CWD.

%% Section 0: Setup
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));
projectRoot = pwd();

logSection("RP01", "Section 0: Setup", "ESOL QSAR + L05 (M-REPRO-AUDIT B2)");
emk.setup.initPython();
molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;
snap = emk.setup.snapshot();
lockTemplatePath = fullfile(thisDir, "lock_template.json");
if isfile(lockTemplatePath)
    lockTemplate = jsondecode(fileread(lockTemplatePath));
else
    lockTemplate = struct();
    logWarn("lock_template.json not found; RF03 criteria will fall back to hardcoded defaults.");
end
rng(42, "twister");  % initialize RNG state for Sections 1-4 (no MATLAB randomness currently consumed here); CV fold assignment is independently re-seeded at the top of Section 5

%% Section 1: Load ESOL Dataset
logSection("RP01", "Section 1: Load ESOL Dataset", "ESOL QSAR + L05");

tbl    = emk.dataset.esol();
nTotal = height(tbl);
logInfo("Loaded %d molecules (ESOL / Delaney 2004)", nTotal);
logInfo("logS range: [%.2f, %.2f] log mol/L", min(tbl.logS), max(tbl.logS));

csvFile     = fullfile(projectRoot, "data", "benchmark", "esol.csv");
datasetHash = computeFileHash_(csvFile);
logInfo("Dataset SHA-256: %s", datasetHash);

%% Section 2: Parse SMILES & Compute Descriptors
logSection("RP01", "Section 2: Parse SMILES & Compute Descriptors", "ESOL QSAR + L05");

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
nValidSmiles = sum(validMask);   % SMILES-valid count before descriptor NaN exclusion
nValid       = nValidSmiles;
validMols    = mols(validMask);
yAll         = tbl.logS(validMask);

allDesc = ["LogP", "MolWt", "NumRotatableBonds", "HeavyAtomCount", ...
           "TPSA", "NumHDonors", "NumHAcceptors", "FractionCSP3"];
descTbl = emk.descriptor.batchCalculate(validMols, allDesc);
[descTbl.AromaticProportion, descTbl.QED] = calcAroPropAndQED_(validMols, descTbl.HeavyAtomCount);

% Exclude molecules where any descriptor is NaN (covers all batchCalculate outputs, AroProp, and QED).
nanDescMask = any(ismissing(descTbl), 2);
if any(nanDescMask)
    logWarn("Excluding %d molecule(s) with NaN in any descriptor from model fitting.", sum(nanDescMask));
    validIndices = find(validMask);
    validMask(validIndices(nanDescMask)) = false;
    descTbl   = descTbl(~nanDescMask, :);
    validMols = validMols(~nanDescMask);
    yAll      = yAll(~nanDescMask);
    nValid    = sum(validMask);
end

logInfo("Parsed %d / %d SMILES.  Descriptors: %d used (-%d NaN excluded).", ...
    nValidSmiles, nTotal, nValid, nValidSmiles - nValid);

%% Section 3: Build Models + Extract Coefficients
logSection("RP01", "Section 3: Build Models + Coefficients", "ESOL QSAR + L05");

y  = yAll;
Xa = [descTbl.LogP, descTbl.MolWt, descTbl.NumRotatableBonds, descTbl.AromaticProportion];
Xb = [Xa, descTbl.TPSA, descTbl.NumHDonors, descTbl.NumHAcceptors, ...
      descTbl.FractionCSP3, descTbl.QED];

featsA = {'LogP','MolWt','RotBonds','AroProp'};
featsB = [featsA, {'TPSA','HBD','HBA','FracCSP3','QED'}];

mdlA = fitlm(Xa, y, "VarNames", [featsA, {'logS'}]);
mdlB = fitlm(Xb, y, "VarNames", [featsB, {'logS'}]);

% Training RMSE = sqrt(RSS/n) -- biased estimator matching README "RMSE (full-dataset fit)".
% sqrt(mdl.MSE) uses RSS/(n-p-1) (unbiased) and is smaller; do not use it here.
rmseTrainA = sqrt(mean(mdlA.Residuals.Raw.^2));
rmseTrainB = sqrt(mean(mdlB.Residuals.Raw.^2));

logInfo("Model A (4-feature) full: RMSE=%.4f  R^2=%.4f", rmseTrainA, mdlA.Rsquared.Ordinary);
logInfo("Model B (9-feature) full: RMSE=%.4f  R^2=%.4f", rmseTrainB, mdlB.Rsquared.Ordinary);

% Extract coefficients for TPSA sign check (B2)
coefB     = mdlB.Coefficients;
coefNames = coefB.Properties.RowNames;
tpsaIdx   = find(strcmp(coefNames, "TPSA"));
tpsaCoef  = coefB.Estimate(tpsaIdx);
tpsaPval  = coefB.pValue(tpsaIdx);
tpsaCorr  = corr(descTbl.TPSA, y);   % marginal Pearson r

logInfo("--- TPSA sign verification (B2) ---");
logInfo("  Marginal Pearson r(TPSA, logS) = %.4f  [expected sign: positive]", tpsaCorr);
logInfo("  Model B partial coef of TPSA   = %+.4f  (p=%.2e)", tpsaCoef, tpsaPval);
if tpsaCoef > 0
    logInfo("  => TPSA partial coef is POSITIVE: consistent with marginal correlation.");
else
    logInfo("  => TPSA partial coef is NEGATIVE: reversal in multivariate context.");
    logInfo("     Likely cause: TPSA is correlated with LogP (r=%.3f); LogP already", ...
        corr(descTbl.TPSA, descTbl.LogP));
    logInfo("     absorbs most of TPSA's polarity signal -- TPSA partial effect becomes residual.");
    logInfo("     TPSA-HBD r=%.3f  TPSA-HBA r=%.3f", ...
        corr(descTbl.TPSA, descTbl.NumHDonors), ...
        corr(descTbl.TPSA, descTbl.NumHAcceptors));
end

%% Section 4: VIF Analysis (B2)
logSection("RP01", "Section 4: VIF Analysis", "ESOL QSAR + L05");

vifA = computeVIF_(Xa);
vifB = computeVIF_(Xb);

logInfo("--- VIF: Model A (4-feature) ---");
for k = 1:numel(featsA)
    flag = "";
    if vifA(k) > 10; flag = " *** HIGH";
    elseif vifA(k) > 5; flag = " * moderate"; end
    logInfo("  %-12s  VIF=%.2f%s", featsA{k}, vifA(k), flag);
end

logInfo("--- VIF: Model B (9-feature) ---");
for k = 1:numel(featsB)
    flag = "";
    if vifB(k) > 10; flag = " *** HIGH";
    elseif vifB(k) > 5; flag = " * moderate"; end
    logInfo("  %-12s  VIF=%.2f%s", featsB{k}, vifB(k), flag);
end

highVIF = featsB(vifB > 10);
modVIF  = featsB(vifB > 5 & vifB <= 10);
logInfo("High VIF (>10): %s", strjoin(highVIF, ", "));
logInfo("Moderate VIF (5-10): %s", strjoin(modVIF, ", "));

%% Section 5: 5-Fold Cross-Validation
logSection("RP01", "Section 5: 5-Fold Cross-Validation", "ESOL QSAR + L05");

rng(42, "twister");  % re-seed immediately before cvpartition; guards against RNG consumed by batchCalculate/pyrun above and any future emk.* internals
cv = cvpartition(nValid, "KFold", 5);

[rmseCVa, r2CVa, rmsesA, r2sA] = runCV_(Xa, y, cv);
[rmseCVb, r2CVb, rmsesB, r2sB] = runCV_(Xb, y, cv);

logInfo("Model A 5-fold CV: RMSE=%.4f +/- %.4f  R^2=%.4f +/- %.4f", ...
    rmseCVa, std(rmsesA), r2CVa, std(r2sA));
logInfo("Model B 5-fold CV: RMSE=%.4f +/- %.4f  R^2=%.4f +/- %.4f", ...
    rmseCVb, std(rmsesB), r2CVb, std(r2sB));
logInfo("Per-fold RMSE -- A: %s", sprintf("%.4f ", rmsesA'));
logInfo("Per-fold RMSE -- B: %s", sprintf("%.4f ", rmsesB'));
logInfo("Per-fold R^2  -- A: %s  (min=%.4f)", sprintf("%.4f ", r2sA'), min(r2sA));
logInfo("Per-fold R^2  -- B: %s  (min=%.4f)", sprintf("%.4f ", r2sB'), min(r2sB));

%% Section 6: Paired t-test Model A vs B (B2)
logSection("RP01", "Section 6: Paired t-test A vs B", "ESOL QSAR + L05");

% H1 (one-sided): Model B has lower RMSE => test A-B > 0 (Tail="right")
[~, p_rmse, ~, st_rmse] = ttest(rmsesA, rmsesB, "Tail", "right");
% H1 (one-sided): Model B has higher R^2 => test A-B < 0 (Tail="left")
[~, p_r2, ~, st_r2]     = ttest(r2sA, r2sB, "Tail", "left");

deltaRmse = rmseCVb - rmseCVa;
deltaR2   = r2CVb   - r2CVa;
logInfo("--- Paired t-test: Model A vs B (n_folds=5) ---");
logInfo("  RMSE: A=%.4f  B=%.4f  delta=%.4f", rmseCVa, rmseCVb, deltaRmse);
logInfo("  t(%d)=%.3f  p(one-sided, B better)=%.3f", st_rmse.df, st_rmse.tstat, p_rmse);
logInfo("  R^2:  A=%.4f  B=%.4f  delta=%.4f", r2CVa, r2CVb, deltaR2);
logInfo("  t(%d)=%.3f  p(one-sided, B better)=%.3f", st_r2.df, st_r2.tstat, p_r2);

if p_rmse < 0.05
    logInfo("  => L05 RMSE improvement is SIGNIFICANT (p=%.3f < 0.05)", p_rmse);
else
    logInfo("  => L05 RMSE improvement is NOT significant (p=%.3f, n=5 folds)", p_rmse);
end
if p_r2 < 0.05
    logInfo("  => L05 R^2  improvement is SIGNIFICANT (p=%.3f < 0.05)", p_r2);
else
    logInfo("  => L05 R^2  improvement is NOT significant (p=%.3f, n=5 folds)", p_r2);
end
logInfo("  NOTE: df=4 (n=5 folds) yields very low power -- non-significant result");
logInfo("        does not establish equivalence; effect may exist but be undetectable.");

%% Section 7: RF03 Verification
logSection("RP01", "Section 7: RF03 Verification", "ESOL QSAR + L05");

if isfield(lockTemplate, "rf03_criteria")
    rf03crit = lockTemplate.rf03_criteria;
else
    rf03crit.rmse_cv = struct("upper", 1.20);
    rf03crit.r2_cv   = struct("lower", 0.75);
end
metA.rmse_cv = rmseCVa;  metA.r2_cv = r2CVa;
metB.rmse_cv = rmseCVb;  metB.r2_cv = r2CVb;
resA = emk.repro.verify(metA, rf03crit);
resB = emk.repro.verify(metB, rf03crit);

logInfo("==> Model A: %s  |  Model B: %s", statusStr_(resA.pass), statusStr_(resB.pass));
disp(resA.report);
disp(resB.report);

%% Section 8: Save Results
logSection("RP01", "Section 8: Save Results", "ESOL QSAR + L05");

runDir = makeRunDir("Prefix", "rp01_esol");
% N3: makeRunDir returns a relative path. Guard against absolute paths.
if startsWith(runDir, '/') || (numel(runDir) >= 2 && runDir(2) == ':')
    absRunDir = char(runDir);
else
    absRunDir = char(fullfile(projectRoot, runDir));
end

yPredA = predict(mdlA, Xa);
yPredB = predict(mdlB, Xb);
outTbl = tbl(validMask, ["SMILES","Name","logS"]);
outTbl.logS_predA = yPredA;
outTbl.logS_predB = yPredB;
outTbl.residual_A = y - yPredA;
outTbl.residual_B = y - yPredB;
outTbl = [outTbl, descTbl];
writetable(outTbl, fullfile(absRunDir, "predictions.csv"));

% Axis range based on measured y: keeps the diagonal (y=x) fully visible.
% Predicted outliers appear as clipped points rather than compressing the view.
loAll = floor(min(y) - 0.3);
hiAll = ceil(max(y) + 0.3);

fig = figure("Name", "RP01 ESOL: Model A vs B");
try
    subplot(1,2,1);
    scatter(y, yPredA, 20, "filled", "MarkerFaceAlpha", 0.5);
    hold on;
    plot([loAll,hiAll],[loAll,hiAll],"k--","LineWidth",1.2);
    xlim([loAll,hiAll]); ylim([loAll,hiAll]);
    xlabel("Measured logS");  ylabel("Predicted logS");
    title(sprintf("Model A (4-feat)  RMSE=%.3f  R^2=%.3f", rmseTrainA, mdlA.Rsquared.Ordinary));
    subplot(1,2,2);
    scatter(y, yPredB, 20, "filled", "MarkerFaceAlpha", 0.5, "MarkerFaceColor","#0072BD");
    hold on;
    plot([loAll,hiAll],[loAll,hiAll],"k--","LineWidth",1.2);
    xlim([loAll,hiAll]); ylim([loAll,hiAll]);
    xlabel("Measured logS");  ylabel("Predicted logS");
    title(sprintf("Model B (9-feat)  RMSE=%.3f  R^2=%.3f", rmseTrainB, mdlB.Rsquared.Ordinary));
    saveas(fig, fullfile(absRunDir, "predicted_vs_actual.png"));
finally
    close(fig);
end

% Coefficient table (for TPSA sign archiving)
termCol  = string(coefNames);
coefBmat = table(termCol, coefB.Estimate, coefB.SE, coefB.tStat, coefB.pValue, ...
    'VariableNames', {'term','estimate','se','tStat','pValue'});
writetable(coefBmat, fullfile(absRunDir, "model_b_coefficients.csv"));

metrics = struct( ...
    "model_a_rmse_train",    rmseTrainA, ...
    "model_a_r2_train",      mdlA.Rsquared.Ordinary, ...
    "model_a_rmse_cv",       rmseCVa, ...
    "model_a_rmse_cv_std",   std(rmsesA), ...
    "model_a_r2_cv",         r2CVa, ...
    "model_a_r2_cv_std",     std(r2sA), ...
    "model_a_fold_rmse",     rmsesA(:)', ...
    "model_a_fold_r2",       r2sA(:)', ...
    "model_b_rmse_train",    rmseTrainB, ...
    "model_b_r2_train",      mdlB.Rsquared.Ordinary, ...
    "model_b_rmse_cv",       rmseCVb, ...
    "model_b_rmse_cv_std",   std(rmsesB), ...
    "model_b_r2_cv",         r2CVb, ...
    "model_b_r2_cv_std",     std(r2sB), ...
    "model_b_fold_rmse",     rmsesB(:)', ...
    "model_b_fold_r2",       r2sB(:)', ...
    "l05_delta_rmse",        deltaRmse, ...
    "l05_delta_r2",          deltaR2, ...
    "dataset_sha256",        datasetHash, ...
    "paired_ttest_rmse", struct( ...
        "t",       st_rmse.tstat, ...
        "df",      st_rmse.df, ...
        "p_one",   p_rmse, ...
        "significant", p_rmse < 0.05), ...
    "paired_ttest_r2", struct( ...
        "t",       st_r2.tstat, ...
        "df",      st_r2.df, ...
        "p_one",   p_r2, ...
        "significant", p_r2 < 0.05), ...
    "vif_model_a", struct( ...
        "features", {featsA}, ...
        "values",   vifA(:)'), ...
    "vif_model_b", struct( ...
        "features", {featsB}, ...
        "values",   vifB(:)'), ...
    "tpsa_b2", struct( ...
        "marginal_r",  tpsaCorr, ...
        "partial_coef",tpsaCoef, ...
        "p_value",     tpsaPval, ...
        "tpsa_logp_r", corr(descTbl.TPSA, descTbl.LogP), ...
        "tpsa_hbd_r",  corr(descTbl.TPSA, descTbl.NumHDonors), ...
        "tpsa_hba_r",  corr(descTbl.TPSA, descTbl.NumHAcceptors)), ...
    "n_molecules",           nValid, ...
    "rf03_criteria",         rf03crit, ...
    "rf03_model_a_pass",     resA.pass, ...
    "rf03_model_b_pass",     resB.pass, ...
    "rf03_pass",             resA.pass && resB.pass);

fid = fopen(fullfile(absRunDir, "metrics.json"), "w");
try
    fprintf(fid, "%s\n", jsonencode(metrics, "PrettyPrint", true));
finally
    fclose(fid);
end

snap.run_date  = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir   = runDir;
snap.rf03_pass = resA.pass && resB.pass;
snap.dataset_sha256 = datasetHash;
emk.setup.lockfile(snap, fullfile(absRunDir, "lock_snapshot.json"));

logInfo("Results saved: %s", runDir);
logInfo("==> RP01 complete.");

% ===========================================================================
% Local helper functions
% ===========================================================================

function [aroProp, qedScores] = calcAroPropAndQED_(mols, heavyAtomCounts)
    n = numel(mols);
    aros = NaN(n, 1);
    qeds = NaN(n, 1);
    pyMolList = py.list();
    for i = 1:n; pyMolList.append(mols{i}); end
    try
        [pyAros, pyQeds] = pyrun([ ...
            "from rdkit.Chem import QED as _QED;" ...
            "aros = [float(sum(a.GetIsAromatic() for a in m.GetAtoms())) for m in mols];" ...
            "qeds = [_QED.qed(m) for m in mols]"], ...
            ["aros","qeds"], mols=pyMolList);
        batchAros = double(py.array.array("d", pyAros));
        batchQeds = double(py.array.array("d", pyQeds));
        assert(numel(batchAros) == n && numel(batchQeds) == n, ...
            "pyrun output length mismatch: aros=%d qeds=%d expected=%d", ...
            numel(batchAros), numel(batchQeds), n);
        aros = batchAros(:);
        qeds = batchQeds(:);
    catch ME
        logWarn("calcAroPropAndQED_: batch pyrun failed (%s); retrying per-molecule.", ME.message);
        for i = 1:n
            try
                [pyAro, pyQed] = pyrun([ ...
                    "from rdkit.Chem import QED as _QED;" ...
                    "aro = float(sum(a.GetIsAromatic() for a in m.GetAtoms()));" ...
                    "qed = _QED.qed(m)"], ...
                    ["aro","qed"], m=mols{i});
                aros(i) = double(pyAro);
                qeds(i) = double(pyQed);
            catch ME2
                logWarn("calcAroPropAndQED_: molecule %d failed (%s); NaN assigned.", i, ME2.message);
            end
        end
    end
    aroProp = NaN(n, 1);
    for i = 1:n
        nH = heavyAtomCounts(i);
        if ~isnan(nH) && nH > 0 && ~isnan(aros(i))
            aroProp(i) = aros(i) / nH;
        end
    end
    qedScores = qeds(:);
end

function vifVals = computeVIF_(X)
% VIF_j = 1 / (1 - R^2_j), where R^2_j is from regressing
% column j against all other columns.  VIF > 5 = moderate,
% VIF > 10 = severe multicollinearity.
    p = size(X, 2);
    vifVals = zeros(p, 1);
    for j = 1:p
        otherCols = [X(:, 1:j-1), X(:, j+1:end)];
        mdlJ = fitlm(otherCols, X(:, j));
        if mdlJ.Rsquared.Ordinary > 0.9999
            logWarn("VIF: feature %d is nearly perfectly collinear; VIF unreliable.", j);
        end
        denom = max(1 - mdlJ.Rsquared.Ordinary, 1e-10);
        vifVals(j) = 1 / denom;
    end
end

function [rmseCV, r2CV, rmses, r2s] = runCV_(X, y, cv)
    K    = cv.NumTestSets;
    rmses = zeros(K, 1);
    r2s   = zeros(K, 1);
    for fold = 1:K
        trainIdx = training(cv, fold);
        testIdx  = test(cv, fold);
        mdl      = fitlm(X(trainIdx,:), y(trainIdx));
        yPred    = predict(mdl, X(testIdx,:));
        res      = y(testIdx) - yPred;
        rmses(fold) = sqrt(mean(res.^2));
        ssTot       = sum((y(testIdx) - mean(y)).^2);  % global mean baseline (OOS R^2); y is the full label vector passed as arg 2
        r2s(fold)   = 1 - sum(res.^2) / ssTot;
    end
    rmseCV = mean(rmses);
    r2CV   = mean(r2s);
end

function s = statusStr_(tf)
    if tf; s = "PASS"; else; s = "FAIL"; end
end

function h = computeFileHash_(fpath)
% Compute SHA-256 hex digest of the dataset file for RF02 traceability.
    fpathChar = char(fpath);
    if ispc()
        fpathSafe = strrep(fpathChar, "'", "''");
        cmd = sprintf( ...
            'powershell -Command "(Get-FileHash -LiteralPath ''%s'' -Algorithm SHA256).Hash"', ...
            fpathSafe);
        [status, out] = system(cmd);
        if status == 0
            h = string(upper(strtrim(out)));
            return;
        end
    else
        fpathSafe = strrep(fpathChar, "'", "''");
        [status, out] = system(sprintf("sha256sum '%s'", fpathSafe));
        if status ~= 0
            [status, out] = system(sprintf("shasum -a 256 '%s'", fpathSafe));
        end
        if status == 0
            parts = strsplit(strtrim(out));
            if ~isempty(parts)
                h = string(upper(parts{1}));
                return;
            end
        end
    end
    warning("emk:rp01:computeFileHash:failed", ...
        "SHA-256 hash unavailable for '%s' (platform: %s).", fpathChar, computer());
    h = "unavailable";
end
