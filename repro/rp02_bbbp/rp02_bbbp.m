% rp02_bbbp.m  RP02: BBBP Classification (sklearn LR, nested CV) + R1 fitclinear
%
% Canonical baseline using sklearn LogisticRegression (lbfgs, nested CV).
% Note: original fitclinear yielded AUC ~0.019 below sklearn lbfgs on
% identical fold splits (M-REPRO-AUDIT A1, 2026-06-21).
%
% R1 (M-REPRO-REFINE Phase 1, 2026-06-26):
%   MATLAB fitclinear re-run with Lambda=1/C_opt (C_opt from Section 2 inner CV).
%   Same outer fold splits as sklearn (shared via res.outer_fold_indices).
%   Purpose: Decompose historical AUC gap into regularization and solver components.
%     regularization_gap = AUC(fitclinear,Lambda=1/C_opt) - AUC(fitclinear,Lambda=1/n)
%     solver_gap         = AUC(sklearn,C=C_opt)           - AUC(fitclinear,Lambda=1/C_opt)
%
%   Model:  sklearn LogisticRegression (L2 ridge, lbfgs, C by inner 3-fold CV)
%   CV:     Outer 5-fold stratified x inner 3-fold for C selection (nested CV)
%   RF03:   ROC-AUC CV >= 0.85
%
%   Output: result/runs/rp02_bbbp_<ts>/
%     metrics.json            -- auc_cv, best C per fold, RF03 pass/fail, R1 gap
%     outer_fold_indices.json -- fold indices for RP03/RP04 fold alignment (A2)
%     roc_curve.png           -- pseudo-ROC (all test-fold predictions concatenated)
%     lock_snapshot.json      -- RF02 version lock
%
% Run: Ctrl+Enter in MATLAB with project root as CWD.

%% Section 0: Setup
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("RP02", "Section 0: Setup", "BBBP sklearn Baseline");
emk.setup.initPython();
% L3: warmup call forces the Python engine to start before timed sections.
% Without this, the first real call absorbs engine startup latency (~3-5 s).
molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;
snap = emk.setup.snapshot();

%% Section 1: Ensure BBBP CSV is cached
logSection("RP02", "Section 1: Ensure BBBP CSV cached", "BBBP sklearn Baseline");

% emk.dataset.bbbp() downloads and caches data/benchmark/bbbp.csv.
tbl = emk.dataset.bbbp();
logInfo("RP02r: %d total molecules in table (CSV cached)", height(tbl));
csvPath = fullfile(pwd(), "data", "benchmark", "bbbp.csv");
if ~isfile(csvPath)
    error("emk:rp02rev:csvNotFound", "bbbp.csv not found at: %s", csvPath);
end

%% Section 2: Run Python nested CV (sklearn LR + inner C selection)
logSection("RP02", "Section 2: Python nested CV (sklearn LR)", "BBBP sklearn Baseline");

coreScript = fullfile(thisDir, "rp02_sklearn_core.py");
logInfo("Running nested CV (outer 5-fold x inner 3-fold C selection) ...");
logInfo("C grid: [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0]");

% H2: Python exceptions are caught inside run_rp02_sklearn and returned as
% JSON with success=false. This ensures the error is always surfaced here
% rather than surfacing as a jsondecode failure on the next line.
resultJson = pyrun([ ...
    "exec(open(script_path).read()); " ...
    "result_json = run_rp02_sklearn(csv_path)"], ...
    "result_json", ...
    script_path=coreScript, ...
    csv_path=csvPath);

res = jsondecode(char(string(resultJson)));

% H2: check for Python-side failure before proceeding
if isfield(res, "success") && ~res.success
    error("emk:rp02:pythonError", ...
        "rp02_sklearn_core.py failed:\n%s\n\nTraceback:\n%s", ...
        res.error, res.traceback);
end

logInfo("n_valid=%d  BBB+=%.0f  BBB-=%.0f", res.n_valid, res.n_bbb_pos, res.n_bbb_neg);
logInfo("");

for k = 1:numel(res.fold_results)
    fr = res.fold_results(k);
    logInfo("  fold %d: AUC=%.4f  best_C=%.3g  (inner_AUC=%.4f, train=%d, test=%d)", ...
        fr.fold, fr.test_auc, fr.best_C, fr.inner_auc, fr.train_size, fr.test_size);
end

aucCV = res.auc_mean;
logInfo("");
logInfo("Nested CV AUC = %.4f +/- %.4f [%s]", aucCV, res.auc_std, res.auc_std_definition);
logInfo("Best C per fold: %s", mat2str([res.fold_results.best_C], 3));

%% Section 2b: MATLAB fitclinear with matched Lambda (R1: solver vs regularization)
logSection("RP02", "Section 2b: fitclinear Lambda=1/C_opt (R1)", "BBBP sklearn Baseline");
%
% Uses same outer fold splits as Section 2 (sklearn res.outer_fold_indices).
% Per fold: Lambda = 1/best_C (best_C from inner-CV selection in Section 2).
%   regularization_gap = AUC(fitclinear,Lambda=1/C_opt) - AUC(fitclinear,Lambda=1/n)
%   solver_gap         = AUC(sklearn,C=C_opt)           - AUC(fitclinear,Lambda=1/C_opt)

FITCLINEAR_AUC_HISTORICAL = 0.8826;

% Compute MATLAB-side ECFP4 fingerprints. Reuses tbl from Section 1 to
% avoid re-downloading. Same SMILES-validation pipeline as Python _load_ecfp4.
nTotalR1 = height(tbl);
molsR1   = cell(1, nTotalR1);
vmaskR1  = false(1, nTotalR1);
for i = 1:nTotalR1
    try
        molsR1{i}  = emk.mol.fromSmiles(tbl.SMILES(i));
        vmaskR1(i) = true;
    catch
    end
end
yR1 = double(tbl.BBB(vmaskR1));
XR1 = batchMorganFP_(molsR1(vmaskR1), 2, 2048);
logInfo("R1 MATLAB ECFP4: %d x %d  (BBB+=%.0f, BBB-=%.0f)", ...
    size(XR1,1), size(XR1,2), sum(yR1==1), sum(yR1==0));

if size(XR1,1) ~= res.n_valid
    error("emk:rp02:r1:sizeMismatch", ...
        "MATLAB n_valid=%d vs Python n_valid=%d; check SMILES parse consistency.", ...
        size(XR1,1), res.n_valid);
end

% Per-fold fitclinear with Lambda=1/C_opt
nFoldsR1 = numel(res.fold_results);
aucsR1   = zeros(nFoldsR1, 1);
lambdaR1 = zeros(nFoldsR1, 1);
for k = 1:nFoldsR1
    fr        = res.fold_results(k);
    bestC     = fr.best_C;
    lambdaOpt = 1.0 / bestC;
    lambdaR1(k) = lambdaOpt;

    % Fold indices from Python are 0-based; +1 for MATLAB 1-based indexing
    trIdx = res.outer_fold_indices.(sprintf("fold%d_train", k)) + 1;
    teIdx = res.outer_fold_indices.(sprintf("fold%d_test",  k)) + 1;

    mdl  = fitclinear(XR1(trIdx,:), yR1(trIdx), ...
        "Learner",        "logistic", ...
        "Regularization", "ridge", ...
        "Lambda",         lambdaOpt);
    pIdx = find(mdl.ClassNames == 1, 1);
    [~, sc] = predict(mdl, XR1(teIdx,:));
    [~,~,~, aucK] = perfcurve(yR1(teIdx), sc(:, pIdx), 1);
    aucsR1(k) = aucK;

    logInfo("  fold %d: AUC=%.4f  Lambda=%.4g (1/C=%.3g)  sklearn=%.4f  solver_delta=%+.4f", ...
        k, aucK, lambdaOpt, bestC, fr.test_auc, fr.test_auc - aucK);
end

aucR1mean = mean(aucsR1);
aucR1std  = std(aucsR1);
logInfo("");
logInfo("R1 fitclinear(Lambda=1/C_opt): AUC = %.4f +/- %.4f", aucR1mean, aucR1std);
logInfo("");
regGap    = aucR1mean - FITCLINEAR_AUC_HISTORICAL;
solverGap = aucCV - aucR1mean;
logInfo("--- R1 Gap Decomposition (historical total gap = %+.4f) ---", ...
    aucCV - FITCLINEAR_AUC_HISTORICAL);
logInfo("  regularization_gap = %+.4f  (Lambda=1/C_opt vs Lambda=1/n, same solver)", regGap);
logInfo("  solver_gap         = %+.4f  (sklearn lbfgs vs fitclinear, same Lambda)", solverGap);
logInfo("  sum check          = %+.4f  (= total gap %.4f)", ...
    regGap + solverGap, aucCV - FITCLINEAR_AUC_HISTORICAL);

%% Section 2c: Scaffold 5-fold nested CV (C: random vs scaffold gap)
logSection("RP02", "Section 2c: Scaffold 5-fold CV (C task)", "BBBP sklearn Baseline");
%
% C (M-REPRO-REFINE Phase 2, 2026-06-26):
% Adds Bemis-Murcko scaffold GroupKFold (outer 5-fold x inner 3-fold for C selection).
% Quantifies split_gap = random_auc - scaffold_auc (random-CV advantage from scaffold leakage)
% and wu_gap = scaffold_auc - WU_ET_AL (residual: 5-fold CV vs single test set + dataset diff).
%
% Requires: coreScript (Section 2), csvPath (Section 1), aucCV (Section 2).

WU_ET_AL_AUC = 0.690;  % Wu et al. (2018) Table 4: Circular FP + LogReg, scaffold split

logInfo("Running scaffold 5-fold nested CV (outer GroupKFold x inner GroupKFold C selection) ...");
resultJsonC = pyrun([ ...
    "exec(open(script_path).read()); " ...
    "result_json = run_rp02_scaffold(csv_path)"], ...
    "result_json", ...
    script_path=coreScript, ...
    csv_path=csvPath);

resC = jsondecode(char(string(resultJsonC)));

if isfield(resC, "success") && ~resC.success
    error("emk:rp02:scaffoldError", ...
        "run_rp02_scaffold failed:\n%s\n\nTraceback:\n%s", ...
        resC.error, resC.traceback);
end

if resC.n_scaffold_errors > 0
    logWarn("  %d molecules raised MurckoScaffoldSmiles errors (merged into acyclic group)", ...
        resC.n_scaffold_errors);
end

logInfo("n_valid=%d  BBB+=%.0f  BBB-=%.0f  unique_scaffolds=%d", ...
    resC.n_valid, resC.n_bbb_pos, resC.n_bbb_neg, resC.n_unique_scaffolds);
logInfo("");

for k = 1:numel(resC.fold_results)
    fr = resC.fold_results(k);
    logInfo("  fold %d: AUC=%.4f  best_C=%.3g  (inner_AUC=%.4f, train=%d, test=%d)", ...
        fr.fold, fr.test_auc, fr.best_C, fr.inner_auc, fr.train_size, fr.test_size);
end

aucScaffold    = resC.auc_mean;
aucScaffoldStd = resC.auc_std;
logInfo("");
logInfo("Scaffold 5-fold AUC = %.4f +/- %.4f [%s]", ...
    aucScaffold, aucScaffoldStd, resC.auc_std_definition);
logInfo("");

splitGap = aucCV - aucScaffold;
wuGap    = aucScaffold - WU_ET_AL_AUC;
logInfo("--- C: Random vs Scaffold Gap ---");
logInfo("  Random 5-fold   AUC = %.4f  (Section 2)", aucCV);
logInfo("  Scaffold 5-fold AUC = %.4f  (Section 2c)", aucScaffold);
logInfo("  Wu et al. AUC       = %.3f  (scaffold split, single test set)", WU_ET_AL_AUC);
logInfo("  split_gap = random - scaffold = %+.4f  (split-strategy effect)", splitGap);
logInfo("  wu_gap    = scaffold - wu     = %+.4f  (5-fold vs single-split + dataset diff)", wuGap);
logInfo("  total     = random - wu       = %+.4f", aucCV - WU_ET_AL_AUC);

%% Section 3: RF03 Verification
logSection("RP02", "Section 3: RF03 Verification", "BBBP sklearn Baseline");

% L1: named constant -- MANUAL UPDATE REQUIRED if RP03 is re-run.
% Source: result/runs/rp03_bbbp_gcn_<ts>/metrics.json "auc_test_mean" (2026-06-21).
% Stale value produces silently incorrect comparison output in the log below.
GCN_AUC_RP03 = 0.9151;
% L5: 0.8826 is the historical fitclinear baseline retained for comparison.
% It appears in metrics.json under comparison.rp02_fitclinear_auc but is NOT
% used as an RF03 criterion -- aucCV (nested CV) is the sole pass/fail signal.
% Also defined in Section 2b; kept here for standalone-section safety.
FITCLINEAR_AUC_HISTORICAL = 0.8826;

rf03crit = struct("auc_cv", struct("lower", 0.85));
% M2: aucCV (nested CV mean) is the ONLY metric used for RF03 pass/fail.
% pseudoAuc (computed in Section 4) is a diagnostic view only and must NOT
% be substituted here -- it is computed on the concatenated test-fold
% predictions, which leaks fold-selection information.
metRP02r = struct("auc_cv", aucCV);
resRP02r = emk.repro.verify(metRP02r, rf03crit);

logInfo("==> ROC-AUC CV = %.4f (criterion: >= %.2f): %s", ...
    aucCV, rf03crit.auc_cv.lower, statusStr_(resRP02r.pass));
disp(resRP02r.report);

logInfo("");
logInfo("--- Comparison with original RP02 ---");
logInfo("  RP02 fitclinear AUC  = %.4f  (historical; M-REPRO-AUDIT A1: biased low due to solver)", ...
    FITCLINEAR_AUC_HISTORICAL);
logInfo("  RP02 sklearn AUC     = %.4f  (fair baseline for RP03/RP04 comparison)", aucCV);
logInfo("  Revised GCN advantage    ~ +%.4f  (vs %.4f reported in RP03)", ...
    GCN_AUC_RP03 - aucCV, GCN_AUC_RP03 - FITCLINEAR_AUC_HISTORICAL);

logInfo("");
logInfo("--- R1 Gap Decomposition (Section 2b) ---");
if exist("aucR1mean", "var")
    logInfo("  fitclinear(Lambda=1/C_opt) AUC = %.4f +/- %.4f", aucR1mean, aucR1std);
    logInfo("  regularization_gap = %+.4f  (Lambda mismatch: 1/C_opt vs 1/n)", ...
        aucR1mean - FITCLINEAR_AUC_HISTORICAL);
    logInfo("  solver_gap         = %+.4f  (sklearn lbfgs vs MATLAB fitclinear, same Lambda)", ...
        aucCV - aucR1mean);
else
    logWarn("  R1 results unavailable -- run Section 2b before Section 3.");
end

% C: standalone-section safety (primary definition in Section 2c)
WU_ET_AL_AUC = 0.690;
logInfo("");
logInfo("--- C: Random vs Scaffold Split Gap (Section 2c) ---");
if exist("aucScaffold", "var")
    logInfo("  Scaffold 5-fold AUC             = %.4f +/- %.4f", aucScaffold, aucScaffoldStd);
    logInfo("  split_gap (random - scaffold)   = %+.4f  (split-strategy effect)", ...
        aucCV - aucScaffold);
    logInfo("  wu_gap    (scaffold - wu_et_al) = %+.4f  (5-fold vs single-split + dataset diff)", ...
        aucScaffold - WU_ET_AL_AUC);
else
    logWarn("  Scaffold CV results unavailable -- run Section 2c first.");
end

%% Section 4: ROC curve (pseudo -- all test-fold predictions concatenated)
logSection("RP02", "Section 4: ROC Curve", "BBBP sklearn Baseline");

runDir = makeRunDir("Prefix", "rp02_bbbp");
% N3: makeRunDir returns a CWD-relative path. Guard against absolute paths
% including Unix (/), Windows drive (C:), and UNC (\\server\share).
if startsWith(runDir, '/') || startsWith(runDir, '\\') || ...
        (numel(runDir) >= 2 && runDir(2) == ':')
    absRunDir = char(runDir);
else
    absRunDir = char(fullfile(pwd(), runDir));
end

yTrue = res.y_true;
yProb = res.y_prob;
% M2: pseudoAuc is the AUC of the concatenated test-fold probability scores.
% It approximates full-data performance but is NOT bias-corrected; use aucCV
% for all quantitative reporting and RF03 decisions.
[rocX, rocY, ~, pseudoAuc] = perfcurve(yTrue, yProb, 1);
fig = figure("Name", "RP02 BBBP: sklearn LR ROC (pseudo hold-out)");
plot(rocX, rocY, "b-", "LineWidth", 1.8);
hold on;
% Overlay scaffold pseudo-ROC if Section 2c results are available
if exist("resC", "var") && isfield(resC, "y_true") && exist("aucScaffold", "var")
    [rocXc, rocYc, ~, pseudoAucC] = perfcurve(resC.y_true, resC.y_prob, 1);
    plot(rocXc, rocYc, "r--", "LineWidth", 1.5);
    plot([0 1], [0 1], "k--", "LineWidth", 1.0);
    xlabel("False Positive Rate");
    ylabel("True Positive Rate");
    title(sprintf("BBBP LR+ECFP4  Random CV=%.3f / Scaffold CV=%.3f  (pseudo-ROC)", ...
        aucCV, aucScaffold));
    legend( ...
        sprintf("Random 5-fold (pseudo AUC=%.3f)", pseudoAuc), ...
        sprintf("Scaffold 5-fold (pseudo AUC=%.3f)", pseudoAucC), ...
        "Random chance", "Location", "SouthEast");
else
    plot([0 1], [0 1], "k--", "LineWidth", 1.0);
    xlabel("False Positive Rate");
    ylabel("True Positive Rate");
    title(sprintf("BBBP LR+ECFP4 (sklearn)  CV=%.3f  pseudo-ROC=%.3f (5-fold)", aucCV, pseudoAuc));
    legend(sprintf("sklearn LR ECFP4 (pseudo AUC=%.3f)", pseudoAuc), "Random", "Location", "SouthEast");
end
saveas(fig, fullfile(absRunDir, "roc_curve.png"));
close(fig);
logInfo("Figure saved: roc_curve.png");

%% Section 5: Save results
logSection("RP02", "Section 5: Save Results", "BBBP sklearn Baseline");

% M4: guard against running this section standalone without Section 4
if ~exist("runDir", "var") || isempty(runDir)
    error("emk:rp02:missingRunDir", ...
        "runDir is not defined. Run Section 4 first to create the output directory.");
end

% Outer fold indices (for A2: sharing with RP03/RP04)
foldIdxPath = fullfile(runDir, "outer_fold_indices.json");
fid = fopen(foldIdxPath, "w");
if fid == -1
    error("emk:rp02:fopenFailed", "Cannot open for writing: %s", foldIdxPath);
end
fprintf(fid, "%s\n", jsonencode(res.outer_fold_indices, "PrettyPrint", true));
fclose(fid);
logInfo("Outer fold indices saved (A2 input): %s", foldIdxPath);

% Metrics JSON
% L2: inner_C_grid (per-fold HP selection evidence) is already embedded in
% fold_results; extracted here at top level for quick inspection.
metrics = struct( ...
    "auc_cv",           aucCV, ...
    "auc_cv_std",       res.auc_std, ...
    "auc_cv_std_definition", res.auc_std_definition, ...
    "auc_per_fold",     res.auc_per_fold, ...
    "auc_pseudo_roc",   pseudoAuc, ...
    "best_C_per_fold",  res.best_C_per_fold, ...
    "n_valid",          res.n_valid, ...
    "n_bbb_pos",        res.n_bbb_pos, ...
    "n_bbb_neg",        res.n_bbb_neg, ...
    "ecfp4_radius",     2, ...
    "ecfp4_nbits",      2048, ...
    "solver",           "sklearn_lbfgs", ...
    "cv_outer",         5, ...
    "cv_inner",         3, ...
    "rf03_criteria",    rf03crit, ...
    "rf03_pass",        resRP02r.pass, ...
    "comparison", struct( ...
        "rp02_fitclinear_auc",    FITCLINEAR_AUC_HISTORICAL, ...
        "rp02_rev_sklearn_auc",   aucCV, ...
        "solver_gap_a1",          0.0191, ...
        "fold_gap_a1",            0.0079, ...
        "gcn_auc_rp03",           GCN_AUC_RP03));

% R1 gap decomposition (Section 2b results)
if exist("aucR1mean", "var")
    metrics.r1_fitclinear_matched_lambda = struct( ...
        "auc_cv",             aucR1mean, ...
        "auc_cv_std",         aucR1std, ...
        "auc_cv_std_definition", "sample_std_ddof1", ...
        "auc_per_fold",       aucsR1', ...
        "lambda_per_fold",    lambdaR1', ...
        "regularization_gap", aucR1mean - FITCLINEAR_AUC_HISTORICAL, ...
        "solver_gap",         aucCV - aucR1mean);
    logInfo("R1 gap metrics added to metrics.json.");
else
    logWarn("R1 results not in workspace -- metrics.json will lack r1_* fields.");
end

% C task: scaffold 5-fold CV metrics (Section 2c)
WU_ET_AL_AUC = 0.690;  % standalone-section safety; primary definition in Section 2c
if exist("aucScaffold", "var")
    metrics.scaffold_cv = struct( ...
        "auc_cv",             aucScaffold, ...
        "auc_cv_std",         aucScaffoldStd, ...
        "auc_cv_std_definition", resC.auc_std_definition, ...
        "auc_per_fold",       resC.auc_per_fold, ...
        "n_unique_scaffolds", resC.n_unique_scaffolds, ...
        "split_gap",          aucCV - aucScaffold, ...
        "wu_gap",             aucScaffold - WU_ET_AL_AUC, ...
        "wu_et_al_auc",       WU_ET_AL_AUC);
    logInfo("Scaffold CV metrics added to metrics.json.");
else
    logWarn("Scaffold CV results not in workspace -- metrics.json will lack scaffold_cv fields.");
end

metricsPath = fullfile(runDir, "metrics.json");
fid = fopen(metricsPath, "w");
if fid == -1
    error("emk:rp02:fopenFailed", "Cannot open for writing: %s", metricsPath);
end
fprintf(fid, "%s\n", jsonencode(metrics, "PrettyPrint", true));
fclose(fid);
logInfo("Metrics saved: metrics.json");

% RF02 lock snapshot
snap.run_date  = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir   = runDir;
snap.rf03_pass = resRP02r.pass;
emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));

if resRP02r.pass
    logInfo("==> RP02 REPRODUCTION: PASS  (run_dir=%s)", runDir);
else
    logWarn("==> RP02 REPRODUCTION: NEEDS REVIEW  (run_dir=%s)", runDir);
end
logInfo("NOTE: outer_fold_indices.json available for RP03/RP04 fold alignment.");

% ===========================================================================
function X = batchMorganFP_(mols, radius, nBits)
% batchMorganFP_  Compute ECFP4 fingerprint matrix for a cell array of RDKit mols.
%   mols   : cell array of py.rdkit.Chem.rdchem.Mol objects
%   radius : Morgan radius (2 = ECFP4)
%   nBits  : fingerprint bit length
%   X      : double matrix [n x nBits]
    n = numel(mols);
    pyMolList = py.list();
    for i = 1:n
        pyMolList.append(mols{i});
    end
    pyAllBits = pyrun([ ...
        "from rdkit.Chem import rdFingerprintGenerator as _rfg;" ...
        "gen = _rfg.GetMorganGenerator(radius=int(fp_r), fpSize=int(fp_nb));" ...
        "fps = ''.join(gen.GetFingerprint(m).ToBitString() for m in mols)"], ...
        "fps", mols=pyMolList, fp_r=int32(radius), fp_nb=int32(nBits));
    allBits = char(string(pyAllBits)) == '1';
    X = double(reshape(allBits, nBits, n)');
end

function s = statusStr_(tf)
    if tf; s = "PASS"; else; s = "FAIL"; end
end
