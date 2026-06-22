% rp02_bbbp.m  RP02: BBBP Classification (sklearn LR, nested CV)
%
% Canonical baseline using sklearn LogisticRegression (lbfgs, nested CV).
% Note: original fitclinear yielded AUC ~0.019 below sklearn lbfgs on
% identical fold splits (M-REPRO-AUDIT A1, 2026-06-21).
%
%   Model:  sklearn LogisticRegression (L2 ridge, lbfgs, C by inner 3-fold CV)
%   CV:     Outer 5-fold stratified x inner 3-fold for C selection (nested CV)
%   RF03:   ROC-AUC CV >= 0.85
%
%   Output: result/runs/rp02_bbbp_<ts>/
%     metrics.json            -- auc_cv, best C per fold, RF03 pass/fail
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
logInfo("Nested CV AUC = %.4f +/- %.4f", aucCV, res.auc_std);
logInfo("Best C per fold: %s", mat2str([res.fold_results.best_C], 3));

%% Section 3: RF03 Verification
logSection("RP02", "Section 3: RF03 Verification", "BBBP sklearn Baseline");

% L1: named constant so that if RP03 is re-run the value can be updated here
GCN_AUC_RP03 = 0.9151;
% L5: 0.8826 is the historical fitclinear baseline retained for comparison.
% It appears in metrics.json under comparison.rp02_fitclinear_auc but is NOT
% used as an RF03 criterion -- aucCV (nested CV) is the sole pass/fail signal.
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

%% Section 4: ROC curve (pseudo -- all test-fold predictions concatenated)
logSection("RP02", "Section 4: ROC Curve", "BBBP sklearn Baseline");

runDir    = makeRunDir("Prefix", "rp02_bbbp");
absRunDir = char(fullfile(pwd(), runDir));

yTrue = res.y_true;
yProb = res.y_prob;
% M2: pseudoAuc is the AUC of the concatenated test-fold probability scores.
% It approximates full-data performance but is NOT bias-corrected; use aucCV
% for all quantitative reporting and RF03 decisions.
[rocX, rocY, ~, pseudoAuc] = perfcurve(yTrue, yProb, 1);
fig = figure("Name", "RP02 BBBP: sklearn LR ROC (pseudo hold-out)");
plot(rocX, rocY, "b-", "LineWidth", 1.8);
hold on;
plot([0 1], [0 1], "k--", "LineWidth", 1.0);
xlabel("False Positive Rate");
ylabel("True Positive Rate");
title(sprintf("BBBP LR+ECFP4 (sklearn)  CV=%.3f  pseudo-ROC=%.3f (5-fold)", aucCV, pseudoAuc));
legend(sprintf("sklearn LR ECFP4 (pseudo AUC=%.3f)", pseudoAuc), "Random", "Location", "SouthEast");
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
function s = statusStr_(tf)
    if tf; s = "PASS"; else; s = "FAIL"; end
end
