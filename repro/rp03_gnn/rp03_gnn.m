% rp03_gnn.m  RP03: GCN BBBP Classification
%
% 3-layer GCNConv on BBBP blood-brain barrier permeability.
% Implements leak-free CV: inner 80/20 val split for early stopping;
% test fold evaluated once only after training completes.
% Outer fold indices aligned with RP02 for paired comparison.
%
%   Comparison baseline: RP02 sklearn LR, AUC=0.9118
%   RF03: ROC-AUC CV (5-fold) >= 0.85
%
%   Output: result/runs/rp03_gnn_<ts>/
%     metrics.json, learning_curves.png, fold_auc_comparison.png
%     learning_curves.csv, lock_snapshot.json
%
% Run: Ctrl+Enter in MATLAB with project root as CWD.
% NOTE: GCN training takes ~8-15 min on CPU (5 folds x 150 epochs).

%% Section 0: Setup
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("RP03", "Section 0: Setup", "GCN BBBP (M-REPRO-AUDIT A3)");
emk.setup.initPython();
molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;
snap = emk.setup.snapshot();

%% Section 1: Resolve paths + A2 fold indices
logSection("RP03", "Section 1: Resolve Paths + A2 fold indices", ...
    "GCN BBBP (M-REPRO-AUDIT A3)");

tbl = emk.dataset.bbbp();
logInfo("BBBP: %d total molecules in table", height(tbl));

csvPath    = fullfile(pwd(), "data", "benchmark", "bbbp.csv");
coreScript = fullfile(thisDir, "rp03_gnn_core.py");

% Auto-discover latest rp02_bbbp run for fold index alignment
rp02Runs = dir(fullfile("result", "runs", "*rp02_bbbp*"));
if isempty(rp02Runs)
    logWarn("No rp02_bbbp run found. Using internal StratifiedKFold (fold alignment not applied).");
    foldIdxPath = "";
else
    [~, sortIdx]  = sort([rp02Runs.datenum], "descend");
    latestDir  = fullfile("result", "runs", rp02Runs(sortIdx(1)).name);
    foldIdxPath   = fullfile(pwd(), latestDir, "outer_fold_indices.json");
    if isfile(foldIdxPath)
        logInfo("A2 fold indices: %s", foldIdxPath);
    else
        logWarn("outer_fold_indices.json not found in %s. Using internal KFold.", latestDir);
        foldIdxPath = "";
    end
end

%% Section 2: GCN Training with leak-free CV (Python)
logSection("RP03", "Section 2: GCN Nested CV (leak-free)", ...
    "GCN BBBP (M-REPRO-AUDIT A3)");
logInfo("Fix: inner 80/20 val split for early stopping; test fold evaluated once only.");
logInfo("Expected runtime: 8-15 min on CPU.");

pyResult = pyrun([ ...
    "exec(open(core_script).read()); " ...
    "result_json = run_rp03(csv_path, fold_indices_path=fold_idx_path)"], ...
    "result_json", ...
    core_script=coreScript, ...
    csv_path=csvPath, ...
    fold_idx_path=foldIdxPath);

res     = jsondecode(char(string(pyResult)));
aucCV   = res.auc_cv;
aucStd  = res.auc_cv_std;
foldAUC = res.fold_aucs;

logInfo("fold_source = %s", res.fold_source);
logInfo("");
for k = 1:numel(res.fold_curves)
    fc = res.fold_curves(k);
    logInfo("  fold %d: test_AUC=%.4f  best_val_AUC=%.4f  epochs=%d  (tr=%d val=%d te=%d)", ...
        fc.fold, fc.test_auc, fc.best_val_auc, fc.n_epochs, ...
        fc.train_size, fc.val_size, fc.test_size);
end
logInfo("");
logInfo("GCN-rev 5-fold CV: AUC = %.4f +/- %.4f", aucCV, aucStd);

%% Section 3: Comparison report
logSection("RP03", "Section 3: Comparison Report", "GCN BBBP (M-REPRO-AUDIT A3)");

auc_rp02_orig = 0.8826;   % original fitclinear (biased)
auc_rp02_rev  = 0.9118;   % sklearn nested CV (fair, M-REPRO-AUDIT A1)
auc_rp03_orig = 0.9151;   % original GCN (leaked)
auc_rp03_rev  = aucCV;    % this run

logInfo("=== RP03 AUDIT COMPARISON ===");
logInfo("  RP02 original  (fitclinear, biased)   : %.4f", auc_rp02_orig);
logInfo("  RP02           (sklearn nested CV)     : %.4f  [fair baseline]", auc_rp02_rev);
logInfo("  RP03 original  (GCN, test leak)        : %.4f  [reported]", auc_rp03_orig);
logInfo("  RP03       (GCN, leak-fixed)       : %.4f  [this run]", auc_rp03_rev);
logInfo("");
logInfo("  Original GCN advantage (biased): %+.4f", auc_rp03_orig - auc_rp02_orig);
logInfo("  Revised  GCN advantage (fair)  : %+.4f", auc_rp03_rev  - auc_rp02_rev);
if abs(auc_rp03_rev - auc_rp02_rev) < aucStd
    logInfo("  >> GCN advantage NOT significant (gap < 1 std=%.4f)", aucStd);
else
    logInfo("  >> GCN advantage: %.1f std", abs(auc_rp03_rev - auc_rp02_rev) / aucStd);
end

%% Section 4: RF03 Verification
logSection("RP03", "Section 4: RF03 Verification", "GCN BBBP (M-REPRO-AUDIT A3)");

rf03crit = struct("auc_cv", struct("lower", 0.85));
metRP03r = struct("auc_cv", aucCV);
resRP03r = emk.repro.verify(metRP03r, rf03crit);
logInfo("==> ROC-AUC CV = %.4f (>= 0.85): %s", aucCV, statusStr_(resRP03r.pass));
disp(resRP03r.report);

%% Section 5: Visualization
logSection("RP03", "Section 5: Visualization", "GCN BBBP (M-REPRO-AUDIT A3)");

if exist("fig1", "var") && isvalid(fig1); close(fig1); end
if exist("fig2", "var") && isvalid(fig2); close(fig2); end

nEpRun  = res.n_epochs_run;
avgLoss = res.avg_train_loss;
avgVAUC = res.avg_val_auc;
epochs  = (1:nEpRun)';

fig1 = figure("Name", "RP03 GCN: Learning Curves (val AUC -- no leak)");
set(fig1, "Position", [100 100 800 380]);
yyaxis left;
plot(epochs, avgLoss, "b-", "LineWidth", 1.5);
ylabel("Mean Train Loss (BCE)");
yyaxis right;
plot(epochs, avgVAUC, "r-", "LineWidth", 1.5);
yline(0.85,         "k--", "RF03>=0.85",     "LineWidth", 1.0, "LabelHorizontalAlignment", "left");
yline(auc_rp02_rev, "m--", "RP02 LR",    "LineWidth", 1.0, "LabelHorizontalAlignment", "left");
ylabel("Mean Val ROC-AUC (inner, no leak)");
aucLoV = min([min(avgVAUC(:)), 0.85, auc_rp02_rev]);
aucHiV = max(avgVAUC(:));
aucMargin = max(0.05, (aucHiV - aucLoV) * 0.15);
ylim([max(0.0, aucLoV - aucMargin), min(1.0, aucHiV + aucMargin)]);
if min(avgVAUC(:)) < 0.5
    logWarn("Val AUC dropped below 0.5 (min=%.4f). Check for training failure.", min(avgVAUC(:)));
end
xlabel("Epoch");
title(sprintf("RP03 GCN (leak-fixed): Learning Curves  AUC CV=%.4f", aucCV), ...
    sprintf("avg curves show first %d epochs (truncated to shortest fold)", nEpRun));
legend(["Train Loss", "Val AUC (inner)", "RF03>=0.85", "RP02-rev"], "Location", "east");
grid on;

fig2 = figure("Name", "RP03 GCN: Fold AUC Comparison");
set(fig2, "Position", [920 100 500 380]);
bar(1:numel(foldAUC), foldAUC, "FaceColor", [0.2 0.5 0.85], "EdgeColor", "none");
hold on;
yline(auc_rp02_rev,  "m--", sprintf("RP02 LR (%.4f)", auc_rp02_rev), "LineWidth", 1.5);
yline(auc_rp02_orig, "g:",  sprintf("RP02 orig  (%.4f)", auc_rp02_orig), "LineWidth", 1.2);
yline(aucCV,         "r-",  sprintf("GCN-rev CV (%.4f)", aucCV),         "LineWidth", 1.5);
hold off;
allFoldVals = [foldAUC(:); auc_rp02_rev; auc_rp02_orig; aucCV];
fig2Margin  = max(0.05, (max(allFoldVals) - min(allFoldVals)) * 0.15);
ylim([max(0.0, min(allFoldVals) - fig2Margin), min(1.0, max(allFoldVals) + fig2Margin)]);
xlabel("Fold"); ylabel("Test ROC-AUC");
xticks(1:numel(foldAUC));
xticklabels(arrayfun(@(k) sprintf("Fold %d", k), 1:numel(foldAUC), "UniformOutput", false));
title("RP03 GCN vs RP02-rev LR+ECFP4: Per-Fold AUC (BBBP)");
grid on; box off;

logInfo("Figures ready (fig1, fig2). Run Section 6 to save all results.");

%% Section 6: Save Results
logSection("RP03", "Section 6: Save Results", "GCN BBBP (M-REPRO-AUDIT A3)");

runDir    = makeRunDir("Prefix", "rp03_gnn");
absRunDir = char(fullfile(pwd(), runDir));

if exist("fig1", "var") && isvalid(fig1)
    saveas(fig1, fullfile(absRunDir, "learning_curves.png"));
    close(fig1);
end
if exist("fig2", "var") && isvalid(fig2)
    saveas(fig2, fullfile(absRunDir, "fold_auc_comparison.png"));
    close(fig2);
end
logInfo("Figures saved.");

if isfield(res, "learning_curve_note")
    lcNote_ = string(res.learning_curve_note);
else
    lcNote_ = sprintf("avg curves truncated to %d epochs (shortest fold).", nEpRun);
end

metrics = struct( ...
    "auc_cv",              aucCV, ...
    "auc_cv_std",          aucStd, ...
    "fold_aucs",           foldAUC, ...
    "fold_source",         res.fold_source, ...
    "n_valid",             res.n_valid, ...
    "n_bbb_pos",           res.n_bbb_pos, ...
    "n_bbb_neg",           res.n_bbb_neg, ...
    "n_epochs_run",        nEpRun, ...
    "learning_curve_note", lcNote_, ...
    "train_size_note",     "RP03 sub-train (~80% of RP02 outer train) differs from RP02 train size. Paired comparison valid on identical outer test folds (A2 alignment) only.", ...
    "audit_comparison", struct( ...
        "rp02_orig_fitclinear",  auc_rp02_orig, ...
        "rp02_rev_sklearn",      auc_rp02_rev, ...
        "rp03_orig_leaked",      auc_rp03_orig, ...
        "rp03_rev_leak_fixed",   aucCV, ...
        "orig_gcn_advantage",    auc_rp03_orig - auc_rp02_orig, ...
        "fair_gcn_advantage",    auc_rp03_rev  - auc_rp02_rev), ...
    "hyperparams",         res.hyperparams, ...
    "rf03_criteria",       rf03crit, ...
    "rf03_pass",           resRP03r.pass);

fid = fopen(fullfile(runDir, "metrics.json"), "w");
fprintf(fid, "%s\n", jsonencode(metrics, "PrettyPrint", true));
fclose(fid);

lcTbl = table(epochs, avgLoss(:), avgVAUC(:), ...
    VariableNames=["Epoch", "AvgTrainLoss", "AvgValAUC"]);
writetable(lcTbl, fullfile(runDir, "learning_curves.csv"));

snap.run_date  = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir   = runDir;
snap.rf03_pass = resRP03r.pass;
emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));

logInfo("Results saved: %s", runDir);
if resRP03r.pass
    logInfo("==> RP03 REPRODUCTION: PASS  (run_dir=%s)", runDir);
else
    logWarn("==> RP03 REPRODUCTION: NEEDS REVIEW  (run_dir=%s)", runDir);
end

% ===========================================================================
function s = statusStr_(tf)
    if tf; s = "PASS"; else; s = "FAIL"; end
end
