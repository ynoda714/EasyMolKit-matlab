% rp04_chemberta.m  RP04: ChemBERTa BBBP Classification (linear probe)
%
% Frozen CLS embeddings from seyonec/ChemBERTa-zinc-base-v1 (44M params)
% classified with LogisticRegression (linear probe) on BBBP.
%
% Baselines loaded dynamically from latest RP02/RP03 metrics.json.
% Fair baselines (M-REPRO-AUDIT): RP02=0.9118 (sklearn nested CV),
%                                  RP03=0.9038 (GCN leak-fixed).
% Token length validation: max_length=128 affects 1.3% of SMILES (26/2039).
%
%   Run: Ctrl+Enter in MATLAB with project root as CWD.

%% Section 0: Setup
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("RP04", "Section 0: Setup", "ChemBERTa BBBP (M-REPRO-AUDIT B4)");
emk.setup.initPython();
molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;
snap = emk.setup.snapshot();

%% Section 1: Load BBBP Dataset & Resolve Paths
logSection("RP04", "Section 1: Load BBBP Dataset", "ChemBERTa BBBP");

tbl     = emk.dataset.bbbp();
nTotal  = height(tbl);
logInfo("Loaded %d molecules (BBB+: %d, BBB-: %d)", ...
    nTotal, sum(tbl.BBB), sum(~tbl.BBB));

root       = resolveProjectRoot();
csvPath    = fullfile(root, "data", "benchmark", "bbbp.csv");
helperPath = fullfile(thisDir, "rp04_chemberta_core.py");

if ~isfile(csvPath)
    error("emk:rp04rev:csvNotFound", "BBBP CSV not found: %s", csvPath);
end

%% Section 2: Load Fair Baselines from Rev Run Metrics (B4)
logSection("RP04", "Section 2: Load Fair Baselines (B4)", "ChemBERTa BBBP");

% RP02-rev (sklearn nested CV, A1-corrected)
rp02Dirs = dir(fullfile("result", "runs", "*rp02_bbbp*"));
if ~isempty(rp02Dirs)
    [~, si]   = sort({rp02Dirs.name}, "descend");  % name sort: robust to same-second runs
    mRP02     = jsondecode(fileread(fullfile("result", "runs", rp02Dirs(si(1)).name, "metrics.json")));
    auc_rp02  = mRP02.auc_cv;
    std_rp02  = mRP02.auc_cv_std;
    logInfo("RP02-rev loaded from: %s", rp02Dirs(si(1)).name);
else
    logWarn("No rp02_bbbp run found -- using hardcoded fallback (0.9118).");
    auc_rp02 = 0.9118;
    std_rp02 = 0.0075;
end

% RP03-rev (GCN leak-fixed, A3-corrected)
rp03Dirs = dir(fullfile("result", "runs", "*rp03_gnn*"));
if ~isempty(rp03Dirs)
    [~, si]   = sort({rp03Dirs.name}, "descend");  % name sort: robust to same-second runs
    mRP03     = jsondecode(fileread(fullfile("result", "runs", rp03Dirs(si(1)).name, "metrics.json")));
    auc_rp03  = mRP03.auc_cv;
    std_rp03  = mRP03.auc_cv_std;
    logInfo("RP03-rev loaded from: %s", rp03Dirs(si(1)).name);
else
    logWarn("No rp03_gnn run found -- using hardcoded fallback (0.9038).");
    auc_rp03 = 0.9038;
    std_rp03 = 0.0203;
end

logInfo("Fair baselines -- RP02-rev: %.4f +/- %.4f | RP03-rev: %.4f +/- %.4f", ...
    auc_rp02, std_rp02, auc_rp03, std_rp03);

%% Section 3: ChemBERTa Embedding + 5-Fold CV (Python)
logSection("RP04", "Section 3: ChemBERTa Embedding + 5-Fold CV", "ChemBERTa BBBP");
logInfo("Loading ChemBERTa (seyonec/ChemBERTa-zinc-base-v1, 44M params) ...");
logInfo("Extracting CLS embeddings for %d molecules (~2-4 min on CPU)...", nTotal);

try
    pyResult = pyrun( ...
        "exec(open(hp).read()); result_json = run_rp04(cp)", ...
        "result_json", hp=helperPath, cp=csvPath);
catch ME
    error("emk:rp04:pyrunFailed", ...
        "ChemBERTa Python execution failed: %s", ME.message);
end

res     = jsondecode(char(string(pyResult)));
aucCV   = res.auc_cv;
aucStd  = res.auc_cv_std;
foldAUC = res.fold_aucs;
tkStats = res.token_length_stats;

logInfo("ChemBERTa CLS dim: %d  |  Model params: %.1fM", res.hidden_size, res.n_params_M);
logInfo("5-fold CV: AUC = %.4f +/- %.4f  (n=%d)", aucCV, aucStd, res.n_valid);
logInfo("Per-fold AUCs: %s", sprintf("%.4f  ", foldAUC));

%% Section 4: Token Length Validation (B4)
logSection("RP04", "Section 4: Token Length Validation (B4)", "ChemBERTa BBBP");

logInfo("SMILES token length distribution (WITHOUT truncation, n=%d):", res.n_valid);
logInfo("  min=%d  max=%d  mean=%.1f", tkStats.min, tkStats.max, tkStats.mean);
logInfo("  p50=%d  p90=%d  p95=%d  p99=%d", tkStats.p50, tkStats.p90, tkStats.p95, tkStats.p99);
logInfo("  Truncated (> max_length=%d): %d / %d (%.1f%%)", ...
    res.hyperparams.max_len, tkStats.n_truncated, res.n_valid, tkStats.frac_truncated * 100);

if tkStats.frac_truncated > 0.01
    logWarn("More than 1%% of SMILES truncated at max_length=%d -- may affect embeddings.", ...
        res.hyperparams.max_len);
else
    logInfo("max_length=%d validated: truncation affects < 1%% of molecules.", ...
        res.hyperparams.max_len);
end

%% Section 5: Method Comparison Report (Fair Baselines)
logSection("RP04", "Section 5: Method Comparison (Fair Baselines)", "ChemBERTa BBBP");

logInfo("--- Method comparison (BBBP random 5-fold CV, fair baselines) ---");
logInfo("RP02-rev LR+ECFP4  (sklearn nested CV): AUC = %.4f +/- %.4f", auc_rp02, std_rp02);
logInfo("RP03-rev GCN       (leak-fixed):        AUC = %.4f +/- %.4f", auc_rp03, std_rp03);
logInfo("RP04     ChemBERTa+LR:                  AUC = %.4f +/- %.4f", aucCV, aucStd);
logInfo("");
logInfo("ChemBERTa advantage vs RP02-rev: %+.4f", aucCV - auc_rp02);
logInfo("ChemBERTa advantage vs RP03-rev: %+.4f", aucCV - auc_rp03);

nStdVsRP02 = abs(aucCV - auc_rp02) / aucStd;
nStdVsRP03 = abs(aucCV - auc_rp03) / aucStd;
logInfo("  vs RP02: %.2f std (RP04 std=%.4f)", nStdVsRP02, aucStd);
logInfo("  vs RP03: %.2f std (RP04 std=%.4f)", nStdVsRP03, aucStd);

%% Section 6: Visualization
logSection("RP04", "Section 6: Visualization", "ChemBERTa BBBP");

methodLabels = {"ECFP4+LR (RP02-rev)", "GCN (RP03-rev)", "ChemBERTa+LR (RP04)"};
aucVals = [auc_rp02, auc_rp03, aucCV];
aucSDs  = [std_rp02, std_rp03, aucStd];

barColors = [0.65 0.65 0.65;   % RP02 gray
             0.20 0.50 0.85;   % RP03 blue
             0.85 0.30 0.20];  % RP04 red

fig1 = figure("Name", "RP04: BBBP Method Comparison (Fair Baselines)");
set(fig1, "Position", [100 100 560 420]);
bh = bar(1:3, aucVals, "FaceColor", "flat", "EdgeColor", "none");
bh.CData = barColors;
hold on;
errorbar(1:3, aucVals, aucSDs, "k.", "LineWidth", 1.5, "CapSize", 8);
yline(0.85, "k--", "RF03 >= 0.85", "LineWidth", 1.0, "LabelHorizontalAlignment", "right");
hold off;
ylim([0.80 1.00]);
xticks(1:3);
xticklabels(methodLabels);
ylabel("ROC-AUC (5-fold random CV)");
title("BBBP Classification: Fair Baseline Comparison (M-REPRO-AUDIT B4)");
subtitle("RP02-rev (sklearn nested CV) | RP03-rev (leak-fixed) | RP04 (ChemBERTa)");
grid on; box off;
for k = 1:3
    text(k, aucVals(k) + aucSDs(k) + 0.004, sprintf("%.4f", aucVals(k)), ...
        "HorizontalAlignment", "center", "FontSize", 9, "FontWeight", "bold");
end

fig2 = figure("Name", "RP04: Per-Fold AUC vs Fair Baselines");
set(fig2, "Position", [680 100 500 380]);
foldNums = 1:numel(foldAUC);
bar(foldNums, foldAUC, "FaceColor", [0.85 0.30 0.20], "EdgeColor", "none");
hold on;
yline(aucCV,    "r-",  sprintf("ChemBERTa CV (%.4f)", aucCV), "LineWidth", 1.5);
yline(auc_rp03, "b--", sprintf("RP03-rev GCN (%.4f)", auc_rp03), "LineWidth", 1.2);
yline(auc_rp02, "m--", sprintf("RP02-rev LR (%.4f)", auc_rp02), "LineWidth", 1.2);
hold off;
ylim([0.80 1.00]);
xlabel("Fold"); ylabel("ROC-AUC");
xticks(foldNums);
xticklabels(arrayfun(@(k) sprintf("Fold %d", k), foldNums, "UniformOutput", false));
title("RP04 ChemBERTa: Per-Fold AUC vs Rev Baselines");
grid on; box off;

%% Section 7: RF03 Verification
logSection("RP04", "Section 7: RF03 Verification", "ChemBERTa BBBP");

rf03crit = struct("auc_cv", struct("lower", 0.85));
metRP04r = struct("auc_cv", aucCV);
resRP04r = emk.repro.verify(metRP04r, rf03crit);
logInfo("==> ROC-AUC CV = %.4f (>= 0.85): %s", aucCV, statusStr_(resRP04r.pass));
disp(resRP04r.report);

%% Section 8: Save Results
logSection("RP04", "Section 8: Save Results", "ChemBERTa BBBP");

runDir    = makeRunDir("Prefix", "rp04_chemberta");
absRunDir = char(fullfile(pwd(), runDir));
saveas(fig1, fullfile(absRunDir, "method_comparison_fair.png"));
saveas(fig2, fullfile(absRunDir, "fold_auc_fair.png"));
close(fig1); close(fig2);
logInfo("Figures saved.");

metrics = struct( ...
    "auc_cv",             aucCV, ...
    "auc_cv_std",         aucStd, ...
    "fold_aucs",          foldAUC, ...
    "n_valid",            res.n_valid, ...
    "n_bbb_pos",          res.n_bbb_pos, ...
    "n_bbb_neg",          res.n_bbb_neg, ...
    "model_name",         res.model_name, ...
    "hidden_size",        res.hidden_size, ...
    "n_params_M",         res.n_params_M, ...
    "token_length_stats", tkStats, ...
    "hyperparams",        res.hyperparams, ...
    "comparison_fair", struct( ...
        "rp02_rev_sklearn",  auc_rp02, ...
        "rp02_rev_std",      std_rp02, ...
        "rp03_rev_gcn",      auc_rp03, ...
        "rp03_rev_std",      std_rp03, ...
        "rp04_chemberta",    aucCV, ...
        "chemberta_vs_rp02", aucCV - auc_rp02, ...
        "chemberta_vs_rp03", aucCV - auc_rp03), ...
    "comparison_biased", struct( ...   % historical record only -- do not use for analysis
        "rp02_fitclinear",   0.8826, ...  % biased: fitclinear solver (corrected to 0.9118 in RP02-rev)
        "rp03_gcn_leaked",   0.9151), ... % biased: GCN test leakage (corrected to 0.9038 in RP03-rev)
    "rf03_criteria",      rf03crit, ...
    "rf03_pass",          resRP04r.pass);

fid = fopen(fullfile(runDir, "metrics.json"), "w");
fprintf(fid, "%s\n", jsonencode(metrics, "PrettyPrint", true));
fclose(fid);
logInfo("Metrics saved: metrics.json");

snap.run_date  = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir   = runDir;
snap.rf03_pass = resRP04r.pass;
emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));

logInfo("RP04 complete.  run_dir=%s", runDir);

% ===========================================================================
function s = statusStr_(tf)
    if tf; s = "PASS"; else; s = "FAIL"; end
end
