% rp04_chemberta.m  RP04: ChemBERTa Molecular Language Model on BBBP
%
% Extracts CLS token embeddings from frozen pre-trained ChemBERTa
% (seyonec/ChemBERTa-zinc-base-v1, 44M params) and trains a logistic
% regression classifier for BBBP BBB permeability classification.
%
%   Paper:  Chithrananda, S. et al. (2020). ChemBERTa: Large-Scale
%           Self-Supervised Pretraining for Molecular Property Prediction.
%           arXiv:2010.09885.
%
%   Model:  seyonec/ChemBERTa-zinc-base-v1 (HuggingFace)
%             RoBERTa-based, 44M params, pre-trained on ZINC SMILES.
%             Approach: frozen CLS embedding (768-dim) + LR (linear probe).
%
%   Task:   Binary classification -- BBB permeability (BBB+ vs BBB-)
%           Dataset: BBBP (MoleculeNet / Wu et al. 2018, 2039 molecules)
%
%   Split:  5-fold stratified random CV (consistent with RP02 / RP03)
%
%   Comparison (same dataset + split):
%     RP02: ECFP4 + LR          AUC CV = 0.8826
%     RP03: GCN (3-layer)       AUC CV = 0.9151
%     RP04: ChemBERTa CLS + LR  AUC CV = ?
%
%   RF01: repro/rp04_chemberta/README.md
%   RF02: emk.setup.snapshot() + emk.setup.lockfile()
%   RF03: ROC-AUC CV (5-fold) >= 0.85
%
%   Prerequisites:
%     emk.setup.install() has been run once.
%     transformers>=5.0 in python_env/ (emk.setup.installExtra("transformers")).
%     Model cached: seyonec/ChemBERTa-zinc-base-v1 (~280MB, auto-downloaded if absent).
%
%   Run: Ctrl+Enter in MATLAB with project root as CWD.
%   NOTE: Embedding extraction takes ~2-4 min on CPU (no training, just forward passes).

%% Section 0: Setup & Environment
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("RP04", "Section 0: Setup & Environment", "ChemBERTa BBBP");
emk.setup.initPython();

molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;

snap = emk.setup.snapshot();
logInfo("RP04 setup complete.");

%% Section 1: Load BBBP Dataset & Resolve Paths
logSection("RP04", "Section 1: Load BBBP Dataset", "ChemBERTa BBBP");

tbl    = emk.dataset.bbbp();
nTotal = height(tbl);
logInfo("Loaded %d molecules (BBB+: %d, BBB-: %d)", ...
    nTotal, sum(tbl.BBB), sum(~tbl.BBB));

root       = resolveProjectRoot();
csvPath    = fullfile(root, "data", "benchmark", "bbbp.csv");
helperPath = fullfile(thisDir, "rp04_chemberta_core.py");

if ~isfile(csvPath)
    error("emk:rp04:csvNotFound", "BBBP CSV not found: %s", csvPath);
end
if ~isfile(helperPath)
    error("emk:rp04:helperNotFound", "Python core not found: %s", helperPath);
end
logInfo("CSV:    %s", csvPath);
logInfo("Helper: %s", helperPath);

%% Section 2: ChemBERTa Embedding Extraction + 5-Fold CV (Python)
logSection("RP04", "Section 2: ChemBERTa Embedding + 5-Fold CV", "ChemBERTa BBBP");
logInfo("Loading ChemBERTa (seyonec/ChemBERTa-zinc-base-v1, 44M params) ...");
logInfo("Model is cached locally -- no download needed.");
logInfo("Extracting CLS embeddings for %d molecules (~2-4 min on CPU)...", nTotal);

pyResult = pyrun( ...
    "exec(open(hp).read()); result_json = run_rp04(cp)", ...
    "result_json", hp=helperPath, cp=csvPath);

res = jsondecode(char(string(pyResult)));

aucCV      = res.auc_cv;
aucStd     = res.auc_cv_std;
nValid     = res.n_valid;
foldAUC    = res.fold_aucs;
hiddenSize = res.hidden_size;
nParamsM   = res.n_params_M;

logInfo("ChemBERTa CLS dim: %d  |  Model params: %.1fM", hiddenSize, nParamsM);
logInfo("5-fold CV: AUC = %.4f +/- %.4f  (n=%d)", aucCV, aucStd, nValid);
logInfo("Per-fold AUCs: %s", sprintf("%.4f  ", foldAUC));
logInfo("--- Method comparison (BBBP random 5-fold CV) ---");
logInfo("RP02 ECFP4+LR:        AUC = 0.8826");
logInfo("RP03 GCN (3-layer):   AUC = 0.9151");
logInfo("RP04 ChemBERTa+LR:    AUC = %.4f  (delta vs RP02: %+.4f, vs RP03: %+.4f)", ...
    aucCV, aucCV - 0.8826, aucCV - 0.9151);

%% Section 3: Comparison Table
logSection("RP04", "Section 3: Method Comparison", "ChemBERTa BBBP");

methods = {"ECFP4+LR (RP02)", "GCN (RP03)", "ChemBERTa+LR (RP04)"};
aucVals = [0.8826, 0.9151, aucCV];
aucSDs  = [0.022,  0.019,  aucStd];   % RP02 std=0.022, RP03 std=0.019

logInfo("%-25s  AUC CV    Std", "Method");
for k = 1:3
    logInfo("%-25s  %.4f    %.4f", methods{k}, aucVals(k), aucSDs(k));
end

%% Section 4: Visualization
logSection("RP04", "Section 4: Visualization", "ChemBERTa BBBP");

% -- Figure 1: Three-way comparison bar chart --
fig1 = figure("Name", "RP04: BBBP Method Comparison");
set(fig1, "Position", [100 100 560 420]);

barColors = [0.65 0.65 0.65;   % RP02 gray
             0.20 0.50 0.85;   % RP03 blue
             0.85 0.30 0.20];  % RP04 red
xPos = 1:3;
bh   = bar(xPos, aucVals, "FaceColor", "flat", "EdgeColor", "none");
bh.CData = barColors;
hold on;
er = errorbar(xPos, aucVals, aucSDs, "k.", "LineWidth", 1.5, "CapSize", 8);
yline(0.85, "k--", "RF03 >= 0.85", "LineWidth", 1.0, "LabelHorizontalAlignment", "right");
hold off;

ylim([0.80 1.00]);
xticks(xPos);
xticklabels(methods);
ylabel("ROC-AUC (5-fold random CV)");
title("BBBP Classification: ECFP4+LR vs GCN vs ChemBERTa");
subtitle("Same dataset (BBBP 2039 mol), same split (random 5-fold stratified)");
grid on; box off;

% Annotate bars with AUC values
for k = 1:3
    text(xPos(k), aucVals(k) + aucSDs(k) + 0.004, sprintf("%.4f", aucVals(k)), ...
        "HorizontalAlignment", "center", "FontSize", 9, "FontWeight", "bold");
end

% -- Figure 2: Per-fold AUC (ChemBERTa) --
fig2 = figure("Name", "RP04: ChemBERTa Per-Fold AUC");
set(fig2, "Position", [680 100 500 380]);

foldNums = 1:numel(foldAUC);
bar(foldNums, foldAUC, "FaceColor", [0.85 0.30 0.20], "EdgeColor", "none");
hold on;
yline(aucCV,   "r-",  sprintf("ChemBERTa CV mean (%.4f)", aucCV), "LineWidth", 1.5);
yline(0.9151,  "b--", "RP03 GCN (0.9151)", "LineWidth", 1.0);
yline(0.8826,  "m--", "RP02 LR+ECFP4 (0.8826)", "LineWidth", 1.0);
hold off;
ylim([0.75 1.00]);
xlabel("Fold");
ylabel("ROC-AUC");
xticks(foldNums);
xticklabels(arrayfun(@(k) sprintf("Fold %d", k), foldNums, "UniformOutput", false));
title("RP04 ChemBERTa: Per-Fold ROC-AUC vs RP02/RP03");
grid on; box off;

logInfo("Figures 1-2 created.");

%% Section 5: RF03 Verification
logSection("RP04", "Section 5: RF03 Verification", "ChemBERTa BBBP");

rf03crit = struct("auc_cv", struct("lower", 0.85));
metRP04  = struct("auc_cv", aucCV);
resRP04  = emk.repro.verify(metRP04, rf03crit);

logInfo("==> ROC-AUC CV = %.4f (>= 0.85): %s", aucCV, statusStr_(resRP04.pass));
disp(resRP04.report);

if resRP04.pass
    logInfo("==> RP04 REPRODUCTION: PASS");
else
    logWarn("==> RP04 REPRODUCTION: NEEDS REVIEW -- see README.md Discussion");
end

%% Section 6: Save Results
logSection("RP04", "Section 6: Save Results", "ChemBERTa BBBP");

runDir    = makeRunDir("Prefix", "rp04_chemberta");
absRunDir = char(fullfile(pwd(), runDir));

% Save figures
saveas(fig1, fullfile(absRunDir, "method_comparison.png"));
saveas(fig2, fullfile(absRunDir, "fold_auc_chemberta.png"));
logInfo("Figures saved.");

% Metrics JSON
metrics = struct( ...
    "auc_cv",             aucCV, ...
    "auc_cv_std",         aucStd, ...
    "fold_aucs",          foldAUC, ...
    "n_valid",            nValid, ...
    "n_bbb_pos",          res.n_bbb_pos, ...
    "n_bbb_neg",          res.n_bbb_neg, ...
    "model_name",         res.model_name, ...
    "hidden_size",        hiddenSize, ...
    "n_params_M",         nParamsM, ...
    "hyperparams",        res.hyperparams, ...
    "comparison",         struct( ...
        "rp02_ecfp4_lr",  0.8826, ...
        "rp03_gcn",       0.9151, ...
        "rp04_chemberta", aucCV), ...
    "rf03_criteria",      rf03crit, ...
    "rf03_pass",          resRP04.pass);
fid = fopen(fullfile(runDir, "metrics.json"), "w");
fprintf(fid, "%s\n", jsonencode(metrics, "PrettyPrint", true));
fclose(fid);
logInfo("Metrics saved: metrics.json");

% RF02 lock snapshot
snap.run_date  = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir   = runDir;
snap.rf03_pass = resRP04.pass;
emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));

logInfo("RP04 complete.  run_dir=%s", runDir);

% ===========================================================================
% Local helper functions
% ===========================================================================

function s = statusStr_(tf)
    if tf; s = "PASS"; else; s = "FAIL"; end
end
