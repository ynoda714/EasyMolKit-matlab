% rp03_gnn.m  RP03: Graph Neural Network -- GCN on BBBP Classification
%
% Trains a 3-layer Graph Convolutional Network (GCN) on the BBBP dataset
% with 5-fold stratified random CV, benchmarking against the ECFP4+LR
% baseline from RP02 (AUC CV = 0.883).
%
%   Paper:  Yang, K. et al. (2019). Analyzing Learned Molecular
%           Representations for Property Prediction.
%           J. Chem. Inf. Model. 59(8):3370-3388.
%           DOI: 10.1021/acs.jcim.9b00237
%
%   Model:  3-layer GCNConv + BatchNorm + GlobalMeanPool + 2xFC
%             Atom features: 25-dim (type/degree/charge/aromaticity/H)
%             BCEWithLogitsLoss with pos_weight for class imbalance
%             Adam lr=1e-3, patience=20 early stopping
%
%   Task:   Binary classification -- BBB permeability (BBB+ vs BBB-)
%           Dataset: BBBP (MoleculeNet / Wu et al. 2018, 2039 molecules)
%
%   Split:  5-fold stratified random CV (consistent with RP02 / RP05)
%
%   RF01: repro/rp03_gnn/README.md
%   RF02: emk.setup.snapshot() + emk.setup.lockfile()
%   RF03: ROC-AUC CV (5-fold) >= 0.85
%
%   RP02 baseline (ECFP4+LR): AUC CV = 0.8826
%
%   Prerequisites:
%     emk.setup.install() has been run once.
%     torch>=2.0 and torch_geometric>=2.0 in python_env/.
%     (emk.setup.installExtra("torch") and installExtra("torch_geometric"))
%
%   Run: Ctrl+Enter in MATLAB with project root as CWD.
%   NOTE: GCN training takes ~5-10 min on CPU (5 folds x 150 epochs).

%% Section 0: Setup & Environment
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("RP03", "Section 0: Setup & Environment", "GCN BBBP Classification");
emk.setup.initPython();

molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;

snap = emk.setup.snapshot();
logInfo("RP03 setup complete.");

%% Section 1: Load BBBP Dataset & Resolve Paths
logSection("RP03", "Section 1: Load BBBP Dataset", "GCN BBBP Classification");

tbl    = emk.dataset.bbbp();
nTotal = height(tbl);
logInfo("Loaded %d molecules (BBB+: %d, BBB-: %d)", ...
    nTotal, sum(tbl.BBB), sum(~tbl.BBB));

root       = resolveProjectRoot();
csvPath    = fullfile(root, "data", "benchmark", "bbbp.csv");
helperPath = fullfile(thisDir, "rp03_gnn_core.py");

if ~isfile(csvPath)
    error("emk:rp03:csvNotFound", "BBBP CSV not found: %s", csvPath);
end
if ~isfile(helperPath)
    error("emk:rp03:helperNotFound", "Python core not found: %s", helperPath);
end
logInfo("CSV:    %s", csvPath);
logInfo("Helper: %s", helperPath);

%% Section 2: GCN Training -- 5-Fold CV (Python)
logSection("RP03", "Section 2: GCN Training + 5-Fold CV", "GCN BBBP Classification");
logInfo("Loading GCN core and running 5-fold CV ...");
logInfo("(Expected runtime: 5-10 min on CPU with early stopping)");

pyResult = pyrun( ...
    "exec(open(hp).read()); result_json = run_rp03(cp)", ...
    "result_json", hp=helperPath, cp=csvPath);

res = jsondecode(char(string(pyResult)));

aucCV   = res.auc_cv;
aucStd  = res.auc_cv_std;
nValid  = res.n_valid;
foldAUC = res.fold_aucs;   % 5-element double vector

logInfo("GCN 5-fold CV: AUC = %.4f +/- %.4f  (n=%d)", aucCV, aucStd, nValid);
logInfo("Per-fold AUCs: %s", sprintf("%.4f  ", foldAUC));
logInfo("RP02 baseline (ECFP4+LR): AUC = 0.8826");
logInfo("GCN vs baseline: %+.4f", aucCV - 0.8826);

%% Section 3: Extract Learning Curves
logSection("RP03", "Section 3: Extract Learning Curves", "GCN BBBP Classification");

avgLoss = res.avg_train_loss;   % avg train loss per epoch across folds
avgVAUC = res.avg_val_auc;      % avg val AUC per epoch across folds
nEpRun  = res.n_epochs_run;     % actual epochs run (shortest fold)
atomDim = res.atom_feat_dim;
hp      = res.hyperparams;

logInfo("Atom feature dim: %d", atomDim);
logInfo("Epochs (shortest fold): %d / %d (patience=%d)", ...
    nEpRun, hp.n_epochs, hp.patience);

%% Section 4: Visualization
logSection("RP03", "Section 4: Visualization", "GCN BBBP Classification");

epochs = (1:nEpRun)';

% -- Figure 1: Learning curves (avg train loss + avg val AUC) --
fig1 = figure("Name", "RP03 GCN: Learning Curves");
set(fig1, "Position", [100 100 800 380]);

yyaxis left;
plot(epochs, avgLoss, "b-", "LineWidth", 1.5);
ylabel("Mean Train Loss (BCE)");

yyaxis right;
plot(epochs, avgVAUC, "r-", "LineWidth", 1.5);
yline(0.85, "k--", "RF03 >= 0.85", "LineWidth", 1.0, "LabelHorizontalAlignment", "left");
yline(0.8826, "m--", "RP02 LR+ECFP4", "LineWidth", 1.0, "LabelHorizontalAlignment", "left");
ylabel("Mean Val ROC-AUC");
ylim([0.5 1.0]);

xlabel("Epoch");
title(sprintf("RP03 GCN BBBP: Learning Curves (5-fold mean, AUC CV=%.4f)", aucCV));
legend(["Train Loss", "Val AUC", "RF03 >=0.85", "RP02 Baseline"], ...
    "Location", "east");
grid on;

% -- Figure 2: Fold AUC comparison + baseline --
fig2 = figure("Name", "RP03 GCN: Fold AUC vs RP02 Baseline");
set(fig2, "Position", [920 100 500 380]);

barColors = repmat([0.2 0.5 0.85], numel(foldAUC), 1);
bar(1:numel(foldAUC), foldAUC, "FaceColor", "flat", "CData", barColors, "EdgeColor", "none");
hold on;
yline(0.8826, "m--", "RP02 LR+ECFP4 (0.8826)", "LineWidth", 1.5);
yline(aucCV,  "r-",  sprintf("GCN CV mean (%.4f)", aucCV), "LineWidth", 1.5);
hold off;
ylim([0.75 1.0]);
xlabel("Fold");
ylabel("ROC-AUC");
xticks(1:numel(foldAUC));
xticklabels(arrayfun(@(k) sprintf("Fold %d", k), 1:numel(foldAUC), ...
    "UniformOutput", false));
title("RP03 GCN vs RP02 LR+ECFP4: Per-Fold ROC-AUC (BBBP)");
grid on;
box off;

logInfo("Figures 1-2 created.");

%% Section 5: RF03 Verification
logSection("RP03", "Section 5: RF03 Verification", "GCN BBBP Classification");

rf03crit = struct("auc_cv", struct("lower", 0.85));

metRP03 = struct("auc_cv", aucCV);
resRP03 = emk.repro.verify(metRP03, rf03crit);

logInfo("==> ROC-AUC CV = %.4f (>= 0.85): %s", aucCV, statusStr_(resRP03.pass));
disp(resRP03.report);

if resRP03.pass
    logInfo("==> RP03 REPRODUCTION: PASS");
else
    logWarn("==> RP03 REPRODUCTION: NEEDS REVIEW -- see README.md Discussion");
end

%% Section 6: Save Results
logSection("RP03", "Section 6: Save Results", "GCN BBBP Classification");

runDir    = makeRunDir("Prefix", "rp03_gnn");
absRunDir = char(fullfile(pwd(), runDir));

% Save figures
saveas(fig1, fullfile(absRunDir, "learning_curves.png"));
saveas(fig2, fullfile(absRunDir, "fold_auc_comparison.png"));
logInfo("Figures saved.");

% Metrics JSON
metrics = struct( ...
    "auc_cv",          aucCV, ...
    "auc_cv_std",      aucStd, ...
    "fold_aucs",       foldAUC, ...
    "n_valid",         nValid, ...
    "n_bbb_pos",       res.n_bbb_pos, ...
    "n_bbb_neg",       res.n_bbb_neg, ...
    "atom_feat_dim",   atomDim, ...
    "n_epochs_run",    nEpRun, ...
    "rp02_baseline_auc", 0.8826, ...
    "gnn_vs_baseline", aucCV - 0.8826, ...
    "hyperparams",     hp, ...
    "rf03_criteria",   rf03crit, ...
    "rf03_pass",       resRP03.pass);
fid = fopen(fullfile(runDir, "metrics.json"), "w");
fprintf(fid, "%s\n", jsonencode(metrics, "PrettyPrint", true));
fclose(fid);
logInfo("Metrics saved: metrics.json");

% Learning curve CSV
lcTbl = table(epochs, avgLoss(:), avgVAUC(:), ...
    VariableNames=["Epoch", "AvgTrainLoss", "AvgValAUC"]);
writetable(lcTbl, fullfile(runDir, "learning_curves.csv"));
logInfo("Learning curves saved: learning_curves.csv");

% RF02 lock snapshot
snap.run_date  = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir   = runDir;
snap.rf03_pass = resRP03.pass;
emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));

logInfo("RP03 complete.  run_dir=%s", runDir);

% ===========================================================================
% Local helper functions
% ===========================================================================

function s = statusStr_(tf)
    if tf; s = "PASS"; else; s = "FAIL"; end
end
