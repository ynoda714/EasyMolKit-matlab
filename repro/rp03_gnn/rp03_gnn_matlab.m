% rp03_gnn_matlab.m  RP03: MATLAB DLT GCN BBBP Classification
%
% Implements a 3-layer GCN using MATLAB Deep Learning Toolbox (DLT) and
% compares it with the latest available Python/PyG reference artifact.
%
% Outer fold indices are shared with rp02_bbbp (A2 alignment) so that
% per-fold AUC values are comparable across RP02 / RP03 Python / RP03 MATLAB.
%
%   Python reference:  loaded from the newest result/runs/*rp03_gnn*/ artifact
%   RF03 criterion:    AUC CV >= 0.85
%
%   Output: result/runs/rp03_gnn_matlab_<ts>/
%     metrics_matlab.json, fold_predictions.csv,
%     fold_auc_matlab.png, featurized.json, lock_snapshot.json
%
% Run: Ctrl+Enter each section in MATLAB with project root as CWD.
% Mini-batch (bs=64) with padded 3D tensors and pagemtimes; val AUC evaluated
% every eval_freq epochs. Expected runtime: ~9 min on GPU, 30-60 min on CPU.

%% Section 0a: Parameters
hidden      = 128;
lr          = 1e-3;   % same as Python; mini-batch (bs=64) gives ~20 Adam steps/epoch
batch_size  = 64;     % mini-batch size (Run-5: added to match Python ~3000 steps total)
use_bn      = false;
n_epochs    = 150;    % matches Python n_epochs=150
patience    = 20;
seed        = 42;
useGPU      = canUseGPU();
weight_decay = 1e-5;
eval_freq   = 5;      % val AUC evaluated every N epochs
n_fold      = 5;

python_auc_cv  = NaN;
python_auc_std = NaN;
python_ref_src = "";

%% Section 0b: Environment setup
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));
addpath(thisDir);           % expose rp03_gnn_matlab_core.m as local function file

logSection("RP03-MATLAB", "Section 0b: Setup", "MATLAB DLT GCN");
emk.setup.initPython();
molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;
snap = emk.setup.snapshot();

runDir = makeRunDir("Prefix", "rp03_gnn_matlab");
runDirStr = char(runDir);
if startsWith(runDirStr, '/') || (numel(runDirStr) >= 2 && runDirStr(2) == ':')
    absRunDir = runDirStr;
else
    absRunDir = char(fullfile(pwd(), runDir));
end
logInfo("Run directory: %s", absRunDir);
logInfo("useGPU = %d", useGPU);

tbl     = emk.dataset.bbbp();
csvPath = fullfile(pwd(), "data", "benchmark", "bbbp.csv");
logInfo("BBBP: %d total molecules in table", height(tbl));

[python_auc_cv, python_auc_std, python_ref_src] = loadLatestPythonRef_(pwd());
if ~isnan(python_auc_cv)
    logInfo("Python reference loaded: AUC CV = %.4f +/- %.4f  (%s)", ...
        python_auc_cv, python_auc_std, python_ref_src);
else
    logWarn("No RP03 Python reference artifact found under result/runs. Comparison section will be partial.");
end

%% Section 1: Load fold indices (A2 alignment with rp02_bbbp)
logSection("RP03-MATLAB", "Section 1: Fold Indices", "MATLAB DLT GCN");

coreScript = fullfile(thisDir, "rp03_gnn_core.py");
if ~isfile(coreScript)
    error("emk:repro:rp03:missingScript", "Core script not found: %s", coreScript);
end

rp02Runs = dir(fullfile("result", "runs", "*rp02_bbbp*"));
if isempty(rp02Runs)
    logWarn("No rp02_bbbp run found. Inner KFold will be used (fold alignment not applied).");
    foldIdxPath = "";
    foldIdx     = struct();
    hasFoldIdx  = false;
else
    [~, sortIdx]  = sort([rp02Runs.datenum], "descend");
    latestRp02    = fullfile("result", "runs", rp02Runs(sortIdx(1)).name);
    foldIdxPath   = fullfile(pwd(), latestRp02, "outer_fold_indices.json");
    if isfile(foldIdxPath)
        foldIdx    = jsondecode(fileread(foldIdxPath));
        hasFoldIdx = true;
        logInfo("A2 fold indices: %s", foldIdxPath);
    else
        logWarn("outer_fold_indices.json not found in %s. Inner KFold used.", latestRp02);
        foldIdx    = struct();
        hasFoldIdx = false;
    end
end

if ~hasFoldIdx
    error("emk:repro:rp03:noFoldIdx", ...
        "A2 fold indices required for RP03-MATLAB. Run rp02_bbbp.m first.");
end

%% Section 2: Python featurization -> molecules struct array
logSection("RP03-MATLAB", "Section 2: Python Featurization", ...
    "MATLAB DLT GCN");
logInfo("Calling featurize_bbbp() via pyrun -> featurized.json ...");

tmpJson = fullfile(absRunDir, "featurized.json");
% isfile guards against missing file; coreScript is hardcoded to a known project-local path.
pyrun(["exec(open(core_script).read()); " ...
       "import json, pathlib; " ...
       "pathlib.Path(tmp).write_text(json.dumps(featurize_bbbp(csv_path)), 'utf-8')"], ...
      core_script=coreScript, tmp=tmpJson, csv_path=csvPath);

data      = jsondecode(fileread(tmpJson));
molecules = data.molecules;
n_valid   = data.n_valid;
n_bbb_pos = data.n_bbb_pos;
n_bbb_neg = data.n_bbb_neg;

logInfo("Featurized: n_valid=%d  n_bbb_pos=%d  n_bbb_neg=%d", n_valid, n_bbb_pos, n_bbb_neg);
logInfo("featurized.json saved to %s", tmpJson);

%% Section 3: MATLAB GCN training (5-fold CV)
logSection("RP03-MATLAB", "Section 3: MATLAB GCN Training", ...
    "MATLAB DLT GCN");
logInfo("Architecture: GCNConv(25->128)x3 + GlobalMeanPool + FC(64) + FC(1)");
logInfo("Hyperparams: hidden=%d lr=%.4f n_epochs=%d patience=%d seed=%d useGPU=%d", ...
    hidden, lr, n_epochs, patience, seed, useGPU);
logInfo("Batch mode: 3D padded tensors + pagemtimes, bs=%d, eval_freq=%d", batch_size, eval_freq);
logInfo("Expected runtime: ~9 min on GPU (5-fold).");

t_train = tic;
mResult = rp03_gnn_matlab_core(molecules, foldIdx, ...
    'Hidden', hidden, 'LR', lr, 'NEpochs', n_epochs, ...
    'Patience', patience, 'Seed', seed, 'UseGPU', useGPU, ...
    'WeightDecay', weight_decay, 'EvalFreq', eval_freq, ...
    'BatchSize', batch_size, 'UseBN', use_bn, ...
    'SaveWeightDir', absRunDir);
elapsed_matlab = toc(t_train);

fold_aucs_matlab = mResult.fold_aucs;
auc_cv_matlab    = mResult.auc_cv;
auc_std_matlab   = mResult.auc_cv_std;

logInfo("");
logInfo("=== MATLAB GCN 5-fold CV Results ===");
for k = 1:n_fold
    logInfo("  Fold %d: AUC = %.4f", k, fold_aucs_matlab(k));
end
logInfo("  Mean AUC = %.4f +/- %.4f", auc_cv_matlab, auc_std_matlab);
logInfo("  Elapsed:  %.1f sec (%.1f min)", elapsed_matlab, elapsed_matlab / 60);

%% Section 4: Comparison report (Python vs MATLAB)
logSection("RP03-MATLAB", "Section 4: Comparison Report", ...
    "MATLAB DLT GCN");

delta_auc = auc_cv_matlab - python_auc_cv;

logInfo("=== Python GCN vs MATLAB DLT GCN ===");
if ~isnan(python_auc_cv)
    logInfo("  Python AUC CV : %.4f +/- %.4f  (PyTorch-Geometric, batch_size=64)", ...
        python_auc_cv, python_auc_std);
    logInfo("    source      : %s", python_ref_src);
else
    logWarn("  Python AUC CV : unavailable");
end
logInfo("  MATLAB AUC CV : %.4f +/- %.4f  (DLT, batch_size=%d, BN=%s)", ...
    auc_cv_matlab, auc_std_matlab, batch_size, mat2str(use_bn));
if ~isnan(delta_auc)
    logInfo("  Delta AUC     : %+.4f", delta_auc);
else
    logWarn("  Delta AUC     : unavailable");
end
logInfo("  MATLAB elapsed: %.1f sec  |  Python elapsed: see %s", elapsed_matlab, python_ref_src);
logInfo("");
if ~isnan(delta_auc) && abs(delta_auc) < max(python_auc_std, auc_std_matlab)
    logInfo("  >> AUC gap < 1 std: implementations are practically equivalent.");
elseif ~isnan(delta_auc)
    logInfo("  >> AUC gap >= 1 std: see README implementation notes for root causes.");
end

%% Section 5: Visualization
logSection("RP03-MATLAB", "Section 5: Visualization", ...
    "MATLAB DLT GCN");

if exist("figM", "var") && isvalid(figM); close(figM); end

figM = figure("Name", "RP03 MATLAB GCN: Fold AUC");
set(figM, "Position", [100 100 560 400]);
fold_ids = (1:n_fold)';
bar(fold_ids, fold_aucs_matlab, "FaceColor", [0.85 0.4 0.1], "EdgeColor", "none");
hold on;
if ~isnan(python_auc_cv)
    yline(python_auc_cv, "b--", sprintf("Python GCN (%.4f)", python_auc_cv), ...
        "LineWidth", 1.5, "LabelHorizontalAlignment", "left");
end
yline(auc_cv_matlab, "r-",  sprintf("MATLAB GCN (%.4f)", auc_cv_matlab), ...
    "LineWidth", 1.5, "LabelHorizontalAlignment", "left");
hold off;
allVals    = [fold_aucs_matlab; auc_cv_matlab];
if ~isnan(python_auc_cv)
    allVals = [allVals; python_auc_cv];
end
figMargin  = max(0.05, (max(allVals) - min(allVals)) * 0.15);
ylim([max(0.0, min(allVals) - figMargin), min(1.0, max(allVals) + figMargin)]);
xlabel("Fold"); ylabel("Test ROC-AUC");
xticks(fold_ids);
xticklabels(arrayfun(@(k) sprintf("Fold %d", k), fold_ids, "UniformOutput", false));
title(sprintf("RP03 MATLAB DLT GCN vs Python GCN: Per-Fold AUC (BBBP)"), ...
    sprintf("MATLAB=%.4f  Python=%s  delta=%s", ...
    auc_cv_matlab, localNumStr_(python_auc_cv), localNumStr_(delta_auc)));
grid on; box off;
logInfo("Figure ready (figM). Run Section 6 to save all results.");

%% Section 6: Save results
logSection("RP03-MATLAB", "Section 6: Save Results", ...
    "MATLAB DLT GCN");

if ~exist("mResult", "var")
    error("emk:repro:rp03:noResult", "Run Section 3 first: mResult not defined.");
end

% Figure
if exist("figM", "var") && isvalid(figM)
    saveas(figM, fullfile(absRunDir, "fold_auc_matlab.png"));
    close(figM);
    logInfo("Figure saved.");
end

% metrics_matlab.json
metrics_matlab = struct( ...
    "model",            "MATLAB_GCN", ...
    "auc_cv",           auc_cv_matlab, ...
    "auc_cv_std",       auc_std_matlab, ...
    "fold_aucs",        fold_aucs_matlab, ...
    "elapsed_sec_total", elapsed_matlab, ...
    "n_valid",          n_valid, ...
    "fold_source",      "a2_shared_rp02rev", ...
    "hyperparams", struct( ...
        "hidden",      hidden, ...
        "n_layers",    3, ...
        "dropout",     0.3, ...
        "lr",          lr, ...
        "n_epochs",    n_epochs, ...
        "patience",    patience, ...
        "batch_size",  batch_size, ...
        "seed",        seed, ...
        "use_gpu",     useGPU), ...
    "implementation_notes", struct( ...
        "batch_mode",          "padded_3d_tensors_pagemtimes", ...
        "batchnorm",           ternaryStr_(use_bn, "f2_learnable_gamma_beta_ema_running_stats", "off"), ...
        "dropout",             "inverted_random_mask", ...
        "kaiming_init",        true, ...
        "a_norm_computed",     "precomputed_3d_tensor_before_fold_loop", ...
        "val_eval_freq",       eval_freq), ...
    "comparison", struct( ...
        "python_auc_cv",      python_auc_cv, ...
        "python_auc_std",     python_auc_std, ...
        "python_reference_source", python_ref_src, ...
        "matlab_auc_cv",      auc_cv_matlab, ...
        "delta_auc",          delta_auc, ...
        "python_elapsed_sec_est", [], ...
        "matlab_elapsed_sec", elapsed_matlab));

fid = fopen(fullfile(absRunDir, "metrics_matlab.json"), "w", "n", "UTF-8");
fprintf(fid, "%s\n", jsonencode(metrics_matlab, "PrettyPrint", true));
fclose(fid);
logInfo("metrics_matlab.json saved.");

% fold_predictions.csv
allPreds = vertcat(mResult.all_preds{:});
predTbl  = array2table(allPreds, "VariableNames", ["fold", "mol_idx", "score", "label"]);
writetable(predTbl, fullfile(absRunDir, "fold_predictions.csv"));
logInfo("fold_predictions.csv saved (%d rows).", height(predTbl));

% lock_snapshot.json
snap.run_date   = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir    = char(runDir);
snap.matlab_gcn_auc_cv = auc_cv_matlab;
snap.rf03_pass  = auc_cv_matlab >= 0.85;
emk.setup.lockfile(snap, fullfile(absRunDir, "lock_snapshot.json"));

logInfo("Results saved: %s", absRunDir);
if auc_cv_matlab >= 0.85
    logInfo("==> RP03-MATLAB: AUC CV = %.4f (>= 0.85) PASS", auc_cv_matlab);
else
    logWarn("==> RP03-MATLAB: AUC CV = %.4f (< 0.85) NEEDS REVIEW", auc_cv_matlab);
end

function out = ternaryStr_(tf, whenTrue, whenFalse)
if tf
    out = whenTrue;
else
    out = whenFalse;
end
end

function [auc_cv, auc_std, sourcePath] = loadLatestPythonRef_(projectRoot)
auc_cv = NaN;
auc_std = NaN;
sourcePath = "";

runDirs = dir(fullfile(projectRoot, "result", "runs", "*rp03_gnn*"));
runDirs = runDirs([runDirs.isdir]);
if isempty(runDirs)
    return;
end

[~, order] = sort([runDirs.datenum], "descend");
runDirs = runDirs(order);

for k = 1:numel(runDirs)
    if contains(runDirs(k).name, "matlab", "IgnoreCase", true)
        continue;
    end
    baseDir = fullfile(runDirs(k).folder, runDirs(k).name);
    candidates = {fullfile(baseDir, "metrics.json"), fullfile(baseDir, "raw_result.json")};
    for j = 1:numel(candidates)
        if ~isfile(candidates{j})
            continue;
        end
        rawText = fileread(candidates{j});
        if ~isempty(rawText) && rawText(1) == char(65279)
            rawText = rawText(2:end);
        end
        try
            data = jsondecode(rawText);
        catch
            continue;
        end
        if isfield(data, "auc_cv") && isfield(data, "auc_cv_std")
            auc_cv = double(data.auc_cv);
            auc_std = double(data.auc_cv_std);
            sourcePath = candidates{j};
            return;
        end
    end
end
end

function s = localNumStr_(x)
if isnan(x)
    s = "N/A";
else
    s = sprintf("%.4f", x);
end
end
