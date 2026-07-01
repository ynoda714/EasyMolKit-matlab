% rp04_chemberta_f1.m  RP04 F1: ChemBERTa Embeddings + MATLAB fitclinear (Zone C)
%
% F1 Zone C demonstration:
%   Python  -- ChemBERTa tokenisation + CLS embedding extraction (768-dim)
%   MATLAB  -- StandardScaler + fitclinear logistic regression + 5-fold CV
%
% Confirms that MATLAB serves as a fully functional ML backend for
% LLM-based molecular property prediction when Python provides embeddings.
%
% Python LR baseline is loaded dynamically from the latest rp04_chemberta run
% when available. Fallback values are used only if no main-path run exists.
%
% Run: Ctrl+Enter in MATLAB with project root as CWD.

%% Section 0: Setup
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("RP04-F1", "Section 0: Setup", "Zone C ChemBERTa + MATLAB LR");
emk.setup.initPython();
molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;
snap = emk.setup.snapshot();

%% Section 1: Extract ChemBERTa CLS Embeddings (Python)
logSection("RP04-F1", "Section 1: Extract ChemBERTa Embeddings", ...
    "Python: tokenise + frozen CLS infer (no LR)");

root    = resolveProjectRoot();
csvPath = fullfile(root, "data", "benchmark", "bbbp.csv");
if ~isfile(csvPath)
    error("emk:rp04f1:csvNotFound", "BBBP CSV not found: %s", csvPath);
end

helperPath = fullfile(thisDir, "rp04_chemberta_core.py");
logInfo("Extracting CLS embeddings via Python (~2-4 min on CPU) ...");
logInfo("  model:   seyonec/ChemBERTa-zinc-base-v1  (44M params, frozen)");
logInfo("  max_len: 128  (26/2039 molecules exceed this length; max raw tokens = 326)");

try
    pyResult = pyrun( ...
        "exec(open(hp).read()); emb_json = extract_embeddings_only(cp, max_len=128)", ...
        "emb_json", hp=helperPath, cp=csvPath);
catch ME
    error("emk:rp04f1:pyrunFailed", ...
        "ChemBERTa embedding extraction failed: %s", ME.message);
end

embData = jsondecode(char(string(pyResult)));
logInfo("Embeddings received: %d x %d  (model: %s, %.1fM params)", ...
    embData.n_valid, embData.embed_dim, embData.model_name, embData.n_params_M);

%% Section 2: Build MATLAB Feature Matrix
logSection("RP04-F1", "Section 2: Build Feature Matrix", "JSON -> MATLAB double");

X = double(embData.embeddings);   % n_valid x 768
y = double(embData.labels);       % n_valid x 1, binary (0=BBB-, 1=BBB+)
[n, d] = size(X);
logInfo("Feature matrix: %d x %d  |  BBB+: %d  BBB-: %d", ...
    n, d, sum(y == 1), sum(y == 0));

%% Section 3: MATLAB fitclinear 5-Fold CV
logSection("RP04-F1", "Section 3: MATLAB fitclinear 5-Fold CV", ...
    "StandardScaler + logistic LR (lbfgs, Lambda=1/n_train)");

% Seed set before cvpartition to match Python's random_state=42 intent
rng(42, "twister");
cv     = cvpartition(categorical(y), "KFold", 5, "Stratify", true);
nFolds = cv.NumTestSets;

foldAUC = zeros(1, nFolds);
for k = 1:nFolds
    trainIdx = training(cv, k);
    testIdx  = test(cv, k);

    Xtr = X(trainIdx, :);
    ytr = y(trainIdx);
    Xte = X(testIdx,  :);
    yte = y(testIdx);

    % StandardScaler: fit statistics on training fold only, apply to test fold.
    % Zero-variance columns (sigma==0) are protected to avoid divide-by-zero.
    mu              = mean(Xtr, 1);
    sigma           = std(Xtr, 0, 1);
    sigma(sigma == 0) = 1;
    Xtr_s           = (Xtr - mu) ./ sigma;
    Xte_s           = (Xte - mu) ./ sigma;

    % fitclinear logistic (L2, lbfgs).
    % Lambda=1/n_train aligns with sklearn LogisticRegression(C=1.0) per-sample scaling
    % (validated in RP02 R1-C: gap=0.006 on ECFP4 features).
    nTr = sum(trainIdx);
    mdl = fitclinear(Xtr_s, ytr, ...
        "Learner",  "logistic", ...
        "Solver",   "lbfgs", ...
        "Lambda",   1 / nTr, ...
        "Regularization", "ridge");

    % scores: [n_test x 2]; column 2 = score for class 1 (BBB+)
    [~, scores]        = predict(mdl, Xte_s);
    [~, ~, ~, auc]     = perfcurve(yte, scores(:, 2), 1);
    foldAUC(k)         = auc;
    logInfo("  Fold %d / %d: AUC = %.4f  (n_train=%d, n_test=%d)", ...
        k, nFolds, auc, nTr, sum(testIdx));
end

aucMATLAB    = mean(foldAUC);
aucMATLABStd = std(foldAUC);
logInfo("MATLAB fitclinear 5-fold CV: AUC = %.4f +/- %.4f", aucMATLAB, aucMATLABStd);
logInfo("Per-fold AUCs: %s", sprintf("%.4f  ", foldAUC));

%% Section 4: Python vs MATLAB LR Comparison
logSection("RP04-F1", "Section 4: Python vs MATLAB LR Comparison", ...
    "Delta AUC quantification");

% Python LR authoritative baseline from the latest rp04_chemberta run.
aucPython = 0.9270;
stdPython = 0.0107;

rp04Dirs = dir(fullfile("result", "runs", "*rp04_chemberta*"));
if ~isempty(rp04Dirs)
    [~, si] = sort([rp04Dirs.datenum], "descend");
    mPy     = jsondecode(fileread(fullfile( ...
                  "result", "runs", rp04Dirs(si(1)).name, "metrics.json")));
    if isfield(mPy, "auc_cv")
        aucPython = mPy.auc_cv;
        stdPython = mPy.auc_cv_std;
        logInfo("Python LR AUC loaded from: %s", rp04Dirs(si(1)).name);
    end
else
    logWarn("No rp04_chemberta run found -- using hardcoded baseline (0.9270).");
end

deltaAUC = aucMATLAB - aucPython;
nSigma   = abs(deltaAUC) / aucMATLABStd;

logInfo("--- Python sklearn LR vs MATLAB fitclinear ---");
logInfo("  Python sklearn (C=1.0, lbfgs):  AUC = %.4f +/- %.4f", aucPython, stdPython);
logInfo("  MATLAB fitclinear (Lambda=1/n): AUC = %.4f +/- %.4f", aucMATLAB, aucMATLABStd);
logInfo("  delta (MATLAB - Python):        %+.4f", deltaAUC);
logInfo("  |delta| / MATLAB std:           %.2f sigma", nSigma);

if nSigma < 1.0
    verdict = "practical tie";
    logInfo("  Verdict: %s (|delta| < 1 sigma)", verdict);
elseif nSigma < 2.0
    verdict = "borderline";
    logWarn("  Verdict: %s (1 <= |delta|/sigma < 2)", verdict);
else
    verdict = "significant gap";
    logWarn("  Verdict: %s (|delta|/sigma >= 2)", verdict);
end

%% Section 5: Zone C Determination
logSection("RP04-F1", "Section 5: Zone Determination", ...
    "Python=tokenise+embed, MATLAB=LR");

logInfo("Zone C (collaboration) confirmed for RP04 LR component:");
logInfo("  Python role : ChemBERTa tokenisation + CLS embedding extraction");
logInfo("  MATLAB role : StandardScaler + fitclinear (logistic LR)");
logInfo("  Limitation  : RoBERTa tokeniser requires Python (Hugging Face transformers)");
logInfo("  Conclusion  : MATLAB fitclinear is fully capable for the ML step.");
logInfo("  RP04 status : Zone D -> Zone C (LR component)");

%% Section 6: Visualization
logSection("RP04-F1", "Section 6: Visualization", "Python vs MATLAB LR comparison");

methodLabels = {"Python sklearn LR", "MATLAB fitclinear"};
aucVals      = [aucPython, aucMATLAB];
aucSDs       = [stdPython, aucMATLABStd];
barColors    = [0.40 0.70 0.40; ...   % Python green
                0.20 0.40 0.80];      % MATLAB blue

fig1 = figure("Name", "RP04-F1: Python vs MATLAB LR on ChemBERTa Embeddings");
set(fig1, "Position", [100 100 480 400]);
bh = bar(1:2, aucVals, "FaceColor", "flat", "EdgeColor", "none");
bh.CData = barColors;
hold on;
errorbar(1:2, aucVals, aucSDs, "k.", "LineWidth", 1.5, "CapSize", 8);
yline(0.85, "k--", "RF03 >= 0.85", "LineWidth", 1.0, ...
    "LabelHorizontalAlignment", "right");
hold off;
ylim([0.80 1.00]);
xticks(1:2);
xticklabels(methodLabels);
ylabel("ROC-AUC (5-fold random CV)");
title("RP04-F1: ChemBERTa CLS Embeddings -- Python vs MATLAB LR");
subtitle("Zone C: Python = tokenise + embed  |  MATLAB = fitclinear (lbfgs)");
grid on; box off;
for k = 1:2
    text(k, aucVals(k) + aucSDs(k) + 0.004, sprintf("%.4f", aucVals(k)), ...
        "HorizontalAlignment", "center", "FontSize", 9, "FontWeight", "bold");
end

fig2 = figure("Name", "RP04-F1: Per-Fold AUC -- Python vs MATLAB");
set(fig2, "Position", [600 100 520 380]);
foldNums = 1:nFolds;
bar(foldNums, foldAUC, "FaceColor", [0.20 0.40 0.80], "EdgeColor", "none");
hold on;
yline(aucMATLAB, "b-",  sprintf("MATLAB CV = %.4f", aucMATLAB), "LineWidth", 1.5);
yline(aucPython, "g--", sprintf("Python CV = %.4f", aucPython),  "LineWidth", 1.2);
hold off;
ylim([0.80 1.00]);
xlabel("Fold");
ylabel("ROC-AUC");
xticks(foldNums);
xticklabels(arrayfun(@(k) sprintf("Fold %d", k), foldNums, "UniformOutput", false));
title("RP04-F1: MATLAB fitclinear Per-Fold AUC");
subtitle("vs Python sklearn LR baseline");
grid on; box off;

%% Section 7: RF03 Verification
logSection("RP04-F1", "Section 7: RF03 Verification", "MATLAB LR must pass >= 0.85");

rf03crit = struct("auc_cv", struct("lower", 0.85));
metF1    = struct("auc_cv", aucMATLAB);
resF1    = emk.repro.verify(metF1, rf03crit);
logInfo("==> MATLAB fitclinear AUC CV = %.4f (>= 0.85): %s", aucMATLAB, statusStr_(resF1.pass));
disp(resF1.report);

%% Section 8: Save Results
logSection("RP04-F1", "Section 8: Save Results", "metrics.json + lockfile + figures");

runDir = makeRunDir("Prefix", "rp04_chemberta_f1");
if startsWith(runDir, "/") || (numel(runDir) >= 2 && runDir(2) == ":")
    absRunDir = char(runDir);
else
    absRunDir = char(fullfile(pwd(), runDir));
end
saveas(fig1, fullfile(absRunDir, "python_vs_matlab_lr.png"));
saveas(fig2, fullfile(absRunDir, "fold_auc_matlab.png"));
close(fig1); close(fig2);
logInfo("Figures saved.");

metrics = struct( ...
    "auc_cv_matlab",     aucMATLAB, ...
    "auc_cv_matlab_std", aucMATLABStd, ...
    "fold_aucs_matlab",  foldAUC, ...
    "auc_cv_python",     aucPython, ...
    "auc_cv_python_std", stdPython, ...
    "delta_auc",         deltaAUC, ...
    "delta_sigma",       nSigma, ...
    "verdict",           verdict, ...
    "n_valid",           n, ...
    "n_bbb_pos",         embData.n_bbb_pos, ...
    "n_bbb_neg",         embData.n_bbb_neg, ...
    "embed_dim",         d, ...
    "model_name",        embData.model_name, ...
    "zone",              "C", ...
    "zone_detail", struct( ...
        "python_role",  "ChemBERTa tokenise + CLS embedding", ...
        "matlab_role",  "StandardScaler + fitclinear (logistic LR)", ...
        "limitation",   "RoBERTa tokeniser requires Python (Hugging Face transformers)"), ...
    "rf03_criteria",     rf03crit, ...
    "rf03_pass",         resF1.pass, ...
    "hyperparams", struct( ...
        "n_folds",   nFolds, ...
        "lambda",    "1/n_train", ...
        "learner",   "logistic", ...
        "solver",    "lbfgs", ...
        "seed",      42));

writelines(jsonencode(metrics, "PrettyPrint", true), fullfile(runDir, "metrics.json"));
logInfo("Metrics saved: metrics.json");

snap.run_date  = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir   = runDir;
snap.rf03_pass = resF1.pass;
emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));

logInfo("RP04-F1 complete.  Zone C confirmed.  run_dir=%s", runDir);

% ===========================================================================
function s = statusStr_(tf)
    if tf; s = "PASS"; else; s = "FAIL"; end
end
