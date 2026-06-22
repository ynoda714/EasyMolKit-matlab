% rp05_shap.m  RP05: Explainable AI -- SHAP for BBBP ECFP4 + Logistic Regression
%
% Applies shap.LinearExplainer to interpret a Logistic Regression model
% trained on BBBP Blood-Brain Barrier data with Morgan ECFP4 fingerprints,
% following the methodology of Rodriguez-Perez & Bajorath (2020).
%
%   Paper:  Rodriguez-Perez, R. & Bajorath, J. (2020). Interpretation of
%           Machine Learning Models Using Shapley Values: Application to
%           Compound Potency and Multi-target Activity Predictions.
%           J. Comput.-Aided Mol. Des. 34:1013-1026.
%           DOI: 10.1007/s10822-020-00314-0
%
%   Model:  sklearn LogisticRegression (L2 ridge, C=1.0) + Morgan ECFP4
%             Radius=2, NBits=2048 -- consistent with RP02.
%
%   SHAP:   shap.LinearExplainer -- exact for linear models.
%             global_imp[j] = mean_i(|coef[j] * (X[i,j] - mean(X[:,j]))|)
%
%   Task:   Feature attribution -- which ECFP4 bits drive BBB permeability?
%
%   Split:  5-fold stratified random CV (consistent with RP02).
%
%   RF01: repro/rp05_shap/README.md
%   RF02: emk.setup.snapshot() + emk.setup.lockfile()
%   RF03: (1) ROC-AUC CV (5-fold) >= 0.85
%         (2) Spearman(global_imp, |coef|*std(X)) >= 0.90
%             (verifies LinearExplainer is consistent with LR weights)
%
%   Prerequisites:
%     emk.setup.install() has been run once.
%     shap>=0.49.1 and scikit-learn>=1.7 installed in python_env/.
%     (Install: python_env/python.exe -m pip install shap scikit-learn)
%
%   Run: Ctrl+Enter in MATLAB with project root as CWD.

%% Section 0: Setup & Environment
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("RP05", "Section 0: Setup & Environment", "SHAP BBBP Explainability");
emk.setup.initPython();

molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;

snap = emk.setup.snapshot();
logInfo("RP05 setup complete.");

%% Section 1: Resolve Paths & Ensure Dataset
logSection("RP05", "Section 1: Resolve Paths & Ensure Dataset", "SHAP BBBP Explainability");

root       = resolveProjectRoot();
csvPath    = fullfile(root, "data", "benchmark", "bbbp.csv");
helperPath = fullfile(thisDir, "rp05_shap_core.py");

% Trigger download/cache if not present; counts are reported from Python
% after SMILES validation (single source of truth, avoids hash mismatch).
if ~isfile(csvPath)
    emk.dataset.bbbp();
end

if ~isfile(csvPath)
    error("emk:rp05:csvNotFound", "BBBP CSV not found: %s", csvPath);
end
if ~isfile(helperPath)
    error("emk:rp05:helperNotFound", "Python core not found: %s", helperPath);
end
logInfo("CSV:    %s", csvPath);
logInfo("Helper: %s", helperPath);

%% Section 2: Fit sklearn LR + 5-Fold CV AUC (Python)
logSection("RP05", "Section 2: Fit sklearn LR + Cross-Validation", ...
    "SHAP BBBP Explainability");

logInfo("Running Python: load data, ECFP4, fit LR, 5-fold CV, SHAP ...");
logInfo("(first run may take ~30 s for numba JIT compilation in shap)");

pyResult = pyrun( ...
    "exec(open(hp).read()); result_json = run_rp05(cp)", ...
    "result_json", hp=helperPath, cp=csvPath);

res = jsondecode(char(string(pyResult)));

aucCV   = res.auc_cv;
aucStd  = res.auc_cv_std;
nValid  = res.n_valid;
nTrain  = res.n_train;
nTest   = res.n_test;
shapRho = res.shap_lr_spearman;

logInfo("Dataset: %d valid (%d train / %d test), BBB+: %d, BBB-: %d", ...
    nValid, nTrain, nTest, res.n_bbb_pos, nValid - res.n_bbb_pos);
logInfo("sklearn LR 5-fold CV: AUC=%.4f +/- %.4f  (full dataset n=%d)", ...
    aucCV, aucStd, nValid);
logInfo("SHAP LinearExplainer: Spearman(global_imp, |coef|*std_train) = %.4f", shapRho);

%% Section 3: Extract SHAP Results for MATLAB Visualization
logSection("RP05", "Section 3: Extract SHAP Results", "SHAP BBBP Explainability");

% Global importance (top 20 bits) -- 0-based Python indices -> 1-based MATLAB
topN    = numel(res.top_n_imp);
topIdx  = res.top_n_idx + 1;           % 1-based MATLAB bit indices
topImp  = res.top_n_imp;               % mean |SHAP|
topCoef = res.top_n_coef;              % LR coefficient

% 3 example molecules (TP / TN / misclassified)
exTypes  = res.ex_types;               % cell of strings: {"TP","TN","MIS"}
exLabels = res.ex_labels;              % true labels (0/1)
exPreds  = res.ex_preds;               % predicted labels (0/1)
exProbs  = res.ex_probs;               % P(BBB+)
exNames  = res.ex_names;
exSmiles = res.ex_smiles;

% SHAP matrix for example molecules: 3 x topN
% jsondecode converts JSON 2D array -> MATLAB matrix (n_ex x n_top)
exTopShap = res.ex_top_shap;           % 3 x topN double
exTopXval = res.ex_top_xval;           % 3 x topN double (bit value 0/1)

logInfo("Examples selected: %s (BBB%s, P=%.3f) | %s (BBB%s, P=%.3f) | %s (BBB%s, P=%.3f)", ...
    exNames{1}, labelStr_(exLabels(1)), exProbs(1), ...
    exNames{2}, labelStr_(exLabels(2)), exProbs(2), ...
    exNames{3}, labelStr_(exLabels(3)), exProbs(3));

%% Section 4: Visualization
logSection("RP05", "Section 4: Visualization", "SHAP BBBP Explainability");

% -- Figure 1: Global feature importance (mean |SHAP|, top 20 bits) --
fig1 = figure("Name", "RP05 BBBP: SHAP Global Importance");
set(fig1, "Position", [100 100 700 480]);

barColors = colormap_(topCoef);    % red for positive LR coef, blue for negative
yPos = topN:-1:1;                  % reverse so highest at top

ax1 = axes(fig1);
for k = 1:topN
    barh(ax1, yPos(k), topImp(k), "FaceColor", barColors(k,:), ...
        "EdgeColor", "none", "BarWidth", 0.75);
    hold(ax1, "on");
end
hold(ax1, "off");

yticklabels(ax1, arrayfun(@(i) sprintf("Bit %d", topIdx(i)), ...
    topN:-1:1, "UniformOutput", false));
xlabel(ax1, "Mean |SHAP Value|");
title(ax1, "BBBP ECFP4 LR: Top 20 Bits by Mean |SHAP| (Global Importance)");
subtitle(ax1, "Red = promotes BBB+, Blue = promotes BBB-");
grid(ax1, "on");
box(ax1, "off");

% -- Figure 2: Local waterfall for 3 example molecules --
fig2 = figure("Name", "RP05 BBBP: SHAP Local Explanations");
set(fig2, "Position", [820 100 1000 560]);
% Derive titles from Python-returned ex_types so fallback (TN2) is reflected.
titleStrs = arrayfun(@(k) sprintf("%s: %s  P(BBB+)=%.3f", ...
    typePrefix_(string(exTypes{k})), exNames{k}, exProbs(k)), ...
    1:3, "UniformOutput", false);

for k = 1:3
    sv  = exTopShap(k, :);    % 1 x topN SHAP values
    xv  = exTopXval(k, :);    % 1 x topN bit values (0 or 1)

    % Sort by |SHAP| descending for this molecule
    [~, sortOrd] = sort(abs(sv), "descend");
    showN = min(10, topN);
    ordIdx  = sortOrd(1:showN);
    sv_plot = sv(ordIdx);
    idx_plot = topIdx(ordIdx);

    posColors = repmat([0.85 0.2 0.1], showN, 1);  % red = positive
    negColors = repmat([0.1 0.3 0.75], showN, 1);  % blue = negative
    barFaceColors = posColors;
    barFaceColors(sv_plot < 0, :) = repmat([0.1 0.3 0.75], sum(sv_plot < 0), 1);

    ax = subplot(1, 3, k, "Parent", fig2);
    yp = showN:-1:1;
    for j = 1:showN
        barh(ax, yp(j), sv_plot(j), "FaceColor", barFaceColors(j,:), ...
            "EdgeColor", "none", "BarWidth", 0.75);
        hold(ax, "on");
    end
    hold(ax, "off");

    xline(ax, 0, "k-", "LineWidth", 0.8);
    yticklabels(ax, arrayfun(@(i) sprintf("Bit %d", idx_plot(i)), ...
        showN:-1:1, "UniformOutput", false));
    xlabel(ax, "SHAP Value");
    title(ax, titleStrs{k}, "FontSize", 8);
    grid(ax, "on");
    box(ax, "off");
end
sgtitle(fig2, "RP05 BBBP ECFP4 LR: Local SHAP Explanations (Top 10 bits)");

logInfo("Figures 1-2 created.");

%% Section 5: RF03 Verification
logSection("RP05", "Section 5: RF03 Verification", "SHAP BBBP Explainability");

rf03crit = struct( ...
    "auc_cv",           struct("lower", 0.85), ...
    "shap_lr_spearman", struct("lower", 0.90));

metRP05 = struct( ...
    "auc_cv",           aucCV, ...
    "shap_lr_spearman", shapRho);

resRP05 = emk.repro.verify(metRP05, rf03crit);

logInfo("==> ROC-AUC CV        = %.4f (>= 0.85): %s", ...
    aucCV, statusStr_(resRP05.details.auc_cv.pass));
logInfo("==> SHAP-LR Spearman  = %.4f (>= 0.90): %s", ...
    shapRho, statusStr_(resRP05.details.shap_lr_spearman.pass));
disp(resRP05.report);

if resRP05.pass
    logInfo("==> RP05 REPRODUCTION: PASS");
else
    logWarn("==> RP05 REPRODUCTION: NEEDS REVIEW -- see README.md Discussion");
end

%% Section 6: Save Results
logSection("RP05", "Section 6: Save Results", "SHAP BBBP Explainability");

runDir    = makeRunDir("Prefix", "rp05_shap");
absRunDir = char(fullfile(pwd(), runDir));

% Save figures
saveas(fig1, fullfile(absRunDir, "shap_global_importance.png"));
saveas(fig2, fullfile(absRunDir, "shap_local_waterfall.png"));
logInfo("Figures saved.");

% Metrics JSON
metrics = struct( ...
    "auc_cv",              aucCV, ...
    "auc_cv_std",          aucStd, ...
    "shap_lr_spearman",    shapRho, ...
    "n_valid",             nValid, ...
    "n_train",             nTrain, ...
    "n_test",              nTest, ...
    "n_bbb_pos",           res.n_bbb_pos, ...
    "ecfp4_radius",        2, ...
    "ecfp4_nbits",         2048, ...
    "n_top_bits_shown",    topN, ...
    "rf03_criteria",       rf03crit, ...
    "rf03_pass",           resRP05.pass);
fid = fopen(fullfile(runDir, "metrics.json"), "w");
fprintf(fid, "%s\n", jsonencode(metrics, "PrettyPrint", true));
fclose(fid);
logInfo("Metrics saved: metrics.json");

% Top-N importance table
impTbl = table(topIdx(:), topImp(:), topCoef(:), ...
    VariableNames=["BitIndex", "MeanAbsSHAP", "LR_Coef"]);
impTbl = sortrows(impTbl, "MeanAbsSHAP", "descend");
writetable(impTbl, fullfile(runDir, "top_bits_shap.csv"));
logInfo("Top bits saved: top_bits_shap.csv");

% RF02 lock snapshot
snap.run_date  = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir   = runDir;
snap.rf03_pass = resRP05.pass;
emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));

logInfo("RP05 complete.  run_dir=%s", runDir);

% ===========================================================================
% Local helper functions
% ===========================================================================

function s = statusStr_(tf)
    if tf; s = "PASS"; else; s = "FAIL"; end
end

function s = labelStr_(lbl)
    if lbl == 1; s = "+"; else; s = "-"; end
end

function s = typePrefix_(t)
    switch t
        case "TP";  s = "True Positive";
        case "TN";  s = "True Negative";
        case "TN2"; s = "True Neg (2nd)";
        case "MIS"; s = "Misclassified";
        otherwise
            error("emk:rp05:unknownExType", "Unexpected ex_type from Python: %s", t);
    end
end

function colors = colormap_(coefs)
% Map LR coefficients to red/blue RGB: positive coef -> red, negative -> blue.
    n = numel(coefs);
    colors = zeros(n, 3);
    for k = 1:n
        if coefs(k) >= 0
            colors(k, :) = [0.85 0.2 0.1];   % red: promotes BBB+
        else
            colors(k, :) = [0.1 0.3 0.75];   % blue: promotes BBB-
        end
    end
end
