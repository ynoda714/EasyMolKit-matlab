% r1c_matlab_nested_cv.m  R1-C: MATLAB nested CV with fitclinear lbfgs
%
% Investigates whether MATLAB fitclinear lbfgs can match or beat sklearn
% when given freedom to select its own optimal Lambda.
%
% Section 3a: Fixed-Lambda sweep (lbfgs, same outer folds as sklearn)
%   Scans Lambda grid [100 → 1/n] to find MATLAB's best Lambda for lbfgs.
%
% Section 3b: MATLAB nested CV (apples-to-apples with sklearn)
%   Outer: same 5-fold splits as sklearn (res.outer_fold_indices).
%   Inner: MATLAB stratified 3-fold (independent seed) for Lambda selection.
%   This is the fairest comparison -- MATLAB chooses its own optimal Lambda.
%
% Output: result/runs/r1c_nested_<ts>/r1c_results.json
%
% Run: Ctrl+Enter in MATLAB with project root as CWD.

%% Section 0: Setup
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("R1C", "Section 0: Setup", "MATLAB Nested CV");
emk.setup.initPython();
molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;

%% Section 1: Load BBBP data and compute ECFP4
logSection("R1C", "Section 1: Load BBBP + ECFP4", "MATLAB Nested CV");

tbl    = emk.dataset.bbbp();
csvPath = fullfile(pwd(), "data", "benchmark", "bbbp.csv");
nTotal = height(tbl);
mols   = cell(1, nTotal);
vmask  = false(1, nTotal);
for i = 1:nTotal
    try
        mols{i}  = emk.mol.fromSmiles(tbl.SMILES(i));
        vmask(i) = true;
    catch
    end
end
y = double(tbl.BBB(vmask));
X = batchMorganFP_(mols(vmask), 2, 2048);
logInfo("ECFP4: %d x %d  (BBB+=%.0f, BBB-=%.0f)", ...
    size(X,1), size(X,2), sum(y==1), sum(y==0));

%% Section 2: Python sklearn CV (outer fold splits)
logSection("R1C", "Section 2: Python sklearn CV (outer fold splits)", "MATLAB Nested CV");

coreScript = fullfile(thisDir, "rp02_sklearn_core.py");
resultJson = pyrun([ ...
    "exec(open(script_path).read()); " ...
    "result_json = run_rp02_sklearn(csv_path)"], ...
    "result_json", ...
    script_path=coreScript, csv_path=csvPath);
res = jsondecode(char(string(resultJson)));

if isfield(res, "success") && ~res.success
    error("emk:r1c:pythonError", "sklearn failed: %s", res.error);
end

sklearnAUC = res.auc_mean;
logInfo("sklearn AUC = %.4f +/- %.4f  (best_C per fold: %s)", ...
    sklearnAUC, res.auc_std, mat2str([res.fold_results.best_C], 3));

if size(X,1) ~= res.n_valid
    error("emk:r1c:sizeMismatch", "MATLAB n_valid=%d vs Python n_valid=%d", ...
        size(X,1), res.n_valid);
end

nFolds = numel(res.fold_results);

%% Section 3a: Fixed-Lambda sweep with lbfgs (same outer folds)
logSection("R1C", "Section 3a: Fixed-Lambda lbfgs sweep", "MATLAB Nested CV");
%
% Lambda grid covers sklearn C grid (Lambda=1/C) plus very weak regularization.
% C grid used in sklearn: [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0]
% Equivalent Lambda:      [100,  20,   10,  2,   1,   0.5, 0.2]
% Plus Lambda=1/n (≈ 0.0005) = sklearn C≈2039, historical default.
n = size(X, 1);
LAMBDA_GRID = [100, 20, 10, 2, 1, 0.5, 0.2, 1/n];

logInfo("Lambda grid: %s  (equiv C: %s)", ...
    mat2str(LAMBDA_GRID, 3), mat2str(1./LAMBDA_GRID, 3));
logInfo("");

sweepAUCs = zeros(numel(LAMBDA_GRID), nFolds);
for li = 1:numel(LAMBDA_GRID)
    lam = LAMBDA_GRID(li);
    for k = 1:nFolds
        trIdx = res.outer_fold_indices.(sprintf("fold%d_train", k)) + 1;
        teIdx = res.outer_fold_indices.(sprintf("fold%d_test",  k)) + 1;
        mdl   = fitclinear(X(trIdx,:), y(trIdx), ...
            "Learner", "logistic", "Regularization", "ridge", ...
            "Lambda",  lam, "Solver", "lbfgs");
        pIdx = find(mdl.ClassNames == 1, 1);
        [~, sc] = predict(mdl, X(teIdx,:));
        [~,~,~, sweepAUCs(li,k)] = perfcurve(y(teIdx), sc(:,pIdx), 1);
    end
    mu  = mean(sweepAUCs(li,:));
    sig = std(sweepAUCs(li,:));
    cEq = 1/lam;
    marker = "";
    if abs(lam - 10) < 1e-9;  marker = "  <- R1-B result";  end
    if abs(lam - 1/n) < 1e-9; marker = "  <- historical default (Lambda=1/n)"; end
    logInfo("  Lambda=%-8.4g (C≈%-6.3g): AUC=%.4f +/-%.4f  gap=%+.4f%s", ...
        lam, cEq, mu, sig, sklearnAUC - mu, marker);
end

[bestSweepAUC, bestLi] = max(mean(sweepAUCs, 2));
bestLambdaFixed = LAMBDA_GRID(bestLi);
logInfo("");
logInfo("Best fixed Lambda: %.4g  AUC=%.4f  gap=%+.4f", ...
    bestLambdaFixed, bestSweepAUC, sklearnAUC - bestSweepAUC);

%% Section 3b: MATLAB nested CV (fitclinear lbfgs, free Lambda selection)
logSection("R1C", "Section 3b: MATLAB nested CV (lbfgs, free Lambda)", "MATLAB Nested CV");
%
% Outer: same splits as sklearn (from res.outer_fold_indices).
% Inner: MATLAB stratified 3-fold on outer-train (independent seed=7).
% Lambda candidates: full LAMBDA_GRID from Section 3a.
%
% This is the fairest comparison: MATLAB chooses its own Lambda without
% any constraint from sklearn's C_opt.

INNER_SEED  = 7;
N_INNER     = 3;

aucsNested  = zeros(nFolds, 1);
bestLamNested = zeros(nFolds, 1);
innerAUCNested = zeros(nFolds, 1);

for k = 1:nFolds
    trIdx = res.outer_fold_indices.(sprintf("fold%d_train", k)) + 1;
    teIdx = res.outer_fold_indices.(sprintf("fold%d_test",  k)) + 1;
    Xtr = X(trIdx,:);
    ytr = y(trIdx);

    % Inner 3-fold stratified CV on outer-train to select Lambda
    rng(INNER_SEED, "twister");
    innerCV = cvpartition(ytr, "KFold", N_INNER, "Stratify", true);

    bestInnerAUC = -Inf;
    bestLam      = LAMBDA_GRID(1);
    for li = 1:numel(LAMBDA_GRID)
        lam = LAMBDA_GRID(li);
        innerAUCs = zeros(N_INNER, 1);
        for ki = 1:N_INNER
            trI = training(innerCV, ki);
            teI = test(innerCV,     ki);
            mi  = fitclinear(Xtr(trI,:), ytr(trI), ...
                "Learner", "logistic", "Regularization", "ridge", ...
                "Lambda",  lam, "Solver", "lbfgs");
            pI = find(mi.ClassNames == 1, 1);
            [~, sci] = predict(mi, Xtr(teI,:));
            [~,~,~, innerAUCs(ki)] = perfcurve(ytr(teI), sci(:,pI), 1);
        end
        muInner = mean(innerAUCs);
        if muInner > bestInnerAUC
            bestInnerAUC = muInner;
            bestLam      = lam;
        end
    end

    % Fit on full outer train with best Lambda
    mdl  = fitclinear(Xtr, ytr, ...
        "Learner", "logistic", "Regularization", "ridge", ...
        "Lambda",  bestLam, "Solver", "lbfgs");
    pIdx = find(mdl.ClassNames == 1, 1);
    [~, sc] = predict(mdl, X(teIdx,:));
    [~,~,~, aucK] = perfcurve(y(teIdx), sc(:,pIdx), 1);

    aucsNested(k)     = aucK;
    bestLamNested(k)  = bestLam;
    innerAUCNested(k) = bestInnerAUC;
    logInfo("  fold %d: AUC=%.4f  best_Lambda=%.4g (C≈%.3g)  inner_AUC=%.4f  sklearn=%.4f", ...
        k, aucK, bestLam, 1/bestLam, bestInnerAUC, res.fold_results(k).test_auc);
end

aucNestedMean = mean(aucsNested);
aucNestedStd  = std(aucsNested);
logInfo("");
logInfo("MATLAB nested CV (lbfgs): AUC = %.4f +/- %.4f", aucNestedMean, aucNestedStd);
logInfo("best_Lambda per fold: %s", mat2str(bestLamNested', 3));

%% Section 4: Summary
logSection("R1C", "Section 4: Summary", "MATLAB Nested CV");

SKLEARN_AUC         = sklearnAUC;
R1B_LBFGS_AUC       = 0.8417;
FITCLINEAR_HIST_AUC = 0.8826;

logInfo("======= R1-C SUMMARY =======");
logInfo("  sklearn lbfgs nested CV                      : AUC = %.4f  [target]", SKLEARN_AUC);
logInfo("  fitclinear lbfgs Lambda=1/C_opt (R1-B)       : AUC = %.4f  (gap=%+.4f)", ...
    R1B_LBFGS_AUC, SKLEARN_AUC - R1B_LBFGS_AUC);
logInfo("  fitclinear lbfgs best-fixed Lambda (3a)       : AUC = %.4f  (gap=%+.4f)", ...
    bestSweepAUC, SKLEARN_AUC - bestSweepAUC);
logInfo("  fitclinear lbfgs nested CV / free Lambda (3b) : AUC = %.4f  (gap=%+.4f)", ...
    aucNestedMean, SKLEARN_AUC - aucNestedMean);
logInfo("  fitclinear SGD   Lambda=1/n (historical RP02) : AUC = %.4f  (gap=%+.4f)", ...
    FITCLINEAR_HIST_AUC, SKLEARN_AUC - FITCLINEAR_HIST_AUC);
logInfo("");

gap = SKLEARN_AUC - aucNestedMean;
if gap <= 0
    logInfo("  >> MATLAB WINS: fitclinear nested CV AUC %.4f > sklearn %.4f", ...
        aucNestedMean, SKLEARN_AUC);
elseif gap < 0.01
    logInfo("  >> PRACTICAL TIE: gap = %.4f (< 0.01)", gap);
elseif gap < 0.03
    logInfo("  >> GAP NARROWS: residual = %.4f -- marginal sklearn advantage", gap);
else
    logInfo("  >> SKLEARN WINS: residual gap = %.4f -- sklearn consistently superior", gap);
end

%% Section 5: Save results
logSection("R1C", "Section 5: Save Results", "MATLAB Nested CV");

runDir = makeRunDir("Prefix", "r1c_nested");
if startsWith(runDir, '/') || startsWith(runDir, '\\') || ...
        (numel(runDir) >= 2 && runDir(2) == ':')
    absRunDir = char(runDir);
else
    absRunDir = char(fullfile(pwd(), runDir));
end

output = struct( ...
    "sklearn_auc",                SKLEARN_AUC, ...
    "r1b_lbfgs_fixed_lambda_auc", R1B_LBFGS_AUC, ...
    "fitclinear_hist_auc",        FITCLINEAR_HIST_AUC, ...
    "section3a", struct( ...
        "lambda_grid",      LAMBDA_GRID, ...
        "auc_per_lambda",   mean(sweepAUCs, 2)', ...
        "std_per_lambda",   std(sweepAUCs, 0, 2)', ...
        "best_lambda",      bestLambdaFixed, ...
        "best_auc",         bestSweepAUC), ...
    "section3b", struct( ...
        "auc_cv",           aucNestedMean, ...
        "auc_cv_std",       aucNestedStd, ...
        "auc_per_fold",     aucsNested', ...
        "best_lambda_per_fold", bestLamNested', ...
        "inner_auc_per_fold",   innerAUCNested', ...
        "sklearn_gap",      SKLEARN_AUC - aucNestedMean));

fid = fopen(fullfile(absRunDir, "r1c_results.json"), "w");
if fid == -1
    error("emk:r1c:fopenFailed", "Cannot open: %s", absRunDir);
end
fprintf(fid, "%s\n", jsonencode(output, "PrettyPrint", true));
fclose(fid);
logInfo("Results saved: r1c_results.json  (run_dir=%s)", runDir);

% ===========================================================================
function X = batchMorganFP_(mols, radius, nBits)
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
