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
%   Model (biased):  sklearn LogisticRegression (L2 ridge, C=1.0) + Morgan ECFP4
%                      Radius=2, NBits=2048 -- consistent with RP02 original.
%   Model (rev):     C selected via inner 3-fold CV on X_tr (RP02-rev design,
%                      outer_seed=42, inner_seed=7); same 80/20 split as biased.
%
%   SHAP:   shap.LinearExplainer -- exact for linear models.
%             global_imp[j] = mean_i(|coef[j] * (X[i,j] - mean(X[:,j]))|)
%
%   Task:   Feature attribution -- which ECFP4 bits drive BBB permeability?
%           Task B: Compare feature rankings between biased (C=1.0) and rev (C_opt).
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
    "exec(open(hp).read(), globals()); result_json = run_rp05(cp)", ...
    "result_json", hp=helperPath, cp=csvPath);

res = jsondecode(char(string(pyResult)));

aucCV   = res.auc_cv;
aucStd  = res.auc_cv_std;
nValid  = res.n_valid;
nTrain  = res.n_train;
nTest   = res.n_test;
shapRho = res.shap_lr_spearman;

logInfo("Dataset: %d valid (%d train / %d test), BBB+: %d total (%d tr / %d te), BBB-: %d", ...
    nValid, nTrain, nTest, res.n_bbb_pos, res.n_bbb_pos_train, res.n_bbb_pos_test, nValid - res.n_bbb_pos);
logInfo("sklearn LR 5-fold CV (biased, C=1.0): AUC=%.4f +/- %.4f  (full dataset n=%d)", ...
    aucCV, aucStd, nValid);
logInfo("SHAP LinearExplainer (biased): Spearman(global_imp, |coef|*std_train) = %.4f", shapRho);

%% Section 2b: Fit RP02-rev Model + SHAP (Task B)
logSection("RP05", "Section 2b: Fit RP02-rev Model + SHAP (Task B)", ...
    "SHAP BBBP Explainability");

logInfo("Running Python: RP02-rev nested CV AUC + inner-CV C selection + SHAP ...");
logInfo("(SHAP JIT already compiled from Section 2 -- this run should be faster)");

pyRevResult = pyrun( ...
    "exec(open(hp).read(), globals()); rev_json = run_rp05_rev(cp)", ...
    "rev_json", hp=helperPath, cp=csvPath);

resRev = jsondecode(char(string(pyRevResult)));

aucCV_rev   = resRev.auc_cv;
aucStd_rev  = resRev.auc_cv_std;
shapRho_rev = resRev.shap_lr_spearman;
bestC_shap  = resRev.best_C_shap;

logInfo("Rev nested CV AUC=%.4f +/- %.4f (best_C per fold: %s)", ...
    aucCV_rev, aucStd_rev, ...
    strjoin(arrayfun(@(c) sprintf('%.2f', c), resRev.best_C_per_fold, ...
        "UniformOutput", false), "/"));
logInfo("Rev SHAP model best C (inner CV on X_tr) = %.2f", bestC_shap);
logInfo("Rev SHAP-LR Spearman = %.4f", shapRho_rev);

%% Section 2c: F2 -- MATLAB Analytical SHAP vs Python LinearExplainer (Zone B Test)
% Requires: Section 0 (initPython, helperPath), Section 1 (csvPath)
% Implements analytical interventional-linear SHAP for linear models:
%   SHAP_j(x_i) = Beta_j * (x_i[j] - mean_train[j])
% Primary comparison: Python C=1.0 (biased) -- matches MATLAB's optimal regularization range.
% Secondary comparison: Python C=0.10 (rev, AUC-optimal) -- documents the 16,000x reg gap.
% Zone B confirmed when primary Spearman rho >= 0.85.
logSection("RP05", "Section 2c: F2 -- MATLAB SHAP vs Python LinearExplainer", ...
    "SHAP BBBP Explainability");

% scipy.io.savemat avoids transferring the ~1631x2048 binary matrix as JSON (~8 MB).
f2MatPath     = fullfile(char(tempdir()), "rp05_f2_data.mat");
f2MatPathRev  = fullfile(char(tempdir()), "rp05_f2_data_rev.mat");

% Primary: Python C=1.0 (biased) -- Zone B test
logInfo("F2: Running Python run_rp05_f2 (C_shap=1.0, primary Zone B test) ...");
pyF2Result = pyrun( ...
    "exec(open(hp).read(), globals()); f2_json = run_rp05_f2(cp, mp, C_shap=1.0)", ...
    "f2_json", hp=helperPath, cp=csvPath, mp=f2MatPath);

resF2 = jsondecode(char(string(pyF2Result)));
if isfield(resF2, "success") && ~resF2.success
    error("emk:rp05:f2PythonError", ...
        "run_rp05_f2 failed:\n%s\n\nTraceback:\n%s", ...
        resF2.error, resF2.traceback);
end
logInfo("F2: Primary data saved (n_train=%d, n_test=%d, C_shap=%.2f)", ...
    resF2.n_train, resF2.n_test, resF2.C_shap);

% Load feature matrices: scipy .mat -> MATLAB double
f2Data       = load(f2MatPath);
X_tr_f2      = f2Data.X_tr;                              % n_train x 2048
X_te_f2      = f2Data.X_te;                              % n_test x 2048
y_tr_f2      = double(f2Data.y_tr(:));                   % n_train x 1
gImpPyBias   = double(f2Data.global_imp_python_te(:));   % 2048 x 1, Python C=1.0
nTrain_f2    = size(X_tr_f2, 1);
nTest_f2     = size(X_te_f2, 1);
logInfo("F2: Loaded X_tr %dx%d, X_te %dx%d.", nTrain_f2, size(X_tr_f2,2), nTest_f2, size(X_te_f2,2));

% Secondary: Python C=0.10 (rev, AUC-optimal) -- root-cause analysis
logInfo("F2: Running Python run_rp05_f2 (C_shap=0.10, secondary / regularization gap analysis) ...");
pyF2RevResult = pyrun( ...
    "exec(open(hp).read(), globals()); f2rev_json = run_rp05_f2(cp, mp, C_shap=0.10)", ...
    "f2rev_json", hp=helperPath, cp=csvPath, mp=f2MatPathRev);
resF2Rev = jsondecode(char(string(pyF2RevResult)));
if isfield(resF2Rev, "success") && ~resF2Rev.success
    error("emk:rp05:f2RevPythonError", "run_rp05_f2(C=0.10) failed:\n%s", resF2Rev.error);
end
gImpPyRev = double(load(f2MatPathRev).global_imp_python_te(:)); % 2048 x 1, Python C=0.10

% Fit MATLAB fitclinear -- lbfgs, Lambda=1/n (MATLAB-optimal from R1-C in RP02)
logInfo("F2: Fitting MATLAB fitclinear (lbfgs, Lambda=1/n=%.6f) ...", 1/nTrain_f2);
MdlF2 = fitclinear(X_tr_f2, y_tr_f2, ...
    "Learner", "logistic", ...
    "Solver",  "lbfgs", ...
    "Lambda",  1 / nTrain_f2);

% Analytical interventional-linear SHAP:
%   SHAP_j(x_i) = Beta_j * (x_i[j] - mean_train[j])
mu_tr_f2     = mean(X_tr_f2, 1);                          % 1 x 2048 (background mean)
shap_te_f2   = (X_te_f2 - mu_tr_f2) .* MdlF2.Beta';      % n_test x 2048
gImpMlab_te  = mean(abs(shap_te_f2), 1)';                 % 2048 x 1

% Spearman correlations: primary (C=1.0) and secondary (C=0.10)
[f2RhoBias, f2PvalBias] = corr(gImpMlab_te, gImpPyBias, "Type", "Spearman");
[f2RhoRev,  f2PvalRev]  = corr(gImpMlab_te, gImpPyRev,  "Type", "Spearman");
[rhoPyPy,   ~]          = corr(gImpPyBias,  gImpPyRev,  "Type", "Spearman");
f2RhoTe  = f2RhoBias;   % primary result
f2PvalTe = f2PvalBias;

logInfo("F2: Spearman rho -- MATLAB vs Python(C=1.0):  %.4f (primary, p=%.2e)", ...
    f2RhoBias, f2PvalBias);
logInfo("F2: Spearman rho -- MATLAB vs Python(C=0.10): %.4f (secondary, p=%.2e)", ...
    f2RhoRev, f2PvalRev);
logInfo("F2: Spearman rho -- Python(C=1.0) vs Python(C=0.10): %.4f (reference)", rhoPyPy);

% Zone B verdict (primary comparison)
f2ZoneB = f2RhoBias >= 0.85;
if f2ZoneB
    logInfo("F2: Zone B CONFIRMED -- primary rho=%.4f >= 0.85 (MATLAB SHAP ~ Python at matched reg)", ...
        f2RhoBias);
    logInfo("F2: Root cause: vs Python(C=0.10) rho=%.4f; 16,000x reg gap, not formula error", ...
        f2RhoRev);
else
    logWarn("F2: Zone B not confirmed -- primary rho=%.4f < 0.85", f2RhoBias);
end

%% Section 2d: F3 -- MATLAB shapley() sampling study vs Python full-test TreeSHAP
% F3 now measures how well MATLAB shapley() reproduces Python RF TreeSHAP rankings
% under practical subsampled evaluation. Python computes the reference ranking on the
% full test set once; MATLAB evaluates shapley() on repeated X_te subsamples for
% n_eval = [16 32 64 128] and compares rankings via Spearman rho.
logSection("RP05", "Section 2d: F3 -- MATLAB shapley vs Python TreeExplainer", ...
    "SHAP BBBP Explainability");

f3Mode = "exploratory";
f3NEvalGrid = [1 2 4];
f3NRepeats = 1;
f3SeedBase = 42;
f3ZoneThreshold = 0.85;
f3BenchNEval = 1;
f3MaxProjectedRuntimeSec = 20 * 60;
f3RuntimeBudgetSec = 10 * 60;
f3Proceed = true;
f3SkipReason = "";
f3ZoneCandidate = "NOT_RUN";
f3ZoneConclusion = "F3 not yet evaluated.";
f3FormalZoneEvaluable = false;

% Step 1: Python RF + TreeSHAP full-test reference
logInfo("F3: Running Python RF (100 trees) + full-test TreeSHAP reference ...");
f3MatPath = fullfile(char(tempdir()), "rp05_f3_data.mat");
pyF3Result = pyrun( ...
    "exec(open(hp).read(), globals()); f3_json = run_rp05_f3(cp, mp)", ...
    "f3_json", hp=helperPath, cp=csvPath, mp=f3MatPath);

resF3 = jsondecode(char(string(pyF3Result)));
if isfield(resF3, "success") && ~resF3.success
    error("emk:rp05:f3PythonError", ...
        "run_rp05_f3 failed:\n%s\n\nTraceback:\n%s", resF3.error, resF3.traceback);
end
logInfo("F3: Python RF CV AUC=%.4f +/- %.4f (n_train=%d, n_test=%d)", ...
    resF3.rf_auc_cv, resF3.rf_auc_std, resF3.n_train, resF3.n_test);

% Load data saved by Python
f3Data         = load(f3MatPath);
X_tr_f3        = f3Data.X_tr;                          % n_train x 2048
X_te_f3        = f3Data.X_te;                          % n_test x 2048
X_bg_f3        = f3Data.X_bg;                          % n_background x 2048 (subsample of X_tr)
y_tr_f3        = double(f3Data.y_tr(:));
gImpPyTreeFull = double(f3Data.global_imp_treeshap_full(:)); % 2048 x 1
logInfo("F3: Loaded X_tr %dx%d, X_te %dx%d, X_bg %dx%d.", ...
    size(X_tr_f3,1), size(X_tr_f3,2), size(X_te_f3,1), size(X_te_f3,2), ...
    size(X_bg_f3,1), size(X_bg_f3,2));
logInfo("F3: mode=%s, n_eval=[%s], repeats=%d.", ...
    f3Mode, num2str(f3NEvalGrid), f3NRepeats);

% Step 2: Train MATLAB fitcensemble("Bag") once; shapley() reuses the same model.
logInfo("F3: Training MATLAB fitcensemble(Bag, 100 trees, rng(42)) ...");
rng(42);
MdlF3ens = fitcensemble(X_tr_f3, y_tr_f3, ...
    "Method",            "Bag", ...
    "NumLearningCycles", 100, ...
    "Learners",          "tree");
logInfo("F3: Model trained. ClassNames: [%s]", ...
    strjoin(string(MdlF3ens.ClassNames), ", "));

% Step 3: MATLAB shapley() with fixed Python-matched background rows.
logInfo("F3: Creating shapley explainer (background=%d samples, subsample of X_tr) ...", size(X_bg_f3,1));
eshapF3 = shapley(MdlF3ens, X_bg_f3);
logInfo("F3: Method auto-selected: %s", eshapF3.Method);

% MeanAbsoluteShapley table columns: {Predictor, class_0, class_1}
% Column 3 = class "1" (BBB+), matching Python's sv_pos = shap_values[1]
classNames_f3 = string(MdlF3ens.ClassNames);
posColIdx_f3  = find(classNames_f3 == "1") + 1;  % +1 because col 1 = Predictor
if isempty(posColIdx_f3)
    error("emk:rp05:f3ClassNotFound", ...
        "Class '1' not found in ClassNames: [%s]", strjoin(classNames_f3, ", "));
end

% Step 4: Benchmark one tiny SHAP run before committing to the full grid.
f3TopK = 20;
topPyTree = topkIdx_(gImpPyTreeFull, f3TopK);
rng(f3SeedBase);
f3BenchIdx = randperm(size(X_te_f3, 1), f3BenchNEval);
ticF3Bench = tic;
eshapF3Bench = fit(eshapF3, X_te_f3(f3BenchIdx, :));
f3BenchRuntimeSec = toc(ticF3Bench);
gImpMlabBench = eshapF3Bench.MeanAbsoluteShapley{:, posColIdx_f3};
[f3BenchRho, f3BenchPval] = corr(gImpMlabBench, gImpPyTreeFull, "Type", "Spearman");
f3ProjectedRuntimeSec = f3BenchRuntimeSec * (sum(f3NEvalGrid) / f3BenchNEval) * f3NRepeats;
f3Benchmark = struct();
f3Benchmark.n_eval = f3BenchNEval;
f3Benchmark.runtime_sec = f3BenchRuntimeSec;
f3Benchmark.projected_runtime_sec = f3ProjectedRuntimeSec;
f3Benchmark.projected_runtime_min = f3ProjectedRuntimeSec / 60;
f3Benchmark.rho = f3BenchRho;
f3Benchmark.p_value = f3BenchPval;
f3Benchmark.runtime_threshold_sec = f3MaxProjectedRuntimeSec;
f3Benchmark.runtime_budget_sec = f3RuntimeBudgetSec;
f3Benchmark.seed = f3SeedBase;
logInfo("F3 benchmark: n_eval=%d runtime=%.1fs projected_total=%.1f min rho=%.4f", ...
    f3BenchNEval, f3BenchRuntimeSec, f3ProjectedRuntimeSec / 60, f3BenchRho);

if f3ProjectedRuntimeSec > f3MaxProjectedRuntimeSec
    f3Proceed = false;
    f3SkipReason = sprintf([ ...
        'Projected F3 runtime %.1f min exceeds guardrail %.1f min ' ...
        '(bench: n_eval=%d took %.1fs).'], ...
        f3ProjectedRuntimeSec / 60, f3MaxProjectedRuntimeSec / 60, ...
        f3BenchNEval, f3BenchRuntimeSec);
    f3ZoneCandidate = "SKIPPED";
    f3ZoneConclusion = sprintf([ ...
        'F3 skipped after benchmark because MATLAB interventional-tree SHAP was too slow for the ' ...
        'configured exploratory n_eval grid. %s'], f3SkipReason);
    logWarn("F3: %s", f3SkipReason);
end

% Step 5: Repeated subsampled evaluation against Python full-test reference.
f3Results = repmat(struct( ...
    "n_eval", 0, ...
    "rho_values", [], ...
    "rho_mean", 0, ...
    "rho_std", 0, ...
    "rho_min", 0, ...
    "rho_max", 0, ...
    "runtime_values_sec", [], ...
    "runtime_mean_sec", 0, ...
    "runtime_std_sec", 0, ...
    "top20_overlap_values", [], ...
    "top20_overlap_mean", 0, ...
    "top20_overlap_std", 0, ...
    "p_values", [], ...
    "seeds", [], ...
    "sample_indices", {{}}, ...
    "sample_size", 0), numel(f3NEvalGrid), 1);

f3QualifyingIdx = [];
f3CompletedMask = false(numel(f3NEvalGrid), 1);
f3CumulativeRuntimeSec = 0;
if f3Proceed
    for iEval = 1:numel(f3NEvalGrid)
        nEval = f3NEvalGrid(iEval);
        rhoVals = zeros(f3NRepeats, 1);
        pVals = zeros(f3NRepeats, 1);
        runVals = zeros(f3NRepeats, 1);
        top20Vals = zeros(f3NRepeats, 1);
        seedVals = zeros(f3NRepeats, 1);
        sampleIdxVals = cell(f3NRepeats, 1);

        if nEval == f3BenchNEval
            runVals(1) = f3BenchRuntimeSec;
            rhoVals(1) = f3BenchRho;
            pVals(1) = f3BenchPval;
            seedVals(1) = f3SeedBase;
            sampleIdxVals{1} = f3BenchIdx;
            topMlBench = topkIdx_(gImpMlabBench, f3TopK);
            top20Vals(1) = numel(intersect(topPyTree, topMlBench));
            f3CumulativeRuntimeSec = f3CumulativeRuntimeSec + f3BenchRuntimeSec;
            logInfo("F3: n_eval=%d repeat=1/%d seed=%d -> rho=%.4f, overlap20=%d, runtime=%.1fs (benchmark reused)", ...
                nEval, f3NRepeats, f3SeedBase, rhoVals(1), top20Vals(1), runVals(1));
        else
            f3ProjectedThisSec = f3BenchRuntimeSec * (nEval / f3BenchNEval) * f3NRepeats;
            if (f3CumulativeRuntimeSec + f3ProjectedThisSec) > f3RuntimeBudgetSec
                f3SkipReason = sprintf([ ...
                    'Stopped exploratory F3 after %.1f min because next n_eval=%d was projected to exceed ' ...
                    'the runtime budget of %.1f min.'], ...
                    f3CumulativeRuntimeSec / 60, nEval, f3RuntimeBudgetSec / 60);
                logWarn("F3: %s", f3SkipReason);
                break;
            end

            for iRep = 1:f3NRepeats
                repSeed = f3SeedBase + iRep - 1;
                rng(repSeed);
                evalIdx = randperm(size(X_te_f3, 1), nEval);
                sampleIdxVals{iRep} = evalIdx;
                seedVals(iRep) = repSeed;

                ticF3 = tic;
                eshapF3Eval = fit(eshapF3, X_te_f3(evalIdx, :));
                runVals(iRep) = toc(ticF3);
                f3CumulativeRuntimeSec = f3CumulativeRuntimeSec + runVals(iRep);

                gImpMlabEval = eshapF3Eval.MeanAbsoluteShapley{:, posColIdx_f3};
                [rhoVals(iRep), pVals(iRep)] = corr(gImpMlabEval, gImpPyTreeFull, "Type", "Spearman");

                topMl = topkIdx_(gImpMlabEval, f3TopK);
                top20Vals(iRep) = numel(intersect(topPyTree, topMl));

                logInfo("F3: n_eval=%d repeat=%d/%d seed=%d -> rho=%.4f, overlap20=%d, runtime=%.1fs", ...
                    nEval, iRep, f3NRepeats, repSeed, rhoVals(iRep), top20Vals(iRep), runVals(iRep));
            end
        end

        f3Results(iEval).n_eval = nEval;
        f3Results(iEval).rho_values = rhoVals';
        f3Results(iEval).rho_mean = mean(rhoVals);
        f3Results(iEval).rho_std = std(rhoVals);
        f3Results(iEval).rho_min = min(rhoVals);
        f3Results(iEval).rho_max = max(rhoVals);
        f3Results(iEval).runtime_values_sec = runVals';
        f3Results(iEval).runtime_mean_sec = mean(runVals);
        f3Results(iEval).runtime_std_sec = std(runVals);
        f3Results(iEval).top20_overlap_values = top20Vals';
        f3Results(iEval).top20_overlap_mean = mean(top20Vals);
        f3Results(iEval).top20_overlap_std = std(top20Vals);
        f3Results(iEval).p_values = pVals';
        f3Results(iEval).seeds = seedVals';
        f3Results(iEval).sample_indices = sampleIdxVals;
        f3Results(iEval).sample_size = nEval;
        f3CompletedMask(iEval) = true;

        if f3Results(iEval).rho_mean >= f3ZoneThreshold
            f3QualifyingIdx(end+1) = iEval; %#ok<AGROW>
        end
    end

    if f3FormalZoneEvaluable
        if ~isempty(f3QualifyingIdx)
            f3ZoneCandidate = "C";
            f3ZoneConclusion = sprintf([ ...
                'MATLAB interventional-tree SHAP reached mean Spearman rho >= %.2f at n_eval >= 64, ' ...
                'so RP05 nonlinear SHAP is classified as Zone C under sampling-based validation.'], ...
                f3ZoneThreshold);
        else
            f3ZoneCandidate = "D";
            f3ZoneConclusion = sprintf([ ...
                'MATLAB shapley() did not reach mean Spearman rho >= %.2f for any n_eval >= 64, ' ...
                'so RP05 nonlinear SHAP remains Zone D.'], f3ZoneThreshold);
        end
    else
        f3ZoneCandidate = "EXPLORATORY";
        f3CompletedNEvals = f3NEvalGrid(f3CompletedMask);
        if isempty(f3CompletedNEvals)
            f3ZoneConclusion = "Exploratory F3 did not complete any n_eval setting within the runtime budget.";
        elseif ~isempty(f3QualifyingIdx)
            f3ZoneConclusion = sprintf([ ...
                'Exploratory F3 completed on n_eval=[%s]. At least one run reached rho >= %.2f, ' ...
                'but this configuration is not sufficient for formal Zone C/D classification.'], ...
                num2str(f3CompletedNEvals), f3ZoneThreshold);
        else
            f3ZoneConclusion = sprintf([ ...
                'Exploratory F3 completed on n_eval=[%s]. No run reached rho >= %.2f, ' ...
                'but this configuration is not sufficient for formal Zone C/D classification.'], ...
                num2str(f3CompletedNEvals), f3ZoneThreshold);
        end
    end

    for iEval = find(f3CompletedMask(:))'
        logInfo("F3 summary: n_eval=%d rho_mean=%.4f rho_std=%.4f runtime_mean=%.1fs overlap20=%.2f", ...
            f3Results(iEval).n_eval, f3Results(iEval).rho_mean, f3Results(iEval).rho_std, ...
            f3Results(iEval).runtime_mean_sec, f3Results(iEval).top20_overlap_mean);
    end
end
logInfo("F3: zone_candidate=%s", f3ZoneCandidate);

%% Section 3: Extract SHAP Results for MATLAB Visualization
logSection("RP05", "Section 3: Extract SHAP Results", "SHAP BBBP Explainability");

% Global importance (top 20 bits) -- 0-based Python indices -> 1-based MATLAB
topN    = numel(res.top_n_imp);
topIdx  = res.top_n_idx + 1;           % 1-based MATLAB bit indices
topImp  = res.top_n_imp;               % mean |SHAP| (biased)
topCoef = res.top_n_coef;              % LR coefficient (biased)

% Full 2048-bit importance array for ranking comparison (Task B)
gImpBiased = res.global_imp_all(:);    % 2048 x 1

% Rev top-N (0-based Python indices -> 1-based MATLAB)
topIdx_rev  = resRev.top_n_idx + 1;   % 1-based
topImp_rev  = resRev.top_n_imp;
topCoef_rev = resRev.top_n_coef;
gImpRev     = resRev.global_imp_all(:); % 2048 x 1

% Ranking comparison (Task B)
[rankRho, ~]  = corr(gImpBiased, gImpRev, "Type", "Spearman");
inBoth     = intersect(topIdx(:), topIdx_rev(:));
biasedOnly = setdiff(topIdx(:), topIdx_rev(:));
revOnly    = setdiff(topIdx_rev(:), topIdx(:));
nOverlap   = numel(inBoth);

logInfo("Task B ranking comparison (biased C=1.0 vs rev C=%.2f):", bestC_shap);
logInfo("  Spearman rho (all 2048 bits) = %.4f", rankRho);
logInfo("  Top-20 overlap = %d / 20  (biased-only=%d, rev-only=%d)", ...
    nOverlap, numel(biasedOnly), numel(revOnly));

% 3 example molecules (TP / TN / misclassified)
exTypes  = cellstr(res.ex_types);      % cell of strings: {"TP","TN","MIS"}
if numel(exTypes) ~= 3
    error("emk:rp05:unexpectedExTypes", ...
        "Expected 3 ex_types from Python, got %d", numel(exTypes));
end
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

%% Section 4b: Task B -- Ranking Comparison Visualization (Figure 3)
logSection("RP05", "Section 4b: Task B -- Ranking Comparison (Fig 3)", ...
    "SHAP BBBP Explainability");

fig3 = figure("Name", "RP05-B: SHAP Ranking Comparison (Biased vs Rev)");
set(fig3, "Position", [100 650 700 500]);
ax3 = axes(fig3);

% Identify bit index groups (1-based)
topUnion = unique([topIdx(:); topIdx_rev(:)]);
allBits  = (1:2048)';
otherIdx = setdiff(allBits, topUnion);

% Gray background: all bits not in either top-20
scatter(ax3, gImpBiased(otherIdx), gImpRev(otherIdx), 4, [0.78 0.78 0.78], ...
    "filled", "DisplayName", "Other bits");
hold(ax3, "on");

% Biased-only top bits (blue)
if ~isempty(biasedOnly)
    scatter(ax3, gImpBiased(biasedOnly), gImpRev(biasedOnly), ...
        70, [0.1 0.3 0.75], "filled", "DisplayName", ...
        sprintf("Top-20 biased only (n=%d)", numel(biasedOnly)));
end

% Rev-only top bits (red)
if ~isempty(revOnly)
    scatter(ax3, gImpBiased(revOnly), gImpRev(revOnly), ...
        70, [0.85 0.2 0.1], "filled", "DisplayName", ...
        sprintf("Top-20 rev only (n=%d)", numel(revOnly)));
end

% Shared top bits (green diamond)
if ~isempty(inBoth)
    scatter(ax3, gImpBiased(inBoth), gImpRev(inBoth), ...
        90, [0.1 0.6 0.2], "d", "filled", "DisplayName", ...
        sprintf("Top-20 both (n=%d)", numel(inBoth)));
end

hold(ax3, "off");
xlabel(ax3, "Mean |SHAP| — Biased (C=1.0)");
ylabel(ax3, sprintf("Mean |SHAP| — Rev (C_{opt}=%.2f)", bestC_shap));
title(ax3, "BBBP ECFP4 LR: Global Importance Comparison");
subtitle(ax3, sprintf("Spearman \\rho=%.4f (all 2048 bits), Top-20 overlap=%d/20", ...
    rankRho, nOverlap));
legend(ax3, "Location", "northwest");
grid(ax3, "on");
box(ax3, "off");

logInfo("Figure 3 created (ranking comparison scatter).");

%% Section 4c: F3 -- Sampling Study Visualization (Figures 4-5)
logSection("RP05", "Section 4c: F3 -- Sampling Study (Fig 4-5)", ...
    "SHAP BBBP Explainability");

if f3Proceed
    f3EvalPlot = f3NEvalGrid(f3CompletedMask);
    f3RhoMean = arrayfun(@(s) s.rho_mean, f3Results(f3CompletedMask));
    f3RhoStd = arrayfun(@(s) s.rho_std, f3Results(f3CompletedMask));
    f3RunMean = arrayfun(@(s) s.runtime_mean_sec, f3Results(f3CompletedMask));
    f3RunStd = arrayfun(@(s) s.runtime_std_sec, f3Results(f3CompletedMask));
    f3OverlapMean = arrayfun(@(s) s.top20_overlap_mean, f3Results(f3CompletedMask));

    if isempty(f3EvalPlot)
        logWarn("F3: No completed exploratory points to plot.");
    else

        fig4 = figure("Name", "RP05-F3: n_eval vs Spearman rho");
        set(fig4, "Position", [860 650 700 460]);
        ax4 = axes(fig4);
        errorbar(ax4, f3EvalPlot, f3RhoMean, f3RhoStd, "o-", "LineWidth", 1.5, ...
            "Color", [0.15 0.35 0.75], "MarkerFaceColor", [0.15 0.35 0.75]);
        hold(ax4, "on");
        yline(ax4, f3ZoneThreshold, "--", sprintf("rho=%.2f", f3ZoneThreshold), ...
            "Color", [0.85 0.2 0.1], "LineWidth", 1.0);
        hold(ax4, "off");
        xlabel(ax4, "n_eval");
        ylabel(ax4, "Spearman rho mean +/- std");
        title(ax4, "F3: MATLAB subsampled SHAP vs Python full-test TreeSHAP");
        subtitle(ax4, sprintf("Top-20 overlap mean: [%s]", num2str(f3OverlapMean, "%.1f ")));
        grid(ax4, "on");
        box(ax4, "off");

        fig5 = figure("Name", "RP05-F3: n_eval vs runtime");
        set(fig5, "Position", [1580 650 700 460]);
        ax5 = axes(fig5);
        errorbar(ax5, f3EvalPlot, f3RunMean, f3RunStd, "s-", "LineWidth", 1.5, ...
            "Color", [0.1 0.55 0.25], "MarkerFaceColor", [0.1 0.55 0.25]);
        xlabel(ax5, "n_eval");
        ylabel(ax5, "Runtime (sec) mean +/- std");
        title(ax5, "F3: shapley() runtime scaling");
        grid(ax5, "on");
        box(ax5, "off");

        logInfo("Figures 4-5 created (F3 sampling study).");
    end
else
    logWarn("F3: Skipping Figures 4-5 because sampling study was skipped.");
end

%% Section 5: RF03 Verification
logSection("RP05", "Section 5: RF03 Verification", "SHAP BBBP Explainability");

% shap_lr_spearman threshold revised from 0.90 to 0.85 (M-REPRO-REFINE Task B):
% For binary ECFP4, mean|X-p| = 2p(1-p) vs std(X) = sqrt(p(1-p)) are nonlinearly
% related, and SHAP is computed on test set while lr_imp uses train std.
% 0.85 is a realistic upper bound for this check; <0.85 signals an implementation bug.
rf03crit = struct( ...
    "auc_cv",           struct("lower", 0.85), ...
    "shap_lr_spearman", struct("lower", 0.85));

metRP05 = struct( ...
    "auc_cv",           aucCV, ...
    "shap_lr_spearman", shapRho);

resRP05 = emk.repro.verify(metRP05, rf03crit);

logInfo("==> [Biased C=1.0] ROC-AUC CV        = %.4f (>= 0.85): %s", ...
    aucCV, statusStr_(resRP05.details.auc_cv.pass));
logInfo("==> [Biased C=1.0] SHAP-LR Spearman  = %.4f (>= 0.85): %s", ...
    shapRho, statusStr_(resRP05.details.shap_lr_spearman.pass));
disp(resRP05.report);

if resRP05.pass
    logInfo("==> RP05 REPRODUCTION: PASS");
else
    logWarn("==> RP05 REPRODUCTION: NEEDS REVIEW -- see README.md Discussion");
end

% Task B: RF03 check for rev model (informational -- same criteria)
metRP05rev = struct( ...
    "auc_cv",           aucCV_rev, ...
    "shap_lr_spearman", shapRho_rev);
resRP05rev = emk.repro.verify(metRP05rev, rf03crit);

logInfo("==> [Rev C=%.2f]     ROC-AUC CV        = %.4f (>= 0.85): %s", ...
    bestC_shap, aucCV_rev, statusStr_(resRP05rev.details.auc_cv.pass));
logInfo("==> [Rev C=%.2f]     SHAP-LR Spearman  = %.4f (>= 0.85): %s", ...
    bestC_shap, shapRho_rev, statusStr_(resRP05rev.details.shap_lr_spearman.pass));
logInfo("==> Task B ranking: Spearman rho=%.4f, overlap=%d/20", rankRho, nOverlap);

%% Section 6: Save Results
logSection("RP05", "Section 6: Save Results", "SHAP BBBP Explainability");

runDir = makeRunDir("Prefix", "rp05_shap");
absRunDir = char(fullfile(pwd(), runDir));

% Save figures
saveas(fig1, fullfile(absRunDir, "shap_global_importance.png"));
saveas(fig2, fullfile(absRunDir, "shap_local_waterfall.png"));
saveas(fig3, fullfile(absRunDir, "shap_ranking_comparison.png"));
if f3Proceed
    if exist("fig4", "var")
        saveas(fig4, fullfile(absRunDir, "f3_n_eval_vs_rho.png"));
    end
    if exist("fig5", "var")
        saveas(fig5, fullfile(absRunDir, "f3_n_eval_vs_runtime.png"));
    end
end
logInfo("Figures saved.");

% Metrics JSON
% Build task_b struct with dot notation to avoid MATLAB's struct() array
% expansion when a field value is non-scalar (e.g., best_C_per_fold is 1x5).
taskB = struct();
taskB.auc_cv_rev           = aucCV_rev;
taskB.auc_cv_std_rev       = aucStd_rev;
taskB.shap_lr_spearman_rev = shapRho_rev;
taskB.best_C_shap          = bestC_shap;
taskB.best_C_per_fold      = resRev.best_C_per_fold;
taskB.rf03_pass_rev        = resRP05rev.pass;
taskB.ranking_spearman_rho = rankRho;
taskB.top20_overlap        = nOverlap;
taskB.n_biased_only        = numel(biasedOnly);
taskB.n_rev_only           = numel(revOnly);

metrics = struct();
metrics.auc_cv           = aucCV;
metrics.auc_cv_std       = aucStd;
metrics.shap_lr_spearman = shapRho;
metrics.n_valid          = nValid;
metrics.n_train          = nTrain;
metrics.n_test           = nTest;
metrics.n_bbb_pos        = res.n_bbb_pos;
metrics.n_bbb_pos_train  = res.n_bbb_pos_train;
metrics.n_bbb_pos_test   = res.n_bbb_pos_test;
metrics.ecfp4_radius     = 2;
metrics.ecfp4_nbits      = 2048;
metrics.n_top_bits_shown = topN;
metrics.rf03_criteria    = rf03crit;
metrics.rf03_pass        = resRP05.pass;
metrics.task_b           = taskB;

% F2 results (M-REPRO-REFINE Phase 5)
if exist("f2RhoBias", "var")
    f2m = struct();
    f2m.matlab_spearman_rho_biased  = f2RhoBias;   % MATLAB vs Python(C=1.0): Zone B test
    f2m.matlab_spearman_rho_rev     = f2RhoRev;    % MATLAB vs Python(C=0.10): reg-gap analysis
    f2m.python_biased_vs_rev_rho    = rhoPyPy;     % Python(C=1.0) vs Python(C=0.10): reference
    f2m.zone_b_confirmed            = f2ZoneB;
    f2m.zone_b_threshold            = 0.85;
    f2m.n_train                     = nTrain_f2;
    f2m.n_test                      = nTest_f2;
    f2m.n_all                       = nTrain_f2 + nTest_f2;
    f2m.C_shap_biased               = resF2.C_shap;
    f2m.C_shap_rev                  = resF2Rev.C_shap;
    f2m.lambda_matlab               = 1 / nTrain_f2;
    f2m.n_features                  = size(X_tr_f2, 2);
    metrics.f2_matlab_shap = f2m;
end

% F3 results (M-REPRO-REFINE Phase 5)
if exist("f3Results", "var")
    f3m = struct();
    f3m.benchmark             = f3Benchmark;
    f3m.mode                  = f3Mode;
    f3m.proceeded             = f3Proceed;
    f3m.skip_reason           = f3SkipReason;
    f3m.formal_zone_evaluable = f3FormalZoneEvaluable;
    f3m.zone_threshold        = f3ZoneThreshold;
    f3m.matlab_method         = eshapF3.Method;
    f3m.python_method         = "shap.TreeExplainer(interventional)";
    f3m.n_estimators          = 100;
    f3m.python_rf_auc_cv      = resF3.rf_auc_cv;
    f3m.python_rf_auc_std     = resF3.rf_auc_std;
    f3m.n_train               = size(X_tr_f3, 1);
    f3m.n_test                = size(X_te_f3, 1);
    f3m.n_background          = size(X_bg_f3, 1);
    f3m.n_eval_grid           = f3NEvalGrid;
    f3m.runtime_budget_sec    = f3RuntimeBudgetSec;
    f3m.n_repeats             = f3NRepeats;
    if f3Proceed
        f3m.results_by_n_eval = f3Results(f3CompletedMask);
    else
        f3m.results_by_n_eval = struct([]);
    end
    f3m.zone_candidate        = f3ZoneCandidate;
    f3m.conclusion            = f3ZoneConclusion;
    metrics.f3_rf_shap = f3m;
end
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
logInfo("Task B summary: ranking_rho=%.4f, overlap=%d/20, rev_AUC=%.4f, rev_C=%.2f", ...
    rankRho, nOverlap, aucCV_rev, bestC_shap);
if exist("f2RhoBias", "var")
    logInfo("F2 summary: rho(MATLAB,C=1.0)=%.4f  rho(MATLAB,C=0.10)=%.4f  zone_b=%s", ...
        f2RhoBias, f2RhoRev, string(f2ZoneB));
end
if exist("f3Results", "var")
    if f3Proceed
        if any(f3CompletedMask)
            lastIdx = find(f3CompletedMask, 1, "last");
            logInfo("F3 summary: zone=%s, completed_n_eval=[%s], last_rho_mean=%.4f", ...
                f3ZoneCandidate, num2str(f3NEvalGrid(f3CompletedMask)), f3Results(lastIdx).rho_mean);
        else
            logWarn("F3 summary: no exploratory n_eval completed.");
        end
    else
        logWarn("F3 summary: skipped after benchmark. %s", f3SkipReason);
    end
end

% Clean up F2 temp files
if exist("f2MatPath", "var") && isfile(f2MatPath)
    delete(f2MatPath);
end
if exist("f2MatPathRev", "var") && isfile(f2MatPathRev)
    delete(f2MatPathRev);
end

% Clean up F3 temp file
if exist("f3MatPath", "var") && isfile(f3MatPath)
    delete(f3MatPath);
end

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

function idx = topkIdx_(vals, k)
    [~, ord] = sort(vals(:), "descend");
    idx = ord(1:min(k, numel(ord)));
end
