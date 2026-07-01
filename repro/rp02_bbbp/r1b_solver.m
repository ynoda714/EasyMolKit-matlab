% r1b_solver.m  R1-B: fitclinear solver/convergence comparison (RP02)
%
% Investigates whether the large solver_gap (+0.1626) in R1 (rp02_bbbp.m
% Section 2b) is due to insufficient convergence or a fundamental limitation.
%
% All experiments use Lambda=1/C_opt and the same fold splits as sklearn.
% Convergence parameter names (from fitclinear help):
%   SGD/ASGD solvers  -> PassLimit      (default 1  for sgd, 10 for dual)
%   dual SGD          -> PassLimit      (default 10)
%   BFGS/LBFGS/SpaRSA -> IterationLimit (default 1000)
%
% Configurations:
%   E1: default solver,  no extra args         (R1 baseline replicate)
%   E2: default solver,  PassLimit=10          (10x more SGD passes)
%   E3: default solver,  PassLimit=50          (50x more SGD passes)
%   E4: Solver='sparsa', IterationLimit=1000
%   E5: Solver='sparsa', IterationLimit=5000
%   E6: Solver='lbfgs',  IterationLimit=1000
%   E7: Solver='lbfgs',  IterationLimit=5000
%   E8: Solver='bfgs',   IterationLimit=1000
%
% Output: result/runs/r1b_solver_<ts>/r1b_results.json
%
% Run: Ctrl+Enter in MATLAB with project root as CWD.

%% Section 0: Setup
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("R1B", "Section 0: Setup", "fitclinear Solver Compare");
emk.setup.initPython();
molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;

%% Section 1: Load BBBP data and compute ECFP4
logSection("R1B", "Section 1: Load BBBP + ECFP4", "fitclinear Solver Compare");

tbl     = emk.dataset.bbbp();
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

%% Section 2: Python sklearn CV (fold indices and best_C per fold)
logSection("R1B", "Section 2: Python sklearn CV (fold splits + best_C)", "fitclinear Solver Compare");

coreScript = fullfile(thisDir, "rp02_sklearn_core.py");
logInfo("Running sklearn nested CV for fold splits ...");
resultJson = pyrun([ ...
    "exec(open(script_path).read()); " ...
    "result_json = run_rp02_sklearn(csv_path)"], ...
    "result_json", ...
    script_path=coreScript, csv_path=csvPath);
res = jsondecode(char(string(resultJson)));

if isfield(res, "success") && ~res.success
    error("emk:r1b:pythonError", "sklearn failed: %s", res.error);
end

sklearnAUC = res.auc_mean;
logInfo("sklearn AUC = %.4f +/- %.4f", sklearnAUC, res.auc_std);
logInfo("best_C per fold: %s", mat2str([res.fold_results.best_C], 3));
logInfo("Lambda per fold: %s", mat2str(1 ./ [res.fold_results.best_C], 3));

if size(X,1) ~= res.n_valid
    error("emk:r1b:sizeMismatch", ...
        "MATLAB n_valid=%d vs Python n_valid=%d", size(X,1), res.n_valid);
end

%% Section 3: fitclinear configurations
logSection("R1B", "Section 3: fitclinear Solver/Convergence Comparison", "fitclinear Solver Compare");

% Experiment definitions (semicolon separates rows).
% PassLimit applies to SGD/ASGD/dual; IterationLimit to BFGS/LBFGS/SpaRSA.
cfgs = { ...
    "E1_default",          "",       {}; ...
    "E2_default_pass10",   "",       {"PassLimit", 10}; ...
    "E3_default_pass50",   "",       {"PassLimit", 50}; ...
    "E4_sparsa_iter1000",  "sparsa", {"IterationLimit", 1000}; ...
    "E5_sparsa_iter5000",  "sparsa", {"IterationLimit", 5000}; ...
    "E6_lbfgs_iter1000",   "lbfgs",  {"IterationLimit", 1000}; ...
    "E7_lbfgs_iter5000",   "lbfgs",  {"IterationLimit", 5000}; ...
    "E8_bfgs_iter1000",    "bfgs",   {"IterationLimit", 1000}; ...
};
nCfgs  = size(cfgs, 1);
nFolds = numel(res.fold_results);

expInfos = {};
for ci = 1:nCfgs
    label    = cfgs{ci, 1};
    solver   = cfgs{ci, 2};
    extra    = cfgs{ci, 3};
    aucs     = zeros(nFolds, 1);
    convFlgs = zeros(nFolds, 1);  % FitInfo.Converged per fold
    errored  = false;

    for k = 1:nFolds
        fr        = res.fold_results(k);
        lambdaOpt = 1.0 / fr.best_C;
        trIdx     = res.outer_fold_indices.(sprintf("fold%d_train", k)) + 1;
        teIdx     = res.outer_fold_indices.(sprintf("fold%d_test",  k)) + 1;

        baseArgs = {"Learner", "logistic", "Regularization", "ridge", ...
                    "Lambda", lambdaOpt};
        if ~isempty(solver)
            baseArgs = [baseArgs, {"Solver", solver}];
        end
        fitArgs = [baseArgs, extra];

        try
            [mdl, fi] = fitclinear(X(trIdx,:), y(trIdx), fitArgs{:});
            pIdx = find(mdl.ClassNames == 1, 1);
            [~, sc] = predict(mdl, X(teIdx,:));
            [~,~,~, aucK] = perfcurve(y(teIdx), sc(:, pIdx), 1);
            aucs(k)     = aucK;
            % FitInfo.Converged is a scalar bool for batch solvers; may be
            % empty for SGD (which uses PassLimit, not convergence criterion).
            if isfield(fi, "Converged") && ~isempty(fi.Converged)
                convFlgs(k) = double(fi.Converged);
            else
                convFlgs(k) = NaN;
            end
        catch ME
            logWarn("  %s fold %d FAILED: %s", label, k, ME.message);
            aucs(k)     = NaN;
            convFlgs(k) = NaN;
            errored     = true;
        end
    end

    aucMean = mean(aucs(~isnan(aucs)));
    aucStd  = std(aucs(~isnan(aucs)));
    delta   = sklearnAUC - aucMean;

    solverDisp = solver;
    if isempty(solverDisp); solverDisp = "default"; end
    convStr = "";
    if ~all(isnan(convFlgs))
        convStr = sprintf("  conv=%s", mat2str(convFlgs', 0));
    end
    logInfo("  %s: AUC=%.4f +/-%.4f  gap=%+.4f  (solver=%s)%s", ...
        label, aucMean, aucStd, delta, solverDisp, convStr);

    expInfos{end+1} = struct( ...
        "label",        label, ...
        "solver",       solverDisp, ...
        "extra_args",   {extra}, ...
        "auc_mean",     aucMean, ...
        "auc_std",      aucStd, ...
        "auc_per_fold", aucs', ...
        "solver_gap",   delta, ...
        "converged",    convFlgs', ...
        "errored",      errored);
end

%% Section 4: Summary report
logSection("R1B", "Section 4: Summary", "fitclinear Solver Compare");

R1_DEFAULT_AUC      = 0.7517;
FITCLINEAR_HIST_AUC = 0.8826;

logInfo("======= R1-B SOLVER COMPARISON SUMMARY =======");
logInfo("  sklearn lbfgs (C=C_opt)                : AUC = %.4f  [target]", sklearnAUC);
logInfo("  fitclinear (Lambda=1/n, historical)    : AUC = %.4f  [original RP02]", FITCLINEAR_HIST_AUC);
logInfo("  fitclinear (Lambda=1/C_opt, R1 default): AUC = %.4f  [R1 baseline]", R1_DEFAULT_AUC);
logInfo("");
logInfo("  %-30s  AUC        gap (vs sklearn)", "Configuration");
logInfo("  %-30s  ---------  ----------------", "-------------");

bestAUC = -Inf; bestIdx = 0;
for ci = 1:numel(expInfos)
    e = expInfos{ci};
    logInfo("  %-30s  %.4f     %+.4f", e.label, e.auc_mean, e.solver_gap);
    if e.auc_mean > bestAUC && ~e.errored
        bestAUC = e.auc_mean;
        bestIdx = ci;
    end
end

logInfo("");
if bestIdx > 0
    bestE = expInfos{bestIdx};
    logInfo("  Best fitclinear config : %s", bestE.label);
    logInfo("  Best AUC               : %.4f  (solver_gap=%+.4f)", bestAUC, bestE.solver_gap);
    if bestE.solver_gap > 0.05
        logInfo("  >> CONCLUSION: gap remains large (>0.05) -- fundamental solver limitation.");
    elseif bestE.solver_gap > 0.01
        logInfo("  >> CONCLUSION: gap narrows but persists -- partial convergence improvement.");
    else
        logInfo("  >> CONCLUSION: gap closed -- convergence settings were the main issue.");
    end
end

%% Section 5: Save results
logSection("R1B", "Section 5: Save Results", "fitclinear Solver Compare");

runDir = makeRunDir("Prefix", "r1b_solver");
if startsWith(runDir, '/') || startsWith(runDir, '\\') || ...
        (numel(runDir) >= 2 && runDir(2) == ':')
    absRunDir = char(runDir);
else
    absRunDir = char(fullfile(pwd(), runDir));
end

output = struct( ...
    "sklearn_auc",          sklearnAUC, ...
    "r1_default_auc_ref",   R1_DEFAULT_AUC, ...
    "fitclinear_hist_auc",  FITCLINEAR_HIST_AUC);
if bestIdx > 0
    output.best_config = expInfos{bestIdx}.label;
    output.best_auc    = bestAUC;
    output.remaining_solver_gap = sklearnAUC - bestAUC;
end
output.experiments = expInfos;

fid = fopen(fullfile(absRunDir, "r1b_results.json"), "w");
if fid == -1
    error("emk:r1b:fopenFailed", "Cannot open for writing: %s", absRunDir);
end
fprintf(fid, "%s\n", jsonencode(output, "PrettyPrint", true));
fclose(fid);
logInfo("Results saved: r1b_results.json  (run_dir=%s)", runDir);

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
