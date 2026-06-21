% rp00_esol_pilot.m  RP00: ESOL Aqueous Solubility -- Pilot Reproduction
%
% Reproduces the linear regression model from Delaney (2004).
%
%   Paper:   Delaney, J.S. (2004). ESOL: Estimating Aqueous Solubility
%            Directly from Molecular Structure.
%            J. Chem. Inf. Comput. Sci. 44(3):1000-1005.
%            DOI: 10.1021/ci034243x
%
%   Model:   logS ~ LogP + MolWt + NumRotatableBonds + AromaticProportion
%            Linear regression, 5-fold cross-validation.
%
%   Deliverables (M-REPRO-PILOT / updated in M-REPRO-FOUND):
%     RF01 -- Standard reproduction template: see README.md and repro/TEMPLATE.md
%     RF02 -- Version lock: emk.setup.snapshot() + emk.setup.lockfile()
%     RF03 -- Verification: emk.repro.verify(), RMSE <= 1.20, R^2 >= 0.75
%
%   Prerequisites:
%     emk.setup.install() has been run once.
%     Statistics and Machine Learning Toolbox (fitlm, cvpartition).
%
%   Run: Ctrl+Enter in MATLAB with this file open (project root as CWD).

%% Section 0: Setup & Environment
% Bootstrap: navigate to project root regardless of launch directory.
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("RP00", "Section 0: Setup & Environment", "Repro Pilot");
emk.setup.initPython();

% Warm up Python/RDKit process
molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;

% Capture environment info for RF02 version lock
snap = emk.setup.snapshot();
logInfo("RP00 setup complete.");

%% Section 1: Load ESOL Dataset
logSection("RP00", "Section 1: Load ESOL Dataset", "Repro Pilot");

tbl    = emk.dataset.esol();
nTotal = height(tbl);
logInfo("Loaded %d molecules (ESOL / Delaney 2004)", nTotal);
logInfo("logS range: [%.2f, %.2f] log mol/L", min(tbl.logS), max(tbl.logS));

%% Section 2: Parse SMILES & Compute Descriptors
logSection("RP00", "Section 2: Parse SMILES & Compute Descriptors", "Repro Pilot");

% Parse all SMILES to RDKit mol objects
mols      = cell(1, nTotal);
validMask = false(1, nTotal);

for i = 1:nTotal
    try
        mols{i}      = emk.mol.fromSmiles(tbl.SMILES(i));
        validMask(i) = true;
    catch
        logWarn("SMILES parse failed [%d]: %s", i, tbl.SMILES(i));
    end
end

nValid    = sum(validMask);
validMols = mols(validMask);
yAll      = tbl.logS(validMask);
logInfo("Parsed %d / %d SMILES successfully", nValid, nTotal);

% Batch-compute Delaney (2004) base descriptors via RDKit
% LogP  = RDKit Crippen MolLogP (proxy for Delaney's clogP)
% MolWt = average molecular weight
% NumRotatableBonds = strict SMARTS definition
% HeavyAtomCount    = denominator for AromaticProportion
descNames = ["LogP", "MolWt", "NumRotatableBonds", "HeavyAtomCount"];
logInfo("Computing %d descriptors for %d molecules ...", numel(descNames), nValid);
descTbl = emk.descriptor.batchCalculate(validMols, descNames);

% AromaticProportion = NumAromaticAtoms / HeavyAtomCount  [Delaney 2004, Eq. 1]
descTbl.AromaticProportion = calcAromaticProportion_(validMols, descTbl.HeavyAtomCount);
logInfo("Descriptor computation complete.");
logInfo("  LogP range: [%.2f, %.2f]", min(descTbl.LogP), max(descTbl.LogP));
logInfo("  AroProp range: [%.3f, %.3f]", ...
    min(descTbl.AromaticProportion), max(descTbl.AromaticProportion));

%% Section 3: Linear Regression Model (Delaney 2004)
logSection("RP00", "Section 3: Linear Regression Model", "Repro Pilot");

% Feature matrix X: [LogP | MolWt | RotBonds | AromaticProportion]
X = [descTbl.LogP, descTbl.MolWt, descTbl.NumRotatableBonds, descTbl.AromaticProportion];
y = yAll;

mdlFull   = fitlm(X, y, "VarNames", {'LogP', 'MolWt', 'RotBonds', 'AroProp', 'logS'});
rmseTrain = sqrt(mdlFull.MSE);
r2Train   = mdlFull.Rsquared.Ordinary;
logInfo("Full model (all %d mols): RMSE=%.4f, R^2=%.4f", nValid, rmseTrain, r2Train);
disp(mdlFull);

%% Section 4: 5-Fold Cross-Validation
logSection("RP00", "Section 4: 5-Fold Cross-Validation", "Repro Pilot");

rng(42, "twister");  % fixed seed for reproducibility
cv    = cvpartition(nValid, "KFold", 5);
rmses = zeros(5, 1);
r2s   = zeros(5, 1);

for fold = 1:5
    trainIdx = training(cv, fold);
    testIdx  = test(cv, fold);

    mdlFold = fitlm(X(trainIdx, :), y(trainIdx));
    yPred   = predict(mdlFold, X(testIdx, :));

    res         = y(testIdx) - yPred;
    rmses(fold) = sqrt(mean(res .^ 2));
    ssTot       = sum((y(testIdx) - mean(y(testIdx))) .^ 2);
    r2s(fold)   = 1 - sum(res .^ 2) / ssTot;

    logInfo("  Fold %d: RMSE=%.4f, R^2=%.4f", fold, rmses(fold), r2s(fold));
end

rmseCV = mean(rmses);
r2CV   = mean(r2s);
logInfo("5-fold CV: RMSE=%.4f +/- %.4f", rmseCV, std(rmses));
logInfo("5-fold CV: R^2 =%.4f +/- %.4f", r2CV,   std(r2s));

%% Section 5: Verification (RF03)
logSection("RP00", "Section 5: Verification", "Repro Pilot");

% RF03 thresholds (confirmed in M-REPRO-FOUND based on RP00 pilot results).
% Reference: Delaney (2004) Table 2 -- training RMSE = 0.996.
% Tolerance +0.20 accounts for: (1) RDKit Crippen-Wildman LogP vs. Delaney clogP;
% (2) dataset version difference (MoleculeNet 1128 vs. paper 1144 mols).
% R^2 >= 0.75 (not 0.80): RDKit LogP-based reproduction upper limit is ~0.76.
DELANEY_RMSE_REF = 0.996;

rf03crit.rmse_cv = struct("upper", 1.20);
rf03crit.r2_cv   = struct("lower", 0.75);

rf03met.rmse_cv  = rmseCV;
rf03met.r2_cv    = r2CV;

rf03result = emk.repro.verify(rf03met, rf03crit);
disp(rf03result.report);

allPass = rf03result.pass;
logInfo("Delaney ref RMSE: %.3f (training set, Table 2)", DELANEY_RMSE_REF);
if allPass
    logInfo("==> RP00 REPRODUCTION: PASS");
else
    logWarn("==> RP00 REPRODUCTION: NEEDS REVIEW -- see README.md Discussion");
end

%% Section 6: Save Results
logSection("RP00", "Section 6: Save Results", "Repro Pilot");

runDir = makeRunDir("Prefix", "rp00_esol_pilot");

% Predictions CSV
yPredFull = predict(mdlFull, X);
outTbl = tbl(validMask, ["SMILES", "Name", "logS"]);
outTbl.logS_pred = yPredFull;
outTbl.residual  = y - yPredFull;
outTbl = [outTbl, descTbl];
writetable(outTbl, fullfile(runDir, "predictions.csv"));
logInfo("Predictions saved: %s/predictions.csv", runDir);

% Scatter plot: measured vs. predicted
fig = figure("Name", "RP00 ESOL: Measured vs Predicted");
scatter(y, yPredFull, 25, "filled", "MarkerFaceAlpha", 0.5);
hold on;
lo = floor(min([y; yPredFull]));
hi = ceil(max([y; yPredFull]));
plot([lo, hi], [lo, hi], "k--", "LineWidth", 1.2);
xlabel("Measured logS (log mol/L)");
ylabel("Predicted logS (log mol/L)");
title(sprintf("RP00 ESOL Full Model  RMSE=%.3f  R^2=%.3f", rmseTrain, r2Train));
saveas(fig, fullfile(runDir, "predicted_vs_actual.png"));
logInfo("Figure saved: %s/predicted_vs_actual.png", runDir);

% Metrics JSON
metrics = struct( ...
    "rmse_train",       rmseTrain, ...
    "r2_train",         r2Train, ...
    "rmse_cv",          rmseCV, ...
    "rmse_cv_std",      std(rmses), ...
    "r2_cv",            r2CV, ...
    "r2_cv_std",        std(r2s), ...
    "n_molecules",      nValid, ...
    "delaney_rmse_ref", DELANEY_RMSE_REF, ...
    "rf03_criteria",    rf03crit, ...
    "rf03_pass",        allPass);
fid = fopen(fullfile(runDir, "metrics.json"), "w");
fprintf(fid, "%s\n", jsonencode(metrics, "PrettyPrint", true));
fclose(fid);
logInfo("Metrics saved: %s/metrics.json", runDir);

% Lock snapshot (RF02 -- actual versions for this run via emk.setup.lockfile)
snap.run_date   = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir    = runDir;
snap.rf03_pass  = allPass;
emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));

logInfo("RP00 ESOL pilot complete. Output: %s", runDir);

% ==========================================================================
% Local helper functions
% ==========================================================================

function prop = calcAromaticProportion_(mols, heavyAtomCounts)
% Compute Delaney (2004) AromaticProportion = NumAromaticAtoms / HeavyAtomCount.
% Uses a single pyrun batch call to count aromatic atoms across all molecules.
    n = numel(mols);
    prop = zeros(n, 1);

    % Build py.list of all mol objects (one IPC round-trip)
    pyMolList = py.list();
    for i = 1:n
        pyMolList.append(mols{i});
    end

    % Batch-compute aromatic atom counts in Python
    pyAros = pyrun( ...
        "aros = [float(sum(a.GetIsAromatic() for a in mol.GetAtoms())) for mol in mols]", ...
        "aros", mols=pyMolList);
    aros = double(py.array.array("d", pyAros));

    for i = 1:n
        nHeavy = heavyAtomCounts(i);
        if ~isnan(nHeavy) && nHeavy > 0
            prop(i) = aros(i) / nHeavy;
        end
    end
end

