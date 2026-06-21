% rp01_esol.m  RP01: ESOL Classical QSAR + L05 Extended Descriptors
%
% Reproduces the Delaney (2004) ESOL linear regression model (Model A) and
% extends it with L05 descriptors (Model B) to evaluate whether advanced
% descriptors improve logS prediction.
%
%   Paper:  Delaney, J.S. (2004). ESOL: Estimating Aqueous Solubility
%           Directly from Molecular Structure.
%           J. Chem. Inf. Comput. Sci. 44(3):1000-1005.
%           DOI: 10.1021/ci034243x
%
%   Model A (Delaney original 4-feature):
%     logS ~ LogP + MolWt + NumRotatableBonds + AromaticProportion
%   Model B (L05 extended 9-feature):
%     logS ~ [Model A] + TPSA + NumHDonors + NumHAcceptors + FractionCSP3 + QED
%
%   RF01: repro/rp01_esol/README.md
%   RF02: emk.setup.snapshot() + emk.setup.lockfile()
%   RF03: Both models -- RMSE_CV <= 1.20, R^2_CV >= 0.75
%
%   Prerequisites:
%     emk.setup.install() has been run once.
%     Statistics and Machine Learning Toolbox (fitlm, cvpartition).
%
%   Run: Ctrl+Enter in MATLAB with project root as CWD.

%% Section 0: Setup & Environment
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("RP01", "Section 0: Setup & Environment", "ESOL Classical QSAR + L05");
emk.setup.initPython();

molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;

snap = emk.setup.snapshot();
logInfo("RP01 setup complete.");

%% Section 1: Load ESOL Dataset
logSection("RP01", "Section 1: Load ESOL Dataset", "ESOL Classical QSAR + L05");

tbl    = emk.dataset.esol();
nTotal = height(tbl);
logInfo("Loaded %d molecules (ESOL / Delaney 2004)", nTotal);
logInfo("logS range: [%.2f, %.2f] log mol/L", min(tbl.logS), max(tbl.logS));

%% Section 2: Parse SMILES & Compute Descriptors
logSection("RP01", "Section 2: Parse SMILES & Compute Descriptors", "ESOL Classical QSAR + L05");

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

% Batch-compute all RDKit descriptors in a single IPC round-trip
allDesc = ["LogP", "MolWt", "NumRotatableBonds", "HeavyAtomCount", ...
           "TPSA", "NumHDonors", "NumHAcceptors", "FractionCSP3"];
logInfo("Computing %d RDKit descriptors for %d molecules ...", numel(allDesc), nValid);
descTbl = emk.descriptor.batchCalculate(validMols, allDesc);

% AromaticProportion (Delaney 2004): NumAromaticAtoms / HeavyAtomCount
% QED (L05): computed via single pyrun batch call for performance
logInfo("Computing AromaticProportion and QED for %d molecules ...", nValid);
[descTbl.AromaticProportion, descTbl.QED] = calcAroPropAndQED_(validMols, descTbl.HeavyAtomCount);

logInfo("Descriptor computation complete.");
logInfo("  LogP range:   [%.2f, %.2f]", min(descTbl.LogP),  max(descTbl.LogP));
logInfo("  TPSA range:   [%.2f, %.2f]", min(descTbl.TPSA),  max(descTbl.TPSA));
logInfo("  QED range:    [%.3f, %.3f]", min(descTbl.QED),   max(descTbl.QED));
logInfo("  AroProp range:[%.3f, %.3f]", min(descTbl.AromaticProportion), max(descTbl.AromaticProportion));

%% Section 3: Build Feature Matrices & Full-Data Models
logSection("RP01", "Section 3: Build Feature Matrices & Full-Data Models", "ESOL Classical QSAR + L05");

y  = yAll;

% Model A: Original Delaney (2004) 4-feature QSAR
Xa = [descTbl.LogP, descTbl.MolWt, descTbl.NumRotatableBonds, descTbl.AromaticProportion];
mdlA = fitlm(Xa, y, "VarNames", {'LogP','MolWt','RotBonds','AroProp','logS'});
rmseTrainA = sqrt(mdlA.MSE);
r2TrainA   = mdlA.Rsquared.Ordinary;
logInfo("Model A (4-feature) full: RMSE=%.4f  R^2=%.4f", rmseTrainA, r2TrainA);
disp(mdlA);

% Model B: L05-extended 9-feature QSAR
Xb = [Xa, descTbl.TPSA, descTbl.NumHDonors, descTbl.NumHAcceptors, ...
      descTbl.FractionCSP3, descTbl.QED];
mdlB = fitlm(Xb, y, "VarNames", ...
    {'LogP','MolWt','RotBonds','AroProp','TPSA','HBD','HBA','FracCSP3','QED','logS'});
rmseTrainB = sqrt(mdlB.MSE);
r2TrainB   = mdlB.Rsquared.Ordinary;
logInfo("Model B (9-feature) full: RMSE=%.4f  R^2=%.4f", rmseTrainB, r2TrainB);
disp(mdlB);

%% Section 4: 5-Fold Cross-Validation
logSection("RP01", "Section 4: 5-Fold Cross-Validation", "ESOL Classical QSAR + L05");

rng(42, "twister");
cv = cvpartition(nValid, "KFold", 5);

[rmseCVa, r2CVa, rmsesA, r2sA] = runCV_(Xa, y, cv);
[rmseCVb, r2CVb, rmsesB, r2sB] = runCV_(Xb, y, cv);

logInfo("Model A 5-fold CV: RMSE=%.4f +/- %.4f  R^2=%.4f +/- %.4f", ...
    rmseCVa, std(rmsesA), r2CVa, std(r2sA));
logInfo("Model B 5-fold CV: RMSE=%.4f +/- %.4f  R^2=%.4f +/- %.4f", ...
    rmseCVb, std(rmsesB), r2CVb, std(r2sB));
deltaRmse = rmseCVb - rmseCVa;
deltaR2   = r2CVb   - r2CVa;
logInfo("L05 effect: dRMSE=%.4f  dR^2=%.4f", deltaRmse, deltaR2);

%% Section 5: RF03 Verification
logSection("RP01", "Section 5: RF03 Verification", "ESOL Classical QSAR + L05");

rf03crit.rmse_cv = struct("upper", 1.20);
rf03crit.r2_cv   = struct("lower", 0.75);

metA.rmse_cv = rmseCVa;  metA.r2_cv = r2CVa;
metB.rmse_cv = rmseCVb;  metB.r2_cv = r2CVb;

resA = emk.repro.verify(metA, rf03crit);
resB = emk.repro.verify(metB, rf03crit);

logInfo("==> Model A (Delaney original): %s", statusStr_(resA.pass));
disp(resA.report);
logInfo("==> Model B (L05 extended):     %s", statusStr_(resB.pass));
disp(resB.report);

allPass = resA.pass && resB.pass;
if allPass
    logInfo("==> RP01 REPRODUCTION: PASS");
else
    logWarn("==> RP01 REPRODUCTION: NEEDS REVIEW -- see README.md Discussion");
end

%% Section 6: Save Results
logSection("RP01", "Section 6: Save Results", "ESOL Classical QSAR + L05");

runDir    = makeRunDir("Prefix", "rp01_esol");
absRunDir = char(fullfile(pwd(), runDir));  % absolute path required for saveas

% Predictions CSV
yPredA = predict(mdlA, Xa);
yPredB = predict(mdlB, Xb);

outTbl = tbl(validMask, ["SMILES","Name","logS"]);
outTbl.logS_predA  = yPredA;
outTbl.logS_predB  = yPredB;
outTbl.residual_A  = y - yPredA;
outTbl.residual_B  = y - yPredB;
outTbl = [outTbl, descTbl];
writetable(outTbl, fullfile(runDir, "predictions.csv"));
logInfo("Predictions saved: predictions.csv");

% Comparison scatter plot
fig = figure("Name","RP01 ESOL: Model A vs B");
subplot(1,2,1);
scatter(y, yPredA, 20, "filled", "MarkerFaceAlpha", 0.5);
hold on;
lo = floor(min([y; yPredA]));  hi = ceil(max([y; yPredA]));
plot([lo,hi],[lo,hi],"k--","LineWidth",1.2);
xlabel("Measured logS");  ylabel("Predicted logS");
title(sprintf("Model A (4-feat)  RMSE=%.3f  R^2=%.3f", rmseTrainA, r2TrainA));

subplot(1,2,2);
scatter(y, yPredB, 20, "filled", "MarkerFaceAlpha", 0.5, "MarkerFaceColor","#0072BD");
hold on;
lo = floor(min([y; yPredB]));  hi = ceil(max([y; yPredB]));
plot([lo,hi],[lo,hi],"k--","LineWidth",1.2);
xlabel("Measured logS");  ylabel("Predicted logS");
title(sprintf("Model B (9-feat)  RMSE=%.3f  R^2=%.3f", rmseTrainB, r2TrainB));

saveas(fig, fullfile(absRunDir, "predicted_vs_actual.png"));
logInfo("Figure saved: predicted_vs_actual.png");

% Metrics JSON
metrics = struct(...
    "model_a_rmse_train",    rmseTrainA, ...
    "model_a_r2_train",      r2TrainA, ...
    "model_a_rmse_cv",       rmseCVa, ...
    "model_a_rmse_cv_std",   std(rmsesA), ...
    "model_a_r2_cv",         r2CVa, ...
    "model_a_r2_cv_std",     std(r2sA), ...
    "model_b_rmse_train",    rmseTrainB, ...
    "model_b_r2_train",      r2TrainB, ...
    "model_b_rmse_cv",       rmseCVb, ...
    "model_b_rmse_cv_std",   std(rmsesB), ...
    "model_b_r2_cv",         r2CVb, ...
    "model_b_r2_cv_std",     std(r2sB), ...
    "l05_delta_rmse",        deltaRmse, ...
    "l05_delta_r2",          deltaR2, ...
    "n_molecules",           nValid, ...
    "rf03_criteria",         rf03crit, ...
    "rf03_model_a_pass",     resA.pass, ...
    "rf03_model_b_pass",     resB.pass, ...
    "rf03_pass",             allPass);
fid = fopen(fullfile(runDir, "metrics.json"), "w");
fprintf(fid, "%s\n", jsonencode(metrics, "PrettyPrint", true));
fclose(fid);
logInfo("Metrics saved: metrics.json");

% RF02 lock snapshot
snap.run_date  = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir   = runDir;
snap.rf03_pass = allPass;
emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));

logInfo("RP01 complete.  run_dir=%s", runDir);

% ===========================================================================
% Local helper functions
% ===========================================================================

function [aroProp, qedScores] = calcAroPropAndQED_(mols, heavyAtomCounts)
% Batch-compute AromaticProportion and QED in a single pyrun call.
    n = numel(mols);
    pyMolList = py.list();
    for i = 1:n
        pyMolList.append(mols{i});
    end

    [pyAros, pyQeds] = pyrun([ ...
        "from rdkit.Chem import QED as _QED;" ...
        "aros = [float(sum(a.GetIsAromatic() for a in m.GetAtoms())) for m in mols];" ...
        "qeds = [_QED.qed(m) for m in mols]"], ...
        ["aros","qeds"], mols=pyMolList);

    aros = double(py.array.array("d", pyAros));
    qeds = double(py.array.array("d", pyQeds));

    aroProp = zeros(n, 1);
    for i = 1:n
        nH = heavyAtomCounts(i);
        if ~isnan(nH) && nH > 0
            aroProp(i) = aros(i) / nH;
        end
    end
    qedScores = qeds(:);
end

function [rmseCV, r2CV, rmses, r2s] = runCV_(X, y, cv)
% 5-fold CV returning mean RMSE, mean R^2, and per-fold arrays.
    K    = cv.NumTestSets;
    rmses = zeros(K, 1);
    r2s   = zeros(K, 1);
    for fold = 1:K
        trainIdx = training(cv, fold);
        testIdx  = test(cv, fold);
        mdl      = fitlm(X(trainIdx,:), y(trainIdx));
        yPred    = predict(mdl, X(testIdx,:));
        res      = y(testIdx) - yPred;
        rmses(fold) = sqrt(mean(res.^2));
        ssTot       = sum((y(testIdx) - mean(y(testIdx))).^2);
        r2s(fold)   = 1 - sum(res.^2) / ssTot;
    end
    rmseCV = mean(rmses);
    r2CV   = mean(r2s);
end

function s = statusStr_(tf)
    if tf; s = "PASS"; else; s = "FAIL"; end
end
