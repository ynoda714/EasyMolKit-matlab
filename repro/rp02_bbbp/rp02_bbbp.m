% rp02_bbbp.m  RP02: MoleculeNet BBBP Classification Baseline (ECFP4 + Logistic Regression)
%
% Reproduces the MoleculeNet (Wu et al. 2018) classification baseline for
% Blood-Brain Barrier Permeability (BBBP) using Morgan fingerprints (ECFP4)
% with logistic regression.
%
%   Paper:  Wu, Z. et al. (2018). MoleculeNet: A Benchmark for Molecular
%           Machine Learning. Chem. Sci. 9:513-530.
%           DOI: 10.1039/C7SC02664A
%
%   Model:  Logistic Regression (L2 ridge) with Morgan ECFP4
%             Radius=2, NBits=2048 (Wu et al. Table 4 "Logreg with Circular FP")
%
%   Task:   Binary classification -- BBB permeability (BBB+ = 1, BBB- = 0)
%
%   Split:  5-fold stratified random CV (NOT scaffold split as in Wu et al.;
%           see README.md Discussion for impact on expected AUC).
%
%   RF01: repro/rp02_bbbp/README.md
%   RF02: emk.setup.snapshot() + emk.setup.lockfile()
%   RF03: ROC-AUC CV (5-fold, random split, stratified) >= 0.85
%
%   Prerequisites:
%     emk.setup.install() has been run once.
%     Statistics and Machine Learning Toolbox (fitclinear, perfcurve, cvpartition).
%
%   Run: Ctrl+Enter in MATLAB with project root as CWD.

%% Section 0: Setup & Environment
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("RP02", "Section 0: Setup & Environment", "BBBP Classification Baseline");
emk.setup.initPython();

molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;

snap = emk.setup.snapshot();
logInfo("RP02 setup complete.");

%% Section 1: Load BBBP Dataset
logSection("RP02", "Section 1: Load BBBP Dataset", "BBBP Classification Baseline");

tbl    = emk.dataset.bbbp();
nTotal = height(tbl);
nPos   = sum(tbl.BBB);
nNeg   = sum(~tbl.BBB);
logInfo("Loaded %d molecules (BBB+: %d [%.1f%%], BBB-: %d [%.1f%%])", ...
    nTotal, nPos, 100*nPos/nTotal, nNeg, 100*nNeg/nTotal);

%% Section 2: Parse SMILES & Compute Morgan Fingerprints (ECFP4)
logSection("RP02", "Section 2: Parse SMILES & Compute Morgan Fingerprints", ...
    "BBBP Classification Baseline");

nBits  = 2048;
radius = 2;  % ECFP4 = Morgan radius 2 (Wu et al. standard)

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
yAll      = tbl.BBB(validMask);   % logical: true=BBB+, false=BBB-
y         = double(yAll);          % 1=BBB+, 0=BBB- (required by fitclinear/perfcurve)
logInfo("Parsed %d / %d SMILES successfully", nValid, nTotal);

% Batch-compute Morgan ECFP4 fingerprints in a single IPC round-trip.
% Returns logical matrix n x nBits via concatenated bit-string approach.
logInfo("Computing Morgan ECFP4 (radius=%d, nBits=%d) for %d molecules ...", ...
    radius, nBits, nValid);
X = batchMorganFP_(validMols, radius, nBits);
density = 100 * sum(X(:)) / numel(X);
logInfo("Fingerprint matrix: %d x %d  (mean bit density: %.2f%%)", ...
    size(X,1), size(X,2), density);

%% Section 3: Build Full-Data Logistic Regression Model
logSection("RP02", "Section 3: Build Full-Data Logistic Regression Model", ...
    "BBBP Classification Baseline");

% fitclinear: efficient logistic regression for high-dimensional sparse data.
% L2 (ridge) regularization with default Lambda=1/n.
% This corresponds to Wu et al. 'Logreg with Circular FP' baseline.
mdlFull = fitclinear(X, y, "Learner", "logistic", "Regularization", "ridge");

posIdx = find(mdlFull.ClassNames == 1, 1);   % index of BBB+ class in ClassNames
[labelsFull, scoresFull] = predict(mdlFull, X);
[~, ~, ~, aucFull] = perfcurve(y, scoresFull(:, posIdx), 1);
accFull  = mean(labelsFull == y);
baFull   = balancedAcc_(y, labelsFull);

logInfo("Full model (train): AUC=%.4f  Acc=%.4f  BalAcc=%.4f", ...
    aucFull, accFull, baFull);

%% Section 4: 5-Fold Stratified Cross-Validation
logSection("RP02", "Section 4: 5-Fold Stratified Cross-Validation", ...
    "BBBP Classification Baseline");

rng(42, "twister");
cv    = cvpartition(y, "KFold", 5);   % stratified: each fold preserves class ratio
nFold = cv.NumTestSets;

aucs  = zeros(nFold, 1);
accs  = zeros(nFold, 1);
bAccs = zeros(nFold, 1);

for fold = 1:nFold
    trainIdx = training(cv, fold);
    testIdx  = test(cv, fold);

    mdlFold  = fitclinear(X(trainIdx,:), y(trainIdx), ...
        "Learner", "logistic", "Regularization", "ridge");

    pIdx = find(mdlFold.ClassNames == 1, 1);
    [labFold, scFold] = predict(mdlFold, X(testIdx,:));

    [~, ~, ~, aucs(fold)] = perfcurve(y(testIdx), scFold(:, pIdx), 1);
    accs(fold)  = mean(labFold == y(testIdx));
    bAccs(fold) = balancedAcc_(y(testIdx), labFold);
end

aucCV  = mean(aucs);
accCV  = mean(accs);
bAccCV = mean(bAccs);

logInfo("5-fold CV: AUC=%.4f +/- %.4f  Acc=%.4f +/- %.4f  BalAcc=%.4f +/- %.4f", ...
    aucCV, std(aucs), accCV, std(accs), bAccCV, std(bAccs));

%% Section 5: RF03 Verification
logSection("RP02", "Section 5: RF03 Verification", "BBBP Classification Baseline");

rf03crit = struct("auc_cv", struct("lower", 0.85));   % explicit init avoids workspace contamination

metRP02 = struct("auc_cv", aucCV);
resRP02 = emk.repro.verify(metRP02, rf03crit);

logInfo("==> ROC-AUC CV = %.4f (criterion: >= %.2f): %s", ...
    aucCV, rf03crit.auc_cv.lower, statusStr_(resRP02.pass));
disp(resRP02.report);

if resRP02.pass
    logInfo("==> RP02 REPRODUCTION: PASS");
else
    logWarn("==> RP02 REPRODUCTION: NEEDS REVIEW -- see README.md Discussion");
end

%% Section 6: Save Results
logSection("RP02", "Section 6: Save Results", "BBBP Classification Baseline");

runDir    = makeRunDir("Prefix", "rp02_bbbp");
absRunDir = char(fullfile(pwd(), runDir));

% Predictions CSV (full model)
outTbl = tbl(validMask, ["SMILES","Name","BBB"]);
outTbl.BBB_pred  = logical(labelsFull);
outTbl.score_pos = scoresFull(:, posIdx);
outTbl.correct   = outTbl.BBB == outTbl.BBB_pred;
writetable(outTbl, fullfile(runDir, "predictions.csv"));
logInfo("Predictions saved: predictions.csv");

% ROC curve (full model, for reference)
[rocX, rocY] = perfcurve(y, scoresFull(:, posIdx), 1);
fig = figure("Name", "RP02 BBBP: ROC Curve");
plot(rocX, rocY, "b-", "LineWidth", 1.8);
hold on;
plot([0 1], [0 1], "k--", "LineWidth", 1.0);
xlabel("False Positive Rate");
ylabel("True Positive Rate");
title(sprintf("BBBP LR+ECFP4  AUC (train)=%.3f  CV=%.3f (5-fold)", aucFull, aucCV));
legend(sprintf("LR+ECFP4 train (AUC=%.3f)", aucFull), "Random", "Location", "SouthEast");
saveas(fig, fullfile(absRunDir, "roc_curve.png"));
logInfo("Figure saved: roc_curve.png");

% Metrics JSON
metrics = struct( ...
    "model_auc_train",    aucFull, ...
    "model_acc_train",    accFull, ...
    "model_balacc_train", baFull, ...
    "auc_cv",             aucCV, ...
    "auc_cv_std",         std(aucs), ...
    "acc_cv",             accCV, ...
    "acc_cv_std",         std(accs), ...
    "balacc_cv",          bAccCV, ...
    "balacc_cv_std",      std(bAccs), ...
    "n_molecules",        nValid, ...
    "n_bbb_pos",          sum(y == 1), ...
    "n_bbb_neg",          sum(y == 0), ...
    "ecfp4_radius",       radius, ...
    "ecfp4_nbits",        nBits, ...
    "rf03_criteria",      rf03crit, ...
    "rf03_pass",          resRP02.pass);
fid = fopen(fullfile(runDir, "metrics.json"), "w");
fprintf(fid, "%s\n", jsonencode(metrics, "PrettyPrint", true));
fclose(fid);
logInfo("Metrics saved: metrics.json");

% RF02 lock snapshot
snap.run_date  = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir   = runDir;
snap.rf03_pass = resRP02.pass;
emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));

logInfo("RP02 complete.  run_dir=%s", runDir);

% ===========================================================================
% Local helper functions
% ===========================================================================

function X = batchMorganFP_(mols, radius, nBits)
% Batch-compute Morgan ECFP4 fingerprints in a single IPC round-trip.
% All bit strings are concatenated in Python and reshaped in MATLAB.
    n = numel(mols);
    pyMolList = py.list();
    for i = 1:n
        pyMolList.append(mols{i});
    end

    % Join all n bit strings into one string of length n*nBits (1 IPC call)
    pyAllBits = pyrun([ ...
        "from rdkit.Chem import rdFingerprintGenerator as _rfg;" ...
        "gen = _rfg.GetMorganGenerator(radius=int(fp_r), fpSize=int(fp_nb));" ...
        "fps = ''.join(gen.GetFingerprint(m).ToBitString() for m in mols)"], ...
        "fps", mols=pyMolList, fp_r=int32(radius), fp_nb=int32(nBits));

    allBits = char(string(pyAllBits)) == '1';    % logical row vector: 1 x (n*nBits)
    X = double(reshape(allBits, nBits, n)');     % n x nBits double (fitclinear requires numeric)
end

function bacc = balancedAcc_(yTrue, yPred)
% Balanced accuracy = (sensitivity + specificity) / 2.
    tp   = sum(yTrue == 1 & yPred == 1);
    tn   = sum(yTrue == 0 & yPred == 0);
    fp   = sum(yTrue == 0 & yPred == 1);
    fn   = sum(yTrue == 1 & yPred == 0);
    sens = tp / max(tp + fn, 1);
    spec = tn / max(tn + fp, 1);
    bacc = (sens + spec) / 2;
end

function s = statusStr_(tf)
    if tf; s = "PASS"; else; s = "FAIL"; end
end
