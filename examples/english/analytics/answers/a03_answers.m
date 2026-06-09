%[text] # A03 Answers: QSAR Regression (LogP Prediction)
%[text] This is a reference answer for the "Try it yourself" exercise in a03_qsar_regression.m.
%[text] First, run a03_qsar_regression.m to prepare the workspace variables.
%[text] Note that this file can also be executed standalone, and each section can be run independently.
addpath(genpath("src"));
emk.setup.initPython();
logInfo("A03 Answers: Setup Complete");

%[text] Reconstruct workspace (for standalone execution)
DATA_FILE  = "data/list/fda_drugs.csv";
FEAT_NAMES = ["MolWt","TPSA","NumHDonors","NumHAcceptors", ...
              "NumRotatableBonds","RingCount","FractionCSP3","HeavyAtomCount"];
N_FEATS    = numel(FEAT_NAMES);

rawTbl  = readtable(DATA_FILE, TextType="string");
nRaw    = height(rawTbl);
mols    = cell(1, nRaw);
valid   = false(1, nRaw);
for k = 1:nRaw
    try
        mols{k} = emk.mol.fromSmiles(rawTbl.SMILES(k));
        valid(k) = true;
    catch
    end
end
validIdx = find(valid);
molNames = rawTbl.Name(validIdx);
nMols    = numel(validIdx);

X = nan(nMols, N_FEATS);
y = nan(nMols, 1);
for j = 1:nMols
    s     = emk.descriptor.calculate(mols{validIdx(j)}, FEAT_NAMES);
    alogp = emk.descriptor.calculate(mols{validIdx(j)}, "LogP");
    for d = 1:N_FEATS; X(j, d) = s.(FEAT_NAMES(d)); end
    y(j) = alogp.LogP;
    logProgress(j, nMols, "Calculating descriptors");
end

Xmean = mean(X, 1);
Xstd  = std(X, 0, 1);
Xstd(Xstd < 1e-12) = 1;
Xs = (X - Xmean) ./ Xstd;
%%
%[text] ## Let's Try 1: Correlation of Each Descriptor with ALogP; Multicollinearity Check

corrVec = corr(X, y);
logInfo("Pearson correlation with ALogP:");
for d = 1:N_FEATS
    logInfo("  %-22s  r = %+.3f", FEAT_NAMES(d), corrVec(d));
end
[~, sortCorr] = sort(abs(corrVec), "descend");
logInfo("Strongest correlated descriptor: %s (r=%+.3f)", ...
    FEAT_NAMES(sortCorr(1)), corrVec(sortCorr(1)));

%[text] Check for pairwise collinearity
C = corr(X);
C_off = C - diag(diag(C));
[maxC, linIdx] = max(abs(C_off(:)));
[r, c] = ind2sub(size(C), linIdx);
logInfo("Most collinear pair: %s and %s  (r = %.3f)", ...
    FEAT_NAMES(r), FEAT_NAMES(c), C(r, c));

%[text] Answer: TPSA (negative correlation r ~-0.6) and NumHAcceptors show the strongest correlation with ALogP.
%[text]    It is chemically natural that polar surface area reduces lipophilicity.
%[text]    MolWt and HeavyAtomCount have strong collinearity (r > 0.95).
%[text]    High collinearity increases the variance of coefficients in ordinary least squares (OLS), making individual coefficient estimates unstable
%[text]    (sign reversal can also occur).
%%
%[text] ## Let's Try 2: Identify High ALogP Outliers; Check Distribution Shape

[sorted_y, ord] = sort(y, "descend");
logInfo("Top 5 Most Lipophilic Molecules:");
for k = 1:5
    logInfo("  %d. %-30s  ALogP = %.2f", k, molNames(ord(k)), sorted_y(k));
end

logInfo("ALogP: Mean=%.2f  Median=%.2f  Skewness=%.2f", ...
    mean(y), median(y), skewness(y));
logInfo("Number of Molecules with ALogP > 6: %d", sum(y > 6));

%[text] Answer: Molecules with ALogP > 6 are highly lipophilic compounds (steroids and lipid-like compounds are typical examples),
%[text]    and because OLS treats all residuals equally as penalties, these outliers can pull the decision hyperplane
%[text]    towards them, potentially reducing the overall fit quality for the majority of molecules.
%[text]    Practical approaches include excluding or down-weighting outliers, or using robust regression (robustfit()).
%[text]    If the skewness is close to 0, the LogP distribution is nearly symmetric and can be modeled without logarithmic transformation.
%%
%[text] ## Let's Try 3: Linear Regression — Large Coefficients; Insignificant Features

ftTbl = array2table(Xs, VariableNames=cellstr(FEAT_NAMES));
ftTbl.ALogP = y;
lmModel = fitlm(ftTbl, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));

coeffTbl = lmModel.Coefficients;
logInfo("Linear model coefficients (standardized features):");
disp(coeffTbl);

%[text] Extract statistically insignificant features with p > 0.05
pVals   = coeffTbl.pValue(2:end);   % Skip the intercept
names   = string(coeffTbl.Properties.RowNames(2:end));
insig   = names(pVals > 0.05);
logInfo("Statistically insignificant (p>0.05): %s", strjoin(insig, ", "));

%[text] Remove HeavyAtomCount and refit
dropIdx  = strcmp(cellstr(FEAT_NAMES), "HeavyAtomCount");
featRed  = FEAT_NAMES(~dropIdx);
trTblRed = array2table(Xs(:, ~dropIdx), VariableNames=cellstr(featRed));
trTblRed.ALogP = y;
lm2 = fitlm(trTblRed, "ALogP ~ " + strjoin(featRed, " + "));
logInfo("Full model RMSE: %.3f   Reduced model RMSE: %.3f", ...
    lmModel.RMSE, lm2.RMSE);

%[text] Check molecules with large residuals
residuals = lmModel.Residuals.Raw;
[~, bigResOrd] = sort(abs(residuals), "descend");
logInfo("Molecules with largest |residual|:");
for k = 1:5
    logInfo("  %-30s  residual = %+.2f", molNames(bigResOrd(k)), residuals(bigResOrd(k)));
end

%[text] Answer: Features with positive coefficients: HeavyAtomCount, NumRotatableBonds
%[text]    (More atoms and bonds tend to increase lipophilicity on average).
%[text]    Negative coefficients: TPSA, NumHDonors, NumHAcceptors (Polar groups decrease lipophilicity).
%[text]    Removing HeavyAtomCount (collinear with MolWt) and refitting shows RMSE remains almost unchanged,
%[text]    confirming feature redundancy.
%%
%[text] ## Let's Try 4: RF Feature Importance; Impact of MinLeafSize; OOB Convergence

rfModel = fitrensemble(Xs, y, ...
    Method="Bag", NumLearningCycles=200, ...
    Learners=templateTree(MinLeafSize=3));

importances = predictorImportance(rfModel);
[sortedImp, sortOrd] = sort(importances, "descend");
logInfo("RF Feature Importance Ranking:");
for d = 1:N_FEATS
    logInfo("  %d. %-22s  MDI = %.4f", d, FEAT_NAMES(sortOrd(d)), sortedImp(d));
end

%[text] Comparison of MinLeafSize
for mls = [1, 10, 20]
    rf_mls = fitrensemble(Xs, y, Method="Bag", NumLearningCycles=100, ...
        Learners=templateTree(MinLeafSize=mls));
    yPredTrain = predict(rf_mls, Xs);
    rmseTrain  = sqrt(mean((y - yPredTrain).^2));
    logInfo("MinLeafSize=%2d -> Training RMSE = %.3f", mls, rmseTrain);
end

%[text] OOB Error Convergence History
oobRmse = oobLoss(rfModel, "Mode", "cumulative");
figure("Name", "A03 OOB RMSE Convergence");
plot(oobRmse, "Color", [0.4 0.7 0.4], LineWidth=1.5);
xlabel("Number of Trees"); ylabel("OOB RMSE");
title("RF OOB Error Convergence");
grid on;
logInfo("OOB RMSE stabilizes around %d trees", ...
    find(abs(diff(oobRmse)) < 0.001, 1, "first"));

%[text] Answer: RF often selects polar surface area (TPSA) or HeavyAtomCount as the most important features.
%[text]    These capture both axes of ALogP variation (polarity and size).
%[text]    MinLeafSize=1 learns deep trees: Training RMSE approaches 0 (overfitting).
%[text]    MinLeafSize=20 prunes trees heavily: Both training and test RMSE are high.
%[text]    MinLeafSize=3−5 is a typical compromise for QSAR datasets.
%[text]    OOB error usually flattens around 150−200 trees.
%%
%[text] ## Let's Try 5: 5-Fold CV Comparison; Overfitting Gap; 10-Fold Sensitivity

K_FOLDS = 5;
rng(42);
cv = cvpartition(nMols, KFold=K_FOLDS);
yCV_lm = nan(nMols, 1);
yCV_rf = nan(nMols, 1);

for fold = 1:K_FOLDS
    trainIdx = training(cv, fold);
    testIdx  = test(cv, fold);
    trTbl = array2table(Xs(trainIdx, :), VariableNames=cellstr(FEAT_NAMES));
    trTbl.ALogP = y(trainIdx);
    lm_f = fitlm(trTbl, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));
    teTbl = array2table(Xs(testIdx, :), VariableNames=cellstr(FEAT_NAMES));
    yCV_lm(testIdx) = predict(lm_f, teTbl);
    rf_f = fitrensemble(Xs(trainIdx, :), y(trainIdx), ...
        Method="Bag", NumLearningCycles=200, Learners=templateTree(MinLeafSize=3));
    yCV_rf(testIdx) = predict(rf_f, Xs(testIdx, :));
    logProgress(fold, K_FOLDS, "CV Fold");
end

r2cv_lm  = 1 - sum((y-yCV_lm).^2) / sum((y-mean(y)).^2);
rmse_lm  = sqrt(mean((y-yCV_lm).^2));
r2cv_rf  = 1 - sum((y-yCV_rf).^2) / sum((y-mean(y)).^2);
rmse_rf  = sqrt(mean((y-yCV_rf).^2));

yTr_rf   = predict(rfModel, Xs);
rmse_rfTr = sqrt(mean((y-yTr_rf).^2));

logInfo("5-Fold CV — Linear Model: R2=%.3f  RMSE=%.3f", r2cv_lm, rmse_lm);
logInfo("5-Fold CV — Random Forest: R2=%.3f  RMSE=%.3f", r2cv_rf, rmse_rf);
logInfo("RF Overfitting Gap: Training RMSE=%.3f  CV RMSE=%.3f  Gap=%.3f", ...
    rmse_rfTr, rmse_rf, rmse_rf - rmse_rfTr);

%[text] Comparison with 10-Fold CV

rng(42);
cv10 = cvpartition(nMols, KFold=10);
yCV_rf10 = nan(nMols, 1);
for fold = 1:10
    trainIdx = training(cv10, fold);
    testIdx  = test(cv10, fold);
    rf10 = fitrensemble(Xs(trainIdx, :), y(trainIdx), ...
        Method="Bag", NumLearningCycles=200, Learners=templateTree(MinLeafSize=3));
    yCV_rf10(testIdx) = predict(rf10, Xs(testIdx, :));
    logProgress(fold, 10, "10-Fold CV");
end
r2cv_rf10 = 1 - sum((y-yCV_rf10).^2) / sum((y-mean(y)).^2);
rmse_rf10 = sqrt(mean((y-yCV_rf10).^2));
logInfo("10-Fold CV — Random Forest: R2=%.3f  RMSE=%.3f", r2cv_rf10, rmse_rf10);

%[text] Answer: For 200 FDA-approved drugs, RF typically achieves a generalized CV R2 > 0.8, while the linear model is around ~0.6-0.7.
%[text]    The RF overfitting gap (CV RMSE - Training RMSE) is usually greater than 0.3, larger than the linear model,
%[text]    indicating RF fits the training data more strongly.
%[text]    With n=200, the R2 for 5-Fold vs 10-Fold CV is similar, but
%[text]    10-Fold has less bias, while the variance of CV estimates increases.
%%
%[text] ## Let's Try 6: Comparison with RDKit Crippen Prediction; Domain Check; Rationale for Retraining
%[text] Compare Ethanol's QSAR prediction with RDKit Crippen value
molEt = emk.mol.fromSmiles("CCO");
sCrippen = emk.descriptor.calculate(molEt, "LogP");
logInfo("RDKit Crippen LogP for Ethanol: %.2f", sCrippen.LogP);

%[text] Final Model — Retraining with All Data
if r2cv_rf >= r2cv_lm
    finalModel = fitrensemble(Xs, y, Method="Bag", NumLearningCycles=200, ...
        Learners=templateTree(MinLeafSize=3));
    predictFn = @(xs) predict(finalModel, xs);
    modelLabel = "RF";
else
    ftAll = array2table(Xs, VariableNames=cellstr(FEAT_NAMES));
    ftAll.ALogP = y;
    finalModel = fitlm(ftAll, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));
    predictFn  = @(xs) predict(finalModel, array2table(xs, VariableNames=cellstr(FEAT_NAMES)));
    modelLabel = "Linear";
end
logInfo("Final Model: %s (Retrained with All Data)", modelLabel);

%[text] Domain Check for Ethanol
sEt = emk.descriptor.calculate(molEt, FEAT_NAMES);
xEt = zeros(1, N_FEATS);
for d = 1:N_FEATS; xEt(d) = sEt.(FEAT_NAMES(d)); end
xEtS = (xEt - Xmean) ./ Xstd;
inDomain = all(xEtS >= min(Xs,[],1) - 2 & xEtS <= max(Xs,[],1) + 2);
logInfo("Ethanol Domain Check: %s", string(ternary_(inDomain, "In Domain", "Out of Domain")));
logInfo("Ethanol QSAR Prediction: %.2f  (Crippen: %.2f)", predictFn(xEtS), sCrippen.LogP);

%[text] Domain Check for Large Molecule (Peptide)
molPep = emk.mol.fromSmiles("NCC(=O)NCC(=O)O");
sPep   = emk.descriptor.calculate(molPep, FEAT_NAMES);
xPep   = zeros(1, N_FEATS);
for d = 1:N_FEATS; xPep(d) = sPep.(FEAT_NAMES(d)); end
xPepS  = (xPep - Xmean) ./ Xstd;
inDomPep = all(xPepS >= min(Xs,[],1) - 2 & xPepS <= max(Xs,[],1) + 2);
logInfo("Ala-Gly Peptide Domain: %s", string(ternary_(inDomPep, "In Domain", "Out of Domain")));

%[text] Answer: The QSAR prediction for Ethanol and the Crippen value should be close (approximately -0.3 ~ -0.7).
%[text]    This is because the learning target ALogP itself is a value calculated by the Crippen method.
%[text]    Ethanol is a small simple molecule, different in chemical class from FDA-approved drugs (training data),
%[text]    but it may fall within the descriptor space.
%[text]    Dipeptides are likely to be judged out of domain by HeavyAtomCount or NumHDonors.
%[text]    Using predictions out of domain carries a risk of large errors.
%[text]    In CV, 20% of the data is held out in each fold, so the final model is retrained with 100% of the data to reduce bias.
logInfo("A03 Answer Completed.");

function out = ternary_(cond, a, b)
    if cond; out = a; else; out = b; end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]%   data: {"layout":"inline","rightPanelPercent":40}
%---
