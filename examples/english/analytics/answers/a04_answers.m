%[text] # A04 Answers: Drug vs Non-drug Classification
%[text] Reference answers for the "Try It Yourself" exercise in a04_drug_classification.m.
%[text] First, run a04_drug_classification.m or execute this file standalone
%[text] (rebuild all necessary variables from scratch).
addpath(genpath("src"));
emk.setup.initPython();
logInfo("A04 Answers: Setup Complete");

%[text] -- Rebuild workspace (standalone) --
FP_NBITS = 2048;
N_PCS    = 20;
TEST_RATIO = 0.20;

rawDrugs = readtable("data/list/fda_drugs.csv",          TextType="string");
rawChem  = readtable("data/list/everyday_chemicals.csv", TextType="string");
logInfo("FDA Drugs: %d rows | Everyday Chemicals: %d rows", ...
    height(rawDrugs), height(rawChem));

function fps = buildFps_(smiles, nBits)
    n   = numel(smiles);
    fps = cell(1, n);
    for k = 1:n
        try
            fps{k} = emk.fingerprint.morgan( ...
                emk.mol.fromSmiles(smiles(k)), NBits=nBits);
        catch
        end
    end
end

fps1 = buildFps_(rawDrugs.SMILES, FP_NBITS);
valid1 = ~cellfun(@isempty, fps1);
fps0 = buildFps_(rawChem.SMILES, FP_NBITS);
valid0 = ~cellfun(@isempty, fps0);

n1     = sum(valid1);  n0 = sum(valid0);
names1 = rawDrugs.Name(valid1);
names0 = rawChem.CommonName(valid0);
nTotal = n1 + n0;

allFps = [fps1(valid1), fps0(valid0)];
bitMat = zeros(nTotal, FP_NBITS, "logical");
for k = 1:nTotal
    bitMat(k, :) = emk.fingerprint.toArray(allFps{k});
end
labels   = [ones(n1, 1); zeros(n0, 1)];
allNames = [names1(:); names0(:)];

logInfo("Total: %d  (%d drugs, %d non-drugs)", nTotal, n1, n0);

%[text] PCA on all data (for visualization only -- not used for train/test)
%[text] Use 60 PCs so the scree curve shows a meaningful decay curve
pca_ws_viz = warning('off', 'stats:pca:ColRankDefX');
[pcCoeffAll, pcScoreAll, ~, ~, explainedAll] = pca(double(bitMat), NumComponents=60);
warning(pca_ws_viz);
%%
%[text] ## Let's Try 1: On-bit Comparison of Drugs vs Non-drugs; No-information Rate

onBits1 = sum(double(bitMat(1:n1, :)), 2);
onBits0 = sum(double(bitMat(n1+1:end, :)), 2);
logInfo("Average on-bits: Drugs=%.1f  Non-drugs=%.1f", mean(onBits1), mean(onBits0));
logInfo("No-information rate (all-drug classifier): %.3f", n1/nTotal);

%[text] Everyday chemicals that are also pharmaceuticals
pharmaList = ["aspirin","caffeine","ibuprofen","paracetamol","ethanol","nicotine"];
pharma_in_chem = intersect(lower(names0), pharmaList);
logInfo("Everyday chemicals that are pharmaceuticals: %s", ...
    strjoin(pharma_in_chem, ", "));

%[text] Answer: FDA drugs have more on-bits (~55-70) than small everyday chemicals
%[text]    (~25-45) because drug molecules are structurally more complex.
%[text]    No-information rate = 87% (160/184 in training) -- a model predicting
%[text]    "always drug" would score 87% accuracy.  This makes accuracy a poor
%[text]    metric: a model with 90% accuracy might just classify everything as drug.
%[text]    Recall for non-drugs (the minority class) is a more informative metric.
%%
%[text] ## Let's Try 2: PC1 vs PC2 Separation; Scree Curve; Position of Aspirin/Caffeine
%[text] Classes in PC space
figure("Name", "A04 Answers PCA");
hold on;
scatter(pcScoreAll(1:n1, 1), pcScoreAll(1:n1, 2), 15, [0.3 0.6 0.9], "filled", ...
    MarkerFaceAlpha=0.4, DisplayName="FDA Approved Drugs (Class 1)");
scatter(pcScoreAll(n1+1:end, 1), pcScoreAll(n1+1:end, 2), 40, [0.9 0.4 0.3], "filled", ...
    MarkerFaceAlpha=0.8, DisplayName="Everyday Chemicals (Class 0)");

%[text] Aspirin and caffeine
aspIdx = find(allNames == "aspirin");
cafIdx = find(allNames == "caffeine");
if ~isempty(aspIdx)
    scatter(pcScoreAll(aspIdx, 1), pcScoreAll(aspIdx, 2), 100, "rx", LineWidth=2.5, ...
        DisplayName="aspirin");
end
if ~isempty(cafIdx)
    scatter(pcScoreAll(cafIdx, 1), pcScoreAll(cafIdx, 2), 100, "g+", LineWidth=2.5, ...
        DisplayName="caffeine");
end
hold off;
xlabel("PC1"); ylabel("PC2"); title("PCA of All Data (20 PCs Used in Model)");
legend(Location="best"); grid on;

%[text] Scree curve
figure("Name", "A04 Scree (60 PCs)");
nShow = min(60, numel(explainedAll));
plot(1:nShow, explainedAll(1:nShow), "-o", Color=[0.3 0.6 0.9], LineWidth=1.5);
xlabel("PC"); ylabel("Explained Variance (%)");
title("Scree Curve -- Fingerprint PCA"); grid on;
cumExp60 = cumsum(explainedAll(1:nShow));
logInfo("Cumulative variance: PC5=%.1f%%  PC20=%.1f%%  PC50=%.1f%%", ...
    cumExp60(5), cumExp60(20), cumExp60(min(50,nShow)));

%[text] Answer: The two classes overlap substantially in the first two PCs because
%[text]    drugs and non-drugs share many common substructures (benzene rings,
%[text]    amide bonds, etc.).  Drug-like space is not linearly separable from
%[text]    everyday chemicals -- hence the need for SVM or RF.
%[text]    Aspirin and caffeine (everyday drugs) plot inside the drug cluster,
%[text]    consistent with their label noise.
%[text]    The scree curve for ECFP4 fingerprints decays very slowly -- 20 PCs
%[text]    typically explain only 30-40% of variance (high-dimensional sparse data).
%%
%[text] ## Let's Try 3: Verification of Stratification; Impact of Non-stratified Split

rng(42);
cvPart  = cvpartition(labels, Holdout=TEST_RATIO);
trainIdx = training(cvPart);
testIdx  = test(cvPart);

drugFracTrain = sum(labels(trainIdx)==1) / sum(trainIdx);
drugFracTest  = sum(labels(testIdx)==1)  / sum(testIdx);
logInfo("Drug fraction: train=%.3f  test=%.3f  (difference=%.3f)", ...
    drugFracTrain, drugFracTest, abs(drugFracTrain - drugFracTest));

%[text] Non-stratified comparison
cvPart2 = cvpartition(nTotal, Holdout=0.2);
nonDrugsInTest2 = sum(labels(test(cvPart2)) == 0);
logInfo("Non-stratified split: non-drugs in test set = %d (can be 0!)", nonDrugsInTest2);

%[text] Seed sensitivity
logInfo("Seed sensitivity (test accuracy may vary by 2-5%% across seeds):");
for seed = [1, 7, 2026]
    rng(seed);
    cv_s = cvpartition(labels, Holdout=TEST_RATIO);
    df   = abs(sum(labels(training(cv_s))==1)/sum(training(cv_s)) - ...
               sum(labels(test(cv_s))==1)/sum(test(cv_s)));
    logInfo("  seed=%d  drug-fraction gap = %.3f", seed, df);
end

%[text] Answer: Stratified split keeps drug-fraction gap < 2-3%.
%[text]    Non-stratified can allocate 0 non-drugs to the test set by chance,
%[text]    making test accuracy 100% trivially (all test samples are drug).
%[text]    Seed sensitivity: single-split test accuracy varies by 2-5% for this
%[text]    small dataset -- a strong argument for CV.
%%
%[text] ## Let's Try 4: Confusion Matrix Cell; No-Information Rate; Cost-Sensitive SVM; RBF Kernel
%[text] Rebuild train/test sets with correct PCA (train-only fit)
trainMean = mean(double(bitMat(trainIdx, :)), 1);
Xtrain_c  = double(bitMat(trainIdx, :)) - trainMean;
pca_ws_tr = warning('off', 'stats:pca:ColRankDefX');
[pcCoeffTr, XTrain, ~, ~, explainedTr] = pca(Xtrain_c, NumComponents=N_PCS);
warning(pca_ws_tr);
Xtest_c = double(bitMat(testIdx, :)) - trainMean;
XTest   = Xtest_c * pcCoeffTr;
yTrain  = labels(trainIdx);
yTest   = labels(testIdx);

%[text] Default SVM
svmModel = fitcsvm(XTrain, yTrain, KernelFunction="linear", BoxConstraint=1, ...
    Standardize=true, ClassNames=[0; 1]);
[yPred_svm, ~] = predict(svmModel, XTest);
cm_svm   = confusionmat(yTest, yPred_svm, Order=[0 1]);
acc_svm  = sum(diag(cm_svm)) / sum(cm_svm(:));
logInfo("SVM: TN=%d FP=%d FN=%d TP=%d  Accuracy=%.3f", ...
    cm_svm(1,1), cm_svm(1,2), cm_svm(2,1), cm_svm(2,2), acc_svm);
logInfo("No-Information Rate: %.3f  SVM exceeds it: %s", n1/nTotal, ...
    string(ternary_(acc_svm > n1/nTotal, "YES", "NO")));

%[text] Cost-sensitive SVM (higher misclassification cost for non-drugs)
costMat = [0, n1/n0; 1, 0];   % cost(true=0, pred=1) = n1/n0
svmCost = fitcsvm(XTrain, yTrain, KernelFunction="linear", BoxConstraint=1, ...
    Standardize=true, ClassNames=[0; 1], Cost=costMat);
[yPred_cost, ~] = predict(svmCost, XTest);
cm_cost = confusionmat(yTest, yPred_cost, Order=[0 1]);
rec0_default = cm_svm(1,1)  / max(cm_svm(1,1)  + cm_svm(1,2),  1);
rec0_cost    = cm_cost(1,1) / max(cm_cost(1,1) + cm_cost(1,2), 1);
logInfo("Non-drug recall: default=%.3f  cost-weighted=%.3f", rec0_default, rec0_cost);

%[text] RBF kernel
svmRBF = fitcsvm(XTrain, yTrain, KernelFunction="rbf", Standardize=true, ClassNames=[0;1]);
[yPred_rbf, ~] = predict(svmRBF, XTest);
acc_rbf  = mean(yPred_rbf == yTest);
logInfo("RBF SVM Accuracy: %.3f", acc_rbf);

%[text] Answer: FN (missed drugs) is more costly in drug screening: you would fail to
%[text]    identify a promising compound.  FP (false alarm: non-drug classified as drug)
%[text]    wastes resources but is recoverable in follow-up assays.
%[text]    Cost weighting improves non-drug recall (less FP) at the cost of overall accuracy.
%[text]    RBF kernel may improve accuracy slightly by capturing non-linear boundaries.
%%
%[text] ## Let's Try 5: RF vs SVM Comparison; Impact of Prior=uniform; Feature Importance PC

rfClf = fitcensemble(XTrain, yTrain, Method="Bag", NumLearningCycles=200, ...
    Learners=templateTree(MinLeafSize=1), Prior="uniform", ClassNames=[0; 1]);
[yPred_rf, ~] = predict(rfClf, XTest);
cm_rf    = confusionmat(yTest, yPred_rf, Order=[0, 1]);
acc_rf   = sum(diag(cm_rf)) / sum(cm_rf(:));
rec_rf0  = cm_rf(1,1) / max(cm_rf(1,1) + cm_rf(1,2), 1);
f1_rf    = 2 * cm_rf(2,2) / max(2*cm_rf(2,2) + cm_rf(1,2) + cm_rf(2,1), 1);

logInfo("RF (uniform prior): Accuracy=%.3f  NonDrug-Recall=%.3f  F1=%.3f", ...
    acc_rf, rec_rf0, f1_rf);

%[text] Without Prior=uniform
rfNoUniform = fitcensemble(XTrain, yTrain, Method="Bag", NumLearningCycles=200, ...
    Learners=templateTree(MinLeafSize=1), ClassNames=[0; 1]);
[yPred_rfN, ~] = predict(rfNoUniform, XTest);
cm_rfN   = confusionmat(yTest, yPred_rfN, Order=[0, 1]);
acc_rfN  = sum(diag(cm_rfN)) / sum(cm_rfN(:));
rec_rfN0 = cm_rfN(1,1) / max(cm_rfN(1,1) + cm_rfN(1,2), 1);
logInfo("RF (no uniform prior): Accuracy=%.3f  NonDrug-Recall=%.3f", acc_rfN, rec_rfN0);

%[text] Feature importance
featureImp = predictorImportance(rfClf);
[~, topPC] = max(featureImp);
logInfo("Highest MDI importance: PC%d", topPC);
logInfo("PC1 importance rank: %d / %d", ...
    find(sort(featureImp,"descend") == featureImp(1), 1), N_PCS);

%[text] Answer: Prior=uniform improves non-drug recall (the minority class gets equal
%[text]    weight during training) at the cost of slightly lower overall accuracy.
%[text]    The most important PC for classification is often NOT PC1: variance
%[text]    and discriminative power are different.  PCs that capture drug-specific
%[text]    substructures (ring systems, H-bond patterns) may have higher MDI
%[text]    even if they explain little total variance.
%%
%[text] ## Let's Try 6: ROC AUC; Optimal Threshold; Null Model Check; 5-Fold CV AUC

[~, scoreSVM] = predict(svmModel, XTest);
[~, scoreRF]  = predict(rfClf,   XTest);

scSVM = scoreSVM(:, 2);
scRF  = scoreRF(:, 2);

[xSVM, ySVM, ~, aucSVM] = perfcurve(yTest, scSVM, 1);
[xRF,  yRF,  ~, aucRF]  = perfcurve(yTest, scRF,  1);

logInfo("AUC: SVM=%.3f  RF=%.3f", aucSVM, aucRF);

%[text] Optimal threshold (closest to top-left corner)
[~, optIdxSVM] = min((1-ySVM).^2 + xSVM.^2);
logInfo("SVM optimal threshold index: FPR=%.3f  TPR=%.3f", ...
    xSVM(optIdxSVM), ySVM(optIdxSVM));

%[text] Null model (shuffled labels)
rng(99);
shuffIdx = randperm(numel(yTest));
[~, ~, ~, aucNull] = perfcurve(yTest, scSVM(shuffIdx), 1);
logInfo("Null model AUC (shuffled): %.3f  (should be ~0.5)", aucNull);

%[text] 5-fold CV AUC
cvSVM = crossval(svmModel, KFold=5);
[yCV_svm, scoresCV] = kfoldPredict(cvSVM);
[~, ~, ~, aucCV_svm] = perfcurve(yTrain, scoresCV(:, 2), 1);
logInfo("SVM 5-fold CV AUC: %.3f  (single-split AUC: %.3f)", aucCV_svm, aucSVM);

%[text] Answer: AUC > 0.5 for both SVM and RF confirms the models have genuine discriminative
%[text]    power.  The null model AUC ~0.5 verifies there is no data leakage.
%[text]    The optimal threshold is typically < 0.5 for imbalanced data (lower threshold
%[text]    = more positives predicted = higher recall for minority class).
%[text]    CV AUC is a more reliable estimate than single-split AUC for n=230.
logInfo("A04 Answer Complete.");

function out = ternary_(cond, a, b)
    if cond; out = a; else; out = b; end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]%   data: {"layout":"inline","rightPanelPercent":40}
%---
