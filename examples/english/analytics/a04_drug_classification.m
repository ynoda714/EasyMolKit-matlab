%[text] # A04: Drug Classification — FDA Approved Drugs vs Everyday Chemicals
%[text] EasyMolKit Analytics — Layer 3
%[text] 
%[text] "Does this white powder look like a drug?" — In forensic settings, instant judgment based solely on molecular structure is sometimes required.
%[text] Convert 200 FDA approved drugs and 30 everyday chemicals into 2048-bit fingerprints, and train "drug-likeness" using SVM and Random Forest.
%[text] Compare accuracy, F1 score, and ROC/AUC to learn how to choose appropriate evaluation metrics for class-imbalanced data.
%[text] In this script, experience preprocessing steps for p \>\> n regime and practical evaluation methods for classification models.
%[text] 
%[text] **What you will learn in this tutorial**
%[text] - Understand the "p \>\> n" problem in fingerprint-based machine learning
%[text] - Learn how to apply PCA as preprocessing for supervised classification
%[text] - Perform binary molecular classification using `fitcsvm()` and `fitcensemble()`
%[text] - Interpret and understand confusion matrix, accuracy, recall, and F1 score
%[text] - Interpret ROC curve and AUC as threshold-independent performance metrics
%[text] - Recognize the impact of class imbalance on accuracy-based evaluation \
%[text] 
%[text] **Prerequisites**
%[text] - Completion of F03 (Fingerprints) and F04 (Similarity) is required
%[text] - Recommended: Understand context with A01 (PCA) and A02 (Clustering)
%[text] - Statistics and Machine Learning Toolbox (`fitcsvm`, `fitcensemble`, `perfcurve`) is required
%[text] - Internet connection is not required \
%[text] 
%[text] Estimated time: 30-45 minutes | Execution method: Run each section with Ctrl+Enter
%[text] 
%[text] **Data:** `data/list/fda_drugs.csv` (200 FDA approved drugs, ChEMBL CC-BY-SA 3.0), `data/list/everyday_chemicals.csv` (30 household chemicals, PubChem CC0)
%[text] 
%[text] **Note on Label Quality**
%[text] Use dataset origin as class labels (drug vs non-drug).
%[text] Some everyday chemicals like Aspirin, Caffeine, and Ibuprofen are actually drugs as well—this is intentional label noise.
%[text] In actual QSAR classification, experimentally measured activity values are used.
%[text] Incomplete labels limit the maximum achievable AUC (see Section 6).
%[text] 
%[text] **References**
%[text] - Cortes & Vapnik (1995) Support-vector networks. *Machine Learning* 20:273-297.
%[text] - Breiman (2001) Random forests. *Machine Learning* 45:5-32.
%[text] - Fawcett (2006) An introduction to ROC analysis. *Pattern Recognition Letters* 27:861-874.
%[text] - Rogers & Hahn (2010) Extended-connectivity fingerprints. *J Chem Inf Model* 50:742-754. \
%%
%[text] ## Section 0: Environment Setup
logSection("A04", "Section 0: Setup", "Analytics L3");
% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
emk.setup.initPython();
%[text] Warm up the Python/RDKit process before executing the main procedure
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;
%%
%[text] ## Section 1: Dataset Loading and Fingerprint Calculation
%[text] 
%[text] Setup is complete. First, load two types of datasets and convert each molecule into a 2048-bit fingerprint.
%[text] This bit vector will be used as input features for machine learning.
%[text] 
%[text] ### Concept: Fingerprints as Binary Feature Vectors for Machine Learning
%[text] Morgan / ECFP4 fingerprints are 2048-bit binary vectors. Each bit encodes the presence of circular substructures (radius 2 from each heavy atom, i.e., up to 2 bonds from the central atom). The "4" in ECFP4 indicates the diameter (= 2×radius).
%[text] 
%[text] **Role in Machine Learning:** Each bit becomes a feature (2048 features per molecule), and molecules are represented as rows in a 0/1 feature matrix. The classifier learns bit patterns that distinguish two classes.
%[text] 
%[text] **Challenge — "p \>\> n" Regime:** With about 230 molecules (n) and 2048 features (p), there are more features than samples. This almost guarantees linear separability but poses a risk of overfitting. Many bits are nearly constant and uninformative, making SVM hyperplane estimation unstable. Dimensionality reduction via PCA is considered as a solution (Section 2).
%[text] 
%[text] **Label Convention:** 1 = FDA-approved drug compound (from ChEMBL) / 0 = Everyday household chemical (from PubChem)
logSection("A04", "Section 1: Dataset Loading and Fingerprint Calculation", "Analytics L3");
FP_RADIUS = 2;
FP_NBITS  = 2048;
FDA_FILE  = "data/list/fda_drugs.csv";
CHEM_FILE = "data/list/everyday_chemicals.csv";

fdaTbl  = readtable(FDA_FILE,  TextType="string");
chemTbl = readtable(CHEM_FILE, TextType="string");
logInfo("FDA Approved Drugs: %d rows  |  Everyday Chemicals: %d rows", ...
    height(fdaTbl), height(chemTbl));
%[text] Parse FDA-approved drugs and calculate fingerprints (Class 1).
fps1 = {}; names1 = string.empty; failed1 = 0;
for k = 1:height(fdaTbl)
    try
        mol = emk.mol.fromSmiles(fdaTbl.SMILES(k));
        fps1{end+1}   = emk.fingerprint.morgan(mol, Radius=FP_RADIUS, NBits=FP_NBITS); %#ok<AGROW>
        names1(end+1) = fdaTbl.Name(k); %#ok<AGROW>
    catch
        failed1 = failed1 + 1;
    end
end
logInfo("FDA Class:     %d parsed successfully, %d failed", numel(fps1), failed1);
%[text] Parse everyday chemicals and calculate fingerprints (Class 0).
fps0 = {}; names0 = string.empty; failed0 = 0;
for k = 1:height(chemTbl)
    try
        mol = emk.mol.fromSmiles(chemTbl.SMILES(k));
        fps0{end+1}   = emk.fingerprint.morgan(mol, Radius=FP_RADIUS, NBits=FP_NBITS); %#ok<AGROW>
        names0(end+1) = chemTbl.CommonName(k); %#ok<AGROW>
    catch
        failed0 = failed0 + 1;
    end
end
logInfo("Non-drug Class:     %d parsed successfully, %d failed", numel(fps0), failed0);
%[text] Assemble the bit matrix and label vector.
n1     = numel(fps1);
n0     = numel(fps0);
nTotal = n1 + n0;

bitMat = zeros(nTotal, FP_NBITS, "single");
allFps = [fps1, fps0];
for j = 1:nTotal
    bitMat(j, :) = single(emk.fingerprint.toArray(allFps{j}));
end

labels   = [ones(n1, 1); zeros(n0, 1)];    % 1=drug, 0=non-drug
allNames = [names1(:); names0(:)];

logInfo("Total: %d molecules (%d drugs, %d non-drugs)", nTotal, n1, n0);
logInfo("Class Balance: Drugs %.1f%% vs Non-drugs %.1f%%", ...
    100*n1/nTotal, 100*n0/nTotal);
%[text] **💡 Observation Point 1**
%[text] Let's check the average number of on-bits for drugs and non-drugs.
%[text] Refer to the following code:
%[text] onBits1 = sum(bitMat(1:n1, :), 2);           % Count per drug
%[text] onBits0 = sum(bitMat(n1+1:end, :), 2);        % Count per non-drug
%[text] \[mean(onBits1), mean(onBits0)\]
%[text] Interpret whether drug molecules have more or fewer structural features.
%[text] Check the accuracy of a naive classifier that always predicts "drug."
%[text] What would that accuracy be? (Answer: n1/nTotal)
%[text] Consider why accuracy is a misleading metric for imbalanced class data.
%[text] Identify everyday chemicals that are also drugs.
%[text] pharma\_in\_chem = intersect(lower(names0), \["aspirin","caffeine","ibuprofen",...
%[text] "paracetamol","ethanol","nicotine"\]);
%[text] These cause intentional label noise.
% ... (Try writing code here)
%%
%[text] ## Section 2: Dimensionality Reduction using PCA
%[text] The fingerprint calculation is complete. With 2048 dimensions, there is a high risk of overfitting, so we will use PCA to reduce dimensions.
%[text] We will also visualize the class distribution.
%[text] 
%[text] ### Concept: Why perform PCA before classification?
%[text] With 2048 binary features and about 230 molecules, we are in a p \>\> n situation. A linear classifier can achieve 100% training accuracy even with random labels, indicating a risk of overfitting. Moreover, many of the 2048 bits are almost constant and carry no information.
%[text] 
%[text] PCA projects the 2048-bit vector into N\_PCS dimensions in the direction of maximum variance (similar to the approach for descriptors in A01). This automatically removes nearly constant bits, retains structural variations that best separate the two classes, and reduces the feature matrix from N × 2048 to N × N\_PCS.
%[text] 
%[text] **Choosing N\_PCS:** In fingerprint-based machine learning, 10~50 principal components (PCs) are recommended in the literature. Here, we use N\_PCS=20.
%[text] 
%[text] **Important:** PCA should be fitted only on the training data, and the same transformation must be applied to test and new molecules. Fitting PCA on the entire dataset before splitting is "data leakage."
logSection("A04", "Section 2: Dimensionality Reduction using PCA", "Analytics L3");
N_PCS = 20;
%[text] Center the features (scaling is unnecessary — binary features are already in the \[0,1\] range).
bitMean   = mean(bitMat, 1);
bitMat_c  = double(bitMat) - bitMean;
%[text] Suppress expected warnings: binary fingerprint columns are linearly dependent.
%[text] (Bits that are almost constant have near-zero variance. The pca() function reduces the number of T^2 components to handle this — the results are correct.)
pca_ws = warning('off', 'stats:pca:ColRankDefX');
[pcCoeff, pcScore, ~, ~, explained] = pca(bitMat_c, NumComponents=N_PCS);
warning(pca_ws);

logInfo("PCA: Top %d PCs explain %.1f%% of the fingerprint variance", ...
    N_PCS, sum(explained(1:N_PCS)));
%[text] Visualize the first two principal components (PCs) color-coded by class.
figure("Name", "A04 Fingerprint PCA");
hold on;
scatter(pcScore(labels==1, 1), pcScore(labels==1, 2), 40, [0.3 0.6 0.9], ...
    "filled", MarkerFaceAlpha=0.6, DisplayName="FDA Approved Drugs (Class 1)");
scatter(pcScore(labels==0, 1), pcScore(labels==0, 2), 70, [0.9 0.4 0.3], ...
    "^", "filled", MarkerFaceAlpha=0.9, DisplayName="Household Chemicals (Class 0)");
hold off;
xlabel(sprintf("PC1 (%.1f%%)", explained(1)));
ylabel(sprintf("PC2 (%.1f%%)", explained(2)));
title("ECFP4 Fingerprint PCA -- Drugs vs Non-drugs");
legend(Location="best");
grid on;
%[text] **💡 Observation Point 2**
%[text] Check if the two classes are visually separated in the PC1 and PC2 plot.
%[text] If there is overlap, it indicates that the two classes share many structural features.
%[text] Check which class has greater variance (spread) in PC space.
%[text] Try N\_PCS = 5, 10, 50 and see how the explained variance changes.
%[text] Plot the scree curve: figure; plot(1:numel(explained), explained, "-o");
%[text] Identify where the elbow (diminishing returns point) is.
%[text] Aspirin and Caffeine are both household chemicals and drugs.
%[text] Find their positions in PC space.
%[text] aspIdx = find(allNames == "aspirin");
%[text] scatter(pcScore(aspIdx,1), pcScore(aspIdx,2), 100, "rx", LineWidth=2);
%[text] Check which cluster, drug or non-drug, they are closer to.
% ... (Try writing code here)
%%
%[text] ## Section 3: Stratified Training/Test Split
%[text] 
%[text] The feature space is organized using PCA (Principal Component Analysis). Next, we will split the data into training and test sets.
%[text] Due to class imbalance, we use stratified splitting to maintain the ratio of each class.
%[text] 
%[text] ### Concept: Stratified Splitting to Address Class Imbalance
%[text] A simple random 80/20 split risks biasing the test samples of minority classes. Stratified splitting divides each class independently, maintaining the same train:test ratio within each class as in the entire dataset.
%[text] 
%[text] **Correct PCA Procedure:** In a production environment, PCA should be fitted only on the training data, and the same transformation should be applied to the test data. Fitting PCA on the entire dataset before splitting is "data leakage."
%[text] 
logSection("A04", "Section 3: Stratified Training/Test Split", "Analytics L3");
TEST_RATIO = 0.20;
rng(42);
cvPart   = cvpartition(labels, Holdout=TEST_RATIO);
trainIdx = training(cvPart);
testIdx  = test(cvPart);
%[text] Correct method: Fit PCA only on the training data and apply the transformation to both sets.
trainMean = mean(double(bitMat(trainIdx, :)), 1);
Xtrain_c  = double(bitMat(trainIdx, :)) - trainMean;
pca_ws2 = warning('off', 'stats:pca:ColRankDefX');
[pcCoeffTr, XTrain, ~, ~, explainedTr] = pca(Xtrain_c, NumComponents=N_PCS);
warning(pca_ws2);
%[text] Apply the same transformation to the test data.
Xtest_c = double(bitMat(testIdx, :)) - trainMean;
XTest   = Xtest_c * pcCoeffTr;   % project onto training PCs

yTrain = labels(trainIdx);
yTest  = labels(testIdx);

logInfo("Training: %d samples (%d drugs, %d non-drugs)", ...
    sum(trainIdx), sum(yTrain==1), sum(yTrain==0));
logInfo("Testing: %d samples (%d drugs, %d non-drugs)", ...
    sum(testIdx), sum(yTest==1), sum(yTest==0));
logInfo("Training PCA: Top %d PCs explain %.1f%% of training variance", ...
    N_PCS, sum(explainedTr(1:N_PCS)));
%[text] **💡 Observation Point 3**
%[text] Verify stratification: Check if the class ratios in the training and test sets are similar.
%[text] Drug ratio (training): sum(yTrain==1)/numel(yTrain)
%[text] Drug ratio (test): sum(yTest==1)/numel(yTest)
%[text] These should be nearly equal (within 2-3%).
%[text] Let's see what happens without stratification.
%[text] cvPart2 = cvpartition(nTotal, Holdout=0.2);  % Non-stratified (integer n)
%[text] sum(labels(test(cvPart2)) == 0)   % Non-drugs in test set -- may become 0!
%[text] Change the rng seed (e.g., rng(1), rng(7)) and rerun to observe test accuracy fluctuations.
%[text] Observe how much the test accuracy fluctuates. This fluctuation is the "split noise" inherent in a single train/test split. CV (Observation Point 5) can reduce this.
% ... (Try writing code here)
%%
%[text] ## Section 4: SVM Classifier (fitcsvm)
%[text] 
%[text] The training-test split is complete. Next, we will perform classification using an SVM that learns the maximum margin boundary.
%[text] Let's examine the prediction details using a confusion matrix.
%[text] 
%[text] ### Concept: Support Vector Machine (SVM) for Binary Classification
%[text] SVM finds the maximum margin hyperplane that separates two classes. Maximizing the margin $= 2 / ||w||$ is equivalent to minimizing $||w||^2$. Points on the margin boundary are "support vectors" and define the decision boundary.
%[text] 
%[text] **Key Parameters:** `BoxConstraint` (C) — Trade-off between margin width and training error (small C for wide margin, large C for narrow margin) / `KernelFunction` — Feature space includes `"linear"` / `Standardize` — Normalize features to mean 0, unit variance
%[text] 
%[text] **Handling Class Imbalance:** If 87% are drugs and 13% are non-drugs, the default SVM may ignore the minority class (non-drugs) as it tries to minimize the total number of misclassifications. This is addressed in 💡 Observation Point 4.
%[text] 
logSection("A04", "Section 4: SVM Classifier (fitcsvm)", "Analytics L3");
svmModel = fitcsvm(XTrain, yTrain, ...
    KernelFunction="linear", ...
    BoxConstraint=1.0, ...
    Standardize=true, ...
    ClassNames=[0; 1]);

[yPred_svm, scores_svm] = predict(svmModel, XTest);
%[text] Classification Evaluation Metrics (Manual Calculation — No Additional Toolbox Functions Required)
cm_svm   = confusionmat(yTest, yPred_svm, Order=[0, 1]);
acc_svm  = sum(diag(cm_svm)) / sum(cm_svm(:));
prec_svm = cm_svm(2,2) / max(cm_svm(1,2) + cm_svm(2,2), 1);  % avoid /0
rec_svm  = cm_svm(2,2) / max(cm_svm(2,1) + cm_svm(2,2), 1);
f1_svm   = 2 * prec_svm * rec_svm / max(prec_svm + rec_svm, eps);

logInfo("SVM (Linear, C=1) -- Test Set:");
logInfo("  Accuracy=%.3f  Precision=%.3f  Recall=%.3f  F1=%.3f", ...
    acc_svm, prec_svm, rec_svm, f1_svm);
logInfo("  Confusion Matrix: TN=%d FP=%d FN=%d TP=%d", ...
    cm_svm(1,1), cm_svm(1,2), cm_svm(2,1), cm_svm(2,2));

figure("Name", "A04 SVM Confusion Matrix");
confusionchart(categorical(yTest, [0 1], ["Non-drug","Drug"]), ...
    categorical(yPred_svm, [0 1], ["Non-drug","Drug"]), ...
    ColumnSummary="column-normalized", RowSummary="row-normalized", ...
    Title="SVM Confusion Matrix");
logInfo("Displayed SVM Confusion Matrix");
%[text] **💡 Observation Point 4**
%[text] Let's confirm what the four cells mean.
%[text] - TP = Drugs correctly labeled as drugs (True Positive)
%[text] - TN = Non-drugs correctly labeled as non-drugs (True Negative)
%[text] - FP = Non-drugs incorrectly labeled as drugs (False Positive — "False Alarm")
%[text] - FN = Drugs incorrectly labeled as non-drugs (False Negative — "Missed Drug") \
%[text] In the context of drug screening, consider which type of error is more costly.
%[text] Consider the no-information rate. A "lazy" classifier always predicts "drug."
%[text] Its accuracy is n1/nTotal. Check if the SVM accuracy exceeds the no-information rate.
%[text] (Otherwise, the SVM adds no value over a naive baseline!)
%[text] Use a cost-sensitive SVM to handle class imbalance.
%[text] w = \[n1/nTotal, n0/nTotal\];   % Inverse frequency weights
%[text] svmW = fitcsvm(XTrain, yTrain, KernelFunction="linear", ...
%[text] Cost=\[0 n1/n0; 1 0\]);   % cost(true=0, pred=1) = n1/n0
%[text] \[yW, ~\] = predict(svmW, XTest);
%[text] Observe how cost weighting changes the recall of the minority class.
%[text] Try KernelFunction="rbf". Check if the nonlinear SVM improves accuracy.
%[text] (Note: The RBF kernel in 20-dimensional PCA space is still nonlinear)
% ... (Try writing code here)
%%
%[text] ## Section 5: Random Forest Classifier (fitcensemble)
%[text] 
%[text] We have built an SVM. Next, let's try Random Forest, an ensemble of decision trees.
%[text] Let's examine the differences from SVM, such as handling class imbalance and feature importance.
%[text] 
%[text] Concept: Random Forest for Classification
%[text] Similar to the tree ensemble approach in A03 regression, but each tree votes for a class.
%[text] 1. Build B trees, each using bootstrap samples and a random subset of features.
%[text] 2. Each tree outputs a class vote (0 or 1).
%[text] 3. The final prediction is determined by majority vote (> 50% of trees agree).
%[text] 4. Class probabilities are represented by the proportion of trees voting for that class. \
%[text] (Smoother than hard thresholds, directly usable in ROC)
%[text] 
%[text] Advantages over SVM:
%[text] - Obtain class probabilities without post-calibration.
%[text] - Less sensitive to feature scaling, no need for Standardize.
%[text] - Prior="uniform" is a simple way to handle class imbalance.
%[text] - Feature importance using split gain is interpretable. \
%[text] 
%[text] Prior="uniform":
%[text] Set prior class probabilities to 0.5/0.5 regardless of actual 87%/13% split.
%[text] This treats both classes as equally important during training, improving recall for non-drugs.
%[text] 
%[text] MinLeafSize=1: Fully grown trees (maximum complexity, slight overfitting).
%[text] For larger datasets, increase MinLeafSize to reduce overfitting.
logSection("A04", "Section 5: Random Forest Classifier (fitcensemble)", "Analytics L3");
rfClf = fitcensemble(XTrain, yTrain, ...
    Method="Bag", ...
    NumLearningCycles=200, ...
    Learners=templateTree(MinLeafSize=1), ...
    Prior="uniform", ...     % Handle class imbalance
    ClassNames=[0; 1]);

[yPred_rf, scores_rf] = predict(rfClf, XTest);

cm_rf   = confusionmat(yTest, yPred_rf, Order=[0, 1]);
acc_rf  = sum(diag(cm_rf)) / sum(cm_rf(:));
prec_rf = cm_rf(2,2) / max(cm_rf(1,2) + cm_rf(2,2), 1);
rec_rf  = cm_rf(2,2) / max(cm_rf(2,1) + cm_rf(2,2), 1);
f1_rf   = 2 * prec_rf * rec_rf / max(prec_rf + rec_rf, eps);

logInfo("Random Forest (200 trees, uniform prior) -- Test Set:");
logInfo("  Accuracy=%.3f  Precision=%.3f  Recall=%.3f  F1=%.3f", ...
    acc_rf, prec_rf, rec_rf, f1_rf);
logInfo("  Confusion Matrix: TN=%d FP=%d FN=%d TP=%d", ...
    cm_rf(1,1), cm_rf(1,2), cm_rf(2,1), cm_rf(2,2));

figure("Name", "A04 RF Confusion Matrix");
confusionchart(categorical(yTest, [0 1], ["Non-drug","Drug"]), ...
    categorical(yPred_rf, [0 1], ["Non-drug","Drug"]), ...
    ColumnSummary="column-normalized", RowSummary="row-normalized", ...
    Title="Random Forest Confusion Matrix");
logInfo("Displayed RF Confusion Matrix");
%[text] Feature Importance (per PCA Component)
featureImp = predictorImportance(rfClf);
figure("Name", "A04 RF Feature Importance (PCA Components)");
bar(featureImp, FaceColor=[0.4 0.7 0.4]);
xlabel("PCA Component"); ylabel("MDI Importance");
title("RF Feature Importance for PC1..PC20");
grid on;
%[text] **💡 Observation Point 5**
%[text] Compare the confusion matrices of RF and SVM,
%[text] and check which has higher recall for the non-drug class (fewer FN).
%[text] Check which has higher precision for the drug class (fewer FP).
%[text] Remove Prior="uniform" from rfClf and refit,
%[text] and check if recall for non-drugs decreases or overall accuracy improves.
%[text] This demonstrates the trade-off between accuracy and recall under class imbalance.
%[text] Check the PC component with the highest RF feature importance,
%[text] and see if it is PC1, which explains the most variance, or a later PC.
%[text] (High variance does not necessarily mean high discriminative power)
%[text] Change MinLeafSize from 1 to 10 and see how test accuracy changes.
%[text] Also check training accuracy.
%[text] \[yTrain\_pred, ~\] = predict(rfClf, XTrain);
%[text] mean(yTrain\_pred == yTrain)   % Training accuracy (should be ~1 with MinLeafSize=1)
% ... (Try writing code here)
%%
%[text] ## Section 6: Comparison of ROC Curves
%[text] 
%[text] We have two classifiers ready. Here, we will fairly compare both models using AUC (Area Under the Curve), a comprehensive metric independent of thresholds.
%[text] In data with class imbalance, AUC is a more reliable evaluation metric than accuracy.
%[text] 
%[text] Concept: Receiver Operating Characteristic (ROC) Curve
%[text] A binary classifier outputs scores (probabilities or decision values).
%[text] The default threshold (0.5) converts scores into binary predictions.
%[text] Changing the threshold creates trade-offs:
%[text] - High threshold -> Fewer predicted positives -> Lower recall, higher precision
%[text] - Low threshold -> More predicted positives -> Higher recall, lower precision \
%[text] 
%[text] ROC Curve: A plot of true positive rate (recall/sensitivity) versus false positive rate (1 - specificity) across all thresholds.
%[text] 
%[text] AUC (Area Under the ROC Curve):
%[text] AUC = 1.0  Perfect classifier (TPR=1 at all FPR)
%[text] AUC = 0.5  Random classifier (diagonal -- no discrimination ability)
%[text] AUC = 0.0  Perfectly inverse (perfect if all predictions are flipped)
%[text] 
%[text] Interpretation of AUC in Cheminformatics:
%[text] AUC > 0.9  Excellent
%[text] 0.7-0.9    Practically good model
%[text] 0.6-0.7    Barely useful
%[text] < 0.6      Poor (check for bugs/label noise)
%[text] 
%[text] AUC is threshold-independent and measures overall discrimination ability.
%[text] It is preferred over accuracy when classes are imbalanced.
%[text] 
%[text] Note on Label Noise:
%[text] If some household chemicals are true drugs (mislabelled as class 0),
%[text] AUC = 1.0 is impossible. The achievable maximum AUC is limited by label quality.
%[text] 
%[text] **Appearance of ROC with Small Samples**
%[text] The test set has only 6 non-drug samples. Therefore, the FPR axis can only take
%[text] 7 steps: 0, 1/6≈0.167, 2/6≈0.333, ..., making the ROC a coarse **step function** rather than a smooth S-curve seen in textbooks.
%[text] This is not an implementation or calculation error but normal behavior for small test sets.
%[text] The curve always starts at (0,0), but if drug scores are higher than non-drug,
%[text] a vertical segment (hard to see) appears at the left with FPR=0.
%[text] Use cross-validation (💡 Observation Point 6) for a smoother ROC.
%[text] 
%[text] SVM: Column 2 of scores = decision values for class 1 (higher means drug)
scoreSVM = scores_svm(:, 2);
%[text] RF: Column 2 = posterior probability of class 1 (drug)
%[text] ClassNames=[0;1] so column 1 = P(class 0), column 2 = P(class 1)
scoreRF  = scores_rf(:, 2);

[xSVM, ySVM, ~, aucSVM] = perfcurve(yTest, scoreSVM, 1);
[xRF,  yRF,  ~, aucRF]  = perfcurve(yTest, scoreRF,  1);

figure("Name", "A04 ROC Curve");
hold on;
plot(xSVM, ySVM, Color=[0.3 0.6 0.9], LineWidth=2, ...
    DisplayName=sprintf("SVM Linear (AUC=%.3f)", aucSVM));
plot(xRF,  yRF,  Color=[0.4 0.7 0.4], LineWidth=2, ...
    DisplayName=sprintf("Random Forest (AUC=%.3f)", aucRF));
plot([0 1], [0 1], "--k", LineWidth=1, DisplayName="Random Baseline (AUC=0.500)");
hold off;
xlabel("False Positive Rate (1 - Specificity)");
ylabel("True Positive Rate (Sensitivity/Recall)");
title("ROC Curve -- Drug vs Non-drug Classification");
legend(Location="southeast");
grid on; axis square;   % Square aspect; unlike axis equal, keeps xlim/ylim at [0,1]
xlim([0 1]); ylim([0 1]);
%[text] Summary Table
summaryTbl = table( ...
    ["SVM (Linear)"; "Random Forest"], ...
    [acc_svm;  acc_rf],  ...
    [prec_svm; prec_rf], ...
    [rec_svm;  rec_rf],  ...
    [f1_svm;   f1_rf],   ...
    [aucSVM;   aucRF],   ...
    VariableNames=["Model","Accuracy","Precision","Recall","F1","AUC"]);
logInfo("Performance Summary:");
disp(summaryTbl);

logInfo("ROC AUC: SVM=%.3f  RF=%.3f", aucSVM, aucRF);
%[text] **💡 Observation Point 6**
%[text] Check which model has a higher AUC. Also, verify if both have AUC > 0.5.
%[text] AUC <= 0.5 means worse than random, a red flag.
%[text] Find the "optimal threshold" (point closest to the top-left corner of the ROC).
%[text] [~, optIdx] = min((1-ySVM).^2 + xSVM.^2);   % Distance to (FPR=0, TPR=1)
%[text] Calculate sensitivity and specificity at this threshold and check if it differs from the default 0.5.
%[text] Is this threshold different from the default 0.5?
%[text] Perform a null model check (shuffle labels and rerun).
%[text] shuffleIdx = randperm(numel(yTest));
%[text] [~, ~, ~, aucNull] = perfcurve(yTest, scoreSVM(shuffleIdx), 1);
%[text] The AUC of shuffled labels should be around 0.5. Ensure the pipeline is
%[text] correctly implemented and there is no data leakage.
%[text] Perform 5-fold cross-validation and calculate the mean AUC.
%[text] mdlCV = crossval(svmModel, KFold=5);
%[text] [yCV, scoresCV] = kfoldPredict(mdlCV);
%[text] [~,~,~,aucCV] = perfcurve(yTrain, scoresCV(:,2), 1);
%[text] Compare CV AUC with single test set AUC.
%[text] (For small datasets, CV AUC is a more reliable estimate)
%[text] Label noise from household drugs (Aspirin, Caffeine) limits maximum AUC.
%[text] Remove these and rerun. Check if AUC improves and by how much.
%[text] **Summary**
%[text] - p>>n 2048 bits 1 230 molecules require PCA before SVM.
%[text] - Class imbalance makes accuracy misleading — use AUC and F1.
%[text] - RF's Prior=uniform and SVM's cost weighting address class imbalance.
%[text] - AUC is a threshold-independent measure of discrimination ability.
%[text] - Label quality limits achievable maximum AUC.
%[text] - Always validate with a null model (shuffled labels). \

%[appendix]{"version":"1.0"}
%[metadata:view]

%---

%   data: {"layout":"inline","rightPanelPercent":40}
%---
