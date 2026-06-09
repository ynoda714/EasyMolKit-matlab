%[text] # A03: QSAR Regression — Predicting LogP from Molecular Descriptors
%[text] EasyMolKit Analytics — Layer 3
%[text] 
%[text] LogP (the logarithm of the octanol-water partition coefficient) is a crucial physicochemical property that determines whether a drug can pass through lipid membranes and be absorbed into the body.
%[text] Experimental measurement requires compound synthesis and liquid-liquid partition assays, which are costly and time-consuming.
%[text] If LogP can be predicted using only 8 types of structural descriptors representing the "shape" of the molecule, thousands of candidates can be virtually screened before synthesis.
%[text] This script builds two models, linear regression (`fitlm`) and random forest (`fitrensemble`), and fairly compares their performance using 5-fold cross-validation.
%[text] 
%[text] **What you will learn in this tutorial**
%[text] - Understand the importance of QSAR regression in early drug discovery
%[text] - Build and interpret a linear regression model using `fitlm()`
%[text] - Perform nonlinear ensemble (random forest) regression using `fitrensemble()`
%[text] - Apply k-fold cross-validation to estimate generalization performance
%[text] - Compare models using R², RMSE, and observed vs predicted plots
%[text] - Understand the concept of the applicability domain \
%[text] 
%[text] **Prerequisites**
%[text] - Completion of F02 (Property Calculation) — Basics of Descriptors
%[text] - Recommended: A01 (Chemical Space PCA) to understand the context of descriptors
%[text] - Statistics and Machine Learning Toolbox (`fitlm`, `fitrensemble`, `cvpartition`)
%[text] - No internet connection required \
%[text] 
%[text] Estimated time required: 30–45 minutes
%[text] 
%[text] **Data:**
%[text] `data/list/fda_drugs.csv` — 200 FDA-approved drugs (ChEMBL, CC-BY-SA 3.0)
%[text] ALogP column: Crippen-Wildman estimates stored in ChEMBL. Used as the regression target variable in this exercise (treated as the "ground truth").
%[text] 
%[text] **References**
%[text] - Wildman SA & Crippen GM (1999) Prediction of physicochemical parameters by atomic contributions. J Chem Inf Comput Sci 39:868-873.  doi:10.1021/ci990307l (Crippen-Wildman LogP / ALogP) \
%[text] - Breiman L (2001) Random forests. Machine Learning 45:5-32. doi:10.1023/A:1010933404324 \
%[text] - Cherkasov A et al. (2014) QSAR modeling: where have you been?  J Med Chem 57:4977-5010. doi:10.1021/jm4004285 \
%[text] - Oprea TI (2000) Property distribution of drug-related chemical databases.  J Comput Aided Mol Des 14:251-264. doi:10.1023/A:1008130001697 \
%[text] 
%[text] How to run: Execute each section with Ctrl+Enter
%%
%[text] ## Section 0: Setup
logSection("A03", "Section 0: Setup", "Analytics L3");
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
%[text] Prepare Python and RDKit processes before main execution
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;
%%
%[text] ## Section 1: Loading FDA-approved drugs and calculating descriptors
%[text] 
%[text] Setup is complete. First, prepare the data and features for prediction.
%[text] Calculate descriptors from SMILES and create 8 types of feature matrices X and target variable y.
%[text] 
%[text] ### Concept: Feature Selection — Why exclude LogP?
%[text] The target variable is ALogP (ChEMBL's Crippen-Wildman estimate).
%[text] Since RDKit's built-in LogP descriptor also uses the same Crippen-Wildman method, including it as a feature would result in a trivial model that "reads LogP to predict LogP," which is not meaningful for QSAR.
%[text] 
%[text] Instead, use 8 complementary descriptors that encode molecular structure without directly representing lipophilicity:
%[text] - Size: MolWt (g/mol), HeavyAtomCount (non-hydrogen atom count)
%[text] - Polarity: TPSA (Å²), NumHDonors, NumHAcceptors
%[text] - Topology: NumRotatableBonds, RingCount
%[text] - Shape: FractionCSP3 (fraction of sp3 carbons) \
%[text] Chemical rationale for LogP prediction using these descriptors:
%[text] - High TPSA / HBD / HBA -> More polar -> Lower LogP (increased solubility in water)
%[text] - More RingCount -> Carbon aromatic rings contribute to higher LogP, but in FDA-approved drugs with many polar heterocycles, the direction can reverse (check model coefficients)
%[text] - High FractionCSP3 -> More 3D shape -> Reduced planarity -> Lower LogP
%[text] - Larger MolWt has a weak positive correlation with LogP within the drug discovery space \
DATA_FILE  = "data/list/fda_drugs.csv";
FEAT_NAMES = ["MolWt", "TPSA", "NumHDonors", "NumHAcceptors", ...
              "NumRotatableBonds", "RingCount", "FractionCSP3", "HeavyAtomCount"];
N_FEATS    = numel(FEAT_NAMES);

rawTbl = readtable(DATA_FILE, TextType="string");
logInfo("Loaded %d molecules from %s", height(rawTbl), DATA_FILE);
%[text] Retrieve ALogP column (may be read as numeric or string by the parser)
if isnumeric(rawTbl.ALogP)
    alogpVec = rawTbl.ALogP;
else
    alogpVec = str2double(rawTbl.ALogP);
end
%[text] Parse SMILES and calculate descriptors
nRaw  = height(rawTbl);
X_all = nan(nRaw, N_FEATS);
y_all = nan(nRaw, 1);
valid = false(nRaw, 1);

for k = 1:nRaw
    if isnan(alogpVec(k)), continue; end      % Skip if ALogP is missing
    try
        mol = emk.mol.fromSmiles(rawTbl.SMILES(k));
        s   = emk.descriptor.calculate(mol, FEAT_NAMES);
        row = zeros(1, N_FEATS);
        for d = 1:N_FEATS
            row(d) = s.(FEAT_NAMES(d));
        end
        X_all(k, :) = row;
        y_all(k)    = alogpVec(k);
        valid(k)    = true;
    catch ME
        logWarn("Skipping %s: %s", rawTbl.Name(k), ME.message);
    end
end
logInfo("Completed descriptor calculation for %d / %d molecules", sum(valid), nRaw);

X        = X_all(valid, :);    % N×8 feature matrix (raw data, non-standardized)
y        = y_all(valid);       % N×1 target variable vector (ALogP)
molNames = rawTbl.Name(valid);
nMols    = sum(valid);

logInfo("Dataset: Analyzed %d / %d molecules", nMols, nRaw);
logInfo("ALogP range: %.2f to %.2f (mean=%.2f, std=%.2f)", ...
    min(y), max(y), mean(y), std(y));
%[text] Check the first 5 rows in table format
prevN = min(5, nMols);
prevTbl = array2table(X(1:prevN, :), VariableNames=cellstr(FEAT_NAMES));
prevTbl.ALogP = y(1:prevN);
prevTbl.Name  = molNames(1:prevN);
prevTbl = movevars(prevTbl, "Name", Before=1);
disp(prevTbl);
%[text] **💡 Observation Point 1**
%[text] Let's identify the descriptor with the strongest Pearson correlation with ALogP.
%[text] Use corr(X, y) to check the 8x1 correlation vector.
%[text] Interpret whether the sign of each correlation is chemically reasonable.
%[text] (Hint: TPSA has a negative correlation with LogP -- high polarity = low lipophilicity)
%[text] Calculate and check pairwise correlations among the 8 features: corr(X)
%[text] Identify the pair with the highest multicollinearity (prediction: MolWt and HeavyAtomCount).
%[text] Consider how strong multicollinearity affects linear regression coefficients.
% ... (Try writing code here)
%%
%[text] ## Section 2: Exploring LogP Distribution
%[text] 
%[text] Descriptor calculation is complete. Next, let's visualize the distribution of the target variable, ALogP, before modeling.
%[text] Pre-checking the data helps in accurately interpreting the subsequent model results.
%[text] 
%[text] ### Concept: Examining the Target Variable Before Modeling
%[text] Linear regression assumes that residuals are approximately normally distributed.
%[text] (y itself does not need to be normally distributed, but it is related).
%[text] What you can learn by examining y first:
%[text] 1. Shape: Skewed distributions may be improved with log or Box-Cox transformations.
%[text] 2. Outliers: Extreme y values can significantly distort OLS fit.
%[text] 3. Range: The model's applicable domain is limited by the training range. \
%[text] Extrapolation beyond this range is less reliable.
%[text] 
%[text] Considering drug-like LogP (Lipinski's Rule of 5: LogP ≤ 5), the distribution of FDA-approved drugs peaks in the range of -2 to +5, with a tail towards lipophilic outliers (including poorly soluble drugs).
figure("Name", "A03 LogP Distribution");
tiledlayout(1, 2);

nexttile;
histogram(y, 20, FaceColor=[0.3 0.6 0.9]);
xlabel("ALogP"); ylabel("Frequency");
title("Distribution of ALogP (FDA-approved drugs)");
xline(mean(y),   "--r", sprintf("mean=%.1f",   mean(y)),   LabelHorizontalAlignment="right");
xline(median(y), "--g", sprintf("median=%.1f", median(y)), LabelHorizontalAlignment="left");
grid on;

nexttile;
scatter(y, X(:, strcmp(cellstr(FEAT_NAMES), "TPSA")), 30, [0.3 0.6 0.9], ...
    "filled", MarkerFaceAlpha=0.6);
xlabel("ALogP"); ylabel("TPSA (A^2)");
title(sprintf("TPSA vs ALogP  r=%.2f", corr(y, X(:, strcmp(cellstr(FEAT_NAMES), "TPSA")))));
grid on;

logInfo("ALogP: Mean=%.2f, Median=%.2f, Std Dev=%.2f, Skewness=%.2f", ...
    mean(y), median(y), std(y), skewness(y));
%[text] Scatter Plot of Features and Target Variable (2×4 Grid)
figure("Name", "A03 Feature-LogP Scatter Plot");
for f = 1:N_FEATS
    subplot(2, 4, f);
    scatter(X(:, f), y, 20, [0.2 0.5 0.8], "filled", MarkerFaceAlpha=0.5);
    xlabel(FEAT_NAMES(f)); ylabel("ALogP");
    r = corr(X(:, f), y);
    title(sprintf("r=%.2f", r));
    grid on;
end
sgtitle("8 Descriptors vs ALogP");
%[text] **💡 Observation Point 2**
%[text] Let's check if there are molecules with ALogP \> 6.
%[text] \[sorted\_y, ord\] = sort(y, "descend");
%[text] table(molNames(ord(1:5)), sorted\_y(1:5), VariableNames=\["Name","ALogP"\])
%[text] Consider whether these outliers might adversely affect the linear model.
%[text] The overall LogP distribution of FDA-approved drugs is said to peak around 2.5 (Lipinski space).
%[text] Let's check if this dataset meets that expectation.
%[text] Look for descriptors that show a clear nonlinear relationship.
%[text] (Nonlinear patterns may indicate the limits of linear regression)
% ... (Try writing code here)
%%
%[text] ## Section 3: Linear Regression (fitlm)
%[text] 
%[text] We have confirmed the distribution characteristics. Here, we begin analysis using the most basic model, linear regression.
%[text] Linear models are easy to interpret in terms of coefficients and are important as a baseline.
%[text] 
%[text] ### Concept: Multiple Regression Analysis for QSAR
%[text] Linear QSAR models predict LogP (the logarithm of the partition coefficient) as a weighted sum of descriptors:
%[text] 
%[text] LogP\_pred = b0 + b1*MolWt + b2*TPSA + ... + b8\*HeavyAtomCount
%[text] 
%[text] Coefficients b1...b8 are estimated by ordinary least squares (OLS):
%[text] Minimize: sum\_i (y\_i - X\_i \* b)^2
%[text] 
%[text] Standardizing features (z-score) before fitting allows for comparison:
%[text] - Large |values| in standardized data coefficients indicate that the descriptor has a strong influence on LogP in this dataset. \
%[text] - Coefficients without standardization are on different scales (e.g., b\_MolWt units are LogP / g mol^-1). \
%[text] 
%[text] Key diagnostics of fitlm:
%[text] - R2      -- Proportion of variance explained (training set; optimistic)
%[text] - RMSE    -- Root mean square error (LogP units)
%[text] - p-value -- Whether this coefficient is significantly different from zero
%[text] - Residuals -- Should be normally distributed without patterns against fitted values \
%[text] 
%[text] Standardize features (normalize to mean 0, variance 1)
Xmean = mean(X, 1);
Xstd  = std(X,  0, 1);
Xstd(Xstd < 1e-12) = 1;      % Prevent division by zero for constant descriptors
Xs    = (X - Xmean) ./ Xstd;  % N×8 standardized feature matrix
%[text] Fit a linear model with all training data
ftTbl = array2table(Xs, VariableNames=cellstr(FEAT_NAMES));
ftTbl.ALogP = y;
lmModel = fitlm(ftTbl, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));

logInfo("Linear model (training): R2=%.3f  RMSE=%.3f LogP units", ...
    lmModel.Rsquared.Ordinary, lmModel.RMSE);
disp(lmModel);
%[text] Plot and check the observed vs. predicted values and residuals
yPred_lm  = predict(lmModel, ftTbl);
residuals = y - yPred_lm;

figure("Name", "A03 Linear Model Diagnostics");
set(gcf, "Position", [100 100 1100 380]);
tiledlayout(1, 3);

nexttile;
scatter(y, yPred_lm, 30, [0.3 0.6 0.9], "filled", MarkerFaceAlpha=0.7);
hold on;
refLim = [min(y)-0.5, max(y)+0.5];
plot(refLim, refLim, "--k", LineWidth=1.5);
hold off;
xlabel("Observed ALogP"); ylabel("Predicted ALogP");
title(sprintf("Linear Model\nR^2=%.3f  RMSE=%.3f", lmModel.Rsquared.Ordinary, lmModel.RMSE));
grid on;

nexttile;
histogram(residuals, 20, FaceColor=[0.9 0.5 0.3]);
xlabel("Residuals (Observed - Predicted)"); ylabel("Frequency");
title("Distribution of Residuals");
xline(0, "--k"); grid on;

nexttile;
[sortCoeff, sortOrd] = sort(lmModel.Coefficients.Estimate(2:end), "descend");
barh(fliplr(sortCoeff), FaceColor=[0.3 0.6 0.9]);
yticks(1:N_FEATS);
yticklabels(fliplr(cellstr(FEAT_NAMES(sortOrd))));
xlabel("Standardized Coefficients"); xline(0, "k");
title("Regression Coefficients (Standardized)"); grid on;

logInfo("Descriptor with the strongest positive influence: %s (Coefficient=%.3f)", ...
    FEAT_NAMES(sortOrd(1)), sortCoeff(1));
logInfo("Descriptor with the strongest negative influence: %s (Coefficient=%.3f)", ...
    FEAT_NAMES(sortOrd(end)), sortCoeff(end));
%[text] **💡 Observation Point 3**
%[text] Check which standardized coefficient has the largest positive descriptor.
%[text] Check which has the largest negative coefficient.
%[text] Read whether these signs match the chemical expectations from Section 1.
%[text] Check if there are features with p-value \> 0.05 (not statistically significant).
%[text] Read whether removing them and refitting improves RMSE.
%[text] For example, remove HeavyAtomCount (multicollinearity with MolWt):
%[text] ftTbl2 = ftTbl; ftTbl2.HeavyAtomCount = \[\]; 
%[text] lm2 = fitlm(ftTbl2, "ALogP ~ MolWt + TPSA + NumHDonors + ...");
%[text] Check if the distribution of residuals appears approximately normal.
%[text] A right-skewed residual distribution indicates potential outliers.
%[text] Check which molecule has the largest |residual|.
%[text] \[~, bigResOrd\] = sort(abs(residuals), "descend");  molNames(bigResOrd(1:5))
% ... (Try writing code here)
%%
%[text] ## Section 4: Random Forest Regression (fitrensemble)
%[text] 
%[text] We have constructed a linear model. Next, let's try Random Forest, which can capture nonlinear relationships. logSection("A03", "Section 4: Random Forest Regression (fitrensemble)", "Analytics L3");%\[text\] In the later section, we will fairly compare the two models through cross-validation.
%[text] 
%[text] ### Concept: Random Forest for QSAR
%[text] Random Forest (Breiman 2001) is an ensemble of B decision trees:
%[text] 1. Extract B bootstrap samples from the training data (with replacement).
%[text] 2. Fit a decision tree to each sample, considering only m features randomly at each split. \
%[text] (mtry = floor(p/3) by default, reducing correlation between trees)
%[text] 1. Predictions are made by averaging the outputs of the B trees. \
%[text] 
%[text] Advantages over linear regression:
%[text] - Can capture nonlinear relationships (tree splits are nonlinear).
%[text] - Can handle interactions between features (splits can depend on two features).
%[text] - Provides built-in feature importance via decrease in impurity (MDI).
%[text] - Robust to outliers (trees use rank-based splits). \
%[text] 
%[text] Drawbacks: Training R2 is very optimistic (trees can memorize data). It is crucial to always evaluate with cross-validation (Section 5).
%[text] 
%[text] Feature importance by MDI (Mean Decrease in Impurity):
%[text] Calculate the sum of variance reduction × sample proportion at each node for all trees and splits using feature j. Larger values indicate greater usefulness.
rfModel = fitrensemble(Xs, y, ...
    Method="Bag", ...
    NumLearningCycles=200, ...
    Learners=templateTree(MinLeafSize=3));

yPred_rfTrain  = predict(rfModel, Xs);
ss_res_rf      = sum((y - yPred_rfTrain).^2);
ss_tot         = sum((y - mean(y)).^2);
r2_rf_train    = 1 - ss_res_rf / ss_tot;
rmse_rf_train  = sqrt(mean((y - yPred_rfTrain).^2));
logInfo("Random Forest (Training, Optimistic): R2=%.3f  RMSE=%.3f", ...
    r2_rf_train, rmse_rf_train);
logInfo("Note: RF training metrics are optimistic -- check with CV in Section 5");
%[text] Visualization of feature importance
importances = predictorImportance(rfModel);
[sortedImp, sortOrdRF] = sort(importances, "descend");

refLim = [min(y)-0.5, max(y)+0.5];   % Redefine to run independently from Section 3
figure("Name", "A03 Random Forest");
tiledlayout(1, 2);

nexttile;
barh(fliplr(sortedImp), FaceColor=[0.4 0.7 0.4]);
yticks(1:N_FEATS);
yticklabels(fliplr(cellstr(FEAT_NAMES(sortOrdRF))));
xlabel("Feature Importance (MDI)");
title("RF Feature Importance");
grid on;

nexttile;
scatter(y, yPred_rfTrain, 30, [0.4 0.7 0.4], "filled", MarkerFaceAlpha=0.7);
hold on;
plot(refLim, refLim, "--k", LineWidth=1.5);
hold off;
xlabel("Observed ALogP"); ylabel("Predicted (RF, Training)");
title(sprintf("RF Training  R^2=%.3f  RMSE=%.3f", r2_rf_train, rmse_rf_train));
grid on;

logInfo("RF Most Important Feature: %s  Least Important: %s", ...
    FEAT_NAMES(sortOrdRF(1)), FEAT_NAMES(sortOrdRF(end)));
%[text] **💡 Observation Point 4**
%[text] Check the features that Random Forest evaluates as most important.
%[text] See if they match the largest coefficients in the linear model.
%[text] (RF captures nonlinear contributions, so they may differ.)
%[text] Check changes when trying MinLeafSize=1 and MinLeafSize=20.
%[text] Consider how training RMSE changes and the risks of MinLeafSize=1.
%[text] (Smaller leaves = deeper trees = increased risk of overfitting to training data.)
%[text] Check if increasing NumLearningCycles to 500 improves training RMSE.
%[text] Consider how many trees stabilize the model (point of diminishing returns).
%[text] oobLoss(rfModel) returns OOB RMSE as a function of the number of trees.
%[text] Try it: plot(oobLoss(rfModel)); xlabel("Number of Trees"); ylabel("OOB RMSE");
% ... (Try writing code here)
%%
%[text] ## Section 5: Comparison by 5-Fold Cross-Validation
%[text] 
%[text] We have two models ready. This is the key point of this session.
%[text] Errors evaluated on training data tend to be optimistic, so we use cross-validation to fairly estimate generalization performance.
%[text] 
%[text] ### Concept: Cross-Validation for Unbiased Model Comparison
%[text] Training error is always optimistic because the model has already seen the data.
%[text] k-fold cross-validation (CV) is a method to estimate generalization error.
%[text] 1. Divide n molecules into k equal folds.
%[text] 2. For each fold i, do the following: \
%[text] - Train on data excluding fold i.
%[text] - Predict on fold i (unused during training).
%[text] - Record residuals.
%[text] - Aggregate out-of-fold predictions to compute CV metrics. \
%[text] 
%[text] k=5 is standard, training on 80% of the data per fold. For small datasets (n<100), use k=10 or LOO (leave-one-out).
%[text] 
%[text] Evaluation metrics are as follows:
%[text] -  - R2 = 1 - SS\_res/SS\_tot (1=perfect, 0=equivalent to mean prediction, negative=below mean)
%[text] -  - RMSE = sqrt(mean((y\_true - y\_pred)^2))  \[LogP units\]
%[text] -  - MAE = mean(|y\_true - y\_pred|) (less sensitive to outliers than RMSE) \
%[text] 
%[text] A good QSAR regression model typically achieves CV R2 > 0.6. CV R2 > 0.8 is considered excellent for drug-like property prediction.
K_FOLDS = 5;
rng(42);
cv = cvpartition(nMols, KFold=K_FOLDS);

yCV_lm = nan(nMols, 1);   % out-of-fold predictions: linear model
yCV_rf = nan(nMols, 1);   % out-of-fold predictions: random forest

for fold = 1:K_FOLDS
    trainIdx = training(cv, fold);
    testIdx  = test(cv, fold);

    % Linear model for this fold
    trTbl = array2table(Xs(trainIdx, :), VariableNames=cellstr(FEAT_NAMES));
    trTbl.ALogP = y(trainIdx);
    lm_f  = fitlm(trTbl, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));
    teTbl = array2table(Xs(testIdx, :), VariableNames=cellstr(FEAT_NAMES));
    yCV_lm(testIdx) = predict(lm_f, teTbl);

    % Random Forest for this fold
    rf_f = fitrensemble(Xs(trainIdx, :), y(trainIdx), ...
        Method="Bag", NumLearningCycles=200, ...
        Learners=templateTree(MinLeafSize=3));
    yCV_rf(testIdx) = predict(rf_f, Xs(testIdx, :));
end
logInfo("%d-Fold CV completed", K_FOLDS);
%[text] Calculate CV metrics (generate structure with anonymous function).
computeMetrics = @(yt, yp) struct( ...
    "R2",   1 - sum((yt-yp).^2) / sum((yt-mean(yt)).^2), ...
    "RMSE", sqrt(mean((yt-yp).^2)), ...
    "MAE",  mean(abs(yt-yp)));

cv_lm = computeMetrics(y, yCV_lm);
cv_rf = computeMetrics(y, yCV_rf);

logInfo("--- %d-Fold CV Results ---", K_FOLDS);
logInfo("Linear Model:     R2=%.3f  RMSE=%.3f  MAE=%.3f", ...
    cv_lm.R2, cv_lm.RMSE, cv_lm.MAE);
logInfo("Random Forest: R2=%.3f  RMSE=%.3f  MAE=%.3f", ...
    cv_rf.R2, cv_rf.RMSE, cv_rf.MAE);
logInfo("Overfitting Gap (RF): Training RMSE=%.3f  CV RMSE=%.3f", ...
    rmse_rf_train, cv_rf.RMSE);
%[text] ### Why did linear regression perform better this time?
%[text] The target variable ALogP in this dataset is an estimate calculated by the Crippen-Wildman method (a **linear sum** of atomic contributions). Since the generation mechanism is linear, linear regression fits well with the data structure, and the flexibility of RF led to overfitting in the small dataset with n=200.
%[text] When using experimentally measured LogP (including nonlinear intermolecular interactions) as the target variable, RF may outperform. This is an important lesson that "the strongest model does not always win." Always use CV for model selection.
%[text] Plot comparison of CV actual and predicted values
figure("Name", "A03 Cross-Validation Comparison");
tiledlayout(1, 2);

nexttile;
scatter(y, yCV_lm, 30, [0.3 0.6 0.9], "filled", MarkerFaceAlpha=0.7);
hold on; plot(refLim, refLim, "--k", LineWidth=1.5); hold off;
xlabel("Actual ALogP"); ylabel("CV Predicted ALogP");
title(sprintf("Linear Model (5-Fold CV)\nR^2=%.3f  RMSE=%.3f", cv_lm.R2, cv_lm.RMSE));
grid on;

nexttile;
scatter(y, yCV_rf, 30, [0.4 0.7 0.4], "filled", MarkerFaceAlpha=0.7);
hold on; plot(refLim, refLim, "--k", LineWidth=1.5); hold off;
xlabel("Actual ALogP"); ylabel("CV Predicted ALogP");
title(sprintf("Random Forest (5-Fold CV)\nR^2=%.3f  RMSE=%.3f", cv_rf.R2, cv_rf.RMSE));
grid on;
%[text] **💡 Observation Point 5**
%[text] Check which model is superior based on CV R2 and RMSE. Check if the R2 difference exceeds 0.1 (practically significant).
%[text] Compare training RMSE (Section 4) and CV RMSE (here) to check the "overfitting gap" = CV RMSE - training RMSE.
%[text] Consider which model has a larger gap and what it implies.
%[text] Check how CV results change if K\_FOLDS is changed to 10.
%[text] For n=200, the bias-variance trade-off differs between 5-fold vs 10-fold CV.
%[text] (More folds reduce bias but increase the variance of CV estimates.)
%[text] Try tuning RF hyperparameters,
%[text] rfTuned = fitrensemble(Xs, y, "OptimizeHyperparameters", "auto", ...
%[text] HyperparameterOptimizationOptions=struct(MaxObjectiveEvaluations=20));
%[text] Check if Bayesian optimization improves CV R2 over default settings.
% ... (Try writing code here)
%%
%[text] ## Section 6: Applying the Model to New Molecules
%[text] 
%[text] We have selected a superior model through cross-validation. Finally, let's apply this model in practice.
%[text] We will predict LogP (the logarithm of the partition coefficient) for new molecules not included in the training data and check the applicability domain.
%[text] 
%[text] ### Concept: Prediction Pipeline and Applicability Domain
%[text] The QSAR model trained on FDA-approved drugs is most reliable for molecules similar to the training set (within the applicability domain).
%[text] 
%[text] Applicability Domain Check (Simple Range Method):
%[text] - Standardize the descriptors of new molecules using the mean and standard deviation of the training set.
%[text] - If the standardized features are more than 2 sigma outside the training range, \
%[text] flag the prediction as "out-of-domain" (less reliable).
%[text] More rigorous domain methods exist (leverage-based Euclidean distance or k-NN distance), but range checks are effective as a practical initial screening.
%[text] 
%[text] Practical Note:
%[text] Before applying to new molecules, retrain the selected model with all data. (The more data, the better the model. CV is only for model selection).
%[text] 
%[text] Select the model based on CV R2.
if cv_rf.R2 >= cv_lm.R2
    finalModel = fitrensemble(Xs, y, Method="Bag", NumLearningCycles=200, ...
        Learners=templateTree(MinLeafSize=3));
    modelLabel = "Random Forest";
    predictFn  = @(xS) predict(finalModel, xS);
else
    finalTbl   = array2table(Xs, VariableNames=cellstr(FEAT_NAMES));
    finalTbl.ALogP = y;
    finalModel = fitlm(finalTbl, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));
    modelLabel = "Linear";
    predictFn  = @(xS) predict(finalModel, ...
        array2table(xS, VariableNames=cellstr(FEAT_NAMES)));
end
logInfo("Final model: %s (retrained with all %d molecules)", modelLabel, nMols);
%[text] New molecules to predict
newSmiles = ["CCO", ...                              % Ethanol
             "c1ccccc1", ...                         % Benzene
             "CC(=O)Oc1ccccc1C(=O)O", ...           % Aspirin
             "CN1C=NC2=C1C(=O)N(C)C(=O)N2C", ...   % Caffeine
             "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34"]; % Testosterone

logInfo("Predicting LogP for %d new molecules:", numel(newSmiles));
for k = 1:numel(newSmiles)
    try
        mol = emk.mol.fromSmiles(newSmiles(k));
        s   = emk.descriptor.calculate(mol, FEAT_NAMES);
        xNew = zeros(1, N_FEATS);
        for d = 1:N_FEATS
            xNew(d) = s.(FEAT_NAMES(d));
        end
        xNewS = (xNew - Xmean) ./ Xstd;   % Standardize with training statistics

        % Applicability Domain: 2σ Range Check (after standardization with training statistics, all features within |z| <= 2)
        inDomain = all(abs(xNewS) <= 2);
        domainStr = "OK";
        if ~inDomain, domainStr = "OUT-OF-DOMAIN"; end

        logpPred = predictFn(xNewS);
        logInfo("  %-46s -> LogP_pred=%5.2f  [%s]", newSmiles(k), logpPred, domainStr);
    catch ME
        logWarn("  Failed for %s: %s", newSmiles(k), ME.message);
    end
end
%[text] **💡 Observation Point 6**
%[text] Compare the predicted LogP with RDKit's built-in Crippen estimate.
%[text] mol = `emk.mol.fromSmiles("CCO")`;
%[text] s = `emk.descriptor.calculate(mol, "LogP")`;  s.LogP
%[text] Check how close the QSAR prediction is to the Crippen value.
%[text] ALogP (the learning target variable) is also calculated by the Crippen method, so they should be similar.
%[text] In other words, this can be considered a "meta-model" of Crippen.
%[text] Try it with very large molecules as well (e.g., peptide Ala-Gly: "NCC(=O)NCC(=O)O").
%[text] Check if the applicability domain triggers "out-of-domain".
%[text] Consider what happens if you use the prediction anyway.
%[text] Think about why the final model should be retrained with all data before deployment.
%[text] (Hint: CV trains with 80% of the data per fold. Using all data yields a better model. The 20% holdout is only for fair evaluation).
%[text] **Summary**
%[text] - Standardize features to z-scores before linear regression to obtain interpretable coefficients.
%[text] - Training error is always optimistic. Always estimate generalization performance with cross-validation.
%[text] - Random Forest captures non-linear effects that linear models miss. However, linear models may outperform with small data.
%[text] - Always check the applicability domain before trusting QSAR predictions. \

%[appendix]{"version":"1.0"}
%[metadata:view]%---

%   data: {"layout":"inline","rightPanelPercent":40}
%---
