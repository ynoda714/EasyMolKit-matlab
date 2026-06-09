%[text] # A05: Molecular Property Prediction Using Neural Networks
%[text] EasyMolKit Analytics — Layer 3
%[text] 
%[text] QSAR (Quantitative Structure-Activity Relationship) models have traditionally relied on linear regression. If the relationship between descriptors and properties is nearly linear, linear models are fast, interpretable, and effective.
%[text] However, flexibility comes at a cost.
%[text] Neural networks offer a more flexible alternative, capturing nonlinear interactions between descriptors without requiring chemists to specify which interactions are important.
%[text] On the other hand, neural networks require more data, are difficult to interpret, and can severely overfit on the typically small datasets in drug discovery.
%[text] In this exercise, we will build a small feedforward neural network to predict ALogP (Crippen-Wildman lipophilicity) from eight structural descriptors of FDA-approved drugs and directly compare it with the linear regression model from A03.
%[text] 
%[text] **What You Will Learn in This Tutorial**
%[text] - Understand why small QSAR datasets are challenging for neural networks
%[text] - Build and train a feedforward regression network using trainnet()
%[text] - Reduce overfitting with dropout and learning rate schedules
%[text] - Interpret learning curves: convergence, overfitting, underfitting
%[text] - Know when to choose simpler linear models over neural networks \
%[text] 
%[text] **Prerequisites**
%[text] - Recommended: A03 (QSAR Regression) for context on the same dataset
%[text] - Deep Learning Toolbox (trainnet, trainingOptions, fullyConnectedLayer)
%[text] - No internet connection required \
%[text] 
%[text] Estimated Time Required: 35–50 minutes | How to Run: Execute sections one by one with Ctrl+Enter
%[text] 
%[text] **Data:** data/list/fda\_drugs.csv — 200 FDA-approved drugs (ChEMBL, CC-BY-SA 3.0). The ALogP column is used as the regression target for Crippen-Wildman lipophilicity estimates.
%[text] 
%[text] **References**
%[text] - Srivastava N, Hinton G, Krizhevsky A, Sutskever I, Salakhutdinov R (2014) Dropout: a simple way to prevent neural networks from overfitting. *J Mach Learn Res* 15:1929-1958.
%[text] - Kingma D & Ba J (2015) Adam: A Method for Stochastic Optimization. arXiv:1412.6980.
%[text] - Ramsundar B et al. (2015) Massively multitask networks for drug discovery. arXiv:1502.02072.
%[text] - Tropsha A (2010) Best practices for QSAR model development, validation, and exploitation. *Mol Inform* 29:476-488. doi:10.1002/minf.201000061. \
%%
%[text] ## Section 0: Setup
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
%[text] Warm up the Python/RDKit process before main execution
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;
logSection("A05", "Section 0: Setup", "Analytics L3");
%%
%[text] ## Section 1: Loading FDA-approved Drugs and Calculating Descriptors
%[text] 
%[text] Setup is complete. First, we will calculate the 8 types of descriptors used in A03 to prepare the dataset.
%[text] By executing this section, we ensure that A05 operates independently.
%[text] 
%[text] ### Concept: Reusing the Descriptor Matrix from A03
%[text] The 8 descriptors used in A03 (MolWt, TPSA, NumHDonors, NumHAcceptors, NumRotatableBonds, RingCount, FractionCSP3, HeavyAtomCount) encode molecular structures in four chemical dimensions: size, polarity, topology, and shape. The target variable is ALogP, which is the Crippen-Wildman lipophilicity estimate stored in ChEMBL.
%[text] 
%[text] If you have already executed A03, you can reuse the workspace data, but executing this section ensures that A05 is fully self-contained.
logSection("A05", "Section 1: Loading FDA-approved Drugs and Calculating Descriptors", "Analytics L3");
DATA_FILE  = "data/list/fda_drugs.csv";
FEAT_NAMES = ["MolWt", "TPSA", "NumHDonors", "NumHAcceptors", ...
              "NumRotatableBonds", "RingCount", "FractionCSP3", "HeavyAtomCount"];
N_FEATS    = numel(FEAT_NAMES);

rawTbl = readtable(DATA_FILE, TextType="string");
logInfo("Loaded %d rows from %s", height(rawTbl), DATA_FILE);

if isnumeric(rawTbl.ALogP)
    alogpVec = rawTbl.ALogP;
else
    alogpVec = str2double(rawTbl.ALogP);
end

nRaw  = height(rawTbl);
X_all = nan(nRaw, N_FEATS);
y_all = nan(nRaw, 1);
valid = false(nRaw, 1);

for k = 1:nRaw
    if isnan(alogpVec(k)), continue; end
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
logInfo("Completed descriptor calculations for %d compounds", sum(valid));

X        = X_all(valid, :);
y        = y_all(valid);
molNames = rawTbl.Name(valid);
nMols    = sum(valid);

logInfo("Dataset: %d / %d molecules (ALogP range: %.2f to %.2f)", ...
    nMols, nRaw, min(y), max(y));
%[text] **💡 Observation Point 1 — Parameter-to-Sample Ratio and Outlier Check**
%[text] Let's check the scale of the feature matrix (nMols x N\_FEATS) compared to the number of learnable parameters defined in the network in Section 4.
%[text] If the number of parameters exceeds the number of learning samples, the risk of overfitting increases. Let's check that ratio here.
%[text] Check if there are molecules with |ALogP| \> 7. You can identify them with the following code:
%[text] \[~,ord\] = sort(abs(y),"descend");
%[text] table(molNames(ord(1:5)), y(ord(1:5)), VariableNames=\["Name","ALogP"\])
%[text] Let's check if these outliers have an unbalanced impact on MSE-based neural network learning.
% ... (Try writing code here)
%%
%[text] ## Section 2: Exploring ALogP Distribution
%[text] 
%[text] Descriptor calculation is complete. Before building the model, let's examine the distribution of the target variable ALogP (a variant of LogP).
%[text] Understanding the presence of outliers and the shape of the distribution will help in correctly interpreting the learning results later.
%[text] 
%[text] ### Concept: Investigating the Target Before Learning
%[text] Neural network regressors minimize mean squared error (MSE). Since MSE penalizes large errors quadratically, outliers in y can disproportionately affect the learned weights. The network might sacrifice the accuracy of most molecules to reduce errors in a few extreme cases.
%[text] 
%[text] By visualizing the target distribution before learning, you can decide whether to apply y transformations (such as clipping or winsorization) to reduce the impact of outliers. In this exercise, we will use raw ALogP.
logSection("A05", "Section 2: Exploring ALogP Distribution", "Analytics L3");
figure("Name", "A05 ALogP Distribution");
histogram(y, 20, FaceColor=[0.2 0.6 0.8]);
xlabel("ALogP"); ylabel("Count");
title("ALogP Distribution -- FDA Approved Drugs");
xline(mean(y),   "r--", sprintf("Mean=%.2f",   mean(y)),   LabelHorizontalAlignment="right",  LabelOrientation="horizontal");
xline(median(y), "g--", sprintf("Median=%.2f", median(y)), LabelHorizontalAlignment="left",   LabelOrientation="horizontal");
grid on;

logInfo("ALogP: Mean=%.2f, Std Dev=%.2f, Skewness=%.2f", mean(y), std(y), skewness(y));
%[text] **💡 Observation Point 2 — Normality Test and Winsorization of ALogP Distribution**
%[text] Let's check if ALogP is approximately Gaussian. Test with the following code:
%[text] \[h, p\] = kstest((y - mean(y)) / std(y))
%[text] If p < 0.05, the distribution can be considered significantly non-Gaussian.
%[text] Is this result important for neural network learning (compared to linear regression)?
%[text] (Hint: OLS assumes Gaussian residuals, but NN using MSE loss does not)
%[text] Let's winsorize the target at the 5th and 95th percentiles:
%[text] y\_clip = min(max(y, prctile(y,5)), prctile(y,95));
%[text] In Sections 3 to 7, replace y with y\_clip and compare R².
% ... (Try writing the code here)
%%
%[text] ## Section 3: Training / Validation / Test Split and Feature Standardization
%[text] 
%[text] We have checked the data distribution. Next, we will perform a three-way split of the data necessary for neural network training.
%[text] The validation set is used to monitor overfitting during training, and the test set is used for final evaluation.
%[text] 
%[text] ### Concept: Three-way Split for Neural Network Training
%[text] Unlike the k-fold cross-validation used in A03, neural networks require a dedicated validation set for (a) real-time monitoring of overfitting / (b) triggering early stopping / (c) tuning hyperparameters without contaminating the test set.
%[text] 
%[text] **Convention:** 70% Training / 15% Validation / 15% Test. For 200 molecules, this is approximately 140 / 30 / 30.
%[text] **Data Leakage Rule:** Statistics for standardization (mean, standard deviation) are calculated only from the training data and applied to the validation and test sets.
logSection("A05", "Section 3: Training / Validation / Test Split and Feature Standardization", "Analytics L3");
MINI_BATCH_SIZE = 32;
VAL_FREQUENCY   = 5;   % validate every 5 training iterations

rng(42);  % fix random seed for reproducibility
cv1    = cvpartition(nMols, "HoldOut", 0.30);
trIdx  = training(cv1);
tmpIdx = test(cv1);

Xtmp = X(tmpIdx, :);
ytmp = y(tmpIdx);
cv2     = cvpartition(sum(tmpIdx), "HoldOut", 0.50);
valIdx  = training(cv2);
tstIdx  = test(cv2);

Xtr  = X(trIdx,  :);  ytr  = y(trIdx);
Xval = Xtmp(valIdx,:); yval = ytmp(valIdx);
Xte  = Xtmp(tstIdx,:); yte  = ytmp(tstIdx);
%[text] Standardize features using only the statistics from the training set.
Xmean = mean(Xtr, 1);
Xstd  = std(Xtr,  0, 1);
Xstd(Xstd < 1e-12) = 1;  % guard: avoid /0 for constant descriptors

Xtr_s  = (Xtr  - Xmean) ./ Xstd;
Xval_s = (Xval - Xmean) ./ Xstd;
Xte_s  = (Xte  - Xmean) ./ Xstd;

logInfo("Split: Training=%d / Validation=%d / Test=%d (Total=%d)", ...
    size(Xtr, 1), size(Xval, 1), size(Xte, 1), nMols);
%[text] **💡 Observation Point 3 — Reliability of Random Seed and Split**
%[text] Let's check why rng(42) is necessary. Remove it and rerun Sections 3-7 five times,
%[text] recording the test R² each time. How much variance is there?
%[text] What does this indicate about the reliability of a single test split on a dataset of about 200 molecules?
%[text] As an alternative, there is nested cross-validation (outer loop: test split; inner loop: validation split for tuning).
%[text] Explain why this provides a less biased performance estimate, albeit at a higher computational cost.
% ... (Try writing code here)
%%
%[text] ## Section 4: Defining the Neural Network Architecture
%[text] 
%[text] The data preparation is complete. Next, we will design the network structure (architecture).
%[text] We will use a shallow 2-layer network with dropout regularization suitable for small datasets.
%[text] 
%[text] ### Concept: Feedforward Network Design for Small Dataset Regression
%[text] The network maps 8 standardized descriptors to a single ALogP (estimated common logarithm) prediction through two hidden layers: `h1 = ReLU(W1×x + b1)` → Dropout(0.20) → `h2 = ReLU(W2×h1d + b2)` → `yhat = W3×h2 + b3` (linear output). The total number of learnable parameters is 2689.
%[text] 
%[text] **ReLU** — Low computational cost, avoids vanishing gradients, and generates sparse activations. **Dropout** (p=0.20) — Randomly zeroes 20% of neurons in each training forward pass to prevent co-adaptation. With about 140 training molecules and 2689 parameters, the parameter-to-sample ratio is approximately 19. Dropout is essential.
%[text] 
logSection("A05", "Section 4: Defining the Neural Network Architecture", "Analytics L3");
layers = [
    featureInputLayer(N_FEATS, Normalization="none")
    fullyConnectedLayer(64)
    reluLayer
    dropoutLayer(0.20)
    fullyConnectedLayer(32)
    reluLayer
    fullyConnectedLayer(1)];

net = dlnetwork(layers);
%[text] Count the learnable parameters
numParams = 0;
for k = 1:height(net.Learnables)
    numParams = numParams + numel(extractdata(net.Learnables.Value{k}));
end
logInfo("Architecture: FC(64)-ReLU-Drop(0.2)-FC(32)-ReLU-FC(1)");
logInfo("Learnable Parameters: %d  (Parameter / Training Sample Ratio = %.1f)", ...
    numParams, numParams / size(Xtr_s, 1));
%[text] **💡 Observation Point 4 — Comparing Architecture Changes and Effects**
%[text] Let's add a third hidden layer (fullyConnectedLayer(16) + reluLayer) before the output layer to observe the effect.
%[text] Will a deeper network improve test R2? How will the cost (training time, overfitting risk) change?
%[text] Replace reluLayer with tanhLayer to observe the impact. How will saturated activation functions affect learning speed and final performance on this small dataset?
%[text] (Hint: Near the saturation region of tanh, gradients approach zero — this is known as the vanishing gradient problem)
%[text] Set dropoutLayer(0.50) instead of 0.20 to observe the impact.
%[text] Will there be more or less overfitting observed in the learning curve (Section 6)? Consider the reasons.
% ... (Try writing the code here)
%%
%[text] ## Section 5: Network Training
%[text] 
%[text] The architecture definition is complete. Next, we will train the network.
%[text] In the progress window displayed during training, let's check the transition of training loss and validation loss.
%[text] 
%[text] ### Concept: Adam Optimizer and Piecewise Learning Rate Schedule
%[text] Adam adaptively adjusts the learning rate for each parameter using the exponential moving average of past gradients. In the piecewise schedule, the learning rate is halved at epoch 80: epochs 1-80 have lr=1e-3, epochs 81-200 have lr=5e-4.
%[text] 
%[text] **Mini-batch SGD:** Each Adam step uses a random subset of 32 molecules. With 140 training molecules and MiniBatchSize=32, 1 epoch consists of 4 iterations. The training progress window is displayed automatically.
%[text] 
logSection("A05", "Section 5: Network Training", "Analytics L3");
opts = trainingOptions("adam", ...
    MaxEpochs=200, ...
    MiniBatchSize=MINI_BATCH_SIZE, ...
    InitialLearnRate=1e-3, ...
    LearnRateSchedule="piecewise", ...
    LearnRateDropFactor=0.5, ...
    LearnRateDropPeriod=80, ...
    ValidationData={Xval_s, yval}, ...
    ValidationFrequency=VAL_FREQUENCY, ...
    Plots="training-progress", ...
    Verbose=false);

logInfo("Training for 200 epochs (approximately 30-60 seconds)...");
[net, trainInfo] = trainnet(Xtr_s, ytr, net, "mse", opts);
%[text] Explanation of version-compatible loss extraction methods.
%[text] R2024a: trainInfo.TrainingLoss / trainInfo.ValidationLoss (direct arrays)
%[text] R2025b: TrainingHistory.Loss / ValidationHistory.Loss (separate tables, column name "Loss")
%[text] R2025a: TrainingHistory (table with TrainingLoss, ValidationLoss columns)
try
    % R2024a style: direct array fields
    trLossAll   = double(trainInfo.TrainingLoss);
    valLossAll  = double(trainInfo.ValidationLoss);
    valIterNums = [];
catch
    try
        % R2025b style: separate TrainingHistory / ValidationHistory tables, column "Loss"
        trLossAll   = double(trainInfo.TrainingHistory.Loss);
        valLossAll  = double(trainInfo.ValidationHistory.Loss);
        valIterNums = double(trainInfo.ValidationHistory.Iteration);
    catch
        % R2025a style: single TrainingHistory table with TrainingLoss/ValidationLoss columns
        h           = trainInfo.TrainingHistory;
        trLossAll   = double(h.TrainingLoss);
        valMask     = ~isnan(h.ValidationLoss);
        valLossAll  = double(h.ValidationLoss(valMask));
        valIterNums = [];
    end
end

finalValLoss = valLossAll(end);
logInfo("Training complete. Final validation MSE=%.4f (RMSE=%.3f)", ...
    finalValLoss, sqrt(max(finalValLoss, 0)));
%%
%[text] ## Section 6: Learning Curve Analysis
%[text] 
%[text] Training is complete. Next, let's analyze the learning curves in detail.
%[text] By examining the difference between training loss and validation loss, you can diagnose whether overfitting is occurring.
%[text] 
%[text] Concept: Diagnose learning dynamics using loss curves
%[text] ### Concept: How to Read Learning Curves
%[text] Learning curves plot the change in loss against training epochs. **Good Fit:** Both training and validation mean squared errors decrease and converge to similar values. **Overfitting:** Training loss continues to decrease, but validation loss stagnates or increases. **Underfitting:** Both losses remain high with slow progress.
%[text] 
logSection("A05", "Section 6: Learning Curve Analysis", "Analytics L3");
% Calculate itersPerEpoch from actual training history (accurate with formula comparison)
nIterTotal    = numel(trLossAll);
nValPoints    = numel(valLossAll);
if isprop(trainInfo, 'TrainingHistory') && ismember('Epoch', trainInfo.TrainingHistory.Properties.VariableNames)
    itersPerEpoch = max(1, round(nIterTotal / max(trainInfo.TrainingHistory.Epoch)));
else
    itersPerEpoch = max(1, ceil(size(Xtr_s, 1) / MINI_BATCH_SIZE));
end

epochsTr  = (1:nIterTotal) / itersPerEpoch;
if ~isempty(valIterNums)
    epochsVal = valIterNums / itersPerEpoch;
else
    epochsVal = (1:nValPoints) * VAL_FREQUENCY / itersPerEpoch;
end
%[text] To facilitate interpretation, convert MSE to RMSE (same units as ALogP).
trRmse  = sqrt(max(trLossAll,  0));
valRmse = sqrt(max(valLossAll, 0));
%[text] Smooth the learning curves using a 1-epoch moving average.
winLen    = max(1, itersPerEpoch);
trSmoothed = movmean(trRmse, winLen, "Endpoints", "fill");

figure("Name", "A05 Learning Curves");
plot(epochsTr,  trRmse,    "Color", [0.7 0.8 1.0], LineWidth=0.8, ...
    DisplayName="Training (Raw)"); hold on;
plot(epochsTr,  trSmoothed, "b-", LineWidth=2.0, DisplayName="Training (Smoothed)");
plot(epochsVal, valRmse,    "r--", LineWidth=2.0, DisplayName="Validation");
xlabel("Epoch"); ylabel("RMSE (ALogP Units)");
title("A05: Learning Curves");
legend(Location="northeast");
grid on;
%[text] Mark the best validation epoch.
[minValRmse, minValIdx] = min(valRmse);
xline(epochsVal(minValIdx), "--k", ...
    sprintf("Best Validation RMSE=%.3f (Epoch %.0f)", minValRmse, epochsVal(minValIdx)), ...
    LabelHorizontalAlignment="left", LabelVerticalAlignment="bottom", ...
    LabelOrientation="horizontal", HandleVisibility="off");

logInfo("Best Validation RMSE=%.3f (Epoch %.0f)", minValRmse, epochsVal(minValIdx));
%[text] **💡 Observation Point 5 — Diagnosing Learning Curves and Early Stopping**
%[text] Do the learning curves indicate overfitting? (Training loss continues to decrease while validation loss increases)
%[text] Check at which epoch the gap between training and validation RMSE starts to widen.
%[text] Try adding early stopping to the trainingOptions call in Section 5.
%[text] ValidationPatience=20 (stop if no improvement after 20 checks)
%[text] Check if early stopping improves, worsens, or has no effect on test R2.
%[text] Consider what happens when MiniBatchSize=nMols (full gradient descent).
%[text] Compare the expected shape of learning curves with mini-batch.
% ... (Try writing code here)
%%
%[text] ## Section 7: Test Set Evaluation
%[text] 
%[text] We have verified the model's behavior using the learning curve. Finally, we will fairly evaluate the model's performance using the test set.
%[text] The test set is used here for the first time. This result will be the final performance indicator.
%[text] 
%[text] ### Concept: One-time Test Set Evaluation
%[text] The test set is evaluated exactly once after all design decisions (architecture, hyperparameters) are finalized. **Key regression metrics:** $R^2 = 1 - SS\_{res}/SS\_{tot}$ / RMSE (same units as y) / MAE (robust to outliers). In the actual vs predicted scatter plot, points on the y=x diagonal indicate perfect predictions.
%[text] 
logSection("A05", "Section 7: Test Set Evaluation", "Analytics L3");
yPred = double(predict(net, Xte_s));   % (nTest x 1) predicted ALogP
res   = yte - yPred;
rmse  = sqrt(mean(res.^2));
mae   = mean(abs(res));
r2    = 1 - sum(res.^2) / sum((yte - mean(yte)).^2);

logInfo("Test Set: R2=%.3f  RMSE=%.3f  MAE=%.3f  (n=%d)", r2, rmse, mae, numel(yte));

figure("Name", "A05 Test Predictions");
tiledlayout(1, 2);

nexttile;
scatter(yte, yPred, 50, "filled", MarkerFaceAlpha=0.7, ...
    MarkerFaceColor=[0.2 0.5 0.9]);
hold on;
lim = [min([yte; yPred]) - 0.5, max([yte; yPred]) + 0.5];
plot(lim, lim, "k--", LineWidth=1.5);
xlabel("Actual ALogP"); ylabel("Predicted ALogP");
title(sprintf("Actual vs Predicted  (R^2=%.3f)", r2));
grid on;

nexttile;
histogram(res, 15, FaceColor=[0.9 0.4 0.3]);
xline(0, "--k", LineWidth=1.5);
xlabel("Residuals (Actual - Predicted)"); ylabel("Count");
title(sprintf("Residuals  (RMSE=%.3f, MAE=%.3f)", rmse, mae));
grid on;
%[text] **💡 Observation Point 6 — Identifying Residual Bias and Difficult-to-Predict Molecules**
%[text] Are the residuals centered around zero? A non-zero mean indicates systematic bias.
%[text] (This means the network consistently over- or under-predicts).
%[text] What causes this in neural networks with small datasets?
%[text] Identify the three worst-predicted molecules (maximum |residual|).
%[text] Are they structural outliers, or do they share a common scaffold?
%[text] \[~,ord\] = sort(abs(res), "descend");
%[text] tstNames = molNames(tmpIdx);  % Narrow down to tmpIdx subset
%[text] table(tstNames(tstIdx)(ord(1:3)), yte(ord(1:3)), yPred(ord(1:3)), ...
%[text] VariableNames=\["Name","Actual","Predicted"\])
% ... (Let's write the code here)
%%
%[text] ## Section 8: Linear Baseline Comparison and Discussion
%[text] 
%[text] The evaluation of the neural network is complete. Finally, we will compare it with linear regression using the same data.
%[text] Let's confirm the important QSAR lesson that "more complex models are not always better."
%[text] 
%[text] ### Concept: Occam's Razor — Choose the simplest and most appropriate model
%[text] In QSAR, the "simplest and most appropriate model" is usually linear regression. Replacing it with a neural network is limited to cases where (a) the data per descriptor is large with N\>1000, (b) there is a theoretical basis to expect non-linearity, or (c) the performance of a cross-validated neural network significantly exceeds that of linear regression.
%[text] 
%[text] In drug datasets of 100-300 molecules, linear regression, random forest, and Gaussian process models usually match or exceed neural networks. GNN is where deep learning truly shines in cheminformatics. For fairness in this comparison, we use the same training/test split.
%[text] 
logSection("A05", "Section 8: Linear Baseline Comparison and Discussion", "Analytics L3");
ftTbl    = array2table(Xtr_s, VariableNames=cellstr(FEAT_NAMES));
ftTbl.ALogP = ytr;
lmMdl    = fitlm(ftTbl, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));

teTbl    = array2table(Xte_s, VariableNames=cellstr(FEAT_NAMES));
yLinPred = predict(lmMdl, teTbl);
r2Lin    = 1 - sum((yte - yLinPred).^2) / sum((yte - mean(yte)).^2);
rmseLin  = sqrt(mean((yte - yLinPred).^2));
maeLin   = mean(abs(yte - yLinPred));

logInfo("Linear Regression: R2=%.3f  RMSE=%.3f  MAE=%.3f", r2Lin, rmseLin, maeLin);
logInfo("Neural Network:    R2=%.3f  RMSE=%.3f  MAE=%.3f", r2, rmse, mae);
logInfo("NN Advantage:      dR2=%.3f  dRMSE=%+.3f  dMAE=%+.3f", ...
    r2 - r2Lin, rmseLin - rmse, maeLin - mae);
%[text] Compare the two models side by side.
figure("Name", "A05 Model Comparison");
set(gcf, "Position", [100 100 1100 400]);
tiledlayout(1, 3);

metrics  = [r2Lin, r2; rmseLin, rmse; maeLin, mae];
ylabels  = ["R^2 (Higher is better)", "RMSE (Lower is better)", "MAE (Lower is better)"];
for m = 1:3
    nexttile;
    bar(metrics(m, :), FaceColor="flat", CData=[0.4 0.6 0.9; 0.9 0.4 0.3]);
    set(gca, "XTickLabel", ["Linear", "Neural Network"]);
    ylabel(ylabels(m)); title(ylabels(m)); grid on;
    ylim([0, max(metrics(m, :)) * 1.3]);
end
sgtitle("A05: Linear Regression vs Neural Network (Same Train/Test Split)");
%[text] **💡 Observation Point 7 — Comparison of Neural Network and Linear Model, and Prospects for GNN**
%[text] Let's check if the neural network outperforms linear regression on this dataset.
%[text] If the linear model is superior, list three specific reasons why.
%[text] Consider what needs to change for the neural network to gain a clear advantage.
%[text] Examine dataset size (N), feature type (fingerprint vs descriptor), and property type (activity vs physicochemical properties).
%[text] Graph Neural Networks (GNN) represent molecules as graphs rather than fixed-length descriptor vectors. How does this change the feature engineering step?
%[text] And explain why better generalization is possible with larger datasets.
%[text] (Reference: Gilmer et al. 2017, Neural Message Passing)
% ... (Try writing code here)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
