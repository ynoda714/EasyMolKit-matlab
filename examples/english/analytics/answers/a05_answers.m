%[text] # A05 Answers: Neural Network Property Prediction
%[text] Reference answers for the "Try It" exercise in a05_neural_net_property.m.
%[text] First, run a05_neural_net_property.m or execute standalone.
%[text]
%[text] Required: Deep Learning Toolbox, Statistics & Machine Learning Toolbox
addpath(genpath("src"));
emk.setup.initPython();
logInfo("A05 Answers: Setup complete");

%[text] Rebuild workspace (for standalone execution)
DATA_FILE  = "data/list/fda_drugs.csv";
FEAT_NAMES = ["MolWt","TPSA","NumHDonors","NumHAcceptors", ...
              "NumRotatableBonds","RingCount","FractionCSP3","HeavyAtomCount"];
N_FEATS    = numel(FEAT_NAMES);
TEST_RATIO = 0.20;
VAL_RATIO  = 0.15;
rng(42);

rawTbl  = readtable(DATA_FILE, TextType="string");
nRaw    = height(rawTbl);
if isnumeric(rawTbl.ALogP)
    alogpVec = rawTbl.ALogP;
else
    alogpVec = str2double(rawTbl.ALogP);
end

X_all = nan(nRaw, N_FEATS);
y_all = nan(nRaw, 1);
valid = false(nRaw, 1);
for k = 1:nRaw
    if isnan(alogpVec(k)), continue; end
    try
        mol = emk.mol.fromSmiles(rawTbl.SMILES(k));
        s   = emk.descriptor.calculate(mol, FEAT_NAMES);
        for d = 1:N_FEATS; X_all(k, d) = s.(FEAT_NAMES(d)); end
        y_all(k) = alogpVec(k);
        valid(k) = true;
    catch
    end
    logProgress(k, nRaw, "Calculating descriptors");
end

X        = X_all(valid, :);
y        = y_all(valid);
molNames = rawTbl.Name(valid);
nMols    = sum(valid);

Xmean = mean(X, 1);
Xstd  = std(X, 0, 1);
Xstd(Xstd < 1e-12) = 1;
Xs = (X - Xmean) ./ Xstd;

cv   = cvpartition(nMols, Holdout=TEST_RATIO);
trIdx = training(cv);
teIdx = test(cv);
cvVal = cvpartition(sum(trIdx), Holdout=VAL_RATIO / (1 - TEST_RATIO));
trSub = find(trIdx);
valIdx_rel = test(cvVal);
trainMask  = true(sum(trIdx), 1);
trainMask(valIdx_rel) = false;

trIdxFull = trSub(trainMask);
valIdxFull = trSub(valIdx_rel);

XTr = Xs(trIdxFull, :);  yTr = y(trIdxFull);
XV  = Xs(valIdxFull, :); yV  = y(valIdxFull);
XTe = Xs(teIdx, :);      yTe = y(teIdx);
nTr = numel(yTr);

logInfo("Dataset: %d training / %d validation / %d test", nTr, numel(yV), numel(yTe));
%%
%[text] ## Let's Try 1: Parameter-to-Sample Ratio; ALogP Outliers
%[text] Number of trainable parameters for baseline FC64-ReLU-Drop0.2-FC32-ReLU-FC1
p1 = N_FEATS * 64 + 64;    % fc1 weights + bias
p2 = 64 * 32 + 32;          % fc2 weights + bias
p3 = 32 * 1 + 1;            % output layer
totalParams = p1 + p2 + p3;
ratio = totalParams / nTr;
logInfo("Baseline network: %d trainable parameters", totalParams);
logInfo("Parameter / Training Sample Ratio: %.1f (Parameters / Training Samples %d)", ratio, nTr);

%[text] ALogP Outliers
[sortedY, ord] = sort(y, "descend");
logInfo("Top 5 Highest ALogP:");
for k = 1:5
    logInfo("  %d. %-30s  ALogP = %.2f", k, molNames(ord(k)), sortedY(k));
end
pct95 = prctile(y, 95);
logInfo("ALogP 95th Percentile Value: %.2f (Above %d molecules)", ...
    pct95, sum(y > pct95));

%[text] Answer: The baseline has about 2700 trainable parameters for approximately 130 training samples
%[text]    (The exact number varies depending on the number of active drugs in the training split).
%[text]    A ratio of ~20 indicates the network is appropriately parameterized for the data,
%[text]    but it is already approaching the "more parameters than samples" region.
%[text]    Adding a third hidden layer increases the ratio above 30, raising the risk of overfitting.
%[text]    Compounds with ALogP > 6 are mainly steroid hormones and cyclic lipid compounds.
%%
%[text] ## Let's Try 2: Distribution Shape; Normality Test; Winsorization

hKS = kstest((y - mean(y)) / std(y));
logInfo("Kolmogorov-Smirnov Normality Test: H=%d  (H=1 -> Non-normal)", hKS);
logInfo("Skewness=%.3f  Kurtosis=%.3f", skewness(y), kurtosis(y));

%[text] Winsorize at 1-99 percentiles
lo = prctile(y, 1);
hi = prctile(y, 99);
yWin = max(min(y, hi), lo);
logInfo("Winsorized range: [%.2f, %.2f] (Original: [%.2f, %.2f])", ...
    lo, hi, min(y), max(y));
logInfo("Values changed: %d", sum(yWin ~= y));

%[text] Compare linear baseline RMSE of original and Winsorized data
ftTbl = array2table(Xs, VariableNames=cellstr(FEAT_NAMES));
ftTbl.y = y;
lm_raw = fitlm(ftTbl, "y ~ " + strjoin(FEAT_NAMES, " + "));
ftTbl2 = ftTbl; ftTbl2.y = yWin;
lm_win = fitlm(ftTbl2, "y ~ " + strjoin(FEAT_NAMES, " + "));
logInfo("Linear RMSE: Original=%.3f  Winsorized=%.3f", lm_raw.RMSE, lm_win.RMSE);

%[text] Answer: K-S test H=1 indicates that the ALogP distribution is non-normal (right-skewed due to high lipophilicity outliers).
%[text]    Neural networks minimize MSE, so outliers impact proportionally to squared residuals.
%[text]    Winsorization caps extreme values: typically only 2-4 molecules are affected.
%[text]    RMSE with Winsorized targets is lower, but predictive ability for extreme values is lost.
%%
%[text] ## Let's Try 3: Purpose of rng(42); Variability without a fixed seed

logInfo("Variance of 5 NN runs without a fixed seed:");
rmsesNoSeed = nan(1, 5);
for rep = 1:5
    netTmp = [featureInputLayer(N_FEATS), ...
              fullyConnectedLayer(64), reluLayer, dropoutLayer(0.2), ...
              fullyConnectedLayer(32), reluLayer, ...
              fullyConnectedLayer(1),  regressionLayer];
    opts = trainingOptions("adam", MaxEpochs=30, MiniBatchSize=32, ...
        ValidationData={XV', yV'}, ValidationFrequency=10, ...
        Verbose=false, Plots="none", ...
        LearnRateSchedule="none", InitialLearnRate=1e-3);
    netT = trainNetwork(XTr', yTr', netTmp, opts);
    pred = predict(netT, XTr')';
    rmsesNoSeed(rep) = sqrt(mean((yTr - pred(:)).^2));
end
logInfo("Training RMSE of 5 runs: Mean=%.3f  Std Dev=%.3f", ...
    mean(rmsesNoSeed), std(rmsesNoSeed));

%[text] Answer: rng(42) fixes the random seed for weight initialization and dropout masks.
%[text]    Without this, RMSE can vary by 0.1 to 0.3 units depending on initial weight configuration.
%[text]    Fixing the seed is essential for reproducible experiments and fair comparison between architectures.
%%
%[text] ## Let's Try 4: Add a 3rd Hidden Layer; tanhLayer; dropoutLayer(0.50)
%[text] Deeper Architecture
rng(42);
netDeep = [featureInputLayer(N_FEATS), ...
           fullyConnectedLayer(64), reluLayer, dropoutLayer(0.2), ...
           fullyConnectedLayer(32), reluLayer, dropoutLayer(0.2), ...
           fullyConnectedLayer(16), reluLayer, ...
           fullyConnectedLayer(1),  regressionLayer];
p_deep = (N_FEATS*64+64) + (64*32+32) + (32*16+16) + (16*1+1);
logInfo("Deep network parameter count: %d (ratio=%.1f)", p_deep, p_deep/nTr);

%[text] Replace ReLU with tanh
rng(42);
netTanh = [featureInputLayer(N_FEATS), ...
           fullyConnectedLayer(64), tanhLayer, dropoutLayer(0.2), ...
           fullyConnectedLayer(32), tanhLayer, ...
           fullyConnectedLayer(1),  regressionLayer];

%[text] Dropout 0.50
rng(42);
netDrop50 = [featureInputLayer(N_FEATS), ...
             fullyConnectedLayer(64), reluLayer, dropoutLayer(0.5), ...
             fullyConnectedLayer(32), reluLayer, ...
             fullyConnectedLayer(1),  regressionLayer];

opts30 = trainingOptions("adam", MaxEpochs=30, MiniBatchSize=32, ...
    ValidationData={XV', yV'}, ValidationFrequency=10, ...
    Verbose=false, Plots="none", InitialLearnRate=1e-3);

logInfo("Training 3 architectures for comparison...");
[netD, tiD] = trainNetwork(XTr', yTr', netDeep,   opts30);
[netT, tiT] = trainNetwork(XTr', yTr', netTanh,   opts30);
[netP, tiP] = trainNetwork(XTr', yTr', netDrop50, opts30);

archNets  = {netD,    netT,    netP};
archNames = {"Deep3", "Tanh",  "Drop0.5"};
for aIdx = 1:3
    yTePredA = predict(archNets{aIdx}, XTe')';
    rmseTE   = sqrt(mean((yTe - yTePredA(:)).^2));
    logInfo("%-8s  Test RMSE = %.3f", archNames{aIdx}, rmseTE);
end

%[text] Answer: Adding a 3rd hidden layer to a small dataset (n~130 training) usually has no effect or worsens performance.
%[text]    Parameters increase but unique learning signals do not, worsening overfitting.
%[text]    Tanh can saturate with large activations (gradient vanishing), sometimes showing smooth behavior on small data.
%[text]    Dropout 0.5 is a strong regularizer: validation RMSE often improves, but training RMSE increases.
%%
%[text] ## Let's Try 5: Early Stopping; ValidationPatience; Full Batch Gradient Descent

rng(42);
netBase = [featureInputLayer(N_FEATS), ...
           fullyConnectedLayer(64), reluLayer, dropoutLayer(0.2), ...
           fullyConnectedLayer(32), reluLayer, ...
           fullyConnectedLayer(1),  regressionLayer];
opts150 = trainingOptions("adam", MaxEpochs=150, MiniBatchSize=32, ...
    ValidationData={XV', yV'}, ValidationFrequency=5, ...
    Verbose=false, Plots="none", InitialLearnRate=1e-3, ...
    ValidationPatience=20);
[netBase, tiBase] = trainNetwork(XTr', yTr', netBase, opts150);

valLoss = tiBase.ValidationLoss;
logInfo("Training stopped at epoch %d out of a maximum of 150 epochs", numel(valLoss)*5);
[minValLoss, minEpIdx] = min(valLoss);
logInfo("Best validation MSE: %.4f (around epoch %d) RMSE=%.3f", ...
    minValLoss, minEpIdx*5, sqrt(max(minValLoss,0)));

%[text] Identifying the start of overfitting: the point where training loss continues to decrease but validation loss plateaus
trainLoss = tiBase.TrainingLoss;
if numel(trainLoss) == numel(valLoss)
    gap = valLoss - trainLoss;
    [~, gapPeak] = max(gap);
    logInfo("Start of overfitting (maximum train-validation gap): around epoch %d", gapPeak*5);
end

%[text] Answer: ValidationPatience=20 stops training when the validation metric does not improve for 20 consecutive checks.
%[text]    Each check is every 5 epochs, so it stops after 100 epochs without improvement.
%[text]    Full batch GD (MiniBatchSize=nTr) results in a smooth loss curve but slow training per epoch.
%[text]    The start of overfitting is the epoch where validation loss begins to diverge from training loss — typically around 30-60 epochs for this architecture.
%%
%[text] ## Let's Try 6: Residual Bias; Worst Predicted Molecules

rng(42);
yTePred = predict(netBase, XTe')';
residuals = yTe - yTePred(:);
logInfo("Test Residuals: Mean=%+.3f  Std Dev=%.3f  RMSE=%.3f", ...
    mean(residuals), std(residuals), sqrt(mean(residuals.^2)));
logInfo("Residual Bias: %s", ...
    string(ternary_(abs(mean(residuals)) > 0.2, "Biased", "Almost Unbiased")));

%[text] Worst Predicted Molecules
[~, worstOrd] = sort(abs(residuals), "descend");
teNames = molNames(find(teIdx));  %#ok<FNDSB>
logInfo("Top 3 Worst Predictions:");
for k = 1:min(3, numel(worstOrd))
    i = worstOrd(k);
    logInfo("  %-28s  actual=%.2f  pred=%.2f  err=%+.2f", ...
        teNames(i), yTe(i), yTePred(i), residuals(i));
end

%[text] Residuals vs Predicted Values Plot
figure("Name", "A05 Residual Analysis");
scatter(yTePred(:), residuals, 40, [0.3 0.6 0.9], "filled", MarkerFaceAlpha=0.7);
hold on; yline(0, "r--"); hold off;
xlabel("Predicted ALogP"); ylabel("Residuals (Actual - Predicted)");
title("Residual Plot -- Test Set"); grid on;

%[text] Answer: A non-zero mean residual indicates systematic bias (the network tends to overpredict or underpredict on average).
%[text]    For ALogP, there is a tendency to underpredict highly lipophilic steroids (too few in training).
%[text]    The worst predicted molecules are often extreme outliers or belong to structural classes that were underrepresented in the training data.
%%
%[text] ## Let's Try 7: NN vs Linear Comparison; When NN Wins; Explanation of GNN

ftTbl  = array2table(Xs(trIdx, :), VariableNames=cellstr(FEAT_NAMES));
ftTbl.ALogP = y(trIdx);
lmFull = fitlm(ftTbl, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));
yTe_lm = predict(lmFull, array2table(XTe, VariableNames=cellstr(FEAT_NAMES)));
rmse_lm = sqrt(mean((yTe - yTe_lm).^2));
rmse_nn = sqrt(mean((yTe - yTePred(:)).^2));
logInfo("Test RMSE: Linear=%.3f  NN=%.3f", rmse_lm, rmse_nn);
logInfo("NN Advantage: %+.3f  (%s)", rmse_lm - rmse_nn, ...
    string(ternary_(rmse_nn < rmse_lm, "NN Wins", "Linear Wins")));

%[text] Answer: In ALogP prediction using tabular descriptors, the advantage of NN over linear regression is slight and can be zero or negative.
%[text]    Descriptors already encode relevant chemical information, making linear combinations a reasonable approximation.
%[text]    NN becomes more powerful in the following cases:
%[text]      (1) When the input is raw molecular graphs or fingerprints (2048 bits),
%[text]      (2) When the relationship is highly nonlinear (e.g., activity cliffs in pIC50),
%[text]      (3) When the training data is large (n > 10,000).
%[text]    GNN directly takes molecular graphs as input (atoms as nodes, bonds as edges) and learns atom embeddings by aggregating neighborhood information.
%[text]    GNNs like MPNN and AttentiveFP achieve SOTA results on the MoleculeNet benchmark.
logInfo("A05 Answer Completed.");

%[text] Local Helper Function
function out = ternary_(cond, a, b)
    if cond; out = a; else; out = b; end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]%   data: {"layout":"inline","rightPanelPercent":40}
%---
