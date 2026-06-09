%[text] # A05 解答: ニューラルネットワークプロパティ予測
%[text] a05_neural_net_property.m の「やってみよう」演習の参照解答。
%[text] 最初に a05_neural_net_property.m を実行するか、スタンドアロンで実行。
%[text]
%[text] 必要: Deep Learning Toolbox、Statistics & Machine Learning Toolbox
addpath(genpath("src"));
emk.setup.initPython();
logInfo("A05 解答: セットアップ完了");

%[text] ワークスペース再構築（スタンドアロン実行用）
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
    logProgress(k, nRaw, "記述子を計算中");
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

logInfo("データセット: %d 学習 / %d 検証 / %d テスト", nTr, numel(yV), numel(yTe));
%%
%[text] ## やってみよう 1: パラメータ対サンプル比; ALogP 外れ値
%[text] ベースライン FC64-ReLU-Drop0.2-FC32-ReLU-FC1 の学習可能パラメータ数
p1 = N_FEATS * 64 + 64;    % fc1 weights + bias
p2 = 64 * 32 + 32;          % fc2 weights + bias
p3 = 32 * 1 + 1;            % output layer
totalParams = p1 + p2 + p3;
ratio = totalParams / nTr;
logInfo("ベースラインネットワーク: %d 学習可能パラメータ", totalParams);
logInfo("パラメータ / 学習サンプル比: %.1f（パラメータ数 / 学習サンプル %d）", ratio, nTr);

%[text] ALogP 外れ値
[sortedY, ord] = sort(y, "descend");
logInfo("最高 ALogP 上位 5:");
for k = 1:5
    logInfo("  %d. %-30s  ALogP = %.2f", k, molNames(ord(k)), sortedY(k));
end
pct95 = prctile(y, 95);
logInfo("ALogP 95パーセンタイル値: %.2f（以上 %d 分子）", ...
    pct95, sum(y > pct95));

%[text] 解答: ベースラインは約 130 学習サンプルに対して約 2700 の学習可能パラメータを持ちます
%[text]    （正確な数は学習分割の有効な薬の数により変わります）。
%[text]    比率 ~20 はデータに対してネットワークが適切にパラメータ化されていることを示しますが、
%[text]    既に「パラメータがサンプル数より多い」領域に近づいています。
%[text]    3 番目の隠れ層を追加すると比率が 30 を超え、過学習のリスクが高まります。
%[text]    ALogP > 6 の化合物は主にステロイドホルモンや環状脂質系化合物です。
%%
%[text] ## やってみよう 2: 分布形状; 正規性検定; ウィンソライゼーション

hKS = kstest((y - mean(y)) / std(y));
logInfo("Kolmogorov-Smirnov 正規性検定: H=%d  (H=1 -> 非正規)", hKS);
logInfo("歪度=%.3f  尖度=%.3f", skewness(y), kurtosis(y));

%[text] 1～99 パーセンタイルでウィンソライズ
lo = prctile(y, 1);
hi = prctile(y, 99);
yWin = max(min(y, hi), lo);
logInfo("ウィンソライズ後範囲: [%.2f, %.2f]（元: [%.2f, %.2f]）", ...
    lo, hi, min(y), max(y));
logInfo("変更された値: %d", sum(yWin ~= y));

%[text] 元データとウィンソライズ済みデータの線形ベースライン RMSE 比較
ftTbl = array2table(Xs, VariableNames=cellstr(FEAT_NAMES));
ftTbl.y = y;
lm_raw = fitlm(ftTbl, "y ~ " + strjoin(FEAT_NAMES, " + "));
ftTbl2 = ftTbl; ftTbl2.y = yWin;
lm_win = fitlm(ftTbl2, "y ~ " + strjoin(FEAT_NAMES, " + "));
logInfo("線形 RMSE: 元データ=%.3f  ウィンソライズ済み=%.3f", lm_raw.RMSE, lm_win.RMSE);

%[text] 解答: K-S 検定 H=1 は ALogP 分布が非正規（高親油性外れ値による右歪り）であることを意味します。
%[text]    ニューラルネットワークは MSE を最小化するため、外れ値は二乗残差に比例して影響します。
%[text]    ウィンソライゼーションは極端値を上限設定します：通常 2～4 分子のみ影響を受けます。
%[text]    ウィンソライズ済みターゲットでの RMSE は低くなりますが、極端値の予測能力は失われます。
%%
%[text] ## やってみよう 3: rng(42) の目的; 固定シードなしの分散

logInfo("シード固定なしの 5 回 NN 実行の分散:");
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
logInfo("5 回実行の学習 RMSE: 平均=%.3f  標準偏差=%.3f", ...
    mean(rmsesNoSeed), std(rmsesNoSeed));

%[text] 解答: rng(42) は重みの初期化とドロップアウトマスクのランダムシードを固定します。
%[text]    これなしでは、初期重みの構成によって RMSE が 0.1～0.3 単位変動することがあります。
%[text]    シードを固定することは再現可能な実験とアーキテクチャ間の公平な比較に不可欠です。
%%
%[text] ## やってみよう 4: 3 番目の隠れ層を追加; tanhLayer; dropoutLayer(0.50)
%[text] より深いアーキテクチャ
rng(42);
netDeep = [featureInputLayer(N_FEATS), ...
           fullyConnectedLayer(64), reluLayer, dropoutLayer(0.2), ...
           fullyConnectedLayer(32), reluLayer, dropoutLayer(0.2), ...
           fullyConnectedLayer(16), reluLayer, ...
           fullyConnectedLayer(1),  regressionLayer];
p_deep = (N_FEATS*64+64) + (64*32+32) + (32*16+16) + (16*1+1);
logInfo("深層ネット パラメータ数: %d（比率=%.1f）", p_deep, p_deep/nTr);

%[text] ReLU を tanh に置き換え
rng(42);
netTanh = [featureInputLayer(N_FEATS), ...
           fullyConnectedLayer(64), tanhLayer, dropoutLayer(0.2), ...
           fullyConnectedLayer(32), tanhLayer, ...
           fullyConnectedLayer(1),  regressionLayer];

%[text] ドロップアウト 0.50
rng(42);
netDrop50 = [featureInputLayer(N_FEATS), ...
             fullyConnectedLayer(64), reluLayer, dropoutLayer(0.5), ...
             fullyConnectedLayer(32), reluLayer, ...
             fullyConnectedLayer(1),  regressionLayer];

opts30 = trainingOptions("adam", MaxEpochs=30, MiniBatchSize=32, ...
    ValidationData={XV', yV'}, ValidationFrequency=10, ...
    Verbose=false, Plots="none", InitialLearnRate=1e-3);

logInfo("比較のために 3 つのアーキテクチャを学習中...");
[netD, tiD] = trainNetwork(XTr', yTr', netDeep,   opts30);
[netT, tiT] = trainNetwork(XTr', yTr', netTanh,   opts30);
[netP, tiP] = trainNetwork(XTr', yTr', netDrop50, opts30);

archNets  = {netD,    netT,    netP};
archNames = {"Deep3", "Tanh",  "Drop0.5"};
for aIdx = 1:3
    yTePredA = predict(archNets{aIdx}, XTe')';
    rmseTE   = sqrt(mean((yTe - yTePredA(:)).^2));
    logInfo("%-8s  テスト RMSE = %.3f", archNames{aIdx}, rmseTE);
end

%[text] 解答: 小規模データセット（n~130 学習）で 3 番目の隠れ層を追加しても、通常効果はないか悪化します。
%[text]    パラメータは増えますがユニークな学習シグナルは増えないため、過学習が悪化します。
%[text]    tanh は大きな活性化で飽和する（勾配消失）ため、小規模データでは滑らかな挙動を示すこともあります。
%[text]    ドロップアウト 0.5 は強い正則化子です：検証 RMSE は改善することが多いですが、学習 RMSE は高くなります。
%%
%[text] ## やってみよう 5: 早期終了; ValidationPatience; フルバッチ勾配降下法

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
logInfo("最大 150 エポック中 エポック %d で学習停止", numel(valLoss)*5);
[minValLoss, minEpIdx] = min(valLoss);
logInfo("最良検証 MSE: %.4f（約エポック %d） RMSE=%.3f", ...
    minValLoss, minEpIdx*5, sqrt(max(minValLoss,0)));

%[text] 過学習開始点の特定: 学習損失は下降し続けるが検証損失が横ばいになる点
trainLoss = tiBase.TrainingLoss;
if numel(trainLoss) == numel(valLoss)
    gap = valLoss - trainLoss;
    [~, gapPeak] = max(gap);
    logInfo("過学習開始点（学習検証ギャップ最大）: 約エポック %d", gapPeak*5);
end

%[text] 解答: ValidationPatience=20 は 20 回連続して検証指標が改善しないとき学習を停止します。
%[text]    各チェックは 5 エポックごとなので、100 エポック改善なしで停止します。
%[text]    フルバッチ GD（MiniBatchSize=nTr）は滑らかな損失曲線になりますが、エポック当たりの学習が遅いです。
%[text]    過学習開始点は検証損失が学習損失から乖離し始めるエポックです — このアーキテクチャでは通常 30～60 エポック頃です。
%%
%[text] ## やってみよう 6: 残差バイアス; 最も予測が悪い分子

rng(42);
yTePred = predict(netBase, XTe')';
residuals = yTe - yTePred(:);
logInfo("テスト残差: 平均=%+.3f  標準偏差=%.3f  RMSE=%.3f", ...
    mean(residuals), std(residuals), sqrt(mean(residuals.^2)));
logInfo("残差バイアス: %s", ...
    string(ternary_(abs(mean(residuals)) > 0.2, "バイアスあり", "ほぼ非偏然")));

%[text] 予測が最も悪い分子
[~, worstOrd] = sort(abs(residuals), "descend");
teNames = molNames(find(teIdx));  %#ok<FNDSB>
logInfo("予測が最も悪い上位 3:");
for k = 1:min(3, numel(worstOrd))
    i = worstOrd(k);
    logInfo("  %-28s  actual=%.2f  pred=%.2f  err=%+.2f", ...
        teNames(i), yTe(i), yTePred(i), residuals(i));
end

%[text] 残差 vs 予測値プロット
figure("Name", "A05 残差分析");
scatter(yTePred(:), residuals, 40, [0.3 0.6 0.9], "filled", MarkerFaceAlpha=0.7);
hold on; yline(0, "r--"); hold off;
xlabel("予測 ALogP"); ylabel("残差（実際値 - 予測値）");
title("残差プロット -- テストセット"); grid on;

%[text] 解答: 非ゼロの平均残差は系統的バイアスを示します（ネットワークが平均的に過大または過小予測）。
%[text]    ALogP では高親油性ステロイド（学習に少なすぎる）を過小予測する側面があります。
%[text]    予測が最も悪い分子は、極端な外れ値または学習データで少数だった構造クラスのことが少なくありません。
%%
%[text] ## やってみよう 7: NN vs 線形比較; NN が勝つ場合; GNN の説明

ftTbl  = array2table(Xs(trIdx, :), VariableNames=cellstr(FEAT_NAMES));
ftTbl.ALogP = y(trIdx);
lmFull = fitlm(ftTbl, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));
yTe_lm = predict(lmFull, array2table(XTe, VariableNames=cellstr(FEAT_NAMES)));
rmse_lm = sqrt(mean((yTe - yTe_lm).^2));
rmse_nn = sqrt(mean((yTe - yTePred(:)).^2));
logInfo("テスト RMSE: 線形=%.3f  NN=%.3f", rmse_lm, rmse_nn);
logInfo("NN 優位性: %+.3f  (%s)", rmse_lm - rmse_nn, ...
    string(ternary_(rmse_nn < rmse_lm, "NN 勝利", "線形勝利")));

%[text] 解答: 表形式記述子を用いた ALogP 予測では、NN の線形回帰に対する優位性はわずかで、ゼロまたは負になる場合もあります。
%[text]    記述子は既に関連する化学情報をエンコードしており、線形結合が合理的な近似となります。
%[text]    NN がより強力になるのは次の場合です:
%[text]      (1) 入力が生の分子グラフまたはフィンガープリント（2048 ビット）の場合、
%[text]      (2) 関係が高度に非線形な場合（例: pIC50 の活性クリフ）、
%[text]      (3) 学習データが大規模な場合（n > 10,000）。
%[text]    GNN は分子グラフを直接入力とし（原子をノード、結合をエッジとして）、近傍情報を集約する
%[text]    原子埋め込みを学習します。MPNN や AttentiveFP など GNN は MoleculeNet ベンチマークで SOTA 結果を達成しています。
logInfo("A05 解答完了。");

%[text] ローカルヘルパー関数
function out = ternary_(cond, a, b)
    if cond; out = a; else; out = b; end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
