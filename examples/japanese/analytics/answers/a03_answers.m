%[text] # A03 解答: QSAR 回帰（LogP 予測）
%[text] a03_qsar_regression.m の「やってみよう」演習の参照解答です。
%[text] まず a03_qsar_regression.m を実行してワークスペース変数を準備してください。
%[text] なおこのファイルはスタンドアロンでも実行でき、各セクションを独立して動かすこともできます。
addpath(genpath("src"));
emk.setup.initPython();
logInfo("A03 解答: セットアップ完了");

%[text] ワークスペースを再構築（スタンドアロン実行用）
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
    logProgress(j, nMols, "記述子を計算中");
end

Xmean = mean(X, 1);
Xstd  = std(X, 0, 1);
Xstd(Xstd < 1e-12) = 1;
Xs = (X - Xmean) ./ Xstd;
%%
%[text] ## やってみよう 1: 各記述子と ALogP の相関; 多重共線性チェック

corrVec = corr(X, y);
logInfo("ALogP とのピアソン相関:");
for d = 1:N_FEATS
    logInfo("  %-22s  r = %+.3f", FEAT_NAMES(d), corrVec(d));
end
[~, sortCorr] = sort(abs(corrVec), "descend");
logInfo("最強相関記述子: %s (r=%+.3f)", ...
    FEAT_NAMES(sortCorr(1)), corrVec(sortCorr(1)));

%[text] ペアワイズ共線性の確認
C = corr(X);
C_off = C - diag(diag(C));
[maxC, linIdx] = max(abs(C_off(:)));
[r, c] = ind2sub(size(C), linIdx);
logInfo("最も多重共線性が高いペア: %s と %s  (r = %.3f)", ...
    FEAT_NAMES(r), FEAT_NAMES(c), C(r, c));

%[text] 解答: TPSA（負の相関 r ~-0.6）と NumHAcceptors が ALogP と最も強い相関を示します。
%[text]    極性表面積は親油性を低下させることが化学的に自然です。
%[text]    MolWt と HeavyAtomCount は強い共線性（r > 0.95）を持ちます。
%[text]    共線性が高いと最小二乗法（OLS）の係数の分散が大きくなり、個々の係数推定が不安定になります
%[text]    （符号の反転も起こりえます）。
%%
%[text] ## やってみよう 2: 高 ALogP 外れ値を特定; 分布形状を確認

[sorted_y, ord] = sort(y, "descend");
logInfo("最も親油性の高い上位 5 分子:");
for k = 1:5
    logInfo("  %d. %-30s  ALogP = %.2f", k, molNames(ord(k)), sorted_y(k));
end

logInfo("ALogP: 平均=%.2f  中央値=%.2f  歪度=%.2f", ...
    mean(y), median(y), skewness(y));
logInfo("ALogP > 6 の分子数: %d", sum(y > 6));

%[text] 解答: ALogP > 6 の分子は高親油性化合物であり（ステロイドなど脂質様化合物が代表的）、
%[text]    OLS はすべての残差を平等にペナルティとして扱うため、これらの外れ値が
%[text]    大多数の分子に対する決定超平面を引き寄せ、全体のフィット品質を低下させる場合があります。
%[text]    実用上の対処としては、外れ値を除外または重みを軽くするか、ロバスト回帰（robustfit()）を利用する方法があります。
%[text]    歪度が 0 に近い場合、LogP 分布は対称に近く、対数変換なしでモデリングできる状態です。
%%
%[text] ## やってみよう 3: 線形回帰 — 大きな係数; 非有意な特徴量

ftTbl = array2table(Xs, VariableNames=cellstr(FEAT_NAMES));
ftTbl.ALogP = y;
lmModel = fitlm(ftTbl, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));

coeffTbl = lmModel.Coefficients;
logInfo("線形モデル係数（標準化特徴量）:");
disp(coeffTbl);

%[text] p > 0.05 の統計的非有意特徴量を抽出
pVals   = coeffTbl.pValue(2:end);   % 切片をスキップ
names   = string(coeffTbl.Properties.RowNames(2:end));
insig   = names(pVals > 0.05);
logInfo("統計的に非有意（p>0.05）: %s", strjoin(insig, ", "));

%[text] HeavyAtomCount を削除して再フィット
dropIdx  = strcmp(cellstr(FEAT_NAMES), "HeavyAtomCount");
featRed  = FEAT_NAMES(~dropIdx);
trTblRed = array2table(Xs(:, ~dropIdx), VariableNames=cellstr(featRed));
trTblRed.ALogP = y;
lm2 = fitlm(trTblRed, "ALogP ~ " + strjoin(featRed, " + "));
logInfo("フルモデル RMSE: %.3f   削減モデル RMSE: %.3f", ...
    lmModel.RMSE, lm2.RMSE);

%[text] 大きな残差を持つ分子を確認
residuals = lmModel.Residuals.Raw;
[~, bigResOrd] = sort(abs(residuals), "descend");
logInfo("最大 |残差| を持つ分子:");
for k = 1:5
    logInfo("  %-30s  residual = %+.2f", molNames(bigResOrd(k)), residuals(bigResOrd(k)));
end

%[text] 解答: 正の係数を持つ特徴量: HeavyAtomCount、NumRotatableBonds
%[text]    （原子数・結合数が多いほど平均的に親油性が高い傾向がある）。
%[text]    負の係数: TPSA、NumHDonors、NumHAcceptors（極性基が親油性を低下させる）。
%[text]    HeavyAtomCount （MolWt と共線性）を削除して再フィットしても RMSE はほぼ変わらず、
%[text]    特徴量の冗長性が確認できます。
%%
%[text] ## やってみよう 4: RF 特徴量重要度; MinLeafSize の影響; OOB 収束

rfModel = fitrensemble(Xs, y, ...
    Method="Bag", NumLearningCycles=200, ...
    Learners=templateTree(MinLeafSize=3));

importances = predictorImportance(rfModel);
[sortedImp, sortOrd] = sort(importances, "descend");
logInfo("RF 特徴量重要度ランキング:");
for d = 1:N_FEATS
    logInfo("  %d. %-22s  MDI = %.4f", d, FEAT_NAMES(sortOrd(d)), sortedImp(d));
end

%[text] MinLeafSize の比較
for mls = [1, 10, 20]
    rf_mls = fitrensemble(Xs, y, Method="Bag", NumLearningCycles=100, ...
        Learners=templateTree(MinLeafSize=mls));
    yPredTrain = predict(rf_mls, Xs);
    rmseTrain  = sqrt(mean((y - yPredTrain).^2));
    logInfo("MinLeafSize=%2d -> 学習 RMSE = %.3f", mls, rmseTrain);
end

%[text] OOB 誤差の収束履歴
oobRmse = oobLoss(rfModel, "Mode", "cumulative");
figure("Name", "A03 OOB RMSE 収束");
plot(oobRmse, "Color", [0.4 0.7 0.4], LineWidth=1.5);
xlabel("木の本数"); ylabel("OOB RMSE");
title("RF OOB 誤差収束");
grid on;
logInfo("OOB RMSE は約 %d 本の木で安定", ...
    find(abs(diff(oobRmse)) < 0.001, 1, "first"));

%[text] 解答: RF は最重要特徴量として極性表面積（TPSA）または HeavyAtomCount を選ぶことが多いです。
%[text]    これらは ALogP の変動の 2 軸（極性と大きさ）を両方捕えるためです。
%[text]    MinLeafSize=1 は深い木を学習: 学習 RMSE は 0 に近づきます（過学習）。
%[text]    MinLeafSize=20 は木を強く刈り込みます: 学習・テスト共に RMSE が高くなります。
%[text]    MinLeafSize=3−5 が QSAR データセットでの典型的な妥協点です。
%[text]    OOB 誤差は通常 150−200 本前後に平坦に達します。
%%
%[text] ## やってみよう 5: 5 分割 CV 比較; 過学習ギャップ; 10 分割感度

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
    logProgress(fold, K_FOLDS, "CV フォールド");
end

r2cv_lm  = 1 - sum((y-yCV_lm).^2) / sum((y-mean(y)).^2);
rmse_lm  = sqrt(mean((y-yCV_lm).^2));
r2cv_rf  = 1 - sum((y-yCV_rf).^2) / sum((y-mean(y)).^2);
rmse_rf  = sqrt(mean((y-yCV_rf).^2));

yTr_rf   = predict(rfModel, Xs);
rmse_rfTr = sqrt(mean((y-yTr_rf).^2));

logInfo("5 分割 CV — 線形モデル: R2=%.3f  RMSE=%.3f", r2cv_lm, rmse_lm);
logInfo("5 分割 CV — Random Forest: R2=%.3f  RMSE=%.3f", r2cv_rf, rmse_rf);
logInfo("RF 過学習ギャップ: 学習 RMSE=%.3f  CV RMSE=%.3f  ギャップ=%.3f", ...
    rmse_rfTr, rmse_rf, rmse_rf - rmse_rfTr);

%[text] 10 分割 CV との比較

rng(42);
cv10 = cvpartition(nMols, KFold=10);
yCV_rf10 = nan(nMols, 1);
for fold = 1:10
    trainIdx = training(cv10, fold);
    testIdx  = test(cv10, fold);
    rf10 = fitrensemble(Xs(trainIdx, :), y(trainIdx), ...
        Method="Bag", NumLearningCycles=200, Learners=templateTree(MinLeafSize=3));
    yCV_rf10(testIdx) = predict(rf10, Xs(testIdx, :));
    logProgress(fold, 10, "10 分割 CV");
end
r2cv_rf10 = 1 - sum((y-yCV_rf10).^2) / sum((y-mean(y)).^2);
rmse_rf10 = sqrt(mean((y-yCV_rf10).^2));
logInfo("10 分割 CV — Random Forest: R2=%.3f  RMSE=%.3f", r2cv_rf10, rmse_rf10);

%[text] 解答: 200 種の FDA 承認薬で RF は汎化 CV R2 > 0.8 を達成することが常で、線形モデルは ~0.6〜0.7 程度です。
%[text]    RF の過学習ギャップ（CV RMSE - 学習 RMSE）は通常 0.3 以上で線形モデルより大きく、
%[text]    RF の方がより強く学習データに適合することを示します。
%[text]    n=200 では 5 分割 vs 10 分割 CV の R2 は近似しますが、
%[text]    10 分割の方がバイアスが少ない一方で、CV 推定値の分散は大きくなります。
%%
%[text] ## やってみよう 6: RDKit Crippen との予測比較; ドメインチェック; 再学習の根拠
%[text] エタノールの QSAR 予測値と RDKit Crippen 値を比較
molEt = emk.mol.fromSmiles("CCO");
sCrippen = emk.descriptor.calculate(molEt, "LogP");
logInfo("エタノールの RDKit Crippen LogP: %.2f", sCrippen.LogP);

%[text] 最終モデル — 全データで再学習
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
logInfo("最終モデル: %s（全データ再学習）", modelLabel);

%[text] エタノールのドメインチェック
sEt = emk.descriptor.calculate(molEt, FEAT_NAMES);
xEt = zeros(1, N_FEATS);
for d = 1:N_FEATS; xEt(d) = sEt.(FEAT_NAMES(d)); end
xEtS = (xEt - Xmean) ./ Xstd;
inDomain = all(xEtS >= min(Xs,[],1) - 2 & xEtS <= max(Xs,[],1) + 2);
logInfo("エタノール ドメインチェック: %s", string(ternary_(inDomain, "ドメイン内", "ドメイン外")));
logInfo("エタノール QSAR 予測: %.2f  (Crippen: %.2f)", predictFn(xEtS), sCrippen.LogP);

%[text] 大分子（ペプチド）のドメインチェック
molPep = emk.mol.fromSmiles("NCC(=O)NCC(=O)O");
sPep   = emk.descriptor.calculate(molPep, FEAT_NAMES);
xPep   = zeros(1, N_FEATS);
for d = 1:N_FEATS; xPep(d) = sPep.(FEAT_NAMES(d)); end
xPepS  = (xPep - Xmean) ./ Xstd;
inDomPep = all(xPepS >= min(Xs,[],1) - 2 & xPepS <= max(Xs,[],1) + 2);
logInfo("Ala-Gly ペプチド ドメイン: %s", string(ternary_(inDomPep, "ドメイン内", "ドメイン外")));

%[text] 解答: エタノールに対する QSAR 予測値と Crippen 値は近い値（およそ -0.3 ~ -0.7）になるはずです。
%[text]    これは学習ターゲットの ALogP 自体が Crippen 法で計算された値なためです。
%[text]    エタノールは小さな単純分子で、FDA 承認薬（学習データ）と化学クラスが異なりますが、
%[text]    記述子空間内に入る可能性があります。
%[text]    ジペプチドは HeavyAtomCount または NumHDonors でドメイン外と判定される可能性が高いです。
%[text]    ドメイン外の予測を使用すると大きな誤差のリスクがあります。
%[text]    CV では各フォールドで 20% のデータを保留するため、本番モデルは 100% のデータで再学習することでバイアスを減らします。
logInfo("A03 解答完了。");

function out = ternary_(cond, a, b)
    if cond; out = a; else; out = b; end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
