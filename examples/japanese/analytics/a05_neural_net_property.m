%[text] # A05: ニューラルネットワークによる分子特性予測
%[text] EasyMolKit アナリティクス — レイヤー 3
%[text] 
%[text] QSAR（定量的構造-活性相関）モデルは、従来から線形回帰に依存してきました。記述子と特性の関係がほぼ線形であれば、線形モデルは高速で解釈可能かつ有効です。
%[text] しかし、柔軟性にはコストが伴います。
%[text] ニューラルネットワークはより柔軟な代替手段を提供し、化学者がどの相互作用が重要かを指定することなく、記述子間の非線形相互作用を捉えることができます。
%[text] 一方で、ニューラルネットワークはより多くのデータを必要とし、解釈が難しく、創薬で典型的な小規模データセットでは深刻な過学習を起こす可能性があります。
%[text] この演習では、ALogP（Crippen-Wildman 親油性）をFDA承認薬の8つの構造記述子から予測する小型フィードフォワードニューラルネットワークを構築し、A03の線形回帰モデルと直接比較します。
%[text] 
%[text] **このチュートリアルで学べること**
%[text] - 小規模QSARデータセットがなぜニューラルネットワークにとって挑戦的かを理解する
%[text] - trainnet() を使ってフィードフォワード回帰ネットワークを構築・学習する
%[text] - ドロップアウトと学習率スケジュールで過学習を削減する
%[text] - 学習曲線を解釈する: 収束、過学習、アンダーフィット
%[text] - いつニューラルネットワークより単純な線形モデルを選ぶべきかを知る \
%[text] 
%[text] **前提条件**
%[text] - 推奨: 同じデータセットの文脈のために A03（QSAR 回帰）
%[text] - Deep Learning Toolbox（trainnet、trainingOptions、fullyConnectedLayer）
%[text] - インターネット接続不要 \
%[text] 
%[text] 推定所要時間: 35～50 分 | 実行方法: Ctrl+Enter でセクションを 1 つずつ実行
%[text] 
%[text] **データ:** data/list/fda\_drugs.csv — 200 FDA 承認薬（ChEMBL、CC-BY-SA 3.0）。ALogP 列は Crippen-Wildman 親油性推定値を回帰ターゲットとして使用します。
%[text] 
%[text] **参考文献**
%[text] - Srivastava N, Hinton G, Krizhevsky A, Sutskever I, Salakhutdinov R (2014) Dropout: a simple way to prevent neural networks from overfitting. *J Mach Learn Res* 15:1929-1958.
%[text] - Kingma D & Ba J (2015) Adam: A Method for Stochastic Optimization. arXiv:1412.6980.
%[text] - Ramsundar B et al. (2015) Massively multitask networks for drug discovery. arXiv:1502.02072.
%[text] - Tropsha A (2010) Best practices for QSAR model development, validation, and exploitation. *Mol Inform* 29:476-488. doi:10.1002/minf.201000061. \
%%
%[text] ## セクション 0: セットアップ
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
%[text] メイン実行前に Python/RDKit プロセスをウォームアップします
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;
logSection("A05", "セクション 0: セットアップ", "アナリティクス L3");
%%
%[text] ## セクション 1: FDA 承認薬の読み込みと記述子計算
%[text] 
%[text] セットアップが完了しました。まず、A03 で使用した 8 種類の記述子を計算してデータセットを準備します。
%[text] このセクションを実行することで、A05 が単独で動作することを保証します。
%[text] 
%[text] ### コンセプト: A03 の記述子行列を再利用する
%[text] A03 で使用した 8 つの記述子（MolWt、TPSA、NumHDonors、NumHAcceptors、NumRotatableBonds、RingCount、FractionCSP3、HeavyAtomCount）は、サイズ、極性、トポロジー、形状の 4 つの化学的次元で分子構造をエンコードします。ターゲット変数は ALogP で、これは ChEMBL に格納された Crippen-Wildman 親油性推定値です。
%[text] 
%[text] A03 を既に実行済みであれば、ワークスペースデータを再利用できますが、このセクションを実行することで A05 が完全に自己完結することを保証します。
logSection("A05", "セクション 1: FDA 承認薬の読み込みと記述子計算", "アナリティクス L3");
DATA_FILE  = "data/list/fda_drugs.csv";
FEAT_NAMES = ["MolWt", "TPSA", "NumHDonors", "NumHAcceptors", ...
              "NumRotatableBonds", "RingCount", "FractionCSP3", "HeavyAtomCount"];
N_FEATS    = numel(FEAT_NAMES);

rawTbl = readtable(DATA_FILE, TextType="string");
logInfo("%s から %d 行を読み込みました", DATA_FILE, height(rawTbl));

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
        logWarn("%s をスキップ: %s", rawTbl.Name(k), ME.message);
    end
end
logInfo("%d 個の記述子計算完了", sum(valid));

X        = X_all(valid, :);
y        = y_all(valid);
molNames = rawTbl.Name(valid);
nMols    = sum(valid);

logInfo("データセット: %d / %d 分子（ALogP 範囲: %.2f 〜 %.2f）", ...
    nMols, nRaw, min(y), max(y));
%[text] **💡 観察ポイント 1 — パラメータ対サンプル比と外れ値の確認**
%[text] 特徴行列（nMols x N\_FEATS）は、セクション 4 で定義されるネットワークの学習可能パラメータ数と比較してどのくらいの規模か確認しましょう。
%[text] パラメータ数が学習サンプル数を超えると過学習のリスクが高まります。ここで、その比率を確認してみましょう。
%[text] |ALogP| \> 7 の分子が存在するか確認しましょう。次のコードで特定できます:
%[text] \[~,ord\] = sort(abs(y),"descend");
%[text] table(molNames(ord(1:5)), y(ord(1:5)), VariableNames=\["Name","ALogP"\])
%[text] これらの外れ値が MSE ベースのニューラルネットワーク学習に不均衡な影響を与えるか確認しましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: ALogP 分布の探索
%[text] 
%[text] 記述子の計算が完了しました。モデルを構築する前に、ターゲット変数である ALogP（LogP の変種）の分布を確認しましょう。
%[text] 外れ値の有無や分布の形状を把握することで、後の学習結果を正しく解釈する助けになります。
%[text] 
%[text] ### コンセプト: 学習前にターゲットを調べる
%[text] ニューラルネットワーク回帰器は平均二乗誤差（MSE）を最小化します。MSE は大きな誤差を二乗で罰するため、y の外れ値が学習済み重みに不均衡な影響を与えます。ネットワークは少数の極端な例の誤差を減らすために、ほとんどの分子の精度を犠牲にする可能性があります。
%[text] 
%[text] 学習前にターゲット分布を可視化することで、外れ値の影響を減らすための y 変換（クリッピングやウィンソライゼーションなど）を適用するかどうかを判断できます。この演習では生の ALogP を使用します。
logSection("A05", "セクション 2: ALogP 分布の探索", "アナリティクス L3");
figure("Name", "A05 ALogP 分布");
histogram(y, 20, FaceColor=[0.2 0.6 0.8]);
xlabel("ALogP"); ylabel("件数");
title("ALogP 分布 -- FDA 承認薬");
xline(mean(y),   "r--", sprintf("平均=%.2f",   mean(y)),   LabelHorizontalAlignment="right",  LabelOrientation="horizontal");
xline(median(y), "g--", sprintf("中央値=%.2f", median(y)), LabelHorizontalAlignment="left",   LabelOrientation="horizontal");
grid on;

logInfo("ALogP: 平均=%.2f、標準偏差=%.2f、歪度=%.2f", mean(y), std(y), skewness(y));
%[text] **💡 観察ポイント 2 — ALogP 分布の正規性検定とウィンソライズ**
%[text] ALogP はほぼガウス分布かどうかを確認しましょう。次のコードでテストします:
%[text] \[h, p\] = kstest((y - mean(y)) / std(y))
%[text] p \< 0.05 なら、分布は有意に非ガウスであると判断できます。
%[text] この結果はニューラルネットワーク学習において重要でしょうか（線形回帰と比較して）?
%[text] （ヒント: OLS はガウス残差を仮定しますが、MSE 損失を用いる NN はそうではありません）
%[text] 5 パーセンタイルと 95 パーセンタイルでターゲットをウィンソライズしてみましょう:
%[text] y\_clip = min(max(y, prctile(y,5)), prctile(y,95));
%[text] セクション 3〜7 で y を y\_clip に置き換えて R² を比較してみましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: 学習 / 検証 / テスト分割と特徴量標準化
%[text] 
%[text] データの分布を確認しました。次に、ニューラルネットワーク学習に必要なデータの三方向分割を行います。
%[text] 検証セットは学習中の過学習を監視するために、テストセットは最終評価のために使用します。
%[text] 
%[text] ### コンセプト: ニューラルネットワーク学習のための三方向分割
%[text] A03で使用したk分割交差検証とは異なり、ニューラルネットワークでは専用の検証セットが必要です。（a）過学習をリアルタイムで監視するため / （b）早期停止をトリガーするため / （c）テストセットを汚染せずにハイパーパラメータをチューニングするためです。
%[text] 
%[text] **規約:** 70% 学習 / 15% 検証 / 15% テスト。200分子の場合、約140 / 30 / 30となります。
%[text] **データリークルール:** 標準化のための統計量（平均、標準偏差）は学習データのみで計算し、検証・テストセットに適用します。
logSection("A05", "セクション 3: 学習 / 検証 / テスト分割と特徴量標準化", "アナリティクス L3");
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
%[text] 学習セットの統計量のみを用いて特徴量を標準化します。
Xmean = mean(Xtr, 1);
Xstd  = std(Xtr,  0, 1);
Xstd(Xstd < 1e-12) = 1;  % guard: avoid /0 for constant descriptors

Xtr_s  = (Xtr  - Xmean) ./ Xstd;
Xval_s = (Xval - Xmean) ./ Xstd;
Xte_s  = (Xte  - Xmean) ./ Xstd;

logInfo("分割: 学習=%d / 検証=%d / テスト=%d（合計=%d）", ...
    size(Xtr, 1), size(Xval, 1), size(Xte, 1), nMols);
%[text] **💡 観察ポイント 3 — ランダムシードと分割の信頼性**
%[text] なぜ rng(42) が必要なのか確認しましょう。これを削除して、セクション 3〜7 を5回再実行し、
%[text] 各回のテスト R² を記録しましょう。分散はどの程度でしょうか？
%[text] これは約200分子のデータセットでの単一テスト分割の信頼性について何を示しているでしょうか？
%[text] 代替案として、ネストされた交差検証（外ループ: テスト分割; 内ループ: チューニング用検証分割）があります。
%[text] なぜこれが計算コストは高いものの、より偏りの少ない性能推定を提供するのかを説明しましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: ニューラルネットワークアーキテクチャの定義
%[text] 
%[text] データの準備が整いました。次にネットワークの構造（アーキテクチャ）を設計します。
%[text] 小規模データセットに適した、ドロップアウトによる正則化を施した2層の浅いネットワークを使用します。
%[text] 
%[text] ### コンセプト: 小規模データセット回帰のためのフィードフォワードネットワーク設計
%[text] ネットワークは8つの標準化された記述子を2つの隠れ層を通して単一のALogP（常用対数の推定値）予測にマッピングします: `h1 = ReLU(W1×x + b1)` → Dropout(0.20) → `h2 = ReLU(W2×h1d + b2)` → `yhat = W3×h2 + b3`（線形出力）。学習可能パラメータの合計は2689です。
%[text] 
%[text] **ReLU** — 計算コストが低く、勾配消失を避け、スパースな活性化を生成します。**ドロップアウト** (p=0.20) — 各学習フォワードパスでニューロンの20%をランダムにゼロ化し、共適応を防ぎます。約140の学習分子と2689のパラメータでは、パラメータ対サンプル比は約19です。ドロップアウトは不可欠です。
%[text] 
logSection("A05", "セクション 4: ニューラルネットワークアーキテクチャの定義", "アナリティクス L3");
layers = [
    featureInputLayer(N_FEATS, Normalization="none")
    fullyConnectedLayer(64)
    reluLayer
    dropoutLayer(0.20)
    fullyConnectedLayer(32)
    reluLayer
    fullyConnectedLayer(1)];

net = dlnetwork(layers);
%[text] 学習可能パラメータを数えます
numParams = 0;
for k = 1:height(net.Learnables)
    numParams = numParams + numel(extractdata(net.Learnables.Value{k}));
end
logInfo("アーキテクチャ: FC(64)-ReLU-Drop(0.2)-FC(32)-ReLU-FC(1)");
logInfo("学習可能パラメータ: %d  （パラメータ / 学習サンプル比 = %.1f）", ...
    numParams, numParams / size(Xtr_s, 1));
%[text] **💡 観察ポイント 4 — アーキテクチャの変更と効果の比較**
%[text] 出力層の前に3番目の隠れ層（fullyConnectedLayer(16) + reluLayer）を追加して、効果を確認しましょう。
%[text] より深いネットワークはテストR2を改善するでしょうか? コスト（学習時間、過学習リスク）はどう変化するでしょうか?
%[text] reluLayerをtanhLayerに置き換えて、影響を確認しましょう。飽和活性化関数はこの小規模データセットでの学習速度と最終性能にどのように影響するでしょうか?
%[text] （ヒント: tanhの飽和領域近くでは勾配がゼロに近づきます — これは勾配消失と呼ばれる問題です）
%[text] 0\.20の代わりにdropoutLayer(0.50)を設定して、影響を確認しましょう。
%[text] 学習曲線（セクション6）で過学習がより多く見られるでしょうか、それとも少なくなるでしょうか? 理由を考えてみましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: ネットワークの学習
%[text] 
%[text] アーキテクチャの定義が完了しました。次にネットワークの学習を行います。
%[text] 学習中に表示される進捗ウィンドウで、学習損失と検証損失の推移を確認しましょう。
%[text] 
%[text] ### コンセプト: Adam オプティマイザとピースワイズ学習率スケジュール
%[text] Adam は過去の勾配の指数減衰平均を用いて、各パラメータの学習率を適応的に調整します。ピースワイズスケジュールでは、エポック 80 で学習率を半分にします: エポック 1-80 は lr=1e-3、エポック 81-200 は lr=5e-4 です。
%[text] 
%[text] **ミニバッチ SGD:** 各 Adam ステップでは、32 分子のランダムサブセットを使用します。140 学習分子で MiniBatchSize=32 の場合、1 エポックは 4 イテレーションとなります。学習進捗ウィンドウが自動的に表示されます。
%[text] 
logSection("A05", "セクション 5: ネットワークの学習", "アナリティクス L3");
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

logInfo("200 エポック学習中（約 30〜60 秒）...");
[net, trainInfo] = trainnet(Xtr_s, ytr, net, "mse", opts);
%[text] バージョン互換の損失抽出方法について説明します。
%[text] R2024a: trainInfo.TrainingLoss / trainInfo.ValidationLoss（直接配列）
%[text] R2025b: TrainingHistory.Loss / ValidationHistory.Loss（分離テーブル、列名 "Loss"）
%[text] R2025a: TrainingHistory（TrainingLoss, ValidationLoss 列のテーブル）
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
logInfo("学習完了。最終検証 MSE=%.4f（RMSE=%.3f）", ...
    finalValLoss, sqrt(max(finalValLoss, 0)));
%%
%[text] ## セクション 6: 学習曲線分析
%[text] 
%[text] 学習が完了しました。次に、学習曲線を詳しく分析してみましょう。
%[text] 学習損失と検証損失の差を見ることで、過学習が発生しているかどうかを診断できます。
%[text] 
%[text] コンセプト: 損失曲線を用いて学習の動きを診断する
%[text] ### コンセプト: 学習曲線の読み方
%[text] 学習曲線は、学習エポックに対する損失の変化をプロットします。**良好なフィット:** 学習と検証の平均二乗誤差がともに減少し、近い値に収束する。**過学習:** 学習損失は下がり続けるが、検証損失は横ばいまたは上昇する。**アンダーフィット:** 両方の損失が高いままで、進展が遅い。
%[text] 
logSection("A05", "セクション 6: 学習曲線分析", "アナリティクス L3");
% itersPerEpoch を実際の学習履歴から算出（計算式比較で正確）
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
%[text] 解釈を容易にするために、MSE を RMSE に変換します（ALogP と同じ単位）。
trRmse  = sqrt(max(trLossAll,  0));
valRmse = sqrt(max(valLossAll, 0));
%[text] 1 エポック移動平均を用いて学習曲線を平滑化します。
winLen    = max(1, itersPerEpoch);
trSmoothed = movmean(trRmse, winLen, "Endpoints", "fill");

figure("Name", "A05 学習曲線");
plot(epochsTr,  trRmse,    "Color", [0.7 0.8 1.0], LineWidth=0.8, ...
    DisplayName="学習（生）"); hold on;
plot(epochsTr,  trSmoothed, "b-", LineWidth=2.0, DisplayName="学習（平滑化）");
plot(epochsVal, valRmse,    "r--", LineWidth=2.0, DisplayName="検証");
xlabel("エポック"); ylabel("RMSE（ALogP 単位）");
title("A05: 学習曲線");
legend(Location="northeast");
grid on;
%[text] 最良の検証エポックをマークします。
[minValRmse, minValIdx] = min(valRmse);
xline(epochsVal(minValIdx), "--k", ...
    sprintf("最良検証 RMSE=%.3f（エポック %.0f）", minValRmse, epochsVal(minValIdx)), ...
    LabelHorizontalAlignment="left", LabelVerticalAlignment="bottom", ...
    LabelOrientation="horizontal", HandleVisibility="off");

logInfo("最良検証 RMSE=%.3f（エポック %.0f）", minValRmse, epochsVal(minValIdx));
%[text] **💡 観察ポイント 5 — 学習曲線の診断と早期停止**
%[text] 学習曲線は過学習を示していますか？（学習損失が下がり続ける一方で、検証損失が上昇する）
%[text] 学習と検証の RMSE の差が広がり始めるのは何エポック目か確認しましょう。
%[text] セクション 5 の trainingOptions 呼び出しに早期停止を追加してみましょう。
%[text] ValidationPatience=20 （20 回チェックで改善がなければ停止）
%[text] 早期停止はテスト R2 を改善するか、悪化させるか、それとも変わらないか確認しましょう。
%[text] MiniBatchSize=nMols（フル勾配降下法）に設定するとどうなるか考えてみましょう。
%[text] 期待される学習曲線の形状をミニバッチと比較して確認しましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 7: テストセット評価
%[text] 
%[text] 学習曲線を用いてモデルの動作を確認しました。最後に、テストセットを用いてモデルの性能を公正に評価します。
%[text] テストセットはここで初めて使用します。この結果が最終的な性能指標となります。
%[text] 
%[text] ### コンセプト: 一回限りのテストセット評価
%[text] テストセットは、全ての設計決定（アーキテクチャ、ハイパーパラメータ）が確定した後に、正確に1回だけ評価します。**回帰の主要指標:** $R^2 = 1 - SS\_{res}/SS\_{tot}$ / RMSE（yと同じ単位）/ MAE（外れ値に対して頑健）。実際値 vs 予測値の散布図で、y=xの対角線上の点が完全予測を示します。
%[text] 
logSection("A05", "セクション 7: テストセット評価", "アナリティクス L3");
yPred = double(predict(net, Xte_s));   % (nTest x 1) predicted ALogP
res   = yte - yPred;
rmse  = sqrt(mean(res.^2));
mae   = mean(abs(res));
r2    = 1 - sum(res.^2) / sum((yte - mean(yte)).^2);

logInfo("テストセット: R2=%.3f  RMSE=%.3f  MAE=%.3f  (n=%d)", r2, rmse, mae, numel(yte));

figure("Name", "A05 テスト予測");
tiledlayout(1, 2);

nexttile;
scatter(yte, yPred, 50, "filled", MarkerFaceAlpha=0.7, ...
    MarkerFaceColor=[0.2 0.5 0.9]);
hold on;
lim = [min([yte; yPred]) - 0.5, max([yte; yPred]) + 0.5];
plot(lim, lim, "k--", LineWidth=1.5);
xlabel("実際の ALogP"); ylabel("予測 ALogP");
title(sprintf("実際値 vs 予測値  (R^2=%.3f)", r2));
grid on;

nexttile;
histogram(res, 15, FaceColor=[0.9 0.4 0.3]);
xline(0, "--k", LineWidth=1.5);
xlabel("残差（実際値 - 予測値）"); ylabel("件数");
title(sprintf("残差  (RMSE=%.3f, MAE=%.3f)", rmse, mae));
grid on;
%[text] **💡 観察ポイント 6 — 残差バイアスと予測が難しい分子の特定**
%[text] 残差はほぼゼロを中心としていますか？ 非ゼロの平均は系統的バイアスを示します。
%[text] （ネットワークが一貫して過大または過小予測していることを意味します）。
%[text] 小規模データセットのニューラルネットワークで、これが起こる原因は何でしょうか？
%[text] 最も予測が悪い3つの分子（最大 |残差|）を特定しましょう。
%[text] それらは構造的外れ値ですか、それとも共通のスキャフォールドを共有していますか？
%[text] \[~,ord\] = sort(abs(res), "descend");
%[text] tstNames = molNames(tmpIdx);  % tmpIdx サブセットに絞り込む
%[text] table(tstNames(tstIdx)(ord(1:3)), yte(ord(1:3)), yPred(ord(1:3)), ...
%[text] VariableNames=\["Name","実際値","予測値"\])
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 8: 線形ベースライン比較と考察
%[text] 
%[text] ニューラルネットワークの評価が完了しました。最後に、同じデータを用いて線形回帰と比較します。
%[text] 「より複雑なモデルが常に優れているわけではない」というQSARの重要な教訓を確認しましょう。
%[text] 
%[text] ### コンセプト: オッカムの剃刀 — 最もシンプルで適切なモデルを選ぶ
%[text] QSARでは「最もシンプルで適切なモデル」は通常、線形回帰です。ニューラルネットワークに置き換えるのは、(a) 記述子あたりのデータがN\>1000と大きい場合、(b) 非線形性を期待する理論的根拠がある場合、(c) 交差検証済みのニューラルネットワークの性能が線形を大幅に上回る場合に限られます。
%[text] 
%[text] 100〜300分子の薬物データセットでは、線形回帰、ランダムフォレスト、ガウス過程モデルが通常、ニューラルネットワークに匹敵するか上回ります。GNNこそがケモインフォマティクスで深層学習の真価が発揮される領域です。この比較では公平を期すために同じ学習/テスト分割を使用します。
%[text] 
logSection("A05", "セクション 8: 線形ベースライン比較と考察", "アナリティクス L3");
ftTbl    = array2table(Xtr_s, VariableNames=cellstr(FEAT_NAMES));
ftTbl.ALogP = ytr;
lmMdl    = fitlm(ftTbl, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));

teTbl    = array2table(Xte_s, VariableNames=cellstr(FEAT_NAMES));
yLinPred = predict(lmMdl, teTbl);
r2Lin    = 1 - sum((yte - yLinPred).^2) / sum((yte - mean(yte)).^2);
rmseLin  = sqrt(mean((yte - yLinPred).^2));
maeLin   = mean(abs(yte - yLinPred));

logInfo("線形回帰:           R2=%.3f  RMSE=%.3f  MAE=%.3f", r2Lin, rmseLin, maeLin);
logInfo("ニューラルネット:   R2=%.3f  RMSE=%.3f  MAE=%.3f", r2, rmse, mae);
logInfo("NN の優位性:  dR2=%.3f  dRMSE=%+.3f  dMAE=%+.3f", ...
    r2 - r2Lin, rmseLin - rmse, maeLin - mae);
%[text] 2つのモデルを並べて比較します。
figure("Name", "A05 モデル比較");
set(gcf, "Position", [100 100 1100 400]);
tiledlayout(1, 3);

metrics  = [r2Lin, r2; rmseLin, rmse; maeLin, mae];
ylabels  = ["R^2（高いほど良い）", "RMSE（低いほど良い）", "MAE（低いほど良い）"];
for m = 1:3
    nexttile;
    bar(metrics(m, :), FaceColor="flat", CData=[0.4 0.6 0.9; 0.9 0.4 0.3]);
    set(gca, "XTickLabel", ["線形", "ニューラルネット"]);
    ylabel(ylabels(m)); title(ylabels(m)); grid on;
    ylim([0, max(metrics(m, :)) * 1.3]);
end
sgtitle("A05: 線形回帰 vs ニューラルネットワーク（同一学習・テスト分割）");
%[text] **💡 観察ポイント 7 — ニューラルネットワークと線形モデルの比較、およびGNNへの展望**
%[text] ニューラルネットワークはこのデータセットで線形回帰を上回っているか確認しましょう。
%[text] 線形モデルが優れている場合、その理由を3つ具体的に挙げてみましょう。
%[text] ニューラルネットワークが明確な優位性を得るには何が変わる必要があるか考えてみましょう。
%[text] データセットサイズ（N）、特徴量タイプ（フィンガープリント vs 記述子）、特性タイプ（活性 vs 物理化学的特性）について検討してみましょう。
%[text] グラフニューラルネットワーク（GNN）は分子を固定長記述子ベクトルではなく、グラフとして表現します。これは特徴量エンジニアリングのステップをどのように変えるでしょうか？
%[text] そして、なぜより大きなデータセットで良い汎化が可能になるのか説明してみましょう。
%[text] （参照: Gilmer et al. 2017, Neural Message Passing）
% ... （ここにコードを書いてみましょう）

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
