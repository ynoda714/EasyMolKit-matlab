%[text] # A04: 薬物分類 — FDA 承認薬 vs 日用化学品
%[text] EasyMolKit アナリティクス — レイヤー 3
%[text] 
%[text] 「この白い粉末は医薬品に見えますか？」—法医学の現場では、分子構造だけで即座に判断を求められることがあります。
%[text] FDA 承認薬 200 種と日用化学品 30 種を 2048 ビットのフィンガープリントに変換し、SVM とランダムフォレストで「薬らしさ」を学習させます。
%[text] 精度、F1 スコア、ROC/AUC を比較し、クラス不均衡データにおける適切な評価指標の選び方を学びます。
%[text] このスクリプトでは、p \>\> n レジームの前処理手順と分類モデルの実践的な評価方法を体験します。
%[text] 
%[text] **このチュートリアルで学べること**
%[text] - フィンガープリントベースの機械学習における「p \>\> n」問題を理解する
%[text] - 教師あり分類の前処理として PCA を適用する方法を学ぶ
%[text] - `fitcsvm()` と `fitcensemble()` を用いて二値分子分類を行う
%[text] - 混同行列、精度、再現率、F1 スコアを読み解き、解釈する
%[text] - ROC 曲線と AUC を閾値に依存しない性能指標として解釈する
%[text] - クラス不均衡が精度ベースの評価に与える影響を認識する \
%[text] 
%[text] **前提条件**
%[text] - F03（フィンガープリント）と F04（類似度）の完了が必要です
%[text] - 推奨: A01（PCA）と A02（クラスタリング）でコンテキストを把握すること
%[text] - Statistics and Machine Learning Toolbox（`fitcsvm`、`fitcensemble`、`perfcurve`）が必要です
%[text] - インターネット接続は不要です \
%[text] 
%[text] 推定所要時間: 30〜45 分 | 実行方法: Ctrl+Enter でセクションを 1 つずつ実行
%[text] 
%[text] **データ:** `data/list/fda_drugs.csv`（200 FDA 承認薬、ChEMBL CC-BY-SA 3.0）、`data/list/everyday_chemicals.csv`（30 種の家庭用化学品、PubChem CC0）
%[text] 
%[text] **ラベル品質に関する注意**
%[text] データセットの出所をクラスラベルとして使用します（薬 vs 非薬）。
%[text] アスピリン、カフェイン、イブプロフェンなど一部の日用化学品は実際に医薬品でもあります—これは意図的なラベルノイズです。
%[text] 実際の QSAR 分類では、実験的に測定された活性値を使用します。
%[text] 不完全なラベルは、達成可能な最大 AUC を制限します（セクション 6 参照）。
%[text] 
%[text] **参考文献**
%[text] - Cortes & Vapnik (1995) Support-vector networks. *Machine Learning* 20:273-297.
%[text] - Breiman (2001) Random forests. *Machine Learning* 45:5-32.
%[text] - Fawcett (2006) An introduction to ROC analysis. *Pattern Recognition Letters* 27:861-874.
%[text] - Rogers & Hahn (2010) Extended-connectivity fingerprints. *J Chem Inf Model* 50:742-754. \
%%
%[text] ## セクション 0: 環境セットアップ
logSection("A04", "セクション 0: セットアップ", "アナリティクス L3");
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
%[text] メイン処理を実行する前に、Python/RDKit プロセスをウォームアップします
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;
%%
%[text] ## セクション 1: データセット読み込みとフィンガープリント計算
%[text] 
%[text] セットアップが完了しました。まず、2種類のデータセットを読み込み、各分子を2048ビットのフィンガープリントに変換します。
%[text] このビットベクトルが機械学習の入力特徴量として使用されます。
%[text] 
%[text] ### コンセプト: 機械学習のためのバイナリ特徴ベクトルとしてのフィンガープリント
%[text] Morgan / ECFP4 フィンガープリントは2048ビットのバイナリベクトルです。各ビットは円形部分構造（各重原子から半径2、すなわち中心原子から2結合先まで）の有無をエンコードします。ECFP4の「4」は直径（= 2×半径）を示します。
%[text] 
%[text] **機械学習における役割:** 各ビットが1つの特徴量となり（分子あたり2048特徴量）、分子は0/1特徴行列の行として表現されます。分類器は2つのクラスを区別するビットパターンを学習します。
%[text] 
%[text] **課題 — 「p \>\> n」レジーム:** 約230分子（n）に対して2048特徴量（p）では、特徴量がサンプルより多くなります。これにより、線形分離可能性はほぼ保証されますが、過学習のリスクがあります。多くのビットがほぼ一定で情報がないため、SVMの超平面推定が不安定になります。解決策として、PCAによる次元削減が考えられます（セクション2）。
%[text] 
%[text] **ラベル規約:** 1 = FDA承認医薬品化合物（ChEMBL由来）/ 0 = 日用家庭用化学品（PubChem由来）
logSection("A04", "セクション 1: データセット読み込みとフィンガープリント計算", "アナリティクス L3");
FP_RADIUS = 2;
FP_NBITS  = 2048;
FDA_FILE  = "data/list/fda_drugs.csv";
CHEM_FILE = "data/list/everyday_chemicals.csv";

fdaTbl  = readtable(FDA_FILE,  TextType="string");
chemTbl = readtable(CHEM_FILE, TextType="string");
logInfo("FDA 承認薬: %d 行  |  日用化学品: %d 行", ...
    height(fdaTbl), height(chemTbl));
%[text] FDA承認薬をパースしてフィンガープリントを計算します（クラス1）。
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
logInfo("FDA クラス:     %d 解析成功、%d 失敗", numel(fps1), failed1);
%[text] 日用化学品をパースしてフィンガープリントを計算します（クラス0）。
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
logInfo("非薬クラス:     %d 解析成功、%d 失敗", numel(fps0), failed0);
%[text] ビット行列とラベルベクトルを組み立てます。
n1     = numel(fps1);
n0     = numel(fps0);
nTotal = n1 + n0;

bitMat = zeros(nTotal, FP_NBITS, "single");
allFps = [fps1, fps0];
for j = 1:nTotal
    bitMat(j, :) = single(emk.fingerprint.toArray(allFps{j}));
end

labels   = [ones(n1, 1); zeros(n0, 1)];    % 1=薬, 0=非薬
allNames = [names1(:); names0(:)];

logInfo("合計: %d 分子（%d 薬、%d 非薬）", nTotal, n1, n0);
logInfo("クラスバランス: 薬 %.1f%% vs 非薬 %.1f%%", ...
    100*n1/nTotal, 100*n0/nTotal);
%[text] **💡 観察ポイント 1**
%[text] 薬と非薬の平均オンビット数はどれくらいかを確認しましょう。
%[text] 以下のコードを参考にしてください:
%[text] onBits1 = sum(bitMat(1:n1, :), 2);           % 薬ごとのカウント
%[text] onBits0 = sum(bitMat(n1+1:end, :), 2);        % 非薬ごとのカウント
%[text] \[mean(onBits1), mean(onBits0)\]
%[text] 医薬品分子の構造的特徴が多いか少ないかを読み取りましょう。
%[text] 無情報分類器が常に「薬」と予測する場合の精度を確認しましょう。
%[text] その精度はどうなるでしょうか?（答え: n1/nTotal）
%[text] なぜ精度はクラス不均衡データに対して誤解を招く指標なのかを考えてみましょう。
%[text] 医薬品でもある日用化学品を特定してみましょう。
%[text] pharma\_in\_chem = intersect(lower(names0), \["aspirin","caffeine","ibuprofen",...
%[text] "paracetamol","ethanol","nicotine"\]);
%[text] これらが意図的なラベルノイズの原因となります。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: PCA による次元削減
%[text] フィンガープリントの計算が完了しました。2048 次元のままでは過学習のリスクが高いため、PCA を用いて次元を削減します。
%[text] クラスの分布も併せて可視化します。
%[text] 
%[text] ### コンセプト: なぜ分類の前に PCA を行うのか?
%[text] 2048 のバイナリ特徴量と約 230 の分子では、p \>\> n の状況にあります。線形分類器はランダムラベルでも 100% の学習精度を達成できるため、過学習の可能性があります。また、2048 ビットの多くはほぼ一定で情報を持ちません。
%[text] 
%[text] PCA は 2048 ビットベクトルを最大分散方向の N\_PCS 次元に射影します（記述子に対する A01 と同様のアプローチです）。これにより、ほぼ一定のビットを自動的に除去し、2 クラスを最もよく分離する構造的変動を保持し、特徴行列を N × 2048 から N × N\_PCS に削減します。
%[text] 
%[text] **N\_PCS の選択:** フィンガープリントベースの機械学習では、10~50 個の主成分（PC）が文献で推奨されています。ここでは N\_PCS=20 を使用します。
%[text] 
%[text] **重要:** PCA は学習データのみでフィットし、テストや新規分子にも同じ変換を適用する必要があります。データを分割する前に全データで PCA をフィットするのは「データリーク」です。
logSection("A04", "セクション 2: PCA による次元削減", "アナリティクス L3");
N_PCS = 20;
%[text] 特徴量を中心化します（スケーリングは不要です — バイナリ特徴量はすでに \[0,1\] の範囲です）。
bitMean   = mean(bitMat, 1);
bitMat_c  = double(bitMat) - bitMean;
%[text] 期待される警告を抑制します: バイナリフィンガープリント列は線形従属です。
%[text] （ほぼ定数のビットはゼロに近い分散を持ちます。pca() 関数は T^2 成分数を減らして対処します — 結果は正しいです。）
pca_ws = warning('off', 'stats:pca:ColRankDefX');
[pcCoeff, pcScore, ~, ~, explained] = pca(bitMat_c, NumComponents=N_PCS);
warning(pca_ws);

logInfo("PCA: 上位 %d PC がフィンガープリント分散の %.1f%% を説明", ...
    N_PCS, sum(explained(1:N_PCS)));
%[text] クラスで色分けした最初の 2 つの主成分（PC）を可視化します。
figure("Name", "A04 フィンガープリントの PCA");
hold on;
scatter(pcScore(labels==1, 1), pcScore(labels==1, 2), 40, [0.3 0.6 0.9], ...
    "filled", MarkerFaceAlpha=0.6, DisplayName="FDA 承認薬（クラス 1）");
scatter(pcScore(labels==0, 1), pcScore(labels==0, 2), 70, [0.9 0.4 0.3], ...
    "^", "filled", MarkerFaceAlpha=0.9, DisplayName="日用化学品（クラス 0）");
hold off;
xlabel(sprintf("PC1 (%.1f%%)", explained(1)));
ylabel(sprintf("PC2 (%.1f%%)", explained(2)));
title("ECFP4 フィンガープリントの PCA -- 薬 vs 非薬");
legend(Location="best");
grid on;
%[text] **💡 観察ポイント 2**
%[text] PC1 と PC2 のプロットで、2 つのクラスが視覚的に分離しているか確認しましょう。
%[text] 重なりがある場合、それは 2 つのクラスが多くの構造的特徴を共有していることを示します。
%[text] PC 空間でどちらのクラスの分散が大きい（広がりが大きい）かを確認しましょう。
%[text] N\_PCS = 5、10、50 を試してみて、説明分散がどのように変わるかを確認しましょう。
%[text] スクリー曲線をプロットします: figure; plot(1:numel(explained), explained, "-o");
%[text] エルボー（収益逓減点）はどこにあるかを確認しましょう。
%[text] アスピリンとカフェインは日用化学品であると同時に医薬品でもあります。
%[text] PC 空間でのそれらの位置を見つけましょう。
%[text] aspIdx = find(allNames == "aspirin");
%[text] scatter(pcScore(aspIdx,1), pcScore(aspIdx,2), 100, "rx", LineWidth=2);
%[text] 薬クラスターと非薬クラスターのどちらに近いかを確認しましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: 層化学習 / テスト分割
%[text] 
%[text] PCA（主成分分析）で特徴空間が整いました。次にデータを学習セットとテストセットに分割します。
%[text] クラス不均衡があるため、層化分割を用いて各クラスの比率を維持します。
%[text] 
%[text] ### コンセプト: クラス不均衡に対応する層化分割
%[text] 単純なランダム 80/20 分割では、少数クラスのテストサンプルが偏る危険があります。層化分割は各クラスを独立して分割し、全データと同じ学習：テスト比率を各クラス内で維持します。
%[text] 
%[text] **正しい PCA の手順:** 本番環境では PCA は学習データのみでフィットし、同じ変換をテストデータに適用する必要があります。分割前に全データで PCA をフィットするのは「データリーク」です。
%[text] 
logSection("A04", "セクション 3: 層化学習 / テスト分割", "アナリティクス L3");
TEST_RATIO = 0.20;
rng(42);
cvPart   = cvpartition(labels, Holdout=TEST_RATIO);
trainIdx = training(cvPart);
testIdx  = test(cvPart);
%[text] 正しい方法: 学習データのみで PCA をフィットし、両セットに変換を適用します。
trainMean = mean(double(bitMat(trainIdx, :)), 1);
Xtrain_c  = double(bitMat(trainIdx, :)) - trainMean;
pca_ws2 = warning('off', 'stats:pca:ColRankDefX');
[pcCoeffTr, XTrain, ~, ~, explainedTr] = pca(Xtrain_c, NumComponents=N_PCS);
warning(pca_ws2);
%[text] テストデータに同じ変換を適用します。
Xtest_c = double(bitMat(testIdx, :)) - trainMean;
XTest   = Xtest_c * pcCoeffTr;   % project onto training PCs

yTrain = labels(trainIdx);
yTest  = labels(testIdx);

logInfo("学習: %d サンプル（%d 薬、%d 非薬）", ...
    sum(trainIdx), sum(yTrain==1), sum(yTrain==0));
logInfo("テスト: %d サンプル（%d 薬、%d 非薬）", ...
    sum(testIdx), sum(yTest==1), sum(yTest==0));
logInfo("学習 PCA: 上位 %d PC が学習分散の %.1f%% を説明", ...
    N_PCS, sum(explainedTr(1:N_PCS)));
%[text] **💡 観察ポイント 3**
%[text] 層化を確認しましょう: 学習セットとテストセットのクラス比率は同程度かを確認します。
%[text] 薬の割合（学習）: sum(yTrain==1)/numel(yTrain)
%[text] 薬の割合（テスト）: sum(yTest==1)/numel(yTest)
%[text] これらはほぼ等しい（2～3% 以内）はずです。
%[text] 層化なしでは何が起きるかを確認しましょう。
%[text] cvPart2 = cvpartition(nTotal, Holdout=0.2);  % 非層化（整数 n）
%[text] sum(labels(test(cvPart2)) == 0)   % テストセットの非薬 -- 0 になる場合があります!
%[text] rng シードを変えて（例: rng(1)、rng(7)）再実行し、テスト精度の変動を確認しましょう。
%[text] テスト精度はどれだけ変動するかを観察します。この変動は単一学習 / テスト分割に
%[text] 内在する「分割ノイズ」です。CV（観察ポイント 5）でこれを削減できます。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: SVM 分類器（fitcsvm）
%[text] 
%[text] 学習・テスト分割が完了しました。次に、最大マージン境界を学習する SVM を用いて分類を行います。
%[text] 混同行列を用いて予測の詳細を確認しましょう。
%[text] 
%[text] ### コンセプト: 二値分類のためのサポートベクターマシン（SVM）
%[text] SVM は 2 クラスを分離する最大マージン超平面を見つけます。マージン $= 2 / ||w||$ を最大化することは、$||w||^2$ を最小化することに等しいです。マージン境界上の点が「サポートベクター」となり、決定境界を定義します。
%[text] 
%[text] **主要パラメータ:** `BoxConstraint` (C) — マージン幅と学習誤差のトレードオフ（小さいCは広いマージン、大きいCは狭いマージン）/ `KernelFunction` — 特徴空間には `"linear"` / `Standardize` — 特徴量を平均0、単位分散に正規化
%[text] 
%[text] **クラス不均衡の処理:** 87% が薬で 13% が非薬の場合、デフォルトの SVM は誤分類の総数を最小化しようとするため、少数クラス（非薬）を無視する可能性があります。これは💡 観察ポイント 4 で対処します。
%[text] 
logSection("A04", "セクション 4: SVM 分類器（fitcsvm）", "アナリティクス L3");
svmModel = fitcsvm(XTrain, yTrain, ...
    KernelFunction="linear", ...
    BoxConstraint=1.0, ...
    Standardize=true, ...
    ClassNames=[0; 1]);

[yPred_svm, scores_svm] = predict(svmModel, XTest);
%[text] 分類評価指標（手動計算 — 追加ツールボックス関数不要）
cm_svm   = confusionmat(yTest, yPred_svm, Order=[0, 1]);
acc_svm  = sum(diag(cm_svm)) / sum(cm_svm(:));
prec_svm = cm_svm(2,2) / max(cm_svm(1,2) + cm_svm(2,2), 1);  % avoid /0
rec_svm  = cm_svm(2,2) / max(cm_svm(2,1) + cm_svm(2,2), 1);
f1_svm   = 2 * prec_svm * rec_svm / max(prec_svm + rec_svm, eps);

logInfo("SVM（線形、C=1）-- テストセット:");
logInfo("  精度=%.3f  適合率=%.3f  再現率=%.3f  F1=%.3f", ...
    acc_svm, prec_svm, rec_svm, f1_svm);
logInfo("  混同行列: TN=%d FP=%d FN=%d TP=%d", ...
    cm_svm(1,1), cm_svm(1,2), cm_svm(2,1), cm_svm(2,2));

figure("Name", "A04 SVM 混同行列");
confusionchart(categorical(yTest, [0 1], ["非薬","薬"]), ...
    categorical(yPred_svm, [0 1], ["非薬","薬"]), ...
    ColumnSummary="column-normalized", RowSummary="row-normalized", ...
    Title="SVM 混同行列");
logInfo("SVM 混同行列を表示しました");
%[text] **💡 観察ポイント 4**
%[text] 4 つのセルが何を意味するかを確認しましょう。
%[text] - TP = 薬として正しくラベル付けされた薬（真陽性）
%[text] - TN = 非薬として正しくラベル付けされた非薬（真陰性）
%[text] - FP = 薬として誤ってラベル付けされた非薬（偽陽性 — 「誤警報」）
%[text] - FN = 非薬として誤ってラベル付けされた薬（偽陰性 — 「見逃した薬」） \
%[text] 薬スクリーニングの文脈で、どのエラータイプがより高コストかを考えてみましょう。
%[text] 無情報率について考えましょう。「怠惰な」分類器は常に「薬」と予測します。
%[text] その精度は n1/nTotal です。SVM の精度が無情報率を上回っているか確認しましょう。
%[text] （そうでなければ、SVM はナイーブベースラインに価値を追加していないことになります！）
%[text] コスト感応型 SVM を使ってクラス不均衡を処理しましょう。
%[text] w = \[n1/nTotal, n0/nTotal\];   % 逆頻度重み
%[text] svmW = fitcsvm(XTrain, yTrain, KernelFunction="linear", ...
%[text] Cost=\[0 n1/n0; 1 0\]);   % cost(true=0, pred=1) = n1/n0
%[text] \[yW, ~\] = predict(svmW, XTest);
%[text] コスト重み付けによって少数クラスの再現率がどう変わるかを確認しましょう。
%[text] KernelFunction="rbf" を試してみましょう。非線形 SVM は精度を改善するか確認しましょう。
%[text] （注: 20 次元 PCA 空間の RBF カーネルは依然として非線形です）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: Random Forest 分類器（fitcensemble）
%[text] 
%[text] SVM が構築できました。次に、決定木のアンサンブルである Random Forest を試してみましょう。
%[text] クラス不均衡への対処法や特徴重要度など、SVM との違いを確認しましょう。
%[text] 
%[text] コンセプト: 分類のための Random Forest
%[text] A03 回帰と同様の木のアンサンブルアプローチですが、各木がクラスに投票します。
%[text] 1. B 本の木を構築し、各木はブートストラップサンプルとランダムな特徴量サブセットを使用します。
%[text] 2. 各木がクラス票を出力します（0 または 1）。
%[text] 3. 最終予測は多数決で決定されます（\> 50% の木が同意）。
%[text] 4. クラス確率は、そのクラスに投票した木の割合で表されます。 \
%[text] （ハード閾値よりも滑らかで、ROC に直接使用可能です）
%[text] 
%[text] SVM に対する優位性:
%[text] - 事後校正なしでクラス確率を得られます。
%[text] - 特徴量スケーリングへの感度が低く、Standardize が不要です。
%[text] - Prior="uniform" はクラス不均衡を扱うシンプルな方法です。
%[text] - 分割ゲインを用いた特徴重要度が解釈可能です。 \
%[text] 
%[text] Prior="uniform":
%[text] 実際の 87%/13% の分割に関係なく、事前クラス確率を 0.5/0.5 に設定します。
%[text] これにより、学習中に両クラスを等しく重要として扱い、非薬の再現率が向上します。
%[text] 
%[text] MinLeafSize=1: 完全成長木（最大の複雑度、わずかな過学習）。
%[text] より大きなデータセットでは、MinLeafSize を増やして過学習を削減します。
logSection("A04", "セクション 5: Random Forest 分類器（fitcensemble）", "アナリティクス L3");
rfClf = fitcensemble(XTrain, yTrain, ...
    Method="Bag", ...
    NumLearningCycles=200, ...
    Learners=templateTree(MinLeafSize=1), ...
    Prior="uniform", ...     % クラス不均衡を処理
    ClassNames=[0; 1]);

[yPred_rf, scores_rf] = predict(rfClf, XTest);

cm_rf   = confusionmat(yTest, yPred_rf, Order=[0, 1]);
acc_rf  = sum(diag(cm_rf)) / sum(cm_rf(:));
prec_rf = cm_rf(2,2) / max(cm_rf(1,2) + cm_rf(2,2), 1);
rec_rf  = cm_rf(2,2) / max(cm_rf(2,1) + cm_rf(2,2), 1);
f1_rf   = 2 * prec_rf * rec_rf / max(prec_rf + rec_rf, eps);

logInfo("Random Forest（200 木、一様事前分布）-- テストセット:");
logInfo("  精度=%.3f  適合率=%.3f  再現率=%.3f  F1=%.3f", ...
    acc_rf, prec_rf, rec_rf, f1_rf);
logInfo("  混同行列: TN=%d FP=%d FN=%d TP=%d", ...
    cm_rf(1,1), cm_rf(1,2), cm_rf(2,1), cm_rf(2,2));

figure("Name", "A04 RF 混同行列");
confusionchart(categorical(yTest, [0 1], ["非薬","薬"]), ...
    categorical(yPred_rf, [0 1], ["非薬","薬"]), ...
    ColumnSummary="column-normalized", RowSummary="row-normalized", ...
    Title="Random Forest 混同行列");
logInfo("RF 混同行列を表示しました");
%[text] 特徴重要度（PCA 成分ごと）
featureImp = predictorImportance(rfClf);
figure("Name", "A04 RF 特徴重要度（PCA 成分）");
bar(featureImp, FaceColor=[0.4 0.7 0.4]);
xlabel("PCA 成分"); ylabel("MDI 重要度");
title("PC1..PC20 の RF 特徴重要度");
grid on;
%[text] **💡 観察ポイント 5**
%[text] RF と SVM の混同行列を比較し、
%[text] 非薬クラスの再現率（FN が少ない）が高いのはどちらかを確認しましょう。
%[text] 薬クラスの適合率（FP が少ない）が高いのはどちらかを確認しましょう。
%[text] rfClf から Prior="uniform" を削除して再フィットし、
%[text] 非薬の再現率が低下するか、全体精度が向上するかを確認しましょう。
%[text] これはクラス不均衡下での精度と再現率のトレードオフを示しています。
%[text] 最も高い RF 特徴重要度を持つ PC 成分を確認し、
%[text] 最大分散を説明する PC1 か、それとも後の PC かを確認しましょう。
%[text] （高い分散は必ずしも高い識別力を意味しません）
%[text] MinLeafSize を 1 から 10 に変更し、テスト精度がどう変わるかを確認しましょう。
%[text] 学習精度も確認しましょう。
%[text] \[yTrain\_pred, ~\] = predict(rfClf, XTrain);
%[text] mean(yTrain\_pred == yTrain)   % 学習精度（MinLeafSize=1 では ~1 のはず）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: ROC 曲線の比較
%[text] 
%[text] 2つの分類器が揃いました。ここでは、閾値に依存しない総合的な指標であるAUC（曲線下面積）を用いて、両モデルを公正に比較します。
%[text] クラス不均衡があるデータでは、精度よりもAUCが信頼できる評価指標です。
%[text] 
%[text] コンセプト: 受信者操作特性（ROC）曲線
%[text] 二値分類器はスコア（確率または決定値）を出力します。
%[text] デフォルト閾値（0.5）はスコアを二値予測に変換します。
%[text] 閾値を変えることでトレードオフが生じます:
%[text] - 高い閾値 -\> 予測陽性が少ない -\> 再現率低下、適合率向上
%[text] - 低い閾値 -\> 予測陽性が多い -\> 再現率向上、適合率低下 \
%[text] 
%[text] ROC 曲線: 全閾値での真陽性率（再現率 / 感度）と
%[text] 偽陽性率（1 - 特異度）のプロットです。
%[text] 
%[text] AUC（ROC 曲線下面積）:
%[text] AUC = 1.0  完全分類器（全FPRでTPR=1）
%[text] AUC = 0.5  ランダム分類器（対角線 -- 識別力なし）
%[text] AUC = 0.0  完全逆（全予測を反転すれば完全）
%[text] 
%[text] ケモインフォマティクスにおけるAUCの解釈:
%[text] AUC \> 0.9  優秀
%[text] 0\.7-0.9    実用的に良いモデル
%[text] 0\.6-0.7    ギリギリ（かろうじて有用）
%[text] \< 0.6      不良（バグ / ラベルノイズを確認）
%[text] 
%[text] AUCは閾値非依存で、全体的な識別能力を測定します。
%[text] クラスが不均衡な場合には、精度よりも優先されます。
%[text] 
%[text] ラベルノイズに関する注意:
%[text] 一部の日用化学品が真の医薬品（クラス0として誤ラベル）である場合、
%[text] AUC = 1.0 は不可能です。達成可能な最大AUCはラベル品質によって制限されます。
%[text] 
%[text] **小サンプルでのROCの見た目について**
%[text] テストセットの非薬は6サンプルのみです。そのため、FPR軸は
%[text] 0, 1/6≈0.167, 2/6≈0.333, ... の7段階しか取れず、
%[text] ROCは教科書の滑らかなS字曲線ではなく、粗い**階段関数**になります。
%[text] これは実装や計算の誤りではなく、小テストセットの正常な挙動です。
%[text] 曲線は常に(0,0)から始まりますが、薬スコアが非薬より高い場合、
%[text] 左端にFPR=0のまま垂直に上がる区間（見えにくい）が現れます。
%[text] 滑らかなROCが必要な場合は、交差検証（💡 観察ポイント6）を使いましょう。
%[text] 
%[text] SVM: スコアの列2 = クラス1の決定値（高いほど薬）
scoreSVM = scores_svm(:, 2);
%[text] RF: 列2 = クラス1（薬）の事後確率
%[text] ClassNames=\[0;1\] なので列1 = P(クラス0)、列2 = P(クラス1)
scoreRF  = scores_rf(:, 2);

[xSVM, ySVM, ~, aucSVM] = perfcurve(yTest, scoreSVM, 1);
[xRF,  yRF,  ~, aucRF]  = perfcurve(yTest, scoreRF,  1);

figure("Name", "A04 ROC 曲線");
hold on;
plot(xSVM, ySVM, Color=[0.3 0.6 0.9], LineWidth=2, ...
    DisplayName=sprintf("SVM 線形 (AUC=%.3f)", aucSVM));
plot(xRF,  yRF,  Color=[0.4 0.7 0.4], LineWidth=2, ...
    DisplayName=sprintf("Random Forest (AUC=%.3f)", aucRF));
plot([0 1], [0 1], "--k", LineWidth=1, DisplayName="ランダムベースライン (AUC=0.500)");
hold off;
xlabel("偽陽性率（1 - 特異度）");
ylabel("真陽性率（感度 / 再現率）");
title("ROC 曲線 -- 薬 vs 非薬分類");
legend(Location="southeast");
grid on; axis square;   % 正方形アスペクト; axis equal と異なり xlim/ylim を [0,1] に維持
xlim([0 1]); ylim([0 1]);
%[text] サマリーテーブル
summaryTbl = table( ...
    ["SVM（線形）"; "Random Forest"], ...
    [acc_svm;  acc_rf],  ...
    [prec_svm; prec_rf], ...
    [rec_svm;  rec_rf],  ...
    [f1_svm;   f1_rf],   ...
    [aucSVM;   aucRF],   ...
    VariableNames=["モデル","精度","適合率","再現率","F1","AUC"]);
logInfo("性能サマリー:");
disp(summaryTbl);

logInfo("ROC AUC: SVM=%.3f  RF=%.3f", aucSVM, aucRF);
%[text] **💡 観察ポイント 6**
%[text] どちらのモデルのAUCが高いかを確認しましょう。両方でAUC \> 0.5かどうかも確認します。
%[text] AUC \<= 0.5 はランダムより良くないことを意味し、危険信号です。
%[text] 「最適閾値」を見つけましょう（ROCの左上角に最も近い点）。
%[text] \[~, optIdx\] = min((1-ySVM).^2 + xSVM.^2);   % (FPR=0, TPR=1) への距離
%[text] この閾値での感度と特異度を計算し、デフォルトの0.5と異なるか確認しましょう。
%[text] この閾値はデフォルトの0.5と異なるか?
%[text] ヌルモデルチェックを行いましょう（ラベルをシャッフルして再実行）。
%[text] shuffleIdx = randperm(numel(yTest));
%[text] \[~, ~, ~, aucNull\] = perfcurve(yTest, scoreSVM(shuffleIdx), 1);
%[text] シャッフルされたラベルのAUCは約0.5のはずです。パイプラインが
%[text] 正しく実装され、データリークがないことを確認します。
%[text] 5分割交差検証を実行して平均AUCを計算しましょう。
%[text] mdlCV = crossval(svmModel, KFold=5);
%[text] \[yCV, scoresCV\] = kfoldPredict(mdlCV);
%[text] \[~,~,~,aucCV\] = perfcurve(yTrain, scoresCV(:,2), 1);
%[text] CV AUCと単一テストセットAUCを比較しましょう。
%[text] （小規模データセットではCV AUCがより信頼性の高い推定値です）
%[text] 日用薬品（アスピリン、カフェイン）のラベルノイズが最大AUCを制限します。
%[text] これらを削除して再実行しましょう。AUCは改善するか、どれくらいか確認します。
%[text] **まとめ**
%[text] - p\>\>nﾈ2048 ビット 1 230 分子ではSVM前にPCAが必要です。
%[text] - クラス不均衡は精度を誤解を招くものにします — AUCとF1を使いましょう。
%[text] - RFのPrior=uniformとSVMのコスト重み付けがクラス不均衡に対処します。
%[text] - AUCは識別能力の閾値非依存な測定値です。
%[text] - ラベル品質が達成可能な最大AUCを制限します。
%[text] - ヌルモデル（シャッフルされたラベル）で常に検証することが重要です。 \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
