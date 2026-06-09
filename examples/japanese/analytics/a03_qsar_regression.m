%[text] # A03: QSAR 回帰 — 分子記述子から LogP を予測する
%[text] EasyMolKit アナリティクス — レイヤー 3
%[text] 
%[text] LogP（オクタノール-水分配係数の常用対数）は、薬物が脂質膜を通過して体内に吸収されるかどうかを決定する重要な物性値です。
%[text] 実験での測定には化合物の合成と液-液分配アッセイが必要で、コストも時間もかかります。
%[text] 分子の「形」を表す8種類の構造記述子だけからLogPを予測できれば、合成前に何千もの候補を仮想スクリーニングできます。
%[text] このスクリプトでは線形回帰（`fitlm`）とランダムフォレスト（`fitrensemble`）の2つのモデルを構築し、5分割交差検証で性能を公正に比較します。
%[text] 
%[text] **このチュートリアルで学べること**
%[text] - 早期創薬におけるQSAR回帰の重要性を理解する
%[text] - `fitlm()`を用いて線形回帰モデルを構築・解釈する
%[text] - `fitrensemble()`を用いて非線形アンサンブル（ランダムフォレスト）回帰を行う
%[text] - k分割交差検証を適用して汎化性能を推定する
%[text] - R²、RMSE、実測値 vs 予測値プロットでモデルを比較する
%[text] - 適用可能ドメインの概念を理解する \
%[text] 
%[text] **前提条件**
%[text] - F02（特性計算）の完了 — 記述子の基礎
%[text] - 推奨: A01（化学空間PCA）で記述子のコンテキストを把握
%[text] - Statistics and Machine Learning Toolbox（`fitlm`、`fitrensemble`、`cvpartition`）
%[text] - インターネット接続不要 \
%[text] 
%[text] 推定所要時間: 30〜45分
%[text] 
%[text] **データ:**
%[text] `data/list/fda_drugs.csv` — 200種のFDA承認薬（ChEMBL、CC-BY-SA 3.0）
%[text] ALogP列: ChEMBLに格納されたCrippen-Wildman推定値。この演習では回帰目的変数として使用します（「正解」として扱います）。
%[text] 
%[text] **参考文献**
%[text] - Wildman SA & Crippen GM (1999) Prediction of physicochemical parameters by atomic contributions. J Chem Inf Comput Sci 39:868-873.  doi:10.1021/ci990307l (Crippen-Wildman LogP / ALogP) \
%[text] - Breiman L (2001) Random forests. Machine Learning 45:5-32. doi:10.1023/A:1010933404324 \
%[text] - Cherkasov A et al. (2014) QSAR modeling: where have you been?  J Med Chem 57:4977-5010. doi:10.1021/jm4004285 \
%[text] - Oprea TI (2000) Property distribution of drug-related chemical databases.  J Comput Aided Mol Des 14:251-264. doi:10.1023/A:1008130001697 \
%[text] 
%[text] 実行方法: Ctrl+Enterでセクションを1つずつ実行
%%
%[text] ## セクション 0: セットアップ
logSection("A03", "セクション 0: セットアップ", "アナリティクス L3");
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
%[text] メイン実行前に Python と RDKit プロセスを事前に準備します
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;
%%
%[text] ## セクション 1: FDA 承認薬の読み込みと記述子計算
%[text] 
%[text] セットアップが完了しました。まず、予測に使用するデータと特徴量を準備します。
%[text] SMILES から記述子を計算し、8 種類の特徴量行列 X と目的変数 y を作成します。
%[text] 
%[text] ### コンセプト: 特徴量選択 — なぜ LogP を除外するのか？
%[text] 目的変数は ALogP（ChEMBL の Crippen-Wildman 推定値）です。
%[text] RDKit の組み込み LogP 記述子も同じ Crippen-Wildman 法を使用しているため、特徴量に含めるとLogP を読んで LogP を予測する」という自明なモデルになり、意味のある QSAR になりません。
%[text] 
%[text] 代わりに、親油性を直接表現せずに分子構造をエンコードする8 種類の相補的記述子を使用します:
%[text] - サイズ: MolWt（g/mol）、HeavyAtomCount（非水素原子数）
%[text] - 極性: TPSA（Å²）、NumHDonors、NumHAcceptors
%[text] - トポロジー: NumRotatableBonds、RingCount
%[text] - 形状: FractionCSP3（sp3 炭素の割合） \
%[text] これらの記述子による LogP 予測の化学的根拠:
%[text] - TPSA / HBD / HBA が高い -\> より極性 -\> LogP 低下（水への溶解度増加）
%[text] - RingCount が多い -\> 炭素芳香族環は LogP 上昇に寄与するが、極性ヘテロ環が多い FDA 承認薬では方向が反転することもある（モデル係数を確認）
%[text] - FractionCSP3 が高い -\> より 3D 形状 -\> 平面性低下 -\> LogP 低下
%[text] - MolWt が大きいと創薬空間内では LogP と弱い正相関がある \
DATA_FILE  = "data/list/fda_drugs.csv";
FEAT_NAMES = ["MolWt", "TPSA", "NumHDonors", "NumHAcceptors", ...
              "NumRotatableBonds", "RingCount", "FractionCSP3", "HeavyAtomCount"];
N_FEATS    = numel(FEAT_NAMES);

rawTbl = readtable(DATA_FILE, TextType="string");
logInfo("%s から %d 個の分子を読み込みました", DATA_FILE, height(rawTbl));
%[text] ALogP 列の取得（パーサーによって数値型または文字列型で読み込まれる場合があります）
if isnumeric(rawTbl.ALogP)
    alogpVec = rawTbl.ALogP;
else
    alogpVec = str2double(rawTbl.ALogP);
end
%[text] SMILES を解析して記述子を計算します
nRaw  = height(rawTbl);
X_all = nan(nRaw, N_FEATS);
y_all = nan(nRaw, 1);
valid = false(nRaw, 1);

for k = 1:nRaw
    if isnan(alogpVec(k)), continue; end      % ALogP が欠損の場合はスキップ
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
logInfo("%d / %d 分子の記述子計算完了", sum(valid), nRaw);

X        = X_all(valid, :);    % N×8 特徴量行列（生データ、非標準化）
y        = y_all(valid);       % N×1 目的変数ベクトル（ALogP）
molNames = rawTbl.Name(valid);
nMols    = sum(valid);

logInfo("データセット: %d / %d 分子を解析しました", nMols, nRaw);
logInfo("ALogP の範囲: %.2f〜%.2f（平均=%.2f、標準偏差=%.2f）", ...
    min(y), max(y), mean(y), std(y));
%[text] 先頭 5 行をテーブル形式で確認します
prevN = min(5, nMols);
prevTbl = array2table(X(1:prevN, :), VariableNames=cellstr(FEAT_NAMES));
prevTbl.ALogP = y(1:prevN);
prevTbl.Name  = molNames(1:prevN);
prevTbl = movevars(prevTbl, "Name", Before=1);
disp(prevTbl);
%[text] **💡 観察ポイント 1**
%[text] ALogP と最も強いピアソン相関を持つ記述子を確認しましょう。
%[text] corr(X, y) を使用して 8x1 の相関ベクトルを確認します。
%[text] 各相関の符号が化学的に妥当かを読み取りましょう。
%[text] （ヒント: TPSA は LogP と負の相関があります -- 極性が高い = 脂溶性が低い）
%[text] 8 つの特徴量間のペアワイズ相関を計算し、確認しましょう: corr(X)
%[text] 最も多重共線性が高いペアを確認しましょう（予想: MolWt と HeavyAtomCount）。
%[text] 強い多重共線性が線形回帰係数にどう影響するかを考察しましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: LogP 分布の探索
%[text] 
%[text] 記述子の計算が完了しました。次に、モデリングの前に目的変数である ALogP の分布を可視化しましょう。
%[text] データを事前に確認することで、後のモデル結果を正確に解釈する助けになります。
%[text] 
%[text] ### コンセプト: モデリング前の目的変数の検査
%[text] 線形回帰は、残差がほぼ正規分布していることを前提としています。
%[text] （y 自体が正規分布である必要はありませんが、関連性があります）。
%[text] y を先に調べることで分かること:
%[text] 1. 形状: 偏った分布は対数変換や Box-Cox 変換で改善できる場合があります。
%[text] 2. 外れ値: 極端な y 値は OLS フィットを大きく歪める可能性があります。
%[text] 3. 範囲: モデルの適用可能ドメインは学習範囲によって制限されます。 \
%[text] この範囲を超えた外挿は信頼性が低くなります。
%[text] 
%[text] 薬物様 LogP（Lipinski の 5 則: LogP ≤ 5）を考慮すると、FDA 承認薬の分布は -2 から +5 の範囲でピークを持ち、親油性の外れ値（溶解性が低い薬も含まれる）に向けた裾野があります。
figure("Name", "A03 LogP 分布");
tiledlayout(1, 2);

nexttile;
histogram(y, 20, FaceColor=[0.3 0.6 0.9]);
xlabel("ALogP"); ylabel("頻度");
title("ALogP の分布（FDA 承認薬）");
xline(mean(y),   "--r", sprintf("mean=%.1f",   mean(y)),   LabelHorizontalAlignment="right");
xline(median(y), "--g", sprintf("median=%.1f", median(y)), LabelHorizontalAlignment="left");
grid on;

nexttile;
scatter(y, X(:, strcmp(cellstr(FEAT_NAMES), "TPSA")), 30, [0.3 0.6 0.9], ...
    "filled", MarkerFaceAlpha=0.6);
xlabel("ALogP"); ylabel("TPSA (A^2)");
title(sprintf("TPSA vs ALogP  r=%.2f", corr(y, X(:, strcmp(cellstr(FEAT_NAMES), "TPSA")))));
grid on;

logInfo("ALogP: 平均=%.2f、中央値=%.2f、標準偏差=%.2f、歪度=%.2f", ...
    mean(y), median(y), std(y), skewness(y));
%[text] 特徴量と目的変数の散布図（2×4 グリッド）
figure("Name", "A03 特徴量-LogP 散布図");
for f = 1:N_FEATS
    subplot(2, 4, f);
    scatter(X(:, f), y, 20, [0.2 0.5 0.8], "filled", MarkerFaceAlpha=0.5);
    xlabel(FEAT_NAMES(f)); ylabel("ALogP");
    r = corr(X(:, f), y);
    title(sprintf("r=%.2f", r));
    grid on;
end
sgtitle("8 種の記述子 vs ALogP");
%[text] **💡 観察ポイント 2**
%[text] ALogP \> 6 の分子が存在するかを確認しましょう。
%[text] \[sorted\_y, ord\] = sort(y, "descend");
%[text] table(molNames(ord(1:5)), sorted\_y(1:5), VariableNames=\["Name","ALogP"\])
%[text] これらの外れ値が線形モデルに悪影響を及ぼす可能性があるかを考えてみましょう。
%[text] FDA 承認薬全体の LogP 分布は 2.5 付近でピークを持つとされています（Lipinski 空間）。
%[text] このデータセットがその期待に合致しているかを確認しましょう。
%[text] 明確な非線形関係を示す記述子があるかを見てみましょう。
%[text] （非線形パターンは線形回帰の限界を示すことがあります）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: 線形回帰（fitlm）
%[text] 
%[text] 分布の特徴を確認しました。ここでは、最も基本的なモデルである線形回帰を用いて分析を始めます。
%[text] 線形モデルは係数の解釈が容易で、ベースラインとして重要です。
%[text] 
%[text] ### コンセプト: QSAR のための重回帰分析
%[text] 線形 QSAR モデルは LogP（分配係数の常用対数）を記述子の重み付き和として予測します:
%[text] 
%[text] LogP\_pred = b0 + b1*MolWt + b2*TPSA + ... + b8\*HeavyAtomCount
%[text] 
%[text] 係数 b1...b8 は最小二乗法（OLS）で推定されます:
%[text] 最小化: sum\_i (y\_i - X\_i \* b)^2
%[text] 
%[text] フィット前に特徴量を標準化（z スコア）することで比較が可能になります:
%[text] - 標準化データで |値| が大きい係数は、その記述子がこのデータセットで LogP に強い影響を持つことを示します。 \
%[text] - 標準化なしの係数は異なるスケール（例: b\_MolWt の単位は LogP / g mol^-1）です。 \
%[text] 
%[text] fitlm の主要な診断指標:
%[text] - R2      -- 説明された分散の割合（学習セット; 楽観的）
%[text] - RMSE    -- 二乗平均平方根誤差（LogP 単位）
%[text] - p 値   -- この係数がゼロと有意に異なるかどうか
%[text] - 残差   -- フィット値に対してパターンなく正規分布すべきです \
%[text] 
%[text] 特徴量を標準化する（平均 0、分散 1 に正規化）
Xmean = mean(X, 1);
Xstd  = std(X,  0, 1);
Xstd(Xstd < 1e-12) = 1;      % 定数記述子でのゼロ除算を防ぐ
Xs    = (X - Xmean) ./ Xstd;  % N×8 標準化済み特徴量行列
%[text] 全学習データで線形モデルをフィットします
ftTbl = array2table(Xs, VariableNames=cellstr(FEAT_NAMES));
ftTbl.ALogP = y;
lmModel = fitlm(ftTbl, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));

logInfo("線形モデル（学習）: R2=%.3f  RMSE=%.3f LogP 単位", ...
    lmModel.Rsquared.Ordinary, lmModel.RMSE);
disp(lmModel);
%[text] 実測値と予測値のプロットおよび残差を確認します
yPred_lm  = predict(lmModel, ftTbl);
residuals = y - yPred_lm;

figure("Name", "A03 線形モデル診断");
set(gcf, "Position", [100 100 1100 380]);
tiledlayout(1, 3);

nexttile;
scatter(y, yPred_lm, 30, [0.3 0.6 0.9], "filled", MarkerFaceAlpha=0.7);
hold on;
refLim = [min(y)-0.5, max(y)+0.5];
plot(refLim, refLim, "--k", LineWidth=1.5);
hold off;
xlabel("実測 ALogP"); ylabel("予測 ALogP");
title(sprintf("線形モデル\nR^2=%.3f  RMSE=%.3f", lmModel.Rsquared.Ordinary, lmModel.RMSE));
grid on;

nexttile;
histogram(residuals, 20, FaceColor=[0.9 0.5 0.3]);
xlabel("残差（実測 - 予測）"); ylabel("頻度");
title("残差の分布");
xline(0, "--k"); grid on;

nexttile;
[sortCoeff, sortOrd] = sort(lmModel.Coefficients.Estimate(2:end), "descend");
barh(fliplr(sortCoeff), FaceColor=[0.3 0.6 0.9]);
yticks(1:N_FEATS);
yticklabels(fliplr(cellstr(FEAT_NAMES(sortOrd))));
xlabel("標準化係数"); xline(0, "k");
title("回帰係数（標準化後）"); grid on;

logInfo("最も正の影響が強い記述子: %s（係数=%.3f）", ...
    FEAT_NAMES(sortOrd(1)), sortCoeff(1));
logInfo("最も負の影響が強い記述子: %s（係数=%.3f）", ...
    FEAT_NAMES(sortOrd(end)), sortCoeff(end));
%[text] **💡 観察ポイント 3**
%[text] 標準化係数が最も大きい正の記述子はどれかを確認しましょう。
%[text] 最も大きい負の係数はどれかを確認しましょう。
%[text] これらの符号がセクション 1 の化学的な期待と一致しているか読み取りましょう。
%[text] p 値 \> 0.05（統計的に有意でない）の特徴量があるか確認しましょう。
%[text] それらを除去して再フィットすると RMSE が改善するかを読み取りましょう。
%[text] 例えば HeavyAtomCount（MolWt と多重共線性あり）を削除します:
%[text] ftTbl2 = ftTbl; ftTbl2.HeavyAtomCount = \[\]; 
%[text] lm2 = fitlm(ftTbl2, "ALogP ~ MolWt + TPSA + NumHDonors + ...");
%[text] 残差の分布がほぼ正規に見えるか確認しましょう。
%[text] 右に歪んだ残差分布は外れ値の可能性を示します。
%[text] 最も大きい |残差| を持つ分子はどれかを確認しましょう。
%[text] \[~, bigResOrd\] = sort(abs(residuals), "descend");  molNames(bigResOrd(1:5))
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: Random Forest 回帰（fitrensemble）
%[text] 
%[text] 線形モデルを構築しました。次に、非線形な関係も捉えられる Random Forest を試してみましょう。logSection("A03", "セクション 4: Random Forest 回帰（fitrensemble）", "アナリティクス L3");%\[text\] 後のセクションで、2つのモデルを交差検証により公正に比較します。
%[text] 
%[text] ### コンセプト: QSAR のための Random Forest
%[text] Random Forest（Breiman 2001）は、B 本の決定木のアンサンブルです:
%[text] 1. 学習データから B 個のブートストラップサンプルを抽出します（重複あり）。
%[text] 2. 各サンプルに決定木をフィットし、各分岐でランダムに m 個の特徴量のみを考慮します。 \
%[text] （mtry = floor(p/3) がデフォルトで、木間の相関を低減します）
%[text] 1. 予測は、B 本の木の出力を平均して行います。 \
%[text] 
%[text] 線形回帰に対する優位性:
%[text] - 非線形関係を捉えることができる（木の分岐は非線形です）。
%[text] - 特徴量間の交互作用を扱える（分岐は2つの特徴量に依存できます）。
%[text] - 分岐純度低下（MDI）による組み込みの特徴重要度が得られます。
%[text] - 外れ値に対して頑健です（木はランク的な分岐を使用します）。 \
%[text] 
%[text] 欠点: 学習 R2 は非常に楽観的です（木はデータを記憶できるため）。常に交差検証（セクション 5）で評価することが重要です。
%[text] 
%[text] MDI（平均不純度減少）による特徴重要度:
%[text] 特徴量 j を使用する全ての木・分岐において、分散減少 × その節点のサンプル割合の和を計算します。値が大きいほど有用です。
rfModel = fitrensemble(Xs, y, ...
    Method="Bag", ...
    NumLearningCycles=200, ...
    Learners=templateTree(MinLeafSize=3));

yPred_rfTrain  = predict(rfModel, Xs);
ss_res_rf      = sum((y - yPred_rfTrain).^2);
ss_tot         = sum((y - mean(y)).^2);
r2_rf_train    = 1 - ss_res_rf / ss_tot;
rmse_rf_train  = sqrt(mean((y - yPred_rfTrain).^2));
logInfo("Random Forest（学習、楽観的）: R2=%.3f  RMSE=%.3f", ...
    r2_rf_train, rmse_rf_train);
logInfo("注意: RF の学習指標は楽観的 -- セクション 5 の CV で確認してください");
%[text] 特徴重要度の可視化
importances = predictorImportance(rfModel);
[sortedImp, sortOrdRF] = sort(importances, "descend");

refLim = [min(y)-0.5, max(y)+0.5];   % セクション 3 から独立して実行できるよう再定義
figure("Name", "A03 Random Forest");
tiledlayout(1, 2);

nexttile;
barh(fliplr(sortedImp), FaceColor=[0.4 0.7 0.4]);
yticks(1:N_FEATS);
yticklabels(fliplr(cellstr(FEAT_NAMES(sortOrdRF))));
xlabel("特徴重要度（MDI）");
title("RF 特徴重要度");
grid on;

nexttile;
scatter(y, yPred_rfTrain, 30, [0.4 0.7 0.4], "filled", MarkerFaceAlpha=0.7);
hold on;
plot(refLim, refLim, "--k", LineWidth=1.5);
hold off;
xlabel("実測 ALogP"); ylabel("予測（RF、学習）");
title(sprintf("RF 学習  R^2=%.3f  RMSE=%.3f", r2_rf_train, rmse_rf_train));
grid on;

logInfo("RF 最重要特徴量: %s  最低重要度: %s", ...
    FEAT_NAMES(sortOrdRF(1)), FEAT_NAMES(sortOrdRF(end)));
%[text] **💡 観察ポイント 4**
%[text] Random Forest が最も重要と評価する特徴量を確認しましょう。
%[text] 線形モデルの最大係数と一致しているかを読み取りましょう。
%[text] （RF は非線形な寄与を捉えるため、異なる場合があります。）
%[text] MinLeafSize=1 と MinLeafSize=20 を試した場合の変化を確認しましょう。
%[text] 学習 RMSE はどう変わるか、MinLeafSize=1 のリスクを考察しましょう。
%[text] （葉が小さい = 深い木 = 学習データへの過学習が増える可能性があります。）
%[text] NumLearningCycles を 500 に増やした場合、学習 RMSE が改善するか確認しましょう。
%[text] 何本の木でモデルが「安定化」するか（収益逓減点）を考察しましょう。
%[text] oobLoss(rfModel) は木の数の関数として OOB RMSE を返します。
%[text] 試してみましょう: plot(oobLoss(rfModel)); xlabel("木の数"); ylabel("OOB RMSE");
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: 5 分割交差検証による比較
%[text] 
%[text] 2 つのモデルが揃いました。ここが今回の重要なポイントです。
%[text] 学習データで評価した誤差は楽観的になりがちなので、交差検証を用いて汎化性能を公正に推定します。
%[text] 
%[text] ### コンセプト: 不偏なモデル比較のための交差検証
%[text] 学習誤差は常に楽観的です。これはモデルがすでにデータを見ているためです。
%[text] k 分割交差検証（CV）は汎化誤差を推定する手法です。
%[text] 1. n 個の分子を k 個の等しいフォールドに分割します。
%[text] 2. 各フォールド i ごとに以下を行います。 \
%[text] - フォールド i を除いたデータで学習します。
%[text] - フォールド i（学習中に未使用）で予測します。
%[text] - 残差を記録します。
%[text] - フォールド外予測を集計して CV 評価指標を計算します。 \
%[text] 
%[text] k=5 が標準的で、1 フォールドあたり 80% のデータで学習します。小規模データセット（n\<100）では k=10 または LOO（1 件抜き）を用います。
%[text] 
%[text] 評価指標は以下の通りです。
%[text] -  - R2 = 1 - SS\_res/SS\_tot（1=完全、0=平均予測と等価、負=平均以下）
%[text] -  - RMSE = sqrt(mean((y\_true - y\_pred)^2))  \[LogP 単位\]
%[text] -  - MAE = mean(|y\_true - y\_pred|)（RMSE より外れ値の影響を受けにくい） \
%[text] 
%[text] 優れた QSAR 回帰モデルは通常 CV R2 \> 0.6 を達成します。CV R2 \> 0.8 は薬物様物性予測において優秀とみなされます。
K_FOLDS = 5;
rng(42);
cv = cvpartition(nMols, KFold=K_FOLDS);

yCV_lm = nan(nMols, 1);   % out-of-fold predictions: linear model
yCV_rf = nan(nMols, 1);   % out-of-fold predictions: random forest

for fold = 1:K_FOLDS
    trainIdx = training(cv, fold);
    testIdx  = test(cv, fold);

    % このフォールドの線形モデル
    trTbl = array2table(Xs(trainIdx, :), VariableNames=cellstr(FEAT_NAMES));
    trTbl.ALogP = y(trainIdx);
    lm_f  = fitlm(trTbl, "ALogP ~ " + strjoin(FEAT_NAMES, " + "));
    teTbl = array2table(Xs(testIdx, :), VariableNames=cellstr(FEAT_NAMES));
    yCV_lm(testIdx) = predict(lm_f, teTbl);

    % このフォールドの Random Forest
    rf_f = fitrensemble(Xs(trainIdx, :), y(trainIdx), ...
        Method="Bag", NumLearningCycles=200, ...
        Learners=templateTree(MinLeafSize=3));
    yCV_rf(testIdx) = predict(rf_f, Xs(testIdx, :));
end
logInfo("%d 分割 CV 完了", K_FOLDS);
%[text] CV 評価指標を計算します（無名関数で構造体を生成）。
computeMetrics = @(yt, yp) struct( ...
    "R2",   1 - sum((yt-yp).^2) / sum((yt-mean(yt)).^2), ...
    "RMSE", sqrt(mean((yt-yp).^2)), ...
    "MAE",  mean(abs(yt-yp)));

cv_lm = computeMetrics(y, yCV_lm);
cv_rf = computeMetrics(y, yCV_rf);

logInfo("--- %d 分割 CV 結果 ---", K_FOLDS);
logInfo("線形モデル:     R2=%.3f  RMSE=%.3f  MAE=%.3f", ...
    cv_lm.R2, cv_lm.RMSE, cv_lm.MAE);
logInfo("Random Forest: R2=%.3f  RMSE=%.3f  MAE=%.3f", ...
    cv_rf.R2, cv_rf.RMSE, cv_rf.MAE);
logInfo("過学習ギャップ（RF）: 学習 RMSE=%.3f  CV RMSE=%.3f", ...
    rmse_rf_train, cv_rf.RMSE);
%[text] ### なぜ今回は線形回帰が優れたのか?
%[text] このデータセットの目的変数 ALogP は Crippen-Wildman 法（原子寄与の**線形和**）で計算された推定値です。生成メカニズムが線形であるため、線形回帰はデータ構造と相性が良く、n=200 の小規模データセットでは RF の柔軟性がかえって過学習に繋がりました。
%[text] 実験で測定した実測 LogP（非線形な分子間相互作用を含む）を目的変数にするとRF が優位になる場合があります。「強力なモデルが常に勝つわけではない」という重要な教訓です。モデル選択には必ず CV を使用しましょう。
%[text] CV 実測値と予測値の比較プロット
figure("Name", "A03 交差検証比較");
tiledlayout(1, 2);

nexttile;
scatter(y, yCV_lm, 30, [0.3 0.6 0.9], "filled", MarkerFaceAlpha=0.7);
hold on; plot(refLim, refLim, "--k", LineWidth=1.5); hold off;
xlabel("実測 ALogP"); ylabel("CV 予測 ALogP");
title(sprintf("線形モデル（5 分割 CV）\nR^2=%.3f  RMSE=%.3f", cv_lm.R2, cv_lm.RMSE));
grid on;

nexttile;
scatter(y, yCV_rf, 30, [0.4 0.7 0.4], "filled", MarkerFaceAlpha=0.7);
hold on; plot(refLim, refLim, "--k", LineWidth=1.5); hold off;
xlabel("実測 ALogP"); ylabel("CV 予測 ALogP");
title(sprintf("Random Forest（5 分割 CV）\nR^2=%.3f  RMSE=%.3f", cv_rf.R2, cv_rf.RMSE));
grid on;
%[text] **💡 観察ポイント 5**
%[text] CV R2 と RMSE からどちらのモデルが優れているかを確認しましょう。R2 の差が 0.1 を超えるか（実践的に有意か）を確認しましょう。
%[text] 学習 RMSE（セクション 4）と CV RMSE（ここ）を比較し、「過学習ギャップ」= CV RMSE - 学習 RMSE を確認しましょう。
%[text] どちらのモデルのギャップが大きいか、それが何を意味するかを考えましょう。
%[text] K\_FOLDS を 10 に変更した場合、CV 結果がどのように変わるかを確認しましょう。
%[text] n=200 の場合、5 分割 vs 10 分割 CV ではバイアス-分散のトレードオフが異なります。
%[text] （フォールドが多いほどバイアスが少なく、CV 推定値の分散が大きくなります。）
%[text] RF のハイパーパラメータチューニングを試し、
%[text] rfTuned = fitrensemble(Xs, y, "OptimizeHyperparameters", "auto", ...
%[text] HyperparameterOptimizationOptions=struct(MaxObjectiveEvaluations=20));
%[text] ベイズ最適化がデフォルト設定より CV R2 を改善するかを確認しましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: 新規分子へのモデル適用
%[text] 
%[text] 交差検証で優れたモデルを選びました。最後に、このモデルを実際に適用してみましょう。
%[text] 学習データに含まれない新規分子に対して LogP（分配係数の常用対数）を予測し、適用可能ドメインも確認します。
%[text] 
%[text] ### コンセプト: 予測パイプラインと適用可能ドメイン
%[text] FDA 承認薬で学習した QSAR モデルは、学習セットと類似した分子（適用可能ドメイン内）で最も信頼性が高くなります。
%[text] 
%[text] 適用可能ドメインチェック（単純な範囲法）:
%[text] - 学習セットの平均と標準偏差を使って新規分子の記述子を標準化します。
%[text] - 標準化された特徴量が学習範囲を 2 シグマ以上外れている場合、 \
%[text] その予測を「適用範囲外」（信頼性が低い）としてフラグを立てます。
%[text] より厳格なドメイン手法もあります（レバレッジベースのユークリッド距離や k-NN 距離）が、範囲チェックは実用的な最初のスクリーニングとして有効です。
%[text] 
%[text] 実用上の注意:
%[text] 新規分子への適用前に、選ばれたモデルを全データで再学習します。（データが多いほどモデルが良くなります。CV はモデル選択のためだけです）。
%[text] 
%[text] CV R2 に基づいてモデルを選択します。
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
logInfo("最終モデル: %s（全 %d 分子で再学習）", modelLabel, nMols);
%[text] 予測する新規分子
newSmiles = ["CCO", ...                              % エタノール
             "c1ccccc1", ...                         % ベンゼン
             "CC(=O)Oc1ccccc1C(=O)O", ...           % アスピリン
             "CN1C=NC2=C1C(=O)N(C)C(=O)N2C", ...   % カフェイン
             "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34"]; % テストステロン

logInfo("%d 種の新規分子の LogP を予測中:", numel(newSmiles));
for k = 1:numel(newSmiles)
    try
        mol = emk.mol.fromSmiles(newSmiles(k));
        s   = emk.descriptor.calculate(mol, FEAT_NAMES);
        xNew = zeros(1, N_FEATS);
        for d = 1:N_FEATS
            xNew(d) = s.(FEAT_NAMES(d));
        end
        xNewS = (xNew - Xmean) ./ Xstd;   % 学習統計量で標準化

        % 適用可能ドメイン: 2σ 範囲チェック（学習統計量で標準化後、全特徴量が |z| <= 2 以内か）
        inDomain = all(abs(xNewS) <= 2);
        domainStr = "OK";
        if ~inDomain, domainStr = "OUT-OF-DOMAIN"; end

        logpPred = predictFn(xNewS);
        logInfo("  %-46s -> LogP_pred=%5.2f  [%s]", newSmiles(k), logpPred, domainStr);
    catch ME
        logWarn("  %s で失敗: %s", newSmiles(k), ME.message);
    end
end
%[text] **💡 観察ポイント 6**
%[text] 予測された LogP を RDKit の組み込み Crippen 推定値と比較してみましょう。
%[text] mol = `emk.mol.fromSmiles("CCO")`;
%[text] s = `emk.descriptor.calculate(mol, "LogP")`;  s.LogP
%[text] QSAR 予測が Crippen 値にどれくらい近いかを確認しましょう。
%[text] ALogP（学習目的変数）も Crippen 法で計算されているため、類似するはずです。
%[text] つまり、これは Crippen の「メタモデル」といえます。
%[text] 非常に大きな分子でも試してみましょう（例: ペプチド Ala-Gly: "NCC(=O)NCC(=O)O"）。
%[text] 適用範囲チェックで「適用範囲外」がトリガーされるかを確認しましょう。
%[text] それでも予測を使った場合、何が起こるかを考えてみましょう。
%[text] なぜ本番展開前に最終モデルを全データで再学習すべきかを考えましょう。
%[text] （ヒント: CV はフォールドごとに 80% のデータで学習します。
%[text] 全データを使うとより良いモデルが得られます。20% のホールドアウトは公正な評価のためだけです）。
%[text] **まとめ**
%[text] - 解釈可能な係数を得るために、線形回帰前に特徴量を z スコア標準化します。
%[text] - 学習誤差は常に楽観的です。常に交差検証で汎化性能を推定しましょう。
%[text] - Random Forest は線形モデルが見逃す非線形効果を捕えます。しかし、小規模データでは線形モデルが勝ることもあります。
%[text] - QSAR 予測を信頼する前に、適用可能ドメインを必ず確認しましょう。 \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
