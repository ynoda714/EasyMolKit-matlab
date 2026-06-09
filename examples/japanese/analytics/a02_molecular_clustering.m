%[text] # A02: 分子クラスタリング — 構造的類似性による化学品のグループ化
%[text] EasyMolKit アナリティクス — レイヤー 3
%[text] 
%[text] 天然物化学者が食品や家庭用品から単離した30種類の化合物があります。
%[text] 生物活性試験の前に「構造的に似たグループ」を把握することで、冗長なテストを省き、コストを削減できます。
%[text] このスクリプトでは、Morgan ECFP4フィンガープリントとタニモト距離を用いて、階層的クラスタリングとk-meansでグループ化する方法を学びます。
%[text] 
%[text] **このチュートリアルで学べること**
%[text] - クラスタリングに記述子ではなくフィンガープリントを使用する理由を理解する
%[text] - `emk.similarity.matrix()` でタニモト行列を構築する方法を学ぶ
%[text] - タニモト類似度を階層的手法の距離指標に変換する方法を理解する
%[text] - Statistics Toolbox の `linkage()` と `dendrogram()` を使用する方法を学ぶ
%[text] - MATLAB の `kmeans()` と `silhouette()` を用いてクラスタの品質を評価する方法を学ぶ
%[text] - クラスタリング結果を化学的文脈で解釈する方法を学ぶ \
%[text] 
%[text] **前提条件**
%[text] - F03（フィンガープリント）とF04（類似度）の完了
%[text] - 推奨: A01（化学空間PCA）でコンテキストを把握すること
%[text] - Statistics and Machine Learning Toolbox（linkage、kmeans、silhouette）の使用
%[text] - インターネット接続は不要 \
%[text] 
%[text] 推定所要時間: 30〜45分
%[text] 
%[text] **データ:**
%[text] `data/list/everyday_chemicals.csv` — 30種類の一般分子（PubChem CC0）
%[text] 
%[text] **参考文献**
%[text] - Willett P, Barnard JM, Downs GM (1998) Chemical similarity searching.  J Chem Inf Comput Sci 38:983-996. doi:10.1021/ci9800211
%[text] - Rogers D & Hahn M (2010) Extended-connectivity fingerprints.  J Chem Inf Model 50:742-754. doi:10.1021/ci100050t
%[text] - Ward JH (1963) Hierarchical grouping to optimise an objective function.  J Am Stat Assoc 58:236-244. doi:10.1080/01621459.1963.10500845
%[text] - Kaufman L & Rousseeuw PJ (1990) Finding Groups in Data: An Introduction  to Cluster Analysis. Wiley. (silhouette coefficient, Chapter 2) \
%[text] 
%[text] 実行方法: Ctrl+Enter でセクションを1つずつ実行
%%
%[text] ## セクション 0: 環境セットアップ
logSection("A02", "セクション 0: セットアップ", "アナリティクス L3");
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

%[text] メイン処理の前に、Python/RDKit プロセスをウォームアップします
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;
%%
%[text] ## セクション 1: 全分子への Morgan フィンガープリント計算
%[text] 
%[text] ### コンセプト: 構造クラスタリングにフィンガープリントを使う理由
%[text] 記述子ベースの手法（PCA、回帰）は物理化学的特性（サイズ、極性など）を捉えます。ただし、2つの分子が類似した特性を持っていても異なるスキャフォールドを持つ場合や、その逆の場合に、構造的類似性を見逃すことがあります。
%[text] 
%[text] フィンガープリントは局所的な部分構造（各原子周辺の円形環境）の有無をバイナリビットベクトルとしてエンコードします。同じスキャフォールドを持つ2つの分子は、置換基が異なっていても多くの共通ビットを持ちます。これは記述子では捉えられない性質です。
%[text] 
%[text] **Morgan / ECFP フィンガープリント:**
%[text] - 半径 2（ECFP4）: 各原子から2結合先までの部分構造を捉える（直径4が名称ECFP**4**の由来）
%[text] - 2048 ビット（NBits=2048）: 衝突確率を最小化
%[text] - 分子の向きに依存しない（3D フィンガープリントとは異なる） \
%[text] 
%[text] **記述子ベースのクラスタリングが適する場合:** 溶解性や吸収性など物理的特性が類似性を定義する場合（製剤／ADMET目的）
%[text] **フィンガープリントベースのクラスタリングが適する場合:** スキャフォールド／部分構造の類似性が重要な場合（SAR、創薬化学）
%[text] 
logSection("A02", "セクション 1: 全分子への Morgan フィンガープリント計算", "アナリティクス L3");
DATA_FILE = "data/list/everyday_chemicals.csv";
FP_RADIUS = 2;     % Morgan radius 2 = ECFP4
FP_NBITS  = 2048;

rawTbl = readtable(DATA_FILE, TextType="string");
nRaw   = height(rawTbl);
logInfo("%s から %d 個の分子を読み込みました", DATA_FILE, nRaw);
%[text] SMILES を解析してフィンガープリントを計算します
fps   = cell(1, nRaw);
valid = false(1, nRaw);
for k = 1:nRaw
    try
        mol     = emk.mol.fromSmiles(rawTbl.SMILES(k));
        fps{k}  = emk.fingerprint.morgan(mol, Radius=FP_RADIUS, NBits=FP_NBITS);
        valid(k) = true;
    catch ME
        logWarn("%s をスキップ: %s", rawTbl.CommonName(k), ME.message);
    end
end

validIdx  = find(valid);
validFps  = fps(validIdx);
molNames  = rawTbl.CommonName(validIdx);
molCats   = rawTbl.Category(validIdx);
nMols     = numel(validIdx);
logInfo("ECFP4 フィンガープリント計算完了: %d / %d 分子（%d ビット）", ...
    nMols, nRaw, FP_NBITS);
%[text] **💡 観察ポイント 1 — カフェインのフィンガープリント密度を調べる**
%[text] カフェインの on-ビット（セットされているビット）は何個あるかを確認しましょう。
%[text] 以下のコードを参考にしてください:
%[text] fp\_caf = validFps{molNames == "caffeine"};
%[text] bits = `emk.fingerprint.toArray(fp_caf)`;
%[text] sum(bits)   % on-ビット数
%[text] sum(bits) / FP\_NBITS  % 密度
%[text] on-ビットが多い分子は必ず他の分子とのタニモト係数が高くなるかを考えてみましょう。
%[text] ヒント: Tanimoto = |A and B| / |A or B|
%[text] on-ビットが多い分子は分母が大きくなりますが、それが Tanimoto にとって有利か不利かを考えてみましょう。
%[text] 
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: タニモト類似度行列の構築
%[text] セクション 1 では、全ての分子の ECFP4 フィンガープリントを計算しました。
%[text] このセクションでは、全ての N×N 分子ペアの類似度を計算し、行列として表現します。
%[text] 
%[text] ### コンセプト: 分子距離の代替指標としてのタニモト類似度
%[text] バイナリフィンガープリントに対するタニモト（Jaccard）係数: 
%[text]{"align":"center"} $T(A,B) = |A \\cap B| / |A \\cup B| = c / (a + b - c)$
%[text]{"align":"center"} （ただし、$a = |A|$、$b = |B|$、$c$ = 共通の on-ビット数）
%[text] 
%[text] **性質:** $T = 1.0$ は同一フィンガープリント、$T = 0.0$ は共通の on-ビットなしを示します。
%[text] $d = 1 - T$（Jaccard 距離）は、三角不等式を満たす正式な距離指標です。
%[text] 
%[text] **業界標準のタニモト閾値の目安:** $T \\geq 0.85$ はほぼ同一、$T \\geq 0.70$ は密接なアナログ、$T \\geq 0.40$ は同じ化学クラス、$T \< 0.40$ は構造的に多様です。
%[text] 
logSection("A02", "セクション 2: タニモト類似度行列の構築", "アナリティクス L3");
simMat  = emk.similarity.matrix(validFps, Metric="tanimoto");  % N x N double
distMat = 1 - simMat;                                          % 距離行列
%[text] 対角要素を正確に 0 に設定します（Python IPC の数値ノイズ対策）。
distMat(logical(eye(nMols))) = 0;

logInfo("類似度行列: %d x %d (範囲: %.3f - %.3f)", ...
    nMols, nMols, min(simMat(:)), max(simMat(~logical(eye(nMols)))));
%[text] \-- 類似度行列をヒートマップとして可視化します --
figure("Name", "A02 タニモト類似度行列");
imagesc(simMat);
colormap(parula);
colorbar;
caxis([0 1]);
title("全対全タニモト類似度（ECFP4）");
xlabel("分子インデックス");
ylabel("分子インデックス");
%[text] カテゴリ区切り線を重ねて、見やすくするためにカテゴリでソートします。
[sortedCats, sortOrd] = sort(molCats);
logInfo("%d x %d 類似度行列のヒートマップを表示しました", nMols, nMols);
%[text] **💡 観察ポイント 2 — 類似度行列を読み取る**
%[text] 全ペアの平均タニモト類似度を確認しましょう。
%[text] ヒント: offDiag = simMat(~logical(eye(nMols))); mean(offDiag)
%[text] 最も類似した分子ペア（非同一で最高タニモト）を見つけましょう。
%[text] ヒント: \[r,c\] = find(simMat == max(simMat(~logical(eye(nMols)))));
%[text] molNames(r(1)), molNames(c(1))
%[text] そのペアが化学的に意味があるか考えてみましょう。
%[text] カテゴリでソートして行列を再プロットし、観察しましょう（sortOrd を使用）。
%[text] 同じカテゴリの分子が明るいブロックにまとまるか確認しましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: 階層的クラスタリングとデンドログラム
%[text] 類似度行列の作成が完了しました。
%[text] 次に、1 - タニモト係数を距離として使用し、階層的クラスタリングで分子をボトムアップでグループ化します。
%[text] 
%[text] ### コンセプト: 階層的クラスタリング（Ward 連結法）
%[text] 階層的クラスタリングはボトムアップ方式でツリー（デンドログラム）を構築します。
%[text] 1. 開始: 各分子が独自のクラスタを形成します。
%[text] 2. 最も近い2つのクラスタを統合します（連結基準に基づいて決定）。
%[text] 3. 全ての分子が1つのクラスタになるまでこのプロセスを繰り返します。 \
%[text] 
%[text] **連結基準:** `"single"` — 最近傍への距離（チェーン効果あり） / `"complete"` — 最遠傍への距離 / `"average"` — 全ペア距離の平均（UPGMA） / `"ward"` — クラスタ内分散を最小化（Ward 1963）
%[text] 
%[text] Ward 連結法は、ケモインフォマティクスでコンパクトで均等なサイズのクラスタを生成するために好まれます。クラスタ内誤差二乗和を最小化する点で k-means と目的関数を共有し、類似した結果が得られることが多いです。
%[text] 
%[text] **デンドログラムの読み方:** 結合の高さは、2つのクラスタが統合されたときの非類似度を示します。高さ $h$ での水平線（「カット」）が $k$ 個のクラスタを定義します。低い高さで統合されるグループは非常に類似しています。
%[text] 
logSection("A02", "セクション 3: 階層的クラスタリングとデンドログラム", "アナリティクス L3");
%[text] 対称距離行列を縮約ベクトル（上三角）に変換します（MATLAB linkage() 関数に必要）。
condDist = squareform(distMat, "tovector");

Z_tree = linkage(condDist, "ward");
%[text] \-- デンドログラムのプロット --
figure("Name", "A02 階層的クラスタリング デンドログラム");
D = dendrogram(Z_tree, 0, ...
    Labels=cellstr(molNames), ...
    Orientation="left", ...
    ColorThreshold=0.65 * max(Z_tree(:, 3)));  % 最大高さの 65%
title("Ward 連結法デンドログラム — 日用化学品（ECFP4）");
xlabel("Ward 距離");
set(gca, FontSize=7);
logInfo("デンドログラムを表示しました");
%[text] \-- 最大高さの 65% で自動カット（約 7〜8 クラスタ）--
CUT_HEIGHT   = 0.65 * max(Z_tree(:, 3));
clusterLabels_hier = cluster(Z_tree, Cutoff=CUT_HEIGHT, Criterion="distance");
nClusters_hier = max(clusterLabels_hier);
logInfo("階層的カット: 閾値=%.3f -> %d クラスタ", CUT_HEIGHT, nClusters_hier);

for c = 1:nClusters_hier
    members = molNames(clusterLabels_hier == c);
    logInfo("  クラスタ H%d (%d 分子): %s", c, numel(members), ...
        strjoin(members, ", "));
end
%[text] **💡 観察ポイント 3 — カットの高さと連結法を変えてみる**
%[text] カットの高さを最大の 20% に下げると、または 60% に上げると、
%[text] クラスタ数はどのように変わるかを確認しましょう。
%[text] 以下のコードを参考にしてください:
%[text] c1 = cluster(Z\_tree, Cutoff=0.2\*max(Z\_tree(:,3)), Criterion="distance");
%[text] c2 = cluster(Z\_tree, Cutoff=0.6\*max(Z\_tree(:,3)), Criterion="distance");
%[text] \[max(c1), max(c2)\]
%[text] 連結法を "average" または "complete" に変更してみましょう。
%[text] 以下のコードを参考にしてください: Z2 = linkage(condDist, "complete"); figure; dendrogram(Z2,0);
%[text] Ward 連結法とツリーの形がどう異なるかを観察しましょう。
%[text] このデータセットで、どちらの方法がよりバランスのとれたクラスタを生成するかを考えてみましょう。
%[text] グルコースとフルクトースはどちらも C6H12O6（構造異性体）です。
%[text] デンドログラムで早い段階（低い高さ）で統合されるかを確認しましょう。
%[text] これは構造異性体に対する ECFP4 の特徴について何を示唆しているかを考えてみましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: フィンガープリントビットベクトルへの K-means クラスタリング
%[text] 階層的クラスタリングは、ツリー全体を一度に構築し、後からカット位置を決める手法です。
%[text] 一方、k-means は k を固定して、直接 k 個のクラスタに分割するアプローチです。
%[text] 
%[text] ### コンセプト: ハミング距離による K-means
%[text] K-means は $N$ 点を $k$ クラスタに分割するために以下を繰り返します:
%[text] 1. 各点を最近傍の重心に割り当てます。
%[text] 2. クラスタの平均として重心を再計算します。
%[text] 3. 割り当てが安定するまで繰り返します。 \
%[text] 
%[text] 標準の k-means はユークリッド距離を使用しますが、バイナリフィンガープリントにはハミング距離（異なるビットの割合）が適しています。これは `kmeans()` の `Distance="hamming"` オプションでサポートされています。
%[text] 
%[text] $k$ **の選択:** ドメイン知識（期待される化学クラス数）、エルボー法（クラスタ内二乗和 vs $k$）、シルエット分析（平均シルエット幅が最適 $k$ でピーク）を使用します。
%[text] 
%[text] フィンガープリントを数値ビット行列（$N \\times$ FP\_NBITS）に変換します。
logSection("A02", "セクション 4: フィンガープリントビットベクトルへの K-means クラスタリング", "アナリティクス L3");
bitMat = zeros(nMols, FP_NBITS, "logical");
for j = 1:nMols
    bitMat(j, :) = emk.fingerprint.toArray(validFps{j});
end
bitMat_double = double(bitMat);
logInfo("フィンガープリントビット行列: %d x %d", size(bitMat_double));
%[text] \-- エルボー法: k = 2..8 で評価します --
K_MAX = 8;
wcss  = zeros(1, K_MAX);     % クラスタ内二乗和

rng(42);   % 再現性のため
for k = 2:K_MAX
    [~, ~, sumd] = kmeans(bitMat_double, k, ...
        Distance="hamming", Replicates=5, Display="off");
    wcss(k) = sum(sumd);
end

figure("Name", "A02 K-Means エルボー曲線");
plot(2:K_MAX, wcss(2:K_MAX), "-o", Color=[0.2 0.5 0.8], LineWidth=1.5);
xlabel("クラスタ数 k");
ylabel("クラスタ内二乗和");
title("エルボー法 — ECFP4 ビットベクトルへの K-Means");
grid on;
logInfo("エルボー曲線を表示しました（急な落下が平坦になる点を探す）");
%[text] **💡 観察ポイント 4 — エルボーで最適 k を探す**
%[text] WCSS 曲線でどの k でエルボー（収益逓減点）が現れるかを確認しましょう。それがこのデータセットの最適な k です。なぜここでユークリッド距離でなくハミング距離を使うのかを考えてみましょう。
%[text] ヒント: フィンガープリントベクトルの各要素は 0 か 1 です。
%[text] 2 ビット間のユークリッド距離は sqrt(0)=0 または sqrt(1)=1 で、ハミング距離と同じです。しかし、ビット密度が変化するバイナリベクトルでは、ユークリッド距離は一般化しにくいです。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: k 選択のためのシルエット分析
%[text] k-means 法では、クラスタ数 k を事前に指定する必要があります。
%[text] シルエット分析を用いて、データに基づいた最適な k の選び方を学びましょう。
%[text] 
%[text] ### コンセプト: シルエット係数
%[text] 各分子 $i$ について、シルエット係数 $s(i)$ は、その分子が自身のクラスにどれだけ適合しているかを示します。これは最近傍クラスとの比較で測定されます: $s(i) = (b(i) - a(i)) / \\max(a(i), b(i))$
%[text] （ここで $a(i)$ はクラス内平均距離、$b(i)$ は最近傍クラスへの平均距離です）
%[text] 
%[text] **解釈:** $s \\approx 1$ の場合、分子はよくクラスタ化されています。$s \\approx 0$ の場合、分子は2つのクラスタの境界にあります。$s \\approx -1$ の場合、誤ってクラスタに割り当てられている可能性があります。
%[text] 
%[text] 全分子の平均シルエットはクラスタ品質の全体的な指標です。平均シルエットを最大化する $k$ が最適なクラスタ数となります。
%[text] 
%[text] ハミング距離を用いてシルエットを計算します（bitMat\_double に Distance="hamming" を指定）。
logSection("A02", "セクション 5: k 選択のためのシルエット分析", "アナリティクス L3");
silScores = zeros(1, K_MAX);
rng(42);
for k = 2:K_MAX
    labels_k = kmeans(bitMat_double, k, ...
        Distance="hamming", Replicates=5, Display="off");
    sil = silhouette(bitMat_double, labels_k, "hamming");
    silScores(k) = mean(sil);
end

figure("Name", "A02 シルエットスコア");
plot(2:K_MAX, silScores(2:K_MAX), "-s", Color=[0.8 0.3 0.2], LineWidth=1.5);
xlabel("クラスタ数 k");
ylabel("平均シルエットスコア");
title("シルエット分析 — ECFP4 ビットベクトルへの K-Means");
grid on;
logInfo("シルエットプロットを表示しました");

[bestSil, bestK_idx] = max(silScores(2:K_MAX));
bestK = bestK_idx + 1;
logInfo("シルエットによる最適 k: k = %d  (平均シルエット = %.3f)", bestK, bestSil);

%[text] 最適な k を用いて最終的な k-means を実行します。
rng(42);
finalLabels_km = kmeans(bitMat_double, bestK, ...
    Distance="hamming", Replicates=10, Display="off");
logInfo("K-Means (k=%d) 最終クラスタ割り当て:", bestK);
for c = 1:bestK
    members = molNames(finalLabels_km == c);
    logInfo("  クラスタ KM%d (%d 分子): %s", c, numel(members), ...
        strjoin(members, ", "));
end
%[text] **💡 観察ポイント 5 — k-means と階層的クラスタリングを比較する**
%[text] k-means クラスタ（KM\*）と階層的クラスタ（H\*）を比較し、どの分子が同じクラスタにまとまるか確認しましょう。2つの手法間でクラスタが入れ替わる分子があるかを確認しましょう。同じ k で階層的クラスタに `silhouette()` を実行し、
%[text] hier\_labels = cluster(Z\_tree, MaxClust=bestK);
%[text] sil\_hier = silhouette(bitMat\_double, hier\_labels, "hamming");
%[text] mean(sil\_hier)
%[text] このデータセットでどちらの手法がより高いシルエットスコアを出すか確認しましょう。化学カテゴリラベル（興奮剤、鎮痛剤など）はクラスタ化の「正解」となります。Adjusted Rand Index を計算するか、混同行列を確認して、
%[text] \[~, catIdx\] = ismember(molCats, unique(molCats));
%[text] k-means クラスタがどの程度化学カテゴリを再現できているかを確認しましょう。
%[text] （ヒント: これは教師なし学習であり、完全な再現は期待しないでください）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: 類似度行列上へのクラスタ可視化
%[text] クラスタ割り当てが完了しました。
%[text] 最後に類似度行列をクラスタ順に並べ替え、クラスタリングの品質を視覚的に確認しましょう。
%[text] 
%[text] ### コンセプト: ソート済みヒートマップ — ブロック構造の発見
%[text] クラスタ割り当てに基づいて類似度行列の行と列を並び替えると、高類似度の分子が対角線沿いに明るい正方形ブロックを形成します。**良いクラスタリング:** タイトなブロック、オフブロック値が低い / **悪いクラスタリング:** 拡散しており、明確なブロックがない
%[text] 
%[text] 階層的クラスタでソート（ヒートマップでは k-means より視覚的に美しい）
logSection("A02", "セクション 6: 類似度行列上へのクラスタ可視化", "アナリティクス L3");
[~, sortByCluster] = sort(clusterLabels_hier);
simMatSorted = simMat(sortByCluster, sortByCluster);
namesSorted  = molNames(sortByCluster);
catsSorted   = molCats(sortByCluster);

figure("Name", "A02 類似度行列（クラスタ順ソート）");
imagesc(simMatSorted);
colormap(parula);
colorbar;
caxis([0 1]);
title("類似度行列（階層的クラスタ順ソート）");
xlabel("分子（クラスタ順）");
ylabel("分子（クラスタ順）");

%[text] 目盛りラベルを追加（省略形）
tick_labels = cellfun(@(n) n(1:min(8,end)), cellstr(namesSorted), ...
    UniformOutput=false);
xticks(1:nMols); xticklabels(tick_labels); xtickangle(90);
yticks(1:nMols); yticklabels(tick_labels);
set(gca, FontSize=6);
%[text] クラスタ境界線を描く
boundaries = [0; find(diff(clusterLabels_hier(sortByCluster))); nMols] + 0.5;
hold on;
for b = boundaries'
    xline(b, "w-", LineWidth=1);
    yline(b, "w-", LineWidth=1);
end
hold off;
logInfo("クラスタ境界付きのソート済み類似度行列を表示しました");
%[text] **💡 観察ポイント 6 — ソート済みヒートマップで代表分子を選ぶ**
%[text] 対角のブロックは一貫して明るい（高類似度）か確認しましょう。
%[text] 薄いブロック（低いクラスタ内類似度）はあるか確認しましょう。
%[text] これはより多くのクラスタが必要か、またはその分子がタイトなグループに属さない真の外れ値である可能性を示します。
%[text] 階層的ラベルの代わりに k-means ラベルでソートした場合を確認しましょう。
%[text] どちらがより明確なブロック構造を示すか確認しましょう。
%[text] 実用的な応用: コスト削減のために各クラスタから 1 つの代表分子をバイオアッセイに送るとしたら、どれを選びますか？
%[text] 一般的な戦略: クラスタ重心に最も近い分子を選ぶ（全クラスタ員への中央値タニモトが最も高い）。
%[text] 
%[text] **まとめ**
%[text] - ECFP4 フィンガープリントは局所的な部分構造をビットベクトルとしてエンコードします。
%[text] - タニモト距離（$1 - T$）が階層的クラスタリングの入力として使われます。
%[text] - Ward 連結法はコンパクトでバランスのとれたクラスタを生成します。
%[text] - シルエット分析はデータ駆動で最適な $k$ を選ぶ手段を提供します。
%[text] - ソート済みヒートマップはブロック構造を視覚的品質チェックとして示します。 \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
