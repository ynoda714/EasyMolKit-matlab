%[text] # A02 解答: 分子クラスタリング
%[text] a02_molecular_clustering.m の「やってみよう」演習の参照解答。
%[text] 最初に a02_molecular_clustering.m を実行してワークスペース変数
%[text] (validFps, molNames, molCats, simMat, distMat, Z_tree, bitMat_double,
%[text] を利用可能にしてください。
addpath(genpath("src"));
emk.setup.initPython();
logInfo("A02 解答: セットアップ完了");

%[text] ── ワークスペース変数の再構築（スタンドアロン実行） ──
DATA_FILE = "data/list/everyday_chemicals.csv";
FP_RADIUS = 2;
FP_NBITS  = 2048;
rawTbl    = readtable(DATA_FILE, TextType="string");
fps       = cell(1, height(rawTbl));
valid     = false(1, height(rawTbl));
for k = 1:height(rawTbl)
    try
        mol = emk.mol.fromSmiles(rawTbl.SMILES(k));
        fps{k} = emk.fingerprint.morgan(mol, Radius=FP_RADIUS, NBits=FP_NBITS);
        valid(k) = true;
    catch
    end
end
validIdx  = find(valid);
validFps  = fps(validIdx);
molNames  = rawTbl.CommonName(validIdx);
molCats   = rawTbl.Category(validIdx);
nMols     = numel(validIdx);

simMat  = emk.similarity.matrix(validFps, Metric="tanimoto");
distMat = 1 - simMat;
distMat(logical(eye(nMols))) = 0;

bitMat = zeros(nMols, FP_NBITS, "logical");
for j = 1:nMols
    bitMat(j, :) = emk.fingerprint.toArray(validFps{j});
end
bitMat_double = double(bitMat);

condDist = squareform(distMat, "tovector");
Z_tree   = linkage(condDist, "ward");
%%
%[text] ## やってみよう 1: カフェインのオンビット; オンビット増加は Tanimoto を改善するか?

fp_caf = validFps{molNames == "caffeine"};
bits   = emk.fingerprint.toArray(fp_caf);
nOn    = sum(bits);
logInfo("カフェイン ECFP4: %d / %d ビット設定（密度 %.1f%%）", ...
    nOn, FP_NBITS, 100 * nOn / FP_NBITS);

%[text] 全分子のオンビット密度を比較する
densities = sum(bitMat_double, 2) / FP_NBITS;
logInfo("オンビット密度: 最小=%.1f%%  平均=%.1f%%  最大=%.1f%%", ...
    min(densities)*100, mean(densities)*100, max(densities)*100);

%[text] **解答:** カフェインは 2048 ビット中 ~40–50 ビットがオン（密度 ~2%）です。
%[text] 二環式プリン骨格が単純な非環式分子よりも多くの局所環境を生成します。
%[text] **疎なベクトルでは、オンビットが多いほど Tanimoto は低下します**: T = c/(a+b−c)。
%[text] 分子 A が 200 ビット、B が 10 ビットを持ち、B の全 10 ビットが A と重なる場合でも
%[text] T = 10/(200+10−10) = 10/200 = 0.05 と非常に低い値になります。
%[text] 多くのオンビットを持つ分子は、同様にオンビットが多い分子以外とは
%[text] 低い Tanimoto スコアを示します。
%%
%[text] ## やってみよう 2: 平均ペアワイズ類似度; 最も類似したペア

offDiag = simMat(~logical(eye(nMols)));
logInfo("平均ペアワイズ Tanimoto: %.3f", mean(offDiag));
logInfo("中央値ペアワイズ Tanimoto: %.3f", median(offDiag));

%[text] 最も類似したペアを探索する
C = simMat;
C(logical(eye(nMols))) = 0;
[maxVal, linIdx] = max(C(:));
[r, c] = ind2sub(size(C), linIdx);
logInfo("最も類似したペア: %s & %s  (T = %.3f)", ...
    molNames(r), molNames(c), maxVal);

%[text] カテゴリ別ソート済み類似度行列の表示
[sortedCats, sortOrd] = sort(molCats);
simSorted = simMat(sortOrd, sortOrd);
figure("Name", "A02 Sorted by Category");
imagesc(simSorted); colormap(parula); colorbar; caxis([0 1]);
title("カテゴリ別ソート済み類似度行列");
xticks(1:nMols); xticklabels(cellstr(molNames(sortOrd))); xtickangle(90);
yticks(1:nMols); yticklabels(cellstr(molNames(sortOrd)));
set(gca, FontSize=6);

%[text] **解答:** 30 種の日用化学品の多様なセットでは、平均ペアワイズ Tanimoto は
%[text] 通常 0.05–0.12 となり、構造的な多様性が確認されます。
%[text] 最も類似したペアは、リモネン＆カルボン（どちらも環式テルペン）または
%[text] カフェイン＆テオブロミン（どちらもメチルキサンチン類）となることが多いです。
%[text] カテゴリ別にソートすると、構造的に一貫したカテゴリ（糖類、覚醒剤）では
%[text] 対角線沿いに明るいブロックが現れます。
%%
%[text] ## やってみよう 3: カット高さの影響; 結合法比較

c1 = cluster(Z_tree, Cutoff=0.2*max(Z_tree(:,3)), Criterion="distance");
c2 = cluster(Z_tree, Cutoff=0.6*max(Z_tree(:,3)), Criterion="distance");
logInfo("最大高さの 20%% でカット -> %d クラスター", max(c1));
logInfo("最大高さの 60%% でカット -> %d クラスター", max(c2));

%[text] 平均連結法・完全連結法との比較
Z_avg  = linkage(condDist, "average");
Z_comp = linkage(condDist, "complete");

figure("Name", "A02 Linkage Comparison");
tiledlayout(1, 3);
nexttile; dendrogram(Z_tree, 0, Orientation="left"); title("Ward");    xlabel("Distance"); set(gca, FontSize=6);
nexttile; dendrogram(Z_avg,  0, Orientation="left"); title("Average"); xlabel("Distance"); set(gca, FontSize=6);
nexttile; dendrogram(Z_comp, 0, Orientation="left"); title("Complete");xlabel("Distance"); set(gca, FontSize=6);
sgtitle("結合法比較 -- 日用化学品");

%[text] グルコースとフルクトースの類似度
gIdx = find(molNames == "glucose");
fIdx = find(molNames == "fructose");
if ~isempty(gIdx) && ~isempty(fIdx)
    logInfo("グルコース-フルクトース Tanimoto: %.3f", simMat(gIdx, fIdx));
end

%[text] **解答:** カットを最大高さの 20% に下げると、多くの小さなクラスター（シングルトンも）が生じます。
%[text] 60% に上げると、ほとんどの分子が 3〘5 個の大きなクラスターに統合されます。
%[text] 平均連結法は Ward より若干コンパクト性が低いですが、外れ値への感度も低いです。
%[text] 完全連結法は最大径のクラスターを強制するため、非常に異なる結果になることがあります。
%[text] グルコースとフルクトース（構造異性体、C6H12O6）の Tanimoto は中程度（~0.3–0.4）です。
%[text] ECFP4 はデフォルトで立体化学を考慮しないため、環状と開鎖型の互変異性体が
%[text] 多くの部分構造を共有することが理由です。
%%
%[text] ## やってみよう 4: エルボー法 — どの k で収益逓減が始まるか?

K_MAX = 8;
wcss  = zeros(1, K_MAX);
rng(42);
for k = 2:K_MAX
    [~, ~, sumd] = kmeans(bitMat_double, k, ...
        Distance="hamming", Replicates=5, Display="off");
    wcss(k) = sum(sumd);
end

%[text] 差分で最大降下点（エルボー）を探す
drops = diff(wcss(2:K_MAX));
[~, elbowRel] = min(drops);   % smallest delta = where drop plateaus
elbowK = elbowRel + 2;        % offset: index 1 -> k=3
logInfo("エルボーの近似 k = %d  (WCSS 低下プラトー)", elbowK);

%[text] **解答:** この 30 分子セットでは、エルボーは通常 k = 4〘6 にあります。
%[text] Hamming 距離はバイナリビットベクトルに自然な指標で、
%[text] 異なる位置の割合 (a+b−2c)/(a+b) をカウントします。
%[text] これは Dice 類似度の補数に相当します。
%%
%[text] ## やってみよう 5: k-means と階層的クラスタリングの比較; 階層的のシルエット

rng(42);
[bestSilScore, bestK_idx] = deal(0, 0);
silScores = zeros(1, K_MAX);
for k = 2:K_MAX
    labels_k = kmeans(bitMat_double, k, Distance="hamming", Replicates=5, Display="off");
    sil = silhouette(bitMat_double, labels_k, "hamming");
    silScores(k) = mean(sil);
    if silScores(k) > bestSilScore
        bestSilScore = silScores(k);
        bestK_idx    = k;
    end
end
bestK = bestK_idx;
rng(42);
finalLabels_km = kmeans(bitMat_double, bestK, Distance="hamming", Replicates=10, Display="off");

%[text] bestK での階層的クラスタリング
hier_labels = cluster(Z_tree, MaxClust=bestK);
sil_hier    = silhouette(bitMat_double, hier_labels, "hamming");
sil_km      = silhouette(bitMat_double, finalLabels_km, "hamming");

logInfo("K-Means シルエット (k=%d): %.3f",       bestK, mean(sil_km));
logInfo("階層的シルエット (k=%d): %.3f",  bestK, mean(sil_hier));
logInfo("優勝: %s", ...
    string(ternary_(mean(sil_km) >= mean(sil_hier), "K-Means", "Hierarchical")));

%[text] カテゴリ回収率の確認
[~, catIdx] = ismember(molCats, unique(molCats));
logInfo("カテゴリラベル vs k-means (k=%d): crosstab() で比較可能", bestK);

%[text] **解答:** 両手法は、キサンチン類（カフェイン・テオブロミン）、テルペン類（リモネン・カルボン）、
%[text] 糖類（スクロース・グルコース・フルクトース）のレベルでは通常一致します。
%[text] この小規模データセットのシルエットスコアは控えめ（~0.2–0.4）です。
%[text] 完全なスコアを得るには非常に密なクラスターが必要です。
%[text] K-Means と階層的クラスタリングは、境界上の分子については異なる結果になる場合があります。
%%
%[text] ## やってみよう 6: クラスター代表選択戦略

CUT_HEIGHT     = 0.4 * max(Z_tree(:, 3));
clusterLabels_hier = cluster(Z_tree, Cutoff=CUT_HEIGHT, Criterion="distance");
nClusters_hier = max(clusterLabels_hier);

logInfo("クラスター代表（最高クラスター内平均類似度）:");
for c = 1:nClusters_hier
    idx = find(clusterLabels_hier == c);
    if numel(idx) == 1
        rep = idx;
    else
        subSim = simMat(idx, idx);
        meanSim = (sum(subSim, 2) - 1) / max(numel(idx) - 1, 1);
        [~, repLocal] = max(meanSim);
        rep = idx(repLocal);
    end
    logInfo("  クラスター H%d: 代表 = %s", c, molNames(rep));
end

%[text] **解答:** 代表分子とは、同一クラスター内の全分子に対して平均 Tanimoto が最も高い分子です。
%[text] これはフィンガープリント空間における「重心」と呼べます。
%[text] 実際のスクリーニングキャンペーンでは、クラスターごとに 1 つの代表分子を選んで生物活性を評価し、
%[text] 活性クラスター周辺の探索を拡大します。
logInfo("A02 解答完了。");

%[text] ローカルヘルパー（分岐を持つ無名関数の代替）
function out = ternary_(cond, a, b)
    if cond; out = a; else; out = b; end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
