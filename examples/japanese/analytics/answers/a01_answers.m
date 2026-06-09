%% A01 解答: PCA による化学空間マッピング
% a01_chemical_space_pca.m の「やってみよう」演習の参照解答。
% このファイルは自己完結しています。CSV 読み込みから記述子計算まで
% 独立して実行できます（main ファイルの事前実行は不要です）。

addpath(genpath("src"));
emk.setup.initPython();
logInfo("A01 解答: セットアップ完了");
%%
%[text] ## やってみよう 1: 最大の絶対範囲を持つ記述子はどれか?
%[text] Q: 範囲の大きさが非標準化 PCA でなぜ問題になるのか説明してみましょう。
DATA_FILE  = "data/list/everyday_chemicals.csv";
DESC_NAMES = ["MolWt", "LogP", "TPSA", "NumHDonors", "NumHAcceptors", ...
              "RingCount", "NumRotatableBonds", "FractionCSP3", "HeavyAtomCount"];

rawTbl = readtable(DATA_FILE, TextType="string");
nMols  = height(rawTbl);
mols   = cell(1, nMols);
valid  = false(1, nMols);
for k = 1:nMols
    try
        mols{k} = emk.mol.fromSmiles(rawTbl.SMILES(k));
        valid(k) = true;
    catch
    end
end
validIdx = find(valid);
descMat  = nan(numel(validIdx), numel(DESC_NAMES));
for j = 1:numel(validIdx)
    s = emk.descriptor.calculate(mols{validIdx(j)}, DESC_NAMES);
    for d = 1:numel(DESC_NAMES)
        descMat(j, d) = s.(DESC_NAMES(d));
    end
end

descRange = max(descMat) - min(descMat);
[maxRange, iMax] = max(descRange);
logInfo("最大範囲: %s  (%.2f)", DESC_NAMES(iMax), maxRange);

%[text] 各記述子の範囲一覧:
for d = 1:numel(DESC_NAMES)
    logInfo("  %-22s  range = %.3f", DESC_NAMES(d), descRange(d));
end

%[text] **解答:** このデータセットでは MolWt が最大範囲（～300 単位）を持ちます。
%[text] 一方、FractionCSP3 は 0–1 の範囲にすぎません。
%[text] 非標準化 PCA では、MolWt のスケールだけで第 1 主成分を支配してしまいます。
%[text] 化学的な情報量とは無関係に、単純に数値スケールで結果が左右されてしまいます。
%[text] z スコア標準化により、各記述子がデータセット内の分散に比例して寄与できるようになります。
%%
%[text] ## やってみよう 2: 標準化を確認; 相関行列を計算

descMean = mean(descMat, 1);
descStd  = std(descMat, 0, 1);
constMask = descStd < 1e-12;
descMat2  = descMat(:, ~constMask);
descMean2 = descMean(~constMask);
descStd2  = descStd(~constMask);
activeDes = DESC_NAMES(~constMask);
Z = (descMat2 - descMean2) ./ descStd2;

colMeans = mean(Z);
colStds  = std(Z);
logInfo("Z の列平均（~0 であるべき）: max |mean| = %.2e", max(abs(colMeans)));
logInfo("Z の列標準偏差（~1 であるべき）: max |std-1| = %.2e", max(abs(colStds - 1)));

%[text] 相関行列の計算:
C = corr(Z);
%[text] 最も相関の高いペア（対角以外）を検出:
C_offdiag = C - diag(diag(C));
[maxCorr, linIdx] = max(abs(C_offdiag(:)));
[r, c] = ind2sub(size(C), linIdx);
logInfo("最も相関のあるペア: %s & %s  (r = %.3f)", ...
    activeDes(r), activeDes(c), C(r, c));

%[text] **解答:** Z の各列の平均は ~1e-16（数値的にゼロ）、標準偏差は 1 になります。
%[text] MolWt と HeavyAtomCount は強い相関（r > 0.95）を持ちます。
%[text] 分子が大きいほど原子数が多いため、この 2 つはほぼ同じ情報を担っています。
%[text] そのため、どちらかを除去しても PCA 結果はほとんど変わらず、冗長性を削減できます。
%%
%[text] ## やってみよう 3: >= 80% の分散に必要な主成分数は?

[~, ~, ~, ~, explained] = pca(Z);

cumVar = cumsum(explained);
nFor80 = find(cumVar >= 80, 1);
logInfo(">= 80%% の分散に必要な PC 数: %d  (cumulative at PC%d: %.1f%%)", ...
    nFor80, nFor80, cumVar(nFor80));

for p = 1:min(5, numel(explained))
    logInfo("  PC%d: %.1f%%  (cumulative: %.1f%%)", p, explained(p), cumVar(p));
end

%[text] **解答:** この 9 記述子・30 分子のデータセットでは、通常 3 つの PC で >= 80% の分散をカバーできます。
%[text] カットオフ閾値を設ける場合は、エルボー法よりも累積閾値法（>= 80%）の方が明確です。
%[text] PCA は最大 min(N, D) = min(30, 9) = 9 個の意味ある PC しか生成できません。
%[text] 30×9 の標準化行列のランクは最大 9 になります。
%%
%[text] ## やってみよう 4: PC 空間で原点から最も遠い分子を特定

[coeff, score, ~, ~, explained] = pca(Z);
validNames = rawTbl.CommonName(validIdx);

dist2 = score(:, 1).^2 + score(:, 2).^2;
[~, outerIdx] = sort(dist2, "descend");

logInfo("PC 原点から最も遠い上位 5 分子:");
for k = 1:5
    logInfo("  %d. %-20s  d = %.3f  (PC1=%.2f, PC2=%.2f)", ...
        k, validNames(outerIdx(k)), sqrt(dist2(outerIdx(k))), ...
        score(outerIdx(k), 1), score(outerIdx(k), 2));
end

%[text] **解答:** スクロース（分子量 342、OH 基多数）およびアスコルビン酸などは、
%[text] 大きなサイズと高い極性により PC1 方向に遠く位置する傾向があります。
%[text] エタノールなどの小さなアルコールは原点付近にクラスタを形成します。
%[text] スクロースは正の PC1（大きい/極性）、リモネンやイブプロフェンは負の PC2（親油性、低 TPSA）側に位置します。
%%
%[text] ## やってみよう 5: PC 負荷量の解釈; HeavyAtomCount の削除
%[text] PC1 と PC2 のローディング上位寄与因子を表示:
[~, ord1] = sort(abs(coeff(:, 1)), "descend");
[~, ord2] = sort(abs(coeff(:, 2)), "descend");
logInfo("PC1 主要寄与因子: %s (loading=%.3f)", ...
    activeDes(ord1(1)), coeff(ord1(1), 1));
logInfo("PC2 主要寄与因子: %s (loading=%.3f)", ...
    activeDes(ord2(1)), coeff(ord2(1), 2));

%[text] HeavyAtomCount を除去して PCA を再実行:
noHAC = activeDes ~= "HeavyAtomCount";
Z_noHAC = Z(:, noHAC);
[~, ~, ~, ~, exp2] = pca(Z_noHAC);
logInfo("HeavyAtomCount あり PC1 説明分散: %.1f%%", explained(1));
logInfo("HeavyAtomCount なし PC1 説明分散: %.1f%%", exp2(1));
logInfo("PC1 分散変化: %+.1f pp", exp2(1) - explained(1));

%[text] **解答:** PC1 は NumHAcceptors、TPSA、MolWt など「サイズ/極性」の記述子が主導します。
%[text] PC2 は LogP と FractionCSP3 など「親油性」の記述子が主導します。
%[text] HeavyAtomCount （MolWt と共線性が高い）を除去しても、
%[text] PC1 の説明分散の変化は 2% 未満です。冗長な記述子は独自の情報をほとんど追加しないことが確認できます。
logInfo("A01 解答完了。");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
