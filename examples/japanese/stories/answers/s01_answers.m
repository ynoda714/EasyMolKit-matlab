%[text] # S01 解答: カフェインの仲間を探せ
%[text] s01_find_caffeine_cousins.m の「やってみよう」演習の参考解答。
%[text] まず s01_find_caffeine_cousins.m を実行してデータベースとフィンガープリントを
%[text] 構築してから、このファイルで答え合わせをすること。
addpath(genpath("src"));
emk.setup.initPython();
logInfo("S01 解答: セットアップ完了");
%%
%[text] ## やってみよう 1: カフェインは BBB 透過基準を満たすか？

mol_caf  = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
desc_caf = emk.descriptor.calculate(mol_caf, ["TPSA", "LogP"]);

logInfo("カフェイン TPSA = %.1f A^2  (BBB 推奨: < 90 A^2)", desc_caf.TPSA);
logInfo("カフェイン LogP = %.2f       (中程度の親油性)",  desc_caf.LogP);

%[text] Q: カフェインは両条件を満たすか？
%[text] A: はい。TPSA ~62 A^2（< 90 -- 良好）、LogP ~-1.0（わずかに負 --
%[text]    非常に水溶性が高い）。
%[text]    LogP が低いにもかかわらず、カフェインは小さなサイズと
%[text]    平面状の芳香族構造のために BBB を効率よく通過する。
%[text]    LogP が負ということは水に溶けやすいことを意味する（コーヒー・お茶）。
%%
%[text] ## やってみよう 2: 最多分子カテゴリ・最重量分子は？

tbl = readtable("data/list/everyday_chemicals.csv", "TextType", "string");

%[text] 最多分子カテゴリ
cats        = unique(tbl.Category);
catCounts   = cellfun(@(c) sum(tbl.Category == c), cellstr(cats));
[~, iMax]   = max(catCounts);
logInfo("最多分子カテゴリ: %s (%d 件)", cats(iMax), catCounts(iMax));

%[text] 最重量分子
[maxMW, iHeavy] = max(tbl.MolecularWeight);
logInfo("最重量分子: %s  (MW = %.1f Da)", tbl.CommonName(iHeavy), maxMW);

%[text] A: 実行して出力を確認する。「flavor（香料）」が最多になることが多い。
%[text]    スクロース（砂糖）や脂質様化合物が最重量になる傾向がある。
%%
%[text] ## やってみよう 3: カフェインのモルガンフィンガープリントの ON ビット数は？

fp_caf = emk.fingerprint.morgan(mol_caf);
bits   = emk.fingerprint.toArray(fp_caf);
nOn    = sum(bits);

logInfo("カフェイン ECFP4 ON ビット: %d / %d  (密度 %.1f%%)", ...
    nOn, numel(bits), 100 * nOn / numel(bits));

%[text] A: 2048 ビット中 ON ビットは約 30〜50 個。
%[text]    正確な数はモルガン半径とビット数によって変わる。
%[text]    カフェインの二環式プリン系は多くの異なる局所環境を生成するため、
%[text]    単環や非環状分子よりも多くのビットが立つ。
%%
%[text] ## やってみよう 4: 1 位ヒット（カフェイン自身を除く）は？
%[text] データベースのフィンガープリントは s01 で計算済み（ワークスペース変数を使用）。
%[text] スタンドアロン実行の場合は再構築する:
tbl         = readtable("data/list/everyday_chemicals.csv", "TextType", "string");
nMols       = height(tbl);
all_fps     = cell(1, nMols);
valid       = false(1, nMols);
for i = 1:nMols
    if emk.mol.isValid(tbl.SMILES(i))
        all_fps{i} = emk.fingerprint.morgan(emk.mol.fromSmiles(tbl.SMILES(i)));
        valid(i)   = true;
    end
end
valid_fps   = all_fps(valid);
valid_names = tbl.CommonName(valid);

result = emk.similarity.rankBy(fp_caf, valid_fps);
logInfo("順位 1（自身）: %s  T=%.4f", valid_names(result.Indices(1)), result.Scores(1));
logInfo("順位 2（最良非自身）: %s  T=%.4f", ...
    valid_names(result.Indices(2)), result.Scores(2));

%[text] A: 順位 1 はカフェイン自身（T = 1.0）。
%[text]    順位 2 は THEOBROMINE（カカオ・チョコレートに含まれる）。
%[text]    テオブロミンは同じキサンチン（プリン）コアを持つが、
%[text]    N-メチル基が 3 個ではなく 2 個。タニモト ~0.53。
%[text]    天然物データベースとしては高い類似度。
%%
%[text] ## やってみよう 5: カフェインとテオブロミンの構造的違い

mol_theo = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)NC(=O)N2C");
figure("Name", "カフェイン（N-CH3 x 3）", "Position", [100 500 440 380]);
emk.viz.draw2d(mol_caf,  Title="カフェイン（N-CH3 × 3）");
figure("Name", "テオブロミン（N-CH3 x 2）", "Position", [560 500 440 380]);
emk.viz.draw2d(mol_theo, Title="テオブロミン（N-CH3 × 2）");

logInfo("カフェイン  分子式: C8H10N4O2  (N-メチル基 3 個: N1, N3, N7)");
logInfo("テオブロミン分子式: C7H8N4O2  (N-メチル基 2 個: N3, N7; N1 は NH)");

%[text] A: 唯一の構造的違いは N1 の置換基:
%[text] - カフェイン:     N1 にメチル基（N-CH3）
%[text] - テオブロミン:   N1 に水素（N-H）
%[text]    この 1 つの置換の違いでフィンガープリントからビット環境が 1 つ失われ、
%[text]    タニモトが 1.0 から ~0.53 に低下する。

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
