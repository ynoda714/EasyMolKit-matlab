%[text] # S04 解答: バーチャルスクリーニング入門
%[text] s04_virtual_screening.m の「やってみよう」演習の参考解答。
%[text] まず s04_virtual_screening.m を実行してから、このファイルで答え合わせをすること。
addpath(genpath("src"));
emk.setup.initPython();
logInfo("S04 解答: セットアップ完了");

%[text] ---- クエリ: イブプロフェン -----------------------------------------------
QUERY_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O";
mol_hit      = emk.mol.fromSmiles(QUERY_SMILES);
fp_hit       = emk.fingerprint.morgan(mol_hit);

%[text] ---- ライブラリ: FDA 承認薬 CSV -------------------------------------------
library      = readtable("data/list/fda_drugs.csv", "TextType", "string");
library.MolecularWeight = double(library.MolecularWeight);
library.ALogP           = double(library.ALogP);

valid   = arrayfun(@(s) emk.mol.isValid(s), library.SMILES);
lib_names_valid  = library.Name(valid);
lib_smiles_valid = library.SMILES(valid);
lib_mw_valid     = library.MolecularWeight(valid);
lib_alogp_valid  = library.ALogP(valid);

logInfo("ライブラリ: %d / %d 件が有効", sum(valid), height(library));

lib_fps_valid = cell(1, sum(valid));
for i = 1:sum(valid)
    lib_fps_valid{i} = emk.fingerprint.morgan(emk.mol.fromSmiles(lib_smiles_valid(i)));
end
vs_result = emk.similarity.rankBy(fp_hit, lib_fps_valid);
%[text] ---
%%
%[text] ## やってみよう 1: イブプロフェンは Ro5 を満たすか？

desc_hit = emk.descriptor.calculate(mol_hit, ["MolWt","LogP","NumHDonors","NumHAcceptors"]);
t_hit = struct2table(desc_hit);
ro5   = emk.filter.lipinski(t_hit);
logInfo("イブプロフェン Ro5: 合格=%d  違反=%d", ro5.Pass_Ro5, ro5.Violations_Ro5);

%[text] Q: 計算した LogP（中性形）と血中の見かけ LogP はなぜ異なるか？
%[text] A: LogP は中性（非イオン化）分子種で測定される。
%[text]    血液 pH 7.4 では、イブプロフェン（pKa ~4.4）は ~99% がカルボキシラートアニオンとしてイオン化。
%[text]    分布係数 D = LogP - log(1 + 10^(pH-pKa)) で生理的 pH では D << LogP。
%[text]    イオン化形は膜透過性が低いが、中性分子のわずかな割合が受動拡散を駆動するため
%[text]    イブプロフェンは依然として吸収される。
%%
%[text] ## やってみよう 2: ライブラリの概要

logInfo("ライブラリ内のユニーク薬物名: %d", numel(unique(library.Name)));

[maxMW, idx] = max(library.MolecularWeight);
logInfo("最重量化合物: %s  (MW = %.1f Da)", library.Name(idx), maxMW);

%[text] イブプロフェンの ALogP vs LogP（ライブラリに存在する場合）
ibIdx = find(strcmpi(lib_names_valid, "IBUPROFEN"), 1);
if ~isempty(ibIdx)
    logInfo("イブプロフェン ALogP（ライブラリ）: %.2f  vs  Wildman-Crippen LogP: %.2f", ...
        lib_alogp_valid(ibIdx), desc_hit.LogP);
end

%[text] A: ALogP（Ghose-Crippen）と Wildman-Crippen LogP は同じ物性を推定するが
%[text]    異なる原子寄与モデルを使用する。
%[text]    典型的な薬物様分子では両者の差は ~0.5 単位以内。
%%
%[text] ## やってみよう 3: イブプロフェン vs カフェインの ON ビット数

mol_caf = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
fp_caf  = emk.fingerprint.morgan(mol_caf);
fp_ibu  = emk.fingerprint.morgan(mol_hit);

nOn_caf = sum(emk.fingerprint.toArray(fp_caf));
nOn_ibu = sum(emk.fingerprint.toArray(fp_ibu));

logInfo("カフェイン  ON ビット: %d", nOn_caf);
logInfo("イブプロフェン ON ビット: %d", nOn_ibu);

%[text] A: カフェインの方が通常 ON ビット数が多い。
%[text]    縮合二環式プリン環系が多くの異なるモルガン環境を生成し、
%[text]    より多くのユニークなビット位置が立つ。
%[text]    イブプロフェンの単環ベンゼン環と短いアルキル鎖は環境の種類が少ない。
%[text]
%[text] 半径 3 の実験:
fp_r3   = emk.fingerprint.morgan(mol_hit, Radius=3);
nOn_r3  = sum(emk.fingerprint.toArray(fp_r3));
logInfo("イブプロフェン ECFP6（半径=3）ON ビット: %d  (vs ECFP4: %d)", nOn_r3, nOn_ibu);

%[text] A: 半径 3 では大きな近傍を考慮する。
%[text]    個々のビット衝突が減る（大きな環境はよりユニーク）可能性があるが、
%[text]    全体密度は通常似たか、わずかに減少する。
%%
%[text] ## やってみよう 4: 上位ヒット -- 他の NSAID も含まれるか？

REPORT_TOP = 15;
for k = 1:min(REPORT_TOP, numel(vs_result.Scores))
    idx = vs_result.Indices(k);
    logInfo("  %2d. %-30s  T=%.4f", k, lib_names_valid(idx), vs_result.Scores(k));
end

%[text] A: 上位 2 位（非自身）の期待値: FLURBIPROFEN（T~0.40）、KETOPROFEN（T~0.39）。
%[text]    両者ともアリルプロピオン酸コア（アリール--CH(CH3)--COOH）を共有する。
%[text]    これらの近接アナログでも T ~0.40 にとどまるのは、
%[text]    フルルビプロフェン/ケトプロフェンのフッ素や追加環が
%[text]    十分な局所 ECFP4 環境の変化をもたらすため。
%[text]
%[text] アスピリンを代替クエリとして使う:
mol_asp = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
fp_asp  = emk.fingerprint.morgan(mol_asp);
res_asp = emk.similarity.rankBy(fp_asp, lib_fps_valid);
logInfo("アスピリン上位 3 位: %s, %s, %s", ...
    lib_names_valid(res_asp.Indices(1)), ...
    lib_names_valid(res_asp.Indices(2)), ...
    lib_names_valid(res_asp.Indices(3)));
%[text] A: アスピリンの上位ヒットはイブプロフェンとは異なる --
%[text]    ECFP4 はアセトキシ安息香酸スキャフォールドを捉えるため。
%%
%[text] ## やってみよう 5: ヒットレポート -- MW 昇順ソート・最高 ALogP

rankVec  = (1:REPORT_TOP)';
nameVec  = strings(REPORT_TOP, 1);
scoreVec = zeros(REPORT_TOP, 1);
mwVec    = zeros(REPORT_TOP, 1);
logpVec  = zeros(REPORT_TOP, 1);
smiVec   = strings(REPORT_TOP, 1);
for k = 1:REPORT_TOP
    if k > numel(vs_result.Scores); break; end
    idx        = vs_result.Indices(k);
    nameVec(k) = lib_names_valid(idx);
    scoreVec(k)= vs_result.Scores(k);
    mwVec(k)   = lib_mw_valid(idx);
    logpVec(k) = lib_alogp_valid(idx);
    smiVec(k)  = lib_smiles_valid(idx);
end
hits_tbl = table(rankVec, nameVec, scoreVec, mwVec, logpVec, smiVec, ...
    'VariableNames', ["Rank","Name","Tanimoto","MW_Da","ALogP","SMILES"]);

logInfo("MW 昇順ソートのヒット:");
sorted_tbl = sortrows(hits_tbl, "MW_Da");
disp(sorted_tbl(:, 1:5));

[~, iMaxLogP] = max(hits_tbl.ALogP);
logInfo("最高 ALogP ヒット: %s  ALogP=%.2f  (薬物様 < 5)", ...
    hits_tbl.Name(iMaxLogP), hits_tbl.ALogP(iMaxLogP));

%[text] A: 上位ヒットはイブプロフェン MW（206 Da）付近に集まる傾向がある。
%[text]    ECFP4 類似度は構造類似度と相関し、このスキャフォールドクラスでは
%[text]    MW とも相関するため。
%%
%[text] ## やってみよう 6: スコア分布とカバレッジ

all_scores = vs_result.Scores;
nHigh   = sum(all_scores >= 0.65);
nModerate = sum(all_scores >= 0.40 & all_scores < 0.65);
logInfo("T >= 0.65 （類似）: %d / %d", nHigh, numel(all_scores));
logInfo("T >= 0.40 （中程度）: %d / %d", nModerate, numel(all_scores));

figure("Name", "スコア分布 -- イブプロフェン vs FDA", "Position", [100 500 560 400]);
histogram(all_scores, 30, "BinEdges", linspace(0, 1, 31));
xlabel("タニモト（ECFP4）");  ylabel("件数");
title("イブプロフェン vs FDA 承認薬 200 件");
xline(0.65, "--r", "類似（0.65）");

%[text] カフェインとの比較
mol_caf2 = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
fp_caf2  = emk.fingerprint.morgan(mol_caf2);
res_caf  = emk.similarity.rankBy(fp_caf2, lib_fps_valid);
figure("Name", "スコア分布 -- カフェイン vs FDA", "Position", [680 500 560 400]);
histogram(res_caf.Scores, 30, "BinEdges", linspace(0, 1, 31));
xlabel("タニモト（ECFP4）");  ylabel("件数");  title("カフェイン vs FDA 承認薬 200 件");

%[text] A: イブプロフェンは FDA ライブラリにある NSAID/プロフェン系がよく代表されるため
%[text]    最大タニモトが高い。
%[text]    カフェインはメチルキサンチン系覚醒剤で FDA 承認薬に構造的アナログが少なく、
%[text]    最大タニモトが低い。
%%
%[text] ## やってみよう 7: 上位 3 ヒットを描画して視覚的に比較する

for k = 1:3
    mol_k = emk.mol.fromSmiles(smiVec(k));
    emk.viz.draw2d(mol_k, Title=sprintf("ヒット %d: %s T=%.2f", k, nameVec(k), scoreVec(k)));
end

%[text] A: 上位 3 ヒット（フルルビプロフェン・ケトプロフェンなどのプロフェン類）は
%[text]    視覚的にアリルプロピオン酸コア（ベンゼン環 -- CH(CH3) -- COOH）を共有する。
%[text]    違いはアリール環上の置換基（F、2 環目の C=O）のみ。
%%
%[text] ## やってみよう 8: ECFP4 vs MACCS キー比較

fp_hit_maccs = emk.fingerprint.maccs(mol_hit);
lib_maccs    = cell(1, sum(valid));
for i = 1:sum(valid)
    lib_maccs{i} = emk.fingerprint.maccs(emk.mol.fromSmiles(lib_smiles_valid(i)));
end
res_maccs = emk.similarity.rankBy(fp_hit_maccs, lib_maccs);

logInfo("ECFP4 上位 3 位: %s, %s, %s", ...
    lib_names_valid(vs_result.Indices(1)), ...
    lib_names_valid(vs_result.Indices(2)), ...
    lib_names_valid(vs_result.Indices(3)));
logInfo("MACCS 上位 3 位: %s, %s, %s", ...
    lib_names_valid(res_maccs.Indices(1)), ...
    lib_names_valid(res_maccs.Indices(2)), ...
    lib_names_valid(res_maccs.Indices(3)));

%[text] A: ECFP4 と MACCS キーは類似するが同一ではないランキングを生成する。
%[text]    MACCS キー（166 ビットの固定サブ構造チェックリスト）は官能基の存在に焦点を当て、
%[text]    ECFP4 はより細かい原子環境を捉える。
%[text]    両者が一致する場合、そのランキングに信頼性が増す。

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
