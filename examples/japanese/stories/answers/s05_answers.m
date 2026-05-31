%[text] # S05 解答: 未知物質の同定チャレンジ
%[text] s05_unknown_substance.m の「やってみよう」演習の参考解答。
%[text] まず s05_unknown_substance.m を実行してから、このファイルで答え合わせをすること。
addpath(genpath("src"));
emk.setup.initPython();
logInfo("S05 解答: セットアップ完了");

%[text] ---- 未知化合物（証拠品 FC-2026-0004）--------------------------------------
UNKNOWN_SMILES = "CN(C)CCCN1c2ccccc2CCc2ccccc21";   % イミプラミン
mol_unk  = emk.mol.fromSmiles(UNKNOWN_SMILES);
fp_unk   = emk.fingerprint.morgan(mol_unk);
desc_unk = emk.descriptor.calculate(mol_unk, ...
    ["MolWt","LogP","TPSA","NumHDonors","NumHAcceptors","NumRotatableBonds","RingCount"]);

%[text] ---- 参照データベース -------------------------------------------------------
DB_FILE = "data/list/forensic_challenge.csv";
refDB   = readtable(DB_FILE, "TextType", "string");
if isnumeric(refDB.is_drug)
    refDB.is_drug = double(refDB.is_drug);
else
    refDB.is_drug = double(str2double(refDB.is_drug));
end

valid_ref = arrayfun(@(s) emk.mol.isValid(s), refDB.SMILES);
ref_smi_valid   = refDB.SMILES(valid_ref);
ref_names_valid = fillmissing(refDB.Name(valid_ref),     "constant", "");
ref_cat_valid   = fillmissing(refDB.Category(valid_ref), "constant", "");
ref_drug_valid  = refDB.is_drug(valid_ref);

ref_fps = cell(1, sum(valid_ref));
for i = 1:sum(valid_ref)
    ref_fps{i} = emk.fingerprint.morgan(emk.mol.fromSmiles(ref_smi_valid(i)));
end
id_result = emk.similarity.rankBy(fp_unk, ref_fps);
%[text] ---
%%
%[text] ## やってみよう 1: Ro5 チェック・CNS 薬物クラス仮説

d = desc_unk;
roPass = d.MolWt < 500 && d.LogP <= 5 && d.NumHDonors <= 5 && d.NumHAcceptors <= 10;
logInfo("未知物質 Ro5 合格: %d", roPass);
logInfo("  MW=%.1f  LogP=%.2f  TPSA=%.1f  HBD=%d  HBA=%d", ...
    d.MolWt, d.LogP, d.TPSA, d.NumHDonors, d.NumHAcceptors);

%[text] A1: はい -- Ro5 基準 4 つすべて合格。
%[text] A2: TPSA < 20 A^2（CNS 透過閾値 90 A^2 を大幅に下回る）。
%[text]     抗うつ薬（TCA）・抗精神病薬・抗ヒスタミン薬はいずれも
%[text]     意図的に低 TPSA に設計されている。
%[text] A3: 環 3 つ + HBD=0（NH/OH なし）は三環系アミンのフィンガープリント。
%[text]     三環系抗うつ薬（イミプラミン・アミトリプチリン・クロミプラミン）は
%[text]     すべてこの構造的青写真を共有する。
%[text] A4: 描画した構造はジベンザゼピン（ベンゼン環 2 個が 7 員窒素含有環に縮合）と
%[text]     ジメチルアミノプロピルチェーンを示す。典型的な TCA（三環系抗うつ薬）構造。
%%
%[text] ## やってみよう 2: ユニークカテゴリ・処方薬数

cats = unique(refDB.Category);
cats = cats(~ismissing(cats));
cats(cats == "") = [];
logInfo("ユニークカテゴリ (%d): %s", numel(cats), strjoin(cats, " | "));

drugNames = refDB.Name(refDB.is_drug == 1);
logInfo("処方薬 (%d 件):", numel(drugNames));
for i = 1:numel(drugNames)
    logInfo("  %s", drugNames(i));
end

%[text] A: 典型的なカテゴリ: stimulant（覚醒剤）, analgesic（鎮痛薬）, flavor（香料）,
%[text]    solvent（溶媒）, food_acid（食品酸）, sugar（糖）, vitamin（ビタミン）, drug（薬）。
%[text]    処方薬は約 30 件以上。
%[text]    環 3 つ・低 TPSA から最も可能性が高いカテゴリは "drug"（CNS 系）。
%%
%[text] ## やってみよう 3: 最多 ON ビット化合物・未知物質の複雑さ

nBitsEach = cellfun(@(fp) sum(emk.fingerprint.toArray(fp)), ref_fps);
[maxBits, maxIdx] = max(nBitsEach);
logInfo("最複雑（ON ビット最多）: %s (%d ビット)", ...
    ref_names_valid(maxIdx), maxBits);

nBits_unk = sum(emk.fingerprint.toArray(fp_unk));
logInfo("未知物質 ON ビット: %d / %d", nBits_unk, 2048);

mol_eth = emk.mol.fromSmiles("CCO");
fp_eth  = emk.fingerprint.morgan(mol_eth);
logInfo("エタノール ON ビット: %d（単純な参考値）", sum(emk.fingerprint.toArray(fp_eth)));

%[text] A: データベース内で最複雑な化合物は通常、最多の環と官能基を持つ。
%[text]    縮合環 3 つの未知物質は単純分子（エタノール: ~3〜5 ビット）より
%[text]    多くの ON ビットを持つが、複雑なマクロライドよりは少ない。
%%
%[text] ## やってみよう 4: タニモトスコア・スコアギャップ・分布

T_top1 = id_result.Scores(1);
T_top2 = id_result.Scores(2);
logInfo("上位候補タニモト: %.4f", T_top1);
logInfo("スコアギャップ（1 位 - 2 位）: %.4f", T_top1 - T_top2);

figure("Name", "同定検索 -- スコア分布");
histogram(id_result.Scores, 20, "BinEdges", linspace(0, 1, 21));
xlabel("未知物質へのタニモト");  ylabel("件数");
title("法科学的同定 -- スコア分布");
xline(0.65, "r--", "類似（0.65）");
xline(1.00, "g-",  "完全一致");

nAbove03 = sum(id_result.Scores > 0.3);
logInfo("T > 0.3 の化合物: %d / %d", nAbove03, numel(id_result.Scores));

%[text] 2 位一致化合物を比較のために描画する
idx2  = id_result.Indices(2);
mol2  = emk.mol.fromSmiles(ref_smi_valid(idx2));
figure("Name", "順位 2 の一致");
emk.viz.draw2d(mol2, Title=sprintf("順位 2: %s T=%.4f", ref_names_valid(idx2), T_top2));

%[text] A: T_top1 = 1.0（データベースの完全一致 -- 同じ化合物）。
%[text]    1 位と 2 位のスコアギャップが大きい（> 0.2）場合、同定の証拠が強化される。
%[text]    2 位ヒットは通常構造的に近い TCA（例: ノルトリプチリン・デシプラミン）。
%[text]    ほとんどのデータベース化合物は T < 0.2 で、ECFP4 の特異性を示す。
%%
%[text] ## やってみよう 5: 上位候補との物性比較・Ro5

TOP_IDX = id_result.Indices(1);
logInfo("上位候補: %s", ref_names_valid(TOP_IDX));
logInfo("  カテゴリ: %s  処方薬: %d", ref_cat_valid(TOP_IDX), ref_drug_valid(TOP_IDX));

smi_top  = ref_smi_valid(TOP_IDX);
mol_top  = emk.mol.fromSmiles(smi_top);
d_top    = emk.descriptor.calculate(mol_top, ["MolWt","LogP","TPSA","NumHDonors","NumHAcceptors"]);
logInfo("MW 差: %.2f Da（未知 %.1f vs 候補 %.1f）", ...
    abs(desc_unk.MolWt - d_top.MolWt), desc_unk.MolWt, d_top.MolWt);

t_top = struct2table(d_top);
ro5   = emk.filter.lipinski(t_top);
logInfo("上位候補 Ro5: 合格=%d  違反=%d", ro5.Pass_Ro5, ro5.Violations_Ro5);

%[text] A: MW 差 ~0（完全一致で同一化合物を確認）。
%[text]    上位候補は処方薬（IsDrug = 1）。
%[text]    処方箋なしにこの物質を所持することは臨床的・法的に重大。
%[text]    Ro5 合格 -- 薬物様経口製剤。
%%
%[text] ## やってみよう 6: イミプラミン vs ノルトリプチリン・カフェイン

smi_nor = refDB.SMILES(refDB.Name == "NORTRIPTYLINE");
if ~isempty(smi_nor)
    mol_nor = emk.mol.fromSmiles(smi_nor(1));
    fp_nor  = emk.fingerprint.morgan(mol_nor);
    T_nor   = emk.similarity.tanimoto(fp_unk, fp_nor);
    logInfo("イミプラミン vs ノルトリプチリン: T=%.4f", T_nor);
else
    logWarn("NORTRIPTYLINE がデータベースに見つからない");
end

mol_caf = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
fp_caf  = emk.fingerprint.morgan(mol_caf);
T_caf   = emk.similarity.tanimoto(fp_unk, fp_caf);
logInfo("イミプラミン vs カフェイン: T=%.4f", T_caf);

%[text] A: イミプラミン vs ノルトリプチリン: T ~0.26（中程度 -- T > 0.7 ではない）。
%[text]    どちらも TCA でも、イミプラミンはジベンズ[b,f]アゼピンコア
%[text]    （7 員環に窒素含有）、ノルトリプチリンはジベンズ[b,f]シクロヘプタジエンコア
%[text]    （環内に C=C、環窒素なし）を持つ。
%[text]    リング系の違いが ECFP4 環境を十分に変えるため類似度は中程度にとどまる。
%[text]    「同じ薬物クラス」は T >> 0.5 を意味しない -- スキャフォールドの違いが ECFP4 を支配する。
%[text]    同定の確認には T=1.0 が必要; T=0.26 では不十分。
figure("Name", "未知物質（イミプラミン）");
emk.viz.draw2d(mol_unk, Title="未知物質（イミプラミン）");
if exist("mol_nor", "var")
    figure("Name", "ノルトリプチリン（代謝物）");
    emk.viz.draw2d(mol_nor, Title="ノルトリプチリン（代謝物）");
end
%%
%[text] ## やってみよう 7: クロスバリデーション -- MACCS キー
%[text] MACCS キーで同定パイプラインを実行し、ECFP4 の結果を確認する.
fp_unk_maccs = emk.fingerprint.maccs(mol_unk);
maccs_fps    = cellfun(@(smi) emk.fingerprint.maccs(emk.mol.fromSmiles(smi)), ...
    cellstr(ref_smi_valid), "UniformOutput", false);
res_maccs    = emk.similarity.rankBy(fp_unk_maccs, maccs_fps);
logInfo("MACCS 1 位: %s T=%.4f", ref_names_valid(res_maccs.Indices(1)), res_maccs.Scores(1));

%[text] A: MACCS キーでも同じ化合物が 1 位に来るはず。
%[text]    異なるフィンガープリントタイプで一致した同定は
%[text]    法科学的証拠としてより強力。
%[text]
%[text] フルコナゾールを別の未知物質として試す:
SMILES2  = "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F";
mol2     = emk.mol.fromSmiles(SMILES2);
fp2      = emk.fingerprint.morgan(mol2);
res2     = emk.similarity.rankBy(fp2, ref_fps);
logInfo("フルコナゾール 1 位: %s T=%.4f", ref_names_valid(res2.Indices(1)), res2.Scores(1));

%[text] 新規化合物（データベースに存在しない）:
SMILES3 = "CN(C)c1ccc(C2=CC=Cc3ccc(CN(C)C)cc32)cc1";
mol3    = emk.mol.fromSmiles(SMILES3);
fp3     = emk.fingerprint.morgan(mol3);
res3    = emk.similarity.rankBy(fp3, ref_fps);
logInfo("新規化合物 1 位 T=%.4f -- データベース外（T < 1.0）", res3.Scores(1));

%[text] A: フルコナゾール（データベースに存在）: T=1.0（完全一致）。
%[text]    上位非自身ヒットはトリアゾール環を共有 -- ECFP4 のトリアゾール断片感度を反映。
%[text]    新規化合物: 最大 T < 1.0 -- データベースに存在しないことを確認。
%[text]    これは「データベースギャップ」問題: 法科学データベースは
%[text]    新規向精神性物質（NPS）に追いつくために継続的更新が必要。

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
