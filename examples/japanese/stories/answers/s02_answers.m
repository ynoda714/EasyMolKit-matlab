%[text] # S02 解答: 薬物フィルタ -- リピンスキーのルール・オブ・ファイブ
%[text] s02_drug_filter_lipinski.m の「やってみよう」演習の参考解答。
%[text] まず s02_drug_filter_lipinski.m を実行してから、このファイルで答え合わせをすること。
addpath(genpath("src"));
emk.setup.initPython();
logInfo("S02 解答: セットアップ完了");
%%
%[text] ## やってみよう 1: アスピリンとイブプロフェンの Ro5 チェック

mol_asp  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
mol_ibu  = emk.mol.fromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O");

for pair = {{"Aspirin",    mol_asp}, {"Ibuprofen", mol_ibu}}
    name = pair{1}{1};  mol = pair{1}{2};
    d = emk.descriptor.calculate(mol, ["MolWt","LogP","NumHDonors","NumHAcceptors"]);
    logInfo("%s: MW=%.1f  LogP=%.2f  HBD=%d  HBA=%d", ...
        name, d.MolWt, d.LogP, d.NumHDonors, d.NumHAcceptors);
end

%[text] A: アスピリン   -- MW ~180, LogP ~1.3, HBD=1, HBA=3 -- 違反 0（全合格）。
%[text]    イブプロフェン -- MW ~206, LogP ~3.1, HBD=1, HBA=1 -- 違反 0（全合格）。
%[text]    どちらも小さく吸収されやすい経口薬 -- Ro5 の予測どおり。
%%
%[text] ## やってみよう 2: FDA データセット内で最重量の薬物は？

rawTbl = readtable("data/list/fda_drugs.csv", "TextType", "string");

%[text] 列名を `emk.filter.lipinski` の期待する名前に変更
rawTbl = renamevars(rawTbl, ...
    ["MolecularWeight", "ALogP",  "HBondDonors",  "HBondAcceptors"], ...
    ["MolWt",           "LogP",   "NumHDonors",   "NumHAcceptors"]);

%[text] 数値列に変換
rawTbl.MolWt         = double(rawTbl.MolWt);
rawTbl.LogP          = double(rawTbl.LogP);
rawTbl.NumHDonors    = double(rawTbl.NumHDonors);
rawTbl.NumHAcceptors = double(rawTbl.NumHAcceptors);

[maxMW, idx] = max(rawTbl.MolWt);
logInfo("最重量薬物: %s  (MW = %.1f Da)", rawTbl.Name(idx), maxMW);
logInfo("SMILES: %s", rawTbl.SMILES(idx));

%[text] A: アンフォテリシン B（MW ~924 Da）-- Streptomyces nodosus 由来のポリエン抗真菌薬。
%[text]    経口投与ではなく静脈内投与される。Ro5 違反と整合。
%%
%[text] ## やってみよう 3: 最も頻繁に違反される Ro5 基準は？

drugTbl = emk.filter.lipinski(rawTbl);

vMW  = sum(rawTbl.MolWt           > 500);
vLP  = sum(rawTbl.LogP             > 5);
vHBD = sum(rawTbl.NumHDonors      > 5);
vHBA = sum(rawTbl.NumHAcceptors   > 10);

logInfo("違反数:");
logInfo("  MW  > 500: %d 薬物", vMW);
logInfo("  LogP > 5 : %d 薬物", vLP);
logInfo("  HBD > 5  : %d 薬物", vHBD);
logInfo("  HBA > 10 : %d 薬物", vHBA);

%[text] A: MW 違反が通常最多（大型の天然物由来薬物）、次いで HBA 違反。
%[text]    ほとんどの薬物は NH/OH 基が少ないため HBD > 5 は最も稀。
%%
%[text] ## やってみよう 4: 違反ちょうど 1 件の薬物; 緩和フィルタの合格数

violCounts = histcounts(drugTbl.Violations_Ro5, 0:6);
logInfo("違反 0 件の薬物: %d", violCounts(1));
logInfo("違反 1 件の薬物: %d", violCounts(2));
logInfo("違反 2 件の薬物: %d", violCounts(3));
logInfo("違反 3 件の薬物: %d", violCounts(4));

drugTblRelaxed = emk.filter.lipinski(rawTbl, MaxViolations=1);
nPassStrict  = sum(drugTbl.Pass_Ro5);
nPassRelaxed = sum(drugTblRelaxed.Pass_Ro5);
logInfo("厳格（違反 0 件）: %d / %d 合格", nPassStrict, height(rawTbl));
logInfo("緩和（違反 1 件）: %d / %d 合格", nPassRelaxed, height(rawTbl));
logInfo("追加回収薬物: %d 件", nPassRelaxed - nPassStrict);

%[text] A: 通常 ~14 薬物が違反ちょうど 1 件。
%[text]    MaxViolations=1 に緩和するとこれらが回収される（大型マクロライド系抗生物質や
%[text]    臨床実績のある天然物が多い）。
%%
%[text] ## やってみよう 5: 化学空間散布図 -- MW vs LogP

figure("Name", "化学空間: MW vs LogP", "Position", [100 100 560 440]);
passIdx = drugTbl.Pass_Ro5;
scatter(rawTbl.LogP(passIdx),  rawTbl.MolWt(passIdx),  40, "b", "filled", ...
    "DisplayName", "Ro5 合格");
hold on;
scatter(rawTbl.LogP(~passIdx), rawTbl.MolWt(~passIdx), 60, "r", "filled", ...
    "DisplayName", "Ro5 不合格");
xline(5,   "--k", "LogP=5");
yline(500, "--k", "MW=500");
xlabel("LogP");  ylabel("MW (Da)");
title("FDA 承認薬の化学空間");
legend("Location", "northwest");  grid("on");

%[text] A: はい -- ほとんどの FDA 承認薬は左下象限（MW < 500, LogP < 5）に集中する。
%[text]    この象限の外にある薬物は注射剤・外用剤・プロドラッグが典型。
%[text]    一部の Ro5 合格薬は LogP ~4〜5、MW ~450〜490 付近に集まる -- 薬物様だが境界に近い。
%%
%[text] ## やってみよう 6: ルール違反薬物を描画する
%[text] 不合格薬物の例 -- リファンピシンが典型的（あれば）
rifIdx = find(strcmpi(rawTbl.Name, "RIFAMPICIN"), 1);
if isempty(rifIdx)
    [~, rifIdx] = max(rawTbl.MolWt);   % フォールバック: 最重量薬物
end
mol_big = emk.mol.fromSmiles(rawTbl.SMILES(rifIdx));
figure("Name", "Ro5 違反薬物: " + rawTbl.Name(rifIdx), "Position", [100 100 440 380]);
emk.viz.draw2d(mol_big, Title=rawTbl.Name(rifIdx) + " (Ro5 違反)");
logInfo("ルール違反薬物: %s  MW=%.1f  違反数=%d", ...
    rawTbl.Name(rifIdx), rawTbl.MolWt(rifIdx), drugTbl.Violations_Ro5(rifIdx));

%[text] A: リファンピシンは大きなマクロ環を持つ（複雑な芳香族系）。
%[text]    ドナー/アクセプターを数えると多くの OH・NH 基が確認できる。
%[text]    それでも経口吸収される -- Ro5 はガイドラインであり絶対ルールではない。
%%
%[text] ## やってみよう 7: MaxViolations=2 フィルタ

drugTbl2 = emk.filter.lipinski(rawTbl, MaxViolations=2);
%[text] MaxViolations=2 で合格するが MaxViolations=1 では不合格の薬物:
newPass = drugTbl2.Pass_Ro5 & ~drugTblRelaxed.Pass_Ro5;
logInfo("違反ちょうど 2 件の薬物: %d", sum(newPass));
disp(rawTbl(newPass, ["Name", "MolWt", "LogP", "NumHDonors", "NumHAcceptors"]));

%[text] A: 違反 2 件の薬物には大型マクロライド（アジスロマイシン・エリスロマイシン）、
%[text]    タキサン（パクリタキセル）、免疫抑制薬（タクロリムス）が含まれる。
%[text]    これらはすべて静脈内投与 -- 天然物由来であることと複雑な薬理学を反映している。

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
