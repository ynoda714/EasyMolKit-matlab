%[text] # F02 解答: 分子特性の計算
%[text] `f02_calculate_properties.m` の演習問題の参照解答です。
%[text] 先に `f02_calculate_properties.m` を実行してから、こちらで答え合わせをしてください。
addpath(genpath("src"));
emk.setup.initPython();
%%
%[text] ## 解答 E1: パラセタモールのリピンスキー確認

mol_para = emk.mol.fromSmiles("CC(=O)NC1=CC=C(C=C1)O");
d = emk.descriptor.calculate(mol_para, ["MolWt","LogP","NumHDonors","NumHAcceptors"]);
logInfo("パラセタモール記述子: MW=%.2f, LogP=%.2f, HBD=%d, HBA=%d", ...
    d.MolWt, d.LogP, d.NumHDonors, d.NumHAcceptors);

pass_mw  = d.MolWt  <= 500;
pass_lp  = d.LogP   <= 5;
pass_hbd = d.NumHDonors    <= 5;
pass_hba = d.NumHAcceptors <= 10;
logInfo("Ro5: MW=%d, LogP=%d, HBD=%d, HBA=%d  (1=PASS)", ...
    pass_mw, pass_lp, pass_hbd, pass_hba);
%[text] 期待値: MW ~151.2、LogP ~0.91、HBD=2、HBA=2 → 4 基準すべて PASS
%%
%[text] ## 解答 E2: MW・LogP の最大・最小と Ro5 違反数

data = readtable("data/list/everyday_chemicals.csv", TextType="string");
mols = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(data.SMILES), "UniformOutput", false);
desc = emk.descriptor.batchCalculate(mols, ["MolWt","LogP"]);
desc.Name = data.CommonName;

[~, iMwMax] = max(desc.MolWt);  [~, iMwMin] = min(desc.MolWt);
[~, iLpMax] = max(desc.LogP);   [~, iLpMin] = min(desc.LogP);
logInfo("MW 最大: %s (%.2f g/mol)", desc.Name(iMwMax), desc.MolWt(iMwMax));
logInfo("MW 最小: %s (%.2f g/mol)", desc.Name(iMwMin), desc.MolWt(iMwMin));
logInfo("LogP 最大: %s (%.2f)",     desc.Name(iLpMax), desc.LogP(iLpMax));
logInfo("LogP 最小: %s (%.2f)",     desc.Name(iLpMin), desc.LogP(iLpMin));
logInfo("MW > 500: %d 件 / LogP > 5: %d 件", sum(desc.MolWt > 500), sum(desc.LogP > 5));
%%
%[text] ## 解答 E3: 4 薬物の記述子テーブル

smiles_list = ["CC(=O)Oc1ccccc1C(=O)O", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", ...
               "CC(=O)NC1=CC=C(C=C1)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"];
drug_names  = ["アスピリン", "イブプロフェン", "パラセタモール", "カフェイン"];

mols_drug = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(smiles_list), "UniformOutput", false);
t = emk.descriptor.batchCalculate(mols_drug, ...
    ["MolWt","LogP","TPSA","NumHDonors","NumHAcceptors","NumRotatableBonds"]);
t.Properties.RowNames = cellstr(drug_names);
disp(t);

%[text] ---
%[text] やってみようの解答
%[text] ---
%%
%[text] ## やってみよう 1: イブプロフェン記述子とリピンスキー Ro5 確認

mol_ibup = emk.mol.fromSmiles("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
d_ibup = emk.descriptor.calculate(mol_ibup, ...
    ["MolWt","LogP","NumHDonors","NumHAcceptors"]);

logInfo("イブプロフェン: MW=%.2f, LogP=%.2f, HBD=%d, HBA=%d", ...
    d_ibup.MolWt, d_ibup.LogP, d_ibup.NumHDonors, d_ibup.NumHAcceptors);

pass_mw  = d_ibup.MolWt          <= 500;
pass_lp  = d_ibup.LogP           <= 5;
pass_hbd = d_ibup.NumHDonors     <= 5;
pass_hba = d_ibup.NumHAcceptors  <= 10;
logInfo("Ro5: MW=%d, LogP=%d, HBD=%d, HBA=%d  (1=PASS)", ...
    pass_mw, pass_lp, pass_hbd, pass_hba);
%[text] MW ~206.3、LogP ~3.52、HBD = 1、HBA = 1 → 4 基準すべて PASS
%[text] イブプロフェンは小さな親油性 NSAID で、Ro5 空間に十分収まります。
%%
%[text] ## やってみよう 2: Ro5 ゾーン外の化合物を特定

data = readtable("data/list/everyday_chemicals.csv", TextType="string");
mols_ec = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(data.SMILES), ...
    "UniformOutput", false);
desc_ec = emk.descriptor.batchCalculate(mols_ec, ["MolWt","LogP","TPSA", ...
    "NumHDonors","NumHAcceptors"]);
desc_ec.Name = data.CommonName;

figure("Name", "MW vs LogP（Ro5 ライン付き）", "Position", [100 100 650 500]);
scatter(desc_ec.LogP, desc_ec.MolWt, 60, desc_ec.TPSA, "filled");
colorbar;
xlabel("LogP  （親油性）");
ylabel("分子量  [g/mol]");
title("日用化学品 -- Ro5 境界確認");
hold on;
xline(5,   "r--", "Ro5: LogP=5",  "LabelVerticalAlignment","bottom");
yline(500, "r--", "Ro5: MW=500",  "LabelHorizontalAlignment","left");
ylim([0, 550]);
hold off;

%[text] Ro5 の MW/LogP ゾーン外の分子を特定します。
outside = desc_ec( desc_ec.MolWt > 500 | desc_ec.LogP > 5, ["Name","MolWt","LogP"]);
logInfo("Ro5 MW/LogP ゾーン外の分子: %d 件", height(outside));
disp(outside);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
