%[text] # F02 Answers: Calculate Molecular Properties (Descriptors)
%[text] This is the reference answer for the exercise problem in `f02_calculate_properties.m`.
%[text] Please run `f02_calculate_properties.m` first, then check your answers here.
addpath(genpath("src"));
emk.setup.initPython();
%%
%[text] ## Answer E1: Lipinski Check for Paracetamol

mol_para = emk.mol.fromSmiles("CC(=O)NC1=CC=C(C=C1)O");
d = emk.descriptor.calculate(mol_para, ["MolWt","LogP","NumHDonors","NumHAcceptors"]);
logInfo("Paracetamol Descriptors: MW=%.2f, LogP=%.2f, HBD=%d, HBA=%d", ...
    d.MolWt, d.LogP, d.NumHDonors, d.NumHAcceptors);

pass_mw  = d.MolWt  <= 500;
pass_lp  = d.LogP   <= 5;
pass_hbd = d.NumHDonors    <= 5;
pass_hba = d.NumHAcceptors <= 10;
logInfo("Ro5: MW=%d, LogP=%d, HBD=%d, HBA=%d  (1=PASS)", ...
    pass_mw, pass_lp, pass_hbd, pass_hba);
%[text] Expected values: MW ~151.2, LogP ~0.91, HBD=2, HBA=2 → All 4 criteria PASS
%%
%[text] ## Answer E2: Maximum and Minimum of MW and LogP, and Number of Ro5 Violations

data = readtable("data/list/everyday_chemicals.csv", TextType="string");
mols = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(data.SMILES), "UniformOutput", false);
desc = emk.descriptor.batchCalculate(mols, ["MolWt","LogP"]);
desc.Name = data.CommonName;

[~, iMwMax] = max(desc.MolWt);  [~, iMwMin] = min(desc.MolWt);
[~, iLpMax] = max(desc.LogP);   [~, iLpMin] = min(desc.LogP);
logInfo("Max MW: %s (%.2f g/mol)", desc.Name(iMwMax), desc.MolWt(iMwMax));
logInfo("Min MW: %s (%.2f g/mol)", desc.Name(iMwMin), desc.MolWt(iMwMin));
logInfo("Max LogP: %s (%.2f)",     desc.Name(iLpMax), desc.LogP(iLpMax));
logInfo("Min LogP: %s (%.2f)",     desc.Name(iLpMin), desc.LogP(iLpMin));
logInfo("MW > 500: %d cases / LogP > 5: %d cases", sum(desc.MolWt > 500), sum(desc.LogP > 5));
%%
%[text] ## Answer E3: Descriptor Table for 4 Drugs

smiles_list = ["CC(=O)Oc1ccccc1C(=O)O", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", ...
               "CC(=O)NC1=CC=C(C=C1)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"];
drug_names  = ["Aspirin", "Ibuprofen", "Paracetamol", "Caffeine"];

mols_drug = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(smiles_list), "UniformOutput", false);
t = emk.descriptor.batchCalculate(mols_drug, ...
    ["MolWt","LogP","TPSA","NumHDonors","NumHAcceptors","NumRotatableBonds"]);
t.Properties.RowNames = cellstr(drug_names);
disp(t);

%[text] ---
%[text] Try It Answers
%[text] ---
%%
%[text] ## Try It 1: Ibuprofen Descriptors and Lipinski Ro5 Check

mol_ibup = emk.mol.fromSmiles("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
d_ibup = emk.descriptor.calculate(mol_ibup, ...
    ["MolWt","LogP","NumHDonors","NumHAcceptors"]);

logInfo("Ibuprofen: MW=%.2f, LogP=%.2f, HBD=%d, HBA=%d", ...
    d_ibup.MolWt, d_ibup.LogP, d_ibup.NumHDonors, d_ibup.NumHAcceptors);

pass_mw  = d_ibup.MolWt          <= 500;
pass_lp  = d_ibup.LogP           <= 5;
pass_hbd = d_ibup.NumHDonors     <= 5;
pass_hba = d_ibup.NumHAcceptors  <= 10;
logInfo("Ro5: MW=%d, LogP=%d, HBD=%d, HBA=%d  (1=PASS)", ...
    pass_mw, pass_lp, pass_hbd, pass_hba);
%[text] MW ~206.3, LogP ~3.52, HBD = 1, HBA = 1 → All 4 criteria PASS
%[text] Ibuprofen is a small lipophilic NSAID and fits well within the Ro5 space.
%%
%[text] ## Try It 2: Identify Compounds Outside the Ro5 Zone

data = readtable("data/list/everyday_chemicals.csv", TextType="string");
mols_ec = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(data.SMILES), ...
    "UniformOutput", false);
desc_ec = emk.descriptor.batchCalculate(mols_ec, ["MolWt","LogP","TPSA", ...
    "NumHDonors","NumHAcceptors"]);
desc_ec.Name = data.CommonName;

figure("Name", "MW vs LogP (with Ro5 lines)", "Position", [100 100 650 500]);
scatter(desc_ec.LogP, desc_ec.MolWt, 60, desc_ec.TPSA, "filled");
colorbar;
xlabel("LogP  (lipophilicity)");
ylabel("Molecular Weight  [g/mol]");
title("Everyday Chemicals -- Ro5 Boundary Check");
hold on;
xline(5,   "r--", "Ro5: LogP=5",  "LabelVerticalAlignment","bottom");
yline(500, "r--", "Ro5: MW=500",  "LabelHorizontalAlignment","left");
ylim([0, 550]);
hold off;

%[text] Identifying molecules outside the Ro5 MW/LogP zone.
outside = desc_ec( desc_ec.MolWt > 500 | desc_ec.LogP > 5, ["Name","MolWt","LogP"]);
logInfo("Molecules outside Ro5 MW/LogP zone: %d cases", height(outside));
disp(outside);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---