%[text] # S05 Answers: Unknown Substance Identification Challenge
%[text] Reference answers for the "Try It" exercise in s05_unknown_substance.m.
%[text] First, run s05_unknown_substance.m, then use this file to check your answers.
addpath(genpath("src"));
emk.setup.initPython();
logInfo("S05 Answers: Setup Complete");

%[text] ---- Unknown Compound (Exhibit FC-2026-0004) --------------------------------------
UNKNOWN_SMILES = "CN(C)CCCN1c2ccccc2CCc2ccccc21";   % Imipramine
mol_unk  = emk.mol.fromSmiles(UNKNOWN_SMILES);
fp_unk   = emk.fingerprint.morgan(mol_unk);
desc_unk = emk.descriptor.calculate(mol_unk, ...
    ["MolWt","LogP","TPSA","NumHDonors","NumHAcceptors","NumRotatableBonds","RingCount"]);

%[text] ---- Reference Database -------------------------------------------------------
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
%[text] ## Let's Try 1: Ro5 Check & CNS Drug Class Hypothesis

d = desc_unk;
roPass = d.MolWt < 500 && d.LogP <= 5 && d.NumHDonors <= 5 && d.NumHAcceptors <= 10;
logInfo("Unknown substance Ro5 pass: %d", roPass);
logInfo("  MW=%.1f  LogP=%.2f  TPSA=%.1f  HBD=%d  HBA=%d", ...
    d.MolWt, d.LogP, d.TPSA, d.NumHDonors, d.NumHAcceptors);

%[text] A1: Yes -- All 4 Ro5 criteria passed.
%[text] A2: TPSA < 20 A^2 (well below CNS penetration threshold of 90 A^2).
%[text]     Antidepressants (TCA), antipsychotics, and antihistamines are all
%[text]     intentionally designed with low TPSA.
%[text] A3: 3 rings + HBD=0 (no NH/OH) is a fingerprint of tricyclic amines.
%[text]     Tricyclic antidepressants (Imipramine, Amitriptyline, Clomipramine) all
%[text]     share this structural blueprint.
%[text] A4: The drawn structure shows dibenzazepine (two benzene rings fused to a 7-membered nitrogen-containing ring) and
%[text]     a dimethylaminopropyl chain. Typical TCA (tricyclic antidepressant) structure.
%%
%[text] ## Let's Try 2: Unique Categories and Prescription Drug Count

cats = unique(refDB.Category);
cats = cats(~ismissing(cats));
cats(cats == "") = [];
logInfo("Unique categories (%d): %s", numel(cats), strjoin(cats, " | "));

drugNames = refDB.Name(refDB.is_drug == 1);
logInfo("Prescription drugs (%d items):", numel(drugNames));
for i = 1:numel(drugNames)
    logInfo("  %s", drugNames(i));
end

%[text] A: Typical categories: stimulant, analgesic, flavor,
%[text]    solvent, food_acid, sugar, vitamin, drug.
%[text]    There are about 30 or more prescription drugs.
%[text]    The most likely category from 3 rings and low TPSA is "drug" (CNS system).
%%
%[text] ## Let's Try 3: Compound with Most ON Bits & Complexity of Unknown Substance

nBitsEach = cellfun(@(fp) sum(emk.fingerprint.toArray(fp)), ref_fps);
[maxBits, maxIdx] = max(nBitsEach);
logInfo("Most complex (most ON bits): %s (%d bits)", ...
    ref_names_valid(maxIdx), maxBits);

nBits_unk = sum(emk.fingerprint.toArray(fp_unk));
logInfo("Unknown substance ON bits: %d / %d", nBits_unk, 2048);

mol_eth = emk.mol.fromSmiles("CCO");
fp_eth  = emk.fingerprint.morgan(mol_eth);
logInfo("Ethanol ON bits: %d (simple reference value)", sum(emk.fingerprint.toArray(fp_eth)));

%[text] A: The most complex compound in the database usually has the most rings and functional groups.
%[text]    An unknown substance with 3 fused rings has more ON bits than a simple molecule (Ethanol: ~3-5 bits),
%[text]    but fewer than a complex macrolide.
%%
%[text] ## Let's Try 4: Tanimoto Score, Score Gap, Distribution

T_top1 = id_result.Scores(1);
T_top2 = id_result.Scores(2);
logInfo("Top candidate Tanimoto: %.4f", T_top1);
logInfo("Score gap (1st - 2nd): %.4f", T_top1 - T_top2);

figure("Name", "Identification Search -- Score Distribution");
histogram(id_result.Scores, 20, "BinEdges", linspace(0, 1, 21));
xlabel("Tanimoto to Unknown");  ylabel("Count");
title("Forensic Identification -- Score Distribution");
xline(0.65, "r--", "Similarity (0.65)");
xline(1.00, "g-",  "Exact Match");

nAbove03 = sum(id_result.Scores > 0.3);
logInfo("Compounds with T > 0.3: %d / %d", nAbove03, numel(id_result.Scores));

%[text] Draw the 2nd matching compound for comparison
idx2  = id_result.Indices(2);
mol2  = emk.mol.fromSmiles(ref_smi_valid(idx2));
figure("Name", "2nd Rank Match");
emk.viz.draw2d(mol2, Title=sprintf("Rank 2: %s T=%.4f", ref_names_valid(idx2), T_top2));

%[text] A: T_top1 = 1.0 (Exact match in database -- same compound).
%[text]    When the score gap between 1st and 2nd is large (> 0.2), the evidence for identification is strengthened.
%[text]    The 2nd hit is usually structurally similar TCA (e.g., Nortriptyline, Desipramine).
%[text]    Most database compounds have T < 0.2, indicating the specificity of ECFP4.
%%
%[text] ## Let's Try 5: Property Comparison with Top Candidate & Ro5

TOP_IDX = id_result.Indices(1);
logInfo("Top Candidate: %s", ref_names_valid(TOP_IDX));
logInfo("  Category: %s  Prescription Drug: %d", ref_cat_valid(TOP_IDX), ref_drug_valid(TOP_IDX));

smi_top  = ref_smi_valid(TOP_IDX);
mol_top  = emk.mol.fromSmiles(smi_top);
d_top    = emk.descriptor.calculate(mol_top, ["MolWt","LogP","TPSA","NumHDonors","NumHAcceptors"]);
logInfo("MW Difference: %.2f Da (Unknown %.1f vs Candidate %.1f)", ...
    abs(desc_unk.MolWt - d_top.MolWt), desc_unk.MolWt, d_top.MolWt);

t_top = struct2table(d_top);
ro5   = emk.filter.lipinski(t_top);
logInfo("Top Candidate Ro5: Pass=%d  Violations=%d", ro5.Pass_Ro5, ro5.Violations_Ro5);

%[text] A: MW Difference ~0 (Confirmed identical compound by exact match).
%[text]    Top candidate is a prescription drug (IsDrug = 1).
%[text]    Possession of this substance without a prescription is clinically and legally significant.
%[text]    Ro5 Pass -- Drug-like oral formulation.
%%
%[text] ## Let's Try 6: Imipramine vs Nortriptyline & Caffeine

smi_nor = refDB.SMILES(refDB.Name == "NORTRIPTYLINE");
if ~isempty(smi_nor)
    mol_nor = emk.mol.fromSmiles(smi_nor(1));
    fp_nor  = emk.fingerprint.morgan(mol_nor);
    T_nor   = emk.similarity.tanimoto(fp_unk, fp_nor);
    logInfo("Imipramine vs Nortriptyline: T=%.4f", T_nor);
else
    logWarn("NORTRIPTYLINE not found in the database");
end

mol_caf = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
fp_caf  = emk.fingerprint.morgan(mol_caf);
T_caf   = emk.similarity.tanimoto(fp_unk, fp_caf);
logInfo("Imipramine vs Caffeine: T=%.4f", T_caf);

%[text] A: Imipramine vs Nortriptyline: T ~0.26 (moderate -- not T > 0.7).
%[text]    Both are TCAs, but Imipramine has a dibenz[b,f]azepine core
%[text]    (7-membered ring with nitrogen), while Nortriptyline has a dibenz[b,f]cycloheptadiene core
%[text]    (C=C in the ring, no ring nitrogen).
%[text]    The difference in ring systems sufficiently alters the ECFP4 environment, resulting in moderate similarity.
%[text]    "Same drug class" does not imply T >> 0.5 -- scaffold differences dominate ECFP4.
%[text]    T=1.0 is required for identification confirmation; T=0.26 is insufficient.
figure("Name", "Unknown Substance (Imipramine)");
emk.viz.draw2d(mol_unk, Title="Unknown Substance (Imipramine)");
if exist("mol_nor", "var")
    figure("Name", "Nortriptyline (Metabolite)");
    emk.viz.draw2d(mol_nor, Title="Nortriptyline (Metabolite)");
end
%%
%[text] ## Let's Try 7: Cross-Validation -- MACCS Keys
%[text] Run the identification pipeline with MACCS keys and check the results of ECFP4.
fp_unk_maccs = emk.fingerprint.maccs(mol_unk);
maccs_fps    = cellfun(@(smi) emk.fingerprint.maccs(emk.mol.fromSmiles(smi)), ...
    cellstr(ref_smi_valid), "UniformOutput", false);
res_maccs    = emk.similarity.rankBy(fp_unk_maccs, maccs_fps);
logInfo("MACCS Rank 1: %s T=%.4f", ref_names_valid(res_maccs.Indices(1)), res_maccs.Scores(1));

%[text] A: The same compound should rank 1st with MACCS keys as well.
%[text]    Identification matches with different fingerprint types
%[text]    are stronger forensic evidence.
%[text]
%[text] Try Fluconazole as another unknown substance:
SMILES2  = "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F";
mol2     = emk.mol.fromSmiles(SMILES2);
fp2      = emk.fingerprint.morgan(mol2);
res2     = emk.similarity.rankBy(fp2, ref_fps);
logInfo("Fluconazole Rank 1: %s T=%.4f", ref_names_valid(res2.Indices(1)), res2.Scores(1));

%[text] Novel compound (not in the database):
SMILES3 = "CN(C)c1ccc(C2=CC=Cc3ccc(CN(C)C)cc32)cc1";
mol3    = emk.mol.fromSmiles(SMILES3);
fp3     = emk.fingerprint.morgan(mol3);
res3    = emk.similarity.rankBy(fp3, ref_fps);
logInfo("Novel compound Rank 1 T=%.4f -- Outside database (T < 1.0)", res3.Scores(1));

%[text] A: Fluconazole (exists in the database): T=1.0 (perfect match).
%[text]    Top non-self hits share a triazole ring -- reflecting ECFP4's sensitivity to triazole fragments.
%[text]    Novel compound: Maximum T < 1.0 -- confirms it is not in the database.
%[text]    This is the "database gap" issue: forensic databases
%[text]    need continuous updates to keep up with novel psychoactive substances (NPS).

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
