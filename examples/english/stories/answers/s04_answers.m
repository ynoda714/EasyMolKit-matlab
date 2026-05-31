%[text] # S04 Answers: Introduction to Virtual Screening
%[text] Reference answers for the "Try It Yourself" exercises in s04_virtual_screening.m.
%[text] First, run s04_virtual_screening.m, then use this file to check your answers.
addpath(genpath("src"));
emk.setup.initPython();
logInfo("S04 Answers: Setup Complete");

%[text] ---- Query: Ibuprofen -----------------------------------------------
QUERY_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O";
mol_hit      = emk.mol.fromSmiles(QUERY_SMILES);
fp_hit       = emk.fingerprint.morgan(mol_hit);

%[text] ---- Library: FDA Approved Drugs CSV -------------------------------------------
library      = readtable("data/list/fda_drugs.csv", "TextType", "string");
library.MolecularWeight = double(library.MolecularWeight);
library.ALogP           = double(library.ALogP);

valid   = arrayfun(@(s) emk.mol.isValid(s), library.SMILES);
lib_names_valid  = library.Name(valid);
lib_smiles_valid = library.SMILES(valid);
lib_mw_valid     = library.MolecularWeight(valid);
lib_alogp_valid  = library.ALogP(valid);

logInfo("Library: %d / %d entries are valid", sum(valid), height(library));

lib_fps_valid = cell(1, sum(valid));
for i = 1:sum(valid)
    lib_fps_valid{i} = emk.fingerprint.morgan(emk.mol.fromSmiles(lib_smiles_valid(i)));
end
vs_result = emk.similarity.rankBy(fp_hit, lib_fps_valid);
%[text] ---
%%
%[text] ## Let's Try 1: Does Ibuprofen satisfy Ro5?

desc_hit = emk.descriptor.calculate(mol_hit, ["MolWt","LogP","NumHDonors","NumHAcceptors"]);
t_hit = struct2table(desc_hit);
ro5   = emk.filter.lipinski(t_hit);
logInfo("Ibuprofen Ro5: Pass=%d  Violations=%d", ro5.Pass_Ro5, ro5.Violations_Ro5);

%[text] Q: Why is the calculated LogP (neutral form) different from the apparent LogP in blood?
%[text] A: LogP is measured for neutral (non-ionized) molecular species.
%[text]    At blood pH 7.4, ~99% of Ibuprofen (pKa ~4.4) is ionized as carboxylate anion.
%[text]    The distribution coefficient D = LogP - log(1 + 10^(pH-pKa)) shows D << LogP at physiological pH.
%[text]    Ionized forms have low membrane permeability, but the small fraction of neutral molecules drives passive diffusion,
%[text]    so Ibuprofen is still absorbed.
%%
%[text] ## Let's Try 2: Overview of the Library

logInfo("Unique drug names in the library: %d", numel(unique(library.Name)));

[maxMW, idx] = max(library.MolecularWeight);
logInfo("Heaviest compound: %s  (MW = %.1f Da)", library.Name(idx), maxMW);

%[text] Ibuprofen ALogP vs LogP (if present in the library)
ibIdx = find(strcmpi(lib_names_valid, "IBUPROFEN"), 1);
if ~isempty(ibIdx)
    logInfo("Ibuprofen ALogP (library): %.2f  vs  Wildman-Crippen LogP: %.2f", ...
        lib_alogp_valid(ibIdx), desc_hit.LogP);
end

%[text] A: ALogP (Ghose-Crippen) and Wildman-Crippen LogP estimate the same property but
%[text]    use different atomic contribution models.
%[text]    In typical drug-like molecules, the difference between them is within ~0.5 units.
%%
%[text] ## Let's Try 3: Number of ON Bits for Ibuprofen vs Caffeine

mol_caf = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
fp_caf  = emk.fingerprint.morgan(mol_caf);
fp_ibu  = emk.fingerprint.morgan(mol_hit);

nOn_caf = sum(emk.fingerprint.toArray(fp_caf));
nOn_ibu = sum(emk.fingerprint.toArray(fp_ibu));

logInfo("Caffeine ON bits: %d", nOn_caf);
logInfo("Ibuprofen ON bits: %d", nOn_ibu);

%[text] A: Caffeine typically has more ON bits.
%[text]    The fused bicyclic purine ring system generates many different Morgan environments,
%[text]    resulting in more unique bit positions being set.
%[text]    Ibuprofen's monocyclic benzene ring and short alkyl chain have fewer types of environments.
%[text]
%[text] Experiment with radius 3:
fp_r3   = emk.fingerprint.morgan(mol_hit, Radius=3);
nOn_r3  = sum(emk.fingerprint.toArray(fp_r3));
logInfo("Ibuprofen ECFP6 (radius=3) ON bits: %d  (vs ECFP4: %d)", nOn_r3, nOn_ibu);

%[text] A: With radius 3, a larger neighborhood is considered.
%[text]    There is a possibility of fewer individual bit collisions (larger environments are more unique),
%[text]    but the overall density is usually similar or slightly reduced.
%%
%[text] ## Let's Try 4: Top Hits -- Are Other NSAIDs Included?

REPORT_TOP = 15;
for k = 1:min(REPORT_TOP, numel(vs_result.Scores))
    idx = vs_result.Indices(k);
    logInfo("  %2d. %-30s  T=%.4f", k, lib_names_valid(idx), vs_result.Scores(k));
end

%[text] A: Expected values for top 2 (excluding self): FLURBIPROFEN (T~0.40), KETOPROFEN (T~0.39).
%[text]    Both share an arylpropionic acid core (aryl--CH(CH3)--COOH).
%[text]    The fact that these close analogs remain at T ~0.40 is because
%[text]    the fluorine or additional rings in flurbiprofen/ketoprofen
%[text]    cause sufficient changes in the local ECFP4 environment.
%[text]
%[text] Use Aspirin as an alternative query:
mol_asp = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
fp_asp  = emk.fingerprint.morgan(mol_asp);
res_asp = emk.similarity.rankBy(fp_asp, lib_fps_valid);
logInfo("Top 3 hits for Aspirin: %s, %s, %s", ...
    lib_names_valid(res_asp.Indices(1)), ...
    lib_names_valid(res_asp.Indices(2)), ...
    lib_names_valid(res_asp.Indices(3)));
%[text] A: Top hits for Aspirin differ from Ibuprofen --
%[text]    ECFP4 captures the acetoxybenzoic acid scaffold.
%%
%[text] ## Let's Try 5: Hit Report -- MW Ascending Sort & Highest ALogP

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

logInfo("Hits sorted by MW in ascending order:");
sorted_tbl = sortrows(hits_tbl, "MW_Da");
disp(sorted_tbl(:, 1:5));

[~, iMaxLogP] = max(hits_tbl.ALogP);
logInfo("Highest ALogP hit: %s  ALogP=%.2f  (Drug-like < 5)", ...
    hits_tbl.Name(iMaxLogP), hits_tbl.ALogP(iMaxLogP));

%[text] A: Top hits tend to cluster around Ibuprofen MW (206 Da).
%[text]    ECFP4 similarity correlates with structural similarity, and in this scaffold class,
%[text]    it also correlates with MW.
%%
%[text] ## Let's Try 6: Score Distribution and Coverage

all_scores = vs_result.Scores;
nHigh   = sum(all_scores >= 0.65);
nModerate = sum(all_scores >= 0.40 & all_scores < 0.65);
logInfo("T >= 0.65 (Similar): %d / %d", nHigh, numel(all_scores));
logInfo("T >= 0.40 (Moderate): %d / %d", nModerate, numel(all_scores));

figure("Name", "Score Distribution -- Ibuprofen vs FDA", "Position", [100 500 560 400]);
histogram(all_scores, 30, "BinEdges", linspace(0, 1, 31));
xlabel("Tanimoto (ECFP4)");  ylabel("Count");
title("Ibuprofen vs 200 FDA Approved Drugs");
xline(0.65, "--r", "Similar (0.65)");

%[text] Comparison with Caffeine
mol_caf2 = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
fp_caf2  = emk.fingerprint.morgan(mol_caf2);
res_caf  = emk.similarity.rankBy(fp_caf2, lib_fps_valid);
figure("Name", "Score Distribution -- Caffeine vs FDA", "Position", [680 500 560 400]);
histogram(res_caf.Scores, 30, "BinEdges", linspace(0, 1, 31));
xlabel("Tanimoto (ECFP4)");  ylabel("Count");  title("Caffeine vs 200 FDA Approved Drugs");

%[text] A: Ibuprofen has a high maximum Tanimoto because it is well represented by NSAID/profen types in the FDA library.
%[text]    Caffeine, as a methylxanthine stimulant, has few structural analogs among FDA approved drugs,
%[text]    resulting in a low maximum Tanimoto.
%%
%[text] ## Let's Try 7: Draw and visually compare the top 3 hits

for k = 1:3
    mol_k = emk.mol.fromSmiles(smiVec(k));
    emk.viz.draw2d(mol_k, Title=sprintf("Hit %d: %s T=%.2f", k, nameVec(k), scoreVec(k)));
end

%[text] A: The top 3 hits (profens such as Flurbiprofen and Ketoprofen)
%[text]    visually share an allyl propionic acid core (benzene ring -- CH(CH3) -- COOH).
%[text]    The differences are only in the substituents on the aryl ring (F, second ring C=O).
%%
%[text] ## Let's Try 8: ECFP4 vs MACCS Key Comparison

fp_hit_maccs = emk.fingerprint.maccs(mol_hit);
lib_maccs    = cell(1, sum(valid));
for i = 1:sum(valid)
    lib_maccs{i} = emk.fingerprint.maccs(emk.mol.fromSmiles(lib_smiles_valid(i)));
end
res_maccs = emk.similarity.rankBy(fp_hit_maccs, lib_maccs);

logInfo("Top 3 ECFP4: %s, %s, %s", ...
    lib_names_valid(vs_result.Indices(1)), ...
    lib_names_valid(vs_result.Indices(2)), ...
    lib_names_valid(vs_result.Indices(3)));
logInfo("Top 3 MACCS: %s, %s, %s", ...
    lib_names_valid(res_maccs.Indices(1)), ...
    lib_names_valid(res_maccs.Indices(2)), ...
    lib_names_valid(res_maccs.Indices(3)));

%[text] A: ECFP4 and MACCS keys generate similar but not identical rankings.
%[text]    MACCS keys (a fixed 166-bit substructure checklist) focus on the presence of functional groups,
%[text]    while ECFP4 captures finer atomic environments.
%[text]    When both agree, the reliability of the ranking increases.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
