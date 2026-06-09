%[text] # A07 Answers: Scaffold Analysis and R Group Decomposition
%[text] Reference answers for the "Try It" exercises in a07_scaffold_analysis.m.
%[text] First, execute a07_scaffold_analysis.m (at least sections 0-5) to
%[text] construct the necessary workspace variables. Then verify with this file.
%[text]
%[text] Variables expected from a07 workspace:
%[text]   scafSmi, scafNAtoms, molNames, molSmiles, validIdx, mols, rawTbl,
%[text]   scafMap, uniqueScafs, scafFreq, scafFreqSorted, uniqueScafsSorted,
%[text]   nMols, nUnique, SDI, nSingletons, familyScafs, familyFreqs,
%[text]   rgroupTbl, scores, explained, ACYCLIC_TAG
addpath(genpath("src"));
emk.setup.initPython();

mol_warmup = emk.mol.fromSmiles("C"); clear mol_warmup;
logInfo("A07 Answers: Setup complete");
%%
%[text] ## Reconstructing Prerequisites (Corresponding to Section 0-4 of a07)

DATA_FILE   = "data/list/fda_drugs.csv";
ACYCLIC_TAG = "<acyclic>";
MIN_FAMILY_SIZE = 3;

rawTbl = readtable(DATA_FILE, TextType="string");
nRaw   = height(rawTbl);

scaffoldSmiles = strings(nRaw, 1);
scaffoldNAtoms = nan(nRaw, 1);
mols           = cell(1, nRaw);
valid          = false(1, nRaw);

logInfo("%d compounds: Extracting Murcko scaffolds ...", nRaw);
for k = 1:nRaw
    try
        mol        = emk.mol.fromSmiles(rawTbl.SMILES(k));
        scaf       = emk.mol.scaffold(mol);
        nScafAtoms = double(scaf.GetNumAtoms());
        if nScafAtoms == 0
            scaffoldSmiles(k) = ACYCLIC_TAG;
        else
            scaffoldSmiles(k) = emk.mol.toSmiles(scaf);
        end
        scaffoldNAtoms(k) = nScafAtoms;
        mols{k}           = mol;
        valid(k)          = true;
    catch; end
end

validIdx   = find(valid);
nMols      = numel(validIdx);
molNames   = rawTbl.Name(validIdx);
molSmiles  = rawTbl.SMILES(validIdx);
scafSmi    = scaffoldSmiles(validIdx);
scafNAtoms = scaffoldNAtoms(validIdx);

%[text] Constructing Scaffold Frequency Map
scafMap = containers.Map("KeyType", "char", "ValueType", "any");
for k = 1:nMols
    key = char(scafSmi(k));
    if isKey(scafMap, key)
        scafMap(key) = [scafMap(key), k];
    else
        scafMap(key) = k;
    end
end

uniqueScafs  = keys(scafMap);
nUnique      = numel(uniqueScafs);
SDI          = nUnique / nMols;
scafFreq     = zeros(1, nUnique);
for i = 1:nUnique; scafFreq(i) = numel(scafMap(uniqueScafs{i})); end
[scafFreqSorted, sortOrd] = sort(scafFreq, "descend");
uniqueScafsSorted = uniqueScafs(sortOrd);
nSingletons  = sum(scafFreq == 1);

familyScafs  = uniqueScafsSorted(scafFreqSorted >= MIN_FAMILY_SIZE);
familyFreqs  = scafFreqSorted(scafFreqSorted >= MIN_FAMILY_SIZE);
nFamilies    = numel(familyScafs);

%[text] Reconstructing R Group Property Table
familyNames    = cell(nFamilies, 1);
familyNMembers = zeros(nFamilies, 1);
familyMeanLogP = zeros(nFamilies, 1);
familyStdLogP  = zeros(nFamilies, 1);
familyMeanTPSA = zeros(nFamilies, 1);
familyMeanMW   = zeros(nFamilies, 1);
familyMeanHBD  = zeros(nFamilies, 1);
familyMeanHBA  = zeros(nFamilies, 1);
for f = 1:nFamilies
    memberIdx = scafMap(familyScafs{f});
    rawRows   = rawTbl(validIdx(memberIdx), :);
    logP = double(rawRows.ALogP);
    familyNames{f}    = familyScafs{f};
    familyNMembers(f) = numel(memberIdx);
    familyMeanLogP(f) = mean(logP);
    familyStdLogP(f)  = std(logP);
    familyMeanTPSA(f) = mean(double(rawRows.TPSA));
    familyMeanMW(f)   = mean(double(rawRows.MolecularWeight));
    familyMeanHBD(f)  = mean(double(rawRows.HBondDonors));
    familyMeanHBA(f)  = mean(double(rawRows.HBondAcceptors));
end
rgroupTbl = table( ...
    string(familyNames), familyNMembers, ...
    round(familyMeanLogP, 2), round(familyStdLogP, 2), ...
    round(familyMeanTPSA, 1), round(familyMeanMW, 1), ...
    round(familyMeanHBD, 2), round(familyMeanHBA, 2), ...
    VariableNames=["Scaffold","N","mean_ALogP","std_ALogP", ...
                   "mean_TPSA","mean_MW","mean_HBD","mean_HBA"]);
rgroupTbl = sortrows(rgroupTbl, "N", "descend");

%[text] PCA (For Try It Yourself 4)
FEAT_COLS = ["MolecularWeight","ALogP","TPSA","HBondDonors", ...
             "HBondAcceptors","RotatableBonds"];
X_pca = zeros(nMols, numel(FEAT_COLS));
for fi = 1:numel(FEAT_COLS)
    X_pca(:,fi) = double(rawTbl.(FEAT_COLS(fi))(validIdx));
end
mu_pca = mean(X_pca,1); sigma_pca = std(X_pca,0,1);
sigma_pca(sigma_pca == 0) = 1;
X_std = (X_pca - mu_pca) ./ sigma_pca;
[~, scores, ~, ~, explained] = pca(X_std);

logInfo("Reconstructing prerequisites: %d molecules, %d unique scaffolds, SDI=%.3f", ...
    nMols, nUnique, SDI);
%%
%[text] ## Let's Try 1: Acyclic Count; Largest Scaffold; Aspirin vs Ibuprofen

logInfo("Let's Try 1 -- Acyclic Molecules:");
nAcyclic    = sum(scafSmi == ACYCLIC_TAG);
acyclicMask = scafSmi == ACYCLIC_TAG;
logInfo("  Acyclic (no ring system): %d / %d  (%.1f%%)", ...
    nAcyclic, nMols, 100*nAcyclic/nMols);
if any(acyclicMask)
    logInfo("  Acyclic drug names: %s", strjoin(molNames(acyclicMask), ", "));
end
%[text] Answer: Typically 5-10 drugs in the FDA set are acyclic (e.g. some amino-acid-
%[text]    derived drugs, simple carboxylic acids, GABA analogues).
%[text]    Acyclic drugs are unusual in modern drug discovery because ring systems
%[text]    provide shape, rigidity, and metabolic stability.
[maxAtoms, maxIdx] = max(scafNAtoms);
logInfo("Let's Try 1 -- Largest Scaffold:");
logInfo("  Compound: %s", molNames(maxIdx));
logInfo("  Scaffold SMILES: %s", scafSmi(maxIdx));
logInfo("  Scaffold Heavy Atom Count: %d", maxAtoms);
%[text] Answer: The compound with the largest scaffold is typically a polycyclic natural-
%[text]    product-derived drug (e.g. a steroid, macrolide antibiotic, or alkaloid).
%[text]    Scaffolds with > 20 heavy atoms represent complex ring systems that would
%[text]    be difficult to synthesise de novo in a drug discovery campaign.
logInfo("Let's Try 1 -- Aspirin vs Ibuprofen Scaffold:");
mol_asp = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
scaf_asp = emk.mol.scaffold(mol_asp);
smi_asp  = emk.mol.toSmiles(scaf_asp);

mol_ibu = emk.mol.fromSmiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1");
scaf_ibu = emk.mol.scaffold(mol_ibu);
smi_ibu  = emk.mol.toSmiles(scaf_ibu);

logInfo("  Aspirin Scaffold: %s (%d atoms)", smi_asp, double(scaf_asp.GetNumAtoms()));
logInfo("  Ibuprofen Scaffold: %s (%d atoms)", smi_ibu, double(scaf_ibu.GetNumAtoms()));
yesNo = ["No", "Yes"];
logInfo("  Same Scaffold?: %s", yesNo(strcmp(char(smi_asp), char(smi_ibu))+1));
%[text] Answer: Both aspirin and ibuprofen are built on the benzene ring scaffold
%[text]    (c1ccccc1, 6 heavy atoms).  Despite very different pharmacology and
%[text]    substituents, they share the same Murcko framework.
%[text]    This illustrates why scaffold frequency alone does not capture the
%[text]    full diversity of drug action -- R-groups matter enormously for activity.
%%
%[text] ## Let's Try 2: Scaffolds with >= 3 Drugs; Shannon Entropy; Most Frequent Scaffold

logInfo("Let's Try 2 -- Privileged Scaffolds (Frequency >= 3):");
nPrivileged = sum(scafFreq >= 3);
pctPriv     = nPrivileged / nUnique * 100;
logInfo("  Scaffolds with >= 3 drugs: %d / %d (%.1f%%)", ...
    nPrivileged, nUnique, pctPriv);
if pctPriv < 10
    logInfo("  (note: < 10%% -- below typical 10-20%% for full FDA set; expected for a %d-compound subset)", nMols);
end
%[text] Answer: Typically 10-20% of scaffolds appear in 3+ drugs.
%[text]    These "privileged scaffolds" (Evans et al.) are frameworks that recur
%[text]    because they provide: (1) metabolic stability of aromatic rings,
%[text]    (2) H-bond acceptors from heteroatoms, (3) appropriate size for
%[text]    protein binding pockets, (4) synthetic accessibility.
%[text]    Benzene, pyridine, piperidine, morpholine and indole are classic examples.
logInfo("Let's Try 2 -- Shannon Entropy of Scaffold Frequency Distribution:");
p = scafFreq / sum(scafFreq);
p_nz = p(p > 0);
H = -sum(p_nz .* log2(p_nz));
H_max = log2(nMols);   % maximum entropy: all singletons, uniform distribution
logInfo("  Shannon Entropy H = %.3f bits", H);
logInfo("  Max entropy (all singletons) = %.3f bits  (log2(%d))", H_max, nMols);
logInfo("  Relative Entropy H / H_max = %.3f", H / H_max);
if H / H_max < 0.85
    logInfo("  (note: below typical 0.85-0.95 for full FDA set; expected for a %d-compound subset with large benzene family)", nMols);
end
%[text] Answer: Shannon entropy measures the "spread" of probability across scaffolds.
%[text]    H near H_max (= log2(nMols)) means the scaffold distribution is nearly
%[text]    uniform -- many singletons, few dominant scaffolds.
%[text]    H much less than H_max means a few scaffolds dominate.
%[text]    FDA drugs typically show H/H_max ~ 0.85-0.95 (moderately diverse).
%[text]    A combinatorial library of one scaffold would give H = 0.
logInfo("Let's Try 2 -- Most Frequent Scaffold:");
topScaf = uniqueScafsSorted{1};
topFreq = scafFreqSorted(1);
topMembers = molNames(scafMap(topScaf));
logInfo("  Top Scaffold: %s  (Frequency=%d)", topScaf, topFreq);
logInfo("  Member drugs: %s", strjoin(topMembers, ", "));
%[text] Answer: The most frequent scaffold in the FDA set is often a simple benzene ring
%[text]    (c1ccccc1) or a closely related aromatic system (e.g. piperidine c1ccncc1).
%[text]    Benzene's prevalence reflects: (a) it is the smallest carbocyclic ring,
%[text]    so many substituent patterns all collapse to the same scaffold;
%[text]    (b) phenyl groups provide flat hydrophobic surface for pi-stacking;
%[text]    (c) they are metabolically stable and synthetically accessible.
%[text]    Not all benzene-scaffold drugs are related -- aspirin, ibuprofen, and
%[text]    paracetamol all share this scaffold but treat different conditions.
%%
%[text] ## Let's Try 3: Maximum ALogP Range; Highest Average TPSA; Most Frequent Family Members

logInfo("Let's Try 3 -- Scaffold family with the widest ALogP range:");
[maxStd, maxStdIdx] = max(rgroupTbl.std_ALogP);
topRow = rgroupTbl(maxStdIdx, :);
logInfo("  Scaffold: %s", topRow.Scaffold);
logInfo("  N=%d  std_ALogP=%.2f  mean_ALogP=%.2f", ...
    topRow.N, topRow.std_ALogP, topRow.mean_ALogP);
%[text] Answer: The scaffold with the widest ALogP range is often a common ring system
%[text]    (benzene, piperidine) that accommodates both highly polar and lipophilic
%[text]    substituents.  High std_ALogP means R-group changes dominate lipophilicity
%[text] -- the scaffold itself is "neutral" with respect to LogP.
%[text]    This is valuable in drug discovery: a "flexible" scaffold lets the chemist
%[text]    tune LogP without changing the core binding geometry.
logInfo("Let's Try 3 -- Scaffold family with the highest average TPSA:");
[maxTPSA, maxTPSAIdx] = max(rgroupTbl.mean_TPSA);
topTPSA = rgroupTbl(maxTPSAIdx, :);
logInfo("  Scaffold: %s", topTPSA.Scaffold);
logInfo("  N=%d  mean_TPSA=%.1f A^2", topTPSA.N, topTPSA.mean_TPSA);
if topTPSA.mean_TPSA > 140
    logInfo("  -> mean TPSA > 140 A^2: Veber rule predicts poor oral absorption");
else
    logInfo("  -> mean TPSA <= 140 A^2: Within Veber rule (predicts good oral absorption)");
end
%[text] Answer: Scaffolds with high mean TPSA often contain multiple H-bond donors/acceptors
%[text]    (e.g. guanidine groups, multiple amide bonds, glycosidic scaffolds).
%[text]    Aminoglycosides (e.g. gentamicin) and heparin-like compounds are often
%[text]    given IV because high TPSA (>140 A^2) prevents oral absorption.
%[text]    The Veber rule: TPSA <= 140 A^2 and RotatableBonds <= 10 for good oral BA.
logInfo("Let's Try 3 -- Members of the most frequent scaffold family:");
key = char(rgroupTbl.Scaffold(1));
memberIdx = scafMap(key);
logInfo("  Family scaffold: %s  (n=%d)", key, numel(memberIdx));
logInfo("  Member drugs:");
for k = 1:numel(memberIdx)
    logInfo("    %s", molNames(memberIdx(k)));
end
%[text] Answer: Members of the most frequent scaffold family are often structurally diverse
%[text]    NSAIDs, antihistamines, or CNS drugs that share a common aromatic core.
%[text]    They rarely belong to the same pharmacological class -- demonstrating that
%[text]    scaffold identity alone does not predict mechanism of action.
logInfo("Let's Try 3 -- SAR significance of singleton scaffolds:");
logInfo("  Singletons: %d / %d scaffolds (%.0f%%)", ...
    nSingletons, nUnique, 100*nSingletons/nUnique);
%[text] Answer: A singleton scaffold (one molecule) has no within-family variation.
%[text]    SAR requires at least 2-3 analogue pairs to establish a trend.
%[text]    A medicinal chemist cannot derive R-group SAR from a single data point.
%[text]    Singleton scaffolds are "orphan frameworks" -- they may represent novel
%[text]    chemical space but cannot be optimised within this dataset alone.
%%
%[text] ## Try It 4: Scaffold Clustering in PCA Space; PC Variance; Intra-family Tanimoto

logInfo("Try It 4 -- Variance captured by PCA:");
logInfo("  PC1 = %.1f%%, PC2 = %.1f%%, cumulative = %.1f%%", ...
    explained(1), explained(2), sum(explained(1:2)));
logInfo("  3D cumulative (PC1+PC2+PC3) = %.1f%%", sum(explained(1:3)));
%[text] Answer: For drug-like descriptor matrices (6 features), PC1+PC2 typically captures
%[text]    55-70% of variance.  PC3 adds another 10-15%.  If PC1+PC2 < 60%, the 2D
%[text]    projection misses significant structure -- check for outliers or skewed
%[text]    descriptors that inflate variance of a single PC.
%[text]    Common result: PC1 ~ "size" (MW, HeavyAtoms), PC2 ~ "polarity" (TPSA, HBD).
logInfo("Try It 4 -- Intra-family Tanimoto vs Global:");
topFamilyKey = uniqueScafsSorted{1};
topFamilyMemberIdx = scafMap(topFamilyKey);
nTopFamily = numel(topFamilyMemberIdx);

if nTopFamily >= 2
    mols_family = mols(validIdx(topFamilyMemberIdx));
    fps_family  = cell(1, nTopFamily);
    for k = 1:nTopFamily
        fps_family{k} = emk.fingerprint.morgan(mols_family{k});
    end

    % Intra-family pairwise Tanimoto
    simMat_family = zeros(nTopFamily, nTopFamily);
    for i = 1:nTopFamily
        for j = 1:nTopFamily
            simMat_family(i,j) = emk.similarity.tanimoto(fps_family{i}, fps_family{j});
        end
    end
    offDiag = simMat_family(~logical(eye(nTopFamily)));
    meanIntra = mean(offDiag);

    % Global mean Tanimoto (sample 50 random pairs for speed)
    rng(42);
    nSample   = min(50, nMols);
    sampleIdx = randperm(nMols, nSample);
    mols_sample = mols(validIdx(sampleIdx));
    fps_sample  = cell(1, nSample);
    for k = 1:nSample
        fps_sample{k} = emk.fingerprint.morgan(mols_sample{k});
    end
    globalSim = zeros(nSample*(nSample-1)/2, 1);
    cnt = 0;
    for i = 1:nSample
        for j = i+1:nSample
            cnt = cnt + 1;
            globalSim(cnt) = emk.similarity.tanimoto(fps_sample{i}, fps_sample{j});
        end
    end
    meanGlobal = mean(globalSim);

    logInfo("Try It 4 -- Tanimoto Similarity:");
    logInfo("  Intra-family mean (top scaffold, n=%d): %.4f", nTopFamily, meanIntra);
    logInfo("  Global random mean (n=%d molecules, %d pairs): %.4f", nSample, cnt, meanGlobal);
    delta = meanIntra - meanGlobal;
    if delta > 0
        logInfo("  Within-family higher than global? yes (delta=+%.4f)", delta);
    else
        logInfo("  Within-family higher than global? no (delta=%.4f)", delta);
    end
    if abs(delta) < 0.10
        logInfo("  (note: small delta -- scaffold '%s' (freq=%d) is small relative to full Morgan fingerprint; R-groups dominate)", ...
            topFamilyKey, scafFreqSorted(1));
    end
else
    logWarn("TRY IT 4: top family has < 2 members -- Tanimoto skipped.");
end
%[text] Answer: Within-scaffold Tanimoto is typically HIGHER than global random pairs,
%[text]    because family members share the same ring scaffold (a substantial portion
%[text]    of their fingerprint bits are identical).  The difference quantifies how
%[text]    much the scaffold contributes to overall structural similarity vs R-groups.
%[text]    If intra-family T is only slightly above global, the scaffold is small
%[text]    relative to R-groups -- R-groups dominate the fingerprint.
%%
%[text] ## Let's Try 5: SDI Comparison; Scaffold Richness Curve; Everyday Chemicals

logInfo("Let's Try 5 -- SDI Comparison:");
logInfo("  Calculated SDI = %.3f  (Langdon 2011 reports ~0.75 for FDA drugs)", SDI);
logInfo("  Mean family size = %.2f  (nMols/nUnique = %d/%d)", ...
    nMols/nUnique, nMols, nUnique);
logInfo("  Largest family = %d  (scaffold: %s)", ...
    scafFreqSorted(1), uniqueScafsSorted{1});
%[text] Answer: Langdon et al. (2011) report SDI ~ 0.75 for their FDA drug set.
%[text]    Minor differences arise from database version (more recent FDA sets include
%[text]    more biologics with complex scaffolds), RDKit version (scaffold canonicalisation
%[text]    may differ), and how acyclic molecules are handled.
%[text]    SDI > 0.7 confirms the FDA drug set is structurally diverse, not dominated
%[text]    by a single scaffold family.  This is expected: approved drugs treat
%[text]    hundreds of different targets and diseases.
logInfo("Let's Try 5 -- Scaffold Richness Curve:");
rng(2026);
nPerms = 3;   % show 3 random orderings
figure("Name","A07 Answer: Scaffold Richness Curve");
hold on;
for p = 1:nPerms
    perm     = randperm(nMols);
    richness = arrayfun(@(n) numel(unique(scafSmi(perm(1:n)))), 1:nMols);
    plot(1:nMols, richness, "-", LineWidth=1.2, ...
        DisplayName=sprintf("Order %d", p));
end
hold off;
xlabel("Number of Molecules Added"); ylabel("Unique Scaffolds");
title("Scaffold Richness Curve (3 Random Orders)");
legend(Location="southeast"); grid on;
%[text] Log a brief quantitative summary
logInfo("  Final unique scaffolds at N=%d: %d (SDI=%.3f)", nMols, nUnique, SDI);
if (nUnique / nMols) > 0.6
    logInfo("  -> Curve is nearly linear (no plateau): high scaffold diversity, library not saturated");
else
    logInfo("  -> Curve flattens before N=%d: saturation detected", nMols);
end
%[text] Answer: The scaffold richness curve typically rises steeply at first (early drugs
%[text]    introduce new scaffolds) then flattens as the library fills in.
%[text]    The "knee" of the curve approximates the coverage saturation point.
%[text]    Different random orderings produce similar curves (not strongly order-
%[text]    dependent), which confirms the statistical robustness of the SDI measure.
logInfo("Let's Try 5 -- Everyday Chemicals Comparison:");
EVERYDAY_FILE = "data/list/everyday_chemicals.csv";
try
    evTbl = readtable(EVERYDAY_FILE, TextType="string");
    nEv   = height(evTbl);
    evScafSmi = strings(nEv, 1);
    evValid   = false(1, nEv);
    for k = 1:nEv
        try
            mol = emk.mol.fromSmiles(evTbl.SMILES(k));
            scaf = emk.mol.scaffold(mol);
            if double(scaf.GetNumAtoms()) == 0
                evScafSmi(k) = ACYCLIC_TAG;
            else
                evScafSmi(k) = emk.mol.toSmiles(scaf);
            end
            evValid(k) = true;
        catch; end
    end
    nEvValid  = sum(evValid);
    nEvUnique = numel(unique(evScafSmi(evValid)));
    SDI_ev    = nEvUnique / nEvValid;
    logInfo("  Everyday Chemicals: n=%d  Unique Scaffolds=%d  SDI=%.3f", ...
        nEvValid, nEvUnique, SDI_ev);
    logInfo("  FDA Drugs SDI=%.3f  vs  Everyday Chemicals SDI=%.3f", SDI, SDI_ev);
catch ME
    logWarn("Everyday chemicals file not found or error: %s", ME.message);
end
%[text] Answer: The everyday chemicals dataset (30 compounds) typically shows a LOWER SDI
%[text]    than the FDA drug set (200 compounds), for two reasons:
%[text]    (1) Smaller absolute size: with only 30 compounds, there are fewer unique
%[text]        scaffolds possible even at maximum diversity (SDI_max = 1.0).
%[text]    (2) The everyday chemicals include many simple benzene-based molecules
%[text]        (ibuprofen, paracetamol, aspirin, caffeine-like structures) that share
%[text]        common scaffolds, making the set less diverse per compound.
%[text]    This illustrates the size-dependence of SDI: raw SDI values should always
%[text]    be compared between libraries of similar size.
logInfo("A07 Answer: Completed.");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]%   data: {"layout":"inline","rightPanelPercent":40}
%---
