%[text] # A07 解答: スキャフォールド分析と R 基分解
%[text] a07_scaffold_analysis.m の「やってみよう」演習の参照解答。
%[text] 最初に a07_scaffold_analysis.m（最低セクション 0〜5）を実行して
%[text] 必要なワークスペース変数を構築してください。その後このファイルで確認。
%[text]
%[text] Variables expected from a07 workspace:
%[text]   scafSmi, scafNAtoms, molNames, molSmiles, validIdx, mols, rawTbl,
%[text]   scafMap, uniqueScafs, scafFreq, scafFreqSorted, uniqueScafsSorted,
%[text]   nMols, nUnique, SDI, nSingletons, familyScafs, familyFreqs,
%[text]   rgroupTbl, scores, explained, ACYCLIC_TAG
addpath(genpath("src"));
emk.setup.initPython();

mol_warmup = emk.mol.fromSmiles("C"); clear mol_warmup;
logInfo("A07 解答: セットアップ完了");
%%
%[text] ## 前提条件の再構築（a07 セクション 0〜4 に対応）

DATA_FILE   = "data/list/fda_drugs.csv";
ACYCLIC_TAG = "<acyclic>";
MIN_FAMILY_SIZE = 3;

rawTbl = readtable(DATA_FILE, TextType="string");
nRaw   = height(rawTbl);

scaffoldSmiles = strings(nRaw, 1);
scaffoldNAtoms = nan(nRaw, 1);
mols           = cell(1, nRaw);
valid          = false(1, nRaw);

logInfo("%d 化合物の Murcko スキャフォールドを抽出中 ...", nRaw);
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

%[text] スキャフォールド頻度マップを構築
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

%[text] R 基特性テーブルを再構築
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

%[text] PCA（やってみよう 4 用）
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

logInfo("前提条件再構築: %d 分子、%d ユニークスキャフォールド、SDI=%.3f", ...
    nMols, nUnique, SDI);
%%
%[text] ## やってみよう 1: 非環式数; 最大スキャフォールド; アスピリン vs イブプロフェン

logInfo("やってみよう 1 -- 非環式分子:");
nAcyclic    = sum(scafSmi == ACYCLIC_TAG);
acyclicMask = scafSmi == ACYCLIC_TAG;
logInfo("  非環式（環系なし）: %d / %d  (%.1f%%)", ...
    nAcyclic, nMols, 100*nAcyclic/nMols);
if any(acyclicMask)
    logInfo("  Acyclic drug names: %s", strjoin(molNames(acyclicMask), ", "));
end
%[text] 解答: Typically 5-10 drugs in the FDA set are acyclic (e.g. some amino-acid-
%[text]    derived drugs, simple carboxylic acids, GABA analogues).
%[text]    Acyclic drugs are unusual in modern drug discovery because ring systems
%[text]    provide shape, rigidity, and metabolic stability.
[maxAtoms, maxIdx] = max(scafNAtoms);
logInfo("やってみよう 1 -- 最大スキャフォールド:");
logInfo("  化合物: %s", molNames(maxIdx));
logInfo("  スキャフォールド SMILES: %s", scafSmi(maxIdx));
logInfo("  スキャフォールド重原子数: %d", maxAtoms);
%[text] 解答: The compound with the largest scaffold is typically a polycyclic natural-
%[text]    product-derived drug (e.g. a steroid, macrolide antibiotic, or alkaloid).
%[text]    Scaffolds with > 20 heavy atoms represent complex ring systems that would
%[text]    be difficult to synthesise de novo in a drug discovery campaign.
logInfo("やってみよう 1 -- アスピリン vs イブプロフェンスキャフォールド:");
mol_asp = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
scaf_asp = emk.mol.scaffold(mol_asp);
smi_asp  = emk.mol.toSmiles(scaf_asp);

mol_ibu = emk.mol.fromSmiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1");
scaf_ibu = emk.mol.scaffold(mol_ibu);
smi_ibu  = emk.mol.toSmiles(scaf_ibu);

logInfo("  アスピリン  スキャフォールド: %s (%d atoms)", smi_asp, double(scaf_asp.GetNumAtoms()));
logInfo("  イブプロフェン スキャフォールド: %s (%d atoms)", smi_ibu, double(scaf_ibu.GetNumAtoms()));
yesNo = ["いいえ", "はい"];
logInfo("  同じスキャフォールド？: %s", yesNo(strcmp(char(smi_asp), char(smi_ibu))+1));
%[text] 解答: Both aspirin and ibuprofen are built on the benzene ring scaffold
%[text]    (c1ccccc1, 6 heavy atoms).  Despite very different pharmacology and
%[text]    substituents, they share the same Murcko framework.
%[text]    This illustrates why scaffold frequency alone does not capture the
%[text]    full diversity of drug action -- R-groups matter enormously for activity.
%%
%[text] ## やってみよう 2: >= 3 薬のスキャフォールド; Shannon エントロピー; 最頻スキャフォールド

logInfo("やってみよう 2 -- 特権的スキャフォールド (頻度 >= 3):");
nPrivileged = sum(scafFreq >= 3);
pctPriv     = nPrivileged / nUnique * 100;
logInfo("  Scaffolds with >= 3 drugs: %d / %d (%.1f%%)", ...
    nPrivileged, nUnique, pctPriv);
if pctPriv < 10
    logInfo("  (note: < 10%% -- below typical 10-20%% for full FDA set; expected for a %d-compound subset)", nMols);
end
%[text] 解答: Typically 10-20% of scaffolds appear in 3+ drugs.
%[text]    These "privileged scaffolds" (Evans et al.) are frameworks that recur
%[text]    because they provide: (1) metabolic stability of aromatic rings,
%[text]    (2) H-bond acceptors from heteroatoms, (3) appropriate size for
%[text]    protein binding pockets, (4) synthetic accessibility.
%[text]    Benzene, pyridine, piperidine, morpholine and indole are classic examples.
logInfo("やってみよう 2 -- スキャフォールド頻度分布の Shannon エントロピー:");
p = scafFreq / sum(scafFreq);
p_nz = p(p > 0);
H = -sum(p_nz .* log2(p_nz));
H_max = log2(nMols);   % maximum entropy: all singletons, uniform distribution
logInfo("  Shannon エントロピー H = %.3f ビット", H);
logInfo("  Max entropy (all singletons) = %.3f bits  (log2(%d))", H_max, nMols);
logInfo("  相対エントロピー H / H_max = %.3f", H / H_max);
if H / H_max < 0.85
    logInfo("  (note: below typical 0.85-0.95 for full FDA set; expected for a %d-compound subset with large benzene family)", nMols);
end
%[text] 解答: Shannon entropy measures the "spread" of probability across scaffolds.
%[text]    H near H_max (= log2(nMols)) means the scaffold distribution is nearly
%[text]    uniform -- many singletons, few dominant scaffolds.
%[text]    H much less than H_max means a few scaffolds dominate.
%[text]    FDA drugs typically show H/H_max ~ 0.85-0.95 (moderately diverse).
%[text]    A combinatorial library of one scaffold would give H = 0.
logInfo("やってみよう 2 -- 最頻スキャフォールド:");
topScaf = uniqueScafsSorted{1};
topFreq = scafFreqSorted(1);
topMembers = molNames(scafMap(topScaf));
logInfo("  上位スキャフォールド: %s  (頻度=%d)", topScaf, topFreq);
logInfo("  Member drugs: %s", strjoin(topMembers, ", "));
%[text] 解答: The most frequent scaffold in the FDA set is often a simple benzene ring
%[text]    (c1ccccc1) or a closely related aromatic system (e.g. piperidine c1ccncc1).
%[text]    Benzene's prevalence reflects: (a) it is the smallest carbocyclic ring,
%[text]    so many substituent patterns all collapse to the same scaffold;
%[text]    (b) phenyl groups provide flat hydrophobic surface for pi-stacking;
%[text]    (c) they are metabolically stable and synthetically accessible.
%[text]    Not all benzene-scaffold drugs are related -- aspirin, ibuprofen, and
%[text]    paracetamol all share this scaffold but treat different conditions.
%%
%[text] ## やってみよう 3: 最大 ALogP 範囲; 最高平均 TPSA; 最頻ファミリーメンバー

logInfo("やってみよう 3 -- ALogP 範囲が最大のスキャフォールドファミリー:");
[maxStd, maxStdIdx] = max(rgroupTbl.std_ALogP);
topRow = rgroupTbl(maxStdIdx, :);
logInfo("  Scaffold: %s", topRow.Scaffold);
logInfo("  N=%d  std_ALogP=%.2f  mean_ALogP=%.2f", ...
    topRow.N, topRow.std_ALogP, topRow.mean_ALogP);
%[text] 解答: The scaffold with the widest ALogP range is often a common ring system
%[text]    (benzene, piperidine) that accommodates both highly polar and lipophilic
%[text]    substituents.  High std_ALogP means R-group changes dominate lipophilicity
%[text] -- the scaffold itself is "neutral" with respect to LogP.
%[text]    This is valuable in drug discovery: a "flexible" scaffold lets the chemist
%[text]    tune LogP without changing the core binding geometry.
logInfo("やってみよう 3 -- 平均 TPSA が最高のスキャフォールドファミリー:");
[maxTPSA, maxTPSAIdx] = max(rgroupTbl.mean_TPSA);
topTPSA = rgroupTbl(maxTPSAIdx, :);
logInfo("  Scaffold: %s", topTPSA.Scaffold);
logInfo("  N=%d  mean_TPSA=%.1f A^2", topTPSA.N, topTPSA.mean_TPSA);
if topTPSA.mean_TPSA > 140
    logInfo("  -> 平均 TPSA > 140 A^2: Veber ルールは経口吸収が困難と予測");
else
    logInfo("  -> 平均 TPSA <= 140 A^2: Veber ルール範囲内（経口吸収良好と予測）");
end
%[text] 解答: Scaffolds with high mean TPSA often contain multiple H-bond donors/acceptors
%[text]    (e.g. guanidine groups, multiple amide bonds, glycosidic scaffolds).
%[text]    Aminoglycosides (e.g. gentamicin) and heparin-like compounds are often
%[text]    given IV because high TPSA (>140 A^2) prevents oral absorption.
%[text]    The Veber rule: TPSA <= 140 A^2 and RotatableBonds <= 10 for good oral BA.
logInfo("やってみよう 3 -- 最頻スキャフォールドファミリーのメンバー:");
key = char(rgroupTbl.Scaffold(1));
memberIdx = scafMap(key);
logInfo("  Family scaffold: %s  (n=%d)", key, numel(memberIdx));
logInfo("  Member drugs:");
for k = 1:numel(memberIdx)
    logInfo("    %s", molNames(memberIdx(k)));
end
%[text] 解答: Members of the most frequent scaffold family are often structurally diverse
%[text]    NSAIDs, antihistamines, or CNS drugs that share a common aromatic core.
%[text]    They rarely belong to the same pharmacological class -- demonstrating that
%[text]    scaffold identity alone does not predict mechanism of action.
logInfo("やってみよう 3 -- シングルトンスキャフォールドの SAR 意義:");
logInfo("  Singletons: %d / %d scaffolds (%.0f%%)", ...
    nSingletons, nUnique, 100*nSingletons/nUnique);
%[text] 解答: A singleton scaffold (one molecule) has no within-family variation.
%[text]    SAR requires at least 2-3 analogue pairs to establish a trend.
%[text]    A medicinal chemist cannot derive R-group SAR from a single data point.
%[text]    Singleton scaffolds are "orphan frameworks" -- they may represent novel
%[text]    chemical space but cannot be optimised within this dataset alone.
%%
%[text] ## やってみよう 4: PCA 空間でのスキャフォールドクラスタリング; PC 分散; ファミリー内 Tanimoto

logInfo("やってみよう 4 -- PCA で捕捉された分散:");
logInfo("  PC1 = %.1f%%, PC2 = %.1f%%, cumulative = %.1f%%", ...
    explained(1), explained(2), sum(explained(1:2)));
logInfo("  3D cumulative (PC1+PC2+PC3) = %.1f%%", sum(explained(1:3)));
%[text] 解答: For drug-like descriptor matrices (6 features), PC1+PC2 typically captures
%[text]    55-70% of variance.  PC3 adds another 10-15%.  If PC1+PC2 < 60%, the 2D
%[text]    projection misses significant structure -- check for outliers or skewed
%[text]    descriptors that inflate variance of a single PC.
%[text]    Common result: PC1 ~ "size" (MW, HeavyAtoms), PC2 ~ "polarity" (TPSA, HBD).
logInfo("やってみよう 4 -- ファミリー内 Tanimoto vs グローバル:");
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

    logInfo("やってみよう 4 -- Tanimoto 類似度:");
    logInfo("  ファミリー内平均（上位スキャフォールド、n=%d）: %.4f", nTopFamily, meanIntra);
    logInfo("  グローバルランダム平均（n=%d 分子、%d ペア）: %.4f", nSample, cnt, meanGlobal);
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
%[text] 解答: Within-scaffold Tanimoto is typically HIGHER than global random pairs,
%[text]    because family members share the same ring scaffold (a substantial portion
%[text]    of their fingerprint bits are identical).  The difference quantifies how
%[text]    much the scaffold contributes to overall structural similarity vs R-groups.
%[text]    If intra-family T is only slightly above global, the scaffold is small
%[text]    relative to R-groups -- R-groups dominate the fingerprint.
%%
%[text] ## やってみよう 5: SDI 比較; スキャフォールド豊富度曲線; 日用化学品

logInfo("やってみよう 5 -- SDI 比較:");
logInfo("  計算 SDI = %.3f  (Langdon 2011 は FDA 薬で ~0.75 と報告)", SDI);
logInfo("  Mean family size = %.2f  (nMols/nUnique = %d/%d)", ...
    nMols/nUnique, nMols, nUnique);
logInfo("  Largest family = %d  (scaffold: %s)", ...
    scafFreqSorted(1), uniqueScafsSorted{1});
%[text] 解答: Langdon et al. (2011) report SDI ~ 0.75 for their FDA drug set.
%[text]    Minor differences arise from database version (more recent FDA sets include
%[text]    more biologics with complex scaffolds), RDKit version (scaffold canonicalisation
%[text]    may differ), and how acyclic molecules are handled.
%[text]    SDI > 0.7 confirms the FDA drug set is structurally diverse, not dominated
%[text]    by a single scaffold family.  This is expected: approved drugs treat
%[text]    hundreds of different targets and diseases.
logInfo("やってみよう 5 -- スキャフォールド豊富度曲線:");
rng(2026);
nPerms = 3;   % show 3 random orderings
figure("Name","A07 解答: スキャフォールド豊富度曲線");
hold on;
for p = 1:nPerms
    perm     = randperm(nMols);
    richness = arrayfun(@(n) numel(unique(scafSmi(perm(1:n)))), 1:nMols);
    plot(1:nMols, richness, "-", LineWidth=1.2, ...
        DisplayName=sprintf("順序 %d", p));
end
hold off;
xlabel("追加した分子数"); ylabel("ユニークスキャフォールド");
title("スキャフォールド豊富度曲線（3 ランダム順序）");
legend(Location="southeast"); grid on;
%[text] Log a brief quantitative summary
logInfo("  N=%d での最終ユニークスキャフォールド数: %d (SDI=%.3f)", nMols, nUnique, SDI);
if (nUnique / nMols) > 0.6
    logInfo("  -> Curve is nearly linear (no plateau): high scaffold diversity, library not saturated");
else
    logInfo("  -> Curve flattens before N=%d: saturation detected", nMols);
end
%[text] 解答: The scaffold richness curve typically rises steeply at first (early drugs
%[text]    introduce new scaffolds) then flattens as the library fills in.
%[text]    The "knee" of the curve approximates the coverage saturation point.
%[text]    Different random orderings produce similar curves (not strongly order-
%[text]    dependent), which confirms the statistical robustness of the SDI measure.
logInfo("やってみよう 5 -- 日用化学品比較:");
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
    logInfo("  日用化学品: n=%d  ユニークスキャフォールド=%d  SDI=%.3f", ...
        nEvValid, nEvUnique, SDI_ev);
    logInfo("  FDA 薬 SDI=%.3f  vs  日用化学品 SDI=%.3f", SDI, SDI_ev);
catch ME
    logWarn("Everyday chemicals file not found or error: %s", ME.message);
end
%[text] 解答: The everyday chemicals dataset (30 compounds) typically shows a LOWER SDI
%[text]    than the FDA drug set (200 compounds), for two reasons:
%[text]    (1) Smaller absolute size: with only 30 compounds, there are fewer unique
%[text]        scaffolds possible even at maximum diversity (SDI_max = 1.0).
%[text]    (2) The everyday chemicals include many simple benzene-based molecules
%[text]        (ibuprofen, paracetamol, aspirin, caffeine-like structures) that share
%[text]        common scaffolds, making the set less diverse per compound.
%[text]    This illustrates the size-dependence of SDI: raw SDI values should always
%[text]    be compared between libraries of similar size.
logInfo("A07 解答: 完了。");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
