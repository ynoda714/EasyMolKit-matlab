%[text] # A07: Scaffold Analysis and R Group Decomposition
%[text] EasyMolKit Analytics — Layer 3
%[text]
%[text] When lining up 200 FDA-approved drugs, you notice that many share the same ring system core.
%[text] This "scaffold" is called the **Bemis-Murcko scaffold**, and medicinal chemists use it to understand which scaffolds support multiple drugs during the design phase.
%[text] Even with the same scaffold, changing only the side chains (R groups) can alter absorption, toxicity, and activity, allowing systematic optimization of R groups. This is the essence of **Structure-Activity Relationship (SAR)** analysis.
%[text] In this script, we automatically extract scaffolds from an FDA-approved drug dataset and explore scaffold diversity using frequency ranking, R group property tables, and PCA visualization.
%[text]
%[text] **Story**
%[text]
%[text] A medicinal chemist is reviewing the structural diversity of 200 FDA-approved drugs before designing a new compound library. She posed three questions:
%[text]
%[text] - (a) What is the most common core framework (scaffold) in approved drugs, and
%[text] how many unique scaffolds are there?
%[text] - (b) Which physicochemical properties (LogP: logarithm of the partition coefficient between water and octanol, TPSA: topological polar surface area, MW: molecular weight) differ among scaffold families? This is the essence of SAR analysis, linking structural changes to property changes.
%[text] - (c) How do the same side chains (R groups) affect drug-likeness within a specific scaffold family?
%[text]
%[text] In this exercise, you will:
%[text]
%[text] 1. Load 200 FDA drugs and extract the Bemis-Murcko scaffold for each drug.
%[text] 2. Map unique scaffold SMILES to identify scaffold families.
%[text] 3. Rank scaffolds by frequency and plot the distribution.
%[text] 4. Construct an R group property table: mean and standard deviation of key descriptors per family.
%[text] 5. Visualize SAR: PCA colored by scaffold class and ALogP box plots.
%[text] 6. Summarize scaffold diversity with the Scaffold Diversity Index (SDI).
%[text]
%[text] **Learning Objectives**
%[text]
%[text] - Understand the definition of Bemis-Murcko scaffolds and their use in SAR
%[text] - Use `emk.mol.scaffold()` to extract canonical scaffolds from SMILES
%[text] - Group molecules by scaffold SMILES using MATLAB `containers.Map`
%[text] - Construct an R group property table summarizing intra-scaffold variance
%[text] - Visualize how scaffolds partition chemical space (PCA, box plots)
%[text] - Calculate and interpret the Scaffold Diversity Index (SDI)
%[text]
%[text] **Prerequisites**
%[text]
%[text] - Completion of F05 (Substructure Search) — the concept of scaffolds has been introduced
%[text] - Recommended: Familiarity with A01 (PCA) and A02 (Clustering) for visualization context
%[text] - Statistics and Machine Learning Toolbox required (`pca`, `boxplot`)
%[text] - No internet connection required
%[text]
%[text] Estimated time: 30–45 minutes
%[text]
%[text] **Data:**
%[text]
%[text] `data/list/fda_drugs.csv` — 200 FDA-approved drugs (ChEMBL CC-BY-SA 3.0)
%[text] Columns: ChEMBLID, Name, SMILES, MolecularWeight, ALogP,
%[text] HBondDonors, HBondAcceptors, TPSA, RotatableBonds
%[text]
%[text] **References**
%[text]
%[text] Bemis & Murcko (1996) *J. Med. Chem.* 39:2887–2893. doi:10.1021/jm9602928
%[text] — Original paper on the concept of scaffolds
%[text]
%[text] Bemis & Murcko (1999) *J. Med. Chem.* 42:5095–5099. doi:10.1021/jm9903996
%[text] — Extension to R group/side chain analysis
%[text]
%[text] Langdon et al. (2011) *J. Chem. Inf. Model.* 51:2174–2185. doi:10.1021/ci200319g
%[text] — Quantitative scaffold diversity analysis of FDA drug set
%[text]
%[text] RDKit `MurckoScaffold` module:
%[text] https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
%[text]
%[text] How to run: Execute each section with Ctrl+Enter

%%
%[text] ## Section 0: Setup

% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
emk.setup.initPython();

%[text] Warm up the Python/RDKit process before the main execution.
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;
logSection("A07", "Section 0: Setup", "Analytics L3");
%%
%[text] ## Section 1: Loading FDA-approved Drugs and Extracting Murcko Scaffolds
%[text]
%[text] Setup is complete. First, load the FDA-approved drug data and extract the Bemis-Murcko scaffolds from each molecule.
%[text] If extracted correctly, drugs with the same ring system core will be automatically grouped.
%[text]
%[text] ### Concept: Bemis-Murcko Scaffold
%[text] A scaffold (molecular framework) is the core ring system of a molecule after removing all terminal substituents (side chains).
%[text]
%[text] **Bemis-Murcko Definition (1996)** retains the following elements:
%[text] - All ring atoms (including heteroatoms like N, O, S)
%[text] - Linker bonds and atoms connecting rings
%[text] - Single-atom linkers between rings
%[text]
%[text] Side chains attached to the ring system are removed. For example:
%[text] - Aspirin `CC(=O)Oc1ccccc1C(=O)O` → Scaffold `c1ccccc1` (Benzene)
%[text] - Ibuprofen `CC(C)Cc1ccc(C(C)C(=O)O)cc1` → Scaffold `c1ccc(cc1)` (Benzene)
%[text]
%[text] Both have a benzene scaffold, but the R groups differ significantly.
%[text] Scaffold analysis reveals which frameworks support multiple drugs.
%[text]
%[text] Molecules without rings (e.g., amino acids, short peptides) have no scaffold.
%[text] RDKit returns an empty Mol (0 atoms) for these.
%[text] Label them as `"<acyclic>"` to retain in the dataset.
logSection("A07", "Section 1: Loading FDA-approved Drugs and Extracting Murcko Scaffolds", "Analytics L3");
DATA_FILE   = "data/list/fda_drugs.csv";
ACYCLIC_TAG = "<acyclic>";   % placeholder for ring-free molecules

rawTbl = readtable(DATA_FILE, TextType="string");
nRaw   = height(rawTbl);
logInfo("Loaded %d rows from %s", nRaw, DATA_FILE);

%[text] Extract scaffolds for each drug
scaffoldSmiles = strings(nRaw, 1);   % Scaffold SMILES for each drug
scaffoldNAtoms = nan(nRaw, 1);       % Number of heavy atoms in the scaffold
mols           = cell(1, nRaw);
valid          = false(1, nRaw);

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
    catch ME
        logWarn("Skipping row %d (%s): %s", k, rawTbl.Name(k), ME.message);
    end
end

validIdx = find(valid);
nMols    = numel(validIdx);
logInfo("Extracted scaffolds for %d / %d molecules", nMols, nRaw);

%[text] Restrict working arrays to valid molecules
molNames    = rawTbl.Name(validIdx);
molSmiles   = rawTbl.SMILES(validIdx);
scafSmi     = scaffoldSmiles(validIdx);
scafNAtoms  = scaffoldNAtoms(validIdx);

%[text] Preview (first 5 rows)
previewTbl = table( ...
    molNames(1:5), molSmiles(1:5), scafSmi(1:5), scafNAtoms(1:5), ...
    VariableNames=["Name", "SMILES", "Scaffold", "ScaffoldAtoms"]);
disp(previewTbl);

%[text] **💡 Observation Point 1 — Let's check the scaffold extraction results**
%[text] Out of 200 drugs, how many are acyclic (ring-free)?
%[text] Run `sum(scafSmi == "<acyclic>")` to check.
%[text] Do you recognize any names among the acyclic drugs? (Hint: simple amino acids, some antiviral drugs)
%[text] Which is the largest scaffold (most heavy atoms)?
%[text] Does the scaffold SMILES appear as a complex ring system?
%[text] The SMILES for Aspirin is `"CC(=O)Oc1ccccc1C(=O)O"`.
%[text] What scaffold is obtained? Next, try with Ibuprofen `"CC(C)Cc1ccc(C(C)C(=O)O)cc1"`.
%[text] Do they share the same scaffold?
% ... (Try writing code here)
%%
%[text] ## Section 2: Scaffold Frequency Analysis
%[text]
%[text] Scaffold extraction is complete. Next, we will examine the frequency distribution of how often each scaffold appears in different drugs.
%[text] Let's check how concentrated the privileged structures are.
%[text]
%[text] ### Concept: Scaffold Diversity and Scaffold Frequency
logSection("A07", "Section 2: Scaffold Frequency Analysis", "Analytics L3");
%[text] By examining the distribution of scaffold frequency, we can determine whether a drug library is "diverse" (many unique scaffolds each appearing only once) or "concentrated" (many compounds sharing a few scaffolds).
%[text]
%[text] Trends in the FDA drug set (Langdon et al. 2011):
%[text] - About 50% of scaffolds are "singletons" (appear in only one drug)
%[text] - A few scaffolds (benzene, pyridine, piperidine rings) appear in 5-10 or more drugs
%[text]
%[text] **Why is this important?**
%[text] - High scaffold diversity = indicates broad coverage of chemical space
%[text] - High-frequency scaffolds are "**privileged structures**",
%[text] evolutionarily selected frameworks that fit many protein binding sites
%[text]
%[text] **Definition of Scaffold Diversity Index (SDI):**
%[text] $\text{SDI} = N_{\text{unique}} / N_{\text{total}}$
%[text] The range of SDI is from $1/N$ (all the same scaffold) to $1.0$ (all unique).
%[text] A random drug-like library typically achieves SDI ~ 0.7-0.9.
%[text]
%[text] Map scaffold SMILES to a list of molecule indices (using containers.Map).
scafMap = containers.Map("KeyType", "char", "ValueType", "any");
for k = 1:nMols
    key = char(scafSmi(k));
    if isKey(scafMap, key)
        scafMap(key) = [scafMap(key), k];
    else
        scafMap(key) = k;
    end
end

uniqueScafs = keys(scafMap);
nUnique     = numel(uniqueScafs);
SDI         = nUnique / nMols;

logInfo("Unique scaffolds: %d / %d molecules  (SDI = %.3f)", ...
    nUnique, nMols, SDI);

%[text] Calculate the frequency of each scaffold.
scafFreq = zeros(1, nUnique);
for i = 1:nUnique
    scafFreq(i) = numel(scafMap(uniqueScafs{i}));
end

%[text] Sort frequencies in descending order.
[scafFreqSorted, sortOrd] = sort(scafFreq, "descend");
uniqueScafsSorted = uniqueScafs(sortOrd);

%[text] Report the top 10 scaffolds.
TOP_N = 10;
logInfo("Top %d scaffolds by frequency:", TOP_N);
for i = 1:min(TOP_N, nUnique)
    memberIdx  = scafMap(uniqueScafsSorted{i});
    memberList = strjoin(molNames(memberIdx), ", ");
    if strlength(memberList) > 80
        memberList = extractBefore(memberList, 81) + "...";
    end
    logInfo("  [%2d] Frequency=%d  Scaffold=%s  Members=%s", ...
        i, scafFreqSorted(i), uniqueScafsSorted{i}, memberList);
end

nSingletons = sum(scafFreq == 1);
logInfo("Singletons (frequency=1): %d / %d scaffolds (%.0f%%)", ...
    nSingletons, nUnique, 100 * nSingletons / nUnique);

%[text] -- Bar Graph: Frequency of Top 10 Scaffolds --
figure("Name", "A07 Scaffold Frequency");
barh(scafFreqSorted(1:TOP_N), FaceColor=[0.3 0.6 0.9]);
yticks(1:TOP_N);
%[text] Use the first 30 characters of scaffold SMILES as labels.
shortLabels = cell(TOP_N, 1);
for i = 1:TOP_N
    s = uniqueScafsSorted{i};
    if numel(s) > 30
        shortLabels{i} = [s(1:30) "..."];
    else
        shortLabels{i} = s;
    end
end
yticklabels(shortLabels);   % YDir="reverse" makes y=1 the top, so no flipud needed
xlabel("Number of FDA drugs with the same scaffold");
title(sprintf("Top %d Most Frequent Murcko Scaffolds (FDA Drugs)", TOP_N));
grid on;
set(gca, FontSize=8, YDir="reverse");

%[text] -- Histogram: Distribution of All Scaffold Frequencies --
figure("Name", "A07 Scaffold Frequency Distribution");
histogram(scafFreq, max(scafFreq), FaceColor=[0.3 0.6 0.9], EdgeColor="none");
xlabel("Number of drugs per scaffold");
ylabel("Number of scaffolds");
title(sprintf("Scaffold Frequency Distribution (Unique N=%d)", nUnique));
xline(mean(scafFreq), "--r", sprintf("Mean=%.1f", mean(scafFreq)), ...
    LabelHorizontalAlignment="right");
grid on;

%[text] **💡 Observation Point 2 ―― Let's Analyze the Frequency Distribution**
%[text] Check the proportion of scaffolds appearing in 3 or more drugs.
%[text] You can calculate this using `sum(scafFreq >= 3) / nUnique * 100`.
%[text] Consider why these are called "privileged" scaffolds.
%[text] Calculate the Shannon entropy of the scaffold frequency distribution.
%[text] A higher H means a more even distribution among scaffolds.
%[text] Compared to a hypothetical "all singletons" library, `H = log2(nMols)`.
%[text] In many FDA analyses, the most frequent scaffold is a simple benzene ring (`c1ccccc1`).
%[text] Check if this is the case in this dataset as well. Check `uniqueScafsSorted{1}`.
%[text] Why is benzene so common?
%[text] (Hint: Aromatic rings provide hydrophobic surfaces, pi-stacking, and metabolic stability)
% ... (Try writing code here)
%%
%[text] ## Section 3: R-Group Property Table
%[text]
%[text] We have understood the distribution of scaffolds. Next, we will quantify how R-groups (side chains) within each scaffold family affect properties.
%[text] You will see how much ALogP (an indicator of lipophilicity) and TPSA (topological polar surface area) can vary even with the same scaffold.
%[text]
%[text] ### Concept: R-Group Analysis and Structure-Activity Relationship (SAR)
%[text] In medicinal chemistry, "R-group" refers to variable substituents added to a fixed scaffold. R-group analysis poses the following question:
%[text]
%[text] With a fixed scaffold, how do different side chains affect properties such as lipophilicity (ALogP), polarity (TPSA), and molecular weight?
%[text]
%[text] This is the core of lead optimization. Drug discovery researchers modify the R-groups of promising scaffolds to improve ADMET (absorption, distribution, metabolism, excretion, toxicity) properties while maintaining activity.
%[text]
%[text] In this section, we will construct a property summary table for each scaffold family with three or more members. The report for each family includes:
%[text] - Number of members (number of approved drugs sharing the same scaffold)
%[text] - Mean and standard deviation of ALogP, TPSA, MolecularWeight, HBondDonors, HBondAcceptors, RotatableBonds
%[text]
%[text] A high `std(ALogP)` within a scaffold indicates that the scaffold allows a wide range of lipophilicity due to R-group changes. This is a "flexible" scaffold.
%[text] A low std indicates that all drugs on that scaffold have similar lipophilicity.
logSection("A07", "Section 3: R-Group Property Table", "Analytics L3");
MIN_FAMILY_SIZE = 3;   % minimum members to include in the R-group table

%[text] Select scaffold families that meet the size threshold.
%[text] Note: Exclude the `<acyclic>` placeholder. Molecules without ring systems do not have a fixed scaffold,
%[text] and thus are not subject to R-group analysis (only ring system drugs can be compared in SAR).
realScafMask = scafFreqSorted >= MIN_FAMILY_SIZE & ~strcmp(uniqueScafsSorted, char(ACYCLIC_TAG));
familyScafs  = uniqueScafsSorted(realScafMask);
familyFreqs  = scafFreqSorted(realScafMask);
nFamilies    = numel(familyScafs);
logInfo("Scaffold families with members >= %d (ring systems only): %d", MIN_FAMILY_SIZE, nFamilies);

%[text] Use property columns present in the CSV (to avoid RDKit overhead).
PROP_COLS = ["ALogP", "TPSA", "MolecularWeight", "HBondDonors", ...
             "HBondAcceptors", "RotatableBonds"];

%[text] Construct the R-group summary table.
familyNames     = cell(nFamilies, 1);
familyNMembers  = zeros(nFamilies, 1);
familyMeanLogP  = zeros(nFamilies, 1);
familyStdLogP   = zeros(nFamilies, 1);
familyMeanTPSA  = zeros(nFamilies, 1);
familyMeanMW    = zeros(nFamilies, 1);
familyMeanHBD   = zeros(nFamilies, 1);
familyMeanHBA   = zeros(nFamilies, 1);

for f = 1:nFamilies
    memberIdx = scafMap(familyScafs{f});   % indices into validIdx
    % Retrieve rows from rawTbl for these members
    rawRows = rawTbl(validIdx(memberIdx), :);

    logP = double(rawRows.ALogP);
    tpsa = double(rawRows.TPSA);
    mw   = double(rawRows.MolecularWeight);
    hbd  = double(rawRows.HBondDonors);
    hba  = double(rawRows.HBondAcceptors);

    familyNames{f}    = familyScafs{f};
    familyNMembers(f) = numel(memberIdx);
    familyMeanLogP(f) = mean(logP);
    familyStdLogP(f)  = std(logP);
    familyMeanTPSA(f) = mean(tpsa);
    familyMeanMW(f)   = mean(mw);
    familyMeanHBD(f)  = mean(hbd);
    familyMeanHBA(f)  = mean(hba);

    % Display member names for the first 5 families
    if f <= 5
        memberStr = strjoin(rawRows.Name, " | ");
        if strlength(memberStr) > 100
            memberStr = extractBefore(memberStr, 101) + "...";
        end
        logInfo("Family %d (n=%d): %s", f, familyNMembers(f), memberStr);
    end
end

%[text] Assemble the summary table.
rgroupTbl = table( ...
    string(familyNames), familyNMembers, ...
    round(familyMeanLogP, 2), round(familyStdLogP, 2), ...
    round(familyMeanTPSA, 1), round(familyMeanMW, 1), ...
    round(familyMeanHBD, 2), round(familyMeanHBA, 2), ...
    VariableNames=["Scaffold","N","mean_ALogP","std_ALogP", ...
                   "mean_TPSA","mean_MW","mean_HBD","mean_HBA"]);
rgroupTbl = sortrows(rgroupTbl, "N", "descend");

logInfo("R-Group Property Table (%d scaffold families, >= %d members):", ...
    nFamilies, MIN_FAMILY_SIZE);
disp(rgroupTbl);

%[text] **💡 Observation Point 3 ―― Let's interpret the R-Group Table**
%[text] Check the scaffold family with the widest ALogP range (highest `std_ALogP`).
%[text] Is the scaffold's SMILES a common aromatic ring?
%[text] What does high LogP variance indicate about R-group flexibility?
%[text] Check the scaffold with the highest mean TPSA and interpret if it is the most polar.
%[text] High TPSA ($> 140\,\text{Å}^2$) is associated with poor oral absorption (Veber rule).
%[text] Can you list drug classes that are usually administered intravenously?
%[text] Display all member drug names of the most frequent scaffold.
%[text] Do these drugs belong to the same pharmacological class?
%[text] Scaffolds with n=1 member (singletons) have no variation within the family.
%[text] What does this mean for SAR? Is R-group optimization possible with a single data point?
%[text] (Hint: SAR requires multiple analogs)
% ... (Try writing code here)
%%
%[text] ## Section 4: SAR Visualization — Property Variation Among Scaffold Families
%[text]
%[text] The R group property table is complete. Next, we will visually examine the structure-activity relationship (SAR) using box plots and principal component analysis (PCA).
%[text] Let's see how scaffold families partition the chemical space.
%[text]
%[text] ### Concept: Structure-Activity Relationship (SAR) Visualization
logSection("A07", "Section 4: SAR Visualization — Property Variation Among Scaffold Families", "Analytics L3");
%[text] SAR plots show how structural changes by different scaffold families affect physicochemical properties.
%[text]
%[text] Here, we use two complementary perspectives.
%[text]
%[text] **(a) ALogP Box Plot for Each Scaffold Family**
%[text] Each box shows the median and interquartile range of ALogP within a family.
%[text] A wide box (large interquartile range) indicates diverse R groups affecting lipophilicity.
%[text]
%[text] **(b) PCA of Descriptor Matrix Colored by Scaffold Class**
%[text] Project all drugs into a 2D chemical space. If scaffold families form distinct clusters, it indicates that scaffold identity predicts the overall physicochemical profile.
%[text] This is characteristic of "scaffold-driven" SAR.
%[text] Overlapping families suggest R groups dominate properties.
%[text]
%[text] PCA Coloring by Scaffold:
%[text] - Assign a unique color to each "top scaffold family."
%[text] - Plot singletons and small families in gray.
%[text] - Separation in PCA space is evidence that scaffolds drive overall properties.
%[text]
%[text] Use only drugs from scaffold families with at least MIN_FAMILY_SIZE members for the box plot.
familyMask   = ismember(scafSmi, string(familyScafs));
nFamilyMols  = sum(familyMask);
logInfo("Molecules in scaffold families (>= %d members): %d / %d", ...
    MIN_FAMILY_SIZE, nFamilyMols, nMols);

if nFamilyMols > 0
    % Construct label vector for grouped box plot
    familyScafStr = string(familyScafs);
    [~, familyLabelIdx] = ismember(scafSmi(familyMask), familyScafStr);

    % Abbreviate scaffold labels for display (12 characters)
    shortScafLabels = cell(nFamilies, 1);
    for f = 1:nFamilies
        s = familyScafs{f};
        shortScafLabels{f} = s(1:min(12, numel(s)));
    end

    % Retrieve ALogP for family molecules
    logP_family = double(rawTbl.ALogP(validIdx(familyMask)));

    % Sort families by mean ALogP for ordered display
    [~, sortByMean] = sort(familyMeanLogP, "descend");
    remappedLabels  = zeros(size(familyLabelIdx));
    for fi = 1:nFamilies
        remappedLabels(familyLabelIdx == sortByMean(fi)) = fi;
    end

    % Note: LabelOrientation="inline" causes graphics timeout in MATLAB Online
    % (renderLabels freezes with many labels). Use boxplot with numeric groups
    % and set XTickLabel separately -- faster rendering.
    figure("Name", "A07 ALogP by Scaffold Family");
    boxplot(logP_family, remappedLabels, Labels=shortScafLabels(sortByMean));
    xtickangle(30);
    xlabel("Scaffold (First 12 SMILES Characters, Sorted by Mean ALogP)");
    ylabel("ALogP");
    title(sprintf("ALogP Distribution by Scaffold Family (Family >= %d Members)", ...
        MIN_FAMILY_SIZE));
    yline(0, "--", "LogP=0", Color=[0.5 0.5 0.5]);
    yline(5, "--", "Ro5 Limit", Color=[0.8 0.2 0.2]);
    grid on;
else
    logWarn("No scaffold families with >= %d members found -- skipping box plot.", ...
        MIN_FAMILY_SIZE);
end

%[text] --- PCA of Descriptor Matrix Colored by Scaffold Family ---
%[text] Use descriptors provided by CSV to avoid per-molecule RDKit calls for speed.
FEAT_COLS = ["MolecularWeight", "ALogP", "TPSA", "HBondDonors", ...
             "HBondAcceptors", "RotatableBonds"];
nFeats    = numel(FEAT_COLS);

%[text] Construct the feature matrix for all valid molecules.
X_pca = zeros(nMols, nFeats);
for fi = 1:nFeats
    X_pca(:, fi) = double(rawTbl.(FEAT_COLS(fi))(validIdx));
end

%[text] Perform standardization (zero mean, unit variance).
mu_pca    = mean(X_pca, 1);
sigma_pca = std(X_pca, 0, 1);
sigma_pca(sigma_pca == 0) = 1;   % Guard against zero variance columns
X_std = (X_pca - mu_pca) ./ sigma_pca;

%[text] Perform PCA.
[~, scores, ~, ~, explained] = pca(X_std);

logInfo("PCA: PC1=%.1f%%, PC2=%.1f%%, Cumulative=%.1f%%", ...
    explained(1), explained(2), sum(explained(1:2)));

%[text] Assign colors based on scaffold family membership.
%[text] Rank families by frequency, assign distinct colors to the top 6, others in gray.
MAX_COLOUR_FAMILIES = 6;
palette = lines(MAX_COLOUR_FAMILIES);   % MATLAB colourmap with distinct colours
colourIdx = zeros(nMols, 1);            % 0 = grey (singleton / small family)

for f = 1:min(MAX_COLOUR_FAMILIES, nFamilies)
    memberIdx = scafMap(familyScafs{f});
    colourIdx(memberIdx) = f;
end

figure("Name", "A07 Chemical Space PCA by Scaffold");
hold on;

%[text] Plot gray "other" molecules first to create a background layer.
greyMask = colourIdx == 0;
scatter(scores(greyMask, 1), scores(greyMask, 2), 30, ...
    [0.75 0.75 0.75], "filled", MarkerFaceAlpha=0.4, DisplayName="Others");

%[text] Then plot colored scaffold families.
legendHandles = gobjects(1, min(MAX_COLOUR_FAMILIES, nFamilies) + 1);
legendHandles(1) = scatter(scores(greyMask, 1), scores(greyMask, 2), 30, ...
    [0.75 0.75 0.75], "filled", MarkerFaceAlpha=0.4);  % Others

for f = 1:min(MAX_COLOUR_FAMILIES, nFamilies)
    memberIdx = scafMap(familyScafs{f});
    col = palette(f, :);
    s = familyScafs{f};
    labelStr = sprintf("Scaf %d (n=%d): %s", f, familyFreqs(f), ...
        s(1:min(20, numel(s))));
    legendHandles(f + 1) = scatter( ...
        scores(memberIdx, 1), scores(memberIdx, 2), ...
        50, col, "filled", MarkerFaceAlpha=0.85, DisplayName=labelStr);
end

hold off;
xlabel(sprintf("PC1 (%.1f%%)", explained(1)));
ylabel(sprintf("PC2 (%.1f%%)", explained(2)));
title("Chemical Space PCA -- Colored by Scaffold Family");
legend(legendHandles, ["Others (Singleton / Small Family)", ...
    arrayfun(@(f) sprintf("Scaffold %d (n=%d)", f, familyFreqs(f)), ...
    1:min(MAX_COLOUR_FAMILIES, nFamilies), UniformOutput=false)], ...
    Location="best", FontSize=7);
grid on;

%[text] **💡 Observation Point 4 — Let's Visualize SAR**
%[text] Do molecules of the same scaffold family form clusters in PCA space?
%[text] If they do, scaffolds strongly determine the physicochemical profile.
%[text] If not, R groups are often dominant, seen in flexible scaffolds.
%[text] How much of the total variance do PC1 and PC2 explain?
%[text] If less than 60%, the 2D projection might miss important structures.
%[text] Try this: Check how much 3D PCA captures with `sum(explained(1:3))`.
%[text] Calculate pairwise Tanimoto similarity within the largest scaffold family,
%[text] and compare with random pairs across the dataset.
%[text] `mean(simMat_family(~logical(eye(numel(memberIdx)))))`
%[text] Is the intra-family Tanimoto similarity higher than the global average in A02?
% ... (Try writing code here)
%%
%[text] ## Section 5: Scaffold Diversity Summary
%[text]
%[text] Visualization is complete. Finally, we will summarize the diversity of this library numerically using the Scaffold Diversity Index (SDI).
%[text] Does the FDA drug dataset cover a broad chemical space?
%[text]
%[text] ### Concept: Scaffold Diversity Index and Library Design
logSection("A07", "Section 5: Scaffold Diversity Summary", "Analytics L3");
%[text] The Scaffold Diversity Index (SDI) is a single number that indicates how much a compound library "spreads" in chemical space:
%[text]
%[text] **SDI = Number of Unique Scaffolds / Total Number of Molecules**
%[text]
%[text] **Interpretation**:
%[text] - If SDI is close to 1.0 → Each molecule has a unique scaffold (maximum diversity)
%[text] - If SDI is close to 0.0 → All molecules share one scaffold (focused library)
%[text]
%[text] **Related Metrics**:
%[text] - Singleton Ratio: Proportion of scaffolds appearing only once
%[text] (High ratio = diverse set, difficult to establish structure-activity relationship (SAR) without analogs)
%[text] - Average Family Size: Total molecules / Unique scaffolds (inverse of SDI)
%[text] - Maximum Family Size: Dominance of a single privileged scaffold
%[text]
%[text] **Library Design Heuristics**:
%[text] - Diverse Screening Library: SDI > 0.7, Singleton Ratio > 50%
%[text] - Focused SAR Library: Target SDI < 0.5, families of 5-20 analogs
%[text] - FDA Approved Drugs: SDI typically 0.7-0.8 (Langdon 2011)
logInfo("=== Scaffold Diversity Summary ===");
logInfo("Total Number of Analyzed Molecules: %d", nMols);
logInfo("Unique Scaffolds: %d", nUnique);
logInfo("Scaffold Diversity Index: %.3f", SDI);
logInfo("Singleton Scaffolds: %d (%.0f%%)", ...
    nSingletons, 100 * nSingletons / nUnique);
logInfo("Average Family Size: %.2f molecules / scaffold", nMols / nUnique);
logInfo("Maximum Scaffold Family: %d molecules (%s)", ...
    scafFreqSorted(1), uniqueScafsSorted{1});
logInfo("Acyclic Molecules: %d (no ring systems)", ...
    sum(scafSmi == ACYCLIC_TAG));

%[text] Scaffold Size Distribution (Number of Heavy Atoms in Scaffolds)
scafAtomsNonAcyclic = scafNAtoms(scafNAtoms > 0);
logInfo("Scaffold Heavy Atom Count: Mean=%.1f, Std Dev=%.1f, Range=[%d,%d]", ...
    mean(scafAtomsNonAcyclic), std(scafAtomsNonAcyclic), ...
    min(scafAtomsNonAcyclic), max(scafAtomsNonAcyclic));

figure("Name", "A07 Scaffold Size Distribution");
histogram(scafAtomsNonAcyclic, "BinWidth", 2, ...
    FaceColor=[0.3 0.6 0.9], EdgeColor="none");
xlabel("Number of Heavy Atoms in Scaffold");
ylabel("Number of Drugs");
title("Distribution of Scaffold Sizes (Murcko Scaffolds)");
xline(mean(scafAtomsNonAcyclic), "--r", ...
    sprintf("Mean=%.1f", mean(scafAtomsNonAcyclic)), ...
    LabelHorizontalAlignment="right");
grid on;

%[text] **Summary**
%[text]
%[text] - Bemis-Murcko scaffolds extract ring system cores and automatically group drug families
%[text] - In FDA approved drugs, a few "privileged scaffolds" are shared by many drugs
%[text] - R group analysis clarifies structure-activity relationships (SAR) within the same scaffold
%[text] - The Scaffold Diversity Index (SDI) allows quantitative evaluation of library chemical space coverage
%[text]

%[text] **💡 Observation Point 5 — Let's Evaluate Scaffold Diversity**
%[text] How does the SDI of the FDA drug set compare to the value expected in a random diversity library? (Langdon 2011 reports SDI ~ 0.75 for 836 FDA approved drugs)
%[text] Hint: Since this exercise uses a subset of 200 drugs, the SDI may be slightly lower.
%[text] What can be said about the relationship between subset size and SDI?
%[text] Let's calculate the "Scaffold Richness" curve:
%[text] Examine how many unique scaffolds appear as you add the first 10, 20, ... 200 drugs.
%[text] `figure; plot(1:nMols, richness, "-");`
%[text] `xlabel("Number of Molecules Added"); ylabel("Unique Scaffolds"); title("Scaffold Richness Curve");`
%[text] How quickly does the curve flatten?
%[text] Let's redo the scaffold analysis using the `everyday_chemicals.csv` dataset.
%[text] Is the SDI of a smaller, curated set high or low?
%[text] What does it tell us about how dataset size affects diversity metrics?
% ... (Try writing code here)

%[appendix]{"version":"1.0"}
%[metadata:view]

%---

%   data: {"layout":"inline","rightPanelPercent":40}
%---
