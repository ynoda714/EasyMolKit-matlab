%[text] # S01: Find Caffeine Cousins
%[text] EasyMolKit Application Story ─ Layer 2
%[text] 
%[text] Why do coffee, tea, and chocolate provide that unique "wakefulness"? 
%[text] The answer is not just caffeine. Plants produce a family of wakefulness-inducing substances that are chemically similar to caffeine. In this script, you will experience finding "caffeine cousins" from a database of 30 everyday chemicals using molecular fingerprints and similarity scores.
%[text] ## Learning Objectives
%[text] - Encode molecular structures as bit vectors with `emk.fingerprint.morgan`
%[text] - Perform bulk searches of chemical databases with `emk.similarity.rankBy`
%[text] - Interpret Tanimoto scores in real-world contexts
%[text] - Understand the relationship between structural similarity and bioactivity \
%[text] ## Prerequisites
%[text] - Completion of F01 (Molecule Drawing) and F04 (Similarity)
%[text] - RDKit installed (execute `emk.setup.install()` once)
%[text] - No additional Toolbox required (runs with MATLAB only) \
%[text] **Duration**: 15-20 minutes | Execution: Run each section with Ctrl+Enter
%[text] **Data**
%[text] - `data/list/everyday_chemicals.csv` — 30 general molecules (PubChem CC0) \
%%
%[text] ## Section 0: Setup
%[text] Initialize paths and Python environment.
%[text] **Always execute this section first.**
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
%[text] Warm up the Python/RDKit process (the first call may take some time).
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;
logInfo("S01: Setup complete");
%%
%[text] ## Section 1: Caffeine
%[text] First, let's understand the chemical properties of **Caffeine**.
%[text] Caffeine (1,3,7-trimethylxanthine) is a plant alkaloid belonging to the methylxanthines. The SMILES encodes a fused bicyclic ring (purine skeleton) with three methyl groups ($&dollar;CH\_3&dollar;$) and two carbonyl oxygens ($&dollar;=O&dollar;$). Here, we will learn how to load the chemical structure and calculate the physicochemical properties of the molecule (molecular weight, lipophilicity, polar surface area, etc.).
%[text] ### Caffeine Content Reference
%[text] - Espresso 80 mg/shot
%[text] - Tea 50 mg/cup
%[text] - Dark Chocolate 20 mg/30 g
%[text] - Decaf Coffee \< 5 mg/cup \
CAFFEINE_SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C";
CAFFEINE_NAME   = "Caffeine";

mol_caffeine = emk.mol.fromSmiles(CAFFEINE_SMILES);
logInfo("Analyzed caffeine. Number of heavy atoms: %d", double(mol_caffeine.GetNumHeavyAtoms()));
desc = emk.descriptor.calculate(mol_caffeine, ["MolWt", "LogP", "TPSA", "RingCount"]);
logInfo("Properties of caffeine:");
logInfo("  Molecular Weight   : %.2f g/mol  (C8H10N4O2)", desc.MolWt);
logInfo("  LogP               : %.2f        (slightly hydrophilic)", desc.LogP);
logInfo("  TPSA               : %.1f A^2   (CNS permeability < 90 is a guideline)", desc.TPSA);
logInfo("  Ring Count         : %d          (purine bicyclic skeleton)", desc.RingCount);
%[text] Draw the structure of caffeine
figure("Name", "Caffeine", "Position", [100 100 440 380]);
emk.viz.draw2d(mol_caffeine, Title="Caffeine (1,3,7-trimethylxanthine)", ...
    Width=350, Height=350);
%[text] **✏️ Try It 1 — Check BBB Permeability**
%[text] Caffeine is a CNS stimulant. The blood-brain barrier (BBB) is more permeable when $&dollar&;TPSA \< 90\\text{ \\AA}^2&dollar&;$ and $&dollar&;LogP&dollar&;$ is not too high (usually $&dollar&;\< 5&dollar&;$). Check the output above to see if caffeine meets these two conditions.
%[text] **Expected Results and Tips**
%[text] - **Expected Values**: $&dollar&;TPSA \\sim 62\\text{ \\AA}^2&dollar&;$ (good as it is \< 90), $&dollar&;LogP \\sim -1.0&dollar&;$ (very high water solubility due to negative value).
%[text] - **Note**: A negative $&dollar&;LogP&dollar&;$ indicates that caffeine dissolves easily in water (coffee or tea). Despite the low $&dollar&;LogP&dollar&;$, caffeine efficiently permeates the BBB due to its small size and planar aromatic structure. \
% ... (Try writing code here)
%%
%[text] ## Section 2: Load the Everyday Chemicals Database
%[text] The CSV contains 30 molecules found in food, pharmaceuticals, and daily products.
%[text] Each row: CommonName, CID, SMILES, MolecularFormula, MolecularWeight, Category, Source
%[text] Load a small database (CSV file) registered with "30 familiar chemicals" to be explored. This table includes components of foods (sugars and amino acids), pharmaceuticals, and daily products that we commonly encounter, each linked with structural information (SMILES). Let's check the distribution of the categories of molecules included.
%[text] - `CommonName`: Common name
%[text] - `CID`: PubChem Compound ID
%[text] - `SMILES`: String describing the structure
%[text] - `Category`: Classification such as "Pharmaceutical" or "Food Additive" \
dataFile = fullfile(projectRoot, "data", "list", "everyday_chemicals.csv");
tbl = readtable(dataFile, "TextType", "string");

logInfo("Database loaded: %d molecules / %d categories", ...
    height(tbl), numel(unique(tbl.Category)));
%[text] Preview of categories
cats = unique(tbl.Category);
logInfo("List of categories:");
for k = 1:numel(cats)
    n = sum(tbl.Category == cats(k));
    logInfo("  %-15s: %d molecules", cats(k), n);
end
%[text] **✏️ Try It 2 — Browse the Database**
%[text] Let's view the entire table in the MATLAB workspace.
%[text] Which category has the most molecules?
%[text] Which molecule has the highest molecular weight?
%[text] Hint: max(tbl.MolecularWeight) and tbl(tbl.MolecularWeight == max(...), :)
% ... (Try writing code here)
%%
%[text] ## Section 3: Create Fingerprints for the Entire Database
%[text] Computers cannot effectively calculate "similarity" using SMILES strings directly. Therefore, it is necessary to convert structures into a mathematically computable representation.
%[text] **Morgan fingerprints** (ECFP4) check each atom in a molecule and its "surrounding connections (radius of 2 atoms)" to convert molecular features into 2048 "0 or 1 flags (bits)." It is like a **digital fingerprint of chemical structures**. Molecules with similar functional groups or ring structures (scaffolds) will have similar fingerprint patterns (positions where 1s appear).
%[text] In this section, we will use a loop to generate fingerprints for all molecules in the previously loaded database.
logInfo("Calculating Morgan fingerprints for %d molecules...", height(tbl));
fps   = cell(1, height(tbl));
valid = true(1, height(tbl));     % Track molecules successfully analyzed

for i = 1:height(tbl)
    smi = tbl.SMILES(i);
    if ~emk.mol.isValid(smi)
        logWarn("  Skipping: %s (invalid SMILES)", tbl.CommonName(i));
        valid(i) = false;
        continue;
    end
    mol    = emk.mol.fromSmiles(smi);
    fps{i} = emk.fingerprint.morgan(mol);
end

nValid = sum(valid);
logInfo("Fingerprint calculation complete: %d / %d molecules", nValid, height(tbl));
%[text] Create a filtered list with valid entries
fps_valid    = fps(valid);
names_valid  = tbl.CommonName(valid);
smiles_valid = tbl.SMILES(valid);
%[text] **✏️ Try It 3 — Count the ON Bits**
%[text] How many bits are ON in the fingerprint of caffeine?
%[text] Hint:
%[text]   fp\_caf = `emk.fingerprint.morgan(mol_caffeine)`;
%[text]   bits   = `emk.fingerprint.toArray(fp_caf)`;
%[text]   sum(bits)
%[text] Expected: Approximately 30 to 50 bits ON out of 2048.
% ... (Try writing the code here)
%%
%[text] ## Section 4: Rank All Molecules by Similarity to Caffeine
%[text] emk.similarity.rankBy utilizes the fast bulk calculation feature (BulkTanimotoSimilarity) of Python/RDKit running in the background.
%[text] Calculating by sending data back and forth between MATLAB and Python one by one incurs communication overhead (time loss) proportional to the number of molecules ($&dollar&;N&dollar&;$). However, this function efficiently processes by sending all $&dollar&;N&dollar&;$ data to Python at once and receiving the results in one go, allowing the entire database search to be completed astonishingly quickly.
%[text] It compares the fingerprint (bit vector) of caffeine with the fingerprints of all molecules in the database and calculates similarity using the Tanimoto Coefficient. The Tanimoto Coefficient is a metric ranging from 0.0 (not similar at all) to 1.0 (identical).
fp_caffeine = emk.fingerprint.morgan(mol_caffeine);

result = emk.similarity.rankBy(fp_caffeine, fps_valid);   % All molecules

logInfo("--- Top 10 Molecules Most Similar to Caffeine ---");
logInfo("%-5s  %-20s  %s", "Rank", "Molecule Name", "Tanimoto");
for k = 1:min(10, numel(result.Scores))
    idx   = result.Indices(k);
    name  = names_valid(idx);
    if strlength(name) == 0; name = "(No Name)"; end
    score = result.Scores(k);
    % Skip caffeine itself (score = 1.0)
    if score >= 1.0
        logInfo("  %2d.  %-20s  %.4f  <-- Caffeine itself", k, name, score);
    else
        logInfo("  %2d.  %-20s  %.4f", k, name, score);
    end
end
%[text] **✏️ Try It 4 — Check the Top Hits**
%[text] Which is the top rank excluding caffeine itself?
%[text] Is that molecule included in the database CSV?
%[text] Expected: Theobromine (found in cocoa/chocolate) should rank high.
%[text] Theobromine has the same xanthine skeleton but only two N-methyl groups.
% ... (Try writing code here)
%%
%[text] ## Section 5: Understanding Top Hits
%[text] Create a summary table of the top 5 most similar molecules (excluding self-match).
%[text] Identify the index of Caffeine in the rank list and extract the "top 5 hits (excluding itself)" from the search results. Format it into MATLAB's `table` for detailed observation to facilitate data analysis. Pay attention to how structural similarity (Tanimoto score) correlates with biological classification (category) and actual use.
%[text] First, identify and exclude the index of Caffeine itself from the rank list.
selfIdx = find(result.Scores >= 1.0, 1);
%[text] Create the top 5 excluding self-match
topK = 5;
topRows = struct("Name", {}, "SMILES", {}, "Tanimoto", {}, "Category", {});
count = 0;
for k = 1:numel(result.Scores)
    if k == selfIdx
        continue;
    end
    count = count + 1;
    idx = result.Indices(k);
    topRows(count).Name     = names_valid(idx);
    topRows(count).SMILES   = smiles_valid(idx);
    topRows(count).Tanimoto = result.Scores(k);
    % Search for category from the original table
    match = tbl.CommonName == names_valid(idx);
    topRows(count).Category = tbl.Category(find(match, 1));
    if count >= topK; break; end
end
%[text] Display as a MATLAB table
topTbl = struct2table(topRows);
disp("Top 5 Caffeine-related compounds:");
disp(topTbl);
%[text] **Chemical Insights**
%[text] Methylxanthines (Caffeine, Theobromine, Theophylline) share a purine bicyclic scaffold.
%[text] Tanimoto \> 0.5 reflects this shared scaffold.
%[text] Scores for non-xanthine molecules are distinctly lower (\< 0.2).
%%
%[text] ## Section 6: Draw the Top 3 Analogues
%[text] Visualize the most similar molecules side by side.
%[text] Let's visually confirm the "actual shape differences" behind the numerical (Tanimoto score). Draw the top 3 molecular structures with the highest scores in 2D and compare them with the structure of caffeine (drawn in Section 1). Which functional groups are common, and where are the changes?
topDraw = min(3, numel(topRows));
figure("Name", "Top Caffeine Analogues", "Position", [100 100 1200 400]);
for k = 1:topDraw
    subplot(1, topDraw, k);
    mol_hit = emk.mol.fromSmiles(topRows(k).SMILES);
    titleStr = sprintf("%s\n(Tanimoto = %.3f)", topRows(k).Name, topRows(k).Tanimoto);
    emk.viz.draw2d(mol_hit, Title=titleStr);
end
logInfo("Displayed the structures of the top %d in one figure.", topDraw);
%[text] **✏️ Try It 5 — Investigate Structural Differences**
%[text] Compare the structure of caffeine from Section 1 with theobromine above.
%[text] What are the structural differences between them? Count the number of N-methyl groups in each molecule.
%[text] Expected: Caffeine has 3 N-CH3 groups (positions 1, 3, and 7 on the xanthine skeleton).
%[text] Theobromine has 2 N-CH3 groups (positions 3 and 7); N1 is NH.
%[text] This single difference lowers the Tanimoto score below 1.0.
% ... (Try writing code here)
%%
%[text] ## Section 7: Similarity Heatmap of Caffeine and Top 4 Relatives
%[text] Compare caffeine and its top 4 nearest neighbors using a heatmap.
%[text] So far, we have compared "Caffeine vs Other Molecules" in a one-to-many manner. But how similar are the "relatives" that were hit? We will create a pairwise similarity matrix (Distance Matrix) that cross-compares the top 5 molecules including caffeine, and visualize it as a heatmap. This allows us to visually identify "chemical clusters" formed by similar molecules.
nHeat   = min(4, numel(topRows));
heatSmi = [CAFFEINE_SMILES, arrayfun(@(r) r.SMILES, topRows(1:nHeat), ...
           "UniformOutput", false)];
heatNames = [CAFFEINE_NAME, arrayfun(@(r) char(r.Name), topRows(1:nHeat), ...
             "UniformOutput", false)];
%[text] Analysis and Fingerprint Calculation
heatFps = cell(1, numel(heatSmi));
for k = 1:numel(heatSmi)
    if iscell(heatSmi)
        smi = heatSmi{k};
    else
        smi = heatSmi(k);
    end
    heatFps{k} = emk.fingerprint.morgan(emk.mol.fromSmiles(smi));
end

S = emk.similarity.matrix(heatFps);

figure("Name", "Caffeine Relatives Heatmap", "Position", [100 100 520 460]);
imagesc(S);
colormap("hot");
colorbar;
clim([0 1]);
axis square;
xticks(1:numel(heatNames));  xticklabels(heatNames);  xtickangle(30);
yticks(1:numel(heatNames));  yticklabels(heatNames);
title("Tanimoto Similarity -- Caffeine and Nearest Relatives");

for r = 1:numel(heatNames)
    for c = 1:numel(heatNames)
        clr = "black";
        if S(r, c) < 0.5
            clr = "white";
        end
        text(c, r, sprintf("%.2f", S(r, c)), ...
            "HorizontalAlignment", "center", "FontSize", 9, "Color", clr);
    end
end
%[text] Note the brightness of the diagonal (self-similarity = 1.0) and the high score (0.53) between caffeine/theobromine. The xanthine cluster (caffeine/theobromine/...) shows high similarity with each other.
%%
%[text] ## Exercise
%[text] E1: Rank using Dice similarity (instead of Tanimoto).
%[text]     Does the order change? Why or why not?
%[text]     Hint: Use `emk.similarity.rankBy(..., Metric="dice")`
result_dice = emk.similarity.rankBy(fp_caffeine, fps_valid, Inf, Metric="dice");
% Compare result_dice.Scores(1:5) with result.Scores(1:5)
%[text] 
%[text] E2: Repeat the search using MACCS key fingerprints (instead of Morgan).
%[text]     Is Theobromine still ranked first?
%[text]     Hint: Replace `emk.fingerprint.morgan` with `emk.fingerprint.maccs`.
fp_caf_maccs = emk.fingerprint.maccs(mol_caffeine);
fps_maccs = cell(1, numel(fps_valid));
for i = 1:numel(fps_valid)
    mol_i = emk.mol.fromSmiles(smiles_valid(i));
    fps_maccs{i} = emk.fingerprint.maccs(mol_i);
end
result_maccs = emk.similarity.rankBy(fp_caf_maccs, fps_maccs);
%[text] 
%[text] E3: Find the molecule least similar to caffeine in the database.
%[text]     What is the Tanimoto score? What is its category?
%[text]     Hint: result.Scores is sorted in descending order, so the last entry is the minimum.
lastIdx   = result.Indices(end);
lastName  = names_valid(lastIdx);
lastScore = result.Scores(end);
logInfo("Least similar: %s (Tanimoto = %.4f)", lastName, lastScore);
% Check its category in tbl

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
