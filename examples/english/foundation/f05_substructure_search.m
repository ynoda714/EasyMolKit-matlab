%[text] # F05: Substructure Search and SMARTS Patterns
%[text] EasyMolKit Basic Tutorial — Layer 1
%[text] 
%[text] Wouldn't it be convenient if you could specify the structural features of a molecule directly as a search condition (query) when you want to "extract all drugs that contain carboxylic acids" or "screen only molecules that contain a benzene ring"? The language that elegantly achieves this with a single line of string is called "SMARTS." While SMILES accurately represents a specific molecule, SMARTS can describe a 'common feature or class' of molecules, such as "a group of molecules that contain carboxylic acids." Because it encodes chemists' intuition — such as benzene rings and amide bonds — directly as search queries, it is routinely used in drug discovery as a powerful filter for bulk screening of compounds with specific functional groups or for eliminating problematic structures (assay-interfering structures) that cause issues in experiments. In this script, let's explore the basics of substructure searching using SMARTS and the workflow of "PAINS filtering," which is essential in practice.
%[text] **Learning Objectives**
%[text] - Understand SMARTS (an extended pattern language for molecular queries)
%[text] - Confirm if a molecule contains a substructure using `emk.mol.hasSubstruct`
%[text] - Apply SMARTS filters to molecular datasets
%[text] - Recognize commonly used functional group patterns \
%[text] **Prerequisites**
%[text] - RDKit installed (run `emk.setup.install()` once)
%[text] - F01 to F04
%[text] - No additional Toolbox required (works with MATLAB only) \
%[text] Estimated time: 10–15 minutes | How to run: Execute sections one by one with Ctrl+Enter
%%
%[text] ## Section 0: Setup
%[text] Initialize the path and Python environment.
%[text] **Always run this section first.**
% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
emk.setup.initPython(); %[output:444a3b00]
%[text] Warm up Python/RDKit process
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;
logInfo("F05: Setup complete"); %[output:186a4c5b]
%%
%[text] ## Section 1: What is SMARTS?
%[text] ### **Referring to "one molecule" or "common features"**
%[text] At first glance, both SMILES and SMARTS appear to be strings of the same alphabet, which can confuse beginners, but their roles are fundamentally different. To illustrate, SMILES is an 'address' that points to a specific house (molecule), such as "123 Main Street, Springfield." In contrast, SMARTS serves as a 'search condition' that captures all properties that match the criteria, such as "a house with a red roof, two stories, and a garden."
%[text] For example, writing `"C"` in SMILES represents the specific molecule known as "methane," while writing `[#6]` in SMARTS means that it matches any molecule containing at least one carbon atom (atomic number 6), whether it be methane, aspirin, or DNA. This powerful pattern matching mechanism allows us to find "specific frameworks" like a puzzle within complex structures.
%[text] **Main differences from SMILES:**
%[text] - Wildcards are usable: `[#6]` = any carbon, `[*]` = any atom
%[text] - Bond types can be specified: `-` single bond, `=` double bond, `:` aromatic bond
%[text] - Logical operators: `[#6,#7]` = carbon OR nitrogen \
%[text] **Commonly used SMARTS patterns:**
%[text] - `c1ccccc1` — Aromatic benzene ring
%[text] - `[OH]` — Hydroxyl group (-OH)
%[text] - `[NH2]` — Primary amine
%[text] - `C(=O)[OH]` — Carboxylic acid
%[text] - `C(=O)[O]` — Ester or acid
%[text] - `C(=O)[N]` — Amide
%[text] - `[F,Cl,Br,I]` — Any halogen \
%%
%[text] ## Section 2: Basic Substructure Check
%[text] ### **Check if "that framework" is hidden within the molecule**
%[text] The core function for performing substructure searches in EasyMolKit is `emk.mol.hasSubstruct`. Its usage is very intuitive; you simply pass the molecule (or an array of multiple molecules) you want to investigate and the SMARTS string that serves as the condition to this function. If it matches the condition, it returns `true` (1), and if not, it returns `false` (0). First, let's have the computer determine whether the familiar aspirin contains a "benzene ring (SMARTS: `c1ccccc1`)."
%[text] Using `emk.mol.hasSubstruct(mol, query)`, you can check if the molecule contains the query substructure with `true`/`false`. The query can be specified as a SMARTS string or a SMILES string.
mol_aspirin     = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
mol_benzene     = emk.mol.fromSmiles("c1ccccc1");
mol_ethanol     = emk.mol.fromSmiles("CCO");
mol_caffeine    = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
mol_acetaminophen = emk.mol.fromSmiles("CC(=O)NC1=CC=C(C=C1)O");
%[text] Does aspirin contain a benzene ring?
has_ring = emk.mol.hasSubstruct(mol_aspirin, "c1ccccc1");
logInfo("Aspirin contains a benzene ring: %d", has_ring); %[output:42b5f0c9]
%[text] Does ethanol contain a benzene ring?
no_ring = emk.mol.hasSubstruct(mol_ethanol, "c1ccccc1");
logInfo("Ethanol contains a benzene ring: %d", no_ring); %[output:65242d3e]
%[text] Does aspirin have a carboxylic acid?
has_cooh = emk.mol.hasSubstruct(mol_aspirin, "C(=O)[OH]");
logInfo("Aspirin contains -COOH: %d", has_cooh); %[output:8b36612b]
%[text] Does acetaminophen have an amide?
has_amide = emk.mol.hasSubstruct(mol_acetaminophen, "C(=O)[NH]");
logInfo("Acetaminophen contains an amide: %d", has_amide); %[output:50c7289a]
%[text] 
%[text] **✏️ Try It 1 — Test SMARTS with caffeine**
%[text] Let's check if caffeine (`"CN1C=NC2=C1C(=O)N(C(=O)N2C)C"`) contains the following substructures.
%[text] - (a) Any nitrogen atom: `"[#7]"`
%[text] - (b) Ketone group: `"[C](=O)[#6]"`
%[text] - (c) Hydroxyl group: `"[OH]"` \
%[text] (a) will be `true`, while (b) and (c) will be `false`. Let's think about why (b) is `false`.
%[text] Hint: The C=O group in caffeine is an amide/imide carbonyl, bonded to nitrogen (N) rather than carbon. Try `"[C](=O)[N]"`.
% ... (let's write code here)
%%
%[text] ## Section 3: Batch Substructure Screening
%[text] ### **Extract target molecules from accumulated data**
%[text] The great thing about `emk.mol.hasSubstruct` is that it can take not just a single molecule but an entire "cell array" containing multiple molecules (batch processing). When you pass a list of molecules to the function, it returns the results of whether each molecule meets the condition as a logical array in MATLAB in one go.
%[text] Using this, you can automatically filter out only those compounds that have hydroxyl groups (`[OH]`) from "the 10 compounds at hand" without writing a loop (for loop). By using the returned logical array directly as MATLAB indices (logical indices), you can write clean code to extract only the molecules that match the condition.
%[text] By passing a cell array of Mol objects to `emk.mol.hasSubstruct`, you can scan the entire dataset in a single call. The result is returned as a logical row vector indicating where the substructure exists as `true`.
mols  = {mol_aspirin, mol_benzene, mol_ethanol, mol_caffeine, mol_acetaminophen};
names = {"Aspirin", "Benzene", "Ethanol", "Caffeine", "Acetaminophen"};

% Note: "c1ccccc1" matches benzene (6-membered all-carbon aromatic ring) only.
% Heterocyclic aromatics (pyridine, imidazole...) are NOT detected by this SMARTS.
% Caffeine has imidazole/pyrimidine rings so its benzene ring result is false -- expected.
patterns = struct( ...
    "description", {"Benzene ring", "Carboxylic acid", "Hydroxyl (-OH)", ...
                    "Amide (-CONH)", "Nitrogen atom"}, ...
    "smarts",      {"c1ccccc1", "C(=O)[OH]", "[OH]", "C(=O)[N]", "[#7]"} ...
);

hit_mat   = false(numel(patterns), numel(mols));
for i = 1:numel(patterns)
    hit_mat(i, :) = emk.mol.hasSubstruct(mols, patterns(i).smarts);
end
col_names = ["Aspirin", "Benzene", "Ethanol", "Caffeine", "Acetaminophen"];
row_names = string({patterns.description})';
fg_tbl    = array2table(hit_mat, "VariableNames", col_names, "RowNames", row_names);
logInfo("Functional group screening results:"); %[output:24f38584]
disp(fg_tbl); %[output:08a7c914]
%%
%[text] ## Section 4: PAINS Filtering
%[text] ### **Detecting Problematic Compounds with PAINS SMARTS Filters**
%[text] Throughout the history of drug discovery, many compounds appeared highly active in cell-based assays only to turn out as false positives — not because they were truly effective, but because they interfered with assay signals or measurement instruments. These problematic substructures are referred to as **PAINS (Pan-Assay Interference Compounds)**.
%[text] In actual drug discovery projects, to avoid wasting experimental resources, it has become standard practice to apply a "PAINS SMARTS filter" early in the workflow, discarding these compounds before they reach costly in vitro testing. The following code identifies whether the compounds at hand are likely to cause assay interference using SMARTS patterns that capture representative problematic scaffolds such as "quinone frameworks" and "azo dye structures."
%[text] The true value of the SMARTS language lies not just in stringing together element symbols but in its advanced expressiveness that allows for precise specification of the "state" of atoms. For example, `[OH]` strictly refers to a "hydroxyl group" where hydrogen is directly attached to oxygen, and writing `C(=O)[OH]` can target the clean structure of a "carboxylic acid."
%[text] Furthermore, using a comma `,` within square brackets means "or (OR)." For example, writing `[F,Cl,Br,I]` means a query for "any of fluorine, chlorine, bromine, or iodine (i.e., any halogen)." Additionally, by using a sharp symbol and atomic number like `#7` instead of element symbols, you can specify "all nitrogen atoms" regardless of whether they are aromatic (lowercase `n`) or aliphatic (uppercase `N`). By combining these symbols, you can pinpoint even the most complex blends of functional groups.
%[text] The accompanying `pains.csv` contains SMARTS patterns from Baell & Holloway 2010.
% NOTE: NumHeaderLines=1, VariableNamesLine=1 must be specified explicitly because
% SMARTS patterns contain commas inside quoted fields which confuse MATLAB's
% auto-detection (detectImportOptions infers 13 columns instead of 4).
pains_data = readtable("data/list/pains.csv", TextType="string", ...
    Delimiter=",", NumHeaderLines=1, VariableNamesLine=1);
logInfo("Loaded %d PAINS SMARTS patterns", height(pains_data)); %[output:48b4b33e]
pains_smarts = pains_data.SMARTS;   % SMARTS column by name (robust vs column-index)

%[text] Test molecules
test_smiles = ["CCO",                             % Ethanol -- clean
               "CC(=O)Oc1ccccc1C(=O)O",          % Aspirin -- clean
               "O=C1C=C(O)C(=CC1=O)c1ccccc1",    % Quinone (known PAINS)
               "c1ccc(cc1)N=Nc2ccccc2"];          % Azo dye (known PAINS)
test_names  = ["Ethanol", "Aspirin", "Quinone", "Azo dye"];

%[text] **Note (Performance)**: This loop makes `numel(test_smiles) × n_patterns` individual IPC calls. The total of 480 patterns × 4 molecules = about 60 seconds (due to Python IPC overhead). For large-scale PAINS filtering, consider using cell array batch calls or SMARTS OR combinations for speedup.
logInfo("--- PAINS Screening (Total %d Patterns) ---", height(pains_data)); %[output:65443a13]
n_patterns = height(pains_data);   % Use all 480 patterns (quinone matches row 205, azo matches row 469)

for i = 1:numel(test_smiles) %[output:group:84d44a38]
    if ~emk.mol.isValid(test_smiles(i)); continue; end
    mol  = emk.mol.fromSmiles(test_smiles(i));
    hits = 0;
    for p = 1:n_patterns
        if emk.mol.hasSubstruct(mol, pains_smarts(p))
            hits = hits + 1;
        end
    end
    status = "Clean";
    if hits > 0; status = sprintf("PAINS (%d hits)", hits); end
    logInfo("  %-14s  %s", test_names(i), status); %[output:2f7f2809]
end %[output:group:84d44a38]
%%
%[text] ## Section 5: Build Custom Functional Group Table
%[text] ### **Building a Custom Functional Group Detection Table**
%[text] The SMARTS-based approach demonstrated above can be applied to any molecular dataset. Here, we scan the `everyday_chemicals` dataset against five common functional group patterns to produce a presence/absence table. This structural summary is useful for exploratory analysis and for quickly characterizing the chemical diversity of a compound collection.
data = readtable("data/list/everyday_chemicals.csv", TextType="string");
fg_patterns = {"c1ccccc1", "[OH]", "C(=O)[OH]", "[NH,NH2]", "[C](=O)[N]"};
fg_names    = ["Benzene ring", "Hydroxyl", "Carboxylic acid", "Amine", "Amide"];

n = height(data);
fg_matrix = false(n, numel(fg_patterns));

for i = 1:n    if ~emk.mol.isValid(data.SMILES(i)); continue; end
    mol = emk.mol.fromSmiles(data.SMILES(i));
    for p = 1:numel(fg_patterns)
        fg_matrix(i, p) = emk.mol.hasSubstruct(mol, fg_patterns{p});
    end
end

fg_tbl = array2table(fg_matrix, "VariableNames", fg_names);
fg_tbl.Name = data.CommonName;
fg_tbl = movevars(fg_tbl, "Name", "Before", "Benzene ring");
disp(fg_tbl); %[output:18f19cfe]
%[text] **✏️ Try It 2 — Add a Halogen Column**
%[text] Use the SMARTS `"[F,Cl,Br,I]"` to add a `"Halogen"` column to `fg_tbl`.
%[text] Which everyday chemicals contain halogen atoms?
%[text] Hint: Add a new entry to `fg_patterns` and `fg_names`, then rerun Section 5 from the top.
% ... (Let's write code here)
%%
%[text] ## Section 6: Summary
%[text] A summary of the substructure search API learned in this tutorial and SMARTS quick reference for common practical use.
%[text] Key API learned in this section:
%[text] - `emk.mol.hasSubstruct(mol, "c1ccccc1")` — Single Mol, SMARTS query
%[text] - `emk.mol.hasSubstruct(mols, "[OH]")` — Cell array, batch mode \
%[text] **SMARTS Quick Reference:**
%[text]     Benzene ring    : `c1ccccc1`  (6-membered all-carbon aromatic ring only; heterocyclic aromatics not detected)
%[text]     Hydroxyl        : `[OH]`
%[text]     Carboxylic acid : `C(=O)[OH]`
%[text]     Amine (1/2)     : `[NH,NH2]`
%[text]     Amide           : `C(=O)[N]`
%[text]     Any halogen     : `[F,Cl,Br,I]`
%[text]     Any nitrogen    : `[#7]`
%[text]     Any carbon      : `[#6]`
%[text] 
%[text] **In next F06**, we will learn how to read and write molecular files (SDF and SMILES lists).
%[text] You will be able to import and export multiple molecules at once.
%%
%[text] ## Exercises
%[text] Try to solve each exercise before referring to `answers/f05_answers.m`.
%[text] 
%[text] **E1.** Write a script that takes a SMILES string and displays a human-readable summary of the functional groups it contains.
%[text] Target functional groups: aromatic ring, hydroxyl, carboxylic acid, amine, amide,
%[text] halogen, ketone (`"C(=O)[#6]"`).
%[text] Test with Nicotine `"CN1CCC[C@H]1C2=CN=CC=C2"` and Morphine
%[text] `"OC1=CC2=C(CC3N(CC23)CC4=CC=CC=C14)C=C1"`.
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:444a3b00]
%   data: {"dataType":"text","outputData":{"text":"[09:36:34][INFO]  initPython: Python 3.10 already active (Status: Loaded) -- skipping.\n","truncated":false}}
%---
%[output:186a4c5b]
%   data: {"dataType":"text","outputData":{"text":"[09:36:34][INFO]  F05: Setup complete\n","truncated":false}}
%---
%[output:42b5f0c9]
%   data: {"dataType":"text","outputData":{"text":"[09:36:34][INFO]  Aspirin contains a benzene ring: 1\n","truncated":false}}
%---
%[output:65242d3e]
%   data: {"dataType":"text","outputData":{"text":"[09:36:34][INFO]  Ethanol contains a benzene ring: 0\n","truncated":false}}
%---
%[output:8b36612b]
%   data: {"dataType":"text","outputData":{"text":"[09:36:34][INFO]  Aspirin contains -COOH: 1\n","truncated":false}}
%---
%[output:50c7289a]
%   data: {"dataType":"text","outputData":{"text":"[09:36:34][INFO]  Acetaminophen contains an amide: 1\n","truncated":false}}
%---
%[output:24f38584]
%   data: {"dataType":"text","outputData":{"text":"[09:36:34][INFO]  Functional group screening results:\n","truncated":false}}
%---
%[output:08a7c914]
%   data: {"dataType":"text","outputData":{"text":"                       Aspirin    Benzene    Ethanol    Caffeine    Acetaminophen\n                       _______    _______    _______    ________    _____________\n\n    Benzene ring        true       true       false      false          true     \n    Carboxylic acid     true       false      false      false          false    \n    Hydroxyl (-OH)      true       false      true       false          true     \n    Amide (-CONH)       false      false      false      false          true     \n    Nitrogen atom       false      false      false      true           true     \n\n","truncated":false}}
%---
%[output:48b4b33e]
%   data: {"dataType":"text","outputData":{"text":"[09:36:34][INFO]  Loaded 479 PAINS SMARTS patterns\n","truncated":false}}
%---
%[output:65443a13]
%   data: {"dataType":"text","outputData":{"text":"[09:36:34][INFO]  --- PAINS Screening (Total 479 Patterns) ---\n","truncated":false}}
%---
%[output:2f7f2809]
%   data: {"dataType":"text","outputData":{"text":"[09:36:42][INFO]    Ethanol         Clean\n[09:36:50][INFO]    Aspirin         Clean\n[09:36:58][INFO]    Quinone         PAINS (2 hits)\n[09:37:05][INFO]    Azo dye         PAINS (1 hits)\n","truncated":false}}
%---
%[output:18f19cfe]
%   data: {"dataType":"text","outputData":{"text":"          Name          Benzene ring    Hydroxyl    Carboxylic acid    Amine    Amide\n    ________________    ____________    ________    _______________    _____    _____\n\n    \"caffeine\"             false         false           false         false    false\n    \"nicotine\"             false         false           false         false    false\n    \"theobromine\"          false         false           false         false    false\n    \"aspirin\"              true          true            true          false    false\n    \"acetaminophen\"        true          true            false         true     true \n    \"ibuprofen\"            true          true            true          false    false\n    \"salicylic acid\"       true          true            true          false    false\n    \"ethanol\"              false         true            false         false    false\n    \"methanol\"             false         true            false         false    false\n    \"isopropanol\"          false         true            false         false    false\n    \"acetone\"              false         false           false         false    false\n    \"glycerol\"             false         true            false         false    false\n    \"sucrose\"              false         true            false         false    false\n    \"glucose\"              false         true            false         false    false\n    \"fructose\"             false         true            false         false    false\n    \"acetic acid\"          false         true            true          false    false\n    \"citric acid\"          false         true            true          false    false\n    \"lactic acid\"          false         true            true          false    false\n    \"benzoic acid\"         true          true            true          false    false\n    \"formic acid\"          false         true            true          false    false\n    \"vanillin\"             true          true            false         false    false\n    \"capsaicin\"            true          true            false         true     true \n    \"menthol\"              false         true            false         false    false\n    \"limonene\"             false         false           false         false    false\n    \"eugenol\"              true          true            false         false    false\n    \"benzaldehyde\"         true          false           false         false    false\n    \"linalool\"             false         true            false         false    false\n    \"ascorbic acid\"        false         true            false         false    false\n    \"urea\"                 false         false           false         true     true \n    \"carvone\"              false         false           false         false    false\n\n","truncated":false}}
%---
