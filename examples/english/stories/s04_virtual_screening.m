%[text] # S04: Virtual Screening — Finding Drug Candidates by Similarity
%[text] EasyMolKit Application Story — Layer 2
%[text] 
%[text] Can we automatically pick out unknown compounds with "similar shapes" to existing drugs from a library of hundreds?
%[text] In drug discovery research, the number of candidate compounds is enormous, and testing each one experimentally is not practical. **Virtual screening** is a method that uses fingerprints and similarity scores to computationally narrow down promising candidates before experimentation.
%[text] In this script, we will rank 200 FDA-approved drugs using Ibuprofen (a COX inhibitor) as a query, experiencing the entire workflow of ligand-based similarity screening.
%[text] ## Learning Objectives
%[text] - Understand the workflow of Ligand-Based Virtual Screening (LBVS)
%[text] - Perform efficient bulk Tanimoto searches with `emk.similarity.rankBy`
%[text] - Interpret Tanimoto score thresholds in the context of actual drug discovery
%[text] - Understand the impact of fingerprint types (ECFP4 vs MACCS) on the hit list \
%[text] ## Prerequisites
%[text] - Completion of F03 (Fingerprints) and F04 (Similarity)
%[text] - Recommended: S01 (Caffeine Cousins) and S02 (Lipinski Filter)
%[text] - RDKit installed (execute `emk.setup.install()` once)
%[text] - No additional Toolbox required (runs with MATLAB only) \
%[text] **Estimated Time**: 25–40 minutes | Execution Method: Run each section with Ctrl+Enter
%[text] **Data**
%[text] - `data/list/fda_drugs.csv` — 200 FDA-approved drugs (ChEMBL, CC-BY-SA 3.0)
%[text] - Columns: ChEMBLID, Name, SMILES, MolecularWeight, ALogP, HBondDonors, HBondAcceptors, TPSA, RotatableBonds, Source \
%[text] **References**
%[text] - Willett P (2006) Similarity-based virtual screening using 2D fingerprints. *Drug Discov Today* 11:1046-1053. [requires institutional access]
%[text] - Rogers D & Hahn M (2010) Extended-connectivity fingerprints. *J Chem Inf Model* 50:742-754. [requires institutional access]
%[text] - Johnson MA & Maggiora GM (1990) *Concepts and Applications of Molecular Similarity*. Wiley. [book]\
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
logInfo("S04: Setup complete");
%%
%[text] ## Section 1: Hit Compounds
%[text] ### Ligand-Based Virtual Screening (LBVS)
%[text]   Johnson & Maggiora (1990) Molecular Similarity Principle:
%[text]     "Structurally similar molecules tend to have similar bioactivity."
%[text]  LBVS leverages this principle:
%[text] \[Query (Hit Compound)\] → \[Fingerprint Calculation\] → \[Similarity Calculation with Library\] → \[Ranked Candidate List\] → \[Experimental Validation (Follow-up)\]
%[text]   Advantages over structure-based docking:
%[text] - No need for protein structure
%[text] - Speed to screen millions of compounds
%[text] - Effective for scaffold hopping within chemical series \
%[text]   Limitations:
%[text] - Cannot discover structurally novel active compounds ("scaffold hopping") that differ greatly from the query. \
%[text] Hit Compound: Ibuprofen
%[text]   Discovered in 1961 by Stewart Adams at Boots Laboratories (UK). Inhibits COX-1 and COX-2 enzymes, blocking prostaglandin synthesis. One of the most consumed drugs worldwide (OTC analgesic/antipyretic). Belongs to the arylpropionic acid (profen) class of NSAIDs.
HIT_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O";
HIT_NAME   = "Ibuprofen";

mol_hit = emk.mol.fromSmiles(HIT_SMILES);
logInfo("Hit Compound: %s", HIT_NAME);
logInfo("  SMILES    : %s", HIT_SMILES);
logInfo("  Heavy Atom Count  : %d", double(mol_hit.GetNumHeavyAtoms()));

desc_hit = emk.descriptor.calculate(mol_hit, ...
    ["MolWt", "LogP", "TPSA", "NumHDonors", "NumHAcceptors", "RingCount"]);

logInfo("Key Properties:");
logInfo("  MW   : %.1f Da   (Ro5 Limit: 500 -- PASS)", desc_hit.MolWt);
logInfo("  LogP : %.2f     (Ro5 Limit: 5   -- PASS)", desc_hit.LogP);
logInfo("  TPSA : %.1f A^2  (Oral Limit: 140 -- PASS)", desc_hit.TPSA);
logInfo("  HBD  : %d        (Ro5 Limit: 5   -- PASS)", desc_hit.NumHDonors);
logInfo("  HBA  : %d        (Ro5 Limit: 10  -- PASS)", desc_hit.NumHAcceptors);
logInfo("  Ring Count : %d        (One benzene ring; arylpropionic acid class)", desc_hit.RingCount);

figure("Name", "Ibuprofen (Hit Compound)", "Position", [100 100 440 380]);
emk.viz.draw2d(mol_hit, Title="Ibuprofen (Hit Compound)");
%[text] **✏️ Try It 1 — Check the meaning of Ro5 and LogP**
%[text] Ibuprofen has one chiral center (alpha carbon adjacent to the carboxyl group).
%[text] Pharmacologically active is only the (S)-enantiomer (COX inhibitor),
%[text] but the (R)-form is converted to the (S)-form by isomerase in the body,
%[text] so a racemic (R/S) mixture is sold as an OTC drug.
%[text] 
%[text] **Q1**: Does Ibuprofen pass Lipinski's Rule of Five?
%[text] Create a one-line table and call `emk.filter.lipinski()` to check.
%[text]  `d = emk.descriptor.calculate(mol_hit, ["MolWt","LogP","NumHDonors","NumHAcceptors"])`
%[text]  `t = struct2table(d)`
%[text]  `emk.filter.lipinski(t)`
%[text] Expected result: Violations 0, Pass_Ro5 = true.
%[text] 
%[text] **Q2**: Why does the calculated LogP (neutral form) differ from the "apparent LogP" in blood? 
%[text] Hint: LogP is measured for neutral species, but Ibuprofen (pKa ~4.4) is about 99% ionized at blood pH 7.4. D (distribution coefficient) accounts for ionization, so D << LogP.
% ... (Try writing code here)
%%
%[text] ## Section 2: Loading the Compound Library
%[text] ### What is a Compound Library?
%[text] In drug discovery research, a "compound library" is a collection of molecules available for screening (searching for drug candidates). The scale and nature of the library can vary greatly depending on the purpose and stage.
%[text] - Focused Library (about 100 to 10,000 compounds): A collection of analogs designed by chemists targeting specific proteins or chemical scaffolds.  \
%[text] - Corporate Library (about 10,000 to 1,000,000 compounds): A vast asset accumulated independently by pharmaceutical companies through past drug projects.  \
%[text] - Commercial Library (about 1 billion to 10 billion compounds): An ultra-large compound space provided by chemical suppliers (e.g., \*Enamine REAL Space\*). \
%[text] Includes those designed on computers and synthesizable upon request.
%[text] 
%[text] Approach in this tutorial: Drug Repurposing (Repositioning) In this script, instead of a huge library, we use 200 FDA-approved drugs as a "proxy library." Screening existing approved drugs is called a "drug repurposing" strategy. Since safety in humans is already confirmed, early safety trials can be bypassed, significantly shortening the development period (reducing drug discovery from the usual 10-15 years to 3-5 years).
%[text] 
%[text] Representative examples of successful repurposing:
%[text] - **Thalidomide**: 1950s "sedative and hypnotic" → Reapproved in 2006 as a treatment for "multiple myeloma (blood cancer)."
%[text] - **Sildenafil**: Originally a treatment for "angina" → Leveraging side effects, became a hit as "erectile dysfunction (Viagra)" in 1998.
%[text] - **Metformin**: A staple "diabetes" treatment → Recently researched worldwide for anticancer effects and anti-aging (longevity research).
LIBRARY_FILE = fullfile(projectRoot, "data", "list", "fda_drugs.csv");
library = readtable(LIBRARY_FILE, "TextType", "string");

logInfo("Library loaded: %d compounds", height(library));
%[text] Library Summary
logInfo("Library Descriptor Summary:");
logInfo("  MW    -- Min: %5.1f  Median: %5.1f  Max: %6.1f  Da", ...
    min(library.MolecularWeight), median(library.MolecularWeight), max(library.MolecularWeight));
logInfo("  ALogP -- Min: %5.2f  Median: %5.2f  Max: %5.2f", ...
    min(library.ALogP), median(library.ALogP), max(library.ALogP));
%[text] Check if the hit compound is in the library
hitInLib = strcmpi(library.Name, HIT_NAME);
if any(hitInLib)
    logInfo("  Note: %s is present in the library (self-match at the top T=1.0 appears)", ...
        HIT_NAME);
else
    logInfo("  Note: %s is not in the library -- all hits are true analogs", HIT_NAME);
end
%[text] **✏️ Try It 2 — Explore the Library Contents**
%[text] **Q1**: How many unique drug names are in the library?
%[text] Hint: `numel(unique(library.Name))`
%[text] 
%[text] **Q2**: What is the heaviest compound (maximum MolecularWeight)?
%[text] Hint: `library(library.MolecularWeight == max(library.MolecularWeight), :)`
%[text] Expectation: A large natural product-derived drug (such as macrolide or peptide).
%[text] 
%[text] **Q3**: Both ALogP (Ghose-Crippen) of the library and LogP (Wildman-Crippen) of `emk.descriptor.calculate` estimate octanol-water partitioning. Are the two values close for Ibuprofen?
%[text] Hint: `desc_hit.LogP` vs ALogP of Ibuprofen in the library (if present).
% ... (Try writing code here)
%%
%[text] ## Section 3: Calculate ECFP4 Fingerprints for the Library
%[text] ### **How ECFP4 (Morgan Fingerprint) Works**
%[text] Computers cannot directly understand the "picture" of a molecule. Therefore, fingerprints are used to convert the molecular structure into a digital code of 0s and 1s (bit vector).
%[text] The representative ECFP4 (Extended Connectivity Fingerprint) explores and records the chemical environment around each atom in the molecule as a starting point, based on "bond radius (number of hops)" (Morgan algorithm).
%[text] - Radius 0: Information about the atom itself (element type, charge, aromaticity, etc.)
%[text] - Radius 1: The atom of interest + the environment of directly adjacent atoms (1 bond away)
%[text] - Radius 2: The atom of interest + the environment up to 2 hops away (equivalent to ECFP4, covering a diameter of 4 bonds) \
%[text] All these explored local environments are mathematically hashed (compressed) into a fixed-length vector of 2048 bits. Two molecules sharing many similar local environments are more likely to have "bits in the same position turned ON (1)" in their fingerprints. As a result, the Tanimoto score, discussed later, is calculated to be high.
%[text] Main Uses of ECFP4
%[text] In the field of drug discovery, ECFP4 is used worldwide for the following purposes.
%[text] 1. Similarity Search: Instantly find compounds with structures similar to the query molecule (the main theme of this exercise).
%[text] 2. Compound Clustering: Classify a vast library into groups with similar structures and evaluate diversity.
%[text] 3. Machine Learning (AI Drug Discovery): Use as input data (descriptors) representing molecular features for AI models. \
logInfo("Calculating ECFP4 fingerprints for %d compounds in the library...", height(library));

lib_fps   = cell(1, height(library));
lib_valid = true(1, height(library));

for i = 1:height(library)
    smi = library.SMILES(i);
    if ~emk.mol.isValid(smi)
        logWarn("  Skipping line %d (%s): Invalid SMILES", i, library.Name(i));
        lib_valid(i) = false;
        continue;
    end
    mol_lib    = emk.mol.fromSmiles(smi);
    lib_fps{i} = emk.fingerprint.morgan(mol_lib);
    logProgress(i, height(library), "FP");
end

nLibValid = sum(lib_valid);
logInfo("Fingerprint preparation complete: %d / %d compounds", nLibValid, height(library));
%[text] Extract Valid Subset
lib_fps_valid    = lib_fps(lib_valid);
lib_names_valid  = library.Name(lib_valid);
lib_smiles_valid = library.SMILES(lib_valid);
lib_mw_valid     = library.MolecularWeight(lib_valid);
lib_alogp_valid  = library.ALogP(lib_valid);
%[text] Calculate Query Fingerprint
fp_hit = emk.fingerprint.morgan(mol_hit);
nOnBits = sum(emk.fingerprint.toArray(fp_hit));
logInfo("Ibuprofen ECFP4: %d ON bits / 2048 total (density %.1f%%)", ...
    nOnBits, 100*nOnBits/2048);
%[text] **✏️ Try It 3 — Compare Fingerprint Density**
%[text] Ibuprofen is a monocyclic compound with 13 heavy atoms, while Caffeine (S01) is a fused bicyclic compound with 14 heavy atoms. Which has more ON bits in ECFP4?
%[text] **Q1**: Calculate and compare the fingerprint of Caffeine.
%[text]  `mol_caf = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")`
%[text]  `fp_caf  = emk.fingerprint.morgan(mol_caf)`
%[text]  `sum(emk.fingerprint.toArray(fp_caf))`
%[text] Expected: Caffeine has more ON bits. The fused bicyclic system creates many different local environments.
%[text] 
%[text] **Q2**: What happens if you increase the radius from 2 to 3?
%[text]  `fp_r3 = emk.fingerprint.morgan(mol_hit, Radius=3)`
%[text]  `sum(emk.fingerprint.toArray(fp_r3))`
%[text] Expected: Usually, the number of ON bits remains the same or decreases. With radius 3, each neighborhood becomes larger and more unique, reducing collisions, but shared bits may be lost even in very similar compounds.
% ... (Try writing code here)

%%
%[text] ## Section 4: Performing Virtual Screening
%[text] ### What is Bulk Tanimoto Similarity Search?
%[text] The most widely used metric to quantify how similar two molecules (bit vectors) are is the Tanimoto Similarity (Tanimoto Similarity / Jaccard Index).
%[text] $&dollar&;&dollar&;T(A, B) = \\frac{|A \\cap B|}{|A \\cup B|} = \\frac{\\text{Number of bits ON in both molecules A and B}}{\\text{Number of bits ON in at least one of molecules A or B}}&dollar&;&dollar&; $
%[text] $&dollar&;T = 1.0&dollar&;$: Fingerprints are identical (same compound, or extremely similar scaffold indistinguishable by ECFP4).
%[text] $&dollar&;T = 0.0&dollar&;$: No common bits (completely different structures).
%[text] ### Practical Interpretation of Scores (Thresholds)
%[text] In medicinal chemistry, the Tanimoto score for ECFP4 is empirically interpreted as follows (Willett 2006). It is an important guideline when reviewing hit lists.
%[text] - $&dollar&;T \\ge 0.85&dollar&; $ (Very Similar): Analog with almost the same scaffold, only slightly different substituents.
%[text] - $&dollar&;T \\ge 0.65&dollar&;$ (Similar): Shares the same chemical class or a common core ring system.
%[text] - $&dollar&;T \\ge 0.40&dollar&; $ (Moderate): Has common functional groups, but different scaffold.
%[text] - $&dollar&;T \< 0.40&dollar&;$ (Low/Non-similar): Few structural commonalities.
%[text] ### Efficiency of Calculations: `emk.similarity.rankBy`
%[text] The `emk.similarity.rankBy` function used in this script internally calls RDKit's `BulkTanimotoSimilarity`. Instead of sending 200 compounds one by one to Python for calculation, it screens the entire library in a single batch process. This minimizes the number of inter-process communications (IPC) between MATLAB and Python to one, allowing it to operate extremely fast even for tens of thousands of compounds.
tic;
vs_result = emk.similarity.rankBy(fp_hit, lib_fps_valid);
tElapsed = toc;

logInfo("Virtual screening completed: Ranked %d compounds in %.3f seconds", ...
    numel(vs_result.Scores), tElapsed);

%[text] Display the top 10 hits with interpretation
TOP_DISPLAY = 10;
logInfo("--- Top %d hits for %s ---", TOP_DISPLAY, HIT_NAME);
logInfo("%-5s  %-28s  %9s  %s", "Rank", "Name", "Tanimoto", "Interpretation");
for k = 1:min(TOP_DISPLAY, numel(vs_result.Scores))
    idx   = vs_result.Indices(k);
    score = vs_result.Scores(k);
    name  = lib_names_valid(idx);

    if score >= 1.0 - 1e-6
        interp = "<-- Self-match (T=1.0)";
    elseif score >= 0.85
        interp = "Very Similar";
    elseif score >= 0.65
        interp = "Similar";
    elseif score >= 0.40
        interp = "Moderate";
    else
        interp = "Low";
    end

    logInfo("  %2d.  %-28s  %9.4f  %s", k, name, score, interp);
end

%[text] **✏️ Try It 4 — Find Other NSAIDs**
%[text] Do the top hits include other NSAIDs (Flurbiprofen, Ketoprofen)?
%[text] These share the same arylpropionic acid core (arene—CH(CH₃)—COOH) as Ibuprofen.
%[text] 
%[text] **Q1**: Check the top 2 hits excluding self-matches.
%[text] Expected: FLURBIPROFEN (T~0.40), KETOPROFEN (T~0.39). The reason they remain around T~0.40 is that fluorine substituents and extra rings alter the local environment, making ECFP4 distinguish them. (T \> 0.30 with ECFP4 still holds meaning as drug-like analogs. 0.65/0.85 are guidelines.)
%[text] **Note**: AMPHETAMINE SULFATE (T~0.35) may appear in 4th place.
%[text] The local environment of a phenyl ring plus a short alkyl chain overlaps with Ibuprofen's substructure, but it lacks carboxylic acid and is a CNS stimulant with no anti-inflammatory activity. This is a typical example of a "structural false positive" with ECFP4. This is why chemists must verify structures before proceeding to experiments.
%[text] 
%[text] **Q2**: What rank is Indomethacin (an acetic acid NSAID)?
%[text] Expected: Lower Tanimoto than the profen class. This is because ECFP4 captures scaffold differences.
%[text] 
%[text] **Q3**: How does the hit list change when querying with Aspirin (`"CC(=O)Oc1ccccc1C(=O)O"`)?
%[text]  `mol_asp = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O")`
%[text]  `fp_asp  = emk.fingerprint.morgan(mol_asp)`
%[text]  `res_asp = emk.similarity.rankBy(fp_asp, lib_fps_valid)`
%[text]  `lib_names_valid(res_asp.Indices(1:5))`
% ... (Try writing the code here)
%%
%[text] ## Section 5: Create a Hit Report
%[text] Create a structured result table: Rank, Name, Tanimoto, MW, ALogP. This will be the report handed to the experimental team.
REPORT_TOP = 15;
rankVec  = (1:REPORT_TOP)';
nameVec  = strings(REPORT_TOP, 1);
scoreVec = zeros(REPORT_TOP, 1);
mwVec    = zeros(REPORT_TOP, 1);
logpVec  = zeros(REPORT_TOP, 1);
smiVec   = strings(REPORT_TOP, 1);

for k = 1:REPORT_TOP
    if k > numel(vs_result.Scores); break; end
    idx         = vs_result.Indices(k);
    nameVec(k)  = lib_names_valid(idx);
    scoreVec(k) = vs_result.Scores(k);
    mwVec(k)    = lib_mw_valid(idx);
    logpVec(k)  = lib_alogp_valid(idx);
    smiVec(k)   = lib_smiles_valid(idx);
end

hits_tbl = table(rankVec, nameVec, scoreVec, mwVec, logpVec, smiVec, ...
    'VariableNames', ["Rank", "Name", "Tanimoto", "MW_Da", "ALogP", "SMILES"]);

logInfo("Top %d hit report:", REPORT_TOP);
disp(hits_tbl(:, 1:5));   % Display without SMILES column

%[text] Tanimoto Threshold Summary
logInfo("Breakdown by threshold (Top %d):", REPORT_TOP);
logInfo("  T >= 0.85 (Very Similar): %d", sum(scoreVec >= 0.85));
logInfo("  T >= 0.65 (Similar)     : %d", sum(scoreVec >= 0.65 & scoreVec < 0.85));
logInfo("  T >= 0.40 (Moderate)    : %d", sum(scoreVec >= 0.40 & scoreVec < 0.65));
logInfo("  T <  0.40 (Low)         : %d", sum(scoreVec < 0.40 & scoreVec > 0));

%[text] **✏️ Try It 5 — Analyze the Hit Report**
%[text] **Q1**: Sort `hits_tbl` by MW and check if the top hits cluster around Ibuprofen (206 Da).
%[text] Hint: `sortrows(hits_tbl, "MW_Da")`
%[text] 
%[text] **Q2**: Which hit has the highest ALogP? Is it still drug-like (ALogP \< 5)?
%[text] High ALogP often leads to poor water solubility, posing formulation challenges despite activity.
%[text] 
%[text] **Q3**: Can you add a Ro5 compliance column to the hit table?
%[text] Change the table column names and call `emk.filter.lipinski()`.
%[text] (MW, LogP, NumHDonors, NumHAcceptors are needed. HBD/HBA are in the CSV.)
% ... (Try writing the code here)
%%
%[text] ## Section 6: Score Distribution
%[text] ### Chemical Space Coverage (Library Diversity Evaluation)
%[text] It is very important to overview the distribution (histogram) of Tanimoto scores for all compounds in the library, not just the top hits. This allows you to see the position of the query (the starting point of the search) in relation to the prepared library (coverage in chemical space).
%[text] 
%[text] If there is a sharp peak near $&dollar&;T = 0&dollar&; $ and no tail on the high score side, it means the query molecule is an "isolated outlier" in that library. Since similar compounds do not exist in the library, it is a difficult case for expansion (follow-up) through screening.
%[text] If the distribution spreads gently towards the high $&dollar&;T&dollar&; $ side (right side), it indicates that the library is rich in "close analogs" related to the query. It is a dense space suitable for tracking activity changes (SAR: Structure-Activity Relationship) when the structure is slightly altered.
allScores = vs_result.Scores;

figure("Name", "Virtual Screening Score Distribution", "Color", "white", "Position", [100 100 580 420]);
histogram(allScores, 30, "FaceColor", [0.2 0.55 0.85], ...
    "EdgeColor", "white", "FaceAlpha", 0.8, ...
    "DisplayName", "Tanimoto Score Distribution");
hold on;
%[text] Mark major thresholds
xline(0.85, "r--", "LineWidth", 1.5, "DisplayName", "Very Similar (0.85)");
xline(0.65, "m--", "LineWidth", 1.5, "DisplayName", "Similar (0.65)");
xline(0.40, "g--", "LineWidth", 1.5, "DisplayName", "Moderate (0.40)");

hold off;
xlabel("Tanimoto Score (ECFP4)");
ylabel("Number of Library Compounds");
title(sprintf("Similarity Distribution: %s vs FDA Drugs 200", HIT_NAME));
legend("Location", "northeast");
grid("on");
%[text] Count by zone
nVSim  = sum(allScores >= 0.85);
nSim   = sum(allScores >= 0.65 & allScores < 0.85);
nMod   = sum(allScores >= 0.40 & allScores < 0.65);
nLow   = sum(allScores < 0.40);

logInfo("%d Overall Score Distribution of Compounds:", nLibValid);
logInfo("  Very Similar (T>=0.85): %3d  (%.0f%%)", nVSim, 100*nVSim/nLibValid);
logInfo("  Similar       (T>=0.65): %3d  (%.0f%%)", nSim,  100*nSim /nLibValid);
logInfo("  Moderate     (T>=0.40): %3d  (%.0f%%)", nMod,  100*nMod /nLibValid);
logInfo("  Low         (T< 0.40): %3d  (%.0f%%)", nLow,  100*nLow /nLibValid);
%[text] **✏️ Try It 6 — Compare Score Distributions**
%[text] **Q1**: What is the proportion of the library of 200 drugs with T \>= 0.65 with Ibuprofen?
%[text] Is Ibuprofen well represented in this library (many neighbors)?
%[text] Or is it an outlier (few neighbors)?
%[text] 
%[text] **Q2**: Change the query to Caffeine and compare the two histograms.
%[text] Which has more high Tanimoto hits?
%[text]  `mol_caf = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")`
%[text]  `fp_caf  = emk.fingerprint.morgan(mol_caf)`
%[text]  `res_caf = emk.similarity.rankBy(fp_caf, lib_fps_valid)`
%[text]  `figure; histogram(res_caf.Scores, 30); title("Caffeine vs FDA Drugs")`
% ... (Try writing code here)
%%
%[text] ## Section 7: Visualize the Top 3 Hits
%[text] Draw the top 3 non-self matches side by side to check structural similarity.
nonSelf = find(vs_result.Scores < 1.0 - 1e-6);
drawK   = nonSelf(1:min(3, numel(nonSelf)));

logInfo("Drawing the top %d non-self matches:", numel(drawK));
figure("Name", "Top Hits", "Position", [100 100 960 340]);
for j = 1:numel(drawK)
    k     = drawK(j);
    idx   = vs_result.Indices(k);
    name  = lib_names_valid(idx);
    smi   = lib_smiles_valid(idx);
    score = vs_result.Scores(k);
    mol_j = emk.mol.fromSmiles(smi);
    subplot(1, numel(drawK), j);
    emk.viz.draw2d(mol_j, Title=sprintf("Rank %d: %s (T=%.3f)", k, name, score));
    logInfo("  Rank %d: %s (T=%.4f, MW=%.1f, ALogP=%.2f)", ...
        k, name, score, lib_mw_valid(idx), lib_alogp_valid(idx));
end

%[text] **✏️ Try It 7 — Visually Compare Structures**
%[text] **Q1**: Do the top 3 hits visually resemble Ibuprofen?
%[text] What are the common structural features?
%[text] Expected: Shared aryl ring + propionic acid side chain ("profen" scaffold).
%[text] 
%[text] **Q2**: Try drawing Ibuprofen and the top hit side by side in a 1×2 figure.
%[text]  `figure("Position", [100 100 750 300])`
%[text]  `subplot(1,2,1); emk.viz.draw2d(mol_hit, Title="Ibuprofen (Query)")`
%[text]  `subplot(1,2,2); emk.viz.draw2d(emk.mol.fromSmiles(smiVec(1)), Title=nameVec(1))`
%[text] (Note: `draw2d` opens its own figure by default. Using subplot allows you to
%[text] compare two structures side by side.)
% ... (Try writing your code here)
%%
%[text] ## Section 8: ECFP4 vs MACCS Keys -- Is the choice of fingerprint important?
%[text] 
%[text] ### Impact of fingerprint choice on hit lists
%[text] There are several types of "fingerprints" that digitize molecules, and the choice can significantly affect screening results (hit lists). Here, we summarize the two characteristics compared this time.
%[text] - ECFP4 (2048 bits / circular descriptors) Encodes the local atomic environment (connectivity) around each atom in detail.
%[text] - MACCS Keys (166 bits / fragment-based descriptors) Encodes the presence or absence of "specific structural fragments (substructures)" defined by medicinal chemists as 0 or 1. Examples: "Is there a carboxylic acid in the molecule?", "Is there an aromatic ring?", "Is there a halogen atom?", etc.
%[text] ### Characteristics and usage of each
%[text] - ECFP4: Very sensitive to subtle structural differences in molecules (such as the position of substituents and changes in ring size).
%[text] - MACCS Keys: Insensitive to changes in structural details (decorations) but suitable for broadly capturing "chemical classes (commonality of skeletons)".
%[text] ### Why does using MACCS Keys tend to result in higher Tanimoto scores?
%[text] When actually calculated, MACCS tends to yield higher (inflated) scores overall. The reasons are as follows:
%[text] 1. Overwhelming difference in the number of bits (166 bits $&dollar&;\\ll&dollar&;$ 2048 bits): Due to the shorter bit length (coarser resolution), structural "mismatches" are less likely to be counted.
%[text] 2. Shared basic fragments: The query ibuprofen and other non-steroidal anti-inflammatory drugs (NSAIDs) in the library share many of the major fragments defined by MACCS, such as "having an aromatic ring" and "having a carboxylic acid".
%[text] ### Interpretation of screening results
%[text] - If both ECFP4 and MACCS Keys rank high: It means that not only the broad chemical class but also the local environment details are very similar, supporting "strong structural similarity (kinship)".
%[text] - If ranked high by MACCS but low by ECFP4: It indicates that while they share the broad class (properties) of "aromatic compounds with carboxylic acid", the detailed skeleton and substituent attachment are completely different.
logInfo("Calculating MACCS fingerprints for the library...");

maccs_fp_hit = emk.fingerprint.maccs(mol_hit);
maccs_lib_fps = cell(1, numel(lib_fps_valid));

for i = 1:numel(lib_fps_valid)
    mol_i = emk.mol.fromSmiles(lib_smiles_valid(i));
    maccs_lib_fps{i} = emk.fingerprint.maccs(mol_i);
    logProgress(i, numel(lib_fps_valid), "MACCS");
end

maccs_result = emk.similarity.rankBy(maccs_fp_hit, maccs_lib_fps);

%[text] Compare the top 10 results with ECFP4 and MACCS
COMPARE_N = 10;
logInfo("--- ECFP4 vs MACCS Keys: Top %d Comparison ---", COMPARE_N);
logInfo("%-5s  %-24s  %9s  %-24s  %9s", ...
    "Rank", "ECFP4 Hit", "T(ECFP4)", "MACCS Hit", "T(MACCS)");
for k = 1:COMPARE_N
    ie = vs_result.Indices(k);
    im = maccs_result.Indices(k);
    logInfo("  %2d.  %-24s  %9.4f  %-24s  %9.4f", k, ...
        lib_names_valid(ie), vs_result.Scores(k), ...
        lib_names_valid(im), maccs_result.Scores(k));
end

%[text] Rank correlation: Agreement of two fingerprints in the top 20
TOP_CORR = 20;
ecfp4_top = arrayfun(@(k) lib_names_valid(vs_result.Indices(k)),    1:TOP_CORR);
maccs_top = arrayfun(@(k) lib_names_valid(maccs_result.Indices(k)), 1:TOP_CORR);
nOverlap  = numel(intersect(ecfp4_top, maccs_top));
logInfo("Overlap of top %d between ECFP4 and MACCS: %d / %d (%.0f%%)", ...
    TOP_CORR, nOverlap, TOP_CORR, 100*nOverlap/TOP_CORR);

%[text] **✏️ Try It 8 — Compare the choice of fingerprints**
%[text] **Q1**: Do ECFP4 and MACCS return the same top 3 compounds?
%[text] If not, what causes the difference?
%[text] (MACCS Keys include generic fragments like "has COOH".
%[text] Any carboxylic acid will score high with MACCS, but ECFP4 only returns a high score if the COOH is
%[text] in a similar ring context.)
%[text] 
%[text] **Q2**: The MACCS Tanimoto scores for top hits tend to be higher than ECFP4.
%[text] Does this mean MACCS is "better"?
%[text] (Hint: The high scores are due to fewer bits and less mismatch, not higher accuracy.
%[text] ECFP4, sensitive to scaffold differences, is preferred for lead optimization.)
%[text] 
%[text] **Q3**: Try calculating similarity using the Dice index instead of Tanimoto.
%[text]  `res_dice = emk.similarity.rankBy(fp_hit, lib_fps_valid, Metric="dice")`
%[text] Dice = 2|A AND B| / (|A| + |B|). For bit vectors, Dice \>= Tanimoto always holds.
%[text] How does the ranking change?
% ... (Try writing the code here)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
