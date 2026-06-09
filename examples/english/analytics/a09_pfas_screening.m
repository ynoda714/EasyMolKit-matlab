%[text] # A09: PFAS and Environmental Chemical Screening
%[text]
%[text] Consider how to determine which chemicals are PFAS, the "forever chemicals," when unknown chemicals are detected in groundwater near a manufacturing facility.
%[text] PFAS are a group of over 9,000 synthetic chemicals used since the 1940s in non-stick coatings and firefighting foams. Due to the C-F bond (about 544 kJ/mol), they have properties that make them resistant to biological and environmental degradation.
%[text] By combining SMARTS patterns and physicochemical scoring, you can automatically narrow down numerous candidates and prioritize them according to their level of concern.
%[text] In this script, you will experience PFAS screening, persistence scoring, weight optimization, and Tanimoto similarity analysis from the perspective of an environmental toxicologist at a water quality management agency.
%[text]
%[text] **Learning Objectives**
%[text] - Understand the structural definition of PFAS (OECD 2021 criteria)
%[text] - Apply multi-pattern SMARTS screening using `emk.mol.hasSubstruct()`
%[text] - Learn how to encode expert knowledge as quantitative concern scores
%[text] - Solve constrained weight optimization problems using fmincon
%[text] - Interpret Tanimoto similarity heatmaps for structural clustering
%[text] - Report screening results as a prioritized concern table
%[text]
%[text] **Workflow (4 Steps):**
%[text] 1. **PFAS Flagging** — Detect perfluoroalkyl chains, CF3 terminals, and sulfonyl head groups using SMARTS
%[text] 2. **Persistence Scoring** — Estimate environmental persistence using three proxies: LogP, number of fluorines, and TPSA
%[text] 3. **Multi-Criteria Prioritization** — Use fmincon to find the weight set that best fits the reference rubric
%[text] 4. **Report** — Create a ranked concern table, scatter plots, and Tanimoto similarity heatmaps
%[text]
%[text] **Prerequisites**
%[text] - Completion of A07 (Scaffold Analysis) — Understanding of SMARTS and substructure concepts
%[text] - Optimization Toolbox (fmincon) — Required for Section 4 (automatic fallback to equal weights if unavailable)
%[text] - Statistics and Machine Learning Toolbox — Required for Section 5 (`clusterdata`, `linkage`, `dendrogram`)
%[text] - Both toolboxes are included in MATLAB Online Basic (free tier)
%[text]
%[text] Data: All molecules are defined inline (SMILES literals), no external files needed (US EPA CompTox + OECD 2021 criteria).
%[text]
%[text] Estimated Time Required: 35-50 minutes | Execution Method: Run each section one by one with Ctrl+Enter
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

%[text] Warm up the Python/RDKit process before main execution
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;

%[text] Check available toolboxes
hasOptTbx   = license("test", "optimization_toolbox");
hasStatsTbx = license("test", "statistics_toolbox");

if hasOptTbx
    logInfo("A09: Optimization Toolbox detected -- weight optimization enabled.");
else
    logWarn("A09: Optimization Toolbox not detected.");
    logWarn("     Section 4 will use equal weight fallback.");
end

if hasStatsTbx
    logInfo("A09: Statistics and ML Toolbox detected -- hierarchical clustering enabled.");
else
    logWarn("A09: Statistics and ML Toolbox not detected.");
    logWarn("     Clustering in Section 5 will be skipped.");
end

logSection("A09", "Section 0: Setup", "Analytics L3");
logInfo("A09: Setup complete.");
%%
%[text] ## Section 1: Definition of Test Chemical Inventory
%[text]
%[text] Setup is complete. First, we define an inventory of 20 chemicals to be analyzed.
%[text] This is a realistic set including four classes of PFAS and non-PFAS controls.
%[text]
%[text] ### Concept: What are PFAS?
%[text]
%[text] OECD (2021) Definition:
%[text] PFAS are substances that contain at least one perfluorinated methyl (-CF3) or methylene (-CF2-) carbon atom, except when all C-F bonds are directly attached to heteroatoms (e.g., -CF2Cl, SF6 are excluded).
%[text]
%[text] Major Structural Classes:
%[text]
%[text] 1. Perfluoroalkyl Carboxylic Acids (PFCAs)
%[text] General formula: F(CF2)_n-COOH
%[text] Examples: PFOA (n=7), PFBA (n=3), PFHxA (n=5)
%[text] SMARTS flag: Perfluoroalkyl chain [C](F)(F)(F) or [CF2] repeat
%[text]
%[text] 2. Perfluoroalkyl Sulfonic Acids (PFSAs)
%[text] General formula: F(CF2)_n-SO3H
%[text] Examples: PFOS (n=8), PFBS (n=4), PFHxS (n=6)
%[text] SMARTS flag: Combined with perfluoroalkyl chain -S(=O)(=O)O
%[text]
%[text] 3. Fluorotelomer Substances (FTS / FTOH)
%[text] Structure: F(CF2)_n-CH2CH2- (polyfluoroalkyl: some C-H bonds present)
%[text] These are precursor chemicals that degrade to PFCA in the environment. In this set, we use fluorotelomer alcohols (FTOH) with a head group of -OH, and the class label is unified as "FTS".
%[text] SMARTS flag: -CF2-CH2- junction
%[text]
%[text] 4. Non-PFAS Fluorinated Compounds
%[text] Includes drugs (Fluoxetine, Ciprofloxacin) and pesticides (Flurbiprofen).
%[text] Contains isolated C-F bonds but no perfluoroalkyl chains.
%[text] SMARTS flag: Fails perfluoroalkyl chain test -> classified as NonPFAS
%[text]
%[text] The test set of 20 chemicals covers all four classes and includes some non-PFAS controls from the US EPA CompTox PFAS Universe list.
%[text]
%[text] Each entry: {Display Name, SMILES, True Class Label}
%[text] Class Labels: "PFCA", "PFSA", "FTS", "NonPFAS"
logSection("A09", "Section 1: Definition of Test Chemical Inventory", "Analytics L3");
CHEMICALS = { ...
    % --- PFCAs (Perfluoroalkyl Carboxylic Acids) ---
    "PFBA",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)F",                           "PFCA";    ...  % C4 PFCA: COOH-CF2-CF2-CF3 (4 carbons, fully fluorinated)
    "PFHxA",   "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",           "PFCA";    ...  % C6 PFCA
    "PFOA",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "PFCA"; ...  % C8 PFCA (OECD Annex B substance)
    "PFNA",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F","PFCA"; ... % C9 PFCA
    "PFDA",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F","PFCA"; ... % C10 PFCA
    % --- PFSAs (Perfluoroalkyl Sulfonic Acids) ---
    "PFBS",    "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",               "PFSA";    ...  % C4 PFSA: SO3H-CF2-CF2-CF2-CF3 (4 carbons, fully fluorinated)
    "PFHxS",   "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F","PFSA";    ...  % C6 PFSA
    "PFOS",    "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F","PFSA"; ... % C8 PFSA (Stockholm Convention Annex B)
    % --- Fluorotelomers and Polyfluoroalkyl Precursors ---
    "6:2FTS",  "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",       "FTS";     ...  % 6:2 Fluorotelomer Sulfonate
    "8:2FTS",  "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F","FTS"; ... % 8:2 Fluorotelomer
    "4:2FTS",  "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                     "FTS";     ...  % 4:2 Fluorotelomer
    % --- Non-PFAS Fluorinated Pharmaceuticals and Pesticides (Controls) ---
    "Fluoxetine", "CNCCC(c1ccccc1)Oc1ccc(cc1)C(F)(F)F",                 "NonPFAS"; ...  % Antidepressant (1 trifluoromethyl)
    "Ciprofloxacin","OC(=O)c1cn(C2CCNCC2)c2cc(F)c(nc12)N1CCCC1",        "NonPFAS"; ...  % Antibiotic (1 C-F)
    "Flurbiprofen","OC(=O)C(C)c1ccc(cc1)-c1cccc(F)c1",                  "NonPFAS"; ...  % NSAID (1 C-F)
    "Diflunisal","OC(=O)c1ccc(cc1O)-c1ccc(F)cc1F",                      "NonPFAS"; ...  % NSAID (2 isolated C-F)
    "Trifluridine","OC1C(CO)OC(n2cc(c(=O)[nH]2)C(F)(F)F)C1F",           "NonPFAS"; ...  % Antiviral (mixed halogenation)
    "Halothane","FC(F)(F)C(Cl)Br",                                       "NonPFAS"; ...  % Anesthetic (3F on 1C, no chain)
    "Sevoflurane","FC(F)(F)C(F)OCC(F)(F)F",                              "NonPFAS"; ...  % Volatile anesthetic
    "Flecainide","OCC(NC(=O)c1cc(OCC(F)(F)F)c(cc1OCC(F)(F)F)OCC(F)(F)F)CC", "NonPFAS"; ... % Antiarrhythmic (OCH2CF3 group)
    "Ethanol",  "CCO",                                                   "NonPFAS"; ...  % Negative control (no F)
};

nChem   = size(CHEMICALS, 1);
names   = string(CHEMICALS(:,1));
smiles  = string(CHEMICALS(:,2));
trueCls = string(CHEMICALS(:,3));

logInfo("Inventory: %d chemicals (PFCA %d, PFSA %d, FTS %d, NonPFAS %d)", ...
    nChem, sum(trueCls=="PFCA"), sum(trueCls=="PFSA"), ...
    sum(trueCls=="FTS"),  sum(trueCls=="NonPFAS"));

%[text] Parsing molecules.
mols  = cell(1, nChem);
valid = false(1, nChem);
for k = 1:nChem
    try
        mols{k} = emk.mol.fromSmiles(smiles(k));
        valid(k) = true;
    catch ME
        logWarn("Cannot parse %s: %s", names(k), ME.message);
    end
end
logInfo("Parsed %d / %d molecules.", sum(valid), nChem);

%[text] **💡 Observation Point 1 — Check Carbon Chain Length and Classification**
%[text] PFOS has 8 fluorinated carbons, while PFHxS has only 6. The numeric prefix in the name (PFHxS = 6, PFOS = 8) indicates the carbon chain length.
%[text] Count the -CF2- units in each SMILES above to verify if the carbon count matches the name.
%[text] Fluoxetine has a -CF3 group but is classified as NonPFAS. Why?
%[text] (Hint: -CF3 alone does not form a perfluoroalkyl chain. Continuous fully fluorinated carbons are required.)
%[text] Halothane (`FC(F)(F)C(Cl)Br`) has 3 fluorines on one carbon. Is it PFAS by OECD definition?
%[text] (Does a single -CF3 with adjacent non-F halogens fit "at least one -CF3 or -CF2-"?)
% ... (Try writing code here)
%%
%[text] ## Section 2: SMARTS-Based PFAS Flagging
%[text]
%[text] The inventory has been defined. Next, we will automatically determine whether each molecule is a PFAS (perfluoroalkyl compound) using SMARTS patterns.
%[text] Let's check which compounds match each pattern.
%[text]
%[text] ### Concept: SMARTS Patterns for PFAS Detection
%[text]
%[text] SMARTS (SMiles ARbitrary Target Specification) patterns encode structural queries with constraints on atoms and bonds. Unlike SMILES, which describes specific molecules, SMARTS describes structural motifs common to many molecules.
%[text]
%[text] Key SMARTS atoms used below:
%[text] [#6]    -- Any carbon atom (atomic number 6)
%[text] [F]     -- Fluorine (element symbol, no charge constraints)
%[text] [#16]   -- Any sulfur atom (atomic number 16)
%[text] [!F]    -- Any atom that is not fluorine
%[text] [$(...)] -- Recursive SMARTS (nested subquery)
%[text]
%[text] Three complementary patterns cover the OECD PFAS structural definition:
%[text]
%[text] (A) Perfluoroalkyl Chain  [#6](F)(F)-[#6](F)(F)
%[text] Matches when a carbon bonded to exactly two Fs is adjacent to another identical carbon.
%[text] Signature of a fully fluorinated carbon chain -CF2-CF2-.
%[text] Requires at least two consecutive CF2 units.
%[text]
%[text] (B) Trifluoromethyl End  [#6](F)(F)F
%[text] Matches -CF3. The end of a perfluoroalkyl chain.
%[text] Also matches isolated CF3 groups (control), so (B) alone is too broad. Reliable PFAS flagging requires combination with (A).
%[text]
%[text] (C) Sulfonyl Head Group  [#16](=O)(=O)
%[text] Matches -S(=O)(=O)- core. Identifies PFSA when combined with (A).
%[text]
%[text] Flagging Logic:
%[text] hasChain   = (A)             --> Presence of a full perfluoroalkyl chain
%[text] hasCF3     = (B)             --> Presence of CF3 group (may be isolated)
%[text] hasSulfonyl = (C)            --> Presence of sulfonic/sulfonamide head group
%[text]
%[text] isPFAS     = hasChain OR (hasCF3 AND NOT aliphatic)
%[text] classGuess: isPFAS AND hasSulfonyl -> "PFSA"
%[text] isPFAS AND NOT hasSulfonyl -> "PFCA_or_FTS"
%[text] else -> "NonPFAS"
%[text]
%[text] Note: This rule-based classifier is for illustration and not regulatory-grade.
%[text] Actual PFAS screening workflows (e.g., EPA DSSTox) use large SMARTS libraries with hundreds of patterns covering PFAS subclasses.
%[text]
%[text] SMARTS Patterns
logSection("A09", "Section 2: SMARTS-Based PFAS Flagging", "Analytics L3");
SMARTS_PFCHAIN  = "[#6](F)(F)-[#6](F)(F)";    % CF2-CF2 backbone
SMARTS_CF3      = "[#6](F)(F)F";              % -CF3 end
SMARTS_SULFONYL = "[#16](=O)(=O)";            % -SO2- head group
SMARTS_POLYFLUORO = "[#6](F)(F)[#6H2]";          % CF2 directly bonded to CH2 (FTS junction only)

validMols = mols(valid);
validNames = names(valid);
validSmiles = smiles(valid);
validCls = trueCls(valid);
nValid = sum(valid);

%[text] Apply substructure flags
hasChain    = emk.mol.hasSubstruct(validMols, SMARTS_PFCHAIN);
hasCF3      = emk.mol.hasSubstruct(validMols, SMARTS_CF3);
hasSulfonyl = emk.mol.hasSubstruct(validMols, SMARTS_SULFONYL);
hasPolyFluoro = emk.mol.hasSubstruct(validMols, SMARTS_POLYFLUORO);

%[text] Classification: isPFAS if perfluoroalkyl chain is present
%[text] Note: `hasCF3 & hasChain` is a subset of `hasChain`, so effectively `isPFAS = hasChain` (redundant terms are kept for clarity).
isPFAS = hasChain | (hasCF3 & hasChain);

%[text] Subclassify PFAS hits
classGuess = repmat("NonPFAS", nValid, 1);
for k = 1:nValid
    if isPFAS(k)
        if hasSulfonyl(k)
            classGuess(k) = "PFSA";
        elseif hasPolyFluoro(k) && ~hasSulfonyl(k)
            classGuess(k) = "FTS";
        else
            classGuess(k) = "PFCA";
        end
    end
end

%[text] Confusion matrix against true labels
nPFAS_hit = sum(isPFAS);
nPFAS_true = sum(validCls ~= "NonPFAS");
correct = sum(classGuess == validCls);

logInfo("PFAS Flagging Results:");
logInfo("  True PFAS: %d    Detected PFAS: %d    Correct Classification: %d / %d", ...
    nPFAS_true, nPFAS_hit, correct, nValid);

%[text] Display results table
flagTbl = table( ...
    validNames, validCls, classGuess, hasChain', hasCF3', hasSulfonyl', isPFAS', ...
    VariableNames=["Name","TrueClass","GuessedClass","HasChain","HasCF3","HasSulfonyl","IsPFAS"]);
disp(flagTbl);

%[text] **💡 Observation Point 2 — Limitations and Improvements of SMARTS Patterns**
%[text] Controls like halothane and sevoflurane also have -CF3 groups. Would the SMARTS_CF3 pattern match them?
%[text] Check the relevant rows for `hasCF3` and consider why `hasCF3` alone is insufficient for PFAS identification.
%[text] What structural features does SMARTS miss that distinguish FTS from PFCA? (Hint: FTS has at least one C-H bond in the fluorinated chain.)
%[text] Design a better SMARTS for FTS. FTS has a -CF2-CH2- junction.
%[text] Try `emk.mol.hasSubstruct(validMols, "[#6](F)(F)-[#6H2]")` to see if you can correctly isolate FTS compounds.
% ... (Try writing code here)
%%
%[text] ## Section 3: Calculating Persistence Proxies from Descriptors
%[text]
%[text] PFAS flagging is complete. Next, we will quantify environmental persistence using three physicochemical proxies (LogP: the logarithm of the partition coefficient between water and octanol, number of fluorines, TPSA: topological polar surface area).
%[text] Let's check which proxies each compound scores highly on.
%[text]
%[text] ### Concept: Three Proxies for Environmental Persistence
%[text]
%[text] Direct measurement of environmental half-life requires expensive laboratory assays, so regulators use physicochemical proxies as early warning indicators.
%[text] Three widely used proxies (Wang et al. 2017, Cousins et al. 2020):
%[text]
%[text] (1) Hydrophobicity (LogP)
%[text] High LogP indicates strong affinity for lipids, leading to bioaccumulation in fatty tissues.
%[text] The bioconcentration factor (BAF) often correlates with LogP for neutral organics (log BAF ~ 0.79 * LogP - 0.40, Veith 1979).
%[text] For PFAS, interpretation of LogP is complicated by ionic head groups (PFCA are acids at neutral pH). RDKit's Wildman-Crippen LogP is used as a structural proxy, not a true partition coefficient.
%[text]
%[text] Concern contribution: score_logP = clamp(LogP / 10, 0, 1)
%[text] (LogP = 10 maps to score 1.0, negative LogP maps to score 0)
%[text]
%[text] (2) Number of Fluorines (F Density)
%[text] Many C-F bonds increase chemical stability, slowing degradation.
%[text] C-F bond dissociation energy (~544 kJ/mol) far exceeds C-Cl (~397 kJ/mol) and C-C (~346 kJ/mol).
%[text] As a simple proxy, normalize the number of fluorine atoms (F count) to a range of 0-1.
%[text]
%[text] Concern contribution: score_F = clamp(nF / 17, 0, 1)
%[text] (17 fluorines are in PFDA, the longest compound in the set)
%[text]
%[text] (3) Water Solubility Proxy (TPSA)
%[text] Highly hydrophilic PFAS (high TPSA) are more likely to leach into water.
%[text] Highly hydrophobic ones (low TPSA) adsorb to sediments but accumulate in organisms.
%[text] Extreme values are of high concern. The bell-shaped concern peak is around TPSA ~ 40 Å², typical for neutral long-chain PFAS.
%[text]
%[text] Concern contribution: score_TPSA = exp(-((TPSA - 40)^2) / (2*30^2))
%[text] (Gaussian centered at TPSA=40 with σ=30 -- heuristic, for explanation)
logSection("A09", "Section 3: Calculating Persistence Proxies from Descriptors", "Analytics L3");
DESCS = ["LogP", "TPSA", "HeavyAtomCount"];
descTbl = emk.descriptor.batchCalculate(validMols, DESCS);

logp_vec    = descTbl.LogP;            % N x 1
tpsa_vec    = descTbl.TPSA;            % N x 1
hatom_vec   = descTbl.HeavyAtomCount;  % N x 1

%[text] Count fluorine atoms from SMILES (fast: no Python call needed).
nF_vec = zeros(nValid, 1);
for k = 1:nValid
    nF_vec(k) = count(validSmiles(k), "F");
end

%[text] Calculate the three persistence proxy scores (each 0-1).
score_logP = min(max(logp_vec / 10, 0), 1);
score_F    = min(max(nF_vec / 17,   0), 1);
score_TPSA = exp(-((tpsa_vec - 40).^2) / (2 * 30^2));

scoresMat = [score_logP, score_F, score_TPSA];  % N x 3

logInfo("Persistence Proxy Statistics (PFAS Hits Only):");
pfasIdx = find(isPFAS);
if ~isempty(pfasIdx)
    logInfo("  LogP Score    : Mean = %.3f  Max = %.3f", ...
        mean(score_logP(pfasIdx)), max(score_logP(pfasIdx)));
    logInfo("  F Count Score: Mean = %.3f  Max = %.3f", ...
        mean(score_F(pfasIdx)),    max(score_F(pfasIdx)));
    logInfo("  TPSA Score  : Mean = %.3f  Max = %.3f", ...
        mean(score_TPSA(pfasIdx)), max(score_TPSA(pfasIdx)));
end

%[text] Visualize the three proxy scores for all chemicals.
figure("Name","A09 Sec3: Persistence Proxy Scores");
X = 1:nValid;
bar(X, scoresMat, "grouped");
xticks(X);
xticklabels(validNames);
xtickangle(45);
ylabel("Score (0-1)");
title("Three Persistence Proxy Scores per Chemical");
legend(["LogP Proxy", "F Count Proxy", "TPSA Proxy"], Location="northeast");
grid on;

%[text] **💡 Observation Point 3 — Explore Persistence Proxies**
%[text] Check which compound has the highest F count score. Investigate `score_F` and `validNames`.
%[text] The longest PFAS chain should have the most fluorines, but do the results match expectations?
%[text] Consider why Ethanol scores near 0 for all three proxies. How many fluorine atoms does it have?
%[text] The TPSA proxy is a Gaussian centered at 40 Å². PFOS has a high TPSA ≈ 115 Å² due to its sulfonate head group.
%[text] Does the Gaussian unfairly underestimate PFOS compared to low-polarity PFCA? How would you redesign this proxy?
% ... (Try writing code here)
%%
%[text] ## Section 4: Weight Optimization Using fmincon
%[text]
%[text] The persistence proxy scores are ready. Next, we use fmincon to find the optimal weights for the three proxies.
%[text] We derive scientifically justified weight distribution using a least squares fit to a reference rubric.
%[text]
%[text] ### Concept: Constrained Weight Optimization for Multi-Criteria Scoring
%[text]
%[text] The final "concern score" for each chemical is a weighted sum of the three persistence proxies:
%[text]
%[text] score(k) = w1 * score_logP(k) + w2 * score_F(k) + w3 * score_TPSA(k)
%[text]
%[text] The weights (w1, w2, w3) represent the relative importance of each proxy according to regulatory guidance.
%[text] Experts may assign w2 > w1 because C-F bond density is a more direct measure of persistence than LogP.
%[text]
%[text] To find the "optimal" weights, we define a reference scoring rubric:
%[text] - Long-chain PFAS (PFOA, PFOS, PFNA, PFDA) -> Reference score >= 0.7
%[text] - Short-chain PFAS (PFBA, PFBS, PFHxA) -> Reference score < 0.5
%[text] - Non-PFAS controls -> Reference score <= 0.2
%[text]
%[text] We minimize the maximum squared deviation from the reference targets under the following constraints:
%[text] w1 + w2 + w3 = 1   (sum of weights is 1)
%[text] w1, w2, w3 >= 0.05  (each proxy contributes at least 5%)
%[text] w1, w2, w3 <= 0.80  (no single proxy dominates more than 80%)
%[text]
%[text] This is a box-constrained quadratic problem solved with fmincon.
%[text] Objective function: sum_k (score(k) - target(k))^2 (least squares fit)
%[text]
%[text] Reference target scores for each chemical
%[text] Long-chain PFAS: High concern (target >= 0.7)
%[text] Short-chain PFAS: Medium concern (target < 0.5)
%[text] Non-PFAS: Low concern (target <= 0.2)
targetScores = zeros(nValid, 1);
for k = 1:nValid
    cls = validCls(k);
    nm  = validNames(k);
    if cls == "NonPFAS"
        targetScores(k) = 0.10;
    elseif any(nm == ["PFOA","PFOS","PFNA","PFDA"])
        targetScores(k) = 0.80;   % Long-chain, high concern
    elseif any(nm == ["PFHxA","PFHxS","8:2FTS"])
        targetScores(k) = 0.60;   % Mid-chain, medium to high concern
    elseif any(nm == ["PFBA","PFBS","6:2FTS","4:2FTS"])
        targetScores(k) = 0.40;   % Short-chain, medium concern
    else
        targetScores(k) = 0.50;   % Default for unassigned PFAS
    end
end

%[text] Objective function: Least squares fit to target scores
%[text] f(w) = sum_k (scoresMat * w - targetScores)^2
logSection("A09", "Section 4: Weight Optimization Using fmincon", "Analytics L3");
objFun = @(w) sum((scoresMat * w - targetScores).^2);

%[text] Constraints
%[text] Aeq * w = beq  -->  w1 + w2 + w3 = 1
Aeq  = [1, 1, 1];
beq  = 1;
lb   = [0.05; 0.05; 0.05];   % Lower bounds
ub   = [0.80; 0.80; 0.80];   % Upper bounds
w0   = [1/3; 1/3; 1/3];       % Initial guess: equal weights

if hasOptTbx
    opts = optimoptions("fmincon", Display="off", Algorithm="sqp");
    [w_opt, fval] = fmincon(objFun, w0, [], [], Aeq, beq, lb, ub, [], opts);
    logInfo("fmincon optimized weights:  w_logP=%.3f  w_F=%.3f  w_TPSA=%.3f  (RSS=%.4f)", ...
        w_opt(1), w_opt(2), w_opt(3), fval);
else
    w_opt = w0;   % equal-weight fallback
    logWarn("Using equal weights [1/3, 1/3, 1/3] (Optimization Toolbox not available).");
end

%[text] Calculate final concern scores with optimized weights.
concernScore = scoresMat * w_opt;

%[text] Display top concerns.
[sortedScores, sortIdx] = sort(concernScore, "descend");
logInfo("Ranked concern table (Top 10):");
for k = 1:min(10, nValid)
    i = sortIdx(k);
    logInfo("  %2d. %-18s  [%s]  Score = %.3f", k, validNames(i), validCls(i), sortedScores(k));
end

%[text] Bar graph: Ranked concern scores
figure("Name","A09 Sec4: Ranked Chemical Concern Scores");
barColors = zeros(nValid, 3);
for k = 1:nValid
    i = sortIdx(k);
    switch validCls(i)
        case "PFCA",    barColors(k,:) = [0.85, 0.20, 0.10];   % red
        case "PFSA",    barColors(k,:) = [0.95, 0.55, 0.10];   % orange
        case "FTS",     barColors(k,:) = [0.25, 0.55, 0.85];   % blue
        case "NonPFAS", barColors(k,:) = [0.70, 0.70, 0.70];   % grey
    end
end
barH = bar(sortedScores, "FaceColor","flat", "HandleVisibility","off");
barH.CData = barColors;
xticks(1:nValid);
xticklabels(validNames(sortIdx));
xtickangle(55);
ylabel("Concern Score (0-1)");
title(sprintf("Environmental Concern Scores (w=[%.2f, %.2f, %.2f])", ...
    w_opt(1), w_opt(2), w_opt(3)));
%[text] Manual legend patches
patch(NaN,NaN,[0.85 0.20 0.10], DisplayName="PFCA");
patch(NaN,NaN,[0.95 0.55 0.10], DisplayName="PFSA");
patch(NaN,NaN,[0.25 0.55 0.85], DisplayName="FTS");
patch(NaN,NaN,[0.70 0.70 0.70], DisplayName="NonPFAS");
legend(Location="northeast");
grid on;
yline(0.5, "--k", "Concern Threshold", LabelHorizontalAlignment="left", HandleVisibility="off");

%[text] **💡 Observation Point 4 — Let's review the results of weight optimization**
%[text] Among w_logP, w_F, and w_TPSA, which proxy received the highest weight?
%[text] Does it align with scientific literature? (Buck et al. 2011: C-F bond count suggested as a major factor)
%[text] If you rerun with the PFOS target score changed from 0.80 to 0.95, how does `w_opt` change?
%[text] Does the weight for F count increase? (PFOS has 17 F atoms.)
%[text] What happens if you add a fourth proxy (-CF2- unit count = chain length)?
% ... (Try writing the code here)
%%
%[text] ## Section 5: Structural Similarity Heatmap of PFAS Hits
%[text]
%[text] The calculation of concern scores is complete. Next, we will visualize the structural similarity between PFAS hits using a heatmap.
%[text] Let's see how similar PFAS of the same class are.
%[text]
%[text] ### Concept: Tanimoto Similarity and Structural Clustering
%[text]
%[text] The Tanimoto coefficient between two molecules is calculated from Morgan (circular) fingerprints.
%[text]
%[text] T(A, B) = |A AND B| / |A OR B|
%[text]
%[text] Here, |A AND B| represents the number of bits set in both fingerprints, and |A OR B| represents the number of bits set in either.
%[text]
%[text] T = 1.0 indicates identical fingerprints (structurally indistinguishable).
%[text] T = 0.0 indicates no shared bits (completely dissimilar).
%[text]
%[text] In PFAS, long-chain homologs (PFOA and PFNA) differ only in the number of -CF2- repeat units.
%[text] The fingerprints are highly similar (T ~ 0.7–0.9), sharing all substructures except for the additional -CF2- repeats.
%[text]
%[text] In contrast, PFCA and PFSA have different head groups: -COOH and -SO3H.
%[text] The Tanimoto similarity is moderate (T ~ 0.3–0.6).
%[text]
%[text] Hierarchical clustering dendrograms (average linkage / UPGMA) group the most similar structures, revealing subfamilies of PFCA / PFSA / FTS.
%[text] Average linkage is used instead of Ward because Tanimoto distance is non-Euclidean. Ward linkage assumes Euclidean space (Willett 1998).
logSection("A09", "Section 5: Structural Similarity Heatmap of PFAS Hits", "Analytics L3");
pfasHitIdx = find(isPFAS);
nHits      = numel(pfasHitIdx);

if nHits > 1
    % Calculate Morgan fingerprints for PFAS hits
    fps = cell(1, nHits);
    for k = 1:nHits
        fps{k} = emk.fingerprint.morgan(validMols{pfasHitIdx(k)});
    end

    % Construct Tanimoto similarity matrix
    simMat = zeros(nHits, nHits);
    for i = 1:nHits
        for j = 1:nHits
            simMat(i,j) = emk.similarity.tanimoto(fps{i}, fps{j});
        end
    end

    hitNames = validNames(pfasHitIdx);

    % Visualize as a heatmap
    figure("Name","A09 Sec5: PFAS Tanimoto Similarity Heatmap");
    imagesc(simMat);
    colormap(hot);
    colorbar();
    clim([0 1]);
    xticks(1:nHits); xticklabels(hitNames); xtickangle(45);
    yticks(1:nHits); yticklabels(hitNames);
    title("Morgan Fingerprint Tanimoto Similarity (PFAS Hits)");
    axis square;

    % Hierarchical clustering (Statistics and ML Toolbox)
    % Use average linkage (UPGMA) because Tanimoto distance is non-Euclidean;
    % Ward linkage requires Euclidean distances and generates false warnings.
    % Average linkage is the standard choice for chemoinformatics fingerprint clustering (Willett 1998).
    if hasStatsTbx
        distVec  = squareform(1 - simMat);   % condensed distance vector
        Z        = linkage(distVec, "average");
        figure("Name","A09 Sec5: PFAS Average Linkage Dendrogram");
        dendrogram(Z, "Labels", hitNames, Orientation="left");
        title("Average Linkage Dendrogram -- PFAS Structural Families");
        xlabel("Distance (1 - Tanimoto)");
    else
        logWarn("Statistics and ML Toolbox not available -- skipping dendrogram.");
    end

    logInfo("Tanimoto Statistics for PFAS Hits:");
    upperTri = simMat(triu(true(nHits), 1));
    logInfo("  Mean = %.3f   Min = %.3f   Max = %.3f", ...
        mean(upperTri), min(upperTri), max(upperTri));
else
    logWarn("Fewer than 2 PFAS hits; skipping heatmap.");
end

%[text] **💡 Observation Point 5 — Delve into Tanimoto Similarity**
%[text] Identify the pair of PFAS compounds with the highest Tanimoto similarity.
%[text] Check the Tanimoto similarity between PFOA (PFCA, C8) and PFOS (PFSA, C8). They have the same carbon chain length but different head groups.
%[text] Determine whether T is closer to 0.5 (moderate) or 0.9 (high).
%[text] Consider why PFOS is placed away from PFBS and PFHxS (other PFSAs) and towards FTS in the dendrogram.
%[text] (Hint: Both PFOS and 8:2FTS have F atom count = 17. Morgan FP may emphasize long-chain similarity over head group differences.)
%[text] Try recalculating with the Morgan fingerprint radius set to `Radius=1`.
%[text] `fps{k} = emk.fingerprint.morgan(validMols{pfasHitIdx(k)}, Radius=1);`
%[text] Consider whether reducing the radius increases or decreases similarity between PFCA and PFSA families. Why?
% ... (Try writing code here)
%%
%[text] ## Section 6: Summary Report and Screening Dashboard
%[text]
%[text] Structural clustering is complete. Finally, we will compile all analysis results into a single screening report.
%[text] We will create a ranked table of concerns and scatter plots for decision-makers.
%[text]
%[text] ### Concept: Communicating Screening Results
logSection("A09", "Section 6: Summary Report and Screening Dashboard", "Analytics L3");
%[text] Construct a summary table for all chemicals.
[~, rankOrder] = sort(concernScore, "descend");
rankVec        = zeros(nValid, 1);
rankVec(rankOrder) = (1:nValid)';

summaryTbl = table( ...
    rankVec, ...
    validNames, ...
    validCls, ...
    classGuess, ...
    round(logp_vec,2), ...
    nF_vec, ...
    round(tpsa_vec,1), ...
    round(concernScore,3), ...
    VariableNames=["Rank","Name","TrueClass","FlaggedClass", ...
                   "LogP","nF","TPSA_A2","ConcernScore"]);
summaryTbl = sortrows(summaryTbl, "Rank");

logInfo("=== PFAS Screening Report ===");
disp(summaryTbl);

%[text] Scatter plot: LogP vs F number colored by concern score
figure("Name","A09 Sec6: LogP vs F Number (Concern Score)");
scatter(logp_vec, nF_vec, 80, concernScore, "filled", MarkerFaceAlpha=0.8);
colormap(turbo);
cb = colorbar();
cb.Label.String = "Concern Score";
clim([0 1]);
for k = 1:nValid
    text(logp_vec(k) + 0.1, nF_vec(k), validNames(k), FontSize=7);
end
xlabel("LogP (Wildman-Crippen)");
ylabel("Number of Fluorine Atoms");
title("Chemical Space View: LogP vs F Number (Color = Concern Score)");
grid on;

%[text] **Summary**
%[text] We learned about automatic PFAS detection using SMARTS patterns, persistence scoring using physicochemical proxies, constrained weight optimization with fmincon, Tanimoto similarity heatmaps, and hierarchical clustering.
%[text] We built a dashboard that meets regulatory screening requirements (ranking, justification, refutability).
logInfo("A09: Complete.");

%[text] **💡 Observation Point 6 — Let's Interpret the Final Screening Results**
%[text] Check for visual clusters of high-concern PFAS in the upper right region of the scatter plot (high LogP and high F number).
%[text] For PFAS compounds, determine which single screening criterion results in the fewest false negatives: LogP threshold or F number threshold.
%[text] Sevoflurane and Halothane have C-F bonds but are volatile anesthetics, differing in environmental persistence from PFAS.
%[text] Consider how to modify the scoring rubric to flag anesthetic gases separately from true PFAS.
%[text] (Consider adding a volatility proxy: MW < 200 g/mol as a "volatility exemption" to lower the final score.)
%[text] Export the summary table to CSV: `writetable(summaryTbl, "pfas_screening_report.csv")`
%[text] Open it in Excel. Consider how to present it to non-expert stakeholders, which columns to keep or remove.
% ... (Try writing the code here)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]%   data: {"layout":"inline","rightPanelPercent":40}
%---
