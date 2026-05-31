%[text] # S06: Search Compounds in PubChem
%[text] EasyMolKit Application Story — Layer 2
%[text] 
%[text] PubChem, the world's largest free chemical database, has over 100 million compounds registered. Checking each one via a browser is inefficient, but by querying directly from MATLAB, you can obtain and compare data for dozens of compounds in seconds.
%[text] In this script, we will retrieve data for 5 antipyretic analgesics via the PubChem PUG REST API and combine it with RDKit to construct a property comparison table.
%[text] ### Story
%[text] You are a researcher at a small pharmaceutical company.
%[text] Your team has just launched a new anti-inflammatory drug project and needs to quickly understand the profiles of existing competing antipyretic analgesics (Aspirin, Ibuprofen, Acetaminophen, Naproxen, Celecoxib).
%[text] Normally, this would involve manually browsing the PubChem website one by one, but here we will use the `emk.db.searchPubchem` function to directly hit the API from MATLAB and obtain and compare structure and property data in bulk.
%[text] **Exercise Flow**
%[text] 1. Search for a single compound by common name
%[text] 2. Search for a compound from a SMILES string (when the structure is known but the name is not)
%[text] 3. Query a list of compounds and create a comparison table
%[text] 4. Visualize structures obtained from PubChem
%[text] 5. Sort and filter the table to find candidates with the best profiles \
%[text] ### Learning Objectives
%[text] - Query PubChem via PUG REST API using `emk.db.searchPubchem`
%[text] - Utilize three query modes: name, SMILES, CID
%[text] - Construct a comparison table for multiple compounds from API results
%[text] - Combine obtained properties with Lipinski Ro5 filter
%[text] - Visualize structures from obtained SMILES strings \
%[text] ### Prerequisites
%[text] - Completion of F01 (Molecule Drawing) and F02 (Property Calculation)
%[text] - **Internet connection required** (PubChem PUG REST API)
%[text] - RDKit installed (execute `emk.setup.install()` once)
%[text] - No additional Toolbox required (works with MATLAB alone) \
%[text] **Estimated Time**: 25–40 minutes | Execution Method: Run each section with Ctrl+Enter
%[text] 
%[text] **References**: 
%[text] - Kim S et al. (2023) PubChem 2023 update. *Nucleic Acids Res* 51:D1373-D1380.
%[text] - doi:10.1093/nar/gkac956 — PubChem PUG REST: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest [Open Access]\
%[text] **Note**: This story sends queries to the live PubChem API. An internet connection is required. Results may vary slightly if PubChem updates its database.
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
%[text] Warm up Python/RDKit process (the first call may take some time).
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;
logInfo("S06: Setup complete");
%%
%[text] ## Section 1: Search for a Single Compound by Common Name
%[text] ### PubChem PUG REST API
%[text] PubChem is the world's largest free chemical information database managed by the National Institutes of Health (NIH). It serves as an infrastructure for drug discovery research with a vast amount of data:
%[text] - Compound Database: Approximately 115 million unique compound structures (standardized)
%[text] - Substance Database: Approximately 300 million unrefined registered structures (including source information)
%[text] - BioAssay Database: Approximately 270 million bioactivity test data \
%[text] The mechanism for direct programmatic access to this database is the "PUG REST API". API requests (URLs) are basically constructed according to the following rules (web patterns):
%[text] [https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/\<input\_type\>/\<query\>/property/\<property\_list\>/JSON](https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/%3Cinput_type%3E/%3Cquery%3E/property/%3Cproperty_list%3E/JSON)
%[text] [Meaning of Each Parameter]
%[text] - \<input\_type\>  : Type of search key (name (common name), smiles (structural formula), cid (compound ID), etc.)
%[text] - \<query\>       : Specific search keyword (e.g., "aspirin" or structural string)
%[text] -  \<property\_list\> : List of desired property items (comma-separated) \
%[text] [Specific Examples of API Requests]
%[text] - Retrieve properties by common name:   /compound/name/aspirin/property/MolecularWeight,IsomericSMILES/JSON
%[text] - Retrieve molecular formula from SMILES: /compound/smiles/CCO/property/MolecularFormula,MolecularWeight/JSON
%[text] - Retrieve IUPAC name from compound ID: /compound/cid/2244/property/IUPACName,MolecularFormula/JSON \
%[text] The `emk.db.searchPubchem` function automatically handles this complex URL generation, HTTP communication, and response (JSON) parsing internally, and returns it wrapped in a MATLAB "table type" for easy handling. Let's first search for "Aspirin (Acetylsalicylic Acid)", a classic antipyretic analgesic discovered in 1897, by its common name and check the API behavior.
logInfo("Querying PubChem for Aspirin by common name...");
tbl_aspirin = emk.db.searchPubchem("aspirin");

logInfo("PubChem results for Aspirin:");
logInfo("  CID              : %d",   tbl_aspirin.CID(1));
logInfo("  IUPAC Name       : %s",   tbl_aspirin.IUPACName(1));
logInfo("  Molecular Formula: %s",   tbl_aspirin.MolecularFormula(1));
logInfo("  Molecular Weight : %.2f g/mol", tbl_aspirin.MolecularWeight(1));
logInfo("  SMILES           : %s",   tbl_aspirin.IsomericSMILES(1));
%[text] **✏️ Try It 1**
%[text] Q: What is the PubChem CID for Aspirin?
%[text]    CID (Compound ID) is an invariant number that uniquely identifies a compound in the database.
%[text]    Hint: Check the value of `tbl\_aspirin.CID(1)`.
%[text]    Expected value: `2244`
%[text] Q: The IUPAC name returned by PubChem is much more systematic than "aspirin".
%[text]    What does it reveal about the structure?
%[text]    Expected value: `2-acetyloxybenzoic acid` (The name indicates that an acetyloxy group is attached at the 2-position (ortho position) of benzoic acid).
%[text] Q: Try searching for "paracetamol" -- this is the British name for "acetaminophen".
%[text]    Do both queries return the same CID?
%[text]    Even if called by different common names worldwide (UK: Paracetamol, US/JP: Acetaminophen), if they refer to the same compound, PubChem returns the same CID.  
%[text]    Hint: Execute `emk.db.searchPubchem("paracetamol")` and `emk.db.searchPubchem("acetaminophen")` in the command window and compare the CIDs.
%[text]    Expected value: Both return CID `1983`.
% ... (Try writing code here)
%%
%[text] ## Section 2: Search with SMILES -- When the structure is known
%[text] ### Structure-based search
%[text] In actual drug discovery workflows, there are cases where you have "structures drawn with ChemDraw" or "new structures designed with computational chemistry," but you don't know their exact chemical names, or you want to check if they are registered as known compounds in a database.
%[text] By querying PubChem with the "SMILES string," which is a structural representation of the compound, you can gain the following powerful advantages:
%[text] 1. Quick determination of whether the structure is "known (registered)"
%[text] 2. Identification of unique identifiers (CID) and official IUPAC names registered in the database
%[text] 3. Acquisition of highly reliable basic physical property data verified and standardized by public institutions \
%[text] **Note**: The SMILES sent as a query must be a valid structure that can be correctly parsed by RDKit, etc. PubChem performs its own "canonicalization" internally, so even if the input SMILES order is slightly different, it will still match correctly.
%[text] Here, let's try a reverse search using the structure (SMILES) of "Ibuprofen (R/S racemate)," which is widely used as a commercial antipyretic analgesic.
IBU_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O";

logInfo("Querying PubChem for Ibuprofen with SMILES...");
tbl_ibu = emk.db.searchPubchem(IBU_SMILES, Type="smiles");

logInfo("PubChem results for Ibuprofen (Type=smiles):");
logInfo("  CID              : %d",   tbl_ibu.CID(1));
logInfo("  IUPAC Name       : %s",   tbl_ibu.IUPACName(1));
logInfo("  Molecular Formula: %s",   tbl_ibu.MolecularFormula(1));
logInfo("  Molecular Weight : %.2f g/mol", tbl_ibu.MolecularWeight(1));
logInfo("  SMILES (PubChem) : %s",   tbl_ibu.IsomericSMILES(1));

%[text] **✏️ Try It 2**
%[text] Q: Does PubChem return the same SMILES as input,
%[text]    or does it return a different (canonicalized) form?
%[text]    SMILES are not unique -- there are many valid SMILES strings representing the same molecule.
%[text]    PubChem always returns its canonical form.
%[text] Q: What is the PubChem CID for Ibuprofen? Expected value: 3672.
%[text] Q: Let's search for the pharmacologically active (S)-Ibuprofen.
%[text]    Does it have a different CID?
%[text]    (S)-Ibuprofen SMILES: "[C@@H](C(=O)O)(Cc1ccc(cc1)CC(C)C)C"
%[text]    Hint: Pass the SMILES `"[C@@H](C(=O)O)(Cc1ccc(cc1)CC(C)C)C"` and `Type="smiles"` to `emk.db.searchPubchem`.
%[text]    Note: The (S) form is registered in PubChem as a structure with different stereochemistry,
%[text]    so it has its own unique CID.
% ... (Try writing the code here)
%%
%[text] ## Section 3: Create a Comparison Table for Antipyretic Analgesic Panel
%[text] ### Drug Panel Comparison
%[text] Comparing a "set of similar and effective drugs (drug panel)" targeting a specific disease or a "series of synthetic compounds (chemical series)" under study side by side is one of the most fundamental tasks in medicinal chemistry. To decipher how slight structural differences affect properties and activity (Structure-Activity Relationship: SAR), it is first necessary to compile the data into a uniform table.
%[text] In this section, we will treat "five famous pharmaceuticals" representing the evolution of antipyretic analgesic and anti-inflammatory drugs over more than a century as a panel and obtain data at once:
%[text] - Aspirin         : Classical non-selective COX inhibitor (introduced in 1897)
%[text] - Ibuprofen       : Representative non-steroidal anti-inflammatory drug (NSAID / introduced in 1961)
%[text] - Acetaminophen   : Antipyretic analgesic mainly acting on the central nervous system (introduced in 1956)
%[text] - Naproxen        : Long-acting non-selective NSAID (introduced in 1976)
%[text] - Celecoxib       : Selective inhibitor targeting COX-2 to reduce gastrointestinal side effects (introduced in 1998) \
%[text] Let's use MATLAB's loop processing to extract and compare how the "evolutionary history" of molecules, from simple small molecules like aspirin to modern molecules optimized for increased target selectivity like celecoxib, is reflected in property data.
PANEL_NAMES = ["aspirin"; "ibuprofen"; "acetaminophen"; "naproxen"; "celecoxib"];
N = numel(PANEL_NAMES);

logInfo("%d compounds being queried in PubChem...", N);

%[text] Pre-allocate result arrays
cids     = zeros(N, 1, 'uint32');
iupac    = strings(N, 1);
formulas = strings(N, 1);
mws      = zeros(N, 1);
panelSmiles = strings(N, 1);

for i = 1:N
    name = PANEL_NAMES(i);
    try
        r = emk.db.searchPubchem(name);
        cids(i)     = r.CID(1);
        iupac(i)    = r.IUPACName(1);
        formulas(i) = r.MolecularFormula(1);
        mws(i)      = r.MolecularWeight(1);
        panelSmiles(i) = r.IsomericSMILES(1);
        logInfo("  [%d/%d] %-15s -> CID %d, MW=%.1f", i, N, name, cids(i), mws(i));
    catch ME
        logWarn("  [%d/%d] %s: Query failed -- %s", i, N, name, ME.message);
    end
end

%[text] Assemble the comparison table
panelTbl = table(PANEL_NAMES, cids, iupac, formulas, mws, panelSmiles, ...
    'VariableNames', {'Name', 'CID', 'IUPACName', 'MolecularFormula', ...
                      'MolecularWeight', 'SMILES'});

logInfo("--- Antipyretic Analgesic Panel (from PubChem) ---");
disp(panelTbl(:, {'Name', 'CID', 'MolecularFormula', 'MolecularWeight'}));

%[text] **✏️ Try It 3**
%[text] Q: Which compound in the panel has the highest molecular weight?
%[text]    Hint: panelTbl(panelTbl.MolecularWeight == max(panelTbl.MolecularWeight), :)
%[text]    Expected: Celecoxib ($&dollar&;MW \\approx 381~\\text{Da}&dollar&;$).
%[text]    Reflects the trend towards larger, more selective molecules in modern drug discovery.
%[text] Q: Both aspirin and acetaminophen are small molecules ($&dollar&;MW \< 180~\\text{Da}&dollar&;$).
%[text]    Despite similar sizes, their mechanisms of action are completely different.
%[text]    What does this indicate about the relationship between MW and activity?
%[text]    (Hint: Size alone does not determine mechanism; functional groups do.)
%[text] Q: Looking at the molecular formulas, which compound contains sulfur (S)?
%[text]    Hint: panelTbl(contains(panelTbl.MolecularFormula, "S"), :)
%[text]    Expected: Celecoxib (C17H14F3N3O2S -- the sulfonamide group is key to COX-2 selectivity).
% ... (Try writing code here)
%%
%[text] ## Section 4: Calculate RDKit Descriptors and Integrate with PubChem Data
%[text] ### Combining API Data and Local Calculations
%[text] External web APIs (such as PubChem) are very useful, but they may not cover all the advanced compound descriptors needed in drug discovery research (e.g., detailed lipophilicity indices or specific 3D descriptors). A commonly used workflow in practice is a hybrid approach: "Obtain basic information from the API and calculate more specialized descriptors locally to combine them."
%[text] [Hybrid Data Construction Flow]
%[text] 1. Obtain the correct "structural information (standardized SMILES)" and "basic information" of compounds from the PubChem API
%[text] 2. Input the obtained SMILES into a local cheminformatics engine (RDKit)
%[text] 3. Use `emk.descriptor.calculate` to quickly calculate descriptors (LogP, TPSA, various atom and ring counts, etc.) that were not available from the API locally
%[text] 4. Use MATLAB's powerful table join function to integrate all data into a single refined data frame (matrix table) \
%[text] Through this section, you will learn the method to construct a "complete compound property profile table" that can be directly applied to data science and machine learning inputs.
PROP_NAMES = ["LogP", "TPSA", "NumHDonors", "NumHAcceptors", ...
              "NumRotatableBonds", "RingCount"];

nProps   = numel(PROP_NAMES);
descMat  = nan(N, nProps);     % Rows: compounds, Columns: descriptors
validIdx = false(N, 1);

logInfo("Calculating RDKit descriptors from PubChem SMILES...");

for i = 1:N
    smi = panelSmiles(i);
    if strlength(smi) == 0
        logWarn("  Skipping %s: Could not obtain SMILES from PubChem", PANEL_NAMES(i));
        continue;
    end
    if ~emk.mol.isValid(smi)
        logWarn("  Skipping %s: Cannot parse SMILES with RDKit", PANEL_NAMES(i));
        continue;
    end
    mol = emk.mol.fromSmiles(smi);
    d   = emk.descriptor.calculate(mol, PROP_NAMES);
    for j = 1:nProps
        descMat(i, j) = d.(PROP_NAMES(j));
    end
    validIdx(i) = true;
end
%[text] Construct descriptor table
descTbl = array2table(descMat, 'VariableNames', cellstr(PROP_NAMES));
%[text] Combine PubChem table and descriptor table
fullTbl = [panelTbl(:, {'Name', 'CID', 'MolecularFormula', 'MolecularWeight'}), ...
           descTbl];

logInfo("--- Complete Property Table (PubChem + RDKit) ---");
disp(fullTbl);
%[text] **✏️ Try It 4**
%[text] Q: Verify that all 5 compounds pass Lipinski's Rule of Five ($&dollar&;MW \< 500&dollar&;$,
%[text]    $&dollar&;\\text{LogP} \\le 5&dollar&;$, $&dollar&;\\text{HBD} \\le 5&dollar&;$, $&dollar&;\\text{HBA} \\le 10&dollar&;$) using
%[text]  `emk.filter.lipinski`.
%[text]    Hint:
%[text]      roTbl = table(fullTbl.MolecularWeight, fullTbl.LogP, ...
%[text]                    fullTbl.NumHDonors,      fullTbl.NumHAcceptors, ...
%[text]                    'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
%[text]      result = `emk.filter.lipinski(roTbl)`;
%[text]      disp(result)
%[text]    Expected: All 5 compounds have Pass\_Ro5 = true.
%[text] Q: Which compound has the highest TPSA?
%[text]    TPSA (Topological Polar Surface Area) correlates with membrane permeability:
%[text]    $&dollar&;\\text{TPSA} \< 90~\\text{A}^2&dollar&;$ predicts good oral absorption.
%[text]    Expected: Celecoxib (sulfonamide + pyrazole with highest TPSA).
%[text]    Still, Celecoxib's TPSA is below $&dollar&;90~\\text{A}^2&dollar&;$ ($&dollar&;77.98~\\text{A}^2&dollar&;$), indicating no issues with oral absorption.
%[text] Q: Which compound has the lowest LogP (most hydrophilic)?
%[text]    Hydrophilicity affects solubility and bioavailability.
%[text]    Expected: Aspirin (RDKit $&dollar&;\\text{LogP} \\approx 1.31&dollar&;$) or Acetaminophen (RDKit $&dollar&;\\text{LogP} \\approx 1.35&dollar&;$)
%[text]    — Both are nearly equivalent and the most hydrophilic in this panel.
%[text]    **Note**: Experimentally, Acetaminophen ($&dollar&;\\approx 0.46&dollar&;$) is the lowest, but
%[text]      RDKit Crippen model may not match experimental values.
% ... (Try writing code here)
%%
%[text] ## Section 5: Visualizing Panel Structures
%[text] ### Visualization of Structures and Properties
%[text] In medicinal chemistry, it is crucial to not only look at property data as numbers (such as molecular weight and LogP) but also to visually observe which parts of the molecule (functional groups or skeleton) these numbers originate from by aligning 2D structural formulas.
%[text] When visualizing the current panel of antipyretic analgesics, let's focus on the link between the following "chemical structural features" and "actual pharmacological properties and side effects":
%[text] - Presence of carboxylic acid group (-COOH): Aspirin, Ibuprofen, and Naproxen all share the carboxylic acid. This indicates acidity and is the main cause of the stomach pain and gastrointestinal disturbances (physical and physiological irritation to the gastric mucosa) typical of NSAIDs.
%[text] - Adoption of amide group (-NHCO-): Acetaminophen does not have a carboxylic acid but has an amide structure. This structural difference is directly linked to its clinical property of being gentle on the stomach (less likely to cause gastrointestinal disturbances).
%[text] - Distinctive functional groups: Visually confirm how the "sulfonamide group (-SO2NH2)" and "trifluoromethyl group (-CF3)" introduced in Celecoxib change the molecular weight and lipophilicity. \
%[text] Pass the high-precision SMILES data obtained from PubChem to `emk.viz.draw2d` and utilize MATLAB's graphics functions (subplot) to draw a neat structural formula grid.
logInfo("Drawing panel structures (from PubChem SMILES)...");

%[text] Draw all compounds as a grid in one figure (same pattern as F01 Section 5)
nCols = ceil(sqrt(N));
nRows = ceil(N / nCols);
figure("Name", "Antipyretic Analgesic Panel -- Structures", "Position", [100 500 1100 560]);
for i = 1:N
    smi = panelSmiles(i);
    subplot(nRows, nCols, i);
    if strlength(smi) == 0 || ~emk.mol.isValid(smi)
        text(0.5, 0.5, "N/A", HorizontalAlignment="center");
        axis("off");
        continue;
    end
    mol = emk.mol.fromSmiles(smi);
    titleStr = sprintf("%s\n[CID %d] MW=%.0f", PANEL_NAMES(i), cids(i), mws(i));
    emk.viz.draw2d(mol, Title=titleStr, Width=320, Height=300);
end

%[text] **✏️ Try It 5**
%[text] Observe the 5 structures drawn above. Identify each of them:
%[text]   (a) Acidic or amide functional groups (related to mechanism and gastrointestinal effects)
%[text]   (b) Ring systems (1 ring vs. 2 rings vs. none)
%[text]   (c) Heteroatoms in the ring (Aspirin/Naproxen: benzene only;
%[text]       Celecoxib: pyrazole ring with 2 nitrogens)
%[text] Q: Celecoxib has a trifluoromethyl group (CF3) and a pyrazole ring.
%[text]    Both are common "bioisosteres" in medicinal chemistry.
%[text]    Why introduce a CF3 group?
%[text]    (Hint: CF3 increases lipophilicity, improves metabolic stability,
%[text]     can mimic a methyl group while inhibiting oxidative metabolism at that site.)
% ... (Try writing code here)
%%
%[text] ## Section 6: Rank and Filter Panels
%[text] ### Property-Based Prioritization
%[text] When narrowing down (or prioritizing) lead compounds from hit compounds or existing drug candidates for the desired treatment, researchers use medicinal chemistry screening criteria to filter and sort (rank) the candidates.
%[text] [Examples of Common Filtering Criteria]
%[text] 1. Evaluation of oral drug-likeness (Lipinski's Rule of Five: Ro5)  \
%[text]              \-Molecular weight (MW) $&dollar&;\< 500&dollar&;$, lipophilicity (LogP) $&dollar&;\\le 5&dollar&;$, number of hydrogen bond donors (HBD) $&dollar&;\\le 5&dollar&;$, number of hydrogen bond acceptors (HBA) $&dollar&;\\le 10&dollar&;$
%[text] 2\. CNS permeability (ease of crossing the blood-brain barrier (BBB))
%[text]             \-Topological polar surface area (TPSA) $&dollar&;\< 90~\\text{A}^2&dollar&;$ (ideally $&dollar&;60~\\text{A}^2&dollar&;$ or less is more ideal)
%[text]             \-Molecular weight (MW) $&dollar&;\< 450&dollar&;$+    3. Avoidance of side effect risks (e.g., reduction of gastrointestinal irritation)
%[text] 3\. Avoidance or masking of strongly acidic functional groups like carboxylic acids
%[text] In this section, we will sort and graph panel compounds predicted to have high CNS permeability (easy to reach the central nervous system) based on locally calculated "TPSA (polar surface area)" values using RDKit. Learn the basics of a logical decision-making process based on property data.
logInfo("Ranking compounds in ascending order of TPSA (ascending = best CNS permeability):");

[~, sortOrder] = sort(fullTbl.TPSA, "ascend");
ranked = fullTbl(sortOrder, :);
ranked.Properties.RowNames = {};

for i = 1:height(ranked)
    cnsFlag = "";
    if ranked.TPSA(i) < 60 && ranked.LogP(i) >= 1 && ranked.LogP(i) <= 3
        cnsFlag = " [Good CNS permeability: TPSA<60 + LogP 1-3]";
    elseif ranked.TPSA(i) < 60
        cnsFlag = " [CNS permeability: Low TPSA]";
    elseif ranked.TPSA(i) < 90
        cnsFlag = " [CNS permeability present]";
    else
        cnsFlag = " [Difficult CNS access]";
    end
    logInfo("  %d. %-15s  TPSA=%.1f A^2  LogP=%.2f%s", ...
        i, ranked.Name(i), ranked.TPSA(i), ranked.LogP(i), cnsFlag);
end

%[text] Bar Graph: TPSA Comparison
figure("Name", "Antipyretic Analgesic Panel -- TPSA Comparison", "Position", [100 100 560 380]);
bar(fullTbl.TPSA(sortOrder));
xticks(1:N);
xticklabels(ranked.Name);
xtickangle(20);
ylabel("TPSA (A^2)");
title("Antipyretic Analgesic Panel -- Topological Polar Surface Area");
yline(90,  "--r", "TPSA=90 (Oral Absorption Threshold)", "LabelHorizontalAlignment", "left");
yline(60,  ":b", "TPSA=60 (CNS Recommended)", "LabelHorizontalAlignment", "left");
grid("on");

%[text] **✏️ Try It 6** — Reading Property Profiles
%[text] Q: Explain why TPSA (topological polar surface area) is used as an indicator of blood-brain barrier (BBB) permeability.
%[text]    Why do highly polar molecules find it difficult to pass through lipid bilayers (cell membranes)?
%[text]    Also, how does the lipophilicity indicated by LogP contribute to CNS permeability?
%[text]    (Hint: CNS permeability requires passing through the phospholipid bilayer.
%[text]    The larger the polar surface area, the higher the water solubility, making membrane permeability difficult.)
%[text] Q: Celecoxib has the highest TPSA in the panel.
%[text]    Do you think it can easily pass through the blood-brain barrier?
%[text]    Why is this acceptable for its intended use (arthritis inflammation)?
%[text]    (Hint: CNS permeability is not required for COX-2 inhibition against inflammation.
%[text]    Celecoxib targets peripheral COX-2 in inflamed tissues.)
%[text] Q: Create a bar graph of molecular weight with a horizontal line at MW=500.
%[text]    Hint:
%[text]      figure; bar(fullTbl.MolecularWeight);
%[text]      xticks(1:N); xticklabels(fullTbl.Name); xtickangle(20);
%[text]      yline(500, "--r", "Lipinski MW Upper Limit"); grid("on");
%[text]    Which compound is closest to the Lipinski MW upper limit (500)?
% ... (Try writing the code here)
%%
%[text] ## Exercise
%[text] 
%[text] E1: Search PubChem by CID.
%[text]     The CID for the antimalarial drug Quinine is 3034034.
%[text]     Retrieve the properties and draw the structure.
%[text]     Hint:
%[text]       tbl\_q = `emk.db.searchPubchem("3034034", Type="cid")`;
%[text]       mol\_q = `emk.mol.fromSmiles(tbl_q.IsomericSMILES(1)`);
%[text]  `emk.viz.draw2d(mol_q, Title="Quinine")`;
%[text]  `emk.descriptor.calculate(mol_q, ["MolWt","LogP","TPSA","RingCount"])`
%[text] 
%[text] E2: Search by InChIKey.
%[text]     The InChIKey for Aspirin is "BSYNRYMUTXBXSQ-UHFFFAOYSA-N".
%[text]     Confirm that the same CID is returned as with the common name search.
%[text]     Hint: `emk.db.searchPubchem("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", Type="inchikey")`
%[text]     Expected value: CID 2244.
%[text] 
%[text] E3: Expand the panel.
%[text]     Add Diclofenac (another NSAID) and Tramadol (an opioid analgesic)
%[text]     to the panel and reconstruct the comparison table and bar graph.
%[text]     How does the TPSA ranking change?
%[text]     Hint: Add names to PANEL\_NAMES and rerun sections 3-6.
%[text]     Diclofenac SMILES: "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl"
%[text]     Tramadol: Search by common name ("tramadol").
%[text] 
%[text] E4: Error handling.
%[text]     What happens when you search for a name not present in PubChem?
%[text]     Hint: Try `emk.db.searchPubchem("xyzzy_notacompound123")`
%[text]     Expected value: An error with ID emk:db:searchPubchem:notFound occurs.
%[text]     To handle the error gracefully:
%[text]       try
%[text]           r = `emk.db.searchPubchem("xyzzy_notacompound123")`;
%[text]       catch ME
%[text]           logWarn("Not found: %s", ME.message);
%[text]       end
%[text] 
%[text] E5: Create a similarity heatmap for the analgesic panel.
%[text]     Since the SMILES for all 5 compounds have been obtained from PubChem,
%[text]     calculate the Morgan fingerprints and construct the Tanimoto similarity matrix.
%[text]     Which two compounds are the most similar?
%[text]     Hint:
%[text]       fps = cell(1, N);
%[text]       for i = 1:N
%[text]           fps{i} = `emk.fingerprint.morgan(emk.mol.fromSmiles(panelSmiles(i)`));
%[text]       end
%[text]       S = `emk.similarity.matrix(fps)`;
%[text]       figure; imagesc(S); colormap("hot"); colorbar; clim(\[0 1\]);
%[text]       xticks(1:N); xticklabels(PANEL\_NAMES); xtickangle(20);
%[text]       yticks(1:N); yticklabels(PANEL\_NAMES);
%[text]       title("Analgesic Panel -- Tanimoto Similarity Matrix");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
