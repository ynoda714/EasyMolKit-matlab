%[text] # S07: Retrieve Bioactivity Data from ChEMBL
%[text] EasyMolKit Application Story — Layer 2
%[text] 
%[text] What if hundreds of thousands of data points indicating whether a drug works or not were freely available?
%[text] ChEMBL is a bioactivity database managed by the European Bioinformatics Institute, providing IC50 measurements extracted from peer-reviewed journals via a REST API. Thousands of bioactivity data points for EGFR inhibitors are included, and this script will guide you through obtaining data from ChEMBL, filtering for highly active compounds, and exploring their relationship with the existing drug Erlotinib through structural similarity.
%[text] ### Story
%[text] You are a graduate student in a computational drug discovery lab.
%[text] Your advisor has tasked you with using the publicly available bioactivity database ChEMBL to investigate the landscape of known inhibitors of EGFR (epidermal growth factor receptor). EGFR is a receptor tyrosine kinase overexpressed in solid tumors such as non-small cell lung cancer (NSCLC). Since 2004, several small molecule EGFR inhibitors have been approved by the FDA:
%[text] - Gefitinib   (Iressa,   2003, 1st generation)
%[text] - Erlotinib   (Tarceva,  2004, 1st generation)
%[text] - Afatinib    (Gilotrif, 2013, 2nd generation, covalent)
%[text] - Osimertinib (Tagrisso, 2015, 3rd generation, T790M selective) \
%[text] ### Your Tasks:
%[text] 1. Identify EGFR in ChEMBL and confirm the target identifier.
%[text] 2. Retrieve 50 IC50 bioactivity data points for EGFR.
%[text] 3. Examine the IC50 distribution and understand the unit of effect size.
%[text] 4. Filter for highly active inhibitors ($&dollar&;\\text{IC}\_{50} \\le 100\\text{ nM}&dollar&;$).
%[text] 5. Save the active compound set to an SDF file.
%[text] 6. Perform similarity screening against Erlotinib and rank the active compounds.
%[text] 7. Visualize the top hits and compare them with another query compound. \
%[text] ### Learning Objectives
%[text] - Search for ChEMBL protein targets using emk.db.searchChemblTarget
%[text] - Download IC50 bioactivity data using emk.db.getChemblActivity
%[text] - Understand activity units (nM) and effect size thresholds in drug discovery
%[text] - Integrate bioactivity data with structural similarity analysis
%[text] - Export curated compound sets using emk.io.writeSdf \
%[text] ### Prerequisites
%[text] - Completion of F03 (fingerprint) and F04 (similarity)
%[text] - Recommended: S04 (virtual screening) and S06 (PubChem search)
%[text] - RDKit installed (run `emk.setup.install()` once)
%[text] - No additional Toolbox required (works with MATLAB only)
%[text] - Internet connection required (ChEMBL REST API) \
%[text] **Estimated Time**: 35-50 minutes | Execution: Run each section with Ctrl+Enter
%[text] 
%[text] ## Data
%[text] - No local data files needed. Compounds are retrieved in real-time from:
%[text] - ChEMBL REST API  https://www.ebi.ac.uk/chembl/api/data/ \
%[text] ## References
%[text] - Mendez D et al. (2019) ChEMBL: towards direct deposition of bioassay data. *Nucleic Acids Res* 47:D930-D940. doi:10.1093/nar/gky1075 [Open Access]
%[text] - Paez JG et al. (2004) EGFR mutations in lung cancer: correlation with clinical response to gefitinib therapy. *Science* 304:1497-1500. doi:10.1126/science.1099314 [requires institutional access]
%[text] - Yun CH et al. (2008) The T790M mutation in EGFR kinase causes drug resistance by increasing the affinity for ATP. *Proc Natl Acad Sci USA* 105:2070-2075. doi:10.1073/pnas.0709662105 [Open Access]
%[text] - Willett P (2006) Similarity-based virtual screening using 2D fingerprints. *Drug Discov Today* 11:1046-1053. doi:10.1016/j.drudis.2006.10.005 [requires institutional access]\
%[text] ## Note
%[text] This story queries the live ChEMBL REST API. An internet connection is required. Results may vary slightly due to updates in the ChEMBL database.
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
logInfo("S07: Setup complete");
%%
%[text] ## Section 1: Target -- EGFR in Oncology Drug Discovery
%[text] ### Protein Kinase Targets and Success Stories of EGFR
%[text] EGFR (Epidermal Growth Factor Receptor / Gene Name: EGFR) is a "transmembrane receptor tyrosine kinase" that controls the switch for cell growth and survival. It is known as one of the most dramatically successful target proteins in cancer treatment, especially in drug discovery for non-small cell lung cancer (NSCLC).
%[text] [Evolution and Generational Shift of EGFR Inhibitors]
%[text] When mutations (such as exon 19 deletion or L858R mutation) occur in the genes of cancer cells, this switch is fixed to "always ON," causing cells to run amok and tumors to expand. Small molecule inhibitors block the signal by fitting into the "ATP binding pocket" at the active site.
%[text] First Generation (Reversible Inhibitors): Erlotinib, Gefitinib
%[text] - Although they work dramatically on mutant EGFR, long-term use faces the challenge of "T790M (gatekeeper mutation)" where the atoms at the pocket entrance change, rendering the drug unable to inhibit (drug resistance). \
%[text] Second Generation (Irreversible Covalent Inhibitors): Afatinib, Dacomitinib
%[text] - By forming a strong "covalent bond" with the amino acid residue (Cys797) in the pocket, they partially overcome resistance due to the T790M mutation. \
%[text] Third Generation (Mutation-Selective Inhibitors): Osimertinib
%[text] - Specifically designed to strongly inhibit the "T790M mutant," the cause of resistance, while successfully reducing side effects on normal wild-type EGFR. It is now the standard first-choice treatment. \
%[text] [Significance of the ChEMBL Database]
%[text] Behind the development of these successive inhibitors is a vast amount of "compound and activity data" accumulated by researchers worldwide. This has been manually curated by the European Bioinformatics Institute (EMBL-EBI) into "ChEMBL."
%[text] As of 2024, over 2.4 million compounds, more than 20 million bioactivities, and 15,000 target proteins are registered, with thousands of valuable activity data linked to the target human EGFR (CHEMBL203) from past experimental papers.
%[text] [Indicator of Bioactivity: What is IC50]
%[text] It refers to the "concentration of a compound required to inhibit the function of a target protein by 50% (block it halfway)" in an assay (test). The unit is "nM (nanomolar: $&dollar&;10^{-9}&dollar&;$ M)," and the smaller the number (effective at a lower concentration), the more potent it is as an inhibitor. The general activity benchmarks in drug discovery research are as follows:
%[text] - $&dollar&;\\text{IC}\_{50} \\le 100\\text{ nM}&dollar&;$: High activity (excellent quality to proceed to the next step, "lead compound")
%[text] - $&dollar&;100\\text{ nM} \< \\text{IC}\_{50} \\le 1000\\text{ nM}&dollar&;$: Moderate activity (requires optimization/modification of structure)
%[text] - $&dollar&;\\text{IC}\_{50} \> 1000\\text{ nM}&dollar&;$: Weak activity (often not adopted as it cannot surpass the 1 $&dollar&;\\mu\\text{M}&dollar&;$ barrier)
%[text] - $&dollar&;\\text{IC}\_{50} \> 10000\\text{ nM}&dollar&;$: Inactive (requires more than 10 $&dollar&;\\mu\\text{M}&dollar&;$, judged ineffective as a drug) \
EGFR_NAME   = "EGFR";
ERLOTINIB_SMILES = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1";
ERLOTINIB_NAME   = "Erlotinib";

logInfo("Story Target: %s", EGFR_NAME);
logInfo("Reference Compound: %s", ERLOTINIB_NAME);
logInfo("  SMILES: %s", ERLOTINIB_SMILES);
%[text] **✏️ Try It 1**
%[text] Q: What does IC50 stand for, and what is its unit?
%[text]    Why is a lower IC50 considered superior in terms of drug efficacy?
%[text]    Expected Answer: IC50 = Half Maximal Inhibitory Concentration; unit is nM.
%[text]    A lower IC50 means achieving 50% inhibition at a lower concentration (= more potent).
%[text] Q: Erlotinib is a first-generation EGFR inhibitor. Identify the following from the SMILES:
%[text]      (a) Alkyne group (-C#C-) -- Pharmacophore for EGFR back pocket binding.
%[text]      (b) Quinazoline core (bicyclic with two nitrogens) -- Hinge binding motif for hinge region.
%[text]      (c) Two methoxyethoxy chains (-OCCOC-) -- Determinants of solubility and selectivity.
%[text]    Hint: mol = `emk.mol.fromSmiles(ERLOTINIB_SMILES)`; `emk.viz.draw2d(mol)`
% ... (Try writing code here)
%%
%[text] ## Section 2: Search for EGFR in ChEMBL
%[text] ### ChEMBL Target Identifier
%[text] Each of the vast number of proteins registered in ChEMBL is assigned a "ChEMBL Target ID (unique identifier)" in the format `CHEMBL{number}` (e.g., `CHEMBL203` for human EGFR, `CHEMBL230` for human COX-2). This ID serves as an immutable master key to retrieve all literature and compound data related to that protein.
%[text] Normally, you would search for this ID on the website, but with EasyMolKit, you can directly query and obtain it from a program using the `emk.db.searchChemblTarget` function with the official name of the protein (Preferred Name) as a keyword. Let's identify the exact target ID linked to the human (Homo sapiens) EGFR we want to analyze.
logInfo("Searching for EGFR target in ChEMBL...");
targetTbl = emk.db.searchChemblTarget( ...
    "Epidermal growth factor receptor", ...
    TargetType="SINGLE PROTEIN", MaxRows=5);

logInfo("ChEMBL target search results:");
disp(targetTbl(:, ["TargetChEMBLID","PreferredName","Organism"]));
%[text] Confirm that CHEMBL203 is human EGFR. The row with Organism == "Homo sapiens" is the target of interest.
humanIdx = find(targetTbl.Organism == "Homo sapiens", 1, "first");
if isempty(humanIdx)
    % Fallback: Use the known ChEMBL ID directly.
    logWarn("Human EGFR not found in results; using CHEMBL203 directly");
    EGFR_CHEMBL_ID = "CHEMBL203";
else
    EGFR_CHEMBL_ID = targetTbl.TargetChEMBLID(humanIdx);
    logInfo("EGFR target confirmed: %s (%s, %s)", ...
        EGFR_CHEMBL_ID, ...
        targetTbl.PreferredName(humanIdx), ...
        targetTbl.Organism(humanIdx));
end
%[text] **✏️ Try It 2**
%[text] Q: What is the ChEMBL target ID for EGFR?
%[text]    Expected: "CHEMBL203"
%[text] Q: Multiple EGFR entries were returned in the search results. Why?
%[text]    Expected: ChEMBL includes targets from multiple species (human, mouse, rat, etc.). The human (Homo sapiens) entry has the most bioactivity data and is therapeutically important.
%[text] Q: Try searching for another target.
%[text]    What is the ChEMBL ID for human COX-2 (cyclooxygenase-2)?
%[text]    Hint: `emk.db.searchChemblTarget("cyclooxygenase-2", MaxRows=10)`
%[text]    Expected: "CHEMBL230"
% ... (Try writing the code here)

%%
%[text] ## Section 3: Download EGFR IC50 Bioactivity Data
%[text] ### ChEMBL Activity Endpoint
%[text] Once the target ID (master key) is identified, the next step is to download the "measurement data (bioactivity records)" from past experiments conducted on that protein. A single record includes information about the tested compound (ChEMBL compound ID and SMILES structure), experimental conditions, measurement values, units, and a "relation operator" indicating the relationship of the numerical value.
%[text] [Types and Meanings of Relation Operators]
%[text] - `'='` : Accurate data where activity was quantitatively measured (e.g., $&dollar&;\\text{IC}\_{50} = 5.2\\text{ nM}&dollar&;$)
%[text] - `'\<'` : Data indicating strong effect beyond the screening upper limit (e.g., $&dollar&;\\text{IC}\_{50} \< 10\\text{ nM}&dollar&;$)
%[text] - `'\>'` : Data where activity was too weak to measure effect even at maximum experimental concentration (e.g., $&dollar&;\\text{IC}\_{50} \> 10000\\text{ nM}&dollar&;$) \
%[text] If these are mixed, subsequent statistical analysis and graph plotting cannot be calculated correctly. Therefore, the EasyMolKit `emk.db.getChemblActivity` function automatically selects only the most reliable `'='` (accurate equality measurement values), unifies units to `nM`, and filters rows with clean structures (normalized SMILES) that RDKit can handle, returning them as a MATLAB table.
%[text] Currently, ChEMBL has over 5,000 entries for EGFR data, but for this exercise, to ensure smooth API communication, let's limit the number to "a maximum of 50 entries (`MaxRows=50`)" and retrieve the data in real-time.
logInfo("Retrieving EGFR IC50 data from ChEMBL (up to 50 entries)...");
actTbl = emk.db.getChemblActivity(EGFR_CHEMBL_ID, ...
    ActivityType="IC50", MaxRows=50);

% Fallback: use ChEMBL ID as display name when preferred name is absent
emptyName = strlength(actTbl.Name) == 0 | actTbl.Name == "";
actTbl.Name(emptyName) = actTbl.MoleculeChEMBLID(emptyName);

logInfo("Activity dataset: IC50 data for %d compounds", height(actTbl));
logInfo("Columns: %s", strjoin(actTbl.Properties.VariableNames, ", "));
logInfo("IC50 Statistics (nM):");
logInfo("  Minimum  : %.1f nM", min(actTbl.Value_nM));
logInfo("  Median   : %.1f nM", median(actTbl.Value_nM));
logInfo("  Maximum  : %.1f nM", max(actTbl.Value_nM));
%[text] Display the first few rows
disp(actTbl(1:min(5, height(actTbl)), ...
    ["MoleculeChEMBLID","Name","Value_nM"]));
%[text] **Data Note**: The same compound (same MoleculeChEMBLID) may be returned in multiple rows.
%[text] Each row in ChEMBL represents an independent assay measurement, and duplicates occur when the same compound is measured multiple times under different assay conditions (cell lines, biochemical, different labs, etc.).
%[text] This is a characteristic of ChEMBL data, and in actual analysis, it is common to aggregate by compound using `unique` or `groupsummary`.
%[text] **✏️ Try It 3**
%[text] Q: How many compounds have IC50 \< 100 nM (high activity)?
%[text]    Hint: sum(actTbl.Value\_nM \<= 100)
%[text] Q: Which compound has the highest activity in the dataset?
%[text]    Hint: actTbl(actTbl.Value\_nM == min(actTbl.Value\_nM), :)
%[text] Q: Try obtaining Ki values instead of IC50. Are there Ki measurements?
%[text]    Hint: kiTbl = `emk.db.getChemblActivity(EGFR_CHEMBL_ID, ActivityType="Ki", MaxRows=25)`
%[text]    Note: Ki (inhibition constant) is a thermodynamic indicator of binding affinity; IC50 is a functional indicator.
%[text]    Both have units of nM.
% ... (Try writing code here)
%%
%[text] ## Section 4: Investigating IC50 Distribution
%[text] ### Effect Size Distribution and Logarithmic Scale Thinking
%[text] Looking at the statistics (minimum and maximum values) of the $&dollar&;\\text{IC}\_{50}&dollar&;$ obtained earlier, we see that the data spans an enormous range of magnitudes, from "$&dollar&;0.5\\text{ nM}&dollar&;$ (super potent)" to "$&dollar&;3,000,000\\text{ nM}&dollar&;$ (almost just water)." Attempting to graph such wide-ranging data in raw numbers ($&dollar&;\\text{nM}&dollar&;$) would crush the differences of highly active compounds around "1-100 $&dollar&;\\text{nM}&dollar&;$" into the zero vicinity, making them indistinguishable.
%[text] Therefore, medicinal chemists standardly use the "pIC50" index, which is the negative common logarithm of the $&dollar&;\\text{IC}\_{50}&dollar&;$ value converted to molar concentration ($&dollar&;\\text{M}&dollar&;$). $&dollar&;\\text{pIC}\_{50} = -\\log\_{10}(\\text{IC}\_{50}\\text{ molar concentration}) = 9 - \\log\_{10}(\\text{Value\\\_nM})&dollar&;$
%[text] [Benefits and Intuitive Correspondence of the $&dollar&;\\text{pIC}\_{50}&dollar&;$ Scale]
%[text]  "Larger numbers indicate stronger potency": Raw $&dollar&;\\text{IC}\_{50}&dollar&;$ values are counterintuitive as smaller numbers indicate stronger potency, but $&dollar&;\\text{pIC}\_{50}&dollar&;$ increases with potency.
%[text]  Treats magnitudes as equidistant (linear): A tenfold difference in activity is neatly plotted as a "1" difference, making statistical analysis and histogram visualization significantly easier.
%[text] - $&dollar&;\\text{pIC}\_{50} \> 7&dollar&;$ ($&dollar&;\\text{IC}\_{50} \< 100\\text{ nM}&dollar&;$): Lead compounds with excellent activity
%[text] - $&dollar&;\\text{pIC}\_{50} \> 8&dollar&;$ ($&dollar&;\\text{IC}\_{50} \< 10\\text{ nM}&dollar&;$): High activity comparable to existing drugs
%[text] - $&dollar&;\\text{pIC}\_{50} \> 9&dollar&;$ ($&dollar&;\\text{IC}\_{50} \< 1\\text{ nM}&dollar&;$): Ultra-potent molecules seen in covalent inhibitors \
%[text] Let's leverage MATLAB's vector operations to instantly add a `pIC50` column to the downloaded table and visualize the effect size distribution with a beautiful histogram.
pIC50 = 9 - log10(actTbl.Value_nM);
actTbl.pIC50 = pIC50;

logInfo("pIC50 Statistics:");
logInfo("  Minimum  : %.2f (IC50 = %.1f nM)", min(pIC50), max(actTbl.Value_nM));
logInfo("  Median   : %.2f (IC50 = %.1f nM)", median(pIC50), median(actTbl.Value_nM));
logInfo("  Maximum  : %.2f (IC50 = %.1f nM)", max(pIC50), min(actTbl.Value_nM));
%[text] Number of Compounds by Effect Size Class
nHigh    = sum(actTbl.Value_nM <= 100);
nMod     = sum(actTbl.Value_nM > 100  & actTbl.Value_nM <= 1000);
nWeak    = sum(actTbl.Value_nM > 1000);

logInfo("Effect Size Distribution:");
logInfo("  High Activity     (IC50 <= 100 nM)    : %d entries", nHigh);
logInfo("  Moderate Activity (100 < IC50 <= 1000): %d entries", nMod);
logInfo("  Weak/Inactive     (IC50 > 1000 nM)   : %d entries", nWeak);
%[text] Histogram of pIC50 Values (Using Base MATLAB Only)
figure("Name", "EGFR IC50 Distribution", "Position", [100 100 580 420]);
histogram(pIC50, 10, "FaceColor", [0.2 0.5 0.8]);
xlabel("pIC50  (-log_{10}[IC50 in M])");
ylabel("Number of Compounds");
title(sprintf("EGFR IC50 Distribution  (n=%d, from ChEMBL)", height(actTbl)));
xline(7, "r--", "LineWidth", 1.5, "Label", "pIC50=7 (100 nM)");
xline(8, "g--", "LineWidth", 1.5, "Label", "pIC50=8 (10 nM)");
grid on;
%[text] **✏️ Try It 4**
%[text] Q: What is the pIC50 corresponding to Erlotinib's IC50 ~2 nM against wild-type EGFR?
%[text]    Hint: Use the formula pIC50 = 9 - log10(IC50\_nM).
%[text] Q: Is the distribution skewed? If so, in which direction and why?
%[text]    Expected: Typically right-skewed (few ultra-potent compounds, many moderately active compounds).
%[text]    The ChEMBL IC50 dataset often shows a log-normal distribution as it reports both active and inactive compounds in the same series.
%[text] Q: The histogram displays pIC50 on the horizontal axis.
%[text]    Would a bar graph with raw IC50 values (in nM) on a linear scale increase or decrease information content?
%[text]    Why do medicinal chemists prefer a logarithmic scale?
%[text]    Expected: On a linear scale, all values compress towards the lower end;
%[text]    pIC50 (logarithmic scale) spreads the distribution, visualizing differences among potent compounds.
% ... (Try writing code here)

%%
%[text] ## Section 5: Filter High-Activity Compounds ($&dollar&;\\text{IC}\_{50} \\le 100\\text{ nM}&dollar&;$)
%[text] ### Setting the Activity Cutoff
%[text] The collected distribution data includes many "inactive compounds" with little effect as drugs. To efficiently advance the project, it is necessary to set an "activity cutoff (threshold)" to narrow down to promising molecules worth investing resources in.
%[text] [Common Cutoff Criteria in Drug Discovery Stages]
%[text] - Primary Screening (filtering tens of thousands): $&dollar&;\\text{IC}\_{50} \< 10,000\\text{ nM}&dollar&;$ ($&dollar&;10\\text{ }\\mu\\text{M}&dollar&;$)
%[text] - Hit-to-Lead (starting point for serious consideration): $&dollar&;\\text{IC}\_{50} \< 1,000\\text{ nM}&dollar&;$ ($&dollar&;1\\text{ }\\mu\\text{M}&dollar&;$)
%[text] - Lead Optimization (aiming for refinement as a drug): $&dollar&;\\text{IC}\_{50} \\le 100\\text{ nM}&dollar&;$
%[text] - Clinical Development Candidate (level to compete as an actual drug): $&dollar&;\\text{IC}\_{50} \< 10\\text{ nM}&dollar&;$ \
%[text] In this section, we set the cutoff criterion as "$&dollar&;\\text{IC}\_{50} \\le 100\\text{ nM}&dollar&;$", the passing standard for lead optimization. Using MATLAB's logical indexing feature, let's swiftly extract (filter) only high-quality compounds that meet this stringent elite standard from 50 miscellaneous data, and further sort them in order of high activity (most potent first).
IC50_CUTOFF_NM = 100;

activeMask = actTbl.Value_nM <= IC50_CUTOFF_NM;
activeTbl  = actTbl(activeMask, :);

logInfo("Activity Cutoff: IC50 <= %.0f nM", IC50_CUTOFF_NM);
logInfo("Active Compounds: %d / %d entries", height(activeTbl), height(actTbl));

%[text] Sort in ascending order of IC50 (most potent first)
activeTbl = sortrows(activeTbl, "Value_nM", "ascend");
logInfo("Top 5 EGFR Inhibitors (by IC50):");
disp(activeTbl(1:min(5, height(activeTbl)), ...
    ["MoleculeChEMBLID","Name","Value_nM","pIC50"]));

%[text] **✏️ Try It 5**
%[text] Q: How many active compounds remain after applying the 100 nM cutoff?
%[text]    If very few, try relaxing it to 1000 nM.
%[text] Q: What is the IC50 of the most active compound?
%[text]    Can you find its structure and name on PubChem or ChEMBL?
%[text]     Hint: disp(activeTbl(1, ["MoleculeChEMBLID","Name","Value_nM","SMILES"]))
%[text] Q: Change IC50\_CUTOFF\_NM to 1000 and rerun this section.
%[text]    How many compounds qualify?
%[text]    Is that number reasonable for a fragment-based or HTS-derived EGFR dataset?
% ... (Try writing some code here)
%%
%[text] ## Section 6: Save Active Compounds to SDF
%[text] ### SDF (Structure-Data File) Format
%[text] A well-curated list of highly active compounds should not remain only on the MATLAB screen if you want to share it with other researchers. Therefore, we will export it locally in the "SDF (Structure-Data File)" format, which is the industry standard in chemistry.
%[text] An SDF file is simply a text file, but internally it contains two important pieces of information as a single record:
%[text] 1. MOL Block: Types of atoms, how they are bonded, and 2D or 3D coordinate information (bond table)
%[text] 2. Data Fields: Text metadata associated with the compound (e.g., ChEMBL ID, $&dollar&;\\text{IC}\_{50}&dollar&;$ values, assay conditions, etc.) \
%[text] This file format is highly versatile and can be directly read by any chemical software worldwide, such as ChemDraw (structure drawing), Schrodinger Maestro and MOE (molecular modeling), and KNIME (data flow construction). Experience the standard workflow of "a data analyst (computational chemist) selecting attractive hit compounds and exporting them to an SDF file to request synthesis or ordering from an experimentalist (medicinal chemist)" that is performed daily in real pharmaceutical settings using the `emk.io.writeSdf` function.
if height(activeTbl) == 0
    logWarn("No active compounds to save -- Skipping SDF export");
else
    % Convert SMILES of active compounds to RDKit Mol objects.
    logInfo("Converting %d active compounds...", height(activeTbl));
    validMols = {};
    validIdx  = [];

    for i = 1:height(activeTbl)
        smi = activeTbl.SMILES(i);
        if ~emk.mol.isValid(smi)
            logWarn("  [%d] Skipping due to invalid SMILES: %s", i, smi);
            continue;
        end
        mol = emk.mol.fromSmiles(smi);
        validMols{end+1} = mol; %#ok<AGROW>
        validIdx(end+1)  = i;   %#ok<AGROW>
    end

    logInfo("Conversion successful: %d / %d", numel(validMols), height(activeTbl));

    if numel(validMols) > 0
        runDir  = makeRunDir("Prefix", "s07_egfr");
        sdfPath = fullfile(runDir, "egfr_actives.sdf");
        emk.io.writeSdf(validMols, sdfPath);
        logInfo("Saved %d active compounds: %s", numel(validMols), sdfPath);

        % Save metadata as CSV next to the SDF
        csvPath = fullfile(runDir, "egfr_actives_metadata.csv");
        writetable(activeTbl(validIdx, :), csvPath);
        logInfo("Saved metadata CSV: %s", csvPath);
    end
end

%[text] **✏️ Try It 6**
%[text] Q: Open the SDF file in a text editor and check the following:
%[text]    (a) MOL Block (ends with "M  END")
%[text]    (b) Record separator ("$$$$")
%[text]    (c) Does the file contain coordinates?
%[text]        RDKit automatically generates 2D coordinates when writing SDF files.
%[text] Q: How to reload the SDF into MATLAB?
%[text]    Hint: reloaded = `emk.io.readSdf(sdfPath)`
%[text]    (if your version of EasyMolKit includes `emk.io.readSdf`)
%[text] Q: Count the number of records in the SDF file:
%[text]    Hint: txt = fileread(sdfPath); sum(contains(strsplit(txt, newline), "$$$$"))
% ... (Try writing code here)
%%
%[text] ## Section 7: Similarity Screening for Erlotinib
%[text] ### Hit Selection by Similarity Guide
%[text] The next question we should ask after obtaining a set of carefully selected high-activity inhibitors is, "Do any of these active compounds resemble an existing drug that has already succeeded in the clinic (Erlotinib)?"
%[text] Finding compounds structurally similar (sharing the scaffold) to already approved drugs has significant advantages:
%[text] - Safety Prediction: Likely to have a scaffold that is less prone to severe toxicity or unexpected side effects (off-target toxicity)
%[text] - Knowledge Transfer: Can leverage vast past research data (Structure-Activity Relationship: SAR) to speed up development
%[text] - Ease of Synthesis: Can easily apply existing manufacturing and synthesis routes (analog expansion) \
%[text] To have the computer calculate this "degree of structural similarity," we compute the "ECFP4 (Morgan radius=2) fingerprint," which scans the local chemical environment within two bonds of each atom into a 2048-bit digital fingerprint. Then, we calculate the "Tanimoto similarity (Tanimoto Metric)" which scores the overlap of these fingerprints in the range of `0.0` (completely different) to `1.0` (fingerprints match).
%[text] Using the `emk.similarity.rankBy` function, let's sort the active compounds in ChEMBL by their similarity scores to the reference "Erlotinib" and create a ranking table. By visually confirming the 2D structures of the compounds ranked at the top using `emk.viz.draw2d`, you will experience the powerful approach of ligand-based drug discovery.
if ~exist("validMols", "var") || numel(validMols) == 0
    logWarn("No valid active mols exist -- please run Section 6 first");
else
    logInfo("Running similarity screening for %s...", ERLOTINIB_NAME);

    % Calculate query fingerprint (Erlotinib)
    mol_erl = emk.mol.fromSmiles(ERLOTINIB_SMILES);
    fp_erl  = emk.fingerprint.morgan(mol_erl, Radius=2, NBits=2048);

    % Calculate fingerprints for all active compounds
    logInfo("Calculating ECFP4 fingerprints for %d active compounds...", numel(validMols));
    dbFps = cell(1, numel(validMols));
    for i = 1:numel(validMols)
        dbFps{i} = emk.fingerprint.morgan(validMols{i}, Radius=2, NBits=2048);
        logProgress(i, numel(validMols), "Fingerprint Calculation");
    end

    % Rank by Tanimoto similarity to Erlotinib
    rankResult = emk.similarity.rankBy(fp_erl, dbFps, Inf, Metric="tanimoto");

    % Construct result table
    rankedIdx    = rankResult.Indices;
    rankedScores = rankResult.Scores;
    nRanked      = numel(rankedIdx);

    metaTbl = activeTbl(validIdx, :);
    rankTbl = table( ...
        (1:nRanked)', ...
        metaTbl.MoleculeChEMBLID(rankedIdx(:)), ...
        metaTbl.Name(rankedIdx(:)), ...
        metaTbl.Value_nM(rankedIdx(:)), ...
        rankedScores(:), ...
        VariableNames=["Rank","ChEMBLID","Name","IC50_nM","Tanimoto_Erlotinib"]);

    logInfo("Top 10 active compounds most similar to Erlotinib (Tanimoto, ECFP4):");
    disp(rankTbl(1:min(10, nRanked), :));
%[text] **Note**: The same ChEMBLID may appear in multiple ranks (corresponding to duplicate rows in Section 3).
%[text] Each row represents an independent IC50 measurement, especially seen in well-studied compounds.

    logInfo("Tanimoto score statistics:");
    logInfo("  Max  : %.3f", max(rankedScores));
    logInfo("  Median: %.3f", median(rankedScores));
    logInfo("  Min  : %.3f", min(rankedScores));

    % Number of compounds with Tanimoto > 0.4 (scaffold similarity threshold)
    nSimilar = sum(rankedScores > 0.4);
    logInfo("Tanimoto > 0.4 with Erlotinib: %d / %d", nSimilar, nRanked);

    % Visualize top hit structures
    logInfo("Visualizing top hit: %s (IC50 = %.1f nM, Tanimoto = %.3f)", ...
        metaTbl.MoleculeChEMBLID(rankedIdx(1)), ...
        metaTbl.Value_nM(rankedIdx(1)), ...
        rankedScores(1));
    emk.viz.draw2d(validMols{rankedIdx(1)}, ...
        Title=sprintf("%s IC50=%.1f nM Tan=%.3f", ...
            metaTbl.MoleculeChEMBLID(rankedIdx(1)), ...
            metaTbl.Value_nM(rankedIdx(1)), ...
            rankedScores(1)));
end
%[text] **✏️ Try It 7**
%[text] Q: What is the Tanimoto score of the compound most similar to Erlotinib?
%[text]    Is it structurally very close (Tanimoto \> 0.7) or moderately similar?
%[text] Q: Investigate the top compounds in ChEMBL or PubChem.
%[text]    Are they known EGFR inhibitors? Does the name suggest a relation to Erlotinib?
%[text]    Hint:
%[text]      topChEMBLID = rankTbl.ChEMBLID(1);
%[text]      disp(rankTbl(1, ["ChEMBLID","Name","Tanimoto"]))
%[text] Q: Compare the query compound with Erlotinib and Gefitinib.
%[text]    Gefitinib SMILES: "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1"
%[text]    Rank the same active set against Gefitinib and compare the top 10:
%[text]    Hint:
%[text]      mol\_gef = emk.mol.fromSmiles("COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1");
%[text]      fp\_gef  = `emk.fingerprint.morgan(mol_gef, Radius=2, NBits=2048)`;
%[text]      res\_gef = `emk.similarity.rankBy(fp_gef, dbFps, 10, Metric="tanimoto")`;
%[text]    Do Erlotinib and Gefitinib identify the same top compounds?
%[text]    What can be said about the chemical diversity of EGFR inhibitors?
% ... (Try writing code here)
%%
%[text] ## Section 8: Summary and Key Points
%[text] ### What We Did in This Story
%[text] 1. Target Identification: Used `emk.db.searchChemblTarget` to identify EGFR (CHEMBL203) from the ChEMBL target database.
%[text] 2. Data Download: Used `emk.db.getChemblActivity` to retrieve 50 $ &dollar&;\\textrm{IC}\_{50}&dollar&;$ records for EGFR from peer-reviewed ChEMBL literature.
%[text] 3. Effect Size Analysis: Examined the distribution of $&dollar&;\\textrm{IC}\_{50}&dollar&;$, converted to $&dollar&;\\textrm{pIC}\_{50}&dollar&;$ for intuitive understanding on a logarithmic scale, and learned about effect size thresholds in drug discovery.
%[text] 4. Activity Filtering: Applied an activity cutoff of $&dollar&;{\\textrm{IC}}\_{50} \\le 100\\;\\textrm{nM}&dollar&;$ to select only high-quality lead compounds.
%[text] 5. Data Saving: Used `emk.io.writeSdf` to save the active compound set in an SDF file, the industry-standard format for chemical structure exchange.
%[text] 6. Similarity Screening: Calculated ECFP4 fingerprints and experienced a ligand-based hit selection strategy by ranking active compounds by Tanimoto similarity to the existing drug Erlotinib. \
%[text] ### Key Concepts Introduced
%[text] - ChEMBL: A curated bioactivity database for drug-like small molecules
%[text] - $&dollar&;\\textrm{IC}\_{50}&dollar&; / &dollar&;\\textrm{pIC}\_{50}&dollar&;$: Indicators of effect size and logarithmic scale thinking
%[text] - Activity Cutoff: Choosing appropriate thresholds according to the drug discovery stage
%[text] - SDF Format: Industry-standard format for managing chemical structures and metadata in bulk
%[text] - Structural Similarity: Ligand-based hit selection using fingerprints \
%[text] ### Next Steps
%[text] - S04 Virtual Screening: Conduct a full ligand-based virtual screening (LBVS) against an FDA-approved drug library
%[text] - A01 Chemical Space Mapping: Visualize the diversity of EGFR inhibitors using Principal Component Analysis (PCA)
%[text] - A03 QSAR Regression: Build a machine learning model to predict $&dollar&;\\textrm{pIC}\_{50}&dollar&;$ from compound descriptors \
logInfo("S07 Complete -- EGFR Bioactivity Analysis Finished.");
if exist("runDir", "var")
    logInfo("File Save Location: %s", runDir);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
