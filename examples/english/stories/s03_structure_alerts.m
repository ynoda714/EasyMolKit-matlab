%[text] # S03: Structural Alerts — Detecting Hazardous Functional Groups with SMARTS
%[text] EasyMolKit Application Story — Layer 2
%[text] 
%[text] If there were "chemical time bombs" hidden among the vast number of new compounds sent to the laboratory, directly reacting with proteins or DNA to cause damage, how would you find them?
%[text] Such reactive sites harmful to living organisms (structural alerts) can be instantly detected and screened on a computer by writing just one line of a special string called **SMARTS pattern**.
%[text] Furthermore, in the pharmaceutical industry, there is a list of very troublesome false positive substances called PAINS (Pan-Assay Interference Compounds), which "do not actually work but trick the device sensors into appearing as drugs" during experiments, becoming a "source of pain" for researchers.
%[text] In this script, you will experience the process of applying 8 types of custom reactivity alerts and the strictest PAINS filter to 200 FDA-approved drugs in bulk, automatically identifying problematic molecules and generating a report.
%[text] ## Learning Objectives
%[text] - Understand the basics of "SMARTS patterns," which are conditional expressions for searching parts of a structure
%[text] - Use `emk.mol.hasSubstruct` to perform bulk substructure searches (vectorized screening) on multiple molecules without loops
%[text] - Understand the concept of "PAINS" in drug discovery and the importance of eliminating assay false positives (triage)
%[text] - Automatically organize screening results and create structured flag tables and CSV reports \
%[text] ## Prerequisites
%[text] - Completion of F05 (Substructure Search)
%[text] - RDKit installed (execute `emk.setup.install()` once)
%[text] - No additional Toolbox required (works with MATLAB alone) \
%[text] **Estimated Time**: 15–20 minutes | Execution: Run each section with Ctrl+Enter
%[text] **Data**
%[text] - `data/list/pains.csv` — 480 PAINS SMARTS (RDKit wehi\_pains, BSD-3)
%[text] - `data/list/fda_drugs.csv` — 200 FDA-approved drugs (ChEMBL, CC-BY-SA 3.0) \
%[text] **References**
%[text] - Baell JB & Holloway GA (2010) *J Med Chem* 53:2719–2740. doi:10.1021/jm901137j [requires institutional access]\
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
logInfo("S03: Setup complete");
%%
%[text] ## Section 1: What are Structural Alerts?
%[text] ### Concept: Chemical groups that "directly react" with biomolecules
%[text] Some chemical groups have the property of directly reacting with biomolecules such as proteins and DNA. In the early stages of drug discovery, these groups are flagged as "structural alerts" and are deprioritized unless there is a clear reason not to.
%[text] Here are some well-known alerts and examples of their risks:
%[text] -   Epoxide          `[C]1CO1`           — Alkylates DNA and proteins
%[text] -   Michael acceptor `[C]=[C]-C=O`       — 1,4-addition with cysteine (-SH)
%[text] -   Aldehyde         `[CH]=O`            — Forms Schiff base with amines
%[text] -   Acyl halide      `[C](=O)[F,Cl]`     — Highly reactive electrophile
%[text] -   Nitro group      `[N+](=O)[O-]`      — Potential mutagenicity (Ames test)
%[text] -   Hydrazine        `[NH]-[NH2]`        — Concerns of hepatotoxicity
%[text] -   Diazo            `[#6]=[N+]=[N-]`    — Highly reactive and unstable
%[text] -   Thiol            `[SH]`              — Can also be a pharmacophore \
%[text] However, the presence of an alert does not mean a compound is "unusable." For example, the ester in Aspirin is designed to be intentionally hydrolyzed in the body (prodrug). "Flagged" is merely a signal for experts to make a judgment.
%[text] Define the alert panel (name → SMARTS correspondence).
ALERT_NAMES  = ["Epoxide",      "Michael_acceptor", "Aldehyde", ...
                "Acyl_halide",  "Nitro",             "Hydrazine", ...
                "Diazo",        "Thiol"];
ALERT_SMARTS = ["[C]1CO1",      "[C]=[C]-C=O",      "[CH]=O", ...
                "[C](=O)[F,Cl]","[N+](=O)[O-]",     "[NH]-[NH2]", ...
                "[#6]=[N+]=[N-]","[SH]"];

logInfo("Alert panel loaded: %d patterns", numel(ALERT_NAMES));
for k = 1:numel(ALERT_NAMES)
    logInfo("  %-18s  %s", ALERT_NAMES(k), ALERT_SMARTS(k));
end
%%
%[text] ## Section 2: Screening a Test Set of Known Problematic Molecules
%[text] We have prepared a test set of "representative problematic molecules" known to contain reactive functional groups. Let's use EasyMolKit to verify if alerts can be correctly detected.
TEST_NAMES  = ["Epichlorohydrin",  "Acrolein",    "Formaldehyde", ...
               "Malondialdehyde",  "Nitrobenzene", "Phenylhydrazine", ...
               "Aspirin",          "Captopril"];
TEST_SMILES = ["ClCC1CO1",          "C=CC=O",       "C=O", ...
               "O=CCC=O",          "c1ccc([N+](=O)[O-])cc1", ...
               "c1ccccc1NN",       "CC(=O)Oc1ccccc1C(=O)O", ...
               "CC(CS)C(=O)N1CCCC1C(=O)O"];

%[text] ### Concept: Vectorized Screening by Passing Cell Arrays
%[text] By passing a cell array of Mol objects as the first argument to `emk.mol.hasSubstruct`, you can receive the results for each molecule as a `logical(1, N)` row vector all at once. You can screen all molecules in bulk without writing explicit loops.
nMols   = numel(TEST_SMILES);
nAlerts = numel(ALERT_NAMES);
flagMat = false(nMols, nAlerts);   % Rows = molecules, Columns = alerts

mols_test = cell(1, nMols);
for i = 1:nMols
    mols_test{i} = emk.mol.fromSmiles(TEST_SMILES(i));
end

for j = 1:nAlerts
    flagMat(:, j) = emk.mol.hasSubstruct(mols_test, ALERT_SMARTS(j))';
end

%[text] Format the results into a readable table.
flagTbl = array2table(flagMat, "VariableNames", cellstr(ALERT_NAMES));
flagTbl.Molecule     = TEST_NAMES';
flagTbl.TotalAlerts  = sum(flagMat, 2);
flagTbl = movevars(flagTbl, ["Molecule", "TotalAlerts"], "Before", "Epoxide");

logInfo("--- Structural Alert Screen (Test Set) ---");
disp(flagTbl(:, ["Molecule", "TotalAlerts", ALERT_NAMES]));

%[text] **✏️ Try It 1 — Read the Table**
%[text] Check the table and answer the following two questions.
%[text] - Which test molecule triggers the most alerts?
%[text] - Did Aspirin trigger any alerts? What does that mean? \
%[text] 
%[text] **Expected Values**:
%[text] - Acrolein (`C=CC=O`) triggers **2 alerts**: Aldehyde and Michael\_acceptor.
%[text] - Malondialdehyde triggers only Aldehyde (1 alert — no Michael acceptor due to lack of C=C bond).
%[text] - Formaldehyde (`C=O`) triggers **0 alerts** — `[CH]=O` requires 1 H on carbon, \
%[text]   but Formaldehyde's carbon has 2 Hs (`[CH2]=O`), so it doesn't match.
%[text]   This is a good example of the specificity of SMARTS. Remember that it is not flagged despite being "aldehyde-like."
%[text] - Aspirin has no alerts in this panel — esters are not flagged.
%[text] - Captopril (an ACE inhibitor) triggers Thiol (intentional pharmacophore). \
%[text] 
%[text] "Flagged" is a trigger for expert judgment, not a basis for automatic exclusion.
% ... (Try writing code here)
%%
%[text] ## Section 3: Drawing Flagged Molecules
%[text] Let's visualize molecules in the test set that triggered at least one alert.
flaggedIdx = find(sum(flagMat, 2) > 0);
nFlagged_test = numel(flaggedIdx);
logInfo("Drawing %d flagged molecules in the test set...", nFlagged_test);

%[text] Draw flagged molecules in a grid (same pattern as Section F01).
nCols_f = ceil(sqrt(nFlagged_test));
nRows_f = ceil(nFlagged_test / nCols_f);
figure("Name", "Flagged Molecules -- Structural Alerts", ...
    "Position", [100 100 nCols_f*280 nRows_f*260]);
for k = 1:nFlagged_test
    i = flaggedIdx(k);
    subplot(nRows_f, nCols_f, k);
    alertList = ALERT_NAMES(flagMat(i, :));
    titleStr = sprintf("%s [%s]", TEST_NAMES(i), strjoin(alertList, ", "));
    emk.viz.draw2d(mols_test{i}, Title=titleStr);
end

%[text] **✏️ Try It 2 — Check the Two Reactive Centers of Acrolein**
%[text] Draw Acrolein (CH2=CH-CHO) and consider which parts are more reactive.
%[text] 
%[text]     mol\_acrolein = emk.mol.fromSmiles("C=CC=O");
%[text]     emk.viz.draw2d(mol\_acrolein, Title="Acrolein — Michael + Aldehyde", ...
%[text]         Width=300, Height=250);
%[text] 
%[text] Does the thiol (-SH) of cysteine prefer to attack the aldehyde carbon or the β carbon?
%[text] Hint: Thiols are "soft" nucleophiles. Recall the HSAB theory.
%[text] 
%[text] **Expected Outcome**: It attacks the β carbon (the second carbon in C=C, Michael acceptor site).
%[text] Soft nucleophiles favor 1,4-addition due to the LUMO extending to the terminal of the conjugated system (β position), making it kinetically favorable.
% ... (Try writing the code here)
%%
%[text] ## Section 4: Load and Explore the PAINS Filter Set
%[text] ### Concept: What is PAINS?
%[text] Baell & Holloway (2010) analyzed 6.8 million assay data points and identified 480 SMARTS patterns that characterize compounds causing recurrent false positives in multiple biochemical assays. These are known as "PAINS."
%[text] The main causes include fluorescence interference, redox cycling, colloidal aggregation, and optical interference. The 480 patterns are classified into three subfilters based on stringency.
%[text] - **PAINS\_A**: 16 patterns (most stringent)
%[text] - **PAINS\_B**: 55 patterns
%[text] - **PAINS\_C**: 409 patterns (widest coverage) \
%[text] 
%[text] **Implementation Note**: The SMARTS strings in the PAINS CSV contain commas within quoted fields.
%[text] Using `textscan` with `%q` (quoted string format) correctly handles embedded delimiters.
fid = fopen(fullfile(projectRoot, "data", "list", "pains.csv"), "r");
textscan(fid, "%s", 1, "Delimiter", newline);  % Consume header line
C = textscan(fid, "%q%q%q%q", "Delimiter", ",");   % 4 columns: Name,SMARTS,FilterSet,Source
fclose(fid);
painsCsv = table(string(C{1}), string(C{2}), string(C{3}), ...
    VariableNames=["Name", "SMARTS", "FilterSet"]);

logInfo("PAINS database: %d patterns loaded", height(painsCsv));
filterSets = unique(painsCsv.FilterSet);
for k = 1:numel(filterSets)
    n = sum(painsCsv.FilterSet == filterSets(k));
    logInfo("  %-10s: %d patterns", filterSets(k), n);
end

%[text] **✏️ Try It 3 — Explore the PAINS Database**
%[text] Open `painsCsv` in the MATLAB workspace and check the first 3 patterns.
%[text] Hint: Display the table with `painsCsv(1:3, ["Name","SMARTS","FilterSet"])`.
%[text] 
%[text] The SMARTS of PAINS are complex recursive queries. It's more accurate to let RDKit's SMARTS parser handle them than to attempt a manual explanation.
% ... (Try writing some code here)
%%
%[text] ## Section 5: Batch Screening of PAINS on FDA Drugs
%[text] Apply PAINS patterns to 200 FDA drugs in bulk.
%[text] Record the following three for each drug:
%[text] - Whether it matches any PAINS pattern (`isPainsFlagged`)
%[text] - How many different patterns it matches (`painsCount`)
%[text] - The name of the first triggered pattern (`firstAlert`) \
%[text] **Performance Note**: 480 patterns × 200 drugs = 96,000 IPC calls, taking 10-30 minutes to complete. This demo uses only the strictest PAINS_A (16 patterns) to complete within 30 seconds. For screening all 480 patterns, see Try It 4.
rawDrugs = readtable(fullfile(projectRoot, "data", "list", "fda_drugs.csv"), "TextType", "string");
nDrugs   = height(rawDrugs);

%[text] For interactive demo, use PAINS_A (16 types, strictest). To try all 480 patterns, change this line to `painsScreen = painsCsv;`.
painsScreen = painsCsv(painsCsv.FilterSet == "PAINS_A", :);
nPains      = height(painsScreen);

logInfo("Screening %d drugs with %d PAINS_A patterns...", nDrugs, nPains);
logInfo("  (PAINS_A = strictest %d patterns; full database is %d types)", nPains, height(painsCsv));

%[text] First, parse the SMILES of all drugs (skip with a warning for invalid SMILES).
mols_drug = cell(1, nDrugs);
validDrug  = false(1, nDrugs);
for i = 1:nDrugs
    if emk.mol.isValid(rawDrugs.SMILES(i))
        mols_drug{i} = emk.mol.fromSmiles(rawDrugs.SMILES(i));
        validDrug(i) = true;
    else
        logWarn("  Skipping: %s (invalid SMILES)", rawDrugs.Name(i));
    end
end
logInfo("  Parsed drug structures: %d / %d.", sum(validDrug), nDrugs);

%[text] Screen all drugs with each PAINS pattern. `hitMat(i, j) = true` means drug i hit pattern j.
hitMat = false(nDrugs, nPains);
mols_valid_list = mols_drug(validDrug);
validIdx        = find(validDrug);

for j = 1:nPains
    logProgress(j, nPains, "PAINS Screen");
    smarts_j = painsScreen.SMARTS(j);
    hits_j   = emk.mol.hasSubstruct(mols_valid_list, smarts_j);   % 1×nValid logical
    hitMat(validIdx, j) = hits_j;
end

%[text] Summarize hit counts and first alert name for each drug.
painsCount = sum(hitMat, 2);          % Number of patterns each drug hits
isPainsFlagged = painsCount > 0;

%[text] Get the name of the first hit pattern (empty if no hits).
firstAlert = repmat("", nDrugs, 1);
for i = 1:nDrugs
    idx = find(hitMat(i, :), 1);
    if ~isempty(idx)
        firstAlert(i) = painsScreen.Name(idx);
    end
end

%[text] Compile results into a table.
resultTbl = table(rawDrugs.Name, rawDrugs.SMILES, isPainsFlagged, painsCount, firstAlert, ...
    VariableNames=["Name", "SMILES", "PAINS_Flagged", "PAINS_Count", "First_Alert"]);

nFlagged = sum(isPainsFlagged);
logInfo("PAINS screening complete:");
logInfo("  Flagged: %d / %d  (%.0f%%)", nFlagged, nDrugs, 100*nFlagged/nDrugs);
logInfo("  Clean   : %d / %d  (%.0f%%)", nDrugs - nFlagged, nDrugs, ...
    100*(nDrugs-nFlagged)/nDrugs);

%[text] **✏️ Try It 4 — Investigate the PAINS Proportion in FDA Drugs**
%[text] - What percentage of FDA drugs triggered at least one PAINS_A pattern?
%[text] - Are the results surprising? Consider the reasons. \
%[text] 
%[text] **Expected Value**: Typically, 3-5% of approved drugs are flagged with PAINS_A (16 patterns). With all 480 patterns, it rises to 15-25% (PAINS_B/C are more lenient patterns). It's natural for approved drugs to have many flags — reactive motifs themselves can be pharmacophores.
%[text] To try screening all 480 patterns (~10-30 minutes), change `painsScreen = painsCsv;` in Section 5.
% ... (Try writing code here)


%%
%[text] ## Section 6: Identify the Most Frequently Triggered PAINS Patterns
%[text] Let's aggregate the trigger counts for each pattern to identify the top alerts.
patternHits = sum(hitMat, 1);   % 1×nPains

%[text] Sort in descending order of trigger counts.
[sortedHits, sortOrder] = sort(patternHits, "descend");
topK = min(10, sum(sortedHits > 0));   % Top 10 or number of triggers

if topK == 0
    logInfo("No PAINS patterns were triggered in this dataset.");
else
    topPatternNames = painsScreen.Name(sortOrder(1:topK));   % already column vector
    topHitCounts    = sortedHits(1:topK)';
    topTbl = table(topPatternNames, topHitCounts, ...
        VariableNames=["PatternName", "FlaggedDrugs"]);
    logInfo("Top %d PAINS patterns:", topK);
    disp(topTbl);

    % Bar graph: Only plot if there is variation in hit counts
    % In the PAINS_A demo (16 patterns, few hits), all pattern counts tend to be uniform,
    % so the bar graph shows no significant difference. Skip it. Significant differences appear with all 480 patterns.
    if max(topHitCounts) - min(topHitCounts) > 0
        figure("Name", "Top PAINS Patterns", "Color", "white", "Position", [100 100 580 380]);
        topNames = painsScreen.Name(sortOrder(1:topK));
        barh(topK:-1:1, sortedHits(1:topK), "FaceColor", [0.85 0.33 0.10]);
        yticks(1:topK);
        yticklabels(flip(cellstr(topNames)));
        xlabel("Number of Flagged FDA Drugs");
        title("Top PAINS Patterns in 200 FDA Drugs");
        grid("on");
    else
        logInfo("  Skipping bar graph as all pattern hit counts are uniform (%d).", topHitCounts(1));
        logInfo("  Modify Section 5 to painsScreen = painsCsv; and run with all 480 patterns");
        logInfo("  to reveal the distribution of hit counts, making the bar graph significant.");
    end
end
%[text] **✏️ Try It 5 — Find the Drug with the Most PAINS Hits**
%[text] Let's investigate the drug that hits the most PAINS patterns.
%[text] Hint:
%[text]     topDrugs = sortrows(resultTbl, "PAINS\_Count", "descend");
%[text]     disp(topDrugs(1:5, \["Name", "PAINS\_Count", "First\_Alert"\]));
%[text] 
%[text] The drugs that rank high may have their "alerts" as the core reason for their efficacy.
%[text] Search for the drug name + mechanism of action to verify.
% ... (Try writing code here)
%%
%[text] ## Section 7: Save Flagged Compound Report
%[text] Save the CSV report summarizing the screening results in `result/runs/`.
runDir = makeRunDir();
outFile = fullfile(runDir, "s03_pains_report.csv");
writetable(resultTbl, outFile);

logInfo("Report saved: %s", outFile);
logInfo("Columns: Name, SMILES, PAINS_Flagged, PAINS_Count, First_Alert");
logInfo("  Extract flagged drugs with resultTbl(resultTbl.PAINS_Flagged, :).");
logInfo("  Extract clean drugs with resultTbl(~resultTbl.PAINS_Flagged, :).");

logInfo("--- S03 Complete ---");
logInfo("Summary:");
logInfo("  Custom alert panel: Screened %d patterns on %d test molecules", ...
    nAlerts, nMols);
logInfo("  PAINS screen: %d patterns on %d FDA drugs -- %d flagged (%.0f%%)", ...
    nPains, nDrugs, nFlagged, 100*nFlagged/nDrugs);
logInfo("  Key point: PAINS flag = trigger for investigation, not automatic exclusion");
%[text] **✏️ Try It 6 — Apply Custom Alerts to FDA Drugs**
%[text] Apply the custom alert panel defined in Section 1 to FDA drugs.
%[text] Hint:
%[text]     drugAlertMat = false(nDrugs, nAlerts);
%[text]     for j = 1:nAlerts
%[text]         hits\_j = emk.mol.hasSubstruct(mols\_drug(validDrug), ALERT\_SMARTS(j));
%[text]         drugAlertMat(validIdx, j) = hits\_j;
%[text]     end
%[text]     alertFlagged = sum(drugAlertMat, 2) \> 0;
%[text]     logInfo("Custom Alerts: %d/%d flagged", sum(alertFlagged), nDrugs);
%[text] 
%[text] Is the flag rate for custom alerts higher or lower than PAINS?
%[text] Which criteria do you think is more stringent?

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
