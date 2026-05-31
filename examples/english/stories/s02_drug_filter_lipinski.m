%[text] # S02: Drug Filter — Lipinski's "Rule of Five"
%[text] EasyMolKit Application Story — Layer 2
%[text] 
%[text] The "drugs" we take need to not only reach targets (such as proteins) in the body but also be absorbed from the intestines, enter the bloodstream, and circulate properly throughout the body.
%[text] No matter how potent a molecule is found in the laboratory, it cannot become a drug if it is not absorbed by the body. In 1997, Christopher Lipinski and colleagues at the pharmaceutical company Pfizer discovered four common physicochemical characteristics of many oral drugs. This is the famous "Lipinski's Rule of 5."
%[text] In this script, we apply Lipinski's rules to a database of 200 FDA-approved drugs obtained from ChEMBL to determine which molecules possess "drug-like oral absorption" using MATLAB's data processing techniques.
%[text] ## Learning Objectives
%[text] - Use the descriptor calculation function (`emk.descriptor.calculate`) to comprehensively obtain molecular properties
%[text] - Implement the criteria of Lipinski's Rule of 5 in code
%[text] - Use MATLAB's conditional judgment and table operations (logical masking) to extract molecules that meet the criteria
%[text] - Calculate a unique "pass score" and screen the entire database \
%[text] ## Prerequisites
%[text] - Completion of F02 (Property Calculation)
%[text] - RDKit installed (execute `emk.setup.install()` once)
%[text] - No additional Toolbox required (works with MATLAB alone) \
%[text] **Estimated Time**: 15-20 minutes | Execution Method: Run each section one by one with Ctrl+Enter
%[text] **Data**
%[text] - `data/list/fda_drugs.csv` — 200 FDA-approved drugs (ChEMBL, CC-BY-SA 3.0) \
%[text] **References**
%[text] - Lipinski CA et al. (1997) *Adv Drug Deliv Rev* 23:3-25. doi:10.1016/S0169-409X(96)00423-1 [requires institutional access]\
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
logInfo("S02: Setup complete");
%%
%[text] ## Section 1: What is the Rule of 5 (Ro5)?
%[text] The Rule of 5 is a guideline to determine if a compound is likely to be well-absorbed when taken orally by permeating biological membranes. It is called the **Rule of 5 (Ro5)** because all criteria involve the number "5" or its multiples.
%[text] **[The Four Rules of Lipinski]**
%[text] 1. Molecular Weight: $&dollar&;\\leq 500 \\text{ g/mol}&dollar&;$ (If the molecule is too large, it cannot pass through the membrane)
%[text] 2. Lipophilicity (LogP): $&dollar&;5&dollar&;$ or less (If it dissolves too much in oil, it remains in the membrane and cannot enter the bloodstream)
%[text] 3. Number of Hydrogen Bond Donors: $&dollar&;5&dollar&;$ or less (If it binds too strongly with water molecules, it dislikes lipid membranes)
%[text] 4. Number of Hydrogen Bond Acceptors: $&dollar&;10&dollar&;$ or less (For the same reason, if too many, it cannot permeate the membrane) \
%[text] Generally, if a compound meets three or more of these four criteria, it is considered a "well-absorbed molecule (a drug-like candidate)."
%[text] Let's use EasyMolKit to calculate these values for Aspirin, the world's most consumed oral analgesic.
ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O";
ASPIRIN_NAME   = "Aspirin";

mol_aspirin = emk.mol.fromSmiles(ASPIRIN_SMILES);
desc_aspirin = emk.descriptor.calculate(mol_aspirin, ...
    ["MolWt", "LogP", "NumHDonors", "NumHAcceptors", "TPSA"]);

logInfo("--- %s (Acetylsalicylic Acid) ---", ASPIRIN_NAME);
logInfo("  MW              : %.2f Da  (Rule: <= 500)", desc_aspirin.MolWt);
logInfo("  LogP            : %.2f     (Rule: <= 5  )", desc_aspirin.LogP);
logInfo("  H-Bond Donors   : %d       (Rule: <= 5  )", desc_aspirin.NumHDonors);
logInfo("  H-Bond Acceptors: %d       (Rule: <= 10 )", desc_aspirin.NumHAcceptors);
logInfo("  TPSA            : %.1f A^2 (Oral < 130 A^2 guideline)", desc_aspirin.TPSA);

figure("Name", "Aspirin", "Position", [100 100 440 380]);
emk.viz.draw2d(mol_aspirin, Title="Aspirin (Acetylsalicylic Acid)");

%[text] **✏️ Try It 1 — Check with another analgesic**
%[text] Look at the numbers above and count how many of the 4 Ro5 criteria Aspirin violates.
%[text] Expected: 0 violations. MW ~180, LogP ~1.3, HBD=1, HBA=3 — all pass.
%[text] Next, try Ibuprofen (`"CC(C)Cc1ccc(cc1)C(C)C(=O)O"`).
%[text]   mol\_ibu  = `emk.mol.fromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")`;
%[text]   desc\_ibu = `emk.descriptor.calculate(mol_ibu, ["MolWt","LogP","NumHDonors","NumHAcceptors"])`;
%[text] Expected: MW ~206, LogP ~3.1, HBD=1, HBA=1 — both pass Ro5.
%[text] Q: Ibuprofen has a higher LogP than Aspirin. Where do you think the chemical structure differs?
% ... (Try writing the code here)
%%
%[text] ## Section 2: Load FDA Drug Database
%[text] Load data for 200 FDA-approved drugs obtained from ChEMBL (CC-BY-SA 3.0).
%[text] The CSV includes precomputed descriptors (e.g., ALogP) from ChEMBL.
%[text] Rename to match the column names expected by `emk.filter.lipinski`.
dataFile = fullfile(projectRoot, "data", "list", "fda_drugs.csv");
rawTbl = readtable(dataFile, "TextType", "string");
logInfo("Loaded %d FDA-approved drugs from ChEMBL.", height(rawTbl));
%[text] Rename to match the required column names for `emk.filter.lipinski`
drugTbl = renamevars(rawTbl, ...
    ["MolecularWeight", "ALogP",  "HBondDonors",  "HBondAcceptors"], ...
    ["MolWt",           "LogP",   "NumHDonors",   "NumHAcceptors"]);
logInfo("Descriptor summary (before filtering):");
logInfo("  MW   -- Min: %5.1f  Median: %5.1f  Max: %6.1f  Da", ...
    min(drugTbl.MolWt), median(drugTbl.MolWt), max(drugTbl.MolWt));
logInfo("  LogP -- Min: %5.2f  Median: %5.2f  Max: %5.2f", ...
    min(drugTbl.LogP), median(drugTbl.LogP), max(drugTbl.LogP));
logInfo("  HBD  -- Min: %d     Median: %g    Max: %d", ...
    min(drugTbl.NumHDonors), median(drugTbl.NumHDonors), max(drugTbl.NumHDonors));
logInfo("  HBA  -- Min: %d     Median: %g    Max: %d", ...
    min(drugTbl.NumHAcceptors), median(drugTbl.NumHAcceptors), max(drugTbl.NumHAcceptors));
%[text] **✏️ Try It 2 — Find the Heaviest Drug**
%[text] Among the 200 drugs loaded, find the one with the largest molecular weight. Use MATLAB's "extract rows that meet conditions" feature to check the name and molecular weight of that drug.
%[text] Hint: Use `drugTbl(drugTbl.MolWt == max(drugTbl.MolWt), :)` to extract the heaviest row.
%[text] Expected: A large natural product-derived drug (e.g., polyene antifungal) with MW over 900 Da.
%[text] Q: Is the drug administered orally or by injection? Check if it matches the number of Ro5 violations.
% ... (Write your code here)
%%
%[text] ## Section 3: Apply Lipinski's Rule of Five
%[text] ### emk.filter.lipinski
%[text] `emk.filter.lipinski` adds 2 columns to the descriptor table.
%[text] - `Pass_Ro5` (logical) — `true` if the number of violations is `MaxViolations` or less
%[text] - `Violations_Ro5` (double) — Number of criteria violated (0–4) \
%[text] Let's first apply a strict filter with zero violations (default).
drugTbl = emk.filter.lipinski(drugTbl);          % MaxViolations=0 (default)

nTotal  = height(drugTbl);
nPass   = sum(drugTbl.Pass_Ro5);
nFail   = nTotal - nPass;

logInfo("Ro5 filter results (strict, MaxViolations=0):");
logInfo("  Passed : %d / %d  (%.0f%%)", nPass,  nTotal, 100*nPass/nTotal);
logInfo("  Failed : %d / %d  (%.0f%%)", nFail,  nTotal, 100*nFail/nTotal);

%[text] Count the number of violations for each criterion
vMW  = sum(drugTbl.MolWt          > 500);
vLP  = sum(drugTbl.LogP           > 5  );
vHBD = sum(drugTbl.NumHDonors     > 5  );
vHBA = sum(drugTbl.NumHAcceptors  > 10 );

logInfo("Number of violations for each criterion:");
logInfo("  MW  > 500: %d cases", vMW);
logInfo("  LogP > 5 : %d cases", vLP);
logInfo("  HBD > 5  : %d cases", vHBD);
logInfo("  HBA > 10 : %d cases", vHBA);

%[text] **✏️ Try It 3 — Which criterion is most frequently violated?**
%[text] Compare the values of vMW, vLP, vHBD, and vHBA displayed above.
%[text] Q: Which Ro5 criterion has the most violations in this dataset?
%[text] Expectation: Drugs derived from natural products tend to be large and highly polar, leading to more violations in MW and HBA.
% ... (Try writing code here)
%%
%[text] ## Section 4: Distribution of Violation Counts
%[text] Let's display specific numbers to identify why the dropped molecules violated Lipinski's rule.
violCounts = zeros(1, 5);
for v = 0:4
    violCounts(v+1) = sum(drugTbl.Violations_Ro5 == v);
end

figure("Name", "Distribution of Ro5 Violation Counts", "Color", "white", "Position", [100 100 560 400]);
bar(0:4, violCounts, "FaceColor", [0.2 0.6 0.9]);
xlabel("Number of Ro5 Violations");
ylabel("Number of FDA Drugs");
title("Lipinski Ro5 -- Violation Count Distribution (200 FDA Drugs)");
xticks(0:4);
xticklabels({"0 (Pass)", "1", "2", "3", "4 (All Violations)"});
grid("on");
for v = 0:4
    if violCounts(v+1) > 0
        text(v, violCounts(v+1) + 0.5, sprintf("%d", violCounts(v+1)), ...
            "HorizontalAlignment", "center", "FontWeight", "bold");
    end
end
logInfo("Plotted the distribution of violation counts.");
%[text] **✏️ Try It 4 — Experiment with Relaxed Filters**
%[text] Look at the bar graph and check how many drugs have exactly 1 violation.
%[text] Hint: `violCounts(2)` gives the number with 1 violation.
%[text] Q: If you relax to `MaxViolations=1`, how many additional drugs will pass?
%[text] Hint: You can find the total passing with relaxed filters using `nPass + violCounts(2)`.
% ... (Try writing code here)
%%
%[text] ## Section 5: Visualization of Chemical Space -- MW vs LogP
%[text] ### What is Drug-like Chemical Space?
%[text] Plotting MW (y-axis) and LogP (x-axis) reveals a "drug-like" region enclosed by the Ro5 boundaries. Most oral drugs cluster below the MW=500 Da line and to the left of the LogP=5 line. Drugs outside these boundaries are often injectables, topicals, or prodrugs.
figure("Name", "MW vs LogP (Drug Space)", "Color", "white", "Position", [100 100 560 440]);
hold on;
%[text] Passed drugs (blue circles)
passIdx = drugTbl.Pass_Ro5;
scatter(drugTbl.LogP(passIdx),  drugTbl.MolWt(passIdx),  40, ...
    [0.15 0.55 0.85], "o", "filled", "MarkerFaceAlpha", 0.6);
%[text] Failed drugs (orange triangles)
scatter(drugTbl.LogP(~passIdx), drugTbl.MolWt(~passIdx), 60, ...
    [0.9 0.45 0.1], "^", "filled", "MarkerFaceAlpha", 0.7);
%[text] Ro5 boundaries
xLim = xlim;  yLim = ylim;
plot([5 5], [0 max(yLim(2), 600)], "r--", "LineWidth", 1.5);   % LogP = 5
yline(500, "r--", "LineWidth", 1.5);                           % MW = 500
%[text] Annotation for Aspirin
aspIdx = find(strcmpi(drugTbl.Name, "ASPIRIN"), 1);
if ~isempty(aspIdx)
    scatter(drugTbl.LogP(aspIdx), drugTbl.MolWt(aspIdx), 120, ...
        "g", "p", "filled");
    text(drugTbl.LogP(aspIdx) + 0.1, drugTbl.MolWt(aspIdx) + 10, ...
        "Aspirin", "FontSize", 9, "Color", [0 0.5 0]);
end

hold off;
legend({"Passed Ro5", "Failed Ro5", "LogP=5 Limit", "MW=500 Limit", "Aspirin"}, ...
    "Location", "northwest");
xlabel("LogP  (Lipophilicity)");
ylabel("Molecular Weight (Da)");
title("FDA Drugs: MW vs LogP — Lipinski Space");
grid("on");

logInfo("Displayed chemical space scatter plot.");
logInfo("  Blue region (MW<500, LogP<5): Typical oral drug-like space");
logInfo("  Drugs outside boundaries: Often injectables, topicals, or prodrugs");
%[text] **✏️ Try It 5 — Read the Scatter Plot**
%[text] Observe the generated scatter plot.
%[text] Q1: Do many blue points (passed Ro5) cluster in the lower left quadrant (MW<500, LogP<5)?
%[text] Q2: Do you see a cluster of passing drugs near LogP~5, MW~480? (Just passing Ro5 line)
%[text] Q3: What kind of drugs do you think the orange points (failed Ro5) in the upper right are?
% ... (Try writing code here)
%%
%[text] ## Section 6: Drugs that Break the Rules -- Ro5 Exceptions
%[text] 
%[text] ### Ro5 is a guideline, not an absolute rule
%[text] Ro5 is a rule for small molecule oral drugs absorbed by passive diffusion.
%[text] Some approved drugs are orally absorbed despite violating Ro5 for the following reasons:
%[text] - Absorbed by active transport (carrier proteins)
%[text] - Administered as prodrugs and converted to active forms in the body \
%[text] Famous examples:
%[text] - Cyclosporin A — MW ~1202 Da, 3 violations (immunosuppressant)
%[text] - Rifampicin — MW ~823 Da, 2 violations (antibiotic)
%[text] - Azithromycin — MW ~749 Da, 2 violations (macrolide antibiotic)
%[text] - Atorvastatin — MW ~559 Da, 1 violation (cholesterol drug) \
failTbl = drugTbl(~drugTbl.Pass_Ro5, :);
logInfo("--- %d FDA drugs failed strict Ro5 ---", height(failTbl));
%[text] Sort by number of violations (descending), then by MW
failTbl = sortrows(failTbl, ["Violations_Ro5", "MolWt"], ["descend", "descend"]);
%[text] Display the top 10 worst offenders
nShow = min(10, height(failTbl));
logInfo("Top %d worst Ro5 offenders (approved oral/systemic drugs):", nShow);
disp(failTbl(1:nShow, ["Name","MolWt","LogP","NumHDonors","NumHAcceptors","Violations_Ro5"]));
%[text] **✏️ Try It 6 — Examine the structures of Ro5 violating drugs**
%[text] Choose one drug from the list above and draw its structure.
%[text] Hint:
%[text]   idx = find(strcmpi(drugTbl.Name, "RIFAMPICIN"), 1);   % You can change to any drug name
%[text]   mol\_big = `emk.mol.fromSmiles(drugTbl.SMILES(idx))`;
%[text]   figure("Name", drugTbl.Name(idx));
%[text]  `emk.viz.draw2d(mol_big, Title=drugTbl.Name(idx))`;
%[text] How many hydrogen bond donors (OH, NH) and acceptors (C=O, N) can you see in the structure?
%[text] Check if it matches the Ro5 violation count.
% ... (Try writing your code here)
%%
%[text] ## Section 7: Relaxed Filter (MaxViolations=1)
%[text] ### "Rule of Five + α"
%[text] In many real-world drug discovery programs, `MaxViolations=1` (allowing one criterion violation) is used.
%[text] The reasons are as follows:
%[text] - Ro5 is based on a specific dataset of passively diffused drugs
%[text] - Active transport can compensate for high polarity or large size
%[text] - One violation can sometimes be corrected in the optimization step \
drugTblRelaxed = emk.filter.lipinski(drugTbl, MaxViolations=1);
nPassRelaxed   = sum(drugTblRelaxed.Pass_Ro5);

logInfo("Relaxed Filter (MaxViolations=1): %d / %d passed  (%.0f%%)", ...
    nPassRelaxed, nTotal, 100*nPassRelaxed/nTotal);
logInfo("Additional drugs passing: %d types", nPassRelaxed - nPass);
%[text] Identify newly passing drugs (strict filter fail → relaxed filter pass)
isNewPass = ~drugTbl.Pass_Ro5 & drugTblRelaxed.Pass_Ro5;
newPassTbl = drugTbl(isNewPass, :);
newPassTbl = sortrows(newPassTbl, "MolWt", "descend");

logInfo("Sample of newly passing drugs (each with 1 violation):");
nSample = min(5, height(newPassTbl));
disp(newPassTbl(1:nSample, ["Name","MolWt","LogP","NumHDonors","NumHAcceptors","Violations_Ro5"]));

logInfo("--- S02 Complete ---");
logInfo("Summary:");
logInfo("  Strict Ro5 (0 violations): %d / %d passed (%.0f%%)", ...
    nPass, nTotal, 100*nPass/nTotal);
logInfo("  Relaxed Ro5 (1 violation): %d / %d passed (%.0f%%)", ...
    nPassRelaxed, nTotal, 100*nPassRelaxed/nTotal);
logInfo("  Most violated criteria: Compare vMW, vLP, vHBD, vHBA in Sections 3-4");
%[text] **✏️ Try It 7 — Extend to MaxViolations=2**
%[text] Reapply the filter with `MaxViolations=2` (allowing up to 2 criteria) and check the newly passing drugs.
%[text] Hint:
%[text]   tbl2 = `emk.filter.lipinski(drugTbl, MaxViolations=2)`;
%[text]   pass2 = tbl2(tbl2.Pass\_Ro5 & ~drugTblRelaxed.Pass\_Ro5, :);
%[text]   disp(pass2(:, \["Name","MolWt","LogP","NumHDonors","NumHAcceptors","Violations\_Ro5"\]))
%[text] Q: Do the additional drugs include large macrolide antibiotics or immunosuppressants?
% ... (Try writing code here)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
