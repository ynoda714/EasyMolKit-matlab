%[text] # S02 Answers: Drug Filter -- Lipinski's Rule of Five
%[text] Reference answers for the "Try It" exercise in s02_drug_filter_lipinski.m.
%[text] First, run s02_drug_filter_lipinski.m, then use this file to check your answers.
addpath(genpath("src"));
emk.setup.initPython();
logInfo("S02 Answers: Setup complete");
%%
%[text] ## Let's Try 1: Ro5 Check for Aspirin and Ibuprofen

mol_asp  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
mol_ibu  = emk.mol.fromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O");

for pair = {{"Aspirin",    mol_asp}, {"Ibuprofen", mol_ibu}}
    name = pair{1}{1};  mol = pair{1}{2};
    d = emk.descriptor.calculate(mol, ["MolWt","LogP","NumHDonors","NumHAcceptors"]);
    logInfo("%s: MW=%.1f  LogP=%.2f  HBD=%d  HBA=%d", ...
        name, d.MolWt, d.LogP, d.NumHDonors, d.NumHAcceptors);
end

%[text] A: Aspirin   -- MW ~180, LogP ~1.3, HBD=1, HBA=3 -- Violations 0 (all pass).
%[text]    Ibuprofen -- MW ~206, LogP ~3.1, HBD=1, HBA=1 -- Violations 0 (all pass).
%[text]    Both are small, easily absorbed oral drugs -- as predicted by Ro5.
%%
%[text] ## Let's Try 2: What is the heaviest drug in the FDA dataset?

rawTbl = readtable("data/list/fda_drugs.csv", "TextType", "string");

%[text] Rename columns to match the expected names for `emk.filter.lipinski`
rawTbl = renamevars(rawTbl, ...
    ["MolecularWeight", "ALogP",  "HBondDonors",  "HBondAcceptors"], ...
    ["MolWt",           "LogP",   "NumHDonors",   "NumHAcceptors"]);

%[text] Convert to numeric columns
rawTbl.MolWt         = double(rawTbl.MolWt);
rawTbl.LogP          = double(rawTbl.LogP);
rawTbl.NumHDonors    = double(rawTbl.NumHDonors);
rawTbl.NumHAcceptors = double(rawTbl.NumHAcceptors);

[maxMW, idx] = max(rawTbl.MolWt);
logInfo("Heaviest drug: %s  (MW = %.1f Da)", rawTbl.Name(idx), maxMW);
logInfo("SMILES: %s", rawTbl.SMILES(idx));

%[text] A: Amphotericin B (MW ~924 Da) -- A polyene antifungal drug derived from Streptomyces nodosus.
%[text]    Administered intravenously, not orally. Consistent with Ro5 violation.
%%
%[text] ## Let's Try 3: Which Ro5 criterion is most frequently violated?

drugTbl = emk.filter.lipinski(rawTbl);

vMW  = sum(rawTbl.MolWt           > 500);
vLP  = sum(rawTbl.LogP             > 5);
vHBD = sum(rawTbl.NumHDonors      > 5);
vHBA = sum(rawTbl.NumHAcceptors   > 10);

logInfo("Number of violations:");
logInfo("  MW  > 500: %d drugs", vMW);
logInfo("  LogP > 5 : %d drugs", vLP);
logInfo("  HBD > 5  : %d drugs", vHBD);
logInfo("  HBA > 10 : %d drugs", vHBA);

%[text] A: MW violations are usually the most common (large natural product-derived drugs), followed by HBA violations.
%[text]    Most drugs have few NH/OH groups, so HBD > 5 is the rarest.
%%
%[text] ## Let's Try 4: Drugs with Exactly 1 Violation; Number Passing Relaxed Filter

violCounts = histcounts(drugTbl.Violations_Ro5, 0:6);
logInfo("Drugs with 0 violations: %d", violCounts(1));
logInfo("Drugs with 1 violation: %d", violCounts(2));
logInfo("Drugs with 2 violations: %d", violCounts(3));
logInfo("Drugs with 3 violations: %d", violCounts(4));

drugTblRelaxed = emk.filter.lipinski(rawTbl, MaxViolations=1);
nPassStrict  = sum(drugTbl.Pass_Ro5);
nPassRelaxed = sum(drugTblRelaxed.Pass_Ro5);
logInfo("Strict (0 violations): %d / %d pass", nPassStrict, height(rawTbl));
logInfo("Relaxed (1 violation): %d / %d pass", nPassRelaxed, height(rawTbl));
logInfo("Additional recovered drugs: %d", nPassRelaxed - nPassStrict);

%[text] A: Typically ~14 drugs have exactly 1 violation.
%[text]    Relaxing to MaxViolations=1 recovers these (many are large macrolide antibiotics or
%[text]    clinically proven natural products).
%%
%[text] ## Let's Try 5: Chemical Space Scatter Plot -- MW vs LogP

figure("Name", "Chemical Space: MW vs LogP", "Position", [100 100 560 440]);
passIdx = drugTbl.Pass_Ro5;
scatter(rawTbl.LogP(passIdx),  rawTbl.MolWt(passIdx),  40, "b", "filled", ...
    "DisplayName", "Ro5 Pass");
hold on;
scatter(rawTbl.LogP(~passIdx), rawTbl.MolWt(~passIdx), 60, "r", "filled", ...
    "DisplayName", "Ro5 Fail");
xline(5,   "--k", "LogP=5");
yline(500, "--k", "MW=500");
xlabel("LogP");  ylabel("MW (Da)");
title("Chemical Space of FDA Approved Drugs");
legend("Location", "northwest");  grid("on");

%[text] A: Yes -- Most FDA approved drugs are concentrated in the lower left quadrant (MW < 500, LogP < 5).
%[text]    Drugs outside this quadrant are typically injectables, topicals, or prodrugs.
%[text]    Some Ro5 pass drugs cluster around LogP ~4-5, MW ~450-490 -- drug-like but near the boundary.
%%
%[text] ## Let's Try 6: Visualizing Rule-Breaking Drugs
%[text] Example of a non-compliant drug -- Rifampicin is typical (if available)
rifIdx = find(strcmpi(rawTbl.Name, "RIFAMPICIN"), 1);
if isempty(rifIdx)
    [~, rifIdx] = max(rawTbl.MolWt);   % Fallback: Heaviest drug
end
mol_big = emk.mol.fromSmiles(rawTbl.SMILES(rifIdx));
figure("Name", "Ro5 Violating Drug: " + rawTbl.Name(rifIdx), "Position", [100 100 440 380]);
emk.viz.draw2d(mol_big, Title=rawTbl.Name(rifIdx) + " (Ro5 Violation)");
logInfo("Rule-breaking drug: %s  MW=%.1f  Number of violations=%d", ...
    rawTbl.Name(rifIdx), rawTbl.MolWt(rifIdx), drugTbl.Violations_Ro5(rifIdx));

%[text] A: Rifampicin has a large macrocycle (complex aromatic system).
%[text]    Counting donors/acceptors reveals many OH and NH groups.
%[text]    Yet it is orally absorbed -- Ro5 is a guideline, not an absolute rule.
%%
%[text] ## Let's Try 7: MaxViolations=2 Filter

drugTbl2 = emk.filter.lipinski(rawTbl, MaxViolations=2);
%[text] Drugs that pass with MaxViolations=2 but fail with MaxViolations=1:
newPass = drugTbl2.Pass_Ro5 & ~drugTblRelaxed.Pass_Ro5;
logInfo("Drugs with exactly 2 violations: %d", sum(newPass));
disp(rawTbl(newPass, ["Name", "MolWt", "LogP", "NumHDonors", "NumHAcceptors"]));

%[text] A: Drugs with 2 violations include large macrolides (Azithromycin, Erythromycin),
%[text]    taxanes (Paclitaxel), and immunosuppressants (Tacrolimus).
%[text]    All are administered intravenously -- reflecting their natural product origins and complex pharmacology.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
