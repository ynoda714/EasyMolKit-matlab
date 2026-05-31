%[text] # S03 Answers: Structural Alerts — Detecting Hazardous Functional Groups with SMARTS
%[text] This is a reference answer for the "Try It Yourself" exercise in `s03_structure_alerts.m`.
%[text] First, run `s03_structure_alerts.m`, then use this file to check your answers.
addpath(genpath("src"));
emk.setup.initPython();
logInfo("S03 Answers: Setup Complete");

%[text] ---- Alert Panel (Same as s03 Section 1)--------------------------------
ALERT_NAMES  = ["Epoxide",      "Michael_acceptor", "Aldehyde", ...
                "Acyl_halide",  "Nitro",             "Hydrazine", ...
                "Diazo",        "Thiol"];
ALERT_SMARTS = ["[C]1CO1",      "[C]=[C]-C=O",      "[CH]=O", ...
                "[C](=O)[F,Cl]","[N+](=O)[O-]",     "[NH]-[NH2]", ...
                "[#6]=[N+]=[N-]","[SH]"];

TEST_NAMES  = ["Epichlorohydrin", "Acrolein",   "Formaldehyde", ...
               "Malondialdehyde", "Nitrobenzene", "Phenylhydrazine", ...
               "Aspirin",         "Captopril"];
TEST_SMILES = ["ClCC1CO1",          "C=CC=O",       "C=O", ...
               "O=CCC=O",          "c1ccc([N+](=O)[O-])cc1", ...
               "c1ccccc1NN",       "CC(=O)Oc1ccccc1C(=O)O", ...
               "CC(CS)C(=O)N1CCCC1C(=O)O"];

nMols   = numel(TEST_SMILES);
nAlerts = numel(ALERT_NAMES);
mols_test = cell(1, nMols);
for i = 1:nMols
    mols_test{i} = emk.mol.fromSmiles(TEST_SMILES(i));
end
flagMat = false(nMols, nAlerts);
for j = 1:nAlerts
    flagMat(:, j) = emk.mol.hasSubstruct(mols_test, ALERT_SMARTS(j))';
end
%[text] ---
%%
%[text] ## Let's Try 1: Most Alerts Molecule & Aspirin's Alerts

totalAlerts = sum(flagMat, 2);
[maxAl, iMax] = max(totalAlerts);
logInfo("Most Alerts: %s (%d items)", TEST_NAMES(iMax), maxAl);
for i = 1:nMols
    alertList = ALERT_NAMES(flagMat(i, :));
    if numel(alertList) > 0
        logInfo("  %s: %s", TEST_NAMES(i), strjoin(alertList, ", "));
    else
        logInfo("  %s: (No alerts)", TEST_NAMES(i));
    end
end

%[text] A: Acrolein (`C=CC=O`) triggers both Aldehyde AND Michael Receptor (2 items).
%[text]    Aspirin has no alerts in this panel.
%[text]    Aspirin's ester group (acetyl of phenol) is not captured by these SMARTS.
%[text]    The ester is an intentional pharmacophore that hydrolyzes in vivo to salicylic acid.
%[text]    "No alerts" does not mean "safe", just the absence of these specific patterns.
%[text]    Captopril (an ACE inhibitor) triggers thiol — expected as it chelates the zinc ion of ACE as a pharmacophore.
%%
%[text] ## Let's Try 2: Michael Addition to the Beta Carbon of Acrolein

mol_acrolein = emk.mol.fromSmiles("C=CC=O");
figure("Name", "Acrolein -- Michael + Aldehyde");
emk.viz.draw2d(mol_acrolein, Title="Acrolein — Michael + Aldehyde");

%[text] Q: Which atom does cysteine thiol attack?
%[text] A: β Carbon (C2 of the C=C double bond — Michael receptor site).
%[text]    Soft nucleophiles (thiols) prefer 1,4 addition (Michael addition) over
%[text]    1,2 addition (direct carbonyl attack) because the LUMO extends to the β position of the conjugated system.
%%
%[text] ## Let's Try 3: The First 3 PAINS Patterns

fid = fopen("data/list/pains.csv", "r");
textscan(fid, "%s", 1, "Delimiter", "\n");
C = textscan(fid, "%q%q%q%q", "Delimiter", ",");
fclose(fid);
painsCsv = table(string(C{1}), string(C{2}), string(C{3}), ...
    VariableNames=["Name", "SMARTS", "FilterSet"]);
painsCsv = painsCsv(strtrim(painsCsv.SMARTS) ~= "", :);

logInfo("The first 3 PAINS patterns:");
disp(painsCsv(1:3, ["Name","SMARTS","FilterSet"]));

%[text] A: PAINS SMARTS are complex recursive queries for substructure matching -- not for human reading.
%[text]    The FilterSet labels (A/B/C) correspond to the three types of HTS assays in the
%[text]    Baell & Holloway 2010 paper.
%%
%[text] ## Let's Try 4: Proportion of PAINS in FDA-Approved Drugs

rawDrugs = readtable("data/list/fda_drugs.csv", "TextType", "string");
nDrugs   = height(rawDrugs);

validDrug = false(nDrugs, 1);
mols_drug = cell(nDrugs, 1);
for i = 1:nDrugs
    if emk.mol.isValid(rawDrugs.SMILES(i))
        mols_drug{i} = emk.mol.fromSmiles(rawDrugs.SMILES(i));
        validDrug(i) = true;
    end
end
validIdx    = find(validDrug);
%[text] Using PAINS_A (the strictest 16 patterns)
painsScreen = painsCsv(painsCsv.FilterSet == "PAINS_A", :);
nPains      = height(painsScreen);

isPainsFlagged = false(nDrugs, 1);
painsCount     = zeros(nDrugs, 1);
firstAlert     = strings(nDrugs, 1);
logInfo("Screening %d drugs with %d PAINS_A patterns...", nDrugs, nPains);
logInfo("  (PAINS_A = strictest subset; %d in the entire database)", height(painsCsv));
for j = 1:nPains
    hits_j = emk.mol.hasSubstruct(mols_drug(validDrug), painsScreen.SMARTS(j));
    for k = 1:numel(hits_j)
        idx = validIdx(k);
        if hits_j(k)
            isPainsFlagged(idx) = true;
            painsCount(idx)     = painsCount(idx) + 1;
            if firstAlert(idx) == ""
                firstAlert(idx) = painsScreen.Name(j);
            end
        end
    end
    logProgress(j, nPains, "PAINS Screening");
end
nFlagged = sum(isPainsFlagged);
logInfo("PAINS_A flagged: %d / %d  (%.0f%%)", nFlagged, nDrugs, 100*nFlagged/nDrugs);

%[text] A: Typically, ~3-5% of approved drugs are flagged by PAINS_A (16 patterns).
%[text]    Using all 480 patterns, this rises to 15-25% (PAINS_B/C are more lenient patterns).
%[text]    Approved drugs may have reactive motifs as pharmacophores
%[text]    (e.g., thiol in captopril, beta-lactam antibiotics).
%[text]    PAINS are a trigger for investigation, not grounds for automatic rejection.
%%
%[text] ## Let's Try 5: Drugs with Most PAINS Patterns
%[text] (isPainsFlagged / painsCount / firstAlert are calculated above)
resultTbl = table(rawDrugs.Name, rawDrugs.SMILES, isPainsFlagged, painsCount, firstAlert, ...
    'VariableNames', ["Name","SMILES","PAINS_Flagged","PAINS_Count","First_Alert"]);
topDrugs = sortrows(resultTbl, "PAINS_Count", "descend");
logInfo("Top 5 PAINS Count:");
disp(topDrugs(1:5, ["Name","PAINS_Count","First_Alert"]));

%[text] A: Drugs with multiple PAINS hits tend to be complex natural products with many functional groups.
%[text]    Let's investigate the mechanism of action of the top drugs —
%[text]    the "alert" might be the key reason why the drug is effective.
%%
%[text] ## Let's Try 6: Applying Custom Alerts to FDA Drugs

%[text] Apply the alert panel from Section 1 to FDA drugs and compare the flag rate with PAINS.
%[text] (`mols_drug`, `validDrug`, `validIdx`, `nDrugs`, `nAlerts`, `ALERT_SMARTS` are pre-calculated above)

drugAlertMat = false(nDrugs, nAlerts);
for j = 1:nAlerts
    hits_j = emk.mol.hasSubstruct(mols_drug(validDrug), ALERT_SMARTS(j));
    drugAlertMat(validIdx, j) = hits_j;
end
alertFlagged = sum(drugAlertMat, 2) > 0;
nAlertFlagged = sum(alertFlagged);
logInfo("Custom Alerts: %d / %d flagged (%.0f%%)", ...
    nAlertFlagged, nDrugs, 100*nAlertFlagged/nDrugs);
logInfo("PAINS_A:      %d / %d flagged", nFlagged, nDrugs);

%[text] A: Custom alerts (8 types) identify specific chemical structure patterns,
%[text]    so the hit rate for FDA drugs differs from PAINS_A.
%[text]    Captopril (Thiol) and beta-lactam antibiotics (acyl halide) may trigger alerts,
%[text]    potentially resulting in a higher flag rate than PAINS.
%[text]    Both approaches mean "flag = trigger for investigation" and are not grounds for automatic exclusion.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
