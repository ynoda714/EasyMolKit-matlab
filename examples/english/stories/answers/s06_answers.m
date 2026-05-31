%[text] # S06 Answers: PubChem Database Search
%[text] Reference answers for the "Try it yourself" exercise in s06_pubchem_search.m.
%[text] First, run s06_pubchem_search.m, then use this file to check your answers.
%[text] Note: This file requires an internet connection (PubChem REST API).
addpath(genpath("src"));
emk.setup.initPython();
logInfo("S06 Answers: Setup complete (Internet connection required)");
%%
%[text] ## Let's Try 1: Aspirin CID & Paracetamol == Acetaminophen?

tbl_aspirin = emk.db.searchPubchem("aspirin");
logInfo("Aspirin CID: %d", tbl_aspirin.CID(1));
logInfo("Aspirin IUPAC Name: %s", tbl_aspirin.IUPACName(1));

tbl_para = emk.db.searchPubchem("paracetamol");
tbl_acet = emk.db.searchPubchem("acetaminophen");
logInfo("Paracetamol CID: %d", tbl_para.CID(1));
logInfo("Acetaminophen CID: %d", tbl_acet.CID(1));
logInfo("Same compound? %d", tbl_para.CID(1) == tbl_acet.CID(1));

%[text] A: Aspirin CID = 2244. This stable identifier has been the record for aspirin since the database began.
%[text]    IUPAC Name: "2-acetyloxybenzoic acid" -- benzoic acid with an acetyl ester at the ortho position.
%[text]    Paracetamol == Acetaminophen: Both have CID 1983 (same compound, different regional names).
%%
%[text] ## Let's Try 2: Search for Ibuprofen with SMILES and (S)-enantiomer CID

tbl_ibu = emk.db.searchPubchem("CC(C)Cc1ccc(cc1)C(C)C(=O)O", Type="smiles");
logInfo("Ibuprofen (SMILES) CID: %d", tbl_ibu.CID(1));
logInfo("PubChem canonical SMILES: %s", tbl_ibu.IsomericSMILES(1));

%[text] (S)-Ibuprofen
tbl_s_ibu = emk.db.searchPubchem("[C@@H](C(=O)O)(Cc1ccc(cc1)CC(C)C)C", Type="smiles");
logInfo("(S)-Ibuprofen CID: %d", tbl_s_ibu.CID(1));
logInfo("Same as racemic? %d", tbl_ibu.CID(1) == tbl_s_ibu.CID(1));

%[text] A: Ibuprofen (racemic) CID = 3672.
%[text]    PubChem returns its canonical SMILES (the input SMILES may have a different atom order but represent the same molecule).
%[text]    (S)-Ibuprofen is registered in PubChem as a structure with different stereochemistry,
%[text]    thus it has a different CID (PubChem distinguishes stereoisomers by default).
%%
%[text] ## Let's Try 3: Antipyretic Analgesic Panel -- Heaviest, Smallest, Sulfur-Containing Compounds

PANEL_NAMES = ["aspirin", "ibuprofen", "acetaminophen", "naproxen", "celecoxib"];
N = numel(PANEL_NAMES);

cids    = zeros(N, 1, 'uint32');
iupac   = strings(N, 1);
formulas= strings(N, 1);
mws     = zeros(N, 1);
panelSmiles = strings(N, 1);

for i = 1:N
    try
        r = emk.db.searchPubchem(PANEL_NAMES(i));
        cids(i)     = r.CID(1);
        iupac(i)    = r.IUPACName(1);
        formulas(i) = r.MolecularFormula(1);
        mws(i)      = r.MolecularWeight(1);
        panelSmiles(i) = r.IsomericSMILES(1);
    catch ME
        logWarn("  %s: Query failed -- %s", PANEL_NAMES(i), ME.message);
    end
end

panelTbl = table(PANEL_NAMES', cids, iupac, formulas, mws, panelSmiles, ...
    'VariableNames', {'Name','CID','IUPACName','MolecularFormula','MolecularWeight','SMILES'});

%[text] Heaviest
[~, iHeavy] = max(panelTbl.MolecularWeight);
logInfo("Heaviest: %s (MW=%.1f)", panelTbl.Name(iHeavy), panelTbl.MolecularWeight(iHeavy));

%[text] Sulfur-Containing
hasSulfur = contains(panelTbl.MolecularFormula, "S");
logInfo("Sulfur-Containing: %s", strjoin(panelTbl.Name(hasSulfur), ", "));

%[text] A: Heaviest: Celecoxib (~381 Da). The latest selective COX-2 inhibitors are generally larger and more complex than traditional NSAIDs.
%[text]    Aspirin and acetaminophen both have MW < 180 Da but completely different mechanisms of action --
%[text]    MW does not predict mechanism.
%[text]    Sulfur-Containing: Celecoxib (C17H14F3N3O2S) -- the sulfonamide (SO2NH2) is key to COX-2 selectivity.
%%
%[text] ## Let's Try 4: Complete Property Table, Ro5 Check, TPSA Ranking
%[text] Construct RDKit descriptors
descNames = ["MolWt","LogP","TPSA","NumHDonors","NumHAcceptors","NumRotatableBonds","RingCount"];
mwArr  = zeros(N,1);  lpArr  = zeros(N,1);  tpArr  = zeros(N,1);
hdArr  = zeros(N,1);  haArr  = zeros(N,1);  rbArr  = zeros(N,1); rcArr  = zeros(N,1);

for i = 1:N
    if panelSmiles(i) ~= ""
        mol = emk.mol.fromSmiles(panelSmiles(i));
        d   = emk.descriptor.calculate(mol, descNames);
        mwArr(i) = d.MolWt;  lpArr(i) = d.LogP;  tpArr(i) = d.TPSA;
        hdArr(i) = d.NumHDonors;  haArr(i) = d.NumHAcceptors;
        rbArr(i) = d.NumRotatableBonds;  rcArr(i) = d.RingCount;
    end
end

fullTbl = table(PANEL_NAMES', cids, formulas, mwArr, lpArr, tpArr, ...
    hdArr, haArr, rbArr, rcArr, ...
    'VariableNames', {'Name','CID','MolecularFormula','MolecularWeight', ...
                      'LogP','TPSA','NumHDonors','NumHAcceptors', ...
                      'NumRotatableBonds','RingCount'});

%[text] Ro5 Check
roTbl  = table(fullTbl.MolecularWeight, fullTbl.LogP, fullTbl.NumHDonors, fullTbl.NumHAcceptors, ...
    'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
result = emk.filter.lipinski(roTbl);
logInfo("Ro5 Results:");
for i = 1:N
    logInfo("  %-15s  Pass=%d  Violations=%d", fullTbl.Name(i), result.Pass_Ro5(i), result.Violations_Ro5(i));
end

%[text] Highest TPSA, Lowest LogP
[~, iMaxTPSA] = max(fullTbl.TPSA);
[~, iMinLogP] = min(fullTbl.LogP);
logInfo("Highest TPSA: %s (%.1f A^2)", fullTbl.Name(iMaxTPSA), fullTbl.TPSA(iMaxTPSA));
logInfo("Lowest LogP: %s (%.2f)", fullTbl.Name(iMinLogP), fullTbl.LogP(iMinLogP));

%[text] A: All 5 pass Ro5 (all widely used oral analgesics).
%[text]    Highest TPSA: Celecoxib (~78 A^2, sulfonamide + pyrazole).
%[text]    Lowest LogP: Aspirin (~1.31) -- most hydrophilic in the panel.
%%
%[text] ## Let's Try 5: Visual Investigation of Structure -- Celecoxib's CF3

mol_cel = emk.mol.fromSmiles(panelSmiles(panelTbl.Name == "celecoxib"));
emk.viz.draw2d(mol_cel, Title="Celecoxib (COX-2 Selective Inhibitor)", ...
    Width=400, Height=400);

%[text] Q: Why introduce CF3?
%[text] A: CF3 (trifluoromethyl) increases lipophilicity (LogP +0.5 per group),
%[text]    improves metabolic stability (prevents CYP oxidation at that site),
%[text]    and can sterically mimic a methyl group without undergoing metabolism.
%[text]    One of the most frequently used substitution patterns in medicinal chemistry.
%%
%[text] ## Let's Try 6: CNS Profile & Celecoxib TPSA & MW Bar Graph
%[text] CNS Target Profile: TPSA < 60, LogP 1-3
cnsFit = fullTbl.TPSA < 60 & fullTbl.LogP >= 1 & fullTbl.LogP <= 3;
logInfo("CNS Profile (TPSA<60, LogP 1-3): %s", strjoin(fullTbl.Name(cnsFit), ", "));

%[text] A: Ibuprofen (TPSA ~37, LogP ~3.1) and Aspirin (TPSA ~64, LogP ~1.31) are
%[text]    closest to the CNS profile. Neither is used as a major CNS drug, but
%[text]    structurally they can permeate the BBB.
%[text]    Celecoxib has a high TPSA ~78 A^2, but is acceptable as it targets peripheral COX-2 at inflammation sites.
%[text]
%[text] MW Bar Graph (with Lipinski Limit)
figure("Name", "Antipyretic Analgesic Panel -- MW and Lipinski Limit", "Position", [100 100 560 380]);
bar(fullTbl.MolecularWeight);
xticks(1:N);  xticklabels(fullTbl.Name);  xtickangle(20);
yline(500, "--r", "Lipinski MW=500 Limit");
ylabel("Molecular Weight (Da)");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
