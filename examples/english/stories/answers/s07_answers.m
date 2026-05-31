%[text] # S07 Answers: ChEMBL Bioactivity Data
%[text] Reference answers for the "Try it yourself" exercise in s07_chembl_activity.m.
%[text] First, run s07_chembl_activity.m, then use this file to check your answers.
%[text] Note: This file requires an internet connection (ChEMBL REST API).
addpath(genpath("src"));
emk.setup.initPython();
set(0, "DefaultFigureWindowStyle", "normal");   % Display figures in pop-up windows
logInfo("S07 Answers: Setup complete (Internet connection required)");

%[text] ---- Reproduce the shared state of s07 ----------------------------------------------
ERLOTINIB_SMILES = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1";
ERLOTINIB_NAME   = "Erlotinib";
EGFR_CHEMBL_ID   = "CHEMBL203";
IC50_CUTOFF_NM   = 100;

actTbl    = emk.db.getChemblActivity(EGFR_CHEMBL_ID, ActivityType="IC50", MaxRows=50);
actTbl.pIC50 = 9 - log10(actTbl.Value_nM);
activeTbl = sortrows(actTbl(actTbl.Value_nM <= IC50_CUTOFF_NM, :), "Value_nM", "ascend");

%[text] Convert SMILES of active compounds
validMols = {};
validIdx  = [];
for i = 1:height(activeTbl)
    smi = activeTbl.SMILES(i);
    if ~emk.mol.isValid(smi); continue; end
    validMols{end+1} = emk.mol.fromSmiles(smi); %#ok<AGROW>
    validIdx(end+1)  = i;                        %#ok<AGROW>
end

mol_erl = emk.mol.fromSmiles(ERLOTINIB_SMILES);
fp_erl  = emk.fingerprint.morgan(mol_erl, Radius=2, NBits=2048);

dbFps = cell(1, numel(validMols));
for i = 1:numel(validMols)
    dbFps{i} = emk.fingerprint.morgan(validMols{i}, Radius=2, NBits=2048);
end
rankResult = emk.similarity.rankBy(fp_erl, dbFps, Inf, Metric="tanimoto");
metaTbl    = activeTbl(validIdx, :);
%[text] ---
%%
%[text] ## Let's Try 1: Definition of IC50 and Structure of Erlotinib
%[text] IC50 = Half maximal inhibitory concentration (units: nM or uM).
%[text] Lower IC50 = More potent (achieves 50% inhibition at lower concentration).
logInfo("IC50 = Half maximal inhibitory concentration");
logInfo("  Lower IC50 = More potent inhibitor");
logInfo("  Erlotinib (wild-type EGFR): IC50 ~2 nM  -->  pIC50 = %.2f", 9 - log10(2));

%[text] Draw Erlotinib and identify pharmacophore elements
emk.viz.draw2d(mol_erl, Title="Erlotinib -- EGFR Inhibitor (1st Generation)");

%[text] A1: IC50 = Half maximal inhibitory concentration (nM / uM).
%[text]     Lower IC50 is better: Achieves 50% inhibition at a lower dose.
%[text]
%[text] A2: Structural features of Erlotinib:
%[text]   (a) Alkyne (-C#C-): Binds to the ribose pocket behind the ATP site.
%[text]       Provides shape complementarity and hydrophobic contact.
%[text]   (b) Quinazoline core (bicyclic, 2 nitrogens): Hinge-binding motif.
%[text]       N1 of quinazoline donates/accepts H-bond with kinase hinge Met769.
%[text]   (c) Two -OCCOC- chains (C6 and C7): Improve water solubility and
%[text]       contribute to selectivity through contact with the gatekeeper region.
%%
%[text] ## Let's Try 2: ChEMBL Target ID & Multiple Species

logInfo("EGFR (Homo sapiens): %s", EGFR_CHEMBL_ID);

%[text] COX-2 Search -- Try multiple queries as ChEMBL preferred names may vary
cox2Queries = ["Cyclooxygenase-2", "cyclooxygenase 2", "prostaglandin"];
humanCox2 = table();
for qi = 1:numel(cox2Queries)
    cox2Tbl = emk.db.searchChemblTarget(cox2Queries(qi), MaxRows=10);
    humanCox2 = cox2Tbl(cox2Tbl.Organism == "Homo sapiens", :);
    if ~isempty(humanCox2); break; end
end
if ~isempty(humanCox2)
    logInfo("COX-2 (Homo sapiens) ChEMBL ID: %s", humanCox2.TargetChEMBLID(1));
else
    logInfo("COX-2 human target not found by name search.");
    logInfo("  Direct reference: CHEMBL230 (known stable ID for human COX-2)");
end

%[text] A: EGFR ChEMBL ID = CHEMBL203.
%[text]    Multiple EGFR entries exist because ChEMBL records targets for each species
%[text]    (human, mouse, etc.) separately. The human (Homo sapiens) entry
%[text]    has the most bioactivity data.
%[text]    COX-2 (Homo sapiens) expected value: CHEMBL230.
%%
%[text] ## Let's Try 3: Number of Compounds Below 100 nM, Strongest Activity, Ki Comparison

nBelow100 = sum(actTbl.Value_nM <= 100);
logInfo("Compounds with IC50 < 100 nM: %d / %d", nBelow100, height(actTbl));

[minIC50, minRow] = min(actTbl.Value_nM);
logInfo("Strongest active compound: %s  IC50=%.2f nM", ...
    actTbl.MoleculeChEMBLID(minRow), minIC50);

%[text] Ki values for EGFR
kiTbl = emk.db.getChemblActivity(EGFR_CHEMBL_ID, ActivityType="Ki", MaxRows=25);
logInfo("Retrieved Ki records: %d", height(kiTbl));
if height(kiTbl) > 0
    logInfo("  Ki range: %.1f - %.1f nM", min(kiTbl.Value_nM), max(kiTbl.Value_nM));
end

%[text] A: The strongest active compound in a 50-entry dataset is usually in the sub-nM range.
%[text]    ChEMBL returns results in random order -- "strongest activity" changes with each run without explicit sorting.
%[text]    Ki (inhibition constant) is a thermodynamic equilibrium indicator independent of enzyme concentration.
%[text]    IC50 depends on assay conditions.
%[text]    For competitive inhibitors: \\$\\text{IC}\_{50} = K_i (1 + [\\text{S}]/K_m)\\$. \\$K_i \\le \\text{IC}\_{50}\\$.
%%
%[text] ## Let's Try 4: Erlotinib's pIC50, Distribution Shape, and Log Scale Rationale

pIC50_erl_2nM = 9 - log10(2);
logInfo("Erlotinib pIC50 (IC50=2 nM): %.2f", pIC50_erl_2nM);

skewness_val = skewness(actTbl.pIC50);
logInfo("Skewness of pIC50 distribution: %.3f", skewness_val);

%[text] A: Erlotinib ~2 nM pIC50 = 9 - log10(2) = 8.70.
%[text]    pIC50 >= 8 is evaluated as high activity (less than 10 nM).
%[text]
%[text]    The IC50 distribution in ChEMBL shows positive skewness in pIC50 space (extremely potent outliers
%[text]    pull the distribution to the right). In raw IC50 nM space, the same data shows right skewness (many weakly active compounds
%[text]    and few strongly active compounds), a classic shape.
%[text]    The exact direction depends on the 50 records returned by ChEMBL.
%[text]    Reasons medicinal chemists use pIC50:
%[text]      (a) Dynamic range spans 6 orders of magnitude -- impractical on a linear axis.
%[text]      (b) Structure-activity relationships (SAR) are additive in log space:
%[text]          A single methyl group typically adds ~0.3 pIC50 units.
%[text]      (c) Binding free energy (\\$\\Delta G = RT \\ln K_i\\$) is proportional to \\$\\text{pK}_i\\$, so
%[text]          pIC50 is closer to a thermodynamic quantity.
%%
%[text] ## Let's Try 5: Number of Active Compounds, Most Potent Activity Name, Cutoff Relaxation

logInfo("Active compounds (IC50 <= %d nM): %d", IC50_CUTOFF_NM, height(activeTbl));
logInfo("Most potent activity: %s  IC50=%.2f nM  pIC50=%.2f", ...
    activeTbl.MoleculeChEMBLID(1), activeTbl.Value_nM(1), activeTbl.pIC50(1));

%[text] Relax the cutoff to 1000 nM
activeTbl_1uM = actTbl(actTbl.Value_nM <= 1000, :);
logInfo("Active compounds (IC50 <= 1000 nM): %d", height(activeTbl_1uM));

%[text] A: Relaxing to 1000 nM adds a layer of "moderate" effect size.
%[text]    In actual drug discovery, 1 uM compounds are often called "hits" in primary screening,
%[text]    and 100 nM compounds are often called "leads" in optimization.
%[text]    The appropriate cutoff depends on the target and assay.
%%
%[text] ## Let's Try 6: SDF File Structure, Reloading, and Record Count
%[text] Re-export (if runDir exists); if not, only demonstrate the concept of reloading
if numel(validMols) > 0
    runDir  = makeRunDir("Prefix", "s07_egfr_ans");
    sdfPath = fullfile(runDir, "egfr_actives_ans.sdf");
    emk.io.writeSdf(validMols, sdfPath);
    logInfo("SDF written: %s", sdfPath);

    % $$$$ delimiter count = record count
    txt       = fileread(sdfPath);
    lines     = strsplit(txt, newline);
    nRecords  = sum(strcmp(lines, "$$$$"));
    logInfo("SDF record count (measured by '$$$$'): %d", nRecords);

    % Confirm it matches numel(validMols)
    logInfo("Expected record count: %d  Match: %d", numel(validMols), nRecords == numel(validMols));
end

%[text] A: An SDF record consists of 3 parts:
%[text] 1. MOL block: Header (3 lines) + Bond table (atom/bond lines) + "M  END"
%[text] 2. Optional SD data fields: Values follow "> <Field Name>"
%[text] 3. Record delimiter: Single line "$$$$"
%[text]
%[text]   RDKit automatically generates 2D coordinates from SMILES when writing SDF,
%[text]   allowing immediate inspection in ChemDraw or Maestro.
%[text]
%[text]   Counting "$$$$ " lines provides the exact number of structure records.
%%
%[text] ## Let's Try 7: Top Tanimoto Compound Matching and Gefitinib Comparison

T_top = rankResult.Scores(1);
top_chembl = metaTbl.MoleculeChEMBLID(rankResult.Indices(1));
top_ic50   = metaTbl.Value_nM(rankResult.Indices(1));
logInfo("Top compound against Erlotinib: %s  T=%.4f  IC50=%.1f nM", ...
    top_chembl, T_top, top_ic50);

%[text] Gefitinib Comparison
GEFITINIB_SMILES = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1";
mol_gef = emk.mol.fromSmiles(GEFITINIB_SMILES);
fp_gef  = emk.fingerprint.morgan(mol_gef, Radius=2, NBits=2048);
res_gef = emk.similarity.rankBy(fp_gef, dbFps, Inf, Metric="tanimoto");

logInfo("Top 3 Gefitinib (ChEMBLID / T / IC50):");
for k = 1:min(3, numel(res_gef.Indices))
    idx = res_gef.Indices(k);
    logInfo("  %d. %s  T=%.4f  IC50=%.1f nM", k, ...
        metaTbl.MoleculeChEMBLID(idx), res_gef.Scores(k), metaTbl.Value_nM(idx));
end

%[text] A: If the Tanimoto score of the top compound and Erlotinib is > 0.7, they are structurally very similar.
%[text]    Erlotinib and Gefitinib are both EGFR inhibitors with a quinazoline core, but
%[text]    due to different substituents, ECFP4 screening does not necessarily identify the same top compounds.
%[text]    Compounds matching both queries are most likely to have a quinazoline core.
%[text]    Different top compounds reflect the diversity of the EGFR inhibitor chemical space.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
