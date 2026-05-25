%[text] # F05 Answers: Substructure Search
%[text] Reference answers for f05_substructure_search.m exercises.
addpath(genpath("src"));
emk.setup.initPython();
%%
%[text] ## Answer E1: Finding Phenolic Compounds in Everyday Chemicals

data = readtable("data/list/everyday_chemicals.csv", TextType="string");

logInfo("Phenolic compounds in everyday chemicals (aromatic ring + hydroxyl):");
count = 0;
for i = 1:height(data)
    if ~emk.mol.isValid(data.SMILES(i)); continue; end
    mol = emk.mol.fromSmiles(data.SMILES(i));
    has_ring = emk.mol.hasSubstruct(mol, "c1ccccc1");
    has_oh   = emk.mol.hasSubstruct(mol, "[OH]");
    if has_ring && has_oh
        count = count + 1;
        logInfo("  %s", data.CommonName(i));
    end
end
logInfo("Total phenolic compounds found: %d", count);
%[text] Phenolic compounds are common in disinfectants (thymol, eugenol), fragrances,
%[text] and analgesics (paracetamol/acetaminophen).
%%
%[text] ## Answer E2: Number of Carboxylic Acids in FDA Drugs

fda = readtable("data/list/fda_drugs.csv", TextType="string");
cooh_count = 0;
logInfo("Scanning for carboxylic acids in %d FDA drugs...", height(fda));
for i = 1:height(fda)
    if ~emk.mol.isValid(fda.SMILES(i)); continue; end
    mol = emk.mol.fromSmiles(fda.SMILES(i));
    if emk.mol.hasSubstruct(mol, "C(=O)[OH]")
        cooh_count = cooh_count + 1;
    end
    logProgress(i, height(fda), "FDA drugs");
end
frac = cooh_count / height(fda) * 100;
logInfo("Carboxylic acid content in FDA drugs: %d / %d  (%.1f%%)", ...
    cooh_count, height(fda), frac);
%[text] Many NSAIDs (Aspirin, Ibuprofen, Naproxen, etc.) and several ACE inhibitors contain
%[text] carboxylic acids. About 20-30% is expected overall.
%%
%[text] ## Answer E3: Functional Group Summary for Any Molecule

smiles_test = {"CN1CCC[C@H]1C2=CN=CC=C2", ...  % Nicotine
               "OC1=CC2=C(CC3N(CC23)CC4=CC=CC=C14)C=C1"};  % Morphine
mol_names   = {"Nicotine", "Morphine"};

fg_patterns = {"c1ccccc1", "[OH]", "C(=O)[OH]", "[NH,NH2]", "C(=O)[N]", ...
               "[F,Cl,Br,I]", "C(=O)[#6]"};
fg_names    = {"Aromatic ring","Hydroxyl","Carboxylic acid","Amine","Amide", ...
               "Halogen","Ketone"};

for m = 1:numel(smiles_test)
    if ~emk.mol.isValid(smiles_test{m}); continue; end
    mol = emk.mol.fromSmiles(smiles_test{m});
    logInfo("--- %s ---", mol_names{m});
    for p = 1:numel(fg_patterns)
        has = emk.mol.hasSubstruct(mol, fg_patterns{p});
        logInfo("  %-14s : %s", fg_names{p}, string(has));
    end
end
%%
%[text] ## Try It 1: Let's Test SMARTS with Caffeine

mol_caf = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");

has_n   = emk.mol.hasSubstruct(mol_caf, "[#7]");
has_ket = emk.mol.hasSubstruct(mol_caf, "[C](=O)[#6]");
has_oh  = emk.mol.hasSubstruct(mol_caf, "[OH]");
has_ami = emk.mol.hasSubstruct(mol_caf, "[C](=O)[N]");

logInfo("Caffeine -- Any nitrogen [#7]: %d", has_n);
logInfo("Caffeine -- Ketone [C](=O)[#6]: %d", has_ket);
logInfo("Caffeine -- Hydroxyl [OH]: %d", has_oh);
logInfo("Caffeine -- Amide [C](=O)[N]: %d", has_ami);
%[text] (a) [#7] = true: Caffeine has 4 nitrogen atoms.
%[text] (b) [C](=O)[#6] = false: The C=O group is an amide/imide carbonyl, adjacent to nitrogen (N), not carbon.
%[text] (c) [OH] = false: Caffeine does not have a hydroxyl group.
%[text] (d) [C](=O)[N] = true: Amide carbonyls (2) are detected.
%%
%[text] ## Try It 2: Let's Add a Halogen Column

data_hal = readtable("data/list/everyday_chemicals.csv", TextType="string");
fg_pat2 = {"c1ccccc1", "[OH]", "C(=O)[OH]", "[NH,NH2]", "[C](=O)[N]", "[F,Cl,Br,I]"};
fg_nm2  = ["Benzene ring", "Hydroxyl", "Carboxylic acid", "Amine", "Amide", "Halogen"];

n_hal = height(data_hal);
fg_mat2 = false(n_hal, numel(fg_pat2));
for i = 1:n_hal
    if ~emk.mol.isValid(data_hal.SMILES(i)); continue; end
    mol = emk.mol.fromSmiles(data_hal.SMILES(i));
    for p = 1:numel(fg_pat2)
        fg_mat2(i, p) = emk.mol.hasSubstruct(mol, fg_pat2{p});
    end
end
fg_tbl2 = array2table(fg_mat2, "VariableNames", fg_nm2);
fg_tbl2.Name = data_hal.CommonName;
fg_tbl2 = movevars(fg_tbl2, "Name", "Before", "Benzene ring");

logInfo("Halogen-containing everyday chemicals:");
halogen_rows = fg_tbl2(fg_tbl2.Halogen, ["Name", "Halogen"]);
disp(halogen_rows);
%[text] Many everyday chemicals do not contain halogens. Chlorine-containing disinfectants
%[text] (such as triclocarban) may be found. 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---