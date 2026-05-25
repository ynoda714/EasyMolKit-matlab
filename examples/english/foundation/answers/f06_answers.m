%[text] # F06 Answers: File I/O
%[text] Reference answers for the exercises in `f06_file_io.m` (E1 to E3).
addpath(genpath("src"));
emk.setup.initPython();
runDir = makeRunDir();
%%
%[text] ## Answer E1: SMILES → SDF Round-Trip Validation

my_smiles = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", ...
             "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"];
my_names  = ["ethanol", "benzene", "aspirin", "caffeine", "ibuprofen"];

%[text] Writing out the SMILES list.
smi_path = fullfile(runDir, "e1_molecules.smi");
lines = strings(numel(my_smiles), 1);
for i = 1:numel(my_smiles)
    lines(i) = sprintf("%s    %s", my_smiles(i), my_names(i));
end
writelines(lines, smi_path);

%[text] Read, write SDF, read back, and compare contents.
mols_in  = emk.io.readSmilesList(smi_path);
sdf_path = fullfile(runDir, "e1_output.sdf");
emk.io.writeSdf(mols_in, sdf_path);
mols_out = emk.io.readSdf(sdf_path);

logInfo("--- Round-Trip Valid SMILES ---");
for i = 1:numel(mols_in)
    smi_before = emk.mol.toSmiles(mols_in{i});
    smi_after  = emk.mol.toSmiles(mols_out{i});
    match = isequal(smi_before, smi_after);
    logInfo("  %-14s  Match=%d  [%s]", my_names(i), match, smi_before);
end
%[text] The valid SMILES should all be the same before and after the SDF round-trip.
%%
%[text] ## Answer E2: Ro5 Filter and Save for Everyday Chemicals

data = readtable("data/list/everyday_chemicals.csv", TextType="string");
valid_mask = cellfun(@(s) emk.mol.isValid(s), cellstr(data.SMILES));
mols_valid = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(data.SMILES(valid_mask)), ...
    "UniformOutput", false);

desc_tbl    = emk.descriptor.batchCalculate(mols_valid, ...
    ["MolWt","LogP","NumHDonors","NumHAcceptors"]);
lipinski_tbl = emk.filter.lipinski(desc_tbl);
pass_mask    = lipinski_tbl.Pass_Ro5;

logInfo("Everyday chemicals passing Ro5: %d / %d", sum(pass_mask), numel(mols_valid));

sdf_out = fullfile(runDir, "everyday_ro5_pass.sdf");
emk.io.writeSdf(mols_valid(pass_mask), sdf_out);
logInfo("Save location: %s", sdf_out);
%[text] Most everyday chemicals are small molecules (solvents, fragrances), 
%[text] so it is expected that most will pass Ro5 (approximately 25-30 out of 30).
%%
%[text] ## Answer E3: Full Mini-Pipeline — Aspirin-Like FDA Drugs
%[text] **Step a**: Load FDA drugs and apply Ro5 filter.
fda = readtable("data/list/fda_drugs.csv", TextType="string");
valid_fda = cellfun(@(s) emk.mol.isValid(s), cellstr(fda.SMILES));
fda_mols  = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(fda.SMILES(valid_fda)), ...
    "UniformOutput", false);
fda_valid = fda(valid_fda, :);

desc_fda     = emk.descriptor.batchCalculate(fda_mols, ...
    ["MolWt","LogP","NumHDonors","NumHAcceptors"]);
lip_fda      = emk.filter.lipinski(desc_fda);
pass_fda     = lip_fda.Pass_Ro5;

passing_mols = fda_mols(pass_fda);
logInfo("FDA molecules passing Ro5: %d", sum(pass_fda));

%[text] **Step b**: Calculate Morgan fingerprints for all passing molecules.
fps_fda = cellfun(@(m) emk.fingerprint.morgan(m), passing_mols, ...
    "UniformOutput", false);

%[text] **Step c**: Rank by Tanimoto similarity to Aspirin.
aspirin_fp = emk.fingerprint.morgan(emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));
result = emk.similarity.rankBy(aspirin_fp, fps_fda, 10);

passing_names = fda_valid.Name(pass_fda);
logInfo("Top 10 FDA drugs most similar to Aspirin:");
for k = 1:10
    idx = result.Indices(k);
    logInfo("  %d. %-25s  Tanimoto=%.4f", k, passing_names(idx), result.Scores(k));
end

%[text] **Step d**: Save the top 10 to SDF.
top_mols = passing_mols(result.Indices(1:10));
sdf_top  = fullfile(runDir, "fda_aspirin_similar_top10.sdf");
emk.io.writeSdf(top_mols, sdf_top);
logInfo("Saved top 10: %s", sdf_top);
%%
%[text] ## Try It 1: Behavior when passing broken SMILES to readSmilesList
%[text] `emk.io.readSmilesList` skips invalid SMILES without crashing and logs 
%[text] the line number and SMILES string to `logWarn`.
%[text] You can create a .smi file containing invalid entries as follows to check the behavior.

invalid_smi_path = fullfile(runDir, "try1_invalid_test.smi");
writelines([ ...
    "CCO  ethanol", ...
    "INVALID_SMILES  bad_entry", ...
    "c1ccccc1  benzene" ...
], invalid_smi_path);

mols_y1 = emk.io.readSmilesList(invalid_smi_path);
logInfo("Try It 1 -- Number of molecules read: %d (1 invalid skipped)", numel(mols_y1));
%[text] **Observation Points**
%[text] - The log output will show `[WARN]` `readSmilesList: line 2, SMILES 'INVALID_SMILES' skipped`
%[text] - The function does not crash and `numel(mols_y1) == 2`
%[text] - Finally, a summary `readSmilesList: loaded 2 molecules ... (1 skipped)` will be output.
%%
%[text] ## Try It 2: How does the number of passes change with relaxed rules?

fda_y2 = readtable("data/list/fda_drugs.csv", TextType="string");
valid_mask_y2 = cellfun(@(s) emk.mol.isValid(s), cellstr(fda_y2.SMILES));
mols_y2 = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(fda_y2.SMILES(valid_mask_y2)), ...
    "UniformOutput", false);
fda_y2_valid = fda_y2(valid_mask_y2, :);

desc_tbl_y2 = emk.descriptor.batchCalculate(mols_y2, ...
    ["MolWt", "LogP", "NumHDonors", "NumHAcceptors"]);
desc_tbl_y2.Name = fda_y2_valid.Name;

lipinski_tbl         = emk.filter.lipinski(desc_tbl_y2);
lipinski_tbl_relaxed = emk.filter.lipinski(desc_tbl_y2, "MaxViolations", 1);

pass_strict  = lipinski_tbl.Pass_Ro5;
pass_relaxed = lipinski_tbl_relaxed.Pass_Ro5;

logInfo("Strict Ro5 passed: %d / %d", sum(pass_strict), numel(pass_strict));
logInfo("Relaxed Ro5 passed (MaxViolations=1): %d / %d", sum(pass_relaxed), numel(pass_relaxed));
logInfo("Newly passed drugs (%d):", sum(pass_relaxed & ~pass_strict));

newly_pass = desc_tbl_y2.Name(pass_relaxed & ~pass_strict);
for k = 1:numel(newly_pass)
    logInfo("  %s", newly_pass(k));
end
%[text] The relaxed filter with MaxViolations=1 typically adds a few to several additional passes.
%[text] These are compounds that slightly exceed one of the Ro5 rules and are "almost orally bioavailable." 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---