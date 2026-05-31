%[text] # F06: Reading and Writing Molecular Files
%[text] EasyMolKit Basic Tutorial — Layer 1
%[text] 
%[text] In actual research projects and drug discovery, you will handle thousands to tens of thousands of compounds at once. It is impossible to input these one by one manually. Therefore, the technology of "file input/output (I/O)" becomes important for directly exchanging data with major chemical databases around the world (such as PubChem and ChEMBL) and external docking simulation software. In this tutorial, you will master reading and writing two major formats: the simplest and lightweight "SMILES list file" and the industry-standard "SDF (Structure-Data File)" that can also pack the 3D shapes of molecules and various internal data together. Furthermore, we will challenge ourselves to build a "mini data pipeline" that automatically processes data and saves results to files by linking all the knowledge learned in F01 to F05 (structure checking, property calculation, substructure filtering, etc.).
%[text] **Learning Objectives**
%[text] - Read SMILES list files with `emk.io.readSmilesList`
%[text] - Read SDF files with `emk.io.readSdf`
%[text] - Write molecules to SDF files with `emk.io.writeSdf`
%[text] - Combine file I/O with descriptor calculation and filtering \
%[text] **Prerequisites**
%[text] - RDKit installed (run `emk.setup.install()` once)
%[text] - F01 to F04
%[text] - No additional Toolbox required (works with MATLAB only) \
%[text] Estimated time: 10–15 minutes | How to run: Execute sections one by one with Ctrl+Enter
%%
%[text] ## Section 0: Setup
%[text] Initialize the path and Python environment.
%[text] **Always run this section first.**
% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
emk.setup.initPython(); %[output:96fa9617]
%[text] Warm up Python/RDKit process
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;

runDir = makeRunDir();   % Create result/runs/<timestamp>/
logInfo("F06: Setup complete -- Output directory: %s", runDir); %[output:8dd385bb]
%%
%[text] ## Section 1: Read SMILES List File
%[text] ### **Simple molecular name list with one text line each**
%[text] The simplest way to save molecular data is to create a text file where "SMILES strings" and "molecule names (IDs)" are separated by tabs or spaces, one per line (commonly with `.smi` or `.txt` extensions). This format is very small in file size and can be directly opened in a text editor for human verification, making it very useful when you want to roughly share a large list of molecules.
%[text] In EasyMolKit, you can read all molecules from this file at once by simply calling the `emk.io.readSmilesList` function. The loaded data is automatically converted into a "cell array of molecular objects," which is the most user-friendly format for MATLAB users. Now, let's read a pre-prepared list of common chemicals (`everyday_chemicals.smi`) and check how the molecules are neatly arranged in the MATLAB workspace.
%[text] 
%[text] `emk.io.readSmilesList` reads a text file with "1 line 1 SMILES." If there is a second column (tab or space-separated), it is treated as a name/label.
%[text] Comment lines starting with `#` and empty lines are automatically skipped.
%[text] 
%[text] The return value is a cell array of RDKit Mol objects. Even if there are invalid SMILES, it will skip them with a warning, so it won't crash.
%[text] The accompanying `everyday_chemicals.csv` is comma-separated, so in this section, we will write out a plain SMILES list for practice.
smiles_file = fullfile(runDir, "sample_molecules.smi");
writelines([ ...
    "# F06 Tutorial Sample SMILES List", ...
    "CCO            ethanol", ...
    "c1ccccc1       benzene", ...
    "CC(=O)Oc1ccccc1C(=O)O  aspirin", ...
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C  caffeine", ...
    "CC(=O)NC1=CC=C(C=C1)O  acetaminophen" ...
], smiles_file);

mols_from_smi = emk.io.readSmilesList(smiles_file); %[output:390d8d09]
logInfo("Read %d molecules from the SMILES list", numel(mols_from_smi)); %[output:2c91ec2c]
%[text] 
%[text] **✏️ Try It 1 — What happens if you mix in broken SMILES?**
%[text] Add an invalid SMILES line (e.g., `"INVALID_SMILES  bad_entry"`) to `sample_molecules.smi` and try running this section again.
%[text] Will the function crash? Or will it skip with a warning?
%[text] Hint: Check if a `logWarn` message appears in the log output.
% ... (Let's write code here)
%%
%[text] ## Section 2: Write Molecules to SDF
%[text] ### **A format that bundles both structure and data**
%[text] While SMILES lists are convenient and easy to use, they cannot store complex information such as "3D coordinates of each atom" or "activity values obtained from experiments (property data)" together. Therefore, the format used as the de facto industry standard in the pharmaceutical and chemical fields around the world is "SDF (Structure-Data File)." SDF can store not only the bonding relationships (connectivity table) between atoms that make up the molecule but also spatial coordinate information and various text data (such as molecular weight and notes) added by the user, all packed into one file.
%[text] 
%[text] ### Concept: Why SDF continues to be used as an industry standard
%[text] SDF (Structure-Data File) is a molecular exchange format introduced by MDL Information Systems in the 1980s. There are three reasons why it has been used as an industry standard for many years.
%[text] 1. **Coordinate retention**: It can store 2D/3D atomic coordinates (not possible with SMILES alone). Essential for docking and 3D pharmacophore analysis.
%[text] 2. **Property inclusion**: Active data and IDs can be stored in one file together with the structure using the `> <PROP_NAME>` field.
%[text] 3. **Safe concatenation**: Multiple molecules can be safely concatenated using delimiters. \
%[text] PubChem, ChEMBL, and Reaxys all provide downloads in SDF format, and docking engines such as AutoDock and Glide accept SDF as input. Reading and writing SDF is the most practical file I/O skill in cheminformatics. `emk.io.writeSdf` writes a cell array of Mol objects to disk.
%[text] The output path can safely use the directory created by `makeRunDir()`.
sdf_out = fullfile(runDir, "sample_output.sdf");
emk.io.writeSdf(mols_from_smi, sdf_out); %[output:72cef0d6]
logInfo("Wrote %d molecules: %s", numel(mols_from_smi), sdf_out); %[output:98ae2446]
%%
%[text] ## Section 3: Read SDF Again
%[text] ### **Verifying the Round-Trip: Read the SDF Back**
%[text] Now let's read back the SDF file we just wrote in Section 2. This "round-trip" test — write then read — confirms that no structural information was lost during the conversion. `emk.io.readSdf` reads an SDF file and returns a cell array of Mol objects, just like `readSmilesList`. The canonical SMILES displayed below may differ from the original input strings, which is expected: RDKit always normalizes the same molecule to the same canonical form.
mols_from_sdf = emk.io.readSdf(sdf_out); %[output:998b8a2d]
logInfo("Reloaded %d molecules from SDF", numel(mols_from_sdf)); %[output:91439fb7]
%[text] 
%[text] Check if the canonical SMILES are retained after the round trip. Note that the displayed SMILES may differ from the input (e.g., Caffeine).
%[text] This is normal behavior as RDKit always converts the same molecule to the same "canonical SMILES."
logInfo("--- Round Trip SMILES Check ---"); %[output:82318383]
for i = 1:numel(mols_from_sdf) %[output:group:2778ade1]
    smi = emk.mol.toSmiles(mols_from_sdf{i});
    logInfo("  Mol %d: %s", i, smi); %[output:66806ecf]
end %[output:group:2778ade1]
%%
%[text] ## Section 4: Use Attached Datasets
%[text] ### **Building an Automated Molecular Processing Pipeline**
%[text] Let's combine all the techniques we've learned so far to build an automated processing system (data pipeline) that professional cheminformaticians use daily. Our goal is to execute the entire flow of **"① Load files ➔ ② Automatically eliminate invalid structures ➔ ③ Calculate property values for all ➔ ④ Filter out non-drug-like molecules ➔ ⑤ Save only the passers to a new file"** automatically with the push of a button.
%[text] Once you set up such a pipeline, even if the number of files increases from 10 to 10,000, the computer will quietly and accurately process them in the background. Let's experience a workflow that combines MATLAB's strengths in matrix and array processing with RDKit's chemical engine, which is the true essence of EasyMolKit.
%[text] EasyMolKit comes with curated lists of molecules under `data/list/`.
%[text] -  `everyday_chemicals.csv`  — 30 types of daily chemicals (PubChem CC0)
%[text] -  `fda_drugs.csv`           — 200 types of FDA-approved drugs (ChEMBL CC-BY-SA 3.0)
%[text] -  `pains.csv`               — PAINS SMARTS filter (BSD-3)
%[text] -  `forensic_challenge.csv`  — Virtual forensic challenge set (labeled) \
%[text] Here, we will read FDA-approved drugs from `fda_drugs.csv` and convert them to Mol objects.
fda = readtable("data/list/fda_drugs.csv", TextType="string");
logInfo("FDA Drug Dataset: %d compounds", height(fda)); %[output:4a25a224]
%[text] 
%[text] Parse the SMILES (there are 200 molecules, so it may take a few seconds). First, check validity with `isValid`, and convert only valid molecules with `fromSmiles`.
%[text] By doing this in order, you can safely skip invalid SMILES.
valid_mask = false(height(fda), 1);
fda_mols   = cell(height(fda), 1);
for i = 1:height(fda)
    if emk.mol.isValid(fda.SMILES(i))
        fda_mols{i}   = emk.mol.fromSmiles(fda.SMILES(i));
        valid_mask(i) = true;
    end
end
logInfo("Valid FDA molecules: %d / %d", sum(valid_mask), height(fda)); %[output:8597961b]
%%
%[text] ## Section 5: Filter and Export Subsets
%[text] Filter drug-like molecules using Lipinski's Ro5 (Rule of Five) and save to SDF.
fda_valid_mols = fda_mols(valid_mask);
fda_valid      = fda(valid_mask, :);

%[text] Calculate descriptor table.
desc_tbl = emk.descriptor.batchCalculate(fda_valid_mols, ...
    ["MolWt", "LogP", "NumHDonors", "NumHAcceptors"]);
%[text] 
%[text] Apply Lipinski filter.
lipinski_tbl = emk.filter.lipinski(desc_tbl);   % Add Pass_Ro5 / Violations_Ro5 %[output:603355c3]

pass_mask = lipinski_tbl.Pass_Ro5;
logInfo("Ro5 Passed: %d / %d FDA molecules", sum(pass_mask), numel(fda_valid_mols)); %[output:28b668da]
%[text] 
%[text] **✏️ Try It 2 — How does the number of passes change with relaxed rules?**
%[text] Try the relaxed filter with the `MaxViolations=1` option to allow one item to exceed Ro5.
%[text] `lipinski_tbl_relaxed = emk.filter.lipinski(desc_tbl, "MaxViolations", 1);`
%[text] How many more pass compared to strict Ro5?
%[text] Let's also check the names of the newly passed drugs.
%[text] (Note: MaxViolations=1 is a relaxed Ro5 that allows "up to 1 exceedance of 4 criteria."
%[text]   This is a different concept from the rules of Veber et al. 2002 (TPSA ≤ 140 Å², number of rotatable bonds ≤ 10).)
%[text] 
%[text] Save the molecules that passed Ro5 to SDF.
passing_mols = fda_valid_mols(pass_mask);
sdf_filtered = fullfile(runDir, "fda_ro5_pass.sdf");
emk.io.writeSdf(passing_mols, sdf_filtered); %[output:4ec9371a]
logInfo("Saved Ro5 passing molecules: %s", sdf_filtered); %[output:3e55f303]
%[text] 
%[text] Also save the descriptor table to CSV for later review.
csv_out = fullfile(runDir, "fda_ro5_pass_descriptors.csv");
writetable([fda_valid(pass_mask, "Name"), desc_tbl(pass_mask, :), ...
    lipinski_tbl(pass_mask, ["Violations_Ro5","Pass_Ro5"])], csv_out);
logInfo("Saved descriptor table: %s", csv_out); %[output:93bcbda8]
%%
%[text] ## Section 6: Summary of Mini Pipeline
%[text] This is a list of functions from the file I/O module learned in this tutorial and a template for a standard pipeline design.
%[text] -  Step 1 -- Read molecules from a file  `emk.io.readSmilesList` / readSdf
%[text] -  Step 2 -- Validity check               `emk.mol.isValid`
%[text] -  Step 3 -- Descriptor calculation       `emk.descriptor.batchCalculate`
%[text] -  Step 4 -- Apply filter                `emk.filter.lipinski`
%[text] -  Step 5 -- Fingerprint + similarity    `emk.fingerprint.morgan`, `emk.similarity.*` (Learned in F04 and F05. This script omits it)
%[text] -  Step 6 -- Export results              `emk.io.writeSdf` / writetable \
%[text] ## Section 7: Summary
%[text] 
%[text:table]
%[text] | Function | Purpose |
%[text] | --- | --- |
%[text] | `emk.io.readSmilesList(path)` | SMILES list → Cell array of Mol |
%[text] | `emk.io.readSdf(path)` | SDF file → Cell array of Mol |
%[text] | `emk.io.writeSdf(mols, path)` | Cell array of Mol → SDF file |
%[text] | `makeRunDir()` | Create output directory with timestamp |
%[text:table]
%[text] 
%[text] **Typical Pipeline**
%[text]     mols -\> isValid -\> batchCalculate -\> filter.lipinski -\> writeSdf
%[text] 
%[text] **Congratulations!** You have completed the Foundation tutorial (F01 to F06).
%[text] **Next Step**: Let's move on to the application story tutorials in `examples/stories/` (Layer 2).
%%
%[text] ## Exercises
%[text] Try to solve each exercise before referring to `answers/f06_answers.m`.
%[text] 
%[text] **E1.** Write out a SMILES list file for 5 molecules, read it back with `readSmilesList`, and
%[text]     validate the round trip. Display the normalized SMILES before and after the SDF round trip
%[text]     and confirm they match.
runDir = makeRunDir();
% Write out SMILES file, read with readSmilesList, write SDF, read back SDF
% Compare toSmiles output before and after round trip
%[text] 
%[text] **E2.** Load `everyday_chemicals.csv`, calculate descriptors, and save the results after applying the Lipinski filter to the SDF file `"everyday_ro5_pass.sdf"`.
%[text]     How many of the 30 types of everyday chemicals pass Ro5?
data = readtable("data/list/everyday_chemicals.csv", TextType="string");
% Analyze -> Descriptors -> lipinski -> writeSdf
%[text] 
%[text] **E3.** Mini pipeline challenge: Starting from `fda_drugs.csv`
%[text]     (a) Filter for Ro5 passing molecules,
%[text]     (b) Calculate Morgan fingerprints for all passing molecules,
%[text]     (c) Rank by similarity to Aspirin (`"CC(=O)Oc1ccccc1C(=O)O"`),
%[text]     (d) Save the top 10 FDA drugs most similar to Aspirin to an SDF file.
%[text]     This is a culmination of skills learned in F02, F03, F04, and F06.
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:96fa9617]
%   data: {"dataType":"text","outputData":{"text":"[09:40:24][INFO]  initPython: Python 3.10 already active (Status: Loaded) -- skipping.\n","truncated":false}}
%---
%[output:8dd385bb]
%   data: {"dataType":"text","outputData":{"text":"[09:40:24][INFO]  F06: Setup complete -- Output directory: result\\runs\\20260525_094024\n","truncated":false}}
%---
%[output:390d8d09]
%   data: {"dataType":"text","outputData":{"text":"[09:40:24][INFO]  readSmilesList: loaded 5 molecules from 'result\\runs\\20260525_094024\\sample_molecules.smi' (0 skipped)\n","truncated":false}}
%---
%[output:2c91ec2c]
%   data: {"dataType":"text","outputData":{"text":"[09:40:24][INFO]  Read 5 molecules from the SMILES list\n","truncated":false}}
%---
%[output:72cef0d6]
%   data: {"dataType":"text","outputData":{"text":"[09:40:24][INFO]  writeSdf: wrote 5 molecules to 'result\\runs\\20260525_094024\\sample_output.sdf'\n","truncated":false}}
%---
%[output:98ae2446]
%   data: {"dataType":"text","outputData":{"text":"[09:40:25][INFO]  Wrote 5 molecules: result\\runs\\20260525_094024\\sample_output.sdf\n","truncated":false}}
%---
%[output:998b8a2d]
%   data: {"dataType":"text","outputData":{"text":"[09:40:25][INFO]  readSdf: loaded 5 molecules from 'result\\runs\\20260525_094024\\sample_output.sdf' (0 skipped)\n","truncated":false}}
%---
%[output:91439fb7]
%   data: {"dataType":"text","outputData":{"text":"[09:40:25][INFO]  Reloaded 5 molecules from SDF\n","truncated":false}}
%---
%[output:82318383]
%   data: {"dataType":"text","outputData":{"text":"[09:40:25][INFO]  --- Round Trip SMILES Check ---\n","truncated":false}}
%---
%[output:66806ecf]
%   data: {"dataType":"text","outputData":{"text":"[09:40:25][INFO]    Mol 1: CCO\n[09:40:25][INFO]    Mol 2: c1ccccc1\n[09:40:25][INFO]    Mol 3: CC(=O)Oc1ccccc1C(=O)O\n[09:40:25][INFO]    Mol 4: Cn1c(=O)c2c(ncn2C)n(C)c1=O\n[09:40:25][INFO]    Mol 5: CC(=O)Nc1ccc(O)cc1\n","truncated":false}}
%---
%[output:4a25a224]
%   data: {"dataType":"text","outputData":{"text":"[09:40:25][INFO]  FDA Drug Dataset: 200 compounds\n","truncated":false}}
%---
%[output:8597961b]
%   data: {"dataType":"text","outputData":{"text":"[09:40:39][INFO]  Valid FDA molecules: 200 \/ 200\n","truncated":false}}
%---
%[output:603355c3]
%   data: {"dataType":"text","outputData":{"text":"[09:40:42][INFO]  lipinski: 176 \/ 200 row(s) pass Ro5 (MaxViolations=0)\n","truncated":false}}
%---
%[output:28b668da]
%   data: {"dataType":"text","outputData":{"text":"[09:40:42][INFO]  Ro5 Passed: 176 \/ 200 FDA molecules\n","truncated":false}}
%---
%[output:4ec9371a]
%   data: {"dataType":"text","outputData":{"text":"[09:40:47][INFO]  writeSdf: wrote 176 molecules to 'result\\runs\\20260525_094024\\fda_ro5_pass.sdf'\n","truncated":false}}
%---
%[output:3e55f303]
%   data: {"dataType":"text","outputData":{"text":"[09:40:47][INFO]  Saved Ro5 passing molecules: result\\runs\\20260525_094024\\fda_ro5_pass.sdf\n","truncated":false}}
%---
%[output:93bcbda8]
%   data: {"dataType":"text","outputData":{"text":"[09:40:47][INFO]  Saved descriptor table: result\\runs\\20260525_094024\\fda_ro5_pass_descriptors.csv\n","truncated":false}}
%---
