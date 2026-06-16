%[text] # R08: Docking Simulation — SARS-CoV-2 Mpro vs 10 Antiviral Drugs
%[text] EasyMolKit Research — Layer 4 (<u>**MATLAB Online Only**</u>)
%[text] Why does Nirmatrelvir (the active ingredient in Paxlovid) strongly inhibit the SARS-CoV-2 main protease (Mpro), while earlier HIV drugs fell short in clinical trials?
%[text] Lopinavir and Ritonavir — touted as early COVID-19 candidates — proved ineffective in controlled trials. Meanwhile, the HCV drug Boceprevir turned out to bind Mpro surprisingly well, a result that has since been documented in the literature.
%[text] These contrasting outcomes come down to a single factor: the shape of the molecule.
%[text] Molecular docking lets us computationally "fit" a small-molecule ligand into a protein's active site and estimate the binding energy — a key step in rational drug design.
%[text] In this script, we dock 10 antiviral drugs against SARS-CoV-2 Mpro (PDB: 6LU7) and explore why some molecules score well and others do not, using numerical data and interactive 3D visualization.
%[text] ## Learning Objectives
%[text] - Understand the molecular docking pipeline (protein preparation → ligand 3D generation → PDBQT conversion → Vina execution).
%[text] - Understand the meaning of binding energy (kcal/mol) and its computational limitations.
%[text] - Consider why off-target drugs "accidentally" bind to Mpro from a structural perspective.
%[text] - Experience protein structure preprocessing with `pdbfixer`.
%[text] - Understand 3D structure generation with RDKit (EmbedMolecule + UFF optimization).
%[text] - Experience interactive 3D ribbon visualization with MATLAB's `uihtml` and 3Dmol.js. \
%[text] ## Overview of SARS-CoV-2 Mpro (3CLpro)
%[text] - Mpro (main protease) is a cysteine protease (EC 3.4.22.69) essential for the replication of SARS-CoV-2.
%[text] - The catalytic residues are His41 (general base) and Cys145 (nucleophile), forming a "catalytic dyad."
%[text] - Key residues involved in substrate recognition are Glu166 (recognizes Gln of P1), Met49, and Met165 (hydrophobic P2 pocket).
%[text] - The PDB structure 6LU7 used is a complex with the N3 covalent inhibitor (Jin et al. 2020). \
%[text] ## 10 Ligands
%[text] **Correct (Mpro Inhibitors):**
%[text] 1. Nirmatrelvir — Active ingredient of Paxlovid. Benchmark
%[text] 2. Ensitrelvir — Mpro inhibitor from Japan (Shionogi). Different design \
%[text] **Cross-Reactivity (Mpro ≠ Original Target, but Strong Binding):**
%[text] 1. Boceprevir — HCV NS3/4A protease inhibitor. Cross-reactivity with Mpro demonstrated in literature \
%[text] **Clinical Failures (Learn Discrepancy Between Structure Score and Trial Results):**
%[text] 1. Lopinavir — HIV protease inhibitor. Early COVID candidate but clinically ineffective
%[text] 2. Ritonavir — HIV PK booster. Component of Paxlovid but weak on its own
%[text] 3. Nelfinavir — HIV protease inhibitor. Limitations in Mpro affinity \
%[text] **Completely Different Targets (Predict Low Scores):**
%[text] 1. Remdesivir — RNA polymerase inhibitor. Non-target for protease
%[text] 2. Favipiravir — RNA polymerase inhibitor. Lightweight (MW=157) and does not fit in the active site
%[text] 3. Oseltamivir — Influenza neuraminidase inhibitor
%[text] 4. Acyclovir — Herpes thymidine kinase inhibitor (nucleic acid analog) \
%[text] ## Prerequisites
%[text] - Completion of A05 (Neural Networks) or R04 (Protein-Ligand) recommended
%[text] - No Toolbox required
%[text] - MATLAB Online required (meeko, vina, pdbfixer need a Linux backend) \
%[text] ## Operating Environment
%[text] This exercise is exclusive to **MATLAB Online** (Linux backend).
%[text] Reason for Desktop non-support: vina lacks a Windows pip wheel (requires Boost C++ build), and `pdbfixer` is blocked by Windows Smart App Control on `_openmm.pyd`.
%[text] In MATLAB Online, all three can be automatically resolved with `pip install`.
%[text] 
%[text] **Duration**: 10–30 minutes (depends on MATLAB Online CPU performance and molecule size. Large molecules like Lopinavir/Ritonavir take 1–2 minutes/ligand, small molecules take a few seconds)
%[text] 
%[text] ## Data Used
%[text] - Protein: PDB 6LU7 (download from RCSB in Section 1)
%[text] - Ligands: Obtain 10 normalized SMILES from PubChem REST API (Section 2) \
%[text] ## Simplification Notes
%[text] - Rigid receptor (no induced fit or flexible side chains)
%[text] - Gasteiger partial charges (not AM1-BCC)
%[text] - exhaustiveness = 4 (educational speed. Recommend 8 or more for publication quality)
%[text] - Nirmatrelvir: Non-covalent approximation (ignores covalent bond with Cys145) \
%[text] ## References
%[text] - Jin Z et al. (2020) Structure of Mpro from SARS-CoV-2 and discovery of its inhibitors. Nature 582:289-293. doi:10.1038/s41586-020-2223-y [6LU7 crystal structure and N3 inhibitor]
%[text] - Owen DR et al. (2021) An oral SARS-CoV-2 Mpro inhibitor clinical candidate. Science 374:1586-1593. doi:10.1126/science.abl4784 [Nirmatrelvir discovery paper]
%[text] - Ma C et al. (2020) Boceprevir, GC-376, and calpain inhibitors II, XII inhibit SARS-CoV-2 viral replication. Cell Res 30:678-692. doi:10.1038/s41422-020-0356-z [Boceprevir's Mpro cross-reactivity]
%[text] - Cao B et al. (2020) A trial of lopinavir-ritonavir in adults hospitalized with severe Covid-19. N Engl J Med 382:1787-1799. doi:10.1056/NEJMoa2001282 [Lopinavir-Ritonavir's COVID-19 clinical failure]
%[text] - Trott O & Olson AJ (2010) AutoDock Vina: improving the speed and accuracy of docking. J Comput Chem 31:455-461. doi:10.1002/jcc.21334 [Vina scoring function]
%[text] - Eberhardt J et al. (2021) AutoDock Vina 1.2.0. J Chem Inf Model 61:3891-3898. doi:10.1021/acs.jcim.1c00203 [Vina Python API] \
%[text] 
%[text] Please execute each section with Ctrl+Enter.
%%
%[text] ## Section 0: Setup — Online Guard and Library Installation
% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif isfolder(fullfile(pwd, 'src'))
    addpath(genpath(fullfile(pwd, 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
%[text] MATLAB Online Exclusive Guard
if ~emk.util.isOnline()
    logWarn("R08: This exercise is exclusive to MATLAB Online.");
    logWarn("     meeko / vina / pdbfixer do not work on Windows Desktop.");
    logWarn("     Reason: vina lacks Windows pip wheel (Boost C++ build required),");
    logWarn("           pdbfixer is blocked by Windows Smart App Control on _openmm.pyd.");
    logWarn("     Please run this script on MATLAB Online (https://matlab.mathworks.com).");
    emk.setup.recipe("docking");
    return
end
%[text] Enable docking library flag and pass to installOnline
cfg = emkLoadConfig();
cfg.useCase.docking = true;   % scipy + meeko + vina + pdbfixer
emk.setup.installOnline(Config=cfg);
%[text] Section 0a: Tuning Parameters (Modify here before execution)
EXHAUSTIVENESS = 4;       % Vina search iterations (4=educational speed, 8=publication quality, 32=high precision)
N_POSES        = 3;       % Number of docking poses to save
%[text] 6LU7 Mpro active site binding box (calculated from N3 inhibitor centroid; Chain A)
%[text] Reference: Jin et al. 2020 / AutoDock Vina Tutorial (widely cited coordinates)
BOX_CENTER = [-26.3, 12.7, 58.6];  % [x, y, z] Angstrom
BOX_SIZE   = [28.0, 28.0, 28.0];   % [dx, dy, dz] Angstrom
%[text] Output Directory
runDir = makeRunDir("Prefix", "r08_docking");
logInfo("R08: Output Directory -> %s", runDir);
%[text] Prepare Directory Structure
receptorDir = fullfile(runDir, "receptor");
ligandDir   = fullfile(runDir, "ligands");
poseDir     = fullfile(runDir, "poses");
mkdir(receptorDir);
mkdir(ligandDir);
mkdir(poseDir);
logInfo("R08: Section 0 Complete");
%%
%[text] ## Section 1: PDB Retrieval and Protein Preprocessing (pdbfixer)
PDB_ID        = "6LU7";
receptorPdb   = fullfile(receptorDir, "6lu7_raw.pdb");
preparedPdb   = fullfile(receptorDir, "6lu7_prepared.pdb");
receptorPdbqt = fullfile(receptorDir, "receptor.pdbqt");
%[text] \--- 1a: Download 6LU7 from RCSB ---
if ~isfile(receptorPdb)
    logInfo("R08: Downloading PDB %s from RCSB...", PDB_ID);
    pdbUrl = sprintf("https://files.rcsb.org/download/%s.pdb", upper(PDB_ID));
    try
        websave(receptorPdb, pdbUrl);
        logInfo("R08: PDB download complete -> %s", receptorPdb);
    catch ME
        error("emk:r08:pdbDownloadFailed", ...
            "Failed to download PDB %s: %s\n" + ...
            "Please check your internet connection.", PDB_ID, ME.message);
    end
else
    logInfo("R08: Using cached PDB -> %s", receptorPdb);
end
%[text] \--- 1b: Prepare protein with pdbfixer ---
%[text] 
%[text] ### Necessity of Protein Preprocessing
%[text] Crystal structures (PDB files) are not immediately suitable for docking — experimental constraints leave them with several well-known artifacts.
%[text] Common issues include missing residues (from regions of weak electron density), absent hydrogen atoms (X-rays detect electron density, not hydrogen positions), and co-crystallized ligands or solvent molecules.
%[text] `pdbfixer` automatically fills in missing residues and atoms, strips heteroatoms, and adds hydrogens, yielding a clean, protonated receptor structure ready for docking.
logInfo("R08: Preparing protein with pdbfixer...");
try
    pdbfixer = py.importlib.import_module("pdbfixer");

    % Create PDBFixer object
    fixer = pdbfixer.PDBFixer(filename=receptorPdb);
    logInfo("R08: PDB loaded");

    % Complete missing residues and heavy atoms
    fixer.findMissingResidues();
    fixer.findMissingAtoms();
    fixer.addMissingAtoms();
    logInfo("R08: Completion of missing atoms done");

    % Remove heteroatoms (ligands, water) to isolate the receptor
    % keepIds: Residues to keep (empty = remove all heteroatoms)
    fixer.removeHeterogens(false);
    logInfo("R08: Removal of heteroatoms (ligands, crystallographic water) done");

    % Protonate at pH 7.4
    fixer.addMissingHydrogens(py.float(7.4));
    logInfo("R08: Hydrogen addition complete (pH 7.4)");

    % Save prepared protein in PDB format
    openmm = py.importlib.import_module("openmm.app");
    fh     = py.open(preparedPdb, "w");
    openmm.PDBFile.writeFile(fixer.topology, fixer.positions, fh);
    fh.close();
    logInfo("R08: Prepared receptor saved -> %s", preparedPdb);

catch ME
    error("emk:r08:pdbfixerFailed", ...
        "Failed to preprocess protein with pdbfixer: %s\n" + ...
        "Please check if pdbfixer is installed: emk.setup.validate()", ...
        ME.message);
end
%[text] \--- 1c: Generate Receptor PDBQT with meeko mk\_prepare\_receptor ---
%[text] 
%[text] ### The PDBQT Format
%[text] AutoDock Vina takes input in PDBQT format, an extension of PDB that adds partial charges and AutoDock atom-type labels to each atom.
%[text] `meeko`'s `mk_prepare_receptor` computes Gasteiger partial charges (semi-empirical values derived from atomic electronegativity) and writes the formatted output automatically.
logInfo("R08: Generating receptor PDBQT...");
try
    % Add ~/.local/bin to PATH (location of scripts when installed with pip --user)
    % Note API change in meeko 0.7.x:
    %   Old: mk_prepare_receptor -i <pdb>  -> Changed to --read_with_prody in meeko 0.7
    %   New: mk_prepare_receptor --read_pdb <pdb> -p <pdbqt>  (prody not needed)
    prepareCmd = sprintf( ...
        "export PATH=$HOME/.local/bin:$PATH && mk_prepare_receptor --read_pdb '%s' -p '%s'", ...
        preparedPdb, receptorPdbqt);
    [status, cmdOut] = system(prepareCmd);
    if status ~= 0
        % Fallback: directly call meeko.cli module with python3 -c
        % (robust as it doesn't depend on shell script presence or PATH settings)
        logWarn("R08: mk_prepare_receptor failed -> Retrying with Python module execution...");
        pySys = py.importlib.import_module("sys");
        pySub = py.importlib.import_module("subprocess");
        % Pass paths via sys.argv (safe for paths with single quotes or spaces)
        % When using python3 -c "code" <pdb> <pdbqt>, sys.argv = ['-c', pdb, pdbqt]
        % Convert to 1×N char with strjoin (char() on string array creates 2D matrix)
        pyInlineCode = char(strjoin([ ...
            "import sys;", ...
            "sys.argv=['mk_prepare_receptor','--read_pdb',sys.argv[1],'-p',sys.argv[2]];", ...
            "from meeko.cli.mk_prepare_receptor import main;main()"], ""));
        result = pySub.run( ...
            py.list({char(pySys.executable), '-c', pyInlineCode, ...
                     char(preparedPdb), char(receptorPdbqt)}), ...
            capture_output=py.True, text=py.True);
        if int32(result.returncode) ~= 0
            error("emk:r08:receptorPdbqtFailed", ...
                "mk_prepare_receptor failed:\n%s\n%s", ...
                char(result.stdout), char(result.stderr));
        end
    end
    logInfo("R08: Receptor PDBQT generation complete -> %s", receptorPdbqt);
catch ME
    error("emk:r08:receptorPdbqtFailed", ...
        "Failed to generate receptor PDBQT: %s\n" + ...
        "Please check if meeko is installed: emk.setup.validate()", ...
        ME.message);
end

logInfo("R08: Section 1 complete");
%%
%[text] ## Section 2: Obtain SMILES for 10 Ligands and Generate 3D Structures
%[text] ### Ligand 3D Structure Generation Pipeline
%[text] Docking requires 3D coordinates of ligands. Starting from SMILES (2D topology), generate 3D structures using the following steps.
%[text] 1. Generate a Mol object with hydrogens added using `AddHs`.
%[text] 2. Calculate initial 3D coordinates using `EmbedMolecule` (ETKDGv3 distance geometry method).
%[text] 3. Perform force field optimization using `UFFOptimizeMolecule` (Universal Force Field) to converge to a strain-free structure.
%[text] 4. Use `meeko MoleculePreparation` to assign Gasteiger charges and convert to PDBQT format (Vina input). \
%[text] 
%[text] \--- 2a: Ligand Definition (PubChem CID + SMILES Notes) ---
%[text] SMILES Source: PubChem REST API (normalized via webread)
%[text] CID is the official PubChem ID (immutable). SMILES will be obtained in Section 2b.
ligandDefs = struct( ...
    "name",    {"Nirmatrelvir",   "Ensitrelvir",  "Boceprevir", ...
                "Lopinavir",      "Ritonavir",    "Nelfinavir", ...
                "Remdesivir",     "Favipiravir",  "Oseltamivir", ...
                "Acyclovir"}, ...
    "pubchemQuery", ...
               {"Nirmatrelvir",   "Ensitrelvir",  "Boceprevir", ...
                "Lopinavir",      "Ritonavir",    "Nelfinavir", ...
                "Remdesivir",     "Favipiravir",  "Oseltamivir", ...
                "Acyclovir"}, ...
    "target",  {"SARS-CoV-2 Mpro",  "SARS-CoV-2 Mpro", "HCV NS3/4A protease", ...
                "HIV protease",      "HIV protease",     "HIV protease", ...
                "RNA polymerase",    "RNA polymerase",   "Neuraminidase", ...
                "Thymidine kinase"}, ...
    "expectation", ...
               {"HIGH",   "HIGH",  "MODERATE-HIGH", ...
                "LOW-MOD", "LOW",  "LOW", ...
                "VERY-LOW", "VERY-LOW", "VERY-LOW", ...
                "VERY-LOW"} ...
);
nLigands = numel(ligandDefs);
logInfo("R08: Processing %d ligands", nLigands);
%[text] \--- 2b: Obtain canonical SMILES via PubChem REST API ---
logInfo("R08: Retrieving canonical SMILES from PubChem...");
ligandSmiles = strings(1, nLigands);
for k = 1:nLigands
    query = ligandDefs(k).pubchemQuery;
    try
        % PubChem PUG REST: Name → canonical SMILES
        apiUrl = sprintf( ...
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/CanonicalSMILES/TXT", ...
            urlencode(query));
        smi = strtrim(webread(apiUrl));
        % Use only the first line if webread returns multiple lines
        lines = strsplit(smi, newline);
        ligandSmiles(k) = strtrim(lines{1});
        logInfo("  [%d/%d] %-16s -> %s", k, nLigands, query, ...
            extractBefore(ligandSmiles(k) + " ", min(50, strlength(ligandSmiles(k))+1)));
    catch ME
        logWarn("  [%d/%d] %-16s -> PubChem retrieval failed: %s", k, nLigands, query, ME.message);
    end
end
%[text] \--- 2c: SMILES → 3D Structure Generation (RDKit EmbedMolecule + UFF) ---
logInfo("R08: Generating 3D coordinates with RDKit...");
%[text] Execute 3D generation as a subprocess (same fallback pattern as Section 1)
%[text] Reason: RDKit's AddHs returns an RWMol (read-write mol). MATLAB fails to map RWMol's Boost.Python metaclass to `py.Boost.Python.class`, so it's more robust to complete within the Python process.
pySys = py.importlib.import_module("sys");
pySub = py.importlib.import_module("subprocess");
pythonExe = char(pySys.executable);
%[text] Python Inline: argv\[1\]=SMILES, argv\[2\]=Output SDF Path
%[text] ETKDGv3 + Fixed Seed + UFF Optimization → Hydrogen Removed SDF
pyGen3dCode = char(strjoin([ ...
    "from rdkit import Chem; from rdkit.Chem import AllChem; import sys;", ...
    "smi,out=sys.argv[1],sys.argv[2];", ...
    "mol=Chem.MolFromSmiles(smi);", ...
    "assert mol is not None,'invalid SMILES: '+smi;", ...
    "mol=Chem.AddHs(mol);", ...
    "p=AllChem.ETKDGv3(); p.randomSeed=42;", ...
    "r=AllChem.EmbedMolecule(mol,p);", ...
    "assert r>=0,'EmbedMolecule failed (r='+str(r)+')';", ...
    "AllChem.UFFOptimizeMolecule(mol,maxIters=500);", ...
    "mol=Chem.RemoveHs(mol);", ...
    "w=Chem.SDWriter(out); w.write(mol); w.close()"], " "));

ligandSdfFiles = strings(1, nLigands);
ligandValid    = false(1, nLigands);

for k = 1:nLigands
    name = ligandDefs(k).name;
    smi  = ligandSmiles(k);
    if smi == ""
        logWarn("  [%02d] %s: No SMILES -- Skipping", k, name);
        continue;
    end

    try
        sdfPath = fullfile(ligandDir, sprintf("%02d_%s.sdf", k, name));

        % Generate 3D structure with Python subprocess
        result = pySub.run( ...
            py.list({pythonExe, '-c', pyGen3dCode, char(smi), char(sdfPath)}), ...
            capture_output=py.True, text=py.True);

        if int32(result.returncode) ~= 0
            error("emk:r08:embedFailed", "3D generation failed:\n%s", char(result.stderr));
        end

        ligandSdfFiles(k) = sdfPath;
        ligandValid(k)    = true;

        mol2d = emk.mol.fromSmiles(smi);   % For MW calculation
        mw = emk.descriptor.molWeight(mol2d);
        logInfo("  [%02d] %-15s: 3D generation complete (MW=%.1f, %s)", ...
            k, name, mw, sdfPath);

    catch ME
        logWarn("  [%02d] %s: Error -> %s", k, name, ME.message);
    end
end

nValid = sum(ligandValid);
logInfo("R08: 3D structure generation complete -- %d / %d successful", nValid, nLigands);
if nValid == 0
    error("emk:r08:noValidLigands", ...
        "No valid ligand 3D structures could be generated. Please check Section 0.");
end
logInfo("R08: Section 2 complete");
%%
%[text] ## Section 3: PDBQT Conversion (meeko MoleculePreparation)
%[text] ### Why PDBQT?
%[text] AutoDock Vina requires each atom to carry an AutoDock atom-type label and a partial charge.
%[text] `meeko` handles this automatically by mapping RDKit atom types and computing Gasteiger charges, writing the result as a valid PDBQT file.
%[text] As an LGPL-2.1 library, meeko does not impose copyleft obligations on EasyMolKit when accessed via dynamic linking in Python (see docs/compliance.md CL-7).
%[text] 
%[text] PDBQT conversion is also executed via subprocess (same pattern as Section 2).
%[text] Reason: Subscript access of `SDMolSupplier` ({}) returns RWMol (Boost.Python) and MATLAB fails. Also, due to API changes in meeko 0.7.x where the return type of `MoleculePreparation.prepare()` changes, it is more robust to complete within the Python process.
pySys2 = py.importlib.import_module("sys");
pySub2 = py.importlib.import_module("subprocess");
pythonExe2 = char(pySys2.executable);
%[text] Python inline: argv\[1\]=SDF path, argv\[2\]=PDBQT output path
%[text] Compatible with both meeko 0.5.x / 0.7.x (check if the return value of prepare / write\_string is a tuple)
pyMeekoCode = char(strjoin([ ...
    "import sys; from rdkit import Chem; from rdkit.Chem import SDMolSupplier;", ...
    "from meeko import MoleculePreparation,PDBQTWriterLegacy;", ...
    "sdf,out=sys.argv[1],sys.argv[2];", ...
    "suppl=SDMolSupplier(sdf,removeHs=False);", ...
    "mol=next(iter(suppl));", ...
    "assert mol is not None,'SDF read failed: '+sdf;", ...
    "mol=Chem.AddHs(mol,addCoords=True);", ...
    "prep=MoleculePreparation();", ...
    "res=prep.prepare(mol);", ...
    "setup=res[0] if isinstance(res,list) else list(prep.setup_dict.values())[0];", ...
    "assert setup is not None,'meeko setup empty';", ...
    "wr=PDBQTWriterLegacy.write_string(setup);", ...
    "pdbqt=wr[0] if isinstance(wr,tuple) else wr;", ...
    "open(out,'w').write(pdbqt)"], " "));

ligandPdbqtFiles = strings(1, nLigands);

logInfo("R08: Converting to PDBQT with meeko...");
for k = 1:nLigands
    if ~ligandValid(k)
        continue;
    end
    name    = ligandDefs(k).name;
    sdfPath = ligandSdfFiles(k);

    try
        pdbqtPath = fullfile(ligandDir, sprintf("%02d_%s.pdbqt", k, name));

        result2 = pySub2.run( ...
            py.list({pythonExe2, '-c', pyMeekoCode, char(sdfPath), char(pdbqtPath)}), ...
            capture_output=py.True, text=py.True);

        if int32(result2.returncode) ~= 0
            error("emk:r08:meekoPdbqtFailed", "meeko failed:\n%s", char(result2.stderr));
        end

        ligandPdbqtFiles(k) = pdbqtPath;
        logInfo("  [%02d] %-15s: PDBQT generation complete", k, name);

    catch ME
        logWarn("  [%02d] %s: meeko error -> %s", k, name, ME.message);
        ligandValid(k) = false;
    end
end

nValid = sum(ligandValid);
logInfo("R08: PDBQT conversion complete -- %d / %d successful", nValid, nLigands);
logInfo("R08: Section 3 complete");
%%
%[text] ## Section 4: Docking Execution (AutoDock Vina Python API)
%[text] ### AutoDock Vina Scoring Function
%[text] Vina explores the ligand's conformational space via Iterated Local Search, minimizing an empirical free-energy function: $\\Delta G \\approx w_1 \\cdot \\mathrm{gauss_1} + w_2 \\cdot \\mathrm{gauss_2} + w_3 \\cdot \\mathrm{repulsion} + w_4 \\cdot \\mathrm{hydrophobic} + w_5 \\cdot \\mathrm{H\_bond}$ (Trott & Olson 2010).
%[text] Scores are reported in kcal/mol; more negative values indicate stronger predicted binding affinity.
%[text] As a rough guide, a score of $-7$ kcal/mol or below is generally considered a promising lead, while $-9$ kcal/mol or below suggests strong binding.
%[text]
%[text] Note the two key approximations made in this exercise:
%[text] - Nirmatrelvir forms a covalent bond with Cys145 in reality; here we use a non-covalent approximation, which underestimates its contribution.
%[text] - Protein flexibility is ignored (rigid receptor); induced-fit effects are not accounted for. \
vina = py.importlib.import_module("vina");

bindingEnergies  = nan(1, nLigands);    % best pose [kcal/mol]
dockingPoseFiles = strings(1, nLigands);

logInfo("R08: Docking started (exhaustiveness=%d)...", EXHAUSTIVENESS);
logInfo("R08: Estimated time required %d ~ %d minutes", nValid * 1, nValid * 3);

for k = 1:nLigands
    if ~ligandValid(k)
        continue;
    end
    name        = ligandDefs(k).name;
    pdbqtPath   = ligandPdbqtFiles(k);
    poseOutPath = fullfile(poseDir, sprintf("%02d_%s_poses.pdbqt", k, name));

    logInfo("  [%02d/%02d] %s -- Docking...", k, nLigands, name);
    tic;

    try
        % Create Vina instance (CPU thread count is auto-configured)
        v = vina.Vina(sf_name="vina", verbosity=py.int(0));

        % Set receptor and ligand
        v.set_receptor(receptorPdbqt);
        v.set_ligand_from_file(pdbqtPath);

        % Set binding box
        % Use py.float() for explicit conversion as num2cell cannot convert to Python scalar
        vinaCenter  = py.list({py.float(BOX_CENTER(1)), py.float(BOX_CENTER(2)), py.float(BOX_CENTER(3))});
        vinaBoxSize = py.list({py.float(BOX_SIZE(1)),   py.float(BOX_SIZE(2)),   py.float(BOX_SIZE(3))});
        v.compute_vina_maps(center=vinaCenter, box_size=vinaBoxSize);

        % Execute docking
        v.dock(exhaustiveness=py.int(EXHAUSTIVENESS), n_poses=py.int(N_POSES));

        % Retrieve results
        % v.energies() returns a numpy ndarray (n_poses x n_terms)
        % Convert to 1D with flatten() and then to MATLAB double
        pyEn     = v.energies(n_poses=py.int(N_POSES));
        energies = reshape(double(pyEn.flatten()), N_POSES, []);
        bindingEnergies(k) = energies(1, 1);   % best pose total score [kcal/mol]

        % Save poses in PDBQT
        v.write_poses(poseOutPath, n_poses=py.int(N_POSES), overwrite=true);
        dockingPoseFiles(k) = poseOutPath;

        elapsed = toc;
        logInfo("    -> Score: %.2f kcal/mol (%.1f seconds)", bindingEnergies(k), elapsed);

    catch ME
        toc;
        logWarn("  [%02d] %s: Docking error -> %s", k, name, ME.message);
    end
end

logInfo("R08: All docking completed");
logInfo("R08: Section 4 completed");
%%
%[text] ## Section 5: Binding Energy Comparison and Discussion
%[text] \--- 5a: Create Result Table ---
names_all     = string({ligandDefs.name})';
targets_all   = string({ligandDefs.target})';
expect_all    = string({ligandDefs.expectation})';
energies_col  = bindingEnergies';
valid_col     = ligandValid';

resultTbl = table(names_all, targets_all, expect_all, energies_col, valid_col, ...
    VariableNames=["Compound", "OriginalTarget", "Expected", "DeltaG_kcal_mol", "Success"]);
%[text] Sort by score (most negative = strongest binding)
validRows = resultTbl.Success;
scored    = resultTbl(validRows, :);
scored    = sortrows(scored, "DeltaG_kcal_mol", "ascend");

%[text] ### Docking Result Summary (by Score)
logInfo("%-20s  %-25s  %-8s  %s", "Compound", "Original Target", "Expected", "Score (kcal/mol)");
logInfo("%s", repmat("-", 1, 75));
for i = 1:height(scored)
    logInfo("%-20s  %-25s  %-8s  %.2f", ...
        scored.Compound(i), scored.OriginalTarget(i), ...
        scored.Expected(i), scored.DeltaG_kcal_mol(i));
end
%[text] \--- 5b: CSV Export ---
csvPath = fullfile(runDir, "r08_docking_results.csv");
writetable(resultTbl, csvPath);
logInfo("R08: Result CSV -> %s", csvPath);
%[text] \--- 5c: Horizontal Bar Graph (barh) ---
%[text] 
%[text] ### Why Covalent Inhibitors Score Lower Than Expected
%[text] Nirmatrelvir and Boceprevir both rank lower here than their experimental potency would suggest.
%[text] The reason is that Vina scores only non-covalent interactions and cannot capture the additional stabilization from covalent bond formation (+2-4 kcal/mol).
%[text] - Nirmatrelvir: the nitrile warhead undergoes Michael addition with Cys145-SH.
%[text] - Boceprevir: the ketoamide forms a nucleophilic addition with Cys145-SH (Ma et al. 2020; Mpro Ki ~9 uM). \
%[text] Conversely, high scores for HIV drugs such as Lopinavir reflect a coincidental shape match with the active site pocket — they showed no benefit in standalone COVID-19 clinical trials (Cao et al. 2020).
fig1 = figure("Name", "R08: Docking Scores vs Mpro (6LU7)", "NumberTitle", "off");
set(fig1, "Position", [100, 100, 800, 500]);

cmap = [0.85, 0.33, 0.10;    % orange-red  = Mpro inhibitors
        0.85, 0.33, 0.10;
        0.93, 0.69, 0.13;    % yellow      = cross-reactive
        0.47, 0.67, 0.19;    % green       = off-target (HIV)
        0.47, 0.67, 0.19;
        0.47, 0.67, 0.19;
        0.30, 0.75, 0.93;    % blue        = RNA pol / other
        0.30, 0.75, 0.93;
        0.30, 0.75, 0.93;
        0.30, 0.75, 0.93];

nScored = height(scored);
scoreVals = scored.DeltaG_kcal_mol;
compNames = scored.Compound;
%[text] Rearrange colors according to rank
barColors = zeros(nScored, 3);
for i = 1:nScored
    idx = find(strcmp(names_all, compNames(i)), 1);
    if ~isempty(idx) && idx <= size(cmap, 1)
        barColors(i, :) = cmap(idx, :);
    else
        barColors(i, :) = [0.5, 0.5, 0.5];
    end
end

bh = barh(scoreVals);
bh.FaceColor = "flat";
for i = 1:nScored
    bh.CData(i, :) = barColors(i, :);
end

yticklabels(compNames);
xlabel("Binding Free Energy (kcal/mol)", "FontSize", 12);
title(sprintf("AutoDock Vina Scores vs SARS-CoV-2 Mpro (PDB: 6LU7)\nexhaustiveness=%d", ...
    EXHAUSTIVENESS), "FontSize", 12);
grid on;
ax = gca;
ax.GridAlpha = 0.3;
%[text] Add general guideline lines
xline(-7.0, "--", "Promising lead (-7)", "Color", [0.5, 0, 0.5], ...
    "LabelHorizontalAlignment", "left", "FontSize", 9);
xline(-9.0, "--", "Strong binding (-9)", "Color", [0.8, 0, 0], ...
    "LabelHorizontalAlignment", "left", "FontSize", 9);
%[text] Add numeric labels to each bar
for i = 1:nScored
    text(scoreVals(i) - 0.2, i, sprintf("%.1f", scoreVals(i)), ...
        "HorizontalAlignment", "right", "FontSize", 9, "Color", "white", "FontWeight", "bold");
end

%[text] Add legend (color and category correspondence)
hold on;
hLeg = [ ...
    patch(NaN, NaN, [0.85, 0.33, 0.10]); ...
    patch(NaN, NaN, [0.93, 0.69, 0.13]); ...
    patch(NaN, NaN, [0.47, 0.67, 0.19]); ...
    patch(NaN, NaN, [0.30, 0.75, 0.93])];
legend(hLeg, ["Mpro inhibitors", "Cross-reactive (HCV)", "HIV protease", "RNA pol / other"], ...
    "Location", "southeast", "FontSize", 9);
hold off;

plotPath = fullfile(runDir, "r08_docking_scores.png");
saveas(fig1, plotPath);
logInfo("R08: Score Comparison Plot -> %s", plotPath);
%[text] \--- 5d: Discussion ---
if height(scored) == 0
    logWarn("R08: No docking results. Check errors in Section 3-4.");
else
    bestCompound  = scored.Compound(1);
    bestScore     = scored.DeltaG_kcal_mol(1);
    worstCompound = scored.Compound(end);
    worstScore    = scored.DeltaG_kcal_mol(end);

    logInfo("R08: Discussion");
    logInfo("  Best Score  : %s (%.2f kcal/mol)", bestCompound, bestScore);
    logInfo("  Worst Score : %s (%.2f kcal/mol)", worstCompound, worstScore);

    % Compare scores of Nirmatrelvir and Boceprevir
    nirmIdx = strcmp(scored.Compound, "Nirmatrelvir");
    bocepIdx = strcmp(scored.Compound, "Boceprevir");
    if any(nirmIdx) && any(bocepIdx)
        nirmScore  = scored.DeltaG_kcal_mol(nirmIdx);
        bocepScore = scored.DeltaG_kcal_mol(bocepIdx);
        logInfo("  Nirmatrelvir vs Boceprevir: %.2f vs %.2f kcal/mol", nirmScore, bocepScore);
        logInfo("  [Note] Both compounds form covalent bonds with Cys145 -- underestimated by Vina's non-covalent approximation.");
        if bocepScore < nirmScore
            logInfo("  -> Boceprevir scores stronger than Nirmatrelvir: consistent with ketoamide cross-reactivity.");
        else
            logInfo("  -> Boceprevir scores lower than Nirmatrelvir: underestimation due to non-covalent approximation.");
            logInfo("     Experimental value: Boceprevir Mpro Ki ~9 uM (Ma et al. 2020) confirms cross-reactivity.");
        end
    end

    % Supplement if Ritonavir scores higher than Nirmatrelvir
    ritoIdx = strcmp(scored.Compound, "Ritonavir");
    if any(ritoIdx) && any(nirmIdx)
        ritoScore = scored.DeltaG_kcal_mol(ritoIdx);
        if ritoScore < nirmScore
            logInfo("  Ritonavir (%.2f) > Nirmatrelvir (%.2f): HIV drug scores higher than Mpro inhibitor.", ...
                ritoScore, nirmScore);
            logInfo("  -> Result of shape coincidence + randomness of exhaustiveness=4.");
            logInfo("     Ritonavir belongs to protease inhibitor family with similar active site shape,");
            logInfo("     but showed no effect in standalone COVID clinical trials (Cao et al. 2020).");
            logInfo("     A typical example of score and clinical activity divergence.");
        end
    end

    % Check MW and score of Favipiravir
    favIdx = strcmp(scored.Compound, "Favipiravir");
    if any(favIdx)
        favScore = scored.DeltaG_kcal_mol(favIdx);
        logInfo("  Favipiravir (MW=157): %.2f kcal/mol -- %s", favScore, ...
            ternary_(favScore > -5.5, "Expected low score (does not fit pocket)", ...
                "Higher than expected score -- possibly partially accommodated in shallow active site region (randomness of exhaustiveness=4 also affects)"));
    end
end
logInfo("R08: Section 5 Complete");
%%
%[text] ## Section 6: 3D Visualization -- Active Site and Binding Pose
%[text] ### Purpose of 3D Visualization
%[text] Numerical scores alone cannot explain *why* a compound binds. 3D visualization lets you see directly how the docking pose sits within the active site — including the spatial relationship to catalytic residues (His41, Cys145, Glu166), how well the ligand fills the hydrophobic P2 pocket (Met49, Met165), and whether the compound is deeply buried or merely perching at the pocket edge.
%[text]
%[text] **Viewer Controls (Interactive):** Left-drag to rotate, scroll to zoom, right-drag to pan.
%[text] The initial view is centered on the docking compound.
%[text] Scroll out to see the full protein ribbon.
%[text]
%[text] ### Interpreting the Visualization: What the Calculation Cannot Tell You
%[text] These poses are computational models; they may differ substantially from experimentally determined binding modes.
%[text] **Rigid Receptor**: The protein is held fixed throughout. Real binding events involve induced-fit movements that can shift key side chains.
%[text] **Score Accuracy**: Vina scores are empirical approximations with an absolute uncertainty of roughly +/-2 kcal/mol; they are most reliable for *relative* comparisons within a single run.
%[text] **exhaustiveness=4**: This setting (chosen for classroom speed) samples conformation space sparsely. Use 8 or higher for publication-quality results.
%[text]
%[text] **Reading the Display:**
%[text] - Ribbons (gold = beta sheet / magenta = alpha helix / gray = loop) show the secondary structure of the Mpro backbone. The gold beta-barrel (Domains I/II) forms the active site; the magenta helix bundle (Domain III) is the regulatory domain.
%[text] - Thin stick clusters (8 groups) highlight the side chains of key active-site residues (His41, Met49, Cys145, His163, Met165, Glu166, Asp187, Gln189).
%[text] - CPK element colors: white/gray = carbon, red = oxygen, blue = nitrogen, yellow = sulfur (Met, Cys).
%[text] - Thick sticks at the center mark the docking compound (best pose only). Check whether the ligand is buried deep in the pocket.
%[text] - Only the top-scored pose (MODEL 1) is displayed; atoms after the first ENDMDL are excluded to prevent the molecule from appearing duplicated. \
%[text] 
%[text] This exercise aims to grasp qualitative trends. For quantitative predictions and actual drug development, experimental validation such as IC50 measurement and crystal structure analysis is necessary.
%[text] 
%[text] \--- 6a: Retrieve Cα Coordinates from Prepared Receptor PDB ---
%[text] Parse the PDB file line by line (Base MATLAB only, Bioinformatics Toolbox not required).
logInfo("R08: Parsing Cα coordinates from PDB file...");

caCoords = [];   % Cα coordinates [N x 3]
caResIDs  = [];  % Residue numbers

fh = fopen(preparedPdb, "r");
if fh < 0
    logWarn("R08: %s not found -- skipping visualization", preparedPdb);
else
    while ~feof(fh)
        line = fgetl(fh);
        if ~ischar(line) || length(line) < 54
            continue;
        end
        recType  = strtrim(line(1:6));
        atomName = strtrim(line(13:16));
        if (strcmp(recType, "ATOM") || strcmp(recType, "HETATM")) && ...
                strcmp(atomName, "CA")
            x = str2double(line(31:38));
            y = str2double(line(39:46));
            z = str2double(line(47:54));
            resNum = str2double(strtrim(line(23:26)));
            if ~isnan(x) && ~isnan(y) && ~isnan(z)
                caCoords = [caCoords; x, y, z];  %#ok<AGROW>
                caResIDs  = [caResIDs;  resNum];  %#ok<AGROW>
            end
        end
    end
    fclose(fh);
    logInfo("R08: Loaded %d Cα atoms", size(caCoords, 1));
end
%[text] \--- 6b: Define Active Site Residues ---
activeSiteResNums = [41, 145, 163, 164, 165, 166, 168, 172, 187, 189];
activeSiteLabels  = ["His41", "Cys145", "Pro163", "Leu164", "Met165", ...
                     "Glu166", "Ala168", "Gln172", "Asp187", "Gln189"];
activeCoords = [];
activeLabels = strings(0);

for i = 1:numel(activeSiteResNums)
    idx = find(caResIDs == activeSiteResNums(i), 1);
    if ~isempty(idx)
        activeCoords = [activeCoords; caCoords(idx, :)];  %#ok<AGROW>
        activeLabels = [activeLabels; activeSiteLabels(i)];  %#ok<AGROW>
    end
end
logInfo("R08: Detected %d / %d active site residues", size(activeCoords, 1), numel(activeSiteResNums));
%[text] \--- 6c: Select One Best Score Ligand (for Ribbon Visualization) ---
displayLigandName = "";
validPoseIdx = find(dockingPoseFiles ~= "");
if ~isempty(validPoseIdx)
    [~, bestIdx]     = min(bindingEnergies(validPoseIdx));
    displayLigandIdx  = validPoseIdx(bestIdx);
    displayLigandName = ligandDefs(displayLigandIdx).name;
    displayPoseFile   = dockingPoseFiles(displayLigandIdx);
    logInfo("R08: Ribbon display target: %s (best score)", displayLigandName);
end
%[text] \--- 6d: 3Dmol.js Ribbon Visualization (Displayed in uifigure) ---
%[text] 
%[text] ### Self-contained HTML Generation
%[text] `py3Dmol`'s `write_html()` loads 3Dmol.js from an external CDN, which fails in MATLAB Online's sandbox browser.
%[text] Therefore, `urllib.request` is used to download 3Dmol.js directly and embed it in the HTML.
logInfo("R08: Generating 3Dmol.js ribbon visualization (self-contained HTML)...");

pySub6     = py.importlib.import_module("subprocess");
pySys6     = py.importlib.import_module("sys");
pythonExe6 = char(pySys6.executable);

htmlPath  = fullfile(runDir, "r08_ribbon.html");
vizScript = fullfile(runDir, "viz_ribbon.py");
%[text] Generate Python Script
%[text] Uses only stdlib (urllib.request + json), no external packages required.
%[text] Downloads 3Dmol.js from CDN and embeds it inline in HTML.
%[text] PDB data is safely escaped with json.dumps() and passed to JS variables.
fid = fopen(char(vizScript), 'w');
fprintf(fid, "import sys, urllib.request, json\n");
fprintf(fid, "rec, lig, out = sys.argv[1], sys.argv[2], sys.argv[3]\n");
fprintf(fid, "rec_pdb = open(rec).read()\n");
fprintf(fid, "lig_pdb = ''\n");
fprintf(fid, "if lig != 'none':\n");
fprintf(fid, "    _ll = []\n");
fprintf(fid, "    for _l in open(lig):\n");
fprintf(fid, "        if _l.startswith('ENDMDL'): break\n");
fprintf(fid, "        if _l[:6] in ('ATOM  ','HETATM'): _ll.append(_l)\n");
fprintf(fid, "    if not _ll:  # fallback: no MODEL/ENDMDL tags\n");
fprintf(fid, "        _ll = [_l for _l in open(lig) if _l[:6] in ('ATOM  ','HETATM')]\n");
fprintf(fid, "    lig_pdb = ''.join(_ll)\n");
fprintf(fid, "js3d = ''\n");
fprintf(fid, "for u in ['https://3dmol.org/build/3Dmol-min.js',\n");
fprintf(fid, "          'https://3dmol.csb.pitt.edu/build/3Dmol-min.js']:\n");
fprintf(fid, "    try:\n");
fprintf(fid, "        req = urllib.request.Request(u, headers={'User-Agent': 'Mozilla/5.0'})\n");
fprintf(fid, "        js3d = urllib.request.urlopen(req, timeout=20).read().decode('utf-8', 'replace')\n");
fprintf(fid, "        break\n");
fprintf(fid, "    except: pass\n");
fprintf(fid, "if not js3d:\n");
fprintf(fid, "    print('ERROR: 3Dmol.js download failed from both CDNs', file=sys.stderr); sys.exit(1)\n");
fprintf(fid, "rec_js = json.dumps(rec_pdb)\n");
fprintf(fid, "lig_js = json.dumps(lig_pdb)\n");
fprintf(fid, "has_lig = 'true' if lig_pdb else 'false'\n");
fprintf(fid, "jv_lines = [\n");
fprintf(fid, "    'var v=$3Dmol.createViewer(document.getElementById(""v""),{backgroundColor:""#1a1a2e""});',\n");
fprintf(fid, "    'v.addModel(' + rec_js + ',""pdb"");',\n");
fprintf(fid, "    'v.setStyle({model:0,chain:""A""},{cartoon:{colorscheme:""ssJmol"",opacity:0.85}});',\n");
fprintf(fid, "    '// Active site residues (Chain A only): element-colored sticks to show binding pocket',\n");
fprintf(fid, "    '[41,49,145,163,165,166,187,189].forEach(function(r){',\n");
fprintf(fid, "    '    v.addStyle({model:0,resi:r,chain:""A""},{stick:{colorscheme:""element"",radius:0.18}});',\n");
fprintf(fid, "    '});',\n");
fprintf(fid, "    'if(' + has_lig + '){',\n");
fprintf(fid, "    '    v.addModel(' + lig_js + ',""pdb"");',\n");
fprintf(fid, "    '    v.setStyle({model:1},{stick:{colorscheme:""element"",radius:0.25}});',\n");
fprintf(fid, "    '    v.zoomTo({model:1});',\n");
fprintf(fid, "    '} else { v.zoomTo({model:0}); }',\n");
fprintf(fid, "    'v.render();'\n");
fprintf(fid, "]\n");
fprintf(fid, "jv = '\\n'.join(jv_lines)\n");
fprintf(fid, "css = ('html,body{margin:0;padding:0;width:100%%;height:100%%;overflow:hidden;background:#1a1a2e;}'\n");
fprintf(fid, "       '#v{width:100%%;height:100%%;display:block;position:absolute;top:0;left:0;}')\n");
fprintf(fid, "html = ('<html><head>'\n");
fprintf(fid, "        '<style>' + css + '</style>'\n");
fprintf(fid, "        '<script type=""text/javascript"">\\n' + js3d + '\\n</script>'\n");
fprintf(fid, "        '</head>'\n");
fprintf(fid, "        '<body>'\n");
fprintf(fid, "        '<div id=""v""></div>'\n");
fprintf(fid, "        '<script type=""text/javascript"">\\n' + jv + '\\n</script>'\n");
fprintf(fid, "        '</body></html>')\n");
fprintf(fid, "open(out, 'w').write(html)\n");
fclose(fid);
%[text] Select Ligand File (Use if Docked, 'none' if Not Completed)
if displayLigandName ~= "" && exist('displayPoseFile', 'var') && isfile(char(displayPoseFile))
    ligArg = char(displayPoseFile);
else
    ligArg = 'none';
end
%[text] Generate HTML (Retrieve 3Dmol.js from CDN and Embed Inline)
vizRes = pySub6.run(py.list({pythonExe6, char(vizScript), ...
    char(preparedPdb), ligArg, char(htmlPath)}), ...
    capture_output=py.True, text=py.True);

if int32(vizRes.returncode) == 0
    logInfo("R08: Ribbon HTML saved -> %s", htmlPath);
    % Display ribbon in MATLAB Figure panel using uifigure + uihtml
    % (uifigure uses the same WebGL-compatible frame as App Designer)
    try
        figW = 960; figH = 780;
        figRibbon = uifigure( ...
            "Name",     sprintf("R08: Ribbon -- %s vs Mpro 6LU7", displayLigandName), ...
            "Position", [50, 50, figW, figH]);
        % Align uihtml Position with uifigure for full display
        uh = uihtml(figRibbon, "HTMLSource", char(fullfile(pwd, htmlPath)));
        uh.Position = [0, 0, figW, figH];
    catch ME2
        logWarn("R08: uifigure display failed (%s) -- falling back to web()", ME2.message);
        try
            web(char(fullfile(pwd, htmlPath)));
        catch
            logWarn("R08: Please open '%s' manually from the Files panel", htmlPath);
        end
    end
else
    logWarn("R08: Ribbon generation failed (returncode=%d)", int32(vizRes.returncode));
    if strlength(string(char(vizRes.stderr))) > 0
        logWarn("R08: stderr:\n%s", char(vizRes.stderr));
    end
    if strlength(string(char(vizRes.stdout))) > 0
        logWarn("R08: stdout:\n%s", char(vizRes.stdout));
    end
end

logInfo("R08: Section 6 complete");
%[text] ### How to Read the Visualization
%[text] As you explore the ribbon diagram, focus on the following:
%[text] 1. **Docking Compound Position**: Are the thick sticks seated deep inside the beta-barrel pocket? A ligand perching at the rim rather than buried inside suggests a shallow or misdocked pose.
%[text] 2. **Active-Site Proximity**: Are the thin sticks for His41 (general base) and Cys145 (nucleophile) close to the docking compound? Covalent inhibitors such as Nirmatrelvir should sit extremely close to Cys145.
%[text] 3. **P2 Pocket Filling**: Are the yellow sulfur sticks (Met49, Met165) adjacent to the hydrophobic portion of the ligand? Good P2 pocket occupancy correlates strongly with a high docking score.
%[text] 4. **Rotate and Explore**: Drag to rotate and examine pocket depth and orientation from multiple angles. \
%%
%[text] ## Section 7: Comprehensive Discussion
%[text] ### Rationale for the Experimental Design
%[text] The three compound groups let us examine the gap between docking scores and clinical outcomes at each level.
%[text] - **Group A — Designed Mpro Inhibitors** (Nirmatrelvir, Ensitrelvir): We expect high (very negative) scores.
%[text] - **Group B — Cross-Reactive and Clinically Failed Agents** (Boceprevir through Nelfinavir): We expect medium to low scores, with interesting exceptions.
%[text] - **Group C — Off-Target Drugs** (Remdesivir through Acyclovir): We expect low scores, reflecting poor geometric complementarity with the Mpro pocket. \
%[text] ### Limitations of Docking Scores
%[text] - Scores approximate binding free energy; they do not always correlate with in vivo activity.
%[text] - The rigid receptor approximation ignores induced-fit effects, which can be significant in flexible active sites.
%[text] - Covalent inhibitors such as Nirmatrelvir are systematically underestimated, because the covalent warhead contribution is not captured by the scoring function. \
%[text] ### Lessons from This Drug Repurposing Case Study
%[text] - Within the cysteine protease family, conserved active-site geometry can produce unexpected cross-reactivity — Boceprevir is a textbook example.
%[text] - Even favorable docking scores do not guarantee success: ADMET liabilities or poor PK can eliminate a compound entirely.
%[text] - Lopinavir's clinical failure reflects a combination of marginal Mpro affinity and CYP3A4 autoinduction that reduces its effective plasma concentration. \
%[text] ### Final Score List (Console Output)
logInfo("%-3s  %-20s  %-8s  %s", "#", "Compound", "Score", "Assessment");
logInfo("%s", repmat("-", 1, 60));
for i = 1:height(scored)
    sc = scored.DeltaG_kcal_mol(i);
    if sc <= -9
        assess = "*** STRONG binding";
    elseif sc <= -7
        assess = "**  Promising lead";
    elseif sc <= -5
        assess = "*   Moderate";
    else
        assess = "    Weak / No binding";
    end
    logInfo("%-3d  %-20s  %6.2f  %s", i, scored.Compound(i), sc, assess);
end
%[text] Well done — you have completed the docking exercise.
%[text] As a next step, consider exploring R04 (Protein-Ligand Analysis) or R09 (GNN Molecular Property Prediction).
%[text] ## Summary
%[text] - We built a complete molecular docking pipeline — pdbfixer for receptor preparation, RDKit for 3D ligand generation, meeko for PDBQT conversion, and AutoDock Vina for scoring — driven entirely from MATLAB via Python subprocess calls.
%[text] - Vina's empirical $\\Delta G$ score cannot capture covalent contributions; as a result, covalent inhibitors (Nirmatrelvir, Boceprevir) are systematically underscored relative to their measured potencies.
%[text] - Active-site geometry is conserved across the cysteine protease family, which explains why the HCV drug Boceprevir binds Mpro with measurable affinity — an outcome confirmed both computationally and experimentally.
%[text] - High docking scores for HIV drugs such as Lopinavir reflect shape coincidence, not Mpro selectivity. Clinical failure results from a combination of weak target affinity and pharmacokinetic issues (CYP3A4 autoinduction).
%[text] - Embedding a self-contained 3Dmol.js page in `uihtml` enables interactive ribbon visualization directly inside MATLAB Online, without requiring external tools.
%[text] - Docking is a hypothesis-generation tool. Experimental follow-up — IC50 measurement, X-ray crystallography, and in vivo pharmacology — remains essential before drawing conclusions about actual drug efficacy. \
%[text] Local helpers
%[text] Ternary Operator Helper (Local Function)
function out = ternary_(cond, trueVal, falseVal)
    if cond
        out = trueVal;
    else
        out = falseVal;
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]%   data: {"layout":"inline","rightPanelPercent":40}
%---
