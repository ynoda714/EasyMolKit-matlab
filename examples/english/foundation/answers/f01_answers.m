%[text] # F01 Answers: Drawing Molecules from SMILES
%[text] Reference solution for the exercise in `f01_draw_molecules.m`.
addpath(genpath("src"));
emk.setup.initPython();
%%
%[text] ## Answer E1: Let's draw Nicotine
mol_nic = emk.mol.fromSmiles("CN1CCC[C@H]1C2=CN=CC=C2");
figure("Name", "Nicotine");
emk.viz.draw2d(mol_nic, Title="Nicotine");
%[text] RDKit draws two rings (pyrrolidine and pyridine).
%[text] The chiral center (@H) is encoded in SMILES, but in 2D drawing,
%[text] it may not be explicitly shown as a wedge bond.
%%
%[text] ## Answer E2: Let's arrange and compare Terpenes

smiles_pair = ["OC1CC(CC(C1)C(C)C)C", "Cc1ccc(cc1O)C(C)C"];
names_pair  = ["Menthol", "Thymol"];

figure("Name", "Menthol vs Thymol", "Position", [100 100 800 400]);
for i = 1:2
    subplot(1, 2, i);
    mol = emk.mol.fromSmiles(smiles_pair(i));
    emk.viz.draw2d(mol, Title=names_pair(i));
    can = emk.mol.toSmiles(mol);
    logInfo("%s's canonical SMILES: %s", names_pair(i), can);
end
%[text] Menthol is a saturated cyclic alcohol (C10H20O).
%[text] Thymol is a phenol (aromatic, C10H14O).
%[text] The molecular formulas differ, so the canonical SMILES are also different.
%%
%[text] ## Answer E3: Let's create a validity filter

candidates = {"CCO", "C(C", "c1ccccc1", "BADSMILES", "CC(=O)O"};
valid_smiles = {};
for i = 1:numel(candidates)
    if emk.mol.isValid(candidates{i})
        valid_smiles{end+1} = candidates{i}; %#ok<SAGROW>
    end
end
logInfo("Valid SMILES (%d / %d):", numel(valid_smiles), numel(candidates));
for i = 1:numel(valid_smiles)
    logInfo("  %s", valid_smiles{i});
end
%[text] Valid: "CCO", "c1ccccc1", "CC(=O)O" (3 out of 5)
%[text] "C(C" has unclosed parentheses and is a syntax error. "BADSMILES" is not a SMILES string.

%[text] ---
%[text] Try It Answers
%[text] ---
%%
%[text] ## Try It 1: Class Check for Ibuprofen and Nicotine

mol_ibupro = emk.mol.fromSmiles("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
mol_nicot  = emk.mol.fromSmiles("CN1CCC[C@H]1C2=CN=CC=C2");
logInfo("Ibuprofen: class = %s", class(mol_ibupro));
logInfo("Nicotine:       class = %s", class(mol_nicot));
%[text] Both return `py.rdkit.Chem.rdchem.Mol`.
%[text] MATLAB only holds references to Python objects,
%[text] so you can use them as inputs to `emk.*` functions without being aware of the RDKit API.
%%
%[text] ## Try It 2: Canonical SMILES for Pyridine and Furan

mol_pyr  = emk.mol.fromSmiles("c1ccncc1");   % Pyridine
mol_fur  = emk.mol.fromSmiles("c1ccoc1");    % Furan

can_pyr = emk.mol.toSmiles(mol_pyr);
can_fur = emk.mol.toSmiles(mol_fur);
logInfo("Canonical SMILES for Pyridine: %s", can_pyr);
logInfo("Canonical SMILES for Furan:   %s", can_fur);
%[text] RDKit canonical forms: Pyridine → `c1ccncc1` (aromatic, retains lowercase n)
%[text] Furan → `c1ccoc1` (aromatic, retains lowercase o)
%[text] RDKit uses Hückel's rule (4n+2 electrons) to determine both as aromatic.
%%
%[text] ## Try It 3: Why "C(C)(C)(C)(C)C" is Invalid

tf_neo  = emk.mol.isValid("C(C)(C)(C)(C)C");
logInfo("C(C)(C)(C)(C)C isValid = %d", tf_neo);
%[text] isValid = 0 (false). While you can write a 5-coordinate carbon in SMILES syntax,
%[text] RDKit's isValid checks for chemical validity (valence rules).
%[text] Since the maximum number of bonds for carbon is 4, RDKit deems this molecule invalid.

tf_bad  = emk.mol.isValid("C(CC");
logInfo("C(CC            isValid = %d", tf_bad);
%[text] isValid = 0. Unclosed parentheses are a syntax error in SMILES.
%%
%[text] ## Try It 4: 5 Molecule Grid (Adaptive Subplot Layout)

smiles_list = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", ...
               "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", ...
               "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"];   % Added Ibuprofen
names       = ["Ethanol", "Benzene", "Aspirin", "Caffeine", "Ibuprofen"];

N    = numel(smiles_list);
ncol = ceil(sqrt(N));
nrow = ceil(N / ncol);      % For N=5, 2 rows x 3 columns

figure("Name", "Molecule Grid (5)", "Position", [100 100 1000 700]);
for i = 1:N
    subplot(nrow, ncol, i);
    mol = emk.mol.fromSmiles(smiles_list(i));
    emk.viz.draw2d(mol, Title=names(i));
end
%[text] ceil(sqrt(5)) = 3 columns, ceil(5/3) = 2 rows → 2×3 grid (1 space empty)

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---