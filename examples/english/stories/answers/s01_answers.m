%[text] # S01 Answer: Find Caffeine Relatives
%[text] Reference answer for the "Try it yourself" exercise in s01_find_caffeine_cousins.m.
%[text] First, run s01_find_caffeine_cousins.m to build the database and fingerprints,
%[text] then use this file to check your answers.
addpath(genpath("src"));
emk.setup.initPython();
logInfo("S01 Answer: Setup Complete");
%%
%[text] ## Let's Try 1: Does Caffeine Meet the BBB Permeability Criteria?

mol_caf  = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
desc_caf = emk.descriptor.calculate(mol_caf, ["TPSA", "LogP"]);

logInfo("Caffeine TPSA = %.1f A^2  (BBB recommendation: < 90 A^2)", desc_caf.TPSA);
logInfo("Caffeine LogP = %.2f       (moderate lipophilicity)",  desc_caf.LogP);

%[text] Q: Does caffeine meet both conditions?
%[text] A: Yes. TPSA ~62 A^2 (< 90 -- good), LogP ~-1.0 (slightly negative --
%[text]    very high water solubility).
%[text]    Despite the low LogP, caffeine efficiently crosses the BBB due to its small size and
%[text]    planar aromatic structure.
%[text]    A negative LogP indicates high solubility in water (coffee, tea).
%%
%[text] ## Let's Try 2: Most Common Molecular Category & Heaviest Molecule?

tbl = readtable("data/list/everyday_chemicals.csv", "TextType", "string");

%[text] Most Common Molecular Category
cats        = unique(tbl.Category);
catCounts   = cellfun(@(c) sum(tbl.Category == c), cellstr(cats));
[~, iMax]   = max(catCounts);
logInfo("Most Common Molecular Category: %s (%d items)", cats(iMax), catCounts(iMax));

%[text] Heaviest Molecule
[maxMW, iHeavy] = max(tbl.MolecularWeight);
logInfo("Heaviest Molecule: %s  (MW = %.1f Da)", tbl.CommonName(iHeavy), maxMW);

%[text] A: Run and check the output. "Flavor" often becomes the most common.
%[text]    Sucrose (sugar) or lipid-like compounds tend to be the heaviest.
%%
%[text] ## Let's Try 3: How many ON bits are in the Morgan fingerprint of Caffeine?

fp_caf = emk.fingerprint.morgan(mol_caf);
bits   = emk.fingerprint.toArray(fp_caf);
nOn    = sum(bits);

logInfo("Caffeine ECFP4 ON bits: %d / %d  (Density %.1f%%)", ...
    nOn, numel(bits), 100 * nOn / numel(bits));

%[text] A: There are about 30-50 ON bits out of 2048.
%[text]    The exact number varies with the Morgan radius and bit count.
%[text]    Caffeine's bicyclic purine structure generates many different local environments,
%[text]    resulting in more ON bits compared to monocyclic or acyclic molecules.
%%
%[text] ## Let's Try 4: What is the top hit (excluding caffeine itself)?
%[text] The database fingerprints were calculated in s01 (use workspace variables).
%[text] Rebuild if running standalone:
tbl         = readtable("data/list/everyday_chemicals.csv", "TextType", "string");
nMols       = height(tbl);
all_fps     = cell(1, nMols);
valid       = false(1, nMols);
for i = 1:nMols
    if emk.mol.isValid(tbl.SMILES(i))
        all_fps{i} = emk.fingerprint.morgan(emk.mol.fromSmiles(tbl.SMILES(i)));
        valid(i)   = true;
    end
end
valid_fps   = all_fps(valid);
valid_names = tbl.CommonName(valid);

result = emk.similarity.rankBy(fp_caf, valid_fps);
logInfo("Rank 1 (self): %s  T=%.4f", valid_names(result.Indices(1)), result.Scores(1));
logInfo("Rank 2 (best non-self): %s  T=%.4f", ...
    valid_names(result.Indices(2)), result.Scores(2));

%[text] A: Rank 1 is caffeine itself (T = 1.0).
%[text]    Rank 2 is THEOBROMINE (found in cocoa and chocolate).
%[text]    Theobromine shares the same xanthine (purine) core,
%[text]    but has 2 N-methyl groups instead of 3. Tanimoto ~0.53.
%[text]    High similarity for a natural product database.
%%
%[text] ## Let's Try 5: Structural Differences Between Caffeine and Theobromine

mol_theo = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)NC(=O)N2C");
figure("Name", "Caffeine (N-CH3 x 3)", "Position", [100 500 440 380]);
emk.viz.draw2d(mol_caf,  Title="Caffeine (N-CH3 × 3)");
figure("Name", "Theobromine (N-CH3 x 2)", "Position", [560 500 440 380]);
emk.viz.draw2d(mol_theo, Title="Theobromine (N-CH3 × 2)");

logInfo("Caffeine  Molecular Formula: C8H10N4O2  (3 N-methyl groups: N1, N3, N7)");
logInfo("Theobromine Molecular Formula: C7H8N4O2  (2 N-methyl groups: N3, N7; N1 is NH)");

%[text] A: The only structural difference is the substituent at N1:
%[text] - Caffeine:     Methyl group at N1 (N-CH3)
%[text] - Theobromine:  Hydrogen at N1 (N-H)
%[text]    This single substitution difference results in the loss of one bit environment from the fingerprint,
%[text]    reducing the Tanimoto from 1.0 to ~0.53.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
