%[text] # F03 Answers: Fingerprints
%[text] Reference answers for the exercise in f03_fingerprints.m.
addpath(genpath("src"));
emk.setup.initPython();
%%
%[text] ## Answer E1: Comparison of ON Bits for Morgan vs MACCS

smiles_list = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", ...
               "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", ...
               "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"};
mol_names   = {"Ethanol","Benzene","Aspirin","Caffeine","Ibuprofen"};

on_bits = zeros(5, 2);
for i = 1:5
    mol = emk.mol.fromSmiles(smiles_list{i});
    fp_morgan = emk.fingerprint.morgan(mol);
    fp_maccs  = emk.fingerprint.maccs(mol);
    on_bits(i, 1) = sum(emk.fingerprint.toArray(fp_morgan));
    on_bits(i, 2) = sum(emk.fingerprint.toArray(fp_maccs));
end

logInfo("%-14s  %8s  %6s", "Molecule Name", "Morgan", "MACCS");
logInfo("  %s", repmat("-", 1, 32));
for i = 1:5
    logInfo("%-14s  %8d  %6d", mol_names{i}, on_bits(i,1), on_bits(i,2));
end
logInfo("Average ON Bits -- Morgan: %.1f, MACCS: %.1f", ...
    mean(on_bits(:,1)), mean(on_bits(:,2)));
%[text] Morgan (2048 bits, Radius 2) typically generates 20 to 60 ON bits for drug-like molecules.
%[text] MACCS (167 bits) typically generates 30 to 60 ON bits.
%[text] MACCS has a higher density (ON/Total) due to its smaller size.
%%
%[text] ## Answer E2: Morgan FP Density vs Heavy Atom Count

data = readtable("data/list/everyday_chemicals.csv", TextType="string");
mols = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(data.SMILES), ...
    "UniformOutput", false);

n = numel(mols);
density    = zeros(n, 1);
heavy_count = zeros(n, 1);

for i = 1:n
    fp   = emk.fingerprint.morgan(mols{i});
    bits = emk.fingerprint.toArray(fp);
    density(i) = sum(bits) / numel(bits);
    d = emk.descriptor.calculate(mols{i}, ["HeavyAtomCount"]);
    heavy_count(i) = d.HeavyAtomCount;
end

figure("Name", "FP Density vs Heavy Atom Count");
scatter(heavy_count, density, 60, "filled");
xlabel("Heavy Atom Count");
ylabel("Morgan FP Density (ON Bits / 2048)");
title("Everyday Chemicals -- Fingerprint Density");
%[text] Larger molecules tend to have more ON bits, so density generally increases with HeavyAtomCount.
%[text] However, it saturates with a fixed size of 2048 bits, so the increase becomes gradual for very large molecules.
%%
%[text] ## Answer E3: Thought Question (No Code)
%[text] A larger Radius encodes a more extended circular environment,
%[text] generating more unique hash values. However, since all hashes are
%[text] folded into a fixed vector of 2048 bits, "collisions" can occur where different environments
%[text] map to the same bit, resulting in fewer ON bits than expected.
%[text] With a very large Radius (e.g., radius=6 for small molecules), most atoms
%[text] view the entire molecule, leading to fewer distinct environments, and the number of ON bits
%[text] may actually decrease compared to radius=2.
logInfo("E3 is a thought question -- refer to the comments in the answer file.");

%[text] ---
%[text] Try It Answers
%[text] ---
%%
%[text] ## Try It 1: Number of ON Bits for Aspirin at Radii (0, 1, 2, 3)

mol_asp = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");

radii   = [0, 1, 2, 3];
on_bits = zeros(1, numel(radii));
for i = 1:numel(radii)
    fp = emk.fingerprint.morgan(mol_asp, Radius=radii(i));
    on_bits(i) = sum(emk.fingerprint.toArray(fp));
end
logInfo("Aspirin Morgan ON Bits:");
for i = 1:numel(radii)
    logInfo("  Radius %d: %d ON Bits", radii(i), on_bits(i));
end
%[text] The number of ON bits increases significantly from Radius 0 to 1,
%[text] and the increase from Radius 2 to 3 is more gradual (due to hash collisions).
%%
%[text] ## Try It 2: Caffeine's Morgan FP; Get ON Bit Indices with find

fp   = emk.fingerprint.morgan(emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C"));
bits = emk.fingerprint.toArray(fp);
on_bits_idx = find(bits);
logInfo("Caffeine's ON Bits: %d", numel(on_bits_idx));
logInfo("ON Bit Indices (first 10): %s", mat2str(on_bits_idx(1:min(10,end))));
%[text] The indices are scattered across the entire 2048-bit space.
%[text] This is because the hash function is designed to "widely disperse adjacent environments."

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---