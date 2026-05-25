%[text] # F04 Answers: Similarity
%[text] Reference answers for the f04_similarity.m exercise.
addpath(genpath("src"));
emk.setup.initPython();
%%
%[text] ## Answer E1: Analgesics Similarity Matrix

smiles = {"CC(=O)Oc1ccccc1C(=O)O", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", ...
          "CC(c1ccc2cc(OC)ccc2c1)C(=O)O", "CC(=O)NC1=CC=C(C=C1)O", ...
          "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"};
names = {"Aspirin","Ibuprofen","Naproxen","Paracetamol","Caffeine"};

fps = cellfun(@(s) emk.fingerprint.morgan(emk.mol.fromSmiles(s)), smiles, ...
    "UniformOutput", false);
S = emk.similarity.matrix(fps);

figure("Name", "Analgesics Similarity", "Position", [100 100 500 440]);
imagesc(S);
colorbar;
clim([0, 1]);
colormap(parula);
xticks(1:5);  xticklabels(names);  xtickangle(30);
yticks(1:5);  yticklabels(names);
title("Tanimoto Similarity (Morgan ECFP4)");
for r = 1:5
    for c = 1:5
        text(c, r, sprintf("%.2f", S(r,c)), ...
            "HorizontalAlignment","center","FontSize",8,"Color","white");
    end
end

%[text] Find the most similar pair off the diagonal
S_masked = S - eye(5);
[max_val, lin_idx] = max(S_masked(:));
[r_idx, c_idx] = ind2sub([5 5], lin_idx);
logInfo("Most similar pair: %s -- %s  (Tanimoto = %.4f)", ...
    names{r_idx}, names{c_idx}, max_val);
%[text] Expected: Both Ibuprofen and Naproxen are propionic acid NSAIDs
%[text] and share a carboxylic acid + aromatic scaffold, resulting in the highest score.
%%
%[text] ## Answer E2: Top 3 most similar to Ibuprofen in everyday_chemicals

query_fp = emk.fingerprint.morgan(emk.mol.fromSmiles( ...
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"));
data    = readtable("data/list/everyday_chemicals.csv", TextType="string");
db_mols = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(data.SMILES), ...
    "UniformOutput", false);
db_fps  = cellfun(@(m) emk.fingerprint.morgan(m), db_mols, "UniformOutput", false);

result = emk.similarity.rankBy(query_fp, db_fps, 3);
logInfo("Top 3 most similar to Ibuprofen:");
for k = 1:3
    logInfo("  %d. %-20s  Tanimoto=%.4f", k, data.CommonName(result.Indices(k)), ...
        result.Scores(k));
end
%%
%[text] ## Answer E3: Count pairs exceeding Tanimoto threshold

threshold = 0.3;
count = 0;
for r = 1:5
    for c = (r+1):5   % Upper triangle only
        if S(r, c) >= threshold
            count = count + 1;
            logInfo("  %s -- %s  (%.4f)", names{r}, names{c}, S(r,c));
        end
    end
end
logInfo("Pairs with Tanimoto >= %.1f: %d", threshold, count);
%[text] Propionic acid NSAIDs (Ibuprofen, Naproxen) tend to exceed 0.3
%[text] due to sharing a carboxylic acid + aromatic ring scaffold.
%[text]
%[text] ---
%[text] Try It Answers
%[text] ---
%%
%[text] ## Try It 1: Tanimoto similarity between Ethanol and Methanol ("CO")

mol_eth  = emk.mol.fromSmiles("CCO");
mol_meth = emk.mol.fromSmiles("CO");
mol_prop = emk.mol.fromSmiles("CCCO");

fp_eth  = emk.fingerprint.morgan(mol_eth);
fp_meth = emk.fingerprint.morgan(mol_meth);
fp_prop = emk.fingerprint.morgan(mol_prop);

t_eth_meth = emk.similarity.tanimoto(fp_eth, fp_meth);
t_eth_prop = emk.similarity.tanimoto(fp_eth, fp_prop);

logInfo("Ethanol vs Methanol  = %.4f", t_eth_meth);
logInfo("Ethanol vs Propanol = %.4f", t_eth_prop);
%[text] Methanol is smaller than Ethanol and shares very little of the Morgan environment.
%[text] Ethanol vs Propanol is likely to score high due to similar chain extension.
%%
%[text] ## Try It 2: Similarity matrix using Dice coefficient

smiles_6 = ["CCO", "CCCO", "CC(C)O", "CC(=O)Oc1ccccc1C(=O)O", ...
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CC(=O)NC1=CC=C(C=C1)O"];
mols_6 = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(smiles_6), "UniformOutput", false);
fps_6  = cellfun(@(m) emk.fingerprint.morgan(m), mols_6, "UniformOutput", false);

S      = emk.similarity.matrix(fps_6);
S_dice = emk.similarity.matrix(fps_6, Metric="dice");

logInfo("All elements with Dice >= Tanimoto: %d", all(S_dice(:) >= S(:)));

S_off      = S      - eye(numel(fps_6));
S_dice_off = S_dice - eye(numel(fps_6));
[~, idx_t] = max(S_off(:));
[~, idx_d] = max(S_dice_off(:));
logInfo("Most similar pair -- Tanimoto idx: %d, Dice idx: %d  (same: %d)", ...
    idx_t, idx_d, isequal(idx_t, idx_d));

%[text] **Most dissimilar pair** (lowest off-diagonal Tanimoto score):
[~, idx_min] = min(S_off(S_off > 0));   % exclude zeros from diagonal
% Reconstruct the true off-diagonal minimum
S_off_pos = S_off;
S_off_pos(S_off_pos <= 0) = Inf;
[min_val, lin_idx_min] = min(S_off_pos(:));
[r_min, c_min] = ind2sub([numel(fps_6), numel(fps_6)], lin_idx_min);
names_6 = {"Ethanol","Propanol","Isopropanol","Aspirin","Caffeine","Paracetamol"};
logInfo("Most dissimilar pair: %s -- %s  (Tanimoto = %.4f)", ...
    names_6{r_min}, names_6{c_min}, min_val);
%[text] In binary bit vectors, Dice >= Tanimoto always holds (Dice = 2J/(A+B), Tanimoto = J/(A+B-J)).
%[text] The ranking of molecular pairs usually does not change even with different metrics.
%%
%[text] ## Try It 3: Search using Ibuprofen as a query

query_mol = emk.mol.fromSmiles("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");   % Ibuprofen
query_fp  = emk.fingerprint.morgan(query_mol);

data    = readtable("data/list/everyday_chemicals.csv", TextType="string");
valid_mask_db = cellfun(@(s) emk.mol.isValid(s), cellstr(data.SMILES));
db_smiles = cellstr(data.SMILES(valid_mask_db));
db_mols   = cellfun(@(s) emk.mol.fromSmiles(s), db_smiles, "UniformOutput", false);
db_fps    = cellfun(@(m) emk.fingerprint.morgan(m), db_mols, "UniformOutput", false);
db_names  = data.CommonName(valid_mask_db);

result = emk.similarity.rankBy(query_fp, db_fps, 5);
logInfo("Top 5 everyday chemicals most similar to Ibuprofen:");
for k = 1:5
    idx   = result.Indices(k);
    label = db_names(idx);
    if strlength(label) == 0; label = "(unknown)"; end
    logInfo("  %d. %-20s  Tanimoto=%.4f", k, label, result.Scores(k));
end
%[text] Compounds with aromatic carboxylic acid scaffolds (such as benzoic acid derivatives)
%[text] tend to rank high.

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---