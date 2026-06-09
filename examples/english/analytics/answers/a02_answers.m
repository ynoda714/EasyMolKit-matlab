%[text] # A02 Answers: Molecular Clustering
%[text] Reference answers for the "Try It Yourself" exercise in a02_molecular_clustering.m.
%[text] First, run a02_molecular_clustering.m to make workspace variables
%[text] (validFps, molNames, molCats, simMat, distMat, Z_tree, bitMat_double,
%[text] available for use.
addpath(genpath("src"));
emk.setup.initPython();
logInfo("A02 Answers: Setup Complete");

%[text] ── Reconstructing Workspace Variables (Standalone Execution) ──
DATA_FILE = "data/list/everyday_chemicals.csv";
FP_RADIUS = 2;
FP_NBITS  = 2048;
rawTbl    = readtable(DATA_FILE, TextType="string");
fps       = cell(1, height(rawTbl));
valid     = false(1, height(rawTbl));
for k = 1:height(rawTbl)
    try
        mol = emk.mol.fromSmiles(rawTbl.SMILES(k));
        fps{k} = emk.fingerprint.morgan(mol, Radius=FP_RADIUS, NBits=FP_NBITS);
        valid(k) = true;
    catch
    end
end
validIdx  = find(valid);
validFps  = fps(validIdx);
molNames  = rawTbl.CommonName(validIdx);
molCats   = rawTbl.Category(validIdx);
nMols     = numel(validIdx);

simMat  = emk.similarity.matrix(validFps, Metric="tanimoto");
distMat = 1 - simMat;
distMat(logical(eye(nMols))) = 0;

bitMat = zeros(nMols, FP_NBITS, "logical");
for j = 1:nMols
    bitMat(j, :) = emk.fingerprint.toArray(validFps{j});
end
bitMat_double = double(bitMat);

condDist = squareform(distMat, "tovector");
Z_tree   = linkage(condDist, "ward");
%%
%[text] ## Let's Try 1: Caffeine On-bits; Does increasing on-bits improve Tanimoto?

fp_caf = validFps{molNames == "caffeine"};
bits   = emk.fingerprint.toArray(fp_caf);
nOn    = sum(bits);
logInfo("Caffeine ECFP4: %d / %d bits set (density %.1f%%)", ...
    nOn, FP_NBITS, 100 * nOn / FP_NBITS);

%[text] Compare the on-bit density of all molecules
densities = sum(bitMat_double, 2) / FP_NBITS;
logInfo("On-bit density: Min=%.1f%%  Mean=%.1f%%  Max=%.1f%%", ...
    min(densities)*100, mean(densities)*100, max(densities)*100);

%[text] **Answer:** Caffeine has ~40–50 bits on out of 2048 bits (density ~2%).
%[text] The bicyclic purine skeleton generates more local environments than simple acyclic molecules.
%[text] **In sparse vectors, more on-bits result in lower Tanimoto**: T = c/(a+b−c).
%[text] Even if molecule A has 200 bits and B has 10 bits, with all 10 bits of B overlapping with A,
%[text] T = 10/(200+10−10) = 10/200 = 0.05, which is very low.
%[text] Molecules with many on-bits show low Tanimoto scores with molecules that do not also have many on-bits.
%%
%[text] ## Let's Try 2: Average Pairwise Similarity; Most Similar Pair

offDiag = simMat(~logical(eye(nMols)));
logInfo("Average Pairwise Tanimoto: %.3f", mean(offDiag));
logInfo("Median Pairwise Tanimoto: %.3f", median(offDiag));

%[text] Search for the most similar pair
C = simMat;
C(logical(eye(nMols))) = 0;
[maxVal, linIdx] = max(C(:));
[r, c] = ind2sub(size(C), linIdx);
logInfo("Most Similar Pair: %s & %s  (T = %.3f)", ...
    molNames(r), molNames(c), maxVal);

%[text] Display similarity matrix sorted by category
[sortedCats, sortOrd] = sort(molCats);
simSorted = simMat(sortOrd, sortOrd);
figure("Name", "A02 Sorted by Category");
imagesc(simSorted); colormap(parula); colorbar; caxis([0 1]);
title("Similarity Matrix Sorted by Category");
xticks(1:nMols); xticklabels(cellstr(molNames(sortOrd))); xtickangle(90);
yticks(1:nMols); yticklabels(cellstr(molNames(sortOrd)));
set(gca, FontSize=6);

%[text] **Answer:** In a diverse set of 30 types of household chemicals, the average pairwise Tanimoto is
%[text] usually 0.05–0.12, indicating structural diversity.
%[text] The most similar pair is often Limonene & Carvone (both cyclic terpenes) or
%[text] Caffeine & Theobromine (both methylxanthines).
%[text] When sorted by category, structurally consistent categories (sugars, stimulants) show
%[text] bright blocks along the diagonal.
%%
%[text] ## Let's Try 3: Effect of Cut Height; Linkage Method Comparison

c1 = cluster(Z_tree, Cutoff=0.2*max(Z_tree(:,3)), Criterion="distance");
c2 = cluster(Z_tree, Cutoff=0.6*max(Z_tree(:,3)), Criterion="distance");
logInfo("Cut at 20%% of maximum height -> %d clusters", max(c1));
logInfo("Cut at 60%% of maximum height -> %d clusters", max(c2));

%[text] Comparison with Average and Complete Linkage
Z_avg  = linkage(condDist, "average");
Z_comp = linkage(condDist, "complete");

figure("Name", "A02 Linkage Comparison");
tiledlayout(1, 3);
nexttile; dendrogram(Z_tree, 0, Orientation="left"); title("Ward");    xlabel("Distance"); set(gca, FontSize=6);
nexttile; dendrogram(Z_avg,  0, Orientation="left"); title("Average"); xlabel("Distance"); set(gca, FontSize=6);
nexttile; dendrogram(Z_comp, 0, Orientation="left"); title("Complete");xlabel("Distance"); set(gca, FontSize=6);
sgtitle("Linkage Method Comparison -- Household Chemicals");

%[text] Similarity between Glucose and Fructose
gIdx = find(molNames == "glucose");
fIdx = find(molNames == "fructose");
if ~isempty(gIdx) && ~isempty(fIdx)
    logInfo("Glucose-Fructose Tanimoto: %.3f", simMat(gIdx, fIdx));
end

%[text] **Answer:** Lowering the cut to 20% of the maximum height results in many small clusters (including singletons).
%[text] Raising it to 60% consolidates most molecules into 3–5 large clusters.
%[text] The average linkage is slightly less compact than Ward but also less sensitive to outliers.
%[text] Complete linkage can yield very different results as it forces the maximum diameter clusters.
%[text] The Tanimoto between glucose and fructose (structural isomers, C6H12O6) is moderate (~0.3–0.4).
%[text] This is because ECFP4, by default, does not consider stereochemistry, leading cyclic and acyclic tautomers to
%[text] share many substructures.
%%
%[text] ## Let's Try 4: Elbow Method — At which k does diminishing returns start?

K_MAX = 8;
wcss  = zeros(1, K_MAX);
rng(42);
for k = 2:K_MAX
    [~, ~, sumd] = kmeans(bitMat_double, k, ...
        Distance="hamming", Replicates=5, Display="off");
    wcss(k) = sum(sumd);
end

%[text] Find the maximum drop point (elbow) using differences
drops = diff(wcss(2:K_MAX));
[~, elbowRel] = min(drops);   % smallest delta = where drop plateaus
elbowK = elbowRel + 2;        % offset: index 1 -> k=3
logInfo("Approximate elbow k = %d  (WCSS drop plateau)", elbowK);

%[text] **Answer:** For this set of 30 molecules, the elbow is usually at k = 4–6.
%[text] Hamming distance is a natural metric for binary bit vectors,
%[text] counting the proportion of differing positions (a+b−2c)/(a+b).
%[text] This corresponds to the complement of the Dice similarity.
%%
%[text] ## Let's Try 5: Comparison of k-means and Hierarchical Clustering; Hierarchical Silhouette

rng(42);
[bestSilScore, bestK_idx] = deal(0, 0);
silScores = zeros(1, K_MAX);
for k = 2:K_MAX
    labels_k = kmeans(bitMat_double, k, Distance="hamming", Replicates=5, Display="off");
    sil = silhouette(bitMat_double, labels_k, "hamming");
    silScores(k) = mean(sil);
    if silScores(k) > bestSilScore
        bestSilScore = silScores(k);
        bestK_idx    = k;
    end
end
bestK = bestK_idx;
rng(42);
finalLabels_km = kmeans(bitMat_double, bestK, Distance="hamming", Replicates=10, Display="off");

%[text] Hierarchical Clustering at bestK
hier_labels = cluster(Z_tree, MaxClust=bestK);
sil_hier    = silhouette(bitMat_double, hier_labels, "hamming");
sil_km      = silhouette(bitMat_double, finalLabels_km, "hamming");

logInfo("K-Means Silhouette (k=%d): %.3f",       bestK, mean(sil_km));
logInfo("Hierarchical Silhouette (k=%d): %.3f",  bestK, mean(sil_hier));
logInfo("Winner: %s", ...
    string(ternary_(mean(sil_km) >= mean(sil_hier), "K-Means", "Hierarchical")));

%[text] Checking Category Recall
[~, catIdx] = ismember(molCats, unique(molCats));
logInfo("Category Labels vs k-means (k=%d): Comparable with crosstab()", bestK);

%[text] **Answer:** Both methods usually agree at the level of xanthines (Caffeine, Theobromine),
%[text] terpenes (Limonene, Carvone), and sugars (Sucrose, Glucose, Fructose).
%[text] The silhouette scores for this small dataset are modest (~0.2–0.4).
%[text] Very tight clusters are needed for perfect scores.
%[text] K-Means and Hierarchical Clustering may yield different results for borderline molecules.
%%
%[text] ## Let's Try 6: Cluster Representative Selection Strategy

CUT_HEIGHT     = 0.4 * max(Z_tree(:, 3));
clusterLabels_hier = cluster(Z_tree, Cutoff=CUT_HEIGHT, Criterion="distance");
nClusters_hier = max(clusterLabels_hier);

logInfo("Cluster Representative (Highest Average Intra-cluster Similarity):");
for c = 1:nClusters_hier
    idx = find(clusterLabels_hier == c);
    if numel(idx) == 1
        rep = idx;
    else
        subSim = simMat(idx, idx);
        meanSim = (sum(subSim, 2) - 1) / max(numel(idx) - 1, 1);
        [~, repLocal] = max(meanSim);
        rep = idx(repLocal);
    end
    logInfo("  Cluster H%d: Representative = %s", c, molNames(rep));
end

%[text] **Answer:** A representative molecule is the one with the highest average Tanimoto similarity to all molecules in the same cluster.
%[text] This can be called the "centroid" in the fingerprint space.
%[text] In actual screening campaigns, one representative molecule per cluster is selected to evaluate bioactivity,
%[text] and exploration around active clusters is expanded.
logInfo("A02 Answer Complete.");

%[text] Local Helper (Alternative to Anonymous Functions with Branching)
function out = ternary_(cond, a, b)
    if cond; out = a; else; out = b; end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]%   data: {"layout":"inline","rightPanelPercent":40}
%---
