%[text] # A02: Molecular Clustering — Grouping Chemicals by Structural Similarity
%[text] EasyMolKit Analytics — Layer 3
%[text] 
%[text] There are 30 compounds isolated by natural product chemists from food and household items.
%[text] By identifying "structurally similar groups" before bioactivity testing, redundant tests can be avoided, reducing costs.
%[text] In this script, you will learn how to group using hierarchical clustering and k-means with Morgan ECFP4 fingerprints and Tanimoto distance.
%[text] 
%[text] **What you will learn in this tutorial**
%[text] - Understand why fingerprints are used instead of descriptors for clustering
%[text] - Learn how to construct a Tanimoto matrix with `emk.similarity.matrix()`
%[text] - Understand how to convert Tanimoto similarity into a distance metric for hierarchical methods
%[text] - Learn how to use `linkage()` and `dendrogram()` from the Statistics Toolbox
%[text] - Learn how to evaluate cluster quality using MATLAB's `kmeans()` and `silhouette()`
%[text] - Learn how to interpret clustering results in a chemical context \
%[text] 
%[text] **Prerequisites**
%[text] - Completion of F03 (Fingerprints) and F04 (Similarity)
%[text] - Recommended: Understanding context from A01 (Chemical Space PCA)
%[text] - Use of Statistics and Machine Learning Toolbox (linkage, kmeans, silhouette)
%[text] - No internet connection required \
%[text] 
%[text] Estimated time required: 30–45 minutes
%[text] 
%[text] **Data:**
%[text] `data/list/everyday_chemicals.csv` — 30 common molecules (PubChem CC0)
%[text] 
%[text] **References**
%[text] - Willett P, Barnard JM, Downs GM (1998) Chemical similarity searching. J Chem Inf Comput Sci 38:983-996. doi:10.1021/ci9800211
%[text] - Rogers D & Hahn M (2010) Extended-connectivity fingerprints. J Chem Inf Model 50:742-754. doi:10.1021/ci100050t
%[text] - Ward JH (1963) Hierarchical grouping to optimise an objective function. J Am Stat Assoc 58:236-244. doi:10.1080/01621459.1963.10500845
%[text] - Kaufman L & Rousseeuw PJ (1990) Finding Groups in Data: An Introduction to Cluster Analysis. Wiley. (silhouette coefficient, Chapter 2) \
%[text] 
%[text] How to run: Execute each section with Ctrl+Enter
%%
%[text] ## Section 0: Environment Setup
logSection("A02", "Section 0: Setup", "Analytics L3");
% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
emk.setup.initPython();

%[text] Before the main process, warm up the Python/RDKit process
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;
%%
%[text] ## Section 1: Morgan Fingerprint Calculation for All Molecules
%[text] 
%[text] ### Concept: Why Use Fingerprints for Structural Clustering
%[text] Descriptor-based methods (PCA, regression) capture physicochemical properties (size, polarity, etc.). However, they may miss structural similarity when two molecules have similar properties but different scaffolds, or vice versa.
%[text] 
%[text] Fingerprints encode the presence of local substructures (circular environment around each atom) as binary bit vectors. Two molecules with the same scaffold will share many common bits even if their substituents differ. This is a property not captured by descriptors.
%[text] 
%[text] **Morgan / ECFP Fingerprint:**
%[text] - Radius 2 (ECFP4): Captures substructures up to two bonds from each atom (the diameter 4 is the origin of the name ECFP**4**)
%[text] - 2048 bits (NBits=2048): Minimizes collision probability
%[text] - Independent of molecular orientation (unlike 3D fingerprints) \
%[text] 
%[text] **When Descriptor-Based Clustering is Suitable:** When physical properties like solubility or absorption define similarity (formulation/ADMET purposes)
%[text] **When Fingerprint-Based Clustering is Suitable:** When scaffold/substructure similarity is important (SAR, medicinal chemistry)
%[text] 
logSection("A02", "Section 1: Morgan Fingerprint Calculation for All Molecules", "Analytics L3");
DATA_FILE = "data/list/everyday_chemicals.csv";
FP_RADIUS = 2;     % Morgan radius 2 = ECFP4
FP_NBITS  = 2048;

rawTbl = readtable(DATA_FILE, TextType="string");
nRaw   = height(rawTbl);
logInfo("Loaded %d molecules from %s", nRaw, DATA_FILE);
%[text] Parse SMILES and calculate fingerprints
fps   = cell(1, nRaw);
valid = false(1, nRaw);
for k = 1:nRaw
    try
        mol     = emk.mol.fromSmiles(rawTbl.SMILES(k));
        fps{k}  = emk.fingerprint.morgan(mol, Radius=FP_RADIUS, NBits=FP_NBITS);
        valid(k) = true;
    catch ME
        logWarn("Skipping %s: %s", rawTbl.CommonName(k), ME.message);
    end
end

validIdx  = find(valid);
validFps  = fps(validIdx);
molNames  = rawTbl.CommonName(validIdx);
molCats   = rawTbl.Category(validIdx);
nMols     = numel(validIdx);
logInfo("ECFP4 fingerprint calculation complete: %d / %d molecules (%d bits)", ...
    nMols, nRaw, FP_NBITS);
%[text] **💡 Observation Point 1 — Investigate Caffeine's Fingerprint Density**
%[text] Let's check how many on-bits (set bits) caffeine has.
%[text] Refer to the following code:
%[text] fp\_caf = validFps{molNames == "caffeine"};
%[text] bits = `emk.fingerprint.toArray(fp_caf)`;
%[text] sum(bits)   % Number of on-bits
%[text] sum(bits) / FP\_NBITS  % Density
%[text] Consider whether molecules with many on-bits always have a high Tanimoto coefficient with other molecules.
%[text] Hint: Tanimoto = |A and B| / |A or B|
%[text] Consider whether having many on-bits is advantageous or disadvantageous for Tanimoto.
%[text] 
% ... (Try writing code here)
%%
%[text] ## Section 2: Construction of Tanimoto Similarity Matrix
%[text] In Section 1, we calculated the ECFP4 fingerprints for all molecules.
%[text] In this section, we will calculate the similarity for all N×N molecular pairs and represent it as a matrix.
%[text] 
%[text] ### Concept: Tanimoto Similarity as an Alternative Measure of Molecular Distance
%[text] Tanimoto (Jaccard) coefficient for binary fingerprints: 
%[text]{"align":"center"} $T(A,B) = |A \\cap B| / |A \\cup B| = c / (a + b - c)$
%[text]{"align":"center"} (where $a = |A|$, $b = |B|$, $c$ = number of common on-bits)
%[text] 
%[text] **Properties:** $T = 1.0$ indicates identical fingerprints, $T = 0.0$ indicates no common on-bits.
%[text] $d = 1 - T$ (Jaccard distance) is a formal distance measure satisfying the triangle inequality.
%[text] 
%[text] **Industry Standard Tanimoto Threshold Guidelines:** $T \\geq 0.85$ is nearly identical, $T \\geq 0.70$ is a close analog, $T \\geq 0.40$ is the same chemical class, $T \< 0.40$ is structurally diverse.
%[text] 
logSection("A02", "Section 2: Construction of Tanimoto Similarity Matrix", "Analytics L3");
simMat  = emk.similarity.matrix(validFps, Metric="tanimoto");  % N x N double
distMat = 1 - simMat;                                          % Distance matrix
%[text] Set diagonal elements to exactly 0 (to counteract numerical noise in Python IPC).
distMat(logical(eye(nMols))) = 0;

logInfo("Similarity matrix: %d x %d (range: %.3f - %.3f)", ...
    nMols, nMols, min(simMat(:)), max(simMat(~logical(eye(nMols)))));
%[text] \-- Visualize the similarity matrix as a heatmap --
figure("Name", "A02 Tanimoto Similarity Matrix");
imagesc(simMat);
colormap(parula);
colorbar;
caxis([0 1]);
title("All-to-All Tanimoto Similarity (ECFP4)");
xlabel("Molecule Index");
ylabel("Molecule Index");
%[text] Overlay category dividing lines and sort by category for clarity.
[sortedCats, sortOrd] = sort(molCats);
logInfo("Displayed %d x %d similarity matrix heatmap", nMols, nMols);
%[text] **💡 Observation Point 2 — Reading the Similarity Matrix**
%[text] Check the average Tanimoto similarity for all pairs.
%[text] Hint: offDiag = simMat(~logical(eye(nMols))); mean(offDiag)
%[text] Find the most similar molecular pair (non-identical with highest Tanimoto).
%[text] Hint: \[r,c\] = find(simMat == max(simMat(~logical(eye(nMols)))));
%[text] molNames(r(1)), molNames(c(1))
%[text] Consider if that pair is chemically meaningful.
%[text] Replot the matrix sorted by category and observe (use sortOrd).
%[text] Check if molecules in the same category cluster into bright blocks.
% ... (Try writing code here)
%%
%[text] ## Section 3: Hierarchical Clustering and Dendrograms
%[text] The similarity matrix creation is complete.
%[text] Next, we will use 1 - Tanimoto coefficient as the distance to group molecules in a bottom-up manner using hierarchical clustering.
%[text] 
%[text] ### Concept: Hierarchical Clustering (Ward's Method)
%[text] Hierarchical clustering builds a tree (dendrogram) in a bottom-up approach.
%[text] 1. Start: Each molecule forms its own cluster.
%[text] 2. Merge the two closest clusters (determined by linkage criteria).
%[text] 3. Repeat this process until all molecules form a single cluster. \
%[text] 
%[text] **Linkage Criteria:** `"single"` — Distance to nearest neighbor (chain effect) / `"complete"` — Distance to farthest neighbor / `"average"` — Average of all pairwise distances (UPGMA) / `"ward"` — Minimize variance within clusters (Ward 1963)
%[text] 
%[text] Ward's method is preferred in cheminformatics for generating compact and evenly sized clusters. It shares the objective function of minimizing within-cluster sum of squares with k-means, often yielding similar results.
%[text] 
%[text] **Reading a Dendrogram:** The height of the linkage indicates the dissimilarity when two clusters are merged. A horizontal line at height $h$ ("cut") defines $k$ clusters. Groups merged at lower heights are very similar.
%[text] 
logSection("A02", "Section 3: Hierarchical Clustering and Dendrograms", "Analytics L3");
%[text] Convert the symmetric distance matrix to a condensed vector (upper triangle) required by MATLAB's linkage() function.
condDist = squareform(distMat, "tovector");

Z_tree = linkage(condDist, "ward");
%[text] \-- Plotting the Dendrogram --
figure("Name", "A02 Hierarchical Clustering Dendrogram");
D = dendrogram(Z_tree, 0, ...
    Labels=cellstr(molNames), ...
    Orientation="left", ...
    ColorThreshold=0.65 * max(Z_tree(:, 3)));  % 65% of maximum height
title("Ward Method Dendrogram — Household Chemicals (ECFP4)");
xlabel("Ward Distance");
set(gca, FontSize=7);
logInfo("Displayed the dendrogram");
%[text] \-- Automatic cut at 65% of maximum height (approximately 7-8 clusters) --
CUT_HEIGHT   = 0.65 * max(Z_tree(:, 3));
clusterLabels_hier = cluster(Z_tree, Cutoff=CUT_HEIGHT, Criterion="distance");
nClusters_hier = max(clusterLabels_hier);
logInfo("Hierarchical cut: threshold=%.3f -> %d clusters", CUT_HEIGHT, nClusters_hier);

for c = 1:nClusters_hier
    members = molNames(clusterLabels_hier == c);
    logInfo("  Cluster H%d (%d molecules): %s", c, numel(members), ...
        strjoin(members, ", "));
end
%[text] **💡 Observation Point 3 — Try Changing Cut Height and Linkage Method**
%[text] Lower the cut height to 20% of maximum or raise it to 60%,
%[text] and observe how the number of clusters changes.
%[text] Refer to the following code:
%[text] c1 = cluster(Z\_tree, Cutoff=0.2\*max(Z\_tree(:,3)), Criterion="distance");
%[text] c2 = cluster(Z\_tree, Cutoff=0.6\*max(Z\_tree(:,3)), Criterion="distance");
%[text] \[max(c1), max(c2)\]
%[text] Try changing the linkage method to "average" or "complete".
%[text] Refer to the following code: Z2 = linkage(condDist, "complete"); figure; dendrogram(Z2,0);
%[text] Observe how the shape of the tree differs from Ward's method.
%[text] Consider which method generates more balanced clusters for this dataset.
%[text] Both glucose and fructose are C6H12O6 (structural isomers).
%[text] Check if they are merged at an early stage (low height) in the dendrogram.
%[text] Consider what this suggests about ECFP4 features for structural isomers.
% ... (Try writing code here)
%%
%[text] ## Section 4: K-means Clustering on Fingerprint Bit Vectors
%[text] Hierarchical clustering builds the entire tree at once and then determines the cut position.
%[text] In contrast, k-means fixes k and directly partitions into k clusters.
%[text] 
%[text] ### Concept: K-means with Hamming Distance
%[text] K-means repeatedly performs the following to partition $N$ points into $k$ clusters:
%[text] 1. Assign each point to the nearest centroid.
%[text] 2. Recalculate the centroid as the mean of the cluster.
%[text] 3. Repeat until the assignment stabilizes. \
%[text] 
%[text] Standard k-means uses Euclidean distance, but Hamming distance (proportion of differing bits) is suitable for binary fingerprints. This is supported by the `kmeans()` option `Distance="hamming"`.
%[text] 
%[text] **Choosing $k$:** Use domain knowledge (expected number of chemical classes), elbow method (within-cluster sum of squares vs $k$), silhouette analysis (peak average silhouette width at optimal $k$).
%[text] 
%[text] Convert fingerprints to a numerical bit matrix ($N \\times$ FP\_NBITS).
logSection("A02", "Section 4: K-means Clustering on Fingerprint Bit Vectors", "Analytics L3");
bitMat = zeros(nMols, FP_NBITS, "logical");
for j = 1:nMols
    bitMat(j, :) = emk.fingerprint.toArray(validFps{j});
end
bitMat_double = double(bitMat);
logInfo("Fingerprint bit matrix: %d x %d", size(bitMat_double));
%[text] \-- Elbow Method: Evaluate for k = 2..8 --
K_MAX = 8;
wcss  = zeros(1, K_MAX);     % Within-cluster sum of squares

rng(42);   % For reproducibility
for k = 2:K_MAX
    [~, ~, sumd] = kmeans(bitMat_double, k, ...
        Distance="hamming", Replicates=5, Display="off");
    wcss(k) = sum(sumd);
end

figure("Name", "A02 K-Means Elbow Curve");
plot(2:K_MAX, wcss(2:K_MAX), "-o", Color=[0.2 0.5 0.8], LineWidth=1.5);
xlabel("Number of Clusters k");
ylabel("Within-cluster Sum of Squares");
title("Elbow Method — K-Means on ECFP4 Bit Vectors");
grid on;
logInfo("Displayed elbow curve (look for the point where the steep drop flattens)");
%[text] **💡 Observation Point 4 — Find Optimal k at the Elbow**
%[text] Check at which k the elbow (point of diminishing returns) appears on the WCSS curve. That is the optimal k for this dataset. Consider why Hamming distance is used instead of Euclidean distance here.
%[text] Hint: Each element of the fingerprint vector is 0 or 1.
%[text] The Euclidean distance between 2 bits is sqrt(0)=0 or sqrt(1)=1, same as Hamming distance. However, for binary vectors with varying bit density, Euclidean distance is hard to generalize.
% ... (Try writing code here)
%%
%[text] ## Section 5: Silhouette Analysis for Selecting k
%[text] In the k-means method, you need to specify the number of clusters k in advance.
%[text] Let's learn how to choose the optimal k based on data using silhouette analysis.
%[text] 
%[text] ### Concept: Silhouette Coefficient
%[text] For each molecule $i$, the silhouette coefficient $s(i)$ indicates how well the molecule fits within its own class. It is measured in comparison to the nearest class: $s(i) = (b(i) - a(i)) / \\max(a(i), b(i))$
%[text] (where $a(i)$ is the average distance within the class, and $b(i)$ is the average distance to the nearest class)
%[text] 
%[text] **Interpretation:** If $s \\approx 1$, the molecule is well-clustered. If $s \\approx 0$, the molecule is on the boundary between two clusters. If $s \\approx -1$, it may be incorrectly assigned to a cluster.
%[text] 
%[text] The average silhouette of all molecules is an overall indicator of cluster quality. The k that maximizes the average silhouette is the optimal number of clusters.
%[text] 
%[text] Calculate the silhouette using Hamming distance (specify Distance="hamming" in bitMat\_double).
logSection("A02", "Section 5: Silhouette Analysis for Selecting k", "Analytics L3");
silScores = zeros(1, K_MAX);
rng(42);
for k = 2:K_MAX
    labels_k = kmeans(bitMat_double, k, ...
        Distance="hamming", Replicates=5, Display="off");
    sil = silhouette(bitMat_double, labels_k, "hamming");
    silScores(k) = mean(sil);
end

figure("Name", "A02 Silhouette Score");
plot(2:K_MAX, silScores(2:K_MAX), "-s", Color=[0.8 0.3 0.2], LineWidth=1.5);
xlabel("Number of Clusters k");
ylabel("Average Silhouette Score");
title("Silhouette Analysis — K-Means on ECFP4 Bit Vectors");
grid on;
logInfo("Displayed silhouette plot");

[bestSil, bestK_idx] = max(silScores(2:K_MAX));
bestK = bestK_idx + 1;
logInfo("Optimal k by silhouette: k = %d  (Average Silhouette = %.3f)", bestK, bestSil);

%[text] Execute the final k-means using the optimal k.
rng(42);
finalLabels_km = kmeans(bitMat_double, bestK, ...
    Distance="hamming", Replicates=10, Display="off");
logInfo("K-Means (k=%d) Final Cluster Assignment:", bestK);
for c = 1:bestK
    members = molNames(finalLabels_km == c);
    logInfo("  Cluster KM%d (%d molecules): %s", c, numel(members), ...
        strjoin(members, ", "));
end
%[text] **💡 Observation Point 5 — Compare k-means and Hierarchical Clustering**
%[text] Compare k-means clusters (KM\*) and hierarchical clusters (H\*) to see which molecules group together. Check if there are molecules that switch clusters between the two methods. Execute `silhouette()` on hierarchical clusters with the same k,
%[text] hier\_labels = cluster(Z\_tree, MaxClust=bestK);
%[text] sil\_hier = silhouette(bitMat\_double, hier\_labels, "hamming");
%[text] mean(sil\_hier)
%[text] Determine which method yields a higher silhouette score for this dataset. Chemical category labels (stimulants, analgesics, etc.) serve as the "correct answer" for clustering. Calculate the Adjusted Rand Index or check the confusion matrix to
%[text] \[~, catIdx\] = ismember(molCats, unique(molCats));
%[text] assess how well k-means clusters reproduce chemical categories.
%[text] (Hint: This is unsupervised learning, so do not expect perfect reproduction)
% ... (Try writing code here)
%%
%[text] ## Section 6: Visualizing Clusters on the Similarity Matrix
%[text] Cluster assignment is complete.
%[text] Finally, let's reorder the similarity matrix by cluster order to visually confirm the quality of clustering.
%[text] 
%[text] ### Concept: Sorted Heatmap — Discovering Block Structures
%[text] By rearranging the rows and columns of the similarity matrix based on cluster assignments, molecules with high similarity form bright square blocks along the diagonal. **Good Clustering:** Tight blocks, low off-block values / **Poor Clustering:** Diffuse, no clear blocks
%[text] 
%[text] Sort with hierarchical clusters (visually more appealing than k-means in heatmaps)
logSection("A02", "Section 6: Visualizing Clusters on the Similarity Matrix", "Analytics L3");
[~, sortByCluster] = sort(clusterLabels_hier);
simMatSorted = simMat(sortByCluster, sortByCluster);
namesSorted  = molNames(sortByCluster);
catsSorted   = molCats(sortByCluster);

figure("Name", "A02 Similarity Matrix (Cluster Order Sorted)");
imagesc(simMatSorted);
colormap(parula);
colorbar;
caxis([0 1]);
title("Similarity Matrix (Hierarchical Cluster Order Sorted)");
xlabel("Molecules (Cluster Order)");
ylabel("Molecules (Cluster Order)");

%[text] Add tick labels (abbreviated)
tick_labels = cellfun(@(n) n(1:min(8,end)), cellstr(namesSorted), ...
    UniformOutput=false);
xticks(1:nMols); xticklabels(tick_labels); xtickangle(90);
yticks(1:nMols); yticklabels(tick_labels);
set(gca, FontSize=6);
%[text] Draw cluster boundaries
boundaries = [0; find(diff(clusterLabels_hier(sortByCluster))); nMols] + 0.5;
hold on;
for b = boundaries'
    xline(b, "w-", LineWidth=1);
    yline(b, "w-", LineWidth=1);
end
hold off;
logInfo("Displayed sorted similarity matrix with cluster boundaries");
%[text] **💡 Observation Point 6 — Selecting Representative Molecules in a Sorted Heatmap**
%[text] Check if the diagonal blocks are consistently bright (high similarity).
%[text] Check if there are faint blocks (low intra-cluster similarity).
%[text] This may indicate the need for more clusters or that the molecule is a true outlier not belonging to a tight group.
%[text] Check the case when sorted by k-means labels instead of hierarchical labels.
%[text] Determine which shows a clearer block structure.
%[text] Practical application: If sending one representative molecule from each cluster for bioassay to reduce costs, which would you choose?
%[text] Common strategy: Choose the molecule closest to the cluster centroid (highest median Tanimoto to all cluster members).
%[text] 
%[text] **Summary**
%[text] - ECFP4 fingerprints encode local substructures as bit vectors.
%[text] - Tanimoto distance ($1 - T$) is used as input for hierarchical clustering.
%[text] - Ward linkage method produces compact and balanced clusters.
%[text] - Silhouette analysis provides a data-driven means to select the optimal $k$.
%[text] - Sorted heatmaps serve as a visual quality check for block structures. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]%   data: {"layout":"inline","rightPanelPercent":40}
%---
