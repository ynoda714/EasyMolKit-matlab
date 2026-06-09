%[text] # A01: Chemical Space Mapping by PCA
%[text] EasyMolKit Analytics — Layer 3
%[text] 
%[text] When 30 types of everyday chemicals (soaps, analgesics, caffeine derivatives, etc.) are quantified by 9 properties, they are represented as points in a 9-dimensional space.
%[text] However, it is difficult for humans to visually understand a 9-dimensional space directly.
%[text] By using PCA (Principal Component Analysis), the nature of the information in the data (variance, i.e., data dispersion) can be compressed into a 2D map (plane) with minimal loss.
%[text] By color-coding chemical categories, you can instantly grasp the relationship between structure and properties, seeing which molecules appear "close" and which are "isolated."
%[text] In this script, you will experience standardizing the descriptor matrix, performing PCA, and interpreting the loadings.
%[text] 
%[text] **What you will learn in this tutorial**
%[text] - Understanding why standardization is necessary before PCA
%[text] - How to use MATLAB's `pca()` on a molecular descriptor matrix
%[text] - Reading a scree plot and determining the number of necessary principal components (PCs)
%[text] - Interpreting principal component loadings and understanding descriptors driving variance
%[text] - Linking the geometric meaning of PCA to actual chemical differences \
%[text] 
%[text] **Prerequisites**
%[text] - Completion of F02 (Property Calculation) — Basic knowledge of descriptors
%[text] - Statistics and Machine Learning Toolbox (`pca` function)
%[text] - No internet connection required \
%[text] 
%[text] Estimated time required: 30–45 minutes | How to run: Execute each section with Ctrl+Enter
%[text] 
%[text] **Data:** `data/list/everyday_chemicals.csv` — 30 types of everyday chemicals (PubChem CC0)
%[text] 
%[text] **References**
%[text] - Wold S, Esbensen K, Geladi P (1987) Principal component analysis. Chemometrics Intell Lab Syst 2:37-52.
%[text] - RDKit: https://www.rdkit.org/docs/GettingStartedInPython.html\#list-of-available-descriptors \
%%
%[text] ## Section 0: Setup
logSection("A01", "Section 0: Setup", "Analytics L3");
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
%[text] Prepare Python and RDKit processes before main execution
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;
%%
%[text] ## Section 1: Loading Everyday Chemicals and Calculating Descriptors
%[text] 
%[text] ### Concept: Molecular Descriptors for Chemical Space
%[text] Molecular descriptors are numerical representations of molecular structures.
%[text] They convert 2D graphs of atoms and bonds into fixed-length vectors.
%[text] The 9 descriptors used in this tutorial cover four property dimensions.
%[text] 
%[text] **Size Dimension**
%[text] - `MolWt` — Molecular weight (g/mol). Larger molecules tend to occupy receptors more easily but are less likely to be absorbed (Ro5: MW ≤ 500).
%[text] - `HeavyAtomCount` — Number of non-hydrogen atoms. Correlates with molecular complexity. \
%[text] 
%[text] **Lipophilicity Dimension**
%[text] - LogP — Logarithm of the partition coefficient between water and octanol. Higher values indicate higher lipophilicity and lower water solubility (Oral bioavailability guideline: Lipinski's Ro5 suggests LogP ≤ 5).
%[text] - `FractionCSP3` — Fraction of sp3 carbons. High Fsp3 indicates 3D shape, leading to improved solubility and selectivity (Beyond-Ro5 concept). \
%[text] **Polarity Dimension**
%[text] - TPSA — Topological polar surface area (Å²). Larger values tend to reduce membrane permeability and oral absorption (Guideline for CNS-active drugs: TPSA < 90).
%[text] - NumHDonors — Number of hydrogen bond donors (e.g., NH, OH groups). Ro5 guideline: HBD ≤ 5
%[text] - NumHAcceptors — Number of hydrogen bond acceptors (e.g., N, O atoms). Ro5 guideline: HBA ≤ 10. \
%[text] **Topology Dimension**
%[text] - `RingCount` — Number of SSSR rings. Higher values in aromatic drugs.
%[text] - `NumRotatableBonds` — Indicator of flexibility. Higher values reduce oral bioavailability. \
logSection("A01", "Section 1: Loading Everyday Chemicals and Calculating Descriptors", "Analytics L3");
DATA_FILE    = "data/list/everyday_chemicals.csv";
DESC_NAMES   = ["MolWt", "LogP", "TPSA", "NumHDonors", "NumHAcceptors", ...
                "RingCount", "NumRotatableBonds", "FractionCSP3", "HeavyAtomCount"];
N_DESCRIPTORS = numel(DESC_NAMES);
%[text] Read CSV
rawTbl = readtable(DATA_FILE, TextType="string");
logInfo("Loaded %d molecules from %s", height(rawTbl), DATA_FILE);
%[text] Convert SMILES to Mol objects
nMols = height(rawTbl);
mols  = cell(1, nMols);
valid = false(1, nMols);
for k = 1:nMols
    try
        mols{k} = emk.mol.fromSmiles(rawTbl.SMILES(k));
        valid(k) = true;
    catch ME
        logWarn("Skipping %s: %s", rawTbl.CommonName(k), ME.message);
    end
end
logInfo("Successfully parsed %d / %d molecules", sum(valid), nMols);
%[text] Calculate descriptors for valid molecules only
validIdx   = find(valid);
validNames = rawTbl.CommonName(validIdx);
validCats  = rawTbl.Category(validIdx);

descMat = nan(numel(validIdx), N_DESCRIPTORS);   % Rows = molecules, Columns = descriptors
for j = 1:numel(validIdx)
    s = emk.descriptor.calculate(mols{validIdx(j)}, DESC_NAMES);
    for d = 1:N_DESCRIPTORS
        descMat(j, d) = s.(DESC_NAMES(d));
    end
end

logInfo("Descriptor matrix: %d molecules x %d descriptors", ...
    size(descMat, 1), size(descMat, 2));
%[text] Display the first 5 rows as a table
previewTbl = array2table(descMat(1:5, :), VariableNames=cellstr(DESC_NAMES));
previewTbl.Name = validNames(1:5);
previewTbl = movevars(previewTbl, "Name", Before=1);
disp(previewTbl);
%[text] **💡 Observation Point 1**
%[text] Before standardizing in the next section, let's check the scale of each descriptor.
%[text] 
%[text] Identify the descriptor with the largest range among 30 molecules. Use `range(descMat)` or `max(descMat) - min(descMat)`.
%[text] 
%[text] Consider why performing PCA on non-standardized data can be problematic due to large range differences. Hint: PCA maximizes variance. A descriptor with a range 10 times larger will dominate the first principal component, regardless of its chemical significance.
%%
%[text] ## Section 2: Standardization of Descriptor Matrix
%[text] In Section 1, we quantified 30 types of molecules using 9 descriptors.
%[text] Before performing PCA, it is necessary to convert the matrix to unify the scale of the descriptors.
%[text] 
%[text] ### Concept: Why z-score standardization is needed before PCA
%[text] PCA seeks the direction of maximum variance in the data.
%[text] If the scales of the descriptors differ significantly (e.g., MolWt: 40–350, FractionCSP3: 0–1), descriptors with larger scales may dominate the first principal component due to unit differences.
%[text] 
%[text] **Formula for z-score standardization (normalization):** $z = (x - \\mu) / \\sigma$ (subtract the mean from each value and divide by the standard deviation)
%[text] 
%[text] **Effects after standardization:**
%[text] - The mean of all descriptors becomes 0, and the standard deviation becomes 1
%[text] - PCA is unaffected by scale and can choose directions that explain chemical variance
%[text] - Loadings become directly comparable, revealing which descriptors are more important \
%[text] 
%[text] **When you can skip standardization:**
%[text] - When all descriptors are on the same scale (e.g., binary bits)
%[text] - When the scale itself has meaning (rare in cheminformatics) \
logSection("A01", "Section 2: Standardization of Descriptor Matrix", "Analytics L3");
descMean = mean(descMat, 1);            % 1 x D
descStd  = std(descMat,  0, 1);        % 1 x D  (normalized by N-1)
%[text] Guard: Exclude constant descriptors (division by zero occurs if std is 0)
constMask = descStd < 1e-12;
if any(constMask)
    logWarn("Removed constant descriptors: %s", ...
        strjoin(DESC_NAMES(constMask), ", "));
    descMat  = descMat(:,  ~constMask);
    descMean = descMean(:, ~constMask);
    descStd  = descStd(:,  ~constMask);
    activeDes = DESC_NAMES(~constMask);
else
    activeDes = DESC_NAMES;
end

Z = (descMat - descMean) ./ descStd;   % Standardized matrix (N x D)
logInfo("Standardization complete: Mean ~ %.2e, Std Dev ~ %.2f", mean(abs(mean(Z))), mean(std(Z)));
%[text] **💡 Observation Point 2**
%[text] Verify that the mean of each column of Z is approximately 0 and the standard deviation is approximately 1. You can check using `mean(Z)` and `std(Z)`.
%[text] 
%[text] Calculate the correlation matrix of Z and identify the most correlated descriptor pairs: `corr(Z)`. For example, MolWt and HeavyAtomCount are expected to be strongly correlated. High correlation suggests what information two descriptors add in PCA.
%%
%[text] ## Section 3: Performing PCA and Checking Explained Variance
%[text] The standardized matrix Z is ready.
%[text] Now, execute `pca()` to compress the 9-dimensional descriptor space into 2-3 dimensions.
%[text] 
%[text] ### Concept: What PCA Computes
%[text] PCA is an algorithm that sequentially finds the directions in which the data spreads the most.
%[text] For the standardized matrix Z, new axes (principal components) are determined in the following order:
%[text] - **PC1 (First Principal Component)**: The direction with the greatest differences between molecules
%[text] - **PC2 (Second Principal Component)**: Orthogonal to PC1, the next direction with the greatest variance
%[text] - Subsequent axes are also orthogonal and maximize variance \
%[text] 
%[text] The percentage of total variance each principal component explains is called the "explained variance ratio," which can be visually confirmed with a scree plot.
%[text] 
%[text] **Return Values of `pca()`**:
%[text] - `coeff` — D × D loading matrix (columns = PC directions in the original space)
%[text] - `score` — N × D score matrix (coordinates of each molecule in PC space)
%[text] - `latent` — D × 1 eigenvalues (variance for each PC)
%[text] - `explained` — D × 1 explained variance (%) \
logSection("A01", "Section 3: Performing PCA and Checking Explained Variance", "Analytics L3");
[coeff, score, latent, ~, explained] = pca(Z);

logInfo("PCA completed. Explained variance of the first 3 PCs:");
logInfo("  PC1: %.1f%%", explained(1));
logInfo("  PC2: %.1f%%", explained(2));
logInfo("  PC3: %.1f%%  (Cumulative: %.1f%%)", explained(3), sum(explained(1:3)));
%[text] \-- Scree Plot --
nPCs = numel(explained);
figure("Name", "A01 Scree Plot");
bar(1:nPCs, explained, FaceColor=[0.3 0.6 0.9]);
hold on;
plot(1:nPCs, cumsum(explained), "-o", Color=[0.8 0.2 0.2], LineWidth=1.5);
yline(80, "--k", "80% Threshold", LabelHorizontalAlignment="left");
hold off;
xlabel("Principal Component");
ylabel("Explained Variance (%)");
title("Scree Plot -- Household Chemicals PCA");
legend("Individual", "Cumulative", Location="east");
grid on;
%[text] **💡 Observation Point 3**
%[text] Check how many PCs are needed to explain more than 80% of the total variance. Consider whether to use the elbow method (scree drop-off point) or the cumulative threshold method (>= 80%).
%[text] 
%[text] The dataset of 30 molecules has 9 types of descriptors (D=9). Consider whether PCA can generate more than 9 meaningful PCs. (Hint: The rank of an $N \\times D$ matrix is at most $\\min(N, D)$)
%%
%[text] ## Section 4: Visualization of 2D Chemical Space
%[text] PCA (Principal Component Analysis) scores have been obtained.
%[text] In this section, we will plot 30 molecules on a 2D plane using PC1 and PC2 as axes and examine the distribution by chemical category.
%[text] 
%[text] ### Concept: Scatter Plot of Chemical Space
%[text] By plotting PC1 and PC2, each molecule can be projected onto the plane of maximum variance.
%[text] Molecules forming clusters have similar descriptor profiles and are located in the same region of chemical space.
%[text] 
%[text] By coloring by chemical category, we can check if chemical classes naturally separate.
%[text] - **Complete Separation** → Descriptors capture the differences between classes well.
%[text] - **Overlap** → Classes share physical properties. \
%[text] 
%[text] **Note:** PCA is a linear method and cannot capture the nonlinear structure of chemical space (e.g., the "shape" of drug-like space).
%[text] Nonlinear alternatives include UMAP and t-SNE (`tsne` is available from R2017b, `umap` from R2024a).
logSection("A01", "Section 4: Visualization of 2D Chemical Space", "Analytics L3");
categories = unique(validCats);
nCats      = numel(categories);
cmap       = lines(nCats);    % Colormap for identification

figure("Name", "A01 Chemical Space (PC1 vs PC2)");
hold on;
hGroups = gobjects(nCats, 1);
for c = 1:nCats
    mask = validCats == categories(c);
    hGroups(c) = scatter(score(mask, 1), score(mask, 2), 60, ...
        "filled", MarkerFaceColor=cmap(c,:), DisplayName=categories(c));
    % Label each point with the molecule's short name
    for k = find(mask)'
        shortName = extractBefore(validNames(k) + "          ", 11);
        text(score(k, 1) + 0.05, score(k, 2) + 0.05, shortName, ...
            FontSize=7, Color=cmap(c,:), Interpreter="none");
    end
end
hold off;

xlabel(sprintf("PC1 (%.1f%% Variance)", explained(1)));
ylabel(sprintf("PC2 (%.1f%% Variance)", explained(2)));
title("2D Chemical Space -- Household Chemicals");
legend(hGroups, Location="bestoutside", Interpreter="none");
grid on;
%[text] **💡 Observation Point 4**
%[text] Check if molecules of the same chemical category (e.g., stimulants, analgesics) form clusters. Which category shows the most densely packed cluster?
%[text] 
%[text] 
%[text] Consider where sucrose (molecular weight 342, many OH groups) and ethanol (molecular weight 46, 1 OH group) are located on the PCA plot.
%%
%[text] ## Section 5: Interpretation of PC Loadings
%[text] We have checked the "position" of molecules on the scatter plot. Next, to understand "why they are in that position," we interpret which descriptors each PC axis corresponds to using loadings.
%[text] 
%[text] ### Concept: How to Read a Loading Plot (Biplot)
%[text] The loading vector of PC1 indicates "which original descriptors point in the direction of maximum variance."
%[text] 
%[text] - **Positive Loading:** High values of that descriptor move molecules in the positive PC1 direction.
%[text] - **Large Absolute Loading:** That descriptor strongly drives PC1. \
%[text] 
%[text] A biplot overlays scores (molecule positions) and loadings (descriptor arrows) in the same space for visual interpretation.
%[text] 
%[text] **Quick Interpretation:**
%[text] - Molecules along the positive PC1 axis have high values of positively loaded descriptors.
%[text] - Descriptors pointing in similar directions are correlated with each other.
%[text] - Orthogonal loading vectors indicate uncorrelated descriptors. \
logSection("A01", "Section 5: Interpretation of PC Loadings", "Analytics L3");
nActive = numel(activeDes);
figure("Name", "A01 PC Loadings (PC1 vs PC2)");
%[text] \-- Loading Bar Graph for PC1 and PC2 --
subplot(1, 2, 1);
barh(coeff(:, 1), FaceColor=[0.3 0.6 0.9]);
yticks(1:nActive);
yticklabels(cellstr(activeDes));
xlabel("Loading on PC1");
title(sprintf("PC1 (%.1f%% Variance)", explained(1)));
xline(0, "k");
grid on;

subplot(1, 2, 2);
barh(coeff(:, 2), FaceColor=[0.9 0.5 0.3]);
yticks(1:nActive);
yticklabels(cellstr(activeDes));
xlabel("Loading on PC2");
title(sprintf("PC2 (%.1f%% Variance)", explained(2)));
xline(0, "k");
grid on;

sgtitle("PC Loadings -- Everyday Chemicals");
%[text] Display top 3 contributors to PC1 and PC2
[~, ord1] = sort(abs(coeff(:, 1)), "descend");
[~, ord2] = sort(abs(coeff(:, 2)), "descend");
logInfo("Top 3 contributors to PC1: %s, %s, %s", ...
    activeDes(ord1(1)), activeDes(ord1(2)), activeDes(ord1(3)));
logInfo("Top 3 contributors to PC2: %s, %s, %s", ...
    activeDes(ord2(1)), activeDes(ord2(2)), activeDes(ord2(3)));
%[text] **💡 Observation Point 5**
%[text] Let's check the descriptors that contribute most to PC1. Consider if they make chemical sense. (Hint: PC1 often corresponds to the "size" axis in drug-like chemical space.)
%[text] 
%[text] Let's check the descriptors that contribute most to PC2. Can you identify the characteristics that separate the positive and negative sides of the PC2 axis? (Hint: Check which categories of molecules are located on the positive and negative sides of PC2.)
%[text] 
%[text] Observe the biplot displaying scores and loadings together: `biplot(coeff(:,1:2), Scores=score(:,1:2), VarLabels=cellstr(activeDes))` Are molecules of different categories pulled in the direction of their characteristic descriptors?
%[text] 
%[text] Try removing `HeavyAtomCount` from `DESC_NAMES` (it is strongly correlated with MolWt). If you rerun PCA, does the explained variance of PC1 change significantly? What does this tell you about redundant descriptors?
%[text] **Summary**
%[text] - By z-score standardizing descriptors before PCA, you can remove scale bias.
%[text] - Use a scree plot to check cumulative explained variance and select the number of PCs covering 80%.
%[text] - By interpreting PC loadings, you can understand which descriptors drive the variance in chemical space.
%[text] - In descriptor-based PCA, molecules of the same chemical category tend to form spatial clusters. \

%[appendix]{"version":"1.0"}
