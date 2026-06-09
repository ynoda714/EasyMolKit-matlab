%% A01 Answer: Chemical Space Mapping with PCA
% Reference answer for the "Try it yourself" exercise in a01_chemical_space_pca.m.
% This file is self-contained. It can be executed independently from CSV reading to descriptor calculation (no need to run the main file beforehand).

addpath(genpath("src"));
emk.setup.initPython();
logInfo("A01 Answer: Setup complete");
%%
%[text] ## Let's Try 1: Which descriptor has the largest absolute range?
%[text] Q: Explain why the size of the range is problematic in non-standardized PCA.
DATA_FILE  = "data/list/everyday_chemicals.csv";
DESC_NAMES = ["MolWt", "LogP", "TPSA", "NumHDonors", "NumHAcceptors", ...
              "RingCount", "NumRotatableBonds", "FractionCSP3", "HeavyAtomCount"];

rawTbl = readtable(DATA_FILE, TextType="string");
nMols  = height(rawTbl);
mols   = cell(1, nMols);
valid  = false(1, nMols);
for k = 1:nMols
    try
        mols{k} = emk.mol.fromSmiles(rawTbl.SMILES(k));
        valid(k) = true;
    catch
    end
end
validIdx = find(valid);
descMat  = nan(numel(validIdx), numel(DESC_NAMES));
for j = 1:numel(validIdx)
    s = emk.descriptor.calculate(mols{validIdx(j)}, DESC_NAMES);
    for d = 1:numel(DESC_NAMES)
        descMat(j, d) = s.(DESC_NAMES(d));
    end
end

descRange = max(descMat) - min(descMat);
[maxRange, iMax] = max(descRange);
logInfo("Largest range: %s  (%.2f)", DESC_NAMES(iMax), maxRange);

%[text] List of ranges for each descriptor:
for d = 1:numel(DESC_NAMES)
    logInfo("  %-22s  range = %.3f", DESC_NAMES(d), descRange(d));
end

%[text] **Answer:** In this dataset, MolWt has the largest range (~300 units).
%[text] On the other hand, FractionCSP3 is only in the range of 0–1.
%[text] In non-standardized PCA, the scale of MolWt alone dominates the first principal component.
%[text] The results are influenced by the numerical scale, unrelated to the chemical information content.
%[text] By z-score standardization, each descriptor can contribute proportionally to the variance in the dataset.
%%
%[text] ## Let's Try 2: Verify Standardization; Calculate Correlation Matrix

descMean = mean(descMat, 1);
descStd  = std(descMat, 0, 1);
constMask = descStd < 1e-12;
descMat2  = descMat(:, ~constMask);
descMean2 = descMean(~constMask);
descStd2  = descStd(~constMask);
activeDes = DESC_NAMES(~constMask);
Z = (descMat2 - descMean2) ./ descStd2;

colMeans = mean(Z);
colStds  = std(Z);
logInfo("Column means of Z (should be ~0): max |mean| = %.2e", max(abs(colMeans)));
logInfo("Column standard deviations of Z (should be ~1): max |std-1| = %.2e", max(abs(colStds - 1)));

%[text] Calculate Correlation Matrix:
C = corr(Z);
%[text] Detect the most correlated pair (excluding diagonal):
C_offdiag = C - diag(diag(C));
[maxCorr, linIdx] = max(abs(C_offdiag(:)));
[r, c] = ind2sub(size(C), linIdx);
logInfo("Most correlated pair: %s & %s  (r = %.3f)", ...
    activeDes(r), activeDes(c), C(r, c));

%[text] **Answer:** The mean of each column of Z is ~1e-16 (numerically zero), and the standard deviation is 1.
%[text] MolWt and HeavyAtomCount have a strong correlation (r > 0.95).
%[text] Larger molecules have more atoms, so these two carry almost the same information.
%[text] Therefore, removing one of them will not significantly change the PCA results, reducing redundancy.
%%
%[text] ## Let's Try 3: How many principal components are needed for >= 80% variance?

[~, ~, ~, ~, explained] = pca(Z);

cumVar = cumsum(explained);
nFor80 = find(cumVar >= 80, 1);
logInfo("Number of PCs needed for >= 80%% variance: %d  (cumulative at PC%d: %.1f%%)", ...
    nFor80, nFor80, cumVar(nFor80));

for p = 1:min(5, numel(explained))
    logInfo("  PC%d: %.1f%%  (cumulative: %.1f%%)", p, explained(p), cumVar(p));
end

%[text] **Answer:** In this dataset of 9 descriptors and 30 molecules, usually 3 PCs cover >= 80% of the variance.
%[text] When setting a cutoff threshold, the cumulative threshold method (>= 80%) is clearer than the elbow method.
%[text] PCA can only generate a maximum of min(N, D) = min(30, 9) = 9 meaningful PCs.
%[text] The rank of the 30×9 standardized matrix is at most 9.
%%
%[text] ## Let's Try 4: Identify the molecule farthest from the origin in PC space

[coeff, score, ~, ~, explained] = pca(Z);
validNames = rawTbl.CommonName(validIdx);

dist2 = score(:, 1).^2 + score(:, 2).^2;
[~, outerIdx] = sort(dist2, "descend");

logInfo("Top 5 molecules farthest from the PC origin:");
for k = 1:5
    logInfo("  %d. %-20s  d = %.3f  (PC1=%.2f, PC2=%.2f)", ...
        k, validNames(outerIdx(k)), sqrt(dist2(outerIdx(k))), ...
        score(outerIdx(k), 1), score(outerIdx(k), 2));
end

%[text] **Answer:** Sucrose (molecular weight 342, many OH groups) and ascorbic acid tend to be far in the PC1 direction due to large size and high polarity.
%[text] Small alcohols like ethanol form a cluster near the origin.
%[text] Sucrose is located on the positive PC1 (large/polar), while limonene and ibuprofen are on the negative PC2 (lipophilic, low TPSA) side.
%%
%[text] ## Let's Try 5: Interpretation of PC Loadings; Removing HeavyAtomCount
%[text] Display top contributing factors for PC1 and PC2 loadings:
[~, ord1] = sort(abs(coeff(:, 1)), "descend");
[~, ord2] = sort(abs(coeff(:, 2)), "descend");
logInfo("PC1 main contributing factor: %s (loading=%.3f)", ...
    activeDes(ord1(1)), coeff(ord1(1), 1));
logInfo("PC2 main contributing factor: %s (loading=%.3f)", ...
    activeDes(ord2(1)), coeff(ord2(1), 2));

%[text] Remove HeavyAtomCount and rerun PCA:
noHAC = activeDes ~= "HeavyAtomCount";
Z_noHAC = Z(:, noHAC);
[~, ~, ~, ~, exp2] = pca(Z_noHAC);
logInfo("PC1 explained variance with HeavyAtomCount: %.1f%%", explained(1));
logInfo("PC1 explained variance without HeavyAtomCount: %.1f%%", exp2(1));
logInfo("PC1 variance change: %+.1f pp", exp2(1) - explained(1));

%[text] **Answer:** PC1 is dominated by descriptors like NumHAcceptors, TPSA, and MolWt, which relate to "size/polarity."
%[text] PC2 is dominated by descriptors like LogP and FractionCSP3, which relate to "lipophilicity."
%[text] Even after removing HeavyAtomCount (highly collinear with MolWt),
%[text] the change in PC1 explained variance is less than 2%. It confirms that redundant descriptors add little unique information.
logInfo("A01 Answer complete.");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]%   data: {"layout":"inline","rightPanelPercent":40}
%---
