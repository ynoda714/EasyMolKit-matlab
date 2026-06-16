%[text] # A08: Mass Spectrometry x Cheminformatics
%[text] EasyMolKit Analytics -- Layer 3
%[text] Five drug samples have arrived at the pharmaceutical QC lab. The labels are damaged by water and unreadable. Can you identify "what drug is this" using only LC-MS data?
%[text] Mass Spectrometry (MS) is an analytical technique that ionizes molecules and measures their precise mass at the ppm level.
%[text] Isotope patterns -- spectral fingerprints of $M$, $M+1$, and $M+2$ -- allow us to narrow multiple mass-matched candidates down to a single compound.
%[text] In this tutorial, we work through all three stages: low-resolution search, high-precision mass, and isotope confirmation.
%[text] Identification unfolds in three stages with increasing confidence.
%[text:table]{"ignoreHeader":true}
%[text] | Stage | Description |
%[text] | --- | --- |
%[text] | Stage 1 -- Low-resolution search (+-0.5 Da) | In real commercial databases, dozens of matches occur; in this 200-compound reference, expect 1-2. |
%[text] | Stage 2 -- High-resolution precise mass (5 ppm) | Only a few candidates survive. |
%[text] | Stage 3 -- Isotope pattern confirmation | The fingerprints of $M$, $M+1$, $M+2$ narrow each unknown substance down to one compound. |
%[text:table]
%[text] ### Exercise Content
%[text:table]{"ignoreHeader":true}
%[text] | Step | Description |
%[text] | --- | --- |
%[text] | 1 | Use `emk.descriptor.calculate()` to build a precise mass reference table from 200 FDA-approved drugs. |
%[text] | 2 | Simulate ESI-MS spectra for 5 drugs (Gaussian peaks, isotope clusters, random noise). |
%[text] | 3 | Detect peaks using `smoothdata` and `findpeaks` from the Signal Processing Toolbox. |
%[text] | 4 | Search the reference table with mass tolerance and compare the number of low and high-resolution candidates. |
%[text] | 5 | Calculate theoretical isotope patterns (relative intensities of $M$, $M+1$, $M+2$) from molecular formulas. |
%[text] | 6 | Score candidates with cosine similarity and create a final ranked list. |
%[text:table]
%[text] ### Learning Objectives
%[text] - Explain why precise mass measurement reduces ambiguity in compound ID.
%[text] - Use `findpeaks` and `smoothdata` on spectral data.
%[text] - Calculate ppm mass error and set appropriate search tolerances.
%[text] - Derive $M+1$ and $M+2$ isotope intensities from natural elemental abundance.
%[text] - Score candidates with cosine similarity and build a ranked hit list. \
%[text] ### Prerequisites
%[text] - Completion of A03 (QSAR Regression) -- Understand the basics of `emk.descriptor.calculate()`.
%[text] - Completion of S05 (Unknown Substance ID) -- Understand the concept of similarity-based identification.
%[text] - Ability to use the Signal Processing Toolbox (`findpeaks`, `smoothdata`). \
%[text] ### Environment
%[text] - Signal Processing Toolbox is required for Sections 3 and 4.
%[text] - Statistics and ML Toolbox (`corrcoef`) is used in Section 5.
%[text] - Estimated time: 35-50 minutes \
%[text] ### Data: 
%[text] - `data/list/fda_drugs.csv` -- 200 FDA-approved drugs (ChEMBL, CC-BY-SA 3.0)
%[text] - Columns: ChEMBLID, Name, SMILES, MolecularWeight, ALogP, HBondDonors, HBondAcceptors, TPSA, RotatableBonds \
%[text] ### References
%[text] - Gross JH (2011) *Mass Spectrometry: A Textbook*, 2nd ed. Springer. ISBN 978-3-642-10709-2 -- Monoisotopic mass, isotope patterns, ESI-MS
%[text] - Kind T & Fiehn O (2006) Metabolomic database annotations via query of elemental compositions. *BMC Bioinformatics* 7:234. doi:10.1186/1471-2105-7-234
%[text] - Claesen J et al. (2012) Efficient method for isotopic distribution calculation. *J Am Soc Mass Spectrom* 23:753-763. doi:10.1007/s13361-011-0326-2
%[text] - Stein SE & Scott DR (1994) Optimization and testing of mass spectral library search algorithms. *J Am Soc Mass Spectrom* 5:859-866. doi:10.1016/1044-0305(94)87009-8
%[text] - Meija J et al. (2016) Atomic weights of the elements 2013 (IUPAC Technical Report). *Pure Appl. Chem.* 88:265-291. doi:10.1515/pac-2015-0305 \
%[text] How to run: Execute each section one by one with Ctrl+Enter.
%%
%[text] ## Section 0: Setup
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
emk.setup.initPython(); %[output:930256b9]
%[text] Warm up Python and RDKit before the main process.
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;
logSection("A08", "Section 0: Setup", "Analytics L3"); %[output:054dd858]
%%
%[text] ## Section 1: Building a Precise Mass Reference Table
%[text] We begin by building a precise mass reference table from 200 FDA-approved drugs -- our "dictionary" for all database searches in this tutorial.
%[text] ### Monoisotopic Mass vs Average Molecular Weight
%[text] Two distinct mass quantities arise in chemistry; only exact monoisotopic mass belongs in an MS database search.
%[text] **Average Molecular Weight (MolWt)** is the weighted average mass of all stable isotopes weighted by their natural abundance.
%[text] For example, the average atomic weight of carbon is 12.011 g/mol (a mixture of 98.9% $^{12}$C and 1.1% $^{13}$C).
%[text] This is the value chemists use for stoichiometry and solution preparation.
%[text] **Exact Monoisotopic Mass (ExactMolWt)** is the mass of a molecule composed only of the most abundant isotopes of each element ($^{12}$C, $^{1}$H, $^{14}$N, $^{16}$O, $^{32}$S, $^{35}$Cl, etc.).
%[text] Mass spectrometers measure the mass-to-charge ratio ($m/z$) of individual ions, so ExactMolWt is used in MS.
%[text] Why is the distinction important?
%[text] For Aspirin (C9H8O4), MolWt = 180.16 g/mol, ExactMolWt = 180.0423 Da, with a difference of 0.12 Da, leading to database mismatches if searched with the wrong mass type.
%[text] 
%[text] In positive mode ESI-MS (electrospray ionization), the $\[M+H\]^+$ ion with an added proton is observed.
%[text] The observed $m/z = \\text{ExactMolWt} + 1.00728$ Da (proton mass).
logSection("A08", "Section 1: Building a Precise Mass Reference Table", "Analytics L3"); %[output:7f13b345]
DATA_FILE = "data/list/fda_drugs.csv";

rawTbl = readtable(DATA_FILE, TextType="string");
nRaw   = height(rawTbl);
logInfo("Read %d rows from %s", nRaw, DATA_FILE); %[output:337aef13]
%[text] Calculate the exact mass and molecular formula for each compound.
exactMass = nan(nRaw, 1);
formula   = strings(nRaw, 1);
valid     = false(1, nRaw);

logInfo("Calculating exact mass (may take 1-2 minutes)..."); %[output:75d7087d]
for k = 1:nRaw
    try
        mol          = emk.mol.fromSmiles(rawTbl.SMILES(k));
        d            = emk.descriptor.calculate(mol, ["ExactMolWt","MolFormula"]);
        exactMass(k) = d.ExactMolWt;
        formula(k)   = d.MolFormula;
        valid(k)     = true;
    catch ME
        logWarn("Row %d (%s): %s", k, rawTbl.Name(k), ME.message);
    end
end

validIdx = find(valid);
refTbl = table( ...
    rawTbl.Name(validIdx), ...
    rawTbl.SMILES(validIdx), ...
    formula(validIdx), ...
    exactMass(validIdx), ...
    VariableNames=["Name","SMILES","Formula","ExactMass"]);
nRef = height(refTbl);
logInfo("Reference table: %d compounds", nRef); %[output:85be30ce]
%[text] Preview the first 5 entries.
disp(refTbl(1:5, ["Name","Formula","ExactMass"])); %[output:1c8b109d]
%%
%[text] ## Section 2: Simulation of 5 Unknown ESI-MS Spectra
%[text] With the reference table built, we now simulate realistic ESI-MS spectra for all 5 unknown drugs.
%[text] ### ESI-MS Spectra and Isotope Clusters
%[text] Mass spectra form isotope clusters around $\[M+H\]^+$, consisting of $M$, $M+1$, and $M+2$.
%[text] An isotope cluster is a group of peaks including the $M$ peak with only the lightest isotopes, and $M+1$ and $M+2$ peaks containing one or two heavy isotopes (such as $^{13}$C, $^{34}$S, $^{37}$Cl).
%[text] It also includes fragment ions from neutral losses like water (18 Da) or CO2 (44 Da), and baseline noise.
%[text] The simulation is reproducible with a fixed seed.
logSection("A08", "Section 2: Simulation of 5 Unknown ESI-MS Spectra", "Analytics L3"); %[output:44614b03]
PROTON_MASS  = 1.007276;   % Da
MZ_SIGMA_LR  = 0.3;        % Da -- Simulated low-resolution peak width (FWHM ~0.7 Da)
NOISE_LEVEL  = 0.03;       % Proportion relative to base peak
N_UNKNOWNS   = 5;
%[text] Select N\_UNKNOWNS drugs evenly distributed across the reference table (fixed selection -- independent of rng).
unknownIdx = round(linspace(ceil(nRef * 0.05), floor(nRef * 0.95), N_UNKNOWNS));

unknownNames    = refTbl.Name(unknownIdx);
unknownFormulas = refTbl.Formula(unknownIdx);
unknownMass     = refTbl.ExactMass(unknownIdx);

logInfo("Selected unknown substances:"); %[output:46addd2a]
for k = 1:N_UNKNOWNS %[output:group:79a78ecc]
    logInfo("  Unknown substance %d: %s  (%.4f Da, molecular formula %s)", ... %[output:0103aab5]
        k, unknownNames(k), unknownMass(k), unknownFormulas(k)); %[output:0103aab5]
end %[output:group:79a78ecc]
%[text] Generate spectra (struct array: .mz, .intensity).
rng(2026);   % Fixed seed for reproducible noise
spectra = cell(1, N_UNKNOWNS);
for k = 1:N_UNKNOWNS
    mH    = unknownMass(k) + PROTON_MASS;  % [M+H]+ m/z
    fc    = parseFormula(unknownFormulas(k));
    iso   = isoPattern(fc);                 % [M, M+1, M+2] relative intensities

    % Peak positions and relative intensities
    pkMz  = [mH,          mH+1,       mH+2, ...   % Isotope cluster
             mH-18,       mH-44];                  % Common neutral losses
    pkInt = [iso(1),      iso(2),     iso(3), ...
             iso(1)*0.35, iso(1)*0.15];

    % Retain only peaks within valid m/z range (above 50 Da)
    keep  = pkMz > 50;
    spectra{k} = simulateSpectrum(pkMz(keep), pkInt(keep), ...
                     [max(50, mH-80), mH+10], MZ_SIGMA_LR, NOISE_LEVEL);
end
%[text] Display the spectrum for unknown substance 1.
figure("Name","A08: Sample Spectrum (Unknown Substance 1)"); %[output:6a938815]
plot(spectra{1}.mz, spectra{1}.intensity, "b-", LineWidth=0.8); %[output:6a938815]
xlabel("m/z"); ylabel("Relative Intensity"); %[output:6a938815]
title(sprintf("Simulated ESI-MS Spectrum -- Unknown Substance 1 (Labels Hidden)")); %[output:6a938815]
grid on; %[output:6a938815]
%%
%[text] ## Section 3: Peak Detection Using Signal Processing Toolbox
%[text] With spectra generated, we now use `smoothdata` and `findpeaks` from the Signal Processing Toolbox to extract the $m/z$ of $\[M+H\]^+$.
%[text] ### Roles of `smoothdata` and `findpeaks`
%[text] `smoothdata(..., "gaussian")` reduces noise with Gaussian weighted averaging, improving detection accuracy without distorting peak shapes.
%[text] The main parameters of `findpeaks` are shown below.
%[text:table]{"ignoreHeader":true}
%[text] | Parameter | Meaning | Setting in this exercise |
%[text] | --- | --- | --- |
%[text] | MinPeakHeight | Minimum peak height (normalized intensity) | 0\.05 (noise removal) |
%[text] | MinPeakProminence | Prominence from surrounding valleys | 0\.04 |
%[text] | MinPeakDistance | Minimum distance between adjacent peaks (Da) | 0\.8 (set below 1 to resolve isotopic peaks 1 Da apart) |
%[text:table]
%[text] Isotopic peaks are 1 Da apart, so MinPeakDistance must be set below 1 Da.
logSection("A08", "Section 3: Peak Detection Using Signal Processing Toolbox", "Analytics L3"); %[output:8253ebe1]
MIN_PEAK_HEIGHT      = 0.05;   % 5% of normalized intensity
MIN_PEAK_PROMINENCE  = 0.04;
MIN_PEAK_DISTANCE    = 0.8;    % Da -- needs to be below 1 to resolve isotopic peaks
%[text] Demonstration with Unknown Substance 1
sp1       = spectra{1};
smoothInt = smoothdata(sp1.intensity, "gaussian", 5);

[pkHeight, pkLoc] = findpeaks(smoothInt, sp1.mz, ...
    "MinPeakHeight",     MIN_PEAK_HEIGHT, ...
    "MinPeakProminence", MIN_PEAK_PROMINENCE, ...
    "MinPeakDistance",   MIN_PEAK_DISTANCE);
%[text] Normalize the peak intensities.
pkHeight = pkHeight / max(pkHeight);

logInfo("Unknown Substance 1: Detected %d peaks", numel(pkLoc)); %[output:140bf7ea]
if numel(pkLoc) > 0 %[output:group:384b4a5e]
    logInfo("  Detected m/z values: %s", strjoin(string(round(pkLoc, 3)), ", ")); %[output:86391e69]
end %[output:group:384b4a5e]
%[text] $\[M+H\]^+$ is the strongest peak within the high $m/z$ cluster (within 3 Da of the highest detected $m/z$).
%[text] For compounds containing Cl or S, the $M+2$ isotopic peak becomes strong, so select the peak with the maximum intensity in the cluster as the $M$ peak.
highMzMask  = pkLoc >= max(pkLoc) - 3;
clusterHts  = pkHeight(highMzMask);
clusterLocs = pkLoc(highMzMask);
[~, iM]     = max(clusterHts);
mH_obs_1    = clusterLocs(iM);
logInfo("Unknown Substance 1: Extracted [M+H]+ = %.4f Da", mH_obs_1); %[output:9db13e61]
%[text] Visualize raw data, smoothed data, and detected peaks.
figure("Name","A08: Peak Detection (Unknown Substance 1)"); %[output:11cc638e]
plot(sp1.mz, sp1.intensity, Color=[0.7 0.7 0.7], DisplayName="Raw Data"); hold on; %[output:11cc638e]
plot(sp1.mz, smoothInt,     "b-", LineWidth=1.2, DisplayName="Smoothed"); %[output:11cc638e]
plot(pkLoc,  pkHeight,      "rv", MarkerSize=8, MarkerFaceColor="r", ... %[output:11cc638e]
     DisplayName="Detected Peaks"); %[output:11cc638e]
xlabel("m/z"); ylabel("Relative Intensity"); %[output:11cc638e]
title("Unknown Substance 1 -- Raw Data, Smoothed, Detected Peaks"); %[output:11cc638e]
legend; grid on; %[output:11cc638e]
%[text] Extract approximate $\[M+H\]^+$ for all unknown substances.
%[text] With a grid step of about 0.018 Da, the accuracy is around 50-100 ppm.
%[text] This is sufficient for low-resolution searches but insufficient for high-resolution searches at 5 ppm.
mH_obs = nan(1, N_UNKNOWNS);
for k = 1:N_UNKNOWNS
    sp = spectra{k};
    sm = smoothdata(sp.intensity, "gaussian", 5);
    [pks, locs] = findpeaks(sm, sp.mz, ...
        "MinPeakHeight",     MIN_PEAK_HEIGHT, ...
        "MinPeakProminence", MIN_PEAK_PROMINENCE, ...
        "MinPeakDistance",   MIN_PEAK_DISTANCE);
    if ~isempty(locs)
        highMzMask  = locs >= max(locs) - 3;   % Molecular ion cluster window
        clusterLocs = locs(highMzMask);
        [~, iM]     = max(pks(highMzMask));     % M peak is stronger than M+2
        mH_obs(k)   = clusterLocs(iM);
    end
end
%[text] Simulate high-resolution mass readings (used in Sections 4-6).
%[text] Actual Orbitrap or Q-TOF reports mass with about 1 ppm rms due to internal lock mass calibration.
%[text] Here, model by adding Gaussian(0, 1 ppm) noise to the true $\[M+H\]^+$.
rng(99, "twister");
mH_hires = (unknownMass + PROTON_MASS)' + ...
           (unknownMass + PROTON_MASS)' .* randn(1, N_UNKNOWNS) * 1e-6;
%%
%[text] ## Section 4: Accurate Mass Database Search
%[text] With peak $m/z$ values extracted, we narrow down candidates by matching against the reference table.
%[text] ### Mass Accuracy in ppm (parts per million)
%[text] ppm error is calculated as $\\frac{(m/z\_{obs} - m/z\_{theo})}{m/z\_{theo}} \\times 10^6$.
%[text] The accuracy of mass spectrometers varies greatly depending on the type of instrument.
%[text:table]{"ignoreHeader":true}
%[text] | Instrument Type | Mass Accuracy | Identification Power |
%[text] | --- | --- | --- |
%[text] | Low Resolution (Single Quadrupole, Ion Trap) | \+-0.2-0.5 Da (~ 700-1700 ppm) | Over 40 matches at nominal mass 300 |
%[text] | High Resolution (Orbitrap, Q-TOF, FT-ICR) | \<5 ppm (at $m/z$ 300 \<0.0015 Da) | Can narrow down to 1-2 candidates |
%[text:table]
%[text] 
%[text] We mimic two types of MS instruments here.
%[text] `mH_obs(k)` is the `findpeaks` centroid (accuracy ~ 50-100 ppm) used for low-resolution Da window search.
%[text] `mH_hires(k)` is a simulated Orbitrap value (true $\[M+H\]^+$ with 1 ppm noise added) used for high-resolution ppm window search.
%[text] This separation is physically accurate.
%[text] Unit resolution detectors report peak positions with about 0.5 Da accuracy, while high-resolution instruments provide precise masses with internal calibration of about 1-5 ppm.
logSection("A08", "Section 4: Accurate Mass Database Search", "Analytics L3"); %[output:9c31587c]
LOW_RES_TOL  = 0.50;   % Da
HIGH_RES_PPM = 5.0;    % ppm

candidateCounts = zeros(N_UNKNOWNS, 2);  % [Low Resolution, High Resolution]

for k = 1:N_UNKNOWNS %[output:group:1b44df07]
    % Low-resolution search (Da window) -- using findpeaks extraction (~50-100 ppm)
    lowResCand = [];
    if ~isnan(mH_obs(k))
        mNeutral_lr = mH_obs(k) - PROTON_MASS;
        daDiff      = abs(refTbl.ExactMass - mNeutral_lr);
        lowResCand  = find(daDiff <= LOW_RES_TOL);
    end

    % High-resolution search (ppm window) -- using calibrated instrument readings (~1 ppm)
    mNeutral_hr = mH_hires(k) - PROTON_MASS;
    ppmDiff     = abs(refTbl.ExactMass - mNeutral_hr) ./ refTbl.ExactMass * 1e6;
    highResCand = find(ppmDiff <= HIGH_RES_PPM);

    candidateCounts(k, :) = [numel(lowResCand), numel(highResCand)];

    logInfo("Unknown substance %d (High-resolution [M+H]+ = %.4f, Neutral = %.4f Da):", ... %[output:7dcbac68] %[output:1b890286] %[output:8bb7d389] %[output:0b035b93] %[output:4bf03f71]
        k, mH_hires(k), mNeutral_hr); %[output:7dcbac68] %[output:1b890286] %[output:8bb7d389] %[output:0b035b93] %[output:4bf03f71]
    logInfo("  Low Resolution (+-%.2f Da): %d candidates", LOW_RES_TOL, numel(lowResCand)); %[output:29e6414f] %[output:53229f52] %[output:725e94cf] %[output:6961240c] %[output:7beb5664]
    logInfo("  High Resolution (%.1f ppm): %d candidates", HIGH_RES_PPM, numel(highResCand)); %[output:5fa2fda3] %[output:4e4fd212] %[output:9ff42df6] %[output:24fc1537] %[output:6491ba3f]
end %[output:group:1b44df07]
%[text] Display a bar graph of the number of candidates.
figure("Name","A08: Number of Candidates by Resolution"); %[output:0cfb32e5]
bar(1:N_UNKNOWNS, candidateCounts); %[output:0cfb32e5]
legend(sprintf("Low Resolution (+-%.2f Da)", LOW_RES_TOL), ... %[output:0cfb32e5]
       sprintf("High Resolution (%.0f ppm)", HIGH_RES_PPM)); %[output:0cfb32e5]
xlabel("Unknown Sample"); ylabel("Number of Database Candidates"); %[output:0cfb32e5]
title("Candidate Reduction: Low vs High Resolution MS"); %[output:0cfb32e5]
grid on; %[output:0cfb32e5]
maxCand = ceil(max(candidateCounts(:)));
yticks(0:maxCand + 1); %[output:0cfb32e5]
ylim([0, maxCand + 1.5]); %[output:0cfb32e5]
%%
%[text] ## Section 5: Isotope Pattern Scoring
%[text] With the candidate list narrowed by exact mass, we now apply the final discriminator: isotope pattern scoring.
%[text] ### Why Isotope Clusters are "Structural Fingerprints"
%[text] All elements have multiple stable isotopes in nature, and their abundance ratios are unique to each element.
%[text] Different elemental compositions in a molecule change the relative intensities of $M$, $M+1$, and $M+2$, functioning as a "fingerprint" of the elemental composition.
%[text] 
%[text] The $M+1$ relative intensity (with $M$ normalized to 100%) follows from the approximation:
%[text] $M+1\\% \\approx 1.10 \\times n(C) + 0.37 \\times n(N) + 0.015 \\times n(H) + 0.04 \\times n(O)$
%[text] where $n(C)$, $n(N)$, $n(H)$, and $n(O)$ denote the number of C, N, H, and O atoms, respectively.
%[text] The largest contribution is from $^{13}$C (natural abundance 1.103%), so the more carbon atoms, the larger the $M+1$ peak.
%[text] 
%[text] The $M+2$ relative intensity arises from two-heavy-isotope combinations:
%[text] $M+2\\% \\approx \\frac{(1.10 \\times n(C))^2}{200} + 0.21 \\times n(O) + 4.25 \\times n(S) + 32.7 \\times n(Cl)$
%[text] Molecules containing sulfur ($^{34}$S: natural abundance 4.25%) or chlorine ($^{37}$Cl: 24.23%) have an unusually large $M+2$.
%[text] Especially for molecules with one chlorine atom, $M : M+2 \\approx 3 : 1$, making it easily distinguishable in the spectrum.
%[text] 
%[text] ### Scoring by Cosine Similarity
%[text] We score each candidate by computing cosine similarity between the observed isotope vector $\[I(M),\\ I(M+1),\\ I(M+2)\]$ (normalized to 1 at the $M$ peak) and the theoretical vector.
%[text] Cosine similarity is $\\frac{a \\cdot b}{|a| \\times |b|}$, where 1 indicates a perfect match and 0 indicates completely orthogonal patterns.
logSection("A08", "Section 5: Isotope Pattern Scoring", "Analytics L3"); %[output:6a8239e9]
isoScoreTbls = cell(1, N_UNKNOWNS);

for k = 1:N_UNKNOWNS
    % Use calibration mass for candidate selection and isotope window targeting
    mNeutral   = mH_hires(k) - PROTON_MASS;
    ppmDiff    = abs(refTbl.ExactMass - mNeutral) ./ refTbl.ExactMass * 1e6;
    candIdx    = find(ppmDiff <= HIGH_RES_PPM);
    if isempty(candIdx), continue; end

    % Observed isotope intensities (M, M+1, M+2) from smoothed spectrum
    sp     = spectra{k};
    smInt  = smoothdata(sp.intensity, "gaussian", 5);
    smInt  = smInt / max(smInt);   % normalise to 1

    obsIso = zeros(1, 3);
    for ishift = 0:2
        targetMz       = mH_hires(k) + ishift;
        [~, bestPt]    = min(abs(sp.mz - targetMz));
        obsIso(ishift+1) = max(smInt(max(1,bestPt-3):min(end,bestPt+3)));
    end
    if obsIso(1) < 1e-6, obsIso(1) = 1; end  % Zero division guard
    obsIso = obsIso / obsIso(1);              % Normalize relative to M peak

    % Scoring each candidate
    nCand  = numel(candIdx);
    scores = nan(nCand, 1);
    for ci = 1:nCand
        fc        = parseFormula(refTbl.Formula(candIdx(ci)));
        theoIso   = isoPattern(fc);   % [1, M+1/M, M+2/M]
        % Cosine similarity
        scores(ci) = dot(obsIso, theoIso) / ...
                     (norm(obsIso) * max(norm(theoIso), 1e-12));
    end

    [scoresSorted, sortOrd] = sort(scores, "descend");
    isoScoreTbls{k} = table( ...
        refTbl.Name(candIdx(sortOrd)), ...
        refTbl.Formula(candIdx(sortOrd)), ...
        refTbl.ExactMass(candIdx(sortOrd)), ...
        ppmDiff(candIdx(sortOrd)), ...
        scoresSorted, ...
        VariableNames=["Name","Formula","ExactMass","ppmError","IsoScore"]);
end
%[text] Display the isotope score ranking for Unknown 1.
logInfo("Unknown 1 -- High-resolution candidates ranked by isotope score:"); %[output:75c7f4a6]
if ~isempty(isoScoreTbls{1}) %[output:group:0de5348d]
    disp(isoScoreTbls{1}(1:min(5,height(isoScoreTbls{1})), :)); %[output:701fdc6c]
end %[output:group:0de5348d]
%%
%[text] ## Section 6: Complete Identification Workflow
%[text] With isotope scores computed, we now integrate all three stages into a complete end-to-end identification workflow.
%[text] ### Integration of 3-Stage Workflow
%[text] The five-step workflow proceeds as follows:
%[text:table]{"ignoreHeader":true}
%[text] | Step | Description |
%[text] | --- | --- |
%[text] | 1 | Extract $\[M+H\]^+$ from raw spectrum (`findpeaks`) |
%[text] | 2 | Back-calculate neutral exact mass (subtract proton 1.00728 Da) |
%[text] | 3 | High-resolution mass search (5 ppm window) |
%[text] | 4 | Rank remaining candidates by isotope score (cosine similarity) |
%[text] | 5 | Select top candidate as identification |
%[text:table]
%[text] 
%[text] **Relation to S05 (Unknown Substance Identification)**
%[text] In S05, unknown compounds were identified using fingerprint similarity (Tanimoto) from a structural information database.
%[text] In A08, the same goal is achieved using physical signals (mass spectra) without prior structural knowledge.
%[text] In practice, MS provides exact mass and isotopic composition, while FP similarity confirms matches when SMILES are available, forming a complementary relationship.
logSection("A08", "Section 6: Complete Identification Workflow", "Analytics L3"); %[output:4e58e823]
logInfo("Executing complete identification workflow..."); %[output:71c00218]
nCorrect = 0;
resultRows = cell(N_UNKNOWNS, 1);

for k = 1:N_UNKNOWNS %[output:group:5a0f4cf5]
    if isempty(isoScoreTbls{k})
        logWarn("Unknown substance %d: Insufficient data for identification", k);
        continue;
    end

    predicted = isoScoreTbls{k}.Name(1);   % top-ranked candidate
    trueID    = unknownNames(k);
    correct   = strcmpi(predicted, trueID);
    if correct, nCorrect = nCorrect + 1; end

    resultRows{k} = {k, trueID, predicted, ...
                     isoScoreTbls{k}.ppmError(1), ...
                     isoScoreTbls{k}.IsoScore(1), correct};

    resultLabel = " Incorrect"; if correct, resultLabel = " Correct"; end
    logInfo("Unknown substance %d:  True=%s  Predicted=%s  [%s]", ... %[output:3f094f67]
        k, trueID, predicted, resultLabel); %[output:3f094f67]
end %[output:group:5a0f4cf5]
%[text] Display summary table.
validRows = find(~cellfun(@isempty, resultRows));
if ~isempty(validRows) %[output:group:18193514]
    resultTbl = cell2table(vertcat(resultRows{validRows}), ...
        VariableNames=["Unknown","TrueID","Predicted","ppmError","IsoScore","Correct"]);
    disp(resultTbl); %[output:7d7b5589]
    logInfo("Identification accuracy: %d / %d = %.0f%%", ... %[output:686bed15]
        nCorrect, numel(validRows), nCorrect/numel(validRows)*100); %[output:686bed15]
end %[output:group:18193514]
%[text] ## Summary
%[text] - Exact monoisotopic mass (ExactMolWt) forms the basis for MS database searches.
%[text] - High-resolution MS (below 5 ppm) narrows candidates from dozens to 1-2.
%[text] - Molecules containing Cl or S have abnormally large $M+2$ peaks, serving as an elemental composition "fingerprint."
%[text] - A two-stage filter of mass + isotopic cosine similarity improves identification accuracy.
%[text] -  \
%%
%[text] ## Overview of Local Functions
function spectrum = simulateSpectrum(peakMz, peakInt, mzRange, sigma, noiseLevel)
% SIMULATESPECTRUM  Generates a synthetic mass spectrum with Gaussian peaks.
%
%   spectrum = simulateSpectrum(peakMz, peakInt, mzRange, sigma, noiseLevel)
%
%   Inputs:
%     peakMz     -- 1xN vector of peak center m/z values (Da)
%     peakInt    -- 1xN vector of relative peak intensities (0 to 1)
%     mzRange    -- Range of the spectrum axis [mzMin, mzMax]
%     sigma      -- Gaussian width (Da); FWHM = 2.355*sigma
%                   sigma = 0.3 Da --> FWHM = 0.71 Da (simulates low-resolution ESI)
%     noiseLevel -- Uniform noise amplitude as a fraction of the base peak
%
%   Output:
%     spectrum   -- struct with fields .mz (1x5000) and .intensity (1x5000)
%
%   About intensity range:
%     Gaussian components are normalized to [0, 1] before noise addition.
%     After uniform noise addition (noiseLevel = 0.03), the maximum intensity can reach up to 1.03.
%     This is physically realistic (random noise rides over peak tops) and intentional.
%     Callers needing strict [0, 1] range should divide by max(spectrum.intensity) later.
    mz        = linspace(mzRange(1), mzRange(2), 5000);
    intensity = zeros(1, 5000);
    for i = 1:numel(peakMz)
        intensity = intensity + peakInt(i) * ...
            exp(-0.5 * ((mz - peakMz(i)) / sigma).^2);
    end
    if max(intensity) > 0
        intensity = intensity / max(intensity);
    end
    intensity = intensity + noiseLevel * rand(1, 5000);
    spectrum  = struct("mz", mz, "intensity", intensity);
end

function counts = parseFormula(formulaStr)
% PARSEFORMULA  Extracts element counts from a molecular formula string.
%
%   counts = parseFormula("C20H25N3O2")
%   Returns a struct with fields C, H, N, O, S, Cl, F, Br, P, I (all double).
%   Non-existent elements are 0.
%
%   Regular expression: ([A-Z][a-z]?)(\d*)
%     [A-Z][a-z]?  Matches element symbols: 1 uppercase + optional lowercase letter
%                  --> Correctly captures 2-letter symbols: Cl, Br, etc.
%                  --> "Cl" matches as token {"Cl",""}, not {"C",""} + {"l",""}
%     (\d*)        Matches count digits; empty string --> treated as 1
%
%   Verified examples (Hill notation returned by RDKit rdMolDescriptors):
%     "CH4"       --> C:1, H:4
%     "C9H8O4"    --> C:9, H:8, O:4 (Aspirin)
%     "C6H5Cl"    --> C:6, H:5, Cl:1 (Chlorobenzene -- 2-letter element OK)
%     "CCl4"      --> C:1, Cl:4 (Carbon tetrachloride)
%     "H2O"       --> H:2, O:1 (No number for last element -- guard OK)
%
%   Edge cases:
%     Empty string ""    --> All counts 0 (tokens = {}; does not enter loop)
%     Unknown elements   --> Implicitly ignored by isfield guard
%     (e.g., "Fe", "Mg" from inorganic impurities)
    counts = struct("C",0,"H",0,"N",0,"O",0,"S",0,"Cl",0,"F",0,"Br",0,"P",0,"I",0);
    tokens = regexp(char(formulaStr), "([A-Z][a-z]?)(\d*)", "tokens");
    for i = 1:numel(tokens)
        elem = tokens{i}{1};
        n    = tokens{i}{2};
        if isempty(n), n = "1"; end
        if isfield(counts, elem)
            counts.(elem) = counts.(elem) + str2double(n);
        end
    end
end

function rel = isoPattern(fc)
% ISOPATTERN  Calculates relative isotopic intensities for M, M+1, M+2.
%
%   rel = isoPattern(fc)
%   fc  -- Element count struct from parseFormula()
%   rel -- [1, M+1/M, M+2/M] (normalized to M = 1.0)
%
%   Approximation (Gross 2011, Chapter 3; valid for <100 carbon atoms):
%     M+1% = 1.103*nC + 0.366*nN + 0.015*nH + 0.038*nO
%     M+2% = (1.103*nC)^2/200 + 0.205*nO + 4.25*nS + 32.7*nCl
%
%   Coefficient derivation (natural abundance: Meija et al. 2016 IUPAC):
%     13C: 1.103% natural abundance --> M+1 contribution = 1.103% per C
%     34S: 4.25%  natural abundance --> M+2 contribution = 4.25% per S
%     37Cl: 24.23% natural abundance
%           M+2 contribution = p(37Cl)/p(35Cl) * 100
%                    = 24.23/75.77 * 100 = 31.98% per Cl
%           Note: 31.98% is exact; 32.7% used here per Gross (2011) Table 3.2.
%           Difference < 2.2%. Reference "37Cl 24.23%" is raw natural abundance;
%           M+2 coefficient is not abundance but 37Cl/35Cl ratio.
%
%   Cross-validated reference values (validated in Claesen et al. 2012):
%     CH4 (C:1, H:4):            rel ~ [1.000, 0.01163, 0.00006]
%       M+1% = 1.103 + 0.060 = 1.163%;  M+2% = 0.006% (barely detectable)
%
%     C9H8O4 Aspirin (C:9, H:8, O:4): rel ~ [1.000, 0.10199, 0.01313]
%       M+1% = 9.927 + 0.120 + 0.152 = 10.199%
%       M+2% = (9.927)^2/200 + 0.205*4 = 0.493 + 0.820 = 1.313%
%
%     C6H5Cl Chlorobenzene (C:6, H:5, Cl:1): rel ~ [1.000, 0.06693, 0.32919]
%       M+1% = 6.618 + 0.075 = 6.693%
%       M+2% = (6.618)^2/200 + 32.7 = 0.219 + 32.7 = 32.919%
%       --> M:M+2 ratio = 1:0.329 confirms "Cl 3:1 rule" (Gross 2011 p.95)
    m1pct = 1.103*fc.C + 0.366*fc.N + 0.015*fc.H + 0.038*fc.O;
    m2pct = (1.103*fc.C)^2 / 200 + 0.205*fc.O + 4.25*fc.S + 32.7*fc.Cl;

    rel = [1.0, m1pct/100, m2pct/100];
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:930256b9]
%   data: {"dataType":"text","outputData":{"text":"[21:18:11][INFO]  initPython: Python 3.10 already active (Status: Loaded) -- skipping.\n","truncated":false}}
%---
%[output:054dd858]
%   data: {"dataType":"text","outputData":{"text":"[21:18:11][INFO]  --- A08 | Section 0: Setup  [Analytics L3] ---\n","truncated":false}}
%---
%[output:7f13b345]
%   data: {"dataType":"text","outputData":{"text":"[21:18:11][INFO]  --- A08 | Section 1: Building a Precise Mass Reference Table  [Analytics L3] ---\n","truncated":false}}
%---
%[output:337aef13]
%   data: {"dataType":"text","outputData":{"text":"[21:18:11][INFO]  Read 200 rows from data\/list\/fda_drugs.csv\n","truncated":false}}
%---
%[output:75d7087d]
%   data: {"dataType":"text","outputData":{"text":"[21:18:11][INFO]  Calculating exact mass (may take 1-2 minutes)...\n","truncated":false}}
%---
%[output:85be30ce]
%   data: {"dataType":"text","outputData":{"text":"[21:18:51][INFO]  Reference table: 200 compounds\n","truncated":false}}
%---
%[output:1c8b109d]
%   data: {"dataType":"text","outputData":{"text":"          Name             Formula       ExactMass\n    ________________    _____________    _________\n\n    \"PRAZOSIN\"          \"C19H21N5O4\"      383.16  \n    \"NICOTINE\"          \"C10H14N2\"        162.12  \n    \"OFLOXACIN\"         \"C18H20FN3O4\"     361.14  \n    \"NALIDIXIC ACID\"    \"C12H12N2O3\"      232.08  \n    \"INDOMETHACIN\"      \"C19H16ClNO4\"     357.08  \n\n","truncated":false}}
%---
%[output:44614b03]
%   data: {"dataType":"text","outputData":{"text":"[21:18:51][INFO]  --- A08 | Section 2: Simulation of 5 Unknown ESI-MS Spectra  [Analytics L3] ---\n","truncated":false}}
%---
%[output:46addd2a]
%   data: {"dataType":"text","outputData":{"text":"[21:18:51][INFO]  Selected unknown substances:\n","truncated":false}}
%---
%[output:0103aab5]
%   data: {"dataType":"text","outputData":{"text":"[21:18:51][INFO]    Unknown substance 1: AMPHETAMINE  (135.1048 Da, molecular formula C9H13N)\n[21:18:51][INFO]    Unknown substance 2: FUROSEMIDE  (330.0077 Da, molecular formula C12H11ClN2O5S)\n[21:18:51][INFO]    Unknown substance 3: MORPHINE  (285.1365 Da, molecular formula C17H19NO3)\n[21:18:51][INFO]    Unknown substance 4: CYPROHEPTADINE  (287.1674 Da, molecular formula C21H21N)\n[21:18:51][INFO]    Unknown substance 5: CAFFEINE  (194.0804 Da, molecular formula C8H10N4O2)\n","truncated":false}}
%---
%[output:6a938815]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAewAAAEpCAYAAABRDxFnAAAQAElEQVR4Aeydz6tc2XXvl42MrsAPugdCFq9Dnwj\/D9IgXFcg+QveyCK4y4OMNeppXAaNnkZvHEGOCWgQcAaZhOBBdpOAY8g8o4jV6IZWcCAGG9KhBS\/9qdNLtevcc+pW3ao6P7\/C6+z1a++99uece5aqdG1\/+\/\/rjwiIgAiIgAiIwOAJfNv0RwREQAREQAREYPAEht2wB49PBYqACIiACIhANwTUsLvhrF1EQAREQARE4CgCati3x6eZIiACIiACItAZATXszlBrIxEQAREQARG4PQE17NuzG\/ZMVScCIiACIjApAmrYk7qdOowIiIAIiMBUCahhT\/XODvtcqk4ERGAEBB49ejSCKqddYn4P1LCnfa91OhEQARG4FQEaxevXr281V5NOR4B7wL1gxZM2bBbNhQ1yIZbbp9TPuXZbnfvsSU6b1Net5+XxPJb7d+nM2RWPWD0PO5fI22fM56HvM2dwOSMoaBfbXbH60Q7Jrc8dgk39uRxaE3MPnTPm\/GPPe+x82B27xjHzd83NY7lOzbnsipF3U5yc28rJGjZF8jeBXPDdtrApzcuZ5Hp+RljlMXR8kYONHiP6PpKv0ZRfj2OzRy74mubWfeTl89Dx1fNki8ApCPBs8Yzlgu8Ua09xjUPYkAvXKXIY45m4F9yTkzRsFmLBOgh8xMKPHXrfY17XEGppYoOvyzrZiz3rPPARq\/tzmzh5uQ8dHzH0rqTr\/bo61zf7aPiaAPeYZ+trdes\/+IhtOc9odLnXMccYS53HnPGUc3mOTrneqdY6ScPeVcxQD76r5iHFjuXH\/LYfVvzE9znvvnlNax0zt2k9+URgFwE9b9fpwAS5HpFnTARO0rB5EHj5I7sOn8dDZ0RiHjoSdoz7+sgnNxd8IeHHznVsJHyM2LngC8n9x+r78jt2n5vmH1PHvnPhRx2MIdi5hJ8x94eOP5e6Hzvi6Ag2I4KOhM6YS8Twhc6I4EPQEfQhS9TIGLKrXnIiHjpjSMRiDD\/j2vfNpclu8pEefsYQ\/Lvk0OctX4s9chsdXwh2SPhirPux22J1P7kI\/hjREexc8OWSx9B3xYj3KXlt6G21EAup54SfsR7LbeK55LFjdNaM+egh4cvHiDHm\/lwnFpL70fHHiI5gN8lJGjYL8wOEsFkI\/l1CHnMQdAQdQd81ty3GPObngi\/yw4+d69jkhY8RGz+Cji8EG\/8+Qm5d6vPydSO3npPbkdM05nmH6ofWka+\/71xqjlxG7FgHHV8IdsQYsSMWIz5iYdd1bIQ8JM\/Df5Pkc9CRWAN913zibcK8thh+4qcQ1op6GbGb1sVPPI+FDz+CHXF0fCHYETt0ZG6sw4h90xrkIeSG3DSnKc5c1gnBJo8xfDHiIxZ2XScesRjxkZcLvogzYkccHV8u+NrieSxy+hqpJa8bHV+9HnzEQrAjBz38jNgRy0f8xHPBl+c06eQ0SVvurvVZZ1ecNU+Vw1ona9gshtxUPDkh5IbOWLfxHSq3XSOg5vuxFv7cFzqx0G8aya1L25w8r21v5uZ5dZ14LsTra2Hjz\/NynVgIuXnsJj3mMTbNxd+0Brn1GDZ+8hmx0XNp8uXx0MlDwt53rM+p27vWIbdNmNcWw0\/8FLLPWoeybcpnH\/zUnOvYdSGPnPDnevj2HZkbwro75jWGmNsYaHDelHtTPJZsy6P+plj4muLE8MfafY7Uss\/+bXmcox7Dxn\/Muvlc1muSPKdNZ15bDH89Tt11HzZ+8kPwhb5rPEnDrm8eG1JEWyxyzjGyZy777pHPCX3fucfksVfT\/K75HVPHMXPzs7NOXfK49NMTgPc5n7VY\/5R7sGYTiVPuEWuxV0jTnnVf5MZYjx9rx7r5eOyap5yf14V+6NrMqUvTGre9P01rndNXPwv2bfc7ScO+7ebnmAcMbmQu++6Tz8n1fecPOY\/zwIYaGbHRhybU1SRDq3NK9cB7Suc55Vlgkws\/O7vWJ57no+\/K34rtabBmk+w5fa801ucseyVnScxhbi5ZeC81n5vrTZPzODr7N+X16aOuJjmkJs7FGidp2CzEgvUC8BGr+09hs\/Y+6+ybR51NuU0+9m3zEztUdu1N7ND1bpvPXk3nwkds17rEyavn4COW+\/HldujkNcXCd1M81jlmjL2OWaOLubtYELtNDczb9\/xNuczFH3uj7+OL\/EPGWLs+p75fU3wfHzmsxXiMHLrGrnNRx01xcoYkbedv8x9yvrY1znn+m\/asxw85zz51n6Rhs1EURsEh+IidQlgr1mXEbloXP\/GQphx8eR42kvtiPr6mWPiJ3SSxVn3M57FePY4vcoihx4h+qLAe8xnb5hIjJxd8bfm5n7x8Hjq+PAcdH7EQbPwIevhjxEcMQQ9\/jPiIhWBHLHxtY57LHOy23KH5qZWac8F3TJ3MZ7191ohc8hHsfeadKof92DcXfPn62Lvi5Lbl1P2sg485Idj4EXy5jQ8b\/yHCHObmgi\/WQM9j6PgifqqRNVn7m\/XeD\/iahISYE3Fs\/HXBHzmM2JGDji8XfBGPEV+eg44v4qcYWY91Q7DzdbEjxoidx9HxEcsFH7F9hHmRf7KGzcYsmgu+XIiFnev4brLJCYncGPHnetj4EOwmIYbkMexc2mL4yWPcJeS0SX1ePS+P57Hcv0tnTj2+j4+cXOpr7LLzeehtucRC6jnhj\/HQOPn1udj4m4RYCHF0RiTX97HJ6VKoL5f63sR2+XbFd8ViTXJCwpePxHIbve6r2005+JqEubnclEOcfEYkdMYQ\/CHhizH8+ViPhc1IXozoSN1u8pGTCzm55DH0PLZLPySXdcinYaAj2G1CHMnjYTOGEEdnDMHOJfwx1mNhRzzG8LeN5O0Ty\/PQQ5iLzhiCHYIPnTEXfLnUY7mNTi4j7EPHPmnDZsFzCYXvszZ5+QH3maMcERABERCBZgKjeJ82lz56b539aBo2hSOjvwMzP4Du4cwfAB1fBETg1gRG07D3PWG9Iez6xH11dWUSMdAzoGdAz8B8noF9e8kJ80621OQadk7mpmb96aef2uXlpUQM9AzoGdAzMJNn4Ic\/\/OH6g1reK8aiT7Zh72rW3Bz+Rv2rX\/3KXrx4Ya9evZKIwWyegWfPnvEjoGdfz\/xsnvl4x\/Ps895f\/wCM8HKWht03h5uadV7f48eP7cmTJxIxmM0zwDPPzwCjnn397M\/pGeCZ59kfq0yuYR\/SrMd601S3CBxD4Hvf+5796Ec\/MsZj1tFcERCBbglMpmHTqAMdel0iZiZNBOZNgEb9ySefzBuCTi8CIyQwuoZd\/y3wYB5+xiaJPI0iIAIiIAIiMEYCo2vYY4R8SM3KFQEREAEREIEmAmrYTVTkEwEREAEREIGBEVDDHtgNGXY5qk4EREAERKAvAmrYfZHXviIgAiIgAiJwAAE17ANgKXXYBFSdCIiACEyZgBr2lO+uziYCIiACIjAZAmrYk7mVOsiwCag6ERABETiOgBr2cfw0WwREQAREQAQ6IaCG3QlmbSICwyag6kRABIZPQA17+PdIFYqACIiACIiAqWHrIRABERg4AZUnAiIAATVsKEhEQAREQAREYOAE1LAHfoNUngiIwLAJqDoR6IqAGnZXpLWPCIiACIiACBxBQA37CHiaKgIiIALDJqDqpkRADXtKd1NnEQEREAERmCwBNezJ3lodTAREQASGTUDVHUZADfswXsoWAREQAREQgV4IqGH3gl2bioAIiIAIDJvA8KpTwx7ePVFFIiACIiACInCNgBr2NSRyiIAIiIAIiMDwCOQNe3jVqSIREAEREAEREIE1ATXsNQZdREAEREAERGDYBMbTsIfNUdWJgAiIgAiIwFkJqGGfFa8WFwEREAEREIHTEFDDPg1HrSICIiACIiACZyWghn1WvFpcBERABERABE5DQA37NByHvYqqEwEREAERGD0BNezR30IdQAREQAROT8DdLKXTr6sVb09gsg370aNHt6eimV0S0F4iIAIDJPDTn5ohAyxttiVNsmGrWc\/2edbBRUAERGCyBCbXsNWsJ\/us9nMw7SoCMybgPuPDD\/Dok2vYr1+\/NmRf1ldXVxby7t07k4iBngE9A3oG3q1foe42+ndivN8Z14ca8WVyDfvQe\/H06VO7vLxcy\/Pnz+3NmzcSMRjLM3CrOnlxvX379lZz9fMxn\/fDb3\/7u\/Xr9Je\/\/GLUz8rLly\/X73fe87zv14ca6WX2DfvFixf26tWrtSyXS3v48KFEDCb9DDx48MDu378\/6TPq5\/j499i\/\/\/vFuq2N\/VnhvR7v+GfPnq3PNNbL7Bv248eP7cmTJ2spisIuLi4kYjDpZ+DevXt29+7d859Rz9GoGV9d3Vn3tbG\/E3mvxzue9\/36UCO9zL5hj\/S+qWwREAEREIGZEVDDntkN13FFQATWBHQRgdERUMMe3S1TwSIgAiIgAnMkMNmGfch\/tWuON15nFgERGDABlSYCDQQm27AbziqXCIiACIiACIyWgBr2aG+dChcBERCBXgho054IqGH3BF7bioAIiIAIiMAhBNSwD6GlXBEQARGYGQH3kR14wuWqYU\/45upoIiACInBbAu63nal55yKghn0uslpXBERABERABLYJHGWpYR+FT5NFQAREQAREoBsCatjdcNYuIiACIiACInAUgbM37KOq02QREAEREAEREIE1ATXsNQZdREAEREAERGDYBGbesId9c1SdCIiACPRNwL3vCrR\/EFDDDhIaRUAEREAERGDABNSwB3xzVJoIiIAIiIAIBAE17CChUQREQAREYE3AfT3oMjACatgDuyHjKUeVioAIiIAIdEmgl4b96NGjLs+ovURABERABERg9AR6adivX782mnbI6CnqAIMjoIJEQAREYGoEemnYQKRph6hxQ0QiAiIgAsMj4D68muZaUW8NOweuxp3TkD59AjqhCIiACBxOYBANOz5h54378KNohgiIgAiIgAhMl0BvDTuaNGM06uli1slEYDwEVKkIuIvBEAn00rDzJk2zroNp8tVzZIuACIiACIjAnAj00rDVkOf0iOmsInBKAlpLBOZLoJeGzSfsJuRt\/qZc+URABERABERgTgR6adhzAqyzioAIzIeATioC5yTQacPmEzTCgRjroq\/KISMRAREQAREQgesEOm3YNGSEMhjrgl8iAiIgAiJwDgKHr1kUh8\/RjPMR6LRhn+8YWlkEREAERODUBIri1CtqvWMI9NKw+WR926Lzr9Hb1shz0Nvy5BcBERABERgGAVVxM4FOGjZNE4ly0NskcppG5tDsQ7DrefgiHiO+ep5sERABERCBmwl8\/vnNOcrohkAnDTsaZxwp7KYxcuojTZf83I+NP\/dJFwEREAEREIHTEhjGap007GEcVVWIgAiIgAjsQ8B9nyzldE2gt4adfzJGR05x+PjUzXoh+NrWvrq6spB3796ZRAz0DOgZmPszkL8vv\/pq3M9DvN8Z83ONUW9r2Gc9C400mmjo2OjHe0dAHgAAEABJREFUbswarJULvrZ1nz59apeXl2t5\/vy5vXnzRiIGk34GeHG9fft20mfUz\/Fx77H\/+I9fr1+Z\/\/VfX9qXX3456mfl5cuX6\/c773ne9+uDjfTSS8MeEqsXL17Yq1ev1rJcLu3hw4cSMZj0M\/DgwQO7f\/\/+pM+on+Pj3mMffPDh+jX9\/e\/fsYuLi1E\/K7zX4x3\/7NkzG\/OfcTbsExJ\/\/PixPXnyZC1FUawfTh5QyYVYfP2imuJzcO\/ePbt7967u70Tv7yme2Tt37qzfst\/5zh1DTrFmX2vwXo93PO\/79cFGehlNw+Yr7vpX29j4R8peZYuACIiACIjA3gR6adg0WZotgk61uY7dJOSSF4IdefjQ8aHngo9YR6JtREAEREAERODkBHpp2JyCJoqgI7mO3SbkheQ5+MJGzyX8GkVABERABERgrAR6a9j5J+BcHyvI0dStQkVABETgAALuByQr9awEemnYNOj8E3Cun\/W0WlwEREAEREAERkqgl4Y9UlYq+\/wEtIMIiIAIiEALATXsFjByi4AIiMCcCRTFnE8\/zLP30rD5CpyvxYeJRFWJQAsBuUVgJgTcZ3LQkR2zl4YdzZqxLiPjp3JFQAREQAREoBMCvTRsPmG3SSen1iYiMD0COpEIiMDECfTSsCfOVMcTAREQAREQgZMT6LVhx9fhnAqdUSICIjBBAjqSCIjA0QR6a9g0aL4WjxOg4wtbowiIgAiIgAiIwIZALw2bxkyD3pQhTQREQAR6I6CNRWAUBHpp2KMgoyJFQAREYOYEPv545gAGdnw17IHdEJUjAiIgAlsEZIjANwR6adh8Hc7X4t\/UsB6w8a8NXURABERABERABLYI9NKwqYDmTJNGZ8RGl4iACIiACAyHgPvOWhTskEBvDZsz0qRDsCUiIAIiIAIiIALNBHpp2Hyibiqnzd+UK58IiIAIiIAI7CQwsWAvDXtiDHUcERABEZgcgaKY3JFGf6BOGzafoBGoMdaFr8eJSURABERABERg4gQOPl6nDZuGjFAlY13wS0RABERABPol8Pnn\/e6v3ZsJdNqwm0uQVwREQAREQARE4CYCnTbsKIZP1qFrFAEREAERGCaBohhmXXOtqpeGDez6v1+HTUwiAiIgAiIgAiKwTaCXhk1z5lN2k2yX16WlvURABERABERguAR6adjDxaHKREAEREAERGCYBNSwh3lfrlUlhwiIgAj0QcC9j121ZxOBXho2X4XztXhTQfKJgAiIgAiIgAhcJ9BLw45mzViX6yXKM3wCqlAEREAERODcBHpp2HzCbpObDpw3+F25++btWkMxERABERABERgKgV4a9m0PTxPOGz1201r498lrmivf9AjoRCIgAocTKIrD52jGeQl00rBpoPscY1ceMZpwvg42\/tyHjT\/31e08Jl0EREAEREAExkCgk4YNCBrpTULeqSTf61Rrah0ROD0BrSgCIiAC+xHorGHzKfcm2a\/km7No1vle2G2zrq6uLOTdu3cmEQM9A3oG5v4M8L786qt3VhRoNur3YrzfGavTjPfaScOmee6DaN+8m9Y6ZJ2nT5\/a5eXlWp4\/f25v3ryRiMGknwFeXG\/fvt3rjPp5mOf74F\/+5Uv78ssv7Ysvvli\/bsf8HLx8+XL9fuc9z\/t+faCRXjpp2ENm8+LFC3v16tValsulPXz4UCIGk34GHjx4YPfv35\/0GfVzfNx77DvfuWMXFxf24YcfGn\/GzJP3erzjnz17xnFGK7Nv2I8fP7YnT56spSiK9UPKgyq5EIuvX1hTfA7u3btnd+\/encD91TN6rufzzp07Fk2b7naufbpYl\/d6vON533OescpoGjZfc9f\/LRobfw4fG3\/uky4CIiACIiACYycwmoYN6GjGNGQEGz+CzYjgxw7Bxi8RAREQga4IaB8RODWBUTVsDk\/zDcEOwRc6I3YItkQEREAEREAExkyg14Ydn4ABiM4oEQEREAERODeBw9Z3Pyxf2ech0FvDpkHzCTiOhY4vbI0iIAIiIAIiIAIbAr00bBozDXpThjQREAEREAERqAjo2kygl4bdXIq8IiACIiACIiACbQTUsNvIyC8CIiACMybw8cczPvzOo\/cX7KVh83U4X4vnx8bGn\/uki4AIiIAIiIAIVAR6adhsTXOmSaMzYqNLREAEREAEREAErhPYp2Ffn3UiD0065ERLahkREAEREIEjCbgfuYCmn4VALw2bT9RnOY0WFQEREAEREIGJEuilYfOpmqYdchRbTRYBERABETgLgaI4y7Ja9JYEemnY1ErTDlHjhohEBERABERABNoJ9Naw85Lyxp37J6DrCCIgAiIwegLuoz\/CJA4wiIYdn7Bp3JOgqkOIgAiIgAiIwIkJ9Nawo0kz0qiRE59Ny91EQHEREAERaCFQFC0BuXsj0EvDzpu0GnVv914bi4AIiIAIjIhALw1bTXpET0h\/pWpnERABERCBjEAnDZtP1Ejsi94mkaNRBERABERABERgQ6CThs0naiS2RW+TyNEoAoMmoOJEQAREoGMCnTTsjs+k7URABERABI4g4F5NLopq1HUYBHpp2Hwd3nT8Nn9TrnwiIAKtBBQQARGYIIFeGvYEOepIIiACIiACInBWAp02bD5BI5yIsS78uzYxiQiIwIQJ6GgiIAK3ItBpw6YhI1TKWBf8EhEQAREQgWERcB9WPXOtptOGPVfIOrcIiMBoCKhQERgsgV4aNp+sB0tEhYmACIiACIjAAAn00rAHyEEliYAIiMDwCajCWRPorWHXf+Es7FnfDR1eBERABAZCoCgGUojKeE+gl4ZNc+ZrcYRK6iM+iQiIgAiIwKgIqNgzE+ilYTediaZNI2+KyScCIiACIiACcycwmIa9742gqYfcNGffvJvWUVwEREAERGDkBCZQ\/qgaNg2YT+Ih2BO4BzqCCIiACAyGgPtgSlEhNQK9NGwabjTb0LHRa\/W9N5vi5ON\/n5Qp+IlnLqkiIAIiIAIHEiiKAyco\/TYE9prTS8OmsryZoiP4u5arqysLeffunUnEQM+AnoE5PwPxDg4G2KGPcYz3OyNnGbN00rD5tLuvnAIme+37F4CnT5\/a5eXlWp4\/f25v3ryRiMGknwFeXG\/fvp30GfVzfNx7jPfwr3\/96\/UzQpP+zW\/+c62PkevLly\/X73fe87zvOdtYpZOGTfOsS5t9LMhDmjV7vXjxwl69erWW5XJpDx8+lIjBpJ+BBw8e2P379yd9Rv0c3\/49xrPBu5ERjnfu3LEPPvhwtM8L7\/V4xz979oyjjVY6adhd06Fph7A3OmOTPH782J48ebKWoijs4uJCIgaTfgbu3btnd+\/enfQZ9XN83HuMdyWNGo65jj024b0e73je95xnrNJrw6aRIsCLEf0YqX9yZy18jPuLMkVABERABIKAe2ga+yTQW8OmQeeNFB1fG4ymOPn42+bILwIiIAIiIAJTIdBLw75to6U5MzcEO24EvtCnPup8IiACpyWQktm3vmXmbvrzDYGv\/4XwG03DUAj00rCPOTxNOiRfB19uh97mj7hGERABERABERgDgdE17DFAnXeNOr0IjJeA+3hrV+XTJ9BLw+ZTb\/0rbGz800euE4qACAyNgHtVkXs16rohoK\/GNyz61npp2Bya5kyTRmfERpeIwDkJaG0RaCLgXnk\/+6wadd0m8Pnn27asfgj01rA5Lk06BFsiAiIgAiLQLwH3fvfX7u0EOm3YfJJG2srZFWubI78ITIeATiICIiAC7QQ6a9g04\/g0jR4loYcQD79GERABERABERCBDYFOGjYNOW\/G6PgQ9JBNWdJEQASGRmDK9ejfaKd8d6dztk4adhsuGnVbTH4REAER6JqAe9c7aj8R2J9Arw17\/zKVKQIiIAK7CCgmAtMnoIY9\/XusE4qACNxAwP2GhBmG9d+\/Ht5NV8Me3j1RRSIgAhMjoOOIwCkIdNaw+QWzXCg+t0PHLxEBERCBPgi497Gr9hSB\/Qh00rD55bJ9Zb+ylSUCIiACInAaArtX0Vfju\/l0Ge2kYXd5IO0lAiIgAiIgAlMkoIY9xbuqM4mACBxEwN1MnyQrZO7VOJSr6tgQUMPesJAmAiIwcwLuMweg4w+agBr2oG+PihMBERCB\/gm491\/DMCvotio17G55azcREIGBEtBX4gO9MSrrPQE17PcopIiACMyRgLuZu5ka9hzv\/rjOfGjDHtfpVK0IiIAIHEjA\/cAJSheBjgioYXcEWtuIgAiIgAiIwDEEptWwjyGhuSIgArMm8PHHsz7+1uH1zwNbOAZjqGEP5laoEBEQAREYHgH9RWY490QNu7t7oZ1EQAQGTECfKqub416N+dU9t6T3RUANuy\/y2lcERGBQBNSwB3U7VEwDATXsBiizdOnQIjBTAu7bB3fftudu8RcZ97lTGMb51bCHcR9UhQiIgAiIgAjsJKCGvROPggMhoDJE4OwE+CR59k1GsoFYDPNGja5hP3r0yEJ2IY0cxl15iomACIhATsA9t6SLwHAIjKph03xfv35tIdhNKPFHDiN2U558InASAlpEBERABDogMJqGTdOl+eZMsPHnPmz8uU+6CIiACLQRcK8ifA2MVNZ8r59\/Pt+zD\/3ko2nYQwep+kRgoARU1h4E8kbtvscEpYhADwQm17Drn65v+sR9dXVlIe\/evTOJGOgZmNcz4G7vf+55B8\/9\/tcZxF9mxsol3u+MnG3MMrmGnd+Mm5o1uU+fPrXLy8u1PH\/+3N68eSMRg0k\/A7y43r59O4wzDuBZ+81v\/pNXwZoHTel3v\/vtWp\/ru+C3v\/2dffXVu\/cMvvjii\/d8xsjk5cuX6\/c773ne9+vDjPQy2Ya9T7Pmnr148cJevXq1luVyaQ8fPpSIwaSfgQcPHtj9+\/cnfcZDfo6\/+93\/ZXfu3FnzYMQ+ZP7Uci8uLuzevYs1D8724Ycf8qq0\/\/7vcb4bea\/HO\/7Zs2frs4z1MsmGvW+z5qY9fvzYnjx5spbi6+9+eFglFyYG02Vw7949u3v3ru7x142J55wm\/fWP\/poH74R\/+7c7a\/3i4mKW43e+c8eCB3wQ+\/oP4xiF93q843nff32U0f5nNA2bf5umEeeksfHf5Mvj0kVABERABNoJuG\/HaN7bHll9ERhNwwYQzZkmHYKNH8HHiKDXBb9EBIZEwN3sxz8eUkXzrGW0\/zWmjm+Xe8cbartrBEbVsKmeJh2CHYIPnbFJiElEYEgEUjIrSzP3IVU171r0aXLe93\/opx9dwx46UNUnAvsScK8yy7Iade2fgBr2ye6BFjoDATXsM0DVkiIgAiIwZgIff7ypXn+J2bDoW1PD7vsOaP\/ZEtC\/nQ7n1udNyX04dQ2pEvchVXNkLSOdroY90hunssdPwL06gxp3xaGvq3tfOw9zX\/fmutyb\/fJ2R0ANuzvW2kkEGgm4N7rlFAERmCeB1lOrYbeiUUAEzkvAvVrfvRp17Y9A\/m+27v3VMdSd+ScDfRPU\/9FwAzYAABAASURBVN1Rw+7\/HqiCmRJwrw7uXo269kPAfbNv3rg3XmkQKEuzskST9EVgEA27r8NrXxHomwCfXPquQfubxX2IUUyaCfzsZ81+ebshoIbdDWftIgJbBFKqTDWIikOfV\/c+dx\/H3p98UtXpXo269kNADftG7koQgfMRUMM+H9vbrux+25nTmdf2XLpP54xjPIka9hjvmmqeDIEf\/GAyRxn1Qdoa1KgPdcLiF4vNYu4bXVq3BNSwu+V98t204DgJuI+z7qlV7T61E53nPIuF2WJRre1ejbp2T0ANu3vm2lEE1gSKwqwo1qr+D0AqDLoOgIB7cxE\/+Unld69GXbsnoIbdPfMZ7aijthH47LMqslhUY0rVqGu\/BIqi2t+9GnW9TsD9uk+ebgioYXfDWbuIQCuBomgNKdARgaLY3sh9256L5d5+0sWiPaZINwTUsLvhrF0GSEAliYC7GDQRKIomr3x9E1DD7vsOaP\/ZEiiKzdHdN7o0EeibQFE0V1AUZvHPOc0Z8p6TgBr2OelqbRHYQaAoNsHr\/zvNm5i08xMoimqPxaIadW0mUBTNfnm7IaCG3Q1n7SICWwTcN+ZisdGldUvAvXk\/92b\/1L3uUz\/huM+nhj3u+6fqRaAXAlPatCimdJrznqUozFI67x5avZ2AGnY7G0VE4GwE3M30\/wx1Nrx7L+x+PbUozNyv++XZPLPuotEHATXsPqhrz1kTSMnMfYOAxu2+saUdS+D4+fqdguMZaoXTE1DDPj1TrSgCBxEoCjP3g6Yo+UQE+I3notherCi27TlZKVWnLYpqrF+LovKkVI26dktADbtb3tpNBN4TWCzeq2rYGxSda0WxvWVRmLlv+05pjXmtxWLM1Y+\/djXs8d9DnWBkBFJqLjilZr+85yPgfr61x7gy\/xRQFGOsfB41q2HP4z7rlAMiwEuxqRz3Jq985yTgvvlFqnwf99yaj+5uVhQ3n\/dnP7s5RxmnJ6CGfXqmWlEEdhJwNysKs8WiSiuKanSvRl27IZBSN\/uMaZeUdldbFGZFYZaS2Y9\/bOa+O1\/R0xJQwz4tz\/eruZuVpZm7mbtZWZqtVlW4LM3czcrSLCWzsjQry+oHoCyrHK7uZimZpWRWlngqcTdLaaOvVhsdrSzNViszdzN3s7I0K0szd7OyNCtLsjaSkpl7ZbubuVd6WVYjV3czd7RK3M3K0qwszcrSLCUzd7OyNCtLs5TM3M3czcrSLCUzd7OyNHM3W62qdcrSrCzNUjIrSzP3yu9ulpKZe2WXpVlZmq1WZimZuZutVmYpmaVklpJZWZqlVOWvVmbuZimZrVZmZWnmbrZamaVklpJZWZqVpZm7WUpm7mZlaeZulpJZWVZruW9Gd7OUzMqy8sXV3Wy1MlutzFYrs7KMiJm7mbuZu1lZmqVkVhSb+GJR6T\/9qVlZmqVU2VxTMitLM3czd7PVyqwsiWzE3awszVIyczdLqYq5m5WlWVmauZu5m7lXsbiWpZm7mbtZSuE1W63MUjJLyawszdyrmLuZu1lKG3u1MitLM3ez1cqsLM1WqyqeX93N3M1SMnM3W63MViuzlMxSMkupynY3K0szdzN3s9XKzN2sLM3K0qwsN3kpbXT363pKZimZuZutVlWca1matX1KjN\/aL0syzcrSbLUyK8vKdjcrS7PVamOnZLZamZWlmXvlT8lstTJzN1utzMrSLCUz9yoe15TMyrKy3M3K0szdrCzN3M3czVYrs7I0czcrSzN3s5TMytLMfXtuWZqVZeUry2p0N0vJzL2yUzIrS7OUzNzN3M3KsorddF0sqoyyNPv93zf7wz80W63MyrLyp2SWkllKZu5m7mYpmbmbuZulZOa+ya2069eUzMrSzL2KlaVZSpU+16sa9pnufEpVA+aBRvjbKC9kHm708KGH8ECik0McGx3B\/61vVT8csQ5+BJt89NWq2jf3MRchzoiEHmsyP9fJQfCvVtUPZujUhU48hPV2+SIedTHGXEbijKyBYCPoCDGEefgZQ7AR4jESYx42OjHs0LFD8JMXNjqCvVpVPLHJQ9CJoTMi6Kwdgg\/BTz4jgm\/XI0ecfCR0RuYirI9NnPu1WlX14UMihxg5+BD8IX\/8x3fs00\/v25\/+6Z31pyT85CLkMrIPI4KPnFzwxx7kksNzETpj5HPelDbPEHOJkYMwFx\/CmuEjB8FmJC8Em\/yw0fGtVpt9WAs\/Qox10BHmUS+1FQXX60JOCHPRY010fKyLsCY2fmwkfOgRw4eNpGSWr4eOP9aIEV\/MR8fPOgg6vtXKjPNgh6xW288G+azDiJCHzcgajPb1n6L4+rLjPz\/4wXYwJbNYh3Xrwtr4GJG6zrnZuy7k4Ys56Ah+zrpdxTysyTbsR48eWciQbmVK29W4b9tYKZm5Vz8E2LmkZBYPa0pm7lXU3cz9+hz3Ks7VnWslKZnFOpWnurpXY0rV6L69ZvxgVtH2q7uZ+\/W4+8bXtD9RdzN3tErcqzG\/xtyUcm+lR6yyrl\/dzdyv+1OqfO7VyJXzpmSWEta2uJuxF7IdqSz87mbulb3rulhsR1MyY\/62d2OlVOlRH5Y7V7OUqtG9GuvXf\/zHO\/bzn393a333Kos9U6r0+tW97jFz3\/ioZWOZuZu5Vw2Ul2wey3X33DKjhm3PdcvdzN3M3Yx89yqnXkPl3VxTMktpYzdpy+XGy9obq9Lcq5GrO9fr4n7dl3vcbf3JNPcdortvZ3NuJPfGNwicIaUqkut4UjJLCW0jn3yy0Q\/VUtqe4b5t55b7xqKuumyiZilVlrtZSrb+C0LlGdL1\/LVMsmHTqF+\/fm0h2OdHOa0d3Kd1nqGdpv4p5ZiX5NDOVq\/Hve4Zll0U2\/UUxbY9Vsv9dpUvFrvnLZdmy+XunHNHi+LcOwxz\/ck1bJozjTrHjY0\/951bdz\/3Dlp\/rASKwmy12q5+udy2ZXVHYLHobq+p7PSTn0zlJOM6xzENe1wnban26urKziFF0bLhBNx37lzZYjGBg\/R0hD\/6o+ZnDn9PJc12W57lpp\/\/5fI4JKx73Ar9zKbuJh51H3n9VFjt+v3vN\/8M1eus29Xs8V5n27A\/+ugje\/z4sT19+tQuLy9PLn\/2Z5fvn4p4uBnr8vDhD+273\/35Wi4u\/mk9hxwURgQ9BBsJm\/GDD\/6fxVzsEPIQbEZyEOx8Dj7s3\/u9y\/U65OZCnDm575\/\/+eeGjb9JOFPMQyeHfHzYjPgQ\/Ag6McZd0jQ39zG3aR32IA\/m5OSCnzg+9MjB17QWecTgxsgcfAg+RoQYwhqMyN\/9XfMz96\/\/erlmSg7C\/Lrgj73QkaacJh\/zyEeIMyLUy0gcf0j48IfOSJwx\/Oj4OCOCzUgcIYbgQ9AR8uoC9\/v3P33\/HJIXkq+Fj7XCxzrY+NFzPzaCj\/XRI485TT\/\/n332yMjlZ4Ic5iB1VvjytcinfkbmEUPIQ0JnpB5G\/KEzxjz8EWfN2Dv8EWOsCzkIc1iPkbVZB5sYc+rjxcWv9n4XckbWZU3WQlgPQUdir4vau40YeQg6gs5ajOQj4WcMId72M9R0L3Mf73ve+7z\/Y70xjdNt2DfcBW7Yixcv7NWrV2eRv\/qr\/2u\/\/OU\/2d\/\/va3lL\/7C7B\/+4Wot+EL\/679+Zr\/4xf9eC\/6\/\/MurdT5z\/\/zPP1rnoxMLYS7rhf23f\/t4PQcf84lHDD0k9+VzfvGLjwybmsmJ\/BjxRQ34sKkZHSGWCz7i5FETOj5sBJsRYV7EqJ1Y+PAj2KyDjsQ8fNQePr6mQ8fPefJ5zCHGCHN0hD0Z8TMi6JGDTk34kViTkb3ZBz95MeKLdbmHxFiDEYFz23NHnHUQ9mAdzoNNDGFf7BBy0BmZg85ILjpjCDZCLj506kVH8LfNjX2ph1yE+QjzOCOCzUgcifXwIcxHiJHLiKDD\/W\/+5v+sn2ds1o2RHOaFsBY+hBzsyKfW2Jd7QJw81kePGGdvuxfkcq9Yl7nMI58RO9YImzzyqZ81sSMvcupziEcMHWEeefg5Dz7WzPfGDwdi5CHMQYhhI8xhPUZyWQcbnThj5DOXGLXvI5yVdVmTWmI9xliT9bDZhxzuC3aTcJ9Yi9wQ5oTOSI2M7L1PjU05vPdvaA+DDc+2YXNHaNpPnjyxc8piYfYHf\/DR+pc0Yh\/s0Ovjn\/zJR+t8\/Mulva9tsTBbLKq1IrZYmC0Wm5zl0izmLxab3Pp+ub1cbvJYd5fk8yLvJt9yWdVHHhLzGBeLKoZOjNrREewQ7OWyqhMfNrJcmi0WmzVWqypnuazGyGHOYrHJw48PYU9GfEhdD5sRIWe5rNZaLKoRH0J8sah8se5yuamFOEJum9TjrLNcVmssFtVIDhJrkIPNGD7GxaLKXyyqEV8IucwJO3T8+BaL6hzoCHEEfblsXo9YmywWm\/WWS7PlcrMG6yJNc6mHGCPx5dJsuTRbLjfrLRYbPfKYQz7jcrmJ4wshFvpNY1su\/pCmNRaLzd7kkcOIoC8Wm\/hisdGJIZwncrERbPzLZcUQGyGGEGPcJeSH7JO\/ay1iy+WmFmzWZH30kOXSbLHYziNGHrJcVudHx48sl1X+YmG2WFyPk3Oo8N7n\/T9G+fYYi55AzTqCCIiACIiACBxEYHINu+k3wvkNcfwHkVGyCIiACIiACAyIwOQaNmxpzjTpEGz8kj0JKE0EREAERGBwBCbZsKFMkw7BloiACIiACIjAmAlMtmGP+aao9p0EFBQBERCBWRKYRcOOr8brY37H81july4CYyTA89xWN7GQek74Gesx2SIwBgJNzy6+JsnPk8dz\/5D0WTRsgMfX4\/mIH+FG5X5s\/BIROJjAACbsen6JtT3ru2IDOJZKEIEbCfAMtyXlz33okcu88DFiR2xI42wadht0bgw3KI9j48990kVgDAR2PbfEeLbzc2DjR9CbYrlPuggMlQDP8G1qY95Ynv3ZNGxuSshtbqrmiMAYCPDiQVpqlVsEJkuA5x5pO2C8\/xnbcobun03D5kaGjPmGDf2BUn0iIAIiMEQC8f5nHGsPmEXD5gblDxD2WG9Yfg7pIjApAjqMCJyJAO\/8fGnsMfaAWTTs\/EZJFwEREAEREIExEphFwx7j36TG+DCpZhGYMAEdbcQEptIDZtGw688ZN4+vRPAzYqOHYOMPW6MITIEAzzTPdn4WbPwIelMs90kXgSkQ4FnnmecsjNjoIdj4wx7KOIuGDXhuQAh2fgOwI8aInceli8BUCPBs84yHYMfZ0MPPiB0xjQMnoPJ2EuBZ5pkOwc4nYEeMETuPD0WfRcMGNjcgBLsuEWOsx2SLwNgI7HqOiYXUzxV+xnpMtgiMgUDbs4s\/pOkcEWNsig\/BN5uGPQTYqkEEREAEZkhARz4RATXsE4HUMiIgAiIgAiJwTgJq2Oekq7VFQAREQASGTWBE1alhj+hmqVQROCcBftnmnOtrbREQgeMIqGEfx0+zRWASBGjWQ\/5lm0lA1iFE4HACWzPUsLdwyBABERCEuBtDAAABK0lEQVQBERCBYRJQwx7mfVFVItAZAX267gy1NhKBowgMrmEfdRpNFgEREAEREIGJElDDnuiN1bHmRYBPyZyYEUFH0BH0fYX8Jtl3vvJEQATOQ0AN+yCuShaB4RKgyfKLYwg6go6gN1WOn3gew86FGDajRAREoD8Catj9sdfOInBSAvWmWrdvs1lTQ7\/NOpojAiJwPAE17OMZDmYFFSIChxC4qRnfFD9kL+WKgAgcT0AN+3iGWkEEJkdAzXpyt1QHmgABNewJ3MRxHEFVDonAroa8KzakM6gWEZgbATXsud1xnVcE9iBA067LHtOUIgIicEYCathnhKulx0Ng7JXWf8Fsl00jrsfz8xNrkjxHugiIQPcE1LC7Z64dRaBXAjTjXgvQ5iIgArci8D8AAAD\/\/3+iH9IAAAAGSURBVAMAaFkzcfj0hScAAAAASUVORK5CYII=","height":237,"width":394}}
%---
%[output:8253ebe1]
%   data: {"dataType":"text","outputData":{"text":"[21:18:51][INFO]  --- A08 | Section 3: Peak Detection Using Signal Processing Toolbox  [Analytics L3] ---\n","truncated":false}}
%---
%[output:140bf7ea]
%   data: {"dataType":"text","outputData":{"text":"[21:18:51][INFO]  Unknown Substance 1: Detected 3 peaks\n","truncated":false}}
%---
%[output:86391e69]
%   data: {"dataType":"text","outputData":{"text":"[21:18:51][INFO]    Detected m\/z values: 92.101, 118.098, 136.102\n","truncated":false}}
%---
%[output:9db13e61]
%   data: {"dataType":"text","outputData":{"text":"[21:18:51][INFO]  Unknown Substance 1: Extracted [M+H]+ = 136.1021 Da\n","truncated":false}}
%---
%[output:11cc638e]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAewAAAEoCAYAAACaU8LCAAAQAElEQVR4Aeydf6xlZXX3F3bs\/HhtAhnG6bRD2U4Vo7YFktoZTLk90\/5hK6+FpqTWSRnuAGK0tRejvlQjcK0TaRnrC9bYlncIqORqbG1itBJrhGs0QZpqgDY2GVPcZK5hNEVNMIMpk\/ruz96umefu2eecfc49Z\/\/83szaz3rWWs\/zrOezn7PXPedeLs\/7sb5EQAREQAREQAQaT+B5pi8REAEREAEREIHGE2h2wW48PiUoAiIgAiIgAtUQUMGuhrNWEQEREAEREIENEVDBnh6fRoqACIiACIhAZQRUsCtDrYVEQAREQAREYHoCKtjTs2v2SGUnAiIgAiLQKQIq2J26ndqMCIiACIhAVwmoYHf1zjZ7X8pOBERABERgQgIq2BMCU7gIiIAIiIAI1EGg8QV7z549Q7mM8uUHTRKbH9uEPvmHMmlOjJ10TJvjN7TfGWyc9YtkBlOXnqLu9csmms+z7Liq48hzlmtOOh\/xeZllPk2Ziz2Oy4WYIhk3blI\/a0w6Zp7xjS\/Y89x8W+bm0DzxxBMWCra25F91nk1hE94v16vOzdf1tur1x9178vHcvMU2blxf\/c7I2z6zcgZh23UeKtgNf+VzADmQ+TSx4cvb59Wvcq2N7KHpec7gvm0ET\/pN3zSMphkzLlHmhEc+Dhu+vL3Kft3rl93rtKzmtb95zTtvHmXnrztOBbvuO7CB9XmxbmB4J4fCBOnk5ma0KfjU\/WAdtxVyHBcjf0YAVk2\/n1mmum6UQGcKth9YWpdRcIhxv+u0Lu7z1u20bqMt6hfZwlj8LthHib8YiR8VV+QfZsOOhPPRD8V9bqMf6mE\/b8eHYPcWHaEfCrZQQh\/6KB\/+tgn74Z563vRDcTstdloX+oj3afN9bJMKc4QSjnc7tlAP++vsOEoKHMqMJYYpaRF0BB1Bzwt2l7yPvvto6bvQR+jTIuihYHMJ7ehup6UfCjaX0D4P3dehDeenj2CjRdBd6Lu4LWzd5637wn6o5\/343OYtNhe3zbL1uWnz82ILJe\/P94l1G3oobp9n25mCDSTg8RBwoY89L9iJCe1uw47Qdz86Nhf67pu0ZazPQ0t\/3BzEIcS6jBtT5Gcs87jQJ47Wbd5iw+f9vI7ffd5iIy4UbO6npe9+dGyhYBvmD30es9GWOYcJcw\/zYcc\/SojJC3v1Mfjoh4LN\/VW0rBeuj47N16aP0KdF0IlBDwUbvknExzPWpWg8vnxs2A\/HhLHE0C\/rJx4hnhZBd2EubC70N+rz8bNqycnzo6Xvc9NH6NMi6Ahx9F3oY3eh7z5vseH3fl6nT4z7aeljR9CxudDHPq0wnrl8vPexIfSH+fJ+j\/OWscTQd52+CzZ885ROFWzAjYMF1KK4IhtzFcUTix1\/qNPPC3HEuD3U3Va2ZawL85Yd53GMdX1cOy52nN\/nHxZH\/kU+txX58WH3uWfRMucwYf5hPuz4RwkxoeRj8eVtYR\/\/qP3iIyYcM6k+7fhpxw3Lj\/lc2Fc+Dl9oy\/fdx9i8jz52MzNa+h5PSx87+jghtiiG8XkffexF8fiK7Bu1sV5+bvrYR82Nn7gwhj52bLT00UMpsoX+onGMwR7GuY7P9XEtc+QlHI8v7DMffeyu05YRxjB2VOw4\/6ixZX2dKtjjNu3QacfFTuNnXm4agj7NHPkxw+aZ5Ro+F2u55PMo6nust0UxG7H5vGG7kfnqHuucwzzCvaGHvlDHx3gEPfRtVGe+UMrOF45BLzvO44aNmccefc15tuwnL\/Ncb9jc+RzoD4sN7cTlJfRPq+fnpD\/tXOE4zkleQj86a+UFu8soXxjDOsS6zfvYXNw3z7ZXBRvI84TZ5rlhEwqHcNR+8Ifx6KPip\/ExZ5FMM1cTx1TBsGjfvi4+10PO2MfJtOPGzdtmf8gw1EvvacpAvxc+PFw71N0\/rA1jQ31YfFl7OFeolx2\/kbhwvVBnTueWt+MLBX\/Ydx17KMznvnm1jS\/YACkCgQ3fNGAYx\/gyY4tiGYvdx6OXsXn8JK3PnR+TX6\/IX8ZGDHPRbkQmnWPUvshjnJ+YrkkRQ+dA6\/tFJ5bWbWXbcePwl50rjJtmHPkXjcOGL5y\/rM44xofx9LFjo6WP7kIfu\/enaRnPPPmxRTZihtnxTSLMw9o+Bh2b970tsrmPdty4cX7mKJJJx43Ls2iNYbZZrx3ON8s8h+VfZG98wSZpBwUkF2z4hglx+LxFD4Xxw3xhHLrHEo\/Qx16VsB7rhoItXJ\/+KD+xw2LydubBxhgX+tgRbGEfG33skwhjGBsKNp8DPfShY3N\/G1vyZx\/k7jp9hD72WQrzhpJfg\/4ov+cSxmEL+4ynjz0U7GG\/SGcccaFgK4ota2P8qPnG+VknjKFfRsIxvj42xtK6jZY+9lCwh\/0inZhQiubBFsagYwvno48dcXtow45g+4n\/9H\/Dj90l9BNH3330kdDmPmxFPrfjm4Uwn6\/pLTbmpnUbLX3so4QYj6UNBd+osbPwtaJgs9EQRqjjQ0IbEOm7eJ+4UPDT9xbdJW+j7+IxYYsv7KPnbfl+UQy2ImFsKONi8BNPi7hO64LdxW3euj1s8z7v0xLnLTqS7xfZiAmFmFBCH3roG6VPEjtqnml9o9YPfegurIVOG0pZWzgGnXF5wZ6XMAYffdq8YEfcju6CDZ12UmFcKPnx+ELbuD6xxLjQz4v7aPM+7+NDwr7r3oZ+bPRDweaSt9N3X5mW+LwMG1cmzmPCOdzmbehDd7u32PJS5HObt+EYt9Fi9xZ9mJSJ8bHEhuJ22rydPnaXfB+722hDwTdvaU3BpuiWgUEcEMvEKkYERGA+BPQ6nIyreE3G66zonhhaU7ApwkhP7ou2KQKtJqDX6mS3T7wm49XX6NYU7LI3KH\/wR33nura2ZhIx0BnQGdAZ6M8ZKFtLZhg3s6k6V7BDMuOK9Tve8Q5bWFiQiIHOgM6AzkBPzsDrX\/\/69I1aWCvaone2YI8q1twcvqN+5JFH7MiRI7aysiIRg96cgaWlJV4COvs687058\/6M5+zz3E9fAC28zKVg181hXLEO89u7d6\/t27dPIga9OQOceV4DtDr7eu336Qxw5jn7bZXOFexJinVbb5ryFoGNEPjZn\/1ZO3jwoNFuZB6NFQERqJZAZwo2hdrRoefFfWbSRKDfBCjU1157bb8haPci0EICrSvY+d8Cd+Zupy0Sj1MrAiIgAiIgAm0k0LqC3UbIk+SsWBEQAREQAREoIqCCXURFNhEQAREQARFoGAEV7IbdkGano+xEQAREQATqIqCCXRd5rSsCIiACIiACExBQwZ4AlkKbTUDZNYsAf5xI0p8\/+VnHvW7WiZ9\/NirY82esFUSgdwR4eOtP\/+rPHs\/7Tz+3+c+MTvNQUMGehprGiMDEBPo1gILNn4DUn\/7Vnz32Pws667btf2Z0mieCCvY01DRGBESgFAH+FGSf\/vSl9lrdn3rlbJU6hB0KUsHu0M3UVkRgWgIaJwIi0HwCKtjNv0fKUAREQAREQARMBVuHQAREoOEEmp9e\/v9d4P1ZZ+7zhu2s19B8zSWggt3ce6PMREAEWkSg6P9hQGGd9Rby65Rdo2zcrPPVfLMjoII9O5aaSQREoIcE6t4yBVzFuO67UM36KtjVcNYqIiACCYGTJ09aFyTZykT\/KKih+GBs43T3l2mZLxQfgw2dFgl1+gg2SbMJqGA3+\/4oOxHoFAH+++wHH3zQ2i7Hjh07675Q9PLCu18CaUMhDvusJVwDfc+ePekS6Ci0SKjTR+aVE2tJZkNABXs2HDWLCIhACQK7d++2yy67rPXCPiz3RdFzwYVO60JBdHGbt9iJp3XbtC1zuIybw+Nox8XKXz8BFez674EyEIHeENi2bZtt37699cI+Rt20fPGlIGJzGTV2I75J1pkkdiM5jRor32QEVLAn46VoERABEZgZAQq4F04mzfexjZNw\/LhY+dtNQAW73fdP2YuACDSUgBdf0nOd4opgm1YYHwpz+1zooc\/ttKEv32cMNklIoHm6Cnbz7okyEgERaBkBimFRyqEdPRSPx+Y6bb6PzQVfXtzn7Si\/+8rEeoza5hBQwW7OvVAmIiACIiACIjCUQFiwhwbJIQIiIAIiIAIiUC8BFex6+Wt1ERABERABEShFoD0Fu9R2FCQCIiACIiAC3SSggt3N+6pdiYAIVEyA37TOS8UprFuOXNYZpuzMap4pl9ewgIAKdgBjA6qGioAI9JgARc1\/AztssVeFpcq1qtqT1llPQAV7PQ\/1REAEREAERKCRBFSwG3lbZpyUphOBhhA4dMjsnHPaL\/fdZ+u+eFfNO1wkdGCn73ZaBBuCjqC70A\/F7bShHR0b4jotgs2FvovbaN3mLTYXt9G6TW39BDpbsHXQ6j9cykAE+kSA4ozw7EHye8eGH0FH0BF04mnph4JtnI94YmgRdISx9F3oF9nxD\/MRL2kGgU4WbD94zUCsLMYQkLtHBO691+zHP26\/LC4Ov2kUPyT\/HMIWjsr3Q9+s9FFrkF8os1pT88yPQOcKNgdwfrg0swiIgAicTaDouUOxLLKfPboeC\/nlpZ5MtGpZAp0r2H4AywJQnAiMJCCnCPSEQJO\/uejJLRi7zc4V7LE7zgWsra2Zy6lTp0wiBjoDGz8DuZdZ57u8UaDg5QX7JJsnftgco3ysEfrpj5Iw1tfDxhhat9Fia7KMe7368522yfsok1vvC\/aBAwdsYWEhlcOHD9vx48clYtCWMzBVnjy4Tpw4MdXYsq8P5i\/zAOpSDIUuL74\/7K7TjurjC4V4l9CO7nZvsSH0vUV3CW3ooXgMbWhHx9ZUeeqpp0ae5aNHj6bPd57zPO+buo8yefW+YB85csRWVlZSWVxctF27dknEoNNnYOfOnbZjx4657vHcc88t8\/xRjAhsmMC4s8xz3Z\/xS0tLG16vzgl6X7D37t1r+\/btSyWKItuyZYtEDDp9BrZu3WqbN2+e+x6tzieb1u4NgXHPbJ7r\/ozned9mML0v2G2+ecpdBHpJYHXVbHXVbHXVbHXVbHXVbHXVbHXVbHXVbHXVbHXVbHXVbHXVbHXVbHXVbHXVbHXVLI57iU2bbj8BFez230PtQAT6R2D\/frNpJY7hJRGB1hHobMFu+i9KtO6kKGERaAqBwcBsMJgum8HAbDCYbuyIUfw2dV5GhDfWxR6KksNeJEWxZW3MVzZWcRmBzhbsbHu6ioAIdJIAfzJt0o1Fkdltt006qnQ8bxJCmaQgnRVbetWzA2c5Vzh7uDfX57VWuK70MwRUsM+wkCYCItAWAlFktrg4WbaDgdlgMNmYDURT1FTQNgBQQ88ioIJ9FhIZREAEWkFgknfZUWQ2SfwcAFC8Q\/ElsKHTIugIeijYXEI7emhHx4agI+ihYHMZZnd\/QVtoCudB9yD0UNyeb4lxG7qL29SaqWDrFIiACLSXQNkiPMePwsvAo\/jwjjsUbIzF5q3r+NBDwUYcbWhHx4YP3VvX8aGHgo04H5aFCwAAEABJREFU2tCObZgQmxfGEo8dPRRs+EIbutvxuWDDR991+gh97BIVbJ0BERCBNhNYXDQbDEbvIIrMFhdHx1TgpfCEMm7JMBZ9VDyFbZSf8aGMih3mY428hLHh\/OjDfKEdnVjmRUfQsaEj9GlLS4cD9Q67wzdXWxOBXhAY9+657LvwOcOi8ORl1JL5WPqj4kf5GJuXUfHT+PLz02ceii+6CzYX99G6jZZYbAh9SUZABTvjoKsIiEBbCQwGZoNBcfaDgdlgUOybs5ViQ+EZtgz+Yb4i+6j4Ub5J5yqKn8ZWJqciPj4OH+L9aXJo4JgNpaSCvSF8GiwCItAIAkXvoqOo0l80o7CEQrFxNuihDx1bkR8bPmJCwTbOl\/fn+z7fsLmIn0aYz+f2Fhtz0bqNFltePAa768Qi9LFL9DNsnQEREIEuEIiis\/8b68HALIoq2R1FJS\/5hcv6fdyo+FE+xrsfHfG+t9hc3Oat28MWX9gv0okJJYwJ7ej4vEVHwj66Cz5JRmDu77CzZXQVAREQgTkTWF42i6JskSiq9N11tqiuIjBfAirY8+Wr2UWg1wT4f29XKu99r61t2mRrN9xgla67tqb1KmbQxxdWzwt2H2+59iwC8yewe\/du439leODAAVtYWKhObr3VFi64wBbuuae6Navcn9Y6fV85W5wxztr8T3QzVnheM9JQFiIgAl0iwEP0yJEjtrKyImkhgze+ccWeeiqTl760ufeQM9al1824vahgjyNUo19Li0CbCVC09+3bZ5L2Mdi5c5\/96EdnpKn3kDPW5tfIpLmrYE9KTPEiIAIiIAIiUAMBFewaoHdjSe1CBERABESgSgK1FGz+Y\/gqN6m1REAEREAEyhOI4\/KxiqyOQC0Fm\/8gnqLtUt12tVJfCGifIiACItA1ArUUbCBStF1UuCEiEQEREAEREIHhBGor2GFKKtwhDendJ6AdioAIiMDkBBpRsP0ddli4J9+KRoiACIiACIhAdwnUVrC9SNN6oe4uZu1MBNpDQJmKgAg0k0AtBTss0hTrPJoiWz5GfREQAREQgfkTiOP5r6EVyhGopWCrIJe7OYoSARHIE1BfBPpLoJaCzTvsIuTD7EWxsomACIiACIhAnwjUUrD7BFh7FQER6A+Bruz0ySe7spNu7aPSgs07aASEtHnRR+WQkYiACIiACIjA2QQqLdgUZIQ0aPOCXSICIiACIjAPApqz7QQqLdhth6X8RUAEREAERKAuArUUbN5ZT7vh8GP0YXOEMejD4mQXAREQARFoBgFlMZ5AJQWbool4OujDxGOKWsZQ7F3o5+Owud9bbPk49UVABERABESgTQQqKdheOB2M94taj8m3FF3iQzt97KFNugiIgAiIgAjMlkAzZqukYDdjq8VZrK2tmcupU6dMIgY6AzoDfT8D4dMyjq3Vz0V\/vtOG+2qjXlvBDt8ZoyOzAOjvupnPBduwuQ8cOGALCwupHD582I4fPy4Rg06fAR5cJ06c6PQe9Tre2HPsmWd+uO6R2WaeR48eTZ\/vPOd53q\/bWMs6wwr2XLdBIfUi6jp99I0uzBzMFQq2YfMeOXLEVlZWUllcXLRdu3ZJxKDTZ2Dnzp22Y8eOTu9Rr+ONPce2bNmy7pHZZp481\/0Zv7S0tG5fbevUUrCbBGnv3r22b9++VKIoMg6qZIs4JA+srp6DrVu32ubNm3WPO3yPN3p2n\/\/8Tese0xudr87xPNf9Gc\/zft3GWtZpZ8FuGWSlKwIiIAIiIAIbJdCags1H3PmPtulj3ygEjRcBERABEThDII7P6NKaQ6CWgk2Rpdgi6OAIdfpFQixxLvQ9Dhs6NvRQsOGrSLSMCIiACIiACMycQC0Fm11QRBF0JNTpDxPiXMIYbN5HD8XtakVABERABESgrQRqK9jhO+BQbyvI1uStREVABERABFpJoJaCTYEO3wGHeispKmkREAEREAERmDOBWgr2nPek6dtLQJmLgAiIgAgMIaCCPQSMzCIgAiIgAiLQJAK1FGw+Audj8SaBUC4iMJaAAkRABESgRgK1FGwv1rR5qZGFlhYBERABESggEMcFRpkqJ1BLweYd9jCpnIAWFIFuENAuREAEOk6gloLdcabangiIgAiIgAjMnECtBds\/DmdX6LQSERCBDhLQlkRABDZMoLaCTYHmY3HfATo276sVAREQAREQARE4Q6CWgk1hpkCfSUOaCIiACNRGQAuLQCsI1FKwW0FGSYqACIhATwnEcU833vBtq2A3\/AYpPREQgZ4T0PZF4CcEainYfBzOx+I\/ySFt6GNPO7qIgAiIgAiIgAisI1BLwSYDijNFGp2WPrpEBERABESgNQSUaIUEaivY7JEi7UJfIgIiIAIiIAIiUEygloLNO+qidIbZi2JlEwEREAEREIGRBDrmrKVgd4yhtiMCIiACnSQQRdm24jhrda2XQKUFm3fQCFumzQsfj+OTiIAIiIAIiEDHCUy8vUoLNgUZIUvavGCXiIAIiIAIiIAInE2g0oJ99vKyiIAIiIAIiIAIlCFQacH2hHhn7bpaERABERABERCB8QRqKdiklf\/5tffxSURABERABERABNYTqKVgU5x5l10k69Orsqe1REAEREAEIBDHXM2iKGt1bQaBWgp2M7auLERABERABESgPQRUsFtyr5SmCIiACIhAvwnUUrD5KJyPxfuNXrsXAREQAREQgfIEainYXqxp81I+dUU2h4AyEQEREAERmDeBWgo277CHybgNhwV+VGzZuFFzyCcCIiACfSYQRdnu4zhrda2XQC0Fe9otU4TDQk+\/aC7sZeKKxsrWPQLakQiIgAh0gUAlBZsCWgbWqDh8FOFwHvrYQxt97KEt3w990kVABERABESgDQQqKdiAoJCOE+JmJeFas5pT84jA7AloRhEQAREoR6Cygs273HFSLuXxURTrcC36w0atra2Zy6lTp0wiBjoDOgN9PwP552WbefjznTa\/r7b1KynYFM8yYMrGjZtrknkOHDhgCwsLqRw+fNiOHz8uEYNOnwEeXCdOnCi1R70e+vk88GfsM8\/8MFV\/8IPvt\/a8HD16NH2+85zneZ9uqKWXSgp2k9kcOXLEVlZWUllcXLRdu3ZJxKDTZ2Dnzp22Y8eOTu9Rr+ONPcf8mb1ly5ZUPffc81p7Xniu+zN+aWkp3U9bL70v2Hv37rV9+\/alEkWRcUAlW8QheVB19Rxs3brVNm\/e3IF7rHM6rzPqBe35z9+Uqps2bWrteeG57s94nvfW4q\/WFGw+5s7\/LJo+9pA\/feyhTboIiIAIiEA5AnGcxSXvXzJF18YQaE3BhpgXYwoyQh87Qp8WwU7fhT52iQiIgAhURUDriMCsCbSqYLN5iq8LfRdsrtPSd6EvEQEREAEREIE2E6i1YPs7YACi00pEQAREQATmTUDzt5FAbQWbAs07YIeGjs37akVABERABERABM4QqKVgU5gp0GfSkCYCIiACIiACGQFdiwnUUrCLU5FVBERABESgSQQuvDDLJo6zVtd6Cahg18tfq4uACIiACLSKQH3J1lKw+Ticj8XDbdPHHtqki4AIiIAIiIAIZARqKdgsTXGmSKPT0keXiIAIiIAIiIAInE2gTME+e9SMLBRplxlNqWlEQAREQAREoJMEainYvKPuJE1tSgREQAREQATmRKCWgs27aoq2y4b2psEiIAIiIAIzIxDH2VT6W+IZhyZdaynYAKBou6hwQ0QiAiIgAiIgAsMJ1Faww5TCwh3aO6BrCyIgAiIgAiIwEwKNKNj+DpvCPZNdaRIREAEREAER6BiB2gq2F2laCjXSMbbN344yFAEREAERaA2BWgp2WKRVqFtzVpSoCIhAzwj4L549+WTPNt7Q7dZSsFWkG3oampWWshEBERABEQgIVFKweUeN+Lrow8Rj1IqACIiACIiACJwhUEnB5h014suiDxOPUSsCjSag5ERABESgYgKVFOyK96TlREAEREAERKBzBGop2HwcXkRymL0oVjYREIGhBOQQARHoIIFaCnYHOWpLIiACIiACIjBXApUWbN5BI+yINi\/8XBufRAREoMMEtDUREIGpCFRasCnICJnS5gW7RAREQAREoD4CcZyt7f8NdtbTtQkEKi3YTdiwchABERCBEQTkEoHGEqilYPPOurFElJgIiIAIiIAINJBALQW7gRyUkgiIgAg0n4Ay7DWB2gp2\/hfOvN\/ru6HNi4AIiECDCPjPseO4QUn1OJVaCjbFmY\/FEdjnW2wSERABERCBVhFQsnMmUEvBLtoTRZtCXuQLbcS4hPYivWxc0VjZREAEREAERKBJBBpTsMtAoQBT2F3olxmnGBEQAREQgZ4T6MD2W1OwKc4U6pA5feyhzXXs+L2vVgREQAREQATaTKCWgk0hpaACznX66NgkIiACIiACItAjAqW2WkvBJrOwOKMj2GchkxT\/tbU1czl16pRJxEBnQGegz2fAn8HPPXfK1VY\/F\/35Tnt6Qy1VKinYFNCyslGOrDNJ8T9w4IAtLCykcvjwYTt+\/LhEDDp9BnhwnThxotN71Ot4+ufYo49+P30Mn3\/+M\/bUU0+lOsW7rUyPHj2aPt95zvO8TzfU0kslBZsCmpdh\/VlwpGi7MB86bZEcOXLEVlZWUllcXLRdu3ZJxKDTZ2Dnzp22Y8eOTu+x8HX8d39nuzYin\/98L5i94AU\/kz4qac8777xU37p1S2v3znPdn\/FLS0vpftp6qaRgVwkn\/40Aa2OjLZK9e\/favn37UomiyLZs2SIRg06fga1bt9rmzZs7vcehr+O\/\/EvbMq088kgvmG3atCl9VNLCMe0kF\/Q2Cs91f8bzvE+20tp\/tRZs3vki0PMWvUgouvkY+tiL4jdm02gREIHOEVheNku+KZ9qX4y7996phmqQCMyKQG0FO19sKbzYRm3MY4hD6Hs8fdfVioAIiEAhgdtuKzSPNU47buzEzQ7g+xQyjGOukroJ1FKwKa5hsZ0EAuNcwnHYwr7rw+zub2OrnEVABKYksLhoNhhMNjiKzBYXJxujaBGYA4FaCvYc9qEpRUAERKAcgUnfLeuj8HJcFTV3AirYc0fctwW0XxFoOIHBwGxxsVySi4tmg0G5WEWJwJwJ1FKw+Ziaj8XDvdHHHtqki4AIiMBcCJR51xxFZpO+G59LsppUBDICtRRslqY4U6TRaemjS0RgngQ0twicJjCuGF97rVkUnQ6XIgJ1E6itYLNxirQLfYkIiIAIVEZgedksioqXiyKz5eVin6wiUBOBSgs276SRYXsd5Rs2RnYR6A4B7aRyAsPeZZf5yLzyZLVg3wlUVrApxv5uGt3Bo7vgd7taERABEZg7gcVFs8Fg\/TKDgdlgsN6mngg0gEAlBZuCHBZjdGwIuksDeCgFERCBIQQ6aw7eZccW9f4XzZ580tKvKEobXRpEoJKCPWy\/FOphPtlFQAREoBICg4HdZ8k7bTNbtYHZYGD6yghEUdbGcdbqWtjBSYIAABAASURBVC+BWgt2vVvX6iIgAt0hsLGdHLJ7jXfX77HbLI5NXyLQSAIq2I28LUpKBESgagIvsm8ZRdv0JQINJaCC3dAbo7REQASqIRDH819HK4jALAhUVrD5BbNQSD7su45dIgIiIAJ1EYjjulbWuiIwmkAlBZtfLisro9OVVwREQAREYLYENFtbCFRSsNsCQ3mKgAiIgAiIQFMJqGA39c4oLxEQgVoIxHEty2rRIQRkPkNABfsMC2kiIAI9JBDHPdy0ttxKAirYrbxtSloEREAE5kMgjrN5oyhrdR1FoFqfCna1vLWaCIiACLSKQBRl6cZx1upaHwEV7PrYa2UREAEREAERKE1g0oJdemIFioAIiEAbCMTx+izjeH1fPRFoCgEV7KbcCeUhAiJQC4E4rmVZLSoCExPoVsGeePsaIAIiIAIiIALtIKCC3Y77pCxFQAREQAR6TkAFu7oDoJVEQAREoHUEoihLOY6zVtf6CKhg18deK4uACDSIQBRlyTz5ZNbqKgJNI6CC3bQ7Ulc+WlcEREAEEgJxnFySf1GUXPSvUQRUsBt1O5SMCIhA1QT8HXUUVb2y1hOByQioYE\/GS9H1ENCqIiACNRGIomzhOM5aXesj0LqCvWfPHnMZhc1jaEfFyScCIiACEIgirmZxnLW6ikDTCLSqYFN8n3jiCXOhXwQUu8fQ0i+Kk00EZkJAk7SaQBy3On0l3yMCrSnYFF2Kb3hv6GMPbfSxhzbpIiACIjCOwIUXjouQXwTqJdCagj0vTGtra+Zy6tQpk4hBx86AzvSY13X+2dL3++88nEO+7\/a2tP58p\/W9tLXtXMHOv7se9477wIEDtrCwkMrhw4ft+PHjEjHo9BngwXXixIlO73GS1\/Gzz\/4ofX7\/4AffT9vnnjvVazZxnGKwn\/qp7Fn4zDM\/TA3wmYRrU2KPHj2aPt95zvO8TzfT0kvnCnZ4H8YVa2KPHDliKysrqSwuLtquXbskYtDpM7Bz507bsWNHM\/bYgLP2\/Odv4lFgv\/RLP5O23\/nOll6zSSEkF38WbtmyJemZnXvuea3kwnPdn\/FLS0vpXtp66WzBLlOsuWl79+61ffv2pRJFkXE4JVvEIXlIzeQc\/MVf2JaNyCc+MfN7sXXrVtu8efPM550Jr1lxn2CetbVNPArsxS\/OWjpt3css8mb\/iM\/lXODktja1PNf9Gc\/znr21VTpZsMsW67beNOXdMgLveY\/ZtPKlL7Vss51KV5sRgUYRaE3B5mfTFOKQHn3s42yhX7oIVEpgedks+eRmqjUZd++9Uw3VoMkJgJtRccxVIgLNI9Cagg06ijNF2oU+dgQbLYKeF+wSEaiFwG23TbfstOOmW02j2kagonz9Gxn\/E64VLatlCgi0qmCTP0Xahb4LNnTaIsEnEYFaCCwumg0Gky0dRWaLi5ONUbQIbJBAHGcTRFHW6tosAq0r2M3Cp2xEoCSBSd8t66PwkmA3HhbH2RxRlLW6zoSAJpkDARXsOUDVlCJwFoHBwGxx8SxzoWFx0WwwKHTJOF8CUZTNH8dZq6sINImACnaT7oZy6TaBMu+ao8hs0nfj3aam3YnA7Am0dEYV7JbeOKXdUgLjivG115pFUUs317604zjLOYqyVtezCURRZovjrNW1PgIq2PWx18p9JLC8bBZFxTuPIrPl5WKfrCIgAn0hMHSfKthD0cghAnMiMOxddpmPzOeUkqYVgZBAFJ3pRVGmf+MbJ9O\/sZ71dK2DgAp2HdS1Zr8JLC6aDQbrGQwGZoPBept6cycQx9kSUbS+jeOs37drHJ+94yjKbN\/97jZ7+umns46utRBoRMGuZedaVARqJLD6G2f+mEpsyRNx2Ltu05cI1E8girIcbr\/9okzRtRYCKti1YNeifSewagO7z5J32maW6rHeXZu+Gk\/gk5\/cZnHc+DQ7m6AK9thbqwARmD0B\/szjIbvXeHf9Hjvzbtv0JQINJPDhD5+0F77wZJrZ6mra6FIDARXsGqBrSRFwAi+ybxlF2\/RVC4E4zpaNovVtHGd9XTMCv\/M722wwyHRd6yOggl0f+5msrEnaSSCOs7wHg6z90peyVlcRaCqBl71sW5paHKeNLjUQUMGuAbqWFAEnEEWuqRWB+gnEcZZDFGVt0ZUf5xTZZZs\/ARXs+TPu8Qra+jACcZx5fuM3sjaOs1ZXEWgqgShqamb9yUsFuz\/3WjttEIE4zpKJoqzVtR4CcZyte+GF69s4zvq6niEQRZkex1mra\/UEVLCrZ64VG0JAaYhAnkAUZZa+fuwbx9n+\/RuYrJddoyhr4zhrda2egAp29cy1ogicJhBFmRrHWatrtQS8MEdRteu2cbUoyrKO46zVtXoCKtjVM9eKInCaQBRlahxn7ZmrNBEQARFYT0AFez0P9URg7gTiOFsiirJW1\/oIxHG2dhStb+M46+u6nkAUZf04zlpdqyWggl0tb60mAmcRiKLMFMdZ24arcuwmAf8RQTd31\/5dqWC3\/x5qBy0j8Oyz2Z94bFnavUg3irJtxnHW9vUaRcU7j6LMHsdZq2u1BFSwq+Wt1UTAvvOd7C9GRVEGI4qyNo6zVteNEig\/Po6z2Cha38Zx1td1PYEoyvpxnLW6VktABbta3lpNBGx1NYPgfzQl6+laB4E4zlaNoqzt+zWOMwJRlLX5q\/\/nXvpTunky1fRVsKvhrFVE4DQB\/zlhFGWmKMraOM5aXashEMfZOlGUtX6NokyL46yd9bXJ88Vxll0UZW3+GkV5i\/pVElDBrpK21hKBhEAcJ5fkXxQlF\/0TgRYRiKIs2W98Q7+HkZGo9qqCXS1vrdZzAnFstrq6HoI+ZlzPo6peHGcrRVHW9v0ax2ZxbBZFZlFUTGMwyOz\/8i\/b7ORJFe2MRnVXFezqWGulCQjM6mEQx+sXjePsobTeWl1vdfXMWlGU6VGUtV\/72tOZomutBKIoW\/4\/\/qNfBSmOs31HUdYOu0ZR5vnEJx62hx9+2J5+Wuc2IzL\/qwr2nBhziL\/4xS8a8uijj6YHG\/3Tn37U3vKWp+3uu4+dbjn0CH7k2LFj6YuAOYjHxxxf+crx03b6x48fN2w33HDc7rrr0XQtYpkD+ed\/PmbEYaNl3o9\/\/GF717uOGWMpitjcz1yMw4YfOzqCjpAT4\/7mb04auRHPfPhYg\/bw4ePGGPzsExvCvl\/3upOnc37ggZNGn5yIRWcO5CtfWUv5ZPpx+5VfedrYp6\/F\/KzNvMTQ\/tVfPZ1yJj\/WZc79+8327j2Z7pmY3\/3dp+1FLzK7\/\/7jact8cCLe56FFmB8OrIUQxxrYmYsYBB0\/+\/D7gB6OxU8+cOPI\/dqvnTz9LmZx0eyFLzxp\/\/Zv2+2ccyzdN3H33WdpnzkZz\/q33\/6wMTfFHRt7RScncnFhbe4nY9kja9OHx+c+d9Kee+45i2NLGcDi5ptPGveNGOYlnrEIOvNip0XIARbkRAzCmtiJZy586OSIeJ8YzgI+8sZOn3lZg3mw04ctOveHlj5CHOMZy9z00TmfjEdnX8QwD4ycg58z5rTkK\/\/Lf1GUGJN\/y8tmhw6ZwYycfV3W+uxnP5ueKeamj4\/8EHRyoGUcbDh3rEdePoZxCHHYyBXd54AJHLknzIMwnhgXYhnLvhFi8GH3ubGjMz8cyI04eHzqU0+nzwLmQE+2bZxN2mEyGGSed77zMrvllt3peWWPzME65Mha5EDu5IPg95Y492Gnzzha7ic6sczDawE782HDTz\/Lol\/XzhbsPXv2mEsdtzSOzfjlos9+drvdccc2W10142Okq666xD70oe32xjdelLbor3rVZXbVVRfbgQOX2e23X2QcUB4qHFri8THH5ZdfYC9\/+db0AfLqV19kN9\/8rGG7554L7KabLknHM+9b35rp11yz2y699BJjfuLf\/OataQxr8CC6667H0uLJi87nJwfWX1tbM95hLC+bXX319jS\/v\/7r7emL8w\/\/8Fl785u3JbZLjLWYj3U\/8hFLX8C33HJBktfuRL8g3Sf535K8sNnrJz+5LR1Lzq95zTajf9NNF6ex6OT79rdvN3Innrx4kFDMGONrLS9byuqTn9xqX\/\/608n43cY41rnhhrV03RtvvMjiOOPu45iH8\/C2t21PfdhhQ\/4HD2b5Mt8DDzxrrPeFLxxL9wID4lgDRuyJ\/WJj7PKyGfvw+4B+xRXbbHnZ0vGw515yBijOd975GGmclsHgtJqeC\/jenBRRrMzJPWHtd73rsvQevuY1W9N5eYBROLkP5E1O3LvHHnvMlpfNGMseWZuzsrpq9qY3bUm\/obr11idSBtdcc0F6Rm9J71t274h\/3\/suSsejs08YcD\/gBAPywc4aCGeKfRPPXPjQ4QBP4ukzlnvL\/YEHdvrce9Ygd9izDmcOHcbskT7nGR6+BvebcazH64b10NkX65EzuTkHzhlz\/\/u\/bwevRVHanL54ASe3++4zW121lDlrsy738vrrfzM9f8xDn\/PA\/OwN\/ebktUk8fdh85jPb0zNJXtyj5WVL79\/ll\/M62Z3yh0e4T5iwx9\/\/\/e3p+szFeObn3FOIn3322fT8U9gR+MOHXO64Y2vyGr04fWZkOV5iv\/qr25NvVteS+3pxmj\/c7r9\/zcjJf\/N769bsPzs8DSSnOJ\/vfnebffGLFyTfQF5knBXu4\/KyWZbjJcaz7\/LkmUW+sMjsF6drUazZF\/uFIVyvSZ5XxPDNDS3jYMJrYXnZ0mcN+8b\/QPLNfi6tXnSf18VdUqifeOIJc6Ff9T45WNdf\/1t2552XpAeaBy2Sz8MLCIcf4QXAGFpeoMRj\/\/jHL0I1dMbQui11JBds+BB0JDGn\/9CZM+0kF3TyIUd0n4s41r8lKbD4sPt86AgPn2SK9B8+FFrGMRd95sHmutvpux0dCWPpE4st03enDNFdGE8exLEmeaLjp8WH7nOgI4yjRfI+bPgR5oMN86ATi50YhD7C2rT4iEUP\/aGdvPC5XHRRdj+9\/7KXnf2Q9PmIYS3mQ0fwIRQx1sZHrsRdnxQTJL8m4xDGsT8vWNhCwU\/fW9dZw9fClhfWC8e4f3XVkm\/StqafILiNlljyYF76CPprX\/u\/09cMOush7nOdfeLHjs11cnAdHxL66SPEIOhRxHW8MA\/rsga5M4I5vO8tPmIRYvLCHPiIQxjnfeZD6Ps4Ylynpc8ccOI+X\/+T5ww2fMQwJ9\/o0EfoY0dgjg0doc9Y1uWbSb7pwT5M\/uAPTtr737\/+Y3DmY40wb+ZkDnzYaRHs5IzOmgg2+nFs6TcBPg4bOuOJW1m5KPUTj71ZMv9sOlewKc4U6hAdfeyhbd46353Pe415zs+LbyPz88LfyHgf6y9Y73ehjSKz7duzd3e+n8VFs1\/+5fUPQfeNanmI5f0wQ\/L2sM84JLTNUx+XzzzXHjU35zSK1kdE0fp+1b2ic0CeRXnMmuvOnc9a0TeP4drbtm0zPqH6rd86Hpor0X2\/119f\/dqVbHDMIp0r2GP2e5bb6IBnAAALhklEQVSbjw\/nIR\/96HF7\/euPrVuPA\/6+9z2cfHz0sK2sPJz6sRGH5F+ovEhDG\/11EyYdxt9zzxfTn4Em3dP\/wljmYH53oiOvfe36AkEc45iTWHTiwr7HfP7zx9ICQww24hiDYHv72\/8p3V\/e5344EIdgC9egzzj30UdnDca5jxYfgo8WITb0+dz4EGI9hhY\/rQsxLthcp6Xv8cxDPgg2+twLBJvHMwYdQX\/rWx+3\/JnbtGkt+Yjyq8aDiLkYf9NNjzLktDAWH+1pY6J439vElP4jlpzSzk8u2PJxuLCFzBjngs9zQXc7LX3m9DlCHR92F+IRX4cWcT\/xCDZaYlmXOUMbbPC5nZZ4Ynwu\/Iyljw+dOOz0Xf+933vcYB\/ejxe\/eM14F+lxjGWeUPDRp2UuWvrcP+7\/u971zeRHKme\/LokhNsyVPuI+1vvMZz6bjidfYu+558H09UQcNmKHCX7i8KMzHp0WQUc8htbtl1762FlnM2QT6u985zEjV8Yzn8+Bji3kgo49L8Rho0V8jjBvxiLuy+K+VzrPMGfWarNspGC3ed+2e\/du27t3b\/KzoQO2sLAwc7n++uvtF3\/xY8nPpv7eXvnKrHhdccU\/JIfsfjt58nP2jW98LvXv2PF\/7Mc\/vi2Va6550HjBvPrVX00+Qvzj5AX7oF155Z3GuDvu+Pvkgf5wOt\/Bgw+mxRL7wYMP2de+9o924423pz4ONsUSYcwLXvCp5Ge7\/2S\/\/dtftT\/7s79Lfv59XZoPuV188Z3pPLwAeGidf\/477Dd\/83q77LK7k5+v3Zr8\/OuuNMeXvvQv7MILH7LrrvtYmg9rffrTd9mhQx+1I0f+IZ2fPWJnDeK+\/vV\/tPPOuyv95oS1DiY5I+eee5f993\/\/v5TB25Oi\/spXfi7NizV50DF2y5avpg8CfsGKFy7z8sA6cOCY\/fzPfzPN4aqr7krXZb13v\/uWNE8e5KzBOBha8nVTUvQOJoxgwVzs9fHHP3g6b+bFzz6YCyEv8kmG2\/e+972UxQc+8H9TVpde+rh5PEwfeui29J5yH8lv06Zv27e+tZrmyZrkQv5XX\/329P6wzvvf\/yeF5+0tb3mLveIV\/5Deb84I9w4m5O25sDY58kCH+VVX3Zly8LXIe1NS\/InhfBDD+EuTBzH7YvzCwq3pOTuY3BPyIj+4cz\/9\/LAX7tt5yT28+up32C\/8wkPJebgzPUfuw39Vch9YkzyZh7PCuoxZSNZ585v\/2LjncH\/uubvtv\/7rT9J7e1NyX8iHNblv+GH1ilf8if30T7\/aGE\/u8L7ppsfSc3RVshZr\/3LySQTsDyb3lTlo2S95cN+vSpiQG2NvTF4XMIclfl5r2K5IXovE\/Ou\/Ft+L\/\/zPhTQH7ivnlRyZe\/\/+67gV6euBNeHG+qxx8cVX2aOP\/p59+MM3Jx\/pvzp9XWJnHc7oddd91H7u576ZsH88vWfkcTC5B8SwX\/bCnF\/4wlF75JFHkm8kvp08n47ZoeR19sEPfjB9RuBnPvLhvnKv4QErfLQwui55rTI3+7388lvTM0X7P\/9zWzLvWvoMIP4jH\/loep5Zg3zgVPZ5+IY3vMG+\/OUb7FWvujudk\/w5h8yDfjC5P\/jIiXPBGTmY7BeA5M0zhTjy4N7fk3xT8ork\/tPn3HG\/bkrOyRXJvTqYzOWx8Hrve68vfA2Ny\/3AgQPpc5\/nP3m0TXpdsI8cOZK8012Zi\/z5n\/958iK5InnYX2Af\/OB2+9M\/\/V5SCC+za665JrVfccUVabu0tGQu\/ILGhz\/8bPrLQsSdf\/75acy73707neflL\/9fafumN22zu+\/+pmG\/KPlZKLHMx\/jDh79t\/DIHOvKFL\/y8EYNceeXFyYPvxjQP4hHmWV215IW3lvqWknyw33HHa+zQof3p+sz\/iU9sTQrlxWkfP7bLLrvs9NzoxF+ZrIEsJfMQ43sgZ+SBB\/Yawhzkxx6uTMbQZ3\/oDz1khn755b+QcmNe5kHYB7HY0Imnj+zff46xBmNh8OUvr9nrXvejNEfWes97zrFPfep7dvToH6U2xvuc6OwBYewHPvAr9rGPraVcYIEdVsvLdnosY7jP7HMp2a\/PRS74WJNc6L\/tbVen9w77hz70oaFnjnUQxiCM5+ywF4Tx+FmLc3XzzfvSeX0t8ibuyoQpcQjjOVf8Ag\/jyYU+rNAPJffZ4\/38MD\/rszeEeViLMe7Dz1ha8qQllntKy9y03O9Pf\/rx5GPUqw1exPp9Yd79+88x\/Njf974bDT6Mw0e+rIccSvKkRbCHQiz+5WUz8iQX\/Nho6SNLyX2iJZ6WtVZWzn4GYCcH4mn37z8nfV2SH3xhf2XCGGF+5G\/\/9s\/S3BnLnIzDfv\/9L0lfN294w0XJ+Xs6nQf7oWQ\/8EQnlnNHXks\/yRE7e8UGN+zo5L1\/\/znpa4NcOJeeDy1jrkxyY+6lZC7m9nvC2WcMcczFGr7uoSQfYsm9rLBXeDMn83EOmYeWufGxFvNyf8mJWIRnCnFXJrniI++lJN8rkz7xzMM5YV7mQvDRsm7ZHPNxPPetpV\/Pa2ne49MuEcF3Wfv27bO+y8tetq2RDF7ykpdsKK9f\/\/Xd68YzX9426t7\/0R+tHz8qdt6+SfIelgv7H+aTvfxzYBb3oou8y56vuvnx3LeWfvW6YLf0niltERABERCBHhLoXMEu+o3wPXv2pP+JV4Pur1IRAREQAREQgYkIdK5gs3sv2hRqhD52iQiIgAiIgAi0lUAnCzY3gyLtQl8yAQGFioAIiIAINI5AZwt240grIREQAREQARHYAIFeFGw+Fi+SkFvoD+3SG0dACZUgwHkeFobPJR\/jdtq8T30RaAOBorOLrUjC\/YT+0N4kvRcFG+D+8XjYYke4UaGdPnaJCLSRwKjzi2\/YWR\/layMH5dw\/ApzhYbsOz73rHss4t9HSd1+T2t4U7GHQuTHcoNBPH3toky4CpQjUHDTq3OLjbIcp0seOoBf5Qpt0EWgqAc7wNLkxri1nvzcFm5viMs1N1RgRaAMBHjxIG3JVjiIwSwKce2TYnP78px0W03R7bwo2N9KlzTes6QdK+TWagJITgd4S8Oc\/bVtrQC8KNjcoPKX023rDwn1IFwEREAERGE+AZ34YRb+NNaAXBTu8UdJFQAQaSkBpiYAIjCTQi4Ldxu+kRt41OUVABERABEoT6EoN6EXBzt9Vbh4fiWCnpY\/uQh+799WKQBcIcKY52+Fe6GNH0It8oa3HurbeIQKcdc48W6Klj+5CH7v3m9L2omADnhvgQj+8AfTdR0s\/9EsXga4Q4Gxzxl3o+97Q3U5L331qRaDNBDjLnGkX+uF+6LuPln7ob4rei4INbG6AC\/28uI8271NfBNpGYNQ5xueS35fbafM+9RtMQKmdJjDs7GJ3OR0cKO6jDcyNUntTsBtFXcmIgAiIgAiIwIQEVLAnBKZwERABERCBiQgoeEYEVLBnBFLTiIAIiIAIiMA8Cahgz5Ou5haBFhHgl21alK5SFYHZEGjRLCrYLbpZSlUE5kWAYt3kX7aZ1741rwi0iYAKdpvulnIVAREQARHoE4F1e1XBXodDHRHoHwG9u+7fPdeO20lABbud901Zi4AIiIAI9IxA4wp2z\/hruyIwEwK8S2YiWgQdQUfQywrxRVJ2vOJEQATmQ0AFez5cNasIVE6AIssvjiHoCDqCXpQQdvyhj34o+OjTSkRABOojoII9EXsFi0BzCeSLar4\/TeZFBX2aeTRGBERg4wRUsDfOUDOIQCsJjCvG4\/yt3LSSFoEWE1DBbvHNy6euvgjMioCK9axIah4RmB0BFezZsdRMItAaAqMK8ihfazaoREWggwRUsDt4U5u5JWXVJgIU7by0KX\/lKgJdJKCC3cW7qj31jkD+F8xG9SnEeX8IDF+RhDHSRUAEqieggl09c63YQAJ9Soli3Kf9aq8i0BUC\/x8AAP\/\/1pmiOwAAAAZJREFUAwDwrB4y\/PC3kgAAAABJRU5ErkJggg==","height":237,"width":394}}
%---
%[output:9c31587c]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]  --- A08 | Section 4: Accurate Mass Database Search  [Analytics L3] ---\n","truncated":false}}
%---
%[output:7dcbac68]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]  Unknown substance 1 (High-resolution [M+H]+ = 136.1122, Neutral = 135.1049 Da):\n","truncated":false}}
%---
%[output:29e6414f]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]    Low Resolution (+-0.50 Da): 1 candidates\n","truncated":false}}
%---
%[output:5fa2fda3]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]    High Resolution (5.0 ppm): 1 candidates\n","truncated":false}}
%---
%[output:1b890286]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]  Unknown substance 2 (High-resolution [M+H]+ = 331.0150, Neutral = 330.0077 Da):\n","truncated":false}}
%---
%[output:53229f52]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]    Low Resolution (+-0.50 Da): 1 candidates\n","truncated":false}}
%---
%[output:4e4fd212]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]    High Resolution (5.0 ppm): 1 candidates\n","truncated":false}}
%---
%[output:8bb7d389]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]  Unknown substance 3 (High-resolution [M+H]+ = 286.1441, Neutral = 285.1369 Da):\n","truncated":false}}
%---
%[output:725e94cf]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]    Low Resolution (+-0.50 Da): 2 candidates\n","truncated":false}}
%---
%[output:9ff42df6]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]    High Resolution (5.0 ppm): 1 candidates\n","truncated":false}}
%---
%[output:0b035b93]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]  Unknown substance 4 (High-resolution [M+H]+ = 288.1742, Neutral = 287.1669 Da):\n","truncated":false}}
%---
%[output:6961240c]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]    Low Resolution (+-0.50 Da): 1 candidates\n","truncated":false}}
%---
%[output:24fc1537]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]    High Resolution (5.0 ppm): 1 candidates\n","truncated":false}}
%---
%[output:4bf03f71]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]  Unknown substance 5 (High-resolution [M+H]+ = 195.0878, Neutral = 194.0805 Da):\n","truncated":false}}
%---
%[output:7beb5664]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]    Low Resolution (+-0.50 Da): 2 candidates\n","truncated":false}}
%---
%[output:6491ba3f]
%   data: {"dataType":"text","outputData":{"text":"[21:18:52][INFO]    High Resolution (5.0 ppm): 1 candidates\n","truncated":false}}
%---
%[output:0cfb32e5]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAewAAAEoCAYAAACaU8LCAAAQAElEQVR4AeydQcsdx5X3GxOwDTODFxaKwAPCi0l2hmwkXoIW+QKBF7yIMjgi8C6NFoMyX0CrEZNPIALxQh9ilmE2wd6YWb2z8EYCy0yG2QRjD4KZ\/B772KVS3Xu7+3bfPt39M\/6rqk6dqjr1q7p9dO99\/Pi1\/\/EfCUhAAhKQgATSE3it8x8JSEACEpCABNITyJ2w0+MzQAlIQAISkMBlCJiwL8PZVSQgAQlIQAJnETBhj8fnSAlIQAISkMDFCJiwL4bahSQgAQlIQALjCZiwx7PLPdLoJCABCUhgUwRM2Js6TjcjAQlIQAJbJWDC3urJ5t6X0UlAAhKQwEACJuyBwHSXgAQkIAEJLEHAhL0A9XfffbcrdYkQWK9cp2737Sv95qgTV0vnrsWcvec4w\/FS6xDiJddivbl1bD9lX1k\/FlNfv3IOxrRU+sxdZ\/1z15hijmMx9J2\/9qNd6tga9r1KwIT9KpNZLVzWzz77rCuFbdZFVzZ5ySbqMlrZIa443LhzZen9e\/VATzGp+2mXTKlje3VmLYcImLAPkZnBzuXkktZTY6Ovts\/ZZs0p5587fuKde42hPFrxEOfQeQb67959Ccas2TrvDIfRiot4M8QWMRBjKyZs9IWf5XECJuzjfCbr5VJyOSeb0IkkIAEJLEiA5xnPtVYI2Olv9dW2vn71uD22TdhJTr2+tFz4UnWY9GGjDNEOhY0ybGVZ2qmHSp+oR1+UYacsbWWdPhQ2StpjxfhjjOp58Q\/VfbTpoyxV22iXCt+w0S7r0aYMRT9l2KIMG2Uo+iixUZ4r5gmVc2Gr2y1b6RP12u87+7vvRvWln9M45P+d84BKPRftWuV0ZV9pH1JnjrH3j7H1WthCdV\/ZxqdsUw8bJQpb1KNNGaIvFLYosVOnDNGeUrCba+4p48w+lwk74QlxsbngpbDVoWJr+Ryy1+Npn\/Kt+1kPG2MRbVTXaeNHX4g29lPCrxZzlOPoxxaiHf3Uw05JO\/r6loxhbClsjA9bXaddCv\/wpaRd9lPHRl+INnaEjfIcMR\/zhGifM1+MZb56LtrY8Yk67RA2+o4Jn5YOjcE35qfEL0rqdT9t7MeET61yTsbSjy1EGzuiHvYosdGHqIedkjb2oWIsYhwlol6L+ekL0R7jU48Z2i7XJwY0dI69+5uwE94ALnafsKb2Y816zrqNTx\/xYqzH0sZ+ajx+pWp\/5qC\/tNPGXtqiTl\/U+5TM0xrTsh2arzUH47GXY7CV7SnrrFXPTxs765R12rXww6e2j233mQufliZas9c09fr1oBYXxmCvfaNNP3V8ok4b0cZOfWoxL\/OX89LGXtvKdt96ay7mxn5oDvpC+B7y0\/4qARP2q0wWsdQXl3apRYL6dtEyDurfmk8W+NY6OajhwIubecou2rXKfuvDCMASzoj6sNGvesc8zBV61WsdlthLGW3sqSyjP\/xbfeGztxIWrT0Hq1aftlcJmLBfZTKLZcjF5HLjX2qWoHpMek4sZfxlvceyJ13K+cr6yYE6TEYA7twPJqSkTT1EuxQ+0TdVyfzMG6I91dzH5mGdlmJM3Ud80XewXGkHe439UdJe6VbSh23CvuARcZG50PWS2Oir7dGmP+pjyiHjT\/me6o\/42E\/Lt2WLMX3LoXP3WbP0GTp\/K+7WHKyBveU\/h421WLOcmzb2sFHvYwv\/viVz9vU9x4912EPonLn6jmUt1q39wxZl3U+7NRZ\/7PSfEr6nfMp+5q3H0MZe+s1dZz3WrdfBRl9tt90mYMJuc5nNyuXkkpbCVi5I+1h\/6duqDxl\/yvdUf6xf+rVssR\/8on9IyTjmiDHRxhbCRj9l2ChpYy+Fjb4Q7WP9+B3zKcdGHX\/GhWhHX5+Scaf88GkpxrFm2U87+qYomY\/5Kcv5aGMvha30maperkF9qnnLeYi9nDva2ELYGEMZtiix0Yeoh52SNvaW6MMnRLv2wxb9dR\/tsh8\/2tgH6qg7c56aO3zwC2E7OrGdLxEwYb+E4zINLmmp1qp1P+3Sr27TV9qoh+q+uh1+lHVftOlD0aasRT8q7bRLlX2tOr4tO7a6j3YpfEK1nXb0RYkthI06ZYh2qbCXZfSHjXbUKWmHaJfCXrapl7ayTl8t+g+p9C19SnvU6Y96lC1b9NXlIV\/spepxdRvf2hbtsi\/qPPSp18LOOOyUpVq2vv31WNql6nkO9eF3qg+fUO1LO\/qixIbKdtQp6QvRLoW9bFNv2bCHWv19bPiUivks+xEwYffjpJcEJJCMAA9+knMt7MlCNZy5CexkfhP2Tg7abUpgiwRIzrW2uE\/3JAEI7DphP3v2rFMy8A54B7wD+7kDJL4La7LldpuweYE+ePCgu3PnjpKBd8A74B3YyR34xS9+cfVGbbIsesGJdp2w\/\/jHP3aPHj3qnjx5kkb379+\/Ov5scWViVMYir+F3V2YyK19Dc9Sz3jHi4rl\/9ZBd4R+zJOw1cbh161Z3+\/btNCIe+FFmiitrLHCS17D7K7NhvLj7MhvGLCuviItnxhq1+4Sd7dB++MMfdh988EFHmS22jPHASV7DTkZmw3jhLTMo9Je8+rMa4rnDhD0Ez+V9uei\/+tWvLr\/wSleU1\/CDk5nMhhMYNsI7NoxXX28Tdl9S+klAAhKQgAQWJGDCXhB+a2ltEpCABCQggRYBE3aLijYJSEACEpBAMgIm7GQHkjsco5OABCQggaUImLCXIu+6EpCABCQggQEE0idsfrF\/7Ic6iralBEoCmer8Jr2sev78eYeyxpcxLnihjLFljAlWaO7YMr3mLxFL6oRNcuYX+wMi6rSpY1MSyEiAh1TmX3v7s5\/9rPvlL3\/ZUfqrefv9amJYyawfK+7UpXit+deMjnl2pU7YYzbkGAksTYCEza8\/fPnXyw7\/dZhz\/MpI5\/QctnIH1v5rRsc8p0zYY6g5RgI9CPBrEPm1lmrYr7WUl7z63AFeXz1ehptyMWFv6jjdTGYCvPPOqiFxZWa8h9j28pXgXvY55M6mTtjxfTUHR52NlXXaSgJrIEBCzPy9Nt879tXevjfMdL\/29Pzjmc9+M\/FfOpbUCRs4HBqijso6bSWBNRAgYfO99pc\/+XX355\/+ZrX66sc\/79hHyXz+h+o3q7FOS9\/0zvMn642d+ZyxY9ccMo740KEx9LVU+kd\/aavr4VOWtY\/tfgTSJ2y2EQcddUolgTUSePH2j7q1a0nu\/IW9Fs+HJWNi7UvEwBrsnfVq0VfbjrXxZy5E\/ZAv\/bXCl3HRRz3srTL8ojzlH3Pg39c3xmy5TJ+wOSwOLQ6BOrZoW0pAAhJYksChtXlWHepb0s7zs4yNOrYhMeHPuBhDHVu0T5VD\/U\/Nt5f+1AmbC8DB7uUw3KcEtkiA13Gp2CO2U\/Xo71syZ6lyXB87PuWYqLfsYSvLqDOurtMO0Y9oR0kd0c4i4gnNFVPMH+Vc62xh3tQJewuA3YME9kyAhzB\/6S6FbSwTxtZibubDTr0UtlbfITtjo49xfcQY\/CgR9VLMh70UtvChfqgvfFol40L0t+rvvvtuF3ZK\/IZoTFxD5se3XIP6mDiZZw8yYe\/hlN2jBBIT4AHd90GNX4gtUacMMVepsOOHvWxHfemS2MbEwLgQ41v1sEWJX18xpvSlXTIs+86tM2\/o3Lm2PD51wm5dEA4V+5YPxb1JQAKnCfAc4HlQemKrFf3Y8Udhs1yWAKtzHpxNCJtqE0idsAmZQ+RAqVPSpq4kIIF1E+C1XL6m6\/YUu2N+5omSNVC06Rujc8ePWXOKuFm3noe9YKOvFPa6HX6UZT91bKX\/sXpf\/75+x9baUl\/6hA1sLkKItpLAWgm89uV\/dq99+afV6hB3Hqy18OV1W9ux0TeFmIv5mSvqtEPYWn2H7IyLPsaFsNEXoh19lLSjj3apsi98sJU+59aHzof\/oViwE0\/pg4029hBt7Ih62FslPqVKf+plX2v8MrZ8q6ZO2BxiC9khe8tXmwQyEHjnnXc6fvfxX\/3rP3V\/8y\/\/uFoRP\/tgP8GVB25Lh\/rDTsk4ylDdDjvlob7STr0U40J97PiEP2XZph6q+6JNP3VU12mH6Ee0KUu1bNFP31TPP+ZCMXeUpY16KPrL8lhf+IVPWUZflGUfdezsM+q0VdelTtgekAS2QoAEt5X\/exf72Mq5rHEfe0lie9nnkDtYJuwh42b15W9WiEUoa3mQkFFrJEDi3oLWyN6YJbB2AikTNgkZAZeyFnYlgTUR4HeJ+z\/\/WNOJ5Y+VNzJZoswUSxYmc8SRMmE3NzrAyOUpNWCorhKYhQAJm\/9pxoOb\/9H98999vlp9cOO\/Xvmff3T+c3ECPN94I3PxhQ8sSCzEdKBb80QEUidsLsHQfXJpGFcK29B59JfAHATe++uvurWr5nLo9VXay3o9nvapfnxC+LYU\/XOUrDd23nPG9lmT+VvqM1afdRFInbBXhNJQJSCBCxIo\/0IedZLWBUNoLnWJGFiDPZcB0K5V9l+izvrEdom19rpG+oTNBWjp0IFxaQ71tex8VBl68eJFp2Rw7h1o3bO924a+LrfEa897v8Q5nnq9xvOd8hLxzLlG6oRNouayIyDUJbZDYmwoxrV879692925c+dKDx8+7J4+fbqouFTPnz+fNoaF9zQn04y8OL\/WXduzjddi7J86oh0ldUQ7RHusYo4oy3nCRnnIXveFX8setrKMOuPqOu0Q\/Yh2lNQR7b7CP3RsTPhEGb7RjrK0U2\/ZsdG3tD7\/\/POjz8vHjx9fPd95zvO8Xzrec9ZPnbBbGyP59rko+IWO+fPflD558qRD9+7d627cuLGorl+\/3l27dm3RGJZmMGT9jLzeeuut1tXdtI3XWK3WhvEpX5fUw6\/sw047+uqSvlqMwQ879VLYWn2H7IyNPsb1EWPwo0TUSzEf9lLYwof6ob7wOVQOGXfI95C9jIs6Cl\/qh2K6lP3U85LnOs93dP\/+\/UuFNcs6q0vYU1Pgtzbdvn27Qzdv3uzeeOONRfXmm292r7\/++qIxXJjBWXvNymvqe5p9vniAl2UrZvp5yCPqpU\/dLvvqOr4h+qhThpi\/VNjxw162o750SWxjYqjH0S73OGbOcgzzHWuXfUvUTz2veK7zfEc875eIcao1N5ewp7yoU0F2HglI4BsCvD5JAOgby\/l\/MhfzljNhqxX92PFHYbOUwBoIpE7Y8cICZNR5kVHHpiQwO4GJF\/ji6x90z\/97vZoCB6\/h0BTz9ZmD9fCLkmcIijZ9Y3Tu+DFr1nHXMdDGZ8zc54xZat1zYl7b2NQJG5jlxaOOsB8S\/VycUtgO+WuXwCUI8OtI+TjuH\/79Rvf3\/\/a3qxXxsw\/2M5Ybr8cQr9Ox85TjmC\/mijrtEDb8KcNGSbtlL\/voD+FPX4h29FHSjj7apcq+8MFW+oypM0fMR0n72Dz4hErfsFGW9mNz2XdZAukT9hgcXLZSY+ZwjASmJECCK3\/AkR+AmUBXPyx56XnYR8mG11rZjnppjzol7fx2XwAAEABJREFUCSFEG\/8oqYdaNvr62PEpxbhQHzs+4U9ZtqmH6r5o008d1XXaIfoRbcpSLVv00wfDso0Nhe1QiU+o9AkbZdjLOrZDbWKp+\/BX0xJIlbA59L6aFoOzSWB+AiRtfvBl7WIf59DiwR46Z569j4VhFgaZYsnCZI44UiVsDr0Wmy5t0aZUEshMgP9GfNd69qxz\/7kY\/OEPf2ieySF75vPL\/NqfK7ZUCbveJO+2SdalnTb20mZdApkI8A6U73n5JQ38sgb1zS8mkoMcprwDvL54nfF6y\/T6nzOW1Al7zo07twTmIsADhO95L\/3dct\/1fve733W\/\/e1vu48++miR78D7xjmz36C9y+ybXy7V90wuxYvX2Vyv44zzmrAznooxrZ4ASTvzd9Xvvffe1S8LyhxjtthkdnvQnbkEL15nq39YDNhA6oQdH3\/zEXgp7AP2qKsEJCCBdREwWgk0CKRO2MRLcq6FXUlAAhKQgAT2RCBVwo530XEA0W6V4WMpAQlIQAIXJeBiCxFIlbDjnXSwiHarDB9LCUhAAhKQwB4IpErYewDuHiUgAQlIYEYCG57ahL3hw3VrEpCABCSwHQKpEnbru+pDtu0cgTuRgAQkIIGdEDhrm6kS9qHvqks7u6VNqSQgAQlIQAJ7IZAqYdfQeXddJ2fa2Gtf2xKQgAQkIIEtE5g9YW8ZnnuTgAQkIAEJXIqACftSpF1HAhKQgAQkcAaB1Ak7Pv7mI\/BS2M\/YczHUqgQkIAEJSGAdBFInbBCSnGthVxKQgAQkIIE9EUifsMt31mV9D4fkHiUgAQlIQAJBIHXCJkHX766jHRuwlIAEJCABCeyBQOqEvYcDWO8ejVwCEpCABC5JwIR9SdquJQEJSEACEhhJIHXC5uNvPhYfuTeH7ZiAW5eABCSwNQKpE3Yka8paWzsI9yMBCUhAAhI4RiB1wuYd9iEd25R9EshNwOgkIAEJDCeQOmEP344jJCABCUhAAtskYMLe5rm6KwmMJuBACUggJ4H0Cbv+7jraOXEalQQkIAEJSGAeAqkTNsk5vsNm+9TLkrqSgAT2RMC9SmC\/BFIn7NaxkLRJ5K0+bRKQgAQkIIGtElhdwt7qQbgvCUhg\/QTcgQTmJGDCnpOuc0tAAhKQgAQmIpA6YZcff0edj8OpT7R\/p5GABCSwEwJuc+0EUids4JbJmTrCriQgAQlIQAJ7IpA+Ye\/pMNyrBCQggb0ScN+nCaRM2HzsjcrwaaPSZl0CEpCABCSwFwLpEjZJmY+9UXkItBH9pd26BCQgAQlIYF4COWZPlbBJxiTlY2jox++Yj30SkIAEJCCBrRFIlbC3Btf9SEACEpCABKYicChhTzW\/80hAAhKQgAQkMAEBE\/YEEJ1CAhKQgAQkMDeBVAm7z\/fTfH\/92Wefzc3F+SUgAQlIQAKpCKRK2JA5lrRN1hBSEpCABCSwRwLpEjaHEEmbBF0KO\/3JZXgSkIAEJCCByQmkTNjskuRcC7uSgAQkIAEJ7JFA2oS9x8O4yJ5dRAISkIAEVknAhL3KYzNoCUhAAhLYGwET9t5OPPd+jU4CEpCABA4QMGEfAKNZAhKQgAQkkImACTvTaRhLbgJGJwEJSGBBAqtI2PGfdsGJOqWSgAQkIAEJ7IlA+oRNguY\/74pDoY4t2pYSkMAVAf+QgAQ2TiB1wiYxk6A3fgZuTwISkIAEJHCSQOqEfTJ6HSQggXUQMEoJSOBsAibssxE6gQQkIAEJSGB+AqkTNh+H87F4iYE29tJmXQISkMAZBBwqgVUQSJ2wIUhyJklTp6RN\/ZTwDZ3ytV8CEpCABCSQnUD6hA1AknSI9imRqMOfkvapMfZLQAISSEnAoCTwLYH0CbtMttTRt7E3C\/pJ0s1OjRKQgAQkIIGVEkidsMvkG3WSMfWV8jZsCUhAAlsi4F4uSCB1wh7DgYRejiO517ay\/9mzZ13oxYsXnZJBxjsQd\/TcMuPejCnPa+7c+8X4bOdJTKHy2b\/G+uYSdnkIp5I1vnfv3u3u3LlzpYcPH3ZPnz5dVFys58+fLxrD0gyGrL8HXh9\/\/HH34YcfXt3RuKtjy\/fff7\/75JNPOu\/YsNf5Xu7Z\/\/m\/\/+\/se5bqjv3lef748ePv9sTznuf+WrXZhN0nWXNojx496p48eXKle\/fudTdu3FhU169f765du7ZoDEszGLL+Hnh9\/fXX3aefftp9+ZNfd3\/+6W9G66sf\/\/xqnrfffts7NvB1vpd79oM\/\/f+z7lnGO8ZzPZ7x9+\/f57G\/WqVO2HyUTeJF1KFc1mm31Mcnxt26dau7ffv2lW7evNm98cYbi+rNN9\/sXn\/99UVjWJrBkPX3wov7+uLtH3Xninng6x0b9jr3nvW\/e9nuGM\/1eMbzvCe+JBocRuqEzW5I1Ig6Kuu0aw1J1vVY2xKQgAQkIIGsBNIn7DHgSNq1xszjGAlIQAISkEAWAhdN2GM2XSfeaB+ai3fgLR3y1y4BCUhAAhJYA4HUCZvkHMkXmNTLkrqSgAQkIAEJ7IFA6oTdOgCSNom81XeezdESkIAEJCCBvARWl7DzojQyCUhAAhKQwHwETNjzsZ10ZieTgAQkIIF9E0idsMuPv6POx+HU931s7l4CEpCABPZGIHXC5jDK5EwdYVeZCBiLBCQgAQnMTSB9wuYdda25oTi\/BCQgAQlIIBuB1AmbRM076lrYs4E0nrwEjEwCEpDAFgikTthbAOweJCABCUhAAlMQMGFPQdE5JDCagAMlIAEJ9COQOmHzUXj98Tdt7P22p5cEJCABCUhgGwRSJWyScS0wl7ZoUyoJSGBeAs4uAQnkIZAqYfPOua\/yIDQSCUhAAhKQwPwEUiXs+bfrChKQwHYIuBMJ7IuACXtf5+1uJSABCUhgpQTSJ+zy++uyvlLehi0BCeyEgNuUwNQEUidsEnR8p83GqZcldSUBCUhAAhLYA4HUCbt1ACRtEnmrT5sEJCABCfQhoM8aCawuYa8RsjFLQAISkIAEziVgwj6XoOMlIAEJSGBSAk7WJpA6YZcff0edj8Opt7ejVQISkIAEJLBNAqkTNsjL5EwdYVcSkIAEJCCByxNYbsXUCZt30y00h+wtX20SkIAEJCCBLRBInbC3ANg9SEACEpCABKYg0CdhT7HOoDl4B40YRFnLj8UhoyQgAQlIYE8EUiZsEjLiIChrYVcSkIAEJCCBPRFImbAHHYDOEpCABCQggR0QSJ2weWe9gzNwixKQgAQkIIGTBFIn7JPR53cwQglIQAISkMAkBNIn7PoHzqI9ye6dRAISkIAEJLASAqkTNsmZj8URPOsSmzqDgEMlIAEJSGA1BFIn7BZFkjaJvNWnTQISkIAEJLBVAqtL2Fs9CPf1CgENEpCABCRQEDBhFzCsSkACEpCABLISSJ2wy4+\/o87H4dSzAjWunRBwmxKQgAQuTCB1woZFmZypI+xKAhKQgAQksCcCaRM276Rr7elg3KsEziDgUAlIYIMEUiZsEjXvpGth3+AZuCUJSEACEpDASQLpEjZJmUTdihw7\/a0+bRKQwEoIGKYEJDCKQKqETTImKR\/bCf34HfOxTwISkIAEJLA1AqkS9tbguh8JSGB1BAxYAmkJmLDTHo2BSUACEpCABL4nYML+noU1CUhAArkJGN2uCaRL2Hw\/fUq7PjE3LwEJSEACuySQKmHzA2V9tcvTctMSkIAE8hIwspkJpErYM+\/V6SUgAQlIQAKrJWDCXu3RGbgEJCABCfQmsAFHE\/YGDtEtSEACEpDA9gmYsLd\/xu5QAhKQgARyE+gVXaqEzU+H94paJwlIQAISkMDOCKRK2DV7E3hNxLYEJCABCeyVwGIJe6\/A3bcEJCABCUhgDIFUCZv\/Bpt31SE2FPW6pE9JQAISkIAE9kIgVcIGOkk7VLfDTknffHJmCUhAAhKQQC4C6RJ2LjxGIwEJSEACEshBIHXC9p10+5JolYAEJCCB\/RFInbDjOMrvr8NmKQEJSEACEtgTgfQJm2TNO+0Q7T4H1Nevz1z6DCGgrwQkIAEJzEEgdcIm6ZKoy43Txl7a6vqp\/trftgQkIAEJSCA7gdQJeww8k\/UYavsZ404lIAEJrJXA5hI278BR3wN59uxZF3rx4kU3RjH+nLK17jnzxdjWvBlsEd85Zb2Pc+Yqx9bzLt3ue5f7+pX7Kfd9Tr2cM0v9nP3E2NZeou+csjXv0ra+96ePX72Xc1jF2HrOvu0YT9kn9sw+qRM2ibd+x0wb+1RQ79692925c+dKDx8+7J4+fTpIH3\/8cffhhx9ejY95xpTvv\/9+x1xcqufPn3effPLJpPMO3dec\/uxzKmZwmpIXZxdnMSeD9tztu8f+prrvzMN8CHZTnENGZlPeMebaw+uSO8H9mELMhbjn8Jvino19XT5+\/Pi75zPP+yn2t9QcqRM2UEjOJOkQbexT6dGjR92TJ0+udO\/eve7GjRuD9PXXX3effvpp9+VPft39+ae\/GaWvfvzzqzmuXbvWXb9+vaP86quvrmwPbv5H989\/9\/kofXDjv67mYL6h+5rTf0pmb7\/99mS84JyR2VtvvdVN+c8emE15x3j97OF1OeU9izvGcyTO4sFCzzKe6\/GMv3\/\/frfmf9InbOCSpEO0p9StW7e627dvX+nmzZvdG2+8MVjE8+LtH3XniDlY+8033+xef\/31qxiwvffXX3XniDmYN5uI6xxejGUO9jUlL1jHvMydRcQ0ldjTGGawOSRiY95MIibuyTliDva0p9clez5XMIs7Rp35Dt2dvnbmYK6h4rkez3ie98yzVq0iYa8VrnFLQAISkIAEpiJgwp6KpPNIQAILEnBpCWyfgAl7+2fsDiUgAQlIYAMENpuw+c57A+fjFiQggQ0QcAsSmIJA6oTNT4ZPsUnnkIAEJCABCaydQOqEvXa4xi8BCUggPwEjXAuB1Ambj7V9l72Wq2ScEpCABCQwJ4HUCXvOjTu3BCQgAQnkJ2CE3xNInbDj3TVlre+3YE0CEpCABCSwfQKpEzYfiR\/S9o\/GHUpAAhKQQG4Cl40udcK+LApXk4AEJCABCeQlsIqEHR+Hg5E6pZKABCQgAQnsicDQhH1xNiRoPhaPhalji7alBCQgAQlIYA8EUidsEjMJeg8H4R4lIAEJSEACxwikTtjHAm\/2aZSABCQgAQlslIAJe6MH67YkIAEJSGBbBFInbD4O52PxEjlt7KVtJXXDlIAEJCABCYwmkDphsyuSM0maOiVt6koCEpCABCSwJwLpEzaHQZIO0VYzEHBKCUhAAhJITSB9wuZdda3URA1OAhKQgAQkMAOB1AmbRB3vrMsS+wwsnDIvASOTgAQksHsCqRP27k9HABKQgAQkIIFvCZiwvwVhIYHRBBwoAQlI4AIEUidsPgavP\/6mjf0CbFxCAhKQgAQkkIZAqoRNMq4FqdIWbUolAQmcJKCDBCSwEQKpEjbvnPtqI\/zdhgQkIAEJSKAXgVQJu1fEOklAAsgOwf4AAAkHSURBVNsh4E4kIIHeBEzYvVHpKAEJSEACEliOQPqEXX5\/XdaXQ+bKEpDATgi4TQmkIpA6YZOgD32nnYqiwUhAAhKQgARmJpA6Yc+8d6eXgAQksF4CRr47AqkTNu+ueZe9u1NxwxKQgAQkIIGKQOqEHbGStGtFn6UEJCABCaQjYEAzEEidsEnSvMtuaQYWTikBCUhAAhJISyB1wk5LzcAkIAEJSGC9BFYaeeqEzTtr3mWvlK1hS0ACEpCABCYjkDphR7KmrDUZASeSgAQkIAEJ5CFwMJLUCZt32Id0cEd2SEACEpCABDZIIHXC3iBvtyQBCUhAAhIYRSBFwj4Uef0xeNk+NEa7BCQgAQlIYIsEUifs1sfhHAJ2SiUBCUhAAhLYC4HUCbt1CCRr3mm3+uaxOasEJCABCUhgeQKrS9jLIzMCCUhAAhKQwOUJmLAvz3zSFZ1MAhKQgAT2QSB1wuaj75b4WHwfx+MuJSABCUhAAt8QSJ2wScwtfRO6f+YnYIQSkIAEJDAVgdQJe6pNOo8EJCABCUhg7QRSJezWx9+HbGsHb\/zLEzACCUhAAmsikCphtz7+rm3AxUapJCABCUhAAnshkCphH4Me77RN1sco2bcdAu5EAhKQwMsEVpGwSdYkavRy+LYkIAEJSEAC+yCQOmGTqJGJeh+X0V2uh4CRSkAClyeQMmGTpBGJGl0eiytKQAISkIAEchFIl7BN1LkuiNFIYH0EjFgC2ySQKmGTrMFMeUr4KQlIQAISkMBeCKRK2Hz83Vd7OSD3KQEJbIuAu5HAWAKpEvbYTThOAhKQgAQksHUCJuytn7D7k4AEJNCbgI6ZCZiwM5+OsUlAAhKQgAS+JWDC\/haEhQQkIAEJ5Caw9+g2m7DLnzLf+yG7fwlIQAISWD+BTSZsknX50+a013JUz58\/737\/+993X3zxxVpCXjROeQ3HLzOZDScwbETcMcphI9fsPX\/sm0vYJGeSdYmONvbSlrXOBf\/oo4+yhpcuLv5iI69hxyKzYbzw9nUJhf7yjvVnNcRzcwl7yObxffbsWXeOmAO99uV\/dq99+adRYjwiDspSX3z9g+75f49TzMO8mRRxvTYBs5grynN4wTnm2Rqv1\/5yN2NvdSmz9us2OHEXoh7lOcxiDubNpIjrtQlfl+wv5s3ALGJZa3lOwl7rnq\/ifuedd7pbt251d+\/e7e7cuTNajGfCv\/rXf+r+5l\/+cZQYyxzMhag\/ePCAovuHf7\/R\/f2\/\/e0oMZZJmPOcPU49lniIi32fyyw4Rcmex\/JiHOOJjRin3vfY+YiFmM7hBWfGM0+wipI9s\/exYjzzEufYPU49jliIiT2z9zFiLHMwF6I+BbOMvOAfe2TfY3gxhrElJ+ZE2Nj30neMWHju8\/wnprVp1wn70aNH3ZMnT5QMvAPeAe\/ATu4Az\/21JeqId7sJO3Z4pORvWbdv3+6UDLwD3gHvwD7uAM\/9I2khddeuE3bqkzE4CUhAAhKQQEFgcwm79RPh\/IQ49mLfS1ddXwISkIAEJDCIwOYSNrsnOZOkQ7SxKwlIQAISkMBaCWwyYXMYJOkQbTWAgK4SkIAEJJCOwGYTdjrSBiQBCUhAAhI4g4AJ+wx4cwzlY\/w55t3QnK9sBWahVzo1vEQgOEX5UqeNowRkdhTPd53BqS6\/c7AymoAJezS66Qdywaefddszwiy++qCkve0dj98dbGBUCtv4GR0pgTaB8o5Fve2pdQgBE\/YQWjP6+uAcDhdmPAyGj5xxhFNvkoB3bZPHurpNmbCTHBmJByUJxzA2SMD7tcFDTbol\/oITShriKsMyYa\/y2AwaAnUC4gFR2\/BT3xG4qsApJK8rJEf\/gJWcjiJ6pRNeIfi94qBhFAET9ihsDspGgIcCD4hscWWMB04huGWMMUtM8IFVlnjWEEfNizYc1xB79hhN2NlPyPhOEuBhwEPhpKMOuQkkjY77FSJE6pRKApcmYMK+NHHXm5QAD0+TdT+ksOrnqVcQ4G6Vwk6bUrUJeM\/aXKawmrCnoOgcixDgweDDcxH0e1zUPY8k4Ot0JLjGMBN2A4qm9RDgYVBrPdFfNlL+clOzwnbZKFxt6wS4U+U9o731PV9qfybsS5HuuY6Xuyeov7jBqqW\/dPnvAQI1rwNumg8QgN+BrmXNyVaHUyhZaKsOx4S96uMzeAlIQAIS2AsBE\/ZeTtp9SkACEliGgKtORMCEPRFIp5GABCQgAQnMScCEPSdd55aABCQggdwEVhSdCXtFh2WolyXAT7oeWvFYXz1miG89NkOb+EstHROxLB2D60tgCQIm7CWou6YEVkKA5Bg\/7RsltpWEb5gSWDuBl+I3Yb+Ew4YEJBAESMwk6WhHiY2+aFtKQAKXIWDCvgxnV5HApgiQtDe1ITcjgRUQSJewV8DMECXwEoF4t0kZesmhauATpqhThqIvyrBTho2y1W7ZSl\/6Q9iPiaTcxzd8oqznxI6NElFH1BH1UmGjDJX9dT18KOs+2xLYEgET9pZO070sRoBkQYIL0W4Fgx2fsi9s2BHt6KeOLUQ7+oaWjI15KGmfmgM\/hG+oHION\/lLYSh\/q2MKHOirb+JQq+\/GjXfZHHTv9IdrRZymBrREwYQ86UZ0l0CZAwmj3fG8lmbT8WjZGtfzxxU5\/WaddCz98wl7Ww9a3ZGyIeWMctqgfK2u\/ul2PPdWPP3HUfrSx068ksDUCJuytnaj7SUmAJDJnMon5p1yDOVsw6zXwK9UaM5etXDfqc63lvBJYmoAJe+kTmHB9p8pLgCSXN7rxkZEk2Vup8bMNH1muW9aHz+QICeQnYMLOf0ZGuBABEgAJqV4eG321vU+bcYwf68tY5ojx1PvYwn9IGXPXY+r1yn76yvbYep95jsU3dl3HSSAzARN25tPZVGzr3EwkBRJICNs5u2E8c\/WZI3zxR7T7jJvKh\/VYtxS2mJ\/6ob7wGVP2nbf2IxZsY9Z0jASyEzBhZz8h41ucAAmgVB0Qfcdsx\/qP9cWc+ITCVpb0lW3qta1ut3ywtcTYUrVP3Ue79BnajrGMC4WNEhtliHapsFtKYGsETNhbO1H3M4qAgyQgAQlkJ\/C\/AAAA\/\/9mXmPiAAAABklEQVQDABP2VHtch6yZAAAAAElFTkSuQmCC","height":237,"width":394}}
%---
%[output:6a8239e9]
%   data: {"dataType":"text","outputData":{"text":"[21:18:53][INFO]  --- A08 | Section 5: Isotope Pattern Scoring  [Analytics L3] ---\n","truncated":false}}
%---
%[output:75c7f4a6]
%   data: {"dataType":"text","outputData":{"text":"[21:18:53][INFO]  Unknown 1 -- High-resolution candidates ranked by isotope score:\n","truncated":false}}
%---
%[output:701fdc6c]
%   data: {"dataType":"text","outputData":{"text":"        Name         Formula     ExactMass    ppmError    IsoScore\n    _____________    ________    _________    ________    ________\n\n    \"AMPHETAMINE\"    \"C9H13N\"      135.1      0.68123     0.99957 \n\n","truncated":false}}
%---
%[output:4e58e823]
%   data: {"dataType":"text","outputData":{"text":"[21:18:53][INFO]  --- A08 | Section 6: Complete Identification Workflow  [Analytics L3] ---\n","truncated":false}}
%---
%[output:71c00218]
%   data: {"dataType":"text","outputData":{"text":"[21:18:53][INFO]  Executing complete identification workflow...\n","truncated":false}}
%---
%[output:3f094f67]
%   data: {"dataType":"text","outputData":{"text":"[21:18:53][INFO]  Unknown substance 1:  True=AMPHETAMINE  Predicted=AMPHETAMINE  [ Correct]\n[21:18:53][INFO]  Unknown substance 2:  True=FUROSEMIDE  Predicted=FUROSEMIDE  [ Correct]\n[21:18:53][INFO]  Unknown substance 3:  True=MORPHINE  Predicted=MORPHINE  [ Correct]\n[21:18:53][INFO]  Unknown substance 4:  True=CYPROHEPTADINE  Predicted=CYPROHEPTADINE  [ Correct]\n[21:18:53][INFO]  Unknown substance 5:  True=CAFFEINE  Predicted=CAFFEINE  [ Correct]\n","truncated":false}}
%---
%[output:7d7b5589]
%   data: {"dataType":"text","outputData":{"text":"    Unknown         TrueID            Predicted        ppmError    IsoScore    Correct\n    _______    ________________    ________________    ________    ________    _______\n\n       1       \"AMPHETAMINE\"       \"AMPHETAMINE\"        0.68123    0.99957      true  \n       2       \"FUROSEMIDE\"        \"FUROSEMIDE\"        0.053748    0.99982      true  \n       3       \"MORPHINE\"          \"MORPHINE\"            1.2995    0.99966      true  \n       4       \"CYPROHEPTADINE\"    \"CYPROHEPTADINE\"      1.8051    0.99963      true  \n       5       \"CAFFEINE\"          \"CAFFEINE\"           0.61124    0.99948      true  \n\n","truncated":false}}
%---
%[output:686bed15]
%   data: {"dataType":"text","outputData":{"text":"[21:18:53][INFO]  Identification accuracy: 5 \/ 5 = 100%\n","truncated":false}}
%---
