%[text] # A08: Mass Spectrometry × Cheminformatics
%[text] EasyMolKit Analytics — Layer 3
%[text]
%[text] Five drug samples have arrived at the pharmaceutical QC lab. The labels are damaged by water and unreadable.
%[text] Can you identify "what drug is this" using only the available LC-MS data, without HPLC?
%[text] Mass spectrometry (MS) measures the "exact mass" of molecules at the ppm level.
%[text] By utilizing isotope patterns (spectral fingerprints of M, M+1, M+2), you can narrow down candidates with similar masses to one.
%[text] In this script, you will experience the process of identifying compounds through three stages: low-resolution search, high-precision mass, and isotope confirmation.
%[text]
%[text] **Story**
%[text]
%[text] A pharmaceutical QC analyst received a batch of five drug samples from the warehouse.
%[text] The labels are damaged by water and unreadable.
%[text] Each sample has already been measured with the lab's LC-MS equipment, and raw mass spectra have been obtained.
%[text] The analyst must confirm the identification of each sample before releasing the batch.
%[text]
%[text] The identification process unfolds in three stages with increasing confidence.
%[text]
%[text] - **Stage 1 — Low-resolution search** (unit mass tolerance 0.5 Da): In a real commercial database (with tens of thousands of entries), dozens of compounds match, but in this exercise's 200 compound DB, it will be 1-2.
%[text] - **Stage 2 — High-resolution exact mass** (5 ppm tolerance): Only a few candidates survive. Promising, but not definitive.
%[text] - **Stage 3 — Isotope pattern confirmation**: The characteristic fingerprints of M, M+1, M+2 narrow each unknown to a single compound with high confidence.
%[text]
%[text] **Exercise Content**
%[text]
%[text] 1. Use `emk.descriptor.calculate()` to build a reference table of exact masses from 200 FDA-approved drugs.
%[text] 2. Simulate five realistic ESI-MS spectra from five drugs (Gaussian peak shapes, isotope clusters, random noise).
%[text] 3. Use the Signal Processing Toolbox (`smoothdata` + `findpeaks`) to detect peaks in the raw spectra.
%[text] 4. Search the reference table with mass tolerance and compare the number of candidates at low and high resolution.
%[text] 5. Calculate theoretical isotope patterns from molecular formulas (relative intensities of M, M+1, M+2 from 13C / 34S / 37Cl abundance).
%[text] 6. Score each candidate by cosine similarity between its isotope pattern and the observed spectrum, creating a final ranked list.
%[text]
%[text] **Learning Objectives**
%[text]
%[text] - Explain why exact mass measurements reduce ambiguity in compound ID.
%[text] - Use `findpeaks` and `smoothdata` (Signal Processing Toolbox) on spectral data.
%[text] - Calculate ppm mass error and set appropriate search tolerances.
%[text] - Derive M+1 and M+2 isotope intensities from natural element abundances.
%[text] - Score candidates by cosine similarity and build a ranked hit list.
%[text] - Understand MS compound ID in connection with the fingerprint-based approach of S05.
%[text]
%[text] **Prerequisites**
%[text]
%[text] - Completion of A03 (QSAR Regression) — Understanding the basics of `emk.descriptor.calculate()`.
%[text] - Completion of S05 (Unknown Substance ID) — Understanding the concept of similarity-based identification.
%[text] - Ability to use the Signal Processing Toolbox (`findpeaks`, `smoothdata`).
%[text] - Ability to use the Statistics and Machine Learning Toolbox (`corrcoef`).
%[text] - No internet connection required.
%[text]
%[text] **Operating Environment**
%[text]
%[text] Sections 3 and 4 require the Signal Processing Toolbox.
%[text] Section 5 uses the Statistics and ML Toolbox (`corrcoef`).
%[text] Both toolboxes are included in MATLAB Online Basic (free tier).
%[text]
%[text] Estimated time required: 35–50 minutes
%[text]
%[text] **Data:**
%[text]
%[text] `data/list/fda_drugs.csv` — 200 FDA-approved drugs (ChEMBL, CC-BY-SA 3.0)
%[text] Columns: ChEMBLID, Name, SMILES, MolecularWeight, ALogP, HBondDonors, HBondAcceptors, TPSA, RotatableBonds
%[text]
%[text] **References**
%[text]
%[text] Gross JH (2011) *Mass Spectrometry: A Textbook*, 2nd ed. Springer.
%[text] ISBN 978-3-642-10709-2 — Single isotope mass, isotope patterns, ESI-MS
%[text]
%[text] Kind T & Fiehn O (2006) Metabolomic database annotations via query of elemental compositions.
%[text] *BMC Bioinformatics* 7:234. doi:10.1186/1471-2105-7-234 — Practical limits of mass accuracy in compound identification
%[text]
%[text] Claesen J et al. (2012) Efficient method for isotopic distribution calculation.
%[text] *J Am Soc Mass Spectrom* 23:753-763. doi:10.1007/s13361-011-0326-2
%[text]
%[text] Stein SE & Scott DR (1994) Optimization and testing of mass spectral library search algorithms.
%[text] *J Am Soc Mass Spectrom* 5:859-866. doi:10.1016/1044-0305(94)87009-8 — Cosine similarity scoring
%[text]
%[text] Meija J et al. (2016) Atomic weights of the elements 2013 (IUPAC Technical Report).
%[text] *Pure Appl. Chem.* 88:265-291. doi:10.1515/pac-2015-0305 — Natural isotope abundances: 13C 1.103%, 34S 4.25%, 37Cl 24.23%
%[text]
%[text] How to run: Execute each section with Ctrl+Enter
%%
%[text] ## Section 0: Setup

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

%[text] Warm up Python and RDKit processes before main processing
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;
logSection("A08", "Section 0: Setup", "Analytics L3");
%%
%[text] ## Section 1: Building a Precise Mass Reference Table
%[text]
%[text] Setup is complete. First, we will build a precise mass reference table from 200 FDA-approved drugs.
%[text] This will serve as the "dictionary" for database searches.
%[text]
%[text] ### Concept: Monoisotopic Mass vs. Average Molecular Weight
%[text]
%[text] The periodic table lists two types of "mass" for each element:
%[text]
%[text] Average Molecular Weight (MolWt):
%[text] A weighted average of all stable isotopes using natural abundance.
%[text] Example: Carbon = 12.011 g/mol (12C is 98.9%, 13C is 1.1%)
%[text] Used by chemists for stoichiometry and solution preparation.
%[text]
%[text] Exact Monoisotopic Mass (ExactMolWt):
%[text] The mass of a molecule composed only of the most abundant isotope of each element.
%[text] Example: 12C, 1H, 14N, 16O, 32S, 35Cl, etc.
%[text] Example: Carbon = exactly 12.000000 Da (the standard for mass scale!)
%[text] Used in mass spectrometry because it measures the mass-to-charge ratio (m/z) of individual ions.
%[text]
%[text] Why is this distinction important?
%[text] Aspirin (C9H8O4): MolWt = 180.16 g/mol, ExactMolWt = 180.0423 Da
%[text] This difference (0.12 Da) is large enough to cause database mismatches if the wrong mass type is searched.
%[text]
%[text] For [M+H]+ ions (protonated molecules in positive mode ESI-MS):
%[text] Observed m/z = ExactMolWt + 1.00728 (proton mass)
%[text] (1.00728 Da = proton mass; electron mass ~0.00055 Da is negligible)
logSection("A08", "Section 1: Building a Precise Mass Reference Table", "Analytics L3");
DATA_FILE = "data/list/fda_drugs.csv";

rawTbl = readtable(DATA_FILE, TextType="string");
nRaw   = height(rawTbl);
logInfo("Read %d rows from %s", nRaw, DATA_FILE);

%[text] Calculate the exact mass and molecular formula for each compound.
exactMass = nan(nRaw, 1);
formula   = strings(nRaw, 1);
valid     = false(1, nRaw);

logInfo("Calculating exact mass (may take 1-2 minutes)...");
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
logInfo("Reference table: %d compounds", nRef);

%[text] Preview
disp(refTbl(1:5, ["Name","Formula","ExactMass"]));

%[text] **💡 Observation Point 1 — Check the Exact Mass of Aspirin and Another Drug**
%[text] Let's check the exact monoisotopic mass of Aspirin.
%[text] Execute: mol = `emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O")`; 
%[text] d   = `emk.descriptor.calculate(mol, "ExactMolWt")`; 
%[text] fprintf("Aspirin ExactMolWt = %.4f Da\n", d.ExactMolWt)
%[text] Check the expected [M+H]+ m/z in the positive mode ESI spectrum.
%[text] (Add 1.00728 Da for the proton.)
%[text] The base peak in the mass spectrum of Ibuprofen appears at m/z 207.14.
%[text] Verify if this matches the [M+H]+ of Ibuprofen (C13H18O2).
%[text] Verify: `emk.descriptor.calculate(emk.mol.fromSmiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1"), "ExactMolWt")`
% ... (Try writing code here)

%%
%[text] ## Section 2: Simulation of 5 Unknown ESI-MS Spectra
%[text]
%[text] The reference table is complete. Next, we simulate realistic ESI-MS spectra for 5 unknown drugs.
%[text] This reproduces the spectra measured by actual instruments, preparing for identification exercises.
%[text]
%[text] ### Concept: ESI-MS Spectra and Isotope Clusters
%[text]
%[text] - Isotope clusters of M, M+1, M+2 are formed around **[M+H]+** ions.
%[text] - Includes baseline noise, H2O, CO2, and other fragment ions.
%[text] - Simulations are reproducible with a fixed seed.
logSection("A08", "Section 2: Simulation of 5 Unknown ESI-MS Spectra", "Analytics L3");
PROTON_MASS  = 1.007276;   % Da
MZ_SIGMA_LR  = 0.3;        % Da -- Simulated low-resolution peak width (FWHM ~0.7 Da)
NOISE_LEVEL  = 0.03;       % Proportion to base peak
N_UNKNOWNS   = 5;

%[text] Select N_UNKNOWNS drugs evenly distributed across the reference table.
%[text] (Fixed selection -- does not depend on rng)
unknownIdx = round(linspace(ceil(nRef * 0.05), floor(nRef * 0.95), N_UNKNOWNS));

unknownNames    = refTbl.Name(unknownIdx);
unknownFormulas = refTbl.Formula(unknownIdx);
unknownMass     = refTbl.ExactMass(unknownIdx);

logInfo("Selected unknown substances:");
for k = 1:N_UNKNOWNS
    logInfo("  Unknown substance %d: %s  (%.4f Da, molecular formula %s)", ...
        k, unknownNames(k), unknownMass(k), unknownFormulas(k));
end

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

    % Keep only peaks within valid m/z range (above 50 Da)
    keep  = pkMz > 50;
    spectra{k} = simulateSpectrum(pkMz(keep), pkInt(keep), ...
                     [max(50, mH-80), mH+10], MZ_SIGMA_LR, NOISE_LEVEL);
end

%[text] Display the spectrum of Unknown Substance 1.
figure("Name","A08: Sample Spectrum (Unknown Substance 1)");
plot(spectra{1}.mz, spectra{1}.intensity, "b-", LineWidth=0.8);
xlabel("m/z"); ylabel("Relative Intensity");
title(sprintf("Simulated ESI-MS Spectrum -- Unknown Substance 1 (Labels Hidden)"));
grid on;

%[text] **💡 Observation Point 2 — Let's observe the spectrum of another unknown substance**
%[text] Let's check the spectrum of Unknown Substance 3. Change the figure index to 3 and
%[text] re-run to see if you can visually identify the molecular ion cluster.
%[text] (Find the highest peak cluster on the high m/z side of the spectrum.)
%[text] The "molecular ion region" of an ESI spectrum is usually the highest m/z
%[text] cluster (for singly charged ions). Check the x-axis range.
%[text] What is the approximate m/z of the [M+H]+ for Unknown Substance 1?
% ... (Try writing code here)
%%
%[text] ## Section 3: Peak Detection Using Signal Processing Toolbox
%[text]
%[text] The spectrum is ready. Next, we will use the Signal Processing Toolbox to remove noise and automatically detect peaks.
%[text] This is preprocessing to accurately read the m/z of [M+H]+.
%[text]
%[text] ### Concept: findpeaks and smoothdata
%[text]
%[text] - Use `smoothdata(..., "gaussian")` to reduce high-frequency noise and improve detection accuracy without distorting peak shapes.
%[text] - Key parameters for `findpeaks`: MinPeakHeight, MinPeakProminence, MinPeakDistance
%[text] - Since isotope peaks are spaced 1 Da apart, MinPeakDistance should be set to less than 1 Da.
logSection("A08", "Section 3: Peak Detection Using Signal Processing Toolbox", "Analytics L3");
MIN_PEAK_HEIGHT      = 0.05;   % 5% of normalized intensity
MIN_PEAK_PROMINENCE  = 0.04;
MIN_PEAK_DISTANCE    = 0.8;    % Da -- needs to be less than 1 to resolve isotope peaks

%[text] Demonstration with Unknown Substance 1
sp1       = spectra{1};
smoothInt = smoothdata(sp1.intensity, "gaussian", 5);

[pkHeight, pkLoc] = findpeaks(smoothInt, sp1.mz, ...
    "MinPeakHeight",     MIN_PEAK_HEIGHT, ...
    "MinPeakProminence", MIN_PEAK_PROMINENCE, ...
    "MinPeakDistance",   MIN_PEAK_DISTANCE);

%[text] Normalize the peak intensities.
pkHeight = pkHeight / max(pkHeight);

logInfo("Unknown Substance 1: Detected %d peaks", numel(pkLoc));
if numel(pkLoc) > 0
    logInfo("  Detected m/z values: %s", strjoin(string(round(pkLoc, 3)), ", "));
end

%[text] The molecular ion [M+H]+ is the strongest peak within the high m/z cluster (within 3 Da of the highest detected m/z).
%[text] Using only max(pkLoc) may select the M+2 isotope peak for Cl/S compounds,
%[text] for example, the M+2 of Furosemide is about 39%, causing a systematic 2 Da error in back-calculating the neutral mass.
highMzMask  = pkLoc >= max(pkLoc) - 3;
clusterHts  = pkHeight(highMzMask);
clusterLocs = pkLoc(highMzMask);
[~, iM]     = max(clusterHts);
mH_obs_1    = clusterLocs(iM);
logInfo("Unknown Substance 1: Extracted [M+H]+ = %.4f Da", mH_obs_1);

%[text] Visualize raw data, smoothed data, and detected peaks.
figure("Name","A08: Peak Detection (Unknown Substance 1)");
plot(sp1.mz, sp1.intensity, Color=[0.7 0.7 0.7], DisplayName="Raw Data"); hold on;
plot(sp1.mz, smoothInt,     "b-", LineWidth=1.2, DisplayName="Smoothed");
plot(pkLoc,  pkHeight,      "rv", MarkerSize=8, MarkerFaceColor="r", ...
     DisplayName="Detected Peaks");
xlabel("m/z"); ylabel("Relative Intensity");
title("Unknown Substance 1 -- Raw Data, Smoothed, Detected Peaks");
legend; grid on;

%[text] Use findpeaks to extract approximate [M+H]+ for all unknown substances.
%[text] A grid step of ~0.018 Da limits accuracy to ~50-100 ppm, sufficient for low-resolution searches at 0.5 Da, but insufficient for high-resolution searches at 5 ppm.
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

%[text] Simulate calibrated high-resolution mass readings (Sections 4-6).
%[text] Actual Orbitrap/Q-TOF reports mass with ~1 ppm rms due to internal lock mass calibration.
%[text] This is fundamentally different from the grid-limited findpeaks centroid above,
%[text] where the mass analyzer (not the spectral image) determines accuracy.
%[text] Modeling: True [M+H]+ with added Gaussian(0, 1 ppm) noise.
rng(99, "twister");
mH_hires = (unknownMass + PROTON_MASS)' + ...
           (unknownMass + PROTON_MASS)' .* randn(1, N_UNKNOWNS) * 1e-6;

%[text] **💡 Observation Point 3 — Let's adjust the peak detection parameters**
%[text] Set MinPeakHeight to 0.01 and observe the changes.
%[text] Observe how many false noise peaks appear, and also check the changes when set to 0.20.
%[text] Try setting MinPeakDistance to 0.1. This will allow findpeaks to resolve isotope peaks individually.
%[text] Check if you can identify the triplet of M, M+1, M+2.
%[text] (Verify if the three peaks [M+H]+, [M+H]+1, [M+H]+2 are consecutive.)
% ... (Try writing code here)


%%
%[text] ## Section 4: Accurate Mass Database Search
%[text]
%[text] We have obtained the m/z values of the peaks. Next, we will narrow down the candidates by matching with a reference table.
%[text] Let's compare the difference in the number of candidates between low and high resolution.
%[text]
%[text] ### Concept: Mass Accuracy in ppm (parts per million)
%[text]
%[text] Mass spectrometers vary greatly in mass accuracy.
%[text]
%[text] Low resolution (single quadrupole, ion trap):
%[text] At m/z 300, the accuracy is about 0.2 to 0.5 Da, which corresponds to approximately 700 to 1700 ppm.
%[text] Many drugs share the same integer nominal mass.
%[text] Example: Nominal mass 300 matches over 40 FDA-approved drugs.
%[text]
%[text] High resolution (Orbitrap, Q-TOF, FT-ICR):
%[text] Accuracy is less than 5 ppm (less than 0.0015 Da at m/z 300).
%[text] Candidates can be narrowed down from dozens to 1 or 2.
%[text]
%[text] ppm error: (observed mass - theoretical mass) / theoretical mass * 1e6
%[text]
%[text] We use two mass variables representing two different instruments for the search.
%[text]
%[text] mH_obs(k) -- findpeaks centroid from simulated spectrum (about 50 to 100 ppm)
%[text] Used for Da window search in low resolution.
%[text] mH_hires(k) -- Simulated Orbitrap reading: true [M+H]+ with 1 ppm noise added
%[text] Used for ppm window search in high resolution.
%[text]
%[text] This separation is physically accurate. Unit resolution detectors report spectral peak positions with about 0.5 Da accuracy,
%[text] while internal mass calibration of high-resolution instruments provides accurate mass within about 1 to 5 ppm.
%[text]
%[text] LOW_RES_TOL = 0.50 Da (simulate unit resolution instrument)
%[text] HIGH_RES_PPM = 5.0 ppm (simulate Orbitrap / Q-TOF)
logSection("A08", "Section 4: Accurate Mass Database Search", "Analytics L3");
LOW_RES_TOL  = 0.50;   % Da
HIGH_RES_PPM = 5.0;    % ppm

candidateCounts = zeros(N_UNKNOWNS, 2);  % [low resolution, high resolution]

for k = 1:N_UNKNOWNS
    % Low resolution search (Da window) -- using findpeaks extraction (~50 to 100 ppm)
    lowResCand = [];
    if ~isnan(mH_obs(k))
        mNeutral_lr = mH_obs(k) - PROTON_MASS;
        daDiff      = abs(refTbl.ExactMass - mNeutral_lr);
        lowResCand  = find(daDiff <= LOW_RES_TOL);
    end

    % High resolution search (ppm window) -- using calibrated instrument reading (~1 ppm)
    mNeutral_hr = mH_hires(k) - PROTON_MASS;
    ppmDiff     = abs(refTbl.ExactMass - mNeutral_hr) ./ refTbl.ExactMass * 1e6;
    highResCand = find(ppmDiff <= HIGH_RES_PPM);

    candidateCounts(k, :) = [numel(lowResCand), numel(highResCand)];

    logInfo("Unknown substance %d (high resolution [M+H]+ = %.4f, neutral = %.4f Da):", ...
        k, mH_hires(k), mNeutral_hr);
    logInfo("  Low resolution (+-%.2f Da): %d candidates", LOW_RES_TOL, numel(lowResCand));
    logInfo("  High resolution (%.1f ppm): %d candidates", HIGH_RES_PPM, numel(highResCand));
end

%[text] Display a bar graph of the number of candidates.
figure("Name","A08: Number of Candidates by Resolution");
bar(1:N_UNKNOWNS, candidateCounts);
legend(sprintf("Low resolution (+-%.2f Da)", LOW_RES_TOL), ...
       sprintf("High resolution (%.0f ppm)", HIGH_RES_PPM));
xlabel("Unknown Sample"); ylabel("Number of Database Candidates");
title("Candidate Reduction: Low vs High Resolution MS");
grid on;
maxCand = ceil(max(candidateCounts(:)));
yticks(0:maxCand + 1);
ylim([0, maxCand + 1.5]);

%[text] **💡 Observation Point 4 — Change the ppm tolerance and observe the change in the number of candidates**
%[text] Let's check how many candidates there are if HIGH_RES_PPM is changed to 50 for unknown substance 1.
%[text] Next, try 1 ppm and see if the true compound remains.
%[text] (mH_hires has about 1 ppm rms noise. With a 1 ppm threshold, 1 to 2 unknown substances may fall just outside the window due to noise.)
%[text] Kind & Fiehn (2006) showed that even at 1 ppm, there are molecular formulas that cannot be distinguished by mass alone.
%[text] Let's see what additional information can be used to further narrow down the candidates.
%[text] (Hint: Try reading ahead to Section 5.)
% ... (Try writing code here)
%%
%[text] ## Section 5: Isotope Pattern Scoring
%[text]
%[text] The candidates have been narrowed down by accurate mass. Next, we perform final identification using isotope patterns.
%[text] Candidates are ranked using cosine similarity scores.
%[text]
%[text] ### Concept: Isotope Clusters as Structural Fingerprints
%[text]
%[text] All elements have characteristic distributions of isotopes in nature.
%[text] The relative intensities of M, M+1, and M+2 peaks indicate the elemental composition of the molecule.
%[text]
%[text] Major contributions (for singly charged organic drug molecules):
%[text]
%[text] M+1 relative abundance (% of M):
%[text] 13C: 1.103% per carbon atom (dominant contribution)
%[text] 15N: 0.366% per nitrogen atom
%[text] 2H:  0.015% per hydrogen atom
%[text] 17O: 0.038% per oxygen atom
%[text] --> M+1% ~ 1.10*nC + 0.37*nN + 0.015*nH + 0.04*nO
%[text]
%[text] M+2 relative abundance (% of M):
%[text] Contribution of 13C^2 (binomial): (1.10*nC)^2 / 200
%[text] 18O: 0.205% per oxygen atom
%[text] 34S: 4.25% per sulfur atom (large -- S is an identifier!)
%[text] 37Cl: ~32.0% per chlorine atom (very large -- Cl is highly distinctive)
%[text] Derivation: M+2/M = p(37Cl)/p(35Cl) = 24.23%/75.77% = 31.98%
%[text] Gross (2011) Table 3.2 approximates to 32.7% (difference < 2%)
%[text] Note: 24.23% is the natural abundance of 37Cl (IUPAC); 32.0% is
%[text] the relative intensity of M+2 (normalized to the 35Cl single isotope peak)
%[text] --> M+2% ~ (1.10*nC)^2/200 + 0.21*nO + 4.25*nS + 32.7*nCl
%[text]
%[text] Molecules containing Cl or S show unusually large M+2 peaks.
%[text] Drugs with one Cl atom show an M:M+2 ratio of about 3:1,
%[text] unmistakable in the spectrum.
%[text]
%[text] Scoring strategy: Cosine similarity of observed isotope vector [I_M, I_{M+1}, I_{M+2}] (normalized) with
%[text] theoretical vectors for each candidate.
%[text]
%[text] For each unknown: Extract observed isotope cluster from spectrum,
%[text] score all high-resolution candidates by isotope cosine similarity.
logSection("A08", "Section 5: Isotope Pattern Scoring", "Analytics L3");
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

%[text] Display isotope score ranking for Unknown 1
logInfo("Unknown 1 -- High-resolution candidates ranked by isotope score:");
if ~isempty(isoScoreTbls{1})
    disp(isoScoreTbls{1}(1:min(5,height(isoScoreTbls{1})), :));
end

%[text] **💡 Observation Point 5 — Check the top candidate by isotope score**
%[text] Look at the isotope score table for Unknown 1 and verify if the top candidate
%[text] matches the true identification.
%[text] (True identification is: unknownNames(1) -- displayed in Section 2.)
%[text] Search the FDA database for compounds containing at least one Cl atom,
%[text] and predict their M+2 relative intensity.
%[text] Next check: Verify if it appears in refTbl.
%[text] (Hint: Search for "Cl" in the Formula column:
%[text] refTbl(contains(refTbl.Formula, "Cl"), :))
% ... (Try writing code here)
%%
%[text] ## Section 6: Full Identification Workflow
%[text]
%[text] Isotope scoring is complete. Finally, integrate all three stages to run the full identification workflow.
%[text] Verify whether all five unknown compounds can be correctly identified.
%[text]
%[text] ### Concept: Compound Identification by Combined Mass + Isotope Evidence
%[text]
%[text] Apply the full workflow in sequence:
%[text] Step 1: Extract [M+H]+ from the raw spectrum (findpeaks)
%[text] Step 2: Back-calculate the neutral exact mass (subtract the proton)
%[text] Step 3: High-resolution mass search (5 ppm window)
%[text] Step 4: Rank remaining candidates by isotope score (cosine similarity)
%[text] Step 5: Select the top-ranked candidate as the identification
%[text]
%[text] The final identification is compared against the known true identifications
%[text] (stored in unknownNames) to compute identification accuracy.
%[text]
%[text] Connection to S05 (Unknown Compound Identification):
%[text] S05 used fingerprint similarity (Tanimoto) to identify unknown compounds
%[text] from a structural information database.
%[text] A08 achieves the same goal using physical signals (mass spectra),
%[text] without requiring prior structural knowledge of the unknowns.
%[text] In practice, both approaches are complementary: MS provides exact mass and
%[text] isotope formula, while FP similarity confirms matches when SMILES from a
%[text] reference library are available.
logSection("A08", "Section 6: Full Identification Workflow", "Analytics L3");
logInfo("Running full identification workflow...");
nCorrect = 0;
resultRows = cell(N_UNKNOWNS, 1);

for k = 1:N_UNKNOWNS
    if isempty(isoScoreTbls{k})
        logWarn("Unknown %d: insufficient data for identification", k);
        continue;
    end

    predicted = isoScoreTbls{k}.Name(1);   % top-ranked candidate
    trueID    = unknownNames(k);
    correct   = strcmpi(predicted, trueID);
    if correct, nCorrect = nCorrect + 1; end

    resultRows{k} = {k, trueID, predicted, ...
                     isoScoreTbls{k}.ppmError(1), ...
                     isoScoreTbls{k}.IsoScore(1), correct};

    resultLabel = " WRONG"; if correct, resultLabel = " CORRECT"; end
    logInfo("Unknown %d:  true=%s  predicted=%s  [%s]", ...
        k, trueID, predicted, resultLabel);
end

%[text] Summary table
validRows = find(~cellfun(@isempty, resultRows));
if ~isempty(validRows)
    resultTbl = cell2table(vertcat(resultRows{validRows}), ...
        VariableNames=["Unknown","TrueID","Predicted","ppmError","IsoScore","Correct"]);
    disp(resultTbl);
    logInfo("Identification accuracy: %d / %d = %.0f%%", ...
        nCorrect, numel(validRows), nCorrect/numel(validRows)*100);
end

%[text] **💡 Observation Point 6 — Explore the Effect of Mass Error and the MS vs FP Approaches**
%[text] The workflow above uses ExactMolWt from RDKit (monoisotopic mass).
%[text] In real experiments, measured m/z values carry instrument uncertainty.
%[text] Try adding a small random mass error to mH_obs before database search:
%[text] mH_obs_noisy = mH_obs + randn(1, N_UNKNOWNS) * 0.002;  % 2 mDa noise
%[text] Then re-run Sections 4–6 and observe how identification accuracy changes.
%[text] The isotope cosine similarity score ranges from 0 to 1.
%[text] What threshold would you apply to reject an identification as
%[text] "not in database" rather than returning a wrong answer?
%[text] (Hint: Run the workflow on a compound absent from the reference table
%[text] and observe the score distribution of false candidates.)
%[text] S05 used Tanimoto fingerprint similarity to identify unknowns from a structural database.
%[text] For compounds truly absent from the database, Tanimoto finds the most structurally similar analogue.
%[text] Discussion: When is the MS approach (A08) preferable, and when is the FP approach (S05) preferable?
%[text] What happens when you combine both?
% ... (Try writing code here)
% ... (Try writing code here)

%[text] **Summary**
%[text]
%[text] - Exact monoisotopic mass (ExactMolWt) is the foundation of MS database search
%[text] - High-resolution MS (<5 ppm) narrows candidates from dozens to 1–2
%[text] - Molecules containing Cl or S show abnormally large M+2 peaks, serving as an elemental "fingerprint"
%[text] - The two-stage filter of mass + isotope cosine similarity improves identification accuracy
%[text]
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
%     After uniform noise addition (noiseLevel = 0.03), max intensity can reach up to 1.03.
%     This is physically realistic (random noise atop peak apex) and intentional.
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
%     [A-Z][a-z]?  Matches element symbols: 1 uppercase + optional lowercase
%                  --> Correctly captures 2-letter symbols: Cl, Br, etc.
%                  --> "Cl" matches as token {"Cl",""}, not {"C",""} + {"l",""}
%     (\d*)        Matches count digits; empty string --> treated as 1
% 
%   Verified examples (Hill notation returned by RDKit rdMolDescriptors):
%     "CH4"       --> C:1, H:4
%     "C9H8O4"    --> C:9, H:8, O:4 (Aspirin)
%     "C6H5Cl"    --> C:6, H:5, Cl:1 (Chlorobenzene -- 2-letter element OK)
%     "CCl4"      --> C:1, Cl:4 (Carbon tetrachloride)
%     "H2O"       --> H:2, O:1 (no digit for last element -- guard OK)
% 
%   Edge cases:
%     Empty string ""    --> All counts 0 (tokens = {}; loop not entered)
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
% ISOPATTERN  Calculates relative isotopic intensities of M, M+1, M+2.
% 
%   rel = isoPattern(fc)
%   fc  -- Element count struct from parseFormula()
%   rel -- [1, M+1/M, M+2/M] (M normalized to 1.0)
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
%           M+2 coefficient is ratio of 37Cl/35Cl, not abundance.
% 
%   Cross-validation reference values (validated in Claesen et al. 2012):
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
%[metadata:view]%---

%   data: {"layout":"inline","rightPanelPercent":40}
%---
