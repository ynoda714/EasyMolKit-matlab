%[text] # A08 Answers: Mass Spectrometry × Cheminformatics
%[text] Reference answers for the "Try It Yourself" exercises in a08_mass_spec_cheminformatics.m.
%[text] First, execute a08_mass_spec_cheminformatics.m (at least sections 0-5)
%[text] to construct the necessary workspace variables. Then verify with this file.
%[text]
%[text] Prerequisites reconstructed here: exact mass reference table, 5 simulated spectra,
%[text] peak detection results, isotope scoring table.
addpath(genpath("src"));
emk.setup.initPython();

mol_warmup = emk.mol.fromSmiles("C"); clear mol_warmup;
logInfo("A08 Answers: Setup complete");
%%
%[text] ## Rebuild prerequisites (mirrors Sections 0-5 of a08)

DATA_FILE    = "data/list/fda_drugs.csv";
PROTON_MASS  = 1.007276;
MZ_SIGMA_LR  = 0.3;
NOISE_LEVEL  = 0.03;
N_UNKNOWNS   = 5;
LOW_RES_TOL  = 0.50;
HIGH_RES_PPM = 5.0;
MIN_PEAK_HEIGHT     = 0.05;
MIN_PEAK_PROMINENCE = 0.04;
MIN_PEAK_DISTANCE   = 0.8;

rawTbl = readtable(DATA_FILE, TextType="string");
nRaw   = height(rawTbl);
exactMass = nan(nRaw, 1);
formula   = strings(nRaw, 1);
valid     = false(1, nRaw);

logInfo("Building exact mass reference table (%d compounds)...", nRaw);
for k = 1:nRaw
    try
        mol          = emk.mol.fromSmiles(rawTbl.SMILES(k));
        d            = emk.descriptor.calculate(mol, ["ExactMolWt","MolFormula"]);
        exactMass(k) = d.ExactMolWt;
        formula(k)   = d.MolFormula;
        valid(k)     = true;
    catch; end
end
validIdx = find(valid);
refTbl = table( ...
    rawTbl.Name(validIdx), rawTbl.SMILES(validIdx), ...
    formula(validIdx), exactMass(validIdx), ...
    VariableNames=["Name","SMILES","Formula","ExactMass"]);
nRef = height(refTbl);

unknownIdx      = round(linspace(ceil(nRef*0.05), floor(nRef*0.95), N_UNKNOWNS));
unknownNames    = refTbl.Name(unknownIdx);
unknownFormulas = refTbl.Formula(unknownIdx);
unknownMass     = refTbl.ExactMass(unknownIdx);

rng(2026);
spectra = cell(1, N_UNKNOWNS);
for k = 1:N_UNKNOWNS
    mH   = unknownMass(k) + PROTON_MASS;
    fc   = parseFormula(unknownFormulas(k));
    iso  = isoPattern(fc);
    pkMz = [mH, mH+1, mH+2, mH-18, mH-44];
    pkInt= [iso(1), iso(2), iso(3), iso(1)*0.35, iso(1)*0.15];
    keep = pkMz > 50;
    spectra{k} = simulateSpectrum(pkMz(keep), pkInt(keep), ...
                     [max(50, mH-80), mH+10], MZ_SIGMA_LR, NOISE_LEVEL);
end

mH_obs = nan(1, N_UNKNOWNS);
for k = 1:N_UNKNOWNS
    sp = spectra{k};
    sm = smoothdata(sp.intensity, "gaussian", 5);
    [pks, locs] = findpeaks(sm, sp.mz, ...
        "MinPeakHeight", MIN_PEAK_HEIGHT, ...
        "MinPeakProminence", MIN_PEAK_PROMINENCE, ...
        "MinPeakDistance", MIN_PEAK_DISTANCE);
    if ~isempty(locs)
        highMzMask = locs >= max(locs) - 3;
        clusterLocs = locs(highMzMask);
        [~, iM]    = max(pks(highMzMask));
        mH_obs(k)  = clusterLocs(iM);
    end
end

rng(99, "twister");
mH_hires = (unknownMass + PROTON_MASS)' + ...
           (unknownMass + PROTON_MASS)' .* randn(1, N_UNKNOWNS) * 1e-6;

%[text] High-res candidates + isotope scoring (needed for TRY IT 5-6)
isoScoreTbls = cell(1, N_UNKNOWNS);
for k = 1:N_UNKNOWNS
    mNeutral   = mH_hires(k) - PROTON_MASS;
    ppmDiff    = abs(refTbl.ExactMass - mNeutral) ./ refTbl.ExactMass * 1e6;
    candIdx    = find(ppmDiff <= HIGH_RES_PPM);
    if isempty(candIdx), continue; end
    sp    = spectra{k};
    smInt = smoothdata(sp.intensity, "gaussian", 5);
    smInt = smInt / max(smInt);
    obsIso = zeros(1,3);
    for ishift = 0:2
        targetMz   = mH_hires(k) + ishift;
        [~, bestPt] = min(abs(sp.mz - targetMz));
        obsIso(ishift+1) = max(smInt(max(1,bestPt-3):min(end,bestPt+3)));
    end
    if obsIso(1) < 1e-6, obsIso(1) = 1; end
    obsIso = obsIso / obsIso(1);
    nCand  = numel(candIdx);
    scores = nan(nCand,1);
    for ci = 1:nCand
        fc      = parseFormula(refTbl.Formula(candIdx(ci)));
        theoIso = isoPattern(fc);
        scores(ci) = dot(obsIso, theoIso) / (norm(obsIso) * max(norm(theoIso),1e-12));
    end
    [scoresSorted, sortOrd] = sort(scores, "descend");
    isoScoreTbls{k} = table( ...
        refTbl.Name(candIdx(sortOrd)), refTbl.Formula(candIdx(sortOrd)), ...
        refTbl.ExactMass(candIdx(sortOrd)), ppmDiff(candIdx(sortOrd)), ...
        scoresSorted, ...
        VariableNames=["Name","Formula","ExactMass","ppmError","IsoScore"]);
end

logInfo("Prerequisites rebuilt (%d reference compounds, %d unknowns).", nRef, N_UNKNOWNS);
%%
%[text] ## Let's Try 1: Aspirin Exact Mass, [M+H]+ m/z; Verify Ibuprofen Base Peak

mol_asp  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
d_asp    = emk.descriptor.calculate(mol_asp, "ExactMolWt");
mH_asp   = d_asp.ExactMolWt + PROTON_MASS;

logInfo("Let's Try 1 -- Aspirin:");
logInfo("  Exact Molecular Weight = %.4f Da  (C9H8O4)", d_asp.ExactMolWt);
logInfo("  [M+H]+ m/z  = %.4f Da  (adding proton 1.00728 Da)", mH_asp);
%[text] Answer: Aspirin ExactMolWt = 180.0423 Da; [M+H]+ = 181.0495 Da.
%[text]    The AVERAGE MolWt (180.16) would produce a search error of ~0.12 Da --
%[text]    large enough to miss a 5 ppm search window at m/z 181.
mol_ibu  = emk.mol.fromSmiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1");
d_ibu    = emk.descriptor.calculate(mol_ibu, "ExactMolWt");
mH_ibu_calc = d_ibu.ExactMolWt + PROTON_MASS;
mH_ibu_obs  = 207.14;   % stated base peak
ppm_ibu = (mH_ibu_obs - mH_ibu_calc) / mH_ibu_calc * 1e6;

logInfo("Let's Try 1 -- Ibuprofen:");
logInfo("  ExactMolWt  = %.4f Da  (C13H18O2)", d_ibu.ExactMolWt);
logInfo("  Theoretical [M+H]+ = %.4f Da", mH_ibu_calc);
logInfo("  Observed Base Peak = %.2f Da", mH_ibu_obs);
logInfo("  ppm error           = %+.1f ppm", ppm_ibu);
%[text] Answer: Ibuprofen ExactMolWt ~ 206.1307 Da => [M+H]+ theoretical = 207.1380 Da.
%[text]    Observed 207.14 => error ~ +9.9 ppm (207.14 has only 2 decimal places;
%[text]    rounding alone accounts for most of the error).  This illustrates why
%[text]    observed m/z values must be reported to >= 4 decimal places for sub-5 ppm
%[text]    accuracy.  Note: ibuprofen commonly shows a base peak at m/z 161 (loss of
%[text]    CO2+H2O) in typical EI spectra; ESI-MS gives predominant [M+H]+.
%%
%[text] ## Let's Try 2: Investigate Unknown Substance 3 Spectrum; Identify Molecular Ion Cluster

figure("Name","A08 Answers: Unknown 3 Spectrum");
plot(spectra{3}.mz, spectra{3}.intensity, "b-", LineWidth=0.8);
xlabel("m/z"); ylabel("Relative Intensity");
title(sprintf("Unknown Substance 3 Spectrum (True Compound: %s)", unknownNames(3)));
grid on;

mH_true3 = unknownMass(3) + PROTON_MASS;
logInfo("Let's Try 2 -- Unknown Substance 3:");
logInfo("  True Compound: %s  (Formula: %s)", unknownNames(3), unknownFormulas(3));
logInfo("  True [M+H]+: %.4f Da", mH_true3);
logInfo("  The molecular ion cluster appears at the high m/z end of the spectrum.");
%[text] Answer: The molecular ion cluster (M, M+1, M+2) is the group of peaks at the
%[text]    rightmost (highest m/z) position on the spectrum.  For singly charged ions
%[text]    each peak in the cluster is separated by ~1 Da.  Peaks at M-18 and M-44
%[text]    (neutral losses) appear as smaller peaks ~18 and ~44 Da below the cluster.
%%
%[text] ## Let's Try 3: MinPeakHeight Sensitivity; MinPeakDistance for Isotope Resolution

sp1       = spectra{1};
smoothInt = smoothdata(sp1.intensity, "gaussian", 5);

logInfo("Let's Try 3 -- MinPeakHeight Sensitivity:");
for mph = [0.01, 0.05, 0.10, 0.20]
    [pks, ~] = findpeaks(smoothInt, sp1.mz, ...
        "MinPeakHeight",     mph, ...
        "MinPeakProminence", 0.01, ...
        "MinPeakDistance",   MIN_PEAK_DISTANCE);
    logInfo("  MinPeakHeight = %.2f  -->  %d peaks detected", mph, numel(pks));
end
%[text] Answer: MinPeakHeight = 0.01 detects many noise spikes (10-20+ peaks in a typical
%[text]    noisy spectrum).  MinPeakHeight = 0.20 may miss the M+2 isotope peak for
%[text]    molecules with few heavy atoms (M+2 often ~1-5% relative intensity).
%[text]    The default 0.05 (5%) is a practical compromise.
logInfo("Let's Try 3 -- MinPeakDistance for Isotope Resolution = 0.1:");
[pkH_iso, pkL_iso] = findpeaks(smoothInt, sp1.mz, ...
    "MinPeakHeight",     MIN_PEAK_HEIGHT, ...
    "MinPeakProminence", MIN_PEAK_PROMINENCE, ...
    "MinPeakDistance",   0.1);
mH_true1 = unknownMass(1) + PROTON_MASS;
%[text] Identify the M, M+1, M+2 triplet: look for consecutive peaks near mH_true1
isoMask = pkL_iso >= mH_true1 - 0.3 & pkL_iso <= mH_true1 + 2.5;
if sum(isoMask) >= 3
    isoLocs = pkL_iso(isoMask);
    isoHts  = pkH_iso(isoMask) / max(pkH_iso(isoMask));
    logInfo("  Isotope triplet for Unknown 1 (near %.4f Da):", mH_true1);
    for ii = 1:min(3, sum(isoMask))
        logInfo("    Peak %d: m/z = %.4f  Relative Intensity = %.3f", ii, isoLocs(ii), isoHts(ii));
    end
else
    logInfo("  Triplet not clearly resolved (smoothing already applied -- peaks overlap at sigma=%.1f Da; narrower Gaussian or higher-res simulation needed).", MZ_SIGMA_LR);
end
%[text] Answer: Setting MinPeakDistance = 0.1 Da resolves peaks 1 Da apart (the isotope
%[text]    cluster spacing).  Three consecutive peaks at [M+H]+, [M+H]+1, [M+H]+2
%[text]    become visible.  Relative intensities match the isoPattern formula.
%[text]    Warning: noise at this low distance threshold also increases false positives
%[text] -- use with smoothdata() pre-processing.
%%
%[text] ## Let's Try 4: HIGH_RES_PPM Sensitivity; Additional Disambiguation Information

HIGH_RES_TESTS = [50, 10, 5, 2, 1];
logInfo("Let's Try 4 -- HIGH_RES_PPM vs. Number of Candidates (Unknown Substance 1):");
allSingle = true;
for ppm = HIGH_RES_TESTS
    mNeutral = mH_hires(1) - PROTON_MASS;
    ppmDiff  = abs(refTbl.ExactMass - mNeutral) ./ refTbl.ExactMass * 1e6;
    nCand    = sum(ppmDiff <= ppm);
    trueRank = find(strcmpi(refTbl.Name(ppmDiff <= ppm), unknownNames(1)));
    if isempty(trueRank), trueStr = "NOT FOUND"; else, trueStr = "found"; end
    candLabel = "candidates"; if nCand == 1, candLabel = "candidate "; end
    logInfo("  %4.0f ppm -> %2d %s  (true compound: %s)", ppm, nCand, candLabel, trueStr);
    if nCand > 1, allSingle = false; end
end
if allSingle
    logInfo("  (note: 1 candidate at all ppm -- %d-compound db is too small to produce ambiguity for m/z %.1f; real FDA db has ~1500 approved small molecules)", nRef, mH_hires(1));
end
%[text] Answer: At 50 ppm many candidates survive (5-15 for typical drug masses).
%[text]    At 5 ppm only 1-3 remain in a 200-compound database.
%[text]    At 1 ppm the true compound may fall outside the window because mH_hires
%[text]    has ~1 ppm rms noise (some unknowns will be > 1 sigma from the true value).
%[text]    Kind & Fiehn (2006) showed that even at 1 ppm, some molecular formulas
%[text]    share the same nominal mass; additional disambiguation requires:
%[text]    (a) isotope pattern scoring (Section 5), (b) retention time from LC,
%[text]    (c) MS/MS fragmentation spectra, or (d) comparison to a spectral library.
%%
%[text] ## Let's Try 5: Top-ranked Isotope Candidates; Cl-containing Compounds M+2

logInfo("Let's Try 5 -- Isotope Score Verification:");
iso5_ids    = zeros(N_UNKNOWNS, 1);
iso5_names  = strings(N_UNKNOWNS, 1);
iso5_scores = nan(N_UNKNOWNS, 1);
iso5_match  = strings(N_UNKNOWNS, 1);
for k = 1:N_UNKNOWNS
    if ~isempty(isoScoreTbls{k})
        top = isoScoreTbls{k}.Name(1);
        topScore = isoScoreTbls{k}.IsoScore(1);
        match = strcmpi(top, unknownNames(k));
        yesNo = ["no","yes"];
        iso5_ids(k)    = k;
        iso5_names(k)  = top;
        iso5_scores(k) = topScore;
        iso5_match(k)  = yesNo(match+1);
    end
end
validRows = iso5_ids > 0;
iso5Tbl = table(iso5_ids(validRows), iso5_names(validRows), iso5_scores(validRows), iso5_match(validRows), ...
    "VariableNames", ["Unknown","TopCandidate","IsoScore","Match"]);
logInfo("Isotope Score Verification Results:");
disp(iso5Tbl);
%[text] Answer: The top-ranked candidate should match the true identity for most unknowns
%[text]    because the synthetic spectra were generated from the exact same isotope
%[text]    formula as used by the scoring function (self-consistent test).
%[text]    In a real experiment the isotope pattern would differ slightly due to
%[text]    instrument calibration, detector noise, and incomplete resolution of
%[text]    overlapping peaks -- reducing the cosine score below 1.0.
%[text]
%[text] Find a Cl-containing compound and predict M+2
clMask = contains(refTbl.Formula, "Cl");
logInfo("  Cl-containing compounds in FDA reference table: %d", sum(clMask));
if any(clMask)
    clTbl = refTbl(clMask, :);
    sampleName = clTbl.Name(1);
    sampleFormula = clTbl.Formula(1);
    fc_cl   = parseFormula(sampleFormula);
    iso_cl  = isoPattern(fc_cl);
    logInfo("  Example: %s  Formula: %s", sampleName, sampleFormula);
    clOnly = fc_cl.Cl * 32.7;
    logInfo("    Theoretical M+2 relative intensity = %.1f%%  (Cl count = %d)", ...
        iso_cl(3)*100, fc_cl.Cl);
    logInfo("    Of which Cl contribution = ~%.1f%%  (Cl 3:1 rule: 32.7%% per Cl)", clOnly);
    logInfo("    Remaining %.1f%% from C/O/S natural isotope abundance", ...
        iso_cl(3)*100 - clOnly);
end
%[text] Answer: Each chlorine atom contributes ~32.7% to the M+2 relative intensity.
%[text]    A drug with 1 Cl shows M+2 ~ 33% (M:M+2 ~ 3:1), clearly visible on-screen.
%[text]    A drug with 2 Cl shows M+2 ~ 66% (M:M+2 ~ 3:2), and a 3:3:1 triplet at M+4.
%[text]    This "chlorine flag" is one of the most recognisable features in MS spectra.
%%
%[text] ## Let's Try 6 (Comprehensive): Adding Mass Error; Score Threshold; MS vs FP Comparison

logInfo("Let's Try 6 -- Impact of 2 mDa Mass Error on Identification:");
%[text] First verify noiseless baseline
nCorrect_base = 0;
for k = 1:N_UNKNOWNS
    if ~isempty(isoScoreTbls{k}) && strcmpi(isoScoreTbls{k}.Name(1), unknownNames(k))
        nCorrect_base = nCorrect_base + 1;
    end
end
logInfo("  Baseline Accuracy (Zero Noise): %d / %d", nCorrect_base, N_UNKNOWNS);
rng(7);
mH_obs_noisy = mH_hires + randn(1, N_UNKNOWNS) * 0.002;   % 2 mDa noise

nCorrect_noisy = 0;
noisy_pred  = strings(N_UNKNOWNS, 1);
noisy_match = strings(N_UNKNOWNS, 1);
noisy_valid = false(N_UNKNOWNS, 1);
for k = 1:N_UNKNOWNS
    mNeutral   = mH_obs_noisy(k) - PROTON_MASS;
    actualPpm  = abs(mH_obs_noisy(k) - mH_hires(k)) / mH_hires(k) * 1e6;
    ppmDiff    = abs(refTbl.ExactMass - mNeutral) ./ refTbl.ExactMass * 1e6;
    candIdx    = find(ppmDiff <= HIGH_RES_PPM);
    if isempty(candIdx)
        logInfo("  Unknown %d: No candidates at 2 mDa noise with %.1f ppm (Noise = %.1f ppm at m/z %.1f)", ...
            k, HIGH_RES_PPM, actualPpm, mH_hires(k));
        continue;
    end
    sp    = spectra{k};
    smInt = smoothdata(sp.intensity, "gaussian", 5);
    smInt = smInt / max(smInt);
    obsIso = zeros(1,3);
    for ishift = 0:2
        targetMz   = mH_obs_noisy(k) + ishift;
        [~, bestPt] = min(abs(sp.mz - targetMz));
        obsIso(ishift+1) = max(smInt(max(1,bestPt-3):min(end,bestPt+3)));
    end
    if obsIso(1) < 1e-6, obsIso(1) = 1; end
    obsIso = obsIso / obsIso(1);
    nCand  = numel(candIdx);
    scores = nan(nCand,1);
    for ci = 1:nCand
        fc      = parseFormula(refTbl.Formula(candIdx(ci)));
        theoIso = isoPattern(fc);
        scores(ci) = dot(obsIso, theoIso) / (norm(obsIso) * max(norm(theoIso),1e-12));
    end
    [~, topci] = max(scores);
    predicted  = refTbl.Name(candIdx(topci));
    correct    = strcmpi(predicted, unknownNames(k));
    if correct, nCorrect_noisy = nCorrect_noisy + 1; end
    yesNo = ["no","yes"];
    noisy_pred(k)  = predicted;
    noisy_match(k) = yesNo(correct+1);
    noisy_valid(k) = true;
end
noisyTbl = table((1:N_UNKNOWNS)', noisy_pred, noisy_match, ...
    "VariableNames", ["Unknown","Predicted","Match"]);
logInfo("Identification Results with 2 mDa Noise:");
disp(noisyTbl(noisy_valid, :));
logInfo("  Accuracy with 2 mDa Noise: %d / %d", nCorrect_noisy, N_UNKNOWNS);
%[text] Answer: 2 mDa noise corresponds to different ppm values depending on m/z:
%[text]    at m/z ~137 (AMPHETAMINE): 2 mDa = ~14.6 ppm => exceeds 5 ppm window
%[text]    at m/z ~195 (CAFFEINE):    2 mDa = ~10.3 ppm => exceeds 5 ppm window
%[text]    at m/z ~286-320:           2 mDa = ~6-7 ppm  => may exceed 5 ppm window
%[text]    Low-mass compounds are disproportionately affected by absolute mass error.
%[text]    Identification accuracy drops to 3/5 (60%) from 5/5 with zero noise.
%[text]    If the window is tightened to 2 ppm, even more unknowns will lose all
%[text]    candidates.  This illustrates why absolute mDa accuracy matters most for
%[text]    small molecules (<200 Da).
logInfo("Let's Try 6 -- Discussion on Isotope Score Rejection Threshold:");
logInfo("  Cosine score >= 0.9 is a practical acceptance threshold.");
logInfo("  Flag identification as uncertain if score < 0.8.");
logInfo("  False candidate score distribution is usually 0.6–0.85.");
%[text] Answer: For a compound truly absent from the database, the best cosine score
%[text]    against any database entry is usually 0.70-0.85 (the closest compound
%[text]    by formula is similar but not identical).  A threshold of >= 0.9 rejects
%[text]    most false positives.  In practice, <0.85 => report "no match found."
%[text]
%[text] MS (A08) vs FP (S05) comparison:
%[text]    MS approach (A08): uses physical mass signal -- works without knowing
%[text]      the structure upfront; resolves formulas accurately; fails when the
%[text]      compound is outside the mass reference library or when fragment ions
%[text]      are complex.  Best for dereplication in known-compound databases.
%[text]    FP approach (S05): uses structural similarity from a SMILES database;
%[text]      can find close analogues even for unknown compounds; requires a
%[text]      structural database and prior SMILES.  Best for lead-finding when
%[text]      exact MS data is unavailable.
%[text]    Combined: exact mass narrows the formula; fingerprint similarity
%[text]      confirms or ranks the structural candidates.  In metabolomics/NPS
%[text]      workflows, both stages are often applied in sequence.
%%
%[text] ## Local Functions (copied from a08_mass_spec_cheminformatics.m)

function spectrum = simulateSpectrum(peakMz, peakInt, mzRange, sigma, noiseLevel)
    mz        = linspace(mzRange(1), mzRange(2), 5000);
    intensity = zeros(1, 5000);
    for i = 1:numel(peakMz)
        intensity = intensity + peakInt(i) * ...
            exp(-0.5 * ((mz - peakMz(i)) / sigma).^2);
    end
    if max(intensity) > 0, intensity = intensity / max(intensity); end
    intensity = intensity + noiseLevel * rand(1, 5000);
    spectrum  = struct("mz", mz, "intensity", intensity);
end

function counts = parseFormula(formulaStr)
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
    m1pct = 1.103*fc.C + 0.366*fc.N + 0.015*fc.H + 0.038*fc.O;
    m2pct = (1.103*fc.C)^2 / 200 + 0.205*fc.O + 4.25*fc.S + 32.7*fc.Cl;
    rel = [1.0, m1pct/100, m2pct/100];
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
