% test_m2_smoke.m  End-to-end M2 smoke test.
%
%   Covers the full M2 pipeline on both Desktop and MATLAB Online.
%   Verifies that all core Chemoinformatics functions work together
%   as an integrated system.
%
%   Prerequisites (Desktop):
%     1. Run emk.setup.install() once to deploy Embedded Python + RDKit
%     2. addpath(genpath("src"))
%     3. run("tests/smoke/test_m2_smoke.m")
%
%   Prerequisites (MATLAB Online):
%     1. Run emk.setup.installOnline() once
%     2. addpath(genpath("src"))
%     3. run("tests/smoke/test_m2_smoke.m")
%
%   Test Cases:
%     TC01  calculate()       - aspirin: all 10 descriptors (struct)
%     TC02  batchCalculate()  - 3 molecules -> MATLAB table with MW column
%     TC03  morgan()          - ethanol: 2048-bit Python FP object
%     TC04  maccs()           - ethanol: 167-bit Python FP object
%     TC05  toArray()         - Morgan FP -> logical(1,2048); ON-bit cross-check
%     TC06  tanimoto() same   - identical FPs -> 1.0
%     TC07  tanimoto() diff   - ethanol vs aspirin -> (0, 1) range; symmetry
%     TC08  dice()            - D >= T for same pair; D=1.0 for same FP
%     TC09  writeSdf/readSdf  - round-trip: 3 molecules, SMILES preserved
%     TC10  readSmilesList()  - temp SMILES file -> 3 mol objects
%     TC11  draw2d()          - ethanol: figure handle; image 300x300 RGB uint8
%     TC12  full pipeline     - library of 5 SMILES; rank by Tanimoto to aspirin

addpath(genpath("src"));

nPass = 0;
nFail = 0;
nSkip = 0;

logInfo("=== M2 Smoke Test start ===");

% ---------------------------------------------------------------------------
% Init Python (required for all TC02+)
% ---------------------------------------------------------------------------
try
    emk.setup.initPython();
    result = emk.setup.verify();
    if ~result.rdkit
        logError("RDKit not available -- cannot run M2 smoke test");
        error("emk:smoke:test_m2_smoke:rdkitUnavailable", ...
            "RDKit not available. Run emk.setup.install() first.");
    end
    logInfo("Python %s + RDKit %s ready", string(pyenv().Version), result.version);
catch ex
    logError("Python/RDKit init failed: %s", ex.message);
    error("emk:smoke:test_m2_smoke:initFailed", ...
        "Python/RDKit init failed: %s", ex.message);
end

% Reference molecules used throughout the test
SMILES_ETHANOL  = "CCO";
SMILES_BENZENE  = "c1ccccc1";
SMILES_ASPIRIN  = "CC(=O)Oc1ccccc1C(=O)O";

molEthanol  = [];
molBenzene  = [];
molAspirin  = [];
fpMorganEth = [];
fpMorganAsp = [];

try
    molEthanol = emk.mol.fromSmiles(SMILES_ETHANOL);
    molBenzene = emk.mol.fromSmiles(SMILES_BENZENE);
    molAspirin = emk.mol.fromSmiles(SMILES_ASPIRIN);
catch ex
    logError("Reference mol creation failed: %s", ex.message);
    error("emk:smoke:test_m2_smoke:refMolFailed", ...
        "Reference mol creation failed: %s", ex.message);
end

% ---------------------------------------------------------------------------
% TC01  calculate() -- all 10 descriptors for aspirin
%   Reference (PubChem CID 2244, IUPAC 2021):
%     MolWt ~180.159 g/mol, HeavyAtomCount=13, RingCount=1, NumHDonors=1
%   algorithm_guide.md sec. 4.2
% ---------------------------------------------------------------------------
try
    desc = emk.descriptor.calculate(molAspirin);

    assert(isstruct(desc), "calculate must return struct");
    required = ["MolWt","ExactMolWt","LogP","TPSA","NumHAcceptors", ...
                "NumHDonors","NumRotatableBonds","RingCount","FractionCSP3","HeavyAtomCount"];
    for k = 1:numel(required)
        assert(isfield(desc, required(k)), ...
            sprintf("calculate result missing field: %s", required(k)));
    end

    assert(abs(desc.MolWt - 180.159) < 0.1, ...
        sprintf("aspirin MolWt=%.4f expected ~180.159", desc.MolWt));
    assert(desc.HeavyAtomCount == 13, ...
        sprintf("aspirin HeavyAtomCount=%d expected 13", desc.HeavyAtomCount));
    assert(desc.RingCount == 1, ...
        sprintf("aspirin RingCount=%d expected 1", desc.RingCount));
    assert(desc.NumHDonors == 1, ...
        sprintf("aspirin NumHDonors=%d expected 1", desc.NumHDonors));

    logInfo("TC01 PASS  calculate(aspirin): MolWt=%.3f HAC=%d rings=%d donors=%d", ...
        desc.MolWt, desc.HeavyAtomCount, desc.RingCount, desc.NumHDonors);
    nPass = nPass + 1;
catch ex
    logError("TC01 FAIL  calculate(aspirin): %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% TC02  batchCalculate() -- 3 molecules -> MATLAB table
%   Verifies: correct row count, MolWt column type, aspirin MW cross-check
% ---------------------------------------------------------------------------
try
    mols = {molEthanol, molBenzene, molAspirin};
    tbl  = emk.descriptor.batchCalculate(mols);

    assert(isa(tbl, "table"), "batchCalculate must return table");
    assert(height(tbl) == 3, ...
        sprintf("batchCalculate table must have 3 rows, got %d", height(tbl)));
    assert(ismember("MolWt", tbl.Properties.VariableNames), ...
        "batchCalculate table must have MolWt column");

    % Aspirin is row 3 -- cross-validate MolWt
    assert(abs(tbl.MolWt(3) - 180.159) < 0.1, ...
        sprintf("batchCalculate aspirin MolWt=%.4f expected ~180.159", tbl.MolWt(3)));

    % All values finite
    for c = 1:width(tbl)
        vals = tbl{:, c};
        assert(all(isfinite(vals)), ...
            sprintf("batchCalculate column %d has non-finite values", c));
    end

    logInfo("TC02 PASS  batchCalculate(3 mols): %d rows x %d cols  aspirinMW=%.3f", ...
        height(tbl), width(tbl), tbl.MolWt(3));
    nPass = nPass + 1;
catch ex
    logError("TC02 FAIL  batchCalculate: %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% TC03  morgan() -- ethanol: 2048-bit Python FP
%   Verifies: Python object returned, bit length = 2048, ON-bits > 0
% ---------------------------------------------------------------------------
try
    fpMorganEth = emk.fingerprint.morgan(molEthanol);
    fpMorganAsp = emk.fingerprint.morgan(molAspirin);

    assert(startsWith(class(fpMorganEth), "py."), ...
        sprintf("morgan must return Python object, got: %s", class(fpMorganEth)));
    assert(double(fpMorganEth.GetNumBits()) == 2048, ...
        sprintf("morgan NBits expected 2048, got %d", double(fpMorganEth.GetNumBits())));
    nOn = double(fpMorganEth.GetNumOnBits());
    assert(nOn > 0, "ethanol morgan ON-bits must be > 0");
    assert(nOn <= 2048, "ethanol morgan ON-bits must not exceed 2048");

    logInfo("TC03 PASS  morgan(ethanol): %d-bit FP, %d ON-bits", ...
        double(fpMorganEth.GetNumBits()), nOn);
    nPass = nPass + 1;
catch ex
    logError("TC03 FAIL  morgan(): %s", ex.message);
    nFail = nFail + 1;
    fpMorganEth = [];
    fpMorganAsp = [];
end

% ---------------------------------------------------------------------------
% TC04  maccs() -- ethanol: 167-bit Python FP
%   Verifies: Python object, exactly 167 bits, ON-bits in range
% ---------------------------------------------------------------------------
fpMaccsEth = [];
fpMaccsAsp = [];
try
    fpMaccsEth = emk.fingerprint.maccs(molEthanol);
    fpMaccsAsp = emk.fingerprint.maccs(molAspirin);

    assert(startsWith(class(fpMaccsEth), "py."), ...
        sprintf("maccs must return Python object, got: %s", class(fpMaccsEth)));
    assert(double(fpMaccsEth.GetNumBits()) == 167, ...
        sprintf("maccs NBits expected 167, got %d", double(fpMaccsEth.GetNumBits())));
    nOn = double(fpMaccsEth.GetNumOnBits());
    assert(nOn > 0, "ethanol maccs ON-bits must be > 0");
    assert(nOn <= 167, "ethanol maccs ON-bits must not exceed 167");

    logInfo("TC04 PASS  maccs(ethanol): 167-bit FP, %d ON-bits", nOn);
    nPass = nPass + 1;
catch ex
    logError("TC04 FAIL  maccs(): %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% TC05  toArray() -- Morgan FP -> logical(1, 2048); ON-bit cross-check
%   algorithm_guide.md sec.5.3: sum(toArray(fp)) == fp.GetNumOnBits()
% ---------------------------------------------------------------------------
if ~isempty(fpMorganEth)
    try
        bits = emk.fingerprint.toArray(fpMorganEth);

        assert(islogical(bits), ...
            sprintf("toArray must return logical, got %s", class(bits)));
        assert(isrow(bits), ...
            "toArray must return a row vector");
        assert(numel(bits) == 2048, ...
            sprintf("toArray length expected 2048, got %d", numel(bits)));

        % Cross-validation: MATLAB sum == Python GetNumOnBits
        nOnMatlab  = sum(bits);
        nOnPython  = double(fpMorganEth.GetNumOnBits());
        assert(nOnMatlab == nOnPython, ...
            sprintf("toArray ON-bit mismatch: MATLAB=%d Python=%d", nOnMatlab, nOnPython));

        logInfo("TC05 PASS  toArray(morgan eth): logical(1,%d), %d ON-bits (MATLAB==Python)", ...
            numel(bits), nOnMatlab);
        nPass = nPass + 1;
    catch ex
        logError("TC05 FAIL  toArray(): %s", ex.message);
        nFail = nFail + 1;
    end
else
    logWarn("TC05 SKIP  Morgan FP unavailable (TC03 failed)");
    nSkip = nSkip + 1;
end

% ---------------------------------------------------------------------------
% TC06  tanimoto() same FP -> 1.0
%   Verifies: identical FPs => Tanimoto = 1.0 (double scalar, isreal, isfinite)
% ---------------------------------------------------------------------------
if ~isempty(fpMorganEth)
    try
        score = emk.similarity.tanimoto(fpMorganEth, fpMorganEth);

        assert(isa(score, "double") && isscalar(score), ...
            "tanimoto must return double scalar");
        assert(isreal(score) && isfinite(score), ...
            "tanimoto must be real and finite");
        assert(abs(score - 1.0) < 1e-9, ...
            sprintf("tanimoto(same) expected 1.0, got %.10f", score));

        logInfo("TC06 PASS  tanimoto(ethanol,ethanol) = %.10f  (==1.0 confirmed)", score);
        nPass = nPass + 1;
    catch ex
        logError("TC06 FAIL  tanimoto same FP: %s", ex.message);
        nFail = nFail + 1;
    end
else
    logWarn("TC06 SKIP  Morgan FP unavailable (TC03 failed)");
    nSkip = nSkip + 1;
end

% ---------------------------------------------------------------------------
% TC07  tanimoto() different FPs: range [0,1], symmetry, < 1.0
%   Cross-validation: T = c/(a+b-c) via toArray (algorithm_guide.md sec.6.1)
% ---------------------------------------------------------------------------
if ~isempty(fpMorganEth) && ~isempty(fpMorganAsp)
    try
        t_ea = emk.similarity.tanimoto(fpMorganEth, fpMorganAsp);
        t_ae = emk.similarity.tanimoto(fpMorganAsp, fpMorganEth);

        assert(t_ea >= 0 && t_ea <= 1, ...
            sprintf("tanimoto(eth,asp) out of range: %.6f", t_ea));
        assert(t_ea < 1.0, ...
            "tanimoto(eth,asp) should be < 1.0 for different molecules");
        assert(abs(t_ea - t_ae) < 1e-9, ...
            sprintf("tanimoto symmetry broken: T(A,B)=%.10f T(B,A)=%.10f", t_ea, t_ae));

        % Cross-validation via toArray
        a = emk.fingerprint.toArray(fpMorganEth);
        b = emk.fingerprint.toArray(fpMorganAsp);
        tManual = sum(a & b) / sum(a | b);
        assert(abs(t_ea - tManual) < 1e-9, ...
            sprintf("tanimoto cross-validation failed: RDKit=%.10f manual=%.10f", t_ea, tManual));

        logInfo("TC07 PASS  tanimoto(eth,asp)=%.6f  symmetry OK  cross-val OK", t_ea);
        nPass = nPass + 1;
    catch ex
        logError("TC07 FAIL  tanimoto different FPs: %s", ex.message);
        nFail = nFail + 1;
    end
else
    logWarn("TC07 SKIP  Morgan FPs unavailable (TC03 failed)");
    nSkip = nSkip + 1;
end

% ---------------------------------------------------------------------------
% TC08  dice() -- D = 1.0 for same FP; D >= T for different FPs
%   algorithm_guide.md sec.6.2: D(A,B) >= T(A,B) for binary vectors
% ---------------------------------------------------------------------------
if ~isempty(fpMorganEth) && ~isempty(fpMorganAsp)
    try
        dSame = emk.similarity.dice(fpMorganEth, fpMorganEth);
        assert(abs(dSame - 1.0) < 1e-9, ...
            sprintf("dice(same) expected 1.0, got %.10f", dSame));

        d_ea = emk.similarity.dice(fpMorganEth, fpMorganAsp);
        t_ea = emk.similarity.tanimoto(fpMorganEth, fpMorganAsp);
        assert(d_ea >= 0 && d_ea <= 1, ...
            sprintf("dice(eth,asp) out of range: %.6f", d_ea));
        assert(d_ea + 1e-9 >= t_ea, ...
            sprintf("D(%.8f) < T(%.8f): D>=T theorem violated", d_ea, t_ea));

        % Cross-validation: D = 2c/(a+b)
        a = emk.fingerprint.toArray(fpMorganEth);
        b = emk.fingerprint.toArray(fpMorganAsp);
        dManual = 2 * sum(a & b) / (sum(a) + sum(b));
        assert(abs(d_ea - dManual) < 1e-9, ...
            sprintf("dice cross-validation failed: RDKit=%.10f manual=%.10f", d_ea, dManual));

        logInfo("TC08 PASS  dice(same)=%.1f  dice(eth,asp)=%.6f >= tanimoto=%.6f  cross-val OK", ...
            dSame, d_ea, t_ea);
        nPass = nPass + 1;
    catch ex
        logError("TC08 FAIL  dice: %s", ex.message);
        nFail = nFail + 1;
    end
else
    logWarn("TC08 SKIP  Morgan FPs unavailable (TC03 failed)");
    nSkip = nSkip + 1;
end

% ---------------------------------------------------------------------------
% TC09  writeSdf / readSdf round-trip -- 3 molecules, SMILES preserved
%   Canonical SMILES are compared after round-trip (algorithm_guide.md sec.7.1)
% ---------------------------------------------------------------------------
tmpSdf = string(tempname()) + ".sdf";
try
    origMols = {molEthanol, molBenzene, molAspirin};
    emk.io.writeSdf(origMols, tmpSdf);

    assert(isfile(tmpSdf), "writeSdf: output file not created");
    assert(dir(tmpSdf).bytes > 0, "writeSdf: output file is empty");

    readMols = emk.io.readSdf(tmpSdf);

    assert(iscell(readMols), "readSdf must return cell array");
    assert(numel(readMols) == 3, ...
        sprintf("readSdf count mismatch: expected 3, got %d", numel(readMols)));
    assert(isrow(readMols), "readSdf must return 1xN row cell");

    % SMILES round-trip verification
    origSmiles = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"];
    for i = 1:3
        readSmi = string(emk.mol.toSmiles(readMols{i}));
        origSmi = string(emk.mol.toSmiles(origMols{i}));
        assert(readSmi == origSmi, ...
            sprintf("SDF round-trip SMILES mismatch [%d]: orig='%s' read='%s'", ...
            i, origSmi, readSmi));
    end

    if isfile(tmpSdf); delete(tmpSdf); end
    logInfo("TC09 PASS  writeSdf/readSdf round-trip: 3 mols, SMILES preserved");
    nPass = nPass + 1;
catch ex
    if isfile(tmpSdf); delete(tmpSdf); end
    logError("TC09 FAIL  writeSdf/readSdf: %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% TC10  readSmilesList() -- temp SMILES file -> 3 mol objects
%   Verifies: cell count, Python objects, SMILES cross-check
% ---------------------------------------------------------------------------
tmpSmiles = string(tempname()) + ".smi";
try
    % Write a SMILES file with comment and name columns
    fid = fopen(tmpSmiles, "w");
    fprintf(fid, "# Smoke test SMILES file\n");
    fprintf(fid, "CCO\tEthanol\n");
    fprintf(fid, "c1ccccc1\tBenzene\n");
    fprintf(fid, "CC(=O)Oc1ccccc1C(=O)O\tAspirin\n");
    fclose(fid);

    molList = emk.io.readSmilesList(tmpSmiles);

    assert(iscell(molList), "readSmilesList must return cell array");
    assert(numel(molList) == 3, ...
        sprintf("readSmilesList count mismatch: expected 3, got %d", numel(molList)));
    assert(isrow(molList), "readSmilesList must return 1xN row cell");

    for i = 1:3
        assert(startsWith(class(molList{i}), "py."), ...
            sprintf("readSmilesList element %d is not a Python object", i));
    end

    % SMILES cross-validation for aspirin (row 3)
    readAspSmi = string(emk.mol.toSmiles(molList{3}));
    assert(readAspSmi == "CC(=O)Oc1ccccc1C(=O)O", ...
        sprintf("readSmilesList aspirin SMILES mismatch: '%s'", readAspSmi));

    if isfile(tmpSmiles); delete(tmpSmiles); end
    logInfo("TC10 PASS  readSmilesList: 3 mols, Python objects, aspirin SMILES verified");
    nPass = nPass + 1;
catch ex
    if isfile(tmpSmiles); delete(tmpSmiles); end
    logError("TC10 FAIL  readSmilesList: %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% TC11  draw2d() -- ethanol: figure handle returned; image 300x300 RGB uint8
%   algorithm_guide.md sec.8.1
% ---------------------------------------------------------------------------
drawnFig = [];
try
    drawnFig = emk.viz.draw2d(molEthanol);

    assert(isa(drawnFig, "matlab.ui.Figure"), ...
        sprintf("draw2d must return matlab.ui.Figure, got %s", class(drawnFig)));

    % Retrieve image data from the axes
    ax = findall(drawnFig, "type", "axes");
    assert(~isempty(ax), "draw2d figure must have at least one axes");
    img = getimage(ax(1));

    assert(~isempty(img), "draw2d axes must contain image data");
    assert(isa(img, "uint8"), ...
        sprintf("draw2d image must be uint8, got %s", class(img)));
    assert(ndims(img) == 3 && size(img, 3) == 3, ...
        sprintf("draw2d image must be HxWx3 RGB, got %s", mat2str(size(img))));
    assert(size(img, 1) == 300 && size(img, 2) == 300, ...
        sprintf("draw2d default size must be 300x300, got %dx%d", size(img,1), size(img,2)));

    close(drawnFig);
    drawnFig = [];
    logInfo("TC11 PASS  draw2d(ethanol): figure returned, 300x300 RGB uint8 verified");
    nPass = nPass + 1;
catch ex
    if ~isempty(drawnFig) && isvalid(drawnFig)
        try; close(drawnFig); catch; end
    end
    logError("TC11 FAIL  draw2d: %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% TC12  Full pipeline: library similarity search
%   Build Morgan FPs for 5 molecules, rank by Tanimoto similarity to aspirin.
%   Expected: aspirin itself has score 1.0 (rank 1); diverse molecules rank lower.
%   Verifies the complete SMILES -> Mol -> FP -> similarity pipeline end-to-end.
% ---------------------------------------------------------------------------
try
    libSmiles = ["CC(=O)Oc1ccccc1C(=O)O", ... % aspirin (query)
                 "CCO", ...                    % ethanol
                 "c1ccccc1", ...               % benzene
                 "CC(=O)O", ...                % acetic acid
                 "c1ccc(cc1)C(=O)O"];          % benzoic acid
    libNames  = ["aspirin","ethanol","benzene","acetic_acid","benzoic_acid"];

    nLib = numel(libSmiles);
    libMols = cell(1, nLib);
    for i = 1:nLib
        libMols{i} = emk.mol.fromSmiles(libSmiles(i));
    end

    queryFP = emk.fingerprint.morgan(libMols{1});  % aspirin
    scores  = zeros(1, nLib);
    for i = 1:nLib
        libFP      = emk.fingerprint.morgan(libMols{i});
        scores(i)  = emk.similarity.tanimoto(queryFP, libFP);
    end

    % Aspirin vs itself must be 1.0
    assert(abs(scores(1) - 1.0) < 1e-9, ...
        sprintf("aspirin self-similarity expected 1.0, got %.10f", scores(1)));

    % All scores in [0, 1]
    assert(all(scores >= 0 & scores <= 1), ...
        "All similarity scores must be in [0, 1]");

    % Sort by similarity (descending)
    [sortedScores, sortIdx] = sort(scores, "descend");

    logInfo("TC12 PASS  similarity search (query: aspirin):");
    for i = 1:nLib
        logInfo("  Rank %d: %-14s  Tanimoto = %.4f", i, libNames(sortIdx(i)), sortedScores(i));
    end
    nPass = nPass + 1;
catch ex
    logError("TC12 FAIL  full pipeline: %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% Summary
% ---------------------------------------------------------------------------
total = nPass + nFail + nSkip;
logInfo("=== M2 Smoke Test Result: %d PASS / %d FAIL / %d SKIP / %d Total ===", ...
    nPass, nFail, nSkip, total);

if nFail > 0
    error("emk:smoke:test_m2_smoke:failures", ...
        "M2 smoke test completed with %d failure(s). See log above.", nFail);
end
