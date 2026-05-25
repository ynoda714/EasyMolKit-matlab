% test_mvp_smoke.m  End-to-end MVP smoke test (M1-10).
%
%   Covers the full M1 MVP pipeline on both Desktop and MATLAB Online.
%   Also serves as M0-9 verification: emk.util.isOnline() is checked on
%   MATLAB Online to confirm the detection heuristic works correctly.
%
%   Prerequisites (Desktop):
%     1. Run emk.setup.install() once to deploy Embedded Python + RDKit
%     2. addpath(genpath("src"))
%     3. run("tests/smoke/test_mvp_smoke.m")
%
%   Prerequisites (MATLAB Online):
%     1. Run emk.setup.installOnline() once
%     2. addpath(genpath("src"))
%     3. run("tests/smoke/test_mvp_smoke.m")
%
%   Test Cases:
%     TC01  isOnline()          - platform detection (M0-9)
%     TC02  initPython()        - Python environment initialization
%     TC03  verify()            - Python + RDKit availability check
%     TC04  fromSmiles("CCO")   - SMILES -> Mol (ethanol)
%     TC05  fromSmiles (invalid)- error path on invalid SMILES
%     TC06  toSmiles(mol)       - Mol -> canonical SMILES
%     TC07  isValid("CCO")      - valid SMILES -> true
%     TC08  isValid (invalid)   - invalid SMILES -> false
%     TC09  molWeight (ethanol) - MW ~46.07 g/mol
%     TC10  full pipeline       - aspirin: fromSmiles -> molWeight -> toSmiles

addpath(genpath("src"));

nPass = 0;
nFail = 0;
nSkip = 0;

logInfo("=== MVP Smoke Test (M1-10) start ===");

% ---------------------------------------------------------------------------
% TC01  Platform detection (M0-9 verification)
% ---------------------------------------------------------------------------
try
    online = emk.util.isOnline();
    assert(islogical(online) && isscalar(online), ...
        "isOnline must return a logical scalar");
    if online
        logInfo("TC01 PASS  isOnline() = true  (MATLAB Online confirmed -- M0-9 verified)");
    else
        logInfo("TC01 PASS  isOnline() = false  (Desktop)");
    end
    nPass = nPass + 1;
catch ex
    logError("TC01 FAIL  isOnline(): %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% TC02  Python initialization
% ---------------------------------------------------------------------------
try
    emk.setup.initPython();
    pe = pyenv();
    % In OutOfProcess mode, Status stays "NotLoaded" until Python is first
    % invoked.  Check Version instead: it is set by pyenv(Version=...) even
    % before the process has started.
    assert(strlength(string(pe.Version)) > 0, ...
        "pyenv Version should be non-empty after initPython()");
    logInfo("TC02 PASS  initPython() OK  Python %s (mode: %s, status: %s)", ...
        string(pe.Version), string(pe.ExecutionMode), string(pe.Status));
    nPass = nPass + 1;
catch ex
    logError("TC02 FAIL  initPython(): %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% TC03  Python + RDKit environment verification
% ---------------------------------------------------------------------------
try
    result = emk.setup.verify();
    assert(result.python, "verify: Python not available");
    assert(result.rdkit,  "verify: RDKit not available");
    logInfo("TC03 PASS  verify() python=%d rdkit=%d version=%s", ...
        result.python, result.rdkit, result.version);
    nPass = nPass + 1;
catch ex
    logError("TC03 FAIL  verify(): %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% TC04  SMILES parsing: ethanol
% ---------------------------------------------------------------------------
mol = [];
try
    mol = emk.mol.fromSmiles("CCO");
    assert(isa(mol, "py.rdkit.Chem.rdchem.Mol"), ...
        "fromSmiles must return py.rdkit.Chem.rdchem.Mol");
    nAtoms = double(mol.GetNumAtoms());
    assert(nAtoms == 3, ...
        sprintf("ethanol (CCO) must have 3 heavy atoms, got %d", nAtoms));
    logInfo("TC04 PASS  fromSmiles(CCO) -> Mol with %d heavy atoms", nAtoms);
    nPass = nPass + 1;
catch ex
    logError("TC04 FAIL  fromSmiles(CCO): %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% TC05  fromSmiles error path: invalid SMILES
% ---------------------------------------------------------------------------
try
    emk.mol.fromSmiles("INVALID_XYZ_999");
    logError("TC05 FAIL  fromSmiles should throw for invalid SMILES but did not");
    nFail = nFail + 1;
catch ex
    if strcmp(ex.identifier, "emk:mol:fromSmiles:invalidSmiles")
        % Correct: RDKit returned None for unparseable SMILES.
        logInfo("TC05 PASS  fromSmiles correctly threw invalidSmiles");
        nPass = nPass + 1;
    elseif strcmp(ex.identifier, "emk:mol:fromSmiles:rdkitError")
        % rdkitError means the Python/RDKit environment is broken, not a
        % parse failure.  This is a genuine environment problem -> FAIL.
        logError("TC05 FAIL  fromSmiles threw rdkitError (RDKit env broken): %s", ex.message);
        nFail = nFail + 1;
    else
        logError("TC05 FAIL  fromSmiles threw unexpected error: %s  msg: %s", ...
            ex.identifier, ex.message);
        nFail = nFail + 1;
    end
end

% ---------------------------------------------------------------------------
% TC06  Canonical SMILES round-trip  (depends on TC04)
% ---------------------------------------------------------------------------
if ~isempty(mol)
    try
        smi = emk.mol.toSmiles(mol);
        smi = string(smi);
        assert(isstring(smi) && strlength(smi) > 0, ...
            "toSmiles must return non-empty string");
        % RDKit canonical SMILES for ethanol is deterministic: "CCO".
        assert(smi == "CCO", ...
            sprintf("ethanol canonical SMILES expected 'CCO', got '%s'", smi));
        logInfo("TC06 PASS  toSmiles(CCO mol) = '%s'  (canonical verified)", smi);
        nPass = nPass + 1;
    catch ex
        logError("TC06 FAIL  toSmiles(): %s", ex.message);
        nFail = nFail + 1;
    end
else
    logWarn("TC06 SKIP  mol unavailable (TC04 failed)");
    nSkip = nSkip + 1;
end

% ---------------------------------------------------------------------------
% TC07  Molecule validation: valid SMILES
% ---------------------------------------------------------------------------
try
    tf = emk.mol.isValid("CCO");
    assert(tf == true, "isValid(CCO) must return true");
    logInfo("TC07 PASS  isValid(CCO) = true");
    nPass = nPass + 1;
catch ex
    logError("TC07 FAIL  isValid(CCO): %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% TC08  Molecule validation: invalid SMILES
% ---------------------------------------------------------------------------
try
    tf = emk.mol.isValid("INVALID_XYZ_999");
    assert(tf == false, "isValid(invalid) must return false");
    logInfo("TC08 PASS  isValid(INVALID_XYZ_999) = false");
    nPass = nPass + 1;
catch ex
    logError("TC08 FAIL  isValid(invalid): %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% TC09  Molecular weight: ethanol  (depends on TC04)
%   Reference: ethanol C2H6O, IUPAC 2021 avg MW = 46.069 g/mol
%   (algorithm_guide.md sec.4.1, PubChem CID 702)
%   Smoke AbsTol = 0.1 g/mol (unit test uses +/-0.01; looser here to tolerate
%   minor RDKit version differences across platforms without false alarms).
% ---------------------------------------------------------------------------
if ~isempty(mol)
    try
        mw = emk.descriptor.molWeight(mol);
        assert(isnumeric(mw) && isscalar(mw), ...
            "molWeight must return a numeric scalar");
        assert(abs(mw - 46.069) < 0.1, ...
            sprintf("ethanol MW=%.4f expected 46.069 +/-0.1 g/mol", mw));
        logInfo("TC09 PASS  molWeight(CCO) = %.4f g/mol  (ref 46.069)", mw);
        nPass = nPass + 1;
    catch ex
        logError("TC09 FAIL  molWeight(CCO): %s", ex.message);
        nFail = nFail + 1;
    end
else
    logWarn("TC09 SKIP  mol unavailable (TC04 failed)");
    nSkip = nSkip + 1;
end

% ---------------------------------------------------------------------------
% TC10  Full pipeline: aspirin (acetylsalicylic acid)
%   SMILES: CC(=O)Oc1ccccc1C(=O)O
%   Reference MW: 180.159 g/mol  (algorithm_guide.md sec.4.1, IUPAC 2021, PubChem CID 2244)
%   Smoke AbsTol = 0.1 g/mol (same rationale as TC09).
%   Canonical SMILES: RDKit deterministically returns "CC(=O)Oc1ccccc1C(=O)O"
%   (verified on RDKit 2024.03.x; update if canonical form changes with major version).
% ---------------------------------------------------------------------------
try
    aspSmiles    = "CC(=O)Oc1ccccc1C(=O)O";
    aspMol       = emk.mol.fromSmiles(aspSmiles);
    assert(isa(aspMol, "py.rdkit.Chem.rdchem.Mol"), ...
        "aspirin Mol type check failed");

    aspMW = emk.descriptor.molWeight(aspMol);
    assert(isnumeric(aspMW) && isscalar(aspMW), ...
        "aspirin molWeight must return a numeric scalar");
    assert(abs(aspMW - 180.159) < 0.1, ...
        sprintf("aspirin MW=%.4f expected 180.159 +/-0.1 g/mol", aspMW));

    aspSmi = string(emk.mol.toSmiles(aspMol));
    assert(strlength(aspSmi) > 0, "aspirin canonical SMILES must be non-empty");
    assert(aspSmi == "CC(=O)Oc1ccccc1C(=O)O", ...
        sprintf("aspirin canonical SMILES mismatch: '%s'", aspSmi));

    logInfo("TC10 PASS  aspirin: MW=%.4f  smi='%s'  (canonical verified)", aspMW, aspSmi);
    nPass = nPass + 1;
catch ex
    logError("TC10 FAIL  aspirin pipeline: %s", ex.message);
    nFail = nFail + 1;
end

% ---------------------------------------------------------------------------
% Summary
% ---------------------------------------------------------------------------
total = nPass + nFail + nSkip;
logInfo("=== Smoke Test Result: %d PASS / %d FAIL / %d SKIP / %d Total ===", nPass, nFail, nSkip, total);

if nFail > 0
    error("emk:smoke:test_mvp_smoke:failures", ...
        "Smoke test completed with %d failure(s). See log above.", nFail);
end
