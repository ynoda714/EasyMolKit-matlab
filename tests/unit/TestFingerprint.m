classdef TestFingerprint < matlab.unittest.TestCase
% TestFingerprint  Unit tests for emk.fingerprint.morgan and
%                  emk.fingerprint.toArray.
%
% Run with:
%   addpath(genpath("src"));
%   suite = testsuite("tests/unit");           % recommended (full suite)
%   runner = matlab.unittest.TestRunner.withNoPlugins;
%   results = runner.run(suite);
%
% NOTE on test-suite discovery:
%   testsuite("tests/unit/TestFingerprint") fails because MATLAB cannot
%   resolve a path-style string for an individual test class file unless
%   "tests/unit" is already on the MATLAB path.
%   Always use testsuite("tests/unit") (directory form) or add the
%   directory first: addpath("tests/unit"); testsuite("TestFingerprint").
%
% ======================================================================
% Coverage (morgan):
%
%   TC1: Non-Mol input => invalidInput (no RDKit required)
%         double, string, char, cell, logical, empty -- all rejected
%         error message contains offending class name
%   TC2: Return type: Python object whose class contains "BitVect" (requires RDKit)
%   TC3: Default NBits=2048 (requires RDKit)
%   TC4: NBits parameter respected: 1024 / 512 (requires RDKit)
%   TC5: ON-bits sanity: > 0 and <= NBits; aspirin > ethanol (requires RDKit)
%   TC6: Idempotency via toArray: two calls => identical bit arrays (requires RDKit)
%   TC7: Different molecules: different bit arrays (requires RDKit)
%   TC8: Radius effect: ECFP2 != ECFP4, ECFP4 != ECFP6 via toArray (requires RDKit)
%   TC9: Radius=0 is valid (requires RDKit)
%
% ======================================================================
% Coverage (toArray):
%
%   TC10: Non-Python input => invalidInput (no RDKit required)
%   TC11: Python object lacking ToBitString (Mol) => invalidInput (requires RDKit)
%   TC12: Return type: logical row vector (requires RDKit)
%   TC13: numel(toArray(fp)) == fp.GetNumBits() (requires RDKit)
%   TC14: sum(toArray(fp)) == fp.GetNumOnBits() cross-validation (requires RDKit)
%   TC15: Idempotency: two toArray calls => identical arrays (requires RDKit)
%   TC16: Custom NBits flows through to toArray length (requires RDKit)
%   TC17: All elements are logical, no NaN (requires RDKit)
%
% Tests requiring RDKit use assumeTrue(tc.rdkitAvailable()) to skip
% gracefully when Python/RDKit is not configured in the current session.
%
% ======================================================================
% Design note on py.rdkit.DataStructs direct calls:
%   Idempotency / distinctness tests use emk.fingerprint.toArray()-based
%   isequal() comparisons rather than py.rdkit.DataStructs.TanimotoSimilarity:
%     1. Avoids direct py.rdkit.* calls outside emk.* wrappers
%     2. Provides exact (bit-by-bit) rather than threshold-based comparison
%     3. Exercises toArray() as a side effect

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: Input validation -- no RDKit required
        % ------------------------------------------------------------------

        function test_morgan_numericInput_throwsInvalidInput(tc)
        % Numeric input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.fingerprint.morgan(42), ...
                "emk:fingerprint:morgan:invalidInput", ...
                "Numeric input must throw invalidInput");
        end

        function test_morgan_numericInput_errorMessage_containsClass(tc)
        % Error message must contain the offending class name for debugging.
            ME = tc.captureError(@() emk.fingerprint.morgan(42));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain the offending class name 'double'");
        end

        function test_morgan_stringInput_throwsInvalidInput(tc)
        % Passing a raw SMILES string must throw invalidInput.
            tc.verifyError(@() emk.fingerprint.morgan("CCO"), ...
                "emk:fingerprint:morgan:invalidInput", ...
                "String input must throw invalidInput");
        end

        function test_morgan_stringInput_errorMessage_containsClass(tc)
        % Error message for string input must contain "string".
            ME = tc.captureError(@() emk.fingerprint.morgan("CCO"));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "string", ...
                "Error message must contain 'string' class name");
        end

        function test_morgan_charInput_throwsInvalidInput(tc)
        % char literal must throw invalidInput.
            tc.verifyError(@() emk.fingerprint.morgan('CCO'), ...
                "emk:fingerprint:morgan:invalidInput", ...
                "char input must throw invalidInput");
        end

        function test_morgan_cellInput_throwsInvalidInput(tc)
        % Cell array input must throw invalidInput.
            tc.verifyError(@() emk.fingerprint.morgan({"CCO"}), ...
                "emk:fingerprint:morgan:invalidInput", ...
                "Cell array input must throw invalidInput");
        end

        function test_morgan_emptyInput_throwsInvalidInput(tc)
        % Empty matrix input must throw invalidInput.
            tc.verifyError(@() emk.fingerprint.morgan([]), ...
                "emk:fingerprint:morgan:invalidInput", ...
                "Empty matrix input must throw invalidInput");
        end

        function test_morgan_logicalInput_throwsInvalidInput(tc)
        % logical input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.fingerprint.morgan(true), ...
                "emk:fingerprint:morgan:invalidInput", ...
                "logical input must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC2: Return type (requires RDKit)
        % ------------------------------------------------------------------

        function test_morgan_returnIsPythonObject(tc)
        % morgan() must return a Python object (class starts with "py.").
        % Confirms ADR-002 Python reference retention.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyTrue(startsWith(class(fp), "py."), ...
                "morgan must return a Python object (class starts with 'py.')");
        end

        function test_morgan_returnClassContainsBitVect(tc)
        % The class name must contain "BitVect", ruling out arbitrary py.* objects.
        % Stable across RDKit versions that use ExplicitBitVect.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyTrue(contains(class(fp), "BitVect"), ...
                sprintf("Return class must contain 'BitVect', got: %s", class(fp)));
        end

        % ------------------------------------------------------------------
        % TC3: Default NBits = 2048 (requires RDKit)
        % ------------------------------------------------------------------

        function test_morgan_defaultNBits_is2048(tc)
        % Default call must produce a 2048-bit fingerprint.
        % 2048 bits is the EasyMolKit default (algorithm_guide.md 5.1).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyEqual(double(fp.GetNumBits()), 2048, ...
                "Default fingerprint must have 2048 bits");
        end

        % ------------------------------------------------------------------
        % TC4: NBits parameter is honored (requires RDKit)
        % ------------------------------------------------------------------

        function test_morgan_NBits1024_produces1024bits(tc)
        % NBits=1024 must produce a 1024-bit fingerprint.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("c1ccccc1");
            fp  = emk.fingerprint.morgan(mol, NBits=1024);
            tc.verifyEqual(double(fp.GetNumBits()), 1024, ...
                "NBits=1024 must produce a 1024-bit fingerprint");
        end

        function test_morgan_NBits512_produces512bits(tc)
        % NBits=512 must produce a 512-bit fingerprint.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("c1ccccc1");
            fp  = emk.fingerprint.morgan(mol, NBits=512);
            tc.verifyEqual(double(fp.GetNumBits()), 512, ...
                "NBits=512 must produce a 512-bit fingerprint");
        end

        % ------------------------------------------------------------------
        % TC5: ON bits sanity check (requires RDKit)
        % ------------------------------------------------------------------

        function test_morgan_ethanol_numOnBits_isPositive(tc)
        % A valid molecule must produce at least one ON bit.
        % Ethanol (CCO) has 3 heavy atoms; its ECFP4 must be non-zero.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol    = emk.mol.fromSmiles("CCO");
            fp     = emk.fingerprint.morgan(mol);
            numOn  = double(fp.GetNumOnBits());
            tc.verifyGreaterThan(numOn, 0, ...
                "Ethanol ECFP4 must have at least one ON bit");
        end

        function test_morgan_ethanol_numOnBits_notExceedNBits(tc)
        % Number of ON bits must never exceed the total bit length.
        % Defensive guard against bit vector overflow.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.morgan(mol);
            numOn = double(fp.GetNumOnBits());
            nBits = double(fp.GetNumBits());
            tc.verifyLessThanOrEqual(numOn, nBits, ...
                "Number of ON bits must not exceed total bits");
        end

        function test_morgan_aspirin_numOnBits_isPositive(tc)
        % Aspirin (multi-ring molecule) must also produce a non-zero FP.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp    = emk.fingerprint.morgan(mol);
            numOn = double(fp.GetNumOnBits());
            tc.verifyGreaterThan(numOn, 0, ...
                "Aspirin ECFP4 must have at least one ON bit");
        end

        function test_morgan_aspirin_numOnBits_greaterThan_ethanol(tc)
        % Aspirin is structurally richer than ethanol; its ECFP4 must set
        % more bits (larger neighborhood diversity).
        % Validates that ON-bit count scales with molecular complexity.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            molE = emk.mol.fromSmiles("CCO");
            molA = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fpE  = emk.fingerprint.morgan(molE);
            fpA  = emk.fingerprint.morgan(molA);
            tc.verifyGreaterThan(double(fpA.GetNumOnBits()), ...
                double(fpE.GetNumOnBits()), ...
                "Aspirin ECFP4 must have more ON bits than ethanol ECFP4");
        end

        % ------------------------------------------------------------------
        % TC6: Idempotency via toArray -- two calls => identical bit arrays
        % (requires RDKit)
        % Note: isequal(bits1, bits2) is stronger than Tanimoto == 1.0 because
        %       Tanimoto is defined as 1.0 for two all-zero vectors, while
        %       isequal would still correctly pass only when bits are identical.
        % ------------------------------------------------------------------

        function test_morgan_ethanol_idempotency_arrayEqual(tc)
        % Calling morgan() twice on the same mol must produce identical arrays.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp1   = emk.fingerprint.morgan(mol);
            fp2   = emk.fingerprint.morgan(mol);
            bits1 = emk.fingerprint.toArray(fp1);
            bits2 = emk.fingerprint.toArray(fp2);
            tc.verifyTrue(isequal(bits1, bits2), ...
                "Same molecule called twice must produce identical bit arrays");
        end

        % ------------------------------------------------------------------
        % TC7: Different molecules produce different bit arrays
        % (requires RDKit)
        % ------------------------------------------------------------------

        function test_morgan_ethanol_vs_benzene_arrayDiffer(tc)
        % Ethanol and benzene are structurally distinct; bit arrays must differ.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            molE  = emk.mol.fromSmiles("CCO");
            molB  = emk.mol.fromSmiles("c1ccccc1");
            bitsE = emk.fingerprint.toArray(emk.fingerprint.morgan(molE));
            bitsB = emk.fingerprint.toArray(emk.fingerprint.morgan(molB));
            tc.verifyFalse(isequal(bitsE, bitsB), ...
                "Ethanol and benzene must produce different bit arrays");
        end

        function test_morgan_ethanol_vs_aspirin_arrayDiffer(tc)
        % Ethanol and aspirin (aromatic + ester) must produce different arrays.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            molE  = emk.mol.fromSmiles("CCO");
            molA  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            bitsE = emk.fingerprint.toArray(emk.fingerprint.morgan(molE));
            bitsA = emk.fingerprint.toArray(emk.fingerprint.morgan(molA));
            tc.verifyFalse(isequal(bitsE, bitsA), ...
                "Ethanol and aspirin must produce different bit arrays");
        end

        % ------------------------------------------------------------------
        % TC8: Radius parameter effect via toArray (requires RDKit)
        % Aspirin is used because its rich structure ensures different radii
        % capture distinct neighborhoods.  See algorithm_guide.md 5.1.
        % ------------------------------------------------------------------

        function test_morgan_radius1_vs_radius2_aspirin_arrayDiffer(tc)
        % ECFP2 (radius=1) and ECFP4 (radius=2) bit arrays must differ.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            bits1 = emk.fingerprint.toArray(emk.fingerprint.morgan(mol, Radius=1));
            bits2 = emk.fingerprint.toArray(emk.fingerprint.morgan(mol, Radius=2));
            tc.verifyFalse(isequal(bits1, bits2), ...
                "ECFP2 and ECFP4 of aspirin must produce different bit arrays");
        end

        function test_morgan_radius2_vs_radius3_aspirin_arrayDiffer(tc)
        % ECFP4 (radius=2) vs ECFP6 (radius=3) bit arrays must differ.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            bits2 = emk.fingerprint.toArray(emk.fingerprint.morgan(mol, Radius=2));
            bits3 = emk.fingerprint.toArray(emk.fingerprint.morgan(mol, Radius=3));
            tc.verifyFalse(isequal(bits2, bits3), ...
                "ECFP4 and ECFP6 of aspirin must produce different bit arrays");
        end

        % ------------------------------------------------------------------
        % TC9: Radius=0 is valid (requires RDKit)
        % ------------------------------------------------------------------

        function test_morgan_radius0_doesNotThrow(tc)
        % Radius=0 computes single-atom environment fingerprints.
        % This is a degenerate but valid use case (atom-count encoding).
        % RDKit accepts Radius=0 without error.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tc.verifyWarningFree(@() emk.fingerprint.morgan(mol, Radius=0), ...
                "Radius=0 must not produce a warning");
        end

        function test_morgan_radius0_returns2048bits(tc)
        % Radius=0 with default NBits must still return a 2048-bit FP.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            fp   = emk.fingerprint.morgan(mol, Radius=0);
            tc.verifyEqual(double(fp.GetNumBits()), 2048, ...
                "Radius=0 must produce a 2048-bit fingerprint");
        end

        % ==================================================================
        % TC10: toArray input validation -- no RDKit required
        % ==================================================================

        function test_toArray_numericInput_throwsInvalidInput(tc)
        % Non-Python numeric input must throw invalidInput without RDKit.
            tc.verifyError(@() emk.fingerprint.toArray(42), ...
                "emk:fingerprint:toArray:invalidInput", ...
                "Numeric input must throw invalidInput");
        end

        function test_toArray_numericInput_errorMessage_containsClass(tc)
        % Error message must name the offending type.
            ME = tc.captureError(@() emk.fingerprint.toArray(42));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain 'double'");
        end

        function test_toArray_stringInput_throwsInvalidInput(tc)
        % String input must throw invalidInput.
            tc.verifyError(@() emk.fingerprint.toArray("CCO"), ...
                "emk:fingerprint:toArray:invalidInput", ...
                "String input must throw invalidInput");
        end

        function test_toArray_charInput_throwsInvalidInput(tc)
        % char literal must throw invalidInput.
            tc.verifyError(@() emk.fingerprint.toArray('abc'), ...
                "emk:fingerprint:toArray:invalidInput", ...
                "char input must throw invalidInput");
        end

        function test_toArray_logicalInput_throwsInvalidInput(tc)
        % logical input must throw invalidInput.
            tc.verifyError(@() emk.fingerprint.toArray(true), ...
                "emk:fingerprint:toArray:invalidInput", ...
                "logical input must throw invalidInput");
        end

        function test_toArray_emptyInput_throwsInvalidInput(tc)
        % Empty matrix must throw invalidInput.
            tc.verifyError(@() emk.fingerprint.toArray([]), ...
                "emk:fingerprint:toArray:invalidInput", ...
                "Empty matrix input must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC11: toArray with wrong Python type (requires RDKit)
        % A Mol object starts with "py." so it passes the first check,
        % but lacks ToBitString() => caught and rethrown as invalidInput.
        % Guards against accidental mol/fp argument confusion.
        % ------------------------------------------------------------------

        function test_toArray_molInput_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tc.verifyError(@() emk.fingerprint.toArray(mol), ...
                "emk:fingerprint:toArray:invalidInput", ...
                "Passing a Mol to toArray must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC12: toArray return type and shape (requires RDKit)
        % ------------------------------------------------------------------

        function test_toArray_returnsLogical(tc)
        % toArray must return a logical array, not double or char.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            fp   = emk.fingerprint.morgan(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyClass(bits, "logical", ...
                "toArray must return a logical array");
        end

        function test_toArray_returnsRowVector(tc)
        % toArray must return a row vector (1 x N), not a column.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            fp   = emk.fingerprint.morgan(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyEqual(size(bits, 1), 1, ...
                "toArray must return a row vector (1 x N)");
        end

        % ------------------------------------------------------------------
        % TC13: toArray length == GetNumBits() (requires RDKit)
        % ------------------------------------------------------------------

        function test_toArray_defaultFP_length2048(tc)
        % Default morgan FP -> toArray length must be 2048.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            fp   = emk.fingerprint.morgan(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyEqual(numel(bits), 2048, ...
                "toArray of 2048-bit FP must have 2048 elements");
        end

        function test_toArray_length_matchesGetNumBits(tc)
        % numel(toArray(fp)) must equal double(fp.GetNumBits()).
        % Cross-validates MATLAB array length against the Python API.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            fp   = emk.fingerprint.morgan(mol, NBits=1024);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyEqual(numel(bits), double(fp.GetNumBits()), ...
                "numel(toArray(fp)) must equal fp.GetNumBits()");
        end

        % ------------------------------------------------------------------
        % TC14: sum(toArray(fp)) == fp.GetNumOnBits() cross-validation
        % (requires RDKit)
        %
        % This is the primary consistency check for toArray correctness:
        % every ON bit in the Python object must appear as a true element
        % in the MATLAB array.  Any mismatch indicates a conversion error.
        % ------------------------------------------------------------------

        function test_toArray_ethanol_sumEqualGetNumOnBits(tc)
        % sum(toArray(fp)) must exactly match fp.GetNumOnBits() for ethanol.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            fp   = emk.fingerprint.morgan(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyEqual(sum(bits), double(fp.GetNumOnBits()), ...
                "sum(toArray(fp)) must exactly match fp.GetNumOnBits() for ethanol");
        end

        function test_toArray_aspirin_sumEqualGetNumOnBits(tc)
        % Same cross-validation for aspirin (more ON bits than ethanol).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp   = emk.fingerprint.morgan(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyEqual(sum(bits), double(fp.GetNumOnBits()), ...
                "sum(toArray(fp)) must exactly match fp.GetNumOnBits() for aspirin");
        end

        function test_toArray_benzene_sumEqualGetNumOnBits(tc)
        % Same cross-validation for benzene (aromatic ring).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccccc1");
            fp   = emk.fingerprint.morgan(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyEqual(sum(bits), double(fp.GetNumOnBits()), ...
                "sum(toArray(fp)) must exactly match fp.GetNumOnBits() for benzene");
        end

        % ------------------------------------------------------------------
        % TC15: toArray idempotency (requires RDKit)
        % ------------------------------------------------------------------

        function test_toArray_idempotency_sameResultTwice(tc)
        % Two toArray calls on the same fp must give identical arrays.
        % Confirms conversion is deterministic and does not mutate the object.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.morgan(mol);
            bits1 = emk.fingerprint.toArray(fp);
            bits2 = emk.fingerprint.toArray(fp);
            tc.verifyTrue(isequal(bits1, bits2), ...
                "Two toArray calls on the same fp must produce identical arrays");
        end

        % ------------------------------------------------------------------
        % TC16: Custom NBits flows through to toArray length (requires RDKit)
        % ------------------------------------------------------------------

        function test_toArray_NBits512_length512(tc)
        % End-to-end: morgan(..., NBits=512) -> fp -> toArray -> numel==512.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccccc1");
            fp   = emk.fingerprint.morgan(mol, NBits=512);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyEqual(numel(bits), 512, ...
                "toArray of 512-bit FP must have 512 elements");
        end

        % ------------------------------------------------------------------
        % TC17: toArray output values are logical {0,1}, no NaN (requires RDKit)
        % ------------------------------------------------------------------

        function test_toArray_allBitsAreLogical(tc)
        % Every element must be true or false (islogical + unique values).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp   = emk.fingerprint.morgan(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyTrue(islogical(bits), ...
                "toArray output must be logical class");
            uniqueVals = unique(double(bits));
            tc.verifyTrue(all(ismember(uniqueVals, [0, 1])), ...
                "toArray output must contain only 0 and 1 values");
        end

        function test_toArray_noNaNValues(tc)
        % toArray must not produce NaN values.
        % double(logical) is always 0 or 1; NaN cannot appear.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            fp   = emk.fingerprint.morgan(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyFalse(any(isnan(double(bits))), ...
                "toArray must not produce NaN values");
        end

        % ==================================================================
        % maccs -- TC-M1: Input validation (no RDKit required)
        % ==================================================================

        function test_maccs_numericInput_throwsInvalidInput(tc)
        % Numeric input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.fingerprint.maccs(42), ...
                "emk:fingerprint:maccs:invalidInput", ...
                "Numeric input must throw invalidInput");
        end

        function test_maccs_numericInput_errorMessage_containsClass(tc)
        % Error message must contain the offending class name.
            ME = tc.captureError(@() emk.fingerprint.maccs(42));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain 'double'");
        end

        function test_maccs_stringInput_throwsInvalidInput(tc)
        % Passing a raw SMILES string must throw invalidInput.
            tc.verifyError(@() emk.fingerprint.maccs("CCO"), ...
                "emk:fingerprint:maccs:invalidInput", ...
                "String input must throw invalidInput");
        end

        function test_maccs_charInput_throwsInvalidInput(tc)
        % char literal must throw invalidInput.
            tc.verifyError(@() emk.fingerprint.maccs('CCO'), ...
                "emk:fingerprint:maccs:invalidInput", ...
                "char input must throw invalidInput");
        end

        function test_maccs_logicalInput_throwsInvalidInput(tc)
        % logical input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.fingerprint.maccs(true), ...
                "emk:fingerprint:maccs:invalidInput", ...
                "logical input must throw invalidInput");
        end

        function test_maccs_emptyInput_throwsInvalidInput(tc)
        % Empty matrix input must throw invalidInput.
            tc.verifyError(@() emk.fingerprint.maccs([]), ...
                "emk:fingerprint:maccs:invalidInput", ...
                "Empty matrix input must throw invalidInput");
        end

        function test_maccs_cellInput_throwsInvalidInput(tc)
        % Cell array input must throw invalidInput.
            tc.verifyError(@() emk.fingerprint.maccs({"CCO"}), ...
                "emk:fingerprint:maccs:invalidInput", ...
                "Cell array input must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC-M2: Return type is a Python ExplicitBitVect (requires RDKit)
        % ------------------------------------------------------------------

        function test_maccs_returnIsPythonObject(tc)
        % maccs() must return a Python object (class starts with "py.").
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.maccs(mol);
            tc.verifyTrue(startsWith(class(fp), "py."), ...
                "maccs must return a Python object (class starts with 'py.')");
        end

        function test_maccs_returnClassContainsBitVect(tc)
        % The class name must contain "BitVect".
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.maccs(mol);
            tc.verifyTrue(contains(class(fp), "BitVect"), ...
                sprintf("Return class must contain 'BitVect', got: %s", class(fp)));
        end

        % ------------------------------------------------------------------
        % TC-M3: Bit length = 167 (requires RDKit)
        % MACCS keys are 166 public keys; RDKit returns 167 bits with
        % bit 0 always 0 (unused).
        % ------------------------------------------------------------------

        function test_maccs_ethanol_numBits_is167(tc)
        % MACCS must always produce a 167-bit fingerprint.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.maccs(mol);
            tc.verifyEqual(double(fp.GetNumBits()), 167, ...
                "MACCS fingerprint must have 167 bits");
        end

        function test_maccs_benzene_numBits_is167(tc)
        % 167 bits is fixed regardless of molecule.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("c1ccccc1");
            fp  = emk.fingerprint.maccs(mol);
            tc.verifyEqual(double(fp.GetNumBits()), 167, ...
                "MACCS must produce 167 bits for benzene");
        end

        % ------------------------------------------------------------------
        % TC-M4: ON bits sanity (requires RDKit)
        % ------------------------------------------------------------------

        function test_maccs_ethanol_numOnBits_isPositive(tc)
        % A valid molecule must produce at least one ON bit.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.maccs(mol);
            numOn = double(fp.GetNumOnBits());
            tc.verifyGreaterThan(numOn, 0, ...
                "Ethanol MACCS must have at least one ON bit");
        end

        function test_maccs_ethanol_numOnBits_notExceedNBits(tc)
        % Number of ON bits must not exceed total bit length.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.maccs(mol);
            numOn = double(fp.GetNumOnBits());
            nBits = double(fp.GetNumBits());
            tc.verifyLessThanOrEqual(numOn, nBits, ...
                "Number of ON bits must not exceed total bits");
        end

        function test_maccs_benzene_numOnBits_isPositive(tc)
        % Benzene (aromatic ring) must produce at least one ON bit.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("c1ccccc1");
            fp    = emk.fingerprint.maccs(mol);
            numOn = double(fp.GetNumOnBits());
            tc.verifyGreaterThan(numOn, 0, ...
                "Benzene MACCS must have at least one ON bit");
        end

        % ------------------------------------------------------------------
        % TC-M5: Aspirin has more ON bits than ethanol (requires RDKit)
        % Aspirin (C9H8O4: aromatic ring, ester, carboxyl, carbonyl) activates
        % many more MACCS structural keys than simple ethanol (aliphatic C-O).
        % Benzene was not used because MACCS key counts are not guaranteed
        % strictly greater for benzene vs ethanol (pattern-based counts depend
        % on which simple aliphatic keys ethanol activates).
        % ------------------------------------------------------------------

        function test_maccs_aspirin_numOnBits_greaterThan_ethanol(tc)
        % Aspirin (multi-functional) must activate more MACCS keys than ethanol.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            molE = emk.mol.fromSmiles("CCO");
            molA = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fpE  = emk.fingerprint.maccs(molE);
            fpA  = emk.fingerprint.maccs(molA);
            tc.verifyGreaterThan(double(fpA.GetNumOnBits()), ...
                double(fpE.GetNumOnBits()), ...
                "Aspirin MACCS must have more ON bits than ethanol MACCS");
        end

        % ------------------------------------------------------------------
        % TC-M6: Idempotency via toArray (requires RDKit)
        % ------------------------------------------------------------------

        function test_maccs_ethanol_idempotency_arrayEqual(tc)
        % Two calls on the same mol must produce identical MACCS arrays.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp1   = emk.fingerprint.maccs(mol);
            fp2   = emk.fingerprint.maccs(mol);
            bits1 = emk.fingerprint.toArray(fp1);
            bits2 = emk.fingerprint.toArray(fp2);
            tc.verifyTrue(isequal(bits1, bits2), ...
                "Same molecule called twice must produce identical MACCS arrays");
        end

        % ------------------------------------------------------------------
        % TC-M7: Different molecules produce different arrays (requires RDKit)
        % ------------------------------------------------------------------

        function test_maccs_ethanol_vs_benzene_arrayDiffer(tc)
        % Ethanol and benzene must produce different MACCS arrays.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            molE  = emk.mol.fromSmiles("CCO");
            molB  = emk.mol.fromSmiles("c1ccccc1");
            bitsE = emk.fingerprint.toArray(emk.fingerprint.maccs(molE));
            bitsB = emk.fingerprint.toArray(emk.fingerprint.maccs(molB));
            tc.verifyFalse(isequal(bitsE, bitsB), ...
                "Ethanol and benzene must produce different MACCS arrays");
        end

        % ------------------------------------------------------------------
        % TC-M8: toArray cross-validation: sum == GetNumOnBits (requires RDKit)
        % ------------------------------------------------------------------

        function test_maccs_ethanol_toArray_sumEqualGetNumOnBits(tc)
        % sum(toArray(fp)) must equal fp.GetNumOnBits() for MACCS FP.
        % This is the primary cross-validation: MATLAB bits match Python object.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            fp   = emk.fingerprint.maccs(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyEqual(sum(bits), double(fp.GetNumOnBits()), ...
                "sum(toArray(maccs(mol))) must equal fp.GetNumOnBits()");
        end

        function test_maccs_aspirin_toArray_sumEqualGetNumOnBits(tc)
        % Cross-validation for aspirin (more complex molecule).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp   = emk.fingerprint.maccs(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyEqual(sum(bits), double(fp.GetNumOnBits()), ...
                "sum(toArray(maccs(mol))) must equal fp.GetNumOnBits() for aspirin");
        end

        % ------------------------------------------------------------------
        % TC-M9: toArray length = 167 (requires RDKit)
        % Verifies end-to-end interoperability of maccs() and toArray().
        % ------------------------------------------------------------------

        function test_maccs_toArray_length167(tc)
        % End-to-end: maccs() -> toArray() -> numel == 167.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccccc1");
            fp   = emk.fingerprint.maccs(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyEqual(numel(bits), 167, ...
                "toArray(maccs(mol)) must have 167 elements");
        end

        function test_maccs_toArray_returnsLogical(tc)
        % toArray(maccs(fp)) must return a logical array, not double or char.
        % Verifies that ExplicitBitVect from maccs() is correctly handled
        % by toArray() -- the same ExplicitBitVect type as morgan().
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            fp   = emk.fingerprint.maccs(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyClass(bits, "logical", ...
                "toArray(maccs(mol)) must return a logical array");
        end

        function test_maccs_toArray_returnsRowVector(tc)
        % toArray(maccs(fp)) must return a row vector (1 x 167).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            fp   = emk.fingerprint.maccs(mol);
            bits = emk.fingerprint.toArray(fp);
            tc.verifyEqual(size(bits, 1), 1, ...
                "toArray(maccs(mol)) must return a row vector (1 x 167)");
        end

    end

    % ======================================================================
    methods (Access = private)

        function tf = rdkitAvailable(~)
        % Return true if Python is configured and RDKit can be imported.
        %
        % In MATLAB OutOfProcess mode, pyenv().Status stays "NotLoaded" until
        % the first Python API call (the process starts on demand).  Checking
        % Status and returning early would always produce false before first use.
        % Instead, ensure pyenv is configured via initPython() when NotLoaded,
        % then attempt the import directly -- matching verify() behaviour.
            tf = false;
            try
                pe = pyenv();
                if strcmp(string(pe.Status), "NotLoaded")
                    emk.setup.initPython();
                end
                py.importlib.import_module("rdkit.Chem");
                tf = true;
            catch
                % Python not configured or RDKit not installed
            end
        end

        function ME = captureError(~, fcn)
        % Call fcn and return the MException if thrown, or [] on success.
            ME = [];
            try
                fcn();
            catch e
                ME = e;
            end
        end

    end
end
