classdef TestSimilarity < matlab.unittest.TestCase
% TestSimilarity  Unit tests for emk.similarity.tanimoto and dice.
%
% Run with:
%   addpath(genpath("src"));
%   suite = testsuite("tests/unit");
%   runner = matlab.unittest.TestRunner.withNoPlugins;
%   results = runner.run(suite);
%
% ======================================================================
% Coverage (tanimoto):
%
%   TC1:  Non-Python fp1 => invalidInput (no RDKit required)
%          double, string, char, logical, empty -- all rejected
%          error message contains offending class name
%   TC2:  Non-Python fp2 => invalidInput (requires RDKit for valid fp1)
%          same types as TC1; all with assumeTrue(rdkitAvailable)
%   TC3:  Return type: double scalar, finite, isreal (requires RDKit)
%   TC4:  Identical fingerprints => score == 1.0 (requires RDKit)
%          Same molecule morgan, same molecule maccs
%   TC5:  (reserved -- all-zero pair not reliably constructible with RDKit)
%   TC6:  Score symmetry: tanimoto(fp1,fp2) == tanimoto(fp2,fp1) (requires RDKit)
%   TC7:  Different molecules => score < 1.0 (requires RDKit)
%   TC8:  Score in [0, 1] range for various molecule pairs (requires RDKit)
%   TC9:  Manual cross-validation: score == c/(a+b-c) via toArray (requires RDKit)
%          Primary correctness test: MATLAB computation vs RDKit API
%  TC10:  MACCS fingerprints also accepted (requires RDKit)
%  TC11:  Mol object passed as fp => rdkitError (requires RDKit)
%  TC12:  Bit-length mismatch (Morgan 2048 vs MACCS 167) => rdkitError (requires RDKit)
%
% Tests requiring RDKit use assumeTrue(tc.rdkitAvailable()) to skip
% gracefully when Python/RDKit is not configured.
%
% ======================================================================
% Coverage (dice):
%
%   TC1d: Non-Python fp1 => invalidInput (no RDKit required)
%         + numeric, string, char, logical, empty [] rejected
%   TC2d: Non-Python fp2 => invalidInput (requires RDKit for valid fp1)
%   TC3d: Return type: double scalar, finite, isreal (requires RDKit)
%   TC4d: Identical fingerprints => score == 1.0 (requires RDKit)
%   TC5d: Score in [0, 1] range — ethanol/aspirin and ethanol/benzene (requires RDKit)
%   TC6d: Symmetry: dice(fp1,fp2) == dice(fp2,fp1) (requires RDKit)
%   TC7d: Different molecules => score < 1.0 (requires RDKit)
%   TC8d: Cross-validation: D = 2c/(a+b) via toArray — Morgan AND MACCS, two pairs (requires RDKit)
%   TC9d: Dice >= Tanimoto — two molecule pairs (requires RDKit)
%  TC10d: Mol object passed as fp => rdkitError (requires RDKit)
%  TC11d: Bit-length mismatch => rdkitError (requires RDKit)
%
% ======================================================================
% Design note on cross-validation (TC9):
%   Manual Tanimoto from toArray(): T = sum(a&b) / sum(a|b)
%   This verifies that the RDKit API and the MATLAB logical-array
%   calculation agree to within floating-point tolerance (AbsTol=1e-10).
%   Any bit-encoding mismatch or RDKit API change will be caught here.

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: fp1 input validation -- no RDKit required
        % ------------------------------------------------------------------

        function test_tanimoto_fp1_numericInput_throwsInvalidInput(tc)
        % Numeric fp1 must throw invalidInput before any RDKit call.
        % Both arguments are non-Python; fp1 check fires first.
            tc.verifyError(@() emk.similarity.tanimoto(42, 42), ...
                "emk:similarity:tanimoto:invalidInput", ...
                "Numeric fp1 must throw invalidInput");
        end

        function test_tanimoto_fp1_numericInput_errorMessage_containsClass(tc)
        % Error message must name the offending type for fp1.
            ME = tc.captureError(@() emk.similarity.tanimoto(42, 42));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain 'double' for numeric fp1");
        end

        function test_tanimoto_fp1_stringInput_throwsInvalidInput(tc)
        % String fp1 must throw invalidInput.
            tc.verifyError(@() emk.similarity.tanimoto("CCO", "CCO"), ...
                "emk:similarity:tanimoto:invalidInput", ...
                "String fp1 must throw invalidInput");
        end

        function test_tanimoto_fp1_charInput_throwsInvalidInput(tc)
        % char fp1 must throw invalidInput.
            tc.verifyError(@() emk.similarity.tanimoto('CCO', 'CCO'), ...
                "emk:similarity:tanimoto:invalidInput", ...
                "char fp1 must throw invalidInput");
        end

        function test_tanimoto_fp1_logicalInput_throwsInvalidInput(tc)
        % logical fp1 must throw invalidInput.
            tc.verifyError(@() emk.similarity.tanimoto(true, true), ...
                "emk:similarity:tanimoto:invalidInput", ...
                "logical fp1 must throw invalidInput");
        end

        function test_tanimoto_fp1_emptyInput_throwsInvalidInput(tc)
        % Empty matrix fp1 must throw invalidInput.
            tc.verifyError(@() emk.similarity.tanimoto([], []), ...
                "emk:similarity:tanimoto:invalidInput", ...
                "Empty matrix fp1 must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC2: fp2 input validation -- no RDKit required
        % (fp1 must be a Python object so fp2 is the failing argument)
        % ------------------------------------------------------------------

        function test_tanimoto_fp2_numericInput_throwsInvalidInput(tc)
        % When fp1 is valid Python but fp2 is numeric, must throw invalidInput.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp1 = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.similarity.tanimoto(fp1, 42), ...
                "emk:similarity:tanimoto:invalidInput", ...
                "Numeric fp2 with valid fp1 must throw invalidInput");
        end

        function test_tanimoto_fp2_numericInput_errorMessage_containsClass(tc)
        % Error message must name the offending type for fp2.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp1 = emk.fingerprint.morgan(mol);
            ME  = tc.captureError(@() emk.similarity.tanimoto(fp1, 42));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain 'double' for numeric fp2");
        end

        function test_tanimoto_fp2_stringInput_throwsInvalidInput(tc)
        % String fp2 with valid fp1 must throw invalidInput.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp1 = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.similarity.tanimoto(fp1, "CCO"), ...
                "emk:similarity:tanimoto:invalidInput", ...
                "String fp2 must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC3: Return type is double scalar (requires RDKit)
        % ------------------------------------------------------------------

        function test_tanimoto_returnsDoubleScalar(tc)
        % tanimoto() must return a double scalar.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.morgan(mol);
            score = emk.similarity.tanimoto(fp, fp);
            tc.verifyClass(score, "double", ...
                "tanimoto must return a double");
            tc.verifySize(score, [1, 1], ...
                "tanimoto must return a scalar");
        end

        function test_tanimoto_returnIsFinite(tc)
        % Return value must be finite (not NaN or Inf).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.morgan(mol);
            score = emk.similarity.tanimoto(fp, fp);
            tc.verifyTrue(isfinite(score), ...
                "tanimoto must return a finite value");
        end

        function test_tanimoto_returnIsReal(tc)
        % Return value must be real-valued (not complex).
        % Complex output would indicate a conversion error in double().
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.morgan(mol);
            score = emk.similarity.tanimoto(fp, fp);
            tc.verifyTrue(isreal(score), ...
                "tanimoto must return a real-valued (non-complex) double");
        end

        % ------------------------------------------------------------------
        % TC4: Identical fingerprints => score == 1.0 (requires RDKit)
        % ------------------------------------------------------------------

        function test_tanimoto_sameFP_morgan_equals1(tc)
        % Tanimoto of a fingerprint with itself must be exactly 1.0.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.morgan(mol);
            score = emk.similarity.tanimoto(fp, fp);
            tc.verifyEqual(score, 1.0, ...
                "Self-similarity (Morgan) must be exactly 1.0");
        end

        function test_tanimoto_sameMol_twoCalls_morgan_equals1(tc)
        % Two separate morgan() calls on same mol must give score = 1.0.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1   = emk.fingerprint.morgan(mol);
            fp2   = emk.fingerprint.morgan(mol);
            score = emk.similarity.tanimoto(fp1, fp2);
            tc.verifyEqual(score, 1.0, ...
                "Identical morgan FPs from same mol must give score=1.0");
        end

        function test_tanimoto_sameFP_maccs_equals1(tc)
        % Tanimoto of a MACCS fingerprint with itself must be exactly 1.0.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("c1ccccc1");
            fp    = emk.fingerprint.maccs(mol);
            score = emk.similarity.tanimoto(fp, fp);
            tc.verifyEqual(score, 1.0, ...
                "Self-similarity (MACCS) must be exactly 1.0");
        end

        % ------------------------------------------------------------------
        % TC5: Score in [0, 1] for various pairs (requires RDKit)
        % ------------------------------------------------------------------

        function test_tanimoto_ethanolVsAspirin_inRange(tc)
        % Tanimoto score must be in [0, 1] for ethanol vs aspirin.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            score = emk.similarity.tanimoto(fp1, fp2);
            tc.verifyGreaterThanOrEqual(score, 0.0, ...
                "Tanimoto score must be >= 0");
            tc.verifyLessThanOrEqual(score, 1.0, ...
                "Tanimoto score must be <= 1");
        end

        function test_tanimoto_ethanolVsBenzene_inRange(tc)
        % Tanimoto score must be in [0, 1] for ethanol vs benzene.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("c1ccccc1");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            score = emk.similarity.tanimoto(fp1, fp2);
            tc.verifyGreaterThanOrEqual(score, 0.0, ...
                "Tanimoto score must be >= 0");
            tc.verifyLessThanOrEqual(score, 1.0, ...
                "Tanimoto score must be <= 1");
        end

        % ------------------------------------------------------------------
        % TC6: Symmetry: tanimoto(fp1, fp2) == tanimoto(fp2, fp1) (requires RDKit)
        % ------------------------------------------------------------------

        function test_tanimoto_symmetry_ethanolVsAspirin(tc)
        % tanimoto(fp1, fp2) must equal tanimoto(fp2, fp1).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1   = emk.mol.fromSmiles("CCO");
            mol2   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1    = emk.fingerprint.morgan(mol1);
            fp2    = emk.fingerprint.morgan(mol2);
            score1 = emk.similarity.tanimoto(fp1, fp2);
            score2 = emk.similarity.tanimoto(fp2, fp1);
            tc.verifyEqual(score1, score2, ...
                "Tanimoto must be symmetric: T(A,B) == T(B,A)");
        end

        function test_tanimoto_symmetry_maccs(tc)
        % Symmetry also holds for MACCS fingerprints.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1   = emk.mol.fromSmiles("CCO");
            mol2   = emk.mol.fromSmiles("c1ccccc1");
            fp1    = emk.fingerprint.maccs(mol1);
            fp2    = emk.fingerprint.maccs(mol2);
            score1 = emk.similarity.tanimoto(fp1, fp2);
            score2 = emk.similarity.tanimoto(fp2, fp1);
            tc.verifyEqual(score1, score2, ...
                "Tanimoto symmetry must hold for MACCS fingerprints");
        end

        % ------------------------------------------------------------------
        % TC7: Different molecules => score < 1.0 (requires RDKit)
        % ------------------------------------------------------------------

        function test_tanimoto_ethanolVsAspirin_lessThan1(tc)
        % Structurally different molecules must have score < 1.0.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            score = emk.similarity.tanimoto(fp1, fp2);
            tc.verifyLessThan(score, 1.0, ...
                "Ethanol vs aspirin Tanimoto must be less than 1.0");
        end

        function test_tanimoto_ethanolVsBenzene_lessThan1(tc)
        % Ethanol and benzene (aliphatic vs aromatic) must have score < 1.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("c1ccccc1");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            score = emk.similarity.tanimoto(fp1, fp2);
            tc.verifyLessThan(score, 1.0, ...
                "Ethanol vs benzene Tanimoto must be less than 1.0");
        end

        % ------------------------------------------------------------------
        % TC8: Score is non-negative (requires RDKit)
        % ------------------------------------------------------------------

        function test_tanimoto_ethanolVsAspirin_nonNegative(tc)
        % Tanimoto score must be >= 0.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            score = emk.similarity.tanimoto(fp1, fp2);
            tc.verifyGreaterThanOrEqual(score, 0.0, ...
                "Tanimoto score must be non-negative");
        end

        % ------------------------------------------------------------------
        % TC9: Manual cross-validation via toArray (requires RDKit)
        %
        % This is the primary correctness test.
        % Compute Tanimoto manually from MATLAB logical arrays and compare
        % with the RDKit API value.  Any mismatch indicates an encoding bug
        % in the Python-to-MATLAB conversion or a wrong API call.
        %
        % T(A,B) = sum(a & b) / sum(a | b)
        % where a = toArray(fp1), b = toArray(fp2)
        % AbsTol = 1e-10 (floating-point arithmetic only, no rounding)
        % ------------------------------------------------------------------

        function test_tanimoto_crossValidation_ethanolVsAspirin(tc)
        % Manual T from toArray must match RDKit TanimotoSimilarity.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            % RDKit API value
            rdkitScore = emk.similarity.tanimoto(fp1, fp2);
            % Manual calculation via toArray
            a = emk.fingerprint.toArray(fp1);
            b = emk.fingerprint.toArray(fp2);
            c = sum(a & b);
            manualScore = c / sum(a | b);
            tc.verifyEqual(rdkitScore, manualScore, "AbsTol", 1e-10, ...
                "Manual Tanimoto from toArray must match RDKit API value");
        end

        function test_tanimoto_crossValidation_ethanolVsBenzene(tc)
        % Cross-validation for ethanol vs benzene.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("c1ccccc1");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            rdkitScore = emk.similarity.tanimoto(fp1, fp2);
            a = emk.fingerprint.toArray(fp1);
            b = emk.fingerprint.toArray(fp2);
            c = sum(a & b);
            manualScore = c / sum(a | b);
            tc.verifyEqual(rdkitScore, manualScore, "AbsTol", 1e-10, ...
                "Manual Tanimoto from toArray must match RDKit API for ethanol/benzene");
        end

        function test_tanimoto_crossValidation_maccs_ethanolVsAspirin(tc)
        % Cross-validation for MACCS fingerprints.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1   = emk.fingerprint.maccs(mol1);
            fp2   = emk.fingerprint.maccs(mol2);
            rdkitScore = emk.similarity.tanimoto(fp1, fp2);
            a = emk.fingerprint.toArray(fp1);
            b = emk.fingerprint.toArray(fp2);
            c = sum(a & b);
            manualScore = c / sum(a | b);
            tc.verifyEqual(rdkitScore, manualScore, "AbsTol", 1e-10, ...
                "MACCS Tanimoto from toArray must match RDKit API value");
        end

        % ------------------------------------------------------------------
        % TC10: MACCS fingerprints are accepted (requires RDKit)
        % ------------------------------------------------------------------

        function test_tanimoto_maccs_returnsDouble(tc)
        % tanimoto() must accept MACCS fingerprints and return double scalar.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("c1ccccc1");
            fp1   = emk.fingerprint.maccs(mol1);
            fp2   = emk.fingerprint.maccs(mol2);
            score = emk.similarity.tanimoto(fp1, fp2);
            tc.verifyClass(score, "double", ...
                "tanimoto with MACCS FPs must return double");
            tc.verifySize(score, [1, 1], ...
                "tanimoto with MACCS FPs must return scalar");
        end

        % ------------------------------------------------------------------
        % TC11: Mol object passed as fp => rdkitError (requires RDKit)
        % A Mol starts with "py." so it passes fp1/fp2 validation,
        % but RDKit's TanimotoSimilarity will reject it.
        % ------------------------------------------------------------------

        function test_tanimoto_molAsFp1_throwsRdkitError(tc)
        % Passing a Mol object as fp1 must throw rdkitError.
        % Mol starts with "py." so it passes the Python-object check,
        % but RDKit's TanimotoSimilarity rejects it.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.similarity.tanimoto(mol, fp), ...
                "emk:similarity:tanimoto:rdkitError", ...
                "Passing Mol as fp1 must throw rdkitError");
        end

        % ------------------------------------------------------------------
        % TC12: Bit-length mismatch => rdkitError (requires RDKit)
        % Morgan default = 2048 bits; MACCS = 167 bits.
        % Mixing fingerprint types must be detected and reported as rdkitError.
        % This guards against silent wrong results from cross-type comparisons.
        % ------------------------------------------------------------------

        function test_tanimoto_morgan2048_vs_maccs167_throwsRdkitError(tc)
        % Morgan (2048-bit) vs MACCS (167-bit) must throw rdkitError.
        % RDKit's TanimotoSimilarity raises when bit lengths differ.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fpM   = emk.fingerprint.morgan(mol);  % 2048-bit
            fpMAC = emk.fingerprint.maccs(mol);   % 167-bit
            tc.verifyError(@() emk.similarity.tanimoto(fpM, fpMAC), ...
                "emk:similarity:tanimoto:rdkitError", ...
                "Morgan 2048-bit vs MACCS 167-bit must throw rdkitError");
        end

        % ==================================================================
        % Dice tests (TC1d - TC11d)
        % ==================================================================

        % ------------------------------------------------------------------
        % TC1d: fp1 input validation -- no RDKit required
        % ------------------------------------------------------------------

        function test_dice_fp1_numericInput_throwsInvalidInput(tc)
        % Numeric fp1 must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.similarity.dice(42, 42), ...
                "emk:similarity:dice:invalidInput", ...
                "Numeric fp1 must throw invalidInput");
        end

        function test_dice_fp1_numericInput_errorMessage_containsClass(tc)
        % Error message must name the offending type for fp1.
            ME = tc.captureError(@() emk.similarity.dice(42, 42));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain 'double' for numeric fp1");
        end

        function test_dice_fp1_stringInput_throwsInvalidInput(tc)
        % String fp1 must throw invalidInput.
            tc.verifyError(@() emk.similarity.dice("CCO", "CCO"), ...
                "emk:similarity:dice:invalidInput", ...
                "String fp1 must throw invalidInput");
        end

        function test_dice_fp1_logicalInput_throwsInvalidInput(tc)
        % logical fp1 must throw invalidInput.
            tc.verifyError(@() emk.similarity.dice(true, true), ...
                "emk:similarity:dice:invalidInput", ...
                "logical fp1 must throw invalidInput");
        end

        function test_dice_fp1_charInput_throwsInvalidInput(tc)
        % char fp1 must throw invalidInput (no RDKit required).
            tc.verifyError(@() emk.similarity.dice('CCO', 'CCO'), ...
                "emk:similarity:dice:invalidInput", ...
                "char fp1 must throw invalidInput");
        end

        function test_dice_fp1_emptyInput_throwsInvalidInput(tc)
        % Empty matrix fp1 must throw invalidInput (no RDKit required).
            tc.verifyError(@() emk.similarity.dice([], []), ...
                "emk:similarity:dice:invalidInput", ...
                "Empty matrix fp1 must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC2d: fp2 input validation -- requires RDKit for valid fp1
        % ------------------------------------------------------------------

        function test_dice_fp2_numericInput_throwsInvalidInput(tc)
        % Valid fp1 + numeric fp2 must throw invalidInput.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp1 = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.similarity.dice(fp1, 42), ...
                "emk:similarity:dice:invalidInput", ...
                "Numeric fp2 with valid fp1 must throw invalidInput");
        end

        function test_dice_fp2_stringInput_throwsInvalidInput(tc)
        % Valid fp1 + string fp2 must throw invalidInput.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp1 = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.similarity.dice(fp1, "CCO"), ...
                "emk:similarity:dice:invalidInput", ...
                "String fp2 must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC3d: Return type is double scalar (requires RDKit)
        % ------------------------------------------------------------------

        function test_dice_returnsDoubleScalar(tc)
        % dice() must return a double scalar.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.morgan(mol);
            score = emk.similarity.dice(fp, fp);
            tc.verifyClass(score, "double", ...
                "dice must return a double");
            tc.verifySize(score, [1, 1], ...
                "dice must return a scalar");
        end

        function test_dice_returnIsFinite(tc)
        % Return value must be finite (not NaN or Inf).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.morgan(mol);
            score = emk.similarity.dice(fp, fp);
            tc.verifyTrue(isfinite(score), ...
                "dice must return a finite value");
        end

        function test_dice_returnIsReal(tc)
        % Return value must be real-valued (not complex).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.morgan(mol);
            score = emk.similarity.dice(fp, fp);
            tc.verifyTrue(isreal(score), ...
                "dice must return a real-valued (non-complex) double");
        end

        % ------------------------------------------------------------------
        % TC4d: Identical fingerprints => score == 1.0 (requires RDKit)
        % ------------------------------------------------------------------

        function test_dice_sameFP_morgan_equals1(tc)
        % Dice of a fingerprint with itself must be exactly 1.0.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.morgan(mol);
            score = emk.similarity.dice(fp, fp);
            tc.verifyEqual(score, 1.0, ...
                "Dice self-similarity (Morgan) must be exactly 1.0");
        end

        function test_dice_sameFP_maccs_equals1(tc)
        % Dice of a MACCS fingerprint with itself must be exactly 1.0.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("c1ccccc1");
            fp    = emk.fingerprint.maccs(mol);
            score = emk.similarity.dice(fp, fp);
            tc.verifyEqual(score, 1.0, ...
                "Dice self-similarity (MACCS) must be exactly 1.0");
        end

        % ------------------------------------------------------------------
        % TC5d: Score in [0, 1] (requires RDKit)
        % ------------------------------------------------------------------

        function test_dice_ethanolVsAspirin_inRange(tc)
        % Dice score must be in [0, 1] for ethanol vs aspirin.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            score = emk.similarity.dice(fp1, fp2);
            tc.verifyGreaterThanOrEqual(score, 0.0, "Dice score must be >= 0");
            tc.verifyLessThanOrEqual(score, 1.0, "Dice score must be <= 1");
        end

        function test_dice_ethanolVsBenzene_inRange(tc)
        % Dice score must be in [0, 1] for ethanol vs benzene.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("c1ccccc1");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            score = emk.similarity.dice(fp1, fp2);
            tc.verifyGreaterThanOrEqual(score, 0.0, "Dice score must be >= 0");
            tc.verifyLessThanOrEqual(score, 1.0, "Dice score must be <= 1");
        end

        % ------------------------------------------------------------------
        % TC6d: Symmetry: dice(fp1,fp2) == dice(fp2,fp1) (requires RDKit)
        % ------------------------------------------------------------------

        function test_dice_symmetry_ethanolVsAspirin(tc)
        % dice(fp1, fp2) must equal dice(fp2, fp1).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1   = emk.mol.fromSmiles("CCO");
            mol2   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1    = emk.fingerprint.morgan(mol1);
            fp2    = emk.fingerprint.morgan(mol2);
            score1 = emk.similarity.dice(fp1, fp2);
            score2 = emk.similarity.dice(fp2, fp1);
            tc.verifyEqual(score1, score2, ...
                "Dice must be symmetric: D(A,B) == D(B,A)");
        end

        % ------------------------------------------------------------------
        % TC7d: Different molecules => score < 1.0 (requires RDKit)
        % ------------------------------------------------------------------

        function test_dice_ethanolVsAspirin_lessThan1(tc)
        % Structurally different molecules must have Dice score < 1.0.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            score = emk.similarity.dice(fp1, fp2);
            tc.verifyLessThan(score, 1.0, ...
                "Ethanol vs aspirin Dice must be less than 1.0");
        end

        % ------------------------------------------------------------------
        % TC8d: Cross-validation D = 2c/(a+b) via toArray (requires RDKit)
        % Primary correctness test: MATLAB computation vs RDKit API.
        % ------------------------------------------------------------------

        function test_dice_crossValidation_ethanolVsAspirin(tc)
        % Manual D = 2c/(a+b) from toArray must match RDKit DiceSimilarity.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            rdkitScore  = emk.similarity.dice(fp1, fp2);
            a = emk.fingerprint.toArray(fp1);
            b = emk.fingerprint.toArray(fp2);
            manualScore = 2 * sum(a & b) / (sum(a) + sum(b));
            tc.verifyEqual(rdkitScore, manualScore, "AbsTol", 1e-10, ...
                "Manual Dice from toArray must match RDKit DiceSimilarity");
        end

        function test_dice_crossValidation_ethanolVsBenzene(tc)
        % Cross-validation for ethanol vs benzene (aliphatic vs aromatic).
        % Tests a different structural class pair to ensure generality.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("c1ccccc1");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            rdkitScore  = emk.similarity.dice(fp1, fp2);
            a = emk.fingerprint.toArray(fp1);
            b = emk.fingerprint.toArray(fp2);
            manualScore = 2 * sum(a & b) / (sum(a) + sum(b));
            tc.verifyEqual(rdkitScore, manualScore, "AbsTol", 1e-10, ...
                "Manual Dice from toArray must match RDKit DiceSimilarity (ethanol/benzene)");
        end

        function test_dice_crossValidation_maccs_ethanolVsAspirin(tc)
        % Cross-validation for MACCS fingerprints.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1   = emk.fingerprint.maccs(mol1);
            fp2   = emk.fingerprint.maccs(mol2);
            rdkitScore  = emk.similarity.dice(fp1, fp2);
            a = emk.fingerprint.toArray(fp1);
            b = emk.fingerprint.toArray(fp2);
            manualScore = 2 * sum(a & b) / (sum(a) + sum(b));
            tc.verifyEqual(rdkitScore, manualScore, "AbsTol", 1e-10, ...
                "MACCS Dice from toArray must match RDKit DiceSimilarity");
        end

        % ------------------------------------------------------------------
        % TC9d: Dice >= Tanimoto (requires RDKit)
        % For binary vectors, D(A,B) >= T(A,B) always holds.
        % This property cross-validates both metrics simultaneously.
        % ------------------------------------------------------------------

        function test_dice_greaterOrEqualTanimoto_ethanolVsAspirin(tc)
        % Dice >= Tanimoto for any pair of fingerprints.
        % Mathematical proof: D >= T iff 2c/(a+b) >= c/(a+b-c)
        % iff a+b >= 2c (since c <= min(a,b) <= (a+b)/2).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            dScore = emk.similarity.dice(fp1, fp2);
            tScore = emk.similarity.tanimoto(fp1, fp2);
            tc.verifyGreaterThanOrEqual(dScore, tScore - 1e-10, ...
                "Dice must be >= Tanimoto for any fingerprint pair");
        end

        function test_dice_greaterOrEqualTanimoto_ethanolVsBenzene(tc)
        % Dice >= Tanimoto also holds for a dissimilar pair (ethanol vs benzene).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1  = emk.mol.fromSmiles("CCO");
            mol2  = emk.mol.fromSmiles("c1ccccc1");
            fp1   = emk.fingerprint.morgan(mol1);
            fp2   = emk.fingerprint.morgan(mol2);
            dScore = emk.similarity.dice(fp1, fp2);
            tScore = emk.similarity.tanimoto(fp1, fp2);
            tc.verifyGreaterThanOrEqual(dScore, tScore - 1e-10, ...
                "Dice must be >= Tanimoto for ethanol vs benzene pair");
        end

        % ------------------------------------------------------------------
        % TC10d: Mol object as fp => rdkitError (requires RDKit)
        % ------------------------------------------------------------------

        function test_dice_molAsFp1_throwsRdkitError(tc)
        % Passing a Mol object as fp1 must throw rdkitError.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.similarity.dice(mol, fp), ...
                "emk:similarity:dice:rdkitError", ...
                "Passing Mol as fp1 to dice must throw rdkitError");
        end

        % ------------------------------------------------------------------
        % TC11d: Bit-length mismatch => rdkitError (requires RDKit)
        % ------------------------------------------------------------------

        function test_dice_morgan2048_vs_maccs167_throwsRdkitError(tc)
        % Morgan (2048-bit) vs MACCS (167-bit) must throw rdkitError.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fpM   = emk.fingerprint.morgan(mol);  % 2048-bit
            fpMAC = emk.fingerprint.maccs(mol);   % 167-bit
            tc.verifyError(@() emk.similarity.dice(fpM, fpMAC), ...
                "emk:similarity:dice:rdkitError", ...
                "Morgan 2048-bit vs MACCS 167-bit must throw rdkitError for dice");
        end

    end

    % ======================================================================
    % rankBy tests
    %
    % Coverage:
    %   RB1:  Non-Python queryFp => invalidQueryFp (no RDKit)
    %   RB2:  Non-cell dbFps => invalidDbFps (no RDKit, queryFp validated first)
    %         NOTE: queryFp must be Python for dbFps check to fire.
    %   RB3:  Empty cell dbFps => invalidDbFps (no RDKit)
    %   RB4:  Cell with non-Python element => invalidDbFps (no RDKit)
    %   RB5:  N = 0 => invalidN (no RDKit)
    %   RB6:  N = -1 => invalidN (no RDKit)
    %   RB7:  N = 1.5 (non-integer) => invalidN (no RDKit)
    %   RB8:  Invalid Metric => invalidMetric (no RDKit, checked first in arguments)
    %   RB9:  Returns struct with Indices, Scores, Metric fields (requires RDKit)
    %   RB10: Indices are double, Scores are double (requires RDKit)
    %   RB11: Scores are sorted descending (requires RDKit)
    %   RB12: Self-similarity: queryFp in dbFps => score=1.0 at top (requires RDKit)
    %   RB13: N=1 returns single result (requires RDKit)
    %   RB14: N > M returns all M results (requires RDKit)
    %   RB15: Default N (omitted) returns all M results (requires RDKit)
    %   RB16: Metric="dice" works and returns results in [0, 1] (requires RDKit)
    %   RB17: Bit-length mismatch => rdkitError (requires RDKit)
    %   RB18: result.Metric field matches requested Metric (requires RDKit)
    %   RB19: Cross-validation (tanimoto): top score matches tanimoto(queryFp, dbFps{idx}) (requires RDKit)
    %   RB20: N=matrix (non-scalar) => invalidN (requires RDKit to reach N check)
    %   RB21: Explicit N=Inf accepted, returns all M results (requires RDKit)
    %   RB22: MACCS fingerprints accepted (requires RDKit)
    %   RB23: All Scores[k] == tanimoto(query, db{Indices[k]}) for every k (requires RDKit)
    %   RB24: Cross-validation (dice): top score matches dice(queryFp, dbFps{idx}) (requires RDKit)
    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % RB1: queryFp validation -- no RDKit required
        % ------------------------------------------------------------------

        function test_rankBy_numericQueryFp_throwsInvalidQueryFp(tc)
        % Non-Python queryFp must throw invalidQueryFp before any RDKit call.
            tc.verifyError(@() emk.similarity.rankBy(42, {42}), ...
                "emk:similarity:rankBy:invalidQueryFp", ...
                "Numeric queryFp must throw invalidQueryFp");
        end

        function test_rankBy_stringQueryFp_throwsInvalidQueryFp(tc)
        % String queryFp must throw invalidQueryFp.
            tc.verifyError(@() emk.similarity.rankBy("CCO", {"CCO"}), ...
                "emk:similarity:rankBy:invalidQueryFp", ...
                "String queryFp must throw invalidQueryFp");
        end

        function test_rankBy_emptyQueryFp_throwsInvalidQueryFp(tc)
        % Empty matrix queryFp must throw invalidQueryFp.
            tc.verifyError(@() emk.similarity.rankBy([], {[]}), ...
                "emk:similarity:rankBy:invalidQueryFp", ...
                "Empty queryFp must throw invalidQueryFp");
        end

        % ------------------------------------------------------------------
        % RB2-RB4: dbFps validation -- requires Python queryFp
        % ------------------------------------------------------------------

        function test_rankBy_nonCellDbFps_throwsInvalidDbFps(tc)
        % When queryFp is valid Python but dbFps is not a cell, must throw.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            tc.verifyError(@() emk.similarity.rankBy(fp, 42), ...
                "emk:similarity:rankBy:invalidDbFps", ...
                "Non-cell dbFps must throw invalidDbFps");
        end

        function test_rankBy_emptyCellDbFps_throwsInvalidDbFps(tc)
        % Empty cell dbFps must throw invalidDbFps.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            tc.verifyError(@() emk.similarity.rankBy(fp, {}), ...
                "emk:similarity:rankBy:invalidDbFps", ...
                "Empty cell dbFps must throw invalidDbFps");
        end

        function test_rankBy_cellWithNonPythonElement_throwsInvalidDbFps(tc)
        % Cell containing a non-Python element must throw invalidDbFps.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            tc.verifyError(@() emk.similarity.rankBy(fp, {fp, 42}), ...
                "emk:similarity:rankBy:invalidDbFps", ...
                "Cell with non-Python element must throw invalidDbFps");
        end

        % ------------------------------------------------------------------
        % RB5-RB7: N validation -- no RDKit required (queryFp/dbFps validated first)
        % Use double as proxy inputs to trigger N check (cannot reach N check
        % with non-Python queryFp; must have valid inputs up to N).
        % NOTE: Since arguments block parses N before code runs, test with
        % valid Python objects to reach the N validation path.
        % ------------------------------------------------------------------

        function test_rankBy_nZero_throwsInvalidN(tc)
        % N = 0 must throw invalidN.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            tc.verifyError(@() emk.similarity.rankBy(fp, {fp}, 0), ...
                "emk:similarity:rankBy:invalidN", ...
                "N=0 must throw invalidN");
        end

        function test_rankBy_nNegative_throwsInvalidN(tc)
        % N = -1 must throw invalidN.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            tc.verifyError(@() emk.similarity.rankBy(fp, {fp}, -1), ...
                "emk:similarity:rankBy:invalidN", ...
                "Negative N must throw invalidN");
        end

        function test_rankBy_nNonInteger_throwsInvalidN(tc)
        % N = 1.5 (non-integer) must throw invalidN.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            tc.verifyError(@() emk.similarity.rankBy(fp, {fp}, 1.5), ...
                "emk:similarity:rankBy:invalidN", ...
                "Non-integer N must throw invalidN");
        end

        % ------------------------------------------------------------------
        % RB8: Metric validation -- fires before Python validation in arguments block
        % ------------------------------------------------------------------

        function test_rankBy_invalidMetric_throwsInvalidMetric(tc)
        % Unknown metric string must throw invalidMetric.
        % This fires at the start of the function before input validation.
            tc.verifyError(@() emk.similarity.rankBy(42, {42}, Inf, Metric="cosine"), ...
                "emk:similarity:rankBy:invalidMetric", ...
                "Unknown Metric must throw invalidMetric");
        end

        % ------------------------------------------------------------------
        % RB9-RB10: Return type and structure (requires RDKit)
        % ------------------------------------------------------------------

        function test_rankBy_returnsStruct(tc)
        % rankBy must return a struct.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp  = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            res = emk.similarity.rankBy(fp, {fp});
            tc.verifyClass(res, "struct", ...
                "rankBy must return a struct");
        end

        function test_rankBy_structHasRequiredFields(tc)
        % Result struct must contain Indices, Scores, and Metric.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp  = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            res = emk.similarity.rankBy(fp, {fp});
            tc.verifyTrue(isfield(res, "Indices"), "Missing field: Indices");
            tc.verifyTrue(isfield(res, "Scores"),  "Missing field: Scores");
            tc.verifyTrue(isfield(res, "Metric"),  "Missing field: Metric");
        end

        function test_rankBy_indicesAreDouble(tc)
        % Indices field must be a double array.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp  = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            res = emk.similarity.rankBy(fp, {fp});
            tc.verifyClass(res.Indices, "double", ...
                "Indices must be double");
        end

        function test_rankBy_scoresAreDouble(tc)
        % Scores field must be a double array.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp  = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            res = emk.similarity.rankBy(fp, {fp});
            tc.verifyClass(res.Scores, "double", ...
                "Scores must be double");
        end

        % ------------------------------------------------------------------
        % RB11: Scores are sorted descending (requires RDKit)
        % ------------------------------------------------------------------

        function test_rankBy_scoresAreSortedDescending(tc)
        % Scores in the result must be sorted in non-increasing order.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mol3 = emk.mol.fromSmiles("c1ccccc1");
            mol4 = emk.mol.fromSmiles("CCCO");
            fp1  = emk.fingerprint.morgan(mol1);
            db   = { emk.fingerprint.morgan(mol2), ...
                     emk.fingerprint.morgan(mol3), ...
                     emk.fingerprint.morgan(mol4) };
            res  = emk.similarity.rankBy(fp1, db);
            tc.verifyTrue(all(diff(res.Scores) <= 0), ...
                "Scores must be sorted in non-increasing (descending) order");
        end

        % ------------------------------------------------------------------
        % RB12: Self-similarity -- queryFp in dbFps gives score=1.0 at top
        % ------------------------------------------------------------------

        function test_rankBy_selfInDb_scoreIsOneAtTop(tc)
        % When queryFp is included in dbFps, the self-similarity score
        % must be 1.0 and appear as the first (highest-ranked) result.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1  = emk.fingerprint.morgan(mol1);
            fp2  = emk.fingerprint.morgan(mol2);
            % queryFp is fp1; dbFps = {fp2, fp1} -- fp1 is at index 2
            res = emk.similarity.rankBy(fp1, {fp2, fp1});
            tc.verifyEqual(res.Scores(1), 1.0, ...
                "Self-similarity must be 1.0 and ranked first");
            tc.verifyEqual(res.Indices(1), 2, ...
                "Self-similarity must be at the position of queryFp in dbFps");
        end

        % ------------------------------------------------------------------
        % RB13: N=1 returns single result (requires RDKit)
        % ------------------------------------------------------------------

        function test_rankBy_n1_returnsSingleResult(tc)
        % N=1 must return exactly one element in Indices and Scores.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mol3 = emk.mol.fromSmiles("c1ccccc1");
            fp1  = emk.fingerprint.morgan(mol1);
            db   = { emk.fingerprint.morgan(mol2), ...
                     emk.fingerprint.morgan(mol3) };
            res  = emk.similarity.rankBy(fp1, db, 1);
            tc.verifySize(res.Indices, [1, 1], "N=1 Indices must have size [1,1]");
            tc.verifySize(res.Scores,  [1, 1], "N=1 Scores must have size [1,1]");
        end

        % ------------------------------------------------------------------
        % RB14: N > M returns all M results (requires RDKit)
        % ------------------------------------------------------------------

        function test_rankBy_nGreaterThanM_returnsAllM(tc)
        % When N > M, all M results must be returned (no error).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1  = emk.fingerprint.morgan(mol1);
            db   = { emk.fingerprint.morgan(mol2) };  % M = 1
            res  = emk.similarity.rankBy(fp1, db, 100);  % N=100 > M=1
            tc.verifySize(res.Indices, [1, 1], "N>M must return M=1 results");
            tc.verifySize(res.Scores,  [1, 1], "N>M must return M=1 scores");
        end

        % ------------------------------------------------------------------
        % RB15: Default N (omitted) returns all M results (requires RDKit)
        % ------------------------------------------------------------------

        function test_rankBy_defaultN_returnsAllM(tc)
        % When N is omitted, all M fingerprints must be returned.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mol3 = emk.mol.fromSmiles("c1ccccc1");
            fp1  = emk.fingerprint.morgan(mol1);
            db   = { emk.fingerprint.morgan(mol2), emk.fingerprint.morgan(mol3) };
            res  = emk.similarity.rankBy(fp1, db);  % N omitted
            tc.verifySize(res.Indices, [1, 2], "Default N must return all M=2 results");
            tc.verifySize(res.Scores,  [1, 2], "Default N must return all M=2 scores");
        end

        % ------------------------------------------------------------------
        % RB16: Metric="dice" works and scores are in [0, 1] (requires RDKit)
        % ------------------------------------------------------------------

        function test_rankBy_metricDice_scoresInRange(tc)
        % Metric="dice" must return scores in [0, 1].
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1  = emk.fingerprint.morgan(mol1);
            db   = { emk.fingerprint.morgan(mol2) };
            res  = emk.similarity.rankBy(fp1, db, Inf, Metric="dice");
            tc.verifyGreaterThanOrEqual(res.Scores(1), 0.0, ...
                "Dice score must be >= 0");
            tc.verifyLessThanOrEqual(res.Scores(1), 1.0, ...
                "Dice score must be <= 1");
        end

        % ------------------------------------------------------------------
        % RB17: Bit-length mismatch => rdkitError (requires RDKit)
        % ------------------------------------------------------------------

        function test_rankBy_bitLengthMismatch_throwsRdkitError(tc)
        % Mixing Morgan (2048-bit) query with MACCS (167-bit) db must throw.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fpM   = emk.fingerprint.morgan(mol);  % 2048-bit
            fpMAC = emk.fingerprint.maccs(mol);   % 167-bit
            tc.verifyError(@() emk.similarity.rankBy(fpM, {fpMAC}), ...
                "emk:similarity:rankBy:rdkitError", ...
                "Bit-length mismatch must throw rdkitError");
        end

        % ------------------------------------------------------------------
        % RB18: result.Metric field matches requested Metric (requires RDKit)
        % ------------------------------------------------------------------

        function test_rankBy_metricField_tanimoto(tc)
        % result.Metric must equal "tanimoto" when Metric="tanimoto".
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp  = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            res = emk.similarity.rankBy(fp, {fp}, Inf, Metric="tanimoto");
            tc.verifyEqual(res.Metric, "tanimoto", ...
                "result.Metric must be ""tanimoto""");
        end

        function test_rankBy_metricField_dice(tc)
        % result.Metric must equal "dice" when Metric="dice".
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp  = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            res = emk.similarity.rankBy(fp, {fp}, Inf, Metric="dice");
            tc.verifyEqual(res.Metric, "dice", ...
                "result.Metric must be ""dice""");
        end

        % ------------------------------------------------------------------
        % RB19: Cross-validation vs tanimoto() (requires RDKit)
        %
        % The top score from rankBy(Metric="tanimoto") must exactly match
        % emk.similarity.tanimoto(queryFp, dbFps{topIdx}).
        % This ensures the Bulk API and the single-pair API agree.
        % ------------------------------------------------------------------

        function test_rankBy_crossValidation_topScoreMatchesTanimoto(tc)
        % Top score from rankBy must match emk.similarity.tanimoto pairwise call.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mol3 = emk.mol.fromSmiles("CCCO");  % n-propanol, structurally close to ethanol
            fp1  = emk.fingerprint.morgan(mol1);
            db   = { emk.fingerprint.morgan(mol2), emk.fingerprint.morgan(mol3) };
            res  = emk.similarity.rankBy(fp1, db);
            topIdx   = res.Indices(1);
            topScore = res.Scores(1);
            pairScore = emk.similarity.tanimoto(fp1, db{topIdx});
            tc.verifyEqual(topScore, pairScore, "AbsTol", 1e-10, ...
                "Top score from rankBy must match tanimoto() pairwise call");
        end

        % ------------------------------------------------------------------
        % RB20: N=matrix (non-scalar) => invalidN (requires RDKit to reach N)
        % N validation fires after queryFp and dbFps checks (both require Python).
        % ------------------------------------------------------------------

        function test_rankBy_nMatrix_throwsInvalidN(tc)
        % Non-scalar N=[1,2] must throw invalidN (fails ~isscalar check).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            tc.verifyError(@() emk.similarity.rankBy(fp, {fp}, [1, 2]), ...
                "emk:similarity:rankBy:invalidN", ...
                "Non-scalar N=[1,2] must throw invalidN");
        end

        % ------------------------------------------------------------------
        % RB21: Explicit N=Inf accepted, returns all M (requires RDKit)
        %
        % Documents that Inf is the canonical "return everything" sentinel
        % and must not be rejected by the N validator.
        % ------------------------------------------------------------------

        function test_rankBy_nInfExplicit_returnsAllM(tc)
        % Passing N=Inf explicitly must not throw and must return all M results.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            fp1  = emk.fingerprint.morgan(mol1);
            db   = { emk.fingerprint.morgan(mol2) };   % M = 1
            res  = emk.similarity.rankBy(fp1, db, Inf);
            tc.verifySize(res.Indices, [1, 1], ...
                "Explicit N=Inf must return all M=1 results");
        end

        % ------------------------------------------------------------------
        % RB22: MACCS fingerprints accepted (requires RDKit)
        %
        % The function spec accepts any ExplicitBitVect; MACCS (167-bit)
        % must work without error as long as all FPs share the same bit length.
        % ------------------------------------------------------------------

        function test_rankBy_maccs_returnsCorrectSize(tc)
        % rankBy with MACCS (167-bit) FPs must return struct of correct size.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mol3 = emk.mol.fromSmiles("CCCO");
            fp1  = emk.fingerprint.maccs(mol1);
            db   = { emk.fingerprint.maccs(mol2), ...
                     emk.fingerprint.maccs(mol3) };
            res  = emk.similarity.rankBy(fp1, db);
            tc.verifyClass(res, "struct", ...
                "rankBy with MACCS FPs must return struct");
            tc.verifySize(res.Indices, [1, 2], ...
                "rankBy with MACCS FPs must return M=2 results");
            tc.verifyClass(res.Scores, "double", ...
                "rankBy with MACCS FPs must return double Scores");
        end

        % ------------------------------------------------------------------
        % RB23: Full index-score correspondence for ALL results (requires RDKit)
        %
        % RB19 only verifies the top-1 result.  This test verifies that
        % for every k, Scores(k) == tanimoto(queryFp, dbFps{Indices(k)}).
        % Any off-by-one in the sort or index mapping will be caught here.
        % ------------------------------------------------------------------

        function test_rankBy_allIndices_correspondToCorrectScores(tc)
        % Every Scores(k) must equal tanimoto(queryFp, dbFps{Indices(k)}).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mol3 = emk.mol.fromSmiles("CCCO");
            mol4 = emk.mol.fromSmiles("c1ccccc1");
            fp1  = emk.fingerprint.morgan(mol1);
            db   = { emk.fingerprint.morgan(mol2), ...
                     emk.fingerprint.morgan(mol3), ...
                     emk.fingerprint.morgan(mol4) };
            res  = emk.similarity.rankBy(fp1, db);
            for k = 1:numel(res.Indices)
                expected = emk.similarity.tanimoto(fp1, db{res.Indices(k)});
                tc.verifyEqual(res.Scores(k), expected, "AbsTol", 1e-10, ...
                    sprintf("Scores(%d) must equal tanimoto(query, db{Indices(%d)})", k, k));
            end
        end

        % ------------------------------------------------------------------
        % RB24: Dice cross-validation (requires RDKit)
        %
        % Analogous to RB19 but for Metric="dice".
        % Top dice score from rankBy must match emk.similarity.dice() pairwise.
        % Ensures BulkDiceSimilarity agrees with DiceSimilarity.
        % ------------------------------------------------------------------

        function test_rankBy_diceMetric_crossValidation_topScoreMatchesDice(tc)
        % Top dice score from rankBy must match emk.similarity.dice() pairwise call.
        % AbsTol=1e-10: only floating-point rounding expected between BulkDice and Dice.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mol3 = emk.mol.fromSmiles("CCCO");   % closest to ethanol
            fp1  = emk.fingerprint.morgan(mol1);
            db   = { emk.fingerprint.morgan(mol2), emk.fingerprint.morgan(mol3) };
            res  = emk.similarity.rankBy(fp1, db, Inf, Metric="dice");
            topIdx    = res.Indices(1);
            topScore  = res.Scores(1);
            pairScore = emk.similarity.dice(fp1, db{topIdx});
            tc.verifyEqual(topScore, pairScore, "AbsTol", 1e-10, ...
                "Top dice score from rankBy must match dice() pairwise call");
        end

    end

    % ======================================================================
    % matrix tests
    %
    % Coverage:
    %   M1:  Non-cell fps => invalidInput (no RDKit required)
    %   M2:  Empty cell fps => invalidInput (no RDKit required)
    %   M3:  Cell with non-Python element => invalidInput (requires RDKit
    %         for fps{1}; fps{2}=42 triggers the check)
    %   M4:  Invalid Metric => invalidMetric (no RDKit, checked first)
    %   M5:  Returns NxN double matrix (requires RDKit)
    %   M6:  Returns symmetric matrix S == S' (requires RDKit)
    %   M7:  Diagonal entries == 1.0 (requires RDKit)
    %   M8:  Off-diagonal scores in [0, 1] (requires RDKit)
    %   M9:  Cross-validation: S(i,j) == tanimoto(fps{i}, fps{j}) (requires RDKit)
    %        Primary correctness test -- guards against Bulk API regression.
    %   M10: Metric="dice" produces valid NxN result (requires RDKit)
    %   M11: N=3 cell of Python objects does not throw (requires RDKit)
    %        Regression test for G13: py.list(cell) with Python objects
    %        must use the element-wise append loop, not py.list(cell) directly.
    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % M1: Non-cell fps => invalidInput (no RDKit required)
        % ------------------------------------------------------------------

        function test_matrix_numericInput_throwsInvalidInput(tc)
        % Numeric fps must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.similarity.matrix(42), ...
                "emk:similarity:matrix:invalidInput", ...
                "Numeric fps must throw invalidInput");
        end

        function test_matrix_stringInput_throwsInvalidInput(tc)
        % String fps must throw invalidInput.
            tc.verifyError(@() emk.similarity.matrix("CCO"), ...
                "emk:similarity:matrix:invalidInput", ...
                "String fps must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % M2: Empty cell fps => invalidInput (no RDKit required)
        % ------------------------------------------------------------------

        function test_matrix_emptyCellInput_throwsInvalidInput(tc)
        % Empty cell array must throw invalidInput.
            tc.verifyError(@() emk.similarity.matrix({}), ...
                "emk:similarity:matrix:invalidInput", ...
                "Empty cell must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % M3: Cell with non-Python element => invalidInput (requires RDKit)
        % ------------------------------------------------------------------

        function test_matrix_cellWithNonPythonElement_throwsInvalidInput(tc)
        % A cell containing a non-Python element must throw invalidInput.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            fp = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
            tc.verifyError(@() emk.similarity.matrix({fp, 42}), ...
                "emk:similarity:matrix:invalidInput", ...
                "Cell with non-Python element must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % M4: Invalid Metric => invalidMetric (no RDKit required)
        % ------------------------------------------------------------------

        function test_matrix_invalidMetric_throwsInvalidMetric(tc)
        % Unknown Metric string must throw invalidMetric.
            tc.verifyError(@() emk.similarity.matrix({42}, Metric="cosine"), ...
                "emk:similarity:matrix:invalidMetric", ...
                "Unknown Metric must throw invalidMetric");
        end

        % ------------------------------------------------------------------
        % M5: Returns NxN double matrix (requires RDKit)
        % ------------------------------------------------------------------

        function test_matrix_returnsSizeNxN(tc)
        % matrix() must return an N x N double matrix for N fingerprints.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smis = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
            fps  = cellfun(@(s) emk.fingerprint.morgan(emk.mol.fromSmiles(s)), ...
                smis, "UniformOutput", false);
            S = emk.similarity.matrix(fps);
            tc.verifyClass(S, "double", "matrix must return double");
            tc.verifySize(S, [3, 3], "3-fp set must produce 3x3 matrix");
        end

        % ------------------------------------------------------------------
        % M6: Symmetric matrix S == S' (requires RDKit)
        % ------------------------------------------------------------------

        function test_matrix_isSymmetric(tc)
        % S must equal its transpose (symmetric within floating-point tolerance).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smis = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
            fps  = cellfun(@(s) emk.fingerprint.morgan(emk.mol.fromSmiles(s)), ...
                smis, "UniformOutput", false);
            S = emk.similarity.matrix(fps);
            tc.verifyEqual(S, S', "AbsTol", 1e-12, ...
                "Similarity matrix must be symmetric (S == S')");
        end

        % ------------------------------------------------------------------
        % M7: Diagonal == 1.0 (requires RDKit)
        % ------------------------------------------------------------------

        function test_matrix_diagonalIsOne(tc)
        % All diagonal entries must equal exactly 1.0 (self-similarity).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smis = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
            fps  = cellfun(@(s) emk.fingerprint.morgan(emk.mol.fromSmiles(s)), ...
                smis, "UniformOutput", false);
            S = emk.similarity.matrix(fps);
            tc.verifyEqual(diag(S), ones(3, 1), ...
                "Diagonal entries must all be 1.0 (self-similarity)");
        end

        % ------------------------------------------------------------------
        % M8: Off-diagonal scores in [0, 1] (requires RDKit)
        % ------------------------------------------------------------------

        function test_matrix_offDiagonalInRange(tc)
        % All off-diagonal entries must be in [0, 1].
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smis = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
            fps  = cellfun(@(s) emk.fingerprint.morgan(emk.mol.fromSmiles(s)), ...
                smis, "UniformOutput", false);
            S = emk.similarity.matrix(fps);
            offDiag = S(~eye(size(S), "logical"));
            tc.verifyGreaterThanOrEqual(min(offDiag), 0.0, ...
                "Off-diagonal entries must be >= 0");
            tc.verifyLessThanOrEqual(max(offDiag), 1.0, ...
                "Off-diagonal entries must be <= 1");
        end

        % ------------------------------------------------------------------
        % M9: Cross-validation S(i,j) == tanimoto(fps{i}, fps{j}) (requires RDKit)
        % Primary correctness test: each matrix entry must match the pairwise API.
        % ------------------------------------------------------------------

        function test_matrix_crossValidation_matchesPairwiseTanimoto(tc)
        % Every S(i,j) must exactly equal tanimoto(fps{i}, fps{j}).
        % Guards against BulkTanimotoSimilarity vs TanimotoSimilarity divergence.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smis = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
            fps  = cellfun(@(s) emk.fingerprint.morgan(emk.mol.fromSmiles(s)), ...
                smis, "UniformOutput", false);
            S = emk.similarity.matrix(fps);
            for i = 1:numel(fps)
                for j = 1:numel(fps)
                    expected = emk.similarity.tanimoto(fps{i}, fps{j});
                    tc.verifyEqual(S(i,j), expected, "AbsTol", 1e-10, ...
                        sprintf("S(%d,%d) must match pairwise tanimoto", i, j));
                end
            end
        end

        % ------------------------------------------------------------------
        % M10: Metric="dice" produces valid NxN result (requires RDKit)
        % ------------------------------------------------------------------

        function test_matrix_metricDice_returnsValidMatrix(tc)
        % Metric="dice" must produce an NxN symmetric matrix with diagonal=1.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smis = {"CCO", "c1ccccc1"};
            fps  = cellfun(@(s) emk.fingerprint.morgan(emk.mol.fromSmiles(s)), ...
                smis, "UniformOutput", false);
            S = emk.similarity.matrix(fps, Metric="dice");
            tc.verifyClass(S, "double", "Dice matrix must be double");
            tc.verifySize(S, [2, 2], "2-fp dice matrix must be 2x2");
            tc.verifyEqual(diag(S), ones(2, 1), ...
                "Dice diagonal must be 1.0");
            tc.verifyEqual(S, S', "AbsTol", 1e-12, ...
                "Dice matrix must be symmetric");
        end

        % ------------------------------------------------------------------
        % M11: N=3 cell of Python objects does not throw (requires RDKit)
        % Regression test for G13: py.list(cell) fails when the cell contains
        % Python objects (ExplicitBitVect).  The fix uses element-wise append.
        % ------------------------------------------------------------------

        function test_matrix_threeFps_noError(tc)
        % A 3-element cell of Python fps must complete without error.
        % Regression: was crashing with "MATLAB cell conversion to Python
        % only supported for 1-N vectors" before G13 fix.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smis = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
            fps  = cellfun(@(s) emk.fingerprint.morgan(emk.mol.fromSmiles(s)), ...
                smis, "UniformOutput", false);
            S = [];
            tc.verifyWarningFree(@() assignin("caller", "S", ...
                emk.similarity.matrix(fps)), ...
                "3-fp matrix must complete without warning or error");
            S = emk.similarity.matrix(fps);
            tc.verifySize(S, [3, 3], "Result must be 3x3");
        end

    end

    % ======================================================================
    methods (Access = private)

        function tf = rdkitAvailable(~)
        % Return true if Python is configured and RDKit can be imported.
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
