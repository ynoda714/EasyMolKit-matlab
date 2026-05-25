classdef TestIo < matlab.unittest.TestCase
% TestIo  Unit tests for emk.io (readSdf, writeSdf, readSmilesList).
%
% Run with:
%   addpath(genpath("src"));
%   suite = testsuite("tests/unit");
%   runner = matlab.unittest.TestRunner.withNoPlugins;
%   results = runner.run(suite);
%
% ======================================================================
% Coverage (readSdf):
%
%   TC1:  Non-string/char filePath => invalidInput (no RDKit required)
%         + double, logical, cell rejected
%   TC2:  Non-existent file => fileNotFound (no RDKit required)
%         + error message contains path
%   TC3:  char filePath accepted (requires RDKit)
%   TC4:  Valid SDF => cell array; each element is a Python object (requires RDKit)
%   TC5:  Result shape is 1xN row cell array (requires RDKit)
%   TC6:  Molecule count matches write count (requires RDKit)
%         Round-trip write+read: numel(mols) preserved
%   TC7:  SMILES round-trip: canonical SMILES preserved (requires RDKit)
%         toSmiles(readSdf{i}) == toSmiles(original mols{i})
%
% Coverage (writeSdf):
%
%   TC8:  Non-cell mols => invalidInput (no RDKit required)
%         + double, string rejected
%   TC9:  Empty cell mols => invalidInput (no RDKit required)
%   TC10: Cell element not a Python object => invalidMol (no RDKit required)
%         + error message contains element index (index 1 and index 2)
%   TC11: Non-string/char filePath => invalidInput (no RDKit required)
%         (filePath check fires before mols-content check)
%   TC12: Non-existent parent directory => dirNotFound (requires RDKit)
%   TC13: char filePath accepted (requires RDKit)
%   TC14: Single molecule => file created and non-empty (requires RDKit)
%   TC15: Multiple molecules => file created (requires RDKit)
%
% Coverage (readSmilesList):
%
%   TC16: Non-string/char filePath => invalidInput (no RDKit required)
%   TC17: Non-existent file => fileNotFound (no RDKit required)
%         + error message contains path
%   TC18: Empty file => empty cell {} (no RDKit required)
%   TC19: Comment-only file => empty cell {} (no RDKit required)
%   TC20: char filePath accepted (requires RDKit)
%   TC21: Valid SMILES => cell array of Python objects (requires RDKit)
%   TC22: Result shape is 1xN row cell array (requires RDKit)
%   TC23: Molecule count matches SMILES count (requires RDKit)
%   TC24: SMILES cross-validation: canonical SMILES preserved (requires RDKit)
%   TC25: Name column (tab-separated) is ignored (requires RDKit)
%   TC26: Comment lines are skipped (requires RDKit)
%   TC27: Invalid SMILES are skipped without throwing (requires RDKit)
%   TC28: Mixed valid/invalid SMILES (requires RDKit)
%
% Design note on round-trip (TC5, TC6):
%   writeSdf + readSdf round-trip is the primary integration test.
%   Canonical SMILES comparison guards against atom reordering / sanitization
%   side effects.  RDKit may re-order atoms during SDF read, so we compare
%   via emk.mol.toSmiles() which returns canonical SMILES.
%
% ======================================================================

    properties
        tmpDir  string  % per-test temporary directory
    end

    methods (TestMethodSetup)
        function setupPath(tc)
            addpath(genpath("src"));
        end

        function createTmpDir(tc)
        % Create a unique temp directory for each test.
            tc.tmpDir = string(tempname());
            mkdir(tc.tmpDir);
        end
    end

    methods (TestMethodTeardown)
        function removeTmpDir(tc)
        % Remove temp directory after each test.
            if isfolder(tc.tmpDir)
                rmdir(tc.tmpDir, "s");
            end
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: readSdf filePath input validation -- no RDKit required
        % ------------------------------------------------------------------

        function test_readSdf_numericPath_throwsInvalidInput(tc)
        % Numeric filePath must throw invalidInput.
            tc.verifyError(@() emk.io.readSdf(42), ...
                "emk:io:readSdf:invalidInput", ...
                "Numeric filePath must throw invalidInput");
        end

        function test_readSdf_logicalPath_throwsInvalidInput(tc)
        % Logical filePath must throw invalidInput.
            tc.verifyError(@() emk.io.readSdf(true), ...
                "emk:io:readSdf:invalidInput", ...
                "Logical filePath must throw invalidInput");
        end

        function test_readSdf_cellPath_throwsInvalidInput(tc)
        % Cell filePath must throw invalidInput.
            tc.verifyError(@() emk.io.readSdf({"path.sdf"}), ...
                "emk:io:readSdf:invalidInput", ...
                "Cell filePath must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC2: readSdf non-existent file -- no RDKit required
        % ------------------------------------------------------------------

        function test_readSdf_nonExistentFile_throwsFileNotFound(tc)
        % Non-existent file path must throw fileNotFound.
            tc.verifyError(@() emk.io.readSdf("nonexistent_file_12345.sdf"), ...
                "emk:io:readSdf:fileNotFound", ...
                "Non-existent file must throw fileNotFound");
        end

        function test_readSdf_fileNotFound_errorMessage_containsPath(tc)
        % Error message must contain the offending file path.
            ME = tc.captureError(@() emk.io.readSdf("no_such_file.sdf"));
            tc.assertNotEmpty(ME, "Expected fileNotFound to be thrown");
            tc.verifySubstring(ME.message, "no_such_file.sdf", ...
                "Error message must contain the missing file path");
        end

        % ------------------------------------------------------------------
        % TC3: readSdf accepts char filePath (requires RDKit)
        % ------------------------------------------------------------------

        function test_readSdf_charPath_accepted(tc)
        % readSdf must accept a char (not just string) file path.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            sdfPath = tc.writeTestSdf({"CCO"});
            % Pass as char array -- the function must not error
            result = emk.io.readSdf(char(sdfPath));
            tc.verifyEqual(numel(result), 1, ...
                "readSdf with char path must return 1 molecule");
        end

        % ------------------------------------------------------------------
        % TC4-5: readSdf returns correct type and shape (requires RDKit)
        % ------------------------------------------------------------------

        function test_readSdf_validSdf_returnsCell(tc)
        % readSdf must return a cell array.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            sdfPath = tc.writeTestSdf({"CCO"});
            result = emk.io.readSdf(sdfPath);
            tc.verifyClass(result, "cell", ...
                "readSdf must return a cell array");
        end

        function test_readSdf_result_isRowCell(tc)
        % readSdf result must be a 1×N row cell (not N×1 column).
        % This is required for consistent iteration with {result{:}} and for(i=1:numel).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            sdfPath = tc.writeTestSdf({"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"});
            result = emk.io.readSdf(sdfPath);
            tc.verifyEqual(size(result, 1), 1, ...
                "readSdf result must have 1 row (row cell array)");
            tc.verifyEqual(size(result, 2), 3, ...
                "readSdf result must have N columns");
        end

        function test_readSdf_validSdf_eachElementIsPythonObject(tc)
        % Each element of the returned cell must be a Python Mol object.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            sdfPath = tc.writeTestSdf({"CCO", "c1ccccc1"});
            result = emk.io.readSdf(sdfPath);
            tc.verifyGreaterThan(numel(result), 0, ...
                "readSdf must return at least one molecule");
            for i = 1:numel(result)
                tc.verifyTrue(startsWith(class(result{i}), "py."), ...
                    sprintf("mols{%d} must be a Python object, got: %s", ...
                    i, class(result{i})));
            end
        end

        % ------------------------------------------------------------------
        % TC5: Round-trip molecule count (requires RDKit)
        % ------------------------------------------------------------------

        function test_readSdf_roundTrip_countMatches_single(tc)
        % Round-trip: writeSdf + readSdf => same count (1 molecule).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            sdfPath = tc.writeTestSdf({"CCO"});
            result = emk.io.readSdf(sdfPath);
            tc.verifyEqual(numel(result), 1, ...
                "Round-trip must preserve molecule count (1)");
        end

        function test_readSdf_roundTrip_countMatches_multiple(tc)
        % Round-trip: writeSdf + readSdf => same count (3 molecules).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smilesList = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
            sdfPath = tc.writeTestSdf(smilesList);
            result = emk.io.readSdf(sdfPath);
            tc.verifyEqual(numel(result), numel(smilesList), ...
                "Round-trip must preserve molecule count (3)");
        end

        % ------------------------------------------------------------------
        % TC6: Round-trip SMILES cross-validation (requires RDKit)
        % Canonical SMILES from readSdf must match original SMILES.
        % This is the primary correctness test for both readSdf and writeSdf.
        % ------------------------------------------------------------------

        function test_readSdf_roundTrip_smilesPreserved_ethanol(tc)
        % Round-trip: canonical SMILES of ethanol is preserved after write+read.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            sdfPath = tc.writeTestSdf({"CCO"});
            result = emk.io.readSdf(sdfPath);
            tc.verifyEqual(numel(result), 1, "Expected 1 molecule");
            origMol   = emk.mol.fromSmiles("CCO");
            origSmiles = emk.mol.toSmiles(origMol);
            rtSmiles   = emk.mol.toSmiles(result{1});
            tc.verifyEqual(rtSmiles, origSmiles, ...
                "Round-trip canonical SMILES must match original (ethanol)");
        end

        function test_readSdf_roundTrip_smilesPreserved_aspirin(tc)
        % Round-trip: canonical SMILES of aspirin is preserved after write+read.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            sdfPath = tc.writeTestSdf({"CC(=O)Oc1ccccc1C(=O)O"});
            result = emk.io.readSdf(sdfPath);
            tc.verifyEqual(numel(result), 1, "Expected 1 molecule");
            origMol    = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            origSmiles  = emk.mol.toSmiles(origMol);
            rtSmiles    = emk.mol.toSmiles(result{1});
            tc.verifyEqual(rtSmiles, origSmiles, ...
                "Round-trip canonical SMILES must match original (aspirin)");
        end

        function test_readSdf_roundTrip_smilesPreserved_allThree(tc)
        % Round-trip for 3 molecules: all canonical SMILES preserved.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            origSmiles = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
            sdfPath = tc.writeTestSdf(origSmiles);
            result = emk.io.readSdf(sdfPath);
            tc.verifyEqual(numel(result), 3, "Expected 3 molecules");
            for i = 1:3
                orig = emk.mol.toSmiles(emk.mol.fromSmiles(origSmiles{i}));
                rt   = emk.mol.toSmiles(result{i});
                tc.verifyEqual(rt, orig, ...
                    sprintf("Round-trip SMILES must match for molecule %d", i));
            end
        end

        % ------------------------------------------------------------------
        % TC7: writeSdf mols input validation -- no RDKit required
        % ------------------------------------------------------------------

        function test_writeSdf_numericMols_throwsInvalidInput(tc)
        % Numeric mols must throw invalidInput.
            sdfPath = fullfile(tc.tmpDir, "out.sdf");
            tc.verifyError(@() emk.io.writeSdf(42, sdfPath), ...
                "emk:io:writeSdf:invalidInput", ...
                "Numeric mols must throw invalidInput");
        end

        function test_writeSdf_stringMols_throwsInvalidInput(tc)
        % String mols must throw invalidInput.
            sdfPath = fullfile(tc.tmpDir, "out.sdf");
            tc.verifyError(@() emk.io.writeSdf("CCO", sdfPath), ...
                "emk:io:writeSdf:invalidInput", ...
                "String mols must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC8: writeSdf empty cell -- no RDKit required
        % ------------------------------------------------------------------

        function test_writeSdf_emptyCell_throwsInvalidInput(tc)
        % Empty cell mols must throw invalidInput.
            sdfPath = fullfile(tc.tmpDir, "out.sdf");
            tc.verifyError(@() emk.io.writeSdf({}, sdfPath), ...
                "emk:io:writeSdf:invalidInput", ...
                "Empty cell mols must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC9: writeSdf non-Python element -- no RDKit required
        % ------------------------------------------------------------------

        function test_writeSdf_nonPythonElement_throwsInvalidMol(tc)
        % Cell with non-Python element must throw invalidMol.
            sdfPath = fullfile(tc.tmpDir, "out.sdf");
            tc.verifyError(@() emk.io.writeSdf({42}, sdfPath), ...
                "emk:io:writeSdf:invalidMol", ...
                "Non-Python element must throw invalidMol");
        end

        function test_writeSdf_nonPythonElement_errorMessage_containsIndex(tc)
        % Error message must mention the offending element index (index=1).
            sdfPath = fullfile(tc.tmpDir, "out.sdf");
            ME = tc.captureError(@() emk.io.writeSdf({42}, sdfPath));
            tc.assertNotEmpty(ME, "Expected invalidMol to be thrown");
            tc.verifySubstring(ME.message, "1", ...
                "Error message must contain the element index '1'");
        end

        function test_writeSdf_secondElementBad_errorMessage_containsIndex2(tc)
        % When mols{2} is not a Python object, error message must say "2".
        % This verifies the per-element index is correctly reported for i>1.
        % mols{1} is a valid Python object (string starts with "py." not tested
        % here, but non-Python check fires for index 2).
        % Note: no RDKit needed; only index detection logic is tested.
        % We pass two non-Python elements; fp1 check fires at index 1, so
        % use a valid-looking non-Python placeholder at index 1 is not possible
        % without RDKit. Instead we pass {"ok_string", "bad"} and verify the
        % error fires at index 1 (first non-Python element found).
        % To specifically test index 2, we need a Python object at index 1.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol     = emk.mol.fromSmiles("CCO");      % valid Python mol
            sdfPath = fullfile(tc.tmpDir, "out.sdf");
            ME = tc.captureError(@() emk.io.writeSdf({mol, 42}, sdfPath));
            tc.assertNotEmpty(ME, "Expected invalidMol at index 2 to be thrown");
            tc.verifyEqual(string(ME.identifier), "emk:io:writeSdf:invalidMol", ...
                "Error must be invalidMol");
            tc.verifySubstring(ME.message, "2", ...
                "Error message must contain the element index '2'");
        end

        % ------------------------------------------------------------------
        % TC11: writeSdf filePath input validation -- no RDKit required
        % filePath is checked FIRST in the implementation (before mols content).
        % Passing numeric filePath fires invalidInput regardless of mols content.
        % ------------------------------------------------------------------

        function test_writeSdf_numericPath_throwsInvalidInput(tc)
        % Numeric filePath must throw invalidInput.
        % filePath validation fires before mols-content validation.
            tc.verifyError(@() emk.io.writeSdf({"mol"}, 42), ...
                "emk:io:writeSdf:invalidInput", ...
                "Numeric filePath must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC12: writeSdf non-existent parent directory (requires RDKit)
        % ------------------------------------------------------------------

        function test_writeSdf_nonExistentParentDir_throwsDirNotFound(tc)
        % Non-existent parent directory must throw dirNotFound.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            badPath = fullfile("nonexistent_dir_99999xyz", "out.sdf");
            tc.verifyError(@() emk.io.writeSdf({mol}, badPath), ...
                "emk:io:writeSdf:dirNotFound", ...
                "Non-existent parent dir must throw dirNotFound");
        end

        % ------------------------------------------------------------------
        % TC13: writeSdf accepts char filePath (requires RDKit)
        % ------------------------------------------------------------------

        function test_writeSdf_charPath_accepted(tc)
        % writeSdf must accept a char (not just string) file path.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol     = emk.mol.fromSmiles("CCO");
            sdfPath = char(fullfile(tc.tmpDir, "char_test.sdf"));
            emk.io.writeSdf({mol}, sdfPath);  % must not throw
            tc.verifyTrue(isfile(sdfPath), ...
                "writeSdf with char path must create the SDF file");
        end

        % ------------------------------------------------------------------
        % TC14: writeSdf single molecule (requires RDKit)
        % ------------------------------------------------------------------

        function test_writeSdf_singleMol_fileIsCreated(tc)
        % writeSdf with a single mol must create the output file.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol     = emk.mol.fromSmiles("CCO");
            sdfPath = fullfile(tc.tmpDir, "single.sdf");
            emk.io.writeSdf({mol}, sdfPath);
            tc.verifyTrue(isfile(sdfPath), ...
                "writeSdf must create the SDF file");
        end

        function test_writeSdf_singleMol_fileIsNonEmpty(tc)
        % Written SDF file must be non-empty (contains data).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol     = emk.mol.fromSmiles("CCO");
            sdfPath = fullfile(tc.tmpDir, "single.sdf");
            emk.io.writeSdf({mol}, sdfPath);
            info = dir(sdfPath);
            tc.verifyGreaterThan(info.bytes, 0, ...
                "Written SDF file must be non-empty");
        end

        % ------------------------------------------------------------------
        % TC15: writeSdf multiple molecules (requires RDKit)
        % ------------------------------------------------------------------

        function test_writeSdf_multipleMols_fileIsCreated(tc)
        % writeSdf with multiple mols must create the output file.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1    = emk.mol.fromSmiles("CCO");
            mol2    = emk.mol.fromSmiles("c1ccccc1");
            mol3    = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            sdfPath = fullfile(tc.tmpDir, "multi.sdf");
            emk.io.writeSdf({mol1, mol2, mol3}, sdfPath);
            tc.verifyTrue(isfile(sdfPath), ...
                "writeSdf (3 mols) must create the SDF file");
        end

        % ==================================================================
        % TC16-28: readSmilesList tests
        % ==================================================================

        % ------------------------------------------------------------------
        % TC16: readSmilesList filePath input validation -- no RDKit required
        % ------------------------------------------------------------------

        function test_readSmilesList_numericPath_throwsInvalidInput(tc)
        % Numeric filePath must throw invalidInput.
            tc.verifyError(@() emk.io.readSmilesList(42), ...
                "emk:io:readSmilesList:invalidInput", ...
                "Numeric filePath must throw invalidInput");
        end

        function test_readSmilesList_logicalPath_throwsInvalidInput(tc)
        % Logical filePath must throw invalidInput.
            tc.verifyError(@() emk.io.readSmilesList(true), ...
                "emk:io:readSmilesList:invalidInput", ...
                "Logical filePath must throw invalidInput");
        end

        function test_readSmilesList_cellPath_throwsInvalidInput(tc)
        % Cell filePath must throw invalidInput.
            tc.verifyError(@() emk.io.readSmilesList({"path.txt"}), ...
                "emk:io:readSmilesList:invalidInput", ...
                "Cell filePath must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC17: readSmilesList non-existent file -- no RDKit required
        % ------------------------------------------------------------------

        function test_readSmilesList_nonExistentFile_throwsFileNotFound(tc)
        % Non-existent file path must throw fileNotFound.
            tc.verifyError(@() emk.io.readSmilesList("no_such_smiles_file.txt"), ...
                "emk:io:readSmilesList:fileNotFound", ...
                "Non-existent file must throw fileNotFound");
        end

        function test_readSmilesList_fileNotFound_errorMessage_containsPath(tc)
        % Error message must contain the offending file path.
            ME = tc.captureError(@() emk.io.readSmilesList("missing_file.txt"));
            tc.assertNotEmpty(ME, "Expected fileNotFound to be thrown");
            tc.verifySubstring(ME.message, "missing_file.txt", ...
                "Error message must contain the missing file path");
        end

        % ------------------------------------------------------------------
        % TC18: readSmilesList empty file -- no RDKit required
        % ------------------------------------------------------------------

        function test_readSmilesList_emptyFile_returnsEmptyCell(tc)
        % An empty file must return an empty cell {}.
        % No RDKit call is made because no SMILES lines exist.
            txtPath = tc.writeSmilesFile_({});
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyClass(result, "cell", ...
                "readSmilesList must return a cell array");
            tc.verifyEmpty(result, ...
                "readSmilesList on empty file must return empty cell");
        end

        % ------------------------------------------------------------------
        % TC19: readSmilesList comment-only file -- no RDKit required
        % ------------------------------------------------------------------

        function test_readSmilesList_commentOnlyFile_returnsEmptyCell(tc)
        % A file containing only comment lines must return empty cell {}.
        % Verifies that '#' prefix detection is correct before any RDKit call.
            txtPath = tc.writeSmilesFile_({}, "# This is a comment", "# Another comment");
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyEmpty(result, ...
                "Comment-only file must yield empty cell");
        end

        % ------------------------------------------------------------------
        % TC20: readSmilesList accepts char filePath (requires RDKit)
        % ------------------------------------------------------------------

        function test_readSmilesList_charPath_accepted(tc)
        % readSmilesList must accept a char (not just string) file path.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            txtPath = tc.writeSmilesFile_({"CCO"});
            result  = emk.io.readSmilesList(char(txtPath));
            tc.verifyEqual(numel(result), 1, ...
                "char filePath must work and return 1 molecule");
        end

        % ------------------------------------------------------------------
        % TC21-22: readSmilesList type and shape (requires RDKit)
        % ------------------------------------------------------------------

        function test_readSmilesList_validFile_returnsCell(tc)
        % readSmilesList must return a cell array.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            txtPath = tc.writeSmilesFile_({"CCO"});
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyClass(result, "cell", ...
                "readSmilesList must return a cell array");
        end

        function test_readSmilesList_result_isRowCell(tc)
        % readSmilesList result must be a 1xN row cell.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            txtPath = tc.writeSmilesFile_({"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"});
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyEqual(size(result, 1), 1, ...
                "readSmilesList result must have 1 row");
            tc.verifyEqual(size(result, 2), 3, ...
                "readSmilesList result must have 3 columns");
        end

        function test_readSmilesList_validFile_eachElementIsPythonObject(tc)
        % Each returned element must be a Python Mol object.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            txtPath = tc.writeSmilesFile_({"CCO", "c1ccccc1"});
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyGreaterThan(numel(result), 0, "Must return at least one mol");
            for i = 1:numel(result)
                tc.verifyTrue(startsWith(class(result{i}), "py."), ...
                    sprintf("mols{%d} must be a Python object", i));
            end
        end

        % ------------------------------------------------------------------
        % TC23: readSmilesList molecule count (requires RDKit)
        % ------------------------------------------------------------------

        function test_readSmilesList_countMatches_single(tc)
        % One SMILES line => one molecule returned.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            txtPath = tc.writeSmilesFile_({"CCO"});
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyEqual(numel(result), 1, "Count must be 1");
        end

        function test_readSmilesList_countMatches_multiple(tc)
        % Three SMILES lines => three molecules returned.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles  = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
            txtPath = tc.writeSmilesFile_(smiles);
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyEqual(numel(result), 3, "Count must be 3");
        end

        % ------------------------------------------------------------------
        % TC24: readSmilesList SMILES cross-validation (requires RDKit)
        % Canonical SMILES must be preserved after loading via readSmilesList.
        % ------------------------------------------------------------------

        function test_readSmilesList_smiles_crossValidation_ethanol(tc)
        % Ethanol canonical SMILES must match after load.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            txtPath  = tc.writeSmilesFile_({"CCO"});
            result   = emk.io.readSmilesList(txtPath);
            tc.verifyEqual(numel(result), 1, "Expected 1 molecule");
            expected = emk.mol.toSmiles(emk.mol.fromSmiles("CCO"));
            actual   = emk.mol.toSmiles(result{1});
            tc.verifyEqual(actual, expected, ...
                "Canonical SMILES for ethanol must be preserved");
        end

        function test_readSmilesList_smiles_crossValidation_aspirin(tc)
        % Aspirin canonical SMILES must match after load.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smi      = "CC(=O)Oc1ccccc1C(=O)O";
            txtPath  = tc.writeSmilesFile_({smi});
            result   = emk.io.readSmilesList(txtPath);
            tc.verifyEqual(numel(result), 1, "Expected 1 molecule");
            expected = emk.mol.toSmiles(emk.mol.fromSmiles(smi));
            actual   = emk.mol.toSmiles(result{1});
            tc.verifyEqual(actual, expected, ...
                "Canonical SMILES for aspirin must be preserved");
        end

        % ------------------------------------------------------------------
        % TC25: readSmilesList ignores tab-separated name column (requires RDKit)
        % ------------------------------------------------------------------

        function test_readSmilesList_tabSeparatedName_isIgnoredAndMolParsed(tc)
        % A tab-separated "SMILES<tab>name" format must parse only the SMILES.
        % The name column must be silently ignored.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            txtPath = tc.writeSmilesFileRaw_("CCO" + sprintf("\t") + "Ethanol");
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyEqual(numel(result), 1, ...
                "Tab-separated name column: must parse 1 molecule");
            expected = emk.mol.toSmiles(emk.mol.fromSmiles("CCO"));
            actual   = emk.mol.toSmiles(result{1});
            tc.verifyEqual(actual, expected, ...
                "Mol from SMILES+name line must match plain SMILES");
        end

        function test_readSmilesList_spaceSeparatedName_isIgnored(tc)
        % A space-separated "SMILES name" format must parse only the SMILES.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            txtPath = tc.writeSmilesFileRaw_("c1ccccc1 Benzene");
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyEqual(numel(result), 1, ...
                "Space-separated name column: must parse 1 molecule");
        end

        % ------------------------------------------------------------------
        % TC26: readSmilesList skips comment lines (requires RDKit)
        % ------------------------------------------------------------------

        function test_readSmilesList_commentLinesSkipped(tc)
        % Comment lines (starting with '#') must not become molecules.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            % Write file with 2 comments + 2 valid SMILES
            lines = "# header" + newline() + ...
                    "CCO" + newline() + ...
                    "# another comment" + newline() + ...
                    "c1ccccc1";
            txtPath = tc.writeSmilesFileRaw_(lines);
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyEqual(numel(result), 2, ...
                "Comments must be skipped; 2 valid SMILES must yield 2 mols");
        end

        % ------------------------------------------------------------------
        % TC27: readSmilesList skips invalid SMILES without throwing
        % (requires RDKit)
        % ------------------------------------------------------------------

        function test_readSmilesList_invalidSmiles_skippedWithoutThrow(tc)
        % An invalid SMILES line must be skipped silently (no exception).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            lines = "CCO" + newline() + "INVALID_SMILES_XYZ" + newline() + "c1ccccc1";
            txtPath = tc.writeSmilesFileRaw_(lines);
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyEqual(numel(result), 2, ...
                "Invalid SMILES must be skipped; 2 valid must remain");
        end

        % ------------------------------------------------------------------
        % TC28: readSmilesList mixed valid/invalid (requires RDKit)
        % ------------------------------------------------------------------

        function test_readSmilesList_mixedContent_correctCount(tc)
        % Mixed file: comments, blank lines, valid, invalid SMILES.
        % Only valid SMILES must appear in the result.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            lines = "# Molecule list" + newline() + ...
                    "" + newline() + ...           % blank line
                    "CCO" + newline() + ...
                    "BADSMILES" + newline() + ...
                    "# comment" + newline() + ...
                    "CC(=O)Oc1ccccc1C(=O)O";
            txtPath = tc.writeSmilesFileRaw_(lines);
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyEqual(numel(result), 2, ...
                "Mixed content: only 2 valid SMILES expected");
        end

        % ------------------------------------------------------------------
        % TC29: readSmilesList all SMILES invalid => allLinesFailed (M3-1)
        % Requires RDKit to attempt parsing (which then fails for all lines).
        % ------------------------------------------------------------------

        function test_readSmilesList_allInvalidSmiles_throwsAllLinesFailed(tc)
        % A file containing only invalid SMILES must throw allLinesFailed.
        % Verifies that partial failure (logWarn) escalates to a hard error
        % when there are no successful parses.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            lines   = "INVALID_XYZ_1" + newline() + "INVALID_XYZ_2";
            txtPath = tc.writeSmilesFileRaw_(lines);
            tc.verifyError(@() emk.io.readSmilesList(txtPath), ...
                "emk:io:readSmilesList:allLinesFailed", ...
                "All invalid SMILES must throw allLinesFailed");
        end

        function test_readSmilesList_allLinesFailed_errorMessage_containsPath(tc)
        % Error message must contain the file path for diagnosis.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            lines   = "BAD_SMILES_XYZ";
            txtPath = tc.writeSmilesFileRaw_(lines);
            ME      = tc.captureError(@() emk.io.readSmilesList(txtPath));
            tc.assertNotEmpty(ME, "Expected allLinesFailed to be thrown");
            tc.verifyEqual(string(ME.identifier), ...
                "emk:io:readSmilesList:allLinesFailed", ...
                "Error ID must be allLinesFailed");
            % Error message must mention the file path (last segment)
            [~, fname, ext] = fileparts(txtPath);
            tc.verifySubstring(ME.message, fname + ext, ...
                "Error message must contain the file name");
        end

        function test_readSmilesList_oneValidOneInvalid_doesNotThrow(tc)
        % A mixed file with at least one valid SMILES must NOT throw.
        % Regression guard: allLinesFailed must only fire when ALL lines fail.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            lines   = "CCO" + newline() + "INVALID_XYZ";
            txtPath = tc.writeSmilesFileRaw_(lines);
            result  = emk.io.readSmilesList(txtPath);
            tc.verifyEqual(numel(result), 1, ...
                "One valid + one invalid: must return 1 mol (no throw)");
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

        function sdfPath = writeTestSdf(tc, smilesList)
        % Helper: write a temporary SDF from a cell array of SMILES strings.
        % Returns the path to the written SDF file.
            mols = cell(1, numel(smilesList));
            for i = 1:numel(smilesList)
                mols{i} = emk.mol.fromSmiles(smilesList{i});
            end
            sdfPath = fullfile(tc.tmpDir, "test_input.sdf");
            emk.io.writeSdf(mols, sdfPath);
        end

        function txtPath = writeSmilesFile_(tc, smilesList, varargin)
        % Helper: write a temporary SMILES list file.
        % smilesList: cell array of SMILES strings (one per line).
        % varargin:   additional header lines (e.g., comment strings).
        %
        % Returns the path to the created text file.
            txtPath = fullfile(tc.tmpDir, "test_smiles.txt");
            fid = fopen(txtPath, "w");
            for i = 1:numel(varargin)
                fprintf(fid, "%s\n", varargin{i});
            end
            for i = 1:numel(smilesList)
                fprintf(fid, "%s\n", smilesList{i});
            end
            fclose(fid);
        end

        function txtPath = writeSmilesFileRaw_(tc, content)
        % Helper: write raw string content to a temp text file.
        % content: scalar string with newlines already embedded.
            txtPath = fullfile(tc.tmpDir, "test_smiles_raw.txt");
            fid = fopen(txtPath, "w");
            fprintf(fid, "%s", content);
            fclose(fid);
        end

    end
end
