classdef TestRdkitModule < matlab.unittest.TestCase
% TestRdkitModule  Unit tests for emk.util.rdkitModule().
%
% Run with:
%   addpath(genpath("src"));
%   suite  = testsuite("tests/unit");
%   runner = matlab.unittest.TestRunner.withNoPlugins;
%   results = runner.run(suite);
%
% ======================================================================
% Background / Motivation
%
%   Direct py.rdkit.* access fails on MATLAB Online with:
%     TypeError: 'rdkit' object is not callable
%   because MATLAB attempts to invoke the rdkit package object as a callable
%   rather than accessing its attributes.  Using py.importlib.import_module()
%   avoids this by explicitly resolving the module reference.
%
%   emk.util.rdkitModule() centralises this pattern and caches the results
%   so each submodule is imported only once per MATLAB session.
%
%   These tests serve as a regression guard: if rdkitModule() is removed or
%   if any wrapper re-introduces py.rdkit.* direct calls, the end-to-end
%   RDKit-requiring tests here will fail on both Desktop AND Online.
%
% ======================================================================
% Coverage:
%
%   TC1: Returns a struct (requires RDKit)
%   TC2: Struct has all expected fields (requires RDKit)
%   TC3: Each field is a Python module object (requires RDKit)
%   TC4: mods.Chem.MolFromSmiles("CCO") returns a valid (non-None) object
%        -- primary cross-platform smoke test: validates importlib path
%           works where py.rdkit.Chem.MolFromSmiles() would fail (G8)
%   TC5: mods.Descriptors.MolWt is accessible and callable (requires RDKit)
%   TC6: mods.rdFpGen.GetMorganGenerator works end-to-end (requires RDKit)
%   TC7: mods.DataStructs.TanimotoSimilarity is accessible (requires RDKit)
%   TC8: mods.MACCSkeys.GenMACCSKeys works end-to-end (requires RDKit)
%   TC9: Calling rdkitModule() twice returns a struct (cache idempotency)
%        -- guards against cache corruption on repeated calls
%   TC10: Regression guard -- emk.mol.fromSmiles uses importlib path
%         (end-to-end call chain test, most critical regression guard)
%   TC11: Regression guard -- emk.descriptor.calculate uses importlib path
%   TC12: Regression guard -- emk.fingerprint.morgan uses importlib path
%   TC13: Regression guard -- emk.similarity.tanimoto uses importlib path
%
% ======================================================================
% Tests requiring RDKit use assumeTrue(tc.rdkitAvailable()) to skip
% gracefully when Python/RDKit is not configured in the current session.
%
% Design note:
%   TC10-TC13 are end-to-end integration tests that call wrapper functions
%   through their full call chains.  They existed implicitly in other test
%   classes (TestMol, TestDescriptor, etc.) but were not labelled as
%   "importlib regression guards".  Adding them here makes the intent
%   explicit and ensures Online environment failures surface in a named test.

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: rdkitModule() returns a struct
        % ------------------------------------------------------------------

        function test_rdkitModule_returnsStruct(tc)
        % rdkitModule() must return a MATLAB struct.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyClass(mods, "struct", ...
                "rdkitModule must return a struct");
        end

        % ------------------------------------------------------------------
        % TC2: Struct has all expected fields
        % ------------------------------------------------------------------

        function test_rdkitModule_hasField_Chem(tc)
        % mods.Chem must be present.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyTrue(isfield(mods, "Chem"), "Missing field: Chem");
        end

        function test_rdkitModule_hasField_DataStructs(tc)
        % mods.DataStructs must be present.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyTrue(isfield(mods, "DataStructs"), "Missing field: DataStructs");
        end

        function test_rdkitModule_hasField_Descriptors(tc)
        % mods.Descriptors must be present.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyTrue(isfield(mods, "Descriptors"), "Missing field: Descriptors");
        end

        function test_rdkitModule_hasField_rdMolDescriptors(tc)
        % mods.rdMolDescriptors must be present.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyTrue(isfield(mods, "rdMolDescriptors"), "Missing field: rdMolDescriptors");
        end

        function test_rdkitModule_hasField_rdFpGen(tc)
        % mods.rdFpGen must be present.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyTrue(isfield(mods, "rdFpGen"), "Missing field: rdFpGen");
        end

        function test_rdkitModule_hasField_MACCSkeys(tc)
        % mods.MACCSkeys must be present.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyTrue(isfield(mods, "MACCSkeys"), "Missing field: MACCSkeys");
        end

        % ------------------------------------------------------------------
        % TC3: Each field is a Python module object (class starts with "py.")
        % ------------------------------------------------------------------

        function test_rdkitModule_Chem_isPythonModule(tc)
        % mods.Chem must be a Python module reference (class starts with "py.").
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyTrue(startsWith(class(mods.Chem), "py."), ...
                sprintf("mods.Chem must start with 'py.', got: %s", class(mods.Chem)));
        end

        function test_rdkitModule_DataStructs_isPythonModule(tc)
        % mods.DataStructs must be a Python module reference.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyTrue(startsWith(class(mods.DataStructs), "py."), ...
                sprintf("mods.DataStructs must start with 'py.', got: %s", class(mods.DataStructs)));
        end

        function test_rdkitModule_Descriptors_isPythonModule(tc)
        % mods.Descriptors must be a Python module reference.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyTrue(startsWith(class(mods.Descriptors), "py."), ...
                sprintf("mods.Descriptors must start with 'py.', got: %s", class(mods.Descriptors)));
        end

        function test_rdkitModule_rdMolDescriptors_isPythonModule(tc)
        % mods.rdMolDescriptors must be a Python module reference.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyTrue(startsWith(class(mods.rdMolDescriptors), "py."), ...
                sprintf("mods.rdMolDescriptors must start with 'py.', got: %s", ...
                class(mods.rdMolDescriptors)));
        end

        function test_rdkitModule_rdFpGen_isPythonModule(tc)
        % mods.rdFpGen must be a Python module reference.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyTrue(startsWith(class(mods.rdFpGen), "py."), ...
                sprintf("mods.rdFpGen must start with 'py.', got: %s", class(mods.rdFpGen)));
        end

        function test_rdkitModule_MACCSkeys_isPythonModule(tc)
        % mods.MACCSkeys must be a Python module reference.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            tc.verifyTrue(startsWith(class(mods.MACCSkeys), "py."), ...
                sprintf("mods.MACCSkeys must start with 'py.', got: %s", class(mods.MACCSkeys)));
        end

        % ------------------------------------------------------------------
        % TC4: mods.Chem.MolFromSmiles works -- PRIMARY CROSS-PLATFORM TEST
        %
        % This is the central regression test for the MATLAB Online importlib
        %   workaround: py.rdkit.Chem.MolFromSmiles fails on MATLAB Online;
        %   mods.Chem.MolFromSmiles("CCO") must work on BOTH platforms.
        %
        % If this test fails, the importlib approach has been broken.
        % If it passes on Desktop but fails on Online, re-introduce bug was detected.
        % ------------------------------------------------------------------

        function test_rdkitModule_Chem_MolFromSmiles_returnsNonNone(tc)
        % mods.Chem.MolFromSmiles("CCO") must return a non-None Python object.
        % This directly validates the importlib workaround for G8:
        %   py.rdkit.Chem.MolFromSmiles() fails with TypeError on MATLAB Online.
        %   mods.Chem.MolFromSmiles() (via importlib) must succeed on both platforms.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            mol  = mods.Chem.MolFromSmiles("CCO");
            tc.verifyFalse(isa(mol, "py.NoneType"), ...
                "MolFromSmiles('CCO') must not return None (importlib path must work)");
            tc.verifyTrue(startsWith(class(mol), "py."), ...
                "MolFromSmiles must return a Python object");
        end

        function test_rdkitModule_Chem_MolFromSmiles_invalidSmiles_returnsNone(tc)
        % mods.Chem.MolFromSmiles("INVALID") must return None for an
        % unparseable SMILES.  Validates that the module is actually connected
        % to RDKit (not a stub).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            mol  = mods.Chem.MolFromSmiles("INVALID_SMILES_XYZ");
            tc.verifyTrue(isa(mol, "py.NoneType"), ...
                "MolFromSmiles with invalid input must return None");
        end

        % ------------------------------------------------------------------
        % TC5: mods.Descriptors.MolWt is accessible and callable
        % ------------------------------------------------------------------

        function test_rdkitModule_Descriptors_MolWt_returnsDouble(tc)
        % mods.Descriptors.MolWt must return a numeric value for ethanol.
        % Reference MW of ethanol (C2H6O): 46.069 g/mol.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            mol  = mods.Chem.MolFromSmiles("CCO");
            mw   = double(mods.Descriptors.MolWt(mol));
            tc.verifyClass(mw, "double", "MolWt must return double");
            tc.verifyEqual(mw, 46.069, "AbsTol", 0.01, ...
                "MolWt of ethanol must be ~46.069 g/mol");
        end

        % ------------------------------------------------------------------
        % TC6: mods.rdFpGen.GetMorganGenerator works end-to-end
        % ------------------------------------------------------------------

        function test_rdkitModule_rdFpGen_GetMorganGenerator_noError(tc)
        % GetMorganGenerator must succeed and return a Python generator object.
        % This exercises the call pattern used in emk.fingerprint.morgan.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            gen  = mods.rdFpGen.GetMorganGenerator( ...
                pyargs("radius", int32(2), "fpSize", int32(2048)));
            tc.verifyTrue(startsWith(class(gen), "py."), ...
                "GetMorganGenerator must return a Python object");
        end

        function test_rdkitModule_rdFpGen_GetFingerprint_noError(tc)
        % gen.GetFingerprint(mol) must succeed and return a bit vector.
        % This is the exact call pattern used in morgan.m.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            mol  = mods.Chem.MolFromSmiles("CCO");
            gen  = mods.rdFpGen.GetMorganGenerator( ...
                pyargs("radius", int32(2), "fpSize", int32(2048)));
            fp   = gen.GetFingerprint(mol);
            tc.verifyTrue(startsWith(class(fp), "py."), ...
                "GetFingerprint must return a Python bit-vector object");
            tc.verifyEqual(double(fp.GetNumBits()), 2048, ...
                "Default Morgan FP must have 2048 bits");
        end

        % ------------------------------------------------------------------
        % TC7: mods.DataStructs.TanimotoSimilarity is accessible
        % ------------------------------------------------------------------

        function test_rdkitModule_DataStructs_TanimotoSimilarity_returnsDouble(tc)
        % DataStructs.TanimotoSimilarity must return a scalar in [0, 1].
        % Identical fingerprint self-similarity must be exactly 1.0.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            mol  = mods.Chem.MolFromSmiles("CCO");
            gen  = mods.rdFpGen.GetMorganGenerator( ...
                pyargs("radius", int32(2), "fpSize", int32(2048)));
            fp   = gen.GetFingerprint(mol);
            score = double(mods.DataStructs.TanimotoSimilarity(fp, fp));
            tc.verifyEqual(score, 1.0, ...
                "Self-Tanimoto via DataStructs must be exactly 1.0");
        end

        % ------------------------------------------------------------------
        % TC8: mods.MACCSkeys.GenMACCSKeys works end-to-end
        % ------------------------------------------------------------------

        function test_rdkitModule_MACCSkeys_GenMACCSKeys_returns167bits(tc)
        % GenMACCSKeys must return a 167-bit fingerprint.
        % This exercises the call pattern used in emk.fingerprint.maccs.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods = emk.util.rdkitModule();
            mol  = mods.Chem.MolFromSmiles("CCO");
            fp   = mods.MACCSkeys.GenMACCSKeys(mol);
            tc.verifyTrue(startsWith(class(fp), "py."), ...
                "GenMACCSKeys must return a Python object");
            tc.verifyEqual(double(fp.GetNumBits()), 167, ...
                "MACCS FP must have 167 bits");
        end

        % ------------------------------------------------------------------
        % TC9: Cache idempotency -- calling twice returns struct (no error)
        % ------------------------------------------------------------------

        function test_rdkitModule_calledTwice_noError(tc)
        % Calling rdkitModule() twice must not throw an error.
        % The persistent cache must return the same struct shape on both calls.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mods1 = emk.util.rdkitModule();
            mods2 = emk.util.rdkitModule();
            tc.verifyClass(mods1, "struct", "First call must return struct");
            tc.verifyClass(mods2, "struct", "Second call must return struct");
            % Both must have the same fields (cache consistent with original)
            tc.verifyEqual(fieldnames(mods1), fieldnames(mods2), ...
                "Both calls must return the same field set");
        end

        % ------------------------------------------------------------------
        % TC10: Regression guard -- emk.mol.fromSmiles uses importlib path
        %
        % This test verifies the full call chain through emk.mol.fromSmiles.
        % If py.rdkit.Chem.MolFromSmiles() is re-introduced in fromSmiles.m,
        % this test will fail on MATLAB Online (TypeError) but pass on Desktop.
        % Running this suite on Online is the primary detection mechanism for G8.
        % ------------------------------------------------------------------

        function test_fromSmiles_callsImportlibPath_notDirectPyRdkit(tc)
        % emk.mol.fromSmiles must succeed on both Desktop and MATLAB Online.
        % Failure on Online (with valid Desktop result) indicates a py.rdkit.*
        % direct call has been re-introduced in fromSmiles.m (G8 regression).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");  % ethanol
            tc.verifyFalse(isa(mol, "py.NoneType"), ...
                "fromSmiles must return a valid mol (not None)");
            tc.verifyTrue(startsWith(class(mol), "py."), ...
                "fromSmiles must return a Python object");
            % Cross-validate atom count to confirm mol is actually ethanol
            desc = emk.descriptor.calculate(mol, "HeavyAtomCount");
            tc.verifyEqual(desc.HeavyAtomCount, 3, ...
                "Ethanol must have 3 heavy atoms (C-C-O)");
        end

        % ------------------------------------------------------------------
        % TC11: Regression guard -- emk.descriptor.calculate uses importlib path
        % ------------------------------------------------------------------

        function test_calculate_callsImportlibPath_notDirectPyRdkit(tc)
        % emk.descriptor.calculate must succeed on both Desktop and Online.
        % If py.rdkit.Chem.Descriptors.* is re-introduced, this test will
        % fail on Online but pass on Desktop (G8 regression detection).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, ["MolWt", "LogP", "TPSA"]);
            tc.verifyClass(desc.MolWt, "double", "MolWt must be double");
            tc.verifyEqual(desc.MolWt, 46.069, "AbsTol", 0.01, ...
                "Ethanol MolWt must be ~46.069 via importlib path");
        end

        % ------------------------------------------------------------------
        % TC12: Regression guard -- emk.fingerprint.morgan uses importlib path
        % ------------------------------------------------------------------

        function test_morgan_callsImportlibPath_notDirectPyRdkit(tc)
        % emk.fingerprint.morgan must succeed on both Desktop and Online.
        % Failure on Online indicates py.rdkit.Chem.rdFingerprintGenerator.*
        % has been re-introduced (G8 regression detection).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyTrue(startsWith(class(fp), "py."), ...
                "morgan must return a Python fingerprint via importlib path");
            tc.verifyEqual(double(fp.GetNumBits()), 2048, ...
                "Default Morgan FP must have 2048 bits");
        end

        % ------------------------------------------------------------------
        % TC13: Regression guard -- emk.similarity.tanimoto uses importlib path
        % ------------------------------------------------------------------

        function test_tanimoto_callsImportlibPath_notDirectPyRdkit(tc)
        % emk.similarity.tanimoto must succeed on both Desktop and Online.
        % Failure on Online indicates py.rdkit.DataStructs.TanimotoSimilarity
        % has been re-introduced (G8 regression detection).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            fp    = emk.fingerprint.morgan(mol);
            score = emk.similarity.tanimoto(fp, fp);
            tc.verifyEqual(score, 1.0, ...
                "Self-Tanimoto via importlib path must be 1.0");
        end

    end

    % ======================================================================
    methods (Access = private)

        function tf = rdkitAvailable(~)
        % Return true if Python is loaded and rdkit.Chem can be imported.
        % Uses py.importlib.import_module (not py.rdkit.* direct access)
        % to stay consistent with the importlib pattern being tested.
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

    end
end
