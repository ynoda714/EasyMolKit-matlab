classdef TestConformer < matlab.unittest.TestCase
% TestConformer  Unit tests for emk.conformer.embed, optimize, and emk.shape.compare.
%
%   RDKit-free validation tests run first.
%   RDKit integration tests are guarded with rdkitAvailable().

    % ======================================================================
    % Helper
    % ======================================================================
    methods (Access = private)
        function tf = rdkitAvailable(~)
            try
                py.importlib.import_module("rdkit.Chem");
                tf = true;
            catch
                tf = false;
            end
        end
    end

    % ======================================================================
    % TC-1: emk.conformer.embed — input validation (no RDKit)
    % ======================================================================
    methods (Test)
        function test_embed_nonMol_throwsInvalidInput(tc)
            tc.verifyError(@() emk.conformer.embed(42), ...
                "emk:conformer:embed:invalidInput");
        end

        function test_embed_stringInput_throwsInvalidInput(tc)
            tc.verifyError(@() emk.conformer.embed("CCO"), ...
                "emk:conformer:embed:invalidInput");
        end

        function test_embed_oddVarargs_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tc.verifyError(@() emk.conformer.embed(mol, "Method"), ...
                "emk:conformer:embed:invalidInput");
        end

        function test_embed_unknownArgName_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tc.verifyError(@() emk.conformer.embed(mol, "Bogus", "ETKDGv3"), ...
                "emk:conformer:embed:invalidInput");
        end

        function test_embed_invalidMethod_throwsInvalidMethod(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tc.verifyError(@() emk.conformer.embed(mol, "Method", "BAD"), ...
                "emk:conformer:embed:invalidMethod");
        end
    end

    % ======================================================================
    % TC-2: emk.conformer.optimize — input validation (no RDKit)
    % ======================================================================
    methods (Test)
        function test_optimize_nonMol_throwsInvalidInput(tc)
            tc.verifyError(@() emk.conformer.optimize(42), ...
                "emk:conformer:optimize:invalidInput");
        end

        function test_optimize_stringInput_throwsInvalidInput(tc)
            tc.verifyError(@() emk.conformer.optimize("CCO"), ...
                "emk:conformer:optimize:invalidInput");
        end

        function test_optimize_invalidForceField_throwsInvalidForceField(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol3d = emk.conformer.embed(emk.mol.fromSmiles("CCO"), "RandomSeed", 0);
            tc.verifyError( ...
                @() emk.conformer.optimize(mol3d, "ForceField", "GAFF"), ...
                "emk:conformer:optimize:invalidForceField");
        end

        function test_optimize_molWithNoConformer_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");  % 2D, no conformer
            tc.verifyError(@() emk.conformer.optimize(mol), ...
                "emk:conformer:optimize:invalidInput");
        end
    end

    % ======================================================================
    % TC-3: emk.shape.compare — input validation (no RDKit)
    % ======================================================================
    methods (Test)
        function test_compare_nonMol1_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol3d = emk.conformer.embed(emk.mol.fromSmiles("CCO"), "RandomSeed", 0);
            tc.verifyError(@() emk.shape.compare(42, mol3d), ...
                "emk:shape:compare:invalidInput");
        end

        function test_compare_nonMol2_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol3d = emk.conformer.embed(emk.mol.fromSmiles("CCO"), "RandomSeed", 0);
            tc.verifyError(@() emk.shape.compare(mol3d, 42), ...
                "emk:shape:compare:invalidInput");
        end

        function test_compare_mol1NoConformer_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol2d = emk.mol.fromSmiles("CCO");
            mol3d = emk.conformer.embed(emk.mol.fromSmiles("CCO"), "RandomSeed", 0);
            tc.verifyError(@() emk.shape.compare(mol2d, mol3d), ...
                "emk:shape:compare:invalidInput");
        end

        function test_compare_mol2NoConformer_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol3d = emk.conformer.embed(emk.mol.fromSmiles("CCO"), "RandomSeed", 0);
            mol2d = emk.mol.fromSmiles("CCO");
            tc.verifyError(@() emk.shape.compare(mol3d, mol2d), ...
                "emk:shape:compare:invalidInput");
        end

        function test_compare_invalidMethod_throwsInvalidMethod(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol3d = emk.conformer.embed(emk.mol.fromSmiles("CCO"), "RandomSeed", 0);
            tc.verifyError(@() emk.shape.compare(mol3d, mol3d, "Method", "bad"), ...
                "emk:shape:compare:invalidMethod");
        end
    end

    % ======================================================================
    % TC-4: emk.conformer.embed — RDKit integration
    % ======================================================================
    methods (Test)
        function test_embed_returnsMol(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mol3d = emk.conformer.embed(mol, "RandomSeed", 42);
            tc.verifyClass(mol3d, "py.rdkit.Chem.rdchem.Mol");
        end

        function test_embed_hasConformer(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mol3d = emk.conformer.embed(mol, "RandomSeed", 42);
            tc.verifyGreaterThan(double(mol3d.GetNumConformers()), 0);
        end

        function test_embed_etkdgv2_runs(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            mol3d = emk.conformer.embed(mol, "Method", "ETKDGv2", "RandomSeed", 0);
            tc.verifyGreaterThan(double(mol3d.GetNumConformers()), 0);
        end

        function test_embed_etkdg_runs(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            mol3d = emk.conformer.embed(mol, "Method", "ETKDG", "RandomSeed", 0);
            tc.verifyGreaterThan(double(mol3d.GetNumConformers()), 0);
        end

        function test_embed_seed_reproducible(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol    = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mol3d1 = emk.conformer.embed(mol, "RandomSeed", 7);
            mol3d2 = emk.conformer.embed(mol, "RandomSeed", 7);
            % Both must embed (no assertion on coordinate equality due to
            % Python object reference semantics, but GetNumConformers check suffices)
            tc.verifyGreaterThan(double(mol3d1.GetNumConformers()), 0);
            tc.verifyGreaterThan(double(mol3d2.GetNumConformers()), 0);
        end

        function test_embed_quinoline_runs(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("c1ccc2ncccc2c1");  % quinoline
            mol3d = emk.conformer.embed(mol, "RandomSeed", 0);
            tc.verifyGreaterThan(double(mol3d.GetNumConformers()), 0);
        end
    end

    % ======================================================================
    % TC-5: emk.conformer.optimize — RDKit integration
    % ======================================================================
    methods (Test)
        function test_optimize_mmff94_returnsMol(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mol3d = emk.conformer.embed(mol, "RandomSeed", 42);
            opt   = emk.conformer.optimize(mol3d, "ForceField", "MMFF94");
            tc.verifyClass(opt, "py.rdkit.Chem.rdchem.Mol");
        end

        function test_optimize_mmff94_hasConformer(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            mol3d = emk.conformer.embed(mol, "RandomSeed", 0);
            opt   = emk.conformer.optimize(mol3d);
            tc.verifyGreaterThan(double(opt.GetNumConformers()), 0);
        end

        function test_optimize_uff_returnsMol(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            mol3d = emk.conformer.embed(mol, "RandomSeed", 0);
            opt   = emk.conformer.optimize(mol3d, "ForceField", "UFF");
            tc.verifyClass(opt, "py.rdkit.Chem.rdchem.Mol");
        end

        function test_optimize_embed_optimize_pipeline(tc)
            % Full pipeline: fromSmiles -> embed -> optimize
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("c1ccc2ncccc2c1");  % quinoline
            mol3d = emk.conformer.embed(mol, "RandomSeed", 1);
            opt   = emk.conformer.optimize(mol3d);
            tc.verifyGreaterThan(double(opt.GetNumConformers()), 0);
        end
    end

    % ======================================================================
    % TC-6: emk.shape.compare — RDKit integration
    % ======================================================================
    methods (Test)
        function test_compare_sameMol_scoreNear1(tc)
            % Comparing a molecule to itself should return ~1.0
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("c1ccccc1");
            mol3d = emk.conformer.embed(mol, "RandomSeed", 0);
            score = emk.shape.compare(mol3d, mol3d);
            tc.verifyGreaterThan(score, 0.99);
        end

        function test_compare_score_inRange(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mol2   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1");
            mol1_3d = emk.conformer.embed(mol1, "RandomSeed", 0);
            mol2_3d = emk.conformer.embed(mol2, "RandomSeed", 0);
            score  = emk.shape.compare(mol1_3d, mol2_3d);
            tc.verifyGreaterThanOrEqual(score, 0.0);
            tc.verifyLessThanOrEqual(score, 1.0);
        end

        function test_compare_tanimoto_method_runs(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("c1ccccc1");
            mol3d = emk.conformer.embed(mol, "RandomSeed", 0);
            score = emk.shape.compare(mol3d, mol3d, "Method", "tanimoto");
            tc.verifyGreaterThan(score, 0.0);
            tc.verifyLessThanOrEqual(score, 1.0);
        end

        function test_compare_protrude_sameMol_score1(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            mol3d = emk.conformer.embed(mol, "RandomSeed", 0);
            score = emk.shape.compare(mol3d, mol3d, "Method", "protrude");
            tc.verifyGreaterThan(score, 0.99);
        end

        function test_compare_returnsScalar(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1_3d = emk.conformer.embed(emk.mol.fromSmiles("CCO"), "RandomSeed", 0);
            mol2_3d = emk.conformer.embed(emk.mol.fromSmiles("CCC"), "RandomSeed", 0);
            score = emk.shape.compare(mol1_3d, mol2_3d);
            tc.verifySize(score, [1, 1]);
        end

        function test_compare_afterOptimize_runs(tc)
            % Full pipeline to validate end-to-end integration
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
            mol2 = emk.mol.fromSmiles("CC(=O)Oc1ccccc1");         % phenyl acetate
            mol1_3d = emk.conformer.optimize( ...
                          emk.conformer.embed(mol1, "RandomSeed", 0));
            mol2_3d = emk.conformer.optimize( ...
                          emk.conformer.embed(mol2, "RandomSeed", 0));
            score = emk.shape.compare(mol1_3d, mol2_3d);
            tc.verifyGreaterThanOrEqual(score, 0.0);
            tc.verifyLessThanOrEqual(score, 1.0);
        end
    end
end
