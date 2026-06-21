classdef TestDiversity < matlab.unittest.TestCase
% TestDiversity  Unit tests for emk.diversity.pick.
%
%   RDKit-free tests run first (input validation, argument parsing).
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
    % TC-1: Input validation — no RDKit required
    % ======================================================================
    methods (Test)
        function test_nonCell_throwsInvalidInput(tc)
            tc.verifyError(@() emk.diversity.pick(42, 1), ...
                "emk:diversity:pick:invalidInput");
        end

        function test_emptyCell_throwsInvalidInput(tc)
            tc.verifyError(@() emk.diversity.pick({}, 1), ...
                "emk:diversity:pick:invalidInput");
        end

        function test_cellWithNonPyObj_throwsInvalidInput(tc)
            tc.verifyError(@() emk.diversity.pick({1, 2, 3}, 2), ...
                "emk:diversity:pick:invalidInput");
        end

        function test_invalidN_zero_throwsInvalidN(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.diversity.pick({fp}, 0), ...
                "emk:diversity:pick:invalidN");
        end

        function test_invalidN_exceedsMolCount_throwsInvalidN(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.diversity.pick({fp}, 2), ...
                "emk:diversity:pick:invalidN");
        end

        function test_invalidN_fractional_throwsInvalidN(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.diversity.pick({fp}, 0.5), ...
                "emk:diversity:pick:invalidN");
        end

        function test_invalidN_nan_throwsInvalidN(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.diversity.pick({fp}, NaN), ...
                "emk:diversity:pick:invalidN");
        end

        function test_invalidMetric_throwsInvalidMetric(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.diversity.pick({fp}, 1, "Metric", "dice"), ...
                "emk:diversity:pick:invalidMetric");
        end

        function test_oddVarargs_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.diversity.pick({fp}, 1, "Seed"), ...
                "emk:diversity:pick:invalidInput");
        end

        function test_unknownArgName_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.diversity.pick({fp}, 1, "Bogus", 1), ...
                "emk:diversity:pick:invalidInput");
        end
    end

    % ======================================================================
    % TC-2: RDKit integration — diversity picking behaviour
    % ======================================================================
    methods (Test)
        function test_returnIsDouble(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"CCO", "c1ccccc1", "CC(=O)O"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            fps    = cellfun(@emk.fingerprint.morgan, mols, "UniformOutput", false);
            indices = emk.diversity.pick(fps, 2, "Seed", 42);
            tc.verifyClass(indices, "double");
        end

        function test_returnSize_equalsN(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"CCO", "c1ccccc1", "CC(=O)O", "c1ccc2ccccc2c1"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            fps    = cellfun(@emk.fingerprint.morgan, mols, "UniformOutput", false);
            indices = emk.diversity.pick(fps, 3, "Seed", 0);
            tc.verifySize(indices, [1, 3]);
        end

        function test_indices_are1Based(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"CCO", "c1ccccc1", "CC(=O)O"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            fps    = cellfun(@emk.fingerprint.morgan, mols, "UniformOutput", false);
            indices = emk.diversity.pick(fps, 2, "Seed", 1);
            tc.verifyGreaterThanOrEqual(min(indices), 1);
            tc.verifyLessThanOrEqual(max(indices), numel(smiles));
        end

        function test_indices_unique(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"CCO", "c1ccccc1", "CC(=O)O", "CCCO", "c1ccc2ccccc2c1"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            fps    = cellfun(@emk.fingerprint.morgan, mols, "UniformOutput", false);
            indices = emk.diversity.pick(fps, 4, "Seed", 0);
            tc.verifyEqual(numel(unique(indices)), 4);
        end

        function test_pickAll_returnsAllIndices(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"CCO", "c1ccccc1", "CC(=O)O"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            fps    = cellfun(@emk.fingerprint.morgan, mols, "UniformOutput", false);
            indices = emk.diversity.pick(fps, 3, "Seed", 0);
            tc.verifyEqual(sort(indices), [1, 2, 3]);
        end

        function test_seed_reproducible(tc)
            % Same seed must produce the same indices
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"CCO", "c1ccccc1", "CC(=O)O", "CCCO", "c1ccc2ccccc2c1"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            fps    = cellfun(@emk.fingerprint.morgan, mols, "UniformOutput", false);
            idx1 = emk.diversity.pick(fps, 3, "Seed", 42);
            idx2 = emk.diversity.pick(fps, 3, "Seed", 42);
            tc.verifyEqual(idx1, idx2);
        end

        function test_pickOne_returnsScalar(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            indices = emk.diversity.pick({fp}, 1, "Seed", 0);
            tc.verifySize(indices, [1, 1]);
            tc.verifyEqual(indices, 1);
        end
    end
end
