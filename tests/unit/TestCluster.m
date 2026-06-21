classdef TestCluster < matlab.unittest.TestCase
% TestCluster  Unit tests for emk.cluster.butina.
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
            tc.verifyError(@() emk.cluster.butina(42), ...
                "emk:cluster:butina:invalidInput");
        end

        function test_emptyCell_throwsInvalidInput(tc)
            tc.verifyError(@() emk.cluster.butina({}), ...
                "emk:cluster:butina:invalidInput");
        end

        function test_cellWithNonPyObj_throwsInvalidInput(tc)
            tc.verifyError(@() emk.cluster.butina({1, 2, 3}), ...
                "emk:cluster:butina:invalidInput");
        end

        function test_oddVarargs_throwsInvalidInput(tc)
            % Need a py object to pass the fp check, so skip if no RDKit
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.cluster.butina({fp}, "Threshold"), ...
                "emk:cluster:butina:invalidInput");
        end

        function test_unknownArgName_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.cluster.butina({fp}, "Bogus", 0.2), ...
                "emk:cluster:butina:invalidInput");
        end

        function test_invalidThreshold_zero_throwsInvalidThreshold(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.cluster.butina({fp}, "Threshold", 0), ...
                "emk:cluster:butina:invalidThreshold");
        end

        function test_invalidThreshold_negative_throwsInvalidThreshold(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.cluster.butina({fp}, "Threshold", -0.1), ...
                "emk:cluster:butina:invalidThreshold");
        end

        function test_invalidThreshold_greaterThanOne_throwsInvalidThreshold(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.cluster.butina({fp}, "Threshold", 1.1), ...
                "emk:cluster:butina:invalidThreshold");
        end

        function test_invalidThreshold_nan_throwsInvalidThreshold(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.cluster.butina({fp}, "Threshold", NaN), ...
                "emk:cluster:butina:invalidThreshold");
        end

        function test_invalidMetric_throwsInvalidMetric(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyError(@() emk.cluster.butina({fp}, "Metric", "dice"), ...
                "emk:cluster:butina:invalidMetric");
        end
    end

    % ======================================================================
    % TC-2: RDKit integration — clustering behaviour
    % ======================================================================
    methods (Test)
        function test_singleMol_oneCluster(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            clusters = emk.cluster.butina({fp}, "Threshold", 0.4);
            tc.verifyEqual(numel(clusters), 1);
            tc.verifyEqual(clusters{1}, 1);  % sole molecule, index=1
        end

        function test_returnIsCell(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            clusters = emk.cluster.butina({fp});
            tc.verifyClass(clusters, "cell");
        end

        function test_clusterMembers_are1Based(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"CCO", "CCCO", "c1ccccc1"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            fps    = cellfun(@emk.fingerprint.morgan, mols, "UniformOutput", false);
            clusters = emk.cluster.butina(fps, "Threshold", 0.4);
            for c = 1:numel(clusters)
                tc.verifyGreaterThanOrEqual(min(clusters{c}), 1);
                tc.verifyLessThanOrEqual(max(clusters{c}), numel(smiles));
            end
        end

        function test_totalMembersEqualsMolCount(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"CCO", "CCCO", "CCCCO", "c1ccccc1", "c1ccc2ccccc2c1"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            fps    = cellfun(@emk.fingerprint.morgan, mols, "UniformOutput", false);
            clusters = emk.cluster.butina(fps, "Threshold", 0.4);
            total = sum(cellfun(@numel, clusters));
            tc.verifyEqual(total, numel(smiles));
        end

        function test_identicalMols_oneCluster(tc)
            % Two identical SMILES must cluster together
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("c1ccccc1");
            fp  = emk.fingerprint.morgan(mol);
            clusters = emk.cluster.butina({fp, fp}, "Threshold", 0.2);
            tc.verifyEqual(numel(clusters), 1);
            tc.verifyEqual(numel(clusters{1}), 2);
        end

        function test_dissimilarMols_separateClusters(tc)
            % Ethanol and benzene are very dissimilar; low threshold -> 2 clusters
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("c1ccccc1");
            fp1  = emk.fingerprint.morgan(mol1);
            fp2  = emk.fingerprint.morgan(mol2);
            clusters = emk.cluster.butina({fp1, fp2}, "Threshold", 0.1);
            tc.verifyEqual(numel(clusters), 2);
        end

        function test_threshold1_allOneCluster(tc)
            % Threshold=1 means any molecule is a neighbour of any other
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"CCO", "c1ccccc1", "CC(=O)O"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            fps    = cellfun(@emk.fingerprint.morgan, mols, "UniformOutput", false);
            clusters = emk.cluster.butina(fps, "Threshold", 1.0);
            total = sum(cellfun(@numel, clusters));
            tc.verifyEqual(total, numel(smiles));
        end

        function test_defaultThreshold_runs(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            fp  = emk.fingerprint.morgan(mol);
            tc.verifyWarningFree(@() emk.cluster.butina({fp}));
        end
    end
end
