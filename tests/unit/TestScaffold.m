classdef TestScaffold < matlab.unittest.TestCase
% TestScaffold  Unit tests for emk.scaffold.* (L01).
%
%   Covers: genericMurcko, brics, rgroup
%   Structure: RDKit-free validation tests first, then RDKit integration tests.
%
%   Run:
%     addpath(genpath("src")); addpath(genpath("tests"));
%     suite = testsuite("tests/unit/TestScaffold");
%     runner = matlab.unittest.TestRunner.withNoPlugins;
%     results = runner.run(suite);

    methods (Access = private)
        function tf = rdkitAvailable(~)
            try
                py.importlib.import_module("rdkit.Chem");
                tf = true;
            catch
                tf = false;
            end
        end

        function mol = aspirin(~)
            mol = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
        end

        function mol = ethanol(~)
            mol = emk.mol.fromSmiles("CCO");
        end

        function mol = imatinib(~)
            % Imatinib (Gleevec): complex multi-ring scaffold
            mol = emk.mol.fromSmiles( ...
                "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1");
        end
    end

    % ======================================================================
    % TC-1: genericMurcko -- input validation (RDKit NOT required)
    % ======================================================================
    methods (Test)
        function test_genericMurcko_nonMol_throwsInvalidInput(tc)
            tc.verifyError(@() emk.scaffold.genericMurcko("CCO"), ...
                "emk:scaffold:genericMurcko:invalidInput");
        end

        function test_genericMurcko_double_throwsInvalidInput(tc)
            tc.verifyError(@() emk.scaffold.genericMurcko(42), ...
                "emk:scaffold:genericMurcko:invalidInput");
        end

        function test_genericMurcko_struct_throwsInvalidInput(tc)
            tc.verifyError(@() emk.scaffold.genericMurcko(struct("a", 1)), ...
                "emk:scaffold:genericMurcko:invalidInput");
        end
    end

    % ======================================================================
    % TC-2: brics -- input validation (RDKit NOT required)
    % ======================================================================
    methods (Test)
        function test_brics_nonMol_throwsInvalidInput(tc)
            tc.verifyError(@() emk.scaffold.brics("CCO"), ...
                "emk:scaffold:brics:invalidInput");
        end

        function test_brics_double_throwsInvalidInput(tc)
            tc.verifyError(@() emk.scaffold.brics(42), ...
                "emk:scaffold:brics:invalidInput");
        end
    end

    % ======================================================================
    % TC-3: rgroup -- input validation (RDKit NOT required)
    % ======================================================================
    methods (Test)
        function test_rgroup_nonCell_throwsInvalidInput(tc)
            tc.verifyError(@() emk.scaffold.rgroup("CC", "c1ccccc1"), ...
                "emk:scaffold:rgroup:invalidInput");
        end

        function test_rgroup_emptyCell_throwsInvalidInput(tc)
            tc.verifyError(@() emk.scaffold.rgroup({}, "c1ccccc1"), ...
                "emk:scaffold:rgroup:invalidInput");
        end

        function test_rgroup_emptyCoreSmiles_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = tc.aspirin();
            tc.verifyError(@() emk.scaffold.rgroup({mol}, ""), ...
                "emk:scaffold:rgroup:invalidInput");
        end

        function test_rgroup_nonStringCore_throwsInvalidInput(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = tc.aspirin();
            tc.verifyError(@() emk.scaffold.rgroup({mol}, 42), ...
                "emk:scaffold:rgroup:invalidInput");
        end
    end

    % ======================================================================
    % TC-4: genericMurcko RDKit integration tests
    % ======================================================================
    methods (Test)
        function test_genericMurcko_aspirin_returnsMol(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            scaf = emk.scaffold.genericMurcko(tc.aspirin());
            tc.verifyTrue(isa(scaf, "py.rdkit.Chem.rdchem.Mol"), ...
                "Return type must be py.rdkit.Chem.rdchem.Mol");
        end

        function test_genericMurcko_aspirin_hasAtoms(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            scaf = emk.scaffold.genericMurcko(tc.aspirin());
            tc.verifyGreaterThan(double(scaf.GetNumAtoms()), 0, ...
                "Aspirin scaffold must have at least 1 atom");
        end

        function test_genericMurcko_aspirin_allCarbonAtoms(tc)
            % Generic scaffold must contain only carbon atoms (atomic number 6)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            scaf = emk.scaffold.genericMurcko(tc.aspirin());
            mods = emk.util.rdkitModule();
            smiles = string(mods.Chem.MolToSmiles(scaf));
            % Only uppercase/lowercase C in SMILES (no heteroatoms N,O,S,etc.)
            heteroCount = sum(ismember(char(smiles), "NOSPHFIB"));
            tc.verifyEqual(heteroCount, 0, ...
                "Generic scaffold SMILES must not contain heteroatoms");
        end

        function test_genericMurcko_acyclic_returnsZeroAtomMol(tc)
            % Ethanol has no rings -> Murcko scaffold is empty (0 atoms)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            scaf = emk.scaffold.genericMurcko(tc.ethanol());
            tc.verifyEqual(double(scaf.GetNumAtoms()), 0, ...
                "Acyclic molecule must yield 0-atom generic scaffold");
        end

        function test_genericMurcko_quinoline_onlyCarbon(tc)
            % Quinoline has a ring N which should be replaced with C
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1cc2ccccc2nc1");  % quinoline
            scaf = emk.scaffold.genericMurcko(mol);
            mods = emk.util.rdkitModule();
            smiles = string(mods.Chem.MolToSmiles(scaf));
            tc.verifyGreaterThan(double(scaf.GetNumAtoms()), 0, ...
                "Quinoline scaffold must have atoms");
            % Verify no nitrogen in SMILES
            tc.verifyFalse(contains(smiles, "N") || contains(smiles, "n"), ...
                "Generic scaffold from quinoline must not contain N");
        end

        function test_genericMurcko_differFromMolScaffold(tc)
            % genericMurcko differs from emk.mol.scaffold for N-heterocycles
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol     = emk.mol.fromSmiles("c1cc2ccccc2nc1");  % quinoline
            bemisScaf   = emk.mol.scaffold(mol);
            genericScaf = emk.scaffold.genericMurcko(mol);
            bemisSmiles   = emk.mol.toSmiles(bemisScaf);
            genericSmiles = emk.mol.toSmiles(genericScaf);
            % Bemis-Murcko retains N; generic replaces with C
            tc.verifyNotEqual(bemisSmiles, genericSmiles, ...
                "Generic and Bemis-Murcko scaffolds must differ for N-heterocycle");
        end
    end

    % ======================================================================
    % TC-5: brics RDKit integration tests
    % ======================================================================
    methods (Test)
        function test_brics_aspirin_returnsStringArray(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            frags = emk.scaffold.brics(tc.aspirin());
            tc.verifyClass(frags, "string", ...
                "brics must return a string array");
        end

        function test_brics_aspirin_atLeastOneFragment(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            frags = emk.scaffold.brics(tc.aspirin());
            tc.verifyGreaterThanOrEqual(numel(frags), 1, ...
                "Must return at least one fragment");
        end

        function test_brics_acyclic_returnsOneFragment(tc)
            % Simple acyclic molecules typically cannot be BRICS-fragmented
            % and return their own SMILES as one fragment
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");
            frags = emk.scaffold.brics(mol);
            tc.verifyGreaterThanOrEqual(numel(frags), 1);
        end

        function test_brics_aspirin_fragmensAreNonEmpty(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            frags = emk.scaffold.brics(tc.aspirin());
            for i = 1:numel(frags)
                tc.verifyGreaterThan(strlength(frags(i)), 0, ...
                    sprintf("Fragment %d must not be empty string", i));
            end
        end

        function test_brics_imatinib_multipleFragments(tc)
            % Complex molecule should yield multiple BRICS fragments
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            frags = emk.scaffold.brics(tc.imatinib());
            tc.verifyGreaterThan(numel(frags), 1, ...
                "Imatinib should yield multiple BRICS fragments");
        end

        function test_brics_isRowVector(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            frags = emk.scaffold.brics(tc.aspirin());
            tc.verifyEqual(size(frags, 1), 1, ...
                "brics must return a row vector (1xN string)");
        end
    end

    % ======================================================================
    % TC-6: rgroup RDKit integration tests
    % ======================================================================
    methods (Test)
        function test_rgroup_returnsTable(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"Cc1ccc(N)cc1", "Clc1ccc(N)cc1"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            [tbl, ~] = emk.scaffold.rgroup(mols, "c1ccccc1");
            tc.verifyClass(tbl, "table", "rgroup must return a table");
        end

        function test_rgroup_hasCoreColumn(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"Cc1ccc(N)cc1", "Clc1ccc(N)cc1"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            [tbl, ~] = emk.scaffold.rgroup(mols, "c1ccccc1");
            tc.verifyTrue(ismember("Core", tbl.Properties.VariableNames), ...
                "Table must have a Core column");
        end

        function test_rgroup_coreColumnIsString(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"Cc1ccc(N)cc1", "Clc1ccc(N)cc1"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            [tbl, ~] = emk.scaffold.rgroup(mols, "c1ccccc1");
            tc.verifyClass(tbl.Core, "string", "Core column must be string type");
        end

        function test_rgroup_unmatchedIdx_isDouble(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"Cc1ccc(N)cc1", "Clc1ccc(N)cc1"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            [~, unmatchedIdx] = emk.scaffold.rgroup(mols, "c1ccccc1");
            tc.verifyClass(unmatchedIdx, "double", ...
                "unmatchedIdx must be a double array");
        end

        function test_rgroup_allMatch_unmatchedEmpty(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"Cc1ccccc1", "Clc1ccccc1"};  % both have benzene core
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            [~, unmatchedIdx] = emk.scaffold.rgroup(mols, "c1ccccc1");
            tc.verifyEmpty(unmatchedIdx, ...
                "Both molecules match benzene core; unmatchedIdx must be empty");
        end

        function test_rgroup_noMatch_throwsNoMatch(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mols = {emk.mol.fromSmiles("CCO")};  % ethanol - no benzene core
            tc.verifyError(@() emk.scaffold.rgroup(mols, "c1ccccc1N"), ...
                "emk:scaffold:rgroup:noMatch");
        end

        function test_rgroup_invalidCore_throwsInvalidCore(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mols = {emk.mol.fromSmiles("CCO")};
            tc.verifyError(@() emk.scaffold.rgroup(mols, "not_valid_smiles!!!"), ...
                "emk:scaffold:rgroup:invalidCore");
        end
    end
end
