classdef TestFilter < matlab.unittest.TestCase
% TestFilter  Unit tests for emk.filter.lipinski.
%
%   All tests that do NOT require RDKit are listed first so that
%   validation-order errors surface early.
%   Tests that require RDKit are guarded with
%       tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");

    properties (Access = private)
        % A minimal table that satisfies all Ro5 criteria (used by many tests)
        passTbl
        % A table where every row violates all 4 Ro5 criteria
        failTbl
        % A mixed table: row 1 passes, row 2 fails
        mixedTbl
    end

    methods (TestMethodSetup)
        function buildFixtures(tc)
            % Build pure-MATLAB tables (no RDKit needed) for RDKit-free tests.

            % All-pass: small drug-like molecules (MW~100, LogP~1, HBD=1, HBA=2)
            tc.passTbl = table( ...
                [150; 200], ...
                [1.5; 2.0], ...
                [1; 1], ...
                [2; 3], ...
                'VariableNames', {'MolWt', 'LogP', 'NumHDonors', 'NumHAcceptors'});

            % All-fail: huge, greasy, polar molecules exceeding all thresholds
            tc.failTbl = table( ...
                [600; 700], ...
                [6; 7], ...
                [6; 8], ...
                [11; 12], ...
                'VariableNames', {'MolWt', 'LogP', 'NumHDonors', 'NumHAcceptors'});

            % Mixed: row 1 passes all, row 2 fails all
            tc.mixedTbl = table( ...
                [150; 600], ...
                [1.5; 6.0], ...
                [1; 6], ...
                [2; 11], ...
                'VariableNames', {'MolWt', 'LogP', 'NumHDonors', 'NumHAcceptors'});
        end
    end

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
    % TC-1: Input validation (RDKit NOT required)
    % ======================================================================
    methods (Test)
        function test_nonTable_throwsInvalidInput(tc)
            % Table is the primary argument; its check fires BEFORE MaxViolations.
            tc.verifyError(@() emk.filter.lipinski(42), ...
                "emk:filter:lipinski:invalidInput");
        end

        function test_structInput_throwsInvalidInput(tc)
            s.MolWt = 100; s.LogP = 1; s.NumHDonors = 1; s.NumHAcceptors = 2;
            tc.verifyError(@() emk.filter.lipinski(s), ...
                "emk:filter:lipinski:invalidInput");
        end

        function test_cellInput_throwsInvalidInput(tc)
            tc.verifyError(@() emk.filter.lipinski({}), ...
                "emk:filter:lipinski:invalidInput");
        end

        function test_missingColumn_MolWt_throwsMissingColumns(tc)
            t = tc.passTbl;
            t = removevars(t, "MolWt");
            tc.verifyError(@() emk.filter.lipinski(t), ...
                "emk:filter:lipinski:missingColumns");
        end

        function test_missingColumn_LogP_throwsMissingColumns(tc)
            t = tc.passTbl;
            t = removevars(t, "LogP");
            tc.verifyError(@() emk.filter.lipinski(t), ...
                "emk:filter:lipinski:missingColumns");
        end

        function test_missingColumn_NumHDonors_throwsMissingColumns(tc)
            t = tc.passTbl;
            t = removevars(t, "NumHDonors");
            tc.verifyError(@() emk.filter.lipinski(t), ...
                "emk:filter:lipinski:missingColumns");
        end

        function test_missingColumn_NumHAcceptors_throwsMissingColumns(tc)
            t = tc.passTbl;
            t = removevars(t, "NumHAcceptors");
            tc.verifyError(@() emk.filter.lipinski(t), ...
                "emk:filter:lipinski:missingColumns");
        end

        function test_missingColumn_errorMessage_containsColumnName(tc)
            t = removevars(tc.passTbl, "MolWt");
            try
                emk.filter.lipinski(t);
                tc.verifyFail("Expected error was not thrown.");
            catch ME
                tc.verifySubstring(ME.message, "MolWt");
            end
        end

        function test_unknownArgName_throwsInvalidInput(tc)
            tc.verifyError(@() emk.filter.lipinski(tc.passTbl, "Bogus", 0), ...
                "emk:filter:lipinski:invalidInput");
        end

        function test_oddNumberOfVarargs_throwsInvalidInput(tc)
            tc.verifyError(@() emk.filter.lipinski(tc.passTbl, "MaxViolations"), ...
                "emk:filter:lipinski:invalidInput");
        end

        function test_maxViol_negative_throwsInvalidMaxViol(tc)
            tc.verifyError( ...
                @() emk.filter.lipinski(tc.passTbl, "MaxViolations", -1), ...
                "emk:filter:lipinski:invalidMaxViol");
        end

        function test_maxViol_five_throwsInvalidMaxViol(tc)
            tc.verifyError( ...
                @() emk.filter.lipinski(tc.passTbl, "MaxViolations", 5), ...
                "emk:filter:lipinski:invalidMaxViol");
        end

        function test_maxViol_fractional_throwsInvalidMaxViol(tc)
            tc.verifyError( ...
                @() emk.filter.lipinski(tc.passTbl, "MaxViolations", 0.5), ...
                "emk:filter:lipinski:invalidMaxViol");
        end

        function test_maxViol_nan_throwsInvalidMaxViol(tc)
            tc.verifyError( ...
                @() emk.filter.lipinski(tc.passTbl, "MaxViolations", NaN), ...
                "emk:filter:lipinski:invalidMaxViol");
        end

        function test_maxViol_inf_throwsInvalidMaxViol(tc)
            tc.verifyError( ...
                @() emk.filter.lipinski(tc.passTbl, "MaxViolations", Inf), ...
                "emk:filter:lipinski:invalidMaxViol");
        end

        function test_nonTable_with_badMaxViol_throwsInvalidInput(tc)
            % Table validation fires BEFORE MaxViolations validation.
            % So lipinski(42, MaxViolations=-1) → invalidInput, not invalidMaxViol.
            tc.verifyError( ...
                @() emk.filter.lipinski(42, "MaxViolations", -1), ...
                "emk:filter:lipinski:invalidInput");
        end
    end

    % ======================================================================
    % TC-2: Output schema (RDKit NOT required)
    % ======================================================================
    methods (Test)
        function test_returnIsTable(tc)
            result = emk.filter.lipinski(tc.passTbl);
            tc.verifyClass(result, "table");
        end

        function test_returnHasPassRo5Column(tc)
            result = emk.filter.lipinski(tc.passTbl);
            tc.verifyTrue(ismember("Pass_Ro5", result.Properties.VariableNames));
        end

        function test_returnHasViolationsRo5Column(tc)
            result = emk.filter.lipinski(tc.passTbl);
            tc.verifyTrue(ismember("Violations_Ro5", result.Properties.VariableNames));
        end

        function test_passRo5Column_isLogical(tc)
            result = emk.filter.lipinski(tc.passTbl);
            tc.verifyClass(result.Pass_Ro5, "logical");
        end

        function test_violationsRo5Column_isDouble(tc)
            result = emk.filter.lipinski(tc.passTbl);
            tc.verifyClass(result.Violations_Ro5, "double");
        end

        function test_rowCount_unchanged(tc)
            result = emk.filter.lipinski(tc.mixedTbl);
            tc.verifyEqual(height(result), height(tc.mixedTbl));
        end

        function test_originalColumns_preserved(tc)
            result = emk.filter.lipinski(tc.passTbl);
            originalCols = tc.passTbl.Properties.VariableNames;
            for k = 1:numel(originalCols)
                tc.verifyTrue(ismember(originalCols{k}, result.Properties.VariableNames));
            end
        end

        function test_passRo5_isColumnVector(tc)
            result = emk.filter.lipinski(tc.passTbl);
            tc.verifySize(result.Pass_Ro5, [height(tc.passTbl), 1]);
        end

        function test_violationsRo5_isColumnVector(tc)
            result = emk.filter.lipinski(tc.passTbl);
            tc.verifySize(result.Violations_Ro5, [height(tc.passTbl), 1]);
        end

        function test_newColumns_appendedAtEnd(tc)
            % Pass_Ro5 and Violations_Ro5 must be the LAST two columns.
            % Downstream code that uses column indices depends on this.
            % VariableNames returns cell-of-char, so compare with char literals.
            result = emk.filter.lipinski(tc.passTbl);
            names = result.Properties.VariableNames;
            n = numel(names);
            tc.verifyEqual(names{n-1}, 'Violations_Ro5');
            tc.verifyEqual(names{n},   'Pass_Ro5');
        end

        function test_columnCount_increasesByTwo(tc)
            % Exactly two new columns must be added.
            result = emk.filter.lipinski(tc.passTbl);
            tc.verifyEqual(width(result), width(tc.passTbl) + 2);
        end
    end

    % ======================================================================
    % TC-3: Filter logic — all-pass fixture (RDKit NOT required)
    % ======================================================================
    methods (Test)
        function test_allPass_passRo5_allTrue(tc)
            result = emk.filter.lipinski(tc.passTbl);
            tc.verifyTrue(all(result.Pass_Ro5));
        end

        function test_allPass_violations_allZero(tc)
            result = emk.filter.lipinski(tc.passTbl);
            tc.verifyEqual(result.Violations_Ro5, zeros(height(tc.passTbl), 1));
        end
    end

    % ======================================================================
    % TC-4: Filter logic — all-fail fixture (RDKit NOT required)
    % ======================================================================
    methods (Test)
        function test_allFail_passRo5_allFalse(tc)
            result = emk.filter.lipinski(tc.failTbl);
            tc.verifyTrue(~any(result.Pass_Ro5));
        end

        function test_allFail_violations_allFour(tc)
            result = emk.filter.lipinski(tc.failTbl);
            tc.verifyEqual(result.Violations_Ro5, 4 * ones(height(tc.failTbl), 1));
        end
    end

    % ======================================================================
    % TC-5: Filter logic — mixed fixture (RDKit NOT required)
    % ======================================================================
    methods (Test)
        function test_mixed_passRo5_correctPattern(tc)
            result = emk.filter.lipinski(tc.mixedTbl);
            tc.verifyEqual(result.Pass_Ro5, [true; false]);
        end

        function test_mixed_violations_correctCounts(tc)
            result = emk.filter.lipinski(tc.mixedTbl);
            tc.verifyEqual(result.Violations_Ro5(1), 0);
            tc.verifyEqual(result.Violations_Ro5(2), 4);
        end
    end

    % ======================================================================
    % TC-6: Boundary conditions (RDKit NOT required)
    % ======================================================================
    methods (Test)
        function test_boundaryExact_MW500_passes(tc)
            % MW == 500 exactly is a pass (criterion is > 500)
            t = table(500, 1.0, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyTrue(result.Pass_Ro5);
            tc.verifyEqual(result.Violations_Ro5, 0);
        end

        function test_boundaryExact_MW501_fails(tc)
            t = table(501, 1.0, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyFalse(result.Pass_Ro5);
            tc.verifyEqual(result.Violations_Ro5, 1);
        end

        function test_boundaryExact_LogP5_passes(tc)
            t = table(300, 5.0, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyTrue(result.Pass_Ro5);
        end

        function test_boundaryExact_LogP5p1_fails(tc)
            t = table(300, 5.1, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyFalse(result.Pass_Ro5);
        end

        function test_boundaryExact_HBD5_passes(tc)
            t = table(300, 1.0, 5, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyTrue(result.Pass_Ro5);
        end

        function test_boundaryExact_HBD6_fails(tc)
            t = table(300, 1.0, 6, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyFalse(result.Pass_Ro5);
        end

        function test_boundaryExact_HBA10_passes(tc)
            t = table(300, 1.0, 1, 10, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyTrue(result.Pass_Ro5);
        end

        function test_boundaryExact_HBA11_fails(tc)
            t = table(300, 1.0, 1, 11, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyFalse(result.Pass_Ro5);
        end

        function test_singleViolation_violationsIsOne(tc)
            % Only MW violates
            t = table(501, 1.0, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyEqual(result.Violations_Ro5, 1);
        end

        function test_singleViolation_logP_isOne(tc)
            % Only LogP violates (MW, HBD, HBA all pass)
            t = table(300, 5.1, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyEqual(result.Violations_Ro5, 1);
            tc.verifyFalse(result.Pass_Ro5);
        end

        function test_singleViolation_hbd_isOne(tc)
            % Only NumHDonors violates
            t = table(300, 1.0, 6, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyEqual(result.Violations_Ro5, 1);
            tc.verifyFalse(result.Pass_Ro5);
        end

        function test_singleViolation_hba_isOne(tc)
            % Only NumHAcceptors violates
            t = table(300, 1.0, 1, 11, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyEqual(result.Violations_Ro5, 1);
            tc.verifyFalse(result.Pass_Ro5);
        end

        function test_twoViolations_isTwo(tc)
            % MW and LogP both violate; HBD and HBA pass
            t = table(501, 5.1, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyEqual(result.Violations_Ro5, 2);
        end

        function test_violations_inRange(tc)
            % Violations_Ro5 must always be in [0, 4] for all fixture rows.
            for tbl = {tc.passTbl, tc.failTbl, tc.mixedTbl}
                result = emk.filter.lipinski(tbl{1});
                tc.verifyGreaterThanOrEqual(min(result.Violations_Ro5), 0);
                tc.verifyLessThanOrEqual(   max(result.Violations_Ro5), 4);
            end
        end
    end

    % ======================================================================
    % TC-7: MaxViolations option (RDKit NOT required)
    % ======================================================================
    methods (Test)
        function test_maxViol0_oneViolation_fails(tc)
            % MW=501 => 1 violation; MaxViolations=0 => fail
            t = table(501, 1.0, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t, "MaxViolations", 0);
            tc.verifyFalse(result.Pass_Ro5);
        end

        function test_maxViol1_oneViolation_passes(tc)
            % MW=501 => 1 violation; MaxViolations=1 => pass
            t = table(501, 1.0, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t, "MaxViolations", 1);
            tc.verifyTrue(result.Pass_Ro5);
        end

        function test_maxViol1_twoViolations_fails(tc)
            % MW=501, LogP=6 => 2 violations; MaxViolations=1 => fail
            t = table(501, 6.0, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t, "MaxViolations", 1);
            tc.verifyFalse(result.Pass_Ro5);
        end

        function test_maxViol4_allViolations_passes(tc)
            % 4 violations; MaxViolations=4 => pass
            result = emk.filter.lipinski(tc.failTbl, "MaxViolations", 4);
            tc.verifyTrue(all(result.Pass_Ro5));
        end

        function test_maxViolDefault_equalsZero(tc)
            % Default MaxViolations=0 must behave like explicit MaxViolations=0
            t = table(501, 1.0, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            r0 = emk.filter.lipinski(t);
            r1 = emk.filter.lipinski(t, "MaxViolations", 0);
            tc.verifyEqual(r0.Pass_Ro5,       r1.Pass_Ro5);
            tc.verifyEqual(r0.Violations_Ro5, r1.Violations_Ro5);
        end

        function test_maxViol2_twoViolations_passes(tc)
            % MW=501, LogP=5.1 => 2 violations; MaxViolations=2 => pass
            t = table(501, 5.1, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t, "MaxViolations", 2);
            tc.verifyTrue(result.Pass_Ro5);
        end

        function test_maxViol2_threeViolations_fails(tc)
            % MW=501, LogP=5.1, HBD=6 => 3 violations; MaxViolations=2 => fail
            t = table(501, 5.1, 6, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t, "MaxViolations", 2);
            tc.verifyFalse(result.Pass_Ro5);
        end

        function test_maxViol3_threeViolations_passes(tc)
            % 3 violations; MaxViolations=3 => pass
            t = table(501, 5.1, 6, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t, "MaxViolations", 3);
            tc.verifyTrue(result.Pass_Ro5);
        end
    end

    % ======================================================================
    % TC-7b: Edge case — 0-row table (RDKit NOT required)
    % ======================================================================
    methods (Test)
        function test_emptyTable_returns0Rows(tc)
            % A 0-row table must be returned unchanged in row count.
            empty = tc.passTbl([], :);   % 0 rows, same columns
            result = emk.filter.lipinski(empty);
            tc.verifyEqual(height(result), 0);
        end

        function test_emptyTable_hasPassRo5Column(tc)
            empty  = tc.passTbl([], :);
            result = emk.filter.lipinski(empty);
            tc.verifyTrue(ismember("Pass_Ro5", result.Properties.VariableNames));
        end

        function test_emptyTable_passRo5_isLogicalEmpty(tc)
            empty  = tc.passTbl([], :);
            result = emk.filter.lipinski(empty);
            tc.verifyClass(result.Pass_Ro5, "logical");
            tc.verifySize(result.Pass_Ro5, [0, 1]);
        end

        function test_emptyTable_violations_isDoubleEmpty(tc)
            empty  = tc.passTbl([], :);
            result = emk.filter.lipinski(empty);
            tc.verifyClass(result.Violations_Ro5, "double");
            tc.verifySize(result.Violations_Ro5, [0, 1]);
        end
    end

    % ======================================================================
    % TC-8: NaN handling (RDKit NOT required)
    % NOTE: MATLAB evaluates NaN > threshold as false.
    %       Therefore NaN is treated as a NON-violation for that criterion.
    %       This is a footgun: NaN molecules may spuriously appear to pass Ro5.
    % ======================================================================
    methods (Test)
        function test_nanMolWt_treatedAsNonViolation(tc)
            % NaN > 500 is false in MATLAB, so MW does NOT count as a violation.
            t = table(NaN, 1.0, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyEqual(result.Violations_Ro5, 0);
            tc.verifyTrue(result.Pass_Ro5);
        end

        function test_nanLogP_treatedAsNonViolation(tc)
            % NaN > 5 is false; LogP NaN does NOT count as a violation.
            t = table(300, NaN, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyEqual(result.Violations_Ro5, 0);
            tc.verifyTrue(result.Pass_Ro5);
        end

        function test_nanHBD_treatedAsNonViolation(tc)
            % NaN > 5 is false; NumHDonors NaN does NOT count as a violation.
            t = table(300, 1.0, NaN, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyEqual(result.Violations_Ro5, 0);
            tc.verifyTrue(result.Pass_Ro5);
        end

        function test_nanHBA_treatedAsNonViolation(tc)
            % NaN > 10 is false; NumHAcceptors NaN does NOT count as a violation.
            t = table(300, 1.0, 1, NaN, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyEqual(result.Violations_Ro5, 0);
            tc.verifyTrue(result.Pass_Ro5);
        end

        function test_allNaN_treatedAsAllNonViolation(tc)
            % All four NaN descriptors => 0 violations (all criteria: NaN > x = false)
            t = table(NaN, NaN, NaN, NaN, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyEqual(result.Violations_Ro5, 0);
            tc.verifyTrue(result.Pass_Ro5);
        end

        function test_nanMixedWithViolation_countsOnlyRealViolation(tc)
            % NaN MolWt does NOT count; LogP=6 DOES count => 1 violation.
            t = table(NaN, 6.0, 1, 2, ...
                'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
            result = emk.filter.lipinski(t);
            tc.verifyEqual(result.Violations_Ro5, 1);
            tc.verifyFalse(result.Pass_Ro5);
        end
    end

    % ======================================================================
    % TC-9: RDKit integration tests (require Python + RDKit)
    % ======================================================================
    methods (Test)
        function test_aspirin_passesRo5(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            % Aspirin: MW=180, LogP=1.19, HBD=1, HBA=4 => all pass
            mol = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            tbl = emk.mol.toTable(mol);
            result = emk.filter.lipinski(tbl);
            tc.verifyTrue(result.Pass_Ro5);
            tc.verifyEqual(result.Violations_Ro5, 0);
        end

        function test_ethanol_passesRo5(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            % Ethanol: MW=46, LogP=-0.31, HBD=1, HBA=1 => all pass
            mol = emk.mol.fromSmiles("CCO");
            tbl = emk.mol.toTable(mol);
            result = emk.filter.lipinski(tbl);
            tc.verifyTrue(result.Pass_Ro5);
            tc.verifyEqual(result.Violations_Ro5, 0);
        end

        function test_largeMolecule_failsRo5(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            % Hexatriacontane (n-C36H74): MW ~ 507 Da => violates MW criterion
            % Simple linear alkane chosen for SMILES validity
            mol = emk.mol.fromSmiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
            tbl = emk.mol.toTable(mol);
            result = emk.filter.lipinski(tbl);
            tc.verifyFalse(result.Pass_Ro5);
            tc.verifyGreaterThan(result.Violations_Ro5, 0);
        end

        function test_violations_crossValidate_withDescriptors(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            % Manually count violations from the descriptor struct and compare
            % to lipinski's Violations_Ro5 column (cross-validation, P3)
            mol = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
            desc = emk.descriptor.calculate(mol);
            expected = double(desc.MolWt         > 500) + ...
                       double(desc.LogP          > 5  ) + ...
                       double(desc.NumHDonors    > 5  ) + ...
                       double(desc.NumHAcceptors > 10 );
            tbl    = emk.mol.toTable(mol);
            result = emk.filter.lipinski(tbl);
            tc.verifyEqual(result.Violations_Ro5, expected);
        end

        function test_passFilter_returnsSubset(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smilesList = { ...
                "CCO", ...                             % ethanol     (pass)
                "CC(=O)Oc1ccccc1C(=O)O", ...          % aspirin     (pass)
                "c1ccccc1"};                           % benzene     (pass)
            mols = cellfun(@emk.mol.fromSmiles, smilesList, ...
                "UniformOutput", false);
            tbl    = emk.mol.toTable(mols);
            result = emk.filter.lipinski(tbl);
            passing = result(result.Pass_Ro5, :);
            % All three small molecules should pass Ro5
            tc.verifyEqual(height(passing), 3);
        end

        function test_passFilter_withSomeFailing(tc)
            % Use hexatriacontane (MW~507, LogP>>5) as the failing molecule.
            % Ethanol passes, hexatriacontane fails => 1 in passing set.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mols = { ...
                emk.mol.fromSmiles("CCO"), ...
                emk.mol.fromSmiles("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")};
            tbl    = emk.mol.toTable(mols);
            result = emk.filter.lipinski(tbl);
            passing = result(result.Pass_Ro5, :);
            tc.verifyEqual(height(passing), 1);
        end

        function test_toTable_lipinski_pipeline_compatible(tc)
            % Verify that emk.mol.toTable output is directly compatible with
            % emk.filter.lipinski without any manual column manipulation.
            % This tests the M3-5 -> M3-6 integration contract.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
            tbl = emk.mol.toTable(mol);
            % Must not throw — required columns are present in default toTable output
            result = emk.filter.lipinski(tbl);
            tc.verifyClass(result, "table");
            tc.verifyTrue(ismember("Pass_Ro5", result.Properties.VariableNames));
        end
    end
end
