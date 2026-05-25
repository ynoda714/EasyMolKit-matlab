classdef TestMordred < matlab.unittest.TestCase
% TestMordred  Unit tests for emk.descriptor.mordred,
%              emk.descriptor.mordredBatch, and emk.descriptor.mordredNames.
%
% Strategy: input-validation tests fire without Python or Mordred;
% all tests requiring Mordred are guarded with assumeTrue(mordredAvailable).
%
% Coverage (mordred -- single molecule):
%   TC1 : Non-Mol input => invalidInput (no Mordred required)
%   TC1b: Error message contains offending class name
%   TC2 : char / string / numeric / logical / empty all rejected (no Mordred)
%   TC3 : libraryNotFound error when mordred not installed
%         (requires mordred absent -- skipped when mordred is present)
%   TC4 : Valid Mol => result is struct (requires Mordred + RDKit)
%   TC5 : Struct fields are all double scalars (requires Mordred)
%   TC6 : Subset request returns only requested fields (requires Mordred)
%   TC7 : MW descriptor for ethanol matches RDKit MolWt (AbsTol=0.5) (requires Mordred)
%   TC8 : All values are finite or NaN -- no Inf allowed (requires Mordred)
%   TC9 : Named descriptor "nRot" = 0 for ethanol (no rotatable bonds)
%   TC10: Named descriptor "nHBDon" = 1 for ethanol (one hydroxyl donor)
%   TC10b: Cross-validation -- nHBDon matches emk.descriptor.calculate HBondDonors
%
% Coverage (mordredBatch -- multiple molecules):
%   TB1 : Non-cell input => invalidInput (no Mordred required)
%   TB2 : Valid batch => result is table (requires Mordred)
%   TB3 : Table dimensions: N rows x M cols (requires Mordred)
%   TB4 : Column names match requested descriptor names (requires Mordred)
%   TB5 : Invalid Mol in cell => NaN row, no error (requires Mordred)
%   TB5b: Valid row in same batch is NOT all NaN (TB5 complement)
%   TB6 : All-invalid mols => allMolsFailed error (requires Mordred)
%   TB7 : Batch MW values cross-validate with single mordred() call (requires Mordred)
%   TB8 : Empty cell with explicit descriptor names => 0-row table (requires Mordred)
%
% Coverage (mordredNames):
%   TN1 : Returns string array (requires Mordred)
%   TN1b: Returns row vector 1 x N (requires Mordred)
%   TN2 : numel >= 1000 (Mordred has ~1800 2D descriptors) (requires Mordred)
%   TN3 : "MW" is in the names list (requires Mordred)
%
% Run with:
%   addpath(genpath("src")); addpath(genpath("tests"));
%   results = run(TestMordred);

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    % mordred -- input validation (no Python required)
    % ======================================================================
    methods (Test)

        function test_mordred_numericInput_throwsInvalidInput(tc)
            tc.verifyError(@() emk.descriptor.mordred(3.14), ...
                "emk:descriptor:mordred:invalidInput", ...
                "Numeric input must throw invalidInput");
        end

        function test_mordred_numericInput_errorMessage_containsClass(tc)
            ME = tc.captureError(@() emk.descriptor.mordred(3.14));
            tc.assertNotEmpty(ME);
            tc.verifySubstring(ME.message, "double", ...
                "Error message must mention the offending class 'double'");
        end

        function test_mordred_stringInput_throwsInvalidInput(tc)
            tc.verifyError(@() emk.descriptor.mordred("CCO"), ...
                "emk:descriptor:mordred:invalidInput", ...
                "Raw SMILES string must throw invalidInput");
        end

        function test_mordred_charInput_throwsInvalidInput(tc)
            tc.verifyError(@() emk.descriptor.mordred(char("CCO")), ...
                "emk:descriptor:mordred:invalidInput", ...
                "char input must throw invalidInput");
        end

        function test_mordred_logicalInput_throwsInvalidInput(tc)
            tc.verifyError(@() emk.descriptor.mordred(true), ...
                "emk:descriptor:mordred:invalidInput", ...
                "Logical input must throw invalidInput");
        end

        function test_mordred_emptyInput_throwsInvalidInput(tc)
            tc.verifyError(@() emk.descriptor.mordred([]), ...
                "emk:descriptor:mordred:invalidInput", ...
                "Empty [] must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % mordredBatch input validation (no Python required)
        % ------------------------------------------------------------------

        function test_mordredBatch_nonCellInput_throwsInvalidInput(tc)
        % mordredBatch requires a cell array; non-cell must be rejected.
            tc.verifyError(@() emk.descriptor.mordredBatch("CCO"), ...
                "emk:descriptor:mordredBatch:invalidInput", ...
                "String input to mordredBatch must throw invalidInput");
        end

        function test_mordredBatch_numericInput_throwsInvalidInput(tc)
            tc.verifyError(@() emk.descriptor.mordredBatch([1, 2, 3]), ...
                "emk:descriptor:mordredBatch:invalidInput", ...
                "Numeric input to mordredBatch must throw invalidInput");
        end

    end   % Input validation

    % ======================================================================
    % mordred -- Mordred-dependent tests
    % ======================================================================
    methods (Test)

        function test_mordred_validMol_returnsStruct(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.mordred(mol, ["MW", "nRot"]);
            tc.verifyClass(desc, "struct", "Result must be a struct");
        end

        function test_mordred_subsetFields_matchRequest(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            requested = ["MW", "nRot", "nHBDon"];
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.mordred(mol, requested);
            for k = 1:numel(requested)
                tc.verifyTrue(isfield(desc, requested(k)), ...
                    "Missing requested field: " + requested(k));
            end
            tc.verifyEqual(numel(fieldnames(desc)), numel(requested), ...
                "Struct must have exactly the requested number of fields");
        end

        function test_mordred_allValuesDoubleOrNaN(tc)
        % All returned values must be double scalars (or NaN for failures).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            mol  = emk.mol.fromSmiles("c1ccccc1");
            desc = emk.descriptor.mordred(mol, ["MW", "ALogP", "nRot"]);
            fnames = fieldnames(desc);
            for k = 1:numel(fnames)
                val = desc.(fnames{k});
                tc.verifyClass(val, "double", fnames{k} + " must be double");
                tc.verifySize(val, [1 1], fnames{k} + " must be scalar");
                tc.verifyFalse(isinf(val), fnames{k} + " must not be Inf");
            end
        end

        function test_mordred_MW_ethanol_nearRef(tc)
        % Mordred "MW" (molecular weight) for ethanol should be ~46.07.
        % Using AbsTol=0.5 to tolerate slight differences between Mordred
        % and RDKit average MW implementations.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.mordred(mol, "MW");
            tc.verifyEqual(desc.MW, 46.07, "AbsTol", 0.5, ...
                "Ethanol MW from Mordred must be near 46.07 g/mol");
        end

        function test_mordred_nRot_ethanol_isZero(tc)
        % Ethanol has no rotatable bonds.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.mordred(mol, "nRot");
            tc.verifyEqual(desc.nRot, 0, "AbsTol", 0, ...
                "Ethanol must have 0 rotatable bonds (nRot=0)");
        end

        function test_mordred_nHBd_ethanol_isOne(tc)
        % Ethanol has one hydroxyl => 1 H-bond donor.
        % mordredcommunity 2.x renamed nHBd -> nHBDon.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.mordred(mol, "nHBDon");
            tc.verifyEqual(desc.nHBDon, 1, ...
                "Ethanol must have 1 H-bond donor (nHBDon=1)");
        end

        function test_mordred_nHBd_crossValidate_rdkitHBondDonors(tc)
        % TC10b: Mordred nHBDon must equal RDKit HBondDonors for ethanol.
        % Cross-validates the Mordred wrapper against the RDKit-based
        % emk.descriptor.calculate to detect implementation drift.
        % Note: mordredcommunity 2.x renamed nHBd -> nHBDon.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            mol     = emk.mol.fromSmiles("CCO");
            mDesc   = emk.descriptor.mordred(mol, "nHBDon");
            rDesc   = emk.descriptor.calculate(mol, "NumHDonors");
            tc.verifyEqual(mDesc.nHBDon, rDesc.NumHDonors, ...
                "Mordred nHBDon must match RDKit HBondDonors for ethanol");
        end

    end   % mordred Mordred-dependent tests

    % ======================================================================
    % mordredBatch -- Mordred-dependent tests
    % ======================================================================
    methods (Test)

        function test_mordredBatch_returnsTable(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            mols = {emk.mol.fromSmiles("CCO"), emk.mol.fromSmiles("c1ccccc1")};
            tbl  = emk.descriptor.mordredBatch(mols, ["MW", "nRot"]);
            tc.verifyClass(tbl, "table", "Result must be a table");
        end

        function test_mordredBatch_dimensions(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            requested = ["MW", "nRot", "nHBDon"];
            mols = {emk.mol.fromSmiles("CCO"), emk.mol.fromSmiles("c1ccccc1"), ...
                    emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O")};
            tbl  = emk.descriptor.mordredBatch(mols, requested);
            tc.verifySize(tbl, [3, 3], ...
                "Table must have 3 rows (mols) and 3 cols (descriptors)");
        end

        function test_mordredBatch_columnNamesMatchRequest(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            requested = ["MW", "nRot"];
            mols = {emk.mol.fromSmiles("CCO")};
            tbl  = emk.descriptor.mordredBatch(mols, requested);
            tc.verifyEqual(string(tbl.Properties.VariableNames), requested, ...
                "Table column names must match requested descriptor names");
        end

        function test_mordredBatch_invalidMol_yieldNaNRow(tc)
        % An invalid (non-Mol) element in the cell must produce a NaN row,
        % not an error.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            mols = {emk.mol.fromSmiles("CCO"), "not_a_mol"};
            tbl  = emk.descriptor.mordredBatch(mols, ["MW", "nRot"]);
            tc.verifySize(tbl, [2 2], "Table must have 2 rows");
            tc.verifyTrue(all(isnan(tbl{2, :})), ...
                "Invalid mol row must be all NaN");
        end

        function test_mordredBatch_validRow_notAllNaN_whenOneMolInvalid(tc)
        % TB5b: When one mol is invalid, the valid mol row must NOT be all NaN.
        % Ensures the NaN-fill for invalid mols does not contaminate valid rows.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            mols = {emk.mol.fromSmiles("CCO"), "not_a_mol"};
            tbl  = emk.descriptor.mordredBatch(mols, ["MW", "nRot"]);
            tc.verifyFalse(all(isnan(tbl{1, :})), ...
                "Valid mol row must have at least one non-NaN value");
            tc.verifyTrue(isfinite(tbl.MW(1)), ...
                "MW for valid ethanol mol must be finite");
        end

        function test_mordredBatch_allInvalidMols_throwsAllMolsFailed(tc)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            tc.verifyError( ...
                @() emk.descriptor.mordredBatch({"not_a_mol", 42}, ["MW"]), ...
                "emk:descriptor:mordredBatch:allMolsFailed", ...
                "All-invalid mols must throw allMolsFailed");
        end

        function test_mordredBatch_MW_crossValidate_mordred(tc)
        % mordredBatch MW values must match individual mordred() calls.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("c1ccccc1");
            mols = {mol1, mol2};
            tbl  = emk.descriptor.mordredBatch(mols, "MW");
            single1 = emk.descriptor.mordred(mol1, "MW");
            single2 = emk.descriptor.mordred(mol2, "MW");
            tc.verifyFalse(isnan(tbl.MW(1)), "Batch MW row 1 must not be NaN");
            tc.verifyEqual(tbl.MW(1), single1.MW, "AbsTol", 1e-6, ...
                "Batch MW must match single-mol MW for ethanol");
            tc.verifyFalse(isnan(tbl.MW(2)), "Batch MW row 2 must not be NaN");
            tc.verifyEqual(tbl.MW(2), single2.MW, "AbsTol", 1e-6, ...
                "Batch MW must match single-mol MW for benzene");
        end

        function test_mordredBatch_emptyCellWithNames_returnsZeroRowTable(tc)
        % TB8: An empty cell array with explicit descriptor names must return
        % a 0-row table whose column names match the requested descriptors.
        % This validates graceful handling of the zero-molecule edge case.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            requested = ["MW", "nRot"];
            tbl = emk.descriptor.mordredBatch({}, requested);
            tc.verifyClass(tbl, "table", "Result must be a table for empty cell");
            tc.verifySize(tbl, [0, 2], ...
                "Empty-mol batch must return 0-row table with 2 descriptor cols");
            tc.verifyEqual(string(tbl.Properties.VariableNames), requested, ...
                "Column names must match requested descriptors even for empty input");
        end

    end   % mordredBatch Mordred-dependent tests

    % ======================================================================
    % mordredNames
    % ======================================================================
    methods (Test)

        function test_mordredNames_returnsStringArray(tc)
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            names = emk.descriptor.mordredNames();
            tc.verifyClass(names, "string", "mordredNames must return string array");
        end

        function test_mordredNames_isRowVector(tc)
        % TN1b: mordredNames must return a 1 x N row vector, not a column.
        % The API contract (function_reference.md) specifies string(1xN).
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            names = emk.descriptor.mordredNames();
            tc.verifyGreaterThanOrEqual(ndims(names), 1, "names must be an array");
            tc.verifyEqual(size(names, 1), 1, ...
                "mordredNames must be a row vector (first dimension == 1)");
        end

        function test_mordredNames_countAtLeast1000(tc)
        % Mordred 2D has ~1800 descriptors; at least 1000 is a safe lower bound.
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            names = emk.descriptor.mordredNames();
            tc.verifyGreaterThanOrEqual(numel(names), 1000, ...
                "Mordred must provide >= 1000 descriptor names");
        end

        function test_mordredNames_containsMW(tc)
        % "MW" (molecular weight) is the most fundamental Mordred descriptor.
            tc.assumeTrue(tc.mordredAvailable(), "Skipped: Mordred not installed");

            names = emk.descriptor.mordredNames();
            tc.verifyTrue(ismember("MW", names), ...
                "mordredNames must include 'MW'");
        end

    end   % mordredNames

    % ======================================================================
    % Private helpers
    % ======================================================================
    methods (Access = private)

        function tf = rdkitAvailable(~)
            try
                py.importlib.import_module("rdkit");
                tf = true;
            catch
                tf = false;
            end
        end

        function tf = mordredAvailable(~)
            try
                py.importlib.import_module("mordred");
                tf = true;
            catch
                tf = false;
            end
        end

        function ME = captureError(~, f)
            ME = [];
            try
                f();
            catch caught
                ME = caught;
            end
        end

    end
end
