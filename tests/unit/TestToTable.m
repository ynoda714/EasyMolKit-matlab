classdef TestToTable < matlab.unittest.TestCase
% TestToTable  Unit tests for src/+emk/+mol/toTable.m
%
% Run with:
%   addpath(genpath("src"));
%   results = run(TestToTable);
%
% Coverage:
%   TC1: Non-cell / non-Mol input => invalidInput error (no RDKit required)
%         + error message contains offending class name
%         + numeric, logical, string-scalar all rejected
%   TC2: Unknown property name => unknownProperty error (no RDKit required)
%         + error message contains the offending property name
%         + partial unknown (mix of valid and invalid) also rejected
%   TC3: Empty Properties argument => emptyProperties error (no RDKit required)
%   TC4: Empty cell array => 0-row table with correct variable names
%         (no RDKit required)
%   TC5: Default properties => all 11 columns (SMILES + 10 descriptors)
%         (requires RDKit)
%   TC6: Single Mol auto-normalised to 1-row table (requires RDKit)
%         + SMILES column is string, descriptor columns are double
%   TC7: Multiple Mols => N-row table (requires RDKit)
%         + height == number of input mols (verifySize)
%   TC8: Subset Properties -- ["SMILES","MolWt","LogP"] => 3-column table
%         (requires RDKit)
%   TC9: "SMILES" only => single-column string table (requires RDKit)
%         + canonical SMILES round-trip vs emk.mol.toSmiles
%   TC10: Descriptor-only (no "SMILES") => double columns, no SMILES col
%         (requires RDKit)
%   TC11: Invalid Mol in cell => SMILES="<invalid>", descriptor=NaN
%         (requires RDKit)
%   TC12: MolWt reference value cross-validation vs emk.descriptor.calculate
%         (requires RDKit)
%   TC13: Properties column order is preserved in output table (requires RDKit)
%   TC14: Output is a table (verifyClass) (requires RDKit)
%
% Tests requiring RDKit use assumeTrue(tc.rdkitAvailable()) to skip
% gracefully when Python/RDKit is not configured in the current session.

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: Input validation -- non-cell / non-Mol input (no RDKit required)
        % ------------------------------------------------------------------

        function test_toTable_numericInput_throwsInvalidInput(tc)
        % Numeric input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.mol.toTable(42), ...
                "emk:mol:toTable:invalidInput", ...
                "Numeric input must throw invalidInput");
        end

        function test_toTable_numericInput_errorMessage_containsClass(tc)
        % Error message must contain the offending class name.
            ME = tc.captureError(@() emk.mol.toTable(42));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain 'double' for numeric input");
        end

        function test_toTable_logicalInput_throwsInvalidInput(tc)
        % logical scalar must throw invalidInput.
            tc.verifyError(@() emk.mol.toTable(true), ...
                "emk:mol:toTable:invalidInput", ...
                "logical input must throw invalidInput");
        end

        function test_toTable_stringScalarInput_throwsInvalidInput(tc)
        % A plain string (not a Mol/cell) must throw invalidInput.
            tc.verifyError(@() emk.mol.toTable("CCO"), ...
                "emk:mol:toTable:invalidInput", ...
                "String scalar input must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC2: Unknown property name (no RDKit required)
        % ------------------------------------------------------------------

        function test_toTable_unknownProperty_throwsUnknownProperty(tc)
        % An unknown property name in a cell must throw unknownProperty.
            tc.verifyError(@() emk.mol.toTable({}, Properties="NotAProperty"), ...
                "emk:mol:toTable:unknownProperty", ...
                "Unknown property must throw unknownProperty");
        end

        function test_toTable_unknownProperty_errorMessage_containsName(tc)
        % Error message must contain the offending property name.
            ME = tc.captureError(@() emk.mol.toTable({}, Properties="BadProp"));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "BadProp", ...
                "Error message must mention the offending property name");
        end

        function test_toTable_partialUnknownProperty_throwsUnknownProperty(tc)
        % Mix of valid and invalid property names must still throw unknownProperty.
            tc.verifyError( ...
                @() emk.mol.toTable({}, Properties=["SMILES","NotAProperty"]), ...
                "emk:mol:toTable:unknownProperty", ...
                "Partial unknown in Properties must throw unknownProperty");
        end

        % ------------------------------------------------------------------
        % TC3: Empty Properties (no RDKit required)
        % ------------------------------------------------------------------

        function test_toTable_emptyProperties_throwsEmptyProperties(tc)
        % An empty Properties array must throw emptyProperties.
            tc.verifyError( ...
                @() emk.mol.toTable({}, Properties=string.empty), ...
                "emk:mol:toTable:emptyProperties", ...
                "Empty Properties must throw emptyProperties");
        end

        % ------------------------------------------------------------------
        % TC4: Empty cell => 0-row table (no RDKit required)
        % ------------------------------------------------------------------

        function test_toTable_emptyCell_returnsEmptyTable(tc)
        % An empty cell must return a 0-row table -- no RDKit call needed.
            tbl = emk.mol.toTable({});
            tc.verifyClass(tbl, "table", ...
                "Empty-cell input must return a table");
            tc.verifyEqual(height(tbl), 0, ...
                "Table from empty cell must have 0 rows");
        end

        function test_toTable_emptyCell_correctVariableNames(tc)
        % Empty-cell table must have all 11 default variable names.
        % (SMILES + 10 descriptor names)
            tbl = emk.mol.toTable({});
            expected = ["SMILES","MolWt","ExactMolWt","LogP","TPSA", ...
                        "NumHAcceptors","NumHDonors","NumRotatableBonds", ...
                        "RingCount","FractionCSP3","HeavyAtomCount"];
            actual = string(tbl.Properties.VariableNames);
            tc.verifyEqual(actual, expected, ...
                "Default variable names must match expected order");
        end

        function test_toTable_emptyCell_smilesColumn_isStringType(tc)
        % The SMILES column of a 0-row table must have type string.
            tbl = emk.mol.toTable({});
            tc.verifyClass(tbl.SMILES, "string", ...
                "SMILES column of empty table must be string type");
        end

        function test_toTable_emptyCell_descriptorColumns_areDoubleType(tc)
        % All descriptor columns of a 0-row table must have type double.
            tbl      = emk.mol.toTable({});
            descCols = ["MolWt","ExactMolWt","LogP","TPSA","NumHAcceptors", ...
                        "NumHDonors","NumRotatableBonds","RingCount", ...
                        "FractionCSP3","HeavyAtomCount"];
            for k = 1:numel(descCols)
                tc.verifyClass(tbl.(descCols(k)), "double", ...
                    descCols(k) + " column of empty table must be double type");
            end
        end

        function test_toTable_emptyCell_customProps_correctNames(tc)
        % Empty-cell table with custom Properties must have matching names.
            tbl = emk.mol.toTable({}, Properties=["SMILES","MolWt"]);
            actual = string(tbl.Properties.VariableNames);
            tc.verifyEqual(actual, ["SMILES","MolWt"], ...
                "Custom Properties must set correct variable names");
        end

        % ------------------------------------------------------------------
        % TC5: Default properties (requires RDKit)
        % ------------------------------------------------------------------

        function test_toTable_defaultProps_has11Columns(tc)
        % Default call must return a table with 11 columns (SMILES + 10 desc).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tbl = emk.mol.toTable({mol});
            tc.verifyEqual(width(tbl), 11, ...
                "Default table must have 11 columns");
        end

        function test_toTable_defaultProps_containsSMILES(tc)
        % Default table must contain a "SMILES" column.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tbl = emk.mol.toTable({mol});
            tc.verifyTrue(ismember("SMILES", string(tbl.Properties.VariableNames)), ...
                "Default table must contain SMILES column");
        end

        % ------------------------------------------------------------------
        % TC6: Single Mol auto-normalised (requires RDKit)
        % ------------------------------------------------------------------

        function test_toTable_singleMol_returns1RowTable(tc)
        % A single (non-cell) Mol must produce a 1-row table.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tbl = emk.mol.toTable(mol);
            tc.verifyClass(tbl, "table", ...
                "Single-Mol input must return a table");
            tc.verifyEqual(height(tbl), 1, ...
                "Single-Mol input must produce a 1-row table");
        end

        function test_toTable_singleMol_smilesColumnIsString(tc)
        % SMILES column must have class "string".
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tbl = emk.mol.toTable(mol, Properties="SMILES");
            tc.verifyClass(tbl.SMILES, "string", ...
                "SMILES column must be a string array");
        end

        function test_toTable_singleMol_molwtColumnIsDouble(tc)
        % MolWt column must have class "double".
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tbl = emk.mol.toTable(mol, Properties="MolWt");
            tc.verifyClass(tbl.MolWt, "double", ...
                "MolWt column must be a double array");
        end

        % ------------------------------------------------------------------
        % TC7: Multiple Mols => N-row table (requires RDKit)
        % ------------------------------------------------------------------

        function test_toTable_threeMols_returns3RowTable(tc)
        % Three Mols must produce a 3-row table.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mols = {emk.mol.fromSmiles("CCO"), ...
                    emk.mol.fromSmiles("c1ccccc1"), ...
                    emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O")};
            tbl = emk.mol.toTable(mols, Properties=["SMILES","MolWt"]);
            tc.verifySize(tbl, [3, 2], ...
                "Three Mols must produce a 3x2 table");
        end

        % ------------------------------------------------------------------
        % TC8: Subset Properties (requires RDKit)
        % ------------------------------------------------------------------

        function test_toTable_subsetProps_correctWidth(tc)
        % ["SMILES","MolWt","LogP"] must produce a 3-column table.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tbl = emk.mol.toTable({mol}, Properties=["SMILES","MolWt","LogP"]);
            tc.verifyEqual(width(tbl), 3, ...
                "Subset Properties must produce correct column count");
        end

        function test_toTable_subsetProps_correctNames(tc)
        % Variable names must match requested Properties exactly.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            req  = ["SMILES","MolWt","LogP"];
            tbl  = emk.mol.toTable({mol}, Properties=req);
            actual = string(tbl.Properties.VariableNames);
            tc.verifyEqual(actual, req, ...
                "VariableNames must match requested Properties");
        end

        % ------------------------------------------------------------------
        % TC9: SMILES-only properties (requires RDKit)
        % ------------------------------------------------------------------

        function test_toTable_smilesOnly_singleColumn(tc)
        % Properties=["SMILES"] must produce a 1-column string table.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tbl = emk.mol.toTable({mol}, Properties="SMILES");
            tc.verifyEqual(width(tbl), 1, ...
                "SMILES-only Properties must produce 1 column");
            tc.verifyClass(tbl.SMILES, "string", ...
                "SMILES-only column must be a string array");
        end

        function test_toTable_smilesOnly_matchesToSmiles(tc)
        % SMILES from toTable must equal emk.mol.toSmiles round-trip.
        % Principle 3: cross-validation between MATLAB and Python API.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol      = emk.mol.fromSmiles("OCC");
            tbl      = emk.mol.toTable({mol}, Properties="SMILES");
            expected = emk.mol.toSmiles(mol);
            tc.verifyEqual(tbl.SMILES(1), expected, ...
                "toTable SMILES must equal emk.mol.toSmiles output");
        end

        function test_toTable_multipleMols_allSmilesMatchToSmiles(tc)
        % Each row's SMILES must individually match emk.mol.toSmiles.
        % Cross-validates all rows, not just the first.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smiles = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
            mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
            tbl    = emk.mol.toTable(mols, Properties="SMILES");
            for i = 1:numel(mols)
                expected = emk.mol.toSmiles(mols{i});
                tc.verifyEqual(tbl.SMILES(i), expected, ...
                    sprintf("Row %d SMILES must match emk.mol.toSmiles", i));
            end
        end

        % ------------------------------------------------------------------
        % TC10: Descriptor-only (no SMILES) (requires RDKit)
        % ------------------------------------------------------------------

        function test_toTable_descriptorOnly_noSMILESColumn(tc)
        % Properties=["MolWt","LogP"] must produce a table without SMILES.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            tbl  = emk.mol.toTable({mol}, Properties=["MolWt","LogP"]);
            actual = string(tbl.Properties.VariableNames);
            tc.verifyFalse(any(actual == "SMILES"), ...
                "Descriptor-only table must not contain SMILES column");
            tc.verifyEqual(width(tbl), 2, ...
                "Descriptor-only Properties must produce 2 columns");
        end

        % ------------------------------------------------------------------
        % TC11: Invalid Mol in cell (requires RDKit)
        % ------------------------------------------------------------------

        function test_toTable_invalidMolInCell_smilesIsInvalidMarker(tc)
        % A non-Mol element in the cell must produce "<invalid>" in SMILES.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            validMol = emk.mol.fromSmiles("CCO");
            mols = {validMol, 42};   % second element is invalid
            tbl  = emk.mol.toTable(mols, Properties=["SMILES","MolWt"]);
            tc.verifyEqual(tbl.SMILES(2), "<invalid>", ...
                "Non-Mol element must produce '<invalid>' in SMILES column");
        end

        function test_toTable_invalidMolInCell_descriptorIsNaN(tc)
        % A non-Mol element in the cell must produce NaN in descriptor columns.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            validMol = emk.mol.fromSmiles("CCO");
            mols = {validMol, 42};
            tbl  = emk.mol.toTable(mols, Properties=["SMILES","MolWt"]);
            tc.verifyTrue(isnan(tbl.MolWt(2)), ...
                "Non-Mol element must produce NaN in descriptor column");
        end

        function test_toTable_invalidMolInCell_validMolUnchanged(tc)
        % The valid mol in the same cell must still produce correct values.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            mols = {mol, 42};
            tbl  = emk.mol.toTable(mols, Properties=["SMILES","MolWt"]);
            expected = emk.mol.toSmiles(mol);
            tc.verifyEqual(tbl.SMILES(1), expected, ...
                "Valid mol must still have correct SMILES when invalid mol is present");
        end

        function test_toTable_allInvalidMols_allSmilesAreInvalidMarker(tc)
        % When every element is non-Mol, every SMILES must be "<invalid>".
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mols = {42, "bad", true};
            tbl  = emk.mol.toTable(mols, Properties="SMILES");
            for i = 1:height(tbl)
                tc.verifyEqual(tbl.SMILES(i), "<invalid>", ...
                    sprintf("Row %d SMILES must be <invalid>", i));
            end
        end

        function test_toTable_allInvalidMols_descriptorRequest_throwsAllMolsFailed(tc)
        % When ALL elements are non-Mol and a descriptor column is requested,
        % batchCalculate has no valid mols to process and throws allMolsFailed.
        % (Contrast with all-invalid + SMILES-only: that path succeeds.)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mols = {42, "bad"};
            tc.verifyError( ...
                @() emk.mol.toTable(mols, Properties=["MolWt","LogP"]), ...
                "emk:descriptor:batchCalculate:allMolsFailed");
        end

        % ------------------------------------------------------------------
        % TC12: MolWt reference cross-validation (requires RDKit)
        % ------------------------------------------------------------------

        function test_toTable_ethanolMolWt_matchesCalculate(tc)
        % toTable MolWt must match emk.descriptor.calculate for the same mol.
        % Principle 3: cross-validation between toTable and calculate.
        % Reference: ethanol C2H6O = 46.069 g/mol (PubChem CID 702)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol      = emk.mol.fromSmiles("CCO");
            tbl      = emk.mol.toTable({mol}, Properties="MolWt");
            expected = emk.descriptor.calculate(mol, "MolWt").MolWt;
            tc.verifyEqual(tbl.MolWt(1), expected, "AbsTol", 1e-9, ...
                "toTable MolWt must exactly match emk.descriptor.calculate");
        end

        function test_toTable_aspirinMolWt_referenceValue(tc)
        % Aspirin (CC(=O)Oc1ccccc1C(=O)O) MW = 180.159 g/mol.
        % AbsTol = 0.01 g/mol (same tolerance as TestDescriptor TC2).
        % Reference: PubChem CID 2244 (Aspirin).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            tbl = emk.mol.toTable({mol}, Properties="MolWt");
            tc.verifyEqual(tbl.MolWt(1), 180.159, "AbsTol", 0.01, ...
                "Aspirin MolWt must be 180.159 +/- 0.01 g/mol");
        end

        % ------------------------------------------------------------------
        % TC13: Column order preserved (requires RDKit)
        % ------------------------------------------------------------------

        function test_toTable_columnOrderPreserved(tc)
        % Output columns must appear in the same order as requested Properties.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            req  = ["LogP","SMILES","MolWt"];
            tbl  = emk.mol.toTable({mol}, Properties=req);
            actual = string(tbl.Properties.VariableNames);
            tc.verifyEqual(actual, req, ...
                "Column order must match Properties request order");
        end

        % ------------------------------------------------------------------
        % TC14: Output type (requires RDKit)
        % ------------------------------------------------------------------

        function test_toTable_outputIsTable(tc)
        % toTable must always return a MATLAB table object.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("c1ccccc1");
            tbl = emk.mol.toTable({mol});
            tc.verifyClass(tbl, "table", ...
                "Return type must be table");
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

    end
end
