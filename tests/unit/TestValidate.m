classdef TestValidate < matlab.unittest.TestCase
% TestValidate  Unit tests for src/+emk/+setup/validate.m
%
% validate() is a non-throwing diagnostic function that returns a table
% describing the installation status of known libraries.
%
% Coverage:
%   TC1 : Return type is a MATLAB table
%   TC2 : Table has exactly 4 columns (Library, Installed, Version, Track)
%   TC3 : Installed column is logical (not numeric)
%   TC4 : Library column is string (not char or cell)
%   TC5 : Version column is string (not char or cell)
%   TC6 : Track column is string containing "1", "2", or "?"
%   TC7 : Non-string Libraries argument -> invalidInput error
%   TC8 : Default call (no Libraries) includes an "rdkit" row
%   TC9 : Default call returns exactly 7 rows (one per known library)
%   TC10: Custom Libraries returns only the specified rows (row count match)
%   TC11: Custom Libraries row order matches the input order
%   TC12: validate() never throws in any Python state (non-throwing design)
%   TC13: string.empty Libraries arg -> returns 0-row table without error
%
% Run with:
%   addpath(genpath("src")); addpath(genpath("tests"));
%   results = run(TestValidate);

    properties (Access = private)
        DefaultResult  % cached emk.setup.validate() result for default (no-arg) call
    end

    methods (TestClassSetup)
        function cacheDefaultResult(tc)
        % Run validate() once per test class and cache the 7-row table.
        % All test methods that verify the default output reuse this result.
        % validate() is read-only (no state change), so sharing is safe.
            addpath(genpath("src"));
            tc.DefaultResult = emk.setup.validate();
        end
    end

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: Return type
        % ------------------------------------------------------------------

        function test_returnType_isTable(tc)
        % validate() must always return a MATLAB table regardless of
        % Python environment state.
            T = tc.DefaultResult;
            tc.verifyClass(T, "table", "validate() must return a table");
        end

        % ------------------------------------------------------------------
        % TC2: Column names
        % ------------------------------------------------------------------

        function test_returnTable_hasColumn_Library(tc)
        % The returned table must have a 'Library' column.
            T = tc.DefaultResult;
            tc.verifyTrue(ismember("Library", string(T.Properties.VariableNames)), ...
                "Table must have column 'Library'");
        end

        function test_returnTable_hasColumn_Installed(tc)
        % The returned table must have an 'Installed' column.
            T = tc.DefaultResult;
            tc.verifyTrue(ismember("Installed", string(T.Properties.VariableNames)), ...
                "Table must have column 'Installed'");
        end

        function test_returnTable_hasColumn_Version(tc)
        % The returned table must have a 'Version' column.
            T = tc.DefaultResult;
            tc.verifyTrue(ismember("Version", string(T.Properties.VariableNames)), ...
                "Table must have column 'Version'");
        end

        function test_returnTable_hasColumn_Track(tc)
        % The returned table must have a 'Track' column.
            T = tc.DefaultResult;
            tc.verifyTrue(ismember("Track", string(T.Properties.VariableNames)), ...
                "Table must have column 'Track'");
        end

        function test_returnTable_hasExactly4Columns(tc)
        % The table must have exactly 4 columns; no extra or missing columns.
            T = tc.DefaultResult;
            tc.verifyNumElements(T.Properties.VariableNames, 4, ...
                "Table must have exactly 4 columns");
        end

        % ------------------------------------------------------------------
        % TC3-TC6: Column types
        % ------------------------------------------------------------------

        function test_installedColumn_isLogical(tc)
        % Installed must be logical, not double or uint8, to prevent
        % silent numeric-comparison bugs in callers.
            T = tc.DefaultResult;
            tc.verifyClass(T.Installed, "logical", ...
                "Installed column must be logical (not numeric)");
        end

        function test_libraryColumn_isString(tc)
        % Library must be a MATLAB string array, not char or cell.
        % verifyClass catches accidental char/cell return.
            T = tc.DefaultResult;
            tc.verifyClass(T.Library, "string", ...
                "Library column must be string array");
        end

        function test_versionColumn_isString(tc)
        % Version must be a string array (empty string "" for not-installed).
        % verifyClass catches accidental char return from pip output parsing.
            T = tc.DefaultResult;
            tc.verifyClass(T.Version, "string", ...
                "Version column must be string array");
        end

        function test_trackColumn_isString(tc)
        % Track must be a string array containing "1", "2", or "?".
        % verifyClass catches accidental numeric track values.
            T = tc.DefaultResult;
            tc.verifyClass(T.Track, "string", ...
                "Track column must be string array");
        end

        function test_trackColumn_valuesAreKnownTokens(tc)
        % All Track values must be "1", "2", or "?".
        % Ensures resolveLibMeta_() covers every library in the default list.
            T = tc.DefaultResult;
            allowed = ["1", "2", "?"];
            allValid = all(ismember(T.Track, allowed));
            tc.verifyTrue(allValid, ...
                "All Track values must be ""1"", ""2"", or ""?""");
        end

        % ------------------------------------------------------------------
        % TC7: Invalid Libraries argument
        % ------------------------------------------------------------------

        function test_invalidLibraries_nonString_throwsInvalidInput(tc)
        % A non-string Libraries argument must throw invalidInput before
        % any pip or system() call.
            tc.verifyError( ...
                @() emk.setup.validate(Libraries=42), ...
                "emk:setup:validate:invalidInput", ...
                "Numeric Libraries must throw invalidInput");
        end

        function test_invalidLibraries_cellArray_throwsInvalidInput(tc)
        % A cell array Libraries argument must also throw invalidInput.
            tc.verifyError( ...
                @() emk.setup.validate(Libraries={"rdkit"}), ...
                "emk:setup:validate:invalidInput", ...
                "Cell array Libraries must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC8-TC9: Default library list
        % ------------------------------------------------------------------

        function test_defaultList_containsRdkitRow(tc)
        % The default call must include a row for "rdkit" (core library).
            T = tc.DefaultResult;
            tc.verifyTrue(any(T.Library == "rdkit"), ...
                "Default validate() must include a row for 'rdkit'");
        end

        function test_defaultList_returns7Rows(tc)
        % The default list covers exactly 11 known libraries:
        %   rdkit, pubchempy, mordred, biopython, scipy, meeko,
        %   vina, pdbfixer, openbabel, mdanalysis, pymol
            T = tc.DefaultResult;
            tc.verifyNumElements(T.Library, 11, ...
                "Default validate() must return exactly 11 rows");
        end

        function test_defaultList_containsAllKnownLibraries(tc)
        % Verify each of the 11 expected library names appears exactly once.
            T = tc.DefaultResult;
            expected = ["rdkit"; "pubchempy"; "mordred"; "biopython"; ...
                        "scipy"; "meeko"; "vina"; "pdbfixer"; ...
                        "openbabel"; "mdanalysis"; "pymol"];
            for i = 1:numel(expected)
                tc.verifyTrue(any(T.Library == expected(i)), ...
                    "Default list must include library: " + expected(i));
            end
        end

        % ------------------------------------------------------------------
        % TC10-TC11: Custom Libraries argument
        % ------------------------------------------------------------------

        function test_customLibraries_rowCountMatchesInput(tc)
        % When Libraries is specified, the returned table must have exactly
        % the same number of rows as the number of requested libraries.
            libs = ["rdkit", "pubchempy"];
            T    = emk.setup.validate(Libraries=libs);
            tc.verifyNumElements(T.Library, numel(libs), ...
                "Row count must match the number of requested libraries");
        end

        function test_customLibraries_rowOrderMatchesInput(tc)
        % Rows must appear in the same order as the Libraries input.
        % Callers rely on row i corresponding to Libraries(i).
            libs = ["pubchempy", "rdkit", "mordred"];
            T    = emk.setup.validate(Libraries=libs);
            tc.verifyEqual(T.Library, libs(:), ...
                "Library column order must match Libraries input order");
        end

        function test_customLibraries_onlyRequestedLibrariesReturned(tc)
        % Only the requested libraries appear; validate() must not append
        % extra rows from the default list.
            T = emk.setup.validate(Libraries="rdkit");
            tc.verifyNumElements(T.Library, 1, ...
                "Single-library request must return exactly 1 row");
            tc.verifyEqual(T.Library(1), "rdkit", ...
                "The single row must be 'rdkit'");
        end

        % ------------------------------------------------------------------
        % TC12: Non-throwing design
        % ------------------------------------------------------------------

        function test_neverThrows(tc)
        % validate() must never throw regardless of Python state.
        % This is the fundamental non-throwing contract.
            ME = [];
            try
                emk.setup.validate();
            catch e
                ME = e;
            end
            tc.verifyEmpty(ME, ...
                "validate() must never throw in any Python state");
        end

        function test_customLibraries_neverThrows(tc)
        % validate() with a custom list must also be non-throwing.
            ME = [];
            try
                emk.setup.validate(Libraries=["rdkit", "pymol"]);
            catch e
                ME = e;
            end
            tc.verifyEmpty(ME, ...
                "validate() with custom Libraries must never throw");
        end

        % ------------------------------------------------------------------
        % TC13: string.empty Libraries -> treated as "use default list"
        % ------------------------------------------------------------------

        function test_emptyStringArray_usesDefaultList(tc)
        % string.empty passed as Libraries is treated as "not specified"
        % and triggers the default 11-library list.
        % isempty(string.empty) == true, so the isempty guard fires.
        % This is consistent with the function reference: "Default: all known
        % libraries when Libraries is omitted or empty."
            T = emk.setup.validate(Libraries=string.empty);
            tc.verifyClass(T, "table", "Must still return a table");
            tc.verifyNumElements(T.Library, 11, ...
                "string.empty Libraries must trigger the default 11-library list");
        end

    end
end
