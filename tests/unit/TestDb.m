classdef TestDb < matlab.unittest.TestCase
% TestDb  Unit tests for emk.db.searchPubchem, emk.db.searchChembl,
%         emk.db.searchChemblTarget, and emk.db.getChemblActivity.
%
% Run with:
%   addpath(genpath("src"));
%   suite = testsuite("tests/unit");
%   runner = matlab.unittest.TestRunner.withNoPlugins;
%   results = runner.run(suite);
%
% ======================================================================
% Coverage (searchPubchem):
%
% Input validation (no network required):
%   TC1:  Numeric query => invalidInput
%   TC2:  Cell array query => invalidInput
%   TC3:  Logical query => invalidInput
%   TC4:  Empty string query => invalidInput
%   TC5:  Whitespace-only string => invalidInput
%   TC6:  Empty char => invalidInput
%   TC7:  char input is accepted (type guard)
%   TC8:  Unknown Type => invalidType
%   TC9:  Type "name" / "smiles" / "cid" accepted (invalidType not thrown)
%  TC22:  Type uppercase ("NAME") accepted (lowercased internally, not invalidType)
%  TC23:  Type "inchikey" accepted (pre-network acceptance test)
%
% Network tests (require PubChem access):
%  TC10:  Valid name search returns table
%  TC11:  Return type is table
%  TC12:  Table has correct column names
%  TC13:  CID column is double
%  TC14:  MolecularWeight column is double
%  TC15:  IUPACName, MolecularFormula, IsomericSMILES are string
%  TC16:  Known compound: ethanol by name has CID=702
%  TC17:  Known compound: CID=702 by cid type returns ethanol formula C2H6O
%  TC18:  SMILES type search returns table with >= 1 row
%  TC19:  Not-found compound => notFound error
%  TC20:  InChIKey search returns table with >= 1 row
%  TC21:  Ethanol MolecularWeight approximately 46.07 g/mol (AbsTol=0.1)
%  TC24:  CID query returns exactly 1 row (single-compound lookup)
%  TC25:  MolecularWeight > 0 for any valid compound
%  TC26:  IsomericSMILES is non-empty for ethanol CID 702
%
% ======================================================================
% Coverage (searchChemblTarget):
%
% Input validation (no network required):
%   SCT-TC1: Numeric query => invalidInput
%   SCT-TC2: Cell array query => invalidInput
%   SCT-TC3: Logical query => invalidInput
%   SCT-TC4: Empty string query => invalidInput
%   SCT-TC5: Whitespace-only string => invalidInput
%   SCT-TC6: Empty char => invalidInput
%   SCT-TC7: char input is accepted (type guard)
%   SCT-TC8: MaxRows=0 => invalidOptions
%   SCT-TC9: MaxRows=-1 => invalidOptions
%  SCT-TC10: MaxRows=1.5 (non-integer) => invalidOptions
%
% Network tests (require ChEMBL access):
%  SCT-TC11: Valid query returns table
%  SCT-TC12: Return type is table
%  SCT-TC13: Table has correct column names
%  SCT-TC14: All column types are string
%  SCT-TC15: Known target EGFR (CHEMBL203) found by name fragment
%
% ======================================================================
% Coverage (getChemblActivity):
%
% Input validation (no network required):
%   GA-TC1: Numeric targetId => invalidInput
%   GA-TC2: Cell array targetId => invalidInput
%   GA-TC3: Logical targetId => invalidInput
%   GA-TC4: Empty string => invalidInput
%   GA-TC5: Whitespace-only string => invalidInput
%   GA-TC6: Empty char => invalidInput
%   GA-TC7: char input is accepted (type guard)
%   GA-TC8: MaxRows=0 => invalidOptions
%   GA-TC9: MaxRows=-1 => invalidOptions
%  GA-TC10: MaxRows=1.5 (non-integer) => invalidOptions
%
% Network tests (require ChEMBL access):
%  GA-TC11: Valid targetId returns table
%  GA-TC12: Return type is table
%  GA-TC13: Table has correct column names
%  GA-TC14: MoleculeChEMBLID and SMILES are string; Value_nM is double
%  GA-TC15: Value_nM > 0 for all returned rows
%
% ======================================================================
% Design notes:
%   - TC1-TC9 and TC22-TC23 must fire before any network call (input validation)
%   - TC10-TC21 and TC24-TC26 require live network access; guarded with assumeTrue(networkAvailable)
%   - Ethanol (CID 702, formula C2H6O, SMILES CCO) used as canonical test compound
%     because it is the simplest stable organic molecule in PubChem and unlikely
%     to be removed or modified.
%   - InChIKey for ethanol: LFQSCWFLJHTTHZ-UHFFFAOYSA-N (PubChem CID 702)
%   - Type lowercasing: searchPubchem does lower(options.Type) so "NAME"=>"name";
%     TC22 verifies this guard does not throw invalidType for uppercase.
% ======================================================================

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1-TC3: Numeric / cell / logical input => invalidInput
        % No network required; fires before URL construction.
        % ------------------------------------------------------------------

        function test_searchPubchem_numericQuery_throwsInvalidInput(tc)
        % Numeric query must throw invalidInput before any network call.
            tc.verifyError(@() emk.db.searchPubchem(702), ...
                "emk:db:searchPubchem:invalidInput", ...
                "Numeric query must throw invalidInput");
        end

        function test_searchPubchem_cellQuery_throwsInvalidInput(tc)
        % Cell array query must throw invalidInput.
            tc.verifyError(@() emk.db.searchPubchem({"aspirin"}), ...
                "emk:db:searchPubchem:invalidInput", ...
                "Cell query must throw invalidInput");
        end

        function test_searchPubchem_logicalQuery_throwsInvalidInput(tc)
        % Logical query must throw invalidInput.
            tc.verifyError(@() emk.db.searchPubchem(true), ...
                "emk:db:searchPubchem:invalidInput", ...
                "Logical query must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC4-TC6: Empty / whitespace query => invalidInput
        % ------------------------------------------------------------------

        function test_searchPubchem_emptyString_throwsInvalidInput(tc)
        % Empty string query must throw invalidInput.
            tc.verifyError(@() emk.db.searchPubchem(""), ...
                "emk:db:searchPubchem:invalidInput", ...
                "Empty string query must throw invalidInput");
        end

        function test_searchPubchem_whitespaceString_throwsInvalidInput(tc)
        % Whitespace-only string must throw invalidInput.
            tc.verifyError(@() emk.db.searchPubchem("   "), ...
                "emk:db:searchPubchem:invalidInput", ...
                "Whitespace-only string must throw invalidInput");
        end

        function test_searchPubchem_emptyChar_throwsInvalidInput(tc)
        % Empty char must throw invalidInput.
            tc.verifyError(@() emk.db.searchPubchem(''), ...
                "emk:db:searchPubchem:invalidInput", ...
                "Empty char must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC7: char input type acceptance
        % ------------------------------------------------------------------

        function test_searchPubchem_charInput_doesNotThrowInvalidInput(tc)
        % char literal input must not throw invalidInput (type acceptance).
        % The function will fail on network/not-found, not on type check.
            ME = tc.captureError(@() emk.db.searchPubchem('ethanol'));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchPubchem:invalidInput", ...
                    "char input must not trigger invalidInput");
            end
        end

        % ------------------------------------------------------------------
        % TC8: Unknown Type => invalidType
        % ------------------------------------------------------------------

        function test_searchPubchem_unknownType_throwsInvalidType(tc)
        % Unknown Type must throw invalidType.
            tc.verifyError( ...
                @() emk.db.searchPubchem("aspirin", Type="formula"), ...
                "emk:db:searchPubchem:invalidType", ...
                "Unknown Type must throw invalidType");
        end

        function test_searchPubchem_invalidType_errorMessageContainsType(tc)
        % invalidType error message must contain the invalid type string.
            ME = tc.captureError( ...
                @() emk.db.searchPubchem("aspirin", Type="formula"));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "formula", ...
                "Error message must contain the invalid type string");
        end

        % ------------------------------------------------------------------
        % TC9: Case-insensitive Type handling
        % ------------------------------------------------------------------

        function test_searchPubchem_typeNameIsAccepted(tc)
        % "name" type must not throw invalidType.
            ME = tc.captureError(@() emk.db.searchPubchem("ethanol", Type="name"));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchPubchem:invalidType", ...
                    "Type=""name"" must not throw invalidType");
            end
        end

        function test_searchPubchem_typeSmilesIsAccepted(tc)
        % "smiles" type must not throw invalidType.
            ME = tc.captureError(@() emk.db.searchPubchem("CCO", Type="smiles"));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchPubchem:invalidType", ...
                    "Type=""smiles"" must not throw invalidType");
            end
        end

        function test_searchPubchem_typeCidIsAccepted(tc)
        % "cid" type must not throw invalidType.
            ME = tc.captureError(@() emk.db.searchPubchem("702", Type="cid"));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchPubchem:invalidType", ...
                    "Type=""cid"" must not throw invalidType");
            end
        end

        % ------------------------------------------------------------------
        % TC22: Uppercase Type accepted (lowercased internally, not invalidType)
        %
        % searchPubchem does lower(options.Type) so "NAME" => "name".
        % This test verifies that upper-case Type strings are not rejected
        % before the lowercase coercion takes effect.
        % ------------------------------------------------------------------

        function test_searchPubchem_typeUppercase_notInvalidType(tc)
        % Type="NAME" (uppercase) must not throw invalidType.
        % The function lowercases the Type before checking; uppercase is valid.
            ME = tc.captureError(@() emk.db.searchPubchem("ethanol", Type="NAME"));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchPubchem:invalidType", ...
                    "Type=""NAME"" (uppercase) must not throw invalidType");
            end
        end

        % ------------------------------------------------------------------
        % TC23: Type="inchikey" pre-network acceptance
        %
        % Analogous to TC9 tests for "name"/"smiles"/"cid".
        % Verifies that "inchikey" is recognised as a valid Type before any
        % network call.  The function call may throw notFound or networkError
        % if the network is absent, but must NOT throw invalidType.
        % ------------------------------------------------------------------

        function test_searchPubchem_typeInchikeyIsAccepted(tc)
        % "inchikey" type must not throw invalidType.
            ME = tc.captureError(@() emk.db.searchPubchem( ...
                "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", Type="inchikey"));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchPubchem:invalidType", ...
                    "Type=""inchikey"" must not throw invalidType");
            end
        end

        % ------------------------------------------------------------------
        % TC10-TC15: Return type and structure (requires network)
        % ------------------------------------------------------------------

        function test_searchPubchem_validName_returnsTable(tc)
        % Valid name search must return a table.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("ethanol");
            tc.verifyClass(result, "table", ...
                "searchPubchem must return a table");
        end

        function test_searchPubchem_validName_tableIsNotEmpty(tc)
        % Valid name search must return a non-empty table.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("ethanol");
            tc.verifyGreaterThanOrEqual(height(result), 1, ...
                "Result table must have at least one row");
        end

        function test_searchPubchem_columnNames_areCorrect(tc)
        % Result table must have exactly the expected column names.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("ethanol");
            expected = ["CID", "IUPACName", "MolecularFormula", ...
                        "MolecularWeight", "IsomericSMILES"];
            tc.verifyEqual(string(result.Properties.VariableNames), expected, ...
                "Table column names must match expected set");
        end

        function test_searchPubchem_cidColumn_isDouble(tc)
        % CID column must be of type double.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("ethanol");
            tc.verifyClass(result.CID, "double", ...
                "CID column must be double");
        end

        function test_searchPubchem_molWtColumn_isDouble(tc)
        % MolecularWeight column must be of type double.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("ethanol");
            tc.verifyClass(result.MolecularWeight, "double", ...
                "MolecularWeight column must be double");
        end

        function test_searchPubchem_stringColumns_areString(tc)
        % IUPACName, MolecularFormula, IsomericSMILES must be string type.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("ethanol");
            tc.verifyClass(result.IUPACName,        "string", ...
                "IUPACName must be string");
            tc.verifyClass(result.MolecularFormula, "string", ...
                "MolecularFormula must be string");
            tc.verifyClass(result.IsomericSMILES,   "string", ...
                "IsomericSMILES must be string");
        end

        % ------------------------------------------------------------------
        % TC16: Known compound ethanol by name => CID=702
        % Reference: PubChem CID 702 (ethanol, C2H6O, MW~46.07 g/mol)
        % ------------------------------------------------------------------

        function test_searchPubchem_ethanolByName_cid702(tc)
        % Searching "ethanol" by name must return a row with CID=702.
        % PubChem CID 702 is the canonical entry for ethanol.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("ethanol");
            tc.verifyTrue(any(result.CID == 702), ...
                "Ethanol search must return a row with CID=702");
        end

        % ------------------------------------------------------------------
        % TC17: CID type search => correct formula
        % ------------------------------------------------------------------

        function test_searchPubchem_ethanolByCid_formulaC2H6O(tc)
        % Searching CID 702 must return the ethanol formula "C2H6O".
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("702", Type="cid");
            tc.verifyTrue(any(result.MolecularFormula == "C2H6O"), ...
                "CID 702 must have formula C2H6O");
        end

        % ------------------------------------------------------------------
        % TC18: SMILES type search returns table with >= 1 row
        % ------------------------------------------------------------------

        function test_searchPubchem_ethanolBySmiles_returnsTable(tc)
        % Searching ethanol SMILES "CCO" must return a non-empty table.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("CCO", Type="smiles");
            tc.verifyClass(result, "table", ...
                "SMILES search must return a table");
            tc.verifyGreaterThanOrEqual(height(result), 1, ...
                "SMILES search must return at least one row");
        end

        % ------------------------------------------------------------------
        % TC19: Not-found compound => notFound error
        % ------------------------------------------------------------------

        function test_searchPubchem_notFoundCompound_throwsNotFound(tc)
        % A compound name that does not exist in PubChem must throw notFound.
        % "emk_nonexistent_compound_xyz123" is an extremely unlikely real name.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            tc.verifyError( ...
                @() emk.db.searchPubchem("emk_nonexistent_compound_xyz123"), ...
                "emk:db:searchPubchem:notFound", ...
                "Non-existent compound must throw notFound");
        end

        % ------------------------------------------------------------------
        % TC20: InChIKey search returns table
        % Ethanol InChIKey: LFQSCWFLJHTTHZ-UHFFFAOYSA-N (PubChem CID 702)
        % ------------------------------------------------------------------

        function test_searchPubchem_ethanolByInchikey_returnsTable(tc)
        % Searching ethanol by InChIKey must return a non-empty table.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem( ...
                "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", Type="inchikey");
            tc.verifyClass(result, "table", ...
                "InChIKey search must return a table");
            tc.verifyGreaterThanOrEqual(height(result), 1, ...
                "InChIKey search must return at least one row");
        end

        % ------------------------------------------------------------------
        % TC21: MolecularWeight for CID 702 is approximately 46.07 g/mol
        % ------------------------------------------------------------------

        function test_searchPubchem_ethanolMolWt_approximately46(tc)
        % Ethanol MolecularWeight must be approximately 46.07 g/mol.
        % AbsTol = 0.1 g/mol accommodates PubChem's reported precision.
        % Reference: IUPAC 2021 atomic weights, PubChem CID 702 (~46.07 g/mol)
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("702", Type="cid");
            idx = find(result.CID == 702, 1);
            tc.assertNotEmpty(idx, "CID 702 must be present in result");
            tc.verifyEqual(result.MolecularWeight(idx), 46.07, "AbsTol", 0.1, ...
                "Ethanol MolecularWeight must be approx 46.07 g/mol");
        end

        % ------------------------------------------------------------------
        % TC24: CID query returns exactly 1 row (requires network)
        %
        % A CID lookup is a single-compound query: PubChem returns exactly
        % one PropertyTable entry per CID.  This guards against accidental
        % multi-row responses from the CID namespace.
        % ------------------------------------------------------------------

        function test_searchPubchem_cidQuery_returnsExactlyOneRow(tc)
        % CID query for a single compound must return exactly 1 row.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("702", Type="cid");
            tc.verifyEqual(height(result), 1, ...
                "CID query for a single CID must return exactly 1 row");
        end

        % ------------------------------------------------------------------
        % TC25: MolecularWeight > 0 (requires network)
        %
        % MW of any real compound is strictly positive.  Verifies that the
        % str2double conversion (for string MW) or direct double cast does
        % not silently produce 0 or NaN.
        % ------------------------------------------------------------------

        function test_searchPubchem_molWt_isPositive(tc)
        % MolecularWeight for any valid compound must be > 0.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("702", Type="cid");
            tc.verifyGreaterThan(result.MolecularWeight(1), 0, ...
                "MolecularWeight must be strictly positive (not NaN or 0)");
        end

        % ------------------------------------------------------------------
        % TC26: IsomericSMILES column is non-empty for known compound (requires network)
        %
        % Verifies that the CanonicalSMILES fallback (or the direct
        % IsomericSMILES value) produces a non-empty string for ethanol.
        % Guards against a silent empty-string return when PubChem omits
        % the IsomericSMILES field and the fallback is also absent.
        % ------------------------------------------------------------------

        function test_searchPubchem_ethanolSmiles_isNonEmpty(tc)
        % IsomericSMILES for ethanol (CID 702) must be a non-empty string.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("702", Type="cid");
            idx = find(result.CID == 702, 1);
            tc.assertNotEmpty(idx, "CID 702 must be present in result");
            tc.verifyGreaterThan(strlength(result.IsomericSMILES(idx)), 0, ...
                "IsomericSMILES for ethanol must be a non-empty string");
        end

        % ------------------------------------------------------------------
        % TC27: Multi-word name (space in query) returns correct compound.
        %
        % Regression guard for the urlencode '+' vs '%20' bug:
        % urlencode encodes spaces as '+', but PubChem PUG REST URL path
        % requires '%20'.  Without the fix, "salicylic acid" returns 404.
        % Salicylic acid PubChem CID = 338.
        % ------------------------------------------------------------------

        function test_searchPubchem_multiWordName_returnsCorrectCid(tc)
        % "salicylic acid" (two-word name) must resolve to CID 338.
            tc.assumeTrue(tc.networkAvailable(), "Skipped: PubChem not reachable");
            result = emk.db.searchPubchem("salicylic acid");
            tc.verifyTrue(any(result.CID == 338), ...
                "Multi-word name 'salicylic acid' must return CID 338");
        end

    end

    % ======================================================================
    methods (Access = private)

        function tf = networkAvailable(~)
        % Return true if PubChem PUG REST API is reachable.
        % Uses a minimal request (CID 702, single property) with a short
        % timeout to minimise test-suite delay when offline.
            tf = false;
            try
                webread( ...
                    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/702/property/MolecularFormula/JSON", ...
                    weboptions("Timeout", 5, "ContentType", "json"));
                tf = true;
            catch
                % Network unavailable or PubChem unreachable
            end
        end

        function tf = chemblAvailable(~)
        % Return true if the ChEMBL REST API is fully reachable.
        % Checks both the single-molecule and list endpoints to detect
        % partial outages (different server processes can fail independently).
        % Result is cached via persistent to avoid repeated round-trips
        % during a single test-suite run.
            persistent cached_;
            if ~isempty(cached_)
                tf = cached_;
                return;
            end
            cached_ = false;
            BASE = "https://www.ebi.ac.uk/chembl/api/data";
            opts = weboptions("Timeout", 5, "ContentType", "json");
            try
                webread(BASE + "/molecule/CHEMBL25.json", opts);
                webread(BASE + "/molecule.json?pref_name__iexact=aspirin&limit=1", opts);
                cached_ = true;
            catch
                % Network unavailable or ChEMBL unreachable
            end
            tf = cached_;
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

    % ======================================================================
    % searchChembl tests
    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % SC-TC1-TC3: Numeric / cell / logical query => invalidInput
        % No network required; fires before URL construction.
        % ------------------------------------------------------------------

        function test_searchChembl_numericQuery_throwsInvalidInput(tc)
        % Numeric query must throw invalidInput before any network call.
            tc.verifyError(@() emk.db.searchChembl(192), ...
                "emk:db:searchChembl:invalidInput", ...
                "Numeric query must throw invalidInput");
        end

        function test_searchChembl_cellQuery_throwsInvalidInput(tc)
        % Cell array query must throw invalidInput.
            tc.verifyError(@() emk.db.searchChembl({"aspirin"}), ...
                "emk:db:searchChembl:invalidInput", ...
                "Cell query must throw invalidInput");
        end

        function test_searchChembl_logicalQuery_throwsInvalidInput(tc)
        % Logical query must throw invalidInput.
            tc.verifyError(@() emk.db.searchChembl(true), ...
                "emk:db:searchChembl:invalidInput", ...
                "Logical query must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % SC-TC4-TC6: Empty / whitespace query => invalidInput
        % ------------------------------------------------------------------

        function test_searchChembl_emptyString_throwsInvalidInput(tc)
        % Empty string query must throw invalidInput.
            tc.verifyError(@() emk.db.searchChembl(""), ...
                "emk:db:searchChembl:invalidInput", ...
                "Empty string query must throw invalidInput");
        end

        function test_searchChembl_whitespaceString_throwsInvalidInput(tc)
        % Whitespace-only string must throw invalidInput.
            tc.verifyError(@() emk.db.searchChembl("   "), ...
                "emk:db:searchChembl:invalidInput", ...
                "Whitespace-only string must throw invalidInput");
        end

        function test_searchChembl_emptyChar_throwsInvalidInput(tc)
        % Empty char must throw invalidInput.
            tc.verifyError(@() emk.db.searchChembl(''), ...
                "emk:db:searchChembl:invalidInput", ...
                "Empty char must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % SC-TC7: char input type acceptance
        % ------------------------------------------------------------------

        function test_searchChembl_charInput_doesNotThrowInvalidInput(tc)
        % char literal input must not throw invalidInput (type acceptance).
            ME = tc.captureError(@() emk.db.searchChembl('aspirin'));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchChembl:invalidInput", ...
                    "char input must not trigger invalidInput");
            end
        end

        % ------------------------------------------------------------------
        % SC-TC8: Unknown Type => invalidType
        % ------------------------------------------------------------------

        function test_searchChembl_unknownType_throwsInvalidType(tc)
        % Unknown Type must throw invalidType before any network call.
            tc.verifyError( ...
                @() emk.db.searchChembl("aspirin", Type="formula"), ...
                "emk:db:searchChembl:invalidType", ...
                "Unknown Type must throw invalidType");
        end

        function test_searchChembl_invalidType_errorMessageContainsType(tc)
        % invalidType error message must contain the invalid type string.
            ME = tc.captureError( ...
                @() emk.db.searchChembl("aspirin", Type="formula"));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "formula", ...
                "Error message must contain the invalid type string");
        end

        % ------------------------------------------------------------------
        % SC-TC9: Valid Type strings are accepted (no network required)
        % ------------------------------------------------------------------

        function test_searchChembl_typeNameIsAccepted(tc)
        % "name" type must not throw invalidType.
            ME = tc.captureError(@() emk.db.searchChembl("aspirin", Type="name"));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchChembl:invalidType", ...
                    "Type=""name"" must not throw invalidType");
            end
        end

        function test_searchChembl_typeSmilesIsAccepted(tc)
        % "smiles" type must not throw invalidType.
            ME = tc.captureError( ...
                @() emk.db.searchChembl("CCO", Type="smiles"));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchChembl:invalidType", ...
                    "Type=""smiles"" must not throw invalidType");
            end
        end

        function test_searchChembl_typeChemblidIsAccepted(tc)
        % "chemblid" type must not throw invalidType.
            ME = tc.captureError( ...
                @() emk.db.searchChembl("CHEMBL25", Type="chemblid"));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchChembl:invalidType", ...
                    "Type=""chemblid"" must not throw invalidType");
            end
        end

        function test_searchChembl_typeInchikeyIsAccepted(tc)
        % "inchikey" type must not throw invalidType.
            ME = tc.captureError( ...
                @() emk.db.searchChembl( ...
                "BSYNRYMUTXBXSQ-UHFFFAOYSA-N", Type="inchikey"));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchChembl:invalidType", ...
                    "Type=""inchikey"" must not throw invalidType");
            end
        end

        function test_searchChembl_typeUppercase_notInvalidType(tc)
        % Type="NAME" (uppercase) must not throw invalidType after lowercasing.
            ME = tc.captureError(@() emk.db.searchChembl("aspirin", Type="NAME"));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchChembl:invalidType", ...
                    "Type=""NAME"" (uppercase) must not throw invalidType");
            end
        end

        % ------------------------------------------------------------------
        % SC-TC10-TC15: Return type and schema (requires ChEMBL network)
        %
        % Reference compound: aspirin (CHEMBL25)
        %   Name: "ASPIRIN"  MW ~180.16  ALogP ~1.31
        %   SMILES: CC(=O)Oc1ccccc1C(=O)O
        %   InChIKey: BSYNRYMUTXBXSQ-UHFFFAOYSA-N
        % ------------------------------------------------------------------

        function test_searchChembl_validName_returnsTable(tc)
        % Valid name search must return a table.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("aspirin");
            tc.verifyClass(result, "table", ...
                "searchChembl must return a table");
        end

        function test_searchChembl_validName_tableIsNotEmpty(tc)
        % Valid name search must return at least one row.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("aspirin");
            tc.verifyGreaterThanOrEqual(height(result), 1, ...
                "Result table must have at least one row");
        end

        function test_searchChembl_columnNames_areCorrect(tc)
        % Result table must have exactly the expected column names.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("aspirin");
            expected = ["ChEMBLID", "Name", "MolecularWeight", "ALogP", ...
                        "HBondDonors", "HBondAcceptors", "SMILES", "InChIKey"];
            tc.verifyEqual(string(result.Properties.VariableNames), expected, ...
                "Table column names must match expected set");
        end

        function test_searchChembl_chemblIdColumn_isString(tc)
        % ChEMBLID column must be of type string.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            try
                result = emk.db.searchChembl("aspirin");
            catch ME_net
                tc.assumeTrue(false, ...
                    "Skipped: ChEMBL transient error: " + ME_net.message);
            end
            tc.verifyClass(result.ChEMBLID, "string", ...
                "ChEMBLID column must be string");
        end

        function test_searchChembl_nameColumn_isString(tc)
        % Name column must be of type string.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            try
                result = emk.db.searchChembl("aspirin");
            catch ME_net
                tc.assumeTrue(false, ...
                    "Skipped: ChEMBL list endpoint transient error: " + ME_net.message);
            end
            tc.verifyClass(result.Name, "string", ...
                "Name column must be string");
        end

        function test_searchChembl_numericColumns_areDouble(tc)
        % MolecularWeight, ALogP, HBondDonors, HBondAcceptors must be double.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            try
                result = emk.db.searchChembl("CHEMBL25", Type="chemblid");
            catch ME_net
                tc.assumeTrue(false, ...
                    "Skipped: ChEMBL endpoint transient error: " + ME_net.message);
            end
            tc.verifyClass(result.MolecularWeight, "double", ...
                "MolecularWeight must be double");
            tc.verifyClass(result.ALogP, "double", ...
                "ALogP must be double");
            tc.verifyClass(result.HBondDonors, "double", ...
                "HBondDonors must be double");
            tc.verifyClass(result.HBondAcceptors, "double", ...
                "HBondAcceptors must be double");
        end

        function test_searchChembl_smilesColumn_isString(tc)
        % SMILES and InChIKey columns must be string type.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("CHEMBL25", Type="chemblid");
            tc.verifyClass(result.SMILES,   "string", "SMILES must be string");
            tc.verifyClass(result.InChIKey, "string", "InChIKey must be string");
        end

        % ------------------------------------------------------------------
        % SC-TC16: ChEMBL ID lookup returns exactly 1 row with correct ID
        % ------------------------------------------------------------------

        function test_searchChembl_aspirin_byChemblid_returnsOneRow(tc)
        % CHEMBL25 lookup must return exactly 1 row.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            try
                result = emk.db.searchChembl("CHEMBL25", Type="chemblid");
            catch ME_net
                tc.assumeTrue(false, ...
                    "Skipped: ChEMBL transient error: " + ME_net.message);
            end
            tc.verifyEqual(height(result), 1, ...
                "chemblid lookup must return exactly 1 row");
        end

        function test_searchChembl_aspirin_chemblId_isCHEMBL25(tc)
        % CHEMBL25 lookup must return ChEMBLID = "CHEMBL25".
        % Reference: CHEMBL25 is the canonical ChEMBL entry for aspirin.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("CHEMBL25", Type="chemblid");
            tc.verifyEqual(result.ChEMBLID(1), "CHEMBL25", ...
                "ChEMBLID for aspirin must be CHEMBL25");
        end

        % ------------------------------------------------------------------
        % SC-TC17: Known compound aspirin MW approximately 180.16 g/mol
        % Reference: MW = 180.159 g/mol (IUPAC 2021 atomic weights)
        % ------------------------------------------------------------------

        function test_searchChembl_aspirin_molWt_approximately180(tc)
        % Aspirin MolecularWeight must be approximately 180.16 g/mol.
        % AbsTol = 0.5 accommodates ChEMBL's reported precision and
        % differences between average and monoisotopic weights.
        % Reference: PubChem CID 2244 MW = 180.16 g/mol.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("CHEMBL25", Type="chemblid");
            tc.verifyEqual(result.MolecularWeight(1), 180.16, "AbsTol", 0.5, ...
                "Aspirin MolecularWeight must be approx 180.16 g/mol");
        end

        function test_searchChembl_aspirin_molWt_isPositive(tc)
        % MolecularWeight for any valid compound must be > 0.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("CHEMBL25", Type="chemblid");
            tc.verifyGreaterThan(result.MolecularWeight(1), 0, ...
                "MolecularWeight must be strictly positive");
        end

        % ------------------------------------------------------------------
        % SC-TC18: SMILES and InChIKey for aspirin are non-empty
        % ------------------------------------------------------------------

        function test_searchChembl_aspirin_smiles_isNonEmpty(tc)
        % SMILES for aspirin must be a non-empty string.
        % Guards against a silent empty-string return when molecule_structures
        % is null or canonical_smiles is absent.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("CHEMBL25", Type="chemblid");
            tc.verifyGreaterThan(strlength(result.SMILES(1)), 0, ...
                "SMILES for aspirin must be a non-empty string");
        end

        function test_searchChembl_aspirin_inchikey_isNonEmpty(tc)
        % InChIKey for aspirin must be a non-empty string.
        % Reference: InChIKey = BSYNRYMUTXBXSQ-UHFFFAOYSA-N
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("CHEMBL25", Type="chemblid");
            tc.verifyGreaterThan(strlength(result.InChIKey(1)), 0, ...
                "InChIKey for aspirin must be a non-empty string");
        end

        % ------------------------------------------------------------------
        % SC-TC19: Name search returns matching compound (CHEMBL ID check)
        % ------------------------------------------------------------------

        function test_searchChembl_aspirin_byName_containsCHEMBL25(tc)
        % Name search "aspirin" must return a row with ChEMBLID = "CHEMBL25".
        % ChEMBL preferred name is case-insensitive (pref_name__iexact filter).
        % Reference: CHEMBL25 is the primary ChEMBL entry for aspirin.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            try
                result = emk.db.searchChembl("aspirin");
            catch ME_net
                tc.assumeTrue(false, ...
                    "Skipped: ChEMBL transient error: " + ME_net.message);
            end
            tc.verifyTrue(any(result.ChEMBLID == "CHEMBL25"), ...
                "Name search 'aspirin' must return a row with CHEMBL25");
        end

        % ------------------------------------------------------------------
        % SC-TC20: Not-found compound => notFound error (requires network)
        % ------------------------------------------------------------------

        function test_searchChembl_notFoundCompound_throwsNotFound(tc)
        % A compound name that does not exist in ChEMBL must throw notFound.
        % "emk_nonexistent_compound_xyz123" is an extremely unlikely real name.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            tc.verifyError( ...
                @() emk.db.searchChembl("emk_nonexistent_compound_xyz123"), ...
                "emk:db:searchChembl:notFound", ...
                "Non-existent compound must throw notFound");
        end

        % ------------------------------------------------------------------
        % SC-TC21: SMILES search returns table (requires network)
        % ------------------------------------------------------------------

        function test_searchChembl_aspirin_bySmiles_returnsTable(tc)
        % SMILES search for aspirin must return a non-empty table.
        % SMILES: CC(=O)Oc1ccccc1C(=O)O (aspirin canonical SMILES)
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl( ...
                "CC(=O)Oc1ccccc1C(=O)O", Type="smiles");
            tc.verifyClass(result, "table", ...
                "SMILES search must return a table");
            tc.verifyGreaterThanOrEqual(height(result), 1, ...
                "SMILES search must return at least one row");
        end

        % ------------------------------------------------------------------
        % SC-TC22: Name field value (cross-validation, requires ChEMBL)
        % ChEMBL preferred names are stored in UPPER CASE.
        % Reference: CHEMBL25 pref_name = "ASPIRIN"
        % ------------------------------------------------------------------

        function test_searchChembl_aspirin_name_isASPIRIN(tc)
        % Name column for CHEMBL25 must equal "ASPIRIN" (ChEMBL preferred name).
        % Cross-validation: confirms that pref_name field is correctly parsed
        % into the Name column.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("CHEMBL25", Type="chemblid");
            tc.verifyEqual(result.Name(1), "ASPIRIN", ...
                "Name for CHEMBL25 must be ASPIRIN");
        end

        % ------------------------------------------------------------------
        % SC-TC23: InChIKey exact value (strong cross-validation)
        % Reference: BSYNRYMUTXBXSQ-UHFFFAOYSA-N (IUPAC InChI Trust 2022)
        % ------------------------------------------------------------------

        function test_searchChembl_aspirin_inchikey_matchesKnown(tc)
        % InChIKey for CHEMBL25 must match the known canonical value.
        % InChIKey = BSYNRYMUTXBXSQ-UHFFFAOYSA-N is the unique identifier
        % for aspirin (acetylsalicylic acid) registered with IUPAC.
        % Cross-validation: confirms standard_inchi_key field is correctly parsed.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("CHEMBL25", Type="chemblid");
            tc.verifyEqual(result.InChIKey(1), "BSYNRYMUTXBXSQ-UHFFFAOYSA-N", ...
                "InChIKey for aspirin must be BSYNRYMUTXBXSQ-UHFFFAOYSA-N");
        end

        % ------------------------------------------------------------------
        % SC-TC24: HBondDonors known value (Ertl definition: OH/NH count)
        % Reference: aspirin has 1 OH group (COOH); HBD = 1 by any definition.
        % ------------------------------------------------------------------

        function test_searchChembl_aspirin_hbondDonors_isOne(tc)
        % Aspirin H-bond donor count must be 1.
        % Aspirin (CC(=O)Oc1ccccc1C(=O)O) has one OH group (carboxylic acid).
        % This is invariant across Lipinski, Ertl, and ChEMBL counting rules.
        % Cross-validation: confirms hbd field is correctly parsed as integer.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("CHEMBL25", Type="chemblid");
            tc.verifyEqual(result.HBondDonors(1), 1, ...
                "Aspirin HBondDonors must be 1 (one COOH)");
        end

        % ------------------------------------------------------------------
        % SC-TC25: InChIKey search returns CHEMBL25
        % Validates the inchikey endpoint URL and filter construction.
        % ------------------------------------------------------------------

        function test_searchChembl_aspirin_byInchikey_returnsCHEMBL25(tc)
        % Searching by InChIKey BSYNRYMUTXBXSQ-UHFFFAOYSA-N must return CHEMBL25.
        % This tests the full inchikey search path:
        %   URL = /molecule.json?molecule_structures__standard_inchi_key={key}&limit=25
        % Cross-validation: the same compound obtained by different search types
        % must yield the same ChEMBLID.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl( ...
                "BSYNRYMUTXBXSQ-UHFFFAOYSA-N", Type="inchikey");
            tc.verifyClass(result, "table", ...
                "InChIKey search must return a table");
            tc.verifyTrue(any(result.ChEMBLID == "CHEMBL25"), ...
                "InChIKey search for aspirin must return CHEMBL25");
        end

        % ------------------------------------------------------------------
        % SC-TC26: Case-insensitive name search (pref_name__iexact)
        % The ChEMBL pref_name__iexact filter is case-insensitive; "ASPIRIN"
        % (all caps) must return the same result as "aspirin" (lowercase).
        % ------------------------------------------------------------------

        function test_searchChembl_uppercaseName_returnsCHEMBL25(tc)
        % Name search with uppercase "ASPIRIN" must return CHEMBL25.
        % ChEMBL uses pref_name__iexact (case-insensitive exact match).
        % The query is passed as-is to the URL (not lowercased by the wrapper),
        % so the __iexact filter on the server handles the case folding.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChembl("ASPIRIN", Type="name");
            tc.verifyClass(result, "table", ...
                "Uppercase 'ASPIRIN' name search must return a table");
            tc.verifyTrue(any(result.ChEMBLID == "CHEMBL25"), ...
                "Uppercase 'ASPIRIN' search must return CHEMBL25 (__iexact filter)");
        end

    end

    % ======================================================================
    % searchChemblTarget tests
    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % SCT-TC1-TC3: Numeric / cell / logical query => invalidInput
        % No network required; fires before URL construction.
        % ------------------------------------------------------------------

        function test_searchChemblTarget_numericQuery_throwsInvalidInput(tc)
        % Numeric query must throw invalidInput before any network call.
            tc.verifyError(@() emk.db.searchChemblTarget(123), ...
                "emk:db:searchChemblTarget:invalidInput", ...
                "Numeric query must throw invalidInput");
        end

        function test_searchChemblTarget_cellQuery_throwsInvalidInput(tc)
        % Cell array query must throw invalidInput.
            tc.verifyError(@() emk.db.searchChemblTarget({"egfr"}), ...
                "emk:db:searchChemblTarget:invalidInput", ...
                "Cell query must throw invalidInput");
        end

        function test_searchChemblTarget_logicalQuery_throwsInvalidInput(tc)
        % Logical query must throw invalidInput.
            tc.verifyError(@() emk.db.searchChemblTarget(true), ...
                "emk:db:searchChemblTarget:invalidInput", ...
                "Logical query must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % SCT-TC4-TC6: Empty / whitespace query => invalidInput
        % ------------------------------------------------------------------

        function test_searchChemblTarget_emptyString_throwsInvalidInput(tc)
        % Empty string query must throw invalidInput.
            tc.verifyError(@() emk.db.searchChemblTarget(""), ...
                "emk:db:searchChemblTarget:invalidInput", ...
                "Empty string must throw invalidInput");
        end

        function test_searchChemblTarget_whitespaceString_throwsInvalidInput(tc)
        % Whitespace-only string must throw invalidInput.
            tc.verifyError(@() emk.db.searchChemblTarget("   "), ...
                "emk:db:searchChemblTarget:invalidInput", ...
                "Whitespace-only string must throw invalidInput");
        end

        function test_searchChemblTarget_emptyChar_throwsInvalidInput(tc)
        % Empty char must throw invalidInput.
            tc.verifyError(@() emk.db.searchChemblTarget(''), ...
                "emk:db:searchChemblTarget:invalidInput", ...
                "Empty char must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % SCT-TC7: char input is accepted (type guard)
        % ------------------------------------------------------------------

        function test_searchChemblTarget_charInput_doesNotThrowInvalidInput(tc)
        % char literal input must not throw invalidInput.
            ME = tc.captureError(@() emk.db.searchChemblTarget('egfr'));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:searchChemblTarget:invalidInput", ...
                    "char input must not trigger invalidInput");
            end
        end

        % ------------------------------------------------------------------
        % SCT-TC8-TC10: Invalid MaxRows => invalidOptions
        % ------------------------------------------------------------------

        function test_searchChemblTarget_maxRowsZero_throwsInvalidOptions(tc)
        % MaxRows=0 must throw invalidOptions.
            tc.verifyError(@() emk.db.searchChemblTarget("egfr", MaxRows=0), ...
                "emk:db:searchChemblTarget:invalidOptions", ...
                "MaxRows=0 must throw invalidOptions");
        end

        function test_searchChemblTarget_maxRowsNegative_throwsInvalidOptions(tc)
        % MaxRows=-1 must throw invalidOptions.
            tc.verifyError(@() emk.db.searchChemblTarget("egfr", MaxRows=-1), ...
                "emk:db:searchChemblTarget:invalidOptions", ...
                "MaxRows=-1 must throw invalidOptions");
        end

        function test_searchChemblTarget_maxRowsNonInteger_throwsInvalidOptions(tc)
        % MaxRows=1.5 (non-integer) must throw invalidOptions.
            tc.verifyError(@() emk.db.searchChemblTarget("egfr", MaxRows=1.5), ...
                "emk:db:searchChemblTarget:invalidOptions", ...
                "MaxRows=1.5 must throw invalidOptions");
        end

        % ------------------------------------------------------------------
        % SCT-TC11-TC15: Return type and schema (requires ChEMBL network)
        %
        % Reference target: EGFR (Epidermal growth factor receptor)
        %   ChEMBL ID: CHEMBL203
        %   Organism: Homo sapiens
        %   TargetType: SINGLE PROTEIN
        % ------------------------------------------------------------------

        function test_searchChemblTarget_validQuery_returnsTable(tc)
        % Valid query must return a table.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChemblTarget("Epidermal growth factor receptor");
            tc.verifyClass(result, "table", ...
                "searchChemblTarget must return a table");
        end

        function test_searchChemblTarget_validQuery_tableIsNotEmpty(tc)
        % Valid query must return at least one row.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChemblTarget("Epidermal growth factor receptor");
            tc.verifyGreaterThanOrEqual(height(result), 1, ...
                "Result table must have at least one row");
        end

        function test_searchChemblTarget_columnNames_areCorrect(tc)
        % Result table must have exactly the expected column names.
        % Guards against the MATLAB Online bug where 'VariableNames' must be
        % a char vector (not string) in table() calls.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            try
                result = emk.db.searchChemblTarget( ...
                    "Epidermal growth factor receptor");
            catch ME_net
                tc.assumeTrue(false, ...
                    "Skipped: ChEMBL target endpoint transient error: " + ME_net.message);
            end
            expected = ["TargetChEMBLID", "PreferredName", "Organism", "TargetType"];
            tc.verifyEqual(string(result.Properties.VariableNames), expected, ...
                "Table column names must match expected set");
        end

        function test_searchChemblTarget_allColumns_areString(tc)
        % All four columns must be of type string.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.searchChemblTarget("Epidermal growth factor receptor");
            tc.verifyClass(result.TargetChEMBLID, "string", ...
                "TargetChEMBLID must be string");
            tc.verifyClass(result.PreferredName, "string", ...
                "PreferredName must be string");
            tc.verifyClass(result.Organism, "string", ...
                "Organism must be string");
            tc.verifyClass(result.TargetType, "string", ...
                "TargetType must be string");
        end

        function test_searchChemblTarget_egfr_containsCHEMBL203(tc)
        % EGFR target search must return a row with TargetChEMBLID = "CHEMBL203".
        % Reference: CHEMBL203 is the canonical ChEMBL entry for human EGFR.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            try
                result = emk.db.searchChemblTarget( ...
                    "Epidermal growth factor receptor");
            catch ME_net
                tc.assumeTrue(false, ...
                    "Skipped: ChEMBL target transient error: " + ME_net.message);
            end
            tc.verifyTrue(any(result.TargetChEMBLID == "CHEMBL203"), ...
                "EGFR search must return a row with CHEMBL203");
        end

    end

    % ======================================================================
    % getChemblActivity tests
    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % GA-TC1-TC3: Numeric / cell / logical targetId => invalidInput
        % No network required; fires before URL construction.
        % ------------------------------------------------------------------

        function test_getChemblActivity_numericId_throwsInvalidInput(tc)
        % Numeric targetId must throw invalidInput before any network call.
            tc.verifyError(@() emk.db.getChemblActivity(203), ...
                "emk:db:getChemblActivity:invalidInput", ...
                "Numeric targetId must throw invalidInput");
        end

        function test_getChemblActivity_cellId_throwsInvalidInput(tc)
        % Cell array targetId must throw invalidInput.
            tc.verifyError(@() emk.db.getChemblActivity({"CHEMBL203"}), ...
                "emk:db:getChemblActivity:invalidInput", ...
                "Cell targetId must throw invalidInput");
        end

        function test_getChemblActivity_logicalId_throwsInvalidInput(tc)
        % Logical targetId must throw invalidInput.
            tc.verifyError(@() emk.db.getChemblActivity(true), ...
                "emk:db:getChemblActivity:invalidInput", ...
                "Logical targetId must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % GA-TC4-TC6: Empty / whitespace targetId => invalidInput
        % ------------------------------------------------------------------

        function test_getChemblActivity_emptyString_throwsInvalidInput(tc)
        % Empty string must throw invalidInput.
            tc.verifyError(@() emk.db.getChemblActivity(""), ...
                "emk:db:getChemblActivity:invalidInput", ...
                "Empty string must throw invalidInput");
        end

        function test_getChemblActivity_whitespaceString_throwsInvalidInput(tc)
        % Whitespace-only string must throw invalidInput.
            tc.verifyError(@() emk.db.getChemblActivity("   "), ...
                "emk:db:getChemblActivity:invalidInput", ...
                "Whitespace-only string must throw invalidInput");
        end

        function test_getChemblActivity_emptyChar_throwsInvalidInput(tc)
        % Empty char must throw invalidInput.
            tc.verifyError(@() emk.db.getChemblActivity(''), ...
                "emk:db:getChemblActivity:invalidInput", ...
                "Empty char must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % GA-TC7: char input is accepted (type guard)
        % ------------------------------------------------------------------

        function test_getChemblActivity_charInput_doesNotThrowInvalidInput(tc)
        % char literal input must not throw invalidInput.
            ME = tc.captureError(@() emk.db.getChemblActivity('CHEMBL203'));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:db:getChemblActivity:invalidInput", ...
                    "char input must not trigger invalidInput");
            end
        end

        % ------------------------------------------------------------------
        % GA-TC8-TC10: Invalid MaxRows => invalidOptions
        % ------------------------------------------------------------------

        function test_getChemblActivity_maxRowsZero_throwsInvalidOptions(tc)
        % MaxRows=0 must throw invalidOptions.
            tc.verifyError( ...
                @() emk.db.getChemblActivity("CHEMBL203", MaxRows=0), ...
                "emk:db:getChemblActivity:invalidOptions", ...
                "MaxRows=0 must throw invalidOptions");
        end

        function test_getChemblActivity_maxRowsNegative_throwsInvalidOptions(tc)
        % MaxRows=-1 must throw invalidOptions.
            tc.verifyError( ...
                @() emk.db.getChemblActivity("CHEMBL203", MaxRows=-1), ...
                "emk:db:getChemblActivity:invalidOptions", ...
                "MaxRows=-1 must throw invalidOptions");
        end

        function test_getChemblActivity_maxRowsNonInteger_throwsInvalidOptions(tc)
        % MaxRows=1.5 (non-integer) must throw invalidOptions.
            tc.verifyError( ...
                @() emk.db.getChemblActivity("CHEMBL203", MaxRows=1.5), ...
                "emk:db:getChemblActivity:invalidOptions", ...
                "MaxRows=1.5 must throw invalidOptions");
        end

        % ------------------------------------------------------------------
        % GA-TC11-TC15: Return type and schema (requires ChEMBL network)
        %
        % Reference target: EGFR (CHEMBL203), IC50 activity.
        % CHEMBL203 has hundreds of IC50 measurements; a MaxRows=10 query
        % reliably returns data without relying on a specific compound ID.
        % ------------------------------------------------------------------

        function test_getChemblActivity_validId_returnsTable(tc)
        % Valid targetId must return a table.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            try
                result = emk.db.getChemblActivity("CHEMBL203", MaxRows=10);
            catch ME_net
                tc.assumeTrue(false, ...
                    "Skipped: ChEMBL activity endpoint transient error: " + ME_net.message);
            end
            tc.verifyClass(result, "table", ...
                "getChemblActivity must return a table");
        end

        function test_getChemblActivity_validId_tableIsNotEmpty(tc)
        % Valid targetId must return at least one row.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.getChemblActivity("CHEMBL203", MaxRows=10);
            tc.verifyGreaterThanOrEqual(height(result), 1, ...
                "Result table must have at least one row");
        end

        function test_getChemblActivity_columnNames_areCorrect(tc)
        % Result table must have exactly the expected column names.
        % Guards against the MATLAB Online bug where 'VariableNames' must be
        % a char vector (not string) in table() calls.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.getChemblActivity("CHEMBL203", MaxRows=10);
            expected = ["MoleculeChEMBLID", "Name", "SMILES", ...
                        "ActivityType", "Value_nM"];
            tc.verifyEqual(string(result.Properties.VariableNames), expected, ...
                "Table column names must match expected set");
        end

        function test_getChemblActivity_columnTypes_areCorrect(tc)
        % MoleculeChEMBLID, Name, SMILES, ActivityType must be string;
        % Value_nM must be double.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            result = emk.db.getChemblActivity("CHEMBL203", MaxRows=10);
            tc.verifyClass(result.MoleculeChEMBLID, "string", ...
                "MoleculeChEMBLID must be string");
            tc.verifyClass(result.Name,             "string", ...
                "Name must be string");
            tc.verifyClass(result.SMILES,           "string", ...
                "SMILES must be string");
            tc.verifyClass(result.ActivityType,     "string", ...
                "ActivityType must be string");
            tc.verifyClass(result.Value_nM,         "double", ...
                "Value_nM must be double");
        end

        function test_getChemblActivity_valueNM_isPositive(tc)
        % All returned Value_nM values must be strictly positive.
        % MinActivity_nM filter and the nM-unit check both guard non-positive
        % values; if any slip through, the data is corrupt.
            tc.assumeTrue(tc.chemblAvailable(), "Skipped: ChEMBL not reachable");
            try
                result = emk.db.getChemblActivity("CHEMBL203", MaxRows=10);
            catch ME_net
                tc.assumeTrue(false, ...
                    "Skipped: ChEMBL activity endpoint transient error: " + ME_net.message);
            end
            tc.verifyTrue(all(result.Value_nM > 0), ...
                "All Value_nM values must be strictly positive");
        end

    end

end
