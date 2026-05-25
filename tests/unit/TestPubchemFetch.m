classdef TestPubchemFetch < matlab.unittest.TestCase
% TestPubchemFetch  Unit tests for emk.db.pubchemFetch.
%
% Coverage strategy: input-validation tests fire without network or Python;
% integration tests (requiring pubchempy + network) are guarded.
%
% Input validation (no Python, no network required):
%   TC1 : Numeric non-scalar array => invalidInput
%   TC2 : Cell array => invalidInput
%   TC3 : Logical input => invalidInput
%   TC4 : Empty string => invalidInput
%   TC5 : Whitespace-only string => invalidInput
%   TC6 : char input accepted (coerced to string)
%   TC7 : Numeric scalar CID accepted (coerced to string)
%   TC8 : Invalid NameSpace => invalidNamespace
%   TC9 : NameSpace case-insensitive ("NAME" treated as "name")
%   TC10: Valid NameSpace values accepted (no error thrown)
%
% PubChemPy availability guard:
%   TC11: libraryNotFound error when pubchempy not installed
%         (runs only in environments where pubchempy is absent)
%
% Integration tests (require pubchempy + network):
%   TC12: Ethanol by name returns struct
%   TC13: Struct has expected field names
%   TC14: CID field equals 702 for ethanol
%   TC14b: CID is double scalar (verifyClass + verifySize)
%   TC15: MolecularWeight is double near 46.07 (AbsTol=0.5)
%   TC16: IsomericSMILES is non-empty string
%   TC17: InChIKey matches LFQSCWFLJHTTHZ-UHFFFAOYSA-N
%   TC18: Synonyms is non-empty string array
%   TC18b: Synonyms is a row vector (1 x N)
%   TC19: SMILES namespace returns same CID as name namespace
%   TC20: CID namespace (numeric) returns correct compound
%   TC21: notFound error for nonsense identifier
%   TC22: MaxSynonyms option limits Synonyms array length
%   TC23: XLogP is a double scalar for ethanol
%   TC24: TPSA is a double scalar for ethanol (reference ~20.23 Ang^2)
%   TC25: HBondDonors = 1 for ethanol
%   TC25b: HBondAcceptors = 1 for ethanol (oxygen as acceptor)
%   TC25c: RotatableBonds = 0 for ethanol
%   TC25d: Charge = 0 for ethanol (neutral molecule)
%   TC25e: Complexity is a finite double scalar for ethanol
%   TC25f: MolecularFormula = "C2H6O" for ethanol
%   TC25g: IUPACName is non-empty string for ethanol
%   TC26: HeavyAtomCount = 3 for ethanol (C2 + O = 3 heavy atoms)
%   TC27: InChI starts with "InChI=" for ethanol
%   TC28: Cross-validation -- pubchemFetch CID matches searchPubchem CID
%
% Run with:
%   addpath(genpath("src")); addpath(genpath("tests"));
%   results = run(TestPubchemFetch);

    properties (TestParameter)
        validNamespace = {"name", "smiles", "cid", "inchi", "inchikey", "formula"}
    end

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    % Input validation tests (no Python / network required)
    % ======================================================================
    methods (Test)

        function test_nonScalarNumeric_throwsInvalidInput(tc)
        % Non-scalar numeric array must be rejected before any Python call.
            tc.verifyError(@() emk.db.pubchemFetch([702, 241]), ...
                "emk:db:pubchemFetch:invalidInput", ...
                "Non-scalar numeric must throw invalidInput");
        end

        function test_cellInput_throwsInvalidInput(tc)
            tc.verifyError(@() emk.db.pubchemFetch({"ethanol"}), ...
                "emk:db:pubchemFetch:invalidInput", ...
                "Cell input must throw invalidInput");
        end

        function test_logicalInput_throwsInvalidInput(tc)
            tc.verifyError(@() emk.db.pubchemFetch(true), ...
                "emk:db:pubchemFetch:invalidInput", ...
                "Logical scalar must throw invalidInput");
        end

        function test_emptyString_throwsInvalidInput(tc)
            tc.verifyError(@() emk.db.pubchemFetch(""), ...
                "emk:db:pubchemFetch:invalidInput", ...
                "Empty string must throw invalidInput");
        end

        function test_whitespaceOnly_throwsInvalidInput(tc)
            tc.verifyError(@() emk.db.pubchemFetch("   "), ...
                "emk:db:pubchemFetch:invalidInput", ...
                "Whitespace-only string must throw invalidInput");
        end

        function test_charInput_acceptedBeforePython(tc)
        % char input must be coerced to string; any subsequent error must NOT
        % be invalidInput (it may be libraryNotFound or networkError).
            err = tc.captureError(@() emk.db.pubchemFetch(char("ethanol")));
            if ~isempty(err)
                tc.verifyNotEqual(err.identifier, "emk:db:pubchemFetch:invalidInput", ...
                    "char input must not trigger invalidInput");
            end
        end

        function test_numericScalarCID_acceptedBeforePython(tc)
        % Numeric scalar CID must not trigger invalidInput.
            err = tc.captureError(@() emk.db.pubchemFetch(702, NameSpace="cid"));
            if ~isempty(err)
                tc.verifyNotEqual(err.identifier, "emk:db:pubchemFetch:invalidInput", ...
                    "Numeric scalar CID must not trigger invalidInput");
            end
        end

        function test_invalidNamespace_throwsInvalidNamespace(tc)
            tc.verifyError( ...
                @() emk.db.pubchemFetch("ethanol", NameSpace="badns"), ...
                "emk:db:pubchemFetch:invalidNamespace", ...
                "Unknown NameSpace must throw invalidNamespace");
        end

        function test_namespaceUppercase_acceptedBeforePython(tc)
        % NameSpace is lowercased internally; "NAME" must not trigger invalidNamespace.
            err = tc.captureError(@() emk.db.pubchemFetch("ethanol", NameSpace="NAME"));
            if ~isempty(err)
                tc.verifyNotEqual(err.identifier, "emk:db:pubchemFetch:invalidNamespace", ...
                    "Uppercase NameSpace must not trigger invalidNamespace");
            end
        end

        function test_validNamespaces_noInvalidNamespaceError(tc, validNamespace)
        % All documented NameSpace values must pass validation.
            err = tc.captureError(@() emk.db.pubchemFetch("x", ...
                NameSpace=string(validNamespace)));
            if ~isempty(err)
                tc.verifyNotEqual(err.identifier, "emk:db:pubchemFetch:invalidNamespace", ...
                    sprintf("Valid NameSpace must not throw invalidNamespace for: %s", validNamespace));
            end
        end

    end   % Input validation

    % ======================================================================
    % TC11: libraryNotFound (runs only when pubchempy is absent)
    % ======================================================================
    methods (Test)

        function test_libraryNotFound_throwsLibraryNotFound(tc)
        % TC11: When pubchempy is NOT installed, pubchemFetch must throw
        % libraryNotFound before any network request is made.
        % This test runs only in environments that lack pubchempy.
            tc.assumeFalse(tc.pubchemPyAvailable(), ...
                "Skipped: pubchempy IS installed (TC11 only runs when absent)");
            tc.verifyError(@() emk.db.pubchemFetch("ethanol"), ...
                "emk:db:pubchemFetch:libraryNotFound", ...
                "Absent pubchempy must throw libraryNotFound");
        end

    end   % TC11

    % ======================================================================
    % Integration tests (require pubchempy + network)
    % ======================================================================
    methods (Test)

        function test_ethanolByName_returnsStruct(tc)
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s, "struct", "Result must be a struct");
        end

        function test_ethanolByName_hasExpectedFields(tc)
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            expectedFields = ["CID","IUPACName","MolecularFormula","MolecularWeight", ...
                              "IsomericSMILES","InChI","InChIKey","XLogP","TPSA", ...
                              "HBondDonors","HBondAcceptors","RotatableBonds", ...
                              "HeavyAtomCount","Charge","Complexity","Synonyms"];
            for k = 1:numel(expectedFields)
                tc.verifyTrue(isfield(s, expectedFields(k)), ...
                    "Missing field: " + expectedFields(k));
            end
        end

        function test_ethanolCID_is702(tc)
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyEqual(s.CID, 702, "Ethanol CID must be 702");
        end

        function test_ethanolMW_nearRef(tc)
        % MolecularWeight must be near 46.07 g/mol (AbsTol=0.5 for rounding).
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s.MolecularWeight, "double", ...
                "MolecularWeight must be double");
            tc.verifyEqual(s.MolecularWeight, 46.07, "AbsTol", 0.5, ...
                "Ethanol MW must be near 46.07 g/mol");
        end

        function test_ethanolSMILES_nonEmpty(tc)
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s.IsomericSMILES, "string", ...
                "IsomericSMILES must be string");
            tc.verifyGreaterThan(strlength(s.IsomericSMILES), 0, ...
                "IsomericSMILES must not be empty");
        end

        function test_ethanolInChIKey_matches(tc)
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyEqual(s.InChIKey, "LFQSCWFLJHTTHZ-UHFFFAOYSA-N", ...
                "Ethanol InChIKey must match reference");
        end

        function test_ethanolSynonyms_nonEmpty(tc)
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s.Synonyms, "string", "Synonyms must be string array");
            tc.verifyGreaterThan(numel(s.Synonyms), 0, ...
                "Synonyms must contain at least one entry");
        end

        function test_smileNamespace_sameCID(tc)
        % SMILES search for ethanol must return CID 702.
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("CCO", NameSpace="smiles");
            tc.verifyEqual(s.CID, 702, ...
                "SMILES fetch must return same CID as name fetch for ethanol");
        end

        function test_cidNamespace_numericInput(tc)
        % Numeric CID 702 must return ethanol.
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch(702, NameSpace="cid");
            tc.verifyEqual(s.CID, 702, ...
                "CID namespace with numeric input must return CID 702");
        end

        function test_notFound_throwsNotFound(tc)
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            tc.verifyError( ...
                @() emk.db.pubchemFetch("xyzzy_no_such_compound_emk_test"), ...
                "emk:db:pubchemFetch:notFound", ...
                "Nonsense identifier must throw notFound");
        end

        function test_maxSynonyms_limitsLength(tc)
        % MaxSynonyms=3 must not return more than 3 synonyms.
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("aspirin", MaxSynonyms=3);
            tc.verifyLessThanOrEqual(numel(s.Synonyms), 3, ...
                "MaxSynonyms=3 must limit Synonyms array to <= 3 entries");
        end

        function test_ethanolXLogP_isDoubleScalar(tc)
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s.XLogP, "double", "XLogP must be double");
            tc.verifySize(s.XLogP, [1 1], "XLogP must be scalar");
        end

        function test_ethanolTPSA_isDoubleScalar(tc)
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s.TPSA, "double", "TPSA must be double");
            tc.verifyEqual(s.TPSA, 20.23, "AbsTol", 0.5, ...
                "Ethanol TPSA reference ~20.23 Ang^2");
        end

        function test_ethanolHBondDonors_equals1(tc)
        % Ethanol has one hydroxyl group => 1 H-bond donor.
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyEqual(s.HBondDonors, 1, ...
                "Ethanol must have exactly 1 H-bond donor");
        end

        function test_ethanolHeavyAtomCount_equals3(tc)
        % Ethanol C2H5OH has 3 heavy atoms (2C + 1O).
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyEqual(s.HeavyAtomCount, 3, ...
                "Ethanol must have 3 heavy atoms");
        end

        function test_ethanolInChI_startsWithPrefix(tc)
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifySubstring(s.InChI, "InChI=", ...
                "InChI must start with 'InChI=' prefix");
        end

        % --- TC14b: CID type + size ---
        function test_ethanolCID_isDoubleScalar(tc)
        % TC14b: CID must be a double scalar, not int or array.
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s.CID, "double", "CID must be double");
            tc.verifySize(s.CID, [1 1], "CID must be a scalar");
        end

        % --- TC18b: Synonyms row vector ---
        function test_ethanolSynonyms_isRowVector(tc)
        % TC18b: Synonyms must be a row vector (1 x N), not a column.
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyEqual(size(s.Synonyms, 1), 1, ...
                "Synonyms must be a row vector (first dim == 1)");
        end

        % --- TC25b: HBondAcceptors ---
        function test_ethanolHBondAcceptors_equals1(tc)
        % TC25b: Ethanol oxygen is an H-bond acceptor => HBondAcceptors = 1.
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s.HBondAcceptors, "double", "HBondAcceptors must be double");
            tc.verifyEqual(s.HBondAcceptors, 1, ...
                "Ethanol must have exactly 1 H-bond acceptor");
        end

        % --- TC25c: RotatableBonds ---
        function test_ethanolRotatableBonds_isZero(tc)
        % TC25c: Ethanol has no rotatable bonds (PubChem CID 702).
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s.RotatableBonds, "double", "RotatableBonds must be double");
            tc.verifyEqual(s.RotatableBonds, 0, ...
                "Ethanol must have 0 rotatable bonds");
        end

        % --- TC25d: Charge ---
        function test_ethanolCharge_isZero(tc)
        % TC25d: Ethanol is neutral; Charge must be 0 (double scalar).
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s.Charge, "double", "Charge must be double");
            tc.verifySize(s.Charge, [1 1], "Charge must be scalar");
            tc.verifyEqual(s.Charge, 0, "Ethanol Charge must be 0 (neutral)");
        end

        % --- TC25e: Complexity ---
        function test_ethanolComplexity_isFiniteDoubleScalar(tc)
        % TC25e: Complexity (Bertz CT) for ethanol is a finite double
        % scalar (~10.4). NaN is allowed only for compounds with no value.
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s.Complexity, "double", "Complexity must be double");
            tc.verifySize(s.Complexity, [1 1], "Complexity must be scalar");
            tc.verifyFalse(isinf(s.Complexity), "Complexity must not be Inf");
            tc.verifyGreaterThan(s.Complexity, 0, ...
                "Ethanol Complexity must be positive (PubChem: ~10.4)");
        end

        % --- TC25f: MolecularFormula ---
        function test_ethanolMolecularFormula_isC2H6O(tc)
        % TC25f: Ethanol molecular formula must be "C2H6O".
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s.MolecularFormula, "string", ...
                "MolecularFormula must be string");
            tc.verifyEqual(s.MolecularFormula, "C2H6O", ...
                "Ethanol MolecularFormula must be 'C2H6O'");
        end

        % --- TC25g: IUPACName ---
        function test_ethanolIUPACName_nonEmpty(tc)
        % TC25g: IUPACName for ethanol must be a non-empty string.
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            s = emk.db.pubchemFetch("ethanol");
            tc.verifyClass(s.IUPACName, "string", "IUPACName must be string");
            tc.verifyGreaterThan(strlength(s.IUPACName), 0, ...
                "IUPACName must not be empty for ethanol");
        end

        % --- TC28: Cross-validation with searchPubchem ---
        function test_crossValidate_withSearchPubchem_CID(tc)
        % TC28: CID returned by pubchemFetch (PubChemPy) must match CID
        % returned by searchPubchem (PUG REST / webread) for ethanol.
        % This cross-validates the two independent PubChem access paths.
            tc.assumeTrue(tc.pubchemPyAvailable(), "Skipped: pubchempy not installed");
            tc.assumeTrue(tc.networkAvailable(), "Skipped: network unavailable");

            sFetch  = emk.db.pubchemFetch("ethanol");
            tSearch = emk.db.searchPubchem("ethanol");
            tc.verifyEqual(sFetch.CID, tSearch.CID(1), ...
                "pubchemFetch CID must match searchPubchem CID for ethanol");
        end

    end   % Integration tests

    % ======================================================================
    % Helper methods
    % ======================================================================
    methods (Access = private)

        function tf = pubchemPyAvailable(~)
        % Return true when pubchempy is importable in the active Python env.
            try
                py.importlib.import_module("pubchempy");
                tf = true;
            catch
                tf = false;
            end
        end

        function tf = networkAvailable(~)
        % Probe PubChem with a short-timeout webread to check connectivity.
            try
                opts = weboptions("Timeout", 5);
                webread("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/702/JSON", opts);
                tf = true;
            catch
                tf = false;
            end
        end

        function ME = captureError(~, f)
        % Run f() and capture any MException; return [] if no error thrown.
            ME = [];
            try
                f();
            catch caught
                ME = caught;
            end
        end

    end
end
