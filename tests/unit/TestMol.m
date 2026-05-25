classdef TestMol < matlab.unittest.TestCase
% TestMol  Unit tests for src/+emk/+mol/fromSmiles.m, toSmiles.m, isValid.m
%
% Run with:
%   addpath(genpath("src"));
%   results = run(TestMol);
%
% Coverage (fromSmiles):
%   TC1: Non-string/char input => invalidInput error (no RDKit required)
%   TC2a: Empty / whitespace-only SMILES => invalidSmiles (no RDKit required;
%         explicit guard in fromSmiles fires before RDKit call)
%   TC2b: Invalid SMILES rejected by RDKit => invalidSmiles (requires RDKit)
%   TC3: Valid SMILES => correct Mol, atom count, bond count (requires RDKit)
%   TC4: char input accepted (requires RDKit)
%   TC5: Complex molecule (aspirin) parsed correctly -- atom count + bond count
%        (requires RDKit)
%
% Coverage (toSmiles):
%   TC6: Non-Mol input => invalidInput error (no RDKit required)
%         + error message contains offending class name
%         + char and logical inputs rejected
%   TC7: Valid Mol => canonical SMILES string + scalar type (requires RDKit)
%         + round-trip preserves atom count AND bond count
%         + complex molecule (aspirin) round-trip -- atom count + bond count
%   TC8: Non-canonical input SMILES is normalised (requires RDKit)
%         + canonical SMILES is idempotent (fixed-point property)
%
% Coverage (isValid):
%   TC9:  Non-string/char inputs return false, no throw (no RDKit required)
%         + non-scalar string array, empty matrix also return false
%         + return type is logical scalar
%   TC9b: Empty / whitespace-only SMILES return false (no RDKit required)
%         + empty char '' return false
%   TC10: Invalid SMILES returns false (requires RDKit)
%         + invalid char literal also returns false
%   TC11: Valid SMILES returns true; char literal accepted (requires RDKit)
%         + return type is logical scalar for valid input
%
% Tests requiring RDKit use assumeTrue(tc.rdkitAvailable()) to skip
% gracefully when Python/RDKit is not configured in the current session.
%
% ======================================================================
% Design note on py.rdkit.* calls:
%   Atom count tests use emk.descriptor.calculate(mol,"HeavyAtomCount")
%   rather than mol.GetNumAtoms() (py.rdkit.* direct call).  Bond count
%   tests call mol.GetNumBonds() directly because no emk.* wrapper
%   exists for bond count; these are labelled as cross-validation tests
%   (principle 3) and verify the Mol object's connectivity beyond what
%   atom count alone can confirm.
%   Return-type checks use startsWith/contains instead of verifyClass
%   with a hardcoded "py.rdkit.Chem.rdchem.Mol" string to stay stable
%   across RDKit versions (mirrors TestFingerprint.m design).

    methods (TestMethodSetup)
        function setupPath(tc) %#ok<MANU>
            addpath(genpath("src"));
        end
    end

    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC1: Input validation -- no RDKit required
        % ------------------------------------------------------------------

        function test_fromSmiles_numericInput_throwsInvalidInput(tc)
        % Numeric input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.mol.fromSmiles(42), ...
                "emk:mol:fromSmiles:invalidInput", ...
                "Numeric input must throw invalidInput");
        end

        function test_fromSmiles_numericInput_errorMessage_containsClass(tc)
        % Error message must contain the actual class name for debugging.
            ME = tc.captureError(@() emk.mol.fromSmiles(42));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain the offending class name");
        end

        function test_fromSmiles_cellInput_throwsInvalidInput(tc)
        % Cell array input must throw invalidInput.
            tc.verifyError(@() emk.mol.fromSmiles({"CCO"}), ...
                "emk:mol:fromSmiles:invalidInput", ...
                "Cell array input must throw invalidInput");
        end

        function test_fromSmiles_stringArrayInput_throwsInvalidInput(tc)
        % Non-scalar string array input must throw invalidInput.
            tc.verifyError(@() emk.mol.fromSmiles(["CCO", "c1ccccc1"]), ...
                "emk:mol:fromSmiles:invalidInput", ...
                "Non-scalar string array must throw invalidInput");
        end

        function test_fromSmiles_logicalInput_throwsInvalidInput(tc)
        % logical input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.mol.fromSmiles(true), ...
                "emk:mol:fromSmiles:invalidInput", ...
                "logical input must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC2a: Empty / whitespace-only SMILES -- no RDKit required
        %       (explicit guard in fromSmiles fires before RDKit call)
        % ------------------------------------------------------------------

        function test_fromSmiles_emptyString_throwsInvalidSmiles(tc)
        % Empty SMILES string must throw invalidSmiles.
        % No RDKit required: the explicit empty-string guard in fromSmiles
        % fires before MolFromSmiles is called, so this test is independent
        % of the installed RDKit version.
            tc.verifyError(@() emk.mol.fromSmiles(""), ...
                "emk:mol:fromSmiles:invalidSmiles", ...
                "Empty SMILES must throw invalidSmiles");
        end

        function test_fromSmiles_whitespaceString_throwsInvalidSmiles(tc)
        % Whitespace-only string must throw invalidSmiles (strtrim collapses to
        % empty before the length check).
        % No RDKit required: same pre-RDKit guard as empty string.
            tc.verifyError(@() emk.mol.fromSmiles("   "), ...
                "emk:mol:fromSmiles:invalidSmiles", ...
                "Whitespace-only SMILES must throw invalidSmiles");
        end

        % ------------------------------------------------------------------
        % TC2b: Invalid SMILES rejected by RDKit (requires RDKit)
        % ------------------------------------------------------------------

        function test_fromSmiles_invalidSmiles_throwsInvalidSmiles(tc)
        % An unparseable SMILES must throw invalidSmiles.
            tc.assumeTrue(tc.rdkitAvailable(), ...
                "Skipped: RDKit not available");
            tc.verifyError(@() emk.mol.fromSmiles("invalid_smiles_xyz"), ...
                "emk:mol:fromSmiles:invalidSmiles", ...
                "Invalid SMILES must throw invalidSmiles");
        end

        function test_fromSmiles_invalidSmiles_errorMessage_containsSmiles(tc)
        % Error message must contain the offending SMILES for debugging.
            tc.assumeTrue(tc.rdkitAvailable(), ...
                "Skipped: RDKit not available");
            badSmiles = "invalid_smiles_xyz";
            ME = tc.captureError(@() emk.mol.fromSmiles(badSmiles));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, char(badSmiles), ...
                "Error message must contain the offending SMILES string");
        end

        % ------------------------------------------------------------------
        % TC3: Valid SMILES => Mol object returned (requires RDKit)
        % ------------------------------------------------------------------

        function test_fromSmiles_ethanol_returnsMolObject(tc)
        % Valid SMILES (ethanol "CCO") must return a Python Mol object.
        % Uses startsWith/contains to remain stable across RDKit versions.
            tc.assumeTrue(tc.rdkitAvailable(), ...
                "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tc.verifyTrue(startsWith(class(mol), "py."), ...
                "fromSmiles must return a Python object (class starts with 'py.')");
            tc.verifyTrue(contains(class(mol), "Mol"), ...
                sprintf("fromSmiles must return a Mol type, got: %s", class(mol)));
        end

        function test_fromSmiles_benzene_returnsMolObject(tc)
        % Valid SMILES (benzene) must return a Python Mol object.
        % Mirrors ethanol type check; confirms the contract holds for
        % aromatic input (not only aliphatic chains).
            tc.assumeTrue(tc.rdkitAvailable(), ...
                "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("c1ccccc1");
            tc.verifyTrue(startsWith(class(mol), "py."), ...
                "fromSmiles must return a Python object for benzene");
            tc.verifyTrue(contains(class(mol), "Mol"), ...
                sprintf("fromSmiles must return a Mol type for benzene, got: %s", class(mol)));
        end

        function test_fromSmiles_ethanol_heavyAtomCount(tc)
        % Ethanol (CCO) has 3 heavy atoms (2C + 1O).
        % HeavyAtomCount via emk.descriptor.calculate avoids mol.GetNumAtoms().
        % verifyClass + verifySize confirm type contract (principle 2).
        % Reference: PubChem CID 702 (Ethanol).
            tc.assumeTrue(tc.rdkitAvailable(), ...
                "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "HeavyAtomCount");
            tc.verifyClass(desc.HeavyAtomCount, "double", ...
                "HeavyAtomCount must be a double scalar");
            tc.verifySize(desc.HeavyAtomCount, [1 1], ...
                "HeavyAtomCount must be a scalar [1 1]");
            tc.verifyEqual(desc.HeavyAtomCount, 3, ...
                "Ethanol (CCO) must have 3 heavy atoms");
        end

        function test_fromSmiles_benzene_heavyAtomCount(tc)
        % Benzene (c1ccccc1) has 6 heavy atoms (6C).
        % Uses emk.descriptor.calculate to avoid mol.GetNumAtoms() direct call.
        % Reference: PubChem CID 241 (Benzene).
            tc.assumeTrue(tc.rdkitAvailable(), ...
                "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccccc1");
            desc = emk.descriptor.calculate(mol, "HeavyAtomCount");
            tc.verifyClass(desc.HeavyAtomCount, "double", ...
                "HeavyAtomCount must be a double scalar");
            tc.verifySize(desc.HeavyAtomCount, [1 1], ...
                "HeavyAtomCount must be a scalar [1 1]");
            tc.verifyEqual(desc.HeavyAtomCount, 6, ...
                "Benzene (c1ccccc1) must have 6 heavy atoms");
        end

        function test_fromSmiles_ethanol_bondCount(tc)
        % Ethanol (CCO) has 2 bonds: C-C and C-O.
        % Cross-validation: mol.GetNumBonds() is called directly because
        % no emk.* wrapper exists for bond count.  Confirms connectivity
        % beyond atom count (principle 3).  Reference: PubChem CID 702.
            tc.assumeTrue(tc.rdkitAvailable(), ...
                "Skipped: RDKit not available");
            mol      = emk.mol.fromSmiles("CCO");
            numBonds = double(mol.GetNumBonds());   % cross-validation: no emk wrapper
            tc.verifyEqual(numBonds, 2, ...
                "Ethanol (CCO) must have 2 bonds (C-C and C-O)");
        end

        % ------------------------------------------------------------------
        % TC4: char input accepted (requires RDKit)
        % ------------------------------------------------------------------

        function test_fromSmiles_charInput_accepted(tc)
        % char literal input must be accepted (ischar || isStringScalar guard).
        % Return type must be a Python Mol object (same as string input).
            tc.assumeTrue(tc.rdkitAvailable(), ...
                "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles('CCO');   % char literal
            tc.verifyTrue(startsWith(class(mol), "py."), ...
                "char input must return a Python object");
            tc.verifyTrue(contains(class(mol), "Mol"), ...
                sprintf("char input must return a Mol type, got: %s", class(mol)));
        end

        % ------------------------------------------------------------------
        % TC5: Complex molecule -- aspirin (requires RDKit)
        % Verifies that aromatic rings and multiple functional groups are
        % handled correctly (not just simple aliphatic chains).
        % ------------------------------------------------------------------

        function test_fromSmiles_aspirin_returnsMolObject(tc)
        % Aspirin SMILES must parse successfully and return a Python Mol object.
        % SMILES: CC(=O)Oc1ccccc1C(=O)O  (PubChem CID 2244)
            tc.assumeTrue(tc.rdkitAvailable(), ...
                "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            tc.verifyTrue(startsWith(class(mol), "py."), ...
                "Aspirin SMILES must return a Python object");
            tc.verifyTrue(contains(class(mol), "Mol"), ...
                sprintf("Aspirin SMILES must return a Mol type, got: %s", class(mol)));
        end

        function test_fromSmiles_aspirin_heavyAtomCount(tc)
        % Aspirin (C9H8O4) has 13 heavy atoms: 9 C + 4 O.
        % Uses emk.descriptor.calculate to avoid mol.GetNumAtoms() direct call.
        % Confirms aromatic ring + multiple functional groups are parsed correctly.
        % Reference: PubChem CID 2244.
            tc.assumeTrue(tc.rdkitAvailable(), ...
                "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            desc = emk.descriptor.calculate(mol, "HeavyAtomCount");
            tc.verifyClass(desc.HeavyAtomCount, "double", ...
                "HeavyAtomCount must be a double scalar");
            tc.verifySize(desc.HeavyAtomCount, [1 1], ...
                "HeavyAtomCount must be a scalar [1 1]");
            tc.verifyEqual(desc.HeavyAtomCount, 13, ...
                "Aspirin (C9H8O4) must have 13 heavy atoms");
        end

        function test_fromSmiles_aspirin_bondCount(tc)
        % Aspirin (C9H8O4) has 13 bonds (ring closure: 6-membered ring
        % with 2 substituents gives 13 bonds for 13 heavy atoms).
        % Cross-validation: mol.GetNumBonds() is called directly because
        % no emk.* wrapper exists for bond count (principle 3).
        % Reference: PubChem CID 2244 (13 atoms, 13 bonds).
            tc.assumeTrue(tc.rdkitAvailable(), ...
                "Skipped: RDKit not available");
            mol      = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            numBonds = double(mol.GetNumBonds());   % cross-validation: no emk wrapper
            tc.verifyEqual(numBonds, 13, ...
                "Aspirin (C9H8O4) must have 13 bonds");
        end

        % ------------------------------------------------------------------
        % TC6: toSmiles -- Input validation (no RDKit required)
        % ------------------------------------------------------------------

        function test_toSmiles_numericInput_throwsInvalidInput(tc)
        % Numeric input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.mol.toSmiles(42), ...
                "emk:mol:toSmiles:invalidInput", ...
                "Numeric input must throw invalidInput");
        end

        function test_toSmiles_stringInput_throwsInvalidInput(tc)
        % Passing a raw SMILES string (not a Mol object) must throw invalidInput.
            tc.verifyError(@() emk.mol.toSmiles("CCO"), ...
                "emk:mol:toSmiles:invalidInput", ...
                "String input must throw invalidInput");
        end

        function test_toSmiles_emptyInput_throwsInvalidInput(tc)
        % Empty matrix input must throw invalidInput.
            tc.verifyError(@() emk.mol.toSmiles([]), ...
                "emk:mol:toSmiles:invalidInput", ...
                "Empty input must throw invalidInput");
        end

        function test_toSmiles_numericInput_errorMessage_containsClass(tc)
        % Error message must contain the actual class name for debugging.
        % Mirrors the analogous fromSmiles test (TC1 error message contract).
            ME = tc.captureError(@() emk.mol.toSmiles(42));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain the offending class name");
        end

        function test_toSmiles_charInput_throwsInvalidInput(tc)
        % char literal (not a Mol object) must throw invalidInput.
        % Ensures the isa() guard rejects all non-Mol types including char.
            tc.verifyError(@() emk.mol.toSmiles('CCO'), ...
                "emk:mol:toSmiles:invalidInput", ...
                "char input must throw invalidInput");
        end

        function test_toSmiles_logicalInput_throwsInvalidInput(tc)
        % logical input must throw invalidInput.
            tc.verifyError(@() emk.mol.toSmiles(true), ...
                "emk:mol:toSmiles:invalidInput", ...
                "logical input must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC7: toSmiles -- Valid Mol returns canonical SMILES (requires RDKit)
        % ------------------------------------------------------------------

        function test_toSmiles_ethanol_returnsString(tc)
        % toSmiles must return a MATLAB string scalar for ethanol.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol    = emk.mol.fromSmiles("CCO");
            smiles = emk.mol.toSmiles(mol);
            tc.verifyClass(smiles, "string", ...
                "toSmiles must return a MATLAB string");
            tc.verifyTrue(isStringScalar(smiles), ...
                "toSmiles must return a scalar string");
        end

        function test_toSmiles_ethanol_notEmpty(tc)
        % toSmiles result for ethanol must be non-empty.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol    = emk.mol.fromSmiles("CCO");
            smiles = emk.mol.toSmiles(mol);
            tc.verifyGreaterThan(strlength(smiles), 0, ...
                "toSmiles must return a non-empty string");
        end

        function test_toSmiles_benzene_returnsString(tc)
        % toSmiles must return a MATLAB string scalar for benzene.
        % Mirrors ethanol type check to confirm the contract holds for
        % aromatic input (not only aliphatic).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol    = emk.mol.fromSmiles("c1ccccc1");
            smiles = emk.mol.toSmiles(mol);
            tc.verifyClass(smiles, "string", ...
                "toSmiles must return a MATLAB string for benzene");
            tc.verifyTrue(isStringScalar(smiles), ...
                "toSmiles must return a scalar string for benzene");
        end

        function test_toSmiles_benzene_notEmpty(tc)
        % toSmiles result for benzene must be non-empty.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol    = emk.mol.fromSmiles("c1ccccc1");
            smiles = emk.mol.toSmiles(mol);
            tc.verifyGreaterThan(strlength(smiles), 0, ...
                "toSmiles must return a non-empty string for benzene");
        end

        function test_toSmiles_roundtrip_ethanolAtomCount(tc)
        % Round-trip: fromSmiles -> toSmiles -> fromSmiles must preserve
        % atom count (3 heavy atoms for ethanol CCO).
        % HeavyAtomCount via emk.descriptor.calculate avoids mol.GetNumAtoms().
        % Reference: PubChem CID 702 (Ethanol).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1   = emk.mol.fromSmiles("CCO");
            smiles = emk.mol.toSmiles(mol1);
            mol2   = emk.mol.fromSmiles(smiles);
            desc   = emk.descriptor.calculate(mol2, "HeavyAtomCount");
            tc.verifyClass(desc.HeavyAtomCount, "double", ...
                "Round-trip HeavyAtomCount must be double");
            tc.verifySize(desc.HeavyAtomCount, [1 1], ...
                "Round-trip HeavyAtomCount must be scalar");
            tc.verifyEqual(desc.HeavyAtomCount, 3, ...
                "Round-trip atom count must be 3 for ethanol");
        end

        function test_toSmiles_roundtrip_ethanolBondCount(tc)
        % Round-trip: fromSmiles -> toSmiles -> fromSmiles must preserve
        % bond count (2 bonds: C-C and C-O).
        % Cross-validation: mol.GetNumBonds() is called directly because
        % no emk.* wrapper exists for bond count (principle 3).
        % Reference: PubChem CID 702 (Ethanol).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1     = emk.mol.fromSmiles("CCO");
            smiles   = emk.mol.toSmiles(mol1);
            mol2     = emk.mol.fromSmiles(smiles);
            numBonds = double(mol2.GetNumBonds());   % cross-validation: no emk wrapper
            tc.verifyEqual(numBonds, 2, ...
                "Round-trip bond count must be 2 for ethanol (C-C and C-O)");
        end

        function test_toSmiles_aspirin_roundtrip_atomCount(tc)
        % Aspirin round-trip: fromSmiles -> toSmiles -> fromSmiles must
        % preserve heavy atom count (13 atoms: 9 C + 4 O).
        % HeavyAtomCount via emk.descriptor.calculate avoids mol.GetNumAtoms().
        % Verifies toSmiles correctly serialises aromatic + functional groups.
        % Reference: PubChem CID 2244.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            smiles = emk.mol.toSmiles(mol1);
            mol2   = emk.mol.fromSmiles(smiles);
            desc   = emk.descriptor.calculate(mol2, "HeavyAtomCount");
            tc.verifyClass(desc.HeavyAtomCount, "double", ...
                "Round-trip HeavyAtomCount must be double");
            tc.verifySize(desc.HeavyAtomCount, [1 1], ...
                "Round-trip HeavyAtomCount must be scalar");
            tc.verifyEqual(desc.HeavyAtomCount, 13, ...
                "Aspirin round-trip atom count must be 13");
        end

        function test_toSmiles_aspirin_roundtrip_bondCount(tc)
        % Aspirin round-trip: fromSmiles -> toSmiles -> fromSmiles must
        % preserve bond count (13 bonds for 13 heavy atoms).
        % Cross-validation: mol.GetNumBonds() is called directly because
        % no emk.* wrapper exists.  Confirms ring closure and ester/carboxyl
        % connectivity survive canonical serialisation (principle 3).
        % Reference: PubChem CID 2244 (13 atoms, 13 bonds).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1     = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            smiles   = emk.mol.toSmiles(mol1);
            mol2     = emk.mol.fromSmiles(smiles);
            numBonds = double(mol2.GetNumBonds());   % cross-validation: no emk wrapper
            tc.verifyEqual(numBonds, 13, ...
                "Aspirin round-trip bond count must be 13");
        end

        % ------------------------------------------------------------------
        % TC8: toSmiles -- Non-canonical input is normalised (requires RDKit)
        %
        % "OCC" is a valid but non-canonical SMILES for ethanol.
        % RDKit's canonical algorithm assigns unique atom indices; the
        % canonical output for all representations of ethanol is "CCO".
        % Reference: Weininger 1989 (SMILES notation), RDKit canonicalisation.
        % ------------------------------------------------------------------

        function test_toSmiles_nonCanonical_isNormalised(tc)
        % Non-canonical SMILES "OCC" must round-trip to the same canonical
        % form as "CCO" (both represent ethanol).
        % Validates that toSmiles always returns the canonical representation.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1 = emk.mol.fromSmiles("CCO");
            mol2 = emk.mol.fromSmiles("OCC");
            tc.verifyEqual(emk.mol.toSmiles(mol1), emk.mol.toSmiles(mol2), ...
                "Both ethanol representations must yield the same canonical SMILES");
        end

        function test_toSmiles_canonical_isIdempotent(tc)
        % Canonical SMILES must be a fixed point (idempotent):
        %   toSmiles(fromSmiles(s)) == toSmiles(fromSmiles(toSmiles(fromSmiles(s))))
        % This is a fundamental mathematical property of canonical form.
        % If canonical SMILES were not idempotent, it would not be a true
        % canonical form.  Starting from non-canonical "OCC" exercises the
        % full normalisation path before the idempotency check.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol1    = emk.mol.fromSmiles("OCC");
            smiles1 = emk.mol.toSmiles(mol1);
            mol2    = emk.mol.fromSmiles(smiles1);
            smiles2 = emk.mol.toSmiles(mol2);
            tc.verifyEqual(smiles1, smiles2, ...
                "Canonical SMILES must be idempotent (fixed-point property)");
        end

        % ------------------------------------------------------------------
        % TC9: isValid -- Non-string inputs return false (no RDKit required)
        % ------------------------------------------------------------------

        function test_isValid_numericInput_returnsFalse(tc)
        % Numeric input must return false without throwing.
            tf = emk.mol.isValid(42);
            tc.verifyFalse(tf, "Numeric input must return false without throw");
        end

        function test_isValid_cellInput_returnsFalse(tc)
        % Cell array input must return false without throwing.
            tf = emk.mol.isValid({"CCO"});
            tc.verifyFalse(tf, "Cell array input must return false without throw");
        end

        function test_isValid_logicalInput_returnsFalse(tc)
        % Logical input must return false without throwing.
            tf = emk.mol.isValid(true);
            tc.verifyFalse(tf, "Logical input must return false without throw");
        end

        function test_isValid_stringArrayInput_returnsFalse(tc)
        % Non-scalar string array must return false without throwing.
        % isStringScalar(["CCO","XYZ"]) == false, so the type guard fires.
            tf = emk.mol.isValid(["CCO", "XYZ"]);
            tc.verifyFalse(tf, "Non-scalar string array must return false without throw");
        end

        function test_isValid_emptyMatrix_returnsFalse(tc)
        % Empty matrix [] must return false without throwing.
        % ischar([]) == false, isStringScalar([]) == false -> false.
            tf = emk.mol.isValid([]);
            tc.verifyFalse(tf, "Empty matrix must return false without throw");
        end

        function test_isValid_numericInput_isLogicalScalar(tc)
        % Return value for non-string input must be a logical scalar.
            tf = emk.mol.isValid(42);
            tc.verifyClass(tf, "logical", "isValid must return logical for non-string input");
            tc.verifyTrue(isscalar(tf), "isValid must return a scalar for non-string input");
        end

        % ------------------------------------------------------------------
        % TC9b: isValid -- Empty / whitespace (no RDKit required)
        % ------------------------------------------------------------------

        function test_isValid_emptyString_returnsFalse(tc)
        % Empty SMILES string must return false without throwing.
        % Mirrors the fromSmiles empty-string guard but returns false
        % instead of throwing (isValid never throws for parse failures).
            tf = emk.mol.isValid("");
            tc.verifyFalse(tf, "Empty string must return false");
        end

        function test_isValid_emptyChar_returnsFalse(tc)
        % Empty char literal '' must return false without throwing.
        % ischar('') == true, so it passes the type guard, then fromSmiles
        % fires the empty-string guard and throws invalidSmiles, which
        % isValid translates to false.
            tf = emk.mol.isValid('');
            tc.verifyFalse(tf, "Empty char must return false");
        end

        function test_isValid_whitespaceString_returnsFalse(tc)
        % Whitespace-only SMILES string must return false without throwing.
            tf = emk.mol.isValid("   ");
            tc.verifyFalse(tf, "Whitespace-only string must return false");
        end

        % ------------------------------------------------------------------
        % TC10: isValid -- Invalid SMILES returns false (requires RDKit)
        % ------------------------------------------------------------------

        function test_isValid_invalidSmiles_returnsFalse(tc)
        % An unparseable SMILES string must return false without throwing.
        % Validates that RDKit rejection is translated to false, not an error.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tf = emk.mol.isValid("invalid_smiles_xyz");
            tc.verifyFalse(tf, "Invalid SMILES must return false without throw");
        end

        function test_isValid_invalidCharInput_returnsFalse(tc)
        % char literal with an unparseable SMILES must return false.
        % Complements test_isValid_charInput_returnsTrue: char is accepted
        % as input type, but an invalid SMILES still returns false.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tf = emk.mol.isValid('invalid_smiles_xyz');
            tc.verifyFalse(tf, "char invalid SMILES must return false without throw");
        end

        % ------------------------------------------------------------------
        % TC11: isValid -- Valid SMILES returns true (requires RDKit)
        % ------------------------------------------------------------------

        function test_isValid_ethanol_returnsTrue(tc)
        % Valid SMILES for ethanol must return true.
        % Reference: PubChem CID 702 (Ethanol).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tf = emk.mol.isValid("CCO");
            tc.verifyTrue(tf, "Ethanol SMILES must return true");
        end

        function test_isValid_benzene_returnsTrue(tc)
        % Valid aromatic SMILES for benzene must return true.
        % Confirms that aromatic notation is accepted.
        % Reference: PubChem CID 241 (Benzene).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tf = emk.mol.isValid("c1ccccc1");
            tc.verifyTrue(tf, "Benzene SMILES must return true");
        end

        function test_isValid_aspirin_returnsTrue(tc)
        % Valid SMILES for aspirin (complex molecule) must return true.
        % SMILES: CC(=O)Oc1ccccc1C(=O)O  (PubChem CID 2244)
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tf = emk.mol.isValid("CC(=O)Oc1ccccc1C(=O)O");
            tc.verifyTrue(tf, "Aspirin SMILES must return true");
        end

        function test_isValid_charInput_returnsTrue(tc)
        % char literal input for a valid SMILES must return true.
        % Confirms that char (not only string) is accepted.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tf = emk.mol.isValid('CCO');   % char literal
            tc.verifyTrue(tf, "char input for valid SMILES must return true");
        end

        function test_isValid_result_isLogicalScalar(tc)
        % Return value for valid SMILES must be a logical scalar.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            tf = emk.mol.isValid("CCO");
            tc.verifyClass(tf, "logical", "isValid must return logical");
            tc.verifyTrue(isscalar(tf), "isValid must return a scalar");
        end

        % ------------------------------------------------------------------
        % TC12: hasSubstruct -- Input validation (no RDKit required)
        % ------------------------------------------------------------------

        function test_hasSubstruct_numericMol_throwsInvalidMol(tc)
        % Numeric mol input must throw invalidMol before any RDKit call.
            tc.verifyError(@() emk.mol.hasSubstruct(42, "c1ccccc1"), ...
                "emk:mol:hasSubstruct:invalidMol", ...
                "Numeric mol must throw invalidMol");
        end

        function test_hasSubstruct_stringMol_throwsInvalidMol(tc)
        % Raw SMILES string as mol must throw invalidMol (not a Mol object).
            tc.verifyError(@() emk.mol.hasSubstruct("CCO", "c1ccccc1"), ...
                "emk:mol:hasSubstruct:invalidMol", ...
                "String mol must throw invalidMol");
        end

        function test_hasSubstruct_emptyCellMol_throwsInvalidMol(tc)
        % Empty cell array as mol must throw invalidMol.
            tc.verifyError(@() emk.mol.hasSubstruct({}, "c1ccccc1"), ...
                "emk:mol:hasSubstruct:invalidMol", ...
                "Empty cell must throw invalidMol");
        end

        function test_hasSubstruct_cellWithNonMolElement_throwsInvalidMol(tc)
        % Cell array where an element is NOT a Mol must throw invalidMol.
        % cellfun in hasSubstruct validates every element; {42} contains a
        % double which fails isa(m, "py.rdkit.Chem.rdchem.Mol").
        % No RDKit required: mol validation fires before any RDKit call.
            tc.verifyError(@() emk.mol.hasSubstruct({42}, "c1ccccc1"), ...
                "emk:mol:hasSubstruct:invalidMol", ...
                "Cell with non-Mol element must throw invalidMol");
        end

        function test_hasSubstruct_invalidMol_errorContainsClass(tc)
        % Error message for invalid mol must contain the offending class name.
            ME = tc.captureError(@() emk.mol.hasSubstruct(42, "c1ccccc1"));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain the offending class name");
        end

        function test_hasSubstruct_numericQuery_throwsInvalidQuery(tc)
        % Numeric query must throw invalidQuery before any RDKit call.
        % mol validation fires first for invalid mol types, so we need a
        % valid Mol for this test; use a known-valid one if available, else skip.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tc.verifyError(@() emk.mol.hasSubstruct(mol, 42), ...
                "emk:mol:hasSubstruct:invalidQuery", ...
                "Numeric query must throw invalidQuery");
        end

        function test_hasSubstruct_emptyStringQuery_throwsInvalidQuery(tc)
        % Empty SMARTS string query must throw invalidQuery.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tc.verifyError(@() emk.mol.hasSubstruct(mol, ""), ...
                "emk:mol:hasSubstruct:invalidQuery", ...
                "Empty SMARTS string must throw invalidQuery");
        end

        function test_hasSubstruct_whitespaceQuery_throwsInvalidQuery(tc)
        % Whitespace-only SMARTS string must throw invalidQuery.
        % Mirrors fromSmiles whitespace guard: strtrim("   ") collapses to empty.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tc.verifyError(@() emk.mol.hasSubstruct(mol, "   "), ...
                "emk:mol:hasSubstruct:invalidQuery", ...
                "Whitespace-only SMARTS must throw invalidQuery");
        end

        function test_hasSubstruct_invalidSmarts_throwsInvalidQuery(tc)
        % An unparseable SMARTS pattern must throw invalidQuery (not rdkitError).
        % '[C' (unclosed bracket) causes RDKit MolFromSmarts to return None,
        % which hasSubstruct translates to invalidQuery via isa(x, "py.NoneType").
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tc.verifyError(@() emk.mol.hasSubstruct(mol, "[C"), ...
                "emk:mol:hasSubstruct:invalidQuery", ...
                "Unclosed-bracket SMARTS must throw invalidQuery");
        end

        function test_hasSubstruct_invalidQuery_errorContainsSmarts(tc)
        % Error message for an invalid SMARTS must contain the offending pattern.
        % Mirrors fromSmiles:invalidSmiles message contract.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol       = emk.mol.fromSmiles("CCO");
            badSmarts = "[C";
            ME = tc.captureError(@() emk.mol.hasSubstruct(mol, badSmarts));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, char(badSmarts), ...
                "Error message must contain the offending SMARTS string");
        end

        % ------------------------------------------------------------------
        % TC13: hasSubstruct -- Single mol, SMARTS query (requires RDKit)
        %
        % Principle 1: emk.* wrapper used throughout (no py.rdkit.* calls).
        % Principle 2: verifyClass + verifySize confirm type contract.
        % Principle 3: cross-validate true/false with chemically meaningful cases.
        % ------------------------------------------------------------------

        function test_hasSubstruct_benzeneRingInBenzylAlcohol_returnsTrue(tc)
        % Benzyl alcohol (c1ccccc1CO) contains a benzene ring (c1ccccc1).
        % SMARTS query: aromatic 6-membered ring.
        % Reference: PubChem CID 244 (Benzyl alcohol, C7H8O).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("c1ccccc1CO");
            tf  = emk.mol.hasSubstruct(mol, "c1ccccc1");
            tc.verifyClass(tf, "logical", ...
                "hasSubstruct must return logical");
            tc.verifySize(tf, [1 1], ...
                "hasSubstruct must return a scalar [1 1]");
            tc.verifyTrue(tf, ...
                "Benzyl alcohol must contain a benzene ring substructure");
        end

        function test_hasSubstruct_amineAbsentInBenzylAlcohol_returnsFalse(tc)
        % Benzyl alcohol (c1ccccc1CO) does NOT contain a primary amine [NH2].
        % Cross-validates against the true case above (same molecule, different query).
        % Reference: PubChem CID 244.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("c1ccccc1CO");
            tf  = emk.mol.hasSubstruct(mol, "[NH2]");
            tc.verifyClass(tf, "logical", ...
                "hasSubstruct must return logical");
            tc.verifySize(tf, [1 1], ...
                "hasSubstruct must return a scalar [1 1]");
            tc.verifyFalse(tf, ...
                "Benzyl alcohol must NOT contain a primary amine substructure");
        end

        function test_hasSubstruct_hydroxylInEthanol_returnsTrue(tc)
        % Ethanol (CCO) contains a hydroxyl group ([OH]).
        % verifyClass + verifySize confirm type contract (principle 2).
        % Reference: PubChem CID 702 (Ethanol).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tf  = emk.mol.hasSubstruct(mol, "[OH]");
            tc.verifyClass(tf, "logical", ...
                "hasSubstruct must return logical");
            tc.verifySize(tf, [1 1], ...
                "hasSubstruct must return a scalar [1 1]");
            tc.verifyTrue(tf, "Ethanol must contain a hydroxyl group [OH]");
        end

        function test_hasSubstruct_benzeneAbsentInEthanol_returnsFalse(tc)
        % Ethanol (CCO) does NOT contain a benzene ring (c1ccccc1).
        % verifyClass + verifySize confirm type contract (principle 2).
        % Reference: PubChem CID 702 (Ethanol).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tf  = emk.mol.hasSubstruct(mol, "c1ccccc1");
            tc.verifyClass(tf, "logical", ...
                "hasSubstruct must return logical");
            tc.verifySize(tf, [1 1], ...
                "hasSubstruct must return a scalar [1 1]");
            tc.verifyFalse(tf, "Ethanol must NOT contain a benzene ring");
        end

        function test_hasSubstruct_carbonylInAspirin_returnsTrue(tc)
        % Aspirin (CC(=O)Oc1ccccc1C(=O)O) contains a carbonyl group [CX3](=O).
        % verifyClass + verifySize confirm type contract (principle 2).
        % Reference: PubChem CID 2244 (Aspirin, C9H8O4).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            tf  = emk.mol.hasSubstruct(mol, "[CX3](=O)");
            tc.verifyClass(tf, "logical", ...
                "hasSubstruct must return logical");
            tc.verifySize(tf, [1 1], ...
                "hasSubstruct must return a scalar [1 1]");
            tc.verifyTrue(tf, "Aspirin must contain a carbonyl group [CX3](=O)");
        end

        function test_hasSubstruct_ketoneAbsentInCaffeine_returnsFalse(tc)
        % Caffeine (CN1C=NC2=C1C(=O)N(C(=O)N2C)C) has TWO carbonyl groups,
        % but BOTH are amide/imide carbonyls bonded to nitrogen (C(=O)-N).
        % The SMARTS [C](=O)[#6] requires the carbonyl carbon to have a
        % carbon neighbor -- caffeine has none, so it must return false.
        % This cross-validates that RDKit distinguishes amide vs ketone
        % carbonyls correctly and that emk wrappers propagate the result.
        % Reference: PubChem CID 2519 (Caffeine, C8H10N4O2).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
            tf  = emk.mol.hasSubstruct(mol, "[C](=O)[#6]");
            tc.verifyClass(tf, "logical", ...
                "hasSubstruct must return logical");
            tc.verifySize(tf, [1 1], ...
                "hasSubstruct must return a scalar [1 1]");
            tc.verifyFalse(tf, ...
                "Caffeine has only amide carbonyls (C=O bonded to N); " + ...
                "[C](=O)[#6] ketone SMARTS must return false");
        end

        function test_hasSubstruct_amideInCaffeine_returnsTrue(tc)
        % Caffeine's C=O groups are both bonded to nitrogen (amide/imide).
        % RDKit perceives the xanthine ring as aromatic, so the carbonyl
        % carbons are aromatic atoms.  SMARTS [#6](=O)[#7] matches any
        % carbon (aromatic or aliphatic) with =O bonded to any nitrogen.
        % Reference: PubChem CID 2519 (Caffeine, C8H10N4O2).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
            tf  = emk.mol.hasSubstruct(mol, "[#6](=O)[#7]");
            tc.verifyClass(tf, "logical", ...
                "hasSubstruct must return logical");
            tc.verifySize(tf, [1 1], ...
                "hasSubstruct must return a scalar [1 1]");
            tc.verifyTrue(tf, ...
                "Caffeine must contain a carbonyl bonded to nitrogen [#6](=O)[#7]");
        end

        % ------------------------------------------------------------------
        % TC14: hasSubstruct -- Mol object as query (requires RDKit)
        % Principle 1: no py.rdkit.* direct calls.
        % ------------------------------------------------------------------

        function test_hasSubstruct_molObjectQuery_benzene(tc)
        % Query as a Mol object (from emk.mol.fromSmiles) must work identically
        % to a SMARTS string query for simple cases.
        % Benzyl alcohol contains benzene (c1ccccc1) as a substructure.
        % verifySize added to confirm scalar contract (principle 2).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("c1ccccc1CO");
            query = emk.mol.fromSmiles("c1ccccc1");
            tf    = emk.mol.hasSubstruct(mol, query);
            tc.verifyClass(tf, "logical", ...
                "hasSubstruct with Mol query must return logical");
            tc.verifySize(tf, [1 1], ...
                "hasSubstruct with Mol query must return scalar [1 1]");
            tc.verifyTrue(tf, ...
                "Benzyl alcohol must contain benzene as Mol object query");
        end

        % ------------------------------------------------------------------
        % TC15: hasSubstruct -- Cell array of mols (requires RDKit)
        % Principle 2: verifyClass + verifySize confirm array contract.
        % ------------------------------------------------------------------

        function test_hasSubstruct_cellInput_returnsLogicalRowVector(tc)
        % Cell array input must return a logical row vector of matching length.
        % Input: 3 mols. Query: benzene ring.
        %   1. Ethanol (CCO)        -- no ring -> false
        %   2. Benzene (c1ccccc1)   -- is benzene -> true
        %   3. Benzyl alcohol       -- contains ring -> true
        % Expected: [false, true, true].
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mols = { emk.mol.fromSmiles("CCO"), ...
                     emk.mol.fromSmiles("c1ccccc1"), ...
                     emk.mol.fromSmiles("c1ccccc1CO") };
            tf = emk.mol.hasSubstruct(mols, "c1ccccc1");
            tc.verifyClass(tf, "logical", ...
                "hasSubstruct on cell must return logical array");
            tc.verifySize(tf, [1 3], ...
                "hasSubstruct on 3-element cell must return 1x3 array");
            tc.verifyEqual(tf, [false, true, true], ...
                "Expected [false, true, true] for benzene ring query");
        end

        function test_hasSubstruct_cellInput_allMatch(tc)
        % When all mols match the query the result must be all-true.
        % Query: hydroxyl group [OH]; mols: ethanol + benzyl alcohol (both have OH).
        % verifyClass + verifySize added for type contract (principle 2).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mols = { emk.mol.fromSmiles("CCO"), ...
                     emk.mol.fromSmiles("c1ccccc1CO") };
            tf = emk.mol.hasSubstruct(mols, "[OH]");
            tc.verifyClass(tf, "logical", ...
                "hasSubstruct all-match must return logical");
            tc.verifySize(tf, [1 2], ...
                "hasSubstruct on 2-element cell must return [1 2]");
            tc.verifyEqual(tf, [true, true], ...
                "Both ethanol and benzyl alcohol must match [OH]");
        end

        function test_hasSubstruct_cellInput_noneMatch(tc)
        % When no mols match the query the result must be all-false.
        % Query: amine [NH2]; mols: ethanol + benzene (neither has amine).
        % verifyClass + verifySize added for type contract (principle 2).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mols = { emk.mol.fromSmiles("CCO"), ...
                     emk.mol.fromSmiles("c1ccccc1") };
            tf = emk.mol.hasSubstruct(mols, "[NH2]");
            tc.verifyClass(tf, "logical", ...
                "hasSubstruct none-match must return logical");
            tc.verifySize(tf, [1 2], ...
                "hasSubstruct on 2-element cell must return [1 2]");
            tc.verifyEqual(tf, [false, false], ...
                "Neither ethanol nor benzene must match [NH2]");
        end

        function test_hasSubstruct_singleElementCell_returnsScalarSizedArray(tc)
        % A cell array with one element must return a 1x1 logical array
        % (not a scalar -- the cell path always returns a row vector).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mols = { emk.mol.fromSmiles("CCO") };
            tf = emk.mol.hasSubstruct(mols, "[OH]");
            tc.verifyClass(tf, "logical", ...
                "Single-element cell must return logical");
            tc.verifySize(tf, [1 1], ...
                "Single-element cell must return [1 1] array");
            tc.verifyTrue(tf, ...
                "Ethanol in a cell must match [OH]");
        end

        function test_hasSubstruct_charQuery_accepted(tc)
        % char literal SMARTS query must be accepted (ischar || isStringScalar guard).
        % verifyClass + verifySize confirm type contract also holds for char input.
        % Reference: PubChem CID 702 (Ethanol).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            tf  = emk.mol.hasSubstruct(mol, '[OH]');   % char literal
            tc.verifyClass(tf, "logical", ...
                "hasSubstruct must return logical for char query");
            tc.verifySize(tf, [1 1], ...
                "hasSubstruct must return scalar [1 1] for char query");
            tc.verifyTrue(tf, "char SMARTS query must be accepted for ethanol [OH]");
        end

        % ------------------------------------------------------------------
        % TC16: hasSubstruct -- Cross-validation (requires RDKit)
        %
        % Principle 3: compare emk.* wrapper result with Python's
        % mol.HasSubstructMatch(query) direct call.
        % mol.HasSubstructMatch is called directly because no emk.* wrapper
        % for raw HasSubstructMatch exists (mirror: mol.GetNumBonds() in TC5).
        % ------------------------------------------------------------------

        function test_hasSubstruct_crossValidate_trueCase(tc)
        % Cross-validation [true case]: emk.mol.hasSubstruct must equal
        % logical(mol.HasSubstructMatch(query)) when the match exists.
        % Reference: PubChem CID 244 (Benzyl alcohol) contains CID 241 (Benzene).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("c1ccccc1CO");  % benzyl alcohol
            query = emk.mol.fromSmiles("c1ccccc1");    % benzene
            % emk.* wrapper result
            wrapperResult = emk.mol.hasSubstruct(mol, query);
            % Python cross-validation: direct method call on mol (no emk wrapper)
            pythonResult  = logical(mol.HasSubstructMatch(query));
            tc.verifyEqual(wrapperResult, pythonResult, ...
                "hasSubstruct must match Python HasSubstructMatch [true case]");
            tc.verifyTrue(wrapperResult, ...
                "Benzyl alcohol must contain benzene substructure");
        end

        function test_hasSubstruct_crossValidate_falseCase(tc)
        % Cross-validation [false case]: emk.mol.hasSubstruct must equal
        % logical(mol.HasSubstructMatch(query)) when no match exists.
        % Reference: PubChem CID 702 (Ethanol) does not contain CID 241 (Benzene).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CCO");         % ethanol
            query = emk.mol.fromSmiles("c1ccccc1");    % benzene
            % emk.* wrapper result
            wrapperResult = emk.mol.hasSubstruct(mol, query);
            % Python cross-validation: direct method call on mol (no emk wrapper)
            pythonResult  = logical(mol.HasSubstructMatch(query));
            tc.verifyEqual(wrapperResult, pythonResult, ...
                "hasSubstruct must match Python HasSubstructMatch [false case]");
            tc.verifyFalse(wrapperResult, ...
                "Ethanol must NOT contain benzene substructure");
        end

    end

    % ======================================================================
    % scaffold tests
    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % SC-TC1: Input validation -- no RDKit required
        % ------------------------------------------------------------------

        function test_scaffold_numericInput_throwsInvalidInput(tc)
        % Numeric input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.mol.scaffold(42), ...
                "emk:mol:scaffold:invalidInput", ...
                "Numeric input must throw invalidInput");
        end

        function test_scaffold_numericInput_errorMessage_containsClass(tc)
        % Error message must contain the actual class name for debugging.
            ME = tc.captureError(@() emk.mol.scaffold(42));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain the offending class name");
        end

        function test_scaffold_stringInput_throwsInvalidInput(tc)
        % Passing a SMILES string (not a Mol object) must throw invalidInput.
            tc.verifyError(@() emk.mol.scaffold("c1ccccc1"), ...
                "emk:mol:scaffold:invalidInput", ...
                "String input must throw invalidInput");
        end

        function test_scaffold_charInput_throwsInvalidInput(tc)
        % char literal input must throw invalidInput.
            tc.verifyError(@() emk.mol.scaffold('c1ccccc1'), ...
                "emk:mol:scaffold:invalidInput", ...
                "char input must throw invalidInput");
        end

        function test_scaffold_logicalInput_throwsInvalidInput(tc)
        % logical input must throw invalidInput.
            tc.verifyError(@() emk.mol.scaffold(true), ...
                "emk:mol:scaffold:invalidInput", ...
                "logical input must throw invalidInput");
        end

        function test_scaffold_emptyInput_throwsInvalidInput(tc)
        % Empty matrix input must throw invalidInput.
            tc.verifyError(@() emk.mol.scaffold([]), ...
                "emk:mol:scaffold:invalidInput", ...
                "Empty matrix input must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % SC-TC2: Return type (requires RDKit)
        % ------------------------------------------------------------------

        function test_scaffold_benzene_returnsPythonObject(tc)
        % scaffold must return a Python Mol object.
        % Uses startsWith/contains to remain stable across RDKit versions.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccccc1");
            scaf = emk.mol.scaffold(mol);
            tc.verifyTrue(startsWith(class(scaf), "py."), ...
                "scaffold must return a Python object (class starts with 'py.')");
            tc.verifyTrue(contains(class(scaf), "Mol"), ...
                sprintf("scaffold must return a Mol type, got: %s", class(scaf)));
        end

        % ------------------------------------------------------------------
        % SC-TC3: Benzene scaffold == benzene (pure ring => identity)
        % Reference: PubChem CID 241 (Benzene, c1ccccc1, 6 heavy atoms)
        % Murcko scaffold of a single ring system is the ring itself.
        % ------------------------------------------------------------------

        function test_scaffold_benzene_hasSixAtoms(tc)
        % Benzene scaffold has 6 heavy atoms (the ring = the scaffold).
        % HeavyAtomCount via emk.descriptor.calculate avoids direct py.* call.
        % Cross-validation: canonical SMILES round-trip also confirms identity.
        % Reference: PubChem CID 241.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccccc1");
            scaf = emk.mol.scaffold(mol);
            desc = emk.descriptor.calculate(scaf, "HeavyAtomCount");
            tc.verifyClass(desc.HeavyAtomCount, "double", ...
                "scaffold HeavyAtomCount must be double");
            tc.verifyEqual(desc.HeavyAtomCount, 6, ...
                "Benzene scaffold must have 6 heavy atoms");
        end

        function test_scaffold_benzene_ringCount_isOne(tc)
        % Benzene scaffold must have RingCount = 1 (the benzene ring).
        % Validates that the ring structure is preserved in the scaffold.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccccc1");
            scaf = emk.mol.scaffold(mol);
            desc = emk.descriptor.calculate(scaf, "RingCount");
            tc.verifyEqual(desc.RingCount, 1, ...
                "Benzene scaffold must have RingCount = 1");
        end

        % ------------------------------------------------------------------
        % SC-TC4: Aspirin scaffold is benzene ring (side chains stripped)
        % Aspirin: CC(=O)Oc1ccccc1C(=O)O (PubChem CID 2244, 13 heavy atoms)
        % Murcko scaffold = benzene (side chains -C(=O)CH3 and -C(=O)O removed)
        % ------------------------------------------------------------------

        function test_scaffold_aspirin_fewerAtomsThanParent(tc)
        % Aspirin scaffold must have fewer heavy atoms than aspirin itself.
        % Aspirin has 13 heavy atoms; the benzene scaffold has 6.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol           = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            scaf          = emk.mol.scaffold(mol);
            descParent    = emk.descriptor.calculate(mol,  "HeavyAtomCount");
            descScaffold  = emk.descriptor.calculate(scaf, "HeavyAtomCount");
            tc.verifyLessThan(descScaffold.HeavyAtomCount, descParent.HeavyAtomCount, ...
                "Aspirin scaffold must have fewer atoms than aspirin (side chains stripped)");
        end

        function test_scaffold_aspirin_ringCount_isOne(tc)
        % Aspirin scaffold (benzene ring) must have RingCount = 1.
        % Cross-validation: confirms ring structure is preserved (principle 3).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            scaf = emk.mol.scaffold(mol);
            desc = emk.descriptor.calculate(scaf, "RingCount");
            tc.verifyEqual(desc.RingCount, 1, ...
                "Aspirin scaffold (benzene) must have RingCount = 1");
        end

        function test_scaffold_aspirin_smiles_isBenzene(tc)
        % Aspirin scaffold SMILES must be canonical benzene "c1ccccc1".
        % This is the primary cross-validation test for scaffold correctness
        % (principle 3): SMILES of the scaffold must match the known Murcko
        % scaffold of aspirin (benzene ring with side chains removed).
        % Reference: PubChem CID 2244 (aspirin) -> scaffold = benzene.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            scaf = emk.mol.scaffold(mol);
            smi  = emk.mol.toSmiles(scaf);
            tc.verifyEqual(smi, "c1ccccc1", ...
                "Aspirin Murcko scaffold SMILES must be c1ccccc1 (benzene)");
        end

        % ------------------------------------------------------------------
        % SC-TC5: Acyclic molecule => empty scaffold (0 atoms)
        % Ethanol (CCO) has no rings; Murcko scaffold is the empty molecule.
        % This is the canonical RDKit behaviour for acyclic inputs.
        % Reference: PubChem CID 702 (Ethanol)
        % ------------------------------------------------------------------

        function test_scaffold_ethanol_returnsEmptyScaffold(tc)
        % Ethanol has no rings; its Murcko scaffold must have 0 atoms.
        % The return value is a valid (non-None) Mol Python object with 0 atoms.
        % Cross-validation: HeavyAtomCount via emk.descriptor.calculate confirms 0.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            scaf = emk.mol.scaffold(mol);
            tc.verifyTrue(startsWith(class(scaf), "py."), ...
                "Empty scaffold must still be a Python object (not None)");
            desc = emk.descriptor.calculate(scaf, "HeavyAtomCount");
            tc.verifyClass(desc.HeavyAtomCount, "double", ...
                "HeavyAtomCount must be double");
            tc.verifyEqual(desc.HeavyAtomCount, 0, ...
                "Ethanol Murcko scaffold must have 0 heavy atoms (acyclic molecule)");
        end

        % ------------------------------------------------------------------
        % SC-TC6: Idempotency -- scaffold(scaffold(mol)) == scaffold(mol)
        % The Murcko scaffold is a fixed-point operation: applying scaffold
        % to a scaffold must return the same SMILES.
        % ------------------------------------------------------------------

        function test_scaffold_idempotency_aspirin(tc)
        % scaffold(scaffold(aspirin)) must equal scaffold(aspirin).
        % Idempotency is a fundamental property of scaffold extraction:
        % a scaffold of a scaffold is the same scaffold.
        % Reference: aspirin scaffold = benzene (6 atoms, 1 ring).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            scaf1 = emk.mol.scaffold(mol);
            scaf2 = emk.mol.scaffold(scaf1);
            smi1  = emk.mol.toSmiles(scaf1);
            smi2  = emk.mol.toSmiles(scaf2);
            tc.verifyEqual(smi1, smi2, ...
                "scaffold(scaffold(mol)) must equal scaffold(mol) (idempotency)");
        end

        function test_scaffold_idempotency_benzene(tc)
        % scaffold(benzene) must equal benzene (pure ring => identity).
        % Benzene is already its own Murcko scaffold (no side chains).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccccc1");
            scaf = emk.mol.scaffold(mol);
            smi1 = emk.mol.toSmiles(mol);
            smi2 = emk.mol.toSmiles(scaf);
            tc.verifyEqual(smi1, smi2, ...
                "scaffold(benzene) must equal benzene (ring is its own scaffold)");
        end

        % ------------------------------------------------------------------
        % SC-TC7: Heteroatom preservation -- Bemis-Murcko retains ring
        % heteroatoms in the framework (ring N/O/S atoms are NOT stripped).
        % Reference: Bemis & Murcko (1996) §2: framework includes ring atoms.
        % ------------------------------------------------------------------

        function test_scaffold_pyridine_smiles_containsNitrogen(tc)
        % Pyridine scaffold SMILES must contain 'n' (ring nitrogen preserved).
        % Bemis-Murcko framework retains ring heteroatoms verbatim.
        % If the nitrogen were stripped, GetScaffoldForMol would return benzene
        % ('c1ccccc1') -- this test guards against that regression.
        % Reference: PubChem CID 1049 (pyridine, c1ccncc1, 6 heavy atoms).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccncc1");   % pyridine
            scaf = emk.mol.scaffold(mol);
            smi  = emk.mol.toSmiles(scaf);
            tc.verifySubstring(smi, "n", ...
                "Pyridine scaffold SMILES must contain 'n' (heteroatom preserved)");
        end

        function test_scaffold_pyridine_hasSixAtoms(tc)
        % Pyridine scaffold (= pyridine itself) has 6 heavy atoms (5C + 1N).
        % Cross-validation: HeavyAtomCount verifies that no atoms were lost
        % when extracting the scaffold of a pure heteroaromatic ring.
        % Reference: PubChem CID 1049 (pyridine).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccncc1");
            scaf = emk.mol.scaffold(mol);
            desc = emk.descriptor.calculate(scaf, "HeavyAtomCount");
            tc.verifyEqual(desc.HeavyAtomCount, 6, ...
                "Pyridine scaffold must have 6 heavy atoms (5C + 1N)");
        end

        % ------------------------------------------------------------------
        % SC-TC8: Fused multi-ring system preserves all rings
        % Naphthalene has 2 fused rings; scaffold = naphthalene (no side chains).
        % RingCount = 2 confirms both rings are retained in the scaffold.
        % ------------------------------------------------------------------

        function test_scaffold_naphthalene_ringCountIsTwo(tc)
        % Naphthalene scaffold must have RingCount = 2.
        % Naphthalene is a fused bicyclic compound; since it has no side chains,
        % its scaffold is itself. This test ensures that fused ring systems are
        % preserved intact in the scaffold.
        % Reference: PubChem CID 931 (naphthalene, c1ccc2ccccc2c1, 10 heavy atoms).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccc2ccccc2c1");   % naphthalene
            scaf = emk.mol.scaffold(mol);
            desc = emk.descriptor.calculate(scaf, "RingCount");
            tc.verifyEqual(desc.RingCount, 2, ...
                "Naphthalene scaffold must have RingCount = 2 (fused bicyclic)");
        end

        function test_scaffold_naphthalene_hasTenAtoms(tc)
        % Naphthalene scaffold (= naphthalene itself) has 10 heavy atoms.
        % Cross-validation: confirms the fused scaffold has no atom loss.
        % Reference: PubChem CID 931 (naphthalene, 10 C atoms).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccc2ccccc2c1");
            scaf = emk.mol.scaffold(mol);
            desc = emk.descriptor.calculate(scaf, "HeavyAtomCount");
            tc.verifyEqual(desc.HeavyAtomCount, 10, ...
                "Naphthalene scaffold must have 10 heavy atoms");
        end

        % ------------------------------------------------------------------
        % SC-TC9: Linker atom preservation
        % Bemis-Murcko retains atoms that LINK two ring systems (linker atoms).
        % Unlike side chains (connected to only one ring), linker atoms are
        % on the path between two ring systems and are part of the framework.
        %
        % Test molecule: 4-methylbibenzyl (Cc1ccc(CCc2ccccc2)cc1)
        %   Side chain: CH3 on ring 1 (removed)
        %   Linker:     CH2-CH2 connecting ring 1 and ring 2 (RETAINED)
        %   Expected scaffold: c1ccc(CCc2ccccc2)cc1 (bibenzyl, 14 heavy atoms)
        %
        % Reference: Bemis & Murcko (1996) §2.1: "Linker atoms connect ring
        %   systems and are retained in the molecular framework."
        %   PubChem CID 241 (benzene) + PubChem CID 12298 (bibenzyl).
        % ------------------------------------------------------------------

        function test_scaffold_linker_isPreserved(tc)
        % 4-methylbibenzyl scaffold must be bibenzyl (linker CH2CH2 preserved).
        % Side chain CH3 is stripped; linker CH2-CH2 between the two phenyl
        % rings is retained. This is the core Bemis-Murcko framework property.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("Cc1ccc(CCc2ccccc2)cc1");
            scaf = emk.mol.scaffold(mol);
            % RingCount = 2: both phenyl rings are retained
            descRing = emk.descriptor.calculate(scaf, "RingCount");
            tc.verifyEqual(descRing.RingCount, 2, ...
                "4-methylbibenzyl scaffold must have RingCount = 2 (two phenyl rings)");
        end

        function test_scaffold_linker_atomCount_is14(tc)
        % 4-methylbibenzyl scaffold (bibenzyl) must have 14 heavy atoms.
        % Parent has 15 heavy atoms (6+6 ring C + 2 linker C + 1 side-chain CH3).
        % After stripping CH3: 14 atoms (6+6 ring C + 2 linker C).
        % Cross-validation: confirms linker CH2CH2 is preserved (not 12 atoms).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("Cc1ccc(CCc2ccccc2)cc1");
            scaf = emk.mol.scaffold(mol);
            desc = emk.descriptor.calculate(scaf, "HeavyAtomCount");
            tc.verifyEqual(desc.HeavyAtomCount, 14, ...
                "4-methylbibenzyl scaffold must have 14 heavy atoms (linker preserved)");
        end

        function test_scaffold_linker_fewerAtomsThanParent(tc)
        % 4-methylbibenzyl scaffold must have fewer atoms than parent (side chain stripped).
        % Parent: 15 heavy atoms.  Scaffold (bibenzyl): 14 heavy atoms.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol          = emk.mol.fromSmiles("Cc1ccc(CCc2ccccc2)cc1");
            scaf         = emk.mol.scaffold(mol);
            descParent   = emk.descriptor.calculate(mol,  "HeavyAtomCount");
            descScaffold = emk.descriptor.calculate(scaf, "HeavyAtomCount");
            tc.verifyLessThan(descScaffold.HeavyAtomCount, descParent.HeavyAtomCount, ...
                "Scaffold must have fewer atoms than parent (side chain removed)");
        end

    end

    % ======================================================================
    methods (Access = private)

        function tf = rdkitAvailable(~)
        % Return true if Python is configured and RDKit can be imported.
        %
        % In MATLAB OutOfProcess mode, pyenv().Status stays "NotLoaded" until
        % the first Python API call (the process starts on demand).  Checking
        % Status and returning early would always produce false before first use.
        % Instead, ensure pyenv is configured via initPython() when NotLoaded,
        % then attempt the import directly -- matching verify() behaviour.
            tf = false;
            try
                pe = pyenv();
                if strcmp(string(pe.Status), "NotLoaded")
                    emk.setup.initPython();
                end
                py.importlib.import_module("rdkit.Chem");
                tf = true;
            catch
                % Python not configured or RDKit not installed
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
