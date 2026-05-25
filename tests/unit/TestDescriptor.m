classdef TestDescriptor < matlab.unittest.TestCase
% TestDescriptor  Unit tests for emk.descriptor.molWeight and
%                 emk.descriptor.calculate.
%
% Run with:
%   addpath(genpath("src"));
%   results = run(TestDescriptor);
%
% Coverage (molWeight):
%   TC1: Non-Mol input => invalidInput error (no RDKit required)
%         + error message contains offending class name for double AND string
%         + char, cell, string, logical, and empty inputs all rejected
%   TC2: Valid Mol => MW value matches reference (requires RDKit)
%         + ethanol    (CCO,                    C2H6O):   46.069 g/mol
%         + benzene    (c1ccccc1,               C6H6):    78.114 g/mol
%         + aspirin    (CC(=O)Oc1ccccc1C(=O)O, C9H8O4): 180.159 g/mol
%         + water      (O,                      H2O):     18.015 g/mol
%         AbsTol = 0.01 g/mol (rationale: docs/algorithm_guide.md 4.1)
%   TC3: Output type and shape: double scalar, real-valued (requires RDKit)
%   TC4: MW is strictly positive for valid molecules (requires RDKit)
%   TC5: MW is finite (isfinite) -- guards against Inf output which
%         isPositive alone cannot detect (Inf > 0 == true in MATLAB)
%         + MW is real-valued (isreal) -- guards against complex output
%
% Coverage (calculate):
%   TC6: Non-Mol input => invalidInput error (no RDKit required)
%         + char, string, numeric, logical, empty all rejected
%   TC7: Unknown descriptor name => unknownDescriptor error (no RDKit required)
%         + error message contains offending name
%         + partial unknown (mix of valid and invalid) also rejected
%
% Coverage additions (batchCalculate -- M3-1 error handling):
%   TB10: All non-Mol elements => allMolsFailed error (no RDKit required)
%         + error message contains molecule count
%         + single-element and multi-element cells both trigger the error
%   TC8: Full calculate (no names arg) => struct with all 10 fields (requires RDKit)
%   TC9: Descriptor values for ethanol (requires RDKit)
%         + MolWt = 46.069 +/- 0.01, TPSA = 20.23 +/- 0.1
%         + NumHAcceptors=1, NumHDonors=1, NumRotatableBonds=0
%         + RingCount=0, FractionCSP3=1.0 +/- 0.001, HeavyAtomCount=3
%   TC10: Subset calculate returns only requested fields (requires RDKit)
%   TC11: All values are finite double scalars (requires RDKit)
%   TC12: Aspirin spot checks: HeavyAtomCount=13, NumHDonors=1, RingCount=1
%
% Coverage additions (calculate -- gaps addressed):
%   TC6+: cell input rejected; invalidInput message contains class name
%   TC7+: error message mentions "Supported" (hints valid names to user)
%   TC13: char descriptor name is accepted (flexibility of input type)
%   TC14: ExactMolWt reference value for ethanol (~46.042 g/mol, AbsTol=0.002)
%         LogP sign for ethanol (negative, hydrophilic) and benzene (positive)
%         LogP bounds for ethanol: in range (-2, 0)
%   TC15: FractionCSP3=0.0 for benzene (all sp2 carbons -- complement to Fsp3=1)
%         NumRotatableBonds=1 for n-butane (non-zero -- complement to 0 for ethanol)
%         RingCount=1 for benzene (pure aromatic; cross-validates aspirin result)
%   TC16: Aspirin NumHAcceptors=3 (acetyl C=O, ester O, carboxyl C=O)
%   TC17: calculate(mol,"MolWt") == molWeight(mol) consistency check
%
% Coverage additions (TC18 -- Python-side cross-validation):
%   TC18: HeavyAtomCount == double(mol.GetNumAtoms()) cross-validation
%         Integer-value check for count descriptors
%         (NumHAcceptors, NumHDonors, NumRotatableBonds, RingCount, HeavyAtomCount)
%
% Coverage additions (TC19 -- MolFormula string descriptor):
%   TC19a: "MolFormula" is accepted by validation (not unknownDescriptor) [no RDKit]
%          Default calculate() does NOT include MolFormula (string/double mix risk)
%   TC19b: MolFormula returns string scalar, non-empty (requires RDKit)
%   TC19c: Reference values: ethanol="C2H6O", aspirin="C9H8O4", benzene="C6H6"
%   TC19d: MolFormula combinable with numeric descriptors in one call
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
        % TC1: Input validation -- no RDKit required
        % ------------------------------------------------------------------

        function test_molWeight_numericInput_throwsInvalidInput(tc)
        % Numeric input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.descriptor.molWeight(42), ...
                "emk:descriptor:molWeight:invalidInput", ...
                "Numeric input must throw invalidInput");
        end

        function test_molWeight_numericInput_errorMessage_containsClass(tc)
        % Error message must contain the actual class name for debugging.
            ME = tc.captureError(@() emk.descriptor.molWeight(42));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain the offending class name");
        end

        function test_molWeight_stringInput_throwsInvalidInput(tc)
        % Passing a raw SMILES string (not a Mol object) must throw invalidInput.
            tc.verifyError(@() emk.descriptor.molWeight("CCO"), ...
                "emk:descriptor:molWeight:invalidInput", ...
                "String input must throw invalidInput");
        end

        function test_molWeight_stringInput_errorMessage_containsClass(tc)
        % Error message for string input must contain "string" class name.
        % Mirrors the analogous double test: the error contract applies to
        % every non-Mol type, not just double.
            ME = tc.captureError(@() emk.descriptor.molWeight("CCO"));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "string", ...
                "Error message must contain the offending class name 'string'");
        end

        function test_molWeight_charInput_throwsInvalidInput(tc)
        % char literal (not a Mol object) must throw invalidInput.
        % The isa() guard rejects all non-Mol types including char.
        % Mirrors test_toSmiles_charInput_throwsInvalidInput in TestMol.m.
            tc.verifyError(@() emk.descriptor.molWeight('CCO'), ...
                "emk:descriptor:molWeight:invalidInput", ...
                "char input must throw invalidInput");
        end

        function test_molWeight_cellInput_throwsInvalidInput(tc)
        % Cell array input must throw invalidInput.
        % Mirrors test_fromSmiles_cellInput_throwsInvalidInput in TestMol.m.
            tc.verifyError(@() emk.descriptor.molWeight({"CCO"}), ...
                "emk:descriptor:molWeight:invalidInput", ...
                "Cell array input must throw invalidInput");
        end

        function test_molWeight_emptyInput_throwsInvalidInput(tc)
        % Empty matrix input must throw invalidInput.
            tc.verifyError(@() emk.descriptor.molWeight([]), ...
                "emk:descriptor:molWeight:invalidInput", ...
                "Empty matrix input must throw invalidInput");
        end

        function test_molWeight_logicalInput_throwsInvalidInput(tc)
        % logical input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.descriptor.molWeight(true), ...
                "emk:descriptor:molWeight:invalidInput", ...
                "logical input must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC2: MW values for reference molecules (requires RDKit)
        %
        % AbsTol = 0.01 g/mol:
        %   RDKit uses IUPAC 2021 atomic weights. 0.01 g/mol tolerance
        %   accommodates rounding in the last digit of reported PubChem
        %   values and minor IUPAC update cycles without masking gross
        %   errors.  See docs/algorithm_guide.md section 4.1 for full
        %   rationale and reference values.
        % ------------------------------------------------------------------

        function test_molWeight_ethanol_value(tc)
        % Ethanol (CCO, C2H6O) average MW = 46.069 g/mol.
        % Reference: PubChem CID 702 (MolecularWeight = 46.07 g/mol).
        % Computed: 2*12.011 + 6*1.008 + 15.999 = 46.069.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            mw  = emk.descriptor.molWeight(mol);
            tc.verifyEqual(mw, 46.069, "AbsTol", 0.01, ...
                "Ethanol MW must be 46.069 +/- 0.01 g/mol");
        end

        function test_molWeight_benzene_value(tc)
        % Benzene (c1ccccc1, C6H6) average MW = 78.114 g/mol.
        % Reference: PubChem CID 241 (MolecularWeight = 78.11 g/mol).
        % Computed: 6*12.011 + 6*1.008 = 78.114.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("c1ccccc1");
            mw  = emk.descriptor.molWeight(mol);
            tc.verifyEqual(mw, 78.114, "AbsTol", 0.01, ...
                "Benzene MW must be 78.114 +/- 0.01 g/mol");
        end

        function test_molWeight_aspirin_value(tc)
        % Aspirin (CC(=O)Oc1ccccc1C(=O)O, C9H8O4) average MW = 180.159 g/mol.
        % Reference: PubChem CID 2244 (MolecularWeight = 180.16 g/mol).
        % Computed: 9*12.011 + 8*1.008 + 4*15.999 = 180.159.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mw  = emk.descriptor.molWeight(mol);
            tc.verifyEqual(mw, 180.159, "AbsTol", 0.01, ...
                "Aspirin MW must be 180.159 +/- 0.01 g/mol");
        end

        function test_molWeight_water_value(tc)
        % Water (O, H2O) average MW = 18.015 g/mol.
        % Reference: PubChem CID 962 (MolecularWeight = 18.015 g/mol).
        % Edge case: single heavy atom with 2 implicit hydrogens.
        % Computed: 2*1.008 + 15.999 = 18.015.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("O");
            mw  = emk.descriptor.molWeight(mol);
            tc.verifyEqual(mw, 18.015, "AbsTol", 0.01, ...
                "Water MW must be 18.015 +/- 0.01 g/mol");
        end

        % ------------------------------------------------------------------
        % TC3: Output type and shape (requires RDKit)
        % ------------------------------------------------------------------

        function test_molWeight_outputIsDoubleScalar(tc)
        % molWeight must return a scalar double.
        % Verifies that the Python float is properly converted via double().
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            mw  = emk.descriptor.molWeight(mol);
            tc.verifyClass(mw, "double", ...
                "molWeight must return a double");
            tc.verifySize(mw, [1 1], ...
                "molWeight must return a scalar");
        end

        function test_molWeight_outputIsRealValued(tc)
        % molWeight must return a real-valued (non-complex) double.
        % Average molecular weight is a real physical quantity; complex
        % output would indicate an error in the Python-to-MATLAB conversion.
        % isreal() returns false for complex doubles even when imag part is 0.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            mw  = emk.descriptor.molWeight(mol);
            tc.verifyTrue(isreal(mw), ...
                "molWeight must return a real-valued double");
        end

        % ------------------------------------------------------------------
        % TC4: MW is strictly positive (requires RDKit)
        % ------------------------------------------------------------------

        function test_molWeight_ethanol_isPositive(tc)
        % MW must be strictly positive for a valid molecule.
        % Defensive test: guards against sign inversion or zero output.
        % Note: NaN > 0 == false in MATLAB, so this also catches NaN output.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            mw  = emk.descriptor.molWeight(mol);
            tc.verifyGreaterThan(mw, 0, "MW must be strictly positive");
        end

        % ------------------------------------------------------------------
        % TC5: MW is finite -- isPositive alone cannot detect Inf
        %      (Inf > 0 == true in MATLAB, so isPositive passes for Inf)
        % ------------------------------------------------------------------

        function test_molWeight_ethanol_isFinite(tc)
        % MW must be finite (not Inf or -Inf).
        % isPositive alone does NOT catch Inf because Inf > 0 == true.
        % isfinite() returns false for both Inf and NaN, making this a
        % stricter guard that subsumes the NaN case covered by isPositive.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CCO");
            mw  = emk.descriptor.molWeight(mol);
            tc.verifyTrue(isfinite(mw), ...
                "MW must be finite (not Inf or NaN)");
        end

        function test_molWeight_aspirin_isFinite(tc)
        % Aspirin MW must also be finite.
        % Uses a structurally complex molecule (aromatic + multiple
        % functional groups) to confirm finiteness is not molecule-specific.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mw  = emk.descriptor.molWeight(mol);
            tc.verifyTrue(isfinite(mw), ...
                "Aspirin MW must be finite (not Inf or NaN)");
        end

    end

    % ======================================================================
    % calculate() tests
    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC6: Input validation -- no RDKit required
        % ------------------------------------------------------------------

        function test_calculate_numericInput_throwsInvalidInput(tc)
        % Numeric input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.descriptor.calculate(42), ...
                "emk:descriptor:calculate:invalidInput", ...
                "Numeric input must throw invalidInput");
        end

        function test_calculate_stringInput_throwsInvalidInput(tc)
        % Passing a raw SMILES string must throw invalidInput.
            tc.verifyError(@() emk.descriptor.calculate("CCO"), ...
                "emk:descriptor:calculate:invalidInput", ...
                "String input must throw invalidInput");
        end

        function test_calculate_charInput_throwsInvalidInput(tc)
        % char input must throw invalidInput.
            tc.verifyError(@() emk.descriptor.calculate('CCO'), ...
                "emk:descriptor:calculate:invalidInput", ...
                "char input must throw invalidInput");
        end

        function test_calculate_logicalInput_throwsInvalidInput(tc)
        % logical input must throw invalidInput.
            tc.verifyError(@() emk.descriptor.calculate(true), ...
                "emk:descriptor:calculate:invalidInput", ...
                "logical input must throw invalidInput");
        end

        function test_calculate_emptyInput_throwsInvalidInput(tc)
        % Empty matrix input must throw invalidInput.
            tc.verifyError(@() emk.descriptor.calculate([]), ...
                "emk:descriptor:calculate:invalidInput", ...
                "Empty matrix input must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TC7: Unknown descriptor names -- no RDKit required
        % ------------------------------------------------------------------

        function test_calculate_unknownName_throwsUnknownDescriptor(tc)
        % A completely unknown descriptor name must throw unknownDescriptor.
            tc.verifyError( ...
                @() emk.descriptor.calculate(struct(), "NotADescriptor"), ...
                "emk:descriptor:calculate:unknownDescriptor", ...
                "Unknown name must throw unknownDescriptor");
        end

        function test_calculate_unknownName_errorMessage_containsName(tc)
        % Error message must contain the offending descriptor name.
            ME = tc.captureError( ...
                @() emk.descriptor.calculate(struct(), "BadName"));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "BadName", ...
                "Error message must contain the unknown descriptor name");
        end

        function test_calculate_partialUnknown_throwsUnknownDescriptor(tc)
        % Mix of valid and invalid names must still throw unknownDescriptor.
        % The validation rejects the entire request when any name is unknown.
            tc.verifyError( ...
                @() emk.descriptor.calculate(struct(), ["MolWt","BadName"]), ...
                "emk:descriptor:calculate:unknownDescriptor", ...
                "Partial unknown must throw unknownDescriptor");
        end

        % ------------------------------------------------------------------
        % TC8: Full calculate returns struct with all 10 fields (requires RDKit)
        % ------------------------------------------------------------------

        function test_calculate_allDescriptors_returnsStruct(tc)
        % calculate with no names argument must return a struct.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol);
            tc.verifyClass(desc, "struct", ...
                "calculate must return a struct");
        end

        function test_calculate_allDescriptors_hasAllTenFields(tc)
        % Default call must produce exactly 10 fields (all supported descriptors).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol      = emk.mol.fromSmiles("CCO");
            desc     = emk.descriptor.calculate(mol);
            expected = ["MolWt","ExactMolWt","LogP","TPSA", ...
                        "NumHAcceptors","NumHDonors","NumRotatableBonds", ...
                        "RingCount","FractionCSP3","HeavyAtomCount"];
            for k = 1:numel(expected)
                tc.verifyTrue(isfield(desc, expected(k)), ...
                    "Missing field: " + expected(k));
            end
            tc.verifyEqual(numel(fieldnames(desc)), 10, ...
                "Struct must have exactly 10 fields");
        end

        % ------------------------------------------------------------------
        % TC9: Descriptor values for ethanol (requires RDKit)
        % Reference: PubChem CID 702 / RDKit 2024.3.x
        % Rationale for tolerances: docs/algorithm_guide.md section 4.2
        % ------------------------------------------------------------------

        function test_calculate_ethanol_MolWt(tc)
        % Ethanol average MW = 46.069 g/mol (cross-check with molWeight.m).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "MolWt");
            tc.verifyClass(desc.MolWt, "double", "MolWt must be a double scalar");
            tc.verifySize(desc.MolWt, [1 1], "MolWt must be a scalar [1 1]");
            tc.verifyEqual(desc.MolWt, 46.069, "AbsTol", 0.01, ...
                "Ethanol MolWt must be 46.069 +/- 0.01 g/mol");
        end

        function test_calculate_ethanol_TPSA(tc)
        % Ethanol TPSA = 20.23 Ang^2.
        % Reference: PubChem CID 702 (TPSA = 20.2 Ang^2).
        % RDKit Descriptors.TPSA uses Ertl (2000) contribution method.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "TPSA");
            tc.verifyClass(desc.TPSA, "double", "TPSA must be a double scalar");
            tc.verifySize(desc.TPSA, [1 1], "TPSA must be a scalar [1 1]");
            tc.verifyEqual(desc.TPSA, 20.23, "AbsTol", 0.1, ...
                "Ethanol TPSA must be 20.23 +/- 0.1 Ang^2");
        end

        function test_calculate_ethanol_NumHAcceptors(tc)
        % Ethanol has 1 H-bond acceptor (the oxygen, Ertl definition).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "NumHAcceptors");
            tc.verifyClass(desc.NumHAcceptors, "double", "NumHAcceptors must be a double scalar");
            tc.verifySize(desc.NumHAcceptors, [1 1], "NumHAcceptors must be a scalar [1 1]");
            tc.verifyEqual(desc.NumHAcceptors, 1, ...
                "Ethanol must have 1 H-bond acceptor");
        end

        function test_calculate_ethanol_NumHDonors(tc)
        % Ethanol has 1 H-bond donor (the OH group, Ertl definition).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "NumHDonors");
            tc.verifyClass(desc.NumHDonors, "double", "NumHDonors must be a double scalar");
            tc.verifySize(desc.NumHDonors, [1 1], "NumHDonors must be a scalar [1 1]");
            tc.verifyEqual(desc.NumHDonors, 1, ...
                "Ethanol must have 1 H-bond donor");
        end

        function test_calculate_ethanol_NumRotatableBonds(tc)
        % Ethanol has 0 rotatable bonds (strict definition: bonds to
        % terminal heavy atoms are excluded; both C-C and C-O have a
        % terminal atom on one end).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "NumRotatableBonds");
            tc.verifyClass(desc.NumRotatableBonds, "double", "NumRotatableBonds must be a double scalar");
            tc.verifySize(desc.NumRotatableBonds, [1 1], "NumRotatableBonds must be a scalar [1 1]");
            tc.verifyEqual(desc.NumRotatableBonds, 0, ...
                "Ethanol must have 0 rotatable bonds (strict definition)");
        end

        function test_calculate_ethanol_RingCount(tc)
        % Ethanol has no rings.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "RingCount");
            tc.verifyClass(desc.RingCount, "double", "RingCount must be a double scalar");
            tc.verifySize(desc.RingCount, [1 1], "RingCount must be a scalar [1 1]");
            tc.verifyEqual(desc.RingCount, 0, ...
                "Ethanol must have 0 rings");
        end

        function test_calculate_ethanol_FractionCSP3(tc)
        % Ethanol has two sp3 carbons out of two total carbons => Fsp3 = 1.0.
        % Reference: Lovering et al. (2009) J. Med. Chem. 52(21):6752-6756.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "FractionCSP3");
            tc.verifyClass(desc.FractionCSP3, "double", "FractionCSP3 must be a double scalar");
            tc.verifySize(desc.FractionCSP3, [1 1], "FractionCSP3 must be a scalar [1 1]");
            tc.verifyEqual(desc.FractionCSP3, 1.0, "AbsTol", 0.001, ...
                "Ethanol FractionCSP3 must be 1.0 (both carbons are sp3)");
        end

        function test_calculate_ethanol_HeavyAtomCount(tc)
        % Ethanol (CCO) has 3 heavy atoms: 2 C + 1 O.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "HeavyAtomCount");
            tc.verifyClass(desc.HeavyAtomCount, "double", "HeavyAtomCount must be a double scalar");
            tc.verifySize(desc.HeavyAtomCount, [1 1], "HeavyAtomCount must be a scalar [1 1]");
            tc.verifyEqual(desc.HeavyAtomCount, 3, ...
                "Ethanol must have 3 heavy atoms");
        end

        % ------------------------------------------------------------------
        % TC10: Subset calculate (requires RDKit)
        % ------------------------------------------------------------------

        function test_calculate_subset_returnsOnlyRequestedFields(tc)
        % Requesting ["MolWt","LogP"] must return a struct with exactly
        % those two fields and no others.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, ["MolWt","LogP"]);
            tc.verifyTrue(isfield(desc, "MolWt"), "MolWt field must be present");
            tc.verifyTrue(isfield(desc, "LogP"),  "LogP field must be present");
            tc.verifyEqual(numel(fieldnames(desc)), 2, ...
                "Subset calculate must return exactly 2 fields");
        end

        function test_calculate_singleDescriptor_returnsStruct(tc)
        % Requesting a single descriptor name (not an array) must still
        % return a struct.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "TPSA");
            tc.verifyClass(desc, "struct", ...
                "Single-descriptor call must return a struct");
            tc.verifyTrue(isfield(desc, "TPSA"), "TPSA field must be present");
        end

        % ------------------------------------------------------------------
        % TC11: All values are finite double scalars (requires RDKit)
        % ------------------------------------------------------------------

        function test_calculate_allDescriptors_allValuesAreDoubleScalar(tc)
        % Every field in the full descriptor struct must be a double scalar.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol    = emk.mol.fromSmiles("CCO");
            desc   = emk.descriptor.calculate(mol);
            fields = fieldnames(desc);
            for k = 1:numel(fields)
                v = desc.(fields{k});
                tc.verifyClass(v, "double", fields{k} + " must be double");
                tc.verifySize(v, [1 1],    fields{k} + " must be scalar");
            end
        end

        function test_calculate_allDescriptors_allValuesAreFinite(tc)
        % Every field value must be finite (not Inf or NaN).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol    = emk.mol.fromSmiles("CCO");
            desc   = emk.descriptor.calculate(mol);
            fields = fieldnames(desc);
            for k = 1:numel(fields)
                tc.verifyTrue(isfinite(desc.(fields{k})), ...
                    fields{k} + " must be finite");
            end
        end

        % ------------------------------------------------------------------
        % TC12: Aspirin spot checks (requires RDKit)
        % Reference: PubChem CID 2244
        % ------------------------------------------------------------------

        function test_calculate_aspirin_HeavyAtomCount(tc)
        % Aspirin (C9H8O4) has 13 heavy atoms: 9 C + 4 O.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            desc = emk.descriptor.calculate(mol, "HeavyAtomCount");
            tc.verifyClass(desc.HeavyAtomCount, "double", "HeavyAtomCount must be a double scalar");
            tc.verifySize(desc.HeavyAtomCount, [1 1], "HeavyAtomCount must be a scalar [1 1]");
            tc.verifyEqual(desc.HeavyAtomCount, 13, ...
                "Aspirin must have 13 heavy atoms");
        end

        function test_calculate_aspirin_NumHDonors(tc)
        % Aspirin has 1 H-bond donor: the carboxylic acid OH.
        % The ester oxygen is not a donor; carbonyl oxygens are acceptors only.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            desc = emk.descriptor.calculate(mol, "NumHDonors");
            tc.verifyClass(desc.NumHDonors, "double", "NumHDonors must be a double scalar");
            tc.verifySize(desc.NumHDonors, [1 1], "NumHDonors must be a scalar [1 1]");
            tc.verifyEqual(desc.NumHDonors, 1, ...
                "Aspirin must have 1 H-bond donor");
        end

        function test_calculate_aspirin_RingCount(tc)
        % Aspirin has 1 ring: the benzene ring.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            desc = emk.descriptor.calculate(mol, "RingCount");
            tc.verifyClass(desc.RingCount, "double", "RingCount must be a double scalar");
            tc.verifySize(desc.RingCount, [1 1], "RingCount must be a scalar [1 1]");
            tc.verifyEqual(desc.RingCount, 1, ...
                "Aspirin must have 1 ring");
        end

    end

    % ======================================================================
    % calculate() -- gap-fill tests (TC6+, TC7+, TC13-TC17)
    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC6+: cell input and error message content -- no RDKit required
        % ------------------------------------------------------------------

        function test_calculate_cellInput_throwsInvalidInput(tc)
        % Cell array input must throw invalidInput.
        % Mirrors test_molWeight_cellInput_throwsInvalidInput -- must be
        % symmetric across all descriptor functions.
            tc.verifyError(@() emk.descriptor.calculate({"CCO"}), ...
                "emk:descriptor:calculate:invalidInput", ...
                "Cell array input must throw invalidInput");
        end

        function test_calculate_invalidInput_errorMessage_containsClass(tc)
        % Error message must contain the class name of the invalid input.
        % Mirrors test_molWeight_numericInput_errorMessage_containsClass.
        % The message format is "mol must be ..., got: <class>".
            ME = tc.captureError(@() emk.descriptor.calculate(42));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "double", ...
                "Error message must contain the offending class name");
        end

        % ------------------------------------------------------------------
        % TC7+: supported names hint in error message -- no RDKit required
        % ------------------------------------------------------------------

        function test_calculate_unknownName_errorMessage_mentionsSupportedNames(tc)
        % Error message must contain the word "Supported" to guide the user.
        % The full message format is:
        %   "Unknown descriptor(s): [<bad>]. Supported: [<list>]"
        % Testing the "Supported" keyword ensures the help text is present.
            ME = tc.captureError( ...
                @() emk.descriptor.calculate(struct(), "BadName"));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "Supported", ...
                "Error message must mention supported descriptor names");
        end

        % ------------------------------------------------------------------
        % TC13: Descriptor name input flexibility -- requires RDKit
        % ------------------------------------------------------------------

        function test_calculate_charDescriptorName_accepted(tc)
        % Passing a char array as descriptorNames must be accepted.
        % calculate.m converts char via ischar() guard before ismember check.
        % Users may pass 'MolWt' instead of "MolWt" from older MATLAB code.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, 'MolWt');  % char, not string
            tc.verifyTrue(isfield(desc, "MolWt"), ...
                "char descriptor name must be accepted and produce MolWt field");
            tc.verifyClass(desc.MolWt, "double", ...
                "MolWt from char name must be double");
        end

        % ------------------------------------------------------------------
        % TC14: ExactMolWt and LogP reference values (requires RDKit)
        % ------------------------------------------------------------------

        function test_calculate_ethanol_ExactMolWt(tc)
        % Ethanol monoisotopic MW: 2*12.000 + 6*1.00783 + 15.99491 = 46.042.
        % RDKit ExactMolWt uses the most abundant isotope of each element
        % (IUPAC 2016 monoisotopic masses).
        % AbsTol = 0.002: monoisotopic mass is more tightly defined than avg MW.
        % Reference: PubChem CID 702 (ExactMassValue = 46.042).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "ExactMolWt");
            tc.verifyEqual(desc.ExactMolWt, 46.042, "AbsTol", 0.002, ...
                "Ethanol ExactMolWt must be 46.042 +/- 0.002 g/mol");
        end

        function test_calculate_ethanol_LogP_isNegative(tc)
        % Ethanol Wildman-Crippen LogP must be negative (hydrophilic molecule).
        % Experimental logP(octanol/water) = -0.31 (NIST Webbook).
        % A negative LogP is unambiguous for ethanol regardless of
        % contribution table version; this is a sign test, not a value test.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "LogP");
            tc.verifyLessThan(desc.LogP, 0, ...
                "Ethanol LogP must be negative (hydrophilic)");
        end

        function test_calculate_ethanol_LogP_isReasonable(tc)
        % Ethanol LogP must be in (-2, 0): negative but not extremely hydrophilic.
        % Lower bound -2 guards against sign-inverted or saturated outputs.
        % Combined with test_calculate_ethanol_LogP_isNegative, the pair
        % tests: -2 < LogP < 0.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "LogP");
            tc.verifyGreaterThan(desc.LogP, -2, ...
                "Ethanol LogP must be > -2 (not extremely hydrophilic)");
        end

        function test_calculate_benzene_LogP_isPositive(tc)
        % Benzene Wildman-Crippen LogP must be positive (lipophilic aromatic).
        % Experimental logP(octanol/water) = 2.13 (NIST Webbook).
        % Positive sign is unambiguous for benzene; complements ethanol test.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccccc1");
            desc = emk.descriptor.calculate(mol, "LogP");
            tc.verifyGreaterThan(desc.LogP, 0, ...
                "Benzene LogP must be positive (lipophilic)");
        end

        % ------------------------------------------------------------------
        % TC15: Cross-molecule validation (requires RDKit)
        % ------------------------------------------------------------------

        function test_calculate_benzene_FractionCSP3(tc)
        % Benzene has 6 aromatic (sp2) carbons, 0 sp3 carbons => Fsp3 = 0.0.
        % Edge case: all carbons are non-sp3.
        % Complements ethanol Fsp3=1.0 test; together they bound the [0,1] range.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccccc1");
            desc = emk.descriptor.calculate(mol, "FractionCSP3");
            tc.verifyEqual(desc.FractionCSP3, 0.0, "AbsTol", 0.001, ...
                "Benzene FractionCSP3 must be 0.0 (all sp2 carbons)");
        end

        function test_calculate_nbutane_NumRotatableBonds(tc)
        % n-Butane (CCCC) has 1 rotatable bond under the strict definition.
        % C1-C2 and C3-C4 are excluded (terminal atoms, degree=1).
        % C2-C3 is the only bond where both atoms have degree > 1 => 1 bond.
        % Complements ethanol test (0 rotatable bonds); verifies non-zero output.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCCC");
            desc = emk.descriptor.calculate(mol, "NumRotatableBonds");
            tc.verifyEqual(desc.NumRotatableBonds, 1, ...
                "n-Butane must have 1 rotatable bond (strict definition)");
        end

        function test_calculate_benzene_RingCount(tc)
        % Benzene has 1 ring.
        % Pure aromatic system cross-validates the aspirin RingCount=1 result;
        % absence of substituents removes any ambiguity about ring counting.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccccc1");
            desc = emk.descriptor.calculate(mol, "RingCount");
            tc.verifyEqual(desc.RingCount, 1, ...
                "Benzene must have 1 ring");
        end

        % ------------------------------------------------------------------
        % TC16: Aspirin NumHAcceptors -- TC12 supplement (requires RDKit)
        % ------------------------------------------------------------------

        function test_calculate_aspirin_NumHAcceptors(tc)
        % Aspirin has 3 H-bond acceptors (Ertl definition):
        %   O1: acetyl C=O     -- not OH => HBA
        %   O2: ester -O-      -- not OH => HBA
        %   O3: carboxyl C=O   -- not OH => HBA
        %   O4: carboxyl OH    -- is OH  => HBD only, NOT HBA
        % Reference: PubChem CID 2244 (HBondAcceptorCount = 3).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            desc = emk.descriptor.calculate(mol, "NumHAcceptors");
            tc.verifyEqual(desc.NumHAcceptors, 3, ...
                "Aspirin must have 3 H-bond acceptors");
        end

        % ------------------------------------------------------------------
        % TC17: Consistency -- calculate vs molWeight (requires RDKit)
        % ------------------------------------------------------------------

        function test_calculate_MolWt_consistentWithMolWeight(tc)
        % calculate(mol, "MolWt").MolWt must equal molWeight(mol) exactly.
        % Both functions call Descriptors.MolWt via emk.util.rdkitModule()
        % (importlib path).  Any divergence indicates a regression in one
        % function.  Uses aspirin (complex molecule) to maximize sensitivity.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            desc = emk.descriptor.calculate(mol, "MolWt");
            mw   = emk.descriptor.molWeight(mol);
            tc.verifyEqual(desc.MolWt, mw, "AbsTol", 1e-9, ...
                "calculate MolWt must equal molWeight for the same molecule");
        end

    end

    % ======================================================================
    % TC18: Python-side cross-validation
    % ======================================================================
    methods (Test)

        function test_calculate_HeavyAtomCount_crossValidates_ethanol(tc)
        % Cross-validation: HeavyAtomCount must equal double(mol.GetNumAtoms()).
        % Bridges the MATLAB wrapper (emk.descriptor.calculate) and the
        % Python Mol object, confirming the conversion is lossless.
        % mol.GetNumAtoms() is called directly here (principle 3: MATLAB
        % wrapper result == Python API value).  Ethanol (CCO) = 3 heavy atoms.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "HeavyAtomCount");
            tc.verifyEqual(desc.HeavyAtomCount, double(mol.GetNumAtoms()), ...
                "HeavyAtomCount must equal mol.GetNumAtoms() for ethanol");
        end

        function test_calculate_HeavyAtomCount_crossValidates_aspirin(tc)
        % Cross-validation: HeavyAtomCount must equal double(mol.GetNumAtoms()).
        % Aspirin (13 heavy atoms) exercises a structurally complex molecule.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            desc = emk.descriptor.calculate(mol, "HeavyAtomCount");
            tc.verifyEqual(desc.HeavyAtomCount, double(mol.GetNumAtoms()), ...
                "HeavyAtomCount must equal mol.GetNumAtoms() for aspirin");
        end

        function test_calculate_countDescriptors_areIntegerValued_ethanol(tc)
        % NumHAcceptors, NumHDonors, NumRotatableBonds, RingCount, HeavyAtomCount
        % must all be integer-valued doubles (floor(v) == v).
        % Guards against incorrect Python float-to-MATLAB conversion for
        % count descriptors that are conceptually integers.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol);
            intFields = ["NumHAcceptors", "NumHDonors", "NumRotatableBonds", ...
                         "RingCount", "HeavyAtomCount"];
            for k = 1:numel(intFields)
                v = desc.(intFields(k));
                tc.verifyEqual(v, floor(v), ...
                    intFields(k) + " must be integer-valued (no fractional part)");
            end
        end

        function test_calculate_countDescriptors_areIntegerValued_aspirin(tc)
        % Same integer-value check for aspirin (more structurally complex).
        % RingCount=1, HeavyAtomCount=13, NumHDonors=1, NumHAcceptors=3.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            desc = emk.descriptor.calculate(mol);
            intFields = ["NumHAcceptors", "NumHDonors", "NumRotatableBonds", ...
                         "RingCount", "HeavyAtomCount"];
            for k = 1:numel(intFields)
                v = desc.(intFields(k));
                tc.verifyEqual(v, floor(v), ...
                    intFields(k) + " must be integer-valued for aspirin");
            end
        end

    end

    % ======================================================================
    % batchCalculate() tests
    %
    % Coverage:
    %   TB1: Non-cell input throws invalidInput (no RDKit required)
    %   TB2: Unknown descriptor name throws unknownDescriptor (no RDKit required)
    %   TB3: N molecules -> N-row table (requires RDKit)
    %   TB4: Output type is table (requires RDKit)
    %   TB5: Column names match requested descriptorNames (requires RDKit)
    %   TB6: Subset of descriptors -> correct column count (requires RDKit)
    %   TB7: Values match single calculate() call (requires RDKit)
    %   TB8: Invalid mol in cell -> NaN row, valid rows unaffected (requires RDKit)
    %   TB9: Empty cell array -> 0-row table with correct columns (no RDKit required)
    %        (descriptorNames must be specified to determine columns without mols)
    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TB1: Non-cell input throws invalidInput -- no RDKit required
        % ------------------------------------------------------------------

        function test_batchCalculate_numericInput_throwsInvalidInput(tc)
        % Numeric input must throw invalidInput before any RDKit call.
            tc.verifyError(@() emk.descriptor.batchCalculate(42), ...
                "emk:descriptor:batchCalculate:invalidInput", ...
                "Numeric input must throw invalidInput");
        end

        function test_batchCalculate_stringInput_throwsInvalidInput(tc)
        % String input must throw invalidInput.
            tc.verifyError(@() emk.descriptor.batchCalculate("CCO"), ...
                "emk:descriptor:batchCalculate:invalidInput", ...
                "String input must throw invalidInput");
        end

        % ------------------------------------------------------------------
        % TB2: Unknown descriptor name -- no RDKit required
        % ------------------------------------------------------------------

        function test_batchCalculate_unknownName_throwsUnknownDescriptor(tc)
        % Unknown descriptor name must throw unknownDescriptor.
            tc.verifyError(@() emk.descriptor.batchCalculate({}, "BadName"), ...
                "emk:descriptor:batchCalculate:unknownDescriptor", ...
                "Unknown descriptor name must throw unknownDescriptor");
        end

        function test_batchCalculate_unknownName_errorMessage_mentionsSupported(tc)
        % Error message must mention "Supported" to guide users.
            ME = tc.captureError( ...
                @() emk.descriptor.batchCalculate({}, "NotAName"));
            tc.assertNotEmpty(ME, "Expected an error to be thrown");
            tc.verifySubstring(ME.message, "Supported", ...
                "Error message must mention Supported descriptors");
        end

        % ------------------------------------------------------------------
        % TB9: Empty cell array with explicit descriptorNames -- no RDKit required
        % ------------------------------------------------------------------

        function test_batchCalculate_emptyMols_returnsZeroRowTable(tc)
        % Empty cell array with explicit names must return a 0-row table.
        % Column names must match the requested descriptorNames.
            tbl = emk.descriptor.batchCalculate({}, ["MolWt", "LogP"]);
            tc.verifyClass(tbl, "table", ...
                "Output must be a table even for empty mols");
            tc.verifySize(tbl, [0 2], ...
                "Empty mols must produce a 0-row, 2-column table");
            tc.verifyEqual(tbl.Properties.VariableNames, {'MolWt', 'LogP'}, ...
                "Column names must match requested descriptorNames");
        end

        % ------------------------------------------------------------------
        % TB3: N rows -- requires RDKit
        % ------------------------------------------------------------------

        function test_batchCalculate_threeSmiles_returnsThreeRows(tc)
        % Three valid molecules must produce a 3-row table.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            smilesList = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"};
            mols = cellfun(@(s) emk.mol.fromSmiles(s), smilesList, ...
                           "UniformOutput", false);
            tbl = emk.descriptor.batchCalculate(mols, ["MolWt", "LogP"]);
            tc.verifyEqual(height(tbl), 3, ...
                "Three molecules must produce 3 rows");
        end

        % ------------------------------------------------------------------
        % TB4: Output type is table -- requires RDKit
        % ------------------------------------------------------------------

        function test_batchCalculate_outputIsTable(tc)
        % batchCalculate must return a MATLAB table.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            tbl  = emk.descriptor.batchCalculate({mol});
            tc.verifyClass(tbl, "table", ...
                "batchCalculate must return a table");
        end

        % ------------------------------------------------------------------
        % TB5: Column names match descriptorNames -- requires RDKit
        % ------------------------------------------------------------------

        function test_batchCalculate_columnNames_matchDefault(tc)
        % Default call must produce 10 columns matching all supported descriptors.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol      = emk.mol.fromSmiles("CCO");
            tbl      = emk.descriptor.batchCalculate({mol});
            expected = {'MolWt','ExactMolWt','LogP','TPSA', ...
                        'NumHAcceptors','NumHDonors','NumRotatableBonds', ...
                        'RingCount','FractionCSP3','HeavyAtomCount'};
            tc.verifyEqual(tbl.Properties.VariableNames, expected, ...
                "Default columns must match all 10 supported descriptor names");
        end

        function test_batchCalculate_columnNames_matchSubset(tc)
        % Subset request must produce exactly those columns in the given order.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            tbl  = emk.descriptor.batchCalculate({mol}, ["TPSA", "MolWt"]);
            tc.verifyEqual(tbl.Properties.VariableNames, {'TPSA', 'MolWt'}, ...
                "Subset columns must match requested names in order");
        end

        % ------------------------------------------------------------------
        % TB6: Column count matches requested descriptors -- requires RDKit
        % ------------------------------------------------------------------

        function test_batchCalculate_subsetDescriptors_correctColumnCount(tc)
        % Requesting 3 descriptors must produce 3 columns.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            tbl  = emk.descriptor.batchCalculate({mol}, ["MolWt","LogP","TPSA"]);
            tc.verifySize(tbl, [1 3], ...
                "3-descriptor request must produce [1 3] table");
        end

        % ------------------------------------------------------------------
        % TB7: Values match single calculate() call -- requires RDKit
        % ------------------------------------------------------------------

        function test_batchCalculate_ethanolMolWt_matchesCalculate(tc)
        % batchCalculate MolWt for ethanol must equal calculate MolWt.
        % Cross-validates batch output against the single-mol reference.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            tbl  = emk.descriptor.batchCalculate({mol}, "MolWt");
            desc = emk.descriptor.calculate(mol, "MolWt");
            tc.verifyEqual(tbl.MolWt(1), desc.MolWt, ...
                "batchCalculate MolWt must equal calculate MolWt");
        end

        function test_batchCalculate_allValues_areFiniteDoubles(tc)
        % All values for valid molecules must be finite doubles.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            tbl  = emk.descriptor.batchCalculate({mol});
            vals = table2array(tbl(1, :));
            tc.verifyClass(vals, "double", ...
                "All table values must be double");
            tc.verifyTrue(all(isfinite(vals)), ...
                "All values for a valid molecule must be finite");
        end

        % ------------------------------------------------------------------
        % TB8: Invalid mol in cell produces NaN row -- requires RDKit
        % ------------------------------------------------------------------

        function test_batchCalculate_invalidMolInCell_nanRow(tc)
        % A non-Mol element in the cell array must yield NaN in its row.
        % Valid neighbors must be unaffected.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            ethanol = emk.mol.fromSmiles("CCO");
            aspirin = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            mols    = {ethanol, "invalid", aspirin};
            tbl     = emk.descriptor.batchCalculate(mols, "MolWt");
            % Row 2 (invalid) must be NaN
            tc.verifyTrue(isnan(tbl.MolWt(2)), ...
                "Invalid mol in cell must produce NaN in the row");
            % Rows 1 and 3 (valid) must be finite
            tc.verifyTrue(isfinite(tbl.MolWt(1)), ...
                "Valid mol before invalid mol must remain finite");
            tc.verifyTrue(isfinite(tbl.MolWt(3)), ...
                "Valid mol after invalid mol must remain finite");
        end

        % ------------------------------------------------------------------
        % TB10: All mols invalid => allMolsFailed (M3-1, no RDKit required)
        % ------------------------------------------------------------------

        function test_batchCalculate_allNonMolInputs_throwsAllMolsFailed(tc)
        % When every element in the cell array is a non-Mol object the
        % function must throw allMolsFailed.
        % No RDKit needed: all failures are detected by type check before
        % any Python descriptor computation is attempted.
            tc.verifyError( ...
                @() emk.descriptor.batchCalculate({"invalid", 42}), ...
                "emk:descriptor:batchCalculate:allMolsFailed", ...
                "All non-Mol inputs must throw allMolsFailed");
        end

        function test_batchCalculate_singleNonMol_throwsAllMolsFailed(tc)
        % Single-element cell with a non-Mol must also throw allMolsFailed.
        % Verifies that nMols=1 is handled correctly (not suppressed by a
        % minimum-element guard).
            tc.verifyError( ...
                @() emk.descriptor.batchCalculate({"not_a_mol"}), ...
                "emk:descriptor:batchCalculate:allMolsFailed", ...
                "Single non-Mol input must throw allMolsFailed");
        end

        function test_batchCalculate_allMolsFailed_errorMessage_containsCount(tc)
        % Error message must contain the total molecule count for diagnosis.
        % Uses a 2-element cell so the count '2' is unambiguous.
            ME = tc.captureError( ...
                @() emk.descriptor.batchCalculate({"bad1", "bad2"}));
            tc.assertNotEmpty(ME, "Expected allMolsFailed to be thrown");
            tc.verifySubstring(ME.message, "2", ...
                "Error message must contain molecule count '2'");
        end

    end

    % ======================================================================
    % calculate() -- TC19: MolFormula descriptor (string-valued)
    % ======================================================================
    methods (Test)

        % ------------------------------------------------------------------
        % TC19a: MolFormula validation -- no RDKit required
        % ------------------------------------------------------------------

        function test_calculate_MolFormula_notUnknown_passesValidation(tc)
        % "MolFormula" must NOT trigger unknownDescriptor.
        % It is a supported but non-default descriptor.  Validation must
        % accept it before any mol-type check occurs.
            ME = tc.captureError( ...
                @() emk.descriptor.calculate(struct(), "MolFormula"));
            if ~isempty(ME)
                tc.verifyNotEqual(ME.identifier, ...
                    "emk:descriptor:calculate:unknownDescriptor", ...
                    "'MolFormula' must not trigger unknownDescriptor");
            end
        end

        function test_calculate_MolFormula_notInDefaultSet(tc)
        % Default calculate() (no names arg) must NOT include MolFormula.
        % MolFormula returns string while all default descriptors return
        % double; mixing types would break code that assumes double fields.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol);
            tc.verifyFalse(isfield(desc, "MolFormula"), ...
                "Default calculate() must not include MolFormula field");
        end

        % ------------------------------------------------------------------
        % TC19b: MolFormula output type and shape (requires RDKit)
        % ------------------------------------------------------------------

        function test_calculate_MolFormula_outputIsStringScalar(tc)
        % MolFormula must return a string scalar (not char, not cell).
        % Verifies that the Python str is converted via string() not char().
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "MolFormula");
            tc.verifyClass(desc.MolFormula, "string", ...
                "MolFormula must return a string (not char or cell)");
            tc.verifySize(desc.MolFormula, [1 1], ...
                "MolFormula must return a scalar [1 1]");
        end

        function test_calculate_MolFormula_outputIsNonEmpty(tc)
        % MolFormula must be a non-empty string for a valid molecule.
        % Guards against empty string returned when CalcMolFormula fails silently.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "MolFormula");
            tc.verifyNotEqual(desc.MolFormula, "", ...
                "MolFormula must be a non-empty string for a valid molecule");
        end

        % ------------------------------------------------------------------
        % TC19c: MolFormula reference values (requires RDKit)
        % Reference: Hill notation (Hill, E.A. (1900) J. Am. Chem. Soc. 22(8):478-494)
        %            C first, H second, then remaining elements alphabetically.
        % ------------------------------------------------------------------

        function test_calculate_MolFormula_ethanol(tc)
        % Ethanol (CCO) molecular formula = "C2H6O".
        % Reference: PubChem CID 702 (MolecularFormula = C2H6O).
        % Includes 6 implicit H (2 on C1, 3 on C2, 1 on O).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, "MolFormula");
            tc.verifyEqual(desc.MolFormula, "C2H6O", ...
                "Ethanol molecular formula must be C2H6O");
        end

        function test_calculate_MolFormula_aspirin(tc)
        % Aspirin (CC(=O)Oc1ccccc1C(=O)O) molecular formula = "C9H8O4".
        % Reference: PubChem CID 2244 (MolecularFormula = C9H8O4).
        % Complex molecule: aromatic ring + ester + carboxylic acid.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
            desc = emk.descriptor.calculate(mol, "MolFormula");
            tc.verifyEqual(desc.MolFormula, "C9H8O4", ...
                "Aspirin molecular formula must be C9H8O4");
        end

        function test_calculate_MolFormula_benzene(tc)
        % Benzene (c1ccccc1) molecular formula = "C6H6".
        % Reference: PubChem CID 241. Edge case: aromatic implicit H.
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("c1ccccc1");
            desc = emk.descriptor.calculate(mol, "MolFormula");
            tc.verifyEqual(desc.MolFormula, "C6H6", ...
                "Benzene molecular formula must be C6H6");
        end

        % ------------------------------------------------------------------
        % TC19d: MolFormula can be combined with numeric descriptors (requires RDKit)
        % ------------------------------------------------------------------

        function test_calculate_MolFormula_withNumericDescriptor(tc)
        % Requesting ["ExactMolWt","MolFormula"] must return a struct with
        % both fields: ExactMolWt (double) and MolFormula (string).
        % This is the primary a08 use-case (exact mass + formula together).
            tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
            mol  = emk.mol.fromSmiles("CCO");
            desc = emk.descriptor.calculate(mol, ["ExactMolWt","MolFormula"]);
            tc.verifyEqual(numel(fieldnames(desc)), 2, ...
                "Struct must have exactly 2 fields");
            tc.verifyTrue(isfield(desc, "ExactMolWt"), ...
                "ExactMolWt field must be present");
            tc.verifyTrue(isfield(desc, "MolFormula"), ...
                "MolFormula field must be present");
            tc.verifyClass(desc.ExactMolWt, "double", ...
                "ExactMolWt must be double");
            tc.verifyClass(desc.MolFormula, "string", ...
                "MolFormula must be string");
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
