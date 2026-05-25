function desc = calculate(mol, descriptorNames)
% calculate  Compute physicochemical descriptors for a molecule.
%
%   desc = emk.descriptor.calculate(mol)
%   desc = emk.descriptor.calculate(mol, ["MolWt","LogP","TPSA"])
%
%   Computes one or more RDKit-based physicochemical descriptors and
%   returns them as a MATLAB struct.  When descriptorNames is omitted all
%   ten supported descriptors are computed.
%
%   Arguments:
%     mol             - py.rdkit.Chem.rdchem.Mol  RDKit molecule object
%     descriptorNames - string array (optional)  Names of descriptors to
%                       compute.  Default: all supported descriptors.
%
%   Returns:
%     desc - struct  Field names = descriptor names, values = double scalar.
%
%   Supported descriptor names:
%     MolWt              Average molecular weight (g/mol)            [double]
%     ExactMolWt         Monoisotopic molecular weight (g/mol)        [double]
%     LogP               Wildman-Crippen LogP (unitless)               [double]
%     TPSA               Topological Polar Surface Area (Ang^2)        [double]
%     NumHAcceptors      H-bond acceptors (Ertl definition)            [double]
%     NumHDonors         H-bond donors (Ertl definition)               [double]
%     NumRotatableBonds  Rotatable bonds (strict SMARTS definition)    [double]
%     RingCount          Total ring count (SSSR)                       [double]
%     FractionCSP3       Fraction of sp3 carbons (Lovering Fsp3)      [double]
%     HeavyAtomCount     Number of heavy (non-H) atoms                 [double]
%     MolFormula         Molecular formula string (e.g., "C2H6O")     [string]
%
%   Errors:
%     emk:descriptor:calculate:invalidInput       - mol is not Mol object
%     emk:descriptor:calculate:unknownDescriptor  - unrecognized name(s)
%     emk:descriptor:calculate:rdkitError         - unexpected Python exception
%
%   Example:
%     mol  = emk.mol.fromSmiles("CCO");
%     desc = emk.descriptor.calculate(mol);
%     desc = emk.descriptor.calculate(mol, ["MolWt","LogP","TPSA"]);
%     mw   = desc.MolWt;   % 46.069 g/mol (ethanol)
%     fml  = emk.descriptor.calculate(mol, "MolFormula").MolFormula;  % "C2H6O"
%
%   Note: MolFormula returns a string scalar, not double.  It is excluded
%   from the default set (no-arg call) to preserve backward compatibility.
%
%   See also: emk.descriptor.molWeight, emk.mol.fromSmiles,
%             emk.descriptor.batchCalculate

    SUPPORTED_NAMES = supportedNames_();

    % --- Input validation: descriptorNames (checked before mol so that
    %     unknown-name errors fire regardless of mol validity) ---
    if nargin < 2 || isempty(descriptorNames)
        descriptorNames = defaultNames_();  % 10 numeric descriptors (no MolFormula)
    else
        if ischar(descriptorNames)
            descriptorNames = string(descriptorNames);
        end
        descriptorNames = reshape(string(descriptorNames), 1, []);
        unknown = descriptorNames(~ismember(descriptorNames, SUPPORTED_NAMES));
        if ~isempty(unknown)
            error("emk:descriptor:calculate:unknownDescriptor", ...
                "Unknown descriptor(s): [%s]. Supported: [%s]", ...
                strjoin(unknown, ", "), strjoin(SUPPORTED_NAMES, ", "));
        end
    end

    % --- Input validation: mol type ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:descriptor:calculate:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    logDebug("calculate: computing %d descriptor(s) for mol (%d heavy atoms)", ...
        numel(descriptorNames), double(mol.GetNumHeavyAtoms()));

    % --- Compute each descriptor ---
    desc = struct();
    for i = 1:numel(descriptorNames)
        name = descriptorNames(i);
        try
            desc.(name) = computeOne_(mol, name);
        catch ME
            error("emk:descriptor:calculate:rdkitError", ...
                "RDKit raised an exception computing '%s': %s", name, ME.message);
        end
    end

    logDebug("calculate: done (%d descriptor(s) computed)", numel(descriptorNames));
end

% -------------------------------------------------------------------------
function names = supportedNames_()
% Return the ordered list of all supported descriptor names (numeric + string).
    names = [defaultNames_(), "MolFormula"];
end

% -------------------------------------------------------------------------
function names = defaultNames_()
% Return the 10 numeric descriptors used when no names argument is given.
% MolFormula is intentionally excluded because it returns string not double.
    names = ["MolWt", "ExactMolWt", "LogP", "TPSA", ...
             "NumHAcceptors", "NumHDonors", "NumRotatableBonds", ...
             "RingCount", "FractionCSP3", "HeavyAtomCount"];
end

% -------------------------------------------------------------------------
function val = computeOne_(mol, name)
% Dispatch a single descriptor computation to RDKit.
% All returned Python scalars are converted to MATLAB double.
    mods = emk.util.rdkitModule();
    switch name
        case "MolWt"
            val = double(mods.Descriptors.MolWt(mol));
        case "ExactMolWt"
            val = double(mods.Descriptors.ExactMolWt(mol));
        case "LogP"
            val = double(mods.Descriptors.MolLogP(mol));
        case "TPSA"
            val = double(mods.Descriptors.TPSA(mol));
        case "NumHAcceptors"
            val = double(mods.rdMolDescriptors.CalcNumHBA(mol));
        case "NumHDonors"
            val = double(mods.rdMolDescriptors.CalcNumHBD(mol));
        case "NumRotatableBonds"
            val = double(mods.rdMolDescriptors.CalcNumRotatableBonds(mol));
        case "RingCount"
            val = double(mods.rdMolDescriptors.CalcNumRings(mol));
        case "FractionCSP3"
            val = double(mods.Descriptors.FractionCSP3(mol));
        case "HeavyAtomCount"
            val = double(mol.GetNumHeavyAtoms());
        case "MolFormula"
            % CalcMolFormula returns a Python str; convert to MATLAB string.
            % Includes implicit H and all elements in Hill order (C, H first).
            % Reference: RDKit rdMolDescriptors.CalcMolFormula (no literature
            % formula definition -- Hill notation is the de facto standard:
            % Hill, E.A. (1900) J. Am. Chem. Soc. 22(8):478-494).
            val = string(mods.rdMolDescriptors.CalcMolFormula(mol));
        otherwise
            % Defensive: unreachable after name validation in calculate().
            error("emk:descriptor:calculate:unknownDescriptor", ...
                "Unrecognized descriptor: %s", name);
    end
end
