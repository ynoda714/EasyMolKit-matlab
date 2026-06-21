function counts = fragmentCount(mol)
% fragmentCount  Ring system and functional group fragment counts.
%
%   counts = emk.descriptor.fragmentCount(mol)
%
%   Computes structural fragment counts for a molecule: ring system counts
%   (total, aromatic, saturated, heterocyclic) and common functional group
%   counts (carbonyl, amine, hydroxyl, halogen, nitrile, sulfonamide, amide).
%
%   Arguments:
%     mol - py.rdkit.Chem.rdchem.Mol  RDKit molecule object
%
%   Returns:
%     counts - struct with fields:
%       NumRings         - int  Total number of SSSR rings
%       NumAromaticRings - int  Aromatic rings
%       NumAliphaticRings - int  Aliphatic (saturated/partially saturated) rings
%       NumHeteroRings   - int  Rings containing at least one non-C atom
%       NumCarbonyl      - int  C=O groups (ketone, aldehyde, acid, ester, amide)
%       NumAmine         - int  N atoms with at least 1 H (primary/secondary amine)
%       NumHydroxyl      - int  OH groups (alcohol, phenol; not acid OH)
%       NumHalogen       - int  F, Cl, Br, I atoms
%       NumNitrile       - int  C#N groups
%       NumSulfonamide   - int  S(=O)(=O)N groups
%       NumAmide         - int  C(=O)N groups
%
%   Errors:
%     emk:descriptor:fragmentCount:invalidInput - mol is not a Mol object
%     emk:descriptor:fragmentCount:rdkitError   - unexpected Python exception
%
%   Example:
%     mol    = emk.mol.fromSmiles("CC(=O)Nc1ccc(O)cc1");  % paracetamol
%     counts = emk.descriptor.fragmentCount(mol);
%     fprintf("Aromatic rings: %d, Amide groups: %d, Hydroxyl: %d\n", ...
%         counts.NumAromaticRings, counts.NumAmide, counts.NumHydroxyl);
%
%   References:
%     RDKit Documentation: rdkit.Chem.rdMolDescriptors
%       CalcNumRings, CalcNumAromaticRings, CalcNumAliphaticRings,
%       CalcNumHeterocycles, CalcNumHeteroatoms
%     RDKit Fragments: rdkit.Chem.Fragments (fr_* SMARTS patterns)
%
%   See also: emk.descriptor.qed, emk.descriptor.bcut,
%             emk.mol.computeDescriptors

    % --- Validate mol ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:descriptor:fragmentCount:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    try
        descMod = py.importlib.import_module("rdkit.Chem.rdMolDescriptors");
        fragMod = py.importlib.import_module("rdkit.Chem.Fragments");

        % Ring counts
        counts.NumRings          = double(descMod.CalcNumRings(mol));
        counts.NumAromaticRings  = double(descMod.CalcNumAromaticRings(mol));
        counts.NumAliphaticRings = double(descMod.CalcNumAliphaticRings(mol));
        counts.NumHeteroRings    = double(descMod.CalcNumHeterocycles(mol));

        % Functional group counts via RDKit Fragments module (SMARTS-based)
        % fr_C_O counts all C=O groups (ketone, aldehyde, acid, ester, amide)
        counts.NumCarbonyl    = double(fragMod.fr_C_O(mol));
        counts.NumAmine       = double(fragMod.fr_NH0(mol)) + ...
                                double(fragMod.fr_NH1(mol)) + ...
                                double(fragMod.fr_NH2(mol));
        counts.NumHydroxyl    = double(fragMod.fr_Al_OH(mol)) + ...
                                double(fragMod.fr_Ar_OH(mol));
        counts.NumHalogen     = double(fragMod.fr_halogen(mol));
        counts.NumNitrile     = double(fragMod.fr_nitrile(mol));
        counts.NumSulfonamide = double(fragMod.fr_sulfonamd(mol));
        counts.NumAmide       = double(fragMod.fr_amide(mol));

    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:descriptor:fragmentCount:rdkitError", ...
            "Fragment count calculation failed: %s", ME.message);
    end

    logDebug("descriptor.fragmentCount: %d rings, %d aromatic", ...
        counts.NumRings, counts.NumAromaticRings);
end
