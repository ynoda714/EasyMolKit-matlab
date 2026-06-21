function score = qed(mol)
% qed  Quantitative Estimate of Drug-likeness (QED) score.
%
%   score = emk.descriptor.qed(mol)
%
%   Computes the QED score, a composite measure of drug-likeness combining
%   eight molecular properties into a single score using a desirability
%   function approach.  Higher values indicate more drug-like properties.
%
%   The eight properties combined by QED:
%     MW    - Molecular weight
%     ALOGP - Calculated LogP
%     HBA   - Hydrogen bond acceptors
%     HBD   - Hydrogen bond donors
%     PSA   - Polar surface area
%     ROTB  - Rotatable bonds
%     AROM  - Aromatic rings
%     ALERTS - Number of structural alerts
%
%   Arguments:
%     mol - py.rdkit.Chem.rdchem.Mol  RDKit molecule object
%
%   Returns:
%     score - double in [0, 1]  QED score.
%             Higher is more drug-like (0.67+ is a common threshold).
%             Typical approved drugs have QED > 0.67.
%
%   Errors:
%     emk:descriptor:qed:invalidInput - mol is not a Mol object
%     emk:descriptor:qed:rdkitError   - unexpected Python exception
%
%   Example:
%     mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
%     score = emk.descriptor.qed(mol);
%     fprintf("QED: %.3f\n", score);  % ~0.55 for aspirin
%
%   References:
%     Bickerton, G.R. et al. (2012). Quantifying the Chemical Beauty of
%       Drugs. Nature Chemistry 4(2):90-98. DOI: 10.1038/nchem.1243
%     RDKit Documentation: rdkit.Chem.QED.qed
%
%   See also: emk.descriptor.saScore, emk.filter.lipinski, emk.filter.veber

    % --- Validate mol ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:descriptor:qed:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    try
        qedMod = py.importlib.import_module("rdkit.Chem.QED");
        score  = double(qedMod.qed(mol));
    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:descriptor:qed:rdkitError", ...
            "QED calculation failed: %s", ME.message);
    end

    logDebug("descriptor.qed: score=%.4f", score);
end
