function score = saScore(mol)
% saScore  Synthetic Accessibility (SA) score.
%
%   score = emk.descriptor.saScore(mol)
%
%   Computes the SA Score, a measure of how easy it is to synthesise a
%   molecule.  The score combines fragment contributions (learned from
%   PubChem) with a complexity penalty based on ring systems, stereocentres,
%   and macrocycle features.
%
%   Score interpretation:
%     1.0 - 3.0  Easy to synthesise
%     3.0 - 6.0  Moderate difficulty
%     6.0 - 10.0 Hard/impossible to synthesise
%
%   NOTE: SA Score requires the RDKit Contrib/SA_Score package to be on
%   the Python path.  In standard RDKit distributions this is available as
%   rdkit.Contrib.SA_Score.sascorer (RDKit >= 2022.03) or via the
%   rdkit-stdinchi / rdkit-utilities path.  A Python RuntimeError is raised
%   if the contrib module is not available.
%
%   Arguments:
%     mol - py.rdkit.Chem.rdchem.Mol  RDKit molecule object
%
%   Returns:
%     score - double in [1, 10]  SA score (lower = easier to synthesise).
%
%   Errors:
%     emk:descriptor:saScore:invalidInput - mol is not a Mol object
%     emk:descriptor:saScore:rdkitError   - SA_Score module not available
%                                           or unexpected Python exception
%
%   Example:
%     mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
%     score = emk.descriptor.saScore(mol);
%     fprintf("SA Score: %.2f\n", score);  % ~2.3 (easy to synthesise)
%
%   References:
%     Ertl, P. & Schuffenhauer, A. (2009). Estimation of Synthetic
%       Accessibility Score of Drug-like Molecules Based on Molecular
%       Complexity and Fragment Contributions. J. Cheminform. 1(8):1-11.
%       DOI: 10.1186/1758-2946-1-8
%     RDKit Contrib: rdkit.Contrib.SA_Score.sascorer
%
%   See also: emk.descriptor.qed, emk.mol.fromSmiles

    % --- Validate mol ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:descriptor:saScore:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    try
        sascorer = py.importlib.import_module("rdkit.Contrib.SA_Score.sascorer");
        score    = double(sascorer.calculateScore(mol));
    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:descriptor:saScore:rdkitError", ...
            "SA Score calculation failed: %s. " + ...
            "Ensure rdkit.Contrib.SA_Score is available in your Python environment.", ...
            ME.message);
    end

    logDebug("descriptor.saScore: score=%.4f", score);
end
