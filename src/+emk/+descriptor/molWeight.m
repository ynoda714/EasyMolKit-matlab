function mw = molWeight(mol)
% molWeight  Compute the average molecular weight of a molecule.
%
%   mw = emk.descriptor.molWeight(mol)
%
%   Returns the average molecular weight (g/mol) using RDKit's
%   Descriptors.MolWt.  This sums the average atomic masses of all atoms
%   including implicit hydrogens, using IUPAC standard atomic weights.
%   For monoisotopic mass use emk.descriptor.exactMolWeight (not yet
%   implemented).
%
%   Arguments:
%     mol  - py.rdkit.Chem.rdchem.Mol  RDKit molecule object (Python ref)
%
%   Returns:
%     mw   - double  Average molecular weight in g/mol
%
%   Errors:
%     emk:descriptor:molWeight:invalidInput  - mol is not Mol object
%     emk:descriptor:molWeight:rdkitError    - unexpected Python exception
%
%   Example:
%     mol = emk.mol.fromSmiles("CCO");
%     mw  = emk.descriptor.molWeight(mol);   % ~46.07 g/mol (ethanol)
%
%   See also: emk.mol.fromSmiles, emk.descriptor.calculate

    % --- Input validation ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:descriptor:molWeight:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    logDebug("molWeight: computing MW for mol with %d heavy atoms", ...
        double(mol.GetNumAtoms()));

    % --- Call RDKit ---
    % Use emk.util.rdkitModule() to get the Descriptors module via importlib.
    % Direct py.rdkit.* access fails (TypeError); importlib avoids this.
    mods = emk.util.rdkitModule();
    try
        pyMW = mods.Descriptors.MolWt(mol);
    catch ME
        error("emk:descriptor:molWeight:rdkitError", ...
            "RDKit Descriptors.MolWt raised an exception: %s", ME.message);
    end

    mw = double(pyMW);

    logDebug("molWeight: MW = %.4f g/mol", mw);
end
