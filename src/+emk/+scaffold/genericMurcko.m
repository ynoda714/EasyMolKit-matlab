function scaf = genericMurcko(mol)
% genericMurcko  Compute the generic Murcko scaffold (all-carbon skeleton).
%
%   scaf = emk.scaffold.genericMurcko(mol)
%
%   Extracts the Bemis-Murcko scaffold and then strips all heteroatom
%   identity by replacing every atom with carbon and every bond with a
%   single bond.  The result is the carbon-skeleton framework used for
%   scaffold diversity comparisons that treat ring-type and heteroatom
%   substitution as equivalent.
%
%   Processing pipeline (two RDKit calls):
%     1. MurckoScaffold.GetScaffoldForMol(mol)  -> Bemis-Murcko scaffold
%     2. MurckoScaffold.MakeScaffoldGeneric(scaf) -> all atoms C, all bonds single
%
%   For acyclic molecules (no ring systems) the Murcko scaffold has 0 atoms;
%   a 0-atom Mol is returned -- consistent with emk.mol.scaffold.
%
%   Arguments:
%     mol  - py.rdkit.Chem.rdchem.Mol  RDKit molecule object (Python ref)
%
%   Returns:
%     scaf - py.rdkit.Chem.rdchem.Mol  Generic Murcko scaffold Mol object.
%            All ring atoms are carbon; all bonds are single bonds.
%            May have 0 atoms for acyclic molecules.
%
%   Errors:
%     emk:scaffold:genericMurcko:invalidInput - mol is not a Mol object
%     emk:scaffold:genericMurcko:rdkitError   - unexpected Python exception
%
%   Example:
%     mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
%     scaf = emk.scaffold.genericMurcko(mol);
%     emk.mol.toSmiles(scaf)   % "C1CCCCC1" (cyclohexane carbon skeleton)
%
%     mol2 = emk.mol.fromSmiles("c1cc2ccccc2nc1");  % quinoline (N-heterocycle)
%     scaf2 = emk.scaffold.genericMurcko(mol2);
%     emk.mol.toSmiles(scaf2)  % "C1CCC2CCCCC2C1" (decalin-like skeleton)
%
%   References:
%     Bemis, G.W. & Murcko, M.A. (1996). The Properties of Known Drugs.
%       1. Molecular Frameworks. J. Med. Chem. 39(15):2887-2893.
%       DOI: 10.1021/jm9602928
%     Bemis, G.W. & Murcko, M.A. (1999). Properties of Known Drugs.
%       2. Side Chains. J. Med. Chem. 42(25):5095-5099.
%       DOI: 10.1021/jm9903996
%
%   See also: emk.mol.scaffold, emk.scaffold.brics, emk.scaffold.rgroup,
%             emk.mol.toSmiles

    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:scaffold:genericMurcko:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    logDebug("genericMurcko: computing generic Murcko scaffold");

    try
        murckoMod = py.importlib.import_module( ...
            "rdkit.Chem.Scaffolds.MurckoScaffold");

        % Step 1: Bemis-Murcko scaffold (retains heteroatoms and bond orders)
        murcko = murckoMod.GetScaffoldForMol(mol);
        if isa(murcko, "py.NoneType")
            error("emk:scaffold:genericMurcko:rdkitError", ...
                "GetScaffoldForMol returned None");
        end

        % Step 2: Generic scaffold -- all atoms -> C, all bonds -> single
        scaf = murckoMod.MakeScaffoldGeneric(murcko);
        if isa(scaf, "py.NoneType")
            error("emk:scaffold:genericMurcko:rdkitError", ...
                "MakeScaffoldGeneric returned None");
        end
    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:scaffold:genericMurcko:rdkitError", ...
            "MurckoScaffold operation failed: %s", ME.message);
    end

    logDebug("genericMurcko: done, scaffold has %d atoms", ...
        double(scaf.GetNumAtoms()));
end
