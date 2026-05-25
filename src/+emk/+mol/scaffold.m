function scaf = scaffold(mol)
% scaffold  Extract the Murcko scaffold from an RDKit Mol object.
%
%   scaf = emk.mol.scaffold(mol)
%
%   Extracts the Bemis-Murcko scaffold from a molecule using RDKit's
%   MurckoScaffold module.  The scaffold is defined as the union of all
%   ring systems and the linker bonds between them, with all side chains
%   stripped.  For acyclic molecules (no ring systems), an empty Mol object
%   (0 atoms) is returned -- this is the canonical RDKit behaviour.
%
%   The returned object is a py.rdkit.Chem.rdchem.Mol Python reference
%   (ADR-002: Python object reference retention).  Use emk.mol.toSmiles to
%   convert to a SMILES string, or pass directly to emk.fingerprint.*,
%   emk.descriptor.*, etc.
%
%   Arguments:
%     mol  - py.rdkit.Chem.rdchem.Mol  RDKit molecule object (Python ref)
%
%   Returns:
%     scaf - py.rdkit.Chem.rdchem.Mol  Murcko scaffold molecule object.
%            May have 0 atoms if mol contains no ring systems.
%
%   Errors:
%     emk:mol:scaffold:invalidInput  - mol is not a py.rdkit.Chem.rdchem.Mol
%     emk:mol:scaffold:rdkitError    - unexpected Python exception
%
%   Example:
%     mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
%     scaf = emk.mol.scaffold(mol);
%     emk.mol.toSmiles(scaf)                                 % "c1ccccc1"
%
%     % Acyclic molecule => empty scaffold (0 atoms)
%     mol2  = emk.mol.fromSmiles("CCO");  % ethanol
%     scaf2 = emk.mol.scaffold(mol2);
%     double(scaf2.GetNumAtoms()) == 0    % true
%
%   Note:
%     RDKit's GetScaffoldForMol preserves ring heteroatoms (Bemis-Murcko
%     framework).  For all-carbon skeletons use GetScaffoldForMol with the
%     makeScaffoldGeneric option inside a custom py.* call.
%
%   References:
%     Bemis, G.W. & Murcko, M.A. (1996). The Properties of Known Drugs.
%       1. Molecular Frameworks. J. Med. Chem. 39(15):2887-2893.
%       DOI: 10.1021/jm9602928
%     RDKit Documentation: MurckoScaffold
%       https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
%
%   See also: emk.mol.fromSmiles, emk.mol.toSmiles, emk.mol.hasSubstruct

    % --- Input validation ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:mol:scaffold:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    logDebug("scaffold: calling MurckoScaffold.GetScaffoldForMol");

    % --- Import MurckoScaffold module and extract scaffold ---
    % py.rdkit.Chem.Scaffolds.MurckoScaffold is not auto-loaded by MATLAB's
    % py.* namespace; py.importlib.import_module ensures reliable loading
    % (same pattern as emk.viz.draw2d for AllChem and Draw modules).
    try
        murckoMod = py.importlib.import_module( ...
            "rdkit.Chem.Scaffolds.MurckoScaffold");
        scaf = murckoMod.GetScaffoldForMol(mol);
    catch ME
        error("emk:mol:scaffold:rdkitError", ...
            "MurckoScaffold.GetScaffoldForMol failed: %s", ME.message);
    end

    % --- None check ---
    % GetScaffoldForMol returns an empty Mol (not None) for acyclic molecules,
    % but guard defensively against None in case of unexpected RDKit behaviour.
    if isa(scaf, "py.NoneType")
        error("emk:mol:scaffold:rdkitError", ...
            "MurckoScaffold.GetScaffoldForMol returned None for the input molecule");
    end

    nAtoms = double(scaf.GetNumAtoms());
    logDebug("scaffold: done, scaffold has %d atoms", nAtoms);
end
