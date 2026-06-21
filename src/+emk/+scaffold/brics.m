function frags = brics(mol)
% brics  Fragment a molecule using the BRICS retrosynthetic rules.
%
%   frags = emk.scaffold.brics(mol)
%
%   Applies BRICS (Breaking Retrosynthetically Interesting Chemical
%   Substructures) fragmentation to identify synthetically accessible
%   building blocks.  The algorithm identifies retrosynthetically
%   interesting bonds as defined by Degen et al. (2008) and breaks them,
%   labelling each attachment point with a numbered dummy atom: [1*],
%   [3*], ... [16*] (BRICS environment types 1-16).
%
%   If the molecule cannot be fragmented (single ring, simple acyclic,
%   or too small), BRICSDecompose returns the original molecule as a
%   1-element set.  This function always returns at least one SMILES.
%
%   Note: The order of returned fragments is non-deterministic because
%   BRICSDecompose uses a Python frozenset internally.  Sort the output
%   if reproducible ordering is required.
%
%   Arguments:
%     mol   - py.rdkit.Chem.rdchem.Mol  RDKit molecule object
%
%   Returns:
%     frags - string(1,N)  SMILES strings of the N BRICS fragments.
%             Attachment point dummy atoms appear as [1*], [3*] etc.
%             N >= 1 (original mol SMILES if no fragmentable bond found).
%
%   Errors:
%     emk:scaffold:brics:invalidInput - mol is not a py.rdkit.Chem.rdchem.Mol
%     emk:scaffold:brics:rdkitError   - unexpected Python exception
%
%   Example:
%     mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
%     frags = emk.scaffold.brics(mol);
%     % frags contains fragment SMILES such as:
%     %   "[1*]OC(C)=O"   (acetyl fragment)
%     %   "[3*]c1ccccc1[3*]"  (aromatic linker fragment)
%     %   etc. (exact SMILES depend on BRICS bond detection)
%
%     % Acyclic molecule: typically returns original molecule
%     mol2  = emk.mol.fromSmiles("CCO");
%     frags2 = emk.scaffold.brics(mol2);
%     % frags2 == "CCO" (no BRICS bonds to cut)
%
%   References:
%     Degen, J., Wegscheid-Gerlach, C., Zaliani, A. & Rarey, M. (2008).
%       On the Art of Compiling and Using Drug-Like Chemical Fragment
%       Spaces. ChemMedChem 3(10):1503-1507. DOI: 10.1002/cmdc.200800178
%     RDKit Documentation: rdkit.Chem.BRICS.BRICSDecompose
%
%   See also: emk.scaffold.genericMurcko, emk.scaffold.rgroup,
%             emk.mol.scaffold, emk.mol.fromSmiles

    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:scaffold:brics:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    logDebug("brics: applying BRICS fragmentation");

    try
        bricsmod = py.importlib.import_module("rdkit.Chem.BRICS");

        % BRICSDecompose returns a Python frozenset of SMILES strings.
        % py.list(frozenset) converts it to a Python list for indexed access.
        pyfragset = bricsmod.BRICSDecompose(mol);
        pylist    = py.list(pyfragset);
        n = double(py.len(pylist));

        frags = strings(1, n);
        for i = 1:n
            frags(i) = string(pylist{i});
        end
    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:scaffold:brics:rdkitError", ...
            "BRICS.BRICSDecompose failed: %s", ME.message);
    end

    logDebug("brics: produced %d fragment(s)", numel(frags));
end
