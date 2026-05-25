function smilesArray = randomSmiles(mol, n)
% randomSmiles  Generate random (non-canonical) SMILES enumerations of a molecule.
%
%   smilesArray = emk.mol.randomSmiles(mol, n)
%
%   Calls RDKit MolToSmiles with doRandom=True to generate n non-canonical
%   SMILES representations of the same molecule.  Returned strings are
%   deduplicated; simple molecules may return fewer than n unique strings.
%
%   This function supports SMILES enumeration data augmentation (Bjerrum 2017):
%   training a language model on multiple valid SMILES per molecule effectively
%   multiplies the corpus size and improves generalisation at small dataset sizes.
%
%   Arguments:
%     mol  - py.rdkit.Chem.rdchem.Mol  RDKit molecule object (Python ref)
%     n    - positive integer  number of random SMILES to attempt
%
%   Returns:
%     smilesArray - string column vector  (up to n unique random SMILES;
%                   may be fewer for small / symmetric molecules)
%
%   Errors:
%     emk:mol:randomSmiles:invalidInput  - bad argument type or value
%     emk:mol:randomSmiles:rdkitError    - unexpected Python exception
%
%   Example:
%     mol    = emk.mol.fromSmiles("c1ccccc1CC(=O)O");
%     rands  = emk.mol.randomSmiles(mol, 8);
%     disp(rands)     % up to 8 distinct traversal-order SMILES
%
%   Reference:
%     Bjerrum EJ (2017) SMILES enumeration as data augmentation for neural
%     network modeling of molecules. arXiv:1703.07076.
%
%   See also: emk.mol.toSmiles, emk.mol.fromSmiles

    % --- Input validation ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:mol:randomSmiles:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end
    if ~isnumeric(n) || ~isscalar(n) || n < 1 || n ~= round(n)
        error("emk:mol:randomSmiles:invalidInput", ...
            "n must be a positive integer, got: %s", mat2str(n));
    end
    n = double(n);

    logDebug("randomSmiles: generating %d random SMILES", n);

    mods = emk.util.rdkitModule();

    collected = strings(n, 1);
    for k = 1:n
        try
            s = mods.Chem.MolToSmiles(mol, ...
                pyargs("canonical", false, "doRandom", true));
            collected(k) = string(s);
        catch ME
            error("emk:mol:randomSmiles:rdkitError", ...
                "RDKit MolToSmiles(doRandom=True) raised: %s", ME.message);
        end
    end

    smilesArray = unique(collected);   % deduplicate (symmetric/simple mols repeat)
end
