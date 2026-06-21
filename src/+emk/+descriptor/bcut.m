function scores = bcut(mol)
% bcut  BCUT2D molecular descriptors.
%
%   scores = emk.descriptor.bcut(mol)
%
%   Computes BCUT2D descriptors, which encode molecular topology and atomic
%   properties via eigenvalues of the Burden matrix (a modified adjacency
%   matrix weighted by atomic properties).  BCUT descriptors capture both
%   connectivity and property information in a compact form.
%
%   The 8 BCUT2D descriptors returned (in order):
%     BCUT2D_MWHI   - Max eigenvalue (atomic mass Burden matrix)
%     BCUT2D_MWLOW  - Min eigenvalue (atomic mass Burden matrix)
%     BCUT2D_CHGHI  - Max eigenvalue (Gasteiger charge Burden matrix)
%     BCUT2D_CHGLO  - Min eigenvalue (Gasteiger charge Burden matrix)
%     BCUT2D_LOGPHI - Max eigenvalue (Crippen LogP Burden matrix)
%     BCUT2D_LOGPLOW - Min eigenvalue (Crippen LogP Burden matrix)
%     BCUT2D_MRHI   - Max eigenvalue (Crippen molar refractivity Burden matrix)
%     BCUT2D_MRLOW  - Min eigenvalue (Crippen molar refractivity Burden matrix)
%
%   Arguments:
%     mol - py.rdkit.Chem.rdchem.Mol  RDKit molecule object
%
%   Returns:
%     scores - double(1,8)  BCUT2D descriptor values in the order above.
%
%   Errors:
%     emk:descriptor:bcut:invalidInput - mol is not a Mol object
%     emk:descriptor:bcut:rdkitError   - unexpected Python exception
%
%   Example:
%     mol    = emk.mol.fromSmiles("c1ccccc1");  % benzene
%     scores = emk.descriptor.bcut(mol);
%     names  = ["BCUT2D_MWHI","BCUT2D_MWLOW","BCUT2D_CHGHI","BCUT2D_CHGLO", ...
%               "BCUT2D_LOGPHI","BCUT2D_LOGPLOW","BCUT2D_MRHI","BCUT2D_MRLOW"];
%     for k = 1:8
%         fprintf("%s: %.4f\n", names(k), scores(k));
%     end
%
%   References:
%     Burden, F.R. (1989). Molecular Identification Number for Substructure
%       Searches. J. Chem. Inf. Comput. Sci. 29(3):225-227.
%       DOI: 10.1021/ci00063a011
%     RDKit Documentation: rdkit.Chem.rdMolDescriptors.BCUT2D
%
%   See also: emk.descriptor.qed, emk.descriptor.fragmentCount,
%             emk.mol.computeDescriptors

    % --- Validate mol ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:descriptor:bcut:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    try
        descMod = py.importlib.import_module("rdkit.Chem.rdMolDescriptors");
        % BCUT2D returns a tuple of 6 floats
        pyResult = descMod.BCUT2D(mol);
        scores   = double(py.array.array("d", pyResult));
    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:descriptor:bcut:rdkitError", ...
            "BCUT2D calculation failed: %s", ME.message);
    end

    % Ensure 1x6 row vector
    scores = reshape(scores, 1, []);

    logDebug("descriptor.bcut: computed %d BCUT2D descriptors", numel(scores));
end
