function mol3d = embed(mol, varargin)
% embed  Embed a molecule in 3D space using the ETKDG algorithm.
%
%   mol3d = emk.conformer.embed(mol)
%   mol3d = emk.conformer.embed(mol, RandomSeed=42)
%   mol3d = emk.conformer.embed(mol, Method="ETKDGv3")
%
%   Generates a single 3D conformer for the input molecule using RDKit's
%   ETKDG (Experimental Torsion Knowledge Distance Geometry) algorithm.
%   Hydrogen atoms are added before embedding and removed from the output
%   (implicit H retained), but the 3D coordinates are optimised in the
%   H-rich environment for better geometry.
%
%   The returned mol3d is the input molecule with one 3D conformer attached.
%   Use emk.conformer.optimize to further refine coordinates with a force
%   field, or emk.shape.compare to compute 3D shape similarity.
%
%   Arguments:
%     mol        - py.rdkit.Chem.rdchem.Mol  2D RDKit molecule (Python ref).
%     Method     - string (optional, default "ETKDGv3")
%                  "ETKDGv3" (recommended), "ETKDGv2", "ETKDG", "KDG".
%     RandomSeed - integer (optional, default -1 = random)
%                  Fixed seed for reproducible embedding.  Use an integer
%                  >= 0 for reproducibility.
%
%   Returns:
%     mol3d - py.rdkit.Chem.rdchem.Mol  Molecule with a 3D conformer.
%             Hydrogen atoms are stripped (RemoveHs).
%             Has exactly one conformer (conformer index 0).
%
%   Errors:
%     emk:conformer:embed:invalidInput     - mol is not a Mol object
%     emk:conformer:embed:invalidMethod    - unsupported embedding method
%     emk:conformer:embed:embeddingFailed  - no conformer could be generated
%     emk:conformer:embed:rdkitError       - unexpected Python exception
%
%   Example:
%     mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
%     mol3d = emk.conformer.embed(mol, "RandomSeed", 42);
%     mol3d = emk.conformer.optimize(mol3d);
%     % Use mol3d in emk.shape.compare or draw with emk.viz.draw2d
%
%   References:
%     Wang, S. et al. (2020). Improving Conformer Generation for Small
%       Molecules: Learning Torsional Distributions from the Cambridge
%       Structural Database. J. Chem. Inf. Model. 60(4):2044-2058.
%       DOI: 10.1021/acs.jcim.0c00025  (ETKDGv3)
%     Riniker, S. & Landrum, G.A. (2015). Better Informed Distance
%       Geometry: Using What We Know To Improve Conformer Generation.
%       J. Chem. Inf. Model. 55(12):2562-2574. DOI: 10.1021/acs.jcim.5b00654
%     RDKit Documentation: AllChem.EmbedMolecule, EmbedParameters
%
%   See also: emk.conformer.optimize, emk.shape.compare, emk.mol.fromSmiles

    VALID_METHODS = ["ETKDGv3", "ETKDGv2", "ETKDG", "KDG"];

    % --- Parse optional name-value arguments ---
    method     = "ETKDGv3";
    randomSeed = -1;
    if nargin > 1
        if mod(numel(varargin), 2) ~= 0
            error("emk:conformer:embed:invalidInput", ...
                "Optional arguments must be name-value pairs.");
        end
        for k = 1:2:numel(varargin)
            argName = string(varargin{k});
            switch lower(argName)
                case "method"
                    method = string(varargin{k+1});
                case "randomseed"
                    randomSeed = varargin{k+1};
                otherwise
                    error("emk:conformer:embed:invalidInput", ...
                        "Unknown option: '%s'. Supported: Method, RandomSeed.", ...
                        varargin{k});
            end
        end
    end

    % --- Validate mol ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:conformer:embed:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    % --- Validate Method (case-insensitive, canonicalise to stored form) ---
    idx = strcmpi(method, VALID_METHODS);
    if ~any(idx)
        error("emk:conformer:embed:invalidMethod", ...
            "Method must be one of [%s], got: %s", ...
            strjoin(VALID_METHODS, ", "), method);
    end
    method = VALID_METHODS(idx);

    logDebug("conformer.embed: embedding with %s (seed=%d)", method, randomSeed);

    try
        allchem = py.importlib.import_module("rdkit.Chem.AllChem");

        % Add explicit Hs for better 3D geometry
        molH = allchem.AddHs(mol);

        % Choose embedding parameters
        switch method
            case "ETKDGv3"
                params = allchem.ETKDGv3();
            case "ETKDGv2"
                params = allchem.ETKDGv2();
            case "ETKDG"
                params = allchem.ETKDG();
            case "KDG"
                params = allchem.KDG();
        end

        % Set random seed for reproducibility
        if randomSeed >= 0
            params.randomSeed = int32(randomSeed);
        end

        % EmbedMolecule modifies molH in-place and returns 0 on success, -1 on fail
        result = double(allchem.EmbedMolecule(molH, params));
        if result == -1
            error("emk:conformer:embed:embeddingFailed", ...
                "EmbedMolecule could not generate a conformer for this molecule. " + ...
                "This can happen for very large, highly rigid, or unusual structures.");
        end

        % Remove Hs (but retain the 3D coordinates on heavy atoms)
        mol3d = allchem.RemoveHs(molH);

    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:conformer:embed:rdkitError", ...
            "3D embedding failed: %s", ME.message);
    end

    nConformers = double(mol3d.GetNumConformers());
    logDebug("conformer.embed: done, mol3d has %d conformer(s)", nConformers);
end
