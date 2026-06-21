function score = compare(mol1, mol2, varargin)
% compare  3D shape similarity between two molecules with conformers.
%
%   score = emk.shape.compare(mol1, mol2)
%   score = emk.shape.compare(mol1, mol2, Method="protrude")
%
%   Computes a 3D shape similarity score between two molecules.  Both
%   molecules must have at least one conformer (see emk.conformer.embed).
%
%   Available methods:
%   - "protrude" (default): Uses rdkit.Chem.rdShapeHelpers.ShapeProtrudeDist.
%     Returns score = 1 - ProtrududeDistance.  The protrude distance
%     measures the fraction of mol1 volume NOT overlapping with mol2.
%     Score range: [0, 1].  1 = identical shapes, 0 = no overlap.
%   - "tanimoto": Uses rdkit.Chem.rdShapeHelpers.ShapeTverskyIndex with
%     alpha=beta=1 (equivalent to shape Tanimoto).
%     Score range: [0, 1].
%
%   NOTE: Shape comparison depends critically on 3D conformer quality.
%   Use emk.conformer.embed + emk.conformer.optimize before comparing.
%   Results may vary with conformer alignment; mol2 is used as the reference.
%
%   Arguments:
%     mol1   - py.rdkit.Chem.rdchem.Mol  Probe molecule (with 3D conformer)
%     mol2   - py.rdkit.Chem.rdchem.Mol  Reference molecule (with 3D conformer)
%     Method - string (optional, default "protrude") "protrude" or "tanimoto"
%
%   Returns:
%     score - double in [0, 1]  Shape similarity score.
%             Higher values indicate greater 3D shape similarity.
%
%   Errors:
%     emk:shape:compare:invalidInput  - mol1/mol2 not Mol objects or
%                                       missing 3D conformers
%     emk:shape:compare:invalidMethod - unsupported method name
%     emk:shape:compare:rdkitError    - unexpected Python exception
%
%   Example:
%     mol1   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
%     mol2   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1");         % phenyl acetate
%     mol1_3d = emk.conformer.optimize(emk.conformer.embed(mol1, "RandomSeed", 0));
%     mol2_3d = emk.conformer.optimize(emk.conformer.embed(mol2, "RandomSeed", 0));
%     score = emk.shape.compare(mol1_3d, mol2_3d);
%     fprintf("Shape similarity: %.3f\n", score);
%
%   References:
%     Ballester, P.J. & Richards, W.G. (2007). Ultrafast Shape Recognition
%       to Search Compound Libraries for Similar Molecular Shapes. J.
%       Comput. Chem. 28(10):1711-1723. DOI: 10.1002/jcc.20681
%     Grant, J.A. & Pickup, B.T. (1995). A Gaussian Description of
%       Molecular Shape. J. Phys. Chem. 99(11):3503-3510.
%       DOI: 10.1021/j100011a016
%     RDKit Documentation: rdkit.Chem.rdShapeHelpers
%
%   See also: emk.conformer.embed, emk.conformer.optimize, emk.mol.fromSmiles

    VALID_METHODS = ["protrude", "tanimoto"];

    % --- Parse optional name-value arguments ---
    method = "protrude";
    if nargin > 2
        if mod(numel(varargin), 2) ~= 0
            error("emk:shape:compare:invalidInput", ...
                "Optional arguments must be name-value pairs.");
        end
        for k = 1:2:numel(varargin)
            argName = string(varargin{k});
            switch lower(argName)
                case "method"
                    method = lower(string(varargin{k+1}));
                otherwise
                    error("emk:shape:compare:invalidInput", ...
                        "Unknown option: '%s'. Supported: Method.", varargin{k});
            end
        end
    end

    % --- Validate mol1 and mol2 ---
    for molName = {"mol1", "mol2"}
        switch molName{1}
            case "mol1"; m = mol1;
            case "mol2"; m = mol2;
        end
        if ~isa(m, "py.rdkit.Chem.rdchem.Mol")
            error("emk:shape:compare:invalidInput", ...
                "%s must be a py.rdkit.Chem.rdchem.Mol, got: %s", ...
                molName{1}, class(m));
        end
        if double(m.GetNumConformers()) == 0
            error("emk:shape:compare:invalidInput", ...
                "%s has no 3D conformer. Call emk.conformer.embed first.", ...
                molName{1});
        end
    end

    % --- Validate Method ---
    if ~ismember(method, VALID_METHODS)
        error("emk:shape:compare:invalidMethod", ...
            "Method must be one of [%s], got: %s", ...
            strjoin(VALID_METHODS, ", "), method);
    end

    logDebug("shape.compare: computing 3D shape similarity (method=%s)", method);

    try
        shapeHelpers = py.importlib.import_module( ...
            "rdkit.Chem.rdShapeHelpers");

        switch method
            case "protrude"
                % ShapeProtrudeDist returns distance in [0,1]: 0=same shape
                dist  = double(shapeHelpers.ShapeProtrudeDist(mol1, mol2, ...
                            pyargs("allowReordering", false)));
                score = 1.0 - dist;

            case "tanimoto"
                % ShapeTverskyIndex(mol1, mol2, alpha, beta)
                % With alpha=beta=1 gives Shape Tanimoto
                score = double(shapeHelpers.ShapeTverskyIndex( ...
                            mol1, mol2, 1.0, 1.0));
        end

    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:shape:compare:rdkitError", ...
            "3D shape comparison failed: %s", ME.message);
    end

    % Clamp to [0,1] for numerical safety
    score = max(0.0, min(1.0, score));

    logDebug("shape.compare: score=%.4f (method=%s)", score, method);
end
