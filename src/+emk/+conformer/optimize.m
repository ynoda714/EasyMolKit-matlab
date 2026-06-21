function mol = optimize(mol, varargin)
% optimize  Optimize a 3D conformer using a molecular force field.
%
%   mol = emk.conformer.optimize(mol3d)
%   mol = emk.conformer.optimize(mol3d, ForceField="MMFF94")
%   mol = emk.conformer.optimize(mol3d, MaxIter=2000)
%
%   Performs force-field minimization on the 3D conformer produced by
%   emk.conformer.embed.  The optimised coordinates are stored in-place
%   in the Python mol object and returned.
%
%   Supported force fields:
%   - MMFF94 (Merck Molecular Force Field 94) -- recommended for drug-like
%     molecules; available for C, H, N, O, S, P, F, Cl, Br, I elements.
%   - UFF   (Universal Force Field) -- fallback for uncommon element types;
%     lower geometric accuracy than MMFF94 for drug-like molecules.
%
%   Arguments:
%     mol        - py.rdkit.Chem.rdchem.Mol  Mol with a 3D conformer
%                  (output of emk.conformer.embed).
%     ForceField - string (optional, default "MMFF94") "MMFF94" or "UFF".
%     MaxIter    - positive integer (optional, default 2000)
%                  Maximum number of minimization steps.
%
%   Returns:
%     mol - py.rdkit.Chem.rdchem.Mol  Same Python object with optimised
%           3D coordinates.  The conformer is modified in-place by RDKit.
%
%   Errors:
%     emk:conformer:optimize:invalidInput      - mol is not a Mol object or
%                                                mol has no 3D conformer
%     emk:conformer:optimize:invalidForceField - unsupported force field name
%     emk:conformer:optimize:optimizeFailed    - force field setup/minimization
%                                                failed (e.g. unsupported element)
%     emk:conformer:optimize:rdkitError        - unexpected Python exception
%
%   Example:
%     mol   = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
%     mol3d = emk.conformer.embed(mol, "RandomSeed", 42);
%     mol3d = emk.conformer.optimize(mol3d);
%     % mol3d now has MMFF94-optimised coordinates
%
%   References:
%     Halgren, T.A. (1996). Merck Molecular Force Field. I. Basis, Form,
%       Scope, Parameterization, and Performance of MMFF94. J. Comput.
%       Chem. 17(5-6):490-519. DOI: 10.1002/jcc.540170510
%     Rappe, A.K. et al. (1992). UFF, a Full Periodic Table Force Field
%       for Molecular Mechanics and Molecular Dynamics Simulations.
%       J. Am. Chem. Soc. 114(25):10024-10035. DOI: 10.1021/ja00051a040
%
%   See also: emk.conformer.embed, emk.shape.compare

    VALID_FF = ["MMFF94", "UFF"];

    % --- Parse optional name-value arguments ---
    forceField = "MMFF94";
    maxIter    = 2000;
    if nargin > 1
        if mod(numel(varargin), 2) ~= 0
            error("emk:conformer:optimize:invalidInput", ...
                "Optional arguments must be name-value pairs.");
        end
        for k = 1:2:numel(varargin)
            argName = string(varargin{k});
            switch lower(argName)
                case "forcefield"
                    forceField = upper(string(varargin{k+1}));
                case "maxiter"
                    maxIter = varargin{k+1};
                otherwise
                    error("emk:conformer:optimize:invalidInput", ...
                        "Unknown option: '%s'. Supported: ForceField, MaxIter.", ...
                        varargin{k});
            end
        end
    end

    % --- Validate mol ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:conformer:optimize:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end
    if double(mol.GetNumConformers()) == 0
        error("emk:conformer:optimize:invalidInput", ...
            "mol has no 3D conformer. Call emk.conformer.embed first.");
    end

    % --- Validate ForceField ---
    if ~ismember(forceField, VALID_FF)
        error("emk:conformer:optimize:invalidForceField", ...
            "ForceField must be one of [%s], got: %s", ...
            strjoin(VALID_FF, ", "), forceField);
    end

    logDebug("conformer.optimize: minimizing with %s (MaxIter=%d)", ...
        forceField, maxIter);

    try
        allchem = py.importlib.import_module("rdkit.Chem.AllChem");

        % Add Hs for accurate force field evaluation, then strip after
        molH = allchem.AddHs(mol, pyargs("addCoords", true));

        switch forceField
            case "MMFF94"
                % MMFFOptimizeMolecule returns 0=success, 1=not converged, -1=error
                status = double(allchem.MMFFOptimizeMolecule(molH, ...
                                    pyargs("maxIters", int32(maxIter))));
            case "UFF"
                status = double(allchem.UFFOptimizeMolecule(molH, ...
                                    pyargs("maxIters", int32(maxIter))));
        end

        if status == -1
            error("emk:conformer:optimize:optimizeFailed", ...
                "%s force field setup failed. " + ...
                "This may be due to unsupported atom types; try ForceField='UFF'.", ...
                forceField);
        end
        if status == 1
            logWarn("conformer.optimize: %s minimization did not converge " + ...
                "(increase MaxIter or check molecule quality)", forceField);
        end

        % Strip Hs again (retain 3D coords on heavy atoms)
        mol = allchem.RemoveHs(molH);

    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:conformer:optimize:rdkitError", ...
            "Force field optimization failed: %s", ME.message);
    end

    logDebug("conformer.optimize: done (%s, status=%d)", forceField, status);
end
