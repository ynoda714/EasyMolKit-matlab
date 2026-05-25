function varargout = draw2d(mol, opts)
% draw2d  Draw a 2D molecular structure in a MATLAB figure.
%
%   emk.viz.draw2d(mol)
%   emk.viz.draw2d(mol, Title="Aspirin")
%   emk.viz.draw2d(mol, Width=400, Height=300)
%   fig = emk.viz.draw2d(mol)
%
%   Generates a 2D depiction of a molecule using RDKit's drawing engine and
%   displays the result in the CURRENT axes (respects figure/subplot context
%   set up by the caller).  If no figure exists, one is created automatically.
%   2D coordinates are computed via AllChem.Compute2DCoords before rendering.
%
%   Arguments:
%     mol     - py.rdkit.Chem.rdchem.Mol  (from emk.mol.fromSmiles)
%     Title   - (string, optional, default "") Figure title text.
%               Empty string means no title is added.
%     Width   - (double, optional, default 300) Image width in pixels.
%     Height  - (double, optional, default 300) Image height in pixels.
%
%   Returns:
%     fig     - (optional) matlab.ui.Figure handle of the current figure.
%
%   Errors:
%     emk:viz:draw2d:invalidInput  - mol is not a py.rdkit.Chem.rdchem.Mol
%     emk:viz:draw2d:rdkitError    - RDKit drawing or I/O failed
%
%   Example:
%     mol = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
%     emk.viz.draw2d(mol, Title="Aspirin");
%
%   See also: emk.mol.fromSmiles, emk.descriptor.calculate

    arguments
        mol
        opts.Title  string = ""
        opts.Width  double = 300
        opts.Height double = 300
    end

    % --- Input validation (RDKit not required) ---
    if ~isa(mol, "py.rdkit.Chem.rdchem.Mol")
        error("emk:viz:draw2d:invalidInput", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol, got: %s", class(mol));
    end

    if opts.Width < 1 || opts.Height < 1
        error("emk:viz:draw2d:invalidInput", ...
            "Width and Height must be >= 1 (got Width=%d, Height=%d)", ...
            opts.Width, opts.Height);
    end

    logDebug("draw2d: generating 2D image (%dx%d)", opts.Width, opts.Height);

    % --- Generate PNG via RDKit ---
    tmpFile = string(tempname()) + ".png";
    try
        allchem = py.importlib.import_module("rdkit.Chem.AllChem");
        draw    = py.importlib.import_module("rdkit.Chem.Draw");

        % Compute 2D coordinates (idempotent; required before rendering)
        allchem.Compute2DCoords(mol);

        % Write PNG directly to temp file via Draw.MolToFile
        imgSize = py.tuple({int32(opts.Width), int32(opts.Height)});
        draw.MolToFile(mol, tmpFile, pyargs("size", imgSize));
    catch ME
        if isfile(tmpFile)
            delete(tmpFile);
        end
        error("emk:viz:draw2d:rdkitError", ...
            "RDKit draw2d failed: %s", ME.message);
    end

    % --- Read PNG with MATLAB imread ---
    try
        img = imread(tmpFile);
    catch ME
        if isfile(tmpFile)
            delete(tmpFile);
        end
        error("emk:viz:draw2d:rdkitError", ...
            "Failed to read generated PNG image: %s", ME.message);
    end
    delete(tmpFile);

    % --- Display in the current axes (caller controls figure/subplot context) ---
    % Using gca() respects figure("Name",...) and subplot() set up by the caller.
    % If no figure exists yet, gca() creates a default one automatically.
    ax = gca();
    imshow(img, "Parent", ax);
    axis(ax, "off");
    if strlength(opts.Title) > 0
        title(ax, opts.Title, "Interpreter", "none");
    end

    if nargout > 0
        varargout{1} = ancestor(ax, "figure");
    end
end
