function [tbl, unmatchedIdx] = rgroup(mols, coreSmiles)
% rgroup  R-group decomposition: decompose molecules into core and R groups.
%
%   [tbl, unmatchedIdx] = emk.scaffold.rgroup(mols, coreSmiles)
%
%   Decomposes a set of molecules sharing a common core into their
%   core scaffold and variable R-group substituents using RDKit's
%   RGroupDecompose algorithm.  The decomposition identifies the
%   common substructure (core) in each molecule and labels the
%   attachment points as R1, R2, ... in SMILES notation.
%
%   Arguments:
%     mols       - 1xN cell array of py.rdkit.Chem.rdchem.Mol objects
%     coreSmiles - string or char  SMILES (or SMARTS) of the common core.
%                  SMARTS is attempted first; SMILES as fallback.
%
%   Returns:
%     tbl          - table  One row per matched molecule.  Columns:
%                    Core (string): core SMILES with R-group labels [*:k]
%                    R1 (string), R2 (string), ...: R-group SMILES
%                    Unmatched rows are excluded (see unmatchedIdx).
%     unmatchedIdx - double(1,K)  1-based indices into mols that did not
%                    match the core.  Empty if all molecules matched.
%
%   Errors:
%     emk:scaffold:rgroup:invalidInput - mols is not a non-empty cell array
%                                        of Mol objects, or coreSmiles is invalid
%     emk:scaffold:rgroup:invalidCore  - coreSmiles could not be parsed
%     emk:scaffold:rgroup:rdkitError   - unexpected Python exception
%     emk:scaffold:rgroup:noMatch      - no molecules matched the core
%
%   Example:
%     smiles = {"Cc1ccc(N)cc1", "Clc1ccc(N)cc1", "Cc1ccc(O)cc1"};
%     mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
%     [tbl, unmatched] = emk.scaffold.rgroup(mols, "c1ccccc1");
%     % tbl has columns Core, R1, R2 with R-group SMILES
%     % tbl.R1 might contain "[*:1]C", "[*:1]Cl", "[*:1]C"
%     % tbl.R2 might contain "[*:2]N", "[*:2]N", "[*:2]O"
%
%   References:
%     RDKit Documentation: rdkit.Chem.rdRGroupDecomposition.RGroupDecompose
%     Degen, J. & Rarey, M. (2006). FlexNovo: Structure-Based Searching
%       in Multidimensional Chemical Space. ChemMedChem 1(8):854-868.
%       DOI: 10.1002/cmdc.200600073
%
%   See also: emk.scaffold.brics, emk.scaffold.genericMurcko,
%             emk.mol.hasSubstruct, emk.mol.toSmiles

    % --- Validate mols ---
    if ~iscell(mols) || isempty(mols)
        error("emk:scaffold:rgroup:invalidInput", ...
            "mols must be a non-empty cell array, got: %s", class(mols));
    end
    for i = 1:numel(mols)
        if ~isa(mols{i}, "py.rdkit.Chem.rdchem.Mol")
            error("emk:scaffold:rgroup:invalidInput", ...
                "mols{%d} must be a py.rdkit.Chem.rdchem.Mol, got: %s", ...
                i, class(mols{i}));
        end
    end

    % --- Validate coreSmiles ---
    if ~(ischar(coreSmiles) || isStringScalar(coreSmiles))
        error("emk:scaffold:rgroup:invalidInput", ...
            "coreSmiles must be a string scalar or char, got: %s", class(coreSmiles));
    end
    coreSmiles = string(coreSmiles);
    if strlength(coreSmiles) == 0
        error("emk:scaffold:rgroup:invalidInput", ...
            "coreSmiles must not be empty");
    end

    logDebug("rgroup: decomposing %d molecule(s) with core '%s'", ...
        numel(mols), coreSmiles);

    mods = emk.util.rdkitModule();

    % --- Parse core: try SMARTS first, then SMILES ---
    try
        core = mods.Chem.MolFromSmarts(coreSmiles);
        if isa(core, "py.NoneType")
            core = mods.Chem.MolFromSmiles(coreSmiles);
        end
        if isa(core, "py.NoneType")
            error("emk:scaffold:rgroup:invalidCore", ...
                "Could not parse coreSmiles as SMARTS or SMILES: '%s'", coreSmiles);
        end
    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:scaffold:rgroup:rdkitError", ...
            "Core parsing failed: %s", ME.message);
    end

    % --- Build Python list of molecules ---
    pyMols = py.list();
    for i = 1:numel(mols)
        pyMols.append(mols{i});
    end
    pyCores = py.list({core});

    % --- Call RGroupDecompose ---
    try
        rgMod = py.importlib.import_module( ...
            "rdkit.Chem.rdRGroupDecomposition");

        % RGroupDecompose returns a 2-tuple (rows, unmatched_indices)
        % rows:    Python list of dicts {colname -> SMILES string}
        % unmatched: Python list of 0-based int indices
        result    = rgMod.RGroupDecompose(pyCores, pyMols, ...
                        pyargs("asSmiles", true));
        pyRows    = result{1};
        pyUnmatched = result{2};
    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:scaffold:rgroup:rdkitError", ...
            "RGroupDecompose failed: %s", ME.message);
    end

    % --- Convert unmatched 0-based indices to 1-based MATLAB indices ---
    nUnmatched = double(py.len(pyUnmatched));
    unmatchedIdx = zeros(1, nUnmatched);
    for i = 1:nUnmatched
        unmatchedIdx(i) = double(pyUnmatched{i}) + 1;
    end

    % --- Check that at least one molecule matched ---
    nRows = double(py.len(pyRows));
    if nRows == 0
        error("emk:scaffold:rgroup:noMatch", ...
            "No molecules matched the core '%s'. " + ...
            "Check coreSmiles or verify that mols share this scaffold.", coreSmiles);
    end

    % --- Determine column names from the first row ---
    firstRow  = pyRows{1};
    keysView  = firstRow.keys();
    keysList  = py.list(keysView);
    nCols     = double(py.len(keysList));
    colNames  = strings(1, nCols);
    for k = 1:nCols
        colNames(k) = string(keysList{k});
    end

    % Sort: Core first, then R1, R2, ... in numeric order
    iCore  = find(colNames == "Core", 1);
    rCols  = sort(colNames(colNames ~= "Core"));
    if ~isempty(iCore)
        colNames = ["Core", rCols];
    else
        colNames = rCols;
    end

    % --- Build MATLAB table from Python rows ---
    % Pre-allocate string columns
    dataCell = cell(nRows, numel(colNames));
    for r = 1:nRows
        rowDict = pyRows{r};
        for c = 1:numel(colNames)
            val = rowDict.get(colNames(c));
            if isa(val, "py.NoneType") || isempty(val)
                dataCell{r, c} = "";
            else
                dataCell{r, c} = string(val);
            end
        end
    end

    % Build table with string columns
    tbl = cell2table(dataCell, "VariableNames", cellstr(colNames));

    logDebug("rgroup: %d matched, %d unmatched, %d R-group column(s)", ...
        nRows, nUnmatched, numel(colNames) - 1);
end
