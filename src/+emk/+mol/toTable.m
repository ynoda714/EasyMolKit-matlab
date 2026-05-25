function tbl = toTable(mols, varargin)
% toTable  Build a MATLAB table from Mol objects with SMILES and descriptors.
%
%   tbl = emk.mol.toTable(mols)
%   tbl = emk.mol.toTable(mols, Properties=["SMILES","MolWt","LogP"])
%
%   Combines SMILES conversion and descriptor calculation into a single
%   MATLAB table.  This is the primary aggregate-table API for end users
%   (ADR-002 rev.4 Aggregate Table Builder pattern).
%
%   A single Mol object is automatically wrapped in a cell array.
%   For each molecule, SMILES is retrieved via emk.mol.toSmiles and
%   descriptors via emk.descriptor.batchCalculate (single IPC round-trip
%   per ADR-002 rev.3).  Invalid Mol entries produce SMILES = "<invalid>"
%   and NaN in all descriptor columns.
%
%   Arguments:
%     mols       - (1xN cell | py.rdkit.Chem.rdchem.Mol)
%                  Cell array of Mol objects, or a single Mol object.
%     Properties - (string array, name-value, optional)
%                  Column names to include.  "SMILES" and any descriptor
%                  name supported by emk.descriptor.calculate are valid.
%                  Default: ["SMILES", all 10 descriptors].
%                  Order is preserved in the output table.
%
%   Returns:
%     tbl  - (table) N-row table.  Columns match Properties.
%            SMILES column: string scalar per row.
%            Descriptor columns: double scalar per row.
%            Rows for invalid mols: SMILES = "<invalid>", descriptors = NaN.
%
%   Supported Properties:
%     "SMILES"           - Canonical SMILES string (from emk.mol.toSmiles)
%     "MolWt"            - Average molecular weight (g/mol)
%     "ExactMolWt"       - Monoisotopic molecular weight (g/mol)
%     "LogP"             - Wildman-Crippen LogP
%     "TPSA"             - Topological polar surface area (A^2)
%     "NumHAcceptors"    - H-bond acceptor count
%     "NumHDonors"       - H-bond donor count
%     "NumRotatableBonds"- Rotatable bond count
%     "RingCount"        - Total ring count
%     "FractionCSP3"     - Fraction of sp3 carbons
%     "HeavyAtomCount"   - Heavy atom count
%
%   Errors:
%     emk:mol:toTable:invalidInput     - mols is not a cell or Mol object
%     emk:mol:toTable:unknownProperty  - unrecognised property name(s)
%     emk:mol:toTable:emptyProperties  - Properties argument is empty
%
%   Example:
%     mols = {emk.mol.fromSmiles("CCO"), emk.mol.fromSmiles("c1ccccc1")};
%     tbl  = emk.mol.toTable(mols);
%     tbl  = emk.mol.toTable(mols, Properties=["SMILES","MolWt","LogP"]);
%     % Single mol:
%     mol  = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");  % aspirin
%     tbl  = emk.mol.toTable(mol);
%
%   See also: emk.mol.fromSmiles, emk.mol.toSmiles,
%             emk.descriptor.calculate, emk.descriptor.batchCalculate

    SUPPORTED_DESCS = ["MolWt", "ExactMolWt", "LogP", "TPSA", ...
                       "NumHAcceptors", "NumHDonors", "NumRotatableBonds", ...
                       "RingCount", "FractionCSP3", "HeavyAtomCount"];
    ALL_PROPS = ["SMILES", SUPPORTED_DESCS];

    % --- Parse optional Properties name-value argument ---
    % Use a simple parser so that validation fires before any RDKit call.
    props = ALL_PROPS;   % default
    if nargin > 1
        if mod(numel(varargin), 2) ~= 0
            error("emk:mol:toTable:invalidInput", ...
                "Optional arguments must be name-value pairs.");
        end
        for k = 1:2:numel(varargin)
            argName = varargin{k};
            if ~(ischar(argName) || isStringScalar(argName))
                error("emk:mol:toTable:invalidInput", ...
                    "Argument name must be a string scalar, got: %s", class(argName));
            end
            if strcmpi(argName, "Properties")
                props = varargin{k+1};
                if ischar(props)
                    props = string(props);
                end
                props = reshape(string(props), 1, []);
            else
                error("emk:mol:toTable:invalidInput", ...
                    "Unknown argument name: '%s'. Valid names: Properties", argName);
            end
        end
    end

    % --- Validate Properties before any RDKit call ---
    if isempty(props)
        error("emk:mol:toTable:emptyProperties", ...
            "Properties must not be empty.");
    end
    unknown = props(~ismember(props, ALL_PROPS));
    if ~isempty(unknown)
        error("emk:mol:toTable:unknownProperty", ...
            "Unknown property name(s): [%s]. Supported: [%s]", ...
            strjoin(unknown, ", "), strjoin(ALL_PROPS, ", "));
    end

    % --- Normalize mols input ---
    % Accept a single Mol object: wrap in a cell for uniform processing.
    if ~iscell(mols)
        if startsWith(class(mols), "py.") && contains(class(mols), "Mol")
            mols = {mols};
        else
            error("emk:mol:toTable:invalidInput", ...
                "mols must be a cell array or a single Mol object, got: %s", class(mols));
        end
    end

    nMols = numel(mols);
    logDebug("toTable: %d mol(s), properties=[%s]", nMols, strjoin(props, ", "));

    % --- Identify which column types are requested ---
    wantSmiles = any(props == "SMILES");
    descProps  = props(props ~= "SMILES");

    % --- Compute SMILES column ---
    if wantSmiles
        smilesCol = strings(nMols, 1);
        for i = 1:nMols
            try
                smilesCol(i) = emk.mol.toSmiles(mols{i});
            catch
                smilesCol(i) = "<invalid>";
            end
        end
    end

    % --- Compute descriptor columns via batchCalculate ---
    if ~isempty(descProps)
        descTbl = emk.descriptor.batchCalculate(mols, descProps);
    end

    % --- Assemble output table in requested column order ---
    if nMols == 0
        % Return an empty table with the correct variable names and types.
        varArrays    = cell(1, numel(props));
        varNames     = cell(1, numel(props));
        for p = 1:numel(props)
            varNames{p} = char(props(p));
            if props(p) == "SMILES"
                varArrays{p} = strings(0, 1);
            else
                varArrays{p} = zeros(0, 1);
            end
        end
        tbl = table(varArrays{:}, 'VariableNames', varNames);
        return;
    end

    varArrays = cell(1, numel(props));
    varNames  = cell(1, numel(props));
    for p = 1:numel(props)
        varNames{p} = char(props(p));
        if props(p) == "SMILES"
            varArrays{p} = smilesCol;
        else
            varArrays{p} = descTbl.(props(p));
        end
    end

    tbl = table(varArrays{:}, 'VariableNames', varNames);

    logDebug("toTable: returned table %dx%d", height(tbl), width(tbl));
end
