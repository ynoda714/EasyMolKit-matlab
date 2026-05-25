function tf = hasSubstruct(mol, query)
% hasSubstruct  Test whether molecule(s) contain a substructure.
%
%   tf = emk.mol.hasSubstruct(mol, query)
%
%   Tests whether one or more molecules contain a substructure defined by a
%   SMARTS pattern string or a query Mol object.  Uses RDKit
%   HasSubstructMatch() internally.
%
%   When mol is a cell array of Mol objects (e.g. from emk.io.readSdf),
%   a logical row vector of length N is returned -- one result per molecule.
%
%   Arguments:
%     mol    - py.rdkit.Chem.rdchem.Mol  |  1-by-N cell array of same
%              Target molecule(s). Cell array returns a 1-by-N logical array.
%     query  - string | char | py.rdkit.Chem.rdchem.Mol
%              Substructure query. Strings are interpreted as SMARTS patterns
%              via MolFromSmarts (supports extended SMARTS notation).
%              Pass a Mol object directly for exact subgraph queries.
%
%   Returns:
%     tf     - logical scalar (single mol) or logical(1,N) row vector (cell)
%              true if the molecule contains the substructure; false otherwise.
%
%   Errors:
%     emk:mol:hasSubstruct:invalidMol    - mol is not a Mol or cell of Mols
%     emk:mol:hasSubstruct:invalidQuery  - query string cannot be parsed
%     emk:mol:hasSubstruct:rdkitError    - unexpected Python exception
%
%   Example:
%     mol  = emk.mol.fromSmiles("c1ccccc1CO");  % benzyl alcohol
%     tf   = emk.mol.hasSubstruct(mol, "c1ccccc1");   % true (benzene ring)
%     tf2  = emk.mol.hasSubstruct(mol, "[NH2]");      % false (no amine)
%
%     mols = {emk.mol.fromSmiles("CCO"), emk.mol.fromSmiles("c1ccccc1")};
%     hits = emk.mol.hasSubstruct(mols, "c1ccccc1"); % [false, true]
%
%   See also: emk.mol.fromSmiles, emk.io.readSdf

    % --- Input validation: mol ---
    isSingleMol = isa(mol, "py.rdkit.Chem.rdchem.Mol");
    isCellOfMol = iscell(mol) && ~isempty(mol) && ...
                  all(cellfun(@(m) isa(m, "py.rdkit.Chem.rdchem.Mol"), mol));

    if ~isSingleMol && ~isCellOfMol
        error("emk:mol:hasSubstruct:invalidMol", ...
            "mol must be a py.rdkit.Chem.rdchem.Mol or cell array of same, got: %s", ...
            class(mol));
    end

    % --- Input validation: query ---
    if ischar(query) || isStringScalar(query)
        % Parse SMARTS string into a query Mol object
        queryStr = string(query);
        if strlength(strtrim(queryStr)) == 0
            error("emk:mol:hasSubstruct:invalidQuery", ...
                "query SMARTS string must not be empty or whitespace-only.");
        end
        mods = emk.util.rdkitModule();
        try
            queryMol = mods.Chem.MolFromSmarts(queryStr);
        catch ME
            error("emk:mol:hasSubstruct:rdkitError", ...
                "RDKit MolFromSmarts raised an exception for '%s': %s", ...
                queryStr, ME.message);
        end
        if isa(queryMol, "py.NoneType")
            error("emk:mol:hasSubstruct:invalidQuery", ...
                "SMARTS pattern could not be parsed: '%s'", queryStr);
        end
        logDebug("hasSubstruct: parsed SMARTS query '%s'", queryStr);
    elseif isa(query, "py.rdkit.Chem.rdchem.Mol")
        queryMol = query;
        logDebug("hasSubstruct: received Mol object as query");
    else
        error("emk:mol:hasSubstruct:invalidQuery", ...
            "query must be a SMARTS string or py.rdkit.Chem.rdchem.Mol, got: %s", ...
            class(query));
    end

    % --- Perform substructure match ---
    if isSingleMol
        tf = matchOne_(mol, queryMol);
    else
        n  = numel(mol);
        tf = false(1, n);
        for i = 1:n
            tf(i) = matchOne_(mol{i}, queryMol);
        end
    end

    logDebug("hasSubstruct: done. %d / %d match(es)", sum(tf), numel(tf));
end

% -------------------------------------------------------------------------
% Local helper
% -------------------------------------------------------------------------

function tf = matchOne_(mol, queryMol)
% matchOne_  Call HasSubstructMatch on a single mol.
    try
        result = mol.HasSubstructMatch(queryMol);
        tf     = logical(result);
    catch ME
        error("emk:mol:hasSubstruct:rdkitError", ...
            "RDKit HasSubstructMatch raised an exception: %s", ME.message);
    end
end
