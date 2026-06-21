function tbl = pains(tbl)
% pains  Annotate a descriptor table with PAINS (pan-assay interference) alerts.
%
%   tbl = emk.filter.pains(tbl)
%
%   Applies the PAINS filter (Baell & Holloway 2010) to each row in the
%   descriptor table.  PAINS (Pan-Assay INterference compoundS) are
%   structural motifs known to cause false positives in biochemical assays
%   due to non-specific reactivity, aggregation, redox cycling, or optical
%   interference.
%
%   Requires Python + RDKit (uses RDKit's built-in FilterCatalog).
%
%   Three columns are appended:
%     HasPains          (logical) - true if any PAINS alert matched
%     NumPainsAlerts    (double)  - number of PAINS substructure matches
%     PainsAlerts       (string)  - comma-joined alert names, or "" if none
%
%   The SMILES column (from emk.mol.toTable) is used to reconstruct RDKit
%   mol objects for SMARTS matching.
%
%   Arguments:
%     tbl - (table) Descriptor table with a SMILES column (string).
%           Obtainable from emk.mol.toTable.
%
%   Returns:
%     tbl - (table) Input table with three new columns appended:
%           NumPainsAlerts (double)  - count of matched PAINS alerts
%           PainsAlerts    (string)  - comma-joined alert names
%           HasPains       (logical) - true if NumPainsAlerts > 0
%
%   Errors:
%     emk:filter:pains:invalidInput   - tbl is not a table
%     emk:filter:pains:missingColumns - SMILES column absent
%     emk:filter:pains:rdkitError     - unexpected Python/RDKit exception
%
%   Example:
%     mol  = emk.mol.fromSmiles("O=C(Nc1ccc(Oc2ccc(NC(=O)c3ccc(F)cc3)cc2)cc1)c1ccc(Cl)cc1");
%     tbl  = emk.mol.toTable(mol);
%     tbl  = emk.filter.pains(tbl);
%     disp(tbl.PainsAlerts)
%
%   References:
%     Baell, J.B. & Holloway, G.A. (2010). New Substructure Filters for
%       Removal of Pan Assay Interference Compounds (PAINS) from Screening
%       Libraries and for Their Exclusion in Bioassays. J. Med. Chem.
%       53(7):2719-2740. DOI: 10.1021/jm901137j
%     Baell, J. & Walters, M.A. (2014). Chemical con artists foil drug
%       discovery. Nature 513:481-483. DOI: 10.1038/513481a
%
%   See also: emk.filter.lipinski, emk.filter.veber, emk.filter.reos,
%             emk.mol.toTable

    % --- Validate input table ---
    if ~istable(tbl)
        error("emk:filter:pains:invalidInput", ...
            "Input must be a MATLAB table, got: %s", class(tbl));
    end

    % --- Check SMILES column ---
    if ~ismember("SMILES", tbl.Properties.VariableNames)
        error("emk:filter:pains:missingColumns", ...
            "Table must have a SMILES column. Use emk.mol.toTable to build the table.");
    end

    logDebug("pains: building PAINS FilterCatalog");

    nRows = height(tbl);

    % --- Build SMILES Python list ---
    smilesPyList = py.list();
    for i = 1:nRows
        smilesPyList.append(char(tbl.SMILES(i)));
    end

    % --- Run PAINS via pyrun ---
    % FilterCatalogParams and FilterCatalog are Boost.Python classes;
    % MATLAB cannot instantiate them via py.module.Class() pattern.
    % pyrun() executes Python code in the interpreter's own namespace,
    % bypassing the Boost.Python instantiation restriction.
    try
        pyCode = ...
            "from rdkit.Chem import MolFromSmiles" + newline + ...
            "from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams" + newline + ...
            "p = FilterCatalogParams()" + newline + ...
            "p.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)" + newline + ...
            "cat = FilterCatalog(p)" + newline + ...
            "n_alerts = []" + newline + ...
            "alert_names = []" + newline + ...
            "for smi in smiles_list:" + newline + ...
            "    mol = MolFromSmiles(smi)" + newline + ...
            "    if mol is not None:" + newline + ...
            "        ms = cat.GetMatches(mol)" + newline + ...
            "        n_alerts.append(len(ms))" + newline + ...
            "        names = list(set([m.GetDescription() for m in ms]))" + newline + ...
            "        alert_names.append(', '.join(names))" + newline + ...
            "    else:" + newline + ...
            "        n_alerts.append(0)" + newline + ...
            "        alert_names.append('')";

        [pyNAlerts, pyNames] = pyrun(pyCode, ["n_alerts", "alert_names"], ...
            smiles_list=smilesPyList);
    catch ME
        error("emk:filter:pains:rdkitError", ...
            "Failed to apply PAINS FilterCatalog: %s", ME.message);
    end

    % --- Convert Python list results to MATLAB arrays ---
    nAlerts    = zeros(nRows, 1);
    alertNames = strings(nRows, 1);
    for i = 1:nRows
        nAlerts(i)    = double(pyNAlerts{i});
        alertNames(i) = string(pyNames{i});
    end

    tbl.NumPainsAlerts = nAlerts;
    tbl.PainsAlerts    = alertNames;
    tbl.HasPains       = logical(nAlerts > 0);

    logInfo("pains: %d / %d row(s) have PAINS alerts", ...
        sum(tbl.HasPains), nRows);
end
