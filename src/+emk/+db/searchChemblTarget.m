function result = searchChemblTarget(query, options)
% searchChemblTarget  Search ChEMBL for protein targets by preferred name.
%
%   result = emk.db.searchChemblTarget(query)
%   result = emk.db.searchChemblTarget(query, TargetType="SINGLE PROTEIN")
%   result = emk.db.searchChemblTarget(query, MaxRows=10)
%
%   Queries the ChEMBL REST API target endpoint for protein targets whose
%   preferred name contains the given query string (case-insensitive).
%   Python is not required.  Returns a MATLAB table.
%
%   ChEMBL REST endpoint used:
%     https://www.ebi.ac.uk/chembl/api/data/target.json
%       ?pref_name__icontains={query}&target_type={TargetType}&limit={MaxRows}
%     Note: spaces in TargetType are encoded as %20 (not +) to avoid HTTP 500.
%
%   Arguments:
%     query      - string | char, partial preferred name (e.g.
%                  "Epidermal growth factor receptor", "cyclooxygenase").
%     TargetType - string, ChEMBL target_type filter.
%                  Default: "SINGLE PROTEIN".
%                  Set to "" to return all target types.
%     MaxRows    - positive integer, max targets to return. Default: 10.
%
%   Returns:
%     result  - MATLAB table with columns:
%       TargetChEMBLID  string  ChEMBL target ID (e.g. "CHEMBL203")
%       PreferredName   string  Target preferred name
%       Organism        string  Organism (e.g. "Homo sapiens")
%       TargetType      string  Target classification (e.g. "SINGLE PROTEIN")
%
%   Errors:
%     emk:db:searchChemblTarget:invalidInput   - query is not a non-empty string/char
%     emk:db:searchChemblTarget:invalidOptions - invalid MaxRows value
%     emk:db:searchChemblTarget:notFound       - no targets match the query
%     emk:db:searchChemblTarget:networkError   - webread failure or bad response
%
%   Example:
%     tbl = emk.db.searchChemblTarget("Epidermal growth factor receptor");
%     tbl = emk.db.searchChemblTarget("cyclooxygenase", MaxRows=5);
%     tbl = emk.db.searchChemblTarget("kinase", TargetType="", MaxRows=20);
%
%   See also: emk.db.getChemblActivity, emk.db.searchChembl

    arguments
        query
        options.TargetType (1,1) string = "SINGLE PROTEIN"
        options.MaxRows    (1,1) double = 10
    end

    % --- Input validation: query ---
    if ~(ischar(query) || isStringScalar(query))
        error("emk:db:searchChemblTarget:invalidInput", ...
            "query must be a string or char, got: %s", class(query));
    end
    query = string(query);
    if strlength(strtrim(query)) == 0
        error("emk:db:searchChemblTarget:invalidInput", ...
            "query must not be empty or whitespace-only");
    end

    % --- Input validation: MaxRows ---
    if ~(isscalar(options.MaxRows) && options.MaxRows > 0 && ...
         floor(options.MaxRows) == options.MaxRows && isfinite(options.MaxRows))
        error("emk:db:searchChemblTarget:invalidOptions", ...
            "MaxRows must be a positive finite integer, got: %g", options.MaxRows);
    end
    maxRows = floor(options.MaxRows);

    % --- Build ChEMBL REST URL ---
    % target_type must be encoded with %20 (not +) for spaces.
    % urlencode() produces '+', which triggers HTTP 500 on ChEMBL API.
    % We replace '+' with '%20' for this parameter only.
    % Client-side filtering is kept as a safety net.
    BASE_URL     = "https://www.ebi.ac.uk/chembl/api/data";
    encodedQuery = string(urlencode(char(query)));
    filterTargetType = strlength(options.TargetType) > 0;
    url = sprintf("%s/target.json?pref_name__icontains=%s&limit=%d", ...
        BASE_URL, encodedQuery, maxRows);

    if filterTargetType
        encodedType = strrep(string(urlencode(char(options.TargetType))), "+", "%20");
        url = sprintf("%s&target_type=%s", url, encodedType);
    end

    logInfo("searchChemblTarget: query='%s', TargetType='%s'", ...
        query, options.TargetType);
    logDebug("searchChemblTarget: URL=%s", url);

    % --- Call ChEMBL REST API via webread ---
    opts = weboptions("Timeout", 30, "ContentType", "json");
    try
        data = webread(url, opts);
    catch ME
        msg = ME.message;
        if contains(lower(msg), "404") || contains(lower(msg), "not found")
            error("emk:db:searchChemblTarget:notFound", ...
                "No targets found in ChEMBL for query: %s", query);
        end
        error("emk:db:searchChemblTarget:networkError", ...
            "ChEMBL request failed: %s", msg);
    end

    % --- Parse response ---
    if ~isfield(data, "targets")
        error("emk:db:searchChemblTarget:networkError", ...
            "Unexpected ChEMBL response (targets field missing)");
    end

    raw = data.targets;
    if isstruct(raw)
        targets = raw;
    elseif iscell(raw)
        if isempty(raw)
            error("emk:db:searchChemblTarget:notFound", ...
                "No targets found in ChEMBL for query: %s", query);
        end
        targets = [raw{:}];
    elseif isempty(raw)
        error("emk:db:searchChemblTarget:notFound", ...
            "No targets found in ChEMBL for query: %s", query);
    else
        error("emk:db:searchChemblTarget:networkError", ...
            "Unexpected targets field type: %s", class(raw));
    end

    nRows = numel(targets);
    targetIds   = strings(nRows, 1);
    prefNames   = strings(nRows, 1);
    organisms   = strings(nRows, 1);
    targetTypes = strings(nRows, 1);

    for i = 1:nRows
        t = targets(i);
        targetIds(i)   = safeStr(t, "target_chembl_id", "");
        prefNames(i)   = safeStr(t, "pref_name",        "");
        organisms(i)   = safeStr(t, "organism",         "");
        targetTypes(i) = safeStr(t, "target_type",      "");
    end

    % --- Client-side TargetType filtering ---
    if filterTargetType
        mask = strcmpi(targetTypes, options.TargetType);
        targetIds   = targetIds(mask);
        prefNames   = prefNames(mask);
        organisms   = organisms(mask);
        targetTypes = targetTypes(mask);
        if isempty(targetIds)
            error("emk:db:searchChemblTarget:notFound", ...
                "No targets of type '%s' found in ChEMBL for query: %s", ...
                options.TargetType, query);
        end
        % Trim to requested MaxRows
        if numel(targetIds) > maxRows
            targetIds   = targetIds(1:maxRows);
            prefNames   = prefNames(1:maxRows);
            organisms   = organisms(1:maxRows);
            targetTypes = targetTypes(1:maxRows);
        end
    end

    result = table(targetIds, prefNames, organisms, targetTypes, ...
        VariableNames=["TargetChEMBLID","PreferredName","Organism","TargetType"]);

    logInfo("searchChemblTarget: found %d target(s)", height(result));
end

% -----------------------------------------------------------------------
function val = safeStr(s, fname, default)
    if isfield(s, fname) && ~isempty(s.(fname)) && ...
       ~isstruct(s.(fname)) && ~iscell(s.(fname))
        val = string(s.(fname));
    else
        val = string(default);
    end
end
