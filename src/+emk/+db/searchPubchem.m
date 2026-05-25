function result = searchPubchem(query, options)
% searchPubchem  Search PubChem for compound information via PUG REST API.
%
%   result = emk.db.searchPubchem(query)
%   result = emk.db.searchPubchem(query, Type="name")
%   result = emk.db.searchPubchem(query, Type="smiles")
%   result = emk.db.searchPubchem(query, Type="cid")
%   result = emk.db.searchPubchem(query, Type="inchikey")
%
%   Queries the PubChem PUG REST API using MATLAB webread only.
%   Python is not required.  Returns compound properties as a MATLAB table.
%
%   PUG REST endpoint used:
%     https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{type}/{query}/
%       property/IUPACName,MolecularFormula,MolecularWeight,IsomericSMILES/JSON
%
%   Arguments:
%     query   - string | char, the compound identifier.  Examples:
%               "aspirin" (Type="name"), "CCO" (Type="smiles"),
%               "702" (Type="cid"), "LFQSCWFLJHTTHZ-UHFFFAOYSA-N" (Type="inchikey")
%     Type    - string, query namespace.  One of:
%               "name" (default), "smiles", "cid", "inchikey"
%
%   Returns:
%     result  - MATLAB table with the following columns:
%       CID              - double - PubChem Compound Identifier
%       IUPACName        - string - IUPAC systematic name
%       MolecularFormula - string - molecular formula (e.g. "C2H6O")
%       MolecularWeight  - double - average molecular weight (g/mol)
%       IsomericSMILES   - string - isomeric SMILES string
%
%   Errors:
%     emk:db:searchPubchem:invalidInput  - query is not a non-empty string/char
%     emk:db:searchPubchem:invalidType   - Type is not a recognised namespace
%     emk:db:searchPubchem:notFound      - no compound found for the query
%     emk:db:searchPubchem:networkError  - webread failure or unexpected response
%
%   Example:
%     tbl = emk.db.searchPubchem("ethanol");
%     tbl = emk.db.searchPubchem("CCO", Type="smiles");
%     tbl = emk.db.searchPubchem("702",  Type="cid");
%
%   See also: emk.mol.fromSmiles, emk.descriptor.calculate,
%             emk.fingerprint.morgan

    arguments
        query
        options.Type (1,1) string = "name"
    end

    % --- Input validation: query must be a non-empty string/char ---
    if ~(ischar(query) || isStringScalar(query))
        error("emk:db:searchPubchem:invalidInput", ...
            "query must be a string or char, got: %s", class(query));
    end
    query = string(query);
    if strlength(strtrim(query)) == 0
        error("emk:db:searchPubchem:invalidInput", ...
            "query must not be empty or whitespace-only");
    end

    % --- Input validation: Type ---
    VALID_TYPES = ["name", "smiles", "cid", "inchikey"];
    qtype = lower(options.Type);
    if ~ismember(qtype, VALID_TYPES)
        error("emk:db:searchPubchem:invalidType", ...
            "Type must be one of: %s. Got: %s", ...
            strjoin(VALID_TYPES, ", "), options.Type);
    end

    % --- Build PUG REST URL ---
    % SMILES and other inputs may contain characters that need URL-encoding.
    % urlencode() is available in MATLAB R2022a+.
    % Request both IsomericSMILES and CanonicalSMILES so that compounds
    % without defined stereochemistry still return a usable SMILES string.
    BASE_URL  = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound";
    PROP_LIST = "IUPACName,MolecularFormula,MolecularWeight,IsomericSMILES,CanonicalSMILES";
    % urlencode encodes spaces as '+', but PubChem PUG REST URL path
    % requires '%20'.  Replace '+' with '%20' after encoding.
    encodedQuery = strrep(string(urlencode(char(query))), "+", "%20");
    url = sprintf("%s/%s/%s/property/%s/JSON", ...
        BASE_URL, qtype, encodedQuery, PROP_LIST);

    logInfo("searchPubchem: querying PubChem (Type=%s, query=%s)", qtype, query);
    logDebug("searchPubchem: URL=%s", url);

    % --- Call PUG REST API via webread ---
    opts = weboptions("Timeout", 15, "ContentType", "json");
    try
        data = webread(url, opts);
    catch ME
        msg = ME.message;
        % Detect HTTP 404 / not-found responses from PubChem.
        % PubChem 404 body: {"Fault":{"Code":"PUGREST.NotFound",...}}
        % MATLAB raises the HTTP status in the error message or identifier.
        if contains(lower(msg), "404") || ...
                contains(lower(msg), "not found") || ...
                contains(lower(msg), "pugrest.notfound")
            error("emk:db:searchPubchem:notFound", ...
                "No compound found in PubChem for query: %s (Type=%s)", ...
                query, qtype);
        end
        error("emk:db:searchPubchem:networkError", ...
            "PubChem request failed: %s", msg);
    end

    % --- Parse JSON response ---
    % webread with ContentType="json" returns a MATLAB struct via jsondecode.
    % PropertyTable.Properties is a struct array (one element per compound).
    try
        rawProps = data.PropertyTable.Properties;
    catch
        error("emk:db:searchPubchem:networkError", ...
            "Unexpected PubChem response format (PropertyTable missing)");
    end

    % Normalize rawProps to a cell array of scalar structs for uniform access.
    % jsondecode returns a struct array when all JSON objects share identical
    % field sets, or a cell array when field sets differ (heterogeneous).
    % Avoid [rawProps{:}] struct concatenation which fails when fields differ.
    if isstruct(rawProps)
        props_cell = num2cell(rawProps);  % struct array -> cell of scalar structs
    elseif iscell(rawProps)
        props_cell = rawProps;
    else
        error("emk:db:searchPubchem:networkError", ...
            "Unexpected Properties type in PubChem response: %s", class(rawProps));
    end

    nRows = numel(props_cell);

    % --- Build output table ---
    cids     = zeros(nRows, 1);
    names    = strings(nRows, 1);
    formulas = strings(nRows, 1);
    mws      = zeros(nRows, 1);
    smiles   = strings(nRows, 1);

    for i = 1:nRows
        p = props_cell{i};
        % CID: PubChem normally returns a JSON integer, but guard against
        % the rare case where jsondecode gives a char/string (same pattern
        % as MolecularWeight below).
        if isnumeric(p.CID)
            cids(i) = double(p.CID);
        else
            cids(i) = str2double(string(p.CID));
        end
        % IUPACName / MolecularFormula may be absent for some records
        if isfield(p, 'IUPACName')
            names(i) = string(p.IUPACName);
        end
        if isfield(p, 'MolecularFormula')
            formulas(i) = string(p.MolecularFormula);
        end
        % MolecularWeight may arrive as string ("46.1") or number; handle both.
        if isfield(p, 'MolecularWeight')
            mwVal = p.MolecularWeight;
            if isnumeric(mwVal)
                mws(i) = double(mwVal);
            else
                mws(i) = str2double(string(mwVal));
            end
        end
        % PubChem PUG REST maps the requested property names to different
        % JSON keys in the response.  Known mappings (verified 2026-04-19):
        %   Requested "IsomericSMILES"  -> returned as "SMILES"
        %   Requested "CanonicalSMILES" -> returned as "ConnectivitySMILES"
        % We probe all four field names in order of preference so that the
        % code remains correct if PubChem corrects the response keys.
        smiVal = "";
        for smiField = {'SMILES', 'IsomericSMILES', 'ConnectivitySMILES', 'CanonicalSMILES'}
            fn = smiField{1};
            if isfield(p, fn) && ~isempty(p.(fn))
                candidate = strtrim(string(p.(fn)));
                if strlength(candidate) > 0
                    smiVal = candidate;
                    break;
                end
            end
        end
        smiles(i) = smiVal;
    end

    % Use char-vector 'VariableNames' + cell array to ensure correct
    % Name-Value pair parsing across all MATLAB versions.
    result = table(cids, names, formulas, mws, smiles, ...
        'VariableNames', ...
        {'CID', 'IUPACName', 'MolecularFormula', 'MolecularWeight', 'IsomericSMILES'});

    logInfo("searchPubchem: found %d compound(s)", nRows);
end
