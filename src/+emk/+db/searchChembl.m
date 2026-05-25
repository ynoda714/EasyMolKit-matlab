function result = searchChembl(query, options)
% searchChembl  Search ChEMBL for molecule information via REST API.
%
%   result = emk.db.searchChembl(query)
%   result = emk.db.searchChembl(query, Type="name")
%   result = emk.db.searchChembl(query, Type="smiles")
%   result = emk.db.searchChembl(query, Type="chemblid")
%   result = emk.db.searchChembl(query, Type="inchikey")
%
%   Queries the ChEMBL REST API using MATLAB webread only.
%   Python is not required.  Returns molecule properties as a MATLAB table.
%
%   ChEMBL REST endpoints used:
%     List:   https://www.ebi.ac.uk/chembl/api/data/molecule.json?{filter}&limit=25
%     Single: https://www.ebi.ac.uk/chembl/api/data/molecule/{chemblid}.json
%
%   Arguments:
%     query   - string | char, the compound identifier.  Examples:
%               "aspirin" (Type="name"),
%               "CC(=O)Oc1ccccc1C(=O)O" (Type="smiles"),
%               "CHEMBL192" (Type="chemblid"),
%               "BSYNRYMUTXBXSQ-UHFFFAOYSA-N" (Type="inchikey")
%     Type    - string, query namespace.  One of:
%               "name" (default), "smiles", "chemblid", "inchikey"
%
%   Returns:
%     result  - MATLAB table with the following columns:
%       ChEMBLID         - string - ChEMBL compound identifier (e.g. "CHEMBL192")
%       Name             - string - preferred molecule name (may be empty)
%       MolecularWeight  - double - full molecular weight (g/mol)
%       ALogP            - double - calculated ALogP
%       HBondDonors      - double - H-bond donor count
%       HBondAcceptors   - double - H-bond acceptor count
%       SMILES           - string - canonical SMILES (may be empty for no-structure)
%       InChIKey         - string - standard InChIKey (may be empty)
%
%   Errors:
%     emk:db:searchChembl:invalidInput  - query is not a non-empty string/char
%     emk:db:searchChembl:invalidType   - Type is not a recognised namespace
%     emk:db:searchChembl:notFound      - no molecule found for the query
%     emk:db:searchChembl:networkError  - webread failure or unexpected response
%
%   Example:
%     tbl = emk.db.searchChembl("aspirin");
%     tbl = emk.db.searchChembl("CC(=O)Oc1ccccc1C(=O)O", Type="smiles");
%     tbl = emk.db.searchChembl("CHEMBL192",             Type="chemblid");
%     tbl = emk.db.searchChembl("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", Type="inchikey");
%
%   See also: emk.db.searchPubchem, emk.mol.fromSmiles,
%             emk.descriptor.calculate, emk.fingerprint.morgan

    arguments
        query
        options.Type (1,1) string = "name"
    end

    % --- Input validation: query must be a non-empty string/char ---
    if ~(ischar(query) || isStringScalar(query))
        error("emk:db:searchChembl:invalidInput", ...
            "query must be a string or char, got: %s", class(query));
    end
    query = string(query);
    if strlength(strtrim(query)) == 0
        error("emk:db:searchChembl:invalidInput", ...
            "query must not be empty or whitespace-only");
    end

    % --- Input validation: Type ---
    VALID_TYPES = ["name", "smiles", "chemblid", "inchikey"];
    qtype = lower(options.Type);
    if ~ismember(qtype, VALID_TYPES)
        error("emk:db:searchChembl:invalidType", ...
            "Type must be one of: %s. Got: %s", ...
            strjoin(VALID_TYPES, ", "), options.Type);
    end

    % --- Build ChEMBL REST URL ---
    BASE_URL = "https://www.ebi.ac.uk/chembl/api/data";
    encodedQuery = string(urlencode(char(query)));

    if qtype == "chemblid"
        % Single-compound lookup by ChEMBL ID.
        % The ChEMBL ID itself is URL-safe (alphanumeric only).
        url = sprintf("%s/molecule/%s.json", BASE_URL, encodedQuery);
        isSingleLookup = true;
    else
        % List endpoint with ORM-style filter parameters.
        switch qtype
            case "name"
                filterParam = sprintf("pref_name__iexact=%s", encodedQuery);
            case "smiles"
                filterParam = sprintf( ...
                    "molecule_structures__canonical_smiles=%s", encodedQuery);
            case "inchikey"
                filterParam = sprintf( ...
                    "molecule_structures__standard_inchi_key=%s", encodedQuery);
            otherwise
                % Should not reach here; covered by VALID_TYPES check above.
                filterParam = "";
        end
        url = sprintf("%s/molecule.json?%s&limit=25", BASE_URL, filterParam);
        isSingleLookup = false;
    end

    logInfo("searchChembl: querying ChEMBL (Type=%s, query=%s)", qtype, query);
    logDebug("searchChembl: URL=%s", url);

    % --- Call ChEMBL REST API via webread ---
    opts = weboptions("Timeout", 15, "ContentType", "json");
    try
        data = webread(url, opts);
    catch ME
        msg = ME.message;
        % Detect HTTP 404 / not-found responses.
        if contains(lower(msg), "404") || contains(lower(msg), "not found")
            error("emk:db:searchChembl:notFound", ...
                "No molecule found in ChEMBL for query: %s (Type=%s)", ...
                query, qtype);
        end
        error("emk:db:searchChembl:networkError", ...
            "ChEMBL request failed: %s", msg);
    end

    % --- Normalise to struct array of molecule records ---
    if isSingleLookup
        % webread returns a single molecule struct directly.
        molecules = data;
        nRows = 1;
    else
        % List response: data.molecules is a struct array, cell array, or empty.
        if ~isfield(data, 'molecules')
            error("emk:db:searchChembl:networkError", ...
                "Unexpected ChEMBL response format (molecules field missing)");
        end
        raw = data.molecules;
        if isstruct(raw)
            molecules = raw;
        elseif iscell(raw)
            if isempty(raw)
                error("emk:db:searchChembl:notFound", ...
                    "No molecule found in ChEMBL for query: %s (Type=%s)", ...
                    query, qtype);
            end
            molecules = [raw{:}];
        elseif isempty(raw)
            % jsondecode may return [] for an empty JSON array.
            error("emk:db:searchChembl:notFound", ...
                "No molecule found in ChEMBL for query: %s (Type=%s)", ...
                query, qtype);
        else
            error("emk:db:searchChembl:networkError", ...
                "Unexpected molecules field type in ChEMBL response: %s", ...
                class(raw));
        end
        nRows = numel(molecules);
        if nRows == 0
            error("emk:db:searchChembl:notFound", ...
                "No molecule found in ChEMBL for query: %s (Type=%s)", ...
                query, qtype);
        end
    end

    % --- Build output table ---
    chemblIds = strings(nRows, 1);
    names     = strings(nRows, 1);
    mws       = zeros(nRows, 1);
    alogps    = zeros(nRows, 1);
    hbds      = zeros(nRows, 1);
    hbas      = zeros(nRows, 1);
    smiles    = strings(nRows, 1);
    inchikeys = strings(nRows, 1);

    for i = 1:nRows
        m = molecules(i);

        % --- ChEMBL ID ---
        if isfield(m, 'molecule_chembl_id') && ~isempty(m.molecule_chembl_id)
            chemblIds(i) = string(m.molecule_chembl_id);
        end

        % --- Preferred name ---
        if isfield(m, 'pref_name') && ~isempty(m.pref_name)
            names(i) = string(m.pref_name);
        end

        % --- Molecule properties sub-struct ---
        % ChEMBL returns molecule_properties as null for some entries (no
        % computed properties available).  jsondecode converts null -> [].
        if isfield(m, 'molecule_properties') && isstruct(m.molecule_properties)
            mp = m.molecule_properties;

            mws(i)    = extractNumericField(mp, 'full_mwt');
            alogps(i) = extractNumericField(mp, 'alogp');
            hbds(i)   = extractNumericField(mp, 'hbd');
            hbas(i)   = extractNumericField(mp, 'hba');
        end

        % --- Molecule structures sub-struct ---
        if isfield(m, 'molecule_structures') && isstruct(m.molecule_structures)
            ms = m.molecule_structures;
            if isfield(ms, 'canonical_smiles') && ~isempty(ms.canonical_smiles)
                smiles(i) = strtrim(string(ms.canonical_smiles));
            end
            if isfield(ms, 'standard_inchi_key') && ~isempty(ms.standard_inchi_key)
                inchikeys(i) = string(ms.standard_inchi_key);
            end
        end
    end

    result = table(chemblIds, names, mws, alogps, hbds, hbas, smiles, inchikeys, ...
        'VariableNames', ...
        {'ChEMBLID', 'Name', 'MolecularWeight', 'ALogP', ...
         'HBondDonors', 'HBondAcceptors', 'SMILES', 'InChIKey'});

    logInfo("searchChembl: found %d molecule(s)", nRows);
end

% ======================================================================
% Local helper
% ======================================================================

function val = extractNumericField(s, fieldName)
% Return double value of a struct field; handles both numeric and string
% JSON types (ChEMBL returns numeric properties as JSON strings).
    val = 0;
    if isfield(s, fieldName) && ~isempty(s.(fieldName))
        raw = s.(fieldName);
        if isnumeric(raw)
            val = double(raw);
        else
            val = str2double(string(raw));
            if isnan(val)
                val = 0;
            end
        end
    end
end
