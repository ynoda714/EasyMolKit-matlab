function result = getChemblActivity(targetId, options)
% getChemblActivity  Retrieve bioactivity data from ChEMBL for a protein target.
%
%   result = emk.db.getChemblActivity(targetId)
%   result = emk.db.getChemblActivity(targetId, ActivityType="IC50")
%   result = emk.db.getChemblActivity(targetId, MaxRows=100)
%   result = emk.db.getChemblActivity(targetId, MinActivity_nM=100)
%
%   Queries the ChEMBL REST API activity endpoint for bioactivity measurements
%   against the specified ChEMBL protein target.  Only rows with a valid
%   canonical SMILES and a numeric standard_value in nM are returned.
%   Python is not required.
%
%   ChEMBL REST endpoint used:
%     https://www.ebi.ac.uk/chembl/api/data/activity.json
%       ?target_chembl_id={targetId}&standard_type={ActivityType}
%       &standard_relation=%3D&limit={MaxRows}
%
%   Arguments:
%     targetId       - string | char, ChEMBL target identifier (e.g. "CHEMBL203").
%     ActivityType   - string, bioactivity measurement type. Default: "IC50".
%                      Common values: "IC50", "Ki", "Kd", "EC50".
%     MaxRows        - positive integer, max activity records to fetch from the
%                      ChEMBL API before filtering. Default: 50.
%     MinActivity_nM - double, upper potency threshold in nM.  Rows with
%                      Value_nM > MinActivity_nM are discarded.
%                      Default: Inf (no threshold; return all valid rows).
%
%   Returns:
%     result  - MATLAB table with columns:
%       MoleculeChEMBLID  string  Compound ChEMBL identifier (e.g. "CHEMBL553")
%       Name              string  Preferred molecule name (may be "")
%       SMILES            string  Canonical SMILES
%       ActivityType      string  Measurement type (e.g. "IC50")
%       Value_nM          double  Activity value in nM
%
%   Errors:
%     emk:db:getChemblActivity:invalidInput   - targetId is not a non-empty string
%     emk:db:getChemblActivity:invalidOptions - invalid MaxRows value
%     emk:db:getChemblActivity:notFound       - no valid activities found for target
%     emk:db:getChemblActivity:networkError   - webread failure or bad response
%
%   Example:
%     tbl = emk.db.getChemblActivity("CHEMBL203");
%     tbl = emk.db.getChemblActivity("CHEMBL203", ActivityType="Ki", MaxRows=100);
%     tbl = emk.db.getChemblActivity("CHEMBL203", MinActivity_nM=100);
%
%   See also: emk.db.searchChemblTarget, emk.db.searchChembl,
%             emk.io.writeSdf, emk.fingerprint.morgan

    arguments
        targetId
        options.ActivityType   (1,1) string = "IC50"
        options.MaxRows        (1,1) double = 50
        options.MinActivity_nM (1,1) double = Inf
    end

    % --- Input validation: targetId ---
    if ~(ischar(targetId) || isStringScalar(targetId))
        error("emk:db:getChemblActivity:invalidInput", ...
            "targetId must be a string or char, got: %s", class(targetId));
    end
    targetId = string(targetId);
    if strlength(strtrim(targetId)) == 0
        error("emk:db:getChemblActivity:invalidInput", ...
            "targetId must not be empty");
    end

    % --- Input validation: MaxRows ---
    if ~(isscalar(options.MaxRows) && options.MaxRows > 0 && ...
         floor(options.MaxRows) == options.MaxRows && isfinite(options.MaxRows))
        error("emk:db:getChemblActivity:invalidOptions", ...
            "MaxRows must be a positive finite integer, got: %g", options.MaxRows);
    end
    maxRows = floor(options.MaxRows);

    % --- Build ChEMBL REST URL ---
    % CHEMBL IDs are alphanumeric; no URL-encoding required.
    % standard_relation=%3D selects only exact-equality measurements ("=").
    BASE_URL    = "https://www.ebi.ac.uk/chembl/api/data";
    encodedType = string(urlencode(char(options.ActivityType)));
    url = sprintf( ...
        "%s/activity.json?target_chembl_id=%s&standard_type=%s" + ...
        "&standard_relation=%%3D&limit=%d", ...
        BASE_URL, targetId, encodedType, maxRows);

    logInfo("getChemblActivity: target=%s, type=%s, maxRows=%d", ...
        targetId, options.ActivityType, maxRows);
    logDebug("getChemblActivity: URL=%s", url);

    % --- Call ChEMBL REST API via webread ---
    opts = weboptions("Timeout", 20, "ContentType", "json");
    try
        data = webread(url, opts);
    catch ME
        msg = ME.message;
        if contains(lower(msg), "404") || contains(lower(msg), "not found")
            error("emk:db:getChemblActivity:notFound", ...
                "No activities found in ChEMBL for target: %s", targetId);
        end
        error("emk:db:getChemblActivity:networkError", ...
            "ChEMBL request failed: %s", msg);
    end

    % --- Parse response ---
    if ~isfield(data, "activities")
        error("emk:db:getChemblActivity:networkError", ...
            "Unexpected ChEMBL response (activities field missing)");
    end

    raw = data.activities;
    if isstruct(raw)
        activities = raw;
    elseif iscell(raw)
        if isempty(raw)
            error("emk:db:getChemblActivity:notFound", ...
                "No activities found in ChEMBL for target: %s", targetId);
        end
        activities = [raw{:}];
    elseif isempty(raw)
        error("emk:db:getChemblActivity:notFound", ...
            "No activities found in ChEMBL for target: %s", targetId);
    else
        error("emk:db:getChemblActivity:networkError", ...
            "Unexpected activities field type: %s", class(raw));
    end

    nRaw = numel(activities);
    logInfo("getChemblActivity: received %d raw activity records", nRaw);

    % --- Extract fields ---
    molIds     = strings(nRaw, 1);
    names      = strings(nRaw, 1);
    smilesList = strings(nRaw, 1);
    actTypes   = strings(nRaw, 1);
    values     = nan(nRaw, 1);
    isNanoMolar = true(nRaw, 1);

    for i = 1:nRaw
        a = activities(i);
        molIds(i)     = safeStr(a, "molecule_chembl_id", "");
        names(i)      = safeStr(a, "molecule_pref_name", "");
        smilesList(i) = safeStr(a, "canonical_smiles",   "");
        actTypes(i)   = safeStr(a, "standard_type",      "");

        % standard_value is returned as a numeric string by ChEMBL JSON.
        if isfield(a, "standard_value") && ~isempty(a.standard_value) && ...
           ~isstruct(a.standard_value) && ~iscell(a.standard_value)
            if isnumeric(a.standard_value)
                values(i) = double(a.standard_value);
            else
                values(i) = str2double(string(a.standard_value));
            end
        end

        % Keep only nM measurements; discard uM, mg/mL, etc.
        units = safeStr(a, "standard_units", "nM");
        if ~strcmpi(units, "nM")
            isNanoMolar(i) = false;
        end
    end

    % --- Filter: require valid SMILES + non-NaN value + positive nM units ---
    keepMask = strlength(smilesList) > 0 & ~isnan(values) & isNanoMolar & values > 0;

    % --- Apply optional potency threshold ---
    if isfinite(options.MinActivity_nM)
        keepMask = keepMask & (values <= options.MinActivity_nM);
    end

    if ~any(keepMask)
        error("emk:db:getChemblActivity:notFound", ...
            "No valid nM activities with SMILES found for target: %s " + ...
            "(ActivityType=%s, MinActivity_nM=%.0f)", ...
            targetId, options.ActivityType, options.MinActivity_nM);
    end

    result = table( ...
        molIds(keepMask), names(keepMask), smilesList(keepMask), ...
        actTypes(keepMask), values(keepMask), ...
        'VariableNames', ...
        {'MoleculeChEMBLID','Name','SMILES','ActivityType','Value_nM'});

    logInfo("getChemblActivity: returning %d records (of %d raw)", ...
        sum(keepMask), nRaw);
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
