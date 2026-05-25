function result = rankBy(queryFp, dbFps, N, options)
% rankBy  Rank database fingerprints by similarity to a query fingerprint.
%
%   result = emk.similarity.rankBy(queryFp, dbFps)
%   result = emk.similarity.rankBy(queryFp, dbFps, N)
%   result = emk.similarity.rankBy(queryFp, dbFps, N, Metric="tanimoto")
%   result = emk.similarity.rankBy(queryFp, dbFps, N, Metric="dice")
%
%   Returns the top-N fingerprints from dbFps ranked by decreasing
%   similarity to queryFp.  All M database fingerprints are scored in a
%   single Python call via BulkTanimotoSimilarity / BulkDiceSimilarity,
%   minimising IPC round-trips (ADR-002 rev.3 principle).
%
%   Arguments:
%     queryFp  - py.rdkit.DataStructs.cDataStructs.ExplicitBitVect
%                Query fingerprint (Python reference from emk.fingerprint.*)
%     dbFps    - 1xM cell array of py.rdkit.DataStructs.ExplicitBitVect
%                Database fingerprints.  All entries must share the same
%                bit length as queryFp.
%     N        - positive integer scalar, number of top results to return.
%                Default: Inf (return full ranking of all M entries).
%                If N > M, all M results are returned.
%     Metric   - string, "tanimoto" (default) or "dice".
%
%   Returns:
%     result   - struct with fields:
%       .Indices  - double (1xK) - 1-based indices into dbFps, sorted by
%                  descending similarity score.  K = min(N, M).
%       .Scores   - double (1xK) - similarity scores in [0, 1], descending.
%       .Metric   - string - the metric used for scoring.
%
%   Errors:
%     emk:similarity:rankBy:invalidQueryFp  - queryFp is not a Python object
%     emk:similarity:rankBy:invalidDbFps    - dbFps is not a non-empty cell
%                                             array of Python objects
%     emk:similarity:rankBy:invalidN        - N is not a positive integer
%                                             scalar (or Inf)
%     emk:similarity:rankBy:invalidMetric   - Metric is not "tanimoto" or "dice"
%     emk:similarity:rankBy:rdkitError      - Python exception (e.g. bit-length
%                                             mismatch between queryFp and dbFps)
%
%   Example:
%     query = emk.fingerprint.morgan(emk.mol.fromSmiles("CCO"));
%     smis  = {"CC(=O)Oc1ccccc1C(=O)O", "c1ccccc1", "CCCO"};
%     db    = cellfun(@(s) emk.fingerprint.morgan(emk.mol.fromSmiles(s)), ...
%                 smis, "UniformOutput", false);
%     res   = emk.similarity.rankBy(query, db, 2);
%     % res.Indices: top-2 indices, res.Scores: corresponding scores
%
%   See also: emk.similarity.tanimoto, emk.similarity.dice,
%             emk.fingerprint.morgan, emk.fingerprint.maccs

    arguments
        queryFp
        dbFps
        N = Inf
        options.Metric (1,1) string = "tanimoto"
    end

    % --- Validate Metric ---
    metric = options.Metric;
    if ~ismember(metric, ["tanimoto", "dice"])
        error("emk:similarity:rankBy:invalidMetric", ...
            "Metric must be ""tanimoto"" or ""dice"", got: %s", metric);
    end

    % --- Validate queryFp: must be Python object (fires before any RDKit call) ---
    if ~startsWith(class(queryFp), "py.")
        error("emk:similarity:rankBy:invalidQueryFp", ...
            "queryFp must be a Python fingerprint object (ExplicitBitVect), got: %s", ...
            class(queryFp));
    end

    % --- Validate dbFps: must be non-empty cell array of Python objects ---
    if ~iscell(dbFps)
        error("emk:similarity:rankBy:invalidDbFps", ...
            "dbFps must be a cell array of Python fingerprint objects, got: %s", ...
            class(dbFps));
    end
    if isempty(dbFps)
        error("emk:similarity:rankBy:invalidDbFps", ...
            "dbFps must be a non-empty cell array");
    end
    for i = 1:numel(dbFps)
        if ~startsWith(class(dbFps{i}), "py.")
            error("emk:similarity:rankBy:invalidDbFps", ...
                "dbFps{%d} must be a Python fingerprint object, got: %s", ...
                i, class(dbFps{i}));
        end
    end

    % --- Validate N ---
    if ~isscalar(N) || ~isnumeric(N) || N <= 0 || ...
            (isfinite(N) && floor(N) ~= N)
        error("emk:similarity:rankBy:invalidN", ...
            "N must be a positive integer scalar or Inf, got: %s", ...
            mat2str(N));
    end

    M = numel(dbFps);
    if isinf(N)
        K = M;
    else
        K = min(double(N), M);
    end

    logDebug("rankBy: M=%d, K=%d, Metric=%s", M, K, metric);

    % --- Compute all similarity scores in a single IPC round trip ---
    % BulkTanimotoSimilarity / BulkDiceSimilarity evaluates one queryFp
    % against a Python list of M database fingerprints at once.
    % Use emk.util.rdkitModule() to get DataStructs via importlib.
    mods = emk.util.rdkitModule();
    try
        % Build Python list from MATLAB cell of Python objects.
        % py.list(cell) only works for 1-N numeric/logical cells; for cells
        % containing Python objects we must append element-by-element.
        pyList = py.list();
        for k = 1:M
            pyList.append(dbFps{k});
        end
        if metric == "tanimoto"
            pyScores = mods.DataStructs.BulkTanimotoSimilarity(queryFp, pyList);
        else
            pyScores = mods.DataStructs.BulkDiceSimilarity(queryFp, pyList);
        end
        % Convert Python list of floats to MATLAB double row vector.
        % py.array.array("d", ...) produces a typed Python array that
        % double() can consume directly.
        scores = double(py.array.array("d", pyScores));
    catch ME
        error("emk:similarity:rankBy:rdkitError", ...
            "RDKit bulk %s similarity raised an exception: %s", ...
            metric, ME.message);
    end

    % --- Sort descending and select top K ---
    [sortedScores, sortedIdx] = sort(scores, "descend");

    result.Indices = sortedIdx(1:K);
    result.Scores  = sortedScores(1:K);
    result.Metric  = metric;

    logDebug("rankBy: top score=%.4f at index=%d", result.Scores(1), result.Indices(1));
end
