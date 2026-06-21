function clusters = butina(fps, varargin)
% butina  Cluster fingerprints using the Butina algorithm.
%
%   clusters = emk.cluster.butina(fps)
%   clusters = emk.cluster.butina(fps, Threshold=0.2)
%   clusters = emk.cluster.butina(fps, Threshold=0.2, Metric="tanimoto")
%
%   Performs Butina (Taylor-Butina) sphere-exclusion clustering on a set
%   of molecular fingerprints.  The algorithm assigns each fingerprint to
%   exactly one cluster; cluster centroids are the molecules with the most
%   neighbours within the given distance threshold.
%
%   Key properties of Butina clustering:
%   - Non-hierarchical; single-pass algorithm (linear in N)
%   - Each molecule belongs to exactly one cluster
%   - Cluster size and count depend heavily on the Threshold parameter
%   - Default Threshold=0.2 corresponds to 80% Tanimoto similarity
%     (molecules within the same cluster are >= 80% similar)
%
%   Arguments:
%     fps        - 1xM cell array of py.rdkit.DataStructs.ExplicitBitVect
%                  All fingerprints must have the same bit length.
%     Threshold  - double in (0, 1] (optional, default 0.2)
%                  Tanimoto DISTANCE threshold for cluster membership.
%                  Lower value -> fewer, larger clusters (stricter similarity).
%                  Higher value -> more, smaller clusters (looser similarity).
%     Metric     - string (optional, default "tanimoto")
%                  Similarity metric for distance calculation.
%                  Only "tanimoto" is supported by RDKit's Butina module.
%
%   Returns:
%     clusters - 1xC cell array (C = number of clusters)
%                Each element is a double(1,K) array of 1-based indices
%                into fps for the K molecules assigned to that cluster.
%                clusters{1} always contains the largest cluster (most members).
%                The first index in each clusters{i} is the centroid.
%
%   Errors:
%     emk:cluster:butina:invalidInput   - fps is not a non-empty cell array
%                                         of Python objects
%     emk:cluster:butina:invalidThreshold - Threshold not in (0,1]
%     emk:cluster:butina:invalidMetric  - Metric is not "tanimoto"
%     emk:cluster:butina:rdkitError     - unexpected Python exception
%
%   Example:
%     smiles = {"CCO", "CCCO", "CCCCO", "c1ccccc1", "c1ccc2ccccc2c1"};
%     mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
%     fps    = cellfun(@emk.fingerprint.morgan, mols, "UniformOutput", false);
%     clusters = emk.cluster.butina(fps, "Threshold", 0.4);
%     fprintf("Number of clusters: %d\n", numel(clusters));
%     fprintf("Largest cluster size: %d\n", numel(clusters{1}));
%
%   References:
%     Butina, D. (1999). Unsupervised Data Base Clustering Based on
%       Daylight's Fingerprint and Tanimoto Similarity: A Fast and
%       Automated Way to Cluster Small and Large Data Sets. J. Chem.
%       Inf. Comput. Sci. 39(4):747-750. DOI: 10.1021/ci9803381
%     Taylor, R. (1995). Simulation Analysis of Experimental Design
%       Strategies for Screening Random Compounds as Potential New Drugs
%       and Agrochemicals. J. Chem. Inf. Comput. Sci. 35(1):59-67.
%       DOI: 10.1021/ci00023a009
%
%   See also: emk.diversity.pick, emk.similarity.tanimoto,
%             emk.similarity.matrix, emk.fingerprint.morgan

    % --- Parse optional name-value arguments ---
    threshold = 0.2;
    metric    = "tanimoto";
    if nargin > 1
        if mod(numel(varargin), 2) ~= 0
            error("emk:cluster:butina:invalidInput", ...
                "Optional arguments must be name-value pairs.");
        end
        for k = 1:2:numel(varargin)
            argName = string(varargin{k});
            switch lower(argName)
                case "threshold"
                    threshold = varargin{k+1};
                case "metric"
                    metric = lower(string(varargin{k+1}));
                otherwise
                    error("emk:cluster:butina:invalidInput", ...
                        "Unknown option: '%s'. Supported: Threshold, Metric.", ...
                        varargin{k});
            end
        end
    end

    % --- Validate Threshold ---
    if ~(isnumeric(threshold) && isscalar(threshold) && isfinite(threshold) && ...
         threshold > 0 && threshold <= 1)
        error("emk:cluster:butina:invalidThreshold", ...
            "Threshold must be a finite scalar in (0, 1], got: %s", mat2str(threshold));
    end

    % --- Validate Metric ---
    if metric ~= "tanimoto"
        error("emk:cluster:butina:invalidMetric", ...
            "Metric must be 'tanimoto' (only RDKit Butina metric supported), got: %s", metric);
    end

    % --- Validate fps ---
    if ~iscell(fps) || isempty(fps)
        error("emk:cluster:butina:invalidInput", ...
            "fps must be a non-empty cell array, got: %s", class(fps));
    end
    for i = 1:numel(fps)
        if ~startsWith(class(fps{i}), "py.")
            error("emk:cluster:butina:invalidInput", ...
                "fps{%d} must be a Python fingerprint object, got: %s", ...
                i, class(fps{i}));
        end
    end

    M = numel(fps);
    logDebug("butina: clustering %d fingerprints, threshold=%.3f", M, threshold);

    % --- Call RDKit Butina clusterer ---
    try
        butinaMod = py.importlib.import_module("rdkit.ML.Cluster.Butina");
        mods = emk.util.rdkitModule();

        % Compute full lower-triangular distance matrix using BulkTanimotoSimilarity.
        % RDKit Butina.ClusterData expects a flat list in the order:
        %   d(1,0), d(2,0), d(2,1), d(3,0), d(3,1), d(3,2), ...
        % (0-based i > j pairs, row-major lower triangle)
        % We build an accumulating list of previously-seen fps for efficiency:
        % each iteration adds exactly one BulkTanimotoSimilarity IPC call.
        distList  = py.list();
        accFpList = py.list();
        for i = 1:M
            if i > 1
                % Similarity of fps{i} to all previous fps{1..i-1}
                sims   = mods.DataStructs.BulkTanimotoSimilarity( ...
                             fps{i}, accFpList);
                simArr = double(py.array.array("d", sims));
                for j = 1:numel(simArr)
                    distList.append(1.0 - simArr(j));
                end
            end
            accFpList.append(fps{i});
        end

        % ClusterData(data, nPts, distThresh, isDistData=True) -> tuple of tuples
        % Each inner tuple contains 0-based indices of cluster members.
        pyClusterResult = butinaMod.ClusterData(distList, int32(M), ...
                                                 threshold, ...
                                                 pyargs("isDistData", true));

    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:cluster:butina:rdkitError", ...
            "Butina clustering failed: %s", ME.message);
    end

    % --- Convert Python result to MATLAB cell array of index arrays ---
    nClusters = double(py.len(pyClusterResult));
    clusters  = cell(1, nClusters);
    for c = 1:nClusters
        pyCluster = pyClusterResult{c};
        nMembers  = double(py.len(pyCluster));
        idxArr    = zeros(1, nMembers);
        for m = 1:nMembers
            idxArr(m) = double(pyCluster{m}) + 1;  % 0-based -> 1-based
        end
        clusters{c} = idxArr;
    end

    logInfo("butina: %d cluster(s) from %d fingerprints (threshold=%.3f)", ...
        nClusters, M, threshold);
end
