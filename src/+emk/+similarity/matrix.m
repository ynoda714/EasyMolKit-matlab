function S = matrix(fps, options)
% matrix  Compute a pairwise similarity matrix for a set of fingerprints.
%
%   S = emk.similarity.matrix(fps)
%   S = emk.similarity.matrix(fps, Metric="tanimoto")
%   S = emk.similarity.matrix(fps, Metric="dice")
%
%   Computes the N x N symmetric pairwise similarity matrix for N
%   fingerprints.  The diagonal is always 1.0.  BulkTanimotoSimilarity /
%   BulkDiceSimilarity is called once per row (N IPC round trips) so that
%   memory usage and IPC overhead scale linearly with N rather than N^2.
%
%   Mathematical definition (Tanimoto):
%
%     S(i, j) = T(fps{i}, fps{j}) = |A_i AND A_j| / |A_i OR A_j|
%
%   Value range: [0, 1].  S is symmetric and S(i,i) = 1.0.
%
%   Arguments:
%     fps     - 1xN cell array of py.rdkit.DataStructs.ExplicitBitVect
%               All fingerprints must have the same bit length.
%     Metric  - string, "tanimoto" (default) or "dice".
%
%   Returns:
%     S       - double (N x N) pairwise similarity matrix.
%               S(i, j) is the similarity between fps{i} and fps{j}.
%               S is symmetric and all diagonals equal 1.0.
%
%   Errors:
%     emk:similarity:matrix:invalidInput   - fps is not a non-empty cell array
%                                            of Python objects
%     emk:similarity:matrix:invalidMetric  - Metric is not "tanimoto" or "dice"
%     emk:similarity:matrix:rdkitError     - bit-length mismatch or Python exc
%
%   Example:
%     smiles = {"CCO", "CC(=O)Oc1ccccc1C(=O)O", "c1ccccc1", "CCCO"};
%     mols = cellfun(@(s) emk.mol.fromSmiles(s), smiles, "UniformOutput", false);
%     fps  = cellfun(@(m) emk.fingerprint.morgan(m), mols, "UniformOutput", false);
%     S    = emk.similarity.matrix(fps);
%     % S is 4x4, symmetric, diagonal = 1.0
%     imagesc(S); colorbar; axis square;
%
%   See also: emk.similarity.tanimoto, emk.similarity.dice,
%             emk.similarity.rankBy, emk.fingerprint.morgan

    arguments
        fps
        options.Metric (1,1) string = "tanimoto"
    end

    % --- Validate Metric ---
    metric = options.Metric;
    if ~ismember(metric, ["tanimoto", "dice"])
        error("emk:similarity:matrix:invalidMetric", ...
            "Metric must be ""tanimoto"" or ""dice"", got: %s", metric);
    end

    % --- Validate fps: must be non-empty cell array of Python objects ---
    if ~iscell(fps)
        error("emk:similarity:matrix:invalidInput", ...
            "fps must be a cell array of Python fingerprint objects, got: %s", ...
            class(fps));
    end
    if isempty(fps)
        error("emk:similarity:matrix:invalidInput", ...
            "fps must be a non-empty cell array");
    end
    for i = 1:numel(fps)
        if ~startsWith(class(fps{i}), "py.")
            error("emk:similarity:matrix:invalidInput", ...
                "fps{%d} must be a Python fingerprint object, got: %s", ...
                i, class(fps{i}));
        end
    end

    N = numel(fps);
    logDebug("matrix: N=%d, Metric=%s", N, metric);

    % --- Build N x N pairwise similarity matrix ---
    % Use BulkTanimotoSimilarity / BulkDiceSimilarity to score each row fp
    % against all N fps in one Python call.  This gives N IPC round trips
    % instead of N^2, at the cost of some redundancy (upper & lower triangles).
    % Use emk.util.rdkitModule() to get DataStructs via importlib.
    S = zeros(N, N);
    % Build Python list from MATLAB cell of Python objects.
    % py.list(cell) only works for 1-N numeric/logical cells; for cells
    % containing Python objects we must append element-by-element.
    pyList = py.list();
    for k = 1:N
        pyList.append(fps{k});
    end
    mods = emk.util.rdkitModule();
    try
        for i = 1:N
            if metric == "tanimoto"
                pyRow = mods.DataStructs.BulkTanimotoSimilarity(fps{i}, pyList);
            else
                pyRow = mods.DataStructs.BulkDiceSimilarity(fps{i}, pyList);
            end
            row = double(py.array.array("d", pyRow));
            S(i, :) = row;
        end
    catch ME
        error("emk:similarity:matrix:rdkitError", ...
            "RDKit bulk %s similarity raised an exception: %s", ...
            metric, ME.message);
    end

    % Enforce exact symmetry and unit diagonal to remove floating-point jitter.
    S = (S + S') / 2;
    S(1:N+1:end) = 1.0;

    logDebug("matrix: done, min=%.4f, max=%.4f", min(S(:)), max(S(:)));
end
