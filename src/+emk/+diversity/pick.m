function indices = pick(fps, N, varargin)
% pick  Select N maximally diverse molecules using the MaxMin algorithm.
%
%   indices = emk.diversity.pick(fps, N)
%   indices = emk.diversity.pick(fps, N, Metric="tanimoto")
%   indices = emk.diversity.pick(fps, N, Seed=1)
%
%   Selects N molecules from a fingerprint library that are maximally
%   diverse using RDKit's MaxMin diversity picker (Kennard-Stone variant
%   adapted for chemical diversity).  The algorithm iteratively adds the
%   molecule most dissimilar to the already-selected set.
%
%   MaxMin is a greedy algorithm: at each step it picks the molecule that
%   maximises the MINIMUM distance to any already-selected molecule.
%   The first molecule is selected randomly (or by Seed for reproducibility).
%
%   Arguments:
%     fps     - 1xM cell array of py.rdkit.DataStructs.ExplicitBitVect
%               All fingerprints must have the same bit length.  M >= N.
%     N       - positive integer  Number of diverse molecules to select.
%               Must satisfy 1 <= N <= M.
%     Metric  - string (optional, default "tanimoto") Similarity metric.
%               Only "tanimoto" is supported by RDKit's MaxMin picker.
%     Seed    - integer (optional, default -1 = random seed)
%               Fixed seed for reproducible first-molecule selection.
%
%   Returns:
%     indices - double(1,N)  1-based indices of selected molecules into fps.
%               indices(1) is the seed molecule; remaining are added in
%               diversity-maximising order.
%
%   Errors:
%     emk:diversity:pick:invalidInput   - fps not a non-empty cell of py objects
%                                         or N is not a positive integer
%     emk:diversity:pick:invalidN       - N > M or N < 1
%     emk:diversity:pick:invalidMetric  - Metric is not "tanimoto"
%     emk:diversity:pick:rdkitError     - unexpected Python exception
%
%   Example:
%     smiles = {"CCO", "CCCO", "c1ccccc1", "c1ccc2ccccc2c1", "CC(=O)O"};
%     mols   = cellfun(@emk.mol.fromSmiles, smiles, "UniformOutput", false);
%     fps    = cellfun(@emk.fingerprint.morgan, mols, "UniformOutput", false);
%     idx    = emk.diversity.pick(fps, 3, "Seed", 42);
%     disp(smiles(idx))  % 3 maximally diverse molecules
%
%   References:
%     Kennard, R.W. & Stone, L.A. (1969). Computer Aided Design of
%       Experiments. Technometrics 11(1):137-148.
%       DOI: 10.1080/00401706.1969.10490666
%     Ashton, M. et al. (2002). Identification of Diverse Database
%       Subsets Using Property-Based and Fragment-Based Molecular
%       Descriptions. QSAR Comb. Sci. 21(8):598-604.
%       DOI: 10.1002/1521-3838(200211)21:8<598::AID-QSAR598>3.0.CO;2-U
%     RDKit Documentation: rdkit.SimDivFilters.rdSimDivPickers.MaxMinPicker
%
%   See also: emk.cluster.butina, emk.similarity.matrix, emk.fingerprint.morgan

    % --- Parse optional name-value arguments ---
    metric = "tanimoto";
    seed   = -1;
    if nargin > 2
        if mod(numel(varargin), 2) ~= 0
            error("emk:diversity:pick:invalidInput", ...
                "Optional arguments must be name-value pairs.");
        end
        for k = 1:2:numel(varargin)
            argName = string(varargin{k});
            switch lower(argName)
                case "metric"
                    metric = lower(string(varargin{k+1}));
                case "seed"
                    seed = varargin{k+1};
                otherwise
                    error("emk:diversity:pick:invalidInput", ...
                        "Unknown option: '%s'. Supported: Metric, Seed.", ...
                        varargin{k});
            end
        end
    end

    % --- Validate fps ---
    if ~iscell(fps) || isempty(fps)
        error("emk:diversity:pick:invalidInput", ...
            "fps must be a non-empty cell array, got: %s", class(fps));
    end
    for i = 1:numel(fps)
        if ~startsWith(class(fps{i}), "py.")
            error("emk:diversity:pick:invalidInput", ...
                "fps{%d} must be a Python fingerprint object, got: %s", ...
                i, class(fps{i}));
        end
    end

    M = numel(fps);

    % --- Validate N ---
    if ~(isnumeric(N) && isscalar(N) && isfinite(N) && N >= 1 && floor(N) == N)
        error("emk:diversity:pick:invalidN", ...
            "N must be a positive integer, got: %s", mat2str(N));
    end
    N = double(N);
    if N > M
        error("emk:diversity:pick:invalidN", ...
            "N (%d) must not exceed the number of fingerprints (%d)", N, M);
    end

    % --- Validate Metric ---
    if metric ~= "tanimoto"
        error("emk:diversity:pick:invalidMetric", ...
            "Metric must be 'tanimoto', got: %s", metric);
    end

    logDebug("diversity.pick: selecting %d / %d molecules (MaxMin, seed=%d)", ...
        N, M, seed);

    % --- Build Python list of fingerprints ---
    pyFps = py.list();
    for i = 1:M
        pyFps.append(fps{i});
    end

    % --- Run MaxMin picker via RDKit ---
    % MaxMinPicker is a Boost.Python extension class; MATLAB cannot instantiate
    % it via pickerMod.MaxMinPicker() directly (raises py.Boost.Python.class
    % not-found error).  pyrun() bypasses this restriction by running Python
    % code in the Python interpreter's own namespace.
    try
        % pyrun creates picker and calls LazyBitVectorPick in one Python scope.
        % LazyBitVectorPick returns py.rdkit.rdBase._vectint (not a tuple/list),
        % so we convert it to a py.list inside pyrun before returning.
        if seed >= 0
            pyResult = pyrun( ...
                "from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker;" + ...
                "picked = list(MaxMinPicker().LazyBitVectorPick(fps, M, N, seed=seed))", ...
                "picked", ...
                fps=pyFps, M=int32(M), N=int32(N), seed=int32(seed));
        else
            pyResult = pyrun( ...
                "from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker;" + ...
                "picked = list(MaxMinPicker().LazyBitVectorPick(fps, M, N))", ...
                "picked", ...
                fps=pyFps, M=int32(M), N=int32(N));
        end

        % Convert Python list of 0-based int to MATLAB 1-based indices
        nPicked = double(py.len(pyResult));
        indices = zeros(1, nPicked);
        for i = 1:nPicked
            indices(i) = double(pyResult{i}) + 1;
        end

    catch ME
        if startsWith(ME.identifier, "emk:")
            rethrow(ME);
        end
        error("emk:diversity:pick:rdkitError", ...
            "MaxMin diversity picking failed: %s", ME.message);
    end

    logInfo("diversity.pick: selected %d molecule(s) from %d", numel(indices), M);
end
