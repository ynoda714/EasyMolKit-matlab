function result = benchmarkBatch(smilesList, descriptorNames)
% benchmarkBatch  Measure batch descriptor calculation throughput.
%
%   result = emk.util.benchmarkBatch(smilesList)
%   result = emk.util.benchmarkBatch(smilesList, descriptorNames)
%
%   Parses molecules from smilesList, then times emk.descriptor.batchCalculate.
%   Returns a struct with timing statistics suitable for performance logging.
%
%   Arguments:
%     smilesList      - (string array | cell of string | char) SMILES strings.
%                       Must not be empty.
%     descriptorNames - (string array, optional) Descriptor names to compute.
%                       Defaults to all 10 supported descriptors.
%
%   Returns:
%     result - struct with fields:
%       .nMols        (double) Number of input SMILES strings
%       .nDescriptors (double) Number of descriptors computed
%       .parseSecSec  (double) Time to parse all molecules (seconds)
%       .batchSec     (double) Time for batchCalculate call (seconds)
%       .totalSec     (double) Total elapsed time (seconds)
%       .secPerMol    (double) Average time per molecule (seconds)
%       .molsPerSec   (double) Throughput (molecules per second)
%       .tbl          (table)  batchCalculate output
%
%   Errors:
%     emk:util:benchmarkBatch:invalidInput  - smilesList is not string/cell
%     emk:util:benchmarkBatch:emptyInput    - smilesList is empty
%
%   Example:
%     smiles = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"];
%     result = emk.util.benchmarkBatch(smiles);
%     logInfo("Throughput: %.1f mol/s", result.molsPerSec);
%
%   See also: emk.descriptor.batchCalculate, emk.mol.fromSmiles

    % --- Input validation ---
    if ~(ischar(smilesList) || isstring(smilesList) || iscell(smilesList))
        error("emk:util:benchmarkBatch:invalidInput", ...
            "smilesList must be a string array, cell array, or char, got: %s", ...
            class(smilesList));
    end

    if ischar(smilesList)
        smilesList = string(smilesList);
    elseif iscell(smilesList)
        smilesList = string(smilesList);
    end

    smilesList = reshape(smilesList, 1, []);

    if isempty(smilesList)
        error("emk:util:benchmarkBatch:emptyInput", ...
            "smilesList must not be empty");
    end

    if nargin < 2 || isempty(descriptorNames)
        descriptorNames = ["MolWt", "ExactMolWt", "LogP", "TPSA", ...
                           "NumHAcceptors", "NumHDonors", "NumRotatableBonds", ...
                           "RingCount", "FractionCSP3", "HeavyAtomCount"];
    end

    nMols = numel(smilesList);
    nDesc = numel(descriptorNames);

    logInfo("benchmarkBatch: %d mol(s) x %d descriptor(s)", nMols, nDesc);

    % --- Parse molecules ---
    t0Parse = tic;
    mols    = cell(1, nMols);
    for i = 1:nMols
        try
            mols{i} = emk.mol.fromSmiles(smilesList(i));
        catch ME
            logWarn("benchmarkBatch: SMILES parse failed for mol %d: %s", i, ME.message);
        end
    end
    parseSec = toc(t0Parse);

    % --- Batch calculate ---
    t0Batch  = tic;
    tbl      = emk.descriptor.batchCalculate(mols, descriptorNames);
    batchSec = toc(t0Batch);

    totalSec  = parseSec + batchSec;
    secPerMol = totalSec / nMols;
    molsPerSec = nMols / totalSec;

    logInfo("benchmarkBatch: parse=%.3fs  batch=%.3fs  total=%.3fs  " + ...
        "%.1f mol/s", parseSec, batchSec, totalSec, molsPerSec);

    % --- Build result struct ---
    result = struct( ...
        "nMols",        nMols, ...
        "nDescriptors", nDesc, ...
        "parseSec",     parseSec, ...
        "batchSec",     batchSec, ...
        "totalSec",     totalSec, ...
        "secPerMol",    secPerMol, ...
        "molsPerSec",   molsPerSec, ...
        "tbl",          tbl);
end
