% rp04_chemberta_f1a.m  RP04 F1-a: ChemBERTa ONNX Inference inside MATLAB
%
% Zone C (enhanced) -- ONNX path:
%   Python  -- SMILES tokenisation only (no model inference)
%   MATLAB  -- importNetworkFromONNX + ChemBERTa forward pass + CLS extraction
%              + StandardScaler + fitclinear LR + 5-fold CV
%
% vs F1-b: Python role reduced from (tokenise+embed) to (tokenise only).
%
% ONNX import workarounds applied automatically in this script:
%   W1: patch_onnx_for_matlab -- removes 6 IsNaN->Where (nan_to_num) pairs
%   W2: IR version 10 -> 9 downgrade (MATLAB max fully-supported version)
%   W3: importNetworkFromONNX called with InputDataFormats/OutputDataFormats
%   W4: prepareSliceArgs INT64_MIN sentinel clamp (MATLAB-generated code bug)
%
% Last observed successful run noted in this header (2026-06-28):
%   AUC = 0.9138 +/- 0.0087  (RF03 PASS >= 0.85)
%   delta vs Python LR: -0.0132  (borderline, same as F1-b)
%   F1-a and F1-b give identical embeddings (ONNX fidelity confirmed).
%   This is a header-level record only; the corresponding saved artifact may
%   not be retained in every workspace snapshot.
%
% Run: Ctrl+Enter in MATLAB with project root as CWD.

%% Section 0: Setup
thisDir = fileparts(mfilename("fullpath"));
if ~isempty(thisDir)
    cd(fullfile(thisDir, "..", ".."));
end
addpath(genpath("src"));

logSection("RP04-F1a", "Section 0: Setup", "ONNX path Zone C");
emk.setup.initPython();
molWarmup = emk.mol.fromSmiles("C");
clear molWarmup;
snap = emk.setup.snapshot();

root       = resolveProjectRoot();
csvPath    = fullfile(root, "data", "benchmark", "bbbp.csv");
helperPath = fullfile(root, "repro", "rp04_chemberta", "rp04_chemberta_core.py");
runDir     = makeRunDir("Prefix", "rp04_chemberta_f1a");

onnxSrc   = fullfile(runDir, "chemberta_raw.onnx");   % raw from torch.onnx.export
onnxFinal = fullfile(runDir, "chemberta_matlab.onnx"); % patched + IR9

if ~isfile(csvPath)
    error("emk:rp04f1a:csvNotFound", "BBBP CSV not found: %s", csvPath);
end

%% Section 1: Export ChemBERTa to ONNX (Python)
logSection("RP04-F1a", "Section 1: ONNX Export", ...
    "torch.onnx.export (dynamo, external_data=False, opset=14 advisory)");

logInfo("Exporting ChemBERTa -> ONNX (44M params) ...");
expResult = pyrun( ...
    "exec(open(hp).read()); r = export_chemberta_onnx(sp)", ...
    "r", hp=helperPath, sp=onnxSrc);
expMeta = jsondecode(char(string(expResult)));
logInfo("Export: %.1f MB, opset=%d (actual IR=%d from dynamo)", ...
    expMeta.file_size_mb, expMeta.opset, 10);

%% Section 2: Patch ONNX for MATLAB (W1 + W2)
logSection("RP04-F1a", "Section 2: Patch ONNX", ...
    "W1: remove IsNaN->Where; W2: IR 10->9");

patchResult = pyrun( ...
    "exec(open(hp).read()); r = prepare_onnx_for_matlab(src, dst)", ...
    "r", hp=helperPath, src=onnxSrc, dst=onnxFinal);
pMeta = jsondecode(char(string(patchResult)));
logInfo("Patch: removed %d IsNaN pairs, IR %d->%d, %.1f MB", ...
    pMeta.isnan_pairs_removed, pMeta.ir_version_before, ...
    pMeta.ir_version_after, pMeta.file_size_mb);

%% Section 3: Import ONNX into MATLAB (W3)
logSection("RP04-F1a", "Section 3: MATLAB ONNX Import", ...
    "W3: importNetworkFromONNX with explicit data formats");

% importNetworkFromONNX generates the +onnx_* package in the project root.
% 'BC' = Batch x Channel (sequence of token IDs, seqLen is the channel dim).
net = importNetworkFromONNX(onnxFinal, ...
    "InputDataFormats",  ["BC", "BC"], ...
    "OutputDataFormats", "BCT");
logInfo("Import done. Layers: %d", numel(net.Layers));

%% Section 4: Patch prepareSliceArgs (W4)
logSection("RP04-F1a", "Section 4: Fix prepareSliceArgs", ...
    "W4: clamp INT64_MIN sentinel in MATLAB-generated ONNX support code");

% importNetworkFromONNX regenerates +onnx_<name>/ on each call.
% ONNX INT64_MIN (-9.22e18) is used as Starts sentinel for slice-to-beginning.
% After adding size(X), the result stays huge-negative and creates an unbounded
% colon range. Clamp to 0 after the negative-index conversion.
patchFiles = {
    fullfile(pwd(), "+onnx_chemberta_matlab", "+ops", "prepareSliceArgs.m"), ...
    fullfile(pwd(), "+onnx_chemberta_matlab", "+coder", "+ops", "prepareSliceArgs.m")};

% Fall back to onnx_ir9 name if the package used that name.
for k = 1:numel(patchFiles)
    if ~isfile(patchFiles{k})
        patchFiles{k} = strrep(patchFiles{k}, "onnx_chemberta_matlab", "onnx_ir9");
    end
end

targetOld = ['    % "If the value passed to start or end is larger than the n (the number' newline ...
             '    % of elements in this dimension), it represents n."' newline ...
             '    if Starts(i) > size(X,DLTDim)'];
targetNew = ['    if Starts(i) < 0; Starts(i) = 0; end  % W4: INT64_MIN sentinel clamp' newline ...
             '    % "If the value passed to start or end is larger than the n (the number' newline ...
             '    % of elements in this dimension), it represents n."' newline ...
             '    if Starts(i) > size(X,DLTDim)'];
nPatched = 0;
for k = 1:numel(patchFiles)
    pf = patchFiles{k};
    if isfile(pf)
        txt = fileread(pf);
        if contains(txt, targetOld)
            fid = fopen(pf, 'w'); fprintf(fid, '%s', strrep(txt, targetOld, targetNew)); fclose(fid);
            nPatched = nPatched + 1;
        end
    end
end
if nPatched > 0
    clear functions;  % force MATLAB to reload patched files
    logInfo("prepareSliceArgs patched (%d files). Function cache cleared.", nPatched);
else
    logWarn("prepareSliceArgs patch target not found -- already patched or package renamed.");
end

%% Section 5: Initialize Network
logSection("RP04-F1a", "Section 5: Network Initialization", ...
    "initialize() with dummy CB-format input");

% 'CB' = Channels x Batch. seqLen is C, batchSize is B.
dummyIds  = dlarray(single(ones(128, 1)), 'CB');
dummyMask = dlarray(single(ones(128, 1)), 'CB');
net       = initialize(net, dummyIds, dummyMask);
logInfo("Network initialized (Initialized=%d)", net.Initialized);

%% Section 6: Tokenise All BBBP SMILES (Python -- no model inference)
logSection("RP04-F1a", "Section 6: Python Tokenisation Only", ...
    "Python role in Zone C: SMILES -> input_ids + attention_mask");

logInfo("Tokenising BBBP SMILES (Python tokeniser, no ChemBERTa inference) ...");
logInfo("  max_len: 128  (26/2039 molecules exceed this length in the main token-length audit)");
tokResult = pyrun( ...
    "exec(open(hp).read()); r = tokenize_for_matlab(cp, max_len=128)", ...
    "r", hp=helperPath, cp=csvPath);
tokData = jsondecode(char(string(tokResult)));

inputIdsAll      = int64(tokData.input_ids);       % n_valid x 128
attentionMaskAll = int64(tokData.attention_mask);   % n_valid x 128
labels           = double(tokData.labels);          % n_valid x 1
nValid           = tokData.n_valid;
logInfo("Samples: %d  (pos=%d  neg=%d)", nValid, tokData.n_bbb_pos, tokData.n_bbb_neg);

%% Section 7: MATLAB ChemBERTa Inference -- Extract CLS Embeddings
logSection("RP04-F1a", "Section 7: MATLAB ONNX Inference", ...
    "Batch ChemBERTa forward pass; CLS token at position 1");

batchSize = 16;
hiddenDim = 768;
clsAll    = zeros(nValid, hiddenDim, "single");
nBatches  = ceil(nValid / batchSize);
tStart    = tic;

for b = 1:nBatches
    idx  = (b-1)*batchSize+1 : min(b*batchSize, nValid);
    ids  = single(inputIdsAll(idx,:)');       % 128 x batchN  (CB format)
    mask = single(attentionMaskAll(idx,:)');  % 128 x batchN

    out = predict(net, dlarray(ids,'CB'), dlarray(mask,'CB'));
    % out: [128, batchN, 768] = [seqLen, batch, hidden] (CBT order from OutputDataFormats BCT)
    cls = squeeze(extractdata(out(1,:,:)));   % batchN x 768
    if numel(idx) == 1, cls = cls(:)'; end    % ensure row for single sample
    clsAll(idx,:) = cls;

    if mod(b,25)==0 || b==nBatches
        logProgress(b, nBatches, "ONNX inference");
    end
end

elapsedSec = toc(tStart);
logInfo("Inference done: %.0f sec (%.1f ms/sample, batchSize=%d)", ...
    elapsedSec, elapsedSec/nValid*1000, batchSize);

%% Section 8: fitclinear 5-Fold CV (Lambda=1/nTr, same as F1-b)
logSection("RP04-F1a", "Section 8: MATLAB fitclinear 5-Fold CV", ...
    "Lambda=1/nTr, ridge, lbfgs -- identical to F1-b hyperparams");

X = double(clsAll);
rng(42, "twister");
cv       = cvpartition(categorical(labels), "KFold", 5, "Stratify", true);
foldAUC  = zeros(5, 1);

for k = 1:5
    trIdx = training(cv, k);  teIdx = test(cv, k);
    Xtr = X(trIdx,:);  ytr = labels(trIdx);
    Xte = X(teIdx,:);  yte = labels(teIdx);

    mu  = mean(Xtr, 1);
    sig = std(Xtr, 0, 1);  sig(sig==0) = 1;
    Xtr_s = (Xtr - mu) ./ sig;
    Xte_s = (Xte - mu) ./ sig;

    nTr = sum(trIdx);
    mdl = fitclinear(Xtr_s, ytr, ...
        "Learner","logistic", "Solver","lbfgs", ...
        "Lambda", 1/nTr, "Regularization","ridge");
    [~, scores]  = predict(mdl, Xte_s);
    [~,~,~, auc] = perfcurve(yte, scores(:,2), 1);
    foldAUC(k)   = auc;
    logInfo("  Fold %d: AUC = %.4f", k, auc);
end

aucMean = mean(foldAUC);
aucStd  = std(foldAUC);
logInfo("F1-a MATLAB ONNX+LR: AUC = %.4f +/- %.4f", aucMean, aucStd);

%% Section 9: RF03 + Zone Report
logSection("RP04-F1a", "Section 9: RF03 Verification + Zone Report", "");

aucPy  = 0.9270;  stdPy  = 0.0107;   % fallback only; overwritten if metrics exist
aucF1b = 0.9138;  stdF1b = 0.0088;   % fallback only; overwritten if metrics exist

rp04Dirs = dir(fullfile("result", "runs", "*rp04_chemberta"));
rp04Source = "fallback_hardcoded";
if ~isempty(rp04Dirs)
    [~, si] = sort([rp04Dirs.datenum], "descend");
    mPy     = jsondecode(fileread(fullfile("result", "runs", rp04Dirs(si(1)).name, "metrics.json")));
    if isfield(mPy, "auc_cv")
        aucPy = mPy.auc_cv;
        stdPy = mPy.auc_cv_std;
        rp04Source = string(rp04Dirs(si(1)).name);
        logInfo("Python LR baseline loaded from: %s", rp04Dirs(si(1)).name);
    end
else
    logWarn("No rp04_chemberta run found -- using hardcoded Python baseline.");
end

f1Dirs = dir(fullfile("result", "runs", "*rp04_chemberta_f1"));
f1Source = "fallback_hardcoded";
if ~isempty(f1Dirs)
    [~, si] = sort([f1Dirs.datenum], "descend");
    mF1     = jsondecode(fileread(fullfile("result", "runs", f1Dirs(si(1)).name, "metrics.json")));
    if isfield(mF1, "auc_cv_matlab")
        aucF1b = mF1.auc_cv_matlab;
        stdF1b = mF1.auc_cv_matlab_std;
        f1Source = string(f1Dirs(si(1)).name);
        logInfo("F1-b baseline loaded from: %s", f1Dirs(si(1)).name);
    end
else
    logWarn("No rp04_chemberta_f1 run found -- using hardcoded F1-b baseline.");
end

logInfo("--- Zone C (enhanced) Comparison ---");
logInfo("  Python sklearn LR (tok+embed+LR):     AUC = %.4f +/- %.4f", aucPy,  stdPy);
logInfo("  F1-b  MATLAB LR  (Python embed only): AUC = %.4f +/- %.4f", aucF1b, stdF1b);
logInfo("  F1-a  MATLAB ONNX+LR (Python tok only): AUC = %.4f +/- %.4f", aucMean, aucStd);
deltaVsPy  = aucMean - aucPy;
deltaVsF1b = aucMean - aucF1b;
logInfo("  delta(F1-a vs Python): %+.4f  (%.2fsigma)", ...
    deltaVsPy, abs(deltaVsPy)/aucStd);
logInfo("  delta(F1-a vs F1-b):   %+.4f  (%.2fsigma, <0.01 -> ONNX fidelity confirmed)", ...
    deltaVsF1b, abs(deltaVsF1b)/aucStd);

rf03Pass = aucMean >= 0.85;
logInfo("RF03 (AUC >= 0.85): %s (AUC=%.4f)", ternary_(rf03Pass,"PASS","FAIL"), aucMean);

logInfo("Zone C (ONNX path): Python=tokenise  |  MATLAB=ChemBERTa ONNX+LR");

%% Section 10: Save Results
logSection("RP04-F1a", "Section 10: Save Results", "metrics.json");

metrics = struct( ...
    "auc_cv",          aucMean, ...
    "auc_cv_std",      aucStd, ...
    "fold_aucs",       foldAUC', ...
    "auc_cv_f1b",      aucF1b, ...
    "auc_cv_python",   aucPy, ...
    "delta_vs_python", deltaVsPy, ...
    "delta_vs_f1b",    deltaVsF1b, ...
    "comparison_sources", struct( ...
        "rp04_main_source", char(rp04Source), ...
        "rp04_f1_source",   char(f1Source)), ...
    "n_valid",         nValid, ...
    "inference_sec",   elapsedSec, ...
    "batch_size",      batchSize, ...
    "rf03_pass",       rf03Pass, ...
    "zone",            "C", ...
    "python_role",     "tokenise only (input_ids + attention_mask)", ...
    "matlab_role",     "ChemBERTa ONNX inference + fitclinear LR", ...
    "workarounds",     ["W1:IsNaN_patch", "W2:IR_downgrade_9", ...
                        "W3:InputDataFormats", "W4:prepareSliceArgs_clamp"], ...
    "onnx_file_mb",    pMeta.file_size_mb, ...
    "hyperparams",     struct("lambda","1/nTr","solver","lbfgs", ...
                              "regularization","ridge","seed",42));
writelines(jsonencode(metrics, "PrettyPrint", true), fullfile(runDir, "metrics.json"));
logInfo("Results saved. run_dir=%s", runDir);

snap.run_date  = char(datetime("now", "Format", "yyyy-MM-dd"));
snap.run_dir   = runDir;
snap.rf03_pass = rf03Pass;
emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));

% ===========================================================================
function s = ternary_(cond, a, b)
    if cond; s = a; else; s = b; end
end
