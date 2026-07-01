function result = runMatlabGcn(molecules, foldIdx, varargin)
% RUNMATLABGCN  5-fold CV GCN on BBBP using MATLAB Deep Learning Toolbox.
%   Molecules are padded to max_atoms and stacked into N x F x B tensors.
%   GCN convolution uses pagemtimes for correct batched matrix multiplication.
%   dlfeval is called once per mini-batch.
%
%   result = runMatlabGcn(molecules, foldIdx, Name, Value, ...)
%
%   molecules -- struct array from jsondecode(featurized.json).
%                Each element has fields: x (N x 25 double), adj (N x N double),
%                label (0/1), n_atoms, smiles.
%   foldIdx   -- struct from jsondecode(outer_fold_indices.json).
%                Fields: fold1_train, fold1_test, ..., fold5_train, fold5_test
%                (0-based Python indices; converted to 1-based internally).
%
%   Name-Value pairs (defaults match Python GCN hyperparams):
%     Hidden       128    hidden channel dimension
%     LR           1e-3   initial learning rate
%     NEpochs      150    max epochs per fold
%     Patience     20     early-stopping patience (val AUC checks)
%     Seed         42     base RNG seed (fold k uses seed+k-1)
%     UseGPU       false  wrap arrays in gpuArray if canUseGPU()
%     WeightDecay  1e-5   L2 regularization coefficient
%     EvalFreq     5      val AUC check every N epochs (Idea 4)
%     BatchSize    64     mini-batch size
%     UseBN        true   enable BN (learnable gamma/beta + EMA running stats)
%
%   Returns struct with fields:
%     fold_aucs, auc_cv, auc_cv_std, elapsed_sec, all_preds, n_precomputed
%
% Architecture: GCNConv(25->128)x3 + GlobalMeanPool + FC(64) + FC(1)
% BN:  applyBNLayer -- learnable gamma/beta + EMA running mean/var.
%      bn_stats managed outside model so adamupdate does not touch them.

p = inputParser;
addParameter(p, 'Hidden',      128,   @isnumeric);
addParameter(p, 'LR',          1e-3,  @isnumeric);
addParameter(p, 'NEpochs',     150,   @isnumeric);
addParameter(p, 'Patience',    20,    @isnumeric);
addParameter(p, 'Seed',        42,    @isnumeric);
addParameter(p, 'UseGPU',      false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'WeightDecay', 1e-5,  @isnumeric);
addParameter(p, 'EvalFreq',    5,     @isnumeric);
addParameter(p, 'BatchSize',   64,    @isnumeric);
addParameter(p, 'UseBN',       true,  @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SaveWeightDir', '', @(x) ischar(x) || isStringScalar(x));
parse(p, varargin{:});
hp = p.Results;

IN_CH    = 25;
DROPOUT  = single(0.3);
BN_MOM   = single(0.9);   % EMA momentum for running stats (PyTorch default: 1-0.1)
N_FOLD   = sum(contains(fieldnames(foldIdx), '_train'));
if N_FOLD == 0; N_FOLD = 5; end
hidden   = hp.Hidden;
lr0      = hp.LR;
nEp      = hp.NEpochs;
pat      = hp.Patience;
seed     = hp.Seed;
useGPU   = logical(hp.UseGPU) && canUseGPU();
wd       = hp.WeightDecay;
evalFreq = hp.EvalFreq;
bs       = hp.BatchSize;
useBN         = logical(hp.UseBN);
saveWeightDir = char(hp.SaveWeightDir);

n_mol  = numel(molecules);
labels = zeros(n_mol, 1);
for i = 1:n_mol
    labels(i) = double(molecules(i).label);
end

% --- Pre-compute 3D padded tensors (N x F x B) once before fold loop ---
n_atoms_vec = arrayfun(@(i) size(molecules(i).x, 1), 1:n_mol)';
max_atoms   = max(n_atoms_vec);
logInfo('Pre-computing 3D tensors: max_atoms=%d, n_mol=%d ...', max_atoms, n_mol);

X_all    = zeros(max_atoms, IN_CH,    n_mol, 'single');  % N x 25 x B
A_all    = zeros(max_atoms, max_atoms, n_mol, 'single'); % N x N  x B
mask_all = zeros(max_atoms, 1,        n_mol, 'single');  % N x 1  x B

for i = 1:n_mol
    n = n_atoms_vec(i);
    X_all(1:n, :, i)    = single(molecules(i).x);
    A_all(1:n, 1:n, i)  = normalizeAdj_plain(molecules(i).adj);
    mask_all(1:n, 1, i) = 1;
end

if useGPU
    X_all    = gpuArray(X_all);
    A_all    = gpuArray(A_all);
    mask_all = gpuArray(mask_all);
end
% -------------------------------------------------------------------------

fold_aucs  = zeros(N_FOLD, 1);
all_preds  = cell(N_FOLD, 1);
pat_checks = pat;

t_start = tic;

for fold_idx = 1:N_FOLD
    rng(seed + fold_idx - 1);

    trval = double(foldIdx.(['fold' num2str(fold_idx) '_train'])) + 1;
    te    = double(foldIdx.(['fold' num2str(fold_idx) '_test']))  + 1;
    trval = trval(:);
    te    = te(:);

    [tr, val] = stratifiedSplit(trval, labels, 0.2, seed + fold_idx - 1);

    n_tr_pos   = sum(labels(tr) == 1);
    n_tr_neg   = sum(labels(tr) == 0);
    pos_weight = single(n_tr_neg / max(n_tr_pos, 1));

    model          = initModel(IN_CH, hidden, seed + fold_idx - 1, useGPU, useBN);
    [avgG, avgGsq] = initAdamState(model);
    bn_stats       = initBNStats(hidden);   % running stats outside model (F2)
    iter           = 1;
    lr             = lr0;               % reset per fold; StepLR(50, 0.5) updates in-place
    best_val_auc   = 0.0;
    best_model     = model;
    best_bn_stats  = bn_stats;
    no_improve     = 0;

    for epoch = 1:nEp
        tr_shuf = tr(randperm(numel(tr)));

        epoch_loss = 0;
        n_batches  = 0;
        for b_start = 1 : bs : numel(tr_shuf)
            b_end  = min(b_start + bs - 1, numel(tr_shuf));
            b_idx  = tr_shuf(b_start:b_end);
            X_b    = X_all(:, :, b_idx);
            A_b    = A_all(:, :, b_idx);
            m_b    = mask_all(:, :, b_idx);
            y_b    = single(labels(b_idx));

            [loss, grads, b_means, b_vars] = dlfeval(@lossGradBatch, ...
                model, X_b, A_b, m_b, y_b, pos_weight, DROPOUT, bn_stats, useBN);
            grads         = addL2Grad(grads, model, wd);
            [model, avgG, avgGsq] = adamupdate(model, grads, avgG, avgGsq, iter, lr, 0.9, 0.999, 1e-8);

            % F2: EMA update for BN running stats (outside dlfeval, plain arithmetic)
            if useBN
                for k = 1:3
                    bm_k = single(gather(extractdata(b_means{k})));  % 1 x F x 1
                    bv_k = single(gather(extractdata(b_vars{k})));   % 1 x F x 1
                    bn_stats.rm{k} = BN_MOM * bn_stats.rm{k} + (single(1) - BN_MOM) * reshape(bm_k, 1, []);
                    bn_stats.rv{k} = BN_MOM * bn_stats.rv{k} + (single(1) - BN_MOM) * reshape(bv_k, 1, []);
                end
            end

            epoch_loss = epoch_loss + double(gather(extractdata(loss)));
            iter       = iter + 1;
            n_batches  = n_batches + 1;
        end
        epoch_loss = epoch_loss / n_batches;

        % Evaluate val AUC every evalFreq epochs
        if mod(epoch, evalFreq) == 0 || epoch == 1
            val_auc = evalAUCBatch(model, X_all, A_all, mask_all, labels, val, single(0), bn_stats, useBN);

            if val_auc > best_val_auc
                best_val_auc  = val_auc;
                best_model    = model;
                best_bn_stats = bn_stats;   % F2: save running stats at best epoch
                no_improve    = 0;
            else
                no_improve    = no_improve + 1;
            end

            if mod(epoch, 10) == 0 || epoch == 1
                logInfo('[Fold %d/%d] Epoch %3d/%d | loss=%.4f | val_AUC=%.3f | no_improve=%d/%d', ...
                    fold_idx, N_FOLD, epoch, nEp, epoch_loss, val_auc, no_improve, pat_checks);
            end

            if no_improve >= pat_checks
                logInfo('[Fold %d/%d] Early stop at epoch %d | best_val_AUC=%.3f', ...
                    fold_idx, N_FOLD, epoch, best_val_auc);
                break;
            end
        end

        % StepLR(step_size=50, gamma=0.5): matches Python scheduler.step() after eval
        if mod(epoch, 50) == 0; lr = lr * 0.5; end
    end

    test_auc            = evalAUCBatch(best_model, X_all, A_all, mask_all, labels, te, single(0), best_bn_stats, useBN);
    fold_aucs(fold_idx) = test_auc;
    logInfo('[Fold %d/%d] test_AUC=%.4f | best_val_AUC=%.4f | tr=%d val=%d te=%d', ...
        fold_idx, N_FOLD, test_auc, best_val_auc, numel(tr), numel(val), numel(te));

    % Save best-model weights (MATLAB shape [in, out]; Python forward: H @ W).
    if ~isempty(saveWeightDir)
        w_out = struct(...
            "W1",   double(gather(extractdata(best_model.W{1}))), ...
            "b1",   double(gather(extractdata(best_model.b{1}))), ...
            "W2",   double(gather(extractdata(best_model.W{2}))), ...
            "b2",   double(gather(extractdata(best_model.b{2}))), ...
            "W3",   double(gather(extractdata(best_model.W{3}))), ...
            "b3",   double(gather(extractdata(best_model.b{3}))), ...
            "Wfc1", double(gather(extractdata(best_model.Wfc1))), ...
            "bfc1", double(gather(extractdata(best_model.bfc1))), ...
            "Wfc2", double(gather(extractdata(best_model.Wfc2))), ...
            "bfc2", double(gather(extractdata(best_model.bfc2))));
        wPath = fullfile(saveWeightDir, sprintf("best_weights_fold%d.json", fold_idx));
        fid = fopen(wPath, "w", "n", "UTF-8");
        fprintf(fid, "%s\n", jsonencode(w_out));
        fclose(fid);
        logInfo("[Fold %d/%d] Weights saved: %s", fold_idx, N_FOLD, wPath);
    end

    % Collect test-fold predictions (inference, no dlfeval needed)
    X_te = dlarray(X_all(:, :, te));
    A_te = dlarray(A_all(:, :, te));
    m_te = mask_all(:, :, te);           % plain array
    [logits_te, ~, ~] = gcnForwardBatch(best_model, X_te, A_te, m_te, false, DROPOUT, best_bn_stats, useBN);
    scores_te = double(gather(extractdata(1 ./ (1 + exp(-logits_te)))));  % n_te x 1
    preds = [fold_idx * ones(numel(te), 1), te - 1, scores_te(:), labels(te)];
    all_preds{fold_idx} = preds;
end

elapsed = toc(t_start);

result.fold_aucs     = fold_aucs;
result.auc_cv        = mean(fold_aucs);
result.auc_cv_std    = std(fold_aucs);
result.elapsed_sec   = elapsed;
result.all_preds     = all_preds;
result.n_precomputed = n_mol;
end


% ===========================================================================
% Local functions
% ===========================================================================

function model = initModel(in_ch, hidden, seed, useGPU, useBN)
% Initialize GCN parameters with Kaiming normal N(0, sqrt(2/fan_in)).
% useBN=true adds learnable bn_gamma/bn_beta (identity init: 1/0); excluded from L2.
rng(seed);
if useGPU
    mkW  = @(r, c) dlarray(gpuArray(randn(r, c, 'single') * sqrt(2 / r)));
    mkB  = @(c)    dlarray(gpuArray(zeros(1, c, 'single')));
    mkG  = @(c)    dlarray(gpuArray(ones(1,  c, 'single')));
    mkBe = @(c)    dlarray(gpuArray(zeros(1, c, 'single')));
else
    mkW  = @(r, c) dlarray(randn(r, c, 'single') * sqrt(2 / r));
    mkB  = @(c)    dlarray(zeros(1, c, 'single'));
    mkG  = @(c)    dlarray(ones(1,  c, 'single'));
    mkBe = @(c)    dlarray(zeros(1, c, 'single'));
end
model.W{1}  = mkW(in_ch,  hidden);
model.b{1}  = mkB(hidden);
model.W{2}  = mkW(hidden, hidden);
model.b{2}  = mkB(hidden);
model.W{3}  = mkW(hidden, hidden);
model.b{3}  = mkB(hidden);
model.Wfc1  = mkW(hidden, hidden / 2);
model.bfc1  = mkB(hidden / 2);
model.Wfc2  = mkW(hidden / 2, 1);
model.bfc2  = mkB(1);
if useBN
    for k = 1:3
        model.bn_gamma{k} = mkG(hidden);   % scale init = 1 (identity)
        model.bn_beta{k}  = mkBe(hidden);  % shift init = 0 (identity)
    end
end
end


function [m, v] = initAdamState(model)
% Zero structs matching model structure for Adam first/second moments.
zW  = @(w) dlarray(zeros(size(extractdata(w)), 'like', extractdata(w)));
m.W    = cellfun(zW, model.W, 'UniformOutput', false);
m.b    = cellfun(zW, model.b, 'UniformOutput', false);
m.Wfc1 = zW(model.Wfc1);
m.bfc1 = zW(model.bfc1);
m.Wfc2 = zW(model.Wfc2);
m.bfc2 = zW(model.bfc2);
if isfield(model, 'bn_gamma')
    m.bn_gamma = cellfun(zW, model.bn_gamma, 'UniformOutput', false);
    m.bn_beta  = cellfun(zW, model.bn_beta,  'UniformOutput', false);
end
v = m;
end


function bn_stats = initBNStats(hidden)
% Initialize BN running stats as plain CPU arrays (outside model to avoid adamupdate).
% rm: running mean (init 0), rv: running variance (init 1).
for k = 1:3
    bn_stats.rm{k} = zeros(1, hidden, 'single');
    bn_stats.rv{k} = ones(1,  hidden, 'single');
end
end


function A_norm = normalizeAdj_plain(adj)
% D^{-1/2}(A+I)D^{-1/2} as plain single array (no dlarray); used in preprocessing.
A_hat  = double(adj) + eye(size(adj, 1));
d      = sum(A_hat, 2);
d_inv  = 1 ./ sqrt(d + 1e-8);
A_norm = single((d_inv .* A_hat) .* d_inv');
end


function [loss, grads, b_means, b_vars] = lossGradBatch(model, X_b, A_b, m_b, labels_b, pos_weight, dropout, bn_stats, useBN)
% Batch BCE loss + gradients. Called via dlfeval.
% Returns b_means/b_vars (cell of dlarray) for EMA running-stat update (F2 BN).
% When useBN=false, b_means/b_vars are empty cells {}.
% Loss: log-sum-exp stable BCE matching PyTorch BCEWithLogitsLoss.
%   -log(sig(x))   = max(0,-x) + log(1+exp(-|x|))
%   -log(1-sig(x)) = max(0, x) + log(1+exp(-|x|))
X_d = dlarray(X_b);    % max_atoms x 25       x n_tr
A_d = dlarray(A_b);    % max_atoms x max_atoms x n_tr
% m_b stays plain array -- wrapping in dlarray accumulates unnecessary graph
[logits, b_means, b_vars] = gcnForwardBatch(model, X_d, A_d, m_b, true, dropout, bn_stats, useBN);
y           = single(labels_b);
softplus_nx = log(single(1) + exp(-abs(logits)));  % log(1+exp(-|x|)), always <= log(2)
loss_v      = pos_weight .* y       .* (max(single(0), -logits) + softplus_nx) ...
            + (single(1) - y)       .* (max(single(0),  logits) + softplus_nx);
loss   = mean(loss_v, 'all');
grads  = dlgradient(loss, model);
end


function [logits, b_means, b_vars] = gcnForwardBatch(model, X_b, A_b, m_b, isTraining, dropout, bn_stats, useBN)
% Full GCN forward for a padded batch.
% X_b : N x 25  x B  dlarray   (N = max_atoms)
% A_b : N x N   x B  dlarray
% m_b : N x 1   x B  plain array (1=valid node, 0=padding)
% Returns logits: B x 1 dlarray
% Returns b_means/b_vars: {mu1,mu2,mu3}/{var1,var2,var3} dlarray (training)
%                         or {} (useBN=false or isTraining=false)
b_means = {};
b_vars  = {};

H = gcnConvBatch(A_b, X_b, model.W{1}, model.b{1});  % N x hidden x B
if useBN
    [H, bm, bv] = applyBNLayer(H, m_b, model, bn_stats, 1, isTraining);
    if isTraining; b_means{1} = bm; b_vars{1} = bv; end
end
H = relu(H);
H = applyDropoutBatch(H, dropout, isTraining);

H = gcnConvBatch(A_b, H,  model.W{2}, model.b{2});
if useBN
    [H, bm, bv] = applyBNLayer(H, m_b, model, bn_stats, 2, isTraining);
    if isTraining; b_means{2} = bm; b_vars{2} = bv; end
end
H = relu(H);
H = applyDropoutBatch(H, dropout, isTraining);

H = gcnConvBatch(A_b, H,  model.W{3}, model.b{3});
if useBN
    [H, bm, bv] = applyBNLayer(H, m_b, model, bn_stats, 3, isTraining);
    if isTraining; b_means{3} = bm; b_vars{3} = bv; end
end
H = relu(H);
H = applyDropoutBatch(H, dropout, isTraining);

g = globalMeanPoolBatch(H, m_b);          % B x hidden

g = g * model.Wfc1 + model.bfc1;          % B x (hidden/2)
g = relu(g);
g = applyDropoutBatch(g, dropout, isTraining);

logits = g * model.Wfc2 + model.bfc2;     % B x 1
end


function [H, bm, bv] = applyBNLayer(H, m_b, model, bn_stats, k, isTraining)
% F2 BN: learnable gamma/beta + running mean/var.
% Training: batch stats over all valid nodes -> normalize -> gamma/beta.
%           Returns bm (1xFx1) and bv (1xFx1) dlarray for EMA update.
% Inference: running stats (bn_stats.rm{k}/rv{k}) -> normalize -> gamma/beta.
%            Returns bm=[], bv=[].
% H: N x F x B dlarray, m_b: N x 1 x B plain array.
F = size(extractdata(H), 2);
if isTraining
    n_valid = max(sum(m_b(:)), single(1));              % scalar: total valid nodes
    mu      = sum(sum(H .* m_b, 1), 3) ./ n_valid;     % 1 x F x 1 dlarray
    Hc      = (H - mu) .* m_b;
    var_    = sum(sum(Hc .^ 2, 1), 3) ./ n_valid;      % 1 x F x 1 dlarray
    sig     = sqrt(var_ + single(1e-5));
    H       = (Hc ./ sig) .* m_b;
    H       = H .* model.bn_gamma{k} + model.bn_beta{k};
    bm      = mu;
    bv      = var_;
else
    % Move running stats to GPU if H is GPU-backed (outside dlfeval only)
    rm_plain = reshape(single(bn_stats.rm{k}), 1, F, 1);
    rv_plain = reshape(single(bn_stats.rv{k}), 1, F, 1);
    Hd = extractdata(H);
    if isa(Hd, 'gpuArray')
        rm_plain = gpuArray(rm_plain);
        rv_plain = gpuArray(rv_plain);
    end
    sig = sqrt(dlarray(rv_plain) + single(1e-5));
    H   = ((H - dlarray(rm_plain)) ./ sig) .* m_b;
    H   = H .* model.bn_gamma{k} + model.bn_beta{k};
    bm  = [];
    bv  = [];
end
end


function H = gcnConvBatch(A, X, W, b)
% Batched GCN conv: H = A*X*W + b via two pagemtimes calls.
% A: N x N x B, X: N x Fin x B, W: Fin x Fout -> H: N x Fout x B
%
% pagemtimes is required instead of reshape(AX, N*B, Fin)*W because MATLAB
% column-major storage causes reshape(N x Fin x B, N*B, Fin) to interleave
% pages incorrectly: row k does NOT correspond to node/batch pair k of AX.
% pagemtimes broadcasts W (Fin x Fout) correctly across all B pages.
AX = pagemtimes(A, X);    % N x Fin x B
H  = pagemtimes(AX, W) + b;  % N x Fout x B  (W broadcast: Fin x Fout x 1)
end


function X = simpleBNBatch(X, mask)
% Simplified BN over node axis, mask-aware (padding excluded from statistics).
% Per-molecule normalization -- kept for reference; AUC~0.5 due to info erasure.
% X: N x F x B, mask: N x 1 x B -> X: N x F x B
n_each = max(sum(mask, 1), single(1));              % 1 x 1 x B
mu     = sum(X .* mask, 1) ./ n_each;              % 1 x F x B
Xc     = (X - mu) .* mask;                         % N x F x B
sigma  = sqrt(sum(Xc .^ 2, 1) ./ n_each + single(1e-5));
X      = (Xc ./ sigma) .* mask;
end


function X = crossBatchBNBatch(X, mask)
% BN over ALL valid nodes across the full mini-batch (PyTorch BatchNorm1d-like).
% Kept for reference; running stats absent -> train/inference stat mismatch.
% X: N x F x B, mask: N x 1 x B -> X: N x F x B
n_valid = max(sum(mask(:)), single(1));
mu      = sum(sum(X .* mask, 1), 3) ./ n_valid;            % 1 x F x 1
Xc      = (X - mu) .* mask;
sigma   = sqrt(sum(sum(Xc .^ 2, 1), 3) ./ n_valid + single(1e-5));
X       = (Xc ./ sigma) .* mask;
end


function g = globalMeanPoolBatch(H, mask)
% Masked mean pooling: N x F x B -> B x F
g  = sum(H .* mask, 1) ./ max(sum(mask, 1), single(1));  % 1 x F x B
Hd = extractdata(H);
F  = size(Hd, 2);   % size(X, dim) returns 1 for missing dim (safe for B=1)
B  = size(Hd, 3);
g  = reshape(g, F, B)';                              % B x F
end


function X = applyDropoutBatch(X, p, isTraining)
% Inverted dropout; handles both 3D node tensors (N x F x B) and 2D FC tensors (B x F).
if ~isTraining || p == 0; return; end
xd = extractdata(X);                   % plain array -- reliable size/ndims/isa
if ndims(xd) == 3
    mask_raw = single(rand(size(xd,1), size(xd,2), size(xd,3)) > p) ./ single(1 - p);
else
    mask_raw = single(rand(size(xd,1), size(xd,2)) > p) ./ single(1 - p);
end
if isa(xd, 'gpuArray')
    mask_raw = gpuArray(mask_raw);
end
X = X .* dlarray(mask_raw);
end


function auc = evalAUCBatch(model, X_all, A_all, mask_all, labels, idx, dropout, bn_stats, useBN)
% Batch inference AUC (no dlfeval; inference uses running stats from bn_stats).
X_b    = dlarray(X_all(:, :, idx));
A_b    = dlarray(A_all(:, :, idx));
m_b    = mask_all(:, :, idx);           % plain array
[logits, ~, ~] = gcnForwardBatch(model, X_b, A_b, m_b, false, dropout, bn_stats, useBN);
scores = double(gather(extractdata(1 ./ (1 + exp(-logits)))));
auc    = computeROCAUC(scores(:), labels(idx(:)));
end


function grads = addL2Grad(grads, model, wd)
% L2 gradient to W and b, matching PyTorch Adam(weight_decay) applied to all parameters.
% bn_gamma/bn_beta excluded per PyTorch convention (scale/shift params not regularized).
for k = 1:numel(model.W)
    grads.W{k} = grads.W{k} + wd * model.W{k};
    grads.b{k} = grads.b{k} + wd * model.b{k};
end
grads.Wfc1 = grads.Wfc1 + wd * model.Wfc1;
grads.bfc1 = grads.bfc1 + wd * model.bfc1;
grads.Wfc2 = grads.Wfc2 + wd * model.Wfc2;
grads.bfc2 = grads.bfc2 + wd * model.bfc2;
end


function [tr_idx, val_idx] = stratifiedSplit(idx, labels, ratio, seed)
% Stratified split preserving class ratio.
rng(seed);
lbl = labels(idx);
pos = idx(lbl == 1);
neg = idx(lbl == 0);

n_vp = max(1, round(numel(pos) * ratio));
n_vn = max(1, round(numel(neg) * ratio));
n_vp = min(n_vp, numel(pos) - 1);
n_vn = min(n_vn, numel(neg) - 1);

pos_s = pos(randperm(numel(pos)));
neg_s = neg(randperm(numel(neg)));

val_idx = [pos_s(1:n_vp);     neg_s(1:n_vn)];
tr_idx  = [pos_s(n_vp+1:end); neg_s(n_vn+1:end)];
end


function auc = computeROCAUC(scores, labels)
% Trapezoidal ROC-AUC. Returns 0.5 when only one class is present.
if numel(unique(labels)) < 2
    auc = 0.5;
    return;
end
[~, ord]   = sort(scores, 'descend');
srt_lbl    = labels(ord);
n_pos      = sum(labels == 1);
n_neg      = sum(labels == 0);
tpr = [0; cumsum(srt_lbl == 1) / n_pos];
fpr = [0; cumsum(srt_lbl == 0) / n_neg];
auc = trapz(fpr, tpr);
end
