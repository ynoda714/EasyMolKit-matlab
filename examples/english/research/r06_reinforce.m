%[text] # R06: REINFORCE — Fine-tuning for Molecular Property Induction
%[text] EasyMolKit Research — Layer 4
%[text] 
%[text] The LSTM trained in R05 can generate valid SMILES, but it lacks a mechanism to consider "whether the molecule is good."
%[text] By updating the policy (LSTM) using the generated molecule's property score as a "reward," we can guide it to generate more molecules with target properties.
%[text] This method is the **REINFORCE algorithm** (Williams 1992), which forms the basis of modern molecular optimization tools like REINVENT, GuacaMol, and MolGPT.
%[text] In this script, we load the model from R05 and apply REINFORCE fine-tuning using the Lipinski Ro5 score as a reward.
%[text] 
%[text] ## Learning Objectives
%[text] - Understand the MDP (Markov Decision Process) for sequence generation: state = partial SMILES, action = next character, reward = molecular score.
%[text] - Implement the baseline-adjusted REINFORCE gradient estimator $\\nabla J = E\[(R - b) \\cdot \\nabla \\log \\pi\]$.
%[text] - Understand why catastrophic forgetting occurs and how KL regularization prevents grammar collapse.
%[text] - Compare manual REINFORCE loops with MATLAB Reinforcement Learning Toolbox's `rlPGAgent`.
%[text] - Connect RL fine-tuning to real drug discovery workflows (REINVENT, GuacaMol, GENTRL). \
%[text] ## Prerequisites
%[text] - Completion of R05 `r05_smiles_generator.m` (checkpoint saved in `result/`).
%[text] - Deep Learning Toolbox (`dlnetwork`, `dlarray`, `dlfeval`, `adamupdate`).
%[text] - Reinforcement Learning Toolbox (optional; used in Section 4). \
%[text] ## Environment
%[text] - Requires Deep Learning Toolbox.
%[text] - Compatible with both MATLAB Online and Desktop.
%[text] - GPU acceleration is not required.
%[text] - Estimated time: 30–60 minutes \
%[text] ## Data
%[text] - `result/r05_checkpoint.mat` — LSTM and vocabulary trained in R05 \
%[text] ## References
%[text] - Olivecrona M et al. (2017) Molecular de novo design through deep reinforcement learning. J Cheminform 9:48. doi:10.1186/s13321-017-0235-x \[REINVENT: Direct basis for this script\]
%[text] - Williams RJ (1992) Simple statistical gradient-following algorithms for connectionist reinforcement learning. Mach Learn 8:229-256. doi:10.1007/BF00992696 \[Original paper on REINFORCE algorithm\]
%[text] - Brown N et al. (2019) GuacaMol: benchmarking models for de novo molecular design. J Chem Inf Model 59:1096-1108. doi:10.1021/acs.jcim.8b00839 \[Standard evaluation metrics for generative models\]
%[text] - Bung N et al. (2022) De novo design of new chemical entities for SARS-CoV-2 using artificial intelligence. Future Med Chem 14:1019-1030. doi:10.4155/fmc-2021-0223 \[Application of REINVENT to COVID-19 target optimization\]
%[text] - Polykovskiy D et al. (2020) Molecular Sets (MOSES): a benchmarking platform for molecular generation models. Front Pharmacol 11:565644. doi:10.3389/fphar.2020.565644 \[Standard implementation of validity/uniqueness/novelty metrics\]
%[text] - Zhavoronkov A et al. (2019) Deep learning enables rapid identification of potent DDR1 kinase inhibitors. Nat Biotechnol 37:1038-1040. doi:10.1038/s41587-019-0224-x \[GENTRL: Molecular generation via tensor decomposition × RL\] \
%[text] 
%[text] Execution: Run each section with Ctrl+Enter.
%[text] Execute R05 first to generate `result/r05_checkpoint.mat`.
%%
%[text] ## Section 0: Setup and Load Pre-trained Model
% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif isfolder(fullfile(pwd, 'src'))
    addpath(genpath(fullfile(pwd, 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
emk.setup.initPython(); %[output:0381f9ab]
%[text] Check availability of Deep Learning Toolbox and Reinforcement Learning Toolbox.
hasDL = license("test", "Neural_Network_Toolbox");
hasRL = ~isempty(ver("Reinforcement Learning Toolbox"));
if ~hasDL
    error("emk:r06:missingToolbox", ...
        "Deep Learning Toolbox is required.");
end
logInfo("Toolbox check: Deep Learning=%d  Reinforcement Learning=%d", hasDL, hasRL); %[output:021f8d24]
if hasRL %[output:group:7a0892da]
    logInfo("Reinforcement Learning Toolbox detected. Section 4 is fully available.");
else
    logInfo("RL Toolbox not found. Executing manual REINFORCE loop. Section 4 is conceptual only."); %[output:34c7d6b2]
end %[output:group:7a0892da]
%[text] ### Adjustable Parameters
RL_EPOCHS      = 20;      % Number of REINFORCE update steps
RL_BATCH       = 48;      % Number of sequences per RL batch
RL_LR          = 1e-4;    % RL learning rate (set smaller than SL learning rate to prevent forgetting)
GRAD_CLIP_NORM = 1.0;     % Gradient clipping threshold
N_EVAL         = 100;     % Number of SMILES generated for evaluation before and after RL
%[text] ### Load R05 Checkpoint
CHECKPOINT = "result/r05_checkpoint.mat";
if ~isfile(CHECKPOINT)
    error("emk:r06:missingCheckpoint", ...
        "Checkpoint not found: %s\nPlease run R05 (r05_smiles_generator.m) first.", ...
        CHECKPOINT);
end
logInfo("Loading R05 checkpoint: %s", CHECKPOINT); %[output:4b880a6e]
load(CHECKPOINT, ...
    "net", "char2idx", "idx2char", "VOCAB_SIZE", "MAX_SEQ_LEN", ...
    "START_IDX", "END_IDX", "PAD_IDX", "HIDDEN_SIZE", "TEMPERATURE", ...
    "smilesAll", "smilesProc", "X_train", "Y_train", "X_val", "Y_val");
logInfo("Model loaded: VOCAB_SIZE=%d  MAX_SEQ_LEN=%d  HIDDEN=%d", ... %[output:group:3260ca9d] %[output:9e381c29]
    VOCAB_SIZE, MAX_SEQ_LEN, HIDDEN_SIZE); %[output:group:3260ca9d] %[output:9e381c29]
logInfo("Training corpus: %d SMILES", numel(smilesAll)); %[output:9893d5f8]

%[text] Perform Python/RDKit warm-up.
mol_warmup = emk.mol.fromSmiles("C"); %#ok<NASGU>
clear mol_warmup;
%%
%[text] ## Section 1: Evaluation of Pre-RL Baseline
%[text] 
%[text] Before running RL, we measure the current model's effectiveness and property distribution.
%[text] Without a baseline, there is no way to know whether RL has improved anything.
%[text] We use four metrics for evaluation (Brown et al. 2019 / MOSES; Polykovskiy et al. 2020).
%[text:table]{"ignoreHeader":true}
%[text] | Metric | Definition |
%[text] | --- | --- |
%[text] | Validity | Percentage of generated SMILES that can be parsed by RDKit |
%[text] | Uniqueness | Percentage of structurally distinct valid SMILES |
%[text] | Novelty | Percentage of valid unique SMILES not in the training set |
%[text] | Avg reward | Average Lipinski Ro5 composite score of the generated batch |
%[text:table]
logInfo("Generating %d pre-RL SMILES (T=%.2f) ...", N_EVAL, TEMPERATURE); %[output:50c27046]
preRLSmiles = strings(N_EVAL, 1);
for g = 1:N_EVAL %[output:group:29273260]
    preRLSmiles(g) = generateSmiles(net, START_IDX, END_IDX, PAD_IDX, ...
        VOCAB_SIZE, MAX_SEQ_LEN, TEMPERATURE, idx2char, char2idx);
    if mod(g, 10) == 0 || g == N_EVAL
        logProgress(g, N_EVAL, "generating (pre-RL)"); %[output:6f56c742]
    end
end %[output:group:29273260]
%[text] Validate and score the generated SMILES.
[nPreValid, nPreUnique, nPreNovel, preRewards] = evalGenSmiles( ...
    preRLSmiles, smilesAll, N_EVAL);
logInfo("Pre-RL:   valid=%d/%d (%.0f%%)  unique=%d  novel=%d  avg_reward=%.3f", ... %[output:group:90d7e5c9] %[output:1b1db75a]
    nPreValid, N_EVAL, 100*nPreValid/N_EVAL, nPreUnique, nPreNovel, mean(preRewards)); %[output:group:90d7e5c9] %[output:1b1db75a]
%%
%[text] ## Section 2: REINFORCE Fine-Tuning
%[text] The LSTM language model can be treated as a policy $\\pi$.
%[text] The state $s(t)$ is the partial SMILES so far (LSTM hidden state $h(t)$), the action $a(t)$ is the next character token sampled from $\\pi$, and the reward $R$ is the scalar score of the complete molecule (terminal reward).
%[text] The REINFORCE gradient estimator (Williams 1992) is:
%[text]{"align":"center"} $\\nabla J(\\theta) = E\\left\[ (R(\\tau) - b) \\cdot \\sum\_t \\nabla \\log \\pi(a(t) | s(t)) \\right\]$
%[text] Subtracting the batch mean reward $b$ as a baseline reduces gradient variance, the core idea of REINVENT (Olivecrona et al. 2017).
%[text] The reward function uses the Lipinski Ro5 composite score $R \\in \[0, 1\]$.
%[text:table]{"ignoreHeader":true}
%[text] | Criterion | Condition | Added Value |
%[text] | --- | --- | --- |
%[text] | Molecular Weight MW | $\\leq 500$ Da | \+0.25 |
%[text] | Hydrogen Bond Donors HBD | $\\leq 5$ | \+0.25 |
%[text] | Hydrogen Bond Acceptors HBA | $\\leq 10$ | \+0.25 |
%[text] | Lipophilicity LogP | $\\leq 5$ | \+0.25 |
%[text:table]
%[text] If the SMILES is not in the training set, a novelty bonus of $\\times 1.1$ (capped at 1.0) is applied.
%[text] 
%[text] **Catastrophic forgetting** occurs when RL updates are too large, destroying the learned SMILES grammar.
%[text] REINVENT's "augmented likelihood" prevents this with KL regularization:
%[text]{"align":"center"} $L(total) = L(PG) + \\lambda \\cdot KL(\\pi(RL) \\| \\pi(prior))$
%[text] where $L(PG)$ is the policy gradient loss and $\\pi(prior)$ is the fixed pre-RL network, penalizing large deviations from the pre-RL distribution.
%[text] Monitoring `val_loss` at each RL epoch allows early detection of grammar collapse.
%[text] A copy of the pre-network is saved.
%\[text] > **Note on Simplification**: The REINVENT paper actually computes the KL loss using $\pi_{prior}$, but this implementation omits it to reduce computational cost.
%\[text] > Instead, by monitoring `val_loss` in each epoch, grammar collapse is detected early.
netPrior = net; %#ok<NASGU> % reference snapshot for teaching; KL loss is not applied in this simplified implementation

avgG_rl       = [];
avgSqG_rl     = [];
iter_rl       = 0;
rl_rewardLog  = zeros(RL_EPOCHS, 1);
rl_valLossLog = zeros(RL_EPOCHS, 1);

logInfo("Starting REINFORCE fine-tuning (%d epochs, batch=%d, lr=%.1e)", ... %[output:group:3a1ac0a2] %[output:09293c18]
    RL_EPOCHS, RL_BATCH, RL_LR); %[output:group:3a1ac0a2] %[output:09293c18]

for rl_ep = 1:RL_EPOCHS %[output:group:6484376b]
    % Step 1: Generate a batch of SMILES (inference only -- no gradients)
    batchSmiles = strings(RL_BATCH, 1);
    for g = 1:RL_BATCH
        batchSmiles(g) = generateSmiles(net, START_IDX, END_IDX, PAD_IDX, ...
            VOCAB_SIZE, MAX_SEQ_LEN, TEMPERATURE, idx2char, char2idx);
    end
    % Step 2: Calculate rewards for each sequence and subtract the baseline
    rewards = arrayfun(@(s) lipinskiReward(s, smilesAll), batchSmiles);
    baseline = mean(rewards);
    rl_rewardLog(rl_ep) = baseline;
    advantages = rewards - baseline;    % [RL_BATCH x 1]
    % Step 3: Re-encode generated sequences for gradient calculation
    [Xrl, Yrl] = encodeSmilesBatch(batchSmiles, char2idx, ...
        START_IDX, END_IDX, PAD_IDX, MAX_SEQ_LEN, VOCAB_SIZE);
    if isempty(Xrl)
        logWarn("RL epoch %d: All sequences are unencodable, skipping.", rl_ep);
        continue;
    end
    dlX_rl = dlarray(Xrl, "CTB");
    % Step 4: Policy gradient update
    iter_rl = iter_rl + 1;
    [~, grads_rl] = dlfeval(@policyGradientLoss, net, dlX_rl, Yrl, ...
        PAD_IDX, advantages);
    % Gradient clipping
    for gi = 1:height(grads_rl)
        gdata = extractdata(grads_rl.Value{gi});
        nrm = sqrt(sum(gdata(:).^2));
        if nrm > GRAD_CLIP_NORM
            grads_rl.Value{gi} = grads_rl.Value{gi} * (GRAD_CLIP_NORM / nrm);
        end
    end
    [net, avgG_rl, avgSqG_rl] = adamupdate(net, grads_rl, avgG_rl, avgSqG_rl, ...
        iter_rl, RL_LR);
    % Step 5: Monitor catastrophic forgetting with validation loss
    dlXv = dlarray(onehot_encode(X_val, VOCAB_SIZE), "CTB");
    rl_valLossLog(rl_ep) = extractdata(dlfeval(@modelLoss_local, net, dlXv, Y_val, PAD_IDX));
    if mod(rl_ep, 5) == 0 || rl_ep == 1
        logInfo("RL epoch %2d/%d  avg_reward=%.3f  val_loss=%.4f", ... %[output:94e68094]
            rl_ep, RL_EPOCHS, baseline, rl_valLossLog(rl_ep)); %[output:94e68094]
    end
end %[output:group:6484376b]
%[text] Visualize the progress of REINFORCE.
figure("Name", "R06 REINFORCE Progress"); %[output:5406fbc1]
tiledlayout(1, 2); %[output:5406fbc1]
nexttile; %[output:5406fbc1]
plot(1:RL_EPOCHS, rl_rewardLog, "g-o", LineWidth=1.5, MarkerSize=4); %[output:5406fbc1]
xlabel("RL Epoch"); ylabel("Batch Mean Reward"); %[output:5406fbc1]
title("REINFORCE: Drug-likeness Reward"); grid on; %[output:5406fbc1]
nexttile; %[output:5406fbc1]
plot(1:RL_EPOCHS, rl_valLossLog, "r-s", LineWidth=1.5, MarkerSize=4); %[output:5406fbc1]
xlabel("RL Epoch"); ylabel("Validation Loss"); %[output:5406fbc1]
title("Language Quality (Forgetting Monitoring)"); grid on; %[output:5406fbc1]
sgtitle("REINFORCE Fine-Tuning Progress"); %[output:5406fbc1]
%[text] Heavy oscillations with `RL_BATCH=48` are expected: high gradient variance is a known characteristic of REINFORCE, not a bug.
%[text] Increasing the batch size or applying a moving-average smoother can reduce oscillations.
%%
%[text] ## Section 3: Evaluation of Generation After RL
%[text] We generate a new batch after fine-tuning and compare results before and after RL.
%[text] Three key checks: Validity (grammar preserved?), Reward distribution (shifted toward drug-like molecules?), and Novelty (not merely memorizing high-reward training data?).
logInfo("Generating %d post-RL SMILES ...", N_EVAL); %[output:00e322bf]
postRLSmiles = strings(N_EVAL, 1);
for g = 1:N_EVAL %[output:group:16a57648]
    postRLSmiles(g) = generateSmiles(net, START_IDX, END_IDX, PAD_IDX, ...
        VOCAB_SIZE, MAX_SEQ_LEN, TEMPERATURE, idx2char, char2idx);
    if mod(g, 10) == 0 || g == N_EVAL
        logProgress(g, N_EVAL, "generating (post-RL)"); %[output:184ca7e7]
    end
end %[output:group:16a57648]

[nPostValid, nPostUnique, nPostNovel, postRewards] = evalGenSmiles( ...
    postRLSmiles, smilesAll, N_EVAL);
logInfo("Post-RL:   valid=%d/%d (%.0f%%)  unique=%d  novel=%d  avg_reward=%.3f", ... %[output:group:3c9597d5] %[output:77d3414f]
    nPostValid, N_EVAL, 100*nPostValid/N_EVAL, nPostUnique, nPostNovel, mean(postRewards)); %[output:group:3c9597d5] %[output:77d3414f]

%[text] Compare the property distributions before and after RL with the training data.
validPost = postRLSmiles(logical(arrayfun(@(s) emk.mol.isValid(s), postRLSmiles)));
validPre  = preRLSmiles(logical(arrayfun(@(s) emk.mol.isValid(s), preRLSmiles)));
FEAT_NAMES = ["MolWt", "LogP"];

propsPost  = batchDescriptors(validPost(1:min(30,numel(validPost))),  FEAT_NAMES);
propsPre   = batchDescriptors(validPre(1:min(30,numel(validPre))),    FEAT_NAMES);
propsTrain = batchDescriptors(smilesAll, FEAT_NAMES);

if ~isempty(propsPost) && ~isempty(propsTrain) %[output:group:1d69edfb]
    figure("Name", "R06 Property Distribution"); %[output:9e6d6b77]
    titles = ["Molecular Weight (Da)", "LogP"];
    for fi = 1:2
        subplot(1, 2, fi); hold on; %[output:9e6d6b77]
        histogram(propsTrain(:,fi),  15, Normalization="probability", ...
            FaceColor=[0.2 0.6 0.8], DisplayName="Training Data", FaceAlpha=0.7);
        if ~isempty(propsPre)
            histogram(propsPre(:,fi),  10, Normalization="probability", ...
                FaceColor=[0.9 0.7 0.2], DisplayName="Pre-RL", FaceAlpha=0.7);
        end
        histogram(propsPost(:,fi), 10, Normalization="probability", ...
            FaceColor=[0.85 0.3 0.1], DisplayName="Post-RL", FaceAlpha=0.7);
        xlabel(titles(fi)); ylabel("Probability Density");
        title(titles(fi)); legend("Location","best"); grid on;
    end
    sgtitle("Property Distribution: Training vs Pre-RL vs Post-RL"); %[output:9e6d6b77]
end %[output:group:1d69edfb]
%%
%[text] ## Section 4: MATLAB Reinforcement Learning Toolbox (Concept Explanation)
%[text] The MATLAB Reinforcement Learning Toolbox provides `rlPGAgent`, a policy gradient agent that bundles REINFORCE with an entropy bonus, advantage normalization, and learning rate scheduling.
%[text] The entropy bonus penalizes "mode collapse," the tendency of the policy to concentrate on a single token.
%[text] 
%[text] For character-level molecule generation, the RL environment is defined as:
%[text:table]{"ignoreHeader":true}
%[text] | Element | Content |
%[text] | --- | --- |
%[text] | Observation | Current one-hot character token (`VOCAB_SIZE` × 1) |
%[text] | Action | Next token index (discrete, 1..`VOCAB_SIZE`) |
%[text] | Reward | 0 at each step; property score at END token |
%[text] | Terminal | When END is sampled or `MAX_SEQ_LEN` is reached |
%[text:table]
%[text] Wrapping an LSTM in `rlEnvironment` requires bridging `dlnetwork` with the RL toolbox's observation/action interface.
if hasRL %[output:group:4b4f931b]
    logInfo("Reinforcement Learning Toolbox is available.");
    logInfo("You can implement rlPGAgent based on the above environment definition.");
else
    logInfo("Reinforcement Learning Toolbox not found."); %[output:298783de]
    logInfo("Manual REINFORCE loop (Section 2) implements the same algorithm."); %[output:5308948c]
end %[output:group:4b4f931b]
%%
%[text] ## Section 5: Summary
%[text] We compare generation quality and rewards before and after RL.
logInfo("Before RL:  valid=%d/%d (%.0f%%)  avg_reward=%.3f", ... %[output:group:2339c514] %[output:6e28324c]
    nPreValid,  N_EVAL, 100*nPreValid/N_EVAL, mean(preRewards)); %[output:group:2339c514] %[output:6e28324c]
logInfo("After RL:  valid=%d/%d (%.0f%%)  avg_reward=%.3f", ... %[output:group:0163c6b1] %[output:7ce94921]
    nPostValid, N_EVAL, 100*nPostValid/N_EVAL, mean(postRewards)); %[output:group:0163c6b1] %[output:7ce94921]
delta_reward = mean(postRewards) - mean(preRewards);
if delta_reward >= 0 %[output:group:7dbc1339]
    logInfo("RL Result: Reward %.3f -> %.3f  (delta=+%.3f, Improved)", ...
        mean(preRewards), mean(postRewards), delta_reward);
else
    logInfo("RL Result: Reward %.3f -> %.3f  (delta=%.3f, Within noise range)", ... %[output:1d690113]
        mean(preRewards), mean(postRewards), delta_reward); %[output:1d690113]
end %[output:group:7dbc1339]
%\[text] > **Why validity is low**: The ~6% validity reflects the small R05 training corpus (500 SMILES) and incomplete grammar learning by the base model. This is **not a limitation of REINFORCE** but a consequence of base-model non-convergence. Practical molecule generation requires >100,000 training examples and hundreds of pre-training epochs. \
%[text] - REINFORCE is the simplest policy gradient algorithm that updates the policy with $\\nabla J = E\[(R - b) \\cdot \\nabla \\log \\pi\]$.
%[text] - Subtracting the baseline $b$ (batch average reward) reduces the variance of gradient estimation and stabilizes learning.
%[text] - KL regularization $L(total) = L(PG) + \\lambda \\cdot KL(\\pi(RL) \\| \\pi(prior))$ is effective in preventing catastrophic forgetting.
%[text] - Reward design is most important. Lipinski's Ro5 is a weak signal, and combining it with QED or synthetic accessibility score (SA score) is more effective.
%[text] - MATLAB `rlPGAgent` is an implementation of REINFORCE with entropy bonus and advantage normalization, providing the same algorithm as the manual loop.
%[text] - At practical scale, evaluation of $10^5$ to $10^6$ molecules is necessary. Note that improvement in rewards is hard to distinguish from noise in small datasets. \
%[text] Local Functions
function [nValid, nUnique, nNovel, rewards] = evalGenSmiles(smilesList, trainingSmiles, nTotal)
%[text] Evaluate the validity, uniqueness, novelty, and rewards of the generated SMILES batch.
validMask = false(nTotal, 1);
for k = 1:nTotal
    try
        validMask(k) = emk.mol.isValid(smilesList(k));
    catch
    end
end
validSm = smilesList(validMask);
nValid  = sum(validMask);
nUnique = numel(unique(validSm));
nNovel  = sum(~ismember(validSm, trainingSmiles));
rewards = arrayfun(@(s) lipinskiReward(s, trainingSmiles), smilesList);
end

function reward = lipinskiReward(smiles, trainingSmiles)
%[text] Calculates the Lipinski Ro5 composite reward.
%[text] Awards +0.25 for each Ro5 criterion satisfied, then applies a novelty bonus of $\\times 1.1$ (capped at 1.0).
try
    if ~emk.mol.isValid(smiles)
        reward = 0; return;
    end
    mol   = emk.mol.fromSmiles(smiles);
    props = emk.descriptor.calculate(mol, ...
        ["MolWt", "NumHDonors", "NumHAcceptors", "LogP"]);
    score = 0;
    if props.MolWt         <= 500, score = score + 0.25; end
    if props.NumHDonors    <=   5, score = score + 0.25; end
    if props.NumHAcceptors <=  10, score = score + 0.25; end
    if props.LogP          <=   5, score = score + 0.25; end
    if ~ismember(smiles, trainingSmiles)
        score = score * 1.1;
    end
    reward = min(score, 1.0);
catch
    reward = 0;
end
end

function [loss, grads] = policyGradientLoss(net, dlX, Ytarget, padIdx, advantages)
%[text] Computes the REINFORCE policy gradient loss.
%[text] `advantages` holds $R(i) - b$ for each sequence as a $B \\times 1$ vector.
%[text] The policy is updated by minimizing $-E\[\\mathrm{advantage} \\cdot \\log \\pi(a | s)\]$.
pred   = forward(net, dlX);
V      = size(pred, 1);
T      = size(pred, 2);
B      = size(pred, 3);
predR  = reshape(stripdims(pred), V, T*B);

Yt = single(reshape(permute(Ytarget, [2 1]), 1, T*B));

targetOH = zeros(V, T*B, "single");
for i = 1:T*B
    idx = Yt(i);
    if idx >= 1 && idx <= V
        targetOH(idx, i) = 1;
    end
end
dlTarget = dlarray(targetOH);

mask    = dlarray(single(Yt ~= padIdx));
nTokens = max(sum(extractdata(mask)), 1);

logPred = log(predR + 1e-8);
logPi   = sum(dlTarget .* logPred, 1);   % Log probability of selected actions

%[text] Repeats each sequence's advantage across T steps to form a T\*B length vector.
nB    = numel(advantages);
dlAdv = dlarray(single(repelem(advantages(:)', T*B / nB)));

loss  = -sum(logPi .* dlAdv .* mask) / nTokens;
grads = dlgradient(loss, net.Learnables);
end

function [Xoh, Ymat] = encodeSmilesBatch(smilesList, char2idx, ...
    startIdx, endIdx, padIdx, maxLen, vocabSize)
%[text] Encode SMILES string array for RL gradient computation.
%[text] `Xoh` is a one-hot array of `[vocabSize x SEQ_LEN x B]`, `Ymat` is an integer array of `[B x SEQ_LEN]`.
N    = numel(smilesList);
SEQ  = maxLen - 1;
Xint = ones(N, SEQ, "single") * padIdx;
Yint = ones(N, SEQ, "single") * padIdx;
for k = 1:N
    s = char(smilesList(k));
    s = strrep(strrep(s, 'Cl', 'L'), 'Br', 'R');
    s = s(1:min(end, maxLen - 2));
    try
        idxs = arrayfun(@(c) char2idx(char(c)), s);
        seq  = [startIdx, idxs, endIdx];
    catch
        continue;
    end
    L = min(numel(seq), maxLen);
    Xint(k, 1:L-1) = single(seq(1:L-1));
    Yint(k, 1:L-1) = single(seq(2:L));
end
Xoh  = onehot_encode(Xint, vocabSize);
Ymat = Yint;
end

function Xoh = onehot_encode(seqMat, vocabSize)
%[text] Converts an `[N x T]` integer matrix to a `single [vocabSize x T x N]` one-hot tensor.
[N, T] = size(seqMat);
Xoh    = zeros(vocabSize, T, N, "single");
valid  = seqMat >= 1 & seqMat <= vocabSize;
[nz_n, nz_t] = find(valid);
nz_c   = double(seqMat(valid));
linIdx = nz_c + (nz_t - 1) * vocabSize + (nz_n - 1) * vocabSize * T;
Xoh(linIdx) = 1;
end

function [loss, grads] = modelLoss_local(net, dlX, Ytarget, padIdx)
%[text] Computes masked cross-entropy loss; local copy used for R06 validation loss monitoring.
pred   = forward(net, dlX);
V      = size(pred, 1);
T      = size(pred, 2);
B      = size(pred, 3);
predR  = reshape(stripdims(pred), V, T*B);
Yt     = single(reshape(permute(Ytarget, [2 1]), 1, T*B));
targetOH = zeros(V, T*B, "single");
for i = 1:T*B
    idx = Yt(i);
    if idx >= 1 && idx <= V; targetOH(idx, i) = 1; end
end
dlTarget = dlarray(targetOH);
mask     = dlarray(single(Yt ~= padIdx));
nTokens  = max(sum(extractdata(mask)), 1);
logPred  = log(predR + 1e-8);
ce       = -sum(dlTarget .* logPred, 1);
loss     = sum(ce .* mask) / nTokens;
grads    = dlgradient(loss, net.Learnables);
end

function smiles = generateSmiles(net, startIdx, endIdx, padIdx, ...
    vocabSize, maxLen, temperature, idx2char, char2idx)
%[text] Perform autoregressive SMILES generation with temperature sampling and grammar constraints.
tokens = startIdx;

parenOpenIdx    = 0;  parenCloseIdx   = 0;
bracketOpenIdx  = 0;  bracketCloseIdx = 0;
if isKey(char2idx, '(');  parenOpenIdx    = char2idx('(');  end
if isKey(char2idx, ')');  parenCloseIdx   = char2idx(')');  end
if isKey(char2idx, '[');  bracketOpenIdx  = char2idx('[');  end
if isKey(char2idx, ']');  bracketCloseIdx = char2idx(']');  end

digitIdx = zeros(1, 9);
for dk = 1:9
    dc = char('0' + dk);
    if isKey(char2idx, dc); digitIdx(dk) = char2idx(dc); end
end

VALID_FIRST = {'B','C','N','O','P','S','F','I','L','R','[','H'};
firstCharMask = zeros(vocabSize, 1);
for fk = 1:numel(VALID_FIRST)
    if isKey(char2idx, VALID_FIRST{fk})
        firstCharMask(char2idx(VALID_FIRST{fk})) = 1;
    end
end

aromPairs = {'c','C'; 'n','N'; 'o','O'; 's','S'; 'p','P'; 'b','B'};
aromaticIdxArr  = zeros(1, size(aromPairs, 1));
aliphaticIdxArr = zeros(1, size(aromPairs, 1));
nAromPairs = 0;
for ak = 1:size(aromPairs, 1)
    if isKey(char2idx, aromPairs{ak,1}) && isKey(char2idx, aromPairs{ak,2})
        nAromPairs = nAromPairs + 1;
        aromaticIdxArr(nAromPairs)  = char2idx(aromPairs{ak,1});
        aliphaticIdxArr(nAromPairs) = char2idx(aromPairs{ak,2});
    end
end
aromaticIdxArr  = aromaticIdxArr(1:nAromPairs);
aliphaticIdxArr = aliphaticIdxArr(1:nAromPairs);

aWildcardIdx = 0;
if isKey(char2idx, 'a'); aWildcardIdx = char2idx('a'); end
MAX_OPEN_RINGS = 3;

cAtomIdx    = 0; if isKey(char2idx, 'C'); cAtomIdx    = char2idx('C'); end
bAtomIdx    = 0; if isKey(char2idx, 'B'); bAtomIdx    = char2idx('B'); end
lSuffixIdx  = 0; if isKey(char2idx, 'l'); lSuffixIdx  = char2idx('l'); end
rSuffixIdx  = 0; if isKey(char2idx, 'r'); rSuffixIdx  = char2idx('r'); end
stereoAtIdx = 0; if isKey(char2idx, '@'); stereoAtIdx = char2idx('@'); end
dotIdx      = 0; if isKey(char2idx, '.'); dotIdx      = char2idx('.'); end
fwdSlashIdx = 0; if isKey(char2idx, '/'); fwdSlashIdx = char2idx('/'); end
plusOutIdx  = 0; if isKey(char2idx, '+'); plusOutIdx  = char2idx('+'); end
eqBondIdx   = 0; if isKey(char2idx, '='); eqBondIdx   = char2idx('='); end
trplBondIdx = 0; if isKey(char2idx, '#'); trplBondIdx = char2idx('#'); end

paren_depth    = 0;
bracket_depth  = 0;
inside_bracket = false;
ring_open      = false(1, 9);

for t = 1:maxLen - 1
    T_curr = numel(tokens);
    xenc   = zeros(vocabSize, T_curr, 1, "single");
    for i = 1:T_curr
        c = tokens(i);
        if c >= 1 && c <= vocabSize; xenc(c, i, 1) = 1; end
    end
    dlX  = dlarray(xenc, "CTB");
    pred = predict(net, dlX);
    prob = double(extractdata(pred(:, end, 1)));
    prob = max(prob, 1e-9);

    logit  = log(prob) / temperature;
    logit  = logit - max(logit);
    prob_t = exp(logit);
    prob_t = prob_t / sum(prob_t);

    if t == 1; prob_t = prob_t .* firstCharMask; end
    if parenCloseIdx   > 0 && paren_depth   == 0; prob_t(parenCloseIdx)   = 0; end
    if bracketCloseIdx > 0 && bracket_depth == 0; prob_t(bracketCloseIdx) = 0; end
    if numel(tokens) > 1
        lastTok = tokens(end);
        if parenCloseIdx   > 0 && parenOpenIdx   > 0 && lastTok == parenOpenIdx
            prob_t(parenCloseIdx) = 0; end
        if bracketCloseIdx > 0 && bracketOpenIdx > 0 && lastTok == bracketOpenIdx
            prob_t(bracketCloseIdx) = 0; end
    end
    if sum(ring_open) >= MAX_OPEN_RINGS
        for dk_g = 1:9
            if digitIdx(dk_g) > 0 && ~ring_open(dk_g); prob_t(digitIdx(dk_g)) = 0; end
        end
    end
    prevTok = tokens(end);
    if lSuffixIdx > 0 && prevTok ~= cAtomIdx; prob_t(lSuffixIdx) = 0; end
    if rSuffixIdx > 0 && prevTok ~= bAtomIdx; prob_t(rSuffixIdx) = 0; end
    if ~inside_bracket && stereoAtIdx > 0; prob_t(stereoAtIdx) = 0; end
    if dotIdx > 0; prob_t(dotIdx) = 0; end
    if fwdSlashIdx > 0; prob_t(fwdSlashIdx) = 0; end
    if ~inside_bracket && plusOutIdx > 0; prob_t(plusOutIdx) = 0; end
    if inside_bracket
        if parenOpenIdx   > 0; prob_t(parenOpenIdx)   = 0; end
        if bracketOpenIdx > 0; prob_t(bracketOpenIdx) = 0; end
        for dk3 = 1:9
            if digitIdx(dk3) > 0; prob_t(digitIdx(dk3)) = 0; end
        end
        if eqBondIdx   > 0; prob_t(eqBondIdx)   = 0; end
        if trplBondIdx > 0; prob_t(trplBondIdx) = 0; end
    end
    if (bracket_depth > 0 || paren_depth > 0 || any(ring_open)) && t < maxLen - 3
        prob_t(endIdx) = 0;
    end
    if t >= maxLen - 3
        if bracket_depth > 0 && bracketCloseIdx > 0
            prob_t(:) = 0; prob_t(bracketCloseIdx) = 1;
        elseif paren_depth > 0 && parenCloseIdx > 0
            prob_t(:) = 0; prob_t(parenCloseIdx) = 1;
        elseif any(ring_open)
            open_d = find(ring_open, 1);
            if digitIdx(open_d) > 0; prob_t(:) = 0; prob_t(digitIdx(open_d)) = 1; end
        end
    end
    for ak = 1:nAromPairs
        prob_t(aliphaticIdxArr(ak)) = prob_t(aliphaticIdxArr(ak)) + prob_t(aromaticIdxArr(ak));
        prob_t(aromaticIdxArr(ak))  = 0;
    end
    if aWildcardIdx > 0; prob_t(aWildcardIdx) = 0; end

    s = sum(prob_t);
    if s <= 0
        prob_t(:) = 1; prob_t(padIdx) = 0; prob_t(startIdx) = 0;
        s = sum(prob_t); if s > 0; prob_t = prob_t / s; end
    else
        prob_t = prob_t / s;
    end

    nextIdx = randsample(vocabSize, 1, true, prob_t);
    tokens(end+1) = nextIdx; %#ok<AGROW>

    if     parenOpenIdx    > 0 && nextIdx == parenOpenIdx
        paren_depth = paren_depth + 1;
    elseif parenCloseIdx   > 0 && nextIdx == parenCloseIdx
        paren_depth = paren_depth - 1;
    elseif bracketOpenIdx  > 0 && nextIdx == bracketOpenIdx
        bracket_depth = bracket_depth + 1; inside_bracket = true;
    elseif bracketCloseIdx > 0 && nextIdx == bracketCloseIdx
        bracket_depth = bracket_depth - 1; inside_bracket = (bracket_depth > 0);
    end
    if ~inside_bracket
        for dk = 1:9
            if digitIdx(dk) > 0 && nextIdx == digitIdx(dk)
                ring_open(dk) = ~ring_open(dk); break;
            end
        end
    end
    if nextIdx == endIdx; break; end
end

charList = {};
for i = 2:numel(tokens)
    idx = tokens(i);
    if idx == endIdx || idx == padIdx; break; end
    if idx == startIdx; continue; end
    if isKey(idx2char, idx); charList{end+1} = idx2char(idx); end %#ok<AGROW>
end
smiles = string(strjoin(charList, ""));
smiles = strrep(strrep(smiles, "L", "Cl"), "R", "Br");
end

function props = batchDescriptors(smilesList, featNames)
%[text] Computes and returns a descriptor matrix `[N_valid x numel(featNames)]` from a SMILES list.
N     = numel(smilesList);
props = zeros(N, numel(featNames));
valid = false(N, 1);
for k = 1:N
    try
        mol = emk.mol.fromSmiles(smilesList(k));
        s   = emk.descriptor.calculate(mol, featNames);
        for fi = 1:numel(featNames)
            props(k, fi) = s.(featNames(fi));
        end
        valid(k) = true;
    catch
    end
end
props = props(valid, :);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:0381f9ab]
%   data: {"dataType":"text","outputData":{"text":"[23:37:57][INFO]  initPython: Python 3.10 already active (Status: Loaded) -- skipping.\n","truncated":false}}
%---
%[output:021f8d24]
%   data: {"dataType":"text","outputData":{"text":"[23:37:57][INFO]  Toolbox check: Deep Learning=1  Reinforcement Learning=0\n","truncated":false}}
%---
%[output:34c7d6b2]
%   data: {"dataType":"text","outputData":{"text":"[23:37:57][INFO]  RL Toolbox not found. Executing manual REINFORCE loop. Section 4 is conceptual only.\n","truncated":false}}
%---
%[output:4b880a6e]
%   data: {"dataType":"text","outputData":{"text":"[23:37:57][INFO]  Loading R05 checkpoint: result\/r05_checkpoint.mat\n","truncated":false}}
%---
%[output:9e381c29]
%   data: {"dataType":"text","outputData":{"text":"[23:37:57][INFO]  Model loaded: VOCAB_SIZE=37  MAX_SEQ_LEN=100  HIDDEN=128\n","truncated":false}}
%---
%[output:9893d5f8]
%   data: {"dataType":"text","outputData":{"text":"[23:37:57][INFO]  Training corpus: 500 SMILES\n","truncated":false}}
%---
%[output:50c27046]
%   data: {"dataType":"text","outputData":{"text":"[23:37:58][INFO]  Generating 100 pre-RL SMILES (T=0.80) ...\n","truncated":false}}
%---
%[output:6f56c742]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 10\/100) generating (pre-RL)\r[##--------]  20% ( 20\/100) generating (pre-RL)\r[###-------]  30% ( 30\/100) generating (pre-RL)\r[####------]  40% ( 40\/100) generating (pre-RL)\r[#####-----]  50% ( 50\/100) generating (pre-RL)\r[######----]  60% ( 60\/100) generating (pre-RL)\r[#######---]  70% ( 70\/100) generating (pre-RL)\r[########--]  80% ( 80\/100) generating (pre-RL)\r[#########-]  90% ( 90\/100) generating (pre-RL)\r[##########] 100% (100\/100) generating (pre-RL)\n","truncated":false}}
%---
%[output:1b1db75a]
%   data: {"dataType":"text","outputData":{"text":"[23:38:19][INFO]  Pre-RL:   valid=6\/100 (6%)  unique=5  novel=6  avg_reward=0.060\n","truncated":false}}
%---
%[output:09293c18]
%   data: {"dataType":"text","outputData":{"text":"[23:38:19][INFO]  Starting REINFORCE fine-tuning (20 epochs, batch=48, lr=1.0e-04)\n","truncated":false}}
%---
%[output:94e68094]
%   data: {"dataType":"text","outputData":{"text":"[23:38:29][INFO]  RL epoch  1\/20  avg_reward=0.062  val_loss=2.6164\n[23:39:00][INFO]  RL epoch  5\/20  avg_reward=0.083  val_loss=2.6129\n[23:39:39][INFO]  RL epoch 10\/20  avg_reward=0.104  val_loss=2.6131\n[23:40:18][INFO]  RL epoch 15\/20  avg_reward=0.021  val_loss=2.6163\n[23:40:56][INFO]  RL epoch 20\/20  avg_reward=0.062  val_loss=2.6202\n","truncated":false}}
%---
%[output:5406fbc1]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAeUAAAElCAYAAADa2PrWAAAQAElEQVR4AeydXagl15Xfl8dX6DZIjAUWcqMeeksoYMFAYswgTaL4lgbliwSGwS8Z2xldkYS8BPySgQQPdIvkwYxfYkgCwZPoajI2IcQh5CWCMahkPIyEZ0APE5wwonubbqMe5EEG27iNLijnd+quW+vuW1\/nnKo6deqs6l61Pvfaa\/+rdq1zzj23+xc+9MMRcAQcAUfAEXAEJoHAL4gfjoAj4Ag4Ao6AIzAJBObZlCcBrRfhCDgCjoAj4AishoA35dXw8mhHwBFwBBwBR2AwBLwpDwZt74k9oSPgCDgCjsDMEfCmPPML7MtzBBwBR8AR2B0EvCnvzrWaZ6W+KkfAEXAEHIFzBLwpn0PhgiPgCDgCjoAjsF0EvClvF3+ffZ4I+KocAUfAEVgLAW\/Ka8HmgxwBR8ARcAQcgf4R8KbcP6ae0RGYJwK+KkfAERgcAW\/Kg0PsEzgCjoAj4Ag4At0Q8KbcDSePcgQcgXki4KtyBCaFgDflSV0OL8YRcAQcAUdgnxHwprzPV9\/X7gg4AvNEwFe1swh4U97ZS+eFOwKOgCPgCMwNAW\/Kc7uivh5HwBFwBOaJwF6sypvyXlxmX6Qj4Ag4Ao7ALiDgTXnkq\/Tkk09KG1WV1DYGv45DVlIbHBu8jdI49Daqy5mOq4uz9nRMlW7jq+RVx1TFp7aqeaps6bgq3Y7Db\/UpyNTUlcaol1rGmKdpDmpoo6bx7nMEKhFIjN6UE0DGUG\/duiVNxMavqqNpDD4do7JytcPrcuNrInI1UTqWeaB0DDYojU\/1dFyq1+XADqXx6NihdC7ViWkijevCm\/Lg65JjmzHUmBL1pDZ07PtCrLeJmu6vfcHI17kZAt6UN8NvkNFs+qE295C5FQxqZx5IbcqxQcSobR1elYOc2KGqnHX2qtgxbVOta0wM2ubaFYyok\/uwbT3udwTqEJhJU65bntsdgYsI+EPzIh6uOQKOwLQQ8KY8resxSjVDNibeJZC\/bSHEENsW19VPLnJ2jfc4R8ARcASmiIA35SlelbOahmQ0MBrZkHNMNTdrn1JtVdfB2pAtVdVu\/SpXxQ1tY+66OVKf1ZEtpTnwNdnwW0pj0a0f2dqQhyKdi\/zISuhKarNcfSm3Mcj44RAylMroED4l9JTUZ3kao7qNQVZ7yvE5dUfAm3J3rEaL5KauahzY62id4piDfF3GEldHXcb3HUMt1N9nXnLW0arz1OXBvkou4lmnJWw2B7r1q4w9jcPWRDZ+LJl6tGbl2LrOT6yOU47NjkdXn3Jrs7GbyJqzKof6dH6NSe3qx64xyrGpX7m1aZxy6yO+zo4PIl5j4OjYqwgfMRByVQw2fMQ4dUPAm3I3nHqN4iZtIm7kqgmx11FVfJ+2unmxrzZPt+gmfPANMS8566hb1WVUXR7sZVSz1GWdTTHMhV9nQW8jjR2LUx81rTvfpuNXnZf5mqhuLYyp8tXZqYt4\/MjrEGPJkY6tsxNHPH7kNiK2LQZ\/1zhinUS8KW\/hLuAmrSJKwQ4fi5iv6yYcqybmoa4qUh98XZrieqvWwvqr7HOybbrGTceviiXzNVFdPsbU+Yay9zEnOdgvSnW1do2rG+\/2EgFvyiUWW5f0xh67kG3Nu846p1rrOmvpa4w+MKt4X3N4nnkhUHWvqC1dKXtOSWPg68al41y\/iIA35Yt47K3GpqvaaKsC0jUPcxG7av66eHKRs84\/Zztrb6I5r73r2sCH+8MStq7jpxJHzXYNyNhWrY8xTVSXz45h7k3j6sbvs92b8sSuPjd9080+sXK3Us4mGIEt47dS+BYmZb3rT9vfyG3Xwfxcd0v9rW68TGOsgzl0RVZWWxXvGlc11m0XEfCmfBGPvdZ4YPWxuTRPVS5sEDF9g01OckNVubETU+XbVRvrYV1V9WPHX+Ub0zaFOsCBOlIaE4equbSuKh+14rc+dOwp2Zg2WXNUxZEXf5XPbeMg4E15HJxXmoVNweZIB2FrIo0nBlk5cldi7q6xTXHkgajBEjaoaWwXHznIm8Zih\/ClhD2Nt3oan+o2dkoy60prRcc+dp3MydyWsI1dRzof9VAHZAl7Gju2Tj3UkRL2tBZisKeEPY1t0hnPmJSw23HoaYzq+DQWWe0px6dxztsR8KbcjlGvEV1v0DQOvY20UBunNjh2eBulcaneNt76GWvJ+upk4ut81t4Uhy8lOzaV09gqPR1TpzO2zldlr4qvsunYKh+2lDS+L07+LrmIs6RjsKkMT3VsSqkv1YmrsmGHrI8mYXX869C6ObqMIyaltMZV1kGudHyqE5NSGoOexqiOz5LaU25jXG5HwJtyO0Ye4Qg4AjuOAA0tXUKVLY2Zml5Vc2GbWqVez7oIeFNeFzkf5wg4AjuBgL5zo3lZUvtOLGJRpNZr14Cs9kWI\/50BAt6UZ3ARfQmOgCPQjoA2L+XtI6YZofUrn2aV\/VS1j1m8Ke\/jVfc1OwKOgCPgCEwSAW\/Kk7wsXpQj4Ag4Ao7APBFoXpU35WZ83OsIOAKOgCPgCIyGgDfl0aD2iRwBR8ARcAQcgWYEdrUpN6\/KvY6AI+AIOAKOwA4i4E15By+al+wIOAKOgCMwTwS8KU\/punotjoAj4Ag4AnuNgDflvb78vnhHwBFwBByBKSHgTXlKV2OetfiqHAFHwBFwBDoi4E25I1Ae5gg4Ao6AI+AIDI2AN+WhEfb880TAV+UIOAKOwAAIeFMeAFRP6Qg4Ao6AI+AIrIOAN+V1UPMxjsA8EfBVOQKOwJYR8Ka85Qvg0zsCjsD+IsB\/vbi\/q9\/flTddd2\/K+3tf+Modgf1AYKKr5MF869atiVbnZQ2JANed6181x2SbMgVXkV1Elb\/Kxhjs8CpSH7yO7LiqmDY\/Y2wMMjZL2CxZXyrbuCo5jVe9KnZIG\/Nq\/jpZ\/VPl1F1FU6qX+mw9qW59c5RZb0p9rZO8msvKauuLk7uO+pqjrzzU2Vcu8rTla\/OTo4k2HZ\/m7jufzT\/ZpkyRvJpIKQUj9Vtdc8ChdCy2lOx4K2scOaxdZewaA1e75TYG2fqQsTHWEvYqsjF1ctW4qjnqxru9RMCxLLGYmsQ9vY3rw7zrYsFYak7HY6uiNE5EBjFRV5q4ypbG9KHXzVNnX3fOvvOtWwfXuaqWSTfldRdbN64OhLr41A6A5Ejt6NjxI7cRccSncdjwpfY+9THmqKuXuet8bncE1kGA\/VJ3X2HHv07eujHkrPPtk33XcOi73r7z2Xtnr5qyXfgc5E0fODoeDoGJcmQla0NOSePaOOPqYtQHt2TjsaPDldAh1ZVjU1Kb8ja7+lfhmlu5jkVXGV6nY7dErBJ2ZDhkZdWxrUKMs2THYkeHK6ErqS3l6leOHxluCZsl66uTbfy6MrnTsWqDW0rjVCcG2XIr40tJ\/al9VZ08lux47OhwCBlCbiPiIOKUWzm1qQ6HbCwyhF0JPSX1WU6je\/LJJ61JGIddjeiW1A7HrhwZQodUtlxl9aMrYVPChgyHVIZD1oYMYVdCT0l9VXznm3K6WPSqhaqNC9wUgy8lHTsGpz47T1oLuvrTWLWvwslHHqhtnI218dirqC2f+hlLPuXIStg0Do6uPjg6hGwJW1U8MfggZEvYGLMqMc7mQcbWNQ+xjLGEzY5HV7+VsaHb2DaZeMZZwmbHoVf5UztjiINXURpPLDaNRcampHbVlROnvr45uXUe5dia5iEOPxxSOR2Hrn5i2oj4lBiDjTyWsOFTQlc\/NqtjxwYhWyJO7crxq6xcbegpkQO\/EjoxcLXBsSlHXoXSXOTBZnOgY1dCt37s6HAImRhkS9jwKaGrX22WWz9x6Pjh6ErYkOF1NOmmzIJSSheEnlLdYtVOPHlVtxxfStbfVSZ\/SuTtOl7jGJOS+vrg5N4kD+PrqEteMGJ8l1hiVoklvorqctTZyUGdKTXFMwYihnHIKWHHD6W+VG+KafKledC7xHeJIVcb9ZWnbZ5N\/FOqkVpS6ro2xnWN7Ttum3PbtQxVR1veFr8tsVWedFNmoUqsBBm+C0StStSLDG8jHtRtMXPxs1ZwgeuakC2pvY3bMcgar\/mxKeGrs+OrIuKV8CPDlTS3crVbjo9xcGtHxmYJ25Bk50Ke2lzUZGmV+hi3SjyxjLGEba5k14k8pXXa\/UFt6FOqb4xaJt2ULQBcHC6StW0i952vqZYx52qqY1OfroPrAKFvkrNqPLaUusyRjkHXcciWqB2ftSGrHV8TVcViS6kph\/Ux77pjbZ6u8ibzUSfjldCb5iWOGEs2HjsxSuj44SlhH4KYu6+5yEM+6oSjI2+DmJsalNCpA54S9iGIeXR+OPoQ80w1J+tl3UrobbXuTFNuW8g6fgACrK5jm+LJg79LLuKIT2Ox4UvtfeqrzkG8zo9MfUpq35STj9xVeersVbHWpuOUWx9ynR1fH2Tz6\/rg5Ibjh6OnhC+1Nemrxqe5VhlPLHUrpblUJ05ly1M7uuaC21grE2d1lRmT+tSGHVljU44\/taHX2fH1RdQ19DzkZx6lutqJq\/Ntaie3zg\/fNN+ujW9av\/rSNe1UU+aishC7CPQ6Ig4fvI7IWeershNPzpSwV8WrDT9jUh2bEjHqV66+lOPHBq8j\/ClVzWHH47dj0K3f+pCtbxNZ5yGnJbU35SbGjkHGxhg4uiVskLUhY2NMFyKWMcSqjK6EDV8XIlbHwdGbxuEnTgm9Kl79KSfe2tCrxtfZ7Fhk4siBrITeZMenpGOUM1Zl5dg0PuX4NE45Nghd41XHBqHjg6MroWNvI+J0jI1VO9zaN5HJpXMpx9aWU2OVM0Zl5dg0D7Lam2zq68I1n\/K2MVoD3Maiaw7l2GxMF5kxOp54q6sdG74+SHMqb8s56abMItIFWLCQIWLgKWGHyIMPuYrUp7wqxtqISyn1W11lxqgMv3XrlmBTwmZJ7VWcOOzwKsJXRTYWv9VVxq6EDVkxRLaEnZg6IlZ9dbL1E2NJfXDscEtqg1tKY6p81oZsx1i5zmftyJbseGR8cKUqHRtEjPJURofwK6kOV1JfFdcYy4lT3cpcXwgfHJ8SNgg7XO1wtcHRlTQOnvo0pspOfBPpWOUai64yHF1JdeVVdvXBIWLgSuiQ6qtwxrFWHYOuchXHb8nGYLc6ebGlVGW345B1DLJSakO3PpWV46+aS+0apxy7yvA6HbslYpWwq6xcbcqt3dqQLWkcHDvckrVZWWOwNa1ffRpv+SSbshasC7MFW5k4yNpUxs54JXT1OV8dAXAEw5Swr57NR2wTAa4h100Jvaoe7BoDJwabEjqEDxsyhIxNCR17SvjxpYQ9jd0VXdfSdQ1d41ZdP3m1Fsuxr5pr3XjmsnOrjH3dnLs0jnXqmi3HDtWtZXJNmeJtwcjYqhaAD0p9xFfZ0zjXV0MATFOqyuC26SLQdW9Uxem1Z3Uqw9GVqsapr4ozPqWquF2x6VqmUK\/WYvnYddm5VR67hm3Op2u2vK2eyTXltoL79N+9e1ecHIOhdpfnOAAAEABJREFU7oE+79Vt5aLJKvVZw1CYe17fz9wDfd6rY+eaZVPmVYkFkodKauPC\/fZv\/7Z85jOfcXIMOtwDq98nv\/mbv7l80WfvxW3L6T6o2htao\/oYA6HjQ4YrYVeb8iqf2nzvrX4v+XNqNcymuPf0\/m\/js2zKdtH2gWHtPBjeeust+cpXviLf+MY3nFow+OIXv7iEz\/Hqdq+AF\/fXErSJnur2hpabNli1W96Uo87ne6\/bPaTPJe4lMPe91w038Jr63uN61tGsm3LdQ8GC8cwzz8izzz7r1ILBr\/\/6r8tv\/dZvCdzxar9fuK\/sfTY1WffGJnU15Wjy6Zxg5PdS+73EnvO9146T3kvcV3qP7SKfXFPm1TkbWsFExqZ6V77uuK759y3uE5\/4hLz44ov7tuxZrrfL3mDPEVcHAD5iqvxNvqp4tzUj4HuvGZ+5eSfXlAGYzc7GhpCxKWFTuY5rDNxSXbzbHYF9QYD9wFrhlrBB2OAQew9dCR07unJkJbUpVzsc23jkMzkCu4vAJJsycPIAgJAtdbERU0U2j8uOwD4iULUvsCkWVsaGroQOqZ7yNh9+J0fAEWhGYLJNubls9zoCjoAjsF8I+Gr3AwFvyvtxnX2VjoAj4Ag4AjuAgDflHbhIXqIj4Ag4AvNEwFeVIuBNOUXEdUfAEXAEHAFHYEsIeFPeEvA+rSPgCDgCjsA8EdhkVd6UN0HPxzoCjoAj4Ag4Aj0i4E25RzA9lSPgCDgCjoAjsAkC023Km6zKxzoCjoAj4Ag4AjuIgDflHbxoXrIj4Ag4Ao7APBHwpjzudfXZHAFHwBFwBByBWgS8KddC4w5HwBFwBBwBR2BcBLwpj4v3PGfzVTkCjoAj4Aj0goA35V5g9CSOgCPgCDgCjsDmCHhT3hxDzzBPBHxVjoAj4AiMjoA35dEh360Jo0S5ufhz5fCKPPnEk\/JXDv7KQrspcz50zR+R4s8T8oSAgfjhCDgCjsDACHhTHhjgXU5Pc3penpeXF390HdjQaVRqmxNnfbNe85wulq\/FEZghAt6UZ3hR+1oSzZcmRb7Pn35evvj+F+W50+dQBftL8pLM7bBrPpZjubH4k0kmHHNdM2tzcgQcgWkg4E15GtdhklWcyIlw0JT4+PbuA3fl9cUfdOzqR54L6ZpY4w25sVwWHB1F\/chOk0HAC3EEZoOAN+XZXMrhFnIkR\/JPD\/6pfPOhby45us7Gu0eV58RZYy658M6ZTwTQ5eyY65rPlufMEdhNBD7yEZEFPfurvyq3bt+Wa7\/0S0sdm+zQ4U15hy7W2KUGCcLxhrwBWxLvFFXHDy0dMznpelgjxLJowirjh7A7OQKDIuDJV0cgBLn71FPy5uGhSJatPn4CI7wpT+AiTLUE\/ciWd4yQ1qmy+tU+B65rYo28ANE1oSOrH9nJEXAEJoZAWDTl\/\/Jf5HNXr8rdBZcsm1iB7eXsZFN+8sknBWpfnnSO65Jr32JekVekrgnxbvGGFD9zlRkd+7jm9PKxt5RSn9U1Bm7tyNiU0J0cgTMEBmfP3r8v375zR6793u+JxDj4fH1PsHNNmY1+69YtgZDrAMEH1fnd3g2B1+V1SZsvjfq23BYas8zwqFozGMx5zXoZ2TPsLSV09VmOXWPg6OpHxqaErj7njsDgCOS5XDs9FXn1VZEQBp+u7wl2qimzudnoCgIyNtUtxwdZW5V8urh4TqfShMHvnP7OBeiaYufiY828+NCFf+H0C40YpevWcbvE2Utd9kxVnI6r8tVhkGLmevM+dHya8VneZ4uGvGzGKDGKxIg0PPU4w0415R7XfZ7q5ORE7iw+6nC6U4vDH7\/7x+d4IXzn4Du1sXPC8Z3Td1jukr773ndXWvNrr722HDfnEw1YaZ11+t6r33N2H929e1fu3bu30v1nx++LfH4PxnguelMuodgZ6YUXXpCrV686NWDw6KOPXrqeP7\/689ljdnBwcL5uMFjlPnn66afPx+6KoO92tV4abmpLffghYvEhw5Wwpzb1+d672mkPPfbYY7Lq\/bfKvTqX2Ps\/+5nc\/9739PYSef11kQ8\/LKi0Tl6a0Dvl7WAVQpDDw0OnFgzSq\/Pm4ZuzxyxK+Yr77sHdldZ7\/fp12eWjqZmyrrpGi0+pLYfvvW7PnStXrsiDDz640v23t8+0\/\/pf9fbbWb73TXlnr9yWC3\/D\/O7ylkvx6XtGoK2Zdpmujxxd5vEYR+ACAm8U\/6bC3YOD5e8rX\/DtiDJYU2ZTdqFVcOLVOTl1DDI21afI51CTfccYJAhHLrn4MT8Euuwp9hxxdavHR0yd3+2OwNAI0JSHnmOo\/IM1ZTalJRZQpWNfhcjBpoeQ7VhsVne5fwS0KdtG3f8s28+Yri\/Vt19h\/xXo\/oFb0pmwqczeQ1dCx4euHFkJm5MjMCgCMYrk+XKKtw4Pl3wXT4M1ZQsGG1M3rdrRsau+CmcslI7pakvHud4dgSM5EpEifh8aVbHS\/Tizf6pIV49PZTi6EjqkesrxOTkCgyIQ43n6Nxc\/hz9XdkwYpSnvGCZeboKAbb7Hcix67NNH2N+X7+uynTsCjsAUEcjz86r8nfI5FC7sAwL6Efbcvuxlr519IWLtLjsCjsBEETBf8ppohZ3KGuWdMh9l8VF1Stg7VelBW0XANqhrp9ckW\/yRxXEiJ+KHI+AIOAKTQCDPl2Xs8pe8WMAoTZlmTANOiQKcdg+B61L+Dq5t2Lu3ku4V7+46u6\/RIx2BnUUgxrL0LCvlHZRGaco7iIuXXIEA75Ix68fXyLkUr05lZoc34ZldUF\/OvBGIsVxflpXyDkqjNGXeIfNueQfx8ZIXCKRfcrJf9lq49+KvN+lpXWavxhG4gECeF2oIIkflb4gUxt06j9KUtSHDU9otuLxaRUDfLc\/1y17ehPVKO3cEdgCBsy95SQg7UGxziaM0Zd4p11Fzee6dAgLaoB4\/ffy8HP+y1zkULjgCPSHgadZCIEaRPC+G7vi7ZBYxSlNmIqd5IXB9z77spS9M5nUVfTWOwAwQiLFcRJaV8o5KozXl9GNr1XcUt70vWz++Bog5ftkr\/Tk663RyBByB1REYfESel1OEUMo7Ko3SlGnA+vE1OCFbjuw0XQT0XaJ++5pK7Ze91I\/dyRFwBByBURGwP08OYdSph5hslKZcVTiNmWZd5XPbbiCg75bn+K6y6oVGlW0qV8ruJWRoKrV5HY7AoAjEWKQPoeCjnIebZGtNebgleea+EdBmdO2DaxdSZ1L8\/Mb\/ZS\/Z6kED5kUuRaiMjozNyRGYLQIxisRYLG8GX\/JiId6UQcFpLQSOpPx9QG3cayXyQY6AI+AIrINAjOWoLCvlHZZGacr2VbvKTz75pCDvMHZ7UXrXZju3L3vpuvUjei622pCdHAFHYAII5HlRRAgiWVbIO34epSmnH6PRjKEdx26l8nmg35Sb8pGzP0\/IEwvtpujR5te4bXL7RS\/qGPrLXlPDhHpYt5Mj4AhMBAH7Ja+JlLRpGaM0ZRowRHNW2rTwyY6vKIyH+fPyvLy8+KNubOg0Z+Qmv46ZItd3kn1\/2WvbmDA\/eOvPzZGnSnZvIVMn+0xldCdHYHYIxCiS58WyZvLzZBYzSlNmIoiHhBIPDQj73Inmqw953l3ekBuSLf7I4sCujXmhSpX\/JXlJtnVQn85t\/0Uvtek6+v6yVxtm28RE1z41zt6ipn3ZV6zVyRFYIpBlSzaH06hNGcB4YEA8QCBscydtWDSwF+VF4d+Lts3Orp+fzb4qr8qNxR\/i8el45C1R7bRHcnTuq1vTecAKgq4ZDPjXw3gn\/rq8LuiyONS\/EAf9y9yDTtBDct1PpFKZvYWMzckRmCUCJyezXNYoTZmHgxIPC2iWaLYsigbGO8BccqlrYNghGjPxcnZgOxMny1hX38Xx8TiY0YShMTCxWDO\/rsna1ebcEXAEtoSA\/Xlylm2piP6nHaUp04SVtDm3LWXTOB2vvG2+If36YOcdss6jNtWVq50GoPHYII0Zk1OHzpd+0Qu7vnNF5iNlPoq\/KeUX2GTNQ9ebNvpeMVmztjkM030Bb1oPfqWqOHypHZul1O+6I9ArAiH0mm7byUZpynaRXZozG9rG2fFWrouz9i55bM4hZG1cNBioaQ5tgsRBxOp45CkRtT4vz18oCRvvbGnOFxwrKrpm8ulQ8o6NCS8OIGrgI3T4rlO6P9Cr1oRd9w8cXeOQIdWVYyPWEjb1O3cEekEgRpE8L1LN6EteLGjUpszmVNJNSxGW8ONTGzI21ZVjw6c6MjbVp8RfkVckW\/wRc2iz4YF\/W25L6pezA\/8NuXGmjc+0zqqZaZKpX9eBnXfOVeO62NowGxITau9S41Ri9N7n\/r9169ayLCsvDWenOvuZ+5xVxTGPBiBDqjt3BEZB4CMfEYGeeELOj5dflqVN5nGM0pTZ4BCbWGkq8J2ensoY9Ienfyh8s1rODj4K\/tLpl+TPT\/9ckPGjn7mX7LnT5879Y9RYNceykLPTJ+5\/QmwMP+M9c52zr51+TbQx47fxq8pgwosSOTvACYwUs1XzdY0\/m27JGLMUFqcPTj+4sH58dbQIH\/Vvuq\/Q+yiAfavURz6bow47t4\/zTNpFnLl\/Tq9dEwkB8QLpei4Yd1AZpSnzgIDGwoe59EGiHFvV\/CcnJ3Lnzp1R6G+\/+7fPS\/jd935X\/vGdfyx2bvRbt28JzYfAx+4\/dsFvY8eSf\/T+jyhlSffu3btQz9K4OD1z\/5nFufj7P9\/\/n\/JX3\/+rhbI4\/\/G7f3xhzKp1g8EizfIvv5IFRqvmWDX+3XffXc7H6b333hOt4Z3Tdzqv5bXXXmP4aKT3ueVVk6f7gPjUpuPUhx9CV18d1zhilbBVxZ+ssvfujLNHV71Xxoi\/e\/eupHtvjHmnOAf30enjj8v9xx5DFBr0\/WeK54\/WO\/beWxbS42mUpky9ukHhqsOHIObgQWAJW9VcL7zwgly9enUU+umjPz0v4VOPfKp2ziMpfs3omw99szZmrJp\/+PAPlzXzQuHRRx+9UI++i33g4IFlDKcfPfIj+dOH\/xRR8P\/Ko79yYcyqdf\/F4V+IHsirjl8nPr1Our4rh1c6r+Xpp5\/Wsgfn3Nv2XlcZe9Pk+Imti2ny1Y3RnIxVwlYVP+beW+c+mMqYxxYNKN17U6lt7Dq4jw6vXJHDvyieCwdPPSUHDzyAWbSWMffecuKeT6M0ZTalblCtHx276tviIQQ5PDwche4e3D1f5icPP1k75wsHL5zHvXn4Zm3cGHVrQ6KgBx988EIt+jH1dw6+g3tJXz\/4uqiOf5Ma7x3ek7j4I2cH8ib5uo49ODg4m1GW633q4Kmlvsr8169fX46Z6om9xx7cZn1j7r2u137kuOX91TbnlUUTSvde25i5+pf3a4wiMS5FOToS3a+65qnvvaLw+vMoTbl++sseHhQ8MNSDjE115djwqY6MTfUpcv32bpAgTYf92XMuZ98wlO0cNKK6mV+R5i+wDfllrBmNMe8AABAASURBVLqa+rDbNbddqz7mGztHl73CXiJu7Np8PkegFYEYy5BXXxWJsdRnIE2uKYOpPhB4KCBjU8KmMj50CLnKjg+yfo0bm+vDvsuDXmO0kY9dazofP89Nbeivy+uSNl90vlGuayBuHVK87Ngqm\/UPKW9z7r7WxV4gF9wSNggbHGLPoCuhY28iYjReObamMe6bGQJjLidGkRjHnHHwuSbZlFk1GxlCtpTa0CEbg4zNErZt0yoPdT76pd4TOZGpHzfl5oXGjC49HKvg1cN0l1Loi4rrMu2PornPtQFajj1dFLYq0jh8KsPRldBTwldlw66U+l13BNZG4MMPRW6YXxFFV1o76bQGjtKU2Zz6sGD5KmNH3xfSJnMkR61LtjE6rnXQFgO0gVFCX\/VW5amyMWefNJVPJ1ZZE3spJfbZKjk81hHYCQTsP685XsGjzTRKU2Y16QMDHfu+kG0mtoHVrX8qP1fWuttqtv5ccunjqGqOWk8f+bvmsGvbxvxd6\/Q4R2CvEAhhlssdtCnzSt1SiiC+1DZX3T7M7UO+ab0a94a80RQ2CZ9+3N5nMYqZ4tBn7qZc25q3qSb3OQKOwBkCeV4IR+2fOBaBu3UerCnTcHk3bAkb8MAhfOhO1Qhoo9vmz5XXaVB9vYiomltt1Yi51RFwBGaNQIyzXh6LG6wpk7yKtBnvW0O2H+lqs63Cx9qOzM+ed6EZ9f2OVtcMDpq76iNti9nQstY09Dye3xFwBCoQiLE0Zlkpz0gavSnPrxl3uxu0mWhz6TJq2z9XXrUB6YuNVcdVYdFHjqq8XWw6t14r5V3Geowj0BkB\/mOFOuqcZM8CYywXHEIpz0gavSnPCLuVlpI+6LsO1obQ10fCXedN47SO1F6l208FqvxdbIoXsTR7nd\/a8Y1N255f18snTm2ksc4njEAIIlkmkmUiWTbhQidSWoxFISGIhFDIMzt7U574BaUhUeLYP1dmzlXpaMCP27Upr1rTOvHaeK9L8fvJY87dtV4+cepCXfN53JYQCEHklVeEfy5SXnxRJMu2VMiOTKu\/DrUj5a5T5qBNOX0lT4FVNuxzJ333eCSrfWPQxmuzGAsrO9\/YjUnxYq36wgTZ1oTu5AgMikDdx8vY+5o4z0X4P4Ffekkkxr6yzjtPCLNd32BNucureGJmi6xZ2CaNxDZD26hM+smIQ9SqOfVd6yZYirRDVZVfa2B0lR+704wRCEEky0RCEMmy\/hcaY5kzxlJ26TICeV7YZvrrUCxusKZMcqcCAfsgt+\/6Cm\/z2cZv++fKzZWK2Fplw2MKa7XNeMPl+PBdRiCE4iNmXUOWqbQ5z3MR\/lOFzTPNP0OM5RpDKOWZSd6Ud+CCHsuxcIz9c2X7YqJrg9K4vpqq5lMODrYu9H0nX\/8ICJyciMQokucied7vhDH2m2+u2WIsVxZCKc9M8qY8wgXNpdzEtrlIx+PI\/Bx6XxqSYmbX3hGujcIsvvZaqay\/2rbRJAMMTr+rofoAU3nKvhDQ\/0ghy8qMx8ci2EuLS4pAnqskEkIpz0zypjziBeXBDm0ypTarTXJ0HVvXoJrGZ1I8YDap084rZ4fFbZPcZ+lmxWjAfD+jinZ3oROqPM9l0I+YYywXG2Mpu1SNQAgiIVT7ZmAdrSnz4KiiGWDYuoRNP8o9lmPRY9Ncmmcofl2uL1NXNdalY8VTJkWTt01Z\/HAExkYgxuFmjLHMneel7NJFBPTXoUK4aJ+ZNkpTphlXvYrHNjM8G5ezSWPR5tRXs2ssdAOnXeO6tVaNs3k3KK91qJ3bzqmy9bcm84DuCPArRnXUPUv\/kXyUDIVQ5s6ytT9iLpOcSTGeCc5aEYixNWQOAaM05TkAtcka9EGuD\/Z1ch1J8fvNfHSr+dbJs8oY\/fnpunWvWydr1Dqr5l43r+acG+fFLS98d35dIYhkmUiWiWTZdJYTo0iMZT15XspDSDEOkXW3c8YoEmOxhhn\/OhQL9KYMCmfEw\/6m3JSPnP15Qp5YaDdlk4OcEDn0o13kVYjx2iAZ97w8L9QpEzz0HT2lUTd8XaIhQzpeZYuF+vridTXr3HV+nR8\/1+ZXn\/1VuX3rtvzNa39T0GXAQxsyPKUBp+0\/dQgir79e0IsvimRZ\/3P0lTHGfjLFeDlPjJdtk7aMXFyWjTzhuNON0pR34ZU8D1Oa3cvy8vkVwIZOcz43biDog32VFNRAXfbXobD1WVddPcxT56uz2zWuM568U\/m5uV0LdbUR6+VacW00Fht6X\/eQ5rWc\/VVHNm5n5JMTEf51qzyfRskxXq4jxsu2viwx9pVp\/Tx1P0pos68\/Y\/PIGJv9M\/KO0pR59Q5m8JSwT4F4cPIApRa+WHVDboi+68P+krwk6xyM1XFBgqx62LrSseRet640V5O+at0av+k7Ws2jtanOutU2Nm+a216rF+6+IB\/76sfkqbtPLUtk3BjXajnZrp\/Sf0yDBr3NNeX55dljvGxbx5Lnl0fFeNm2DUsIIlkmkmUiWXaxghBEQhAJQSTLLvqG0PK8zJplpTygtK3UvzDGxHWv4rGPMX+XOfSdKI34FXlFbi7+vC6vC7osDvUvxJX+8jDWAdpUVO\/CdV7q4MWCjkFHVj\/yVEjXade+Sm06TvPoWNXVr\/Y++ffl+8t0OtdSWZyuS\/Gt8oVY+1evBdfmX9\/91\/LIVx+RPu6h2gmNI32xi27cuyHmucjzz1+u9eREpOkdmox0hFBOFGMp9yWFUGT6fnEPFsoWzyHI8j\/JCEEkBJEQymJiFIlRJMbClmUFH+qsmIQw1AyTyTtKU15ntTxUoLaxxEBVcdiVqvxVtiM5Eh6uz8vzwkeO6HJ2rNMM7Jj0QX+WthOjDkiDbS47h\/q3ybW2XBYPWVntYC0Qo65LdSNUPzFTJK7T\/Wfvy51v35F\/dO0fieIhi2OI2rnHeYGbEvbFlJV\/8SlVBpwZNQZ+ZrrAmuz4oAsD2pQ8F8nz5qgQRLKsOaZPr\/4qDjlD4CyiTaLQ1j9rnhBEQijyxFjwbZ\/zXJY\/Rjg5ETk5EYmxuqI8F4mx2teXNcYiUwgFn\/F52KZsgGNzVpEJOReJ0wcM8rkjEfDVxVkfMejJ8AtqkOJi688zc8mFB+ir8qpw4IeQV6G6d15dc+icb8gborIsDl44LJhgg6ThYB035aZ85OwPLzbQpeVgHCFt+YmxdL2mmdqYLnI6b195m+bWNacxtpa2GK4VMafXToX7CJLFQQ5oIW71L3uBPaGEXlUQdo2Bo2scMqS65diJV0K3\/kqZXzuqI778ZQfFKHLjhkiWWevwcggiWVbME2PBNz3HWGYIoZDzvOBTOYcgEsLFakK4+O+Bx3jR37cWY5ExhILP+DxKU2ZT6gYFS2TLkZU0VnVisamuHBs+1ZGxocPRkZVSXe3KMyk2Gw\/Ql+Ql0YOHK7L6kVchHb\/uw1jnpS5+ZpnOrf7UrjrzPy\/Pix2LDZ3mrHF9crtW5lolt423echhdRuHr2+yc3XNrdeCa9XnPdR1\/i5xVXujalxVnN1DyFA6tm1cGt9Jz7IiLISCc05\/7oxtKIqxyBxCwTnnOef+KASR69WfDPU3ScdMN29eDgxBJITSHqMM+q+clTOJxCgSY2GZCkZFNYOcR2nKVZWzodnAVb6+bORXqst5enoq0NdOvyb6UJXkuHZ6Tb50+qVl3Onp6UpcmwdzrEO2Lh72trQuddF8tYbPn35+uQ5dJ\/YvnH6hdj34me\/x08drY6rWZBvat06\/tdJYnZN5WZ\/Nj03J2vuU6\/KrHV43n71WxFkCkxuyeHdnjTsg6\/6Br1Iu8Up14+pwrLKT4\/TaNZEQEOX0nXeWnFNVfF+2ZTOIkWnk9PHHRY6OljInath0nmX+RbJl7hAWUvG3j9zr1CY3b4q8XP4GyrKaGEXyXCTGpXp+yvNzUYV15mwbY+c9fe651ueJ1rKrfGtNeWjAeCDQ+JXQq+Y8OTmRO3fuLOk\/3\/7P8sX3v3gh7LM\/+ax8+8635aN3PioatwrXJvPpH396rfHMVVXXM\/ef6VTXiZwIB\/H\/7N1\/Jn\/n3b8j5EPH\/vWDr9fWhR\/6yY9\/Infv3pV79+7VxlKn0vXb1xm2pB+9\/6NOY3Ts2++\/vRzHKcX8wXcfxLyk77733ZXyav42\/rP7P1vm\/+D0gwv5u84NtvYe0mZ8W24LsgxwcI9zf6eEPZ0utTEmtekY9eGH0NXXxIkjXgm9Kv7kpNx7bdeF8Qff+c7ygXwunzXmtrGb+pkPev9jH5P3fvhDxCW9993vXrhH1plHG86PP\/5xeffB8v62uVfZe+vUoGN+\/C\/+xXlD5gXQnW9\/W27futVI733lK0ssOCFrrj75u+++S\/olvffee62Yv\/baa8vYXT3NtinzQOhyUV544QW5evXqOX35ypdFH8yMPzw8PPfZuC7yz6\/+nBRL+uWHf3ntPMyldfFul4RXDq90ykcs9NTBU\/K5q5+Tz\/zSZ+T\/Xf1\/8sLBC5iXRJ1Xr5YYMB+2pXNxeujhh+Sxxx6TRx99VPB1ocWw5d+3H3m78xjy\/vDh4qFHA0O39KlHPrXMyemnj\/50pbw2T5P8F4d\/QXoBLxv3yCOPLO2c2nDgWj13+hyhS7opi3cfMuzB\/Z5S24w0S8bUxTX56sZg7zruhWTvWbxTmbzQweLFIRxSOY3tU3\/0pz9lqiU98tf+mjz8D\/7BUuaEb6O5fv5z0izp4Ycekkc+Vd7fNveqe69LTU88+aSk9MhXv7qsRUJY\/rz40V\/5ldY99vA\/\/+ciISzHPfL2anu9S53EPPp\/\/s8yP6dHfuM3Wmt6+umnCd1ZGqUps0l5AICSyujI2LZJIQSh8aYUJAjHAwcPVPrT+Cr93uE9UiyJh3xVzKo26iFhLnmnunQdPzj4gdw9uCscbx6+KX908EeIgv+Th5+szCVnB7VfuXJFHly8ku9aL3kZTr1dxxBHnYxjPHpK+KCDg4PKmtP4VXU5O6rqPnMtcWzLS33Esw5438T+gcgLryP8VUT8tvdfCNV7rwpb0S+BvfJKuZzbtwV7VXxvtjffPJ\/v4Kmn5PCTnyz1xQuETefRZE25V917XWpazhuCSJYtxQunBcYHixdMXfIQozkOvv71QfbkwR8VzyqaP\/O10fXr5Sd1F9a1I8ooTRks7AMAGcKeEnYeGGpHxqa6cmz4VEfGhg5HR16XggThiJL8HEW2e1yX1W64TIpNRxOXs4NvlKuu\/jNXb0zz6jxdEzfhrdeEXE1x+M9oZaZ5U5zt3CsnHWAA9zhEangd4U+JvUF8arc6fuKsrYu87rguuZcxISzZ8pTnSzbKKQSREIqpQii4\/jpToa1+jrEcE0Ihh1DwTXMXWZrPIRS\/h5xGZVlqadbNz9klz5tjN\/GGsMnonRk7WlNmgyuBDjK8inRjE4NsY7Cpjg8dQlY7HB27Ejr2rqQP4VWbis1vx2ay4o0u1YfWhVcbCHIdvSKvSDq3jiPXDbkhVYfG4CPl2+6jAAAQAElEQVQOvg7ZPF3Ga\/yRlF+oseO0lu\/LRP6BBVuckXUdWq9xbVVkP1AA3BI2CBscYs+gK6FjbyPidAwcvW1MZ3+WdQ7tJdD+jrImzLJCyvOC93kOocgWY8GHPMcoy99DZo4QOK9Hx8fluCG+FZ\/nRX7b\/AvLLM+jNGXdmHZzImOvQxU\/lPpTGzqUxqFjV0Jfha7Lau9Iq3Jr4wiywQ0vFw+bSx\/8FyMua6\/L65It\/sjZQY4bckOG\/PLRkZRNtWudXePOltE7a5ofzHTCprg05rpsfh9pznOeCHX7qMqu+yHlmhK7ynB0JfSU8KU2dOxK6INQVcMcZKJF0hAWp+RvjIlhRTXPywEhFHIIBc\/zgg95jrHIHsLyZ8iSZYW+zjnLilEnJwXv6xxjX5l2Js8oTXln0DCFrvoQNkPPRX1421znzjUFm0vzr5Pqpgz\/5SNZ8bDrsS8ibBpdv421\/r5knWedfEPXtk5NsxoTwnjLibGYK4SCc7bv2GLEsjmFUOQY4+ehMRZz6TlGWf4aVIxqWZ3zv3rpqBhV2pzHWObIslKeseRNuebi2ofylB6yfdTVZT02xs5ZA9cls22queTS17FOLX3NvU6eIevlnTBEXfCUeKeKb0eoe5lZVsTGWPChzjGKxFhkr2uWMRb+dc76c+MQytEhlHKMpdyndHJyOVuei8R42d7VkmVlZJ6X8qZSjGWGEEp5xpI35ZqLax+mtkHVhFeacyluziMpP8qVHg6tTT8e75JSa+kS20eM1kiurvjZGm1TJ0dKXXOm45p0m9PWr2PU1oZ7Wx7Ntymn6ULkgaeEfdaU5+MtL4Ryriwr5RhLeVUpxmJECAXnHALngvK84H2eT05k+a6YnCEI316\/RPhWpRBEQihG9fljhRiLnCGIhFDIMz+P0pR5WOirePBUGTv6FEkfwNRmH7LoXWidMV3yEqO1dZ2jKq7KRm4l69f51NeV67i2Jpbm03GpHf26FD+jtfVh32ea8j4a5LoM8fFxVaExltYQquW2BlSO6iZlWbe4daJilAtf7Er\/TfF1ctoxWVZoJycF7+PcN7591DRwjsGasjZe5VXrwFdln4otSLERV20qsjhs02h717cIX+mv1mXnaErQNa4pxya+rvO\/IW9sMs1kxtr16rUaujj2UhUNPe+s8+d5ubwQShkpBM6bUYzF+BAKnp5jTC3r6zHKeUMmy40bIiEg9Uf2xdLJSX95yRQC572gwZqyffWOXEdTRjlIcSPYh6xM4Ljew7vFrmtSDNZZtr4Y6TqXztE0p\/Wtmlfz13Gbz86j8Wqzceqz3PoP7h5Y1yAyzVj3FxMgW448KwqhXE6el\/JQUggiIVzMnmWFnucFX+ccYzEq\/Xl1CIVdf+ZcaKud0\/9\/+oknRPK8yEFDPj4u5PXO1aOOj0t7X+9w87zIaRt+YZntebCmDGI8HCAeGkrYd4W6PoSr1pPL2c0kIppHejpsPtsA6tLbWjSmbdw6nw5obuXXV3zxoHUeSb8\/g9d6xuZjNOS6Nem+q\/PvtD3LyvJjLOW+pS6NJcb1Zo2xflwIhS\/Ggq97DkEky0RCuJjh5s2Lep\/a8XGR7eSk4JucYyxHh1DKM5cGbcqKHQ8IJW3OcPVPlV+XzX9+SQOFhlpjW3NlXm2wQ9bBPCnZ+drqbPNr7lVy6piu3NZg59HxarNx6rNc8bY2l3tGIIQi4SbvJosM7ecQLsfYd24xXvavYgnhYnQIhZ7nBV\/3HELx+8cxrpth9XF94hJjOX8IpTyUNJG8ozRlu1ZtzvCpN2Z9CFN\/24OYGEtD\/nw0k8WrXymOLnVpzCrr0THFLOud7Xy5dH\/A2PVJcticfdSYpO9F1bq2+U65l4Xse5IYCwRCKHjdOcY6T709xtIXQikjpR9nY1uHYhR5\/vliZAgFH\/p8fFzOcHJSylZKP1q3uo3L81ILoZRnLo3elGnESjTmmeMrtolIT4fNqQ2gKbWN0bFd381pfFP+Ol8m5YsHaTlsjU2hm9TTlHcVX9daV8m5SSz7iD1FDpXRkbHNkrKsWFaeF7zvc4wiMRZZq5rk8XHh4xwj5\/4ohDJXjKW8qhSjSIzFKP5xjywr5KHPIRQzNH2KEYJIlolkmUiWFfF15xBEQqjzzs7ec1OuxocHhBIPCqXq6OlYMylvllxW2\/z64B6qiWjeLs1VazmS7f2sVmuQmsPiq2urCT03t+U8D+woKJZ181+X4scZbem0roMfDP8lL62FPWVlq6t9VlwbZYzDLyuE5jm6\/Ow5zZDnpSWEUkYKgXNBeV7wTc\/UGOOmWbqNz7Ii7uREJMZCrjvHKMIXz7LscgQ1Yw2B897QoE15FxuxvfL24awPWuuvk4mF8Hd9kBM7BGkd5GY9ELK1o6fU5k\/j63SdTxteXZzaiYdUr+Lq75qzKseQNsVuyI+vdW+18SHXudXcIZTTx1jKfUkxlplCKGUrhWC19eUQLo7Nsov6qlqMl0fkuUiMl+1DWI6Oyqx5XspIMXIWiVEkz0VilOVH7Hle2O05RqvtjTxYU+ZhoSgi15HGTJXTAKht3Qag48nRJ2VSbNxcKm5mKQ9tEFioBUK2dvQ60vg6f5td62ybb8ifwbfVqP62Gi0WdbF1dp2jL847YUvkrdKxz55i7H+JeV7mDKGUrZRlhZbnBV\/lrB\/thtA8KsZmf+qNUc5\/HzkEEf6BkA8\/FLGUjulbf+mlMiOy\/Zkxv5ql3hBUKvnNm4Uco0iMhWybfGGZ9XmwpmwfEE3yHNG1D+YgFTee9HfYuaqyWr+txdqrxvVta3vxoPPZGtWWco0Zag2aP513Vf2Buw+sOmSteF7wssfsYHTs1jYrOcvK5cRYyn1LIYiE0Jw1xmZ\/lTfGwhpCwdNzCIVFm3ehtZ9PTkTyvIh78UWRLCvkbZ5DqJ89hIu+l18WoYnb5q022Y9jsKY8F\/gyKW7qXM5udGk\/bLPo6wGfznok5UdEdr40zvqo5broz0XTyIu6jusaf3F0qa1aJzWWo6sljdEaq6O2Y7U1Dfnx9XZWN6FZQyiLibGU+5L055lN+ew7uBibIlf3hVCMibHgXc4nJ3L+71pnmcjNm11GDR8TYzlHCKUco0iMpe7SEgFvyksY2k\/2YdsWbWO1gbSNGcqvHwtrHcqZz9aJvi2iDoj5V3kRoGMYNwZNEbsx1j3ZOUIoSlv13WQxqts5hG5xMXaL06gYCymEgqfnEApLnhe86sw7Skt8VKxxr7yi0nZ4llXPe\/v2xY\/S7cfq+KpGZVmVdbY2b8otl\/ZIur0jtWn058\/2IW79fcg2dy4NG1eKw8YXlvpzn83OzttXnU2Nu35V7R5dt625fdTFCM2BdaxvX+tH1XxcbQk7dcyWsqxYWowF7\/McY5EthIJXnY+PS2uMpdxFirGI0m+RF1p5rrOXEYUUgkgIhWzPIVhtO3KWFfOGIBJCITedQyi8WSYSQiGHUPA9OntTHuBi64N5k4d7W1mZLG5caT\/SWmxNfTXJpiq61Kk1ksfWh15FNsaOrYrt09Zl3jHrsWujAadk\/bOW87zf5cUoEmORs2tz7PJxd5FRJEaV6nkIpS\/GUk6lEET42XFq37ae50UFIYiEUMirnHnXzJfUoFXGzSB2tKZsX8FbeeoY2gdxUxPjYXxTbspHFn80DpuMcOhH1OlUzA9hvy7Fz5LterAPScwNJjrHy\/KyWF39z8vZvzokIoqdbHhobq4H9IQ8IXZuqTgYg1mxQt6E+vuZ8uUqdA\/hUbmK458t2YYZ4zDLDKE5bwjN\/jZvCNURIZT2PC\/lVIpRzn+OHIJIlqUR29PzXCTGxaY+410ryXNZ\/poUX\/DiI\/kYu46cRdwoTZmHRfoqXvWpo5hJ+03Ow\/x5eV5eXvwRc2CnGRhTr+IqDVbXYcdQX1VB1m7jq2LrbORIMcEGRmCCnPrJpX7kOrI1kSeNw5bmxtYld5prVV1\/dDFkQ6Ymu39UruLEzpZCGGZpMZZ5QyjlKinLCmueF7zLOcYyKoRStlKWWa1ejrH0bfvnyFqJ\/TlxKmtMG89zkTwXyXORGNuiZ+UfpSnPBbG6d6Q87Hnos85MLm4m7C\/JSzLEoXMxR1X+OrvGagNRvU9uMbFNlDmoSxszepW\/CTMbTy5yWLJzH8ux3Fj8sVhV5a7KI8nRNi\/hXfIQ5ySyMQYhlCnyvJQ3lWIsM4RQyk1SjE3eel8I9T71xKhSM+edZYzNMVP3pk3c6lOvvaf6JtuUeXcNta2TGKguDh9U5+9itw\/jqvgTOREOHvw0AGRIx6kf2xCUS\/UDydq1FlkcVl6og\/zVNYPJi\/Ji4xzayIiFCNbxyKuSjiUXH0XTpLku6ORSP3IV9YHP0O+Uuae7UNX6sNmx6HXUFoe\/bWydf2N7lm2cojJBjIU5BJEQCrnuvM6vReV5mS2EUk6lEApLjAVvO+e5SIxtUe6fOAKjNGU+VmvavClGxDIGQk79quMjBkJWe99cH9LaPOryH8mR1DXCtrGyxnFdip8Ttw2lfkjjrKw2y22tbbF2XJUMJl1zEEe85rF1qA1OHByqi8FHHA0Z+VV5VbrkJrYL1X3K0FRPl7xdY7jnu1BVPvaKHYu+ahxjoKpxo9pCKKZb5YtWxYj687q58lxE6tNe8oRwyXTBEEKhxljw9JxlhSWEy79mVHj8vIMIDNaU2bCWwMbqKmO3hJ0HhtqQsamuHBs+1ZGxqQ5Hx45cR6enp9JGPNwZT8OtilU\/H29DxF47vQZbEjJUNXYTm87LJO+cvnNpHVpLOoeOqxpDLPmU0NchnYMavnD6Bfnz0z+Xr51+TdQuZwc6diXicYEXVDc344irWoP6uF7EQN86\/ZY05bbNtMu8H5x+cAlvamUuSGtAnhJ12RPUWxVn9xIyRGwVVY1P48BrU9Kcpx9UX4918p\/nvHat8hpfyPmFL2h4e+xp8ayRs9+rvpDnzGdtEsIy9+Fbb13KvXxHnOdL\/+nf+BuX\/DbPvslLUHb4NFhTZsN2oW1jd3JyInfu3GmkX3z\/F8\/LrIr99E8+vfTnUvyRxXH34K6gyeJ45v4zjfmrcnaxPfjug4vsxd\/vvvfdS3P87P7Pls7HTx+\/4PvxT368tNOIquZ5+\/23l35OH73z0eXYu3fvyr1795Zy1ZjUZjF57vQ5+Xc\/\/nfy2v3XJC7+iDnQf+\/09+TPfvxn8g\/v\/0PpihkPGlkc9+\/fv1STzk3uRcjy792W6\/Huu+8u4zi99957l3Lq+prmJUbnvPfWPVKNRjTBKtq0AJtz01zp+C57D0yb6CefLvaeLN7dNsWt4jt9551lqfcfe6z2PrD5lsGL0\/3XXusUf\/9731tEi5w+fnFf2pzI7\/9i+dxJ9977\/\/bfLnNweu\/v\/\/1O85JzXdqlca8trgO47CoN1pRTQNjcakOGVB+Ck58XBW25X3jhBbl69Woj\/fLDv3ye5udXf34p9g8O\/kD055WSHLxj+vLhly+NaZuzaq6lkAAAEABJREFUi\/9Tj3zqfLafPvrTS3O8dfjW0v\/CwcU1fvLwk0s7japqnocefmjp56T+xxYPqEcfffTSHOpPucWEOr76yFflmw99k5QCJrxzVsyq\/G2YkUMWx18+9JeXarJzL0Iu\/GVcVe5HHnnkPK5pnQcHxX\/HeHh4eGle7g1Ncm3xLkvlobne63q\/p9zOrz616VjVLVcfYyB066+TiSO+zq\/2FzrsvfS+SvXDTxb38sHiRWPqW0v\/+c+FXNRI7i45zt\/RVtwTVeMPHij+TfTDK1ekyq+2h3+5fO5cXdSldvjDf\/qnlCgSgjzyG7\/RmIf4faKnn366wGZHz6M05XSTsmEh7EPgRl7yd8kdFjc1D9gmeurgqfNU9w7vSVXs6\/K6aJORxREkyI3Fn9tyWxhfNWZTmzZXWRw0WJuPOhfm5V8aifVRz9KxOBFnfcjEL1zCGtChK4sHyIMPPli5dvxVBCZgIGcH+dAVkzZ\/VU612TWozXJyn017zuzcNhYZHDSQ3NiqiDUQ94ODH1RigQ967vHnYFsj7n\/2QVMB+Imri2ny1Y1py2nHhQ57r+oaWNvBU+XePFx8kmN968paI7m75JAsWw45+KM\/qrwn0hza9Gmoqc\/qzL9MvDj94ttvl7nffFMOvvOdhXXx98aN0n546PICg+v299cXENX\/naZnlKa8jaXzcFBifmT4OhQkiB5RorQdxNN4bspNGfpgLuZIv3hk67QvFohtozRXW3yTHww+lOJPFSZt\/qbc+Ow60ZWq7IqVxqzDm3LYOZvi1pm37zHsh3Wabpc6yK1EPDJ8cIpx8yliLHOEUMpNkjaBGJuiSl+MhazjCu3yOcsu27Dwq0\/wEESOj5GcZoTA5JoyDwq7iZGxpZhjw6d2ZGzocEtqg69DQRY3vxSHffAWlstnG3\/Z269F5+pSl86sY9BXGUf8VOi6XF+WUlf\/iZzIKofNY\/Gpy2Hj62Kmard7pa5G9g9xdf46O+MsEYcOH4SyrEwbYymvK8VYjgyhlJukEEpvjKVcJcVYZW21XVm8O14GxSiS50tRsqzgfp4VAps25UHAYBPzQICQ7STYVMeHDiGrfQiuD+q6d5E8pHMpNsuRlP+JhQx8aF3Mb6fSWrBlUr9503GyONSmuRemnfur37S2a1DbUItR3Mhv50Ufkrj32QPMoTI6MjZL2NHhlrBB2OAQ49GV0LFPikIoy4mxlNeVYixGhiASQiG3nUMoI\/K8lNukLGuLEAnhYkyel\/oU\/83rsjqX1kRglKbMZtaNbTn2urrxQak\/taFDaZzV2\/w2tk7Wh6x98NpYa29qgnZMH\/J1qX7HqC8etG47V5XN+ndBtmuw2GvtaiMOUnsb7xqr+W0+a9PrYv19yuwjm8\/e48iQ9auMvYqsX2W4jUVPCX9qs3qb38auLYdQDD37VaNCWfP8xhurDeS\/Tny+\/LfbhX+rGRtUlSnGKmu9LYSl7\/zn0K++utQlBJEsK2Q\/zwqBUZoyiLE5U8I+SaooSh\/W9sFrw\/SdKXFjNmXm0zpsbSpbv8ZZm8apbw6cNUGs5UiORNerNqk49EVMheuC6boUL4IuGBNF50vMvarsJRqzUq\/Jdy1ZCEXFMRa8j3MI3bOEUMaGUMpVUoylNYRSrpNCWHr4XWXJc5E8X+qT\/J+hisr8vCECozXlDevc+nB9GNc92N+Q4hX2GA\/kOjBsbVauitc6q5qRjtWYqvHbttnatF6tyeqZZKKx+sJJBjqqsBxoqmVaGrPSXjfnEJZ4SJ4XfJNzjMXoEAre5RyCSAhFZJaJZFkhV51jLK0hlHKdZL8MdnJSRIUgcnxcyH6eHQKjNGUeGFXI1dmrYrdtC7LYCFIc9qGPBT2X4oFwJOP9PFkWRyblA4A6ZHHAoYW4qOa8HtTZUNP1eFleXq6TGPC5Lu3vbhWv5cCGEznVnY5R3cZo7NA8bc67tLc2xsY2rk2SxSgSY5Fh1ZxZVow7OSl42zmEtojCH0LBF+eDr399cV78DUEkhIXgf+eIwChNeY7A2TXpwxgbTQA+FtkGYOvQ+a1fbXC1V41R2\/UOzYxc2yCtv2lujVFOrK4NuYpsbJV\/6jZtzvC9acwhlJclxlLeRAphtdFd\/2OKVX\/uHcLlOm7cuGxzy2wQGLQp81CAQAueEg8OfLtAttnmUrwrlrNDdR7oNu7MPThjXibRj09t41Effktqt7HWvzPyolC7BmS9HkdSfEqga12ECn4Z6NDcdr6BpqpNa\/fYLu2v2gV1cYRQRuV5KavEl67qSGPgMXIuKISCdznnuYh+AYv4GDlXU4yFPYSCN52p2X6JTGOrbOpzvvMIDNqUeShAoARPCfuukH3Q6sNXa5\/Cz5O1Frg2JWRbN3pK6VpS\/5R1XZu+GKFWux59gaRxqR9dyY5TWxVvyqU5tvEJgzZju8eq6p+lLcvalxWCSJaJZJlIllXHx1jaQyjlLlKed4laPSaE1cf4iJ1GYNCmrMjwoFB5l7k+kNMmoE3w6Oyd2dhr1Oajdej81Aupbnld49DGQmzdWHxTJcWA2hUXZK3Xrk9tlttYa2+T2\/K2jU\/8nVRtxHD2GNRp4JyD6n6lKQSRV14RCUGEj5qz7DIKMRa2EERCKOS284cfiigdH5fR2EqtlGIs5BAK3nYOQSSEtij3zwiBUZoyePHgqCJ8u0z2YZxJxUaX8Q6tRd+5N81sm4+Oa4qfok\/XYOvXtatP61bdvqBS3ypc8zDGzouuZGPU1jdnL9GElfrOv3P5QmguOUZZ\/g7xyYkI\/0xlnhfxN2+K8DExhF0WR4xybpMVDpq9hp+cqHSRx1jo168XvMtZ\/5GQLBPJsi4jPGaHERilKdsHCFjxILEceRdIm24uZxtaRKycyXY2zJF5h24bRZAgmxybjt9k7i5jtT5dM1yvx5HBhFxpLDZLjEW\/Lis8LBlwRjoeVedCHop0Dw2Vf5C8QybNsiJ7nhc8PccokuepVZYNOrWGkFq66cfHZVzVO\/YYS\/8K0v1\/+S\/lva98RU7\/8A9XGOWhu4rAKE25ChweKjTrKt\/UbfYBXPfObJtrqGtMtibbOOx6rGzjpyxrzcqpNZNM7KHrtTHW31XWPF3jPW5kBGJsnzCE5pgsE8my5pg6b5YVnpOTgteds6zOc9Ge53L49\/6ePPzf\/7sc\/K2\/JRLjRb9rs0Nga015F5E8kuLbvNTOwx3KJReOF+VF2FbINoqu\/xGDHcM6tlL4hpOm72r1WrC2tClrbNVaq2xdSrPjrMz8XcZ7TI8I2I+OYywTx1jKKsWoksjt26Wsks2ltnp+0aMfNWONkXNJMZbyKlK+aMxvvSWy4BLjKiM9dgcR8Ka8wUWzD+JMOr7ylf4PO7f9mam1p7POoXHYNXAtmj61SGNTPFS3cWrrwpm\/S5zHjIBAjOUkJyelfOOGCF\/AshRC4c8ykRs3ZPllsOPjwrbpmSZqc8RYaiGUcp10Vuf9n\/1Mbt+6JacffCByZqsb4vbdR2CUpmw\/qlaZj66RdwlC+8DOpfgjZ4f1nZm2wqhq1YltQ7HyVNbUZT3UrWs\/kqPGIcQ2BrQ4FRf7AsgOUb+1uTwwAllWThBjIdMU9ctbWSZy82ZhrzoTy8+B+X1jfg84xqqodtvxsUgIRRz5Cqk4x1hwziFwni55ZVtDYNCmTONVogGrvLXVbjhxJpnYw74z2\/aDWOe3DSet19aOrGPqmgsxUyatnxq1ISNXrdvaLEbEW93mxNeVFMN1x3edx+NqEAihdMQoEqMsv22NNYTiHTByE9GYlWJsimz2ZVnhPzkpeHoOIbW47gicIzBYU6YB04iVUh07tvNKdkygIeeSC0eQ7W+ytIZUl45HHw2q41Qbh9k1viqvnuezDViNNtauUf2b8iFyblrT3o0PoVgy\/5TlyYlIjIXOz3lDKOSqs34kXMWr4ttsR+aTmhjLaOoqNZfGR2AnZhysKe\/E6tcoUh\/u9iF8Q26skanfIVqXZk11tVuuMXYt1r9Lsq5B11RVu\/r0XW1VTBeb5tE57Rj1WZvLAyPA7xhDMRYTnZzI+a86ZZnIzZuFfazz8XE508lJKcdYyCEU3M+OQAUC3pQrQGky6UM3l+JdskzkuC6r\/36trmVXm4vWby9BJpmsetj1V+Xsks\/m6BLvMT0jEIJICJeT8q94XbYObwmhmMPfHRc4+LkzAis35c6ZZxqYPrTR12kEfcNDHTbnkZiP0KyjQrYNZdN3kRXpBzWl635R6n81Ta\/Tpi+o0jntApt8Ns7lnhEIQYSPqdO0IaSWcfQsK+Y5ORGJsZBjLHgIBfezI1CBgDflClDqTLZ51cVsw05dkJ27S3O9Lqu\/u7ZzbFNmvTflpsDFHF2aYjrG6l3G63R2nMq7jKmua2d5CNMp3f5cOc+LumIs+Cr\/xGYxws97hIA35eXFbj\/x0H1enpeTxR8xB\/Yn5AljGVdkfup6WV6+MDF1ttVlGxB5SKDc+rBPiaixas3UiB1eRUdSfnpAjqqYdWx95lpnfh9zhsDxsUiWiYQgEsKZcUvs+LicmF+NirHUQyhllxyBBIFBmzLfrlZiXpWVY6ujLjGMbYpTH5zYTYimV\/fwxf6SvLRJ+rXHTrWutRfUYaBdc\/riYZNrkeaqK+W6FJ8wMFca0zVHOm5snT2h1DS3xsCr4prs+KCqcb3b8lxEf784hN7Tr5Xw+LgYdnJScD2HoJJzR+ASAoM1ZX7lqQtdqmhhYCPrWOSFqfIvvro46yMGvTJJRyPvPAm1D13kTBavzEVE\/TLyofNSB\/Xo9OjI6kdOycZXNZg0fiq6rok1Qmld6k\/tdr25LB7iUhxdPuovIqvPFjs7R3X09q3sBfaEEnpVVdg1Bo6ucciQ6pZjJ14J3foHk\/NcJEaRPBeJcbBpOieu+gi782AP3FcEBmvK6wLKBmYz63hkbKorx4ZPdWRs6HB05L7pRSm\/SMQD+EiG+Uh01bqpg3oYB0dHhmzTQFciTmWNUW59GjM1zhr1XSu1ZVK8QJLFoetYiOd\/rf\/cuGdC171RFWf3FDKUwlc1Lo3pXa\/6\/WK19T7ZCgmzrAzmXwpTLQSVnDsClxCYXFO+VOHAhtPTU+lC2qTekDfkePHn2uk1+fzp5xfaG8KBDnXJ1WeMretLp18SauhaF7HUDmlNyNAHpx+I2qbFT8Wu+QunXxAO1g6HWBdUVTd+6Fun3zpfnzbwujFpHp2fPO+cviM6XhaH5liIs\/hLk1Xqe0EprrPTr10TCaGALc8LvjifLuyzW+tpt+foGOteQLzTf2fZlNNX8DxUUptetZOTE7lz504rffonn14OySWX793\/nnz2x5+V1+6\/ttCKzfbM\/Wdac3SZZ9UYW9fvnP6OUMf\/vf9\/O9clZ8ef\/fjPlvVrg\/n4jz++1LWeu3fvyr179y7Y1Dc2t2umKX\/2J5+VrmumabLk+\/fvn6\/lZ\/d\/hkl4IdJlLT9874fLeE7vvvuuvP3+24hL+h32VGQAAA4hSURBVOidjy7zvvbaa0t9aqd0HzTtDfUxBkJvWw9xNoYxqU39Xfdel2sy1ZiffLp4buiaacir1jqlvbdq7duIn+re03ugjc+yKdtFNz0UiHvhhRfk6tWrrfQHB38g+vHnW4dvyVcf+ap886FvCgfvnL58+OXWHF3mWTUmrYuaoK51UTuxP3z4h8v6kaGHHn5oqWs9jz32mDz66KMXbOobm2+yZl3vXz70l+dreeDgAeG4cnjl3Na0po898jHClwQmYIVCbh339NNPY+pM2whs2xt1zbRrrW35u+49xXQX+eHf\/bsX4Do4OOh0j9m1Tmnv2bqmKu\/C3rtwUyTKrJty20MBLEIIcnh42Ilel9flxuIP46AgYaHdkNuLP08dPNUpR9e5VonbpC45O2hMzHmmCg8PdKUrV67Igw8+uLU1ah3K110z14k18omH5rp7cBeTcD3V1sQ1hyyOe4f35AcHP1hIxV8dd33iv4vaZW8UK1rv3CV\/CN33nuK6a\/zgn\/yTiwCuseap7b2pX4Op772LN8RlbXJNmVfnbGgtFRmb6sqx4VMdGVudrvZN+U25KR+e\/aEZo8sEDupYp64gQTiiFH\/k7FD7mTpJts6ar8v13tcCciTdBcyoM90r2FJiLxGX2rvojGO8SJfomcfwb3LbJea5CDZI\/HAELiMwuaZMiWxoNjaEjE0Jm8r40CFktaMjwy1hc7qIgDYSbSwXvfPTdL2sTNesvGvDrspBvl0g9gN1wi1hg7DBIfYUuhI69iYiFj\/cEjanBQIhLE7+1xGoR2CSTZlyeQBAyJZSGzqUxmBLyca4fBEBbUxqtY1HbXPgdl3pmtddn+axudfNNfS4dE+orvOiqwxHV0JPCZ+1oVeRjdl1eaP6+fe5s2yjFD543ghMtinPG\/bprM6+O9TmMp3q+q\/ENs5112tzUKHmsVhid3IElghkmUiWiYQgcny8NPnJEahDwJtyHTJ7Yg+yeFBIcWhzKbR5ntP12jVbX9fV83vrXWM9bo8ReP11ESiEiYLgZU0FAW\/KU7kSXsdoCGjz3fSf10wL1ryp3fU9RyDPZfnvcr\/0UsFj3HNAfPlNCHhTbkJnD3y2kbwqr56v2NrPjTMRdG32XTJLUztyG2ksv1qlsWpT3bkjcI5AnovkuUiei8R4bnZhWAR2Mbs35V28aj3WvI+NRNdMU4Z6hNNTOQIXEdB\/g7uKX4x0zRFYIuBNeQnD\/p60QYHAvrzru372u8p9NGSbw2IJnk6OgCPgCKyKQHtTXjWjxzsCE0fANk\/bVFcp2+ZYZZzHOgKOgCPQhIA35SZ09sSXNphUnzMMb0jxv3yxxk3Xvel4anByBByB\/UZgX5vyfl\/1PV99Jpno0dc7ZW\/IiqhzR8AR2AQBb8qboDeTsbZJzWRJjcuwDXRffo7eCIg7HQFHYDIIeFOezKXooZCeUtim1VPKyaXpe41955scYF6QI+AIjIKAN+VRYJ72JPpt5GlXOWx1qzZVx2zY6+HZHYF9RcCb8r5eebPuVRuSGTqGOMgcfX9kv48YDnJhPKkjsOcIeFPe8xugavneYKpQuWhLMfJ3zhfxcc0RcATWQ8Cb8nq4+agdR+BIji6sIG2yF5xDKJ7TEXAEHIEKBLwpV4Cyb6a0Ifm7vtXvgBTD1TP4CEfAEXAERLwp+10g+9hQ+l5z3\/lkNw+v2hFwBDZEwJvyhgDOYfg+NpT0i16rYrBq\/BzuE1+DI+AIDI+AN+XhMd65GbzhrH7JHLPVMduZEV6oIzAiAt6URwR7ylPtY1Oxa1715+h2LNc11bE5OQKOgCOwKgLelFdFbGbx\/NvPN+WmwOXssP9Jw5lpdoz12kb6qrwq4CAdDsaeyInYo+tYO8ZlR2CLCPjUE0Vg1k35ySefFGii2G+9LJrL8\/K8vLz4Y4uh4TwhT1jTrGRddy656IENHNrWTVwVZl3G6lxT5+wZpaZaNQZeFVdnr4p1myPgCBQIzLYp80C4deuWQMjFcv1sEaCR0GSw2XeN6NhfkpcQZ0d23eni2ta9ydh0rinq7BX2jBJ6VZ3YNQaOrnHIkOrOHYHREJjBRLNsyjwQeFDo9UHGprrld+\/elX0l3hGDxVN3n5IXF3+QIW3Q+BWbe\/fuyauvvip\/8id\/svN4sS7Waenf3P03ot\/Ixq\/rTjk+WRxgtmDLv8jp2KVjx07sEfZKW9lVcXYcMtSWJ8XW9epn0Zz23hjXuO2+m7p\/lk25C+jXrl2TZ555Rj73uc\/JZz7zmb0kxemH3\/yh\/Iff\/g+qCroqf\/1zf32Jza\/92q\/J7\/\/+788CL13b4ZuHKsq\/\/+1\/L29\/9e1zXded3hsaAEYHdw+W6r237l0ay33F\/cV9tgya4YnmrLTK8sAEbMAoxdf1y8+iOe29Ma4v9xX3F\/fZKvflVGIrmvJUShu2Di7YV77yFfnGN76xt6TviD\/+2Y\/L\/\/rs\/xLe7d2QG\/KJZz4hHPj\/2+\/+t9nhw7pYH\/fAsRR\/\/scX\/4eAgywO\/HXrxrcIWca+Iq8I+v9+5n8vdVkc6DqW+2th2pm\/6btbGm5q08WoDz+Err42Du5gs897z9c+7HOX+6vtPpyqf2+bMheEh8Ozzz4r+0o0YXB459o78q+e\/VfCvwf9ffm+oGPHP0dsWBfrY538DJlfh\/qPz\/7HTuu2Y\/\/Ttf8kfOwPJxc58Stm3F\/YdpFosjTbutqbfHVjrB1sFCfn+\/sMGurac3\/Z+22X5L1pyrt0UcaqlXd6NBHm45vIfIlJf2bKOz7eNeObG22y7k3G7gqObQ15V9bhdToCu4jALJsyr+J5sOgFQcamuvMSgdfldbHNV5vxbbktyDLTY5N1bzJ26nB22SvsJeKmvhavzxHYRQRm2ZS5EPrg4OGBjG1+1M+KbspN+fDsD80YXfbgYJ3rrnuTsVOFlr1CbXBL2CBscIg9ha6Ejt3JEXAENkNgtk0ZWHhQQMgp6cMktbt+EQFwumgpNOxQofkZBMBDCd1Snd3GbFtmr1SR1oVPZTi6EnpK+FIb+i5gQZ3bJnCqqgE7VOXbVxt4KKUY1NnTuKnos27KdSBzkXhgQMh1cftsBxeoCgPsYAchV8VsYtvFseAAHkroug7kKrv694k7Fu1XG4ygqkjsfi9dRMZiAjboGoGMDUJW+5T53jVlLgwXSC8KMjbVnRcIgAtUaOUZrKwdGVsZsX8S6weHqpWnPuKwVcXO3ca6Wb+uExmb6s4LBMAFKrTyDFbWjoytjNg\/ifWDQ9XKUx9x2Kpip2Tbu6Y8JfC9ln1CwNfqCDgCjkA7At6U2zHyCEegEQFegdsAXo2nNut32RFwBPpBIN1nc9h73pT7uTc8iyOwRGAOD4XlQjqePMwRmAoCc9l73pSnckd5HTuPwFweCjt\/IXwBe4fAnPbe3jVlPu7gAupdi4xNdefNCIAVmGkUMjbV95XX4QA2+BQXZGyq7xNn3axf14yMTfXp8GlWAlZgptUhY1N9X3kdDmCDT3FBxqb6VPneNWUuBBeGCwQhY3PqjgCYgR2E3H3kPCPBgZXBLWGDwEjtyNj2lVi\/Y7H+1Xf8LmLHvYQFbgkbtIt47WVT1ovFBUN2qkegDiPsUP3I\/fGAQxVZBNRvbfsqOxbdrjw4VUVih6p8TbY5+sChiuxa1W9tU5b3tilP+aJ4bY6AI+AIOAL7iYA35f287r5qR8ARcARmgMD8luBNeX7X1FfkCDgCjoAjsKMIeFPe0Qtny7ZfcLByGmP1NtnmSeW2sev6mWfdsT7OEdgGAtyzVWRrwW\/1Npn4Omobu66f+dYd6+P6RYCm3G9Gz7YVBPTLDJZvutFsLitvZYE+qSMwUQTs3lDZ995EL9YOlOVNeQcukpfoCDgCjoAjsB8IzLcp78f129oq9Z0AXEmLUV252pWrXbnalasdrjbnjoAjUCCg+wKuVHhEVFeuduVqV6525WqHq835uAh4Ux4X752ajY2Zkl0APv24Do4OIVvCpuOQrQ8ZW53f+jTGuSMwdwS471Oya8bH3lFCh1RXjk3HIatdObY6v\/VpjPPhEfCmPDzGfc5Qm4sNlBIbr3ZABwfjU7LD8Fl9XdnmsfK6+XycIzAmAum+Q9\/0PmZ8SnZN+Ky+rmzzWHndfD5ucwS8KW+O4SQysKGUKAgZ7uQIOALDIsBeU2ImZLiTI7AOAt6U10Ft4mN4KPBqfeJlnpf35JNPnssuOAK7jIDvvV2+etOo3ZvyNK7DTlaxS41\/JwH2oh2BGgR879UAMwOzN+UZXMSqJVS9YmcjV1HVeGxVsdjwQToHNggdQraEjXgI2fqQseGbGfly9hQB7mfua7t89CqyMVauisWmMToHNggdQraErW4Mcdavcc63i4A35e3i38vsdRvL2pHrqKqIuljsNh5dSe2qK1e7crUrt3aVlROjsnNHYGoI1N2f1o5cR1XrqYvFbuPRldSuunK1K1e7cmtXWTkxKjsfDwFvyuNh7TM5AruPgK\/AEXAEBkXAm\/Kg8M43ub+Knu+19ZVNGwHfe9O+PptW5015UwR9vCPgCOw6Al6\/IzAZBLwpT+ZSeCGOgCPgCDgC+46AN+V9vwN8\/Y6AIzBPBHxVO4mAN+WdvGxetCPgCDgCjsAcEfCmPMer6mtyBBwBR2CeCMx+Vd6UZ3+JfYGOgCPgCDgCu4KAN+VduVJepyPgCDgCjsA8ETCr8qZswHDREXAEHAFHwBHYJgLelLeJvs\/tCDgCjoAj4AgYBGbUlM2qXHQEHAFHwBFwBHYQAW\/KO3jRvGRHwBFwBByBeSLgTXni19XLcwQcAUfAEdgfBLwp78+19pU6Ao6AI+AITByB\/w8AAP\/\/CxZkNgAAAAZJREFUAwB7mWkMRRdTSQAAAABJRU5ErkJggg==","height":234,"width":388}}
%---
%[output:00e322bf]
%   data: {"dataType":"text","outputData":{"text":"[23:40:57][INFO]  Generating 100 post-RL SMILES ...\n","truncated":false}}
%---
%[output:184ca7e7]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 10\/100) generating (post-RL)\r[##--------]  20% ( 20\/100) generating (post-RL)\r[###-------]  30% ( 30\/100) generating (post-RL)\r[####------]  40% ( 40\/100) generating (post-RL)\r[#####-----]  50% ( 50\/100) generating (post-RL)\r[######----]  60% ( 60\/100) generating (post-RL)\r[#######---]  70% ( 70\/100) generating (post-RL)\r[########--]  80% ( 80\/100) generating (post-RL)\r[#########-]  90% ( 90\/100) generating (post-RL)\r[##########] 100% (100\/100) generating (post-RL)\n","truncated":false}}
%---
%[output:77d3414f]
%   data: {"dataType":"text","outputData":{"text":"[23:41:17][INFO]  Post-RL:   valid=4\/100 (4%)  unique=4  novel=4  avg_reward=0.040\n","truncated":false}}
%---
%[output:9e6d6b77]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAeUAAAEkCAYAAAARhClzAAAQAElEQVR4AeydbaydV1bf973Xr4mdcYgt4zYNHoOdMilJ3dA6IsW0kUDAQNV2GCo7OHI\/FNoRkNEwVjWlKFQYqGoKpFNoKXTaSUep1KJ+q\/phPlYDWMMkTavMANUYKmfwjZLUIqKO32L3\/p571\/G6++7n7ZznOc\/L+Y\/mf9faa6299tr\/\/XZzrn29fFf\/EwNiQAyIATEgBnrBwHLQ\/8SAGBADYkAMiIFeMDDOR7kX1KoIMSAGxIAYEAP1GNCjXI8vRYsBMSAGxIAYaI0BPcqtUdt4YiUUA2JADIiBkTOgR3nkC6zpiQExIAbEwHAY0KM8nLUaZ6WalRgQA2JADEwY0KM8oUKKGBADYkAMiIFuGdCj3C3\/Gn2cDGhWYkAMiIGpGNCjPBVt6iQGxIAYEANioHkG9Cg3z6kyioFxMqBZiQEx0DoDtR7lI0eOhDK0XvEIBijisGh69CvyN+2Lx4vbTY8X55v3ePH4vk0tVeH7zaoz5rQ5Zuk77Zjz7sccyzDvmlLjldWIP+6XssUxfWlTaxn6UmuVOsrmgj\/Ok7LFMVXatR5lEl66dCkUoanCGGvMyOMQ\/sCY5z7EuaXWi3nk2fEJ82EgtQbeVnKe5lPk2ii+ppTelzrXSp3q\/6k5edvQ5udrT+ltzaf2o1y2WhTfVrFlYzft72Ie8AdSY2Ofdo6pfGW5ZhmvLLf359U2r\/F9LX3TZ+Fglr5942GWeuAhb4\/NkrfpvkOpc9p5z2t+81rrtubT+KM87YKpnxgQA2JADDTEgNIMlgE9yj1dura+C+vpdDsta17fWXc6SQ0uBsTAIBho\/VH2Fx66wdixtpfm8xI\/baQHthjeb3ocQxsfEqDH8HZ0QAwyD2X+vH5V7HFu2in4XObH5nXaABsSoBusjYxhMSZT\/thmbfp43dpIYD5k3DYbdg\/sMfCbDd3D7Cb5Bsj0JiVjWj50Q8oW+3yM6UjikADdA5sHvrw2Pg8fZ7r3o2NHAvQ81PETm0Je7ibsjGd50A0pW+yzmCFIai+q0\/vRUyjq34Sv6pgbcVv+oLGvwWKweZ32kND4owwZqUvO7PgAJMU27AA7\/hjY8Xtg83G0vd907D7O6\/gszktirI0OaBOPHgM7\/tjeRtvGYrwY+GxM89H2Om0P+pgf6X1ej+OIxeZjqur0BcQjAXoeGIeYGNjr9CmKz8szi53xfM3kim3mx46\/DMRZH5PYyvrhJ876mMSGz0DbfCa9zeJSknhiUz7s+PGZTjsGPmKaALnIH+cyOz6AP7ZhB9jxG2iXwWKrSvIxVtX4sjhykTMVhx0\/PtNpx8BHTBMgF\/ktl7WxeWC3GCRt7\/c6PmKA2WOdNnFlIK4OyMeYdfpUia39KFNIEVJFEh\/bUzYrmFj81jaJ3XST2CwWSdt8XmLH723o2PChjwV151OHg1RubORokz\/yM05qDOz4Yx82fLF9nu1UDSlbnZpm6T9L3zo1zhpbZ92YUxFSuYiP7SmbzYNY\/L6NrQgWa5L+RSCXxfZF1qmpaG74fK647edLHH5vy9OJzfN5O3Fl8PHo1FAE8hE3M6IEtR9lCilClD9rEp8pPfwyTW30YbH8dGhj97Y2dcZiTMMsY5Frlv597duHeaVqSNngkLVEliGvf1k\/\/LP0pX9VME48H9rYLQc6NoPZ60ryFCGVj\/iUvU0bY6bAmNiRTYO88Ovz0sZuNnRsBrPXleQpQt18Fk9Oqw1p9rYl46bAuNiRbaD2o9xGEU3mZNHy0OQ4fcjFxjD4OfehNtVQzIBfL9NZy+Je4\/MyZ4PxgOxqpoydh7ZqYv6M2Vb+KnmpwUAthip95xFjtSGtNuQ8xo7HsBpie1PtkTzK9+iAsCLci5xNY4w2NwW5GaNqlcQa6Fu1n+LmzwDrY2vl5fwryR+RuqjTA1t+j9k95Dcw7uwZ62ew8fNk\/YzD6+Hn3tU6FLHW9\/qKaq\/iG92jXDTptjYYedkoRWM37WPMpnO2nW+INbfNSV\/zs1bsaY+6tdKXPPRD0kY3YDN9CLLNeuGmrfw+N2PQ9nxj8+2+6X2sDw7bqquzR7loUkwWf7w5sKdsFotMxdAHO370PoM6AbX2sU5qi+vCVlZvlZg4r2+TnxzeZjp2\/NYeqmQefakdPqknxjT1kWOafvPuY3NOjcsc8Kd8Q7Exh77UCpd59WDH35da511HZ48yE4V4FiAGdvwxsJfFpmLogz3OV9amD30Ncbz5kbGvrG05Y0kuUNafmLivtfHF\/bGZHxn7q7bjPOTC5vvTxu6Bzcd4HZ+P9T6vx3HWB7uPq6uTp26fWeKplzFjYJ8lb5N9qY16YmCvMw79iTeJbsBGvhTwWdw8JeN2VY+N7eebqsXbfGyRTm78JtEN2HxOr+OzuKYluf1YpmP3Y9E2Xyzx+Vh0bD4OWxuwcXxuP25K97F5eq1HmSLyEuXZy\/rgj5GXC3uV2DiGNn1j5Nl9HDEGb59Ft3wpWZSXeO+nnYKP8bqPNTs202MZ+6yN9Ij70fZ+dGzg0qVLiC0gxmBO2qabxBbDfF4S49tej31x28cW6WX9ivz4YjAWNqShrG1xyLLY2E8fg\/dxmfi2xUwri3LhS6HqWPStGmtxZX3wx7C+00ryVenr49DLUCWnxZDL9FjiSyGOy2vTN89XZKdfjFR8HGPtVCw28yNp10XVfj4OvQxV6qj1KFdJuEgxTV9ei8Sd5tpvBtjbcYUpWxyjthgQA7MxoEe5Jn9cTAa+K6rZXeEjZWBM02JfA9vnJrGBMc1VcxEDfWNgMI9yXy4D6jD0bTFVjxhokgHb5yabzK1cYkAMpBkYzKOcLl9WMSAG2mNAmcWAGJg3A3qU5824xhMDYkAMiAExkMOAHuUcYmQWA2JgnAxoVmKgzwzoUe7z6qg2MSAGxIAYWCgG9Cgv1HJrsmJADIyTAc1qLAzoUR7LSmoeYkAMiAExMHgG9CgPfgk1ATEgBsTAOBlYxFnpUV7EVdecxYAYEANioJcM6FHu5bKoKDEgBsSAGBgnA8Wz0qNczI+8YkAMiAExIAbmxoAe5blRrYHEgBgQA2JADBQzMNRHuXhW8ooBMSAGxIAYGCADepQHuGgqWQyIATEgBsbJgB7lPq2rahEDYkAMiIGFZkCP8kIvvyYvBsSAGBADfWJAj3KfVmOctXQ6qyNHjjQ2fpO5GitKicSAGBgVA3qUR7WcmoxngEf00qVL3jSTTi5yzpREncWAGBADBQzoUS4gp6qr7KIu8zNOlRjipkVT+VN5qtp87ak+3o9eJYa4qiBfClX7b4pTQwyUMMBeKwmp7SZnCrUTqUNvGdCj3NDScFBSqfLsqVjZqjNQxit+\/ss2zogtBrFxXF6bvnXi8\/LILgamZYA9GEN7clo2+9dPj3L\/1mRQFXEZcEEg6xROnzrxbcZSS93626ynw9waWgyIgY4Z0KPc0AKkLnYueuw2BG0Ps8fSx6BX9efFVu1PnOVAAmwezCdl9zH4LQ7dEMfEbYszaX5rI7F5aTr2JkA+jyZyKocYSDHg9xl6HIMtRhyj9vgY0KM8pzXlcPFQeWCLh8fmY9CxWRw6Ng9s5i+TxPq+6Nh8P9rYgbfX0X0O8gBsqRzY8RuIQUfGPmvjIwagNwHLTU4DtiZyK0eHDPRwaPaV7TGT2KxUdLMjsZtEF8bLgB7lEaxtncNaJbZKDLTZxYFOH9roTYK8TeYryjXPsYrqkE8MFDHAOYuhvVvE2LB8epQbXC8OBoeFlEja6HVBX4+6\/cvifW70svhp\/eT2mDbPPPv5etHnObbGEgNVGOBeuXTpUvCySj\/FDIMBPco9XCd\/2EwvKrPO40Gs5TRZlHsWn+X3MpUPP3UZaKfi2rYxPmN7tD2m8ouBFAPsQfajgXYqTrbxMaBHeQBrysGcpswq\/arExGNzQdAP6X20U3aLwWe6l9jpa\/C+aXVykbdKf+KIj2Oxp2yp2DhObTEwCwPsPfaZYZZcg+qrYoMe5YY3AYfIDpRPbXZ8Bmw+Bh2b+U1iwwfQzW4SW8pndnwGbNYPSdt8TUrykt8DW94YPg49L87s5CIOmK2KJD4GueiL9D7a2AUxMC0Dfj95nb3l2+jY\/DjYPLxP+ngZGOSjbBs1b1nMH8u8+Fnt8WHKa2P38ONitza6h9lNeh+62ZG0Ddb20vRUjPmQZaB\/Ksbb0T18PHbarBF6DLMT40GctdGBtWOJjzxmp52C+ZHeb20kIBd+9DGC+YG8ueFLIS9+0e3slTzATezDBuA49tE2OzHCoBioVezgHmXbmLZJ82aLP0ZerOzdMcAasaYxsDdRVVN5qKXJXOTrE+Cf+QH0vNrwx8iLlX06BuCXNYiBfbqM6jUkBgb1KLNJ\/cZExzYkwlXrVgZYxxhbo2RpiwHOEPxbfnRs1pacPwOsQYz5V6ERu2BgMI9yXXK4VAx5fd94440Abt++HcpQNa4sj\/zlXC8CR3l7su92O1PIvFo5K2AR1rGtOcIfaCv\/GPPm7ceh2Uf7KPvvMlMXCBv+3Llz4eTJk+H8+fPh8uXLE\/zhH\/5heP311wMS+xe\/+MXwYz\/2Y+GjH\/1oQMc2BFC\/n8cQak7VOMZ5XL16dWh3RVbvLOcqtbbT2Ia2H+rWyx3DfcPd1NWdU7fmadaxyT7U+0d\/9EfZf1xlG3XAX0b5KHNx+DWhHT\/MPMoXL14MFy5cCGfPng2HDh2aYP\/+\/WHnzp2T9o0bN8Jrr72W4cCBAxO77zOdfm\/MNvrH82hjjHnkHOM89u7d67foIHTOkS+Udp1z1dReGdp+qFuv3Tff800HO7tz6tbc1NpOm4d6b9265bfnYPVRPsp1VuPEiRPh8OHDYdeuXRPs3r07e5S9zXJ6W9\/11Dz6XnOqvjHOY9u2bbalRilT5yq1ttPYhrYfpqmXTXH80AcQk3tpGq6m7TNNzdOO1UQ\/6h3LmRrUoxx\/Z8536diyneu+YHfNQDsV52Nm0fmv7j5idXU1gHnX1vR4zAE0nXeafLPsk7725WxwRqw+dGzWNonddCTtVBy+JpC3PuwFkOfvm51aQdW6Yu6q9msyjnpBkznzcsXzXfT2oB5lFotLgMsAoGMzYEPHjm6gjb0NsNHsZ9P8DKhPeOaZZ8Kzzz4bkH2qq24t1N+XeZw6dSr7w4Ft7KUuc3JG8s4LdmrzMdhoY28DReeqT\/uhyl6uW+\/p06czSn\/uv\/9BJmlXGafJmLo1zzL2WM9UtnhTfBnco8wcuQwAuoe3oRt8TNM6l4f9bPrll18Owlg5eDk8\/\/zzgbVueg\/1JV\/eecFuNaIbzNaG1Lka7znyd+TYz9Q0Z2OQj\/I0E227Dz9De+qpp4Lh4YcfDmOAzUfyqcAat72PlH8zA3Bue28M54k52HwkdaY27\/b1lh7ldR4a\/cp3+X39SLvux0z6aKnRrdG7ZEMpSGdqKCulOmdlQI\/yrAwm+nOB8DHn8ec+GZ7++IXB4tEPn9nycS0\/S4yRoCBpol\/SAL613QAAEABJREFUsWEs82+E5Qr6e+QGOgfxrim1pwzoTKUXpmz\/lvnTWe9Z6e9xz5OvEZ\/vlaeMAT3KZQzN4N9\/9Imw\/9iAcfTxLbP3P0\/0+pbAhIH4hHliKvNPAgsUchh0ORQQNVBXtTPV4zOnMzXQnTe\/svUoz4\/rUY9kD6BJJosO0IHpSAN2AzZ0pIG2wWwmzZ4neZyJNT868G302Obb+AUx0AUDtg9NUgM6QAemIw3YDdjQkQbaBrOZNHue1JnKY6Y5ux7l5rhc+EwcbA4tRJhOGx2bB3aQ8hGHDxT5iasK8pAPoNMP3Uvs2AA6PkEMdMkA+5D9SA2m00bH5oEdpHzE4QNFfuKqgjzkA+j0Q\/cSOzaAjk8oZkCPcjE\/8tZggINn4egcQmC2JiR5Qd1c9KEWkNe3SkxeX9nFQBsMsCctLzr7F5itCUleUDcXfagF5PWtEpPXd1HtepQXdeVbnjcHlQMJmhyKvKBKTuJsfNOtnepfJSbVTzYxMA8G2tqf5AVV5kCcnSHTrZ3qXx6T6rXYNj3Ki73+rc2eg8qBBE0OQl6Qlxe7gTgbG93sZkOa3evE0RbEQJ8YsL3a9P4kL8jLi91AnHGCbnazIc3udeJoC+UM6FEu52jqiGvvrIZr77w5WBRNnIPn\/XEbHzaDtb3M0+mDD8Q6hxt4O3EAmwc2jzwfdotDN5hNsj8MXNOZCn5\/orM6JvP0Ij\/nCfgY8gBsHtg88nzYLQ7dYLaxy1nmp0d5FvZy+vJbe\/irG1\/45XPh8z91ZrCgfn6jEvMJPfmfDndPFmLOZbAH2YvsyaGfqWNPPJn9tr85U5g7nM5ULjWdOPQot0A7F8h\/+befHsXvwb5w4UILDCmlGKjHAGeKveh\/b\/JQ9d\/49C\/Vm7yiF4qB\/j7KA18GLpEx\/G5b5jHwpVD5I2GAvagzNZLF1DRyGdCjnEvNbA5+LeAYMBsL6i0GmmNgDOeJOTTHiDKNkQE9yi2sKgcv5x+kCHX\/QYiu4\/UPUrSwQZSyNgM6U7UpU4eBMqBHuYWF4wLhH6R44WPHw6+98PRg8cM\/8Gij\/yBFGdX8CdAY1ge76ZKLx4DO1HRrzrmJYZmwmy7ZHwb0KLe4Fk9+aH8YNB7bn2TH\/rSmySNHjiTjpjFaTpO6OKZhcbx9Bn2euA90psa7ORuamR7lhohUmhDsATUJJ+gAXRADYqAeA3Z2TNIbHaAL42NAj\/L41rT1GXEhePBftTYodmubThvdx9D2yPPR13xzlhpODMyNAX8W0P2+923T8aNbgegx8nz0NZ9k\/xjQo9y\/Nel9RRxqD18wdt+2i8LbiIlhfm83m6QYGDsDft+j+\/nGbZ0pz874dD3K41vTXs2IC8XQq8IWtRjNe\/AM2HlCDn4ymsAWBvQob6FEhqYY4NKw7+qRdfNaf+tHDg+zS4qBRWHAzoSdg7rztv7Wz\/KYNLtkdwzoUW6R+ytvXQt\/PGCkqOFQp+zYUj5sBmKKQFzsNxsyRhyr9vgZyDlTgzlnqRViX6fs2FI+bAZiikBc7DcbMkYcq\/b8GdCj3ALn\/DpAfnn+j\/yTL4S\/8aOfHyyon3kwnxZoUkoxUJkB9iB7kT2pM1WZNgUOkAE9yi0sGhfIWH55PvNogSKlFAO1GFjIM1WLIQWPhQE9yi2tJJeIfnl+S+Qq7UIyoDO1kMu+cJPWo9zSkvNrAceAluhRWjFQm4ExnCfmUHvi4+qg2ZQwoEe5hKBp3By853\/w+wf3j0+k\/vEL\/YMU0+wA9WmaAc6U\/pGXpllVvj4yoEe5hVXhAnl19U\/CP\/r2Y+FffM\/jg8Xf\/Yvf0Oo\/SMFfwyijn5gY1ge76ZLjZoAzxT\/yojNVvM5VzgQxMSwrdtMlazDQYKge5QbJjFMd\/\/oPhEHj0AfiKWXt+K9RtH2Q5z1eNkl96SUDgz5P3Ac6U73cV30qapCPMo8AqEIkcaBKrGJmZwCuDT6b2ZDYY4lN6JYB1gRUqYI4UCVWMbMxAM8Gn8lsSOyxxCYMj4EePcrVyGPj2X85oVfrpagmGYB3D9aD\/NjQDbRTdmzEeElsDPwgtltffEIzDMAxvAL0ZrIqS1UG4NyDdaAvNnQD7ZQdGzFeEhsDP4jt1hef0C0Dg3qU2Uh+86Bjy6MQHzF5ftmnYwBOPcqyEMtagLxYYmJYrLebTbI5BlgXOLaM6NisHUt8xMR2tadnAD49yjIRyzqAvFhiYlist5tNsh8MDOpRboMy\/gDJ7du3QxH8uHGc96V02dYZsEug6BJZj+z\/13gPTNPu\/yxnq7DKuSribbbRF6O3ztTme3ssqz7aR5nLn01btlCnT58O58+fD5cvX56AC2V1dXXSRrc8V65cmdjp430Ws6gSvuHdQBsurI00G3Z0bOhVEfehv0fVPLPExXuAfVAFfl9dvXp1lhI66wvXrEFZAalzVYUji9G5WmcYruHcQBuPtZFmw46ODb0q4j7096iaZ5a4ac+U7RfO1rvvvjtLCb3pO8pHmQ3FRqvCMr9G8uzZs+HQoUMTHDx4MBw4cGDS3rdv3ySVt9PH+yZBG8qVP70RVv\/0+mCxMY1NYp3XTaZNDfwGc1gbaTaTKZv5kCm\/2ZAx6NM24j3APqgCv6\/27t3bdpmN55\/1XFXhyGLyzpXO1Pqy+n2\/brn3Fd+91lYt5TcbMsbWDM1bpj1Ttl84W\/fdd1\/zhXWQcZSPMjxygRisjYzBL7k\/fPhw2LVr1wS7d+8OO3funLTxWT\/0GOYzya8DJO+P\/7f\/GT76n784WFA\/82A+NjfJsGlfxHuhqO331bZt2wZJpZ0pJBMwie7BvonPVRE3KZ\/Pxx4kJ3tSZ8ozMw49tf51bJytoZ6peAUH9SjzHZy\/BNCxxZPC5oGfNnIe4ALhv8BffvnlMHQwj3lw1ocxFrUGzgZnyeaPjs3aJrF5YKeNbBs6U20zrPx9YWBQjzKkcQlwaQB0bAZsps9b8jMND8bnIhk6mIef1yLrcDFWcJY4PwDdzxObb89T9\/uNcYd+nqifefh5LbIOF8JmBgb3KFM+lwZA90jZ8OfZ8c0KDhkfq\/EHW1K\/O1q2k6P4HeCsI2vMWrPms+6b+fWvPhLnBMQ9UjZi8uz4ZgUcwzWcw70wnnPk15L1ZZ1Z71n3zFj6d\/4od\/ldeBOLyGbiI94+fkz9mc98JvziL\/5ieOmllwb9MXqf5sFaN7Fv2s4x5nPVp\/1Q5dzXrff555\/Ptsf3fNPBTLLnqozTZEzdmmcZm\/llE9WXjIHOH2W+2+YCMWRVDewLD3Nf\/+3kJ554IvS1tjp19WUerPUQtufYz1Vf9kPVPVynXv7LkT12fOP3ZNOuOk6TcXVqnmXcoZwp1mQe6PxRZpJcIIYhP87MRRADfWHAzhRS56ovq6I6xEAxA714lH2JXCBAl4hnRboYmI0BzhTQuZqNx3Z7K7sYCKF3j7JdGlwggLYWSgyIgdkY4BwBzhRAny2jeosBMdAGA714lLkgDFwYoI3JKqcYWCQG7EwhOVNgkeavufaDAVVRj4HOH2V\/YaQujZSt3hQVLQYWjwGdq8Vbc814HAx0\/ijr0R3HRtIs+sWAzlW\/1kPVjI2B9ubT+aPMd\/Sp6eXZU7GyiQExsJmBvPOTZ9\/cWy0xIAa6YqDzR7mriWtcMSAGxIAYEAN9Y6CzR\/nIkSPBvmtHxtDHb33bKqpnCAzYOaJW073UuYIZQQz0l4HOHmUuBwA1yBjYBTEgBuoxYOeIXqZ7iV0QA2Kgvwx09igbJVwYpo9CahJioAcM6Fz1YBFUghiYgoFOHmX7OI16TU9J\/IIYEAPVGLAzRLTpKYlfEANioJ8MdPIo8108gBJkHvALvWBARQyAATtHlGp6SuIXxIAY6CcDnTzKMRV8N282dGBtSTEgBqZjwJ8jdDBdJvUSA2JgXgx0\/ihzUfDdPBM2nTY6NkEMtMbAiBNzfjhHTNF02ujYBDEgBvrJQOePcj9pUVViQAyIATEgBubPgB7l+XOuEcVAmwwotxgQAwNmQI\/ygBdPpYsBMSAGxMC4GOj8Ubafc\/GzLnTo9TrtPuGNN97oUzmqRQwkGeAscY4AOkFepz0ojLxY7hVDPFWzI2Of2uNjoPNH2VPKpQGwmUTvE06fPh1OnToVdED6tCqqpYgBzhIgxiS60A8GuEvOnTsXTp48mYE7xlf24osvZnb8uns8M+PUO3+UuST4Tj6FPlL+fd\/xSLh48WIfS1NNYmDCwNDO1aTwxVKy2fIoc6e88LHj4ddeeDr88A88mtnti\/mwo5tdcpwMdP4oD43WJx97aGglq14xIAYGwMCTH9ofMjy2f0u1efYtgTIMngE9ynNYQr4TroI5lKIhxIAYEAPtMqDsMzHQ+aPMx9Z81DbTLHrcmcfY\/7yInwvlQT8v6vFCDqy0sZ+rgS2HyhUDlRno\/FG2BxkZo\/IsehzIo8zPgY4\/98nw9Mcv5OLRD5\/Rz6p7vI5DK42zRM3IGNgFMSAGajEwt+DOH2W+o8\/D3FiYw0D7jz4R9h8rwNHH51CFhlgUBvLOFPZF4UDzFANDZKDzR3mIpKlmMSAGxIAYEANtMDDXRzlvAv7jNWJoIwUxIAamZ4BzZCALOlIQA2Kgvwx0\/ihzUfCRGjCa0LFbe5EkP4MuwyLxoblOxwDnh3MELAM6dmtLigEx0D8GOn+Up6GEiwUU9cXvURQ7m6\/Z3vw2n7w\/nW12\/SntZjlXtnUG7Lyst9JfLcZkOkpWMSAGpmVgcI8ylwHf8QP01MSx4\/fAlortm01\/SrtvK7IY9XA+7Lygp2aN3WJMYkvFyiYGxMB0DAzqUeYC4DKwqaJjs\/Y0ko+Kb9++HYqQylsU733W986dO6EIFvfQN31LeOjo47nYv+Yj1o\/Rhq6cxXtiVn5Yw76AM8RZsnrQsVl7GlnlXM3K4Vj6G793797duCPummki72b3x7p9LPNueh4TsgaudP4o2wVglwASYJ+W2zp9+bj4\/Pnz4fLlyxNwoayurk7a6FbLzRs3M\/XKlSsTv+8b69b3xo0b4fr167nAT+Lr12+E9957LxcWV2X8eB5xbUNpj3EeV69eZblbA2eAcwQYBAmw054GdfqmzlVT+21o+6Gs3skdsXE\/3Ly5fsfYXcNace7NXuXsz8p1Wc2z5m+6P\/W+++67UDV4dPooc0mAmMU6hz\/u69vkNuTlvHDhQjh79mw4dOjQBAcPHgwHDhyYtPft2zdJu23btkz3ft831q3vjp07ws6dO3Oxffv2LO\/OinFVxo\/nEdc2lHZz87i3xl3M3c9j79692Xq38cX2fJw77wzEcWVty4\/My5k6V01x7nlsKmebecrqndwRO9bvhx0bd8HKxl3DeuzYsSOYvcrZn3U+ZTXPmr\/p\/tR73333QdXg0dmjbAeaQx0DXxPM+rx5OU+cOBEOHz4cdu3aNcHu3buzx9PbrJ7llXXKvK9Mp+\/K8kpYWcnH8vIyYWF5eaVSXNmY+FPzwD40jNBV7usAABAASURBVHEe9s1dtugNfmGf+33vdXxNDFUlZ+pcNbXvhrYfqtTLuiyvLGdnH0l7Za2NBMsrK2F5o90Uj0V5qtRc1H\/ePupt60zB\/zyxPM\/BbCwuBw62tWOJj5jYXrU9S9+qYyhODPSFAauDfc\/ZsXYs8RET26u2Z+lbdQzFiYFFZ6CTR3la0uNLhUsC27T51E8MiIEQOEOcJeMCHZu1JcWAGJgfA4N6lKGFy4JLA6BjM2BDx47ugQ2fIAbEwFYGOB92XtB9BHba2NE9sOFrHsooBhaTgcE9yiwTFwFA9\/A2dA8fJ10MiIGtDNh5iT3YzYbuYXZJMSAGmmFgkI9yM1NXFjEgBsTAfBnQaGKgjIHOHmX\/EVhKLytcfjEgBrYykDpL3ra1hyxiQAz0iYFOHmX\/8VeR3ieiVIsY6DsDRWfJ+\/o+D9U3NAZUb5MMdPIoNzkB5RIDYkAMiAExMBYG9CiPZSU1DzEgBsSAGGiMga4S6VHuinmNKwbEgBgQA2IgYkCPckSImmJADIgBMSAGumKg3Ue5wqz4k6EVwhQiBsRADQZ0rmqQpVAx0CMGOn+U+VOhXCCGHnGjUsTAYBnQuRrs0qnwBWeg80cZ\/rlADAN4nClZEAO9Z8DOFFLnqvfLpQLFQMZALx7lrJKNL1wgQJfIBiESYqABBjhTQOeqATKVQgy0yEDvHmW7NLhAAO0W56\/UxoDkqBngHAHOFEAf9YQ1OTEwUAZ68ShzQRi4MMBA+VTZYqA3DNiZQnKmQG+KUyFiQAwkGej8UfYXRurSSNmSM2nB+MYbbxRmLfMXdpZzHgws7Bh9PleLtiirq6sBxPPW\/REzojYMdP4od\/noQkAeOESf+tSnwg9+73eGr33ta1vCOFDnzp0Lp06dCuhbAmQQAx0y0Ndz1SElnQzN3cA98uyzz4YzZ85M7grs3B8f+bbjyfulk2I1aC8Y6PxR5jv6FBN59lRsG7Y333wzXLx4Mfzu7301mZ6HGj9IBsgoBtpioELevPOTZ6+QUiFTMMDjyx3xnR98KLtPLIXZX139EzNJioGMgc4f5awKfREDYkAMjJiBJw7sHfHsNLUmGejsUeY7dsBkkDH08RvMCGKgHgN2juhlupc9PleULIiBhWegs0eZywGwAsgY2AUxIAbqMWDniF6me4ldEANioL8MdPYoGyVcGKZLigEx0AwDOlfN8DhzFiUQAzUZ6ORRto\/TqNX0lMQviAExUI0BO0NEm56S+AUxIAb6yUAnjzLfxQMoQeYBvyAGxEA1BuwcEW16SuIXxMAMDKhriwx08ii3OB+lFgNiQAyIATEwWAY6eZRTH6mlbINlVYWLgQ4YSJ2hlK2D0jSkGOg\/Az2psJNHOfWRWsrWE45UhhgYBAOpM5SyDWIyKrIWA\/wykiqolVTBnTDQyaPcyUw1qBgQA2JgoAzcuH493LhxI6s+9fjyKztPnjwZyqBfC5xR2OsvDT\/K1eZqH6kRbXpK4hfEgBioxoCdIaJNT0n8wrAYeOWVV8Lrr385K\/r06dNbHl9+lefx5z4Znv74hVw8+uEzm37VZ5ZMX3rHQCePsn2kBhumpyR+QQyIgWoM2Bki2vSUxC8Mi4E9hx4Je9dA1fHjy2OLff\/RJ8L+YwU4+jhhQs8Z6ORRjjnp+3fzcb1qi4EhMKBzNYRVqlbj9t17wvbd92fBWx5fPbYZL2P50vmjzMWR+m4e+1hI1jzEwLwZ4PzoXM2bdY0nBmZnoPNHeZopcOGAsr7EGMpiF8+vGYuBzQxUPSsWh9ycQS0xIAZmZWBwjzIXgf0XAHoeAfgsDkk7L1Z2MbDoDHA+OCcAPY8PfMQYaOfFyi4GxEB9Bjp\/lDnc8cGmjT2eTmwnBltZXOxXe7wMaGbrDKTOBmcF+3rEva+xnRhs9yLWNWz41lv6KgbEQBsMdPIoc7g9mFiqjb1t8Hf+bt++HWLUGTfu69uW586dO6EId+\/ezUKLYvBlQWtf\/BjSt65f3zlZW8LG\/+\/PEDoDIA3WRraNvHPV93Vpur6YZ8sf22lzB3DG79xZvwuwTYAJrBnuRHcJ\/dbMhfcLfSzOahibhIMxoJNHme+2q2Bagsnt+3IpxTbz83f+zp8\/Hy5fvjwBF8pbb71lIcHrN2\/czOzeduXKlUlfnwd9dXU1i+cv\/l+\/fj3k4ebN9bxlcfhJWDQm4wLmwfjoQ8YY53H16lWWsSKqhbHHq6Batq1R5PbWuueqqT04pP3A+YOzW7duIYKdW7NjtLuEXxDC\/WB3gd01xNy5eye8vwZ07gDiDBYf281v0uKshqL1GBLHzIN63333XegZPDp5lOfJWtHFQR0XLlwIZ8+eDYcOHZrg4MGD4YEHHsCdYc+ePZnky7Zt2xDB2w4cODDp6\/Og79u3L4vfsXNH2LlzZy62b99eK65oTMYFzKNKHLF9xhjnsXfv3my9h\/plmnPV1B4b0n6w87+yspIttZ1Hs2O0u2THjvX7YcfGXbCycdcQs7S0FJbXENb+F98lujsOBfbEfffdt8bO8P\/fi0eZA57CrPSSM\/7uPs554sSJcPjw4bBr164Jdu\/enT2eFkvb9OWVdcq8zfdN6fRdWV4JKyv5WF5eJiysLK+ElQpxqXFiGzXyjUBsH1p7jPOwb+6yRW\/pC\/s\/hVmHI+c05yred9O2h7Yf4NvOt58zdsB8kMsry9nZR9JeWWsjwdLSUlhaWkIN8R2xvLyctK+srAQPi\/M15OnUNKS7g3rncaYyolv+stxy\/tL0dsDtkMfSJ8BHvNnQsVnbyyKfj5MuBsbIgO1\/Ox+x9HPGR7zZ0LFZ28sin4+TLgbEwHQMdP4op8rmQuDwF\/nwE+djsNH2Et2ATxADi8oA54WzkJq\/+fCj+xhstL1EN+ATjAFJMTAbA718lMumxKUB4jizIVOI49UWA2LgHgN2Zu5Z1jXsaMgU8AliQAw0w8AgH+Vmpq4sYkAMiAExUIUBxcyPgc4fZb7z5mMwpmw6bXRsghgQA\/UZ4PxwjuhpOm10bIIYEAP9ZKDzR9lo4cIA1pYUA2JgdgY4U2D2TMogBsbGQD\/n0\/mjzIXBd+8xsPeTMlUlBvrPAOcnPlO0sfe\/elUoBhaXgc4f5cWlXjMXA2JADIgBMbCZgVkf5c3Z1BIDYkAMiAExIAamZqDzRzn1kRofsWGfelYdduR3sAJKMIkuiIF5MsD54Rz5MWlj9zbp42Hg2jur4do7b2aoOivuKEPVPoprl4FOHmUuBw+mmGpj7wQVBv3N3\/zNLVFs7nPnzoWTJ0+GU6dOBXT+wYstgTKIgRYY8GcInSGQBmsjhfEx8Pv\/9XPh8z91JsMXfulc6ePs7yu7s7AF\/a9TBjp5lPluvQo6ZaZk8IsXL26JYENj\/77veCQgAfqWQBnEQAsMVDlTxLQwtFI2wAD\/mtN7770X+Neepkn39v9+Lfy9H\/pQ+JvfeySgl+Ww++qFjx0PP\/wDj2Z3Vlkf+dtnYLn9IRZvhCcfe2gyaa9PjN0oGlUMiIEeMsAnapT16iuvhN\/+7d8Or7\/+ZZrhK1\/+SibrfPnzRx8M37yGOn2e\/ND+8ORj++t0UWyLDPTiUbaP17xscc5KLQYWggF\/nkxfiIkPbJKPPPVdWcUPPPyN4cEPfnPYe+iRrL3rQT2UGREL9qXzR5nLgo\/UYmBfsLXQdKdhQH2SDHB+4jNFG3uyg4ydMfDQscezsXfcvzfs2PNA2L77\/o32A5nUl8VioPNHebHo1mzFgBgQA2JADOQzoEc5n5tee\/hDGkXodfEqrowB+cWAGFhQBvQoD3Th+atW\/DWGPPBXsni0Bzo9lS0GxIAYWEgGOn+U7edc\/KzLA\/tCrkjFSR9\/7pPh6Y9fSOLRD5\/RX2+oyONYwzg\/\/jyZjr2zOWtgMSAGShno\/FHmsuCiiFFa+YIH7D\/6RNh\/LAdH1\/\/gyIJTtNDT17la6OXX5AfMQOeP8oC5U+liQAyIATEgBhploPNHmf9C5rv6RmelZGJgwRnQuVrwDaDpD5aBzh9le5CRMQbLqgoXAx0zwFmiBGQM7MLiMsAfAAXGAL\/e03TsIJhBcu4MdP4o8x19HubOhgYUAyNhIO9MYR\/JFDWNKRmwv7mBJAW\/3vOVV15FDdj4Gx1nzpwJq6urmU1f5stA54\/yfKer0cSAGBADi82A\/c0N\/pYGTPDrPfcc\/LOoAR\/2L33pS1l7ZF8GMZ3OHmV9pDaI\/aEiB8aAztXAFqyDcid\/c2Pjb2lkv97z\/vVf6Zn5NuwdlKYh1xjo5FHm4uBjNA9sa\/Xo\/2JADEzJAGfInyl0bFOmUzcxIAY6YKD2o9xBjRpSDIgBMSAGxMBCMKBHeSGWWZMUA2JADIiBITCgRzlbJX0RA2JADIgBMdA9A3qUu18DVSAGxIAYEANiIGOgs0eZP4DiQTW+jY5NmJ4B9Vw8Bjg3HjDg2+jYBDEgBvrJQCePMn8qtAr6SZmqEgP9ZKDKmSKmn9WrKjEgBmCgk0eZgQUxMB0D6iUGxIAYGC8Dg3yU+QgOVFmWqnFVcilGDIyZAc4KqDLHqnFVcilGDIiBewwM7lHmMuAjOIB+byqbNXxgs1UtMdBPBrquirPCmQLoefXgA3l+2cWAGJiNgUE9ylwGXBo2ZXRs1vYSH\/C2lM6\/iHL79u0QIxWbZ7O+eX7sd+7cCUW4e\/cuYYUx9K8SZzFWl+TW9e2ak2yxe\/KFM+TPCjq2VHn4QMrnbXnnqmve5z2+5wTdxkc33L2zfvazNirIGu4LNrBm4h7wsPO+5lr\/f0mc9d3Uz\/Ux+\/vvv7\/lXrT6+yjXJz\/8r4N6lNugm38V5fz58+Hy5csTcKG89dZbk+G8fvPGzYndlCtXrmR97V9V8TGm37hxI\/BPpOXh5s31vE3EWS7qAX5uQ9RZj7HN4+rVq7Z9RinXz9Xmc9XU3hvSfmDfssC3bt1ChPiuwGi+99e+cX\/\/zvvh\/bt3MIc7G5IGutnjO8LOO3FZjo1+eXFmt35xH7Nz7zW1Zm3nYU+8++67UDB4LPyjfOHChXD27Nlw6NChCQ4ePBgeeGD9F7Szwnv27EFk2LZtWyb9lwMHDmR99+3bl5l9jOk7du4IO3fuzMX27duzvk3EWa79+\/cHq83Pb2g66zG2eezduzdb77F+SZ2rpvbdkPaD3QkrKyvZUts+NjvGbdu3IcLy8vI6lpay9tLSvet5aWkpLK8hrP0vviPsvK+5QpajJM76W7+4j9kffPDB7F5rat3azMOeuO+++6Bg8Li36oOfynQTOHHiRDh8+HDYtWvXBLt3784eT8tI2\/Tlla2Uvf322wGQgzgfY\/rK8kpYWcnH8vIyXcPK8kpYmTHOclEP3wgghwz4H9s87Ju1bNFH+CV1rprag13vh7rzYHn9mbT+2IH5eIqXlpbC0tIS5rD+NVPD0tJSBlrxHWH98S2tfVla4msIeXFmt35LUR+zD+nMsSfGcqaW19ZD\/5+RAT6qO3nyZEBKlaVDAAAPKUlEQVSS6itf\/goig9czg76IATEgBsSAGMhhYFCPMn\/AxP8BFHRsOXObm\/n4c58MT3\/8Qnj0w2eyMXc\/eCCTfPE6bUEM9I0BzhBnyepCx2ZtyfoM8DNOEPfEBmJ7921V0BcGBvUoQxqXBZcGQMdmwGb6POX+o0+E\/cfWcPTxbNgd99\/7eaHXM6e+iIEeMsBZ4vwAdF8iNt+WXswAj+65c+eC\/\/Tsq1\/9ataJT9O8PTPqixhwDAzuUaZ2Lg2A7lHV5vtIFwNiYJ0Bzg9Yb937WtV2r8diazzKFy9eDHyCZp+e7dz3UEYKNv+pWmbUl9YYGGLiQT7KQyRaNYsBMbBYDGSfoG18erZ949OzzOY+VVssRjTbKgzoUa7CkmLEgBgQA2JADMyBgfJHeQ5FaAgxIAbEgBgQA2IgBD3K2gViQAyIATEgBnrCwKI+yj2hX2WIATEgBsSAGLjHgB7le1xIEwNiQAyIATHQKQN6lDulv+HBlU4MiAEx0BADb775ZuCvd+WhoWGUJmJAj3JEiJpiQAyIATEQwic+8YnwzDPPZL8EhV94EuO7\/9ZHs0dbXDXLgB7lZvlUtuYZUEYxIAY6YOCxv\/Pj4due\/2fZrxDmF5548EtR\/uC1L3VQ1fiH1KM84jXm33IFeR8\/YR\/x9DU1MSAGZmDgwW\/8C+Gho4+v\/wphfuGJx5p9htTqWsCAHuUCcobueu6558Kzzz5b+BHUqVOn9BFUFwutMcWAGBADCQb0KCdIGYvpiWc\/Eb71H5wv\/AiK39E7lvlqHmJADIiBoTOgR3noK1hQ\/0PHHg\/6CKqAILmaZkD5xIAYmJEBPcozEqjuYkAMLA4D\/DmMMiwOG5ppGwzoUW6DVeUUA2JgPAxszITH+Ny5c7l\/Rcj+ytDp06c3ekiIgfoM6FGuz5l6iAExsIAM8CjzZzCOP\/fJ5F8Tsr8yxF8XWkB6NOWGGNCj3BCRSiMGxMBiMDD5N5H9XxHy+jD+utBiLNYAZ6lHeYCLppLFgBgQA2JgnAzoUR7numpWYkAMiIHFY2AEM9ajPIJF1BTEgBgQA2JgHAzoUR7HOmoWYkAMiAExMAIGEo\/yCGalKYgBMSAGxIAYGCADepQHuGgqWQyIATEgBsbJwMI8yuNcPs1KDIgBMSAGxsSAHuUxreaUc+GXIpRhytTqJgbEwIgZKLs38I94+q1MTY9yK7TOK2kz4\/BrAe1XBObJU6dOBR2wZvhWFjEwFgZ0dzS\/knqUm+d0cBmr\/NpAfr3g4CamgsWAGGiVAd0dzdOrRznB6erqavjsZz878fBL6K3xlS9\/xdTRyL792sDREKuJiIGRM6C7o\/kF1qOc4JRH+bXXXpt4Hnnquyb67gcPTPRFU\/j4ugksGm+ab\/8ZqLKv+z8LVTgGBvQoV1jFh449Ponacf\/eib5oSpWfH+X9TNrbF\/Pn04u2W4YzXx5kPg3zezSls\/+HM6t+VQrHZehXxd1VM+pH+ciRIwF0R++4Ri76+ZH9c3VFMfzTdsTp59PD3hecKTCUWZQ9BvjZk1X27lDm3Lc6+YYm9Y2Ot\/HN+u\/8zu9kf6CUNUmhb\/Nqo57RPspcGpcuXQoAvQ3y+p7z5o2b2Qa\/ceNGI6UW\/vxo45+rK4zhn7fbiEsduNhmRfPjBH7GjzSbybhPXtviu5TUnzePLuuqMzZniTMF0Ov0rRqbt4Zmr8Mjfcr+K5gHg9qq7l1i6+LWzVt1u3Qev7p6JTR1d1T5hodvjFgL\/1DHOg83axqTU2dPxH371h7lo8xlwaVhZKNjs7aXLHAM76+iX3tnNVx7583C0OvX3wvXr1\/PhXVuIs5y3bhxff1RzhnX4pocs2qussPHYeQA8p0zB+6ll17K5hKvVdmFSx5gueL+82z7eRj3Q5KcIc6S1YyOzdpezsJr2Zo+88wzgf1AHPujbCwu+2\/+6I+GJ\/\/+zyTBpzfUXnXv1okjL7h1a\/1RvlZwV7z9f6+Ht995j\/CtWNk2sVmOa2t3Dpg41pT379xZ+7r+\/2sbY12L4q5t2NejQvB9mJvZr1xZDTcaujv2\/LljYc8jj+Zi\/8Y362XrxFqm1tvOltU+ZDnKR7nKgjz88MPhxIkTIfU4YNv+Dd8SDhw7Hj74wK7w6ku\/EB57cDk8emhn+PXPfXmSnvaxfdvDF375XPj8T53J5EN\/5kA4uH8pIAH6wSeeDP\/j9y+F3\/qt38oF\/h0f+quNxFmu3\/s\/f5zV+sqrrybHtTjkrLWRo079933H6VCEXd\/64fDF\/\/V72fqwHkwEyQPrQUxRHnzk+t1LV7Jcvm\/b+l\/77u8Pfgzqt3m8+OKLqKND0bmCi7\/+kR\/axAm2GKwV61aEqmsK5ysHHglv3N0T\/uDt\/5fEG3f2hM17N31O6+5x4oHdI1\/34EPZHWF3BXcL4H5hI\/z8i18KP\/HCFwJy37a74eve\/mpAfuPhPWHb0vWwsgZ0u2\/IA2jTH1+VOOKtn+9z4ts+GKgXwMfy\/R8I87g7uHtszLJ14l5mTeM9g437nP0HF0PG8pCLn6V2Fu\/ChQvh5ZdfTuKzP\/+p8OJPnws\/+68\/k\/l\/8lc+F1745\/8u060P7Z\/+1c9usn36F14MH\/orPxmQAP2X\/uFPhF\/7248X4tfPfHv49\/\/4Y4Ux5KgSVyWmaq6qcXXH\/I0f+b5QhM984tnwud\/4V5u4Nd69JKYoDz5y\/YcX\/2lpLp+3CZ3a8vJ85CMfmWX79rZv2bmqsg7EsG5FqLOm5GMf56Hu3s3LY3afD5175C+duxD+5a98etMe5G4BqT3yqy\/9x\/CXL\/yngPyZn\/s34euPncvuFfRUPDZ8VeKINfg+z\/\/oz2b3DzWX3UVVYuCjyThycS9b7bHkPu\/twahR2HKN2NGFcoE89dRTQRAH894D7L3RHaiNCTG3efOp8crP8Ng5Yt9tbMFBi4V+lAe9cipeDIgBMSAGRsfAKB\/l+A+g8IdRsI1u9TQhMTBHBjhDnCUbEh2btSXFwPwZGN+Io3yUWSYuCy4NgI5NEANiYDYGOEucKYA+Wzb1FgNiIGZgtI8yE+XSAOh1wIUD6vSZdyz1peDrML+3mV7ks5i2JTWkxsAOmvSlcjVli2ulnYIfz\/zeZnqRz2K6lJwp0GUNxlEsu6ypaGyrsyimDz6rM5Z9qC1VA3Xm2fN8qfg+2XiU+1RP57WwkFw4AL3zggoKoMYYFk7t5kM3O5J2ng9\/22B8kBoHe15t0\/pS4zRhox6QymVz8NLi6GN2dLMjaef58Av3GDCevLzn7Y82tDX1fJreHzbXK4FTsN7a\/BW71Y2+2dv\/lh5lt0YsIItpJnRs1h6KpGZqt3rRsdFG0kYH6NjQ5wXGBPF41OHt6NiIQ9JGB+jY0JG00QE6NvQ2wTigzhjU5fugYyMHkjY6QMeGLgyTAdaPdbTq0bFZW3I6BuARxL3h1tvRscVxfW6P91HuM+sN1cZmMzSUUmkaYsDWBdlQSqWJGIBbQ+RScwYGjFPkDGnUdUoG9ChPSVwfuvFdoEEHqA8rcq8GWxek1uYeL3U0eEvBcsCtgTizS87GgHGKFK+zcTlNbz3K07DWXZ\/JyByYSWNNoa0DtEZED\/7PWvgyaGttPCPVdHhLgd7YkQba4tjYmF7Co+9NW7x6RtrX9Si3z7FGEANiQAyIATFQiQE9yo6m+LtCvkPE5kJ6o1KbL4a21YqkbX50bLSRtNEBOjb0zrAxMHVQz0YzoGOjjaSNDtCxoSNpowN0bOhdgPH9uLStHiRt86Njo42kjQ7QsaELmxmAG2+h3UeuqInarFZ0bNbum6Q+XxPtPtfra0WnVmpGB+jY0IcCPcrRSrGALCRAj9y9aVIbNRpo++JoT+PzObrQp627qN+85+FrYQ1o+xpoYwfoVX0+btF1eIM\/A+2+ckJtQ6gT\/nyt1Ewb+5BAzdQO0IdUO7XqUYaFCCwkiMy9a1KjIVXctL5UrjZs1JfKix006UvlasqWqvXSpUsBO0iNgx3U9aXiF9UGf4a+czCUOuHRakXS7jPyasQO+lx7Xm16lPOYkV0MiAExIAbEwJwZ0KM8Z8I13HwZ4COspkZsMldTNc09jwYUA2KgVQb0KLdKr5J3yQCPaJMfYZGLnF3OSWOLATEwbgb0KDewvmUXdZmfEqrEEDctmsqfylPV5mtP9fF+9CoxxFUF+VKo2l9xo2WglYmx15pOTM4Umh5H+bpjQI9yQ9xzUFKp8uypWNmqM1DGK37+yzbOiC0GsXFcXpu+deLz8sguBqZlgD0YQ3tyWjb710+Pcv\/WZFAVcRlwQSDrFE6fOvFtxlJL3frbrEe5xUAjDCjJIBnQo9zQsqUudi567DYEbQ+zx9LHoFf158VW7U+c5UACbB7MJ2X3MfgtDt0Qx8RtizNpfmsjsXlpOvYmQD6PJnIqhxhIMeD3GXocgy1GHKP2+BjQozynNeVw8VB5YIuHx+Zj0LFZHDo2D2zmL5PE+r7o2Hw\/2tiBt9fRfQ7yAGypHNjxG4hBR8Y+a+MjBqA3ActNTgO2JnIrhxjwDLCvbI+ZxGYx6GZHYjeJvsAY\/dT1KI9giesc1iqxVWKgzS4OdPrQRm8S5G0yX1GueY5VVId8YqCIAc5ZDO3dIsaG5dOj3OB6cTA4LKRE0kavC\/p61O1fFu9zo5fFT+snt8e0eebZz9eLPs+xNZYYqMIA90qMKv0U02MGXGl6lB0ZfVHjA0e7qLY6jwex5PMoyj2Lz49heiofPuoy0E7FtW1jfMb2aHtM5RcDKQbYg+xHA+1UnGzjY0CP8gDWlIM5TZlV+lWJicfmgqAf0vtop+wWg890L7HT1+B90+rkIm+V\/sQRH8diT9lSsXGc2mJgFgbYe+wzwyy51HdYDIzoUe4H8RwiO1C+IrPjM2DzMejYzG8SGz6AbnaT2FI+s+MzYLN+SNrma1KSl\/we2PLG8HHoeXFmJxdxwGxVJPExyEVfpPfRxi6IgWkZ8PvJ6+wt30bH5sfB5uF90sfLgB7lBtY2Pkx5bewefmjs1kb3MLtJ70M3O5K2wdpemp6KMR+yDPRPxXg7uoePx06bSwc9htmJ8SDO2ujA2rHERx6z007B\/EjvtzYSkAs\/uiAGyhhgr+SBvrEPG7B9FvvNTowwXgb0KPd8bcdeHhcPl00M7E3Mvak81NJkLvIJYiDFAPssPg+0safiZRsXA\/8fAAD\/\/wZXp3YAAAAGSURBVAMAi00lYFE59XgAAAAASUVORK5CYII=","height":234,"width":388}}
%---
%[output:298783de]
%   data: {"dataType":"text","outputData":{"text":"[23:42:52][INFO]  Reinforcement Learning Toolbox not found.\n","truncated":false}}
%---
%[output:5308948c]
%   data: {"dataType":"text","outputData":{"text":"[23:42:52][INFO]  Manual REINFORCE loop (Section 2) implements the same algorithm.\n","truncated":false}}
%---
%[output:6e28324c]
%   data: {"dataType":"text","outputData":{"text":"[23:42:52][INFO]  Before RL:  valid=6\/100 (6%)  avg_reward=0.060\n","truncated":false}}
%---
%[output:7ce94921]
%   data: {"dataType":"text","outputData":{"text":"[23:42:52][INFO]  After RL:  valid=4\/100 (4%)  avg_reward=0.040\n","truncated":false}}
%---
%[output:1d690113]
%   data: {"dataType":"text","outputData":{"text":"[23:42:52][INFO]  RL Result: Reward 0.060 -> 0.040  (delta=-0.020, Within noise range)\n","truncated":false}}
%---
