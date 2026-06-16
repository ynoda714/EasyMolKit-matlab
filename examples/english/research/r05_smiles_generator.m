%[text] # R05: Molecular Language Model -- SMILES Generation
%[text] EasyMolKit Research -- Layer 4
%[text] 
%[text] ## Story
%[text] Molecules have "grammar."
%[text] Atomic symbols like C, N, O form branches with parentheses and close rings with numbers—strings following this rule are SMILES. Can a Markov chain, which has learned this "grammar" from training data, generate entirely new SMILES? In this script, we build a character vocabulary from a corpus of FDA-approved drug SMILES and generate new molecules using two models: a bigram Markov chain and a character-level LSTM.
%[text] We compare both models using evaluation metrics such as validity, diversity, and novelty, and explore the challenges and potential of generative models with small-scale data.
%[text] ## Learning Objectives
%[text] - Understand SMILES as a formal language and perceive language model training as probabilistic model learning of chemical structures.
%[text] - Construct a bigram Markov chain from training data (no gradient descent needed) and understand why inductive bias is crucial with small datasets.
%[text] - Diagnose SMILES generation failures by structural error types (parentheses, ring closure, valency) and relate them to memory requirements.
%[text] - Build and train a character-level LSTM with a custom training loop (backpropagation through time BPTT).
%[text] - Implement autoregressive decoding with temperature sampling.
%[text] - Master four standard metrics for generative models (validity, diversity, novelty, property distribution overlap).
%[text] - Understand REINFORCE (Williams 1992) as a policy gradient algorithm and grasp the connection between reinforcement learning and molecular generation (REINVENT) (details in r06\_reinforce.m).
%[text] - Recognize the scaling challenge where memory, not generalization, occurs with small corpora (production models use $10^6 \\sim 10^7$ molecules). \
%[text] ## Prerequisites
%[text] - Completion of A05 (Neural Network Property Prediction) -- Custom DL Loop
%[text] - Deep Learning Toolbox (dlnetwork, lstmLayer, adamupdate, dlfeval) \
%[text] ## Environment
%[text] - Deep Learning Toolbox is required for Sections 7-9.
%[text] - Compatible with both MATLAB Online and Desktop.
%[text] - No GPU needed; all training completes on a modern CPU. \
%[text] Estimated time: 45-90 minutes (REINFORCE RL section in r06\_reinforce.m adds 30-60 minutes)
%[text] ## Data
%[text] - data/list/fda\_drugs.csv -- 200 FDA-approved drugs (ChEMBL, CC-BY-SA 3.0) \
%[text] ## References
%[text] - Gomez-Bombarelli R et al. (2018) Automatic chemical design using a data-driven continuous representation of molecules. ACS Cent Sci 4:268-276. doi:10.1021/acscentsci.7b00572 \[Pioneering paper on VAE-based SMILES generation model\]
%[text] - Segler MHS et al. (2018) Generating focused molecule libraries for drug discovery with recurrent neural networks. ACS Cent Sci 4:120-131. doi:10.1021/acscentsci.7b00512 \[LSTM-based SMILES generation and transfer learning\]
%[text] - Olivecrona M et al. (2017) Molecular de novo design through deep reinforcement learning. J Cheminform 9:48. doi:10.1186/s13321-017-0235-x \[REINVENT: REINFORCE for molecular generation -- Foundation for Section 7\]
%[text] - Williams RJ (1992) Simple statistical gradient-following algorithms for connectionist reinforcement learning. Mach Learn 8:229-256. doi:10.1007/BF00992696 \[Original paper on REINFORCE algorithm\]
%[text] - Brown N et al. (2019) GuacaMol: benchmarking models for de novo molecular design. J Chem Inf Model 59:1096-1108. doi:10.1021/acs.jcim.8b00839 \[Metrics for generative model validity/diversity/novelty/KL divergence\]
%[text] - Polykovskiy D et al. (2020) Molecular sets (MOSES): a benchmarking platform for molecular generation models. Front Pharmacol 11:565644. doi:10.3389/fphar.2020.565644 \[MOSES benchmark: standard evaluation protocol for generative models\]
%[text] - Hochreiter S & Schmidhuber J (1997) Long short-term memory. Neural Comput 9:1735-1780. doi:10.1162/neco.1997.9.8.1735 \[Original paper on LSTM\] \
%%
%[text] ## Section 0: Setup and Configuration
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
rehash toolboxcache;   % flush package cache -- ensures newly added src/ files are visible
emk.setup.initPython(); %[output:1b68a8eb]
%[text] Check toolbox availability.
hasDL = license("test", "Neural_Network_Toolbox");
if ~hasDL
    error("emk:r05:missingToolbox", ...
        "Deep Learning Toolbox is required for Sections 7-9.");
end
logInfo("Toolbox: Deep Learning=%d", hasDL); %[output:0902b7b6]
%[text] ### Section 0a: Adjustable Parameters (Change Here Before Execution)
HIDDEN_SIZE   = 64;     % Number of LSTM hidden units -- auto-scaled based on CORPUS_TARGET below
DROPOUT_RATE  = 0.30;   % Dropout probability
LEARN_RATE    = 8e-4;   % Initial learning rate for Adam -- auto-scaled based on CORPUS_TARGET below
N_EPOCHS      = 150;    % Number of SL training epochs
                        %   Plateau of val loss becomes clear at 150 (key teaching point);
                        %   Smoke test can be shortened to 20
BATCH_SIZE    = 32;     % Mini-batch size (must be less than n_train)
MAX_SEQ_LEN   = 100;    % Maximum sequence length after padding (including START + END tokens; 4.5% of corpus truncated at 90)
TEMPERATURE   = 0.85;   % Sampling temperature -- auto-scaled based on CORPUS_TARGET below
                        %   Lower to 0.6-0.75 only after sustained efficacy rate >20%
N_GENERATED   = 300;    % Number of generated SMILES for evaluation
                        %   300 results in binomial SE ~2.6% (compared to ~3.5% for 100),
                        %   allowing reliable comparison of Markov vs LSTM; generation takes ~4.5 min with max sequence length 100
CORPUS_TARGET = 500;    % Number of target molecules for training
                        %   200: Local CSV only (no network required; ~1 min; Markov ~13%, LSTM ~10%)
                        %   500: Local CSV + ChEMBL Phase-4 retrieval (initial ~50 sec; total ~5 min)
                        %        --> Markov ~18%, LSTM ~4%: contrast is a teaching point
                        %   Change to 200 if network is unstable
AUGMENT_FACTOR = 1;     % Set to 1: Augmentation adds only noise at this data scale

%[text] Automatically scale architecture based on corpus size.
%[text] For 200 molecules, HIDDEN=64 results in parameters/tokens $\\approx 4$ (memorization regime).
%[text] For 500+ molecules, it becomes $\\approx 1.7$ leading to underfitting, so expand HIDDEN\_SIZE to 128.
%[text] Rule of thumb: Keep parameters/tokens around 3-6 for LSTM character LM.
EMBED_DIM = 16;    % Character embedding dimension (VOCAB_SIZE → EMBED_DIM → LSTM)

if CORPUS_TARGET >= 500
    HIDDEN_SIZE = 128;     % 128 hidden units for 500 mol corpus
    N_EPOCHS    = 150;
    TEMPERATURE = 0.80;
end
if CORPUS_TARGET >= 1000
    HIDDEN_SIZE = 256;
    N_EPOCHS    = 100;
    TEMPERATURE = 0.75;
end
%[text] Gradient clipping threshold: Limit L2 norm per tensor to 1.0 to stabilize LSTM training.
GRAD_CLIP_NORM = 1.0;
%[text] Define special token characters (must never appear in SMILES strings).
START_TOKEN = "^";
END_TOKEN   = "$";
PAD_TOKEN   = "_";

logInfo("R05: Setup complete (HIDDEN=%d  EPOCHS=%d  TEMP=%.2f  AUG=%dx)", ... %[output:group:10372ce6] %[output:334b2162]
    HIDDEN_SIZE, N_EPOCHS, TEMPERATURE, AUGMENT_FACTOR); %[output:group:10372ce6] %[output:334b2162]
%[text] Warm up Python/RDKit (preempt initial delay before training).
mol_warmup = emk.mol.fromSmiles("C");
clear mol_warmup;
%%
%[text] ## Section 1: Loading SMILES Corpus
%[text] ### SMILES as a Language
%[text] SMILES (Simplified Molecular Input Line Entry System) encodes molecular structures as compact ASCII strings: atoms, bonds, rings, and stereocenters each map to one or two characters following a formal grammar.
%[text] Canonical SMILES give a unique string per molecule; non-canonical SMILES can describe the same structure in multiple ways. Because SMILES follow a formal grammar, a language model trained on them can:
%[text] - Generate syntactically valid strings.
%[text] - Implicitly encode chemical rules (ring closure balance, valency).
%[text] - Interpolate between structures in chemical space. \
%[text] The corpus here consists of 200--500 FDA-approved drugs from ChEMBL.
%[text] We deliberately use a small dataset so training completes in minutes on a CPU and overfitting can be observed directly.
%[text] Production models (REINVENT, MolGPT) are trained on $10^6$--$10^7$ molecules from ZINC/ChEMBL.
%[text] 
%[text] \---- Loading Multi-Source Corpus ----
%[text] Indicates source priority (trimmed to CORPUS\_TARGET after deduplication).
%[text] CORPUS\_TARGET ≤ 200: Use only fda\_drugs.csv (200 curated FDA drugs).
%[text] CORPUS\_TARGET 201-500: Use only fda\_drugs.csv + ChEMBL Phase-4 approved drugs (everyday\_chemicals / forensic\_challenge are intentionally excluded to prevent domain mixing from degrading model quality at small N).
%[text] CORPUS\_TARGET \> 500: Use all local CSV + ChEMBL.
smilesCombined = strings(0, 1);

t1 = readtable("data/list/fda_drugs.csv", TextType="string");
smilesCombined = [smilesCombined; t1.SMILES];
logInfo("Local fda_drugs.csv:          %d rows", height(t1)); %[output:3fd1b0d0]

if CORPUS_TARGET > 500 %[output:group:16007640]
    % Mix heterogeneous sources only for large N
    t2 = readtable("data/list/everyday_chemicals.csv", TextType="string");
    smilesCombined = [smilesCombined; t2.SMILES];
    logInfo("Local everyday_chemicals.csv: %d rows added", height(t2));

    t3 = readtable("data/list/forensic_challenge.csv", TextType="string");
    smilesCombined = [smilesCombined; t3.SMILES];
    logInfo("Local forensic_challenge.csv: %d rows added", height(t3));
elseif CORPUS_TARGET > height(t1)
    logInfo("Corpus 201-500: Using only FDA drugs + ChEMBL approved drugs"); %[output:525fc4be]
    logInfo("  (Excluding everyday/forensic -- domain mixing degrades small N models)"); %[output:751f8380]
end %[output:group:16007640]
logInfo("Local corpus: %d SMILES (before ChEMBL fetch)", numel(smilesCombined)); %[output:422f2707]

CHEMBL_CACHE = "data/list/chembl_extended.csv";
if CORPUS_TARGET > numel(smilesCombined) %[output:group:1a5839c5]
    nFetch = CORPUS_TARGET - numel(smilesCombined);
    cacheInsufficient = false;
    if isfile(CHEMBL_CACHE)
        cacheT = readtable(CHEMBL_CACHE, TextType="string");
        useN   = min(height(cacheT), nFetch + 150);  % Buffer window for deduplication loss
        if height(cacheT) < nFetch + 100
            % Cache is too small (e.g., fewer than needed were fetched in the last run)
            logWarn("ChEMBL cache has only %d rows. ~%d rows needed. Refetching ...", ...
                height(cacheT), nFetch + 150);
            cacheInsufficient = true;
        else
            logInfo("ChEMBL cache hit: %d / %d rows loaded", useN, height(cacheT)); %[output:956d762e]
            smilesCombined = [smilesCombined; cacheT.SMILES(1:useN)];
        end
    end
    if ~isfile(CHEMBL_CACHE) || cacheInsufficient
        logInfo("Fetching ~%d SMILES from ChEMBL (first time ~50 seconds) ...", nFetch);
        BASE_URL_C  = "https://www.ebi.ac.uk/chembl/api/data/molecule.json";
        FILTER_C    = "max_phase=4&molecule_type=Small+molecule&withdrawn_flag=0";
        PAGE_SIZE_C = 100;
        opts_c      = weboptions("Timeout", 30, "ContentType", "json");
        fetchedSmiles = strings(0, 1);
        offset_c = 200;   % Skip 200 overlaps with fda_drugs.csv
        nWant    = nFetch + 150;  % Large buffer to ensure enough after deduplication
        while numel(fetchedSmiles) < nWant
            lim_c = min(PAGE_SIZE_C, nWant - numel(fetchedSmiles));
            url_c = sprintf("%s?%s&limit=%d&offset=%d", BASE_URL_C, FILTER_C, lim_c, offset_c);
            try
                data_c = webread(url_c, opts_c);
            catch ME
                logWarn("Stopped ChEMBL fetch at offset=%d: %s", offset_c, ME.message);
                break;
            end
            mols_c = data_c.molecules;
            if isempty(mols_c); break; end
            for ki = 1:numel(mols_c)
                m = mols_c(ki);
                if isstruct(m.molecule_structures) && ...
                        isfield(m.molecule_structures, "canonical_smiles") && ...
                        strlength(string(m.molecule_structures.canonical_smiles)) > 0
                    fetchedSmiles(end+1) = string(m.molecule_structures.canonical_smiles); %#ok<AGROW>
                end
            end
            logInfo("  ChEMBL offset=%d -> %d fetched", offset_c, numel(fetchedSmiles));
            offset_c = offset_c + numel(mols_c);
            pause(0.3);
            if numel(mols_c) < lim_c; break; end
        end
        if ~isempty(fetchedSmiles)
            [~, ~] = mkdir(fileparts(CHEMBL_CACHE));
            writetable(array2table(fetchedSmiles(:), 'VariableNames', {'SMILES'}), CHEMBL_CACHE);
            logInfo("%d ChEMBL SMILES cached -> %s", numel(fetchedSmiles), CHEMBL_CACHE);
            smilesCombined = [smilesCombined; fetchedSmiles(:)];
        else
            logWarn("ChEMBL fetch resulted in 0 SMILES -- check network. Using local data only.");
        end
    end
end %[output:group:1a5839c5]
%[text] After deduplication, trim to CORPUS\_TARGET.
smilesCombined = unique(smilesCombined);
logInfo("After deduplication: %d unique SMILES", numel(smilesCombined)); %[output:55200a7a]
if numel(smilesCombined) > CORPUS_TARGET
    smilesCombined = smilesCombined(1:CORPUS_TARGET);
end
logInfo("Training corpus: %d (target=%d)", numel(smilesCombined), CORPUS_TARGET); %[output:0617c68c]
%[text] Validate SMILES with RDKit and exclude unprocessable rows.
nSmilesCombined = numel(smilesCombined);
validMask = false(nSmilesCombined, 1);
for k = 1:nSmilesCombined %[output:group:8a4fbe84]
    try
        validMask(k) = emk.mol.isValid(smilesCombined(k));
    catch
        validMask(k) = false;
    end
    if mod(k, max(1, round(nSmilesCombined/10))) == 0 || k == nSmilesCombined
        logProgress(k, nSmilesCombined, "Validating SMILES"); %[output:6cc88e94]
    end
end %[output:group:8a4fbe84]
smilesAll = smilesCombined(validMask);
N_MOLS    = numel(smilesAll);
logInfo("Valid SMILES: %d / %d", N_MOLS, numel(smilesCombined)); %[output:1340b2ac]
%[text] Check the distribution of string lengths. Sequences exceeding MAX\_SEQ\_LEN are truncated.
lengths = cellfun(@strlength, cellstr(smilesAll));
logInfo("SMILES length: min=%d  median=%.0f  max=%d  mean=%.1f", ... %[output:group:6e60e7e4] %[output:73633f29]
    min(lengths), median(lengths), max(lengths), mean(lengths)); %[output:group:6e60e7e4] %[output:73633f29]

pctTruncated = 100 * mean(lengths > MAX_SEQ_LEN - 2);
logInfo("Truncated sequences (> %d chars): %.1f%%", MAX_SEQ_LEN-2, pctTruncated); %[output:9a34ceb7]

figure("Name", "R05 SMILES Length Distribution"); %[output:45463a39]
histogram(lengths, 20, FaceColor=[0.2 0.5 0.8]); %[output:45463a39]
xlabel("SMILES Length (chars)"); ylabel("Count"); %[output:45463a39]
title(sprintf("Distribution of SMILES Lengths (N=%d)", N_MOLS)); %[output:45463a39]
xline(MAX_SEQ_LEN - 2, "r--", "Max Sequence Length", LabelHorizontalAlignment="left"); %[output:45463a39]
grid on; %[output:45463a39]
%[text] Save the set of canonical SMILES before performing SMILES enumeration data augmentation (Bjerrum 2017: arXiv:1703.07076).
%[text] This is used for novelty checks (Section 4, 9) to avoid mistakenly counting non-canonical duplicates as "new molecules". It is also used in property distribution plots (Section 9) to reflect the original data distribution.
smilesOriginal = smilesAll;
N_MOLS_ORIG    = N_MOLS;
if AUGMENT_FACTOR > 1 %[output:group:76149936]
    logInfo("Augmenting corpus: %d mol x up to %d random SMILES ...", N_MOLS, AUGMENT_FACTOR);
    aug_mods = emk.util.rdkitModule();  % cached; no extra IPC cost
    augBuf = strings(0, 1);
    for k = 1:N_MOLS_ORIG
        mol_k = emk.mol.fromSmiles(smilesAll(k));
        rands = strings(AUGMENT_FACTOR, 1);
        for ri = 1:AUGMENT_FACTOR
            try
                rands(ri) = string(aug_mods.Chem.MolToSmiles(mol_k, ...
                    pyargs("canonical", false, "doRandom", true)));
            catch
                rands(ri) = smilesAll(k);  % Fallback to canonical SMILES on error
            end
        end
        augBuf = [augBuf; unique(rands)]; %#ok<AGROW>
        if mod(k, max(1, round(N_MOLS_ORIG/10))) == 0 || k == N_MOLS_ORIG
            logProgress(k, N_MOLS_ORIG, "Augmenting SMILES");
        end
    end
    augBuf  = unique(augBuf);
    augBuf  = augBuf(~ismember(augBuf, smilesAll));  % Remove canonical duplicates
    smilesAll = [smilesAll; augBuf];
    N_MOLS    = numel(smilesAll);
    logInfo("Post-augmentation corpus: %d sequences (%.1fx)", N_MOLS, N_MOLS / N_MOLS_ORIG);
else
    logInfo("AUGMENT_FACTOR=1: Augmentation disabled (corpus = %d unique SMILES)", N_MOLS); %[output:233a680b]
end %[output:group:76149936]
%%
%[text] ## Section 2: Building Character Vocabulary
%[text] ### Tokenization and Special Tokens
%[text] In a character-level language model, each unique character in the corpus becomes a vocabulary token. The vocabulary also includes three special tokens that never appear in SMILES strings.
%[text] - START ("^"): Added at the beginning of the entire sequence. It serves as the starting point for generation.
%[text] - END ("\\$"): Added at the end to indicate molecule completion. Training stops loss calculation at this position.
%[text] - PAD ("\_"): Fills beyond the actual sequence end to create fixed-length matrices within a mini-batch. PAD positions are masked from the loss function. \
%[text] Additionally, Cl and Br are replaced with single-character placeholders (L and R) before building the vocabulary. Without replacement, 'l' and 'r' could be sampled after any token, generating meaningless sequences.
%[text] In the output of generateSmiles, L→Cl and R→Br are restored.
%[text] 
%[text] Preprocessing: Replace two-character atoms with single-character placeholders. L = Cl (Chlorine), R = Br (Bromine).
smilesProc = strrep(strrep(smilesAll, "Cl", "L"), "Br", "R");

%[text] num2cell converts a 1×M char row into a cell array of individual characters.
%[text] The three special token cells can be directly concatenated with horzcat.
allCharsCell = num2cell(unique(char(strjoin(smilesProc, ""))));
specialChars  = {char(PAD_TOKEN), char(START_TOKEN), char(END_TOKEN)};
inSpecial     = cellfun(@(c) any(strcmp(c, specialChars)), allCharsCell);
vocabChars    = [specialChars, sort(allCharsCell(~inSpecial))];
VOCAB_SIZE    = numel(vocabChars);

%[text] Warn if characters outside the SMILES token set are included in the vocabulary.
%[text] char(92) is a backslash ('' in E/Z stereochemistry notation). Referenced numerically to avoid JSON escape issues.
%[text] 'a' is an aromatic wildcard atom used in some ChEMBL canonical SMILES, valid in RDKit.
SMILES_ALLOWED = {'B','C','N','O','P','S','F','I','L','R', ...
    'a','b','c','n','o','p','s', ...
    '#','=','-','+','.', ...
    '(',')', '[',']', ...
    '0','1','2','3','4','5','6','7','8','9', ...
    '@','/', char(92), '%','H'};
unexpectedChars = vocabChars(~ismember(vocabChars, [specialChars, SMILES_ALLOWED]));
if ~isempty(unexpectedChars) %[output:group:415def2d]
    logWarn("Unexpected characters included in vocabulary: [%s] -- Generation quality may degrade.", ... %[output:6ccf3410]
        strjoin(string(unexpectedChars), " ")); %[output:6ccf3410]
end %[output:group:415def2d]

char2idx = containers.Map(vocabChars, num2cell(1:VOCAB_SIZE));
idx2char  = containers.Map(num2cell(1:VOCAB_SIZE), vocabChars);

PAD_IDX   = char2idx(char(PAD_TOKEN));
START_IDX = char2idx(char(START_TOKEN));
END_IDX   = char2idx(char(END_TOKEN));

logInfo("Vocabulary size: %d characters", VOCAB_SIZE); %[output:7c8775f2]
logInfo("Special tokens: PAD=%d  START=%d  END=%d", PAD_IDX, START_IDX, END_IDX); %[output:329ce97e]
%%
%[text] ## Section 3: Bigram Markov Chain -- Simple Generation Baseline
%[text] ### n-gram Language Model
%[text] Before deep learning, sequence generation relied on n-gram models.
%[text] The bigram ($n=2$) model estimates the conditional probability:
%[text]{"align":"center"} $P(c(t) \\mid c(t-1)) = \\frac{\\mathrm{count}(c(t-1),\\, c(t))}{\\mathrm{count}(c(t-1))}$
%[text] The probability of each character depends only on the immediately preceding character.
%[text] In code, this is the `transCount` matrix: each observed bigram $(c(t-1), c(t))$ increments `transCount(c(t-1), c(t))` by 1.
%[text] Despite its simplicity, the bigram model works surprisingly well on small SMILES data because many local transitions are chemically constrained: after '(', an atom (C, N, O, S...) almost always follows; after '=' or '\#', an atom always follows (bond order); after '\[', an atomic symbol always follows.
%[text] About 200 molecules is enough to learn these local patterns reliably.
%[text] Laplace smoothing ($\\alpha = 0.1$) prevents zero probabilities; each row is then normalized to a proper distribution. For comparison: the Markov model needs no training and reaches 20--35% validity from 200 molecules, while LSTM on the same data achieves only 5--15% (underfitting). With the default 500-molecule corpus, increased vocabulary diversity reduces the Markov validity rate to 5--15%.
%[text] 
%[text] Construct a bigram transition count matrix from training SMILES (using smilesProc, with Cl/Br replaced in Section 2).
transCount = zeros(VOCAB_SIZE, VOCAB_SIZE, "double");
for mk = 1:N_MOLS
    s = char(smilesProc(mk));
    % Encode sequence: START + characters + END
    nChar  = numel(s);
    idxSeq = zeros(1, nChar + 2);
    idxSeq(1) = START_IDX;
    for ci = 1:nChar
        if isKey(char2idx, s(ci))
            idxSeq(ci + 1) = char2idx(s(ci));
        else
            idxSeq(ci + 1) = PAD_IDX;   % Unknown character -> PAD
        end
    end
    idxSeq(end) = END_IDX;
    % Count all consecutive pairs
    for ti = 1:numel(idxSeq) - 1
        c1 = idxSeq(ti);
        c2 = idxSeq(ti + 1);
        if c1 >= 1 && c1 <= VOCAB_SIZE && c2 >= 1 && c2 <= VOCAB_SIZE
            transCount(c1, c2) = transCount(c1, c2) + 1;
        end
    end
end
%[text] Laplace smoothing ($\\alpha = 0.1$): Add pseudo-counts to each cell to avoid zero probabilities.
%[text] PAD is an absorbing state (self-loop) and does not transition to other tokens.
LAPLACE_ALPHA    = 0.1;
transCountSmooth = transCount + LAPLACE_ALPHA;
transCountSmooth(PAD_IDX, :)          = 0;   % PAD does not transition
transCountSmooth(PAD_IDX, PAD_IDX)    = 1;   % Absorbing state for normalization
%[text] Normalize each row to obtain a probability distribution.
rowSums   = sum(transCountSmooth, 2);
rowSums(rowSums == 0) = 1;          % Guard against all-zero rows
transProb = transCountSmooth ./ rowSums;

logInfo("Markov model construction complete: %dx%d transition table (from %d SMILES, %d bigrams)", ... %[output:group:4b3ee340] %[output:18cc59a8]
    VOCAB_SIZE, VOCAB_SIZE, N_MOLS, round(sum(transCount(:)))); %[output:group:4b3ee340] %[output:18cc59a8]
%[text] Visualize the transition heatmap of major SMILES characters. Rows correspond to "current character" and columns to "next character".
HEAT_CHARS = {'C','N','O','S','F','(',')','+','-','1','2','3','=','#','[',']'};
heatIdx = zeros(1, numel(HEAT_CHARS));  heatN = 0;
for hi = 1:numel(HEAT_CHARS)
    if isKey(char2idx, HEAT_CHARS{hi})
        heatN = heatN + 1;
        heatIdx(heatN) = char2idx(HEAT_CHARS{hi});
    end
end
heatIdx    = heatIdx(1:heatN);
heatLabels = cellfun(@(i) idx2char(i), num2cell(heatIdx), UniformOutput=false);
subMatrix  = transProb(heatIdx, heatIdx);

figure("Name", "R05 Bigram Transition Probabilities"); %[output:2f235c0c]
imagesc(subMatrix); %[output:2f235c0c]
colormap("hot");  colorbar; %[output:2f235c0c]
xticks(1:heatN);  xticklabels(heatLabels);  xtickangle(45); %[output:2f235c0c]
yticks(1:heatN);  yticklabels(heatLabels); %[output:2f235c0c]
xlabel("Next Character");  ylabel("Current Character"); %[output:2f235c0c]
title(sprintf("Bigram Transition Probabilities (Top %d Characters)", heatN)); %[output:2f235c0c]
%%
%[text] ## Section 4: Generation and Evaluation by Markov Model
%[text] ### GuacaMol / MOSES Evaluation Protocol
%[text] We use three standard metrics for molecular generative models (Brown et al. 2019; Polykovskiy et al. 2020):
%[text] - **Validity**: Fraction of generated SMILES parseable by RDKit.
%[text] - **Uniqueness**: Fraction of valid SMILES that are structurally distinct.
%[text] - **Novelty**: Fraction of unique valid SMILES absent from the training set. \
%[text] Production targets: Validity \> 90%, Uniqueness \> 85%, Novelty \> 60%.
%[text] The Markov baseline achieves 20--35% validity with 200 molecules, dropping to 5--15% at 500 (increased vocabulary diversity).
%[text] It is a strong zero-training-time baseline but falls short of production targets -- Sections 7--8 show that LSTM is also insufficient at small scale.
N_MARKOV = N_GENERATED;   % Same number for direct comparison with LSTM
logInfo("Generating %d SMILES with bigram Markov model ...", N_MARKOV); %[output:235d5d03]
markovSmiles = strings(N_MARKOV, 1);
for mg = 1:N_MARKOV %[output:group:8d9a6cf7]
    markovSmiles(mg) = generateSmilesMarkov(transProb, START_IDX, END_IDX, ...
        PAD_IDX, VOCAB_SIZE, MAX_SEQ_LEN, TEMPERATURE, idx2char, char2idx);
    if mod(mg, max(1, round(N_MARKOV/10))) == 0 || mg == N_MARKOV
        logProgress(mg, N_MARKOV, "Generating (Markov)"); %[output:4a72d9ff]
    end
end %[output:group:8d9a6cf7]

%[text] Validate the generated SMILES with RDKit.
logInfo("Validating Markov SMILES with RDKit ..."); %[output:353a9a4a]
markovValid = false(N_MARKOV, 1);
for mg = 1:N_MARKOV %[output:group:345dfb99]
    try
        if strlength(markovSmiles(mg)) > 0
            markovValid(mg) = emk.mol.isValid(markovSmiles(mg));
        end
    catch
        markovValid(mg) = false;
    end
    if mod(mg, max(1, round(N_MARKOV/10))) == 0 || mg == N_MARKOV
        logProgress(mg, N_MARKOV, "Validating (Markov)"); %[output:7b36f722]
    end
end %[output:group:345dfb99]

nMValid  = sum(markovValid);
mValidSm = markovSmiles(markovValid);
nMUnique = numel(unique(mValidSm));
nMNovel  = sum(~ismember(mValidSm, smilesOriginal));  % Novelty against the regular training set

logInfo("=== Markov Baseline Results (%d generated) ===", N_MARKOV); %[output:8098b4ea]
logInfo("  Validity:   %d/%d  (%.1f%%)", nMValid, N_MARKOV, 100*nMValid/N_MARKOV); %[output:0fd21ee6]
logInfo("  Uniqueness: %d/%d  (%.1f%%)", nMUnique, max(nMValid,1), ... %[output:group:3ed5b7d2] %[output:4f3967e5]
    100*nMUnique/max(nMValid,1)); %[output:group:3ed5b7d2] %[output:4f3967e5]
logInfo("  Novelty:    %d/%d  (%.1f%%)", nMNovel, max(nMValid,1), ... %[output:group:104abef8] %[output:9d3c70f5]
    100*nMNovel/max(nMValid,1)); %[output:group:104abef8] %[output:9d3c70f5]
logInfo("  (Production targets: Validity>90%%  Uniqueness>85%%  Novelty>60%%)"); %[output:06b087e6]

%[text] Display the first 5 sample SMILES for diagnostics, regardless of validity.
logInfo("Sample Markov SMILES (first 5):"); %[output:8f440e59]
for mg_s = 1:min(5, N_MARKOV) %[output:group:2366cbc7]
    logInfo("  [%2d] %s", mg_s, markovSmiles(mg_s)); %[output:7147febe]
end %[output:group:2366cbc7]

if nMValid > 0 %[output:group:022431b9]
    logInfo("First %d valid Markov SMILES:", min(5, nMValid)); %[output:2a8966f4]
    for mg = 1:min(5, nMValid)
        logInfo("  [%2d] %s", mg, mValidSm(mg)); %[output:067fae5f]
    end
end %[output:group:022431b9]
%%
%[text] ## Section 5: Limitations of Markov Chains -- Motivation for LSTM
%[text] ### Long-Distance Dependency Problem
%[text] Markov models reproduce local transitions accurately but fail systematically for structural constraints that span more than one character.
%[text] **Problem 1: Parenthesis Matching (Branch Notation)**
%[text] Drug SMILES contain nested branches: in CC(=O)O (acetic acid), each '(' must be closed by a matching ')'. A bigram model has no memory of open branches.
%[text] Grammar rules 2 and 5 in this script enforce parenthesis matching deterministically, substantially reducing mismatch errors at runtime.
%[text] **Problem 2: Ring Closure Matching**
%[text] In C1CCCCC1 (cyclohexane), the digit 1 after the first C opens a ring bond; it must appear again after five more atoms to close it. This dependency spans 7 characters -- far beyond a bigram's one-character lookahead.
%[text] Grammar rules 2c and 5 track open ring-closure numbers and force closure, suppressing ring-closure mismatches at runtime.
%[text] **Problem 3: Valency Consistency**
%[text] Carbon always forms exactly 4 bonds. Bigrams cannot count how many bonds a carbon has already formed from previous tokens.
%[text] **Problem 4: Ring-Bond Accumulation**
%[text] Drug SMILES routinely open ring bonds with digits 1--6. Without a depth limit, unclosed bonds accumulate and cannot all be resolved within the forced-closure window.
%[text] Grammar rule 2c caps simultaneously open ring bonds at three to mitigate this.
%[text] **Solution: LSTM**
%[text] LSTM (Hochreiter & Schmidhuber 1997) maintains a hidden state $h(t)$ and cell state $c(t)$ that carry context across arbitrary distances.
%[text] At 10,000+ molecules, LSTM can implicitly track open parentheses and ring bonds in its hidden state, achieving validity rates above 90%.
%[text] In summary: Markov is data-efficient (bigram counting works well at 200 molecules), but LSTM underfits at 500 molecules and needs far more data to learn long-distance patterns. Production language models (MolGPT, ChemGPT) use $10^6$--$10^7$ molecules.
%[text] 
%[text] We categorize and count four failure modes in Markov-generated SMILES: parenthesis mismatch, bracket mismatch, ring-closure mismatch, and other.
%[text] Grammar constraints (rules 2, 2c, 5) suppress the first two categories, keeping total mismatches below 10%; "other" (valency errors, aromaticity errors, invalid atom symbols) accounts for over 90% of the remaining failures.
%[text] This is the heart of Problem 3. Parenthesis and ring-closure "memory" can be handled with a simple depth counter -- a deterministic rule over a single integer. Valency tracking, by contrast, requires recording how many bonds every atom in the growing sequence currently holds; the state space explodes with branches and rings. Valency correctness must therefore be learned as chemical knowledge from data, not patched with grammar rules.
markovInvalid = markovSmiles(~markovValid);
logInfo("Markov Failure Analysis (%d invalid / %d generated):", ... %[output:group:7ebbdda9] %[output:8b7375fb]
    numel(markovInvalid), N_MARKOV); %[output:group:7ebbdda9] %[output:8b7375fb]
if numel(markovInvalid) > 0 %[output:group:9b729d14]
    nParenErr   = sum(cellfun(@(s) count(char(s), '(') ~= count(char(s), ')'), ...
        cellstr(markovInvalid)));
    nBracketErr = sum(cellfun(@(s) count(char(s), '[') ~= count(char(s), ']'), ...
        cellstr(markovInvalid)));
    % Ring closure mismatch: Numbers appearing an odd number of times remain open (or vice versa).
    nRingErr = sum(cellfun(@(s) any(arrayfun(@(d) mod(count(char(s), char('0'+d)), 2) ~= 0, 1:9)), ...
        cellstr(markovInvalid)));
    nOtherErr = numel(markovInvalid) - nParenErr - nBracketErr - nRingErr;
    logInfo("  Parenthesis Mismatch '()':  %d / %d  (%.0f%%)", ... %[output:5276edff]
        nParenErr,   numel(markovInvalid), 100*nParenErr/numel(markovInvalid)); %[output:5276edff]
    logInfo("  Bracket Mismatch '[]':     %d / %d  (%.0f%%)", ... %[output:83699159]
        nBracketErr, numel(markovInvalid), 100*nBracketErr/numel(markovInvalid)); %[output:83699159]
    logInfo("  Ring Closure Number Mismatch:  %d / %d  (%.0f%%)", ... %[output:14d3854d]
        nRingErr,    numel(markovInvalid), 100*nRingErr/numel(markovInvalid)); %[output:14d3854d]
    logInfo("  Others (Valency/Atoms/Invalid): %d / %d  (%.0f%%)", ... %[output:9f6c2440]
        max(nOtherErr,0), numel(markovInvalid), 100*max(nOtherErr,0)/numel(markovInvalid)); %[output:9f6c2440]
end %[output:group:9b729d14]
logInfo("Next Steps: Build LSTM character language model in Sections 6-8."); %[output:8ce954da]
%%
%[text] ## Section 6: Encoding, Padding, and Data Splitting
%[text] ### Training Objective: Next-Token Prediction
%[text] Training a character-level language model reduces to next-token prediction.
%[text] For each sequence $c(1), c(2), \\ldots, c(T)$, the input is $\[\\mathrm{START},\\, c(1), \\ldots, c(T-1)\]$ and the target is $\[c(1), c(2), \\ldots, c(T),\\, \\mathrm{END}\]$, both of length $T$.
%[text] At each position $t$, the model sees $c(1) \\ldots c(t)$ and predicts the next character $c(t+1)$.
%[text] This teacher-forcing setup is equivalent to maximizing the log-likelihood $L = -\\sum\_t \\log P(c(t+1) \\mid c(1) \\ldots c(t))$.
%[text] In code, `X_seqMat` holds the input tokens and `Y_seqMat` the targets (shifted by one). Sequences shorter than MAX\_SEQ\_LEN are padded with PAD; the loss function masks PAD positions so they contribute no gradient.
encodedSeqs = cell(N_MOLS, 1);
for k = 1:N_MOLS
    s    = char(smilesProc(k));                % Use L/R placeholder format for vocabulary reference
    s    = s(1:min(end, MAX_SEQ_LEN - 2));     % Truncate if necessary
    idxs = arrayfun(@(c) char2idx(char(c)), s);
    encodedSeqs{k} = [START_IDX, idxs, END_IDX];
end
seqLengths = cellfun(@numel, encodedSeqs);
MAX_T      = max(seqLengths);

%[text] Pad to MAX\_T. seqMatrix is an integer matrix of size $\[N \\times \\mathrm{MAX\\\_T}\]$.
seqMatrix = ones(N_MOLS, MAX_T, "single") * PAD_IDX;
for k = 1:N_MOLS
    L = seqLengths(k);
    seqMatrix(k, 1:L) = single(encodedSeqs{k});
end

%[text] Shift input X and target Y by one token to form the teacher-forcing format.
X_seqMat = seqMatrix(:, 1:end-1);   % [N x SEQ_LEN]  Input
Y_seqMat = seqMatrix(:, 2:end);     % [N x SEQ_LEN]  Target (shifted by one token)
SEQ_LEN  = size(X_seqMat, 2);

logInfo("Sequence matrix: %d x %d  (vocabulary=%d)", N_MOLS, SEQ_LEN, VOCAB_SIZE); %[output:07691b4a]
%[text] Split into training and validation sets at an 80/20 ratio (fix random seed to 42 for reproducibility).
rng(42);
perm    = randperm(N_MOLS);
n_train = round(0.8 * N_MOLS);
trainIdx = perm(1:n_train);
valIdx   = perm(n_train+1:end);

X_train = X_seqMat(trainIdx, :);
Y_train = Y_seqMat(trainIdx, :);
X_val   = X_seqMat(valIdx, :);
Y_val   = Y_seqMat(valIdx, :);

logInfo("Split: %d training / %d validation", n_train, N_MOLS - n_train); %[output:0c181dc6]
%%
%[text] ## Section 7: Character-level LSTM -- Architecture and Training
%[text] ### BPTT (Backpropagation Through Time)
%[text] The LSTM reads sequences left to right, accumulating context in hidden state $h(t)$ and cell state $c(t)$.
%[text] During training, BPTT unrolls the LSTM over $T$ time steps and propagates gradients back through all of them.
%[text] The architecture:
%[text] 1. Input: one-hot tensor $\[\\mathrm{VOCAB\\\_SIZE} \\times T \\times B\]$
%[text] 2. Embedding layer (VOCAB\_SIZE → EMBED\_DIM): a learned character embedding that places similar tokens (C/N/O, c/n/o, etc.) near each other in a dense feature space.
%[text] 3. lstmLayer(HIDDEN\_SIZE): outputs $h(t)$ at each step $t$.
%[text] 4. dropoutLayer(DROPOUT\_RATE): regularization.
%[text] 5. fullyConnectedLayer(VOCAB\_SIZE) + softmaxLayer: outputs the next-token probability distribution $\[\\mathrm{VOCAB\\\_SIZE} \\times T \\times B\]$. \
%[text] The loss is masked cross-entropy, excluding PAD positions.
%[text] With only ~160 training samples (80% of a 200-molecule corpus), the train-validation loss gap is minimal (~0.06 at 150 epochs).
%[text] FDA-approved drugs span highly heterogeneous ring systems and molecular weights from 3 to 2000+, so the LSTM converges toward a uniform character distribution rather than memorizing specific structures.
%[text] This is underfitting, not overfitting.
%[text] Key lesson: even when the parameter/token ratio exceeds 1 ("overparameterization"), overfitting does not automatically follow if the data distribution is highly diverse.
%[text] 
%[text] Architecture: one-hot(V) → Embed(EMBED\_DIM) → LSTM(HIDDEN) → Dropout → FC(V) → Softmax
%[text] The embedding fully connected layer (VOCAB\_SIZE → EMBED\_DIM) functions as a learned character embedding.
%[text] It groups similar characters (atoms C/N/O, aromatic c/n/o, bond characters) in a dense feature space, improving gradient flow by having the LSTM receive dense EMBED\_DIM vectors instead of sparse 37-dimensional one-hot vectors.
%[text] As a result, 'C' and 'c' are learned as similar representations as "both carbon."
layers = [
    sequenceInputLayer(VOCAB_SIZE, Name="input")
    fullyConnectedLayer(EMBED_DIM, Name="embed")    % Learned character embedding (no activation)
    lstmLayer(HIDDEN_SIZE, OutputMode="sequence", Name="lstm1")
    dropoutLayer(DROPOUT_RATE, Name="drop")
    fullyConnectedLayer(VOCAB_SIZE, Name="fc")
    softmaxLayer(Name="softmax")
];
net = dlnetwork(layers);

nParams = sum(cellfun(@numel, net.Learnables.Value));
logInfo("Network: Embed(%d->%d) + LSTM(%d) + FC -> softmax  |  Number of parameters: %d", ... %[output:group:43ceed7d] %[output:3b3443ea]
    VOCAB_SIZE, EMBED_DIM, HIDDEN_SIZE, nParams); %[output:group:43ceed7d] %[output:3b3443ea]
logInfo("Training token to parameter ratio: %.2f  (>1 => overparameterization)", ... %[output:group:81631103] %[output:5f85d3e4]
    nParams / (n_train * mean(seqLengths))); %[output:group:81631103] %[output:5f85d3e4]
%[text] Execute custom training loop. Shuffle mini-batches each epoch and update with Adam optimizer.
numBatches   = ceil(n_train / BATCH_SIZE);
avgG_sl      = [];
avgSqG_sl    = [];
iter_sl      = 0;
trainLossLog = zeros(N_EPOCHS, 1);
valLossLog   = zeros(N_EPOCHS, 1);

for epoch = 1:N_EPOCHS %[output:group:914dc138]
    % Shuffle training order each epoch
    perm_e = randperm(n_train);
    X_shuf = X_train(perm_e, :);
    Y_shuf = Y_train(perm_e, :);
    epochLoss = 0;

    for b = 1:numBatches
        i1 = (b-1)*BATCH_SIZE + 1;
        i2 = min(b*BATCH_SIZE, n_train);
        Xb = X_shuf(i1:i2, :);   % [B x T]
        Yb = Y_shuf(i1:i2, :);   % [B x T]

        % One-hot encode: [VOCAB_SIZE x T x B]
        dlX = dlarray(onehot_encode(Xb, VOCAB_SIZE), "CTB");

        iter_sl = iter_sl + 1;
        [lossVal, grads] = dlfeval(@modelLoss, net, dlX, Yb, PAD_IDX);
        % Gradient clipping: Limit L2 norm per tensor to prevent loss spikes
        for gi = 1:height(grads)
            gdata = extractdata(grads.Value{gi});
            nrm = sqrt(sum(gdata(:).^2));
            if nrm > GRAD_CLIP_NORM
                grads.Value{gi} = grads.Value{gi} * (GRAD_CLIP_NORM / nrm);
            end
        end
        epochLoss = epochLoss + extractdata(lossVal);
        % Cosine LR annealing: Decay from LEARN_RATE to LEARN_RATE/20 over N_EPOCHS
        lr_epoch = LEARN_RATE * (0.05 + 0.95 * 0.5 * (1 + cos(pi * (epoch-1) / N_EPOCHS)));
        [net, avgG_sl, avgSqG_sl] = adamupdate(net, grads, avgG_sl, avgSqG_sl, ...
            iter_sl, lr_epoch);
    end
    trainLossLog(epoch) = epochLoss / numBatches;

    % Validation loss (no gradients)
    dlXv = dlarray(onehot_encode(X_val, VOCAB_SIZE), "CTB");
    valLossLog(epoch) = extractdata(dlfeval(@modelLoss, net, dlXv, Y_val, PAD_IDX));

    if mod(epoch, 10) == 0 || epoch == 1
        logInfo("Epoch %3d/%d  train=%.4f  val=%.4f", ... %[output:81098cf3]
            epoch, N_EPOCHS, trainLossLog(epoch), valLossLog(epoch)); %[output:81098cf3]
    end
end %[output:group:914dc138]
figure("Name", "R05 Training Curve"); %[output:84c89218]
plot(1:N_EPOCHS, trainLossLog, "b-", LineWidth=1.5, DisplayName="Training"); %[output:84c89218]
hold on; %[output:84c89218]
plot(1:N_EPOCHS, valLossLog,   "r--", LineWidth=1.5, DisplayName="Validation"); %[output:84c89218]
xlabel("Epoch"); ylabel("Cross-Entropy Loss"); %[output:84c89218]
title("Character-level LSTM: val loss plateau as key evidence of data scale limit"); %[output:84c89218]
legend; grid on; %[output:84c89218]
%[text] The "plateau" where val loss stops improving around epochs 50-70 is a key learning point.
%[text] The numerical difference in efficacy rates (Markov about 9%, LSTM about 5%) is about 3-4 points at $N=300$ samples and is not statistically significant ($z=1.6$, $p \\approx 0.10$, 95% CI overlap).
%[text] Avoid overinterpreting the raw efficacy rate gap and judge by the loss curve.
logInfo("Training complete. Final train/val: %.4f / %.4f", ... %[output:group:3700040c] %[output:7a5fd55e]
    trainLossLog(end), valLossLog(end)); %[output:group:3700040c] %[output:7a5fd55e]
%%
%[text] ## Section 8: SMILES Generation by Temperature Sampling
%[text] ### Autoregressive Decoding and Temperature Sampling
%[text] After training, new SMILES are generated autoregressively one character at a time.
%[text] 1. Input the START token and perform the forward pass of the LSTM.
%[text] 2. Sample the next character from the output distribution.
%[text] 3. Use the sampled character as the next input.
%[text] 4. Repeat (1) to (3) until sampling the END token or reaching MAX\_SEQ\_LEN. \
%[text] Temperature $T$ adjusts the balance between exploration and exploitation: $P\_T(c) = \\mathrm{softmax}\\!\\left(\\log P(c) / T\\right)$
%[text] For $T \< 1$, the distribution sharpens, reducing invalid characters with conservative sampling.
%[text] For $T \> 1$, the distribution flattens, increasing creativity but also invalid SMILES.
%[text] Implementation note: Input the entire growing prefix to the LSTM at each generation step (not just the last token).
%[text] This results in $O(T^2)$ computation but is accurate (can be $O(T)$ with stateful LSTM).
%[text] **Grammar Constraints**: generateSmiles() applies 5 masking rules to eliminate invalid SMILES.
%[text] 1. Leading Position Mask: Only atom characters and '\[' can be the first token.
%[text] 2. Parenthesis Depth Guard: ')' is prohibited when paren\_depth = 0.
%[text] 3. Inside Parenthesis Mask: Prohibit '(' and '\[' inside '\[...\]'.
%[text] 4. END Block: Prohibit END while there are unclosed parentheses or ring bonds.
%[text] 5. Forced Closure: Force closure of open delimiters if remaining steps are insufficient. \
%[text] Structural issues like aromatic valency require a larger training corpus.
logInfo("%d SMILES being generated (T=%.2f) ...", N_GENERATED, TEMPERATURE); %[output:2d058d5c]
generatedSmiles = strings(N_GENERATED, 1);
for g = 1:N_GENERATED %[output:group:525e00c9]
    generatedSmiles(g) = generateSmiles(net, START_IDX, END_IDX, PAD_IDX, ...
        VOCAB_SIZE, MAX_SEQ_LEN, TEMPERATURE, idx2char, char2idx);
    if mod(g, max(1, round(N_GENERATED/10))) == 0 || g == N_GENERATED
        logProgress(g, N_GENERATED, "Generating (pre-RL)"); %[output:3bd50a65]
    end
end %[output:group:525e00c9]
logInfo("Sample of generated SMILES (first 8):"); %[output:74269a7c]
for k = 1:min(8, N_GENERATED) %[output:group:4081c67c]
    logInfo("  [%3d] %s", k, generatedSmiles(k)); %[output:42d269fd]
end %[output:group:4081c67c]
%%
%[text] ## Section 9: Efficacy, Diversity, Novelty, and Property Analysis
%[text] ### Four Standard Evaluation Metrics for Generative Models (Brown et al. 2019)
%[text] **Validity**: The percentage of generated SMILES that RDKit can parse. Reflects how well the model has learned SMILES grammar.
%[text] **Uniqueness**: The percentage of non-duplicate valid generated SMILES. Lower diversity indicates model collapse to a few structures.
%[text] **Novelty**: The percentage of non-duplicate valid SMILES not present in the training set. Lower novelty indicates memorization and lack of generalization.
%[text] **Property KL Divergence**: The KL divergence between property distributions of generated and training molecules. Lower KL indicates better reproduction of the training corpus "style".
%[text] Expected values for default settings (500 molecules, HIDDEN=128, ~80K parameters): Validity 5-15% (with grammar constraints), Diversity 60-90%, Novelty 80-100%. For 200 molecules, HIDDEN=64 (~28K parameters), Validity 5-20%, Novelty around 40-70%.
%[text] Low validity is expected, indicating data scale challenges.
%[text] Production models (REINVENT, MolGPT) achieve Validity \> 95% with over $10^6$ molecules.
validGenMask = false(N_GENERATED, 1);
for k = 1:N_GENERATED %[output:group:354fa23a]
    try
        validGenMask(k) = emk.mol.isValid(generatedSmiles(k));
    catch
        validGenMask(k) = false;
    end
    if mod(k, max(1, round(N_GENERATED/10))) == 0 || k == N_GENERATED
        logProgress(k, N_GENERATED, "Validating (Generated)"); %[output:4c71865b]
    end
end %[output:group:354fa23a]
validGenSmiles = generatedSmiles(validGenMask);
nValid   = sum(validGenMask);
nUnique  = numel(unique(validGenSmiles));
nNovel   = sum(~ismember(validGenSmiles, smilesOriginal));  % novelty vs canonical training set

logInfo("Pre-RL Generation -- Total %d", N_GENERATED); %[output:58037d9b]
logInfo("  Validity:    %d/%d  (%.1f%%)", nValid,  N_GENERATED, 100*nValid/N_GENERATED); %[output:68c1252a]
logInfo("  Uniqueness:  %d/%d  (%.1f%%)", nUnique, max(nValid,1), 100*nUnique/max(nValid,1)); %[output:1307a96b]
logInfo("  Novelty:     %d/%d  (%.1f%%)", nNovel,  max(nValid,1), 100*nNovel/max(nValid,1)); %[output:54f8ce50]
%[text] Compare property distributions of generated and training molecules (molecular weight and LogP).
FEAT_NAMES = ["MolWt", "LogP"];
propsGen   = batchDescriptors(validGenSmiles(1:min(50,nValid)), FEAT_NAMES);
propsTrain = batchDescriptors(smilesOriginal, FEAT_NAMES);  % property dist of original molecules

if ~isempty(propsGen) && ~isempty(propsTrain) %[output:group:1f959db9]
    figure("Name", "R05 Property Distribution (pre-RL)"); %[output:9f971cbd]
    titles = ["Molecular Weight", "LogP (Lipophilicity)"];
    units  = ["Da", "log units"];
    for fi = 1:numel(FEAT_NAMES)
        subplot(1, 2, fi); hold on; %[output:9f971cbd]
        histogram(propsTrain(:,fi), 20, Normalization="probability", ...
            FaceColor=[0.2 0.6 0.8], DisplayName="Training");
        histogram(propsGen(:,fi),   20, Normalization="probability", ...
            FaceColor=[0.9 0.4 0.2], FaceAlpha=0.7, DisplayName="Generated");
        xlabel(sprintf("%s (%s)", FEAT_NAMES(fi), units(fi)));
        ylabel("Probability Density");
        title(titles(fi)); legend; grid on;
    end
    sgtitle("Generated vs Training: Property Distribution (Pre-RL)"); %[output:9f971cbd]
end %[output:group:1f959db9]
%%
%[text] ## Section 10: Checkpoint Saving and Summary
%[text] ### Serialization of Trained Model
%[text] By saving the trained network and vocabulary, you can resume generation in a new session without retraining.
%[text] You can also transfer it to r06\_reinforce.m for REINFORCE fine-tuning.
%[text] Save to result/r05\_checkpoint.mat (result/ is managed by .gitignore, so the file will not be committed to the repository).
CHECKPOINT = "result/r05_checkpoint.mat";
if ~isfolder("result"); mkdir("result"); end
save(CHECKPOINT, ...
    "net", "char2idx", "idx2char", "VOCAB_SIZE", "MAX_SEQ_LEN", ...
    "START_IDX", "END_IDX", "PAD_IDX", "HIDDEN_SIZE", "TEMPERATURE", ...
    "smilesAll", "smilesProc", "X_train", "Y_train", "X_val", "Y_val");
logInfo("Checkpoint saved -> %s", CHECKPOINT); %[output:6816dab6]
logInfo("  Load in R06: load(""result/r05_checkpoint.mat"")"); %[output:4dca65b3]
%[text] ### --- Summary ---
logInfo("Corpus:         %d molecules Valid FDA drug SMILES (trained with %d molecules after expansion)", N_MOLS_ORIG, N_MOLS); %[output:4752e9cd]
logInfo("Vocabulary:     %d characters  |  Max sequence length: %d", VOCAB_SIZE, MAX_T); %[output:86afdc49]
logInfo("Loss plateau:   val=%.4f (epoch %d+)  --  Random baseline=%.4f  (Perplexity=%.1f)", ... %[output:group:44d2b1de] %[output:5038645d]
    min(valLossLog), find(valLossLog == min(valLossLog), 1), ... %[output:5038645d]
    log(VOCAB_SIZE), exp(min(valLossLog))); %[output:group:44d2b1de] %[output:5038645d]
logInfo("--- Generation Metrics (N=%d each; Binomial SE~%.1f%%)", N_GENERATED, ... %[output:group:22e3defc] %[output:134573c2]
    100*sqrt(0.10*0.90/N_GENERATED)); %[output:group:22e3defc] %[output:134573c2]
if exist("nMValid", "var") %[output:group:4103d63f]
    logInfo("Markov Model:   Validity=%.1f%%  Diversity=%.1f%%  Novelty=%.1f%%", ... %[output:6104ce8f]
        100*nMValid/N_MARKOV, 100*nMUnique/max(nMValid,1), 100*nMNovel/max(nMValid,1)); %[output:6104ce8f]
else
    logInfo("Markov Model:   (Please run Sections 3-5 first)");
end %[output:group:4103d63f]
logInfo("LSTM Model:     Validity=%.1f%%  Diversity=%.1f%%  Novelty=%.1f%%", ... %[output:group:08457e7b] %[output:97f77c7d]
    100*nValid/N_GENERATED, 100*nUnique/max(nValid,1), 100*nNovel/max(nValid,1)); %[output:group:08457e7b] %[output:97f77c7d]
logInfo("Warning: With N=%d, the validity gap is not statistically significant (p~0.10)", N_GENERATED); %[output:7a8ab334]
logInfo("Architecture:   LSTM(%d) -> Dropout(%.2f) -> FC(%d) -> Softmax", ... %[output:group:650dc90a] %[output:198ebf94]
    HIDDEN_SIZE, DROPOUT_RATE, VOCAB_SIZE); %[output:group:650dc90a] %[output:198ebf94]
logInfo("Number of Parameters:     %d  (Overparameterized relative to %d training tokens)", ... %[output:group:3f37985d] %[output:17de43a5]
    nParams, n_train * round(mean(seqLengths))); %[output:group:3f37985d] %[output:17de43a5]
%[text] ### Key Lessons
%[text] 1. **Main Evidence**: The val loss plateaus early. LSTM cannot learn SMILES grammar with a few hundred molecules. It's a data scale bottleneck, not a code bug.
%[text] 2. **Validity Gap** (Markov vs LSTM) is a reference value but not statistically significant at $N=300$. Judge by loss curves, not validity numbers.
%[text] 3. **Inductive Bias** (Markov's local transitions) surpasses expressiveness with small $N$.
%[text] 4. Production models (REINVENT, MolGPT) use $10^6 \\sim 10^7$ molecules.
%[text] 5. The next step is r06\_reinforce.m. Learn how to apply REINFORCE policy gradient RL to guide molecule generation towards target properties. \
%[text] \--- Local Functions ---
function smiles = generateSmilesMarkov(transProb, startIdx, endIdx, padIdx, ...
    vocabSize, maxLen, temperature, idx2char, char2idx)
% BIGRAM MARKOV CHAIN SMILES generation with grammar constraints.
% transProb: [vocabSize x vocabSize] row-normalised transition probability matrix.
% Each row transProb(i, :) is the distribution over the next token
% given that the current token is i.
% All 6 grammar constraint rules are identical to generateSmiles() --
% only the next-token distribution differs (table lookup vs LSTM forward pass).
tokens     = startIdx;
currentTok = startIdx;

% Pre-lookup grammar-critical indices (0 = absent in vocab).
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

% Aromatic wildcard 'a' cannot be generated correctly without ring-aromaticity
% context tracking.  Pre-compute its index so it can be suppressed each step.
aWildcardIdx = 0;
if isKey(char2idx, 'a'); aWildcardIdx = char2idx('a'); end

% Maximum simultaneous open ring bonds.  Drug SMILES rarely exceed 3;
% without this limit the Markov model accumulates 5-6 open ring bonds
% that cannot all close before maxLen, producing 0% validity.
MAX_OPEN_RINGS = 3;

% Additional context-sensitive token indices (Rules 2d-2h)
cAtomIdx    = 0; if isKey(char2idx, 'C'); cAtomIdx    = char2idx('C'); end
bAtomIdx    = 0; if isKey(char2idx, 'B'); bAtomIdx    = char2idx('B'); end
lSuffixIdx  = 0; if isKey(char2idx, 'l'); lSuffixIdx  = char2idx('l'); end
rSuffixIdx  = 0; if isKey(char2idx, 'r'); rSuffixIdx  = char2idx('r'); end
stereoAtIdx = 0; if isKey(char2idx, '@'); stereoAtIdx = char2idx('@'); end
dotIdx      = 0; if isKey(char2idx, '.'); dotIdx      = char2idx('.'); end
fwdSlashIdx = 0; if isKey(char2idx, '/'); fwdSlashIdx = char2idx('/'); end
bkSlashIdx  = 0; if isKey(char2idx, char(92)); bkSlashIdx  = char2idx(char(92)); end
plusOutIdx  = 0; if isKey(char2idx, '+'); plusOutIdx  = char2idx('+'); end
eqBondIdx   = 0; if isKey(char2idx, '='); eqBondIdx   = char2idx('='); end
trplBondIdx = 0; if isKey(char2idx, '#'); trplBondIdx = char2idx('#'); end

paren_depth    = 0;
bracket_depth  = 0;
inside_bracket = false;
ring_open      = false(1, 9);

for t = 1:maxLen - 1
    % Bigram table lookup: distribution over next token given current token.
    prob_t = double(transProb(currentTok, :)');   % [vocabSize x 1]

    % Temperature scaling (same formula as generateSmiles)
    prob_t = max(prob_t, 1e-9);
    logit  = log(prob_t) / temperature;
    logit  = logit - max(logit);
    prob_t = exp(logit);
    prob_t = prob_t / sum(prob_t);

    % ---- Grammar constraints (same 6 rules as generateSmiles) -----------
    % 1. First position: restrict to valid SMILES-start characters
    if t == 1
        prob_t = prob_t .* firstCharMask;
    end

    % 2. Closing tokens forbidden when no matching opener exists
    if parenCloseIdx   > 0 && paren_depth   == 0; prob_t(parenCloseIdx)   = 0; end
    if bracketCloseIdx > 0 && bracket_depth == 0; prob_t(bracketCloseIdx) = 0; end

    % 2b. Prevent empty delimiters () and []
    if numel(tokens) > 1
        lastTok = tokens(end);
        if parenCloseIdx   > 0 && parenOpenIdx   > 0 && lastTok == parenOpenIdx
            prob_t(parenCloseIdx)   = 0;
        end
        if bracketCloseIdx > 0 && bracketOpenIdx > 0 && lastTok == bracketOpenIdx
            prob_t(bracketCloseIdx) = 0;
        end
    end

    % 2c. Ring bond depth guard: block NEW ring-opening digits when too many
    %     ring bonds are already open.  Drug SMILES rarely need more than 3
    %     simultaneous ring bonds; without this guard the Markov model
    %     accumulates open ring bonds that cannot close before maxLen.
    if sum(ring_open) >= MAX_OPEN_RINGS
        for dk_guard = 1:9
            if digitIdx(dk_guard) > 0 && ~ring_open(dk_guard)
                prob_t(digitIdx(dk_guard)) = 0;  % block only NEW opens
            end
        end
    end

    % 2d. Multi-char halogen suffix: 'l' (2nd char of Cl) valid ONLY after 'C';
    %     'r' (2nd char of Br) valid ONLY after 'B'.  All other placements create
    %     unrecognised atom symbols that RDKit rejects immediately.
    prevTok = tokens(end);
    if lSuffixIdx > 0 && prevTok ~= cAtomIdx; prob_t(lSuffixIdx) = 0; end
    if rSuffixIdx > 0 && prevTok ~= bAtomIdx; prob_t(rSuffixIdx) = 0; end

    % 2e. '@' (stereochemistry) is valid ONLY inside '[ ]' atom blocks.
    if ~inside_bracket && stereoAtIdx > 0; prob_t(stereoAtIdx) = 0; end

    % 2f. '.' creates disconnected fragments; target is single connected molecules.
    if dotIdx > 0; prob_t(dotIdx) = 0; end

    % 2g. '/' and '\' (E/Z bond direction) require adjacent '=' context; block both.
    if fwdSlashIdx > 0; prob_t(fwdSlashIdx) = 0; end
    if bkSlashIdx  > 0; prob_t(bkSlashIdx)  = 0; end

    % 2h. '+' (charge) outside '[ ]' is a SMILES syntax error.
    if ~inside_bracket && plusOutIdx > 0; prob_t(plusOutIdx) = 0; end

    % 3. Inside bracket: '(' and '[' are not valid (brackets do not nest).
    %    CRITICAL: block ring closure digits inside brackets.  Without this,
    %    the model generates '[N1]'-style invalid bracket atoms AND the
    %    ring_open state tracker toggles incorrectly, causing systematic
    %    ring mismatch that makes every subsequent ring closure invalid.
    %    Also block bond descriptors '=' and '#' (not part of bracket atom syntax).
    if inside_bracket
        if parenOpenIdx   > 0; prob_t(parenOpenIdx)   = 0; end
        if bracketOpenIdx > 0; prob_t(bracketOpenIdx) = 0; end
        for dk3 = 1:9
            if digitIdx(dk3) > 0; prob_t(digitIdx(dk3)) = 0; end
        end
        if eqBondIdx   > 0; prob_t(eqBondIdx)   = 0; end
        if trplBondIdx > 0; prob_t(trplBondIdx) = 0; end
    end

    % 4. Block END while any delimiter or ring bond is unclosed.
    if bracket_depth > 0 || paren_depth > 0 || any(ring_open)
        prob_t(endIdx) = 0;
    end

    % 5. Force-close open delimiters / ring bonds when insufficient steps remain.
    %    nOpenTokens = number of closing tokens still needed (1 per open bracket/paren/ring).
    %    If remaining steps <= nOpenTokens, start force-closing immediately.
    %    Replaces hardcoded maxLen-3 budget which fails when multiple rings are open.
    nOpenTokens = bracket_depth + paren_depth + sum(ring_open);
    if (maxLen - 1 - t) <= nOpenTokens
        if bracket_depth > 0 && bracketCloseIdx > 0
            prob_t(:) = 0;  prob_t(bracketCloseIdx) = 1;
        elseif paren_depth > 0 && parenCloseIdx > 0
            prob_t(:) = 0;  prob_t(parenCloseIdx) = 1;
        elseif any(ring_open)
            open_d = find(ring_open, 1);
            if digitIdx(open_d) > 0
                prob_t(:) = 0;  prob_t(digitIdx(open_d)) = 1;
            end
        end
    end

    % 6. Redirect aromatic atom probability to aliphatic equivalents.
    %    Also block aromatic wildcard 'a': valid only in SMARTS / specific
    %    ring-aromatic contexts that the Markov model cannot track.
    for ak = 1:nAromPairs
        prob_t(aliphaticIdxArr(ak)) = prob_t(aliphaticIdxArr(ak)) + prob_t(aromaticIdxArr(ak));
        prob_t(aromaticIdxArr(ak))  = 0;
    end
    if aWildcardIdx > 0; prob_t(aWildcardIdx) = 0; end

    % Renormalise; fallback to uniform over non-special tokens if all masked.
    % NOTE: re-apply END block AFTER fallback -- the fallback resets prob_t(:)=1
    % which would undo Rule 4 and allow END to be sampled inside open brackets/parens.
    s = sum(prob_t);
    if s <= 0
        prob_t(:) = 1;
        prob_t(padIdx)   = 0;
        prob_t(startIdx) = 0;
        % Re-apply END block so fallback cannot override grammar constraint
        if bracket_depth > 0 || paren_depth > 0 || any(ring_open)
            prob_t(endIdx) = 0;
        end
        s = sum(prob_t);
        if s > 0; prob_t = prob_t / s; end
    else
        prob_t = prob_t / s;
    end

    nextIdx = randsample(vocabSize, 1, true, prob_t);
    tokens(end+1) = nextIdx; %#ok<AGROW>
    currentTok = nextIdx;

    % Update grammar state
    if     parenOpenIdx    > 0 && nextIdx == parenOpenIdx
        paren_depth = paren_depth + 1;
    elseif parenCloseIdx   > 0 && nextIdx == parenCloseIdx
        paren_depth = paren_depth - 1;
    elseif bracketOpenIdx  > 0 && nextIdx == bracketOpenIdx
        bracket_depth  = bracket_depth + 1;
        inside_bracket = true;
    elseif bracketCloseIdx > 0 && nextIdx == bracketCloseIdx
        bracket_depth  = bracket_depth - 1;
        inside_bracket = (bracket_depth > 0);
    end
    % Ring closure state: toggle ONLY outside brackets.
    % Inside brackets, digits are part of isotope/H-count notation, NOT ring closures.
    % Toggling inside brackets would corrupt ring_open and cause every subsequent
    % ring closure to fail -- the dominant cause of 0% Markov validity.
    if ~inside_bracket
        for dk = 1:9
            if digitIdx(dk) > 0 && nextIdx == digitIdx(dk)
                ring_open(dk) = ~ring_open(dk); break;
            end
        end
    end

    if nextIdx == endIdx; break; end
end

% Decode: skip START; stop at END/PAD; restore Cl/Br placeholders.
charList = {};
for i = 2:numel(tokens)
    idx = tokens(i);
    if idx == endIdx || idx == padIdx; break; end
    if idx == startIdx; continue; end
    if isKey(idx2char, idx)
        charList{end+1} = idx2char(idx); %#ok<AGROW>
    end
end
smiles = string(strjoin(charList, ""));
smiles = strrep(strrep(smiles, "L", "Cl"), "R", "Br");
end

function Xoh = onehot_encode(seqMat, vocabSize)
% ONE-HOT ENCODE a [N x T] integer matrix to single [vocabSize x T x N].
% Vectorised via linear indexing: no nested loops, scales to large N.
[N, T] = size(seqMat);
Xoh = zeros(vocabSize, T, N, "single");
valid = seqMat >= 1 & seqMat <= vocabSize;   % [N x T] logical mask
[nz_n, nz_t] = find(valid);                  % row (N) and col (T) of valid entries
nz_c = double(seqMat(valid));                % character indices (column-major order)
% Linear index into [vocabSize x T x N]: Xoh(c, t, n)
linIdx = nz_c + (nz_t - 1) * vocabSize + (nz_n - 1) * vocabSize * T;
Xoh(linIdx) = 1;
end

function [loss, grads] = modelLoss(net, dlX, Ytarget, padIdx)
% MASKED CROSS-ENTROPY LOSS for sequence-to-sequence training.
% dlX:     dlarray [VOCAB_SIZE x T x B] 'CTB'
% Ytarget: single  [B x T]  integer target token indices
% padIdx:  scalar  PAD token index to mask
pred   = forward(net, dlX);          % [V x T x B] 'CTB'
V      = size(pred, 1);
T      = size(pred, 2);
B      = size(pred, 3);
predR  = reshape(stripdims(pred), V, T*B);   % [V x T*B]

% Target: [B x T] -> [T x B] column-major -> [1 x T*B]
Yt = single(reshape(permute(Ytarget, [2 1]), 1, T*B));

% One-hot target matrix [V x T*B] -- vectorised via linear indexing (same as onehot_encode).
% Avoids T*B loop iterations (~3000 per batch) that dominated modelLoss runtime.
targetOH  = zeros(V, T*B, "single");
valid_t   = Yt >= 1 & Yt <= V;              % [1 x T*B] logical mask
validCols = find(valid_t);                   % column indices of valid tokens
validRows = double(Yt(valid_t));             % row indices (character indices)
if ~isempty(validCols)
    linIdx_t = validRows(:) + (double(validCols(:)) - 1) * V;
    targetOH(linIdx_t) = 1;
end
dlTarget = dlarray(targetOH);

% Mask: 1 where target != PAD
mask    = dlarray(single(Yt ~= padIdx));   % [1 x T*B]
nTokens = max(sum(extractdata(mask)), 1);

% Cross entropy: -sum_v target_OH(v,i) * log(pred(v,i)), summed per pos
logPred = log(predR + 1e-8);                       % [V x T*B]
ce      = -sum(dlTarget .* logPred, 1);            % [1 x T*B]
loss    = sum(ce .* mask) / nTokens;               % scalar

grads = dlgradient(loss, net.Learnables);
end

function smiles = generateSmiles(net, startIdx, endIdx, padIdx, ...
    vocabSize, maxLen, temperature, idx2char, char2idx)
% AUTOREGRESSIVE SMILES GENERATION with temperature sampling and extended
% grammar constraints.  Feeds growing prefix to LSTM; O(T^2) but correct.
% 
% Grammar constraints (6 rules -- see Section 5 concept note for details):
% 1. First-position mask  2. Paren/bracket depth  3. Ring-depth guard
% 4. Inside-bracket mask  5. END blocked while delimiters open
% 6. Force-close near maxLen
tokens = startIdx;

% Pre-lookup grammar-critical token indices (0 = absent in vocabulary).
parenOpenIdx    = 0;  parenCloseIdx   = 0;
bracketOpenIdx  = 0;  bracketCloseIdx = 0;
if isKey(char2idx, '(');  parenOpenIdx    = char2idx('(');  end
if isKey(char2idx, ')');  parenCloseIdx   = char2idx(')');  end
if isKey(char2idx, '[');  bracketOpenIdx  = char2idx('[');  end
if isKey(char2idx, ']');  bracketCloseIdx = char2idx(']');  end

% Pre-lookup ring closure digit indices (digits '1'-'9').
digitIdx = zeros(1, 9);  % digitIdx(d) = vocab index for char d (1=>'1', ..., 9=>'9')
for dk = 1:9
    dc = char('0' + dk);
    if isKey(char2idx, dc); digitIdx(dk) = char2idx(dc); end
end

% Valid SMILES-start character mask: aliphatic atom letters + '['.
% Aromatic atoms excluded: Rule 6 redirects them to aliphatic equivalents,
% so they must also be excluded from valid first-position characters.
% Digits, '/', '\', '=', '#', '.', '(', ')', '@', '+', '-' are illegal starts.
VALID_FIRST = {'B','C','N','O','P','S','F','I','L','R','[','H'};
firstCharMask = zeros(vocabSize, 1);
for fk = 1:numel(VALID_FIRST)
    if isKey(char2idx, VALID_FIRST{fk})
        firstCharMask(char2idx(VALID_FIRST{fk})) = 1;
    end
end

% Pre-compute aromatic -> aliphatic index redirects (Rule 6).
% Aromatic atoms (c, n, o, s, p, b) are only valid inside aromatic ring closures.
% Tracking aromatic ring context is complex; instead we redirect their probability
% mass to the aliphatic equivalent at generation time.  This keeps ring-notation
% SMILES generatable (C1CCCCC1 for cyclohexane) while avoiding illegal 'c' outside rings.
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

% Aromatic wildcard 'a' -- suppress in generation (same reason as Markov).
aWildcardIdx = 0;
if isKey(char2idx, 'a'); aWildcardIdx = char2idx('a'); end

% Maximum simultaneous open ring bonds (same constant as generateSmilesMarkov).
MAX_OPEN_RINGS = 3;

% Additional context-sensitive token indices (Rules 2d-2h, same as Markov)
cAtomIdx    = 0; if isKey(char2idx, 'C'); cAtomIdx    = char2idx('C'); end
bAtomIdx    = 0; if isKey(char2idx, 'B'); bAtomIdx    = char2idx('B'); end
lSuffixIdx  = 0; if isKey(char2idx, 'l'); lSuffixIdx  = char2idx('l'); end
rSuffixIdx  = 0; if isKey(char2idx, 'r'); rSuffixIdx  = char2idx('r'); end
stereoAtIdx = 0; if isKey(char2idx, '@'); stereoAtIdx = char2idx('@'); end
dotIdx      = 0; if isKey(char2idx, '.'); dotIdx      = char2idx('.'); end
fwdSlashIdx = 0; if isKey(char2idx, '/'); fwdSlashIdx = char2idx('/'); end
bkSlashIdx  = 0; if isKey(char2idx, char(92)); bkSlashIdx  = char2idx(char(92)); end
plusOutIdx  = 0; if isKey(char2idx, '+'); plusOutIdx  = char2idx('+'); end
eqBondIdx   = 0; if isKey(char2idx, '='); eqBondIdx   = char2idx('='); end
trplBondIdx = 0; if isKey(char2idx, '#'); trplBondIdx = char2idx('#'); end

% Grammar state
paren_depth    = 0;      % unclosed '(' count
bracket_depth  = 0;      % unclosed '[' count (0 or 1 in valid SMILES)
inside_bracket = false;  % true when currently inside a '[ ]' atom block
ring_open      = false(1, 9);  % ring_open(d): true when digit d has been seen an odd
                               % number of times (ring bond opened but not yet closed)

for t = 1:maxLen - 1
    T_curr = numel(tokens);
    % Build one-hot prefix: [V x T_curr x 1]
    xenc = zeros(vocabSize, T_curr, 1, "single");
    for i = 1:T_curr
        c = tokens(i);
        if c >= 1 && c <= vocabSize
            xenc(c, i, 1) = 1;
        end
    end
    dlX  = dlarray(xenc, "CTB");
    pred = predict(net, dlX);                       % [V x T_curr x 1]
    prob = double(extractdata(pred(:, end, 1)));    % last time step: [V x 1]
    prob = max(prob, 1e-9);

    % Temperature sampling
    logit  = log(prob) / temperature;
    logit  = logit - max(logit);                    % numerical stability
    prob_t = exp(logit);
    prob_t = prob_t / sum(prob_t);

    % ---- Grammar constraints (applied in priority order) ---------------
    % 1. First-position: restrict to valid SMILES-start characters only
    if t == 1
        prob_t = prob_t .* firstCharMask;
    end

    % 2. Closing tokens forbidden when no matching opener exists
    if parenCloseIdx   > 0 && paren_depth   == 0; prob_t(parenCloseIdx)   = 0; end
    if bracketCloseIdx > 0 && bracket_depth == 0; prob_t(bracketCloseIdx) = 0; end

    % 2b. Prevent empty delimiters: ')' immediately after '(' or ']' after '['
    %     Empty branches () and empty brackets [] are never valid SMILES.
    if numel(tokens) > 1
        lastTok = tokens(end);
        if parenCloseIdx   > 0 && parenOpenIdx   > 0 && lastTok == parenOpenIdx
            prob_t(parenCloseIdx)   = 0;
        end
        if bracketCloseIdx > 0 && bracketOpenIdx > 0 && lastTok == bracketOpenIdx
            prob_t(bracketCloseIdx) = 0;
        end
    end

    % 2c. Ring bond depth guard: block NEW ring-opening digits when too many
    %     are already open.  Limits simultaneous open ring bonds to MAX_OPEN_RINGS.
    if sum(ring_open) >= MAX_OPEN_RINGS
        for dk_guard = 1:9
            if digitIdx(dk_guard) > 0 && ~ring_open(dk_guard)
                prob_t(digitIdx(dk_guard)) = 0;  % block only NEW opens
            end
        end
    end

    % 2d. Multi-char halogen suffix: 'l' valid ONLY after 'C'; 'r' ONLY after 'B'.
    prevTok = tokens(end);
    if lSuffixIdx > 0 && prevTok ~= cAtomIdx; prob_t(lSuffixIdx) = 0; end
    if rSuffixIdx > 0 && prevTok ~= bAtomIdx; prob_t(rSuffixIdx) = 0; end

    % 2e. '@' (stereochemistry) valid ONLY inside '[ ]' atom blocks.
    if ~inside_bracket && stereoAtIdx > 0; prob_t(stereoAtIdx) = 0; end

    % 2f. '.' (disconnect) blocked; target is single connected molecules.
    if dotIdx > 0; prob_t(dotIdx) = 0; end

    % 2g. '/' and '\' (E/Z bond direction) require adjacent '=' context; block both.
    if fwdSlashIdx > 0; prob_t(fwdSlashIdx) = 0; end
    if bkSlashIdx  > 0; prob_t(bkSlashIdx)  = 0; end

    % 2h. '+' (charge) outside '[ ]' is invalid SMILES syntax.
    if ~inside_bracket && plusOutIdx > 0; prob_t(plusOutIdx) = 0; end

    % 3. Inside bracket: several tokens are invalid in bracket atom content.
    %    CRITICAL: block ring closure digits inside brackets.  Without this,
    %    the model generates '[N1]'-style invalid bracket atoms AND the
    %    ring_open state tracker toggles incorrectly, causing systematic
    %    ring mismatch that makes every subsequent ring closure invalid.
    %    Also block bond descriptors '=' and '#' (not part of bracket atom syntax).
    if inside_bracket
        if parenOpenIdx   > 0; prob_t(parenOpenIdx)   = 0; end
        if bracketOpenIdx > 0; prob_t(bracketOpenIdx) = 0; end
        for dk3 = 1:9
            if digitIdx(dk3) > 0; prob_t(digitIdx(dk3)) = 0; end
        end
        if eqBondIdx   > 0; prob_t(eqBondIdx)   = 0; end
        if trplBondIdx > 0; prob_t(trplBondIdx) = 0; end
    end

    % 4. Block END while any delimiter or ring bond is unclosed.
    if bracket_depth > 0 || paren_depth > 0 || any(ring_open)
        prob_t(endIdx) = 0;
    end

    % 5. Force-close open delimiters / ring bonds when insufficient steps remain.
    %    nOpenTokens = number of closing tokens still needed (1 per open bracket/paren/ring).
    %    If remaining steps <= nOpenTokens, start force-closing immediately.
    %    Replaces hardcoded maxLen-3 budget which fails when multiple rings are open.
    nOpenTokens = bracket_depth + paren_depth + sum(ring_open);
    if (maxLen - 1 - t) <= nOpenTokens
        if bracket_depth > 0 && bracketCloseIdx > 0
            prob_t(:) = 0;  prob_t(bracketCloseIdx) = 1;
        elseif paren_depth > 0 && parenCloseIdx > 0
            prob_t(:) = 0;  prob_t(parenCloseIdx) = 1;
        elseif any(ring_open)
            open_d = find(ring_open, 1);
            if digitIdx(open_d) > 0
                prob_t(:) = 0;  prob_t(digitIdx(open_d)) = 1;
            end
        end
    end

    % 6. Redirect aromatic atom probability to aliphatic equivalents.
    %    Aromatic atoms without a closed aromatic ring produce invalid SMILES.
    %    Transferring probability mass (c->C, n->N, ...) ensures the model's
    %    learned frequency for aromatic atoms still contributes to generation,
    %    just via their aliphatic equivalents.  Ring SMILES like C1CCCCC1 remain
    %    fully generatable; only the aromatic shorthand is suppressed.
    %    Also block 'a' (aromatic wildcard) -- same reason as generateSmilesMarkov.
    for ak = 1:nAromPairs
        prob_t(aliphaticIdxArr(ak)) = prob_t(aliphaticIdxArr(ak)) + prob_t(aromaticIdxArr(ak));
        prob_t(aromaticIdxArr(ak))  = 0;
    end
    if aWildcardIdx > 0; prob_t(aWildcardIdx) = 0; end

    % Renormalise; fallback to uniform over non-special tokens if all masked.
    % NOTE: re-apply END block AFTER fallback -- the fallback resets prob_t(:)=1
    % which would undo Rule 4 and allow END to be sampled inside open brackets/parens.
    s = sum(prob_t);
    if s <= 0
        prob_t(:) = 1;
        prob_t(padIdx)   = 0;
        prob_t(startIdx) = 0;
        % Re-apply END block so fallback cannot override grammar constraint
        if bracket_depth > 0 || paren_depth > 0 || any(ring_open)
            prob_t(endIdx) = 0;
        end
        s = sum(prob_t);
        if s > 0; prob_t = prob_t / s; end
    else
        prob_t = prob_t / s;
    end

    nextIdx = randsample(vocabSize, 1, true, prob_t);
    tokens(end+1) = nextIdx; %#ok<AGROW>

    % Update grammar state counters
    if     parenOpenIdx    > 0 && nextIdx == parenOpenIdx
        paren_depth = paren_depth + 1;
    elseif parenCloseIdx   > 0 && nextIdx == parenCloseIdx
        paren_depth = paren_depth - 1;
    elseif bracketOpenIdx  > 0 && nextIdx == bracketOpenIdx
        bracket_depth  = bracket_depth + 1;
        inside_bracket = true;
    elseif bracketCloseIdx > 0 && nextIdx == bracketCloseIdx
        bracket_depth  = bracket_depth - 1;
        inside_bracket = (bracket_depth > 0);
    end
    % Update ring closure state: toggle ONLY outside brackets.
    % Inside brackets, digits are isotope/H-count notation, NOT ring closures.
    % Toggling inside brackets would corrupt ring_open and cause systematic
    % ring mismatch failure in all subsequent ring closures.
    if ~inside_bracket
        for dk = 1:9
            if digitIdx(dk) > 0 && nextIdx == digitIdx(dk)
                ring_open(dk) = ~ring_open(dk);
                break;
            end
        end
    end

    if nextIdx == endIdx; break; end
end

% Decode: skip START token; stop at END/PAD.
% Restore two-char atom placeholders: L -> Cl,  R -> Br.
charList = {};
for i = 2:numel(tokens)
    idx = tokens(i);
    if idx == endIdx || idx == padIdx; break; end
    if idx == startIdx; continue; end           % skip spurious embedded START tokens
    if isKey(idx2char, idx)
        charList{end+1} = idx2char(idx); %#ok<AGROW>
    end
end
smiles = string(strjoin(charList, ""));
smiles = strrep(strrep(smiles, "L", "Cl"), "R", "Br");
end

function props = batchDescriptors(smilesList, featNames)
% COMPUTE DESCRIPTOR MATRIX [N_valid x numel(featNames)] from SMILES list.
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
    catch; end
end
props = props(valid, :);
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:1b68a8eb]
%   data: {"dataType":"text","outputData":{"text":"[23:30:55][INFO]  initPython: Python 3.10 already active (Status: Loaded) -- skipping.\n","truncated":false}}
%---
%[output:0902b7b6]
%   data: {"dataType":"text","outputData":{"text":"[23:30:55][INFO]  Toolbox: Deep Learning=1\n","truncated":false}}
%---
%[output:334b2162]
%   data: {"dataType":"text","outputData":{"text":"[23:30:55][INFO]  R05: Setup complete (HIDDEN=128  EPOCHS=150  TEMP=0.80  AUG=1x)\n","truncated":false}}
%---
%[output:3fd1b0d0]
%   data: {"dataType":"text","outputData":{"text":"[23:30:55][INFO]  Local fda_drugs.csv:          200 rows\n","truncated":false}}
%---
%[output:525fc4be]
%   data: {"dataType":"text","outputData":{"text":"[23:30:55][INFO]  Corpus 201-500: Using only FDA drugs + ChEMBL approved drugs\n","truncated":false}}
%---
%[output:751f8380]
%   data: {"dataType":"text","outputData":{"text":"[23:30:55][INFO]    (Excluding everyday\/forensic -- domain mixing degrades small N models)\n","truncated":false}}
%---
%[output:422f2707]
%   data: {"dataType":"text","outputData":{"text":"[23:30:55][INFO]  Local corpus: 200 SMILES (before ChEMBL fetch)\n","truncated":false}}
%---
%[output:956d762e]
%   data: {"dataType":"text","outputData":{"text":"[23:30:55][INFO]  ChEMBL cache hit: 450 \/ 865 rows loaded\n","truncated":false}}
%---
%[output:55200a7a]
%   data: {"dataType":"text","outputData":{"text":"[23:30:55][INFO]  After deduplication: 649 unique SMILES\n","truncated":false}}
%---
%[output:0617c68c]
%   data: {"dataType":"text","outputData":{"text":"[23:30:55][INFO]  Training corpus: 500 (target=500)\n","truncated":false}}
%---
%[output:6cc88e94]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 50\/500) Validating SMILES\r[##--------]  20% (100\/500) Validating SMILES\r[###-------]  30% (150\/500) Validating SMILES\r[####------]  40% (200\/500) Validating SMILES\r[#####-----]  50% (250\/500) Validating SMILES\r[######----]  60% (300\/500) Validating SMILES\r[#######---]  70% (350\/500) Validating SMILES\r[########--]  80% (400\/500) Validating SMILES\r[#########-]  90% (450\/500) Validating SMILES\r[##########] 100% (500\/500) Validating SMILES\n","truncated":false}}
%---
%[output:1340b2ac]
%   data: {"dataType":"text","outputData":{"text":"[23:31:32][INFO]  Valid SMILES: 500 \/ 500\n","truncated":false}}
%---
%[output:73633f29]
%   data: {"dataType":"text","outputData":{"text":"[23:31:32][INFO]  SMILES length: min=3  median=38  max=260  mean=45.7\n","truncated":false}}
%---
%[output:9a34ceb7]
%   data: {"dataType":"text","outputData":{"text":"[23:31:32][INFO]  Truncated sequences (> 98 chars): 3.8%\n","truncated":false}}
%---
%[output:45463a39]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdwAAAEeCAYAAAAgg6RKAAAQAElEQVR4AeydXWhl2ZXft\/uq60pkJsgpy2qBCXLTnpcgpukXKQyomX5JnKYghjikVI5bD8M8ZNwjYtL4IQwJQQOJNaQJDgMBQbAhMs5DHqYIkzwkEAjGhYdATyUw5KE8IKWkgs48GOvDbpU093e7t3vfXeee74+9z\/k3vbT3Wnvtvdf6nXP20r1XKr10q\/9EQAREQAREQAQaJ\/CS0X8iIAIiIAIiIAKNE1DBTUOsMREQAREQARGoiYAKbk0gtYwIiMBwCbz66qvRJ9+HHEK\/CCq4oV+hcONTZCIgAhMCFKonT55MenH\/Tw7kEncWYUevghv29VF0IiACAROgQFGoAg6xUGjkQk6FJsk5NwEV3NyounHk5vdlXiT4zRurw+6v7+t17GHXaHJtu0eZlrispM1\/dfIWo\/Wj9X2TbL4PuvWjtYIdQadNE3zmSdI83zfJx7cxx7eFprcZY9peaWNJzPBPkizfpHFs7lroknYJqOC2y7vUbnzX6QoPTamFNKkSAbi712HeYr4fc7D5\/kk218cdZw3GbEs\/rzAnSfz57Of7YfP9pH9MADbw+lib\/crYrKW8xh6++Kuxn+uDXsaHOayTNJ8xSTUCKrjV+HUye94Dgb1oQEUerDLr54knKYam9soTTxUfckmKHRtjVdZuci6xEaO\/BzbGfHtFvdHpIcTbJjfyZT8XKjp2a6OPzeq06NjpS9ohoILbDufad9HDUjvSxhfkmrmboM878LAz7vp31Q8ljq7yb2pfrvE8cff0fdwx9eMioIIb1\/VKjZYH03VAdyVpDJv1oY+g0yL0Edun9YVxK+4YNlen79roI9Zu+1andYVxK67d9hmjT2sFPY9Yf1rfHxuCnRahnyQUJ8aRpPFQbW3EDRMrSRwYw05rBd0XO+a2+Fjd76NbsT601mZbbK5Ye9nWMp03n\/FE8X7i2fchxnlryh42ARXcsK9P6eh4KNMeVDvGBm4fHWE+kjTGuBXXB190O5bV4o\/gR4vQTxLWZdwKehU\/dy5r2XVp0d1xbAg2WoT+PGEcYR0r83zbsNsY3DZpX2JGsvyS5mbZWJO1raAnzcFufWjRXT907FYYo2\/bpD5jiD8XHTtCn7muYGOsSyGeLvfX3vUSUMGtl2fQqxV5ePFFshLK45O1RtY4B5+\/Dzp2fy5235ams4Y\/Bx172rw8Y6xjZd56jPtj6Njz7JHHh7V8SZvn+hJLmm+eMdZgTdcXHbtro4+dtgkpunZR\/6SYWSMpzyRf2XIRiNpJBTfqyzc\/ePug87Bbme+tkboIwDppLXs9ksZCsDUdN+v7EkLeNgZ7fdwY7VjV1q7tr+Pu5fd9X+n9IKCCG+l15AHlQU4Ln3FXmJPmr7FuCHCN7LWhRe8mkuZ2JackKboja8DICnrRNeb5s5Yr7DHPFzu+WT74zRPmzxM7p8r6do2iLXsSV9F58s8m0HjBzQ5BHkUJ5Hkg8Cm6bhn\/IvsU8XVj4eH356Jjd\/3K9FmDtdy56NhdW94+85jv+2NjzLeHohMbMfrxYGPMtxfRmc86\/pwkm+\/j68xhPSv+eFmddcvOzTOPeJvaI2lt9sJuY6OPzeq06NjpS9ohoILbDudKu\/BguJLnIcHHnUMfmx8INsYQfyyP7s5nDXQ7jz42K+h2zG2xWx\/X7vZdH3zR3fEqfdZiTSvoda7HullrMp7mxxgx2ZY+gp4kjFlJGsdmx21rY2DMCjY7ntZaf7+1c1jHH8Nmx4u0\/jr+XNa1Pv7YPN2dY+dim+dv7fjgb\/W0Ft+08aQx5rC+K9h8X2x1+LAGa\/nrS6+HQDAFlwttJSk1O0ZbZjxpTve27Ai4+X2ZNws\/dwzdFXfM7Vsfa0O3fbf17VanteL607d2WqvT+sI4Yu1u37VhR6zNbZPsSTZ3ju3jZ8Xa\/JZx3zZPx9cV34+xIjb8rdh5Vk9qi\/hYX1p\/LWxZ4s9xdXeua6fvjtHPsvHs4+MLdua7Yn2sDd32bevb0F2xflktc9wY0OfNSRtLm8M8K3n8yviQA3vMmyt7dQJBFFx7obnYCLqbGjp2K+hFxl1f9UVABOIkwPPPs+8L9q4zCiGGqgz6kENVBk3P77zg8vCkXeikcfyxA4cWnb4VdOxWVxsnAUUtAj4Bnm1ffB\/pIhAqgc4LbhdgTk5OjEQMdA\/oHtA9EO890EXtqLpn5wWX71bdJHhl6tvc8ap9HrD33nvPbG9vS8Sg8D3wld\/4DfMHb7xhaLu9h3T\/iv+w74H79+9PXzhVrQltzu+84LrJNl1s2YuC++jRI3NwcGCOjo5Kyd7eHktVWqPs3qHMGyqD7\/ze75l\/9mu\/Zr72ta\/pHhj4czDUZ8A9g7piwL6c49OHMKIvwRTcNoqte102NzfN1tZWKWEua9GWXSP2eeQ+SAaT++azn\/2sGWz+zjMTMoM2nq+h5w\/jrhjYfTmDYpIgCm7bxTamC6RYAyOwsWHMw4fGfP7zgQWmcERABEIn0HnBzSq2fJ6LjwsSHTs2WnT6VtCxW12tCIiACIiAJaC2KwKdF1wSp0D6gt0KxdMdR7djtOhp4\/hIRCAXgc98xpgkMfpPBERABKoR6LzgUiyTxE\/L9fHH0LPG8alLXnnlFfP1r3\/d0Na1ZmzrkHvvGFBob2+NSRLvAvUyfy\/HLHXoDIaeP\/dHHxmQV1PSecFtKrEm1+Ume+edd5rcIvi1h85g6Plzgw6dwdDz1z0AgWKigluMl7xFQAREQAREoBSB+AtuqbQ1SQQSCPBWMm8rJwzJJAIiIAJVCajgViWo+XEToMC6QjaubvvYkcePjbl3z5hnz9AkIiACIpCbgApublRROiroLAK8qs0jWetoXAREQAQyCKjgZgDSsAiIgAiIgAjUQUAFtw6KWiNOAn7UvH3s29Dn2RmTiIAIiEBOAiq4OUHJTQREQAREQASqEFDBrUJPc\/tBgFewCNnQ+sJnvIwNS5StCIhAzQRUcGsG2vZy\/LnBqtJ2zMHtR0FFCIzWF+wSERABEahIQAW3IsAup1Nof\/sb3zRV\/xB3jH\/IuXXuvOptfVNtGCwBBSYCJQio4JaAFsoUCu6f\/emfmDcffMu8\/e77peSNL79jHj16FEpK3cdBYbVio0HnVS+ttakVAREQgYIEgim4\/LWfebEz5kqSX9Z40py+2F750utmray89npfMNSTB4XVigpsPUy1ytAIKN85BIIouBTLOfEZxty\/BEQfm+uPjt0KujuuvgiIgAiIgAh0TaDzglu1ODKfQuuCRMfu2tQXAREQAREQgS4JTAtulwFQHJEuY9DeIjAlwFvJvI1sxerTQX0RAREQgWoEOi+4WeFTjHm16gq2rHlZ4\/zAkS\/X19cmJrE53t7cmJuSYteIKe9GY\/3oI3Nthfvhk76ZFGHst5MijFw\/fx7VvdIoMzhJdD80dA\/45zS6Pbdia4MvuBRaCqwr2KqC3tnZeeHXafb3983x8bEnL+pc8LOzs0y\/PGtV8SEGOFxcXJrz8\/NScnF5wRLm9PS0UD6hMKjCr8jcnzx5MuXz57\/6q+bPv\/Mdc\/KLXxj4F1mjb75Duwf86zf0\/OHRBoPDw8MXzmrO7+nBFdmX4AtuUzwPDg7M0dHRjOzu7pq1tbVMWV1dNSsrK5l+edaq4rO8vDzFs7i4aJaWlkrJ+M54ukbRfEJhUIVf0twvvvqqSRLft6\/5+3mm6UNnMPT8uTfaYMC57J\/Ve3t703Mrti+DLbibm5tma2trRtbX1w3FK0sobuPxOJdv1lpVx7nhRgsjs7CwUEpGoxFLFM7FMqgaf1DzJ9+0mMlbxknix9nL\/CffuPl5pulDZzD0\/Lk32mDAueyf1Zzf04Mrsi\/RF1zeavbfYkbHHtm1ULgiIAIiIAI9JhB8waVwUkBdweZeE\/S0cddXfRFoj4B2EgEREIFPCQRTcCman4Y122PMldnRj7Ws8Y+99FUEUgjwdvJnPpPioCEREAERKE8gmIJbPgXNFIGaCNhiS+tLTVtomU8JqCcCQyOggju0K6585xPgFe48sbMePzbm3j1jnj2zFrUiIAIikIuACm4uTHISAREQgTYJaK8+ElDB7eNVVU4iIAIiIALBEVDBDe6SKKBOCbif3RIIOq1EBEQgGAKxBqKCG+uVU9z1E6C42s9w7ero2K2uVgREQARKElDBLQlO00RABERABESgCIF2Cm6RiOQrAiIgAiIgAj0koILbw4uqlERABERABMIjoILb\/TVRBKEQsJ\/X2s9saRHsocSoOERABKIloIIb7aVT4I0QoLj6QtFtZDMtKgIiMCQCKrhDutox5qqYRUAERKAnBIIpuPy1nzSmjFtJ8rNjtEnjsomACIiACIhAlwSCKLhZRZJx968BobvQ0NPGXV\/1RaBHBJSKCIhARAQ6L7gUyzRejFNMXR9XnzeO3Z2jvgiIgAiIgAh0SaDzgkvxRLIgUECtZPlqXAQaI7CxYczDh8asrja2hRauiYCWEYHACHRecPPwoNBSlK2g55mX5nNycmJ8ub6+NjGJze\/25sbclBS7Rkx51xmr4SeQs2QCqc49tVZcz5muV7fXyz+n0SePZJT\/R1FwKbR1093Z2THb29szsr+\/b46PjzOFC352dpbpl2etKj7EAJeLi0tzfn5eSi4uL1jCnJ6eFsonFAZV+DH3J0+emDyCryt9yd\/NqWh\/6Awiz7\/Q8z7v3miDweHh4cw5zbnN+T09uCL7EkXBbYLpwcGBOTo6mpHd3V2ztraWKauTtxNXVlYy\/fKsVcVneXl5imZxcdEsLS2VkvGd8XSNovmEwqAKvypzh54\/7IbOYOj5t3UPcC77Z\/Xe3t703Irty2AL7ubmptna2pqR9fV1Q\/HKEorbeDzO5Zu1VtVxbrjRwsgsLCyUktFoxBKFcwmJQVWGZeYPPX+YDZ3B0PNv6x7gXPbPas7v6cHV1JeG1g2+4PJ2ctpntknj+GNviJmWFQEREAEREIHCBIIvuGRE8aSIWkHHbgXdjtGi2zG1IiACIiACIhACgWAKblaRZNzKi+CMsWO0SeOyiUAtBB4\/NubePWOePatlOS0iAiIwHALBFNzhIFemIiACIiACQySggjuAq64URUAEREAEuieggtv9NVAEIRFw\/xEM4kKnlYiACIhARQIquBUBanrsBJz4Ka63t8Yg1kwfu9XVioAIiEBJAiq4JcFpmgiIgAiIgAgUIaCCW4SWfEVgYASUrgiIQH0EVHDrY6mVREAEREAERGAuARXcuWg0MDgC9vNa+5ktLYJ9cDCUcDYBeYhAMQIquMV4ybvvBCiuvvQ9Z+UnAiLQCgEV3FYwa5NoCPCK1gZLH7G6WhEQgdwE5PgiARXcF5nIMlQCFFde3dr86SPYrU2tCIiACJQkoIJbEpymiYAIiIAIiEARAp8W3CKz5CsCIiACIiACIlCIQDAFlz+rlxU5PkiSH3YrSeOyiUAtBDY2jHn40JjV1VqW0yIiIALDIRBEwaVQthGHNwAAEABJREFUVkHOfP4snxX0KuslzJVpCATs57V8ZusK9iHkrxxFQAQaJdB5wc1bHPGjoPo0kuz4Yfd9pYtAJgGKqy+Zk+QgAiIgAtkEOi+4FEckO9R6PU5OTowv19fXJiaxRG5vbsxNSbFrVMo7Mm7KNa77XNdr2NfLP6fR7bkVW9t5wc0DjFerdRflnZ0ds729PSP7+\/vm+Pg4U7jgZ2dnmX551qriQwzwu7i4NOfn56Xk4vKCJczp6WmhfEJhUIVf0tyFl182SeL79jV\/P880fegMhp4\/90YbDA4PD2fOac5tzu\/pwRXZl+ALbhPFlmt0cHBgjo6OZmR3d9esra1lyurqqllZWcn0y7NWFZ\/l5WVSMYuLi2ZpaamUjO+Mp2sUzScUBlX4+XO\/+Oqr5uryMlF83wL5d36f+LHXpQ+dwdDz5z5qgwHnsn9W7+3tTc+t2L4EX3ABStG1YnXaKrK5uWm2trZmZH19fVq8KGBpQnEbj8e5fNPWqWMMBqOFkVlYWCglo9GIJQrnEhKDOjiyBiBo80gf88+Tt+szdAZDz597oQ0GnMv+Wc35zfMamwRfcHkr2RUAo9NKRKB1Ao8fG3PvnjHPnrW+dS83VFIiMCACwRfcrGtB8eXVr+uHjt21qS8CmQT46WR+HSjTUQ4iIAIiUJxA9AWXlCmuFFkr6NglIlCIgC22tL4UWkjOIlALAS3SMwLBFNy8RXKeH3YrPbtGSqctArzCnSdtxaB9REAEeksgmILbW8JKTAREQAREoH4CEa6oghvhRVPIIiACIiAC8RFQwY3vminiJgm4n92yDzqtRAREQAQqEmix4FaMVNNFoGkCFFf7Ga7dCx271dWKgAiIQEkCKrglwWmaCIiACIiACBQhoIJbhFaDvlpaBERABESg3wRUcPt9fZWdCIiACIhAIARUcAO5EAojjUBLY\/bzWvuZLS2CvaUQtI0IiEB\/Cajg9vfaKrMyBCiuvpRZR3NEQAREwCOggusBkTpgAryanScBY1FoIiACcRBQwY3jOinKNgj4r2zR2de29Dc2jHn40JjVVTSJCIiACOQmoIKbG5UcB0mAYsur3kEm34eklYMIhEMgmILLX\/pJw8K4lSQ\/O0abNC6bCIiACIiACHRJIIiCm1UkGbd\/CYgW3YWGjt0KujuuvgiIgAiIwIsEZGmXQOcFN6s4Mk4hnYclaRx\/7PPmyC4CiQR46zhJeFs5cYKMIiACIpCfQOcFl+KI5A+5Hs+TkxPjy\/X1tYlJLInbmxtzU1LsGjHl3VisH31krpMksvuiMT7iENX50Jf7wD+n0e251Wxb\/+qdF9yslPxizCtX35a1RtL4zs6O2d7enpH9\/X1zfHycKVzws7OzTL88a1XxIQZyu7i4NOfn56Xk4vKCJczp6WmhfEJhUIWfP3fh5ZdNltg5fczf5pa3HTqDoefPfdIGg8PDw5lzmnOb83t6cEX2JfiC6\/Ksq9iy5sHBgTk6OpqR3d1ds7a2limrq6tmZWUl0y\/PWlV8lpeXScUsLi6apaWlUjK+M56uUTSfUBhU4efPvbq8nLKgtYLB9mnXPvzQ\/PXf+R2zOnnruSgzf7\/Y9T7eA0WuydDzh1UbDDiX\/bN6b2+PRzM6iabg5ii2heBvbm6ara2tGVlfX58WLwpYmlDcxuNxLt+0deoYI+nRwsgsLCyUktFoxBKFcwmJQR0cp2tMvmkxk89rp\/3JNzG0U31ip28FZkuT8VDuARtX220v74HJdc3Lcej5w6kNBpzL\/lnN+T09uCL7UlvBpSAm5T7PnuQ7z8YadbyNPG992UVABERABESgaQK1FdymAs0qthRifNz90bG7tsH3BUAEREAERKBTApULLsUNIQtaX+oofP6a6OxnhT2wWUG3Y2pFIDeBydvJZvLZ7AuCPfcichQBERCBZAKVCy7FDWF5Wl+w5xHmJflhTxLf1\/Xxx6SLQAaBT4cprr58OqqeCIiACJQmULnglt5ZE0VABERABERgQARqK7i8whwQN6UqAsMhoExFQARqIVBbwa0lGi0iAl0TcD\/DJRZ0WokIiIAIVCRQa8G1P7TktxVj1HQRaIcAxdV+fmt3RMdudbUiMEtAmgjkJlBbwaXI8rZykuSORo4iIAIiIAIi0FMCtRXcnvJRWiIgAiIgAmUJaN4MARXcGRxSREAEREAERKAZArUVXN5K5m3lZsLUqiLQAgH7ea39zJYWwd7C9tpCBESg3wS8gls+WVtsaX0pv6pmikDLBCiuvrghbGwY8\/ChMaurrlV9ERABEcgkUFvB5RXuPMmMQg4iIAIiIAIi0HMCtRXcnnOaptfnLycnJ6aInJ2dGYQ50XPhbWPETQTdimtXXwREQARKElDBLQmub9N2dnbM9vZ2bnnrrbfMgwcPDO39+\/enxTpKJhRV+xayTcC1MYZux9SKgAiIQEkCtRVc\/3NbV88TG\/5pfoxbSfKzY7RJ47KlE3jzwbfM2+++n1v+9j\/6A\/Obv\/X75vW\/9XXz6NEjY0z6+hoVAREQgaETqK3gJn1+C1zstGmSVSQZZx0r6O566HaMFt0dVz+bwCtfet2sFZBXXvt18\/kvbkznZK8ekQevZnlVG1HIClUERCAOArUV3KR08xS\/rOLIOOu466Njx0aLTt8KOnarqxWBrglofxEQARFotODmwUtxRPL4ykcEOifw+LEx9+4Z8+xZ56EoABEQgbgIdF5wu8LFT9f6cn19bWISy+725sbclJVP3j4tusbtZB572hhi4jYT60cfGcPbyJ\/I9US349bu2sgbuX7+PKp7xeYUZhvXcyeG7V4v\/5xGN5H+V1vB5S3cJAn11WvST+Xu7++b4+PjTOGC8ysxeXyb9CEG7ruLi0tzfn5eSq6uLlnClFnj8vLSXF5cTOefnp5mcmuSRZW1f\/LkibHirpNke\/r0qSHv00kLf9d\/aP1QnoOuuA89f7i3weDw8PCF357g\/J4ePJF9qa3gUliTJFQeBwcH5ujoaEZ2d3fN2tpapqyurpqVlZVMvzxrVfFZXl6e4l1cXDRLS0ulZHxnXHoN9r0zvjOdHwKPKizzzr17964Zj8fmc5PrP5Sc57EJ5TmYF1\/T9rbybzqPKuu3wYBz2T+r9\/b2pudObF9qK7ixJb65uWm2trZmZH193VBEsoTixqGb5dfGONxHCyOzsLBQSkajEUuYUcE1RpN5H8vCdH4buYayB3kvTb7JCeUe6IpLSM9BFwyGnj\/M22DAueyf1Zzf04Mnsi\/BF1xeNfNWtcsVHTs2WnT6VtCxW12tCIiACIhAHwnElVPtBZdiZ6UuFBRPuyYturs2OnYr6O64+iIgAiIgAiLQNYFaCy4Fj2JnBT1vgsxJ82XcSpKfHaNNGpdNBERABERABLokUFvBpbj6xQ4du5OguiIQLgF+NSgpunn2JF\/ZREAERGAOgdoK7pz1ZRaBeAjc3prp794a5z+KLXbHpK4IiIAIlCGggluGWlNztG73BCiuFFkr6N1HpQhEQAR6QKC2gpv09jFvJ2PvASelIAIiIAIiIAKVCNRWcImC4kqRtYKOXSICNRBoZwle2fKq1gp6OztrFxEQgZ4TqLXgwooiawVdIgLREKC4UmjdgNGxuzb1RUAERKAEgdoLbokYNEUEwiBAcU2KxLVvbBjz8KExq6tJnt3ZtLMIiEDwBCoXXPv2sZtpks0dV18EgiXAq9kkCTZgBSYCIhALgUoFl8Ka9PaxtTEeCwjFKQLTXwni1SwCDr\/FJomRgGIWgSAIlC64FFMKa1oWjOOX5qMxEQiaAEWXV7xBB6ngREAEYiBQuuDGkJxiFAEREAERyCCg4dYIqOC2hlobiYAIiIAIDJmACu6Qr75ynyXgvn1s+7ydTH\/WU5oIiMAwCNSaZemCm+fzWT6\/xa9qxKzjStJ6WeNJc2QTgRcIuMWVPuI6PX5szL17xjx75lrVFwEREIFMAqULLitTTCl09H3BzrhvL6rbdVjLCjZ3HXQ7RovujqsvAnMJ8Ao2j8xdQAMiIAIikI9ApYLLFrbAUeRcwc5408Ke7l7sh46dvkQEUgm4r2Dpz5PURTQoAiIgAtkEKhdctqDA+YI9ZDk5OTG+XF9fm5jE8r29uTE3ZYUCM1mo6Bq3k3nsOZk6\/T8mbi\/E+tFH5noi09\/D\/eTV7gs+n9wb5I1cP38e1b0yLx\/Z43rmh3i9\/HMafXroRPilloLbZN4Ucl6tuoKt6p47Oztme3t7Rvb3983x8XGmcMHPzs4y\/fKsVcWHGOBwcXFpzs\/Pc8iLPldXlyxhyqxxeXlpLi8upvNPT08751GFJXN\/8uSJsbLw8svGCmNWnj59asj7dNLC39qH2IbyHHTFfuj5w70NBoeHhzPnNOc25\/f04InsS\/AFl0JLgXUFW1XOBwcH5ujoaEZ2d3fN2tpapqyurpqVlZVMvzxrVfFZXl6eYlhcXDRLS0ulZHxnXHoN9r0zvjOdHwKPKiz9uVeTbyasfPHVV395re\/evWvG47H53OT69y1nn0GWHspzkBVnU+NDzx+ubTDgXPbP6r29vem5E9uX4AtuU0A3NzfN1tbWjKyvrxuKSJZQ3Dh0s\/zaGIfPaGFkFhYWSsloNGIJMyq4xmgy72NZmM5vI9em95hZf\/INzOInYiZvn7tj5L00+SYnlHvAja3NfkjPQZt5272Gnj8c2mDAueyf1Zzf04Mnsi+DLbiRXSeF2waBTz6\/nX6WOymyFNqptLG39hABEeg9gegLLm81+28xo2Pv\/dVTgvUQsIVWRTaFp4ZEQASqEgi+4FI4KaCuYHMTR08bd33VF4EZAhRba6A\/T6yPWhEQAREoSSD4gkteFFRXsPmSNe77SxeBKQH3VW1af+qsLyKQTEBWEchDIIqCmycR+YiACIiACIhAyARUcEO+OopNBERABKInoAQsARVcS0KtCIiACIiACDRIQAW3QbhauocENjaMefjQmNXVHianlERABJokkFRwm9xPa4uACIiACIjAIAmo4A7ysitpERABERCBtgmo4BYlLn8REAEREAERKEFABbcENE0RAREQAREQgaIEVHCLEpN\/GgGNiYAIiIAIzCGggjsHjMwiIAIiIAIiUCcBFdw6aQ54Lf4QdRWJBt3jx8bcu2fMs2fFQ9YMERCBQRNQwR305a8v+Z2dHbO9vV1a7t+\/byjY9UWklURABEQgLALRFFz3rwElIcwaT5ojW30E3nzwLfP2u++Xkje+\/I559OhRfcFopRgJKGYR6D2BKAouxdT9a0Do7pVBTxt3fdVvhsArX3rdrJWV115vJiitKgIiIAIBEQi+4Npi6jKjuFp93jh266NWBERABKImoOB7QSD4gmspU0CtWFuVls8Lfbm+vjZtir9\/Ud3mf3tzY27Kyu3tdJmia9xO5rEn81iAFr2UTNZijTbZl92LvJHr589bvVfKxqt57T7T4l0\/76RzkfMiRomi4FJoeVVrBb0q7KQf8t93fbcAABAASURBVNnf3zfHx8eZwg1wdnaW6Ze21o9\/\/GPz7rvvlv4hI35AiRzgcHFxac7Pz0vJ1dUlS5gya1xeXpqrq6vS823MNobT09NKTNN41zX29OlTQ96nk7bqPVBXTF2tU8dz0FXsdew79Pxh2AaDw8PDF85Je\/ZNDx9jommiKLgU2rqJHhwcmKOjoxnZ3d01a2trmbK6umpWVlYy\/dLW+vnPf24++OADU\/WHjeCyuLholpaWSsn4zpglTJk1mDNeLD\/fxmxjqMo0jXddY3fv3jXj8dh8bnL9Y4i3rryT1qnjOUhaNxbb0PPnOrXBgHPZP6v39vam51ZsX6IouE1A3dzcNFtbWzOyvr5uKCJZQqHg0M3yyxonrzp+2Gi0MDILCwulZDQaEYYZFVxjNJmHvPRSufluvKPJWgSRxSuUceJdmnyTU8c9EEpOZeKo6zkos3cIc4aeP9egDQacy\/5ZzfnNmRGbdFJwY4OkeEVABERABESgKoHgCy5vJ6d9Zps0jj\/2qnA0XwREQAREQATqIhB8wSVRiidF1Ao6divodowW3Y7F1ypiERABERCBPhKIouACniJqBd0XO0brj0kXAREQAREQga4JRFNwuwal\/cMgoChEQAREIFYCKrixXjnF3Q2BjQ1jHj40ZnW1m\/21qwiIQLQEVHCjvXQKXAR8AtJFQARCJqCCG\/LVUWwiIAIiIAK9IaCC25tLqUREQATSCGhMBLomoILb9RXQ\/iIgAiIgAoMgoII7iMusJEVABEQgjYDG2iCggtsGZe0hAiIgAiIweAIquIO\/BQSgEIHHj425d8+YZ88KTZOzCIhAvATqilwFty6SWkcEREAEREAEUgio4KbA0ZAIiIAIiIAI1EWgnwW3LjpaRwREQAREQARqIhBVweUvASFJuWO3kjQumwiIgAiIgAh0SSCqgjsPFIWWvxJkBX2er+xGCERABERABDogEE3BpYhSUH1GSXb8sPu+0kVABERABESgKwLRFNyuAGnfgRFQuiIgAiLQEIEoCi6vVnnVWieDk5MT48v19bVpS2wutzc35qas3N5Ol+lijdvJ3sTN3gRBi15KJmuxRlvsq+xD3sj18+et3StV4tXc9p5psW6GtX9Oo3NexCjBF9wmii0Xamdnx2xvb8\/I\/v6+OT4+zhQu+NnZWaZf2lrMJ46Li0tzfn5eSq6uLlnCdLXG5eWlubq6qi2G09PTSkzTeNc19vTpU0Pep5OWa1jXupGsM3N96ngOYszbxjz0\/OHQBoPDw8OZc5pzm\/N7evBE9iX4ggtPiq4Vq9NWkYODA3N0dDQju7u7Zm1tLVNWV1fNyspKpl\/aWsvLy9PwFxcXzdLSUikZ3xl3ugaxjxfri6Eq0zTedY3dvXvXjMdj87nJ9Y8h3rryTlqnjucgad1YbEPPn+vUBgPOZf+s3tvbm559sX0JvuDyVrIrAEanrSKbm5tma2trRtbX1w1FJEsokBy6WX5Z48Q\/WhiZhYWFUjIajVjCjDpYYzTZG3nppfpiyOIVyvhokvvS5BulOu6BUHIqE0ddz0GZvUOY80L+k3sihLjajKENBpzL\/lnN+T09\/CL7EnzBzeJJ8eXVr+uHjt21qS8CIiACIiACXRKIvuACj+JKkbWCjl0iAiIgAiIQBAEFMSEQXcGdV0yxW5nkpf9FoBkCGxvGPHxozORz\/GY20KoiIAJ9JRBdwe3rhVBeIiACIiAC\/SYwt+D2O21lJwIiIAIiIALtElDBbZe3dkshwO\/0VZGUpTUkAiIgAp0TUMEtdQk0qQkC\/DI7v9ReVu7fvz\/918OaiE1rioAIiEBVAiq4VQlqfm0E3nzwLfP2u++Xkje+\/I559OhRbbFoIREQARGom4AKbt1EtZ4pi+CVL71u1srKa6+X3VbzREAERKAVAiq4rWDWJr0h8PixMffuGfPsWW9SUiIiIALtEFDBbYezdhGBTwioEQERGCoBFdyhXnnlLQIiIAIi0CoBFdxWcWszERCBNAIaE4E+E1DB7fPVVW4iIAIiIALBEFDBDeZSKBAREAERSCOgsdgJRFNw7V8Cok2Cjt1K0rhsIiACIiACItAlgSgKLoXU\/iUgWnQXGjp2K+juuPoiIAIiIAL9JhBDdsEXXIonhXQezKRx\/LHPmyO7CIiACIiACLRNIPiC2zYQ7ScCIiACIiACTRDoruDmzIZXq64rr1x9mzuet5\/0V2mur69NW2LjvL25MTdl5fZ2ukwXa9xO9iZu9iYIWvRSMlmr6hrMR5q+fuSNXD9\/3tq90nROWr+9516si7NOOqt51mOU4AuuC7WuYsuaSX+ZZn9\/3xwfH2cKN8DZ2VmmX9pazCeOi4tLc35+Xkquri5ZwnS1xuXlpbm6uuo0Bsvu4vJiGsfp6Wml65J2zRh7+vSpIe\/TScs1xDZUqeM5iJnd0PPn2rXB4PDw0Ph\/QYzze\/rAR\/YlmoJbZ7HlGh0cHJijo6MZ2d3dNWtra5myurpqVlZWMv3S1lpeXiYMs7i4aJaWlnzJpY\/vjDtdg9jHi93GYNlZFlWvS9o1Y+zu3btmPB6bz02uf9N7sV\/IUsdzEHJ+WbENPX\/4tMGAc9k\/q\/f29qZnX2xfoii4dRdbLtLm5qbZ2tqakfX1dUMRyRIOeQ7dLL+sceIYLYzMwsJCKRmNRixhRh2sMZrsjbz0UncxuNxGk3iA8eGHH5oqknXNGGevpck3SnXcA6wXq9T1HCj\/xVznToic2rgHOJf9s5rzm+c9Ngm+4GYVWz7PxccFj47dtanfIwIpqfBWk\/\/2UxH9\/v37hrfJ5m6xsWHMw4fGTN7lmOujAREQARFIIBB8wSVmCqgv2K1QXN1xdDumdlgE3nzwLfP2u++Xkje+\/I559OjRsIApWxEQgdYIBF9wKZ5J4hNyffwx6cMhoD9ib4ZzsZWpCERGIPiCGxlPhSsCIiACIiACiQRUcBOxyCgCItBLAkpKBDokoILbIXxtLQIiIAIiMBwCKrjDudbKVAREQATSCGisYQIquA0D1vI9I\/D4sTH37hnz7FnPElM6IiACTRNQwW2asNYXAREQARGIn0ANGajgloTIv6PLP5BQVkpuq2ktEEi7pvzb0cjZ5BVu0j3QQnjaQgREIFICKrglLhwH7be\/\/W3z1ltvvfCPauf9V434F5FKbK0pLRDg2sy7jl\/5ylfMD3\/4Q\/PVr37VPHjw4IV74O\/83b+f\/i9VtRC\/thABEQiTQI8LbnPAKbgffPCBqfqvGjUXoVauQiDtuv7Nv\/e75vPrf8O8+Q\/\/qfnN3\/p98+Vv\/Otf\/qtW\/EtVf\/anf1Jla80VARHoMQEV3AoXV\/+qUQV4AU9Nva6vvW4W\/8pfNWtfet18\/osbZu21Xzf0pzIZCzgthSYCItAxgZc63l\/bd0RA24qACIiACLRLQAW3Xd7aTQREQAREYKAEVHAHeuGrpv2zv3hm\/vd\/+7752V+cVV0qwPnZIZH3kPOHED\/L8N3vftfQog9NyHvI+XO9xQAK+aU3Bdf983z505dnWQIUnP\/z379fdnr08342+YZjXv5pv1aUZywUOFmxcth+73vfm\/5UdpZv2fFQWCTFYfNPGhuKTQyKXeleFFyKrfvn+dCLYZC3CNRHIO3Xiub9upFr54\/g\/+hHP6pUyKpmQ4H87W98M\/HX3mys5Mk+tNZWd8uvWXXNghwlIlAHgegLLsWVYuvCQMfu2tQXgbYIpP1a0dvvvv\/LXyNK6vOrRY8ePTJVixhFm6JZNmfm8itOabkQK+vP88kaT8rftTGfGOpgUbVow8MXckd8e1M6e4Ugbn42HteWp2\/nDa2NvuCWvWB5bop5PnbPn\/3\/M\/PTD5+Wkpvb2+kyP\/3wtNT8n0727XKNW\/Nx\/LeR5wFHK3l4fvSLS\/PRR1cmKW87\/1f+2ivmV+6Wk7VPfrWIYkMhKyPMpWjPu3\/z2Kc35+RLWi421nk+WeNZjOx88inDgTnMhUXVop30yp01J4i8b462U98VSFonr41vopr4xiHP\/eD6vPfee7\/MsSyDqrnAPUYZXMH9whe+YDY3Nys9JPYm+8\/f+cfmP\/6LB6Xkj\/\/tN6f3C22Ma\/zXP3xvGv\/\/\/P7BtI01D5c9OZAMrWt3+\/\/+3\/yu+Xf\/7\/+a\/3T0L3E1\/+UP\/8kvrz\/zMHJf\/OCf3zdlhLms8b\/++Lvmf\/yHf1VKmMsa3Kd5D3Pfj7msQTzz8mAszSdrfN661m7nk09VFsQZuzT1jYN\/7bN04qjKkjW4x7L2mjfOXM5xzvOqsbQ5f5AF9+DgwBwdHUnEoPA98J0f\/MD8gz\/6o8LzdL\/pedM9cGTqZMA53maxrGOvwRVcoPFd0dbWlpGIge4B3QO6B+K8BzjHOc9jkkEW3JgukGIVAREQARHoB4H0ghtBjkk\/kcxPKGOPIHyFKAIiIAIiMBAC0RdcrhPFlSJrBR27RAREQAREQARCIdCLggtMiqwV9BZEW4iACIiACIhAbgK9Kbi5M5ajCIiACIiACHRAQAW3A+iD2FJJioAIiIAIzBBQwZ3BIUUEREAEREAEmiGggluCq\/3hLNoS06OZQn5J4ifg+vhjserklBY741aS\/OwYbcJ48KZ5cWNPEj8h18cfi0XPyqHqeOgc5uXn2t2+n0\/amO87FF0Ft+CV5iayP5xFi15wiajcydEXNwHyd8fR3fEY+1k5MJ6Wc9Z46EyIPy1GN3fbd\/2Zb+206O54DH1iJnYr6G7c6HaMFr3IuOsbYp98yMsKuhuntbutO46\/O4bujg+1r4Jb4Mpz03ATuVPQsbu2ofTJm\/zdfNGxu7aY+lmxM06Obk7o2LHRotO3go7d6iG3VeNkPvm6OaJjd22p\/Y4HiZWY54WRNI4\/dubQotO3go7d6iG3xEm8ZWNMms962Muu2Zd5Krh9uZIN5cFDYqWhLYJaloMBCSqoFoMhdyRtS3s\/0Kb5aay\/BLj2VvqbZf2ZqeDWz7RXK3L4WuEB61VySqYUAXs\/0PbxniAvFww5+jZ3vOV+49v5uSblj48VxhsPqicbqOD25EI2kQYPlLsuuh4ul8jw+twDbtbofb4nyI0c3ZyH1E\/K3+eBjt+QuJTNVQW3LDnNEwER6DUBigjFpNdJpiQXZf4p+YQwpIIbwlUINAYeuEBDU1gdERjKPUGeKrZPEu8y2CQOyJhJQAU3E9GnDjyA\/s2Gjv1Tr\/72\/FzJG5ubMTp219anPrmRo5sTOnZstOj0raBjt3qfWj838sTm5oiO3bWF3M+Kl1zwcXNAx46NFp2+FXTsVg+5LRqr70+e2Nwc0bG7tiH2Oy648SHnpuHmsYIeXxb5IiY3myctuj8TG2NW0H2fvunkaPOlRXdzRMduBd0dj7lPLjYvWnQ\/H2yMWUH3fULXbexu68ZMTu4YepFx1zfEvpub7ds4ydXaaNHtmG2xMWYF3Y4NuVXBLXH1uXmslJge1RSbJ+28wBmzMs8nNjv5pMXMuJUkPztGmzQeui254WBQAAABtElEQVQtbsaszMvDjtPO8wnVTsxJ4sfr+vhj6Fnj+IQobtxu3411nr2oj+s\/hL4KbsBXWaGJgAiIgAj0h4AKbn+upTIRAREQAREImIAKbsAXR6GlEdCYCIiACMRFQAU3ruulaEUgWgL8AE0MwccSZwwsFeMsARXcWR7SRKAXBEJLgiLGD9qEFldSPMRJvEljsolAFQIquFXoaa4IiEAmAYoXRSzTMSAH4iXugEJSKD0goILbg4vY1xQ48Fzx82TMtyXp1o\/WivVDt\/15LT7zJGmO75vk49uY49tC09uOsbn9QiOreIZCQAV3KFc6sjw5bHmV4Qo2P40km+vjjrMWY7aln1eYkyT+fPbz\/bD5fkPRyR0eMeZL3MQfY+yKOUwCKrhhXpdBR8Uhx2HnQ8DGmG8PRSc2YvTjwcaYbw9Zjy3ekFnGFpvibY6ACm5zbLVyAwQoXu6y6POKA3bGXf+u+qHE0VX+2lcERMAYFVzdBcERoDhRLJHggksJqI24YWIlKRTGsNNaQffFjrktPlb3++hWrA+ttZVpme\/KvDXSfNwx+v4aro0+gg+tK9gkIpBNoJrHXwIAAP\/\/EGDvgQAAAAZJREFUAwBs2c0sCUboigAAAABJRU5ErkJggg==","height":229,"width":381}}
%---
%[output:233a680b]
%   data: {"dataType":"text","outputData":{"text":"[23:31:32][INFO]  AUGMENT_FACTOR=1: Augmentation disabled (corpus = 500 unique SMILES)\n","truncated":false}}
%---
%[output:6ccf3410]
%   data: {"dataType":"text","outputData":{"text":"[23:31:32][WARN]  Unexpected characters included in vocabulary: [K] -- Generation quality may degrade.\n","truncated":false}}
%---
%[output:7c8775f2]
%   data: {"dataType":"text","outputData":{"text":"[23:31:32][INFO]  Vocabulary size: 37 characters\n","truncated":false}}
%---
%[output:329ce97e]
%   data: {"dataType":"text","outputData":{"text":"[23:31:32][INFO]  Special tokens: PAD=1  START=2  END=3\n","truncated":false}}
%---
%[output:18cc59a8]
%   data: {"dataType":"text","outputData":{"text":"[23:31:33][INFO]  Markov model construction complete: 37x37 transition table (from 500 SMILES, 23229 bigrams)\n","truncated":false}}
%---
%[output:2f235c0c]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdwAAAEeCAYAAAAgg6RKAAAQAElEQVR4AeydDZglR13u6+wOSUiWMCGsYciOhFzgBi5gQHT3guxzEVFEUWB9UBZQBOFRETeioCJqRBQlggYEPxBElNUL7iPII0HkSxEERQhcRO4DhK\/lDmEMO0lmySSZsLd+vVO7tb39XdXndPV59znvVtW\/\/vWvqrf79NvV3XN62zH9EwNiQAyIATEgBnpnYJvRPzEgBsSAGBADYqB3BiS4VRSrTgyIATEgBsRAJAYkuJGIVBgxMDYGLr744rFNSfOJzID2kXaESnDb8SXvkwwoN2IGOJBec801I56hphaDAfYR9pUYseYhhgR3Hray5igGWjDAAZQDaYsmcp1jBthX2GfmmILGUx+c4LLhilA0I\/yK7EOzMc46DGHMjLHJOJr6NYlV50NfRahr16ae+G38876F7T2nfH1d2Wt6Sjbf7pTKHgv0W4ceuz8lNOM4xZArUO+Qq6otunYuLWpAXZFdNjHQhIHBCS6D5owpj5R39PxciuaITShmIM8f5ZT3h+JZDsMKr\/Drj4ayD+r8MnlsfYOxVfVBPWNxoFzl79fh69q5FJvvo3w5A3Amvsr5cTWDFFw3OD8t2qDYfB\/lwxgo4rPoS1TkF9Zz+9aMoWhs7SP134KxVvVSVF80tyK\/qrgDrms9tCI+\/CDUd+WnrC3xqPP7SSGf4phT4DXGGJMR3BiTVQwxIAbSZADxA21G38QfcWri16Zf+YqBMgaSEdyiLwa2\/MSw5eF8sPt5Vyb14XxcSh15UkAekAfkY8GPRx4Qm9QHNh\/UUSZ1oOzD2V3q15HHTgrIg3zelUl94Ovg212eOvKkDpRjgZguFnmQL\/s2V+en1DuU2an36\/w8dQ55u1\/O52njbOQBZVJAHvh5ygCbA+U8XJ1L8\/V9lF1fpEXxnZ3Uocivi83FI+3SPt+mTIyJ75Bv4+wuLap3Nt\/H5S+++GJD3vn4KXYfRXXYnA95B2cjdTY\/9e3kAfWkPrAJ3RgYpOD6G9fly3Z8f9r44udAHXlSB3wAduDnKQNszt+l2KgD5AF5QN75xUiJB4gN\/DxlgC3fFzbqHCg7H\/LO7lJsrj6fOh\/sfp5yHsRxPqSU8z6UsVPvQBl7F9CWOH5bbAA7oM4vY6OMPQ\/s1DtQxofU2VyKjTof2Fw9KWW\/vmmetgB\/UkC+CPRBvQNl34+yq3MpNt8ndp74ri9SykV9YKfegXKRXxsbMVw8Uspt2jf1JS7xHSi7tuSd3aXYXL1LsQHfx+VdSr3zJ6Xs6lyKjTrgbPk8ZfxcPSll7HlgB\/gAP08ZYMu3U7kZA4MUXDZqHrE2sovr6KHs8lVp3i9frmrbto7YwLXz885WlDb1c23b+rt2fsp2ycehjN33I4+dtC2IlUdRLGzAxaeNX8ZOGTt5H9j9ssuX2V29S5v6Of\/QlDnk+6SMvSo2PlX1IXX0nY9PGXs+Lva8LbTcR8yiMVX1U1Xnx8IPOJufdzY\/hcMinyKb3458UVvaYafeB3bg2\/L5uvq8f+Ry0uEGKbhFjLKRi3aQIt+2NuL6aNt+Gv7++Mi37dPxR1uHtjFm5c\/Y85jmWBxfLp1m31V9ufH4qe8PZ34deb9e+XoGunBGGx\/1vRz38NuQP26N8z\/x8mgSWftQE5aa+yQjuE2mlN85KNe1YyfEz0ddm2nXxxqjP0fyxJ32XFLrD47gysdQ5uCPyc\/74\/Pt5JmPX5\/PN\/HJt0m5HHu+8EtMH0346dquSWx8\/PH4eerq4PuTZ6z5Ntioy9tVPpWB3gX31O66l5psUOfDhgddeiNGl3bTbNNljF3aNJkTPOdjU8bepH2fPoyBsfh9UMbu28hjJ61DmV+ZvS5e13rmUNSnb\/PzXftp065oTIwBez4O9rwtpEwfXWOWtSUedSHjIkaX9vl2jCNvI26RDbuPkLZN4vt9KV\/NwCAFl42cBztN9VSO1+bbHbeW\/09cvw3lcu\/Z1DCm0DHmYxAPW92M8MEXlPn6PvhRLvOdtp2xMCYHykVjwO58SCnjR0rZgTL2PLA7H1LKeZ82ZdoTB5S1833wA9icP3lsPrC5+rIUH9qU1VfZXVvaA8pF\/tipd6Bc5NfGRgwXj5Ry0\/b40sYHtqbtnR9tusRo0i7vQz\/YXN8uxUYdKLJhB\/i5+rIUH3x9YPP9qcvb\/HrlTzIwOMFlwxXh5JBP5vBzJbfRsfnA7nywu7yfYnfATp7UoW3ZtStLT8Y73aOsDrsDrciTOuTL2PM2yj7w8UGdX3Z57MAvu7xLqXdwNj+lzi+TL7Jh99HEB\/8qP+oc8M2DOmykDpQdnI0Um0vJA1cmdcDugM3lSevK+AD8AHng5ykDbD6w+fDryPt1VXl8\/e+O70udX87nqXfI1\/ll50Pq26vydb7UO1TFKapz7Vxa5pO34+\/bKDtgJ0\/qkC\/7duoANpeSd8Dmw9nzqfPx7c7mUr+OPHbSPLD78OvZR6jzbcqXMzA4wS0fanUNG52Nnwf26paqFQNioIgBfXeKWJHNZ0D7iM9GfX40gstU2fh5YBfSZECjHicDfEfHOTPNSgxUMzAqwc1P9fDhw0YQB9oHtA9oHyjfB\/LHzaJyLP6KYs+TbbSCyw7ynOc8x+zdu1cQByPcB7Rf67sdZx94whOekC1MyoQv5rG0qi\/\/dmDZWJr4lLUdgn3UgvvBD37Q\/Pjqqnn+ykop9h05km2HOj8X40wbqw4LWzHPsH3X+R6z8WqxFe+YjVfra+OtNsANWzGP2JhN\/K+1MeuwthXzOhuzzvdRNl4dHrAVb6+NV+dL\/f+wMeuwayvmPWzMOl\/q6+bRZs51sVx925hN9oljW\/Nuug9db7mswte24t1oeazyc3X\/ZePV4catmGs2Zp0v9TfYmHW4aSvmURuzzpf6uu9C2+\/NdjvGOmzbGuM2O8Y6X+rdflKVNt2H8OM4mR0ES\/5DcPG54oorzMGDBzvjwIEDhjhF3SCk3GpwoJz3w+bqSSnnfYZeHq3gOuLvvbFh7lMB6vElrfJzddtsrDpstz7ErPOj\/pj1rYOxPsSr83P1N1v\/JiDmRkNf\/OpAny5mne+S7bcJiNfED59zbcw63NH6ELPOz9XXzaPNnOtiufq2Md12r0rb7kO3Wp7qAI+bDfyIc4v1q4ObN2mdL\/X0XQZnp2\/GSepsVSl914F4bKs6P+ondt5NQMwmfvjQdx3om5hN\/fCtw+7dDzJ79nQH7Yv6QDgRUL+OMnbfls838cm3mXU5CcGFeB+zJk39iwExIAbmj4Gb7ZQ3AnCLbTvfn8ELLkLLmYwPbPO92TR7MSAGxEBXBrq227QNu+Hw4S\/YtrdZxPukqAODFlwIRWjzmwgbdXl7l\/LOzU3DfVzSLu2L2kxsTO7jkhbVt7Uds\/Gye3A2bdu2zH\/TxuJ+1G02LfNpaycm94RI27Yt8r+DHRv3cXfYtKi+i+1MG4v7uKRd2ufbMNeYcyZ+HzFj70NftzxyHzfm\/kMs7uOSwkMMME7u45LGiMe2if29MZZL7uPGOl4wT8YZe780BsHsJriHDr3Z7N\/\/DIbWGe64z7EfUO4cbEYNBy240+AEod23tha1K744t1tbixfTfiGPra3Fi2cjcVC7YW3N5uJ9+JJfv7YWLSBC+8C1tWjxCITQLq+tkY2C2HNmUH3E5KAecx9i\/\/na2hrDjQZi3ri2Fi0egRDam9bWyEYBY4z9veF4sW1tLcr4XJBe9iHTTWyNbbdv3yPNgQNPc8PrnCKyDp2D1DTsszpZwYX0JsSsLiwYhyb+8hEDYkAMjJWBBXs89NFunt1XuLt27TS7d9+\/XXc5b1a1OVNyxWQFtynTL1xaMgeWlzMcWlxs2kx+YkAMiIHRMXDOjh3mQns8dLjAHh+bT7L7CpdVrskuSZ\/eG4unvJhSxn6690lLE5+T3sPIJSu4kJ1RWPOf\/\/e1e9fXa7xVLQbEgBgYLwNH7THQ\/\/td7vM2n233Fa6xl5VNieDSP+LKMd2BMnaAjRQbeQfK2FNCsoLblGT\/72u5X9u0nfzEgBgQA2NjgHu7\/t\/mur\/XbTbPfla4rm8E1MHZSLGRAvIOlFPDoAUXYjmbyZOKjbq8XeXTGJBBDIgBMRCJgf5WuJEGOPgwgxZc2ENYEVgf2KgTxIAYEANiYFoM9LvCndYsZtnP4AUXchBYH9gEMRDMgAKIATHQggH90lQLsgpdkxDcwpHLKAbEgBgQA1NkQJeUQ8keveDyg2KfsSzFwidtrJj4oo0XG2fbmLFhQ0b9nGujxQZ\/9BUbO+04Y2K7jRcbD7ExY+MiY8xFxphYOMsYExvn2JixgaTEBAfY2LDTntFHl5RDiWdfCI2h9mJADIgBMTB6BjgVCRFd2o+epMoJSnAr6VGlGBADYqCAgbk0hYgtbSW4gxdcnk4u27er6srayC4GxIAYEANdGEAwEc6uoH2XfsfTZvCCC9USVlgQxIAYEAOzZKCx0NpBFvlKcJMQXP4kSKJr92F9xIAYEAMzYwDBLBLSpjbaz2zwg+g4CcENYerIwoJZ20JIHLUVA2JADKTOgP+mIPLt5tNUWMv8JLjbIDwFdF3lvmZpybxkeTnDuxYXU5iqxigGxIAY6IWBsLcFIZhlYtrETvteppVM0GQEF0a7iO7jVlfNU1dWMjxwXW8LgkdBDIiB+WQg7G1B+qWp0L0mKcHtMtm7b2wYh8VNzsLaRpG\/GBADYmAcDAz5bUHjYLh6FskJbpdVbjUFqhUDYkAMiIF6BrgkzKKlK2hf38uYPZITXDaGRBcWhgGNQgyIgXlhoKvQunbVgstfojiUMerqXVrmN1R7koILmYguqSAGxIAYEAPTYADBdOLZJaV98TgRUI7pDpTznthcvUux5f2GXB684EJsGYFVdWVtZBcD02NAPYmBMTHQRWT9NsWCi2jmj+WUsY+JPeYyeMFlkIIYEANiQAzMmgEE0xfQ5vnDh9fs4Glvkzn+jF5wP2A37nsiwoYa\/OdGO8LYsCGjfi6z0WLjFZ81JjZW7Thj4nY2Xmy80caMjZivoCRW0f4Yaou5XVwsS2XUDzITG1EH2CoYgtkNhw79X7N\/\/1Wtess7u1UvK18HbHm\/IZdHL7hDJl9jEwNiQAykwwBi23xVa8xJ33377m4OHLhf0FQRWQTWB7agoFNunITgQmoeU+ZJ3YkBMSAGpsjAELs6KaC+mDbJ79p1ptm9+85DnNRUxzR4wUVo\/TMal8c+VabUmRgQA2JgrhnovsI9Lsq0n2sCzaAFF1FFYIs2EXbqi+pkEwNiQAyIgdgM9PPTjkXHco7t2MtmkKp90IKbKqkatxgQA2JgfAx0v6Rct8JFXBFZB8qOP2zksZH3gY26VDB6wb1pYcFsbCGVjaJxigExIAb6YIBX8vlo1weXhENEl\/blPSKeDr4XNlcm78PZU0mnI7gzZOPDS0vmfcvLGa5ZXJzhSNS1GBADYmC2DIS9ni9EbGlbLbizZWY6vY9ecO+zumoeuLKS4a56PiR\/BwAAEABJREFUPd909ir1IgbEwCAZCHs9H4KJcHYF7QdJy9QGNXrBPW9jwzicNczX801tY6sjMSAG5psBvZ5vttt\/0ILLtXpukBdRhJ36ojrZxIAYEANiIDYDrFC7rm5pR\/vYY0or3qAFFyoRVcQ1D+zUCyNnQNMTA2JgIAwgmiGQ4A5ecNnTENc8sAtiQAyIATEwLQYQTAluCNtJCG7IBNVWDIyYAU1NDEyRgRCxpS2CPcXhDrCr0QvuuZZ0\/hgoFmy4wX9iv42GeLEnfdujjImNS+9uTGxcaCceExs2Xmws2phDhx1i9A+H79iIPcj724CxYUPO6NPPL03NaDIz6Xb0gjsTVtWpGBADs2dAI4jMAKc3rFS7gvaRh5RYOAluYhtMwxUDYkAMzIaBrkLr2klwkxDc\/BPKlGezw6lXMSAGxMAoGOgwCQTTiWeXlPYduh1Rk8ELLuKaf0KZMvYRbQdNRQyIATEwcAa6iKzfRoI7eMEd+B6o4YkBMSAG5oQBBNMX0LZ52idCVU\/DTFZwWeU24eTowoL52haa+MtHDIgBMTBWBvw3BZFvN8+2Apv3l+D2Lrihl34RVmI4tNtBjPmnpSVz1fJyhk8sLrZtLn8xIAbEwGgYCHtbEIKZF9E2ZdqXU+mO8aRFXtiLUOQ7VFvvghtj4oiugyP81LjlpQetrpq9KysZLlpfL3dUjRgQA2Jg5AyEvS2ojbgW+ZYLLsd1d4wnpZzfFNjzyPsMvdy74EJQEXldiSEeaBpz58aGcThbbwvqSrvaiQExMAIGwt4WhGAWCWlTG+1PJ5FjOcd0v4Yydt+Wz1OPX94+5HLvghs6eUgNjTHv7TV\/MSAGxEA4A\/qlqVAOexdcJ5ikeYQOXu3FgBgQA2JgWgw0Xcme7nf48NftIItXuLZibj69Cy5L\/jI0YZm2eaGmjL1Je\/mIgWoGVCsGxEAzBhDM08XUmHrboUMLZv\/+OzTrpoFXqhrQu+A24K7WBXHNo7aRHMSAGBADYiAiA\/XCWia++\/YdNQcO3BhxLGmGSkJw06RWoxYD6TOgGYiBkwx0X+Hu2nWz2b37ppOh5jQ3NcHlEgCAZ5eS7xt\/bjt4ZUTYUIP\/3GpHGBs2ZNTP9rcaExtXH7vJxMaX7KxjwoaL\/rnORhw67BDn8vNOO+vYsCFn9Om+wj2+8kWwTx86Vy\/zmkAZ++nexlTVFfkPyTYVwc0TBJHYhkSExiIGxIAYaMfAnHm7557QzS6gfQllThPQBUDZuVJ2+dTT3gUXsnzyuhJGnDJ0jal2YkAMiAEx0JABBDMEx6r7QSccfE9sVWW\/buj53gU3JgEQX4SYfSiWGBADYkAMFDDQclVrivwLws6TKSnBnacNo7mKATEgBgbFAKvbIhFtaqP9oCY0\/cH0LrisSLkU7E+NMnbfprwYEANiQAwMmIHQJzF55mrA05vG0E4Kbo+9Ia6ILF2QUibfFrTNoy7G9oUF41Dnq3oxIAbEwJgZ4JV8PlrNlRVq09VskR\/tW3U4PuepCC60IbIOlLvAtffTujjnLS2ZOy8vZzhncbHOXfViQAyIgdEyEPR6PgQzBDUPTY2WdG9ivQsuK1KvvxPZMvsJh0iZG1ZXzZGVlQwb3V\/PF2k0CiMGxIAYmB0DQa\/nK1q1trXNbuqD6Ll3wZ31LG\/Z2DAOt+n1fLPeHOpfDIiBGTIQ9Ho+VrdtBdb3p\/0M5z6ErnsTXFawgEmS5sFlYeqEETCgKYgBMTB+BhDMEOiSsulNcBFUwF5Imgd2QQyIATEgBhJhwF+tds0nMtW+htmb4MYeMIIdO6biiYEpMKAuxMA4GGB121VoaUf7cTDReRa9C66EsvO2UUMxIAbEwHAYQDBDoEvK\/V1SHs5eopGIATEwWAY0sHQYYJUainRm28tIe1\/hMmoemCJ1yJedvU3aNMYn7AyPRsQeO8iYeLqNFxu3szFj43wbMyZsuOif7ZPbm9iIPcjY24V4sceYSryP2IHGhg0Z9XOOjRYbNuRsPvxSVMivTSHWsxn5YHq1UtTvWBDG\/GVlyti79kxbYnRtr3ZiQAyIgQQYGNYQEcxQDGtGUx9N74I79RmpQzEgBsSAGIjPQMj9W9rW3MNlIeVQNXjnQ1rlN8S65AQXkrW6HeKupDGJATEwagYQzZAVLu1LCHLHdY7tgHKRK3bqM1xzjaFc5DdUW3KCO1QiNS4xIAbEwKgZCBFb17aAIEQTAfWrKGP3bZSx+7Z82a8bYr53wYUQiPInTxm7b+srf3hhwTj01YfiigExIAZSYMB\/UxD5VmNmhRqCmkvKTceCfjg0bTMUv94Fl4kirhdffHG2\/Icoytingf0XLJm9Fy5nuPKOelvQNDhXH2JADAyTgeC3BbmVasv08LULxiDWgbQ4\/UBDAOXAkFNtvm1avUGOw7T6pJ8rrls1B69dybDv6DomQQyIATEwlwzM6m1Bh\/5ph9n\/oqVgztGQ4CAzDDA1wZ3VHHdvbJg9W9g14LcFzYof9SsGxMD8MBD8tiBWqR2w73+umwPfd2R+iC6Z6VQEl2V\/EUrGJLMYEANiQAwMjQGEtuWlZLPlv+u8TbP7nhtDm9HUx9O74CK0XAYAzC6fYmsD2hOzTRv5ps6Axi8GxMDMGQj5lSna8ktVBZMoOqZzjMfuu1PG7ttSy\/cuuEWEjIG4onnJJgbEgBgYLQMBK9xspUv7EnKcJiCogLJzpezy2Ck7UHZ1KaQzEdxQYiAZwkPjqL0YGAMDmoMYmAoDW5eHM\/Hsmq8YKMd1B98NW76MDfj2FPJJCi7Epkg24xbEgBgQA0kywAo1BJH+DjdJ7rYG3bvgIoxuNerylMlvjaF1QvumjT5gd5D3RAR\/yRsTX7QTiQ073exP3mKma3acMfEQGy82bMjBf7iVFRuDn3RPA\/w+G7cexrTxsSGjfh5qo8WGDTmbDweUritb2tF+NiMfTK+9Cy7i6IsreRDCAO2JGxJDbcWAGBADYqAFA4hmKFp0N0bX3gV3jKRpTmJADIiBMTDQag6sUEOgS8qmd8HtazXaV9xWO6CcxYAYEAPzwgBiG7LCpf28cFUyz94Ft6RfmcWAGBADYiAlBkLE1rVNab4m\/mB7F1x3r5U0j\/jTUUQxIAbEgBjohQFWqCHQJeXpXFLm8m8RetkpckH\/a2HBOOSqVBQDYkAMzBUDvJLPR6vJ83j9LbZFV9DeNp\/nT+8r3CmSW9jVi5eWzHOXlzO8eXGx0EdGMSAGxMA8MBD8ej6tcIN2k9EL7tNWV81zV1YyPGRdr+cL2lvUWAyIgaQZCHo9H2Lr7sV2SWmfNHvhg5+K4Obv3bpy+PDrI\/z3jQ1zyRbuPM+v56unSh5iQAyMnIGg1\/N1Edl8m5HzWze93gUXcXX3bxkMeT8l3wUubpe2aiMGxIAYEAMtGWCFGgI9NNX\/Q1NFmxTRRTCL6mQTAzNgQF2KATFQxwBim1+xtinTvq6Pkdf3vsLtgz\/EGtHuI7ZiigExIAbEQAEDbcS1zLcg7DyZkhRcie087aKa68wZ0ADEAAywQg2BLin3f0kZcWRFyvZyecrksQliQAyIATGQAAOIbdnKtYmd9hXTRBccitxcXT4t8h2qbSorXF9cyYNpEfJE29HDI+JtNtbQ0WTfn7XP+yyPsWFD6jNHDMR+rSXxOtBX2aSPY8Vtx46ZWHjnu99dOf5TKhHMUJwS8GQBEUUXHCifrD2Zc\/V+erJ2+LmpCO7wadAIxYAYEANioJKB0LN0xLqgA8QVAfWrKGP3bWPIT0VwIa4IYyBQcxADYkAMzAUDm3aW\/DxjG2z5Hl5bMIb2NkTIx9eRkDizatu74EIQZytF6DJp4nVppzZiQAyIATEQwEDACvfQNTvM\/ncvBXR+vKmvIylqQe+Ce5wm\/S8GxIAYEANJM8Al4Y6iu2953Ry45EjQ9BFbPwDl1EQ3J7j+dMaR7\/xmjHFMX7MQA2JADJxg4PDhw8bHiYomGQS3I3adtWl2n7\/RpJdR+\/QuuDHOQjiLcWBruDwp5SpcsLRkLlxeznDHxcUqV9WJATEgBkbNwKFDh8zevXtPYP\/+\/c3n23F1a1w7xLp5b6d5Njnen9ZoYIboggspeTDnvI0y9iZAtB3wd3lSylW4bnXVXLuykoE3ZVT51tWpXgyIATGQMgP79u0zBw8ePIEDBw6Yxv8QTCeeXVLaF3TGcTyvB5SxF7ifMDXxOeE8kEx0wYWkppgGBxsbG8aBN2VMo0\/1IQbEgBgYIgO7du0ye\/bsOYHdu3ebxv8QzFCUdIZmIKAOlJ0rNvLYyDtQxp4SogtuSpPXWGMyoFhiQAyMmoEuq1q\/DWJdQRAC6uC7YXNl8g7OllKanOBCdkoEa6xiQAyIgVEwgGD6Ato2T\/tRENF9Er0Jrlv2+0Mrsvn1youBsTKgeYmB5BlAMEORPAlhE+hFcBFWVqLAHx5lQL1vV14MiAExIAYGzgC\/GnWLHWNXRPilKdt70p\/ogouYIqpVrFCPX5WP6sSAGJgXBjTPJBgIXd3SPomJ9jfI6ILb31AVWQyIgSIGzrHG2DjfxoyJ2OMjnh3ifH72ToyJhQMPa85h23u2eX8Jbv\/vw22+NeUpBsSAGBADeQYGU0Yw8yLapkz7wUxmNgPRCnc2vKtXMSAGxEBaDCCYoUhrxtFHG11wm9yf5f4tftFno4BiQAyIATHQDwNtVrNFvoh19JGlFTC64DJ9xBRRJZ8HdurzdpXFgBgQA2JgwAwgmEVC2tRG+wFPbxpD60VwGTiiirjmgZ16QQyIATEgBhJiAMEMRULT7WOovQkug0Vcc8A8Vej1fFOlW52JATEwYAYOby4YH62G2nQlW+aHWLfqcHzOvQpuLLryq2S\/XNeHXs9Xx5DqxYAYmBcGDn1th9n75eUT2L+61HzqCGaZmDax0755b6P03JbCrPKrZL9cN\/6kXs9XNxnViwExIAYCGNh39ro5uHPlBA6ce6R5NH5pKgT6panx\/x2uezUfqV7P1\/y7JU8xIAbGx8CuhU2z58yNE9ht841nyQq1yUq2zIf2jTsbp2MSK9xxUq9ZtWRA7mJADMySAQQzFBXjb3Or0PlWhBtklQR3kJtFgxIDYkAMDIyBspVrUztiXTIlBNS\/VUi5xDVpc++CW0ZcmT1pNjV4MTArBtSvGOibAQSzqbgW+dG+YIxoAWLrV1HG7ttcHjv1rpxS2rvgpkSGxioGxIAYEAMlDCCYoSgJPS\/m3gSXsxAAkaR5pHqGwnwEMSAGkmJAg43BQNGqtaHt8K0LxiDWgeNAR1LWjt4EF1IA\/JLmgX0aeLbt5AURcZaNFRMX2nixsd3GjA0bMupn2UaLjT02ZmzYkFE\/sbcL8Y7aEcZGzH2cWBt2jLFxkY0ZGzZk1A+vEYyN7e81JhYe\/rEW0z1mfRHNDjh06w6z\/6YWf\/Nru8p\/Uhdb5tOb4BJcEANiQAyIgYEz0HR4DVezpsBv32TdHNje4m9+S8aE6AOWXHAAABAASURBVDrgQp40FfQuuKxsY5GRGrmx5q04YkAMiIGZM1AgpEXiWmTb9fVNs9twjaP7LNASH0SiTJoKehfcWERIbGMxqThiQAyIgQ4M8EtRPfzSFKKZP75Txt5hlLGbRI03FcGFvCI0nQltm\/rKTwyIATEgBnpgIGCFm616ufdbMizEleO8A2Xnis3lU097F1zIgrwiNCXPtW3qLz8xIAbEgBiIzEDAQ1PZE8oVgstI3XGelLJDvlxnd\/VDTHsX3GlPOt\/fjQsLZn0L+TqVxYAYEAPzxID\/ulLyrebe4wq31TgSdh694L51acn87+XlDB9eXEx4U2noYkAMiIEwBs7ZscNcaI+HDry+tHFECW5jqsocexdcLgdwWblsAH3b966umketrGS41\/p6390NPL6GJwbEwDwzcNQeA6+1x0OHtSMt\/lSn50vK87BdehdcJ7akeUyD4KWNDeOwY5PH7KbRq\/oQA2JADAyPAV5RyqtKHW62x8fGo9QKtzFVZY69Cy4r3DKUDUp2MTALBtSnGBADFQzw0FOI6NK+Ivw8VPUuuPNAouYoBsSAGBg9AwhmKEZPUvUEkxJcVsrV01GtGBAD\/TCgqHPPQMjqlraI9ZyTODXBdfdv4Zs8qSAGxIAYEAOJMMAjMD380lQis48yzKkILgLrr07JY4syg5ogN9j6tYjgRC0mGF9snG3nGxu8lSYmbrRjjI3\/ZWPGhg0Z9cMXLjYusSOMjS\/ZmDHBG4Ni4zo7xtiwITt9yhrttBWxYUPO5hN64NMK1\/Dd73XjIawIbK+dKLgYEANiQAz0ygB6GaK5tO91gAkE711wE+BAQxQDYkAMiIEaBhDMbjDZLzvyZ7w1XYy+OinBZbXsMPotowmKATEgBgbEQMjq1rUd0HRmMpTeBZfLyYikPzvK2H1bXd61oR2gXNdG9WJADIgBMRCHAVa3Tji7pLSPM5J0oxQJbvTZ+AKJUFJu00mXNm3iy1cMiAExIAaqGUAwQ6BLymY6D02xGRFZB8rTwtGFBfO1LUyrT\/UjBsSAGBgiA7whyEebMXZZ1ebbtOlvjL5TWeGGEodQ+zHarHj\/aWnJXLW8nOETi4t+mG55tRIDYkAMJMpAyNuCWN3mBbRNmfZVtHFcdyjzc\/UuLfMbqr13wUUsIScWAcQiZtN4D1pdNXtXVjJctL7etJn8xIAYEAOjYyDkbUEIZgiqLim74zrHdkA5Tz426nxgy\/sNudy74MacPORCdpuYOzc2jMPZeltQG+q6+KqNGBADA2Yg5G1B\/MjULXZuXUF72\/y0T9FxneM89tOcEzf0LriONNI82nBHWzZCmzbyFQNiQAyIgTgMhKxu\/9\/Cgqla4TYZ4RiO\/70LLiSVoQnJ+EhsYUFIngFNQAwkzECb+7V537fu2GEOLC1FmT164IC2RAk6pSC9Cy7ExJgLcfKIEVcxxIAYEANioJ4BVrh5IW1afsT6unnykSP1nTTwQGQd0IQGTQbj0rvgxpipIzefxoitGGJADAyCAQ1i4AwguF2xc3PT3H9jY+Az7H94vQsuIpnaWUj\/tKsHMSAGxEBaDDRdzVb5hcx4DDrSu+CGEByj7S4b5KKIsKGifs6w0WKDpwFjo+pL1KUu9qvaiBdzO7tYdvNE\/cTeLsT7lB1hbMR+pdxRO8bYYO6xYYcZ9RPtlZZ2VC6Wzc7k03V169qVPTRVtChDXLHPZKI9dtq74EIc4yfNA7sgBsSAGBADw2cA4exycu3a0L5sloirrw+UnS928tjI+8BGXSroXXAhpAxtSPJJJt+mrXzFgBgQA2IgjAEnnCFp1Qh8nfD9sLsyeR9b9mSS3gU3BhMIrE8yeWwxYiuGGBADYkAM1DPACjUEZZeU63sej0cSgjseujUTMSAGxECaDHC\/vOuvTNGO9mnOPN6oexdcVqJ5uHLTabCibeorPzEgBsSAGIjPQMjqlrZa4U7h9XyIZR7sCthI28AJNWnT9tcvLJgbttCmL\/mKATEgBsbGgP9qPvJt5odohty\/pX2b\/sbo2\/sKt4g0xBLRLKqrstHOoWn7NywtmVctL2d4\/+JiVfiB1GkYYkAMiIF+GAh5PV+I2Lq2\/cwqnagzEdxp0vPI1VXz+JWVDPdd1+v5psm9+hIDYmBYDAz19XzDYqm\/0SQhuE1Xs0U0LW9sGIdzNzeLXGRLiAENVQyIge4MhLyej0vCbqXaJaV995GPo2XvgotYFoFLw+OgULMQA2JADIyfgS4im28zfpaqZ9i74CKsRage1qm1tM+LNrZTvVQSA\/POgOYvBvpjgBVqCPSU8hSeUo61+RFYH7HiKo4YEANiQAzUM4DY5lesbcq0r+9l3B69r3DHTZ9mJwbEQCoMaJxhDLQR1zLfsBGk37o3wXWXgH2Kimx+fQr58+0gY+Irxx5uYqNsZw+x22lH\/XzZRouNH7cxY8OGHPwnZLuWtV21s44JGy76h7erxkbsQX7GBowNG3ImHx455deiuoJ9bSYDH1CnvQguwuou\/\/pzdTbqfbvyYkAMiAExMEsG6vtGMENR38u4PaILLmKKsFbRRj1+VT6qEwNiQAyIgeEwwD3YEOihqYQemhrObqeRiAExIAbmjwHENmSFS\/sq1liEOcTwq4rRti6Wf\/QVbqyBKY4YEANiQAwMh4EQsXVty2aD0HLl04FykS9250NKuchvqLYkBBdSyzBUYjUuMSAGxMCYGGCFGoKyS8oc2xFPnyvK2H0bZey+LbV8dMGFEIipIoJ6\/Kp8\/Dp8y+D7nch7Gb0tyCNDWTEgBuaaAd4Q5KMNGYitW6m2Tb+6sGBo36a\/MfpGF1xIQhwRVfJ5YKc+b++rrLcF9cWs4ooBMZAaA7N6W9BHd+wwr19aCqIrrxvT1pKgwW817kVwiQ05EJIHduqnBb0t6DSmZRADYmBOGZjV24Lus75uHnzkSDTW0ZVpa0mMwfcmuAwOQvLA3haQW4a6WO5NQaR6W1AdW6oXA2JgzAzM6m1BOzY3zYUb\/ExJOLtoAboSHmn6EXoV3FjTgdwyxOpDccRAxoD+EwNioJCBtvdti\/wLA7cwpiy2TDMJwWWgghgQA2JADMyOgb5+2pHFFELqz4wy9jqbX59CXoKbwlbSGMXAMBjQKOaYgaIVa1tbGX2IKyLrQNn5YvPzlH24uhRSCW4KW0ljFANiQAzMmAH+rCcEZX+H66aFyDo4Gyk2l5LPg7pUIMFNZUtpnGJADAybgZGPDrFtu6L1\/Wk\/copqpzd6wX3yY4155pPi4Ut3MyYmzpi808QGO3ZsbLe7kmCMOBAHY9oH2ghAjGOKPYzM9acN33NNlCYvBsSAGJhnBvzVaoe8fmnK7jxJCS43yu2Ys4+fzwz6TwyIATEgBnpjgBVuF6F1bWjf2+ASCZyU4CbCqYYpBsSAGBgdAwhmKEZHSssJlQpuyzhyFwNiQAyIgREz4FaqXVPEesT0NJqaBLcRTXISA2JADMw3AwhmV7GlHe3nm0FjkhBc7tcCNhYpcHnSKhzeWDCHbz6OKr92dfIWA2JADKTHwPaFBeOjzQxutc4h4JeqbIi5\/iQhuO4PndlS+Ty2Kuz\/+JLZ+6HlDFd+YbHKVXViQAyIgVEzwOv5lpaXjcPOpeavzGOFykq1K2g\/anIbTG5bA59BuLCqRWwZjJ+nXIUr7rlqDt53JcO+C9arXFUXiQGFEQNiYJgM8Hq+1ZUV43BDi1fmIZihGCYr0xtVEoKLwEIJKcjnKZdh9x03zJ4t7DpTFzXKeJJdDIiB8TNw2+amuXlj4xQ0nXXXla1rh1g37WusfkkIrlvZkgI2BikgL4iBdBjQSMVAmgwgmE48u6S0T3Pm8UadhODGm64iiQExIAbEQBcGEMxQdOl3TG0kuGPampqLGEicAQ1\/uAx0WdX6bRDr4c5uOiNLRnD9y8d+fjo0qRcxIAbEwHwzgGD6Ato2T\/sqBnk+x6HOr6p+yHXJCO6QSdTYxIAYEAP9MzDbHhDMUJTNAKFlIeVAuci3zF7kO0RbkoIL6T4qiX3pq43587dEw70+b0xM8JfBsVHJR8fKtmez8jdGHIiDoe8DCGjTQ0LoXMr64liO0PrjoIzdt+XLfl0q+SQFF3LZIA6UBTEgBsSAGOiPAf6ocpa\/NFV3vO9v5vEiJyu48ShQJDEgBsSAGKhjIGSFe9PCgt6HawmW4FoS9BEDYkAMiIFqBrgk3FV0v7pjh\/lsi5+RrB5JurWzE9x0OdPIxYAYEANzxwCC2xXnrq+b81v8jORYyR294B4+fJ05fPirGca6ETUvMSAGxEATBhbspV0fTdo4n9tspiu2bW6aszY2bIT5\/oxecPfvf7nZu\/fyDFdeeVUqW1vjFANiQAxEZ4C3BV24vGwcLmhxmZfVbVfBpR3tiybEw1D5J5ApYy\/yT9k2esG94oonmoMHn5Vh377dKW8rjV0MiAExEMQAbwu6dmXFOKy1uMyLYIaibPCIKyLrQNn5YnP51NPRC+7u3fc0e\/Ycx65dd0p9e2n8MCCIATHQiYFNe2l3w17adeDNQU0DsUoNAWJd1Rci6+D7Yasq+3VDz49ecIe+ATQ+MSAGxEAKDCCYfQpuChyEjlGCG8qg2ouBYTGg0YiBXhhAcEPRy8ASCirBTWhjaahiQAyIgVkxwK9M3WI77wp+qco2n+tPsoLLjXSHud6CmrwYEAPNGZBnZwZCV7e079z5SBomKbjcRPcxkm2haYgBMSAGBstAyP1b2kpwjUlOcN2q1qV1e+d\/u\/vTzPbJo6PhM7bDmLjOxosNdu7YsMOM+rntocbExvl2hLFhQw7+82o7wth4no0ZE5faeLFhQ+oTl4HKaMdsLaIZAhtirj\/JCa5WtnO9v2ryYkAMzIiB0JN4hHpGQx9Mt8kJ7mCY00DEgBgQA3PEwNwLboRtLcGNQKJCiAExIAbGzoAuKYdvYQluOIeKIAbEgBgYPQNa4YZv4hELbjg5iiAGxIAYEAPHGZDgHuch5P\/RC27XV1GFkKq2YkAMiIEhMuAfD8m3GaMuKbdhq9h39ILL66fcq6juuLhYzMIcWjVlMSAG5o+BkNfzdf2FKddOvzSV4N\/htv2KXLe6euJVVLyaqm17+YsBMSAGxsIAx0D3aj7SNq\/n0wo3fC8Y\/QrXvYaKlFdT5SnjsgorX9J8XdcysWLGjB2PeaUQ8\/DmgrnyhkVDyphjYLKwYM5YXDSk5fGa16TAI7NZs\/P+x8VFQ0o5FEdtvI\/ZeKShsVz7r9uYN9mYpM4WmsbePrHjMb9pxuQYyLHQYUiv54OLsWPb2CdYNz929sXzzqtza1UfO2bseEwmhZiHb0Nw426bbfagfkbE7Z0Cj2xvhPYfI8573fL4fyLGY4wI7UbkmLG3T+x4zDuVmPxwRciDU7RnvvOMuRfced74mrsY6MKA2swnAwhmKOaTuZOzTk5w3W8ok56chnJiQAyIATHQJwMhq1vaItZ9ji+F2MkJrv9byuTrcf9wAAAPPElEQVTrSOZyTRVc+yqftnWxY8aOx3xmHZP7snVwY+TScp0v9dyXrYOLyaXlOl\/q4aoKLl6VT9u6tjG5XFwHF7POz9Vzb7YKLh6Xlqv8XB2Xi+vgYtb5ufomvLqYTXyb+NTHWzBN4vg+s4zp+m6S8n2wkzNdIcEd8VPKu3btMrt37zb+nwW5Pw\/yU+rZ2Uh9e0ieWDFjxo7H3GYdc++Xl00d9q8uQaMhrfOl\/pzlZVOH2y8dj0la50s9XFVh1jwytpfZedfhdVvzJq3zpf5NNmYV3rEVj7TKz9Vdb+PV4catmKR1vtQz9zrE3j6x4zH+Wcakb46THC+zL1vBf9Ths3NpySzZ7dgVO5eWzF3uchdDvIJu5sK0bayzZKNeccUV5uDBg4I40D6gfUD7QMk+wHGySgdiHkvf8IY3VHU1+rrRCi5bjh1lz549RhAH2ge0D2gfKN4HOE5yvKwCPjH4I05VP2OvG7Xgjn3jaX5iQAyIATGQDgPVgpvOPDRSMSAGxIAYEAODZkCCO+jNo8GJATEgBsTAWBiQ4HbfkmopBsSAGBADYqAxAxLcxlTJUQyIATEgBsRAdwYkuN25U8sqBlQnBsRAYwZuvPFG89u\/\/dvmc5\/7XOM2VY633XabOXbsWJWL6mbAgAR3BqSPuctbb73V\/M7v\/I55y1veEmWa\/\/zP\/2x+4zd+I0osBREDQ2QAsX3e855nzjzzTHPhhRcGD5F4l112mbnyyisN38fggAoQjQEJbjQqTwa66aabzK\/+6q+av\/7rvzZDPsvkLHh9ff3kwANzfLn5I\/r3vOc95oEPfGBgtOPN73CHO0Q5aDDX1dVVE3O+x0fY6X81SpiBW265Jft+X3311cGzcMeKc8891zzzmc80t7vd7YJiIraI9zd8wzeYD37wg+YVr3hFlO9P0KDU+AQDcy24N9xwg3nBC15g7ne\/+5lLLrnEPOtZzzJf+tKXTpDTNXPWWWdlZ6q\/8iu\/Yg4dOhRFdBGKr3zlKwbRQDy6jo12COOrXvUqc+mll5r73\/\/+5tGPfrT58pe\/TFVnEBOxff\/732\/+6I\/+KJt\/52Bew\/POO89cd911wUL51a9+Ndu+3\/3d322++Zu\/2fz0T\/+0+c\/\/\/M8o28Yb7iCzR48eNX\/wB39gOBjHGiD74y\/+4i+ae97zntk+9Bd\/8RfBXPoxH\/SgB2VXSWKesPLd+fjHPx5MwRlnnJH9mM7v\/\/7vm0984hPmIx\/5SOeYt7\/97c2P\/MiPZOL43ve+t3Mc1\/Av\/\/IvzT3ucQ\/zC7\/wC+aP\/\/iPzYc\/\/GGJriNnAOncCi4Hn5\/6qZ\/Kviy\/\/uu\/bp7znOeYf\/\/3fzcf+9jHgjfLZDIxT3\/60w2XdUJFl4PQr\/3ar5kHPOAB2Zec3zR9yEMeYt7xjnd0OsAhjKy+\/+qv\/sq87nWvM3\/\/939vnv\/855sLLrig87yJ2YfYMqCdO3dmZ\/2hJ0LEYc4c1P71X\/\/VPPnJTzaXX365+fmf\/\/lgMUcUEHTGC\/JlbKEIiXnOOeeYz3zmM9kKiv0+xlhe+tKXms9+9rPmzW9+s\/mhH\/oh85u\/+Zvmne98Z+fQzI+YN998s3n7299uHve4x5lf+qVfMmyrzkG9huw\/z3jGM7Ix05dX1TjLiQsn6IcPHzaPfOQjzTd90zdlJ6v\/9m\/\/1um76Domzkte8pKMw3e9613O3ClljgcOHMi+M1wdYoUr0e1EZS+N5lZwX\/3qVxsuDSE6j3nMY8zTnvY0w+qMFZCx\/zY3N+3\/3T+TSbjoImQcyFiJ8UXk7UhcJvrO7\/zO7OD5J3\/yJ62\/6B\/4wAcM90VZ8SDirFAQ8clk0mmyrLa5V\/S3f\/u32T0jdw+KgxqX1BFzDqKdgttGrAC+7du+Lerl+e3bt5tv+ZZvyX5b9973vnd2YhQiRDzosn\/\/fsM2Yt5c1eBkyxdhO5XOHxfz2c9+tuGqTNtA7EcbGxvZKorLliFzpW+E55Of\/KRhZXaf+9wnO1ndu3evedOb3mS6fm84sSTm4x\/\/eHP3u9\/d\/MzP\/Ex2WyJk9chYAWLLyfVTn\/pU873f+71mMum2r3Piwn5zuT1RY6X8d3\/3d4ZLwcvLy51jMj4QU3SJ5yDRdUwMI51LweXLzZkzl1LZIfObgi8Twtv0Hg33YThD5WEh2jiBmUzCRJfV9oc+9CHzohe9yPClZpys1Fjx\/uzP\/qz53d\/93dariv\/4j\/8wHCQ5qBEPcEBC2DmYcIm5zeVBxOsRj3hEJvyMh4O5EwjKP\/ADP5A9DEI\/XUF8VlNve9vbuoYobMfYn\/KUp5hv\/dZvzR706ioWcPlbv\/Vb5oUvfKHh\/hknGpTvdKc7Ffbbxui4JCbbiAN8m\/b4vvvd787EiwfZ2Nahost9xrPPPjs7wUDMKXN7ghMMhJ0+22Lbtm2GV9a9733vy\/YlYu7YscNwYswJ4Ytf\/GLD98y0\/Md882LLpWVuKbBPtQyXrWw5SeUk\/fu\/\/\/uzFf5rX\/tac9VVV5nHPvaxhtVu25jOv6HoOvfGKcc4rXQb09Wr41wKbh2jd7vb3cz5559v\/vRP\/7TRGTuX67hcyarxCU94grnvfe9rHvWoR2UH8U996lPZSoAVT9vLy1dffbW5853vbBBZf8yTycT86I\/+qOFLT5+InF9flWcFyiUmDsJvfetbs4PEwx72sGz1wz3sJz3pSdlKFZ+qOH4dBwruF3Gw4d7Ra17zmuxkgC\/5pfY+se\/bJc8Bg5MMDmyszrvEKGszmUyyy8tra2vmox\/9aJlbrR0OnvjEJ2b37Envete71rapc\/DFlsutXWNyRYTV3b3udS\/zh3\/4h9lzCiGiy9O0rG7\/5V\/+xXz605\/OpsF9w2uvvbaTKBKA1SP7ntvXiTWZTMzLX\/7ybLX7+te\/3nDbAk7wb4Iysf25n\/s5w4r8oosuahLmFB+uZvDd4cSH9t\/4jd9oGPdll12WXTUJ3d\/Zjzh55+SKKyandB5Q4DvE95HvNSknSgHh1LQjA3MpuHy5ER4uCbHazXPHzsl9qS984QuGy2f5+nwZgeVAxoqGs95\/+Id\/MD\/2Yz9m+HLu27cvE2DEDfHkEiuXHJscOHj4irN6Ltvm+2QFgLizqvjiF7+Yry4tc++JA\/BP\/MRPZA+McbmWAwj34n74h3\/Y\/ORP\/mT2IAziWRqkoIIDBaLLipxLzC972cvMpZdeWuDZzcT24kDEgyrcb2d10oTDJr1x2Zp7htzPbuKf92EcbFO2O3NnZR96sHQxWdmGiG1+rKzO2FcRoxDRZb9hjlySd31whYBxUyb+7\/3e77W6P85+yT44mUyyh+4Q2z179mT3iH\/5l385u5rD6pT4TXDkyJHs4Th3GZm2iC1XtvheTibtLy3zfeP+P6\/9fK1d2XKLBgHjKkeTY0WTcfNdYl\/Piy4i2eWWguuT4xpjHbXouskONJ1LwZ1MJoZLnVw65ovjDhL5bcQOymWuvL2ozJOvHMh4cISndLkkjThwWZgHdRDHiy66yPCEIweCohh52yWXXGJ4QKNM\/HjQictuRYKcj+XKCDWXPlmVc7Dg0jRi5uoReL7U\/oHU1dWlHCgQXQSMqwNtVt51sann1V6sdOD2uc99bnYpmJUKlwdDD3YIEVxff\/31dNUKn\/vc5wxiizAyNtI\/+7M\/M023c1FnfsyuK9uiuNiYK\/sqothVdCeTSbYvEy8P4vJMBHZWw6Sh4KoT+\/pk0lwkORHm6s1kMsme7g8VW+bA95wrGHw\/WNmy0mcf5MSI706s2x58l3zRRWxZ4XMy1+b7zph9cEyT6PqMTDe\/bbrdDac37ttxb4edmh2QHdqNjgMGB6THPOYxhtWws9elfBlpxxPEPM3IF3AymZi73OUu5gd\/8Aezp4IRYJ5gnkzqDxzcT+X+JfdwuTSd7\/\/zn\/+84X7x4uJivqpxmTG6VT4pfXG2\/uAHP7hxDN+RAwWiy0kCfzYSW3S57\/rt3\/7t5o1vfKPhhIGV5eWXX264p+iPo22eS4TEgM+2beGLkywnjKzsEVz+pKltLOefj+nssVInuisrK+Zv\/uZvgsO6S8HcBkFsuZKCmHOC1yU4+w23OPguse\/zPeVeLleJ2sbjZIwVcsjKNt8nJ1PcOmE7sc+cddZZBkHnu1h2Ap+PUVfmu8S8Wemy8uf7yYkm34G6tlX1jJWTef7MkKsSVb6qi8vAjAU37mTaRJtMJtmTydxzZCXKAeKVr3xl9vNqCC1nrXxB28TEt0h0sXcBByue1uRgxtk0K2X3ZWY1xhfxu77ruwwrvy7xeUqbe6Pf8R3fYXgC9nu+53uyVRlx6btLTNpwoOhTdOkDcOCBGy7lTyb1JzC0KcNkMsmeOPVPvMp8x2JHdLmtwiXS0DlxssIDU4gk36UQsWUsk8kke36BWx88D8HJJ7cSJpP225mTZoSr62VkxpMHJ1M80MUDhpygc4WDE66HPvShwU8s+33xXWLs\/EUBJw1cPfLr2+bZv1kpI94x4rXtf97951Zw2fAcsHn4iC\/Lwx\/+8OzPGvgTHFa8rM66ik5M0UVQeDCKFflTnvKU7GlT7p8xXu5vcTl1Mml\/EGL+XN7mqVrEnMt1\/H0uf2rEZSfqQ8CBAtGFC86oQ2JNoy0HZR5q8y+vT6PfWffBPjCZdNt\/\/LFfZG+XcOKGQIaKLXHd\/sjqFrCfhogN23cyCZ8nY3PgRI8VNyetXDp3J8OuPlbKd4mVbcj8\/bFw3JDY+oxMLz\/Xguto5k9u+EK\/\/e1vNzwIgbhNJmFfToSGy8s8WMKfDLi+2qTOl7NpHiDhEipnpwgDcfnRhq4nBS42l1K5tM4lcC7VchLi6kJTDhSc0ISOMXQcat8\/A2xj9iH3owv99zicHjhJ47szmYQdM6YxI7YTY40l3tMY85j6kOD2uDURXYSRv8+L0Q2rXVYRPM3JlzxGTMUQA7EY4GRtMhm+6MSar+KIgbYMSHDbMtbSv49LWS2HMFJ3TUsMiAExkBYDEty0tpdGKwbEgBgQA4kyIMFNdMNp2GKgigHViQExMDwGJLjD2yYakRgQA2JADIyQAQnuCDeqpiQGxEAVA6oTA7NhQII7G97VqxgQA2JADMwZAxLcOdvgmq4YEANioIoB1fXHgAS3P24VOYCBiy++uLR1VV1pox4qGIePoi6oL7LLJgbEwPwxIMGdv22ezIyHLFaM7ZprrjE+sCVDrgYqBsRABwbCmvx\/AAAA\/\/+YSj0wAAAABklEQVQDADpZVC3RkaKxAAAAAElFTkSuQmCC","height":229,"width":381}}
%---
%[output:235d5d03]
%   data: {"dataType":"text","outputData":{"text":"[23:31:33][INFO]  Generating 300 SMILES with bigram Markov model ...\n","truncated":false}}
%---
%[output:4a72d9ff]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 30\/300) Generating (Markov)\r[##--------]  20% ( 60\/300) Generating (Markov)\r[###-------]  30% ( 90\/300) Generating (Markov)\r[####------]  40% (120\/300) Generating (Markov)\r[#####-----]  50% (150\/300) Generating (Markov)\r[######----]  60% (180\/300) Generating (Markov)\r[#######---]  70% (210\/300) Generating (Markov)\r[########--]  80% (240\/300) Generating (Markov)\r[#########-]  90% (270\/300) Generating (Markov)\r[##########] 100% (300\/300) Generating (Markov)\n","truncated":false}}
%---
%[output:353a9a4a]
%   data: {"dataType":"text","outputData":{"text":"[23:31:35][INFO]  Validating Markov SMILES with RDKit ...\n","truncated":false}}
%---
%[output:7b36f722]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 30\/300) Validating (Markov)\r[##--------]  20% ( 60\/300) Validating (Markov)\r[###-------]  30% ( 90\/300) Validating (Markov)\r[####------]  40% (120\/300) Validating (Markov)\r[#####-----]  50% (150\/300) Validating (Markov)\r[######----]  60% (180\/300) Validating (Markov)\r[#######---]  70% (210\/300) Validating (Markov)\r[########--]  80% (240\/300) Validating (Markov)\r[#########-]  90% (270\/300) Validating (Markov)\r[##########] 100% (300\/300) Validating (Markov)\n","truncated":false}}
%---
%[output:8098b4ea]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]  === Markov Baseline Results (300 generated) ===\n","truncated":false}}
%---
%[output:0fd21ee6]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]    Validity:   27\/300  (9.0%)\n","truncated":false}}
%---
%[output:4f3967e5]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]    Uniqueness: 18\/27  (66.7%)\n","truncated":false}}
%---
%[output:9d3c70f5]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]    Novelty:    25\/27  (92.6%)\n","truncated":false}}
%---
%[output:06b087e6]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]    (Production targets: Validity>90%  Uniqueness>85%  Novelty>60%)\n","truncated":false}}
%---
%[output:8f440e59]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]  Sample Markov SMILES (first 5):\n","truncated":false}}
%---
%[output:7147febe]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]    [ 1] C(=C)C(N)CC(N(Br-C(NCCC(C(O)C(=O)CN1)CCCCCC4(=CC)CN)CC(=O)[CCCC@H]2N)C(CCN(C21)CO)CC)[NCCCC@H](=O)4C\n[23:31:45][INFO]    [ 2] N(CC)N1CCC(N)C(CC(C(C(C1N(C(OC)CCCC)=OCC)CCNCCC(=O)(=O)CCCC(C)CCCCCC(C)C[CNCCCCC@H](C(OCN)(O)[]))))\n[23:31:45][INFO]    [ 3] C(CO)C(C[CCCCCC@H]2CN1CCCCCN)C(C(CC)C(=O)CCC(C(=O)CC[C@H](=C[CCC)C@H]1C[CCCC@H](=C(NCC(CC(C))))))2C\n[23:31:45][INFO]    [ 4] C(=O)COCCCCCOC[C@@@@@H](N(N)CCCNCC)CCC[CCCCNCCCC@H]2CN1(O)C[NCCCCCCCCC@H]1CN1=CCC(=O)OC(=O)CC[O]12C\n[23:31:45][INFO]    [ 5] CCO\n","truncated":false}}
%---
%[output:2a8966f4]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]  First 5 valid Markov SMILES:\n","truncated":false}}
%---
%[output:067fae5f]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]    [ 1] CCO\n[23:31:45][INFO]    [ 2] CC(=O)O\n[23:31:45][INFO]    [ 3] CCCCC\n[23:31:45][INFO]    [ 4] C=O\n[23:31:45][INFO]    [ 5] C(N)O\n","truncated":false}}
%---
%[output:8b7375fb]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]  Markov Failure Analysis (273 invalid \/ 300 generated):\n","truncated":false}}
%---
%[output:5276edff]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]    Parenthesis Mismatch '()':  8 \/ 273  (3%)\n","truncated":false}}
%---
%[output:83699159]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]    Bracket Mismatch '[]':     6 \/ 273  (2%)\n","truncated":false}}
%---
%[output:14d3854d]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]    Ring Closure Number Mismatch:  17 \/ 273  (6%)\n","truncated":false}}
%---
%[output:9f6c2440]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]    Others (Valency\/Atoms\/Invalid): 242 \/ 273  (89%)\n","truncated":false}}
%---
%[output:8ce954da]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]  Next Steps: Build LSTM character language model in Sections 6-8.\n","truncated":false}}
%---
%[output:07691b4a]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]  Sequence matrix: 500 x 99  (vocabulary=37)\n","truncated":false}}
%---
%[output:0c181dc6]
%   data: {"dataType":"text","outputData":{"text":"[23:31:45][INFO]  Split: 400 training \/ 100 validation\n","truncated":false}}
%---
%[output:3b3443ea]
%   data: {"dataType":"text","outputData":{"text":"[23:31:46][INFO]  Network: Embed(37->16) + LSTM(128) + FC -> softmax  |  Number of parameters: 79621\n","truncated":false}}
%---
%[output:5f85d3e4]
%   data: {"dataType":"text","outputData":{"text":"[23:31:46][INFO]  Training token to parameter ratio: 4.39  (>1 => overparameterization)\n","truncated":false}}
%---
%[output:81098cf3]
%   data: {"dataType":"text","outputData":{"text":"[23:31:48][INFO]  Epoch   1\/150  train=3.4883  val=3.1569\n[23:31:53][INFO]  Epoch  10\/150  train=2.6523  val=2.6592\n[23:31:59][INFO]  Epoch  20\/150  train=2.6141  val=2.6346\n[23:32:05][INFO]  Epoch  30\/150  train=2.6089  val=2.6267\n[23:32:11][INFO]  Epoch  40\/150  train=2.6055  val=2.6201\n[23:32:16][INFO]  Epoch  50\/150  train=2.5971  val=2.6182\n[23:32:22][INFO]  Epoch  60\/150  train=2.5970  val=2.6155\n[23:32:28][INFO]  Epoch  70\/150  train=2.5932  val=2.6101\n[23:32:34][INFO]  Epoch  80\/150  train=2.5960  val=2.6174\n[23:32:39][INFO]  Epoch  90\/150  train=2.5968  val=2.6149\n[23:32:45][INFO]  Epoch 100\/150  train=2.5947  val=2.6161\n[23:32:50][INFO]  Epoch 110\/150  train=2.5913  val=2.6131\n[23:32:56][INFO]  Epoch 120\/150  train=2.5963  val=2.6128\n[23:33:01][INFO]  Epoch 130\/150  train=2.5929  val=2.6152\n[23:33:07][INFO]  Epoch 140\/150  train=2.5956  val=2.6149\n[23:33:13][INFO]  Epoch 150\/150  train=2.5966  val=2.6132\n","truncated":false}}
%---
%[output:84c89218]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdwAAAEeCAYAAAAgg6RKAAAQAElEQVR4Aeydv4tkSXbvv70MaNaTDCEtzGOiB9aW2S2nFSPQXyBLI9BEI8keS9YuVIzdIHi2hDp2EeM9W2bHChlPoD9gQTB9hm1jxDPeM54h0ECrvnn7VEbdujczK3\/cvDfzW8zJOL\/iRMQnMvN0VnXX\/Oi9vkRABERABERABE5O4EfQlwiIgAiIgAiIwMkJqOFuQqyYCIiACIiACByJgBrukUCqjAiIgAiIgAhsIqCGu4mOYpsIKCYCIiACIvAIAmq4j4ClVBEQAREQARHYl4Aa7r7kNE8ENhFQTAREQAR6BNRwe0BkioAIiIAIiMApCKjhnoKqaoqACGwioJgIXCUBNdyrvHYdWgREQAREYGoCq4b72WefoZWhTTA+5L8235Qctq21Ld7eDXNb2RRr81xnPnWO22TXvG11tsW3rbMtvq3+qeJz3depzjtWd2kcNu13U2zs\/KP+Mwd4FpddtsLcXfKWlnPoudr5rv+IyrfffotW6FsaHO13nADvs71f6vT5DNqt0N\/a1OmjtPNo92VbvJ8vWwREYD4E+Prl691lPju7jJ38iGD7R6GP4Pv+udpL2uvUDMmG99lflz7G+n7Z8yOge5rHnfA1M4+dLH8XBzynF3l4f+6svqW8yBNo0wcT8CfBYwpxztiLhX7GH1NPuSIgAiJwLQRGG+7YGyffVF1aSO7zsY1Rp58jhTql1WlT6OsL\/a143H20W502xX0cabfS+qhT2viuOue5tHPoa23X+37aLp5zzJH3eMr6m\/bKdcfiY7HWT72VsVq7+jfVamPU25q0W2ljQzpz6efoQntIPO5jm9P6Wt1z3MfRfT7S14r7fWTMdR+HfB7jyHgr9LXSxqi3sV30\/hzaLu18+lrb9TG\/xzkyx4W2C32ut2Pr7+u0KW1+qzPm0vqp0+8jdQrtVuhrpY1R3xRjfEjG5rifc1qddl88zrEfo01\/K\/S5uJ92q7d238\/YmHiuj\/089\/vYxt3nYxsb0j2P41B8k8\/njDbcocmcxDdwF9rM4+g+H+ljrBX6KG2O6z4y3p\/jMR89x23mtzpt5riPI236W6GPwjiF+pi081xnLue50PYYfa1NP236qVPcpo9Cm\/5jC2tTWN\/l2GsM1eOaQ376GONeqLvQpp+267Rd6GNsH+Fcr8ORttehTl8r9DHOsfVTp4+xTcIc5rrQ7ufT53Ef6fM899FuddrMcx9H2vRTqNPXCn2M7Suc39ajTp\/Xo05fK\/R5fNvIXM71PLfpo9D2WN+mn3H6qY+J5zCPQttz+zb9jNNPvRX3M0ah3cap08eYC236W6HP4xxpe5w6fa3QNxZvY57TH5kzVs\/9nNPqtFvZVIN5\/Thr0ccYhTalrzOH\/lboY96YMN7mU6fP86nT1wp9jHNs\/dTpY2xIGGOOC+2hvG2+RzVcLjZUcMzfz2Uexf2t7r525KGGcoZ82+ZxDuu1efRR3Ed9TDzHR9ZirtscadNPfZswj\/ltHm36W98xddZ3OWQd1ujPp03\/sfZ7zFpDe2N9+sf2y\/g+MZ+zaf5jcjy3Hbnvfn3a9DOPOsdjyj41d53Dfbe5fZvnYJx+6vsI57JGO5c2\/a1vH5112nms2ffRpr\/No6+1XWfeUMx9Q3HG6Pca\/ZEx5rR+2vS3vsfqrNHO6dttbJO+77x+Ta\/Dc7ne5rjPxzY2pg\/V4nz6x+aM+Ucb7mOLMb+VsQX7\/nYO9X58X5u1+rJvrbF5\/fq0x3KH\/J99dv+fYz12\/lDNvm+s5r5PmH79U9rceyunWstZDK21KXaM\/bRrUt+1JnP70s7dFGvzHqNvqrkvJ9b0ue1e6O9LG99H79ejvU+dXeawdl92mbdrTr827V3nnjqPe2ll1\/XaOdS3zfPnDXNdts1p4z7HxzbW1z2nHfs5u9ijDXeXyZ7DTfDwrXhs07jvvE01PdbupdU9foyxrdvqXps+npE2R9rUXWgPicfnPHLfPBP3yJE29ccI53Au53CkTZ3iNn0u9J9KfA0fub6v5T4f25jn7DOyjtf0cdc6nt8fOf+Qupw\/JLvU7O+Fc4ZqtT7OaW3X6R+SNu71OTLXY2Mjc4bE8xljLdocaVPfRzh3SPapNTRnqDZ9Q7lT+pwb9+Kyy\/r7zvM1fGSdU6zn9fvjLmu1Oat\/h9s6qHPTLEx9H+H8Y8zjHoZqDfna9fad19bYph+6xqHzt+3P45vWYczz5ji2e9p2521uX+c5+\/Np089c6hyHZFNsKN99+8zbdQ73PZQ75ON+xvyMueyS47kc+\/l9mzmPkfZMrd7WOGSNU9T0vfX3deha2+Zvi\/u+2nFoDvdNf5v3WJ01Ns3ZFh+bu8u8TTk811B8yMc9jPkZe2wtzhmT1b\/D5WKtcIGxCUN+5u8zf5d5\/RyuQ1+7D9r0t8J4a1NnHv3HEtZj3Vboa+vTZpxj66dOH2Ot0MfYrtLObfV2Pmu2Mer0tTn76KyxqRZj2+qO1XA\/a1Bob6u1Kc75rONC2\/Opu99H+hjn6D4f6WNskzDH8znS7ufTx5gL7X4Obfo9hzal9XmMvqGY+xlzoc\/ncaTtsaGRcea50G7zaHvMR\/ranG068zmXea7TdqGPMRfajHF036aRecxvhb52Dm3GObb+VmeMOS602zh1+jzuI32M7SLM9Xk+0udzqbvfR\/o8PjQy7rkcaQ\/lbfJxDue60G7zaXuMI+027jr9jFPoa236aNO\/SZjD3Fbo8znU2xh1+hjnSNuFNv1jwrjn+kjfWP6Yf\/UtZU5sZSiZ8b6\/9VF3YR51ji59u\/UzRqHPR+oUHo6+VuhjrC9tzpDe5jPe2rvq\/Xm0Wxmqw\/iQnz7GWqHPhX7Xh0bGx6Sf38\/rx1ubua3t+pB\/k28o5rXacSyPfhfmU+fo0rfd72M\/TtvFc3x0v4\/u5+g+H+nbRTyfY5vf2tRdmEOdY1\/op6z9uPfb4TbFOKcfdx\/9FLc5jgnzXJhDnaML7VbcPzYytx9rfdRb6efuY7f1qA\/V2MXPHBfWoM6xFfpa6cdamzpzObrQbsX9PrYx6u7fNDLPZSiPsSF\/62OOC\/3UObrQdqGPOse+0E9xP3UX+qhz3CTMaaWf28aot3HaLvRT5+gyZNPn4nljI\/M85vqq4bpTowiIgAgsgQD\/0E3xN7Il7Fl7FIHZN9zHvKD4AnTR1YrAuQg85jl7rj0ufV0ypoydQ34RmCOB2Tdch+aNlOPYC41+F+b53KHx3bt3kIiBngN6Dug5sMznwND7+tx9i2m43kg5DjVT+lvYtIfymMMX2N\/8zd\/gxYsXEjHQc0DPgSt7DlzG+96f\/dmfrT408T19KbKYhntMoGy4\/\/qv\/4pXr17hm2++kVwwg6+++mr11NFdX\/7zXHd9+Xfs79e8a76Hr17cC3qYfcMd+5TaZ7xrXjvv179+hl\/\/+jn+8z+f4\/lzySUyePbs2erKOV7i+XSm9euWd8zL5iguay6XyIJ3zLs+hZyy5uwb7r6HZwPmt5U3zS8FePkS+PrrTVmKLZnA7\/\/+7+Mv\/uIvwHHJ59DetxPgHeuut3NSxvkIzL7hsmmyebZCnyOjnzp91F1o0y+5bgJ8E\/7yyy+vG8KVnF53fSUXveBjzr7hki2bZyv0uXz77beu3vtFAHdOKSIgAiIgAiIwAwKLaLgz4KQtiIAIiIAIiMBBBK664YZwELslTNYeRUAEREAEZkLgqhuu34GZaxpFQAREQARE4DQE1HBPw1VVl0BAexQBERCBCQmo4U4IW0uJgAgcToC\/uGZIvv\/+e1CGYvLN+9c3Hv6sWEYFNdxl3JN2KQJTE5jlemycY7+W9Y\/\/+I\/x53\/+5+CoX9u6rF\/fuMRf07jPC0QNdx9qmiMCInAWAmy4\/JV++lWdl\/NrHJf6axr3eQGo4e5DTXNEQATOSoC\/2u+sv7JQvwr2aL8Ol3d51ifThItfdcMNoSNt1o16FAEREAEREIFTEbjqhnsqqKorAiIgAldMQEcfIaCGOwJGbhEQARE4NwH\/3fBD4z57Y51d5u2at0st5awJqOGuWUgTAREQgVkR6P8O+b49q81qM1sJrBru1iwliIAIiMCMCZgBZoAZYAaYAWaAGWAGmAFmgBlgBpgBZoAZYAaYAWaAGWAGmAFmgBlgBpgBZoAZYAaYAWaAGWAGmAFmgBlgBpgBZoAZYAaYAWaAGWAGmAFmgBlgBpgBZqcHzIa9yyq75u1SSzlrAmq4axbSREAEFkrg6VPgEsTscRfg3\/rlSOFsjq3Q50I\/9XakTqHfxe12pE7xHI60Ka5zlIwTUMMdZ\/MhokEEREAE5kuADY+fSCncJcdWGKe\/L\/QfktfOd72\/huz7BK664YawhmG21qWJgAgsi8Dr18AlSAiP586m2Z\/FBujSj7k9NM9j7TiWR7+vQb2dI32YwFU33GEk8j6GgHJFYA4EUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSAlICUgJSOg5JNkE2QJfjVH1YpV3nYVSeIQJquENU5BMBERABEdhKgE3XZWuyEqCGqyeBCJyMgAqLwPQE+MnWmyDHU+6Aa7mceq1TnmOq2mq4U5HWOiIgAiJwAAE2tv70IR9z6G+FPgp97UjdxWO0XfeRPhf3cWSTdaHtORqHCajhDnORVwRE4MQEVH75BNhkXZZ\/mtOfQA339Iy1ggiIgAiIgAhc989wP\/10\/QwwW+vSREAEROC8BLT6JRLQJ9xLvFWdSQREQAREYHYE1HBndyXakAiIgAh0BPgXkjrt\/uOYv81qc1p9LKf1b9PH6m2bd6z4Uuuo4S715rRvERABEZiAgJrr8SCr4R6PpSqJgAici0ApQClAKUApQClAKUApQClAKUApQCm77bAUoBSgFKAUoBSgFKAUoBSgFKAUoJTj1utV498A7jc82vT3Ujeaj83fWOw2eOx6tyWv4r9pGu5VoNQhRUAEzkbg5UtgV9llk7vWYt4x65ntUu0uh823lbtAT2GOu6i7DPkYa\/3U6aNQp\/R12i6MU2j7SJ1C+5rlqhtuCOurN1vr0kRABERgCQT4SbOVbU2N8Tbfz9j6qDOPMeo+uk7bhXn0t0LfWLyNec41jVfdcGdy0dqGCIjAoQTevgV2lV3W2rUW845ZL4QH1djMvFFxpN0m0efS+h+rew2Oj507lt\/f61jetfjVcK\/lpnVOEbhkAiEAIQAhACEAIQAhACEAIQAhACEAIexGIQQgBCAEIAQgBCAEIAQgBCAEIAQghOPW263aXRabI5uay13gkcqx6jxy2atLV8O9uitf2IG1XREQARG4EAJquBdykTqGCIjA5RLgJ1j\/FNqe0v2MUdrYkD6WP+ZnjTZGu5U2xvUp9LU50tcE1HDXLKSJwNIIaL9XRGCskdHfiiOhb0xnzKXNcR9HhtH7zQAAEABJREFU93OkTaFO6eu0XRin0ObYypCvjV+6ftUNN4RLv16dTwREQAREYC4Errrhtpdg1lrSRUAE5kzg3bt32CrKWQSjOT\/Pjr03NdxjE1U9ERCBkxH45JNP8OzZM3zxxRd48eKF5AIY8C55p7zbkz1xZlJYDXcmF6FtiIAIbCfAN+VXr17hm2++eSD\/8A\/\/gL\/927\/FL3\/5ywexofwr982KEe90++0vP0MNd\/l3qBOIwFURYNN9\/vw5huQP\/uAPBv1DufINMzwHF97pNTyJJ2+4\/GvjDpY6xW2NIiACIiACInB2AifawKQNl83V\/1q467Spn+h8KisCIiACIiACsyAwacOdxYl7mwihc3z3XTfqUQREQAREQAROQeBCGu4p0KimCIiACIiACByPgBru8ViqkgiIgAiIgAiMEpi04frPa\/kzW+rcVavTlhyfgCqKgAiIgAicn8CkDZfHZaOlUKe0Om2JCIiACIiACFwigckbLj\/ROkjqFLc1isD0BLSiCIiACExDYNKGy+bqn2hdp019muNqFREQAREQARE4D4FJG+55jrh51RC6uFk36lEERGBNQJoIiMDxCFx9wz0eSlUSAREQAREQgXECarjjbBQRAREQgQ0EFBKBxxGYtOH6z2v5M1vq3Gqr05aIgAiIgAiIwCUSmLThEiAbLYU6pdVpS0RABERABJZPQCd4SGDyhsst8FOtC22JCIiACIiACFw6gckbLhstP9W60D4n5BC61c26UY8iIAIiIAIicAoC64Z7iuq9mmyubLStmzb9rU+6CIiACIiACFwagUkb7qXB03lEQAREQAREYFcCari7kVKWCIiACIiACBxEYNKGO\/TtY347mf6DTqHJIiACIiACIjBzApM2XLJgc2WTdaFNv2TBBLR1ERABERCBrQQmb7jcEZusC+05iNkcdqE9iIAIiIAIXCqBszTcPkx+2u37prI\/\/XSqlbTOlRLQsUVABERgRWAWDXe1Ez2IgAiIgAiIwAUTUMO94MvV0URg9gS0QRG4IgJquFd02TqqCIiACIjA+Qio4Z6PvVYWAREQgU0EFLswAmq4F3ahOo4IiIAIiMA8CZy04fJvH+8i50QTwnp1s7UuTQREQAREYMYEFri1kzZc\/7e2u4wLZKcti4AIiIAIiMDOBE7acHfehRJFQAREQARE4MIJTNhwL5ykjicCIiACIiACGwhcd8N9+RLx68\/xGi83IFJIBERABERABA4ncN0N1wzBKgIM5\/7S+iIgAiIgApdN4Lob7mXfrU4nAiIgAiIwIwLX3XBDWF2Ff8I1W5l6mB0BbUgEREAElk9g0obLf5O7fGQ6gQiIgAiIgAg8nsCkDZf\/HpdN1+Xx29UMERCBPgHZIiACyyAwacMlEjZdFzVeEpGIgAiIgAhcA4HJG24L9eyNN4TVdvxnuCtDDyIgAhdEQEcRgfkQOGvD9U+4beOdDxrtRAREQAREQASOR2DyhutNlqM32uMd57BKZofN12wREAERWBIB7XVaApM23LbJstn2jzrk6+fIFgEREAEREIElEpi04c6uod7cwN6+xxO8X+Ldac8iIAIiIAInI3D8wpM2XG6fn3L7Qr9EBERABERABC6ZwKQNl42Wn3L7Qv8lQ9bZREAEREAERGDShnti3CovAiIgAiIgArMlcPUNN4TZ3o02JgIiIAIicEEEJm24\/FZy\/9vHtOmfA1OzOeziRHtQWREQAREQgbMSmLThsrnytBxdhmz6JCIgAiIgAiJwSQQmbbj8JLtNLgmuzrIYAtqoCIiACJycwKQN9+Sn0QIiIAIiIAIiMFMCZ2m4\/u1kjmflYgY84a+9eIKEAn2JgAgMEJBLBETgKAQmb7hssu23lWkf5SQqIgIiIAIiIAIzJjBpw2VzZbNtedCmv\/WdS\/\/uu3OtrHVFQAQWSkDbFoGdCUzacHfe1VSJIUy1ktYRAREQARG4cgLX3XCby\/8Utz\/Phb5EQAREQASORkCF7hGYtOEOffuY306m\/96uZIiACIiACIjAhRGYtOGSHZsrm6wLbfolIiACIiACInDJBHoN97RHZZPlCmyyLrTPKiGslg\/Q35hagdCDCIiACIjASQhM2nBPcgIVFQEREAEREIEFEJi04fJTrX\/KnSMbs827UlQEREAEREAE9iUwacPdd5O7zmMzb2XXecoTAREQAREQgVMTmLThshnyQBz7Qv8hwnr8BN0KfVtrvn2Lz+N7vMTrralK2ERAMREQAREQgU0EJm24bTPs65s2qZgIiIAIiIAILJ3ApA33lLDYwE9ZX7VF4BACmisCIiACkzbcsW\/xjvn3uR7WctnWhN+9eweKr\/PDDz9AIgZ6Dug5oOfAPJ8DfL+m+Hv20sZJG+4UcNhkXdh4N635xRdf4MWLF\/i3f\/tfq7T\/+q8f8Jvf\/EZyQQz44vz+++91p7O+0+O85nTXx+E45\/fAv\/\/7v1+9Z\/O9e\/WmvbCHSRouGx+FbDj2hQ2Ssanl1atX+Oabb\/D8+bPV0v\/xHx\/jJz\/5ieSCGPze7\/0efvd3f1d3ekF3OvYa1V1f\/ntXSmn1nv3VV1+t3rOX9jBJw2VDpRAOx77Qf6iwiT+2xrNnz26b7XN88sknd1M\/\/vhjSC6HwY9\/\/GP81m\/9lu70Cp7Xl3rXej9avx+FEFbv2XzvvnvTXpAyScNdEA9tVQREQAREQAROQmDShstPtic5xW1R1uan3Fbouw3pPxEQAREQgYsksKxDTdpwiaZtiK3O2KHCBtvKTvVyxuvyBO\/xZKd0JYmACIiACIjAPgQmbbhssG1DbPV9Nq85IiACIiACIrAUApM2XACz5mI26+1pcyIgAiIgAgsmoIa74MvT1kVABERABJZDYNKGy28h89vKs8ITwt12As78ERf6EgEREAERuFQCkzZcb7Yc+3KpgHUuERABERABESCBSRsuP+GOCTcjEYENBBQSAREQgUUTmLThLpqUNi8CIiACIiACBxA4acPlt4132duuebvUenROCHdT9DPcOxRSlkZA+xUBEZg9gZM23KHTn7W5Dm2o5zPrOWSKgAiIgAiIwBEITN5wj7BnlRABERCBxxBQrgjMgoAabowor\/mLHd+jIkJfIiACIiACInAKAmq4p6CqmiIgAiKwFALa52QE1HAnQ62FREAEREAErpnAyRsu\/5JUK4Td2tTpk4iACIiACIjAzAgcdTsnbbhjv+RiyH\/UUx1QzOyAyZoqAiIgAiIgAiMETtpwR9acnTuE2W1JGxIBERABEbgwApM33FN\/C3mf+wlhPctsrUsTAREQAREQgWMRmLzhHmvjqiMCIiACIiACSyKghnt7WyHcPnz4z+yDcpGDDiUCIiACInAuAmq4JP\/0KfirL17jJS2JCIiACIiACBydgBpuD+l33\/UcMq+GgA4qAiIgAqckoIZ7SrqqLQIiIAIiIAIfCKjhfgDhg5lrGkVABNYEpImACBxKYPKGy196ceimjz4\/hFXJAHVb6EsEREAEROAkBCZvuO2\/w6VOOcnJ9ixqtudETRMBEbhaAjq4COxCYNKGy+bqn3Bdp019l80qRwREQAREQASWSmDShrsESGZL2KX2KAIiIAJLIaB9OgE1XJIIgY\/Qz3ChLxEQAREQgRMRUMMdAGs24JRLBERABERABA4gMNRwDyi3ear\/vJY\/s6XO7FanLREBERABERCBSyQwacMlQDZaCnVKq9M+i7x+jfrmPZ7i7Wp5s9WgBxEQAREQARE4GoHJGy4\/0fruqVPcXsSoTYqACIiACIjAHgQmbbhsrv6J1nXa1PfY+1GnhLAuZ7bWpYmACIiACIjAMQhM2nCPsWHVmDUBbU4EREAERGCEgBruBzAhfFBuB7PbB\/0nAiIgAiIgAkckoIZ7RJgqJQIbCSgoAiJw1QQmbbj+81r+zJY6ybc67XNKCN3q333XjXoUAREQAREQgWMRmLThctNstBTqlFanLREBEbhKAjq0CFw8gckbLonyU60L7blICN1OzLpRjyIgAiIgAiJwLAKTN1w2Wn6qdaF9rMPsXacU4MkTvKlPkHCrQ18iIAIiMCMC2spFEJi04bK5stG25GjT3\/om11O6t6TZPVOGCIiACIiACBxMYNKGe\/BuJyjwR\/jVBKtoCREQAREQgSMRWEwZNVy\/qhBcg9mdKkUEREAEREAEjkJg0oY79O1jfjuZ\/qOc5pAiMa5mR1ToSwREQAREQASOTWDShsvNs7myybrQpv\/s8umnqy0EdB9vzaAvERABERABETgagUkbLpssd84m60J7FhLC3Ta86d45pIiACIiACIjAgQQmbbgH7vW000O4qx9Rz\/hzXOhLBERABETgAglM2nD5qdY\/5c6OZYyz25I2JAIiIAIicDkEJm24s8cWwmqL\/KdBZitVDzMjoO2IgAiIwFIJTNpw\/dMtx77MCeBLvJ7TdrQXERABERCBCyAwacPlt5THZBYs377FE7xfbcVsNehBBBZEQFsVARGYM4FJG+6cQWhvIiACIiACInBKAmq4PbohdA79P3E7DnoUgUshoHOIwLkJTNJw\/ee17WGHfG1cugiIgAiIgAhcEoGTN1w2Vv+5bQvOfYy3\/nPrMXY7KAUw63Q9ioAIiMBlE9DppiBw0obLZsrGuukgjDNvU86UsS+\/XK9W61qXJgIiIAIiIAKHEDhpwz1kY+eaGyMQQrf6r\/R\/6utA6FEEREAErpjAsY6uhjtAMsbOWQpg1ul6FAEREAEREIFDCKjhDtDjt5Xfr\/5F7hPUOpAglwiIgAiIgAg8ksBJG+4uP5\/lz2+Z98h9b04\/MBrt9qOt1\/hFo7tPowiIgAiIgAg8ksBJGy73wmbKpkq9L\/Qz3vef3U4JCGG1jVRf6tvKKxJ6EAEREAEROITAyRsuN8emyubaF\/oZn6Xc3Nxty76+qE+5d+eSIgIiIAIiMB2BSRouj8Pm2hf6ZyspASGstheLPuWuQOhBBERABERgbwKTNdy9d3jOic2n3PD1y3PuRGtPRUDriIAIiMCJCKjhbgKbEhBjl1EKkDNQa2frUQREQAREQAQeQUANdwss+3L9s1x8\/TXw+efA06dbZiksAhdJQIcSARE4gIAa7hZ4IUWU0DRd3Pbc8BZmt4r+EwEREAEREIEdCajh7gAqvsnIN+\/xFG9XUitWH3Jz3mGyUkRABK6DgE4pAlsIqOFuAcRwCEDOwJu3AV\/e3Bp03op\/h\/lW1X8iIAIiIAIisJGAGu5GPPeDIQA5A2\/fAjF2sVqx+rFuZ90+mgG13ir6TwREQARE4AMBDbcE1HBvITz2vxBuP+2+AWLsZtbafYuZf59q9b1mKk+eAC78S1Y5Q18iIAIiIALXS0AN94C7f9M0XTPgpn4+XM0Md3\/D2Ww4R14REAEREIGLJjDacC\/61Ec83OvXt432BogR+Dq+wed4g5d4jYybO4F\/1dp9FOYnXvdpFAEREAERuAoCargHXnMIQM4AP+1SXnFzs+cAABAASURBVL+N+PQm4WvkO+HfbmYDhn+ZubZ9NANevgT47WmO0JcIiIAIiMASCajh7nVr45NCAHIG3r\/v\/nIV\/4IV\/2bzr2KGN96KuPpRb864+zK7U9eKWddsS+l8pXR2Z+lRBERABERgQQTUcE94WSEAIQA5f\/gE\/DaAjZffdjbD6se6\/O4yP7y2I\/\/OFdxZ6\/0dlgIwQZ9273ORJQIiIAIzJ6CGO+EFhdA1Xv7cN4RuYbNu9EczoFa3PowpdR+XY+wctQKldPq2RzZndvOcgVq3ZR8lriIiIAIiIAIPCajhPmRyck9KXf+8uQFiBFIC2IRpU7gB\/szX5Wl9jVwCcniNggSEAP8yA0oBanVPb6wVMMPq4zSbr39y5idk2hQ2ZEpv6gPTDMgZq58pcz4lZyBnoNYH6XKIgAiIgAisCajhrllMruW8\/sSbEpAzkDPAn\/\/iJsP\/4pUZVv3y6xLAvwH9xN6ivnm\/6qOldP2PfZO9lGMpzVFYLMbGcauaAaUAtQK1AmaA2W1gy39syvz1WqUApQCloNvY11h9m5txb8JbSq3CZkCtQCndIXxuzoDZKuVRD2aAGZDz3X4++sd\/7EqYdeMhj2ZAzkDOB1TRVBEQgWsloIY705vPGavGe3MDhPBwk2ys7H2UNlpr17u896162B+9QU63TfrmDXBzWzBGIAQgRiBGICUgpbbMZj0EIAQghPt5ZkApWDXh+5GHVs5Y\/c0xHoSbLAUopZvLQ\/EAjD2c+dDD3FY4v1bADB\/99V\/j6Wef4aOf\/vThvDFPzkDOuPdJnnvkGqxNoT423\/3cf85Azvdr5QzkvPbl7DO2j7dn4rlQK1BKV4N\/0uJ+ct4+3zNKAUoBSgFy7urwjDkDOQO1euZuo8+tdTyfex+PriM8CyXnte8YGtevFagVKAUoZX1u7v8xa5SynpszUCtQK8A1XB5Tb1Mu622Kb4txfq1AKUDOQCnbZnRxn9dZ3WPfR7uLHP54zFqH7+YkFdRwT4L1eEVz7r79zA+qlDe3PdOrl9JpIWD1z5LYSzsPwOduKUApWPU\/fjr+\/OuIp7\/IyH\/0BuXmLV6G5t8Nf\/p6NcfnD47cAL\/3\/fZttymO9FHoT2lw2qDzV7+67w4BCOG+r9b79pjFw1LaeAit1elsgJ02\/sjmxYZKKQUoBSgFKOX+HJ79vuehVesH+F8DpQClAKU89HGth7Pve3i+nAE2IgrP8vIlUEqXxzjr0Nd5Nj8yz4XzSgFKWe+N9cmi1s11PFpKN9fn+R5pU2ctjp4\/Mq7cPAuF+\/J5nEu9Fe5\/NWHLA\/N8Hxwp9FFKAUoBStlS5EOYc7mXdi73ST+FMZcPUzYOrEPhXBfaOQM5Azlj9d2jjUU+BJ0N12ct1uHY2twr\/R+mbBx8Huu67qP7aFPfWOhDkHnMp3Bf3IcLbfo958OUTcM\/\/+Y3+OSHHzalzC6mhju7K9m8oRi75upZIXR2jEDOXR+8uf0QGwIQAhCCZ3ajv4\/xeV4KUCtQClZNmc93Ss5AKd0f4P11QD\/1bOkuxho5f6gbE+zmNeqb9wAbcOdeNXGzD0Z\/SKn74TWblwvnUr+56bL5Auy08UfmpgSk1NXzGrfjDz\/7Gf7\/n\/4pfvi7vwOYN16li6TUjXwMAYgRCAEIAUgJqz\/Z8A8XjG+TGNcZIQAhACE89K0941oIWF1SPyMEICUghC5SSjduehxiGgIQwsNZvPSH3u0eM6BWoFbAbJ1vttbHtJTuR8wAs\/s+WmZ83C6lALUO54UAhACEMBzve2sFzDpvCN14yGMpQClArUCtQK1AKd1dszlSzIChO8PIlxlQK1AKUOtI0iPdZkCt9yeZ3bd3scwAM6BWoBSgFKAUoNb1bL7+19aoxmb77qOPRuNzDKjhzvFWtuwpxq6x8n3\/ze0n3hDWE0IAcu7ifN5SbvvO6n+4wH4Twv3cGNc2NTOs3tfZTEsBagXMADOg1vuxUjqb7wVsyBS+P3PMGas\/mFOnsB7rU8yA8uUbrH4e\/TLh6ecBOQM5A2a3GSEAOXeH4AFvXfzP7Hbe7Zo502ok567REkhKTQD44ec\/x\/959QpICYjxXmzQYI0WGtenTbmNWYhY\/cU17PDFue8\/\/AGE8136Pto7lFulhNCdlbW93u2+7i44xlXaxgef5yPX7+uszydMShtL3QXbGpwXIxAjECOQEpBSt+9a76aMKjyP12MtSkpASt0fmmhTvvxytMSDQAhASgDP5eJr+NkfTBpwxAjE2O2D87yG1+TeXQamP3CFAIQAhADECMQIhHA\/LYRuPdx3P7BS6vLIJkYgBCClzsf9cb8uDyYPOHg2zuvXcx\/9FNofppt9UPpDKd0+mB8jECMQQpcVApBS9\/xgvPNuffzfH3+8NWduCWq4c7uRHfcTApASEMJuE0IAcu56GF9zfI34yNcVdT7XQ1jXCwGIEUgJSGnt36aZYdW0a11nloLVj0RzBtiU2YBL6eJmXT7\/MM\/mfCe3jfjpy3iXTz\/neV7OXQMuZXjMGfiXfxn\/E7BZt5ecu33cPYawavw5AzkDtQJmQM5YfUeXe+BecgZKAWq9m3lPMQNqBUoBcl6f\/17SYwy\/qJSAGIEQ7s\/OuWsot14zwAwoBcgZqPXW6f+FAIQAhACE4N71GAIQI5Bz9ya4jjzQzICcgZwBQwBCAHKGvX7T\/aGKP7aw18i3P7KoIQEpYdOXGVAKcFsCCAGgQnn9Gqu\/yk\/dJSWYAaUAdLV1zQAzgP7P43vkL99282MEYgRibNMHdTOgFKDWJvzmTcf4trDZB38IQIxAjEBKQEpASh+Cmwd78xZ38vpDbb4Y\/a7f3Ppo3643VskMKOU2+vqWEfMoPs99MQIhACEAIdwm7\/hfjEC\/XowAfS4x3rHm64LC0L0VUgLopHBvFJ7Lz3m7T4sJNebVWUrB6JdZF\/riJz\/plAU9quEu6LKOtdUQgBjvVwsByHndkP11wNfF7Wth9V7X+qhT+Jpho06pez9jPm1WDwFICYiRFlAKVo3YrLNDAFICYuxsfzQDzAAzwAyoFSjFo91ohlUtNr8xYWP+kz\/5CJ999hQ\/\/elHq2aZM1ArkDNWdq1Y1eEfAsyAUrBq8HzT4HwKY253qwNm3Tyu7XGOOWP1B4snT7r69DGHdUrp5rBWzkDOXS5zmE9pY6UAtQK1ArX6yoAZUEq3T87lHArXofR9XNt9HLkOhXP6wnjOQM5AKUCtD9fNeb1vzmd9CvWcuxh17qUUoFbcMe6vy\/VcOIfCeUP1GGMu4xTq9FFn\/o9\/\/DFevPgf4J3TT6G\/1vX669Pc18yAWoFagVqBnLv7Y+12nVIAMyDnLu7nYR6FuZScgVKAnDsejOUM5AyUApTS+bnHVliPQl+1AMS42qgZUCuQM1bPz5wBMyDnbh+szzk5r9JXD2ZAKV0+90TxnFpXKaMPZkDOXe2cATPADDADagVqBWoFagVK6dYgaxY0w+q+uRb3lTNgxkgnZkDO3Rzfk5+ZNudQOD9noFbADCilm0P\/H\/7h89tv6nzbFfzwuIRBDXcJtzTxHkMYXzAEIMZ1PAQg567ZpgTECOQMeDP2Zh3j\/Tk3N11zZ5xN2vPpTwlICUgJiBGIEQgBSAlgLiWEdb1dNDPADKs3Ar6o\/c3B59aKVQPmC71W9z4cU1qftY2aAbV29UtpIw91sy6PeygFqHWdY7aOcS\/cq4u\/KfENh7FagVoBM8AMKAUoBah1Xa\/VzIBa1x4zwAwwA8wAM6DWh+tzXYqvO7Rvr+oxt0MAYnSrG80AM8AMqBWoFagVMOvi\/miG1X2VApgBZkCtQClAKUCtnrke3737aNDPjFoBnoPCs1DIliOFugvPgebLDCgFqz9MMbeNmwGlAKUAtQK1dvvmHTGvFKCUzkebfkopzQI91Qyrxsq1uCcfOb\/WrhZ9tH2q2drPGIXr1ArUCtQKmHU5rOkcmNPa1DmXtc26fNoujLfC+WbdLkIAQuh0M6CU+\/OZyzqsXStQK2DW5fcfzbq5XItzOLfWftaybDXcZd3XIncbQtco2UzZYPmpOOeHRwkByLlraMyjsLlSOId2jECMXbOmb0y8gf\/TP\/0nvvrq\/+JnP\/thNa9dNYRuXzG2XiAErP5+la\/LkXvnyD2k1O3R16CPca8SApBSl8M53CNzKW0e80MAYgRS6taMkd5xMVvHQgBiBFICYrzvT+n++lw3RiBGIKVuLfpSAlICUgJiBGJc19mkhQDE2NXh+XjOELoZIQApdWw9xpHruaQExAjECMQIxAik1O3Zc7tqQAhASkCMQIxACEAIQIxASut1eMeffPLDKpZSV4vM+\/VY1wwwA2oFzOi5LyEAKXW1fc\/3M4CUuvOHAIQAhADE2M8CQnjooycEIKVun3wOUfprmQG1MruTELqxfQyhq+E+M8DMLSAEIEYgRiDGtZ+aGVAKUCstwAyotdP5GAIft0sIHQuypvAcMd6fZwaUsvaFAMQIpASk1J2BDPhc4vwQ1rmuhdCt81d\/9Q6\/\/dv\/092LGc\/XcBeDSBs9FoGcgZSOVQ0IAQgBCAEIAQgBCAEIoVsjBCBG3Dbc\/4ecuzdPfzPwF3aMnZ8v8Bi7Fz1zcgZiBEIAYgRyBmLs6raPIQApATkD\/ubO+ayfEhAjEMJ6Rs5dHt9UPJ8683Pu9uJ+1mHM5eame7NJqcvzOOcyh\/Mo9NOXEhAjEAKQczeHeYzlDOTcnZc2hTGK12Ad2r6u53icsZyBEIAYuz8EeQ7HGNfnDgHIGcgZyLlbl\/Nb4ZyUgBCAnLt6jHMfjFGn0KZQpz9GIATg5z\/\/Af\/8z7\/Bv\/\/7D6sfgaTUrR8CkHNXz8+SEpASECOQ0vh+YgRyBnJez0+pY8m1c+783A+Fe2r5uM6RwhwK8ziyRkpASkBKQM5Azl1N7jUEIEYgpW6PnON1GKfQlxJA\/83N+jmS0noO13NhPtdlboxACECMQEpASkCMQIzdGZlLYW5KXT2v4yPjlJw73nzMuZvPPTHG+fRTQuhi9LMG90JJCUgJiBHIuWPQ5jCXds7AX\/7lO\/zO76jhkqdEBGZLIAQgZyAlIIT1NnPu3gRSWvv20ULYbVaMm\/NCAEIAYgRiBGIEcgZy7t70Ytw8\/9BoCEAIQIxAzkDOQEpAjJsrpwSEsDln12gIQIy7Zm\/PCwHIGci5Y8g3eb6Jc0wJiBGIEYgRCOFhvRCAnLu5MT6Mt54QgBhbT6eHAIQAxNjZY48hADl3Tafdo+eHAOQM5OyebswZyLnbo5+ri6wfQwBSAnLunvNsYr4G51CnxNjNCQHIuauZEhAjECMQIxAjEEKXN\/YYApDz+ixcL8ax7Pv+EIAYgZSAGO\/HlmjpE+48b027EgEREIGLIhACEONFHenRh1HDfTQyTRABERABERCBxxNQw308M804NwGtLwIiIAILJKCGu8BL05ZFQAREQASWR0ANd3l3ph2LwCYCiomACMyUgBruTC9G2xKo9+AiAAAF+ElEQVQBERABEbgsAmq4l3WfOo0IiMAmAoqJwBkJqOGeEb6WFgEREAERuB4CarjXc9c6qQiIgAhsIqDYiQmo4Z4YsMqLgAiIgAiIAAmo4ZKCRAREQAREQAQ2EThCTA33CBBVQgREQAREQAS2EVDD3UZIcREQAREQARE4AoHFNNzPPvsMLkPn9th67PKHcuUTAREQAREQgakJLKLhsol+++23cKE9BMrj7TiUJ9\/1EPj+++\/xi1\/8Ahyv59TXeVLese76Ou9+KaeefcNlc2UDXQrQpezzWvbJN+Ff\/vKX13Lcqz6n7vqqr38Rh599w30MRTZnl13mvXv3DpLLZuDPA93zZd8z71d3ffl3zHum+F0vbZx9w+1\/umVD7fscOv0uzHN\/f\/zkk0\/w7NkzfPHFF3jx4oXkghnwjnn\/HHe\/az0nlsiKd6y7vo7nLu+a7+F8L+edL0Vm33BbkGyibKitz\/W+nzbzPd6OvKRXr17hm2++kYiBngN6Dug5sMDnAN\/D2\/f1JeiLabhsnmyix4LKpvv8+XNIxEDPgcc9B8RLvObwHOB7+LH6wVR1FtFwd2m2zJkKmtYRAREQAREQgccSmH3DZSPd55PtvvMeC1D5IiACIrAmIE0ExgnMvuFy62yefaGfQj9HNmXqLrTpl4iACIiACIjAHAjMvuGycQ6Jw2Os1WlT3KdRBERABERgHgSufRezb7jXfkE6vwiIgAiIwGUQuMqG69925ngZ16hT8C6HpCXTxlu\/9GUR4D2O7Zgxl36O+zn2Y7LnSWDorugbkvYEbbz1n1vf3HDPvbsTrM+L4LecXWifYBmVPAMBv9N29G3wnls\/bY9pXA6BTffG2Ngdb4ot5\/TXtVPe2diJ23t23XM5z30caXvs3ONVNVyC5wW00GnT3\/qkXxYB3i\/vuT0Vbfpbn\/R5E9h0X4zxTtsT0KafQn0o1vqkz4cA72yf3XDenO\/6qhruPhe4YY5CMyPAF5vLzLam7RyBAN9IKUcopRIzJ8B7poxt01\/nHMdy5uhXw53jrWhPexHgC9RlaS\/EvQ6sSSJwpQT8dc5xSa91NdwrfcKe\/NgTL8AXXrsk7SW9ENu9SxcBERgnwNd2G6W9lNe6Gm57c9JFQAREQARE4EQE1HBPBFZlpyWwlD\/hfqCiQQREYE8CS36tX1XDHfrWAy+P\/j3vXtNmSqC9V94v7XartOlvfdKXS4B3yTttT0Cbfgr1oVjrk75MArxb3jF3z5E2dRfa9Lt9zvGqGi5BEzwvwIU2\/ZJlE+A9+p1ypN2eiDb9LrTbuPQZEdhzK7xTv1+OtL0UdfpcaHtM47II8O78HjnSbk9Am34X2m38nPrVNVzC5gW40JZcBgG\/U45DJ6LfZSgu3zII8A7HdsqYSz\/H\/Rz7MdnzJDB2V\/S7DO3cYxyH4ufyXWXDPRdsrSsCIiACBxLQ9AUTUMNd8OVp6yIgAiIgAsshoIa7nLvSTkVABERABDYRmHlMDXfmF6TtiYAIiIAIXAYBNdzLuEedQgREQAREYOYEztxwZ05H2xOBiQn4P2UYGqfaCteeai2tIwLXREAN95puW2ddBAH+U4YhWcTmtUkREIFRAmq4o2jOH9AOREAEREAELoeAGu7l3KVOIgIiIAIiMGMCargzvhxtbROB6435z1g5uvRpuJ9jP0ab\/lbo68u2eD9ftgiIwGYCarib+SgqApMTaBud6\/1N0N\/+nJe251AfizGnH2cufYy50KbfhbbHNIqACOxHQA13P26aJQInI+BNrh37izHW99FmY2SMugtt+mlzpE29lb6vb7e50kVABPYjoIa7HzfNEgEREAEREIFHEVDDfRQuJYuACCyfgE4gAuchoIZ7Hu5aVQQOIsBvDQ8V4LeC+zHa9DOfI23qrQz52rh0ERCBwwmo4R7OUBVE4KgE2PyGpF3EG6fn0fY4dfdzpO0xjrTpb4U+xiQiIAKnI6CGezq2qiwCjybAxjcm\/WJt3mNizG3nUqfPpW\/TP+SjXyICIrA7ATXc3VkpUwREQARE4KoJHHb4\/wYAAP\/\/VgB6fAAAAAZJREFUAwDOR6PVJragbwAAAABJRU5ErkJggg==","height":229,"width":381}}
%---
%[output:7a5fd55e]
%   data: {"dataType":"text","outputData":{"text":"[23:33:13][INFO]  Training complete. Final train\/val: 2.5966 \/ 2.6132\n","truncated":false}}
%---
%[output:2d058d5c]
%   data: {"dataType":"text","outputData":{"text":"[23:33:14][INFO]  300 SMILES being generated (T=0.80) ...\n","truncated":false}}
%---
%[output:3bd50a65]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 30\/300) Generating (pre-RL)\r[##--------]  20% ( 60\/300) Generating (pre-RL)\r[###-------]  30% ( 90\/300) Generating (pre-RL)\r[####------]  40% (120\/300) Generating (pre-RL)\r[#####-----]  50% (150\/300) Generating (pre-RL)\r[######----]  60% (180\/300) Generating (pre-RL)\r[#######---]  70% (210\/300) Generating (pre-RL)\r[########--]  80% (240\/300) Generating (pre-RL)\r[#########-]  90% (270\/300) Generating (pre-RL)\r[##########] 100% (300\/300) Generating (pre-RL)\n","truncated":false}}
%---
%[output:74269a7c]
%   data: {"dataType":"text","outputData":{"text":"[23:33:56][INFO]  Sample of generated SMILES (first 8):\n","truncated":false}}
%---
%[output:42d269fd]
%   data: {"dataType":"text","outputData":{"text":"[23:33:56][INFO]    [  1] CCCOCCOC((CCCCCC1CCCO1CHC=)C=C)CCCCC((CCCC1CCCCCC1(NCCC4CC)CCCCN1O)CCC)C2C=O=Cl(COOCC=H1)HCCC1C=1241\n[23:33:56][INFO]    [  2] C2CCCCCCCCC2C(CC(HC[CCC)C)@COCCCNCNCHCCCCCCCNCCCCCOCCCCCCNCNCCCCCCCOCNCHCCCCCNOCCCC@HSCOOCCO]CC(C)[\n[23:33:56][INFO]    [  3] CCC(1CCCC12)2C=O(H)[CCNCCCCCOCCCOCNOOCCOOCCNCC]C(C3)C(CCC(CCCC2=CNCO2CCC)CNCCCNC2[C)CCCCC-CCCOC]23[\n[23:33:56][INFO]    [  4] Br=H2CN(CC(C1213OC)NCOOC12ONC)COCC3C2[CCCCCCCNCCNC@COCCCCCC@CCCNCCCCOCCCC]CCCCCCCCCC=CCOC(NCOOCCC)12\n[23:33:56][INFO]    [  5] CC(=(CO=)CC)CCNCCC=(C)1CCCOCCNC=C1CCCCC=C(CCCOCNOCC(2CC1(CCC(CCCCCCCCCO(C5CCCC(1CCCC=CC1C))))))125=\n[23:33:56][INFO]    [  6] CCNO(C3C3CCC((ONCO(C))CCCCCCCCCC=C(C)CC=CCCCC)1HCCC1)(C(CC(=CC(3)CCC)C22(CCC(CC(CCCC(C1(C)))))))13C\n[23:33:56][INFO]    [  7] O((NOCSCC(CCCCCC1CC2OOCNC(C(C1C)(CC(C)CC(O((C)CC((C[NCCC)CCCCNCCCCC)]CN1=CCl))COC=C1C1CCCOCC)))))12(\n[23:33:56][INFO]    [  8] C(C(CCC)CCC2C(CC(COCCC)((=CCCCC[C)C)CC]C=CCNC(2)C(C32CCC))(C(H(CS(CCCCNNCN(C(C(C(CN)1)31C())))))))2\n","truncated":false}}
%---
%[output:4c71865b]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 30\/300) Validating (Generated)\r[##--------]  20% ( 60\/300) Validating (Generated)\r[###-------]  30% ( 90\/300) Validating (Generated)\r[####------]  40% (120\/300) Validating (Generated)\r[#####-----]  50% (150\/300) Validating (Generated)\r[######----]  60% (180\/300) Validating (Generated)\r[#######---]  70% (210\/300) Validating (Generated)\r[########--]  80% (240\/300) Validating (Generated)\r[#########-]  90% (270\/300) Validating (Generated)\r[##########] 100% (300\/300) Validating (Generated)\n","truncated":false}}
%---
%[output:58037d9b]
%   data: {"dataType":"text","outputData":{"text":"[23:34:05][INFO]  Pre-RL Generation -- Total 300\n","truncated":false}}
%---
%[output:68c1252a]
%   data: {"dataType":"text","outputData":{"text":"[23:34:05][INFO]    Validity:    16\/300  (5.3%)\n","truncated":false}}
%---
%[output:1307a96b]
%   data: {"dataType":"text","outputData":{"text":"[23:34:05][INFO]    Uniqueness:  11\/16  (68.8%)\n","truncated":false}}
%---
%[output:54f8ce50]
%   data: {"dataType":"text","outputData":{"text":"[23:34:05][INFO]    Novelty:     16\/16  (100.0%)\n","truncated":false}}
%---
%[output:9f971cbd]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAdwAAAEfCAYAAADr33fvAAAQAElEQVR4AeydX6xdx3Xeh7xXIimbDCuZpq5qGYRqN1LcwIqF1i6aUqnABzdyW8NBXJBqrIfWSR+KqIpBGwIKtGndGoFSC27dAoETl1YRAe5LWyAGCrQpmocCYey0dgBbhutIBqiEDOiCqGz+EcWS4W9fLnKdubP\/nr33mb3PJ+g7a81aa2bWfLNnr3vPPdxn9w39JwbEgBgQA2JADAzOwO6g\/8SAGBADYkAMiIHBGVDBraJYPjEgBsSAGBADPTGggtsTkRpGDIgBMSAGxEAVAyq4VezIV8WAfGJADIgBMdCCARXcFmQpVAyIATEgBsRAVwZUcLsyp35ioIoB+cSAGBADEQMquBEhaooBMSAGxIAYGIIBFdwhWNWYYkAMVDEgnxhYSwZUcNdy27VoMSAGxMByDDz00EPLDTDT3lW8qODOdNO1LDEgBibKwATSpqi88sorE8h0\/BThBX5SM6vgpliRTQyIATEgBpIMUEwoKkmnjAUD8ANPRcO9NCq4dIzhxpDagQH47NBtkC7k0hR9JsCcXcdbpm\/XOcv6kUsdyvrKfoeBKg7vRO3U6LfTOpwlni9uDzfz9shjz7c9a\/0redWhfpTaiKUD6nLEH0+SssUxTdqVBZdJANU6BnbQZBLF5M1AvLe0yRgZA7uwk4GYp7its7KTs5Ql5s3a8AdSfWQbjwH2gD0pmxFfFehf1ndMe1WO+PrIMzVOacFlQjqAFBFl9lTsVG1wMNXcp5D3MtfQMn1XwQ35zuV6WsU64A+k5sbedU9T49WNtcx8dWN7f1luY83vc+lLJ\/eydfU1Rx\/jDJVnacFtknSTpJqMoxgxIAbEgBgQA3NnIFlw+QmEYjr3xWt9YqCKAc5BlV++8RngvqR9GZ93zdgPA8mC22ZoDkAcz4GIEcfQJgYJ0D2wxfB+0+MY2viQAN1AG1jbS+wGs9P2Om2D2b00Xyx9DPodf7VWF+v96ClUz9Ddy1zWG92QssU+H2M6kjgkQPfA5oGvrI3Pw8eZ7v3o2JEAHaSubex9wM+DbrCxre2l+bzETxvpgS2G95sex9DGhwToMbwdHRCDLEOdv6xfE3s8Nu0U\/Fjmx+Z12gAbEqAbrI2MYTEmU\/7YZm36eN3aSGA+ZNw2G3YP7DHwmw3dw+xjSua3+dANKVvss5gpyaULbrxYSOFGFQN7HGttfHXxqRj6YLdxYomPGAP+2GY+7PiB2WKdNiDWYrzEjt8Dm49Bx+ZjyvSqWMbAT1\/TacfAR8xQYHw\/J\/PENvNjx18H4qyPSWx1\/fATZ31MYsNnoG0+k95mcX1KGz8e0+yWB\/7YZj7s+GNgtxiT2HwcbfN5id3HeR2fjzWdGK9bm3j0GNiJj+1DtG0u5ouBz+Y0H22v0\/agj\/mR3uf1OI5YbD6mqU5fQDwSoJeBeYiJgb1Nn6r4snHa2pmDPON+ZscH8Mc27AA7fgPtOlhsU8l4zNU0vmlcp4JLMh42GbayJLHjt1iT2PBZOyWrYuiLP+6HDZ+3p2ze30SvGoP58Ns46NisbTJlM1\/fcsi5UutL2dqsyfdv04\/YZfrSfxkwdxVS+0B8bE\/ZLC9i8VvbJHbTTWKzWCRt83mJHb+3oWPDhz4XtF1PGw5SY2NjjCH5Y3zmSc2BHX\/sw4YvtvfRZuwqpOYlPranbJYfsfh9G1sVLNYk\/avAWBbbp+xUcEnGsGwyjLPsGKn+qXFTNvpCPDJHkHOcH23sli86NoPZh5bMG8+RshFDbsg6lPWv64d\/mb70XwbMXYXU2MSn7DnYuuRGn3ifaWMfa03MxZyGZeZlrGX659p3yHUxdhVSnBCfsg9pY84UmBM7cgh0KrhVidiFnpJV\/ep8qfHMVtfX+62Pl20J9n1j3c81ps4aDD6nMXOwufz8ppOb+SX7YKDdGLYPKdlupPyjudYMfr35Z64MPQN+72Ldx\/Wpc90wV59j+rGSBXeZSelbBT95G71qTHxNxoJIYmM06etj4v5x28euQvf5sOYxc2A+P7\/pY+ZQNxc5kacHtrp+U\/azvir0tTbmgNe+xovHYWzmiO1lbWIN9C2Lkz0\/BmzfymR+GddnlCy49d22I9pewG3jt2epfx1q3PqZd0YMkQsXnI2LpO1nxubb0qsZgC849KjuMW8vfAyxQsaF4yHGLhuTOfFNCVPMedX8DskZ1+xQ45cWXJu0bGLsxHjiaWP3NtOx47d2G0k\/+qf6YMef8jW1MUbTWOYqi8eO38ZCx2Ztkymb+apk135VY47lyyl32xdy8hiLi6p5LLdUDLnij33YUzaLRaZi6IMdP3rOIE+Qa67kFvOHrS7fJjHxuL7N+IzhbaZjx2\/tPiTjMW4fY1WNUTUP8+Ov6p+DL5VnacElYRYF6BgDOzExsMextLHHsW3a9GecGNibjkNs3J829tQY2PEbLCa2mx+7xZjEZn6T2MzfVFofk74fNhs7lvh87NA688U50MY+9NxNx7d8yMkDux8jbnvfkDo5MXcM7Kl5sdfFpmLogz01ZpWNPvQ1xLHmR8a+uraNGUvGAnX9iYn7Whtf3B+b+ZGxv2k7HoexsPn+tLF7YPMxXsfnY73P63Gc9cHu4+7o09DI39biJfahV8AczOnnoV0FH1umVxZc68TkMbzPdJNxLG3zeVlmJyblwxaD2BjExDZr44uBDxsyBnaD95nNS+\/3uo9Bx2cSvSmq+uBLoenYcRxjxTbfrvLji0FfbEhDXdvikHWxsZ8+Bu\/jwPi2xaRk33F+jrqx8cfw\/WO9SWwcQzseh3aZHZ+BGIPZlpU2XkpWjU2899NOwcd43ceaHZvpsYx91kZ6xP1oez86NuB12gbsBm8z3aTFeGk+L\/H7tterfHEc58jb0Jv2J9ZQ1wd\/DOvbVTJek74+Dr0Ofkz4Id7b0BsVXAIFMTAXBjgM8VpStjhG7WYMwGXqZtOst6KmwAD7yz5PIde2OS4bDy\/wkxpHBTfFimyzZYCDADgUHtjAbBc+8MJiLgeeTsNnwIDOS3oTqnhRwU1zJuvMGeBQeEx1uawhh9zJw5BDPspBDOTIwPQLbo6sKicxIAbEgBgQAxEDKrgRIWqKATEgBsSAGBiCARXcIVjNZ0xlIgbEgBgQA5kwoIKbyUYoDTEgBsSAGJg3Ayq4895fra6KAfnEgBgQAyMyoII7ItmaSgyIATEgBtaXARXc9d17rVwMVDEgnxgQAz0zoILbM6EaTgyIATEgBsRAigEV3BQrsokBMSAGqhiQTwx0YEAFtwNp6iIGxIAYEANioC0DKrhtGVO8GBADYkAMVDEgXwkDKrglxMgsBsSAGBADYqBPBlRw+2RTY4kBMSAGxIAYKGGgKLglPpnFgBgQA2JADIiBnhhQwe2JSA0jBsSAGBADYqCKARXcKnYKn17EgBgQAyE89NBDoqGGAXFUTZAKbjU\/8ooBMSAGimL7yiuviIkaBuBIRbecJBXccm4G89RdkHV+EmsSQ1xXNB2\/bvzUOE1tfuxUH+9HbxJDnDAdBobYU8ZMoQ0r9K+Kr\/NX9W3iqxrf+8r0qjl8n6o4+dozoILbnrNeepRd1GX2XiZd40HE6xpvfmLp\/CYWo+wawU5sYpjZm1h720XCVZd+beeZYrwK7hR3bcI5cxC7HEj6TG\/ZynhdGVjl9Vo2d5k93qOmcXE\/tesZUMGt52iQCC5qio8fnDZ2s9H2MHssfQx6U39ZbNP+xNkYSIDNg\/Wk7D4Gv8WhG+KYuG1xJs1vbSQ2L03HLkyfAfbTI16R95kexwzRZi4b13SkwfvMhjQ70tpIA3YPsyPN7nWzIWM7bQ9iALZYxjb8Hub3Nuk7GVDB3clJFhYuYIqQB7Y4OWw+Bh2bxaFj88Bm\/jpJrO+Ljs33o40deHsb3Y\/BOABbagzs+A3EoCNjn7XxEQPQhdUzsGwGtrfsqQGbjYtudiR2k+hjIs6FNiAfD2w+L9pl\/iqfHyOlx32ZA5uPxUYbCUyP42ibnxihnAEV3HJuZutpcziaxDaJgUx\/MOlDG3ufYNw+x9NY82SAay\/GkNdO17Gr+lX5uuxa3+N1yWHufVRwV7jDXOAcelJA0kZvC\/p6tO1fF+\/HRq+L7+pnbI+u46ifGKhjgLMWo67PuH7NNkcGVHBnsKvxjYN21bIoalV+7yOW8Ty8v0\/dz2F6anx85GWgnYqTbb0Z4LqwawRJe70ZqV49HFVHbHvh0WKRtLc9eq1jQAW3jqEJ+jkEXdJu0q9JTDw3B5J+SO+jnbJbDD7TvcROX4P3SRcDxkBf1wnXGWPZuF1k1\/5d+3XJsa8+5AxnfY2XGmeqNhXcFe8cF2bqAjU7PgO2OF1s5jeJzeLQzW4SG36k2ZC0sXtgw2eg7f196Yxrc5jEVja+xZgsizM7YzWNtT6SeTBg+xZLv6fmw+azNrtJ71tWtzFjmRqXvHwcbeBt6Nh8f9rYDbS9v6vOODamSWzxeNjM731mR3q79GoGsi24qU2Ol2IxJs1v7Viaf9UyvkjL2tg9fN7YrY3uYXaT3odudiRtg7W9ND0VYz5kHeifivF2dA8fj502e4oew+zEeBBnbXRg7XWQ8ALq1kqMwWKtHUvzjyHZrzIwf+zDBsg59tE2OzFtYH2tD+0yEIMP6YHNYHZrmzS7l+ZDmt3rsc37ynT64PPABrAhDbSBtatkV36rxpyTb5yC25Ix2zQ2GT3VHTt+D2wW6+2mm09yugywl+xzDOzTXdUwmcMRvAD0slnwEWOgbbFm89J8OUvyZR0xsHfNe5m+XefMsZ9xmuIjZctxDavKKbuCy2b6TUPHtiqCNG9+DHBNxMgvy9VmxJmBI8sCHZu1TWLDZ+05SdYVYxXrI4cu83bt12WuNn3IC7Tpo9htBrIruNtp1b\/WbTg3EkPZaK+99lq4du3aqqH512APyq7B3O12hpBlueocLd5D4APo3rLISx98lF2DU7FPtuAawdwIDL4IoxvwW7xJDsTJkyfDpz\/96XDmzJkFvPrqq+Gb3\/xmQMa+KbTJe8r5w\/HU1xDnf+HCBbv0spKcEZ8QZ8Xb0A34fCx61TliH5si5qtpv1XHxXl\/9atfDb\/4dz8Wjh49Gn72Z3820F51jvH8cc6xP9c2eX\/ve98rfkHh2psiJl9w7WaAtBsCut8M2uYzOzeK06dPh2PHjoWtra0FvO1tbwt79uxZsMUxObennj\/c3l5DtDf4poA4\/\/3799ull63kjHBWLEGvY6NNDLqh6hy12aeYrzZ9Vxkb5\/3GG2+Er\/+f74UPveed4Rvf+EY4dOhQdveROOdV8tdmbvJ+88037dKbpJxswY0Pflf2jxw5Evbu3buAffv2FQU3tk+lPfX84Xnqa4jz39zc7HqJjtKP80RB7TpZ6hyxj00R89W036rjUnnD4WMP3odYuK+sOlebP5Wz+XKW5J37OSo2veIlu4LLoefwW87o2KzdRNLHx9FuO4bvL10MZMpAaVpc71z3FoCOzdpelvmwJbT9xwAAEABJREFUN4nzMdLFgBgoZyC7gkuq3Bg47AAdmwEbOnZ0D2wpn9nxCWJgXRjgurfzge7XjZ22l+gGfPSxNpI2dkEMiIFuDGRZcFkKhxuge3gbukccZz5vly4G1omBsjOAHR6QKeAD3kd7UlCyYiAzBrItuJnxpHTEgBgQA2JADCzFgAruUvSpsxgQA2Jgkgwo6RUwoIK7AtI1pRgQA2JADKwfAyq467fnWrEYEANiQAxUMTCQTwV3IGLbDssDBPrCuXPnAuhrvFWMQ\/5gFXN3mbPtfit+GAa67N0yfbhGgY0Rr8rsOUnyBTnlZLnE\/M2trYKbwY5ysfGYSR4H1weeeOKJ8NRTTwVkH+OtYgxyn9Iajh8\/HtjHDC6ntU0B\/vs8R02u+\/g6PXHiRMH\/L\/+X\/11I2k3GGTMmznnMuevmmvs5mknBLa7tyb5wo+Axk88\/\/3x46aWXhIlx8MwzzwT2b7IX4EwS1zma9r1jHc6RCm5GN5v3v\/\/94QMf+IAwMQ7Yt4wuo7VPhf3QOZrefYR9m\/vFq4Kb8Q7zE3sfWPUYGVOs1NaAgVVf\/33NvwZbNfslquBmusUc0rH\/HlX395Wu\/rn\/XSbTS0hp3WRA5+gmCfo\/GwZUcLPZisVEuFHwd8HHP\/6p8ORzL0wW7\/vw0zv+vslzeWMsrr68Rb9ybwh1\/p19ZZkzAzpH6d2tOyd1\/vSostYxoIJbx9CK\/fc\/\/GjYmjIeeXQHg\/Z8Xhxep10H4qti6vxVfeWbLwM6R4t7W3dO6vyLo6nVlAEV3KZMKW5wBuynapNMiA7QgelIA3YDNnSkgbbBbCbNLplmQNbpMcC1TdYmTY\/b3u59ZjeJD9A20PYwu2Q1Ayq41fzIOzIDHGL76dp02uhxKthBykcsPlDlJ04QA3NjgGuea591mU4bHZsHdpDyEYcPVPmJE+oZUMGt4Ii\/\/1S45RqAAQ62DYvOIQdm60MyLuhjLI2xzgzku3Z\/faNzhkCfGTMu6HPMuY+lgluywzz67OTTHw36hG0JQSOYuUFwoEGf0zEu6HNMjSUGcmWAa50zBPrMkXFBn2POfSwV3JId\/pM\/+ZPw+2e+v+MTtiXhMg\/AADcIDjToc3jGBX2P22eOGksM9MWAXet9X++MC8rG7Sv\/OY2jgpv5bv7w++fCDyaMKno5rN4ft\/FhM1jbyzKdPvhArHODAN5OnDBfBnSOXglc74BdjqW3ed3ivM10zhDwMfiEcgYmWXDZZFC+rO1\/j0mMoSo2R9873vGOwKPOvvKZZ8OXP3F8siB\/1sF6cuGZGwTIJR\/lMRwDXHdcf1yHOkf98swZAv2OOu\/R7hTciayTAsomA\/RU2tjxe2BLxeZq40Yxly8zYB258qy85s2AztG893dqq5tUwaVoUkSNZHRs1p6b5GYxh4ews4657Y3WMx0GuP50jqazX3POdFIFt+lGUIibxl67di2k4Ptfu3YtGZPq18Xm5\/I6\/yxpDvBrmrNetfdzXnfua5vDGWINufOs\/OoZmGXBtWXz26+hrAifOnUqnDlzZgFc3OfPn7dhwtmzZxf8cfyybf4J0u3Jbink8OxTHwldvzAgp37r8k+r\/HXC\/rGvdm1cuHDh1s5KjMkA+6BzNCbjmquKgVkXXIqsgcKbIuLYsWNha2trAYcPHw4HDhy4HX7o0KEFfxy\/bPvgwYO35zKFGwX\/LOkff\/Anwq999K\/kjYr8fv4v\/2jyn1axHx627rElObSZsyreXydcQ769f\/\/+NtMoticGdI56IrJmmKpzkeraNj41xhRtsyy4bTbzyJEjYe\/evQvYt29f2LNnz+39jP1DtG9PFimPPfi2MHVESyq+0cd+EDLZZs\/i8XJp++vCriGzbW5u5pLmWuYx9TNE\/vHGcWbs\/JjEFsepnQ8Dkyq4XFT+gkLHlg+dyqQrA34f2Vfgx6Jt8Hb0lB2b+ZAAG0AHpps0m29jA9gAeksoXAyMxoDO0WhUd5poUgWXFXJBceMD6NgM2NCxo3tgwyfkxYDfIzKjzV4BdGwGbMDb0bEBdItF0sbuddrYsaF7iR0bQMcH0LEB2oIYyI0BrlEDuaFzvQJ0bAZswNvRsQF0i0XSxu512tixoXuJHRtAxwfQsQHa64jJFVw2iQ0D6B7ehu7h46Tnw4Dtkc+Igwm8rUonFsQxjG02dGKA2VISP0j5ZBuAAQ3ZCwNc38APxnUMvK1KJxbEMX5cdGJAHOfb+IG3SQ9hkgVXGzdvBjjUhiYrtVhkWTyHHz8oi8GO30BbEANTZcCuY2STNRBnKIvXOSpjppldBbcZT4rqkQEONQfXw4aPfWYvk03jfZwfy+zYTLe8sAFvpz0XxOssW5fFIctiZB+EgcpB\/XXJ3gDrEPvMXiabxvs4P5bZsZlOPgAb8Hba6wgV3Mx3\/ez\/uxT++PXpooxeDl8Mi\/V2b0vp2OriiQGpOLMjQVWM+YibOrgR2nrQy9aDz+KQtMtic7brHG3vDnu4rYXiywxMR+Iz0Aa0kR7YDLHd2uZHmg1J20B73aCCm+mO8zg6Hrr+C\/\/hf4a\/+YX\/OlmQP+tgPZlSvXZpUTS56dnC0bFZ2yQ2fNaeouS64\/rjOtQ5muIOVuQ8QZcKbqabxo2Ch\/6\/9NJLYepgHZnSrLR6YqDqsZZj+VJL0TlKsZK3rep6yTvz+uxUcOs5WlkENws9dH1l9K\/9xPFvt1W\/8aYekWqPtWwieSKUfxRmkz5xDP1Tm6ZzlGIlX5t\/RKrfY66R119\/Pd\/EG2Q2YsFtkM2ah3BBCa+FqXGwDpdtVbFl\/alHpLZ55Gn8KMw2fS3WHpE6tetH+W6fea4j4B+JanuL5Bq55557CJksVHAz2Dp+AufvTCdOnJjFlxXk9MUJY+TCvrF\/7GMGl1PvKdQVWyZMPSLVHmvZRMaPwmzSJ45517veFdgH9mOMfdccR3u9X7Fv7B\/7GO8tba6RqT8iVQWXu8WKwY2av3P29bfaL37xi+Gzn\/1sePHFFyf799+prYH9W\/Fl1Hh63iqmiFoHdGzW9rLK5+Ny0Ps+R03OY3ydPvPMMwUVH3rPOwvJddFknDFj4pzHnLtuLvgqiJvpyygFl0M7U\/5aLavqrSNuFn3+vfa9731v6HO8VYw1pTWwf60uhg7BfZ4jCizjAXSfDjbaXqIb8OUK9mHsa9Vfp\/yGBjePPXgfoviNe+x8msznc24SP1YM+1cQN9OXUQouB9oOK3KmXFYui2J78uTJ0rdgjh8\/XvztsnKQtXVq4TDQ9zliPMDYHmZDpuBjpYsBMdCcgVEKLun4g0vRBdjXBRTc06dPh8c\/\/qnw5HMvLOB9H346+Z2x68KN1tmcgXU\/R82ZUqQYyI+B0QquX7rdNCi6wPvmrt\/\/8KNhK8Yjj8592VrfAAzYOeIMgQGm0JBiQAz0yMBKCi43B+BvGD2uSUOJgbVggDMEdI7WYru1yBkwMFrB5cZgsBvEDPjTEsTAqAzYGULqHDWhfvUx\/DkJrD4TZbBqBkYpuP7mwE0iXnTKFseoLQbWnQGdo+ldARRa+7Dk3\/7QB\/XByOltYa8Zj1JwVVB73TMNtqYM6BxNb+P\/6I\/+qPhAJP8u96vf+k52C1BC4zIwSsHlJ\/PUssrsqVjZxMC6M1B2Xsrs685XTuu3f5ebU07KZXwGRim4XZbFTQTU9SXGYLHWjqX5JcWAGBADYkAMVDPQv3fQgmsFj7RN97LsLTJi8AF0+qeAjxgDbYszm5fmkxQDU2KA6xqQMzIG1zg+QQyIgbwZGLTgciMAUICMgT0GNxPizI6OzdomseGztqQYmCsDXOeA9SFjYBfEgBjIn4FBC64tnxuE6QPKhaEpyIYFR9S4du1aSMGHpfxtbTbejevXw\/UduFG4246p+PTe5chLscFLvqziHC2ZsrqLATHgGBis4PpiZ3pKulxaqfHNh7G9Dd2Ar2zwU6dOhTNnziyAj\/KfP3\/+dpezZ88u+OP4Ju1z584V4126dDlcvHhxAVeuXC58fcxDLuTPfOhTxdTXEOd\/4cKFYo\/bvnDtAvohy4BfEANiIG8GBiu4VuxYvukpiX9ZcBNibBvH69hoE4MeI\/XF2XzR8YEDB26Hln0hMl+K3BT25dj2vY58t6Nhz549xVx9zEM+5L9jrK2tgG8qmPoa4vz3799f7HHbF65dQD9kGfALYkAM5M3AYAU3XrYveOggjunSZhxuQl360ufIkSOBIuhBIbQiSIz3LaMz1sbmRtjc3FzAxsYGrh15dJ3L8u\/aP4d+U19DnD97Xmzyki9c7zYEOrC2pBgQA3kzMErB5aZgRdF02ugxPbGdGGxxHO0yH3b8BtplY1iM5FozMInF++vYdK5r9EksQEmKgTVnYJSC25Zju4lwI0H3\/bHR9hLdgI8+1kbSxi6IATEgBsSAGFgVA1kWXMigSAJ0D7MhU7BY7zNb7pIP2qSQe97Kb+YMaHliQAz0wkC2BbeX1U1skBMnToSjR4\/uwPHjx\/XQ84ntpdIVA2JADMQMjFJw+W2Tt3YBOkl4nbYQwuMf\/1R48rkXFvC+Dz9dPPxc\/IgBzg7nBqDDiNdpC6MzoAnFQGMGRim4PhtuEACbSXQhhPsffjRsxXjkUVEjBnYwwNkBOEyiC2JADOTLwCgFlxsCP5GnkC81ykwM5MWAzlFe+6FsGjCgkAUGRim4CzOqIQbEgBgQA5UMpD48abbKjnJmzYAKbtbbo+TEgBhYNwYorCdPntzx4Un7QKU+RDndKyIquMMshLeSeTtsmNE1qhhYDwZ0jtZjnym4p0+f1ocoZ7jdoxRcK7bIGDPkVEsSA4MwwNlhYGQM7MK8GNCHKOe1n6xmlILLT+ZlIImpQHmKgVUyUHaGsK8yL80tBsRAMwZGKbjNUlGUGBADYkAMNGGAt51TaNJXMatjYLSC698CY7m0kcJcGNA6xmCAc2NgPnSkMC0GUsXSbE1WoqfSNWEpv5hRCi43Bd72AkYBOnZrS4oBMVDNAOeFcwMsEh27tSXzZ4DCWvUpZIpp3Sr0VLo6hvL0j1Jw81y6shID4zGwqpkoxqDJ\/HEc7RSajKWYcgYouGWfQubRrjzOtbz3tkcfqNrmYWqvKrhT2zHlKwYaMkCx5DdggF7WDR9I+ekbIxUnW3sGkkWTR7vqca7tyZxIDxXciWyU0hQDbRiggFIorQ86Nmt7iQ9427i6ZhMD68HAKAWXw8xhB9CKBNhpC2JADNQzwHnh3ACikQA77SHA+Iaq8a9duxaEa1UUFT44KpSbLzeuXw\/Xk7hx0xtC2l\/vY465oiBmwi+DF9yywzrkTWLC+6HUxUCSgVWdI86pgRySyd00njp1Kpw5c6Yz+LvmuXPnOvdfZu5l+sZ5s4abdISrV68iwvnz5wtp7bNnzwaLuXTpcrh48eIOXLlyuehjfh\/TxMccVWuKc66KzclH3q+\/\/nrBzVRfBi24HFA7rLHEN1XSlLcYGJMBzkp8fqyNb6hcmMOPTXwjV2cAABAASURBVLtsvmPHjoWtra3OOHz4cDh06FDn\/svMvUzfOO+DBw8WlG1ubBbyrW99ayGtzRotZu++vWHfvn07sGfPnqLP3r07\/U18zFG1pjjnqticfOR9zz33FNxM9WWwgsvB5ICWEYOPmDJ\/lZ1+oCoGHzEG2oIYmBoDXL+clbK88RFT5h\/LfuTIkUCB6AoKD8Wka\/9V9UvlDee7N7Zvrfh92\/LEtrGxETY3N3cAe+Hf3Olv4rM5yiQ5TZVr+IKbO5iWtn1VTChnbi7cZAB6Wer4iDHQLouVXQzMjQGue3\/No2Nrs076+Hjabcfw\/aWLgXVnYFIFNz7wHH5s8SZiwxfb1RYD68QAZ4CzAND92rH5dkqnD3EG2qk42cSAGGjGwNgFt1lWI0aVfZrPp1AW08Zu47X95KH1azOXYvP7xKzt49iSIgniedvYiAXxGGqLATHQjoFZFtz45sBP6LHNaEp9upJPw9mnC4mr+9Rfk0\/yVX0yseqTh5cuXSKF0CYH8me+JnnlGjP1NcT5X7hwodhHvYgBMbC+DAxacCl0VRiDduYvK7bMn\/p0JZ+GO3DgAO4CdZ\/6a\/JJvtufTKz65GHKd\/f2Jxbb5ED+beKb5D92zNTXEOe\/f\/\/+4lrq8sI1XIUuY6qPGBAD4zMwWMGlyDXBkEvmJkUOVXOkPl1pn+KzfmWf9mtrZ7yNtp88vBlPvzZzWf5t+uQWO\/U1xPl3\/XQl128TcI0IYkAM5M3AYAV3iGVz46GI2tjo2KztZZXPx0mfDANKVAyIATEwaQYmVXBhmgJLMQXo2AzY0L1EN+ATxIAYEANiQAysgoHJFVxIotACdA+zIVPwsdLFwKwY0GLEgBjInoFJFtzsWVWCYkAMiAExIAYiBkYpuLylG807qebJkycD\/8xjUkkr2dkxMPVztMIN0dRiIAsGRim4vL3LzcKQxcpbJHH69OkW0QoVA8MwMPVzNAwrGlUMTIeBUQoudHCzMEy18LIOQQyskgE7Q0ido1XuxIzm1lJGY2C0gutXxM0C6IbhWZEuBtoxwBkCOkfteFO0GFgVAyspuHaD4GYBaK+KAM0rBqbKAOcGcIYA+lTXorzFQKYM9JrWaAWXm4GBmwPodSUaTAysAQN2hpCcIbAGy9YSxcAsGBil4PqbQ+oGkbLNgl0tQgz0yIDOUY9kNhyKf50AmoQ3jbOxfvj9c6ZKrgkDoxTcMQvqXPeNw5zCXNerde1kQOdoJydDWjhvn\/j7fy8cPXo0HD9+vPKfBhL77FMfCSdOnGiU0kNv3Qhf+cyzBX6gwtuIszkEjVJw+ck8RVaZPRW77jYOMgc\/Rt2NYN15m9P6y85LmX1Oa1\/FWiiiX\/3Wd8KH3vPOUPdPA4n9\/TPfL2Kb5PrKD\/9\/ePxdD4Sz3\/56k3DFzISBUQruTLha6TIe\/\/inwpPPvbCA93346dobwWLSaokBMdCWgccevK9xlzaxP7Z1b+NxFTgPBgYtuPzkDaAKGUNvkcFMM9z\/8KNhK8YjjzbrrKhJM2DnhkWY7qXOEcwIYiB\/BgYtuNwIADQgY2AXxEAuDOSah50b8jPdS+yCGBAD+TMwaMG15XNzMF1SDIiBbgzoHHXjTb3EQC4MDFZw7S0vFmp6SuIXxIAYSDNgZwav6SmJf1hodDEgBpZlYLCCy0\/jgASRZcAviAExkGbAzg1e01MSvyAGxEDeDAxWcPNetrITA2JADPTHgEYSA00YGKzgpt72StmaJKkYMbCuDKTOTMq2rvxo3WJgSgwMVnBTb3ulbGVk2U2lzO\/txMZtbDF8jHQxMAUGUmcmZZvCWpRjqHxa1Xz50cqMgcEKrk3QRVIo7aaCXjYGPpDyW38vU3GyiQExIAaGYuDkyZPF0C+\/\/HIheWIcKBp6WTsGBiu4FEIAo8gy4PcgjiJpNnRs1vYSH\/A26WJgTgxw7QPWhCwDfiE\/Bt79kx8sknrLvYcKyRPjeEJc0dDL2jGQKri9kEAhBAyGLAP+IeBvTFXjX7t2LaQQ90nFtLHZeDeuXw\/Xd+BG4e7qa5OHYtP7PTQvxQZ3eLFzQ1fTUxK\/kB8DDzzy3iKpvW\/9kUIWT4zTE+IKLtbxZbCCG5PpC6DpcUyfbX9TYr6ysU+dOhXOnDmzAB5Efv78+YUuZ8+eXYiJ+9S1z53b\/iquS5cuh4sXLy7gypXLxVxdfXFu5M98dTnl7J\/6GuL8L1y4UOzxsi9cyzGWHVP9xYAYGIeBUQouNwhfAE3HPsQyGd+PS7tsrmPHjoWtra0FHD58OBw4cMAPEQ4dOrQdE8XGfcvaBw8eLMbbu3dv2Ldv3wL27NmzlC\/OjfxjW1leudqnvoY4\/\/379xd7vMwL1zDXcgzsy4yrvmJADIzDwCgFd5yldJvlyJEjgSLoQUG0Imijen9XnbE2NjfC5ubmAjY2NnCFjY6+OB\/LP7ZPqT31NcT5s+fFJmf0QqEGTVJqGtdkLMWIgXVlILuCy0\/v\/nCjY2uzQfTx8bTbjuH7S2\/MgAInwoCdCc4Felna+ECZX3YxIAaaMzBKwU0dag4x9lSq2PEDdB+DzbdTOn2IM9BOxckmBqbEANcx17TPmTZ2b6vT4z70x5bqhw+kfN429IfOVjG+Xx96VQ74wY0b2x+ARAfXr2+3zb79wciUzX+Y8pZ\/x4crian3VeU5dR+cThmDFVwOsAckpdrYU+CQg9jXxkYsiMdQWwyshIEOk\/ozg84QSIO1katE6sOHbT6UF3\/IrE3foWL54CGcXr16FRHiDycyr+VtsW+++WYRa33eeOONon39ZvFE4YOR9iFJb\/MfpDQ\/sd6O3sSXypNcDZaztaciyfv111+HxslisIJLoWuCyTKnxMXACAw0OUPEjJBK5RSpDx+2+UBe\/CGzNn2HirUPOm5ubBZrT30Q0fK22Ls2t2Otz91331303b17+1bL5yjs8yHext\/8DeYn1mwmm\/hSeXqOLGdvm4JO3vfcc0\/B51Rftq+CqWavvMWAGMiCgdSHDykYTbFv375AMWkaP1Yc5O7e2L5NpuakEFrexO66VVitz8atvrt27cId+GDkxsb2hyR37bpj40N1BvMTazaTTXypPL3N5+ztuevkDQ8FkRN92b6SRkje3gKL5QhTawoxMBsG4vNj7dksUAsRAzNmYJSCy02Bt70AXMYSmyAGxEA1A32dI84fY9ls6NisLZkhAy1S4m+dKbQYQqEDMTBKwU3lzgHnoKd8sokBMdCMga7nyPpxBtH9bNh8W\/q0GODLEY4ePRpiHD9+XN9WtOKtXFnBXfG6Nb0YWHsGKLQgJqKpLe6ndh4M8AUJTz73QvDgCxNOnz6dR4L9ZzGZEVVwJ7NVSlQMiIHcGOCfAwHL682r2\/8syNqrkMUXJDz8aNjy0BcmrGIrdsw5SsHlJ2Z7m8p02ug7MsrUwN9EmqRGXApN+ipGDFQxwHnh3BBjOm10bMK4DHDOn3vuufDUU0+Fj33sY8Xk3\/3D7xbSvv+2aOhFDNxiYJSCe2uuwM0BWHtKkr+L1P0NhAN48uTJHX874W8p9J\/SepVrvgxwhkC+Ga5HZpx33qb9S3\/nmfATf2u74B54+wPF4t9y76FC6kUMeAZGKbjcHPgpPAZ2n0zO+ofe887A4arK0Q5g2d9QqvrKJwbqGOC8xGeINva6vvIPx8Db3\/3jgbdxmeGuvdsPZrDvv8UmiAFjYJSCa5NNWT724H2N07\/f\/+3E9MZ\/Q2k8jQLFgBgQA2JgQgyo4E5os5SqGBADYkAMTJeBUQpu6m0v3gbDPl3qlPkqGFjnOTkvnBvPAW3s3iZdDIiBPBkYrOByI\/Bg+ak2dkEMiIE0A\/7MoBOFNFgbKYgBMZA3A4MVXH7qboK86VF2YmC1DDQ5Q8RsZ6lXMSAGcmZgsIKb86KVmxgQA2JADIiBsRkYreDaW2Bejr1YzScGps6APz+mT31NY+WvecTAqhkYpeByY+BtrxjYV02A5hcDU2GA8xKfIdrYp7IG5SkG1pmBUQpu3wRzgwFNxm0a12QsxYgBMSAG5smAVjUGA5MruBRQfqoH6GUk4QNlftnFgBgQA2JADIzJwKQKLgWUQmsEoWOztpf4gLel9GvXroUUUrHYUrFmww9uXL8eru\/ADVxhCJ\/NL5neyxx4KTZfL2JADEySgb6SHqXgUvgojDGw97WQruOcOnUqnDlzZgE8E\/n8+fMLQ169erVonz17diHW97Wv6bp06XK4ePHiAq5cuVz0H8IX50T+5OJzm5o+9TXE+V+4cKHY\/2VeOC\/xGaKNfZlx1VcMiIFxGBil4NpNgRuDxzhLrJ7l2LFjYWtrawGHDx8OBw4cWOi4ubFZtA8dOrQQ6\/sePHiwiNm7b2\/Yt2\/fAvbs2bPt29u\/L86J\/GObz3MK+tTXEOe\/f\/\/+Yv+Xecn5HC2zLvUVA+vCwCgFd3QyW0x45MiRsPdmEfSgWFqBtKF2b2xT5eNSOvEbGxthc3NzAdgK32b\/vjgPyz+2T6k99TXE+XM9sP+CGBAD68vAdhUZeP38VstP5wNPo+HFwKwZ0Dma9fZqcWvAwCgF14otMkYbjuMbDmNhazOGYoMomCgDXO+kjoyBXRADYiBvBkYpuBTFMrSlh3HsZoPu+2P3beliYE4McL2XYU7r1FqGY8A+zMeHKtENw82okT0DoxRcP2Efut104rGwN7HFMWqLgdsMSBEDE2Dgh98\/F34QoUnaJ06cCE888UR46qmnCnn06NEAjh8\/Hii+TcZQTHcGBi24\/Mbp0T1N9RQD68uAP0Po68uEVm4MfOUzz4Yvf+L4An7nC79i7lL5+Mc\/FT74yV8Nf+2Zfx7++qc+G5587oXwvg8\/HU6fPl3aR47+GBis4HJj4DdOD2z9pa6RxMD8GeDM+DOEjm1FK9e0mTDw+LseKDLxkt96C2PFy\/0PPxru\/9H3hre\/+8fD1sPvvYlHw9Yjj1b0kKtPBgYruH0mqbHEgBgQA2LgDgM\/tnVv0YhlYdRLtgyo4Ga7NUpMDIiBSTGgZMVADQMquDUEyS0GxIAYEANioA8GVHD7YFFjiAExIAbEQBUD8t1kYNCCy4c7PG7OF3wbHZsgBsRAOQOcEw8ifRsdmyAGxEDeDAxWcPk0ZRPkTY+yEwOrZaDJGSKmLEuKMSjzm50YYG0k7RTwCWJADLRnoLTgth9KPcSAGMiJAYolxRigl+WGjxiA7uOwxfB+6WJADDRnQAW3OVeKFAOTYYDCSaG0hNGxWdskNnzWRsdm7aby2rVrYW6I1x6vz\/w3blwPgPaNGzcQN9vb8vr1bXnbfv16SNuwG271KWLNZvKWL5rHxmfyG1X9bvqIvX5T3sH2mPH6cmyzvilDBbfT7qmTGFgPBii+hqoVnzp1Kpw5c6YzeKwgz\/ddZoy++5IPa7569SoinD17dmF95r985Uq4cvlKEUMRQzH5xhtv0LxZZK8X8tKly+HKlcuFbjHYLl68GAzmj+34zWd9Y8k321JYAAAQAElEQVTAVf3Md\/ny5R3zxevrm89lx+Maef3111niZKGCO9mtU+JiYHgG+I3XQOEtm\/HYsWNha2urMw4fPhwOHTrUuf8yc5f1PXjwYLHczY3NQsb5mZ\/vzt6zd08Rs3v39i3V5N13371g5zupicdoMdj4\/mSD+WM7fvNZ31gyblU\/8yHvueee4MeM11fGy6rsXCPkzBqniu2rY6rZK+8sGVBS82CAQutXQrus6B45ciRwE+8Ku\/F37T9UP9a\/e2P7NpmaA\/\/G7o2w+ybQd+3ahQi7dm3LjY3di+3NjbCxsbHDtrm5GQzm37gZazaT5tu1a3v8XbsWJQNX9cPHGAbGRadfan052bhGyJdcp4rdU018VXnztoafm7bB28t0nnfKt3yU+WUXA2JADIiBeTKggttyX\/l6K\/sqKwrtyZMni6+34iuuPvnJT9aO9r\/+45eKb\/jg2z4ovABbbceKAPLw4G9LoKKLXCtjYJyJ499G+c0UWzw7NnxmR8dGGx1poG0+s0mKATHQnAEV3OZc3Y60r7KiyKF\/6D3vLHy0C6Xi5ey3vx74hg8kYfzGazrtLuCHAAq+wb7v8ud+7uf0HZddCJ1JH4ojRRKg+2VhszY+2gA9ZY99FiMpBsRAcwZUcJtzVRr52IP3lfpSDvuGj5Svi43vuHzyuReK77ZE8n2X7\/np4\/qOyy5kzqwPBRTEy4pttEEqDjuIfUO0NebqGOAXhhRWl9H8ZlbBncGe8h2XWw8\/Ggzb33f5F2awMi1BDIiBsRiI3ymzd8zsT2hj5THnebItuLyFBerIJwb4ONop+BjpYkAMiIFpMTBstvE7Zbxb9r4PP613ynqkPcuCS7HkLSyAXrZefMQAdB+HLYb3SxcDYkAMiIE7DMTvlBXvmD3y6J0AaUszkF3BpXBSKG1l6NisbRIbPmujY7N2U1n2+LK4P49D8zb6+bbXqx6tRpyNtR23\/Vi1wr7wuLXoUW4JH33A9jgWf\/32Y+Xwkaew+scOshdCvgzob5f57k3TzKYQl13B7Ys0iq+haszUI+k4fOfPn1\/o9uabby60eQya\/dMbe\/SbBdjj03gUm8EeyUaMPY6NOG+nbfEmzZ\/yXbp0ieHCpUt3HtN2p9\/2o+bIc9lHqq2iP3sAv6uYu4854\/wvXLhQ7JVe8mOAvfL\/vM\/+dsnfNMn21VdfRQhiYGkGZltw+Y3XQOEtYyr1SDoeIXbgwIGFLndtbj\/ezYw8Bs0e7WaPfjPf3n17i0em8WQUgz2SjZjdu7dp5yku3k7b4k2aP\/bxiLM9d28\/Ti720ffuu+5mquwel9f0sXDsARw3jc8tLs5\/\/\/79xX7oJT8GKLj88774b5j8\/ZJsr1zZ\/uGVOA98ghhow8D2nb9Nj75iBxyHQuuHp11WdI8kHklHwbJCZ+PsulUkrU2RA7R3byzSyKPSeASZBzZiwa5d249js8esYQO0fR906xf7sGNL9cNnOZHjFGF7MMXcyTnOn71kr4R8GdjxN8zo75f8xmu\/\/SJp57saZZYjA4uVIscMlZMYEANiIAMGyn4DziA1pTARBrIruPFvo\/xmii3mExs+s6Njo42ONNA2n9kyl72l598C83pvE2ggMbAmDNT9BrwmNGiZSzCQXcFlLRRHiiRAx2bAZjo+2gA9ZY99FrMukre9ePsrxk9\/5KN67OO6XARapxgQA1kwkGXBhRkKKED3iG20gY9Bx2agva6I3wZ78rkXAh8G+fbXvzZdSpS5GBADYmCCDGRbcCfIZZYp73gbjEdARh8GyTJxJSUGxIAYmBkDKrgz21AtZ+0ZEAFiQAxkyoAKbqYbo7TEgBgQA2JgXgyo4M5rP7UaMSAGqhiQTwyskAEV3BWSr6nFgBgQA2JgfRhQwe241\/yb1o5d1U0MiAExkCMDymlgBlRwOxLMv28FdH\/55ZcRghgQA2JADIiBUgZUcEupqXbw71v596xEveXeQwhBDIgBMSAG5spAD+tSwe1IYvHvW2\/9e9a9b\/2RjqOomxgQA2JADKwLAyq467LTWqcYEANioCMDfGYlhY7DrW23GRfctd3TxgtPHSBsjQdQoBgQA2vBAJ9XiZ\/HTvv48eN6JnuLK0AFtwVZcwvVIZrbjmo9YmAYBvjMCs9h9+AzLKdPnx5mwpmOqoI7042tWxZ+HSJYEMSAGKhjoPjMCs9h97j1GZa6vvLfYUAF9w4Xa6fpEK3dlmvBYkAMrJABFdwVkj\/HqfkbcBmms15lOkcGdF0Os6vitTmvKrgNuXrz6psNI9c3jIP387\/4S4EPU6SgD1is77Wx6pVzbZ48eTJ5bfJZhlXnN+X54U\/nvdkOquAmeDp37lz4\/K+fWvB89w+\/u9Cee4MbVApV6yb+21\/\/WtDfhqtYmr4v5xVwDZaBD\/iUXZs5ryn33Mo4he+y3Mv2CHtZnznYVXATu0jBfeVbf7DgOfD2Bxbac2+U\/dT60x\/5aPjd3\/3d4p8CcDg8jBP9bdiYkPTXR6x3ZScex7frfovVtdmV9fJ+bTllv8r2id+UuccQUz7jdD2zLrgPPfRQAH1sz1177+ljmMmMUfZTK7\/BlhVj7HUL5CClUNevzJ8ay2xlfWRvxwBnCLTrFYofyqpurF3+xMDeVo3Jb1Vl127b\/LvFq1cdA+xh1T5xj6kbY6r+2RZcbhCvvPJKAOhT3aAueV++fDm8dua1Ll1v97nff\/zf9Fv\/DCB1Q+Pf5\/Hv8m4PUKJQlPkpNkZ88+Vdhi996UsBWTLUIDf0srnMzs2iDBaDJO+6\/InLHZwdzhBAb5MvPFXdWPERA4yvr33ta8W+ls1DLP1S16Bdf1XXbtm4y9hfffWV8MYbV5YZYtS+3B++853vhEuXLo06bzxZ232yawQZjzWV9iwLLjcGbhC2CejYrO0lBziG95fpP\/z+uTJXwPeDm36PVDBx3n7x4sUQw\/yxnXaZjwN15rUzhZu4GIXj5ktsp33TXPyPHqNw3HzZ\/+C7w4EEtm4V5Lgf7Zvdiv8fO\/4Pwk\/94qcXwI2Sm6jfBw7Viy++WNx8vT3W6dd0zLhvl3bVb1f80GBvt8f5F4uf2AtnhrNjaaNjs7aXKS7Nn7pe7FqxH8CeeOKJwH5b23Ppx24yJtdbDOsX22mX+cxukvOaOtN\/ePPPT9evX7ewSskYcQC21LixnZi4b6rNmmJYHHbuD1ZwaXsfbY9lfX7v0EGTMYmLYWfK+k9RzrLgNtmId7zjHeH9739\/sAPuf+PCtnHvVviLDxy8jV2Xf3hb\/3N\/5i3hK595tsCuffuD+bAzN74vf+J48MCGz+KJxQbu3XdXePcDh8N\/\/+3fDv\/5y7+5gN++adv8s3++le+\/\/dZ\/YqrQth9zd5mPfqBJ35f\/+P+G3\/vGNxfwrZs2eIF32wd0FoE0Wyzx0a9sTPaQmLjfMm0K\/J4f\/6kQ4653PRZ+\/7ULt68n5vX5f+5zn6M5OzQ5R1XXdcwjbbj8vT\/41m0u\/X7BK3teNWZfPq5nzuV9B\/YHO6+pM33l934r\/I9\/9Y8C59juBUj6\/M4XfiWQr7U578Bi77v33uI+khqXOG9Hx2Z9uT8xrpdNzjz3hxuXfxCQTc9tF07Lzh97yEGoGpMYv+\/o2Lhnc83Rf4pY64L7\/PPPh5deeimJf\/\/5fxme\/dV\/exs\/+Uv\/4rb+z\/7NF273+c3f+LVgPm8vG9fifeznf+NL4Zd\/9YXw67\/wN3bgi\/\/wRHjx+X+yw05sTr5l84GXMs6q7PRj7hhwwx5W9e3q+3fP\/XyI8aV\/+olin8rG\/Jmf+Zkp3h9qc+bmV3eO4r2hzf5wXcc80oZL9rWMS3yMEcPGjO20u\/jow7n8q899LvjzWpYX59juBUjrQ76+TX+L\/def\/\/ztewn2JrC+3J8Y10s4Zb0xWMvYvrrzF+dI2\/Is44FrrfaizDiguuBmnHgfqXGz+MAHPhAEcTD0NcC11sc1m+MYrG1o\/jS+zijXANdajmegaU5rXXCbkqQ4MSAGxIAYEAPLMjDLght\/uIMPemBblqyov5piYNYMcGY4O7ZIdGzWlhQDYqAdA7MsuFDAjYEbBEDHJogBMdCOAc4OZwigt+utaDEgBjwDsy24LJIbBEBvA24uoE2fMWPJLQWfg\/m9zfQqn8UsLUsGYO6UCzvo05caqw9bnCftFPxc5vc206t8FrNKyRkCY+RgXMRyjLm7zGF5dum7ij6WbyxXkUuTOckzFYcdpHw522ZdcLsQzyZycwHoXcYYow\/5xbB5ydt86GZH0i7z4R8KzAtS42Mvy6mrLzXPsjZyAalxLH8vLY4+Zkc3O5J2mQ\/\/OsL48DJHHqa6d55X03PjF25BKi\/sljd6KiZXmwqu2xk2j400Ezo2a09Bki95W67o2GgjaaMDdGzoQ4O5QDwP83s7OjbikLTRATo2dCRtdICODX0oMAdoMz45+T7oDz30UDFEla8I0Eu2DGjvht0azgmIZ5k67yq48Y5OpM2FZ5hIyrNO0\/YCOeuFjrA4ODSMMN3aTWHcItdu8StesAruijeg6\/T89GfQwenKYn\/9bC+Q2o9qXuEnBev1yq1noItLY6RfCa8G9qHf0TVaFQMquFXsZOrjsPjUaOvgeEbG1eHfz0hb++EZWdThJwWisCMNtMWlsbG8hE8\/Cm3x6xkZVlfBHZZfjS4GxIAY6JMBjTVhBlRw3ebFP+3xkx82F5KFSl4+EdqWJ5K2+dGx0UbSRgfo2NBXBeYnD5sfHRttJG10gI4NHUkbHaBjQx8bzO3npG25IGmbHx0bbSRtdICODX1dAQd+7bRz5IScyM1yRcdm7VwlefrcaE8hb8uZXMnZ2ujYrJ27VMGNdojNYxMBeuTOokle5Geg7ROj3cXnxxhT75pvVb9V5Q\/v5OXnp40doDf1+bh10eEHngy0c107uU0hT8+fz5ncaXv\/FHRyJneAvpBz5g0V3MQGsYkg4crGRH6GVFJdfamx+rSRV2o87KBPX2qsPmypPLEZUnN09aXGmrvNuELmvlZyBLnn6fMjX4O356iTZyov7CDly9mmgpvz7ig3MSAGxIAYmA0DKy64s+FRC0kwwFs+CXMnU59jdUpAncRACwbmdL3OaS0ttnCQUBXcQWjVoBzSPt\/yYSzGFLNiIHcGuE65XnPPs2l+rIU1NY1XXDkDKrjl3KzcU5dA3SGo8zO+xZjEZmhqs\/gqyVgpVPWRTwy0YYDrq018k1jGTKFJX4uhv+mrlLnksUoOVj23Cu6qd2DJ+csOUZl9yemK7nVj4+en4iLYvWCLQawLqVTp2ya+cjA5xUBDBrjuYpRdh9iJbTj0ysLIs83krKltnzbjr0usCu667HTLdXK4xjhk3edouSCFi4E1Z4CztuYUrHz5Krgr34LlEuAQURz9KLSxm422h9m9JJ4Yb4t170cHccwybcbzWGYs9RUDMQP+2kKv86di4j5t24zpEff3PtNTMWU2+uBDGmgD2rE0m9lpG7AJ\/TKgtIk1fgAAAhtJREFUgtsvn9mNxuGhmHpg65oo49AXCdD7ADkxnge2PsZexzG05kUGuJb8tYWOzaLQsRmwoyP7QjwH42Oz8dGxGbCjI9sgHoe2729jIgE+YtA9sOET+mNABbc\/LmczEgeNg8eCkLTRhwTzDDm+xhYDXRng+o+R8\/XaV259jdOV9zn2U8Gdwa5yMLghsBQkbfSpgdw9ppa\/8p0KA+3y5DzFaDdC\/tGsz5899Pyznl6GKrjT27NZZswB59B7zHKhWlSWDHDdcQ0aaI+dKHPa\/EjaY+bAfB7kMOb86zCXCu467HLDNXLYOGRI34V2yu5jvG7x3laml42LPe6DjbFju9piYFkG7Nri+gJdx6MvY3XpTz\/6G5qOQb+msWVxfoxUDH7ySvlka86ACm5zrrKO5DCkDoXZ8RmwLbMY+ttYbcaxPl4yFmMgU3Z8ghhoyoC\/hrweX1\/4sPlxsXl4Xx868\/nx0bH5sbF5eJ\/p9PExtM3XVNLHxqCPb5sdGz6hPwZUcPvjcvSR4gNR1sbu4RPFXtU2XyoutlksEh8HFx3QTgGfwfux0UYCxvJtbIIY8AxwfZSBuNiHDdi1FfvNTkwbMA59rQ9tr9M2mJ14s3mJ3WK89DHYaXuJbijzYQc+jrbB7OSAzdrrLZdbvQrucvypdwUDfR7SPseqSFmuNWSAa4uiEgN7Vzra9iU+np829q459NUvhxz6Wsuqx\/lTAAAA\/\/+RYgeYAAAABklEQVQDABGZe7ZCzG5vAAAAAElFTkSuQmCC","height":229,"width":381}}
%---
%[output:6816dab6]
%   data: {"dataType":"text","outputData":{"text":"[23:35:34][INFO]  Checkpoint saved -> result\/r05_checkpoint.mat\n","truncated":false}}
%---
%[output:4dca65b3]
%   data: {"dataType":"text","outputData":{"text":"[23:35:34][INFO]    Load in R06: load(\"result\/r05_checkpoint.mat\")\n","truncated":false}}
%---
%[output:4752e9cd]
%   data: {"dataType":"text","outputData":{"text":"[23:35:34][INFO]  Corpus:         500 molecules Valid FDA drug SMILES (trained with 500 molecules after expansion)\n","truncated":false}}
%---
%[output:86afdc49]
%   data: {"dataType":"text","outputData":{"text":"[23:35:34][INFO]  Vocabulary:     37 characters  |  Max sequence length: 100\n","truncated":false}}
%---
%[output:5038645d]
%   data: {"dataType":"text","outputData":{"text":"[23:35:34][INFO]  Loss plateau:   val=2.6101 (epoch 70+)  --  Random baseline=3.6109  (Perplexity=13.6)\n","truncated":false}}
%---
%[output:134573c2]
%   data: {"dataType":"text","outputData":{"text":"[23:35:34][INFO]  --- Generation Metrics (N=300 each; Binomial SE~1.7%)\n","truncated":false}}
%---
%[output:6104ce8f]
%   data: {"dataType":"text","outputData":{"text":"[23:35:34][INFO]  Markov Model:   Validity=9.0%  Diversity=66.7%  Novelty=92.6%\n","truncated":false}}
%---
%[output:97f77c7d]
%   data: {"dataType":"text","outputData":{"text":"[23:35:34][INFO]  LSTM Model:     Validity=5.3%  Diversity=68.8%  Novelty=100.0%\n","truncated":false}}
%---
%[output:7a8ab334]
%   data: {"dataType":"text","outputData":{"text":"[23:35:34][INFO]  Warning: With N=300, the validity gap is not statistically significant (p~0.10)\n","truncated":false}}
%---
%[output:198ebf94]
%   data: {"dataType":"text","outputData":{"text":"[23:35:34][INFO]  Architecture:   LSTM(128) -> Dropout(0.30) -> FC(37) -> Softmax\n","truncated":false}}
%---
%[output:17de43a5]
%   data: {"dataType":"text","outputData":{"text":"[23:35:34][INFO]  Number of Parameters:     79621  (Overparameterized relative to 18000 training tokens)\n","truncated":false}}
%---
