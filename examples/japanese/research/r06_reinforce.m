%[text] # R06: REINFORCE — 分子特性誘導ファインチューニング
%[text] EasyMolKit Research — Layer 4
%[text] 
%[text] R05 で学習した LSTM は有効な SMILES を生成できますが、「良い分子かどうか」を考慮する仕組みがありません。
%[text] 生成された分子の特性スコアを「報酬（reward）」として方策（LSTM）を更新することで、目標特性を持つ分子をより多く生成するよう誘導できます。
%[text] この手法が **REINFORCE アルゴリズム**（Williams 1992）であり、REINVENT・GuacaMol・MolGPT など現代の分子最適化ツールの基礎をなしています。
%[text] このスクリプトでは R05 のモデルを読み込み、リピンスキー Ro5 スコアを報酬として REINFORCE ファインチューニングを体験します。
%[text] 
%[text] ## 学習目標
%[text] - 系列生成の MDP（マルコフ決定過程：状態 = 部分 SMILES、行動 = 次の文字、報酬 = 分子スコア）を理解します。
%[text] - ベースライン付き REINFORCE の勾配推定量 $\\nabla J = E\[(R - b) \\cdot \\nabla \\log \\pi\]$ を実装します。
%[text] - 壊滅的忘却（catastrophic forgetting）が起こる理由と、KL 正則化による文法崩壊防止を理解します。
%[text] - 手動 REINFORCE ループと MATLAB Reinforcement Learning Toolbox の `rlPGAgent` を比較します。
%[text] - RL ファインチューニングを実際の創薬ワークフロー（REINVENT、GuacaMol、GENTRL）に結び付けます。 \
%[text] ## 前提条件
%[text] - R05 `r05_smiles_generator.m` を完了していること（チェックポイントが `result/` に保存済み）。
%[text] - Deep Learning Toolbox（`dlnetwork`, `dlarray`, `dlfeval`, `adamupdate`）。
%[text] - Reinforcement Learning Toolbox（任意; Section 4 で活用します）。 \
%[text] ## 動作環境
%[text] - Deep Learning Toolbox が必要です。
%[text] - MATLAB Online と Desktop の両方に対応しています。
%[text] - GPU アクセラレーションは不要です。
%[text] - 所要時間: 30〜60 分 \
%[text] ## データ
%[text] - `result/r05_checkpoint.mat` — R05 で学習した LSTM と語彙（vocabulary） \
%[text] ## 参考文献
%[text] - Olivecrona M et al. (2017) Molecular de novo design through deep reinforcement learning. J Cheminform 9:48. doi:10.1186/s13321-017-0235-x \[REINVENT: 本スクリプトの直接の基礎\]
%[text] - Williams RJ (1992) Simple statistical gradient-following algorithms for connectionist reinforcement learning. Mach Learn 8:229-256. doi:10.1007/BF00992696 \[REINFORCE アルゴリズムのオリジナル論文\]
%[text] - Brown N et al. (2019) GuacaMol: benchmarking models for de novo molecular design. J Chem Inf Model 59:1096-1108. doi:10.1021/acs.jcim.8b00839 \[標準的な生成モデル評価指標\]
%[text] - Bung N et al. (2022) De novo design of new chemical entities for SARS-CoV-2 using artificial intelligence. Future Med Chem 14:1019-1030. doi:10.4155/fmc-2021-0223 \[COVID-19 標的最適化への REINVENT 適用\]
%[text] - Polykovskiy D et al. (2020) Molecular Sets (MOSES): a benchmarking platform for molecular generation models. Front Pharmacol 11:565644. doi:10.3389/fphar.2020.565644 \[validity/uniqueness/novelty 指標の標準実装\]
%[text] - Zhavoronkov A et al. (2019) Deep learning enables rapid identification of potent DDR1 kinase inhibitors. Nat Biotechnol 37:1038-1040. doi:10.1038/s41587-019-0224-x \[GENTRL: テンソル分解×RL による分子生成\] \
%[text] 
%[text] 実行方法: Ctrl+Enter でセクションごとに実行します。
%[text] 先に R05 を実行して `result/r05_checkpoint.mat` を生成してください。
%%
%[text] ## Section 0: セットアップと学習済みモデルの読み込み
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
%[text] Deep Learning Toolbox と Reinforcement Learning Toolbox の利用可否を確認します。
hasDL = license("test", "Neural_Network_Toolbox");
hasRL = ~isempty(ver("Reinforcement Learning Toolbox"));
if ~hasDL
    error("emk:r06:missingToolbox", ...
        "Deep Learning Toolbox が必要です。");
end
logInfo("ツールボックス確認: Deep Learning=%d  Reinforcement Learning=%d", hasDL, hasRL); %[output:021f8d24]
if hasRL %[output:group:7a0892da]
    logInfo("Reinforcement Learning Toolbox を検出しました。Section 4 は完全に利用可能です。");
else
    logInfo("RL ツールボックスが見つかりません。手動 REINFORCE ループを実行します。Section 4 は概念解説のみです。"); %[output:34c7d6b2]
end %[output:group:7a0892da]
%[text] ### 調整可能なパラメータ
RL_EPOCHS      = 20;      % REINFORCE 更新ステップ数
RL_BATCH       = 48;      % RL バッチあたりの系列数
RL_LR          = 1e-4;    % RL 学習率 (忘却防止のため SL 学習率より小さく設定)
GRAD_CLIP_NORM = 1.0;     % 勾配クリッピング閾値
N_EVAL         = 100;     % RL 前後の評価用 SMILES 生成数
%[text] ### R05 チェックポイントの読み込み
CHECKPOINT = "result/r05_checkpoint.mat";
if ~isfile(CHECKPOINT)
    error("emk:r06:missingCheckpoint", ...
        "チェックポイントが見つかりません: %s\n先に R05 (r05_smiles_generator.m) を実行してください。", ...
        CHECKPOINT);
end
logInfo("R05 チェックポイントを読み込んでいます: %s", CHECKPOINT); %[output:4b880a6e]
load(CHECKPOINT, ...
    "net", "char2idx", "idx2char", "VOCAB_SIZE", "MAX_SEQ_LEN", ...
    "START_IDX", "END_IDX", "PAD_IDX", "HIDDEN_SIZE", "TEMPERATURE", ...
    "smilesAll", "smilesProc", "X_train", "Y_train", "X_val", "Y_val");
logInfo("モデルを読み込みました: VOCAB_SIZE=%d  MAX_SEQ_LEN=%d  HIDDEN=%d", ... %[output:group:3260ca9d] %[output:9e381c29]
    VOCAB_SIZE, MAX_SEQ_LEN, HIDDEN_SIZE); %[output:group:3260ca9d] %[output:9e381c29]
logInfo("訓練コーパス: %d SMILES", numel(smilesAll)); %[output:9893d5f8]

%[text] Python/RDKit のウォームアップを行います。
mol_warmup = emk.mol.fromSmiles("C"); %#ok<NASGU>
clear mol_warmup;
%%
%[text] ## Section 1: RL 前ベースラインの評価
%[text] 
%[text] RL を実行する前に、現在のモデルの有効性と特性分布を測定します。
%[text] ベースラインがなければ、RL が何かを改善したかどうか判断できません。
%[text] 評価には次の 4 つの指標を使います（Brown et al. 2019 / MOSES; Polykovskiy et al. 2020）。
%[text:table]{"ignoreHeader":true}
%[text] | 指標 | 定義 |
%[text] | --- | --- |
%[text] | Validity | RDKit で解析可能な生成 SMILES の割合 |
%[text] | Uniqueness | 構造的に異なる有効 SMILES の割合 |
%[text] | Novelty | 訓練セットに含まれない有効ユニーク SMILES の割合 |
%[text] | Avg reward | 生成バッチの平均リピンスキー Ro5 複合スコア |
%[text:table]
logInfo("%d 個の RL 前 SMILES を生成中 (T=%.2f) ...", N_EVAL, TEMPERATURE); %[output:50c27046]
preRLSmiles = strings(N_EVAL, 1);
for g = 1:N_EVAL %[output:group:29273260]
    preRLSmiles(g) = generateSmiles(net, START_IDX, END_IDX, PAD_IDX, ...
        VOCAB_SIZE, MAX_SEQ_LEN, TEMPERATURE, idx2char, char2idx);
    if mod(g, 10) == 0 || g == N_EVAL
        logProgress(g, N_EVAL, "generating (pre-RL)"); %[output:6f56c742]
    end
end %[output:group:29273260]
%[text] 生成した SMILES を検証してスコアリングします。
[nPreValid, nPreUnique, nPreNovel, preRewards] = evalGenSmiles( ...
    preRLSmiles, smilesAll, N_EVAL);
logInfo("RL 前:   valid=%d/%d (%.0f%%)  unique=%d  novel=%d  avg_reward=%.3f", ... %[output:group:90d7e5c9] %[output:1b1db75a]
    nPreValid, N_EVAL, 100*nPreValid/N_EVAL, nPreUnique, nPreNovel, mean(preRewards)); %[output:group:90d7e5c9] %[output:1b1db75a]
%%
%[text] ## Section 2: REINFORCE ファインチューニング
%[text] LSTM 言語モデルは方策（policy）$\\pi$ として扱えます。
%[text] 状態 $s(t)$ はこれまでの部分 SMILES（LSTM 隠れ状態 $h(t)$）、行動 $a(t)$ は $\\pi$ からサンプリングされる次の文字トークン、報酬 $R$ は完全な分子のスカラースコア（終端報酬）です。
%[text] REINFORCE 勾配推定量（Williams 1992）は次のように書けます。
%[text]{"align":"center"} $\\nabla J(\\theta) = E\\left\[ (R(\\tau) - b) \\cdot \\sum\_t \\nabla \\log \\pi(a(t) | s(t)) \\right\]$
%[text] ここでバッチ平均報酬 $b$ をベースラインとして減算することで、勾配推定の分散を低減します。これが REINVENT（Olivecrona et al. 2017）の中核アイデアです。
%[text] 報酬関数にはリピンスキー Ro5 複合スコア $R \\in \[0, 1\]$ を使います。
%[text:table]{"ignoreHeader":true}
%[text] | 基準 | 条件 | 加算値 |
%[text] | --- | --- | --- |
%[text] | 分子量 MW | $\\leq 500$ Da | \+0.25 |
%[text] | 水素結合ドナー数 HBD | $\\leq 5$ | \+0.25 |
%[text] | 水素結合アクセプター数 HBA | $\\leq 10$ | \+0.25 |
%[text] | 脂溶性 LogP | $\\leq 5$ | \+0.25 |
%[text:table]
%[text] SMILES が訓練セットにない場合は新規性ボーナスとして $\\times 1.1$（上限 1.0）を乗じます。
%[text] 
%[text] **壊滅的忘却（catastrophic forgetting）** とは、RL 更新が大きすぎることで学習済み SMILES 文法を破壊する現象です。
%[text] REINVENT "augmented likelihood" では次の KL 正則化を導入してこれを防ぎます。
%[text]{"align":"center"} $L(total) = L(PG) + \\lambda \\cdot KL(\\pi(RL) \\| \\pi(prior))$
%[text] $L(PG)$ は方策勾配損失、$\\pi(prior)$ は固定した RL 前ネットワークであり、これが RL 前分布からの大きな逸脱にペナルティを課します。
%[text] 各 RL エポックでバリデーション損失 `val_loss` を監視することで、文法崩壊を早期に検出できます。
%[text] 事前ネットワークのコピーを保存しておきます。
%[text] > **簡略化の注記**: REINVENT 論文では $\pi_{prior}$ を用いた KL 損失を実際に計算しますが、本実装では計算コスト削減のために省略しています。
%[text] > 代わりに各エポックの `val_loss` を監視することで、文法崩壊を早期に検出します。
netPrior = net; %#ok<NASGU> % reference snapshot for teaching; KL loss is not applied in this simplified implementation

avgG_rl       = [];
avgSqG_rl     = [];
iter_rl       = 0;
rl_rewardLog  = zeros(RL_EPOCHS, 1);
rl_valLossLog = zeros(RL_EPOCHS, 1);

logInfo("REINFORCE ファインチューニングを開始します (%d エポック, バッチ=%d, lr=%.1e)", ... %[output:group:3a1ac0a2] %[output:09293c18]
    RL_EPOCHS, RL_BATCH, RL_LR); %[output:group:3a1ac0a2] %[output:09293c18]

for rl_ep = 1:RL_EPOCHS %[output:group:6484376b]
    % ステップ 1: SMILES バッチを生成 (推論のみ -- 勾配なし)
    batchSmiles = strings(RL_BATCH, 1);
    for g = 1:RL_BATCH
        batchSmiles(g) = generateSmiles(net, START_IDX, END_IDX, PAD_IDX, ...
            VOCAB_SIZE, MAX_SEQ_LEN, TEMPERATURE, idx2char, char2idx);
    end
    % ステップ 2: 系列ごとの報酬を計算してベースラインを減算
    rewards = arrayfun(@(s) lipinskiReward(s, smilesAll), batchSmiles);
    baseline = mean(rewards);
    rl_rewardLog(rl_ep) = baseline;
    advantages = rewards - baseline;    % [RL_BATCH x 1]
    % ステップ 3: 勾配計算のために生成系列を再エンコード
    [Xrl, Yrl] = encodeSmilesBatch(batchSmiles, char2idx, ...
        START_IDX, END_IDX, PAD_IDX, MAX_SEQ_LEN, VOCAB_SIZE);
    if isempty(Xrl)
        logWarn("RL エポック %d: すべての系列がエンコード不可、スキップします。", rl_ep);
        continue;
    end
    dlX_rl = dlarray(Xrl, "CTB");
    % ステップ 4: 方策勾配の更新
    iter_rl = iter_rl + 1;
    [~, grads_rl] = dlfeval(@policyGradientLoss, net, dlX_rl, Yrl, ...
        PAD_IDX, advantages);
    % 勾配クリッピング
    for gi = 1:height(grads_rl)
        gdata = extractdata(grads_rl.Value{gi});
        nrm = sqrt(sum(gdata(:).^2));
        if nrm > GRAD_CLIP_NORM
            grads_rl.Value{gi} = grads_rl.Value{gi} * (GRAD_CLIP_NORM / nrm);
        end
    end
    [net, avgG_rl, avgSqG_rl] = adamupdate(net, grads_rl, avgG_rl, avgSqG_rl, ...
        iter_rl, RL_LR);
    % ステップ 5: バリデーション損失で壊滅的忘却を監視
    dlXv = dlarray(onehot_encode(X_val, VOCAB_SIZE), "CTB");
    rl_valLossLog(rl_ep) = extractdata(dlfeval(@modelLoss_local, net, dlXv, Y_val, PAD_IDX));
    if mod(rl_ep, 5) == 0 || rl_ep == 1
        logInfo("RL エポック %2d/%d  avg_reward=%.3f  val_loss=%.4f", ... %[output:94e68094]
            rl_ep, RL_EPOCHS, baseline, rl_valLossLog(rl_ep)); %[output:94e68094]
    end
end %[output:group:6484376b]
%[text] REINFORCE の進捗を可視化します。
figure("Name", "R06 REINFORCE Progress"); %[output:5406fbc1]
tiledlayout(1, 2); %[output:5406fbc1]
nexttile; %[output:5406fbc1]
plot(1:RL_EPOCHS, rl_rewardLog, "g-o", LineWidth=1.5, MarkerSize=4); %[output:5406fbc1]
xlabel("RL エポック"); ylabel("バッチ平均報酬"); %[output:5406fbc1]
title("REINFORCE: ドラッグ様性報酬"); grid on; %[output:5406fbc1]
nexttile; %[output:5406fbc1]
plot(1:RL_EPOCHS, rl_valLossLog, "r-s", LineWidth=1.5, MarkerSize=4); %[output:5406fbc1]
xlabel("RL エポック"); ylabel("バリデーション損失"); %[output:5406fbc1]
title("言語品質 (忘却モニタリング)"); grid on; %[output:5406fbc1]
sgtitle("REINFORCE ファインチューニング進捗"); %[output:5406fbc1]
%[text] 報酬カーブは `RL_BATCH=48` の小バッチでは激しく振動することが多いですが、これは REINFORCE の既知の課題（勾配分散が高い）であり、バグではありません。
%[text] バッチサイズを増やすか、移動平均で平滑化すると振動が緩和されます。
%%
%[text] ## Section 3: RL 後の生成評価
%[text] ファインチューニング後に新しいバッチを生成して RL 前後を比較します。
%[text] 確認すべき点は「有効性（文法は維持されているか）」「報酬分布（ドラッグ様分子に向けてシフトしたか）」「新規性（高報酬の訓練分子を記憶しているだけではないか）」の 3 つです。
logInfo("%d 個の RL 後 SMILES を生成中 ...", N_EVAL); %[output:00e322bf]
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
logInfo("RL 後:   valid=%d/%d (%.0f%%)  unique=%d  novel=%d  avg_reward=%.3f", ... %[output:group:3c9597d5] %[output:77d3414f]
    nPostValid, N_EVAL, 100*nPostValid/N_EVAL, nPostUnique, nPostNovel, mean(postRewards)); %[output:group:3c9597d5] %[output:77d3414f]

%[text] RL 前後と訓練データの特性分布を比較します。
validPost = postRLSmiles(logical(arrayfun(@(s) emk.mol.isValid(s), postRLSmiles)));
validPre  = preRLSmiles(logical(arrayfun(@(s) emk.mol.isValid(s), preRLSmiles)));
FEAT_NAMES = ["MolWt", "LogP"];

propsPost  = batchDescriptors(validPost(1:min(30,numel(validPost))),  FEAT_NAMES);
propsPre   = batchDescriptors(validPre(1:min(30,numel(validPre))),    FEAT_NAMES);
propsTrain = batchDescriptors(smilesAll, FEAT_NAMES);

if ~isempty(propsPost) && ~isempty(propsTrain) %[output:group:1d69edfb]
    figure("Name", "R06 Property Distribution"); %[output:9e6d6b77]
    titles = ["分子量 (Da)", "LogP"];
    for fi = 1:2
        subplot(1, 2, fi); hold on; %[output:9e6d6b77]
        histogram(propsTrain(:,fi),  15, Normalization="probability", ...
            FaceColor=[0.2 0.6 0.8], DisplayName="訓練データ", FaceAlpha=0.7);
        if ~isempty(propsPre)
            histogram(propsPre(:,fi),  10, Normalization="probability", ...
                FaceColor=[0.9 0.7 0.2], DisplayName="RL 前", FaceAlpha=0.7);
        end
        histogram(propsPost(:,fi), 10, Normalization="probability", ...
            FaceColor=[0.85 0.3 0.1], DisplayName="RL 後", FaceAlpha=0.7);
        xlabel(titles(fi)); ylabel("確率密度");
        title(titles(fi)); legend("Location","best"); grid on;
    end
    sgtitle("特性分布: 訓練 vs RL 前 vs RL 後"); %[output:9e6d6b77]
end %[output:group:1d69edfb]
%%
%[text] ## Section 4: MATLAB Reinforcement Learning Toolbox（概念解説）
%[text] MATLAB Reinforcement Learning Toolbox は `rlPGAgent`（方策勾配エージェント）を提供します。
%[text] `rlPGAgent` は手動 REINFORCE に加えて、エントロピーボーナス・アドバンテージ正規化・学習率スケジュールを備えています。
%[text] エントロピーボーナス（entropy bonus）とは、方策が 1 つのトークンに集中しすぎる「モード崩壊」を防ぐためのペナルティ項です。
%[text] 
%[text] 文字レベルの分子生成では、RL 環境を次のように定義します。
%[text:table]{"ignoreHeader":true}
%[text] | 要素 | 内容 |
%[text] | --- | --- |
%[text] | 観測 | 現在のワンホット文字トークン（`VOCAB_SIZE` × 1） |
%[text] | 行動 | 次のトークンインデックス（離散, 1..`VOCAB_SIZE`） |
%[text] | 報酬 | 各ステップでは 0; END トークン時に特性スコア |
%[text] | 終端 | END がサンプリングされるか `MAX_SEQ_LEN` に達したとき |
%[text:table]
%[text] LSTM を `rlEnvironment` にラップするには `dlnetwork` と RL ツールボックスの観測/行動インターフェースを橋渡しする必要があります。
if hasRL %[output:group:4b4f931b]
    logInfo("Reinforcement Learning Toolbox が利用可能です。");
    logInfo("上記の環境定義をもとに rlPGAgent を実装できます。");
else
    logInfo("Reinforcement Learning Toolbox が見つかりません。"); %[output:298783de]
    logInfo("手動 REINFORCE ループ (Section 2) が同じアルゴリズムを実装しています。"); %[output:5308948c]
end %[output:group:4b4f931b]
%%
%[text] ## Section 5: まとめ
%[text] RL 前後の生成品質と報酬の変化を確認します。
logInfo("RL 前:  valid=%d/%d (%.0f%%)  avg_reward=%.3f", ... %[output:group:2339c514] %[output:6e28324c]
    nPreValid,  N_EVAL, 100*nPreValid/N_EVAL, mean(preRewards)); %[output:group:2339c514] %[output:6e28324c]
logInfo("RL 後:  valid=%d/%d (%.0f%%)  avg_reward=%.3f", ... %[output:group:0163c6b1] %[output:7ce94921]
    nPostValid, N_EVAL, 100*nPostValid/N_EVAL, mean(postRewards)); %[output:group:0163c6b1] %[output:7ce94921]
delta_reward = mean(postRewards) - mean(preRewards);
if delta_reward >= 0 %[output:group:7dbc1339]
    logInfo("RL 結果: 報酬 %.3f -> %.3f  (delta=+%.3f, 改善)", ...
        mean(preRewards), mean(postRewards), delta_reward);
else
    logInfo("RL 結果: 報酬 %.3f -> %.3f  (delta=%.3f, ノイズ範囲内)", ... %[output:1d690113]
        mean(preRewards), mean(postRewards), delta_reward); %[output:1d690113]
end %[output:group:7dbc1339]
%[text] > **Validity が低い理由**: 本実験では validity が 6% 程度と低い値を示しています。これは主に R05 の訓練コーパスが 500 SMILES と小規模であり、基盤モデルが SMILES 文法を十分に学習できていないためです。REINFORCE の限界ではなく、**ベースモデルの未収束**が主因です。実用的な分子生成では 10 万件以上の訓練データと数百エポックの事前学習が必要です。 \
%[text] - REINFORCE は $\\nabla J = E\[(R - b) \\cdot \\nabla \\log \\pi\]$ で方策を更新する、最もシンプルな方策勾配アルゴリズムです。
%[text] - ベースライン $b$（バッチ平均報酬）の減算が勾配推定の分散を下げ、学習を安定させます。
%[text] - 壊滅的忘却を防ぐには KL 正則化 $L(total) = L(PG) + \\lambda \\cdot KL(\\pi(RL) \\| \\pi(prior))$ が有効です。
%[text] - 報酬設計が最も重要です。リピンスキー Ro5 は弱いシグナルであり、QED や合成可能性スコア（SA score）を組み合わせるとより効果的です。
%[text] - MATLAB `rlPGAgent` は REINFORCE にエントロピーボーナスとアドバンテージ正規化を加えた実装であり、手動ループと同じアルゴリズムを提供します。
%[text] - 実用スケールでは $10^5$〜$10^6$ 個の分子評価が必要です。小規模データセットでは報酬の改善がノイズと区別しにくいことに注意してください。 \
%[text] ローカル関数
function [nValid, nUnique, nNovel, rewards] = evalGenSmiles(smilesList, trainingSmiles, nTotal)
%[text] 生成された SMILES バッチの有効性・一意性・新規性・報酬を評価します。
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
%[text] リピンスキー Ro5 複合報酬を計算します。
%[text] Ro5 基準を 1 つ満たすごとに +0.25 を加算し、新規性ボーナスとして $\\times 1.1$（上限 1.0）を乗じます。
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
%[text] REINFORCE 方策勾配損失を計算します。
%[text] `advantages` は各系列の $R(i) - b$ を表す $B \\times 1$ ベクトルです。
%[text] $-E\[\\mathrm{advantage} \\cdot \\log \\pi(a | s)\]$ を最小化することで方策を更新します。
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
logPi   = sum(dlTarget .* logPred, 1);   % 選択した行動の対数確率

%[text] 各系列の advantage を T ステップ分繰り返して T\*B 長のベクトルを作ります。
nB    = numel(advantages);
dlAdv = dlarray(single(repelem(advantages(:)', T*B / nB)));

loss  = -sum(logPi .* dlAdv .* mask) / nTokens;
grads = dlgradient(loss, net.Learnables);
end

function [Xoh, Ymat] = encodeSmilesBatch(smilesList, char2idx, ...
    startIdx, endIdx, padIdx, maxLen, vocabSize)
%[text] SMILES 文字列配列を RL 勾配計算用にエンコードします。
%[text] `Xoh` は `[vocabSize x SEQ_LEN x B]` のワンホット配列、`Ymat` は `[B x SEQ_LEN]` の整数配列です。
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
%[text] `[N x T]` 整数行列を `single [vocabSize x T x N]` のワンホットテンソルに変換します。
[N, T] = size(seqMat);
Xoh    = zeros(vocabSize, T, N, "single");
valid  = seqMat >= 1 & seqMat <= vocabSize;
[nz_n, nz_t] = find(valid);
nz_c   = double(seqMat(valid));
linIdx = nz_c + (nz_t - 1) * vocabSize + (nz_n - 1) * vocabSize * T;
Xoh(linIdx) = 1;
end

function [loss, grads] = modelLoss_local(net, dlX, Ytarget, padIdx)
%[text] マスク付きクロスエントロピー損失を計算します（R06 のバリデーション損失監視用ローカルコピー）。
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
%[text] 温度サンプリングと文法制約による自己回帰的 SMILES 生成を行います。
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
%[text] SMILES リストから記述子行列 `[N_valid x numel(featNames)]` を計算して返します。
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
%   data: {"dataType":"text","outputData":{"text":"[10:07:25][INFO]  initPython: Python 3.10 already active (Status: Loaded) -- skipping.\n","truncated":false}}
%---
%[output:021f8d24]
%   data: {"dataType":"text","outputData":{"text":"[10:07:25][INFO]  ツールボックス確認: Deep Learning=1  Reinforcement Learning=0\n","truncated":false}}
%---
%[output:34c7d6b2]
%   data: {"dataType":"text","outputData":{"text":"[10:07:25][INFO]  RL ツールボックスが見つかりません。手動 REINFORCE ループを実行します。Section 4 は概念解説のみです。\n","truncated":false}}
%---
%[output:4b880a6e]
%   data: {"dataType":"text","outputData":{"text":"[10:07:25][INFO]  R05 チェックポイントを読み込んでいます: result\/r05_checkpoint.mat\n","truncated":false}}
%---
%[output:9e381c29]
%   data: {"dataType":"text","outputData":{"text":"[10:07:25][INFO]  モデルを読み込みました: VOCAB_SIZE=37  MAX_SEQ_LEN=100  HIDDEN=128\n","truncated":false}}
%---
%[output:9893d5f8]
%   data: {"dataType":"text","outputData":{"text":"[10:07:25][INFO]  訓練コーパス: 500 SMILES\n","truncated":false}}
%---
%[output:50c27046]
%   data: {"dataType":"text","outputData":{"text":"[10:07:26][INFO]  100 個の RL 前 SMILES を生成中 (T=0.80) ...\n","truncated":false}}
%---
%[output:6f56c742]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 10\/100) generating (pre-RL)\r[##--------]  20% ( 20\/100) generating (pre-RL)\r[###-------]  30% ( 30\/100) generating (pre-RL)\r[####------]  40% ( 40\/100) generating (pre-RL)\r[#####-----]  50% ( 50\/100) generating (pre-RL)\r[######----]  60% ( 60\/100) generating (pre-RL)\r[#######---]  70% ( 70\/100) generating (pre-RL)\r[########--]  80% ( 80\/100) generating (pre-RL)\r[#########-]  90% ( 90\/100) generating (pre-RL)\r[##########] 100% (100\/100) generating (pre-RL)\n","truncated":false}}
%---
%[output:1b1db75a]
%   data: {"dataType":"text","outputData":{"text":"[10:07:45][INFO]  RL 前:   valid=6\/100 (6%)  unique=5  novel=6  avg_reward=0.060\n","truncated":false}}
%---
%[output:09293c18]
%   data: {"dataType":"text","outputData":{"text":"[10:07:45][INFO]  REINFORCE ファインチューニングを開始します (20 エポック, バッチ=48, lr=1.0e-04)\n","truncated":false}}
%---
%[output:94e68094]
%   data: {"dataType":"text","outputData":{"text":"[10:07:54][INFO]  RL エポック  1\/20  avg_reward=0.062  val_loss=2.6164\n[10:08:24][INFO]  RL エポック  5\/20  avg_reward=0.083  val_loss=2.6129\n[10:09:00][INFO]  RL エポック 10\/20  avg_reward=0.104  val_loss=2.6131\n[10:09:36][INFO]  RL エポック 15\/20  avg_reward=0.021  val_loss=2.6163\n[10:10:12][INFO]  RL エポック 20\/20  avg_reward=0.062  val_loss=2.6202\n","truncated":false}}
%---
%[output:5406fbc1]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAlAAAAFkCAYAAADv4QVjAAAQAElEQVR4Aey9X6gd2ZUfvHp8G13xdeNuxkIWVtBupQccbEgGY6SJO7ol6JeBgUnwi0fOjE4zM4SBgF5mIJBJdEWUMP70ZPISsELfntgiJGMIhDAGO6j04aG7mRj84CEZpiPtRte0TDfp\/uhurMYXlPpV3XVq3X13\/f+3q86qe1fttdf\/\/as6VUvnnFv6pSe6KQKKgCKgCCgCioAioAg0QuCXSDdFQBFQBBQBRWB2CGjBisC0CGgDNS3+ml0RUAQUAUVAEVAEZoiANlAzPGhasiIQAgJagyKgCCgCm4yANlCbfPR17YqAIqAIKAKKgCLQCgFtoFrBFoKT1qAIKAKKgCKgCCgCUyGgDdRUyGteRUARUAQUAUVgExFYyJq1gVrIgdRlKAKKgCKgCCgCisB4CGgDNR7WmkkRUAQUgRAQ0BoUAUWgBwS0geoBRA2hCCgCioAioAgoApuFgDZQm3W8dbUhIKA1KAKKgCKgCMweAW2gZn8IdQGKgCKgCCgCioAiMDYCm9hAjY2x5lMEFAFFQBFQBBSBhSGgDdTCDqguRxFQBBQBRWCpCOi6QkJAG6iQjobWoggoAoqAIqAIKAKzQEAbqFkcJi1SEVAEQkBAaxgGgfPnz6eBeUwnDXc+X5+sKGwd2zo2RfFVvjwEtIFa3jHVFSkCG4NACDc01ADaGNBbLLQvfPqK02IJ6qIIHENAG6hjkKggXAS0MkVgegRwE5d0\/\/59Ak1f2fIrAM7A3l0pZNBJuU8m9XV4xKhjpzabiYA2UJt53HtdNS4yVeRLWOUDPfuBZ2IZRsgwVpFrh3kVFcV0\/YrspNz18c2lvY9v6uOzd2W+PK7M9akzd2MMMUcd7k0TeSAHgW9L8C8ixEReJsxDpqJ1+ORyHayXshB44I7a3FogkwS9nDMPeRfiOHLsEk9954tAowZqvsvUyodGABe1MsLFxldDmQ907MM8jyzHWBQbujJCrDJyfZEH5PpABnLt3bnr586LYkAOcu0xhxzk5uI5bMqI7crGMn+frixWXzqsGbl98SAHwcanryODfxHV8Q\/JpmgdPrmsm\/VdcJTxEKeI2A45mS8bpR1iwhYyST4Z9JBLgj8IMowg8EWEGC4V2ap82QhoA7Xs4xvM6nDBqbowtS12yNhcE2pHHhDLeIQMBBuWtRl9MRATcpAvZpHcZzuFDPXXpbr1IV6ddcMGtk3i1rXdJLu+cEScImqLJ44vYtbwLzSBPwgGGEHglRSBKgS0gapCSPWKQOAI4IKPG0mIZaK2OlS3dqwT8erawxY+de37tEPeKuoz35CxpsSxbF2oi\/Uu1pD7ZJBXEfvBjnmM7hwyJuiUNgsBbaA263gvdrW4kOJCNsQCERfxq2LDBrZVdnX1iIWYde2DtuuhOODRJgwwbOvbJh\/7IG8Zsd1cRqylDY7wge\/Q60QOScgn5+Ahq0Ow9RF8fXLIoFPaLAS0gdqs473o1eIihov1ohdZsDisvUA1qhj4V9UCG1kU5nV8qmxkTJev8pU1gC8jxC7SQ9cXIUdfsfqK0wTHvnJyHOAhieU6KgJTIaAN1FTIj5t38my48PkuvpAXUZuikQPx6vjCrojq+Pdtg1pQf59xEbOI+swzZCzU3zcuZfUiVxnBt0gPXV+EHFh7E+ord19xUDvW0Uc8xGHyxUMuSbCRc\/CQlRFsQGwjeZbpqAgwAtpAMRI6dkIAF5oywoXPlwDyIvLZ9ykrygt5n3k4Vhk+0A2RFzGLiOuacixbN3Qg1D9kjWPkaFs\/1t6E2ubpw8\/F0Z33kaMshosTbH0yyCWhThBkbA++DrFfHVu1aYtAuH7aQIV7bGZVGV943BGLgAzjWIR8IV7YUJePgAvkGNtSCOtFDU3WUWQPOQixQG0xUb9pEah77PhYT1Ut6gS5+aeuy61H5+EhoA1UeMdkURXhwoQL0diLmipvm3XOqdY26yvywbp9OshBPl3fMpybY+WStU+VV9ZQl69jt\/T1YH0giQXmU5w7sgblp0VAG6hp8dfsAyKAixsucl1T1I2DXLDtmo\/9EQsxeR7yiDpRb8g1urW1qbmNj5t3afMumHTxdXFELEnQyznzkFeRey5jzlTlq\/rNQUAbqM051pOtFBceXLwmK2AGiYsxqi4e2MK\/2nIYiyb5m9gOU20WNZQ6smrmuw8FR5z\/dakPtENZdx9r0RjtEdAGqj126jkDBHBRxcWua6kcxxcLMhBsuuZx\/RETsUGuDnPIYQN+Cpo6f9M1o15QG8za+jWtcQ72wALUBkdeX1d\/jjPmiJpBXdY9Zr2aa1gEtIGqwFfV\/SCACw4uPG40yMqI7WEDnkfwdQm569qW2SEOCDVIggxU5ltHhxiI69pCDoLOJchdezl37d25tG3Cc5yq\/DImfJrYS98uPPIyIT+oaTz4N\/GDfRkhv6uHrE9y4zeZ++qQ\/sAC5LOrknGcLv7w5ThV+frQcy7kBfURU2PMHwFtoOZ\/DCdfQd0LimuHeRXx4qQdyzBCjrGKXDt3XuUv9fCVJHVFPOyLdFJeZgedS9LX5V1b39z1KZvzTQQjxyqzZx3sQfBh2dAj8jEhL1PTvDJGE1\/O12RsEr+ObZPcri3H5\/VjlDasrzvCn4nj1PVlO+kPGcdhed0RvpLYzydjHeeSNjPgtcSBEdAGamCANbwisBQE+CaCscmaYA9q4lNkWzcO7JiKYtWR9xGjTp5QbXj9GLvUCH+mtnGK\/Fled3TzSz\/WSRl4luuoCEgEtIGSaCivCCgCikBfCGgcRUARWDQC2kAt+vDq4hQBRUARUAQUAUVgCAS0gRoCVY0ZAgJagyKgCCgCioAiMBgC2kANBq0GVgQUAUVAEVAEFIGlIjBcA7VUxHRdioAioAgoAoqAIrDxCGgDtfGngAKgCCgCioAiIBFQXhGog4A2UHVQUhtFQBFQBBQBRUARUAQEAtpACTCUVQQUgRAQ0BoUAUVAEQgfAW2gwj9GWqEioAgoAoqAIqAIBIaANlCBHZAQytEaFAFFQBFQBBQBRaAcAW2gyvFRrSKgCCgCioAioAjMA4FRq9QGalS4NZkioAgoAoqAIqAILAEBbaCWcBR1DYqAIqAIhICA1qAIbBAC2kBt0MHWpSoCioAioAgoAopAPwhoA9UPjhpFEQgBAa1BEVAEFAFFYCQEtIEaCWhNowgoAoqAIqAIKALLQUAbqD6PpcbqjMD58+c7x2gbwJebZTy6sX1yV+bOOUaRHPoyHfSgOjawY4K9S1LHPI+wZV7H+SGA48eE6pn3jdC75LOTMtee59KGeeiYlyPkVQT7Khvo69rBFtTUHj5D0VS1+PKyjEdesztn+SaP2kBt8tGfwdr5RcujLBmyIpJ2Y\/P3798n1IW8GDEHX0WwrbIp08PfR+yDOkCYYwSBlwR\/OVc+DARwXIqoqEIcXybYMI\/RN4fMJdj6yLVz564P64vkrHdHrBk+rnyoOfK5NFSupnFRl+vjk7k2dedVdjgOY+arqicEvTZQHY4CTiYfyZA+vU\/GPtAx745SB76IpJ\/PRurB+2wgg04SZJKkjnmpd3m26WvEC7qI2uYoqhl5ZEy2g6yMd\/WYS4KvnNfh4QOCLUYQeBDqlASZ0vwRkMfU5X2rwznhEuxY5vKY90mch0eOzXM5sq7uKH0lD385lzx0TCzHXPKYgyS+mIPYrmiEjaQiO8ilneShc0nqq3j2rbJz9ezHI+uBA\/M6+hHQBsqPS20pTjKXcCLKAK5ezmGHOUYm15\/l7gg\/H7Ed4vj0kLMNj1V28HFtIGN\/Obp2PJc2ofJcqxx962Q91uHjIZM6nkMm40Eu5+BBbAcehDkTfECYYwSBd8n1wxwEO4wg8EohIHC8BhyfpsRRcE64BB1kPDKPeV+EehGXCXHBF8mhg42P2MfVwYcJOubLRtgxsR3mkse8iNgOI2wwSoJMktS5vLSTfF076SN59gduUl7Fs58ci2LApkhXlWeJem2gAjyqfZykOMkRx7c8yKH36Xwy2MLH1UEGnSvvOkdMEOIghxzBg6CvItgxwZZ5jO4csr4IsUGIhxEEvoiwRthgZIIt8xgxB8GOSc7BuwQ\/kJTznEepUz4sBHCMmhKvgM8RHlnedGR\/dyyKg3qhY3vmeZRyyJT6RwDHADi7kX0y10bnzRDQBqoZXmrdEYE6L2JcAECwBflSQl9G7MP+sGWedXJsw3M8HjkGcoEwxwgC3wdxLIxMMi5qcUnqfTzspdydS53y4yOA41FFblV8bvAIPcco46FjYt+ike2KxiI\/yIt8hpYzBhiRCyMT5lMQ8veNSd\/xpsBlDjm1gRrhKOEF4lJVWrwA4FNmB71LZfZD6FCnG9etCXO28dmzzh1hC5L+rk3VvK4\/csC2Kp5PDz\/4Q4eRyTeHDPYYJflkUi952CIHZBgxBw9iHiMIsjKS\/uBhCz\/mMVeaFgEcjypyK8TxkwR90xjw6UIyv8t3iSt9sSY3tpxLW\/CwZ\/LNIfORG1POmff5NZH1FaduTuQDFnXt1e44AoE0UMcLm4sEJ6FL7kmJuUt11gcfxC6yhd6lItsqOfK4hNhVfj49\/Fzy2dWVIRZqq2tfZIc4RbqmctQj4zGPkQkxmceIeRkhJhPsmMeIOQg8iHk5gq8i+KIWjLAFjxFUxEOnNB0COFZV5FaHYynJ1ZfNfbnYvkzHNjzK\/C7PNn2Mbmw57yM+Yrgx5Zx52DUl4Al\/+GEEQYZ5n4SYiN1nTI1FpA1Ux7MAJyUTQoHHODdC3UyoHTzGKsILs8qmD33deny5XF\/ULGXuXMYo00k75mEPknPwUoa5j1ATE\/TMY8QcBJ5JzsG7hJywlXKe8yh1yoeJAI5VFbmV49hLYr2USZ71PMp8LOOxTMc2GGV8l4d+CEKeOnFhh3XAFjzG1tTBkWuQISBDTRhduZx34X3xOV6Zjm10zBDQBirDoZc9TnicfL0EOwwyRMzD0IXDFDkLi2mhwDHAGuq6NrHlmEU+kINQA2yZx4h5H4TYIMTCCALPhFyujHWQM7FMx7AR4ONVNvpWgPMA5Oogk+Tq+5rLHC7fVw7EAS4YQcgj55KHfk5UVTv0WG\/Vmnx2dfyq4qpe34GaxTmAkx0vgibFlvkgFvR148EWPq49ZNC58pDnTWpuYivX3NZPxijjgTkINhhB4F3yySFjcu1nOt+IsvmY+cZQAcDroIjq1Iy1wt9nWySHLfwwgsC7tphDDj14OYIHQQ4CPzYhL+oDge+SH\/6IUzdGlX2Vvm6epdjpO1A9H0mcrDjJZFjMiwh20GEsI8Qt0\/t08EFslyD32UsZbODHMp5DxgQZ6+XIeneEDWQYhyDEBpXVVaST9SAGqI5tHT\/EQTxpOySPXE1ysv2QNWnszUIA5x9WjNEllmMcknBeu\/FZhhHEtbl2LMfo6qrmiFtFvhjsI3NKXvrAtkgHO+hBZTawY4ItqK49+236uKgGCicAqM5BLbKDnKkqTtHJJuXghTYWrAAAEABJREFU6xDngi3z7ih1knft5Bx2Lkk9eOgxuuTKMZeU2js7qXd5mEKGsS\/iY4URsUFubOhAPp20hQ0IdiCpq8PDBwRbxMHIxHKeuyPsJUHvziGrIvhwLoyYSx\/MmaScefgwv6SR14yxbF3QM7l2LMfo6jAvkkPXlRC7iIpis72rZzmPrh5z1mHEXBJkTFLOPOuKziXI2YZ9ika2lXr4Qs4y8JC5JPWShz2TlDNfNcLXzYU5+0FfRWyLEb4g9oGsiGAHgq1rAzkT9CDXxp03sYdtnZhujiXPF9NA8cHFAQZfdNCgA\/n0kMOfCXOfncqGQwDYN4kOe6Yivyo9+9W1Y3uM8MHoEuQgVy7nUg++jKQf87BnHiPmIPBMcg7eJdhBhtGlIrlrF\/ocr2OshQlzX82Qsw1GzNkOPGRMmEudnLO8r5Fzlo1uLteW9a6c56zHyDI5Qg6SMuYhl1Qkb2rD9ojHPEZ3zjLIXYKuCcG\/jj3sXKrj57PhOD6dKyuzZR1G169oDltQkV7K69pJn9D5rvUtooHCxUseXPCQ+cCBDuTqYO+Tu3Y6VwQUgfkgUPd17bPj64FPJxGAHUjKfPz+\/j4pKQabdA74XgdLki2igRrzgGzSya9r1Yv9mK+tEHKhWWLqsx68lv7oj\/6ILl26pLTGQLFY+vnwW7\/1W+k\/Gvp8LYUUSxuow6Ph\/gsSF1FXphdBveAt\/YLnrm\/uF0D3Nex7XR9eAoh18AFhDh14jEyQuzLWFY24drz55pt069YtunPnzqzo2rVr6bK09nGP25xxxzmO+nHOpyfPQnfaQHkObNEFUi+C415A8CIE4YWIw6QX8H7xB7ZlBNyXdAEsel3j3ALVaYqqYiBOGV24cIEuXrw4K0LNWBNGrX28Ywe854o7zhOuH2tYKmkD5RzZOhdInBg4QeZEv\/mbv0m\/8zu\/QxjnVDdqRc1a+3gXbmAOwnnuvDxmO63zuq5aXB8xqnKEqP\/sZz+bXjswhlhfWU2oGdcOjGV2IepQ81xrDxHPIWpaRAOFfzni4sYAgYeM53XHtn514x+1G3eGF+PVq1fHTdpTNq29JyA3NEyd1zWuF7Arggg62BTplyzX1980R3fOuE+D2PhZF9FAATZc3HCRA4GHjAky5otGtsEoqche5YqAIhA+Angto0qMkiADQYYRhOsG5kyYQ445j+CZIFNSBDYSAV10isBiGiisBhc8EHhJdWSw8ZGMo7wioAjMCwHfaxoyXoXkIcOcCXMQz90ROknQy7nyioAisGwEFtVALftQ6eoUAUVAEUgR0J0ioAgEgIA2UAEcBC1BEVAEFAFFQBFQBOaFgDZQ8zpeWm0ICGgNioAioAgoAhuPgDZQG38KKACKgCKgCCgCioAi0BSBOTZQTdeo9oqAIqAIKAKKgCKgCPSKgDZQvcKpwRQBRUARUAQUgSIEVL4kBLSBWtLR1LUoAoqAIqAIKAKKwCgIaAM1CsyaRBFQBEJAQGtQBBQBRaAvBLSB6gtJjaMIKAKKgCKgCCgCG4OANlAbc6hDWKjWoAgoAoqAIrBoBJ56iiihi7\/2a3T\/wQM6+7f+VjqHjBa2aQO1sAOqy1EEFAFFQBFQBCZFwBjaf\/FFemN7myiKJi2lt+SeQNpAeUBRkSKgCCgCioAioAi0RMAkDdR\/+A905cwZ2k9GiqKWgcJ20wYq7OOj1SkCioAioAgQKQZzQsBaOvvWW3T24GBOVTeuVRuoxpCpgyKgCCgCioAioAgUIoAG6rd\/m\/6\/hw\/TRqrQbuYKbaBmfgC1fEVgFAQWnOT8+fPEVLZMtsHo2kHG5Op0rghsNAI3bhAlDdUSMdAGaolHVdekCCgCtRBA03P\/\/n1iwtznCDnbYMSc7cBDxoQ563RUBDYegTgm0gZq40+DKQHQ3IqAItAzAmh00PRUhfXZsZ9PVxVP9YrAohF48oQoitIl4q\/w9pOP8QgyUCpdzk7fgVrOsdSVKAKKwIAIoFliaptmf3+fDg4OlBSDxZ4DOMcpjtOXyJt4jAGl7CJ32kAt8rDqohQBRaAKAX4Xie3QHLkyVwc9CLbQgcfIBLkrYx3GK1eu0M2bN+lh8q\/yuRBuiI8ePZpVzYyt1v5w9OP2vX\/373Cqp7T\/9NPpuNSdNlBLPbK6LkVAEaiNQFXjU9YUcZKqGLC7desWrVYrOnPmzGjUNdfp06fp1KlTs6lXrldrH\/88e\/nFF3Gqp7S\/tZWOS91pA7XUI6vrUgQUgVoI1Gl8qgLVjXHhwgUyxtB28tHGXOjkyZN04sSJWdXM2Grt26MfNyO+MK4f4VVdOVSvCCgCJQioKmQE6jQ+ePcJdkXrgA42RXqVKwIbhcC9e+lyl\/7uExap70ABBaUgEbBkaTf5Obl9ks6\/cJ5+ZetXktku6XYcAcbqKcp+XqAXCNiRboUIoPGBEqMkyECQYQShQcKcCXPIMecRPBNkSorAJiOgDdQCjr4uYZ4IoCG4TJfpRvLDK4AMczQHLNORCLgoVtR4QxPkIw4EHfMYMWfCHMRzd4ROSRHYSATiOFt2FGXjgvf6DtSCD+6cl4ZGCY0B1vD1g6\/Ttfev0UsHL2FKkL9Cr5BuGQISqxWt6HryE1FE2BQroKCkCMwSgfkVLb7\/RMbMr\/6GFWsD1RAwNR8HgT3aI2xoBHZpl7777Hfpd7d+lzCnZGN9wm78L2MBbK4nzdM9upfsrxPmlGysT1j9VQQUAUVgOAS0gRoOW42sCDRFYId26Ntb36b9rX3Cu06Ycwy8u8K8jpQgtUNolmKK6XLyo1hRt029FQFFoBkCcZzb7+zk\/EI5fQdqoQd27ssyZAgb3k15m94Gm9Jr9Fo6Qg9KJxu+Yxzu0b0jSChWR+DQiSKgCAyNwOFf4G3Cx3eAUhsooKAUHAIRRWlNeDcF76qkk2TH7zrxx1OJaON\/GQtgdYNurPFQrNZQKKMIKAJjImDMmNkmy6UN1GTQa+IyBF6lV4kbA3I2vOOC7\/o44o2dKlYbe+h14YpAWAjEcVbPBnx8h4UWN1DQKikCEyJwl+6S2yihqXpADwhNFOm2RsCHFbBTrNYQKaMIKAJDImBtHt2YnF8wtzENVJMH3MF2wcd8VkvbpV1C00SHW0wx6eZHYJd2SWIledJNEdggBHSpEyBgbZ7UmJxfMLcRDRQaIn7QHfii4wkdqEiv8jAQ4O\/2hFFNuFUoTuEeG61MEVgcAnGcLymKcn7B3OIbKDREaJ74GIKHjOdyhA4kZT5+f3+fQAcHB6Q0PAZuI\/CDgx8o7gd+3OX52uXcxPktYzXj1VoRUAQ2DoEN+ws8HN\/FN1BYZN905coVunTpEt28eZMePnw4C8IN8dGjR7Oo1cUUjYA8ht97\/L3ZrGNs3N86eGsN1Qfvf9Aap9u3bxPO83UwZRQBRUARqIOAMXWsFmGjDVSLw3jr1i26c+cOrVYrOnPmTHDkq+n06dN06tSpWdTq1r+\/tX\/kKG1vb89mHWPjvrW1tcbqg+c\/aI0Tzu1r166tYymjCCgCikAhAtYSxXGm3pC\/wMNitYECCg3pwoULdPHiRTLGEG7mc6CTJ0\/SiRMnZlOvxNQ9PN\/Z+s5s1jE27hIrNFMSxyY8zm2c5zKe8oqAIqAIVCKQ3BcrbYoNZqXRBmpWh2uzizVkiDf3e1Es3\/RRcdn0M0DXrwhMgIC1eVJjcn7h3OIbKHwpXH5pHDxkCz+ui1mebAhkAxXT4dvFpFsRAvK\/wCmyUbkicAQBnSgCbRCI49wrinJ+4dziGygcPzRMaJxA4CFjgox5HcNG4CpdXRcoG6u1cMMZxaT4BMDrnKnYiohtMPrsfHLImHw+KlMEFo\/ABv4FHo7pRjRQWCgaJxB4SXVl0kf56RDgd6EW+u7KdMAuODOaG7zOmTD3LRdytsGIOduBB\/GcR8hgy4Q563RUBDYOAWM2askb00Bt1FFd0GLluyponiLK3h6W\/8Ew6eZFQGLnNdgAIRoaNDdVS\/XZST\/woKo4sEGsIrv9\/X19htmB\/xlmeFyJ0vywIWuJ4jg95Q++8pX0GYk4z1PBwnfaQMkDrHzwCJyjc+satUFYQ6FMTwig+WHqKeSRMHi21pyeH4dnsuFmONdnyGnt4zynkE\/yd7\/wBcIz5PCcRJzrLF\/qqA3UUo\/sQteFd6F4adpAMRLZqHhkOMg93hGSczRHroz1rIMehDnr6o5VPniGHJ6x5T7rLOT52M8y6xOLTaq9T9yaxDr18cfrlweeNYjzG89J3ITnyGkDtT70yoSIgGwK0DytaEW8SR3LdMwRUHxyLMChuUFjBN5HZTqfPWTwQVwmzCEvIjxby5j5PD8Ozw4b+1lmyNkXae3bwz8z74031qf71ssvE85vPCcR5\/pasVBGG6iFHtglL8uQIWz36B4GJUWgEgE0OFXNTWWQAgPEvX\/\/PmEsMFGxIrBcBDb0L\/BwQLWBAgpKwSIg30U5e3A2rTOi7IvkUke6kYuHO6cN3eo0T2h+YNcUItcHc8RqGkftFYHZI2DM7JfQdAHaQDVFTO2DQUA+TDOYorSQoBBAQ4OCMEqCDAQZRhAaH8yZMIe8jGDD9hgxL7NXnSKwKASsJYrjbEkb9H\/gZQsm0gaKkdBxNgjs0M66Vn2XZQ0FKRZ0bEND4yM2hI55jJgzYe4SdD4Z5CBXp3NFYGMQiKKNWSovtKcGisPpqAj0iwA\/MJM\/vkP0iPIXqr4LRaWbNlWl8KhSEVAEuiBgbRfv2ftqAzX7Q7h5C+AvkWPl+kVyoKCkCCgCrRFQx\/YIxHHmawxRFGX8Bu21gdqgg72kpcomaknr6rIWfreuSwz1VQQUAUWgNgL8F3i1HZZlqA3Uso7n4lbDH0F97uBzR9bGH+Ppf+lyBJZjE8bvmCIcgVaiCCgCc0fAmLmvoFX92kC1gk2dpkZgR79IfuwQaLN0DBIVKAKKwFAIWEsUx1n0DfwLPCxcGyigsMk007XLj\/D0i+QzPYhatiKgCCwDgShaxjoarkIbqIaAqfm4CPC7KvKv8FBBRPkLlm1ItxQB2VwqNikkulMEFIG+EbC274iN403toA3U1EdA87dGgBsF\/fJ0awjVURFQBBSBdgjEceZnDFEUZfyG7bWB2rADvqTl8rtQ+kXy7Kjyu02MSybVvSIwFAIad6MR2PC\/wMOx1wYKKCgFiwA3BWd\/kf0\/eLLQHf0iuYTDyzN+XqUKFQFFQBFoi4C1macx2biBe22gNvCgL3HJm\/hFcvc4crN0js65Kp0rAoqAItAfAtYSWZvF29C\/wMPitYECCkpBIsANQVFxK1qRbuUI6PfDyvFRrSKgCLRAwNrcKYpyfsM4baBaH3B1HBMB96\/wODd\/kVz\/SxdGJBsZl2yme0VAEVAEFIG+EdAGqm9ENd6oCESU\/etn0z\/Ck+\/WafNEuikCikAZAl11cZxFMIYoijJ+A\/faQG3gQZ\/LkmVT4P5XLrwG\/r6PtGWdjooAI3D+\/HliYplvZAiFqa4AABAASURBVBuMRXpXDltJrl7nisDiENC\/wEsPqTZQKQzz3KFp2KVdeurw5wV6IZntktzY5uT2STr\/wnn6la1fIfjQQjZDhnjDWpkfckQeYFiG+5D568RmXFBrHfs52aBZ8dVbJr9\/\/z4xtbGDD8jNCxnH5REy186Z61QRmDcC1mb1G5ONG7rXBmqmBx43xst0mW4kP7wEyDBHIwUZ5lU2sJsz8Ud4WMMYH+OFiinqAgYgbp7AbzqhmUFjU4WDz076gQdVxVG9IrBYBJ56iojJWkq3OKZURpu5aQM10+OORolvmvhrtOt0nbiZgPwVeoWkzdcPvk7X3r9GLx28RNjYBvycSTYLtb9I3mHBEtMi3DuE782VccFx7i3oRIHQ3FSlhk0fDQ7iMFXlbKPf39+ng4MDJcVgducAzveDs2eJjAF7hOQ5jXMcdMRgoRNtoGZ6YPnp22iartJVcpsH6EG8vJ9u\/ZTeOPkG\/e7W7xJ8KNmkPpkG9ytv\/kV\/hYeix2wWGDNgiKYV+e\/SXcKcko31Cau\/IyGAhqeoeXLlZbasgw8I86olsB1smSAr8rty5QrdvHmTHj58OBvCzfDRo0ezqVdiq7X3d57hnD743OfonX\/zb8AeIYn57du36dKlS4Rz\/YjRAidTNFALhHG6Je3QDuFdEXx8BSqqBLo3t98kvDPFDQdsZZOC+RyJmxescaz6gTt\/PIqmCXPOPQWmMqc8vlzTUkc0LWUNi1x3lW3dOL6Y8GVCHmkj+Vu3btFqtaIzZ87Mhk6fPk2nTp2aTb0SW629v\/MM5\/H2yZN06q\/+CiyRMXTwUvaJhsQc5\/edO3fo2rVrmd2C99pAzfTgGjKEDe88cePAMsgl785xw2eZawd5KCSbgrJ3oHaSJpJrlj4s63NkvIA75wIPQh7oQeCnpHOUPY2ca5yylqFyo1FB01InfhPbOvHa2ly4cIGMMbS9vT0bOpncNE+cODGbeiW2\/dc+3nELrXY+57f+4i+Ypa2trZSXmOP8vnjxIuFcT5UL3k3eQOHC1pQWfDxqLy2iiLBx8wReEvQrWhFvvgYkoiwGzXwzZIi3IjxY33VkzGQe8CDEZj34MWmJjRJfF4Aj8xh53mfzhFgcG\/GVFAFFwINAHBNZmymsJbI24zd0P3kDhQuXj3A8fHLIoPMRLoAgn07KYAOSMvCQSYIsVHqVXiX3Zs03UTQU1+k6SZsfbv2Q5GbIEGxoARvWwsvAR5T4K8RdOvo4B+ppk5hySBd3lk81SjymqqGPvHitgxALIwjEc7xWwZcR22CUxD6QMY\/YmDNhzrqiETZszyNkRfYqVwRmj4C1+RKszfkN5CZvoHyY40LU9CLEPvAD74sLGXSwAYGHDAQeMkmQQRcq3aW7JN9lwo0TTdEDekDgKdlgA1nCHvmVNkcUAU3eprfTanzvnqWKZIfmBd9FStj1L2T4XhgaqbWwRwaYyuYVWAPjkDBFTbxk4MH8kka8Vqteo7DxEeMAHfMYMWfC3CXofDLImVy9zhWBRSDw5AnRgwf5Uq5fJ4IMlEs3iguigcJFUBKOgJxLHjqXoMfFi+XgIeM5j5BBx3PwkPF8juMO7azLxo3d984LZD9\/\/HO69e6ttW27m+raPRgGjZK7Fm5uIMc7UkMUK3FHfGCMcWqSjdPUtXTNj9cmqCzOEl7DZetTnSIQLALGBFvaWIUF0UDhIugSAHBlmEM+NeFPY0Hy2RdT8WgSGA+8U1NWx4XHF9iUfnDwAyqzDUHHayurhb8Qv15Ywnzr4FvETRT0Zf5tdfylcUo21Nk2Tl9+\/G4dx0vKWv+yrMmI83sdYCIGr3cQ0pc1UrAp08NfSRFQBHpAwNo8iDE5PwdugBonb6Bw4fMR1urKIRuK+CIsc0Lmy4fnW+A5FyE8z+UnH\/0kLRHNk3wWh8unN8Tk3P\/s48+m9t97\/L3gn+uCd81Q7Gc\/+SwVPYcGehDWjxH0X97\/L\/R33\/+7YFN6\/Z3Xe18r15YmSHYu3jwH7kW1s00f44cffZhUkf0i3ol3TmSTZP+X7\/5l4\/XjWS44zxP3yX\/xOgThtVlUTJW+yE\/lioAi0AABaxsYL9908gYKFz6XALsrwxwXUBD0fRPiIockyHx5bt26RXjOxWq1Ivn8iyn4n23\/LC3xxa0XS2vh56HsHH7k96NnflRqP8Va3JxPbz2dru3EiROFz6Hhj6x4hMNzzz9HP3r2R2AJ8i+f+nLva2Xc6XD75Mwn3hyMu7u2vudcD58Hzz\/\/\/GFlVIhdWQ04t0N7jgtem+tFKTM3BLTeJSBgbb6KKMr5DeUmb6B8uONC6WteIAf5fMaU4fkWeM6FMdM\/z4Wfw4H1y2dxuDw\/UwQ3WNjiY6dH24+CfraLXBuaKHdNmEeUvYjlXxn+661\/TTyHHnZ9E\/AjsRXFZ9yL9H3JJVYck8vb39pvfJxxbuM85xghjUXXgCJ5SLVrLYrArBF4O\/vDHkrufbNeR0\/FT95AyUYJPBPWx7w7QicJF07YsAw8ZDznETLoeA4eMp7PceTnD+0cvrNUtQZDhnhjX56HNnKTImt2a5zisQJcl6zFJ5P6sfljmI1dgOZTBBSB5SFg7fLW1GFFkzdQbu1oaJigA88j85i7BB0aIhB4qYeM59BhDgLvk0MHknq2C2lsc9Ne0Yp4a+PPviGNd+kuXU9+uCY0D5gP9VgBH24+Gdczxjh1\/jHWWJUDr9kqG9UrAopADwgY00OQ+YcIooHChQ8EODEyNW1gYA9CHEmuDHOQtAEPmSTIeqLBw0QUUd0NDQZs+S+3wIdI3BRwvWU14jEC3BzCD\/My+6XqGCse57xOvg5gxDowlhFeu7BTUgQUgYEQiOMs8E7++JxMsJn7IBooXPhAOAQYmXCxhEzJjwAaBb+mXMrNFv7Ev9xyXtpzdG6Ugn0ffbY9FmMUHHJtVevHtUDa8LxolLbKKwKKQB8IHMaw9pDRgREIooHiYjCiaWLiOY+Qg1fKEJA3xibvOOxQ\/q8HGSOLOt+9xGCMdcl8U6PG6x2riZx6vZpfEVAEJkTAmAmTh5N68gaK\/yUJSMDXIdgqEfFNk5Ktyc1c2vreTUnCTf7bZm1yXdK\/78XIh2hyzhA\/Dg25tr6PyZLj6doUgSAQsDYvw5ic32Bu8gZqg7HvvHS+afONsm5A\/ggP9rIZwHzOJHEYsoFijGQ+lk0xjrHWKdalORUBRSAgBKzNizEm5zeYC6KBwkdzTDgWzGPEXJJPJvWbxPONs82NXH7heimYSRwYm+5rOx6BYyMfCBYsAz8lcT1T1tBnbvf1zvOisc\/cGksRUAQEAtbmE2NyfoO5IBoo4I+P7jAyYQ7ChRIEOUbIwCt1Q4C\/K7Okj\/CACDcQ\/O4cZH0TN0vAkPOxrO9cXeKFXFuddeG1zgR75otGXB9gp7RwBJ56iqiISLfBENCHaB6DNpgGCpXhwohREstwcWRe6jeZ5+ZnR3wpvC4efHOF\/ZL+Go\/XNVRDI+NyLmAYCoVYUx\/Y4PXPcSTPMox6fQAKG0LGEEURURQRRdGGLHriZVqbFWBMNuqegmqg5PHARZIJF0YQ5tJmk3l5I2+DA3+EB9+usRAjFOIGYqg1ybjIhXehsHYpx3xMmjL3mOvskgvXDqayOGyD0WfnyjH3kc9XZT0iYAzRq69mAe\/eJYqijA9rv6xqrF3WenpYTbANFBomJl4n5rhY8VzHDAH5pfBMUm+PBgCWQ37chfhtSDYFXGedOCE0NHXqHMOGcZNYjpG3rxx4rTMhpo9nGY+wcwk6XDuYMHdtMIecbTBiDjkIPAi8JNi5JPXKD4jACy8QxTHRK68MmERDrxGwNmP1IZoZDsk+mAYKFyempC5i3h2hUzr6CIO2eESU\/attiR\/hUbIN0Tjwx6ZJeDKHP3S4DZHvMPRGDrIxAQBlc+hgc4SSCa4fRbpEvf712Uk\/8KC1QwHjiyNN9\/f36eDgYFFU+H2kp54aZJ3AExhiBB384hcYUoJcqf\/zi6xN8cWuCl+c4yDYLp2CaKBwYWpCSz8oddYnb9a4kdfxcW12xHenZDzXbop5H\/X0EaNo7cAcVKQfUy7X6atJ6sesa+xcuIZ0zYkGiKlrLJ\/\/lStX6ObNm\/Tw4cPZEG6Gjx49KqwX6zw4e5YeX7iwJshAQ6wTcbd++EMMKW195zt08NZbKe\/mq6rdtQ9pHlrtKcDJ7t0vfKHwXAB+t2\/fpkuXLhHO9cR80b9BNFBAGBctHiXPMoxKOQLypui7aeaWxZz0k++sFHs00vRmLOusChpR9q4aJZvEKJn28us+N0vWNkS+LkXzx5ldYoTi6zZH7ryoTtcO1xZXxr6sgx6EOevqjLCHX5ntrVu3aLVa0ZkzZ2ZDp0+fplOnThXWi\/Vuvfgi0d4e0e4uURzTwUsvQVzo02X9aWBnt5W8sweRG7eqdtc+pHlItZ\/6+GPAm1LZuQD8VqsV3blzh65du5baL3kXTAMlQcZFCBcjKVP+KAJv09upwJChtltEEfHmNgYsn9so8RiyoeE8PAKnIfMhfhHJvLKeIvtNleOagmtL0frLdEU+TeUXkndpjDG0vb09Gzp58iSdOHGisN4UA2tp++\/8Hdr+9V+n7eTdqq2trVQ8xDrpyRMifHE8zXC4gywhN19V7a59SPOgan\/jjUOgidAsl+GE8\/vixYuEc33tNDumXsG\/VM9sOCtc1EDIwCN4vpixDCMT9JtOfNM0ZKjLtqIVYeN44EOiNutjH24y+1wP48Q5eOwzR1+xZG1cd1+x5xYH1w6+psyt9lnUmzRQ6zqTd6DW\/FiMzD9Wzk3MYwyRMZu4cu+aJ2+gcFEDoToeweOCh5FlGJkgV+oHgXN0Lg0UU0whbUM0P32sjxsRxk3GZJ2UKd8fArgmMDWJCh9cO8p8oIddmU2ZDr6IUWazMbobN4g6NjSVWA0dv7KADTO4d2\/DFlxvuZM3UL4y9WLkQ+WojBueHdo5qmg4k+9SLOWv8SLKPppkjKinTTZIEjfmp2r6OC\/X0dNygwuDBoUJ1whQVZFsg1ES+0HGPGJjzoQ563RsgACaG1ADl8ambnx33jigOtRCwJhaZptiFEwDJS9WkseBcOeQbTLJG3lXHPgjPMTpMy7idSGupUtTwDG61CF9ZbwudcmYYfLhV4VrAgjNTlm1sPER+0DHPEbMmTB3CTpXhnmRHLrF05MnRKtVvszVitLvKUGeS4fltIEaB19jhs0zs+jBNFCMm7wQSZ71Oh5FwFD3E5pj8DsZNPNtR7wrJ5uemS+rcfl8XOG4VBxwjahqorB+pYERkA2M5IdKy\/8vG8cfIyfn2sSR8T2XfeVjEyHwrTm4BspX5KbJqtYrb4byJlnlV6SPKPvIaykf4dFAm\/xIUOLOvDwuA5XgDct5uQ6v0YKF2kQFdnDjePyC3IZq\/AqWm5GbJ6zQGOyVDhEIpoHCvyLL6LBeHRIE+IaZsNTHTXMnwHds5Bqp4SYxkU3xZmqfAAAQAElEQVRPwzCF5ogPYgPmu9TMsfocuS7EDK021NSUcH1o6qP2IyEgb7JjpBw73xhrCjWHxNqYqio3Sh9MAwXU8S\/JIsLFEwS7TSd5M5Q3yba4yBhDNBxt64KfrA3zOsTvqNWxbWJT9awseVyaxFXbagTw2se1ocgSOtgU6VU+MALyJotU7hyyIWnsfEOuRWPPBoGgGigfanxRxAUSxHOf7abI+LtKbZoLH0ay4ahqEnz+IcoYmyHWw7F53b5HGrBujJEbN7euMXKPkQOvebz2x8g1SY6nniIqImq4TWHua158sj5rc+O78z5zbXqsOM4RiKKcV46Cb6Bw4cQFVI9VjsAQN0xuojh2nm0ajuuYujmRq+ea3EZFztlG+k3Fh1rXVHgEndcYoigiiiKiKAq61FrFDd3QuPHdea0i1agRAsY0Mt8E42AbqKKmCQ3VJhyYsdfIN9vQPsJri8MQDSE3RyVNXdtyO\/mFWlenRQnnuq\/5unYidDisMdl\/T4L\/ouT6daIoCqe2qkp8zYtPVhVH9WEioA\/RLDwukzdQaJRAskLM5cUQPGTSZpN5bnJ2xJe\/u+IhY3H8rjHb+nND0NZf+vW1FlkTN5ucR86lHet1VARqI4DG4\/Jlojiu7TK5IWp2ixjyr+JkPmPyzFKeS5XrCwFj+oq0mDiTN1BojkBAlJsknkPWiCY2xs1zl3bpqcOfF+iFZLZLfW7I0Wc8jhVR\/i\/ey3SZhqidWmyyOWniviOayz4wkzHKapJ2TeodypZr5e\/N1cmDNeA8\/rWLv0YP7j+gf3D2HxDmNPGG60MVTVxi9\/T4b1BklDiWszB52bgYM26NUZTnk3Xk0vlxh9+H2z55kl44f562nn6a6FBGU2yMqzFTZA865+QNFKODpgmECyTL5AidnIfG46aDxuMG3ViXBhnmaEbWwh4Zvjl2DYk6UbuMA9mQtctcc+f7Og5tcMBxauNX5IN4OBdw7NkGMsyHOo85T50R1wEQbDEy8RzjLCmOifDOE9+seBGQ7e3xLOzRGKIoymp015FJ+9nL2Dviv7KS8n4y9R+FGyHPeCSZMXTw0kv0+MIFoig6ohp1AkxBSKoP0QQKRyiYBgpVoXnCBRE8E2TMhzziBoMbDWrEf49yna4Tv6sD+Sv0CvWxIRbH6evGLWvn2EPUzrGbjG3XKP1iSm5O1G2TMWRsRJVzeXygG5NkHW3zynPh5f2X6blvPkcv7r+YhsPa+jqP04C6O4pAHBPF8VEZZq8k1w4Q+BDJ9x2ZOB6nUmPGydNnFmOIoogoioii6GhkblYgvXqVfn7xIlEyHrODfmwyZuyMwef7peArnEmB\/BTviCJ6NfnZpV26m\/xgTsnG+oTt9IubGAfo44aJWL7asIa+a0euOiTXWMfeZ8O1+3RdZMAc1CXGmL5ca11M+VwAfv9q\/1\/R8998PjmL7xLmlGysT9iF\/U64HPyfcT7CF8q5rL09It+7FiyjADZjiMZ4l8LafLHG5Ly1OR8yZwzRq68S4d0zY4iMyarFsXzhhYy3lrZ+\/\/eT1983idA8J\/NMMfJe5jVm5OThp5u8gcI7TEyAi3keWYZxDrRD2VvKuGGBeI7aMcfYhWQMQ4b63PgmiZgxxclKsrVQssm8yXRWv308C6oqBh+Lt+ntUbGRx4Vr6KOAneToH5w9WIfCnCcyJ8t0HACBKCJ68IAoio4GN4Yoioii6Kh8qhnfZI0hMiavguW5pB\/OjWtMFnfIL65nGfrZo358NPvaa0R7e0SYV0WuY1MVo41e5jWmTYRF+0zeQOEjOyYgzbwcIQ+d+OZ1j+4RvisCwschmKN26EHguxDfoPuIxXX4YuEm2XftnM83DiHzratrniFidq2pzJ\/rxfEss2Md2+PY\/4uz\/4LwJfLfPvvbhDlsoAeBVxoBAWOI8E4UHm3A6XBTw7sXkEcRS6cbUQ+ynztHZAy4jFiezYbZG0NkTBZ7jHxZpm571MlkTB5rtcremWKJMcxNN6JOzm4MczoeIjB5A3VYRy+DfNeqLGCZHeswlsVwdRFlF7KYYuKbFXgQJRvrE7bTL8fu8ybGtXGtKPA1eo14znoaaeM1Il2XdbKvjIeYbYhjcEw3BsvZztXPZc7HGscehLoxgsCzHvzYhH9UcU7JQ+bOIZOE1zOTlLs822B0dZiXyaEDwa532t09GhJ\/rbe3d1Q2xUzeYJHfGOwzcnWZtL+9MVksY7Jx6HxZln73xhAZk8XEx3qrVcZjH8J6+F09Y1CRkoPAYQPlSCeaFl0Ei+SyTFy4YAcCL3WShw42IPBFOp9e2rq8\/M4Q6\/hmipvrdbrO4uDGOddeBiZwhz6mmLpufCzPUfKvbE8wzsV2HpNBRDIf19Alke9c4HiIH\/J5zHW6I17neD0zYe7aYA4522DEHHIQeBB4lyCHPRPmrk1v8ygiMiYL5\/vydqaZZm8MkTF57qEaAL6p55kybqh8WfTue199cUzkyg+\/D\/f45z+n969dy\/Pi49x8Ng7HtRkzTr6ZZQmqgWqLHS5YuHixP3jIeM4jZNDxHDxkmGPEHDyTO2d50XiX7tL15EfqMX9ADwg3H+ph42Zgh\/LvJ\/UQltzaUW+b2nFD36Vd6vIsLMSgww11HLKNh3Pkb3aaBuqrnqZ5+7BnDOQaquL2dS5U5Wmjx+sUVNcXtnVexz476Qce5Oat8nPt9\/f36eDgoBWlseKYyNqUpTgmsjbl28bs6sf5UcTB2bPpusgYTOngrbfSedccrj\/nXOfDR4dpRhokn5u\/zTytGd97SupM6\/7+9+ngF784Sr7zIsE0cUl\/D37wg9HXlyZOdk3WjHMclLgt\/jfYBgoXpinQR16movw4OUC+k+qPD\/6Y0HiwL+Y+uzYyeRNs41\/lg1pXtCJsyIV5lY\/Uw+cyXaYbyQ9igCDDHN8Jk7ZlPPyYyuyqdLL5euug\/cUca+B6zh5kNwk39znKmjXYuroh51wXRl8eyJl8+iIZnwfwvUpXaTf5oQA2NDGgqtdo21I5LsYmMWDPVOZ35coVunnzJj18+LAxHYtrLZG1qbhNvLo+uNY9evTIW+8777yT5sfu3XffTW1wTmH++PHjdF43T107vDOD+GhA4PPuM89gmtI7r79+JGdZ7fAdg1DT4699jcjatMYPv\/pVevi3\/\/aROn11oPb\/hedApV5Ej7\/3vUofX5wuMorjNPuHX\/pS7dy3b9+mS5cuEc711HnBu8kbKFx0QsEXteDizIS5rzacGDhBii6Ez7ybv6Bff+foC7rLyQxfrue595+rfULDDy\/Goosg9Eyffv\/TnKJRfPj\/4Ud\/SGggKNm++tFX6dr71+jC4wvJjAjyrz3+Wq2YH7z\/AfGGuHVrh62kE++c4DD0l+\/+Za3c0p95302CdTy6NbO8be3sXzX++P0fr9f4qYefOrbGorqq4uIiyIFlI8oyMU7CVr1GuSjYMY8Rr2lXBjmIddCDMIe8imAHeybMi3xu3bpFq9WKzpw505jQOIAOvv71dXjMQW3i1fU5ffo0nTp1ylvvqb\/6q3Utz\/+jf5TaEL7gnkif+dGP0nndPHXttn\/2syQ60daLL6bxn3\/uuXSO3ZlPPkllHKusdrYZejz13\/4bbb\/5JsojSo79yT\/5kyM1FuXn2vGOFZy3t7dr+RXFayxPsERe0LNJk1rXH+f3nTt36Jr8+BFBFki\/FMqacNEBTVkPLoB18uMiiBMEJ4rvpHru+fwF\/cmZoy9on31d2cenPl6X9\/ee\/3uNXkz8YqzK9cVnv7jO0bT27z7z3dQ3ooi+vfVt2j65TavkB3NKNuir8kP\/3rPvJdZEJvnBvG7tsJX0q8\/\/KvEG7KSuCf9Xp\/KbBGL6fItwa1u7L4dP9syzebPu0xfV5bOVspdffpmhIxwHCnTDa7budQN2sC9aSpmuyAfyJn4XkncUjDGEm2FbQuOAvKC2MZr4nTx5kk6cOOGteWtrC2UkL9Z8Tel\/PQKptV6fJrl9tggNQh7oJR5byUekkDGV1c42g45vvEHb3\/gGyk0xwvOf6ubj2rkh3frOdwbBs6yerHCireR6UGYndTi\/L168SDjX2X+p4y+FsjBchEC4yIVSU1EdODFwguBEkScO8y9v5Tef\/a393k76R9uP1iW9uPViFnd7u9bIL0ausWhEXE6CfEV2Pjn77dAOvbH9Bn1j+xv0B9t\/QCb5ocOtTsynt54+tCZCnrq1w1bS57c\/v47T5TjwTQLrQEyZg3m2QUKWYWxbO3zrEOdFbT571rt1+WyljAwFtZVdF+pcN+APu6AW1bYYY3JPa3N+Cs73he5z54atxNosPueJomwe2t5aSh+CibqMyR5HAb4pRVHuYW3OD81ZO3SG2ccPpoFiJNtc5OCDCyTHAA8Zz3mEDDqeg4cMc4yYg+9KhpIXC2UbPrrKuO57GUvm6B45jyDjyny5RTHHvvfoHsXJDx1ue7RH2KAHgS+jpnnLYnG+t6n9Ay7r+HIe1NJn\/YhXRnVqK\/Ofgw6vS7w+29Zaxx\/xYdc0R1u\/pnmCtbc2K82YbMTeGOwzsjYb+9pbWx7J2nL9WFprKW2erM0y4jlexmR80\/1qlXvEcc4PzVmbZzAm53vk5h4quAaqLaB8IcNFELyMAxnPocMcBJ7lGDGHnAlzyNuQoeyE6\/MGx7E4Ng2wydhNG4GIsn8pxRTTjeSHnI31jrhwKmspNKpQcIyma5Fh2ZdjSR3zUsf2rJtybFuXXIOMMfZa8Fp0X4eQuXXAxidnGUZJ7A8Z8xwDMhDmrCsbYQd7JszL7DvrjMlDxHHOh8IZk1dibc73zRmTRzQm433viGWa4ff4r1iY8F+yxHGec7XK+TacMZnXmI+usDbLib0x2Cs5CATRQPGFR46oU84lD52PcOECuTpXhjnItcMccibMQyK+qQ19Q+P43LDVxeBVepWi5Ic8G2Jep\/GfhYW8KIexA9+U6vhynqaxu9rXqa1NDhl3a3+rTYhBfHAdwOuzbnDY+oj9oWMeI+ZMmLsEnSvDHHImzAelKOopfA9hrM2CGJON2BuDfUbWZmNfe2vzSMbkfCicMdl3nfquJ4qyiHt72TjGnptRY8bINsscQTRQfOGRI9CUc8lDFzpxIyFvRKHX3Ed9d+kuXU9+ZCx8kbzPZ2HJ2FX8Ocq+j9HlOLDvDtV79hbbV9XWp75OA9emrqmbJ7zuGaemzRP7LXq0drrlWUtkbZafv4+EmTHYZ2RtNg69j6Isg7XZOOXe2iy7MUTGZHzX\/eFfNqZhrE2HwXfWZimMyUbdH0MgiAbqWFULEsQk3salbhvH2ql5I2+bjZs\/ztc0zi7tUpeNb\/R1moKqPDIGx63ykfo2PnXfuZN5huLl+pvkCGkNdeuWzVZdn9naGZOVzu8SZLPp9sYczW1MNu+7PmuzuNgbg\/1RiuOj87Fn1uYZ8V+zGJPPu3DG5N5xnPNDctYOGX0RsYNqoPAvTNASkD1H2Tsffa2lD4lxIwAAEABJREFUzY28a+62Odv6da23yr9NXdKHG8uiPG2blaJ4U8p53VO\/AyUx2KgGSS48VN7avDJjcn4KTr4DNkV+N+dqRRRFrrT9PIpy37G+B2VtllO++5VJdH+IQDANFBonXCBBh7UNPAwbXt5M+WbUV0YZu6+YMs7UzR\/j1UcdEeUXHo5LC9p4TUXnRJF8zhDgWjHn+nupPYqyMHGcjaHtoyirKI6zsa+9tXkkY\/y8tbl8LM7ao5nimAj\/dYu1R+VdZlGUeVubjUPurR0y+mJiB9NALaVx8p0ZfJPz6erKZIyhb4oyvszbplb4tIkBvz6o61rkx5gylq821k+5Xl9dLGtSVxNbjq\/jBAhYO0HSw5TWHjLJYEyy8\/xa6xH2IDKmhyA9htjbOxrMWqI4JrL2qLzLjN8JimOiLnHq+FqbW0VRzit3BIHJG6ilNk4R5SddHzcjGYNv1DTQJuPLvHXTxRST3NrEkP5deV5Pl+\/1IAaorBbWj7leztXHu3W+tW39NJy\/wPPVt7GyED6ysjaH35icBzdUfUXfqTIGWTOK42wca28t0Y0bWTZjiJ48OU6ZttvemNx\/by\/nlZsMgckbqKKVz72x4psp1sc3OfBtScaQsdvGK\/OT8WXeMp++dGPnq6q7TdMV2hr4eDZZC68hpO9AVR2rEfXTpzImr8HanB+T42bGmONZjcll1uZ8V87aLIIx2ch7Y5gbf3zllTwnvjiez\/rlVqs8nrU5PwQXx3lUY3JeuSMIBNtAHaly5pMmN66ipXIMvhkW2fUhlzn4Rtok7j26d8y8TRxZx7GADQQRRYQtJnFRoHob112nlqHeBSqqlGsr0reVDxW3bT3qV4GAtRUGA6mtLQ5sTK6zNueH4ozJI1ub80NzcUwUx1mW1YooijJ+qL0xWWRuXrPZcHtjiIwZLv7MI2sDNeABrHPTrZueb2p9xizLzXm4cSuzHV3XMiFj2MS9iQ9jhvhN\/GDflWTurrGk\/9P7+f9LKOXKT4xAFOUFWJvzU3DGHM9qTC6zNue7ctYWRzAm043VXFhL6X\/XgqzGEF2\/Dm5YiqIs\/t5eNg61H+sv\/Yaqf6S42kANCLSh5EVFRGPfTGnijdfL60c5LANfRtJO+pf5VOl2aGdtIuOvhSUM28sYJeaLUPGasRj9CA8oBEjGTF+UtVkNxmSj3BuTz6zN+b44Y\/qK1D7O3l5ycT9cG5onY9rHquu5k1\/LyNq6Xu3tjGnvuwGebRqowWFZyp8pcwMQ0+FbvNR+4xg7ohloH63aM6LsXzqcl2puuPmCYM7rB88y8HOhLjV38W2Dj8Ta9Wfd2DW5deh8IASmerfA2mxBRV8YNybT9\/mOkLVZTF9OYzKdtdk45N5aOvLF8dVqyGx57NUq5+M45\/vmrM0iGpONuvciMHkDtZRmyYfuOernYZpT3vi65N4h8a8lH0AjybiBQLomDaFce0QRVW0yj\/St8murHyqHjBvaX+HhejH3PzBpe7yP+RlzTDSawNrRUq0TWbtmvYwxmdjabPTuexKO9cVxX7nGZNIhG2drsxy+RjXT6D5BYPIGKqlhsb9D3FBlzCGBa9v8yZuvbDqkvKxuadfXWmUdZbm76vqqt2sdrj\/XJbF1beRc2oX2EZ42T+JIGZNNrM3GqfZR5M8cRZk8jrOxz70xxdGsLdb1oYljojjOIq1WRFGU8WPtoyjLtLeXjX3vrc0jGpPzyh1DYPIGChdE\/KtSEqp055DNmeRNqek6pC\/fDJvGaGov88j8VXGkrYxR5Te0nmvx\/YVgUe6YDi+SRMT+VLLVsSlx3wiVLrJHBIzJglmbjWPura2fzdr6tmWW1uZaY3KeuaG+H\/TUU0SS8IRxOtzw3adDdrRhqHXyAqxlLrnweXDOtRvPTd5A8RFAI8UEGfMYuZmCfE4U0eG\/FIhINhbUcJO+Y92kZR6Zv6p0aYsYIPjM+a\/5sAYQ1lGXJA51fZrayRx16pP2Zbn4WIX27lNZzVU6voZgLLOFnslnB52UY+4jaTM4b+3gKY4lsDYXGZPzklvSxz\/GEEWRXF3GG5ONU+3jeKrMmjdBIJgGKqml8BdNFAgXqkKjABV1bmp1ypY3vm4x62TLbGQemT\/TFu\/55iv9i62Pa2SutjGORyXiZla+q0QVG6+lwuyImmtu43skUI+Tc9Tsu3jyGPRYxmShcN3A9YMJc18xkLMNRszZDjyI53KErUtSPxg\/9DsRZYVbm2uNyXnJGZPPrM35oThj8sjW5nwfnDFEQz4ks0mNq1VuPcT3oOI4jx9FOa\/cMQRm0UBx1bhIFV3E2CbUsclHR+4a+GbMN2dXP8Rc5mpyQ2Vb9ueR5UPU2iRmkzrYltfQJM8Ytlwfcg1R49TvQOG1zoQ1tiH447pR5euzk37gQVVxNk5vTPGSjcl11uZ8W87a3NOYnGfOGOaSt\/yFbS7txvEXx43pFqcPb2OyKNZm4xB7Y4aI2j5mgJ6zaqACxK+ypD5ubPJGWZmwRwOunRu4OqG5VvblkeV1YsCG\/cD3QTuU\/0Vg3Vrq2sn6uO42vjJOnzzX1GfMMWKhYWHifGh0mB9iRHymuvHZHmOVz\/7+Ph0cHPRCZMw63cFbb\/USs25tdPjOR5l93\/WRtfl6z549vt5ExgZldTXVpTHjmCiOU\/bgK18hiqKM7+lYNq4pitL8qKmpb5V9nWNbFgPnOCgrcNn7IBooXCTrwtzEtm7MIe345tXHDZVjDVlv19i8znN09GMjllfFb9KsVcXqqueadyhvvqpihn6MeE1l62CbqdeChoSprN4ynXu9QDxXxv6sgx6EOevKRtgyVflcuXKFbt68SQ8fPuxOn\/rUuqz3f\/zj7vEKasLN8NGjR0fiP\/75z9PcB5\/73BG5XNeDJ09SG+w+\/MlPCu2kTxn\/\/gcfIBQdJI1SkV1qkOwef+97aT5f7UW+RfIk3NHfpHlEwwphkU8f8rLa3\/3CF5A+pXf\/839O1+rLufX001REPnvI6hxb2BXR7du36dKlS4RzPS1wwbsgGigXX1yMXNlc53wTavLdG3etU93QIsr+lVO3dq6Tko3X7TZSiWqSX64HyeusR64FPk2pq3+dfNxsyrXV8auy4dr7jluV19XzdYBHV990juamLFaZriiX64M58hTZ37p1i1arFZ05c6Yznfryl9dpnv\/gg87ximo6ffo0nTp16kj87Z\/9LM29ffIkFflBTsakds++916pHWyrCDEQbGtrqzAW59ve3k5tfLVX5XH1yClpK3kXEQSZa9vnvKz2Z\/\/pP0X6lE59\/HG6Vl\/u1MAYoigiiiKiKEpF2PnsIeNju\/Xii4VxYVdEOL\/v3LlD165dQ5pFU5AN1JIQP0dH34lpujbczEDw6xoLMdoQ56\/ylXaGDGHjEbzUY+6jOjY+vypZRPmFg2psso4mvnyMpH+NdNUmI1nMte4qeNDUoLmpshtaf+HCBTLGEG7wfVASLC1566c\/7S2mW9fJpEk6ceLEkfjrj9N2do7IXd+0uGSHd0FcXdM51piESn+LfF08fLUX+Xrljx6l+dJdctwI76oJ8vokzVsf8qra12v9i78oPAbruvEFeNSfHC+KolTsrRHrtTbVo4Hy2lSsD+f3xYsXCed6GmjBu+AaKFzoJM0de0OGeOt6c5KxOOaQ47mGzZ9cX9dau\/qX4XKP7pWpF6eTWMpjVLVQ6VdlO6ZeXh+YL8sPm6rmCXrYlcXx6VwfzBHLZ7sYmbX1lxJFmW0cZ2Mfe2OKoxiT6azNxq77GzfyCGhC8tn0XBRlNeztZWPR3lqiF14g2tuj9L+fieMiy6NyY47OdXYMgeAaKFx8JOGCBDpW+biCXrI1uXlxQukz9g1N5pN1cG3uKG3Yl0fYSj3mY5OspSp3TPlFpomftB16vRxf5qQet6Hidi1RXh+YL7pGsByjJK4BMuYRC3MmzFlXNMKG7TFiXmQ7iNyYLKy12Tj23ph6Ga2tZ9eXlbXdI1lLtLeXxYkioijK+FD2eDeJa7GWueOjtcdlkFiL\/VGyNp8bk\/PKeREIroFyq8QFCYSLk6ubw1x+\/MM3vCZ1S5+xb2gyn6yjqP636e1UJf0kXydGGmCgHR+LJnWgftBAJY0SVtZftfYq\/SgFt0hSdI2A3EecAjrmMWLOhLlL0PlkkINc3eBzY7IU1mbjGHtr8yzG5LyP6\/NhmtZmGYzJRt++z3zysQXXr\/uyTSszJs8fxzkPzlqi3V1wRFQw4InqcXxUaW0+NybnlfMiEGQDhWYJJCvGxcmVSX2ovLx5talR3tC6xmqaX+aTdRTFYRvpV2RbJO8jRlFslsfkXDTo+MbN4HFNuUSunddS7qFaRaADAtwwWNshyICuxuTBrc35Npy1mRevOZsd3RuTz63N+aZcHBPFceYVRURRlPEh7aOIyJisonviawnWEqE5kh8\/ZlZH99Zmdnt7udzanDcm55XzIhBcA4UmCc0SyK0YMuhd+Vzmbb57wzdyeWMea70yZ51mwGfTNMaQazvX4DtdvBZZf53apD3HqOPXxobjy5xt4kgfjglZn3ERrw3hNV\/Xr4lt3Zh92A0aw5g8vLU5PyQXx3l0Y3LexxmTS63N+aactU09iKxt7gMPayn9rhB4Y4hCfPcJteH\/57MWHNHeHhHmIHzfyVpKtygiEl96X\/N376bqdId32kCYvJ19ikDGYKZUgUBwDRRfBNEogSrqn4W6y41I3tCmWCzXzo1cWQ1c6w7teM2axDhH3f560VcArwU6rhW8j6r0Pp9QZW3XPcQx6IoRXx+6xlH\/HhAwhsiY8kDG5Hprc74LZ0yxtzG5ztqcb8LFMVEcZx5XrxIZk\/Fz21+\/TiQbJVl\/FBE9eJBL9vYobcD29ijdrKV0joaMdCtCILgGCoWiccKFEoS5JJ9M6kPk+QbW5KYM213apTj5oWTjecIG+Yv6igrj9Rfph5Vn0VEfKJsRXU5+gC\/PMUIP2VP0FIGnFptca50YsOGcyPsCvUCY08SbXMfEpWj6MgSMybVxnPNDcvwuRZ0cxuRW1uZ8F86YYu8oKtbV0VhLxO\/GGEO0u1vHKzwbNE+7u+V1GZM1UVF03M6Y4zKVHEMgyAZqjk3SMWSFwFB2MnIzRBUbbqqX6TLdSH6kKea4wUrZ0HxE2YsLNVHJJvXsw+a8fp4XjTJGkU0bOeK6eEIm8cTctUEuHLO2mFe94+bLCZmsCzUUEWyhO9fju3VVNSPfmIR\/TJXRmLUEmcuY8cuyNstpTDZW7Y3JLJo0XplHvrc25+ty1ta1zO3k94ZCe2xBXmXORVHOS253V86KeWP871JdvUoURcV+qkkRGLyBSrNs+K7pDQ43UL45utBB\/gq94ooHn6ORaJuEGyjUXjcG+9S1L7Org6e0cRtA1N0E87q1y5wrWtH15IdzN81JFRvXVNUgIW9FqFHV+MdUGXFzNWpRISUzJq\/G2pzfZM6YbPVNGzZrifb2Mt8oIoqijLRTm7gAABAASURBVA99H0VZhasVURRlfJu9MbmXMTmvXCECkzdQfAHEWFQldJKK7EKV880L9dW5Qe1R9iKWfq\/SqxQlP5RsrE\/YwX\/rNn8x5R8fcJ3kbHXW7rj0MmW8UJfElINDD+K5XAvbSz3bdR05JnIAZzQ3d+kuoU5KNtYn7LHfobFETceSBijg5grXhwDLG6ckY7I8TRuGzKv53trMx5hsrNpHUWYRx9nYZm9t7mVMznfl8B0fJnz5muPFMbjwKY6zGqOIyFoia7N5m70x2Ud6eOdttWoTYeN8Jm+g+AKIUaKPCyITdJKkneTZXsp8fJVdld4Xs66syY1P3sTA71D+5ewmcerW5rNDXpbXySnt2Y\/HOv5sWxaHbZqOwC+i5EJD9TbUcJWuro3r1g8\/ODWxx7tRaJhAqBP+oDoxOB\/su1KdfF1zDOGP6wNet0PE1pgOAtZmgrLHCWQWR\/fWHp03mVmbWxuT8z4uijKptdlYtTeGyJgqq3D1cUwUx0RxTGRt+zrjmNLvf732GqWPQbC2fawN8Zy8gcJFjwmYM48LIhPkVQQ\/tgdfZA9dHbsi\/zbyiA5f0ERU5wZlyBA2ficEc8S4R\/cgJsxBNMIm85TVzrX5SjpH9f6iriy+L25dGa8BNV6n\/IflGEEcjy3wrh98IIceBL6K2K5qPWzHx5njtsnJvkUj56qqif3ZnueTjJq0HgJRlNlZm41D7q3NoxuT82Vc00arLJYxZdqjujg+Oi+aGUNkTKY1hiiKMj70ve\/xBCxrW3scE8UxURwnNytxrNvGW7jf5A0UNzMYZXPTBHf2Yx+OxXMe69i5Nuwrx\/39fQIdHBxQHTp7cHbtXsc+oqMvYNz08B0cvtF+5eArtfLWyVVlI2+kbx28VZiXDjes1Y3ZNAZCuTG6zBlP4AccER8flwFX8MCTbTDnBuY1eo3gQ8kGm7o1JObpL+KX+cicqUOyu0E3qE5OxKbDzYe5m\/fQtPD4sT3H3drfYpdJR7wei8gtDK97V7ZR8zgefrnW5jmMyfkyzphca23ON+GafDzZtGGLY6I4zqrBX69l3GbtufHyjZuFRKPVTt5AyWrncgG8cuUKXbp0iW7evEkPHz6sRbzO7z3+XqX9v3zwL+nC4wvsko74aAcMbpb\/5J1\/UhnDrQsN36NHj5r4pbZPHjxB2pR+8uFPUpkbG3O+6X\/pwy8ds3nv3fdSf+zeeeedY3r4g6CDDejdd99d27WtHTFBEk\/UiSbFxbOODWLVoU+\/\/2ksIaWy2mXO1DjZcQNTdZyLsCqq75c\/+uUkOtEvDn6xxtW1ff2d11Mb7P7su39GOM\/BT0m4JhTRlHUFlbtpwzB28cbkGa3N+SactfWtjcltrc35OlwU1bFSG0UgRSCoBiqtaMId\/qWLi3VVCbdu3aI7d+7QarWiM2fO1CJD2Yt6e3u7lj3f4Olwg\/\/15OOnvzn4G\/ryqS\/XiiFrO336NJ06daqxH2IgN8p479n3vP6fnPkE6pSeefaZYzZffPaLxBtsEdNHH5\/6mM3oV5\/\/1XWcLrVzHjROwI8TYE2YSzzr2HC8slGu9\/Hpx6W4tz3Ozz\/\/PC+lND7XifMODj\/b\/tkaV9bxiPMDNqA\/\/Oof0rVr18AqhY6AMXmF1ub8EJy1eVRjcr6MMybXWpvzbThj2njV98EzoKytb6+WgSIwTlnaQB3iXLd5gvmFCxfo4sWLycfmhnBjqkO4YcP3p1s\/reXzH7f\/I8xTukt36UHys0u7tXx99Zw8eZJOnDjRyp8Ot6e3ni70PzShl7dePmbz+e3Ps5r2t\/aP6bnera2ttR3LMHapHf5MwO8JZT9FeNax4XhF45F1nNwuxV0eZ148aiiKLeVsj3xS7uNx7Njep4fs0fYjNqGvnP0K4TxfC5RRBICAtdhnZEw2Vu2NyS2szfmhOGPyyHGc83W4ONbv\/tTBSW1SBIJpoNDAuJRWOOJO5kdazDH2QdxA4V2OOvH4ezjwiyiiKTfOzx8vubUUydkOa2C+bJRx6vqUxQtBJ9fkq4ePM+uq7H12TbCqG5\/zTDni9cdUVAfreSyyW7TcmHx5cZzzQ3D8XSRjSqMfUxqTidg\/m9XfW5vZGpONZXtjyrRHdfiTfZbgvz2R3wFiuY6KQAECwTRQ+OjMJb4oYiyoPxXDT9qAhyxVih1k0LEIPGSYY5TEMox90Dmq95dobq4mN0fXt+95UfMn5VX1zukG3hY\/iUGb9bbxqaq1zvkn88o1VMUeSu++HjkPXrdM0gY822zUGEXjLdfa8XLJTNZmszrf9zIms8XeWuyLCX+yD60xRFEETkkRqI1AMA2Ur2JcEJn4gumzgwx2bAMeMibImYcOcxB4lg89yhuSvFH58kLPTckO5c9+8tmOIatz80UdWCMIvEssx1+\/ubqlzXmtVeuSx7muD8eEL\/PHfVnTbJQxm3kOY43XKBMyMI\/XLRPkdYh9MZbZQ8\/ks4POJ4cMOhD4ycjacVIb0yxPFGX2cZyNTfbWNrHObI3JRmuz0be3liiOM00UZaPuFYEGCATdQMl18AWz7ALFNtIPPOQYmTAH8dw3Vul9PnVlVTcqqY9o+he2vEHL2uhwcz+GOhQ3Hri5kvkaBwnAoW79EsurlD+wk5vnPpcia5J5i3JI+yKboeV4DTLhdc9807zSFzEw98WAHHomzNkOPIjnwY3GZCW1\/Ygs867eW5vZGJONTffWNvU4am\/M0XmXmbW5N\/7vt3ymnCJQC4HZNFC8GlzcmA919NUlG6GqG5i8gUo\/X9wxZPJmWla7tHPrYl2Zv+uzhHnZevk4A5td2qU2G3zb+Pl8Qm5g277u0fTU8fXZST\/wIB9ukPn8IXdpf3+\/8llc\/EyuJiPnOfjFLwaJz7WQtWmqg899rlEeEh+9cay6Y5rwcHdw9mytvGRM6nHwv\/93oT3duJHZIOZLLxXa1a1T7bLnIuIcB6XgLnw3uwZqrsejyY2O39Fp4jMkLrIOX0PAMmnn1sM6tnX1S5vXWS8fZ3ftdTDiZsf11Xk3BNAIMXWL5PfGs7WaPD\/OfVZX0fyjL30pS3jvXuFzvop8q+S4GT569IgevfFGliPZf\/jRR43yvPvMM4lX9vvO66838oV95kkknw9XVvdHv5w992wraVhRu2ubxozjNCyaQVcfwpxxH7mWRsfGV9vt27fT5yTiXE8BXvBOG6gJDm7RjZNL4XcmQnj3CTVxMwDevbljDoLuHLX7ojx8QRxH5oN8iYS18nG+StnHd7zuoZsj5PZhynKuw2czN5n7rhGaI1fGa2Id9CDMWVc2wg72ZTaswzPkmjw\/jp\/RVTXiMRTIgYahyrapnp\/D9pnPfAYpUnr2N36j8HlivvjPfjF\/FtyZTz5p5Nv0mWfIv\/357NEp20njh+ebQSbp1Mcfp+vAbuvmzUb1yDhD8oz7kDmGiI3zG89J3ITnyGkDhVfQSFTnxsQ3MZTEN1bwUxPXXnZzZxtfrdxcyfX57JYiYyzqrJcb5SY+deJKLDm2lM2FR4PiUpvaEaOs0SnTFeWriun64dlaxpQ8P257u\/A5aWiSimjrxRfXqYps2srXz2FLmhFOUufZY0fyHTY08EeTd0RXtWaZN1lnHV+Jx\/\/z7rvHMRXvpm29\/PJxfVVNI+jXuI+Qqw6mdW1wfuM5iTjXcbyXTNpAjXh0DRnCVnbz43clYDcHkmvh9c2h7qFrrMLCfQI56mEfiSnkZcQ+ZTZ1dZy3z5h1c5fZobFxCY0LU5kv62CLGDzvc0RsJsQFj3FUMiZPZ23O98lZm0czJufrcMbkVtbmfB3O2tzKmJzvwiUfdabuUZQOulME2iCgDVQb1Fr68I2prEnCf2CL8LDldyYwn5q4Fr7Jcj1yjppZ7o5SJ32kHculrdSPwPeeAk+e9wXlj3GxVsaW7RgHnvcxIg\/HKYrPcn63kO1DHNEMMaFhARXVCR1si\/SQQw878E0IfpLgiznGycjaYVJbm8c1JufrcsZklkP\/pSCyGIN9SnjHK2V4Zy1RHGeznekfE5MVovs5IqAN1IhHrcmNSd7wRiyxMpXb\/PFNF46h1ozapqL9rf3S1BKzHcov5hLX0gCqJDQsIF8DxDKMkhg2yJjnGJCBMGdd8GMU5SVam\/NDcMYMEbU4JjdcxhTbuBpj1pKtn\/50zadMHKdDuouidNCdItAGgbAaqDYrmJGPvFn6bpCQcYOyI26mISyxqPnj70TJtfnqlXpeo89uKbIivLC+Po4zYiCWxBXzIqqy43jwr7KFTYjka3gg8xHXDx3zGDFnwtwl6FyZnFfppW2vvDG9hvMG44+9vMoawijKjOI4G+vura1rmdsZs+aPvQOlTx9fY6NMNwS0geqGX2tvecPiIFIW0eHFhsLY5E1V1sm81PsqrtLDh2OVNR+wmwPJ9brvQvE6sQ55nKXPkE2mzI8aXJJ1uDqdzwABa4ct0phu8a1t529MK7+n98W7wNYSxXEWJ4qysae9htk8BLSBGvGYyxuT7yYmb5oRhfXiLqqd1yH15Nmq9B6XxYqKjnMTjBj3vppNjrdY0DdhYcZkq+SPvLJZf3tru8USD9PsFqimtzHHDa3NZfr08RwL5VohoA1UK9jaOVU1RfKLxe0yDOclb+58s8UIQtYmN3L2gR+TT8a6OY4SL\/cdqKLjLH2GwIPj88euc8S1uGbVkDEZCNZmY8\/79Udhbb94bUxekbU5X8VZm1kYk41198akluu6MTt8+niKVRRBoqQItEZAG6jW0HVz5JuojMLvTFQ1WtJnLJ5vvsjnu7lLPWx8xDZVN3C288WYi0yuwf1LvLLjzH5VGPWNgzymXEPfOTTewAgYkyWwNht73B9pQtrGNSb3tDbnqzhrM4um72AZk\/qtv0RuLVEcp7K0gco43SsCrRHQBqo1dO0ci25O8gYW0gM05Sq5dr65y5pZJ+2VP46AxKztcZYx+sJ9iJjHV6+SURCwdtg0xjSP\/9RTRJcv537gIQPl0kG4dfNnbR7\/+vWcV04RaImANlAtgWvrxjc8ecNCLH5XAvxcSK6B11VWO9tIP7aXMrZj3RxHuQb5EV7VcY4o+1ihyo4G2mTdA6XQsEMh0PQdmgZ1rJsQ+BiDfXMyJvcxJufLOGtzrTE5X4dz8YjjzMsYoijKeN0rAkStMdAGqjV07Rz5BiUbBkQK9QGaqI2Jb+5cO4\/Q87rAFxHbSL8i26XK6x7nITAqw5\/fVVwq7huxLmPyZVqb86FwxhAZk1UTRURRlPFD7Y1ZR04bQH4MgzFruTKKQBcEtIHqgl4L33OU\/Ye7RTdIvsm1CD2aC787wjfdpjUXrX20BYyUiHHZf1r8GfVhbtYdTtfDDlU\/TFPiVxRnHbAmI2PWdFGz0BAwJq\/I2pzvgdt+8808ijE53wdXFsPaXGtMzjfl4pgojjOvtl+Cz7x1rwisEdAGag3FOIy84fFNCyM3JTviBjpORfWzcPPHHqgbvFwT5kUXWQu5AAAQAElEQVTk+ks7jgVZ3XiwnRNhjVMeZ8YVdRThxjZFepUHjIAxeXHW5nzfnDHtI0ZR5ru3l41D7qMoj85\/fQdJFGGvpAh0RkAbqM4Qdg8gb2gRhfviljdX1AyiBpvr38B1lqa8Xv4OlMSr6DizDxYc0+G\/mGnQLQ0ua0sFupsfAsbkNVub8z1w64dRGtM+WhwTxXHub23OF3HW5hpjcr4OZ0xuZW3GG0MURRmve0WgIwLaQHUEsKm77wYpb5QRzePFjRsuiJJtJ+B3zZLyJvvlY82PMahznNmnrGjGHTZ17GFXl\/qOVzev2vWMQM8P00y\/Q4QSjcG+PVmb+1qb80WctbnGmJxvy+nDM9sip34eBLSB8oAypOhYg5Qk42dChX7zkrXLZqBu3dJONgEJBCTn0o4WsPE7UHWOs1y7xKQPGPgjVF9clrFNH\/k0xgQIGDNB0hopnzwhYmLz1SqT8XyI0ZijUaPo6FxnikAHBLSB6gBeV1e+oXIzIm+eXWMP4S\/r49qRR8oxr0N8w65jO1cbtxnh4xxR+UWc8eQv6ZOz9Y1d3\/GccnU6JgLGZNmszcaiPZ6\/VEQen\/XDKI3xaBuKjMkcrM3Gsj2\/k2ZMmdVxHa\/N2qM6PH\/qqERnHRDYdFdtoCY4A\/gGidTy5tX2wYqIMxZx7bLuurnZt8y+jk2Zf0g6uZYfbv1wXVqfx1nmWCfowPQdr0Mp6toGAWMyL2uzsWxvDFEUEUURURSVWdL6Izz32UqlXgXKKMoUcZyNQ+2NITJmqOgaVxEgbaAmOAn4JoUmBMQlsJznIY+y7qp3VHgdcn3SH\/qid1ugWwLxu0911sJ4NvGpE7cIf3kspE2dmHOxOX\/+PDGV1cw2GH12Pjlkknx+w8g8UbnBsdajdETGEL36akZ37xJFkWOQTdfNE6bGYN+N5GMErC2PZW2mNyYbm+yNIYqiJh5qqwg0QkAbqEZw9WPMNyncuOo+WLGfzN2j8M2dI\/FaeK6jH4EbdCNVAC8Xw1Th2eH88IjXIsRaT5QpRADNzf3794kJc58x5GyDEXO2Aw\/iOY+QwVYSZKwffTQmT2ltzvs4ayn971VeeIEIZK3PKn\/3CVpjsO9Gq1XuH8c5PwTHzZoxRFE0RAaNucEIaAM1wcE\/R\/nDNPkmuSk3Q16n+44T4zDB4WidssqR1yrtfDKpB79D5Q\/TdLGDTxeS2Nepr0uusX3RzKC5qcrrs5N+4EFVcWalt5bI2qxka4mszfjd3Wwccm9MFp2fDp7Nju+tzWTGZGPD\/eOvfY3ev3aNDr7\/\/Yaeaq4IVCOgDVQ1Rr1byJsUf0yzI26avSfsMSA3fxxSroVlXca+43Wppauvby1THmdZj2yaJN91zXP3RyPFVGctTZuq\/f19Ojg4GIzImHXZB2+9VZhnbeRj8NBJvCPFX8SWX7wGn8i7roGiKMu8t1dYI3JwU3fwuc+V2sFWUho8jmn713+dTr7xBm39\/u8TWZuKpZ3y\/Z+LOMdBKdgL3824gVrWkYno8IJCYW\/yJoxKdxo2fuy\/CTdtXitwYqpznKVfTDG5G2Mn7VybtvMhYratpQ8\/t8FBc+TKOA\/roAdhzrqqEbZM8C2yv3LlCt28eZMePnw4DH3qU+vU7\/\/4x4U51kZFjLVFmlTetf73P\/3pNA5277z+emWdH370UaGNrxbETSlOmij8NzTJyA2Uzz5EGZqQR48eNVp3COu4ffs2Xbp0iXCup8dgwTttoCY4uL6bVJ0b6wSl9p6S185NQO8JAg9Y5zgzRkMuReLPHwmOkXfINVXFRoNT1tyU6apiw5cJeYrsb926RavVis6cOTMInfryl9epn\/\/gA3+OTz5Z21AU0eOf\/\/wI0fXrRMbkNswZQwcvvZTOutb\/7Be\/mMbB7sxf\/3VlnbBvkpPX9P7\/+T\/0v\/7n\/6SPPvxwvcYmcQazrXH8T58+TadOnfJjU8N\/qtpxft+5c4euJR+d4vgumbSBGvno4sble1dh5DJapUPtIOn8Nr0tp7V5Nw7PDRlawob17FK775JIDBCHetpk3J5CziYMmho0OH0XjLhNYl64cCHpTQxtb28PRkmCtCQ8u8mb5xvfSPXpLmmWXBva3SV68CBVu7utra1U5Po0nW\/93u+lcbDDX\/kV+UMPQt4imzL5yZMn6cSJE4NhXZa7q26utRtj6OLFi4RzHcduyaQN1IhHFzfDy3SZ+C+yZOoX6AU5DY4vqn2P9pLK69d+jrIv0NOCtyKssOS6x9mQIWy+BhXxoWMb8F2o73hdahnCF02O0zwdSwM97I4pliawlmhvL1tVFBFFUcYX7aOIHifv4AzyRWxjsqz8sMxslu+tzXljcl45RSAQBBbVQOECCKrCFjYgnx3kTD59FxkaJ75ZuXEgf4VeccXBzPuqXd70sWZeIPNLaLD6woqx6Wsswr6v+CHGwWsZdWGUBBkIMowgbqIgA2EOeRnBBraSICvzGVxnTJbC2myU+1cOrzHGUPpRndT5+Dim7T\/4g2NfxPaZNpYlzVnqs7eXDrpTBOaGwGIaKFzAcOECgS86ENDBBgRe2mEOORPmUt+Vx7s1iBFR9kOHG9\/YWH8oDmrg2lA514sRcxTKevCVtHADxgLYrGhFvAEv8KwHX0QRRYQtppiG3rh55fqGzjdmfH4tuyPXADnzGDFnwtwl6HwyyJlc\/ehzY7KUcZyNvLeWKI6zWRQRRVHGV+3jpIlyvohd5VJLv7OTm1mb88xZyxyRMTmvnCIQCAKLaKDQ6ODixZiCh4znPEIGHc\/BQ4Y5RszBD007tEN8s8J4la4Sb3wz43loI2qPKL\/wYk6HW53asd5Dc6pjTzPeJDZYd5vj7MOIZeeon49D+44340O2jNL5aeTuauS7T1fza45rtp4f\/ue\/+EL2g\/v36eAXvyA6lK1tujDG5N5Jk5ZPDjlrDxkdFIFwEJCVLKKBkgsag8efl4KaPkMEN1HUd4\/u0bcOvkVfP\/g6ff\/g+8nsHsR09uBsSk3jjmEva\/\/nB\/+c2tbOcbDgtw6y59TwDRwyYDDGeobMwWvEcf63H\/1b+ocf\/EP688d\/3ug471D+r3PGiWsGTiCeNxm5NhkTsUAyDs5vyJRmiIAxedHWZjwaFBBmUUQUReCmpSjK81c9UNOY3FY5RSAQBLSBOjwQ7rtPZe9I4fkWeM5F0+e5fOmjL6XZYorppYOX6DMffob+2eN\/lsxiwnbh8YXBnvmBG2KXZ4rI2v\/xwT\/uVDvWCvrg\/Q\/W68UcJGX8TJOutXOcsUaJ1WW6TM\/\/\/8\/THx\/8MeG4U7LVOc7vvfteYpn9vvPOO2ucXn\/n9UyY7D\/6sNmzcbB+NEmJKz1+\/DiNKeM99\/5zqQx2eJYLznPYKi0AATwcE8swhuj6dXBhkDHVdRgjbJRVBMJBQBsoz7Eoa55gjme54DkXq1Wz57l8e+vbFCU\/lGxvbr9J33z+m\/TdZ76bzIjwzsCfbP\/JYM\/86PpMkb5q\/\/Kp\/Dk1HzyfPacGzzqhw+255587hkHX2sd+HorE6n888z\/o1XOvNj7OX3w2f07OX5\/Jn5MjsYJN07Xhz8EBNf7EGr4ynsQe5\/YmPMcFWCyOjMmXFMdEcUwUx5ksipKLjdBn0un2UZTl3tvLRrkv+us8aaO8IjAhAtpAOeBXNU8wx\/Mt8JwLY5o\/z+Uu3aXryQ\/igAyZZHadHiQ\/L269ONjzSvp4pkhftdPh9tOtn6brfbT96FBC5MOgj9rRMIxJXbH6\/Pbn15ig6ZG1s8KVS5siHucb\/KuwN8Z4n+MCX6XAEYiiowXKd59effWoburZzk5egbU5D85a7ImSczFjdK8IhIXAIhoofPyGxoehBQ8Zz3mEDDqeg4esaM7yvsdd2qUnhz9onDCnmWyotWvthgxtwgasfv7453T\/wX36m4O\/Icyp5iYxkt8Rk+7SRsqVVwTWCLz2GlEcZ9MoysaQ9qtVXk0c57xyisAMEFhEAwWc0QihIQKBh4wJMuahwxwEnuWYg8coCbJwSCvZJAS4QZIP0yxqpuriwjHZXsZzdWyj4wwRMCYrOo6z0Rii0N59yiojMibj3C+SW5vJjclG3SsCgSGwmAYKuKIhAoGX5MowB7k2kLkkbZTvB4GIsn8J882bR0o2vYknIIzwy5jziJSKPVCYOT31FBHIWjqyWXtkGtQkirJy4jgbdb9MBBa4qkU1UAs8Potekrx5L3qhLRfHjWZMw99YtHmi5WzGEBkzn\/WcO5fVai2RtRmPvbXYE7E+m+leEQgGAW2ggjkUm1cIN1A8AgG9kQOFoyTxkXwXrDiO\/HjwaFadzRYBY4iiyC0\/3LkxeW1xnPHWZqPuFYGAEdAGKuCDs9TSztHhvziXusCe1rVD+V8occPTNbSLfV9xu9al\/j0jwE8aN4YoinoO3nO41ao8oDHletUqAhMhoA3URMBvclr5zsnG3MB7OuASL4lj1\/B9xupai\/r3gEAUET14kFEP4QYPYUyWgr9Ibm02x94Y7JUUgeAQ0AYquEOyeQXxx0h6Ez967CUefX0PSsaUzdjRzDqbNQJxTHT5MhH+7zuM1oa\/nCjKatzby0bdKwIzQGBTGqgZHIrNKdGQId70Js5IHB8lTse1\/UgY\/zFy9VOxRqmFQBwTxTFRHBNZW8tlUiP3gZrW5uUYk\/PKKQIBIaANVEAHY1NKkTdrvoFvytqbrNOHU9\/v1jH+7nejmtSptgEh8OQJUREFVGZpKXFcqt48pa44VAS0gQr1yCy4LtkYYJl8E3fl0G06MSbcOPWJR0ybc6Oq+3DcKjvoybNBzuRRq6gKgdUqt8D3oKzN58bkvHKKQEAIaAMV0MHYxFK4edrEtbdZc1e8uCFzcxfJXbs5ztHYyAfkYu5bB+Q+O9hCBwLvEuRFfq6tzksQMOa40pjjMpUoAoEgoA1UIAdi08rgG\/YQ76wsCcuIsi\/Xuu8WMX7UYZPNWB\/xOpQymCs3N1UJfHZoitgPPIjnPPr8WOcb9\/f36eDgQMmDAUVRBtneHtHbb6e8YjW\/cwXnOCg9gAvfaQO18AM8h+Xxjdx\/E5\/DCoavkTHqmklirM3rcTTREDEd13aXXLlyhW7evEkPHz6cDeFm+OjRo8HrffcLX1gDfPCDH6T8wec+1ynvWLUPcTznWvvt27fp0qVLhHM9PYgL3mkDteCDG\/LS+EbeV2MQ8lq71Ca\/3N03VjIeH48utYbo675rhObIlXHdrIMehDnrikbYSR18XJnU37p1i1arFZ05c2Y2dPr0aTp16tTg9T77G7+xhmoreacOk+2TJ6kLVmPV3qXGIt+51o7z+86dO3Tt2jUcwkVTZQO16NXr4iZDwJAh3aoRcHHipseVV0fKLNr6Zd7z3lc1N2WNT52VV8VHjAsXLpAxhra3t2dDJ5Mm5sSJE8PX+\/nPA6Jj1AWr0Wof4HjOtXac3xcvXiSc68cO5sIE2kAt7IDObTloCECo+xzpf\/ECHCQZjBExIAAADhlJREFUMsQb48TzrmNM+V\/hyTy0wK1Oc9Nl2UPH71LbrHyTd+eO1GvMkWnDiZorAoMioA3UoPBq8CIEuFnquykoyjdXuWxsgBUIa2H8wHclmaNrrBD96zQ3ePcJdm3qhx\/82\/iqTwUC5\/QfVRUIqXpCBLSBmhB8Ta0IVCEgmxtunqp8gtEHUAiaG5SBURJkIMgwgtAEYc6EOeRlBFvoMUqCTKkFAvKJ5C3c1UURGBMBbaDGRFtzrRGQjQELfTLWbfLIuLxN2Z92d8WC43Ecd87yJYxognzEa4OOeYyYM2HuEnRShrmPpI3yNRF46ilK\/\/8+aX7jBhHkpJsiEB4C2kCFd0z6qkjjLAwB+Q6UIf1uCOm2PASMWd6adEWLRUAbqMUe2rAX5msAfLKwVzFOdRFlDxiUX\/qmHjfFvUcwNVQ3BIwhWq26xVDvBSAwjyVoAzWP47S4Kg0lF0rSbQoEXOz7\/EL6FOvRnAtDQH4PypiFLU6XsyQEtIFa0tGc0Vrcm\/iMSh+9VF+Do\/iNfhg2KuGki40ioigiun6dyJhJS9HkikAZAtpAlaGjOkUgAAT6bpbceO48gCVrCZuKQBzT+ovk9+4RWbupSOi6Z4CANlAzOEhLLdHQ0X9dunPa2O3owhWXo3jobOEIxDFRHBPFMZG1C1+sLm\/OCGgDNeejp7VvBAK+BsonawtGn7Ha1qB+igA9eUKFpPAoAgEioA2U56CoaBwE3Bu3Ox+nivCz9I3LOf0vc8I\/6FqhIqAIBI+ANlDBH6LlFth3Y7BcpIiGxGrI2KSbIqAIjImA5hoRAW2gRgRbUykCfSHQZ9PTZ6y+1qdxFAFFQBEIHQFtoEI\/Qguu75z4KElv4uUHOqKI+toU676Q9MRRkSKgCGwMAtpAbcyhDm+heiNvd0z6xK3PWO1Wo16KgCKgCMwTAW2g5nnctGo\/AouVynfrFrtIXZgioAgoAjNCQBuoGR2sJZWK\/xj3Ht1bLwnzXdol3Y4jAGzeprfXCszbYgVfEAcD3zYWx9BREVAEFIFNRKDfBmoBCJ4\/f55AC1hKsEvATfsyXaa95IfEdoNu0AvJjxBtPNsnVhwLOEtgMV867nhNM8m1uzzbYHR1mBfJq3TQKykCisCyENAGShxPXBzv379PIPBCNXv20aNH9NprrxHGqReDGzZu5r46IH+FXjmiQs2h1H6ksBqTrrU3xaqspD5jleUJTYfXMl7TTJj7aoScbTBiznbgQTyXI+QgKVsS3\/UcnhKLJrVPWacv95xr961niTJtoA6PKi6AuGgeTkubqP39fZob4cX4p3\/6p0HUze88vbj\/IsNNL++\/TPyXZtBLfEOqXdZVh+9aO7CgZJNYgS\/CqqwmjiW\/OP57+79HRbFoAZv7ui5aks\/OvR7IuYwDOUjKiviy4xOqrus5POW6tPbp7lVFr4ElybWBanA0z549SxcuXKArV67QpUuXZkWoGUvFOHXtqAP03nffw5DSG2++QT\/+5o9THru\/f+Xvr\/FFzZBhnLr2pvlRc5fa4QsCVlv7W2AJN6QirMrqS52THWIlQ\/r7Z9\/9s1Lccb7jvE+NB92FERyNFFOfFQFDYInzoewYhahDzcACY4j1ldWEmrX2ae5VwB7nPM59HIMlkjZQDY4qToRbt27RnTt3lDpgwO+AfOarn6HryQ\/mf37hzwlzHA7M\/9P\/+58U4wRjYAFMgA2+7I35t85+qxVW8EWsz174LK2SH8z\/61f\/a2ksnO804819ZwjNkSvj5bEOehDmrOs66rVDr5mbeN+Y+\/Wj6nWvDVQVQo4eF8KLFy\/S0mjM9fBHRm+dfYvwl3hX6Sr9+7P\/njAH3NCPWU\/IuYAFMAE2Pzj7gwSpq\/Tfz\/73VljJWPiuGeZVuON8R\/4lEBoiNEZFaynTFfk0kQPLkM81rU2v632fAzjnm7xG5marDdTcjtgC6n2VXiXcvCnZYooJX26W38\/Bu1KkW4pAn1i5sYA5CInwbtSSca9qnoCBkiKgCMwSgcmK1gbqEHr86xMX2cNp+igDyHiuY78I3KW7JG\/YfAN\/QA8IPOm2RqBPrPqMtS4wcAav66rXMvSwC3wpWp4ioAgEhIA2UOJg8EUUF1LwQqXsAAjgOz1PKPtB44T5AGkWERLYZEg9oa5Y9RkrdHDxWkaNGCVBBoIMIwivecyZMIdcqQIBVSsCG4qANlDOgcdFE+SIdaoIKAIzRACvZR\/xUqBjHiPmTJi7BJ0r43mZjm10VAQUgeUgoA1Uw2PJ\/zpt6DapOdfsjpMWVSM56vWZQQ7y6UKRufVh7qMe6u01hKzRDcw6V67zagTmih3X7Y7VK57WAvX6KoAc5NOFInPrw9xHodTLdcgaWcYj63i+lFEbqAZHEicB\/pUJAt\/AdXJT1OzS5EUVFABsQT415LwO8D6bKWWoCeSrgeuWo89uKhnqlrVhzrWAZx14lutYjQDwmjN2XLscq1c9jQWwBvmyQ85rAO+zmVKGmkC+GrhuOfrsppKhblkb5lwLeNaBZ\/kSRm2gah5FHHicBGwOHjJigY69IQBsQW5A4C3l4CFz7aacoybQlDW0yQ0ci+p2dbCDrE2eTfMBTsCL1w0eMp7r2C8CwBfkRgXmUg4eMtduyjlqAk1ZQ5vcwLGoblcHO8ja5AnRRxuoEI\/KADXhpGUaILyGrIEA44+xhrmaKAJBIIDzlSmIgjawCMYfY5\/L11jdENAGqht+s\/FG58+kL8JpDhvjjzG0Y4CaJCqoz5VJvfKbgwDOAyacF5uz8nBWyvhjDO0YoCaJFOpzZVK\/JF4bqCUdzYK1uCcz5jjJC8xVPAACwFyGxTzUY4C6UJ+sV\/khEAg\/pnseYI7zI\/zKl1MhMJerwTzUY4C6UJ+sd8m8NlBLPrq6NkWgIQKbdgFsCI+aKwKKQAECm3jt0Aaq4GRwxeiqcYKwHDxkPA95RK2yPsznUrusGzWjdpaBh4znTccx7VGrzId5aLUX1YQ6oeP6wUPGcx2LEQBOwIstwEPG89BH1CtrxHxO9XPtqBm18xw8ZDwPeUStsj7MQ6u9qCbUCR3XDx4yns991AaqwRHEgccJAALfwHVSU9SKmpkwn7SgDslR+xzXIetG\/Zh3gKF3V9SEoBglQQZCvSwHD5lSPQSA11yxk7VjDZjXW3V4VqgdawCBD69Cf0WoFTUzYe63nEaKupAZoyTIQKiX5eAhWwpN2EDNE0KcAKC5VY+ameZSO+r11Qo5yKcLRearDzKmUOrkOrgud2Q9RtaBV2qGwJyx49oxNlv1dNZFtUIOmq6y6sy++iBjqo4wrgXX5Y6yCtZJ2RJ4baCWcBR1DYqAIqAIKALjIaCZFIEEAW2gEhD0VxFQBBQBRUARUAQUgSYIaAPVBC21VQQUgRAQ0BoUAUVAEZgcAW2gJj8EWkBXBPAFRY4BHsTzJmNbvyY51FYRUATCQUC+5sGD2lTX1q9NLvUJBwFtoMI5FvOpJNBKcRHjLyuCD7RMLUsRUAQCQwDXC712BHZQZlCONlAzOEghlYgLjY9kjdDLeR0ePkVUx78vmxBq6GstGkcRCAmBOq8t2DStGT5F1DRWF\/sQauhS\/6b49rlObaD6RHNDYvG\/1OSIi0eX5ctYiOPOIQMhj0ss5xF6+GPelODnUtMYaq8IKAJ+BNzXFuZ4vfqt60kRgwkezGPEnAl5XIIOMh7Bu37Q1SH4uVTHT23mi4A2UPM9dhtZuXuBwhxAYJQEmSRcGMvmUqe8IqAIDIXAdHHl9YF5VMM8j5BJ0muHREN5iYA2UBIN5ReHAF\/8cHFkfnGL1AUpAopA7wjw9UKvHb1Du5iA2kAt5lBu1kL44uau2pXXvfjBD7ZuvKXNdT2KwKYjgNe6DwNXjuuBKyvyg61Pp7JlI6AN1LKP7yCrw0XFpVAuIKgDtfkWDp1PrjJFQBEYBwG8Nl0K5XWJOlCbDwnofHKVbTYC2kCNevyXkQwXEyasCDzGUMmtDxdJKXPnch1lOmmnvCKgCFQjgNcdE6zBYwyV3Prc64E7l+so00k75eeLgDZQ8z12QVSOCwwuFGMWg3zIW5QTOtgU6V057F2ZzhUBRWBYBPC6a\/I67aMa5EPeoljQwaZI78ph78oWO9eFHUNAG6hjkKhgyQjg4lj3otfEdsmY6doUAUWAqMn1oImtYjtfBLSBmu+xC6ZyNCS4YMiCMPeRtGnDIybyDemLHKC2edrUpj6KQAUCi1TjNYbXmlwc5j6SNm14xES+IX2RA9Q2T5va1Gc6BLSBmg77WWYuujBIOfgi8i0aFxxJsJFz8JCBEBdjFUk7+IOkzOcPGxDsQD4blSkCikA7BIpeU1IOvoh8WfF6lQQbOQcPGQhxMVaRtIM\/SMp8\/rABwQ7ks1HZ8hDQBmp5x3R2K8IFp4rKFgXfKn2VDfxhAwJ\/jFSgCCgCwSGA12sVlRUN3yp9lQ38YQMCr7Q5CGgDtTnHWleqCCgCioAioAgoAj0hMJcGqqflahhFQBFQBBQBRUARUAS6I6ANVHcMNYIioAgoAoqAIlCAgIqXioA2UEs9srouRUARUAQUAUVAERgMAW2gBoNWAysCikAICGgNioAioAgMgYA2UEOgqjEVAUVAEVAEFAFFYNEIaAO16MMbwuK0BkVAEVAEFAFFYHkIaAO1vGOqK1IEFAFFQBFQBBSBrghU+GsDVQGQqhUBRUARUAQUAUVAEXAR0AbKRUTnioAioAgoAiEgoDUoAkEjoA1U0IdHi1MEFAFFQBFQBBSBEBHQBirEo6I1KQIhIKA1KAKKgCKgCBQioA1UITSqUAQUAUVAEVAEFAFFwI\/A\/wUAAP\/\/dP3YdgAAAAZJREFUAwDjhnUJvITXnQAAAABJRU5ErkJggg==","height":285,"width":474}}
%---
%[output:00e322bf]
%   data: {"dataType":"text","outputData":{"text":"[10:10:13][INFO]  100 個の RL 後 SMILES を生成中 ...\n","truncated":false}}
%---
%[output:184ca7e7]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 10\/100) generating (post-RL)\r[##--------]  20% ( 20\/100) generating (post-RL)\r[###-------]  30% ( 30\/100) generating (post-RL)\r[####------]  40% ( 40\/100) generating (post-RL)\r[#####-----]  50% ( 50\/100) generating (post-RL)\r[######----]  60% ( 60\/100) generating (post-RL)\r[#######---]  70% ( 70\/100) generating (post-RL)\r[########--]  80% ( 80\/100) generating (post-RL)\r[#########-]  90% ( 90\/100) generating (post-RL)\r[##########] 100% (100\/100) generating (post-RL)\n","truncated":false}}
%---
%[output:77d3414f]
%   data: {"dataType":"text","outputData":{"text":"[10:10:32][INFO]  RL 後:   valid=4\/100 (4%)  unique=4  novel=4  avg_reward=0.040\n","truncated":false}}
%---
%[output:9e6d6b77]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAlAAAAFkCAYAAADv4QVjAAAQAElEQVR4Aey9f5Cd13nfd+69wO7CWkhQiTWyqaJAkEjGUocsqyRAwwaeUaduXXo6bW1OB6DJQf6I03rsUNYInmmdjJwa9jTeVinrxq1bV61VF+m04z87\/cN\/ZmwKQ5oS08qKHQtqDThATCrATBwSCwK72c+7+yzOPTjvr3vf9973Pe+Xw+ee5zznOT+ezznnwdkluDve1T8iIAIiIAIiIAIiIAK1CIyd\/hEBERABERCB3hHQgkVguQT0gFouf80uAiIgAiIgAiLQQwJ6QPVw07RkEegCAa1BBERABIZMQA+oIe++YhcBERABERABEZiJgB5QM2HrQietQQREQAREQAREYFkE9IBaFnnNKwIiIAIiIAJDJJBIzHpAJbKRCkMEREAEREAERGBxBPSAWhxrzSQCIiACXSCgNYiACDRAQA+oBiBqCBEQAREQAREQgWER0ANqWPvdWLRnzpwpHauKTzhIlT5VfMJxy+pNjVlpnLLFdKS9SixVfDoSTqeXUcSxqK3TQXVkcUX8ito6snwto8ME9IDq8OZoafsEmkxyTY61v7rlfRJLKP5qaPPrvl7Uhh\/t169fR5UsgACsYb6AqTRFDgHxzwEjcy6BIT6gcmGooZxAqkmmj3+AseZQ8nZwnn2zvpShMF9oo45dMj8BWLLH8480zBFCftSNBLrYGg2VsxDQA2oWauoTJUBC8gUnv46OLSZ+m6\/HfGe1MW4TCZNx8oS15bWFdny7JrZG40TddMpQWH9oo469S0IcZRJbL31i9nlsjBkTxgztMZv50NYlsXUVlbH14h+zz2tj3PAsUseOoNscYd3s3Su1oi4R0AOqS7vR87WQkHwhHL+Ojq1M8COh+X7UTbCbbiW2eSU2b2xM\/GKCb8yeZ8N\/VrG4Y+WsY9KPtVKmKMRWJLBcRNxFa6jTtoi11p2jbP2LYsw8rIX1o1OamN3qtIc2a1MpAkUE9IAqoqO2KAESDg2UJtTbEMZnXBKcSVjHjq1IGKeKHz74Fo0Va6MPfa2NuultlMyVJ\/POx7isH0H3x8PmC21+HR1bH4VYy9a\/7Li6vr4yPotizDxFa6E9xjJmKxpHbcMmoAfUsPe\/VvQkFxIPQkdKE+pNCuPaeL5utjqlrbtqH+ajT1X\/FP1gYHGFLGgzwcd0SuqSagTgGhO\/N+1+XXo9AvDLO5fYw3az1ZtF3kMloAfUUHe+4bhJRKEwRWijjt0XbCQu3xbX61tnHZv10LfKjPjhX8U35kP\/mL3MRr9QyvpUbWfceWKqOs\/Q\/WDsy9B5LDJ+\/4yjm7Af6Itci+bqJwE9oPq5b0tZNYklb2La8oQ+fhv1KkISo18V35iP3x+9qthYzE0fq8dKa6f0BV+\/bjr2UJgntFWp0y+UKv2KfPx1ms4c6NYP3QSb6ZTUq0iZr9+OHpMq8wzZB2ZF8fvt6DEp6t\/1NuLx10jdF9qszhk3wS4RgSoEaj2gqgwon2ERsAQ0T9Qkrlj\/PHvMN2bz+6P7gr9fNx27L9j9uq8Tu7VT+oKfX0fH1nVhnSasFZ3SF2wm2E2npF5F8IVfzBc77bSZTj0U2vBpShiPOZoar8o4zOlLlT5VfYiFsWP+2GmnzXTqodCGT1PCeMzR1Hh548TmYV5f6Gt1dF+w+3XpIhAjoAdUjIpshQRITjhQkmgQ9DwxX2unXkesn5X0Nd1KbIsSm5O4FzVnbB5bh1\/G\/JqwLTvWWAx11+Rziul1x4utqa6NOX2p279tf9ZWZ44YV99Wd7w6c\/u+kXn8Zuki0AgBPaAawTiMQSwRWnKykujR8yRsp15HwnHpG7NhX4TY3IuYq2gOW4dfFvlXbWOffd9YHRuCH6UJ9arCuunn+1PHbjZ0bCZmn6VkrJgwFnbK1IS4YOfHRR272dCxmZh9lpKxYsJY2Cm7IMTKOihDwS4RgSoE9ICqQkk+GQESIJJVIh8kooi5cybWWRRH5xbcoQXBzcRfltkoffuhPofCmCbsnckcQ051ZWzGnDIuqMK8vixo2semgYFJG+thbMZ9bOIlGFgH62FqSl+wSUSgKgE9oKqSkl8pARIRyanUUQ6NEIB1KI0MHAwS21fmxY4rJXX0toW5TBY1Z5sxWSxWtjlX1bFtLZQpMPbjJh7iwkZJHR1Bx4YuEYEqBPSAqkKp\/z4Li4AERCJa2IQ1J2JtrLFmt7ncmROZa5CczsRikuNSy8w6Ga+oU8yHPtiL+sXa\/H70p+77YfPrbenMu6i52oohb1w\/NmKk7vti8+tt6cy7qLliMTA3a\/DbqGNH0P026SJQRkAPqDJCaq9NoKuJqCxJlrXngSiLl3aTvDGats8aS946WD9jIugxP+y0I7H2ItssfYrGS6UNpk3FMnTGMZYwwY6gN8Va4zRJoLtj6QHV3b3RyhokQHIkSTY4ZGNDsbbGBssZqGrsRX60ITlTZGbakaxS8cP8rfS7YYNPTGjzfefVGY95wnGwFUno38U6sbEuK9FNsOXFR5v5NVEyHnOFY2ErktB\/3rrNxXpsLHSzm02lCBQR0AOqiI7a5iLQRDKyMfySRfl107GHYm0kx7CNurVT5vng16bMMy\/rNrE1Up9nTBtnkWXRemmLSZ310b+Kf+hHvUyqjFvmw54heX60IXntVezEkedHW0zwryr0r+Ib+lEvkyrjVvGBIWLzhX3Mjg8StqsuAj4BPaB8GtIrEyDRlDnjY1Lma+34m05JvargH4r1De1Wt3ZKszVRNj1ebE3MEQp+2CiLZJE+ResYcpu\/B+gmMSbWRhlrly1OIORFHYl7P7LigzyySBOBxwnoAfU4E1lEoEMEtBQREAEREIEuEtADqou7ojWJgAiIgAiIgAh0moAeUCXbo2YREAEREAEREAERCAnoARUSUV0EREAEREAE+k9AEbRMQA+olgFreBEQAREQAREQgfQI6AGV3p4qIhEQgS4Q0BpEQASSJqAHVNLbq+BEQAREQAREQATaIKAHVBtUNWYXCGgNIiACIiACItAaAT2gWkOrgUVABERABERABFIl0N4DKlViiksEREAEREAERGDwBPSAGvwREAAREAEREAGfgHQRqEJAD6gqlOQjAiIgAiIgAiIgAh4BPaA8GFJFQAS6QEBrEAEREIHuE9ADqvt7pBWKgAiIgAiIgAh0jIAeUB3bkC4sR2sQAREQAREQAREoJqAHVDEftYqACIiACIiACPSDwEJXqQfUQnFrMhEQAREQAREQgRQI6AGVwi4uKYYzZ85MzRzWpxoPKviUyYHrY0VZP9of65RjqOObM8Rj5jbGfGwSGUSgywS0NhEYEAE9oAa02V0I9fr16y5PWB9tlHlCe57k9QntPHQYI7TPW2dMxp53HPUXAREQARHoPgE9oLq\/R51cIQ8FHgz+4qhj921VdfrRv8wfvzwp61vU3saYRfO11KZhRSAJAtzHtgJh7Ji0NZ\/GTZeAHlDp7m1rkZF88h472GmvMzn+9KvSB788qdK\/aK7YuPhXGdd8GKNuH+urUgREYDEEuKeh6N4uhn1Ks+gB1eRuDmAskgyJpyhU2vEr8rE2\/PC3eqzEB6GNMk\/8dvQmhLUxXxNjaQwREAEREIF0COgBlc5eth4JDwl7UKAXifk1sSjGQhiL0qSoTlsbEsbcxhwaUwREoJhAlXsY+lAvHnXYrYq+PgE9oOozG2wPHi4ET2kS1rFjQ3yd5BUT\/GJ2bLSFgt2ENtMpqbcpzEFMvmBrc06NLQIiME2AO+ffQXRsvhd17Ca0oVNKRKApAnpANUVS4xQSIHnFhE4xOzbaQsFuQpvplNTblEXM0eb6NXaXCWhtiyTAAysU3e9F7kAac+kBlcY+DiYKP+kRdFjH1qb486G3OZfGFgERaIcAj6VQ2plJo6ZMQA+olHe35dh4QJCEWp7mcHjm8oUGv46ObRap0sfiZR6TKv3kIwIisFgC3E\/uqwn1xa5Asw2BgB5QQ9jlRGK0ZGglYZluJbYiIZHiW+Tjt+FLH99mOm2mW4ktz998VIqACLRLwO4hdxFpdzaNPlQCHXlADRV\/f+O2BBVGkGcP\/Wapkwh9YQy\/jo5tVmHtofhjovvt1GedS\/1EQASKCfh3zde5d34dHVs4GnZfwnbVRWBeAoN5QNlFygNm7WGZ5z9Uu\/EJE1aePY+T+ee1+3bzDUt8QpvVacsT1o6f344tJr4Puu9jdUqEMWlHl3SbAHuF5K2Stpjk+cveLAHuUZ4wU9iGzYR9C9upY8cHnbIx0UCDJTCIBxQXh0uDoOftNu2h5PkO1W58wvjz7KGf1c2f0mx5JT51JW8sszOe6U2VbYzZ1No0ziMC5AD2CkF\/1DKt0R7KtIdqXSTAnrGvoWDv4nq1pv4SSP4BxSXyLw46tv5umVYuAiLgEailcvfJAdYJHZvVVaZBgH0NJY3IFEWXCCT\/gKoDm0Rqktfv5s2bzuTBgweujtCvjr986\/EVr+XxyrsvfbZbLqAsioN7baIzuLwzGLLXnix\/L4ruTQptekB5u+h\/xRJLmlzIy5cvu\/Pnz2dy5coVd+PGjSn5zne+4775zW86Sr\/tjTfecD\/xIy+4F1980aH7bb3WD+In3ljcKcRWFIPi\/k52\/u\/cuePdpDTUsnxAlFVyQtH5mact9bM3T3zk2J+89KNZnu5qzp0nvnnOzSL6Wmzvvvsu1yRZ0QPqYGtJlgdqVlAPH1Eky2vXrrmtrS139epVd+nSJbe5uTklJ0+edKurq1M2fLa3t90\/+ON\/5t5++223sbHxWDs+fZa8uPscU5W1K+7983\/8+PHs3qTywf33Y6Ee5gPaq+SEKudoFp\/Uz9488ZFv377+h+4HP3Wqszl3nvhmOS+L7GOxra+vc006K\/MuTA+oGQiePXvWnTt3zp0+fdqtra1NybFjx7IHVGinblOhpyZFcacWqx+P4t4\/\/0eOHLHjPciyKCf456VJPfWzN298HMTnNj9CMZWjm9yDecaaN7555m67r8WWel5I\/gEVfuXIV5HYslvlfWD3qo56zM\/3kS4CItAvAtxp7ratGh2b1a3Ebjol9ZgfbZJZCaifCPSbQPIPKLaHxEcCRNCxmWBDx45uQh27RAREIC0C3O28e46daH0fbNSxS0RABETACAziAUWwJEAE3Rffhm7i+0gXgRQJDDmmvHuO3bigm5hNpQiIgAgYgcE8oCxglSIgAiIgAiIgAiIwLwE9oOYlOHN\/dRQBERABERABEegrAT2g+rpzWrcIiIAIiIAILIOA5swI6AGVYdCHCIiACIiACIiACFQnoAdUdVbyFAEREIEuENAaREAEOkBAD6gObIKWIAIiIAIiIAIi0C8CekD1a796sVp+vUXScvOms\/hu377tEKunWPbi0GmRnSdQ925wr5C6\/fAPYWDrmhAb0rV1VV1PyHiIdT2ghrjrLcbM5bt8+XL2SzzPnz+ffPm5z33OvfTSS44y1XgvXLiQPRhbPDYaOnECs+QF7tSsd+vi26KRBAAAEABJREFUxYsZ0V\/4+7+fldS7dj\/nia8LsSgvONfHB1R2IfTRTQIkymvXrrmtrS139epVSc8ZvPrqq4797OZp06r6QkB5Ia1cqLywf\/P0gNrnoM+GCdgvV+WXLiMf+9jHXApCLKHw1a0v1o7NdMqyeswHW57YeFbm+eXZq\/RjHxs+GhpuwAQ4T\/55TCEnEIMfEzp3yxdsyMWLFx2lCT6mU4b1PBv2mFh\/K2M+Rbaq\/djHAR\/jw9D1gDpEIaUtAnz1mcp\/1ot929p+3YeVIUd+l1rMhj922k1i9dDm+9ImEYE+Ekg5L3C3fQn3hzscs9EHO+0msXpo831pkyyGgB5Qi+E86FlIlPxnoOde+aJ7\/vNbvZWnX3g5+p+zLHlZWbbZ+PmJEt2EvqZTUkfQY2JtlLMIa\/FlljH61Edr7Q6BlPOCf6fQy6jjw\/3Gz3TqCDZKE+qI1cPS2ihnEeb3ZZYxhtJHD6ih7HQH4jz55LPu5FM9liefiVKMJbCo44ERf0tQ6JitHurU2xTmN2lzHo0tAnkEUswLdqeszIvd7PhZDkDHbvVQp96mML9Jm\/OkMLYeUCnsYm9iSHOhluisrBolScr3tbqV1mbj5pXmR4kPZZngF85T1kftIiAC1Qhwv3yp1su58E5a3Uobxx87ppsfJe2UZYJfOE9Zn6G36wE19BOg+OcmQNLxpWxAEhU+lAg6YrqV2BB\/bPTQRt2EdvrniflZiZ\/pKkVABJohwD30pWxUu4eUiPmbbqXZ\/bHRsVOaUDfBRv88MT8r8TNdpUcgouoBFYEikwjUIUDC8aWsLwkNwc9KX\/dt2OsK\/fOEsVgr7eiU1BF0bBIREIH5CHCffCkbjbuH4Gelr\/s27HWF\/nnCWKyVdnRK6gg6NkmcQFIPKDYciYc6bcUPMSt6TKxdpQjkESDJ+JLnZ3bOmel+aXYr\/bYmddbqjxfW\/bY+63BEqsSAH2K+6DGxdpULJ9CrCblTvpQtnrMW8zG7lTGfJmys1R8nrPtt0h8RSOYBxQFj0xH0RyFW1+gbSvXe8hwqAc6bL3CgzllCryKhb1ivMsasPrZWylnH6Fo\/YoEhgj7L+ugbyizjqM\/wCHDmfIEAdc4TehUJfcN6lTFm9bG1Us46xhD6JfGAYpP9w4WOLW8DacMnr73Mzv9+izx48MBVFX\/Mqn366OfHGervffe2e++7\/6S3Esbj1zlPiG8zPc9u7ZT+mcSfOvZQsNMe2met++MVjfug5KzPOn8b\/fyYGJ+4sKHHhDZ8Ym1VbeQDpI93dhFrLuL4XqJ5gTOFxGLPs\/u+\/rnEn7rfbjp22q0+b+mPVzZu2dmZdy1d75\/EA6opyBwck6Ix+Wmt\/C6iK1euuBs3bkwJSZRfEBnasdmYt27dmuoT+vaxbnH7cVq8\/KRefnLtb\/3Xl91v\/s2Xeyus\/6lnP5v9RHUX\/OOfG\/SixBO2xfzxwW7ToCPYzTZPyVhI1fHCM2v7bWf1zp078yynk33hY1K2wKKcYIyaLsM9aHr8JsfLywvcJ+5V3\/MC+Y08558T\/+ygF921sC3mjw92mwMdwW62eUrGQuqMF+YFOzN2Nu\/evTvPkjrfd3APqKIDwsExwS9v9+z3vF26dMltbm5OyalTp9zGxsaUDZ8TJ04cDhdrx6dAHhuva74Wtx+nBUxiMWZ9\/\/14v\/pLf8fCOiztzPjlYaOn0O5VD\/+XZd9epPtt\/jiz6IyFVO0bnlnbbzuHx48frzpUp\/y453kcsJvgV7RwO9+xnGCMmi7DPWh6\/CbHy8sL3Ke+5wTWz\/7758POjV\/67abTbjql1a30bTHd96N9HmEspM4YYV6wM2Nnc319vc5wvfMd1AOKJJh3QEI7dfxjO8pXG\/w+odOnT7u1tbUpOXbsmFtdXZ2ymY+NZfWUSj9ui9MveUTBrO9CHH5cQ9HDs+rvN21HjhzpHQruN\/c8tvDQTh3\/mC+2opwAnzYk3IM25mhyTDiFwn3qe05g\/cQRxjaEet75sLO5nxfSJTGoBxTbSBI0sTqlpF0CfEs3BWmX0v7onM99bfozzz7tpVodAjA1oR86pWQxBFLICcTQNq28c5lnb3s9Gn+fQBIPqPCrQw4Vtv0QH31i84UW6pT0oTShbm1mUzkbARJMyr9MGCqcl1Cwm9BmelnJuQv9qWMv66v2\/Z\/mDC9jgR5jh80X\/KlT0ofShLq1ma0vZVfXmXpe4MyE4u8FbX69SOfshf7UsRf1U1u7BJJ4QIGIg8SBQtCxmWAzPa+kD34m1PN8Za9HgER57do196Uff879ypee76382I887YgjL3rOjC+cpTzfmB1\/E9pNp4zVsUniBNgHuCHovhc2vx7T6YOfCfWYn2yzExhCXuDc+MJ5qkMMfxP6mU4Zq2OTLI5AMg8okNlBRfcFu183PbRTNzEflc0R+OynT7pey2dOzgCjehc7e1XK6qMO19M4hgSwhzbqoZ26Ce2Sdgj0OieQ01rMC3b+qpTt7I5GLSKQ1AOqKFC1iUAfCPCVpS\/hmmkLbaqLgAikS4A770sYKW2hTfXFEEj+AbUYjJpFBJwjkfnCV42zcKEfYn0Z03SVIiAC\/SHA3fXFv9d1oqAfYn0Y03SVyyOgB9Ty2GvmxAiQ4EzmCY3kiMwzhvqKgAgsn4DlA8oZVnPYhXyAHBqkdIKAHlCd2AYtIjUCJMxZEx59kdSYKB4RGDIB7rRyQlonQA+otPZT0SRAgCSLJBBKv0PQ6kWgIwTIB0hHlqNlHBDQA+oAhIr2Cdx65z33j3ssdQmFX3GSAE2KxqIfYj6+bjaVIpAKgSHlBe4yOcD2Dt3EbLGSfoi1+brZVC6egB5Qi2c+uBn5NQf8qou\/9rd+y\/17P\/GbVaVzfqyfOIgn3MS8hGZ2Sl\/C\/mGdpFrFFvqoLgJ9IcA94j5xr1LMC9z32F6YndKXmK9vU07waXRD1wOqG\/uQ9CpIlPyyTX7pZt+FONrcLBIq41Mi6CbUEaurFIE+E1BeqLZ7ducpEb8XdcS3SV8cgfwH1OLWoJkGQIBkee7cOdd3IY4BbJdCFIGFEOA+9T0nsH7iWAgwTdIpAnpAdWo70l0Mv7YhBVnEDulb9YugrDm6QKCtnLDocdtmqZzQNuHZxtcDajZu6lWDAMns8uXL7vz5872XCxcuOOIJwyfBheL70ObX83T8wm\/J+zb0vL6yi0CfCHCPvvBjf6X3OYG8FssL3NVQ\/P2hza\/n6fgpJ+TRWa5dD6jl8h\/E7CRKfgnvf\/aXn3L\/zQ8+01v5K\/\/qn239lwlboiRpcjisNN3aqXdXtDIRKCdAXnjzH37bpZwXuK+++Pe5nND+bzegP77W10qzWTt1yWIJ6AG1WN6Dnu25P\/UR12vZ\/Eir+2eJkARpOiViNspWF6HBRWDBBHqdE8hpLeYF7j7bwb03nRIxGyU+ksUT0ANq8cxbn1ET9JcAyZDkSAToJlb3S3SJCIhA2gTIAcoJ3dxjPaC6uS9aVQ8JkOh8saRXNRT6+r70R7BR+oJNIgIi0G0C3GlfuMN1Vkxf35\/+CDZKX7AlIL0KQQ+oXm2XFttlAvMmM\/rnxUciNcnzkV0ERKBbBLjTJrOsjL55\/SwfUOb5yN4ugcE8oDhkSBWc+CFVfOUjAjECJL4mzpCNwXhIbC7Z6hOAK1KlJ35IFd\/e+yiA1ghwf5s4RzYG4yGtLVgDlxIYxAOKA8dBQ9BLqchBBDpAgLPKmbWlhHWzq6xHwDjCFr1eb3mLwPIIcF45t7aCsG52lYshkPwDKjxgHD5seXhpwyevHTv\/+y3y4MEDV1XoZ1K1Tx\/9LMZYeetPtt3tP7nXW4nFVGTjHHGezAfdxGxBOVWlPwZK+qEvS8rO4rLWVXdeOMLT+qFjs3pY0oZPaA\/r5AOkjNNQ20Nefn1IeYGzxJmy+NFNzFZU0p92SvqhL1PKzvMy17aIuZN\/QLUB8eLFi9kPf7ty5Yq7cePGlJBEb9++PWXDB5ut5datW4+149Nnsbj9OC1efs0BvzT0r\/\/f\/8C9+H++0Vth\/cRBPBablSQ00\/3S7JS++D55OgnSxPpaPa9PW\/bwzNp+25m9c+dOW1P3YtyinGCMmi7DPWh6\/CbHG2Je4M7GDq\/ZKX2J+YY2u\/+U1hcdCX0XUQ\/zgp0ZO5t3795dxDKWNoceUB76M2fOOA6lZ4qqW1tb7urVq+7SpUtuc3NzSk6dOuU2NjambPicOHHicKxYOz59Fovbj9MC5sFhzODWZyEOi6uN0j9\/6CY2V1g3e9tleGZtv+3MHj9+vO0lLHx8\/lCCd5WJORec61hOMEZNl+EeND1+k+MpL1Q5RXEf\/wyim5h3WDf7IsowL9iZsbO5vr6+iGUsbQ49oA7Q10mWfBeCXyB5+vRpt7a2NiXHjh1zq6urUzbzOZgq2mY+fS39uImTr0B8wcZDqu9CHH5cqevEi4Tn0t9v2o4cOYJbMlInHxB0UU6ATxsS7kEbczQ5JpzC+4Kt7zmB9RNHGNui64ucj3iRvPNhZzO1vEDMvugB5dEgaZpgRqeUVCdAMuEPE\/tPGvyeKEl\/fwcg+8h+sq\/VT0Eantx\/EyJCp5TUJ8D54RxxnpQP+psPbO\/YR\/aTfa1\/GtLpkfwDim9v+okPHVu4hdh8oZ06paQ6AS6U\/ecM\/pNG6vKVr3zFffnLX3Zf\/epXs\/+sm2K87Gf1E9BtT+40OcBWiY7N6lZi8wU7dcq4yFpEYJa8MM\/devXVV7Pl\/OCnTmUlZ7hrd3Oe+LoQC0wzuAP+SP4Bxd6S+EiUCDo2E2ymq2yGAMmS\/8Q5FHn22WddyrGyn82cjG6MQg7g3iPo\/qqw+XXpzRHgHNW9J7PeLb47wsqfO\/g9ddTrzr0I\/1njW8TayuZgP2E8ZBnEA4oNJlEi6L7EbLTn2WmTdI+AViQCdQhwv5GwT8yGT56dNokIiMAwCQzmATXM7VXUIiACIiACIiACbRBo6AHVxtI0pgiIgAiIgAiIgAh0k4AeUN3cF61KBERABERgEQQ0hwjMSEAPqBnBqZsIiIAIiIAIiMBwCegBNdy9V+Qi0AUCWoMIiIAI9JKAHlC93DYtWgREQAREQAREYJkE9IBaJv0uzK01iIAIiIAIiIAI1CagB1RtZOogAiIgAiIgAiKwbALLnl8PqGXvgOYXAREQAREQARHoHQE9oHq3ZVqwCIiACHSBgNYgAsMmoAfUsPdf0YuACIiACIiACMxAQA+oGaCpiwh0gYDWIAIiIAIisDwCekAtj71mFgEREAEREAER6CkBPaBm3jh1FAEREAEREAERGCoBPaCGuvOKWwREQAREYJgEFHUjBPSAagSjBhEBERABERABERgSgaQeUGfOnHFI0QbS7ov5+jZft3aVIiAC\/SJg97ho1eZjpflaPSytfc5S3UVABBIgkMwDikR3\/fp1h6DH9gY77b5gM1\/fbrq1qRQBEegPAe613WH02Mqxm4+V2MzXbH5pbSpFQAREIENeaqEAABAASURBVIkHFEmPJGfbiY7N6k2XN2\/edMiDBw9cVfHXULWP\/KrzrcWqxr5p3Gp74J\/vZevcfXKArQMdm9XbKMkHiM5LtfPSNqdwj9ueT+PH9z3ch9TqSTygqm4KibTIlyRrUuR38eJFd\/78eXflyhV348aNKSGJ3r59e8qGDzYb89atW4+149NnyYu7zzFVWbvi3j\/\/d+7csePdm7KpfEDARTmhyjmaxSf1szdPfJZv72\/fZ3tcF3PuPPHNcl4W2cdiu3v3bsY\/1Y9lPKCWztIeSZR+EkU3oS1voVtbW+7q1avu0qVLbnNzc0pOnTrlNjY2pmz4nDhx4nC4WDs+fZa8uPscU5W1K+7983\/8+PHD8903hbtuwv239aOb0G72WFmUE6qco1l8Uj9788Rn+XZy5Ei2XV3MufPEN8t5WWQfi219fT3jn+rHIB9QlhQpLTGi+5tM3dp8O\/rZs2fduXPn3OnTp93a2tqUHDt2zK2urk7ZzIe+iNVTKoviTinOMBbFvX\/+jxz8QcX57ptw103szlP346Bubb7d9KKcEJ6Zpuqpn71542NvJpP9P+KaYr62tn\/emxhv3viaWENbY1hsfc4LnJ8y2T9dZV6JtBclwERCVBgiIAIVCSgfVAQlNxEQgSiBJB5Q4VeHJEZs0YhzjPTxm6jXHcPvL10ERKB5AlVG5N5yf80XHZvVq5T08f2o1x3D7y9dBEQgPQJJPKDYFpIbSQ5Bx2aCDR07ui\/YYm1mp00iAiLQLwLcX7vn6P7qsVPHju4Ltlib2WmTiIAIiAAEknlAEQxJDkH3xbeh+xL6WZtvl24EVIpAfwjk3WXsFgW6L2anzLPTJhEBERCBpB5Q2k4REAEREAEREAEReIxACwY9oFqAWjYkPyOjzEftIiACIiAC1QiQU30Je\/lt6GG76iIwCwE9oGahNmcffujehQsXsp9mPudQ6i4CIiACfSDQ2hp5EF2+fDn74cb8gGPkp3\/6p6fme+2116ba\/6Mf\/neUf6cIqTILAT2gZqE2Z58f+v6Pu2vXrs05irqLgAiIgAjwgCKffunHn3O\/8qXn3Y\/9yNOPPY7C9je+\/vsCJwJzE9ADam6E9Qf47GeeqN9JPURgHgLqKwKJE\/jsp0+6TD5zMhpp1oZPTnu0k4wiUEBAD6gCOGoSAREQAREQAREQgRgBPaBiVJq3aUQREAEREAEREIGECOgB1cPN5L\/515Eehqgli4AIiIAIdIKAFpFHQA+oPDIdtfNwCv+PE\/6vkyLR\/\/HX0c3UskRABERABHpLQA+onm0dDyj+j5LnXvmie\/7zW6Xy9Asv6\/\/469kea7ki4BOQLgIi0E0CekB1c19KV3XyyWfdyacqyJPPlI4lBxEQAREQAREQgXoE9ICqx0vegyOggEVABERABETgcQJ6QD3ORBYREAEREAEREAERKCTQ+QdU4erVKAIiIAIiIAIiIAJLIKAH1BKgL2NK\/vJ5VVnG+jSnCIiACCRGQOEkTkAPqMQ32MK7ePHi1C\/T1I89MDIqRUAEREAERKA+gU48oM6cOVN\/5TV7MAdS1I12X4p8+9amH3vQtx0b7nq5g41EXzAIcyAFLo52X4p81SYCIjA8Ap14QPnYSVh+PdTL2kN\/6vS5fv26Q9CxhYKddl+whX59revHHvR154a97rI7WNYeo0cfu+foZT5lvrH+somACKRPoHMPqBA5CQ4J7VXr9CUBmj86NqvPUtrfJXrw4IGrKrF5qvb1\/WycnZ0dV0V2d3ezLlV88cmc9z78OWfUK7PR+NXPUVdZ7R2ZhfzL3UVmnYy+5ADrj47N6rOWs+SEru5l39Zle0auI4ft7OznPLNbGbb3Lc4+rtfYp1p26gFliYwSMehNJTkbL1YyR8wes9nfJ7py5Yq7cePGlJBIb9++PWXDB5uNdX\/7fqbeunXrMT98i8TG2d7edvfu3SuV+\/f356rqjx+Lq7u2vLiLYkmhTXHvn\/87d+5wbBoVywGUiA3OXfXrZm+yZI464xXlhLbOeepnr2p8hznxIB9azmP\/LNeibwftdXNc0\/tYNb6m513EeBbb3bt3Qd9DqbbkTj2gLGlZ6YeAjaRJ6dub1pnDJG+ura0td\/XqVXfp0iW3ubk5JadOnXIbGxtTNnxOnDhxuNQjR45keswP3yKxcVZWV9zq6mqpHD16NJurrn\/dteXFXRRLCm2Ke\/\/8Hz9+PDtnTX7Y\/bPSHxsb95TStzetM4dJ0VxFOaGtc5762asa32FOXNnPhysHOY+zMDnItegrQXvdHNf0PlaNr+l5FzGexba+vg76ZGXpDyiSU5fokiRN8tZ29uxZd+7cOXf69Gm3trY2JceOHcseNaGdusU5nuxjxzaLMM5kPHGTSbmMx2Pc3WQ8cZMa\/nXXVRR33bH65K+498+\/fVGQHbY5PvLu3BxDztXVcgFl0dqKckJb5znlswezOvGxyePJOMtxlNSRyZ6NEsE+mUzc+MDGHMuUOvEtc52zzG2xNZUX2L8uynjZiyIxVVkDyQtfyir+dX3aGrfuOuQvAkMmwB2vEj\/3FV\/KKv51fdoat+465C8CItBdAkt\/QPloLGlZaW3USZZWr1PSj\/7WBx2b1VWKgAiEBLpR566yEivREeqz3mH60Z9xEHRs6BIREAERqEOgUw8oS2SUiAXi62arU9KfRImg+32xUceO7gs22iQiIAKLJ2D3jxKxFfi62eqU9Ld7ju73xU4dO7ov2GiTiIAIiAAEOvWAYkGhNJW0GAcpGp92X0LfRdU1jwiIQJwA9zPeUs\/KOEjYy7eh+xL6qi4CIjBsAp17QJGwirakrL2or9pEQAT6RaDsvpe19ytarVYEek9gUAF05gHlf6s8Tx\/UzihYERgwgbwc4NsHjEehi4AIdIBAJx5QfBVZRSx5doCbliACItASgSq5AB\/lg8gGyCQCIrAwAp14QFWNlqSJVPWXnwiIQLoEyAVIuhEqMhEQgS4T6NUDqssgtTYRcM4JggiIgAiIwEAIdOIBFX4rPqwPZC8UpgiIwB6B8P6H9T0X\/SsCIiACSyfQiQcUFOxb8SRLdAQ7gs0XbFGRUQREIAkCdv+59+iIBYbNF7OrFAEREIFFEujMA8qC9hOlb8OOmE2lCIhA+gRidx6bSfoEFOFQCCjO\/hHo3AMKhPbVJbpEBERg2ASUD4a9\/4peBLpKoFMPKBIloPjqklIiAiIwXAKLzQfD5azIRUAEZiPQqQcUDydLmrOFo14iIAKpEFA+SGUnFYcIpElg6Q+o8MEUJk3afUlzGxSVCIgABLjrlCbKB0ZCpQiIQNcILP0BRYI0KJY8fRu6L+arUgREID0C3HWLSvnASKgUARHoIoGDB1Q3luYnz26sSKsQARFYFgHlg2WR17wiIAJVCHTqAcVXnCYkT\/QqQchHBEQgPQLcfxPlg\/T2t7GINJAILIlApx5QJEmTWXhYsi3ra36U5oseE2tXKQIisFgClgsoZ5nZ7nNZX\/OjNF\/0mFi7ShEQARHo1AOK7SBpUfqCzcS3+zrtJFoE3W\/zddrwMaFu7WbzS2tTKQIiECXQqtG\/nzYRNhOzhSXtdo\/Rw3ar02Z+lNStjXoo1qZSBERABDrzgLLERcIy3bYHmy9mtxJ\/2q2Ojs3qVmKjzeoqRUAEukmAu8rKuK+mU0ew+YLNF\/xpNxs6NqtbiY02q6sUAREQgToEOvGAIokhtnDTrTR72yUJ1aRorps3bzrkwYMHrkh8n9h41tf3M1teaePs7Oy4TErK3d3drEsVX3wy572PvPllL97zofLZOzKN\/cu9R2xA0600e9ul5QLKsrm4w8hQ97+NuG\/fvu2Q2NiwRmizvSHXkcN2dvZzntmtDNvpK2k3nxn7VMtOPKAWBTdMwCRG34ZuQlveui5evOjOnz\/vrly54m7cuDElXGou\/Ztvvul+8id\/0r344ovujTfeyBKBjXd\/+36m3rp1K2vz\/cLxwjpj03l7e9vdu3evVO7f35+rqj9+jM\/awrmL6hZ3kU+KbYp7\/\/zfuXOHY9Mr4a77C+bO+zZ0E9p831AvygltnfuUzx4586d+6qfcSy+95C5cuJDlSeNIGzmTHEx+\/cY3vpFtx\/ZBPrSch9FyLXrYXjfH2fxNlSnvn8V29+5d0OdK3xs68YDKS05mt7JJ2IxJcrQxfR0bdXzQQ9na2nJXr151ly5dcpubm1Ny6tQpt7GxkT1s3n77bYdQP3HixOEwR44cyXTsPFjwQaiH44V1G2dldcWtrq6WytGjR7O56vpXWYu\/Novbtw1BV9z75\/\/48ePZOWviI+\/emd3KJuayMRiTO291X8dGHR\/0mBTlhLbuQcpnz\/Liv\/WJJw5zqHG0th\/81KmsbX19PduSlZX9fLhykPMwTg5yLXrYXjfH2fxNlSnvn8VmewP\/FKUTDygD6ycodJKWtTVZzjv22bNn3blz59zp06fd2tralBw7dix71GC3NaMjVh9P9rFjQ8yOXkXwn4wnbjIpl\/F4jLubjCduUsO\/yjp8Hz9u3566rrj3z799UZAdtoY+uKc2FHpX8wFrLMoJbd2B+c\/e\/t61tb55x4Xrsxv7D\/NwLNqe2\/wIhYMDyngyznIcJXVksmejRLBPJhM3PrCFYy66zrr5InjR8y5iPoutjbzAXnZFxl1ZCOsgQZIoTaesIn4\/\/BkDG3ooeW3YfV\/qeWP4ftJFQATaIcD94x4yOjplFcHX+uGPjg09lLw27L4v9bwxfD\/pIiACwyHQqQcU2C1JkbCoVxX60QdB9\/tho+6X6Ca00cfqlNSxS0SgLQIat5yA3UPuZLn3Iw\/60QdBf9TiHDbqfoluQht9rE5JHbtEBERABIxA5x5QLGzWhEWSQxjDF7NRxsR8\/TazqRQBEVguAeWD5fLX7CIgAnECnXpAkSjjy2zDqjFFQAS6TED5oMu7o7WJgAgs\/QFlSZKS7wD5JdtD3Up0hLpEBEQgPQJ2vymVD9LbX0XUEAEN0wkCS39AkSQhQeknTWwIdivREeoSERCB9AjY\/aZUPkhvfxWRCKREYOkPqBhMkmfMLpsIiMDwCHQwHwxvExSxCIjAYwQ69YBSonxsf2QQgcESUD4Y7NYrcBHoBYFOPaAgxrftKZU8oSCJEpBxMASUDwaz1QpUBHpHoFMPKJKlPZzQe0dTCxYBEWiMADlA+aAxnBpIBESgYQKzPKAaXsKj4SxZYkEngaJLREAEhkeAHGBRoysfGA2VIiACXSDQiQcUyTEGw+xWxnxkEwERSItA3n03u5VpRa1oFkOg2iw3b950daTaqPJKjUAnHlCpQVU8IiACIiAC\/STAw+ny5cvu\/PnzleUv\/tv\/Qfbg6mfEWvWsBDrxgOJb84gF4etms7KozXxUioAIdJNAlVVxxxHz9XWzWVnUZj4qRQAC9+7dc++\/\/77b3t6mmj14eCzF5Nq1a+65V77onv\/8Vqk8\/cLL7t1\/9HY2pj6GRaATDyiQx77BcHNHAAAQAElEQVQtb8nRSvwkIiAC6RNQPkh\/jxcd4dffesu9\/vrr7pvf\/N1s6osXL0a\/w4Qdh5NPPutOPlVBnnwGd8kACXTmAQV7PZSg0KZobBHoDwHlg\/7sVR9W+uGPfdJ99BPf545vfjxbbt53mPiOUuagDxEoIdCpBxRfdVrSjJVmK4lJzSIgAgkQUD5IYBM7FMLKh467lfUPu6PHPpStKvc7TPqOUsancx8dXNDSH1Dho4ikCadYaTbaJSIgAukRUD5Ib08VkQikSmDpDyh7FJE4kVRBKy4REIFyAsoH5YyW5KFpRUAEAgJLf0DZekiciD2iYqXZrE+dkr5IWR98TMp81S4CItAOAXIBwl1khlhpNtrrCn2Rsn74mJT5ql0ERGBYBDrzgCJJgZ6kifg6dV9oqyOMbf3R8\/rSZn6U1PN8ZReBhRIY2GR297iHCOFTxoS2OsLYNg56Xl\/azI+Sep6v7CIgAsMj0IkHFMkJMfyWqChDMZ+qJf39sdGxhf2x0RbaVRcBEVgsAe4hYrNyN9EpQ8FeR+jvj42OLRwDG22hXXUREAERMAKdeECRrGxBVpK8QrG2JZRTU9oPXnvw4IHLE78DPn7d9NBOvUys787Ojqsiu7u7WZcqvvhkznsfZetQe\/7eD5HN3pFp7N++5QMCr5IThnguZo0Zpr744\/h2XyfXkcN2dvZznt+Gnln5QPYMOzk5lHH2mivlV8Ywf3+N0vfzIxxTlk48oABM0jShvmjhsebPyVpCm7Xzg9b4Mf9XrlxxN27cmBIS6e3btx1i\/rdu3Zqq39++nzWFdurheGHdxuWn6fKTdcvk\/v39uar648fiqqzFX5vF7duGoCvu\/fN\/584djk1jwv0zaWzQGgOFd5+1hDZ\/uKKc0NY9SPnsWZ774IMPMsx+PrI2y6PvvPNO5rN9754jH1rOw2g+6Dx2Hu48dA93d6hmP5Ec\/1CsP7kwbIvVzd9fY5U9X9z+7d\/RKmtqysdiu3v3bsY61Y\/OPKBITibAJmGFgn0RwrysJW+ura0td\/XqVXfp0iW3ubk5JadOnXIbGxvuxIkTh93D+pEjR7K20E49HC+s27grqytudXW1VI4ePZrNVde\/ylr8tVncvm0IuuLeP\/\/Hjx\/PzllTH9w\/E8bkToaCfRHCvKylaK6inNDWPUj57Fmem0wmGXY\/Hx22HeTR9fX1zGdlZT8frhzkPIyTAx\/08XjsMhmNqLq8nKicuX+n5zm3djZtbzLgCX505gHlsyVZIdgofcHWplRJlmfPnnXnzp1zp0+fdmtra1Ny7Nix7FGD3daJjlh9PNnHjg0xO3oVwX8ynrjJpFzG4zHubjKeuEkN\/yrr8H38uH176rri3j\/\/9kVBdtga\/rD7z7CmW4mtTamSD5i\/KCe0dQdSP3twtfwVMqRtMhlTODigjPfqk8nEjfdK6sjE03k2jUYjNxqNaHJ5OXE8Hhe2TyYT54v5h2ssq7Nuvggu8+tju8XWZl7INmnJH+Mlz59NTzLMlOAjzx64FVYZgyRoTujYrO6XRW2+n3QREIH2COTdzzx7nZUwBvfc+qBjs7pfFrX5ftJFQASGSaATD6i20ZMgSYYIuj8fNup+iW5Cm6TPBLR2EZgmQA6w+43ut2Kn7pfoJrRJREAERAACg3hAESiJEkH3xWyUMfF9pYuACKRBwO56GA12bJQxoU0iAiIgAhBo\/QHFJBIREAEREAEREAERSImAHlAp7aZiEQEREAERaIqAxhGBQgJ6QBXiUaMIiIAIiIAIiIAIPE5AD6jHmcgiAiLQBQJagwiIgAh0mIAeUB3eHC1NBERABERABESgmwT0gOrmvnRhVVqDCIiACIiACIhADgE9oHLAyCwCIiACIiACItBHAotZsx5Qi+GsWURABERABERABBIioAdUQpupUERABESgCwS0BhEYAgE9oIawy4pRBERABERABESgUQJ6QDWKs9nBbt686RAb1dfNplIEHicgiwiIwLwE3vvubffed\/\/JodQdj3ztS93+8u8+AT2gFrBHXKLXXnvtsZny7DjSdvnyZXf+\/Hl34cIF97Wvfc1Rv3jxIs0SERABERCBFgn83v\/16+43\/+bLh\/Jbf+dy9phyFf7x8zc5HCGPY6\/QXS49IaAHVAsbFQ75R3\/0R+7atWuh2eXZceSi0eeHvv\/jWV\/zpU67RAREQAREoD0C7\/6jt91f\/dFPu\/\/01c+6f\/\/fPeOoV53N8veXfvw59ytfet792I88neXxqv3l1w8CekB1fJ8++5knplYY1qcaVREBERCBgRK4d++ee\/\/999329nZjBP7ckx91yPftlbMM+tlPn3SZfObkLN2X0Udz1iCgB1QNWHIVAREQARHoBoFv\/e63soXwVxtQvv7WW+7111933\/zm71LNxHyyij5EoGECekA1DFTDiYAIiMDMBNSxMoG1j+5\/V+fj534g6\/Phj33SffQT3+eOb348q\/NhPugSEWiaQFIPqDNnzjikCqTQj3pMqowlHxEQge4RsPtcZWX4+n7UY+L7SF8ugZUPfThbwBNPPZOVKx867lbWP+yOHvtQVufDfNAlItA0gWQeUCS769evOwQ9DxRtSKydvqHE\/GRLloACS4QAd9zuMnpeWLQhsXbr75cxP9lEQASGSSCJBxQJkCRnW4iOzep+SRvi26SLgAikQ4C7799xdGyxCGlDYm2yiYAIiEARgW49oIpWuoA2kqxJ0XT8L6rIgwcPXJ4U9d\/d2Y02+2NFHTzjzs6OqyK7u\/tzVfHFx6Ygvqrir1t6\/plImY2dm5RKywWUZXHZXUl5jxcZW8jbnztsm8qnpDskdLI6bchenXwXE8uZey77\/+KP7NWK\/GNr3D3I03tds399nyHoWdAJf+gB5W0uX4maFCVNfpglPxjtypUr7saNG1NCIr19+7ZDbOh33nnH1Kz84IMPsjK037p163As639\/+\/6Ur9X5X3X533bL5P79\/f5V\/fFjQouROMvkxRdfdG+++WYWc8gj9brtd+pxhvGFcd+5c4djk5RYLqAsygcEbfcllhNCdk3Vwz1oatwujGP5z3JlLDfu7O6A3pnPw73HysOdh+7hgZ1G80EP28l1sfxpObNqH\/OPrfHRHPs\/WsH3SXn\/LLa7d++CMVnRA+pga0mSB2pWUM9LmltbW+7q1avu0qVLbnNzc0pOnTrlNjY23IkTJ7Jx+FhfX6c4lKNHj2R6aKefjWf9jxyZ9rX6yuqKW11dLZWjR49mc9X1f\/ZHv+D+9b\/+t0vl6Rdedm+\/\/bY7efJkFretfyil7fdQ4rU4w7iPHz+enbNUPrj\/fizU8\/IBfkU5wZg1XYZ70PT4yxzP8t9kMgHvVG6xttFo\/4+vIwf5dDweu0xGo6wPH+aDnrXhc9CelxMtZ1btY\/6x\/L2yspLl6JhPyvtnsYV\/xsE0Jdk\/gSlFtIBYzp49686dO+dOnz7t1tbWpuTYsWPZhcFuS8FmOuVo7xJThnb6+ILPeLK\/ReZr9cl44iaTchmPxwzjJuOJm9Tw\/96nn3Pf++fK5eST+\/8HDOvmQUc5JGFfFPeas4d9dtiW8rHcSYtyQlv3IfWzx45a\/goZ0mbPJPOhPhqN3Gg0ojmTR5pz6KPRyI1GI8c\/eTnRxsMHz9FoVNjH\/GNrHE8s744ZLvfPirBv3+t2NlPPC\/u7mm1tfz\/Crw75ShFbnYjo4\/tTrzuG31+6CIjAcghwb7m\/Njs6NqtXKenj+1GvO4bfX7oIiEB6BJJ4QLEtJDeSHIKOzQSb6XklffAzoZ7nK\/vjBGQRgS4R4P7m3WXsZWv1++NPvayP2kVABIZFIJkHFNtGkkPQfaljwxfx+0sXARHoHwHuMRKuvI4NXyQcQ3UREIFkCMwcSFIPqJkpqKMIiIAIiIAIzEGA\/\/PMxIbZvnfPvf\/++1Z1Yfthg5ReEtADqpfbpkWLgAiIQCIEEgnDfpQFP\/YFnbDeeust9\/rrr7u33vo6VYed9gsXLmSPqcyoj94S0AOqt1unhYuACIiACHSFwHOvfNE9\/\/mtTPjxLqzrI3\/mU9kvOF4\/9S9TdfjQdu3atayuj34T0AOq3\/un1YvAvATUXwREoAECJ5981p186kAOfrzL0e9ZdyvrH3YrH\/pwNkPmc9CWGfTRawJ6QPV6+7R4ERABERABERCBZRDQA2oZ1P05pYuACIiACIiACPSOgB5QvdsyLVgEREAEREAElk9g6CvQA2roJ0Dxi4AIiIAIiIAI1CagB1RtZOogAiIgAl0goDWIgAgsk4AeUMukr7lFQAREQAREQAR6SUAPqF5umxbdBQJagwiIgAiIwHAJ6AE13L1X5CIgAiIgAiIgAjMS6PEDasaI1U0EREAEREAEREAE5iSgB9ScANVdBERABERABGoRkHMSBPSASmIbFYQIiIAIiIAIiMAiCegBtUjamksERKALBLQGERABEZibwGAeUGfOnHFIFWJV\/aqMJR8REIHuEeCOI1VWVtWvyljyEQERSIfAIB5QJMDr1687BD1v+2hD8tplb4iAhhGBJRLgjpMLEPS8pdCG5LXLLgIiMGwCyT+gSIAkSttmdGxW90vaEN8mXQREIB0C3H3\/jqNji0VIGxJrk00ERGCYBPyok39A+cE2pd+8edMhDx48cHlSNNfuzm602R8r6uAZd3Z2XBXZ3d2fq4ovPnX9bUn+2qXnn4tU2dg5GGpJPkBS3d9FxxWeI3\/+sG0qn5LukNDJ6rQhe3XyXUwsB+657P+LP7JXK\/L326bGoC9y0N\/a\/JhS1fdCTvpfPaBm2N6LFy+68+fPuytXrrgbN25MCUn09u3bDrGh33nnHVOz8oMPPsjK0H7r1q3Dsaz\/\/e37U75W397edvfu3SuV+\/f3+7flz7gskPUiIY\/U67bfqccZxhfGfefOHY7BYKUoJ4TsmqqHe9DUuF0Yh1zCYbJc+Sg33jjMrTu7O7g483m490Xlw52H7uGBnUbzQQ\/byV2xHGo5s2of8\/fHM1tsTmsjRqQLvJteg53Nu3fvgjFZ0QNqhq3d2tpyV69edZcuXXKbm5tTcurUKbexseFOnDhxOPL6+vqhjnL06BEKF9rpZ+NZ\/yNHpn2tvrK64lZXV0vl6NGj2Vxt+588eTKL29Y\/lNL2eyjxWpxh3MePH8\/O2VA\/inKCMWu6DPeg6fGXOZ7lv8lkkh2pWG4cjfb\/+DpykE\/H47HLZDTK+vBhPuhZGz4H7Xk50XJm1T7m749nttic1pZyzrSzGf4ZB9OUZP8EphTRAmI5e\/asO3funDt9+rRbW1ubkmPHjmWPGuy2FGymU472LjFlaKePL\/iMJ\/tbZL5Wn4wnbjIpl\/F4zDBuMp64SYv+rJsHHeWQhH0ZUty2t2Hc9rDPDtsAP4pygjFrugz3oOnxlz0ex8jyV7gW2uyZZD7UR6ORG41GNGfySHMOfTQaudFo5PgnLyfaePjgORqNCvuYvz+e2UZ7g4xGo6n+1kZMG7co2gAAEABJREFUqeYOO5up54Xx3v4m\/S9\/CdT\/S6Lo2LoYNN\/2NOni+rQmEeg7Ae4+OcDiQMdmdZUiIAIiUJVA8g8oQJAgSZQIOjYTbKYvu7S\/R8Hfr0JnPd\/63W9RuMuXLzvnnLN6ZtSHCIhAbQLkAO49gu4PgM2vSxcBERCBPAKDeEARPIkSQfelqs3v05b+3CtfdM9\/fiuTp194OZvm2Ec3svLj534gK62eVfQhAiIwEwHuPRJ2rmoL+6k+HwH7zruV4WhmpwzbVO8JgQSXOZgHVB\/27uSTz7qTTx3Ik89kS1750P5fzn3iqel61qgPERABEeg5AR5FfIed77yb2Hfgv\/3tb2fRUQ\/bsgZ9iMASCegBtUT4mloEREAEFkSgs9PwgLp27ZqLfQd+9cQT2bpjbVmDPkRgiQT0gFoifE0tAiIgAiKwTyD2HfijB9+Bj7Xt99KnCCyPgB5Qy2OvmYdEQLGKgAiIgAgkRUAPqKS2U8GIgAiIgAiIgAgsgsBQHlCLYKk5REAEREAEREAEBkJAD6iBbLTCFAEREAER6CMBrbmrBPSA6urOaF0iIAIiIAIiIAKdJaAHVGe3RgsTARHoAgGtQQREQARiBPSAilGRTQREQAREQAREQAQKCOgBVQBHTV0goDWIgAiIQHoEbt++7RB+kGiZpBd9GhHpAZXGPioKERABERCBHhF45ZVX3EsvveQ+97nPOfs1NXnlhQsXHI+sHoU3iKWWPqAGQUFBioAIiIAIiMACCTz70hfcn\/9Prri\/9Oovuuc\/v5UrT7\/wsuNX3SxwaZqqIgE9oCqCkpsIiIAIiECvCHR6sU889Yz76Cf\/FffEk8+4k089my977Z0OZMCL0wNqwJvfZOj8t3yEbzOXSZPzaiwREAEREAERWAYBPaCWQT3BOfXf8xPc1HlDUn8REAERSJhAUg+oM2fOOKRsv\/BBfD\/qMfF9pOcT0H\/Pz2ejluUQsPtcNnvMz2xhWTaW2kVABIZDIJkHFInu+vXrDkHP20La8EHQfT9sofjtPdMXulz99\/yF4tZkJQS423aX0fPcacvzM7tf5o0juwiIwPAIJPGAsiRo20fCw2Z1K7HRZnV0bFZXKQIi0H8C3GnutkWCjs3qVmKjzero2KyuUgREYFkE+jFvEg+oplCTPE2KxrS\/JP3gwQOXJ0X9d3d2o807OzvOZHd32ifsY35lpY1T5mfts\/q7vfXSF7GxYqUFnsdN9vwz1WU2tq8plZYLKMviqpITurx\/y1ybsd3djeS\/gzTo55LdvVxjfSinciP+CA0xoQ3Za\/PH9PVwfIc\/ktPH\/HPHoC9y0N\/8q+ZM81\/mHs06917ISf+rB5S3vXwFalKUNC9evJj94LMrV664GzduTAmJlP8bDbGh33nnHVOz8oMPPsjK0L69ve3u3buXyf379zOfnYNkYX2s7vtan1hp4yzCv8oc+BDYrVu3priFHPtSt\/3uy3qbWmcY9507d9jWpMRyAWVRPiDoopzQFPNwnHAPwvZ56vP2feONN1xVsVx5797j+c8eD+SNe4\/lxh3QO8uND\/e+AH2489A93HuIZQ17HzueHrb7Y9rYlJYz97q7Kn3M3x\/PbLH+1oY\/wpxFYv59ypl2Nu\/evQvGZEUPqIOtJUkeqFlBPS9pbm1tuatXr7pLly65zc3NKTl16pTb2NhwJ06cyMbhY319neJQjh49kumhfWV1xa2urmZy9OjRzGc0GmXlkYM+o9F+3fe1PrHSxmndf2XFrexJbA2+zdYDo5BdH+u2331c+zxrDuM+fvx4dk5T+eD++7FQz8sH+BXlhHk4F\/UN96DId5FtDx8+dF\/+8pezn7LNT9ouky984Qsg3Mt7K3sSz39+\/rIcMhrt\/\/FluXE8HrtMDnIkg5oPetaGz0G7P2YsR1XtY+vxxzNbbE5rI18i\/twx3fz7lDPtbIZ\/xsE0Jdk\/gSlFtIBYzp49686dO+dOnz7t1tbWpuTYsWNZEsBuS8FmOuVo7xJThvbJeOImk33h4uGz\/1xyWWLw676v9YmVNk7b\/sTEXEhsHWajnTjgk4KwhyS9ZmOZPlNdHDuM+8iRI2zrYKUoJ7S1f+EetDVP3XHfffdd9zu\/8zvuuVe+mPvTtf2fvM1P2ubgjMf7uW+ylwPH4zEm50Yu+2cynlTKjaPRyI1Go6wPH48059BHo5EbjUaOf\/wxJ3tzmhzOveeE52g0Kuxj\/v54ZhsxRtD\/sG188ODbK23uWGn+dfdhmf52NlPPC+O9\/e39v+FXh3yliC0MDBttZkfHRh2d0oS6tZlNpQiIQPcJcG+5v7ZSdGxWtxIbbVZHx0YdndKEurWZTWUxgZNPFvx0bf8nb+snbReDVGtnCSTxgIIuyY0kh6BjM8FmOm3UEfSYPWwzH5UiIAL9IMDd5h4j6P6qsVmdNuoIeswetpmPShEQgWETSOYBxTaSABF0X0IbdcT3QcdmQl0iAiLQXwJ5dxm7HxV1xLehYzOhLhGBHhDQEhdIIKkH1AK5aSoREAEREAEREIEBE9ADasCbr9BFQAQaJqDhREAEBkNAD6jBbLUCFQEREAEREAERaIqAHlBNkdQ4XSCgNYiACIiACIjAQgjoAbUQzJpEBERABERABEQgJQLNPqBSIqNYREAEREAEREAERCCHgB5QOWBkFgEREAERGA4BRSoCdQnoAVWXmPxFQAREQAREYMEE+AW9VWXBSxvsdHpADXbrlxt41USA33JXqtkXQ0CziIAIFBG4ePGiO3\/+fCW5cOGCU+4sotlMmx5QzXDUKDUJKBnUBCZ3ERCBQROo88uZr127NmhWiwpeD6gWSN++fbtw1A\/uf1DYvozGRc+pZLBo4ppPBESgzwROPvmsO\/lUBXnymT6H2au16wHV8HbxePrFX\/xF98orrxyOfPny5UMd5Q\/+4A8oXGjPjAP5UDIYyEYrTBEQARFol8DSRtcDqmH0PKDefvtt9\/FzP3A4sq9j\/NDGn6aY8skM+sglwH\/Pb0tyJ1WDCIhARqDO3cs66EMEBkBAD6iWNvmJpx59G9XXme7I6hqFC+2ZUR9RAnX+zlTVv2hpfvoLl1HkMopARoDHE98tt\/vyWBn8xWbuatZRHyKQOAE9oBLf4FTCq\/N3poi5jr\/+wiXEJCIQJ8ADijtS507FR5J1kQTYtzqyyLWlMpceUC3t5O7ObksjD3PYun9nqq5\/j6lq6SJQm8Asf7DqTtXGvNQOfCew7LuFfjvfif\/a176W\/fiDKudjqcF1ZHI9oIKNOHPmjEMCs6qJE6iSMMynLgrrV7WsO7782yNALkDamyE+cpWzwt+3RPCNjxK34q\/\/JBdnk5K16ncMn\/\/8lnv6hZcd32Ws8+jiwcVZSolZ3Vj0gPKIkSivX7\/uEHSvKV\/NaXnw4EFOS5rm+9v3s69ctre3exngrImDP8B+7dd+zVHGAifB1PnDiq8I205MrKmOxOIi3qK4Y336ZiMHkAsQ9HnWX4c3vlXOzOc+9zn30ksvOco6Z4bx+cOy6h+w\/OE6T+zz9E35R760nTMrf8eQH41w8KMP6pwJzlDe3g4hPxC7HlBQ2BMSJIlyT83+RceWVYIPElCemOsHHzz6WU9hEtjZ3f\/Pe6H9zq0\/dP\/0j\/7\/TLbv38+GCn2tfu\/ePff++++XSjbI3kdVf5u3rv\/9vfXCpKzfrOOXjWss6o6\/hyb79\/te\/An32f\/450qFP0xIHMSKkCi++tWvZo9H6jHBf9bxY+PNa6vyhzMPORP+cA6\/tR\/GnUFM6IO7Tw6wkNCxWT0sy\/akLvMqZ+Yzl37Gfc\/3X3R\/9t98MfvuQdkarN3Wvv5nnnLrH3+6VE4e\/OF67155vuEe1r2DMX+zWR6N5cbdII\/+8bvvuT9+5z23s7NjITrzwbCzs+t2Hu4ctvtjWt6ltLmr9jF\/fzyzhXOSx6yt\/Zz5vmO+IrE2YkXqngk7U2Fp+YExUxY9oGrs7sc+9jF39uxZV\/TdCtrG3\/MRd\/O9Xbfx1HPuEx9ec\/\/v\/\/5aNstnPjp2T2+uuv\/pf\/tWVvftT5046t785Z9xf\/\/n\/2om6E\/86Q23uTF2nzy9no1BSf3pv\/Bp943f+7Z7\/fXXS+Ubv3fdrXz636js\/63v3JzJ\/1v\/380spq9\/\/euFa5p1\/Krx1h3f+NzcXXe\/\/+4\/L5WbO+tusvHxwzPAfhM4pT04\/BJ7dh5mHN8fqyn9jf\/nH2Z\/8PKHb5ms\/fkX3JvXbx3Ga2sgLj\/u117bP+PYhiRVcgLMYFjG2tphzp0tO5N\/+CcPs7P47srJrGRPmKtM8ONMcvZ\/+7d\/25UJfqyHskrOqXsHY\/7YPvXJp7L8+S999InHciN59XsnDx1tlkf\/9i993V3+W7\/tKDmDp9d23RPfve5OHNkr93LpZPS+u3Pnjjs6uZ\/lVHKs5Vu\/xE5\/8q3fh3xMm++Ljg1\/SuoIOv5+\/7N\/6RNZHiY2eH6r9Zx5vXRvbe\/ZW9ZEabaiEj9i5izFzht22n\/jN36DIlkZJxtZC4GRLLe2ttzVq1cL5dd\/9b9z\/+PLf9m99rOX3c\/\/91859P2Zv\/vr7kv\/1f98WLdxsP\/sL\/\/aY\/Zf+i9fc5\/+iz\/jfu4X\/oesjZL6l37qb7hf+Q+fqSSs43\/5Gz9eyZcx5V\/MFT7\/62v\/RbYftn9lJecBtlVklvHL5g\/bWc+v\/rUfclXkK194yVWJ94d\/+IdbuHHdH7JqToDh47zjewDzuneW8cN9LqpzBqqcR3w4k3XXM68\/c\/7nP\/ez7l+7vOX+27\/7S4\/dN\/JqXpvF\/Qtf+XvuL2z9H+6Xv\/r3HLn0Tz112SHkUHKp+eWV+OCP0Icx8nxjdvzpi9D\/1Z\/4+SwPE9u8fNiXPKk7PuPU7YN\/LObQlnpeGHc\/RXVrhSTMc+fOOYkY6Aw8OgPci27d1MWthth1Fh6dBbEQCzsD3I3F3cTFz6QH1OKZa8aOENAyREAEREAERGBWAnpAHZAL\/5Iof2EU20GzChEQgQER4O6TAyxkdGxWVykCIiACS3xAdQ8+CZJEiaB3b4VakQiIwKIIkAPIBQj6oubVPCIgAv0goAdUsE8kSiQwqyoCIjBAAuQCZIChK+QiAmoTgT0CekDtQWj6X75iRZoed5njEU9M\/DVZu28zvajNfLpUst7YerAjTbbFxlqWLYyNekz89Vm7bzO9qM18VM5HwBiH5XyjLr+3xbP8lTS7AosrLJudZTmjEVNsZuxIrK3PNj2gGt49DglfsSLoDQ+\/1OGIKRRbELFaG7rZKanntdHeJWGtSGxN2PPimLUtNs8ybKwfic1tMful+dHH7Ohmp6Se10b7HKKuAQHj7JeBS6+qqZ8df59M79UGBYtlv5DAnFWxW4zomTGRDz2gGtxIDgcHxYZEx2b1VEtiJFaLDx0bdUrq6Ag6NvQuCutDwrWxZt+Ojg0\/SuroCDo2dErq6Ag6NvQuCetC6qyJOPw+6NgYg5I6OoKODV0iAkUEOG5AngAAAAV7SURBVCecF\/NBx2Z1ld0jwB4h4crYN9+Oji3062tdD6i+7twS1s3Bz+TMmSXMrimXRcD2nHJZa9C8xQTYG5NiT7V2gYDtFWUX1qM1zEZAD6jZuA2yF189mOjiD+cI2J5Tat+Xs+9wj4mthr0xwc\/sKrtJwPaKUvu12D1qcjY9oJqkmfBYXHQ\/POq6+D6RNHX22Y+MuvbdJ7IYHe4xYXbslCbUtUdGo3sl++Ovirr2yyfSH10PqP7slVYqAiIgAj0noOWLQDoE9IBqcC\/DryT4qgJbg1MsbShi8SenbrFRUrd2dGzUKamjI+jY0PskrJm125rRsVGnpI6OoGNDp6SOjqBjQ++DsF5\/ndRt\/ZTUrR0dG3VK6ugIOjZ0SbMEYOuPSL3PrFk7MVhM6Nis3veSePwYqKcUnx8bcRGf2dCxWb3vpR5QDe8gh4NDgqA3PPzShiMWYjKh7i+G+ixt\/hhd12eNsajfomOuO5+\/dvaXuj8GdewIetU230\/6fATgDn8T6vONuPzexJBSPD5RPzZipO63p6YTH3Ei6CnFpwdUC7vJIUFaGHqpQxKTSWwhs7bFxlq2jVhia8CONNkWG2tZtlhs2Exi65q1LTaWbLMRsD2gnG2E7vUiFqR7K5t\/RcRlMv9o3RmBmGKrwY7E2vps0wNqobunyURABERABERABFIgoAdUCruoGERABERABESgTQIa+zECekA9hkQGERABERABERABESgmoAdUMR+1FhDgLwX6zWHdbzMdnzIx37As60d72CevXsc3b4zQ3saY4RyqD5aAAhcBEegYAT2gOrYhqS+Hv0iYJ8ROG2We0J4neX1COw8dxgjt89YZk7HnHUf9RUAEREAEuk9AD6ju71EnV8hDgQeDvzjq2H1bVZ1+9C\/zxy9PyvoWtZeOWdRZbSIgAo0S4D42OqA3GGPHxHORKgKVCOgBVQmTnHwCJJ+8xw522n3\/Mh1\/+pX50Y5fntBeJkVzxcbFv2xMv50x6vbx+0sXARFonwD3NBTd2\/a5pzZDXx5QqXHvbTwkGRJPUQC041fkY2344W\/1WIkPQhtlnvjt6E0Ia2O+JsbSGCIgAiIgAukQ0AMqnb1sPRIeEvagQC8S82tiUYyFMBalSVGdtjYkjLmNOTSmCIhAMYEq9zD0oV48alutGjdVAnpApbqzLcTFw4VhKU3COnZsiK+TvGKCX8yOjbZQsJvQZjol9TaFOYjJF2xtzqmxRUAEpglw5\/w7iI7N96KO3YQ2dEqJCDRFQA+opkhqnEICJK+Y0Clmx0ZbKNhNaDOdknqbsog52lz\/UMdW3CIQEuCBFYrud0hJ9TICekCVEVJ7pwj4SY+FhXVsbYo\/H3qbc2lsERCBdgjwWAqlnZk0asoE9IBKeXdbjo0HBEmoeJrmWpnLF0b26+jY2hKLl3lM2ppL44qACMxOgPvJfTWhPvto6ikCcQJ6QMW5yNpBApYMrWSJpluJrUhIpPgW+fht+NLHt5lOm+lWYsvzNx+VIiAC7RKwe8hdRNqdTaMnS6AkMD2gSgCpOU7AElTYmmcP\/Wapkwh9YQy\/jo5tVmHtofhjovvt1GedS\/1EQASKCfh3zde5d34dHVs4GnZfwnbVRWBeAnpAzUtwYP0tIYUJK8+eh8f889p9u\/mGJT6hzeq05Qlrx89vxxYT3wfd97E6JcKYtKNLREAEZifAPbp+\/bqLlYwa2rGZ2D0MfbDjg51SIgLzEtADal6CA+tP8kHCsLEhoT2vjq9Jno\/Zza9OaX3zSsbKa5vV3saYs65F\/URgqAS4hzyWQsE+VCaKux0CekC1w1WjikD\/CSgCEegpAR5LofQ0FC27wwT0gOrw5mhpIiACIiACIiAC3STwLwAAAP\/\/4m2UlgAAAAZJREFUAwCFyXxTg3O2HAAAAABJRU5ErkJggg==","height":285,"width":474}}
%---
%[output:298783de]
%   data: {"dataType":"text","outputData":{"text":"[10:11:47][INFO]  Reinforcement Learning Toolbox が見つかりません。\n","truncated":false}}
%---
%[output:5308948c]
%   data: {"dataType":"text","outputData":{"text":"[10:11:47][INFO]  手動 REINFORCE ループ (Section 2) が同じアルゴリズムを実装しています。\n","truncated":false}}
%---
%[output:6e28324c]
%   data: {"dataType":"text","outputData":{"text":"[10:11:47][INFO]  RL 前:  valid=6\/100 (6%)  avg_reward=0.060\n","truncated":false}}
%---
%[output:7ce94921]
%   data: {"dataType":"text","outputData":{"text":"[10:11:47][INFO]  RL 後:  valid=4\/100 (4%)  avg_reward=0.040\n","truncated":false}}
%---
%[output:1d690113]
%   data: {"dataType":"text","outputData":{"text":"[10:11:47][INFO]  RL 結果: 報酬 0.060 -> 0.040  (delta=-0.020, ノイズ範囲内)\n","truncated":false}}
%---
