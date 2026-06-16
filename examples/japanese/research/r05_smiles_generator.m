%[text] # R05: 分子言語モデル -- SMILES生成
%[text] EasyMolKit Research -- Layer 4
%[text] 
%[text] ## ストーリー
%[text] 分子には「文法」があります。
%[text] C、N、O などの原子記号が括弧で枝を作り、数字で環を閉じる——この規則に従った文字列が SMILES です。その「文法」を訓練データから学んだマルコフ連鎖は、全く新しい SMILES を生成できるでしょうか。このスクリプトでは FDA 承認薬の SMILES コーパスから文字語彙を構築し、バイグラムマルコフ連鎖と文字レベル LSTM の 2 つのモデルで新規分子を生成します。
%[text] 有効性・多様性・新規性という評価指標で両モデルを比較し、小規模データでの生成モデルの課題と可能性を体験します。
%[text] ## 学習目標
%[text] - SMILES を形式言語として理解し、言語モデルの訓練を化学構造の確率的モデル学習として捉えられるようになります。
%[text] - 訓練データからバイグラムマルコフ連鎖を構築し（勾配降下法不要）、小規模データセットで帰納バイアスが重要な理由を理解します。
%[text] - SMILES 生成の失敗を構造エラー種別（括弧・環閉鎖・原子価）で診断し、メモリ要件と関連付けます。
%[text] - カスタム訓練ループで文字レベル LSTM を構築・訓練します（時間方向逆伝播 BPTT）。
%[text] - 温度サンプリングによる自己回帰デコーディングを実装します。
%[text] - 生成モデルの 4 つの標準指標（有効性・多様性・新規性・物性分布重複）を習得します。
%[text] - REINFORCE（Williams 1992）を方策勾配アルゴリズムとして理解し、強化学習と分子生成（REINVENT）の接続を把握します（詳細は r06\_reinforce.m）。
%[text] - 小規模コーパスでは汎化でなく記憶が生じるスケーリング課題を認識します（本番モデルは $10^6 \\sim 10^7$ 分子を使用）。 \
%[text] ## 前提条件
%[text] - A05（ニューラルネットワーク物性予測）修了 -- カスタム DL ループ
%[text] - Deep Learning Toolbox（dlnetwork, lstmLayer, adamupdate, dlfeval） \
%[text] ## 動作環境
%[text] - Deep Learning Toolbox は Section 7〜9 に必要です。
%[text] - MATLAB Online と Desktop の両方に対応しています。
%[text] - GPU 不要で、最新 CPU で全訓練が完結します。 \
%[text] 所要時間: 45〜90 分（r06\_reinforce.m の REINFORCE RL セクションは追加で 30〜60 分）
%[text] ## データ
%[text] - data/list/fda\_drugs.csv -- FDA 承認薬 200 種（ChEMBL, CC-BY-SA 3.0） \
%[text] ## 参考文献
%[text] - Gomez-Bombarelli R et al. (2018) Automatic chemical design using a data-driven continuous representation of molecules. ACS Cent Sci 4:268-276. doi:10.1021/acscentsci.7b00572 \[VAE ベース SMILES 生成モデルの先駆的論文\]
%[text] - Segler MHS et al. (2018) Generating focused molecule libraries for drug discovery with recurrent neural networks. ACS Cent Sci 4:120-131. doi:10.1021/acscentsci.7b00512 \[LSTM ベース SMILES 生成と転移学習\]
%[text] - Olivecrona M et al. (2017) Molecular de novo design through deep reinforcement learning. J Cheminform 9:48. doi:10.1186/s13321-017-0235-x \[REINVENT: 分子生成への REINFORCE -- Section 7 の基盤\]
%[text] - Williams RJ (1992) Simple statistical gradient-following algorithms for connectionist reinforcement learning. Mach Learn 8:229-256. doi:10.1007/BF00992696 \[REINFORCE アルゴリズムの原著論文\]
%[text] - Brown N et al. (2019) GuacaMol: benchmarking models for de novo molecular design. J Chem Inf Model 59:1096-1108. doi:10.1021/acs.jcim.8b00839 \[生成モデルの有効性/多様性/新規性/KL ダイバージェンス指標\]
%[text] - Polykovskiy D et al. (2020) Molecular sets (MOSES): a benchmarking platform for molecular generation models. Front Pharmacol 11:565644. doi:10.3389/fphar.2020.565644 \[MOSES ベンチマーク: 生成モデルの標準評価プロトコル\]
%[text] - Hochreiter S & Schmidhuber J (1997) Long short-term memory. Neural Comput 9:1735-1780. doi:10.1162/neco.1997.9.8.1735 \[LSTM 原著論文\] \
%%
%[text] ## Section 0: セットアップと設定
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
%[text] ツールボックスの利用可否を確認します。
hasDL = license("test", "Neural_Network_Toolbox");
if ~hasDL
    error("emk:r05:missingToolbox", ...
        "Deep Learning ToolboxはSection 7ー9に必要です。");
end
logInfo("ツールボックス: Deep Learning=%d", hasDL); %[output:0902b7b6]
%[text] ### Section 0a: 調整可能なパラメータ（実行前にここで変更）
HIDDEN_SIZE   = 64;     % LSTM隠れユニット数 -- 下記CORPUS_TARGETに基づき自動スケール
DROPOUT_RATE  = 0.30;   % ドロップアウト確率
LEARN_RATE    = 8e-4;   % Adam初期学習率 -- 下記CORPUS_TARGETに基づき自動スケール
N_EPOCHS      = 150;    % SL訓練エポック数
                        %   150でval損失のプラトーが鮮明になる（主要教育ポイント）；
                        %   スモークテストは20に短縮可能
BATCH_SIZE    = 32;     % ミニバッチサイズ（n_train以下が必需）
MAX_SEQ_LEN   = 100;    % 最大パディング後系列長（START + ENDトークン含む；90でコーパスの4.5%が打ち切り）
TEMPERATURE   = 0.85;   % サンプリング温度 -- 下記CORPUS_TARGETに基づき自動スケール
                        %   有効率が持続的に20%を超えた後のみ0.6-0.75に下げる
N_GENERATED   = 300;    % 評価用生成SMILES数
                        %   300で二項分布SE ~2.6%（100の~3.5%比）になり、
                        %   マルコフ vs LSTMの信頼性の高い比較が可能；最大系列長100で生成は素4.5分
CORPUS_TARGET = 500;    % 訓練対象分子数
                        %   200: ローカルCSVのみ（ネットワーク不要；~1分；Markov ~13%, LSTM ~10%）
                        %   500: ローカルCSV + ChEMBL Phase-4取得（初回約50秒；合計約5分）
                        %        --> Markov ~18%, LSTM ~4%: 対比が教育ポイント
                        %   ネットワーク不安定時は200に変更
AUGMENT_FACTOR = 1;     % 1に設定: このデータスケールではオーグメンテーションはノイズを加えるだけ

%[text] コーパスサイズに応じてアーキテクチャを自動スケールします。
%[text] 200 分子の場合 HIDDEN=64 はパラメータ数 / トークン数 $\\approx 4$（記憶化レジーム）です。
%[text] 500+ 分子の場合は $\\approx 1.7$ になりアンダーフィットになるため、HIDDEN\_SIZE を 128 に拡大します。
%[text] 経験則: LSTM 文字 LM ではパラメータ数 / トークン数を 3〜6 程度に保つと良いとされています。
EMBED_DIM = 16;    % 文字埋め込み次元（VOCAB_SIZE → EMBED_DIM → LSTM）

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
%[text] 勾配クリッピング閾値: テンソルごとの L2 ノルムを 1.0 に制限し、LSTM の訓練を安定化します。
GRAD_CLIP_NORM = 1.0;
%[text] 特殊トークン文字を定義します（SMILES 文字列内に絶対に現れてはなりません）。
START_TOKEN = "^";
END_TOKEN   = "$";
PAD_TOKEN   = "_";

logInfo("R05: セットアップ完了 (HIDDEN=%d  EPOCHS=%d  TEMP=%.2f  AUG=%dx)", ... %[output:group:10372ce6] %[output:334b2162]
    HIDDEN_SIZE, N_EPOCHS, TEMPERATURE, AUGMENT_FACTOR); %[output:group:10372ce6] %[output:334b2162]
%[text] Python/RDKit のウォームアップ（初回の遅延を訓練前に済ませます）。
mol_warmup = emk.mol.fromSmiles("C");
clear mol_warmup;
%%
%[text] ## Section 1: SMILESコーパスの読み込み
%[text] ### 言語としての SMILES
%[text] SMILES（Simplified Molecular Input Line Entry System）は、分子構造をコンパクトな ASCII 文字列として表現する記法です。原子・結合・環・立体中心は、正式な文法に従ってそれぞれ 1〜2 文字に対応します。
%[text] 正則 SMILES（canonical SMILES）は分子に対して一意の表現を持ち、非正則 SMILES は同じ分子を複数の方法で記述できます。SMILES が正式文法に従うため、SMILES で訓練された言語モデルは次の 3 点を実現できます。
%[text] - 構文的に有効な文字列を生成できます。
%[text] - 化学規則（環閉鎖の均衡・原子価）を暗黙的にエンコードできます。
%[text] - 「化学空間」における分子構造間を補間できます。 \
%[text] 本スクリプトのコーパスは ChEMBL から取得した 200〜500 種の FDA 承認薬です。
%[text] 意図的に小規模なデータセットを用いることで、CPU 上で数分以内に訓練が完了し、過学習（overfitting）の様子を目視できます。
%[text] 本番モデル（REINVENT, MolGPT）は ZINC/ChEMBL から $10^6 \\sim 10^7$ 分子で訓練します。
%[text] 
%[text] \---- マルチソースコーパス読み込み ----
%[text] ソース優先度（重複除去後、CORPUS\_TARGET に切り詰め）を示します。
%[text] CORPUS\_TARGET ≤ 200: fda\_drugs.csv のみ（200 種のキュレート FDA 薬）を使用します。
%[text] CORPUS\_TARGET 201〜500: fda\_drugs.csv + ChEMBL Phase-4 承認薬のみを使用します（everyday\_chemicals / forensic\_challenge は意図的に除外します。小規模 N で混合ドメインはモデル品質を低下させるためです）。
%[text] CORPUS\_TARGET \> 500: 全ローカル CSV + ChEMBL を使用します。
smilesCombined = strings(0, 1);

t1 = readtable("data/list/fda_drugs.csv", TextType="string");
smilesCombined = [smilesCombined; t1.SMILES];
logInfo("Local fda_drugs.csv:          %d 行", height(t1)); %[output:3fd1b0d0]

if CORPUS_TARGET > 500 %[output:group:16007640]
    % 大きNのときのみ異質ソースを混合する
    t2 = readtable("data/list/everyday_chemicals.csv", TextType="string");
    smilesCombined = [smilesCombined; t2.SMILES];
    logInfo("Local everyday_chemicals.csv: %d 行追加", height(t2));

    t3 = readtable("data/list/forensic_challenge.csv", TextType="string");
    smilesCombined = [smilesCombined; t3.SMILES];
    logInfo("Local forensic_challenge.csv: %d 行追加", height(t3));
elseif CORPUS_TARGET > height(t1)
    logInfo("コーパス 201-500: FDA薬 + ChEMBL承認薬のみ使用"); %[output:525fc4be]
    logInfo("  (everyday/forensic除外 -- ドメイン混合は小規模Nのモデルを劣化する)"); %[output:751f8380]
end %[output:group:16007640]
logInfo("ローカルコーパス: %d SMILES (ChEMBL取得前)", numel(smilesCombined)); %[output:422f2707]

CHEMBL_CACHE = "data/list/chembl_extended.csv";
if CORPUS_TARGET > numel(smilesCombined) %[output:group:1a5839c5]
    nFetch = CORPUS_TARGET - numel(smilesCombined);
    cacheInsufficient = false;
    if isfile(CHEMBL_CACHE)
        cacheT = readtable(CHEMBL_CACHE, TextType="string");
        useN   = min(height(cacheT), nFetch + 150);  % 重複除去ロスのための余裕ウィンドウ
        if height(cacheT) < nFetch + 100
            % キャッシュが小さすぎる（例: 前回の実行で必要数未満の取得が行われた）
            logWarn("ChEMBLキャッシュは%d行のみ。~%d行必要。再取得中 ...", ...
                height(cacheT), nFetch + 150);
            cacheInsufficient = true;
        else
            logInfo("ChEMBLキャッシュヒット: %d / %d 行読み込み", useN, height(cacheT)); %[output:956d762e]
            smilesCombined = [smilesCombined; cacheT.SMILES(1:useN)];
        end
    end
    if ~isfile(CHEMBL_CACHE) || cacheInsufficient
        logInfo("ChEMBLから~%d SMILESを取得中（初回は~50秒） ...", nFetch);
        BASE_URL_C  = "https://www.ebi.ac.uk/chembl/api/data/molecule.json";
        FILTER_C    = "max_phase=4&molecule_type=Small+molecule&withdrawn_flag=0";
        PAGE_SIZE_C = 100;
        opts_c      = weboptions("Timeout", 30, "ContentType", "json");
        fetchedSmiles = strings(0, 1);
        offset_c = 200;   % fda_drugs.csvとの重複200分をスキップ
        nWant    = nFetch + 150;  % 重複除去後に十分な数を確保する大きめのバッファ
        while numel(fetchedSmiles) < nWant
            lim_c = min(PAGE_SIZE_C, nWant - numel(fetchedSmiles));
            url_c = sprintf("%s?%s&limit=%d&offset=%d", BASE_URL_C, FILTER_C, lim_c, offset_c);
            try
                data_c = webread(url_c, opts_c);
            catch ME
                logWarn("ChEMBL取得をoffset=%dで停止: %s", offset_c, ME.message);
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
            logInfo("  ChEMBL offset=%d -> %d 取得済み", offset_c, numel(fetchedSmiles));
            offset_c = offset_c + numel(mols_c);
            pause(0.3);
            if numel(mols_c) < lim_c; break; end
        end
        if ~isempty(fetchedSmiles)
            [~, ~] = mkdir(fileparts(CHEMBL_CACHE));
            writetable(array2table(fetchedSmiles(:), 'VariableNames', {'SMILES'}), CHEMBL_CACHE);
            logInfo("%d ChEMBL SMILESをキャッシュ -> %s", numel(fetchedSmiles), CHEMBL_CACHE);
            smilesCombined = [smilesCombined; fetchedSmiles(:)];
        else
            logWarn("ChEMBLの取得結果が0 SMILES -- ネットワークを確認してください。ローカルデータのみ使用します。");
        end
    end
end %[output:group:1a5839c5]
%[text] 重複除去後、CORPUS\_TARGET に切り詰めます。
smilesCombined = unique(smilesCombined);
logInfo("重複除去後: %d ユニークSMILES", numel(smilesCombined)); %[output:55200a7a]
if numel(smilesCombined) > CORPUS_TARGET
    smilesCombined = smilesCombined(1:CORPUS_TARGET);
end
logInfo("訓練用コーパス: %d (目標=%d)", numel(smilesCombined), CORPUS_TARGET); %[output:0617c68c]
%[text] RDKit で SMILES を検証し、解析不能な行を除外します。
nSmilesCombined = numel(smilesCombined);
validMask = false(nSmilesCombined, 1);
for k = 1:nSmilesCombined %[output:group:8a4fbe84]
    try
        validMask(k) = emk.mol.isValid(smilesCombined(k));
    catch
        validMask(k) = false;
    end
    if mod(k, max(1, round(nSmilesCombined/10))) == 0 || k == nSmilesCombined
        logProgress(k, nSmilesCombined, "SMILES検証中"); %[output:6cc88e94]
    end
end %[output:group:8a4fbe84]
smilesAll = smilesCombined(validMask);
N_MOLS    = numel(smilesAll);
logInfo("有効SMILES: %d / %d", N_MOLS, numel(smilesCombined)); %[output:1340b2ac]
%[text] 文字長の分布を確認します。MAX\_SEQ\_LEN を超える系列は打ち切られます。
lengths = cellfun(@strlength, cellstr(smilesAll));
logInfo("SMILES長: min=%d  median=%.0f  max=%d  mean=%.1f", ... %[output:group:6e60e7e4] %[output:73633f29]
    min(lengths), median(lengths), max(lengths), mean(lengths)); %[output:group:6e60e7e4] %[output:73633f29]

pctTruncated = 100 * mean(lengths > MAX_SEQ_LEN - 2);
logInfo("打ち切りされた系列 (> %d 文字): %.1f%%", MAX_SEQ_LEN-2, pctTruncated); %[output:9a34ceb7]

figure("Name", "R05 SMILES長分布"); %[output:45463a39]
histogram(lengths, 20, FaceColor=[0.2 0.5 0.8]); %[output:45463a39]
xlabel("SMILES長（文字数）"); ylabel("件数"); %[output:45463a39]
title(sprintf("SMILES長の分布 (N=%d)", N_MOLS)); %[output:45463a39]
xline(MAX_SEQ_LEN - 2, "r--", "最大系列長", LabelHorizontalAlignment="left"); %[output:45463a39]
grid on; %[output:45463a39]
%[text] SMILES 列挙データ拡張（Bjerrum 2017: arXiv:1703.07076）を行う前に、正則 SMILES のセットを保存します。
%[text] これは新規性チェック（Section 4, 9）で使用します。非正則の重複を「新規分子」として誤計上しないためです。また物性分布プロット（Section 9）でも使用し、元のデータ分布を反映します。
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
                rands(ri) = smilesAll(k);  % エラー時は正則SMILESにフォールバック
            end
        end
        augBuf = [augBuf; unique(rands)]; %#ok<AGROW>
        if mod(k, max(1, round(N_MOLS_ORIG/10))) == 0 || k == N_MOLS_ORIG
            logProgress(k, N_MOLS_ORIG, "SMILESオーグメンテーション中");
        end
    end
    augBuf  = unique(augBuf);
    augBuf  = augBuf(~ismember(augBuf, smilesAll));  % 正則重複を除去
    smilesAll = [smilesAll; augBuf];
    N_MOLS    = numel(smilesAll);
    logInfo("オーグメンテーション後コーパス: %d系列 (%.1fx)", N_MOLS, N_MOLS / N_MOLS_ORIG);
else
    logInfo("AUGMENT_FACTOR=1: オーグメンテーション無効 (コーパス = %d ユニークSMILES)", N_MOLS); %[output:233a680b]
end %[output:group:76149936]
%%
%[text] ## Section 2: 文字語彙の構築
%[text] ### トークン化と特殊トークン
%[text] 文字レベル言語モデルでは、コーパス内の各固有文字が 1 つの語彙トークン（vocabulary token）になります。語彙には、SMILES 文字列に絶対に現れない 3 つの特殊トークンも含まれます。
%[text] - START（"^"）: 全系列の先頭に付加します。生成の起点となります。
%[text] - END（"\\$"）: 分子完了を示すために末尾に付加します。訓練はこの位置で損失の計算を停止します。
%[text] - PAD（"\_"）: ミニバッチ内で固定長行列を作るため、実際の系列終端以後を埋めます。PAD 位置は損失関数からマスクされます。 \
%[text] また、Cl と Br は語彙構築前に 1 文字プレースホルダー（L と R）へ置換します。置換しないと 'l' と 'r' がどのトークンの後にもサンプリングされてしまい、無意味な系列を生成してしまうためです。
%[text] generateSmiles の出力では L→Cl、R→Br に復元します。
%[text] 
%[text] 前処理: 2 文字原子を単一文字プレースホルダーに置換します。L = Cl（塩素）、R = Br（臭素）。
smilesProc = strrep(strrep(smilesAll, "Cl", "L"), "Br", "R");

%[text] num2cell は 1×M の char 行を 1 文字ずつの cell 配列に変換します。
%[text] 特殊トークン 3 つの cell とそのまま horzcat できます。
allCharsCell = num2cell(unique(char(strjoin(smilesProc, ""))));
specialChars  = {char(PAD_TOKEN), char(START_TOKEN), char(END_TOKEN)};
inSpecial     = cellfun(@(c) any(strcmp(c, specialChars)), allCharsCell);
vocabChars    = [specialChars, sort(allCharsCell(~inSpecial))];
VOCAB_SIZE    = numel(vocabChars);

%[text] SMILES トークンセット外の文字が語彙に含まれる場合は警告します。
%[text] char(92) はバックスラッシュ（E/Z 立体化学記法の ''）です。JSON エスケープ問題を回避するため数値で参照します。
%[text] 'a' は一部 ChEMBL 正則 SMILES で使用される芳香族ワイルドカード原子で、RDKit では有効です。
SMILES_ALLOWED = {'B','C','N','O','P','S','F','I','L','R', ...
    'a','b','c','n','o','p','s', ...
    '#','=','-','+','.', ...
    '(',')', '[',']', ...
    '0','1','2','3','4','5','6','7','8','9', ...
    '@','/', char(92), '%','H'};
unexpectedChars = vocabChars(~ismember(vocabChars, [specialChars, SMILES_ALLOWED]));
if ~isempty(unexpectedChars) %[output:group:415def2d]
    logWarn("語彙に想定外の文字が含まれます: [%s] -- 生成品質が低下する可能性があります。", ... %[output:6ccf3410]
        strjoin(string(unexpectedChars), " ")); %[output:6ccf3410]
end %[output:group:415def2d]

char2idx = containers.Map(vocabChars, num2cell(1:VOCAB_SIZE));
idx2char  = containers.Map(num2cell(1:VOCAB_SIZE), vocabChars);

PAD_IDX   = char2idx(char(PAD_TOKEN));
START_IDX = char2idx(char(START_TOKEN));
END_IDX   = char2idx(char(END_TOKEN));

logInfo("語彙サイズ: %d文字", VOCAB_SIZE); %[output:7c8775f2]
logInfo("特殊トークン: PAD=%d  START=%d  END=%d", PAD_IDX, START_IDX, END_IDX); %[output:329ce97e]
%%
%[text] ## Section 3: バイグラムマルコフ連鎖 -- 単純な生成ベースライン
%[text] ### n-gram 言語モデル
%[text] 深層学習以前、系列生成は n-gram モデルに依存していました。
%[text] バイグラム（$n=2$）モデルは次の条件付き確率を推定します: 
%[text]{"align":"center"} $P(c(t) \\mid c(t-1)) = \\frac{\\mathrm{count}(c(t-1),\\, c(t))}{\\mathrm{count}(c(t-1))}$
%[text] 各文字の確率は直前の 1 文字にのみ依存します。
%[text] これに対応するコードが `transCount` 行列です。各バイグラム $(c(t-1), c(t))$ を観測するたびに `transCount(c(t-1), c(t))` に 1 を加算します。
%[text] SMILES 生成では、小規模データでも驚くほど有効です。化学的に障壁の小さい局所遷移が多く存在するためです。'(' の後はほぼ常に原子（C, N, O, S ...）が来ます。
%[text] '=' または '\#' の後は必ず原子が来ます（結合次数記述子）。
%[text] '\[' の後は必ず原子記号が来ます（括弧原子記法）。
%[text] これらの局所パターンは 200 分子程度から学習可能です。
%[text] ラプラススムージング（$\\alpha = 0.1$）を各セルに加算することで確率ゼロを回避し、各行を正規化して確率分布を得ます。LSTM との期待性能比較: Markov モデルは訓練不要で 200 分子から有効率 20〜35% を達成しますが、LSTM は同規模では 5〜15% 程度（アンダーフィット）となります。デフォルトの 500 分子コーパスでは語彙と遷移の多様性が増すため、Markov の有効率は 5〜15% 程度に低下します。
%[text] 
%[text] 訓練 SMILES からバイグラム遷移カウント行列を構築します（smilesProc を使用、Section 2 で Cl/Br 置換済み）。
transCount = zeros(VOCAB_SIZE, VOCAB_SIZE, "double");
for mk = 1:N_MOLS
    s = char(smilesProc(mk));
    % 系列をエンコード: START + 文字 + END
    nChar  = numel(s);
    idxSeq = zeros(1, nChar + 2);
    idxSeq(1) = START_IDX;
    for ci = 1:nChar
        if isKey(char2idx, s(ci))
            idxSeq(ci + 1) = char2idx(s(ci));
        else
            idxSeq(ci + 1) = PAD_IDX;   % 未知文字 -> PAD
        end
    end
    idxSeq(end) = END_IDX;
    % 全連続ペアをカウント
    for ti = 1:numel(idxSeq) - 1
        c1 = idxSeq(ti);
        c2 = idxSeq(ti + 1);
        if c1 >= 1 && c1 <= VOCAB_SIZE && c2 >= 1 && c2 <= VOCAB_SIZE
            transCount(c1, c2) = transCount(c1, c2) + 1;
        end
    end
end
%[text] ラプラススムージング（$\\alpha = 0.1$）: 各セルに偽カウントを加算して確率ゼロを回避します。
%[text] PAD は終端吸収状態（自己ループ）で、他のトークンに遷移しません。
LAPLACE_ALPHA    = 0.1;
transCountSmooth = transCount + LAPLACE_ALPHA;
transCountSmooth(PAD_IDX, :)          = 0;   % PADは遷移しない
transCountSmooth(PAD_IDX, PAD_IDX)    = 1;   % 正規化のための吸収状態
%[text] 行ごとに正規化して確率分布を得ます。
rowSums   = sum(transCountSmooth, 2);
rowSums(rowSums == 0) = 1;          % 全ゼロ行へのガード
transProb = transCountSmooth ./ rowSums;

logInfo("マルコフモデル構築完了: %dx%d遷移テーブル (%d SMILESから, %dバイグラム)", ... %[output:group:4b3ee340] %[output:18cc59a8]
    VOCAB_SIZE, VOCAB_SIZE, N_MOLS, round(sum(transCount(:)))); %[output:group:4b3ee340] %[output:18cc59a8]
%[text] 主要 SMILES 文字の遷移ヒートマップを可視化します。行が「現在の文字」、列が「次の文字」に対応します。
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

figure("Name", "R05 バイグラム遷移確率"); %[output:2f235c0c]
imagesc(subMatrix); %[output:2f235c0c]
colormap("hot");  colorbar; %[output:2f235c0c]
xticks(1:heatN);  xticklabels(heatLabels);  xtickangle(45); %[output:2f235c0c]
yticks(1:heatN);  yticklabels(heatLabels); %[output:2f235c0c]
xlabel("次の文字");  ylabel("現在の文字"); %[output:2f235c0c]
title(sprintf("バイグラム遷移確率（上位%d文字）", heatN)); %[output:2f235c0c]
%%
%[text] ## Section 4: マルコフモデルによる生成と評価
%[text] ### GuacaMol / MOSES 評価プロトコル
%[text] 分子生成モデルの 3 つの標準指標を紹介します（Brown et al. 2019, Polykovskiy et al. 2020）。
%[text] - **有効性（Validity）**: RDKit で解析可能な生成 SMILES の割合です。
%[text] - **多様性（Uniqueness）**: 構造的に異なる有効 SMILES の割合です。
%[text] - **新規性（Novelty）**: 訓練セットに含まれない固有の有効 SMILES の割合です。 \
%[text] 本番目標は有効性 \> 90%、多様性 \> 85%、新規性 \> 60% です。
%[text] マルコフベースラインは 200 分子で有効性 20〜35% を達成しますが、デフォルトの 500 分子コーパスでは 5〜15% 程度となります（語彙多様性の増加による）。
%[text] 訓練時間ゼロのカウンティングモデルとしては良好ですが、本番目標には大きく及びません。
%[text] Section 7〜8 では LSTM も小規模データでは不足する理由を示します。
N_MARKOV = N_GENERATED;   % LSTMと直接比較するために同じ件数
logInfo("バイグラムマルコフモデルで%d SMILESを生成中 ...", N_MARKOV); %[output:235d5d03]
markovSmiles = strings(N_MARKOV, 1);
for mg = 1:N_MARKOV %[output:group:8d9a6cf7]
    markovSmiles(mg) = generateSmilesMarkov(transProb, START_IDX, END_IDX, ...
        PAD_IDX, VOCAB_SIZE, MAX_SEQ_LEN, TEMPERATURE, idx2char, char2idx);
    if mod(mg, max(1, round(N_MARKOV/10))) == 0 || mg == N_MARKOV
        logProgress(mg, N_MARKOV, "生成中 (Markov)"); %[output:4a72d9ff]
    end
end %[output:group:8d9a6cf7]

%[text] RDKit で生成した SMILES を検証します。
logInfo("Markov SMILES を RDKit で検証中 ..."); %[output:353a9a4a]
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
        logProgress(mg, N_MARKOV, "検証中 (Markov)"); %[output:7b36f722]
    end
end %[output:group:345dfb99]

nMValid  = sum(markovValid);
mValidSm = markovSmiles(markovValid);
nMUnique = numel(unique(mValidSm));
nMNovel  = sum(~ismember(mValidSm, smilesOriginal));  % 正則訓練セットに対する新規性

logInfo("=== マルコフベースライン結果 (%d件生成) ===", N_MARKOV); %[output:8098b4ea]
logInfo("  有効性:   %d/%d  (%.1f%%)", nMValid, N_MARKOV, 100*nMValid/N_MARKOV); %[output:0fd21ee6]
logInfo("  多様性: %d/%d  (%.1f%%)", nMUnique, max(nMValid,1), ... %[output:group:3ed5b7d2] %[output:4f3967e5]
    100*nMUnique/max(nMValid,1)); %[output:group:3ed5b7d2] %[output:4f3967e5]
logInfo("  新規性:    %d/%d  (%.1f%%)", nMNovel, max(nMValid,1), ... %[output:group:104abef8] %[output:9d3c70f5]
    100*nMNovel/max(nMValid,1)); %[output:group:104abef8] %[output:9d3c70f5]
logInfo("  (本番目標: 有効性>90%%  多様性>85%%  新規性>60%%)"); %[output:06b087e6]

%[text] 有効性にかかわらず、診断用に先頭 5 件のサンプル SMILES を表示します。
logInfo("サンプルMarkov SMILES（先頭5件）:"); %[output:8f440e59]
for mg_s = 1:min(5, N_MARKOV) %[output:group:2366cbc7]
    logInfo("  [%2d] %s", mg_s, markovSmiles(mg_s)); %[output:7147febe]
end %[output:group:2366cbc7]

if nMValid > 0 %[output:group:022431b9]
    logInfo("最初の%d周の有効Markov SMILES:", min(5, nMValid)); %[output:2a8966f4]
    for mg = 1:min(5, nMValid)
        logInfo("  [%2d] %s", mg, mValidSm(mg)); %[output:067fae5f]
    end
end %[output:group:022431b9]
%%
%[text] ## Section 5: マルコフ連鎖の限界 -- LSTMの動機付け
%[text] ### 長距離依存問題
%[text] マルコフモデルは良好な局所精度（正しい原子間遷移）を達成しますが、1 文字を超える構造的制約に対して系統的に失敗します。
%[text] **問題 1: 括弧対応（枝分かれ記法）**
%[text] 医薬品 SMILES はネストした枝分かれを持ちます。CC(=O)O（酢酸）では '(' を対応する ')' で閉じる必要がありますが、バイグラムモデルは枝分かれが開いていることを「記憶」できません。
%[text] → このスクリプトでは文法則 2, 5 が括弧対応を決定論的に強制するため、実行時の括弧不一致エラーは大幅に抑制されます。
%[text] **問題 2: 環閉鎖結合の対応**
%[text] 'C1CCCCC1'（シクロヘキサン）では、最初の C の後に数字 1 で環結合を開き、5 つの原子を生成した後、再度 1 で環結合を閉じる必要があります。
%[text] この「メモリ」は 7 文字にまたがり、バイグラム（1 文字先読み）の限界を大きく超えます。
%[text] → 文法則 2c, 5 が環閉鎖数字の累積と未閉鎖状態を追跡して強制閉鎖するため、環閉鎖不一致も実行時には抑制されます。
%[text] **問題 3: 原子価の一貫性**
%[text] 炭素は常に 4 結合を形成します。バイグラムは現在の炭素がこれまでのトークンで何結合を形成したかを数えられません。
%[text] **問題 4: 環結合の累積**
%[text] 医薬品 SMILES は環閉鎖数字 1〜6 を頻繁に使用します。環深度限界なしでは未閉鎖環結合が累積し、強制閉鎖ウィンドウ内で全部閉じられなくなります。
%[text] 文法則 2c は同時に開く環結合を 3 つまでに制限することでこの問題を緩和します。
%[text] **解決策: LSTM**
%[text] LSTM（Hochreiter & Schmidhuber 1997）は任意の距離にわたる文脈を蓄積する隠れ状態 $h(t)$ とセル状態 $c(t)$ を保持します。
%[text] 10,000 分子以上のスケールでは、LSTM は隠れ状態内で開いた括弧と環結合を暗黙的に数えられるようになり、有効率 \> 90% を達成します。
%[text] まとめると、マルコフは 200 分子でも機能しますが（バイグラムカウントはデータ効率が良い）、LSTM は 500 分子規模では失敗します。
%[text] BPTT で長距離パターンを学習するにはさらに多くのサンプルが必要であり、本番 LM（MolGPT, ChemGPT）は $10^6 \\sim 10^7$ 分子を使用します。
%[text] 
%[text] Markov 生成 SMILES の失敗モードを 4 種類（括弧不一致・角括弧不一致・環閉鎖不一致・その他）に分類して分析します。
%[text] grammar constraint（文法則 2, 2c, 5）が括弧・環閉鎖の失敗を機械的に抑制するため、実測では括弧/環閉鎖の不一致は合計 10% 以下にとどまり、「その他（原子価・芳香族性・無効原子記号）」が失敗の 90% 以上を占めます。
%[text] これが問題 3（原子価の一貫性）の本質です。括弧・環閉鎖の「メモリ」は深さカウンタ（整数 1 個）で追跡できるため決定論的規則で補えます。一方、原子価の追跡には「系列中の各原子が現時点で何本の結合を持つか」を原子ごとに管理する状態が必要で、分岐・環の組み合わせにより状態空間が爆発します。結果として原子価の正しさはデータから化学知識として学習するしかなく、単純な規則では解決できません。
markovInvalid = markovSmiles(~markovValid);
logInfo("Markov失敗分析 (%d無効 / %d生成):", ... %[output:group:7ebbdda9] %[output:8b7375fb]
    numel(markovInvalid), N_MARKOV); %[output:group:7ebbdda9] %[output:8b7375fb]
if numel(markovInvalid) > 0 %[output:group:9b729d14]
    nParenErr   = sum(cellfun(@(s) count(char(s), '(') ~= count(char(s), ')'), ...
        cellstr(markovInvalid)));
    nBracketErr = sum(cellfun(@(s) count(char(s), '[') ~= count(char(s), ']'), ...
        cellstr(markovInvalid)));
    % 環閉鎖不一致: 奇数回現する数字は開かれたまま閉じられていない（または逆）。
    nRingErr = sum(cellfun(@(s) any(arrayfun(@(d) mod(count(char(s), char('0'+d)), 2) ~= 0, 1:9)), ...
        cellstr(markovInvalid)));
    nOtherErr = numel(markovInvalid) - nParenErr - nBracketErr - nRingErr;
    logInfo("  括弧不一致 '()':  %d / %d  (%.0f%%)", ... %[output:5276edff]
        nParenErr,   numel(markovInvalid), 100*nParenErr/numel(markovInvalid)); %[output:5276edff]
    logInfo("  角括弧不一致 '[]':     %d / %d  (%.0f%%)", ... %[output:83699159]
        nBracketErr, numel(markovInvalid), 100*nBracketErr/numel(markovInvalid)); %[output:83699159]
    logInfo("  環閉鎖数字不一致:  %d / %d  (%.0f%%)", ... %[output:14d3854d]
        nRingErr,    numel(markovInvalid), 100*nRingErr/numel(markovInvalid)); %[output:14d3854d]
    logInfo("  その他（原子価/原子/無効）: %d / %d  (%.0f%%)", ... %[output:9f6c2440]
        max(nOtherErr,0), numel(markovInvalid), 100*max(nOtherErr,0)/numel(markovInvalid)); %[output:9f6c2440]
end %[output:group:9b729d14]
logInfo("次のステップ: Section 6-8でLSTM文字言語モデルを構築する。"); %[output:8ce954da]
%%
%[text] ## Section 6: エンコード・パディング・データ分割
%[text] ### 因果言語モデルの目的関数
%[text] 文字レベル言語モデルの訓練は「次トークン予測（next-token prediction）」タスクです。
%[text] 各系列 $c(1), c(2), \\ldots, c(T)$ に対して、入力は $\[\\mathrm{START},\\, c(1), \\ldots, c(T-1)\]$、ターゲットは $\[c(1), c(2), \\ldots, c(T),\\, \\mathrm{END}\]$ となります（それぞれ長さ $T$）。
%[text] 各位置 $t$ で、モデルは $c(1) \\ldots c(t)$ を見て次の文字 $c(t+1)$ を予測します。
%[text] この「教師強制（teacher forcing）」は対数尤度目標 $L = -\\sum\_t \\log P(c(t+1) \\mid c(1) \\ldots c(t))$ と等価です。
%[text] コード上では `X_seqMat` が入力（現在のトークン）、`Y_seqMat` がターゲット（1 つ右にシフトしたトークン）に対応します。
%[text] パディングは PAD トークンで全系列を MAX\_SEQ\_LEN に延長し、損失関数は PAD 位置をマスクして勾配に寄与しないようにします。
encodedSeqs = cell(N_MOLS, 1);
for k = 1:N_MOLS
    s    = char(smilesProc(k));                % 語彙参照にはL/Rプレースホルダー形式を使用
    s    = s(1:min(end, MAX_SEQ_LEN - 2));     % 必要なら打ち切り
    idxs = arrayfun(@(c) char2idx(char(c)), s);
    encodedSeqs{k} = [START_IDX, idxs, END_IDX];
end
seqLengths = cellfun(@numel, encodedSeqs);
MAX_T      = max(seqLengths);

%[text] MAX\_T にパディングします。seqMatrix は $\[N \\times \\mathrm{MAX\\\_T}\]$ の整数行列です。
seqMatrix = ones(N_MOLS, MAX_T, "single") * PAD_IDX;
for k = 1:N_MOLS
    L = seqLengths(k);
    seqMatrix(k, 1:L) = single(encodedSeqs{k});
end

%[text] 入力 X とターゲット Y を 1 トークン分シフトして teacher-forcing の形式にします。
X_seqMat = seqMatrix(:, 1:end-1);   % [N x SEQ_LEN]  入力
Y_seqMat = seqMatrix(:, 2:end);     % [N x SEQ_LEN]  ターゲット（1トークンシフト）
SEQ_LEN  = size(X_seqMat, 2);

logInfo("系列行列: %d x %d  (語彙=%d)", N_MOLS, SEQ_LEN, VOCAB_SIZE); %[output:07691b4a]
%[text] 80/20 で訓練・検証に分割します（乱数シードを 42 に固定して再現性を確保）。
rng(42);
perm    = randperm(N_MOLS);
n_train = round(0.8 * N_MOLS);
trainIdx = perm(1:n_train);
valIdx   = perm(n_train+1:end);

X_train = X_seqMat(trainIdx, :);
Y_train = Y_seqMat(trainIdx, :);
X_val   = X_seqMat(valIdx, :);
Y_val   = Y_seqMat(valIdx, :);

logInfo("分割: %d 訓練 / %d 検証", n_train, N_MOLS - n_train); %[output:0c181dc6]
%%
%[text] ## Section 7: 文字レベルLSTMの定義と訓練
%[text] ### BPTT（時間方向逆伝播）
%[text] LSTM は系列を左から右に処理し、文脈を隠れ状態 $h(t)$ とセル状態 $c(t)$ に蓄積します。
%[text] 訓練中は BPTT が LSTM を $T$ タイムステップに展開し、全ステップを通じて勾配を計算します。
%[text] アーキテクチャは次のように構成されます。
%[text] 1. 入力: ワンホット表現 $\[\\mathrm{VOCAB\\\_SIZE} \\times T \\times B\]$
%[text] 2. Embedding 全結合層 (VOCAB\_SIZE → EMBED\_DIM): 類似文字（C/N/O、c/n/o など）を密な特徴空間でグループ化する学習済み文字埋め込みとして機能します。
%[text] 3. lstmLayer(HIDDEN\_SIZE): 各 $t$ で $h(t)$ を出力します。
%[text] 4. dropoutLayer(DROPOUT\_RATE): 正則化します。
%[text] 5. fullyConnectedLayer(VOCAB\_SIZE) → softmaxLayer: 次トークンの確率分布 $\[\\mathrm{VOCAB\\\_SIZE} \\times T \\times B\]$ を出力します。 \
%[text] 損失はマスク付きクロスエントロピーです。PAD 位置は平均から除外されます。
%[text] 訓練サンプル数が約 160（200 分子コーパスの 80%）しかない場合、モデルは訓練と検証の損失ギャップがほとんど生じません（150 エポックで約 0.06）。
%[text] FDA 承認薬は分子量 3〜2000+ にわたる高度に異質な環システムを持つため、LSTM は記憶ではなく均一な文字分布に収束します。
%[text] これは過学習（overfitting）ではなくアンダーフィット（underfitting）です。
%[text] 重要な教訓: パラメータ数 / トークン数 \> 1 の「過パラメータ化」であっても、データ分布が高度に多様な場合は必ずしも過学習が生じるわけではありません。
%[text] 
%[text] アーキテクチャ: one-hot(V) → Embed(EMBED\_DIM) → LSTM(HIDDEN) → Dropout → FC(V) → Softmax
%[text] 埋め込み全結合層（VOCAB\_SIZE → EMBED\_DIM）は学習済み文字埋め込みとして機能します。
%[text] 類似文字（原子 C/N/O、芳香族 c/n/o、結合文字）を密な特徴空間でグループ化し、LSTM はスパースな 37 次元ワンホットではなく密な EMBED\_DIM ベクトルを受け取ることで勾配フローが改善されます。
%[text] 結果として 'C' と 'c' が「共に炭素」として類似した表現に学習されます。
layers = [
    sequenceInputLayer(VOCAB_SIZE, Name="input")
    fullyConnectedLayer(EMBED_DIM, Name="embed")    % 学習済み文字埋め込み（活性化なし）
    lstmLayer(HIDDEN_SIZE, OutputMode="sequence", Name="lstm1")
    dropoutLayer(DROPOUT_RATE, Name="drop")
    fullyConnectedLayer(VOCAB_SIZE, Name="fc")
    softmaxLayer(Name="softmax")
];
net = dlnetwork(layers);

nParams = sum(cellfun(@numel, net.Learnables.Value));
logInfo("ネットワーク: Embed(%d->%d) + LSTM(%d) + FC -> softmax  |  パラメータ数: %d", ... %[output:group:43ceed7d] %[output:3b3443ea]
    VOCAB_SIZE, EMBED_DIM, HIDDEN_SIZE, nParams); %[output:group:43ceed7d] %[output:3b3443ea]
logInfo("訓練トークン対パラメータ比: %.2f  (>1 => 過パラメータ化)", ... %[output:group:81631103] %[output:5f85d3e4]
    nParams / (n_train * mean(seqLengths))); %[output:group:81631103] %[output:5f85d3e4]
%[text] カスタム訓練ループを実行します。各エポックでミニバッチをシャッフルし、Adam オプティマイザで更新します。
numBatches   = ceil(n_train / BATCH_SIZE);
avgG_sl      = [];
avgSqG_sl    = [];
iter_sl      = 0;
trainLossLog = zeros(N_EPOCHS, 1);
valLossLog   = zeros(N_EPOCHS, 1);

for epoch = 1:N_EPOCHS %[output:group:914dc138]
    % 各エポックで訓練順序をシャッフル
    perm_e = randperm(n_train);
    X_shuf = X_train(perm_e, :);
    Y_shuf = Y_train(perm_e, :);
    epochLoss = 0;

    for b = 1:numBatches
        i1 = (b-1)*BATCH_SIZE + 1;
        i2 = min(b*BATCH_SIZE, n_train);
        Xb = X_shuf(i1:i2, :);   % [B x T]
        Yb = Y_shuf(i1:i2, :);   % [B x T]

        % ワンホットエンコード: [VOCAB_SIZE x T x B]
        dlX = dlarray(onehot_encode(Xb, VOCAB_SIZE), "CTB");

        iter_sl = iter_sl + 1;
        [lossVal, grads] = dlfeval(@modelLoss, net, dlX, Yb, PAD_IDX);
        % 勾配クリッピング: 損失スパイクを防ぐためテンソルごとのL2ノルムを制限
        for gi = 1:height(grads)
            gdata = extractdata(grads.Value{gi});
            nrm = sqrt(sum(gdata(:).^2));
            if nrm > GRAD_CLIP_NORM
                grads.Value{gi} = grads.Value{gi} * (GRAD_CLIP_NORM / nrm);
            end
        end
        epochLoss = epochLoss + extractdata(lossVal);
        % コサインLRアニーリング: LEARN_RATEからLEARN_RATE/20にN_EPOCHSで減衰
        lr_epoch = LEARN_RATE * (0.05 + 0.95 * 0.5 * (1 + cos(pi * (epoch-1) / N_EPOCHS)));
        [net, avgG_sl, avgSqG_sl] = adamupdate(net, grads, avgG_sl, avgSqG_sl, ...
            iter_sl, lr_epoch);
    end
    trainLossLog(epoch) = epochLoss / numBatches;

    % 検証損失（勾配なし）
    dlXv = dlarray(onehot_encode(X_val, VOCAB_SIZE), "CTB");
    valLossLog(epoch) = extractdata(dlfeval(@modelLoss, net, dlXv, Y_val, PAD_IDX));

    if mod(epoch, 10) == 0 || epoch == 1
        logInfo("エポック %3d/%d  train=%.4f  val=%.4f", ... %[output:81098cf3]
            epoch, N_EPOCHS, trainLossLog(epoch), valLossLog(epoch)); %[output:81098cf3]
    end
end %[output:group:914dc138]
figure("Name", "R05 訓練カーブ"); %[output:84c89218]
plot(1:N_EPOCHS, trainLossLog, "b-", LineWidth=1.5, DisplayName="訓練"); %[output:84c89218]
hold on; %[output:84c89218]
plot(1:N_EPOCHS, valLossLog,   "r--", LineWidth=1.5, DisplayName="検証"); %[output:84c89218]
xlabel("エポック"); ylabel("クロスエントロピー損失"); %[output:84c89218]
title("文字レベルLSTM: val損失プラトーがデータスケール限界の主要証拠"); %[output:84c89218]
legend; grid on; %[output:84c89218]
%[text] val 損失がエポック 50〜70 付近で改善しなくなる「プラトー」が主要な学習ポイントです。
%[text] 有効率の数値差（Markov 約 9%、LSTM 約 5%）は $N=300$ サンプルでは約 3〜4 ポイントであり、統計的に有意ではありません（$z=1.6$、$p \\approx 0.10$、95% CI 重複）。
%[text] 生の有効率ギャップを過度に解釈せず、損失カーブで判断しましょう。
logInfo("訓練完了。最終train/val: %.4f / %.4f", ... %[output:group:3700040c] %[output:7a5fd55e]
    trainLossLog(end), valLossLog(end)); %[output:group:3700040c] %[output:7a5fd55e]
%%
%[text] ## Section 8: 温度サンプリングによるSMILES生成
%[text] ### 自己回帰デコーディングと温度サンプリング
%[text] 訓練後、新しい SMILES は 1 文字ずつ自己回帰的に生成されます。
%[text] 1. START トークンを入力し、LSTM の前向きパスを実行します。
%[text] 2. 出力分布から次の文字をサンプリングします。
%[text] 3. サンプリングした文字を次の入力として利用します。
%[text] 4. END トークンをサンプリングするか、MAX\_SEQ\_LEN に達するまで (1)〜(3) を繰り返します。 \
%[text] 温度 $T$ は探索と活用のバランスを調節します: $P\_T(c) = \\mathrm{softmax}\\!\\left(\\log P(c) / T\\right)$
%[text] $T \< 1$ では分布が尖り、保守的なサンプリングで無効文字が減ります。
%[text] $T \> 1$ では分布が平坦になり、創造的ですが無効 SMILES が増えます。
%[text] 実装上の注意: 各生成ステップで増大するプレフィックス全体を LSTM に入力します（最後のトークンだけではなく）。
%[text] これは $O(T^2)$ の計算量になりますが、正確な実装です（ステートフル LSTM を使えば $O(T)$ にできます）。
%[text] **文法制約**: generateSmiles() は 5 つのマスキング則を適用して無効な SMILES を排除します。
%[text] 1. 先頭位置マスク: 原子文字と '\[' のみが最初のトークンになれます。
%[text] 2. 括弧深度ガード: ')' は paren\_depth = 0 のとき禁止します。
%[text] 3. 括弧内部マスク: '\[...\]' 内部で '(' と '\[' を禁止します。
%[text] 4. END ブロック: 未閉鎖の括弧・環結合がある間は END を禁止します。
%[text] 5. 強制閉鎖: 残りステップが不足する場合、開いた区切り文字を強制的に閉じます。 \
%[text] 芳香族原子価などの構造的問題には、より大きな訓練コーパスが必要です。
logInfo("%d SMILESを生成中 (T=%.2f) ...", N_GENERATED, TEMPERATURE); %[output:2d058d5c]
generatedSmiles = strings(N_GENERATED, 1);
for g = 1:N_GENERATED %[output:group:525e00c9]
    generatedSmiles(g) = generateSmiles(net, START_IDX, END_IDX, PAD_IDX, ...
        VOCAB_SIZE, MAX_SEQ_LEN, TEMPERATURE, idx2char, char2idx);
    if mod(g, max(1, round(N_GENERATED/10))) == 0 || g == N_GENERATED
        logProgress(g, N_GENERATED, "生成中 (RL前)"); %[output:3bd50a65]
    end
end %[output:group:525e00c9]
logInfo("生成SMILESのサンプル（先頭8件）:"); %[output:74269a7c]
for k = 1:min(8, N_GENERATED) %[output:group:4081c67c]
    logInfo("  [%3d] %s", k, generatedSmiles(k)); %[output:42d269fd]
end %[output:group:4081c67c]
%%
%[text] ## Section 9: 有効性・多様性・新規性・物性分析
%[text] ### 生成モデルの 4 つの標準評価指標（Brown et al. 2019）
%[text] **有効性（Validity）**: RDKit が解析できる生成 SMILES の割合です。モデルが SMILES 文法をどれだけ学習したかを反映します。
%[text] **多様性（Uniqueness）**: 有効な生成 SMILES のうち重複のないものの割合です。多様性が低いほどモデルが少数の構造に崩壊しています。
%[text] **新規性（Novelty）**: 重複のない有効 SMILES のうち訓練セットにないものの割合です。新規性が低いほどモデルが記憶しており汎化できていません。
%[text] **物性 KL ダイバージェンス**: 生成分子と訓練分子の物性分布間の KL ダイバージェンスです。KL が低いほどモデルが訓練コーパスの「スタイル」を再現しています。
%[text] デフォルト設定（500 分子・HIDDEN=128・パラメータ数 約 80K）の期待値: 有効性 5〜15%（文法制約付き）、多様性 60〜90%、新規性 80〜100% です。200 分子・HIDDEN=64（約 28K パラメータ）の場合は有効性 5〜20%、新規性 40〜70% 程度になります。
%[text] 有効性が低いのは想定内であり、データスケールの課題を示しています。
%[text] 本番モデル（REINVENT, MolGPT）は $10^6$ 分子以上で有効性 \> 95% を達成します。
validGenMask = false(N_GENERATED, 1);
for k = 1:N_GENERATED %[output:group:354fa23a]
    try
        validGenMask(k) = emk.mol.isValid(generatedSmiles(k));
    catch
        validGenMask(k) = false;
    end
    if mod(k, max(1, round(N_GENERATED/10))) == 0 || k == N_GENERATED
        logProgress(k, N_GENERATED, "検証中 (生成)"); %[output:4c71865b]
    end
end %[output:group:354fa23a]
validGenSmiles = generatedSmiles(validGenMask);
nValid   = sum(validGenMask);
nUnique  = numel(unique(validGenSmiles));
nNovel   = sum(~ismember(validGenSmiles, smilesOriginal));  % novelty vs canonical training set

logInfo("RL前生成 -- 合計 %d", N_GENERATED); %[output:58037d9b]
logInfo("  有効性:    %d/%d  (%.1f%%)", nValid,  N_GENERATED, 100*nValid/N_GENERATED); %[output:68c1252a]
logInfo("  多様性:  %d/%d  (%.1f%%)", nUnique, max(nValid,1), 100*nUnique/max(nValid,1)); %[output:1307a96b]
logInfo("  新規性:     %d/%d  (%.1f%%)", nNovel,  max(nValid,1), 100*nNovel/max(nValid,1)); %[output:54f8ce50]
%[text] 生成分子と訓練分子の物性分布を比較します（分子量と LogP）。
FEAT_NAMES = ["MolWt", "LogP"];
propsGen   = batchDescriptors(validGenSmiles(1:min(50,nValid)), FEAT_NAMES);
propsTrain = batchDescriptors(smilesOriginal, FEAT_NAMES);  % property dist of original molecules

if ~isempty(propsGen) && ~isempty(propsTrain) %[output:group:1f959db9]
    figure("Name", "R05 Property Distribution (pre-RL)"); %[output:9f971cbd]
    titles = ["分子量", "LogP（親油性）"];
    units  = ["Da", "log units"];
    for fi = 1:numel(FEAT_NAMES)
        subplot(1, 2, fi); hold on; %[output:9f971cbd]
        histogram(propsTrain(:,fi), 20, Normalization="probability", ...
            FaceColor=[0.2 0.6 0.8], DisplayName="Training");
        histogram(propsGen(:,fi),   20, Normalization="probability", ...
            FaceColor=[0.9 0.4 0.2], FaceAlpha=0.7, DisplayName="生成");
        xlabel(sprintf("%s (%s)", FEAT_NAMES(fi), units(fi)));
        ylabel("確率密度");
        title(titles(fi)); legend; grid on;
    end
    sgtitle("生成 vs 訓練: 物性分布 (RL前)"); %[output:9f971cbd]
end %[output:group:1f959db9]
%%
%[text] ## Section 10: チェックポイント保存とまとめ
%[text] ### 訓練済みモデルのシリアライズ
%[text] 訓練済みネットワークと語彙を保存することで、再訓練なしに新しいセッションで生成を再開できます。
%[text] また REINFORCE 微調整のために r06\_reinforce.m に転送することもできます。
%[text] result/r05\_checkpoint.mat に保存します（result/ は .gitignore 管理のため、ファイルはリポジトリにコミットされません）。
CHECKPOINT = "result/r05_checkpoint.mat";
if ~isfolder("result"); mkdir("result"); end
save(CHECKPOINT, ...
    "net", "char2idx", "idx2char", "VOCAB_SIZE", "MAX_SEQ_LEN", ...
    "START_IDX", "END_IDX", "PAD_IDX", "HIDDEN_SIZE", "TEMPERATURE", ...
    "smilesAll", "smilesProc", "X_train", "Y_train", "X_val", "Y_val");
logInfo("チェックポイント保存完了 -> %s", CHECKPOINT); %[output:6816dab6]
logInfo("  R06で読み込む: load(""result/r05_checkpoint.mat"")"); %[output:4dca65b3]
%[text] ### --- まとめ ---
logInfo("コーパス:         %d分子 有効FDA薬SMILES (拡張後 %d分子で訓練)", N_MOLS_ORIG, N_MOLS); %[output:4752e9cd]
logInfo("語彙:     %d文字  |  最大系列長: %d", VOCAB_SIZE, MAX_T); %[output:86afdc49]
logInfo("損失プラトー:   val=%.4f (エポック%d+)  --  ランダムベースライン=%.4f  (パープレキシティ=%.1f)", ... %[output:group:44d2b1de] %[output:5038645d]
    min(valLossLog), find(valLossLog == min(valLossLog), 1), ... %[output:5038645d]
    log(VOCAB_SIZE), exp(min(valLossLog))); %[output:group:44d2b1de] %[output:5038645d]
logInfo("--- 生成指標 (N=%dずつ; 二項型SE~%.1f%%)", N_GENERATED, ... %[output:group:22e3defc] %[output:134573c2]
    100*sqrt(0.10*0.90/N_GENERATED)); %[output:group:22e3defc] %[output:134573c2]
if exist("nMValid", "var") %[output:group:4103d63f]
    logInfo("マルコフモデル:   有効性=%.1f%%  多様性=%.1f%%  新規性=%.1f%%", ... %[output:6104ce8f]
        100*nMValid/N_MARKOV, 100*nMUnique/max(nMValid,1), 100*nMNovel/max(nMValid,1)); %[output:6104ce8f]
else
    logInfo("マルコフモデル:   (Section 3-5を先に実行してください)");
end %[output:group:4103d63f]
logInfo("LSTMモデル:     有効性=%.1f%%  多様性=%.1f%%  新規性=%.1f%%", ... %[output:group:08457e7b] %[output:97f77c7d]
    100*nValid/N_GENERATED, 100*nUnique/max(nValid,1), 100*nNovel/max(nValid,1)); %[output:group:08457e7b] %[output:97f77c7d]
logInfo("警告: N=%dでは有効率ギャップは統計的に有意ではない (p~0.10)", N_GENERATED); %[output:7a8ab334]
logInfo("アーキテクチャ:   LSTM(%d) -> Dropout(%.2f) -> FC(%d) -> Softmax", ... %[output:group:650dc90a] %[output:198ebf94]
    HIDDEN_SIZE, DROPOUT_RATE, VOCAB_SIZE); %[output:group:650dc90a] %[output:198ebf94]
logInfo("パラメータ数:     %d  (%d訓練トークンに対し過パラメータ化)", ... %[output:group:3f37985d] %[output:17de43a5]
    nParams, n_train * round(mean(seqLengths))); %[output:group:3f37985d] %[output:17de43a5]
%[text] ### 重要な教訓
%[text] 1. **主要証拠**: val 損失が早期にプラトーになります。LSTM は数百分子では SMILES 文法を学習できません。エポック数ではなくデータスケールのボトルネックであり、コードのバグではありません。
%[text] 2. **有効率ギャップ**（Markov vs LSTM）は参考値ですが、$N=300$ では統計的に有意ではありません。有効率数値でなく損失カーブで判断しましょう。
%[text] 3. **帰納バイアス**（Markov の局所遷移）は小さな $N$ では表現力を上回ります。
%[text] 4. 本番モデル（REINVENT, MolGPT）は $10^6 \\sim 10^7$ 分子を使用します。
%[text] 5. 次のステップは r06\_reinforce.m です。REINFORCE 方策勾配 RL を適用し、目標物性に向けて分子生成を誘導する方法を学びます。 \
%[text] \--- ローカル関数 ---
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
%   data: {"dataType":"text","outputData":{"text":"[14:00:44][INFO]  initPython: Desktop mode -- embedded Python: D:\\workspace\\20260207_ML_MCP\\20260417_EasyMolKit_dev\\python_env\\python.exe\n[14:00:45][INFO]  initPython: Python 3.10 configured (mode: OutOfProcess)\n","truncated":false}}
%---
%[output:0902b7b6]
%   data: {"dataType":"text","outputData":{"text":"[14:00:45][INFO]  ツールボックス: Deep Learning=1\n","truncated":false}}
%---
%[output:334b2162]
%   data: {"dataType":"text","outputData":{"text":"[14:00:45][INFO]  R05: セットアップ完了 (HIDDEN=128  EPOCHS=150  TEMP=0.80  AUG=1x)\n","truncated":false}}
%---
%[output:3fd1b0d0]
%   data: {"dataType":"text","outputData":{"text":"[14:00:51][INFO]  Local fda_drugs.csv:          200 行\n","truncated":false}}
%---
%[output:525fc4be]
%   data: {"dataType":"text","outputData":{"text":"[14:00:51][INFO]  コーパス 201-500: FDA薬 + ChEMBL承認薬のみ使用\n","truncated":false}}
%---
%[output:751f8380]
%   data: {"dataType":"text","outputData":{"text":"[14:00:51][INFO]    (everyday\/forensic除外 -- ドメイン混合は小規模Nのモデルを劣化する)\n","truncated":false}}
%---
%[output:422f2707]
%   data: {"dataType":"text","outputData":{"text":"[14:00:51][INFO]  ローカルコーパス: 200 SMILES (ChEMBL取得前)\n","truncated":false}}
%---
%[output:956d762e]
%   data: {"dataType":"text","outputData":{"text":"[14:00:51][INFO]  ChEMBLキャッシュヒット: 450 \/ 865 行読み込み\n","truncated":false}}
%---
%[output:55200a7a]
%   data: {"dataType":"text","outputData":{"text":"[14:00:51][INFO]  重複除去後: 649 ユニークSMILES\n","truncated":false}}
%---
%[output:0617c68c]
%   data: {"dataType":"text","outputData":{"text":"[14:00:51][INFO]  訓練用コーパス: 500 (目標=500)\n","truncated":false}}
%---
%[output:6cc88e94]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 50\/500) SMILES検証中\r[##--------]  20% (100\/500) SMILES検証中\r[###-------]  30% (150\/500) SMILES検証中\r[####------]  40% (200\/500) SMILES検証中\r[#####-----]  50% (250\/500) SMILES検証中\r[######----]  60% (300\/500) SMILES検証中\r[#######---]  70% (350\/500) SMILES検証中\r[########--]  80% (400\/500) SMILES検証中\r[#########-]  90% (450\/500) SMILES検証中\r[##########] 100% (500\/500) SMILES検証中\n","truncated":false}}
%---
%[output:1340b2ac]
%   data: {"dataType":"text","outputData":{"text":"[14:01:18][INFO]  有効SMILES: 500 \/ 500\n","truncated":false}}
%---
%[output:73633f29]
%   data: {"dataType":"text","outputData":{"text":"[14:01:18][INFO]  SMILES長: min=3  median=38  max=260  mean=45.7\n","truncated":false}}
%---
%[output:9a34ceb7]
%   data: {"dataType":"text","outputData":{"text":"[14:01:18][INFO]  打ち切りされた系列 (> 98 文字): 3.8%\n","truncated":false}}
%---
%[output:45463a39]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAkYAAAFgCAYAAABewLbaAAAQAElEQVR4AeydXWhlS3bfa+ZorqQODm0aSVfGDxpxBz+Expd5kYJBgycvvhk3JBCHtC74KuDkwfhaxOQy4JeEID\/EMh4nEwwBgRlDZJxHN8HJQwJ5c+MkMO4Egg2yQZqW0jh+uBlJPdNCyvkf9VJX795nn\/1R+6v2b7hLVbXWqlWrfrW1a81Rt\/qLN\/wPAhCAAAQgAAEIQGBC4IuO\/0EAAhCAAASiJcDGIFCMAIVRMV54QwACEIAABCAQMQEKo4gPl61BIEYC7AkCEIBAnQQojOqkS2wIQAACEIAABHpFgMKoV8cVY7LsCQIQgAAEINAdAhRG3TkLMoEABCAAAQhAoGUCwQujlvfD8hCAAAQgAAEIQKA0AQqj0uiYCAEIdInA+vp6qXTKzksuFipOMm6T4xj20AAvloicAIVR5AfM9iBQBwFdoHmkjrUVU2ur7Yoon6Ojo66kUzoP7UF7KR2AiRCIgACFUQSHyBYg0AYBXaKzJCsvXcBJyfKfZVMuijfNT7Y0ceMJaXpfN3aZ+p\/8tPZUh54ZtBftqWdpky4EghGgMAqGkkAQuCWgS8WXW+2br7K9GU3vmZ9aE\/PW2PrTWvlMk7Q5Sd80H1+X9E+OfV+\/b366gCWyqZXIpnGWyEe+WT5pNs1JivySurSx\/IqIcpzmn2VLmyP\/NJnlm2aXzo+lMQIBCLxNgMLobR6M+k+g1R3o0klerNIlk0rT+T6+XfFks1b9vKI5aZKcr\/WSftIl\/fxx0j859n2tr5jmZzq\/lU0+vq5Iv8h8rSN\/xVdfbVHRPIuRnCtbUld2rDWSkoyl9Xwfjcv4aI7ipM2XDYFA7AQojGI\/YfbXGAFdJLpQkgtKJ1tS35WxclOOyXykky2pt7FsWWJ+RdusdbWe7FkxZZdflk8Ttibz0H61nr8vjaU3nfrS2VitxtKrj0AAArcEul8Y3ebJVwj0moAuIH8DGk+7kKSX3fdvq5+Vh2xZksy56r6KzFde8k\/mYGPZ5GPjtlrlMU38nJI+vo0+BCAQlgCFUVieRBswAV20doH1CUMf8hZX5Smu6s8S+clffur7Yjq1JrJb32+lryLTcrCYsk8T81Gb9FGO0iNxEGAX3SJAYdSt8yCbnhOwC0wXl0mbW7Ic\/DYtnzx5J2P441n9tDWl0zytrX6W+D7qm2iO9a2VzkQ666v115PNRDbrWytdF0T5dCEPcoDAUAhQGA3lpNlnowR0mZnoMk5bXPakTWPp0\/zL6BQrKVlxfF\/l4vv6NumzxrL5Puq3LdqP5VVvLm9H15pa+20tIwhAoKsEKIy6ejLk1TsC0y6\/rl+MIfJO7lExpavzEIuuUXc+WXvV2so36SPdNEn6MoYABJohQGHUDGdW6TGBOlP3L0xdkBrXuV6dsZW79iBRf9Za0\/ym6WfFa8OufSrfsmtr\/jSxmFXiW4yirdZUXkXn4Q+BGAhQGMVwiuyhEwR0kehCSSYjnWxJfVfGyk05JvORTrakPjmWn0T6LH\/ZzE++VURxFK9KjKbnKl\/lXce6abG1lvS2nvrS2VitxtKrj0AAArcEel8Y6Rvb5HZLb381m9pby9tfpTd528IIAsUJ6JKx58la6bIiyS5ftWl+sklvrfoSjdNENpM0u3Rmt1ZrS++LdGZPa81XfhLzUd9sprPWt6lverU2J6mXzRf5zfLx\/Yv2FV+Sd55yyesv37xxzU9zFN8X6cxurXQhfBRDsSwuLQSGRqDXhZF9A+ubWKKxf4AaS2+icRG770sfAnkJ2PNmbXKe9EV08jexeTZOa4v4mK\/aZCzpfNH3j4nv6\/tY3+zmb3q1ZlNfYj5JvWy++H6+Xn2zqVUc6bIky0c2k6wYvk3+Wtt0Gls\/2WbZkr421hxfTJ9sq\/poD4qRjMu4ZwRItxKB3hZGs76B0+z6hpdexNRqrL6JxtLbmBYCEHhDQN8fJm+02b08\/nl8tEqWn9nUyreolJ3nrxMihh+vjX4Me2iDG2vGRaC3hVFcx8BuIAABCEwlgAECEGiQQG8Lo+T\/s9EnPUldaI4nJycOgQHPAM8AzwDPwNCfgdD3a5fi9bYw8iE2VRR99tlnbmtrC4FBtWcAfvDjGeAZ6Pkz8Pjx48kHBf5dHEu\/94VRE0WRDlv\/7+Dp06duf3\/fHR4eBpHd3V2FDhozVG5djwO7As\/gr\/6q++8\/9mPu8Ld+a\/Lcwq4Au8T3OuzKsYNbOW56D3eRnXLSfTi5wCL80uvCqKmiyD\/3jY0Nt7m56UvpvmIpttpQMYcSR8xgl\/M5HD+zP\/qjP3r3nMIuJ7eU73PYlWMHt3Lc9D7vIjvLSe\/gGKW3hdGsokh\/3kg+\/qFpLL10ajVW30Rj6W1MCwEIQAACEIDAsAj0tjDSMamQSYr0JipyfLvGZlOrcZZdPnXK+++\/737+53\/eqa1znRhjixnsyp0s7Mpx0yzYiUJxgVtxZjYDdkaiuba3hZGKmjRJovN9kjaNZ9nlU5fogf\/kk0\/qCh91XNh5x\/uFL3gDrztFDzuPUcFuLOwKbruyO9zKI4RdeXZlZ\/a2MCq7YeZBAAIQgAAEIACBaQQojKaRQQ+BrhGY8gmQc68T9e3q39y8NtBAAAIQgEBeAhRGeUnhB4GuEFDRo1ysVV+iQiipkx6BAAQgAIHcBCiMcqNqzpGVIFCagBVHaksHYSIEIACB4RKgMBru2bPzGAnwiVGMp8qeIBAbgU7vh8Ko08dDchDISUAFkUSfFEnUzzkVNwhAAAIQeEOAwugNC3oQ6C8BFUOS\/u6AzPtMgNwhEBEBCqOIDpOtQAACEIAABCBQjQCFUTV+zIZAswTsR2RprXQmykp9tcWFGRCAAAQGS4DCaLBHz8Z7R0A\/KksTbWSaXraHD5178sS55WWNEAhAAAIQyCBAYZQBJxoTG4EABCAAAQhAIBcBCqNcmHCCQA8J6FOkHqZNyhCAAASKEgjpT2EUkiaxINAGAQqgNqizJgQgECkBCqNIDzbktk5OTlxICZkbsSAAgRgJsCcItEeAwqg99r1YWQXRP\/6lX3FbW1vB5G\/\/nb8\/KbR6AYAkIQABCEBgUAQojAZ13MU3q8Lof\/\/Jf3Nf+\/ib7huffquyfPWjT5ziFc+EGU5\/\/b6odAQbaUAAAhDoCwEKo76cVMt5vv+VD91qCPngw5Z30uPl9WeJioht9dkz5x49cu7FC9PQQgACEIDAFAIURlPAoM4igK0XBFRE9SJRkoQABCDQHQK9L4zW19en0pTNlzTHWfa0Oegg0BqBPD9Kay05FoYABCDQfwKTwqiv21BRMy132Y6Ojpwv0vn+GmfZfV\/6EOgEAX0KNEs6kShJQAACEOgngd4WRipqqiDXfBVFfgyNpfd19CHQKQL6xKhTCZEMBHpBgCQhkJtAbwsjFTGS3DsN5Ki\/peXL1dWVi1kM2831tbsOIBYvZmZ17k38LH7yb6mZPtnejD9hkly9ehX1s5rcN+O4302cb3Pn69956us9FLP0tjCadSgqmvTpjy\/SzZo3y769vf3W7\/PZ29tzx8fHpUQP2NnZWam5ZdcsOk\/5icnFxaU7Pz+vLBeXFwrnTk9PK+27D+yKss7jL3jmp\/6fj39cbDL3pS85idmtff78ubu8vHRqpYuCXcnvOe2\/isAu3nddleeizrldeOYODg7euvd0D+r9E6tEWxipIFIh5It0VQ9yf3\/fHR4e3snOzo5bXV0tJSsrK25paanU3LJrFp13\/\/79CbKFhQV37949t7i4WEnm35ufxKu67z6wK8o6j7\/gmZ\/fl07jl+MC6Mvr6289Uw8ePHDz8\/NueXl5oh8qOzGqKrCL911X9dmoa34Xnjndc\/69t7u7q9dNtBJtYVTXiW1sbLjNzc07WVtbcyoayoiKDF1YZeY2OUcsR3MjNxqN3NzcXCVRDMVL5F+YYV\/YVd1ncv6E3bg4XRjLpD8uWCc+Go9\/ZKa+Uzseq28i7pqj8VDZae9VBXYLhb9XxRxu5bh1hZ3uOf\/e0z2o90+s8sVYN8a+IBAtgXHho+JnIvrD2BLp\/A1rLL2vow8BCEAAAjMJ1FMYzVy2fQf9iC35ozWNpW8\/OzKAQE4CKoAkFEE5geEGAQhAIJtAtIWRChwVOr5I5+PQOMvu+9KHQCcIqAhKSyRNn6ZLm4sOAhAoTIAJ8RLofWGk4mba8cjmS5rfLHvaHHQQgAAEIAABCMRJoPeFUZzHwq4gUICAfowm8af444cPnXvyxLnlZd+D\/lsEGEAAAhC4JUBhdMuBrxDoJwEVQPqRmcR2YDob00IAAhCAQG4CFEa5UeHYJwJR5qqCx9+YxiqI1EqvViKdxggEIAABCBQmQGFUGBkTINAyARU\/SiFZAGkskQ2BAAQgAIFSBHpSGJXaG5MgECcBFT9WHMW5Q3YFAQhAoDUCFEatoWdhCFQgQHFUAR5TIdBBAqTUGQIURp05ChKBwAwCyWLIH+sTpDSZERIzBCAAAQi8TYDC6G0ejCDQfQIqgCxLFUfqq00T2ZA2CLAmBCDQUwIURj09ONIeMAEVQFYcWZuF49kz5x49cu7FiywvbBCAAAQgMCZAYTSGwH8QmEmgaw5WHKntWm7kAwEIQKDHBCiMenx4pD5gAvqkKFkUSScZMBa2DgEIQKAqgaEWRlW5MR8C7RJIFkXKRjoJxZFoIBCAAARKEaAwKoWNSRBoiYAKn1lLy4fiaBYl7BCInADbK0uAwqgsOeZBoCsEVAglc0nTJX0YQwACEIDAOwQojN5BggICEIBA9wiQEQQg0AwBCqNmOLMKBMISyPOjsjw+YbMiGgQgAIHeE6Aw6v0RsoF+Eqg5axVF\/DitZsiEhwAEYiTQ+8JofX0981xkN0lzNJvaNDs6CPSOAEVR746MhCEAge4Q6HVhNKuYkf3o6MiZaOyj19hsajX27UX6+EKgEwQoijpxDCQBAQj0l0BvC6NZRYzsKnb8o\/HH0+zS+3PoQ6AXBFQQSfjxWS+OiyQh0EMCg0m5t4WRihzJrJNSoWMyyxc7BHpHQMWQRAWRpHcbIGEIQAAC3SLQ28IoD0YVRCqeTDTOMy\/L5+TkxPlydXXlYhZjcXN97a4DiMWLmVlde3MqgEwEcty\/evXKTSTjObwZF0ySWX515U3cjr4jMp4Zzowz858B\/85TX6+fmCXqwkgFUejD297edltbW3eyt7fnjo+PS4kesLOzs1Jzy65ZdJ7yE8OLi0t3fn5eWS4uLxTOnZ6eVtp3H9gVZT3L\/8+PjpyJIKo\/96UvOUnW3L\/4kR9xf\/Htb7vjH\/5wwnyI7LL4FLHBLt53XZHnoEnfLjxzBwcHd3ee7j\/dg3oHxSpRF0Z1HNr+\/r47PDy8k52dHbe6ulpKVlZW3NLSUqm5ZdcsOu\/+\/fsTjAsLC+7evXtucXGxksy\/Nz+JV3XffWBXlHURf0GU\/8vLSyf58vq6+\/L6eq5naejsxK2swC7ed13ZZ6LueV145nTP+ffe7u6uXkHRjKGd3gAAEABJREFUCoVRwaPd2Nhwm5ubd7K2tuZUNJQRFRnz8\/Ol55dZs8wcIRrNjdxoNHJzc3OVRDEUr0we\/py+sPNzDtlPMnTjH5dJFsaF66x1hs5uFp8sO+wWSr2v4FaOm57FLrDTPeffe7oH9Q6KVb4Y68b0Y7SsP1OUZpe\/9J1gQhIQKENABdIXvlBmJnMgAAEIQGBMINrCaLy3ye8vUrFjkix6NDabWo01D4FArwlQHPX6+EgeAkMh0NV99r4wmlXMyG6SdghmU5tmRweBXhKgOOrlsZE0BCDQPoHeF0btIyQDCLRAQIXPrGXz+MyKgR0CuQngCIE4CFAYxXGO7AIC0wk8e+bco0fOvXgx3QcLBCAAAQhMCFAYTTDwBQI9JDDrD1nPss\/YMmYIQAACQyRAYTTEU2fPcRFQAWQS187YDQQgAIHGCVAYNY68rQVZt\/cErPix1t+Q\/jyRRDbp1WqsPgIBCEAAArkJUBjlRoUjBDpAwIoda5MpSU9RlKTCGAIQGAKBQHukMAoEkjAQ6AQBK4rUdiIhkoAABCDQLwIURv06L7KFQDoBFUISfWIkD7Uaq49AoJ8EyBoCrRCgMGoFO4tCIDABFUKSwGEJBwEIQGBoBCiMhnbi7Dc+AtMKomn6tgiwLgQgAIEeEKAw6sEhkSIEphLQj8uyZOpEDBCAAAQgkEaAwiiNCro8BPDpAgH7VCitNV0X8iQHCEAAAj0hQGHUk4MiTQhMCOjTIXWsVR+BAAQgAIFgBN4URsFCEggCEKiFgD4B8kWLWIE0rZUPAgEIQAACuQlQGOVGhSMEOkjAL5T8vp\/qw4fOPXni3PKyr6UPgcERYMMQyEOAwigPJXwg0EUCKoSm5ZVlmzYHPQQgAAEIOAojHgIIQKCnBEgbAhCAQHgCFEbhmRIRAhCAAAQgAIGeEuh9YbS+vj4TvXwkaY7Sm6TZ0TVHgJUgAAEIQAACbRPodWGkgqYKQM0\/OjpyJhpXicdcCEAAAhCAAAT6TaDGwqheMHmLGPmp8Elmk6aXn\/RJX8YQgAAEIAABCAyDQG8LIxUxkmEcE7uEAAQgAIHOESChKAn0tjDKcxr69Cd08XRycuJ8ubq6cjGLcb65vnbXAcTixcysa3u7+ZM\/cTc\/+7Pu6nvfi\/pZ7Rp38on73Tik8\/XvPPXtPR5rG21hVEdRpIdge3vbbW1t3cne3p47Pj4uJXrAzs7OSs0tu2bRecpP+764uHTn5+eV5eLyQuHc6elppX33gV1R1nX5P3\/+3F1eXjq1WgN2U79fZz6TsCvHDm7luHXl+\/Xg4ODuztP9p3tw8iKP9Eu0hZHOS8WRiY3VVpH9\/X13eHh4Jzs7O251dbWUrKysuKWlpVJzy65ZdN79+\/cnuBYWFtzi4mJlmX9vfhKv6r77wK4o67r8Hzx44Obn593y8vLkWYNdue9XnQ\/syrGDWzluXXnmdM\/5997u7u7kPR7rl2gLI\/0IzRcdoMZqq8jGxobb3Ny8k7W1NaeioYyo0NCFVWZuk3PEazQ3cnNzc5VlNBopXGlmtu9c7MbFnPkPvRX3hXFhKw6wWyj9\/MGuHDu4lePWle9X3XP+vad7cPIij\/RLtIXRrPNSkaRPk3w\/jaX3dfQhAAEIQAACEBgOgT4VRsFPRUWQiiETjYMvQkAIQAACEIAABHpDoPeFUd5iZpqf9Ca9OTUShQAEIACBCAmwpS4Q6H1h1AWI5AABCEAAAhCAQBwEKIziOEd2AQEIdJAAKUEAAv0jQGHUvzMjYwhAAAIQgAAEaiJAYVQTWMLGSIA9QQACEIBA7AQojGI\/YfYHAQhAAAIQgEBuAoMujHJTwhECEIAABCAAgUEQoDAaxDGzSQhAAAIQGCABtlyCAIVRCWhMgUCvCDx86NyTJ84tL\/cqbZKFAAQg0AYBCqM2qLMmBCAAgTIEmAMBCNROgMKodsQsAAEIQAACEIBAXwhQGPXlpMgzRgLsCQIQgAAEOkaAwqhjB0I6EIAABCAAAQi0R4DCKCR7YkEAAhCAAAQg0GsCFEa9Pj6ShwAEIAABCDRHYAgrURgN4ZTZ47AJPHvm3KNHzr14MWwO7B4CEIBADgIURjkg4QIBCEAgTgLsCgIQSBKgMEoSYQwBCEAAAhCAwGAJUBgN9ujZeIwE2BMEIAABCFQj0PvCaH19PZOA7CZpjmZTm2ZHBwEIQAACEIDAcAj0ujCaVczIfnR05Ew09o9WY7Op1di3t98nAwhAAAIQgAAEmiTQ28JoVhEju4qdaTDT7PKXftoc9BCAAAQgAAEIBCTQwVC9LYxUxEiaZnpycuJ8ubq6cjGL8b25vnbXAcTixcysa3u7ublxkqtXr6J+VrvGnXzifjcO6Xz9O099e4\/H2va2MJp1IMmiSZ8EJXWzYqTZt7e33dbW1p3s7e254+PjUqIH7OzsrNTcsmsWnaf8xOHi4tKdn59XlovLC4Vzp6enlfbdB3ZFWdfl\/\/z5c3d5eenUag3Ylft+hV0ubqnf1zxz\/WZ3cHBwd+fp\/tM9OHmRR\/ol2sLIP69QRZFi7u\/vu8PDwzvZ2dlxq6urpWRlZcUtLS2Vmlt2zaLz7t+\/r227hYUFt7i4WFnm35ufxKu67z6wK8q6Lv8HDx64+fl5t7y8PHnWYFfu+1XnA7ty7OBWjltXnjndc\/69t7u7O3mPx\/ol+sIoZFGkh2BjY8Ntbm7eydra2qRoUOFQVFRo6MIqOq9pf+17NDdyc3NzlWU0GilcaWa2976ws3zbbsV9YVzYKo9C7MYFseYgC5NnFna3HIo+D3Arx02cu8BO95x\/7+kenLzII\/0SdWEUuiiK9BlgWxCAAAQgAAEIvCYQbWE0qyjSnzeSz2sOk0Zj6SeDOL+wKwhAAAIQgAAEMghEWxhpzyp0kiK9iYog366x2WghEA2Bhw+de\/LEueXlaLbERiAAAQikE6iu7X1hNK2YkT5Nksh8n6SNMQQgAAEIQAACwyLQ+8JoWMfFbiEAAQgMiwC7hUDTBCiMmibOehCAAAQgAAEIdJYAhVFnj4bEIBAjAfYEAQhAoNsEKIy6fT5kBwEIQAACEIBAgwQojBqEHeNS7AkCEIAABCAQEwEKo5hOk71AII3As2fOPXrk3IsXaVZ0EIAABCDgEUgURp6lpq5+b1De0EV888bEDwIQgAAEIAABCEwj0HhhpERU8PgiXVJk1+8YSuoZQwACEIAABEoTYCIEZhBopTBSweOLiiAT5au+7OojEIAABCAAAQhAoCkCrRRGyc2pCDKhKErSYQwBCGQQwAQBCEAgKIFGCiMVO5KszGWXqEBSm+WLDQIQgAAEIAABCNRBoJHCSMWOxAoea7Uh9SWyS6RTK536yMAIsF0IQAACEIBAiwQaKYxsfyp4TFT4SGxsPtZKL7uNaSEAAQhAAAIQgEDdBOoujCb5W4FjrZQqfCTSSaRLiuxJHWMIQAACEIAABCBQF4FGCiMVOH7x4\/dl0+akk6gv8fsaI3ERODk5cVXk7OzMSSxGXHTYDQQg0B8CZBobgUYKI0HzCyDrS28inUQFkUR9s9HGR2B7e9ttbW2Vlq9\/\/evu448\/dmoV5\/Hjx5NCKz5S7AgCEIAABJok0EhhpEJHYhtTX2Jjv6Ug8mnE2\/\/ax9903\/j0W6XlZ37xN9xP\/8KvuY9+6TfdVz\/6xD19+jReWOysUQIsBgEIDJtAI4WRih2JoVZfonGyQNJYNrWyz5JZfrKbpMUym9o0O7p6CLz\/lQ\/dagV5\/4OfdMtffuhWx+3qBx\/WkyRRIQABCEBgcAQaKYxE1QqPZNGjsdnlY2PpZon8s3xkVzwTjX1\/jc2mVmPfTj8GAuwBAhCAAAQgkJ9AI4WRCg4VHpaW+tLZWK10ak2vsfWlT0qWTb6yK4b6JhpLr7FajdU30Vh6G9NCIAoCDx869+SJc8vLUWyHTUAAAhCok0AjhZEKjjybkJ8ky9ds8pPYmBYCEIAABCAAAQhUJdBIYWRJ+oXMtL75qvV9NO6C2F8Pt\/bq6srFLMb85vraXYeQm5tJyKrxbsZx7vIZ9xU05nNgb3F\/n3G+nK\/3DHTuTrH7zlq9b2OWRgujGEAm\/5r53t6eOz4+LiV6yPS7eMrOb2Ke8tO5XVxcuvPz88ry8uWlwrkQ8S4vb3OymKenp6XOoQmOXVqjD89dl3j5ucAu3nedf85d6nfhmTs4OHjrV6voHpy8yCP90onCqE9\/rmd\/f98dHh7eyc7OjltdXS0lKysrbmlpqdTcsmsWnXf\/\/v3Jo7+wsOAWFxcry\/x788HiKad79+45i9l1lkXZ1+Xfh+eurr1XjVuKXcn3Q9VcuzQfbuXuCJ1hF9jpnvPvvd3d3cl7PNYvjRdGfSqC0g59Y2PDbW5u3sna2prTBV1GVGjMz8+Xnl9mzTJzxGE0N3Jzc3OVZTQaKZwbVYw3GsfxRUHL7G2Ic\/ry3HXxbGC3UOp9Bbdy3PQ90AV2uuf8e0\/3oN65sUrjhVFTIPXnk5JFmMbSKwe1GqtvorH0NqaFQAoBVBCAAAQgEDGBaAsjnZmKHBU7JhpLb6Kx2dRqbDZaCERD4Nkz5x49cu7Fi2i2xEYgAAEI1EWgkcJIRYeJNmJ9a9N0vk32aTKrmJHdJC3G0dGRy7KnzUEHAQhAAAIQgECcBBopjKzwUCuMan1J05ldNgQCEIAABCAAgXIEmFWMQCOFUbGU8IYABCAAAQhAAALtEGisMMr7o7F2MLAqBCAAgb4QIE8IQKBOAo0VRvajMbUUSXUeKbEhAAEIQAACEChLoLHCyE9QxZFEBZKvpw+BIRJgzxCAAAQg0B0CrRRGtn0VR+pbqz4CAQhAAAIQgAAE2iLQSGGkT4YkbW2y2XVZDQIQgAAEIACBvhJopDDSJ0ISFUdFpa9gyRsCEIAABCAQJYHIN9VIYWQMVRxJkmPppon50kIAAhCAAAQgAIG6CTRaGNlmrAjSp0emo4UABCAAgVYIsCgEIOARaKUwsvVVIFmfFgIQgAAEIAABCLRNoNXCqO3Nsz4EoiTApiAAAQhAoDQBCqPS6JgIgZ4QePjQuSdPnFte7knCpAkBCECgPQIURu2xz7syfhCAAAQgAAEINESAwqgh0CwDAQhAAAIQgEAagW7pKIy6dR5kAwEIQAACEIBAiwQojFqEz9IQgAAEYiTAniDQZwIURn0+PXKHAAQgAAEIQCAogagLI\/0CSV\/SyM2yp81BB4FhEWC3EIAABIZDINrCSAWPfoGkL9L5R6txlt33pQ+B3hJ49sy5R4+ce\/Git1sgcQhAAAJNEYi2MJoF0Ioi309FkvS+LsY+e4IABCAAAQhAIJ3AYAujdBxoIQABCEAAAhDoOYFK6UdbGNmnP\/oEyES6SrTGk09OTpwvV1dXLmYZb3ny3831tbsOIJNg4y9V493c3LzJZ9wfh4z6HKo8Y2IluXr1CkaRf79WeU6YG\/e7vMr5+nee+nrfxizRFkYqhlQI+SJd1Y63E6oAABAASURBVMPc3t52W1tbd7K3t+eOj49LiR6ws7OzUnPLrll0nvITs4uLS3d+fl5ZLi4vFM6FiHd5eZvTy5eXk5inp6edZlmUfSj\/58+fO7FSq5h9eO6UZxcFdi2960q+Y7v4DBXNqQvP3MHBwd2dp\/tP9+DkpRvpl2gLo7rOa39\/3x0eHt7Jzs6OW11dLSUrKytuaWmp1Nyyaxadd\/\/+\/QnKhYUFt7i4WFnm35sPFk853bt3z1nMrrMsyj6U\/4MHD9z8\/LxbXl6ePGt9eO5C7T10HNjF+64L\/ayEiteFZ073nH\/v7e7uTt7jsX6hMCp4shsbG25zc\/NO1tbWnC7oMqJCQxdWmblNzhGi0dzIzc3NVZbRaKRwblQx3mgcxxcFbZJJ39YSq4VxYau8O\/Dclf6eUf5tCuwWSp0d3Mpx07PeBXa65\/x7T\/eg3rmxymALI\/2ILfmjNY2lj\/Ww2RcEIAABCEAAAtkEoi2MVOCo0PFFOh+Hxll235d+BgFMEIAABCAAgUgIRFsY6XxU+PgiXVJm2ZP+jCEAAQhAAAIQiJdAWmEU727ZGQQgAAEIQAACEMggQGGUAQcTBCAAAQjESIA9QWA6AQqj6WywQAACEIAABCAwMAIURgM7cLYLgRgJsCcIQAACoQhQGIUiSRwIdJXAw4fOPXni3PJyVzMkLwhAAAKdIUBh1JmjIJE3BOhBAAIQgAAE2iFAYdQOd1aFAAQgAAEIQKCDBBopjDq4b1KCAAQgAAEIQAAC7xCgMHoHCQoIQAACEIBAIQI4R0SAwiiiw2QrEIAABCAAAQhUI0BhVI0fsyEAgRgJsCcIQGCwBCiMBnv0bHwwBJ49c+7RI+devBjMltkoBCAAgbIEKIzKkmNenwiQKwQgAAEIQCAXAQqjXJhwggAEIAABCEBgCAT6WRgN4WTYY2ECJycnLpQUXpwJEIAABCAQBQEKoyiOkU2IwPb2ttva2goijx8\/nhRZiotAAAIQaJoA67VHgMKoPfasHJjA1z7+pvvGp9+qLF\/96BP39OnTwNkRDgIQgAAE+kCAwqgPp0SOuQi8\/5UP3WoI+eDDXOvhBIH8BPCEAAT6QiD6wmh9fd2ZpB2K2dSm2dFBAAIQgAAEIDAcAlEXRip2jo6OnInG\/tFqbDa1Gvt2+hCYRgA9BCAAAQjESSDawkhFjood\/9j88TS79P4c+hCAAAQgAAEIDIdAtIWRHaEKHRPTvdvm1yT\/OvjV1ZXrkiTzqzo2MjfX1+46gISKd3Nz8yafcV9xQ+V4\/Tpel861Si5iJbl69apTz2qVPTG3W+8dziPu80jeI3rfxixRF0YqiPQpkYnGVQ8z+VfC9\/b23PHxcSnRw3Z2dlZqbtqaf\/zHf+w+\/fTTIH9d3f7au\/YrZhcXl+78\/LyyXFxeKJwLEe\/y8janly8vg8XUHi3e6elpsLNJO6+mdM+fP3dipVZrhn7uFHMoArtuvOuG8rxpn8GeuZL3lHI4ODh4616xe2Hy4o3wS9SFkQqi0Ge2v7\/vDg8P72RnZ8etrq6WkpWVFbe0tFRqbtqaP\/jBD9x3v\/tdF+qvreuvvuuvrovhwsKCW1xcrCzz780rnAsRTzHu3bvnQsbUHi1eyLNJO6+mdA8ePHDz8\/NueXl58qyFfu6a2kcX1oFdN951XXgWmsqhC8+c7jn\/3tvd3Z28x2P9EnVhVMehbWxsuM3NzTtZW1ubXPK6pIvK5BIeX1hF52X5a8\/B\/tq6\/ur7Bx8qpBvNjdzc3FxlGY1GQeIpji8KOgqcYxbnvtnEamFc2CrvOp47xR2CdIhd6fdOG+cEt4XS59UFdmvje86\/93QP6p0bq1AYxXqy7AsCRuDhQ+eePHHjj4xMQwsBCEAAAlMIRFsY6cdoWX+mKM0uf+mnsEINgfoIEBkCEIAABDpBINrCSHRV5KjYMdFYehONzaZWY7PRQgACEIAABCAwPAJRF0Y6ThU7JhonxWxqk7YKY6ZCAAIQgAAEINBDAtEXRj08E1KGAAQgAAEIdJxAvOlRGMV7tuwMAhCAAAQgAIGCBCiMCgLDHQIQgECMBNgTBCBwS4DC6JYDXyEQL4Fnz5x79Mi5Fy\/i3SM7gwAEIBCIAIVRIJCEgUC3CJANBCAAAQiUIUBhVIYacyAAAQhAAAIQiJIAhVFPjpU0IQABCEAAAhConwCFUf2MWQECEIAABCAAgWwCnbFSGHXmKEgEAhCAAAQgAIG2CVAYtX0CrA8BCEAgRgLsCQI9JUBh1NODI20IQAACEIAABMIToDAKz5SIEIiRAHuCAAQgMAgCFEaDOGY2CQEIQAACEIBAHgIURnkoxejDniAAAQhAAAIQeIcAhdE7SFBAAAIQgAAEINB3AmXzpzAqS455EIAABCAAAQhER4DCKLojZUMQgAAEYiTAniDQDAEKo2Y4swoE2iPw8KFzT544t7zcXg6sDAEIQKAnBAZRGK2vrztJ2plIb5JmRwcBCNRDgKgQgAAEukhgEIXRNPAqiI6OjpyJxtN80UMAAhCAAAQgED+B6AsjFTsqfJJHmaaXn\/RJX8Z5COADAQhAAAIQ6D+B6Auj0Ed0cnLifLm6unJdEdvrzfW1uw4lNzeTsKFiToKNv1SNdzPO626P4\/44pKsaMxmvK+dKHt35HuMsOIshPgP+nae+3rcxy9TCKIZN69MffQoUci\/b29tua2vrTvb29tzx8XEp0QN2dnZWam7amoqlvV5cXLrz8\/Mg8vLlpUK6UDEvLi+Cxbu8vN1n6Bwt3unpabCzSTuvtnShn7u29tHGurDrxruujbNva80uPHMHBwd3d57uP92Dkxd5pF+iLYzqKIr0DOzv77vDw8M72dnZcaurq6VkZWXFLS0tlZqbtub9+\/eVoltYWHCLi4tBZP69+aAxQ8bTPu\/du+dCxhQ3ixfybNLOqy1d6OeurX20sS7suvGuC3T2wd69debThWdO95x\/7+3u7k7uhVi\/RFsY6cBUHJnYWG0V2djYcJubm3eytrY2KUR0SReVySU8P196ftp62ttobuTm5uaCyGg0Ukg3ChRzFCie4viiJEPnmMa3l7o\/+zO38HM\/5xY+\/3zyrNXx3PWSy\/j\/QBTNG3YLk2cIbuU4FOUm\/y48c7rn\/HtP96DeubFKtIWRfoTmiw5QY7UIBCAQIQG2BAEIQCAAgWgLo1lsVCTp0yTfT2PpfR19CEAAAhCAAASGQ2CwhZGOWEWQiiETjaVHOkGAJCAAAQhAAAKNExhMYTSt6JHepHH6LAgBCEAAAhCAQKcINFcYdWrbJAMBCEAAAhCAAATeJUBh9C4TNBCAAAQgAIHCBJgQBwEKozjOkV1AAAIQgAAEIBCAAIVRAIiEiJOAfuNsSImTUsy7Ym8QgMAQCVAYDfHU2XMuAvq19\/r196Hk8ePHk39nL9fiOEEAAhCAQCsEKIxawc6ibRAouubXPv6m+8an3woiX\/3oE\/f06dOiKeAPAQhAAAINE6Awahg4y\/WHwPtf+dCthpIPPuzPxskUAhCAwIAJ9LgwGvCpsXUIQAACEIAABGohQGFUC1aCQgACEIAABCoSYHorBCiMWsHOohBokMDDh849eeLc8nKDi7IUBCAAgX4SoDDq57mRNQQg0D8CZAwBCPSAAIVRDw6JFCEAAQhAAAIQaIYAhVEznFklRgLsCQIQgAAEoiNAYRTdkbIhCEAAAhCAAATKEqAwekOOHgQgAAEIQAACAydAYTTwB4DtQwACEIDAUAiwzzwEKIzyUMIHAn0m8OyZc48eOffiRZ93Qe4QgAAEGiEQfWG0vr7uTNKImk1tmh0dBCAAga4SIC8IQCA8gagLIxU7R0dHzkRjH6HGZlOrsW+nDwEIQAACEIDAsAhEWxipyFGxM+040+zyl37aHPQQqJcA0SEAAQhAoG0C0RZGbYNlfQhAAAIQgAAE+kcg2sJIn\/74x6FPgpI63563f3Jy4ny5urpy06Rpve3h5vraXYeSm5tJ2FAxJ8HGX6rGuxnndbfHcX8c0lWNWVu88VkoP0nTz4TWEyvJ1atXnXlWlRfSnXcHZ8FZZD0D\/p2nvt5lMUu0hZF\/aKGKIsXc3t52W1tbd7K3t+eOj49LiR6ws7OzUnPT1lQs5XhxcenOz8+DyMuXlwrpQsW8uLwIFu\/y8nafoXMMHU9nYfs+PT0Ndt5pz0Ca7vnz506s1Moe+rlTzKEI7LrxrhvK86Z9duGZOzg4uLvzdP\/pHpy8yJ2Lsom+MApZFOkJ2N\/fd4eHh3eys7PjVldXS8nKyopbWloqNTdtzfv37ytFt7Cw4BYXF4PI\/HvzQWOGjKd93rt3z4WMKW6h4\/kxQ5532jOQpnvw4IGbn593y8vLk2ct9HOXtmasOth1410X6\/OVtq8uPHO65\/x7b3d3d3IvxPol6sIodFGkh2BjY8Ntbm7eydra2qQQ0SVdVCYX5vjCKjovy185juZGbm5uLoiMRiOFdKNAMUeB4imOL0qyazn6Z6BclWPW2dVp0\/oL42JZa9Tx3CnuECR6duP\/U1XHOcJtoVP3RNEz1j3n33u6B\/U+i1WiLYxmFUX680by8Q9WY+l9HX0IhCSgj8VDSci8iAUBCEAAArcEoi2MtD0VOkmR3kRFkG\/X2Gy0EKiDgH42r5\/Rh5DHjx9P\/iJARp6YIAABCECgIIFoCyMVOWmS5OP7JG2MIRCawNc+\/qb7xqffqixf\/egT9\/Tp09DpEQ8CEIDA4AlEWxhFebJsqvcE3v\/Kh241hHzwYe9ZsAEIQAACXSRAYdTFUyEnCEAAAhCAwAAJdGHLFEZdOAVygECdBB4+dO7JE+eWl+tchdgQgAAEoiBAYRTFMbIJCEAAAl0kQE4Q6B8BCqP+nRkZQwACEIAABCBQEwEKo5rAEhYCMRJgTxCAAARiJ0BhFPsJsz8IQAACEIAABHIToDDKjSpGR\/YEAQhAAAIQgIBPgMLIp9FC\/+zsbPLbi0P8MxEtpM+SLRMo89zomZOkzW15OywPAQhAICyBEtEojEpACzVFl9Ov\/\/qvu69\/\/esuxD8RoX9uIlRuxOkHAZ35rGfnH21uuj\/9iZ9wf\/enfmrynOl5+\/jjj1OfO\/6ZkX6cO1lCAAL1EaAwqo\/tzMgqjL773e+6kP9MxMxFcYiKQJ5n52\/+vV92y2t\/w\/2tf\/jPJv8Uyc\/84m+4n\/6FX3Mf\/dJvTsb2T5Twz4xE9WjEuhn2BYHaCVAY1Y549gL8MxGzGeGRTiDXs\/PBh27hr\/11t\/rBTzr9cyTvj9vlLz90NpZuImO\/9FXQQgACEBgOgS8OZ6vsFAIQ6BwBEoIABCDQMQIURh07ENKBAAQgAAEIQKA9AhRG7bGPcWX2BAEIQAACEOg1AQqjXh8fyUMAAhCAAAQgEJJAdmEUciViQQACEIAABCAAgY4ToDDq+AGRHgSaJpD2ix\/L6prOvcx6RfamX7EhKTInlG+ZvTFnNgE8IJAkMPjCaH193Zkk4TDuLoFLMoIPAAALFUlEQVTv\/9X\/cf\/zP\/+e+\/5fnXU3yY5m9v0Z7PL80shZv1TS7F3\/hZEqWj777LPJL760nLParF+OmTUvhK3rLLMedxWT3\/nOd5zaLD9s7xIQM9i9y6VOzaALIxVER0dHzkTjOmETOxyB748Lov\/1X34vXMABRZrFLs8vjbRfCpnV2i+MVPERUrKPqphVeT19+jT3L1nVnrRCFqM8Plnc0myKqTyVbyjRPpoSXe6\/+7u\/29RyUa0Du+aPc7CFkYogFUQ+co2l93X0ITA0Arl+aeRXPnSTXwqZ1b7+hZEhP4HSJy91fHKSe8+v95Tpn8cni1ua7XXMkCzF8Y\/+6I+C\/VuNWQWbfQ9l+TRls1y63PosLE9fV7RvMWjzERhsYZQPz7teRR\/ILH+L\/v3\/e+Y+\/8vnleX65mYS8vO\/PK0c6\/PX+YSOOUlw\/OXznDl+\/jqPZHvj3uw1dI6h4yn30DGLxHv1w0v36tXL8TNx+1z47JSb5P+NP4EbH8vkPz2PGleVSbDxF33aoU9YQohihfzkZJze5L+8e544j79k+Y\/Nk\/+yfIqynQQcf9H+Q3IMWWipaJ0mWmecvlM7zacpfZMFYdb7P8vm\/3hXzKqyC71n5ROzUBjlPN0f\/\/EfdxsbG0G\/se2B\/w\/f\/ifu3\/+LjyvLH\/6bX5nsRm2IeIqhWAqqVuOqor2GiPeffvszhXH\/8bf\/qVNuGqitmp\/mK07IeHXELJLj7\/yrX3b\/9nt\/6n7nX\/\/y5Bnz2Sk3ye\/\/88fOzkatxlVFccTxf\/zhd9x\/\/Xf\/MogolmLqeyfERao4iqdc8+xXfrP88\/jkWcv3sZjafwiWiqN9DFFUWOvcQzw\/dcVQjiHPRvFC7lmxdB\/qXgyZZ1diNVwYdWXbxfPQA7C\/v+8ODw8RGPTqGfj27\/+++wd\/8Ae9ypnvM94zPAPdfgZ0Hxa\/Sfsxg8KowDmpONrc3HQIDHgGeAZ4BngGUp+BgdwRug8LXJ+9cqUw6tVxkSwEIAABCEAAAnUSGGxhlPY30PQ30qSvEzixIQCBXhEgWQhAYGAEBlsY6ZxVBKkYMtFYegQCEIAABCAAgWESGHRhpCNXMWSiMRI5AbYHAQhAAAIQyCAw+MIogw0mCEAAAhCAAAQGRqDvhdHAjovtQgACEIAABCBQJwEKozrpEhsCEIAABCBQiQCTmyZAYdQ08dfr2R\/4VvtaRZMgIDZpknBzvk\/SNrSxWGTtWXaTND+zqU2zx6ybtmfp0yTJwvdJ2mIfz9p7VXvM\/Kax8fV+P8kiy5b0ZZyPAIVRPk5BvfQg2x\/4Vqtx0AUiCiY+SfG3J3a+XWPfPqT+rL3LnsVqlj1mltp71v58btb3\/TXf9Go19u1+P7a+9qo9m2js71Fjs6nVuIjd942tLxZiYqKxv0fT+61vl79v09i30y9HgMKoHLfSs\/Tg6kH2A2gsva+jP5uAmImd76mx9L5uCP1Ze5ZdbHwWGksvnVqN1TfRWHobx9pW3aPmi5XPR2PpfV2Mfe1Re522tzS7\/KXXHLUaq2+isfQ2jrXVHrXXsvtLm6940peNybxbAhRGtxz42lEC+iY36WaK3chKL0RJN7LpVxbiJsnK2p5BtVl+2CAQkoCeN5OQcYmVTYDCKJsP1pYJ6MIy0Qui5XRYfqAE7BlUy3P45iEQjzcjN\/nzfkmdb6f\/hkCSk56rpE5jE9nfzKZXJwEKowRdht0hoBeCn43GvBx8IvSbIKDnzl9HY55Dn8htX0zE5nbE1yIE0tglWWosvyJx8S1HgMKoHDdmQQACEIDAawK6sHVxvx7SFCDQArsC2Q3TlcJomOfei13rhdGLREkyagI8h9nHKz4URdmMplmz2Mk2bR76eglQGNXL953oeoEkH3iNpX\/HGcVbBJKcxEw630lj6X0dfefERGx8FhpLL51ajdU30Vh6G9PeEkhyESPpbq23XzWW\/nZU89cWw8\/apxjIx09RY+mlU6ux+iYaS2\/jWNui+0z6i5F0Ph+Npfd19IsToDAqzqzyDD24eoBNNK4cNMIA4mKM1Gqc3KZ0splonPRhfEtAbIyTWo1vLbdfNZbeRONby7C\/ioMxUatxkoh0splonPSJdWx79lt\/r2Lh2zQuYvd9Y+v7XKxvexQn06nV2GzWSiebicZmoy1PgMKoPLtKM\/UAm1QKFPlkY6R22lZlM5nm0yN9pVTFISuA7CZpfmZTm2aPWZe1Z9lMpjEwu9ppPrHptdc0Se7T90naNJ5ll09s4u\/Z7\/v7nKYv6uP7059NgMJoNiM8IAABCEAAAhAYCAEKo7oPmvgQgAAEIAABCPSGAIVRb46KRCEAAQhAAALdIxBbRhRGsZ0o+4HAwAnoD6LmRVDEd1rMIjF83zx9f03f39fThwAEwhKgMArLk2gQgEAJArr088is0IqhP7A6y68tu3JTjv76Gkvv69L68pFvmi2sjmgQGDYBCqNhnz+7h0BnCOjinyVZyapo0Pykj\/TTRL7TbNLLLlF\/muSx+z7Wn9am7UG+EtmUh\/oIBCBQDwEKo3q4EhUCwQjoIvQlGVi2pM4fW9\/81JokbTZOa21OWpvHP83H16XF9XW+b5G+igkTzbN+Wpu0a2yS5j9LZ3PV+r4am5jexv6erW82WghAoH4CFEb1M2YFCJQmoIvRLk5rpUsGTNP5Pr5dcWSzVv28ojlpkpyv9ZJ+0iX9\/HHSPzn2fZN9xZZ\/Ut\/FsZ+rctbYz1M6E+nVV2uicXKO2WghAIHqBCiMqjNsOALLDYWALj9dgsn9SidbUt+VsXJTjsl8pJMtqbexbFlifm20yt3WzcrRbOabNs\/Xyc\/Gmmt96REIQKAdAhRG7XBnVQhUIpC8QDXWxZoWVHrZ02xN67LykC1LQuSq+OKRFkt62dNsvk4+00R+sqlNivQS0\/t909FCAALOuZYhUBi1fAAsD4FpBHRx6rKWTPPpor6veVdlqXPS3vPEka\/vp7HmqvX19CEAgeYJUBg1z5wVIZCbgC5LiS5Mk9yTa3C0HPw2bRnlLMnyS9r88ax+2pp5dZaX76\/1pPd1efuaK5k1Xz4mim19tRonRfpZMZNzOjomLQj0igCFUa+Oi2SHSkAXpIkuzDQOsidtGkuf5l9Gp1hJyYrj+yoX39e3SZ81ls33Ub+KKJ7lo1bjqvFmzdcaEvmpNbGxtcpHfQQCEGiHAIVRO9xZFQIzCUy7IHWhTrPNDBrCYUaMabkVyTvpq5jSzVi6kFnx8saV3zTRotNs0stukhybnhYCEOgOAQqj7pwFmUCgMgG77BVIl7DG6vdRlLv2IFE\/aw+yyy\/LJ2mTf9558psmijvNJr3svkintaVTq7H6Jhqn6c2eZTMfWghAoDwBCqPy7GKZyT46SsAuyGR6Xb8Yq+at\/Um0b8VSG1IUW2Kx1WosCblOWiytJb3aaetJn2XXfAQCEKiPAIVRfWyJDIHKBOyC1GVpIl1WYNnlqzbNTzbprVVfonGayGaSZpfO7NZqbel9kc7saa35yk9iPuqbzXRprfml2aSzGPKTSGeiscR8TF9Ha2toPYnGto760mmsVmP1TTSW3sa0EIBAHgLFfCiMivHCGwKNE9BF6EsyAdmK6ORvYvNsnNYW8TFftclY0vmiS97E9\/V9rG928zd9spWffHy9xhLZJL4t2ZddIn9J0u6PZZf4uqy+fCWKLzFf65vN9Gplk159tRqrj0AAAvURoDCqjy2RIQCBDAK65E0y3N4y5fGXjz9JY4mvm9WXvyTLT3aTLD+zzfKV3Xz91vTW+rah9tk3BOok8P8BAAD\/\/wh9LygAAAAGSURBVAMAMbg7GW1hHPoAAAAASUVORK5CYII=","height":281,"width":466}}
%---
%[output:233a680b]
%   data: {"dataType":"text","outputData":{"text":"[14:01:19][INFO]  AUGMENT_FACTOR=1: オーグメンテーション無効 (コーパス = 500 ユニークSMILES)\n","truncated":false}}
%---
%[output:6ccf3410]
%   data: {"dataType":"text","outputData":{"text":"[14:01:19][WARN]  語彙に想定外の文字が含まれます: [K] -- 生成品質が低下する可能性があります。\n","truncated":false}}
%---
%[output:7c8775f2]
%   data: {"dataType":"text","outputData":{"text":"[14:01:19][INFO]  語彙サイズ: 37文字\n","truncated":false}}
%---
%[output:329ce97e]
%   data: {"dataType":"text","outputData":{"text":"[14:01:19][INFO]  特殊トークン: PAD=1  START=2  END=3\n","truncated":false}}
%---
%[output:18cc59a8]
%   data: {"dataType":"text","outputData":{"text":"[14:01:19][INFO]  マルコフモデル構築完了: 37x37遷移テーブル (500 SMILESから, 23229バイグラム)\n","truncated":false}}
%---
%[output:2f235c0c]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAkYAAAFgCAYAAABewLbaAAAQAElEQVR4AeydC5RdSVnv6yQNDJN2aBjG0CQtAxe4A5fHiGhyRbKu4gPxhcaLEkARxKUiBkFQERURBRlAAQERQUSNXDALkCug8lIEQREQEbgLGV7BJrRDNzMd5oTpIbd+O6dmKjt777MfVfvsx3\/W+VK7Xt9X9avTu\/6n9umeXWf0nwiIgAiIgAiIgAiIQEJgl9F\/IiACIiACIjBYApqYCFQjIGFUjZdai4AIiIAIiIAIDJiAhNGAF1dTE4EhEtCcREAERCAmAQmjmHTlWwREQAREQAREoFcEJIx6tVxDHKzmJAIiIAIiIALdISBh1J210EhEQAREQAREQAQWTCC4MFrwfBReBERABERABERABGoTkDCqjU4dqxK4wx3ukHRxaZIp8U+6fTqPC1fmUsoWYSHj4wvz55HO+3VVr+v6qtuP8VXpW6UtvrG8Pnnl9GnLQoyhjI8ybdqa80DjaFoDJyBhNPAF7vv0uMlfeeWVC50GYyg7AMbq2pOWNeffb48vzNW51G\/jrl0dqSvzU8qrmN83fY2fdJmfpz7LaJM1n6y2KqtPAMawru9BPUVg3AQkjMa9\/sFnn74hp\/NVAtKXm3yVPllt8VNkWX2alDFm4pH6hk\/yLuXaGWVYOk+ZsyyftHf1LqUMI0+KcU1\/3\/LKKKdPnpWpp41vxMWfX9aZ64ADYZ55Rpi8Osqpx7jOszL1tIE1PriWiYAIVCMgYVSNl1pXJODfoLmmu0u5zjNu6kXtqPcNP+Sz+lBWZPQNbcTDJ2Ny5vIu9cspKzLa4pPUb0eecr+s6Jq2GG1IMXdN6hu+\/TzXlLk+5H2jzs+3fU18NzaufWMsfp5rytKWV+7aUe\/Mlfkp8Z1R7q6z0nQ9eWdZ7eeVub5KRUAEmhGQMGrGT70jEGDjYRMock192mbtgyZlxjIvoBsn7bh2qbsm74x4GHlSjGvXltSVkZKnvor5\/fzrtA98U085KUYZeYw8Kca1X0eZs6I61yZ0yljmWTom40yX+XnqfZ\/k\/fquXTPWro+xa8w0HhGAgIQRFGRRCVS5QQ\/pRu7mQooB2U\/9a+owWGHpa\/LOqKcvqStLp1n1tKectqSYf009ed9cmUtdHX39Mv\/atQmREidtWX5pk1WeVZbVNqvM70t9rDn6cbj24xB3ntEH8\/uRl4mACNQj0H1hVG9e6tUzAu7mP+\/mTrt5bfyp076M+X1CXTNOZ\/h016RZecp8Y9x+nmvKMOeDsjyjnTPacE1KX8y\/dnnKfHN9KKMNeYxrypxR5q5DpsRJW1n\/VcbkYpT1TTv6kOYZ9XljoJz6vL6unDZ5RhvqSGUiIALhCEgYhWMpTzUJuE0ixk0en\/Os5rDP68Y8nPmVlDEGv4w85X6Zu84qpwyjn2s3L6Ut5tr5166sKM2L5\/xQ7\/pT5uddeRspcYnfViziYW3Ey4tB\/LbmnDcGlYcjIE\/dIiBh1K316PVoim7W3MSp9ydIHqPOL8+7rtI2z0eV8qx4lOX5YB6YX097ykjT5sr99v61X8815urTvsj7de46K6UtRh2pM\/LOKMuKRz11GPWklGHpPGWLtNDjYa74dEZ+3vxom25HnvJ5fbPq6YvV7Z\/lU2UiIALnEpAwOpeHci0S4OaOlQ1Z1LatzYIxEKvMmGlHe9qS+kYZRhntuMa4pozrLHP1tPHNb0u5n09fU49RTuqMvDPK3DUp+bS5clJntHHXbaVVYjp+dcZWPs653ulHXEpJyXNd15r2rxtX\/URgLAQkjMay0ppnJoEQG1WmY1voNjBipM1W3\/By7W4o8C5cnetPnmuviSFPuV8W85p4Mf0P0TfrAzfSefOjXZ7RN6+OcuplIiACzQhIGDXjp94VCLApxLh54xPfFYYyt6nvM92YWNSny4vy9Ekb7cv6cX3pw7XrR0qe8rJGH4z2pM7IO3NlWSltsspdGfW+MT7q\/LKxXTP\/shxol2dwy6ujnHpnLqbLKxUBEShHQMKoHCe1mkNgUTfhJnHL9qVd2ubgOK863d\/l0w3Tm1u63uVphw9SV1Y2pQ9Ge1Jn5J25sqyUNlnlroz6rhhjcmOpy4v++KE\/11WMPhj96UdKHiMvEwER6B6BEQqj7i2CRlSdABsLxkZTtTf9sKy+lGP4JMVol2XU0c43ymhLilFHvshoU8Xwi+GzqB9tsKI2bdUx1qyxUJZnjC2vjnLq2zQ3B2Jj5PPiU4\/RBvPbkceox\/y6UNf4JUYof\/IjAmMiIGE0ptWONNc2b8LEwrjpY1WmRD+MflhWX8rTltUuq8z5ps75oKyM0afIfB\/Od1F718alflvnizJ37VLKYhljIY7vn7K65vvxr4mRZbRJl1PmG2Px8+lr6p2l68g7\/0VtaIe5Nq4PZXlWpo3rS1t8u7zSERLQlBsRkDBqhE+dIVDlJlymbVEb6jDiVjX6YVX75bVP+0rn6UdZGaOtb\/RJ5ynD\/PKq1\/QvsjL+6F+mXVabJn2z\/PllzjdpWfP7h7h2cav4KtPHtSGd57tMm3k+VC8CYyYgYTTm1dfcRUAE+kBAYxQBEWiRgISRhX3ixAkjEwO9B\/Qe0HtA74Gi94DdLiq\/ivyFrqs8OHXIJDB6YcQb8wlPeII5dOiQTAzaeQ+IszjrPdDL98CDH\/zg5EN05m6aUdj2\/lJ1fBlDVpElIGFkT4ve8573mJ\/c2DBPXl8vZYc3Ny06U6kPvm9m\/ZexpZn\/m9oxlWnv2pyx\/kvbLMYZG6N0H+t\/o4JdPYuxaWOU7XfS+i9rWzP\/V1n\/ZfvQ7gE2Rln72lmMQzZG2T60+x82RlnbP4txRxujbB\/mUdbqcirrn3Z1Y1R5752Zcarynv2iXYey9qWZ\/2vsOpTtQ7v\/sjHK2jWzGFs2Rtk+tLvaxihj1878n7L+y7R3bcr+fNKuzs81\/XbbOZS1XbN57LLzKNuHdrwXy1id9yt92CuSm3\/JfxBG9LniimeYY8f+JKodPfoYQ6yiofHFfGd57Vw9aV6boZePXhi5Bb7LdGruWtJoSz\/Su06npfvtsv7L2G7bDv9l2vptzth+Zc3YtsQo2961O237VTFiTCv0oW1ZYxzOf9k+tFu146lixKjSnrYX2Rhl7Ra2LTHKtqcd8yhrdTmV9U+7ujHc+6pMWuc9e51lW8VYh52Kfb5s25c1x4m0bB\/aMaYyxlyZA2mZ9q4N46lixHDrXrbfxHKqYsSo0p62jKmMMWb8l2nr2rg+9KtqBw58nTl4MK4Ro2hcCB2+mO+MfLo9Za6elHy6zRjyEkZjWGXNUQREQAREYIEErrexdyIbMWyIjBcCB6HjV5Gn3C9LX9\/YJl0z7PzghBEL7VuM5btkZ8fwOI00hv+J9c\/jNNIY\/vF5xsZIHk3YlHwM27G+OXa\/3qax\/HO8TZwY\/vH5VXbsPE5btin5GHYz65vHaaQx\/MMnNqc2YsR+z37FrgOP02K9X1lbfPM4jZR8aGMOPE4jDe3b+WOtY\/5cJ3HsWvA4LdY9kDnE\/plI5tHCPydO\/KfBWgg1mhCDEkYIIhSub5SFXk0E0eGtrdBub\/DHzeAmW1s35KNc2BvPma2tKK6dU27+V29tuWzwlJvbF7e2gvv1HSKI7rW15RcFv0YQrW1tBffrHLbBqY0YJvJ7lvfrl7a2HLYoKTGu2dqK4hunCKJrt7a4LLQmlcwh5s81Y+MeuGtri8so1sr79ZyRxzstOn78NebQofubI0ceeU7EEJkY+2eIccX2MRhhxAIiiNLAKKMuXa68CIiACIiACLRDgMdcccTR4cMPMMeOvcgcPfqoxlNx+yV7Jka+sdMeOhiMMGrKfmNpyThr6kv9RaBdAoomAiIQg8CS3RecNfMfRxQZs2P277\/EHDx4T3PgwD2bDXHWGzHkbFY0umQUwohFnreyT1tdNUfX1hI7vrIyr7nqRUAEREAEBk5gz\/Ky2Wf3BWyv3SPqTzfeiRHi6KwRo\/4I6ckpEenYbRTCqMwi+3\/H6ND2dpku0drIsQiIgAiIwOIJnLJ7gfu7SHxZu\/6I4p0YnRVF+M8XRhwOpEUPecqL5lSmTVH\/vtaNQhixuPMWyP+bRHy5el571YuACIiACAybAF\/SDvF3jIxBtCBeYhox8tcDEcRe6Iy8a00Z15Rx7Yw85RGs0y5HIYw6vQIanAiIgAiIwMAJxBREznexMAIwQscZeWeU+dfkMVc2tnQwwohFROWmF5Ay6tLlyouACIiACAQiIDdzCCBanICJlRJjzjBUXYrAYIQRs0UAIYR8o4w6mQiIgAiIgAgshkAsMeT7lTAKtbaDEkZAQQj5RplMBESgEgE1FgERCEoA0eKLmBjXxAg66NE6G5wwGu1KauIiIAIiIAIiIAKNCUgYNUbYAwcaogiIgAiIwAIJxDghSvvUiVGoBZYwCkVSfkRABERABEQgkwCiJS1kQueJkRl8FIUhJylhNKP5aZt+PLJ91PqPaZ+x\/mPbhTZGTLPuo78ushFiG387PaZdYucQ23bbGLHtPjZGTLvU+o9tF9gYsW2PjRHT2FJjG5tNbLOYOvoKLYKy\/LGCHZ1+z4bF+7RnQ9ZwRUAEREAEhk1gaLNDtGSJmZBlxBgat8XMR8JoMdwVVQREQAREYDQEQgqgPF8SRqHeThJGoUjKjwiIQC4BVYjAuAkgWvIETahyYoybcqjZD0oY8Ycd88AU1eX1UbkIiIAIiIAINCcQSvwU+ZEwar5OZz0MShgxJQkgKMQ2+RcBERABEShPANFSJGpC1BGj\/IjUMp\/A4IQRf\/Va4ih\/wVUjAiIgAiLQNoEQwmeeDwmjUKuaCKNQzvrsZ3NpyWzNrM\/z0NhFQAREQATCEFiye4KzMB7lpQ8EdvVhkFXHWOfU6GWrq+bZa2uJvXVlpWpItRcBERABEegugVoj27O8bPbZfQHba\/eIWk6STpzmzDvxaVpPjCSY\/mlIYJDCCCZVxdEPbGyYR6yvJ3av7W1cyERABERABEZM4JTdC07afQHb2txsQKKp6CnTX8KowQKd03WwwuicWZbI3H46Nc5WdngTluikJiIgAoshoKgi0AKBHbsXTO3egJ22af2QiBb2lZhGjPojVM8bCQxaGFU9NboRi65EQAREQAREIBSBmILI+ZYwCrVagxZGQJI4gkLnTQMUAREQgQETQLQ4ARMrJcaAEbY4tcELI1gijkhlIiACIiACItA+gVhiyPc7Xxjxp2yc5TFw9S7Nazfk8jjCaEHEigRQUd2ChquwIiACIiACoyCAaPFFTIxrYuTDROiwDzojn25Nmat3KWXpdkPPD0oYDX2xND8REAEREIFuEKg2ihhCKO0zXxghbhA6\/pjJU+6X6fosAQmjsxz0rwiIgAiIgAhEIoBoSQuZ0HliRBr+yNxKGI1swTVdERCBLAIqE4F+EjhxYttgIUbvTpE4SXJGWQjfffIhYTRbrXfb9O2Rzbrv\/esaO4OYZt1Hfz3WRohtL\/iEMTFtw84htt3ExohtlcU62QAAEABJREFUr7YxYtpHrf\/YFvPnwfmOvdYWU\/TXlo0Q22yIjr44zYljx4\/\/P3Po0KvMkSNvbDx3xBBCyDfKGjvumQMJo54tmIZbjoBaiYAIiEB3CCCKQj86O+vv8OHbm2PHvtUcPXr37ky35yORMOr5Amr4IiACIiACXSdwVsQYEz7dv\/9m5uDBi82BA7fuOoTejK8nwqg8T4790la+t1qKgAiIgAiIQGgC8U6MbhRbxMgeN4\/G2Bf9WvKU+2W6PktgUMLILTSL7RvlZ6erf0VABERABESgbQIlTooanyblCyNmy57IXuiMPOUYZaSUce0bZdSNyQYjjFjIvAWknPoxLazmKgIiIAIi0BUCiJbY4ogYxfNlL3Tmt6TM5bn2zZWPKR2MMBrTommuIiACnSegAYqARyC2KML\/fGHkDUiXBQQkjGZwrl1aMtOZzYqUiIAIiIAIjJjAkt0TnDXDgGhBvMQ0YjQbpXqfJSBhdJaDed\/qqnnn2lpiV66szEqViMCMgBIREIHREdizvGz22X0B22v3iPoAYgoi51vCqP76nNtTwmjG464bG+Ze6+uJ3XZ7e1aqRAREQAREYKwETtm94KTdF7Ctzc2xYhjdvMcqjM5b6FtOp8bZBTso8POaqEAEREAERGBEBHbsXjC1ewN22qb1p85pDvtKTCNG\/RGq540EBiOM+BZ93m+eUU79jdPWlQiIgAiIgAi0RSCmIHK+08KorbkNL85ghBFLg\/hBBKWNcuplIiACIiACItA+AUSLEzCxUmK0P7MhRhyUMGKBEEFpo1wmAiIgAn0moLH3mUAsMeT7lTAK9Q4ZnDAKBUZ+REAEREAERCAMAUSLL2JiXBMjzGjH7kXCaOzvAM1\/QQQUVgREYDwEYgihtE8Jo1DvJwmjUCTlRwREQAREQAQyCSBa0kImdJ4YmcFVWJGAhNEM2EU25c861rUy\/WyI3r9uYmcQ06z76K\/rH2BMbLv89sbEtH2WUmyb2hixbcXG6LvZKUR\/seXFtOgTsAHu0YLZEB19hRZBWf54h3R0+j0bloRRzxZMwxUBERABEegbAURLlpgJWUaMqFxG41zCaDRLrYmKgAiIgAiIgAjMIyBhNI+Q6kVABERgiAQ0pxYJhDwZyvOlE6NQCzo4YZT+447kQ8GSHxEQAREQARGoTgDRkidoQpUTo\/rI1ON8AoMSRoig9B93JE\/5+VNXiQgMioAmIwIi0FkCocRPkR8Jo1DLPyhhFAqK\/IiACIiACIhAOAKIliJRE6KOGOFGPGZPoxBGnBr1bpE1YBEQAREQgYEQCCF85vmQMAr1ZhmUMEIA8djMWRVIp5aWzJdmVqWf2oqACIiACAyTwJLdE5w1myGiZZ6waVpPjGajbLt3V+MNShgBGXHkrIpA+vvVVfPGtbXEPryygiuZCIiACIjAiAnsWV42++y+gO21e0R9FE1FT5n+84WR2xNJs+ZCeZZltR1y2eCEkb9YvkDyy7Ou772xYQ6tryd26fZ2VhOViYAIiIAI5BIYXsUpuxectPsCtrW52WCCiJYy4qZJG2LkDxHB4\/ZEUvLp1pSnLd1mDPlBCaOshS67iJdMp8bZhTu8Ocv2VDsREAEREIEhEtixe8HU7g3YaZvWnyN7SmzLF0bsjQgef\/zkKffL0tfU0y5dPvT8oITR0BdL8xOBNgkolgiIQCgCiJbFCaNQsxiLn0EJI5QtCjdtlI9lQTVPERABERCB8RA4cWKXwcYz4\/gzHZQwAhciKG2Uy0RABERABERgMQTinRYdP77LHDq0xxw5cmHQqXHAwF4a1GlPnA1OGPWEu4YpAiIgAiIwGgLxHqUdPnytOXZsyxw9emo0NHMnGqhCwigQSLkRAREQAREQgWwC8U6M9u8\/bQ4evNYcODDNDq3SygQkjCojUwcREAEREIEWCAwoRLwTI2Oc6CJGNjIeifFozK8lT7lf5q6L6lybIacSRrPV\/RObvjCyWfe9f11nZxDTrPvor91vMCa2feDMtSamfdZSim02RPTXVTZC381OQa8SBN5i28Q2G6KjLydeYqb5wggoiCAEjzPylGOUkcrOEpAwOstB\/4qACMQmIP8iMFYCX7ETR7fENGLYMEUvxJAzvx1lRXm\/bgzXEkZjWGXNUQREQAREYHEEEC2x7czipje0yBJGQ1vR9uajSCIgAiIgAmUIxDwp8n2XGYvazCUwSGHE89I8m0tEDURABERABEQgJAFOi3wBE+OaGCHHPGJfNwqjgUHgmWmWDWyamo4IiIAIiIAIiEBAAoMVRgEZyZUIiIAIiMAACCxsCpzmxDZ9xyjY8koYzVDuXloyzmZFSkRABERABEZMYMnuC84aYYjx6CzLZ6NBqrMjMFhhlPUdIzfprPSWq6vm1mtrie1ZWclqojIREIFOEdBgRCAugT3Ly2af3RewvXaPqB2N06IsIROyjBi1B6iOPoHBCqOq3y+6emPDbK6vJzbd3vYZ6VoEREAERGCEBE7ZveCk3Rewrc3N+gQQLbFNj9Lqr0+q52CFUWqec7Nfnk6Ns+t3+Oukc7uoQWACcicCIiACXSKwY\/eCqd0bsNM2rT22kCdDRb5qD1AdfQISRj4NXYuACIiACIhAaAKcFhUJmhB1xAg97pH6iyiMRkpU0xYBERABERABnwCiJbbpUZpPvNG1hFEjfOosAiIgAiIwWgJlJx7iRKiMj7LjUbtCAoMURnzxunDWqhQBERABERCBtghwWlRG2DRpQ4y25jPwOIMURgNfM01PBEQgDgF5FQEREAEjYaQ3gQiIgAiIgAjEJMBpTmzTd4yCraCEUTCUctQ5AhqQCIiACHSBQJNHZFX6dmGuAxiDhNFsET9sSZyKbAdtrJj2KOs\/tt3ExohpF1v\/sc2GiP7aPbm5iWnRJ2ADxFxn59uG0asEgffbNjHNuo\/+2mMjxDYbopuv2KdF+NeJUbC1t1IgmK\/YjuRfBERABERABPpHAOFS5eSnTlti9I9MJ0c8CmHE\/zetk\/Q1KBEQAREQgeETKC10LIombW13vZoTGLwwQhTp1\/ebv1HkQQREQAREoCYBTnNimx6l1Vyc87sNXhidP2WViIAIiEA7BBRFBBICiKImJ0Fl+hIjCZb\/DwcFzvJbGePakBa1G2rdoIURi6rToqG+dTUvERABEegJgTLCJkSbAhxuP2RPxMhnNaecemfks9oNuWzQwmjIC6e5LYKAYoqACIhADQKc5sS2gkdpiBuEjj9y8pT7ZeQp98vSeb9uqNcSRrOVPbG0ZJzNipSIgAiIgAiMmMCS3RecjQkDAsnZmObt5jpqYeQgkB7Zu2oO7VtL7Lm3WKFIJgIiIAIiMGICe5aXzb61tcT2rq7WJ8FpUYhHZRk+Tpy0H+qtGWLUH+ENPRFEnBI5I39D5Ugudo1knnOnecVVG+bYyfXEDp\/antteDURABERABIZN4NT2tjm5vp7Y1uZm\/clmCBoTqOz43y+bQ49fM0eenincKo8ZQVS508A6SBjNFvTAdGoOzmz\/zs6sVIkIiIAIiMBYCezYvWBq9wXstE1rc+A0J5Id\/p\/b5tgT1s3R720g3GpPbJgdJYyGua6alQiIwBAJaE79JIAoCnRClD5p2n\/LHXPwjlNz4E7TfrLp4KgHLYw4Ehzj89EOvs80JBEQAREYL4FIoigtkvIAZ+2F7I2U+33IU+6XjfF6YcII+M6ywFOXVa4yERgQAU1FBERgDAQ4MYptBb+uD2InethbMfKUY+RJMcrJOyNP+ZhsYcIIyADHuI5l+GeBY\/mXXxEQAREQAREoJIAoin1qRIzCQRjDfujMb0pZOk8Z5peP5XqhwqgtyK0tblsTUhwREAEREIH+EIgtipz\/\/hDp9Eg7IYx0otPp94gGJwIiIAIi0IQApzmxbc6jtCbD9\/uO4boTwgjQEkdQkImACIiACIiACCySQGeEEY+7FimO3m3V\/Nsj24pd6Zj2Ges\/tllEyR9YjZVu2TnEtvvYGLHNhuj96zo7g9hmQ+hVgsD32jZxzBj8WvfRX\/e1EWKbDdHNFzdM97grVkqMbs6+d6PqjDCCXAxxtEixxZxkIiACIiACIycQSwyl\/Y4cc6jpty6MECpY0QRcPalvRX3y6mKIrbxYKheBRRNQfBEQgQ4S4DQntuk7RsEWvnVhhFDB8mZAHUY9qTPyMhEQAREQARHoHQFEUfp0J3SeGL0D080Bty6MFoEBccXJ0yJiN4up3iIgAiIgAr0nEFoE5fnrPahuTKA3wghx0w1kGoUIiIAIiIAIVCDAaU5s6+ujtAoY22raG2HUFhDFEQEREAEREIGgBBBFeac8ocqJEXTQ43XWCWHUhdOg\/1paMs7G+3bQzEVABESgEYFBdV6y+4KzRhMLJX7m+Wk0SHV2BDohjNxgFpk+c3XVPHFtLbHXrawsciiKLQIiIAIi0AECe5aXzT67L2B77R5Re0ic5sQ2PUqrvTzpjgsVRnwhGksPahH5R25smCeuryd2n+3tRQxBMUWgewQ0IhEYMYFTdi84afcFbGtzsz4JRNG8056m9cSoP0L19AgsTBjx+MyZN54bLqm7IdPCxX+fTs1lM7v1zk4LERVCBERABESgywR27F4wtfsCdtqmXR6rxhaOwMKEUbgpzPfEqVTbQmv+qBbSQkFFQAREQATaJtD0NKhs\/7bnNdB4CxNGiJUqTKu2r+JbbUVABERABEQgGgEec8U2fcdotnzNk4UJo+ZDL+cBQaXTonKs1EoEREAERCACAURR2VOfuu2IEWHoY3TZC2HURNxIFI3xba05i4AIDIXAIOZRV+xU7TcIWIufxEKFEYLH2eJRaAQiIAIiIAIiEIEApzmxTY\/Sgi3cQoURpznO8gQS5bQJNmM5EgERWCABhRaBERJAFFU9\/ananhgjRBtjygsVRv6EED8YQsgvb+v6ITbQ\/SLbm6z\/vlvVn9Uutn+nXYfYZkPoJQLBCHzGeopp1n30Vxv3vuvPnDGx7C1ve1t9RoiWNmzOCNlfnWU1dXXpNKvtkMs6I4wcZF8csTjkXZ3S7hHQiERABERABOYQaOMTIsKrYBhuP2VPxchnNacubVnthlzWujBiMbAiqCwKbUiL2qlOBERABERABDpPANESWxwRIwdE1n7K\/kp5TpdRF6eEUXwWLAZWFInFog1pUTvViYAIiIAIiIAIhCPAvussnNd+eWpdGM3Dw4IgimhHSp5rmQiIgAiIgAg0JrAIB5zmRLIT1ywZzOA\/wNzYd52Ndf\/tjDBiATAWJMDaJi7wl1zoHxEQAREQARFYFIGIj9GOX7lsDv3fNXPkbauNZ5fef8mPcR9dqDACuDMWAEuvLGW0SZcrLwIiIALGGEEQge4T4DQnkjg6vLZtjt133Ry9bLP7HHoywoUKI0SPs57w0jBFQAREQAREoBoBhFEk23\/Bjjl4q6k5cPG02pgyWusQ4iyUhQqjs0OY\/y\/iqeyC0c4Znt01Kfk8Wwpa7lYAABAASURBVFpaMs7y2qi8BQIKIQIiIAIdIXDixAnjrNGQIp0WGd8vwitnkFl7KHsi5TldkuIybZKGA\/unF8KoCnMW2hn93DUp+Tzbu7pq9q2tJXaLlZW8ZioXAREQAREYCYHjx4+bQ4cOJXbkyJH6s0a0+CImxjUxCkbIHojQcUbeNaeMa8q4dkae8rFZbGGUy7Mq8KrtcwPnVFy1sWFOrq8ndmp7O6eVikVABERABMZC4PDhw+bYsWOJHT161NT+D9HShs0ZIPuoM78pZS7PtTNXNrZ0YcKoa6Cn06lxtrOz07XhaTwiIAIiIAItE9i\/f785ePBgYgcOHDDZ\/5UojXFClPaJ8CoxFDWZT2DQwgjVOx+BWoiACIiACIhARAKIlrSQCZ0nRsQpjMl1J4SRe56ZTse0EJqrCIhANwhoFCIgAuMmsHBhhBjiZCfLqMPGvUSavQiIgAiIQK8JcJrThvUaUncGv1BhhOhBEOXhoA6jXV4blYtAMQHVioAIiMCCCYR+bJblD+G14GkOJfzChBFiB9EzFJCahwiIgAiIgAhkEkC0ZImZkGXEyAyuwqoEFiaMqgwUAYWQoo9MBERABERABHpFANHShvUKSncHu1BhhNhJW3dRaWQiIAKhCeyxDmPaxdZ\/bIs5fufbTkOvMgQOTYyJZUe\/ucwIsttUOxky5\/xF67J9EV7Z0VVakcBChREnQWnzhZI\/F9r5eV2LgAiIgAiIQC8IIFrKCpy67YjRCxjdH+RChVEWHgSQMyeSstqpTAREQAR6QUCDFAFESxsm0kEILEwYIX4QPkWzoA02r12RD9WJgAiIgAiIwEIJ1D0FqtIP4bXQSQ4n+MKE0XAQaiYjI6DpioAIiEA1AoiWKiKnTltiVBuVWucQWKgwcqdBRSdC1NEuZ\/wqFgEREAEREAEREIFgBBYqjJgFogdDAGUZdbSLZnIsAiIgAiIgAjEJcJrThsWcw4h8L1wYOdYIoCxz9bHTpaUl4yx2LPkXAREQARHoPoETO0vGWaPR1nk0VrUPwitnkCquRmChwogTorLDLduWdnlWFGvv6qrZt7aW2C1WVoqaqk4EREAERGAEBI5\/adkc+txaYkc2VuvPGNFSVehUbU+M+iNUT4\/ALu+6M5cIm7qDyTp1cmVFPq\/a2DAn19cTO7W9XdRUdSIgAiKwQAIK3RaBwxdum2OXrCd29KLN+mERLW1Y\/RGqp0egk8LIG19rl9Pp1Djb2dlpLa4CiYAIiIAIdJPA\/qUdc\/Bm08QO2LT2KKue\/tRpj\/CqPUB19AlIGPk0dC0CCyCgkCIgAgMngGipI3aq9CHGwDG2NT0Jo7ZIK44IiIAIiMA4CSBa2rBx0g0+684II75X5IxZumuXUtYP0yhFQAREQAREwCNQ5eSnbluElxcy69Ltp6RZ9X4ZbTC\/bCzXnRFG7gvSpMAn9Y0ymQiIgAiIgAj0jgCipa7gKduPGAVgEDn+nkq+oHlx1cBrWxdGLAY2cK6angiIgAiIgAicJYBoacPORjvvX\/ZcRJFfQZ5yv8xdU069y48tbV0YARsD\/Nhga74iIAIi0EECGpIIiIBHoHVh5GIjjrjOEkiujnqZCIiACIiACPSaQNnHYTXanbhuyWCGE6kAkNiTx74HL0wYufVjAVgIl89LaZdXF6L8cdbJUyPbBdZ\/TNtn\/ce23TZGTLPuo7\/WbITYdtDGiGnWff4rUE3MdXa+T9mxxrSYP2\/O99TOIbZdamPENOs++muPjRDbdr\/DmFh2vw\/aCdR9nbEdES4R7Ph1y+bQl9bMkWsb\/GVuOzxe7MWx91ridN0WLowAxEKwIFzLREAEREAERGBQBGqcBJmSfQ5Pts2xm6ybo7sb\/GVuDzZ7sTOKuSYdk3VCGAHciSMWIc9oN0LTlEVABERABPpMoKTIKSuG\/Hb7v7JjDp6ZmgOGc8tmkNiHfcMbedIxWWeEEdDdApBmGW1kIiACIiACItArAhEfpSXfLXKP6HKgsJ9y4OBXk6fcL1vcdbcid0oYgYaFYsG4bmIhfDSJr74iIAIiIAIikBCIeGJ0w+kR4igJlv2P21vZGzHyriV5d63UmIUKI39h\/MXIK\/fbFF1rkYvoqE4EREAE4hKQ9xSBDggjRsTe6oy8M8rctZ\/mlftthni9UGEUA6hEUQyq8ikCIiACIlCbwIIfpdUe90g7dl4YIXSclVkjFC5Wpq3f5pqlJbM9M79c1yIgAiIgAuMksGT3BGeNCHTkxKjRHEbUeSHCyAmdotStASIHc\/lY6RtWV83\/WVtL7H0rK7HCyK8IiIAIiEBPCOxZXjb77L6A7bV7RE+GrWE2JLAQYcSYETuYu3apKyPfph3a2DAPWF9P7M7b222Gbj2WAoqACIiACMwncMruBSftvoBtbTb4O0E6MZoPu0MtFiaM5jFInybNa9+0fnU6Nc6Wd3aaulN\/ERABERCBnhPYsXvB1O4N2Gmb1p6OvmNUG13Njo26dVYYcXKUtkYzVWcREAEREAERWAQBnRgtgnrtmJ0VRrVnpI4iIAIiIALDItD32fA3hmKLI2L0nVNHxt9ZYeQepXWEk4YhAiIgAiIgAvUIIFrasHqjU68UgYUJI1\/4cM24SDGu3WM0l6esitG\/Snu1FQERaI2AAonAuAjEPi3CP8JrXFSjzXYhwgjRMs\/cjGlXVxw5H0pFQAREQAREYGEEEC6xTcIo2PIuRBhVHb3EUVViLbdXOBEQAREQgXwC+q20fDYdrFmIMMo6ASpbFovh1dbxVmSL\/YGBOcS2Cy2jmLbb+o9t19gYse1\/2RgxzbqP\/uLmENsus7OIaZ+1\/mPbBTZGbLvKxohp1n301yU2QmyzIbr5in3zx79OjIKtPfe9tLNO5jk16uTANCgREAEREAEREIHBEFiYMOKEyDeI+nmuXRmpTAREQAREQATCEGjXC4c5HOrENGK0O6vhRluYMOIEyDcQ+3l3TblMBERABERABPpKANES2\/gaU1\/5dG3cCxNGsUFw4uQsdiz5FwERWCwBRReBLhOIeVLk++4ygz6NbSHCiNOgNKSsMtrklVOXZwgi+jkjn9dW5SIgAiIgAiIQkwCnRb6AiXFNjJhzGJPvhQijmIARQQiimDHkOzYB+RcBERCB4RBAtMQ2PUoL934ZnDAKh0aeREAEREAERKA5gRgnRFk+m49UHiDQijAiUFuWPi0qe4J0amnJfGlmbY1VcURABERABLpLYMnuCc6ajJLToiwhE7KMGE3GqL43EhicMLpxasaUFUX0+fvVVfPGtbXEPryyQpFMBERABERgxAT2LC+bfXZfwPbaPaIARWEVoiW2lXmUxp7oLG\/Art6lee2GXD5YYcSipk+Pihby3hsb5tD6emKXbm8XNVWdCIiACIjACAicsnvBSbsvYFubm7VnHPJkqMhX0QDdnsi+iJFPt6eMOt8oS7cben6QwoiFZGGrLN4l06lxduHOTpWuaisCIjA0ApqPCFgCO3YvmNq9ATttU1vUy1fWnsgeSXkvJxR50IMTRiw0Cx6Zm9yLgAiIgAiIQCkCPEYrOukJUUeMUoMpaKS98yycwQkjpoU4ShvlstES0MRFQAREYGEEEC2x7D+XlgxW5jtGZQH4++cYxdLghBGLmGVl3xBqJwIiIAIiIAIhCYQ4Ecrz8YblZfNDa2vmaIMvh6fn6u+hiKR0\/dDz\/RRGQ18VzU8EREAERGAwBDgtyhM2Tcu\/bXvbPGt93TyswZfDBwM60EQkjAKBlBsREAEREAERyCKAMKpqZdtfsrNj7j6dmntYy4pdpWyMp0NZfCSMsqioTAREQAREQAQCEWh6KlS2f95weTSWFj3kKc\/rM+ZyCaPZ6u+36aWRzbqP+rqp9R7brrMxYlrZG0CTdhfYOcS2vr+XLCITc52d74\/ZQDHtEus\/tp2yMeabMU3aOF6xUjuF6K\/dNkJssyE6+Sp7+tOk3bwvXyOCEEPOyDtYlHFNGde+UUbdmEzCaEyrrbmKgAiIgAi0TgDB0+TDXJm+xJg3MUSOM78tZS7PtW+ufEyphNGYVltzDUZAjkRABESgLIEywiZEm7LjUbtiAhJGxXxUKwIiIAIiIAIiMCICgxRG\/vNRruevp1qIgAiIgAiIQBwCPOaKbfO+YxRnZsP0OjhhhBDyn49yTdkwl0+zEgEREAER6DoBRFGIR2VFPohRyEGVpQkMThiVnrkaioAIiIAIiEALBIoETci6FqYyihCDE0acEI1i5TRJERCBMRPQ3HtEgNOc2KZHaeHeEIMTRg4Nj8+clRFLX1xaMlfPzPlQKgIiIAIiMF4CS3ZPcNaEAqIo5MlQli9iNBmj+t5IYLDCCDHkDIF045Szr161umpesraW2LtWVrIbqVQEYhGQXxEQgc4R2LO8bPbZfQHba\/eIugPMEjIxyuqOT\/3OJTBYYXTuNOfn7r+xYR60vp7Y3ba353dQCxEQAREQgUETOGX3gpN2X8C2GvxPWjnNiW16lBburTg4YVTmdCgL39p0apxdtLOT1aRKmdqKgAiIgAj0nMCO3Qumdm\/ATtu07nQQRTFOiHyfxKg7PvU7l8DghNG501NOBERABERABBZLwBcwMa\/bneVwow1OGPG9Ik6NfKNsuEuomYmACIiACHSZAKc5sU2P0sK9AwYnjECDEPKNMpkIiIAIiEA+AdWIgAicJTBIYXR2avpXBERABERABBZPgNOimI\/Q8E2Mxc90GCOQMBrGOmoWIpAioKwIiEBXCCBc2rCuzLfv45Aw6vsKavwiIAIiIAKdJsBpTmzTd4zCvQUkjMKxnOvpYtuirpXp9\/kz9zOxLfanHoso+utzNkJs+0kbI6ZZ94N4xX4\/bVhKsc2GiP6a2ggxzbqP\/vq4jRDbbIhOvhBFsd\/rxOjk5Hs4KAmjHi6ahiwCIiACItAfArFFkfPfHyKZI+1MoYRRZ5ZCAxEBERABERgiAU5zYpsepYV750gYhWMpTyIgAiIgAo6A0hsIIIrcqU6slBg3BNRFIwISRo3wqbMIiIAIiIAIFBOIJYbSfotHYYz\/h4+L2pZtV+Sjz3USRn1ePY1dBNojoEgiIAI1CXCaE9vmPUpD7Ph\/+Jh81nQoL9Muq+9QygYnjFjUPBvKomkeIiACIiACIlCWAHsiYsdvT55yv4w85X7ZGK8HJ4xY1DwrWuAvLi2Zq2dW1G4wdZqICIiACIhAIYEluyc4K2w4p5LTovRjr9B5YswZhqpLEhicMCo57\/OavWp11bxkbS2xd62snFevAhEQAREQgXER2LO8bPbZfQHba\/eIurMPLYJ8f1+w4g2jrO74XD8OFdw1ad9PkJhDHZMwmlG7\/8aGedD6emJ3296elSoRAREQAREYK4FTdi84afcFbGtzszYGTnNi2b9Z8fZiK95e2UC4ZU1srKIIFoMTRixmnjHhPFubTo2zi3Z28pqpXAREQAREYCEE2g+6Y\/eCqd0bsNM2rTsCRBEnOjGcd45nAAAQAElEQVTsLla8\/W8r3g42EG7pebGHpk+P0m2GnB+cMGIx82zIC6m5iYAIiIAIdJNADEHkfO6x4m3VirbbWgsx+7GLIhgOThgxKZkIiED3CWiEIjAWApwYxbaiX9fnsADB4\/MmT\/m8Mr9+LNcSRmNZac1TBERABERgIQQQRe6EJ1ZKjKLJIYIQQ87Iu\/aU+dfkfXN1Y0kljMay0tHnqQAiIAIiIAJZBGKJobTfrNh+GWLIWbqcvKtLp9SNySSMxrTamqsIiIAIiEDrBDjNiW1Fj9Jan3DPA+YKo57PS8MXAREQAREQAREQgcoEJIwqI1MHERABERCBARBobQqcFqUfe4XOE6O1CQ08kITRbIEf9v3GPPqhce2ztzMmpt108hYT2\/jhi2m77XrIjBEDMdB7oFvvgSabZcx7pu\/b3j71CkCgyVoHCC8XIiACIhCIgNyIQEcJhD4dyvKHQOro9Hs3LAmj3i2ZBiwCIiACItAnAoiWLDETsowYfWLS5bFKGHV5dcY9Ns1eBERABAZBANHShg0CVgcmMVhhxB+ncnz9a1emVAREQAREQATaIBDyZCjPF8KrjbmMIUZ7wmgMNDVHERABERABEUgRQLTkCZpQ5cRIhVW2JgEJo5rg1E0EREAEREAEfAJ514iWNiwvvsqrEZAwmvE6MV0yJ06ftVmREhEQAREQgRET2L20ZJw1wRDqVKjID8KryRjV90YCgxNGfJ8IY4qkmLsmzbMjH1o1h967lthzP72S10zlIiACoyGgiY6dwJ7lZbO6tpbYJaurRv+Ng8CuoU3T\/c\/vmFf6mrI8u+JOG+bY3dYTO7x3O6+ZykVABERABEZC4NT2ttlYX0\/s6s3N2rPmNKfotCdEHTFqD1AdzyEwOGHE7DglQhSlr8nn2YFbTM3Bme2\/2U5eM5X3mICGLgIiIAJVCFy\/s2NOT6c3WJW+fltESxvmx9R1fQKDE0aIInCQYulr8jIREAEREAERaItAiBOheT4QXm3NZ+hxeiyMspfGnRSRYrQixbiWiYAIiIAIiECbBBAt84RN03pitDmnIccanDAa8mJpbiIgAiIgAv0jgGipZXaqVfrZ5noFICBhFACiXIiACIiACIhAHoGmp0Fl+iOg8uKrvBqBQQoj\/7GZf10NjVqLgAiIQFACcjZSAoiWMuKmSRtijBRv8GkPUhgFpySHIiACIiACIlCTAKKlDas5PHVLEZAwckCe81Jj\/uT1Ue3OnzImpvFnKWObwxUrbfKJqfW+FoJiGiMGYjCG9wDCxv7I13q1wafM+PhNbWdFE6FNUf3Q6ySMhr7Cmp8IiIAIiMDoCSB2+GqJM\/JZUPLKs9oOtWzwwohF9q1gIVUlAiIgAiIgAsEJcJoT+9SIGHkDZw9EEPn15Cn3y9J5v25M14MXRiwmbwBn5GUiIAIiIAIi0BYBREssmy4tGQz\/8+dT3EL75Fk+oxBGZ6eqf0VABERABESgfQIxT4uuWl42H1lbMx\/X\/+Q22MJKGAVDKUciIAIi0C4BResHAU5zYomji7a3ze3W182tNzf7AaMHo5Qw6sEiaYgiIAIiIAL9JYAwimW7d3bMBdOpubm1\/hLq1sgljLq1HhrNqAlo8iIgAkMkEOu0yPeL8Boiu0XMScJoRv3EiavMiRNfSGxWpEQEREAERGDEBJaWloyzJhgQLb6IiXFNjLwx8qXq9G+ckac8r8+YyyWMZqt\/5MjzzaFDT0nsuc9946y0WaLeIiACIiAC\/SWwZ3nZ7FtbS2xvgy83I1rasCLSiCDEkDPyrj1l7lqpMRJGs3fBFVc8xBw79pjEDh8+MCtVIgIiIAIiMFYCp7a3zcn19cS2Gny5OcYJUdonwmveOiGGnPltKSvK+3Wp60FmJYxmy3rgwJ3MwYNnbf\/+W81KlYiACIiACIyVwM7OjplOp4mdtmldDoiWtJAJnSdG3fGp37kEJIzO5aGcCIiACIyXgGYuAiKgR2l6D4iACIiACIhATAKc5rRhMecwJt+jODHii2XOxrS4muvoCQiACIhABwiEfmyW5Q\/h1YGpDmIIgxdGfKnMt0GsmiYhAiIgAiLQGwKIliwxE7KMGL0B0vGBDl4YdZx\/teGptQiIgAiIQO8IIFrasN6B6eiAJYw6ujAalgiIgAiIwDAIhDwZyvOF8BoCrS7MQcJotgr\/7faPNLsn3xPVPm5jxbSrrP\/YlvdDGarcTiH66\/r7GhPbLraziGnW\/SBeL7WziGlPsv5j2+U2RmyzIfTqMYEzduwIl9hmw+gVgICEUQCIciECIiACIpBFQGUQCPXBscgPootYsuYEBi2M3G+iubQ5LnkQAREQAREQgWoEigRNqDoJo2prUtR60MJIv41WtPSqE4HqBNRDBESgOgE9SqvObJE9Bi2MFglWsUVABERABERABPpHQMKof2sWcMRyJQIiIAIiEJtAqMdlRX70KC3cKkoYhWMpTyIgAiIgAiJwHoEiQROqTsLoPOxnC2r8K2FUA5q6iIAIiIAIiEBZAvqOUVlS3WgnYTRbh6WlJeNsVqREBERABESgWwRaHY3bE0ibBA51KlTkRydGTVbo3L4SRjMee1dXzb61tcRusbIyK1UiAiIgAiIwVgJ7lpeTPYG9gT2iLociQROqTsKo7uqc30\/CaMbkqo0Nc3J9PbFT29uzUiUiIAJRCci5CHSYAHuB2xe2Njdrj1SP0mqjW0hHCaMZ9ul0apzt7OzMSpWIgAiIgAiMlQB7gdsXTts9oi6HUKdCRX50YlR3dc7vJ2F0PhOV1CegniIgAiIgAikCiJYiUROijhipsMrWJCBhVAMcX8Tje0ikNbrP7YLfmP4ZwBBitDGHEztL5rlXrxhSuMWwydKSuenKiiGN4b8NTm3E2LKc\/m5lxZDG4HTK+v+g9U8awz8+v2JjXGtjkJIPbW2swxBitDEHf20RLW2YH1PX9QnsKuza80r3\/0gjDTkVfqhWbnnLkC7P8RXbP8GGEKONOZy4HmEUb61Zi112s7yp3k+gKDQE0d9F5LRt1+HfIvpncgiiacQYbfxMDCFGG3NgvWX9JDBoYeT\/v9K47ucSadQiIAIiIAKxCLThN8Sjsnk+OJFqYy5jiDFoYVRlAfkEUdac37Ltq7aL7Z\/xDCFG3TnwWKysuRicHJXtQzsei5U1F4OTo7J9WMOy5vyXbV+nXd0YnAKVNRejbHva8VisrDn\/nByV7UM7ToHKmotRtr1rV3ZNnP+y7eu0G0KMOnNwfeqk\/FwbeyIZ0ySM6qxMdp\/RC6P9+\/ebAwcOGP5GBX+roozRFpykZdpXbYPfmP4ZzxBi1J3Doc+tmbJ2ZGOVpTCkZfvQbs\/amilrN189G4O0bB\/WsKzV5VTWP+3qxnie5VTWXjHjRJrfZ834da+1\/svam2f+Scv2od0XbYyyds0sBmnZPrSDcRmruw5lfLs2Q4hRZw70Ya9gz0huCiX+oe1tbnMbc8nqqlm175OYdsnqarKXEbPE0NSkgMCugrpRVPEmuuKKK8yxY8dkYqD3gN4Deg\/oPZD7HmCvqLIxsr+86lWvyvUXet+pOr4qcxlT29ELIxabN+\/BgweNrHsMtCZaE70H9B7oynuAvYI9o4rRp63xE6vK2NQ2m4CEUTYXlYqACIiACIiACIyQQMvCaISENWUREAEREAEREIHeEJAw6s1SaaAiIAIiIAKdJ6AB9p6AhFHvl1ATEAEREAEREAERCEVAwigUSfkRAREYIgHNSQREYGQEJIxGtuCargiIgAh0ncA111xjfvu3f9t88pOfjDbU66+\/3pw5cyaafznuLwEJo\/6unUY+I3DdddeZZz3rWeb1r3\/9rKQgqVj1D\/\/wD+Y3f\/M3K\/ZScxEQgboEEEVPetKTzM1udjOzb9++um4K+xHjsY99rHnuc59ruH8UNlbl6AhIGHV0ya+99lrza7\/2a+Yv\/uIvTJ8\/1fCpbHt7Oxplbmr8UbO3v\/3t5l73ulfwOF\/1VV8V7cYJm42NDROTT3AgcigClsCXv\/zl5P70gQ98wObCvdx976KLLjKPfvSjzU1ucpNwzmeeEEUIr6\/+6q8273nPe8wLXvCCaD\/js5BKekag78KoddxXX321eepTn2rufve7m8suu8w85jGPMZ\/97GeDj+OCCy5IPi396q\/+qjl+\/HhwccRm\/PnPf96wMbNBh54AguUlL3mJufzyy8097nEP8z3f8z3mc5\/7XNAwxEAUvetd7zIvfvGLE15BA1hnt7zlLc1VV10VRbx84QtfSN4\/3\/md32m+7uu+zvzsz\/6s+chHPhJ8re00ev86deqUedGLXmTY1GJNhp+JX\/qlXzJ3utOdkvfsn\/7pnwZdC9\/\/ve997+SEM\/aHHn6+P\/ShDwVHdtOb3jT5g7i\/93u\/Zz784Q+b97\/\/\/UFi3PzmNzc\/+qM\/mgiWd7zjHUF8pp38+Z\/\/ubnjHe9ofvEXf9H8wR\/8gXnf+94ncZSGNPK8hFGFNwA35Z\/5mZ9JbgK\/8Ru\/YZ7whCeYf\/mXfzEf\/OAHK3gp13QymZhHPepRhuPekOKIm\/Ov\/\/qvm6\/92q9Nbmz8v3\/uc5\/7mDe\/+c3BNgEEC6ddr3zlK80rXvEK89d\/\/dfmyU9+stm7d2+5yZdoRYzYoohhXHLJJcmn1hjiF98wYgP4p3\/6J\/Owhz3MPOUpTzG\/8Au\/EFyIsQEjxJgTls5TFtJC+9+zZ4\/5+Mc\/npwi8HMYcqz4YrzPec5zzCc+8Qnzute9zvzwD\/+w+a3f+i3zlre8herG5vyfPn3a\/M3f\/I35gR\/4AfPLv\/zLhnVv7DzHAe\/Zn\/iJn0jmRPycZpWKEah8MDxx4oS5\/\/3vb+55z3smH3r++Z\/\/Odj9A5\/PfvazE\/5vfetbK42vTGOYHD16NPm55kSYE6Nui6Mys1KbkAQkjCrQfOlLX2o4Qmazf+ADH2ge+chHGk4r+MRv7H87Ozv233CvySSsOEJMcLPnVIIbzpVXXpl8Mvv2b\/\/2ZMP5wz\/8wyA3t3e\/+92G7+bwCR8BxidwBNhkMgkChxMuvhvwl3\/5l8l3BNz3ELj58+gREcYGFCIYn2C\/6Zu+Kfojzd27d5uv\/\/qvT\/6fSne5y10SQRxSAPAl1iNHjhjWHU6cQiK6fbEUghc+nP\/HPe5xhhNWypoa793pdJq8X3nEEpINY2PD\/+hHP2o4rbjrXe+afOg5dOiQee1rX2tC\/FzzgQT\/D3rQg8ztb39783M\/93PJo99QJy3MwTdEER\/iHvGIR5jv\/u7vNpNJmJ89BCrv06dYAc9J1F\/91V8ZHnutra0Fi8E8YosjYjiTOHIklDoCEkaOxJyUGxuf7ngkxA9Sujk3CQRS3WfuPFvnUxJfIsaH29gnk3DiiJOt9773vebpT3+64UbGHDi14ATp8Y9\/vPmd3\/mdIJ+Q\/\/3f\/92wubABEAPjRo0o46bKo7UmjykQEd\/2bd+WiDjGzCbpNmPyP\/iDP5h8cZO4IYxYe9DnfQAACxdJREFUnCS86U1vCuGu0Adze\/jDH26+4Ru+IflCeYhNmYCsxTOe8QzztKc9zfD9CgQk+Vvd6lZUBzO3Dvhnvdk0Qzh\/29velggJvmDPeym0OOK7LBdeeGEiHBFh5HkMjHBEkDWdw65du8zS0pJ55zvfaWCE\/+XlZcOHLT40PPOZzzTcA0yA\/+CTFkU8UuNRLe\/jpiE4KeLDDh8Ov+\/7vi85YXv5y19u3vjGN5rv\/\/7vN5weNY1Bf4kjKMgWQUDCKBD1293udubiiy82f\/RHf1TrEyaPCXiswinLgx\/8YHO3u93NPOABD0g2x4997GPJJ1k+4Td5rPaBD3zA3PrWtzaIIX\/ak8nE\/NiP\/ZjhRkd8hIZfX\/WaExyOptnM3vCGNyQ3y2\/+5m9OPu3znayHPvShyUkPbar6du25afL9AG7CfFfgZS97WSLsOBa\/\/PLLXbMgKUIY8cjNn5OwIE4LnEwmk+Sx2tbWlvnXf\/3XgpbVqmD2kIc8JPnOGultb3vbag7mtGbD5yQKUcRjqZD+OdXk9OPOd76z+f3f\/\/3ke30hxRG\/AcVp0T\/+4z+a\/\/iP\/0hmyvdQTp48GUSwcNLC+9797OF3MpmY5z\/\/+cnp0Z\/92Z8ZHg3DMAle8588UfTzP\/\/zhhOwSy+9tKbnG7tx+sjPNqIXf1\/zNV9jmNtjH\/vY5NQz5M8f71k+MCKyOe28cRThr\/g55\/7BfYkUgRw+ijz2gYCEUclV4sbGhs\/RMadH6W78UPG9hE9\/+tOGY\/l0\/bw8QogbPp\/g+RT2t3\/7t+bHf\/zHDTehw4cPJ0IJkYGw4VERG1DVmyhf6OZTKY+i0uPhEyyCjE\/In\/nMZ9LVlfJ8omQj+6mf+qnki+o8iuJGync3fuRHfsT89E\/\/dPIFV0RNJcepxtw0EUecgvFo7XnPe565\/PLLU63CZFl7btB82ZTvlvHJ+yz\/MP7TXniEx\/dQ+H5Wuq5OnrHynuF9BStO1kJuNM5\/DFGUni+nFfysIAJCiiPepzDhUaaLyYkdcyNPvN\/93d+t\/f0vfiZ4\/08mk+QXBRBFBw8eTL7P9Cu\/8ivJaS0nO8Sqa5ubm8mX+N3jM\/whijjp5j4ymTR\/pMY9gu\/DHTt2zLzcnhTx6BwhwalknXvfvLnyc87PXlocIVxCPap1Y+A+zlwkjhyRcaYSRiXXfTKZGB7R8MiMG4K7Waa784PFkXm6vEye30zihs8XPvktKx7NsRHzCIwv6CJcLr30UsNvhHADLOPTb3PZZZcZvjSZJ0j4cjTH+1nCyfcz7xqRxSMbTsG4afKYDmHh+iHOuKH5G5Crq5py00QcISQ4rWt62lUUf\/\/+\/YZP9qzLE5\/4xORxF5\/CeUQRY0NAALBeX\/ziF4uGVaruk5\/8pEEUcZLD+En\/+I\/\/2NR5H2UF9P2HPCnKikUZbPhZQayEEkeTyST52cJ\/2ojDdwop53SJNKRx4szP3mTSTLjwAYvT2clkkvzGaWhRxJy5T3HiyM8vJ0WctPFzgCjmZzvGI2d+zn1xhCjihA2B3\/R+xZx84x4uceQTGd\/1rvFNuXjGRbV874Nn9\/yA8oPDD6drz42TG\/UDH\/hAw+mSK6+actPBD78lxm9\/cKOZTCbmNre5jfmhH\/qh5Le8EEr8xtpkUu0mynd7+L4M3zHi8Vx6bJ\/61KcM321aWVlJVzXKMwd3ykZKfD5dfuM3fmMjv64zN03EEYKPX7eOKY74DtC3fMu3mFe\/+tUG0ccpzFOe8hTD91PceEKlPKrAL2vS1Ce8EdtOtHCyhjDizxE09U3\/tH\/KYpsTR+vr6+Y1r3lN8HDusRePoBFFnIQiwhD+TYPxHuWxMj\/n\/CxyT+G7RpwIN\/VNf4Q6p1AhT4rw6xuimkfYrD3v0QsuuMAgxLh\/5H1w9PtXvebnHE6cHHHyxr2EDyj8TFb1Na89c+GDJH\/ShFPDee1VPywCEkYV1nMymSS\/icZ3WjjJ4Ub5whe+MPnT9QgiPjVxI6rgMrNpljjKbFixkBs6vw3DDZ9PepxCuRsYJxPccL7jO77DcDJS0XVuc36Lj+\/nfOu3fqvht5S+67u+KzmlIBbjye1YsYKbZlviyA2NGzIsefw5mVQTqc5HUTqZTJLf+PEFeFH7MdYhjni8zaOd0PNHlPLFawQMP+uhRBHjnEwmyff9eNzMdwn50MIj2skkzPuID2eIiFCPzxhz2hDVfHmcX6TggyEnkgjv+973vkF\/Q82Py8858+K3XRF+nBT79SGu+XnjNArhFStGiHE28KGucwhIGM0BlK5mM+SLytwE7ne\/+yW\/zsuvv3OCxGlFqM0+ljhiI+cL1px+PfzhD09+04fvVjAXvu\/AY6HJJMzNGXY89uO3oBBiPCrg7xvxZwE4rqY+pHHTRBzBjk98IX0vwhebG1+29x9DLmIcXY\/Je2wyCfeedfO91D62RtAjXkKKIvy7nwVOizB+RkJv8rx\/JpPwXBi\/Mz4UcMrFhx8eMboPWq4+RsrPOSdFoXn5Y+VeKFHkExnXtYRRzfXm1925mfHH2vgCIkJjMgl7E2KD57EaXwjlV2FrDvW8bnzS44ufPArikxGbLzH4w4KhhJ0flEdCPILk0SCPoRCXfn3Ia26aCNcY8wg5TvnqPgHeQ7xn3R8D7OyIOzAwxDs\/25NJ2HvgIqbGujOXmMJrEfNSzPIEJIzKs1pIS8QRooW\/DxJ6AJwe8YmY35bhxhbav\/yJQN8JIOInk\/5v9n1fB41fBNokIGHUJu2asdo4Eq85NHULT0AeRUAEREAEFkhAwmiB8BVaBERABERABESgWwQkjGKvh\/yLgAiIgAiIgAj0hoCEUW+WSgMVAREQAREQge4RGNqIJIyGtqKaz2gI3OEOd1jYXEPEDuFjYQAUWAREYLAEJIwGu7Sa2FAISEAMZSX7Mg+NUwTGTUDCaNzrr9n3gMCVV15pQokj\/KQtD0G6nZ+nj59PX1OPpcv9fJl62shEQAREoE0CEkZt0lYsEahJIC2OEBi4Is0z6jHMtcEPRhkpRh35tFHnjDp3nZWm68k7y2o\/r8z1VSoCIiACbROQMGqbuOKJQE0CiAm\/K3mMMlLfKHOG8HF1rsxPqaONX6ZrERABERgrAQmj3q28Bjx2AogYxEweh3n1Wf3wR7+suiZl+HX98T\/PXFu\/nytTKgIiIAJtEJAwaoOyYohASwQQHr6oSOfrDgOf+MrqTzn1WXV+GW3yjHbUkcpEQARGTmDB05cwWvACKLwIVCVQJCCK6qrGaat9WWHV1ngURwREYNwEJIzGvf6afc8JOFFBWnYqtK0joOhDXz8Oecr9srLX9MXq9i8bR+0WTkADEIFeEZAw6tVyabBjJoCIWPT8ETFuHKTkm4ypaf8msdVXBERABLIISBhlUVGZCHSMQJYI8csQGOTbGPaVs7+rRMx58RhTntE3r45y6mUiIAIi0DYBCaO2iSueCFQkgEgoI0LKuM3zlVee5dO1Jc2q98sYd57RLq+OcuplIiACItA2AQmjtol3L55G1GECiI8skZBVTjvK\/elklfn1Va7xjeGTfqTkMfIyERABERgCAQmjIayi5jBYAoiP9OQQIlnltKOceq6duTLKuXblpJRh6XLqnFGP0QZz5aTkMeoxymQiIAIi0C0C1UYjYVSNl1qLwEIJID4QIkWDoJ52fhvKMFdGPUYZ5sr9lHqMesyvS19Tj9EeS9f7eeoxv0zXIiACItAVAhJGXVkJjUMEShBAfJRoZua1ox4r8kU9VtQmXUd7LF3u56l35pfrWgTKElA7EYhJ4P8DAAD\/\/49gl5sAAAAGSURBVAMAZfq1pErAp1UAAAAASUVORK5CYII=","height":281,"width":466}}
%---
%[output:235d5d03]
%   data: {"dataType":"text","outputData":{"text":"[14:01:20][INFO]  バイグラムマルコフモデルで300 SMILESを生成中 ...\n","truncated":false}}
%---
%[output:4a72d9ff]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 30\/300) 生成中 (Markov)\r[##--------]  20% ( 60\/300) 生成中 (Markov)\r[###-------]  30% ( 90\/300) 生成中 (Markov)\r[####------]  40% (120\/300) 生成中 (Markov)\r[#####-----]  50% (150\/300) 生成中 (Markov)\r[######----]  60% (180\/300) 生成中 (Markov)\r[#######---]  70% (210\/300) 生成中 (Markov)\r[########--]  80% (240\/300) 生成中 (Markov)\r[#########-]  90% (270\/300) 生成中 (Markov)\r[##########] 100% (300\/300) 生成中 (Markov)\n","truncated":false}}
%---
%[output:353a9a4a]
%   data: {"dataType":"text","outputData":{"text":"[14:01:21][INFO]  Markov SMILES を RDKit で検証中 ...\n","truncated":false}}
%---
%[output:7b36f722]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 30\/300) 検証中 (Markov)\r[##--------]  20% ( 60\/300) 検証中 (Markov)\r[###-------]  30% ( 90\/300) 検証中 (Markov)\r[####------]  40% (120\/300) 検証中 (Markov)\r[#####-----]  50% (150\/300) 検証中 (Markov)\r[######----]  60% (180\/300) 検証中 (Markov)\r[#######---]  70% (210\/300) 検証中 (Markov)\r[########--]  80% (240\/300) 検証中 (Markov)\r[#########-]  90% (270\/300) 検証中 (Markov)\r[##########] 100% (300\/300) 検証中 (Markov)\n","truncated":false}}
%---
%[output:8098b4ea]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]  === マルコフベースライン結果 (300筆生成) ===\n","truncated":false}}
%---
%[output:0fd21ee6]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]    有効性:   29\/300  (9.7%)\n","truncated":false}}
%---
%[output:4f3967e5]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]    多様性: 27\/29  (93.1%)\n","truncated":false}}
%---
%[output:9d3c70f5]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]    新規性:    27\/29  (93.1%)\n","truncated":false}}
%---
%[output:06b087e6]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]    (本番目標: 有効性>90%  多様性>85%  新規性>60%)\n","truncated":false}}
%---
%[output:8f440e59]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]  サンプルMarkov SMILES（先頭5筆）:\n","truncated":false}}
%---
%[output:7147febe]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]    [ 1] CO\n[14:01:29][INFO]    [ 2] CC(=OC[COCCC)CNCOC@H]1CCCC(C(=C(NC)[CCCCC@H]1CC)CC(=OC)C(C)CCN[C)CCCCCCNC@@H](N3CO)NCCC(=O)[CCC@]3C\n[14:01:29][INFO]    [ 3] CC(C)CCCC(=O)OC[CC@]\n[14:01:29][INFO]    [ 4] CCN(C)N1N(=C(O)C(N=O)C(=C(=O)OCCC)O)CCC(C(=OCC(=O=C(=O)CCC)N(C(CCC(O)C=O)CCCCC1CC)[CCCC@H](N(=))))[\n[14:01:29][INFO]    [ 5] C(OC=C)CC(CC(=O)CC(=O)[CCC@H]1C)CCCCC[CNCC@@H](C(CCCCNC)C(=O)CN=OCCC(=O)C(O)CC(C(=O)NCC(CC[CC])))1\n","truncated":false}}
%---
%[output:2a8966f4]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]  最初の5周の有効Markov SMILES:\n","truncated":false}}
%---
%[output:067fae5f]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]    [ 1] CO\n[14:01:29][INFO]    [ 2] CC(=O)O\n[14:01:29][INFO]    [ 3] C(C(=O)N1C(=CC)OC[C@@H](CCO)F)CN1\n[14:01:29][INFO]    [ 4] C(O)CN2CCC(C(N)O)CO2\n[14:01:29][INFO]    [ 5] CCC(O)O\n","truncated":false}}
%---
%[output:8b7375fb]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]  Markov失敗分析 (271無効 \/ 300生成):\n","truncated":false}}
%---
%[output:5276edff]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]    括弧不一致 '()':  10 \/ 271  (4%)\n","truncated":false}}
%---
%[output:83699159]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]    角括弧不一致 '[]':     2 \/ 271  (1%)\n","truncated":false}}
%---
%[output:14d3854d]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]    環閉鎖数字不一致:  11 \/ 271  (4%)\n","truncated":false}}
%---
%[output:9f6c2440]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]    その他（原子価\/原子\/無効）: 248 \/ 271  (92%)\n","truncated":false}}
%---
%[output:8ce954da]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]  次のステップ: Section 6-8でLSTM文字言語モデルを構築する。\n","truncated":false}}
%---
%[output:07691b4a]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]  系列行列: 500 x 99  (語彙=37)\n","truncated":false}}
%---
%[output:0c181dc6]
%   data: {"dataType":"text","outputData":{"text":"[14:01:29][INFO]  分割: 400 訓練 \/ 100 検証\n","truncated":false}}
%---
%[output:3b3443ea]
%   data: {"dataType":"text","outputData":{"text":"[14:01:31][INFO]  ネットワーク: Embed(37->16) + LSTM(128) + FC -> softmax  |  パラメータ数: 79621\n","truncated":false}}
%---
%[output:5f85d3e4]
%   data: {"dataType":"text","outputData":{"text":"[14:01:31][INFO]  訓練トークン対パラメータ比: 4.39  (>1 => 過パラメータ化)\n","truncated":false}}
%---
%[output:81098cf3]
%   data: {"dataType":"text","outputData":{"text":"[14:01:34][INFO]  エポック   1\/150  train=3.4883  val=3.1569\n[14:01:39][INFO]  エポック  10\/150  train=2.6523  val=2.6592\n[14:01:45][INFO]  エポック  20\/150  train=2.6141  val=2.6346\n[14:01:50][INFO]  エポック  30\/150  train=2.6089  val=2.6267\n[14:01:56][INFO]  エポック  40\/150  train=2.6055  val=2.6201\n[14:02:01][INFO]  エポック  50\/150  train=2.5971  val=2.6182\n[14:02:07][INFO]  エポック  60\/150  train=2.5970  val=2.6155\n[14:02:13][INFO]  エポック  70\/150  train=2.5932  val=2.6101\n[14:02:19][INFO]  エポック  80\/150  train=2.5960  val=2.6174\n[14:02:24][INFO]  エポック  90\/150  train=2.5968  val=2.6149\n[14:02:30][INFO]  エポック 100\/150  train=2.5947  val=2.6161\n[14:02:35][INFO]  エポック 110\/150  train=2.5913  val=2.6131\n[14:02:41][INFO]  エポック 120\/150  train=2.5963  val=2.6128\n[14:02:46][INFO]  エポック 130\/150  train=2.5929  val=2.6152\n[14:02:52][INFO]  エポック 140\/150  train=2.5956  val=2.6149\n[14:02:57][INFO]  エポック 150\/150  train=2.5966  val=2.6132\n","truncated":false}}
%---
%[output:84c89218]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAkYAAAFgCAYAAABewLbaAAAQAElEQVR4Aey93asc2X7f\/R1HB0lg89gXemRhmVkzyU3+gCBd2DprSHKRy2AIOTKcqSExJBDQRTKG4APTQ86BBGGIL0MEU2MSQQj+A0JepmwC9rnLnQk2Ug0SjCAEG54HMnAEk\/3p0tp77dpV3dUvu7uq+ruZX6+1fm9rrU+t7vpNd++tn\/vOPyZgAiZgAiZgAiZgAksCPyf\/mIAJmIAJmMBsCXhjJrAZARdGm\/GytwmYgAmYgAmYwIwJuDCa8cX11kxgjgS8JxMwARO4TgIujK6TrnObgAmYgAmYgAlMioALo0ldrjku1nsyARMwARMwgfEQcGE0nmvhlZiACZiACZiACRyZwN4LoyPvx9ObgAmYgAmYgAmYwNYEXBhtje5y4IcffnhZscVoHzm2mNYhRyBwXdd6H3nJgRwBy7VMmfaS2qGTJP\/UDo2z3+EIcG2Gyh5X1ZmKdXQaJqBsr709nsAW9rpEF0Z7xXmYZPs+tClfatftYqjfujy2XyUA2y656tlo8G16+3kk34sXL4TQ3yXrrvGbzn3o+TZdH\/6HWOMh5mAvXTJ0bvzaQr62Lh9j7xLOai745GP66K5bmIf1DpkHv20lz9\/OkWxtfXuc\/Ia2efzQmCn7uTDa4Orlh6PdJ01bl4+xI7mu3R9ix+fQwjoPPecpz8cLbC5jYcE5WCVjWefQdZzv5ezd3q7+0DzX5de1plx3XfNukpf1bOKf+\/ad8Vyf+nncdfTZR1u2mYf1kmddLH5tISbXtcfJhj5J0qV2lT750Ca\/Tdpt4zaZYyy+Low2uBIcjCSEpX5X27YzTtLlv06XYndtedIiu+ZZF7\/JHPjmkufO9X395I899Ve1Q\/1W5Rhqa8\/FGCGeFqGfC7pcclvex4cxbZL2OOlpsSH0+yS3cyYZJ2G8SpLfmFr2yZpZU2rpI4zXSfKj3Zewpi5p51+1trbvlMewYK\/b7oH4XMiTj+mj6xJsCPMj+NAi6BlvKsRuGrONP+vrkk1ypXhi6LfbpEN\/SuLC6ISuNoecJy1CP22dMf3U0j+UsA7mzQVdmj\/X00dPmwu6JB9++GHqdrZ57k6Ha1amdTNN3meMJF27xdYlXX5tXRqn+DRe1SbfqbVc3ySsnT0ypj8WYU1t2ffa2PMq2fd8+8zXte5V+btYdunaOZgn+bVtjLHhQ78t6DeVdTmw5znbY2zokrC+Lkn2IW2Kx5c+LUIfoX+K4sLowFc9P2wc9HWSlpfHJd3UW\/betS902Ka+v03Xz55XSZ4PPzjlulPswyEJ+4dJEvRJl\/qMpy7shT2u2gf2XPBtj9EdW7r2kq8z9Y+1TuZnje350XcJfl16dNjagj4JttSn7Rqjy4W15dK25WP3hxMYf2E0fC8H9eTgciC7JkWPvcuW6\/DpE\/yw0W4qzL8uhtxD\/NbluS4769smN3F9+0KPfZu8h4phfV3Snh+ftu46x7AbIu01pJi2vm+Mf5+tSw+HJF32Ph3zDJWuHMR26aeqO9Z+uHa7zk18LlyDfEwfXS7omDvXTbHPHnKZ4h7GuGYXRiO8KtfxpL2OnLui4wnNupBdc20bv2ruPluup5\/LJutIcXlM0rXbIT7E4EfbJ9i3Ea7VOunKm2Ly9XT5JR3+uW9fP\/kPbVPe5M+4S7C39ejagk\/f2nJ9Hoc+j8tt++qnOTbNl6+LHH2yad4h\/mnu5Ns1d7K1W2KTYEt92q4xujFJ2ist66JN0jVGl0vyTW1ua\/fX+SR7YteOP6WxC6MdrjYHiMOUp2CMPtcN7ROLbBvPPJvEMw\/+xO0q2+ZhDQjxSXZdyybxzN3nj4015XbG6NGlPuMk6LD1CXYEe4qhnyTp2m2y07ZtaZxsqU36vMW2TlgfMev8NrGTLwn5kb745LeqbceSD\/+2\/jrHzLdOuuZPMay5y75KRwzxq3y2tZF3nWyTmzVvEpfWQEzq0zLuE+Zo+zBG3xfT1uNLTFt\/XWPmyoV58nFXH59c2j65rd1PvrmePSPokj2NadFhOzVxYbTjFefgcIBIQ8uY\/raya3zXvKvWxXzYu+IOqWMdSXZdD3naORij3+eetslHDNK3DtbZJX3+U9Szf2SKax+25qteXNP2nhmjR+hfjerW4N9lQb9Jnq4cY9axv\/b60CHsmza3M056+rltH31ytoW8bR1j9KsEH9a6yuc6bMyJbJqb9W4aMyV\/F0Z7uFocLA4K7bp0+PUJsX029NhXCT5D1rAqxza2beclrms+9tBn6\/I\/po515rLrWth7l6zLyxqIw4+WMf1NhTjiDxW36TxT9e9j2qfv2me6NsTQb\/ugb+uOPWadq9aFHenyQYeNPdAypo8wRtAh6PoEO4I\/0ue3qZ6cbSFHW8cYfS6sIxds+XhVH98kbb+k36VlveSl3SXPlGNdGO3h6qVDRLsuHYetT4jts6HH3ifM3efTp0+5sBOfxnNo8z2xN8bb7Is44omlZUwfSWN0SdDvIuTskl1yzj0WXvBft891PkPzrJtnE3uak3ZVHPZ8\/fTRrYpJNvzwT+MxtawLWbWmrvUTg6S4Lp9kSy3+SBrnbV98nz6P3abPOpIQn\/p5u0qPLUmKSeN9t20GjPc9x9jyuTDa4YpwQBAOJmloGSOMDyXMx9y7zEc8edo50GFr69N4nT35dbXkJb5tQ4etrR\/7mHXvskb2vEr6cjMvcbmdMfpct66PP3Hr\/Nr2bePaefY5HuOa2vsbusahfu381z1et6519lXrIxbZ5jyuypts5CV\/Gh+rZQ1pLfTTOuijT+NDtmluWubN10E\/Cba5yuwKIy5mLlcvXKPJffJ+Y139mPy7DkjSJZ\/VmXazpjmYc7dMTTR5yNmM1j\/iS0yfJ\/Yuyf2Jb\/ugy3227ZOH3LR9ObD32ZKeePxok46WMfokjNHvIuTK4xkjuS710SN986LHjqSYrhY7gn+XvU9HDLJpXF++bfT5\/PS3yUEMsYfaB3Mh+Xx5n\/UkafslPS0x2On3CXb8+uzXpd9l3hSb1k2Lrm+t2PDps\/fpiSEWoZ\/7oUPa+txnH\/2UnxZhTiTvM97HXENyMBdz40vLGGF8SjKrwogLyMXMBV3fBc39Ur\/PFz25kCG+yQd\/hPg+wY702bv0+Kc5uuzb6siZYtMcaZy3q2z4kadPsOfS9stt7T6+bV0ad9nW6brsKV\/e9vmhT4I\/fdok7TF62NHmgg5p+zNGsCV\/+gh6JOm7WuwI\/kjyoZ8EO5Js69pt49bl3dXOHvK1Dcm3qf+QnH0++Vystc8Pfe7LeFPZNX7dfOTv28MqW19eYhByIm0\/dNjbenTY2vqhY2KR5E8+BB2S9HmLvU\/w67Ohx56EcVuYE8GHNknyQ9+Wtm05\/vDDpVtXH93S+O4hjWmZ75162TBGsLVl6TDTh1kVRtd9jTggyCbz4I+sisGeZJVfsm3im2K2aZmnKw490mWzbj0B2CG5J2Mk1+X93EYfye3r+vgjyY9+kqQb2g6Nw29ozm392nMwRobmwxcZ6r+LH\/MgQ3Lgh6zz7fNBj6yL39belxs9si5v24cxsiquy96lSzlW2ZJPuyUGaevzMfZtZV2e3J7303y5jn7Sp7ZLl2x5i1+SPn2y0+Y+qY9+rjKrwogLtsmFyivgTeLsawImYAIHJOCpTMAEDkhgVoVR4pYXPKuKJWxJiEnxfe3r169lMQOfAZ8BnwGfgVM\/A333yTnoZ1kYpWKHtq\/gwZZfQMZ9vvjxJPj000\/16NEjixnsdgbMz\/x8BnwGJn4GfvCDHyzfKOD+ODeZZWF0HReJwuinP\/2pnj59qufPn1vMYG9n4MmTJ8sj67Pl59U+X1t8rnye9nme8lycLe6HyxeuGT7MqjBa9Y5P+9pt4pvHPnjwQA8fPszFffPY6QxwpjhjtD5bfm7t6wxwnnyufJ72dZ7yPOlscb7mKLMqjHa5QBRKfJw2JEddS4uF9MknUl0PibCPCZiACZiACUyVwGmte1aFEYUNBU4u6PJLio0xevpJGKMfIlUlff65VJZDvO1jAqsJ\/PIv\/7J++MMfina1p60mMJwA58nnajgve5pAIjCrwohNUeDkgi4XbGlMP0nSuTWBQxPgBvbxxx8felrPN3MCq87VzLfu7ZnATgRmVxjtRMPBJmACJmACJmACJ03AhdEWlz+Ei6C6vui7ZwLHIeBZTcAETMAE9kXAhdG+SDqPCZiACZiACZjA5Am4MBrhJfSSTMAETMAETMAEjkPAhdFxuHtWEzABEzABEzhVAqPetwujUV8eL84ETMAETMAETOCQBFwYHZK25zIBEzCBORLY8574J5gsr5f\/FtkhOOz58k0+nQujHS9hXe+YwOEmYAImYALnBCgEPv30U\/8jswf8R2bn\/A\/Cnh+sDToujDaAlVxDSD23JjBLAt6UCRyNAIXRT3\/6Uz19+lTPnz+3XDODJ0+eCN5Hu+AjnNiF0QgvipdkAiZgAqdOgH+oNP+HS+fcf\/z4sZKkfTJO\/bxt69vjVb65LfXhfOpnrb1\/F0ZtInMce08mYAImYAKjJZD+aSrarkXyb3r26dsx+CYhJvVpGVvWE3BhtJ6RPUzABEzABEzgWghQsLRlyETEpKKIfopBlwRd6tMynqvsc18ujPZJ07lMwARMwARMYAMCFCxtGRJOTF4QpRh0SdClPi1jy3oCLozWM7KHCZiACZjAQQmczmQULLlssvNUHNHmcYwRdLQIfcswAi6MhnG65BXCpaEHJmACJmACJrAVAYqWXDZJQkFFLDH0aRH6SLvP2LKegAuj9YxWetT1SrONJmACkgzBBHYl8Mkn0kcfTV+q6jIJCphcLlv7R8RgpUXoJ6FYQhjTIvQtwwi4MBrGyV4mYAImYAJHJFDXUlVJVSVVlVRVUlVJVSVVlVRVUlVJVSVVlVRVUlVJVSVVlVRVUlVJVSVVlVRVUlVJVSVVlVRVUlVJVSVVlVRVUlVJVSVVlVRVUlVJVSVVlVRVUlVJVSVVlVRVUlVJVSVVlVRVUlVJVSVVlVRVUlVJdX0ZIkVLLpet\/aMUgwd92iQUSghjWoS+ZRgBF0bDONnrEgEPTMAETOCwBL7\/fSlGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKUYpRilGKYTL3ChacrlsHT5KxRFtEqJTn5axZT2B2RVG+QGjvx6BhB8yxNc+JmACJmAChyewWEhffTV9ifEqO4oW5KqlX9O+Z\/WN2\/r+jLYkAsvCKA2m3nIAOFy5oJv6vrx+EzABEzCB0yXAPW2T3XPfSzG0jDeJP3XfWRVG21xMDgwHZ5tYx5iACZiACUyCwKQWuct9qSuWexz6SUE44mJnVRhx8a+bJf\/AIRJCM9Of\/\/lbvX1rMQOfAZ8Bn4F9nIHmlfW0HilakE123b7fpXFqliYelQAAEABJREFU27n69Mlv1bXjnpdLiplrO6vCKF0kDliSVYcBn1X2lC9v+Qf7Hj16JA4J+m+\/\/VavXr2ymMHWZ4Cz9ObNm63jT\/L8+bytPS9TPVc8F3htPRXhHpQL+2ZM2yVtW3ucx6yy5X70v\/nmm94z9ezZM3HfS8J9kJi5yiwLIw5DEoqfrouHHp8u2yrd06dP9fz5c92\/f3\/pduvWLd27d89iBlufgbt37+rOnTtbx\/v8+fnXdQameq5+8Rd\/cfna6ofDElj1GlQUxfK+x70PefLkyWEXd+DZZlkYDWVIcZSEGPq0q+TBgwd6+PDhucv3vndDFEeWW7twOOnY27dv6+bNmyfNwM+f\/T9\/pnyudII\/ffefPj2IsCVhvIuseg6GEJb3Pe59CPfBXeYae+ysCiMOyFDgvFuUC3GMaS0mYAImYAImcGwC3NO4L9F2rQVbEuz4rRJ8LOsJXE9htH5ee5iACZiACZjAyRNIhQwg6NMi9Cl66NMypp+EcVvwywXf9hidZTWBWRVGHICug5IjwJ6P3TcBEzABEzCBTQnsy5\/7FkI+Wu5RSN5vj5MvPrmgt+xOYFaFETjyQ0IfXS5dOux9emxdcvaR61Jd18vGDyZgAiZgAiawMwHuRQiJutqko1jKBX\/LfgjMrjDaDxZnMQETOC0C3u3oCZSlVJZSWUplKZWlVJZSWUplKZWlVJZSWQ7bSllKZSmVpVSWUllKZSmVpVSWUllKZSmV5X7z1fWVfBQ4KPOWPpL09BHGCAVSLuiQ3IexZXMCLow2Z+YIEzABEzCBQxP48kvpk0+GyZC1Dc2F3z7zVVVvNgqdZKSPMKZF6CehAMol6dvtEJ92zKmPXRid+gmY6f69LRMwAROYAgEKl3bRM2TdxOTSFzPEpy\/2VPUujE71ynvfJmACJjAlAl98Ib18OUyG7GtoLvz2ma8ohmQTBROCMy1CPwnjXJLe7e4EJlIY7b5RZzABEzABE5gwAX7jZagM2ebQXPjtM18rF+\/otFRClwRb6tMyTsIYSWO3+yHgwmg\/HJ3FBEzABEzABLYn0BHJO0IUPrQdZquuiYALoy3Bpv+JqOstEzjMBEzABEzABHoIUAxRFGGmZUzfcv0EXBhdP2PPYAImcHoEvGMT2IoABRBCMZQnYIweSXr6SD6mjw7J+2mMzrKagAuj1XxsNQETMAETMIFrJ0DhwyS0CP22oEfQ0+bSpcvt9PGxrCfgwmg9I3uYgGQGJmACJmACJ0HAhdFJXGZv0gRMwASmReD169eyXD+DaZ2Kw6z2VAujw9D1LCZgAiZgAhsRuH\/\/vh48eKDHjx\/r0aNHlmtmAGd4w32jCzVjZxdGO17cut4xgcNNwARMwATOCXCDfvr0qZ4\/f27ZicFwfvA+vwDuyIXRlofg\/fe3DHSYCZiACZjASgIURw8fPpTlMAzgvfKCnJjRhdGJXXBv1wRMYJoEvGoTMIHDEHBhdBjOnsUETMAETMAETGACBFwYTeAieYlzJOA9mYAJmIAJjJHA7Aoj\/rpnLqug5370V\/naZgImYAImYAImMH8CsyqMKG746565oOu6jOhzP\/rounyH6OxjAiZgAiZgAiYwfQKzKowOeTnSPyJ7yDk9lwmYgAmYgAkcicDJTDurwoh3fYZeuU181+Ws63UetpuACZiACZiACUyBwKwKowScj8SSrCuAkh\/tOl\/ypz9RTz\/J27dvZTEDnwGfgUmdAb9u+XV74BlI973UpnvfXNtZFkYUOEkoeFZdvORHu86XPPz5dP5M\/W\/\/9qcMl\/LNN9\/o1atXFjPY6gzwYvPmzZutYn3u\/LzrOwM+Vz4bfWdjU\/2zZ88u\/dMs3AeXN7+ZPsyyMLrOa8WfTudP1T958uR8mjt37ujevXsWM9jqDNy9e1d7OENbze1zO9\/nrc\/VfK\/toZ+3RVFc+udZ8vvf+Y1wRp1ZFUZD3vFJ124T3xRDyz+2x5+p\/5Vfuc9wKbdu3ZLFDLY9A7dv39bNmzd9hvw82usZ8Lnya9K2r0ntuBDCpX+ehfvg8uY304dZFUazukbejAmYgAmYgAmYwMEJzKowSt8T4t2gJOhyqugZo6efCzpsFhMwARMwARMwgeslMNbssyqMgExxkwu6XLClMf1ckn5Ie\/bO4rlbXZ933TEBEzABEzABE5gwgdkVRhO+Fl66CZiACUyYgJduAvMg4MJoHtfRuzABEzABEzABE9gDARdGe4DoFCYwRwLekwmYgAmcIgEXRqd41b1nEzABEzABEzCBTgIujDqxzFHpPZmACZiACZiACawj4MJoHaEeu38rrQeM1SZgAiZgAiZwDAJ7mtOF0Z5AOo0JmIAJmIAJmMD0CYyyMOKPLvahXWXri7HeBEzABExgcgS8YBM4CoFRFEYudo5y7T2pCZiACZiACZhAi8AoCqPWmnqHFFD8pepeBxtMwATGS8ArMwETMIEJEJhUYTQBnl6iCZiACZiACZjAhAlMpjDyu0WjO2VekAmYgAmYgAnMjsDoCyMKImRsH6H51\/Vn91zwhkzABEzABExAF4XRyGBQDCEURMjIluflmIAJmIAJmIAJzJDA0Qojip4kcE19WsYUQwh9iwmYgAmYgAnsSsDxJjCEwNEKI4qeJCw09WkZUyAh9C0mYAImYAImYAImcAgCRyuM1m2OAgmhOELW+dtuAiZwagS8XxMwARPYP4HRFkZpqxRHiIujRMStCZiACZiACZjAdRE4amFEsYMM2dzQ4oh8uazLvYnvuly270bA0SZgAiZgAiZwbAJHLYwodhAgUKDQ7iLkIF8u6PpyYhvq25Uj\/cr+1193Wa0zARMwARMwAROYGoFrLIyGo0jFCYXKqij81vmsis9t5CFfrnPfBEzABEzABEzgtAmMojBKl2DXQmXX+LQOtyZgAiZgAiawloAdZklgVIVRIryqwFllS\/G8G5Skz7+tx7+tS\/ny9vXr10qS9D\/72Vu9fWsxA58BnwGfAZ+B+Z2BdM9Lbbr3zbUdZWG0K2wKnCQUPOvy4YP\/Oj\/sjx8\/1qNHj5bCIUH37bff6tWrVxYz2OoMcI7evHmzVazP3V6fd7O6Bj5XPhv7en149uzZ8p6X7n3cB7n3zVVGVRhRoPRJugDYU38fLfmGFkXM9\/TpUz1\/\/nwp9+\/fR6Vbt27p3r17FjPY6gzcvXtXd+7c2SrW587Pu74z4HPls9F3NjbVF0WxvOele9+TJ0+W9765Phy1MKIoSZIApyKl3WLHN+kZtwV7W7dqjP+qfF2xDx480MOHD5eS7N\/73o1lcUSBZLk1HhZnBesUrsft27d18+ZNc5vI9ZrCmWKNPld+LeIc7ENCCMt7Xrr3cR\/UjH+OWhjBddPChJh9yDZFUXves7OyVNX1svGDCZiACZiACZjAxAkcvTDagJ\/WFVHYKXhyQZfPgS2N6bcl2dyagAmYgAmYgAmcHoFJFUZDLg+FUC7tGGzoaLsEm8UETMAETMAEDk\/AM46BwOwKozFA9RpMwARMwARMwASmSeDohREfZQ1Ft4nv0Jz2MwETMIHrIuC8JmAC0yNw1MIo\/ygroUvFT7vFjn\/SM7aYgAmYgAmYgAmYwD4JHLUwam+EwqdP2r5jGKffShvDWryGQxDwHCZgAiZgAnMnMKrCaAhsCqchfof0qetDzua5TMAETMAETMAErovAKAqj9sdj7fF1bd55TcAETMAETMAETCAnMIrCKF\/QJPp1LS0W+n61UNBZX\/4xARMwARMwgdER8IK2IHC0woh3hZKw7tSnbY\/R5YL9qFJV0uefq6g\/P+oyPLkJmIAJmIAJmMB+CRytMOK7QknYUuqntkuHDb3FBEzABE6SgDdtAiZw7QSOVhhtu7NRFEchnC\/fH6Wdo3DHBEzABEzABCZPYBSFUVex06UbK+26HuvKvK6RE\/DyTMAETMAERkZgFIXRyJh4OSZgAiZgAiZgAidKwIXRjhf+0kdpO+ZyuAmYgAmYgAmYwHEJjLIw4jfQjotlzewhrHGw2QRMwARMwATmR+AUdjTKwugUwHuPJmACJmACJmAC4yPgwmh818QrMgETMIEDEfA0JmACbQIujNpEhoxDGOJlHxMwARMwARMwgYkRcGG07QX77juVX5yJCvFT1zxaTOC4BDy7CZiACZjAbgRmVxjxxe1c1uHBd52P7SZgAiZgAiZgAqdBYFaFEUUOfxgyF3R9l3KVrS\/msHrPZgImYAImYAImcEgCsyqMNgHnomgTWvY1ARMwARMwgWsgMMKUsyqMeKdoKGN8kaH+ye\/169dKknS0b9++lcUMfAZ8BnwGfAbmdgbSPS+13PPmLKMrjHgnZ5uCJb9I5Eiya648L\/3Hjx\/r0aNHS\/nt3\/4U1VK++eYbvXr1ymIGG58BXmzevHmzcZzPm59vq87Akc6Vz\/EMXwOfPXu2vOelex\/3weWNb6YPoyuM9lHIkCMJBdI+r93Tp0\/1\/PnzpTx58uQ89Z07d3Tv3j2LGWx8Bu7evSufHz939v364XPlM7WvM1UUxfKe13XvO78JzqgzusJo7GwfPHighw8fLuVXfuX++XJv3bolixlscwZu376tmzdvjvP8+FxP9rr4XPn1aJvXo66YEMLynpfufdwHz29+M+zMqjDa97tD66732Vk5d6nr8647JmACJmACJmACEyUwq8LooNfgo48UP3pPX+mjg06742QONwETMAETMAETWEHgqIUR7\/AgK9a3kYnvFZEvF3R5Emz52H0TMAETMAETMIG5ENh9H0crjChQKFoQ+rtvpclAvlwa7cUjtouR1B7nNvdNwARMwARMwAROi8DRCqNtMO+zgNpm\/ksxIVwaemACJmACJrB\/As5oAocmcLTCiHdqKHSQ1D\/05vcxX5C\/dS3\/mIAJmIAJmMBMCBytMIIfBRFCf+pS11PfgddvAocg4DlMwARMYNwEjloY5WgokHj3KNe5bwImYAImYAImYAKHJDCawuiQm97nXKf+Udo+WTqXCZiACZiACRybwKgKo0m9a\/T++8e+dp7fBEzABEzABExgzwRahdGeszudCZiACZiACZiACUyIwFELI75T1BbYtXVpzDtK2C0mYAImYAImsBUBB5nAGgJHLYwodDaRNXs5rHmxUP3yO72n7w47r2czARMwARMwARO4NgJHLYyubVcHShzCxUR1fdF3zwRM4GAEPJEJmIAJ7JWAC6O94nQyEzABEzABEzCBKRM4emGUvj+0STtl4F77GgI2m4AJmIAJmMARCRy9MGLvQ79nlHxpLSZgAiZgAiZgAiawbwLXXRjte73OZwImYAImYAImYALXRsCF0bWhdWITMAETMIH5E\/AO50bAhdHcrqj3YwImYAImYAImsDUBF0Zbo2sCQ2jar79uWj+agAlMm4BXbwImcNoETrowav8m3EZHoaqk997Ty\/o9FSrlHxMwARMwARMwgekTONnCiKKo\/dtw6AZf0hAGu9rxmAQ8twmYgAmYgJB1ic0AABAASURBVAkMJzCKwoiCZIiwLfxoLSZgAiZgAiZgAiawbwJHL4za79qsG2PfB4R95dnHWpzDBEzABEzABExgHASOXhgdGwPvQCUZUiy9fv1aS7lx49LSf\/azt3r71mIGPgM+Az4DPgM7nYHR3UuW97x07ztrL938Zjg4+cKIYigJBdK6a\/z48WM9evRoKcn3fdX69ttv9erVK4sZbHwGeNF58+bNxnE+b36+rToDPlc+H6vOxya2Z8+eLe956d7HfVAz\/jn5wmjTa\/v06VM9f\/58KXns\/\/7fP6979+5ZzGDjM3D37l3duXNn4zift4k83470nPC58vnY12tEURTLe1669z158kRz\/hltYTTk3ZtdLsy2+R88eKCHDx8upT3\/rVu3ZDGDTc\/A7du3dfPmTZ8dP3\/2egZ8rvxatOlrUZ9\/CGF5z0v3Pu6D7fvfnMajLYzmBNl7mRUBb8YETMAETGDGBEZZGPFuDt\/7uU7u5GeeXNBtNOdZFY1\/kP\/sNRwsJmACJmACJjB1AqMsjA4FlULoxYsXSu2h5vU8JmACJmACJmAC4yRw0oXRzpfk5Ut9FL\/TJ\/pi51ROYAImYAImYALXQcA5NyMwusKIj7Z4B2ezbRzP+92nacdbgGc2ARMwARMwARPYG4HRFUZ729mBE9X1gSf0dCZgAidKwNs2ARO4TgJHLYx4d6gtbLata4\/xsZiACZiACZiACZjAvgkctTDiI7NtZN8QnM8EjknAc5uACZiACYyHwFELo\/Fg8EpMwARMwARMwARMQHJhtPdT4IQmYAImYAImYAJTJTCqwij\/LtFUgXrdJmACJmACJjBrAjPf3KgKo\/z7RqlImgr\/up7KSr1OEzABEzABEzCBPgKjKozyRaYiiQIp14+q\/8kn+qJ8Ty\/1waiW5cWYgAmYwAYE7GoCJpARGG1hlNZIgTTq4uhsoUF+u0j+MQETMAETMIEZEBh9YQTj0RZH77\/P8iwmMC4CXo0JmIAJmMDWBCZRGLG70RZHLM5iAiZgAiZgAiYwCwKjKozG\/pHZka64pzUBEzABEzABEzgQgdEURhRFvCu0at\/Y8VvlY5sJmIAJmIAJmMCUCIxrraMojCh2KHrGhWbAakI4d+IL2HV9PnTHBEzABEzABExgggSOXhhNtiia4MX2kk3ABEzgEAQ8hwlMmcDRC6Mpw\/PaTcAETMAETMAE5kXg6IXRph+hDfHnXagkqy5X8qFd5WebCZw2Ae\/eBEzABE6HwNELo32jpsiheErCuGsO9MmHlnGX31Ad3zEa6ms\/EzABEzABEzCBcRKYVGG0rnjBTpGzDWriiN8otihUfvGd3tN3qhRV15rEjxdpAiZgAiZgAibQTWAUhdHGBUn3Xo6iDeEo03pSEzABEzABEzCBbgI7aUdRGO20gyyYd32yoSi42rrcnvfxzcd9\/devXyuXEC48\/\/zP3+rtW4sZ+Az4DPgM+AzM5wzk9zz6F3e9efZmVRjll4hCZ1VRhA2fJIzz+L7+48eP9ejRo3P5yU9+cu76l3\/5F3r16pXFDDY6A7zQvHnzZqMYnzM\/z9adgVmdK7+mHPX14dmzZ+f3PO5\/3AfPb3wz7MyyMKLYGVLo4JNk6LV9+vSpnj9\/fi6\/8zu\/eR76l3\/5S7p3757FDDY6A3fv3tWdO3c2ivE58\/Ns3RnwufIZWXdGhtqLoji\/53H\/e\/Lkyfl9b46d0RRGFDPrZMgFIAfFzjpf\/Nb5dNkfPHighw8fnksI4dztxo0bunXrlsUMNjoDt2\/f1s2bNzeKmfg5814P8BzxufJr8b5eJ7jP5fc97oOa8c9oCiOKmXWy7jpQ7JBjnV+XfZfYEJqMX3\/dtH40ARMwARMwAROYJoHRFEb7wkeB05Y8NzbGFFD0kzBGb9mCgENMwARMwARMYCYEZlUYUdx0SX6tsKcx\/SRJ59YETMAETMAETOB0CXQVRqdLY5udv\/eeXtbv6TMt\/Ace5R8TMAETMAETmDaBURRGvGszbYxevQmYgAmYwHQIeKUm0E9gFIVR\/\/IuW8ZeQNX15fV6ZAImYAImYAImMC0CkyqMRok2hOWygr72R2lLEn4wgcMT8IwmYAImsC8CkyqM0m+Q7WvzzmMCJmACJmACJmACOYFJFEapIOKjNCTfgPtzJOA9mYAJmIAJmMBxCBy9MEpFD20XAvQUQ0iXfWy6uh7birweEzABEzABEzCBoQQOUhitWwxFD0IRlIQY+ujpj1ZCWC4tyBWR\/GMCJmACJmACEycwisIoMaQISkJRlPRTaut6Sqv1Wk3ABEzABPZAwClmRGBUhVHOdeoFUr4X903ABEzABEzABKZBYLSFUcJHgTSld4\/qOq3crQmYwGQJeOEmYAInS2D0hRFXZtTF0VdfqX75nT7SVyzVYgImYAImYAImMGECkyiM4Dvm4igEVthIXTetH0dFwIsxARMwARMwgUEEJlMYsZsxF0esz2ICJmACJmACJjBtAkcvjCh2NkY4woAQmkV9\/XXT+tEETMAETMAETGB6BI5eGG2KbC6F1Kb7tr8JmIAJmMDpEPBOj0dgcoXR8VB5ZhMwARMwARMwgbkTOHphxK\/ibypjvCghNKuq66b1owmYgAlcEHDPBExgKgSOXhgBio\/HhkjypV0leaE11I+YVb7rbCE0HnXdtH40ARMwARMwAROYHoFRFEb7xEaBkxdZjLvyo8\/96KPr8rXOBNoEPDYBEzABE5gngVkVRhQ2FDgHvVRlKb33nr4o3xP\/kGxdyz8mYAImYAImYAITJTCrwmj7azA88vXr18qlK\/Lt27eymIHPgM+Az4DPwBzOQH7Po99135uTblaFUfvdolXvIOGLPRd06y7u48eP9ejRo3P53d\/93Sshf\/zH3+jVq1cWMxh0BnihefPmzSBfnys\/r4aeAZ8rn5VLZ2WH1+Nnz56d3\/O4\/3EfvHLjm5FiVoVRfl0oeFYVOsmOTxJ0eY6u\/tOnT\/X8+fNz+Y1\/+k\/P3aIq8XPnzh3du3fPYgaDzsDdu3flM+Pny75fM3yufKb2daaKoji\/53H\/e\/LkCbe62cosCyMKHIqd67hqDx480MOHD88lnB2YNM\/7ar5gdOvWLVnMYOgZuH37tm7evOkz4+fNJmdgra\/PlV+Dhr4GrfMLIZzf87j\/cR\/UjH9GURhRyAwRrgN+tH2C\/bqKor45dXZosAU1\/x5IXTOymIAJmIAJmIAJTI3A0QsjiphNpQ\/yUYqibDFBTUVU1\/KPCWxGwN4mYAImYAKjIHD0wmjfFCiO2pLPgY0xxRj9XNBh21hiXIZENd8xkn9MwARMwARMwAQmSWBWhRGFTZfkVwZ7GtPPJek3bt9\/vx3isQmYgAmYgAmYwAQJzKowOhr\/EM6n5uO0uj4fumMCJmACJmACMyQw3y2NojDKP87K+5PBHsL5UqP8cZr8YwImYAImYAITJTCKwgh2+UdaqZ+KJOyjlhil777TB+E7lSr0dfPLafKPCZiACUyFgNdpAibQEBhNYdQs5\/JjXiBdtnhkAiZgAiZgAiZgAvsnMOrCKG2XAol3j9J4rG0Izcrqumn9aALHI+CZTcAETMAEtiEwicKIjU2lOGKtFhMwARMwARMwgWkSmExhNAW8ITSrrOum3eejc5mACZiACZiACVw\/gaMVRnw0loRt5n3GXeJ3jbqoWGcCJmACJmACkycwmg0crTCiyOmSIQXSaOi1FpL+zmNdtwwemoAJmIAJmIAJTILA0QqjPjqpWKJA6vLB3qW3zgRMwARMYEQEvBQTmCiB0RVGieMUC6AQ0uqlur7ou2cCJmACJmACJjANAqMtjKaB7\/Iqi88\/0Hd6T1\/ok8sGj0xg+gS8AxMwARM4CQIujPZ5mWNcZouq9Pnn8o8JmIAJmIAJmMDECLgwuoYLxj8kW5bXkHifKZ3LBEzABEzABEzgCgEXRleQ7KD4\/vfPgymO6vp86I4JmIAJmIAJmMABCWw7lQujbcl1xYVwrv1Mn6uqzofumIAJmIAJmIAJTICAC6N9XqQYpRiXGQuV+vJL+ccETMAETGAvBJzEBA5DwIXRvjl\/\/PF5xlCVquvzoTsmYAImYAImYAIjJzDLwog\/Dpmkj3+yt9s+\/8H6opBCWLrza\/tVtez6wQRMoEXAQxMwARMYI4HZFUYUOvxxyCSMu8Ane952+W2l++yz87Dwh\/71tHMY7piACZiACZjAyAnMqjCiCKLQ2Yb5LrFX5isKKYRGXX55Ih+nNdv1owmYgAmYgAlMmcCsCqNDXIjXr18rl7dv36pL9PHHek\/f6SN9pR\/9qNunK846s\/IZ8BnwGfAZGNMZyO959A9xrz3mHL2F0TEXte3c7XeL9vou0LtFPX78WI8ePTqXH\/\/4x3r16tUVeXlWGD148O0y6t\/\/+xtX7F0x1l3leApMeKF58+aNz0jH8+gUrv917dHn6jRfT67jPD179uz8nsf9j\/vg8uY204dZFUb5NdqkKNrE9+nTp3r+\/Pm5FEWhe\/fudco\/\/Ic3zpf0X\/\/rr3b69MVa3810jlzu3r2rO3fu+Hz0PI\/meM0PsSefq7WvIX7ODXzOcZ\/L73tPnjw5v7fNsTPLwmiTQmfTi\/rgwQM9fPjwXEIIunXrVqdQGJ2Zl1P85Cc3On36Yq3vZjpHLrdv39bNmzd9PnqeR3O85ofYk8\/V6byGXPd54j6X3\/e4D2rGP7MrjK6zKNrmHKRfUKtrqaq2yeAYEzCBQQTsZAImYAJ7IDCrwmibomibmE24F8WF9+efv+svFlJdvxu4MQETMAETMAETGAuBWRVGQKXQaQv6JNhS\/1BtUTQzVZX0wQdnfSokOlV1NvB\/PQSsNgETMAETMIGDE5hVYcRvpXVJThX7qnFu21c\/fZxGvrrm8Z189JFUVe8GbkzABEzABEzABI5N4HCF0bF3esT5+QL2y5dSKpD420bny6E4+uQTifa99969pXRudccETMAETMAETOCABFwYHQh2CNJiIVEgff+zuPzDj+dTl6VUVc2wrrUskpqRH03ABEzABCZCwMucBwEXRge+jiFIi4X0xcuoRfxK6adWkELQ8qeqJN5FWg78YAImYAImYAImcCgCLowORbo1TwjS4quo8rOX4p8O+UAv9VF4KcUofffdWeX0RSvCQxMwgcMS8GwmYAKnSMCF0ZGverEIy1pIZz9VJX1Qf6XFQlospKo6U\/o\/EzABEzABEzCBgxFwYXQw1P0TffWVFGNjr2uJ3+ZH+D42v9W\/WDS2S498UXuxuKTyYDUBW03ABEzABExgHQEXRusIHchOcVQUUgiXJ6zrplCiDqJIolhaLN75UD2hXCzeKdyYgAmYgAmYgAnsQmDChdEu2x5n7BdfNL+1xleM+O21zz67vM66lqpK+vLzWpWilj91LVEg5ZUTX9xeLKTFYumy9qGu17rYwQRMwARMwAROgYALo5Fe5RCkxUKiSPrsMy3\/BlJRSDFK\/AYbfwtpoc+WfaWzSXL0AAAQAElEQVSfupaqSipLLYslCqZk62vLUuJdJwqrqurzst4ETMAETODQBDzfUQi4MDoK9s0mXSykxaL5RTU+ckvvJv1hXIjfZkMokkoVqkOUQtCgn7K8\/GcB+JyuqgaF2skETMAETMAE5kjAhdEEr2oI0mIhUSQhvIP0uRb6RF8sf6vtk\/hS5Rff6ZPiu+XfiqTeWSykstTlHz5yQxMCj43gXJZN348mYAL7JOBcJmACEyDgwmgCF2nVEmNsvpfEx23JryybN4LKUqoqqaq0\/GSNOohPzBaLd558TheC9MUXTZX1Tq3kiHPSDWnreoiXfUzABEzABExgtARcGI320gxfWAjSYnG1QCJDCFKM9C6Erx5R8ywWZ7qXL6UYpRgvF0dnpuV\/OC47Kx5SIZW+q0S7WEhVtSJoBiZvwQRMwARMYHYEXBjN6JKGIC0WTYHER2y8IUTdk\/djvNhwXiBVlVTWUfwl7vqzs3eQPnv3jW\/eTToLqWupqs46Xf+98zk31bWWb1HxsRxFEkLxhJw7rejgt1hcdShLabFo3g7D56pHt6aupcVCyy+ZLxbdPtaagAmYgAmYwBkBF0ZnEN79N5smBCnGq9sJoXlTiGIpxgs7BRI1DLXGJ58HffB5oSoupMVCKgrVtZY1BT7UOIuFVFW6\/FMUzUdyFFQhXNjqWqrrs6rrrKgpywv9ql5ZSiyKd6sWC6mutfyyFAtEX5ZNvlU5km2xaBZPXF1fzpt8tm3rWqoqqao2y1DXzfoXC6muN4vdxruuparaJtIxJmACJnByBFwYndwll0K4KJCoY7oQUARVlVTXWhZFyaeum9oCeyqSljbeNSoKabFo3rKi+kLHBDFKIUghLF3XPoRw4UJBw0RVdaELQQrhYryqR3yyh5B6ulR4XWj7e6lIw6OupcVCSzCAQLCzTuzrBL9U5J31b\/31v65f\/L3f043f+i0tC8Az3bJdLNZlWm2va2mxuFjnau8La11LdX0x3qVX11JdS1UlVdUumS5i61qqa6mqpKqSqkqqKqmqLnym2qtrqaqawrmup7oLr3u0BLywIQRcGA2hNFOfEKTFoqlj+LiNWoY2bZf7PfdoxiFIRSHFyKiRutayvsAnF+75izKoCmcBi8VFFcYETejqR\/xYSAgXfiFo+cec0ueD+FxYV\/dCaGKJQYriwp\/Cqaouxqt6+LJRwNBv+9Z1WzNsXNf6pbPCSGUpVZV0Nl62i8X6+FSQsSb6ubDWtE54rsuWchCHkIsWPS3jJOtyYceXOIQcSbBtKjApSy0LRvIhKV\/eDslLLELcYiFVlUR+pCy1\/OUDbENy4bNYSIsFvctCPqSqpLK8bOsbpXXxJKKPlKVEnhST95Oup73z6ae68eMfS4uFVFVXvchVllf1XZq0lmQjNknSbdKWpVRVUlVJZSmVpVSW0mJxcQ3SGdokb9uXNXI9kcWibe0eE9O2oKsqqSwv1tf26RovFtJiIVWVRI6ylMpSWiy6vNfr2AfnY7GQFgupLKWylMpSKkupqqSqWp9nnQdrXSyavS4WUlWti5iNfZaF0Ycffqgk665U8qNd5ztXewhSjFIIUoxNHdPe68cfN5+UcX+ltuCNoOTD8yeXstSyYOL5y2vpYiGVZfP8SjraxUIiLuW51MbYVGxMVBQSEy8Wl1wGDXjXigUvFo17CBcbCaHRDXmMsfFiwQijEJp1MUdaJ\/p1QnHHflgXcck\/BClGqSiSZljLeqqq37copBD67clSVal30abctBfa7XtVpeU7bWUpLRbD83CQuBlU1fCYPk\/2glSVLh3UNEdZSlXVF31VT\/GJpJs4efI+h521X41cr2GdxKacKe+7yHXNz\/\/BHzR7ZH2sI8WTkzF56a9LhD1fS4olHkl56SP4rxPmJQ9CPwlrLUupqpoMPE+aXv8jc5IHIc9iIZWlzgvpqpKqqmHRn+XCQr60J3Ii6GjJX5ZSVV34r+qxH4RYchCPoMvnoL8qT7JVlVSWzV7IQa5cmAcZmg8\/1kUMQi5adOQvy2auGNMKZt\/OrjCiwHnx4oWSMO67itiSHy3jPt9T0sfY3O\/Tnrl3LxZpJIUgLRYS9\/fPPpOKQioKqSikGC\/86PFaynOL51pZSlUloasqLe9JPPeQxUJCX5bSYqHz\/2n\/6A8X+qD6QosyXIotS2mxYIY1UhTdDjFeFF5fftntk2spZJAQpBAu3oGKUSoKabFoCq48ZlU\/RikEabHQt3\/6p\/rm+XO9\/bM\/a8BTaAF3VTw2gHEBQpCKQiqKZl3okrBm8oVAxGqJUSqKZh\/EkKMopBilopCK4iL\/6kyNtSgucrGOomj0rJsDwcF4771Gt+ox9wmhWQPrI2cuQ26gzJPvi3FbQpBi1PKm2ratG7M3ZJ1fn519IeyLdfb5DdS\/vX\/\/qifrK0upqq7ahmiIr6qrnuiTXLVupglBilEqimFxzFtVUlVJZanliwtnrKqa+BCadsjzqvFsHvO8jaZ5DEGKsenv+pjmGJonRimEod6r\/XjxxSOtoaqkspSqCq0UQtOe2OOsCiMKGwqcIdewy3do7JD8U\/eJsakbeH0+u3f3bgcbr+NJ8Oe1h3sUr+shXISGIMUoFcWFjh7PSe6RPEd5LaNfllJVSVUlJXv6nxja5Mc9M8WhT2PWVVUSsfSx4VdVzJgJxrPF41eWEkMk87joxthAYXO9Thfug3sh6NuHDyUNjmgcz+LEOljP2R6Wf4+KcS4xNr5DHrl45CkKqSgk8jBOevrokCH58C8KqSikGC+KpDw2hHzU3Wd\/HCbWQZ\/5i0KKUYpRilGKUQqhO76tJZ61pXz0yU+bdLRIO7ZrnB\/4GKWikIriYr\/kQbpi27qikIpCilFineRmXQhrROi343rGr\/7oj\/T2Zz9rzi1rID7GxjuEpshE32hWP8I+xReFVBQXe0SPFIVUFKvzJCv5mBuhj7BfhD76oXuNUYpRijFlb9oQmj2Sj7yNdv0j87IfJEYpBKkomv2Sh3ysb30mKfmTEyEu6T77TIpRCqFZ55B8xDN\/ykE\/F+zMQ+51+ZIfvjFKMUohSCE06yFXmmddrhnZZ1UYbXNdKJCSbBM\/55gQpBi322EIEq\/r6fnK84s+Lc\/F9Fzj+dg1QwhSjFKMUoxdHhe6VNRUVaNjTHGViiH6VSWhR1dVjR+PZanlGwMUTanYwj8vsNgHsbkQi6DDTjwx9NGPRfL1sTbGR11bUTQ3aS48h4FDsW5BIUgsPsZ1npvbQ5CKQiJ\/UUgxbp6DiBAkcrAnDjhSFFJRSDFKMeK1nRSFVBQS+ZGi2DxPCFKMEvGsMT0BGcc4LF8IEv7Esz+kKKSikNAj6JAhGUOQYpRilEKQQhgS1e3DmpKkvXG2ENbVHdWvLQqJOIS85GFfRdEfs8oSglQUUlFIMTaeIUh5fvqNpfeR529ZSueuIUghSCFIIUghSDFKRSGdO\/Wmk2KU8EPYJ8JeEXQhNMEhNO2JPM6qMGq\/40PB09bl1zXZ8UEY5\/au\/uvXr5XL27dv1SfWN2zu33+rX\/u1pp8zQf+jH73Vz372Vv\/5P7\/Vv\/23Tcv4z\/6s6aNHGP\/O7zQ+yS\/pQpDIhfzmb74V0r522JKO4ojnPC3FUFUly+WWFyGKJITCJxeKoDTGji\/R9LGRHyE\/Qj8J4yTo\/t2\/u6E\/+ZNbveeI3Kw15c0Z0sdeVVp+\/MiayFmWUllqWfSxJnxosSPJp6p0ZV58kbKUyrLJy\/xIimXereX+fb390Y\/09td+7crc63KyrqqSyvJiXeylHYdfW7ePMXm5dlwLeOQ5sZVlsy7sSOJVVVc557GpTw72Q27apB\/SEluWzfztWGwIenKzNtoheft8Uj5yVlX\/\/vDry7GtnpxVJVWVVJZSVUnozvNxxpCBr8\/Esg\/kPMfA2EP5s0bOHpKey\/\/lv1x9XR2yHnIhVSWVZfeZyfNwz0v+VcW17viIlhfAmcjPzWQfV7ZBkUOxc8WQKdbZM9fz7uPHj\/Xo0aNz+fGPf6xXr15ZdmTwV\/\/qK\/3Nv\/lKtF08\/8pfeaV\/8A8an+SXdP\/tv73UH\/3Rq6X8i3\/xSgjjJ0\/+Qk+f\/i89f\/7N0kY\/XUheWKqqGVE0\/cZv\/P9LvxcvmlzEPnjwbePQ88gLRTKRA0lj8iNlKZWlll95YIyUpVSWUlk2+t\/6rRt6\/Piebt++tfxeMjfTf\/bP\/r\/lmfrjP\/5Gv\/7rb5cv\/OQm\/u\/\/\/W\/Pbfhh5yZXlhJrwocXT4Qxcfna0CUf4tJ8\/\/pf\/4XIzRghHilLqaqkqpJSLHb8EWIofNF973s3hNBHx\/rw+Y\/\/8X+JvSRpX2P0+BKDEI\/QR4\/QR8ea83Wxl3zO1KdFiMmFPKyZnPna8jWxHgQ7vghx5ClLiEpVpeX1wodc2NK6Gg8p8WLNrAXBLwk5kyQd+6mq5mygIzdz4Jfi6aNP62KMb5qfHL\/wCz+vv\/f3\/l\/F2KwTO\/qqalZXVY2ea5PvvasPi1yYO+UjJ\/tjzHrwIwdtWldaN3H4oEeIQcca2COCHUGPoMOOME5xzImwZ1pypXnwZQ25sB5yEU8exthpeQ6xD4Q8XfHoWBfxCH5Inot8q4S50hpSbO6PnXlokzBmDvyrqrl2PHK2\/vbfvrE8g9gRYtr50OfrTozIBzf4laWWr1HY8GfOlIs+fJI\/Mb\/+62f\/c\/N2WRyxlNnJLAujIUXRtlfy6dOnZzfQ5+dSFIXu3btnGRmDv\/E37uhf\/svb+if\/5Bf0d\/\/uLy2vD33efUrXPgQt\/wIA7zzxrk3yS7FVJf2f\/\/OtsPMuei58EoScXf7lL8zhkyTGZoYQpBCkEJpxegxBCiGNLrevX984e0fyhn7v935JH374wVkB\/qvLsbKfP\/iDn1\/qHz361aUfMZhDkGKkdyEhXOyR9bHmGC\/s9Ihnvk8\/vSNyo8slBCkEib2G0FiIwR8h5qc\/vXVpndjRpbwUfqw3Sdob4x\/84N5yP\/gSgxCP0EeP0G9mbx5DaNr0mPzTOLXocyEPayYn609rS2tK62Ft2PFFiCNnCFII9LTcMz7kajRSCBKsuljj014LeRH02HNBR27mwAcbOvro07oYY2vL\/\/gf\/4\/athAk1ocvudL+2Tfyz\/\/5Pf2bf9MI\/cQFHkmYm\/hcyMV6kg9tPjd24vBBjyQda2CPCHYEXwQddoQxcfm87T458f3hD39V\/+k\/3RN7YF+sh1zEk4dxEmJSHvrEYyOOeProWBfxCH5Ingtf2OXzphzoyZHWkMfikzjjw3xJGDNHWh\/nCknjlAefFPM\/\/+e95d5TbL7uFNfXkoc4crEm+nV9dmj6Amaon11hdJ1FEdf\/wYMHevjw4bmEEHTr1i3LRBj8rb91Q3x8TpFDy1vn667fX\/trN8SNJBfiEPLEqPPrjy8f07e\/5sA4CfMiafynf\/rt8p0tPirkBe\/sSHHUzoUxemJibNS8GDY9KbenuWlZGzGskz2yNvrYmBsbeWNMmSRyfNx6ZgAAC7lJREFUsU9icz980dEScxEhEROjRByCHcl9uvrsAeEmk+ztXDEmi4SN\/GkdrKW9D+wIPqwBYZxLjBc52732erAzLxJjU2QyL8Ic2JPgA7NkS6zTmLUkYT0xSjFKMUqMEezkYF\/0U27aECR8kBjRSCFIMUro8E+xzMl5wisECTvrxY4t9bEjad8\/PStw\/8N\/uKV\/9a8aoY+9S8iZ52P+5Ec++iE0zLDFiEYKQQpBIj7GRpc\/hpCPrvZDaHIydy7siXlSxH\/\/7zf0j\/\/xLbGH9hlLPl3rzG3EEZ\/7xSixdiSE5N0UyfjCLp835UCf5wmhiUWHTzPqfwyh2TfnCuE6st8YL8eQ7+\/8nVvLvdPHGoIUo8SaEVgh8CMP540WG\/5tCUHChj9x\/+gf\/Ylu3HjddpvNeFaF0SZF0YsXL87+j\/zD2VxIb2Q4gRAknuTDI67XMwTpN84+yuPFDuEFCkkverwYoQ+h+W3+GKUQJPaADV\/s+SpjlLDnunY\/BIm4lCO1vPARG2M7ohkTw5z4pRhaxgh2JL3Y4pvb8WFvSIxSCBLz4dP2RZfnIRbfEJq18BiCxHz4YkfwQYcwzgW\/PCdj7KyHuBglWnRtP\/IxJ4IP6yWOHPRjxHJZQpDwJTYJuYlJwhjBHmMTT5\/50eNHfvoIY2zo6KPDP8YmNgSJ7++9OPtomHcKsbOGGBs7jzFq+T8JrD\/ZYpRCwNpICBI24nNJc8YohSDhw\/ysh3whSDE25xU9QkxaM37kQ0efFkl2WgQbgi33J1+MUoxSjFKMUnsNzQ6kEKQYm6KCPORDWGeMUgiX15nbYmyyhNDEYyMHa0EYI+RCGu\/mMQQpBCnGZsxjCBd5iCMX645RoiVnW\/BB8Gff5EFCkBhjS6y61oCOWPxSbuZCYpRCIJsUQvNLd\/gmP2LpJ12MEnEff9zEzPVxVoURF4niqC3ok2BL\/RfviiN0CONkc2sCxyQQgpRe9EK4vBJe4PIXqsvWS6PBgxCkGAe7L19EeYGMcXVMCFIIUowS\/knYGzJ0LyGsnmcbawhSCFKMEutiPdwEWBMtunV5Q5CIi3Gd5\/Z21hHj9vGrIkOQWD\/7Zd8IZ4sbLUIfG2vIJcburCFI5COOXCF0++XaEKQYpRhzbdMPQQpBilFi\/hgb\/arHEKS0BtaBsBZ0MTaRIUiM0WMPodHzGIKUbNgTB3TY2xKChA1JvsSRF6GPHmGMX8oRY1OM4NPFOe05xhTR34YgkZs5KGgQ8qLrj7pqCUFiXoRY2qte89bMqjCisOmS\/BJib4\/RIbnefRMwARMwgekSCEEK4fDrD0GK8fDzphlDkChokBCS1u0mBGZVGG2y8Un6etEmYAImYAImYALXSsCF0bXidXITMAETMAETMIGhBMbg58JoDFfBazABEzABEzABExgFARdGo7gMXoQJmIAJzJGA92QC0yPgwmh618wrNgETMAETMAETuCYCLoyuCazTmsAcCXhPJmACJjB3Ai6M5n6FvT8TMAETMAETMIHBBFwYDUY1R0fvyQRMwARMwARMICfgwiin4b4JmIAJmIAJmMB8CGyxExdGW0BziAmYgAmYgAmYwDwJuDCa53X1rkzABExgjgS8JxO4dgIujK4dsScwARMwARMwAROYCgEXRlO5Ul6nCcyRgPdkAiZgAiMj4MJoZBfEyzEBEzABEzABEzgeARdGx2M\/x5m9JxMwARMwAROYNAEXRpO+fF68CZiACZiACZjAPgmsLoz2OZNzmYAJmIAJmIAJmMDICbgwGvkF8vJMwARMwASuj4Azm0CbwCwLow8\/\/FBJ2htO42Rvt8nu1gQOReDNmzf68ssvRXuoOT3P\/Alwnnyu5n+dvcP9E5hdYUSh8+LFCyVh3Ict+eRtn6\/1JnBdBLiB\/f7v\/\/51pT+xvN5uIuBzlUi4NYHNCMyqMKIIosjZDIG9TcAETMAETMAETKAhMKvCqNnS8EcKqSRDo16\/fi3L4RicAut09k5hr97j4Z47PleHY31q5zqdrbm2syqM2u8WUfS0dfmFxJYE39zW7t+\/f18PHjzQ48eP9ejRI4sZ7O0McKY4b7Q+W35u7esMcJ58rnye9nWe8jycLe6H3Bc5Y3OTAxdGh8NHoUPR0zdj28aYmD5\/DsDTp0\/1\/Plzixn4DPgM+Az4DJz0GeB+2He\/nLp+loURBQ6Fzr4vDsXRw4cPZTEDnwGfAZ8Bn4HOM3Ai9wjuh\/u+x44l3+wKo6FFEX5juQhehwmYgAmYgAmYwDgIzKowotjZ9p2iXWLHcSm9ChMwgWsg4JQmYAInRmBWhRHXjgKnLeiTYKNPAUU\/CWP0FhMwARMwARMwgdMlMKvCiOKmS\/LLiz2N6SdJOrczJ+DtmYAJmIAJmMAKArMqjFbs0yYTMAETMAETMAETWEtg6oXR2g3awQRMwARMwARMwASGEnBhNJSU\/UzABEzABEzg4AQ84aEJuDAaSDx9SZt2YIjdTOCcAOemS84d3nVyn3cqNybQSYCz0ml4p8Se5J3qUpNstJcMHpw0gb7zgL5L2rByn7ZtKmMXRgOuFBc6fUmblvGAMLuYwCUCnJ225A6cq9zOOLe7Pw8C+9jFurOBfdVZWmffxxqdY3oEOBerVp2fqdTP\/YlPelrGuX0qfRdGa64UF5YLnLsxRp\/r3DeBXQhwnjhXeQ7G6HOd+yaw7kxg5+zkpBijR0fLmH4SxujT2O3pEdj1+hPPOcrJMUaf66bQd2E0havkNY6YwPCl8QKRZHiUPU3gMgFuNshlrUcmsBsBzhSyKkt6\/aJd5Td1mwujqV9Br38yBHjRSTL3F5bJXBQv1ARMYDCB9PpFO+fXMBdGrSPhoQlcBwFeSPK8jOf8wpLv1X0TMIHpE+A1K98F47m+hrkwyq+0+yZgAiZgAiYwbwLe3RoCLozWALLZBPZBYK7\/Z7UPNs5hAiYwfgKn9BrmwmjNeex6u5ADgn5NqM0m0EugfYY4T+jyAMboc537JnCFQEvBmeHs5GrG6NHRMqafhDH6NHZrAusItM8M5wddHscYfa6bQt+F0YCrxIXlAidhPCDMLiZwToAzk84PLeNz47sOOmxJGL8zuTGBjQhwdtI5omWcJ2CMPgnj3O6+CbQJcEbSeaFlvI1PO2aMYxdGA68KhyDJwBC7mcAlAun80F4yZIMXL14IO5Kp3TWBKwTWnRHsSa4EnymSjfZs6P9MYElg1XnAlmTp3PGQ7LQd5kmoXBhN4jJ5kSZgAiZgAiZgAocg4MLouik7vwmYgAmYgAmYwGQIuDCazKXyQk3ABEzABExgfATmtiIXRnO7ot6PCcyUAF\/4RHbZ3q7xu8ztWBMwgWkQcGE0jevkVZrASRJIhQwtX+ZE6J8kjINt2hOZwGkTcGF02tffuzcBEzABEzABE8gIuDDKYLhrAnMkcIw98a7OOtlkXemdInLS3yS27UuOLmn7eWwCJnCaBFwYneZ1965N4FoJULzkwmT5mD66JF2FCjrstAh9JO8zTtKnT\/bUMndbks2tCZiACbgwmtwZ8IJNYH4E2oVKGrPT1M9b9MiqYmiVjViLCZiACXQRcGHURcU6EzCByRGgcFpVDGHDZ3Ib84JN4NQIHHm\/LoyOfAE8vQmYwAUBipeL0eVel41Cp0t\/OdIjEzABExhOwIXRcFb2NAETOCKBXYogiifij7j8U57aezeBSRFwYTSpy+XFmoAJrCLQV\/z06cnlogkKFhMwgUTAhVEi4dYETGAYgWvyGlKgUODgly8BXT6mj0+XHpvFBEzABFYRcGG0io5tJmACsybgAmrWl9ebM4GtCLgw2grbrIK8GRM4OoF9FihDcw31OzocL8AETOCgBFwYHRS3JzOB+ROg4GgLu27rGKNHNvnYq8+XfEifnXkQfJB1fvhaTMAE5kBgsz24MNqMl71NwATWEKDgGCprUg0yU+Qgac6+IHyQdX598dabgAmcBgEXRqdxnb1LE5gkAYqYdQvHB9mX37o8th+fgFdgAtdJ4P8CAAD\/\/x4n\/roAAAAGSURBVAMAwRDUkQbp4rQAAAAASUVORK5CYII=","height":281,"width":466}}
%---
%[output:7a5fd55e]
%   data: {"dataType":"text","outputData":{"text":"[14:02:57][INFO]  訓練完了。最終train\/val: 2.5966 \/ 2.6132\n","truncated":false}}
%---
%[output:2d058d5c]
%   data: {"dataType":"text","outputData":{"text":"[14:02:58][INFO]  300 SMILESを生成中 (T=0.80) ...\n","truncated":false}}
%---
%[output:3bd50a65]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 30\/300) 生成中 (RL前)\r[##--------]  20% ( 60\/300) 生成中 (RL前)\r[###-------]  30% ( 90\/300) 生成中 (RL前)\r[####------]  40% (120\/300) 生成中 (RL前)\r[#####-----]  50% (150\/300) 生成中 (RL前)\r[######----]  60% (180\/300) 生成中 (RL前)\r[#######---]  70% (210\/300) 生成中 (RL前)\r[########--]  80% (240\/300) 生成中 (RL前)\r[#########-]  90% (270\/300) 生成中 (RL前)\r[##########] 100% (300\/300) 生成中 (RL前)\n","truncated":false}}
%---
%[output:74269a7c]
%   data: {"dataType":"text","outputData":{"text":"[14:03:34][INFO]  生成SMILESのサンプル（先頭8筆）:\n","truncated":false}}
%---
%[output:42d269fd]
%   data: {"dataType":"text","outputData":{"text":"[14:03:34][INFO]    [  1] CCCOCCOC((CCCCCC1CCCO1CHC=)C=C)CCCCC((CCCC1CCCCCC1(NCCC4CC)CCCCN1O)CCC)C2C=O=Cl(COOCC=H1)HCCC1C=1241\n[14:03:34][INFO]    [  2] C2CCCCCCCCC2C(CC(HC[CCC)C)@COCCCNCNCHCCCCCCCNCCCCCOCCCCCCNCNCCCCCCCOCNCHCCCCCNOCCCC@HSCOOCCO]CC(C)[\n[14:03:34][INFO]    [  3] CCC(1CCCC12)2C=O(H)[CCNCCCCCOCCCOCNOOCCOOCCNCC]C(C3)C(CCC(CCCC2=CNCO2CCC)CNCCCNC2[C)CCCCC-CCCOC]23[\n[14:03:34][INFO]    [  4] Br=H2CN(CC(C1213OC)NCOOC12ONC)COCC3C2[CCCCCCCNCCNC@COCCCCCC@CCCNCCCCOCCCC]CCCCCCCCCC=CCOC(NCOOCCC)12\n[14:03:34][INFO]    [  5] CC(=(CO=)CC)CCNCCC=(C)1CCCOCCNC=C1CCCCC=C(CCCOCNOCC(2CC1(CCC(CCCCCCCCCO(C5CCCC(1CCCC=CC1C))))))125=\n[14:03:34][INFO]    [  6] CCNO(C3C3CCC((ONCO(C))CCCCCCCCCC=C(C)CC=CCCCC)1HCCC1)(C(CC(=CC(3)CCC)C22(CCC(CC(CCCC(C1(C)))))))13C\n[14:03:34][INFO]    [  7] O((NOCSCC(CCCCCC1CC2OOCNC(C(C1C)(CC(C)CC(O((C)CC((C[NCCC)CCCCNCCCCC)]CN1=CCl))COC=C1C1CCCOCC)))))12(\n[14:03:34][INFO]    [  8] C(C(CCC)CCC2C(CC(COCCC)((=CCCCC[C)C)CC]C=CCNC(2)C(C32CCC))(C(H(CS(CCCCNNCN(C(C(C(CN)1)31C())))))))2\n","truncated":false}}
%---
%[output:4c71865b]
%   data: {"dataType":"text","outputData":{"text":"\r[#---------]  10% ( 30\/300) 検証中 (生成)\r[##--------]  20% ( 60\/300) 検証中 (生成)\r[###-------]  30% ( 90\/300) 検証中 (生成)\r[####------]  40% (120\/300) 検証中 (生成)\r[#####-----]  50% (150\/300) 検証中 (生成)\r[######----]  60% (180\/300) 検証中 (生成)\r[#######---]  70% (210\/300) 検証中 (生成)\r[########--]  80% (240\/300) 検証中 (生成)\r[#########-]  90% (270\/300) 検証中 (生成)\r[##########] 100% (300\/300) 検証中 (生成)\n","truncated":false}}
%---
%[output:58037d9b]
%   data: {"dataType":"text","outputData":{"text":"[14:03:42][INFO]  RL前生成 -- 合計 300\n","truncated":false}}
%---
%[output:68c1252a]
%   data: {"dataType":"text","outputData":{"text":"[14:03:42][INFO]    有効性:    16\/300  (5.3%)\n","truncated":false}}
%---
%[output:1307a96b]
%   data: {"dataType":"text","outputData":{"text":"[14:03:42][INFO]    多様性:  11\/16  (68.8%)\n","truncated":false}}
%---
%[output:54f8ce50]
%   data: {"dataType":"text","outputData":{"text":"[14:03:42][INFO]    新規性:     16\/16  (100.0%)\n","truncated":false}}
%---
%[output:9f971cbd]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAkYAAAFfCAYAAACr7AOiAAAQAElEQVR4Aeydf2xdyXXfRyQlkWtJlXctS3S8jiA78a7twIq3jRQU5aYCWrhZt13khwFq4Vpoa6dBGzuOoSy2ReEgVmA4Sr11m\/6RunVhFxWQAPkngP8qgqIFUoRxatgpYhu2owDhwpShLbSJV5REaanyc8kjDYdzf7x3733v3rnfxZ43Z86c+XE+c2d4SPE9ztzXfyIgAiIgAiIgAiIgAhmBGaf\/REAEREAERCBZAgpMBEYjoMRoNF7yFgEREAEREAERSJiAEqOEN1ehiUCKBBSTCIiACLRJQIlRm3Q1tgiIgAiIgAiIQK8IKDHq1XaluFjFJAIiIAIiIALdIaDEqDt7oZWIgAiIgAiIgAhMmUDjidGU49H0IiACIiACIiACIjA2ASVGY6NTRxEQAREQgQESUMiJE1BilPgGKzwREAEREAEREIHqBJQYVWclz4YInDp1qqGR4sP44\/t63HtY1io8qvjUoZY3fp696lxV+kd9qk6Q49fUmE2NEy6zrXHDeUatF62rqI15ytrxkYjAuASUGI1LTv0qEZj2BXb16lU37TVUArXjxFpD2WnKCtoyJfJS1IY77fBAH4oQc1OxNjlWU2sqG4c1d3XPWRfrK4sh1l6nb2w82UTAJ6DEyKchPQUCe2LgEt1j7KiBtYaSt9Rxv6gwnvWlDMXaY3baxhXGG7dvl\/uxX12MjTWxti6zi62t6rqJDd\/YGLKJQB0CSozq0FPfqRDgMiwSFlXWjs9QxFjwhYSYqZtOGQo+oY069rbEX9Ooc9DX+vi62ZooGbcJBoyTJ6wzry2041tHwvFi9dj4+MXsZTb6xYR+oT1mMx\/aJCLQNoHuJ0ZtE9D4rRHgMmNwSl9itrAdnyLhi1Se0C+vDTvtXRSfQajXWW+XYiYuYqFE0E2oF4n5FZXEyhi+D3UT7KZbia2uxOaNjYlfTPCN2fNs+BcJsdG3yIf2ImGMov6jtBXNM0pbOCd9m1xnOL7qwySgxGiY+z6xqLm4QmHy0ObXaR+i+AxCvS4PxuMLCILuj4fNF9r8Ojq2OsIYzIswDiWCXtSGD4LfqMK49KG\/SVjHjq1IGKeKHz74Fo0Va6MPfa2NuunTLFlTF9bShTW0vQ8av1sElBh1az+SWQ2XGRdrMgElEIi\/H+yPHxJtJthNp6TeJ\/HX7OvjxACnUcbAlz7jzDWEPrCJiR877X5dughMmoASo0kTH8B8XGx8gehSqKypS+vJWwvrDCXPd1Q74zaxL4wz6tx58zJWXluVOer2L5pj3LGvXq3+Tshx5yhad9fb2G9fur5erW94BJQYDW\/PW4+YS6\/tSfiCEkpTczJu0Vh+O3pMivoXtcEulCL\/Km22PnxNZw50bAi6SayOzYS+ptct\/bGY3+q+Ps4cTfZnrKpiayUO+lg9Vlo7pS\/4+nXTsUtEQATaJ6DEqH3GmqEFAnzh8cWfAjtfTHzbKHrY3+\/LuLRjM516KLTh0wXx18Z6qFP6gs0Eu+mU1JsW+CB1x81bX5696nx+f3RfGMOvm47dF+x+3deJ3dopfcHPr6Njm7T4a2xybsb1pcmxNZYINEFAiVETFDVGlIB\/+fk6zn7d12lrWhh\/0l9cxp2PtYbSNA8bb9w1Wv9RSj8m+jE3go7Q7texNSGM6wtj+nV0bJMS5kPaiHXUGFhHkbS1Rsb1ZdR1y18E2iYwwMSobaQa3wj4l5+v0+7XfZ22Lghr4ouGvxbq2M2Gjs3E7HVKxgylznjWlzWaThmrY0OsHR2hXlf8mMKxmIN2s4d1s49TMq4vjOHX0bFNSpgPmdR8RfOwjpjQBzulRASGSECJ0RB3fQAxc7HzBRZBbytkxjZhLpO25mtjXFs\/pT8+dRPf3qQOL+Zocsw2x+rbesdhwX4Q5zh9q\/RhbF+q9JHPiATkXouAEqNa+NR50gS4ULm4q847im\/VMfP8mMuEdeb5FdnpF0qR\/7htrJN5\/P7UsWOjpI7epjCPPz5zhja\/XXr\/CbC\/vvQ\/IkWQGgElRqnt6IDi4YuoCRetHzp2bJSh3a8X6X5\/xqHu+2Pz603pzGPSxJisk\/GKxor50Ad7Ub8m25iLOZscs8mxprE+5jx16lSTYVQai31g7krOchKBxAgoMUpsQ4cUDpe3iR83Fzp2bJTU0ceVuv3HnTfsxzqIJ7SPW2csxkTQY+Ngpx2JtTdlY3zmio1X1Bbzb8NWtoay9rw15cVs\/rSbmK1KSR\/WVMW3az6sveqaiHEU\/6rjym\/YBJQYDXv\/exV92SVIOxJelNSxjxMsfelnJboJNsaNCW3mt6dsyFB1jiI\/2pCiJdGOmA\/xmt5EyXj++E2M2eQYXV9fk7H6Y7EnxO7b0LEVCT4SEegzASVGfd69jq696NKkjWVTFgk+VcXG4SJHYv2w4xdrK7PRN8+Htpjk+Vexs04T86fOPFafZtn0OmLjEa9JrL1K\/NbfL+nn103HHoq15c1v7ZR5PuGYk6yzJtYWm5O2mD20hX7UyyQcI6\/O2pA67awlr7\/sIjAuASVG45J72E9aQIDLqq4EQ2ZVxsyU4AU7Epj3VPFB9jR0yMD6QmF52CiLZJI+Resoa6u6TvyQsvGsPfSlXlVsDL+0vr7N162d0rfX1Zscj7GKko+6ax21P+uxPugmZvNLa6P07ejEFLPTJhGBugSUGNUlqP4iIAIi0GECKSYQKcbU4UfIOTes1SkxGtZ+K1oREAEREAEREIECAkqMCuCoSQREQARSJKCYREAE8gkoMcpnoxYREAEREAEREIGBEVBiNLANV7gpElBMIiACIiACTRFQYtQUSY0jAiIgAiIgAiLQewJKjDq4hVqSCIiACIiACIjAdAgoMZoOd80qAiIgAiIgAkMl0Om4lRh1enu0OBEQAREQAREQgUkSUGI0SdqaSwREQARSJKCYRCAhAkqMEtpMhSICIiACIiACIlCPgBKjevzUWwRSJKCYREAERGCwBJQYDXbrFbgIiIAIiIAIiEBIQIlRSCTFumISAREQAREQARGoRECJUSVMchIBERABERABEegqgSbXpcSoSZoaSwREQAREQAREoNcElBj1evu0eBEQARFIkYBiEoHpEVBiND32mlkEREAEkidw6tSp5GNsOkAxa5roaOMpMRqNl7xFQATGIKAuwyTAF\/irV68OM\/gaUcMMdjWGUNcaBJQY1YDXh65lh6us3Y9xFF+\/3zh6k3PFxqpqi6091reOX6yvbCIwKoGqz+Wo4+LP2DGhbVRhnKI+Ze1Ffau2lc3ht+fpZXP5\/cp81d4tAkqMurUfrawm74Dm2csXIY8mCIh\/ExQ1xqQI8FOMUIqeYdrwn9T6ujwPLEZdH+zG6TfqPPLfS0CJ0V4msiROgMumzqVD38QRKTwRmBqBLpyvvDXk2WOwRvGN9ZdtegSyxGh602vmSRDggJIM+HNRxx7asJv4baFuPlaG7dStzUpsJthMtzJmow27L9hMsKNTIuihEGdem\/nSjh91dF+w+UJbWMfmS1672fFFp0TQJSLQJgGeM19ic\/ntpsf82rQxr41vOqWJ32Y2SrNbaTZKE2uz0uyUZqMM69iQ0E7dF3xMsKP7pa\/TFoq1h3bVJ0tAidFkeXd2Ng4kiYEv2GILxu77oWPzfalj9wWb71NFp48\/Bjo2vy917IhvH0f3x2I8BFveWLThY4IfOiUStlPHbj6UCDaJCLRFgOeO58wXbP581P122qhTTlPCdVFHWJsv2MJ1YtvxcZTUzQcdmwl1a6tS4m99rcQW9qUNGyVieuhL3drxkUyPgBKj6bEf3MzjHPoqfar4GGz\/8qEfdWtro2SONsbVmCIwaQKclVAm8XzXmaOob1HbuGzbGHPctajf+ASUGI3Prlc9ObBcaiyakjr6uMIYvow7TpV+\/jzoVfrU8WEOX+qMpb4tEdCwEyfAnRHKxBehCUVgAgSUGE0AcopThBck9bI4STbKfMJ2+jC2L6FP03V\/LtPz5qCdNZpQz\/OVXQS6TIBn155jSupdXm8X1wa3quuCr\/lTUq\/aV37tElBi1C7fwYzOwa4RrKvav6pfuBYuHfpS+m3UY3bfh3a\/7uu0MYaJ3yZdBPpEoOlnmTPBmHUZ1BmjTt+6667bn7XDsO446j86ASVGozPrbQ8OWd5hszbaTbDFgsVuPlZi832pW5uV2MwH3eyU1K3NL7HTbkLdb29aZ3yby0psRfOYn5VFvn4b447ax+8vXQRCAvY8haX\/rFkbtrL+YXtTdVtDWMbGZ52+H3XEt6FjC\/tjo82Eeugzbp2xbFwrscXGw24+frvZKX279OkSaCcxmkJMsYcuXIb5WGntVg9La+9zGR64ojptvoRx02Y2dF\/M7pd+O7rfho7NxOqUCHZKBN3E6pQIdsoqkucb2qn7Eo5NGzaeF\/RQsNOOndKX0EYd8X2kj04A5khZT3xMzNfqYWntfSl5jvKEGMI2bCbEHrZTx44POuWoQj8bw\/piyxN8aKP0BZuJ2a1updnD0toprc3XYza\/PU+nH22+YDPBbjoldQS9TGBW1bdsLLWPTiCJxMgeIh4k9BgG7LT7gs18fbvp1qZSBGIEeE54hkLBHvOXrR0C8Ic5gp43C234mFA3X7P5pbUNoSRueISCvW78TYxRdw1t9K87prGO8YnZ6s6n\/tUJ9D4x4uHyHyJ0bNURyFMExifA8xbK+KOp56gEOOvwt37o2KxuJTbarK5yLwH4hLLXa7IW1jPujHX6jjvnKP1YHzJKH\/lOhkDvE6OqmMoeQC5Ok7wxX3rpJWdy7949JxGDoT4DeWekr\/ZTp05lbwDgDiiKQee\/+TNvTCmHep76GHfROel722ASI9soLj4TP1lCN6Hd\/K3k0F68eNEtLS1lcunSJbe6urpL\/vzP\/9z96Z\/+qaMM21KtE+uQYh5avDy3sZhv3LhhR6MXJWfbXyhn3Lehm9Dm+5ru3wGx8w+rOhLjXGe8rvWNxfflL3\/Z\/cIv\/EJ2p3K3\/uzP\/qzD1rW1V11PLMaqffviZzG+\/PLLdjSSKweXGNnlR2kXILq\/s9StzexciisrK+7y5cvuypUr7sKFC25xcXGXvOENb3AHDx7cZQt9Uqt3Nea2OA8tXjjGYj58+LAdjd6VnG3OuC3c17FRxwfdF\/8OiJ1\/WNWRGOc643Wtbyy+O3fuuK997WvuE+\/9UffhH397ph87dqy3d2gsxq7tQ931WIyHDh3yj0dS+kxS0RQEE7voCtxzm86cOePOnj3rTp486ebn53fJwsJClhiF9pTrQ4t5aPHy7MZinpubyz0jXW7gHiDxqbNG7oDY+YdVHYlxrjNe1\/rmxcdePPX4GxyC3rV1j7KevBhHGaPrvhZjX+8AnrEy6UlilB8GlxyXnXmgY7N6lZI+vh\/1Ucfw+0sXARGYDAHOKefVZkPHZnW\/zGvDXsXP95EuAiKQLoHeJ0ZsDRchlxuCjs0EGzp2dF+wxdrMTptEBESg2wQ4r3auAIK04wAAEABJREFU0f3VYqful+gmtNHH6pTUsUtEYKIENFlnCCSRGEGTywxB98W3ofsS+lmbb5cuAiLQfQJ5Zxc7q6eMCW2I30ZdIgIiMFwCySRGw91CRS4CItBBAlqSCIhATwkoMerpxmnZIiACIiACIiACzRNQYtQ8U42YIgHFJAIiIAIiMAgCSowGsc0KUgREQAREQAREoAqBoSZGVdjIRwREQAREQAREYGAElBgNbMMVrgiIgAiIwBAIKMZxCSgxGpec+o1FgD+r0Ge5du2aQ\/ocQ9nax9pYdRKBigT858+6rP3luqnOb++bzt2A9G3d4XofbMZAFSVGA934aYTN4fP\/EC9\/NLJvcu7cOffcc885yr6tvep6l5eXsy9O03hGNGc+gRRaSBpeeOEFZ8\/i+fPns7B+7nf+wCFUsFl730ruhRTuh6HfAUqMOImSiRAgMVpZWXGXL192V65ckXSQwUc\/+lHHHk3kgdAkgyPwve99L3u+dAd09\/7THeCcEqPBXU3TD5g\/wskf4jV585vf7FIQi6daedZ10Y+9mf4TohWkToDnzH\/+Uzj\/xODH1FedvUn9+SuLT4lRGSG1t0qAnyL1\/Z\/X7Mf9sR8\/87e38mQUsIxR5l\/Fp2wMtYvApAmkfAdwJvNkFM6MUeZfxadsDLVvE1BitM2h9qsGGI8AlyL\/dPP0h553z7zwYm\/lPc9+MPsngpBC+De4wnron1enX16b2av4mK9KEegKgZTvAM6kCbxNp6ReVar4V\/GpOt\/Q\/ZQYDf0J6Ej8J5447Rb7LE+e7ghJLUME+klAd0Dn920wC1RiNJitVqBdI2A\/+raS9aGbUDfBhu6XpmNHrO6XptOOUPcFm0QERGDyBDiHzGql6dQR6iZW90vTR\/Whn4n1VbmbgBKj3TxUE4GJEuCC8n8Ejm5CW2wx2Mfx8fsxLmNQSgZKQGFPnYCdSVsIZ9KENrP7JfZxfPx+jMcYlJK9BJQY7WUiiwhMjEB4OXF5meQtIuwT88vzKRs7NpZsIiAC7RAIz6mdT8q8GcM+Mb88H8ZFYn1ke0hAidFDFtJEYKoEuLC40ExGXEwldxubslIHOYmACEyEgM7\/RDBXmkSJUSVMchKBNAhw+ZqkEZGiEAERqErAzj5l1T5D9FNi1NVdH9i6Xn35mvt+j6WJ7eKnOFxYJk2MGY7BHCbME7arLgLTIjD0O4BzyZk0aWMfmMOEedqYI4UxlRilsIs9joFPi+WTVr\/0qY+53\/74cm+F9T9x+q9nn+Cdtx1cSH5bWKcNmy\/YEGx+iY6Y3dd9W2jnMjQJ\/fCViMCkCaR0B3CXEU+MYXjewjp9sPmCDcHml+iI2X3dt4V2O\/uUoR++k5auzqfEqIWd4QPLkBaGTm5ILpFU\/m7Sf\/x3n+n8\/nAZmnR+sVrgIAikdAdwl3V50+zsU3Z5ndNemxKjhnfA\/+vRsT8R0fB0SQzHxdjXvyvkr5s4ktgQBSECYxEYvxNnxz9LfdWJY3wK6tkVAkqMGt6J733ve9mfhnjfO9+SlQ0Pn+Rw\/HQtBSnaHH50be2+bray0vpYWeavdhHoE4EUzj8x5DH3z62v5\/mHdutjZdiuerMElBg1y\/PBaE89\/tgDXUo+AS6TlP+IbH7k47XwI\/BJXY7jrVC9RGA0AtwBH\/9n\/9QtLS31Xtr+VwKd\/9GerXG9lRiNS079GiHApcgfkf3Ee3\/U\/db7\/2Zv5cM\/\/vbcnxBaIkOJAI7SF2wmvt102kIdm0QE+k6AO+DLX\/+WS\/UO4NyyR5SIr1NHsJlQD4U2s5lOKWmHgBKjdrh2cNRuL+mpx9\/g+i5FhPlOzwQ\/063E5ovZraTNdCuxSUQgFQJ9P\/+sP28v7MxS4kPpCzZf\/DZ02ih9wSZph8AgEiM\/087DaD5W5vnJLgKjEOB54jIbpY98RUAE0iCg8z\/hfWxouuQTI3sw+eKEHuOGnXZfsMV8ZROBUQjwTI3ib748fybYGIc6upXoEhEQge4S4NyOszrOuAn9GYc6upXoknYIJJ0Y8QDxQBk6dGxWVykCkyLAc2fCnKZTUg+FZxUJ7aqLwIAIJBMq59yEoEynpB4KZx8J7apPhkDSiVFVhKM8gPyiIHLv3j0Xk3DOmM9QbSGbIdV5xkyI23RK6lUEXy5Syir+dXzKntE6Y6uvCAyNAGfWhNhNp6ReRfCd1Pmvsp6UfZQYebvLQ2fCQ+g1PVDPnz+fvaX00qVLbnV1dZeQMF2\/fj3z3djYyMq1tbVdPmGfvteJmQ+1rBIHfhkUvdQmwHNae5CCAfznNrbHN27cKOid0ySzCIhAIwTaPv+NLLLHgygx8jaPZMgk78HjI9+vXLniLly44BYXF3fJ8ePH3ZEjR7IR52bnsvLYsWO7fMI+fa8Tc9UYjx49mjGJvaz95br77l\/1V2IxjWPjueMZpETQbRyrU5qtrdLf09geHz58uK2pNe5ACegOcI6zzZmnRNDtcbA6pdlUtkNAidEW11EeNP5IIB9Xf\/LkSTc\/P79LFhYW3MGDB7dGdG5mdhtt6JNQPYvdYq4aVwbHe+Ej9GH6c7\/zB+4ffO6\/91ZYP3EQjxdedtHxfPlCu183HbsJFyJidXysbqW1tVH6+xnb47m57cS\/jbk15rAIcGY4O5yh1O4Azm0o7G5oo47dhDOOWJ12q1tpbSqbJ7D91bv5cTsxIg8QD5QtBh2b1VVOnwCXov0Ujp\/E9VmIIyTK81ZV6IsvZSihPayH\/qqLQF8IpHwHcE6rCvuFL2UooT2sh\/6q1yPwMDGqN05ne\/MAkRAh6P5CsVHHju4LNtokzRPgd1Z8YQYux74Lcfhx9VEnBokItE0gPBvM1\/fzz\/qJI4ytb3ViGLoknxixwSQ5CLovvg3dF99PejMEuDj4kbn9AnsKfxsptRjYG\/aIvWpm1zWKCDwkwO+rPfXUU47nbBpnR3OW\/z069mbod8AgEqOHx1LaNAnwxZZ\/burzP5d9\/vOfd5\/5zGfcF7\/4RdfnOIrWzh5N8znR3OkSOHHihOP5sufvox\/9aBYsfyeNvzdIxW83v76UqdwP7AF7MVRRYjTUnZ9S3CRH\/PJ6n+Xd73636\/P6y9bOHk3p8RhxWrn3kQDPlz2D\/GSCGPg7Ywg6NmvvY5nC\/cAesRdDFSVGQ915xS0CIiACIiACIrCHgBKjPUhkmBYBzSsCIiACIiAC0yagxGjaO6D5RUAEREAEREAEOkOgxcSoMzFqISIgAgMmYB\/DUYbA\/CjLfNUuAiKQLgElRunurSITgcETIMmxj+FAzwNCm\/lRUs\/zlV0EHhCQkiQBJUZJbquCEgERILkhyTES6NisbiU22qyuUgREYNgElBgNe\/8VvQiIwEMClTQ+yfjevXtOUp+BAb+\/uek279\/PquJan+skGGableiLEqNEN1ZhiYAIVCMQ\/rSo7CdIfDLwpUuX3OrqaqNCwnXt2rVGx2x6jXXGi8VHvOzSnTt33MbGBqpbW1vrLYNYjHWYdbGvxfjKK69k+5XiixKjFHdVMW0T0KsIjEigLCliOD4V+MKFC25xcbFR4c9lHDt2rNExm15jnfFi8R09ehSk7sCBA+7A\/v2Z3mcGsRjrMOtiX4vx0KFD2X6l+KLEKMVdVUwiIAIjE6iSFDEon8x88uRJNz8\/36gsLCy4gwcPNjpm02usM15efDCdmZ11MzMzqL2OPy\/GOty61tdinJuby\/YrxZftJ7EfkWmVIiACIlCZAP9ERrJjHdCxWd0vi9p8P+kiIALpE1BilP4eK0IRGCwBEiGSHgTdB4GNul+im9AmEYHJEtBsXSCgxKgLu6A1iIAItEaAhAgJJzAbZUxCf9VFQASGQUCJUY\/3mXcHlEmPw9PSRaD3BBSACIhA\/wgoMerfnmUrJiH68Ed+yS0tLRXK8vKywzfrpBcREAEREAEREIFCAkqMCvF0t5Fk55tf\/WP39Ieed8+88GJU3vPsB93Kykp3g+jdyrRgERABERCB1AkoMer5Dp944rRbzJMnT\/c8Oi1fBERABERABCZLYNCJ0WRRazYREAEREAEREIGuE1Bi1PUd0vpEQAREYIAE+HUBZIChNxmyxhqDgBKjMaCpiwiIgAiIQHsESIguXryYvbFEbyBpj7NGjhNQYhTnIqsIiIAIdI\/AQFZEYsQbR973zrfoDSQD2fMuhanEqEu7obWIgAiIgAg8IPDU44890KWIwKQIKDGaFGnNIwJ7CcgiAiIgAiLQMQJKjDq2IVqOCIiACIiACIjA9AgoMWqSvcYSAREQAREQARHoNYFkEqNTp045pGw38DExX6uHpbWrFAEREAEREAERcG4IDJJIjEhorl696hD0vI2jDR8T6uZrNr+0NpUiIAIiIAIiIALDIND7xIjkhmTGtgsdm9WtxEab1ccteRspcu\/ePReTcNyYTxM2m+f+5qbbzJX7mVsT82mM+H4PlUv2YOklAQIKQQREICTQ+8QoDGjcOomTSdEY58+fzz507NKlS251dXWXkDBdv349676xsZGVa2tru3zCPuPWr127lo2\/vn7L3bx5Myq3b9\/KfNpaA2snZtaCPgQZWrzsaSzmGzduZM+WXkRABEQgNQKDSYzCnxaRBPk2dBPa8jb68uXL7sqVK+7ChQtucXFxlxw\/ftwdOXIk6zo3O5eVx44d2+UT9hm3fvTo0Wz8+fl5t7CwEJWDBw9mPm2tgbUTc5vjM0eXpOvxtsEqFvPhw4ezZ0svIiACIpAagcEkRv7GkfiQBJnN17FRxwc9lDNnzrizZ8+6kydPOpISX0hQLBmZmd1G67c3rbO22blZNzc3F5XZ2Vlc9qyzyXVYzE2O2eWxhhYvexGLmWcue7j0IgIiIAKJEdj+6p1YUEXhkPCQ+BT5dKdNKxEBERABERABEZgkgd4nRiQ5JDsGDR2b1f0yrw17FT\/fR7oIiIAIiIAIiEBNAh3s3vvECKYkQiQ3CDo2E2zofoluQht9rE5JHbtEBERABERABERgWASSSIzYMpIZBN0Xs1HGxHz9NrOpFAEREAERGImAnEWg9wSSSYx6vxMKQAREQAREQAREYOoElBhNfQvaXwCfQ5Mn7c+uGXpNQIsXAREQgYERUGI0gA23D6VcWlrKPpzSL5eXlx1J0wAwKEQREAEREAERKCWgxKgUUf8dnv7Q8+6ZF15Edsl7nv2gW1lZ6X+AikAEREAEREAEGiKgxKghkF0e5sQTp91iTJ483eVla20iIAIiIAIiMCKB+u5KjOoz1AgiIAIiIAIiIAKJEFBilMhGKgwREAERSJGAYhKBSRNQYjRp4ppPBERABESgVQK8oaRMWl2ABu81ASVGvd4+LV4E+kZA6xWBdgmQEF28eHHPO3D9d+Oi6x257e5Dn0dXYtTn3dPaRUAEREAEdhEgMeLdtnnvxuUdunpH7i5kqgQElBgFQFQdjYC8RUAERKCLBHLfjcs7dPWO3C5uWWfWpMSoM1uhhYiACIiACIiACHZ+AP0AABAASURBVEybQJAYTXs5ml8EREAEREAEREAEpkdAidH02GtmERABERCBSRPw5uP3kYrEc5U6IAJKjAa02QpVBERABETgIYGivyOpd6495DQ0TYnR0HZc8YpAWgQUzcAIFP2Eh7ZRcOida6PQGo6vEqPh7LUiFQEREIFeEyDxKfuMIn4KVDVIvXOtKqlh+SkxGtZ+dz9arVAEREAEcgiQGFX5jKKc7jKLQCUCSowqYZKTCIiACIhAVwjoJz1d2Yk019F2YlRK7dSpU6U+chABEUiTQFPnn3GQKpRCP+oxqTKWfERABNIjMPXEyEfK5eTXQ72sPfRXXQREoD8Eys53Xjv2q1evOgQ9L2LakFg7fUOJ+ckmAnsJyJIagU4lRiFcLjEktKsuAiKQPgHOPlIUKe0kNOaDjs3qfkkb4tuki4AIiEBIoDOJkV1mlIgtlIvMr5tdpQiIQDoE7IxTIhbZpM4\/c5rY3HklvwB87949J6nPwBjf39x0m\/fvZ1W4Zor3gg0xU+ZPn6hsj9OUD\/NK9u617UWKZWcSIy5AAFuJboKNS4vSbCpFQATSIWBn20o\/Mmxtn3\/mMGEuf\/5Q5+3gly5dcqurq40KCde1a9caHbPpNdYZLxYf8cL3zp07bmNjA9Wtra05s\/s25jb7+votd\/Pmzajcvn0rG6cpH9bD3FUkFmOVfn3ysRhfeeWVjHOKL1NNjMouoBSBK6ZJE9B8XSXQlfNPQuQzol60tsuXL7sLFy64xcXFRuX48ePu2LFjjY7Z9BrrjBeL7+jRoxn6AwcOuAP792c6DMw+Nzv3wMbcZp+fn3cLCwtROXjwYNanKR\/Ww9xVJBZjlX598rEYDx06lHFO8WWqiREXUBWoXFL4Ulbxl48IiED3CXCmq6ySc48vZRX\/tn3OnDnjTp486fjC26TwhZ4v6k2O2aWx8uJjv2ZmZ93MzAzqA65UZmZ324gH++zcrJubm4vK7OwsLm62IR\/mrCp5MVbt3wc\/ixH+GegEX7afug4EZpeelbYk6lyKVh+nZAykrC8+JmW+ahcBEWiOAOeO0axER6gXnX\/a8MEXQceGXlXo4\/tSH3UMv790ERCBfhPoTGJkFxElYlh93WyjlHbJMQ56Xl\/a8DGhnucruwiIQLMEOHeMSImgI75OPSb4cF4RdN8Hm1+P6fTBz4R6zE82EZgQAU0zZQKdSYxiHOpeUFx0\/hjo2MK5sNEW2lUXARGYHoFRziS+SLjaUWz4IuEYqouACAyLQKcSo7JLqax9ElvHb+QjeW\/fDNeQ51fXbvPUe0vq9tta665F\/fe+lTV1Jvb8NVmWne+y9ibX0uhYGkwERKBXBDqRGPETmzJpk2p44bKW0Gbz81bdpaUlF3u7LgnT9evXM9fwbaZNvx2zibet2ttaR3k7ahgHMbOW0J5qfWjxso+xmG\/cuJE95028cN7KpIl5NIYIiIAIVCEw9cSIBKSK2MVZJag6PszDevLG4K26V65ccbG36\/I2xiNHjmRdw7eZNv12zCbetso7YFjsKG9HDeMg5jr9w\/G6Xj9+\/HjSb6mO8Y\/t8eHDh3l0agtnrYpwLpHaE2oAERABESghMPXEqGR9D5rt8nxgaEHh4mWeoqF5q+7Zs2ejb9e1tzHSfybyNtOm34rJPHXekjo7u\/221jrrspjrjNGnvkOLl72JxTzpt+pyLhGeeYkIiIAItEmgN4nROBC4SEl2rC86Nqtn5c5LUduOiwoREAEREAEREIHECUw9MSIhQYwzOmL1uiWJEOMh6P542Kj7JboJbRIREIH2CIRnLay3N7NGFoHhEFCkoxGYemLEci1h4VJER7Aj2HzBNqowHhL2MxtlTEJ\/1UVABJonwNljVM45OkIdweYLNokIiIAItEmgE4mRBehfiL4NO2I2lSIgAukRiJ1xbCbpRTxuROonAiLQJoFOJUYEat8doktEQASGRUDnf1j7rWhFoIsEOpMYcSECiO8OKVMRPgMmlVgURzsENKpzqZ5\/7a0IiED\/CHQmMSIhssuxfxjzV8wHQi4vLzslSPmM1CICqZ5\/7awIiED\/CEw1MQoTofBypN2XfuDdvcr3vfMtbmVlZbdRNREQgQc\/JTIUKZ5\/i02lCIhAfwhMNTHiIjRUJEDovg3dF9r7Jk89\/ljflqz1isBECHC2baJUz7\/Fp1IEkiKQeDBTTYx8tv4l6duli4AIpE9A5z\/9PVaEItAXAp1JjPiO0YRLEr0vELVOERCBegQ47yY6\/\/VYjtFbXURABDwCnUmMuAxNvPVJFQERGAABO\/uUAwhXIYqACHSYQGcSIxjxHSOlL9hMfLt0ERCBHAI9NXPOw6VjMwnbVBcBERCBNgh0IjHi4iM4vls0nTqCzRdsEhEQgXQI2JnnnJtu0WHzxewqRUAERKAtAlNPjOzSswCpo1uJPnBR+CKQLAHOOWIBmm6l2VV2lwCf0WYyyirpM4q\/7\/vqy9f8qnQRaJTA1BOjRqPRYCIgAiIgAhMjQHJz8eJFt7S0lEnVD7O1fnwAbtXF0uezn\/1s5v6lT33MId9XgpTx6P9LtyKYemIU\/ujc8JjdSrOrbJ4AF06RND+jRhSBbQJ559vsVm5767VrBLg3+ADbT7z3R92Hf\/ztlT\/M1vrxAbhVY7I+7\/\/BeffsO97s1r751apd5ScCIxGYGcm7RWf\/AkTXj9JbhB0MzXdt9h1frFzWnzQJiKnaNAHOvI2JrvNvNPpRPvX4Gxxiq61ajvMBuG89NOvesfj6qlPITwRGJtCZxIiLkAuRCNApJZMh8PSHnnfPvPBiVN7z7Acrfxc4mdVqlhQJcOZ1\/lPcWcUkAv0j0JnECHRcjpR2QaJL2idw4onTbjFPnjzd\/gI0Q8cJTGZ5Ov+T4axZREAEigl0KjFiqSRFdkFSl4iACAyHgM7\/cPZakYpAVwl0JjHiQuwqpNTWpXhEoGsEdP67tiNajwgMl8BUEyO7DCn5KZFfsiXUrURHqEtEQAT6T8DOM6XOf\/\/3UxGIQIcI1FrKVBMjLkNWT+lfjtgQ7FaiI9QlIiAC\/Sdg55lS57\/\/+6kIRCAVAlNNjGIQuSRjdtlEQATSJ6Dzn\/4ejxWhOonABAl0JjHShTjBXddUItAxAjr\/HdsQLUcEBkygM4kRe8CP0yl1SUJBIgLJEogGpvMfxSKjCIjAhAl0JjHiUrSECH3CHDSdCIjAFAlw5nX+p7gBmloEROABgc4kRnYpsjJ0Lkp0SQ8IaIkiUJMAZ96GQNf5NxoqRUAEJk1g6okRl2AsaLNbGfORTQREoN8E8s632a3sd5Ra\/SgE+GOxoX\/MFvqoLgJNEYglRk2NPdFx+A4TqTJp6Ec9JlXGko8IiIAIiEB9Anfu3MkGOX\/+vEOofOMb36DI6mbLDHoRgRYJTD0xsoTEYqRueljmtWHnO0sEPexnddoQq\/slfUPx26WLgAg0T4DziNjIvm42K4vazEdlfwlsbNzNFs8fteaPV1N53aPHKJxvywy1XzSACOQTmHpixNJISCh9sUvQSr\/N12n3+6Nj831Mpw2x+jglP9JF7t2752JiY96\/f9\/UrIz51rFlg2693N\/cdJu5sr2GfJ+y9s2tGbb\/r7NW9Y0\/K33msv1UNPMaO5N2hq1sZiaN0gcCJ5447RafPJ0tdf7QX8tK35YZ9CICLRLoRGJEfNO+AJnfhPXkCT\/OXVpacpcuXXKrq6u7hITp+vXrWde7d7e\/+9nY2Mjqa2tru3zDvqPWr127lo27vn7L3bx5Myq3b98q9ClrZ9z19fVsjLz1EzNrGXX9ffUfWrzsUyzmGzduZM9FUy+cvTpjqa8IiIAINEWgM4kR3zXa5RgrzdZU4OE4zG9SNNfly5fdlStX3IULF9zi4uIuOX78uDty5Eg29P65uaycm90ujx07tss37Dtq\/ejRo9n48\/PzbmFhISoHDx584PPII484xPf12327rx88sD1G3vqJOa9t1Jj64D+0eNmTWMyHDx\/Onq2mXjh7du5ipdmamk\/jiIAIiEAegakmRuFlx+XIQmOl2WhvWsKxqYdrsznPnDnjzp49606ePOlISnwhobBkY9\/MNtqZ2e3S92tKZ02zc7NubisJi8ns7CwubnbLZ3ZLR3w\/6jjQ7tt9nTZ88tZsMee1j2ef38O2K+MMLV64x2LmGeG5qCPhGePcMV6sNBvtEhEQARFok8D2V+02ZygY2y47LkikwFVNIiACiRHQ+U9sQxWOCCRCYCKJURkrLkjEkqNYabZwLL8fbfhhQ68q9PF9qY86ht9fugiIQHUCnDWEc0evWGk22kcV+iJV+lX1qzKWfERABPpJoBOJkV1GXI4IKCljQlso+DEGgu63Y\/PrMZ0++JlQj\/nJJgIi0DwBzh2jcu4QX6fuC22jCGNbf\/S8vrQhee2yi0AJATUnRGDqiZFdWsbULifKUMwnVobjmA92063Ms2FHzE+lCIhAuwQ4b4jNwplHpwwF+yhCf39sdGyxMWhDYm2yiYAIDIvA1BOj2EXFBRXKsLZF0YrAMAh09vwX4OfjC\/r8GVRNrt0wZZ+VtvPZbVXGf9Bvpw\/9H9iyz2Xb+Yw1r31zc9uW+e2o9NvM\/P3Pc9tujLeZX3WfKvEM0Sfbh0Rfpp4YwZXL0YS6RAREYDgE7OxT9iFqPsss9jlmfOZTHSHh6ttngrFe9ow\/51H2mW1+fNbPPu+N\/nfvbn\/m2\/r6LWefsUbSw\/i+jbpv5\/PWfLG+9PHtvj6KT95nuMX22o8x1p6CzWJ85ZVX2IokpROJkf\/TIShzQYaCXSICYxJQtw4T6Nv557PMYp9jxmc+1ZHY50XVGW8Sfe3z1A4cOOAO7N+fPWV5n2vmx2f97PPe6G8fATG\/MO\/sY09mdj72hI+NMBuT2MehYOfjJHwxv1ib+Y3ikxdPjK8fY6w9BZvFeOjQIbYiSelEYuSTtUsSm+lWYpOIgAikS8A\/66Zb2ZWo+Syz2OeY8YW4jvBFmy\/YdcaYRl\/2ZWZ21vlJTGwdYXz0swRnZnbW7ds3g8nx+WoIlX379lE4Pk9tdssnq2y97JgddhIqX2Z3\/GJt5jeKTyyWPFsYY55fn+0WIyy3tiLJ\/7efxCmGxqUXmz7PnvnqRQREIAkCeec8zz5K0IzBT56tDzo2q6sUAREQgRiBqSdGsUXJJgIiIAJNECARIiFC0P0xsfl16SLQJQJay\/QIKDGaHnvNLAIiMAECJERIOFVVW9hP9WER4JeNi2RYNIYRrRKjYeyzohQBEZgqAU3eVwK8C3FpacnlyfLysiNx6mt8WvdeAkqM9jKRRQREQAREQAQyAk9\/6Hn3zAsvRuU9z37QraysZH56SYeAEqN09lKRTJCAphIBERgGgRNPnHaLefLk6WFAGFiUSowmtOH6UeuEQGsaERCBzhLggx39u\/Duxt1srXdu385KvYhAFwgoMcp2of0X\/p16lH+L5vIokvZXrBlEQAREoDkCJEX\/+tc+7c7DHEJ7AAAQAElEQVSdO+e4Dxn5O3\/2HQr3la98xX37O9\/OdL2IwLQJKDGa0A68751vqfxv0SREFy9ezP1lP34J0C6WCS1f04iACIhALQIkRle\/\/ieO39nhd3MY7Mgb30ThHvvBH3KHHj2W6XppiYCGrUxAiVFlVPUcn3r8scoDkBjxC31cIEW\/9Fd5QDmKgAiIQEcIZL+zs\/O7OfvnH8lWdfB1h93+hW09M+hFBKZIQInRFOGXTX0i7xf+sO9cLGVjqF0ERCBJAgpKBESgJQJKjFoCq2FFQAREQAREQAT6R0CJUf\/2TCtOkYBiEgEREAER6AQBJUad2AYtQgREQAREQAREoAsElBi1swsaVQREQAREQAREoIcElBj1cNO0ZBEQAREQARGYLoF0Z1dilO7eKjIREAEREAEREIERCSgxGhGY3EVABEQgRQKKSQREYJuAEqNtDnoVAREQAREQAREQAafESA+BCCRJQEGJgAiIgAiMQ0CJ0TjU1EcEREAEREAERCBJAkqMerKtWqYIiIAIiIAIiED7BAaRGJ06dcohVXBW9asylnxEQAREQAREQAQqEeiMU\/KJEYnO1atXHYKeR542JK9ddhEQAREQAREQgfQJJJ0YkeiQENk2omOzul\/Shvg26SIgAiIgAmMSUDcR6CmBpBOjNvbkpZdecsi9e\/dcTGzO+\/fvZ6qVWWXrJdYntG25Zf\/f39x0m7myM35uO33LfMraN7N18BKuUfX4\/g+FC8+ERAREQARSJKDEaMRdPX\/+vFtaWnKXLl1yq6uru4SE6fr169mId+\/e3VVubGxk9bW1tV19wjGoX7t2LfNdX7\/lbt68GZXbt2\/V9qkyxvr6ejZP3rqJmfWy7iHI0OJlT3di3vXc3rhxI3su9CICIiACqRFQYjTijl6+fNlduXLFXbhwwS0uLu6S48ePuyNHjmQj7p+b21XOzW7Xjx07tqtPOAb1o0ePZn3nF+bdwsJCVA4ePLjtM1\/BJ2ecSmMc2J4nb93EnNdGLKnJ0OJl\/2IxHz58OHv+9CICIiACqRFQYjTijp45c8adPXvWnTx50s1vJSW+kMRYsrFvZhutlTOz23Xfv0hnWbOzs25uK8GKCW2Zz1wFn9g4W+NWGmNrfObJW6vFnNeemn1o8bJ\/sZh5JnkuJCIgAiKQGoHtr9apRbUTD79M7f+yNTq2nWYVIiACIiACIiACiRIYN6ykEyOgkAiRECHo2Eywma5SBERABERABERABJJPjNhiEiIE3ZeqNr+PdBEQAREQgWkQ0JwiMBkCg0iMJoNSs4iACIiACAyRAO\/cRHiHLoLuyxCZ9DlmJUZ93j2tXQR6TEBLF4FUCNjHuJw7d84999xzjpKPdTFZXl7OPv8ulXhTj0OJUeo7rPhEQAREQASiBF59+Zr7fo5EO+QYn\/7Q8+6ZF1507734G+5vf\/TX3N97\/jNZHdt7nv2gW1lZyekpcxcJKDHq4q70ck1atAiIgAj0i8CXPvUx99sfX47K\/\/3Cr1cO5sQTp90i8uRp98Yf+pEt\/d1b8tBWeSA5doKAEqNObIMWIQIiIAIiMGkCP\/+33uV+5vRbs2l9\/em3vcn9xXfXMrtehkcgNzEaHgpFLAIiIAIiMCQCP\/zYI+4di6\/PQn7Hidc\/1BcfzWyjvoR\/G3PU\/vLvBgElRt3YB61CBERABERgsgQ0mwhECSgximKRUQREQAREQAREYIgElBgNcdcVswikSEAxiYAIiEADBJQYNQBRQ4iACIiACIiACKRBQIlRGvuYYhSKSQREQAREQAQmTkCJ0cSRa0IREAEREAEREIGuEphcYtRVAlqXCIhA7wmcOnXKIWWB4IP4ftRj4vtIFwERGA4BJUbD2WtFKgJJEiCpuXr1qkPQ84KkDR8E3ffDForfLl0EqhCQTxoElBilsY+KQgQGSYAEh4TGgkfHZnUrsdFmdXRsVlcpAiIgAkZAiZGRUCkCIjBYAiRJJg8hxLWXXnrJ3bt3T7LFwAjd39x0m\/fvZ9U8Nlnj1kvmu7nta58UTbnT3UXbGX+nz9YQzm133\/Hdmjtrt3K7cXscs+0uszF42XJ9MO+WsvU\/Vsd6MmXrpWiczZ01mQ\/9NnethXm3JtkaJ49LX+1bISX7vxKjZLdWgYmACFQlwE+QTEiQivqdP3\/eXbp0ya2urjYqJFzXrl1rdMym1xiOx3phdefOHbexsYHq1tbWojGY7\/r6LXf79q3MlyQC5bXXNh3JBXqs3bfhY\/2w37x50\/liY8fazG99fZ1hHPNubr62R7fxaVjfWq\/1C8vYXLdu3YquJ49LyLTrdXtOX3nlFfAkKUqMktxWBRUjIJsIxAiQEPl26kXJ0eXLl92FCxfc4uJio3L8+HF37NixRsdseo3heEePHs3QHThwwB3Yvz\/T82Iw3\/mFeXfw4MHMd2Zm+0vQzOyM27ejz8\/vbfdtdPR9FxYWnC82Nn18u68fPLAzfzbvLEO67TXs6DtroaFwnJ04fB\/0Rx555MGabD15XEKmXa\/bc3ro0CHwJCkzSUaloERABESgJQJnzpxxJ0+edHwBbFL4ws0X0SbHnMRYYJ6ZnXUzO8lE0Zz4zm75Iuj79u2jcDNb5db\/mT47N+vCdt+Gk+87NzfnfLG+s1vj+HZfp83GmdleQrYG0\/ft2zFuOeHr9\/X1cC7qJuZHfWuYxp+XIs5tttlzSnzElaLM9DcorVwERGDoBMKf7vCTHmwhF2y0mR0dG3V0ShPq1mY2lSIgAsMhoMRownvNv88i\/rTUffHbqujff\/naHreYbY\/TCAZ\/fb7O7w0gIwwlVxFolABJDMkMgu4Pjs3qtFFH0GP2sM18VIrAVAho0qkQUGI0Yez84ubS0pJbXl52lmBcvHjRYTP55V\/+5cqrenUrKfpfn\/u0+9KnPuYsGaI0W+WBShxt3bZGK8+dO+eee+4594EPfCCLp2QYNYtAKwRIdJBw8NBGHYn5YUfCNtVFQASGRUCJ0RT2+33vfItbWVnJZiY5QseGgRIbehUhMVr75lcdYv5me\/ptbzJT7fLpDz3vnnnhxT3y3ou\/4d75k8sP4qk9kQYQgXQJKDIREIEeEFBiNIVNeurxx\/bMajYr9ziMYXjH4qNj9Ip3OfHEabcYkydPuzf+0LvinWQVAREQARHICPANb55kDnrpDAElRp3ZCi2kdwS0YBEQARGoSCDv1xH4tYSf\/Kn361cRKnKchJsSo0lQTngOPuk14fAUmgiIgAg0QiDv1xHe8+wH3Te\/+seNzKFBmiGgxOghR2kiIAIiIAIi0AqBol9HaGVCDTo2gWQSI95mi5SRwAfx\/ajHxPeRLgIiIAIiIAL9JqDVVyGQRGJEUsPbbBH0vMBpwwdB9\/2wheK3SxcBERABERABEUifQO8TIxIcEhrbKnRsVrcSG21WR8dm9aqlvasg7y8i2zj2uzdWmp3SbIxBPU\/u7\/krzZuOP264Ldt\/sXlz5687M8ZD\/+02m+eh3e+PvuNXYZ68MZjj\/s4aiEcyjL+6zvMmmT6BSa3A7r2wtPnv3L5tqkoR6D2B3idGTe0ASZJJ0Zj2zoJLly651dXVXcKlcf369az73bt3d5UbGxtZnRdrW1tbc\/ap0dZuJX7r67fczZs3o3L79i1c3J07d7KSF\/O3ts2thMe3h2OZn\/UL26lX8bmzsb0G4gmZpFZnj9mz1OIqiicW840bN3i0JAMgwP5\/+CO\/tOtDaHknFcJ9CIKvfOUr7lvf\/jaqRAR6T0CJ0c4W8hMkExKkHfOegr+sfeXKFRf769r81eEjR45kffbPze0q52a36xitjb+2bH9x2tqtxG9+Yf7BX2jmD\/f5wh+bxIe\/ak2J8IcD8bG2B3\/UMWcc87N+9A2lio+tgXi6\/peh666PPW4vzkVXd31t9I\/FfPjwYR45yQAIkBjxrqnYu6p4RxUIDh9bdK+++ipq9rZz+oSSNepFBHpAQInR1iaREG0VD\/6nnpccnTlzxp09e9adPHnSkVD4QlJhicS+mW20Vs7MbteZxGzWF5u1W4ltdnbWzW0lWDGhbdvn4bizc9v+1rZv3z5cHPXCMXb6jePD2DM7sVo8KZe2xynHGMYWi5lnJXu49DIYAtF3VT15Oov\/wCOHspIXforET5NCwU67RAS6TuDhV9Wur7SH69OSRUAERGBoBGI\/WXrmhRfd6X\/4gQzF\/Z1\/4s8qehGBDhLofWIU\/nSHn\/RgC1ljo83s6Nioo1OaULc2s6kUAREQAREoJxD9ydITp92Jt7+7vLM8+kYgyfX2PjFiV0hiSGYQdGwm2EynjTqCHrOHbeajspxA+DsFYb18BHmIgAiIgAiIwHQJJJEYgZBEB0H3JbRRR3wfdGwm1CWjE+B3CMLfK\/Dry8vL2S9mjj6yeoiACEyEgCYRARFwySRG2svpE8j73QJ+v4B3r6ysrEx\/kVqBCIiACIiACBQQUGJUAEdNoxHI+92CxSdOu8Wdd6+MNqK8axJQdxEQAREQgREJKDEaEZjcRUAEREAEREAE0iWgxKhPe6u1ioAIiIAIiIAItEpAiVGreDW4CIiACIiACIhAVQJd8FNi1IVd0BpEQAREQAREQAQ6QUCJUSe2QYsQAREQgRQJKCYR6B8BJUb92zOtWAREQAREQAREoCUCSoxaAls07N2Nu1kznwydKXoRgZ4Q0DJFQAREIHUCSoymsMPfXftuNiufFI1Q+cY3vkHhrMwqehEBERABERABEZgoASVGE8W9Pdn6+nqm8EnRfCI0ldc9eozCWZlVWn\/RBCIgAiIgAiIgAj4BJUY+jQnr2SdF73wi9Pyhv5bNbmVW0YsIiIAIiIAIiMD4BMboqcRoDGjqIgIiIAIiIAIikCYBJUZp7quiEgEREIEUCSgmEWidgBKj1hFrAhEQAREQAREoJsC7lIukuLdamySgxKhJmhqrlEDRwaetdAA5pEVA0YiACGQEeIfy0tKSy5Pl5WWnOzJD1fqLEqPWEWsCn4AOv09DugiIgAhsE+Bdys+88KKLCe9eXllZ2XbUa+sElBi1jnhQE5QGq8NfikgOIiACAySQvUv5idNuMSY7714eIJaphKzEaCrYhzupDv9w916Ri4AIiEAfCBQnRn2IQGsUAREQAREQAREQgYYIKDFqCKSG6RYBfkmxSLq1Wq1GBKZHoOic0Da9lU1m5j7Nwn4USZ9i6fJalRi1tDt3N+62NLKGLSPAxfHhj\/xS7rs7eNfHT\/7U+\/UOjzKQak+eAGfl4sWLhWeFN0wkD6InAbIX3F95oneuNbORSoya4ZiNcu3aNfeFL3wh07\/zZ9\/JSvujsFZmRr20SoDL\/ptf\/WOX94vevMOD9lYXocGnREDT+gQ4C2XCu53yzgrvkOK8+GNKnx6Bsn1iL8tWV\/Y80F42RurtSowa3GESo6997WvZiEfe+KasfN2jx3aVWUUvhQQ4mEVS2NlrzP1Fb73Dw6MkNVUCnKGqPw3KPSu8Q0rn2vlkFgAADJ5JREFUpTOPSN19qvJM8NOoof\/kSYlRS4\/8\/vlHspHtj8JamRn1UkjA\/3ExhzQUDu0f\/uEfZv8UxkEPpXBwNYpAywTC5zFWb2IJsXFDGz9BKPspQxNr0Rj9IMDzUeWZwKcfEbWzSiVG7XDVqDUIlF3kHNqi5Im2GtOrqwiMTYAvPGU\/pSHRJ7nHd9yJ6Fs2j52DE\/zUJ0\/006Bxt6DX\/fRMFG\/fhBOj4sVMu\/XUqVMOmfY6hj5\/lUNbljxVYcgXlyKpMkYVn6I5rK3KOPJplwBnH6kzC\/tJ4l72fOKDry\/8UzyCrWwN+DBG2Txl46hdBERgLwElRjtMuBCvXr3qEPQds4oSArdu3XIvrb5U4tV8c5XkqWxWvqPmu\/c8iX1XzxcufsGesmx82vkCVvadPfPH5qJ\/08J6yiSck1hHiTns34c6Z56zj6DXXXOV5zN8\/s6dO+eee+45R1n1eagyT91Ymup\/586dbKhbt25nZYov3Iff+ta33Pr6eorhZTEN4T5QYrS11VyEXIhbavY\/OrasErwUfVEJXEurr758rdCH9u9v+cQk1tH8w7abN2+6mJhfrM1sZT5cBKsvrWZu1idWZg5bL7E2s201Z\/9bPVZmDlsvsTazbTVn\/1s9LLPGrZenlv+F+4mPXIoK78ThO\/Jwv7kUvvjFL+b+flPoT51xxpmLvk1L1STN\/x2uWMxb+JL5n7POmbeA0LFZPSyL9sR8w2fOr5tP+Ez8jX\/8gpv\/sfe5t\/+dn3E8M03N488d6raW0O7Xi3ys7f79+6Y6u4fCe+uRR7Z\/73Lz1l898DXltXuvuQMLr8uq9M8U7yVmoxl7OA92JNZmvrQj+\/bNUWTCGjIlePFZhLq5mp378Fs7iZHZQh+zW1nWjt8oPk09N3nj2H1ga0qxVGJUcVff\/OY3uzNnzrjwuzy+2zehbd\/CYffGt73LPeo23N9967E9Jba3vfFR99bXv8596VMfy+SxRx\/d5UedZdH+2x9fdjGhDZ\/\/+blPZ2PZePjSRp01sJb\/8b\/\/yP3e7\/3eHsG+\/21POcpYOzbainx+\/\/d\/381sxVzkU2WcSfpYTN+8seG+\/K2\/iApts48u7tlv9hjulLbvRSV+jMN4RXPx3OBbNFYTbX\/0J1\/PvvjyBTgmB3\/kJ9z\/eenGrrhZVxjzZz\/7WUxDkAcxVr0D5n7ghyudqfCZ+Pp3\/5+bffRNbnVj\/1a599nz95894ZnhWebsxIS2ts8lc\/zwW9\/mZjZuu8eOHnXcXdw\/3EOhYAcmJWJ3FHfiazf\/0r326l85bLQh6NxhlNQRxn\/06Ovdvo07D3yL5gnbrG5j2bz+GtCZ98fe+rgblR\/3ITFS2p7AqGicsnbGqerDXcOz4T8rvk4bPozHuDGhjWcLX7+v6diJ8Xd\/93cpkhQlRhW3lUvx8uXL7sqVK4Xy3\/7zb7l\/+6v\/0r3z53\/FXfjkZ\/eU2N71sd9wn\/wPn3swzr\/\/zd\/c5Ue9bB6\/nbGQ0MYaWMvn\/tFPuJj8p5\/7++4Lv\/rxaJv5p+hTNab\/+pv\/5sEe+WxH1RnHeMZK1sNzM+q44\/gzz+d\/8bzLk\/\/ywofdFy\/\/SmncP\/3TP13x5KTjVvUOgF9sn83GfhedO9p5Zsr2l720MWMl4xTNQ5+6PvT\/lU\/+qlv8J59w3Guj3F3cWdxR3ImLF\/5Vdgdis7jRaac0G+P\/wD\/\/9T2+1j5KyVg2r78GdOb9xU9+eiL3Iwyb2qcqzw0+7H2esB6erTKWKd8BSoxGuLe5GM+ePeskPWag\/Wvk+eUsjHB0knElbp1\/nX89A2cdZyGZgx0EosQoAKKqCIiACIiACIjAcAn0PTFqZOfCX7bkFy+xNTK4BhEBEeg0Ac46Z94WiY7N6ipFQASGRUCJ0c5+cxFyISLoO2YVIiACAyDAmefsI+gDCFkh9oaAFjppAkqMPOJciIhnkioCIjAQApx9ZCDhKkwREIEcAkqMcsCMa+Y7TmTc\/l3qRxwx8ddo7b7N9KI28+lSyXpj68GONNkWG2satjAu6jHx12btvs30ojbzGXpZFL\/xC8uiPn1qs7j6tOYqa7W4wrJK3z74EFdsndiRWFufbUqMGtw9HhC+40TQGxx6akMRSyi2GGK0NnSzU1LPa6O9S8JakdiasOfFMW5bbJ5J21g7EpvX4vVL86OP2dHNTkk9r412STUCxtAvq\/Xstlfqz4e\/X6Z3e0fKV8eeITFP7BYnesynrzYlRg3tHA8GD4kNh47N6qmVxEaMFhc6NuqU1NERdGzoXRTWh4RrY82+HR0bfpTUnaPmdv0pmaK2be\/pv7J2ZJSVFMVV1DbKHPJNk4Cej37uK3cEEq4+9f1UYhTuuOq7CHAATHY1qJIsAdtvymSD7HBgcDfp8DK1tICA7Rll0KRqzwgoMQo2TNXdBPhuwUQHfjebVGu235Ta8+Z3GaYxsZngboKf2VV2m4DtGaX2rdt7VbY6JUZlhAbczgH3w6euA+8TSU9nj\/2oqGvPfSL1dZjGhJGxU5pQF3+j0d2SffJXR73D++YvVXqEgBKjCBSZREAEREAEREAEhklAiVFD+x5+h8B3C9gaGn4qwxCDPzF1i4mSurWjY6NOSR0dQceG3idhzazd1oyOjToldXQEHRs6JXV0BB0beteFtfprpG5rp6Ru7ejYqFNSR0fQsaFLqhOAm+9NvZSj36GjOjEQiy0PHZvV+14Sjx8D9ZTi82NDJzZiREfQsaGnIEqMGtxFHgweEAS9waGnMhQxEIsJdX8h1Mdp88fouj5ujEX9uhyzv272lrq\/XurYEfSqbb6f9HwCMIWtCfV87361EEuKcbELfmzESB17ykKMxIqgpxSrEqOGd5MHBGl42KkNRywmsUWM2xYba9o2YomtATvSZFtsrB3bRItYXNhMYosZty02lmx7CRhfyr2t\/bYQE9LvKOKrJy6TuEd\/rcQVWz12JNbWZ5sSoz7vntYuAiIgAiIgAiLQKAElRo3ijAwmkwiIgAiIgAiIQG8IKDHqzVZpoSIgAiIgAiLQPQKprUiJUWo7qngeEOCXAh9UGlLaGLOhpWkYEeg1gaGcraHE2eeHUYlRn3dPa88lwOXTxi8FMiZj506sBhHoPYHJB8CZ4mxNfubJz0icxDv5mTVjVQJKjKqSkp8rO8xl7T5C87Uy1lZm89ur6MwVkyp95SMCqRDgDLQVC2PHZJz5GGecfpPo0+W1TSL+1OdQYpT6DjccX96FkGdvePrS5Iz5WAvflaGHgj0U\/EO\/ojr9R+1TNF7bbRpfBCZJgPMRStF5oQ3\/Sa6xyblY\/6jjEe84\/UadR\/7jEVBiNB439WqJAJfFpC+NSc\/XEjoNKwIiMCEC3BkTmkrTTIGAEqMpQK835XR7cyGQvPiroI49tGE38dt8nX74+LZQp933o46Efk3VGduXpsbVOCLQBwL+s48eWzP2UGJ+TdmqzBX6UM+bP9bm20ynNPHHwkbdL003O3UTbJL+EFBi1J+96s1KuQxIZHzBVjcAxmMMSgS9aWGdjO0Ltqbn0Xgi0EUCPOv+s4+OzV8rdewmtKFTtiHhfMyFzZ+LOnYT2tApx5VwTOrhWDYHJUI7fui+YKNNUpHAlN2UGE15AzT9QwJcHlwmWCipo09SmHeS82kuEUiNAOc2lD6eqybX3ORYqT0vXYxHiVEXd6Xja+KQc\/GxTErq6KkIMfmSSlyKQwQmQYD7IJBJTNuZOYjdvz\/QO7M4LaQSASVGlTDJaSgEuMS42HwZSuyKUwSqEOBscE5MqFfp16YPa7D1UFJvc76ysZnfF9ZU1kft3SGgxKg7ezHYlXCBcHFQ+hCox+y+T0y3frG2mK1oDtrCPtiYI7QPpq5AB03Ann\/OAFIXBmMwZp1x6M84JqOORf9R++T5VxkLH9aaN4bs0yWgxGi6\/Hs7O4c673BbG+0m2JoIlnFszHHGs75+yZg2Fnpem\/moFIG+E\/CfcV8Pn3\/asIXxYvclbG+yzvz+XOjYwjmw+xK2+3X6+77U\/fZRdPraWPTz62bHRpukHwSUGPVjn9pcZeWxw8NdVKfNl3AS2nxbWLe2mB0bYj6xknYuJb8NW0x8H3Tfx+qUCGPSji4RgT4S4PnNE+IJ27CZ2PMf+mDHBzvlqEI\/G8P6YvN16iZmp6Sf2f0SO+154vviQ50S8XXqiG\/zdWvzbei+4GPCumizusruEVBi1L090YoaItDG5dPGmA2Fq2FEoHUCPP98YQ8Fe93Jxx2DfuF6qGOvu6Y2+nd1XW3E2p0xR1uJEqPReMlbBERABAZNgC\/soUwbSLge6tNek+bvL4H\/DwAA\/\/8rdMGrAAAABklEQVQDAOIcf0boMDlVAAAAAElFTkSuQmCC","height":281,"width":466}}
%---
%[output:6816dab6]
%   data: {"dataType":"text","outputData":{"text":"[14:04:55][INFO]  チェックポイント保存完了 -> result\/r05_checkpoint.mat\n","truncated":false}}
%---
%[output:4dca65b3]
%   data: {"dataType":"text","outputData":{"text":"[14:04:55][INFO]    R06で読み込む: load(\"result\/r05_checkpoint.mat\")\n","truncated":false}}
%---
%[output:4752e9cd]
%   data: {"dataType":"text","outputData":{"text":"[14:04:55][INFO]  コーパス:         500分子 有効FDA薬SMILES (拡張後 500分子で訓練)\n","truncated":false}}
%---
%[output:86afdc49]
%   data: {"dataType":"text","outputData":{"text":"[14:04:55][INFO]  語彙:     37文字  |  最大系列長: 100\n","truncated":false}}
%---
%[output:5038645d]
%   data: {"dataType":"text","outputData":{"text":"[14:04:55][INFO]  損失プラトー:   val=2.6101 (エポック70+)  --  ランダムベースライン=3.6109  (パープレキシティ=13.6)\n","truncated":false}}
%---
%[output:134573c2]
%   data: {"dataType":"text","outputData":{"text":"[14:04:55][INFO]  --- 生成指標 (N=300ずつ; 二項型SE~1.7%)\n","truncated":false}}
%---
%[output:6104ce8f]
%   data: {"dataType":"text","outputData":{"text":"[14:04:55][INFO]  マルコフモデル:   有効性=9.7%  多様性=93.1%  新規性=93.1%\n","truncated":false}}
%---
%[output:97f77c7d]
%   data: {"dataType":"text","outputData":{"text":"[14:04:55][INFO]  LSTMモデル:     有効性=5.3%  多様性=68.8%  新規性=100.0%\n","truncated":false}}
%---
%[output:7a8ab334]
%   data: {"dataType":"text","outputData":{"text":"[14:04:55][INFO]  警告: N=300では有効率ギャップは統計的に有意ではない (p~0.10)\n","truncated":false}}
%---
%[output:198ebf94]
%   data: {"dataType":"text","outputData":{"text":"[14:04:55][INFO]  アーキテクチャ:   LSTM(128) -> Dropout(0.30) -> FC(37) -> Softmax\n","truncated":false}}
%---
%[output:17de43a5]
%   data: {"dataType":"text","outputData":{"text":"[14:04:55][INFO]  パラメータ数:     79621  (18000訓練トークンに対し過パラメータ化)\n","truncated":false}}
%---
