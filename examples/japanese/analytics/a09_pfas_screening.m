%[text] # A09: PFAS と環境化学物質スクリーニング
%[text]
%[text] 製造施設近くの地下水から未知の化学物質が検出された場合、どの化学物質が「永遠の化学物質」であるPFASなのかをどのように判断するかを考えてみましょう。
%[text] PFASは、ノンスティックコーティングや消火泡剤に1940年代から使用されてきた9,000種類以上の合成化学物質群です。C-F結合（約544 kJ/mol）のため、生物学的・環境的に分解されにくい特性を持っています。
%[text] SMARTSパターンと物理化学的スコアリングを組み合わせることで、多数の候補を自動的に絞り込み、懸念度に応じて優先順位を付けることができます。
%[text] このスクリプトでは、水質管理機関の環境毒性学者の視点から、PFASのスクリーニング、持続性スコアリング、重み最適化、Tanimoto類似度分析を体験します。
%[text]
%[text] **学習目標**
%[text] - PFASの構造的定義を理解する（OECD 2021基準）
%[text] - `emk.mol.hasSubstruct()`を使用してマルチパターンSMARTSスクリーニングを適用する
%[text] - 専門知識を定量的懸念スコアとしてエンコードする方法を学ぶ
%[text] - fminconを用いて制約付き重み最適化問題を解く
%[text] - 構造クラスタリングのためのTanimoto類似度ヒートマップを解釈する
%[text] - スクリーニング結果を優先順位付き懸念テーブルとして報告する
%[text]
%[text] **ワークフロー（4段階）:**
%[text] 1. **PFASフラグ付け** — SMARTSを用いてパーフルオロアルキル鎖、CF3末端、スルホニル頭部基を検出する
%[text] 2. **持続性スコアリング** — LogP、フッ素数、TPSAの3つのプロキシを用いて環境持続性を推定する
%[text] 3. **多基準優先順位付け** — fminconを用いて参照ルーブリックに最も適した重みセットを求める
%[text] 4. **レポート** — ランク付き懸念テーブル、散布図、Tanimoto類似度ヒートマップを作成する
%[text]
%[text] **前提条件**
%[text] - A07（スキャフォールド分析）修了 — SMARTSとサブ構造の概念を理解していること
%[text] - Optimization Toolbox（fmincon）— セクション4に必要（なければ等重みフォールバックを自動使用）
%[text] - Statistics and Machine Learning Toolbox — セクション5に必要（`clusterdata`、`linkage`、`dendrogram`）
%[text] - 両ツールボックスはMATLAB Online Basic（無料枠）に含まれています
%[text]
%[text] データ: 全分子はインライン（SMILESリテラル）で定義されており、外部ファイルは不要です（US EPA CompTox + OECD 2021基準）。
%[text]
%[text] 推定所要時間: 35〜50分 | 実行方法: Ctrl+Enterでセクションを1つずつ実行
%%
%[text] ## セクション 0: セットアップ

% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
emk.setup.initPython();

%[text] メイン実行に先立ち、Python/RDKit プロセスをウォームアップします
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;

%[text] 使用可能なツールボックスを確認します
hasOptTbx   = license("test", "optimization_toolbox");
hasStatsTbx = license("test", "statistics_toolbox");

if hasOptTbx
    logInfo("A09: Optimization Toolbox を検出 -- 重み最適化が有効。");
else
    logWarn("A09: Optimization Toolbox が検出されない。");
    logWarn("     セクション 4 は等重みフォールバックを使用。");
end

if hasStatsTbx
    logInfo("A09: Statistics and ML Toolbox を検出 -- 階層クラスタリングが有効。");
else
    logWarn("A09: Statistics and ML Toolbox が検出されない。");
    logWarn("     セクション 5 のクラスタリングはスキップ。");
end

logSection("A09", "セクション 0: セットアップ", "アナリティクス L3");
logInfo("A09: セットアップ完了。");
%%
%[text] ## セクション 1: テスト化学物質インベントリの定義
%[text]
%[text] セットアップが完了しました。まず、分析対象となる20種類の化学物質インベントリを定義します。
%[text] PFASの4つのクラスと非PFASコントロールを含む現実的なセットです。
%[text]
%[text] ### コンセプト: PFASとは何か？
%[text]
%[text] OECD（2021）の定義:
%[text] PFASとは、少なくとも1つのパーフルオロ化メチル（-CF3）またはメチレン（-CF2-）炭素原子を含む物質です。ただし、すべてのC-F結合がヘテロ原子に直接結合している場合を除きます（例: -CF2Cl、SF6は除外）。
%[text]
%[text] 主要な構造クラス:
%[text]
%[text] 1. パーフルオロアルキルカルボン酸（PFCAs）
%[text] 一般式: F(CF2)_n-COOH
%[text] 例: PFOA（n=7）、PFBA（n=3）、PFHxA（n=5）
%[text] SMARTSフラグ: パーフルオロアルキル鎖 [C](F)(F)(F) または [CF2] 繰り返し
%[text]
%[text] 2. パーフルオロアルキルスルホン酸（PFSAs）
%[text] 一般式: F(CF2)_n-SO3H
%[text] 例: PFOS（n=8）、PFBS（n=4）、PFHxS（n=6）
%[text] SMARTSフラグ: パーフルオロアルキル鎖と組み合わせた -S(=O)(=O)O
%[text]
%[text] 3. フルオロテロマー物質（FTS / FTOH）
%[text] 構造: F(CF2)_n-CH2CH2-（ポリフルオロアルキル: 一部のC-H結合あり）
%[text] 環境中でPFCAに分解する前駆体化学物質です。本セットでは頭部基が-OHのフルオロテロマーアルコール（FTOH）を使用し、クラスラベルは"FTS"で統一しています。
%[text] SMARTSフラグ: -CF2-CH2- 接合
%[text]
%[text] 4. 非PFASフッ素化化合物
%[text] 薬物（フルオキセチン、シプロフロキサシン）と農薬（フルビプロフェン）を含みます。
%[text] 孤立したC-F結合を含むが、パーフルオロアルキル鎖は持ちません。
%[text] SMARTSフラグ: パーフルオロアルキル鎖テストに失敗 -> 非PFASに分類
%[text]
%[text] 20種類の化学物質のテストセットは、4つのクラスすべてとUS EPA CompTox PFASユニバースリストからのいくつかの非PFASコントロールを網羅しています。
%[text]
%[text] 各エントリ: {表示名, SMILES, 真のクラスラベル}
%[text] クラスラベル: "PFCA"、"PFSA"、"FTS"、"NonPFAS"
logSection("A09", "セクション 1: テスト化学物質インベントリの定義", "アナリティクス L3");
CHEMICALS = { ...
    % --- PFCAs (パーフルオロアルキルカルボン酸) ---
    "PFBA",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)F",                           "PFCA";    ...  % C4 PFCA: COOH-CF2-CF2-CF3（4 炭素、完全フッ素化）
    "PFHxA",   "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",           "PFCA";    ...  % C6 PFCA
    "PFOA",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "PFCA"; ...  % C8 PFCA（OECD 附属書 B 物質）
    "PFNA",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F","PFCA"; ... % C9 PFCA
    "PFDA",    "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F","PFCA"; ... % C10 PFCA
    % --- PFSAs (パーフルオロアルキルスルホン酸) ---
    "PFBS",    "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",               "PFSA";    ...  % C4 PFSA: SO3H-CF2-CF2-CF2-CF3（4 炭素、完全フッ素化）
    "PFHxS",   "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F","PFSA";    ...  % C6 PFSA
    "PFOS",    "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F","PFSA"; ... % C8 PFSA（ストックホルム条約 附属書 B）
    % --- フルオロテロマーとポリフルオロアルキル前駆体 ---
    "6:2FTS",  "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",       "FTS";     ...  % 6:2 フルオロテロマースルホナート
    "8:2FTS",  "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F","FTS"; ... % 8:2 フルオロテロマー
    "4:2FTS",  "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                     "FTS";     ...  % 4:2 フルオロテロマー
    % --- 非 PFAS フッ素化医薬品・農薬（コントロール）---
    "Fluoxetine", "CNCCC(c1ccccc1)Oc1ccc(cc1)C(F)(F)F",                 "NonPFAS"; ...  % 抗うつ薬（1 トリフルオロメチル）
    "Ciprofloxacin","OC(=O)c1cn(C2CCNCC2)c2cc(F)c(nc12)N1CCCC1",        "NonPFAS"; ...  % 抗生物質（1 C-F）
    "Flurbiprofen","OC(=O)C(C)c1ccc(cc1)-c1cccc(F)c1",                  "NonPFAS"; ...  % NSAID（1 C-F）
    "Diflunisal","OC(=O)c1ccc(cc1O)-c1ccc(F)cc1F",                      "NonPFAS"; ...  % NSAID（孤立 C-F 2 個）
    "Trifluridine","OC1C(CO)OC(n2cc(c(=O)[nH]2)C(F)(F)F)C1F",           "NonPFAS"; ...  % 抗ウイルス薬（混合ハロゲン化）
    "Halothane","FC(F)(F)C(Cl)Br",                                       "NonPFAS"; ...  % 麻酔薬（3F だが 1C 上、鎖なし）
    "Sevoflurane","FC(F)(F)C(F)OCC(F)(F)F",                              "NonPFAS"; ...  % 揮発性麻酔薬
    "Flecainide","OCC(NC(=O)c1cc(OCC(F)(F)F)c(cc1OCC(F)(F)F)OCC(F)(F)F)CC", "NonPFAS"; ... % 抗不整脈薬（OCH2CF3 基）
    "Ethanol",  "CCO",                                                   "NonPFAS"; ...  % 陰性コントロール（F なし）
};

nChem   = size(CHEMICALS, 1);
names   = string(CHEMICALS(:,1));
smiles  = string(CHEMICALS(:,2));
trueCls = string(CHEMICALS(:,3));

logInfo("インベントリ: %d 化学物質（PFCA %d、PFSA %d、FTS %d、NonPFAS %d）", ...
    nChem, sum(trueCls=="PFCA"), sum(trueCls=="PFSA"), ...
    sum(trueCls=="FTS"),  sum(trueCls=="NonPFAS"));

%[text] 分子をパースします。
mols  = cell(1, nChem);
valid = false(1, nChem);
for k = 1:nChem
    try
        mols{k} = emk.mol.fromSmiles(smiles(k));
        valid(k) = true;
    catch ME
        logWarn("%s をパースできない: %s", names(k), ME.message);
    end
end
logInfo("%d / %d 分子をパースしました。", sum(valid), nChem);

%[text] **💡 観察ポイント 1 — 炭素鎖長とクラス分けを確認しましょう**
%[text] PFOSは8個のフッ素化炭素を持つのに対し、PFHxSは6個のみです。名前の数字プレフィックス（PFHxS = 6、PFOS = 8）は炭素鎖長を示しています。
%[text] 上記の各SMILESで-CF2-ユニットを数えて、炭素数が名前と一致するか確認しましょう。
%[text] フルオキセチンは-CF3基を持ちますが、NonPFASに分類されます。なぜでしょうか？
%[text] （ヒント: -CF3単独ではパーフルオロアルキル鎖を形成しません。連続した完全フッ素化炭素が必要です。）
%[text] ハロタン（`FC(F)(F)C(Cl)Br`）は1つの炭素上に3個のフッ素を持ちます。OECD定義ではPFASでしょうか？
%[text] （隣接する非Fハロゲンを持つ単一の-CF3は「少なくとも1つの-CF3または-CF2-」に適合しますか？）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: SMARTS ベースの PFAS フラグ付け
%[text]
%[text] インベントリが定義されました。次に、SMARTS パターンを用いて各分子が PFAS（パーフルオロアルキル化合物）かどうかを自動判定します。
%[text] 各パターンがどの化合物にマッチするかを確認しましょう。
%[text]
%[text] ### コンセプト: PFAS 検出のための SMARTS パターン
%[text]
%[text] SMARTS（SMiles ARbitrary Target Specification）パターンは、原子と結合の制約を持つ構造クエリをエンコードします。特定の分子を記述する SMILES とは異なり、SMARTS は多くの分子に共通する構造モチーフを記述します。
%[text]
%[text] 以下で使用する主要な SMARTS 原子:
%[text] [#6]    -- 任意の炭素原子（原子番号 6）
%[text] [F]     -- フッ素（元素記号、電荷制約なし）
%[text] [#16]   -- 任意の硫黄原子（原子番号 16）
%[text] [!F]    -- フッ素ではない任意の原子
%[text] [$(...)] -- 再帰 SMARTS（ネストされたサブクエリ）
%[text]
%[text] 3 つの補完的なパターンが OECD PFAS 構造定義をカバーします:
%[text]
%[text] (A) パーフルオロアルキル鎖  [#6](F)(F)-[#6](F)(F)
%[text] 正確に 2 つの F に結合した炭素が別の同じ炭素に隣接する場合にマッチします。
%[text] 完全フッ素化炭素鎖の -CF2-CF2- シグネチャです。
%[text] 少なくとも 2 つの連続する CF2 ユニットが必要です。
%[text]
%[text] (B) トリフルオロメチル末端  [#6](F)(F)F
%[text] -CF3 にマッチします。パーフルオロアルキル鎖の末端です。
%[text] 孤立した CF3 基（コントロール）にもマッチするため、(B) 単独では広すぎます。信頼できる PFAS フラグ付けには (A) との組み合わせが必要です。
%[text]
%[text] (C) スルホニル頭部基  [#16](=O)(=O)
%[text] -S(=O)(=O)- コアにマッチします。(A) と組み合わせて PFSA を特定します。
%[text]
%[text] フラグ付けロジック:
%[text] hasChain   = (A)             --> 完全パーフルオロアルキル鎖が存在
%[text] hasCF3     = (B)             --> CF3 基が存在（孤立している可能性あり）
%[text] hasSulfonyl = (C)            --> スルホン酸/スルホンアミド頭部基が存在
%[text]
%[text] isPFAS     = hasChain OR (hasCF3 AND NOT 非環式)
%[text] classGuess: isPFAS AND hasSulfonyl -> "PFSA"
%[text] isPFAS AND NOT hasSulfonyl -> "PFCA_or_FTS"
%[text] else -> "NonPFAS"
%[text]
%[text] 注意: このルールベース分類器は説明用であり、規制グレードではありません。
%[text] 実際の PFAS スクリーニングワークフロー（例: EPA DSSTox）は、PFAS サブクラスをカバーする数百のパターンを持つ大規模 SMARTS ライブラリを使用します。
%[text]
%[text] SMARTS パターン
logSection("A09", "セクション 2: SMARTS ベースの PFAS フラグ付け", "アナリティクス L3");
SMARTS_PFCHAIN  = "[#6](F)(F)-[#6](F)(F)";    % CF2-CF2 バックボーン
SMARTS_CF3      = "[#6](F)(F)F";              % -CF3 末端
SMARTS_SULFONYL = "[#16](=O)(=O)";            % -SO2- 頭部基
SMARTS_POLYFLUORO = "[#6](F)(F)[#6H2]";          % CH2 に直接結合した CF2（FTS 接合のみ）

validMols = mols(valid);
validNames = names(valid);
validSmiles = smiles(valid);
validCls = trueCls(valid);
nValid = sum(valid);

%[text] サブ構造フラグを適用
hasChain    = emk.mol.hasSubstruct(validMols, SMARTS_PFCHAIN);
hasCF3      = emk.mol.hasSubstruct(validMols, SMARTS_CF3);
hasSulfonyl = emk.mol.hasSubstruct(validMols, SMARTS_SULFONYL);
hasPolyFluoro = emk.mol.hasSubstruct(validMols, SMARTS_POLYFLUORO);

%[text] 分類: パーフルオロアルキル鎖があれば isPFAS
%[text] 注: `hasCF3 & hasChain` は `hasChain` の部分集合であるため、実質的に `isPFAS = hasChain` と等価です（冗長項は明示性のために残しています）。
isPFAS = hasChain | (hasCF3 & hasChain);

%[text] PFAS ヒットをサブ分類
classGuess = repmat("NonPFAS", nValid, 1);
for k = 1:nValid
    if isPFAS(k)
        if hasSulfonyl(k)
            classGuess(k) = "PFSA";
        elseif hasPolyFluoro(k) && ~hasSulfonyl(k)
            classGuess(k) = "FTS";
        else
            classGuess(k) = "PFCA";
        end
    end
end

%[text] 真のラベルに対する混同行列
nPFAS_hit = sum(isPFAS);
nPFAS_true = sum(validCls ~= "NonPFAS");
correct = sum(classGuess == validCls);

logInfo("PFAS フラグ付け結果:");
logInfo("  真の PFAS: %d    検出 PFAS: %d    正しい分類: %d / %d", ...
    nPFAS_true, nPFAS_hit, correct, nValid);

%[text] 結果テーブルを表示
flagTbl = table( ...
    validNames, validCls, classGuess, hasChain', hasCF3', hasSulfonyl', isPFAS', ...
    VariableNames=["Name","TrueClass","GuessedClass","HasChain","HasCF3","HasSulfonyl","IsPFAS"]);
disp(flagTbl);

%[text] **💡 観察ポイント 2 — SMARTS パターンの限界と改良**
%[text] ハロタンやセボフルランなどのコントロールも -CF3 基を持ちます。SMARTS_CF3 パターンはそれらにマッチするでしょうか？
%[text] `hasCF3` の該当行を確認し、なぜ `hasCF3` 単独では PFAS 特定に不十分かを考えてみましょう。
%[text] SMARTS が見逃す、FTS を PFCA から区別する構造的特徴は何でしょうか？（ヒント: FTS はフッ素化鎖に少なくとも 1 つの C-H 結合を持ちます。）
%[text] FTS のより良い SMARTS を設計しましょう。FTS は -CF2-CH2- 接合を持ちます。
%[text] `emk.mol.hasSubstruct(validMols, "[#6](F)(F)-[#6H2]")` を試して、FTS 化合物を正しく分離できるか確認してみましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: 記述子から持続性プロキシを計算
%[text]
%[text] PFAS のフラグ付けが完了しました。次に、3 つの物理化学的プロキシ（LogP: 水/オクタノール間の分配係数の常用対数、フッ素数、TPSA: 位相的極性表面積）を用いて環境持続性を数値化します。
%[text] 各化合物がどのプロキシで高いスコアを示すかを確認しましょう。
%[text]
%[text] ### コンセプト: 環境持続性のための 3 つのプロキシ
%[text]
%[text] 環境半減期の直接測定には高価な実験室アッセイが必要なため、規制当局は物理化学的プロキシを早期警告指標として使用します。
%[text] 広く使用される 3 つのプロキシ（Wang et al. 2017、Cousins et al. 2020）:
%[text]
%[text] (1) 疎水性（LogP）
%[text] 高い LogP は脂質への強い親和性を示し、脂肪組織への生物蓄積を引き起こします。
%[text] 生物蓄積係数（BAF）は中性有機物の LogP と相関することが多いです（log BAF ~ 0.79 * LogP - 0.40、Veith 1979）。
%[text] PFAS では、イオン性頭部基（PFCA は中性 pH で酸）のため LogP の解釈が複雑です。RDKit の Wildman-Crippen LogP は真の分配係数ではなく、構造的プロキシとして使用されます。
%[text]
%[text] 懸念寄与: score_logP = clamp(LogP / 10, 0, 1)
%[text] （LogP = 10 がスコア 1.0 にマップされ、負の LogP はスコア 0 になります）
%[text]
%[text] (2) フッ素数（F 密度）
%[text] C-F 結合が多いと化学的安定性が高まり、分解が遅くなります。
%[text] C-F 結合解離エネルギー（~544 kJ/mol）は C-Cl（~397 kJ/mol）や C-C（~346 kJ/mol）を大きく超えます。
%[text] 単純なプロキシとして、フッ素原子数（F 数）を 0〜1 の範囲に正規化します。
%[text]
%[text] 懸念寄与: score_F = clamp(nF / 17, 0, 1)
%[text] （17 フッ素は PFDA で、セット内で最長の化合物です）
%[text]
%[text] (3) 水溶性代替指標（TPSA）
%[text] 非常に親水性の高い PFAS（高 TPSA）は水により溶出しやすいです。
%[text] 非常に疎水性のもの（低 TPSA）は堆積物に吸着しますが、生物に蓄積します。
%[text] 極端な値で高い懸念があります。ベル形の懸念最大値は TPSA ~ 40 Å² 付近で、中性形の長鎖 PFAS に典型的です。
%[text]
%[text] 懸念寄与: score_TPSA = exp(-((TPSA - 40)^2) / (2*30^2))
%[text] （σ=30 で TPSA=40 を中心とするガウス -- ヒューリスティック、説明用）
logSection("A09", "セクション 3: 記述子から持続性プロキシを計算", "アナリティクス L3");
DESCS = ["LogP", "TPSA", "HeavyAtomCount"];
descTbl = emk.descriptor.batchCalculate(validMols, DESCS);

logp_vec    = descTbl.LogP;            % N x 1
tpsa_vec    = descTbl.TPSA;            % N x 1
hatom_vec   = descTbl.HeavyAtomCount;  % N x 1

%[text] SMILES からフッ素原子を数えます（高速: Python 呼び出し不要）。
nF_vec = zeros(nValid, 1);
for k = 1:nValid
    nF_vec(k) = count(validSmiles(k), "F");
end

%[text] 3 つの持続性プロキシスコアを計算します（各 0〜1）。
score_logP = min(max(logp_vec / 10, 0), 1);
score_F    = min(max(nF_vec / 17,   0), 1);
score_TPSA = exp(-((tpsa_vec - 40).^2) / (2 * 30^2));

scoresMat = [score_logP, score_F, score_TPSA];  % N x 3

logInfo("持続性プロキシ統計（PFAS ヒットのみ）:");
pfasIdx = find(isPFAS);
if ~isempty(pfasIdx)
    logInfo("  LogP スコア    : 平均 = %.3f  最大 = %.3f", ...
        mean(score_logP(pfasIdx)), max(score_logP(pfasIdx)));
    logInfo("  F 数スコア: 平均 = %.3f  最大 = %.3f", ...
        mean(score_F(pfasIdx)),    max(score_F(pfasIdx)));
    logInfo("  TPSA スコア  : 平均 = %.3f  最大 = %.3f", ...
        mean(score_TPSA(pfasIdx)), max(score_TPSA(pfasIdx)));
end

%[text] 全化学物質の 3 つのプロキシスコアを可視化します。
figure("Name","A09 Sec3: 持続性プロキシスコア");
X = 1:nValid;
bar(X, scoresMat, "grouped");
xticks(X);
xticklabels(validNames);
xtickangle(45);
ylabel("スコア（0〜1）");
title("化学物質ごとの 3 つの持続性プロキシスコア");
legend(["LogP プロキシ", "F 数プロキシ", "TPSA プロキシ"], Location="northeast");
grid on;

%[text] **💡 観察ポイント 3 — 持続性プロキシを探索しましょう**
%[text] F 数スコアが最も高い化合物を確認しましょう。`score_F` と `validNames` を調べてみてください。
%[text] 最長の PFAS 鎖が最も多くのフッ素を持つはずですが、結果は期待と一致していますか？
%[text] エタノールがなぜ 3 つのプロキシ全てで 0 に近いスコアになるのかを考えてみましょう。フッ素原子はいくつありますか？
%[text] TPSA プロキシは 40 Å² を中心とするガウスです。PFOS はスルホナート頭部基のため TPSA ≈ 115 Å² と高い値を持ちます。
%[text] ガウスは低極性 PFCA と比べて PFOS を不当に低評価していませんか？このプロキシをどのように再設計しますか？
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: fmincon による重み最適化
%[text]
%[text] 持続性プロキシスコアが揃いました。次に fmincon を使って3つのプロキシに対する最適な重みを求めます。
%[text] 参照ルーブリックへの最小二乗フィットを用いて、科学的根拠に基づいた重み配分を導出します。
%[text]
%[text] ### コンセプト: 多基準スコアリングの制約付き重み最適化
%[text]
%[text] 各化学物質の最終「懸念スコア」は3つの持続性プロキシの加重和です:
%[text]
%[text] score(k) = w1 * score_logP(k) + w2 * score_F(k) + w3 * score_TPSA(k)
%[text]
%[text] 重み（w1、w2、w3）は規制ガイダンスに従って各プロキシの相対的重要度を表します。
%[text] 専門家は、C-F 結合密度が LogP よりも持続性の直接的な尺度であるため、w2 > w1 を割り当てることがあります。
%[text]
%[text] 「最適」な重みを求めるために参照スコアリングルーブリックを定義します:
%[text] - 長鎖 PFAS（PFOA、PFOS、PFNA、PFDA）-> 参照スコア >= 0.7
%[text] - 短鎖 PFAS（PFBA、PFBS、PFHxA）-> 参照スコア < 0.5
%[text] - 非 PFAS コントロール -> 参照スコア <= 0.2
%[text]
%[text] 以下の制約のもとで参照ターゲットからの最大二乗偏差を最小化します:
%[text] w1 + w2 + w3 = 1   （重みの和が1）
%[text] w1, w2, w3 >= 0.05  （各プロキシが少なくとも5%寄与）
%[text] w1, w2, w3 <= 0.80  （単一プロキシが80%以上を支配しない）
%[text]
%[text] これは fmincon で解く箱制約付き二次問題です。
%[text] 目的関数: sum_k (score(k) - target(k))^2（最小二乗フィット）
%[text]
%[text] 化学物質ごとの参照ターゲットスコア
%[text] 長鎖 PFAS: 高懸念（ターゲット >= 0.7）
%[text] 短鎖 PFAS: 中懸念（ターゲット < 0.5）
%[text] 非 PFAS: 低懸念（ターゲット <= 0.2）
targetScores = zeros(nValid, 1);
for k = 1:nValid
    cls = validCls(k);
    nm  = validNames(k);
    if cls == "NonPFAS"
        targetScores(k) = 0.10;
    elseif any(nm == ["PFOA","PFOS","PFNA","PFDA"])
        targetScores(k) = 0.80;   % 長鎖、高懸念
    elseif any(nm == ["PFHxA","PFHxS","8:2FTS"])
        targetScores(k) = 0.60;   % 中鎖、中〜高懸念
    elseif any(nm == ["PFBA","PFBS","6:2FTS","4:2FTS"])
        targetScores(k) = 0.40;   % 短鎖、中懸念
    else
        targetScores(k) = 0.50;   % 未割り当て PFAS のデフォルト
    end
end

%[text] 目的関数: ターゲットスコアへの最小二乗フィット
%[text] f(w) = sum_k (scoresMat * w - targetScores)^2
logSection("A09", "セクション 4: fmincon による重み最適化", "アナリティクス L3");
objFun = @(w) sum((scoresMat * w - targetScores).^2);

%[text] 制約
%[text] Aeq * w = beq  -->  w1 + w2 + w3 = 1
Aeq  = [1, 1, 1];
beq  = 1;
lb   = [0.05; 0.05; 0.05];   % 下限
ub   = [0.80; 0.80; 0.80];   % 上限
w0   = [1/3; 1/3; 1/3];       % 初期推定: 等重み

if hasOptTbx
    opts = optimoptions("fmincon", Display="off", Algorithm="sqp");
    [w_opt, fval] = fmincon(objFun, w0, [], [], Aeq, beq, lb, ub, [], opts);
    logInfo("fmincon 最適化重み:  w_logP=%.3f  w_F=%.3f  w_TPSA=%.3f  (RSS=%.4f)", ...
        w_opt(1), w_opt(2), w_opt(3), fval);
else
    w_opt = w0;   % equal-weight fallback
    logWarn("等重み [1/3、1/3、1/3] を使用（Optimization Toolbox が利用できない）。");
end

%[text] 最適化した重みで最終懸念スコアを計算します。
concernScore = scoresMat * w_opt;

%[text] トップ懸念を表示します。
[sortedScores, sortIdx] = sort(concernScore, "descend");
logInfo("ランク付き懸念テーブル（上位 10）:");
for k = 1:min(10, nValid)
    i = sortIdx(k);
    logInfo("  %2d. %-18s  [%s]  スコア = %.3f", k, validNames(i), validCls(i), sortedScores(k));
end

%[text] 棒グラフ: ランク付き懸念スコア
figure("Name","A09 Sec4: ランク付き化学物質懸念スコア");
barColors = zeros(nValid, 3);
for k = 1:nValid
    i = sortIdx(k);
    switch validCls(i)
        case "PFCA",    barColors(k,:) = [0.85, 0.20, 0.10];   % red
        case "PFSA",    barColors(k,:) = [0.95, 0.55, 0.10];   % orange
        case "FTS",     barColors(k,:) = [0.25, 0.55, 0.85];   % blue
        case "NonPFAS", barColors(k,:) = [0.70, 0.70, 0.70];   % grey
    end
end
barH = bar(sortedScores, "FaceColor","flat", "HandleVisibility","off");
barH.CData = barColors;
xticks(1:nValid);
xticklabels(validNames(sortIdx));
xtickangle(55);
ylabel("懸念スコア（0〜1）");
title(sprintf("環境懸念スコア (w=[%.2f, %.2f, %.2f])", ...
    w_opt(1), w_opt(2), w_opt(3)));
%[text] 手動凡例パッチ
patch(NaN,NaN,[0.85 0.20 0.10], DisplayName="PFCA");
patch(NaN,NaN,[0.95 0.55 0.10], DisplayName="PFSA");
patch(NaN,NaN,[0.25 0.55 0.85], DisplayName="FTS");
patch(NaN,NaN,[0.70 0.70 0.70], DisplayName="NonPFAS");
legend(Location="northeast");
grid on;
yline(0.5, "--k", "懸念閾値", LabelHorizontalAlignment="left", HandleVisibility="off");

%[text] **💡 観察ポイント 4 — 重み最適化の結果を確認しましょう**
%[text] w_logP、w_F、w_TPSA のうち、どのプロキシが最も高い重みを受け取りましたか？
%[text] 科学文献と一致していますか？（Buck et al. 2011: C-F 結合数が主要な要因と示唆）
%[text] PFOS のターゲットスコアを 0.80 から 0.95 に変えて再実行した場合、`w_opt` はどう変わりますか？
%[text] F 数の重みは増加しますか？（PFOS は 17 個の F 原子を持ちます。）
%[text] 4 番目のプロキシ（-CF2- ユニット数 = 鎖長）を追加したらどうなるでしょうか？
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: PFAS ヒットの構造類似度ヒートマップ
%[text]
%[text] 懸念スコアの計算が完了しました。次に、PFAS ヒット間の構造類似度をヒートマップで視覚化します。
%[text] 同じクラスの PFAS はどれほど似ているかを確認しましょう。
%[text]
%[text] ### コンセプト: Tanimoto 類似度と構造クラスタリング
%[text]
%[text] 2 分子間の Tanimoto 係数は Morgan（円形）フィンガープリントから計算されます。
%[text]
%[text] T(A, B) = |A AND B| / |A OR B|
%[text]
%[text] ここで |A AND B| は両フィンガープリントでセットされているビット数、|A OR B| はいずれかでセットされているビット数を表します。
%[text]
%[text] T = 1.0 は同一フィンガープリント（構造的に区別不可）を示します。
%[text] T = 0.0 は共有ビットがない（完全に非類似）ことを示します。
%[text]
%[text] PFAS では、長鎖ホモログ（PFOA と PFNA）は -CF2- 繰り返し単位の数のみが異なります。
%[text] フィンガープリントは高度に類似しており（T ~ 0.7〜0.9）、追加の -CF2- 繰り返し以外の全サブ構造を共有しています。
%[text]
%[text] 対照的に、PFCA と PFSA は頭部基が異なります: -COOH と -SO3H。
%[text] Tanimoto 類似度は中程度（T ~ 0.3〜0.6）です。
%[text]
%[text] 階層クラスタリングデンドログラム（平均連結 / UPGMA）は最も類似した構造をグループ化し、PFCA / PFSA / FTS のサブファミリーを明らかにします。
%[text] Tanimoto 距離は非ユークリッドであるため、Ward の代わりに平均連結を使用します。Ward 連結はユークリッド空間を仮定します（Willett 1998）。
logSection("A09", "セクション 5: PFAS ヒットの構造類似度ヒートマップ", "アナリティクス L3");
pfasHitIdx = find(isPFAS);
nHits      = numel(pfasHitIdx);

if nHits > 1
    % PFAS ヒットの Morgan フィンガープリントを計算
    fps = cell(1, nHits);
    for k = 1:nHits
        fps{k} = emk.fingerprint.morgan(validMols{pfasHitIdx(k)});
    end

    % Tanimoto 類似度行列を構築
    simMat = zeros(nHits, nHits);
    for i = 1:nHits
        for j = 1:nHits
            simMat(i,j) = emk.similarity.tanimoto(fps{i}, fps{j});
        end
    end

    hitNames = validNames(pfasHitIdx);

    % ヒートマップとして可視化
    figure("Name","A09 Sec5: PFAS Tanimoto 類似度ヒートマップ");
    imagesc(simMat);
    colormap(hot);
    colorbar();
    clim([0 1]);
    xticks(1:nHits); xticklabels(hitNames); xtickangle(45);
    yticks(1:nHits); yticklabels(hitNames);
    title("Morgan フィンガープリント Tanimoto 類似度（PFAS ヒット）");
    axis square;

    % 階層クラスタリング（Statistics and ML Toolbox）
    % Tanimoto 距離は非ユークリッドのため平均連結（UPGMA）を使用;
    % Ward 連結はユークリッド距離を必要とし偽の警告を生成する。
    % 平均連結はケモインフォマティクスフィンガープリントクラスタリングの
    % 標準的な選択（Willett 1998）。
    if hasStatsTbx
        distVec  = squareform(1 - simMat);   % condensed distance vector
        Z        = linkage(distVec, "average");
        figure("Name","A09 Sec5: PFAS 平均連結デンドログラム");
        dendrogram(Z, "Labels", hitNames, Orientation="left");
        title("平均連結デンドログラム -- PFAS 構造ファミリー");
        xlabel("距離（1 - Tanimoto）");
    else
        logWarn("Statistics and ML Toolbox が利用できない -- デンドログラムをスキップ。");
    end

    logInfo("PFAS ヒットの Tanimoto 統計:");
    upperTri = simMat(triu(true(nHits), 1));
    logInfo("  平均 = %.3f   最小 = %.3f   最大 = %.3f", ...
        mean(upperTri), min(upperTri), max(upperTri));
else
    logWarn("PFAS ヒットが 2 未満; ヒートマップをスキップ。");
end

%[text] **💡 観察ポイント 5 — Tanimoto 類似度を深堀りしましょう**
%[text] 最も高い Tanimoto 類似度を持つ PFAS 化合物のペアを確認しましょう。
%[text] PFOA（PFCA、C8）と PFOS（PFSA、C8）の Tanimoto 類似度を確認しましょう。同じ炭素鎖長ですが、頭部基が異なります。
%[text] T は 0.5（中程度）に近いか、0.9（高）に近いかを確認しましょう。
%[text] デンドログラムで PFOS が PFBS・PFHxS（他の PFSA）と離れ、FTS 側に配置される理由を考えましょう。
%[text] （ヒント: PFOS と 8:2FTS はどちらも F 原子数 = 17。Morgan FP は頭部基の差より長鎖の類似性を強調することがあります。）
%[text] Morgan フィンガープリントの半径を `Radius=1` に変えて再計算してみましょう。
%[text] `fps{k} = emk.fingerprint.morgan(validMols{pfasHitIdx(k)}, Radius=1);`
%[text] 半径を小さくすると PFCA と PFSA ファミリー間の類似度は増加するか、減少するかを考えましょう。なぜでしょうか？
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: サマリーレポートとスクリーニングダッシュボード
%[text]
%[text] 構造クラスタリングが完了しました。最後に、すべての分析結果を一つのスクリーニングレポートにまとめます。
%[text] 意思決定者向けにランク付けされた懸念テーブルと散布図を作成します。
%[text]
%[text] ### コンセプト: スクリーニング結果の伝達
logSection("A09", "セクション 6: サマリーレポートとスクリーニングダッシュボード", "アナリティクス L3");
%[text] すべての化学物質のサマリーテーブルを構築します。
[~, rankOrder] = sort(concernScore, "descend");
rankVec        = zeros(nValid, 1);
rankVec(rankOrder) = (1:nValid)';

summaryTbl = table( ...
    rankVec, ...
    validNames, ...
    validCls, ...
    classGuess, ...
    round(logp_vec,2), ...
    nF_vec, ...
    round(tpsa_vec,1), ...
    round(concernScore,3), ...
    VariableNames=["Rank","Name","TrueClass","FlaggedClass", ...
                   "LogP","nF","TPSA_A2","ConcernScore"]);
summaryTbl = sortrows(summaryTbl, "Rank");

logInfo("=== PFAS スクリーニングレポート ===");
disp(summaryTbl);

%[text] 散布図: 懸念スコアで色付けした LogP 対 F 数
figure("Name","A09 Sec6: LogP vs F 数（懸念スコア）");
scatter(logp_vec, nF_vec, 80, concernScore, "filled", MarkerFaceAlpha=0.8);
colormap(turbo);
cb = colorbar();
cb.Label.String = "懸念スコア";
clim([0 1]);
for k = 1:nValid
    text(logp_vec(k) + 0.1, nF_vec(k), validNames(k), FontSize=7);
end
xlabel("LogP（Wildman-Crippen）");
ylabel("フッ素原子数");
title("化学空間ビュー: LogP vs F 数（色 = 懸念スコア）");
grid on;

%[text] **まとめ**
%[text] SMARTS パターンによる PFAS 自動検出、物理化学的プロキシを用いた持続性スコアリング、fmincon による制約付き重み最適化、Tanimoto 類似度ヒートマップと階層クラスタリングを学びました。
%[text] 規制スクリーニングの要件（ランク付け、正当化、反証可能）を満たすダッシュボードを構築しました。
logInfo("A09: 完了。");

%[text] **💡 観察ポイント 6 — 最終スクリーニング結果を解釈しましょう**
%[text] 散布図の右上領域（高 LogP かつ高 F 数）に高懸念 PFAS の視覚的なクラスターがあるか確認しましょう。
%[text] PFAS 化合物に対して偽陰性が最も少ない単一のスクリーニング基準は、LogP 閾値と F 数閾値のどちらかを確認しましょう。
%[text] セボフルランとハロタンは C-F 結合を持ちますが、揮発性麻酔薬であり、環境残留性は PFAS とは異なります。
%[text] 麻酔ガスを真の PFAS とは別にフラグ付けするには、スコアリングルーブリックをどのように変更するか考えてみましょう。
%[text] （揮発性プロキシの追加を検討: MW < 200 g/mol を最終スコアを下げる「揮発性免除」として。）
%[text] サマリーテーブルを CSV にエクスポートしましょう: `writetable(summaryTbl, "pfas_screening_report.csv")`
%[text] Excel で開いてみましょう。非専門家のステークホルダーにどのように提示するか、どの列を維持または削除するか考えてみましょう。
% ... （ここにコードを書いてみましょう）

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
