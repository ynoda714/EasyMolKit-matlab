%[text] # F05: 部分構造検索と SMARTS パターン
%[text] EasyMolKit 基礎チュートリアル — レイヤー 1
%[text] 
%[text] 「カルボン酸を持っている薬だけをすべて取り出したい」あるいは「ベンゼン環を含む分子だけをスクリーニングしたい」と思ったとき、分子の形をそのまま検索の条件（クエリ）に指定できたら便利だと思いませんか？それを1行の文字列でエレガントに実現してくれるのが「SMARTS（スマート）」という言語です。特定の1つの分子を正確に表す SMILES に対し、SMARTS は「カルボン酸を持つ分子のグループ」というように、分子の『共通する特徴やクラス』を記述することができます。ベンゼン環やアミド結合といった化学者の直感をそのまま検索クエリにできるため、創薬の現場では特定の官能基を持つ化合物の一括スクリーニングや、実験で悪さをする問題児（毒性構造）を排除するための強力なフィルターとして日常的に使われています。このスクリプトでは、SMARTS を使った部分構造検索の基本と、実務で必須となる「PAINS フィルタリング」のワークフローを体験してみましょう。
%[text] **学習目標**
%[text] - SMARTS（分子クエリのための拡張パターン言語）を理解する
%[text] - `emk.mol.hasSubstruct` で分子が部分構造を含むか確認する
%[text] - SMARTS フィルタを分子データセットに適用する
%[text] - よく使う官能基パターンを認識する \
%[text] **前提条件**
%[text] - RDKitインストール済み（`emk.setup.install()` を一度だけ実行しておく）
%[text] - F01～F04
%[text] - 追加Toolbox不要（MATLAB だけで動きます） \
%[text] 所要時間: 10〜15 分 | 実行方法: Ctrl+Enter でセクションを一つずつ実行
%%
%[text] ## セクション 0: セットアップ
%[text] パスと Python 環境を初期化します。
%[text] **常にこのセクションを最初に実行してください。**
% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
emk.setup.initPython(); %[output:444a3b00]
%[text] Python/RDKit プロセスのウォームアップ
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;
logInfo("F05: セットアップ完了"); %[output:186a4c5b]
%%
%[text] ## セクション 1: SMARTS とは何か？
%[text] ### **「1つの分子」を指すか、「共通の特徴」を指すか**
%[text] 一見すると SMILES も SMARTS も同じアルファベットの文字列に見えるため、初学者は混乱しがちですが、その役割は根本的に異なります。例えるなら、SMILES は「東京都新宿区西新宿1-1-1」という特定の家（分子）を指す『住所』です。これに対して SMARTS は、「赤い屋根で、2階建ての、庭がある家」という条件に合う物件をすべて引っかけるための『検索条件』にあたります。
%[text] たとえば、SMILES で `"C"` と書けばそれは「メタン」という特定の分子そのものを表しますが、SMARTS で `[#6]` と書けば、それはメタンだけでなく、アスピリンだろうがDNAだろうが「分子の中のどこかに炭素（原子番号6番）が1つでも含まれていればすべてマッチする」というルールになります。この強力なパターンマッチングの仕組みがあるからこそ、私たちは複雑な構造式の中から「特定の骨格だけ」をパズルのように見つけ出すことができるのです。
%[text] **SMILES との主な違い:**
%[text] - ワイルドカード使用可能: `[#6]` = 任意の炭素、`[*]` = 任意の原子
%[text] - 結合指定可能: `[-]` 単結合、`[=]` 二重結合、`[:]` 芳香族結合
%[text] - 論理演算子: `[#6,#7]` = 炭素 OR 窒素 \
%[text] **よく使う SMARTS パターン:**
%[text] - `c1ccccc1` — 芳香族ベンゼン環
%[text] - `[OH]` — ヒドロキシル基（-OH）
%[text] - `[NH2]` — 第一級アミン
%[text] - `C(=O)[OH]` — カルボン酸
%[text] - `C(=O)[O]` — エステルまたは酸
%[text] - `C(=O)[N]` — アミド
%[text] - `[F,Cl,Br,I]` — 任意のハロゲン \
%%
%[text] ## セクション 2: 基本的な部分構造チェック
%[text] ### **分子の中に「あの骨格」が隠れているか調べる**
%[text] EasyMolKit で部分構造検索を行うコア関数が `emk.mol.hasSubstruct` です。使い方はとても直感的で、調べたい分子（または複数の分子の配列）と、条件となる SMARTS 文字列をこの関数に渡すだけです。条件にマッチすれば `true`（1）、含まれていなければ `false`（0）が返ってきます。まずは、お馴染みのアスピリンの中に「ベンゼン環（SMARTS: `c1ccccc1`）」が含まれているかどうかを、コンピューターに判定させてみましょう。
%[text] `emk.mol.hasSubstruct(mol, query)` を使うと、分子がクエリ部分構造を含むかどうかを `true`/`false` で確認できます。クエリには SMARTS 文字列または SMILES 文字列を指定します。
mol_aspirin     = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");
mol_benzene     = emk.mol.fromSmiles("c1ccccc1");
mol_ethanol     = emk.mol.fromSmiles("CCO");
mol_caffeine    = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
mol_acetaminophen = emk.mol.fromSmiles("CC(=O)NC1=CC=C(C=C1)O");
%[text] アスピリンはベンゼン環を含むか？
has_ring = emk.mol.hasSubstruct(mol_aspirin, "c1ccccc1");
logInfo("アスピリンはベンゼン環を含む: %d", has_ring); %[output:42b5f0c9]
%[text] エタノールはベンゼン環を含むか？
no_ring = emk.mol.hasSubstruct(mol_ethanol, "c1ccccc1");
logInfo("エタノールはベンゼン環を含む: %d", no_ring); %[output:65242d3e]
%[text] アスピリンにカルボン酸はあるか？
has_cooh = emk.mol.hasSubstruct(mol_aspirin, "C(=O)[OH]");
logInfo("アスピリンは -COOH を含む: %d", has_cooh); %[output:8b36612b]
%[text] アセトアミノフェンにアミドはあるか？
has_amide = emk.mol.hasSubstruct(mol_acetaminophen, "C(=O)[NH]");
logInfo("アセトアミノフェンはアミドを含む: %d", has_amide); %[output:50c7289a]
%[text] 
%[text] **✏️ やってみよう 1 — カフェインで SMARTS を試してみましょう**
%[text] カフェイン（`"CN1C=NC2=C1C(=O)N(C(=O)N2C)C"`）が以下の部分構造を含むか確認してみましょう。
%[text] - (a) 任意の窒素原子: `"[#7]"`
%[text] - (b) ケトン基: `"[C](=O)[#6]"`
%[text] - (c) ヒドロキシル基: `"[OH]"` \
%[text] (a) は `true`、(b)(c) は `false` になります。なぜ (b) が `false` なのか考えてみましょう。
%[text] ヒント: カフェインの C=O 基はアミド/イミドカルボニルで、炭素ではなく窒素（N）に結合しています。`"[C](=O)[N]"` で試してみましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: バッチ部分構造スクリーニング
%[text] ### **溜まったデータから目的の分子をまとめてあぶり出す**
%[text] `emk.mol.hasSubstruct` の素晴らしいところは、1つの分子だけでなく、複数の分子を詰め込んだ「セル配列（Cell array）」を丸ごと放り込める（バッチ処理）という点です。関数に分子のリストを渡すと、それぞれの分子が条件を満たしているかどうかの結果が、MATLAB の論理配列（`logical`）として一撃で返ってきます。
%[text] これを使えば、「手元にある10個の化合物の中から、水酸基（ヒドロキシル基: `[OH]`）を持っているものだけを自動で抜き出す」といったフィルター処理が、ループ文（for文）を書くことなくスマートに実行できます。返ってきた論理配列をそのまま MATLAB のインデックス（論理インデックス）として使えば、条件にマッチした分子だけを抽出するコードもすっきりと記述できます。
%[text] `emk.mol.hasSubstruct` に Mol オブジェクトのセル配列を渡すと、データセット全体を 1 回の呼び出しでスキャンできます。結果は部分構造が存在する位置が `true` の論理行ベクトルで返されます。
mols  = {mol_aspirin, mol_benzene, mol_ethanol, mol_caffeine, mol_acetaminophen};
names = {"アスピリン", "ベンゼン", "エタノール", "カフェイン", "アセトアミノフェン"};

% Note: "c1ccccc1" matches benzene (6-membered all-carbon aromatic ring) only.
% Heterocyclic aromatics (pyridine, imidazole...) are NOT detected by this SMARTS.
% Caffeine has imidazole/pyrimidine rings so its benzene ring result is false -- expected.
patterns = struct( ...
    "description", {"ベンゼン環", "カルボン酸", "ヒドロキシル(-OH)", ...
                    "アミド(-CONH)", "窒素原子"}, ...
    "smarts",      {"c1ccccc1", "C(=O)[OH]", "[OH]", "C(=O)[N]", "[#7]"} ...
);

hit_mat   = false(numel(patterns), numel(mols));
for i = 1:numel(patterns)
    hit_mat(i, :) = emk.mol.hasSubstruct(mols, patterns(i).smarts);
end
col_names = ["Aspirin", "Benzene", "Ethanol", "Caffeine", "Acetaminophen"];
row_names = string({patterns.description})';
fg_tbl    = array2table(hit_mat, "VariableNames", col_names, "RowNames", row_names);
logInfo("官能基スクリーニング結果:"); %[output:24f38584]
disp(fg_tbl); %[output:08a7c914]
%%
%[text] ## セクション 4: PAINS フィルタリング
%[text] ### **化学者のこだわりを表現する特殊な記号たち**
%[text] SMARTS 言語の真価は、単に元素記号を並べるだけでなく、原子の「状態」を細かく指定できる高度な表現力にあります。たとえば、`[OH]` は酸素に水素が直接くっついている「水酸基」を厳密に指しますし、`C(=O)[OH]` と書けば「カルボン酸」の綺麗な形だけを狙い撃ちできます。
%[text] さらに、角括弧の中にカンマ `,` を使うと「または（OR）」という意味になります。たとえば、`[F,Cl,Br,I]` と書けば「フッ素、塩素、臭素、ヨウ素のどれか（＝任意のハロゲン）」という意味のクエリになります。また、元素記号の代わりに `#7` のようにシャープと原子番号を使うことで、芳香族（アルファベット小文字の `n`）か脂肪族（大文字の `N`）かを問わず、「とにかくすべての窒素原子」を指定することも可能です。これらの記号を組み合わせることで、どんなに複雑な官能基のブレンドであっても、ピンポイントで表現できるようになります。
%[text] 付属の `pains.csv` には Baell & Holloway 2010 の SMARTS パターンが収録されています。
% NOTE: NumHeaderLines=1, VariableNamesLine=1 must be specified explicitly because
% SMARTS patterns contain commas inside quoted fields which confuse MATLAB's
% auto-detection (detectImportOptions infers 13 columns instead of 4).
pains_data = readtable("data/list/pains.csv", TextType="string", ...
    Delimiter=",", NumHeaderLines=1, VariableNamesLine=1);
logInfo("PAINS SMARTS パターン %d 件を読み込んだ", height(pains_data)); %[output:48b4b33e]
pains_smarts = pains_data.SMARTS;   % SMARTS column by name (robust vs column-index)

%[text] テスト分子
test_smiles = ["CCO",                             % エタノール -- クリーン
               "CC(=O)Oc1ccccc1C(=O)O",          % アスピリン -- クリーン
               "O=C1C=C(O)C(=CC1=O)c1ccccc1",    % キノン（既知 PAINS）
               "c1ccc(cc1)N=Nc2ccccc2"];          % アゾ染料（既知 PAINS）
test_names  = ["エタノール", "アスピリン", "キノン", "アゾ染料"];

%[text] **補足（パフォーマンス）**: このループは `numel(test_smiles) × n_patterns` 回の個別 IPC 呼び出しを行います。全 479 パターン × 4 分子 = 約 60 秒かかります（Python IPC オーバーヘッドのため）。大規模な PAINS フィルタリングには、セル配列バッチ呼び出しまたは SMARTS OR 結合による高速化を検討しましょう。
logInfo("--- PAINS スクリーニング（全 %d パターン）---", height(pains_data)); %[output:65443a13]
n_patterns = height(pains_data);   % Use all 479 patterns (quinone matches row 205, azo matches row 469)

for i = 1:numel(test_smiles) %[output:group:84d44a38]
    if ~emk.mol.isValid(test_smiles(i)); continue; end
    mol  = emk.mol.fromSmiles(test_smiles(i));
    hits = 0;
    for p = 1:n_patterns
        if emk.mol.hasSubstruct(mol, pains_smarts(p))
            hits = hits + 1;
        end
    end
    status = "クリーン";
    if hits > 0; status = sprintf("PAINS（%d ヒット）", hits); end
    logInfo("  %-14s  %s", test_names(i), status); %[output:74566d7e]
end %[output:group:84d44a38]
%%
%[text] ## セクション 5: カスタム官能基テーブルを構築する
%[text] ### **実験を台無しにする「問題児」をあらかじめ排除する**
%[text] 創薬の歴史において、細胞を使った実験で「ものすごく効いているように見えるのに、いざ詳しく調べてみると、ただ細胞を傷つけていただけだったり、測定機器の光を邪魔していただけだったりする偽物の化合物」が数多く発見され、研究者を悩ませてきました。こうした、実験（アッセイ）において高確率で偽陽性のエラーを引き起こすお騒がせな特定の部分構造のことを **PAINS（Pan-Assay Interference Compounds：パンス）** と呼びます。
%[text] 実際の創薬プロジェクトでは、無駄な実験コストや時間を省くため、コンピューター上で化合物ライブラリに対して「PAINSのSMARTSフィルター」をかけ、これらの問題児を最初の段階でバッサリとゴミ箱に捨てる処理（カウンタースクリーニング）が必須の常識となっています。ここでは、そうしたお騒がせ構造の代表格である「キノン骨格」や「ロダニン骨格」を模した SMARTS パターンを使って、手元の化合物が安全なものかどうかを厳格にチェックする実務さながらのワークフローを試してみましょう。
data = readtable("data/list/everyday_chemicals.csv", TextType="string");
fg_patterns = {"c1ccccc1", "[OH]", "C(=O)[OH]", "[NH,NH2]", "[C](=O)[N]"};
fg_names    = ["ベンゼン環", "ヒドロキシル", "カルボン酸", "アミン", "アミド"];

n = height(data);
fg_matrix = false(n, numel(fg_patterns));

for i = 1:n
    if ~emk.mol.isValid(data.SMILES(i)); continue; end
    mol = emk.mol.fromSmiles(data.SMILES(i));
    for p = 1:numel(fg_patterns)
        fg_matrix(i, p) = emk.mol.hasSubstruct(mol, fg_patterns{p});
    end
end

fg_tbl = array2table(fg_matrix, "VariableNames", fg_names);
fg_tbl.Name = data.CommonName;
fg_tbl = movevars(fg_tbl, "Name", "Before", "ベンゼン環");
disp(fg_tbl); %[output:73316fe1]
%[text] **✏️ やってみよう 2 — ハロゲン列を追加してみましょう**
%[text] SMARTS `"[F,Cl,Br,I]"` を使って `fg_tbl` に `"ハロゲン"` 列を追加してみましょう。
%[text] ハロゲン原子を含む日用化学品はどれですか？
%[text] ヒント: `fg_patterns` と `fg_names` に新しいエントリを追加して、セクション 5 の先頭から再実行します。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: まとめ
%[text] 今回のチュートリアルで学んだ部分構造検索の API と、実務でよく使われる代表的な SMARTS クイックリファレンスです。
%[text] このセクションで学んだ主な API:
%[text] - `emk.mol.hasSubstruct(mol, "c1ccccc1")` — 単一 Mol、SMARTS クエリ
%[text] - `emk.mol.hasSubstruct(mols, "[OH]")` — セル配列、バッチモード \
%[text] **SMARTS クイックリファレンス:**
%[text]     ベンゼン環   : `c1ccccc1`  ← 6員全炭素芳香環のみ（ヘテロ芳香環は不検出）
%[text]     ヒドロキシル : `[OH]`
%[text]     カルボン酸   : `C(=O)[OH]`
%[text]     アミン（1/2）: `[NH,NH2]`
%[text]     アミド       : `C(=O)[N]`
%[text]     任意のハロゲン: `[F,Cl,Br,I]`
%[text]     任意の窒素   : `[#7]`
%[text]     任意の炭素   : `[#6]`
%[text] 
%[text] **次回 F06** では、分子ファイル（SDF・SMILES リスト）の読み書き方法を学びます。
%[text] 複数の分子をまとめてインポート・エクスポートできるようになります。
%%
%[text] ## 演習
%[text] 各演習は `answers/f05_answers.m` を参照する前に自分で解いてみましょう。
%[text] 
%[text] **E1.** SMILES 文字列を受け取り、含まれる官能基を人間が読みやすいサマリーとして表示するスクリプトを書きましょう。
%[text] 対象官能基: 芳香環、水酸基、カルボン酸、アミン、アミド、
%[text] ハロゲン（`"[F,Cl,Br,I]"`）、ケトン（`"C(=O)[#6]"`）。
%[text] ニコチン `"CN1CCC[C@H]1C2=CN=CC=C2"` と モルヒネ
%[text] `"OC1=CC2=C(CC3N(CC23)CC4=CC=CC=C14)C=C1"` でテストしてください。
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:444a3b00]
%   data: {"dataType":"text","outputData":{"text":"[09:25:55][INFO]  initPython: Python 3.10 already active (Status: Loaded) -- skipping.\n","truncated":false}}
%---
%[output:186a4c5b]
%   data: {"dataType":"text","outputData":{"text":"[09:25:55][INFO]  F05: セットアップ完了\n","truncated":false}}
%---
%[output:42b5f0c9]
%   data: {"dataType":"text","outputData":{"text":"[09:25:56][INFO]  アスピリンはベンゼン環を含む: 1\n","truncated":false}}
%---
%[output:65242d3e]
%   data: {"dataType":"text","outputData":{"text":"[09:25:56][INFO]  エタノールはベンゼン環を含む: 0\n","truncated":false}}
%---
%[output:8b36612b]
%   data: {"dataType":"text","outputData":{"text":"[09:25:56][INFO]  アスピリンは -COOH を含む: 1\n","truncated":false}}
%---
%[output:50c7289a]
%   data: {"dataType":"text","outputData":{"text":"[09:25:56][INFO]  アセトアミノフェンはアミドを含む: 1\n","truncated":false}}
%---
%[output:24f38584]
%   data: {"dataType":"text","outputData":{"text":"[09:25:56][INFO]  官能基スクリーニング結果:\n","truncated":false}}
%---
%[output:08a7c914]
%   data: {"dataType":"text","outputData":{"text":"                       Aspirin    Benzene    Ethanol    Caffeine    Acetaminophen\n                       _______    _______    _______    ________    _____________\n\n    ベンゼン環           true       true       false      false          true     \n    カルボン酸           true       false      false      false          false    \n    ヒドロキシル(-OH)     true       false      true       false          true     \n    アミド(-CONH)        false      false      false      false          true     \n    窒素原子             false      false      false      true           true     \n\n","truncated":false}}
%---
%[output:48b4b33e]
%   data: {"dataType":"text","outputData":{"text":"[09:25:56][INFO]  PAINS SMARTS パターン 479 件を読み込んだ\n","truncated":false}}
%---
%[output:65443a13]
%   data: {"dataType":"text","outputData":{"text":"[09:25:56][INFO]  --- PAINS スクリーニング（全 479 パターン）---\n","truncated":false}}
%---
%[output:74566d7e]
%   data: {"dataType":"text","outputData":{"text":"[09:26:04][INFO]    エタノール           クリーン\n[09:26:12][INFO]    アスピリン           クリーン\n[09:26:20][INFO]    キノン             PAINS（2 ヒット）\n[09:26:28][INFO]    アゾ染料            PAINS（1 ヒット）\n","truncated":false}}
%---
%[output:73316fe1]
%   data: {"dataType":"text","outputData":{"text":"          Name          ベンゼン環    ヒドロキシル    カルボン酸    アミン    アミド\n    ________________    _________    __________    _________    _____    _____\n\n    \"caffeine\"            false        false         false      false    false\n    \"nicotine\"            false        false         false      false    false\n    \"theobromine\"         false        false         false      false    false\n    \"aspirin\"             true         true          true       false    false\n    \"acetaminophen\"       true         true          false      true     true \n    \"ibuprofen\"           true         true          true       false    false\n    \"salicylic acid\"      true         true          true       false    false\n    \"ethanol\"             false        true          false      false    false\n    \"methanol\"            false        true          false      false    false\n    \"isopropanol\"         false        true          false      false    false\n    \"acetone\"             false        false         false      false    false\n    \"glycerol\"            false        true          false      false    false\n    \"sucrose\"             false        true          false      false    false\n    \"glucose\"             false        true          false      false    false\n    \"fructose\"            false        true          false      false    false\n    \"acetic acid\"         false        true          true       false    false\n    \"citric acid\"         false        true          true       false    false\n    \"lactic acid\"         false        true          true       false    false\n    \"benzoic acid\"        true         true          true       false    false\n    \"formic acid\"         false        true          true       false    false\n    \"vanillin\"            true         true          false      false    false\n    \"capsaicin\"           true         true          false      true     true \n    \"menthol\"             false        true          false      false    false\n    \"limonene\"            false        false         false      false    false\n    \"eugenol\"             true         true          false      false    false\n    \"benzaldehyde\"        true         false         false      false    false\n    \"linalool\"            false        true          false      false    false\n    \"ascorbic acid\"       false        true          false      false    false\n    \"urea\"                false        false         false      true     true \n    \"carvone\"             false        false         false      false    false\n\n","truncated":false}}
%---
