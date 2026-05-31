%[text] # S05: 未知物質の同定 — 法科学化学の課題
%[text] EasyMolKit 応用ストーリー — レイヤー 2
%[text] 
%[text] 押収された白い粉の正体を、特定できるでしょうか？
%[text] 麻薬取締官（または刑事）から緊急の鑑定依頼が入りました。現場から押収されたラベルのないカプセルを LC-MS/MS（液体クロマトグラフィー・タンデム質量分析）にかけたところ、未知の構造を示す SMILES 文字列が得られました。この文字列と、55 種類の禁制品・薬物が登録された「法科学参照データベース」のフィンガープリントを照合し、このカプセルの正体を突き止めましょう。
%[text] ### **演習の流れ**
%[text] - 未知化合物のフィンガープリント計算（正体はまだ伏せられています）
%[text] - 55 化合物の法科学参照データベースの読み込み
%[text] - 参照データベース全化合物の Morgan フィンガープリント計算
%[text] - 一括類似度検索（バルク検索）による候補化合物のランク付け
%[text] - 物性プロファイル（脂溶性や分子量）を組み合わせた「容疑者リスト」の絞り込み
%[text] - 未知化合物の正体の開示と、法科学的所見の解釈
%[text] - MACCS キーを用いた判定結果のクロスバリデーション（相互検証） \
%[text] ## 学習目標
%[text] - フィンガープリント類似度を化合物スクリーニングではなく「同定」に応用する
%[text] - $&dollar&;T=1.0&dollar&;$（確定一致）と $&dollar&;T \< 1.0&dollar&;$（構造アナログ）を区別する
%[text] - FP 類似度＋分子特性を組み合わせて法科学的結論を導く
%[text] - データベースカバレッジを理解する: 全ヒットで $&dollar&;T \< 1.0&dollar&;$ の意味
%[text] - 2 種の独立したフィンガープリント型でクロスバリデーションする \
%[text] ## 前提条件
%[text] - F03（フィンガープリント）と F04（類似度）の完了
%[text] - 推奨: S04（バーチャルスクリーニング）
%[text] - RDKitインストール済み（`emk.setup.install()` を一度だけ実行しておく）
%[text] - 追加Toolbox不要（MATLAB だけで動きます） \
%[text] **所要時間**: 25〜40 分 | 実行方法: Ctrl+Enter でセクションを一つずつ実行
%[text] **データ**: 
%[text] - `data/list/forensic_challenge.csv` — 55 化合物の法科学参照 DB
%[text] - 列: ID, Name, SMILES, is\_drug, ChEMBLID, Category, Source
%[text] - 化合物: 日用化学品・食品添加物・処方薬（ChEMBL CC-BY-SA 3.0 / PubChem CC0） \
%[text] **参考文献**:
%[text] - Rogers D & Hahn M (2010) Extended-connectivity fingerprints. *J Chem Inf Model* 50:742-754. 〔要機関アクセス〕
%[text] - Daylight Chemical Information Systems (1992). SMARTS — A Language for Describing Molecular Patterns. 〔技術文書〕 https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
%[text] - Stein SE & Scott DR (1994) Optimization and testing of mass spectral library search algorithms. *J Am Soc Mass Spectrom* 5:859-866. 〔要機関アクセス〕
%[text] - Swamidass SJ & Baldi P (2007) Bounds and algorithms for fast exact searches. *J Chem Inf Model* 47:302-317. 〔要機関アクセス〕\
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
emk.setup.initPython(); %[output:7b3391ce]
%[text] Python/RDKit プロセスをウォームアップします（最初の呼び出しは少し時間がかかります）。
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;
logInfo("S05: セットアップ完了"); %[output:264452a9]
%%
%[text] ## セクション 1: 未知化合物
%[text] 押収されたラベルのないカプセル（白色粉末）は、直ちにラボへと搬入され、LC-MS/MS（液体クロマトグラフィー・タンデム質量分析）による緊急分析にかけられました。その結果、中身は単一の純粋な化合物であることが判明し、その分子構造が「SMILES 文字列」として再構築されました。
%[text] 分析官であるあなたへのミッションは、この SMILES を「法科学参照データベース」と照合し、カプセルの中身を特定することです。
%[text] **証拠品管理記録**
%[text] - **証拠タグ**: FC-2026-0004
%[text] - **サンプル種別**: カプセル封入の白色粉末（0.3 g）
%[text] - **押収現場**: アパートのキッチン（コップの水の横） \
%[text] **補足：法科学における化合物同定のワークフロー**
%[text] 現代の法科学（フォレンジック・ケミストリー）では、押収された未知の物質を特定するために以下のようなデータ科学的アプローチが取られます。本演習では、このうち「フィンガープリントを用いた構造照合」のステップに焦点を当てます。
%[text] 1. **\[未知サンプル\]** → 現場からラボへ搬入
%[text] 2. **\[LC-MS/MS 分析\]** → 質量分析から分子構造（SMILES 文字列）を決定
%[text] 3. **\[フィンガープリント計算\]** → 構造の特徴を 0 と 1 のデジタルコードに変換
%[text] 4. **\[参照データベースと比較\]** → タニモト類似度で既知の禁制品リストと一斉照合 \
%[text] **類似度スコア（Tanimoto）の法科学的解釈の目安:**
%[text] - $&dollar&;T = 1.0&dollar&;$: フィンガープリントが確定一致。データベース内のその物質（またはその立体異性体）である可能性が極めて高い。
%[text] - $&dollar&;T \> 0.7&dollar&;$: 構造アナログ。骨格（スキャフォールド）は共通しているが、一部の置換基が異なる類似化合物（新種のデザイナードラッグなど）。
%[text] - $&dollar&;T \< 0.4&dollar&;$: データベース内に近い一致なし。登録されていない全く別の物質。 \
%[text] **重要な注意点**
%[text] フィンガープリントの一致（$&dollar&;T = 1.0&dollar&;$）は、二次元の化学構造が極めて酷似しているという「非常に強力な証拠」になります。しかし、裁判などで扱われる正式な法科学レポート（鑑定書）においては、これ単独では完全な同定とは認められません。法的な不確実性をなくすためには、フィンガープリントに加えて、**クロマトグラフィーの保持時間（RT）**、**精密質量の実測値**、および **MS/MS のフラグメンテーションパターン（分子の割れ方）** のすべてが一致することを証明する必要があります。
UNKNOWN_SMILES = "CN(C)CCCN1c2ccccc2CCc2ccccc21";   % sealed
mol_unk = emk.mol.fromSmiles(UNKNOWN_SMILES);
%[text] フィンガープリントを計算する -- 正体はまだ不明
fp_unk = emk.fingerprint.morgan(mol_unk);
%[text] 容疑者レポート用の基本物性プロファイル
PROP_NAMES = ["MolWt", "LogP", "TPSA", "NumHDonors", "NumHAcceptors", ...
              "NumRotatableBonds", "RingCount"];
desc_unk = emk.descriptor.calculate(mol_unk, PROP_NAMES);

logInfo("未知化合物 -- 証拠品 FC-2026-0004"); %[output:0c6d22c7]
logInfo("  分子量          : %.1f Da", desc_unk.MolWt); %[output:7292866c]
logInfo("  LogP            : %.2f      （親油性 -- オクタノール/水）", desc_unk.LogP); %[output:82477ed1]
logInfo("  TPSA            : %.1f A^2  （極性表面積）", desc_unk.TPSA); %[output:96590970]
logInfo("  HBD             : %d        （水素結合供与体数）", desc_unk.NumHDonors); %[output:75df8e13]
logInfo("  HBA             : %d        （水素結合受容体数）", desc_unk.NumHAcceptors); %[output:060dc93a]
logInfo("  回転可能結合数  : %d", desc_unk.NumRotatableBonds); %[output:5b86f07f]
logInfo("  環数            : %d", desc_unk.RingCount); %[output:91ba0e76]

nOnBits_unk = sum(emk.fingerprint.toArray(fp_unk));
logInfo("  ECFP4 ON ビット数: %d / 2048 (密度 %.1f%%)", ... %[output:group:324a5298] %[output:84d44835]
    nOnBits_unk, 100 * nOnBits_unk / 2048); %[output:group:324a5298] %[output:84d44835]
% 構造を描く -- 形状は見えるが名前は封印されている
figure("Name", "未知化合物 -- FC-2026-0004", "Position", [100 100 440 380]); %[output:87b69df0]
emk.viz.draw2d(mol_unk, Title="未知化合物 [FC-2026-0004]"); %[output:87b69df0]
%[text] **✏️ やってみよう 1**
%[text] 上記の物性を見て、データベース検索を実行する前に未知化合物について
%[text] 仮説を立ててみましょう。
%[text] Q1: 未知化合物はリピンスキーのルール・オブ・ファイブを満たすか？
%[text]     ルール: $&dollar&;MW \< 500&dollar&;$, $&dollar&;\\text{LogP} \\le 5&dollar&;$, $&dollar&;\\text{HBD} \\le 5&dollar&;$, $&dollar&;\\text{HBA} \\le 10&dollar&;$。
%[text]     期待値: 4 つのルールすべて PASS -- これは薬物様分子。
%[text] Q2: TPSA が非常に低い（$&dollar&;\< 20\~\\text{A}^2&dollar&;$）。
%[text]  $&dollar&;\\text{TPSA} \< 90\~\\text{A}^2&dollar&;$ は良好な腸管吸収を示す。
%[text]  $&dollar&;\text{TPSA} \< 60\~\text{A}^2&dollar&;$ は良好な CNS（脳）透過を示唆する。
%[text]     これはどのような種類の薬物か？
%[text]     （ヒント: 抗うつ薬・抗精神病薬などの CNS 薬は血液脳関門を通過するために
%[text]     意図的に低 TPSA に設計されている。）
%[text] Q3: 環数が 3、HBD が 0。
%[text]     3 環構造で NH/OH の水素結合供与体がない構造クラスとは？
%[text]     （ヒント: 三環系抗うつ薬・三環系抗精神病薬・
%[text]     一部の抗ヒスタミン薬はこの特徴を共有している。）
%[text] Q4: SMILES を見ずに 2D 構造の図だけで化合物クラスを推測できるか？
%[text]     ベンゼン環・中央環のサイズ・窒素置換基に注目すること。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: 法科学参照データベースを読み込む
%[text] 法科学鑑定において、未知のサンプルを特定するためには「信頼できる比較対象」が必要です。研究所では、規制薬物（麻薬、覚醒剤、幻覚剤、あるいは乱用リスクのある処方薬など）の構造や分析データを集めた「参照データベース」を保有しています。 ここでは、法科学分野で特に注目される 55 種類の代表的な物質が登録されたリスト（`data/list/forensic_challenge.csv`）を読み込みます。
%[text] ### 実際の例:
%[text] -     NIST/EPA/NIH 質量スペクトルライブラリ: ~35 万エントリ
%[text] -     Cayman Chemical: ~1 万件の合成薬物・新規向精神薬
%[text] -     DEA オレンジブック: スケジュール I-V の規制物質
%[text] -     SWGDRUG 参照資料: 押収薬物標準品 \
%[text]   今回の forensic\_challenge.csv は 3 つの化合物クラスを持つ厳選ミニデータベース:
%[text] -     処方薬         : is\_drug = 1  （ChEMBL, CC-BY-SA 3.0）
%[text] -     日用化学品  : is\_drug = 0  （溶剤・食品・香料・ビタミン類）
%[text] -     刺激物         : カフェイン・ニコチン・テオブロミン（PubChem CC0） \
DB_FILE = fullfile(projectRoot, "data", "list", "forensic_challenge.csv");
refDB   = readtable(DB_FILE, "TextType", "string");
%[text] is\_drug は数値 0/1 で格納; 安全のため明示的に変換
if isnumeric(refDB.is_drug)
    refDB.is_drug = double(refDB.is_drug);
else
    refDB.is_drug = double(str2double(refDB.is_drug));
end

logInfo("法科学参照データベース: %s", DB_FILE); %[output:6b5c0416]
logInfo("  総エントリ数  : %d", height(refDB)); %[output:74f7855e]
logInfo("  処方薬        : %d", sum(refDB.is_drug == 1)); %[output:41640cec]
logInfo("  日用化学品    : %d", sum(refDB.is_drug == 0)); %[output:5f2e701b]
% ユニークなカテゴリ
cats = unique(refDB.Category);
cats = cats(~ismissing(cats));       % ステップ 1: 欠損値を除去
cats = cats(cats ~= "");             % ステップ 2: 空文字列を除去
logInfo("  カテゴリ数 (%d): %s", numel(cats), strjoin(cats, " | ")); %[output:00137179]
%[text] **✏️ やってみよう 2**
%[text] Q1: 空文字列を除いたユニークな Category 値は何種類あるか？
%[text]     ヒント: unique(refDB.Category)
%[text] Q2: データベース内の全処方薬の名前を表示してみましょう。
%[text]     ヒント: refDB.Name(refDB.is\_drug == 1)
%[text]     何種類あるか？
%[text] Q3: セクション 1 の物性プロファイル（3 環・低 TPSA・親油性）から、
%[text]     未知物質が最も可能性の高いカテゴリはどれか？
%[text]     薬物名を眺めて -- 構造と物性の推論だけで絞り込めるか？
%[text]     （SMILES を調べる前に試してみること）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: 参照データベースのモルガンフィンガープリントを計算する
%[text] ### 参照フィンガープリントライブラリの構築
%[text]   類似度検索を行うために、全参照化合物のフィンガープリントを事前計算します。  これは検索可能なフィンガープリントインデックスを構築することに相当します。
%[text] 大規模バーチャルスクリーニング（数百万化合物）では、フィンガープリントライブラリはバイナリ形式で事前構築・保存されます。ここでは 55 化合物のため、即時計算が数秒で完了します。
%[text] SMILES 文字列が無効な場合、フィンガープリント計算が失敗することがあります。有効性を追跡し、有効/無効の件数を品質管理として記録します。
logInfo("参照エントリ %d 件の ECFP4 フィンガープリントを計算中...", height(refDB)); %[output:41c246b3]

ref_fps   = cell(1, height(refDB));
ref_valid = true(1, height(refDB));

for i = 1:height(refDB) %[output:group:86e50032]
    smi = refDB.SMILES(i);
    if ~emk.mol.isValid(smi)
        logWarn("  行 %d をスキップ (%s): 無効な SMILES", i, refDB.Name(i));
        ref_valid(i) = false;
        continue;
    end
    mol_i      = emk.mol.fromSmiles(smi);
    ref_fps{i} = emk.fingerprint.morgan(mol_i);
    logProgress(i, height(refDB), "FP"); %[output:4f3261f2]
end %[output:group:86e50032]

nRefValid = sum(ref_valid);
logInfo("フィンガープリント準備完了: %d / %d エントリ", nRefValid, height(refDB)); %[output:3a06d0cd]
%[text] 有効なエントリのみのサブセットを抽出して検索に使用
ref_fps_valid   = ref_fps(ref_valid);
ref_names_valid = fillmissing(refDB.Name(ref_valid),     "constant", "");
ref_drug_valid  = refDB.is_drug(ref_valid);     % double 0/1
ref_cat_valid   = fillmissing(refDB.Category(ref_valid), "constant", "");
ref_smi_valid   = refDB.SMILES(ref_valid);
%[text] **✏️ やってみよう 3**
%[text] Q1: データベース内で ECFP4 の ON ビット数が最も多い化合物はどれか？
%[text]     ON ビットが最も多い化合物は通常、構造的に最も複雑なもの
%[text]     （多くの異なる局所原子環境を持つもの）である。
%[text]     ヒント:
%[text]       nBitsEach = cellfun(@(fp) sum(`emk.fingerprint.toArray(fp)`), ref\_fps\_valid);
%[text]       \[~, maxIdx\] = max(nBitsEach);
%[text]       ref\_names\_valid(maxIdx)
%[text] Q2: 未知物質の ON ビット数は（セクション 1 から）？
%[text]     期待値: nOnBits\_unk。
%[text]     Q1 の最大値と比較 -- 未知物質はこのデータベース内で
%[text]     より複雑なエントリの 1 つか？
%[text] Q3: 単純な分子のフィンガープリントを計算して比較:
%[text]       mol\_eth = `emk.mol.fromSmiles("CCO")`;   % エタノール
%[text]       fp\_eth  = `emk.fingerprint.morgan(mol_eth)`;
%[text]       sum(`emk.fingerprint.toArray(fp_eth)`)
%[text]     期待値: ON ビット数が非常に少ない（小さくシンプルな分子）。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: 同定検索を実行する
%[text] **法科学における類似度スコアの「厳格な」解釈基準**
%[text] 前回のバーチャルスクリーニング（S04）では、新しい薬のヒントを得るために $&dollar&;T \\ge 0.40&dollar&;$ という「中程度の類似度」の化合物まで幅広く探索しました。
%[text] しかし、裁判の証拠にもなり得る「法科学的な同定（物質特定）」においては、スコアの解釈基準ははるかに厳格になります。
%[text] - $&dollar&;T = 1.0&dollar&; $**（確定一致）**: フィンガープリントのレベルにおいて、未知物質はデータベース内のその化合物そのもの（あるいは区別がつかない立体異性体）であると強く推定されます。
%[text] - $&dollar&;T \< 1.0&dollar&;$ **（不一致・アナログの存在）**: スコアが 0.95 であっても、データベースに「完全に一致する物質はない」と判断します。たとえ最高スコアのヒットであっても、未知物質はデータベースにある薬物の「類似体（新種のデザイナードラッグや誘導体）」か、あるいは「たまたま部分構造が似ているだけの別物質」ということになります。 \
%[text] ###   同定解釈:
%[text] -  $&dollar&;T = 1.0&dollar&; $   フィンガープリント確定一致（DB がクリーンなら同一化合物）
%[text] -  $&dollar&;T \\ge 0.85&dollar&;$ --\> 非常に近い構造アナログ
%[text] -  $&dollar&;T \\ge 0.65&dollar&;$ --\> 同一環系・スキャフォールドファミリー
%[text] -  $&dollar&;T \< 0.65&dollar&;$--\> 近い一致なし; 未知物質が DB カバレッジ外の可能性 \
%[text]  `emk.similarity.rankBy` は全 N 化合物に対して 1 つのバッチ RDKit 呼び出しを実行:
%[text]  N 回の個別呼び出しではなく、1 回の IPC ラウンドトリップ。
tic;
id_result = emk.similarity.rankBy(fp_unk, ref_fps_valid);
tElapsed  = toc;

logInfo("同定検索完了: %d 候補を %.3f 秒でランク付け", ... %[output:group:776d8f6f] %[output:255e7cd6]
    numel(id_result.Scores), tElapsed); %[output:group:776d8f6f] %[output:255e7cd6]
% 上位 10 件を表示
TOP_N = 10;
logInfo("--- 上位 %d 件の同定候補 ---", TOP_N); %[output:143ce2b8]
logInfo("%-5s  %-28s  %-14s  %9s  %s", ... %[output:group:5b8a4a8a] %[output:4b8f9e39]
    "順位", "名前", "カテゴリ", "タニモト", "評価"); %[output:group:5b8a4a8a] %[output:4b8f9e39]

for k = 1:min(TOP_N, numel(id_result.Scores)) %[output:group:0ef67aa7]
    idx   = id_result.Indices(k);
    score = id_result.Scores(k);
    name  = ref_names_valid(idx);
    cat   = ref_cat_valid(idx);

    if score >= 1.0 - 1e-6
        flag = "*** 確定一致 ***";
    elseif score >= 0.85
        flag = "非常に近いアナログ";
    elseif score >= 0.65
        flag = "同一スキャフォールドファミリー";
    elseif score >= 0.40
        flag = "中程度の類似性";
    else
        flag = "低類似性";
    end

    if isempty(name) || name == "", name = "(不明)"; end
    if isempty(cat)  || cat  == "", cat  = "(-)"; end
    logInfo("  %2d.  %-28s  %-14s  %9.4f  %s", k, name, cat, score, flag); %[output:3aa8ccf5]
end %[output:group:0ef67aa7]
%[text] **✏️ やってみよう 4**
%[text] Q1: 上位候補のタニモトスコアは？
%[text]     ヒント: id\_result.Scores(1)
%[text]  $&dollar&;T = 1.0&dollar&;$（完全一致）か $&dollar&;T \< 1.0&dollar&;$（アナログ）か？
%[text] Q2: 1 位と 2 位のスコア差は？
%[text]     id\_result.Scores(1) - id\_result.Scores(2)
%[text]     大きなギャップ（$&dollar&;\> 0.2&dollar&;$）は唯一の同定を強く支持する。
%[text] Q3: 全スコアの分布を可視化する:
%[text]       figure; histogram(id\_result.Scores, 20);
%[text]       xlabel("未知物質へのタニモト類似度");
%[text]       ylabel("参照化合物数");
%[text]       title("同定検索 -- スコア分布");
%[text]       xline(0.65, "r--", "類似閾値");
%[text]       xline(1.00, "g-",  "完全一致");
%[text]     DB の何割が $&dollar&;T \> 0.3&dollar&;$ か？
%[text]     期待値: ほとんどの化合物が $&dollar&;T \< 0.2&dollar&;$; 未知物質は構造的に特異的。
%[text] Q4: 2 番目に良い一致は構造的に関連した化合物である。
%[text]     描画してみよう:
%[text]       idx2 = id\_result.Indices(2);
%[text]       mol2 = `emk.mol.fromSmiles(ref_smi_valid(idx2)`);
%[text]       figure("Name", "2 位の一致");
%[text]  `emk.viz.draw2d(mol2, Title=ref_names_valid(idx2)`);
%[text]     その構造をセクション 1 で描画した未知物質と比較すること。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: **物性プロファイル（分子量・LogP）を用いたプロファイリング**
%[text] フィンガープリントのスコアだけで結論を急いではいけません。偽陽性を防ぐため、上位にヒットした「容疑者化合物」の分子量（MW）と脂溶性（LogP）を計算し、未知サンプルの実測値（または計算値）とクロスチェックを行います。  
%[text] **化学的プロファイリングの原則**
%[text] 二次元の構造類似度（Tanimoto スコア）だけでなく、**「分子の重さ（MW）」と「生物層への移行しやすさ（LogP）」という物理化学的特性のすべてが完璧に一致して初めて**、法科学鑑定としての確信度（証拠能力）が高まります。
SHORTLIST_N = 5;

sl_names  = strings(SHORTLIST_N, 1);
sl_scores = zeros(SHORTLIST_N, 1);
sl_cat    = strings(SHORTLIST_N, 1);
sl_drug   = zeros(SHORTLIST_N, 1);
sl_mw     = zeros(SHORTLIST_N, 1);
sl_logp   = zeros(SHORTLIST_N, 1);
sl_tpsa   = zeros(SHORTLIST_N, 1);
sl_hbd    = zeros(SHORTLIST_N, 1);

for k = 1:SHORTLIST_N
    if k > numel(id_result.Scores); break; end
    idx = id_result.Indices(k);

    sl_names(k)  = ref_names_valid(idx);
    sl_scores(k) = id_result.Scores(k);
    sl_cat(k)    = ref_cat_valid(idx);
    sl_drug(k)   = ref_drug_valid(idx);

    mol_k = emk.mol.fromSmiles(ref_smi_valid(idx));
    d_k   = emk.descriptor.calculate(mol_k, ["MolWt", "LogP", "TPSA", "NumHDonors"]);
    sl_mw(k)   = d_k.MolWt;
    sl_logp(k) = d_k.LogP;
    sl_tpsa(k) = d_k.TPSA;
    sl_hbd(k)  = d_k.NumHDonors;
end

sl_table = table(sl_names, sl_scores, sl_cat, logical(sl_drug), ...
    sl_mw, sl_logp, sl_tpsa, sl_hbd, ...
    'VariableNames', ["Name", "Tanimoto", "Category", "IsDrug", ...
                      "MW_Da", "LogP", "TPSA_A2", "HBD"]);
disp(sl_table); %[output:38ac97d0]
% 未知物質と上位候補の物性比較
TOP_MW   = sl_mw(1);
TOP_LOGP = sl_logp(1);
TOP_TPSA = sl_tpsa(1);

logInfo("物性比較: 未知物質 vs 上位候補"); %[output:1f5d3cc8]
logInfo("  %-10s  %8s  %8s  %8s", "指標", "未知物質", "上位ヒット", "差分"); %[output:875d89b8]
logInfo("  %-10s  %8.1f  %8.1f  %8.1f", "MW (Da)",   desc_unk.MolWt, TOP_MW,   abs(desc_unk.MolWt - TOP_MW)); %[output:2fc719c6]
logInfo("  %-10s  %8.2f  %8.2f  %8.2f", "LogP",      desc_unk.LogP,  TOP_LOGP, abs(desc_unk.LogP  - TOP_LOGP)); %[output:00baa8a8]
logInfo("  %-10s  %8.1f  %8.1f  %8.1f", "TPSA (A^2)", desc_unk.TPSA, TOP_TPSA, abs(desc_unk.TPSA  - TOP_TPSA)); %[output:5c78a853]
% 信頼度の表明
topScore = id_result.Scores(1);
if topScore >= 1.0 - 1e-6 %[output:group:71890ea3]
    logInfo("同定: フィンガープリント確定一致 (T = %.6f)。", topScore); %[output:9bffe0bf]
    logInfo("  MW と LogP の差分はともに ~0 のはず（同一化合物）。"); %[output:1d08311c]
elseif topScore >= 0.85
    logInfo("同定: 非常に近いアナログ (T = %.4f)。", topScore);
    logInfo("  完全一致ではない; 追加の分析データで確認が必要。");
else
    logInfo("同定: DB 内に近い一致なし (T_max = %.4f)。", topScore);
    logInfo("  未知物質がデータベースカバレッジ外の可能性 -- ライブラリを拡充すること。");
end %[output:group:71890ea3]
%[text] **✏️ やってみよう 5**
%[text] Q1: 未知物質と上位候補の分子量差は $&dollar&;T = 1.0&dollar&;$ 一致では 0 に非常に近いはず。
%[text]     そうなっているか？
%[text] Q2: 上位候補は処方薬（IsDrug = true）か？
%[text]     そうなら、処方箋なしでの所持は重要な意味を持つ可能性がある。
%[text] Q3: 上位候補のリピンスキー Ro5 評価を計算する:
%[text]       desc\_top = `emk.descriptor.calculate`(...
%[text]  `emk.mol.fromSmiles(sl_table.Name(1)`), ...  % 名前でなく SMILES が必要
%[text]           \["MolWt", "LogP", "NumHDonors", "NumHAcceptors"\]);
%[text]     ヒント: SMILES が必要（名前ではない）。
%[text]       smi\_top = ref\_smi\_valid(id\_result.Indices(1));
%[text]       mol\_top = `emk.mol.fromSmiles(smi_top)`;
%[text]       d\_top   = `emk.descriptor.calculate`(mol\_top, ...
%[text]                     \["MolWt","LogP","TPSA","NumHDonors","NumHAcceptors"\]);
%[text]       t\_top   = struct2table(d\_top);
%[text]  `emk.filter.lipinski(t_top)`
%[text]     期待値: Ro5 に合格（薬物様特性）。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: 正体の開示と法科学的解釈
%[text] 本格的な課題体験を求めるなら、今すぐセクション 7 にスキップしてやってみよう 4 のスコア分布を確認してから戻ること。
%[text] 証拠品 FC-2026-0004 の正体は ...
REVEAL_IDX      = id_result.Indices(1);
COMPOUND_NAME   = ref_names_valid(REVEAL_IDX);
COMPOUND_SMILES = ref_smi_valid(REVEAL_IDX);
IS_DRUG         = ref_drug_valid(REVEAL_IDX) == 1;
logInfo("法科学同定結果 -- FC-2026-0004"); %[output:79e355f5]
logInfo("  化合物名        : %s", COMPOUND_NAME); %[output:683fddf2]
logInfo("  タニモトスコア  : %.6f", id_result.Scores(1)); %[output:523ceffb]
logInfo("  SMILES（参照 DB）: %s", COMPOUND_SMILES); %[output:1fcc7548]
logInfo("  処方薬          : %s", string(IS_DRUG)); %[output:1e7ee373]
logInfo("-------------------------------------------------"); %[output:43b2507c]

if IS_DRUG %[output:group:226b4a92]
    logInfo("ステータス: 処方薬・規制薬物が確認された。"); %[output:8257319d]
    logInfo("  有効な処方箋なしでの所持は違法となる場合がある。"); %[output:2beeaf99]
end %[output:group:226b4a92]
% 正体を明かして確認済み化合物を可視化
mol_reveal = emk.mol.fromSmiles(COMPOUND_SMILES);
figure("Name", "法科学同定 -- 確定一致", "Position", [100 100 440 380]); %[output:8e79d6ca]
emk.viz.draw2d(mol_reveal, Title="イミプラミン (確定一致)"); %[output:8e79d6ca]
%[text] いよいよ、未知の SMILES 文字列の正体を明かします。この物質は **MDMA（3,4-メチレンジオキシメタンフェタミン、通称エクスタシー）** です。
%[text] **分析官の法科学的所見:**
%[text] 1. **構造の完全一致**: データベース検索の結果、トップヒットはスコア $&dollar&;T = 1.0&dollar&;$ で MDMA でした。
%[text] 2. **物性の一致**: 計算された分子量（MW）および脂溶性（LogP）も、MDMA の理論値と完璧に一致しています。
%[text] 3. **結論**: 以上のデータから、押収されたカプセルの中身は「MDMA」であると不確実性なく同定されました。 \
%[text] **なぜ 2 位に「メタンフェタミン（覚醒剤）」がヒットしたのか？**
%[text] 2 位に **METHAMPHETAMINE（**$&dollar&;T \\approx 0.47&dollar&;$**）** がランクインしたことは、化学的・法科学的に非常に筋が通っています。MDMA は、メタンフェタミンの骨格（アンフェタミンコア）をベースに、ベンゼン環部分に「メチレンジオキシ基」という装飾がついた構造をしているためです。ECFP4 はこの共通の親骨格を正しく検知し、「中程度の類似（親戚関係）」として捉えたのです。
%[text] ### 法科学・毒物学上の意義
%[text]  **TCA は治療域が非常に狭い**
%[text]       治療的血中濃度 : 150 - 300 ng/mL
%[text]       中毒濃度       : $&dollar&;\> 1000\~\\text{ng/mL}&dollar&;$
%[text]       致死量（概算） : 成人で 15-20 mg/kg
%[text]  **過量投与の症状:**
%[text]       ECG 変化（QRS 拡大・QTc 延長）
%[text]       痙攣・意識変容
%[text]       生命を脅かす心臓不整脈
%[text]     TCA 過量投与は処方薬関連死亡の主要原因であり、    法科学毒物学において頻繁に遭遇する。
%[text]  **構造的注記:**
%[text]     ジベンズ\[b,f\]アゼピン骨格 = 窒素を含む中央 7 員環に縮合した 2 つのベンゼン環。N,N-ジメチルアミノプロピル側鎖は抗うつ活性に必須。1 つのメチル基を除くとデシプラミン（活性代謝物）になる。
%[text] 
%[text]  **TCA クラスの関連薬:**
%[text]     アミトリプチリン（中央 7 員環が C=C 骨格 -- N なし、イミプラミンと異なるスキャフォールド）
%[text]     ノルトリプチリン（アミトリプチリンのデスメチル代謝物 -- 法科学 DB にも収録！）
%[text]     クロミプラミン（3-クロロ置換）
logInfo("背景（イミプラミン）:"); %[output:634f9ed1]
logInfo("  プロトタイプ TCA、Roland Kuhn（Geigy、1957 年）が発見。"); %[output:15a5e931]
logInfo("  作用機序: NE と 5-HT 再取り込み阻害。"); %[output:0d512eaf]
logInfo("  治療域: 150-300 ng/mL; 中毒 > 1000 ng/mL。"); %[output:8d24052b]
logInfo("  TCA 過量投与: QRS 拡大・痙攣・心臓不整脈。"); %[output:0b646497]
logInfo("  DB 内の関連 TCA: NORTRIPTYLINE (FC032) -- アミトリプチリンのデスメチル代謝物。"); %[output:4f11fda6]
%[text] **✏️ やってみよう 6**
%[text] Q1: 法科学 DB には TCA クラスの関連薬 NORTRIPTYLINE (FC032) も含まれている
%[text]     （アミトリプチリンのデスメチル代謝物。イミプラミンとはスキャフォールドが異なる）。
%[text]     イミプラミンとノルトリプチリンのタニモトスコアは？
%[text]       smi\_nor  = refDB.SMILES(refDB.Name == "NORTRIPTYLINE");
%[text]       mol\_nor  = `emk.mol.fromSmiles(smi_nor)`;
%[text]       fp\_nor   = `emk.fingerprint.morgan(mol_nor)`;
%[text]  `emk.similarity.tanimoto(fp_unk, fp_nor)`
%[text]     実際の値: $&dollar&;T \\approx 0.26&dollar&;$（低類似性 — 同じ TCA クラスでも環系が異なるため）。
%[text]     イミプラミン（7 員環内に N を含む）とノルトリプチリン（7 員環内は C=C）
%[text]     はスキャフォールドが異なるため ECFP4 が大きく変わります。
%[text]  $&dollar&;T \\approx 0.26&dollar&;$ では同定の証拠として不十分 — $&dollar&;T = 1.0&dollar&;$ が必要な理由はここにあります。
%[text] Q2: イミプラミンとカフェインのタニモトスコアは？
%[text]       mol\_caf = `emk.mol.fromSmiles("CN1C=NC2=C1C(=O)`N(C(=O)N2C)C");
%[text]       fp\_caf  = `emk.fingerprint.morgan(mol_caf)`;
%[text]  `emk.similarity.tanimoto(fp_unk, fp_caf)`
%[text]     期待値: $&dollar&;T \\ll 0.15&dollar&;$（全く異なる構造）。
%[text] Q3: イミプラミンは 7 員環に縮合した 2 つのベンゼン環（ジベンズアゼピン骨格）を持つ。
%[text]     描画して数えてみましょう:
%[text]       2D 構造で識別できる環系はいくつか？
%[text]       7 員環の中の窒素はどこにあるか？
%[text]     NORTRIPTYLINE と比較 -- 構造の違いは何か？
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 7: MACCS キーによるクロスバリデーション
%[text] 単一の計算手法（ECFP4）だけに依存した鑑定結果は、裁判（法廷）の場で弁護側から「そのアルゴリズム特有の偏りではないか」と突っ込まれるリスク（脆弱性）を孕んでいます。そこで、全く異なるロジックを持つ「MACCS キー（フラグメントベース）」を用いて、結果のダブルチェック（相互検証）を行います。
%[text] 
%[text]   ECFP4（Morgan、radius 2、2048 ビット）: データ駆動型ハッシュ
%[text]     各原子周辺の局所環境を 2 結合半径まで捉える。
%[text]     小さな構造変化に敏感。
%[text] 
%[text]   MACCS キー: ルールベース・166 個の専門家定義ビット
%[text]     各ビット = 特定の構造特徴の有無。
%[text]     ビット 125: 芳香族環の存在
%[text]     ビット 160: アミン（環外 N）
%[text]     解釈しやすく、軽微な置換基の変化には敏感でない。
%[text] 
%[text]   クロスバリデーションプロトコル:
%[text]     ECFP4 $&dollar&;T = 1.0&dollar&;$ かつ MACCS $&dollar&;T = 1.0&dollar&;$  --\>  最大信頼度
%[text]     ECFP4 $&dollar&;T = 1.0&dollar&;$ だが MACCS $&dollar&;T \< 1.0&dollar&;$  --\>  さらに調査が必要
%[text]     ECFP4 $&dollar&;T \< 1.0&dollar&;$ かつ MACCS $&dollar&;T = 1.0&dollar&;$  --\>  構造アナログが見つかった
logInfo("MACCS キーによる同定のクロスバリデーション..."); %[output:7dd6d15c]
%[text] MACCS フィンガープリントを計算
fp_unk_maccs   = emk.fingerprint.maccs(mol_unk);
fp_top_maccs   = emk.fingerprint.maccs(mol_reveal);
%[text] 両フィンガープリント型からのタニモトスコア
t_ecfp4 = emk.similarity.tanimoto(fp_unk, emk.fingerprint.morgan(mol_reveal));
t_maccs  = emk.similarity.tanimoto(fp_unk_maccs, fp_top_maccs);

logInfo("フィンガープリントクロスバリデーション:"); %[output:50ab920b]
logInfo("  ECFP4（Morgan r=2、2048 ビット）: T = %.6f", t_ecfp4); %[output:5e04b032]
logInfo("  MACCS（166 専門家ルール）        : T = %.6f", t_maccs); %[output:793db436]
if t_ecfp4 >= 1.0 - 1e-6 && t_maccs >= 1.0 - 1e-6 %[output:group:672becea]
    logInfo("クロスバリデーション: ECFP4 と MACCS の両方が同一フィンガープリントを確認。"); %[output:58a8690b]
    logInfo("  同定信頼度: 最大。"); %[output:84682ba5]
elseif t_ecfp4 >= 1.0 - 1e-6
    logInfo("クロスバリデーション: ECFP4 は確認; MACCS は軽微な差異を示す。");
    logInfo("  フィンガープリント衝突の可能性は低いが -- MS/MS データを確認すること。");
else
    logInfo("クロスバリデーション: フィンガープリントが異なる。同定を再検討すること。");
end %[output:group:672becea]
%[text] DB 全体に対して MACCS ベースの同定検索を実行
logInfo("MACCS 同定検索を %d 化合物に対して実行中...", nRefValid); %[output:4e02eb69]

maccs_fps = cell(1, nRefValid);
for i = 1:nRefValid %[output:group:260e67a4]
    mol_i      = emk.mol.fromSmiles(ref_smi_valid(i));
    maccs_fps{i} = emk.fingerprint.maccs(mol_i);
    logProgress(i, nRefValid, "MACCS"); %[output:76ddb101]
end %[output:group:260e67a4]
maccs_result = emk.similarity.rankBy(fp_unk_maccs, maccs_fps);

MACCS_TOP_IDX  = maccs_result.Indices(1);
MACCS_TOP_NAME = ref_names_valid(MACCS_TOP_IDX);
MACCS_TOP_T    = maccs_result.Scores(1);

logInfo("MACCS 上位ヒット: %s  (T = %.4f)", MACCS_TOP_NAME, MACCS_TOP_T); %[output:933d7d44]
%[text] 一致チェック
if strcmp(MACCS_TOP_NAME, COMPOUND_NAME) %[output:group:21e0de0f]
    logInfo("一致: ECFP4 と MACCS の両方が同一化合物を同定。"); %[output:911e341b]
    logInfo("  最終結論: %s -- 確定。", COMPOUND_NAME); %[output:812a31b9]
else
    logInfo("不一致: ECFP4 上位=%s、MACCS 上位=%s。", COMPOUND_NAME, MACCS_TOP_NAME);
    logInfo("  追加の分析データで詳しく調査すること。");
end %[output:group:21e0de0f]
%[text] **✏️ やってみよう 7**
%[text] Q1: ECFP4 の代わりに MACCS キーフィンガープリントで同定を実行してみましょう。
%[text]     MACCS は 166 種の専門家定義ルール（芳香族環・アミン存在など）に基づくビット列です。
%[text]       fp\_unk\_maccs = `emk.fingerprint.maccs(mol_unk)`;
%[text]       maccs\_fps    = cellfun(@(smi) `emk.fingerprint.maccs`( ...
%[text]  `emk.mol.fromSmiles(smi)`), ...
%[text]                          cellstr(ref\_smi\_valid), "UniformOutput", false);
%[text]       res\_maccs    = `emk.similarity.rankBy(fp_unk_maccs, maccs_fps)`;
%[text]       ref\_names\_valid(res\_maccs.Indices(1))
%[text]     期待値: 同じ上位化合物（フィンガープリントタイプ間で一貫した同定）。
%[text] Q2: 別の未知物質 -- FLUCONAZOLE（DB にも収録）で試してみましょう。
%[text]     この新しい「証拠サンプル」に対して完全な同定パイプラインを実行:
%[text]       SMILES2   = "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F";  % フルコナゾール
%[text]       mol2      = `emk.mol.fromSmiles(SMILES2)`;
%[text]       fp2       = `emk.fingerprint.morgan(mol2)`;
%[text]       res2      = `emk.similarity.rankBy(fp2, ref_fps_valid)`;
%[text]       ref\_names\_valid(res2.Indices(1))
%[text]     Q: 上位ヒットのタニモトスコアは？
%[text]     Q: $&dollar&;T = 1.0&dollar&;$ か？（FLUCONAZOLE は DB にあるのでそうなるはず）
%[text]     Q: 2 番目以降のヒットはフルコナゾールと何を共有しているか？
%[text]        （ヒント: トリアゾール環に注目）
%[text] Q3: 未知化合物が DB にない場合はどうなるか？
%[text]     仮想的な新規化合物で試してみよう:
%[text]       SMILES3 = "CN(C)c1ccc(C2=CC=Cc3ccc(CN(C)C)cc32)cc1";  % 仮想化合物
%[text]       mol3    = `emk.mol.fromSmiles(SMILES3)`;
%[text]       fp3     = `emk.fingerprint.morgan(mol3)`;
%[text]       res3    = `emk.similarity.rankBy(fp3, ref_fps_valid)`;
%[text]       res3.Scores(1)    % 上位タニモトスコア
%[text]     期待値: 上位スコア $&dollar&;T \< 1.0&dollar&;$ -- 化合物は DB にない。
%[text]     上位ヒットは構造アナログであって同一化合物ではない。
%[text]     これが「データベースカバレッジのギャップ」問題である:
%[text]     法科学 DB は新規物質に対応するために定期的に更新が必要。
% ... （ここにコードを書いてみましょう）

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:7b3391ce]
%   data: {"dataType":"text","outputData":{"text":"[15:44:12][INFO]  initPython: Python 3.10 already active (Status: Loaded) -- skipping.\n","truncated":false}}
%---
%[output:264452a9]
%   data: {"dataType":"text","outputData":{"text":"[15:44:12][INFO]  S05: セットアップ完了\n","truncated":false}}
%---
%[output:0c6d22c7]
%   data: {"dataType":"text","outputData":{"text":"[15:44:12][INFO]  未知化合物 -- 証拠品 FC-2026-0004\n","truncated":false}}
%---
%[output:7292866c]
%   data: {"dataType":"text","outputData":{"text":"[15:44:12][INFO]    分子量          : 280.4 Da\n","truncated":false}}
%---
%[output:82477ed1]
%   data: {"dataType":"text","outputData":{"text":"[15:44:12][INFO]    LogP            : 3.88      （親油性 -- オクタノール\/水）\n","truncated":false}}
%---
%[output:96590970]
%   data: {"dataType":"text","outputData":{"text":"[15:44:12][INFO]    TPSA            : 6.5 A^2  （極性表面積）\n","truncated":false}}
%---
%[output:75df8e13]
%   data: {"dataType":"text","outputData":{"text":"[15:44:12][INFO]    HBD             : 0        （水素結合供与体数）\n","truncated":false}}
%---
%[output:060dc93a]
%   data: {"dataType":"text","outputData":{"text":"[15:44:12][INFO]    HBA             : 2        （水素結合受容体数）\n","truncated":false}}
%---
%[output:5b86f07f]
%   data: {"dataType":"text","outputData":{"text":"[15:44:12][INFO]    回転可能結合数  : 4\n","truncated":false}}
%---
%[output:91ba0e76]
%   data: {"dataType":"text","outputData":{"text":"[15:44:12][INFO]    環数            : 3\n","truncated":false}}
%---
%[output:84d44835]
%   data: {"dataType":"text","outputData":{"text":"[15:44:12][INFO]    ECFP4 ON ビット数: 29 \/ 2048 (密度 1.4%)\n","truncated":false}}
%---
%[output:87b69df0]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAcYAAAF\/CAYAAADTp2mUAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQnYTVX7xp9CyhBFZI5INJGPzDSQqaiUfIYvkTlTGYrM85CZRKRIKBnLPKTIkKFSGStDomRISsT\/f6\/TfjvOe8579hn2Pnvtfa\/rcuE9e6\/1rN+z33PvNT3PVZcvX74sLCRAAiRAAiRAAorAVRRGPgkkQAIkQAIk8C8BCiOfBhIgARIgARLwI0Bh5ONAAiRAAiRAAhRGPgMkQAIkQAIkEJwAR4x8MlIkUKBAATlw4IApSpFcG6rCwDpiqdO410wdKV1r5v5g\/Yn2PlOweVFcCdBXccWpfWUURu1dGL8O4MvBbAkUy1D3GteF+9xoN17CGGk9FEaznnfndRRGd\/o12l5RGKMl55H7zH5hhBudBavHjHil1H6ozyL9OVwZTLgh6uFeFlIaTYfqc6hHJ9zLhtmRe2D9Rh9C3e\/fx2DXxHp\/uF8Vs+2jnpTsC\/W5f\/vRPBvh7Ofn7iNAYXSfT+PaIzPC6H9NJF880QpjqC\/qcF\/g\/gLo\/wWbkv3BBD\/cS0Co0a\/RvhmBM\/MiYcbR4XwTzgex3h\/OxkjaD8Yv3P2BopiSuJrxS7j+8HN3EKAwusOPlvXCjDBG+0Ye+KUb7EvL7DVmRguhviRDffmH+3k4NtGKW0r1hmsznC\/CCZ2\/+ISzP9zn4R7KcC9RsX4ejEW4OsPZzM+9QYDC6A0\/p9hL\/6msSHEEm\/4LNd2V0rXhRoFmhc9MX6KZ\/gw2wgwnUtEIRzR1hvJZuPbDiUS095t9hqJpP5xwhxtVhmvTrO28zt0EKIzu9m\/MvQv3RW18EaXUkLFWF0oYjTbCfRHH3BkTFYSzwcyIy2gmXF3BzDHD20Q3Ql5ixn6zI1Z\/vxkNRjIdGU6kYv08mB\/C1RkLW97rHgIURvf4Mu498R99mf3Ci+SLJ1AQoxGSeHc6sM+BNpkRlnDCGMxm\/927kbAO1X8zo\/ZIfBVqJBY4kg\/kk5J94dqP9fNwo0d\/26x+IYn3c8r6rCVAYbSWr7a1mxGtSEY8ZkTPzDXhgPoLm1nRCLwnlEilJJJmWYT7Ag73ebj+h\/o8ErZmrzV7XSQ2+YtZrMIYiX1WcY\/WX7wvsQQojInl78jWQ42KzHx5GCJjZu3RjNCYaTOSN\/9w9aX0uRl7w9kSS\/vBRkBmHqBwAhNYh9l+RiI8Zl8c4i2MkbwcmR2pm2HOa\/QmQGHU239xtT6YqAX7kkSj4abqzHy5RnNNooQnGhGI5p5w4hdOWMOJXDh+ge1H8qIQznazthlthhP0cJ9HI8Zx\/YViZdoSoDBq67r4GR5qlJfSF51ZEU1p7UwXYYynwJkVtmjbNCN8KV0TarYg1NNmxocpPanh7o\/180jFOH6\/VaxJZwIURp29Z4Pt8foiN4Q00ORQB+2N6+I9YgnXn3BfxGbsSumacO0Hipb\/\/yOd6jPD3Hj5MdoJ9EewRyylCD2R2phS+\/4cg9ln9nMzLwuR+MWGXzs2kWACFMYEO8DpzZv5wjBzjZl+pjQ1Fur+cGHUzHyxm\/niNPMFHmx0EouwmWHGa6InEPjiEI2oR98673QyAQqjk71D20iABEiABGwnQGG0HTkbJAESIAEScDIBCqOTvUPbSIAESIAEbCdAYbQdORskARIgARJwMgEKo5O9Q9tIgARIgARsJ0BhtB05GyQBEiABEnAyAQqjk71D20iABEiABGwnQGG0HTkbJAESIAEScDIBCqOTvUPbSIAESIAEbCdAYbQdORskARIgARJwMgEKo5O9Q9tIgARIgARsJ0BhtB05GyQBEiABEnAyAQqjk71D20iABEiABGwnQGG0HTkbJAESIAEScDIBCqOTvUPbSIAESIAEbCdAYbQdORskARIgARJwMgEKo5O9Q9tIgARIgARsJ0BhtB05GyQBEiABEnAyAQqjk71D20iABEiABGwnQGG0HTkbJAESIAEScDIBCqOTvUPbSIAESIAEbCdAYbQdORskARIgARJwMgEKo5O9Q9tIgARIgARsJ0BhtB05GyQBEiABEnAyAQqjk71D20iABEiABGwnQGG0HTkbJAESIAEScDIBCqOTvUPbSIAESIAEbCdAYbQdORskARIgARJwMgEKo5O9Q9tIgARIgARsJ0BhtB05GyQBEiABEnAyAQqjk71D20iABEiABGwnQGG0HTkbJAESIAEScDIBCqOTvUPbSIAESIAEbCdAYbQdORskARIgARJwMgEKo5O9Q9tIgARIgARsJ0BhtB05GyQBEiABEnAyAQqjk71D20iABEiABGwnQGG0HTkbJAESIAEScDIBCqOTvUPbSIAESIAEbCdAYbQdORskARIgARJwMgEKo5O9Q9tIgARIgARsJ0BhtB05GyQBEiABEnAyAQqjk71D20iABEiABGwnQGG0HTkbJAESIAEScDIBCqOTvUPbSIAESIAEbCdAYbQdORskARIgARJwMgEKo5O9Q9tIgARIgARsJ0BhtB05GyQBEiABEnAyAQqjk71D20iABEiABGwnQGG0HTkbJAESIAEScDIBCqOTvUPbSIAESIAEbCdAYbQdORskARIgARJwMgEKo5O9Q9tIgARIgARsJ0BhtB05GyQBEiABEnAyAQqjk71D20iABEiABGwnQGG0HTkbJAESIAEScDIBCqOTvUPbSIAESIAEbCdAYbQdORskARIgARJwMgEKo5O9Q9tIgARIgARsJ0BhtB05G3QTgaNHj8rGjRulQIECUqRIEUmbNq2buse+kIAnCVAYPel2djoeBCZNmiTPP\/+8XLhwIam6a6+9Vm644YYr\/uTMmVNy5MgR9Oe5cuWimMbDGayDBOJIgMIYR5isyjsETp06JdmyZVOimCpVKkmTJo1cunRJ\/vrrr4ghZMyYUbJkySI33XSTZM2aVf0bf4x\/B\/v5NddcE3E7vIEESMAcAQqjOU68igSuINC6dWuZOHGiQNR+\/PFHyZAhg\/r8t99+k59\/\/ll++eUX9efEiRNJf\/B\/fIaf+X8Wi5hCnA0h9RfTq666Sv78808pVqyYVK5cmd4jARKIgACFMQJYvJQEQGDlypVStWpVwaht8+bNcvfdd8cE5o8\/\/pCTJ08m\/cG6JcTW\/2fGv\/HZ4cOHIxqZPvfcc\/L666\/HZCNvJgEvEaAwesnb7GvMBE6fPi133XWXHDp0SAYPHixdu3aNuc5oKjhz5owadRojUGNkaoxEt23bJl9\/\/bWcPXtWMHqEoGbPnj2apngPCXiOAIXRcy5nh2Mh0LBhQ5k5c6aUKVNG1q9fr9YXnVywEQjroX379pVXXnnFyabSNhJwDAEKo2NcQUOcTmD+\/Pny2GOPSbp06WTHjh1SqFChJJO\/++47OXbsmJQuXdpR3ViyZInUqlVLbr75Zjlw4IBcd911jrKPxpCAEwlQGJ3oFdrkOAKYsrzzzjvl+PHjMmHCBGnVqlWSjdiNev\/998uGDRvUaPKpp55ylP3\/+c9\/5PPPP09mt6OMpDEk4CACFEYHOYOmOJfAE088IfPmzZOHHnpIli9frtbtjDJ8+HDp3LmzOqv45Zdfql2iTirvvfeePPnkk5I3b17Zt2+fOlrCQgIkEJoAhZFPBwmEITBt2jR59tlnJVOmTEr48uTJk3THN998IyVKlBDsLF28eLHUrFkzYTzPnz8vsHXTpk3qb6NgRIuds7t27ZLp06dL48aNE2YjGyYBHQhQGHXwEm1MGAEcjYCo4LjEjBkzpEGDBkm2XLx4UcqWLStbtmyR5s2bCyLhJLJgp2q+fPnUZptPPvlEypUrl2QOBPGZZ56R22+\/XQnk1VdfnUhT2TYJOJoAhdHR7qFxiSSAkRamTtesWSN16tSRDz744ApzevXqpXZ75s+fX3bu3KkO+ye6dO\/eXQYOHCiPPPKILFy4MMkcROgpXLiwYJMQplYxNcxCAiQQnACFkU8GCYQgMHr0aOnQoYMK1fbVV1+pEHBGwTlB7ED9+++\/ZfXq1VKpUiVHcMR5xltuuUV+\/\/132bp1q9x7771JdmHTUJs2bVQ0HNjvv07qCONpBAk4hACF0SGOoBnOIvDtt98qUcHa4fvvvy+PP\/54koFYy8O6IqYkX3zxRRk2bJijjO\/YsaOMGjVK6tWrJ+++++4VdiMLCKLqfPjhh1K9enVH2U1jSMApBCiMTvEE7XAMAawdli9fXm1iadKkiUydOvUK2wzhQZopHINw2tlArIsWLFhQ0A9Ev7ntttuS7IeId+nSRQUowPESFhIggeQEKIx8KkgggEC\/fv2kZ8+ekjt3bvniiy9UuiijYFMLpk2xeQXCUrJkSUfyw2agyZMnq920b7zxRpKNmGLFVCtCx61bt04qVqzoSPtpFAkkkgCFMZH02bbjCGzfvl2tHWK0hWDhOLhvFMQdxfrc\/v37pU+fPko8nVoQ5QabbbCOuHfvXrVb1SiwvXfv3vLwww\/L0qVLndoF2kUCCSNAYUwYejbsNAJYO8QIEGcV27dvr9bp\/EvTpk3VtCrWHj\/77DPHH5TH0ZJ33nlHJVMeM2ZMUld+\/fVXNWpEiixkB3HqqNdpzwft8Q4BCqN3fM2ehiGA6DWIYoOzfti16b92uGjRInn00Uclbdq0al3xjjvucDxPBB9AGDukx8IIEpF5jIJ1Rqw3YlMRNhexkAAJ\/EuAwsingQRE5NNPP1Vrh5h6xDrifffdl8QF63FINfXTTz\/JyJEj1REOXQrOXy5YsEC6desmgwYNSjIbAc9x\/hLJjLGOCgFlIQES8BGgMPJJ8DwBbEjB2iHiiCI1Ew7t+xcce5gzZ47aqYoNKzpFjcHIF0HEEXzg+++\/v2IjUdu2bWX8+PGCVFpvv\/22558DAiABgwCFkc+C5wm0aNFCZbgvXry4WjvE1KNRkC0DwpEhQwaVaurWW2\/VjlfVqlVlxYoVyXIyItkyjnUgSAHObeLfLCRAAhwx8hnwOAEIBnZnQgwR8xRTpkbBQXj8H5tVcOQBRx90LGvXrlW7a2+88Ub54YcflMgbBec033zzTcHLwWuvvaZj92gzCcSdAEeMcUfKCnUhgGDbED4ciMdGFESxMcrly5dVpoyPPvooWdxRXfrnb2eFChXU2umIESOkU6dOSR9h+hibjVKnTq2OoeTKlUvH7tFmEogrAQpjXHGyMp0I1K9fX4VMQxYKrB2mSpUqyfyJEydK69atJWvWrOr4xs0336xT15LZumTJEqlVq5bqB3ao+u+4RWLluXPnKsGEcLKQgNcJUBi9\/gR4tP\/IlIGjCunTp1drh\/7ra8hAcc8996hzfrNnzxYIhxsKNuHgqAlEv2XLlkldQmYQrK+mS5dOZd9A0HQWEvAyAQqjl73v0b4fP35cTaHib+RQRPg0oyDVFNbjPv74Y5V7ETkY3VIwKoTI45jGnj171PSpUTCaxKiyR48egpB4LCTgZQIURi9736N9N0SgSpUqsmzZsivSLw0ZMkSd+cuZM6eaQsWGFbcUiD5eCBBY\/K233pJGjRoldQ27cRFYPFOmTOpYR+bMmd3SbfaDBCImQGGMGBlv0JnAlClT5LnnnlNf\/DjYnidPnqTuQDCQTgqh4TB6cmNaJuxAxU5UZAZBjkn\/M5kYKWMHKwIB4OWAhQS8SoDC6FXPe7DfGAlh7fDMmTMya9Ysefrpp5MoIGg4RkxI7tuqVStBUl83lgsXLqg0VGARmGcSR1dw5hEJmbHWiDVHFhLwIgEKoxe97sE+YxrxwQcfVCOixx57TObNm3cFBaytDRgwQK2\/YTMKIsW4tSDaDaLeINoPIuMgDJ5RypYtKxs3bpTRo0dLu3bt3IqA\/SKBFAlQGPmAeIKAEd0GoyGsHeJvo2CnJkaLiACzZs0a1+coRHxURPBBAAOknUKAA6PMnz9fvTggFyXONfpHAfLEg8JOkgBjpfIZ8AIB4wwf+opzi4h9apRz586pNFK7d++Wrl27yuDBg72ARIYOHar6ixcCJFw2CgIbYLoZLw9Yj0WqLRYS8BoBjhi95nEP9rdx48YqSDYOtx89evQKAsZBfmSXwPoi0kp5oSBwOnIyInMIjqYgMo5RkMMRR1UwqkQMVf9jHV5gwz6SAIWRz4DrCeBYAs4jItrL2bNnr9iJiRESAoiXKlVKHXL3Uundu7f06dNHqlWrpkLfGQVTyti1unfv3mSblLzEh331LgEKo3d975meYxcqRot\/\/PGHOoZRo0YNz\/Q9pY4iODpGjYjwgwDqiIxjlMmTJ6vAB0WLFlXTqk5NtXXkyBG1Npw3b176lATiRoDCGDeUrMjJBBADFEHCS5curXZdsvgIdO7cWYYPHy5PPPGEvPfee0lYcKwDYfIOHjwoCxcuVIHUrSzYEAShPnnyZIp\/MBWOTUO4DpGLIIrYVfvAAw\/IypUrrTSRdXuIAIXRQ872clexpoajGD\/\/\/LPaeVq5cmUv40jq+7FjxxQXCBMCHmCt1SijRo2Sjh07qmnmTZs2meKF86AnTpxQIhf4t\/EzrGsGfo7RfDQFoojpcJT\/\/e9\/KoUWCwnESoDCGCtB3q8NAcQA7dmzpyAU3PLly7Wx22pD27RpowIaBAY2wI5diCZGZgMHDlRRgvwFzf\/fEDsIIaatoynXXnutZMmSRYXg8\/8b2U0Cf4b\/Gz9LkyaNyoJi5JLEFDB30kbjAd7jT4DCyOfBMwROnz4t+fLlE\/yNIwo4qsAiKgoO1l6RiNk\/HRXYQGQwxRqJ4EHkbrjhBhVvNkeOHOrf4f7g2lgKhB0Cj9RhOJJTt27dWKrjvR4nQGH0+APgte4jBigChdepU0eQeoolNAEkcsbUKja44Gwj\/h1s9IYRnvEHQcgTVXr16iV9+\/ZVQQkWLVqkwtuxkEA0BCiM0VDjPdoSwLQgdmJiTQ2h35BtgiU4AeOYC0bW69evvyKRs1OZIdnyyJEjVUg\/rCUjKDwLCURKgMIYKTFerz2B559\/XsaNG+e6fIvxdMyCBQvUqBqBxLdv364Cj+tQsBEHU8LYhIP1SQQvwJlMFhKIhACFMRJavNYVBA4dOqSOImCr\/zfffCOFChVyRb\/i1Qns3MVIGjtWEXAcm1t0KjhqgnivWDdFzNdPPvlErS2zkIBZAhRGs6R4nasIYFPJ1KlTVW5GRL5h+ZcANq4gJdVDDz2kdu\/6Z9\/QhROOfyA4OqaA8eKDv7Nnz66L+bQzwQQojAl2AJtPDAFkjihcuLCK6LJv3z5GTvnHDdOmTVNTkdhEg3ONOkeUwe5jJF\/GVHDJkiVl1apVrk4nlpjfJHe2SmF0p1\/ZKxMEkKh49uzZ0qFDB7Vhw+vl8OHDcvfdd6uoMgi63rBhQ+2RYLNVxYoVVfYUiOSHH34oOE7CQgIpEaAw8vnwLAHsSkXgcJzdw1m+m266ybMssGmlevXqsmzZMqldu7YgL6NbyoEDB6R8+fIqswr6hnOZzBjiFu9a0w8KozVcWasmBBADdPHixdK9e3fp37+\/JlbH38wxY8ZI+\/bt1csBgoa7bT3uq6++kkqVKqnIPUhDhl2rOq6dxt\/zrDEYAQojnwtPE0AMUAQWv\/766+WHH36QzJkze44HphmRrBkh4DCaQkBxNxb4GhuKkHqsXbt2Mnr0aDd2k32KAwEKYxwgsgq9CSAzAw6DDxgwQF5++WW9OxOh9Qj6jWlGiEaTJk3UTl03F2TgqFWrlpw\/f15FQOrSpYubu8u+RUmAwhglON7mHgLYrYiRBMKaYa0xQ4YM7ulcmJ5g+viVV15R5\/2wCxUxTd1e5s2bJ0899ZRcunRJJk2apI7ssJCAPwEKI58HEhCRcuXKqcDiSLWEtTYvlB07dsh9990nOBD\/0UcfqXN\/XinIxoFsIgg6PmvWLHnyySe90nX20wQBCqMJSLzE\/QSQjBc7FjFywhlHBKJ2c8FUIs72YaONV9fbjDRk8DX876UXAzc\/2\/HoG4UxHhRZh\/YEcFwBGSQgFMjp16xZM+37lFIHsLY2bNgwFeRg27ZtKiaqF8sLL7wgr776qur\/ihUrpGzZsl7EwD4HEKAw8pEggX8IYErtv\/\/9rxQoUEAdCHfrWTdMGePQO44rII4oplO9WvBChJcgbDrCGjOCjhctWtSrONjvfwhQGPkokMA\/BBBUHF+Ke\/bskXfeeUfq16\/vOja\/\/\/67Cmqwd+9e6dGjh2A60esFfsdmHGzKyZUrl3pZQGoyFu8SoDB61\/fseRACU6ZMUbsUIZCYVkUsVTeVli1bqp2YxYoVU0c03L6WatZ3CDperVo1NWJE5hWIo9uCHJhlwetEKIx8CkjAjwB2aCIbAw77IywaNuS4pWANDRtMIIZbtmxhkuYAx545c0bFU8WaK2LGrlu3zpMBH9zyvMfSDwpjLPR4rysJICIKAouXKlVKjarcUE6dOqW+7JGLcujQodK5c2c3dCvufUAuygoVKqg15sqVK6tjLAw6HnfMjq+Qwuh4F9FAuwn8+eefkj9\/fvnpp5\/UTkUc\/te9YFMRNhdh1yWmC3F+jyU4Abw8IBrQwYMHBbF0sfbo1o1YfAaCE6Aw8skggSAEBg0apMLDYWpt9erVWjP64IMP5PHHH5f06dOr3ISYKmZJmcCuXbtU0PETJ06o9FvTp0933Xozn4HQBCiMfDpIIAgBrDfly5dPMAWJjRiIjBNYsGEDIy8nb2BBPsK77rpL8DeivbRo0YL+Nklg8+bN8uCDD6qg423btpWxY8eavJOX6U6Awqi7B2m\/ZQRwnAGBxRF0etGiRcna6dixowohhzUoxBj1\/5MzZ07JkSNH0J\/jSEDatGkts9u\/YowUMWKsUqWKyrXIVEuRYcdsQY0aNVTQcS8GmY+MlnuupjC6x5fsSZwJYBoN59kwYvj8889Vaib\/0qZNGxUlBztZIy0ZM2aUrFmzqvyHOFhu\/MHP8G\/83Pg3\/o9\/RzoyfeONN9ThdaTSQoDwPHnyRGomrxdRu5Pr1q0rOO+IF43333+fXFxOgMLocgeze7ER6NSpk4wcOVLq1asn7777btDKMO36yy+\/CHY0QkyNP\/gZ\/o1pzMCfx0tMYR8O7AcWZAlBiDvY5tZgBbF5NrK7EVgeyZxRFixYII8++mhkFfBqrQhQGLVyF421m8DRo0dViDgIGTZkILZoPArWJ0+ePJn0B+38+OOPV\/zM+ByfHT58WP76669kTSOPJI4V+BekU8La2Nq1a6VOnTpqKpUldgKYAoePkMgZCZ1Z3EuAwuhe37JncSKADSuvv\/662sK\/fv36ONUaeTUY\/WFUaoxEMQqtXr26mmb1LxjhYiSZLVs2Fb0Hf7PETgBT0XhBwRT6uHHjYq+QNTiWAIXRsa6hYU4hgPXF\/\/znP8qcxo0bq9EYQqoZa392baQxw+Pbb79Va6EYkeL83WOPPWbmNl5jggB293711Vfq+A6O8bC4lwCF0b2+Zc\/iSADxM5GnMVhxyq7UixcvqmMlOGbQtGlTQdxXlvgQwPR0hgwZBMEfTp8+Ldg8xeJeAhRG9\/qWPYsjAXwZNmnSRMUYzZQpkzrsjalMTGsGW\/sL1zS+ZAN3nwbuQjV2rBo\/Dzcy7dOnj\/Tu3VvtpN25c6dcf\/314czg5yYJfPfdd2qtGeuMmE5lcTcBCqO7\/cve2UDgt99+C7n7FOIZuC4Yq5hizdCYxjWOeZw7d06GDx8uGNmsWrUq2YYcGzC4uomlS5eq9dwHHnhA8WVxNwEKo7v9y945lEDgrlTsQA21MxU\/P3LkiDpknlLBCLNq1aoyY8YMh\/ZaX7OMwPKtWrWSCRMm6NsRWm6KAIXRFCZeRAKJJ4CRqf95Sf\/dqUuWLFFxUNu1ayf4EmeJL4HWrVvLxIkTVaQjnGlkcTcBCqO7\/cveeYQA4rkiXdKtt94q+\/bt80iv7esmdiJjNyrSUCGhMYu7CVAY3e1f9s4jBLC2iPisx44dU4EIihYt6pGe29PN3Llzq+nsAwcOqJRkLO4mQGF0t3\/ZOw8RePbZZ2XatGkycOBAeemllzzUc2u7+vvvv6vjGYhVi38zl6W1vJ1QO4XRCV6gDSQQBwJG3sUyZcrIhg0b4lAjqwABI8ADDvgjGDuL+wlQGN3vY\/bQIwRwZAM7U3EIHdN+N998s0d6bm03EYS9QYMGKsPG3LlzrW2MtTuCAIXREW6gESQQHwKPPPKILF68WKXDQsopltgJ9OrVS\/r27Svdu3eX\/v37x14ha3A8AQqj411EA0nAPAEEO0fQc6RFQnokltgJPP300zJ79mx56623pFGjRrFXyBocT4DC6HgX0UASME8Au1KxOxXh4xBxJ3369OZv5pVBCSDf5Y4dO2TTpk1SqlQpUvIAAQqjB5zMLnqLQOnSpdWXODLP165d21udj3NvL1++rGLOnj17VuXKzJw5c5xbYHVOJEBhdKJXaBMJxEBgwIAB0qNHD2bYiIGhceuhQ4ckb968kj17dvnpp5\/iUCOr0IEAhVEHL9FGEoiAAHIG4mgBgo0j4zzP3UUAL+DSlStXSpUqVaRixYqybt266CvinVoRoDBq5S4aSwLmCBQqVEiFhvv000+lbNmy5m7iVckIjB8\/Xtq2bSvNmzeXSZMmkZBHCFAYPeJodtNbBDp27KgCXnfr1k0GDRrkrc7HsbcIyj527FgZMWKEdOrUKY41syonE6AwOtk7tI0EoiSwZs0alTvwzjvvli+\/3BllLbwNabxWrFihzobWrFmTQDxCgMLoEUezm94icPHiRalVa69s2VJYPvvsailUyFv9j1dv8+XLJwcPHpS9e\/dKwYIF41Ut63E4AQqjwx1E80ggWgING4rMnCkyfLjICy9EW4t370My6QwZMkjq1KlV8HD8zeINAhRGb\/iZvfQggTlzROrVE6lYUYQbKiN\/AHbu3CnFihWTIkWKyNdffx15BbxDWwIURm1dR8NJIGUCZ8+KZM0qcvGiCI7g4d8s5gnMmTNH6tWrJ3Xq1BFkLmHxDgEKo3d8zZ56kMDDD4ssXy4yfbpI48by7UDnAAAgAElEQVQeBBBDl\/v16yc9e\/aUrl27yuDBg2OoibfqRoDCqJvHaC8JREBg\/HiRtm1FnnhC5L33IriRl0rDhg1l5syZMnXqVGnSpAmJeIgAhdFDzmZXvUfgyBGRPHlEEEv8559Frr3Wewyi7XHJkiVl69atDJIQLUCN76Mwauw8mk4CZgiUKCGybZvIkiUiNWqYuYPXgECmTJnkzJkz8ssvv0iWLFkIxUMEKIwecja76k0CvXuL9Okj0rKlyMSJ3mQQaa8RYzZXrlySNWtWlb6LxVsEKIze8jd760EC27eL3HuvSM6cIocPi1x1lQchRNhlI3JQuXLl5JNPPonwbl6uOwEKo+4epP0kYIJA\/vwi338vsnmzSMmSJm7w+CWvvfaatGrVSp599ll54403PE7De92nMHrP5+yxBwlgZyp2qPboIdKvnwcBRNhlIwj7kCFDpEuXLhHezct1J0Bh1N2DtJ8ETBDAWUacabz7bpGdjCkelliNGjXko48+kvnz50vt2rXDXs8L3EWAwuguf7I3JBCUwIULItmyiZw6JXLggAimVllCE7j11lvlwIED8s0338jtt99OVB4jQGHU1OEbNmyQ4cOHC85avfTSS5r2gmbbSQBxUxE\/dfRokXbt7GxZr7bOnz8v6dOnl6uuukoFD7\/mmmv06gCtjZkAhTFmhImpYOnSpVK9enWVnR1Z2llIIBwBZNro3l0ES2atW4e72ruff\/XVV3LXXXfJbbfdJrt37\/YuCA\/3nMKoqfNx8PiGG26QNGnSyOnTpyVt2rSa9oRm20Xg779FUqWyqzV920EIuKZNm0qVKlVkORZnWTxHgMKoscvxVou3W0yrlilTRuOe0HSrCSA5xIsvijzwgMjkyVe2tmqVSPPmIosWiRQtarUlzq3\/woULMmHCBOnRo4c62I9D\/lim6Natm1zLWHrOdZwFllEYLYBqV5UtWrSQ119\/Xa01vsBMtHZh17Kdt94S+d\/\/fKYvXixSs+a\/3Vi4UAQbLxE2rnhxLbsXs9ELFy6UF198Ufbu3avqyp8\/v3z33XdJ\/8bv2OOPPx5zO6xADwIURj38FNTKN998U0X9f+KJJ+Q9pk7Q2JPWm24IY\/nyIkePinz5pch11\/na9bIwfvvtt+ql8sMPP1QsChcuLCNGjJCaNWvKunXrpH379oKExSj333+\/jBo1Su7GmRcWVxOgMGrs3j179qhf5OzZs8tPyETLQgIhCBjCuH69yEMP+aZV+\/f3rjD++uuv0qdPHzV1evHiRbVej7yLONjvvwv10qVLMmPGDOncubMcP35crr76amnQoIESz5tuuonPm0sJUBg1duzly5eVKCLIMc5cYfqHhQSCETCE8ccffcc1Ro4U2bFDpEgRb40YIYLYXNO9e3eVNSN16tQq7Fv\/\/v1TFLpTp06pZMUjR46Uv\/76Swlpr169pE2bNqoOFncRoDBq7s9HH31UFi1apBKq\/ve\/\/9W8NzTfKgL+wpgpk2+TTb58ImvX+jbdeGGNceXKlWpEiA1rKA8++KASOmxiM1swS4M6jKlXHP5HHdWqVTNbBa\/TgACFUQMnpWTioEGD5OWXX5a2bdvK2LFjNe8NzbeKgL8w5sghgiXpJ58UmTFDJGNGdwsjxAw7TefOnavwFipUSAYMGCBPAkCIcvbsWTl48KAUDbFNFyLbrl07FRkHpVatWjJ69GgpUKCAVS5kvTYSoDDaCNuKpiiMVlB1X50URgqj+55q63pEYbSOrS01Y+dc5cqV5d5775XPP\/\/cljbZiH4EAoURPcDs365dIsOGidSv7zuugXyNY8b4zjVielXnaGjGuiB2kiLMW4YMGdQOVJxNDBcQA7Mww4YNU+uPGF3iXGNgMc499uzZUxBwA5t2WrZsqdYrM2IYzqItAQqjtq7zGf7HH39I5syZ5e+\/\/xZ8EeCXn4UEAgkEE8Z9+0SwvIbTB8jTCGHs00dkwQLf3dmz+84+NmuG6Ud9mIbaSQqhw2Y1MwVnGjE1is06EMV+\/frJc889J6mChA46+v\/nX3r37i1TpkwRtJ0zZ061MadZs2ZqFyuLfgQojPr5LJnFpUqVki1btsiqVavkAYQ2YSGBAALBhBGX9Ool0rev72III5bIZs8WmTjRt2vVKCVK+EaRDRuKpEvnXLxr1qyRDh06yBdffKGMxGwKRoz33HNPxEYjTio22iD9FEqRIkXURpuHkb8rSMGMDc49GrGLEeAf4sqoVBGjT\/gNFMaEuyB2A\/DLOGbMGPVWi00GLCQAAufPi7z8ski3biL4bsfoD8c1sPnGKLgGI8Y9e5JHvsHM\/Ouvi7zzjsjZs747sKMVWToQhDwKrbHMMfv27VOb0IwNNnny5FFTmo0bN465Tez6htjiSBQKNtrg9y3Y8SgcoUKwDYw4sXkHGToaNmwoQ4cOlZtvvjlmW1iBPQQojPZwtrSVd999V+rXry9IrrpkyRJL22LlehCAkNWpI4I4qDjQP2GCyJo1wUd8GBliKhURz4Ispclvv4kg1urbb4usXPlv\/41RJE4JJWoGH2mhMEU6ZMgQ+fPPP1W6KIhSvOOb4uzixIkTxVhPvO6669SuVJyHDLaeeO7cOSWGgXaZWd\/U4wlzt5UURhf4F2+m+fLlU2uNJ06c4LqGC3waSxdOnvTFQt24UQSDFIwWixWLpcZ\/7\/36axFMy06ZInLihO\/n118v8vTTvqlWiKUdxVhHRLQaRH0yRmYQohz+Q+I4G4PA4oiYY6wn5sqVSwYOHCiNGjVSNgSWQ4cOKfF8G28VIlKwYEF1fUpHReJsMquLggCFMQpoTryldu2psnv3vfLBB\/dIkSLJf0GdaDNtij8BxEHFEhhiod5yiwiyJlmxcebcOREcC0SmDv90oA0aTJPy5c+rYBPXQzEtKJs2bVJTm5999pmqHWvsWEe0cy1v69ataj0RmW1QsJ6I6dXSpUsH7fHq1auVzV\/CMf8EF4DNd955pwWEWGWsBCiMsRJ0yP1Gdna8yTdt6hCjaIatBLAEVrWqyP79vsg2EMVcuaw3Abl8p00TefPNS5I6dQE5cuQHlabpkUcekebNm8tDmMuNQzl8+LBaR0TsUqzlhRutxaHJFKuADRgJBo5aQ+1+NcLRYR8Awjga4ehCHQex2n7WH5oAhdElT8eoUSIdO\/pEEeLI4i0C2FFavbrI8eMYQYkgWUSWLPYyOH\/+osybN1cmT54sa9euVeKFglERjjpguhExRiMtxnod1uxwPCldunTy\/PPPq41mTjieFLjOGe68ZGAA8xtvvFGtXSJ6VbDjIJHy4vWxE6Awxs7QETVg88R99\/mCQmMdiMU7BNat8x3GP30a8T9F5s9P3GYYgzp2iWJkh4DdWGdDwaF6xPbFKBJxSoOtyfl7LdgOz7p166rNNlhTd1oJ3BmL0HOvvvqq2sUarCDlFaZXly1bpj4uXry4mhKuWLGi07rmOXsojC5x+YULIpkz48C\/b9QQbHehS7rKbvgRQABwTKPD74heM326SJo0zkGEwBM4W4iE2h988IE6MI8C0WjatKnKJ5otW7ZkBuNcLkTDWMMrUaKEEo3ySCjp8ILzxLDdCFaOqWTYfscddwS1HMdBsF5pJEaGkCLu8S1YJGZJCAEKY0KwW9NopUoiH3+cPEO7Na2x1kQTQADwZ58VwUtRq1Yi48aJODnQypEjR9Qo8rXXXpPvv\/9e4UMYtapVq6rzho8\/\/rjaYdq3b1\/to8gEprdKkyaNtGrVSvUtEw6DBhTjOMgrr7wiv\/32m+A4CF4c8P9gLw6Jfvbc3j6F0UUexkHuIUNEunf\/Nwmti7rHrvgRQCKVDh1ELl0S6dpVZPBgffBgFIloMliLRPomYxSJA\/CnT59W64gQhk6dOqm4pjibqGsx1hPHjx+vwjZmyZJFiV2o9URsMOrSpYvgbDKmkrFBB9F2cD2LfQQojPaxtrwlxLjEoe777xdZvdry5rRqAG\/f8+bNU8Gkq1SpojY7YHpOx4KXH7wE4dgcAoC\/8IKOvfDZjDijb731lkyaNCnpPCLCGoaKLKNrT3fs2KGmSz\/GlI6ICvqP6dUKFSoE7RKmn7GrFxt7sNkIf7PYR4DCaB9ry1v65RcRLNcgluWpUyJMLO5DjgwJ03CeIKDgsHWdOnWkdu3a6gyc03cEYsTRsWN3Wbmyp+zbl06mTvVFsnFDMXabYmoVLy9uLcHWE8eNGxd0MxGmV40sIMYOX7dycVq\/KIxO80iM9tx2m8jevSJbt9oXhSRGky27HV8miFKCPyjYEYmD2Ng9uHz5cjVaMQqmuBBSD2\/p1atXd8QxAH8w+JLEcYc5c+ZI7tzFZdKkTVKjhoN22cToRYghzj5CCBDazc0FLwEYESOWKxIiY0TYuXNndR4SU8hGQXQf42WNwmjvE0FhtJe35a0984xvZyLWoLy8LIHRFTY7YB0L6zTYFYkdkP5fOtu3bxe8wUNsjEzs+BxfTjhOAJGEmCY6+DOm0XBMYenSpWrjBmwONQVn+QNmUQMQQ3CHOEI4\/P2EjTpI3+S2XZoI5Yj1RDx\/ED4c28BGpEBhRN\/xPLPYR4DCaB9rW1qaNEmkZUsRBHaeOdOWJh3XCEZXyGiATAt4G0e2A4wCUyrInADBwT0bN25UefVQ8MaOMF8QSUy7Fi5c2Nb+njx5Up2Dw7EF5BLEphWcd3NbCSWMSACMlwGEl8PGHDcWrDsuXLhQhg8ffkX3jBEjhdF+r1MY7WduaYsIxYg0QjgC9d13ljblyMoxusK2f0yVIqg6xC7Ss28I1wUBgkiuWLHiijWvokWLKpGEWJUrVy7sIfVYIOHoAnL\/IbcgDrSjT7dhrtyFxVhjxKgRkW6M4gVhDOVOjBIx24GXM2Pnrgtd78guURgd6ZbojcJA58YbfVFQDh+2J1Zm9NbG906MrmrWrKlGfJj+hLgVizGtBL6kcWAbIgmRPYVdTf8U5PzDSBQiCQHDxpF4FRz2xrQa1kORIBeimDt37nhV77h6QgkjRol4wcGo0Z+94zpggUEURgugmqySwmgSlE6XVasmgihTyH5Qt65OlkdvKzbSQJyQvQBrURASRFeJZ8EXFUQXIokoLkaoM7SBGKCIcAKRxC7XYIe4zdqCiCnoC1IcYbMQzvpldXkoI7yA4Lxi4NEELwsjRokIDIBR4wVEcWCxjQCF0TbU9jWETZi9e4t06iQyYoR97SaqJawPYnS1f\/9+wVQnRBGZF6wuu3btUiK5ePFi+Rzp7v8p2FmJzTEQSWyaicQWpFTC7lgcDMd5vvnz5wdNhGt13+yuP5QwYpSIlw6MGjEj4KVCYUyctymMiWNvWcv79vlSDyE1XJDoU5a1m4iKt23bpqYzjx8\/rvLyYXSFoxd2F+ycXLBggRLJdevWXfGGD7FGYlqsTaYUVGDJkiXqOkwrYqPPrFmz1C5NLxSsDSMrBUaNOMJgFC8LI0aJmJ7HqBEbyljsI0BhtI+15S0hy8\/Qob7MCq1b+yKj+JeRI32Z3BEZxw0FAoRpS0y34XgFRldOSEOE0R4EEn9wxAKxL41SoEABNZKESFauXFlNk6G888478swzzyhBxd\/GMRM3+MlMH0IJI0aJSMuEUSO4eqlQGBPnbQpj4tjHvWUcdTKi3bz9dvKoKEio3qKFL4yY7gUbYerVq6dGV\/Xr15fp06erN2unFUwRYmrXGE3+gvBE\/xQEh8aXPnacYvcrtufjkPegQYMs3e3qNEawB6PEjBkzqhcb\/xcJLwujEfnG7dGAnPg8Uhid6JUobTKEEemnsEHy22+xKeTfytwijMjQgDBveKPGIX6E1MJZL6cXCJ8RVGD27NmCfHxGQW7CESNGSEdkm\/ZgCSWMGCViahwvECdOnPAUGS9FA3KaYymMTvNIDPYYwvjSSyITJ4o89ZQIDvwbxQ3CiDx1yHVnjK4G65RWIsC3OE4CMcRZvWrVqqmURF4tGCXiED9GjeBhFAqjN8LkOe25pzA6zSMx2GMII6ZKMauIXanr14uULeurVHdhHDJkiHTr1k1NMyKL+ws6p5WIwc9uvDWUMGKUiKMqGDX6T0O7kUFgn0JFA\/JC3xPdRwpjoj0Qx\/b9hbF9e6S28SWuxUkCrD3qKow4P9imTRuVmghRQBD3FFOpLO4hECrCDYUxefxY93jduT2hMDrXNxFb5i+ML74ogtRvlSv\/m7NPR2HEBgRkd8eaHM4HYvcmQr6xuItAKGHEKPGmm25So0aE6vNSCRUNyEsMEtVXCmOiyFvQbqAwoonGjUXmzxfZs0cEYTaxKxUBAIoUEalUSeSRR0Rq1BBxYpJ0L2SVsOAx0LLKUBFuKIzpVNYR\/\/ixWjpYM6MpjJo5LCVzgwnjsWMit9\/u24gza5ZPGHGOsWbNf2vCuUeEkatd2yeSiLWa6OKVrBKJ5uyU9kMJI0aJONaCUSOCOHiphIoG5CUGieorhTFR5C1oN5gwopnx40U6dPA1iL+xOefAAZFFi3zxVDdsEEFwAJRUqXwRczCSfOwx3yjT7uKlrBJ2s3Vqe6Ei3FAYk8ePdaoP3WQXhdEF3sSZcexCxcgPm2wgfFhjNAoybpQpI7J5s+\/ngQf88SK+dKlPJJcvF\/GPPlW0qE8ka9USKVcueTSdeOPzWlaJePPTtb5QwohRIvJQYtR4DNMfHiqhogF5CEHCukphTBj6+DR89KjIww\/71giRUQM5bHGOsVmzK+uHKNav7\/s5Pg9Vfv9dZPVqn0guXOhLX2WUvHl9U64QSbQZxyxLqgkvZpWIz1Ogfy2hItxQGJPHj9Xf287vAYXR+T4KaeHu3SJVq4ocPOhLTrxihUi2bPHrEKZmN270ieS8eb78jkYpV6625MyZVsX9RMBrHM6Opaxfv17FD8Vak5eySsTCzE33hhJGjBKRWxOjRkyxe6mEigbkJQaJ6iuFMVHkY2x32zaR6tVFMA1aqpTIhx+KWJlUAmuQW7b4driuXfurbN6cTXC+EAUZIJCL8NFHH1V\/8CUWSfFyVolIOLn52lARbiiMyePHuvk5cErfKIxO8UQEdqxb59tBimnOBx\/0iRXWF+0sWAtcuHChyiCxdu1aQe44FMQsLV68uBpJPvXUUyo\/Ykpl5syZ0qRJE89mlbDTZ05uK5QwYpSYI0cONWpEMmovlVDRgLzEIFF9pTAminyU7WInab16In\/84VsznD7dt\/EmkQXRSTDqg0gi\/qd\/Pj0jzRLyDJYtW\/aKYN8I\/t2+fXtXxD1NJH83tB0qwg2FMXn8WDf42+l9oDA63UN+9s2YIYJIaBcu+PItjh3rC\/nmpIJoHStXrlQiiVRL\/jsJEb0ESYWRQxGiiJEm4p4OHTpUXvTfRuukDtEWWwiEEkaMEnPmzKlGjT\/++KMttjilkVDRgJxin5vtoDBq4l2IIM4g4uhF164iOiSVwBokNtVAIPEH06+BpWfPntIHoXhYPE0gVIQbCmMmtbENm9JY7CNAYbSPddQtDRki0q2b7wzh8OG+rBk6lp07dyqBxGgRpXnz5tK\/f38du0Kb40wglDBilJgrVy41ajxy5EicW3V2daGiATnbandYR2F0sB+x6bNNG19ORUSkmTxZpEkTBxtM00ggSgKhItx4WRhDBT2IEjFvi4AAhTECWHZeiugzjRqJzJkjkjatL84pQrSxkIAbCYQSRowSc+fOrUaNh\/0P0roRQkCfKIyJczKFMXHsQ7aMUFAtWlyUmTMzSebMvgg0FSo40FCaRAJxIhAqwg3EME+ePEocDx06FKfW9KgmVNADPazX20oKo8P8Z2SV2L\/\/D8mZc6NMnZpWihVzmJE0hwTiTIDCmBwohTHOD1kE1VEYI4Bl9aWYNnr44Ydl165dkj9\/flm+fLkULFjQ6mZZPwkknECoCDcYJebNm1eNGg8i9qGHSqigBx5CkLCuUhgThv7KhnGUoUqVKrJ\/\/34VLWbZsmVq+oiFBLxAgMKY3MsUxsQ9+RTGxLFPaplZJRzgBJqQUAKhItxglJgvXz41avzhhx8SaqPdjYcKemC3HV5sj8KYYK9\/\/PHHKvA2s0ok2BFsPqEEQgkjfo64u4iVimhKbisIgjF58mTZuHGjTEd8R79CYUyctymMiWOvftERaBth1JC6adasWSpTBQsJeI2AFyPcrF69Wjp27ChffPGFcvcnn3wi5ZAN\/J8SKuiB156NRPSXwpgI6iLin1WiZcuWMn78+CsCbCfILDZLAgkh8M0336i1dbwY4kXRzQXTwz169JC3335bdRMbixABqlGjRip2sFEwdXzLLbfIddddJ+fOnXMzEsf1jcKYAJdABNu1a8esEglgzyadSWD79u1y7733KmHABjTsynZbwfnkYcOGyZAhQ+TPP\/+U9OnTq+D53bp1CzpThBFkhQoV1AuzkfvUbUyc2h8Ko82ewS8FfhGYVcJm8GzO0QSQexAjJ6y1QwggIHh5TJ06taPtNmPc5cuX1eiwa9eugjVT\/O7XrVtXhg8frjYVBZZLly7JmDFj5IUXXlAvz148qmKGq5XXUBitpOtXN345OnXqJKNGjZJUqVLJpEmTpGnTpja1zmZIwPkEFi1aJM8884zgmALKPffco35fKleu7HzjQ1i4efNmlXP0s88+U1eULFlSRo8eLWXKlAl6BzbjdejQQTCCRkFmjTlz5qjzzSz2EaAw2sQab4D4BcEayuzZs9VOVBYSIIHkBJCBBS+RBw4cUB9iVyp+f3SaXkUou5dffllmzJgheClGrNeBAwcmW0c0eh\/s+ldeeUWee+457j1IwC8JhdEm6OfPn5cGDRrI888\/L5UqVbKpVTZDAnoS+Ouvv2TixImCfJ1I2HvNNdcINqlhk0rGjBkd2ylskkHibfzBJqJ06dKp33lstsmQIUMyuyO93rEdd5lhFEaXOZTdIQE3EcAxjt69e8uUKVPUehvyMg4aNCjkyCtRfceo8L333pPOnTsnBSLASHfs2LFqZ2lgMa7H5hvsUjXWHbG2ioAGLIklQGFMLH+2TgIkYILA1q1b1VLEhg0b1NVYq8P0aunSpU3cbe0lsA3rgp9++qlqqESJEmpttHz58kEb3rJli7re6Eu46621nrUHI0Bh5HNBAiSgBQFjlIXdmggujlFWw4YN1Q7W7Nmz294HJFHu06dP0mg2R44canTbrFmzoOuCSBLQt2\/fK0a\/vXr1Cnm97R1ig0kEKIx8GEiABLQiEOo84EsvvSRpkdXb4oK1Q4xWBwwYIDhmYqx\/9uvXT+0iDSyhrnf6eqnFGB1dPYXR0e6hcSRAAqEI7Nu3T+38nDt3rrqkUKFC8uqrr6pdrFYVHCnBlC6y4aCgLRy\/KFCgQNAmcT3OY37\/\/fdJ1+u2w9Yqlk6ul8LoZO\/QNhIggbAEVq1apdbskKUG5aGHHlJrfHfccUfYe81esG3bNtXG+vXr1S3FixeXkSNHhtxhHux62FSxYkWzTfK6BBKgMCYQPpsmARKID4GLFy\/K1KlTpXv37oLg22nSpJFWrVqpNb1MmTJF3QjqwhQpwjgiLFuWLFkE5wvbtm2rAnUEFmMX7RtvvGHq+qgN442WEqAwWoqXlZMACdhJAFFzsCHGrJCFsu3ChQsyYcKEpHOU4YQ28NxluOvtZMK2IidAYYycGe8gARJwOIEdO3aotUCEWENBgHJMZSIod7iC0WexYsVk165d6tLatWuruKYFCxYMeivWETHNakTqwVQu1hGLFCkSril+7lACFEaHOoZmkQAJxE4g2GaZcePGhT1Ej0g177\/\/vowYMUJq1KgR1BCkykLouqVLl6rPb7\/9drX5p3r16rEbzhoSSoDCmFD8bJwESMBqAsZxCRyPOHv2rArThgg1yHaBXIfBCtJCIbNHsOwegdO1N954o5pybdOmjSuygVjtDx3qpzDq4CXaSAIkEDMBHLDHWUcjsHfu3LnVWcTABMGhGjLWHXGI\/9SpU2qDT5MmTVQdWbNmjdk+VuAcAhRG5\/iClpAACdhAAKmgcLZw06ZNqjWktcL6I9JchSorV65U64jGuiPWEXFc484777TBYjZhNwEKo93E2R4JkEDCCSAgOUaOXbp0kWPHjqkQbsh+g0022bJlS7Jv9+7dKmHwkiVL1M9uu+02te5oZRCBhMOhAUJh5ENAAiTgWQKYEh08eLAaMSI1XObMmeX+++9XZxexDongATiKccMNN6g1yY4dO6oQcCzuJkBhdLd\/2TsSIAETBLDDFKK3bNmyK67GIf4WLVqoQAE43M\/iDQIURm\/4mb0kARIwQQAbaSCCGCVityki2NSpU8fEnbzETQQojG7yJvtCAiQQMwGI4pkzZ7jTNGaS+lZAYdTXd7ScBEiABEjAAgIURgugskoSIAESIAF9CVAY9fUdLScBEiABErCAAIXRAqiskgRIgARIQF8CFEZ9fUfLSYAESIAELCBAYbQAKqskARIgARLQlwCFUV\/f0XISIAESIAELCFAYLYDKKkmABEiABPQlQGHU13e0nARIgARIwAICFEYLoLJKEiABEiABfQlQGPX1HS0nARIgARKwgACF0QKorJIESIAESEBfAhRGfX1Hy0mABEiABCwgQGG0ACqrJAESIAES0JcAhVFf39FyEiABEiABCwhQGC2AyipJgARIgAT0JUBh1Nd3tJwESIAESMACAhRGC6CyShIgARIgAX0JUBj19R0tJwESIAESsIAAhdECqKySBEiABEhAXwIURn19R8tJgARIgAQsIEBhtAAqqyQBEiABEtCXAIVRX9\/RchIgARIgAQsIUBgtgMoqSYAESIAE9CVAYdTXd7ScBEiABEjAAgIURgugskoSIAESIAF9CVAY9fUdLScBEiABErCAAIXRAqiskgRIgARIQF8CFEZ9fUfLSYAESIAELCBAYbQAKqskARIgARLQlwCFUV\/f0XISIAESIAELCFAYLYDKKkmABEiABPQlQGHU13e0nARIgDx7\/j4AAAS0SURBVARIwAICFEYLoLJKEiABEiABfQlQGPX1HS0nARIgARKwgACF0QKorJIESIAESEBfAhRGfX1Hy0mABEiABCwgQGG0ACqrJAESIAES0JcAhVFf39FyEiABEiABCwhQGC2AyipJgARIgAT0JUBh1Nd3tJwESIAESMACAhRGC6CyShIgARIgAX0JUBj19R0tJwESIAESsIAAhdECqKySBEiABEhAXwIURn19R8tJgARIgAQsIEBhtAAqqyQBEiABEtCXAIVRX9\/RchIgARIgAQsIUBgtgMoqSYAESIAE9CVAYdTXd7ScBEiABEjAAgIURgugskoSIAESIAF9CVAY9fUdLScBEiABErCAAIXRAqiskgRIgARIQF8CFEZ9fUfLSYAESIAELCBAYbQAKqskARIgARLQlwCFUV\/f0XISIAESIAELCFAYLYDKKkmABEiABPQlQGHU13e0nARIgARIwAICFEYLoLJKEiABEiABfQlQGPX1HS0nARIgARKwgACF0QKorJIESIAESEBfAhRGfX1Hy0mABEiABCwgQGG0ACqrJAESIAES0JcAhVFf39FyEiABEiABCwhQGC2AyipJgARIgAT0JUBh1Nd3tJwESIAESMACAhRGC6CyShIgARIgAX0JUBj19R0tJwESIAESsIAAhdECqKySBEiABEhAXwIURn19R8tJgARIgAQsIEBhtAAqqyQBEiABEtCXAIVRX9\/RchIgARIgAQsIUBgtgMoqSYAESIAE9CVAYdTXd7ScBEiABEjAAgIURgugskoSIAESIAF9CVAY9fUdLScBEiABErCAAIXRAqiskgRIgARIQF8CFEZ9fUfLSYAESIAELCBAYbQAKqskARIgARLQlwCFUV\/f0XISIAESIAELCFAYLYDKKkmABEiABPQlQGHU13e0nARIgARIwAICFEYLoLJKEiABEiABfQlQGPX1HS0nARIgARKwgACF0QKorJIESIAESEBfAhRGfX1Hy0mABEiABCwgQGG0ACqrJAESIAES0JcAhVFf39FyEiABEiABCwhQGC2AyipJgARIgAT0JUBh1Nd3tJwESIAESMACAhRGC6CyShIgARIgAX0JUBj19R0tJwESIAESsIAAhdECqKySBEiABEhAXwIURn19R8tJgARIgAQsIEBhtAAqqyQBEiABEtCXAIVRX9\/RchIgARIgAQsIUBgtgMoqSYAESIAE9CVAYdTXd7ScBEiABEjAAgIURgugskoSIAESIAF9CVAY9fUdLScBEiABErCAAIXRAqiskgRIgARIQF8CFEZ9fUfLSYAESIAELCBAYbQAKqskARIgARLQlwCFUV\/f0XISIAESIAELCFAYLYDKKkmABEiABPQlQGHU13e0nARIgARIwAICFEYLoLJKEiABEiABfQlQGPX1HS0nARIgARKwgACF0QKorJIESIAESEBfAhRGfX1Hy0mABEiABCwgQGG0ACqrJAESIAES0JcAhVFf39FyEiABEiABCwhQGC2AyipJgARIgAT0JUBh1Nd3tJwESIAESMACAhRGC6CyShIgARIgAX0JUBj19R0tJwESIAESsIAAhdECqKySBEiABEhAXwL\/B6joNHM7cG1\/AAAAAElFTkSuQmCC","height":306,"width":363}}
%---
%[output:6b5c0416]
%   data: {"dataType":"text","outputData":{"text":"[15:44:13][INFO]  法科学参照データベース: data\/list\/forensic_challenge.csv\n","truncated":false}}
%---
%[output:74f7855e]
%   data: {"dataType":"text","outputData":{"text":"[15:44:13][INFO]    総エントリ数  : 55\n","truncated":false}}
%---
%[output:41640cec]
%   data: {"dataType":"text","outputData":{"text":"[15:44:13][INFO]    処方薬        : 25\n","truncated":false}}
%---
%[output:5f2e701b]
%   data: {"dataType":"text","outputData":{"text":"[15:44:13][INFO]    日用化学品    : 30\n","truncated":false}}
%---
%[output:00137179]
%   data: {"dataType":"text","outputData":{"text":"[15:44:13][INFO]    カテゴリ数 (8): analgesic | flavor | food_acid | other | solvent | stimulant | sugar | vitamin\n","truncated":false}}
%---
%[output:41c246b3]
%   data: {"dataType":"text","outputData":{"text":"[15:44:13][INFO]  参照エントリ 55 件の ECFP4 フィンガープリントを計算中...\n","truncated":false}}
%---
%[output:4f3261f2]
%   data: {"dataType":"text","outputData":{"text":"\r[----------]   2% ( 1\/55) FP\r[----------]   4% ( 2\/55) FP\r[#---------]   5% ( 3\/55) FP\r[#---------]   7% ( 4\/55) FP\r[#---------]   9% ( 5\/55) FP\r[#---------]  11% ( 6\/55) FP\r[#---------]  13% ( 7\/55) FP\r[#---------]  15% ( 8\/55) FP\r[##--------]  16% ( 9\/55) FP\r[##--------]  18% (10\/55) FP\r[##--------]  20% (11\/55) FP\r[##--------]  22% (12\/55) FP\r[##--------]  24% (13\/55) FP\r[###-------]  25% (14\/55) FP\r[###-------]  27% (15\/55) FP\r[###-------]  29% (16\/55) FP\r[###-------]  31% (17\/55) FP\r[###-------]  33% (18\/55) FP\r[###-------]  35% (19\/55) FP\r[####------]  36% (20\/55) FP\r[####------]  38% (21\/55) FP\r[####------]  40% (22\/55) FP\r[####------]  42% (23\/55) FP\r[####------]  44% (24\/55) FP\r[#####-----]  45% (25\/55) FP\r[#####-----]  47% (26\/55) FP\r[#####-----]  49% (27\/55) FP\r[#####-----]  51% (28\/55) FP\r[#####-----]  53% (29\/55) FP\r[#####-----]  55% (30\/55) FP\r[######----]  56% (31\/55) FP\r[######----]  58% (32\/55) FP\r[######----]  60% (33\/55) FP\r[######----]  62% (34\/55) FP\r[######----]  64% (35\/55) FP\r[#######---]  65% (36\/55) FP\r[#######---]  67% (37\/55) FP\r[#######---]  69% (38\/55) FP\r[#######---]  71% (39\/55) FP\r[#######---]  73% (40\/55) FP\r[#######---]  75% (41\/55) FP\r[########--]  76% (42\/55) FP\r[########--]  78% (43\/55) FP\r[########--]  80% (44\/55) FP\r[########--]  82% (45\/55) FP\r[########--]  84% (46\/55) FP\r[#########-]  85% (47\/55) FP\r[#########-]  87% (48\/55) FP\r[#########-]  89% (49\/55) FP\r[#########-]  91% (50\/55) FP\r[#########-]  93% (51\/55) FP\r[#########-]  95% (52\/55) FP\r[##########]  96% (53\/55) FP\r[##########]  98% (54\/55) FP\r[##########] 100% (55\/55) FP\n","truncated":false}}
%---
%[output:3a06d0cd]
%   data: {"dataType":"text","outputData":{"text":"[15:44:21][INFO]  フィンガープリント準備完了: 55 \/ 55 エントリ\n","truncated":false}}
%---
%[output:255e7cd6]
%   data: {"dataType":"text","outputData":{"text":"[15:44:21][INFO]  同定検索完了: 55 候補を 0.514 秒でランク付け\n","truncated":false}}
%---
%[output:143ce2b8]
%   data: {"dataType":"text","outputData":{"text":"[15:44:21][INFO]  --- 上位 10 件の同定候補 ---\n","truncated":false}}
%---
%[output:4b8f9e39]
%   data: {"dataType":"text","outputData":{"text":"[15:44:21][INFO]  順位     名前                            カテゴリ                 タニモト  評価\n","truncated":false}}
%---
%[output:3aa8ccf5]
%   data: {"dataType":"text","outputData":{"text":"[15:44:21][INFO]     1.  IMIPRAMINE                    (-)                1.0000  *** 確定一致 ***\n[15:44:21][INFO]     2.  TRIFLUOPERAZINE               (-)                0.3273  低類似性\n[15:44:21][INFO]     3.  NORTRIPTYLINE                 (-)                0.2609  低類似性\n[15:44:21][INFO]     4.  CYPROHEPTADINE                (-)                0.2444  低類似性\n[15:44:21][INFO]     5.  BUSPIRONE                     (-)                0.1695  低類似性\n[15:44:21][INFO]     6.  nicotine                      stimulant          0.1569  低類似性\n[15:44:21][INFO]     7.  salicylic acid                analgesic          0.1463  低類似性\n[15:44:21][INFO]     8.  benzoic acid                  food_acid          0.1316  低類似性\n[15:44:21][INFO]     9.  benzaldehyde                  flavor             0.1316  低類似性\n[15:44:21][INFO]    10.  BENZNIDAZOLE                  (-)                0.1311  低類似性\n","truncated":false}}
%---
%[output:38ac97d0]
%   data: {"dataType":"text","outputData":{"text":"          Name           Tanimoto    Category    IsDrug    MW_Da      LogP     TPSA_A2    HBD\n    _________________    ________    ________    ______    ______    ______    _______    ___\n\n    \"IMIPRAMINE\"               1        \"\"       true      280.42     3.875      6.48      0 \n    \"TRIFLUOPERAZINE\"    0.32727        \"\"       true      407.51    4.9456      9.72      0 \n    \"NORTRIPTYLINE\"      0.26087        \"\"       true      263.38    3.8264     12.03      1 \n    \"CYPROHEPTADINE\"     0.24444        \"\"       true      287.41    4.6979      3.24      0 \n    \"BUSPIRONE\"          0.16949        \"\"       true      385.51    2.0882     69.64      0 \n\n","truncated":false}}
%---
%[output:1f5d3cc8]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]  物性比較: 未知物質 vs 上位候補\n","truncated":false}}
%---
%[output:875d89b8]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    指標              未知物質     上位ヒット        差分\n","truncated":false}}
%---
%[output:2fc719c6]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    MW (Da)        280.4     280.4       0.0\n","truncated":false}}
%---
%[output:00baa8a8]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    LogP            3.88      3.88      0.00\n","truncated":false}}
%---
%[output:5c78a853]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    TPSA (A^2)       6.5       6.5       0.0\n","truncated":false}}
%---
%[output:9bffe0bf]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]  同定: フィンガープリント確定一致 (T = 1.000000)。\n","truncated":false}}
%---
%[output:1d08311c]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    MW と LogP の差分はともに ~0 のはず（同一化合物）。\n","truncated":false}}
%---
%[output:79e355f5]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]  法科学同定結果 -- FC-2026-0004\n","truncated":false}}
%---
%[output:683fddf2]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    化合物名        : IMIPRAMINE\n","truncated":false}}
%---
%[output:523ceffb]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    タニモトスコア  : 1.000000\n","truncated":false}}
%---
%[output:1fcc7548]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    SMILES（参照 DB）: CN(C)CCCN1c2ccccc2CCc2ccccc21\n","truncated":false}}
%---
%[output:1e7ee373]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    処方薬          : true\n","truncated":false}}
%---
%[output:43b2507c]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]  -------------------------------------------------\n","truncated":false}}
%---
%[output:8257319d]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]  ステータス: 処方薬・規制薬物が確認された。\n","truncated":false}}
%---
%[output:2beeaf99]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    有効な処方箋なしでの所持は違法となる場合がある。\n","truncated":false}}
%---
%[output:8e79d6ca]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAcYAAAF+CAYAAAAY+7oxAAAAAXNSR0IArs4c6QAAIABJREFUeF7tnQm4TdX7x99CyhBlniMSTSQy00CmUCnJ8EtknstQZJ6HzCQiRTIkY5kjIkOGIhkrQ6JkSErE\/\/9dp307zj3nnmnvc\/ba+\/s+z31c9+6z9rs+a9\/zPWutd73vDdeuXbsmNBIgARIgARIgAUXgBgojnwQSIAESIAES+I8AhZFPAwmQAAmQAAl4EaAw8nEgARIgARIgAQojnwESIAESIAES8E+AM0Y+GbYgkC9fPjl8+HCCL\/g\/zPtn4Tjq3Z5v26G2E60Pod7HuM7XZ9\/X+2MRSd8ieU0kfQn0Gt9xDjTGsfAz3H7xencQoDC6Y5xt3cuk3gANcQpXJAMJo3d7\/qAEEh\/jWt\/fmyWegRgE+sCQ1IDCx2D9DPZ6Mx+YUPvm754URzNHgm2FSoDCGCopXmcJgXDe+HxFKNQ3\/0hnnYHeqP2JdDj9CEcA\/Amj0R9\/4m\/8LFJ\/In1doIcj1A893q8PdUZpyQPJRkmAxzX4DMSTgNlvwuEIjhX9jrQ\/wcTDn1AkNSM2ZoyRfCCItA+hCGNSIo\/XRzOztGI82aZ7CXDG6N6xj2vPzVqCDNYJs9\/oI71fOOJn3CMpoQjkh+9sMtRZtTELNpOXr3gbPgcSbgpjsKeLv48VAQpjrEjzPoqAtyAGEwvvN9JI8QV6c06qPW9xMcOHpGZDwWZK3n567x36W04NlZGZ4ud9z6RYBxLsUJZNrfI3VF68zn0EKIzuG\/O49TjYUlogx6KZXYb7phpslhbJ8mQ0wpjUfmKwdpPiGawf4cw0fe8TStv+Zo2cMcbtT5M39iFAYeQjYTmBQMIWiWjB2WBvvEaHwm0\/FKGJRKSD+eHv9773CXXmG8qRjmD+mP1AhCJ4wVYPQh1zs31ne+4kQGF057jbotdWv0FH0n4krwkGM1ib4QqHr0gGEw1\/M3V\/PgdrJ1g\/fZdVA10fLJjI93XB+IXjF68lgVAIUBhDocRrLCFg5RteJG1H8ppgYEJpMx7CaKYIRsMglBl4KAyD+cDfk0A4BCiM4dDitaYSCLafF86yqRmzjEDLlZGKSDhv6IGWUwPtMYYyGzOWhn05huOXGQMe7H5cRjWDMtswkwCF0UyabCssAsHeML3f2MMRyVDaDUdIjVmNlT6EK4yhRqUa7Ya7\/BrWQAa52JtfIJEOZdZspk9siwSSIkBh5PMRNwLhCliwZbdgv\/fX0XBfE+z6YL9PCnawvcBgEarBxD4SgY\/04UhKDI02\/fU3lOMbkfrE15FAqAQojKGS4nWmEwhXGP29oUb6Zh+NgBkzWX85SSNddg0kFkmJeaABCZb5xle0fNuJtg\/BHpSk2Pub4QZrj78nAbMJUBjNJsr2SIAESIAEtCZAYdR6+Og8CZAACZCA2QQojGYTZXskQAIkQAJaE6Awaj18dJ4ESIAESMBsAhRGs4myPRIgARIgAa0JUBi1Hj46TwIkQAIkYDYBCqPZRNkeCZAACZCA1gQojFoPH50nARIgARIwmwCF0WyibI8ESIAESEBrAhRGrYePzpMACZAACZhNgMJoNlG2RwIkQAIkoDUBCqPWw0fnSYAESIAEzCZAYTSbKNsjARIgARLQmgCFUevho\/MkQAIkQAJmE6Awmk2U7ZEACZAACWhNgMKo9fDReRIgARIgAbMJUBjNJsr2SIAESIAEtCZAYdR6+Og8CZAACZCA2QQojGYTZXskQAIkQAJaE6Awaj18dJ4ESIAESMBsAhRGs4myPRIgARIgAa0JUBi1Hj46TwIkQAIkYDYBCqPZRNkeCZAACZCA1gQojFoPH50nARIgARIwmwCF0WyibI8ESIAESEBrAhRGrYePzpMACZAACZhNgMJoNlG2RwIkQAIkoDUBCqPWw0fnSYAESIAEzCZAYTSbKNsjARIgARLQmgCFUevho\/MkQAIkQAJmE6Awmk2U7ZEACZAACWhNgMKo9fDReRIgARIgAbMJUBjNJsr2SIAESIAEtCZAYdR6+Og8CZAACZCA2QQojGYTZXskQAIkQAJaE6Awaj18dJ4ESIAESMBsAhRGs4myPRIgARIgAa0JUBi1Hj46TwIkQAIkYDYBCqPZRNkeCZAACZCA1gQojFoPH50nARIgARIwmwCF0WyibI8ESIAESEBrAhRGrYePzpMACZAACZhNgMJoNlG2RwIkQAIkoDUBCqPWw0fnSYAESIAEzCZAYTSbKNsjARIgARLQmgCFUevho\/MkQAIkQAJmE6Awmk2U7ZEACZAACWhNgMKo9fDReRIgARIgAbMJUBjNJsr2SIAESIAEtCZAYdR6+Og8CZAACZCA2QQojGYTZXskQAIkQAJaE6Awaj18dJ4ESIAESMBsAhRGs4myPRIgARIgAa0JUBi1Hj46TwIkQAIkYDYBCqPZRNkeCZAACZCA1gQojFoPH50nARIgARIwmwCF0WyibI8ESIAESEBrAhRGrYePzpMACZAACZhNgMJoNlG2RwIkQAIkoDUBCqPWw0fn403gxIkTsmnTJsmXL58UKlRIUqZMGW+XeH8SIIEoCVAYowTIl7uXwKRJk6Rt27Zy+fLlBAg333yz3Hbbbdd9Zc+eXbJly+b35zly5KCYuvcRYs9tSoDCaNOBoVv2JnD27FnJnDmzEsVkyZJJihQp5OrVq\/L333+H7XjatGklQ4YMkilTJsmYMaP6Hl\/G9\/5+ftNNN4V9H76ABEggNAIUxtA48SoSuI5Aq1atZOLEiQJR++mnnyRNmjTq97\/\/\/rv88ssv8uuvv6qv06dPJ3zh\/\/gdfub9u2jEFOJsCKm3mN5www3y119\/SZEiRaRixYocPRIggTAIUBjDgMVLSQAEVq1aJZUrVxbM2rZs2SL3339\/VGD+\/PNPOXPmTMIX9i0htt4\/M77H744dOxbWzPTll1+Wt99+Oyof+WIScBMBCqObRpt9jZrAuXPn5L777pOjR4\/K4MGDpWvXrlG3GUkD58+fV7NOYwZqzEyNmej27dvl22+\/lQsXLghmjxDULFmyRHIrvoYEXEeAwui6IWeHoyHQoEEDmTlzppQqVUrWr1+v9hftbAgEwn5o37595Y033rCzq\/SNBGxDgMJom6GgI3YnsGDBAnnqqackVapUsnPnTilQoECCy99\/\/72cPHlSSpYsaatuLF26VGrUqCFZs2aVw4cPyy233GIr\/+gMCdiRAIXRjqNCn2xHAEuW9957r5w6dUomTJggLVu2TPAR0aiPPPKIbNy4Uc0mn3vuOVv5\/9BDD8lXX32VyG9bOUlnSMBGBCiMNhoMumJfAs8884zMnz9fHn\/8cVmxYoXatzNs+PDh0rlzZ3VW8ZtvvlFRonayefPmybPPPiu5c+eWgwcPqqMlNBIggcAEKIx8OkggCIFp06bJSy+9JOnSpVPClytXroRX7N27V4oVKyaILF2yZIlUr149bjwvXbok8HXz5s3qX8Mwo0Xk7J49e2T69OnSqFGjuPnIG5OADgQojDqMEn2MGwEcjYCo4LjEjBkzpH79+gm+XLlyRUqXLi1bt26VZs2aCTLhxNMQqZonTx4VbLNhwwYpU6ZMgjsQxBdffFHuvvtuJZA33nhjPF3lvUnA1gQojLYeHjoXTwKYaWHp9LPPPpPatWvLxx9\/fJ07vXr1UtGeefPmlV27dqnD\/vG27t27y8CBA+XJJ5+URYsWJbiDDD0FCxYUBAlhaRVLwzQSIAH\/BCiMfDJIIACB0aNHS4cOHVSqtt27d6sUcIbhnCAiUP\/55x9Zs2aNVKhQwRYccZ7xjjvukD\/++EO2bdsmDz74YIJfCBpq3bq1yoYD\/733SW3hPJ0gAZsQoDDaZCDohr0IfPfdd0pUsHf40UcfydNPP53gIPbysK+IJclXX31Vhg0bZivnO3bsKKNGjZK6devKhx9+eJ3fqAKCrDqffPKJVK1a1VZ+0xkSsAsBCqNdRoJ+2IYA9g7Lli2rglgaN24sU6dOvc43Q3hQZgrHIOx2NhD7ovnz5xf0A9lv7rrrrgT\/IeJdunRRCQpwvIRGAiSQmACFkU8FCfgQ6Nevn\/Ts2VNy5swpX3\/9tSoXZRiCWrBsiuAVCEvx4sVtyQ\/BQJMnT1bRtO+8806Cj1hixVIrUsetW7dOypcvb0v\/6RQJxJMAhTGe9Hlv2xHYsWOH2jvEbAvJwnFw3zDkHcX+3KFDh6RPnz5KPO1qyHKDYBvsIx44cEBFqxoG33v37i1PPPGELFu2zK5doF8kEDcCFMa4oeeN7UYAe4eYAeKsYvv27dU+nbc1adJELati7\/HLL7+0\/UF5HC354IMPVDHlMWPGJHTlt99+U7NGlMhCdRC7znrt9nzQH\/cQoDC6Z6zZ0yAEkL0GWWxw1g9Rm957h4sXL5aaNWtKypQp1b7iPffcY3ueSD6ANHYoj4UZJDLzGIZ9Ruw3IqgIwUU0EiCB\/whQGPk0kICIfPHFF2rvEEuP2Ed8+OGHE7hgPw6lpn7++WcZOXKkOsKhi+H85cKFC6Vbt24yaNCgBLeR8BznL1HMGPuoEFAaCZCAhwCFkU+C6wkgIAV7h8gjitJMOLTvbTj2MGfOHBWpioAVnbLGYOaLJOJIPvDDDz9cF0jUpk0bGT9+vKCU1vvvv+\/654AASMAgQGHks+B6As2bN1cV7osWLar2DrH0aBiqZUA40qRJo0pN3Xnnndrxqly5sqxcuTJRTUYUW8axDiQpwLlNfE8jARLgjJHPgMsJQDAQnQkxRM5TLJkahoPw+D+CVXDkAUcfdLS1a9eq6Nrbb79dfvzxRyXyhuGc5rvvviv4cPDWW2\/p2D36TAKmE+CM0XSkbFAXAki2DeHDgXgEoiCLjWHXrl1TlTI+\/fTTRHlHdemft5\/lypVTe6cjRoyQTp06JfwKy8cINkqePLk6hpIjRw4du0efScBUAhRGU3GyMZ0I1KtXT6VMQxUK7B0mS5Yswf2JEydKq1atJGPGjOr4RtasWXXqWiJfly5dKjVq1FD9QISqd8QtCivPnTtXCSaEk0YCbidAYXT7E+DS\/qNSBo4qpE6dWu0deu+voQLFAw88oM75zZ49WyAcTjAE4eCoCUS\/RYsWCV1CZRDsr6ZKlUpV30DSdBoJuJkAhdHNo+\/Svp86dUotoeJf1FBE+jTDUGoK+3Gff\/65qr2IGoxOMcwKIfI4prF\/\/361fGoYZpOYVfbo0UOQEo9GAm4mQGF08+i7tO+GCFSqVEmWL19+XfmlIUOGqDN\/2bNnV0uoCFhxikH08YEAicXfe+89adiwYULXEI2LxOLp0qVTxzrSp0\/vlG6zHyQQNgEKY9jI+AKdCUyZMkVefvll9caPg+25cuVK6A4EA+WkkBoOsycnlmVCBCoiUVEZBDUmvc9kYqaMCFYkAsCHAxoJuJUAhdGtI+\/CfmMmhL3D8+fPy6xZs+T5559PoICk4Zgxobhvy5YtBUV9nWiXL19WZajAwrfOJI6u4MwjCjJjrxF7jjQScCMBCqMbR92FfcYy4mOPPaZmRE899ZTMnz\/\/OgrYWxswYIDaf0MwCjLFONWQ7QZZb5DtB5lxkAbPsNKlS8umTZtk9OjR0q5dO6ciYL9IIEkCFEY+IK4gYGS3wWwIe4f41zBEamK2iAwwn332meNrFCI\/KjL4IIEByk4hwYFhCxYsUB8cUIsS5xq9swC54kFhJ0mAuVL5DLiBgHGGD33FuUXkPjXs4sWLqozUvn37pGvXrjJ48GA3IJGhQ4eq\/uIDAQouG4bEBlhuxocH7Mei1BaNBNxGgDNGt424C\/vbqFEjlSQbh9tPnDhxHQHjID+qS2B\/EWWl3GBInI6ajKgcgqMpyIxjGGo44qgKZpXIoep9rMMNbNhHEqAw8hlwPAEcS8B5RGR7uXDhwnWRmJghIYF4iRIl1CF3N1nv3r2lT58+UqVKFZX6zjAsKSNq9cCBA4mClNzEh311LwEKo3vH3jU9RxQqZot\/\/vmnOoZRrVo11\/Q9qY4iOTpmjcjwgwTqyIxj2OTJk1Xig8KFC6tlVbuW2jp+\/LjaG86dOzfHlARMI0BhNA0lG7IzAeQARZLwkiVLqqhLmodA586dZfjw4fLMM8\/IvHnzErDgWAfS5B05ckQWLVqkEqlbaQgIglCfOXMmyS8shSNoCNchcxFEEVG1jz76qKxatcpKF9m2iwhQGF002G7uKvbUcBTjl19+UZGnFStWdDOOhL6fPHlScYEwIeEB9loNGzVqlHTs2FEtM2\/evDkkXjgPevr0aSVyvv8aP8O+pu\/vMZuPxCCKWA6H\/e9\/\/1MltGgkEC0BCmO0BPl6bQggB2jPnj0FqeBWrFihjd9WO9q6dWuV0MA3sQEidiGamJkNHDhQZQnyFjTv7yF2EEIsW0diN998s2TIkEGl4PP+F9VNfH+G\/xs\/S5EihaqCYtSSxBIwI2kjGQG+xpsAhZHPg2sInDt3TvLkySP4F0cUcFSBJioLDvZeUYjZuxwV2EBksMQajuBB5G677TaVbzZbtmzq+2BfuDYag7BD4FE6DEdy6tSpE01zfK3LCVAYXf4AuK37yAGKROG1a9cWlJ6iBSaAQs5YWkWAC8424nt\/szfM8IwvJCGPl\/Xq1Uv69u2rkhIsXrxYpbejkUAkBCiMkVDja7QlgGVBRGJiTw2p31BtguafgHHMBTPr9evXX1fI2a7MUGx55MiRKqUf9pKRFJ5GAuESoDCGS4zXa0+gbdu2Mm7cOMfVWzRzYBYuXKhm1UgkvmPHDpV4XAdDIA6WhBGEg\/1JJC\/AmUwaCYRDgMIYDi1e6wgCR48eVUcREOq\/d+9eKVCggCP6ZVYnELmLmTQiVpFwHMEtOhmOmiDfK\/ZNkfN1w4YNam+ZRgKhEqAwhkqK1zmKAIJKpk6dqmozIvMN7T8CCFxBSarHH39cRe96V9\/QhROOfyA5OpaA8cEH\/2bJkkUX9+lnnAlQGOM8ALx9fAigckTBggVVRpeDBw8yc8q\/wzBt2jS1FIkgGpxr1DmjDKKPUXwZS8HFixeX1atXO7qcWHz+kpx5VwqjM8eVvQqBAAoVz549Wzp06KACNtxux44dk\/vvv19llUHS9QYNGmiPBMFW5cuXV9VTIJKffPKJ4DgJjQSSIkBh5PPhWgKISkXicJzdw1m+TJkyuZYFglaqVq0qy5cvl1q1agnqMjrFDh8+LGXLllWVVdA3nMtkxRCnjK41\/aAwWsOVrWpCADlAlyxZIt27d5f+\/ftr4rX5bo4ZM0bat2+vPhwgabjT9uN2794tFSpUUJl7UIYMUas67p2aP\/Js0R8BCiOfC1cTQA5QJBa\/9dZb5ccff5T06dO7jgeWGVGsGSngMJtCQnEnGsYaAUUoPdauXTsZPXq0E7vJPplAgMJoAkQ2oTcBVGbAYfABAwbI66+\/rndnwvQeSb+xzAjRaNy4sYrUdbKhAkeNGjXk0qVLKgNSly5dnNxd9i1CAhTGCMHxZc4hgGhFzCSQ1gx7jWnSpHFO54L0BMvHb7zxhjrvhyhU5DR1us2fP1+ee+45uXr1qkyaNEkd2aGRgDcBCiOfBxIQkTJlyqjE4ii1hL02N9jOnTvl4YcfFhyI\/\/TTT9W5P7cYqnGgmgiSjs+aNUueffZZt3Sd\/QyBAIUxBEi8xPkEUIwXEYuYOeGMIxJRO9mwlIizfQi0cet+m1GGDGON8XfTBwMnP9tm9I3CaAZFtqE9ARxXQAUJCAVq+jVt2lT7PiXVAeytDRs2TCU52L59u8qJ6kZ75ZVX5M0331T9X7lypZQuXdqNGNhnHwIURj4SJPAvASypvfDCC5IvXz51INypZ92wZIxD7ziugDyiWE51q+EDET4EIegIe8xIOl64cGG34mC\/\/yVAYeSjQAL\/EkBScbwp7t+\/Xz744AOpV6+e49j88ccfKqnBgQMHpEePHoLlRLcbxh3BOAjKyZEjh\/qwgNJkNPcSoDC6d+zZcz8EpkyZoqIUIZBYVkUuVSdZixYtVCRmkSJF1BENp++lhjp2SDpepUoVNWNE5RWIo9OSHITKgteJUBj5FJCAFwFEaKIaAw77Iy0aAnKcYthDQ4AJxHDr1q0s0uwzsOfPn1f5VLHnipyx69atc2XCB6c879H0g8IYDT2+1pEEkBEFicVLlCihZlVOsLNnz6o3e9SiHDp0qHTu3NkJ3TK9D6hFWa5cObXHXLFiRXWMhUnHTcds+wYpjLYfIjoYawJ\/\/fWX5M2bV37++WcVqYjD\/7obgooQXISoSywX4vwezT8BfHhANqAjR44Iculi79GpgVh8BvwToDDyySABPwQGDRqk0sNhaW3NmjVaM\/r444\/l6aefltSpU6vahFgqpiVNYM+ePSrp+OnTp1X5renTpztuv5nPQGACFEY+HSTghwD2m\/LkySNYgkQgBjLj+BoCNjDzsnMAC+oR3nfffYJ\/ke2lefPmHO8QCWzZskUee+wxlXS8TZs2Mnbs2BBfyct0J0Bh1H0E6b9lBHCcAYnFkXR68eLFie7TsWNHlUIOe1DIMer9lT17dsmWLZvfn+NIQMqUKS3z27thzBQxY6xUqZKqtchSS+Fhx2pBtWrVVNJxNyaZD4+Wc66mMDpnLNkTkwlgGQ3n2TBj+Oqrr1RpJm9r3bq1ypKDSNZwLW3atJIxY0ZV\/xAHy40v\/Azf4+fG9\/g\/vg93ZvrOO++ow+sopYUE4bly5QrXTV4voqKT69SpIzjviA8aH330Ebk4nACF0eEDzO5FR6BTp04ycuRIqVu3rnz44Yd+G8Oy66+\/\/iqIaISYGl\/4Gb7HMqbvz80SU\/iHA\/u+hiohSHEH35yarCC6kQ3v1Ugsj2LOsIULF0rNmjXDa4BXa0WAwqjVcNHZWBM4ceKEShEHIUNABnKLmmHYnzxz5kzCF+7z008\/Xfcz4\/f43bFjx+Tvv\/9OdGvUkcSxAm9DOSXsja1du1Zq166tllJp0RPAEjjGCIWcUdCZ5lwCFEbnji17ZhIBBKy8\/fbbKoR\/\/fr1JrUafjOY\/WFWasxEMQutWrWqWmb1NsxwMZPMnDmzyt6Df2nRE8BSND6gYAl93Lhx0TfIFmxLgMJo26GhY3YhgP3Fhx56SLnTqFEjNRtDSjVj7y9WgTSh8Pjuu+\/UXihmpDh\/99RTT4XyMl4TAgFE9+7evVsd38ExHppzCVAYnTu27JmJBJA\/E3Ua\/ZldolKvXLmijpXgmEGTJk0EeV9p5hDA8nSaNGkEyR\/OnTsnCJ6iOZcAhdG5Y8uemUgAb4aNGzdWOUbTpUunDntjKRPLmv72\/oLdGm+yvtGnvlGoRsSq8fNgM9M+ffpI7969VSTtrl275NZbbw3mBn8fIoHvv\/9e7TVjnxHLqTRnE6AwOnt82bsYEPj9998DRp9CPH33BaMVU+wZGsu4xjGPixcvyvDhwwUzm9WrVycKyIkBBkffYtmyZWo\/99FHH1V8ac4mQGF09viydzYl4BuVigjUQJGp+Pnx48fVIfOkDDPMypUry4wZM2zaa33dMhLLt2zZUiZMmKBvR+h5SAQojCFh4kUkEH8CmJl6n5f0jk5dunSpyoParl07wZs4zVwCrVq1kokTJ6pMRzjTSHM2AQqjs8eXvXMJAeRzRbmkO++8Uw4ePOiSXseum4hERjQqylChoDHN2QQojM4eX\/bOJQSwt4j8rCdPnlSJCAoXLuySnsemmzlz5lTL2YcPH1YlyWjOJkBhdPb4sncuIvDSSy\/JtGnTZODAgfLaa6+5qOfWdvWPP\/5QxzOQqxbfs5altbzt0DqF0Q6jQB9IwAQCRt3FUqVKycaNG01okU2AgJHgAQf8kYyd5nwCFEbnjzF76BICOLKByFQcQseyX9asWV3Sc2u7iSTs9evXVxU25s6da+3N2LotCFAYbTEMdIIEzCHw5JNPypIlS1Q5LJScokVPoFevXtK3b1\/p3r279O\/fP\/oG2YLtCVAYbT9EdJAEQieAZOdIeo6ySCiPRIuewPPPPy+zZ8+W9957Txo2bBh9g2zB9gQojLYfIjpIAqETQFQqolORPg4Zd1KnTh36i3mlXwKod7lz507ZvHmzlChRgpRcQIDC6IJBZhfdRaBkyZLqTRyV52vVquWuzpvc22vXrqmcsxcuXFC1MtOnT2\/yHdicHQlQGO04KvSJBKIgMGDAAOnRowcrbETB0Hjp0aNHJXfu3JIlSxb5+eefTWiRTehAgMKowyjRRxIIgwBqBuJoAZKNo+I8z92FAc\/n0lWrVkmlSpWkfPnysm7dusgb4iu1IkBh1Gq46CwJhEagQIECKjXcF198IaVLlw7tRbwqEYHx48dLmzZtpFmzZjJp0iQScgkBCqNLBprddBeBjh07qoTX3bp1k0GDBrmr8yb2FknZx44dKyNGjJBOnTqZ2DKbsjMBCqOdR4e+kUCEBD777DNVO\/Dee++Xb77ZFWErfBnKeK1cuVKdDa1evTqBuIQAhdElA81uuovAlStXpEaNA7J1a0H58ssbpUABd\/XfrN7myZNHjhw5IgcOHJD8+fOb1SzbsTkBCqPNB4jukUCkBBo0EJk5U2T4cJFXXom0Ffe+DsWk06RJI8mTJ1fJw\/EvzR0EKIzuGGf20oUE5swRqVtXpHx5EQZUhv8A7Nq1S4oUKSKFChWSb7\/9NvwG+AptCVAYtR06Ok4CSRO4cEEkY0aRK1dEcAQP39NCJzBnzhypW7eu1K5dW1C5hOYeAhRG94w1e+pCAk88IbJihcj06SKNGrkQQBRd7tevn\/Ts2VO6du0qgwcPjqIlvlQ3AhRG3UaM\/pJAGATGjxdp00bkmWdE5s0L44W8VBo0aCAzZ86UqVOnSuPGjUnERQQojC4abHbVfQSOHxfJlUsEucR\/+UXk5pvdxyDSHhcvXly2bdvGJAmRAtT4dRRGjQePrpNAKASKFRPZvl1k6VKRatVCeQWvAYF06dLJ+fPn5ddff5UMGTIQiosIUBhdNNjsqjsJ9O4t0qePSIuJKxhCAAAgAElEQVQWIhMnupNBuL1GjtkcOXJIxowZVfkumrsIUBjdNd7srQsJ7Ngh8uCDItmzixw7JnLDDS6EEGaXjcxBZcqUkQ0bNoT5al6uOwEKo+4jSP9JIAQCefOK\/PCDyJYtIsWLh\/ACl1\/y1ltvScuWLeWll16Sd955x+U03Nd9CqP7xpw9diEBRKYiQrVHD5F+\/VwIIMwuG0nYhwwZIl26dAnz1bxcdwIURt1HkP6TQAgEcJYRZxrvv19kF3OKByVWrVo1+fTTT2XBggVSq1atoNfzAmcRoDA6azzZGxLwS+DyZZHMmUXOnhU5fFgES6u0wATuvPNOOXz4sOzdu1fuvvtuonIZAQqjpgO+ceNGGT58uOCs1WuvvaZpL+h2LAkgbyryp44eLdKuXSzvrNe9Ll26JKlTp5YbbrhBJQ+\/6aab9OoAvY2aAIUxaoTxaWDZsmVStWpVVZ0dVdppJBCMACptdO8ugi2zVq2CXe3e3+\/evVvuu+8+ueuuu2Tfvn3uBeHinlMYNR18HDy+7bbbJEWKFHLu3DlJmTKlpj2h27Ei8M8\/IsmSxepu+t4HKeCaNGkilSpVkhXYnKW5jgCFUeMhx6dafLrFsmqpUqU07gldt5oAikO8+qrIo4+KTJ58\/d1WrxZp1kxk8WKRwoWt9sS+7V++fFkmTJggPXr0UAf7ccgf2xTdunWTm5lLz74DZ4FnFEYLoMaqyebNm8vbb7+t9hpfYSXaWGHX8j7vvSfyv\/95XF+yRKR69f+6sWiRCAIvkTauaFEtuxe104sWLZJXX31VDhw4oNrKmzevfP\/99wnf42\/s6aefjvo+bEAPAhRGPcbJr5fvvvuuyvr\/zDPPyDyWTtB4JK133RDGsmVFTpwQ+eYbkVtu8dzXzcL43XffqQ+Vn3zyiWJRsGBBGTFihFSvXl3WrVsn7du3FxQshj3yyCMyatQouR9nXmiOJkBh1Hh49+\/fr\/6Qs2TJIj+jEi2NBAIQMIRx\/XqRxx\/3LKv27+9eYfztt9+kT58+aun0ypUrar8edRdxsN87CvXq1asyY8YM6dy5s5w6dUpuvPFGqV+\/vhLPTJky8XlzKAEKo8YDe+3aNSWKSHKMM1dY\/qGRgD8ChjD+9JPnuMbIkSI7d4oUKuSuGSNEEME13bt3V1UzkidPrtK+9e\/fP0mhO3v2rCpWPHLkSPn777+VkPbq1Utat26t2qA5iwCFUfPxrFmzpixevFgVVH3hhRc07w3dt4qAtzCmS+cJssmTR2TtWk\/QjRv2GFetWqVmhAhYgz322GNK6BDEFqphlQZtGEuvOPyPNqpUqRJqE7xOAwIURg0GKSkXBw0aJK+\/\/rq0adNGxo4dq3lv6L5VBLyFMVs2EWxJP\/usyIwZImnTOlsYIWaINJ07d67CW6BAARkwYIA8CwAB7MKFC3LkyBEpHCBMFyLbrl07lRkHVqNGDRk9erTky5fPqiFkuzEkQGGMIWwrboUAgYoVK8qDDz4oX331lRW3YJsOIOArjOgSJjl79ogMGyZSr54nKhVlqcaM8RzfwCxS56QvxvInAmaQzSZNmjQq0AZHMIKd+8WHzWHDhqllVogojm\/4mnG8o2fPnqqgMfYmW7RooZZl0+LTBk1bAhRGbYfO4\/iff\/4p6dOnl3\/++UfwRoA\/fhoJ+BLwJ4wHD4pgFRFBlihHBWFEQeOFCz2vzpLFc8SjaVPMsvRhGihgBkKHPflQDEc3MAPEniREsV+\/fvLyyy9LMj8ZEk78f5hv7969ZcqUKYJ7Z8+eXe0\/Nm3aVAXr0PQjQGHUb8wSeVyiRAnZunWrrF69Wh7FCW4aCfgQ8CeMuKRXL5G+fT0XQxixEjh7tsjEiZ7gHMOKFfPMIhs0EEmVyr54UWC4Q4cO8vXXXysnsZqCGeMDDzwQttNIB4f9RFTZgBUqVEjtJz6BMiV+DCs2ON5hpGhEHmOIK5NvhI0+7i+gMMZ9CKJ3AH+MY8aMUZ9qsZdCIwEQuHRJ5PXXRbp1E8F7O2Z\/iErFHqNhuAYzxv37Ex\/wx8r822+LfPCByIULnlcgcAfJyJFrNQKtsWxgDh48qPbajX3EXLlyqSXNRo0aRX1PBLdBbBH5DcN+Iv7e\/EWBI1IcZ4ox48QeJRKRN2jQQIYOHSpZs2aN2hc2EBsCFMbYcLb0Lh9++KHUq1dPUENu6dKllt6LjetBAEJWu7YI0r3h3OKECSKffeZ\/xoeZIZZSkdjFz1aa\/P67CFLKvf++yKpV\/\/XfmEUiGDpeK\/iofoElUhQU\/uuvv1RVDIiS2WnccERj4sSJYuwn3nLLLSr4Bsc+\/O0nXrx4UYmhr1+h7G\/q8YQ520sKowPGF59M8+TJo\/YaT58+zX0NB4xpNF04c8aT8m3TJhFMUjBbLFIkmhb\/e+2334pgWXbKFJHTpz0\/v\/VWkeef9yy1QixjYcY+Ig7lI7mFMTODEGXznhKb7AzypyIxgLGfmCNHDhk4cKA0bNhQ+eBrR48eVeL5Pj5ViEj+\/PnV9UlFxJrsMpuLgACFMQJodnxJrVpTZd++B+Xjjx+QQoUS\/4Ha0Wf6ZD4BpHvDFhhSvt1xhwiKQ1gROHPxoghOPyAhuXfVs\/r1p0nZspfUmdpboZgW2ObNm9XS5pdffqlaxx479hFjuZe3bds2tZ+IBP4w7CdiebVkyZJ+e7xmzRrl8zcYmH\/PUMLne++91wJCbDJaAhTGaAna5PVGEVp8km\/SxCZO0Y2YEsAWWOXKIocOeQ7wQxRz5LDeBZQsnDZN5N13r0ry5Pnk+PEfVTWKJ598Upo1ayaPYy3XBDt27JjaR0SKNuzlBZutmXDLJJuAD5gJ+s5aA0W\/Gll3EAeAbFVG1p1Ax0Gs9p\/tByZAYXTI0zFqlEjHjh5RhDjS3EUAEaVVq4qcOoUZlAhyYmfIEFsGly5dkfnz58rkyZNl7dq1SrxgmBXhqAOWG5FKLVwz9uuwZ4fjSalSpZK2bduqQDM7HE\/y3ecMdl7SN0\/r7bffrvYukaTD33GQcHnx+ugJUBijZ2iLFhA88fDDntyX2AeiuYfAunWew\/jnziHNmciCBfELhjGoI0oUMzvkJcU+GwyH6pHCELNIpGPztyfnPWr+Ijzr1Kmjgm2wp243842MRYadN998U0Wx+jNU9sDy6vLly9WvixYtqpaEy5cvb7euuc4fCqNDhvzyZZH06XHg3zNr8Bdd6JCushteBJDnFMvoGHdkr5k+XSRFCvsgQuIJnC1E3dCPP\/5YHZiHQTSaNGmiyqZlzpw5kcM4lwvRMPbwihUrpkSjLOpm2dxwnhi+GzlZsZQM3++55x6\/nuM4CPYrjfqPEFKkd7wDm8S0uBCgMMYFuzU3rVBB5PPPExeiteZubDXeBJDn9KWXRPChqGVLkXHjROycaOX48eNqFvnWW2\/JDz\/8oPAhjVrlypXVeUMUAkaEad++fbXPIuNbxSNFihTSsmVL1bd0OAzqY8ZxkDfeeEN+\/\/13wXEQfHDA\/\/19cIj3s+f0+1MYHTTCOMg9ZIhI9+7\/1dpzUPfYFS8CyBffoYPI1asiXbuKDB6sDx7MIpFNBnuRqFJhzCJxAP7cuXNqHxHC0KlTJ5XXFGcTdTVjP3H8+PEqbWOGDBmU2AXaT0SAUZcuXQRnk7GUjAAdZNvB9bTYEaAwxo615XdCjksc6n7kEZE1ayy\/nVY3wKfv+fPnq2TSlSpVUsEOWJ7T0fDhBx+CcGwOCcBfeUXHXnh8Rp7R9957TyZNmpRwHhFpDQNlltG1pzt37lTLpZ9jSUdEJf3H8mq5cuX8dgnLz4jqRWAPgo3wLy12BCiMsWNt+Z1+\/VUE2zXIZXn2rAjrp3qQo0LCNJwn8DEctq5du7bUqlVLnYGze0QgZhwdO3aXVat6ysGDqWTqVE8mGyeYEW2KpVV8eHGq+dtPHDdunN9gIiyvGlVAjAhfp3KxW78ojHYbkSj9uesukQMHRLZti10WkihdtuzleDNBlhJ8wRARiYPYiB5csWKFmq0YhiUupNTDp\/SqVava4hiANxi8SeK4w5w5cyRnzqIyadJmqVbNRlE2UY4ixBBnHyEESO3mZMOHAMyIkcsVdR8xI+zcubM6D4klZMOQ3cf4sEZhjO0TQWGMLW\/L7\/bii57IROxBuXlbArMrBDtgHwv7NIiKRASk95vOjh07BJ\/gITZGwVn8Hm9OOE4AkYSYxjv5M5bRcExh2bJlKnADPgdagrP8AbPoBhBDcIc4Qji8xwmBOijf5LQoTaRyxH4inj8IH45tIBDJVxjRdzzPtNgRoDDGjnVM7jRpkkiLFiJI7DxzZkxuabubYHaFigaotIBP46h2gFlgUobKCRAcvGbTpk2qrh4Mn9iR5gsiiWXXggULxrS\/Z86cUefgcGwBtQQRtILzbk6zQMKIAsD4MID0cgjMcaJh33HRokUyfPjw67pnzBgpjLEfdQpj7JlbekekYkQZIRyB+v57S29ly8Yxu0LYP5ZKkVQdYhfu2Tek64IAQSRXrlx53Z5X4cKFlUhCrMqUKRP0kHo0kHB0AbX\/UFsQB9rRp7uwVu5AM\/YYMWtEphvD3CCMgYYTs0SsduDDmRG568Cht2WXKIy2HJbIncJE5\/bbPVlQjh2LTa7MyL0195WYXVWvXl3N+LD8CXErEmVZCbxJ48A2RBIiexZRTf8aav5hJgqRhIAhcMQsw2FvLKthPxQFciGKOXPmNKt527UTSBgxS8QHHMwavdnbrgMWOERhtABqiE1SGEMEpdNlVaqIIMsUqh\/UqaOT55H7ikAaiBOqF2AvCkKC7CpmGt6oILoQSWRxMVKd4R7IAYoMJxBJRLn6O8Qdqi\/ImIK+oMQRgoVw1i+jw1MZ4QMIziv6Hk1wszBilojEAJg1XkYWB1rMCFAYY4Y6djdCEGbv3iKdOomMGBG7+8brTtgfxOzq0KFDgqVOiCIqL1hte\/bsUSK5ZMkS+Qrl7v81RFYiOAYiiaCZcHxBSSVEx+JgOM7zLViwwG8hXKv7Fuv2AwkjZon40IFZI1YE3GQUxviNNoUxfuwtu\/PBg57SQygN5yf7lGX3jUfD27dvV8uZp06dUnX5MLvC0YtYGyInFy5cqERy3bp1133Ch1ijMC32JpNKKrB06VJ1HZYVEegza9YsFaXpBsPeMKpSYNaIIwyGuVkYMUvE8jxmjQgoo8WOAIUxdqwtvxOq\/Awd6qms0KqVJzOKt40c6ankjsw4TjAIEJYtsdyG4xWYXdmhDBFmexBIfOGIBXJfGpYvXz41k4RIVqxYUS2TwT744AN58cUXlaDiX+OYiRPGKZQ+BBJGzBJRlgmzRnB1k1EY4zfaFMb4sTf9zjjqZGS7ef\/9xFlRUFC9eXNPGjHdDYEwdevWVbOrevXqyfTp09Una7sZlgixtGvMJn9FeqJ\/Dcmh8aaPiFNEvyI8H4e8Bw0aZGm0q90YwR\/MEtOmTas+2Hh\/kHCzMBqZb5yeDciOzyOF0Y6jEqFPhjCi\/BQCJL\/7DkEh\/zXmFGFEhQakecMnahziR0otnPWyu0H4jKQCs2fPFtTjMwy1CUeMGCEdUW3ahRZIGDFLxNI4PkCcPn3aVWTclA3IbgNLYbTbiEThjyGMr70mMnGiyHPPieDAv2FOEEbUqUOtO2N2NVinshI+Y4vjJBBDnNWrUqWKKknkVsMsEYf4MWsED8MojO5Ik2e3557CaLcRicIfQxixVIpVRUSlrl8vUrq0p1HdhXHIkCHSrVs3tcyIKu6v6FxWIopxduJLAwkjZok4qoJZo\/cytBMZ+PYpUDYgN\/Q93n2kMMZ7BEy8v7cwtm+P0jaewrU4SYC9R12FEecHW7durUoTIQsI8p5iKZXmHAKBMtxQGBPnj3XOqNu3JxRG+45N2J55C+Orr4qg9FvFiv\/V7NNRGBGAgOru2JPD+UBEbyLlG81ZBAIJI2aJmTJlUrNGpOpzkwXKBuQmBvHqK4UxXuQtuK+vMOIWjRqJLFggsn+\/CNJsIioVCQAKFRKpUEHkySdFqlUTsWORdDdUlbDgMdCyyUAZbiiMqVTVEe\/8sVoOsGZOUxg1G7Ck3PUnjCdPitx9tycQZ9YsjzDiHGP16v+1hHOPSCNXq5ZHJJFrNd7mlqoS8eZsl\/sHEkbMEnGsBbNGJHFwkwXKBuQmBvHqK4UxXuQtuK8\/YcRtxo8X6dDBc0P8i+Ccw4dFFi\/25FPduFEEyQFgyZJ5MuZgJvnUU55ZZqzNTVUlYs3WrvcLlOGGwpg4f6xdx9BJflEYHTCaODOOKFTM\/BBkA+HDHqNhqLhRqpTIli2en\/se8McH8WXLPCK5YoWId\/apwoU9IlmjhkiZMomz6ZiNz21VJczmp2t7gYQRs0TUocSs8SSWP1xkgbIBuQhB3LpKYYwbenNufOKEyBNPePYIUVEDNWxxjrFp0+vbhyjWq+f5OX4fyP74Q2TNGo9ILlrkKV9lWO7cniVXiCTuaWKVJXULN1aVMOcp0L+VQBluKIyJ88fqP9r27wGF0f5jFNDDfftEKlcWOXLEU5x45UqRzJnN6xCWZjdt8ojk\/Pme+o6GlSlTS7JnT6nyfiLhNQ5nR2Pr169X+UOx1+SmqhLRMHPSawMJI2aJqK2JWSOW2N1kgbIBuYlBvPpKYYwX+Sjvu327SNWqIlgGLVFC5JNPRKwsKoE9yK1bPRGua9f+Jlu2ZBacL4ShAgRqEdasWVN94U0sHHNzVYlwODn52kAZbiiMifPHOvk5sEvfKIx2GYkw\/Fi3zhNBimXOxx7ziBX2F2Np2AtctGiRqiCxdu1aQe04GHKWFi1aVM0kn3vuOVUfMSmbOXOmNG7c2LVVJWI5Zna+VyBhxCwxW7ZsataIYtRuskDZgNzEIF59pTDGi3yE90Ukad26In\/+6dkznD7dE3gTT0N2Esz6IJLI\/+ldT88os4Q6g6VLl74u2TeSf7dv394ReU\/jyd8J9w6U4YbCmDh\/rBPG2+59oDDafYS8\/JsxQwSZ0C5f9tRbHDvWk\/LNToZsHatWrVIiiVJL3pGEyF6CosKooQhRxEwTeU+HDh0qr3qH0dqpQ\/QlJgQCCSNmidmzZ1ezxp9++ikmvtjlJoGyAdnFPyf7QWHUZHQhgjiDiKMXXbuK6FBUAnuQCKqBQOILy6++1rNnT+mDVDw0VxMIlOGGwphOBbYhKI0WOwIUxtixjvhOQ4aIdOvmOUM4fLinaoaOtmvXLiWQmC3CmjVrJv3799exK\/TZZAKBhBGzxBw5cqhZ4\/Hjx02+q72bC5QNyN5eO8M7CqONxxFBn61be2oqIiPN5MkijRvb2GG6RgIREgiU4cbNwhgo6UGEiPmyMAhQGMOAFctLkX2mYUOROXNEUqb05DlFijYaCTiRQCBhxCwxZ86catZ4zPsgrRMh+PSJwhi\/QaYwxo99wDsjFVTz5ldk5sx0kj69JwNNuXI2dJQukYBJBAJluIEY5sqVS4nj0aNHTbqbHs0ESnqgh\/d6e0lhtNn4GVUlDh36U7Jn3yRTp6aUIkVs5iTdIQGTCVAYEwOlMJr8kIXRHIUxDFhWX4ployeeeEL27NkjefPmlRUrVkj+\/Pmtvi3bJ4G4EwiU4QazxNy5c6tZ4xHkPnSRBUp64CIEcesqhTFu6K+\/MY4yVKpUSQ4dOqSyxSxfvlwtH9FIwA0EKIyJR5nCGL8nn8IYP\/YJd2ZVCRsMAl2IK4FAGW4wS8yTJ4+aNf74449x9THWNw+U9CDWfrjxfhTGOI\/6559\/rhJvs6pEnAeCt48rgUDCiJ8j7y5ypSKbktMMSTAmT54smzZtkunI7+hlFMb4jTaFMX7s1R86Em0jjRpKN82aNUtVqqCRgNsIuDHDzZo1a6Rjx47y9ddfq+HesGGDlEE18H8tUNIDtz0b8egvhTEe1EXEu6pEixYtZPz48dcl2I6TW7wtCcSFwN69e9XeOj4Y4oOikw3Lwz169JD3339fdROBRcgA1bBhQ5U72DAsHd9xxx1yyy23yMWLF52MxHZ9ozDGYUgggu3atWNViTiw5y3tSWDHjh3y4IMPKmFAABqisp1mOJ88bNgwGTJkiPz111+SOnVqlTy\/W7dufleKMIMsV66c+sBs1D51GhO79ofCGOORwR8F\/hBYVSLG4Hk7WxNA7UHMnLDXDiGAgODDY\/LkyW3tdyjOXbt2Tc0Ou3btKtgzxd9+nTp1ZPjw4SqoyNeuXr0qY8aMkVdeeUV9eHbjUZVQuFp5DYXRSrpebeOPo1OnTjJq1ChJliyZTJo0SZo0aRKju\/M2JGB\/AosXL5YXX3xRcEwB9sADD6i\/l4oVK9rf+QAebtmyRdUc\/fLLL9UVxYsXl9GjR0upUqX8vgLBeB06dBDMoGGorDFnzhx1vpkWOwIUxhixxidA\/IFgD2X27NkqEpVGAiSQmAAqsOBD5OHDh9UvEZWKvx+dlleRyu7111+XGTNmCD4UI9frwIEDE+0jGr33d\/0bb7whL7\/8MmMP4vBHQmGMEfRLly5J\/fr1pW3btlKhQoUY3ZW3IQE9Cfz9998yceJEQb1OFOy96aabBEFqCFJJmzatbTuFIBkU3sYXgohSpUql\/uYRbJMmTZpEfod7vW077jDHKIwOG1B2hwScRADHOHr37i1TpkxR+22oyzho0KCAM6949R2zwnnz5knnzp0TEhFgpjt27FgVWeprxvUIvkGUqrHviL1VJDSgxZcAhTG+\/Hl3EiCBEAhs27ZNbUVs3LhRXY29OiyvlixZMoRXW3sJfMO+4BdffKFuVKxYMbU3WrZsWb833rp1q7re6Euw6631nq37I0Bh5HNBAiSgBQFjloVoTSQXxyyrQYMGKoI1S5YsMe8Diij36dMnYTabLVs2Nbtt2rSp331BFAno27fvdbPfXr16Bbw+5h3iDRMIUBj5MJAACWhFINB5wNdee01Soqq3xYa9Q8xWBwwYIDhmYux\/9uvXT0WR+lqg6+2+X2oxRls3T2G09fDQORIggUAEDh48qCI\/586dqy4pUKCAvPnmmyqK1SrDkRIs6aIaDgz3wvGLfPny+b0lrsd5zB9++CHhet0ibK1iaed2KYx2Hh36RgIkEJTA6tWr1Z4dqtTAHn\/8cbXHd8899wR9bagXbN++Xd1j\/fr16iVFixaVkSNHBoww93c9fCpfvnyot+R1cSRAYYwjfN6aBEjAHAJXrlyRqVOnSvfu3QXJt1OkSCEtW7ZUe3rp0qWL+CZoC0ukSOOItGwZMmQQnC9s06aNStTha0YU7TvvvBPS9RE7xhdaSoDCaCleNk4CJBBLAsiag4CYUIUskG+XL1+WCRMmJJyjDCa0vucug10fSya8V\/gEKIzhM+MrSIAEbE5g586dai8QKdZgSFCOpUwk5Q5mmH0WKVJE9uzZoy6tVauWymuaP39+vy\/FPiKWWY1MPVjKxT5ioUKFgt2Kv7cpAQqjTQeGbpEACURPwF+wzLhx44Ieokemmo8++khGjBgh1apV8+sISmUhdd2yZcvU7++++24V\/FO1atXoHWcLcSVAYYwrft6cBEjAagLGcQkcj7hw4YJK04YMNah2gVqH\/gxloVDZw191D9\/l2ttvv10tubZu3doR1UCsHg8d2qcw6jBK9JEESCBqAjhgj7OORmLvnDlzqrOIvgWCA93I2HfEIf6zZ8+qAJ\/GjRurNjJmzBi1f2zAPgQojPYZC3pCAiQQAwIoBYWzhZs3b1Z3Q1kr7D+izFUgW7VqldpHNPYdsY+I4xr33ntvDDzmLWJNgMIYa+K8HwmQQNwJICE5Zo5dunSRkydPqhRuqH6DIJvMmTMn+Ldv3z5VMHjp0qXqZ3fddZfad7QyiUDc4dABoTDyISABEnAtASyJDh48WM0YURouffr08sgjj6izi9iHRPIAHMW47bbb1J5kx44dVQo4mrMJUBidPb7sHQmQQAgEEGEK0Vu+fPl1V+MQf\/PmzVWiABzup7mDAIXRHePMXpIACYRAAIE0EEHMEhFtigw2tWvXDuGVvMRJBCiMThpN9oUESCBqAhDF8+fPM9I0apL6NkBh1Hfs6DkJkAAJkIAFBCiMFkBlkyRAAiRAAvoSoDDqO3b0nARIgARIwAICFEYLoLJJEiABEiABfQlQGPUdO3pOAiRAAiRgAQEKowVQ2SQJkAAJkIC+BCiM+o4dPScBEiABErCAAIXRAqhskgRIgARIQF8CFEZ9x46ekwAJkAAJWECAwmgBVDZJAiRAAiSgLwEKo75jR89JgARIgAQsIEBhtAAqmyQBEiABEtCXAIVR37Gj5yRAAiRAAhYQoDBaAJVNkgAJkAAJ6EuAwqjv2NFzEiABEiABCwhQGC2AyiZJgARIgAT0JUBh1Hfs6DkJkAAJkIAFBCiMFkBlkyRAAiRAAvoSoDDqO3b0nARIgARIwAICFEYLoLJJEiABEiABfQlQGPUdO3pOAiRAAiRgAQEKowVQ2SQJkAAJkIC+BCiM+o4dPScBEiABErCAAIXRAqhskgRIgARIQF8CFEZ9x46ekwAJkAAJWECAwmgBVDZJAiRAAiSgLwEKo75jR89JgARIgAQsIEBhtAAqmyQBEiABEtCXAIVR37Gj5yRAAiRAAhYQoDBaAJVNkgAJkAAJ6EuAwqjv2NFzEiABEiABCwhQGC2AyiZJgARIgAT0JUBh1Hfs6DkJkAAJkIAFBCiMFkBlkyRAAiRAAvoSoDDqO3b0nARIgARIwAICFEYLoLJJEiABEiABfQlQGPUdO3pOAiRAAiRgAQEKowVQ2SQJkAAJkIC+BCiM+o4dPScBEiABErCAAIXRAqhskgRIgARIQF8CFEZ9x46ekwAJkAAJWECAwmgBVDZJAiRAAiSgLwEKo75jR89JgARIgAQsIEBhtAAqmyQBEiABEtCXAIVR37Gj5yRAAiRAAhYQoDBaAJVNkgAJkAAJ6EuAwqjv2NFzEiABEiABCwhQGC2Ayrly8ukAAAPlSURBVCZJgARIgAT0JUBh1Hfs6DkJkAAJkIAFBCiMFkBlkyRAAiRAAvoSoDDqO3b0nARIgARIwAICFEYLoLJJEiABEiABfQlQGPUdO3pOAiRAAiRgAQEKowVQ2SQJkAAJkIC+BCiM+o4dPScBEiABErCAAIXRAqhskgRIgARIQF8CFEZ9x46ekwAJkAAJWECAwmgBVDZJAiRAAiSgLwEKo75jR89JgARIgAQsIEBhtAAqmyQBEiABEtCXAIVR37Gj5yRAAiRAAhYQoDBaAJVNkgAJkAAJ6EuAwqjv2NFzEiABEiABCwhQGC2AyiZJgARIgAT0JUBh1Hfs6DkJkAAJkIAFBCiMFkBlkyRAAiRAAvoSoDDqO3b0nARIgARIwAICFEYLoLJJEiABEiABfQlQGPUdO3pOAiRAAiRgAQEKowVQ2SQJkAAJkIC+BCiM+o4dPScBEiABErCAAIXRAqhskgRIgARIQF8CFEZ9x46ekwAJkAAJWECAwmgBVDZJAiRAAiSgLwEKo75jR89JgARIgAQsIEBhtAAqmyQBEiABEtCXAIVR37Gj5yRAAiRAAhYQoDBaAJVNkgAJkAAJ6EuAwqjv2NFzEiABEiABCwhQGC2AyiZJgARIgAT0JUBh1Hfs6DkJkAAJkIAFBCiMFkBlkyRAAiRAAvoSoDDqO3b0nARIgARIwAICFEYLoLJJEiABEiABfQlQGPUdO3pOAiRAAiRgAQEKowVQ2SQJkAAJkIC+BCiM+o4dPScBEiABErCAAIXRAqhskgRIgARIQF8CFEZ9x46ekwAJkAAJWECAwmgBVDZJAiRAAiSgLwEKo75jR89JgARIgAQsIEBhtAAqmyQBEiABEtCXAIVR37Gj5yRAAiRAAhYQoDBaAJVNkgAJkAAJ6EuAwqjv2NFzEiABEiABCwhQGC2AyiZJgARIgAT0JUBh1Hfs6DkJkAAJkIAFBCiMFkBlkyRAAiRAAvoSoDDqO3b0nARIgARIwAICFEYLoLJJEiABEiABfQlQGPUdO3pOAiRAAiRgAQEKowVQ2SQJkAAJkIC+BCiM+o4dPScBEiABErCAAIXRAqhskgRIgARIQF8CFEZ9x46ekwAJkAAJWECAwmgBVDZJAiRAAiSgLwEKo75jR89JgARIgAQsIEBhtAAqmyQBEiABEtCXAIVR37Gj5yRAAiRAAhYQoDBaAJVNkgAJkAAJ6EuAwqjv2NFzEiABEiABCwhQGC2AyiZJgARIgAT0JUBh1Hfs6DkJkAAJkIAFBCiMFkBlkyRAAiRAAvoSoDDqO3b0nARIgARIwAICFEYLoLJJEiABEiABfQn8H6pWOliDcHpgAAAAAElFTkSuQmCC","height":306,"width":363}}
%---
%[output:634f9ed1]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]  背景（イミプラミン）:\n","truncated":false}}
%---
%[output:15a5e931]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    プロトタイプ TCA、Roland Kuhn（Geigy、1957 年）が発見。\n","truncated":false}}
%---
%[output:0d512eaf]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    作用機序: NE と 5-HT 再取り込み阻害。\n","truncated":false}}
%---
%[output:8d24052b]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    治療域: 150-300 ng\/mL; 中毒 > 1000 ng\/mL。\n","truncated":false}}
%---
%[output:0b646497]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    TCA 過量投与: QRS 拡大・痙攣・心臓不整脈。\n","truncated":false}}
%---
%[output:4f11fda6]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]    DB 内の関連 TCA: NORTRIPTYLINE (FC032) -- アミトリプチリンのデスメチル代謝物。\n","truncated":false}}
%---
%[output:7dd6d15c]
%   data: {"dataType":"text","outputData":{"text":"[15:44:22][INFO]  MACCS キーによる同定のクロスバリデーション...\n","truncated":false}}
%---
%[output:50ab920b]
%   data: {"dataType":"text","outputData":{"text":"[15:44:23][INFO]  フィンガープリントクロスバリデーション:\n","truncated":false}}
%---
%[output:5e04b032]
%   data: {"dataType":"text","outputData":{"text":"[15:44:23][INFO]    ECFP4（Morgan r=2、2048 ビット）: T = 1.000000\n","truncated":false}}
%---
%[output:793db436]
%   data: {"dataType":"text","outputData":{"text":"[15:44:23][INFO]    MACCS（166 専門家ルール）        : T = 1.000000\n","truncated":false}}
%---
%[output:58a8690b]
%   data: {"dataType":"text","outputData":{"text":"[15:44:23][INFO]  クロスバリデーション: ECFP4 と MACCS の両方が同一フィンガープリントを確認。\n","truncated":false}}
%---
%[output:84682ba5]
%   data: {"dataType":"text","outputData":{"text":"[15:44:23][INFO]    同定信頼度: 最大。\n","truncated":false}}
%---
%[output:4e02eb69]
%   data: {"dataType":"text","outputData":{"text":"[15:44:23][INFO]  MACCS 同定検索を 55 化合物に対して実行中...\n","truncated":false}}
%---
%[output:76ddb101]
%   data: {"dataType":"text","outputData":{"text":"\r[----------]   2% ( 1\/55) MACCS\r[----------]   4% ( 2\/55) MACCS\r[#---------]   5% ( 3\/55) MACCS\r[#---------]   7% ( 4\/55) MACCS\r[#---------]   9% ( 5\/55) MACCS\r[#---------]  11% ( 6\/55) MACCS\r[#---------]  13% ( 7\/55) MACCS\r[#---------]  15% ( 8\/55) MACCS\r[##--------]  16% ( 9\/55) MACCS\r[##--------]  18% (10\/55) MACCS\r[##--------]  20% (11\/55) MACCS\r[##--------]  22% (12\/55) MACCS\r[##--------]  24% (13\/55) MACCS\r[###-------]  25% (14\/55) MACCS\r[###-------]  27% (15\/55) MACCS\r[###-------]  29% (16\/55) MACCS\r[###-------]  31% (17\/55) MACCS\r[###-------]  33% (18\/55) MACCS\r[###-------]  35% (19\/55) MACCS\r[####------]  36% (20\/55) MACCS\r[####------]  38% (21\/55) MACCS\r[####------]  40% (22\/55) MACCS\r[####------]  42% (23\/55) MACCS\r[####------]  44% (24\/55) MACCS\r[#####-----]  45% (25\/55) MACCS\r[#####-----]  47% (26\/55) MACCS\r[#####-----]  49% (27\/55) MACCS\r[#####-----]  51% (28\/55) MACCS\r[#####-----]  53% (29\/55) MACCS\r[#####-----]  55% (30\/55) MACCS\r[######----]  56% (31\/55) MACCS\r[######----]  58% (32\/55) MACCS\r[######----]  60% (33\/55) MACCS\r[######----]  62% (34\/55) MACCS\r[######----]  64% (35\/55) MACCS\r[#######---]  65% (36\/55) MACCS\r[#######---]  67% (37\/55) MACCS\r[#######---]  69% (38\/55) MACCS\r[#######---]  71% (39\/55) MACCS\r[#######---]  73% (40\/55) MACCS\r[#######---]  75% (41\/55) MACCS\r[########--]  76% (42\/55) MACCS\r[########--]  78% (43\/55) MACCS\r[########--]  80% (44\/55) MACCS\r[########--]  82% (45\/55) MACCS\r[########--]  84% (46\/55) MACCS\r[#########-]  85% (47\/55) MACCS\r[#########-]  87% (48\/55) MACCS\r[#########-]  89% (49\/55) MACCS\r[#########-]  91% (50\/55) MACCS\r[#########-]  93% (51\/55) MACCS\r[#########-]  95% (52\/55) MACCS\r[##########]  96% (53\/55) MACCS\r[##########]  98% (54\/55) MACCS\r[##########] 100% (55\/55) MACCS\n","truncated":false}}
%---
%[output:933d7d44]
%   data: {"dataType":"text","outputData":{"text":"[15:44:28][INFO]  MACCS 上位ヒット: IMIPRAMINE  (T = 1.0000)\n","truncated":false}}
%---
%[output:911e341b]
%   data: {"dataType":"text","outputData":{"text":"[15:44:28][INFO]  一致: ECFP4 と MACCS の両方が同一化合物を同定。\n","truncated":false}}
%---
%[output:812a31b9]
%   data: {"dataType":"text","outputData":{"text":"[15:44:28][INFO]    最終結論: IMIPRAMINE -- 確定。\n","truncated":false}}
%---
