%[text] # S04: バーチャルスクリーニング — 類似度で薬物候補を探す
%[text] EasyMolKit 応用ストーリー — レイヤー 2
%[text] 
%[text] 既存薬と「似た形」の未知化合物を数百件のライブラリから自動で拾い出せるでしょうか？
%[text] 創薬研究では候補化合物が膨大にのぼり、ひとつひとつ実験で試すのは現実的ではありません。**バーチャルスクリーニング**はフィンガープリントと類似度スコアを使って、実験の前に有望な候補を計算機で絞り込む手法です。
%[text] このスクリプトでは、イブプロフェン（COX 阻害薬）をクエリとして FDA 承認薬 200 件をランク付けし、リガンドベース類似度スクリーニングの全ワークフローを体験します。
%[text] ## 学習目標
%[text] - リガンドベース バーチャルスクリーニング（LBVS）のワークフローを理解する
%[text] - `emk.similarity.rankBy` で効率的なバルクタニモト検索を行う
%[text] - タニモトスコア閾値を実際の創薬文脈で解釈する
%[text] - フィンガープリント種類（ECFP4 vs MACCS）がヒットリストに与える影響を理解する \
%[text] ## 前提条件
%[text] - F03（フィンガープリント）と F04（類似度）の完了
%[text] - 推奨: S01（カフェインの仲間）と S02（リピンスキーフィルタ）
%[text] - RDKitインストール済み（`emk.setup.install()` を一度だけ実行しておく）
%[text] - 追加Toolbox不要（MATLAB だけで動きます） \
%[text] **所要時間**: 25〜40 分 | 実行方法: Ctrl+Enter でセクションを一つずつ実行
%[text] **データ**
%[text] - `data/list/fda_drugs.csv` — FDA 承認薬 200 種（ChEMBL、CC-BY-SA 3.0）
%[text] - 列: ChEMBLID, Name, SMILES, MolecularWeight, ALogP, HBondDonors, HBondAcceptors, TPSA, RotatableBonds, Source \
%[text] **参考文献**
%[text] - Willett P (2006) Similarity-based virtual screening using 2D fingerprints. *Drug Discov Today* 11:1046-1053. 〔要機関アクセス〕
%[text] - Rogers D & Hahn M (2010) Extended-connectivity fingerprints. *J Chem Inf Model* 50:742-754. 〔要機関アクセス〕
%[text] - Johnson MA & Maggiora GM (1990) *Concepts and Applications of Molecular Similarity*. Wiley. 〔書籍〕\
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
emk.setup.initPython();
%[text] Python/RDKit プロセスをウォームアップします（最初の呼び出しは少し時間がかかります）。
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;
logInfo("S04: セットアップ完了");
%%
%[text] ## セクション 1: ヒット化合物
%[text] ### リガンドベース バーチャルスクリーニング（LBVS）
%[text]   Johnson & Maggiora（1990）の分子類似性原理:
%[text]     「構造的に類似した分子は類似した生物活性を持つ傾向がある。」
%[text]  LBVS はこの原理を活用する:
%[text] \[クエリ(ヒット化合物)\] → \[フィンガープリント計算\] → \[ライブラリとの類似度算出\] → \[順位付けされた候補リスト\] → \[実験による検証(フォローアップ)\]
%[text]   構造ベースドッキングに対する利点:
%[text] - タンパク質構造が不要
%[text] - 数百万化合物をスクリーニングできる速さ
%[text] - 化学シリーズ内のスキャフォールドホッピングに有効 \
%[text]   限界:
%[text] - クエリと大きく異なる構造的に新規な（「スキャフォールドホッピング」）活性化合物を発見できない。 \
%[text] ヒット化合物: イブプロフェン
%[text]   1961 年に Boots Laboratories（英国）の Stewart Adams が発見。  COX-1 と COX-2 酵素を阻害し、プロスタグランジン合成をブロック。 世界で最も多く消費される薬の 1 つ（OTC 鎮痛薬/解熱薬）。アリールプロピオン酸（プロフェン）クラスの NSAID に属する。
HIT_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O";
HIT_NAME   = "Ibuprofen";

mol_hit = emk.mol.fromSmiles(HIT_SMILES);
logInfo("ヒット化合物: %s", HIT_NAME);
logInfo("  SMILES    : %s", HIT_SMILES);
logInfo("  重原子数  : %d", double(mol_hit.GetNumHeavyAtoms()));

desc_hit = emk.descriptor.calculate(mol_hit, ...
    ["MolWt", "LogP", "TPSA", "NumHDonors", "NumHAcceptors", "RingCount"]);

logInfo("主要性質:");
logInfo("  MW   : %.1f Da   (Ro5 上限: 500 -- PASS)", desc_hit.MolWt);
logInfo("  LogP : %.2f     (Ro5 上限: 5   -- PASS)", desc_hit.LogP);
logInfo("  TPSA : %.1f A^2  (経口上限: 140 -- PASS)", desc_hit.TPSA);
logInfo("  HBD  : %d        (Ro5 上限: 5   -- PASS)", desc_hit.NumHDonors);
logInfo("  HBA  : %d        (Ro5 上限: 10  -- PASS)", desc_hit.NumHAcceptors);
logInfo("  環数 : %d        (ベンゼン環 1 個; アリールプロピオン酸クラス)", desc_hit.RingCount);

figure("Name", "イブプロフェン（ヒット化合物）", "Position", [100 100 440 380]);
emk.viz.draw2d(mol_hit, Title="イブプロフェン（ヒット化合物）");
%[text] **✏️ やってみよう 1 — Ro5 と LogP の意味を確認しましょう**
%[text] イブプロフェンにはキラル中心が 1 つあります（カルボキシル基に隣接するアルファ炭素）。
%[text] 薬理学的に活性なのは (S)-エナンチオマー（COX 阻害薬）のみですが、
%[text] (R)-体は体内でイソメラーゼにより (S)-体に変換されるため、
%[text] ラセミ（R/S）混合物が OTC 薬として販売されています。
%[text] 
%[text] **Q1**: イブプロフェンはリピンスキーのルール・オブ・ファイブに合格しますか？
%[text] 1 行のテーブルを作り `emk.filter.lipinski()` を呼んで確認してみましょう。
%[text]  `d = emk.descriptor.calculate(mol_hit, ["MolWt","LogP","NumHDonors","NumHAcceptors"])`
%[text]  `t = struct2table(d)`
%[text]  `emk.filter.lipinski(t)`
%[text] 期待値: 違反 0、Pass\_Ro5 = true。
%[text] 
%[text] **Q2**: 計算 LogP（中性形）と血液中の「見かけの LogP」はなぜ異なるのでしょう？ 
%[text] ヒント: LogP は中性種で測定されますが、イブプロフェン（pKa ~4.4）は 血液 pH 7.4 で約 99% がイオン化しています。 D（分布係数）はイオン化を考慮するため D \<\< LogP になります。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: 化合物ライブラリの読み込み
%[text] ### 化合物ライブラリとは？
%[text] 創薬研究において「化合物ライブラリ」とは、スクリーニング（薬の候補の探索）に利用可能な分子のコレクション（集合体）のことです。目的やステージに応じて、ライブラリの規模や性質は大きく異なります。 
%[text] - フォーカスライブラリ（約100 〜 10,000化合物）: 特定の標的タンパク質や、特定の化学骨格に狙いを絞って化学者が設計したアナログ（類似体）の集まり。  \
%[text] - 企業ライブラリ（約10,000 〜 1,000,000化合物）: 製薬企業などが過去の製薬プロジェクトを通じて独自に蓄積してきた膨大な資産。  \
%[text] - 商業ライブラリ（約10億 〜 100億化合物）: 化学品サプライヤーが提供する超巨大な化合物空間（例: \*Enamine REAL Space\* など）。 \
%[text] 計算機上でデザインされ、注文に応じて合成可能なものも含まれます。
%[text] 
%[text] 本チュートリアルでのアプローチ：薬物再利用（リポジショニング）このスクリプトでは、巨大なライブラリの代わりに 200種のFDA承認薬を「プロキシ（代替）ライブラリ」として使用します。 既存の承認薬を対象にスクリーニングを行う手法は、「薬物再利用（ドラッグ・リポジショニング）」戦略と呼ばれます。すでに人での安全性が確認されているため、初期の安全性試験をバイパスでき、開発期間を大幅に短縮（通常10〜15年かかる創薬を3〜5年に短縮）できる強力なアプローチです。 
%[text] 
%[text] 代表的なリポジショニングの成功例:
%[text] - **サリドマイド**: 1950年代の「鎮静薬・睡眠薬」 → 2006年に「多発性骨髄腫（血液がん）」の治療薬として再承認。
%[text] - **シルデナフィル**: もともとは「狭心症」の治療薬 → 副作用の機序を逆手に取り、1998年に「勃起不全（バイアグラ）」として大ヒット。
%[text] - **メトホルミン**: 定番の「糖尿病」治療薬 → 近年、抗がん作用や抗老化（長寿研究）への転用が世界中で研究中。
LIBRARY_FILE = fullfile(projectRoot, "data", "list", "fda_drugs.csv");
library = readtable(LIBRARY_FILE, "TextType", "string");

logInfo("ライブラリ読み込み完了: %d 化合物", height(library));
%[text] ライブラリのサマリー
logInfo("ライブラリ記述子サマリー:");
logInfo("  MW    -- 最小: %5.1f  中央値: %5.1f  最大: %6.1f  Da", ...
    min(library.MolecularWeight), median(library.MolecularWeight), max(library.MolecularWeight));
logInfo("  ALogP -- 最小: %5.2f  中央値: %5.2f  最大: %5.2f", ...
    min(library.ALogP), median(library.ALogP), max(library.ALogP));
%[text] ヒット化合物がライブラリ内にあるか確認
hitInLib = strcmpi(library.Name, HIT_NAME);
if any(hitInLib)
    logInfo("  注: %s はライブラリ内に存在する（上位に自己マッチ T=1.0 が現れる）", ...
        HIT_NAME);
else
    logInfo("  注: %s はライブラリ内にない -- 全ヒットは真のアナログ", HIT_NAME);
end
%[text] **✏️ やってみよう 2 — ライブラリの中身を調べましょう**
%[text] **Q1**: ライブラリ内のユニークな薬物名は何件ありますか？
%[text] ヒント: `numel(unique(library.Name))`
%[text] 
%[text] **Q2**: 最も重い化合物（MolecularWeight 最大）は何でしょう？
%[text] ヒント: `library(library.MolecularWeight == max(library.MolecularWeight), :)`
%[text] 期待値: 大型の天然物由来薬（マクロライドまたはペプチド系など）です。
%[text] 
%[text] **Q3**: ライブラリの ALogP（Ghose-Crippen）と `emk.descriptor.calculate` の LogP（Wildman-Crippen）は
%[text] どちらもオクタノール-水分配を推定します。イブプロフェンで 2 つの値は近いですか？
%[text] ヒント: `desc_hit.LogP` vs ライブラリ内のイブプロフェンの ALogP（存在する場合）。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: ライブラリの ECFP4 フィンガープリントを計算する
%[text] ### **ECFP4（Morgan フィンガープリント）の仕組み**
%[text] 計算機は分子の「絵」をそのまま理解できません。そこで分子の構造を 0 と 1 のデジタル符号（ビットベクトル）に変換したものがフィンガープリントです。 
%[text] その代表格である ECFP4（拡張連結性フィンガープリント） は、分子の中の各原子を出発点として、周囲の化学的な環境を「結合の半径（ホップ数）」ごとに探索し、記憶していきます（Morgan アルゴリズム）。 
%[text] - 半径 0: 注目した原子そのものの情報（元素の種類、電荷、芳香族性など）
%[text] - 半径 1: 注目した原子 ＋ 直接隣接する原子（1結合先）の環境
%[text] - 半径 2: 注目した原子 ＋ 2ホップ先までの環境（これが ECFP4 に相当。直径4結合分をカバー）  \
%[text] これらの探索されたすべての局所環境が、数式的に 2048ビットの固定長ベクトル にハッシュ（圧縮）されます。 多くの類似した局所環境を共有している 2 つの分子は、フィンガープリント上でも「同じ位置のビットが ON（1）になる」確率が高くなります。その結果、後述するタニモトスコアが高く算出される仕組みです。 
%[text] ECFP4 の主な用途
%[text] 創薬の現場において、ECFP4 は以下の用途で世界中で標準的に使われています。 
%[text] 1. 類似度検索: クエリ分子に似た構造の化合物を一瞬で見つけ出す（この演習のメインテーマ）。
%[text] 2. 化合物クラスタリング: 膨大なライブラリを構造の似たグループごとに分類し、多様性を評価する。
%[text] 3. 機械学習（AI創薬）: 分子の特徴量を表すインプットデータ（記述子）としてAIモデルに入力する。 \
logInfo("ライブラリ %d 化合物の ECFP4 フィンガープリントを計算中...", height(library));

lib_fps   = cell(1, height(library));
lib_valid = true(1, height(library));

for i = 1:height(library)
    smi = library.SMILES(i);
    if ~emk.mol.isValid(smi)
        logWarn("  スキップ 行 %d (%s): 無効な SMILES", i, library.Name(i));
        lib_valid(i) = false;
        continue;
    end
    mol_lib    = emk.mol.fromSmiles(smi);
    lib_fps{i} = emk.fingerprint.morgan(mol_lib);
    logProgress(i, height(library), "FP");
end

nLibValid = sum(lib_valid);
logInfo("フィンガープリント準備完了: %d / %d 化合物", nLibValid, height(library));
%[text] 有効なサブセットを抽出
lib_fps_valid    = lib_fps(lib_valid);
lib_names_valid  = library.Name(lib_valid);
lib_smiles_valid = library.SMILES(lib_valid);
lib_mw_valid     = library.MolecularWeight(lib_valid);
lib_alogp_valid  = library.ALogP(lib_valid);
%[text] クエリフィンガープリントを計算
fp_hit = emk.fingerprint.morgan(mol_hit);
nOnBits = sum(emk.fingerprint.toArray(fp_hit));
logInfo("イブプロフェン ECFP4: %d ON ビット / 2048 合計（密度 %.1f%%）", ...
    nOnBits, 100*nOnBits/2048);
%[text] **✏️ やってみよう 3 — フィンガープリント密度を比べましょう**
%[text] イブプロフェンは単環・重原子 13 個、カフェイン（S01）は縮合二環・重原子 14 個です。ECFP4 の ON ビット数はどちらが多いでしょう？
%[text] **Q1**: カフェインのフィンガープリントを計算して比べてみましょう。
%[text]  `mol_caf = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")`
%[text]  `fp_caf  = emk.fingerprint.morgan(mol_caf)`
%[text]  `sum(emk.fingerprint.toArray(fp_caf))`
%[text] 期待値: カフェインの方が ON ビット数が多いです。縮合二環系が多くの異なる局所環境を作り出すためです。
%[text] 
%[text] **Q2**: 半径を 2 から 3 に増やすとどうなりますか？
%[text]  `fp_r3 = emk.fingerprint.morgan(mol_hit, Radius=3)`
%[text]  `sum(emk.fingerprint.toArray(fp_r3))`
%[text] 期待値: 通常は ON ビット数が同じか少なくなります。半径 3 では各隣接域が
%[text] 大きくなりよりユニークになるため衝突が減りますが、
%[text] 非常に似た化合物でも共有ビットが失われることがあります。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: バーチャルスクリーンを実行する
%[text] ### バルクタニモト類似度検索とは？
%[text] 2つの分子（ビットベクトル）がどれくらい似ているかを数値化する指標として、最も広く使われているのがタニモト類似度（Tanimoto Similarity / Jaccard 指数）です。 
%[text] $&dollar&;&dollar&;T(A, B) = \\frac{|A \\cap B|}{|A \\cup B|} = \\frac{\\text{分子AとBの両方で ON になっているビット数}}{\\text{分子AまたはBの少なくとも一方で ON になっているビット数}}&dollar&;&dollar&; $
%[text] $&dollar&;T = 1.0&dollar&;$: フィンガープリントが完全一致（同一化合物、あるいは ECFP4 で区別できない極めて近いスキャフォールド）。
%[text] $&dollar&;T = 0.0&dollar&;$: 共通するビットがゼロ（完全に異なる構造）。
%[text] ### スコア（閾値）の現実的な解釈
%[text] 創薬化学において、ECFP4 のタニモトスコアは経験的に以下のように解釈されます（Willett 2006）。ヒットリストを眺める際の重要なガイドラインです。 
%[text] - $&dollar&;T \\ge 0.85&dollar&; $（非常に類似）: ほぼ同じ骨格（スキャフォールド）を持ち、置換基がごくわずかに変化しただけのアナログ。
%[text] - $&dollar&;T \\ge 0.65&dollar&;$ （類似）: 同じ化学クラス、または共通のコア環系を共有している。
%[text] - $&dollar&;T \\ge 0.40&dollar&; $（中程度）: 共通の官能基はあるが、骨格自体は異なっている。
%[text] - $&dollar&;T \< 0.40&dollar&;$ （低 / 非類似）: 構造的な共通点はほぼ見当たらない。 \
%[text] ### 計算の効率化: `emk.similarity.rankBy`
%[text] 本スクリプトで使用する `emk.similarity.rankBy` 関数は、裏側で RDKit の `BulkTanimotoSimilarity` を呼び出しています。 これは、200件の化合物を1件ずつバラバラに Python 側に投げて計算するのではなく、ライブラリ全体を1回のバッチ（一括）処理でスクリーニングします。MATLAB と Python 間の通信（IPC）の往復回数を 1 回に抑えるため、何万件もの化合物に対しても極めて高速に動作します。
tic;
vs_result = emk.similarity.rankBy(fp_hit, lib_fps_valid);
tElapsed = toc;

logInfo("バーチャルスクリーン完了: %d 化合物を %.3f 秒で順位付け", ...
    numel(vs_result.Scores), tElapsed);

%[text] 上位 10 件のヒットを解釈付きで表示
TOP_DISPLAY = 10;
logInfo("--- %s の上位 %d 件ヒット ---", HIT_NAME, TOP_DISPLAY);
logInfo("%-5s  %-28s  %9s  %s", "順位", "名前", "タニモト", "解釈");
for k = 1:min(TOP_DISPLAY, numel(vs_result.Scores))
    idx   = vs_result.Indices(k);
    score = vs_result.Scores(k);
    name  = lib_names_valid(idx);

    if score >= 1.0 - 1e-6
        interp = "<-- 自己マッチ (T=1.0)";
    elseif score >= 0.85
        interp = "非常に類似";
    elseif score >= 0.65
        interp = "類似";
    elseif score >= 0.40
        interp = "中程度";
    else
        interp = "低";
    end

    logInfo("  %2d.  %-28s  %9.4f  %s", k, name, score, interp);
end

%[text] **✏️ やってみよう 4 — 他の NSAID を探しましょう**
%[text] 上位ヒットには他の NSAID（フルルビプロフェン、ケトプロフェン）が含まれますか？
%[text] これらはイブプロフェンと同じアリールプロピオン酸コア（アレーン—CH(CH₃)—COOH）を持ちます。
%[text] 
%[text] **Q1**: 自己マッチ以外の上位 2 件を確認しましょう。
%[text] 期待値: FLURBIPROFEN（T~0.40）、KETOPROFEN（T~0.39）あたりです。T~0.40 にとどまる理由は、フッ素置換基や余分な環が局所環境を変えて ECFP4 が区別するためです。（ECFP4 で T \> 0.30 は薬物様アナログとして依然意味があります。0.65/0.85 はガイドラインです。）
%[text] **注目**: 4 位に AMPHETAMINE SULFATE（T~0.35）が現れる場合があります。
%[text] フェニル環＋短いアルキル鎖の局所環境がイブプロフェンの部分構造と重なるためですが、カルボン酸をもたず CNS 刺激薬であり抗炎症活性はありません。これは ECFP4 の「構造的偽陽性」の典型例です。スクリーニング結果は必ず化学者が構造を確認してから実験に進む必要がある理由がここにあります。
%[text] 
%[text] **Q2**: インドメタシン（酢酸系 NSAID）は何位でしょう？
%[text] 期待値: プロフェンクラスよりタニモトが低くなります。これは ECFP4 がスキャフォールドの差を捕えるためです。
%[text] 
%[text] **Q3**: アスピリン（`"CC(=O)Oc1ccccc1C(=O)O"`）をクエリにするとヒットリストはどう変わりますか？
%[text]  `mol_asp = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O")`
%[text]  `fp_asp  = emk.fingerprint.morgan(mol_asp)`
%[text]  `res_asp = emk.similarity.rankBy(fp_asp, lib_fps_valid)`
%[text]  `lib_names_valid(res_asp.Indices(1:5))`
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: ヒットレポートを作成する
%[text] 構造化結果テーブルを作成: 順位、名前、タニモト、MW、ALogP。これが実験チームに渡すレポートとなる。
REPORT_TOP = 15;
rankVec  = (1:REPORT_TOP)';
nameVec  = strings(REPORT_TOP, 1);
scoreVec = zeros(REPORT_TOP, 1);
mwVec    = zeros(REPORT_TOP, 1);
logpVec  = zeros(REPORT_TOP, 1);
smiVec   = strings(REPORT_TOP, 1);

for k = 1:REPORT_TOP
    if k > numel(vs_result.Scores); break; end
    idx         = vs_result.Indices(k);
    nameVec(k)  = lib_names_valid(idx);
    scoreVec(k) = vs_result.Scores(k);
    mwVec(k)    = lib_mw_valid(idx);
    logpVec(k)  = lib_alogp_valid(idx);
    smiVec(k)   = lib_smiles_valid(idx);
end

hits_tbl = table(rankVec, nameVec, scoreVec, mwVec, logpVec, smiVec, ...
    'VariableNames', ["Rank", "Name", "Tanimoto", "MW_Da", "ALogP", "SMILES"]);

logInfo("上位 %d 件ヒットレポート:", REPORT_TOP);
disp(hits_tbl(:, 1:5));   % SMILES 列なしで表示

%[text] タニモト閾値サマリー
logInfo("閾値別内訳（上位 %d 件）:", REPORT_TOP);
logInfo("  T >= 0.85 (非常に類似): %d", sum(scoreVec >= 0.85));
logInfo("  T >= 0.65 (類似)      : %d", sum(scoreVec >= 0.65 & scoreVec < 0.85));
logInfo("  T >= 0.40 (中程度)    : %d", sum(scoreVec >= 0.40 & scoreVec < 0.65));
logInfo("  T <  0.40 (低)        : %d", sum(scoreVec < 0.40 & scoreVec > 0));

%[text] **✏️ やってみよう 5 — ヒットレポートを分析しましょう**
%[text] **Q1**: `hits_tbl` を MW で並び替えて、上位ヒットがイブプロフェン（206 Da）付近に
%[text] クラスターしているか確認しましょう。
%[text] ヒント: `sortrows(hits_tbl, "MW_Da")`
%[text] 
%[text] **Q2**: ALogP が最も高いヒットはどれですか？それでも薬物様（ALogP \< 5）ですか？
%[text] ALogP が高いと水への溶解性が低くなりやすく、活性があっても製剤上の課題となります。
%[text] 
%[text] **Q3**: ヒットテーブルに Ro5 合格列を追加できますか？
%[text] テーブル列名を変更して `emk.filter.lipinski()` を呼びましょう。
%[text] （MW, LogP, NumHDonors, NumHAcceptors が必要。CSV に HBD/HBA があります。）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: スコア分布
%[text] ### 化学空間カバレッジ（ライブラリの多様性評価）
%[text] 上位のトップヒットだけでなく、ライブラリに含まれる全化合物のタニモトスコアの分布（ヒストグラム）を俯瞰することは非常に重要です。これにより、用意したライブラリに対して、クエリ（探索の起点となる分子）がどのような立ち位置にあるのか（化学空間でのカバレッジ）が見えてきます。 
%[text] 
%[text] $&dollar&;T = 0&dollar&; $付近に鋭いピークが立ち、高スコア側に裾野が全くない場合 クエリ分子がそのライブラリの中で「孤立した外れ値」であることを意味します。類似化合物がライブラリ内に存在しないため、スクリーニングによる展開（フォローアップ）が難しいケースです。 
%[text] 高い $&dollar&;T&dollar&; $側（右側）に向かって分布が緩やかに広がっている場合  ライブラリ内にクエリと親戚関係にあるような「近似アナログ」が豊富に含まれていることを示します。構造を少しずつ変えた際の活性変化（SAR: 構造活性相関）を追うのに適した、密度の高い空間と言えます。
allScores = vs_result.Scores;

figure("Name", "バーチャルスクリーンスコア分布", "Color", "white", "Position", [100 100 580 420]);
histogram(allScores, 30, "FaceColor", [0.2 0.55 0.85], ...
    "EdgeColor", "white", "FaceAlpha", 0.8, ...
    "DisplayName", "タニモトスコア分布");
hold on;
%[text] 主要閾値をマーク
xline(0.85, "r--", "LineWidth", 1.5, "DisplayName", "非常に類似 (0.85)");
xline(0.65, "m--", "LineWidth", 1.5, "DisplayName", "類似 (0.65)");
xline(0.40, "g--", "LineWidth", 1.5, "DisplayName", "中程度 (0.40)");

hold off;
xlabel("タニモトスコア (ECFP4)");
ylabel("ライブラリ化合物数");
title(sprintf("類似度分布: %s vs FDA 薬物 200 種", HIT_NAME));
legend("Location", "northeast");
grid("on");
%[text] ゾーン別件数
nVSim  = sum(allScores >= 0.85);
nSim   = sum(allScores >= 0.65 & allScores < 0.85);
nMod   = sum(allScores >= 0.40 & allScores < 0.65);
nLow   = sum(allScores < 0.40);

logInfo("%d 化合物全体のスコア分布:", nLibValid);
logInfo("  非常に類似 (T>=0.85): %3d  (%.0f%%)", nVSim, 100*nVSim/nLibValid);
logInfo("  類似       (T>=0.65): %3d  (%.0f%%)", nSim,  100*nSim /nLibValid);
logInfo("  中程度     (T>=0.40): %3d  (%.0f%%)", nMod,  100*nMod /nLibValid);
logInfo("  低         (T< 0.40): %3d  (%.0f%%)", nLow,  100*nLow /nLibValid);
%[text] **✏️ やってみよう 6 — スコア分布を比べましょう**
%[text] **Q1**: 200 件の薬物ライブラリのうち、イブプロフェンとの T \>= 0.65 の割合はどのくらいですか？
%[text] イブプロフェンはこのライブラリでよく代表されていますか（近傍が多い）？
%[text] それとも外れ値ですか（近傍が少ない）？
%[text] 
%[text] **Q2**: クエリをカフェインに変えて 2 つのヒストグラムを比較しましょう。
%[text] どちらが高タニモトのヒットをより多く持ちますか？
%[text]  `mol_caf = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")`
%[text]  `fp_caf  = emk.fingerprint.morgan(mol_caf)`
%[text]  `res_caf = emk.similarity.rankBy(fp_caf, lib_fps_valid)`
%[text]  `figure; histogram(res_caf.Scores, 30); title("カフェイン vs FDA 薬物")`
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 7: 上位 3 件のヒットを可視化する
%[text] 自己マッチ以外の上位 3 件を並べて描き、構造的類似度を確認する。
nonSelf = find(vs_result.Scores < 1.0 - 1e-6);
drawK   = nonSelf(1:min(3, numel(nonSelf)));

logInfo("自己マッチ以外の上位 %d 件を描画中:", numel(drawK));
figure("Name", "上位ヒット", "Position", [100 100 960 340]);
for j = 1:numel(drawK)
    k     = drawK(j);
    idx   = vs_result.Indices(k);
    name  = lib_names_valid(idx);
    smi   = lib_smiles_valid(idx);
    score = vs_result.Scores(k);
    mol_j = emk.mol.fromSmiles(smi);
    subplot(1, numel(drawK), j);
    emk.viz.draw2d(mol_j, Title=sprintf("順位 %d: %s (T=%.3f)", k, name, score));
    logInfo("  順位 %d: %s (T=%.4f, MW=%.1f, ALogP=%.2f)", ...
        k, name, score, lib_mw_valid(idx), lib_alogp_valid(idx));
end

%[text] **✏️ やってみよう 7 — 構造を視覚的に比べましょう**
%[text] **Q1**: 上位 3 件のヒットはイブプロフェンに視覚的に似ていますか？
%[text] 共通する構造的特徴は何でしょう？
%[text] 期待値: アリール環 + プロピオン酸側鎖（「プロフェン」スキャフォールド）が共有されています。
%[text] 
%[text] **Q2**: イブプロフェンと上位ヒットを 1×2 の figure に並べて描いてみましょう。
%[text]  `figure("Position", [100 100 750 300])`
%[text]  `subplot(1,2,1); emk.viz.draw2d(mol_hit, Title="イブプロフェン（クエリ）")`
%[text]  `subplot(1,2,2); emk.viz.draw2d(emk.mol.fromSmiles(smiVec(1)), Title=nameVec(1))`
%[text] （注: `draw2d` はデフォルトで自前の figure を開きます。subplot を使うと
%[text] 2 つの構造を並べて比較できます。）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 8: ECFP4 vs MACCS キー -- フィンガープリントの選択は重要か？
%[text] 
%[text] ### フィンガープリントの選択がヒットリストに与える影響
%[text] 分子をデジタル化する「フィンガープリント」にはいくつかの種類があり、どれを選ぶかによってスクリーニングの結果（ヒットリスト）は大きく変わります。ここでは、今回比較した2つの特性を整理します。 
%[text] - ECFP4（2048 ビット / 円構造記述子） 各原子の周囲にあるローカルな原子環境（つながり方）を詳細に符号化します。
%[text] - MACCS キー（166 ビット / フラグメントベース記述子） 創薬化学の専門家が定義した「特定の構造フラグメント（部分構造）が存在するかどうか」を 0 と 1 で符号化します。   例: 「分子内にカルボン酸はあるか？」「芳香環はあるか？」「ハロゲン原子はあるか？」など \
%[text] ### それぞれの特徴と使い分け
%[text] - ECFP4: 分子の微細な構造の違い（置換基の位置や環のサイズ変化など）に対して非常に敏感です。
%[text] - MACCS キー: 構造の細部（装飾）の変化には鈍感な反面、大まかな「化学クラス（骨格の共通性）」を広く捉える\*\*のに適しています。  \
%[text] ### なぜ MACCS キーを使うとタニモトスコアが高くなりやすいのか？
%[text] 実際に計算してみると、MACCS の方が全体的にスコアが高く（インフレして）出やすくなります。理由は以下の2点です。 
%[text] 1. ビット数の圧倒的な差（166 ビット $&dollar&;\\ll&dollar&;$ 2048 ビット）: ビット長が短い（解像度が粗い）ため、構造の「不一致（ミスマッチ）」がカウントされにくいため。
%[text] 2. 基本フラグメントの共有: クエリであるイブプロフェンと、ライブラリ内の他の非ステロイド性抗炎症薬（NSAIDs）は、「芳香環がある」「カルボン酸がある」といった MACCS 定義の主要なフラグメントの多くを共通して持っているため。  \
%[text] ### スクリーニング結果の解釈
%[text] - ECFP4 と MACCS キーの両方で上位にヒットした場合: 大まかな化学クラスだけでなく、細部の局所環境までよく似ていることを意味し、「強力な構造類似性（親戚関係）」が裏付けられます。
%[text] - MACCS では上位だが、ECFP4 では下位に落ちる場合: 「カルボン酸を持つ芳香族化合物」という大枠のクラス（性質）は共有しているものの、「詳細な骨格や置換基の付き方はまったく異なる」ことを示しています。 \
logInfo("ライブラリの MACCS フィンガープリントを計算中...");

maccs_fp_hit = emk.fingerprint.maccs(mol_hit);
maccs_lib_fps = cell(1, numel(lib_fps_valid));

for i = 1:numel(lib_fps_valid)
    mol_i = emk.mol.fromSmiles(lib_smiles_valid(i));
    maccs_lib_fps{i} = emk.fingerprint.maccs(mol_i);
    logProgress(i, numel(lib_fps_valid), "MACCS");
end

maccs_result = emk.similarity.rankBy(maccs_fp_hit, maccs_lib_fps);

%[text] 上位 10 件を ECFP4 と MACCS で並べて比較
COMPARE_N = 10;
logInfo("--- ECFP4 vs MACCS キー: 上位 %d 件比較 ---", COMPARE_N);
logInfo("%-5s  %-24s  %9s  %-24s  %9s", ...
    "順位", "ECFP4 ヒット", "T(ECFP4)", "MACCS ヒット", "T(MACCS)");
for k = 1:COMPARE_N
    ie = vs_result.Indices(k);
    im = maccs_result.Indices(k);
    logInfo("  %2d.  %-24s  %9.4f  %-24s  %9.4f", k, ...
        lib_names_valid(ie), vs_result.Scores(k), ...
        lib_names_valid(im), maccs_result.Scores(k));
end

%[text] 順位相関: 上位 20 件で 2 種のフィンガープリントの一致度
TOP_CORR = 20;
ecfp4_top = arrayfun(@(k) lib_names_valid(vs_result.Indices(k)),    1:TOP_CORR);
maccs_top = arrayfun(@(k) lib_names_valid(maccs_result.Indices(k)), 1:TOP_CORR);
nOverlap  = numel(intersect(ecfp4_top, maccs_top));
logInfo("ECFP4 と MACCS の上位 %d 件重複: %d / %d (%.0f%%)", ...
    TOP_CORR, nOverlap, TOP_CORR, 100*nOverlap/TOP_CORR);

%[text] **✏️ やってみよう 8 — フィンガープリントの選択を比較しましょう**
%[text] **Q1**: ECFP4 と MACCS は同じ上位 3 件の化合物を返しますか？
%[text] 一致しない場合、何がその違いを生み出しているのでしょう？
%[text] （MACCS キーには「COOH がある」などの汎用フラグメントが含まれます。
%[text] どんなカルボン酸でも MACCS スコアは高くなりますが、ECFP4 はその COOH が
%[text] 類似の環コンテキストにある場合のみ高いスコアを返します。）
%[text] 
%[text] **Q2**: 上位ヒットの MACCS タニモトスコアは ECFP4 より高い傾向があります。
%[text] これは MACCS が「より優れている」ことを意味しますか？
%[text] （ヒント: スコアが高い理由はビット数が少なく不一致になりにくいためで、精度の高さではありません。
%[text] リード最適化にはスキャフォールドの違いに敏感な ECFP4 が好まれます。）
%[text] 
%[text] **Q3**: タニモトの代わりにダイス指標で類似度を計算してみましょう。
%[text]  `res_dice = emk.similarity.rankBy(fp_hit, lib_fps_valid, Metric="dice")`
%[text] Dice = 2|A AND B| / (|A| + |B|)。ビットベクトルでは Dice \>= Tanimoto が常に成立します。
%[text] 順位はどう変わりますか？
% ... （ここにコードを書いてみましょう）

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
