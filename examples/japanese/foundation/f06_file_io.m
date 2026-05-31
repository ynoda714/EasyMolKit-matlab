%[text] # F06: 分子ファイルの読み書き
%[text] EasyMolKit 基礎チュートリアル — レイヤー 1
%[text] 
%[text] 実際の研究プロジェクトや創薬の現場では、数千から数万件におよぶ大量の化合物を一括して取り扱うことになります。これらを1つずつ手入力していくのは不可能です。そこで重要になるのが、世界の主要な化学データベース（PubChem や ChEMBL など）や、外部のドッキングシミュレーションソフトとデータを直接やり取りするための「ファイル入出力（I/O）」の技術です。このチュートリアルでは、最もシンプルで軽量な「SMILES リストファイル」と、分子の3次元形状や様々な社内データも一緒に詰め込める業界標準の「SDF（Structure-Data File）」という2大フォーマットの読み書きをマスターします。さらに、これまでに学んだF01〜F05の知識（構造チェック、物性計算、部分構造フィルタリングなど）をすべて数珠つなぎに連結し、自動でデータを処理して結果をファイルに保存する「ミニデータパイプライン」の構築に挑戦してみましょう。
%[text] **学習目標**
%[text] - `emk.io.readSmilesList` で SMILES リストファイルを読み込む
%[text] - `emk.io.readSdf` で SDF ファイルを読み込む
%[text] - `emk.io.writeSdf` で分子を SDF ファイルに書き出す
%[text] - ファイル I/O と記述子計算・フィルタリングを組み合わせる \
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
emk.setup.initPython(); %[output:96fa9617]
%[text] Python/RDKit プロセスのウォームアップ
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;

runDir = makeRunDir();   % result/runs/<タイムスタンプ>/ を作成
logInfo("F06: セットアップ完了 -- 出力ディレクトリ: %s", runDir); %[output:8dd385bb]
%%
%[text] ## セクション 1: SMILES リストファイルを読み込む
%[text] ### **テキスト1行ずつのシンプルな分子名簿**
%[text] 分子のデータを最もシンプルに保存する方法が、1行ごとに「SMILES 文字列」と「分子の名前（ID）」をタブやスペースで区切って並べたテキストファイル（一般に `.smi` や `.txt` 拡張子）です。この形式はファイルサイズが非常に小さく、中身をメモ帳などで直接開いて人間の目で確認できるため、大量の分子のリストを大雑把に共有したいときにとても重宝されます。
%[text] EasyMolKit では、`emk.io.readSmilesList` 関数を呼び出すだけで、このファイルからすべての分子を一括で読み込むことができます。読み込まれたデータは、MATLAB ユーザーが最も扱いやすい「分子オブジェクトのセル配列（Cell array）」へと自動で変換されます。それでは、あらかじめ用意されている身近な化学物質のリスト（`everyday_chemicals.smi`）を読み込んで、MATLAB のワークスペースに分子たちが綺麗に並ぶ様子を確認してみましょう。
%[text] 
%[text] `emk.io.readSmilesList` は「1 行 1 SMILES」のテキストファイルを読み込みます。第 2 列（タブまたはスペース区切り）がある場合は名前/ラベルとして扱われます。
%[text] `#` で始まるコメント行と空行は自動でスキップされます。
%[text] 
%[text] 戻り値は RDKit Mol オブジェクトのセル配列です。無効な SMILES があっても警告を出してスキップするので、クラッシュしません。
%[text] 付属の `everyday_chemicals.csv` はカンマ区切りなので、このセクションでは練習用にプレーン SMILES リストを自前で書き出します。
smiles_file = fullfile(runDir, "sample_molecules.smi");
writelines([ ...
    "# F06 チュートリアル用サンプル SMILES リスト", ...
    "CCO            ethanol", ...
    "c1ccccc1       benzene", ...
    "CC(=O)Oc1ccccc1C(=O)O  aspirin", ...
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C  caffeine", ...
    "CC(=O)NC1=CC=C(C=C1)O  acetaminophen" ...
], smiles_file);

mols_from_smi = emk.io.readSmilesList(smiles_file); %[output:390d8d09]
logInfo("SMILES リストから %d 分子を読み込んだ", numel(mols_from_smi)); %[output:2c91ec2c]
%[text] 
%[text] **✏️ やってみよう 1 — 壊れた SMILES を混ぜたらどうなる？**
%[text] `sample_molecules.smi` に無効な SMILES 行（例: `"INVALID_SMILES  bad_entry"`）を
%[text] 追加して、このセクションをもう一度実行してみましょう。
%[text] 関数はクラッシュしますか？それとも警告を出してスキップしますか？
%[text] ヒント: ログ出力に `logWarn` のメッセージが現れるか確認してください。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: 分子を SDF に書き出す
%[text] ### **分子の「形」も「データ」も丸ごとパッキングするリッチな箱**
%[text] SMILES リストは手軽で便利ですが、「各原子の3次元的な位置座標」や「実験で得られた活性値（物性データ）」といった複雑な情報を一緒に保存することができません。そこで、世界の製薬・化学の現場で事実上の業界標準（デファクトスタンダード）として使われているのが「SDF（Structure-Data File）」というフォーマットです。SDF は、分子を構成する原子同士の結合関係（接続表）だけでなく、立体的な座標情報、さらにはユーザーが独自に追加した様々なテキストデータ（分子量や注記など）を1つのファイルの中に何千個分もパックにして格納できます。
%[text] EasyMolKit では `emk.io.readSdf` を使うことで、この立体情報や付加データを含んだリッチな分子たちを、一瞬で MATLAB のセル配列としてロードできます。ファイルの中にどんなユニークなデータ（プロパティ）が隠されているかも含めて、中身をスキャンしてみましょう。
%[text] 
%[text] ### コンセプト: SDF が業界標準として使われ続ける理由
%[text] SDF（Structure Data File）は 1980 年代に MDL Information Systems が導入した分子交換フォーマットです。長年にわたって業界標準として使われてきた理由は 3 つあります。
%[text] 1. **座標保持**: 2D/3D 原子座標を格納できます（SMILES 単独では不可）。ドッキングや 3D ファーマコフォア解析に必須です。
%[text] 2. **プロパティ付き**: `> <PROP_NAME>` フィールドで活性データや ID を構造と一緒に 1 ファイルに格納できます。
%[text] 3. **連結安全**:  区切りにより複数分子を安全に連結できます。 \
%[text] PubChem・ChEMBL・Reaxys はいずれも SDF でダウンロード提供しており、AutoDock・Glide 等のドッキングエンジンも SDF を入力として受け付けます。SDF の読み書きはケモインフォマティクスで最も実用的なファイル I/O スキルです。`emk.io.writeSdf` は Mol オブジェクトのセル配列をディスクに書き出します。
%[text] 出力パスは `makeRunDir()` で作成したディレクトリを使えば安全です。
sdf_out = fullfile(runDir, "sample_output.sdf");
emk.io.writeSdf(mols_from_smi, sdf_out); %[output:72cef0d6]
logInfo("%d 分子を書き出した: %s", numel(mols_from_smi), sdf_out); %[output:98ae2446]
%%
%[text] ## セクション 3: SDF を読み直す
%[text] ### **実験データを記録し、成果物を上書きから守る**
%[text] MATLAB 上で計算した新しい記述子（LogP など）や、フィルタリングの結果を外部へ出力したいときは、`emk.io.writeSdf` 関数を使って分子を新しい SDF ファイルとして書き出します。このとき実務で極めて重要になるのが、「過去の計算結果をうっかり上書きして消してしまわないこと（再現性の確保）」です。
%[text] これをスマートに解決するため、EasyMolKit には `makeRunDir()` という便利な補助関数が用意されています。この関数を実行すると、実行した瞬間の日付と時刻をもとに、たとえば `result/runs/20260523_001500/` といった「タイムスタンプ付きの専用フォルダ」を自動で作成してくれます。この綺麗な使い捨てフォルダの中に出力ファイル（SDF など）を保存するようにコードを組むことで、いつでも過去の実験データをクリーンに振り返ることができます。
%[text] `emk.io.readSdf` は SDF ファイルを読み込み、Mol オブジェクトのセル配列を返します。
mols_from_sdf = emk.io.readSdf(sdf_out); %[output:998b8a2d]
logInfo("SDF から %d 分子を再読み込みした", numel(mols_from_sdf)); %[output:91439fb7]
%[text] 
%[text] ラウンドトリップ後に正規 SMILES が保持されているか確認します。なお、入力した SMILES と表示が異なる場合があります（例: カフェイン）。
%[text] これは RDKit が同一分子を常に同じ「正規 SMILES」に変換するためで、正常な動作です。
logInfo("--- ラウンドトリップ SMILES 確認 ---"); %[output:82318383]
for i = 1:numel(mols_from_sdf) %[output:group:2778ade1]
    smi = emk.mol.toSmiles(mols_from_sdf{i});
    logInfo("  Mol %d: %s", i, smi); %[output:66806ecf]
end %[output:group:2778ade1]
%%
%[text] ## セクション 4: 付属データセットを使う
%[text] ### **バケツリレーで大量の分子を自動処理する**
%[text] これまでに学んできたすべてのテクニックを組み合わせ、プロのケモインフォマティシャンが日常的に行っている自動処理の仕組み（データパイプライン）を構築しましょう。私たちが目指すのは、**「①ファイルの読み込み ➔ ②不正な構造の自動排除 ➔ ③全員の物性値を一括計算 ➔ ④薬らしくない分子を足切り ➔ ⑤合格者だけを新しいファイルに保存」** という一連の流れを、ボタン一つで全自動で実行するバケツリレーのシステムです。
%[text] こうしたパイプラインを1度組んでおけば、中身のファイルが10個から1万個に増えたとしても、コンピューターが裏で黙々と正確に処理を進めてくれるようになります。MATLAB の得意とする行列・配列処理のパワーと RDKit の化学エンジンが融合する、EasyMolKit の真骨頂とも言えるワークフローを体験してみましょう。
%[text] EasyMolKit には `data/list/` 以下にキュレートされた分子リストが付属しています。
%[text] -  `everyday_chemicals.csv`  — 30 種の日用化学品（PubChem CC0）
%[text] -  `fda_drugs.csv`           — 200 種の FDA 承認薬（ChEMBL CC-BY-SA 3.0）
%[text] -  `pains.csv`               — PAINS SMARTS フィルタ（BSD-3）
%[text] -  `forensic_challenge.csv`  — 仮想法科学チャレンジセット（ラベル付き） \
%[text] ここでは `fda_drugs.csv` から FDA 承認薬を読み込み、Mol オブジェクトに変換します。
fda = readtable("data/list/fda_drugs.csv", TextType="string");
logInfo("FDA 薬物データセット: %d 化合物", height(fda)); %[output:4a25a224]
%[text] 
%[text] SMILES を解析します（200 分子あるので数秒かかる場合があります）。`isValid` で先に有効性を確認し、有効な分子だけ `fromSmiles` で変換します。
%[text] この順序にすることで、無効 SMILES があっても安全にスキップできます。
valid_mask = false(height(fda), 1);
fda_mols   = cell(height(fda), 1);
for i = 1:height(fda)
    if emk.mol.isValid(fda.SMILES(i))
        fda_mols{i}   = emk.mol.fromSmiles(fda.SMILES(i));
        valid_mask(i) = true;
    end
end
logInfo("有効な FDA 分子: %d / %d", sum(valid_mask), height(fda)); %[output:8597961b]
%%
%[text] ## セクション 5: サブセットをフィルタしてエクスポートする
%[text] リピンスキーの Ro5（Rule of Five）で薬物様分子をフィルタし、SDF に保存します。
fda_valid_mols = fda_mols(valid_mask);
fda_valid      = fda(valid_mask, :);

%[text] 記述子テーブルを計算します。
desc_tbl = emk.descriptor.batchCalculate(fda_valid_mols, ...
    ["MolWt", "LogP", "NumHDonors", "NumHAcceptors"]);
%[text] 
%[text] リピンスキーフィルタを適用します。
lipinski_tbl = emk.filter.lipinski(desc_tbl);   % Pass_Ro5 / Violations_Ro5 を追加 %[output:603355c3]

pass_mask = lipinski_tbl.Pass_Ro5;
logInfo("Ro5 通過: %d / %d FDA 分子", sum(pass_mask), numel(fda_valid_mols)); %[output:28b668da]
%[text] 
%[text] **✏️ やってみよう 2 — 緩和ルールで通過数はどう変わる？**
%[text] `MaxViolations=1` オプションで Ro5 を 1 項目まで超過できる緩和フィルタを試してみましょう。
%[text] `lipinski_tbl_relaxed = emk.filter.lipinski(desc_tbl, "MaxViolations", 1);`
%[text] 厳格な Ro5 と比べて何件多く通過しますか？
%[text] 新たに通過した薬物の名前も確認してみましょう。
%[text] （補足: MaxViolations=1 は「4 基準のうち 1 つまで超過を許容する」緩和 Ro5 です。
%[text]   Veber et al. 2002 の規則（TPSA ≤ 140 Å²、回転可能結合数 ≤ 10）とは別概念です。）
%[text] 
%[text] Ro5 を通過した分子を SDF に保存します。
passing_mols = fda_valid_mols(pass_mask);
sdf_filtered = fullfile(runDir, "fda_ro5_pass.sdf");
emk.io.writeSdf(passing_mols, sdf_filtered); %[output:4ec9371a]
logInfo("Ro5 通過分子を保存した: %s", sdf_filtered); %[output:3e55f303]
%[text] 
%[text] 記述子テーブルを CSV にも保存して、後で確認できるようにします。
csv_out = fullfile(runDir, "fda_ro5_pass_descriptors.csv");
writetable([fda_valid(pass_mask, "Name"), desc_tbl(pass_mask, :), ...
    lipinski_tbl(pass_mask, ["Violations_Ro5","Pass_Ro5"])], csv_out);
logInfo("記述子テーブルを保存した: %s", csv_out); %[output:93bcbda8]
%%
%[text] ## セクション 6: ミニパイプラインのまとめ
%[text] 今回のチュートリアルで学んだファイル I/O モジュールの関数一覧と、標準的なパイプライン設計のテンプレートです。
%[text] -  ステップ 1 -- ファイルから分子を読み込む  `emk.io.readSmilesList` / readSdf
%[text] -  ステップ 2 -- 有効性確認               `emk.mol.isValid`
%[text] -  ステップ 3 -- 記述子計算               `emk.descriptor.batchCalculate`
%[text] -  ステップ 4 -- フィルタ適用              `emk.filter.lipinski`
%[text] -  ステップ 5 -- フィンガープリント＋類似度  `emk.fingerprint.morgan`, `emk.similarity.*`（F04・F05 で学習済み。このスクリプトでは省略）
%[text] -  ステップ 6 -- 結果のエクスポート        `emk.io.writeSdf` / writetable \
%[text] ## セクション 7: まとめ
%[text] 
%[text:table]
%[text] | 関数 | 用途 |
%[text] | --- | --- |
%[text] | `emk.io.readSmilesList(path)` | SMILES リスト → Mol のセル配列 |
%[text] | `emk.io.readSdf(path)` | SDF ファイル → Mol のセル配列 |
%[text] | `emk.io.writeSdf(mols, path)` | Mol のセル配列 → SDF ファイル |
%[text] | `makeRunDir()` | タイムスタンプ付き出力ディレクトリを作成 |
%[text:table]
%[text] 
%[text] **典型的なパイプライン**
%[text]     mols -\> isValid -\> batchCalculate -\> filter.lipinski -\> writeSdf
%[text] 
%[text] **おめでとうございます！** Foundation チュートリアル（F01〜F06）を完走しました。
%[text] **次のステップ**: `examples/stories/` の応用ストーリーチュートリアル（Layer 2）に進んでみましょう。
%%
%[text] ## 演習
%[text] 各演習は `answers/f06_answers.m` を参照する前に、自分で解いてみましょう。
%[text] 
%[text] **E1.** 5 分子の SMILES リストファイルを書き出し、`readSmilesList` で読み直して
%[text]     ラウンドトリップを検証しましょう。SDF 往復の前後でそれぞれの正規 SMILES を
%[text]     表示して、一致することを確認してください。
runDir = makeRunDir();
% SMILES ファイル書き出し、readSmilesList で読み込み、SDF 書き出し、SDF 読み直し
% ラウンドトリップ前後で toSmiles の出力を比較する
%[text] 
%[text] **E2.** `everyday_chemicals.csv` を読み込み、記述子を計算して Lipinski フィルタを適用した後、
%[text]     通過した分子を SDF ファイル `"everyday_ro5_pass.sdf"` に保存しましょう。
%[text]     日用化学品 30 種類のうち、Ro5 を通過するのは何個ですか？
data = readtable("data/list/everyday_chemicals.csv", TextType="string");
% Analyze -> Descriptors -> lipinski -> writeSdf
%[text] 
%[text] **E3.** ミニパイプラインチャレンジ: `fda_drugs.csv` から始めて
%[text]     (a) Ro5 を通過する分子をフィルタリングし、
%[text]     (b) 通過した分子すべてについて Morgan フィンガープリントを計算し、
%[text]     (c) アスピリン（`"CC(=O)Oc1ccccc1C(=O)O"`）との類似度でランキングし、
%[text]     (d) アスピリンに最も類似した上位 10 個の FDA 薬を SDF ファイルに保存してください。
%[text]     F02・F03・F04・F06 で学んだスキルの集大成です。
%[text] 

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:96fa9617]
%   data: {"dataType":"text","outputData":{"text":"[09:27:01][INFO]  initPython: Python 3.10 already active (Status: Loaded) -- skipping.\n","truncated":false}}
%---
%[output:8dd385bb]
%   data: {"dataType":"text","outputData":{"text":"[09:27:01][INFO]  F06: セットアップ完了 -- 出力ディレクトリ: result\\runs\\20260525_092701\n","truncated":false}}
%---
%[output:390d8d09]
%   data: {"dataType":"text","outputData":{"text":"[09:27:01][INFO]  readSmilesList: loaded 5 molecules from 'result\\runs\\20260525_092701\\sample_molecules.smi' (0 skipped)\n","truncated":false}}
%---
%[output:2c91ec2c]
%   data: {"dataType":"text","outputData":{"text":"[09:27:01][INFO]  SMILES リストから 5 分子を読み込んだ\n","truncated":false}}
%---
%[output:72cef0d6]
%   data: {"dataType":"text","outputData":{"text":"[09:27:02][INFO]  writeSdf: wrote 5 molecules to 'result\\runs\\20260525_092701\\sample_output.sdf'\n","truncated":false}}
%---
%[output:98ae2446]
%   data: {"dataType":"text","outputData":{"text":"[09:27:02][INFO]  5 分子を書き出した: result\\runs\\20260525_092701\\sample_output.sdf\n","truncated":false}}
%---
%[output:998b8a2d]
%   data: {"dataType":"text","outputData":{"text":"[09:27:02][INFO]  readSdf: loaded 5 molecules from 'result\\runs\\20260525_092701\\sample_output.sdf' (0 skipped)\n","truncated":false}}
%---
%[output:91439fb7]
%   data: {"dataType":"text","outputData":{"text":"[09:27:02][INFO]  SDF から 5 分子を再読み込みした\n","truncated":false}}
%---
%[output:82318383]
%   data: {"dataType":"text","outputData":{"text":"[09:27:02][INFO]  --- ラウンドトリップ SMILES 確認 ---\n","truncated":false}}
%---
%[output:66806ecf]
%   data: {"dataType":"text","outputData":{"text":"[09:27:02][INFO]    Mol 1: CCO\n[09:27:02][INFO]    Mol 2: c1ccccc1\n[09:27:02][INFO]    Mol 3: CC(=O)Oc1ccccc1C(=O)O\n[09:27:02][INFO]    Mol 4: Cn1c(=O)c2c(ncn2C)n(C)c1=O\n[09:27:02][INFO]    Mol 5: CC(=O)Nc1ccc(O)cc1\n","truncated":false}}
%---
%[output:4a25a224]
%   data: {"dataType":"text","outputData":{"text":"[09:27:02][INFO]  FDA 薬物データセット: 200 化合物\n","truncated":false}}
%---
%[output:8597961b]
%   data: {"dataType":"text","outputData":{"text":"[09:27:13][INFO]  有効な FDA 分子: 200 \/ 200\n","truncated":false}}
%---
%[output:603355c3]
%   data: {"dataType":"text","outputData":{"text":"[09:27:16][INFO]  lipinski: 176 \/ 200 row(s) pass Ro5 (MaxViolations=0)\n","truncated":false}}
%---
%[output:28b668da]
%   data: {"dataType":"text","outputData":{"text":"[09:27:16][INFO]  Ro5 通過: 176 \/ 200 FDA 分子\n","truncated":false}}
%---
%[output:4ec9371a]
%   data: {"dataType":"text","outputData":{"text":"[09:27:21][INFO]  writeSdf: wrote 176 molecules to 'result\\runs\\20260525_092701\\fda_ro5_pass.sdf'\n","truncated":false}}
%---
%[output:3e55f303]
%   data: {"dataType":"text","outputData":{"text":"[09:27:21][INFO]  Ro5 通過分子を保存した: result\\runs\\20260525_092701\\fda_ro5_pass.sdf\n","truncated":false}}
%---
%[output:93bcbda8]
%   data: {"dataType":"text","outputData":{"text":"[09:27:21][INFO]  記述子テーブルを保存した: result\\runs\\20260525_092701\\fda_ro5_pass_descriptors.csv\n","truncated":false}}
%---
