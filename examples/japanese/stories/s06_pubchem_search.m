%[text] # S06: PubChem で化合物を検索する
%[text] EasyMolKit 応用ストーリー — レイヤー 2
%[text] 
%[text] 世界最大の無料化学データベース PubChem には 1 億件以上の化合物が登録されています。ブラウザで 1 件ずつ調べるのは非効率ですが、MATLAB から直接クエリを投げれば数十件の化合物データを数秒で取得して比較できます。
%[text] このスクリプトでは PubChem PUG REST API を通じて解熱鎮痛薬 5 種のデータを取得し、RDKit と組み合わせて物性比較テーブルを構築します。
%[text] ### ストーリー
%[text] あなたは小規模な製薬会社の研究者です。
%[text] あなたのチームは新しい抗炎症薬のプロジェクトを立ち上げたばかりで、競合となる既存の解熱鎮痛薬（アスピリン、イブプロフェン、アセトアミノフェン、ナプロキセン、セレコキシブ）のプロファイルを迅速に把握する必要があります。
%[text] 通常なら PubChem のウェブサイトを1件ずつ手動でブラウジングする作業になりますが、ここでは`emk.db.searchPubchem` 関数を使い、MATLAB から直接 API を叩いて一括で構造と物性データを取得・比較してみましょう。
%[text] **演習の流れ**
%[text] 1. 一般名で単一の化合物を検索する
%[text] 2. SMILES 文字列から化合物を検索する（構造はわかるが名前がわからない場合）
%[text] 3. 化合物リストをクエリして比較テーブルを作成する
%[text] 4. PubChem から取得した構造を可視化する
%[text] 5. テーブルを並べ替え・フィルタして最良のプロファイルを持つ候補を見つける \
%[text] ### 学習目標
%[text] - `emk.db.searchPubchem` で PUG REST API 経由で PubChem に問い合わせる
%[text] - 3 つのクエリモードを使い分ける: name・SMILES・CID
%[text] - API 結果から複数化合物の比較テーブルを構築する
%[text] - API 取得物性と Lipinski Ro5 フィルタを組み合わせる
%[text] - 取得した SMILES 文字列から構造を可視化する \
%[text] ### 前提条件
%[text] - F01（分子の描画）と F02（物性計算）の完了
%[text] - **インターネット接続が必要**（PubChem PUG REST API）
%[text] - RDKitインストール済み（`emk.setup.install()` を一度だけ実行しておく）
%[text] - 追加Toolbox不要（MATLAB だけで動きます） \
%[text] **所要時間**: 25〜40 分 | 実行方法: Ctrl+Enter でセクションを一つずつ実行
%[text] 
%[text] **参照文献**: 
%[text] - Kim S et al. (2023) PubChem 2023 update. *Nucleic Acids Res* 51:D1373-D1380.
%[text] - doi:10.1093/nar/gkac956 — PubChem PUG REST: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest 〔Open Access〕\
%[text] **注意**: このストーリーはライブ PubChem API にクエリを送信します。インターネット接続が必要です。PubChem がデータベースを更新すると、結果がわずかに異なる場合があります。
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
logInfo("S06: セットアップ完了");
%%
%[text] ## セクション 1: 一般名で単一化合物を検索する
%[text] ### PubChem PUG REST API
%[text] PubChem は米国国立医学図書館（NIH）が管理する世界最大の無料化学情報データベースです。膨大なデータが収録されており、創薬研究のインフラとなっています：
%[text] - Compound データベース：約 1 億 1500 万件のユニークな化合物構造（標準化済み）
%[text] - Substance データベース：約 3 億件の未精製の登録構造（由来情報を含む）
%[text] - BioAssay データベース：約 2 億 7000 万件の生物活性試験データ \
%[text] このデータベースにプログラムから直接アクセスするための仕組みが「PUG REST API」です。API のリクエスト（URL）は基本的に以下のルール（Webパターン）に従って構築されます：
%[text] [https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/\<input\_type\>/\<query\>/property/\<property\_list\>/JSON](https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/%3Cinput_type%3E/%3Cquery%3E/property/%3Cproperty_list%3E/JSON)
%[text] 【各パラメータの意味】
%[text] - \<input\_type\>  : 検索キーの種類（name（一般名）, smiles（構造数式）, cid（化合物ID）など）
%[text] - \<query\>       : 具体的な検索キーワード（"aspirin" や構造文字列など）
%[text] -  \<property\_list\> : 取得したい物性項目のリスト（カンマ区切り） \
%[text] 【APIリクエストの具体例】
%[text] - 一般名から物性を取得：   /compound/name/aspirin/property/MolecularWeight,IsomericSMILES/JSON
%[text] - SMILESから分子式を取得： /compound/smiles/CCO/property/MolecularFormula,MolecularWeight/JSON
%[text] - 化合物IDからIUPAC名を取得：/compound/cid/2244/property/IUPACName,MolecularFormula/JSON \
%[text]  `emk.db.searchPubchem` 関数は、この複雑なURL生成とHTTP通信、およびレスポンス（JSON）の解析を内部で自動的に行い、MATLAB で扱いやすい「テーブル型」にラップして返してくれます。まずは、1897年に発見された古典的な解熱鎮痛薬である「アスピリン（アセチルサリチル酸）」を一般名で検索し、APIの挙動を確認してみましょう。
logInfo("一般名でアスピリンを PubChem に照会中...");
tbl_aspirin = emk.db.searchPubchem("aspirin");

logInfo("アスピリンの PubChem 結果:");
logInfo("  CID              : %d",   tbl_aspirin.CID(1));
logInfo("  IUPAC 名         : %s",   tbl_aspirin.IUPACName(1));
logInfo("  分子式           : %s",   tbl_aspirin.MolecularFormula(1));
logInfo("  分子量           : %.2f g/mol", tbl_aspirin.MolecularWeight(1));
logInfo("  SMILES           : %s",   tbl_aspirin.IsomericSMILES(1));
%[text] **✏️ やってみよう 1**
%[text] Q: アスピリンの PubChem CID は？
%[text]    CID（Compound ID）はデータベース内で化合物を一意に識別する不変の番号です。
%[text]    ヒント: `tbl\_aspirin.CID(1)` の値を確認してみましょう。
%[text]    期待値: `2244`
%[text] Q: PubChem が返す IUPAC 名は「aspirin」よりずっと系統的。
%[text]    構造について何がわかるか？
%[text]    期待値: `2-acetyloxybenzoic acid`（安息香酸の2位（オルト位）にアセチルオキシ基が結合した構造であることが名前から読み取れます）。 
%[text] Q: "paracetamol" を検索してみよう -- これは「acetaminophen」の英国名。
%[text]    両方のクエリが同じ CID を返すか？
%[text]    世界中で異なる一般名（英名：パラセタモール、米名/日局名：アセトアミノフェン）で呼ばれていても、同じ化合物を指していれば PubChem は同一の CID を返します。  
%[text]    ヒント: コマンドウィンドウで `emk.db.searchPubchem("paracetamol")` と `emk.db.searchPubchem("acetaminophen")` を実行し、双方の CID を比較してください。
%[text]    期待値: 両方とも CID `1983` を返す。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: SMILES で検索する -- 構造がわかる場合
%[text] ### 構造ベースの検索
%[text] 実際の創薬ワークフローでは、「ChemDrawなどで描いた構造」や「計算化学でデザインした新規構造」はあるものの、その正確な化学名がわからないケースや、そもそも既知の化合物としてデータベースに登録されているか確かめたいケースが頻発します。
%[text] 化合物の構造表現である「SMILES 文字列」をクエリとして PubChem に問い合わせることで、以下の強力なメリットが得られます：
%[text] 1. その構造が「既知（登録済み）」かどうかのクイックな判定
%[text] 2. データベースに登録されている一意の識別子（CID）や正規の IUPAC 名の特定
%[text] 3. 公的機関によって検証・標準化された信頼性の高い基本物性データの取得 \
%[text] ※ クエリとして送信する SMILES は、RDKit 等で正しく解析できる有効な構造である必要があります。PubChem は内部で独自の「正規化（Canonical化）」を行った上でデータベースと照合するため、入力したSMILESの表記順が少し異なっていても正しくヒットします。
%[text] ここでは、市販の解熱鎮痛薬として広く普及している「イブプロフェン（R/Sラセミ体）」の構造（SMILES）を使って、逆引き検索を試してみましょう。
IBU_SMILES = "CC(C)Cc1ccc(cc1)C(C)C(=O)O";

logInfo("SMILES でイブプロフェンを PubChem に照会中...");
tbl_ibu = emk.db.searchPubchem(IBU_SMILES, Type="smiles");

logInfo("イブプロフェンの PubChem 結果（Type=smiles）:");
logInfo("  CID              : %d",   tbl_ibu.CID(1));
logInfo("  IUPAC 名         : %s",   tbl_ibu.IUPACName(1));
logInfo("  分子式           : %s",   tbl_ibu.MolecularFormula(1));
logInfo("  分子量           : %.2f g/mol", tbl_ibu.MolecularWeight(1));
logInfo("  SMILES（PubChem）: %s",   tbl_ibu.IsomericSMILES(1));

%[text] **✏️ やってみよう 2**
%[text] Q: PubChem は入力した SMILES と同じものを返すか、
%[text]    それとも別の（正規化された）形を返すか？
%[text]    SMILES は一意でない -- 同じ分子を表す多くの有効な SMILES 文字列がある。
%[text]    PubChem は常に自身の正規形を返す。
%[text] Q: イブプロフェンの PubChem CID は？ 期待値: 3672。
%[text] Q: 薬理学的に活性な（S）体イブプロフェンを検索してみましょう。
%[text]    異なる CID を持つか？
%[text]    (S)-イブプロフェン SMILES: "[C@@H](C(=O)O)(Cc1ccc(cc1)CC(C)C)C"
%[text]    ヒント: `emk.db.searchPubchem` に SMILES `"[C@@H](C(=O)O)(Cc1ccc(cc1)CC(C)C)C"` と `Type="smiles"` を渡します。
%[text]    注: (S) 体は立体化学が異なる構造として PubChem に登録されているため
%[text]    独自の CID を持つ。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: 解熱鎮痛薬パネルの比較テーブルを作成する
%[text] ### 薬物パネル比較
%[text] 特定の疾患をターゲットにした「同種・同効薬のセット（薬物パネル）」や、研究中の「合成化合物シリーズ（化学シリーズ）」を横並びで一括比較することは、医薬化学（メディシナルケミストリー）におけるもっとも基本的な作業です。構造のわずかな違いが物性や活性にどう影響するか（構造活性相関：SAR）を読み解くために、まずはデータを均一なテーブルにまとめる必要があります。
%[text] 本セクションでは、1世紀以上にわたる解熱鎮痛・抗炎症薬の進化を代表する「5つの有名医薬品」をパネルとして扱い、一挙にデータを取得します：
%[text] - アスピリン         : 古典的な非選択的 COX 阻害薬（1897年登場）
%[text] - イブプロフェン     : 代表的な非ステロイド性抗炎症薬（NSAID / 1961年登場）
%[text] - アセトアミノフェン : 主に中枢神経系に作用するとされる解熱鎮痛薬（1956年登場）
%[text] - ナプロキセン       : 長時間作用型の非選択的 NSAID（1976年登場）
%[text] - セレコキシブ       : 胃腸障害の副反応を抑えるため、COX-2 を特異的に標的とした選択的阻害薬（1998年登場） \
%[text] アスピリンのような単純な小分子から、標的への選択性を高めるために最適化されたセレコキシブのような現代的分子まで、分子の「進化の歴史」が物性データにどう現れるかを MATLAB のループ処理を用いて一括抽出・比較してみましょう。
PANEL_NAMES = ["aspirin"; "ibuprofen"; "acetaminophen"; "naproxen"; "celecoxib"];
N = numel(PANEL_NAMES);

logInfo("%d 種の化合物を PubChem に照会中...", N);

%[text] 結果配列の事前確保
cids     = zeros(N, 1, 'uint32');
iupac    = strings(N, 1);
formulas = strings(N, 1);
mws      = zeros(N, 1);
panelSmiles = strings(N, 1);

for i = 1:N
    name = PANEL_NAMES(i);
    try
        r = emk.db.searchPubchem(name);
        cids(i)     = r.CID(1);
        iupac(i)    = r.IUPACName(1);
        formulas(i) = r.MolecularFormula(1);
        mws(i)      = r.MolecularWeight(1);
        panelSmiles(i) = r.IsomericSMILES(1);
        logInfo("  [%d/%d] %-15s -> CID %d, MW=%.1f", i, N, name, cids(i), mws(i));
    catch ME
        logWarn("  [%d/%d] %s: クエリ失敗 -- %s", i, N, name, ME.message);
    end
end

%[text] 比較テーブルを組み立てる
panelTbl = table(PANEL_NAMES, cids, iupac, formulas, mws, panelSmiles, ...
    'VariableNames', {'Name', 'CID', 'IUPACName', 'MolecularFormula', ...
                      'MolecularWeight', 'SMILES'});

logInfo("--- 解熱鎮痛薬パネル（PubChem より）---");
disp(panelTbl(:, {'Name', 'CID', 'MolecularFormula', 'MolecularWeight'}));

%[text] **✏️ やってみよう 3**
%[text] Q: パネル内で最も分子量が高い化合物はどれか？
%[text]    ヒント: panelTbl(panelTbl.MolecularWeight == max(panelTbl.MolecularWeight), :)
%[text]    期待値: セレコキシブ（$&dollar&;MW \\approx 381~\\text{Da}&dollar&;$）。
%[text]    現代の創薬では大きくより選択的な分子への傾向を反映している。
%[text] Q: アスピリンとアセトアミノフェンはともに小分子（$&dollar&;MW \< 180~\\text{Da}&dollar&;$）。
%[text]    サイズが似ているのに、作用機序が全く異なる。
%[text]    これは MW と活性の関係について何を示しているか？
%[text]    （ヒント: サイズだけでは機序を決定しない; 官能基が決める。）
%[text] Q: 分子式を見ると、硫黄（S）を含む化合物はどれか？
%[text]    ヒント: panelTbl(contains(panelTbl.MolecularFormula, "S"), :)
%[text]    期待値: セレコキシブ（C17H14F3N3O2S -- スルホンアミド基が COX-2 選択性の鍵）。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: RDKit 記述子を計算して PubChem データと統合する
%[text] ### API データとローカル計算の組み合わせ
%[text] 外部のWeb API（PubChemなど）は非常に便利ですが、創薬研究で必要となるすべての高度な化合物記述子（例えば詳細な脂溶性指標や特定の三次元記述子など）を網羅しているとは限りません。そこで実務でよく使われるのが、「基本情報は API から取得し、より専門的な記述子は手元のローカル環境で計算して組み合わせる」というハイブリッドなワークフローです。
%[text] 【ハイブリッド・データ構築の流れ】
%[text] 1. PubChem API から化合物の正しい「構造情報（標準化されたSMILES）」と「基本情報」を取得する
%[text] 2. 取得した SMILES をローカルのケモインフォマティクスエンジン（RDKit）にインプットする
%[text] 3. `emk.descriptor.calculate` を使って、APIからは得られなかった記述子（LogP, TPSA, 各種原子・環カウントなど）をローカルで高速に一括計算する
%[text] 4. MATLAB の強力なテーブル結合機能を用いて、すべてのデータを1つの洗練されたデータフレーム（行列テーブル）に統合する \
%[text] このセクションを通じて、データサイエンスや機械学習インプットにもそのまま応用できる「完全な化合物物性プロファイルテーブル」の構築手法を学びます。
PROP_NAMES = ["LogP", "TPSA", "NumHDonors", "NumHAcceptors", ...
              "NumRotatableBonds", "RingCount"];

nProps   = numel(PROP_NAMES);
descMat  = nan(N, nProps);     % 行: 化合物、列: 記述子
validIdx = false(N, 1);

logInfo("PubChem SMILES から RDKit 記述子を計算中...");

for i = 1:N
    smi = panelSmiles(i);
    if strlength(smi) == 0
        logWarn("  %s をスキップ: PubChem から SMILES を取得できなかった", PANEL_NAMES(i));
        continue;
    end
    if ~emk.mol.isValid(smi)
        logWarn("  %s をスキップ: RDKit で SMILES を解析できない", PANEL_NAMES(i));
        continue;
    end
    mol = emk.mol.fromSmiles(smi);
    d   = emk.descriptor.calculate(mol, PROP_NAMES);
    for j = 1:nProps
        descMat(i, j) = d.(PROP_NAMES(j));
    end
    validIdx(i) = true;
end
%[text] 記述子テーブルを構築
descTbl = array2table(descMat, 'VariableNames', cellstr(PROP_NAMES));
%[text] PubChem テーブルと記述子テーブルを結合
fullTbl = [panelTbl(:, {'Name', 'CID', 'MolecularFormula', 'MolecularWeight'}), ...
           descTbl];

logInfo("--- 完全物性テーブル（PubChem + RDKit）---");
disp(fullTbl);
%[text] **✏️ やってみよう 4**
%[text] Q: 5 種すべてがリピンスキーのルール・オブ・ファイブ（$&dollar&;MW \< 500&dollar&;$、
%[text]    $&dollar&;\\text{LogP} \\le 5&dollar&;$、$&dollar&;\\text{HBD} \\le 5&dollar&;$、$&dollar&;\\text{HBA} \\le 10&dollar&;$）に合格することを
%[text]  `emk.filter.lipinski` で確認してみましょう。
%[text]    ヒント:
%[text]      roTbl = table(fullTbl.MolecularWeight, fullTbl.LogP, ...
%[text]                    fullTbl.NumHDonors,      fullTbl.NumHAcceptors, ...
%[text]                    'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
%[text]      result = `emk.filter.lipinski(roTbl)`;
%[text]      disp(result)
%[text]    期待値: 5 種すべてが Pass\_Ro5 = true。
%[text] Q: TPSA が最も高い化合物はどれか？
%[text]    TPSA（位相的極性表面積）は膜透過性と相関する:
%[text]    $&dollar&;\\text{TPSA} \< 90~\\text{A}^2&dollar&;$ で良好な経口吸収が予測される。
%[text]    期待値: セレコキシブ（スルホンアミド＋ピラゾールで最高 TPSA）。
%[text]    それでもセレコキシブの TPSA は $&dollar&;90~\\text{A}^2&dollar&;$ 未満（$&dollar&;77.98~\\text{A}^2&dollar&;$）であり、経口吸収に問題はない。
%[text] Q: LogP が最も低い（最も親水性の高い）化合物はどれか？
%[text]    親水性は水溶性とバイオアベイラビリティに影響する。
%[text]    期待値: アスピリン（RDKit $&dollar&;\\text{LogP} \\approx 1.31&dollar&;$）またはアセトアミノフェン（RDKit $&dollar&;\\text{LogP} \\approx 1.35&dollar&;$）
%[text]    — 両者はほぼ同等で、このパネル内では最も親水性が高い。
%[text]    ※ 実験測定値ではアセトアミノフェン（$&dollar&;\\approx 0.46&dollar&;$）が最低だが、
%[text]      RDKit Crippen モデルは実験値と一致しないことがある。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: パネル構造を可視化する
%[text] ### 構造と物性の可視化
%[text] 数値としての物性データ（分子量やLogPなど）を眺めるだけでなく、それらの数字が「分子のどの部分（官能基や骨格）」に由来しているのかを、2次元構造式を横に並べて視覚的に観察することは、医薬化学における非常に重要なアプローチです。
%[text] 今回の解熱鎮痛薬パネルを可視化するにあたり、以下の「化学構造上の特徴」と「実際の薬理特性・副作用」のリンクに着目してみましょう：
%[text] - カルボン酸基（-COOH）の有無 : アスピリン、イブプロフェン、ナプロキセンは共通してカルボン酸を持ちます。これが酸性を示し、NSAIDs特有の胃痛や胃腸障害（胃粘膜への物理的・生理的刺激）の主な原因となります。
%[text] - アミド基（-NHCO-）の採用     : アセトアミノフェンはカルボン酸を持たず、アミド構造を有しています。この構造の違いが、胃に優しい（胃腸障害が出にくい）という臨床特性に直結しています。
%[text] - 特徴的な機能性官能基        : セレコキシブに導入された「スルホンアミド基（-SO2NH2）」や「トリフルオロメチル基（-CF3）」が、分子量や脂溶性をどう変化させているかを視覚的に確認します。 \
%[text] PubChem から取得した高精度な SMILES データを `emk.viz.draw2d` に渡し、MATLAB のグラフィックス機能（subplot）を活かして、綺麗な構造式グリッドを描画してみましょう。
logInfo("パネル構造を描画中（PubChem SMILES より）...");

%[text] 全化合物を 1 つの figure にグリッドとして描画（F01 セクション 5 と同じパターン）
nCols = ceil(sqrt(N));
nRows = ceil(N / nCols);
figure("Name", "解熱鎮痛薬パネル -- 構造", "Position", [100 500 1100 560]);
for i = 1:N
    smi = panelSmiles(i);
    subplot(nRows, nCols, i);
    if strlength(smi) == 0 || ~emk.mol.isValid(smi)
        text(0.5, 0.5, "N/A", HorizontalAlignment="center");
        axis("off");
        continue;
    end
    mol = emk.mol.fromSmiles(smi);
    titleStr = sprintf("%s\n[CID %d] MW=%.0f", PANEL_NAMES(i), cids(i), mws(i));
    emk.viz.draw2d(mol, Title=titleStr, Width=320, Height=300);
end

%[text] **✏️ やってみよう 5**
%[text] 上で描画した 5 つの構造を観察してみましょう。それぞれについて識別してみましょう:
%[text]   (a) 酸性またはアミド官能基（機序・胃腸への影響に関わる）
%[text]   (b) 環系（1 環 vs. 2 環 vs. なし）
%[text]   (c) 環内のヘテロ原子（アスピリン/ナプロキセン: ベンゼンのみ;
%[text]       セレコキシブ: 窒素 2 個を持つピラゾール環）
%[text] Q: セレコキシブにはトリフルオロメチル基（CF3）とピラゾール環がある。
%[text]    これらはどちらも医薬化学の一般的な「バイオアイソスター」。
%[text]    なぜ CF3 基を導入するのか？
%[text]    （ヒント: CF3 は親油性を高め、代謝安定性を改善し、
%[text]     その部位での酸化的代謝を阻害しながらメチル基を模倣できる。）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: パネルをランク付け・フィルタする
%[text] ### 物性ベースの優先順位付け
%[text] ヒット化合物や既存薬の候補群から、目的の治療法に最適な「リード化合物」を絞り込む（あるいは優先順位をつける）際、研究者は医薬化学的なスクリーニング基準を用いて候補をフィルタリング・ソート（ランク付け）します。
%[text] 【一般的なフィルタリング基準の例】
%[text] 1. 経口医薬品らしさの評価（Lipinskiのルール・オブ・ファイブ：Ro5）  \
%[text]              \-分子量（MW）$&dollar&;\< 500&dollar&;$、脂溶性（LogP）$&dollar&;\\le 5&dollar&;$、水素結合供与体数（HBD）$&dollar&;\\le 5&dollar&;$、水素結合受容体数（HBA）$&dollar&;\\le 10&dollar&;$
%[text] 2\. 中枢神経系への移行性（CNS透過性 / 血液脳関門（BBB）の突破しやすさ）
%[text]             \-位相的極性表面積（TPSA）$&dollar&;\< 90~\\text{A}^2&dollar&;$（できれば $&dollar&;60~\\text{A}^2&dollar&;$ 以下がより理想的）
%[text]             \-分子量（MW）$&dollar&;\< 450&dollar&;$+    3. 副作用リスクの回避（例：胃腸刺激の低減）
%[text] 3\. カルボン酸などの強酸性官能基を回避、またはマスクできているか
%[text] このセクションでは、RDKitでローカル計算した「TPSA（極性表面積）」の値に基づき、脳への移行性が高い（中枢に到達しやすい）と予測される順にパネル化合物をソートし、グラフ化します。物性データに基づいたロジカルな意思決定プロセスの基礎を学びましょう。
logInfo("TPSA 昇順で化合物をランク付け（昇順 = CNS 透過性が最良）:");

[~, sortOrder] = sort(fullTbl.TPSA, "ascend");
ranked = fullTbl(sortOrder, :);
ranked.Properties.RowNames = {};

for i = 1:height(ranked)
    cnsFlag = "";
    if ranked.TPSA(i) < 60 && ranked.LogP(i) >= 1 && ranked.LogP(i) <= 3
        cnsFlag = " [CNS 透過性良好: TPSA<60 + LogP 1-3]";
    elseif ranked.TPSA(i) < 60
        cnsFlag = " [CNS 透過性: 低 TPSA]";
    elseif ranked.TPSA(i) < 90
        cnsFlag = " [CNS 透過性あり]";
    else
        cnsFlag = " [CNS 到達困難]";
    end
    logInfo("  %d. %-15s  TPSA=%.1f A^2  LogP=%.2f%s", ...
        i, ranked.Name(i), ranked.TPSA(i), ranked.LogP(i), cnsFlag);
end

%[text] 棒グラフ: TPSA 比較
figure("Name", "解熱鎮痛薬パネル -- TPSA 比較", "Position", [100 100 560 380]);
bar(fullTbl.TPSA(sortOrder));
xticks(1:N);
xticklabels(ranked.Name);
xtickangle(20);
ylabel("TPSA (A^2)");
title("解熱鎮痛薬パネル -- 位相的極性表面積");
yline(90,  "--r", "TPSA=90 （経口吸収閾値）", "LabelHorizontalAlignment", "left");
yline(60,  ":b", "TPSA=60 （CNS 推奨）", "LabelHorizontalAlignment", "left");
grid("on");

%[text] **✏️ やってみよう 6** — 物性プロファイルを読む
%[text] Q: TPSA（位相的極性表面積）が血液脳関門（BBB）透過性の指標として使われる理由を説明しましょう。
%[text]    極性の大きな分子はなぜ脂質二重層（細胞膜）を通りにくいのでしょうか？
%[text]    また、LogP が示す脂溶性は CNS 透過にどのように寄与しますか？
%[text]    （ヒント: CNS 透過にはリン脂質の二重膜を通る必要があります。
%[text]    極性表面積が大きいほど水溶性が上がり、膜透過が難しくなります。）
%[text] Q: セレコキシブはパネル内で最も高い TPSA を持ちます。
%[text]    血液脳関門を容易に通過できると思いますか？
%[text]    なぜ目的の用途（関節炎症）ではそれが許容されるのでしょうか？
%[text]    （ヒント: 炎症に対する COX-2 阻害に CNS 透過は不要です。
%[text]    セレコキシブは炎症組織の末梢 COX-2 を標的とします。）
%[text] Q: MW=500 の水平線を追加した分子量の棒グラフを作成しましょう。
%[text]    ヒント:
%[text]      figure; bar(fullTbl.MolecularWeight);
%[text]      xticks(1:N); xticklabels(fullTbl.Name); xtickangle(20);
%[text]      yline(500, "--r", "リピンスキー MW 上限"); grid("on");
%[text]    Lipinski MW 上限（500）に最も近い化合物はどれですか？
% ... （ここにコードを書いてみましょう）
%%
%[text] ## 演習
%[text] 
%[text] E1: CID で PubChem を検索する。
%[text]     抗マラリア薬キニーネの CID は 3034034。
%[text]     物性を取得して構造を描画してみましょう。
%[text]     ヒント:
%[text]       tbl\_q = `emk.db.searchPubchem("3034034", Type="cid")`;
%[text]       mol\_q = `emk.mol.fromSmiles(tbl_q.IsomericSMILES(1)`);
%[text]  `emk.viz.draw2d(mol_q, Title="Quinine")`;
%[text]  `emk.descriptor.calculate(mol_q, ["MolWt","LogP","TPSA","RingCount"])`
%[text] 
%[text] E2: InChIKey で検索する。
%[text]     アスピリンの InChIKey は "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"。
%[text]     一般名検索と同じ CID が返ることを確認してみましょう。
%[text]     ヒント: `emk.db.searchPubchem("BSYNRYMUTXBXSQ-UHFFFAOYSA-N", Type="inchikey")`
%[text]     期待値: CID 2244。
%[text] 
%[text] E3: パネルを拡張する。
%[text]     ジクロフェナク（別の NSAID）とトラマドール（オピオイド鎮痛薬）を
%[text]     パネルに追加して比較テーブルと棒グラフを再構築してみましょう。
%[text]     TPSA ランキングはどう変わるか？
%[text]     ヒント: PANEL\_NAMES に名前を追加してセクション 3〜6 を再実行。
%[text]     ジクロフェナク SMILES: "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl"
%[text]     トラマドール: 一般名で検索（"tramadol"）。
%[text] 
%[text] E4: エラーハンドリング。
%[text]     PubChem に存在しない名前を検索するとどうなるか？
%[text]     ヒント: `emk.db.searchPubchem("xyzzy_notacompound123")` を試す
%[text]     期待値: ID emk:db:searchPubchem:notFound のエラーが発生する。
%[text]     エラーを穏やかに処理するには:
%[text]       try
%[text]           r = `emk.db.searchPubchem("xyzzy_notacompound123")`;
%[text]       catch ME
%[text]           logWarn("見つからない: %s", ME.message);
%[text]       end
%[text] 
%[text] E5: 解熱鎮痛薬パネルの類似度ヒートマップを作成する。
%[text]     PubChem から全 5 化合物の SMILES を取得済みなので、
%[text]     モルガンフィンガープリントを計算してタニモト類似度行列を構築してみましょう。
%[text]     最も類似している 2 化合物はどれか？
%[text]     ヒント:
%[text]       fps = cell(1, N);
%[text]       for i = 1:N
%[text]           fps{i} = `emk.fingerprint.morgan(emk.mol.fromSmiles(panelSmiles(i)`));
%[text]       end
%[text]       S = `emk.similarity.matrix(fps)`;
%[text]       figure; imagesc(S); colormap("hot"); colorbar; clim(\[0 1\]);
%[text]       xticks(1:N); xticklabels(PANEL\_NAMES); xtickangle(20);
%[text]       yticks(1:N); yticklabels(PANEL\_NAMES);
%[text]       title("解熱鎮痛薬パネル -- タニモト類似度行列");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
