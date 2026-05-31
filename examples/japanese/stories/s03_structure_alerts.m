%[text] # S03: 構造アラート — 危険な官能基を SMARTS で検出する
%[text] EasyMolKit 応用ストーリー — レイヤー 2
%[text] 
%[text] 実験室に送られてくる膨大な新化合物の中に、タンパク質やDNAに直接反応してダメージを与えてしまう「化学的な時限爆弾」が潜んでいるとしたら、あなたはどうやって見つけますか？
%[text] こうした生体にとって有害な反応性を持つ部位（構造アラート）は、**SMARTS（スマート）パターン**と呼ばれる特殊な文字列をたった1行書くだけで、コンピュータ上で瞬時に検出し、スクリーニング（足切り）することができます。
%[text] さらに、製薬業界にはPAINS（ペインズ：Pan-Assay Interference Compounds）と呼ばれる、実験において「本当は効いていないのに、装置のセンサーを誤作動させて薬のように見せかける」という非常に厄介な（まさに研究者の"痛みの種"となる）偽陽性物質のリストが存在します。
%[text] このスクリプトでは、8種のカスタム反応性アラートと、最も厳格なPAINSフィルターをFDA承認薬200種に対して一括適用し、問題のある分子を自動で見つけ出してレポートを生成するプロセスを体験します。
%[text] ## 学習目標
%[text] - 構造の一部を検索するための条件式である「SMARTS（スマート）パターン」の基本を理解する
%[text] - `emk.mol.hasSubstruct` を使い、ループを使わずに複数分子の一括部分構造検索（ベクトル化スクリーニング）を行う
%[text] - 創薬における「PAINS」の概念と、アッセイ偽陽性を排除する重要性（トリアージ）を理解する
%[text] - スクリーニング結果を自動で整理し、構造化されたフラグテーブルやCSVレポートを作成する \
%[text] ## 前提条件
%[text] - F05（部分構造検索）の完了
%[text] - RDKitインストール済み（`emk.setup.install()` を一度だけ実行しておく）
%[text] - 追加Toolbox不要（MATLAB だけで動きます） \
%[text] **所要時間**: 15〜20 分 | 実行方法: Ctrl+Enter でセクションを一つずつ実行
%[text] **データ**
%[text] - `data/list/pains.csv` — PAINS SMARTS 480 種（RDKit wehi\_pains、BSD-3）
%[text] - `data/list/fda_drugs.csv` — FDA 承認薬 200 種（ChEMBL、CC-BY-SA 3.0） \
%[text] **参考文献**
%[text] - Baell JB & Holloway GA (2010) *J Med Chem* 53:2719–2740. doi:10.1021/jm901137j 〔要機関アクセス〕\
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
logInfo("S03: セットアップ完了");
%%
%[text] ## セクション 1: 構造アラートとは何か？
%[text] ### コンセプト: 生体分子と「直接反応」する化学基
%[text] 一部の化学基は、タンパク質や DNA などの生体分子と直接反応する性質があります。創薬の初期段階では、こうした部位を「構造アラート」としてフラグを立て、明確な理由がない限り優先度を下げます。
%[text] よく知られたアラートとそのリスクの一例です:
%[text] -   エポキシド          `[C]1CO1`           — DNA・タンパクをアルキル化
%[text] -   マイケル受容体      `[C]=[C]-C=O`       — システイン(-SH)と 1,4-付加
%[text] -   アルデヒド          `[CH]=O`            — アミンと Schiff 塩基を形成
%[text] -   アシルハライド      `[C](=O)[F,Cl]`     — 高反応性の求電子剤
%[text] -   ニトロ基            `[N+](=O)[O-]`      — 潜在的変異原性（Ames 試験）
%[text] -   ヒドラジン          `[NH]-[NH2]`        — 肝毒性の懸念
%[text] -   ジアゾ              `[#6]=[N+]=[N-]`    — 高反応性・不安定
%[text] -   チオール            `[SH]`              — ファーマコフォアになる場合も \
%[text] ただし、アラートがあるだけで化合物が「使えない」わけではありません。たとえばアスピリンのエステルは体内で意図的に加水分解される設計（プロドラッグ）です。「フラグあり」はあくまでも専門家が判断するための合図です。
%[text] アラートパネルを定義します（名前 → SMARTS の対応）。
ALERT_NAMES  = ["Epoxide",      "Michael_acceptor", "Aldehyde", ...
                "Acyl_halide",  "Nitro",             "Hydrazine", ...
                "Diazo",        "Thiol"];
ALERT_SMARTS = ["[C]1CO1",      "[C]=[C]-C=O",      "[CH]=O", ...
                "[C](=O)[F,Cl]","[N+](=O)[O-]",     "[NH]-[NH2]", ...
                "[#6]=[N+]=[N-]","[SH]"];

logInfo("アラートパネル読み込み完了: %d パターン", numel(ALERT_NAMES));
for k = 1:numel(ALERT_NAMES)
    logInfo("  %-18s  %s", ALERT_NAMES(k), ALERT_SMARTS(k));
end
%%
%[text] ## セクション 2: 既知の問題分子のテストセットをスクリーニングする
%[text] あらかじめ反応性のある官能基を含んでいることが分かっている「代表的な問題分子」を集めたテストセットを用意しました。EasyMolKitを使って、正しくアラートを検出できるか検証してみましょう。
TEST_NAMES  = ["エピクロロヒドリン",  "アクロレイン",    "ホルムアルデヒド", ...
               "マロンアルデヒド",    "ニトロベンゼン", "フェニルヒドラジン", ...
               "アスピリン",          "カプトプリル"];
TEST_SMILES = ["ClCC1CO1",          "C=CC=O",       "C=O", ...
               "O=CCC=O",          "c1ccc([N+](=O)[O-])cc1", ...
               "c1ccccc1NN",       "CC(=O)Oc1ccccc1C(=O)O", ...
               "CC(CS)C(=O)N1CCCC1C(=O)O"];

%[text] ### コンセプト: セル配列を渡してベクトル化スクリーニング
%[text] `emk.mol.hasSubstruct` の第 1 引数に Mol オブジェクトのセル配列を渡すと、分子ごとの結果を `logical(1, N)` 行ベクトルとして一度に受け取れます。明示的なループを書かなくても、全分子を一括スクリーニングできます。
nMols   = numel(TEST_SMILES);
nAlerts = numel(ALERT_NAMES);
flagMat = false(nMols, nAlerts);   % 行 = 分子、列 = アラート

mols_test = cell(1, nMols);
for i = 1:nMols
    mols_test{i} = emk.mol.fromSmiles(TEST_SMILES(i));
end

for j = 1:nAlerts
    flagMat(:, j) = emk.mol.hasSubstruct(mols_test, ALERT_SMARTS(j))';
end

%[text] 結果を読みやすいテーブルに整形します。
flagTbl = array2table(flagMat, "VariableNames", cellstr(ALERT_NAMES));
flagTbl.Molecule     = TEST_NAMES';
flagTbl.TotalAlerts  = sum(flagMat, 2);
flagTbl = movevars(flagTbl, ["Molecule", "TotalAlerts"], "Before", "Epoxide");

logInfo("--- 構造アラートスクリーン（テストセット）---");
disp(flagTbl(:, ["Molecule", "TotalAlerts", ALERT_NAMES]));

%[text] **✏️ やってみよう 1 — テーブルを読んでみましょう**
%[text] テーブルを確認して、次の 2 つの問いに答えてみてください。
%[text] - テスト分子の中で最も多くのアラートをトリガーするのはどれですか？
%[text] - アスピリンはアラートをトリガーしましたか？それは何を意味しますか？ \
%[text] 
%[text] **期待値**:
%[text] - アクロレイン（`C=CC=O`）は Aldehyde と Michael\_acceptor の **2 件**をトリガーします。
%[text] - マロンアルデヒドは Aldehyde のみ（1 件 — C=C 結合がないためマイケル受容体は非該当）。
%[text] - ホルムアルデヒド（`C=O`）は **0 件** — `[CH]=O` は炭素に H を 1 個要求しますが、 \
%[text]   ホルムアルデヒドの炭素は H を 2 個持つため（`[CH2]=O`）マッチしません。
%[text]   これは SMARTS の特異性を示す良い例です。「アルデヒドらしい」のにフラグされない点を覚えておきましょう。
%[text] - アスピリンはこのパネルでアラートなし — エステルはフラグされません。
%[text] - カプトプリル（ACE 阻害薬）は Thiol をトリガー（意図的なファーマコフォア）。 \
%[text] 
%[text] 「フラグあり」は専門家判断のトリガーであり、自動除外の根拠ではありません。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: フラグ付き分子を描く
%[text] テストセットで少なくとも 1 つのアラートをトリガーした分子を可視化してみましょう。
flaggedIdx = find(sum(flagMat, 2) > 0);
nFlagged_test = numel(flaggedIdx);
logInfo("テストセットのフラグ付き分子を %d 件描画中...", nFlagged_test);

%[text] フラグ付き分子をグリッドに描画します（F01 セクション 5 と同じパターン）。
nCols_f = ceil(sqrt(nFlagged_test));
nRows_f = ceil(nFlagged_test / nCols_f);
figure("Name", "フラグ付き分子 -- 構造アラート", ...
    "Position", [100 100 nCols_f*280 nRows_f*260]);
for k = 1:nFlagged_test
    i = flaggedIdx(k);
    subplot(nRows_f, nCols_f, k);
    alertList = ALERT_NAMES(flagMat(i, :));
    titleStr = sprintf("%s [%s]", TEST_NAMES(i), strjoin(alertList, ", "));
    emk.viz.draw2d(mols_test{i}, Title=titleStr);
end

%[text] **✏️ やってみよう 2 — アクロレインの 2 つの反応中心を確認しましょう**
%[text] アクロレイン（CH2=CH-CHO）を描いて、どこが反応しやすいか考えてみましょう。
%[text] 
%[text]     mol\_acrolein = emk.mol.fromSmiles("C=CC=O");
%[text]     emk.viz.draw2d(mol\_acrolein, Title="アクロレイン — マイケル + アルデヒド", ...
%[text]         Width=300, Height=250);
%[text] 
%[text] システインのチオール（-SH）は、アルデヒドの炭素と β 炭素のどちらを攻撃しやすいですか？
%[text] ヒント: チオールは「軟らかい」求核剤です。HSAB 理論を思い出してみましょう。
%[text] 
%[text] **期待値**: β 炭素（C=C の 2 番目の炭素、マイケル受容体部位）を攻撃します。
%[text] 軟求核剤は LUMO が共役系の末端（β 位）に広がるため、1,4-付加が速度論的に有利です。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: PAINS フィルタセットを読み込んで探索する
%[text] ### コンセプト: PAINS とは何か？
%[text] Baell & Holloway（2010）は 680 万件のアッセイデータを分析し、複数の生化学アッセイで偽陽性を繰り返し引き起こす化合物を特徴づける480 種の SMARTS パターンを特定しました。これが「PAINS」です。
%[text] 主な原因は蛍光干渉・酸化還元サイクリング・コロイド状凝集・光学干渉などです。480 パターンは厳格度に応じて 3 つのサブフィルタに分類されています。
%[text] - **PAINS\_A**: 16 パターン（最も厳格）
%[text] - **PAINS\_B**: 55 パターン
%[text] - **PAINS\_C**: 409 パターン（最も広いカバレッジ） \
%[text] 
%[text] **実装メモ**: PAINS CSV の SMARTS 文字列は引用符フィールド内にカンマを含みます。
%[text] `textscan` の `%q`（引用文字列形式）を使うと、埋め込みデリミタを正しく処理できます。
fid = fopen(fullfile(projectRoot, "data", "list", "pains.csv"), "r");
textscan(fid, "%s", 1, "Delimiter", newline);  % ヘッダー行を消費
C = textscan(fid, "%q%q%q%q", "Delimiter", ",");   % 4 列: Name,SMARTS,FilterSet,Source
fclose(fid);
painsCsv = table(string(C{1}), string(C{2}), string(C{3}), ...
    VariableNames=["Name", "SMARTS", "FilterSet"]);

logInfo("PAINS データベース: %d パターン読み込み完了", height(painsCsv));
filterSets = unique(painsCsv.FilterSet);
for k = 1:numel(filterSets)
    n = sum(painsCsv.FilterSet == filterSets(k));
    logInfo("  %-10s: %d パターン", filterSets(k), n);
end

%[text] **✏️ やってみよう 3 — PAINS データベースを認めてみましょう**
%[text] MATLAB ワークスペースで `painsCsv` を開いて、最初の 3 パターンを確認してみましょう。
%[text] ヒント: `painsCsv(1:3, ["Name","SMARTS","FilterSet"])` でテーブルを表示できます。
%[text] 
%[text] PAINS の SMARTS は複雑な再帰クエリです。人が目で解説しようとするより、
%[text] RDKit の SMARTS パーサに任せる方が正確です。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: FDA 薬物への PAINS バッチスクリーニング
%[text] PAINS パターンを 200 種の FDA 薬物に一括適用します。
%[text] 各薬物について以下の 3 つを記録します。
%[text] - いずれかの PAINS パターンに一致するか（`isPainsFlagged`）
%[text] - いくつの異なるパターンに一致するか（`painsCount`）
%[text] - 最初にトリガーされたパターン名（`firstAlert`） \
%[text] **パフォーマンス注記**: 480 パターン × 200 薬物 = 96,000 回の IPC 呼び出しになり、
%[text] 完了まで 10～30 分かかります。このデモでは最も厳格な PAINS\_A（16 パターン）だけを
%[text] 使用して 30 秒以内に完了させています。全 480 パターンのスクリーニング方法は
%[text] やってみよう 4 を参照してください。
rawDrugs = readtable(fullfile(projectRoot, "data", "list", "fda_drugs.csv"), "TextType", "string");
nDrugs   = height(rawDrugs);

%[text] インタラクティブデモ用に PAINS\_A！16 種・最も厳格）を使用します。全 480 パターンを試す場合は、この行を `painsScreen = painsCsv;` に変更してください。
painsScreen = painsCsv(painsCsv.FilterSet == "PAINS_A", :);
nPains      = height(painsScreen);

logInfo("%d 種の薬物を %d 件の PAINS_A パターンでスクリーニング中...", nDrugs, nPains);
logInfo("  （PAINS_A = 最も厳格な %d パターン; 完全データベースは %d 種）", nPains, height(painsCsv));

%[text] まず全薬物の SMILES を解析します（無効な SMILES は警告を出してスキップします）。
mols_drug = cell(1, nDrugs);
validDrug  = false(1, nDrugs);
for i = 1:nDrugs
    if emk.mol.isValid(rawDrugs.SMILES(i))
        mols_drug{i} = emk.mol.fromSmiles(rawDrugs.SMILES(i));
        validDrug(i) = true;
    else
        logWarn("  スキップ: %s（無効な SMILES）", rawDrugs.Name(i));
    end
end
logInfo("  薬物構造を解析した: %d / %d 件。", sum(validDrug), nDrugs);

%[text] 各 PAINS パターンで全薬物をスクリーニングします。`hitMat(i, j) = true` は薬物 i がパターン j にヒットしたことを意味します。
hitMat = false(nDrugs, nPains);
mols_valid_list = mols_drug(validDrug);
validIdx        = find(validDrug);

for j = 1:nPains
    logProgress(j, nPains, "PAINS スクリーン");
    smarts_j = painsScreen.SMARTS(j);
    hits_j   = emk.mol.hasSubstruct(mols_valid_list, smarts_j);   % 1×nValid logical
    hitMat(validIdx, j) = hits_j;
end

%[text] 薬物ごとにヒット数と最初のアラート名をまとめます。
painsCount = sum(hitMat, 2);          % 各薬物がヒットするパターン数
isPainsFlagged = painsCount > 0;

%[text] 最初にヒットしたパターン名を取得します（ヒットなしなら空文字のまま）。
firstAlert = repmat("", nDrugs, 1);
for i = 1:nDrugs
    idx = find(hitMat(i, :), 1);
    if ~isempty(idx)
        firstAlert(i) = painsScreen.Name(idx);
    end
end

%[text] 結果をテーブルにまとめます。
resultTbl = table(rawDrugs.Name, rawDrugs.SMILES, isPainsFlagged, painsCount, firstAlert, ...
    VariableNames=["Name", "SMILES", "PAINS_Flagged", "PAINS_Count", "First_Alert"]);

nFlagged = sum(isPainsFlagged);
logInfo("PAINS スクリーン完了:");
logInfo("  フラグ付き: %d / %d  (%.0f%%)", nFlagged, nDrugs, 100*nFlagged/nDrugs);
logInfo("  クリーン   : %d / %d  (%.0f%%)", nDrugs - nFlagged, nDrugs, ...
    100*(nDrugs-nFlagged)/nDrugs);

%[text] **✏️ やってみよう 4 — FDA 薬物の PAINS 割合を調べましょう**
%[text] - FDA 薬物の何割が少なくとも 1 つの PAINS\_A パターンをトリガーしましたか？
%[text] - 結果は驚くべきことですか？その理由を考えてみましょう。 \
%[text] 
%[text] **期待値**: PAINS\_A（16 パターン）では通常〜3～5% の承認薬がフラグされます。全 480 パターン時は 15～25% に上昇します（PAINS\_B/C はより寛容なパターン）。承認薬にフラグが多いのは自然なことです — 反応性モチーフ自体が薬効団の場合もあるからです。
%[text] 全 480 パターンのスクリーニング（～10～30 分）を試す場合は、セクション 5 の`painsScreen = painsCsv;` への変更で実行できます。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: 最も頻繁にトリガーされる PAINS パターンを特定する
%[text] パターンごとの発火回数を集計して、上位のアラートを特定しましょう。
patternHits = sum(hitMat, 1);   % 1×nPains

%[text] 発火回数の高い順（降順）でソートします。
[sortedHits, sortOrder] = sort(patternHits, "descend");
topK = min(10, sum(sortedHits > 0));   % 上位 10 件または発火数

if topK == 0
    logInfo("このデータセットでトリガーされた PAINS パターンはありませんでした。");
else
    topPatternNames = painsScreen.Name(sortOrder(1:topK));   % already column vector
    topHitCounts    = sortedHits(1:topK)';
    topTbl = table(topPatternNames, topHitCounts, ...
        VariableNames=["PatternName", "FlaggedDrugs"]);
    logInfo("上位 %d 件の PAINS パターン:", topK);
    disp(topTbl);

    % 棒グラフ: ヒット数にばらつきがある場合のみ描画
    % PAINS_A デモ（16 パターン・少数ヒット）では全パターンのカウントが均一になりやすく
    % 棒グラフに有意差が出ないためスキップする。全 480 パターン時に有意差が現れる。
    if max(topHitCounts) - min(topHitCounts) > 0
        figure("Name", "上位 PAINS パターン", "Color", "white", "Position", [100 100 580 380]);
        topNames = painsScreen.Name(sortOrder(1:topK));
        barh(topK:-1:1, sortedHits(1:topK), "FaceColor", [0.85 0.33 0.10]);
        yticks(1:topK);
        yticklabels(flip(cellstr(topNames)));
        xlabel("フラグ付き FDA 薬物数");
        title("FDA 薬物 200 種の上位 PAINS パターン");
        grid("on");
    else
        logInfo("  全パターンのヒット数が均一 (%d 件) のため棒グラフをスキップ。", topHitCounts(1));
        logInfo("  セクション 5 を painsScreen = painsCsv; に変更して全 480 パターンで実行すると");
        logInfo("  ヒット数の分布が現れ、棒グラフが有意になります。");
    end
end
%[text] **✏️ やってみよう 5 — PAINS ヒット数が最多の薬物を探しましょう**
%[text] 最も多くの PAINS パターンにヒットする薬物を調べてみましょう。
%[text] ヒント:
%[text]     topDrugs = sortrows(resultTbl, "PAINS\_Count", "descend");
%[text]     disp(topDrugs(1:5, \["Name", "PAINS\_Count", "First\_Alert"\]));
%[text] 
%[text] 上位に来た薬物には、その「アラート」が実際にその薬が効く理由の核心になっている
%[text] 可能性があります。薬物名 + 作用機序で検索して確認してみましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 7: フラグ付き化合物レポートを保存する
%[text] スクリーニング結果をまとめた CSV レポートを `result/runs/` に保存します。
runDir = makeRunDir();
outFile = fullfile(runDir, "s03_pains_report.csv");
writetable(resultTbl, outFile);

logInfo("レポートを保存した: %s", outFile);
logInfo("列: Name, SMILES, PAINS_Flagged, PAINS_Count, First_Alert");
logInfo("  resultTbl(resultTbl.PAINS_Flagged, :) でフラグ付き薬物を抽出する。");
logInfo("  resultTbl(~resultTbl.PAINS_Flagged, :) でクリーンな薬物を抽出する。");

logInfo("--- S03 完了 ---");
logInfo("サマリー:");
logInfo("  カスタムアラートパネル: %d パターンを %d 件のテスト分子でスクリーニング", ...
    nAlerts, nMols);
logInfo("  PAINS スクリーン: %d パターンを %d 種の FDA 薬物 -- %d 件フラグ（%.0f%%）", ...
    nPains, nDrugs, nFlagged, 100*nFlagged/nDrugs);
logInfo("  要点: PAINS フラグ = 調査のトリガー、自動除外ではない");
%[text] **✏️ やってみよう 6 — FDA 薬物にカスタムアラートも適用してみましょう**
%[text] セクション 1 で定義したカスタムアラートパネルを FDA 薬物に適用しましょう。
%[text] ヒント:
%[text]     drugAlertMat = false(nDrugs, nAlerts);
%[text]     for j = 1:nAlerts
%[text]         hits\_j = emk.mol.hasSubstruct(mols\_drug(validDrug), ALERT\_SMARTS(j));
%[text]         drugAlertMat(validIdx, j) = hits\_j;
%[text]     end
%[text]     alertFlagged = sum(drugAlertMat, 2) \> 0;
%[text]     logInfo("カスタムアラート: %d/%d 件フラグ", sum(alertFlagged), nDrugs);
%[text] 
%[text] カスタムアラートのフラグ率は PAINS より高いですか、低いですか？
%[text] どちらの基準がより厳格だと思いますか？

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
