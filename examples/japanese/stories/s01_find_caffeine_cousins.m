%[text] # S01: カフェインの仲間を探せ
%[text] EasyMolKit 応用ストーリー ─ レイヤー 2
%[text] 
%[text] コーヒー・お茶・チョコレートは、なぜあの独特な「覚醒感」をもたらすのでしょうか？
%[text] 答えはカフェインだけではありません。植物はカフェインと化学的によく似た覚醒物質のファミリーをまとめて産生しています。このスクリプトでは、分子フィンガープリントと類似度スコアを使って「カフェインの仲間」を 30 種の日用化学品データベースから探し出す体験ができます。
%[text] ## 学習目標
%[text] - `emk.fingerprint.morgan` で分子構造をビットベクトルとしてエンコードする
%[text] - `emk.similarity.rankBy` で化学品データベースを一括検索する
%[text] - タニモトスコアを実世界の文脈で解釈する
%[text] - 構造的類似性と生物活性の関連を理解する \
%[text] ## 前提条件
%[text] - F01（分子の描画）と F04（類似度）の完了
%[text] - RDKitインストール済み（`emk.setup.install()` を一度だけ実行しておく）
%[text] - 追加Toolbox不要（MATLAB だけで動きます） \
%[text] **所要時間**: 15〜20 分 | 実行方法: Ctrl+Enter でセクションを一つずつ実行
%[text] **データ**
%[text] - `data/list/everyday_chemicals.csv` — 30 種の一般分子（PubChem CC0） \
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
logInfo("S01: セットアップ完了");
%%
%[text] ## セクション 1: カフェイン
%[text] まずは、**カフェイン**の化学的性質を理解しましょう。
%[text] カフェイン（1,3,7-トリメチルキサンチン）はメチルキサンチン類に属する植物アルカロイドです。SMILESは、3つのメチル基（$&dollar;CH\_3&dollar;$）と2つのカルボニル酸素（$&dollar;=O&dollar;$）を持つ縮合二環式環（プリン骨格）をエンコードしています。ここでは、化学構造の読み込みと、分子の物理化学的プロパティ（分子量、脂溶性、極性表面積など）の計算方法を学びます。
%[text] ### カフェイン含有量の目安
%[text] - エスプレッソ 80 mg/ショット
%[text] - 紅茶 50 mg/カップ
%[text] - ダークチョコレート 20 mg/30 g
%[text] - デカフェコーヒー \< 5 mg/カップ \
CAFFEINE_SMILES = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C";
CAFFEINE_NAME   = "Caffeine";

mol_caffeine = emk.mol.fromSmiles(CAFFEINE_SMILES);
logInfo("カフェインを解析した。重原子数: %d", double(mol_caffeine.GetNumHeavyAtoms()));
desc = emk.descriptor.calculate(mol_caffeine, ["MolWt", "LogP", "TPSA", "RingCount"]);
logInfo("カフェインの性質:");
logInfo("  分子量             : %.2f g/mol  (C8H10N4O2)", desc.MolWt);
logInfo("  LogP               : %.2f        (わずかに親水性)", desc.LogP);
logInfo("  TPSA               : %.1f A^2   (CNS 透過は < 90 が目安)", desc.TPSA);
logInfo("  環数               : %d          (プリン二環式骨格)", desc.RingCount);
%[text] カフェインの構造を描く
figure("Name", "カフェイン", "Position", [100 100 440 380]);
emk.viz.draw2d(mol_caffeine, Title="カフェイン（1,3,7-トリメチルキサンチン）", ...
    Width=350, Height=350);
%[text] **✏️ やってみよう 1 — BBB 透過性を確認しましょう**
%[text] カフェインはCNS刺激薬です。血液脳関門（BBB）は $&dollar&;TPSA \< 90\\text{ \\AA}^2&dollar&;$ かつ $&dollar&;LogP&dollar&;$ があまり高くない（通常 $&dollar&;\< 5&dollar&;$）と透過しやすくなります。上の出力を確認して、カフェインはこの2つの条件を満たしていますか？
%[text] **期待される結果とヒント**
%[text] - **期待値**: $&dollar&;TPSA \\sim 62\\text{ \\AA}^2&dollar&;$ （\< 90 で良好）、$&dollar&;LogP \\sim -1.0&dollar&;$（負の値で非常に水溶性が高い）。
%[text] - **注釈**: 負の $&dollar&;LogP&dollar&;$ は、カフェインが水（コーヒーやお茶）に容易に溶けることを示します。低い $&dollar&;LogP&dollar&;$ にもかかわらず、カフェインはその小さいサイズと平面的な芳香族構造のおかげで、BBBを効率よく透過できます。 \
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: 日用化学品データベースを読み込む
%[text] CSV には食品・医薬品・日用品に含まれる 30 種の分子が収録されています。
%[text] 各行: CommonName, CID, SMILES, MolecularFormula, MolecularWeight, Category, Source
%[text] 探索対象となる「30種類の身近な化学物質」が登録された小規模なデータベース（CSVファイル）を読み込みます。このテーブルには、私たちが普段目にする食品（糖類やアミノ酸）、医薬品、日用品の成分が含まれており、それぞれに構造情報（SMILES）が紐づいています。どのようなカテゴリの分子が含まれているか、分布を確認してみましょう。
%[text] - `CommonName`: 一般名
%[text] - `CID`: PubChem 化合物 ID
%[text] - `SMILES`: 構造を記述した文字列
%[text] - `Category`: 「医薬品」「食品添加物」などの分類 \
dataFile = fullfile(projectRoot, "data", "list", "everyday_chemicals.csv");
tbl = readtable(dataFile, "TextType", "string");

logInfo("データベース読み込み完了: %d 分子 / %d カテゴリ", ...
    height(tbl), numel(unique(tbl.Category)));
%[text] カテゴリのプレビュー
cats = unique(tbl.Category);
logInfo("カテゴリ一覧:");
for k = 1:numel(cats)
    n = sum(tbl.Category == cats(k));
    logInfo("  %-15s: %d 分子", cats(k), n);
end
%[text] **✏️ やってみよう 2 — データベースを一覧してみましょう**
%[text] MATLAB ワークスペースでテーブル全体を閲覧してみましょう。
%[text] 最も多い分子を持つカテゴリはどれですか？
%[text] 最も分子量が高い分子はどれですか？
%[text] ヒント: max(tbl.MolecularWeight)  と  tbl(tbl.MolecularWeight == max(...), :)
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: データベース全体のフィンガープリントを作成する
%[text] コンピュータはSMILES文字列のままでは「形が似ているかどうか」をうまく計算できません。そこで、構造を計算可能な数学的表現に変換する必要があります。
%[text] **Morgan フィンガープリント**（ECFP4）は、分子を構成するそれぞれの原子と、その「周囲のつながり（半径2原子分）」をチェックして、分子の特徴を2048個の「0か1のフラグ（ビット）」に変換する技術です。いわば、**化学構造のデジタル指紋**のようなものです。似たような官能基や環構造（スキャフォールド）を持つ分子同士は、この指紋のパターン（1が立つ位置）がそっくりになります。
%[text] このセクションでは、先ほど読み込んだデータベースの全分子に対して、ループ処理を用いてフィンガープリントを一括生成します。
logInfo("%d 分子の Morgan フィンガープリントを計算中...", height(tbl));
fps   = cell(1, height(tbl));
valid = true(1, height(tbl));     % 正常に解析できた分子を追跡

for i = 1:height(tbl)
    smi = tbl.SMILES(i);
    if ~emk.mol.isValid(smi)
        logWarn("  スキップ: %s（無効な SMILES）", tbl.CommonName(i));
        valid(i) = false;
        continue;
    end
    mol    = emk.mol.fromSmiles(smi);
    fps{i} = emk.fingerprint.morgan(mol);
end

nValid = sum(valid);
logInfo("フィンガープリント計算完了: %d / %d 分子", nValid, height(tbl));
%[text] 有効なエントリでフィルタ済みリストを作成
fps_valid    = fps(valid);
names_valid  = tbl.CommonName(valid);
smiles_valid = tbl.SMILES(valid);
%[text] **✏️ やってみよう 3 — ON ビット数を数えてみましょう**
%[text] カフェインのフィンガープリントで ON になっているビット数は何個ですか？
%[text] ヒント:
%[text]   fp\_caf = `emk.fingerprint.morgan(mol_caffeine)`;
%[text]   bits   = `emk.fingerprint.toArray(fp_caf)`;
%[text]   sum(bits)
%[text] 期待値: 2048 ビット中 約 30〜50 ビットが ON。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: カフェインへの類似度で全分子をランク付け
%[text] emk.similarity.rankBy は、裏側で動いている Python/RDKit の高速な一括計算機能（BulkTanimotoSimilarity）を利用しています。
%[text] MATLAB と Python の間でデータを1件ずつ往復させて計算すると、分子の数（$&dollar&;N&dollar&;$ 件）だけ通信のオーバーヘッド（タイムロス）が発生してしまいます。しかし、この関数は「全$&dollar&;N&dollar&;$ 件のデータを1回でまとめて Python に送り、1回で結果を受け取る」という効率的な処理を行うため、データベース全体の検索を驚くほど一瞬で完了できます。
%[text] カフェインの指紋（ビットベクトル）と、データベース内すべての分子の指紋を比較し、タニモト係数（Tanimoto Coefficient）を用いて類似度を算出します。タニモト係数は 0.0（全く似ていない）から 1.0（完全一致）の範囲の値をとる指標です。
fp_caffeine = emk.fingerprint.morgan(mol_caffeine);

result = emk.similarity.rankBy(fp_caffeine, fps_valid);   % 全分子

logInfo("--- カフェインに最も類似した上位 10 分子 ---");
logInfo("%-5s  %-20s  %s", "順位", "分子名", "タニモト");
for k = 1:min(10, numel(result.Scores))
    idx   = result.Indices(k);
    name  = names_valid(idx);
    if strlength(name) == 0; name = "(名称なし)"; end
    score = result.Scores(k);
    % カフェイン自身（スコア = 1.0）はスキップ
    if score >= 1.0
        logInfo("  %2d.  %-20s  %.4f  <-- カフェイン自身", k, name, score);
    else
        logInfo("  %2d.  %-20s  %.4f", k, name, score);
    end
end
%[text] **✏️ やってみよう 4 — 上位ヒットを確認しましょう**
%[text] （カフェイン自身を除く）第 1 位はどれですか？
%[text] その分子はデータベース CSV に含まれていますか？
%[text] 期待値: テオブロミン（カカオ/チョコレートに含まれる）が上位に来るはずです。
%[text] テオブロミンは同じキサンチン骨格を持つが、N-メチル基が 2 つしかない。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: 上位ヒットを理解する
%[text] 最も類似した上位 5 分子（自己一致を除く）のサマリーテーブルを作ります。
%[text] カフェイン自身のインデックスをランクリストで特定し、検索結果の「上位5件（自分自身を除く）」を抽出。データ分析がしやすいよう MATLAB の `table` 形式に整形して詳しく観察してみましょう。構造の類似性（タニモトスコア）が、生物学的な分類（カテゴリ）や実際の用途とどのように連動しているか注目してください。
%[text] まず、カフェイン自身のインデックスをランクリストから特定して除外します。
selfIdx = find(result.Scores >= 1.0, 1);
%[text] 自己一致を除く上位 5 件を作成
topK = 5;
topRows = struct("Name", {}, "SMILES", {}, "Tanimoto", {}, "Category", {});
count = 0;
for k = 1:numel(result.Scores)
    if k == selfIdx
        continue;
    end
    count = count + 1;
    idx = result.Indices(k);
    topRows(count).Name     = names_valid(idx);
    topRows(count).SMILES   = smiles_valid(idx);
    topRows(count).Tanimoto = result.Scores(k);
    % 元のテーブルからカテゴリを検索
    match = tbl.CommonName == names_valid(idx);
    topRows(count).Category = tbl.Category(find(match, 1));
    if count >= topK; break; end
end
%[text] MATLAB テーブルとして表示
topTbl = struct2table(topRows);
disp("上位 5 件のカフェイン類縁化合物:");
disp(topTbl);
%[text] **化学的洞察**
%[text] メチルキサンチン類（カフェイン、テオブロミン、テオフィリン）はプリン二環式骨格を共有しています。
%[text] タニモト \> 0.5 はこの共有スキャフォールドを反映しています。
%[text] 非キサンチン分子のスコアはくっきり低く（\< 0.2）なります。
%%
%[text] ## セクション 6: 上位 3 件の類縁体を描く
%[text] 最も類似した分子を並べて可視化します。
%[text] 数値（タニモトスコア）の背後にある「実際の形の違い」を目で見て確認しましょう。スコアが高かった上位3件の分子構造を 2D 描画し、カフェインの構造（セクション1で描画したもの）と見比べます。どの官能基が共通していて、どこが変化しているでしょうか。
topDraw = min(3, numel(topRows));
figure("Name", "カフェインの上位類縁体", "Position", [100 100 1200 400]);
for k = 1:topDraw
    subplot(1, topDraw, k);
    mol_hit = emk.mol.fromSmiles(topRows(k).SMILES);
    titleStr = sprintf("%s\n(タニモト = %.3f)", topRows(k).Name, topRows(k).Tanimoto);
    emk.viz.draw2d(mol_hit, Title=titleStr);
end
logInfo("上位 %d 件の構造を 1 つの図に表示した。", topDraw);
%[text] **✏️ やってみよう 5 — 構造の違いを調べましょう**
%[text] セクション 1 のカフェインと上のテオブロミンの構造を比較してみましょう。
%[text] 両者の構造的な違いは何でしょうか？各分子の N-メチル基の数を数えてみましょう。
%[text] 期待値: カフェインは N-CH3 基が 3 個（キサンチン骨格の 1、3、7 位）。
%[text] テオブロミンは N-CH3 が 2 個（3 位と 7 位）；N1 は NH。
%[text] この 1 か所の違いがタニモトスコアを 1.0 未満に下げます。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 7: カフェインと上位 4 類縁体の類似度ヒートマップ
%[text] カフェインとその最近傍 4 件をヒートマップで比較します。
%[text] ここまでは「カフェイン vs 他の分子」という1対多の比較でした。しかし、ヒットした「類縁体同士」はどのくらい似ているのでしょうか？カフェインを含む上位5分子をすべて互いに交差比較した総当たり類似度マトリクス（Distance Matrix）を作成し、ヒートマップとして可視化します。これにより、似た分子同士が形成する「化学的クラスター（集団）」を視覚的に特定できます。
nHeat   = min(4, numel(topRows));
heatSmi = [CAFFEINE_SMILES, arrayfun(@(r) r.SMILES, topRows(1:nHeat), ...
           "UniformOutput", false)];
heatNames = [CAFFEINE_NAME, arrayfun(@(r) char(r.Name), topRows(1:nHeat), ...
             "UniformOutput", false)];
%[text] 解析とフィンガープリント計算
heatFps = cell(1, numel(heatSmi));
for k = 1:numel(heatSmi)
    if iscell(heatSmi)
        smi = heatSmi{k};
    else
        smi = heatSmi(k);
    end
    heatFps{k} = emk.fingerprint.morgan(emk.mol.fromSmiles(smi));
end

S = emk.similarity.matrix(heatFps);

figure("Name", "カフェイン類縁体ヒートマップ", "Position", [100 100 520 460]);
imagesc(S);
colormap("hot");
colorbar;
clim([0 1]);
axis square;
xticks(1:numel(heatNames));  xticklabels(heatNames);  xtickangle(30);
yticks(1:numel(heatNames));  yticklabels(heatNames);
title("タニモト類似度 -- カフェインと最近傍類縁体");

for r = 1:numel(heatNames)
    for c = 1:numel(heatNames)
        clr = "black";
        if S(r, c) < 0.5
            clr = "white";
        end
        text(c, r, sprintf("%.2f", S(r, c)), ...
            "HorizontalAlignment", "center", "FontSize", 9, "Color", clr);
    end
end
%[text] 対角線（自己類似度 = 1.0）の明るさと、カフェイン/テオブロミン間の高いスコア（0.53）に注目してみましょう。キサンチンクラスター（カフェイン/テオブロミン/...）は互いに高い類似度を示しています。
%%
%[text] ## 演習
%[text] E1: Dice 類似度でランク付けする（タニモトの代わりに）。
%[text]     順序は変わるか？なぜそうなるか/ならないか？
%[text]     ヒント: `emk.similarity.rankBy(..., Metric="dice")` を使う
result_dice = emk.similarity.rankBy(fp_caffeine, fps_valid, Inf, Metric="dice");
% result_dice.Scores(1:5) と result.Scores(1:5) を比較する
%[text] 
%[text] E2: MACCS キーフィンガープリントを使って検索を繰り返す（Morgan の代わりに）。
%[text]     テオブロミンは依然として第 1 位か？
%[text]     ヒント: `emk.fingerprint.morgan` を `emk.fingerprint.maccs` に置き換える。
fp_caf_maccs = emk.fingerprint.maccs(mol_caffeine);
fps_maccs = cell(1, numel(fps_valid));
for i = 1:numel(fps_valid)
    mol_i = emk.mol.fromSmiles(smiles_valid(i));
    fps_maccs{i} = emk.fingerprint.maccs(mol_i);
end
result_maccs = emk.similarity.rankBy(fp_caf_maccs, fps_maccs);
%[text] 
%[text] E3: データベースでカフェインに最も類似しない分子を探す。
%[text]     タニモトスコアはいくつか？そのカテゴリは？
%[text]     ヒント: result.Scores は降順ソートなので、最後のエントリが最小値。
lastIdx   = result.Indices(end);
lastName  = names_valid(lastIdx);
lastScore = result.Scores(end);
logInfo("最も類似しない: %s (タニモト = %.4f)", lastName, lastScore);
% tbl でそのカテゴリを調べる

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
