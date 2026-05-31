%[text] # S02: 薬物フィルタ — リピンスキーの「ルール・オブ・ファイブ」
%[text] EasyMolKit 応用ストーリー — レイヤー 2
%[text] 
%[text] 私たちが飲む「薬」は、体内の標的（タンパク質など）に届くだけでなく、腸から吸収されて血液に乗り、適切に全身を巡る必要があります。
%[text] どんなに実験室で強力に効く分子が見つかっても、体に吸収されなければ薬にはなり得ません。 1997年、製薬会社ファイザーの科学者クリストファー・リピンスキーらは、経口医薬品（飲み薬）の多くが共通して持つ「4つの物理化学的な特徴」を発見しました。これが有名な「リピンスキーの規則（Rule of 5）」です。 
%[text] このスクリプトでは、ChEMBL から取得した FDA 承認薬 200 種のデータベースに対してリピンスキーの規則を適用し、どの分子が「医薬品らしい経口吸収性」を備えているかを MATLAB のデータ処理テクニックを使って一括で判定します。 
%[text] ## 学習目標
%[text] - 記述子計算機能（`emk.descriptor.calculate`）を使って分子の性質を網羅的に取得する
%[text] - リピンスキーの規則（Rule of 5）の基準をコードで実装する
%[text] - MATLABの条件判定とテーブル操作（論理マスキング）を使って、条件を満たす分子を抽出する
%[text] - 独自の「合格スコア」を算出し、データベース全体をスクリーニングする \
%[text] ## 前提条件
%[text] - F02（特性計算）の完了
%[text] - RDKitインストール済み（`emk.setup.install()` を一度だけ実行しておく）
%[text] - 追加Toolbox不要（MATLAB だけで動きます） \
%[text] **所要時間**: 15〜20 分 | 実行方法: Ctrl+Enter でセクションを一つずつ実行
%[text] **データ**
%[text] - `data/list/fda_drugs.csv` — FDA 承認薬 200 種（ChEMBL, CC-BY-SA 3.0） \
%[text] **参考文献**
%[text] - Lipinski CA et al. (1997) *Adv Drug Deliv Rev* 23:3-25. doi:10.1016/S0169-409X(96)00423-1 〔要機関アクセス〕\
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
logInfo("S02: セットアップ完了");
%%
%[text] ## セクション 1: リピンスキーの規則（Rule of 5）とは？
%[text] リピンスキーの規則は、経口摂取したときに「生体膜を透過してうまく吸収されやすいか」を見分けるためのガイドラインです。すべての基準に「5」またはその倍数の数字が登場するため、\*\*Rule of 5（Ro5）\*\* と呼ばれています。
%[text] **【リピンスキーの4つの規則】**
%[text] 1. 分子量（Molecular Weight）:$&dollar&;\\leq 500 \\text{ g/mol}&dollar&;$ 以下（分子が大きすぎると膜を通り抜けられません）
%[text] 2. 脂溶性（LogP）: $&dollar&;5&dollar&;$ 以下（油に溶けすぎると膜に留まってしまい、血液に抜けられません）
%[text] 3. 水素結合供与体数（H-Bond Donors）:$&dollar&;5&dollar&;$ 個以下（水分子と強く結びつきすぎると、脂質の膜を嫌います）
%[text] 4. 水素結合受容体数（H-Bond Acceptors）: $&dollar&;10&dollar&;$ 個以下（同上の理由で、多すぎると膜を透過できません） \
%[text] ※一般的に、この4つのうち3つ以上を満たしていれば「吸収性が良好な分子（医薬品らしい候補）」とみなされます。
%[text] 今回はEasyMolKitを使って、これらの値を一括計算するため、世界で最も消費される経口鎮痛薬アスピリンで確認してみましょう。
ASPIRIN_SMILES = "CC(=O)Oc1ccccc1C(=O)O";
ASPIRIN_NAME   = "Aspirin";

mol_aspirin = emk.mol.fromSmiles(ASPIRIN_SMILES);
desc_aspirin = emk.descriptor.calculate(mol_aspirin, ...
    ["MolWt", "LogP", "NumHDonors", "NumHAcceptors", "TPSA"]);

logInfo("--- %s（アセチルサリチル酸）---", ASPIRIN_NAME);
logInfo("  MW              : %.2f Da  （規則: <= 500）", desc_aspirin.MolWt);
logInfo("  LogP            : %.2f     （規則: <= 5  ）", desc_aspirin.LogP);
logInfo("  水素結合供与体  : %d       （規則: <= 5  ）", desc_aspirin.NumHDonors);
logInfo("  水素結合受容体  : %d       （規則: <= 10 ）", desc_aspirin.NumHAcceptors);
logInfo("  TPSA            : %.1f A^2 （経口 < 130 A^2 の目安）", desc_aspirin.TPSA);

figure("Name", "アスピリン", "Position", [100 100 440 380]);
emk.viz.draw2d(mol_aspirin, Title="アスピリン（アセチルサリチル酸）");

%[text] **✏️ やってみよう 1 — 別の鎮痛薬で確認してみましょう**
%[text] 上の数値を見ながら、アスピリンが 4 つの Ro5 基準をいくつ違反しているか数えてみてください。
%[text] 期待値: 違反 0 件。MW ~180、LogP ~1.3、HBD=1、HBA=3 — すべて通過。
%[text] 次に、イブプロフェン（`"CC(C)Cc1ccc(cc1)C(C)C(=O)O"`）も試してみましょう。
%[text]   mol\_ibu  = `emk.mol.fromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O")`;
%[text]   desc\_ibu = `emk.descriptor.calculate(mol_ibu, ["MolWt","LogP","NumHDonors","NumHAcceptors"])`;
%[text] 期待値: MW ~206、LogP ~3.1、HBD=1、HBA=1 — どちらも Ro5 通過。
%[text] Q: イブプロフェンはアスピリンより LogP が高いです。化学構造のどこに違いがあると思いますか？
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: FDA 薬物データベースを読み込む
%[text] ChEMBL（CC-BY-SA 3.0）から取得した FDA 承認薬 200 種のデータを読み込みます。
%[text] CSV には ChEMBL で事前計算された記述子（ALogP 等）が含まれています。
%[text] `emk.filter.lipinski` が期待する列名に合わせてリネームします。
dataFile = fullfile(projectRoot, "data", "list", "fda_drugs.csv");
rawTbl = readtable(dataFile, "TextType", "string");
logInfo("ChEMBL から FDA 承認薬 %d 種を読み込んだ。", height(rawTbl));
%[text] `emk.filter.lipinski` の必要な列名に合わせてリネーム
drugTbl = renamevars(rawTbl, ...
    ["MolecularWeight", "ALogP",  "HBondDonors",  "HBondAcceptors"], ...
    ["MolWt",           "LogP",   "NumHDonors",   "NumHAcceptors"]);
logInfo("記述子サマリー（フィルタ前）:");
logInfo("  MW   -- 最小: %5.1f  中央値: %5.1f  最大: %6.1f  Da", ...
    min(drugTbl.MolWt), median(drugTbl.MolWt), max(drugTbl.MolWt));
logInfo("  LogP -- 最小: %5.2f  中央値: %5.2f  最大: %5.2f", ...
    min(drugTbl.LogP), median(drugTbl.LogP), max(drugTbl.LogP));
logInfo("  HBD  -- 最小: %d     中央値: %g    最大: %d", ...
    min(drugTbl.NumHDonors), median(drugTbl.NumHDonors), max(drugTbl.NumHDonors));
logInfo("  HBA  -- 最小: %d     中央値: %g    最大: %d", ...
    min(drugTbl.NumHAcceptors), median(drugTbl.NumHAcceptors), max(drugTbl.NumHAcceptors));
%[text] **✏️ やってみよう 2 — 最重量の薬物を探してみましょう**
%[text] 読み込んだ 200 種の中で、最も分子量が大きい薬物を探してみましょう。MATLAB の「条件に合う行だけを抜き出す機能」を使って、その薬物の名前と分子量を確認してみましょう。
%[text] ヒント: `drugTbl(drugTbl.MolWt == max(drugTbl.MolWt), :)` で最重量の行を抽出できます。
%[text] 期待値: MW が 900 Da を超える大型天然物由来の薬物（ポリエン系抗真菌薬など）。
%[text] Q: その薬物は経口投与ですか？注射剤ですか？Ro5 違反数と一致するか確認してみましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: リピンスキーのルール・オブ・ファイブを適用する
%[text] ### emk.filter.lipinski
%[text] `emk.filter.lipinski` は記述子テーブルに 2 列を追加します。
%[text] - `Pass_Ro5` (logical) — 違反数が `MaxViolations` 以下なら `true`
%[text] - `Violations_Ro5` (double) — 違反した基準の数（0〜4） \
%[text] まず違反ゼロの厳格フィルタ（デフォルト）から適用してみましょう。
drugTbl = emk.filter.lipinski(drugTbl);          % MaxViolations=0（デフォルト）

nTotal  = height(drugTbl);
nPass   = sum(drugTbl.Pass_Ro5);
nFail   = nTotal - nPass;

logInfo("Ro5 フィルタ結果（厳格、MaxViolations=0）:");
logInfo("  通過 : %d / %d  (%.0f%%)", nPass,  nTotal, 100*nPass/nTotal);
logInfo("  失敗 : %d / %d  (%.0f%%)", nFail,  nTotal, 100*nFail/nTotal);

%[text] 各基準の違反回数を数える
vMW  = sum(drugTbl.MolWt          > 500);
vLP  = sum(drugTbl.LogP           > 5  );
vHBD = sum(drugTbl.NumHDonors     > 5  );
vHBA = sum(drugTbl.NumHAcceptors  > 10 );

logInfo("各基準の違反件数:");
logInfo("  MW  > 500: %d 件", vMW);
logInfo("  LogP > 5 : %d 件", vLP);
logInfo("  HBD > 5  : %d 件", vHBD);
logInfo("  HBA > 10 : %d 件", vHBA);

%[text] **✏️ やってみよう 3 — どの基準が最もよく破られますか？**
%[text] 上に表示された vMW、vLP、vHBD、vHBA の値を比較してみてください。
%[text] Q: このデータセットで最も違反件数が多い Ro5 基準はどれですか？
%[text] 期待値: 天然物由来の薬物は大型で極性が高いため、MW と HBA の違反が多い傾向があります。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: 違反数の分布
%[text] 見つかった脱落分子が、なぜリピンスキーの規則に引っかかってしまったのか、具体的な数値を表示して原因を突き止めましょう。
violCounts = zeros(1, 5);
for v = 0:4
    violCounts(v+1) = sum(drugTbl.Violations_Ro5 == v);
end

figure("Name", "Ro5 違反数分布", "Color", "white", "Position", [100 100 560 400]);
bar(0:4, violCounts, "FaceColor", [0.2 0.6 0.9]);
xlabel("Ro5 違反数");
ylabel("FDA 薬物数");
title("リピンスキー Ro5 -- 違反数分布（FDA 薬物 200 種）");
xticks(0:4);
xticklabels({"0（通過）", "1", "2", "3", "4（全違反）"});
grid("on");
for v = 0:4
    if violCounts(v+1) > 0
        text(v, violCounts(v+1) + 0.5, sprintf("%d", violCounts(v+1)), ...
            "HorizontalAlignment", "center", "FontWeight", "bold");
    end
end
logInfo("違反数分布をプロットした。");
%[text] **✏️ やってみよう 4 — 緩和フィルタを試してみましょう**
%[text] 棒グラフを見て、違反ちょうど 1 件の薬物が何種あるか確認してください。
%[text] ヒント: `violCounts(2)` で違反 1 件の数が分かります。
%[text] Q: `MaxViolations=1` に緩和すると、何種の薬物が新たに通過しますか？
%[text] ヒント: `nPass + violCounts(2)` で緩和フィルタ通過の合計を求められます。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: 化学空間の可視化 -- MW vs LogP
%[text] ### 薬物様化学空間とは？
%[text] MW（y 軸）と LogP（x 軸）をプロットすると、Ro5 の境界線で囲まれた「薬物様」領域が浮かび上がります。経口薬の多くは MW=500 Da ライン以下、LogP=5 ライン左側に集まります。境界の外にある薬物は、注射剤・局所適用・プロドラッグであることがほとんどです。
figure("Name", "MW vs LogP（薬物空間）", "Color", "white", "Position", [100 100 560 440]);
hold on;
%[text] 通過した薬物（青い丸）
passIdx = drugTbl.Pass_Ro5;
scatter(drugTbl.LogP(passIdx),  drugTbl.MolWt(passIdx),  40, ...
    [0.15 0.55 0.85], "o", "filled", "MarkerFaceAlpha", 0.6);
%[text] 失敗した薬物（オレンジの三角）
scatter(drugTbl.LogP(~passIdx), drugTbl.MolWt(~passIdx), 60, ...
    [0.9 0.45 0.1], "^", "filled", "MarkerFaceAlpha", 0.7);
%[text] Ro5 境界線
xLim = xlim;  yLim = ylim;
plot([5 5], [0 max(yLim(2), 600)], "r--", "LineWidth", 1.5);   % LogP = 5
yline(500, "r--", "LineWidth", 1.5);                           % MW = 500
%[text] アスピリンの注釈
aspIdx = find(strcmpi(drugTbl.Name, "ASPIRIN"), 1);
if ~isempty(aspIdx)
    scatter(drugTbl.LogP(aspIdx), drugTbl.MolWt(aspIdx), 120, ...
        "g", "p", "filled");
    text(drugTbl.LogP(aspIdx) + 0.1, drugTbl.MolWt(aspIdx) + 10, ...
        "Aspirin", "FontSize", 9, "Color", [0 0.5 0]);
end

hold off;
legend({"Ro5 通過", "Ro5 失敗", "LogP=5 限界", "MW=500 限界", "アスピリン"}, ...
    "Location", "northwest");
xlabel("LogP  （親油性）");
ylabel("分子量（Da）");
title("FDA 薬物: MW vs LogP — リピンスキー空間");
grid("on");

logInfo("化学空間散布図を表示した。");
logInfo("  青い領域（MW<500、LogP<5）: 典型的な経口薬物様空間");
logInfo("  境界外の薬物: 注射剤・局所適用・プロドラッグが多い");
%[text] **✏️ やってみよう 5 — 散布図を読んでみましょう**
%[text] 生成された散布図を眺めてください。
%[text] Q1: 青い点（Ro5 通過）の多くは左下象限（MW\<500、LogP\<5）に集まっていますか？
%[text] Q2: LogP~5、MW~480 付近に通過薬物のクラスターが見えますか？（Ro5 ぎりぎり通過ライン）
%[text] Q3: 右上に離れているオレンジの点（Ro5 不合格）はどんな薬物だと思いますか？
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: ルールを破る薬物 -- Ro5 の例外
%[text] 
%[text] ### Ro5 はガイドラインであり、絶対ルールではありません
%[text] Ro5 は受動拡散で吸収される小分子経口薬向けのルールです。
%[text] 一部の承認済み薬物は次の理由で Ro5 に違反しても経口吸収されます。
%[text] - 能動輸送（キャリアタンパク）によって吸収される
%[text] - プロドラッグとして投与され、体内で活性体に変換される \
%[text] 有名な例:
%[text] - シクロスポリン A — MW ~1202 Da、3 違反（免疫抑制薬）
%[text] - リファンピシン — MW ~823 Da、2 違反（抗生物質）
%[text] - アジスロマイシン — MW ~749 Da、2 違反（マクロライド系抗生物質）
%[text] - アトルバスタチン — MW ~559 Da、1 違反（コレステロール薬） \
failTbl = drugTbl(~drugTbl.Pass_Ro5, :);
logInfo("--- 厳格な Ro5 を失敗した %d 種の FDA 薬物 ---", height(failTbl));
%[text] 違反数（降順）次に MW でソート
failTbl = sortrows(failTbl, ["Violations_Ro5", "MolWt"], ["descend", "descend"]);
%[text] 上位 10 件の最悪違反者を表示
nShow = min(10, height(failTbl));
logInfo("上位 %d 件の Ro5 最悪違反者（承認済み経口/全身薬物）:", nShow);
disp(failTbl(1:nShow, ["Name","MolWt","LogP","NumHDonors","NumHAcceptors","Violations_Ro5"]));
%[text] **✏️ やってみよう 6 — Ro5 違反薬物の構造を見てみましょう**
%[text] 上の一覧から気になる薬物を 1 つ選んで、構造を描いてみてください。
%[text] ヒント:
%[text]   idx = find(strcmpi(drugTbl.Name, "RIFAMPICIN"), 1);   % 任意の薬物名に変えられます
%[text]   mol\_big = `emk.mol.fromSmiles(drugTbl.SMILES(idx))`;
%[text]   figure("Name", drugTbl.Name(idx));
%[text]  `emk.viz.draw2d(mol_big, Title=drugTbl.Name(idx))`;
%[text] 構造の中に水素結合供与体（OH、NH）や受容体（C=O、N）がいくつ見えますか？
%[text] それが Ro5 違反数と一致しているか確認してみましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 7: 緩和フィルタ（MaxViolations=1）
%[text] ### 「ルール・オブ・ファイブ＋α」
%[text] 実際の多くの創薬プログラムでは `MaxViolations=1`（1 基準の違反を許容）を使います。
%[text] その理由は次のとおりです。
%[text] - Ro5 は受動拡散薬物という特定のデータセットに基づいている
%[text] - 能動輸送が高い極性や大きなサイズを補うことがある
%[text] - 最適化ステップで 1 つの違反を後から修正できることもある \
drugTblRelaxed = emk.filter.lipinski(drugTbl, MaxViolations=1);
nPassRelaxed   = sum(drugTblRelaxed.Pass_Ro5);

logInfo("緩和フィルタ（MaxViolations=1）: %d / %d 通過  (%.0f%%)", ...
    nPassRelaxed, nTotal, 100*nPassRelaxed/nTotal);
logInfo("追加で通過する薬物: %d 種", nPassRelaxed - nPass);
%[text] 新たに通過した薬物を特定（厳格フィルタ失敗 → 緩和フィルタ通過）
isNewPass = ~drugTbl.Pass_Ro5 & drugTblRelaxed.Pass_Ro5;
newPassTbl = drugTbl(isNewPass, :);
newPassTbl = sortrows(newPassTbl, "MolWt", "descend");

logInfo("新たに通過した薬物のサンプル（各 1 違反）:");
nSample = min(5, height(newPassTbl));
disp(newPassTbl(1:nSample, ["Name","MolWt","LogP","NumHDonors","NumHAcceptors","Violations_Ro5"]));

logInfo("--- S02 完了 ---");
logInfo("サマリー:");
logInfo("  厳格 Ro5（0 違反）: %d / %d 通過（%.0f%%）", ...
    nPass, nTotal, 100*nPass/nTotal);
logInfo("  緩和 Ro5（1 違反）: %d / %d 通過（%.0f%%）", ...
    nPassRelaxed, nTotal, 100*nPassRelaxed/nTotal);
logInfo("  最も多く違反される基準: セクション 3〜4 で vMW、vLP、vHBD、vHBA を比較する");
%[text] **✏️ やってみよう 7 — MaxViolations=2 まで広げてみましょう**
%[text] `MaxViolations=2`（2 基準まで許容）でフィルタをかけ直して、新たに通過する薬物を確認してください。
%[text] ヒント:
%[text]   tbl2 = `emk.filter.lipinski(drugTbl, MaxViolations=2)`;
%[text]   pass2 = tbl2(tbl2.Pass\_Ro5 & ~drugTblRelaxed.Pass\_Ro5, :);
%[text]   disp(pass2(:, \["Name","MolWt","LogP","NumHDonors","NumHAcceptors","Violations\_Ro5"\]))
%[text] Q: 追加された薬物の中に、大型マクロライド系抗生物質や免疫抑制薬は含まれていますか？
% ... （ここにコードを書いてみましょう）

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
