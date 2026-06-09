%[text] # A10: リード最適化 — 多目的プロパティ最適化
%[text] EasyMolKit アナリティクス — レイヤー 3
%[text]
%[text] **ストーリー**
%[text] ある創薬化学者が、新しい抗炎症鎮痛薬を開発するための構造リードとしてアスピリンを選定しました。アスピリンは有効ですが、分子量（180 g/mol）が現代の医薬品候補としては低く、特異性を高めるための置換基を追加する余地があります。
%[text] 課題: ターゲット結合を改善するために官能基を追加すると、複数の物理化学的性質が同時に変化することがあります。
%[text] - 親油性基を追加すると LogP（オクタノール/水分配係数の常用対数）が上がり、膜透過性が向上しますが、分子量（MW）も上がり水溶性が低下します。
%[text] - 溶解性を改善するために極性基を追加すると、位相的極性表面積（TPSA）が上がり、経口吸収性が低下する可能性があります。
%[text] これらのトレードオフがリード最適化の「多目的」な性質を定義します。
%[text] 単一の分子が全ての望ましい性質を同時に最適化することは難しいため、化学者は「最も悪くない」妥協点を見つける必要があります。
%[text] この演習では以下を行います:
%[text]
%[text] 1. 3つの経口薬設計ターゲットに対してアスピリンを評価します:
%[text] LogP = 2.0、TPSA = 70 Å²、MW = 350 g/mol。
%[text] 2. 200のFDA承認薬をアナログ参照ライブラリとしてロードし、その3Dプロパティ空間を探索します。
%[text] 3. Derringer-Suich法を使って生のプロパティを無次元望ましさスコアに変換します。
%[text] 4. fgoalattain（Optimization Toolbox）で目標達成問題を定式化・解き、全ターゲットを同時に最もよく満たすプロパティ空間のブレンドを求めます。
%[text] 5. epsilon制約法（linprog）でLogPとTPSAのトレードオフのパレートフロントを追跡します。
%[text] 6. 上位のリード最適化候補を選択し、可視化します。
%[text]
%[text] **学習目標**
%[text] - なぜリード最適化が本質的に多目的なのかを理解する
%[text] - Derringer-Suich望ましさとして薬理学的直感をエンコードする
%[text] - fgoalattain（Optimization Toolbox）で全プロパティゴールにわたる最大スラックを最小化する
%[text] - 達成因子gammaを進捗指標として解釈する
%[text] - epsilon制約法とlinprogを適用してパレートフロントを追跡する
%[text] - パレート最適候補と劣位候補を区別する
%[text]
%[text] **前提条件**
%[text] - S02（ドラッグフィルター）修了 — LipinskiのRo5のコンテキスト
%[text] - 推奨: A03（QSAR回帰）QSARモデリングのコンテキスト
%[text] - Optimization Toolbox（fgoalattain、linprog）— 利用できない場合、セクション4〜5は手動フォールバックを使用するため全ての概念は学習可能。
%[text]
%[text] **動作環境**
%[text] 互換性サマリー:
%[text] 学生ライセンス / キャンパスワイド — フルOptimization Toolbox
%[text] MATLAB Online 無料版 — Optimization Toolbox 非対応
%[text] MATLAB Online（個人 / キャンパス） — Optimization Toolbox 利用可
%[text] ツールボックスがない場合、セクション4〜5は網羅的探索による手動代替に切り替わります。
%[text]
%[text] 推定所要時間: 45〜90分
%[text]
%[text] **データ:**
%[text] data/list/fda_drugs.csv — 200 FDA承認薬（ChEMBL、CC-BY-SA 3.0）
%[text]
%[text] **参考文献**
%[text] Derringer G (1980) Simultaneous optimization of several response variables. J Quality Technology 12:214-219. doi:10.1080/00224065.1980.11980968
%[text] Veber DF et al. (2002) Molecular properties that influence the oral bioavailability of drug candidates. J Med Chem 45:2615-2623. doi:10.1021/jm020017n
%[text] Lipinski CA et al. (2001) Experimental and computational approaches to estimate solubility and permeability in drug discovery and development. Adv Drug Deliv Rev 46:3-26. doi:10.1016/S0169-409X(00)00129-0
%[text] Cohon JL & Marks DH (1975) A review and evaluation of multiobjective programming techniques. Water Resour Res 11:208-220. doi:10.1029/WR011i002p00208 (epsilon制約法)
%[text] Charnes A & Cooper WW (1977) Goal programming and multiple objective optimization. Eur J Oper Res 1:39-54. doi:10.1016/S0377-2217(77)81007-2 (目標達成の基礎)
%[text]
%[text] 実行方法: Ctrl+Enter でセクションを1つずつ実行
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

%[text] メインの実行前に、Python/RDKit プロセスをウォームアップします
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;

%[text] Optimization Toolbox が使用可能かどうかを確認します
hasOptTbx = license("test", "optimization_toolbox");
if hasOptTbx
    logInfo("A10: Optimization Toolbox を検出 -- フル最適化が有効。");
else
    logWarn("A10: Optimization Toolbox が検出されない。");
    logWarn("     セクション 4〜5 は手動グリッドサーチフォールバックを使用。");
    logWarn("     fgoalattain / linprog には Optimization Toolbox をインストール。");
end

logSection("A10", "セクション 0: セットアップ", "アナリティクス L3");
logInfo("A10: セットアップ完了。");
%%
%[text] ## セクション 1: リード化合物と多目的ジレンマ
%[text]
%[text] セットアップが完了しました。まず、リード化合物であるアスピリンの特性を計算し、目標からどの程度ずれているかを確認しましょう。
%[text] LogP（親油性の常用対数）、TPSA（位相的極性表面積）、MW（分子量）の3軸で多目的トレードオフの構造を理解しましょう。
%[text]
%[text] ### コンセプト: 経口薬らしさの3つの柱
%[text]
%[text] 現代の経口薬創薬では、3つの主要な物理化学的閾値に基づいて分子を評価します（Veber 2002 および Lipinski 2001 に基づく）。
%[text]
%[text] (1)  LogP（親油性の常用対数）
%[text] 範囲: 1.0 <= LogP <= 3.5（膜透過性と水溶性のバランスに最適）
%[text] 低すぎる場合: 膜透過性が低い
%[text] 高すぎる場合: 水溶性が低く、CYP 代謝リスクがある
%[text]
%[text] (2)  TPSA（位相的極性表面積、Å²）
%[text] ターゲット: <= 90 Å²（良好な経口吸収の Veber ルール）
%[text] ターゲット: <= 60 Å²（中枢神経系透過 -- 血液脳関門）
%[text] 高すぎる場合: TPSA > 130 Å² -> 一般に受動吸収が低い
%[text]
%[text] (3)  MW（分子量、g/mol）
%[text] ターゲット: <= 400 g/mol（より厳格な「beyond-Ro5」ガイドライン）
%[text] Lipinski: <= 500 g/mol
%[text] 高すぎる場合: 吸収が低下し、排出トランスポーターに認識される可能性がある
%[text]
%[text] 目標達成ターゲット（セクション 3〜5 で使用）:
%[text] LogP* = 2.0、TPSA* = 70.0 Å²、MW* = 350.0 g/mol
%[text]
%[text] リード化合物
logSection("A10", "セクション 1: リード化合物と多目的ジレンマ", "アナリティクス L3");
LEAD_SMILES = "CC(=O)Oc1ccccc1C(=O)O";   % アスピリン
LEAD_NAME   = "Aspirin";

%[text] プロパティターゲット
TARGET_LOGP = 2.0;    % 最適親油性
TARGET_TPSA = 70.0;   % A^2、バランスのとれた極性表面積
TARGET_MW   = 350.0;  % g/mol、厳格な薬らしさ

mol_lead  = emk.mol.fromSmiles(LEAD_SMILES);
d_lead    = emk.descriptor.calculate(mol_lead);

logInfo("リード化合物: %s", LEAD_NAME);
logInfo("  LogP = %+.2f   ターゲット = %.1f  偏差 = %+.2f", ...
    d_lead.LogP, TARGET_LOGP, d_lead.LogP - TARGET_LOGP);
logInfo("  TPSA = %5.1f A^2   ターゲット = %.1f  偏差 = %+.1f", ...
    d_lead.TPSA, TARGET_TPSA, d_lead.TPSA - TARGET_TPSA);
logInfo("  MW   = %5.1f g/mol  ターゲット = %.1f  偏差 = %+.1f", ...
    d_lead.MolWt, TARGET_MW, d_lead.MolWt - TARGET_MW);

%[text] **💡 観察ポイント 1**
%[text] アスピリンの MW（180 g/mol）は 350 g/mol のターゲットをかなり下回っています。
%[text] これは、アナログを開発する創薬化学者にとってどのような意味を持つでしょうか？
%[text] （ヒント: 分子量が低いリードは、置換基を追加する「余地」が多いです。）
%[text] イブプロフェンの特性を計算してアスピリンと比較しましょう:
%[text] mol_ibu = emk.mol.fromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
%[text] d_ibu   = emk.descriptor.calculate(mol_ibu);
%[text] どちらの分子が3つのターゲット全てに同時に近いかを確認しましょう。
%[text] ここでのターゲットは仮想的な中枢神経系（CNS）薬のためのものです。
%[text] 末梢性抗炎症薬のために TARGET_TPSA をどのように変更するべきでしょうか？
%[text] （ヒント: CNS 薬は TPSA < 60 が必要ですが、末梢薬は < 90 を許容できます。）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: 参照ライブラリの構築とプロパティ計算
%[text]
%[text] リード化合物のプロパティを確認しました。次に、FDA 承認薬 200 種を参照ライブラリとしてロードし、化学空間を探索します。
%[text] 実際の薬がどのプロパティ空間に分布するかを確認しましょう。
%[text]
%[text] ### コンセプト: FDA 薬を「参照化学空間」として使用
%[text]
%[text] リード最適化において、参照ライブラリは次の2つの目的を果たします:
%[text] (a)  現実的にアクセス可能な化学空間を定義すること --
%[text] 類似の薬が既に存在し、臨床試験を通過していれば、近接するアナログが成功する可能性が高まります。
%[text] (b)  プロパティ最適化の「ターゲット分布」を提供すること:
%[text] ライブラリは、実際の薬が達成する LogP（分配係数の常用対数）/ TPSA（位相的極性表面積）/ MW（分子量）値を示します。
%[text]
%[text] ここでは、ChEMBL からの 200 FDA 承認薬を参照として使用します。
%[text] 実際のプロジェクトでは、商業的に入手可能なビルディングブロックライブラリや計算的に生成されたアナログのセットを使用します。
logSection("A10", "セクション 2: 参照ライブラリの構築とプロパティ計算", "アナリティクス L3");
DATA_FILE = "data/list/fda_drugs.csv";
logInfo("%s から %d エントリをロード", DATA_FILE, height(rawTbl));

%[text] 分子をパースします。
nLib  = height(rawTbl);
mols  = cell(1, nLib);
valid = false(1, nLib);

for k = 1:nLib
    try
        mols{k} = emk.mol.fromSmiles(rawTbl.SMILES(k));
        valid(k) = true;
    catch
        % パース不可エントリをスキップ
    end
end

validIdx = find(valid);
mols     = mols(validIdx);
libNames = rawTbl.Name(validIdx);
libSmiles = rawTbl.SMILES(validIdx);
nValid   = numel(mols);
logInfo("%d / %d 分子をパースしました。", nValid, nLib);

%[text] 記述子を計算します: LogP、TPSA、MW
DESCS   = ["LogP", "TPSA", "MolWt"];
descTbl = emk.descriptor.batchCalculate(mols, DESCS);

logp_vec = descTbl.LogP;    % N x 1 double
tpsa_vec = descTbl.TPSA;    % N x 1 double
mw_vec   = descTbl.MolWt;   % N x 1 double

%[text] プロパティ行列: N x 3  [LogP、TPSA、MW]
propMat = [logp_vec, tpsa_vec, mw_vec];

logInfo("ライブラリプロパティ範囲:");
logInfo("  LogP : [%.2f, %.2f]  平均 = %.2f", ...
    min(logp_vec), max(logp_vec), mean(logp_vec));
logInfo("  TPSA : [%.1f, %.1f] A^2  平均 = %.1f", ...
    min(tpsa_vec), max(tpsa_vec), mean(tpsa_vec));
logInfo("  MW   : [%.1f, %.1f] g/mol  平均 = %.1f", ...
    min(mw_vec), max(mw_vec), mean(mw_vec));

%[text] --- 3D プロパティ空間プロット ---
figure("Name", "A10 Sec2: FDA 薬プロパティ空間（LogP vs TPSA vs MW）");
scatter3(logp_vec, tpsa_vec, mw_vec, 30, mw_vec, "filled", ...
    MarkerFaceAlpha=0.6);
hold on;
%[text] リード化合物をマークします。
scatter3(d_lead.LogP, d_lead.TPSA, d_lead.MolWt, 120, "r", ...
    "^", "filled", DisplayName=LEAD_NAME);
%[text] ターゲット点をマークします。
scatter3(TARGET_LOGP, TARGET_TPSA, TARGET_MW, 150, "k", ...
    "pentagram", LineWidth=2, DisplayName="Target");
xlabel("LogP");
ylabel("TPSA (A^2)");
zlabel("MW (g/mol)");
title("FDA 薬プロパティ空間");
cb = colorbar();
cb.Label.String = "MW (g/mol)";
legend(Location="best");
grid on;

%[text] **💡 観察ポイント 2**
%[text] 3D 散布図を回転させて、LogP と TPSA の間に目に見える相関があるか確認しましょう。
%[text] LogP と TPSA の間に目に見える相関はありますか？
%[text] （ヒント: 親油性化合物は極性基が少ない傾向があります。）
%[text] ピアソン相関を計算してみましょう: corr(logp_vec, tpsa_vec)
%[text] 3 つのターゲット全てを既に満たすライブラリ分子は何個あるか確認しましょう:
%[text] 1.0 <= LogP <= 3.5 AND TPSA <= 90 AND MW <= 400 の条件を満たす分子の数を数えます。
%[text] 数える: sum((logp_vec >= 1.0 & logp_vec <= 3.5) & ...
%[text] tpsa_vec <= 90 & mw_vec <= 400)
%[text] Lipinski の Ro5 境界を透明なボックスとして 3D プロットに追加してみましょう:
%[text] [LogP <= 5、TPSA（Ro5 基準ではない）、MW <= 500] の条件を考慮します。
%[text] Veber / 厳格な「350 g/mol」ガイドラインとどのように比較できるか確認しましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: 望ましさ関数 -- プロパティからスコアへ
%[text]
%[text] 参照ライブラリのプロパティ空間が明らかになりました。次に、Derringer-Suich 望ましさ関数を用いて各プロパティを [0,1] のスコアに変換します。
%[text] 複数のプロパティを同時に最適化し、単一の指標で表現する方法を学びましょう。
%[text]
%[text] ### コンセプト: Derringer-Suich 望ましさ（1980）
%[text]
%[text] 望ましさ関数は生のプロパティ値を [0, 1] の無次元スコア d_i にマッピングします。
%[text] d_i = 1.0  --  プロパティが理想ターゲットに完全に一致する場合
%[text] d_i = 0.0  --  プロパティが許容範囲外（下限未満または上限超過）の場合
%[text]
%[text] 「両側」（ターゲット値）目標の標準形は次の通りです。
%[text]
%[text] d_i(p) = 0                                p < L_i の場合
%[text] = ((p - L_i) / (T_i - L_i))^s     L_i <= p <= T_i の場合
%[text] = ((U_i - p) / (U_i - T_i))^t     T_i < p <= U_i の場合
%[text] = 0                                p > U_i の場合
%[text]
%[text] ここで、L_i は許容下限、T_i は理想ターゲット、U_i は許容上限、s と t は曲率を制御します。
%[text]
%[text] 複合望ましさ D は全ての個別 d_i を組み合わせて計算します。
%[text]
%[text] D = (d_1 * d_2 * ... * d_k)^(1/k)    （幾何平均）
%[text]
%[text] D = 1 の場合、全プロパティが同時に理想的です。
%[text] D = 0 の場合、少なくとも 1 つのプロパティが許容できません。
%[text] 幾何平均は、単一の非常に低いスコアに対して算術平均よりも強いペナルティを課します。
%[text]
%[text] プロパティ望ましさの定義（この演習）:
%[text]
%[text] LogP:   L = 0.5、T = 2.0、U = 4.0   （s=t=1、線形ランプ）
%[text] TPSA:   片側上限（下限は常に許容）:
%[text] d_TPSA = 1              TPSA <= 60 の場合
%[text] = (90 - TPSA)/30 60 < TPSA <= 90 の場合
%[text] = 0              TPSA > 90 の場合
%[text] MW:     L = 150、T = 350、U = 500   （s=1、t=1）
%[text]
logSection("A10", "セクション 3: 望ましさ関数 — プロパティからスコアへ", "アナリティクス L3");
%[text] --- 個別望ましさ関数（ベクトル化）---
d_logp = @(p) max(0, min(1, ...
    (p >= 0.5 & p <= 2.0) .* ((p - 0.5) ./ 1.5) + ...
    (p >  2.0 & p <= 4.0) .* ((4.0 - p) ./ 2.0) + ...
    (p == 2.0)));
%[text] 注意: 上の式は p=2 で重複する可能性があります。区分的に使用してください。
d_logp = @(p) ...
    ((p >= 0.5) & (p <  2.0)) .* ((p - 0.5) ./ 1.5) + ...
    ((p >= 2.0) & (p <= 4.0)) .* ((4.0 - p) ./ 2.0) + ...
    (p == 2.0);

d_tpsa = @(p) ...
    (p <= 60) .* 1.0 + ...
    (p >  60 & p <= 90) .* ((90 - p) ./ 30) + ...
    (p >  90) .* 0;

d_mw   = @(p) ...
    ((p >= 150) & (p <  350)) .* ((p - 150) ./ 200) + ...
    ((p >= 350) & (p <= 500)) .* ((500 - p) ./ 150) + ...
    (p == 350);

%[text] 複合望ましさ（3 成分の幾何平均）
D_composite = @(logp, tpsa, mw) ...
    (d_logp(logp) .* d_tpsa(tpsa) .* d_mw(mw)) .^ (1/3);

%[text] ライブラリのスコアを計算します。
d1 = d_logp(logp_vec);
d2 = d_tpsa(tpsa_vec);
d3 = d_mw(mw_vec);
D  = (d1 .* d2 .* d3) .^ (1/3);

%[text] リード化合物のスコアを計算します。
D_lead = D_composite(d_lead.LogP, d_lead.TPSA, d_lead.MolWt);
logInfo("リード（%s）望ましさ: D=%.3f (d_LogP=%.2f、d_TPSA=%.2f、d_MW=%.2f)", ...
    LEAD_NAME, D_lead, d_logp(d_lead.LogP), d_tpsa(d_lead.TPSA), d_mw(d_lead.MolWt));

%[text] ライブラリの統計情報を表示します。
logInfo("ライブラリ望ましさ: 平均=%.3f、最大=%.3f、中央値=%.3f", ...
    mean(D), max(D), median(D));

%[text] 望ましさ上位 10 ライブラリ分子を表示します。
[D_sorted, D_rank] = sort(D, "descend");
topN_des = min(10, nValid);
logInfo("複合望ましさ上位 %d ライブラリ分子:", topN_des);
for k = 1:topN_des
    idx = D_rank(k);
    logInfo("  %2d. D=%.3f  LogP=%5.2f  TPSA=%5.1f  MW=%5.1f  %s", ...
        k, D_sorted(k), logp_vec(idx), tpsa_vec(idx), mw_vec(idx), libNames(idx));
end

%[text] --- 望ましさ関数プロット ---
p_range_logp = linspace(-1, 6, 300);
p_range_tpsa = linspace(0, 150, 300);
p_range_mw   = linspace(100, 700, 300);

figure("Name", "A10 Sec3: 望ましさ関数");
subplot(1, 3, 1);
plot(p_range_logp, d_logp(p_range_logp), "b-", LineWidth=2);
xline(TARGET_LOGP, "k--", Label="Target", LabelVerticalAlignment="bottom");
xline(d_lead.LogP, "r:", LineWidth=1.5, Label=LEAD_NAME);
xlabel("LogP");  ylabel("d_{LogP}");
title("LogP 望ましさ");  ylim([0, 1.1]);  grid on;

subplot(1, 3, 2);
plot(p_range_tpsa, d_tpsa(p_range_tpsa), "g-", LineWidth=2);
xline(TARGET_TPSA, "k--", Label="Target", LabelVerticalAlignment="bottom");
xline(d_lead.TPSA, "r:", LineWidth=1.5, Label=LEAD_NAME);
xlabel("TPSA (A^2)");  ylabel("d_{TPSA}");
title("TPSA 望ましさ");  ylim([0, 1.1]);  grid on;

subplot(1, 3, 3);
plot(p_range_mw, d_mw(p_range_mw), "m-", LineWidth=2);
xline(TARGET_MW, "k--", Label="Target", LabelVerticalAlignment="bottom");
xline(d_lead.MolWt, "r:", LineWidth=1.5, Label=LEAD_NAME);
xlabel("MW (g/mol)");  ylabel("d_{MW}");
title("MW 望ましさ");  ylim([0, 1.1]);  grid on;

sgtitle("Derringer-Suich 望ましさ関数");

%[text] **💡 観察ポイント 3**
%[text] 幾何平均 D = (d1*d2*d3)^(1/3) は、単一の悪いスコアに対してペナルティを課します。
%[text] 算術平均に置き換えて結果を確認しましょう。
%[text] D_arith = (d1 + d2 + d3) / 3;
%[text] 上位 10 のランキングはどのように変わるかを確認し、どちらの定式化が低いプロパティを持つ化合物に対して保守的かを考察しましょう。
%[text] CNS 薬を好む LogP 望ましさに変更する場合（TPSA < 60、LogP 1〜3）:
%[text] d_tpsa のターゲットを 45、上限を 60 に調整します。
%[text] この条件でライブラリ分子がいくつ残るかを確認しましょう。
%[text] 最高の d_LogP スコアを持ち、最低の d_TPSA スコアを持つ分子を見つけましょう。
%[text] 低い TPSA スコアを説明する化学的特徴を考察しましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: 目標達成最適化（fgoalattain）
%[text]
%[text] 候補を望ましさスコアでランク付けしました。次に、fgoalattain を使用して3つのプロパティゴールを同時に満たす最適なブレンドを求めます。
%[text] 達成因子 gamma の役割を理解しましょう。
%[text]
%[text] ### コンセプト: 目標達成と達成因子
%[text]
%[text] 目標達成法（Charnes & Cooper 1977）は次の問題を解きます:
%[text]
%[text] min  gamma
%[text] s.t. f_j(x) - gamma * w_j <= goal_j    全 j = 1..m に対して
%[text] x in 実行可能集合 X
%[text]
%[text] ここで:
%[text] x = 決定変数
%[text] f_j(x) = j 番目の目的関数
%[text] goal_j = j 番目の目的の望ましいターゲット
%[text] w_j = 重み（正、スラックをスケーリング）
%[text] gamma = 達成因子（最適化変数）
%[text]
%[text] gamma の解釈:
%[text] gamma <= 0  =>  全てのゴールが同時に達成される
%[text] gamma >  0  =>  いくつかのゴールが過達成され、「最小」のゴールが gamma * w ユニット未達
%[text]
%[text] 問題定式化（この演習）:
%[text]
%[text] 「仮想化合物」をライブラリ分子の凸混合として表現します。
%[text] R^N（N = ライブラリサイズ）の決定変数 x は次を満たします:
%[text]
%[text] x >= 0、sum(x) = 1    （混合制約）
%[text]
%[text] ブレンドされたプロパティは次のようになります:
%[text]
%[text] f(x) = [sum_i x_i * LogP_i、sum_i x_i * TPSA_i、sum_i x_i * MW_i]
%[text] = PropMat' * x     （3 x 1 ベクトル）
%[text]
%[text] ここで PropMat は N x 3 のプロパティ行列です。
%[text]
%[text] fgoalattain は PropMat'*x* がミニマックスの意味でゴールベクトル [TARGET_LOGP、TARGET_TPSA、TARGET_MW] に最も近くなる x* を求めます。
%[text]
%[text] 最適ブレンド点に最も近い実分子:
%[text] 最適ブレンド点 p* = PropMat' * x* を計算した後、正規化プロパティ空間（ゼロ平均、単位分散）で p* に最小ユークリッド距離を持つライブラリ分子を見つけます。
%[text]
%[text] 注意: ライブラリの凸包がターゲット点を含まない可能性があります。
%[text] 達成因子 gamma はターゲットが凸包のどれだけ外側にあるかを定量化します。
logSection("A10", "セクション 4: 目標達成最適化（fgoalattain）", "アナリティクス L3");
n = nValid;
propGoal   = [TARGET_LOGP; TARGET_TPSA; TARGET_MW];
propWeight = [1.0; 1.0; 1.0];   % 全目的に等重み

%[text] 混合制約: sum(x) = 1、x >= 0
Aeq_mix = ones(1, n);
beq_mix = 1;
lb_mix  = zeros(n, 1);
ub_mix  = ones(n, 1);

%[text] 目的関数: ブレンドされたプロパティベクトル（x に線形）
fun_blend = @(x) propMat' * x(:);   % returns 3 x 1

if hasOptTbx
    % --- OPTIMIZATION TOOLBOX パス ---
    x0_mix   = ones(n, 1) / n;   % 一様混合から開始
    opts_ga  = optimoptions("fgoalattain", ...
        Display="off", MaxFunctionEvaluations=n*50, ...
        OptimalityTolerance=1e-6, ConstraintTolerance=1e-6);

    logInfo("fgoalattain を実行中（%d 変数、3 目的関数）...", n);
    tic;
    [x_opt, fval_blend, attainfactor, exitflag] = fgoalattain( ...
        fun_blend, x0_mix, propGoal, propWeight, ...
        [], [], Aeq_mix, beq_mix, lb_mix, ub_mix, [], opts_ga);
    t_opt = toc;

    logInfo("fgoalattain が %.2f 秒で完了（exitflag=%d）。", t_opt, exitflag);
    logInfo("達成因子 gamma = %.4f", attainfactor);
    if attainfactor <= 0
        logInfo("  => 全プロパティゴールが凸包内で達成。");
    else
        logInfo("  => ゴールが部分的に未達 -- ターゲットが凸包の外側にある。");
        logInfo("     最近接ブレンド点: LogP=%.2f、TPSA=%.1f、MW=%.1f", ...
            fval_blend(1), fval_blend(2), fval_blend(3));
    end

else
    % --- 手動フォールバック: 加重ユークリッド距離による最近接分子 ---
    logWarn("fgoalattain が利用できない -- 加重最近傍フォールバックを使用。");
    propSigma = std(propMat);           % normalisation scale
    propGoalN = propGoal(:)' ./ propSigma;
    propMatN  = propMat ./ propSigma;
    dist_all  = sqrt(sum((propMatN - propGoalN).^2, 2));
    [~, nn_idx] = min(dist_all);
    fval_blend = propMat(nn_idx, :)';  % 3 x 1
    attainfactor = max((fval_blend - propGoal) ./ propWeight);
    logInfo("最近傍: %s  (D_composite=%.3f)", libNames(nn_idx), D(nn_idx));
end

%[text] --- 最適ブレンド点に最も近い実分子を見つける ---
%[text] ライブラリ標準偏差で正規化（単位分散距離）
propSigma = std(propMat);
propMatN  = propMat ./ propSigma;       % N x 3 normalised
fvalN     = fval_blend(:)' ./ propSigma;  % 1 x 3

dist_to_blend = sqrt(sum((propMatN - fvalN).^2, 2));
[~, nearestRank] = sort(dist_to_blend);

%[text] 上位 5 の最近接実分子を表示
TOP_GA = 5;
logInfo("最適ブレンド点に最近接の上位 %d 実分子:", TOP_GA);
logInfo("  （最適ブレンド: LogP=%.2f、TPSA=%.1f、MW=%.1f）", ...
    fval_blend(1), fval_blend(2), fval_blend(3));
for k = 1:min(TOP_GA, nValid)
    idx = nearestRank(k);
    logInfo("  %d. 距離=%.3f  D=%.3f  LogP=%5.2f  TPSA=%5.1f  MW=%5.1f  %s", ...
        k, dist_to_blend(idx), D(idx), logp_vec(idx), tpsa_vec(idx), ...
        mw_vec(idx), libNames(idx));
end

%[text] --- レーダー（スパイダー）チャート: リード vs 上位望ましさ候補 ---
%[text] 注意: gamma << 0 のとき、fgoalattain ブレンド最近傍（nearestRank(1)）は LogP~-3.5 に到達します（fgoalattain がゴール未達の最大ギャップを最小化するため）。
%[text] 意味のある比較のために上位望ましさ分子を代わりに使用します。
[~, topGaIdx]  = max(D);
radarProps = ["LogP", "TPSA", "MW"];
radarLead  = [d_lead.LogP, d_lead.TPSA, d_lead.MolWt];
radarCand  = [logp_vec(topGaIdx), tpsa_vec(topGaIdx), mw_vec(topGaIdx)];
radarGoal  = [TARGET_LOGP, TARGET_TPSA, TARGET_MW];

%[text] 固定の薬らしい空間境界を使って [0,1] に正規化します。
%[text] ライブラリの min/max を使うと極端な外れ値（MW=924、TPSA=319）が範囲を伸ばすため、チャートが圧縮されます。固定境界は 3 化合物を比較可能で人間が読みやすい領域に保ちます。
radarMin   = [-1,   0, 100];   % [LogP_lo, TPSA_lo, MW_lo]
radarMax   = [ 5, 120, 500];   % [LogP_hi, TPSA_hi, MW_hi]
radarNorm  = @(v) (v - radarMin) ./ (radarMax - radarMin);
radarNorm  = @(v) max(0, min(1, (v - radarMin) ./ (radarMax - radarMin)));  % clip to [0,1]

figure("Name", "A10 Sec4: プロパティプロファイル -- リード vs 候補 vs ゴール");
theta    = linspace(0, 2*pi, numel(radarProps) + 1);
nLead    = radarNorm(radarLead);   % 1x3
nCand    = radarNorm(radarCand);   % 1x3
nGoal    = radarNorm(radarGoal);   % 1x3
rLead    = [nLead, nLead(1)];      % 1x4  (close the polygon)
rCand    = [nCand, nCand(1)];      % 1x4
rGoal    = [nGoal, nGoal(1)];      % 1x4

polarplot(theta, rLead, "r-o", LineWidth=2, DisplayName=LEAD_NAME);
hold on;
polarplot(theta, rCand, "b-s", LineWidth=2, DisplayName=libNames(topGaIdx));
polarplot(theta, rGoal, "k--^", LineWidth=1.5, DisplayName="Goal");
legend(Location="southoutside");
title(sprintf("プロパティプロファイル（正規化）\n%s vs %s vs ゴール", ...
    LEAD_NAME, libNames(topGaIdx)));
pax = gca;
pax.ThetaTick      = [0, 120, 240];
pax.ThetaTickLabel = {"LogP", "TPSA (A^2)", "MW (g/mol)"};

%[text] **💡 観察ポイント 4**
%[text] 達成因子 gamma は何かを確認しましょう。
%[text] gamma <= 0 なら全ゴールが達成可能で、gamma > 0 なら一部が未達です。
%[text] TARGET_MW を 350 から 250 に変えて再実行し、gamma の変化を確認しましょう。
%[text] gamma はどう変わるか、どのゴールが最も達成しにくくなるかを読み取りましょう。
%[text] 注意: gamma << 0（例: -5.5）のとき、fgoalattain は全ゴールを同時にできるだけ下回るようにブレンドを駆動します。つまり「最適ブレンド」点がゴールから遠く離れる可能性があります（例: LogP = 2.0 の代わりに -3.5）。そのブレンド点への上位 5 の最近接分子は D = 0 になることが多い（望ましさゾーンの外）。
%[text] 最近接の実行可能化合物を見つけるには fgoalattain ブレンド最近傍ではなく上位 D ランキング（BETAXOLOL）を使用します。
%[text] 重みを変えて TPSA より LogP を優先させた場合を確認しましょう:
%[text] propWeight = [0.5; 2.0; 1.0];  % TPSA を 2 倍重み
%[text] fgoalattain に「TPSA ゴールを逃すコストは 2 倍」と伝えます。
%[text] どの分子が上位 1 位に移動するかを確認しましょう。
%[text] 「混合」の解釈は物理的に何を意味するかを考えましょう。
%[text] x_opt は各ライブラリ分子に重みを割り当てます。x_opt が重み [0.3、0.4、0.3] の 3 つの非ゼロエントリを持つなら、「理想的なアナログ」が 3 分子全ての構造特徴を組み合わせることを意味します。
%[text] 実際のリード最適化では、どの置換基をマージするかを導きます。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: パレートフロント -- LogP vs TPSA トレードオフ
%[text]
%[text] 目標達成最適化の結果が得られました。次に LogP（分配係数の常用対数）と TPSA（位相的極性表面積）のパレートフロントを描画し、トレードオフの全体像を把握しましょう。
%[text] どの化合物がパレート最適で、どの化合物が劣位かを判別しましょう。
%[text]
%[text] ### コンセプト: パレート最適性と epsilon 制約法
%[text]
%[text] 解 x がパレート最適（非劣位）とは、全ての目的において同時により良い他の実行可能解が存在しない場合を指します。
%[text]
%[text] 2 目的コンテキスト（|LogP - 2.0| と TPSA を最小化）において:
%[text] 分子 A が分子 B を支配する条件:
%[text] A が |LogP - target| と TPSA の両方で B より低い。
%[text]
%[text] パレートフロントは全非劣位解の集合です。
%[text] これは2つの目的間の真のトレードオフを表します。
%[text]
%[text] EPSILON 制約法（Cohon & Marks 1975）:
%[text] 目的を1つ（加重和）に結合する代わりに、もう一方を制約しながら一方の目的を最適化します:
%[text]
%[text] min  TPSA(x) = tpsa' * x
%[text] s.t. LogP(x) = logp' * x >= epsilon      （LogP 下限）
%[text] sum(x) = 1、x >= 0
%[text]
%[text] epsilon をライブラリの LogP_min から LogP_max まで掃引することで完全なパレートフロントを追跡します。各 epsilon 値での線形計画（LP）は linprog（Optimization Toolbox）で解きます。
%[text]
%[text] フォールバック（Optimization Toolbox なし）:
%[text] 非劣位性を確認してライブラリから直接経験的パレートフロントを抽出します:
%[text] 分子 i は、TPSA_j < TPSA_i かつ |LogP_j - 2| < |LogP_i - 2| を満たす分子 j が存在しない場合に非劣位とされます。
%[text]
%[text] 混合変数 x に関する目的関数の定義
%[text] obj1 = |LogP(x) - TARGET_LOGP|  （一方向使用: LogP(x) >= eps）
%[text] obj2 = TPSA(x)                  （TPSA を最小化）
%[text]
%[text] linprog 混合アプローチに関する注意:
%[text] 200 分子混合に対する linprog の epsilon 制約法は数学的には有効ですが、実際には誤解を招くパレートフロントを与えることがあります。
%[text] 凸結合は常に最低 TPSA 分子（TPSA=3.2）に ~100% の重みを割り当てるため、最適 TPSA はどの LogP でも自明にほぼゼロになります。フロントはプロット下端のフラットラインに収縮し、実際の薬から遠ざかります。
%[text] 代わりに経験的パレートフロント（非劣位実分子）を使用します。
%[text]
%[text] --- 経験的パレートフロント: ライブラリからの非劣位セット ---
logSection("A10", "セクション 5: パレートフロント — LogP vs TPSA トレードオフ", "アナリティクス L3");
obj1 = abs(logp_vec - TARGET_LOGP);   % minimise deviation from target LogP
obj2 = tpsa_vec;                        % minimise TPSA

nondominatedMask = true(nValid, 1);
for ii = 1:nValid
    for jj = 1:nValid
        if jj ~= ii && obj1(jj) <= obj1(ii) && obj2(jj) <= obj2(ii) && ...
           (obj1(jj) < obj1(ii) || obj2(jj) < obj2(ii))
            nondominatedMask(ii) = false;
            break;
        end
    end
end
pIdx        = find(nondominatedMask);
[~, sortP]  = sort(logp_vec(pIdx));
pIdx        = pIdx(sortP);
pareto_logp = logp_vec(pIdx);
pareto_tpsa = tpsa_vec(pIdx);
logInfo("経験的パレートフロント: %d 非劣位分子。", numel(pIdx));

if hasOptTbx
    % linprog concept demo (NOT used for the plot -- see NOTE above).
    % Uncomment to observe the mixture-front collapse:
    %   N_LP = 20;
    %   eps_sw = linspace(min(logp_vec)+0.1, max(logp_vec)-0.1, N_LP);
    %   opt_lp = optimoptions("linprog", Display="off");
    %   for k = 1:N_LP
    %       [~, tv, fk] = linprog(tpsa_vec, -logp_vec', -eps_sw(k), ...
    %           ones(1,n), 1, zeros(n,1), ones(n,1), opt_lp);
    %       if fk==1, fprintf("eps=%.2f TPSA_mix=%.1f\n",eps_sw(k),tv); end
    %   end
end

%[text] --- パレートフロントプロット ---
figure("Name", "A10 Sec5: LogP vs TPSA パレートフロント");

%[text] ライブラリ背景
scatter(logp_vec, tpsa_vec, 25, D, "filled", ...
    MarkerFaceAlpha=0.5, DisplayName="Library (colour = D)");
colormap("cool");
cb2 = colorbar();
cb2.Label.String = "複合望ましさ D";
hold on;

%[text] パレートフロント曲線
plot(pareto_logp, pareto_tpsa, "k-o", LineWidth=2, ...
    MarkerFaceColor="k", MarkerSize=5, DisplayName="Pareto front");

%[text] 理想（ターゲット）点
plot(TARGET_LOGP, TARGET_TPSA, "rp", MarkerSize=18, ...
    LineWidth=2, DisplayName="Ideal target");

%[text] リード化合物
plot(d_lead.LogP, d_lead.TPSA, "rv", MarkerSize=12, ...
    MarkerFaceColor="r", DisplayName=LEAD_NAME);

%[text] Veber TPSA 閾値
yline(90, "g--", LineWidth=1.5, Label="Veber TPSA=90", DisplayName="Veber TPSA=90");
yline(60, "b--", LineWidth=1.5, Label="CNS TPSA=60",  DisplayName="CNS TPSA=60");

xlabel("LogP");
ylabel("TPSA (A^2)");
title("LogP vs TPSA: パレートフロント（epsilon 制約）");
legend(Location="northeast");
grid on;

%[text] **💡 観察ポイント 5**
%[text] パレートフロント曲線を観察し、左から右へ（LogP が増加）動くときに TPSA がどのように変化するかを確認しましょう。トレードオフは単調でしょうか？
%[text] これはライブラリの多様性について何を示しているでしょうか？
%[text] 理想点 [LogP=2、TPSA=70] は星で示されています。
%[text] この点がパレートフロント上にあるかどうか、もしない場合は上か下かを確認しましょう。
%[text] （フロントより下の点は実行不可能です。ライブラリ分子や混合がその組み合わせを達成できないことを示します。）
%[text] パレート掃引を MW（分子量）を2番目の目的として使うように変更してみましょう。
%[text] tpsa_vec を mw_vec に、TARGET_TPSA を TARGET_MW に置換します。
%[text] linprog を再実行し、LogP-MW トレードオフが LogP-TPSA トレードオフより強いか弱いかを確認しましょう。（ヒント: まずピアソン相関を確認してください。）
%[text] epsilon 制約法は LogP 空間でパレートフロントを一様にサンプリングします。
%[text] フロントの「平坦な」領域の点は、より高い LogP を受け入れても利益が少ないことを示しています。
%[text] 「膝点」を特定できるか確認しましょう -- LogP の増加が最大の TPSA 削減をもたらすパレート点です。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: リード最適化サマリーと候補可視化
%[text]
%[text] パレートフロントが得られました。最終的に望ましさ、目標距離、パレート最適性を組み合わせて上位候補を選びます。
%[text] 最終候補のテーブルと構造式を表示し、演習を完了します。
%[text]
%[text] ### コンセプト: 候補の選択 — 望ましさと近接性の組み合わせ
%[text]
%[text] 最終選択は次の3つの基準を組み合わせます:
%[text] (a)  複合望ましさ D（セクション 3 参照）
%[text] (b)  目標達成最適値への近接性（セクション 4 参照）
%[text] (c)  LogP vs TPSA のパレート非劣位性（セクション 5 参照）
%[text]
%[text] 実際のプロジェクトでは、以下の追加フィルターが適用されます:
%[text] - 合成アクセシビリティスコア（SAS）
%[text] - hERG 心毒性予測
%[text] - PAINS / 構造アラートフィルター（S03 参照）
%[text] - ADMET モデリング（吸収、分布、代謝、排泄、毒性）
%[text] 例: pkCSM または SwissADME ウェブサーバー
%[text]
%[text] --- 最終候補テーブルを構築 ---
logSection("A10", "セクション 6: リード最適化サマリーと候補可視化", "アナリティクス L3");
TOP_FINAL = 10;

%[text] 複合望ましさに基づいてランク付け
[D_sorted_f, rank_by_D] = sort(D, "descend");

candidateTbl = table( ...
    (1:TOP_FINAL)', ...
    libNames(rank_by_D(1:TOP_FINAL)), ...
    round(logp_vec(rank_by_D(1:TOP_FINAL)), 2), ...
    round(tpsa_vec(rank_by_D(1:TOP_FINAL)), 1), ...
    round(mw_vec(rank_by_D(1:TOP_FINAL)), 1), ...
    round(D_sorted_f(1:TOP_FINAL), 3), ...
    round(dist_to_blend(rank_by_D(1:TOP_FINAL)), 3), ...
    VariableNames=["Rank", "Name", "LogP", "TPSA_A2", "MW_gmol", ...
                   "Desirability", "Dist2Goal"]);

logInfo("--- 複合望ましさ上位 %d 候補 ---", TOP_FINAL);
disp(candidateTbl);

%[text] --- リード vs 上位候補の比較 ---
bestIdx = rank_by_D(1);

logInfo("--- リード vs 最良候補 ---");
logInfo("               %-25s  %-25s  ターゲット", LEAD_NAME, libNames(bestIdx));
logInfo("  LogP       : %8.2f                      %8.2f          %.1f", ...
    d_lead.LogP, logp_vec(bestIdx), TARGET_LOGP);
logInfo("  TPSA (A^2) : %8.1f                      %8.1f          %.1f", ...
    d_lead.TPSA, tpsa_vec(bestIdx), TARGET_TPSA);
logInfo("  MW (g/mol) : %8.1f                      %8.1f          %.1f", ...
    d_lead.MolWt, mw_vec(bestIdx), TARGET_MW);
logInfo("  望ましさ: %7.3f                      %7.3f          1.000", ...
    D_lead, D(bestIdx));

%[text] --- 上位 5 候補を可視化 ---
TOP_VIZ = 5;
vizMols  = cell(1, TOP_VIZ);
vizNames = strings(1, TOP_VIZ);
for k = 1:TOP_VIZ
    idx          = rank_by_D(k);
    vizMols{k}   = mols{idx};
    vizNames(k)  = sprintf("%s\nD=%.3f LogP=%.2f MW=%.0f", ...
        libNames(idx), D(idx), logp_vec(idx), mw_vec(idx));
end

logInfo("上位 %d 候補を描画中...", TOP_VIZ);
for k = 1:TOP_VIZ
    try
        figure();  % create a fresh figure so draw2d does not reuse the Pareto axes
        fig = emk.viz.draw2d(vizMols{k}, Title=vizNames(k));
        set(fig, "Name", sprintf("A10 Sec6: Candidate %d -- %s", ...
            k, libNames(rank_by_D(k))));
    catch ME
        logWarn("候補 %d の draw2d が失敗: %s", k, ME.message);
    end
end

%[text] --- サマリー棒グラフ: 望ましさ内訳 ---
barData = [d_logp(d_lead.LogP),  d_tpsa(d_lead.TPSA),  d_mw(d_lead.MolWt);   ...
           d_logp(logp_vec(rank_by_D(1:TOP_VIZ-1))), ...
           d_tpsa(tpsa_vec(rank_by_D(1:TOP_VIZ-1))), ...
           d_mw(mw_vec(rank_by_D(1:TOP_VIZ-1)))];

barLabels = [LEAD_NAME; libNames(rank_by_D(1:TOP_VIZ-1))];
%[text] 表示用に長い名前を短縮
for k = 1:numel(barLabels)
    if strlength(barLabels(k)) > 18
        barLabels(k) = extractBefore(barLabels(k), 18) + "...";
    end
end

figure("Name", "A10 Sec6: 望ましさ内訳 -- リード vs 上位候補");
b = bar(barData, "grouped");
b(1).FaceColor = [0.2 0.6 0.9];   % d_LogP
b(2).FaceColor = [0.3 0.8 0.4];   % d_TPSA
b(3).FaceColor = [0.9 0.5 0.2];   % d_MW
legend(["d_{LogP}", "d_{TPSA}", "d_{MW}"], Location="northeast");
xticklabels(barLabels);
xtickangle(20);
ylabel("個別望ましさスコア");
title("望ましさ内訳: リード vs 上位候補");
ylim([0, 1.05]);
yline(1.0, "k--", HandleVisibility="off");
grid on;

%[text] **💡 観察ポイント 6**
%[text] 最高の d_TPSA スコアと最高の d_LogP スコアを持つ分子をそれぞれ確認しましょう。同じ分子ですか?
%[text] これは多目的選択の課題について何を示していますか?
%[text] S03 の PAINS フィルターを候補リストに追加してみましょう:
%[text] pains_smarts = emk.filter.loadPainsSmarts();  % （利用可能な場合）
%[text] hasPains = cellfun(@(m) emk.mol.hasSubstruct(m, pains_smarts), mols);
%[text] 上位 10 リストから PAINS 陽性候補を削除します。
%[text] いくつ残りますか?
%[text] 実際のリード最適化では合成アクセシビリティも考慮します。
%[text] 合成の複雑さに基づいて上位候補を手動でランク付けしてみましょう。
%[text] 高い望ましさは容易な合成と相関していますか?
%[text] 自分のリード化合物とプロパティターゲットのセットを選びましょう。
%[text] セクション 1 からフルパイプラインを再実行します。
%[text] 推奨リード:
%[text] カフェイン: "Cn1cnc2c1c(=O)n(C)c(=O)n2C"   — CNS 刺激薬
%[text] イブプロフェン: "CC(C)Cc1ccc(cc1)C(C)C(=O)O"   — NSAID
%[text] メトホルミン: "CN(C)C(=N)NC(=N)N"              — 糖尿病薬
%[text]
%[text] **まとめ**
%[text] Derringer-Suich 望ましさ関数で薬理学的直感を無次元スコアに変換し、fgoalattain で全プロパティゴールを同時に最小化する目標達成最適化を適用しました。
%[text] epsilon 制約法で LogP vs TPSA のパレートフロントを追跡し、劣位候補を除外しました。
logInfo("A10: 完了。");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
