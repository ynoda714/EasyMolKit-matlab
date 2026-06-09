%[text] # A10 解答: リード最適化 -- 多目的プロパティ最適化
%[text] a10_lead_optimization.m の「やってみよう」演習の参照解答。
%[text] 最初に a10_lead_optimization.m（最低セクション 0〜3）を実行して
%[text] 必要なワークスペース変数を構築してください。その後このファイルで確認。
%[text] addpath(genpath("src"));
%[text] emk.setup.initPython();

%[text] mol_warmup = emk.mol.fromSmiles("C"); clear mol_warmup;
%[text] hasOptTbx  = license("test", "optimization_toolbox");
%[text] logInfo("A10 解答: セットアップ完了  (OptTbx=%d)", hasOptTbx);
%%
%[text] ## 前提条件の再構築（a10 セクション 0〜3 の再現）

LEAD_SMILES = "CC(=O)Oc1ccccc1C(=O)O";
LEAD_NAME   = "Aspirin";
TARGET_LOGP = 2.0;
TARGET_TPSA = 70.0;
TARGET_MW   = 350.0;

mol_lead = emk.mol.fromSmiles(LEAD_SMILES);
d_lead   = emk.descriptor.calculate(mol_lead);

DATA_FILE = "data/list/fda_drugs.csv";
rawTbl    = readtable(DATA_FILE, TextType="string");
nLib      = height(rawTbl);
mols      = cell(1, nLib);
valid     = false(1, nLib);
for k = 1:nLib
    try; mols{k} = emk.mol.fromSmiles(rawTbl.SMILES(k)); valid(k) = true; catch; end
end
validIdx  = find(valid);
mols      = mols(validIdx);
libNames  = rawTbl.Name(validIdx);
nValid    = numel(mols);

DESCS   = ["LogP", "TPSA", "MolWt"];
descTbl = emk.descriptor.batchCalculate(mols, DESCS);
logp_vec = descTbl.LogP;
tpsa_vec = descTbl.TPSA;
mw_vec   = descTbl.MolWt;
propMat  = [logp_vec, tpsa_vec, mw_vec];

d_logp = @(p) ...
    ((p >= 0.5) & (p <  2.0)) .* ((p - 0.5) ./ 1.5) + ...
    ((p >= 2.0) & (p <= 4.0)) .* ((4.0 - p) ./ 2.0) + ...
    (p == 2.0);
d_tpsa = @(p) (p <= 60) .* 1.0 + (p > 60 & p <= 90) .* ((90 - p) ./ 30);
d_mw   = @(p) ...
    ((p >= 150) & (p <  350)) .* ((p - 150) ./ 200) + ...
    ((p >= 350) & (p <= 500)) .* ((500 - p) ./ 150) + ...
    (p == 350);

d1 = d_logp(logp_vec);
d2 = d_tpsa(tpsa_vec);
d3 = d_mw(mw_vec);
D  = (d1 .* d2 .* d3) .^ (1/3);

logInfo("A10 前提条件再構築完了（%d 有効分子）。", nValid);
%%
%[text] ## やってみよう 1: アスピリン MW の「余裕」、イブプロフェン比較、CNS vs 末梢ターゲット
%[text] Q: アスピリン MW = 180 g/mol は 350 g/mol ターゲットを大きく下回る。
%[text] 解答: MW が低いと創薬化学者に「分子量バジェット」（~170 Da）が生まれ、
%[text]    Lipinski Ro5（MW <= 500）を逸脱せずに効力・選択性・PK を改善する置換基を
%[text]    追加できる。現代の創薬ではリード化合物を小さくスタートすることが
%[text]    強く推奨される。
mol_ibu = emk.mol.fromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
d_ibu   = emk.descriptor.calculate(mol_ibu);

logInfo("やってみよう 1 -- アスピリン vs イブプロフェン vs ターゲット:");
logInfo("  %12s: LogP=%+5.2f TPSA=%5.1f MW=%5.1f", ...
    LEAD_NAME, d_lead.LogP, d_lead.TPSA, d_lead.MolWt);
logInfo("  %12s: LogP=%+5.2f TPSA=%5.1f MW=%5.1f", ...
    "Ibuprofen", d_ibu.LogP, d_ibu.TPSA, d_ibu.MolWt);
logInfo("  %12s: LogP=%+5.2f TPSA=%5.1f MW=%5.1f", ...
    "Target", TARGET_LOGP, TARGET_TPSA, TARGET_MW);

%[text] 解答: イブプロフェン（MW=206、LogP~3.1、TPSA~37）はアスピリン（MW=180、LogP=1.3、TPSA=63）
%[text]    より LogP ターゲットに近い。TPSA は 70 A^2 ターゲットを大幅に下回る。
%[text]    ただし、どちらも MW=350 には届かない ── 両者とも「リーン」な NSAIDs。
%[text]
%[text] 末梢性抗炎症薬の TARGET_TPSA: < 90 A^2（Veber ルール）。
%[text] CNS 薬: < 60 A^2（血液脳関門基準）。
%[text] アスピリンとイブプロフェンはいずれも末梢基準を満たす。CNS 基準では
%[text] アスピリンは境界線上（TPSA=63、60 A^2 限界に近い）。
%%
%[text] ## やってみよう 2: LogP-TPSA 相関; 3 基準全て満たす分子数; Ro5 ボックス

r_logp_tpsa = corr(logp_vec, tpsa_vec);
logInfo("やってみよう 2 -- ピアソン r(LogP, TPSA) = %+.3f", r_logp_tpsa);

%[text] 解答: r ~ -0.3〜-0.5（弱い負の相関）。親油性化合物（高 LogP）は
%[text]    極性基が少ない傾向（低 TPSA）。ただし相関は弱く、例外も多い
%[text]    （例: 塩素化芳香族は高 LogP だが、塩素置換により TPSA は中程度）。
n_all_three = sum((logp_vec >= 1.0 & logp_vec <= 3.5) & ...
                   tpsa_vec <= 90 & mw_vec <= 400);
logInfo("やってみよう 2 -- 3 つの厳格基準全てを満たす分子: %d / %d", ...
    n_all_three, nValid);

%[text] 解答: 3 基準全てを同時に満たすのはこの FDA ライブラリで 75 / 200（37.5%）。
%[text]    「経口薬スイートスポット」は承認薬に色濃く反映されている ──
%[text]    ADMET で失敗した分子は FDA 承認にほとんど届かないため。
%[text]
%[text] Lipinski vs Veber 比較:
n_lipinski = sum(logp_vec <= 5 & mw_vec <= 500);
n_veber    = sum(tpsa_vec <= 90 & mw_vec <= 400);
logInfo("やってみよう 2 -- Lipinski (MW<=500, LogP<=5): %d 分子", n_lipinski);
logInfo("やってみよう 2 -- Veber (TPSA<=90, MW<=400):   %d 分子", n_veber);

%[text] 解答: ほとんどの FDA 薬は Lipinski を通過（~90%）。より厳しい
%[text]    Veber / 400 Da ガイドラインを通過するのは少なめ（~60-70%）。
%[text]    FDA ライブラリは承認薬に偏っており、定義上ヒトで許容可能な吸収を持つ。
%%
%[text] ## やってみよう 3: 算術平均 vs 幾何平均; CNS TPSA ターゲット; 最大 d_LogP & 最小 d_TPSA
%[text] デフォルト（幾何）と算術平均の望ましさ
D_arith = (d1 + d2 + d3) / 3;
[D_geo_top10,  idx_geo]   = maxk(D,       10);
[D_arith_top10,idx_arith] = maxk(D_arith, 10);

overlap_top10 = numel(intersect(idx_geo, idx_arith));
logInfo("やってみよう 3 -- 上位 10 重複（幾何 vs 算術平均）: %d / 10", overlap_top10);

%[text] 解答: 算術平均は d_i = 0 が 1 個でも他の 2 つが高ければ許容する。
%[text]    幾何平均はいずれか 1 つの d_i = 0 で全体が 0 に崩壊する。
%[text]    したがって算術平均は「1 プロパティが不合格でも他の 2 つが優秀な
%[text]    分子」を昇格させる ── 創薬では通常望ましくない挙動。
%[text]    この FDA ライブラリでは上位 10 の重複は 9/10 ──
%[text]    上位候補は 3 プロパティ全てが範囲内（d_i = 0 なし）なので
%[text]    ランキングはほぼ変わらない。
%[text]    「1 プロパティ外れ値」（超高 MW や極端な LogP）が多いライブラリでは
%[text]    重複が顕著に低下（5-7）する。
%[text]
%[text] CNS TPSA ターゲット: [0, 60] に絞り込む
d_tpsa_cns = @(p) ...
    (p <= 45) .* 1.0 + ...
    (p > 45 & p <= 60) .* ((60 - p) ./ 15) + ...
    (p > 60) .* 0;

D_cns = (d_logp(logp_vec) .* d_tpsa_cns(tpsa_vec) .* d_mw(mw_vec)) .^ (1/3);
n_cns_survive = sum(D_cns > 0);
logInfo("やってみよう 3 -- D_CNS > 0 の分子（TPSA < 60）: %d / %d", ...
    n_cns_survive, nValid);

%[text] 解答: このライブラリで CNS 基準（TPSA < 60）を通過するのは 50 / 200（25%）。
%[text]    FDA 薬の多くは末梢組織をターゲットとするため、
%[text]    CNS 透過性に最適化されていない。
%[text]
%[text] 最高 d_LogP だが最低 d_TPSA の分子
[~, i_top_logp] = max(d1);
[~, i_low_tpsa] = min(d2);
logInfo("やってみよう 3 -- 最高 d_LogP: %s (LogP=%.2f, TPSA=%.1f)", ...
    libNames(i_top_logp), logp_vec(i_top_logp), tpsa_vec(i_top_logp));
logInfo("やってみよう 3 -- 最低 d_TPSA:  %s (LogP=%.2f, TPSA=%.1f)", ...
    libNames(i_low_tpsa), logp_vec(i_low_tpsa), tpsa_vec(i_low_tpsa));

%[text] 解答: 最低 d_TPSA 分子は PRAZOSIN（LogP=1.78、TPSA=107、d_TPSA=0）。
%[text]    LogP は中程度だが、ピペラジン環・フラン環・2 つのアミド基を持ち、
%[text]    各々が TPSA に寄与する。TPSA の高さは親油性ではなく
%[text]    極性基の数によって決まることを示す。
%[text]    最高 d_LogP 分子は METOCLOPRAMIDE（LogP=2.00、d_LogP=1.0）。
%[text]    LogP がターゲットのピーク値に完全一致する。両者は異なる分子であり、
%[text]    d_LogP と d_TPSA が薬らしさの独立した軸であることを確認する。
%%
%[text] ## やってみよう 4: fgoalattain / gamma の解釈; MW ターゲット変更; 重み変更

if ~hasOptTbx
    logWarn("やってみよう 4 -- Optimization Toolbox が利用不可。概念のみ表示。");
    logInfo("やってみよう 4 -- gamma <= 0 は全ゴールが凸包内にあることを意味する");
    logInfo("           gamma > 0 はターゲットが外にあることを意味する。");
    logInfo("           確認方法: Optimization Toolbox で a10 セクション 4 を実行。");
else
    % --- MW ターゲット厳格化の影響 ---
    propGoal_tight = [TARGET_LOGP; TARGET_TPSA; 250.0];  % MW = 250 instead of 350
    propWeight_def = [1.0; 1.0; 1.0];

    fun_blend = @(x) propMat' * x(:);
    n = nValid;
    x0  = ones(n, 1) / n;
    opt = optimoptions("fgoalattain", Display="off");

    [~, ~, gamma_350, ~] = fgoalattain(fun_blend, x0, ...
        [TARGET_LOGP; TARGET_TPSA; 350], propWeight_def, ...
        [], [], ones(1,n), 1, zeros(n,1), ones(n,1), [], opt);
    [~, ~, gamma_250, ~] = fgoalattain(fun_blend, x0, ...
        propGoal_tight, propWeight_def, ...
        [], [], ones(1,n), 1, zeros(n,1), ones(n,1), [], opt);

    logInfo("やってみよう 4 -- gamma (MW ターゲット=350): %.4f", gamma_350);
    logInfo("やってみよう 4 -- gamma (MW ターゲット=250): %.4f", gamma_250);

    % 解答: TARGET_MW を 350 から 250 に狭めても gamma は変わらない。
    %    ライブラリの凸包が MW = [46, 924] g/mol の範囲にまたがるため。
    %    凸混合はその範囲内の任意の MW を自明に達成できるため、
    %    両ターゲットとも凸包の深い内側にある（gamma << 0）。
    %    gamma がゼロに近づく様子を観察するには、ライブラリの最小値を
    %    下回る TARGET_MW（例: 40 g/mol）または最大値を超える値（例: 1000 g/mol）に設定。

    % --- TPSA を 2 倍重み付けした場合の影響 ---
    propWeight_tpsa = [0.5; 2.0; 1.0];
    [~, fval_wtpsa, gamma_wtpsa, ~] = fgoalattain(fun_blend, x0, ...
        [TARGET_LOGP; TARGET_TPSA; TARGET_MW], propWeight_tpsa, ...
        [], [], ones(1,n), 1, zeros(n,1), ones(n,1), [], opt);

    propSigma = std(propMat);
    propMatN  = propMat ./ propSigma;
    fvalN     = fval_wtpsa(:)' ./ propSigma;
    dist_tw   = sqrt(sum((propMatN - fvalN).^2, 2));
    [~, top_tw] = min(dist_tw);
    logInfo("やってみよう 4 -- 上位 1 位（TPSA 重み 2 倍）: %s", libNames(top_tw));

    % 解答: With double TPSA weight, fgoalattain penalises deviations from
    %    TARGET_TPSA more strongly.  The optimal blend shifts toward
    %    molecules with TPSA close to 70.  The top-1 real molecule may
    %    change if there is a drug with good TPSA but only average LogP/MW.
end
%%
%[text] ## やってみよう 5: パレートフロントの解釈（概念的質問 -- コード不要）

logInfo("やってみよう 5 -- パレートフロント概念解答:");
logInfo("  Q1: TPSA は LogP の増加とともに減少（負のトレードオフ）。");
logInfo("      整った分布では単調; 多様なスキャフォールドからは例外あり。");
logInfo("      フラット領域 = 同じ LogP 帯で類似 TPSA を持つ多くの薬。");
logInfo("      ");
logInfo("  Q2: 理想 [LogP=2, TPSA=70] はおそらくパレートフロント上（達成困難）。");
logInfo("      ライブラリは通常、その複合プロパティを持つ薬を含まない。");
logInfo("      ");
logInfo("  Q3: LogP-MW 相関はこのライブラリで +0.337（弱い正の相関）。");
logInfo("      LogP-TPSA (-0.333) と同程度の大きさだが逆符号。");
logInfo("      大きく重いスキャフォールドはより親油性が高い傾向（アルキル鎖/環付加）。");
logInfo("      極性基は LogP を下げるが MW への影響は TPSA 上昇より小さい。");
logInfo("      ");

%[text] 確認用のピアソン相関
r_logp_mw = corr(logp_vec, mw_vec);
logInfo("  ピアソン r(LogP,TPSA)=%.3f  r(LogP,MW)=%.3f", ...
    corr(logp_vec,tpsa_vec), r_logp_mw);

%[text] Q4: 膝点 -- |Δ TPSA| / |Δ LogP| の最大
%[text] セクション 5 の pareto_logp / pareto_tpsa ベクトルが必要。
%[text] a10 セクション 5 を先に実行してから以下を実行:
%[text]   dTPSA = diff(pareto_tpsa);
%[text]   dLogP = diff(pareto_logp);
%[text]   [~, knee] = max(abs(dTPSA ./ dLogP));
%[text]   logInfo("Knee at LogP = %.2f, TPSA = %.1f", pareto_logp(knee), pareto_tpsa(knee));
%%
%[text] ## やってみよう 6: d_TPSA vs d_LogP チャンピオン; PAINS フィルター; 合成複雑性
%[text] 最高 d_TPSA vs 最高 d_LogP の分子
[~, i_best_tpsa] = max(d2);
[~, i_best_logp] = max(d1);
same_molecule = isequal(i_best_tpsa, i_best_logp);

logInfo("やってみよう 6 -- 最高 d_TPSA: %s (TPSA=%.1f, LogP=%.2f, D=%.3f)", ...
    libNames(i_best_tpsa), tpsa_vec(i_best_tpsa), logp_vec(i_best_tpsa), D(i_best_tpsa));
logInfo("やってみよう 6 -- 最高 d_LogP: %s (LogP=%.2f, TPSA=%.1f, D=%.3f)", ...
    libNames(i_best_logp), logp_vec(i_best_logp), tpsa_vec(i_best_logp), D(i_best_logp));
logInfo("やってみよう 6 -- 同じ分子か? %d  (多目的テンションを示す)", same_molecule);

%[text] 解答: d_TPSA = 1.0 は TPSA <= 60 の全分子に当てはまるため、max(d2) は
%[text]    配列順で最初の該当分子（NICOTINE、TPSA=16.1、LogP=1.85）を返す。
%[text]    「最高 d_TPSA」は一意でなく ── 50 分子が d_TPSA=1 で同率首位。
%[text]    LogP チャンピオンは METOCLOPRAMIDE（LogP=2.00、d_LogP=1.0）。
%[text]    LogP がターゲット値 2.0 に完全一致し、最大の d_LogP スコアを得る。
%[text]    両者は異なる分子（same_molecule=0）であり、
%[text]    リード最適化の多目的テンションの本質を示す。
%[text]
%[text] PAINS フィルター統合（emk.filter.loadPainsSmarts が利用可能な場合に実行）
try
    pains_smarts = emk.filter.loadPainsSmarts();
    [D_sorted_f, rank_by_D] = sort(D, "descend");
    hasPains_top10 = false(1, 10);
    for k = 1:10
        idx_k = rank_by_D(k);
        hasPains_top10(k) = emk.mol.hasSubstruct(mols{idx_k}, pains_smarts);
    end
    n_pains_free = sum(~hasPains_top10);
    logInfo("やってみよう 6 -- 上位 10 候補中の PAINS フリー: %d / 10", n_pains_free);
catch
    logWarn("やってみよう 6 -- PAINS フィルター利用不可 (emk.filter.loadPainsSmarts)。");
    logInfo("            PAINS SMARTS データをインストールするか S03 ワークフローを使用。");
end

%[text] Q: 自分のリード化合物 ── テンプレート
%[text] カフェインをリードとした CNS ターゲット例（LogP 1-3、TPSA < 60、MW 200-400）:
%[text]   LEAD_SMILES = "Cn1cnc2c1c(=O)n(C)c(=O)n2C";
%[text]   TARGET_LOGP = 1.5; TARGET_TPSA = 45.0; TARGET_MW = 280.0;
%[text]   セクション 1 から再実行。
logInfo("A10 解答: 完了。");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
