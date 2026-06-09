%[text] # A09 解答: PFAS と環境化学物質スクリーニング
%[text] a09_pfas_screening.m の「やってみよう」演習の参照解答。
%[text] 最初に a09_pfas_screening.m（最低セクション 0〜5）を実行して
%[text] 必要なワークスペース変数を構築してください。その後このファイルで確認。
addpath(genpath("src"));
emk.setup.initPython();

mol_warmup = emk.mol.fromSmiles("C"); clear mol_warmup;
hasOptTbx   = license("test", "optimization_toolbox");
hasStatsTbx = license("test", "statistics_toolbox");
logInfo("A09 解答: セットアップ完了  (OptTbx=%d  StatsTbx=%d)", hasOptTbx, hasStatsTbx);
%%
%[text] ## 前提変数の再構築（a09 セクション 0～5 の再現）

CHEMICALS = { ...
    "PFBA",       "OC(=O)C(F)(F)C(F)(F)C(F)(F)F",                                   "PFCA"; ...
    "PFHxA",      "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                    "PFCA"; ...
    "PFOA",       "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",      "PFCA"; ...
    "PFNA",       "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "PFCA"; ...
    "PFDA",       "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "PFCA"; ...
    "PFBS",       "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                        "PFSA"; ...
    "PFHxS",      "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",         "PFSA"; ...
    "PFOS",       "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "PFSA"; ...
    "6:2FTS",     "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",               "FTS"; ...
    "8:2FTS",     "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "FTS"; ...
    "4:2FTS",     "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                              "FTS"; ...
    "Fluoxetine", "CNCCC(c1ccccc1)Oc1ccc(cc1)C(F)(F)F",                            "NonPFAS"; ...
    "Ciprofloxacin","OC(=O)c1cn(C2CCNCC2)c2cc(F)c(nc12)N1CCCC1",                   "NonPFAS"; ...
    "Flurbiprofen","OC(=O)C(C)c1ccc(cc1)-c1cccc(F)c1",                             "NonPFAS"; ...
    "Diflunisal",  "OC(=O)c1ccc(cc1O)-c1ccc(F)cc1F",                               "NonPFAS"; ...
    "Trifluridine","OC1C(CO)OC(n2cc(c(=O)[nH]2)C(F)(F)F)C1F",                      "NonPFAS"; ...
    "Halothane",   "FC(F)(F)C(Cl)Br",                                               "NonPFAS"; ...
    "Sevoflurane", "FC(F)(F)C(F)OCC(F)(F)F",                                        "NonPFAS"; ...
    "Flecainide",  "OCC(NC(=O)c1cc(OCC(F)(F)F)c(cc1OCC(F)(F)F)OCC(F)(F)F)CC",     "NonPFAS"; ...
    "Ethanol",     "CCO",                                                            "NonPFAS"; ...
};

nChem   = size(CHEMICALS, 1);
names   = string(CHEMICALS(:,1));
smiles  = string(CHEMICALS(:,2));
trueCls = string(CHEMICALS(:,3));

mols  = cell(1, nChem);
valid = false(1, nChem);
for k = 1:nChem
    try; mols{k} = emk.mol.fromSmiles(smiles(k)); valid(k) = true; catch; end
end

validMols  = mols(valid);
validNames = names(valid);
validSmiles = smiles(valid);
validCls   = trueCls(valid);
nValid     = sum(valid);

SMARTS_PFCHAIN    = "[#6](F)(F)-[#6](F)(F)";
SMARTS_CF3        = "[#6](F)(F)F";
SMARTS_SULFONYL   = "[#16](=O)(=O)";
SMARTS_POLYFLUORO = "[#6](F)(F)[#6H2]";

hasChain     = emk.mol.hasSubstruct(validMols, SMARTS_PFCHAIN);
hasCF3       = emk.mol.hasSubstruct(validMols, SMARTS_CF3);
hasSulfonyl  = emk.mol.hasSubstruct(validMols, SMARTS_SULFONYL);
hasPolyFluoro = emk.mol.hasSubstruct(validMols, SMARTS_POLYFLUORO);

isPFAS     = hasChain | (hasCF3 & hasChain);
classGuess = repmat("NonPFAS", nValid, 1);
for k = 1:nValid
    if isPFAS(k)
        if hasSulfonyl(k),      classGuess(k) = "PFSA";
        elseif hasPolyFluoro(k), classGuess(k) = "FTS";
        else,                    classGuess(k) = "PFCA";
        end
    end
end

DESCS       = ["LogP", "TPSA", "HeavyAtomCount"];
descTbl     = emk.descriptor.batchCalculate(validMols, DESCS);
logp_vec    = descTbl.LogP;
tpsa_vec    = descTbl.TPSA;
nF_vec      = zeros(nValid, 1);
for k = 1:nValid; nF_vec(k) = count(validSmiles(k), "F"); end

score_logP  = min(max(logp_vec / 10, 0), 1);
score_F     = min(max(nF_vec / 17,   0), 1);
score_TPSA  = exp(-((tpsa_vec - 40).^2) / (2 * 30^2));
scoresMat   = [score_logP, score_F, score_TPSA];

targetScores = zeros(nValid, 1);
for k = 1:nValid
    cls = validCls(k); nm = validNames(k);
    if cls == "NonPFAS",                              targetScores(k) = 0.10;
    elseif any(nm == ["PFOA","PFOS","PFNA","PFDA"]),  targetScores(k) = 0.80;
    elseif any(nm == ["PFHxA","PFHxS","8:2FTS"]),     targetScores(k) = 0.60;
    elseif any(nm == ["PFBA","PFBS","6:2FTS","4:2FTS"]),targetScores(k) = 0.40;
    else,                                             targetScores(k) = 0.50;
    end
end

if hasOptTbx
    objFun = @(w) sum((scoresMat * w - targetScores).^2);
    Aeq = [1,1,1]; beq = 1;
    lb = [0.05;0.05;0.05]; ub = [0.80;0.80;0.80]; w0 = [1/3;1/3;1/3];
    opts = optimoptions("fmincon", Display="off", Algorithm="sqp");
    [w_opt, ~] = fmincon(objFun, w0, [], [], Aeq, beq, lb, ub, [], opts);
else
    w_opt = [1/3;1/3;1/3];
end
concernScore = scoresMat * w_opt;

pfasHitIdx = find(isPFAS);
nHits      = numel(pfasHitIdx);

logInfo("前提条件再構築: %d 化学物質、%d PFAS ヒット。", nValid, nHits);
%%
%[text] ## やってみよう 1: PFOS vs PFHxS の炭素数; フルオキセチン非 PFAS; ハロタン OECD

logInfo("やってみよう 1 -- PFOS と PFHxS の炭素鎖長:");
for nm = ["PFOS", "PFHxS"]
    idx  = find(validNames == nm);
    fc   = parseSMILES_FCount(validSmiles(idx));
    mol  = validMols{idx};
    d    = emk.descriptor.calculate(mol, ["HeavyAtomCount","MolFormula"]);
    nCF2 = count(validSmiles(idx), "C(F)(F)");
    logInfo("  %s: formula=%s  heavy=%d  CF2 units in SMILES=%d", ...
        nm, d.MolFormula, d.HeavyAtomCount, nCF2);
end
%[text] 解答: PFOS（パーフルオロオクタンスルホン酸）はスルホニル頭部基に加えて 8 個の CF2 ユニットを持ち、
%[text]    フッ素化炭素数は計 8 個。PFHxS は 6 個。
%[text]    名前のプレフィックス: PF = パーフルオロ、Hx = ヘキサ（6）、S = スルホナート、
%[text]    O = オクタ（8）。この数字が炭素鎖長を表す。
logInfo("やってみよう 1 -- フルオキセチン: なぜ非 PFAS か?");
idx_fluox = find(validNames == "Fluoxetine");
mol_fluox = validMols{idx_fluox};
d_fluox   = emk.descriptor.calculate(mol_fluox, "MolFormula");
hasCF3_fluox  = emk.mol.hasSubstruct({mol_fluox}, SMARTS_CF3);
hasChain_fluox = emk.mol.hasSubstruct({mol_fluox}, SMARTS_PFCHAIN);
logInfo("  フルオキセチン式: %s  hasCF3=%d  hasPFchain=%d", ...
    d_fluox.MolFormula, hasCF3_fluox, hasChain_fluox);
%[text] 解答: フルオキセチンは孤立した -CF3 基（トリフルオロメチル）を 1 つ持つが、
%[text]    OECD PFAS 定義では少なくとも 1 つの -CF2-CF2- ユニット（連続した完全フッ素化炭素 2 個）が必要。
%[text]    CF3 単独ではパーフルオロアルキル鎖を形成しない。hasCF3 は true だが hasPFchain は false
%[text]    => NonPFAS に分類。
%[text]    重要な区別: 多くの医薬品は孤立した CF3 や CF2 基を含むが、
%[text]    PFAS のような環境残留性への懸念対象ではない。
logInfo("やってみよう 1 -- ハロタン (FC(F)(F)C(Cl)Br): OECD による PFAS か?");
idx_halo  = find(validNames == "Halothane");
mol_halo  = validMols{idx_halo};
hasCF3_halo   = emk.mol.hasSubstruct({mol_halo}, SMARTS_CF3);
hasChain_halo = emk.mol.hasSubstruct({mol_halo}, SMARTS_PFCHAIN);
logInfo("  ハロタン: hasCF3=%d  hasPFchain=%d  -> isPFAS=%d", ...
    hasCF3_halo, hasChain_halo, isPFAS(idx_halo));
%[text] 解答: ハロタンは 1 つの炭素上に CF3 を持つが、隣接炭素は Cl と Br を持ち F ではない。
%[text]    -CF2-CF2- 鎖は存在しない。OECD 2021 定義ではハロタンは PFAS ではない。
%[text]    OECD は「全ての C-F 結合がヘテロ原子に直接結合している物質」を除外するが、
%[text]    ハロタンは C-C 結合（CHClBr へ）を持つためこの除外は適用されない。
%[text]    正しい判断: パーフルオロアルキル鎖がない =>
%[text]    -CF3 が存在しても PFAS ではない。
%%
%[text] ## やってみよう 2: コントロールの hasCF3; SMARTS_CF3 の不十分さ; 改良 FTS SMARTS

logInfo("やってみよう 2 -- 非 PFAS コントロールの hasCF3:");
nonPFASidx = find(validCls == "NonPFAS");
for k = nonPFASidx'
    if hasCF3(k)
        logInfo("  %-18s: hasCF3=1  isPFAS=%d", validNames(k), isPFAS(k));
    end
end
%[text] 解答: フルオキセチン、トリフルリジン、ハロタン、セボフルラン、フレカイニドなどのコントロールは
%[text]    hasCF3 = true を返す。しかしパーフルオロアルキル鎖（SMARTS_PFCHAIN）を持たないため
%[text]    全て NonPFAS に分類される。hasCF3 単独の使用は PFAS 検出で
%[text]    非常に高い偽陽性率を引き起こす。
logInfo("やってみよう 2 -- 改良 FTS SMARTS  [#6](F)(F)-[#6H2]:");
hasFTSjunction = emk.mol.hasSubstruct(validMols, "[#6](F)(F)-[#6H2]");
ftsIdx = find(validCls == "FTS");
logInfo("  真の FTS 化合物:");
for k = ftsIdx'
    logInfo("    %-10s: hasFTSjunction=%d  classGuess=%s", ...
        validNames(k), hasFTSjunction(k), classGuess(k));
end
logInfo("  FTS ジャンクション SMARTS にヒットした非 FTS 化合物:");
for k = 1:nValid
    if hasFTSjunction(k) && validCls(k) ~= "FTS"
        logInfo("    %-18s (true=%s)", validNames(k), validCls(k));
    end
end
%[text] 解答: ジャンクションパターン [#6](F)(F)-[#6H2]（CF2 が CH2 に結合）は
%[text]    フルオロテロマー物質の構造的シグネチャー。
%[text]    6:2FTS、8:2FTS、4:2FTS を正しく同定し、PFCA（フッ素化鎖に CH2 なし、全体が -CF2-CF2-）や
%[text]    PFSA にはヒットしない。
%[text]    このより特異的なパターンにより FTS の PFCA への誤分類が減少する。
%%
%[text] ## やってみよう 3: 最高 F 数; エタノールプロキシスコア; TPSA プロキシ再設計

[maxFscore, maxFidx] = max(score_F);
logInfo("やってみよう 3 -- 最高 F 数スコア:");
logInfo("  化合物: %s  (nF=%d  score_F=%.3f)", ...
    validNames(maxFidx), nF_vec(maxFidx), maxFscore);
%[text] 解答: score_F は nF >= 17 で 1.0 に上限固定（正規化定数）。
%[text]    PFNA（nF=17）と PFDA（nF=19）は両方 score_F = 1.000 に達する。max() は最初の
%[text]    出現を返すため PFNA が勝者として表示される。
%[text]    実際にフッ素数が最も多いのは PFDA（19F: CF2 ユニット 9 個 + CF3 1 個 = 19F）。
%[text]    最長鎖の化合物が最も分解しにくい -- Buck et al.（2011）と整合。
%[text]    飽和を避けるため実際にはデータセットの nF_max で正規化することを推奨
%[text]    （例: nF/19）。
ethanolIdx = find(validNames == "Ethanol");
if ~isempty(ethanolIdx)
    logInfo("やってみよう 3 -- エタノールプロキシスコア:");
    logInfo("  nF=%d  score_logP=%.3f  score_F=%.3f  score_TPSA=%.3f", ...
        nF_vec(ethanolIdx), score_logP(ethanolIdx), ...
        score_F(ethanolIdx), score_TPSA(ethanolIdx));
end
%[text] 解答: エタノール（CCO）: フッ素原子なし => score_F = 0; LogP ~ -0.31 =>
%[text]    score_logP = 0（下限固定）; TPSA ~ 20 A^2 => 中心 40 のガウスで ~0.64
%[text]    -- これが唯一の非ゼロスコア。エタノールは持続性への懸念がなく陰性コントロールとして機能する。
%[text]    3 つのプロキシ全てが 0 に近いはずだが、TPSA = 20 がガウス中心（40 A^2）に近いため
%[text]    TPSA プロキシが小さなスコアを返す。
%[text]    これはこのプロキシ設計の既知の弱点。
logInfo("やってみよう 3 -- TPSA プロキシ再設計の議論:");
pfosIdx = find(validNames == "PFOS");
if ~isempty(pfosIdx)
    logInfo("  PFOS: TPSA=%.1f  score_TPSA=%.3f", tpsa_vec(pfosIdx), score_TPSA(pfosIdx));
end
%[text] 解答: PFOS は TPSA ~ 115 A^2（スルホナート頭部基）だが環境持続性は高い。
%[text]    中心 40 A^2 のガウスは PFOS に、短鎖 PFCA（TPSA ~ 37 A^2）より低い TPSA スコアを
%[text]    割り当てる。これは懸念評価として直感に反する。
%[text]    再設計案:
%[text]    (a) シグモイドを使用: score_TPSA = 1 / (1 + exp(-(TPSA - 60) / 10))
%[text]        -> 両極端で高懸念。
%[text]    (b) TPSA の代わりに溶解度予測（`emk.descriptor`: MolLogS）を使用。
%[text]    (c) TPSA をイオン性頭部基（-COOH / -SO3H）のバイナリフラグに置き換える。
%%
%[text] ## やってみよう 4: 最大重みプロキシ; PFOS ターゲット変更の影響; 第 4 プロキシ

logInfo("やってみよう 4 -- 最適化プロキシ重み:");
[maxW, maxWidx] = max(w_opt);
proxyNames = ["LogP", "F-count", "TPSA"];
logInfo("  w = [%.3f, %.3f, %.3f]  ->  最大: %s (%.3f)", ...
    w_opt(1), w_opt(2), w_opt(3), proxyNames(maxWidx), maxW);
%[text] 解答: 通常 w_F（F 数）が最大の重みを受け取る理由:
%[text]    (1) 長鎖 PFAS は nF が高く目標スコアも高い（0.7-0.8）。
%[text]    (2) 短鎖 PFAS と FTS は nF が中程度で目標スコアも中程度（0.4-0.6）。
%[text]    (3) 非 PFAS コントロールは nF = 0 で目標スコアが低い（0.1）。
%[text]    F 数プロキシは 3 つの懸念層を完全に分離するため、fmincon は最大の重みをここに割り当てる。
%[text]    これは Buck et al.（2011）と整合 -- C-F 結合密度が環境持続性の
%[text]    主要な要因として特定されている。
if hasOptTbx
    logInfo("やってみよう 4 -- PFOS ターゲットを 0.80 から 0.95 に上げる影響:");
    targetScores_mod          = targetScores;
    pfos_idx_local            = find(validNames == "PFOS");
    targetScores_mod(pfos_idx_local) = 0.95;
    objFun_mod = @(w) sum((scoresMat * w - targetScores_mod).^2);
    [w_opt_mod, ~] = fmincon(objFun_mod, w0, [], [], Aeq, beq, lb, ub, [], opts);
    logInfo("  元の重み: [%.3f, %.3f, %.3f]", w_opt(1), w_opt(2), w_opt(3));
    logInfo("  変更後の重み: [%.3f, %.3f, %.3f]", ...
        w_opt_mod(1), w_opt_mod(2), w_opt_mod(3));
    % 解答: PFOS のターゲットを上げると最適化器は F 数プロキシにさらに多くの重みを割り当てる
    %    （PFOS は 17 個の F 原子を持つ）。LogP の重みは通常ほぼ変わらないか減少する
    %    （PFOS の LogP はすでに高いため）。
    %    この感度分析は、規制上の懸念ルーブリックが直接プロキシ重みを決定することを示す
    %    —— 専門知識のエンコードの一形態。

    logInfo("やってみよう 4 -- 第 4 プロキシ: CF2 鎖長:");
    nCF2_vec = zeros(nValid, 1);
    for k = 1:nValid; nCF2_vec(k) = count(validSmiles(k), "C(F)(F)"); end
    score_CF2 = min(max(nCF2_vec / 10, 0), 1);  % normalise to ~10 CF2 units
    scoresMat4 = [score_logP, score_F, score_TPSA, score_CF2];
    Aeq4 = [1,1,1,1]; beq4 = 1;
    lb4 = [0.05;0.05;0.05;0.05]; ub4 = [0.80;0.80;0.80;0.80];
    w0_4 = [0.25;0.25;0.25;0.25];
    objFun4 = @(w) sum((scoresMat4 * w - targetScores).^2);
    [w_opt4, fval4] = fmincon(objFun4, w0_4, [], [], Aeq4, beq4, lb4, ub4, [], opts);
    logInfo("  4 プロキシ重み: LogP=%.3f  F=%.3f  TPSA=%.3f  CF2=%.3f  (RSS=%.4f)", ...
        w_opt4(1), w_opt4(2), w_opt4(3), w_opt4(4), fval4);
    % 解答: 専用の CF2 鎖長プロキシを追加すると通常残差（RSS）は減少する。
    %    鎖長プロキシと F 数プロキシは相関しているため、最適化器はそれらの間で
    %    重みを分散する。これは多重共線性があるとき最適な重み付けが一意でないことを示す。
end
%%
%[text] ## やってみよう 5: 最高 Tanimoto ペア; PFOA vs PFOS 類似度; 半径の影響

if nHits > 1
    fps = cell(1, nHits);
    for k = 1:nHits
        fps{k} = emk.fingerprint.morgan(validMols{pfasHitIdx(k)});
    end
    simMat = zeros(nHits, nHits);
    for i = 1:nHits
        for j = 1:nHits
            simMat(i,j) = emk.similarity.tanimoto(fps{i}, fps{j});
        end
    end
    hitNames = validNames(pfasHitIdx);

    logInfo("やってみよう 5 -- 最も類似した PFAS ペア:");
    triMask  = triu(true(nHits), 1);
    triVals  = simMat(triMask);
    [maxSim, maxLinIdx] = max(triVals);
    [ii, jj] = ind2sub([nHits, nHits], find(triMask));
    logInfo("  最大 Tanimoto = %.4f  %s と %s の間", ...
        maxSim, hitNames(ii(maxLinIdx)), hitNames(jj(maxLinIdx)));
    % 解答: 最も類似したペアは PFHxA（C6 PFCA）と PFOA（C8 PFCA）で T = 1.0。
    %    これはエラーではなく、Morgan/ECFP4 が繰り返し鎖でどう機能するかを反映する。
    %    両分子は -CF2- の繰り返しに -COOH を 1 つ加えた構造。
    %    半径 2 のビット列挙後、PFHxA のユニーク環状部分構造のセットは PFOA の厳密な部分集合
    %    （PFOA の余分な -CF2- ユニットは鎖が十分長くなると新しいビットパターンを追加しない）。
    %    Tanimoto = |A AND B| / |A OR B| = |PFHxA bits| / |PFOA bits| = 1.0
    %    PFHxA の全ビットが PFOA に現れ、PFOA が新しいビットを追加しないため。
    %    このフィンガープリント飽和は同族列の既知の制限（Maggiora et al. 2014）。
    %    より高い半径やカウントベースのフィンガープリント（ECFP6 等）で部分的に軽減できる。

    pfoa_idx  = find(hitNames == "PFOA");
    pfos_idx  = find(hitNames == "PFOS");
    if ~isempty(pfoa_idx) && ~isempty(pfos_idx)
        t_pfoa_pfos = simMat(pfoa_idx, pfos_idx);
        logInfo("やってみよう 5 -- PFOA vs PFOS Tanimoto = %.4f", t_pfoa_pfos);
    end
    % 解答: PFOA（C8 PFCA）と PFOS（C8 PFSA）は同じ炭素鎖長だが
    %    頭部基が異なる（-COOH vs -SO3H）。Tanimoto は通常 0.3～0.5:
    %    鎖の部分構造は重なるが、頭部基の部分構造（カルボン酸 vs スルホナート）が
    %    異なり類似度を下げる。
    %    半径 2 は 2 結合の近傍を捕える。半径 1 では頭部基の差異の影響が小さく
    %    T は 0.5～0.7 に上昇する。

    logInfo("やってみよう 5 -- 縮小半径（radius=1）の影響:");
    fps_r1 = cell(1, nHits);
    for k = 1:nHits
        fps_r1{k} = emk.fingerprint.morgan(validMols{pfasHitIdx(k)}, Radius=1);
    end
    if ~isempty(pfoa_idx) && ~isempty(pfos_idx)
        t_r1 = emk.similarity.tanimoto(fps_r1{pfoa_idx}, fps_r1{pfos_idx});
        logInfo("  PFOA vs PFOS T(radius=1) = %.4f  vs T(radius=2) = %.4f", ...
            t_r1, t_pfoa_pfos);
    end
    % 解答: 半径 1 ではアルゴリズムは 1 結合の近傍のみを考慮する。
    %    共有される -CF2-CF2- 鎖が強く寄与し、頭部基の差異は 1～2 原子に
    %    収まるため区別力が低下する。
    %    半径を小さくすると類似度は通常増加する。これは Morgan 半径2（ECFP4）が
    %    標準的な選択である理由を示す: 局所（結合）と中距離（部分構造）の
    %    構造情報をバランスよく捕えるため。
end
%%
%[text] ## やってみよう 6: 高懸念クラスター; 麻酔薬免除; CSV エクスポート

logInfo("やってみよう 6 -- 高懸念クラスター（LogP vs F 数）:");
highConcernMask = concernScore > 0.5 & isPFAS';
logInfo("  懸念 > 0.5 かつ isPFAS の化合物:");
for k = find(highConcernMask)'
    logInfo("    %-12s: LogP=%.2f  nF=%d  score=%.3f  [%s]", ...
        validNames(k), logp_vec(k), nF_vec(k), concernScore(k), validCls(k));
end
%[text] 解答: 長鎖 PFCA と PFSA（PFOA、PFOS、PFNA、PFDA、PFHxS）は LogP vs F 数の散布図の
%[text]    右上領域（高 LogP かつ高 nF）にクラスター化する。
%[text]    単純な 2 値スクリーンとして nF >= 10 の閾値はこのデータセットで偽陰性ゼロで
%[text]    全長鎖 PFAS を捕捉する。ただし F 数のみでは FTS
%[text]    （F は少ないが PFCA の前駆体）を見逃す。
logInfo("やってみよう 6 -- 麻酔薬免除設計:");
volatileMask  = zeros(nValid, 1);
for k = 1:nValid
    d_mw = emk.descriptor.calculate(validMols{k}, "MolWt");
    if d_mw.MolWt < 200
        volatileMask(k) = 1;
    end
end
exemptIdx = find(volatileMask & ~isPFAS');
logInfo("  データセット内の揮発性（MW<200）非 PFAS 化合物:");
for k = exemptIdx'
    logInfo("    %-14s: nF=%d  isPFAS=%d", validNames(k), nF_vec(k), isPFAS(k));
end
%[text] 解答: ハロタン（MW=197）とセボフルラン（MW=200）は揮発性麻酔ガス。
%[text]    「揮発性免除」フラグ（MW < 200）により、F 数スコアが中程度であっても
%[text]    高懸念層から除外される。これは現実の規制慣行を模倣:
%[text]    麻酔ガスは環境残留性フレームワークではなく職業暴露限値で規制される。
%[text]    免除ルーブリック: volatilityScore = (MW < 200) ? 0 : 1
%[text]    適用: finalScore = concernScore .* volatilityScore
%[text]
%[text] CSV エクスポート
[~, rankOrder] = sort(concernScore, "descend");
rankVec = zeros(nValid, 1);
rankVec(rankOrder) = (1:nValid)';
summaryTbl = table( ...
    rankVec, validNames, validCls, classGuess, ...
    round(logp_vec,2), nF_vec, round(tpsa_vec,1), round(concernScore,3), ...
    VariableNames=["Rank","Name","TrueClass","FlaggedClass","LogP","nF","TPSA_A2","ConcernScore"]);
summaryTbl = sortrows(summaryTbl, "Rank");

runDir = makeRunDir("Prefix", "a09_answers");
csvPath = fullfile(runDir, "pfas_screening_report.csv");
writetable(summaryTbl, csvPath);
logInfo("やってみよう 6 -- サマリーテーブルをエクスポート先: %s", csvPath);
disp(summaryTbl);
%[text] 解答: エクスポートした CSV は Excel や任意の分析ツールで読み込める。
%[text]    非専門家のステークホルダー向けに推奨する保持列:
%[text]    Rank、Name、FlaggedClass、ConcernScore。
%[text]    脚注で説明が必要な列: LogP（親油性）、nF（フッ素数）。
%[text]    TrueClass（内部検証用）と TPSA_A2 は不要な場合は削除。
%[text]    レポートの導入文例: "スコア > 0.5 の 5 種の化学物質が長鎖 PFAS として
%[text]    優先審査にフラグ付けされました（表 1）。"
logInfo("A09 解答: 完了。");
%%
%[text] ## ローカルヘルパー関数

function nF = parseSMILES_FCount(smi)
%[text] SMILES 文字列中の "F" の出現数を返す。
    nF = count(smi, "F");
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
