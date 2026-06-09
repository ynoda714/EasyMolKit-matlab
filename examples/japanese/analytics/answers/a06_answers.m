%[text] # A06 解答: 用量反応曲線フィッティング
%[text] a06_dose_response.m の「やってみよう」演習の参照解答。
%[text]
%[text] 注意: このファイルは Curve Fitting Toolbox のみ使用 -- RDKit / Python
%[text]       依存なし。`emk.setup.initPython()` は不要。
%[text]
%[text] `fittype()` と `fitoptions()` の expression 引数・パラメータ名引数には
%[text] すべて char リテラル（二重引用符の string ではなく）が必要です。
addpath(genpath("src"));
logInfo("A06 解答: セットアップ完了（Python 不要）");

%[text] 共有定数とヘルパー関数を（スタンドアロン用に）再定義します
TRUE_EMAX  = 100;
TRUE_EC50  = 1.0;
TRUE_N     = 1.5;
TRUE_EMIN  = 5;
NOISE_STD  = 8.0;
N_CONC     = 8;
rng(42);

hillModel = fittype( ...
    'Emin + (Emax - Emin) ./ (1 + (EC50 ./ x).^n)', ...
    'independent', 'x', ...
    'coefficients', {'Emax', 'n', 'EC50', 'Emin'});

hillFit_ = @(conc, resp) fit(conc(:), resp(:), hillModel, ...
    fitoptions(hillModel, ...
        'Lower',       [50,   0.1,  0.001, -10], ...
        'Upper',       [150,  10,   100,    50], ...
        'StartPoint',  [90,   1.0,  1.0,    10], ...
        'TolFun',      1e-8, ...
        'MaxIter',     2000, ...
        'Robust',      'off'));

conc = logspace(-3, 2, N_CONC)';   % 8 concentrations 0.001 to 100 uM
trueResp = @(c) TRUE_EMIN + (TRUE_EMAX - TRUE_EMIN) ./ ...
    (1 + (TRUE_EC50 ./ c).^TRUE_N);
resp = trueResp(conc) + NOISE_STD * randn(N_CONC, 1);
%%
%[text] ## やってみよう 1: Hill スロープの生物学; x=10*EC50 と x=0.1*EC50 での y; 中点の証明

logInfo("=== やってみよう 1: Hill スロープ n 値 ===");
logInfo("n = 1  : simple Michaelis-Menten (enzyme kinetics)");
logInfo("n = 2.8: ヘモグロビン O2 結合（協同的）");
logInfo("n = 1.5: 弱い協同性（我々の真のモデル）");

%[text] Hill 傾き (n) によるシグモイド形状の比較をプロットします
concFine = logspace(-2, 2, 200)';
figure("Name", "A06 Hill Slope Effect");
hold on;
for nVal = [0.5, 1, 1.5, 2, 4]
    yFine = TRUE_EMIN + (TRUE_EMAX - TRUE_EMIN) ./ (1 + (TRUE_EC50./concFine).^nVal);
    semilogx(concFine, yFine, LineWidth=1.5, DisplayName=sprintf("n=%.1f", nVal));
end
hold off;
xlabel("Concentration (uM)"); ylabel("Response (%)");
title("Hill 方程式の形状 vs Hill スロープ n");
legend(Location="southeast"); grid on;

%[text] x=10*EC50 と x=0.1*EC50 での応答を計算します
y10  = TRUE_EMIN + (TRUE_EMAX-TRUE_EMIN) / (1 + (TRUE_EC50/(10*TRUE_EC50))^TRUE_N);
y01  = TRUE_EMIN + (TRUE_EMAX-TRUE_EMIN) / (1 + (TRUE_EC50/(0.1*TRUE_EC50))^TRUE_N);
logInfo("y at x=10*EC50 : %.1f%%  (expected ~%d%% for Emax=%d)", ...
    y10, 100*(TRUE_EMAX-TRUE_EMIN)/(TRUE_EMAX-TRUE_EMIN+1)*TRUE_EMIN, TRUE_EMAX);
logInfo("y at x=0.1*EC50: %.1f%%  (close to Emin=%.1f)", y01, TRUE_EMIN);

%[text] x=EC50 で中点になることの代数的証明
%[text] Hill 方程式: y = Emin + (Emax-Emin)/(1+(EC50/x)^n)
%[text] x=EC50 のとき: 分母 = 1+(EC50/EC50)^n = 1+1 = 2
%[text] => y = Emin + (Emax-Emin)/2 = (Emax+Emin)/2
midpoint_analytical = (TRUE_EMAX + TRUE_EMIN) / 2;
midpoint_formula    = TRUE_EMIN + (TRUE_EMAX - TRUE_EMIN) / 2;
logInfo("Midpoint proof: Emin+(Emax-Emin)/2 = %.1f  ==  (Emax+Emin)/2 = %.1f", ...
    midpoint_formula, midpoint_analytical);

%[text] 解答: 大きな n（急峻なシグモイド）は、応答が狭い濃度範囲で急激に切り替わることを意味します。
%[text]    ヘモグロビン（n~2.8）は肺（高 pO2）での効率的な O2 結合と
%[text]    組織（低 pO2）での放出を可能にします。
%[text]    x=10*EC50（非常に高い用量）では、y は Emax に近づきますが完全には達しません。
%[text]    x=0.1*EC50（有効用量以下）では、y ≈ Emin になります。
%[text]    代数的に: x=EC50 のとき (EC50/x)^n = 1、分母 = 2 となり、
%[text]    y = Emin + (Emax-Emin)/2 = (Emax+Emin)/2 -- これがちょうど中点です。
%%
%[text] ## やってみよう 2: ノイズの影響; SEM を用いた三連デザイン

logInfo("=== やってみよう 2: ノイズと測定デザイン ===");
for noiseLevel = [0, 4, 8, 15, 25]
    rng(42);
    respN = trueResp(conc) + noiseLevel * randn(N_CONC, 1);
    try
        fN = hillFit_(conc, respN);
        logInfo("noiseStd=%-3d  -> fitted EC50=%.3f (true=%.3f)  error=%+.3f", ...
            noiseLevel, fN.EC50, TRUE_EC50, fN.EC50 - TRUE_EC50);
    catch ME
        logWarn("noiseStd=%-3d  -> fit failed: %s", noiseLevel, ME.message);
    end
end

%[text] 三連設計（濃度ごとに 3 回測定）
rng(42);
nReps = 3;
respTri = trueResp(conc) * ones(1, nReps) + NOISE_STD * randn(N_CONC, nReps);
semPerConc = std(respTri, 0, 2) / sqrt(nReps);
logInfo("三連 SEM（濃度間平均）: %.2f", mean(semPerConc));
logInfo("単一点ノイズ標準偏差: %.2f", NOISE_STD);

%[text] 三連の平均値でフィットします
respMean = mean(respTri, 2);
fTri = hillFit_(conc, respMean);
logInfo("Fit on triplicate means: EC50=%.3f (true=%.3f)", fTri.EC50, TRUE_EC50);

%[text] 解答: ノイズが大きいほど EC50 の信頼区間が広がります。noiseStd=25 では、
%[text]    フィットが収束しないか、生物学的に非合理なパラメータ
%[text]    （EC50 が境界値に当たる）が返る場合があります。
%[text]    三連測定は有効ノイズを sqrt(3) ≈ 1.73 倍だけ低減します。
%[text]    SEM = std/sqrt(n) は各濃度点の不確実性を定量化します。
%[text]    測定回数を増やすとフィット品質は向上しますが、試薬コストが 3 倍になります。
%%
%[text] ## やってみよう 3: 代替開始点; 無制約フィットの危険性

logInfo("=== やってみよう 3: 収束感度 ===");
startPoints = {[90, 1.0, 1.0, 10], [60, 2.0, 5.0, 0], [99, 0.5, 0.1, 1]};
for k = 1:numel(startPoints)
    sp = startPoints{k};
    try
        fk = fit(conc, resp, hillModel, ...
            fitoptions(hillModel, ...
                'Lower',      [50,  0.1, 0.001, -10], ...
                'Upper',      [150, 10,  100,    50], ...
                'StartPoint', sp, ...
                'TolFun',     1e-8, 'MaxIter', 2000));
        logInfo("StartPoint [Emax=%.0f n=%.1f EC50=%.1f Emin=%.0f] -> EC50=%.3f", ...
            sp(1), sp(2), sp(3), sp(4), fk.EC50);
    catch ME
        logWarn("StartPoint failed: %s", ME.message);
    end
end

%[text] 無制約フィットのリスクを確認します
try
    fUC = fit(conc, resp, hillModel, ...
        fitoptions(hillModel, 'MaxIter', 1000));
    logInfo("Unconstrained fit: EC50=%.3f  n=%.3f", fUC.EC50, fUC.n);
catch ME
    logWarn("Unconstrained fit failed: %s", ME.message);
end

%[text] 解答: 真値に近い適切な開始点は確実に収束します。
%[text]    不適切な開始点（真値が 1.0 のときに EC50=0.1 など）は局所最小値に
%[text]    収束するか、MaxIter 内に収束しない場合があります。
%[text]    無制約フィットは負の EC50（物理的に無意味）や n > 10（そのような急峻な
%[text]    シグモイドを生む既知の生物学的メカニズムは存在しない）を返す場合があります。
%[text]    fittype/fit には常に生理学的に意味のある境界を与えてください。
%[text]    Hill 方程式の導出: リガンド-受容体の急速平衡を仮定すると
%[text]    EC50^n = [L]^n * (1 - y)/y （y は占有率分率）
%[text]    => y = [L]^n / ([L]^n + EC50^n)
%%
%[text] ## やってみよう 4: N_CONC=5 での CI 幅; ゼロノイズ CI; predint の意味

logInfo("=== やってみよう 4: 信頼区間幅 ===");
for nPoints = [5, 8, 12, 20]
    rng(42);
    cTmp = logspace(-3, 2, nPoints)';
    rTmp = trueResp(cTmp) + NOISE_STD * randn(nPoints, 1);
    try
        fTmp = hillFit_(cTmp, rTmp);
        ciTmp = confint(fTmp, 0.95);
        width = ciTmp(2,3) - ciTmp(1,3);   % column 3 = EC50
        logInfo("N_CONC=%2d -> EC50 CI 幅 = %.3f", nPoints, width);
    catch
        logWarn("N_CONC=%2d -> fit failed", nPoints);
    end
end

%[text] ゼロノイズの場合を確認します
rng(42);
respZero = trueResp(conc);   % no noise
fZero  = hillFit_(conc, respZero);
ciZero = confint(fZero, 0.95);
logInfo("Zero noise CI for EC50: [%.4f, %.4f]", ciZero(1,3), ciZero(2,3));
logInfo("(Note: CI may be degenerate when residuals=0 -- tool-dependent)");

%[text] predint の比較: observation（観測値 PI）vs functional（平均曲線 PI）
fMain = hillFit_(conc, resp);
concPred = linspace(1e-3, 100, 50)';
predObs  = predint(fMain, concPred, 0.95, "observation");
predFun  = predint(fMain, concPred, 0.95, "functional");
logInfo("EC50 での予測区間幅（観測値 vs 関数値）:");
[~, idxEC50] = min(abs(concPred - TRUE_EC50));
logInfo("  Observation: %.3f  |  Functional: %.3f", ...
    predObs(idxEC50,2) - predObs(idxEC50,1), ...
    predFun(idxEC50,2) - predFun(idxEC50,1));

%[text] 解答: 濃度点が多いほど EC50 CI が狭まります（データがシグモイド形状をより強く制約するため）。
%[text]    8 点から 16 点に倍増するだけで CI 幅は概ね半分になります。
%[text]    ゼロノイズでは残差がゼロになり、アルゴリズムが雑音分散を推定できないため
%[text]    縮退した CI が返されます。
%[text]    "functional" predint は狭め（真の平均曲線の CI を与える）です。
%[text]    "observation" predint は測定雑音分散を加算するため広くなり、
%[text]    新しい単一測定値がどこに落ちるかを正確に表します。
%%
%[text] ## やってみよう 5: バッチ EC50 プロファイリング; 最悪の推定; 部分アゴニスト境界

logInfo("=== やってみよう 5: バッチ EC50 プロファイリング ===");

%[text] 4 つの仮想化合物（注意: この解答ファイル独自のデータ。教材の化合物 A-D とは EC50 値が異なります）
cpdNames    = ["Compound A", "Compound B", "Compound C", "Compound D"];
cpdEC50True = [0.5,  2.0, 0.1,  8.0];
cpdEmax     = [100,  100, 100,  60];   % Compound D is a partial agonist
cpdEmin     = [5,    5,   5,    5];
nCpd        = numel(cpdNames);
CPD_EMIN_CONST = 5;

rng(42);
ec50_fitted = nan(1, nCpd);
for c = 1:nCpd
    cpdConc = logspace(-3, 2, N_CONC)';
    cpdResp = CPD_EMIN_CONST + ...
              (cpdEmax(c) - CPD_EMIN_CONST) ./ ...
              (1 + (cpdEC50True(c) ./ cpdConc).^TRUE_N) + ...
              NOISE_STD * randn(N_CONC, 1);
    % For partial agonist, allow higher Upper bound on EC50 (less constrained)
    upperEC50 = ternary_(cpdEmax(c) < 100, 200, 100);
    try
        fc = fit(cpdConc, cpdResp, hillModel, ...
            fitoptions(hillModel, ...
                'Lower',      [max(cpdEmax(c)*0.5, 30),  0.1,  0.001, -10], ...
                'Upper',      [min(cpdEmax(c)*1.5, 150), 10,   upperEC50, 50], ...
                'StartPoint', [cpdEmax(c)*0.9, 1.0, 1.0, CPD_EMIN_CONST], ...
                'TolFun',     1e-8, 'MaxIter', 2000));
        ec50_fitted(c) = fc.EC50;
        logInfo("%-12s  true EC50=%.2f  fitted EC50=%.3f  error=%+.3f", ...
            cpdNames(c), cpdEC50True(c), fc.EC50, fc.EC50 - cpdEC50True(c));
    catch ME
        logWarn("%-12s  fit failed: %s", cpdNames(c), ME.message);
    end
end

%[text] 真値から最も外れた EC50 推定値を特定します
validFit = ~isnan(ec50_fitted);
[~, worstRel] = max(abs(ec50_fitted(validFit) - cpdEC50True(validFit)));
idxValid = find(validFit);
logInfo("真値から最も遠い EC50: %s", cpdNames(idxValid(worstRel)));

%[text] pEC50 を一括計算します（EC50: uM → M 変換後 pEC50 = -log10）
pec50 = -log10(ec50_fitted * 1e-6);
logInfo("pEC50 ランキング（高いほど強力）:");
[~, potOrd] = sort(pec50, "descend", "MissingPlacement", "last");
for k = 1:nCpd
    c = potOrd(k);
    if ~isnan(pec50(c))
        logInfo("  %d. %-12s  pEC50=%.2f", k, cpdNames(c), pec50(c));
    end
end

%[text] 解答: 化合物 D（部分アゴニスト、Emax=60%）は EC50 フィット誤差が大きくなる可能性があります。
%[text]    シグモイドの上部プラトーが不明瞭なため、上限漸近線を制約するデータ点が少なくなります。
%[text]    部分アゴニストには Emax の Upper 境界を広くとる必要があります
%[text]    （最適化器が早期にフィットをクリップするのを防ぐため）。
%[text]    pEC50 = -log10(EC50 in mol/L): pEC50 が高いほど EC50 が低く、効力が高いです。
%[text]    変換例: EC50=1.0 uM = 1e-6 M -> pEC50=6.0。
%[text]    化合物 C（EC50=0.1 uM = 1e-7 M）-> pEC50=7.0（このセット中で最も高効力）。
logInfo("A06 解答完了。");
%%
%[text] ## やってみよう 6: 可視化と単位変換（概念問題）
%[text]
%[text] やってみよう 6 は概念的・描画的質問のみです（新規コードは不要）。
%[text]
%[text] - **Q1**: フィット済み EC50（uM）を IC50（nM）に変換するには `* 1000` で十分です。
%[text]   セクション 4 の CI も同じ係数で変換できます（スケールが変わるだけで比率は不変）。
%[text]
%[text] - **Q2**: セクション 6 の全化合物曲線プロットと pEC50 バーチャートを参照してください。
%[text]   EC50 CI のエラーバーが大きい化合物ほど順位付けの信頼性が低くなります。
%[text]
%[text] - **Q3**: 協同結合（cooperative binding）は Hill 傾き n > 1 で数値化されます。
%[text]   アロステリック調節（allosteric modulation）は結合部位から離れた部位への
%[text]   リガンド結合が標的タンパク質の形状を変化させる機構です。
%[text]   ヘモグロビンの O2 結合（n≈2.8）が代表例として知られています。

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
