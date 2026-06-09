%[text] # A06: 用量反応曲線フィッティング
%[text] EasyMolKit アナリティクス — レイヤー 3
%[text] 
%[text] 「この薬の濃度を 10 倍にすると効果はどう変わるか」――創薬の現場では日常的に問われる質問です。その答えを数値化するのが **用量反応曲線** と **EC50**（半最大有効濃度）です。
%[text] Hill 方程式（4PL モデル）を用いることで、ノイズを含む実験データから EC50 を統計的に推定できます。
%[text] このスクリプトでは、4 つの合成化合物のアッセイデータを生成し、フィットして、`fit()` を用いて EC50 を 95% 信頼区間付きで求める過程を体験します。
%[text] 
%[text] **ストーリー**
%[text] ある薬理学者が 4 つの薬物候補の酵素阻害活性を試験しています。各化合物を様々な濃度で適用し、標的酵素の阻害割合を測定します。目的は各化合物の EC50（半最大有効濃度）を決定することです。これは創薬における最も重要な効力指標のひとつです。
%[text] Hill 方程式（4 パラメータロジスティック、4PL、またはシグモイド用量反応モデルとも呼ばれます）は、薬物濃度に伴う生物学的応答の変化を記述します。
%[text] 3 つの主要特性を捉えます:
%[text] - ベースライン応答（薬物なし）
%[text] - 最大応答（飽和薬物量）
%[text] - 遷移の急峻さ（Hill 傾き / 協同性） \
%[text] この演習では以下を行います:
%[text] 1.  異なる EC50 値（1000 倍にわたる）を持つ 4 化合物の現実的な用量反応データセットをシミュレートし、実際のアッセイ変動を模倣するための測定ノイズを加えます。
%[text] 2.  4PL モデルを使ったカスタム `fittype` で単一化合物をフィットします。
%[text] 3.  95% 信頼区間付きの EC50 と Hill 傾きを抽出します。
%[text] 4.  全化合物にわたるバッチフィッティングを実行し、結果を表にまとめます。
%[text] 5.  対数濃度軸で用量反応曲線を可視化し、効力で化合物を順位付けします。 \
%[text] 
%[text] **学習目標**
%[text] - Hill 方程式とその薬理学的解釈を理解する
%[text] - カスタム非線形モデル（`fittype`）で MATLAB の `fit()` を使用する
%[text] - 曲線フィッティングに物理的に意味のある下限 / 上限を設定する
%[text] - EC50 信頼区間を抽出・解釈する（`confint`）
%[text] - 標準的な用量反応プロット（対数スケール x 軸）を読み・作成する
%[text] - 化合物を効力で順位付けし、不確実性を議論する \
%[text] 
%[text] **注意**: この演習では Python / RDKit は不要です。
%[text] 純粋な MATLAB + Curve Fitting Toolbox です。
%[text] 
%[text] **前提条件**
%[text] - 基本的な統計知識（フィッティング、残差）
%[text] - Curve Fitting Toolbox（`fit`、`fittype`、`fitoptions`、`confint`、`predint`）
%[text] - インターネット接続不要 \
%[text] 
%[text] 推定所要時間: 30〜45 分
%[text] 
%[text] **参考文献**
%[text] - Hill AV (1910) The possible effects of the aggregation of the molecules of haemoglobin on its dissociation curves.*J Physiol* 40:iv-vii. (original Hill equation)
%[text] - Motulsky H & Christopoulos A (2004) Fitting Models to Biological Data Using Linear and Nonlinear Regression. Oxford University Press.(standard reference for dose-response fitting)
%[text] - Sebaugh JL (2011) Guidelines for accurate EC50/IC50 estimation.*Pharm Stat* 10:128-134. doi:10.1002/pst.426
%[text] - Ritz C, Baty F, Streibig JC, Gerhard D (2015) Dose-response analysis using R. *PLoS ONE* 10:e0146021. doi:10.1371/journal.pone.0146021 \
%[text] 
%[text] 実行方法: Ctrl+Enter でセクションを 1 つずつ実行
%%
%[text] ## セクション 0: セットアップ
%[text] **注意**: このセクションでは `emk.setup.initPython()` を使用しません。この演習は MATLAB のみで完結し、RDKit は使用しません。
% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
logSection("A06", "セクション 0: セットアップ", "アナリティクス L3");
logInfo("A06: セットアップ完了（Curve Fitting Toolbox 必須、Python 不要）");
%%
%[text] ## セクション 1: Hill 方程式と EC50
%[text] セットアップが完了しました。まず、この演習の中心となるHill方程式とEC50（半最大有効濃度）の意味を理解しましょう。形状の違いをグラフで確認し、直感を養いましょう。
%[text] ### コンセプト: 4パラメータロジスティック（4PL）Hill方程式
%[text] Hill方程式はシグモイド用量反応曲線を記述します。
%[text]{"align":"center"}   y = Emin + (Emax - Emin) / (1 + (EC50 / x)^n)
%[text] 
%[text] ### パラメータ:
%[text] - Emax  -- 最大応答（上部プラトー、例: 100%阻害）
%[text] - Emin  -- 最小応答（下部プラトー、例: ゼロ用量で0%）
%[text] - EC50  -- 半最大有効濃度 \[xと同じ単位\]
%[text] - x = EC50のとき: y = (Emax + Emin) / 2（中点）  n     -- Hill傾き（協同性係数）
%[text] -  n = 1: 単純な二分子結合（協同性なし）
%[text] - n \> 1: 正の協同性（急峻な遷移、例: イオンチャネル）
%[text] -  n \< 1: 負の協同性 / 混合機構（緩やかな遷移） \
%[text] ### EC50の解釈:
%[text] 低いEC50 = より高い効力（半最大効果に必要な薬物量が少ない）。  EC50はノイズのあるデータから推定されるため、常に95%信頼区間と共に報告します。 EC50 = 1 uM ± 0.5 uMの化合物はEC50 = 1 uM ± 10 uMと意味的に異なります。
%[text] ### IC50 vs EC50:
%[text] -   IC50: 50%阻害のための濃度（アンタゴニスト、酵素阻害剤）
%[text] -   EC50: 50%最大効果のための濃度（アゴニスト、活性化剤）
%[text] -   数学的には同一ですが、ラベルはアッセイの文脈によります。
%[text] -   この演習では酵素阻害を使用するため、EC50 = IC50です。 \
%[text] 
%[text] 濃度軸: 常にlog10スケールでプロットします。シグモイド形状（S字曲線）は対数スケールでのみ明確に現れます。線形スケールでは曲線は単一点での急激なステップとして現れます。
%[text] 直感を養うため、異なるHill傾き(n)でHill方程式の形状を比較します。
xDemo  = logspace(-2, 2, 200);   % concentration / EC50 (dimensionless)
nVals  = [0.5, 1.0, 2.0, 4.0];
colors = [0.8 0.2 0.2; 0.2 0.6 0.2; 0.2 0.4 0.9; 0.8 0.5 0.0];

figure("Name", "A06 Hill 方程式の形状");
for k = 1:numel(nVals)
    y_demo = 100 ./ (1 + (1 ./ xDemo).^nVals(k));
    semilogx(xDemo, y_demo, "-", LineWidth=2.0, Color=colors(k,:), ...
        DisplayName=sprintf("n = %.1f", nVals(k)));
    hold on;
end
xline(1.0, "k--", "EC50", LabelHorizontalAlignment="right", LineWidth=1.5, ...
    HandleVisibility="off");
yline(50,  "k:",  "50%",  LabelVerticalAlignment="bottom",  LineWidth=1.0, ...
    HandleVisibility="off");
xlabel("濃度 / EC50（対数スケール）"); ylabel("応答（%）");
title("異なる Hill 傾き (n) による Hill 方程式の形状");
legend(Location="northwest"); grid on; ylim([0 110]);

logInfo("A06: Hill 方程式デモをプロット。対数スケールで S 字曲線を確認しよう。");
%[text] **💡 観察ポイント 1**
%[text] n = 4 では遷移が非常に急峻（オール・オア・ナッシングに近い）ことを確認しましょう。
%[text] Hill傾き \> 3 を生じさせる生物学的メカニズムを考えてみましょう。
%[text] （ヒント: ヘモグロビンO2結合はn ~ 2.8; 電位依存性Na+チャネルの実効nは\> 4になり得る）
%[text] n = 1 の場合、x = 10*EC50での応答とx = 0.1*EC50での応答を確認しましょう。
%[text] 解析的に確認してみましょう: y = 100 / (1 + (EC50/x)^1)
%[text] なぜEC50は「半最大」と呼ばれるのかを考え、x = EC50のときに
%[text] 任意のnの値に対してy = (Emin + Emax)/2を代数的に示してみましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: 合成用量反応データの生成
%[text] 
%[text] Hill 方程式の理解が深まりました。次に、実際のアッセイを模倣した合成データを作成します。
%[text] ノイズを加えることで、実験におけるばらつきを再現します。
%[text] 
%[text] ### コンセプト: 現実的なアッセイ測定のシミュレーション
logSection("A06", "セクション 2: 合成用量反応データの生成", "アナリティクス L3");
%[text] 実際の用量反応アッセイでは、化合物ごとに 8〜12 濃度で生物活性を測定し、EC50 周辺の 3〜4 桁をカバーするために通常対数間隔（半対数または全対数ステップ）でサンプリングします。
%[text] 
%[text] 細胞ベースアッセイの測定ノイズは、通常、応答単位での加法的ガウスノイズ（~3～8% 標準偏差）として表現されます。
%[text] 酵素アッセイでは、5% の標準偏差が一般的です。
%[text] 
%[text] 4 つの合成化合物は、3 対数単位にわたる効力範囲を表します。
%[text] 
%[text] 化合物 A -- 高効力    （EC50 = 0.5  uM、n = 1.2）
%[text] 化合物 B -- 中程度の効力（EC50 = 5.0  uM、n = 1.0）
%[text] 化合物 C -- 低効力    （EC50 = 50.0 uM、n = 2.0）
%[text] 化合物 D -- 非常に低い（EC50 = 500  uM、n = 0.8、部分アゴニスト）
%[text] 
%[text] 部分アゴニスト（化合物 D）: Emax = 80% は、飽和濃度でも完全な阻害を達成できないことを意味します。これは、活性部位を部分的にしかブロックしない化合物によく見られます。
rng(2026);  % fix seed for reproducibility
%[text] 各化合物の真のパラメータを設定します。
compNames = ["Compound A", "Compound B", "Compound C", "Compound D"];
EC50_true = [0.5, 5.0, 50.0, 500.0];   % uM
n_true    = [1.2, 1.0,  2.0,   0.8];   % Hill slope
Emax_true = [100, 100, 100,   80.0];   % % response (top)
Emin_true = [  0,   0,   0,    5.0];   % % response (bottom)
noiseStd  = 5.0;                        % assay noise (% response units)

N_CONC    = 10;   % concentrations per compound
%[text] バッチフィッティング用に、全データを構造体に格納します。
data = struct();
for c = 1:numel(compNames)
    % Log-spaced concentrations: 2 log units below to 2 log units above EC50
    logMin = log10(EC50_true(c)) - 2;
    logMax = log10(EC50_true(c)) + 2;
    conc   = logspace(logMin, logMax, N_CONC)';  % uM, column vector

    % True Hill response + noise
    yTrue  = Emin_true(c) + (Emax_true(c) - Emin_true(c)) ./ ...
             (1 + (EC50_true(c) ./ conc).^n_true(c));
    yNoisy = yTrue + noiseStd * randn(N_CONC, 1);
    yNoisy = max(0, min(100, yNoisy));  % clamp to [0, 100]%

    data(c).name     = compNames(c);
    data(c).conc     = conc;     % observed concentrations (uM)
    data(c).response = yNoisy;  % observed responses (% inhibition)
    data(c).EC50     = EC50_true(c);
    data(c).n        = n_true(c);
    data(c).Emax     = Emax_true(c);
    data(c).Emin     = Emin_true(c);
end

logInfo("%d 化合物の用量反応データを生成（各 %d 点、ノイズ=%.0f%%）", ...
    numel(compNames), N_CONC, noiseStd);
%[text] データをざっと確認します（全 4 化合物を 1 つの対数スケールグラフに表示）。
figure("Name", "A06 生データ");
set(gcf, "Position", [100 100 900 520]);
tiledlayout(2, 2, "TileSpacing", "loose");
for c = 1:numel(compNames)
    nexttile;
    semilogx(data(c).conc, data(c).response, "o", MarkerSize=8, ...
        MarkerFaceColor=[0.3 0.5 0.9], MarkerEdgeColor="k");
    xlabel("濃度（uM、対数スケール）"); ylabel("阻害率（%）");
    title(sprintf("%s\nEC50 = %.1f uM（真値）", data(c).name, data(c).EC50));
    ylim([-5 110]); grid on;
end
sgtitle("A06: 生の用量反応データ（合成、ノイズ = 5%）");
%[text] **💡 観察ポイント 2**
%[text] noiseStd = 15.0（高ノイズアッセイ）を試して再プロットし、散布の変化を確認しましょう。
%[text] どの化合物の EC50 が最も正確に推定しにくいか、その理由を読み取りましょう。
%[text] 実際のアッセイでは、各濃度で 3 回の測定を行うことが多いです。
%[text] 各濃度に 3 つの反復測定（N\_CONC x 3）があるようにコードを修正しましょう。
%[text] 濃度ごとの平均 ± SEM をどのように計算し、表示するかを考えましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: 4PL Hill 方程式で単一化合物をフィットする
%[text] 
%[text] データが準備できました。まずは化合物 A を用いてフィッティングの流れを学びましょう。
%[text] `fit()` の使い方と初期値（StartPoint）および制約の設定方法を確認します。
%[text] 
%[text] ### コンセプト: `fit()` による非線形最小二乗フィッティング
%[text] MATLAB の `fit()`（Curve Fitting Toolbox）は、モデルと観測値の二乗残差和を最小化します。
%[text] {Emax, n, EC50, Emin} について最小化: sum\_i (y\_i - y\_model(x\_i))^2
%[text] 
%[text] Hill 方程式はパラメータに対して非線形であるため、最小化には反復アルゴリズム（デフォルトは Levenberg-Marquardt）を使用します。
%[text] このアルゴリズムは初期値に敏感です。
%[text] 
%[text] StartPoint: 各パラメータの初期推定値
%[text] Lower/Upper: 物理的制約を強制するボックス制約
%[text] 
%[text] 阻害率アッセイの物理的境界条件:
%[text] Emax  \[50, 120\]  -- 上部プラトー（~100% に近い; わずかな超過を許容）
%[text] n     \[0.1, 5.0\] -- Hill 傾き（実際に \> 4 になることはまれ）
%[text] EC50  \[1e-4, 1e6\] -- アッセイの濃度範囲（uM）
%[text] Emin  \[-10, 30\]  -- 下部プラトー（~0% に近い; ベースラインノイズを許容）
%[text] 
%[text] ワーキングサンプルとして化合物 A（インデックス 1）をフィットします。
c = 1;
xFit = data(c).conc;
yFit = data(c).response;

%[text] 4 パラメータ Hill モデルを fittype として定義します。
%[text] 係数の順序: {Emax, n, EC50, Emin}
%[text] fittype には式とパラメータ名に char リテラル（string ではなく）が必要です。
hillModel = fittype( ...
    'Emin + (Emax - Emin) ./ (1 + (EC50 ./ x).^n)', ...
    'independent', 'x', ...
    'coefficients', {'Emax', 'n', 'EC50', 'Emin'});

fitOpts = fitoptions(hillModel);
fitOpts.Lower      = [50,  0.1, 1e-4, -10];    % [Emax, n, EC50, Emin]
fitOpts.Upper      = [120, 5.0, 1e6,   30];
fitOpts.StartPoint = [100, 1.0, median(xFit), 0];

[fitResult, gof] = fit(xFit, yFit, hillModel, fitOpts);

%[text] フィットしたパラメータを抽出します。
ec50Fit  = fitResult.EC50;
nFit     = fitResult.n;
EmaxFit  = fitResult.Emax;
EminFit  = fitResult.Emin;

logInfo("化合物 A フィット: EC50=%.3f uM（真値=%.1f）  n=%.2f（真値=%.1f）", ...
    ec50Fit, data(c).EC50, nFit, data(c).n);
logInfo("               Emax=%.1f%%  Emin=%.1f%%  R2=%.4f  RMSE=%.2f%%", ...
    EmaxFit, EminFit, gof.rsquare, gof.rmse);

%[text] データにフィット曲線を重ねてプロットします。
xPlot  = logspace(log10(min(xFit)) - 0.5, log10(max(xFit)) + 0.5, 500);
yPlot  = EminFit + (EmaxFit - EminFit) ./ (1 + (ec50Fit ./ xPlot).^nFit);

figure("Name", "A06 単一フィット -- 化合物 A");
semilogx(xFit,  yFit,  "o", MarkerSize=10, MarkerFaceColor=[0.2 0.5 0.9], ...
    MarkerEdgeColor="k", DisplayName="Observed"); hold on;
semilogx(xPlot, yPlot, "b-", LineWidth=2.5, DisplayName="4PL Hill fit");
xline(ec50Fit, "--r", sprintf("EC50 = %.2f uM", ec50Fit), ...
    LineWidth=1.5, LabelHorizontalAlignment="right", HandleVisibility="off");
yline(50, ":k", LineWidth=1.0, HandleVisibility="off");
xlabel("濃度（uM、対数スケール）"); ylabel("阻害率（%）");
title(sprintf("A06: %s -- 4PL Hill フィット  (R^2=%.3f)", data(c).name, gof.rsquare));
legend(Location="northwest"); ylim([-5 110]); grid on;

%[text] **💡 観察ポイント 3**
%[text] StartPoint = \[100, 3.0, median(xFit), 0\] でフィットを実行した結果を確認しましょう。
%[text] アルゴリズムは同じ EC50 に収束しましたか？収束しなかった場合、その理由を考えてみましょう。
%[text] （ヒント: Levenberg-Marquardt アルゴリズムは局所最小値にトラップされる可能性があります）
%[text] Lower/Upper 制約を削除して再フィットした結果を確認しましょう。無制約の場合、結果はどう変わりますか？
%[text] Emax が非現実的に大きくなる可能性はありますか（例: 500%）？
%[text] なぜ生物学的用量反応フィッティングに制約が重要なのか考えてみましょう。
%[text] Hill 方程式はもともとヘモグロビンへの酸素結合のために導出されました。
%[text] その導出について調べ、Hill 傾きがリガンド結合部位の数と正確に等しくなる
%[text] 仮定について考えてみましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: EC50 の信頼区間
%[text] 
%[text] 単一化合物のフィッティングが完了しました。フィット値だけでなく、その不確かさ（95% 信頼区間）も必ず報告する必要があります。
%[text] EC50 の信頼区間を可視化し、推定の信頼性を評価します。
%[text] 
%[text] ### コンセプト: 非線形フィッティングにおける不確かさの定量化
logSection("A06", "セクション 4: EC50 の信頼区間", "アナリティクス L3");
%[text] `confint()` は最小値での二乗和曲面の曲率に基づいて、各フィットパラメータの 95% 信頼区間を計算します。
%[text] 狭い CI はよく決定されたパラメータを示します。一方、広い CI はデータがパラメータを正確に固定できないことを意味します。
%[text] 
%[text] CI の幅は以下の要因に依存します:
%[text] - データ点数（多いほど狭い CI）
%[text] - ノイズレベル（低いほど狭い CI）
%[text] - データのカバレッジ（EC50 付近に測定点が多い場合、EC50 の CI が狭くなる）
%[text] - モデルの識別可能性（全パラメータが独立して推定可能かどうか） \
%[text] 
%[text] 目安: log10(EC50) の 95% CI が 1 単位を超える場合（10 倍の不確実範囲）、EC50 推定値は信頼できません。
%[text] アッセイを繰り返すか、濃度範囲を拡張しましょう。
ci = confint(fitResult, 0.95);   % 2 x 4 matrix: [lower; upper] for [Emax, n, EC50, Emin]
ec50Lo = ci(1, 3);
ec50Hi = ci(2, 3);
nLo    = ci(1, 2);
nHi    = ci(2, 2);

logInfo("化合物 A -- 95%% CI:");
logInfo("  EC50 = %.3f uM  [%.3f, %.3f]  (不確実性 %.1f 倍)", ...
    ec50Fit, ec50Lo, ec50Hi, ec50Hi / max(ec50Lo, 1e-9));
logInfo("  n    = %.2f     [%.2f, %.2f]", nFit, nLo, nHi);

%[text] 95% 予測区間バンド付きでプロットします。
yPI = predint(fitResult, xPlot, 0.95, "functional", "off");  % 曲線の 95% PI

figure("Name", "A06 信頼区間バンド -- 化合物 A");
fill([xPlot, fliplr(xPlot)], [yPI(:,1)', fliplr(yPI(:,2)')], ...
    [0.7 0.8 1.0], FaceAlpha=0.4, EdgeColor="none", DisplayName="95% PI band");
hold on;
semilogx(xFit,  yFit,  "o", MarkerSize=10, MarkerFaceColor=[0.2 0.5 0.9], ...
    MarkerEdgeColor="k", DisplayName="観測値");
semilogx(xPlot, yPlot, "b-", LineWidth=2.5, DisplayName="4PL フィット");
xline(ec50Fit, "--r", sprintf("EC50 = %.2f [%.2f, %.2f] uM", ec50Fit, ec50Lo, ec50Hi), ...
    LineWidth=1.5, LabelHorizontalAlignment="right", ...
    LabelVerticalAlignment="bottom", HandleVisibility="off");
xlabel("濃度（uM、対数スケール）"); ylabel("阻害率（%）");
title(sprintf("A06: %s -- 95%% 信頼区間バンド", data(c).name));
legend(Location="northwest"); ylim([-10 115]); grid on;
set(gca, "XScale", "log");

logInfo("95%% 予測区間バンドを表示。");

%[text] **💡 観察ポイント 4**
%[text] N\_CONC = 5（データ点数の半分）でシミュレーションを再実行した場合、
%[text] EC50 の 95% CI の幅がどのように変わるかを確認しましょう。
%[text] 通常推奨される最小濃度数は何かを読み取りましょう。
%[text] （参考: Sebaugh 2011 は信頼できる 4PL フィッティングに \>= 10 を推奨しています）
%[text] noiseStd = 0 を設定して再フィットした場合、CI がほぼゼロに縮小するか確認しましょう。
%[text] なぜ実際には CI が正確にゼロにならないのかを考察しましょう。
%[text] "functional" の predint は平均曲線の CI を提供します。
%[text] "observation" に変更して単一の新しい観測値の CI を取得し、その違いを確認しましょう。
%[text] なぜ観測 PI が常に機能的 PI より広いのかを説明しましょう。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: 全化合物のバッチフィッティング
%[text] 
%[text] 単一化合物のフィッティング手順を学びました。これを全4化合物に自動適用し、効力を比較します。結果はサマリーテーブルで確認しましょう。
%[text] 
%[text] ### コンセプト: 系統的な効力プロファイリング
%[text] 実際の創薬キャンペーンでは、同じフィッティング手順を数十から数百の化合物に適用します。
%[text] バッチフィッティングにより、各化合物のEC50値とその信頼区間を自動的に抽出します。
%[text] 
%[text] 結果は効力に基づいて化合物を順位付けしたサマリーテーブルで報告します（最低EC50 = 最高効力）。
%[text] 信頼区間（CI）の倍率幅が信頼性を示します。
%[text] 
%[text] CI 倍率幅 = EC50\_upper / EC50\_lower
%[text] \< 2 倍:   優れた精度
%[text] 2〜5 倍:  許容範囲（ノイズの多い細胞アッセイで典型的）
%\[text]   > 10 倍:  信頼できません。濃度範囲を拡張するか、アッセイを繰り返しましょう。
fitSummary = struct();

for c = 1:numel(compNames)
    xc = data(c).conc;
    yc = data(c).response;

    % 適応的 StartPoint: EC50 を濃度の幾何平均から開始
    startEC50 = exp(mean(log(xc)));

    fopts = fitoptions(hillModel);
    fopts.Lower      = [50,  0.1, 1e-4, -10];
    fopts.Upper      = [120, 5.0, 1e6,   30];
    fopts.StartPoint = [100, 1.0, startEC50, 0];

    try
        [fr, gf] = fit(xc, yc, hillModel, fopts);
        ci_c     = confint(fr, 0.95);

        fitSummary(c).Name      = data(c).name;
        fitSummary(c).EC50      = fr.EC50;
        fitSummary(c).EC50_lo   = ci_c(1, 3);
        fitSummary(c).EC50_hi   = ci_c(2, 3);
        fitSummary(c).HillSlope = fr.n;
        fitSummary(c).Emax      = fr.Emax;
        fitSummary(c).R2        = gf.rsquare;
        fitSummary(c).RMSE      = gf.rmse;
        fitSummary(c).fitObj    = fr;

        logInfo("%-12s  EC50=%7.2f uM [%6.2f, %7.1f]  n=%.2f  R2=%.3f", ...
            data(c).name, fr.EC50, ci_c(1,3), ci_c(2,3), fr.n, gf.rsquare);
    catch ME
        logWarn("%s のフィット失敗: %s", data(c).name, ME.message);
        fitSummary(c).Name = data(c).name;
        fitSummary(c).EC50 = NaN;
    end
end

%[text] サマリーテーブルを構築します。
ec50Vec = [fitSummary.EC50];
loVec   = [fitSummary.EC50_lo];
hiVec   = [fitSummary.EC50_hi];
nVec    = [fitSummary.HillSlope];
r2Vec   = [fitSummary.R2];

summTbl = table(compNames(:), EC50_true(:), ec50Vec(:), loVec(:), hiVec(:), ...
    n_true(:), nVec(:), r2Vec(:), ...
    VariableNames=["Compound","EC50_true_uM","EC50_fit_uM", ...
                   "CI_lower","CI_upper","n_true","n_fit","R2"]);
disp(summTbl);

%[text] **💡 観察ポイント 5**
%[text] どの化合物でフィットしたEC50が真値から最も離れているかを確認しましょう。
%[text] 真値は常に95%信頼区間内にあるかを確認しましょう。
%[text] そうでない場合、これはバグを示すかを考えてみましょう。（ヒント: 95% CIはすべての実験での
%[text] 包含を保証しない -- 平均的に95%の時間のみ）
%[text] 化合物DはEmax = 80%（部分阻害）。部分アゴニストのフィッティングを
%[text] 適切に許可するために下限をどう修正するかを考えましょう。
%[text] （現在の制約はEmax \>= 50%を強制し、EC50推定にバイアスをかける可能性があります）
%[text] バッチフィッティングループを拡張してpEC50 = -log10(EC50\_uM)も計算し、
%[text] pEC50が効力を線形スケールに変換するため薬化学でよく使われることを確認しましょう。
%[text] （pEC50 = 6はEC50 = 1 uMを意味します）。
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: 効力比較と多パネル可視化
%[text] 
%[text] 全化合物の EC50（半数効果濃度）が得られました。最後に、全化合物の用量反応曲線を重ねて視覚的に比較し、効力で順位付けを行います。
%[text] pEC50 スケールでのランキングバーチャートも作成します。
%[text] 
%[text] ### コンセプト: 化合物の効力の可視化と順位付け
%[text] 用量反応曲線は対数スケールの濃度軸でプロットします。このスケールでは、よくフィットした 4PL Hill 曲線は対称な S 字形（シグモイド）として現れ、EC50 は変曲点に対応します。
%[text] 
%[text] **効力の順位付け:** EC50 は昇順で報告します。EC50 分布は対数正規なので、線形スケールでは報告しません。報告例: `化合物 A: EC50 = 0.5 uM（95% CI: 0.3 — 0.9 uM）`
%[text] 
logSection("A06", "セクション 6: 効力比較と多パネル可視化", "アナリティクス L3");
%[text] 全化合物の用量反応曲線を 1 つのグラフに重ね合わせてプロットします。
cmap = lines(numel(compNames));

figure("Name", "A06 全用量反応曲線");
xPlotGlobal = logspace(-3, 4, 500);  % cover full concentration range (uM)

for c = 1:numel(compNames)
    if isnan(fitSummary(c).EC50), continue; end
    fr = fitSummary(c).fitObj;
    yc = fr.Emin + (fr.Emax - fr.Emin) ./ ...
         (1 + (fr.EC50 ./ xPlotGlobal).^fr.n);
    semilogx(xPlotGlobal, yc, "-", LineWidth=2.5, Color=cmap(c, :), ...
        DisplayName=sprintf("%s  EC50=%.1f uM", compNames(c), fr.EC50));
    hold on;
    % EC50 を曲線上にマーク（中点 = (Emin+Emax)/2 が正確な y 座標）
    semilogx(fr.EC50, (fr.Emax + fr.Emin) / 2, "v", MarkerSize=12, MarkerFaceColor=cmap(c, :), ...
        MarkerEdgeColor="k", HandleVisibility="off");
    % 生データ点（MarkerFaceAlpha は R2020b+; 互換性のために set() を使用）
    h_pts = semilogx(data(c).conc, data(c).response, "o", MarkerSize=6, ...
        MarkerFaceColor=cmap(c, :), MarkerEdgeColor="k", ...
        HandleVisibility="off");
    try; set(h_pts, "MarkerFaceAlpha", 0.6); catch; end
end
yline(50, "--k", "50%（EC50 閾値）", LabelHorizontalAlignment="left", LineWidth=1.0, HandleVisibility="off");
xlabel("濃度（uM、対数スケール）"); ylabel("阻害率（%）");
title("A06: 全化合物 -- フィットした用量反応曲線");
legend(Location="east"); ylim([-5 110]); grid on;

%[text] 効力順位棒グラフ（pEC50 スケール）を作成します。
pEC50_true = -log10(EC50_true * 1e-6);  % convert uM to M first
pEC50_fit  = -log10(ec50Vec  * 1e-6);

[~, sortIdx] = sort(pEC50_fit, "descend");  % rank by potency (highest pEC50 first)

figure("Name", "A06 効力ランキング");
colors_bar = cmap(sortIdx, :);
bh = bar(pEC50_fit(sortIdx), FaceColor="flat", DisplayName="フィット pEC50");
bh.CData = colors_bar;
hold on;
%[text] CI エラーバーを追加します（対数スケールで非対称 → 対数空間では対称）。
pLo = -log10(hiVec(sortIdx) * 1e-6);  % note: EC50_hi -> lower pEC50
pHi = -log10(loVec(sortIdx) * 1e-6);  % note: EC50_lo -> higher pEC50
errorbar(1:numel(compNames), pEC50_fit(sortIdx), ...
    pEC50_fit(sortIdx) - pLo, pHi - pEC50_fit(sortIdx), ...
    ".k", LineWidth=1.5, DisplayName="95% CI");
plot(1:numel(compNames), pEC50_true(sortIdx), "r^", MarkerSize=10, ...
    MarkerFaceColor="r", DisplayName="真の pEC50");
set(gca, "XTickLabel", compNames(sortIdx), "XTickLabelRotation", 20);
ylabel("pEC50 = -log_{10}[EC50 (M)]"); title("A06: 化合物効力ランキング");
legend(Location="northeast"); grid on;
logInfo("効力ランキング完了。pEC50 が高いほど効力が高い。");

%[text] **💡 観察ポイント 6**
%[text] フィットした EC50 値（uM）を IC50（nM）に変換して確認しましょう。
%[text] （1 uM = 1000 nM）どの化合物が「ナノモル効力」と呼ばれるかを確認しましょう。
%[text] 薬物候補は通常 IC50 \< 1 uM（pEC50 \> 6）であれば前進します。
%[text] 化合物 D は小さい Emax（部分阻害）を持つ。これは pEC50 推定にどう影響するかを考えましょう。
%[text] 完全阻害剤と並べて化合物 D を pEC50 で順位付けするのは公平かを考えましょう。
%[text] 追加で報告すべき指標は何かを考えましょう。
%[text] 実際のキャンペーンでは化合物はしばしば繰り返し実験で再試験されます。
%[text] 化合物 A がある実行で EC50 = 0.3 uM、別の実行で EC50 = 0.9 uM を与えた場合、
%[text] これらの結果は一致しているかを確認しましょう。
%[text] セクション 4 で計算した 95% CI と比較してみましょう。
%[text] Hill 方程式は n = 1 の場合に単一の結合部位を仮定します。
%[text] 化合物 C が n ~ 2 を持つ理由を説明できる構造情報は何かを考えましょう。
%[text] 「協同結合」と「アロステリック調節」を調べてみましょう。
% ... （ここにコードを書いてみましょう）

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
