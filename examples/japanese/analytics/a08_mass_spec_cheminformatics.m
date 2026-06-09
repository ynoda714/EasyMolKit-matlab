%[text] # A08: 質量分析 × ケモインフォマティクス
%[text] EasyMolKit アナリティクス — レイヤー 3
%[text]
%[text] 製薬 QC ラボに 5 つの薬サンプルが届きました。ラベルが水で傷んで読めなくなっています。
%[text] HPLC を使わずに、手元の LC-MS データだけで「これは何の薬か」を特定できるでしょうか？
%[text] 質量分析（MS）は分子の「精密な質量」を ppm レベルで測定します。
%[text] 同位体パターン（M、M+1、M+2 のスペクトル指紋）を活用すると、質量が似た候補を 1 つに絞り込むことができます。
%[text] このスクリプトでは、低分解能検索、高精密質量、同位体確認の 3 段階で化合物を同定するプロセスを体験します。
%[text]
%[text] **ストーリー**
%[text]
%[text] ある製薬 QC アナリストが倉庫から 5 つの薬サンプルのバッチを受け取りました。
%[text] ラベルが水で傷んで読めなくなっています。
%[text] 各サンプルはすでに実験室の LC-MS 装置で測定されており、生の質量スペクトルが得られています。
%[text] アナリストはバッチをリリースする前に、各サンプルの同定を確認しなければなりません。
%[text]
%[text] 同定プロセスは、確信度が上昇する 3 段階で展開されます。
%[text]
%[text] - **段階 1 — 低分解能検索**（単位質量許容差 0.5 Da）: 実際の商用データベース（数万件以上）では数十の化合物がマッチしますが、本演習の 200 化合物 DB では 1〜2 件になります。
%[text] - **段階 2 — 高分解能精密質量**（5 ppm 許容差）: わずかな候補のみが生き残ります。有望ですが、決定的ではありません。
%[text] - **段階 3 — 同位体パターン確認**: 特徴的な M・M+1・M+2 の指紋が各未知物質を高い信頼度で 1 つの化合物に絞り込みます。
%[text]
%[text] **演習内容**
%[text]
%[text] 1. `emk.descriptor.calculate()` を使って、200 の FDA 承認薬から精密質量参照テーブルを構築します。
%[text] 2. 5 つの薬物から 5 つのリアルな ESI-MS スペクトルをシミュレートします（ガウスピーク形状、同位体クラスター、ランダムノイズ）。
%[text] 3. Signal Processing Toolbox（`smoothdata` + `findpeaks`）を使って、生スペクトルのピークを検出します。
%[text] 4. 質量許容差で参照テーブルを検索し、低分解能と高分解能での候補数を比較します。
%[text] 5. 分子式から理論的同位体パターンを計算します（13C / 34S / 37Cl 存在比からの M・M+1・M+2 相対強度）。
%[text] 6. 各候補をその同位体パターンと観測スペクトルのコサイン類似度でスコアリングし、最終ランク付きリストを作成します。
%[text]
%[text] **学習目標**
%[text]
%[text] - 精密質量測定が化合物 ID の曖昧さを減らす理由を説明できます。
%[text] - スペクトルデータに `findpeaks` と `smoothdata`（Signal Processing Toolbox）を使用できます。
%[text] - ppm 質量誤差を計算し、適切な検索許容差を設定できます。
%[text] - 自然元素存在比から M+1 と M+2 の同位体強度を導き出せます。
%[text] - コサイン類似度で候補をスコアリングし、ランク付きヒットリストを構築できます。
%[text] - MS 化合物 ID を S05 のフィンガープリントベースアプローチと結びつけて理解できます。
%[text]
%[text] **前提条件**
%[text]
%[text] - A03（QSAR 回帰）修了 ―― `emk.descriptor.calculate()` の基本を理解していること。
%[text] - S05（未知物質 ID）修了 ―― 類似度ベース同定の考え方を理解していること。
%[text] - Signal Processing Toolbox（`findpeaks`、`smoothdata`）を使用できること。
%[text] - Statistics and Machine Learning Toolbox（`corrcoef`）を使用できること。
%[text] - インターネット接続は不要です。
%[text]
%[text] **動作環境**
%[text]
%[text] セクション 3・4 には Signal Processing Toolbox が必要です。
%[text] セクション 5 では Statistics and ML Toolbox（`corrcoef`）を使用します。
%[text] 両ツールボックスは MATLAB Online Basic（無料枠）に含まれています。
%[text]
%[text] 推定所要時間: 35〜50 分
%[text]
%[text] **データ:**
%[text]
%[text] `data/list/fda_drugs.csv` ―― 200 FDA 承認薬（ChEMBL、CC-BY-SA 3.0）
%[text] 列: ChEMBLID、Name、SMILES、MolecularWeight、ALogP、HBondDonors、HBondAcceptors、TPSA、RotatableBonds
%[text]
%[text] **参考文献**
%[text]
%[text] Gross JH (2011) *Mass Spectrometry: A Textbook*, 2nd ed. Springer.
%[text] ISBN 978-3-642-10709-2 ―― 単一同位体質量、同位体パターン、ESI-MS
%[text]
%[text] Kind T & Fiehn O (2006) Metabolomic database annotations via query of elemental compositions.
%[text] *BMC Bioinformatics* 7:234. doi:10.1186/1471-2105-7-234 ―― 化合物同定における質量精度の実際的限界
%[text]
%[text] Claesen J et al. (2012) Efficient method for isotopic distribution calculation.
%[text] *J Am Soc Mass Spectrom* 23:753-763. doi:10.1007/s13361-011-0326-2
%[text]
%[text] Stein SE & Scott DR (1994) Optimization and testing of mass spectral library search algorithms.
%[text] *J Am Soc Mass Spectrom* 5:859-866. doi:10.1016/1044-0305(94)87009-8 ―― コサイン類似度スコアリング
%[text]
%[text] Meija J et al. (2016) Atomic weights of the elements 2013 (IUPAC Technical Report).
%[text] *Pure Appl. Chem.* 88:265-291. doi:10.1515/pac-2015-0305 ―― 自然同位体存在比: 13C 1.103%、34S 4.25%、37Cl 24.23%
%[text]
%[text] 実行方法: Ctrl+Enter でセクションを 1 つずつ実行
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

%[text] メイン処理を行う前に、PythonとRDKitのプロセスをウォームアップします
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;
logSection("A08", "セクション 0: セットアップ", "アナリティクス L3");
%%
%[text] ## セクション 1: 精密質量参照テーブルの構築
%[text]
%[text] セットアップが完了しました。まず、200のFDA承認薬から精密質量参照テーブルを構築します。
%[text] これがデータベース検索の「辞書」として機能します。
%[text]
%[text] ### コンセプト: 単一同位体質量 vs 平均分子量
%[text]
%[text] 周期表には各元素に対して2種類の「質量」が記載されています:
%[text]
%[text] 平均分子量（MolWt）:
%[text] 自然存在比を用いた全安定同位体の加重平均です。
%[text] 例: 炭素 = 12.011 g/mol（12Cは98.9%、13Cは1.1%）
%[text] 化学者が化学量論や溶液調製に使用します。
%[text]
%[text] 精密単一同位体質量（ExactMolWt）:
%[text] 各元素の最も存在比の高い同位体のみで構成された分子の質量です。
%[text] 例: 12C、1H、14N、16O、32S、35Cl、など
%[text] 例: 炭素 = 正確に12.000000 Da（質量スケールの基準！）
%[text] 質量分析装置は個々のイオンの質量電荷比（m/z）を測定するため、
%[text] 質量分析で使用されます。
%[text]
%[text] なぜこの区別が重要なのでしょうか？
%[text] アスピリン（C9H8O4）: MolWt = 180.16 g/mol、ExactMolWt = 180.0423 Da
%[text] この差（0.12 Da）は、誤った質量タイプで検索するとデータベース不一致を
%[text] 引き起こすのに十分な大きさです。
%[text]
%[text] [M+H]+ イオン（ポジティブモードESI-MSでのプロトン化分子）の場合:
%[text] 観測されるm/z = ExactMolWt + 1.00728（プロトン質量）
%[text] （1.00728 Da = プロトン質量; 電子質量 ~0.00055 Da は無視できる）
logSection("A08", "セクション 1: 精密質量参照テーブルの構築", "アナリティクス L3");
DATA_FILE = "data/list/fda_drugs.csv";

rawTbl = readtable(DATA_FILE, TextType="string");
nRaw   = height(rawTbl);
logInfo("%s から %d 行を読み込みました", DATA_FILE, nRaw);

%[text] 各化合物の精密質量と分子式を計算します。
exactMass = nan(nRaw, 1);
formula   = strings(nRaw, 1);
valid     = false(1, nRaw);

logInfo("精密質量を計算中（1〜2 分かかる場合があります）...");
for k = 1:nRaw
    try
        mol          = emk.mol.fromSmiles(rawTbl.SMILES(k));
        d            = emk.descriptor.calculate(mol, ["ExactMolWt","MolFormula"]);
        exactMass(k) = d.ExactMolWt;
        formula(k)   = d.MolFormula;
        valid(k)     = true;
    catch ME
        logWarn("行 %d（%s）: %s", k, rawTbl.Name(k), ME.message);
    end
end

validIdx = find(valid);
refTbl = table( ...
    rawTbl.Name(validIdx), ...
    rawTbl.SMILES(validIdx), ...
    formula(validIdx), ...
    exactMass(validIdx), ...
    VariableNames=["Name","SMILES","Formula","ExactMass"]);
nRef = height(refTbl);
logInfo("参照テーブル: %d 化合物", nRef);

%[text] プレビュー
disp(refTbl(1:5, ["Name","Formula","ExactMass"]));

%[text] **💡 観察ポイント 1 — アスピリンと別の薬の精密質量を確認しましょう**
%[text] アスピリンの精密単一同位体質量を確認しましょう。
%[text] 実行: mol = `emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O")`; 
%[text] d   = `emk.descriptor.calculate(mol, "ExactMolWt")`; 
%[text] fprintf("アスピリン ExactMolWt = %.4f Da\n", d.ExactMolWt)
%[text] ポジティブモードESIスペクトルで期待される [M+H]+ m/z を確認しましょう。
%[text] （プロトンに1.00728 Daを加えます。）
%[text] イブプロフェンの質量スペクトルでベースピークがm/z 207.14に現れます。
%[text] これはイブプロフェン（C13H18O2）の [M+H]+ と一致するか確認しましょう。
%[text] 確認: `emk.descriptor.calculate(emk.mol.fromSmiles("CC(C)Cc1ccc(C(C)C(=O)O)cc1"), "ExactMolWt")`
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: 5 つの未知 ESI-MS スペクトルのシミュレーション
%[text]
%[text] 参照テーブルが完成しました。次に、5 つの未知薬物のリアルな ESI-MS スペクトルをシミュレートします。
%[text] 実際の装置が測定するスペクトルを再現し、同定演習の準備を整えます。
%[text]
%[text] ### コンセプト: ESI-MS スペクトルと同位体クラスター
%[text]
%[text] - **[M+H]+** イオンを中心に M、M+1、M+2 の同位体クラスターが形成されます。
%[text] - ベースラインノイズ、H2O、CO2、およびその他のフラグメントイオンが含まれます。
%[text] - シミュレーションは固定シードで再現可能です。
logSection("A08", "セクション 2: 5 つの未知 ESI-MS スペクトルのシミュレーション", "アナリティクス L3");
PROTON_MASS  = 1.007276;   % Da
MZ_SIGMA_LR  = 0.3;        % Da -- シミュレートした低分解能ピーク幅（FWHM ~0.7 Da）
NOISE_LEVEL  = 0.03;       % ベースピークに対する割合
N_UNKNOWNS   = 5;

%[text] 参照テーブル全体に均等に分布した N_UNKNOWNS 薬物を選択します。
%[text] （固定選択 -- rng に依存しません）
unknownIdx = round(linspace(ceil(nRef * 0.05), floor(nRef * 0.95), N_UNKNOWNS));

unknownNames    = refTbl.Name(unknownIdx);
unknownFormulas = refTbl.Formula(unknownIdx);
unknownMass     = refTbl.ExactMass(unknownIdx);

logInfo("選択した未知物質:");
for k = 1:N_UNKNOWNS
    logInfo("  未知物質 %d: %s  (%.4f Da、分子式 %s)", ...
        k, unknownNames(k), unknownMass(k), unknownFormulas(k));
end

%[text] スペクトルを生成します（struct 配列: .mz、.intensity）。
rng(2026);   % 再現可能なノイズのための固定シード
spectra = cell(1, N_UNKNOWNS);
for k = 1:N_UNKNOWNS
    mH    = unknownMass(k) + PROTON_MASS;  % [M+H]+ m/z
    fc    = parseFormula(unknownFormulas(k));
    iso   = isoPattern(fc);                 % [M, M+1, M+2] relative intensities

    % Peak positions and relative intensities
    pkMz  = [mH,          mH+1,       mH+2, ...   % 同位体クラスター
             mH-18,       mH-44];                  % 一般的な中性脱離
    pkInt = [iso(1),      iso(2),     iso(3), ...
             iso(1)*0.35, iso(1)*0.15];

    % 有効な m/z 範囲（50 Da 以上）のピークのみを保持
    keep  = pkMz > 50;
    spectra{k} = simulateSpectrum(pkMz(keep), pkInt(keep), ...
                     [max(50, mH-80), mH+10], MZ_SIGMA_LR, NOISE_LEVEL);
end

%[text] 未知物質 1 のスペクトルを表示します。
figure("Name","A08: サンプルスペクトル（未知物質 1）");
plot(spectra{1}.mz, spectra{1}.intensity, "b-", LineWidth=0.8);
xlabel("m/z"); ylabel("相対強度");
title(sprintf("シミュレート ESI-MS スペクトル -- 未知物質 1（ラベル非表示）"));
grid on;

%[text] **💡 観察ポイント 2 — 別の未知物質のスペクトルを観察してみましょう**
%[text] 未知物質 3 のスペクトルを確認しましょう。フィギュアインデックスを 3 に変えて
%[text] 再実行し、分子イオンクラスターを視覚的に識別できるか確認しましょう。
%[text] （スペクトルの高 m/z 側の最も高いピーククラスターを見つけましょう。）
%[text] ESI スペクトルの「分子イオン領域」は通常、最高 m/z の
%[text] クラスターです（一価電荷イオンの場合）。x 軸の範囲を確認しましょう。
%[text] 未知物質 1 の [M+H]+ はおよそ何 m/z でしょうか？
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: Signal Processing Toolbox によるピーク検出
%[text]
%[text] スペクトルが準備できました。次に Signal Processing Toolbox を使用してノイズを除去し、ピークを自動検出します。
%[text] これは [M+H]+ の m/z を正確に読み取るための前処理です。
%[text]
%[text] ### コンセプト: findpeaks と smoothdata
%[text]
%[text] - `smoothdata(..., "gaussian")` を用いて高周波ノイズを低減し、ピーク形状を歪めずに検出精度を向上させます。
%[text] - `findpeaks` の主要なパラメータ: MinPeakHeight、MinPeakProminence、MinPeakDistance
%[text] - 同位体ピークは 1 Da 間隔であるため、MinPeakDistance は 1 Da 未満で設定する必要があります。
logSection("A08", "セクション 3: Signal Processing Toolbox によるピーク検出", "アナリティクス L3");
MIN_PEAK_HEIGHT      = 0.05;   % 正規化強度の 5%
MIN_PEAK_PROMINENCE  = 0.04;
MIN_PEAK_DISTANCE    = 0.8;    % Da -- 同位体ピークを分解するには 1 未満が必要

%[text] 未知物質 1 でのデモンストレーション
sp1       = spectra{1};
smoothInt = smoothdata(sp1.intensity, "gaussian", 5);

[pkHeight, pkLoc] = findpeaks(smoothInt, sp1.mz, ...
    "MinPeakHeight",     MIN_PEAK_HEIGHT, ...
    "MinPeakProminence", MIN_PEAK_PROMINENCE, ...
    "MinPeakDistance",   MIN_PEAK_DISTANCE);

%[text] ピークの強度を正規化します。
pkHeight = pkHeight / max(pkHeight);

logInfo("未知物質 1: %d ピークを検出", numel(pkLoc));
if numel(pkLoc) > 0
    logInfo("  検出された m/z 値: %s", strjoin(string(round(pkLoc, 3)), ", "));
end

%[text] 分子イオン [M+H]+ は高 m/z クラスター内で最も強いピークです（最高検出 m/z の
%[text] 3 Da 以内）。max(pkLoc) のみを使用すると、Cl/S 化合物では M+2 同位体ピークが選択され、
%[text] 例としてフロセミドの M+2 は約 39% であり、中性質量の逆算に
%[text] 系統的な 2 Da の誤差が生じます。
highMzMask  = pkLoc >= max(pkLoc) - 3;
clusterHts  = pkHeight(highMzMask);
clusterLocs = pkLoc(highMzMask);
[~, iM]     = max(clusterHts);
mH_obs_1    = clusterLocs(iM);
logInfo("未知物質 1: 抽出された [M+H]+ = %.4f Da", mH_obs_1);

%[text] 生データ、平滑化データ、検出されたピークを可視化します。
figure("Name","A08: ピーク検出（未知物質 1）");
plot(sp1.mz, sp1.intensity, Color=[0.7 0.7 0.7], DisplayName="生データ"); hold on;
plot(sp1.mz, smoothInt,     "b-", LineWidth=1.2, DisplayName="平滑化");
plot(pkLoc,  pkHeight,      "rv", MarkerSize=8, MarkerFaceColor="r", ...
     DisplayName="検出ピーク");
xlabel("m/z"); ylabel("相対強度");
title("未知物質 1 -- 生データ、平滑化、検出ピーク");
legend; grid on;

%[text] findpeaks を使用して、全未知物質の概算 [M+H]+ を抽出します。
%[text] グリッドステップ ~0.018 Da は精度を ~50-100 ppm に制限しますが、
%[text] 0.5 Da の低分解能検索には十分で、5 ppm の高分解能検索には不十分です。
mH_obs = nan(1, N_UNKNOWNS);
for k = 1:N_UNKNOWNS
    sp = spectra{k};
    sm = smoothdata(sp.intensity, "gaussian", 5);
    [pks, locs] = findpeaks(sm, sp.mz, ...
        "MinPeakHeight",     MIN_PEAK_HEIGHT, ...
        "MinPeakProminence", MIN_PEAK_PROMINENCE, ...
        "MinPeakDistance",   MIN_PEAK_DISTANCE);
    if ~isempty(locs)
        highMzMask  = locs >= max(locs) - 3;   % 分子イオンクラスターウィンドウ
        clusterLocs = locs(highMzMask);
        [~, iM]     = max(pks(highMzMask));     % M ピークは M+2 より強い
        mH_obs(k)   = clusterLocs(iM);
    end
end

%[text] キャリブレーションされた高分解能質量読み取りをシミュレートします（セクション 4〜6）。
%[text] 実際の Orbitrap/Q-TOF は内部ロックマスキャリブレーションにより ~1 ppm rms で
%[text] 質量を報告します。これは、グリッド制限された上記の findpeaks セントロイドとは根本的に異なり、
%[text] 質量アナライザー（スペクトル画像ではない）が精度を決定します。
%[text] モデル化: 真の [M+H]+ に Gaussian(0, 1 ppm) ノイズを加えたものです。
rng(99, "twister");
mH_hires = (unknownMass + PROTON_MASS)' + ...
           (unknownMass + PROTON_MASS)' .* randn(1, N_UNKNOWNS) * 1e-6;

%[text] **💡 観察ポイント 3 — ピーク検出パラメータを調整してみましょう**
%[text] MinPeakHeight を 0.01 に設定した場合、どのような変化があるか確認しましょう。
%[text] 偽のノイズピークがいくつ現れるかを観察し、0.20 に設定した場合の変化も確認しましょう。
%[text] MinPeakDistance を 0.1 に設定してみましょう。これにより findpeaks は同位体ピークを
%[text] 個別に分解します。M、M+1、M+2 のトリプレットを識別できるか確認しましょう。
%[text] （[M+H]+、[M+H]+1、[M+H]+2 の 3 つのピークが連続しているかを確認します。）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: 精密質量データベース検索
%[text]
%[text] ピークの m/z 値が取得できました。次に、参照テーブルと照合して候補を絞り込みます。
%[text] 低分解能と高分解能での候補数の違いを比較してみましょう。
%[text]
%[text] ### コンセプト: ppm（百万分率）での質量精度
%[text]
%[text] 質量分析計は質量精度において大きく異なります。
%[text]
%[text] 低分解能（シングル四重極、イオントラップ）:
%[text] m/z 300 での精度は約0.2〜0.5 Daで、これは約700〜1700 ppmに相当します。
%[text] 多くの薬物が同じ整数の名目質量を共有しています。
%[text] 例: 名目質量300は40以上のFDA承認薬に一致します。
%[text]
%[text] 高分解能（Orbitrap、Q-TOF、FT-ICR）:
%[text] 精度は5 ppm未満（m/z 300で0.0015 Da未満）です。
%[text] 候補を数十から1〜2つに絞り込むことができます。
%[text]
%[text] ppm誤差: （観測質量 - 理論質量） / 理論質量 * 1e6
%[text]
%[text] 検索には2つの異なる装置を表す2つの質量変数を使用します。
%[text]
%[text] mH_obs(k) -- シミュレートされたスペクトルからのfindpeaksセントロイド（約50〜100 ppm）
%[text] 低分解能のDaウィンドウ検索に使用します。
%[text] mH_hires(k) -- シミュレートされたOrbitrap読み取り: 真の[M+H]+に1 ppmのノイズを加えたもの
%[text] 高分解能のppmウィンドウ検索に使用します。
%[text]
%[text] この分離は物理的に正確です。単位分解能検出器はスペクトルピーク位置を約0.5 Daの精度で報告しますが、
%[text] 高分解能装置の内部質量キャリブレーションは約1〜5 ppmの正確な質量を提供します。
%[text]
%[text] LOW_RES_TOL = 0.50 Da（単位分解能装置をシミュレート）
%[text] HIGH_RES_PPM = 5.0 ppm（Orbitrap / Q-TOFをシミュレート）
logSection("A08", "セクション 4: 精密質量データベース検索", "アナリティクス L3");
LOW_RES_TOL  = 0.50;   % Da
HIGH_RES_PPM = 5.0;    % ppm

candidateCounts = zeros(N_UNKNOWNS, 2);  % [低分解能, 高分解能]

for k = 1:N_UNKNOWNS
    % 低分解能検索（Da ウィンドウ）-- findpeaks 抽出を使用（~50〜100 ppm）
    lowResCand = [];
    if ~isnan(mH_obs(k))
        mNeutral_lr = mH_obs(k) - PROTON_MASS;
        daDiff      = abs(refTbl.ExactMass - mNeutral_lr);
        lowResCand  = find(daDiff <= LOW_RES_TOL);
    end

    % 高分解能検索（ppm ウィンドウ）-- キャリブレーション済み装置読み取りを使用（~1 ppm）
    mNeutral_hr = mH_hires(k) - PROTON_MASS;
    ppmDiff     = abs(refTbl.ExactMass - mNeutral_hr) ./ refTbl.ExactMass * 1e6;
    highResCand = find(ppmDiff <= HIGH_RES_PPM);

    candidateCounts(k, :) = [numel(lowResCand), numel(highResCand)];

    logInfo("未知物質 %d（高分解能 [M+H]+ = %.4f、中性 = %.4f Da）:", ...
        k, mH_hires(k), mNeutral_hr);
    logInfo("  低分解能 (+-%.2f Da): %d 候補", LOW_RES_TOL, numel(lowResCand));
    logInfo("  高分解能 (%.1f ppm): %d 候補", HIGH_RES_PPM, numel(highResCand));
end

%[text] 候補数の棒グラフを表示します。
figure("Name","A08: 分解能別候補数");
bar(1:N_UNKNOWNS, candidateCounts);
legend(sprintf("低分解能 (+-%.2f Da)", LOW_RES_TOL), ...
       sprintf("高分解能 (%.0f ppm)", HIGH_RES_PPM));
xlabel("未知サンプル"); ylabel("データベース候補数");
title("候補削減: 低分解能 vs 高分解能 MS");
grid on;
maxCand = ceil(max(candidateCounts(:)));
yticks(0:maxCand + 1);
ylim([0, maxCand + 1.5]);

%[text] **💡 観察ポイント 4 — ppm 許容差を変えて候補数の変化を確認しましょう**
%[text] 未知物質1でHIGH_RES_PPMを50に変更した場合、候補はいくつになるか確認しましょう。
%[text] 次に1 ppmを試して、真の化合物が残るか確認しましょう。
%[text] （mH_hiresは約1 ppmのrmsノイズを持ちます。1 ppmの閾値では、ノイズの影響で
%[text] 1〜2つの未知物質がウィンドウのギリギリ外に落ちる場合があります。）
%[text] Kind & Fiehn（2006）は、1 ppmでも質量だけでは区別できない分子式があると示しました。
%[text] 候補をさらに絞り込むために使用できる追加情報は何か確認しましょう。
%[text] （ヒント: セクション5を先読みしてみましょう。）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: 同位体パターンスコアリング
%[text]
%[text] 精密質量によって候補が絞られました。次に、同位体パターンを用いて最終的な同定を行います。
%[text] コサイン類似度スコアを用いて候補をランク付けします。
%[text]
%[text] ### コンセプト: 構造フィンガープリントとしての同位体クラスター
%[text]
%[text] すべての元素は、自然界に存在する同位体の特徴的な分布を持っています。
%[text] M、M+1、M+2 ピークの相対強度は分子の元素組成を示します。
%[text]
%[text] 主要な寄与（一価電荷有機薬物分子の場合）:
%[text]
%[text] M+1 相対存在比（M の %）:
%[text] 13C: 炭素原子 1 個あたり 1.103%（支配的な寄与）
%[text] 15N: 窒素原子 1 個あたり 0.366%
%[text] 2H:  水素原子 1 個あたり 0.015%
%[text] 17O: 酸素原子 1 個あたり 0.038%
%[text] --> M+1% ~ 1.10*nC + 0.37*nN + 0.015*nH + 0.04*nO
%[text]
%[text] M+2 相対存在比（M の %）:
%[text] 13C^2 の寄与（二項式）: (1.10*nC)^2 / 200
%[text] 18O: 酸素原子 1 個あたり 0.205%
%[text] 34S: 硫黄原子 1 個あたり 4.25%（大きい -- S は識別子！）
%[text] 37Cl: 塩素原子 1 個あたり ~32.0%（非常に大きい -- Cl は高度に識別的）
%[text] 導出: M+2/M = p(37Cl)/p(35Cl) = 24.23%/75.77% = 31.98%
%[text] Gross（2011）表 3.2 では 32.7% に近似（差 < 2%）
%[text] 注意: 24.23% は 37Cl の自然存在比（IUPAC）; 32.0% は
%[text] M+2 の相対強度（35Cl 単一同位体ピークに正規化）
%[text] --> M+2% ~ (1.10*nC)^2/200 + 0.21*nO + 4.25*nS + 32.7*nCl
%[text]
%[text] Cl または S を含む分子は異常に大きい M+2 ピークを示します。
%[text] Cl 原子 1 つを持つ薬物は M:M+2 比が約 3:1 を示し、
%[text] スペクトルで間違えようがありません。
%[text]
%[text] スコアリング戦略: 観測同位体ベクトル [I_M, I_{M+1}, I_{M+2}]（正規化）と
%[text] 各候補の理論ベクトルのコサイン類似度。
%[text]
%[text] 各未知物質: スペクトルから観測同位体クラスターを抽出し、
%[text] 全高分解能候補を同位体コサイン類似度でスコアリングします。
logSection("A08", "セクション 5: 同位体パターンスコアリング", "アナリティクス L3");
isoScoreTbls = cell(1, N_UNKNOWNS);

for k = 1:N_UNKNOWNS
    % 候補選択と同位体ウィンドウターゲティングにキャリブレーション質量を使用
    mNeutral   = mH_hires(k) - PROTON_MASS;
    ppmDiff    = abs(refTbl.ExactMass - mNeutral) ./ refTbl.ExactMass * 1e6;
    candIdx    = find(ppmDiff <= HIGH_RES_PPM);
    if isempty(candIdx), continue; end

    % 平滑化スペクトルからの観測同位体強度（M、M+1、M+2）
    sp     = spectra{k};
    smInt  = smoothdata(sp.intensity, "gaussian", 5);
    smInt  = smInt / max(smInt);   % normalise to 1

    obsIso = zeros(1, 3);
    for ishift = 0:2
        targetMz       = mH_hires(k) + ishift;
        [~, bestPt]    = min(abs(sp.mz - targetMz));
        obsIso(ishift+1) = max(smInt(max(1,bestPt-3):min(end,bestPt+3)));
    end
    if obsIso(1) < 1e-6, obsIso(1) = 1; end  % ゼロ除算ガード
    obsIso = obsIso / obsIso(1);              % M ピークに相対して正規化

    % 各候補をスコアリング
    nCand  = numel(candIdx);
    scores = nan(nCand, 1);
    for ci = 1:nCand
        fc        = parseFormula(refTbl.Formula(candIdx(ci)));
        theoIso   = isoPattern(fc);   % [1, M+1/M, M+2/M]
        % コサイン類似度
        scores(ci) = dot(obsIso, theoIso) / ...
                     (norm(obsIso) * max(norm(theoIso), 1e-12));
    end

    [scoresSorted, sortOrd] = sort(scores, "descend");
    isoScoreTbls{k} = table( ...
        refTbl.Name(candIdx(sortOrd)), ...
        refTbl.Formula(candIdx(sortOrd)), ...
        refTbl.ExactMass(candIdx(sortOrd)), ...
        ppmDiff(candIdx(sortOrd)), ...
        scoresSorted, ...
        VariableNames=["Name","Formula","ExactMass","ppmError","IsoScore"]);
end

%[text] 未知物質 1 の同位体スコアランキングを表示します
logInfo("未知物質 1 -- 同位体スコアでランク付けした高分解能候補:");
if ~isempty(isoScoreTbls{1})
    disp(isoScoreTbls{1}(1:min(5,height(isoScoreTbls{1})), :));
end

%[text] **💡 観察ポイント 5 — 同位体スコアで最上位候補を確認しましょう**
%[text] 未知物質 1 の同位体スコアテーブルを見て、最上位の候補が
%[text] 真の同定に一致するか確認しましょう。
%[text] （真の同定は: unknownNames(1) -- セクション 2 で表示されています。）
%[text] FDA データベースで Cl 原子を少なくとも 1 つ含む化合物を探し、
%[text] その M+2 相対強度を予測してみましょう。
%[text] 次に確認: それが refTbl に現れるか確認しましょう。
%[text] （ヒント: Formula 列で "Cl" を探す:
%[text] refTbl(contains(refTbl.Formula, "Cl"), :)）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 6: 完全同定ワークフロー
%[text]
%[text] 同位体スコアの計算が完了しました。最後に、全3段階を統合して完全な同定ワークフローを実行します。
%[text] 5つの未知物質を正しく同定できるか確認しましょう。
%[text]
%[text] ### コンセプト: 質量 + 同位体証拠の組み合わせによる化合物ID
%[text]
%[text] 完全ワークフローを順番に適用します:
%[text] ステップ 1: 生スペクトルから [M+H]+ を抽出（findpeaks）
%[text] ステップ 2: 中性精密質量を逆算（プロトンを差し引く）
%[text] ステップ 3: 高分解能質量検索（5 ppm ウィンドウ）
%[text] ステップ 4: 残る候補を同位体スコアでランク付け（コサイン類似度）
%[text] ステップ 5: 最上位候補を同定として選択
%[text]
%[text] 最終同定は同定精度を計算するために、既知の真の同定
%[text] （unknownNames に格納）と比較します。
%[text]
%[text] S05（未知物質同定）との関連:
%[text] S05 は構造情報データベースから未知化合物を同定するために
%[text] フィンガープリント類似度（Tanimoto）を使用しました。
%[text] A08 は物理的シグナル（質量スペクトル）を使って、未知物質の
%[text] 構造知識を事前に必要とせずに同じ目標を達成します。
%[text] 実際には、両アプローチは補完的です: MS は精密質量と同位体式を提供し、
%[text] FP 類似度は参照ライブラリの SMILES が利用可能な場合に一致を確認します。
logSection("A08", "セクション 6: 完全同定ワークフロー", "アナリティクス L3");
logInfo("全同定ワークフローを実行中...");
nCorrect = 0;
resultRows = cell(N_UNKNOWNS, 1);

for k = 1:N_UNKNOWNS
    if isempty(isoScoreTbls{k})
        logWarn("未知物質 %d: 同定に十分なデータがない", k);
        continue;
    end

    predicted = isoScoreTbls{k}.Name(1);   % top-ranked candidate
    trueID    = unknownNames(k);
    correct   = strcmpi(predicted, trueID);
    if correct, nCorrect = nCorrect + 1; end

    resultRows{k} = {k, trueID, predicted, ...
                     isoScoreTbls{k}.ppmError(1), ...
                     isoScoreTbls{k}.IsoScore(1), correct};

    resultLabel = " 不正解"; if correct, resultLabel = " 正解"; end
    logInfo("未知物質 %d:  真=%s  予測=%s  [%s]", ...
        k, trueID, predicted, resultLabel);
end

%[text] サマリーテーブル
validRows = find(~cellfun(@isempty, resultRows));
if ~isempty(validRows)
    resultTbl = cell2table(vertcat(resultRows{validRows}), ...
        VariableNames=["Unknown","TrueID","Predicted","ppmError","IsoScore","Correct"]);
    disp(resultTbl);
    logInfo("同定精度: %d / %d = %.0f%%", ...
        nCorrect, numel(validRows), nCorrect/numel(validRows)*100);
end

%[text] **💡 観察ポイント 6 — 質量誤差の影響と MS vs FP アプローチを考察しましょう**
%[text] 上記ワークフローは RDKit からの ExactMolWt（単一同位体質量）を使用します。
%[text] 実際の実験では、測定された m/z には装置の不確実性があります。
%[text] データベース検索前に mH_obs に小さなランダム質量誤差を追加してみましょう:
%[text] mH_obs_noisy = mH_obs + randn(1, N_UNKNOWNS) * 0.002;  % 2 mDa ノイズ
%[text] 次にセクション 4〜6 を再実行して、同定精度の変化を確認しましょう。
%[text] 同位体コサイン類似度スコアは 0 から 1 の範囲です。
%[text] 誤った答えを返すのではなく「データベースにない」として
%[text] 同定を拒否するためにどの閾値を適用するとよいでしょうか？
%[text] （ヒント: 参照テーブルにない化合物でワークフローを実行し、
%[text] 偽候補のスコア分布を観察してみましょう。）
%[text] S05 は構造データベースから未知物質を同定するために Tanimoto
%[text] フィンガープリント類似度を使用しました。データベースに真にない化合物には、
%[text] Tanimoto が最も構造的に類似したアナログを見つけます。
%[text] 議論: MS アプローチ（A08）はいつ好ましく、FP アプローチ（S05）は
%[text] いつ好ましいでしょうか？ 両方を組み合わせるとどうなるでしょうか？
% ... （ここにコードを書いてみましょう）

%[text] **まとめ**
%[text]
%[text] - 精密単一同位体質量（ExactMolWt）が MS データベース検索の基礎です
%[text] - 高分解能 MS（<5 ppm）は候補を数十から 1〜2 に絞り込みます
%[text] - Cl・S を含む分子は M+2 ピークが異常に大きく、元素組成の「指紋」になります
%[text] - 質量 + 同位体コサイン類似度の 2 段階フィルタで同定精度が向上します
%[text]
%%
%[text] ## ローカル関数の概要

function spectrum = simulateSpectrum(peakMz, peakInt, mzRange, sigma, noiseLevel)
% SIMULATESPECTRUM  ガウスピークを持つ合成質量スペクトルを生成する。
% 
%   spectrum = simulateSpectrum(peakMz, peakInt, mzRange, sigma, noiseLevel)
% 
%   入力:
%     peakMz     -- ピーク中心 m/z 値の 1xN ベクトル（Da）
%     peakInt    -- 相対ピーク強度の 1xN ベクトル（0〜1）
%     mzRange    -- スペクトル軸の範囲 [mzMin, mzMax]
%     sigma      -- ガウス幅（Da）; FWHM = 2.355*sigma
%                   sigma = 0.3 Da --> FWHM = 0.71 Da（低分解能 ESI をシミュレート）
%     noiseLevel -- ベースピークの割合としての均一ノイズ振幅
% 
%   出力:
%     spectrum   -- フィールド .mz（1x5000）と .intensity（1x5000）を持つ struct
% 
%   強度範囲について:
%     ガウス成分はノイズ追加前に [0, 1] に正規化される。
%     均一ノイズ追加後（noiseLevel = 0.03）、最大強度は最大 1.03 に達する。
%     これは物理的にリアルです（ランダムノイズがピーク頂点の上に乗る）し
%     意図的なものです。厳密な [0, 1] 範囲が必要な呼び出し元は
%     後で max(spectrum.intensity) で割ること。
    mz        = linspace(mzRange(1), mzRange(2), 5000);
    intensity = zeros(1, 5000);
    for i = 1:numel(peakMz)
        intensity = intensity + peakInt(i) * ...
            exp(-0.5 * ((mz - peakMz(i)) / sigma).^2);
    end
    if max(intensity) > 0
        intensity = intensity / max(intensity);
    end
    intensity = intensity + noiseLevel * rand(1, 5000);
    spectrum  = struct("mz", mz, "intensity", intensity);
end

function counts = parseFormula(formulaStr)
% PARSEFORMULA  分子式文字列から元素数を抽出する。
% 
%   counts = parseFormula("C20H25N3O2")
%   フィールド C、H、N、O、S、Cl、F、Br、P、I（全て double）を持つ struct を返す。
%   存在しない元素は 0。
% 
%   正規表現: ([A-Z][a-z]?)(\d*)
%     [A-Z][a-z]?  元素記号にマッチ: 大文字 1 文字 + 任意の小文字 1 文字
%                  --> 2 文字記号を正しく捕捉: Cl、Br など
%                  --> "Cl" はトークン {"Cl",""} としてマッチ、{"C",""} + {"l",""} ではない
%     (\d*)        カウント数字にマッチ; 空文字列 --> 1 として扱う
% 
%   検証済みの例（RDKit rdMolDescriptors が返す Hill 表記）:
%     "CH4"       --> C:1, H:4
%     "C9H8O4"    --> C:9, H:8, O:4（アスピリン）
%     "C6H5Cl"    --> C:6, H:5, Cl:1（クロロベンゼン -- 2 文字元素 OK）
%     "CCl4"      --> C:1, Cl:4（四塩化炭素）
%     "H2O"       --> H:2, O:1（最後の元素に数字なし -- ガード OK）
% 
%   エッジケース:
%     空文字列 ""    --> 全カウント 0（tokens = {}; ループに入らない）
%     未知元素       --> isfield ガードで暗黙的に無視
%     （例: 無機不純物からの "Fe"、"Mg"）
    counts = struct("C",0,"H",0,"N",0,"O",0,"S",0,"Cl",0,"F",0,"Br",0,"P",0,"I",0);
    tokens = regexp(char(formulaStr), "([A-Z][a-z]?)(\d*)", "tokens");
    for i = 1:numel(tokens)
        elem = tokens{i}{1};
        n    = tokens{i}{2};
        if isempty(n), n = "1"; end
        if isfield(counts, elem)
            counts.(elem) = counts.(elem) + str2double(n);
        end
    end
end

function rel = isoPattern(fc)
% ISOPATTERN  M、M+1、M+2 の相対同位体強度を計算する。
% 
%   rel = isoPattern(fc)
%   fc  -- parseFormula() からの元素カウント struct
%   rel -- [1, M+1/M, M+2/M]（M = 1.0 に正規化）
% 
%   近似式（Gross 2011、第 3 章; 炭素原子 <100 個に有効）:
%     M+1% = 1.103*nC + 0.366*nN + 0.015*nH + 0.038*nO
%     M+2% = (1.103*nC)^2/200 + 0.205*nO + 4.25*nS + 32.7*nCl
% 
%   係数の導出（自然存在比: Meija et al. 2016 IUPAC）:
%     13C: 1.103% 自然存在比 --> M+1 寄与 = C 1 個あたり 1.103%
%     34S: 4.25%  自然存在比 --> M+2 寄与 = S 1 個あたり 4.25%
%     37Cl: 24.23% 自然存在比
%           M+2 寄与 = p(37Cl)/p(35Cl) * 100
%                    = 24.23/75.77 * 100 = 31.98% per Cl
%           注意: 31.98% は正確な値; 32.7% は Gross（2011）表 3.2 に従いここで使用。
%           差 < 2.2%。参考文献の「37Cl 24.23%」は生の自然存在比;
%           M+2 係数は存在比ではなく 37Cl/35Cl の比です。
% 
%   相互検証参照値（Claesen et al. 2012 で検証）:
%     CH4（C:1、H:4）:            rel ~ [1.000, 0.01163, 0.00006]
%       M+1% = 1.103 + 0.060 = 1.163%;  M+2% = 0.006%（ほとんど検出不可）
% 
%     C9H8O4 アスピリン（C:9、H:8、O:4）: rel ~ [1.000, 0.10199, 0.01313]
%       M+1% = 9.927 + 0.120 + 0.152 = 10.199%
%       M+2% = (9.927)^2/200 + 0.205*4 = 0.493 + 0.820 = 1.313%
% 
%     C6H5Cl クロロベンゼン（C:6、H:5、Cl:1）: rel ~ [1.000, 0.06693, 0.32919]
%       M+1% = 6.618 + 0.075 = 6.693%
%       M+2% = (6.618)^2/200 + 32.7 = 0.219 + 32.7 = 32.919%
%       --> M:M+2 比 = 1:0.329 が「Cl 3:1 ルール」を確認（Gross 2011 p.95）
    m1pct = 1.103*fc.C + 0.366*fc.N + 0.015*fc.H + 0.038*fc.O;
    m2pct = (1.103*fc.C)^2 / 200 + 0.205*fc.O + 4.25*fc.S + 32.7*fc.Cl;

    rel = [1.0, m1pct/100, m2pct/100];
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
