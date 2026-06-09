%[text] # A07: スキャフォールド分析と R 基分解
%[text] EasyMolKit アナリティクス — レイヤー 3
%[text]
%[text] FDA 承認薬を200種類並べると、多くの薬が同じ環系コアを共有していることに気づきます。
%[text] この「骨格」を **Bemis-Murcko スキャフォールド** と呼び、医薬化学者は設計段階でどのスキャフォールドが複数の薬を支えているかを把握するために活用します。
%[text] 同じ骨格でも側鎖（R 基）を変えるだけで吸収、毒性、活性が変わる場合、R 基を系統的に最適化できます。これが **構造-活性相関（SAR）** 解析の核心です。
%[text] このスクリプトでは、FDA 承認薬データセットから骨格を自動抽出し、頻出度ランキング、R 基特性テーブル、PCA 可視化を用いてスキャフォールド多様性を探ります。
%[text]
%[text] **ストーリー**
%[text]
%[text] ある医薬化学者が新しい化合物ライブラリを設計する前に、200 の FDA 承認薬の構造多様性を確認しています。彼女は3つの問いを立てました:
%[text]
%[text] - (a) 承認薬で最も多いコアフレームワーク（スキャフォールド）は何で、
%[text] ユニークなスキャフォールドはいくつあるか？
%[text] - (b) スキャフォールドファミリー間でどの物理化学的特性（LogP: 水/オクタノール間の分配係数の常用対数、TPSA: 位相的極性表面積、MW: 分子量）が
%[text] 異なるか？これが SAR 分析の本質であり、構造変化を特性変化にリンクします。
%[text] - (c) 同じ側鎖（R 基）は、特定のスキャフォールドファミリー内で
%[text] 薬物様性にどう影響するか？
%[text]
%[text] この演習では以下を行います:
%[text]
%[text] 1. 200 の FDA 薬を読み込み、各薬物の Bemis-Murcko スキャフォールドを抽出します。
%[text] 2. ユニークなスキャフォールド SMILES をマッピングしてスキャフォールドファミリーを特定します。
%[text] 3. スキャフォールドを頻度でランク付けし、分布をプロットします。
%[text] 4. R 基特性テーブルを構築します: ファミリーごとの主要記述子の平均と標準偏差。
%[text] 5. SAR を可視化します: スキャフォールドクラスで色分けした PCA と ALogP ボックスプロット。
%[text] 6. スキャフォールド多様性指数（SDI）でスキャフォールド多様性を要約します。
%[text]
%[text] **学習目標**
%[text]
%[text] - Bemis-Murcko スキャフォールドの定義と SAR での使用を理解します
%[text] - `emk.mol.scaffold()` を使って SMILES から正準スキャフォールドを抽出します
%[text] - MATLAB `containers.Map` でスキャフォールド SMILES ごとに分子をグループ化します
%[text] - スキャフォールド内分散を要約する R 基特性テーブルを構築します
%[text] - スキャフォールドが化学空間をどう分割するかを可視化します（PCA、ボックスプロット）
%[text] - スキャフォールド多様性指数（SDI）を計算し、解釈します
%[text]
%[text] **前提条件**
%[text]
%[text] - F05（部分構造検索）修了 ―― スキャフォールドの概念が導入済みです
%[text] - 推奨: 可視化の文脈のために A01（PCA）と A02（クラスタリング）を学習済みであること
%[text] - Statistics and Machine Learning Toolbox（`pca`、`boxplot`）が必要です
%[text] - インターネット接続は不要です
%[text]
%[text] 推定所要時間: 30〜45 分
%[text]
%[text] **データ:**
%[text]
%[text] `data/list/fda_drugs.csv` ―― 200 FDA 承認薬（ChEMBL CC-BY-SA 3.0）
%[text] 列: ChEMBLID、Name、SMILES、MolecularWeight、ALogP、
%[text] HBondDonors、HBondAcceptors、TPSA、RotatableBonds
%[text]
%[text] **参考文献**
%[text]
%[text] Bemis & Murcko (1996) *J. Med. Chem.* 39:2887–2893. doi:10.1021/jm9602928
%[text] ―― スキャフォールド概念の原著論文
%[text]
%[text] Bemis & Murcko (1999) *J. Med. Chem.* 42:5095–5099. doi:10.1021/jm9903996
%[text] ―― R 基 / 側鎖分析への拡張
%[text]
%[text] Langdon et al. (2011) *J. Chem. Inf. Model.* 51:2174–2185. doi:10.1021/ci200319g
%[text] ―― FDA 薬物セットの定量的スキャフォールド多様性分析
%[text]
%[text] RDKit `MurckoScaffold` モジュール:
%[text] https://www.rdkit.org/docs/source/rdkit.Chem.Scaffolds.MurckoScaffold.html
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

%[text] メインの実行前に、Python/RDKit プロセスをウォームアップしておきます。
mol_warmup = emk.mol.fromSmiles("C");   % メタン -- 軽量
clear mol_warmup;
logSection("A07", "セクション 0: セットアップ", "アナリティクス L3");
%%
%[text] ## セクション 1: FDA 承認薬の読み込みと Murcko スキャフォールド抽出
%[text]
%[text] セットアップが完了しました。まず、FDA 承認薬データを読み込み、各分子から Bemis-Murcko スキャフォールドを抽出します。
%[text] 正しく抽出できると、同じ環系コアを持つ薬物が自動的にグループ化されます。
%[text]
%[text] ### コンセプト: Bemis-Murcko スキャフォールド
%[text] スキャフォールド（分子フレームワーク）は、末端置換基（側鎖）をすべて除去した後の分子のコア環系です。
%[text]
%[text] **Bemis-Murcko 定義（1996 年）** が保持するのは以下の要素です:
%[text] - すべての環原子（N、O、S などのヘテロ原子を含む）
%[text] - 環を接続するリンカー結合と原子
%[text] - 環間の一原子リンカー
%[text]
%[text] 環系に付いた側鎖は除去されます。例えば:
%[text] - アスピリン `CC(=O)Oc1ccccc1C(=O)O` → スキャフォールド `c1ccccc1`（ベンゼン）
%[text] - イブプロフェン `CC(C)Cc1ccc(C(C)C(=O)O)cc1` → スキャフォールド `c1ccc(cc1)`（ベンゼン）
%[text]
%[text] 両方ともベンゼンスキャフォールドを持ちますが、R 基は大きく異なります。
%[text] スキャフォールド分析により、どのフレームワークが複数の薬物を担っているかがわかります。
%[text]
%[text] 環を持たない分子（例: アミノ酸、短いペプチド）はスキャフォールドがありません。
%[text] RDKit はこれらに対して空の Mol（0 原子）を返します。
%[text] データセットに保持するために `"<acyclic>"` とラベル付けします。
logSection("A07", "セクション 1: FDA 承認薬の読み込みと Murcko スキャフォールド抜出", "アナリティクス L3");
DATA_FILE   = "data/list/fda_drugs.csv";
ACYCLIC_TAG = "<acyclic>";   % placeholder for ring-free molecules

rawTbl = readtable(DATA_FILE, TextType="string");
nRaw   = height(rawTbl);
logInfo("%s から %d 行を読み込みました", DATA_FILE, nRaw);

%[text] 各薬物のスキャフォールドを抽出
scaffoldSmiles = strings(nRaw, 1);   % 薬物ごとのスキャフォールド SMILES
scaffoldNAtoms = nan(nRaw, 1);       % スキャフォールドの重原子数
mols           = cell(1, nRaw);
valid          = false(1, nRaw);

for k = 1:nRaw
    try
        mol        = emk.mol.fromSmiles(rawTbl.SMILES(k));
        scaf       = emk.mol.scaffold(mol);
        nScafAtoms = double(scaf.GetNumAtoms());
        if nScafAtoms == 0
            scaffoldSmiles(k) = ACYCLIC_TAG;
        else
            scaffoldSmiles(k) = emk.mol.toSmiles(scaf);
        end
        scaffoldNAtoms(k) = nScafAtoms;
        mols{k}           = mol;
        valid(k)          = true;
    catch ME
        logWarn("行 %d (%s) をスキップ: %s", k, rawTbl.Name(k), ME.message);
    end
end

validIdx = find(valid);
nMols    = numel(validIdx);
logInfo("%d / %d 分子のスキャフォールドを抽出しました", nMols, nRaw);

%[text] 有効な分子に作業配列を制限
molNames    = rawTbl.Name(validIdx);
molSmiles   = rawTbl.SMILES(validIdx);
scafSmi     = scaffoldSmiles(validIdx);
scafNAtoms  = scaffoldNAtoms(validIdx);

%[text] プレビュー（最初の 5 行）
previewTbl = table( ...
    molNames(1:5), molSmiles(1:5), scafSmi(1:5), scafNAtoms(1:5), ...
    VariableNames=["Name", "SMILES", "Scaffold", "ScaffoldAtoms"]);
disp(previewTbl);

%[text] **💡 観察ポイント 1 ―― スキャフォールドの抽出結果を確認しましょう**
%[text] 200 薬のうち、非環式（環系なし）のものは何個ありますか？
%[text] `sum(scafSmi == "<acyclic>")` を実行して確認しましょう。
%[text] 非環式薬物名に見覚えはありますか？（ヒント: 単純なアミノ酸、一部の抗ウイルス薬）
%[text] 最大のスキャフォールド（重原子が最多）はどれですか？
%[text] スキャフォールド SMILES は複雑な環系に見えますか？
%[text] アスピリンの SMILES は `"CC(=O)Oc1ccccc1C(=O)O"` です。
%[text] どのスキャフォールドが得られますか？次にイブプロフェン `"CC(C)Cc1ccc(C(C)C(=O)O)cc1"` でも試してみましょう。
%[text] 同じスキャフォールドを共有していますか？
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 2: スキャフォールド頻度分析
%[text]
%[text] スキャフォールドの抽出が完了しました。次に、各スキャフォールドがどの薬物にどれだけ出現するか、頻度分布を調べます。
%[text] 特権的構造がどの程度少数に集中しているかを確認しましょう。
%[text]
%[text] ### コンセプト: スキャフォールド多様性とスキャフォールド頻度
logSection("A07", "セクション 2: スキャフォールド頻度分析", "アナリティクス L3");
%[text] スキャフォールド頻度の分布を確認することで、薬物ライブラリが「多様」（多くのユニークなスキャフォールドがそれぞれ1回だけ登場）か「集中」（少数のスキャフォールドを多くの化合物が共有）かを判断できます。
%[text]
%[text] FDA 薬物セット（Langdon et al. 2011）の傾向:
%[text] - スキャフォールドの約 50% が「シングルトン」（1つの薬物にのみ出現）
%[text] - 少数のスキャフォールド（ベンゼン、ピリジン、ピペリジン環系）が
%[text] 5〜10種以上の薬物に出現します
%[text]
%[text] **なぜ重要でしょうか？**
%[text] - 高スキャフォールド多様性 = 化学空間の広いカバレッジを示します
%[text] - 高頻度スキャフォールドは「**特権的構造**」であり、
%[text] 多くのタンパク質結合部位に適合する進化的に選択されたフレームワークです
%[text]
%[text] **スキャフォールド多様性指数（SDI）** の定義:
%[text] $\text{SDI} = N_{\text{unique}} / N_{\text{total}}$
%[text] SDI の範囲は $1/N$（全て同じスキャフォールド）から $1.0$（全てユニーク）です。
%[text] ランダムな薬物様ライブラリは通常 SDI ~ 0.7〜0.9 を達成します。
%[text]
%[text] スキャフォールド SMILES を分子インデックスのリストにマップします（containers.Map を使用）。
scafMap = containers.Map("KeyType", "char", "ValueType", "any");
for k = 1:nMols
    key = char(scafSmi(k));
    if isKey(scafMap, key)
        scafMap(key) = [scafMap(key), k];
    else
        scafMap(key) = k;
    end
end

uniqueScafs = keys(scafMap);
nUnique     = numel(uniqueScafs);
SDI         = nUnique / nMols;

logInfo("ユニークなスキャフォールド: %d / %d 分子  (SDI = %.3f)", ...
    nUnique, nMols, SDI);

%[text] 各スキャフォールドの頻度を計算します。
scafFreq = zeros(1, nUnique);
for i = 1:nUnique
    scafFreq(i) = numel(scafMap(uniqueScafs{i}));
end

%[text] 頻度を降順でソートします。
[scafFreqSorted, sortOrd] = sort(scafFreq, "descend");
uniqueScafsSorted = uniqueScafs(sortOrd);

%[text] 上位 10 スキャフォールドを報告します。
TOP_N = 10;
logInfo("頻度順上位 %d スキャフォールド:", TOP_N);
for i = 1:min(TOP_N, nUnique)
    memberIdx  = scafMap(uniqueScafsSorted{i});
    memberList = strjoin(molNames(memberIdx), ", ");
    if strlength(memberList) > 80
        memberList = extractBefore(memberList, 81) + "...";
    end
    logInfo("  [%2d] 頻度=%d  スキャフォールド=%s  メンバー=%s", ...
        i, scafFreqSorted(i), uniqueScafsSorted{i}, memberList);
end

nSingletons = sum(scafFreq == 1);
logInfo("シングルトン（頻度=1）: %d / %d スキャフォールド（%.0f%%）", ...
    nSingletons, nUnique, 100 * nSingletons / nUnique);

%[text] -- 棒グラフ: 上位 10 スキャフォールドの頻度 --
figure("Name", "A07 スキャフォールド頻度");
barh(scafFreqSorted(1:TOP_N), FaceColor=[0.3 0.6 0.9]);
yticks(1:TOP_N);
%[text] スキャフォールド SMILES の最初の 30 文字をラベルとして使用します。
shortLabels = cell(TOP_N, 1);
for i = 1:TOP_N
    s = uniqueScafsSorted{i};
    if numel(s) > 30
        shortLabels{i} = [s(1:30) "..."];
    else
        shortLabels{i} = s;
    end
end
yticklabels(shortLabels);   % YDir="reverse" で y=1 が上端になるため flipud 不要
xlabel("同じスキャフォールドを持つ FDA 薬の数");
title(sprintf("最も頻出する Murcko スキャフォールド上位 %d（FDA 薬）", TOP_N));
grid on;
set(gca, FontSize=8, YDir="reverse");

%[text] -- ヒストグラム: 全スキャフォールド頻度の分布 --
figure("Name", "A07 スキャフォールド頻度分布");
histogram(scafFreq, max(scafFreq), FaceColor=[0.3 0.6 0.9], EdgeColor="none");
xlabel("スキャフォールドあたりの薬物数");
ylabel("スキャフォールド数");
title(sprintf("スキャフォールド頻度分布（ユニーク N=%d）", nUnique));
xline(mean(scafFreq), "--r", sprintf("平均=%.1f", mean(scafFreq)), ...
    LabelHorizontalAlignment="right");
grid on;

%[text] **💡 観察ポイント 2 ―― 頻度分布を分析しましょう**
%[text] 3つ以上の薬物に出現するスキャフォールドの割合を確認しましょう。
%[text] `sum(scafFreq >= 3) / nUnique * 100` を用いて計算できます。
%[text] これらが「特権的」スキャフォールドと呼ばれる理由を考えてみましょう。
%[text] スキャフォールド頻度分布のシャノンエントロピーを計算してみましょう。
%[text] H が高いほど、スキャフォールド間でより均等な分布を意味します。
%[text] 仮想的な「全シングルトン」ライブラリと比較すると、`H = log2(nMols)` です。
%[text] 多くの FDA 分析で最も頻出するスキャフォールドは単純なベンゼン環（`c1ccccc1`）です。
%[text] このデータセットでも同様か確認しましょう。`uniqueScafsSorted{1}` をチェックします。
%[text] ベンゼンがこれほど一般的なのはなぜでしょうか？
%[text] （ヒント: 芳香環は疎水性表面、パイスタッキング、代謝安定性を提供します）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 3: R 基特性テーブル
%[text]
%[text] 頼度分布を把握しました。次に、各スキャフォールドファミリー内で R 基（側鎖）が物性にどのように影響するかを定量化します。
%[text] 同じ骨格でも ALogP（親油性の指標）や TPSA（位相的極性表面積）がどれほど変わるかが見えてきます。
%[text]
%[text] ### コンセプト: R 基分析と構造-活性相関（SAR）
%[text] 医薬化学では、「R 基」は固定されたスキャフォールドに付加された可変置換基を指します。R 基分析では次の問いを立てます。
%[text]
%[text] 固定されたスキャフォールドで、異なる側鎖が親油性（ALogP）、極性（TPSA）、分子量などの特性にどのように影響するでしょうか？
%[text]
%[text] これはリード最適化の核心です。創薬研究者は活性を維持しつつ、ADMET（吸収・分布・代謝・排泄・毒性）特性を改善するために有望なスキャフォールドの R 基を修飾します。
%[text]
%[text] このセクションでは、3 つ以上のメンバーを持つスキャフォールドファミリーごとに特性サマリーテーブルを構築します。各ファミリーについて報告する内容は以下の通りです。
%[text] - メンバー数（同じスキャフォールドを共有する承認薬の数）
%[text] - ALogP、TPSA、MolecularWeight、HBondDonors、HBondAcceptors、
%[text] RotatableBonds の平均と標準偏差
%[text]
%[text] スキャフォールド内の高い `std(ALogP)` は、スキャフォールドが R 基変化により広い親油性範囲を許容することを意味します。これは「柔軟な」スキャフォールドです。
%[text] 低い std は、そのスキャフォールド上のすべての薬物が同様の親油性を持つことを示します。
logSection("A07", "セクション 3: R 基特性テーブル", "アナリティクス L3");
MIN_FAMILY_SIZE = 3;   % minimum members to include in the R-group table

%[text] サイズ閾値を満たすスキャフォールドファミリーを選択します。
%[text] 注: `<acyclic>` プレースホルダーは除外します。環系を持たない分子には固定スキャフォールドが
%[text] 存在しないため、R 基分析の対象外となります（環系薬物のみが SAR で比較可能です）。
realScafMask = scafFreqSorted >= MIN_FAMILY_SIZE & ~strcmp(uniqueScafsSorted, char(ACYCLIC_TAG));
familyScafs  = uniqueScafsSorted(realScafMask);
familyFreqs  = scafFreqSorted(realScafMask);
nFamilies    = numel(familyScafs);
logInfo("メンバー >= %d のスキャフォールドファミリー（環系のみ）: %d", MIN_FAMILY_SIZE, nFamilies);

%[text] CSV に存在する特性列を使用します（RDKit のオーバーヘッドを避けるために直接使用）。
PROP_COLS = ["ALogP", "TPSA", "MolecularWeight", "HBondDonors", ...
             "HBondAcceptors", "RotatableBonds"];

%[text] R 基サマリーテーブルを構築します。
familyNames     = cell(nFamilies, 1);
familyNMembers  = zeros(nFamilies, 1);
familyMeanLogP  = zeros(nFamilies, 1);
familyStdLogP   = zeros(nFamilies, 1);
familyMeanTPSA  = zeros(nFamilies, 1);
familyMeanMW    = zeros(nFamilies, 1);
familyMeanHBD   = zeros(nFamilies, 1);
familyMeanHBA   = zeros(nFamilies, 1);

for f = 1:nFamilies
    memberIdx = scafMap(familyScafs{f});   % indices into validIdx
    % これらのメンバーの rawTbl から行を取得
    rawRows = rawTbl(validIdx(memberIdx), :);

    logP = double(rawRows.ALogP);
    tpsa = double(rawRows.TPSA);
    mw   = double(rawRows.MolecularWeight);
    hbd  = double(rawRows.HBondDonors);
    hba  = double(rawRows.HBondAcceptors);

    familyNames{f}    = familyScafs{f};
    familyNMembers(f) = numel(memberIdx);
    familyMeanLogP(f) = mean(logP);
    familyStdLogP(f)  = std(logP);
    familyMeanTPSA(f) = mean(tpsa);
    familyMeanMW(f)   = mean(mw);
    familyMeanHBD(f)  = mean(hbd);
    familyMeanHBA(f)  = mean(hba);

    % 最初の 5 ファミリーのメンバー名を表示
    if f <= 5
        memberStr = strjoin(rawRows.Name, " | ");
        if strlength(memberStr) > 100
            memberStr = extractBefore(memberStr, 101) + "...";
        end
        logInfo("ファミリー %d (n=%d): %s", f, familyNMembers(f), memberStr);
    end
end

%[text] サマリーテーブルを組み立てます。
rgroupTbl = table( ...
    string(familyNames), familyNMembers, ...
    round(familyMeanLogP, 2), round(familyStdLogP, 2), ...
    round(familyMeanTPSA, 1), round(familyMeanMW, 1), ...
    round(familyMeanHBD, 2), round(familyMeanHBA, 2), ...
    VariableNames=["Scaffold","N","mean_ALogP","std_ALogP", ...
                   "mean_TPSA","mean_MW","mean_HBD","mean_HBA"]);
rgroupTbl = sortrows(rgroupTbl, "N", "descend");

logInfo("R 基特性テーブル（%d スキャフォールドファミリー、>= %d メンバー）:", ...
    nFamilies, MIN_FAMILY_SIZE);
disp(rgroupTbl);

%[text] **💡 観察ポイント 3 ―― R 基テーブルを読み解きましょう**
%[text] 最も広い ALogP 範囲（最高の `std_ALogP`）を持つスキャフォールドファミリーを確認しましょう。
%[text] スキャフォールドの SMILES はよくある芳香環でしょうか？
%[text] 高い LogP 分散は R 基の柔軟性について何を示していますか？
%[text] 最高の平均 TPSA を持つスキャフォールドを確認し、それが最も極性が高いかどうかを読み取りましょう。
%[text] 高 TPSA（$> 140\,\text{Å}^2$）は経口吸収不良と関連します（Veber ルール）。
%[text] 通常 IV 投与される薬物クラスを挙げることができますか？
%[text] 最も頻出するスキャフォールドの全メンバー薬物名を表示してみましょう。
%[text] これらの薬物は同じ薬理クラスに属していますか？
%[text] n=1 メンバー（シングルトン）のスキャフォールドはファミリー内変動がありません。
%[text] これは SAR に何を意味するでしょうか？1 つのデータポイントで R 基最適化は可能でしょうか？
%[text] （ヒント: SAR は複数のアナログが必要です）
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 4: SAR 可視化 — スキャフォールドファミリー間の特性変動
%[text]
%[text] R 基特性テーブルが完成しました。次に、ボックスプロットと主成分分析（PCA）を用いて、構造-活性相関（SAR）を視覚的に確認します。
%[text] スキャフォールドファミリーが化学空間をどのように分割しているかを確認しましょう。
%[text]
%[text] ### コンセプト: 構造-活性相関（SAR）可視化
logSection("A07", "セクション 4: SAR 可視化 — スキャフォールドファミリー間の特性変動", "アナリティクス L3");
%[text] SAR プロットは、異なるスキャフォールドファミリーによる構造変化が物理化学的特性にどのように影響するかを示します。
%[text]
%[text] ここでは、2 つの補完的な視点を用います。
%[text]
%[text] **(a) スキャフォールドファミリーごとの ALogP ボックスプロット**
%[text] 各ボックスは、1 つのファミリー内の ALogP の中央値と四分位範囲を示します。
%[text] 広いボックス（大きな四分位範囲）は、親油性に影響する多様な R 基を示します。
%[text]
%[text] **(b) スキャフォールドクラスで色分けした記述子行列の PCA**
%[text] 全薬物を 2D 化学空間に射影します。スキャフォールドファミリーが明確なクラスターを形成する場合、スキャフォールドの同一性が全体的な物理化学的プロファイルを予測することを示します。
%[text] これが「スキャフォールド駆動」SAR の特徴です。
%[text] 重複するファミリーは、R 基が特性を支配していることを示します。
%[text]
%[text] スキャフォールドによる PCA 色付け:
%[text] - 各「上位スキャフォールドファミリー」にユニークな色を割り当てます。
%[text] - シングルトンと小さなファミリーはグレーでプロットします。
%[text] - PCA 空間での分離は、スキャフォールドが全体的な特性を駆動している証拠です。
%[text]
%[text] ボックスプロットには、スキャフォールドファミリーが MIN_FAMILY_SIZE 以上の薬物のみを使用します。
familyMask   = ismember(scafSmi, string(familyScafs));
nFamilyMols  = sum(familyMask);
logInfo("スキャフォールドファミリー内の分子（>= %d メンバー）: %d / %d", ...
    MIN_FAMILY_SIZE, nFamilyMols, nMols);

if nFamilyMols > 0
    % グループ化ボックスプロットのラベルベクトルを構築
    familyScafStr = string(familyScafs);
    [~, familyLabelIdx] = ismember(scafSmi(familyMask), familyScafStr);

    % 表示のためにスキャフォールドラベルを短縮（12 文字）
    shortScafLabels = cell(nFamilies, 1);
    for f = 1:nFamilies
        s = familyScafs{f};
        shortScafLabels{f} = s(1:min(12, numel(s)));
    end

    % ファミリー分子の ALogP を取得
    logP_family = double(rawTbl.ALogP(validIdx(familyMask)));

    % 順序表示のために平均 ALogP でファミリーをソート
    [~, sortByMean] = sort(familyMeanLogP, "descend");
    remappedLabels  = zeros(size(familyLabelIdx));
    for fi = 1:nFamilies
        remappedLabels(familyLabelIdx == sortByMean(fi)) = fi;
    end

    % 注意: LabelOrientation="inline" は MATLAB Online でグラフィックスタイムアウトを引き起こす
    % （renderLabels が多くのラベルで固まる）。数値グループで boxplot を使い
    % XTickLabel を別途設定する -- より高速なレンダリング。
    figure("Name", "A07 スキャフォールドファミリーごとの ALogP");
    boxplot(logP_family, remappedLabels, Labels=shortScafLabels(sortByMean));
    xtickangle(30);
    xlabel("スキャフォールド（SMILES 最初の 12 文字、平均 ALogP でソート）");
    ylabel("ALogP");
    title(sprintf("スキャフォールドファミリーごとの ALogP 分布（ファミリー >= %d メンバー）", ...
        MIN_FAMILY_SIZE));
    yline(0, "--", "LogP=0", Color=[0.5 0.5 0.5]);
    yline(5, "--", "Ro5 上限", Color=[0.8 0.2 0.2]);
    grid on;
else
    logWarn(">= %d メンバーのスキャフォールドファミリーが見つからない -- ボックスプロットをスキップ。", ...
        MIN_FAMILY_SIZE);
end

%[text] --- スキャフォールドファミリーで色分けした記述子行列の PCA ---
%[text] CSV 提供の記述子を使用し、速度向上のために分子ごとの RDKit 呼び出しを避けます。
FEAT_COLS = ["MolecularWeight", "ALogP", "TPSA", "HBondDonors", ...
             "HBondAcceptors", "RotatableBonds"];
nFeats    = numel(FEAT_COLS);

%[text] 全有効分子の特徴行列を構築します。
X_pca = zeros(nMols, nFeats);
for fi = 1:nFeats
    X_pca(:, fi) = double(rawTbl.(FEAT_COLS(fi))(validIdx));
end

%[text] 標準化（ゼロ平均、単位分散）を行います。
mu_pca    = mean(X_pca, 1);
sigma_pca = std(X_pca, 0, 1);
sigma_pca(sigma_pca == 0) = 1;   % ゼロ分散列に対するガード
X_std = (X_pca - mu_pca) ./ sigma_pca;

%[text] PCA を実行します。
[~, scores, ~, ~, explained] = pca(X_std);

logInfo("PCA: PC1=%.1f%%、PC2=%.1f%%、累積=%.1f%%", ...
    explained(1), explained(2), sum(explained(1:2)));

%[text] スキャフォールドファミリーメンバーシップに基づいて色を割り当てます。
%[text] ファミリーは頻度でランク付けし、上位 6 つには明確な色を、残りはグレーにします。
MAX_COLOUR_FAMILIES = 6;
palette = lines(MAX_COLOUR_FAMILIES);   % MATLAB colourmap with distinct colours
colourIdx = zeros(nMols, 1);            % 0 = grey (singleton / small family)

for f = 1:min(MAX_COLOUR_FAMILIES, nFamilies)
    memberIdx = scafMap(familyScafs{f});
    colourIdx(memberIdx) = f;
end

figure("Name", "A07 スキャフォールドによる化学空間 PCA");
hold on;

%[text] 最初にグレーの「その他」分子をプロットし、背景レイヤーを作成します。
greyMask = colourIdx == 0;
scatter(scores(greyMask, 1), scores(greyMask, 2), 30, ...
    [0.75 0.75 0.75], "filled", MarkerFaceAlpha=0.4, DisplayName="その他");

%[text] 次に、色付きスキャフォールドファミリーをプロットします。
legendHandles = gobjects(1, min(MAX_COLOUR_FAMILIES, nFamilies) + 1);
legendHandles(1) = scatter(scores(greyMask, 1), scores(greyMask, 2), 30, ...
    [0.75 0.75 0.75], "filled", MarkerFaceAlpha=0.4);  % その他

for f = 1:min(MAX_COLOUR_FAMILIES, nFamilies)
    memberIdx = scafMap(familyScafs{f});
    col = palette(f, :);
    s = familyScafs{f};
    labelStr = sprintf("Scaf %d (n=%d): %s", f, familyFreqs(f), ...
        s(1:min(20, numel(s))));
    legendHandles(f + 1) = scatter( ...
        scores(memberIdx, 1), scores(memberIdx, 2), ...
        50, col, "filled", MarkerFaceAlpha=0.85, DisplayName=labelStr);
end

hold off;
xlabel(sprintf("PC1 (%.1f%%)", explained(1)));
ylabel(sprintf("PC2 (%.1f%%)", explained(2)));
title("化学空間 PCA -- スキャフォールドファミリーで色分け");
legend(legendHandles, ["その他（シングルトン / 小ファミリー）", ...
    arrayfun(@(f) sprintf("スキャフォールド %d (n=%d)", f, familyFreqs(f)), ...
    1:min(MAX_COLOUR_FAMILIES, nFamilies), UniformOutput=false)], ...
    Location="best", FontSize=7);
grid on;

%[text] **💡 観察ポイント 4 ―― SAR を可視化してみましょう**
%[text] 同じスキャフォールドファミリーの分子は PCA 空間でクラスターを形成していますか？
%[text] 形成する場合、スキャフォールドが物理化学的プロファイルを強く決定しています。
%[text] 形成しない場合、R 基が支配的であることが多く、柔軟なスキャフォールドに見られます。
%[text] PC1 と PC2 は合計分散のどれくらいを説明していますか？
%[text] 60% 未満の場合、2D 射影が重要な構造を見逃している可能性があります。
%[text] 試してみましょう: `sum(explained(1:3))` で 3D PCA がどれくらいを捉えるか確認します。
%[text] 最大スキャフォールドファミリー内のペアワイズ Tanimoto 類似度を計算し、
%[text] データセット全体のランダムペアと比較してみましょう。
%[text] `mean(simMat_family(~logical(eye(numel(memberIdx)))))`
%[text] ファミリー内 Tanimoto 類似度は A02 のグローバル平均より高いですか？
% ... （ここにコードを書いてみましょう）
%%
%[text] ## セクション 5: スキャフォールド多様性サマリー
%[text]
%[text] 可視化が完了しました。最後に、スキャフォールド多様性指数（SDI）を用いてこのライブラリの多様性を数値でまとめます。
%[text] FDA 薬物データセットは広範な化学空間をカバーしているでしょうか？
%[text]
%[text] ### コンセプト: スキャフォールド多様性指数とライブラリ設計
logSection("A07", "セクション 5: スキャフォールド多様性サマリー", "アナリティクス L3");
%[text] スキャフォールド多様性指数（SDI）は、化合物ライブラリが化学空間に
%[text] どれだけ「広がっている」かを示す単一の数値です:
%[text]
%[text] **SDI = ユニークなスキャフォールド数 / 分子の総数**
%[text]
%[text] **解釈**:
%[text] - SDI が 1.0 に近い場合 → 各分子がユニークなスキャフォールドを持つ（最大多様性）
%[text] - SDI が 0.0 に近い場合 → 全分子が 1 つのスキャフォールドを共有（集中ライブラリ）
%[text]
%[text] **関連指標**:
%[text] - シングルトン割合: 1 回だけ出現するスキャフォールドの割合
%[text] （高い割合 = 多様なセット、アナログなしでは構造-活性相関（SAR）が難しい）
%[text] - 平均ファミリーサイズ: 全分子 / ユニークスキャフォールド（SDI の逆数）
%[text] - 最大ファミリーサイズ: 単一の特権的スキャフォールドの支配性
%[text]
%[text] **ライブラリ設計のヒューリスティクス**:
%[text] - 多様なスクリーニングライブラリ: SDI > 0.7、シングルトン割合 > 50%
%[text] - 集中 SAR ライブラリ: 目標 SDI < 0.5、5〜20 アナログのファミリー
%[text] - FDA 承認薬: SDI 通常 0.7〜0.8（Langdon 2011）
logInfo("=== スキャフォールド多様性サマリー ===");
logInfo("分析分子数合計          : %d", nMols);
logInfo("ユニークスキャフォールド : %d", nUnique);
logInfo("スキャフォールド多様性指数: %.3f", SDI);
logInfo("シングルトンスキャフォールド: %d（%.0f%%）", ...
    nSingletons, 100 * nSingletons / nUnique);
logInfo("平均ファミリーサイズ     : %.2f 分子 / スキャフォールド", nMols / nUnique);
logInfo("最大スキャフォールドファミリー: %d 分子 (%s)", ...
    scafFreqSorted(1), uniqueScafsSorted{1});
logInfo("非環式分子               : %d（環系なし）", ...
    sum(scafSmi == ACYCLIC_TAG));

%[text] スキャフォールドサイズ分布（スキャフォールドの重原子数）
scafAtomsNonAcyclic = scafNAtoms(scafNAtoms > 0);
logInfo("スキャフォールド重原子数 : 平均=%.1f、標準偏差=%.1f、範囲=[%d,%d]", ...
    mean(scafAtomsNonAcyclic), std(scafAtomsNonAcyclic), ...
    min(scafAtomsNonAcyclic), max(scafAtomsNonAcyclic));

figure("Name", "A07 スキャフォールドサイズ分布");
histogram(scafAtomsNonAcyclic, "BinWidth", 2, ...
    FaceColor=[0.3 0.6 0.9], EdgeColor="none");
xlabel("スキャフォールド重原子数");
ylabel("薬物数");
title("スキャフォールドサイズの分布（Murcko スキャフォールド）");
xline(mean(scafAtomsNonAcyclic), "--r", ...
    sprintf("平均=%.1f", mean(scafAtomsNonAcyclic)), ...
    LabelHorizontalAlignment="right");
grid on;

%[text] **まとめ**
%[text]
%[text] - Bemis-Murcko スキャフォールドは環系コアを抽出し、薬物ファミリーを自動でグループ化します
%[text] - FDA 承認薬では少数の「特権的スキャフォールド」が多くの薬物に共有されています
%[text] - R 基分析により同一スキャフォールド内の構造-活性相関（SAR）が明確になります
%[text] - スキャフォールド多様性指数（SDI）を用いてライブラリの化学空間カバレッジを定量評価できます
%[text]

%[text] **💡 観察ポイント 5 ―― スキャフォールド多様性を評価してみましょう**
%[text] FDA 薬物セットの SDI は、ランダム多様性ライブラリに期待される値と
%[text] どのように比較されますか？（Langdon 2011 は 836 FDA 承認薬に SDI ~ 0.75 を報告しています）
%[text] ヒント: この演習では 200 薬のサブセットを使用しているため、SDI がやや低くなることがあります。
%[text] サブセットサイズと SDI の関係について何が言えますか？
%[text] 「スキャフォールド豊富さ」曲線を計算してみましょう:
%[text] 最初の 10、20、... 200 薬を追加するにつれてユニークなスキャフォールドがいくつ現れるかを調べます。
%[text] `figure; plot(1:nMols, richness, "-");`
%[text] `xlabel("追加した分子数"); ylabel("ユニークスキャフォールド"); title("スキャフォールド豊富度曲線");`
%[text] 曲線がどれくらい速く平坦になりますか？
%[text] `everyday_chemicals.csv` データセットを使ってスキャフォールド分析をやり直してみましょう。
%[text] より小さく厳選されたセットの SDI は高いですか、低いですか？
%[text] データセットサイズが多様性指標にどのように影響するかについて何を伝えていますか？
% ... （ここにコードを書いてみましょう）

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
