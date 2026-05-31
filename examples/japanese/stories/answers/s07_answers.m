%[text] # S07 解答: ChEMBL 生物活性データ
%[text] s07_chembl_activity.m の「やってみよう」演習の参考解答。
%[text] まず s07_chembl_activity.m を実行してから、このファイルで答え合わせをすること。
%[text] 注意: このファイルはインターネット接続が必要（ChEMBL REST API）。
addpath(genpath("src"));
emk.setup.initPython();
set(0, "DefaultFigureWindowStyle", "normal");   % 図をポップアップウィンドウで表示
logInfo("S07 解答: セットアップ完了（インターネット接続必要）");

%[text] ---- s07 の共有状態を再現する ----------------------------------------------
ERLOTINIB_SMILES = "C#Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1";
ERLOTINIB_NAME   = "Erlotinib";
EGFR_CHEMBL_ID   = "CHEMBL203";
IC50_CUTOFF_NM   = 100;

actTbl    = emk.db.getChemblActivity(EGFR_CHEMBL_ID, ActivityType="IC50", MaxRows=50);
actTbl.pIC50 = 9 - log10(actTbl.Value_nM);
activeTbl = sortrows(actTbl(actTbl.Value_nM <= IC50_CUTOFF_NM, :), "Value_nM", "ascend");

%[text] 活性化合物の SMILES を変換する
validMols = {};
validIdx  = [];
for i = 1:height(activeTbl)
    smi = activeTbl.SMILES(i);
    if ~emk.mol.isValid(smi); continue; end
    validMols{end+1} = emk.mol.fromSmiles(smi); %#ok<AGROW>
    validIdx(end+1)  = i;                        %#ok<AGROW>
end

mol_erl = emk.mol.fromSmiles(ERLOTINIB_SMILES);
fp_erl  = emk.fingerprint.morgan(mol_erl, Radius=2, NBits=2048);

dbFps = cell(1, numel(validMols));
for i = 1:numel(validMols)
    dbFps{i} = emk.fingerprint.morgan(validMols{i}, Radius=2, NBits=2048);
end
rankResult = emk.similarity.rankBy(fp_erl, dbFps, Inf, Metric="tanimoto");
metaTbl    = activeTbl(validIdx, :);
%[text] ---
%%
%[text] ## やってみよう 1: IC50 の定義・エルロチニブの構造
%[text] IC50 = 半数阻害濃度（単位: nM または uM）。
%[text] 低い IC50 = より強力（より低い濃度で 50% 阻害を達成）。
logInfo("IC50 = 半数阻害濃度");
logInfo("  IC50 が低いほど = より強力な阻害剤");
logInfo("  エルロチニブ（野生型 EGFR）: IC50 ~2 nM  -->  pIC50 = %.2f", 9 - log10(2));

%[text] エルロチニブを描画してファーマコフォア要素を特定する
emk.viz.draw2d(mol_erl, Title="エルロチニブ -- EGFR 阻害剤（第 1 世代）");

%[text] A1: IC50 = 半数阻害濃度（nM / uM）。
%[text]     IC50 が低いほど優れている: より低い用量で 50% 阻害を達成できる。
%[text]
%[text] A2: エルロチニブの構造的特徴:
%[text]   (a) アルキン（-C#C-）: ATP 部位の裏側のリボースポケットに結合。
%[text]       形状相補性と疎水性接触を提供する。
%[text]   (b) キナゾリンコア（二環、窒素 2 個）: ヒンジ結合モチーフ。
%[text]       キナゾリンの N1 がキナーゼヒンジ Met769 と H 結合を供与/受容する。
%[text]   (c) -OCCOC- チェーン 2 本（C6 と C7）: 水溶性を改善し、
%[text]       ゲートキーパー領域との接触を通じて選択性に寄与する。
%%
%[text] ## やってみよう 2: ChEMBL ターゲット ID・複数生物種

logInfo("EGFR（Homo sapiens）: %s", EGFR_CHEMBL_ID);

%[text] COX-2 検索 -- ChEMBL の優先名が異なる場合があるため複数クエリを試す
cox2Queries = ["Cyclooxygenase-2", "cyclooxygenase 2", "prostaglandin"];
humanCox2 = table();
for qi = 1:numel(cox2Queries)
    cox2Tbl = emk.db.searchChemblTarget(cox2Queries(qi), MaxRows=10);
    humanCox2 = cox2Tbl(cox2Tbl.Organism == "Homo sapiens", :);
    if ~isempty(humanCox2); break; end
end
if ~isempty(humanCox2)
    logInfo("COX-2（Homo sapiens）ChEMBL ID: %s", humanCox2.TargetChEMBLID(1));
else
    logInfo("COX-2 ヒトターゲットが名前検索で見つからない。");
    logInfo("  直接参照: CHEMBL230（COX-2 ヒトの既知安定 ID）");
end

%[text] A: EGFR ChEMBL ID = CHEMBL203。
%[text]    複数の EGFR エントリが存在するのは ChEMBL が各生物種（ヒト・マウスなど）の
%[text]    ターゲットを個別に収録するため。ヒト（Homo sapiens）エントリが
%[text]    最多のバイオアクティビティデータを持つ。
%[text]    COX-2（Homo sapiens）期待値: CHEMBL230。
%%
%[text] ## やってみよう 3: 100 nM 未満の化合物数・最強活性・Ki 比較

nBelow100 = sum(actTbl.Value_nM <= 100);
logInfo("IC50 < 100 nM の化合物: %d / %d", nBelow100, height(actTbl));

[minIC50, minRow] = min(actTbl.Value_nM);
logInfo("最強活性化合物: %s  IC50=%.2f nM", ...
    actTbl.MoleculeChEMBLID(minRow), minIC50);

%[text] EGFR の Ki 値
kiTbl = emk.db.getChemblActivity(EGFR_CHEMBL_ID, ActivityType="Ki", MaxRows=25);
logInfo("取得した Ki レコード: %d", height(kiTbl));
if height(kiTbl) > 0
    logInfo("  Ki 範囲: %.1f - %.1f nM", min(kiTbl.Value_nM), max(kiTbl.Value_nM));
end

%[text] A: 50 件データセットの最強活性化合物は通常サブ nM 範囲。
%[text]    ChEMBL は結果をランダム順で返す --「最強活性」は明示的ソートなしでは実行ごとに変わる。
%[text]    Ki（阻害定数）は酵素濃度に依存しない熱力学的平衡指標。
%[text]    IC50 はアッセイ条件に依存する。
%[text]    競合阻害剤では: \\$\\text{IC}\_{50} = K_i (1 + [\\text{S}]/K_m)\\$。\\$K_i \\le \\text{IC}\_{50}\\$。
%%
%[text] ## やってみよう 4: エルロチニブの pIC50・分布形状・対数スケールの根拠

pIC50_erl_2nM = 9 - log10(2);
logInfo("エルロチニブ pIC50（IC50=2 nM）: %.2f", pIC50_erl_2nM);

skewness_val = skewness(actTbl.pIC50);
logInfo("pIC50 分布の歪度: %.3f", skewness_val);

%[text] A: エルロチニブ ~2 nM の pIC50 = 9 - log10(2) = 8.70。
%[text]    pIC50 >= 8 は高活性（10 nM 未満）として評価される。
%[text]
%[text]    ChEMBL の IC50 分布は pIC50 空間で正の歪みを示す（超強力な外れ値が
%[text]    分布を右に引く）。生の IC50 nM 空間では同データが右偏り（弱活性化合物が多く
%[text]    強活性化合物が少ない）という古典的な形状になる。
%[text]    正確な方向は ChEMBL が返す 50 件のレコードに依存する。
%[text]    医薬化学者が pIC50 を使う理由:
%[text]      (a) ダイナミックレンジが 6 桁にわたる -- 線形軸では実用的でない。
%[text]      (b) 構造活性相関（SAR）は対数空間で加算的:
%[text]          メチル基 1 つが通常 ~0.3 pIC50 単位を追加する。
%[text]      (c) 結合自由エネルギー（\\$\\Delta G = RT \\ln K_i\\$）は \\$\\text{pK}_i\\$ に比例するため、
%[text]          pIC50 は熱力学的量により近い。
%%
%[text] ## やってみよう 5: 活性化合物数・最強活性名・カットオフ緩和

logInfo("活性化合物（IC50 <= %d nM）: %d", IC50_CUTOFF_NM, height(activeTbl));
logInfo("最強活性: %s  IC50=%.2f nM  pIC50=%.2f", ...
    activeTbl.MoleculeChEMBLID(1), activeTbl.Value_nM(1), activeTbl.pIC50(1));

%[text] カットオフを 1000 nM に緩和する
activeTbl_1uM = actTbl(actTbl.Value_nM <= 1000, :);
logInfo("活性化合物（IC50 <= 1000 nM）: %d", height(activeTbl_1uM));

%[text] A: 1000 nM に緩和すると「中程度」効果量の層が加わる。
%[text]    実際の創薬では 1 uM 化合物を一次スクリーニングの「ヒット」、
%[text]    100 nM 化合物を最適化の「リード」と呼ぶことが多い。
%[text]    適切なカットオフはターゲットとアッセイに依存する。
%%
%[text] ## やってみよう 6: SDF ファイル構造・再読み込み・レコード数
%[text] 再エクスポート（runDir があれば）; なければ再読み込みのコンセプトのみ示す
if numel(validMols) > 0
    runDir  = makeRunDir("Prefix", "s07_egfr_ans");
    sdfPath = fullfile(runDir, "egfr_actives_ans.sdf");
    emk.io.writeSdf(validMols, sdfPath);
    logInfo("SDF を書き込み: %s", sdfPath);

    % $$$$ 区切り数 = レコード数
    txt       = fileread(sdfPath);
    lines     = strsplit(txt, newline);
    nRecords  = sum(strcmp(lines, "$$$$"));
    logInfo("SDF レコード数（'$$$$' で計測）: %d", nRecords);

    % numel(validMols) と一致することを確認
    logInfo("期待レコード数: %d  一致: %d", numel(validMols), nRecords == numel(validMols));
end

%[text] A: SDF レコードは 3 つのパートで構成される:
%[text] 1. MOL ブロック: ヘッダー（3 行）+ 結合テーブル（原子/結合行）+ "M  END"
%[text] 2. オプションの SD データフィールド: "> <フィールド名>" の後に値
%[text] 3. レコード区切り: 単独行の "$$$$"
%[text]
%[text]   RDKit は SDF 書き込み時に SMILES から 2D 座標を自動生成するため、
%[text]   出力ファイルは ChemDraw や Maestro で即座に確認できる。
%[text]
%[text]   "$$$$ " 行を数えることで構造レコードの正確な件数が得られる。
%%
%[text] ## やってみよう 7: 上位タニモト・化合物照合・ゲフィチニブ比較

T_top = rankResult.Scores(1);
top_chembl = metaTbl.MoleculeChEMBLID(rankResult.Indices(1));
top_ic50   = metaTbl.Value_nM(rankResult.Indices(1));
logInfo("エルロチニブに対する上位化合物: %s  T=%.4f  IC50=%.1f nM", ...
    top_chembl, T_top, top_ic50);

%[text] ゲフィチニブ比較
GEFITINIB_SMILES = "COc1cc2ncnc(Nc3ccc(F)c(Cl)c3)c2cc1OCCCN1CCOCC1";
mol_gef = emk.mol.fromSmiles(GEFITINIB_SMILES);
fp_gef  = emk.fingerprint.morgan(mol_gef, Radius=2, NBits=2048);
res_gef = emk.similarity.rankBy(fp_gef, dbFps, Inf, Metric="tanimoto");

logInfo("ゲフィチニブ上位 3 位（ChEMBLID / T / IC50）:");
for k = 1:min(3, numel(res_gef.Indices))
    idx = res_gef.Indices(k);
    logInfo("  %d. %s  T=%.4f  IC50=%.1f nM", k, ...
        metaTbl.MoleculeChEMBLID(idx), res_gef.Scores(k), metaTbl.Value_nM(idx));
end

%[text] A: 上位化合物とエルロチニブのタニモトスコアが > 0.7 であれば構造的に非常に近い。
%[text]    エルロチニブとゲフィチニブはともにキナゾリンコアを持つ EGFR 阻害剤だが、
%[text]    置換基が異なるため ECFP4 スクリーニングで必ずしも同一の上位化合物を特定しない。
%[text]    両クエリで一致する化合物はキナゾリンコアを持つ可能性が最も高い。
%[text]    上位化合物が異なる場合は EGFR 阻害剤ケミカルスペースの多様性を反映している。

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
