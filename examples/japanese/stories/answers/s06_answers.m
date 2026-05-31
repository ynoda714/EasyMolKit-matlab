%[text] # S06 解答: PubChem データベース検索
%[text] s06_pubchem_search.m の「やってみよう」演習の参考解答。
%[text] まず s06_pubchem_search.m を実行してから、このファイルで答え合わせをすること。
%[text] 注意: このファイルはインターネット接続が必要（PubChem REST API）。
addpath(genpath("src"));
emk.setup.initPython();
logInfo("S06 解答: セットアップ完了（インターネット接続必要）");
%%
%[text] ## やってみよう 1: アスピリン CID・パラセタモル == アセトアミノフェン？

tbl_aspirin = emk.db.searchPubchem("aspirin");
logInfo("アスピリン CID: %d", tbl_aspirin.CID(1));
logInfo("アスピリン IUPAC 名: %s", tbl_aspirin.IUPACName(1));

tbl_para = emk.db.searchPubchem("paracetamol");
tbl_acet = emk.db.searchPubchem("acetaminophen");
logInfo("パラセタモル    CID: %d", tbl_para.CID(1));
logInfo("アセトアミノフェン CID: %d", tbl_acet.CID(1));
logInfo("同じ化合物? %d", tbl_para.CID(1) == tbl_acet.CID(1));

%[text] A: アスピリン CID = 2244。この安定な識別子はデータベース開始以来アスピリンのレコード。
%[text]    IUPAC 名: "2-acetyloxybenzoic acid" -- オルト位にアセチルエステルを持つ安息香酸。
%[text]    パラセタモル == アセトアミノフェン: 両方とも CID 1983（同じ化合物、地域名の違い）。
%%
%[text] ## やってみよう 2: SMILES でイブプロフェン検索・(S) 体 CID

tbl_ibu = emk.db.searchPubchem("CC(C)Cc1ccc(cc1)C(C)C(=O)O", Type="smiles");
logInfo("イブプロフェン（SMILES）CID: %d", tbl_ibu.CID(1));
logInfo("PubChem 正規 SMILES: %s", tbl_ibu.IsomericSMILES(1));

%[text] (S)-イブプロフェン
tbl_s_ibu = emk.db.searchPubchem("[C@@H](C(=O)O)(Cc1ccc(cc1)CC(C)C)C", Type="smiles");
logInfo("(S)-イブプロフェン CID: %d", tbl_s_ibu.CID(1));
logInfo("ラセミ体と同じ? %d", tbl_ibu.CID(1) == tbl_s_ibu.CID(1));

%[text] A: イブプロフェン（ラセミ）CID = 3672。
%[text]    PubChem は自身の正規 SMILES を返す（入力 SMILES と原子順序が異なる場合があるが同一分子）。
%[text]    (S)-イブプロフェンは立体化学が異なる構造として PubChem に登録されるため
%[text]    異なる CID を持つ（PubChem はデフォルトで立体異性体を区別）。
%%
%[text] ## やってみよう 3: 解熱鎮痛薬パネル -- 最重量・最小・硫黄含有化合物

PANEL_NAMES = ["aspirin", "ibuprofen", "acetaminophen", "naproxen", "celecoxib"];
N = numel(PANEL_NAMES);

cids    = zeros(N, 1, 'uint32');
iupac   = strings(N, 1);
formulas= strings(N, 1);
mws     = zeros(N, 1);
panelSmiles = strings(N, 1);

for i = 1:N
    try
        r = emk.db.searchPubchem(PANEL_NAMES(i));
        cids(i)     = r.CID(1);
        iupac(i)    = r.IUPACName(1);
        formulas(i) = r.MolecularFormula(1);
        mws(i)      = r.MolecularWeight(1);
        panelSmiles(i) = r.IsomericSMILES(1);
    catch ME
        logWarn("  %s: クエリ失敗 -- %s", PANEL_NAMES(i), ME.message);
    end
end

panelTbl = table(PANEL_NAMES', cids, iupac, formulas, mws, panelSmiles, ...
    'VariableNames', {'Name','CID','IUPACName','MolecularFormula','MolecularWeight','SMILES'});

%[text] 最重量
[~, iHeavy] = max(panelTbl.MolecularWeight);
logInfo("最重量: %s (MW=%.1f)", panelTbl.Name(iHeavy), panelTbl.MolecularWeight(iHeavy));

%[text] 硫黄含有
hasSulfur = contains(panelTbl.MolecularFormula, "S");
logInfo("硫黄含有: %s", strjoin(panelTbl.Name(hasSulfur), ", "));

%[text] A: 最重量: セレコキシブ（~381 Da）。最新の選択的 COX-2 阻害薬は
%[text]    従来の NSAID より一般に大きく複雑。
%[text]    アスピリンとアセトアミノフェンはともに MW < 180 Da だが作用機序が全く異なる --
%[text]    MW は機序を予測しない。
%[text]    硫黄含有: セレコキシブ（C17H14F3N3O2S）-- スルホンアミド（SO2NH2）が COX-2 選択性の鍵。
%%
%[text] ## やってみよう 4: 完全物性テーブル・Ro5 確認・TPSA ランク付け
%[text] RDKit 記述子を構築する
descNames = ["MolWt","LogP","TPSA","NumHDonors","NumHAcceptors","NumRotatableBonds","RingCount"];
mwArr  = zeros(N,1);  lpArr  = zeros(N,1);  tpArr  = zeros(N,1);
hdArr  = zeros(N,1);  haArr  = zeros(N,1);  rbArr  = zeros(N,1); rcArr  = zeros(N,1);

for i = 1:N
    if panelSmiles(i) ~= ""
        mol = emk.mol.fromSmiles(panelSmiles(i));
        d   = emk.descriptor.calculate(mol, descNames);
        mwArr(i) = d.MolWt;  lpArr(i) = d.LogP;  tpArr(i) = d.TPSA;
        hdArr(i) = d.NumHDonors;  haArr(i) = d.NumHAcceptors;
        rbArr(i) = d.NumRotatableBonds;  rcArr(i) = d.RingCount;
    end
end

fullTbl = table(PANEL_NAMES', cids, formulas, mwArr, lpArr, tpArr, ...
    hdArr, haArr, rbArr, rcArr, ...
    'VariableNames', {'Name','CID','MolecularFormula','MolecularWeight', ...
                      'LogP','TPSA','NumHDonors','NumHAcceptors', ...
                      'NumRotatableBonds','RingCount'});

%[text] Ro5 チェック
roTbl  = table(fullTbl.MolecularWeight, fullTbl.LogP, fullTbl.NumHDonors, fullTbl.NumHAcceptors, ...
    'VariableNames', {'MolWt','LogP','NumHDonors','NumHAcceptors'});
result = emk.filter.lipinski(roTbl);
logInfo("Ro5 結果:");
for i = 1:N
    logInfo("  %-15s  合格=%d  違反=%d", fullTbl.Name(i), result.Pass_Ro5(i), result.Violations_Ro5(i));
end

%[text] 最高 TPSA・最低 LogP
[~, iMaxTPSA] = max(fullTbl.TPSA);
[~, iMinLogP] = min(fullTbl.LogP);
logInfo("最高 TPSA: %s (%.1f A^2)", fullTbl.Name(iMaxTPSA), fullTbl.TPSA(iMaxTPSA));
logInfo("最低 LogP: %s (%.2f)", fullTbl.Name(iMinLogP), fullTbl.LogP(iMinLogP));

%[text] A: 5 種すべて Ro5 合格（いずれも広く使われる経口鎮痛薬）。
%[text]    最高 TPSA: セレコキシブ（~78 A^2、スルホンアミド＋ピラゾール）。
%[text]    最低 LogP: アスピリン（~1.31）-- パネル内で最親水性。
%%
%[text] ## やってみよう 5: 構造の視覚的調査 -- セレコキシブの CF3

mol_cel = emk.mol.fromSmiles(panelSmiles(panelTbl.Name == "celecoxib"));
emk.viz.draw2d(mol_cel, Title="セレコキシブ（COX-2 選択的阻害薬）", ...
    Width=400, Height=400);

%[text] Q: なぜ CF3 を導入するのか？
%[text] A: CF3（トリフルオロメチル）は親油性を高め（1 つで LogP +0.5）、
%[text]    代謝安定性を改善し（その部位での CYP 酸化を阻止）、
%[text]    代謝を受けずにメチル基を立体的に模倣できる。
%[text]    医薬化学で最も頻用される置換パターンの一つ。
%%
%[text] ## やってみよう 6: CNS プロファイル・セレコキシブ TPSA・MW 棒グラフ
%[text] CNS ターゲットプロファイル: TPSA < 60、LogP 1-3
cnsFit = fullTbl.TPSA < 60 & fullTbl.LogP >= 1 & fullTbl.LogP <= 3;
logInfo("CNS プロファイル（TPSA<60, LogP 1-3）: %s", strjoin(fullTbl.Name(cnsFit), ", "));

%[text] A: イブプロフェン（TPSA ~37, LogP ~3.1）とアスピリン（TPSA ~64, LogP ~1.31）が
%[text]    CNS プロファイルに最も近い。どちらも主要 CNS 薬として使用されていないが、
%[text]    構造的には BBB を透過できる。
%[text]    セレコキシブは TPSA ~78 A^2 で高いが、炎症部位の末梢 COX-2 を標的とするため許容される。
%[text]
%[text] MW 棒グラフ（リピンスキー上限付き）
figure("Name", "解熱鎮痛薬パネル -- MW とリピンスキー上限", "Position", [100 100 560 380]);
bar(fullTbl.MolecularWeight);
xticks(1:N);  xticklabels(fullTbl.Name);  xtickangle(20);
yline(500, "--r", "リピンスキー MW=500 上限");
ylabel("分子量 (Da)");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
