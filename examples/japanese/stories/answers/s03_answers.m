%[text] # S03 解答: 構造アラート — 危険な官能基を SMARTS で検出する
%[text] `s03_structure_alerts.m` の「やってみよう」演習の参考解答です。
%[text] まず `s03_structure_alerts.m` を実行してから、このファイルで答え合わせをしてください。
addpath(genpath("src"));
emk.setup.initPython();
logInfo("S03 解答: セットアップ完了");

%[text] ---- アラートパネル（s03 セクション 1 と同じ）--------------------------------
ALERT_NAMES  = ["Epoxide",      "Michael_acceptor", "Aldehyde", ...
                "Acyl_halide",  "Nitro",             "Hydrazine", ...
                "Diazo",        "Thiol"];
ALERT_SMARTS = ["[C]1CO1",      "[C]=[C]-C=O",      "[CH]=O", ...
                "[C](=O)[F,Cl]","[N+](=O)[O-]",     "[NH]-[NH2]", ...
                "[#6]=[N+]=[N-]","[SH]"];

TEST_NAMES  = ["エピクロロヒドリン", "アクロレイン",   "ホルムアルデヒド", ...
               "マロンアルデヒド",   "ニトロベンゼン", "フェニルヒドラジン", ...
               "アスピリン",         "カプトプリル"];
TEST_SMILES = ["ClCC1CO1",          "C=CC=O",       "C=O", ...
               "O=CCC=O",          "c1ccc([N+](=O)[O-])cc1", ...
               "c1ccccc1NN",       "CC(=O)Oc1ccccc1C(=O)O", ...
               "CC(CS)C(=O)N1CCCC1C(=O)O"];

nMols   = numel(TEST_SMILES);
nAlerts = numel(ALERT_NAMES);
mols_test = cell(1, nMols);
for i = 1:nMols
    mols_test{i} = emk.mol.fromSmiles(TEST_SMILES(i));
end
flagMat = false(nMols, nAlerts);
for j = 1:nAlerts
    flagMat(:, j) = emk.mol.hasSubstruct(mols_test, ALERT_SMARTS(j))';
end
%[text] ---
%%
%[text] ## やってみよう 1: 最多アラート分子・アスピリンのアラート

totalAlerts = sum(flagMat, 2);
[maxAl, iMax] = max(totalAlerts);
logInfo("最多アラート: %s (%d 件)", TEST_NAMES(iMax), maxAl);
for i = 1:nMols
    alertList = ALERT_NAMES(flagMat(i, :));
    if numel(alertList) > 0
        logInfo("  %s: %s", TEST_NAMES(i), strjoin(alertList, ", "));
    else
        logInfo("  %s: （アラートなし）", TEST_NAMES(i));
    end
end

%[text] A: アクロレイン（`C=CC=O`）はアルデヒド AND マイケル受容体の両方を発火（2 件）。
%[text]    アスピリンはこのパネルではアラートなし。
%[text]    アスピリンのエステル基（フェノールのアセチル）はこれらの SMARTS には捕捉されない。
%[text]    エステルはインビボで加水分解されてサリチル酸になる意図的なファーマコフォア。
%[text]    「アラートなし」は「安全」を意味するわけではなく、これらの特定パターンがないだけ。
%[text]    カプトプリル（ACE 阻害薬）はチオールを発火 — ACE の亜鈑イオンをキレートするファーマコフォアなので想定内。
%%
%[text] ## やってみよう 2: アクロレインとベータ炭素へのマイケル付加

mol_acrolein = emk.mol.fromSmiles("C=CC=O");
figure("Name", "アクロレイン -- マイケル + アルデヒド");
emk.viz.draw2d(mol_acrolein, Title="アクロレイン — マイケル + アルデヒド");

%[text] Q: システインチオールはどの原子を攻撃するか？
%[text] A: β 炭素（C=C 二重結合の C2 — マイケル受容体サイト）。
%[text]    軟求核剤（チオール）は LUMO が共役系のβ 位に広がるため、
%[text]    1,4 付加（マイケル付加）を 1,2 付加（直接カルボニル攻撃）より好みます。
%%
%[text] ## やってみよう 3: 最初の 3 つの PAINS パターン

fid = fopen("data/list/pains.csv", "r");
textscan(fid, "%s", 1, "Delimiter", "\n");
C = textscan(fid, "%q%q%q%q", "Delimiter", ",");
fclose(fid);
painsCsv = table(string(C{1}), string(C{2}), string(C{3}), ...
    VariableNames=["Name", "SMARTS", "FilterSet"]);
painsCsv = painsCsv(strtrim(painsCsv.SMARTS) ~= "", :);

logInfo("最初の 3 つの PAINS パターン:");
disp(painsCsv(1:3, ["Name","SMARTS","FilterSet"]));

%[text] A: PAINS SMARTS は部分構造マッチング用の難解な再帰クエリ -- 人間が読むためではない。
%[text]    FilterSet ラベル（A/B/C）は Baell & Holloway 2010 論文の
%[text]    3 種類の HTS アッセイタイプに対応する。
%%
%[text] ## やってみよう 4: FDA 承認薬中の PAINS 割合

rawDrugs = readtable("data/list/fda_drugs.csv", "TextType", "string");
nDrugs   = height(rawDrugs);

validDrug = false(nDrugs, 1);
mols_drug = cell(nDrugs, 1);
for i = 1:nDrugs
    if emk.mol.isValid(rawDrugs.SMILES(i))
        mols_drug{i} = emk.mol.fromSmiles(rawDrugs.SMILES(i));
        validDrug(i) = true;
    end
end
validIdx    = find(validDrug);
%[text] PAINS_A（最も厳しい 16 パターン）を使用
painsScreen = painsCsv(painsCsv.FilterSet == "PAINS_A", :);
nPains      = height(painsScreen);

isPainsFlagged = false(nDrugs, 1);
painsCount     = zeros(nDrugs, 1);
firstAlert     = strings(nDrugs, 1);
logInfo("%d 薬物を %d 件の PAINS_A パターンでスクリーニング中...", nDrugs, nPains);
logInfo("  （PAINS_A = 最厳格サブセット; データベース全体では %d 件）", height(painsCsv));
for j = 1:nPains
    hits_j = emk.mol.hasSubstruct(mols_drug(validDrug), painsScreen.SMARTS(j));
    for k = 1:numel(hits_j)
        idx = validIdx(k);
        if hits_j(k)
            isPainsFlagged(idx) = true;
            painsCount(idx)     = painsCount(idx) + 1;
            if firstAlert(idx) == ""
                firstAlert(idx) = painsScreen.Name(j);
            end
        end
    end
    logProgress(j, nPains, "PAINS スクリーニング");
end
nFlagged = sum(isPainsFlagged);
logInfo("PAINS_A フラグ: %d / %d  (%.0f%%)", nFlagged, nDrugs, 100*nFlagged/nDrugs);

%[text] A: PAINS_A（16 パターン）では通常 ~3〜5% の承認薬がフラグ対象。
%[text]    全 480 パターン使用時は 15〜25% に上昇（PAINS_B/C はより寛容なパターン）。
%[text]    承認薬は反応性モチーフ自体がファーマコフォアのこともある
%[text]    （カプトプリルのチオール、ベータラクタム系抗生物質など）。
%[text]    PAINS は調査の引き金であって自動却下の根拠ではない。
%%
%[text] ## やってみよう 5: 最多 PAINS パターン薬物
%[text] （isPainsFlagged / painsCount / firstAlert は上で計算済み）
resultTbl = table(rawDrugs.Name, rawDrugs.SMILES, isPainsFlagged, painsCount, firstAlert, ...
    'VariableNames', ["Name","SMILES","PAINS_Flagged","PAINS_Count","First_Alert"]);
topDrugs = sortrows(resultTbl, "PAINS_Count", "descend");
logInfo("PAINS カウント上位 5 件:");
disp(topDrugs(1:5, ["Name","PAINS_Count","First_Alert"]));

%[text] A: 複数 PAINS ヒットを持つ薬物は官能基が多い複雑な天然物になりがちです。
%[text]    上位薬物の作用機序を調べてみましょう —
%[text]    その「アラート」が実際にその薬が効く理由の核心かもしれません。
%%
%[text] ## やってみよう 6: FDA 薬物へのカスタムアラート適用

%[text] セクション 1 のアラートパネルを FDA 薬物に適用し、フラグ率を PAINS と比較します。
%[text] （`mols_drug`・`validDrug`・`validIdx`・`nDrugs`・`nAlerts`・`ALERT_SMARTS` は上で計算済み）

drugAlertMat = false(nDrugs, nAlerts);
for j = 1:nAlerts
    hits_j = emk.mol.hasSubstruct(mols_drug(validDrug), ALERT_SMARTS(j));
    drugAlertMat(validIdx, j) = hits_j;
end
alertFlagged = sum(drugAlertMat, 2) > 0;
nAlertFlagged = sum(alertFlagged);
logInfo("カスタムアラート: %d / %d 件フラグ (%.0f%%)", ...
    nAlertFlagged, nDrugs, 100*nAlertFlagged/nDrugs);
logInfo("PAINS_A:      %d / %d 件フラグ", nFlagged, nDrugs);

%[text] A: カスタムアラート（8 種）は決まった化学構造パターンをちょうど特定するため、
%[text]    FDA 薬物へのヒット率は PAINS_A と異なります。
%[text]    カプトプリル（Thiol）やベータラクタム系抗生物質（アシルハライド）がアラートを発火するため、
%[text]    フラグ率は PAINS より高くなる可能性があります。
%[text]    どちらのアプローチも「フラグ = 調査のトリガー」であり、自動除外の根拠ではありません。

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
