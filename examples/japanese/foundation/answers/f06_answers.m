%[text] # F06 解答: ファイル I/O
%[text] `f06_file_io.m` の演習問題（E1〜E3）の参照解答です。
addpath(genpath("src"));
emk.setup.initPython();
runDir = makeRunDir();
%%
%[text] ## 解答 E1: SMILES → SDF ラウンドトリップ検証

my_smiles = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", ...
             "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"];
my_names  = ["ethanol", "benzene", "aspirin", "caffeine", "ibuprofen"];

%[text] SMILES リストを書き出します。
smi_path = fullfile(runDir, "e1_molecules.smi");
lines = strings(numel(my_smiles), 1);
for i = 1:numel(my_smiles)
    lines(i) = sprintf("%s    %s", my_smiles(i), my_names(i));
end
writelines(lines, smi_path);

%[text] 読み込み、SDF 書き出し、読み直し、内容を比較します。
mols_in  = emk.io.readSmilesList(smi_path);
sdf_path = fullfile(runDir, "e1_output.sdf");
emk.io.writeSdf(mols_in, sdf_path);
mols_out = emk.io.readSdf(sdf_path);

logInfo("--- ラウンドトリップ 正規 SMILES ---");
for i = 1:numel(mols_in)
    smi_before = emk.mol.toSmiles(mols_in{i});
    smi_after  = emk.mol.toSmiles(mols_out{i});
    match = isequal(smi_before, smi_after);
    logInfo("  %-14s  一致=%d  [%s]", my_names(i), match, smi_before);
end
%[text] SDF ラウンドトリップの前後で正規 SMILES はすべて同一になるはずです。
%%
%[text] ## 解答 E2: 日用化学品の Ro5 フィルタと保存

data = readtable("data/list/everyday_chemicals.csv", TextType="string");
valid_mask = cellfun(@(s) emk.mol.isValid(s), cellstr(data.SMILES));
mols_valid = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(data.SMILES(valid_mask)), ...
    "UniformOutput", false);

desc_tbl    = emk.descriptor.batchCalculate(mols_valid, ...
    ["MolWt","LogP","NumHDonors","NumHAcceptors"]);
lipinski_tbl = emk.filter.lipinski(desc_tbl);
pass_mask    = lipinski_tbl.Pass_Ro5;

logInfo("Ro5 を通過した日用化学品: %d / %d", sum(pass_mask), numel(mols_valid));

sdf_out = fullfile(runDir, "everyday_ro5_pass.sdf");
emk.io.writeSdf(mols_valid(pass_mask), sdf_out);
logInfo("保存先: %s", sdf_out);
%[text] ほとんどの日用化学品は小分子（溶剤、香料）であるため、
%[text] 大部分が Ro5 を通過すると予想されます〔30 件中約 25〜30 件が目安です。
%%
%[text] ## 解答 E3: フルミニパイプライン — アスピリン類似 FDA 薬物
%[text] **ステップ a**: FDA 薬物を読み込んで Ro5 フィルタを適用します。
fda = readtable("data/list/fda_drugs.csv", TextType="string");
valid_fda = cellfun(@(s) emk.mol.isValid(s), cellstr(fda.SMILES));
fda_mols  = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(fda.SMILES(valid_fda)), ...
    "UniformOutput", false);
fda_valid = fda(valid_fda, :);

desc_fda     = emk.descriptor.batchCalculate(fda_mols, ...
    ["MolWt","LogP","NumHDonors","NumHAcceptors"]);
lip_fda      = emk.filter.lipinski(desc_fda);
pass_fda     = lip_fda.Pass_Ro5;

passing_mols = fda_mols(pass_fda);
logInfo("Ro5 通過 FDA 分子: %d 件", sum(pass_fda));

%[text] **ステップ b**: すべての通過分子に Morgan フィンガープリントを計算します。
fps_fda = cellfun(@(m) emk.fingerprint.morgan(m), passing_mols, ...
    "UniformOutput", false);

%[text] **ステップ c**: アスピリンとの Tanimoto 類似度でランク付けします。
aspirin_fp = emk.fingerprint.morgan(emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O"));
result = emk.similarity.rankBy(aspirin_fp, fps_fda, 10);

passing_names = fda_valid.Name(pass_fda);
logInfo("アスピリンに最も類似した上位 10 件の FDA 薬物:");
for k = 1:10
    idx = result.Indices(k);
    logInfo("  %d. %-25s  タニモト=%.4f", k, passing_names(idx), result.Scores(k));
end

%[text] **ステップ d**: 上位 10 件を SDF に保存します。
top_mols = passing_mols(result.Indices(1:10));
sdf_top  = fullfile(runDir, "fda_aspirin_similar_top10.sdf");
emk.io.writeSdf(top_mols, sdf_top);
logInfo("上位 10 件を保存した: %s", sdf_top);
%%
%[text] ## やってみよう 1: 壊れた SMILES を readSmilesList に渡した場合の挙動
%[text] `emk.io.readSmilesList` は無効な SMILES をクラッシュせずにスキップし、
%[text] その行番号と SMILES 文字列を `logWarn` に記録します。
%[text] 以下のように無効エントリを含む .smi ファイルを作成して動作を確認できます。

invalid_smi_path = fullfile(runDir, "y1_invalid_test.smi");
writelines([ ...
    "CCO  ethanol", ...
    "INVALID_SMILES  bad_entry", ...
    "c1ccccc1  benzene" ...
], invalid_smi_path);

mols_y1 = emk.io.readSmilesList(invalid_smi_path);
logInfo("やってみよう 1 -- 読み込んだ分子数: %d（無効 1 件はスキップ済み）", numel(mols_y1));
%[text] **観察ポイント**
%[text] - ログ出力に `[WARN]` `readSmilesList: line 2, SMILES 'INVALID_SMILES' skipped` が現れる
%[text] - 関数はクラッシュせず `numel(mols_y1) == 2` になる
%[text] - 最後に `readSmilesList: loaded 2 molecules ... (1 skipped)` と要約が出力される
%%
%[text] ## やってみよう 2: 緩和ルールで通過数はどう変わる？

chembl = readtable("data/list/chembl_phase3plus.csv", TextType="string");
valid_mask_y2 = cellfun(@(s) emk.mol.isValid(s), cellstr(chembl.SMILES));
mols_y2 = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(chembl.SMILES(valid_mask_y2)), ...
    "UniformOutput", false);
chembl_valid = chembl(valid_mask_y2, :);

desc_tbl = emk.descriptor.batchCalculate(mols_y2);
desc_tbl.Name = chembl_valid.Name;

lipinski_tbl         = emk.filter.lipinski(desc_tbl);
lipinski_tbl_relaxed = emk.filter.lipinski(desc_tbl, "MaxViolations", 1);

pass_strict  = lipinski_tbl.Pass_Ro5;
pass_relaxed = lipinski_tbl_relaxed.Pass_Ro5;

logInfo("厳格 Ro5 通過: %d / %d", sum(pass_strict), numel(pass_strict));
logInfo("緩和 Ro5 通過 (MaxViolations=1): %d / %d", sum(pass_relaxed), numel(pass_relaxed));
logInfo("新たに通過した薬物 (%d 件):", sum(pass_relaxed & ~pass_strict));

newly_pass = desc_tbl.Name(pass_relaxed & ~pass_strict);
for k = 1:numel(newly_pass)
    logInfo("  %s", newly_pass(k));
end
%[text] MaxViolations=1 の緩和フィルタにより、通常は数件〜数十件が追加で通過します。
%[text] これらは Ro5 のルールを 1 つだけ僅かに超えた「ほぼ経口投与可能」な化合物です。

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
