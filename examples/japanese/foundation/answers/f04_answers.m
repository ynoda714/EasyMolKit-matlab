%[text] # F04 解答: 類似度
%[text] f04_similarity.m の演習の参照解答。
addpath(genpath("src"));
emk.setup.initPython();
%%
%[text] ## 解答 E1: 鎮痛薬の類似度行列

smiles = {"CC(=O)Oc1ccccc1C(=O)O", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O", ...
          "CC(c1ccc2cc(OC)ccc2c1)C(=O)O", "CC(=O)NC1=CC=C(C=C1)O", ...
          "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"};
names = {"アスピリン","イブプロフェン","ナプロキセン","パラセタモール","カフェイン"};

fps = cellfun(@(s) emk.fingerprint.morgan(emk.mol.fromSmiles(s)), smiles, ...
    "UniformOutput", false);
S = emk.similarity.matrix(fps);

figure("Name", "鎮痛薬類似度", "Position", [100 100 500 440]);
imagesc(S);
colorbar;
clim([0, 1]);
colormap(parula);
xticks(1:5);  xticklabels(names);  xtickangle(30);
yticks(1:5);  yticklabels(names);
title("タニモト類似度（Morgan ECFP4）");
for r = 1:5
    for c = 1:5
        text(c, r, sprintf("%.2f", S(r,c)), ...
            "HorizontalAlignment","center","FontSize",8,"Color","white");
    end
end

%[text] 対角外で最も類似したペアを求める
S_masked = S - eye(5);
[max_val, lin_idx] = max(S_masked(:));
[r_idx, c_idx] = ind2sub([5 5], lin_idx);
logInfo("最も類似したペア: %s -- %s  (タニモト = %.4f)", ...
    names{r_idx}, names{c_idx}, max_val);
%[text] 期待値: イブプロフェンとナプロキセンはどちらもプロフェン系 NSAID で
%[text] カルボン酸 + 芳香族スキャフォールドが類似しており、最高スコアになる。
%%
%[text] ## 解答 E2: everyday_chemicals でイブプロフェンに最も類似した上位 3 件

query_fp = emk.fingerprint.morgan(emk.mol.fromSmiles( ...
    "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"));
data    = readtable("data/list/everyday_chemicals.csv", TextType="string");
db_mols = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(data.SMILES), ...
    "UniformOutput", false);
db_fps  = cellfun(@(m) emk.fingerprint.morgan(m), db_mols, "UniformOutput", false);

result = emk.similarity.rankBy(query_fp, db_fps, 3);
logInfo("イブプロフェンに最も類似した上位 3 件:");
for k = 1:3
    logInfo("  %d. %-20s  タニモト=%.4f", k, data.CommonName(result.Indices(k)), ...
        result.Scores(k));
end
%%
%[text] ## 解答 E3: タニモト閾値を超えるペアをカウント

threshold = 0.3;
count = 0;
for r = 1:5
    for c = (r+1):5   % 上三角のみ
        if S(r, c) >= threshold
            count = count + 1;
            logInfo("  %s -- %s  (%.4f)", names{r}, names{c}, S(r,c));
        end
    end
end
logInfo("タニモト >= %.1f のペア: %d 件", threshold, count);
%[text] プロフェン系薬物（アスピリン、イブプロフェン、ナプロキセン）は
%[text] カルボン酸 + 芳香族環スキャフォールドを共有するため 0.3 を超える傾向がある。
%[text]
%[text] ---
%[text] やってみようの解答
%[text] ---
%%
%[text] ## やってみよう 1: エタノールとメタノール（"CO"）のタニモト類似度

mol_eth  = emk.mol.fromSmiles("CCO");
mol_meth = emk.mol.fromSmiles("CO");
mol_prop = emk.mol.fromSmiles("CCCO");

fp_eth  = emk.fingerprint.morgan(mol_eth);
fp_meth = emk.fingerprint.morgan(mol_meth);
fp_prop = emk.fingerprint.morgan(mol_prop);

t_eth_meth = emk.similarity.tanimoto(fp_eth, fp_meth);
t_eth_prop = emk.similarity.tanimoto(fp_eth, fp_prop);

logInfo("エタノール vs メタノール  = %.4f", t_eth_meth);
logInfo("エタノール vs プロパノール = %.4f", t_eth_prop);
%[text] メタノールはエタノールより小さく、Morgan 環境を少ししか共有しない。
%[text] エタノール vs プロパノールは類似した鎖延長で高スコアになる可能性が高い。
%%
%[text] ## やってみよう 2: ダイス係数による類似度行列

smiles_6 = ["CCO", "CCCO", "CC(C)O", "CC(=O)Oc1ccccc1C(=O)O", ...
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "CC(=O)NC1=CC=C(C=C1)O"];
mols_6 = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(smiles_6), "UniformOutput", false);
fps_6  = cellfun(@(m) emk.fingerprint.morgan(m), mols_6, "UniformOutput", false);

S      = emk.similarity.matrix(fps_6);
S_dice = emk.similarity.matrix(fps_6, Metric="dice");

logInfo("すべての要素で Dice >= Tanimoto: %d", all(S_dice(:) >= S(:)));

S_off      = S      - eye(numel(fps_6));
S_dice_off = S_dice - eye(numel(fps_6));
[~, idx_t] = max(S_off(:));
[~, idx_d] = max(S_dice_off(:));
logInfo("最も類似したペア -- タニモト idx: %d、ダイス idx: %d  (同一: %d)", ...
    idx_t, idx_d, isequal(idx_t, idx_d));
%[text] 2 値ビットベクトルでは常に Dice >= Tanimoto が成立します（Dice = 2J/(A+B)、Tanimoto = J/(A+B-J)）。
%[text] 分子ペアのランキングはメトリックが異なっても通常は変わりません。
%%
%[text] ## やってみよう 3: イブプロフェンをクエリとして検索

query_mol = emk.mol.fromSmiles("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");   % イブプロフェン
query_fp  = emk.fingerprint.morgan(query_mol);

data    = readtable("data/list/everyday_chemicals.csv", TextType="string");
valid_mask_db = cellfun(@(s) emk.mol.isValid(s), cellstr(data.SMILES));
db_smiles = cellstr(data.SMILES(valid_mask_db));
db_mols   = cellfun(@(s) emk.mol.fromSmiles(s), db_smiles, "UniformOutput", false);
db_fps    = cellfun(@(m) emk.fingerprint.morgan(m), db_mols, "UniformOutput", false);
db_names  = data.CommonName(valid_mask_db);

result = emk.similarity.rankBy(query_fp, db_fps, 5);
logInfo("イブプロフェンに最も類似した日用化学品 上位 5 件:");
for k = 1:5
    idx   = result.Indices(k);
    label = db_names(idx);
    if strlength(label) == 0; label = "(unknown)"; end
    logInfo("  %d. %-20s  タニモト=%.4f", k, label, result.Scores(k));
end
%[text] 芳香族カルボン酸スキャフォールド（安息香酸誘導体など）を持つ化合物が
%[text] 上位に来る傾向があります。

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
