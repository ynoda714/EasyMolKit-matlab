%[text] # F03 解答: フィンガープリント
%[text] f03_fingerprints.m の演習の参照解答。
addpath(genpath("src"));
emk.setup.initPython();
%%
%[text] ## 解答 E1: Morgan vs MACCS ON ビット数の比較

smiles_list = {"CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", ...
               "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", ...
               "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"};
mol_names   = {"エタノール","ベンゼン","アスピリン","カフェイン","イブプロフェン"};

on_bits = zeros(5, 2);
for i = 1:5
    mol = emk.mol.fromSmiles(smiles_list{i});
    fp_morgan = emk.fingerprint.morgan(mol);
    fp_maccs  = emk.fingerprint.maccs(mol);
    on_bits(i, 1) = sum(emk.fingerprint.toArray(fp_morgan));
    on_bits(i, 2) = sum(emk.fingerprint.toArray(fp_maccs));
end

logInfo("%-14s  %8s  %6s", "分子名", "Morgan", "MACCS");
logInfo("  %s", repmat("-", 1, 32));
for i = 1:5
    logInfo("%-14s  %8d  %6d", mol_names{i}, on_bits(i,1), on_bits(i,2));
end
logInfo("平均 ON ビット -- Morgan: %.1f、MACCS: %.1f", ...
    mean(on_bits(:,1)), mean(on_bits(:,2)));
%[text] Morgan（2048 ビット、Radius 2）は薬物様分子で通常 20～60 の ON ビットを生成します。
%[text] MACCS（167 ビット）は通常 30～60 の ON ビットを生成します。
%[text] MACCS はサイズが小さいため密度（ON/合計）が高くなります。
%%
%[text] ## 解答 E2: Morgan FP 密度 vs 重原子数

data = readtable("data/list/everyday_chemicals.csv", TextType="string");
mols = cellfun(@(s) emk.mol.fromSmiles(s), cellstr(data.SMILES), ...
    "UniformOutput", false);

n = numel(mols);
density    = zeros(n, 1);
heavy_count = zeros(n, 1);

for i = 1:n
    fp   = emk.fingerprint.morgan(mols{i});
    bits = emk.fingerprint.toArray(fp);
    density(i) = sum(bits) / numel(bits);
    d = emk.descriptor.calculate(mols{i}, ["HeavyAtomCount"]);
    heavy_count(i) = d.HeavyAtomCount;
end

figure("Name", "FP 密度 vs 重原子数");
scatter(heavy_count, density, 60, "filled");
xlabel("重原子数");
ylabel("Morgan FP 密度 (ON ビット / 2048)");
title("日用化学品 -- フィンガープリント密度");
%[text] 大きな分子ほど ON ビットが多い傾向があるため、密度は HeavyAtomCount と
%[text] ともに一般的に増加します。ただし 2048 ビットの固定サイズで飽和するため、それほど大きい分子では増加が緩やかになります。
%%
%[text] ## 解答 E3: 思考問題（コードなし）
%[text] 大きな Radius はより拡張された円形環境をエンコードし、
%[text] より多くの固有ハッシュ値を生成する。ただし、すべてのハッシュは
%[text] 2048 ビットの固定ベクトルに折り畳まれるため、異なる環境が同じビットに
%[text] 写像される「衝突」が発生し、ON ビット数が予想より少なくなることがある。
%[text] 非常に大きな Radius（例: 小分子で radius=6）では、ほとんどの原子が
%[text] 分子全体を見るため、異なる環境が少なくなり、radius=2 より
%[text] ON ビット数が実際に減少することがある。
logInfo("E3 は思考問題 -- 解答ファイルのコメントを参照。");

%[text] ---
%[text] やってみようの解答
%[text] ---
%%
%[text] ## やってみよう 1: Radius (0, 1, 2, 3) に対するアスピリンの ON ビット数

mol_asp = emk.mol.fromSmiles("CC(=O)Oc1ccccc1C(=O)O");

radii   = [0, 1, 2, 3];
on_bits = zeros(1, numel(radii));
for i = 1:numel(radii)
    fp = emk.fingerprint.morgan(mol_asp, Radius=radii(i));
    on_bits(i) = sum(emk.fingerprint.toArray(fp));
end
logInfo("アスピリン Morgan ON ビット数:");
for i = 1:numel(radii)
    logInfo("  Radius %d: %d ON ビット", radii(i), on_bits(i));
end
%[text] ON ビット数は Radius 0 から 1 にかけて大幅に増加し、
%[text] その後 Radius 2→3 での増加は緩やかになります（ハッシュ衝突のため）。
%%
%[text] ## やってみよう 2: カフェインのモルガン FP; ON ビットのインデックスを find で取得

fp   = emk.fingerprint.morgan(emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C"));
bits = emk.fingerprint.toArray(fp);
on_bits_idx = find(bits);
logInfo("カフェインの ON ビット数: %d", numel(on_bits_idx));
logInfo("ON ビットインデックス（最初の 10 個）: %s", mat2str(on_bits_idx(1:min(10,end))));
%[text] インデックスは 2048 ビットの空間全体に散らばっています。
%[text] これはハッシュ関数が「隣接する環境を広く分散させる」ように設計されているためです。

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
