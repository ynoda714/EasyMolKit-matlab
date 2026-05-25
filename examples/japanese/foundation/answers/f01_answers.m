%[text] # F01 解答: SMILES から分子を描画する
%[text] `f01_draw_molecules.m` の演習の参照解答。
addpath(genpath("src"));
emk.setup.initPython();
%%
%[text] ## 解答 E1: ニコチンを描いてみましょう
mol_nic = emk.mol.fromSmiles("CN1CCC[C@H]1C2=CN=CC=C2");
figure("Name", "ニコチン");
emk.viz.draw2d(mol_nic, Title="Nicotine");
%[text] RDKit は 2 つの環（ピロリジン環とピリジン環）を描きます。
%[text] 不斉中心（@H）は SMILES にエンコードされていますが、2D 描画では
%[text] くさび形結合として明示されない場合があります。
%%
%[text] ## 解答 E2: テルペンを並べて比較しましょう

smiles_pair = ["OC1CC(CC(C1)C(C)C)C", "Cc1ccc(cc1O)C(C)C"];
names_pair  = ["Menthol", "Thymol"];

figure("Name", "メントール vs チモール", "Position", [100 100 800 400]);
for i = 1:2
    subplot(1, 2, i);
    mol = emk.mol.fromSmiles(smiles_pair(i));
    emk.viz.draw2d(mol, Title=names_pair(i));
    can = emk.mol.toSmiles(mol);
    logInfo("%s の正規 SMILES: %s", names_pair(i), can);
end
%[text] メントールは飽和環を持つシクロヘキサノール（C10H20O）。
%[text] チモールはフェノール（芳香族、C10H14O）。
%[text] 分子式は異なるため正規 SMILES も異なります。
%%
%[text] ## 解答 E3: 有効性フィルターを作ってみましょう

candidates = {"CCO", "C(C", "c1ccccc1", "BADSMILES", "CC(=O)O"};
valid_smiles = {};
for i = 1:numel(candidates)
    if emk.mol.isValid(candidates{i})
        valid_smiles{end+1} = candidates{i}; %#ok<SAGROW>
    end
end
logInfo("有効な SMILES (%d / %d):", numel(valid_smiles), numel(candidates));
for i = 1:numel(valid_smiles)
    logInfo("  %s", valid_smiles{i});
end
%[text] 有効なもの: "CCO"、"c1ccccc1"、"CC(=O)O"（5 件中 3 件）
%[text] "C(C" は括弧が閉じておらず構文エラー。"BADSMILES" は SMILES 文字列ではありません。

%\[text] ---
%[text] やってみようの解答
%\[text] ---
%%
%[text] ## やってみよう 1: イブプロフェンとニコチンのクラス確認

mol_ibupro = emk.mol.fromSmiles("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O");
mol_nicot  = emk.mol.fromSmiles("CN1CCC[C@H]1C2=CN=CC=C2");
logInfo("イブプロフェン: class = %s", class(mol_ibupro));
logInfo("ニコチン:       class = %s", class(mol_nicot));
%[text] どちらも `py.rdkit.Chem.rdchem.Mol` が返ります。
%[text] MATLAB は Python オブジェクトへの参照を保持するだけなので、
%[text] RDKit の API を意識せずに `emk.*` 関数の入力として使えます。
%%
%[text] ## やってみよう 2: ピリジンとフランの正規 SMILES

mol_pyr  = emk.mol.fromSmiles("c1ccncc1");   % ピリジン
mol_fur  = emk.mol.fromSmiles("c1ccoc1");    % フラン

can_pyr = emk.mol.toSmiles(mol_pyr);
can_fur = emk.mol.toSmiles(mol_fur);
logInfo("ピリジンの正規 SMILES: %s", can_pyr);
logInfo("フランの正規 SMILES:   %s", can_fur);
%[text] RDKit 正規形: ピリジン → `c1ccncc1`（芳香族、小文字 n を保持）
%[text] フラン → `c1ccoc1`（芳香族、小文字 o を保持）
%[text] RDKit はヒュッケル則（4n+2 電子）を使って両方を芳香族と判定します。
%%
%[text] ## やってみよう 3: "C(C)(C)(C)(C)C" が無効な理由

tf_neo  = emk.mol.isValid("C(C)(C)(C)(C)C");
logInfo("C(C)(C)(C)(C)C  isValid = %d", tf_neo);
%[text] isValid = 0（false）。SMILES の構文としては 5 配位炭素を書くことができますが、
%[text] RDKit の isValid は化学的妥当性（原子価ルール）まで確認します。
%[text] 炭素の最大結合数は 4 であるため、RDKit はこの分子を無効と判定します。

tf_bad  = emk.mol.isValid("C(CC");
logInfo("C(CC            isValid = %d", tf_bad);
%[text] isValid = 0。閉じていない括弧は SMILES の構文エラーです。
%%
%[text] ## やってみよう 4: 5 分子グリッド（適応的 subplot レイアウト）

smiles_list = ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", ...
               "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", ...
               "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"];   % イブプロフェンを追加
names       = ["エタノール", "ベンゼン", "アスピリン", "カフェイン", "イブプロフェン"];

N    = numel(smiles_list);
ncol = ceil(sqrt(N));
nrow = ceil(N / ncol);      % N=5 のとき 2 行 x 3 列

figure("Name", "分子グリッド (5)", "Position", [100 100 1000 700]);
for i = 1:N
    subplot(nrow, ncol, i);
    mol = emk.mol.fromSmiles(smiles_list(i));
    emk.viz.draw2d(mol, Title=names(i));
end
%[text] ceil(sqrt(5)) = 3 列、ceil(5/3) = 2 行 → 2×3 グリッド（1 マス空き）

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
