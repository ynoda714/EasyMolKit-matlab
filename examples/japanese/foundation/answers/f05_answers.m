%[text] # F05 解答: 部分構造検索
%[text] f05_substructure_search.m の演習の参照解答。
addpath(genpath("src"));
emk.setup.initPython();
%%
%[text] ## 解答 E1: 日用化学品でフェノール類を見つける

data = readtable("data/list/everyday_chemicals.csv", TextType="string");

logInfo("日用化学品のフェノール類（芳香族環 + ヒドロキシル）:");
count = 0;
for i = 1:height(data)
    if ~emk.mol.isValid(data.SMILES(i)); continue; end
    mol = emk.mol.fromSmiles(data.SMILES(i));
    has_ring = emk.mol.hasSubstruct(mol, "c1ccccc1");
    has_oh   = emk.mol.hasSubstruct(mol, "[OH]");
    if has_ring && has_oh
        count = count + 1;
        logInfo("  %s", data.CommonName(i));
    end
end
logInfo("見つかったフェノール類の合計: %d 件", count);
%[text] フェノール類は消毒薬（チモール、オイゲノール）、香料、
%[text] 鎮痛薬（パラセタモール/アセトアミノフェン）に一般的。
%%
%[text] ## 解答 E2: FDA 薬物のカルボン酸数

fda = readtable("data/list/fda_drugs.csv", TextType="string");
cooh_count = 0;
logInfo("FDA 薬物 %d 件のカルボン酸をスキャン中...", height(fda));
for i = 1:height(fda)
    if ~emk.mol.isValid(fda.SMILES(i)); continue; end
    mol = emk.mol.fromSmiles(fda.SMILES(i));
    if emk.mol.hasSubstruct(mol, "C(=O)[OH]")
        cooh_count = cooh_count + 1;
    end
    logProgress(i, height(fda), "FDA 薬物");
end
frac = cooh_count / height(fda) * 100;
logInfo("FDA 薬物で -COOH 含有: %d / %d  (%.1f%%)", ...
    cooh_count, height(fda), frac);
%[text] NSAID（アスピリン、イブプロフェン、ナプロキセンなど）と多くの ACE 阻害薬が
%[text] カルボン酸を含む。全体の約 20〜30% が期待される。
%%
%[text] ## 解答 E3: 任意の分子の官能基サマリー

smiles_test = {"CN1CCC[C@H]1C2=CN=CC=C2", ...  % ニコチン
               "OC1=CC2=C(CC3N(CC23)CC4=CC=CC=C14)C=C1"};  % モルヒネ
mol_names   = {"ニコチン", "モルヒネ"};

fg_patterns = {"c1ccccc1", "[OH]", "C(=O)[OH]", "[NH,NH2]", "C(=O)[N]", ...
               "[F,Cl,Br,I]", "C(=O)[#6]"};
fg_names    = {"芳香族環","ヒドロキシル","カルボン酸","アミン","アミド", ...
               "ハロゲン","ケトン"};

for m = 1:numel(smiles_test)
    if ~emk.mol.isValid(smiles_test{m}); continue; end
    mol = emk.mol.fromSmiles(smiles_test{m});
    logInfo("--- %s ---", mol_names{m});
    for p = 1:numel(fg_patterns)
        has = emk.mol.hasSubstruct(mol, fg_patterns{p});
        logInfo("  %-14s : %s", fg_names{p}, string(has));
    end
end
%%
%[text] ## やってみよう 1: カフェインで SMARTS を試してみましょう

mol_caf = emk.mol.fromSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");

has_n   = emk.mol.hasSubstruct(mol_caf, "[#7]");
has_ket = emk.mol.hasSubstruct(mol_caf, "[C](=O)[#6]");
has_oh  = emk.mol.hasSubstruct(mol_caf, "[OH]");
has_ami = emk.mol.hasSubstruct(mol_caf, "[C](=O)[N]");

logInfo("カフェイン -- 任意の窒素 [#7]: %d", has_n);
logInfo("カフェイン -- ケトン [C](=O)[#6]: %d", has_ket);
logInfo("カフェイン -- ヒドロキシル [OH]: %d", has_oh);
logInfo("カフェイン -- アミド [C](=O)[N]: %d", has_ami);
%[text] (a) [#7] = true: カフェインは 4 個の窒素原子を持つ。
%[text] (b) [C](=O)[#6] = false: C=O 基はアミド/イミドカルボニルで、隣接するのは炭素ではなく窒素（N）。
%[text] (c) [OH] = false: カフェインにヒドロキシル基はない。
%[text] (d) [C](=O)[N] = true: アミドカルボニル（2 個）が検出される。
%%
%[text] ## やってみよう 2: ハロゲン列を追加してみましょう

data_hal = readtable("data/list/everyday_chemicals.csv", TextType="string");
fg_pat2 = {"c1ccccc1", "[OH]", "C(=O)[OH]", "[NH,NH2]", "[C](=O)[N]", "[F,Cl,Br,I]"};
fg_nm2  = ["ベンゼン環", "ヒドロキシル", "カルボン酸", "アミン", "アミド", "ハロゲン"];

n_hal = height(data_hal);
fg_mat2 = false(n_hal, numel(fg_pat2));
for i = 1:n_hal
    if ~emk.mol.isValid(data_hal.SMILES(i)); continue; end
    mol = emk.mol.fromSmiles(data_hal.SMILES(i));
    for p = 1:numel(fg_pat2)
        fg_mat2(i, p) = emk.mol.hasSubstruct(mol, fg_pat2{p});
    end
end
fg_tbl2 = array2table(fg_mat2, "VariableNames", fg_nm2);
fg_tbl2.Name = data_hal.CommonName;
fg_tbl2 = movevars(fg_tbl2, "Name", "Before", "ベンゼン環");

logInfo("ハロゲン含有日用化学品:");
halogen_rows = fg_tbl2(fg_tbl2.ハロゲン, ["Name", "ハロゲン"]);
disp(halogen_rows);
%[text] 日用化学品の多くはハロゲンを含まない。塩素含有の殺菌剤
%[text] （triclocarban など）が見つかる場合があります。

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
