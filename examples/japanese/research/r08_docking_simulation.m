%[text] # R08: ドッキングシミュレーション — SARS-CoV-2 Mpro vs 10 種の抗ウイルス薬
%[text] EasyMolKit Research — Layer 4（<u>**MATLAB Online 専用**</u>）
%[text] なぜ Nirmatrelvir（Paxlovid の有効成分）は SARS-CoV-2 の主要プロテアーゼ（Mpro）を強力に阻害し、HIV 薬は臨床で失敗したのでしょうか。
%[text] COVID-19 初期に期待された Lopinavir・Ritonavir は臨床試験で無効でしたが、HCV 薬の Boceprevir は Mpro との交叉反応性（本来の標的とは異なるタンパク質への結合）が文献で実証されています。
%[text] この違いは「分子の形」から説明できます。
%[text] 分子ドッキングシミュレーション（molecular docking）とは、リガンド（低分子化合物）をタンパク質の活性部位に仮想的に「はめ込み」、結合エネルギーを計算する計算化学ツールです。
%[text] このスクリプトでは 10 種の抗ウイルス薬を SARS-CoV-2 Mpro（PDB: 6LU7）にドッキングし、「意外な正解・意外な不正解」を数値と 3D 可視化で体験します。
%[text] ## 学習目標
%[text] - 分子ドッキングのパイプライン（タンパク質準備 → リガンド 3D 生成 → PDBQT 変換 → Vina 実行）を理解できます。
%[text] - 結合エネルギー（kcal/mol）の意味と計算上の限界を理解できます。
%[text] - なぜ標的外の薬が Mpro に「偶然」結合するのかを構造から考察できます。
%[text] - `pdbfixer` によるタンパク質構造前処理を体験できます。
%[text] - RDKit による 3D 構造生成（EmbedMolecule + UFF 最適化）を理解できます。
%[text] - MATLAB の `uihtml` と 3Dmol.js でインタラクティブな 3D リボン可視化を体験できます。 \
%[text] ## SARS-CoV-2 Mpro（3CLpro）の概要
%[text] - Mpro（主要プロテアーゼ）は SARS-CoV-2 の複製に必須のシステインプロテアーゼ（EC 3.4.22.69）です。
%[text] - 触媒残基は His41（一般塩基）と Cys145（求核剤）で、これら 2 残基が「触媒ダイアド」を形成します。
%[text] - 基質認識に関わる主要残基は Glu166（P1 の Gln を認識）、Met49 および Met165（疎水性 P2 ポケット）です。
%[text] - 使用する PDB 構造 6LU7 は N3 共有結合阻害剤との複合体です（Jin et al. 2020）。 \
%[text] ## 10 種のリガンド
%[text] **正解（Mpro 阻害剤）:**
%[text] 1. Nirmatrelvir — Paxlovid の有効成分。基準値
%[text] 2. Ensitrelvir — 日本発 Mpro 阻害剤（Shionogi）。設計が異なる \
%[text] **交叉反応（Mpro ≠ 本来の標的、だが結合は強い）:**
%[text] 1. Boceprevir — HCV NS3/4A プロテアーゼ阻害剤。文献で Mpro 交叉反応実証済み \
%[text] **臨床失敗薬（構造スコアと試験結果の乖離を学ぶ）:**
%[text] 1. Lopinavir — HIV プロテアーゼ阻害剤。初期 COVID 候補だが臨床無効
%[text] 2. Ritonavir — HIV PK ブースター。Paxlovid の成分だが単体では弱い
%[text] 3. Nelfinavir — HIV プロテアーゼ阻害剤。Mpro 親和性の限界 \
%[text] **標的が全く異なる（低スコアを予測）:**
%[text] 1. Remdesivir — RNA ポリメラーゼ阻害剤。プロテアーゼ非標的
%[text] 2. Favipiravir — RNA ポリメラーゼ阻害剤。軽量（MW=157）で活性部位に収まらない
%[text] 3. Oseltamivir — インフルエンザ ノイラミニダーゼ阻害剤
%[text] 4. Acyclovir — ヘルペスチミジンキナーゼ阻害剤（核酸アナログ） \
%[text] ## 前提条件
%[text] - A05（ニューラルネット）または R04（タンパク質-リガンド）修了を推奨
%[text] - Toolbox 不要
%[text] - MATLAB Online 必須（meeko、vina、pdbfixer は Linux バックエンドが必要） \
%[text] ## 動作環境
%[text] この演習は **MATLAB Online**（Linux バックエンド）専用です。
%[text] Desktop 非対応の理由: vina は Windows pip wheel がなく（Boost C++ ビルドが必要）、`pdbfixer` は Windows Smart App Control が `_openmm.pyd` をブロックします。
%[text] MATLAB Online では `pip install` で 3 つすべてを自動解決できます。
%[text] 
%[text] **所要時間**: 10〜30 分（MATLAB Online の CPU 性能と分子サイズに依存。Lopinavir/Ritonavir のような大分子は 1〜2 分/リガンド、小分子は数秒）
%[text] 
%[text] ## 使用データ
%[text] - タンパク質: PDB 6LU7（セクション 1 で RCSB からダウンロード）
%[text] - リガンド: PubChem REST API から 10 種の正規化 SMILES を取得（セクション 2） \
%[text] ## 簡略化の注記
%[text] - 剛体受容体（誘起適合・柔軟なサイドチェーンは非対応）
%[text] - Gasteiger 部分電荷（AM1-BCC ではない）
%[text] - exhaustiveness = 4（教育速度。論文品質には 8 以上を推奨）
%[text] - Nirmatrelvir: 非共有結合近似（Cys145 との共有結合は無視） \
%[text] ## 引用文献
%[text] - Jin Z et al. (2020) Structure of Mpro from SARS-CoV-2 and discovery of its inhibitors. Nature 582:289-293. doi:10.1038/s41586-020-2223-y \[6LU7 結晶構造と N3 阻害剤\]
%[text] - Owen DR et al. (2021) An oral SARS-CoV-2 Mpro inhibitor clinical candidate. Science 374:1586-1593. doi:10.1126/science.abl4784 \[Nirmatrelvir 発見論文\]
%[text] - Ma C et al. (2020) Boceprevir, GC-376, and calpain inhibitors II, XII inhibit SARS-CoV-2 viral replication. Cell Res 30:678-692. doi:10.1038/s41422-020-0356-z \[Boceprevir の Mpro 交叉反応性\]
%[text] - Cao B et al. (2020) A trial of lopinavir-ritonavir in adults hospitalized with severe Covid-19. N Engl J Med 382:1787-1799. doi:10.1056/NEJMoa2001282 \[Lopinavir-Ritonavir の COVID-19 臨床失敗\]
%[text] - Trott O & Olson AJ (2010) AutoDock Vina: improving the speed and accuracy of docking. J Comput Chem 31:455-461. doi:10.1002/jcc.21334 \[Vina スコアリング関数\]
%[text] - Eberhardt J et al. (2021) AutoDock Vina 1.2.0. J Chem Inf Model 61:3891-3898. doi:10.1021/acs.jcim.1c00203 \[Vina Python API\] \
%[text] 
%[text] Ctrl+Enter でセクションを 1 つずつ実行してください。
%%
%[text] ## セクション 0: セットアップ — Online ガードとライブラリインストール
% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif isfolder(fullfile(pwd, 'src'))
    addpath(genpath(fullfile(pwd, 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
%[text] MATLAB Online 専用ガード
if ~emk.util.isOnline()
    logWarn("R08: この演習は MATLAB Online 専用です。");
    logWarn("     meeko / vina / pdbfixer は Windows Desktop では動作しません。");
    logWarn("     理由: vina は Windows pip wheel なし（Boost C++ ビルド必須）、");
    logWarn("           pdbfixer は Windows Smart App Control が _openmm.pyd をブロック。");
    logWarn("     MATLAB Online (https://matlab.mathworks.com) でこのスクリプトを実行してください。");
    emk.setup.recipe("docking");
    return
end
%[text] ドッキングライブラリのフラグを有効化して installOnline に渡す
cfg = emkLoadConfig();
cfg.useCase.docking = true;   % scipy + meeko + vina + pdbfixer
emk.setup.installOnline(Config=cfg);
%[text] セクション 0a: チューニングパラメータ（ここで変更してから実行）
EXHAUSTIVENESS = 4;       % Vina 探索回数 (4=教育速度, 8=論文品質, 32=高精度)
N_POSES        = 3;       % 保存するドッキングポーズ数
%[text] 6LU7 の Mpro 活性部位バインディングボックス（N3 阻害剤重心から算出; Chain A）
%[text] 参考: Jin et al. 2020 / AutoDock Vina チュートリアル（広く引用される座標値）
BOX_CENTER = [-26.3, 12.7, 58.6];  % [x, y, z] Angstrom
BOX_SIZE   = [28.0, 28.0, 28.0];   % [dx, dy, dz] Angstrom
%[text] 成果物ディレクトリ
runDir = makeRunDir("Prefix", "r08_docking");
logInfo("R08: 成果物ディレクトリ -> %s", runDir);
%[text] ディレクトリ構造の準備
receptorDir = fullfile(runDir, "receptor");
ligandDir   = fullfile(runDir, "ligands");
poseDir     = fullfile(runDir, "poses");
mkdir(receptorDir);
mkdir(ligandDir);
mkdir(poseDir);
logInfo("R08: Section 0 完了");
%%
%[text] ## セクション 1: PDB 取得とタンパク質前処理 (pdbfixer)
PDB_ID        = "6LU7";
receptorPdb   = fullfile(receptorDir, "6lu7_raw.pdb");
preparedPdb   = fullfile(receptorDir, "6lu7_prepared.pdb");
receptorPdbqt = fullfile(receptorDir, "receptor.pdbqt");
%[text] \--- 1a: RCSB から 6LU7 をダウンロード ---
if ~isfile(receptorPdb)
    logInfo("R08: PDB %s を RCSB からダウンロード中...", PDB_ID);
    pdbUrl = sprintf("https://files.rcsb.org/download/%s.pdb", upper(PDB_ID));
    try
        websave(receptorPdb, pdbUrl);
        logInfo("R08: PDB ダウンロード完了 -> %s", receptorPdb);
    catch ME
        error("emk:r08:pdbDownloadFailed", ...
            "PDB %s のダウンロードに失敗しました: %s\n" + ...
            "インターネット接続を確認してください。", PDB_ID, ME.message);
    end
else
    logInfo("R08: キャッシュ済み PDB を使用 -> %s", receptorPdb);
end
%[text] \--- 1b: pdbfixer でタンパク質を準備 ---
%[text] 
%[text] ### タンパク質前処理の必要性
%[text] 結晶構造（PDB ファイル）は実験上の制約から、そのままではドッキングに使えない不完全さを持っています。
%[text] 主な問題点: 欠損残基（電子密度が弱い領域）の存在、水素原子の欠如（X 線は電子密度を測定するため）、共結晶リガンドや結晶水の混在です。
%[text] `pdbfixer` は欠損残基・原子の補完、ヘテロ原子の除去、水素付加を自動で行い、「タンパク質のみ・水素付加済み」の状態を準備します。
logInfo("R08: pdbfixer でタンパク質を準備中...");
try
    pdbfixer = py.importlib.import_module("pdbfixer");

    % PDBFixer オブジェクトを作成
    fixer = pdbfixer.PDBFixer(filename=receptorPdb);
    logInfo("R08: PDB 読み込み完了");

    % 欠損残基・欠損重原子を補完
    fixer.findMissingResidues();
    fixer.findMissingAtoms();
    fixer.addMissingAtoms();
    logInfo("R08: 欠損原子補完完了");

    % ヘテロ原子（リガンド・水）を除去してレセプターのみにする
    % keepIds: 保持するレジデュー（空 = 全ヘテロ原子削除）
    fixer.removeHeterogens(false);
    logInfo("R08: ヘテロ原子（リガンド・結晶水）除去完了");

    % 水素原子を pH 7.4 でプロトン化
    fixer.addMissingHydrogens(py.float(7.4));
    logInfo("R08: 水素付加完了（pH 7.4）");

    % 準備済みタンパク質を PDB 形式で保存
    openmm = py.importlib.import_module("openmm.app");
    fh     = py.open(preparedPdb, "w");
    openmm.PDBFile.writeFile(fixer.topology, fixer.positions, fh);
    fh.close();
    logInfo("R08: 準備済みレセプター保存 -> %s", preparedPdb);

catch ME
    error("emk:r08:pdbfixerFailed", ...
        "pdbfixer によるタンパク質前処理に失敗しました: %s\n" + ...
        "pdbfixer がインストールされているか確認してください: emk.setup.validate()", ...
        ME.message);
end
%[text] \--- 1c: meeko mk\_prepare\_receptor でレセプター PDBQT を生成 ---
%[text] 
%[text] ### PDBQT 形式とは
%[text] AutoDock Vina は PDBQT（PDB + 電荷 + 原子型）形式を入力として使用します。
%[text] `meeko` の `mk_prepare_receptor` は Gasteiger 電荷（原子の電気陰性度から計算する半経験的な部分電荷）を自動計算して変換します。
logInfo("R08: レセプター PDBQT 生成中...");
try
    % PATH に ~/.local/bin を追加（pip --user インストール時のスクリプト位置）
    % meeko 0.7.x の API 変更に注意:
    %   旧: mk_prepare_receptor -i <pdb>  -> meeko 0.7 で --read_with_prody に変更
    %   新: mk_prepare_receptor --read_pdb <pdb> -p <pdbqt>  (prody 不要)
    prepareCmd = sprintf( ...
        "export PATH=$HOME/.local/bin:$PATH && mk_prepare_receptor --read_pdb '%s' -p '%s'", ...
        preparedPdb, receptorPdbqt);
    [status, cmdOut] = system(prepareCmd);
    if status ~= 0
        % フォールバック: python3 -c で meeko.cli モジュールを直接呼び出す
        % （シェルスクリプトの有無・PATH 設定に依存しないため堅牢）
        logWarn("R08: mk_prepare_receptor 失敗 -> Python モジュール実行で再試行...");
        pySys = py.importlib.import_module("sys");
        pySub = py.importlib.import_module("subprocess");
        % パスは sys.argv 経由で渡す（シングルクォートや空白を含むパスに安全）
        % python3 -c "code" <pdb> <pdbqt> のとき sys.argv = ['-c', pdb, pdbqt]
        % strjoin で 1×N char に変換（string 配列のまま char() すると 2D 行列になる）
        pyInlineCode = char(strjoin([ ...
            "import sys;", ...
            "sys.argv=['mk_prepare_receptor','--read_pdb',sys.argv[1],'-p',sys.argv[2]];", ...
            "from meeko.cli.mk_prepare_receptor import main;main()"], ""));
        result = pySub.run( ...
            py.list({char(pySys.executable), '-c', pyInlineCode, ...
                     char(preparedPdb), char(receptorPdbqt)}), ...
            capture_output=py.True, text=py.True);
        if int32(result.returncode) ~= 0
            error("emk:r08:receptorPdbqtFailed", ...
                "mk_prepare_receptor 失敗:\n%s\n%s", ...
                char(result.stdout), char(result.stderr));
        end
    end
    logInfo("R08: レセプター PDBQT 生成完了 -> %s", receptorPdbqt);
catch ME
    error("emk:r08:receptorPdbqtFailed", ...
        "レセプター PDBQT 生成に失敗しました: %s\n" + ...
        "meeko がインストールされているか確認してください: emk.setup.validate()", ...
        ME.message);
end

logInfo("R08: Section 1 完了");
%%
%[text] ## セクション 2: リガンド 10 種の SMILES 取得と 3D 構造生成
%[text] ### リガンド 3D 構造生成パイプライン
%[text] ドッキングにはリガンドの 3D 座標が必要です。SMILES（2D トポロジー）から始まり、次の手順で 3D 構造を生成します。
%[text] 1. `AddHs` で水素を付加した Mol オブジェクトを生成します。
%[text] 2. `EmbedMolecule`（ETKDGv3 距離ジオメトリ法）で初期 3D 座標を計算します。
%[text] 3. `UFFOptimizeMolecule`（Universal Force Field）で力場最適化を行い、歪みのない構造に収束させます。
%[text] 4. `meeko MoleculePreparation` で Gasteiger 電荷を付与し、PDBQT 形式（Vina 入力）に変換します。 \
%[text] 
%[text] \--- 2a: リガンド定義 (PubChem CID + SMILES 注記) ---
%[text] SMILES ソース: PubChem REST API (webread 経由で正規化済み)
%[text] CID は PubChem 公式 ID（不変）。SMILES はセクション 2b で取得します。
ligandDefs = struct( ...
    "name",    {"Nirmatrelvir",   "Ensitrelvir",  "Boceprevir", ...
                "Lopinavir",      "Ritonavir",    "Nelfinavir", ...
                "Remdesivir",     "Favipiravir",  "Oseltamivir", ...
                "Acyclovir"}, ...
    "pubchemQuery", ...
               {"Nirmatrelvir",   "Ensitrelvir",  "Boceprevir", ...
                "Lopinavir",      "Ritonavir",    "Nelfinavir", ...
                "Remdesivir",     "Favipiravir",  "Oseltamivir", ...
                "Acyclovir"}, ...
    "target",  {"SARS-CoV-2 Mpro",  "SARS-CoV-2 Mpro", "HCV NS3/4A protease", ...
                "HIV protease",      "HIV protease",     "HIV protease", ...
                "RNA polymerase",    "RNA polymerase",   "Neuraminidase", ...
                "Thymidine kinase"}, ...
    "expectation", ...
               {"HIGH",   "HIGH",  "MODERATE-HIGH", ...
                "LOW-MOD", "LOW",  "LOW", ...
                "VERY-LOW", "VERY-LOW", "VERY-LOW", ...
                "VERY-LOW"} ...
);
nLigands = numel(ligandDefs);
logInfo("R08: %d 種のリガンドを処理します", nLigands);
%[text] \--- 2b: PubChem REST API で canonical SMILES を取得 ---
logInfo("R08: PubChem から canonical SMILES を取得中...");
ligandSmiles = strings(1, nLigands);
for k = 1:nLigands
    query = ligandDefs(k).pubchemQuery;
    try
        % PubChem PUG REST: 名前 → canonical SMILES
        apiUrl = sprintf( ...
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/CanonicalSMILES/TXT", ...
            urlencode(query));
        smi = strtrim(webread(apiUrl));
        % webread が複数行返した場合は先頭行のみ使用
        lines = strsplit(smi, newline);
        ligandSmiles(k) = strtrim(lines{1});
        logInfo("  [%d/%d] %-16s -> %s", k, nLigands, query, ...
            extractBefore(ligandSmiles(k) + " ", min(50, strlength(ligandSmiles(k))+1)));
    catch ME
        logWarn("  [%d/%d] %-16s -> PubChem 取得失敗: %s", k, nLigands, query, ME.message);
    end
end
%[text] \--- 2c: SMILES → 3D 構造生成 (RDKit EmbedMolecule + UFF) ---
logInfo("R08: RDKit で 3D 座標を生成中...");
%[text] 3D 生成を subprocess で実行（セクション 1 と同じフォールバックパターン）
%[text] 理由: RDKit の AddHs は RWMol（read-write mol）を返します。RWMol の Boost.Python メタクラスを MATLAB が `py.Boost.Python.class` に型マップしようとして失敗するため、Python プロセス内で完結させる方が堅牢です。
pySys = py.importlib.import_module("sys");
pySub = py.importlib.import_module("subprocess");
pythonExe = char(pySys.executable);
%[text] Python インライン: argv\[1\]=SMILES, argv\[2\]=出力SDF パス
%[text] ETKDGv3 + 固定シード + UFF 最適化 → 水素除去 SDF
pyGen3dCode = char(strjoin([ ...
    "from rdkit import Chem; from rdkit.Chem import AllChem; import sys;", ...
    "smi,out=sys.argv[1],sys.argv[2];", ...
    "mol=Chem.MolFromSmiles(smi);", ...
    "assert mol is not None,'invalid SMILES: '+smi;", ...
    "mol=Chem.AddHs(mol);", ...
    "p=AllChem.ETKDGv3(); p.randomSeed=42;", ...
    "r=AllChem.EmbedMolecule(mol,p);", ...
    "assert r>=0,'EmbedMolecule failed (r='+str(r)+')';", ...
    "AllChem.UFFOptimizeMolecule(mol,maxIters=500);", ...
    "mol=Chem.RemoveHs(mol);", ...
    "w=Chem.SDWriter(out); w.write(mol); w.close()"], " "));

ligandSdfFiles = strings(1, nLigands);
ligandValid    = false(1, nLigands);

for k = 1:nLigands
    name = ligandDefs(k).name;
    smi  = ligandSmiles(k);
    if smi == ""
        logWarn("  [%02d] %s: SMILES なし -- スキップ", k, name);
        continue;
    end

    try
        sdfPath = fullfile(ligandDir, sprintf("%02d_%s.sdf", k, name));

        % Python サブプロセスで 3D 構造生成
        result = pySub.run( ...
            py.list({pythonExe, '-c', pyGen3dCode, char(smi), char(sdfPath)}), ...
            capture_output=py.True, text=py.True);

        if int32(result.returncode) ~= 0
            error("emk:r08:embedFailed", "3D 生成失敗:\n%s", char(result.stderr));
        end

        ligandSdfFiles(k) = sdfPath;
        ligandValid(k)    = true;

        mol2d = emk.mol.fromSmiles(smi);   % MW 計算用
        mw = emk.descriptor.molWeight(mol2d);
        logInfo("  [%02d] %-15s: 3D 生成完了 (MW=%.1f, %s)", ...
            k, name, mw, sdfPath);

    catch ME
        logWarn("  [%02d] %s: エラー -> %s", k, name, ME.message);
    end
end

nValid = sum(ligandValid);
logInfo("R08: 3D 構造生成完了 -- %d / %d 成功", nValid, nLigands);
if nValid == 0
    error("emk:r08:noValidLigands", ...
        "有効なリガンド 3D 構造が 1 つも生成できませんでした。Section 0 を確認してください。");
end
logInfo("R08: Section 2 完了");
%%
%[text] ## セクション 3: PDBQT 変換 (meeko MoleculePreparation)
%[text] ### PDBQT 形式の必要性
%[text] AutoDock Vina は各原子に「AutoDock 原子型」と「部分電荷」を要求します。
%[text] `meeko` はこれを RDKit 原子型から自動マッピングし、Gasteiger 電荷を計算して PDBQT ファイルを生成します。
%[text] LGPL-2.1 ライブラリであるため、Python の動的リンク利用の範囲で EasyMolKit への copyleft 汚染はありません（docs/compliance.md CL-7 参照）。
%[text] 
%[text] PDBQT 変換も subprocess で実行します（セクション 2 と同じパターン）。
%[text] 理由: `SDMolSupplier` の添字アクセス（{}）が RWMol（Boost.Python）を返し MATLAB が失敗します。また meeko 0.7.x では `MoleculePreparation.prepare()` の戻り値の型が変わる API 変更があるため、Python プロセス内で完結させる方が堅牢です。
pySys2 = py.importlib.import_module("sys");
pySub2 = py.importlib.import_module("subprocess");
pythonExe2 = char(pySys2.executable);
%[text] Python インライン: argv\[1\]=SDF パス, argv\[2\]=PDBQT 出力パス
%[text] meeko 0.5.x / 0.7.x 両対応 (prepare / write\_string の戻り値がタプルかどうかを確認)
pyMeekoCode = char(strjoin([ ...
    "import sys; from rdkit import Chem; from rdkit.Chem import SDMolSupplier;", ...
    "from meeko import MoleculePreparation,PDBQTWriterLegacy;", ...
    "sdf,out=sys.argv[1],sys.argv[2];", ...
    "suppl=SDMolSupplier(sdf,removeHs=False);", ...
    "mol=next(iter(suppl));", ...
    "assert mol is not None,'SDF read failed: '+sdf;", ...
    "mol=Chem.AddHs(mol,addCoords=True);", ...
    "prep=MoleculePreparation();", ...
    "res=prep.prepare(mol);", ...
    "setup=res[0] if isinstance(res,list) else list(prep.setup_dict.values())[0];", ...
    "assert setup is not None,'meeko setup empty';", ...
    "wr=PDBQTWriterLegacy.write_string(setup);", ...
    "pdbqt=wr[0] if isinstance(wr,tuple) else wr;", ...
    "open(out,'w').write(pdbqt)"], " "));

ligandPdbqtFiles = strings(1, nLigands);

logInfo("R08: meeko で PDBQT 変換中...");
for k = 1:nLigands
    if ~ligandValid(k)
        continue;
    end
    name    = ligandDefs(k).name;
    sdfPath = ligandSdfFiles(k);

    try
        pdbqtPath = fullfile(ligandDir, sprintf("%02d_%s.pdbqt", k, name));

        result2 = pySub2.run( ...
            py.list({pythonExe2, '-c', pyMeekoCode, char(sdfPath), char(pdbqtPath)}), ...
            capture_output=py.True, text=py.True);

        if int32(result2.returncode) ~= 0
            error("emk:r08:meekoPdbqtFailed", "meeko 失敗:\n%s", char(result2.stderr));
        end

        ligandPdbqtFiles(k) = pdbqtPath;
        logInfo("  [%02d] %-15s: PDBQT 生成完了", k, name);

    catch ME
        logWarn("  [%02d] %s: meeko エラー -> %s", k, name, ME.message);
        ligandValid(k) = false;
    end
end

nValid = sum(ligandValid);
logInfo("R08: PDBQT 変換完了 -- %d / %d 成功", nValid, nLigands);
logInfo("R08: Section 3 完了");
%%
%[text] ## セクション 4: ドッキング実行 (AutoDock Vina Python API)
%[text] ### AutoDock Vina のスコアリング関数
%[text] Vina は反復局所探索（Iterated Local Search）でリガンドの配座空間を探索し、経験的スコアリング関数 $\\Delta G \\approx w(1) \\cdot \\mathrm{gauss1} + w(2) \\cdot \\mathrm{gauss2} + w(3) \\cdot \\mathrm{repulsion} + w(4) \\cdot \\mathrm{hydrophobic} + w(5) \\cdot \\mathrm{HBond}$ を最小化します（Trott & Olson 2010）。
%[text] スコアの単位は kcal/mol で、より負の値ほど強い結合親和性を示します。
%[text] 一般的な目安: $-7$ kcal/mol 以下は有望なリード化合物、$-9$ kcal/mol 以下は強い結合とみなされます。
%[text] 
%[text] この演習での主な簡略化:
%[text] - Nirmatrelvir は Cys145 と共有結合しますが、ここでは非共有結合近似を使用します。
%[text] - 剛体受容体近似のため誘起適合効果は無視されます。 \
vina = py.importlib.import_module("vina");

bindingEnergies  = nan(1, nLigands);    % best pose [kcal/mol]
dockingPoseFiles = strings(1, nLigands);

logInfo("R08: ドッキング開始 (exhaustiveness=%d)...", EXHAUSTIVENESS);
logInfo("R08: 推定所要時間 %d ~ %d 分", nValid * 1, nValid * 3);

for k = 1:nLigands
    if ~ligandValid(k)
        continue;
    end
    name        = ligandDefs(k).name;
    pdbqtPath   = ligandPdbqtFiles(k);
    poseOutPath = fullfile(poseDir, sprintf("%02d_%s_poses.pdbqt", k, name));

    logInfo("  [%02d/%02d] %s -- ドッキング中...", k, nLigands, name);
    tic;

    try
        % Vina インスタンス生成（CPU スレッド数は自動設定）
        v = vina.Vina(sf_name="vina", verbosity=py.int(0));

        % レセプターと配位子をセット
        v.set_receptor(receptorPdbqt);
        v.set_ligand_from_file(pdbqtPath);

        % バインディングボックス設定
        % num2cell では Python がスカラーに変換できないため py.float() で明示変換
        vinaCenter  = py.list({py.float(BOX_CENTER(1)), py.float(BOX_CENTER(2)), py.float(BOX_CENTER(3))});
        vinaBoxSize = py.list({py.float(BOX_SIZE(1)),   py.float(BOX_SIZE(2)),   py.float(BOX_SIZE(3))});
        v.compute_vina_maps(center=vinaCenter, box_size=vinaBoxSize);

        % ドッキング実行
        v.dock(exhaustiveness=py.int(EXHAUSTIVENESS), n_poses=py.int(N_POSES));

        % 結果取得
        % v.energies() は numpy ndarray (n_poses x n_terms) を返す
        % flatten() で 1D にしてから MATLAB double に変換
        pyEn     = v.energies(n_poses=py.int(N_POSES));
        energies = reshape(double(pyEn.flatten()), N_POSES, []);
        bindingEnergies(k) = energies(1, 1);   % best pose total score [kcal/mol]

        % ポーズを PDBQT で保存
        v.write_poses(poseOutPath, n_poses=py.int(N_POSES), overwrite=true);
        dockingPoseFiles(k) = poseOutPath;

        elapsed = toc;
        logInfo("    -> スコア: %.2f kcal/mol (%.1f 秒)", bindingEnergies(k), elapsed);

    catch ME
        toc;
        logWarn("  [%02d] %s: ドッキングエラー -> %s", k, name, ME.message);
    end
end

logInfo("R08: 全ドッキング完了");
logInfo("R08: Section 4 完了");
%%
%[text] ## セクション 5: 結合エネルギー比較と考察
%[text] \--- 5a: 結果テーブル作成 ---
names_all     = string({ligandDefs.name})';
targets_all   = string({ligandDefs.target})';
expect_all    = string({ligandDefs.expectation})';
energies_col  = bindingEnergies';
valid_col     = ligandValid';

resultTbl = table(names_all, targets_all, expect_all, energies_col, valid_col, ...
    VariableNames=["Compound", "OriginalTarget", "Expected", "DeltaG_kcal_mol", "Success"]);
%[text] スコア順（最も負の値 = 最強結合）でソート
validRows = resultTbl.Success;
scored    = resultTbl(validRows, :);
scored    = sortrows(scored, "DeltaG_kcal_mol", "ascend");

%[text] ### ドッキング結果サマリー（スコア順）
logInfo("%-20s  %-25s  %-8s  %s", "Compound", "Original Target", "Expected", "Score (kcal/mol)");
logInfo("%s", repmat("-", 1, 75));
for i = 1:height(scored)
    logInfo("%-20s  %-25s  %-8s  %.2f", ...
        scored.Compound(i), scored.OriginalTarget(i), ...
        scored.Expected(i), scored.DeltaG_kcal_mol(i));
end
%[text] \--- 5b: CSV エクスポート ---
csvPath = fullfile(runDir, "r08_docking_results.csv");
writetable(resultTbl, csvPath);
logInfo("R08: 結果 CSV -> %s", csvPath);
%[text] \--- 5c: 水平棒グラフ (barh) ---
%[text] 
%[text] ### 共有結合阻害剤のスコア過小評価に注意
%[text] このシミュレーションでは Nirmatrelvir と Boceprevir が予想より低スコアになります。
%[text] 理由: Vina は非共有結合スコアリングのため、共有結合の寄与（+2〜4 kcal/mol）を計上できません。
%[text] - Nirmatrelvir: ニトリル基が Cys145-SH と Michael 付加（covalent warhead）
%[text] - Boceprevir: ketoamide が Cys145-SH と求核付加（Ma et al. 2020, Mpro Ki ~9 µM） \
%[text] Lopinavir 等の HIV 薬が高スコアになっても、それは「形状の偶然一致」であり実際の Mpro 活性（臨床無効: Cao et al. 2020）とは無関係です。
fig1 = figure("Name", "R08: Docking Scores vs Mpro (6LU7)", "NumberTitle", "off");
set(fig1, "Position", [100, 100, 800, 500]);

cmap = [0.85, 0.33, 0.10;    % orange-red  = Mpro inhibitors
        0.85, 0.33, 0.10;
        0.93, 0.69, 0.13;    % yellow      = cross-reactive
        0.47, 0.67, 0.19;    % green       = off-target (HIV)
        0.47, 0.67, 0.19;
        0.47, 0.67, 0.19;
        0.30, 0.75, 0.93;    % blue        = RNA pol / other
        0.30, 0.75, 0.93;
        0.30, 0.75, 0.93;
        0.30, 0.75, 0.93];

nScored = height(scored);
scoreVals = scored.DeltaG_kcal_mol;
compNames = scored.Compound;
%[text] 色を順位に合わせて並べ替え
barColors = zeros(nScored, 3);
for i = 1:nScored
    idx = find(strcmp(names_all, compNames(i)), 1);
    if ~isempty(idx) && idx <= size(cmap, 1)
        barColors(i, :) = cmap(idx, :);
    else
        barColors(i, :) = [0.5, 0.5, 0.5];
    end
end

bh = barh(scoreVals);
bh.FaceColor = "flat";
for i = 1:nScored
    bh.CData(i, :) = barColors(i, :);
end

yticklabels(compNames);
xlabel("Binding Free Energy (kcal/mol)", "FontSize", 12);
title(sprintf("AutoDock Vina Scores vs SARS-CoV-2 Mpro (PDB: 6LU7)\nexhaustiveness=%d", ...
    EXHAUSTIVENESS), "FontSize", 12);
grid on;
ax = gca;
ax.GridAlpha = 0.3;
%[text] 一般目安ラインを追加
xline(-7.0, "--", "Promising lead (-7)", "Color", [0.5, 0, 0.5], ...
    "LabelHorizontalAlignment", "left", "FontSize", 9);
xline(-9.0, "--", "Strong binding (-9)", "Color", [0.8, 0, 0], ...
    "LabelHorizontalAlignment", "left", "FontSize", 9);
%[text] 数値ラベルを各バーに追加
for i = 1:nScored
    text(scoreVals(i) - 0.2, i, sprintf("%.1f", scoreVals(i)), ...
        "HorizontalAlignment", "right", "FontSize", 9, "Color", "white", "FontWeight", "bold");
end

%[text] 凡例を追加（色とカテゴリの対応）
hold on;
hLeg = [ ...
    patch(NaN, NaN, [0.85, 0.33, 0.10]); ...
    patch(NaN, NaN, [0.93, 0.69, 0.13]); ...
    patch(NaN, NaN, [0.47, 0.67, 0.19]); ...
    patch(NaN, NaN, [0.30, 0.75, 0.93])];
legend(hLeg, ["Mpro inhibitors", "Cross-reactive (HCV)", "HIV protease", "RNA pol / other"], ...
    "Location", "southeast", "FontSize", 9);
hold off;

plotPath = fullfile(runDir, "r08_docking_scores.png");
saveas(fig1, plotPath);
logInfo("R08: スコア比較プロット -> %s", plotPath);
%[text] \--- 5d: 考察 ---
if height(scored) == 0
    logWarn("R08: ドッキング結果がありません。Section 3-4 のエラーを確認してください。");
else
    bestCompound  = scored.Compound(1);
    bestScore     = scored.DeltaG_kcal_mol(1);
    worstCompound = scored.Compound(end);
    worstScore    = scored.DeltaG_kcal_mol(end);

    logInfo("R08: 考察");
    logInfo("  最高スコア  : %s (%.2f kcal/mol)", bestCompound, bestScore);
    logInfo("  最低スコア  : %s (%.2f kcal/mol)", worstCompound, worstScore);

    % Nirmatrelvir と Boceprevir のスコアを比較
    nirmIdx = strcmp(scored.Compound, "Nirmatrelvir");
    bocepIdx = strcmp(scored.Compound, "Boceprevir");
    if any(nirmIdx) && any(bocepIdx)
        nirmScore  = scored.DeltaG_kcal_mol(nirmIdx);
        bocepScore = scored.DeltaG_kcal_mol(bocepIdx);
        logInfo("  Nirmatrelvir vs Boceprevir: %.2f vs %.2f kcal/mol", nirmScore, bocepScore);
        logInfo("  【注記】両化合物は Cys145 と共有結合 -- Vina 非共有結合近似で過小評価されます。");
        if bocepScore < nirmScore
            logInfo("  -> Boceprevir が Nirmatrelvir より強いスコア: ketoamide 交叉反応性と一致。");
        else
            logInfo("  -> Boceprevir スコアが Nirmatrelvir 以下: 非共有結合近似による過小評価。");
            logInfo("     実験値: Boceprevir Mpro Ki ~9 uM (Ma et al. 2020) で交叉反応性は実証済み。");
        end
    end

    % Ritonavir が Nirmatrelvir より高スコアの場合の補足
    ritoIdx = strcmp(scored.Compound, "Ritonavir");
    if any(ritoIdx) && any(nirmIdx)
        ritoScore = scored.DeltaG_kcal_mol(ritoIdx);
        if ritoScore < nirmScore
            logInfo("  Ritonavir (%.2f) > Nirmatrelvir (%.2f): HIV 薬が Mpro 阻害剤より高スコア。", ...
                ritoScore, nirmScore);
            logInfo("  -> 形状偶然一致 + exhaustiveness=4 の乱択性による結果。");
            logInfo("     Ritonavir はプロテアーゼ阻害剤族で活性部位形状が類似するが、");
            logInfo("     単体 COVID 臨床試験では効果なし（Cao et al. 2020）。");
            logInfo("     スコアと臨床活性が乖離する典型例です。");
        end
    end

    % Favipiravir の MW とスコアを確認
    favIdx = strcmp(scored.Compound, "Favipiravir");
    if any(favIdx)
        favScore = scored.DeltaG_kcal_mol(favIdx);
        logInfo("  Favipiravir (MW=157): %.2f kcal/mol -- %s", favScore, ...
            ternary_(favScore > -5.5, "予想通り低スコア（ポケットに収まらない）", ...
                "想定以上のスコア -- 小分子でも活性部位の浅い領域に部分収容された可能性（exhaustiveness=4 の乱択性も影響）"));
    end
end
logInfo("R08: Section 5 完了");
%%
%[text] ## セクション 6: 3D 可視化 -- 活性部位と結合ポーズ
%[text] ### 3D 可視化の目的
%[text] 数値スコアだけでは「なぜ」結合するかがわかりません。
%[text] 3D 可視化を使うと、活性部位残基（His41、Cys145、Glu166）とリガンドの空間的関係、疎水性ポケット（P2 サイト: Met49、Met165）への充填状況、ドッキング化合物がポケット深部に収まっているか（pocket burial）を直感的に把握できます。
%[text] 
%[text] **ビューワー操作（インタラクティブ）:** 左ドラッグで回転、スクロールでズーム、右ドラッグで平行移動できます。
%[text] 初期ビューはドッキング化合物を中心にズームしています。
%[text] タンパク質全体を見るにはスクロールアウトしてください。
%[text] 
%[text] ### 可視化の解釈に関する注意（近似計算）
%[text] このドッキング結果は計算上のモデルであり、実際の結合構造（結晶構造）とは異なる場合があります。
%[text] **剛体受容体**: タンパク質構造を固定しています。実際の結合では誘起適合（induced fit）によりタンパク質側の残基も動きます。
%[text] **スコア誤差**: Vina スコアは結合自由エネルギーの近似値で、同一手法内での相対比較には有用ですが、絶対精度は ±2 kcal/mol 程度です。
%[text] **exhaustiveness=4**: 教育速度設定のため、コンフォーメーション空間を完全には探索しません（論文品質には 8 以上を推奨）。
%[text] 
%[text] **表示要素の読み方:**
%[text] - リボン（金色 = β シート / マゼンタ = α ヘリックス / 灰色 = ループ）は Mpro 主鎖の二次構造を示します。金色の β バレル（上部・Domain I/II）が活性部位を形成し、マゼンタのヘリックス束（Domain III）が調節ドメインです。
%[text] - 細いスティック群（8 群）は活性部位アミノ酸 8 残基のサイドチェーンを示します。各残基（His41・Met49・Cys145・His163・Met165・Glu166・Asp187・Gln189）が独立した細線クラスタとして見えるのは正常です。
%[text] - 元素色コード（CPK 配色）: 白/灰=炭素、赤=酸素、青=窒素、黄=硫黄（Met・Cys のみ）
%[text] - 太いスティック（中央）はドッキング化合物（最良ポーズのみ）を示します。ポケット中心に収まっているか確認してください。
%[text] - PDBQT ファイルには N\_POSES 個のポーズが格納されています。可視化は最良スコアのポーズ（MODEL 1）のみを表示します。全ポーズを重ねると同一分子が分裂したように見えるため、最初の ENDMDL 以降の原子は除外しています。 \
%[text] 
%[text] 本演習では定性的な傾向の把握を目的としています。定量的予測や実際の医薬品開発には、IC50 測定・結晶構造解析などの実験的検証が必要です。
%[text] 
%[text] \--- 6a: 準備済みレセプター PDB から Cα 座標を取得 ---
%[text] PDB ファイルを行単位でパースします（Base MATLAB のみ、Bioinformatics Toolbox 不要）。
logInfo("R08: PDB ファイルから Cα 座標をパース中...");

caCoords = [];   % Cα 座標 [N x 3]
caResIDs  = [];  % 残基番号

fh = fopen(preparedPdb, "r");
if fh < 0
    logWarn("R08: %s が見つからない -- 可視化をスキップ", preparedPdb);
else
    while ~feof(fh)
        line = fgetl(fh);
        if ~ischar(line) || length(line) < 54
            continue;
        end
        recType  = strtrim(line(1:6));
        atomName = strtrim(line(13:16));
        if (strcmp(recType, "ATOM") || strcmp(recType, "HETATM")) && ...
                strcmp(atomName, "CA")
            x = str2double(line(31:38));
            y = str2double(line(39:46));
            z = str2double(line(47:54));
            resNum = str2double(strtrim(line(23:26)));
            if ~isnan(x) && ~isnan(y) && ~isnan(z)
                caCoords = [caCoords; x, y, z];  %#ok<AGROW>
                caResIDs  = [caResIDs;  resNum];  %#ok<AGROW>
            end
        end
    end
    fclose(fh);
    logInfo("R08: %d 個の Cα 原子を読み込み", size(caCoords, 1));
end
%[text] \--- 6b: 活性部位残基を定義 ---
activeSiteResNums = [41, 145, 163, 164, 165, 166, 168, 172, 187, 189];
activeSiteLabels  = ["His41", "Cys145", "Pro163", "Leu164", "Met165", ...
                     "Glu166", "Ala168", "Gln172", "Asp187", "Gln189"];
activeCoords = [];
activeLabels = strings(0);

for i = 1:numel(activeSiteResNums)
    idx = find(caResIDs == activeSiteResNums(i), 1);
    if ~isempty(idx)
        activeCoords = [activeCoords; caCoords(idx, :)];  %#ok<AGROW>
        activeLabels = [activeLabels; activeSiteLabels(i)];  %#ok<AGROW>
    end
end
logInfo("R08: %d / %d 個の活性部位残基を検出", size(activeCoords, 1), numel(activeSiteResNums));
%[text] \--- 6c: ベストスコアのリガンドを 1 つ選択（ribbon 可視化対象）---
displayLigandName = "";
validPoseIdx = find(dockingPoseFiles ~= "");
if ~isempty(validPoseIdx)
    [~, bestIdx]     = min(bindingEnergies(validPoseIdx));
    displayLigandIdx  = validPoseIdx(bestIdx);
    displayLigandName = ligandDefs(displayLigandIdx).name;
    displayPoseFile   = dockingPoseFiles(displayLigandIdx);
    logInfo("R08: ribbon 表示対象: %s (best score)", displayLigandName);
end
%[text] \--- 6d: 3Dmol.js ribbon 可視化（uifigure 内で表示） ---
%[text] 
%[text] ### 自己完結 HTML 生成
%[text] `py3Dmol` の `write_html()` は外部 CDN から 3Dmol.js を読み込むため、MATLAB Online のサンドボックスブラウザではロードに失敗します。
%[text] そのため `urllib.request` で 3Dmol.js を直接ダウンロードして HTML に埋め込みます。
logInfo("R08: 3Dmol.js ribbon 可視化（自己完結 HTML）を生成中...");

pySub6     = py.importlib.import_module("subprocess");
pySys6     = py.importlib.import_module("sys");
pythonExe6 = char(pySys6.executable);

htmlPath  = fullfile(runDir, "r08_ribbon.html");
vizScript = fullfile(runDir, "viz_ribbon.py");
%[text] Python スクリプト生成
%[text] stdlib のみ使用（urllib.request + json）で外部パッケージは不要です。
%[text] 3Dmol.js を CDN からダウンロードして HTML に inline 埋め込みます。
%[text] PDB データは json.dumps() で安全にエスケープして JS 変数に渡します。
fid = fopen(char(vizScript), 'w');
fprintf(fid, "import sys, urllib.request, json\n");
fprintf(fid, "rec, lig, out = sys.argv[1], sys.argv[2], sys.argv[3]\n");
fprintf(fid, "rec_pdb = open(rec).read()\n");
fprintf(fid, "lig_pdb = ''\n");
fprintf(fid, "if lig != 'none':\n");
fprintf(fid, "    _ll = []\n");
fprintf(fid, "    for _l in open(lig):\n");
fprintf(fid, "        if _l.startswith('ENDMDL'): break\n");
fprintf(fid, "        if _l[:6] in ('ATOM  ','HETATM'): _ll.append(_l)\n");
fprintf(fid, "    if not _ll:  # fallback: no MODEL/ENDMDL tags\n");
fprintf(fid, "        _ll = [_l for _l in open(lig) if _l[:6] in ('ATOM  ','HETATM')]\n");
fprintf(fid, "    lig_pdb = ''.join(_ll)\n");
fprintf(fid, "js3d = ''\n");
fprintf(fid, "for u in ['https://3dmol.org/build/3Dmol-min.js',\n");
fprintf(fid, "          'https://3dmol.csb.pitt.edu/build/3Dmol-min.js']:\n");
fprintf(fid, "    try:\n");
fprintf(fid, "        req = urllib.request.Request(u, headers={'User-Agent': 'Mozilla/5.0'})\n");
fprintf(fid, "        js3d = urllib.request.urlopen(req, timeout=20).read().decode('utf-8', 'replace')\n");
fprintf(fid, "        break\n");
fprintf(fid, "    except: pass\n");
fprintf(fid, "if not js3d:\n");
fprintf(fid, "    print('ERROR: 3Dmol.js download failed from both CDNs', file=sys.stderr); sys.exit(1)\n");
fprintf(fid, "rec_js = json.dumps(rec_pdb)\n");
fprintf(fid, "lig_js = json.dumps(lig_pdb)\n");
fprintf(fid, "has_lig = 'true' if lig_pdb else 'false'\n");
fprintf(fid, "jv_lines = [\n");
fprintf(fid, "    'var v=$3Dmol.createViewer(document.getElementById(""v""),{backgroundColor:""#1a1a2e""});',\n");
fprintf(fid, "    'v.addModel(' + rec_js + ',""pdb"");',\n");
fprintf(fid, "    'v.setStyle({model:0,chain:""A""},{cartoon:{colorscheme:""ssJmol"",opacity:0.85}});',\n");
fprintf(fid, "    '// Active site residues (Chain A only): element-colored sticks to show binding pocket',\n");
fprintf(fid, "    '[41,49,145,163,165,166,187,189].forEach(function(r){',\n");
fprintf(fid, "    '    v.addStyle({model:0,resi:r,chain:""A""},{stick:{colorscheme:""element"",radius:0.18}});',\n");
fprintf(fid, "    '});',\n");
fprintf(fid, "    'if(' + has_lig + '){',\n");
fprintf(fid, "    '    v.addModel(' + lig_js + ',""pdb"");',\n");
fprintf(fid, "    '    v.setStyle({model:1},{stick:{colorscheme:""element"",radius:0.25}});',\n");
fprintf(fid, "    '    v.zoomTo({model:1});',\n");
fprintf(fid, "    '} else { v.zoomTo({model:0}); }',\n");
fprintf(fid, "    'v.render();'\n");
fprintf(fid, "]\n");
fprintf(fid, "jv = '\\n'.join(jv_lines)\n");
fprintf(fid, "css = ('html,body{margin:0;padding:0;width:100%%;height:100%%;overflow:hidden;background:#1a1a2e;}'\n");
fprintf(fid, "       '#v{width:100%%;height:100%%;display:block;position:absolute;top:0;left:0;}')\n");
fprintf(fid, "html = ('<html><head>'\n");
fprintf(fid, "        '<style>' + css + '</style>'\n");
fprintf(fid, "        '<script type=""text/javascript"">\\n' + js3d + '\\n</script>'\n");
fprintf(fid, "        '</head>'\n");
fprintf(fid, "        '<body>'\n");
fprintf(fid, "        '<div id=""v""></div>'\n");
fprintf(fid, "        '<script type=""text/javascript"">\\n' + jv + '\\n</script>'\n");
fprintf(fid, "        '</body></html>')\n");
fprintf(fid, "open(out, 'w').write(html)\n");
fclose(fid);
%[text] リガンドファイルを選択（ドッキング済みなら使用、未完なら 'none'）
if displayLigandName ~= "" && exist('displayPoseFile', 'var') && isfile(char(displayPoseFile))
    ligArg = char(displayPoseFile);
else
    ligArg = 'none';
end
%[text] HTML を生成（3Dmol.js を CDN から取得して inline 埋め込み）
vizRes = pySub6.run(py.list({pythonExe6, char(vizScript), ...
    char(preparedPdb), ligArg, char(htmlPath)}), ...
    capture_output=py.True, text=py.True);

if int32(vizRes.returncode) == 0
    logInfo("R08: ribbon HTML 保存 -> %s", htmlPath);
    % uifigure + uihtml で MATLAB Figure パネル内に ribbon を表示
    % （uifigure は App Designer と同じ WebGL 対応フレームを使用する）
    try
        figW = 960; figH = 780;
        figRibbon = uifigure( ...
            "Name",     sprintf("R08: Ribbon -- %s vs Mpro 6LU7", displayLigandName), ...
            "Position", [50, 50, figW, figH]);
        % uihtml の Position を uifigure と同サイズに揃えて全面表示
        uh = uihtml(figRibbon, "HTMLSource", char(fullfile(pwd, htmlPath)));
        uh.Position = [0, 0, figW, figH];
    catch ME2
        logWarn("R08: uifigure 表示失敗 (%s) -- web() にフォールバック", ME2.message);
        try
            web(char(fullfile(pwd, htmlPath)));
        catch
            logWarn("R08: Files パネルから '%s' を手動で開いてください", htmlPath);
        end
    end
else
    logWarn("R08: ribbon 生成失敗 (returncode=%d)", int32(vizRes.returncode));
    if strlength(string(char(vizRes.stderr))) > 0
        logWarn("R08: stderr:\n%s", char(vizRes.stderr));
    end
    if strlength(string(char(vizRes.stdout))) > 0
        logWarn("R08: stdout:\n%s", char(vizRes.stdout));
    end
end

logInfo("R08: Section 6 完了");
%[text] ### 可視化の読み取り方
%[text] 表示されたリボン図では次の点を確認してください。
%[text] 1. **ドッキング化合物の位置**: 太いスティック（中央）は β バレルドメインのポケット奥に収まっていますか。ポケット外縁（ループ付近）にある場合はドッキングが浅い可能性があります。
%[text] 2. **活性部位との近接**: His41（触媒一般塩基）と Cys145（触媒求核剤）の細いスティックはドッキング化合物の近くに見えますか。共有結合阻害剤（Nirmatrelvir）は Cys145 に極めて近接するはずです。
%[text] 3. **疎水性ポケット充填**: 黄色のスティック（Met49・Met165 の硫黄）はドッキング化合物の疎水性部位と近接していますか。P2 ポケットへの充填が高スコアの鍵です。
%[text] 4. **全体を回転させて多角的に確認**: ドラッグで回転し、ポケット深さや配向を 3 次元的に観察してください。 \
%%
%[text] ## セクション 7: 総合考察
%[text] ### 実験設計の意図
%[text] 次の 3 グループを比較することで、スコアと臨床的意義のズレを体験します。
%[text] - **A — Mpro 設計阻害剤**（Nirmatrelvir, Ensitrelvir）: 高スコアを期待します。
%[text] - **B — 交叉反応・臨床失敗**（Boceprevir〜Nelfinavir）: 中〜低スコアを示します。
%[text] - **C — 全く別標的**（Remdesivir〜Acyclovir）: 低スコアを期待します。 \
%[text] ### ドッキングスコアの限界
%[text] - スコアは結合自由エネルギーの近似値であり、in vivo 活性とは相関しないことがあります。
%[text] - 剛体受容体（Rigid receptor）近似は誘起適合（induced fit）効果を無視します（特に柔軟な活性部位に注意）。
%[text] - Nirmatrelvir は共有結合阻害剤であり、非共有結合近似では結合エネルギーが過小評価されます。 \
%[text] ### 薬物再利用（Drug Repurposing）の教訓
%[text] - 同じシステインプロテアーゼ族では active site geometry が類似するため、意外な交叉反応が起きます（Boceprevir の例）。
%[text] - スコアが良くても、ADMET（吸収・代謝・排泄・毒性）や PK が問題になることがあります。
%[text] - Lopinavir の臨床失敗は「Mpro への低親和性」＋「CYP3A4 自己誘導による血中濃度低下」の複合要因です。 \
%[text] ### 最終スコア一覧（コンソール出力）
logInfo("%-3s  %-20s  %-8s  %s", "#", "Compound", "Score", "Assessment");
logInfo("%s", repmat("-", 1, 60));
for i = 1:height(scored)
    sc = scored.DeltaG_kcal_mol(i);
    if sc <= -9
        assess = "*** STRONG binding";
    elseif sc <= -7
        assess = "**  Promising lead";
    elseif sc <= -5
        assess = "*   Moderate";
    else
        assess = "    Weak / No binding";
    end
    logInfo("%-3d  %-20s  %6.2f  %s", i, scored.Compound(i), sc, assess);
end
%[text] 演習完了です、おつかれさまでした。
%[text] 次のステップとして R04（タンパク質-リガンド解析）または R09（GNN 分子性質予測）に進んでみましょう。
%[text] ## まとめ
%[text] - 分子ドッキングのパイプライン（pdbfixer によるタンパク質前処理 → RDKit による 3D 構造生成 → meeko による PDBQT 変換 → AutoDock Vina 実行）を MATLAB から Python ライブラリを連携して一貫して実装できました。
%[text] - Vina のスコアリング関数 $\\Delta G$ は経験的な近似値で、共有結合阻害剤（Nirmatrelvir・Boceprevir）の Cys145 との共有結合寄与は非共有結合近似では計上されず、スコアが過小評価されます。
%[text] - 同じシステインプロテアーゼ族では活性部位の形状が類似するため、HCV 薬 Boceprevir のような「意図しない交叉反応」が計算・実験の両面で確認されます。
%[text] - HIV 薬（Lopinavir など）が高ドッキングスコアを示しても臨床で無効だった理由は、ADMET・CYP3A4 自己誘導による血中濃度低下などの複合要因であり、スコアと臨床活性は必ずしも一致しません。
%[text] - 3Dmol.js の自己完結 HTML を `uihtml` で表示することで、MATLAB Online 上でインタラクティブな 3D リボン可視化が実現できます。
%[text] - ドッキングシミュレーションは仮説生成ツールであり、実際の医薬品開発には IC50 測定・結晶構造解析・in vivo 試験による実験的検証が不可欠です。 \
%[text] Local helpers
%[text] 三項演算子ヘルパー（ローカル関数）
function out = ternary_(cond, trueVal, falseVal)
    if cond
        out = trueVal;
    else
        out = falseVal;
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
