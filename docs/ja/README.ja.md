# EasyMolKit

🇬🇧 [English](../../README.md)

**MATLAB から RDKit を簡単に使えるケモインフォマティクス統合環境**

> リポジトリをクローンしてセットアップコマンドを 1 本実行するだけで、RDKit が MATLAB 上で使えるようになります

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=ynoda714/EasyMolKit-matlab)

## なぜ EasyMolKit？

ケモインフォマティクス分野では Python + RDKit が事実上の標準ツールチェーンです。
しかし、その環境構築には多くのハードルがあります。

- Python バージョンや仮想環境の管理
- conda / pip による複雑な RDKit インストール
- 商用ツール（PyMOL など）との環境競合
- MATLAB ユーザーには馴染みのない Python エコシステム

**EasyMolKit はこれらのハードルをすべて解消します。** MATLAB の `pyenv` ベースの Python 連携を活用し、
RDKit の機能を通常の MATLAB 関数として提供します — Python の知識は一切不要です。

## 特徴

- **ゼロ設定**: `emk.setup.install()` を 1 回呼ぶだけで Python + RDKit を自動配備
- **MATLAB ネイティブ**: 結果は `table` / `struct` / `double` で返却 — ワークスペースで即利用可能
- **Desktop & Online 対応**: Windows Desktop および MATLAB Online をサポート（macOS / Linux Desktop は未検証）
- **RDKit ラッパー群**: 分子解析・記述子計算・フィンガープリント・類似度検索をすべて MATLAB 関数として提供
- **将来拡張**: PyMOL Open-Source による 3D 構造可視化を将来のリリースで予定

## 対象ユーザー

- **化学・薬学・医学分野の研究者**で MATLAB をメインの研究環境として使用している方
- **ケモインフォマティクスを学ぶ学生**（MATLAB Online Basic 無料枠で Layer 1〜3 が動作）
- **Python 環境構築に時間をかけたくない MATLAB ユーザー**

## 動作要件

| 項目 | Desktop | MATLAB Online |
|---|---|---|
| MATLAB | R2025b 以降 | R2025b 以降 |
| Python | 自動配備（Embedded Python） | システム Python が使用される |
| RDKit | 自動配備 | `emk.setup.installOnline()` でインストール |
| OS | Windows | — |

## クイックスタート

```matlab
% 1. リポジトリをクローン
%    git clone https://github.com/ynoda714/EasyMolKit-matlab.git
%    cd EasyMolKit-matlab

% 2. MATLAB で開き、初回セットアップを実行（1 回だけ）
addpath(genpath("src"));
emk.setup.install();          % Desktop
% emk.setup.installOnline();  % MATLAB Online

% 3. 試してみよう
mol = emk.mol.fromSmiles("CCO");          % エタノール
mw  = emk.descriptor.molWeight(mol);      % 分子量
fp  = emk.fingerprint.morgan(mol);        % Morgan フィンガープリント

% 4. 分子間の類似度を比較
mol2  = emk.mol.fromSmiles("CCCO");       % プロパノール
fp2   = emk.fingerprint.morgan(mol2);
score = emk.similarity.tanimoto(fp, fp2); % Tanimoto 係数
```

詳細は [docs/quickstart.md](../quickstart.md) を参照してください。

## 追加ライブラリ（Track 1 & Track 2）

EasyMolKit はアドオンライブラリを 2 つのトラックで管理します。

| トラック | ライブラリ | インストール方法 | ライセンス |
|---|---|---|---|
| **Track 1** | pubchempy, mordred, biopython, torch, torch_geometric, transformers, datasets 等 | `emk.setup.installExtra()` — Embedded Python に直接追加 | MIT / BSD-3 / Apache-2.0 |
| **Track 2** | Open Babel, MDAnalysis, PyMOL OSS | 別途 CPython 環境が必要; `emk.setup.useExternal()` で接続 | GPLv2 / GPLv2+ / BSD |

### Track 1: Embedded Python への追加インストール

```matlab
% インストール前にレシピとライセンスを確認
emk.setup.recipe("pubchempy")          % インストール手順とライセンスを表示
emk.setup.installExtra("pubchempy")    % Embedded Python にインストール
emk.setup.installExtra("mordred")      % 1800 種以上の記述子ライブラリ
emk.setup.installExtra("biopython")    % PDB / 配列解析

% PyTorch + HuggingFace スタック（R09 / R10 で必要）
emk.setup.installExtra("torch")           % CPU のみ、約 800 MB（最初にインストール）
emk.setup.installExtra("torch_geometric") % GNN ライブラリ（torch が必要）
emk.setup.installExtra("transformers")    % HuggingFace Transformers
emk.setup.installExtra("datasets")        % HuggingFace Datasets

% インストール状況を確認
T = emk.setup.validate()
```

> 全ライブラリの一括インストールは `main_setup_extra.m` を参照。

### Track 2: 外部 CPython 環境が必要なライブラリ

GPL ライセンスや技術的制約により、Open Babel・MDAnalysis・PyMOL は Embedded Python に
バンドルできません。**別途 CPython 3.10+ 環境**にインストールし、`emk.setup.useExternal()` で接続します。

> ⚠️ **重要な制約**: `emk.setup.useExternal()` は Python がロードされる **前**に呼び出してください。
> 一度ロードされると同一セッション内で変更できません（MATLAB `pyenv` の制限）。
> Track 2 使用時はスクリプト先頭の `addpath` 直後に呼び出してください。

#### Step 1: 外部 Python 環境のセットアップ

ライブラリごとに独立した環境を推奨します。例:

**MDAnalysis（MD トラジェクトリ解析）**
```powershell
python -m venv C:\envs\mdenv
C:\envs\mdenv\Scripts\pip install MDAnalysis
```

**Open Babel（化学ファイル形式変換 & 3D 座標生成）**
```
1. https://github.com/openbabel/openbabel/releases から Windows インストーラー（Python バインディング付き）をダウンロード
2. CPython 3.10+ 環境にインストール（python_env/ には使用不可）
3. 以下と同様に useExternal() を呼び出す
```

**PyMOL OSS（3D 可視化）**
```powershell
python -m venv C:\envs\pymolenv
C:\envs\pymolenv\Scripts\pip install pymol-open-source
```

#### Step 2: MATLAB セッション開始時に接続

```matlab
addpath(genpath("src"));

% Python がロードされる前に useExternal() を呼び出す
emk.setup.useExternal("C:\envs\mdenv\python.exe")

% 以降は通常どおり emk.* 関数を使用
mol = emk.mol.fromSmiles("CCO");
```

> 毎回呼び出すのを省略するには `config/settings.json` に設定してください:
>
> ```json
> {
>   "python": {
>     "external_path": "C:\\envs\\mdenv\\python.exe"
>   }
> }
> ```

#### Step 3: 接続確認

```matlab
T = emk.setup.validate()

emk.setup.recipe("openbabel")
emk.setup.recipe("mdanalysis")
emk.setup.recipe("pymol")
```

#### Track 2 利用可能ライブラリ

| ライブラリ | 用途 | pip コマンド | ライセンス |
|---|---|---|---|
| **openbabel** | 化学ファイル変換（110 以上の形式）& 3D 座標生成 | Windows インストーラー必須（上記参照） | GPLv2 |
| **mdanalysis** | MD トラジェクトリ解析（GROMACS / AMBER / NAMD 等） | `pip install MDAnalysis` | GPLv2+ |
| **pymol** | 3D 分子可視化（PyMOL Open-Source） | `pip install pymol-open-source` | Python-2.0/BSD |

> ⚠️ **GPL ライセンスに関する注意**: Open Babel と MDAnalysis は GPLv2 / GPLv2+ ライブラリです。
> EasyMolKit 本体（MIT）には影響しませんが、これらを利用するスクリプトは GPL 条件に
> 従う場合があります。商用利用前に [docs/compliance.md](../compliance.md) を確認してください。

## チュートリアル & サンプル（4 層構造）

EasyMolKit は `examples/` 配下に段階的な学習コンテンツを提供します。

| 層 | 対象者 | 必要ツールボックス | コンテンツ | リリース状況 |
|---|---|---|---|---|
| **L1 Foundation** | 全ユーザー | なし | API を 1 概念ずつ学ぶ（6 モジュール、各 5〜15 分） | ✅ v1.0.0 |
| **L2 Application Stories** | Foundation 修了後 | なし | 複数機能を組み合わせた実践ワークフロー（7 モジュール、各 20〜40 分） | ✅ v1.1.0 |
| **L3 Analytics** | 全ユーザー | 要統計 ML 等 | QSAR・クラスタリング・MS 解析・最適化など（A01〜A10、各 30〜60 分） | ✅ v1.2.0 |
| **L4 Research** | 全ユーザー | 要並列計算等 | 研究レベルの応用（R01〜R10、各 30〜90 分） | ✅ v1.3.0 |

*L1〜L3 は MATLAB Online Basic（無料枠）で完全動作します。*

### Layer 1: Foundation（Base MATLAB のみ）

| # | タイトル | 必要Toolbox | 学習内容 | Desktop | Online |
|---|---|---|---|:---:|:---:|
| F01 | SMILES で分子を描く | なし | 分子表現、SMILES 構文 | ✔ | ✔ |
| F02 | 分子物性を計算する | なし | MW / LogP / TPSA の意味と計算 | ✔ | ✔ |
| F03 | フィンガープリント入門 | なし | ビットベクトル表現、Morgan vs MACCS | ✔ | ✔ |
| F04 | 類似度で分子を比較する | なし | Tanimoto / Dice 類似度の定量化 | ✔ | ✔ |
| F05 | 部分構造検索 | なし | SMARTS パターンマッチング | ✔ | ✔ |
| F06 | ファイルから分子を読み込む | なし | SDF / SMILES ファイル操作 | ✔ | ✔ |

### Layer 2: Application Stories（Base MATLAB のみ）

| # | タイトル | 必要Toolbox | ドメイン | Desktop | Online |
|---|---|---|---|:---:|:---:|
| S01 | カフェインの仲間を探す | なし | 身近な化学 | ✔ | ✔ |
| S02 | 創薬フィルター：リピンスキーの Rule of Five | なし | 薬理学 | ✔ | ✔ |
| S03 | 危険化合物の構造アラート | なし | 安全性 | ✔ | ✔ |
| S04 | バーチャルスクリーニング入門 | なし | 創薬 | ✔ | ✔ |
| S05 | 未知化合物同定チャレンジ | なし | 法科学 | ✔ | ✔ |
| S06 | PubChem で化合物を検索する | なし | データベース | ✔ | ✔ |
| S07 | ChEMBL の活性データを解析する | なし | 創薬 | ✔ | ✔ |

### Layer 3: Analytics

| # | タイトル | 必要Toolbox | 内容 | Desktop | Online |
|---|---|---|---|:---:|:---:|
| A01 | PCA による化学空間マッピング | Statistics and Machine Learning Toolbox | 次元削減・化学空間可視化 | ✔ | ✔ |
| A02 | 分子クラスタリング | Statistics and Machine Learning Toolbox | 階層クラスタリング・構造類似性 | ✔ | ✔ |
| A03 | QSAR 回帰 | Statistics and Machine Learning Toolbox | LogP 予測・回帰モデル評価 | ✔ | ✔ |
| A04 | 薬物分類 | Statistics and Machine Learning Toolbox | SVM / ランダムフォレスト・ROC 曲線 | ✔ | ✔ |
| A05 | ニューラルネットワーク物性予測 | Deep Learning Toolbox | フィードフォワード NN・物性回帰 | ✔ | ✔ |
| A06 | 用量反応カーブフィット | Curve Fitting Toolbox | Hill 式・EC50 推定 | ✔ | ✔ |
| A07 | Scaffold 分析と R グループ分解 | Statistics and Machine Learning Toolbox | 創薬化学分析 | ✔ | ✔ |
| A08 | 質量分析×ケモインフォマティクス | Signal Processing Toolbox + Statistics and Machine Learning Toolbox | 同位体パターンマッチング・MS アノテーション | ✔ | ✔ |
| A09 | PFAS・環境スクリーニング | Optimization Toolbox + Statistics and Machine Learning Toolbox | SMARTS スクリーニング・Pareto 最適化 | ✔ | ✔ |
| A10 | リード最適化 | Optimization Toolbox | 多目的最適化・Derringer-Suich 法 | ✔ | ✔ |

### Layer 4: Research

| # | タイトル | 必要ツールボックス | Desktop | Online |
|---|---|---|:---:|:---:|
| R01 | 大規模類似度スクリーニング（GPU）| Parallel Computing Toolbox (GPU) | ✔ | △（CPU Only）|
| R02 | PK/PD シミュレーション | SimBiology | ✔ | ✔ |
| R03 | 法科学ケモメトリクス | Statistics and Machine Learning Toolbox + Parallel Computing Toolbox | ✔ | ✔ |
| R04 | タンパク質-リガンド解析 † | Bioinformatics Toolbox | ✔ | ✔ |
| R05 | 分子言語モデル：SMILES 生成 | Deep Learning Toolbox | ✔ | ✔ |
| R06 | REINFORCE 分子設計 | Deep Learning Toolbox + Reinforcement Learning Toolbox | ✔ | ✔ |
| R07 | メタボロミクス † | Bioinformatics Toolbox + SimBiology | ✔ | ✔ |
| R08 | タンパク質-リガンド ドッキングシミュレーション ‡ | なし（Track 1: meeko + vina + pdbfixer）| ✕ | ✔ |
| R09 | GNN 分子性質予測 § | Deep Learning Toolbox | ✔ | ✔ |
| R10 | ChemBERTa 転移学習 § | Deep Learning Toolbox | ✔ | ✔ |

> **†** 初回実行前に `emk.setup.installExtra("biopython")` が必要です（Track 1 追加ライブラリ。MATLAB ライセンスとは別要件）。
>
> **‡ MATLAB Online 限定**（Windows Desktop 非対応: vina は Windows PyPI ホイールなし・pdbfixer の openmm は Smart App Control でブロック）。
> セットアップ: `main_rdkit.m` で `cfg.optionalLibraries.meeko/vina/pdbfixer = true` に設定 → `installOnline(Config=cfg)` で一括導入。
>
> **§** PyTorch + HuggingFace スタックが必要です。以下の順でインストールしてください:
> `emk.setup.installExtra("torch")` → `emk.setup.installExtra("torch_geometric")` → `emk.setup.installExtra("transformers")` → `emk.setup.installExtra("datasets")`。
> R10 は R09 の torch 環境が前提です。

## ディレクトリ構成

```
EasyMolKit/
├─ main_rdkit.m               # RDKit セットアップ & 基本操作（セクション実行）
├─ config/
│   └─ settings.example.json  # 設定テンプレート
├─ examples/
│   ├─ japanese/              # 配布教材（日本語版、plain-text Live Code）
│   │   ├─ foundation/        #   L1: API 基礎（Base MATLAB のみ）
│   │   ├─ stories/           #   L2: 応用ストーリー（Base MATLAB のみ）
│   │   ├─ analytics/         #   L3: 統計 & ML 統合（A01〜A10）
│   │   └─ research/          #   L4: 研究レベル（リリース予定）
│   └─ english/               # 配布教材（英語版、コメントのみ異なる）
│       ├─ foundation/
│       ├─ stories/
│       ├─ analytics/
│       └─ research/
├─ src/
│   ├─ +emk/                  # メインパッケージ
│   │   ├─ +setup/            # Python 環境セットアップ
│   │   ├─ +mol/              # 分子オブジェクト操作
│   │   ├─ +descriptor/       # 記述子計算
│   │   ├─ +fingerprint/      # フィンガープリント生成
│   │   ├─ +similarity/       # 類似度計算
│   │   ├─ +filter/           # 分子フィルタリング（Lipinski 等）
│   │   ├─ +io/               # ファイル入出力
│   │   ├─ +viz/              # 可視化（将来 PyMOL 連携）
│   │   └─ +util/             # パッケージユーティリティ
│   ├─ config/                # 設定ローダー
│   └─ util/                  # ログヘルパー & 共通ユーティリティ
├─ result/                    # 実行成果物（Git 非追跡）
├─ tests/
│   ├─ unit/                  # matlab.unittest クラスベーステスト
│   └─ smoke/                 # スモークテスト
├─ data/                      # キュレーション済みサンプルデータ
└─ docs/                      # ドキュメント
```

## API 概要

| モジュール | 主な関数 | 説明 |
|---|---|---|
| `emk.setup` | `install()`, `verify()`, `initPython()` | Python 環境の自動配備 & 初期化 |
| `emk.mol` | `fromSmiles()`, `toSmiles()`, `isValid()`, `hasSubstruct()` | 分子オブジェクトの作成 & 変換 |
| `emk.descriptor` | `molWeight()`, `calculate()`, `batchCalculate()` | 分子記述子の計算 |
| `emk.fingerprint` | `morgan()`, `maccs()`, `toArray()` | フィンガープリント生成 |
| `emk.similarity` | `tanimoto()`, `dice()` | 分子類似度計算 |
| `emk.io` | `readSdf()`, `writeSdf()`, `readSmilesList()` | SDF / SMILES ファイル入出力 |
| `emk.viz` | `draw2d()` | 2D 構造描画（※） |

> **※ 描画に関する注意**: `emk.viz.draw2d()` は RDKit（Python）で PNG を生成し MATLAB に転送します。
> 描画には **1 分子あたり 0.5〜2 秒**かかります。MATLAB Online では
> プロセス間通信のオーバーヘッドがさらに加わります（構造上の制約であり改善不可）。

全 API の詳細は [docs/function_reference.md](../function_reference.md) または [docs/function_catalog.ja.md](../function_catalog.ja.md) を参照してください。

## 主要な規約

- ロジックはすべて `src/` 配下に配置（`main_<feature>.m` エントリポイントを除く）
- `.m` ファイルは英語のみ（コメント・ログ・エラーメッセージ）
- `string` 型（`"..."` リテラル）を優先; R2025b 以降が必要
- `py.rdkit.*` を直接呼び出さない — 必ず `emk.*` ラッパー経由
- 出力はすべて `logInfo` / `logWarn` / `logError` を使用（`fprintf` は禁止）
- 実行成果物は `result/runs/<YYYYMMDD_HHMMSS>/` に保存（Git 非追跡）

## ライセンス

EasyMolKit: [MIT License](../../LICENSE)

### サードパーティライセンス

| ライブラリ | ライセンス | 用途 |
|---|---|---|
| RDKit | BSD-3-Clause | ケモインフォマティクスコア |
| Python (CPython) | PSF License | ランタイム環境 |
| PyMOL Open-Source | Python/BSD-like | 3D 可視化（将来リリース） |

詳細は [THIRD_PARTY_NOTICES.md](../../THIRD_PARTY_NOTICES.md) および [docs/compliance.md](../compliance.md) を参照してください。

## コントリビューション

バグ報告・機能リクエスト・プルリクエストを歓迎します。
[CONTRIBUTING.md](../../CONTRIBUTING.md) のガイドラインを参照してください。

## 免責事項

本ソフトウェアは研究・教育目的で提供されています。

- 本ソフトウェアは明示的または黙示的な保証なしに「現状のまま」提供されます
- 開発者は本ソフトウェアの使用に起因するいかなる損害についても責任を負いません
- 外部データソース（PubChem、ChEMBL 等）の利用はそれぞれの利用規約に従います
- 予測・計算結果は研究目的のみ; 医療・安全判断への直接適用には専門家のレビューが必要です

## ドキュメント

| ファイル | 説明 |
|---|---|
| [docs/quickstart.md](../quickstart.md) | セットアップ手順 & FAQ |
| [docs/function_reference.md](../function_reference.md) | 全関数シグネチャリファレンス |
| [docs/function_catalog.ja.md](../function_catalog.ja.md) | コンパクト関数カタログ（日本語） |
| [docs/test_catalog.ja.md](../test_catalog.ja.md) | テストクラスカタログ（日本語） |
| [docs/python_integration.md](../python_integration.md) | Python 連携アーキテクチャ |
| [docs/platform_support.md](../platform_support.md) | Desktop / Online プラットフォーム対応 |
| [docs/compliance.md](../compliance.md) | ライセンス & コンプライアンス |
