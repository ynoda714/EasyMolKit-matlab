# クイックスタートガイド — EasyMolKit

> 英語版 → [../en/quickstart.md](../en/quickstart.md)

---

## 前提条件

| 項目 | Desktop | MATLAB Online |
|---|---|---|
| MATLAB | R2025b 以降 | R2025b 以降 |
| ネットワーク | 初回セットアップに必要 | 必要 |
| Python | 不要（自動配備） | システム Python が使用される |
| `config/settings.json` | `settings.example.json` からコピー | 同様 |

---

## 初回セットアップ

### MATLAB Desktop

```matlab
% 1. リポジトリをクローン
%    git clone https://github.com/ynoda714/EasyMolKit-matlab.git

% 2. MATLAB で EasyMolKit フォルダを開く

% 3. セットアップを実行（1 回のみ — Python + RDKit を自動配備）
addpath(genpath("src"));
emk.setup.install();
```

> `emk.setup.install()` は以下を自動実行します:
> 1. プラットフォーム検出（Windows Desktop / MATLAB Online）
> 2. Embedded Python 3.10 を `python_env/` にダウンロード・展開
> 3. pip で RDKit をインストール
> 4. `pyenv`（OutOfProcess モード）を設定
> 5. RDKit のインポートを検証

### MATLAB Online

```matlab
addpath(genpath("src"));
emk.setup.installOnline();
```

> MATLAB Online では `get-pip.py` ブートストラップ + `!~/.local/bin/pip install rdkit-pypi` を
> プリインストール済み Python に対して実行します。

### MCP サーバー設定（開発者のみ）

`.vscode/mcp.json` を開き、2 つのパスエントリをプロジェクトの絶対パスに更新:

```json
"--initial-working-folder", "D:\\workspace\\EasyMolKit",
"--log-folder", "D:\\workspace\\EasyMolKit\\logs\\matlab-mcp",
```

### 設定ファイルの準備

```powershell
copy config\settings.example.json config\settings.json
# settings.json を編集してオプションを設定
```

---

## 基本的な使い方

### Step 1. セクション実行で体験する

```matlab
addpath(genpath("src"));
main_emk    % Run Section（Ctrl+Enter）でセクションを 1 つずつ実行
```

> ⚠️ **Run Section（Ctrl+Enter）** を使用してください — "Run"（F5）は全セクションを一括評価します。

### Step 2. API を直接呼び出す

```matlab
% 分子を作成
mol = emk.mol.fromSmiles("CCO");          % エタノール

% 記述子を計算
mw   = emk.descriptor.molWeight(mol);     % 分子量
logP = emk.descriptor.calculate(mol, "LogP");

% フィンガープリント
fp = emk.fingerprint.morgan(mol);

% 類似度計算
mol2  = emk.mol.fromSmiles("CCCO");
fp2   = emk.fingerprint.morgan(mol2);
score = emk.similarity.tanimoto(fp, fp2);

% SDF ファイル入出力
mols = emk.io.readSdf("data/sample.sdf");
```

### Step 3. 出力成果物を確認する

```
result/runs/<YYYYMMDD_HHMMSS>/
  ├─ run_meta.json
  └─ *.csv / *.mat
```

---

## テストの実行

```matlab
% ユニットテスト
addpath(genpath("src"));
suite   = testsuite("tests/unit");
runner  = matlab.unittest.TestRunner.withNoPlugins;
results = runner.run(suite);
fprintf("RESULT: %d PASS / %d FAIL / %d Total\n", ...
    sum([results.Passed]), sum([results.Failed]), numel(results));

% スモークテスト
addpath(genpath("src")); addpath("tests/smoke");
test_mvp_smoke();
```

---

## 3 層 UX

| 層 | 入口 | 対象者 | 特徴 |
|---|---|---|---|
| **Tier 1: 今すぐ体験** | `main_emk.m` | 初回ユーザー | Ctrl+Enter でセクションを順番に実行; プリセット SMILES で即体験 |
| **Tier 2: カスタマイズ** | `main_<feature>.m` | 一般ユーザー | パラメータ制御、セクション実行 |
| **Tier 3: API 直接利用** | `emk.*` 関数群 | 開発者・上級者 | バッチ処理、カスタムスクリプト |

---

## FAQ

**Q: `emk.setup.install()` が失敗する**
A: ネットワーク接続を確認してください。プロキシ環境の場合は `config/settings.json` の `python.proxy` を設定してください。

**Q: パスが認識されない**
A: `addpath(genpath("src"))` を実行してください。`main_emk.m` の Section 0a で自動実行されます。

**Q: 設定ファイルが見つからない**
A: `config/settings.example.json` を `config/settings.json` にコピーしてください。

**Q: `pyenv` からエラーが出る**
A: `pyenv` で別の Python がすでに設定されている場合は、MATLAB を再起動して `emk.setup.initPython()` を実行してください。
`pyenv` の Version は MATLAB セッション中に 1 度しか設定できません。

**Q: MATLAB Online で RDKit が使えない**
A: `emk.setup.installOnline()` を実行してください。内部で `get-pip.py` ブートストラップ + `!~/.local/bin/pip install rdkit-pypi` が実行されます。

**Q: MCP サーバーに接続できない**
A: `.vscode/mcp.json` のパスが正しいか確認してください。MATLAB が起動している必要があります。

**Q: 出力成果物が見つからない**
A: `result/runs/` 配下の最新タイムスタンプフォルダを確認してください。

---

## Track 2: 外部 CPython 環境

GPL ライセンスや技術的制約により、Open Babel・MDAnalysis・PyMOL は Embedded Python に
バンドルできません。**別途 CPython 3.10+ 環境**にインストールし、`emk.setup.useExternal()` で接続します。

> ⚠️ **重要な制約**: `emk.setup.useExternal()` は Python がロードされる **前**に呼び出してください。
> 一度ロードされると同一セッション内で変更できません（MATLAB `pyenv` の制限）。
> Track 2 使用時はスクリプト先頭の `addpath` 直後に呼び出してください。

### Step 1: 外部 Python 環境のセットアップ

ライブラリごとに独立した環境を推奨します。

**MDAnalysis（MD トラジェクトリ解析）**
```powershell
python -m venv C:\envs\mdenv
C:\envs\mdenv\Scripts\pip install MDAnalysis
```

**Open Babel（化学ファイル形式変換 & 3D 座標生成）**
```
1. https://github.com/openbabel/openbabel/releases から Windows インストーラー（Python バインディング付き）をダウンロード
2. CPython 3.10+ 環境にインストール（python_env/ には使用不可）
3. 以下の Step 2 で useExternal() を呼び出す
```

**PyMOL OSS（3D 可視化）**
```powershell
python -m venv C:\envs\pymolenv
C:\envs\pymolenv\Scripts\pip install pymol-open-source
```

### Step 2: MATLAB セッション開始時に接続

```matlab
addpath(genpath("src"));

% Python がロードされる前に useExternal() を呼び出す
emk.setup.useExternal("C:\envs\mdenv\python.exe")

% 以降は通常どおり emk.* 関数を使用
mol = emk.mol.fromSmiles("CCO");
```

> 毎回呼び出すのを省略するには `config/settings.json` に設定してください:
> ```json
> { "python": { "external_path": "C:\\envs\\mdenv\\python.exe" } }
> ```

### Step 3: 接続確認

```matlab
T = emk.setup.validate()
emk.setup.recipe("mdanalysis")   % ライブラリごとの詳細手順
```

### Track 2 利用可能ライブラリ

| ライブラリ | 用途 | ライセンス |
|---|---|---|
| **openbabel** | 化学ファイル変換（110 以上の形式）& 3D 座標生成 | GPLv2 |
| **mdanalysis** | MD トラジェクトリ解析（GROMACS / AMBER / NAMD 等） | GPLv2+ |
| **pymol** | 3D 分子可視化（PyMOL Open-Source） | Python-2.0/BSD |

> ⚠️ Open Babel と MDAnalysis は GPLv2 / GPLv2+ ライブラリです。これらを利用するスクリプトは
> GPL 条件に従う場合があります。商用利用前に [docs/compliance.md](../compliance.md) を確認してください。

---

## ログ出力形式

```
[HH:MM:SS][INFO]  処理完了（42 分子）
[HH:MM:SS][WARN]  5 行目の無効 SMILES — スキップ
[HH:MM:SS][ERROR] RDKit インポート失敗: モジュールが見つかりません
[####------]  40% ( 4/10) molecules
```
