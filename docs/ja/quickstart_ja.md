# クイックスタートガイド — EasyMolKit

🇬🇧 [English](../quickstart.md)

> Version: v0.2.0 / Updated: 2026-04-19

---

## 前提条件

| 項目 | Desktop | MATLAB Online |
|---|---|---|
| MATLAB | R2025b 以降 | R2025b 以降 |
| ネット接続 | 初回セットアップ時に必要 | 必要 |
| Python | 不要（自動配備） | 事前インストール済み |
| `config/settings.json` | `settings.example.json` からコピー | 同左 |

---

## 初回セットアップ

### MATLAB Desktop

```matlab
% 1. リポジトリをクローン
%    git clone https://github.com/ynoda714/EasyMolKit-matlab.git

% 2. MATLAB で EasyMolKit フォルダを開く

% 3. セットアップを実行（1 回のみ — Embedded Python + RDKit を自動配備）
addpath(genpath("src"));
emk.setup.install();
```

> `emk.setup.install()` は以下を自動実行します:
> 1. プラットフォーム検出（Windows Desktop / MATLAB Online）
> 2. Embedded Python 3.10 を `python_env/` にダウンロード・展開
> 3. pip で RDKit をインストール
> 4. `pyenv` を OutOfProcess モードで設定
> 5. RDKit インポートを検証

### MATLAB Online

```matlab
addpath(genpath("src"));
emk.setup.installOnline();
```

> MATLAB Online では、プリインストール済み Python に対して `get-pip.py` ブートストラップ
> `!~/.local/bin/pip install rdkit-pypi` を実行します。

### MCP サーバー設定（開発者向け）

`.vscode/mcp.json` を開き、2 箇所のパスエントリをプロジェクトの絶対パスに更新してください:

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
main_rdkit    % セクションを 1 つずつ Ctrl+Enter で実行
```

> ⚠️ **Run Section（Ctrl+Enter）** を使用してください。「Run」（F5）は押さないでください。
> F5 を押すとすべてのセクションが一度に評価されます。

### Step 2. API を直接使う

```matlab
% 分子の作成
mol = emk.mol.fromSmiles("CCO");          % エタノール

% 記述子の計算
mw   = emk.descriptor.molWeight(mol);     % 分子量
logP = emk.descriptor.calculate(mol, "LogP");

% フィンガープリント
fp = emk.fingerprint.morgan(mol);

% 類似度計算
mol2  = emk.mol.fromSmiles("CCCO");
fp2   = emk.fingerprint.morgan(mol2);
score = emk.similarity.tanimoto(fp, fp2);

% SDF ファイルの入出力
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

> テストクラス一覧 → [テストカタログ](test_catalog_ja.md)

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

## 3 層 UX 構成

| 層 | エントリーポイント | 対象ユーザー | 特徴 |
|---|---|---|---|
| **Tier 1: すぐ試す** | `main_rdkit.m` | 初めての方 | Ctrl+Enter でセクションを順に実行。プリセット SMILES でそのまま体験 |
| **Tier 2: カスタマイズ** | `main_<feature>.m` | 通常ユーザー | パラメーター制御、セクション実行 |
| **Tier 3: 直接 API** | `emk.*` 関数 | 開発者・上級者 | バッチ処理、カスタムスクリプト |

---

## よくある質問

**Q: `emk.setup.install()` が失敗する**  
A: ネットワーク接続を確認してください。プロキシ環境では `config/settings.json` の `python.proxy` を設定してください。

**Q: パスが認識されない**  
A: `addpath(genpath("src"))` を実行してください。`main_rdkit.m` のセクション 0a で自動実行されます。

**Q: 設定ファイルが見つからない**  
A: `config/settings.example.json` を `config/settings.json` にコピーしてください。

**Q: `pyenv` でエラーが発生する**  
A: すでに別の Python が `pyenv` に設定されている場合は、MATLAB を再起動して `emk.setup.initPython()` を実行してください。
`pyenv` の Version は 1 MATLAB セッション内で 1 回のみ設定できます。

**Q: MATLAB Online で RDKit が使えない**  
A: `emk.setup.installOnline()` を実行してください。内部的に `get-pip.py` ブートストラップ `!~/.local/bin/pip install rdkit-pypi` を実行します。

**Q: MCP サーバーに接続できない**  
A: `.vscode/mcp.json` のパスが正しいか確認してください。MATLAB が起動している必要があります。

**Q: 出力成果物が見つからない**  
A: `result/runs/` 以下の最新のタイムスタンプ付きフォルダを確認してください。

---

## ログ出力フォーマット

```
[HH:MM:SS][INFO]  処理完了 (42 分子)
[HH:MM:SS][WARN]  行 5 の無効な SMILES — スキップ
[HH:MM:SS][ERROR] RDKit インポート失敗: モジュールが見つかりません
[####------]  40% ( 4/10) molecules
```
