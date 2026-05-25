# テストカタログ — EasyMolKit

🇬🇧 [English](../test_catalog.md)

> テストクラス一覧（25 クラス、869 テスト）。テストスイートの実行方法は [クイックスタートガイド（日本語）](quickstart_ja.md#テストの実行) を参照。

```matlab
% テスト実行
addpath(genpath("src"));
addpath(genpath("tests"));
suite   = testsuite("tests/unit");
runner  = matlab.unittest.TestRunner.withNoPlugins;
results = runner.run(suite);
fprintf("RESULT: %d PASS / %d FAIL / %d Total\n", ...
    sum([results.Passed]), sum([results.Failed]), numel(results));
```

---

## テストクラス一覧

| クラス名 | テスト数 | 対象モジュール | 主なテスト内容 |
|---|:---:|---|---|
| `TestMol` | 97 | `emk.mol` | SMILES 解析、正準化、部分構造検索、Mol シリアライズ、スキャフォールド抽出 |
| `TestDb` | 95 | `emk.db` | PubChem / ChEMBL 検索・活性データ取得（モック＋実サーバー） |
| `TestSimilarity` | 93 | `emk.similarity` | Tanimoto / Dice 係数、ランキング検索、N×N 類似度行列 |
| `TestDescriptor` | 83 | `emk.descriptor` | 10 標準記述子、バッチ計算、無効 Mol 処理 |
| `TestFilter` | 72 | `emk.filter` | Lipinski Ro5 フィルタ、違反数計算、NaN 処理 |
| `TestFingerprint` | 64 | `emk.fingerprint` | Morgan / MACCS 生成、`toArray()` 変換、ビット数検証 |
| `TestIo` | 49 | `emk.io` | SDF 読み書き、SMILES ファイル読み込み、無効 Mol スキップ |
| `TestPubchemFetch` | 36 | `emk.db` | PubChemPy 詳細取得（CID / 名前 / SMILES / InChIKey） |
| `TestToTable` | 34 | `emk.mol` | `toTable()` 列構成・型検証、無効 Mol の NaN 置換 |
| `TestMordred` | 27 | `emk.descriptor` | Mordred 2D 記述子計算（~1800 種）、バッチ処理 |
| `TestInstallExtra` | 24 | `emk.setup` | Track 1 ライブラリのインストール・検証・バージョン管理 |
| `TestRdkitModule` | 25 | `emk.setup` | RDKit サブモジュール importlib 確認、バージョン取得 |
| `TestRecipe` | 23 | `emk.setup` | `recipe()` 表示、未知ライブラリエラー |
| `TestValidate` | 22 | `emk.setup` | `validate()` テーブル形式・インストール状態検証 |
| `TestLoadConfig` | 15 | `src/config` | 設定ロード優先順位（環境変数 > JSON > デフォルト） |
| `TestInstallTrack2` | 16 | `emk.setup` | Track 2 venv 作成・接続・設定書き込み |
| `TestLogHelpers` | 12 | `src/util` | ログ出力フォーマット（INFO / WARN / ERROR / DEBUG） |
| `TestVerify` | 12 | `emk.setup` | `verify()` の非スロー動作、Python / RDKit 状態判定 |
| `TestInitPython` | 10 | `emk.setup` | `initPython()` 冪等性、pyenv 設定、プラットフォーム検出 |
| `TestUseExternal` | 10 | `emk.setup` | 外部 CPython パス指定、無効パスエラー |
| `TestIsOnline` | 6 | `emk.util` | MATLAB Online / Desktop 環境判定 |
| `TestInstall` | 7 | `emk.setup` | `install()` デスクトップ配備フロー |
| `TestInstallOnline` | 7 | `emk.setup` | `installOnline()` オンライン配備フロー |
| `TestViz` | 25 | `emk.viz` | `draw2d()` 描画・Figure 返却・入力バリデーション |
| `TestMakeRunDir` | 5 | `src/util` | `makeRunDir()` ディレクトリ生成・タイムスタンプ形式 |

**合計: 25 クラス / 869 テスト**

---

## テスト種別とタグ

各テストは `assumeTrue` で **RDKit が不要なもの** から先に実行される設計です。

| タグ | 意味 |
|---|---|
| RDKit 不要 | 入力バリデーション、設定ロード、ログ、ユーティリティ |
| RDKit 必要 | 分子操作、記述子計算、フィンガープリント、類似度 |
| ネットワーク必要 | PubChem / ChEMBL 実サーバーテスト（`emk.db` の一部） |

> RDKit 未インストール環境では `assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available")` により
> 該当テストが自動的にスキップされます。

---

## 関連ドキュメント

- [function_reference.md](../function_reference.md) — 完全なシグネチャ・エラー仕様
- [function_reference_ja.md](function_reference_ja.md) — 関数リファレンス（日本語コンパクト版）
- [algorithm_guide.md](../algorithm_guide.md) — アルゴリズム根拠・テスト検証戦略
