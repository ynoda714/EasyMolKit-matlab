# テストカタログ — EasyMolKit

> ユニットテスト 25 クラス + スモークテスト 2 スクリプト。
> English version → [test_catalog.md](test_catalog.md)

## テストの実行方法

```matlab
addpath(genpath("src"));
addpath(genpath("tests"));   % クラス名ルックアップに必要

% ユニットテストを全て実行
suite   = testsuite("tests/unit");
runner  = matlab.unittest.TestRunner.withNoPlugins;
results = runner.run(suite);
fprintf("RESULT: %d PASS / %d FAIL / %d Total\n", ...
    sum([results.Passed]), sum([results.Failed]), numel(results));

% スモークテストを実行
run("tests/smoke/test_mvp_smoke.m");
run("tests/smoke/test_m2_smoke.m");
```

---

## ユニットテスト（`tests/unit/`）— 25 クラス

| クラス | 対象モジュール | テスト内容 |
|---|---|---|
| `TestDb` | `emk.db` | PubChem / ChEMBL 検索・フェッチ・エラーハンドリング |
| `TestDescriptor` | `emk.descriptor` | 10 種標準記述子・バッチ計算・無効 Mol の NaN |
| `TestFilter` | `emk.filter` | Lipinski Ro5 違反カウント・合否判定・エッジケース |
| `TestFingerprint` | `emk.fingerprint` | Morgan（radius / nBits）・MACCS・`toArray` の型と形状 |
| `TestInitPython` | `emk.setup.initPython` | プラットフォーム検出・冪等な `pyenv` 設定 |
| `TestInstall` | `emk.setup.install` | Embedded Python 配備・パス長ガード |
| `TestInstallExtra` | `emk.setup.installExtra` | Track 1 ライブラリのインストールとインポート検証 |
| `TestInstallOnline` | `emk.setup.installOnline` | MATLAB Online ブートストラップ・バージョン一致時スキップ |
| `TestInstallTrack2` | `emk.setup.installTrack2` | venv 作成・`settings.json` 書き込み・useExternal 呼び出し |
| `TestIo` | `emk.io` | SDF ラウンドトリップ・SMILES リスト読み込み・無効 Mol スキップ |
| `TestIsOnline` | `emk.util.isOnline` | Desktop / Online 検出ヒューリスティック |
| `TestLoadConfig` | `emkLoadConfig` | 環境変数オーバーライド・JSON 読み込み・デフォルト値 |
| `TestLogHelpers` | `logInfo / logWarn / logError / logDebug` | 出力フォーマット・verboseフラグ |
| `TestMakeRunDir` | `makeRunDir` | ディレクトリ作成・タイムスタンプ形式・プレフィックス |
| `TestMol` | `emk.mol` | fromSmiles / toSmiles / isValid / hasSubstruct / scaffold |
| `TestMordred` | `emk.descriptor.mordred*` | ~1800 記述子名・バッチ処理・NaN ハンドリング |
| `TestPubchemFetch` | `emk.db.pubchemFetch` | 拡張フェッチ・全 NameSpace 型・notFound エラー |
| `TestRdkitModule` | `emk`（全般） | RDKit 利用可能確認・バージョン文字列 |
| `TestRecipe` | `emk.setup.recipe` | 全既知ライブラリ名の出力・unknownLibrary エラー |
| `TestSimilarity` | `emk.similarity` | tanimoto / dice / rankBy（N・Metric）/ matrix の対称性 |
| `TestToTable` | `emk.mol.toTable` | 列名・型・無効 Mol の NaN・空入力 |
| `TestUseExternal` | `emk.setup.useExternal` | パス検証・fileNotFound・冪等スキップ |
| `TestValidate` | `emk.setup.validate` | テーブル構造・全ライブラリスキャン・個別名指定 |
| `TestVerify` | `emk.setup.verify` | struct フィールド・logical 値・例外なし保証 |
| `TestViz` | `emk.viz` | draw2d の Figure 返却・無効入力エラー |

### テスト設計原則

- RDKit 不要テストは、各クラス内で RDKit 必要テストより **先に** 記述する。
- RDKit 必要な全テストメソッドの先頭に以下を記述:
  ```matlab
  tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
  ```
- 値の検証に加え `verifyClass` / `verifySize` で型・形状のバグを明示的に検出する。
- テストコードで `py.rdkit.*` を直接呼び出さない — 必ず `emk.*` ラッパー経由。

---

## スモークテスト（`tests/smoke/`）— 2 スクリプト

| スクリプト | 目的 | 主なステップ |
|---|---|---|
| `test_mvp_smoke.m` | コア機能のエンドツーエンド検証 | `install` → `fromSmiles` → `calculate` → `tanimoto` → `draw2d` |
| `test_m2_smoke.m` | Milestone 2 ワークフロー検証 | `readSdf` / `readSmilesList` → `batchCalculate` → `lipinski` → `searchPubchem` / `searchChembl` |
