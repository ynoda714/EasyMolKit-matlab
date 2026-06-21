# テストカタログ — EasyMolKit

> ユニットテスト 30 クラス + スモークテスト 2 スクリプト。
> 英語版 → [../en/test_catalog.md](../en/test_catalog.md)

## テストの実行

```matlab
addpath(genpath("src"));
addpath(genpath("tests"));   % クラス名ルックアップに必要

% 全ユニットテストを実行
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

## ユニットテスト（`tests/unit/`）— 30 クラス

| クラス | 対象モジュール | テスト範囲 |
|---|---|---|
| `TestCluster` | `emk.cluster` | Butina 球面排除クラスタリング・閾値・バリデーション |
| `TestConformer` | `emk.conformer` / `emk.shape` | embed/optimize/compare: バリデーション・統合・シード再現性 |
| `TestDataset` | `emk.dataset` | ESOL / FreeSolv / BBBP / Tox21 バリデーション・ダウンロードガード |
| `TestDb` | `emk.db` | PubChem / ChEMBL 検索・フェッチ・エラーハンドリング |
| `TestDescriptor` | `emk.descriptor` | 標準 10 記述子 + qed / saScore / bcut(1×8) / fragmentCount; バッチ・無効分子の NaN |
| `TestDiversity` | `emk.diversity` | MaxMin 多様性選択・シード再現性・バリデーション |
| `TestFilter` | `emk.filter` | Lipinski / Veber / REOS / PAINS: 違反数・合否・エッジケース |
| `TestFingerprint` | `emk.fingerprint` | Morgan（半径・ビット数）、MACCS、`toArray` 型/形状 |
| `TestInitPython` | `emk.setup.initPython` | プラットフォーム検出・冪等な `pyenv` 設定 |
| `TestInstall` | `emk.setup.install` | Embedded Python 配備・パス長ガード |
| `TestInstallExtra` | `emk.setup.installExtra` | Track 1 ライブラリのインストールとインポート検証 |
| `TestInstallOnline` | `emk.setup.installOnline` | MATLAB Online ブートストラップ・既存バージョンのスキップロジック |
| `TestInstallTrack2` | `emk.setup.installTrack2` | venv 作成・`settings.json` 書き込み・useExternal 呼び出し |
| `TestIo` | `emk.io` | SDF ラウンドトリップ・SMILES リスト読み込み・無効分子スキップ |
| `TestIsOnline` | `emk.util.isOnline` | Desktop / Online 検出ヒューリスティック |
| `TestLoadConfig` | `emkLoadConfig` | 環境変数オーバーライド・JSON ロード・デフォルト値 |
| `TestLogHelpers` | `logInfo / logWarn / logError / logDebug` | 出力フォーマット・verbose フラグ |
| `TestMakeRunDir` | `makeRunDir` | ディレクトリ作成・タイムスタンプ形式・プレフィックス |
| `TestMol` | `emk.mol` | fromSmiles, toSmiles, isValid, hasSubstruct, scaffold |
| `TestMordred` | `emk.descriptor.mordred*` | 約 1800 記述子名・バッチ・NaN ハンドリング |
| `TestPubchemFetch` | `emk.db.pubchemFetch` | 拡張フェッチ・全 NameSpace タイプ・notFound エラー |
| `TestRdkitModule` | `emk`（全般） | RDKit 利用可否チェック・バージョン文字列 |
| `TestRecipe` | `emk.setup.recipe` | 既知ライブラリ名の出力・unknownLibrary エラー |
| `TestScaffold` | `emk.scaffold` | genericMurcko / brics / rgroup: バリデーション・SMARTS マッチング |
| `TestSimilarity` | `emk.similarity` | tanimoto, dice, rankBy（N, Metric）、matrix 対称性 |
| `TestToTable` | `emk.mol.toTable` | 列名・型・無効分子の NaN・空入力 |
| `TestUseExternal` | `emk.setup.useExternal` | パスバリデーション・fileNotFound・冪等スキップ |
| `TestValidate` | `emk.setup.validate` | テーブル構造・全ライブラリスキャン・指定名 |
| `TestVerify` | `emk.setup.verify` | 構造体フィールド・論理値・非スロー保証 |
| `TestViz` | `emk.viz` | draw2d が Figure を返す・無効入力エラー |

### テスト規約

- RDKit が不要なテストは各クラス内で **RDKit 必要テストより前に**配置する。
- RDKit 必要テストのメソッドは必ず以下で始める:
  ```matlab
  tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
  ```
- 値チェックに加えて `verifyClass` / `verifySize` で型・形状バグを明示的に検出する。
- テストコードでも `py.rdkit.*` を直接呼び出さない — 必ず `emk.*` ラッパー経由。

---

## スモークテスト（`tests/smoke/`）— 2 スクリプト

| スクリプト | 目的 | 主要ステップ |
|---|---|---|
| `test_mvp_smoke.m` | コアのエンドツーエンドワークフロー | `install` → `fromSmiles` → `calculate` → `tanimoto` → `draw2d` |
| `test_m2_smoke.m` | Milestone 2 ワークフロー | `readSdf` / `readSmilesList` → `batchCalculate` → `lipinski` → `searchPubchem` / `searchChembl` |
