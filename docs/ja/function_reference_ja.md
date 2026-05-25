# 関数リファレンス — EasyMolKit

🇬🇧 [English](../function_reference.md)

> コンパクト版関数一覧。完全なシグネチャ・オプション詳細は [function_reference.md (English)](../function_reference.md) を参照。

全 **43 関数** / 10 モジュール

---

## emk.setup — Python 環境セットアップ（8 関数）

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `install()` | — | Embedded Python 3.10 + RDKit を自動配備（Desktop） |
| `installOnline()` | — | MATLAB Online 用 RDKit セットアップ |
| `initPython()` | — | プラットフォーム検出 → `pyenv` 設定（冪等） |
| `verify()` | `struct` | Python / RDKit 状態確認（例外なし） |
| `installExtra(name)` | — | Track 1 オプションライブラリを Embedded Python に追加 |
| `useExternal(path)` | — | Track 2 用外部 CPython に切り替え |
| `validate()` | `table` | インストール済みライブラリの診断テーブル |
| `recipe(name)` | — | 指定ライブラリのインストール手順を表示 |
| `installTrack2(name)` | — | Track 2 用 venv 作成・接続 |

**installExtra / recipe の対応 name**: `"pubchempy"` / `"mordred"` / `"biopython"` / `"meeko"` / `"vina"` / `"pdbfixer"` / `"torch"` / `"torch_geometric"` / `"transformers"` / `"datasets"`

---

## emk.mol — 分子オブジェクト操作（8 関数）

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `fromSmiles(smiles)` | `py.Mol` | SMILES → RDKit Mol オブジェクト |
| `toSmiles(mol)` | `string` | 正準 SMILES に変換 |
| `isValid(smiles)` | `logical` | SMILES 妥当性チェック（例外なし） |
| `hasSubstruct(mol, query)` | `logical` | SMARTS 部分構造マッチング |
| `toStruct(mol)` | `struct` | Mol → MATLAB 構造体（molblock / pickle） |
| `fromStruct(s)` | `py.Mol` | MATLAB 構造体 → Mol に復元 |
| `toTable(mols)` | `table` | Mol 配列 → 記述子付き table |
| `scaffold(mol)` | `py.Mol` | Bemis-Murcko スキャフォールド抽出 |

---

## emk.descriptor — 分子記述子計算（5 関数）

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `molWeight(mol)` | `double` | 平均分子量 [g/mol] |
| `calculate(mol)` | `struct` | 10 標準記述子を一括計算 |
| `batchCalculate(mols)` | `table` | 複数分子の記述子バッチ計算 |
| `mordred(mol)` | `struct` | Mordred 2D 記述子（~1800 種） |
| `mordredBatch(mols)` | `table` | Mordred バッチ計算 |
| `mordredNames()` | `string(1×N)` | 利用可能な Mordred 記述子名一覧 |

**標準 10 記述子**: `MolWt` `ExactMolWt` `LogP` `TPSA` `NumHAcceptors` `NumHDonors` `NumRotatableBonds` `RingCount` `FractionCSP3` `HeavyAtomCount`

---

## emk.fingerprint — フィンガープリント生成（3 関数）

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `morgan(mol, Radius=2, NBits=2048)` | `py.ExplicitBitVect` | Morgan (ECFP4) フィンガープリント |
| `maccs(mol)` | `py.ExplicitBitVect` | 167 ビット MACCS keys |
| `toArray(fp)` | `logical(1,N)` | フィンガープリント → MATLAB 論理配列 |

---

## emk.similarity — 類似度計算（4 関数）

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `tanimoto(fp1, fp2)` | `double` ∈ [0,1] | Tanimoto 係数 |
| `dice(fp1, fp2)` | `double` ∈ [0,1] | Dice 係数 |
| `rankBy(queryFp, dbFps, N=Inf)` | `struct` | データベースに対するランキング検索 |
| `matrix(fps)` | `double(N×N)` | N×N 対称類似度行列 |

---

## emk.filter — 分子フィルタリング（1 関数）

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `lipinski(tbl, MaxViolations=0)` | `table` | Lipinski の Ro5 フィルタ（違反数・合否を追記） |

**Ro5 閾値**: MW > 500、LogP > 5、HBD > 5、HBA > 10（各 1 違反）

---

## emk.io — ファイル入出力（3 関数）

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `readSdf(filePath)` | `1×N cell` | SDF ファイルから Mol 配列を読み込み |
| `writeSdf(mols, filePath)` | — | Mol 配列を SDF ファイルに書き出し |
| `readSmilesList(filePath)` | `1×N cell` | 1 行 1 SMILES 形式ファイルの読み込み |

---

## emk.viz — 可視化（1 関数）

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `draw2d(mol, Title="", Width=300, Height=300)` | `Figure` | 2D 構造を MATLAB Figure に描画 |

> 注: 1 分子あたり 0.5〜2 秒の描画時間が必要です（RDKit → PNG → MATLAB 転送）

---

## emk.db — 外部データベース検索（5 関数）

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `searchPubchem(query, Type="name")` | `table` | PubChem 化合物検索（Python 不要） |
| `pubchemFetch(identifier)` | `struct` | PubChemPy を使った詳細情報取得 |
| `searchChembl(query, Type="name")` | `table` | ChEMBL 化合物検索 |
| `searchChemblTarget(query)` | `table` | ChEMBL ターゲット検索 |
| `getChemblActivity(targetId)` | `table` | ChEMBL 活性データ取得 |

---

## emk.util — ユーティリティ（2 関数）

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `isOnline()` | `logical` | MATLAB Online 環境の検出 |
| `benchmarkBatch(smilesList)` | `struct` | バッチ処理スループット計測 |

---

## src/util — ログヘルパー（4 関数）

| 関数 | 説明 |
|---|---|
| `logInfo(msg, ...)` | 情報ログ出力（`fprintf` の代替） |
| `logWarn(msg, ...)` | 警告ログ出力 |
| `logError(msg, ...)` | エラーログ出力 |
| `logDebug(msg, ...)` | デバッグログ（`EMK_LOG_VERBOSE=1` 時のみ） |
| `logProgress(i, n, label)` | ループ進捗バー |
| `logSection(scriptId, label, layer)` | チュートリアルセクションバナー |

---

## src/config — 設定ロード（1 関数）

| 関数 | 戻り値 | 説明 |
|---|---|---|
| `loadConfig()` | `struct` | 設定ロード（環境変数 > settings.json > デフォルト値） |
| `makeRunDir(Prefix="")` | `char` | 実行成果物ディレクトリ作成・パス返却 |
