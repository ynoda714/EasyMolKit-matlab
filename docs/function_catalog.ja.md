# 関数カタログ — EasyMolKit

> 全モジュール 51 関数のコンパクトリファレンスです。
> 完全なシグネチャ・エラー ID → [function_reference.md](function_reference.md)
> English version → [function_catalog.md](function_catalog.md)

---

## emk.setup（9 関数）

| 関数 | 概要 |
|---|---|
| `install(PythonVersion="3.10")` | Embedded Python + RDKit を `python_env/` に自動配備（Desktop 専用） |
| `installOnline(Config=struct())` | pip ブートストラップ → MATLAB Online に RDKit をインストール |
| `initPython()` | プラットフォーム検出 → `pyenv` 設定（冪等スキップあり） |
| `verify()` | `.python` / `.rdkit` / `.version` を持つ struct を返す（例外なし） |
| `installExtra(name)` | Track 1 アドオンを Embedded Python に pip インストール |
| `useExternal(pythonPath)` | 外部 CPython に `pyenv` を切り替える（Python ロード前に呼ぶこと） |
| `validate(Libraries=string.empty)` | ライブラリのインストール状況テーブルを返す（例外なし） |
| `recipe(name)` | インストール手順とライセンス情報をコマンドウィンドウに表示 |
| `installTrack2(name, BasePython="")` | Track 2 ライブラリ用 venv を作成し `settings.json` を更新 |

---

## emk.mol（8 関数）

| 関数 | 概要 |
|---|---|
| `fromSmiles(smiles)` | SMILES 文字列 → `py.rdkit.Chem.Mol` |
| `toSmiles(mol)` | Mol → Canonical SMILES 文字列 |
| `isValid(smiles)` | `logical` を返す（`rdkitError` 以外は例外なし） |
| `hasSubstruct(mol, query)` | SMARTS / Mol で部分構造マッチ; cell 入力 → `logical(1,N)` |
| `toStruct(mol, Format="molblock")` | Mol → MATLAB struct（molblock または pickle 形式） |
| `fromStruct(s)` | `toStruct` の逆変換; molblock / pickle を自動検出 |
| `toTable(mols, Properties=…)` | Mol の cell 配列 → table（SMILES + 選択した記述子） |
| `scaffold(mol)` | Bemis-Murcko スキャフォールドを Mol で返す |

---

## emk.descriptor（6 関数）

| 関数 | 概要 |
|---|---|
| `molWeight(mol)` | 平均分子量 [g/mol] |
| `calculate(mol, descriptorNames=<all10>)` | 最大 11 種の標準記述子を計算 → struct |
| `batchCalculate(mols, descriptorNames=<all10>)` | 複数 Mol → table（無効 Mol の行は NaN） |
| `mordred(mol, descriptorNames=[])` | Mordred 2D 記述子 ~1800 種 → struct（要 `installExtra("mordred")`） |
| `mordredBatch(mols, descriptorNames=[])` | Mordred 一括計算 → table（要 `installExtra("mordred")`） |
| `mordredNames()` | 利用可能な Mordred 記述子名の昇順リストを返す |

**標準記述子（all10）**: `MolWt`, `ExactMolWt`, `LogP`, `TPSA`, `NumHAcceptors`, `NumHDonors`,
`NumRotatableBonds`, `RingCount`, `FractionCSP3`, `HeavyAtomCount`。`MolFormula` は明示指定時のみ。

---

## emk.fingerprint（3 関数）

| 関数 | 概要 |
|---|---|
| `morgan(mol, Radius=2, NBits=2048)` | Morgan（ECFP4）フィンガープリント → `py.rdkit.DataStructs.ExplicitBitVect` |
| `maccs(mol)` | 167 ビット MACCS キーフィンガープリント |
| `toArray(fp)` | フィンガープリント → `logical(1,N)` MATLAB 配列に変換 |

---

## emk.similarity（4 関数）

| 関数 | 概要 |
|---|---|
| `tanimoto(fp1, fp2)` | Tanimoto 係数 ∈ [0, 1] |
| `dice(fp1, fp2)` | Dice 係数 ∈ [0, 1] |
| `rankBy(queryFp, dbFps, N=Inf, Metric="tanimoto")` | データベース FP を類似度でランキング → `.Indices` / `.Scores` / `.Metric` |
| `matrix(fps, Metric="tanimoto")` | N×N 対称ペアワイズ類似度行列 |

---

## emk.filter（1 関数）

| 関数 | 概要 |
|---|---|
| `lipinski(tbl, MaxViolations=0)` | `Violations_Ro5`（double）と `Pass_Ro5`（logical）列を table に追加 |

必須列: `MolWt`, `LogP`, `NumHDonors`, `NumHAcceptors`

---

## emk.io（3 関数）

| 関数 | 概要 |
|---|---|
| `readSdf(filePath)` | SDF ファイル → `1×N cell` の Mol 配列（失敗分子は警告でスキップ） |
| `writeSdf(mols, filePath)` | Mol cell 配列を SDF ファイルに書き出す |
| `readSmilesList(filePath)` | 1 行 1 SMILES ファイル → `1×N cell` の Mol 配列 |

---

## emk.viz（1 関数）

| 関数 | 概要 |
|---|---|
| `draw2d(mol, Title="", Width=300, Height=300)` | 2D 構造を MATLAB figure に描画 |

---

## emk.db（5 関数）

| 関数 | 概要 |
|---|---|
| `searchPubchem(query, Type="name")` | PubChem 検索 → 5 列 table（Python 不要） |
| `pubchemFetch(identifier, NameSpace="name", MaxSynonyms=10)` | PubChemPy で拡張プロパティを取得 → struct |
| `searchChembl(query, Type="name")` | ChEMBL 検索 → 8 列 table（Python 不要） |
| `searchChemblTarget(query, TargetType="SINGLE PROTEIN", MaxRows=10)` | ChEMBL ターゲット検索 → 4 列 table |
| `getChemblActivity(targetId, ActivityType="IC50", MaxRows=50, MinActivity_nM=Inf)` | ChEMBL バイオアクティビティデータ → 5 列 table |

`pubchemFetch` は `installExtra("pubchempy")` が必要。

---

## emk.util（2 関数）

| 関数 | 概要 |
|---|---|
| `isOnline()` | MATLAB Online なら `true` を返す `logical` |
| `benchmarkBatch(smilesList, descriptorNames=<all10>)` | スループット計測 → タイミングフィールドを持つ struct |

---

## src/config（1 関数）

| 関数 | 概要 |
|---|---|
| `emkLoadConfig()` | 設定を読み込む（環境変数 > `settings.json` > デフォルト）→ struct |

---

## src/util（8 関数）

| 関数 | 概要 |
|---|---|
| `logInfo(msg, ...)` | `[HH:MM:SS][INFO]  メッセージ` をコンソールに出力 |
| `logWarn(msg, ...)` | `[HH:MM:SS][WARN]  メッセージ` をコンソールに出力 |
| `logError(msg, ...)` | `[HH:MM:SS][ERROR]  メッセージ` をコンソールに出力 |
| `logDebug(msg, ...)` | `EMK_LOG_VERBOSE=1` のときのみ出力 |
| `logProgress(i, n, label)` | ループの進捗バーを表示 |
| `logSection(scriptId, label, layer)` | チュートリアルスクリプトのセクションバナーを表示 |
| `makeRunDir(Prefix="", BaseDir="result/runs")` | タイムスタンプ付き実行ディレクトリを作成（`mkdir` 直書き禁止） |
| `resolveProjectRoot()` | 3 段階フォールバックでプロジェクトルートを特定し `cd` する |
