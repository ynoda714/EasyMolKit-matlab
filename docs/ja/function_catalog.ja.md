# 関数カタログ — EasyMolKit

> 全 76 公開関数のコンパクトリファレンス。
> 詳細シグネチャ・エラー ID → [function_reference.ja.md](function_reference.ja.md)
> 英語版 → [../en/function_catalog.md](../en/function_catalog.md)

---

## emk.setup（12 関数）

| 関数 | 概要 |
|---|---|
| `install(PythonVersion="3.10")` | Embedded Python + RDKit を `python_env/` に配備（Desktop のみ） |
| `installOnline(Config=struct())` | MATLAB Online で pip をブートストラップして RDKit をインストール |
| `initPython()` | プラットフォームを検出して `pyenv` を設定（冪等） |
| `verify()` | `.python` `.rdkit` `.version` を持つ struct を返す; 例外なし |
| `installExtra(name)` | Track 1 追加ライブラリを Embedded Python に pip インストール |
| `useExternal(pythonPath)` | `pyenv` を外部 CPython に切り替える（Track 2; Python ロード前に呼ぶ） |
| `validate(Libraries=string.empty)` | 全または指定ライブラリのインストール状況テーブルを返す; 例外なし |
| `recipe(name)` | 指定ライブラリのインストール手順・ライセンス情報をコマンドウィンドウに表示 |
| `installTrack2(name, BasePython="")` | Track 2 ライブラリの venv を作成して `settings.json` に書き込む |
| `snapshot()` | RF02 環境スナップショット struct を返す（matlab/python/rdkit/toolboxes/timestamp） |
| `lockfile(snap, filePath)` | スナップショット struct を RF02 ロック JSON に保存 |
| `lockfile(filePath)` | RF02 ロック JSON を読み込んで struct を返す |
| `verifyLock(lockRef)` | 現在の環境と保存済みロックを比較 → `.pass` `.details` `.warnings` |

---

## emk.mol（8 関数）

| 関数 | 概要 |
|---|---|
| `fromSmiles(smiles)` | SMILES 文字列 → `py.rdkit.Chem.Mol` |
| `toSmiles(mol)` | Mol オブジェクトから正規化 SMILES 文字列を返す |
| `isValid(smiles)` | `logical` を返す; 例外なし（`rdkitError` のみ再スロー） |
| `hasSubstruct(mol, query)` | SMARTS/Mol 部分構造マッチ; cell 入力 → `logical(1,N)` |
| `toStruct(mol, Format="molblock")` | Mol → MATLAB struct（molblock または pickle 形式） |
| `fromStruct(s)` | `toStruct` の逆変換; molblock/pickle を自動検出 |
| `toTable(mols, Properties=…)` | Mol の cell → table（SMILES + 指定記述子） |
| `scaffold(mol)` | Bemis-Murcko スキャフォールドを Mol として返す |

---

## emk.scaffold（3 関数）

| 関数 | 概要 |
|---|---|
| `genericMurcko(mol)` | Generic Murcko スキャフォールド（全原子 C・全結合単結合）→ Mol |
| `brics(mol)` | BRICS フラグメント分解 → `string(1,N)` SMILES 配列（非決定論的順序） |
| `rgroup(mols, coreSmiles)` | R グループ分解 → `[table, unmatchedIdx]`（Core, R1, R2… 列） |

---

## emk.dataset（4 関数）

| 関数 | 概要 |
|---|---|
| `esol(CacheDir="", ForceDownload=false)` | ESOL 水溶性データセット → table（約 1128 化合物、`logS` 列） |
| `freesolv(CacheDir="", ForceDownload=false)` | FreeSolv 溶媒和自由エネルギー → table（約 643 化合物） |
| `bbbp(CacheDir="", ForceDownload=false)` | BBBP 血液脳関門透過性 → table（約 2050 化合物、`BBB` logical） |
| `tox21(CacheDir="", ForceDownload=false)` | Tox21 毒性 → table（約 7831 化合物、12 double エンドポイント列） |

---

## emk.descriptor（10 関数）

| 関数 | 概要 |
|---|---|
| `molWeight(mol)` | 平均分子量 [g/mol] |
| `calculate(mol, descriptorNames=<all10>)` | 最大 11 種の標準記述子を計算 → struct |
| `batchCalculate(mols, descriptorNames=<all10>)` | 複数 Mol → table（無効行は NaN） |
| `mordred(mol, descriptorNames=[])` | Mordred 2D 記述子約 1800 種 → struct; `installExtra("mordred")` が必要 |
| `mordredBatch(mols, descriptorNames=[])` | バッチ Mordred 計算 → table; `installExtra("mordred")` が必要 |
| `mordredNames()` | 利用可能な Mordred 記述子名の `string(1×N)` リストを返す |
| `qed(mol)` | QED スコア ∈ [0,1] — 薬物様性複合スコア（8 物性デザイラビリティ関数） |
| `saScore(mol)` | SA スコア ∈ [1,10] — 合成容易性（低いほど合成しやすい） |
| `bcut(mol)` | BCUT2D 記述子 → `double(1,8)` （MWHI/MWLOW/CHGHI/CHGLO/LOGPHI/LOGPLOW/MRHI/MRLOW） |
| `fragmentCount(mol)` | 環系・官能基フラグメント数 → struct（11 フィールド） |

**標準記述子**（`all10`）: `MolWt`, `ExactMolWt`, `LogP`, `TPSA`, `NumHAcceptors`, `NumHDonors`,
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
| `rankBy(queryFp, dbFps, N=Inf, Metric="tanimoto")` | DB フィンガープリントを類似度でランキング; `.Indices` `.Scores` `.Metric` を返す |
| `matrix(fps, Metric="tanimoto")` | N×N 対称ペアワイズ類似度行列 |

---

## emk.filter（4 関数）

| 関数 | 概要 |
|---|---|
| `lipinski(tbl, MaxViolations=0)` | `Violations_Ro5`（double）と `Pass_Ro5`（logical）列を入力テーブルに追加 |
| `veber(tbl)` | `Violations_Veber`（double, 0-2）と `Pass_Veber`（logical）を追加 — 純 MATLAB |
| `pains(tbl)` | `NumPainsAlerts`、`PainsAlerts`、`HasPains`（logical）を追加 — RDKit FilterCatalog |
| `reos(tbl)` | `Violations_REOS`（double, 0-6）と `Pass_REOS`（logical）を追加 — 純 MATLAB |

Lipinski 必要列: `MolWt`, `LogP`, `NumHDonors`, `NumHAcceptors`
Veber 必要列: `NumRotatableBonds`, `TPSA`
PAINS 必要列: `SMILES`
REOS 必要列: `MolWt`, `LogP`, `NumHDonors`, `NumHAcceptors`, `NumRotatableBonds`, `HeavyAtomCount`

---

## emk.cluster（1 関数）

| 関数 | 概要 |
|---|---|
| `butina(fps, Threshold=0.2, Metric="tanimoto")` | Butina 球面排除クラスタリング → `cell(1,C)`（各要素は 1 始まりインデックス `double(1,K)`） |

---

## emk.diversity（1 関数）

| 関数 | 概要 |
|---|---|
| `pick(fps, N, Metric="tanimoto", Seed=-1)` | MaxMin 多様性選択 → `double(1,N)` 1 始まりインデックス |

---

## emk.conformer（2 関数）

| 関数 | 概要 |
|---|---|
| `embed(mol, Method="ETKDGv3", RandomSeed=-1)` | ETKDG 3D コンフォーマー生成 → コンフォーマー付き Mol |
| `optimize(mol, ForceField="MMFF94", MaxIter=2000)` | 力場最小化（MMFF94/UFF）→ 3D 座標最適化済み Mol |

---

## emk.shape（1 関数）

| 関数 | 概要 |
|---|---|
| `compare(mol1, mol2, Method="protrude")` | 3D 形状類似度スコア ∈ [0,1]（protrude/tanimoto; 3D コンフォーマーが必要） |

---

## emk.io（3 関数）

| 関数 | 概要 |
|---|---|
| `readSdf(filePath)` | SDF ファイル → Mol オブジェクトの `1×N cell`（失敗分子は警告付きスキップ） |
| `writeSdf(mols, filePath)` | Mol cell 配列を SDF ファイルに書き込む |
| `readSmilesList(filePath)` | 1 行 1 SMILES ファイル → Mol オブジェクトの `1×N cell` |

---

## emk.viz（1 関数）

| 関数 | 概要 |
|---|---|
| `draw2d(mol, Title="", Width=300, Height=300)` | 2D 構造を MATLAB フィギュアウィンドウに描画 |

---

## emk.db（5 関数）

| 関数 | 概要 |
|---|---|
| `searchPubchem(query, Type="name")` | PubChem 検索 → 5 列テーブル; Python 不要 |
| `pubchemFetch(identifier, NameSpace="name", MaxSynonyms=10)` | PubChemPy による拡張 PubChem フェッチ → struct |
| `searchChembl(query, Type="name")` | ChEMBL 検索 → 8 列テーブル; Python 不要 |
| `searchChemblTarget(query, TargetType="SINGLE PROTEIN", MaxRows=10)` | ChEMBL ターゲット検索 → 4 列テーブル |
| `getChemblActivity(targetId, ActivityType="IC50", MaxRows=50, MinActivity_nM=Inf)` | ChEMBL バイオアクティビティデータ → 5 列テーブル |

`pubchemFetch` には `installExtra("pubchempy")` が必要。

---

## emk.repro（1 関数）

| 関数 | 概要 |
|---|---|
| `verify(metrics, criteria)` | RF03 数値検証 — メトリクスと受入基準を比較 → `.pass` `.details` `.report` |

---

## emk.util（2 関数）

| 関数 | 概要 |
|---|---|
| `isOnline()` | `logical` を返す — MATLAB Online で `true` |
| `benchmarkBatch(smilesList, descriptorNames=<all10>)` | スループットベンチマーク → タイミングフィールドを持つ struct |

---

## src/config（1 関数）

| 関数 | 概要 |
|---|---|
| `emkLoadConfig()` | 設定をロード（環境変数 > `settings.json` > デフォルト）→ struct |

---

## src/util（8 関数）

| 関数 | 概要 |
|---|---|
| `logInfo(msg, ...)` | `[HH:MM:SS][INFO]  メッセージ` をコンソールに出力 |
| `logWarn(msg, ...)` | `[HH:MM:SS][WARN]  メッセージ` をコンソールに出力 |
| `logError(msg, ...)` | `[HH:MM:SS][ERROR]  メッセージ` をコンソールに出力 |
| `logDebug(msg, ...)` | `EMK_LOG_VERBOSE=1` のときのみ出力 |
| `logProgress(i, n, label)` | ループ反復の進捗バー |
| `logSection(scriptId, label, layer)` | チュートリアルスクリプトの各 `%%` セクション先頭のバナーログ |
| `makeRunDir(Prefix="", BaseDir="result/runs")` | タイムスタンプ付き実行ディレクトリを作成; `mkdir` 直接使用禁止 |
| `resolveProjectRoot()` | 3 段階フォールバックでプロジェクトルートを特定して `cd` |
