# 関数リファレンス — EasyMolKit

> 関数シグネチャとオプションの完全リスト。開発時の主要参照。
> アルゴリズム根拠 → [../algorithm_guide.md](../algorithm_guide.md)
> 英語版 → [../en/function_reference.md](../en/function_reference.md)

```
src/+emk/+setup/       src/+emk/+mol/        src/+emk/+descriptor/
src/+emk/+fingerprint/ src/+emk/+similarity/  src/+emk/+filter/
src/+emk/+scaffold/    src/+emk/+dataset/     src/+emk/+cluster/
src/+emk/+diversity/   src/+emk/+conformer/   src/+emk/+shape/
src/+emk/+io/          src/+emk/+viz/         src/+emk/+util/
src/+emk/+repro/       src/config/            src/util/
```

---

## emk.setup

### `emk.setup.install(PythonVersion="3.10")`
Desktop 用に Embedded Python + RDKit を `python_env/` に自動配備する。
パス長チェック → zip 展開 → pip インストール → 検証を 1 回の呼び出しで実行。
エラー ID: `notDesktop`（Online 環境）、`pathTooLong`（240 文字超）、`downloadFailed`、`installFailed`。

### `emk.setup.installOnline(Config=struct())`
MATLAB Online 用。`get-pip.py` → `pip install rdkit-pypi==<ver> --user` → `py.sys.path` 挿入 → 検証を実行。
既に一致するバージョンが検出された場合はスキップ。
`Config` を渡すと `cfg.optionalLibraries.<name>=true` のライブラリを追加で pip インストール（Desktop の `install(Config=cfg)` と対称な API）。
エラー ID: `notOnline`（Desktop 環境）、`bootstrapFailed`、`installFailed`。

### `emk.setup.initPython()`
プラットフォームを検出 → `pyenv(Version=..., ExecutionMode="OutOfProcess")` を呼び出す。
冪等: Python がすでに Loaded の場合は何もしない。
エラー ID: `notInstalled`（Desktop + `python_env/python.exe` なし）、`pyenvFailed`。

### `emk.setup.verify()`
戻り値: `struct` — `.python`（logical）、`.rdkit`（logical）、`.version`（string）。
例外なし（non-throwing 設計）。

### `emk.setup.installExtra(name)`
Track 1 オプションライブラリを Embedded Python（`python_env/`）に pip インストールする。
`name`: `"pubchempy"` / `"mordred"` / `"biopython"` / `"meeko"` / `"vina"` / `"pdbfixer"` /
`"torch"` / `"torch_geometric"` / `"transformers"`。
バージョンは `config/settings.json` の `extraLibraries.<name>` で管理。インストール後に `import` で検証。
`"torch_geometric"` は特殊フロー: インストール済み torch バージョンを自動検出し、PyG ホイール URL
（`https://data.pyg.org/whl/torch-<X.Y.Z>+cpu.html`）を動的に構築し、
`torch_scatter` / `torch_sparse` / `torch_cluster` も併せてインストール。
`"pdbfixer"` は openmm >= 8.2 を依存として自動取得（約 70 MB）。
エラー ID: `notDesktop`（Online 環境）、`unknownLibrary`、`installFailed`、`importVerifyFailed`。

### `emk.setup.useExternal(pythonPath)`
Track 2 用。`pyenv` を `pythonPath` の外部 CPython 実行ファイルに切り替える。
セッション内で冪等（Python がすでに Loaded の場合は警告のみ）。
`pythonPath`: string または char（絶対パス）。呼び出し後に `validate()` を推奨。
エラー ID: `invalidInput`（非 string/char または空）、`fileNotFound`、`pyenvFailed`。

### `emk.setup.validate(Libraries=string.empty)`
Python の `importlib.metadata` を使って指定ライブラリの存在とバージョンを診断する。
`Libraries` 省略時は RDKit + 現在インストール済みの全 Track 1 ライブラリを確認。
戻り値: `table` — 列: `Library`（string）、`Installed`（logical）、`Version`（string）、`Track`（string）。
例外なし（non-throwing 設計）。
エラー ID: `invalidInput`（非 string 配列）。

### `emk.setup.recipe(name)`
指定ライブラリのインストールレシピ（トラック分類・コマンド・ライセンス・備考）を
MATLAB コマンドウィンドウに表示する。
`name`: `"pubchempy"` / `"mordred"` / `"biopython"` / `"meeko"` / `"vina"` / `"pdbfixer"` /
`"docking"`（meeko+vina+pdbfixer 一括）/ `"openbabel"` / `"mdanalysis"` / `"pymol"`。
戻り値なし（表示のみ）。
エラー ID: `unknownLibrary`。

### `emk.setup.installTrack2(name, BasePython="")`
Track 2 ライブラリ用に `python_env_t2/<name>/` に venv を作成し、pip インストールを実行、
`settings.json` に `python.external_path` を書き込み、`useExternal()` で接続する（ADR-007）。
`name`: `"mdanalysis"`（GPLv2+）/ `"pymol"`（PSF/BSD）。Open Babel は MSI 必須のため除外（`recipe("openbabel")` 参照）。
`BasePython`: 省略時は PATH から `py` / `python` を自動検出。
次回セッション以降、`initPython()` が `external_path` を自動検出して `useExternal()` を呼び出す。
エラー ID: `notDesktop`、`unknownLibrary`、`basePythonNotFound`、`venvFailed`、`installFailed`、`importVerifyFailed`、`settingsWriteFailed`。

### `emk.setup.snapshot()`
戻り値: `struct` — `.matlab` `.python` `.rdkit` `.commit` `.toolboxes` `.timestamp`。
RF02 バージョンロック用に現在の環境をキャプチャする。例外なし。
関連: `emk.setup.lockfile`、`emk.setup.verifyLock`。

### `emk.setup.lockfile(snap, filePath)` / `emk.setup.lockfile(filePath)`
保存モード（2 引数）: `snap` struct を整形 JSON として `filePath` に保存。
読み込みモード（1 string 引数）: ロック JSON ファイルを読み込んで struct を返す。
エラー ID: `invalidInput`、`writeError`、`fileNotFound`。

### `emk.setup.verifyLock(lockRef)`
現在の環境を RF02 ロックと比較する。
`lockRef`: JSON ファイルパス（string/char）またはロック struct。
戻り値: `struct` — `.pass`（logical）、`.details`（フィールドごと: `.expected .actual .match`）、`.warnings`（string）。
重要フィールド: `matlab`、`python`、`rdkit`。非重要（警告のみ）: `commit`。
エラー ID: `invalidInput`。

---

## emk.mol

### `emk.mol.fromSmiles(smiles)`
戻り値: `py.rdkit.Chem.rdchem.Mol` — RDKit Mol オブジェクト（Python 参照、ADR-002）。入力: `string|char`。
エラー ID: `invalidInput`（型エラー）、`invalidSmiles`（空または RDKit が None を返す）。

### `emk.mol.toSmiles(mol)`
戻り値: `string` — 正規化 SMILES。Morgan アルゴリズムで正規化; 入力 SMILES と異なる場合がある。
エラー ID: `invalidInput`（非 Mol オブジェクト）、`rdkitError`。

### `emk.mol.isValid(smiles)`
戻り値: `logical` — SMILES が有効なら `true`、無効なら `false`。例外なし（`rdkitError` のみ再スロー）。
エラー ID: `invalidInput`（非 string/char）、`rdkitError`（RDKit が未起動の場合のみ）。

### `emk.mol.hasSubstruct(mol, query)`
戻り値: `logical` または `logical(1,N)` — 部分構造マッチ結果。`query` は SMARTS 文字列または Mol。
`mol` が cell 配列の場合は 1×N logical ベクトルを返す。
エラー ID: `invalidMol`、`invalidQuery`（無効な SMARTS）、`rdkitError`。

### `emk.mol.toStruct(mol, Format="molblock")`
戻り値: `struct` — `.smiles` `.formula` `.numAtoms` `.numBonds` `.molblock`（または `.pickle`）。
高忠実度な保存には `Format="pickle"` を使用。
エラー ID: `invalidInput`、`rdkitError`。

### `emk.mol.fromStruct(s)`
戻り値: `py.rdkit.Chem.Mol` — `toStruct` の逆変換。molblock/pickle フィールドを自動検出して復元。
エラー ID: `invalidInput`（struct でない、または必須フィールドなし）、`rdkitError`。

### `emk.mol.toTable(mols, Properties=["SMILES","MolWt",...all10])`
戻り値: `table`（N×M）— SMILES（string）+ 記述子（double）。
無効 Mol 行は SMILES=`"<invalid>"` で記述子=NaN。
エラー ID: `invalidInput`、`unknownProperty`、`emptyProperties`。

### `emk.mol.scaffold(mol)`
戻り値: `py.rdkit.Chem.rdchem.Mol` — Bemis-Murcko スキャフォールド。非環状分子は 0 原子 Mol を返す（例外なし）。
エラー ID: `invalidInput`（非 Mol）、`rdkitError`。

---

## emk.descriptor

### `emk.descriptor.molWeight(mol)`
戻り値: `double` — 平均分子量 [g/mol]。`Descriptors.MolWt` 使用
（全原子の IUPAC 平均原子量の合計。implicit H を含む）。
エラー ID: `invalidInput`（非 Mol）。

### `emk.descriptor.calculate(mol, descriptorNames=<all10>)`
戻り値: `struct` — フィールド名=記述子名、値=double（`MolFormula` のみ string）。
デフォルト（引数なし）で以下の 10 種の数値記述子を計算。`MolFormula` は明示指定時のみ返す。

| 名前 | 説明 | RDKit API | 型 |
|---|---|---|---|
| `MolWt` | 平均分子量（g/mol） | `Descriptors.MolWt` | double |
| `ExactMolWt` | 同位体正確分子量 | `Descriptors.ExactMolWt` | double |
| `LogP` | Wildman-Crippen LogP | `Descriptors.MolLogP` | double |
| `TPSA` | 位相極性表面積（Å²） | `Descriptors.TPSA` | double |
| `NumHAcceptors` | 水素結合受容体数 | `rdMolDescriptors.CalcNumHBA` | double |
| `NumHDonors` | 水素結合供与体数 | `rdMolDescriptors.CalcNumHBD` | double |
| `NumRotatableBonds` | 回転可能結合数 | `rdMolDescriptors.CalcNumRotatableBonds` | double |
| `RingCount` | 全環数（SSSR） | `rdMolDescriptors.CalcNumRings` | double |
| `FractionCSP3` | sp3 炭素比率（Fsp3） | `Descriptors.FractionCSP3` | double |
| `HeavyAtomCount` | 重原子数 | `mol.GetNumHeavyAtoms()` | double |
| `MolFormula` | 分子式文字列（例: "C2H6O"） | `rdMolDescriptors.CalcMolFormula` | **string** |

エラー ID: `invalidInput`、`unknownDescriptor`。

### `emk.descriptor.batchCalculate(mols, descriptorNames=<all10>)`
戻り値: `table`（N×M）— 複数分子のバッチ記述子計算。
列名=記述子名、値=double。無効 Mol 行は NaN。
エラー ID: `invalidInput`（非 cell）、`unknownDescriptor`、`allMolsFailed`。

### `emk.descriptor.mordred(mol, descriptorNames=[])`
Mordred 2D 記述子を計算して struct を返す。`descriptorNames` 省略時は約 1800 種を計算。
失敗した記述子は NaN。
必要条件: `emk.setup.installExtra("mordred")` が必要。
エラー ID: `invalidInput`（非 Mol）、`libraryNotFound`、`pythonError`。

### `emk.descriptor.mordredBatch(mols, descriptorNames=[])`
戻り値: `table`（N×M）— バッチ Mordred 記述子計算。`run_mordred.py` ヘルパーで
1 回の IPC ラウンドトリップを実現（ADR-002 rev.3）。無効 Mol 行は NaN。
必要条件: `emk.setup.installExtra("mordred")` が必要。
エラー ID: `invalidInput`（非 cell）、`libraryNotFound`、`allMolsFailed`、`pythonError`。

### `emk.descriptor.mordredNames()`
戻り値: `string(1×N)` — 利用可能な Mordred 2D 記述子名のソート済みリスト（約 1800 件）。
`mordred()` / `mordredBatch()` の `descriptorNames` 引数として使用。
必要条件: `emk.setup.installExtra("mordred")` が必要。
エラー ID: `libraryNotFound`、`pythonError`。

### `emk.descriptor.qed(mol)`
戻り値: `double` ∈ [0,1] — QED（薬物様性定量推定）スコア。8 種の分子物性を
デザイラビリティ関数で合成したスコア。スコア ≥ 0.67 が薬物様性の一般的な閾値。
`rdkit.Chem.QED.qed` 使用。
エラー ID: `invalidInput`（非 Mol）、`rdkitError`。

### `emk.descriptor.saScore(mol)`
戻り値: `double` ∈ [1,10] — SA（合成容易性）スコア。低いほど合成しやすい
（1-3: 容易、3-6: 中程度、6-10: 困難/不可能）。フラグメント寄与と環/立体複雑度ペナルティを組み合わせる。
必要条件: `rdkit.Contrib.SA_Score.sascorer`（標準 RDKit に含まれる）。
エラー ID: `invalidInput`、`rdkitError`（モジュール未検出を含む）。

### `emk.descriptor.bcut(mol)`
戻り値: `double(1,8)` — BCUT2D 記述子。MW・部分電荷・LogP・モル屈折率で重み付けした
Burden 行列の固有値。順序: `BCUT2D_MWHI`、`BCUT2D_MWLOW`、`BCUT2D_CHGHI`、`BCUT2D_CHGLO`、
`BCUT2D_LOGPHI`、`BCUT2D_LOGPLOW`、`BCUT2D_MRHI`、`BCUT2D_MRLOW`。`rdMolDescriptors.BCUT2D` 使用。
エラー ID: `invalidInput`、`rdkitError`。

### `emk.descriptor.fragmentCount(mol)`
戻り値: `struct` — 環系・官能基フラグメント数。フィールド: `NumRings`、
`NumAromaticRings`、`NumAliphaticRings`、`NumHeteroRings`、`NumCarbonyl`、`NumAmine`、
`NumHydroxyl`、`NumHalogen`、`NumNitrile`、`NumSulfonamide`、`NumAmide`。
`rdMolDescriptors` + `Fragments`（fr_* SMARTS）使用。
エラー ID: `invalidInput`、`rdkitError`。

---

## emk.scaffold

### `emk.scaffold.genericMurcko(mol)`
戻り値: `py.rdkit.Chem.rdchem.Mol` — Generic Murcko スキャフォールド（全原子を C に設定、全結合を単結合に設定）。
2 ステップ: `GetScaffoldForMol` → `MakeScaffoldGeneric`。非環状分子はエラーなしで 0 原子 Mol を返す。
注意: 原子/結合タイプを保持する `emk.mol.scaffold` とは異なる。
エラー ID: `invalidInput`、`rdkitError`。

### `emk.scaffold.brics(mol)`
戻り値: `string(1,N)` — BRICS フラグメント SMILES 配列（非決定論的順序; frozenset 由来）。
`rdkit.Chem.BRICS.BRICSDecompose` 使用。
エラー ID: `invalidInput`、`rdkitError`。

### `emk.scaffold.rgroup(mols, coreSmiles)`
戻り値: `[tbl, unmatchedIdx]` — R グループ分解結果。
`tbl`: 列 `Core`、`R1`、`R2`、...（各 R グループの string SMILES）を持つ `table`。
`unmatchedIdx`: コアにマッチしなかった分子の 1 始まりインデックス `double(1,K)`。
`coreSmiles` には `[*:1]` スタイルの R マップドコアを推奨。
エラー ID: `invalidInput`、`invalidCore`、`noMatch`（全分子アンマッチ）、`rdkitError`。

---

## emk.dataset

### `emk.dataset.esol(CacheDir="", ForceDownload=false)`
戻り値: `table` — ESOL（Delaney）水溶性データセット（約 1128 化合物）。列: `SMILES`、
`Name`（string）、`logS`、`logS_Delaney`、`MolWt`（double）。DeepChem S3 から `delaney-processed.csv` を
ダウンロード; `data/benchmark/esol.csv` にキャッシュ。
エラー ID: `invalidInput`、`downloadFailed`、`parseFailed`。

### `emk.dataset.freesolv(CacheDir="", ForceDownload=false)`
戻り値: `table` — FreeSolv 水和自由エネルギーデータセット（約 643 化合物）。列: `SMILES`、
`Name`（string）、`DeltaG_exp`、`DeltaG_calc`、`DeltaG_exp_sem`（double）。
エラー ID: `invalidInput`、`downloadFailed`、`parseFailed`。

### `emk.dataset.bbbp(CacheDir="", ForceDownload=false)`
戻り値: `table` — BBBP 血液脳関門透過性データセット（約 2050 化合物）。列: `SMILES`、
`Name`（string）、`BBB`（logical、true=透過性あり）。
エラー ID: `invalidInput`、`downloadFailed`、`parseFailed`。

### `emk.dataset.tox21(CacheDir="", ForceDownload=false)`
戻り値: `table` — Tox21 多標的毒性データセット（約 7831 化合物）。列: `SMILES`、
`MolID`（string）+ 12 毒性エンドポイント（double、0/1/NaN）。`.csv.gz` をダウンロードして `gunzip` で解凍。
エンドポイント: `NR_AR`、`NR_AR_LBD`、`NR_AhR`、`NR_Aromatase`、`NR_ER`、`NR_ER_LBD`、
`NR_PPAR_gamma`、`SR_ARE`、`SR_ATAD5`、`SR_HSE`、`SR_MMP`、`SR_p53`。
エラー ID: `invalidInput`、`downloadFailed`、`parseFailed`。

---

## emk.filter

### `emk.filter.lipinski(tbl, MaxViolations=0)`
戻り値: `table` — 入力テーブルに `Violations_Ro5`（double）と `Pass_Ro5`（logical）列を追加。
必要列: `MolWt`、`LogP`、`NumHDonors`、`NumHAcceptors`。

Ro5 閾値: MW>500、LogP>5、HBD>5、HBA>10 の各違反を 1 としてカウント。
`MaxViolations` の有効範囲は [0,4]。
> NaN 記述子は `NaN > 500 = false` と評価されるため「違反なし」扱いになる（偽陰性リスク）。
> 前処理として `rmmissing(tbl)` を推奨。
エラー ID: `invalidInput`（非テーブル）、`invalidMaxViol`（[0,4] 範囲外 / NaN / Inf）、`missingColumns`。

### `emk.filter.veber(tbl)`
戻り値: `table` — 入力テーブルに `Violations_Veber`（double, 0-2）と `Pass_Veber`（logical）を追加。
必要列: `NumRotatableBonds`、`TPSA`。純 MATLAB; RDKit 不要。
Veber 基準: `NumRotatableBonds > 10` または `TPSA > 140 Å²` を各 1 違反としてカウント。
エラー ID: `invalidInput`、`missingColumns`。

### `emk.filter.pains(tbl)`
戻り値: `table` — 入力テーブルに `NumPainsAlerts`（double）、`PainsAlerts`（string, カンマ区切り）、
`HasPains`（logical）列を追加。必要列: `SMILES`。
RDKit `FilterCatalog(PAINS)` 使用。
エラー ID: `invalidInput`、`missingColumns`、`rdkitError`。

### `emk.filter.reos(tbl)`
戻り値: `table` — 入力テーブルに `Violations_REOS`（double, 0-6）と `Pass_REOS`（logical）を追加。
必要列: `MolWt`、`LogP`、`NumHDonors`、`NumHAcceptors`、`NumRotatableBonds`、`HeavyAtomCount`。
純 MATLAB; RDKit 不要。
REOS 範囲（6/7 基準; FormalCharge は除外）: MW [200,500]、LogP [-5,5]、HBD [0,5]、
HBA [0,10]、RotBonds [0,8]、HeavyAtoms [15,50]。
エラー ID: `invalidInput`、`missingColumns`。

---

## emk.cluster

### `emk.cluster.butina(fps, Threshold=0.2, Metric="tanimoto")`
戻り値: `cell(1,C)` — C クラスター。各要素は `fps` への 1 始まりインデックスの `double(1,K)` 配列。
`clusters{1}` が最大クラスター; 各クラスターの最初のインデックスがセントロイド。
Butina 球面排除クラスタリング。`Threshold` は Tanimoto **距離**閾値
（デフォルト 0.2 は 80% 類似度カットオフに相当）。
`Metric`: `"tanimoto"` のみ。
エラー ID: `invalidInput`、`invalidThreshold`（(0,1] 範囲外）、`invalidMetric`、`rdkitError`。

---

## emk.diversity

### `emk.diversity.pick(fps, N, Metric="tanimoto", Seed=-1)`
戻り値: `double(1,N)` — `fps` から選択した N 個の最大多様性分子の 1 始まりインデックス。
MaxMin（Kennard-Stone 変形）アルゴリズム（`MaxMinPicker.LazyBitVectorPick`）使用。
`Seed=-1`: ランダム開始; `Seed≥0`: 再現可能。
`Metric`: `"tanimoto"` のみ。
エラー ID: `invalidInput`、`invalidN`（[1,M] 範囲外）、`invalidMetric`、`rdkitError`。

---

## emk.conformer

### `emk.conformer.embed(mol, Method="ETKDGv3", RandomSeed=-1)`
戻り値: `py.rdkit.Chem.rdchem.Mol` — 3D コンフォーマー付き分子（水素除去済み）。
パイプライン: `AddHs` → `EmbedMolecule(params)` → `RemoveHs`。
`Method`: `"ETKDGv3"`（推奨）、`"ETKDGv2"`、`"ETKDG"`、`"KDG"`。
`RandomSeed=-1`: ランダム; `RandomSeed≥0`: 再現可能。
エラー ID: `invalidInput`、`invalidMethod`、`embeddingFailed`（RDKit が -1 を返した）、`rdkitError`。

### `emk.conformer.optimize(mol, ForceField="MMFF94", MaxIter=2000)`
戻り値: `py.rdkit.Chem.rdchem.Mol` — 力場最小化済み分子。
パイプライン: `AddHs(addCoords=true)` → 最小化 → `RemoveHs`。
`ForceField`: `"MMFF94"`（薬物様分子に推奨）または `"UFF"`（特殊元素のフォールバック）。
非収束（status=1）は警告のみ; FF セットアップ失敗（status=-1）は `optimizeFailed`。
エラー ID: `invalidInput`（非 Mol またはコンフォーマーなし）、`invalidForceField`、`optimizeFailed`、`rdkitError`。

---

## emk.shape

### `emk.shape.compare(mol1, mol2, Method="protrude")`
戻り値: `double` ∈ [0,1] — 3D 形状類似度スコア（高いほど類似）。両分子に
3D コンフォーマー（`emk.conformer.embed` による）が必要。
`Method="protrude"`: `1 - ShapeProtrudeDist(mol1, mol2, allowReordering=false)`。
`Method="tanimoto"`: `ShapeTverskyIndex(mol1, mol2, 1.0, 1.0)`（形状 Tanimoto）。
浮動小数点安全のため [0,1] にクランプ。
エラー ID: `invalidInput`（非 Mol またはコンフォーマーなし）、`invalidMethod`、`rdkitError`。

---

## emk.fingerprint

### `emk.fingerprint.morgan(mol, Radius=2, NBits=2048)`
戻り値: `py.rdkit.DataStructs.ExplicitBitVect` — Morgan（ECFP）フィンガープリント。Radius=2 ≈ ECFP4。
エラー ID: `invalidInput`、`rdkitError`。

### `emk.fingerprint.maccs(mol)`
戻り値: `py.rdkit.DataStructs.ExplicitBitVect` — 167 ビット MACCS キーフィンガープリント（ビット 0 は未使用）。
エラー ID: `invalidInput`、`rdkitError`。

### `emk.fingerprint.toArray(fp)`
戻り値: `logical(1,N)` — フィンガープリントを MATLAB logical 配列に変換。N = `fp.GetNumBits()`。
エラー ID: `invalidInput`（Python オブジェクトでない、または `ToBitString()` を持たない）。

---

## emk.similarity

### `emk.similarity.tanimoto(fp1, fp2)`
戻り値: `double` ∈ [0,1] — Tanimoto 係数。同一フィンガープリントなら 1.0 を返す。
ビット長不一致は `rdkitError`。
エラー ID: `invalidInput`（非 Python オブジェクト）、`rdkitError`。

### `emk.similarity.dice(fp1, fp2)`
戻り値: `double` ∈ [0,1] — Dice 係数。2 値ベクトルでは常に Dice ≥ Tanimoto が成立。
エラー ID: `invalidInput`、`rdkitError`。

### `emk.similarity.rankBy(queryFp, dbFps, N=Inf, Metric="tanimoto")`
戻り値: `struct` — `.Indices(1×K)` `.Scores(1×K)` `.Metric(string)`。スコアの降順でソート。
`BulkTanimotoSimilarity`（1 回の IPC ラウンドトリップ）使用。
エラー ID: `invalidQueryFp`、`invalidDbFps`、`invalidN`（0/負/非整数）、`invalidMetric`、`rdkitError`。

### `emk.similarity.matrix(fps, Metric="tanimoto")`
戻り値: `double(N×N)` — 対称類似度行列。対角=1.0。浮動小数点誤差は `(S+S')/2` で除去。
エラー ID: `invalidInput`（非空 cell）、`invalidMetric`、`rdkitError`。

---

## emk.io

### `emk.io.readSdf(filePath)`
戻り値: `1×N cell` — SDF ファイルから読み込んだ Mol cell 配列。
`SDMolSupplier(removeHs=True, sanitize=True)` 使用。失敗分子は `logWarn` でスキップ。
エラー ID: `invalidInput`、`fileNotFound`、`rdkitError`。

### `emk.io.writeSdf(mols, filePath)`
SDF ファイルに分子を書き込む。親ディレクトリが存在する必要がある。既存ファイルは上書き。
エラー時も `writer.close()` は必ず呼ばれる。
エラー ID: `invalidInput`、`invalidMol`（非 Python オブジェクト要素）、`dirNotFound`、`rdkitError`。

### `emk.io.readSmilesList(filePath)`
戻り値: `1×N cell` — 1 行 1 SMILES ファイルを読み込む。`#` コメントと空行をスキップ。
タブ/スペース以降の名前列は無視。失敗行は `logWarn` でスキップ。
エラー ID: `invalidInput`、`fileNotFound`、`allLinesFailed`。

---

## emk.viz

### `emk.viz.draw2d(mol, Title="", Width=300, Height=300)`
戻り値: `matlab.ui.Figure` — MATLAB フィギュアに 2D 構造を描画する。
パイプライン: `Compute2DCoords` → `MolToFile` → `imread` → `imshow`（一時ファイル経由）。
エラー ID: `invalidInput`（非 Mol または Width/Height < 1）、`rdkitError`。

---

## emk.db

### `emk.db.searchPubchem(query, Type="name")`
戻り値: `table`（5 列: `CID` double、`IUPACName` `MolecularFormula` `IsomericSMILES` string、`MolecularWeight` double）。
`Type`: `"name"` / `"smiles"` / `"cid"` / `"inchikey"`。`webread` のみ使用（Python 不要）。
エラー ID: `invalidInput`、`invalidType`、`notFound`（HTTP 404）、`networkError`。

### `emk.db.pubchemFetch(identifier, NameSpace="name", MaxSynonyms=10)`
PubChemPy を使って PubChem から拡張プロパティセット（同義語、InChI/InChIKey、XLogP、TPSA、HBD/HBA 等）を取得し struct を返す。
`searchPubchem` の強化版（Python 必要）。
`NameSpace`: `"name"` / `"smiles"` / `"cid"` / `"inchi"` / `"inchikey"` / `"formula"`。数値スカラー CID も可。
主な struct フィールド: `CID`（double）、`IUPACName` `MolecularFormula` `IsomericSMILES` `InChI` `InChIKey`（string）、
`MolecularWeight` `XLogP` `TPSA` `HBondDonors` `HBondAcceptors` `RotatableBonds` `HeavyAtomCount` `Charge` `Complexity`（double）、
`Synonyms`（string 配列）。
必要条件: `emk.setup.installExtra("pubchempy")` が必要。
エラー ID: `invalidInput`、`invalidNamespace`、`libraryNotFound`、`notFound`、`pythonError`。

### `emk.db.searchChembl(query, Type="name")`
戻り値: `table`（8 列: `ChEMBLID` `Name` `MolecularFormula` `SMILES` `InChIKey` string、
`MolecularWeight` `ALogP` `HBondDonors` `HBondAcceptors` double）。
`Type`: `"name"` / `"chemblid"` / `"smiles"` / `"inchikey"`。`webread` のみ使用。
エラー ID: `invalidInput`、`invalidType`、`notFound`、`networkError`。

### `emk.db.searchChemblTarget(query, TargetType="SINGLE PROTEIN", MaxRows=10)`
戻り値: `table`（4 列: `TargetChEMBLID` `PreferredName` `Organism` `TargetType` string）。
`pref_name__icontains` による部分一致検索。`webread` のみ使用（Python 不要）。
エラー ID: `invalidInput`、`invalidOptions`、`notFound`、`networkError`。

### `emk.db.getChemblActivity(targetId, ActivityType="IC50", MaxRows=50, MinActivity_nM=Inf)`
戻り値: `table`（5 列: `MoleculeChEMBLID` `Name` `SMILES` `ActivityType` string、`Value_nM` double）。
`standard_relation="="` かつ単位が nM の行のみ返す。`Value_nM <= MinActivity_nM` の行のみ保持（デフォルト: 全行）。
`webread` のみ使用。
エラー ID: `invalidInput`、`invalidOptions`、`notFound`、`networkError`。

---

## emk.repro

### `emk.repro.verify(metrics, criteria)`
RF03 受入基準に対して再現メトリクスを検証する。
`metrics`: 数値メトリクスフィールドを持つ struct（例: `metrics.rmse_cv = 1.017`）。
`criteria`: 受入境界の struct; 各フィールドは `"upper"` および/または `"lower"` キーを持つ sub-struct。
戻り値: `struct` — `.pass`（logical）、`.details`（メトリクスごと: `.value .pass .criteria`）、`.report`（char、整形済みサマリー）。
エラー ID: `invalidInput`（非 struct）、`missingMetric`（criteria フィールドが metrics に存在しない）。

---

## emk.util

### `emk.util.isOnline()`
戻り値: `logical` — MATLAB Online 上なら `true`。検出順序: `ismatlabonline()` → `MATLAB_ONLINE` 環境変数 → Linux x64 ヒューリスティック。

### `emk.util.benchmarkBatch(smilesList, descriptorNames=<all10>)`
戻り値: `struct` — `.nMols` `.parseSec` `.batchSec` `.totalSec` `.molsPerSec` `.tbl`。スループット計測用。
エラー ID: `invalidInput`、`emptyInput`。

---

## src/config

### `emkLoadConfig()`
戻り値: `struct` — `.python` `.rdkit` `.runtime` `.output` `.run`。
優先順位: `EMK_<SECTION>_<KEY>` 環境変数 > `config/settings.json` > デフォルト値。

---

## src/util

### `logInfo(msg, ...)` / `logWarn(msg, ...)` / `logError(msg, ...)` / `logDebug(msg, ...)`
`[HH:MM:SS][LEVEL]  メッセージ` 形式で出力。`logDebug` は `EMK_LOG_VERBOSE=1` 時のみ出力。
`fprintf` の直接使用は禁止。

### `logProgress(i, n, label)`
`i/n` の割合でループ進捗バーを表示する。

### `logSection(scriptId, label, layer)`
チュートリアルスクリプトの各 `%%` セクション先頭で呼び出すバナーログ。
出力例: `[11:36:40][INFO]  --- R01 | Section 0: Setup  [Research L4] ---`
引数: `scriptId` = "F01"/"S01"/"R01"/"A01"、`label` = セクションヘッダーテキスト、
`layer` = "Foundation L1" / "Stories L2" / "Analytics L3" / "Research L4"。

### `makeRunDir(Prefix="", BaseDir="result/runs")`
戻り値: `char` — `result/runs/yyyyMMdd_HHmmss[_Prefix]` を作成して返す。
`mkdir` の直接使用は禁止。
