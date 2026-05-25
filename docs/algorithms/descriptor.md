# Algorithm Guide — emk.descriptor

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## 4.1 `molWeight`

**設計意図**: 平均分子量（Average Molecular Weight）を返す最初の記述子関数。
M1 MVP で「SMILES → mol → MW」のパイプラインが端から端まで機能することを実証する。

**アルゴリズム概要**:
1. `isa(mol, "py.rdkit.Chem.rdchem.Mol")` で入力型を検証。Mol 以外は即エラー
2. `py.rdkit.Chem.Descriptors.MolWt(mol)` を呼び出す
3. Python 例外 → `rdkitError` に変換（`try/catch`）
4. 返り値（Python `float`）を MATLAB `double` に変換して返す

**数学的定義**:

$$
MW = \sum_{i=1}^{N_{\text{heavy}}} m_i^{\text{avg}} + \sum_{j=1}^{N_H} m_H^{\text{avg}}
$$

ここで $m_i^{\text{avg}}$ は重原子 $i$ の IUPAC 2021 平均原子量、$N_H$ は明示的・暗黙的水素の総数。  
RDKit は暗黙の水素（Implicit H）も含めて合計するため、SMILES で H を省略しても正しい値を返す。

**平均原子量 vs. モノアイソトピック質量の選択根拠**:  
`Descriptors.MolWt` が返す**平均分子量**は、同位体組成を考慮した実験室スケールの質量であり、
薬化学・QSAR での一般的な記述子として広く利用される。
モノアイソトピック質量（`Descriptors.ExactMolWt`）は高分解能質量分析の同位体ピーク帰属に用いられるが、
汎用記述子として平均分子量を優先する。

**参考値（IUPAC 2021 原子量）**:

| 分子 | 分子式 | 計算値 | PubChem CID |
|---|---|---|---|
| エタノール | C₂H₆O | 46.069 g/mol | 702 |
| ベンゼン | C₆H₆ | 78.114 g/mol | 241 |
| アスピリン | C₉H₈O₄ | 180.159 g/mol | 2244 |
| 水 | H₂O | 18.015 g/mol | 962 |

**引用文献**:
- RDKit Documentation: [Descriptors.MolWt](https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-descriptors)
- IUPAC 2021 Standard Atomic Weights: [IUPAC CIAAW](https://ciaaw.org/atomic-weights.htm)

**テスト戦略** (`tests/unit/TestDescriptor.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_molWeight_numericInput_throwsInvalidInput` | 数値入力 → `invalidInput` エラー | 不要 |
| `test_molWeight_numericInput_errorMessage_containsClass` | エラーメッセージに `"double"` を含む | 不要 |
| `test_molWeight_stringInput_throwsInvalidInput` | 生 SMILES 文字列 → `invalidInput` エラー | 不要 |
| `test_molWeight_charInput_throwsInvalidInput` | `char` 入力 → `invalidInput` エラー | 不要 |
| `test_molWeight_emptyInput_throwsInvalidInput` | `[]` → `invalidInput` エラー | 不要 |
| `test_molWeight_ethanol_value` | エタノール MW = 46.069 ± 0.01 g/mol | 要 |
| `test_molWeight_benzene_value` | ベンゼン MW = 78.114 ± 0.01 g/mol | 要 |
| `test_molWeight_aspirin_value` | アスピリン MW = 180.159 ± 0.01 g/mol | 要 |
| `test_molWeight_water_value` | 水 MW = 18.015 ± 0.01 g/mol（暗黙 H の処理確認） | 要 |
| `test_molWeight_outputIsDoubleScalar` | 戻り値が `double` スカラー | 要 |
| `test_molWeight_outputIsRealValued` | `isreal(mw) == true` | 要 |
| `test_molWeight_ethanol_isPositive` | MW > 0 | 要 |
| `test_molWeight_ethanol_isFinite` | `isfinite(mw)` | 要 |
| `test_molWeight_aspirin_isFinite` | 複雑分子でも `isfinite(mw)` | 要 |

**AbsTol = 0.01 g/mol の根拠**:  
PubChem の公称値は小数第 2 位まで（例: 46.07 g/mol）。RDKit バージョン間の差異（~0.001 g/mol 程度）を吸収しつつ、
重大な計算誤り（例: 水素を無視した MW）は確実に検出できる。

---

## 4.2 `calculate`

**設計意図**: 10 種の代表的な physicochemical 記述子を一括計算し、MATLAB struct で返す。

**アルゴリズム概要**:
1. `isa(mol, "py.rdkit.Chem.rdchem.Mol")` で入力型を検証
2. `descriptorNames` が省略された場合は全 10 記述子を対象とする
3. 不明な記述子名 → `unknownDescriptor` エラー（RDKit 呼び出し前に弾く）
4. `for` ループで各記述子を `computeOne_()` で計算
5. 各記述子は RDKit API → Python 型 → `double()` 変換のパターンで実装
6. Python 例外は `rdkitError` に変換して返す

**記述子の定義・数学的根拠**:

| 記述子 | RDKit API | 定義・引用文献 |
|---|---|---|
| `MolWt` | `Descriptors.MolWt` | IUPAC 2021 平均原子量の総和（暗黙 H 含む） |
| `ExactMolWt` | `Descriptors.ExactMolWt` | モノアイソトピック質量（最も豊富な同位体の総和） |
| `LogP` | `Descriptors.MolLogP` | Wildman-Crippen 分配係数。Wildman & Crippen (1999) *J. Chem. Inf. Comput. Sci.* 39(5):868–873 |
| `TPSA` | `Descriptors.TPSA` | Ertl (2000) 寄与法による位相的極性表面積。*J. Med. Chem.* 43(20):3714–3717 |
| `NumHAcceptors` | `rdMolDescriptors.CalcNumHBA` | Ertl 定義 HBA (O, N で特定 SMARTS に合致する原子) |
| `NumHDonors` | `rdMolDescriptors.CalcNumHBD` | Ertl 定義 HBD (N-H, O-H を持つ原子) |
| `NumRotatableBonds` | `rdMolDescriptors.CalcNumRotatableBonds` | 末端重原子への単結合・環内結合を除いた単結合数 (strict definition) |
| `RingCount` | `rdMolDescriptors.CalcNumRings` | SSSR（最小スパニングリングセット）に基づく環数 |
| `FractionCSP3` | `Descriptors.FractionCSP3` | sp3 炭素数 / 全炭素数。Lovering et al. (2009) *J. Med. Chem.* 52(21):6752–6756 |
| `HeavyAtomCount` | `mol.GetNumHeavyAtoms()` | 水素を除く原子数 |
| `MolFormula` | `rdMolDescriptors.CalcMolFormula` | Hill 表記の分子式文字列 (例: "C2H6O")。戻り値は **string**（他は double）。デフォルトセット外。Hill (1900) *J. Am. Chem. Soc.* 22(8):478–494 |

**引用文献**:
- Wildman, S.A. & Crippen, G.M. (1999). *JCICS* 39(5):868–873
- Ertl, P. et al. (2000). *J. Med. Chem.* 43(20):3714–3717
- Lovering, F. et al. (2009). *J. Med. Chem.* 52(21):6752–6756
- RDKit Documentation: [List of Available Descriptors](https://www.rdkit.org/docs/GettingStartedInPython.html#list-of-available-descriptors)

**エタノール参考値（RDKit 2024.3.x / IUPAC 2021）**:

| 記述子 | 値 | 根拠 |
|---|---|---|
| MolWt | 46.069 g/mol ± 0.01 | IUPAC 原子量の総和、PubChem CID 702 |
| TPSA | 20.23 Å² ± 0.1 | Ertl 2000、PubChem CID 702 (20.2 Å²) |
| NumHAcceptors | 1 | O 原子 1 個（Ertl 定義） |
| NumHDonors | 1 | OH 基 1 個 |
| NumRotatableBonds | 0 | strict: 末端原子への結合を除外 |
| RingCount | 0 | 環なし |
| FractionCSP3 | 1.0 | 2 C / 2 C = 1.0（両炭素とも sp3） |
| HeavyAtomCount | 3 | 2 C + 1 O |

**テスト戦略** (`tests/unit/TestDescriptor.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_calculate_numericInput_throwsInvalidInput` | 数値入力 → `invalidInput` | 不要 |
| `test_calculate_unknownName_throwsUnknownDescriptor` | 不明な記述子名 → `unknownDescriptor` | 不要 |
| `test_calculate_unknownName_errorMessage_containsName` | エラーメッセージに不明名を含む | 不要 |
| `test_calculate_allDescriptors_returnsStruct` | 戻り値が struct | 要 |
| `test_calculate_allDescriptors_hasAllTenFields` | 全 10 フィールドを持つ | 要 |
| `test_calculate_ethanol_MolWt` | エタノール MolWt = 46.069 ± 0.01 g/mol | 要 |
| `test_calculate_ethanol_TPSA` | エタノール TPSA = 20.23 ± 0.1 Å² | 要 |
| `test_calculate_ethanol_NumHAcceptors` | NumHAcceptors = 1 | 要 |
| `test_calculate_ethanol_NumHDonors` | NumHDonors = 1 | 要 |
| `test_calculate_ethanol_NumRotatableBonds` | NumRotatableBonds = 0 (strict) | 要 |
| `test_calculate_ethanol_RingCount` | RingCount = 0 | 要 |
| `test_calculate_ethanol_FractionCSP3` | FractionCSP3 = 1.0 ± 0.001 | 要 |
| `test_calculate_ethanol_HeavyAtomCount` | HeavyAtomCount = 3 | 要 |
| `test_calculate_subset_returnsOnlyRequestedFields` | ["MolWt","LogP"] → 2 フィールドのみ | 要 |
| `test_calculate_allDescriptors_allValuesAreDoubleScalar` | 全値が double スカラー | 要 |
| `test_calculate_allDescriptors_allValuesAreFinite` | 全値が有限（NaN/Inf なし） | 要 |
| `test_calculate_aspirin_HeavyAtomCount` | アスピリン HeavyAtomCount = 13 | 要 |
| `test_calculate_aspirin_NumHDonors` | アスピリン NumHDonors = 1 | 要 |
| `test_calculate_aspirin_RingCount` | アスピリン RingCount = 1 | 要 |
| `test_calculate_MolFormula_notUnknown_passesValidation` | "MolFormula" は unknownDescriptor を発生させない | 不要 |
| `test_calculate_MolFormula_notInDefaultSet` | デフォルトセットに MolFormula を含まない | 要 |
| `test_calculate_MolFormula_outputIsStringScalar` | 戻り値が string スカラー | 要 |
| `test_calculate_MolFormula_outputIsNonEmpty` | 空文字でない | 要 |
| `test_calculate_MolFormula_ethanol` | エタノール = "C2H6O" | 要 |
| `test_calculate_MolFormula_aspirin` | アスピリン = "C9H8O4" | 要 |
| `test_calculate_MolFormula_benzene` | ベンゼン = "C6H6" | 要 |
| `test_calculate_MolFormula_withNumericDescriptor` | ["ExactMolWt","MolFormula"] → 2 フィールド（型混在可） | 要 |

---

## 4.3 `batchCalculate`

**設計意図**: 複数分子の記述子を一括計算し MATLAB table で返す。
ADR-002 rev.3 の IPC 最小化原則に従い、Python 側ヘルパー `batch_descriptors.py` を
使って分子数 N の全計算を 1 回の IPC 往復で完了する。ヘルパーが利用できない場合はループ fallback を使用する。

**アルゴリズム概要**:
1. `mols` が cell array でない場合は `invalidInput` エラー
2. `descriptorNames` のバリデーション（`calculate` と同一の許容名リストを使用）
3. `NaN` で初期化した `nMols × nDesc` 行列を確保
4. Python バッチパス（M3-2 最適化）:
   - MATLAB 側で `isa(mol, "py.rdkit.Chem.rdchem.Mol")` を全要素に適用（Python 呼び出し不要）
   - 有効 mol を `py.list` に集約し `batch_descriptors.batch_calculate(pyMols, pyNames)` を 1 回呼び出し
   - 返値 `list[list[float]]` を MATLAB double 行列に変換
   - ヘルパーロード失敗時はループ fallback へ移行
5. ループ fallback: 各 `mols{i}` に `emk.descriptor.calculate` を適用（従来挙動）
6. **M3-1 エラー**: 全分子が失敗した場合 (`nSkip == nMols && nMols > 0`) → `allMolsFailed` を throw
7. `array2table` で MATLAB table に変換

**IPC 最小化の数値的根拠** (ADR-002 rev.3):
- OutOfProcess モードの IPC 往復コスト: ~10–50 ms/call（マシン依存）
- N=100 分子、10 記述子の場合:
  - 旧実装（ループ）: 100 × 10 記述子計算 × IPC コスト ≈ 1–5 s
  - 新実装（バッチ）: 1 IPC 往復 + Python 内ループ ≈ 数十 ms〜数百 ms
- `emk.util.benchmarkBatch()` で実環境計測が可能

**テスト戦略** (`tests/unit/TestDescriptor.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| TB1: 非 cell 入力 | `invalidInput` エラーが RDKit 不要で発火 | なし |
| TB2: 未知記述子名 | `unknownDescriptor` エラー。エラーメッセージに "Supported" を含む | なし |
| TB9: 空 cell | 0 行 × N 列の table。列名が一致 | なし |
| **TB10: 全 mol 無効 → allMolsFailed** | 全非 Mol 要素で `allMolsFailed` | なし |
| TB3: 3 分子 → 3 行 | `height(tbl) == 3` | RDKit |
| TB4: 戻り値型 | `verifyClass(tbl, "table")` | RDKit |
| TB5: 列名（デフォルト） | 10 列すべての列名が一致 | RDKit |
| TB7: 値の一致 | `batchCalculate MolWt == calculate MolWt`（cross-validation） | RDKit |
| TB8: 無効 mol → NaN 行 | 無効 mol の行が NaN。隣接有効行は有限 | RDKit |

---

## 4.4 `benchmarkBatch` (`emk.util`)

**設計意図**: `batchCalculate` の実環境スループットを測定し、最適化効果を定量化する（M3-2 検証用）。

**アルゴリズム概要**:
1. SMILES リストをパース（`emk.mol.fromSmiles` ループ）し、parse 時間を計測
2. `emk.descriptor.batchCalculate(mols, descriptorNames)` を呼び出し、batch 時間を計測
3. スループット統計（mol/s, sec/mol）を計算して struct で返す

**計測の根拠**:
- `tic`/`toc` は MATLAB の高分解能タイマー（通常 1 µs 以下の精度）
- 初回実行は Python 起動コストを含むため、ウォームアップ呼び出し後に計測することを推奨
- `parseSec`（SMILES→Mol）と `batchSec`（記述子計算）を分離することで、ボトルネックを特定できる

---

## 4.5 `mordred`

**設計意図**: Mordred（mordredcommunity フォーク）の 2D 記述子計算を MATLAB 関数として提供する。
RDKit が提供する 10 種の記述子（`calculate`）では不足する QSAR 前処理（~1800 記述子）に対応する。

**前提条件**: `emk.setup.installExtra("mordred")` で `mordredcommunity` を Embedded Python に pip install 済みであること。

**アルゴリズム概要**:
1. `isa(mol, "py.rdkit.Chem.rdchem.Mol")` で入力型を検証（Python 不要）
2. `py.importlib.import_module("mordred")` で Mordred の可用性を確認 → 不可なら `libraryNotFound`
3. `run_mordred.py` ヘルパーを `src/+emk/+util/python/` からロード（`sys.path` 挿入、1 回のみ）
4. `mordred_calculate_list(mol, names)` を 1 IPC 往復で呼び出し → `(names_list, values_list)` を取得
5. MATLAB 側で struct に変換（`colNames(k)` → `values(k)` のフィールド代入）

**IPC 設計（ADR-002 rev.3）**:
- `mordred_calculate_list` は Python tuple `(names_list: list[str], values_list: list[float])` を返す
- MATLAB 側では `cell(pyTuple)` → `string(cell(names))` + `double(py.array.array("d", vals))` の 1 回変換のみ
- 全 ~1800 記述子の場合でも IPC 往復は 1 回（従来の dict lookback N × M 回 → 1 回に削減）

**Mordred 記述子体系**:
- Moriwaki et al. (2018) が提案した Python ライブラリ。2D/3D 記述子を含む
- `ignore_3D=True` で 2D のみ計算（RDKit の 3D 座標生成不要）
- 記述子数: 2D のみで ~1800 種（バージョンにより変動）
- 代表的な記述子: `MW`（分子量）, `ALogP`（Crippen LogP）, `nRot`（回転結合数）, `nHBd`（HBD）, `nHBa`（HBA）, `TPSA`（極性表面積）

**Calculator キャッシュ**:
`run_mordred.py` はモジュールレベルの `_CALC_ALL` 変数で Calculator を 1 回だけ初期化する。
MATLAB セッション中は同一 Python プロセスを使い回すため、2 回目以降の呼び出しはキャッシュから返す。

**引用文献**:
- Moriwaki H. et al. (2018). Mordred: a molecular descriptor calculator. *Journal of Cheminformatics* 10:4. https://doi.org/10.1186/s13321-018-0258-y
- mordredcommunity: https://github.com/mordred-descriptor/mordredcommunity

**テスト戦略** (`tests/unit/TestMordred.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| TC1 非 Mol → `invalidInput` | 数値入力 → 即時エラー（Python 不要） | なし |
| TC1b エラーメッセージにクラス名 | `"double"` を含む | なし |
| TC2 string/char/logical/empty 全拒否 | 非 Mol 全型が `invalidInput` | なし |
| TC4 valid Mol → struct | struct 型を返す | Mordred + RDKit |
| TC5 全フィールドが double scalar | `verifyClass` + `verifySize` + `verifyFalse(isinf)` | Mordred + RDKit |
| TC6 subset request 返却フィールド一致 | 要求数と一致するフィールド数 | Mordred + RDKit |
| TC7 MW ≈ 46.07（AbsTol=0.5） | エタノール分子量の参照値 ★★★ | Mordred + RDKit |
| TC9 nRot = 0 for ethanol | 回転結合なし ★★★ | Mordred + RDKit |
| TC10 nHBd = 1 for ethanol | OH 基 1 個 = HBD 1 ★★★ | Mordred + RDKit |
| TC10b nHBd == RDKit HBondDonors | emk.descriptor.calculate とのクロスバリデーション ★★★ | Mordred + RDKit |

**★★★ 参照データの根拠**:
- `MW = 46.07 g/mol`: Mordred は RDKit と同じ平均原子量ベースだが実装差で ±0.5 g/mol 程度ずれる場合があるため AbsTol=0.5 を設定
- `nRot = 0`: 厳密定義（strict）では末端原子への結合を除外するため、エタノールの C-O 結合も回転結合カウントされない
- `nHBd = 1`: エタノールの OH 基のみ。CH₂ 結合の H は HBD に含まれない（Ertl 定義）
- TC10b クロスバリデーションの意義: Mordred ラッパーが RDKit の同等記述子と整合することを保証し、ラッパー実装のドリフトを検出する

---

## 4.6 `mordredBatch`

**設計意図**: `mordred` の複数分子版。`run_mordred.py` の `mordred_batch_matrix` を介して
N 分子の ~1800 記述子を 1 IPC 往復で計算し MATLAB table に変換する。

**アルゴリズム概要**:
1. `iscell(mols)` でバリデーション（Python 不要）
2. MATLAB 側で `isa(mols{i}, "py.rdkit.Chem.rdchem.Mol")` を全要素に適用してタイプチェック
3. 無効要素は `py.None` に置換した `py.list` を Python に渡す（NaN 行として返る）
4. `mordred_batch_matrix(pyMols, pyNames)` を 1 IPC 往復で呼び出し → `(names_list, matrix)` を取得
5. 各行を `double(py.array.array("d", pyRowCell{i}))` で変換
6. 全行が NaN なら `allMolsFailed` をスロー
7. `array2table` で table に変換
8. **空 cell の扱い**: `numel(mols) == 0` の場合、`allMolsFailed` はスローせず 0 行テーブルを返す（`nMols > 0 && nSkip == nMols` 条件が偽になる）

**IPC 往復数（ADR-002 rev.3）**:
- `mordred_batch_matrix` は Python tuple `(names: list[str], matrix: list[list[float]])` を返す
- MATLAB からの IPC 往復 = **O(N)**（行数分の `py.array.array` 変換）
- `mordred_calculate` を N 回ループした場合: **O(N × M)** IPC 往復（M ≈ 1800）

**テスト戦略** (`tests/unit/TestMordred.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| TB1 非 cell → `invalidInput` | string / numeric / struct → 即時エラー（Python 不要） | なし |
| TB2 valid batch → table | 2 分子 → table 型 | Mordred + RDKit |
| TB3 テーブル次元 | N rows × M cols | Mordred + RDKit |
| TB4 列名一致 | 要求名と VariableNames が一致 | Mordred + RDKit |
| TB5 無効 Mol → NaN 行 | 無効要素行が全 NaN、エラーなし | Mordred + RDKit |
| TB5b 有効行が NaN でない | 同バッチの valid 行に有限値あり ★★★ | Mordred + RDKit |
| TB6 全無効 → `allMolsFailed` | 全非 Mol → 確実にスロー | Mordred + RDKit |
| TB7 batch vs single クロスバリデーション | batch MW == single MW（両者とも NaN でないことも検証） ★★★ | Mordred + RDKit |
| TB8 空 cell + 明示名 → 0 行テーブル | `{}` → 0×M table、列名一致 ★★★ | Mordred + RDKit |

**★★★ 補強テストの根拠**:
- TB5b: TB5（NaN 行確認）だけでは「全行が NaN になるバグ」を見逃す。valid 隣接行の有限値確認が必須
- TB7 NaN ガード: `verifyEqual(NaN, NaN)` は MATLAB では失敗するが `abs(NaN - NaN)` が `NaN` のため AbsTol チェックも失敗する。NaN 状態を先に `verifyFalse(isnan(...))` で確認してから値比較することで、Python 計算失敗を明確に検出できる
- TB8: 0 分子バッチは実パイプラインで filter 後の残余 0 件ケースとして発生し得る。エラーではなく空テーブルを返すことが契約であり、テストで保証する

---

## 4.7 `mordredNames`

**設計意図**: `mordred` / `mordredBatch` の `descriptorNames` 引数に渡す名前リストをユーザーが
インタラクティブに探索できるようにする。記述子のフィルタリングや subset 計算の起点として利用する。

**アルゴリズム概要**:
1. `py.importlib.import_module("mordred")` で可用性を確認
2. `run_mordred.py` をロードし `mordred_list_names()` を呼び出す
3. Python 側の `_get_calc_all()` でキャッシュ済み Calculator から `sorted(str(d) for d in ...)` を返す
4. MATLAB 側で `string(cell(...))` に変換して返す（Python list → 1×N cell → 1×N string）

**テスト戦略** (`tests/unit/TestMordred.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| TN1 string 型 | `verifyClass(names, "string")` | Mordred |
| TN1b 行ベクトル | `size(names, 1) == 1`（1×N 保証）★★★ | Mordred |
| TN2 件数 >= 1000 | Mordred 2D ≈ 1800 → 安全下限は 1000 | Mordred |
| TN3 "MW" を含む | 基本的な MW 記述子が存在する | Mordred |

**★★★ TN1b 行ベクトル保証の根拠**:
`function_reference.md` の API 仕様は `string(1×N)` を明示している。
`mordred()` / `mordredBatch()` の `descriptorNames` 引数に行ベクトルを渡すことを前提としているため、
列ベクトルが返ると下流ユーザーの `for k = 1:numel(names)` ループが正しく動作するが、
`reshape(descriptorNames, 1, [])` を経由しているため問題が潜在する可能性があり、型保証として検証する。

**引用文献**: → 4.5 `mordred` と共通
