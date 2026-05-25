# Algorithm Guide — emk.fingerprint

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## 5.1 `morgan`

**設計意図**: Morgan（円形）フィンガープリントを生成する。

**アルゴリズム概要**:
- Extended-Connectivity Fingerprint (ECFP) の実装
- 各原子の近傍を半径 $r$ まで探索し、ハッシュ化して固定長ビットベクトルに折りたたむ

$$
\text{ECFP}_{2r} = \text{fold}\left(\bigcup_{a \in \text{atoms}} h(N_r(a))\right)
$$

ここで $N_r(a)$ は原子 $a$ の半径 $r$ の近傍、$h$ はハッシュ関数。

**パラメータ選択の根拠**:
- デフォルト半径 = 2（ECFP4 相当）: 医薬品化学で最も広く使用される設定
- デフォルトビット数 = 2048: 衝突率と計算コストのバランス

**引用文献**:
- Rogers, D. & Hahn, M. (2010). "Extended-Connectivity Fingerprints". *J. Chem. Inf. Model.* 50(5), 742-754. DOI: 10.1021/ci100050t
- Morgan, H.L. (1965). "The Generation of a Unique Machine Description for Chemical Structures". *J. Chem. Doc.* 5(2), 107-113.

**RDKit API**: `py.rdkit.Chem.rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)` → `gen.GetFingerprint(mol)`  
Generator API（rdkit 2022+ 推奨）を使用。非推奨の `GetMorganFingerprintAsBitVect` は使用しない。  
返り値は `py.rdkit.DataStructs.cDataStructs.ExplicitBitVect`（ADR-002: Python 参照保持）。  
Python キーワード引数は `pyargs("radius", int32(r), "fpSize", int32(n))` で渡す。

**テスト戦略** (`tests/unit/TestFingerprint.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_morgan_numericInput_throwsInvalidInput` | double 入力 → `invalidInput` | 不要 |
| `test_morgan_stringInput_throwsInvalidInput` | string 入力 → `invalidInput` | 不要 |
| `test_morgan_charInput_throwsInvalidInput` | char 入力 → `invalidInput` | 不要 |
| `test_morgan_emptyInput_throwsInvalidInput` | 空行列 → `invalidInput` | 不要 |
| `test_morgan_returnsPythonObject` | 戻り値の class が `"py."` で始まる | 要 |
| `test_morgan_returnHasGetNumBitsMethod` | `fp.GetNumBits()` が 2048 を返す | 要 |
| `test_morgan_defaultNBits_is2048` | デフォルト NBits=2048 | 要 |
| `test_morgan_NBits1024_produces1024bits` | `NBits=1024` → 1024 bits | 要 |
| `test_morgan_ethanol_numOnBits_isPositive` | エタノールの ON ビット数 > 0 | 要 |
| `test_morgan_ethanol_idempotency_tanimoto1` | 同一 mol 2 回呼び出し → Tanimoto=1.0 | 要 |
| `test_morgan_ethanol_vs_benzene_tanimotoLessThan1` | エタノール vs ベンゼン → Tanimoto < 1.0 | 要 |
| `test_morgan_radius1_vs_radius2_aspirin_differ` | アスピリンの ECFP2 ≠ ECFP4 | 要 |
| `test_morgan_radius2_vs_radius3_aspirin_differ` | アスピリンの ECFP4 ≠ ECFP6 | 要 |
| `test_morgan_radius0_doesNotThrow` | `Radius=0` でエラーなし | 要 |

**冪等性テストの根拠**: 同一 mol オブジェクトを 2 回渡した場合、RDKit の決定論的ハッシュアルゴリズムにより常に同一ビット列が生成される。Tanimoto = 1.0 が exact equality の代わりになる（ExplicitBitVect の MATLAB 直接比較が困難なため）。

**半径差分テストの根拠**: アスピリン（C9H8O4、芳香環 + 2 つのカルボニル）を使用。Radius=1 では直接結合原子のみを捕捉、Radius=2 では 2 ホップ先も捕捉するため FP が必ず異なる。エタノールのような小さな分子では半径が大きくても差異が出ない場合があるため、構造的に豊富なアスピリンを選択。

---

## 5.2 `maccs`

**設計意図**: MACCS (Molecular ACCess System) 166 公開鍵ベースのフィンガープリント。
Morgan FP がハッシュベースの「円形エンコーディング」であるのに対し、MACCS は定義済み
SMARTS パターンとの照合結果をそのままビットに格納するため、各ビットが化学的意味を持つ
（解釈可能性が高い）。スキャフォールドレベルの類似度比較・薬効団フィーチャー比較に広く使用される。

**アルゴリズム概要**:
1. `isa(mol, "py.rdkit.Chem.rdchem.Mol")` で入力型を検証
2. `py.rdkit.Chem.MACCSkeys.GenMACCSKeys(mol)` を呼び出す
3. Python 例外 → `rdkitError` に変換（`try/catch`）
4. 返り値（167-bit ExplicitBitVect）を Python 参照のまま返す（ADR-002）

**ビット長 167 の根拠**:
RDKit の MACCS 実装では bit 0 が未使用（常に 0）で、bits 1–166 が MACCS structural keys に対応する。
そのため `fp.GetNumBits()` は 167 を返す。`numel(bits) == 167` という点は認知しておく必要がある。

**Morgan との比較**:

| | MACCS | Morgan (ECFP4) |
|---|---|---|
| ビット長 | 167（固定） | 2048（デフォルト、変更可） |
| エンコード方式 | SMARTS パターン照合 | 原子近傍のハッシュ折りたたみ |
| 解釈可能性 | 高い（各ビットが定義済み） | 低い（ハッシュ衝突あり） |
| 汎用性 | やや低い（定義済みパターンに限定） | 高い（任意の部分構造を捕捉） |

**引用文献**:
- MDL Information Systems. MACCS Structural Keys. MDL Information Systems, Inc.
- Durant, J.L. et al. (2002). *J. Chem. Inf. Comput. Sci.* 42(6), 1273–1280. DOI: 10.1021/ci010132r
- RDKit Documentation: [MACCSkeys](https://www.rdkit.org/docs/GettingStartedInPython.html#maccs-keys)

**テスト戦略** (`tests/unit/TestFingerprint.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_maccs_numericInput_throwsInvalidInput` | double 入力 → `invalidInput` | 不要 |
| `test_maccs_stringInput_throwsInvalidInput` | string 入力 → `invalidInput` | 不要 |
| `test_maccs_returnIsPythonObject` | 戻り値の class が `"py."` で始まる | 要 |
| `test_maccs_ethanol_numBits_is167` | エタノールの GetNumBits() = 167 | 要 |
| `test_maccs_ethanol_numOnBits_isPositive` | ON ビット数 > 0 | 要 |
| `test_maccs_benzene_numOnBits_greaterThan_ethanol` | ベンゼン > エタノール（芳香環キーで差異） | 要 |
| `test_maccs_ethanol_idempotency_arrayEqual` | 同一 mol 2 回 → 同一配列（冪等性） | 要 |
| `test_maccs_ethanol_toArray_sumEqualGetNumOnBits` | `sum(toArray(fp)) == fp.GetNumOnBits()` エタノール ★ | 要 |
| `test_maccs_aspirin_toArray_sumEqualGetNumOnBits` | 同上、アスピリン | 要 |
| `test_maccs_toArray_length167` | `numel(toArray(maccs(mol))) == 167`（E2E） | 要 |
| `test_maccs_toArray_returnsLogical` | `class(toArray(maccs(mol))) == "logical"` — 型検証 | 要 |
| `test_maccs_toArray_returnsRowVector` | `size(bits,1) == 1` — 行ベクトル検証 | 要 |

**★ クロス検証の重要性**: `sum(toArray(fp)) == fp.GetNumOnBits()` は MACCS の MATLAB 配列変換とPython オブジェクトが一致することを確認する最強テスト。toArray の変換バグはここで検出される。

---

## 5.3 `toArray`

**設計意図**: Python `ExplicitBitVect` を MATLAB `logical(1, N)` 配列に変換し、MATLAB ネイティブの行列演算・`isequal` 比較・型アサーションを可能にする。

**アルゴリズム概要**:
1. 入力が Python オブジェクトでなければ即座に `invalidInput` を投げる（RDKit 不要）
2. `fp.ToBitString()` を 1 回 IPC 呼び出しでビット文字列（例: `"01001..."`) を取得
3. `char(string(...))` で MATLAB `char` 配列に変換後、`== '1'` で `logical` に変換
4. `ToBitString()` を持たない Python オブジェクト（例: `Mol`）は `try/catch` で `invalidInput` に変換

**ADR-002 適合**: Python ↔ MATLAB IPC は `ToBitString()` の 1 回のみ。ビット数分の個別アクセスは行わない。

**テスト戦略** (`tests/unit/TestFingerprint.m`):

| テスト名 | 内容 | RDKit 要否 |
|---|---|---|
| `test_toArray_numericInput_throwsInvalidInput` | `double` 入力 → `invalidInput` | 不要 |
| `test_toArray_stringInput_throwsInvalidInput` | `string` 入力 → `invalidInput` | 不要 |
| `test_toArray_molInput_throwsInvalidInput` | `py.Mol` 入力（`ToBitString` なし）→ `invalidInput` | 必要 |
| `test_toArray_returnsLogical` | 戻り値クラスが `logical` | 必要 |
| `test_toArray_returnsRowVector` | 戻り値が行ベクトル `(1 x N)` | 必要 |
| `test_toArray_defaultFP_length2048` | デフォルト FP の長さ = 2048 | 必要 |
| `test_toArray_length_matchesGetNumBits` | `numel == fp.GetNumBits()` | 必要 |
| `test_toArray_ethanol_sumEqualGetNumOnBits` | `sum(bits) == fp.GetNumOnBits()` エタノール ★ | 必要 |
| `test_toArray_aspirin_sumEqualGetNumOnBits` | 同上、アスピリン | 必要 |
| `test_toArray_idempotency_sameResultTwice` | 同一 FP に 2 回呼び出し → 同一配列 | 必要 |
| `test_toArray_allBitsAreLogical` | 全要素 `∈ {0, 1}`、`islogical` | 必要 |
| `test_toArray_noNaNValues` | NaN 要素なし | 必要 |

**★ TC14 クロス検証の重要性**: `sum(toArray(fp)) == fp.GetNumOnBits()` は MATLAB 変換と Python オブジェクトの ON ビット数が一致することを直接確認する最強のテスト。変換バグ（ビット反転・ずれ等）は必ずここで検出される。

**バリデーション順序の根拠**: Python オブジェクトかどうかのチェックを先に行い RDKit 不要テストを独立させる。これは `emk.descriptor.calculate()` と同一の設計方針。
