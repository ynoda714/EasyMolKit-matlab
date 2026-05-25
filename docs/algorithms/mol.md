# Algorithm Guide — emk.mol

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## 3.1 `fromSmiles`

**設計意図**: SMILES 文字列から RDKit Mol オブジェクトを生成し、MATLAB から化学解析パイプラインへの入口となる関数。

**アルゴリズム概要**:
1. `ischar(smiles) || isStringScalar(smiles)` を確認（RDKit 不要）
2. `strlength(strtrim(smiles)) == 0` で空白 SMILES チェック（RDKit 不要）
3. `py.rdkit.Chem.MolFromSmiles(char(smiles))` で RDKit に渡す
4. `isa(mol, "py.NoneType")` で失敗を検出 → `invalidSmiles` エラー
5. 有効な `py.rdkit.Chem.rdchem.Mol` オブジェクトを返す

**`py.NoneType` 検出の根拠**:
`MolFromSmiles` は SMILES が無効な場合に `None` を返す（例外を投げない RDKit の仕様）。
`isa(mol, "py.NoneType")` はこれを検出する。`isempty(mol)` では正常な Mol も None も
どちらも空として判定される場合があるため `isa` を使用する。

**引用文献**:
- Weininger, D. (1988). SMILES, a chemical language and information system. *JCICS* 28(1):31–36.
- RDKit: Mol object (py.rdkit.Chem.rdchem.Mol) — Python reference, retained per ADR-002

**テスト戦略** (`tests/unit/TestMol.m`):

| テスト ID | 検証内容 | RDKit 要否 |
|---|---|---|
| TC1: `test_fromSmiles_charInput_accepted` | char 型入力を受け入れる ★ | 要 |
| TC2: `test_fromSmiles_stringInput_accepted` | string 型入力を受け入れる | 要 |
| TC3: `test_fromSmiles_numericInput_throwsInvalidInput` | double → `invalidInput` | 不要 |
| TC4: `test_fromSmiles_emptyString_throwsInvalidSmiles` | `""` → `invalidSmiles` | 不要 |
| TC5: `test_fromSmiles_whitespaceOnly_throwsInvalidSmiles` | `"   "` → `invalidSmiles` | 不要 |
| TC6: `test_fromSmiles_validSmiles_returnsPyMol` | 有効 SMILES → `py.rdkit.*` | 要 |
| TC7: `test_fromSmiles_invalidSmiles_throwsInvalidSmiles` | `"XYZ123"` → `invalidSmiles` | 要 |
| TC8: `test_fromSmiles_ethanol_atomCountIsThree` | エタノール原子数 = 3（重原子）★★ | 要 |
| TC9: `test_fromSmiles_benzene_atomCountIsSix` | ベンゼン原子数 = 6 ★★ | 要 |
| TC10: `test_fromSmiles_aspirin_atomCountIsThirteen` | アスピリン原子数 = 13 ★★ | 要 |
| TC11: `test_fromSmiles_emptyChar_throwsInvalidSmiles` | `''` (空 char) → `invalidSmiles` | 不要 |
| TC12: `test_fromSmiles_cellInput_throwsInvalidInput` | cell 配列 → `invalidInput` | 不要 |
| TC13: `test_fromSmiles_pyObjInput_throwsInvalidInput` | Mol オブジェクト → `invalidInput` | 不要 |
| TC14: `test_fromSmiles_validSmiles_classCheck` | 戻り値クラスが `py.rdkit.Chem.rdchem.Mol` | 要 |
| TC15: `test_fromSmiles_returnIsNotPyNone` | 戻り値が `py.NoneType` でない | 要 |
| TC16: `test_fromSmiles_extraWhitespace_parsedSuccessfully` | `"  CCO  "` → 正常パース（trim 動作確認）★ | 要 |
| TC17: `test_fromSmiles_aspirin_getNumAtoms_crossValidation` | MATLAB 計算 == py.GetNumAtoms() ★★★ | 要 |

**★★ 参照原子数（重原子）**:
- エタノール (`CCO`): 3
- ベンゼン (`c1ccccc1`): 6
- アスピリン (`CC(=O)Oc1ccccc1C(=O)O`): 13

**★★★ TC17 クロス検証の重要性**:
MATLAB 側の `mol.GetNumAtoms()` 結果と Python API の結果を一致確認。
型変換バグ（例: `int32` が `double` に化ける）を検出する。

---

## 3.2 `toSmiles`

**設計意図**: RDKit Mol オブジェクトを Canonical SMILES 文字列（MATLAB string 型）に変換する。

**アルゴリズム概要**:
1. `isa(mol, "py.rdkit.Chem.rdchem.Mol")` で型検証（RDKit 不要）
2. `py.rdkit.Chem.MolToSmiles(mol)` で Canonical SMILES を生成
3. `string(...)` で MATLAB string に変換して返す

**Canonical SMILES の定義**:
RDKit の Morgan アルゴリズムに基づいた一意な SMILES 表現。
同一分子は入力 SMILES の書き方によらず同一の Canonical SMILES を返す。

**テスト戦略**:

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_toSmiles_numericInput_throwsInvalidInput` | double → `invalidInput` | 不要 |
| `test_toSmiles_ethanol_returnsString` | 戻り値が string 型 | 要 |
| `test_toSmiles_ethanol_roundTrip` | `fromSmiles(toSmiles(mol))` がエラーなし | 要 |
| `test_toSmiles_idempotency` | `toSmiles(fromSmiles(toSmiles(m))) == toSmiles(m)` ★ | 要 |

**★ 冪等性テストの重要性**: SMILES → Mol → SMILES の変換が安定していることを確認する。

---

## 3.3 `isValid`

**設計意図**: SMILES が有効な分子を表すかを論理値で返す。例外を投げない。

**アルゴリズム概要**:
1. `ischar(s) || isStringScalar(s)` チェック（RDKit 不要）。false の場合は `invalidInput`
2. `try ... emk.mol.fromSmiles(s); result = true; catch me ... ` でエラーを catch
3. `me.identifier` が `"emk:mol:fromSmiles:rdkitError"` の場合 → 再スロー（RDKit エラーは透過）
4. `"emk:mol:fromSmiles:invalidSmiles"` の場合 → `false` を返す

**RDKit エラーを再スローする根拠**:
`rdkitError` は「SMILES が不正」ではなく「RDKit が動作していない」を意味する。
`isValid` は SMILES の妥当性を判定する関数であり、RDKit の可用性を隠蔽すべきでない。

**テスト戦略**:

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_isValid_numericInput_throwsInvalidInput` | double → `invalidInput` | 不要 |
| `test_isValid_validSmiles_returnsTrue` | "CCO" → true | 要 |
| `test_isValid_invalidSmiles_returnsFalse` | "XYZ123" → false | 要 |
| `test_isValid_emptyString_returnsFalse` | "" → false（空文字は fromSmiles でエラー → false） | 不要 |
| `test_isValid_returnIsLogicalScalar` | 戻り値が logical スカラー | 要 |

---

## 3.4 `hasSubstruct`

**設計意図**: 分子がサブ構造（SMARTS パターン）を含むかを判定する。スクリーニング・フィルタリングに使用。

**アルゴリズム概要**:
1. `isa(mol, "py.rdkit.Chem.rdchem.Mol")` で mol の型を確認（RDKit 不要）
2. smarts が string/char かつ非空を確認（RDKit 不要）
3. `py.rdkit.Chem.MolFromSmarts(char(smarts))` でパターンを生成
4. `isa(pattern, "py.NoneType")` で無効な SMARTS パターンを検出 → `invalidSmarts` エラー
5. `mol.HasSubstructMatch(pattern)` で一致を判定
6. MATLAB `logical` スカラーを返す

**SMARTS パターンの `py.NoneType` 検出**:
`MolFromSmarts` も `MolFromSmiles` と同様に無効入力で `None` を返す。
`py.NoneType` チェックは SMARTS のバリデーションに必要。

**テスト戦略**:

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_hasSubstruct_numericMol_throwsInvalidInput` | double mol → `invalidInput` | 不要 |
| `test_hasSubstruct_numericSmarts_throwsInvalidInput` | double smarts → `invalidInput` | 不要 |
| `test_hasSubstruct_emptySmarts_throwsInvalidSmarts` | "" smarts → `invalidSmarts` | 不要 |
| `test_hasSubstruct_invalidSmarts_throwsInvalidSmarts` | 不正 SMARTS → `invalidSmarts` | 要 |
| `test_hasSubstruct_aspirin_hasCarbonyl` | アスピリン に `"[C](=O)"` → true | 要 |
| `test_hasSubstruct_ethanol_hasOH` | エタノールに `"[OH]"` → true | 要 |
| `test_hasSubstruct_benzene_hasAromaticC` | ベンゼンに `"c"` → true | 要 |
| `test_hasSubstruct_ethanol_noAromaticC` | エタノールに `"c"` → false | 要 |
| `test_hasSubstruct_returnsLogicalScalar` | 戻り値が logical スカラー | 要 |

---

## 3.5 `toStruct` / `fromStruct`

**設計意図**: Mol オブジェクトを MATLAB セッションを超えて永続化するための変換関数（ADR-004）。

**永続化フォーマット**:
- **MolBlock (V2000)**: `MolToMolBlock()` → string。人間可読、SDF 相互運用可
- **Pickle**: Python ネイティブシリアライゼーション。高速・情報完全。MATLAB Online 非推奨

**ADR-004 の選択根拠**: MolBlock は ASCII テキストのため MATLAB の `.mat` ファイルに string として保存可能。Pickle バイト列は MATLAB 環境をまたいだ互換性が低い。

---

## 3.6 `toTable`

**設計意図**: 分子リストを MATLAB `table` に変換する Aggregate Table Builder（ADR-002 rev.4）。

**アルゴリズム概要**:
1. `mols` が非空 cell かつ全要素が Python オブジェクトであることを確認
2. `batchCalculate` でバッチ記述子計算（MW, LogP, TPSA, HBD, HBA, RotBonds）
3. 各 mol に `toSmiles` を適用して SMILES 列を追加
4. MATLAB `table()` で組み立てて返す

**テスト戦略** (`tests/unit/TestToTable.m`):
- 出力が `table` 型であることの確認
- 必須カラム（SMILES, MolWeight, LogP, TPSA, HBD, HBA, RotBonds）の存在確認
- 行数が mol 数と一致することの確認
- 型検証（各カラムが string/double の期待型）

---

## 3.7 `scaffold`

**設計意図**: Bemis-Murcko スキャフォールドを抽出し、分子骨格の SAR 分析・クラスタリングに使用する。

**アルゴリズム概要**:
1. `isa(mol, "py.rdkit.Chem.rdchem.Mol")` で型確認（RDKit 不要）
2. `py.importlib.import_module("rdkit.Chem.Scaffolds.MurckoScaffold")` で動的インポート
3. `murcko.MurckoScaffold.GetScaffoldForMol(mol)` でスキャフォールドを取得
4. `scaffold.GetNumAtoms() == 0` の場合（非環状分子）→ `toSmiles` で空 SMILES を返す（エラーにしない）
5. `py.rdkit.Chem.rdchem.Mol` オブジェクトを返す

**Bemis-Murcko スキャフォールドの定義**:
分子から側鎖・リンカーを除去し、環系とそれらを繋ぐリンカーだけを残したコアフレームワーク。

**参照テストデータ**:

| 入力分子 | 期待スキャフォールド | 根拠 |
|---|---|---|
| アスピリン (`CC(=O)Oc1ccccc1C(=O)O`) | ベンゼン環 (`c1ccccc1`) | 側鎖（アセテート・カルボキシル）を除去 |
| ピリジン (`c1ccncc1`) | ピリジン自身 | ヘテロ芳香環は窒素ごと保持 |
| ナフタレン (`c1ccc2ccccc2c1`) | ナフタレン自身（2環保持）| 縮合環は分割しない |
| エタノール (`CCO`) | `""` 空スキャフォールド | 非環状分子 → 0原子 → 空 SMILES |

**引用文献**:
- Bemis, G.W. & Murcko, M.A. (1996). *J. Med. Chem.* 39(15):2887–2893. DOI: 10.1021/jm9602928

**テスト戦略** (`tests/unit/TestMol.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_scaffold_numericInput_throwsInvalidInput` | double → `invalidInput` | 不要 |
| `test_scaffold_emptyInput_throwsInvalidInput` | `[]` → `invalidInput` | 不要 |
| `test_scaffold_aspirin_scaffoldIsNotNoneType` | アスピリン → py.rdkit.* | 要 |
| `test_scaffold_aspirin_returnsValidMol` | `isa(s, "py.rdkit.Chem.rdchem.Mol")` | 要 |
| `test_scaffold_aspirin_scaffoldSmiles_isBenzene` | アスピリン骨格 SMILES = `"c1ccccc1"` ★ | 要 |
| `test_scaffold_pyridine_nitrogenPreserved` | ピリジン窒素が骨格に保持される ★ | 要 |
| `test_scaffold_naphthalene_twoRingsPreserved` | ナフタレン 2 環ともスキャフォールドに含まれる ★ | 要 |
| `test_scaffold_ethanol_zeroAtomMol` | エタノール（非環状）→ GetNumAtoms=0 | 要 |
| `test_scaffold_linker_preserved` | リンカー接続分子 → リンカーが保持される | 要 |
| `test_scaffold_aspirin_getNumAtoms_crossValidation` | MATLAB 計算 == `py.GetNumAtoms()` ★★ | 要 |

**★ スキャフォールド SMILES クロス検証の重要性**:
`test_scaffold_aspirin_scaffoldSmiles_isBenzene` は `toSmiles(scaffold(mol))` の結果を
`"c1ccccc1"` と正確に比較する。MurckoScaffold の API 引数が入れ替わる場合（`GetGenericScaffold` の誤用等）をこのテストが捕捉する。
