# Algorithm Guide — emk.db

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## 11.1 `searchPubchem`

**設計意図**: PubChem PUG REST API を使って化合物名・CID・SMILES で化合物情報を検索し、
MATLAB table として返す。ユーザーが API の詳細（URL 構築・JSON パース・null 処理）を
意識せずに PubChem データベースを利用できる薄いラッパー。

**アルゴリズム概要**:
1. `query` が string/char であることを確認（RDKit 不要）
2. `SearchBy` の値確認 (`"name"`, `"cid"`, `"smiles"` のいずれか）（RDKit 不要）
3. `urlencode(char(query))` で URL エンコード
4. `webread(url)` で JSON を取得（1 回の HTTP GET のみ。`webwrite` は使用しない）
5. HTTP 404 の検出: `MException.message` に `"404"` が含まれる場合 → `notFound` エラー
6. JSON パースと null → `NaN`/`""` 変換
7. 5 カラムの table を返す

**出力スキーマ（5 カラム）**:

| カラム名 | 型 | ソースフィールド |
|---|---|---|
| `CID` | `double` | `cid` |
| `IUPACName` | `string` | `iupacname` |
| `MolecularFormula` | `string` | `molecularformula` |
| `MolecularWeight` | `double` | `molecularweight` |
| `IsomericSMILES` | `string` | `isomericsmiles` |

**PUG REST API エンドポイント設計**:
```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{SearchBy}/{encodedQuery}/JSON?
  properties=IUPACName,MolecularFormula,MolecularWeight,IsomericSMILES
```

**`webread` のみ使用の根拠**:
HTTP GET 操作に `webread` で十分。`webwrite` は POST/PUT 用であり、
データ取得に使うと予期しないペイロードが生成される可能性がある。

**`urlencode` の必要性**:
SMILES クエリに `+`, `#`, `%` などの特殊文字が含まれる場合（例: `CC(=O)Oc1ccccc1C(=O)O`）、
URL に直接埋め込むと API が不正なリクエストと解釈する。`urlencode` でパーセントエンコードする。

**404 検出ロジックの根拠**:
MATLAB の `webread` は HTTP 404 を `MException` としてスローするが、
identifier は `"MATLAB:webservices:HTTP404StatusCodeError"` など実装依存の場合がある。
`MException.message` に `"404"` が含まれるかを文字列検索で判定する方が、MATLAB バージョン間の互換性が高い。

**引用文献**:
- PubChem PUG REST API Documentation: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest
- Kim, S. et al. (2021). PubChem in 2021: new data content and improved web interfaces. *Nucleic Acids Research* 49(D1):D1388–D1395.

**テスト戦略** (`tests/unit/TestDb.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_searchPubchem_numericQuery_throwsInvalidInput` | double → `invalidInput` | 不要 |
| `test_searchPubchem_invalidSearchBy_throwsInvalidSearchBy` | 未知 SearchBy → `invalidSearchBy` | 不要 |
| `test_searchPubchem_nonExistentCompound_throwsNotFound` | `"xyznonexistent99999"` → `notFound` ★★ | 要 |
| `test_searchPubchem_returnsTable` | 戻り値が table | 要 |
| `test_searchPubchem_ethanol_hasCIDColumn` | `CID` カラム存在 | 要 |
| `test_searchPubchem_ethanol_hasIUPACNameColumn` | `IUPACName` カラム存在 | 要 |
| `test_searchPubchem_ethanol_hasMolecularFormulaColumn` | `MolecularFormula` カラム存在 | 要 |
| `test_searchPubchem_ethanol_hasMolecularWeightColumn` | `MolecularWeight` カラム存在 | 要 |
| `test_searchPubchem_ethanol_hasIsomericSMILESColumn` | `IsomericSMILES` カラム存在 | 要 |
| `test_searchPubchem_ethanol_CIDIsDouble` | `CID` が double 型 ★ | 要 |
| `test_searchPubchem_ethanol_MW_isDouble` | `MolecularWeight` が double 型 ★ | 要 |
| `test_searchPubchem_ethanol_CIDis702` | エタノール CID = 702 ★★★ | 要 |
| `test_searchPubchem_ethanol_MW_inExpectedRange` | MW ∈ [40, 50] ★★★ | 要 |
| `test_searchPubchem_ethanol_formulaIsC2H6O` | MolecularFormula = `"C2H6O"` ★★★ | 要 |
| `test_searchPubchem_charQuery_accepted` | char クエリも受け入れ ★ | 要 |
| `test_searchPubchem_searchByCID_ethanol_CIDis702` | SearchBy="cid", query="702" → CID=702 | 要 |
| `test_searchPubchem_searchBySmiles_ethanol_CIDis702` | SearchBy="smiles", query="CCO" → CID=702 | 要 |

**★★ notFound テストの重要性**:
`test_searchPubchem_nonExistentCompound_throwsNotFound` は HTTP 404 の検出ロジックを直接検証する。
`MException.message` に `"404"` が含まれるかの判定が MATLAB バージョン間でブレがないかを確認する。

**★★★ 参照データの根拠**:
- エタノール CID = 702 は PubChem の公式 CID（安定した well-known ID）
- 分子式 `"C2H6O"` は IUPAC 表記（順序: C → H → その他のアルファベット順）
- MW ∈ [40, 50] の範囲チェックはエタノール MW = 46.07 を含む緩い検証（API の有効数字のばらつきを許容）

---

## 11.2 `searchChembl`

**設計意図**: ChEMBL REST API を使って化合物名・ChEMBL ID で化合物情報を検索し、
MATLAB table として返す。PubChem より詳細な薬理活性データが付属する化合物データベース。

**アルゴリズム概要**:
1. `query` が string/char、`SearchBy` が `"name"` または `"chemblid"` であることを確認（RDKit 不要）
2. SearchBy で URL パターンを切り替え:
   - `"name"`: `/molecule.json?pref__icontains={encoded}&limit=25`（部分一致検索）
   - `"chemblid"`: `/molecule/{CHEMBL_ID}.json`（完全一致検索）
3. `webread(url)` で JSON を取得
4. `molecules` 配列から各フィールドを抽出。null フィールドは `extractNumericField` ヘルパーで `NaN`/`""` に変換
5. 8 カラムの table を返す

**出力スキーマ（8 カラム）**:

| カラム名 | 型 | ChEMBL フィールド |
|---|---|---|
| `ChemblID` | `string` | `molecule_chembl_id` |
| `PreferredName` | `string` | `pref_name` |
| `MolecularFormula` | `string` | `molecule_properties.full_molformula` |
| `MolecularWeight` | `double` | `molecule_properties.full_mwt` |
| `AlogP` | `double` | `molecule_properties.alogp` |
| `HBD` | `double` | `molecule_properties.hbd` |
| `HBA` | `double` | `molecule_properties.hba` |
| `SMILES` | `string` | `molecule_structures.canonical_smiles` |

**`/molecule/{id}.json` vs `/molecule.json?filter` の切り替え根拠**:
ChEMBL ID（例: `"CHEMBL25"`）は完全一致エンドポイントを使う方が高速かつ確実。
名前検索は部分一致が必要なため、フィルタ API（`pref__icontains`）を使用する。
`limit=25` は名前検索の結果数上限（デフォルトでは 20）を明示的に設定する。

**`extractNumericField` ヘルパーの役割**:
ChEMBL API の数値フィールドは文字列として返される場合がある（例: `"329.37"` as string）。
`str2double` で変換し、失敗時は `NaN` を返す。型が `py.NoneType` の場合も `NaN` に変換する。

**引用文献**:
- Zdrazil, B. et al. (2024). The ChEMBL database in 2023: a drug discovery platform spanning multiple bioactivity complementary data sources. *Nucleic Acids Research* 52(D1):D1180–D1192.
- ChEMBL REST API Documentation: https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services

**テスト戦略** (`tests/unit/TestDb.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_searchChembl_numericQuery_throwsInvalidInput` | double → `invalidInput` | 不要 |
| `test_searchChembl_invalidSearchBy_throwsInvalidSearchBy` | 未知 SearchBy → `invalidSearchBy` | 不要 |
| `test_searchChembl_nonExistentChemblId_throwsNotFound` | `"CHEMBL999999999"` → `notFound` ★★ | 要 |
| `test_searchChembl_returnsTable` | 戻り値が table | 要 |
| `test_searchChembl_aspirin_hasChemblIDColumn` | `ChemblID` カラム存在 | 要 |
| `test_searchChembl_aspirin_hasMolecularFormulaColumn` | `MolecularFormula` カラム存在 | 要 |
| `test_searchChembl_aspirin_hasMolecularWeightColumn` | `MolecularWeight` カラム存在 | 要 |
| `test_searchChembl_aspirin_hasAlogPColumn` | `AlogP` カラム存在 | 要 |
| `test_searchChembl_aspirin_hasSMILESColumn` | `SMILES` カラム存在 | 要 |
| `test_searchChembl_aspirin_ChemblIDis25` | ChemblID = `"CHEMBL25"` ★★★ | 要 |
| `test_searchChembl_aspirin_MW_inExpectedRange` | MW ∈ [170, 185] ★★★ | 要 |
| `test_searchChembl_aspirin_formula_isC9H8O4` | MolecularFormula = `"C9H8O4"` ★★★ | 要 |
| `test_searchChembl_searchByChemblid_aspirin` | SearchBy="chemblid", "CHEMBL25" → ChemblID="CHEMBL25" | 要 |
| `test_searchChembl_MolecularWeight_isDouble` | `MolecularWeight` が double 型 ★ | 要 |

**★★★ 参照データの根拠**:
- アスピリン ChEMBL ID = `"CHEMBL25"` は ChEMBL データベースの安定した well-known ID
- 分子式 `"C9H8O4"` は CAS RN 50-78-2 のアセチルサリチル酸の IUPAC 表記
- MW ∈ [170, 185] の範囲チェックはアスピリン MW = 180.16 を含む緩い検証（`extractNumericField` の str2double 変換精度を許容）

---

## 11.3 `searchChemblTarget`

**設計意図**: ChEMBL REST API の target エンドポイントを使ってタンパク質ターゲットを名前で検索し、
`TargetChEMBLID` を含む MATLAB table として返す。ユーザーは ChEMBL ID を事前知識なしに
発見できる。Python 不要（`webread` のみ）。

**アルゴリズム概要**:
1. `query` が string/char であることを確認
2. `pref_name__icontains={encoded}` フィルタで ChEMBL target エンドポイントに GET リクエスト
3. `TargetType` が非空の場合、`target_type={encoded}` フィルタを追加（デフォルト: `"SINGLE PROTEIN"`）
4. `webread` で JSON を取得し、`targets` 配列を struct array に正規化
5. 4 カラムの table を返す

**出力スキーマ（4 カラム）**:

| カラム名 | 型 | ChEMBL フィールド |
|---|---|---|
| `TargetChEMBLID` | `string` | `target_chembl_id` |
| `PreferredName`  | `string` | `pref_name` |
| `Organism`       | `string` | `organism` |
| `TargetType`     | `string` | `target_type` |

**ChEMBL target エンドポイント**:
```
GET https://www.ebi.ac.uk/chembl/api/data/target.json
  ?pref_name__icontains={query}&target_type={TargetType}&limit={MaxRows}
```

**`pref_name__icontains` を選択した根拠**:
ChEMBL の `target_synonym` エンドポイントは関係テーブルの JOIN が必要で REST API では直接サポートされない。
`pref_name__icontains` は preferred name の部分一致で十分な発見性を提供する。
遺伝子シンボル（例: "EGFR"）は preferred name（例: "Epidermal growth factor receptor erbB1"）に
含まれないため、ユーザーには full name または部分 name で検索するよう案内する。

**引用文献**:
- ChEMBL REST API: https://chembl.gitbook.io/chembl-interface-documentation/web-services/chembl-data-web-services

**テスト戦略** (`tests/unit/TestDb.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_searchChemblTarget_numericQuery_throwsInvalidInput` | double → `invalidInput` | 不要 |
| `test_searchChemblTarget_emptyQuery_throwsInvalidInput` | `""` → `invalidInput` | 不要 |
| `test_searchChemblTarget_invalidMaxRows_throwsInvalidOptions` | `MaxRows=0` → `invalidOptions` | 不要 |
| `test_searchChemblTarget_returnsTable` | 戻り値が table | 要（ネット） |
| `test_searchChemblTarget_EGFR_hasTargetChEMBLIDColumn` | `TargetChEMBLID` カラム存在 | 要 |
| `test_searchChemblTarget_EGFR_containsCHEMBL203` | "Epidermal growth factor receptor" → CHEMBL203 | 要 |

---

## 11.4 `getChemblActivity`

**設計意図**: ChEMBL REST API の activity エンドポイントを使って、指定ターゲットの生物活性データ（IC50 等）
をダウンロードし、有効 SMILES を持つ nM 単位のレコードのみを MATLAB table として返す。
Python 不要（`webread` のみ）。SDF 保存・類似度解析との連携を想定した設計。

**アルゴリズム概要**:
1. `targetId` が string/char、`MaxRows` が正整数であることを確認
2. `standard_relation=%3D`（完全測定値のみ）でフィルタした URL を構築
3. `webread` で JSON を取得し、`activities` 配列を正規化
4. 各レコードから `molecule_chembl_id`, `molecule_pref_name`, `canonical_smiles`, `standard_type`, `standard_value`, `standard_units` を抽出
5. `standard_value` は JSON では文字列で返されるため `str2double` で変換
6. `canonical_smiles` が有効、`standard_value` が非 NaN、`standard_units == "nM"` の行のみ保持
7. `MinActivity_nM` 閾値を適用（`Value_nM <= MinActivity_nM`）
8. 5 カラムの table を返す

**出力スキーマ（5 カラム）**:

| カラム名 | 型 | ChEMBL フィールド |
|---|---|---|
| `MoleculeChEMBLID` | `string` | `molecule_chembl_id` |
| `Name`             | `string` | `molecule_pref_name` |
| `SMILES`           | `string` | `canonical_smiles` |
| `ActivityType`     | `string` | `standard_type` |
| `Value_nM`         | `double` | `standard_value`（nM 変換済み） |

**ChEMBL activity エンドポイント**:
```
GET https://www.ebi.ac.uk/chembl/api/data/activity.json
  ?target_chembl_id={targetId}&standard_type={ActivityType}
  &standard_relation=%3D&limit={MaxRows}
```

**`standard_relation=%3D` フィルタの根拠**:
`<` や `>` の関係記号は「測定値の上界・下界」を意味する（例: IC50 < 1000 nM）。
これらは正確な測定値でないため QSAR 解析に不適。`=`（URL エンコード: `%3D`）
のみを取得することで測定値の信頼性を確保する。

**nM フィルタの根拠**:
ChEMBL の `standard_value` は ChEMBL が標準化した nM 値であるが、
uM（マイクロモル）や mg/mL など変換されていないレコードが混在する場合がある。
`standard_units == "nM"` でフィルタすることで単位の一貫性を保証する。

**引用文献**:
- Mendez D et al. (2019) ChEMBL: towards direct deposition of bioassay data. *Nucleic Acids Res* 47:D930-D940.
- IC50 概念: Cheng Y & Prusoff WH (1973) Relationship between the inhibition constant (Ki) and the concentration of inhibitor which causes 50% inhibition (I50) of an enzymatic reaction. *Biochem Pharmacol* 22:3099-3108.

**テスト戦略** (`tests/unit/TestDb.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_getChemblActivity_numericTargetId_throwsInvalidInput` | double → `invalidInput` | 不要 |
| `test_getChemblActivity_emptyTargetId_throwsInvalidInput` | `""` → `invalidInput` | 不要 |
| `test_getChemblActivity_invalidMaxRows_throwsInvalidOptions` | `MaxRows=0` → `invalidOptions` | 不要 |
| `test_getChemblActivity_returnsTable` | CHEMBL203 → table | 要（ネット） |
| `test_getChemblActivity_hasRequiredColumns` | 5 カラム全存在 | 要 |
| `test_getChemblActivity_Value_nM_isDouble` | `Value_nM` が double 型 | 要 |
| `test_getChemblActivity_SMILES_nonEmpty` | 全行 SMILES 非空 | 要 |
| `test_getChemblActivity_MinActivity_nM_filters` | `MinActivity_nM=10` → 全行 <= 10 | 要 |

---

## 11.5 `pubchemFetch`

**設計意図**: `searchPubchem`（webread のみ）では取得できない拡張プロパティ（同義語・InChI/InChIKey・
XLogP・TPSA・HBD/HBA・Complexity 等）を PubChemPy Python ライブラリ経由で取得し MATLAB struct で返す。

**前提条件**: `emk.setup.installExtra("pubchempy")` で `pubchempy` を Embedded Python に pip install 済みであること。

**アルゴリズム概要**:
1. 識別子の型チェック（string/char/numeric scalar → string 変換）
2. `NameSpace` バリデーション（"name"/"smiles"/"cid"/"inchi"/"inchikey"/"formula"）
3. `py.importlib.import_module("pubchempy")` で可用性確認 → 不可なら `libraryNotFound`
4. `pcp.get_compounds(identifier, ns)` を呼び出し → PubChemPy が PubChem PUG REST に HTTP リクエスト
5. 結果が空なら `notFound`。先頭要素（最高スコア）の属性をプライベートヘルパーで変換:
   - `pyStr_()`: Python str/None → MATLAB string（None → `""`）
   - `pyDbl_()`: Python numeric/None → MATLAB double（None → NaN）
   - `pyStrDbl_()`: PubChemPy の `molecular_weight` は文字列 → `str2double` で変換
6. `synonyms` は `cell(c.synonyms)` → 先頭 `MaxSynonyms` 件を string 配列として返す

**PubChemPy vs. PUG REST（`searchPubchem`）の比較**:

| 特性 | `searchPubchem` | `pubchemFetch` |
|---|---|---|
| Python 要否 | 不要 (webread) | 要 (PubChemPy) |
| 返値型 | table | struct |
| 対応識別子 | name/smiles/cid/inchikey | name/smiles/cid/inchi/inchikey/formula |
| 同義語 | × | ✓ (up to MaxSynonyms) |
| InChI/InChIKey | × | ✓ |
| XLogP | × | ✓ |
| TPSA | × | ✓ |
| HBD/HBA | × | ✓ |
| Complexity | × | ✓ |

**Python None の変換方針**:
PubChemPy の一部プロパティ（`xlogp`, `iupac_name` 等）は化合物によって `None` が返る。
MATLAB 側で `double(None)` を試みると例外が発生するため、`try/catch` で NaN/"" に変換する。
`molecular_weight` は PubChemPy が文字列（`"46.07"`）として返すため `str2double` を使用。

**引用文献**:
- PubChemPy Documentation: https://pubchempy.readthedocs.io/
- Kim, S. et al. (2021). PubChem in 2021. *Nucleic Acids Research* 49(D1):D1388–D1395.

**テスト戦略** (`tests/unit/TestPubchemFetch.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| TC1 非スカラー数値 → `invalidInput` | `[702, 241]` → 即時エラー | なし |
| TC2 セル配列 → `invalidInput` | `{"ethanol"}` → 即時エラー | なし |
| TC3 論理値 → `invalidInput` | `true` → 即時エラー | なし |
| TC4 空文字列 → `invalidInput` | `""` → 即時エラー | なし |
| TC5 空白のみ → `invalidInput` | `"   "` → 即時エラー | なし |
| TC6 char 入力受容 | char は `invalidInput` にならない | なし |
| TC7 数値スカラー CID 受容 | 数値スカラーは `invalidInput` にならない | なし |
| TC8 無効 NameSpace → `invalidNamespace` | 未知 NameSpace → 即時エラー | なし |
| TC9 NameSpace 大文字小文字非区別 | `"NAME"` は `invalidNamespace` にならない | なし |
| TC10 全有効 NameSpace 値 (TestParameter) | 6 種すべてが `invalidNamespace` にならない | なし |
| TC11 `libraryNotFound` (absent のみ) | pubchempy 不在時のみ実行 | なし（absent 環境） |
| TC12 名前検索 → struct | `"ethanol"` → struct 型 | pubchempy + ネット |
| TC13 全フィールド存在 | 16 フィールドすべて存在 | pubchempy + ネット |
| TC14 CID = 702 | エタノール CID が正確に 702 | pubchempy + ネット |
| TC14b CID 型・サイズ | `double` スカラーであること | pubchempy + ネット |
| TC15 MW 型・参照値 | `double` かつ 46.07 ± 0.5 g/mol ★★★ | pubchempy + ネット |
| TC16 IsomericSMILES 非空 | 非空 string | pubchempy + ネット |
| TC17 InChIKey 完全一致 | `"LFQSCWFLJHTTHZ-UHFFFAOYSA-N"` ★★★ | pubchempy + ネット |
| TC18 Synonyms 非空 | `numel >= 1` | pubchempy + ネット |
| TC18b Synonyms 行ベクトル | `size(s.Synonyms, 1) == 1` | pubchempy + ネット |
| TC19 SMILES 名前空間で CID 一致 | `"CCO"` SMILES → CID = 702 | pubchempy + ネット |
| TC20 CID 名前空間（数値入力） | `702` numeric → CID = 702 | pubchempy + ネット |
| TC21 存在しない化合物 → `notFound` | ランダム文字列 → `notFound` | pubchempy + ネット |
| TC22 MaxSynonyms 上限制御 | `MaxSynonyms=3` → `numel <= 3` | pubchempy + ネット |
| TC23 XLogP double スカラー | 型・サイズ検証 | pubchempy + ネット |
| TC24 TPSA 参照値 | 20.23 ± 0.5 Å² ★★★ | pubchempy + ネット |
| TC25 HBondDonors = 1 | エタノール OH → 1 ★★★ | pubchempy + ネット |
| TC25b HBondAcceptors = 1 | エタノール O → 1 ★★★ | pubchempy + ネット |
| TC25c RotatableBonds = 0 | エタノール → 0 ★★★ | pubchempy + ネット |
| TC25d Charge = 0 (double scalar) | 中性分子 → 0 | pubchempy + ネット |
| TC25e Complexity 有限 double scalar | Bertz CT > 0 for ethanol | pubchempy + ネット |
| TC25f MolecularFormula = "C2H6O" | IUPAC 式 ★★★ | pubchempy + ネット |
| TC25g IUPACName 非空 | 非空 string | pubchempy + ネット |
| TC26 HeavyAtomCount = 3 | 2C + 1O ★★★ | pubchempy + ネット |
| TC27 InChI 接頭辞 | `"InChI="` で始まる | pubchempy + ネット |
| TC28 searchPubchem とのクロスバリデーション | CID が PUG REST と一致 ★★★ | pubchempy + ネット |

**★★★ 参照データの根拠**:
- エタノール CID = 702、InChIKey = LFQSCWFLJHTTHZ-UHFFFAOYSA-N は PubChem 公式安定 ID
- MW ≈ 46.07 g/mol、TPSA ≈ 20.23 Å²、HBD = 1、HBA = 1、RotatableBonds = 0、HeavyAtomCount = 3 は PubChem CID 702 の公式値
- MolecularFormula = "C2H6O" は IUPAC Hill 表記（C → H → アルファベット順）
- TC28 クロスバリデーションの意義: PubChemPy（Python 経由）と PUG REST（webread 直接）の 2 つの独立した経路が同一 CID を返すことを保証する
- TC18b 行ベクトル保証: Python list → MATLAB `cell()` は 1×N セル配列を返す。`string(1×N cell)` は 1×N string 配列になることを API 仕様として検証
