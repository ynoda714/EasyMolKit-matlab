# Algorithm Guide — emk.io

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## 7.1 `readSdf` / `writeSdf`

**設計意図**: SDF（Structure-Data File）ファイルの読み書きを提供する。
RDKit の `SDMolSupplier` / `SDWriter` を使用し、Python 参照を MATLAB cell 配列として返す（ADR-002 準拠）。

**アルゴリズム概要 — `readSdf`**:
1. filePath が string/char であることを確認（RDKit 不要）
2. `isfile(filePath)` でファイル存在確認（RDKit 不要）
3. MATLAB `fileread(filePath)` で SDF 全体を文字列として読み込む
4. `strsplit(sdfText, '$$$$')` で mol ブロック列に分割する
5. 空ブロック（末尾の空要素）は `strtrim` 後にスキップ
6. 各ブロックを `Chem.MolFromMolBlock(char(blk))` で Mol に変換（IPC-safe な 1 回呼び出し）
7. `py.NoneType` の要素（パース失敗分子）をフィルタアウトし `logWarn` を出力
8. 有効な分子の **1×N** row cell 配列を返す

**OutOfProcess 制約と MolFromMolBlock 方式の採用根拠**:
`SDMolSupplier` は Boost.Python C++ 拡張型のイテラブルを返す。MATLAB OutOfProcess IPC がこのオブジェクトを
デシリアライズしようとして pickle に失敗し `クラス 'py.Boost.Python.class' が見つかりません` エラーが発生する。
MATLAB 側でテキスト分割してから `MolFromMolBlock()` を 1 分子ずつ呼ぶことで問題を回避する。
1 分子返し関数は `fromSmiles` と同じパターンで IPC-safe。

**RDKit デフォルト動作 (`readSdf`)**:
`MolFromMolBlock` は `sanitize=True, removeHs=True` で呼び出される（RDKit のデフォルト値）。
サニタイゼーションが通らない分子は `None` として返され、`logWarn` でスキップされる。

**アルゴリズム概要 — `writeSdf`**:
1. filePath が string/char であることを確認（RDKit 不要）
2. mols が cell かつ非空であることを確認（RDKit 不要）
3. 各 `mols{i}` が `py.*` で始まることを確認（RDKit 不要）
4. `fileparts(filePath)` で親ディレクトリを確認。存在しない場合は `dirNotFound` エラー
5. MATLAB `fopen(filePath, 'w')` でファイルを開く
6. 各 mol を `Chem.MolToMolBlock(mols{i})` で mol ブロック文字列に変換（→ IPC 通過可な plain string）
7. `fwrite(fid, [mb, "$$$$\n"], 'char')` で SDF レコードを書き出す
8. 例外発生時はファイルをクローズしてから `rdkitError` に変換して再スロー

**OutOfProcess 制約と MolToMolBlock 方式の採用根拠**:
`SDWriter.write(mol)` を MATLAB の OutOfProcess モードから呼び出す場合、MATLAB は mol オブジェクトを
IPC 境界越しに渡すために pickle シリアライズを試みる。しかし RDKit の `rdchem.Mol` は Boost.Python
C++ 拡張型のため pickle 不可能であり、MATLAB が `クラス 'py.Boost.Python.class' が見つかりません`
エラーを送出する。  
`MolToMolBlock(mol)` は Python 文字列（pickle 可能）を返すためこの問題を回避できる。

**パラメータ選択の根拠**:
- `char()` 変換: MATLAB string → Python str はカーネルバージョンによって自動変換されない場合があるため明示的な変換を行う
- SDF レコード終端 `$$$$\n`: MDL SDF フォーマット仕様に従い各 mol ブロックの末尾に付加する

**引用文献**:
- MDL/Dassault SDF Format Specification (Accelrys, 2011). *CTfile Formats*.
- RDKit Documentation: [SDMolSupplier](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.SDMolSupplier)
- RDKit Documentation: [SDWriter](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html#rdkit.Chem.rdmolfiles.SDWriter)

**テスト戦略** (`tests/unit/TestIo.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_readSdf_numericPath_throwsInvalidInput` | numeric filePath → `invalidInput` | 不要 |
| `test_readSdf_nonExistentFile_throwsFileNotFound` | 存在しないファイル → `fileNotFound` | 不要 |
| `test_readSdf_fileNotFound_errorMessage_containsPath` | エラーメッセージにパスが含まれる | 不要 |
| `test_readSdf_charPath_accepted` | char パスを受け入れる ★ | 要 |
| `test_readSdf_validSdf_returnsCell` | 有効 SDF → cell 配列 | 要 |
| `test_readSdf_result_isRowCell` | 戻り値が 1×N row cell ★ | 要 |
| `test_readSdf_roundTrip_countMatches_single` | round-trip: 分子数保持（1件） | 要 |
| `test_readSdf_roundTrip_countMatches_multiple` | round-trip: 分子数保持（3件） | 要 |
| `test_readSdf_roundTrip_smilesPreserved_ethanol` | round-trip: エタノール SMILES 一致 ★ | 要 |
| `test_readSdf_roundTrip_smilesPreserved_aspirin` | round-trip: アスピリン SMILES 一致 ★ | 要 |
| `test_readSdf_roundTrip_smilesPreserved_allThree` | round-trip: 3分子全 SMILES 一致 ★ | 要 |
| `test_writeSdf_emptyCell_throwsInvalidInput` | 空 cell → `invalidInput` | 不要 |
| `test_writeSdf_nonPythonElement_throwsInvalidMol` | 非 Python 要素 → `invalidMol` | 不要 |
| `test_writeSdf_nonExistentParentDir_throwsDirNotFound` | 存在しない親ディレクトリ → `dirNotFound` | 要 |
| `test_writeSdf_singleMol_fileIsCreated` | 1分子 write → ファイル存在 | 要 |
| `test_writeSdf_multipleMols_fileIsCreated` | 3分子 write → ファイル存在 | 要 |

**★ ラウンドトリップ検証テストの重要性**:
`test_readSdf_roundTrip_smilesPreserved_*` は `writeSdf` + `readSdf` の統合検証であり、
`SDWriter.write()` が不正な MolBlock を出力する場合、`SDMolSupplier` が読み込んだ分子を
`toSmiles()` で変換した際に SMILES が変化する場合、`char()` 変換漏れによるパス処理エラーを確実に検出する。

---

## 7.2 `readSmilesList`

**設計意図**: SMILES リストファイル（1 行 1 分子のプレーンテキスト）を読み込み、Mol オブジェクトの
cell 配列として返す。SDF に比べて軽量なフォーマットで文献・データベースからの相互運用性が高い。

**アルゴリズム概要**:
1. filePath が string/char であることを確認（RDKit 不要）
2. `isfile(filePath)` でファイル存在確認（RDKit 不要）
3. `readlines(filePath)` で全行を MATLAB string 配列として読み込む
4. 各行を `strtrim` して先頭・末尾空白を除去
5. 空行（`strlength == 0`）・コメント行（`startsWith(raw, "#")`）をスキップ（RDKit 不要）
6. `strsplit(raw, {" ", "\t"}, "CollapseDelimiters", true)` で最初のトークンを SMILES として抽出
7. `emk.mol.fromSmiles(smi)` を呼び出して Mol 変換
8. 変換失敗時は `logWarn` を出力してスキップ（例外を再スローしない）
9. 有効な分子の **1×N** row cell 配列を返す

**コメント行スキップを RDKit 呼び出し前に行う根拠**:
`#` で始まる行を先にフィルタすることで、不正な入力が `fromSmiles` に渡って `rdkitError` を発生させるのを防ぐ。
「RDKit 不要テストを先に書く」原則に沿った設計。

**`strsplit` による名前列無視の根拠**:
SMILES + 名前列（例: `"CCO\tEthanol"`, `"CCO Ethanol"`）は化学情報データベースの標準的な出力形式であるため、名前列を自動的に除去する設計が利便性を高める。

**引用文献**:
- Weininger, D. (1988). SMILES, a chemical language and information system. *JCICS* 28(1):31–36.
- Daylight SMILES format: [Daylight Theory Manual](https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html)

**テスト戦略** (`tests/unit/TestIo.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_readSmilesList_numericPath_throwsInvalidInput` | numeric filePath → `invalidInput` | 不要 |
| `test_readSmilesList_nonExistentFile_throwsFileNotFound` | 存在しないファイル → `fileNotFound` | 不要 |
| `test_readSmilesList_emptyFile_returnsEmptyCell` | 空ファイル → `{}` (RDKit 呼び出しなし) | 不要 |
| `test_readSmilesList_commentOnlyFile_returnsEmptyCell` | コメントのみ → `{}` (RDKit 呼び出しなし) | 不要 |
| `test_readSmilesList_validFile_returnsCell` | 有効ファイル → cell | 要 |
| `test_readSmilesList_result_isRowCell` | 1×N row cell ★ | 要 |
| `test_readSmilesList_countMatches_multiple` | 3 行 → 3 分子 | 要 |
| `test_readSmilesList_smiles_crossValidation_ethanol` | エタノール SMILES 一致 ★ | 要 |
| `test_readSmilesList_tabSeparatedName_isIgnoredAndMolParsed` | タブ区切り名前列を無視してパース ★ | 要 |
| `test_readSmilesList_commentLinesSkipped` | コメント行が分子に含まれない | 要 |
| `test_readSmilesList_invalidSmiles_skippedWithoutThrow` | 不正 SMILES → 例外なし・スキップ | 要 |
| `test_readSmilesList_mixedContent_correctCount` | コメント+空行+有効+無効の混在 → 2分子 ★ | 要 |
| TC29: `test_readSmilesList_allLinesFailed` | 有効行ありで全 SMILES パース失敗 → `allLinesFailed` | 要 |

**★ 混在コンテンツテストの重要性**:
`test_readSmilesList_mixedContent_correctCount` は 6 種類の入力パターン（ヘッダー・空行・有効×2・
無効・コメント）を一度に検証する統合テスト。フィルタリング・スキップロジックの全パスを単一テストで確認できる。
