# Algorithm Guide — emk.viz

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## 8.1 `draw2d`

**設計意図**: RDKit の 2D 描画エンジンを使って分子の 2D 構造画像を生成し、MATLAB figure に表示する。
ユーザーが Python の PIL や BytesIO を一切意識せずに `emk.mol.fromSmiles("CCO")` で取得した
Mol オブジェクトをそのまま描画できる薄いラッパーを提供する。

**アルゴリズム概要**:
1. `isa(mol, "py.rdkit.Chem.rdchem.Mol")` で入力型を検証（RDKit 不要）
2. `opts.Width < 1 || opts.Height < 1` で画像サイズを検証（RDKit 不要）
3. `py.importlib.import_module("rdkit.Chem.AllChem")` で AllChem モジュールを動的インポート
4. `allchem.Compute2DCoords(mol)` で 2D 座標を付加（冪等。既存の場合は上書き）
5. `py.importlib.import_module("rdkit.Chem.Draw")` で Draw モジュールを動的インポート
6. `draw.MolToFile(mol, tmpFile, pyargs("size", imgSize))` で一時 PNG ファイルに直接書き出す
7. `imread(tmpFile)` で MATLAB uint8(H, W, 3) 配列として読み込む
8. `figure()` → `imshow(img)` → `axis off` → 任意で `title()` を設定
9. `delete(tmpFile)` で一時ファイルを削除

**`importlib` による動的インポートを使う根拠**:
`AllChem` と `Draw` は `rdkit.Chem` のサブモジュールであり、`py.rdkit.Chem.AllChem` のような
直接 MATLAB 名前空間アクセスでは MATLAB の Python インターフェースがサブモジュールを
自動的にロードしない場合がある。
`py.importlib.import_module("rdkit.Chem.AllChem")` は Python の明示的 import と同等であり、確実にロードできる。

**2D 座標付加の必要性**:
`emk.mol.fromSmiles()` が返す Mol オブジェクトには 2D/3D 座標が含まれない。
RDKit の `MolToFile` は座標がない場合、全原子を原点に配置して描画するため、
`Compute2DCoords` が必須となる。`Compute2DCoords` は冪等であるため同一 mol を
2 回 `draw2d` に渡しても安全（TC5c で検証）。

**MolToFile を MolToImage+BytesIO より優先する根拠**:
`MolToImage` + `py.io.BytesIO` + `uint8(buf.read())` の変換パスは、`py.bytes` オブジェクトの
MATLAB `uint8` への変換で不正な型変換が発生し PNG バイト列が破損することが判明。
`MolToFile(mol, filePath, ...)` はファイルシステム経由で PNG を受け渡しするため、確実に動作する。

**一時ファイルの使用理由**:
MATLAB の `imread` は直接バイト配列を受け付けない。
`string(tempname()) + ".png"` のパターンは MATLAB の標準的な一時ファイル生成方法であり、
ユニークなパスを保証する。例外発生時も try/catch で `delete(tmpFile)` を実行してハンドルリークを防ぐ。

**引用文献**:
- RDKit Documentation: [Drawing Molecules](https://www.rdkit.org/docs/GettingStartedInPython.html#drawing-molecules)
- RDKit Documentation: [AllChem.Compute2DCoords](https://www.rdkit.org/docs/GettingStartedInPython.html#working-with-2d-molecules)
- Riniker, S. & Landrum, G.A. (2015). *J. Chem. Inf. Model.* 55(12):2562–2574.

**テスト戦略** (`tests/unit/TestViz.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_draw2d_numericInput_throwsInvalidInput` | double 入力 → `invalidInput` | 不要 |
| `test_draw2d_stringInput_throwsInvalidInput` | string 入力 → `invalidInput` | 不要 |
| `test_draw2d_emptyInput_throwsInvalidInput` | `[]` 入力 → `invalidInput` | 不要 |
| `test_draw2d_zeroWidth_throwsInvalidInput` | Width=0 → `invalidInput`（サイズ範囲チェック）★ | 要 |
| `test_draw2d_negativeHeight_throwsInvalidInput` | Height=-1 → `invalidInput`（負値ガード）★ | 要 |
| `test_draw2d_ethanol_doesNotThrow` | エタノール → エラーなし | 要 |
| `test_draw2d_benzene_doesNotThrow` | ベンゼン（芳香族）→ エラーなし | 要 |
| `test_draw2d_aspirin_doesNotThrow` | アスピリン（複雑分子）→ エラーなし | 要 |
| `test_draw2d_returnsFigureHandle` | 戻り値が `matlab.ui.Figure` ★ | 要 |
| `test_draw2d_defaultSize_imageIsCorrectHeight` | デフォルト Height=300 → 300行 ★ | 要 |
| `test_draw2d_defaultSize_imageIsCorrectWidth` | デフォルト Width=300 → 300列 ★ | 要 |
| `test_draw2d_defaultSize_imageIsRGB` | 3チャンネル (RGB) ★ | 要 |
| `test_draw2d_defaultSize_imageIsUint8` | uint8 型 ★ | 要 |
| `test_draw2d_customWidth_imageWidthMatches` | Width=400 → 400列 ★ | 要 |
| `test_draw2d_customHeight_imageHeightMatches` | Height=200 → 200行 ★ | 要 |
| `test_draw2d_nonSquare_bothDimensionsCorrect` | Width=400, Height=150 → 150×400 両軸同時検証 ★★ | 要 |
| `test_draw2d_idempotency_calledTwice_doesNotThrow` | 同一 mol を 2 回 draw2d → エラーなし ★ | 要 |
| `test_draw2d_titleOption_setsAxesTitle` | `Title="Ethanol"` → axes title に設定 ★ | 要 |

**★ Width/Height range validation テストの重要性**:
`Width < 1` / `Height < 1` のチェックは RDKit より前に発火する（TC1b）。
無効なサイズが RDKit に渡って `rdkitError` として誤分類されるのを防ぐ。

**★★ 非正方形 cross-validation テストの重要性**:
`test_draw2d_nonSquare_bothDimensionsCorrect` は Width と Height を同時に異なる値で検証する。
`pyargs("size", py.tuple({int32(W), int32(H)}))` で W/H の引数順序が入れ替わるバグを検出できる。

**合成データの科学的根拠**:
- エタノール (`"CCO"`): 3 重原子・単純鎖構造 — PubChem CID 702
- ベンゼン (`"c1ccccc1"`): 芳香族 SMARTS 表記。平面環状の描画を確認 — PubChem CID 241
- アスピリン (`"CC(=O)Oc1ccccc1C(=O)O"`): 複数の官能基・環・エステル結合を含む複雑分子 — PubChem CID 2244
