# Algorithm Guide — emk.scaffold

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## `genericMurcko`

**設計意図**: Bemis-Murcko スキャフォールドの「Generic」形式（原子タイプ・結合次数を消去した最小骨格）を生成する。骨格多様性の定量比較や、原子/結合タイプに依存しない clustering に使用。

**アルゴリズム概要**:
1. `MurckoScaffold.GetScaffoldForMol(mol)` — 環系と連結リンカー部分のみ抽出（側鎖除去）
2. `MurckoScaffold.MakeScaffoldGeneric(scaffold)` — 全原子を炭素、全結合を単結合に変換

**`emk.mol.scaffold` との違い**:

| 関数 | 返り値の特性 |
|---|---|
| `emk.mol.scaffold` | Bemis-Murcko スキャフォールド（原子タイプ・結合次数保持） |
| `emk.scaffold.genericMurcko` | Generic スキャフォールド（全原子 C、全結合単結合） |

Generic スキャフォールドは骨格トポロジーのみを比較したい場合に使用する。

**設計上の決定（PLAN.md）**: `emk.scaffold.murcko()` は作成しない（`emk.mol.scaffold` と重複するため）。

**引用文献**:
- Bemis, G.W. & Murcko, M.A. (1996). The Properties of Known Drugs. 1. Molecular Frameworks. *J. Med. Chem.* 39(15):2887-2893. DOI: 10.1021/jm9602928
- RDKit: `rdkit.Chem.Scaffolds.MurckoScaffold`

---

## `brics`

**設計意図**: BRICS (Breaking of Retrosynthetically Interesting Chemical Substructures) アルゴリズムで分子を合成的に意味のある断片に分解する。フラグメントライブラリ構築・de novo 設計に使用。

**アルゴリズム概要**:
1. `BRICS.BRICSDecompose(mol)` — Python `frozenset` でフラグメント SMILES を返す
2. `py.list(frozenset)` で MATLAB から反復可能なリストに変換
3. SMILES を `string(1,N)` 配列として返す

**重要な特性**:
- 返り値の順序は**不定**（frozenset 起源）
- 結合切断点は `[4*]`, `[8*]` 等のダミー原子で表現される
- 非切断可能分子（小さすぎる、または BRICS ルール非適合）は入力 SMILES のまま返ることがある

**引用文献**:
- Degen, J. et al. (2008). On the Art of Compiling and Using 'Drug-Like' Chemical Fragment Spaces. *ChemMedChem* 3(10):1503-1507. DOI: 10.1002/cmdc.200800178
- RDKit: `rdkit.Chem.BRICS.BRICSDecompose`

---

## `rgroup`

**設計意図**: R-group 分解（コアに対する側鎖の体系的列挙）により構造活性相関 (SAR) 解析のための table を生成する。

**アルゴリズム概要**:
1. `coreSmiles` を `Chem.MolFromSmiles` でパース（`[*:1]` 等の R-map 付き SMARTS コアを推奨）
2. `RGroupDecomposition.RGroupDecompose([core], mols, asSmiles=true)` を呼出
   - `asSmiles=true`: 結果を SMILES 文字列で返す（Mol オブジェクトより軽量）
3. Python `dict` の key 抽出: `py.list(firstRow.keys())` → `string(keysList{k})` で列名取得
4. `Core`, `R1`, `R2`, ... 列を持つ table を組み立て
5. アンマッチインデックスを 0-based → 1-based に変換して返す

**R-map 付きコアの例**:
```
"[*:1]c1ccccc1[*:2]"  — ベンゼン環に2つのR基
"[*:1]N1CCN(CC1)[*:2]"  — ピペラジン環
```

**引用文献**:
- Sauer, W.H.B. & Schwarz, M.K. (2003). Molecular Shape Diversity of Combinatorial Libraries. *J. Chem. Inf. Comput. Sci.* 43(3):987-1003. DOI: 10.1021/ci025599w
- RDKit: `rdkit.Chem.rdRGroupDecomposition.RGroupDecompose`

**テスト戦略** (`tests/unit/TestScaffold.m`):

| カテゴリ | 内容 |
|---|---|
| TC-1〜3: バリデーション | 非Mol入力・無効コア・全アンマッチ (RDKit 不要) |
| TC-4: genericMurcko 統合 | アスピリン・キノリン・非環状分子 |
| TC-5: brics 統合 | アスピリン・イマチニブ・非環状分子 |
| TC-6: rgroup 統合 | ベンゼンコア・無効コア |
