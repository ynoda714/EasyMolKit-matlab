# Algorithm Guide — emk.conformer / emk.shape

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## `emk.conformer.embed`

**設計意図**: ETKDG アルゴリズムで分子の 3D コンフォーマーを生成する。後続の力場最適化 (`emk.conformer.optimize`) や 3D 形状比較 (`emk.shape.compare`) の前処理として使用。

**アルゴリズム概要**:
1. `AllChem.AddHs(mol)` — 明示的 H 追加（3D 幾何精度向上のため）
2. `AllChem.<Method>()` でパラメータオブジェクト取得
3. `params.randomSeed = seed` （seed≥0 の場合）
4. `AllChem.EmbedMolecule(molH, params)` — 0=成功, -1=失敗
5. `AllChem.RemoveHs(molH)` — H を削除して重原子のみの 3D 座標を保持

**サポートメソッド**:

| Method | アルゴリズム | 推奨用途 |
|---|---|---|
| `ETKDGv3` (デフォルト) | Cambridge Structural Database 由来のトーション分布を使用 | 一般的な創薬分子 |
| `ETKDGv2` | トーションライブラリの初期版 | レガシー互換 |
| `ETKDG` | 実験的トーション + 距離幾何 | 小分子 |
| `KDG` | 基本的な距離幾何 | 参照用 |

**埋め込み失敗への対処**:
- 高剛性・大環状・特殊元素の分子では埋め込み失敗率が高い
- 失敗時は `emk:conformer:embed:embeddingFailed` エラー（-1 を RDKit が返した場合）
- 異なる `RandomSeed` でリトライするか、より単純な `Method` を試す

**引用文献**:
- Wang, S. et al. (2020). Improving Conformer Generation for Small Molecules: Learning Torsional Distributions from the Cambridge Structural Database. *J. Chem. Inf. Model.* 60(4):2044-2058. DOI: 10.1021/acs.jcim.0c00025 (ETKDGv3)
- Riniker, S. & Landrum, G.A. (2015). Better Informed Distance Geometry. *J. Chem. Inf. Model.* 55(12):2562-2574. DOI: 10.1021/acs.jcim.5b00654 (ETKDG)
- RDKit: `rdkit.Chem.AllChem.EmbedMolecule`, `ETKDGv3`

---

## `emk.conformer.optimize`

**設計意図**: 力場最小化で 3D コンフォーマーの歪みを解消し、エネルギー的に安定な構造を得る。`emk.conformer.embed` 後に呼ぶことで形状比較精度が向上する。

**アルゴリズム概要**:
1. `AllChem.AddHs(mol, addCoords=true)` — 3D 座標付きで H 追加（位置を最適化に含める）
2. `AllChem.MMFFOptimizeMolecule(molH, maxIters=MaxIter)` または `UFFOptimizeMolecule`
   - status 0 = 収束成功
   - status 1 = 未収束（logWarn のみ）
   - status -1 = FF 設定失敗（`optimizeFailed` エラー）
3. `AllChem.RemoveHs(molH)` — H 削除（重原子の 3D 座標を保持）

**力場の選択**:

| ForceField | 適用対象 | 精度 |
|---|---|---|
| `MMFF94` (デフォルト) | C, H, N, O, S, P, F, Cl, Br, I | 高（創薬分子に推奨） |
| `UFF` | 全周期表元素 | 中（非標準元素のフォールバック） |

MMFF94 が FF 設定に失敗した場合（非標準元素など）は `optimizeFailed` が発生する。`ForceField="UFF"` へのフォールバックを検討すること。

**引用文献**:
- Halgren, T.A. (1996). Merck Molecular Force Field. I. *J. Comput. Chem.* 17(5-6):490-519. DOI: 10.1002/jcc.540170510 (MMFF94)
- Rappe, A.K. et al. (1992). UFF, a Full Periodic Table Force Field. *J. Am. Chem. Soc.* 114(25):10024-10035. DOI: 10.1021/ja00051a040 (UFF)

---

## `emk.shape.compare`

**設計意図**: 2 分子の 3D 形状を数値的に比較し、類似度スコアを返す。仮想スクリーニング・コンフォーマー品質評価に使用。

**サポートメソッド**:

### `protrude` (デフォルト)

`ShapeProtrudeDist(mol1, mol2, allowReordering=false)` を使用。

$$\text{score} = 1 - d_\text{protrude}$$

$d_\text{protrude}$ は mol1 体積のうち mol2 と重ならない割合（非対称）。`allowReordering=false` により mol1 を常に probe として固定。

### `tanimoto`

`ShapeTverskyIndex(mol1, mol2, α=1, β=1)` を使用。

$$\text{score} = \frac{|A \cap B|}{|A \cup B|}$$

形状 Tanimoto（Gaussian ボリューム重複による）。対称指標。

**数値クランプ**: スコアは浮動小数点誤差に備えて `max(0, min(1, score))` でクランプ済み。

**重要な注意**:
- 形状比較の精度はコンフォーマーの配向に大きく依存する
- 本関数は分子のアライメント（重ね合わせ最適化）を行わない
- 精度が必要な場合は外部ツール（ROCS, OpenEye Omega等）の使用を推奨

**引用文献**:
- Ballester, P.J. & Richards, W.G. (2007). Ultrafast Shape Recognition to Search Compound Libraries. *J. Comput. Chem.* 28(10):1711-1723. DOI: 10.1002/jcc.20681
- Grant, J.A. & Pickup, B.T. (1995). A Gaussian Description of Molecular Shape. *J. Phys. Chem.* 99(11):3503-3510. DOI: 10.1021/j100011a016
- RDKit: `rdkit.Chem.rdShapeHelpers.ShapeProtrudeDist`, `ShapeTverskyIndex`

**テスト戦略** (`tests/unit/TestConformer.m`):

| カテゴリ | 内容 |
|---|---|
| TC-1: embed バリデーション | 非Mol入力・無効Method |
| TC-2: optimize バリデーション | 非Mol・コンフォーマーなし・無効FF |
| TC-3: compare バリデーション | 非Mol・コンフォーマーなし・無効Method |
| TC-4: embed 統合 | Mol 返り値・コンフォーマー数・seed 再現性 |
| TC-5: optimize 統合 | MMFF94/UFF 動作・パイプライン |
| TC-6: compare 統合 | 同一分子→スコア≈1・スコア範囲[0,1] |
