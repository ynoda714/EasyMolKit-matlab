# Algorithm Guide — emk.cluster / emk.diversity

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## `emk.cluster.butina`

**設計意図**: Taylor-Butina 球排除クラスタリングで分子ライブラリを類似グループに分割する。各分子は唯一のクラスタに所属し、クラスタ重心（セントロイド）は最も多くの近傍を持つ分子。

**アルゴリズム概要**:
1. 下三角距離リストを `BulkTanimotoSimilarity` で構築（O(N²) の N²/2 ペア）
   ```
   distList には d(1,0), d(2,0), d(2,1), d(3,0), ... の順で格納
   ```
   累積 Python リスト (`accFpList`) を使って `i` 番目のFP と先行 `1..i-1` の類似度を一括取得
2. `ML.Cluster.Butina.ClusterData(distList, N, threshold, isDistData=True)` を呼出
3. 返り値の Python tuple of tuples を MATLAB cell of double arrays に変換（0-based → 1-based）

**Threshold パラメータ**:
- `Threshold` = **Tanimoto 距離**閾値（= 1 − Tanimoto 類似度）
- `Threshold=0.2` → 80% 以上類似の分子が同クラスタ（デフォルト）
- 小さいほどクラスタが大きく・数が少ない（厳格な類似性要件）
- 大きいほどクラスタが小さく・数が多い（緩い類似性要件）

**計算複雑度**: O(N²) — N=10,000 では数秒、N=100,000 では数分かかることがある。

**引用文献**:
- Butina, D. (1999). Unsupervised Data Base Clustering Based on Daylight's Fingerprint and Tanimoto Similarity. *J. Chem. Inf. Comput. Sci.* 39(4):747-750. DOI: 10.1021/ci9803381
- Taylor, R. (1995). Simulation Analysis of Experimental Design Strategies for Screening Random Compounds. *J. Chem. Inf. Comput. Sci.* 35(1):59-67. DOI: 10.1021/ci00023a009
- RDKit: `rdkit.ML.Cluster.Butina.ClusterData`

---

## `emk.diversity.pick`

**設計意図**: MaxMin アルゴリズム (Kennard-Stone の化学情報学的変形) で分子ライブラリから最大多様性サブセットを選択する。スクリーニングライブラリ設計・訓練セット選択に使用。

**アルゴリズム概要**:
1. Python `list` に全 FP を格納
2. `MaxMinPicker.LazyBitVectorPick(pyFps, M, N, seed=seed)` を呼出
   - `M`: ライブラリサイズ
   - `N`: 選択する多様分子数
   - `seed`: ランダム初期分子の固定（-1 = ランダム）
3. 返り値の Python tuple（0-based）を MATLAB double 配列（1-based）に変換

**MaxMin アルゴリズムの動作**:
1. ランダム（or seed 固定）で最初の分子を選択
2. 残り分子のうち、**既選択集合への最小距離が最大**の分子を次に選択
3. N 分子が選ばれるまで繰り返す

**LazyBitVectorPick の効率化**:
RDKit の `LazyBitVectorPick` は BitVector に最適化されており、距離関数の Python コールバックより高速。内部で `BulkTanimotoSimilarity` を活用している。

**引用文献**:
- Kennard, R.W. & Stone, L.A. (1969). Computer Aided Design of Experiments. *Technometrics* 11(1):137-148. DOI: 10.1080/00401706.1969.10490666
- Ashton, M. et al. (2002). Identification of Diverse Database Subsets Using Property-Based and Fragment-Based Molecular Descriptions. *QSAR Comb. Sci.* 21(8):598-604. DOI: 10.1002/1521-3838(200211)21:8<598::AID-QSAR598>3.0.CO;2-U
- RDKit: `rdkit.SimDivFilters.rdSimDivPickers.MaxMinPicker`

**テスト戦略**:

| ファイル | カテゴリ | 内容 |
|---|---|---|
| `TestCluster.m` | TC-1: バリデーション | 非cell入力・空cell・閾値範囲・Metric検証 |
| `TestCluster.m` | TC-2: 統合 | 1分子→1クラスタ・全メンバー数=N・同一分子→同クラスタ |
| `TestDiversity.m` | TC-1: バリデーション | 非cell・N範囲・Metric検証 |
| `TestDiversity.m` | TC-2: 統合 | 返りサイズ=N・1-based・unique・seed再現性 |
