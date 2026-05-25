# Algorithm Guide — emk.similarity

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## 6.1 `tanimoto`

**設計意図**: 2 つのフィンガープリント間の Tanimoto 係数を計算する。化学類似度検索のデファクト
スタンダードであり、0 = 完全不一致、1 = 完全一致。Morgan FP および MACCS FP の両方に適用可能。

**アルゴリズム概要**:
1. `startsWith(class(fp1), "py.")` で Python オブジェクトかチェック（RDKit 不要で発火するバリデーション）
2. `py.rdkit.DataStructs.TanimotoSimilarity(fp1, fp2)` を呼び出す
3. Python 例外（ビット長不一致等）→ `rdkitError` に変換（`try/catch`）
4. 返り値（Python float）を MATLAB `double` に変換して返す

**数学的定義**:

$$
T(A, B) = \frac{|A \cap B|}{|A \cup B|} = \frac{c}{a + b - c}
$$

ここで $a = |A|$（FP A の ON ビット数）, $b = |B|$（FP B の ON ビット数）, $c = |A \cap B|$（両者でオンのビット数）。
両方のベクトルが全ゼロのとき、RDKit は慣例として $T = 0.0$ を返す。

**バリデーション前置の根拠**: `py.*` チェックを RDKit 呼び出しより前に行うことで、RDKit 未起動セッションでも型エラーを確実に検出できる。

**引用文献**:
- Tanimoto, T.T. (1958). An Elementary Mathematical Theory of Classification and Prediction. IBM Internal Report.
- Willett, P., Barnard, J.M. & Downs, G.M. (1998). *J. Chem. Inf. Comput. Sci.* 38(6), 983–996. DOI: 10.1021/ci9800211

**テスト戦略** (`tests/unit/TestSimilarity.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_tanimoto_fp1_numericInput_throwsInvalidInput` | double fp1 → `invalidInput` | 不要 |
| `test_tanimoto_fp1_stringInput_throwsInvalidInput` | string fp1 → `invalidInput` | 不要 |
| `test_tanimoto_fp1_emptyInput_throwsInvalidInput` | `[]` fp1 → `invalidInput` | 不要 |
| `test_tanimoto_fp2_numericInput_throwsInvalidInput` | 有効 fp1 + double fp2 → `invalidInput` | 要 |
| `test_tanimoto_returnsDoubleScalar` | 戻り値が double スカラー | 要 |
| `test_tanimoto_returnIsFinite` | `isfinite(score)` | 要 |
| `test_tanimoto_sameFP_morgan_equals1` | 同一 Morgan FP → score = 1.0 | 要 |
| `test_tanimoto_sameFP_maccs_equals1` | 同一 MACCS FP → score = 1.0 | 要 |
| `test_tanimoto_ethanolVsAspirin_inRange` | スコア ∈ [0, 1] | 要 |
| `test_tanimoto_symmetry_ethanolVsAspirin` | T(A,B) = T(B,A)（対称性） | 要 |
| `test_tanimoto_ethanolVsAspirin_lessThan1` | 異なる分子 → score < 1.0 | 要 |
| `test_tanimoto_crossValidation_ethanolVsAspirin` | 手動 T = c/(a+b-c) と一致（AbsTol=1e-10）★ | 要 |
| `test_tanimoto_crossValidation_ethanolVsBenzene` | 同上、エタノール vs ベンゼン | 要 |
| `test_tanimoto_crossValidation_maccs_ethanolVsAspirin` | MACCS FP で手動計算と一致 | 要 |
| `test_tanimoto_morgan2048_vs_maccs167_throwsRdkitError` | Morgan(2048bit) vs MACCS(167bit) → `rdkitError` | 要 |
| `test_tanimoto_returnIsReal` | `isreal(score) == true` | 要 |

**★ クロス検証テスト (TC9) の重要性**:
手動計算 `sum(a & b) / sum(a | b)` と RDKit API の結果を AbsTol=1e-10 で比較する。
API 呼び出しミス・`toArray()` のビットエンコーディングエラー・`double()` 変換時の精度損失を検出する。

---

## 6.2 `dice`

**設計意図**: 2 つのフィンガープリント間の Dice 係数を計算する。ON ビット数が少ない分子対でより高い類似度を示す傾向があり、反応物・生成物間の比較や官能基ベースの類似度評価（MACCS keys）で補完的に使われる。

**数学的定義**:

$$
D(A, B) = \frac{2|A \cap B|}{|A| + |B|} = \frac{2c}{a + b}
$$

**Dice と Tanimoto の関係**:
バイナリベクトルに対して $D(A,B) \geq T(A,B)$ が常に成立する。

**証明**: $c \leq \min(a, b) \leq (a+b)/2$ より $a + b \geq 2c$。これを用いると:

$$
\frac{2}{a+b} \geq \frac{1}{a+b-c} \iff 2(a+b-c) \geq a+b \iff a+b \geq 2c \quad \checkmark
$$

**引用文献**:
- Sørensen, T. (1948). *Kongelige Danske Videnskabernes Selskab* 5(4), 1–34.

**テスト戦略** (`tests/unit/TestSimilarity.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_dice_fp1_numericInput_throwsInvalidInput` | double fp1 → `invalidInput` | 不要 |
| `test_dice_fp1_emptyInput_throwsInvalidInput` | 空 double fp1 → `invalidInput` | 不要 |
| `test_dice_returnsDoubleScalar` | 戻り値が double スカラー | 要 |
| `test_dice_sameFP_morgan_equals1` | 同一 Morgan FP → 1.0 | 要 |
| `test_dice_symmetry_ethanolVsAspirin` | D(A,B) = D(B,A) | 要 |
| `test_dice_crossValidation_ethanolVsAspirin` | 手動 D = 2c/(a+b) と一致（AbsTol=1e-10）★ | 要 |
| `test_dice_crossValidation_ethanolVsBenzene` | 脂肪族 vs 芳香族ペアで手動計算と一致 ★ | 要 |
| `test_dice_crossValidation_maccs_ethanolVsAspirin` | MACCS でも手動計算と一致 | 要 |
| `test_dice_greaterOrEqualTanimoto_ethanolVsAspirin` | D ≥ T（aspirin ペア）★ | 要 |
| `test_dice_greaterOrEqualTanimoto_ethanolVsBenzene` | D ≥ T（benzene ペア）★ | 要 |
| `test_dice_morgan2048_vs_maccs167_throwsRdkitError` | Morgan(2048bit) vs MACCS(167bit) → `rdkitError` | 要 |

**★ 一貫性テスト (TC9d) の重要性**:
`test_dice_greaterOrEqualTanimoto_*`: $D \geq T$ の数学的性質を利用して Dice と Tanimoto 双方の
結果が整合しているかを交差検証する。片方の実装に問題があればこのテストが検出する。

---

## 6.3 `rankBy`

**設計意図**: クエリフィンガープリントとデータベース全体の類似度を一括計算し、
上位 N 件を降順で返す最上位 API。仮想スクリーニング・Top-N 候補選出の中核関数。

**アルゴリズム概要**:
1. Metric バリデーション（`"tanimoto"` または `"dice"`。arguments block 内で即チェック）
2. `queryFp` の型バリデーション（`startsWith(class(fp), "py.")` — RDKit 不要）
3. `dbFps` の型バリデーション（cell 配列 + 全要素が Python オブジェクト — RDKit 不要）
4. `N` のバリデーション（正整数スカラーまたは Inf）
5. `py.list(dbFps)` で MATLAB cell → Python list に変換（IPC 1 往復に集約）
6. `py.rdkit.DataStructs.BulkTanimotoSimilarity(queryFp, pyList)` または `BulkDiceSimilarity` を呼び出し
7. `double(py.array.array("d", pyScores))` で MATLAB double 行ベクトルに変換
8. `sort(..., "descend")` で降順ソート → 上位 K = min(N, M) 件を抽出
9. `result.Indices`, `result.Scores`, `result.Metric` を struct にまとめて返す

**Bulk API の数値的根拠（ADR-002 IPC 最小化）**:

| 実装 | IPC 往復数 | N=100 での概算時間 |
|---|---|---|
| `tanimoto()` ループ | M 回 | 1–5 s |
| `BulkTanimotoSimilarity()` | 1 回 | 数十 ms〜数百 ms |

**`py.array.array("d", ...)` 変換の根拠**:
Python list of float を MATLAB double に直接変換するには `py.array.array("d", list)` を経由する必要がある。
型コード `"d"` は C double（64-bit IEEE 754 浮動小数点）を指定する。

**ソート結果のインデックス仕様**:
`sort(scores, "descend")` の第 2 出力は MATLAB 1-based インデックス（dbFps の元の位置）を返す。
呼び出し元は `dbFps{result.Indices(k)}` で k 番目の分子を参照できる。

**テスト戦略** (`tests/unit/TestSimilarity.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_rankBy_numericQueryFp_throwsInvalidQueryFp` | 数値 queryFp → `invalidQueryFp` | 不要 |
| `test_rankBy_emptyQueryFp_throwsInvalidQueryFp` | `[]` → `invalidQueryFp` | 不要 |
| `test_rankBy_invalidMetric_throwsInvalidMetric` | 未知 Metric → `invalidMetric` | 不要 |
| `test_rankBy_nZero_throwsInvalidN` | N=0 → `invalidN` | 要 |
| `test_rankBy_nNegative_throwsInvalidN` | N=-1 → `invalidN` | 要 |
| `test_rankBy_nNonInteger_throwsInvalidN` | N=1.5 → `invalidN` | 要 |
| `test_rankBy_returnsStruct` | 戻り値が struct | 要 |
| `test_rankBy_structHasRequiredFields` | Indices/Scores/Metric フィールド存在 | 要 |
| `test_rankBy_scoresAreSortedDescending` | `diff(Scores) <= 0` | 要 |
| `test_rankBy_selfInDb_scoreIsOneAtTop` | queryFp in dbFps → score=1.0 が 1 位 | 要 |
| `test_rankBy_n1_returnsSingleResult` | N=1 → Indices/Scores が [1,1] | 要 |
| `test_rankBy_nGreaterThanM_returnsAllM` | N > M → M 件返す | 要 |
| `test_rankBy_defaultN_returnsAllM` | N 省略 → M 件返す | 要 |
| `test_rankBy_metricDice_scoresInRange` | Metric="dice" → scores ∈ [0,1] | 要 |
| `test_rankBy_crossValidation_topScoreMatchesTanimoto` | topScore == tanimoto(fp, dbFps{topIdx})（AbsTol=1e-10）★ | 要 |
| `test_rankBy_nInfExplicit_returnsAllM` | N=Inf 明示 → M 件全返す | 要 |
| `test_rankBy_allIndices_correspondToCorrectScores` | 全 k に対し Scores(k) == tanimoto(query, dbFps{Indices(k)})（AbsTol=1e-10）★ | 要 |
| `test_rankBy_diceMetric_crossValidation_topScoreMatchesDice` | dice top score == pairwise dice(query, dbFps{topIdx})（AbsTol=1e-10）★ | 要 |

**★ クロス検証テスト (RB19) の重要性**:
`test_rankBy_crossValidation_topScoreMatchesTanimoto` は Bulk API の結果と単一ペア API の結果を
AbsTol=1e-10 で比較する。BulkTanimotoSimilarity と TanimotoSimilarity が一致しないバグや
`py.array.array` 変換でのインデックスずれを確実に検出する。

---

## 6.4 `matrix`

**目的**: N 個のフィンガープリントに対して N×N の全対型類似度行列を計算する。
ヒートマップ可視化、クラスタリング前処理、SAR 解析などで使用する。

**アルゴリズム概要**:
1. 入力バリデーション（fps が非空 cell、全要素が Python オブジェクト、Metric チェック）
2. `pyList = py.list(fps)` で全 FP を Python リストに変換
3. i=1..N のループで `BulkTanimotoSimilarity(fps{i}, pyList)` を 1 回呼び出し（N 回の IPC 往復）
4. `S = (S + S') / 2` で浮動小数点誤差による非対称性を除去
5. `S(1:N+1:end) = 1.0` で対角を厳密に 1.0 に固定

**計算量**: O(N²) の類似度評価、O(N) の IPC 往復。N > 1000 の場合は GPU Bulk 実装を推奨。

**対称性の強制**:
`(S + S') / 2` の操作は浮動小数点の往復誤差（± 1e-15 程度）による `S(i,j) ≠ S(j,i)` を防ぐ。
クラスタリングアルゴリズム（linkage 等）が対称性を仮定するため、強制しておく方が安全。

**テスト検証戦略** (`tests/unit/TestSimilarity.m`):

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_matrix_nonCellInput_throwsInvalidInput` | 数値 fps → `invalidInput` | 不要 |
| `test_matrix_emptyCellInput_throwsInvalidInput` | 空 cell → `invalidInput` | 不要 |
| `test_matrix_invalidMetric_throwsInvalidMetric` | 未知 Metric → `invalidMetric` | 不要 |
| `test_matrix_returnsDouble` | 戻り値が double | 要 |
| `test_matrix_isSquare` | size(S,1) == size(S,2) == N | 要 |
| `test_matrix_diagonalIsOne` | diag(S) がすべて 1.0 | 要 |
| `test_matrix_isSymmetric` | max(abs(S - S')) < 1e-12 | 要 |
| `test_matrix_valuesInRange` | S の全要素 ∈ [0, 1] | 要 |
| `test_matrix_crossValidation_vs_tanimoto` | S(1,2) == tanimoto(fps{1},fps{2})（AbsTol=1e-10）★ | 要 |
| `test_matrix_singleFP_returns1x1One` | N=1 → S = [1.0] | 要 |

**★ クロス検証テストの重要性**:
`test_matrix_crossValidation_vs_tanimoto` は行列の個別要素を `tanimoto()` で検証し、
Bulk API と単一ペア API の乖離バグを検出する。

---

## 6.5 GPU Matrix Tanimoto（R01 研究コンテンツ）

**設計意図**: 100k〜1M 規模のライブラリを GPU で高速スクリーニングするための行列演算実装。
`emk.similarity.rankBy` の BulkTanimotoSimilarity（Tier 2）と同等の結果を MATLAB/GPU で得る。

**アルゴリズム**:

$$
T(i) = \frac{\mathbf{f}_i \cdot \mathbf{q}}{|\mathbf{f}_i| + |\mathbf{q}| - \mathbf{f}_i \cdot \mathbf{q}}
$$

行列形式（$F \in \mathbb{R}^{N \times B}$ single、$\mathbf{q} \in \mathbb{R}^{1 \times B}$ single）:

```
intersection = F * q'               % N x 1  (GEMV)
onBitsA      = sum(F, 2)            % N x 1
onBitsQ      = sum(q)               % scalar
T            = intersection ./ (onBitsA + onBitsQ - intersection)
```

**GPU 実装の要点**:
- `single` 精度（float32）を使用する。GPU GEMM は float32 で最大スループットを達成する
- `gpuArray(F_single)` でホスト→デバイス転送（100k × 2048 × 4 B = 819 MB）
- `gather()` で最終スコアのみデバイス→ホスト転送（N × 8 B = 0.8 MB）
- バッチ戦略: `gpuDevice().AvailableMemory * 0.8 / (B * 4)` でバッチサイズを動的決定

**精度について**:
float32 の累積誤差により RDKit 整数ビットカウントとの最大差は $\approx 10^{-6}$（2048 ビットの場合）。
スクリーニング用途ではトレラブル（最終精度確認は RDKit Bulk API 推奨）。

**スループット（参考値）**:

| 手法 | 概算スループット | 100k mol でのコスト |
|---|---|---|
| 単一 IPC ループ（Tier 1） | ~500 mol/s | ~200 s |
| RDKit BulkTanimoto（Tier 2） | ~200,000 mol/s | ~0.5 s |
| GPU GEMV（Tier 3） | >5,000,000 mol/s | ~0.05 s |

**引用文献**:
- Cao Y et al. (2008). ChemFP. *J Chem Inf Model* 48:2208–2215. DOI: 10.1021/ci8001854
- MATLAB Parallel Computing Toolbox: `gpuArray`, `gpuDevice`, `gputimeit`

---

## 許容差一覧

| 関数 | AbsTol | 根拠 |
|---|---|---|
| `tanimoto` | 1e-10 | クロス検証用（手動計算 vs RDKit API） |
| `dice` | 1e-10 | 同上 |
| `matrix` crossval | 1e-10 | 同上 |
| GPU matrix vs RDKit | 1e-4 | float32 累積誤差（2048 bit GEMV） |
