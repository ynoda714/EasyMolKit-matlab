# Algorithm Guide — EasyMolKit

> 全関数のアルゴリズム根拠・数学的定義・テスト検証戦略。  
> 関数シグネチャ → [function_reference.md](function_reference.md)

**本ドキュメントの役割**: ソースコードの「なぜそう動くか」を記述する。
各モジュールの詳細は下表リンク先を参照のこと。

---

## 1. 処理パイプライン概観

```
config/settings.json ──→ loadConfig() ──→ params (struct)
                                              │
                              emk.setup.initPython()
                                  │ pyenv 設定・RDKit 検証
                                  ▼
SMILES / SDF ─→ emk.mol.fromSmiles() ──→ mol (py.rdkit obj)
                                              │
                              emk.descriptor.calculate(mol)
                              emk.fingerprint.morgan(mol)
                              emk.similarity.tanimoto(fp1, fp2)
                                              │
                              ──→ MATLAB table / struct / double
                                              │
                              makeRunDir() → result/runs/<ts>/
```

**設計原則**:
- 中間オブジェクト（Mol, Fingerprint）は Python 参照のまま保持（ADR-002）
- 最終出力段階で `emk.util.pyToMatlab()` を経由して MATLAB ネイティブ型に変換
- バッチ処理では Python 側でループを集約し IPC 往復を最小化（ADR-002 rev.3）
- 永続化は MolBlock 基本、Pickle は高忠実度オプション（ADR-004）

---

## モジュール別アルゴリズム詳細

| モジュール | ファイル | 主要関数 |
|---|---|---|
| emk.setup | [algorithms/setup.md](algorithms/setup.md) | install, installOnline, verify, initPython |
| emk.mol | [algorithms/mol.md](algorithms/mol.md) | fromSmiles, toSmiles, isValid, hasSubstruct, toStruct, fromStruct, toTable, scaffold |
| emk.descriptor | [algorithms/descriptor.md](algorithms/descriptor.md) | molWeight, calculate, batchCalculate, benchmarkBatch |
| emk.fingerprint | [algorithms/fingerprint.md](algorithms/fingerprint.md) | morgan, maccs, toArray |
| emk.similarity | [algorithms/similarity.md](algorithms/similarity.md) | tanimoto, dice, rankBy, matrix |
| emk.io | [algorithms/io.md](algorithms/io.md) | readSdf, writeSdf, readSmilesList |
| emk.viz | [algorithms/viz.md](algorithms/viz.md) | draw2d |
| emk.util / src/util / src/config | [algorithms/util.md](algorithms/util.md) | pyToMatlab, isOnline, loadConfig, logHelpers, makeRunDir |
| emk.filter | [algorithms/filter.md](algorithms/filter.md) | lipinski |
| emk.db | [algorithms/db.md](algorithms/db.md) | searchPubchem, searchChembl |

---

## 11. サンプルスクリプト: Analytics Layer 3

### A07 Scaffold 分析と R-group 分解 (`examples/analytics/a07_scaffold_analysis.m`)

**設計意図**: FDA 承認薬 200 件を対象に Bemis-Murcko 骨格ごとに分子をグループ化し、
R-group（側鎖）の物性バリエーションを可視化することで SAR（構造活性相関）の基礎を学ぶ。

**アルゴリズム概要**:
1. `emk.mol.fromSmiles` → `emk.mol.scaffold` → `emk.mol.toSmiles` で各薬物の骨格 SMILES を取得
2. 0 原子骨格（非環状分子）を `"<acyclic>"` タグに変換して除外せず追跡
3. `containers.Map` で骨格 SMILES → メンバー分子インデックスのマッピングを構築
4. 骨格頻度を降順ソートし、上位 N 骨格をバーチャートで表示
5. 3 件以上のメンバーを持つ骨格ファミリーごとに ALogP / TPSA / MW の mean・std を集計（R-group table）
6. CSV 既存物性値（MW, ALogP, TPSA 等）を使って記述子行列を組み立て、`pca()` で 2D 投影し骨格ファミリー別着色
7. `boxplot()` で骨格ファミリーごとの ALogP 分布を可視化
8. Scaffold Diversity Index (SDI) = 一意骨格数 / 総分子数 を算出

**数学的定義**:

$$\text{SDI} = \frac{|\{\text{unique scaffolds}\}|}{N}$$

**引用文献**:
- Bemis & Murcko (1996) *J. Med. Chem.* 39(15):2887-2893. doi:10.1021/jm9602928（骨格定義の原典）
- Bemis & Murcko (1999) *J. Med. Chem.* 42(25):5095-5099. doi:10.1021/jm9903996（R-group 分析）
- Langdon et al. (2011) *J. Chem. Inf. Model.* 51(9):2174-2185. doi:10.1021/ci200319g（FDA 薬物骨格多様性の定量解析）
- Ertl et al. (2006) *J. Med. Chem.* 49(15):4568-4573. doi:10.1021/jm060217p（環系多様性と特権骨格）

**テスト戦略**: サンプルスクリプトのためユニットテストは対象外。
動作検証は MATLAB Online で Section 単位実行（Ctrl+Enter）で確認する。

### A08 Mass Spectrometry × Cheminformatics

**設計意図**: 5 種の合成 ESI-MS スペクトルから未知化合物を同定するパイプラインを実装し、
正確質量測定・同位体パターン一致スコアリングの化学情報学的意義を段階的に学ぶ。

#### `simulateSpectrum`

1. `linspace(mzMin, mzMax, 5000)` で 5000 点の m/z 軸を生成
2. 各ピーク位置にガウス関数 $\exp(-\frac{1}{2}(\frac{m/z - m/z_0}{\sigma})^2)$ を重畳
3. ガウス成分を `[0, 1]` に正規化後、一様乱数ノイズ (`noiseLevel × rand`) を加算
4. ノイズ加算後の最大強度は `1 + noiseLevel` 未満（sigma = 0.3 Da → FWHM ≈ 0.71 Da で低分解能 ESI-MS を模擬）

#### `parseFormula`

- 正規表現 `([A-Z][a-z]?)(\d*)` で Hill 記法の分子式をトークン化
- `[A-Z][a-z]?` は 1 文字大文字 + 任意 1 文字小文字 → `Cl`, `Br` 等の 2 文字元素を正しく捕捉
- 空の桁文字列（例: 末尾 `O`）は `"1"` に補完

**検証済みエッジケース:**

| 入力 | 期待出力 | 確認事項 |
|---|---|---|
| `"CH4"` | C:1, H:4 | 単純な 1 桁カウント |
| `"C9H8O4"` | C:9, H:8, O:4 | 複数桁カウント（アスピリン） |
| `"C6H5Cl"` | C:6, H:5, Cl:1 | 2 文字元素（C と l ではなく Cl として解析） |
| `"CCl4"` | C:1, Cl:4 | 四塩化炭素 |
| `""` | 全フィールド 0 | 空文字列は安全にスキップ |

#### `isoPattern`

**近似式** (Gross 2011, Chapter 3; C 原子数 < 100 で有効):

$$M+1\% = 1.103 n_C + 0.366 n_N + 0.015 n_H + 0.038 n_O$$

$$M+2\% = \frac{(1.103 n_C)^2}{200} + 0.205 n_O + 4.25 n_S + 32.7 n_{Cl}$$

**係数の導出根拠** (自然存在比: Meija et al. 2016 IUPAC):

| 元素 | 自然存在比 | 使用する係数 | 導出 |
|---|---|---|---|
| 13C | 1.103% | 1.103% per C (M+1) | 直接適用（12C 基準） |
| 34S | 4.25% | 4.25% per S (M+2) | 直接適用（32S 基準） |
| 37Cl | **24.23%** | **32.7% per Cl (M+2)** | 37Cl/35Cl 比 = 24.23/75.77 = **31.98%** を Gross (2011) 表に従い 32.7% に近似 |

> **24.23% と 32.7% の違い**: 24.23% は 37Cl の天然存在比（IUPAC）。
> 32.7%（≈ 32.0% 厳密値）は M+2 ピーク強度の M ピーク基準の比（= 37Cl 確率 ÷ 35Cl 確率）。
> 単純比では 31.98%; Gross (2011) は 32.7% に丸めている（差 < 2.2%）。

**交差検証: 既知化合物での確認値** (Claesen et al. 2012 と整合):

| 化合物 | 分子式 | M+1% (計算値) | M+2% (計算値) | 検証根拠 |
|---|---|---|---|---|
| Methane | CH4 | 1.163% | 0.006% | 13C 1 原子分のみ; M+2 はほぼ検出不可 |
| Aspirin | C9H8O4 | 10.199% | 1.313% | 9C の寄与が支配的; Claesen et al. 誤差 < 0.5% |
| Chlorobenzene | C6H5Cl | 6.693% | 32.919% | M:M+2 ≈ 3:1 の「Cl 規則」(Gross 2011 p.95) と一致 |

**ppm 検索許容差の根拠**:
- `LOW_RES_TOL = 0.50 Da`: 単一四重極・イオントラップの典型仕様
- `HIGH_RES_PPM = 5.0 ppm`: Orbitrap/Q-TOF の実運用精度（合成スペクトルでは誤差 0 ppm; Kind & Fiehn 2006 は 1 ppm でも曖昧性が残ると報告）

**コサイン類似度スコアリング** (Stein & Scott 1994):

$$\text{IsoScore} = \frac{\vec{I}_{\text{obs}} \cdot \vec{I}_{\text{theo}}}{\|\vec{I}_{\text{obs}}\| \cdot \|\vec{I}_{\text{theo}}\|}$$

$\vec{I} = [I_M,\, I_{M+1},\, I_{M+2}]$ を M ピーク基準で正規化。スコア 1.0 = 完全一致、0.0 = 直交（全く異なる同位体プロファイル）。

**引用文献**:
- Gross JH (2011) *Mass Spectrometry*, 2nd ed. Springer. ISBN 978-3-642-10709-2
- Kind T & Fiehn O (2006) *BMC Bioinformatics* 7:234. doi:10.1186/1471-2105-7-234
- Claesen J et al. (2012) *J Am Soc Mass Spectrom* 23:753-763. doi:10.1007/s13361-011-0326-2
- Stein SE & Scott DR (1994) *J Am Soc Mass Spectrom* 5:859-866. doi:10.1016/1044-0305(94)87009-8
- Meija J et al. (2016) *Pure Appl. Chem.* 88:265-291. doi:10.1515/pac-2015-0305

**テスト戦略**: サンプルスクリプトのためユニットテストは対象外。
ローカル関数の科学的正確性は上表の cross-validation 値で担保する。
動作検証は MATLAB Online で Section 単位実行（Ctrl+Enter）で確認する。

### R03 Forensic Chemometrics

**設計意図**: 55 件の法科学参照データベースを用いたケモメトリクスパイプラインを実装し、
PCA + LDA による多クラス化学分類と Parallel Computing Toolbox (parfor) の適用限界を学ぶ。

**アルゴリズム概要**:
1. ECFP4 フィンガープリント（2048 ビット）行列を構築（55×2048）
2. PCA で 15 主成分に圧縮（`pca()`, ColRankDefX 警告抑制）
3. LDA (`fitcdiscr`, linear) を 8 カテゴリラベルで学習
4. LOOCV (`crossval` + `kfoldPredict`) で分類精度を推定
5. 200 件 FDA 薬物ライブラリを**参照 DB 平均で中心化**してから参照 PCA 空間に射影し LDA で予測
6. 未知 5 件を LDA 事後確率 + Tanimoto 最近傍の 2 源証拠で分類

**LDA 数学定義** (Fisher 1936):

$$J(\mathbf{w}) = \frac{\mathbf{w}^\top S_B \mathbf{w}}{\mathbf{w}^\top S_W \mathbf{w}}$$

$S_W$: クラス内散布行列、$S_B$: クラス間散布行列。
一般化固有値問題 $S_B \mathbf{w} = \lambda S_W \mathbf{w}$ の上位 $C-1$ 固有ベクトルが判別軸。
`fitcdiscr` (linear) は $S_W$ が特異に近い場合に疑似逆行列にフォールバックする。

**パラメータ選択根拠**:

| パラメータ | 値 | 根拠 |
|---|---|---|
| k（PCA 主成分数） | 15 | n=55 に対し p=2048 は `p >> n`。LDA の前処理として `k < n` が必須条件。k=15 でフィンガープリント分散の ~70–80% を説明（スクリープロットで確認）。文献的範囲: フィンガープリント PCA では k=10〜50 が標準 (Varmuza & Filzmoser 2009 §5.3) |
| T（Tanimoto 閾値） | 0.85 | Willett et al. (1998) が T ≥ 0.85 を「非常に類似した骨格」の実用的閾値として提唱した化学情報学の慣習。法的標準ではなく教育目的の指標として使用 |
| LDA prior | 頻度比例（default） | `fitcdiscr` デフォルト。クラス不均衡（drug: ~50%）を事前確率に反映。均一事前確率は `Prior="uniform"` で変更可（TRY IT 3） |

**PCA 中心化の設計意図**:
FDA ライブラリを参照 DB の PCA 空間に射影する際は、**参照 DB の列平均 `bitMean`** で中心化する
（FDA ライブラリ自身の平均ではない）。訓練時と同一の座標系を保つために必須。

**LOOCV 精度の期待範囲**: n=55・8 カテゴリでは LOOCV 精度 55〜80% が現実的。
単一誤分類が精度値を ~2% 変化させる高分散推定であることに注意（結果は「目安」として解釈）。

**Figure 4 設計注記（科学的正確性）**:
`predict(fitcdiscr)` の第 2 戻り値は **事後確率行列** $P(\text{class}|x)$ であり、
LDA 判別座標（canonical variate scores）ではない。
Figure 4 は P(drug | ECFP4) を真カテゴリ別ジッター散布図で可視化することで、
モデルの「薬物様度」識別能力を示す。
LDA 判別空間を直接可視化するには $S_W^{-1} S_B$ の固有ベクトルを手動計算する必要があり、
高度な拡張課題として TRY IT コメントで言及する。

**parfor の適用限界**:
RDKit フィンガープリント計算は Python IPC（アウトプロセス）がボトルネックで parfor で高速化不可。
`toArray()` 変換（純粋 MATLAB）は parfor による並列化で高速化できる。

**引用文献**:
- Fisher RA (1936) *Ann Eugenics* 7:179-188. doi:10.1111/j.1469-1809.1936.tb02137.x（LDA 原典）
- Willett P, Barnard JM & Downs GM (1998) Chemical similarity searching. *J Chem Inf Comput Sci* 38:983-996. doi:10.1021/ci9800211（T ≥ 0.85 閾値慣習）
- Varmuza K & Filzmoser P (2009) *Introduction to Multivariate Statistical Analysis in Chemometrics*. CRC Press.（PCA 主成分数選択・LDA 前処理）
- Rogers D & Hahn M (2010) *J Chem Inf Model* 50:742-754. doi:10.1021/ci100050t（ECFP4）

**テスト戦略**: サンプルスクリプトのためユニットテストは対象外。
LOOCV 精度（55〜80% 範囲）を動作検証の主指標とし、MATLAB Desktop で Section 単位実行で確認する。

---

### R05 分子言語モデル — SMILES 生成

**設計意図**: SMILES 文字列をシーケンスとして扱い、(1) バイグラム Markov Chain で「ゼロ訓練・局所精度」を実演し、(2) 文字レベル LSTM でその限界（長距離依存）を克服する構成で、分子生成 AI の根幹概念を実装ベースで体得する。

#### バイグラム Markov Chain（Sections 3〜5）

**遷移テーブル構築**:

$$P(c_t \mid c_{t-1}) = \frac{\text{count}(c_{t-1}, c_t) + \alpha}{\sum_{c'} \left[\text{count}(c_{t-1}, c') + \alpha\right]}, \quad \alpha = 0.1\ \text{(Laplace平滑化)}$$

訓練データ（200 分子 ≈ 7,000 バイグラム）から遷移確率行列 `transProb [V×V]` を構築。訓練時間 ≈ 0 秒。

**期待性能** (200 mol): Validity 20〜35%。局所遷移（'(' の後は原子）は正確だが括弧対応・環閉鎖は記憶不足で失敗。

**限界（Section 5）**:
- **括弧マッチング**: ブランチ開閉は数文字離れ → 1 ステップ記憶では追跡不能
- **環クロージャ**: 数字 '1' の開閉は 6〜10 文字離れ → Markov の射程外
- **原子価整合**: 現炭素の結合数は複数文字にまたがる

#### キャラクタ言語モデル（Sections 6〜8）

**語彙化**: SMILES 全文字の sorted-unique + 3 特殊トークン（START=`^`, END=`$`, PAD=`_`）。語彙サイズ: 通常 40〜50 文字。

**教師強制 (Teacher Forcing)**:

$$L = -\sum_{t=1}^{T} \log P_\theta(c_{t+1} \mid c_1, \ldots, c_t)$$

PAD トークンはマスクして損失から除外（`mask = Yt ~= PAD_IDX`）。

**温度サンプリング**:

$$P_T(c) = \text{softmax}\!\left(\frac{\log P(c)}{T}\right), \quad T < 1 \Rightarrow \text{保守的}, \; T > 1 \Rightarrow \text{探索的}$$

**生成方式**: 増加プレフィックスを毎ステップ渡す ($O(T^2)$ だが実装が単純で教育的)。本番では `HasStateInputs=true` の LSTM で $O(T)$ にする。

**スケーリング則**: 200〜500 mol では損失がユニグラムエントロピー（≈2.6）付近に収束し、長距離文法学習に必要なデータ（≈10k+）が不足。Markov > LSTM の逆転現象が示す帰納的バイアスの重要性。

#### REINFORCE（R06: `r06_reinforce.m` Section 2）

**政策勾配定理**:

$$\nabla_\theta J = \mathbb{E}_{\tau \sim \pi_\theta}\!\left[R(\tau) \sum_{t=1}^{T} \nabla_\theta \log \pi_\theta(a_t \mid s_t)\right]$$

バッチ平均を baseline として分散低減:

$$\hat{g} = \frac{1}{N} \sum_{i=1}^{N} (R_i - b) \sum_t \nabla_\theta \log \pi_\theta(a_t^i \mid s_t^i)$$

**報酬関数**: Lipinski Ro5 複合スコア（MW ≤ 500 / HBD ≤ 5 / HBA ≤ 10 / LogP ≤ 5 で各 +0.25、新規性ボーナス ×1.1）。最大値 1.0。

**Catastrophic forgetting 対策** (TRY IT 10): KL 正則化項
$\beta \cdot \mathrm{KL}(\pi_\theta \| \pi_{\text{prior}})$
を損失に加算して事前言語モデルからの逸脱を抑制（REINVENT 拡張尤度: Blaschke et al. 2020）。

**評価指標** (Brown et al. 2019 / Polykovskiy et al. 2020):
- **Validity**: RDKit が解析できる SMILES の割合
- **Uniqueness**: 有効 SMILES のうち重複しない割合
- **Novelty**: 訓練セット非包含の有効 SMILES の割合

**引用文献**:
- Gomez-Bombarelli et al. (2018) *ACS Cent Sci* 4:268-276. doi:10.1021/acscentsci.7b00572（VAE-SMILES 生成モデル）
- Segler et al. (2018) *ACS Cent Sci* 4:120-131. doi:10.1021/acscentsci.7b00512（LSTM SMILES 生成・転移学習）
- Olivecrona et al. (2017) *J Cheminform* 9:48. doi:10.1186/s13321-017-0235-x（REINVENT: 分子生成に REINFORCE を適用）
- Williams RJ (1992) *Mach Learn* 8:229-256. doi:10.1007/BF00992696（REINFORCE アルゴリズム原典）
- Brown et al. (2019) *J Chem Inf Model* 59:1096-1108. doi:10.1021/acs.jcim.8b00839（GuacaMol 評価指標）
- Polykovskiy et al. (2020) *Front Pharmacol* 11:565644. doi:10.3389/fphar.2020.565644（MOSES ベンチマーク）
- Hochreiter & Schmidhuber (1997) *Neural Comput* 9:1735-1780. doi:10.1162/neco.1997.9.8.1735（LSTM 原典）

**テスト戦略**: サンプルスクリプトのためユニットテストは対象外。
動作検証は Markov Validity ≥ 15%（200 mol 訓練時の最低期待値）で確認。
REINFORCE 評価指標（Pre/Post-RL 平均報酬比較）は R06 (`r06_reinforce.m`) を参照。

Section 5 の期待値（事後検証に使用）:

| 未知サンプル | 真化合物 | 期待 NN Tanimoto | 期待 LDA カテゴリ |
|---|---|---|---|
| U001 | caffeine | T=1.0（DB 内: FC054） | stimulant（DB ラベルが is_drug=0） |
| U002 | ibuprofen | T=1.0（DB 内: FC050） | analgesic |
| U003 | nicotine | T=1.0（DB 内: FC042） | stimulant |
| U004 | acetaminophen | T=1.0（DB 内: FC034） | analgesic |
| U005 | imipramine | T=1.0（DB 内: FC013） | drug |

> U001 の LDA カテゴリは「stimulant」または「drug」のいずれかになりうる（カフェインはラベルノイズ）。
> T=1.0 が得られない場合は SMILES の正規化差異が原因（RDKit は入力 SMILES を正規化する）。

### R04 Protein-Ligand Analysis

**設計意図**: FDA 承認薬 200 件を 5 種の承認済みキナーゼ阻害剤（イマチニブ・エルロチニブ・ゲフィチニブ・ダサチニブ・ソラフェニブ）を参照とするリガンドベース仮想スクリーニング（LBVS）でランク付けし、Bioinformatics Toolbox + Biopython による配列解析を組み合わせたタンパク質-リガンド解析の教育的パイプラインを提供する。3D ドッキングは使用しない。

**アルゴリズム概要**:
1. 参照キナーゼ阻害剤 5 件の ECFP4 フィンガープリントを計算
2. FDA ライブラリ全件のフィンガープリントを計算し、参照 5 件に対する最大 Tanimoto 類似度 (`maxSim`) を `emk.similarity.rankBy` で計算
3. 複合スコアリング（0.4 × 類似度 + 0.4 × 物性マッチ + 0.2 × LE プロキシ）で候補をランク付け
4. Biopython で PKA 触媒サブユニット（UniProt P17612）の FASTA を解析し、DFG モチーフ・残基組成を注釈
5. Bioinformatics Toolbox (`swalign`) で P-loop モチーフのローカルアライメントを実施（不在時は文字列マッチにフォールバック）
6. 上位 10 件の 2D 構造図を `emk.viz.draw2d` で描画し、CSV + JSON レポートを出力

**複合スコア数学定義**:

$$\text{Score} = 0.4 \cdot S_{\text{sim}} + 0.4 \cdot S_{\text{prop}} + 0.2 \cdot S_{\text{LE}}$$

**物性マッチスコア**（ガウシアンペナルティ）:

$$S_{\text{prop}} = \frac{1}{3} \left[ e^{-\frac{(\text{MW} - \mu_{\text{MW}})^2}{2\sigma_{\text{MW}}^2}} + e^{-\frac{(\text{LogP} - \mu_{\text{LogP}})^2}{2\sigma_{\text{LogP}}^2}} + e^{-\frac{(\text{TPSA} - \mu_{\text{TPSA}})^2}{2\sigma_{\text{TPSA}}^2}} \right]$$

$\mu$ はキナーゼ阻害剤参照 5 件の中央値、$\sigma$ は薬物化学的許容幅（$\sigma_{\text{MW}}=80$, $\sigma_{\text{LogP}}=1.5$, $\sigma_{\text{TPSA}}=40$）。

**LE プロキシ**（Hopkins et al. 2004 をベースに教育用に単純化）:

$$S_{\text{LE}} = \min\left(1,\ \frac{\exp(-|\text{LogP} - 3.0| / 2)}{\text{MW} / 300}\right)$$

理想値 (MW=300, LogP=3.0) で $S_{\text{LE}} = 1.0$。大型分子は $\text{MW}/300 > 1$ で自動的にペナルティを受ける。
キナーゼ阻害剤は bRo5 領域（MW > 500）が多いため、このコンポーネントは意図的に低ウェイト（0.2）に設定している。

**パラメータ選択根拠**:

| パラメータ | 値 | 根拠 |
|---|---|---|
| W_SIM = 0.4 | 類似度ウェイト | 活性クラス予測能が最も高いが、単独では ADMET を無視するため 0.4 に抑制 |
| W_PROP = 0.4 | 物性ウェイト | 類似度と相補的。経口吸収性・ADMET リスクをカバー |
| W_LE = 0.2 | LE プロキシウェイト | 二次的補助指標。高 MW に軽いペナルティを与える |
| maxSim 閾値 0.40 | 「同一スキャフォールド」の実用的下限 | Willett et al. (1998) が T ≥ 0.85 を厳密な閾値とし、T ≥ 0.4 は化学的近傍の目安として使用 |

**引用文献**:
- Bleicher et al. (2003) *Nat Rev Drug Discov* 2:369-378. doi:10.1038/nrd1086（LBVS と複合スコアリング）
- Hopkins AL & Groom CR (2002) *Nat Rev Drug Discov* 1:727-730. doi:10.1038/nrd892（LE の概念）
- Willett P et al. (1998) *J Chem Inf Comput Sci* 38:983-996. doi:10.1021/ci9800211（Tanimoto 類似度閾値の実用的解釈）
- Cock PJA et al. (2009) *Bioinformatics* 25:1422-1423. doi:10.1093/bioinformatics/btp163（Biopython）
- Knighton DR et al. (1991) *Science* 253:407-414. doi:10.1126/science.1862342（1ATP PKA 構造; DFG モチーフの文脈）

**テスト戦略**: サンプルスクリプトのためユニットテストは対象外。
動作検証は Section 単位実行（Ctrl+Enter）で確認する。
Biopython 不在時の fallback（文字列マッチ）と Bioinformatics Toolbox 不在時の fallback（手動モチーフスキャン）が両方とも正常完了することを確認する。

---

### R06 REINFORCE 分子設計

**設計意図**: R05 で訓練した LSTM 言語モデルを REINFORCE 政策勾配法で微調整し、創薬標的（Lipinski Ro5 合致）分子の生成確率を上げる。REINVENT（Olivecrona et al. 2017）の教育実装。

**政策勾配定理**（Williams 1992）:

$$\nabla_\theta J = \mathbb{E}_{\tau \sim \pi_\theta}\!\left[(R(\tau) - b) \sum_{t=1}^{T} \nabla_\theta \log \pi_\theta(a_t \mid s_t)\right]$$

$b = \bar{R}$（バッチ平均報酬）の baseline 減算で勾配分散を低減。

**MDP 定式化**:
- 状態 $s_t$: LSTM 隠れ状態 $h_t$（部分 SMILES たどり）
- 行動 $a_t$: 次文字トークン（離散, $|\mathcal{A}| = V_\text{vocab}$）
- 報酬 $R$: 完全な SMILES の Lipinski Ro5 複合スコア（終端報酬）

**報酬関数**（Lipinski Ro5 複合スコア, $R \in [0, 1]$）:
- MW $\leq 500$: +0.25
- HBD $\leq 5$: +0.25
- HBA $\leq 10$: +0.25
- LogP $\leq 5$: +0.25
- 訓練セット非包含（新規性）: ×1.1（上限 1.0）

**壊滅的忘却 (Catastrophic Forgetting)**:
RL 更新により SMILES 文法が崩壊するリスク。Validation loss 監視で検知。
KL 正則化（REINVENT 拡張層利益, TRY IT 2）:

$$\mathcal{L}_\text{total} = \mathcal{L}_{PG} + \beta \cdot \mathrm{KL}(\pi_{\theta} \| \pi_{\text{prior}})$$

**評価指標** (Brown et al. 2019 / MOSES Polykovskiy et al. 2020):
- **Validity**: RDKit 解析可能な SMILES の割合
- **Uniqueness**: 有効 SMILES のうち重複しない割合
- **Novelty**: 訓練セット非包含の有効 SMILES の割合
- **Avg reward**: バッチ平均 Lipinski 報酬

**引用文献**:
- Olivecrona M et al. (2017) *J Cheminform* 9:48. doi:10.1186/s13321-017-0235-x（REINVENT 原典）
- Williams RJ (1992) *Mach Learn* 8:229-256. doi:10.1007/BF00992696（REINFORCE 原典）
- Brown N et al. (2019) *J Chem Inf Model* 59:1096-1108. doi:10.1021/acs.jcim.8b00839（GuacaMol 評価指標）
- Blaschke T et al. (2020) *J Chem Inf Model* 60:5918-5922. doi:10.1021/acs.jcim.0c00915（KL 正則化 / 拡張層利益）

**テスト戦略**: サンプルスクリプトのためユニットテストは対象外。
動作検証: Pre-RL 有効性 ≥ 5%、Post-RL 平均報酬 > Pre-RL 平均報酬、Validation loss の暗費な坦上がないこと。

---

### R09 GNN 分子性質予測

**設計意図**: 分子グラフ（原子=ノード、結合=エッジ）を入力とする Graph Convolutional Network (GCN) で ALogP を予測し、手作り記述子ベース手法（A03/A05）との比較を通じて GNN が有利な条件・不利な条件を学ぶ。Python GNN パイプライン (`gnn_property.py`) を MATLAB の `system()` 経由で呼び出す実装で、PyG フレームワークを MATLAB ワークフローに組み込むパターンを実証する。

**GCN 更新則** (Kipf & Welling 2017):

$$H^{(l+1)} = \sigma\!\left(\tilde{D}^{-1/2} \tilde{A} \tilde{D}^{-1/2} H^{(l)} W^{(l)}\right)$$

- $\tilde{A} = A + I$（自己ループ付き隣接行列）
- $\tilde{D} = \mathrm{diag}(\tilde{A}\mathbf{1})$（次数行列）
- $H^{(l)} \in \mathbb{R}^{N \times d}$（第 $l$ 層のノード特徴行列）
- $W^{(l)}$（学習可能重み行列）
- $\sigma$（ReLU 活性化関数）

**グラフレベル集約** (global mean pooling):

$$z = \frac{1}{|V|} \sum_{v \in V} h_v^{(L)}, \quad \hat{y} = W_\text{out} z + b$$

**原子特徴ベクトル** (`gnn_property.py::atom_features`, 33 次元):

| 次元 | 特徴 | 根拠 |
|---|---|---|
| 0–9 | 元素 one-hot (C,N,O,S,F,Cl,Br,I,P,他) | 電子配置・電気陰性度を間接的に符号化 |
| 10–16 | 次数 one-hot (0–5, 他) | 局所結合パターン（分岐点 vs 鎖状）|
| 17–22 | 形式電荷 one-hot (−2〜+2, 他) | 荷電状態 → 溶解度・LogP に直接影響 |
| 23–29 | 混成 one-hot (S,SP,SP2,SP3,SP3D,SP3D2,他) | π系（芳香族・カルボニル）を捕捉 |
| 30 | Is aromatic (0/1) | 疎水性 π スタッキング |
| 31 | Is in ring (0/1) | 骨格 vs 側鎖の区別 |
| 32 | 暗黙的 H 数 / 4（正規化）| LogP の水素結合ドナー能と相関 |

**モデルアーキテクチャ**:
```
GCNConv(33→64) → BatchNorm → ReLU
GCNConv(64→64) → BatchNorm → ReLU
GCNConv(64→64) → BatchNorm → ReLU
global_mean_pool
Linear(64→1)                         → 予測 ALogP
```
学習可能パラメータ数: ~17,000（比較: A05 feedforward NN ~10,000; A03 線形回帰 9 パラメータ）

**Python–MATLAB 連携設計**:
- MATLAB → `system("python gnn_property.py --input smiles.csv --output results.json")` → Python GNN 訓練 → JSON 結果読み込み
- `pyenv().Executable` でセッション中の Python 実行ファイルパスを取得（private 関数不要）
- JSON シリアライズにより PyG テンソルの py.* IPC 変換問題を回避（ADR-002 rev.3 の IPC 最小化原則の適用）

**データセット・学習設定**:

| 項目 | 値 | 根拠 |
|---|---|---|
| データセット | FDA drugs (ChEMBL, ~175 valid) | A03/A05 と同一。公平な比較 |
| 分割 | 70 / 15 / 15 (train/val/test) | A05 と同一の乱数シード 42 |
| エポック | 150 | Sec 3 のデフォルト。TRY IT で変更可 |
| バッチサイズ | 32 | GPU なし CPU 学習での安定性 |
| オプティマイザ | Adam (lr=1e-3, wd=1e-4) | 分子 GNN のデファクト標準 |
| LR スケジューラ | StepLR(step=50, γ=0.5) | 後半の学習率低下で収束を安定化 |

**期待性能** (CPU, R2025b, 約 3 分):
ALogP 予測 test RMSE ~0.65–0.85 程度（A03 線形回帰 ~0.77, A05 NN ~0.65 と同程度）。
小データセットでは GNN は descriptor NN に対して有意な優位を持ちにくい（Hu et al. 2020 の知見）。
GNN が真価を発揮するのは ~10k 分子以上のスケール（TRY IT 4 で参照文献を提示）。

**引用文献**:
- Kipf TN & Welling M (2017) Semi-supervised classification with GCNs. ICLR 2017. arXiv:1609.02907（GCN 層）
- Gilmer J et al. (2017) Neural message passing for quantum chemistry. ICML 2017. arXiv:1704.01212（MPNN 統一フレームワーク）
- Duvenaud D et al. (2015) Convolutional networks on graphs for learning molecular fingerprints. NeurIPS 2015. arXiv:1509.09292（分子グラフ畳み込みの原典）
- Yang K et al. (2019) Analyzing learned molecular representations (Chemprop). J Chem Inf Model 59:3370-3388. doi:10.1021/acs.jcim.9b00237（D-MPNN ベンチマーク）
- Hu W et al. (2020) Open Graph Benchmark. NeurIPS 2020. arXiv:2005.00687（GNN ベンチマーク: 大規模データでの優位性）
- Fey M & Lenssen JE (2019) Fast graph representation learning with PyG. arXiv:1903.02428（PyTorch Geometric）

---

### R10 ChemBERTa 転移学習

**設計意図**: HuggingFace Transformers の ChemBERTa（PubChem/ZINC で事前学習済み SMILES BERT）を ALogP 予測に転移学習し、「スクラッチ学習（R09 GCN）vs 転移学習（R10 ChemBERTa）」の比較を通じて事前学習の教育的価値を伝える。Python ヘルパー `chemberta.py` を MATLAB の `system()` 経由で呼び出す設計（`gnn_property.py` と同パターン）。

**モデルアーキテクチャ**:

ChemBERTa は RoBERTa 系 Transformer エンコーダ。SMILES 文字列を BPE トークン列として処理する。

```
SMILES文字列
  → BPE トークナイザ（SMILES専用語彙）
  → [CLS] t1 t2 ... tn [SEP]  (max_length=128)
  → TransformerEncoder (12 layer, 12 head, hidden=768)  ← RoBERTa-base 準拠
  → [CLS] トークンの最終埋め込み (768 次元)
  → Dropout(p=0.1)
  → Linear(768 → 1)                   → 予測 ALogP
```

事前学習パラメータ数: ~86M（`seyonec/ChemBERTa-zinc-base-v1`, RoBERTa-base ベース, hidden=768）

**転移学習戦略比較**:

| 戦略 | 学習可能パラメータ | 適合データ量 | リスク |
|---|---|---|---|
| 凍結エンコーダ（線形プロービング）| ~769（ヘッドのみ）| < 100 分子 | 事前学習埋め込みに ALogP 情報が不足の場合アンダーフィット |
| 全層ファインチューニング | ~86M | > 500 分子 | 小データでは壊滅的忘却リスク。LR 2e-5 + weight decay で緩和 |

**損失関数・学習設定**:

$$\mathcal{L} = \frac{1}{N}\sum_{i=1}^N (\hat{y}_i - y_i)^2 \quad \text{(MSE)}, \quad \text{optimizer: Adam}$$

| 項目 | デフォルト値 | 根拠 |
|---|---|---|
| LR | 2e-5 | Transformer ファインチューニングの標準値（大 LR は壊滅的忘却） |
| バッチサイズ | 16 | CPU メモリ ~4 GB 以内 |
| エポック | 30 | 小データ収束の目安。TRY IT で 60 に増やして確認可 |
| LR スケジューラ | StepLR(step=epochs//3, γ=0.5) | 後半の LR 低下で収束安定化 |
| max_length | 128 | FDA drug SMILES の最大長カバー。計算コストとのトレードオフ |

**Python–MATLAB 連携設計**:
- `chemberta.py` は `--mode embed`（埋め込み抽出）と `--mode finetune`（ファインチューニング）の 2 モードを持つ CLI
- SMILES CSV → `system("python chemberta.py --mode finetune ...")` → JSON（学習曲線 + predictions）→ MATLAB で可視化
- HuggingFace Hub からのモデルダウンロード（初回のみ）はキャッシュに保存（`~/.cache/huggingface/`）
- JSON 経由で PyTorch テンソルの py.* 型変換を回避（ADR-002 rev.3 の IPC 最小化原則）

**期待性能** (CPU, 初回モデルダウンロード込み ~5 分):
小データセット（~175 分子）では全層ファインチューニングは GCN（R09）と同程度か若干劣ることが多い。
凍結エンコーダは正則化された設定として機能し、RMSE ~0.80–1.0 程度を期待。
大規模データ（数千分子以上）では ChemBERTa の転移学習が GCN を上回る（Ahmad et al. 2022 参照）。

**引用文献**:
- Chithrananda S, Grand G & Ramsundar B (2020) ChemBERTa: large-scale self-supervised pretraining for molecular property prediction. arXiv:2010.09885（ChemBERTa 原著）
- Ahmad W et al. (2022) ChemBERTa-2: towards chemical foundation models. arXiv:2209.01712（ChemBERTa-2 拡張）
- Devlin J et al. (2019) BERT: pre-training of deep bidirectional transformers for language understanding. NAACL 2019. arXiv:1810.04805（BERT アーキテクチャ）
- Liu Y et al. (2019) RoBERTa: a robustly optimised BERT pretraining approach. arXiv:1907.11692（ChemBERTa が採用した RoBERTa 改良）
- Hu W et al. (2020) Strategies for pre-training graph neural networks. ICLR 2020. arXiv:1905.12265（GNN 事前学習との比較軸）

**テスト戦略**: サンプルスクリプトのためユニットテストは対象外。
動作検証: test RMSE < 1.5（ベースライン: 全分子の平均値で予測した場合の RMSE ≈ std(ALogP) ≈ 1.3 よりも改善できているか）。JSON 出力の必須キー存在確認（`num_train`, `test_rmse`, `predictions`）。

---

## 12. 付録: 許容差一覧

| 関数 | AbsTol | RelTol | 根拠 |
|---|---|---|---|
| `molWeight` | 0.01 | — | PubChem 公称値の有効数字に基づく |
| `tanimoto` | 0.001 | — | 浮動小数点精度（RDKit C++ → Python → MATLAB 変換経路） |
| `dice` | 0.001 | — | 同上 |

### 合成データ生成テンプレート

```matlab
% Reference molecules
smiles = ["CCO", "c1ccccc1", "CC(=O)OC1=CC=CC=C1C(=O)O"];
names  = ["Ethanol", "Benzene", "Aspirin"];

% Expected molecular weights (PubChem)
expectedMW = [46.069, 78.112, 180.159];
```
