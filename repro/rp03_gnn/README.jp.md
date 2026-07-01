# RP03: Graph Neural Network - GCN on BBBP Classification

[English README](./README.md)

> **目的**: Yang et al. (2019) の Chemprop 系分子予測を参考に、BBBP の血液脳関門透過性分類を
> Graph Convolutional Network (GCN) で再現する。あわせて RP02 の ECFP4+LR ベースライン
> と比較し、固定長フィンガープリントに対してグラフ学習が優位かを確認する。
> さらに MATLAB Deep Learning Toolbox (DLT) 実装を Python/PyG 実装と比較する。

---

## Overview

| 項目 | 内容 |
|---|---|
| 論文 | Yang, K. et al. (2019). Analyzing Learned Molecular Representations for Property Prediction. *J. Chem. Inf. Model.* 59(8):3370-3388. |
| DOI | [10.1021/acs.jcim.9b00237](https://doi.org/10.1021/acs.jcim.9b00237) |
| タスク | 二値分類: BBB 透過性（BBB+ = permeable, BBB- = non-permeable） |
| モデル | Python: 3-layer GCNConv + BatchNorm1d + ReLU + Dropout(0.3) + GlobalMeanPool + 2FC。MATLAB: DLT で同トポロジを実装し、妥当化済み Run-12 では BN を off。 |
| データ | BBBP / MoleculeNet (Wu et al. 2018, 2039 molecules) |
| 論文報告値 | Yang et al. scaffold split: Chemprop 0.919, ECFP 0.877 |
| RP02 ベースライン | ECFP4+sklearn LR random 5-fold CV: AUC = 0.9118 |

---

## Environment (RF02 Version Lock)

実際のバージョンは各 run 後に `result/runs/<timestamp>/lock_snapshot.json` に記録される。

### Python GCN (`rp03_gnn.m`)

| 項目 | 要件 |
|---|---|
| MATLAB | R2025a 以降 |
| Python | 3.10 (EasyMolKit Embedded Python) |
| RDKit | 2022.03 以降 |
| PyTorch | 2.0 以降 (CPU-only) |
| torch_geometric | 2.0 以降 |
| Toolbox | なし。計算本体は `pyrun` 経由で Python 実行 |

### MATLAB DLT GCN (`rp03_gnn_matlab.m`)

| 項目 | 要件 |
|---|---|
| MATLAB | R2025a 以降 + Deep Learning Toolbox |
| Python | 3.10（特徴量生成のみ） |
| RDKit | 2022.03 以降 |
| GPU | 推奨（約 9 分）。CPU でも実行可能（約 30-60 分） |

### 特徴量定義 (RF01 必須)

| 特徴量 | ツール | バージョン | 定義 |
|---|---|---|---|
| 原子特徴量 (25-dim) | RDKit | 実行時記録 | 原子種 onehot x12 + degree onehot x6 + formal charge onehot x5 + `is_aromatic` x1 + H count (0-4 normalized) x1 |

> **注記**: GCNConv は edge feature を直接使わない。結合情報はグラフ構造としてのみ表現される。
> edge feature を明示的に使う MPNN や GATConv は将来 RP 候補。

---

## Data

- **ソース**: MoleculeNet (Wu et al. 2018) / DeepChem distribution（RP02 と同一ファイル）
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv`
- **ライセンス**: MoleculeNet benchmark（research and educational use）
- **キャッシュ**: `data/benchmark/bbbp.csv`（`emk.dataset.bbbp()` により自動取得）
- **件数**: 2039 valid SMILES / 2050 total（RP02 と一致）
- **データハッシュ**: `result/runs/<ts>/lock_snapshot.json` を参照

---

## Scripts

```text
repro/rp03_gnn/rp03_gnn.m              Python GCN オーケストレーション（MATLAB -> Python）
repro/rp03_gnn/rp03_gnn_core.py        Python コア（グラフ構築 + GCN 学習 + 5-fold CV）
repro/rp03_gnn/rp03_gnn_matlab.m       MATLAB DLT GCN オーケストレーション
repro/rp03_gnn/rp03_gnn_matlab_core.m  MATLAB DLT GCN コア（runMatlabGcn + ローカル関数）
repro/rp03_gnn/rp03_gnn_matlab_smoke.m 50 分子・5 epoch のスモークテスト
repro/rp03_gnn/diag1_matlab_forward.m  MATLAB 単分子 forward 診断
repro/rp03_gnn/diag1_python_forward.py diag1 の Python 側 forward 照合
repro/rp03_gnn/diag3_weight_transfer.py MATLAB best-model 検証用の Python 重み転送診断
```

**実行方法**: MATLAB を project root を CWD にして開き、Ctrl+Enter で section ごとに実行する。  
**注記**: 5-fold 全体学習は GPU で約 9 分、CPU で約 30-60 分。

**Python GCN** (`rp03_gnn.m`):

| Section | 内容 |
|---|---|
| Section 0 | Setup (`emk.setup.snapshot()` による RF02 記録) |
| Section 1 | BBBP CSV の解決と RP02 fold indices の自動検出 |
| Section 2 | Python: SMILES->graph 変換、GCN 学習、5-fold CV AUC 評価 |
| Section 3 | RP02 orig / RP02-rev / RP03 の比較レポート |
| Section 4 | RF03 検証 (`emk.repro.verify()`) |
| Section 5 | 可視化（学習曲線、fold AUC 比較、RP02 baseline 線） |
| Section 6 | 結果保存（`makeRunDir()` -> `metrics.json` / `learning_curves.csv` / lockfile） |

**MATLAB DLT GCN** (`rp03_gnn_matlab.m`):

| Section | 内容 |
|---|---|
| Section 0a | パラメータ定義 (`hidden=128`, `lr=1e-3`, `n_epochs=150`, ...) |
| Section 0b | Setup と run directory 作成 |
| Section 1 | RP02 fold indices の自動検出（A2 alignment） |
| Section 2 | Python featurization (`featurize_bbbp` -> `featurized.json`) |
| Section 3 | MATLAB GCN 5-fold CV 学習 |
| Section 4 | Python vs MATLAB 比較レポート |
| Section 5 | 可視化（Python 基準線付き per-fold AUC） |
| Section 6 | 結果保存（`metrics_matlab.json` / `fold_predictions.csv`） |

---

## Results

### Python GCN (PyG / torch_geometric)

| 指標 | 値 | ステータス |
|---|---|---|
| **ROC-AUC CV (5-fold)** | **0.9038 +/- 0.0203** | **PASS (>= 0.85)** |
| vs RP02 sklearn LR（差分） | -0.0080 | 差分 < 1 std（実質比較可能） |

**参照注記**: これは MATLAB 比較で固定参照とする Python 値
（`R2_SPEC.md` の Run-12 Ref。`seed=42`, `batch_size=64`, `dropout=0.3`, `StepLR`, `n_epochs=150`）。
後日の direct rerun ではわずかに値が変わりうるため、役割を明示せずにこの固定参照値と置き換えない。
この固定参照値に対応する per-fold AUC は、後日の rerun の fold 値と混線させないため、
README では再掲しない。

**後日の direct rerun 環境例（上記の固定参照値そのものではない）**:

| 項目 | 値 |
|---|---|
| MATLAB | direct `python_env` rerun では未記録 |
| Python | 3.10.11 |
| RDKit | 2024.03.6 |
| PyTorch | 2.11.0+cpu |
| torch_geometric | 2.7.0 |
| Commit | 4033717 |

### MATLAB DLT GCN

| 指標 | 値 | ステータス |
|---|---|---|
| **ROC-AUC CV (5-fold)** | **0.8872 +/- 0.0151** | **PASS (>= 0.85)** |
| vs Python GCN（差分） | -0.0166 | 差分 < 1 std（実質等価） |
| Runtime (all 5 folds) | ~9 min (GPU) | - |

Per-fold test AUC: 0.8994 / 0.8641 / 0.8796 / 0.8962 / 0.8967

---

## Evaluation Design

### Early Stopping と inner validation split

各 outer fold（RP02 と共有）をさらに 80/20 に分割し、sub-train と inner validation を作る。

- Early stopping は inner validation AUC のみを監視し、test fold は学習中に参照しない
- validation AUC が最大のモデルを `copy.deepcopy` で保存し、early stopping 後に復元する
- 各 test fold は学習完了後に 1 回だけ評価する
- validation split には fold ごとに異なる seed を使う: `train_test_split(random_state=seed+fold_idx)`

> **学習サイズ注記**: RP03 の sub-train は RP02 outer-train の約 80% であり、RP02 の train set より小さい。
> 比較は共有 outer test fold 上では妥当だが、学習サンプル数は同一ではない。
> 正確な件数は `metrics.json -> fold_curves[].train_size` と `val_size` を参照。

### Reproducibility

- モデル初期化 RNG: `torch.manual_seed(seed + fold_idx)`
- DataLoader shuffle も fold-specific generator (`fold_gen`) を使用
- `pos_weight` は sub-train のみから計算し、validation label balance を loss scale に混入させない

---

## Verification (RF03 Numerical Criteria)

> RF03 category: **Cat A (absolute threshold)**。Cat B / Cat C はなし。  
> RF04: compliant（RF01 / RF02 / RF03 Cat A を満たす）

| 指標 | 基準 | 根拠 |
|---|---|---|
| ROC-AUC CV | >= 0.85 | RP02 と同じデータ・同じ split 方式で、GCN は LR+ECFP4 (0.8826) と少なくとも同程度を期待する。Yang et al. の scaffold-split GNN gain (+0.042) を踏まえつつ、保守的に設定。 |

**結果**: 固定参照の Python GCN（0.9038）と MATLAB DLT GCN（0.8872）はどちらも **PASS**。

**許容根拠**:
1. 本再現は GCNConv（undirected, no edge features）であり、論文の directed MPNN より表現力が弱い
2. Random split は scaffold split より高い AUC が出やすい
3. 小規模データ（2039 molecules）では GNN が LR+ECFP4 を下回る場合がある（Hu et al. 2020）

---

## Discussion

### Yang et al. との差異

| 差異 | 内容 |
|---|---|
| Architecture | 論文: directed MPNN (D-MPNN / Chemprop)。本再現: 3-layer GCNConv（undirected, no edge features）。 |
| Split | 論文: Bemis-Murcko scaffold split（train 80% / val 10% / test 10%, 20 runs）。本再現: random 5-fold CV。 |
| Metric | 論文: 20 run の test-set 平均。本再現: 5-fold CV 平均。 |
| Implementation | 論文: Chemprop。本再現: torch_geometric GCNConv。 |

### GNN vs Fingerprints (RP02 との比較)

GCN は分子を原子グラフとして扱い、ECFP4 が捨てる局所環境や接続関係を保持する。
一見すると GCN >= LR+ECFP4 を期待しやすいが、以下の反論がある。

- ECFP4 は n=2039 では過学習しにくく、強く正則化された低次元表現である
- GCN はパラメータ数が多く、小規模データでは汎化で不利になりうる
- MoleculeNet の既報でも ECFP ベースラインは依然として競争力がある

**結果**: GCN 参照 AUC = 0.9038、LR-rev AUC = 0.9118、delta = -0.0080（< 1 std = 0.0203）。  
**BBBP random 5-fold CV（n = 2039）では「FP < GNN」は支持されない（FP approx. GNN）。**

### MATLAB DLT GCN vs Python/PyG

MATLAB DLT GCN（0.8872）は固定参照の Python/PyG GCN（0.9038）から 1 std 以内であり、
この RP では実質等価とみなせる。

主な実装差異:

- Mini-batch 表現: 3D padded tensors（MATLAB） vs PyG DataLoader（Python）
- BatchNorm: MATLAB 側は必須ではない（`BN=off` でも 0.8872 に到達）
- Runtime: MATLAB 約 9 分（GPU） vs Python は後日の direct `python_env` rerun で約 109 秒
- この Python 実行時間は固定参照値と同条件の厳密な比較用タイミングではない。`rp03_gnn.m` には
  CPU 実行の想定時間として約 8-15 分の記載があり、実行環境差を前提に読む必要がある

### Lessons Learned

- **「FP < GNN」は支持されない**: BBBP random CV では固定参照 GCN（0.9038）と LR-rev（0.9118）は近接したまま
- **Early stopping は重要**: 最短の fold は 150 epoch 中 30 epoch で収束した
- **`BCEWithLogitsLoss` と `pos_weight`**: 76.5% / 23.5% の class imbalance 下でも AUC を保てる
- **MATLAB `table()`**: `VariableNames` は文字列リテラルではなく named argument で渡す
- **PyTorch scalar conversion**: `float(loss)` ではなく `loss.item()` を使い `requires_grad` 警告を避ける

---

## MATLAB DLT Implementation Notes

### `pagemtimes` と `reshape`

これは軽微な実装メモではない。Run-1〜Run-11 では `gcnConvBatch` 内の
`reshape(AX, N*B, Fin) * W` が `R2_SPEC.md` の BUG-5 に該当していた。
MATLAB は column-major のため、この reshape が batch/node-to-row 対応を壊し、
初期 MATLAB run を無効化していた。`pagemtimes` へ修正後、AUC は Run-11 の `0.5995`
から妥当化済み Run-12 の `0.8872` へ改善した。

MATLAB DLT GCN 実装では、隣接行列による伝播後の線形変換に `pagemtimes` を使う:

```matlab
AX = pagemtimes(A, X);       % N x Fin x B
H  = pagemtimes(AX, W) + b;  % N x Fout x B
```

`reshape(AX, N*B, Fin) * W` はここでは不適切である。MATLAB は column-major のため、
`reshape(N x Fin x B, N*B, Fin)` はページを interleave し、node/batch と行の対応を崩す。  
`pagemtimes` は `W`（`Fin x Fout`）を全ページに正しく broadcast する。

### `dlarray` 制約

| 制約 | 内容 | 対応 |
|---|---|---|
| sparse 非対応 | `dlarray` は sparse matrix を受け取れない | 前処理で `full()` に変換（BBBP 規模では許容） |
| Graph BN の組み込みなし | `batchNormalizationLayer` は画像/系列向け | 自前 `applyBNLayer` を使用。`BN=off` でも PASS 到達 |

---

## Related Files

| ファイル | 内容 |
|---|---|
| `README.md` | English canonical README |
| `README.jp.md` | 日本語 companion README |
| `rp03_gnn.m` | Python GCN オーケストレーション |
| `rp03_gnn_core.py` | Python コア（validation split + 5-fold CV） |
| `rp03_gnn_matlab.m` | MATLAB DLT GCN オーケストレーション |
| `rp03_gnn_matlab_core.m` | MATLAB DLT GCN コア（runMatlabGcn + local functions） |
| `diag1_matlab_forward.m` | `R2_SPEC.md` Section 9 Step 1 の MATLAB 単分子 forward 診断 |
| `diag1_python_forward.py` | diag1 forward 等価確認の Python 側対応 |
| `diag3_weight_transfer.py` | MATLAB best-model 重み検証と BUG-5 露見に使った Python 診断 |
| `lock_template.json` | RF02 version-lock schema |
| `R2_SPEC.md` | MATLAB DLT 実装の詳細仕様と診断履歴 |
| `result/runs/<ts>/lock_snapshot.json` | 実行時バージョン情報 |
| `result/runs/<ts>/metrics.json` | Python GCN metrics（AUC CV / fold AUCs） |
| `result/runs/<ts>/metrics_matlab.json` | MATLAB GCN metrics |
| `result/runs/<ts>/fold_predictions.csv` | 全 fold の test prediction（`fold / mol_idx / score / label`） |
| `result/runs/<ts>/featurized.json` | MATLAB 用の Python featurization cache |
| `result/runs/<ts>/learning_curves.csv` | epoch ごとの train loss / validation AUC |
| `result/runs/<ts>/fold_auc_matlab.png` | Python 参照線付き per-fold AUC bar chart |
| `repro/rp02_bbbp/` | 比較対象の RP02 ECFP4+LR |
