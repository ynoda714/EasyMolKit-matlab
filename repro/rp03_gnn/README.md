# RP03: Graph Neural Network — GCN on BBBP Classification

[日本語版 README](./README.jp.md)

> **Goal**: Reproduce Blood-Brain Barrier (BBB) permeability classification using a
> Graph Convolutional Network (GCN) following Yang et al. (2019) Chemprop, and compare
> against the ECFP4+LR baseline from RP02 (sklearn LR AUC CV = 0.912).
> This RP tests whether graph learning improves over fixed-length fingerprints on BBBP.
> It also compares a MATLAB Deep Learning Toolbox (DLT) GCN implementation against
> the Python/PyG reference.

---

## Overview

| Item | Content |
|---|---|
| Paper | Yang, K. et al. (2019). Analyzing Learned Molecular Representations for Property Prediction. *J. Chem. Inf. Model.* 59(8):3370–3388. |
| DOI | [10.1021/acs.jcim.9b00237](https://doi.org/10.1021/acs.jcim.9b00237) |
| Task | Binary classification — BBB permeability (BBB+ = permeable, BBB- = non-permeable) |
| Model | Python: 3-layer GCNConv + BatchNorm1d + ReLU + Dropout(0.3) + GlobalMeanPool + 2FC. MATLAB: same topology in DLT, with BN optional and disabled in the validated Run-12 pass. |
| Data | BBBP / MoleculeNet (Wu et al. 2018, 2039 molecules) |
| Published metric | Yang et al. scaffold split: Chemprop 0.919, ECFP 0.877 |
| RP02 baseline | ECFP4+sklearn LR random 5-fold CV: AUC = 0.9118 |

---

## Environment (RF02 Version Lock)

Actual versions are recorded in `result/runs/<timestamp>/lock_snapshot.json` after each run.

### Python GCN (`rp03_gnn.m`)

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later |
| Python | 3.10 (EasyMolKit Embedded Python) |
| RDKit | 2022.03 or later |
| PyTorch | 2.0 or later (CPU-only) |
| torch_geometric | 2.0 or later |
| Toolbox | None (all computation runs in Python via `pyrun`) |

### MATLAB DLT GCN (`rp03_gnn_matlab.m`)

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later + Deep Learning Toolbox |
| Python | 3.10 (featurization only) |
| RDKit | 2022.03 or later |
| GPU | Recommended (~9 min); CPU also works (~30–60 min) |

### Feature Definitions (RF01 Required)

| Feature | Tool | Version | Definition |
|---|---|---|---|
| Atom features (25-dim) | RDKit | logged | Atom type onehot x12 + degree onehot x6 + formal charge onehot x5 + `is_aromatic` x1 + H count (0-4 normalized) x1 |

> **Note**: GCNConv does not use edge features; bond information is represented only
> through graph topology. Extensions to MPNN or GATConv with explicit edge features
> are future RP candidates.

---

## Data

- **Source**: MoleculeNet (Wu et al. 2018) / DeepChem distribution (same file as RP02)
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv`
- **License**: MoleculeNet benchmark (research and educational use)
- **Cache path**: `data/benchmark/bbbp.csv` (shared with RP02; auto-downloaded by `emk.dataset.bbbp()`)
- **Count**: 2039 molecules (valid SMILES) / 2050 total (matches RP02)
- **Data hash**: see `result/runs/<ts>/lock_snapshot.json`

---

## Scripts

```text
repro/rp03_gnn/rp03_gnn.m              Python GCN orchestration (MATLAB -> Python)
repro/rp03_gnn/rp03_gnn_core.py        Python core (graph construction + GCN + 5-fold CV)
repro/rp03_gnn/rp03_gnn_matlab.m       MATLAB DLT GCN orchestration
repro/rp03_gnn/rp03_gnn_matlab_core.m  MATLAB DLT GCN core (runMatlabGcn + local functions)
repro/rp03_gnn/rp03_gnn_matlab_smoke.m Quick smoke test (50 molecules, 5 epochs)
repro/rp03_gnn/diag1_matlab_forward.m  MATLAB single-molecule forward diagnostic
repro/rp03_gnn/diag1_python_forward.py Python counterpart for diag1 weight-aligned forward check
repro/rp03_gnn/diag3_weight_transfer.py Python weight-transfer diagnostic for MATLAB best-model verification
```

**How to run**: Open MATLAB with the project root as CWD and execute section by section with Ctrl+Enter.  
**Note**: Full 5-fold training takes about 9 minutes on GPU or 30–60 minutes on CPU.

**Python GCN** (`rp03_gnn.m`):

| Section | Content |
|---|---|
| Section 0 | Setup (`emk.setup.snapshot()` for RF02 version capture) |
| Section 1 | Resolve BBBP CSV path and auto-discover RP02 fold indices |
| Section 2 | Python: SMILES->graph conversion, GCN training, and 5-fold CV AUC |
| Section 3 | Comparison report (RP02 orig / RP02-rev / RP03 audit values) |
| Section 4 | RF03 verification (`emk.repro.verify()`) |
| Section 5 | Visualization (learning curves, fold AUC comparison, RP02 baseline reference line) |
| Section 6 | Save results (`makeRunDir()` -> `metrics.json` / `learning_curves.csv` / lockfile) |

**MATLAB DLT GCN** (`rp03_gnn_matlab.m`):

| Section | Content |
|---|---|
| Section 0a | Parameter definitions (`hidden=128`, `lr=1e-3`, `n_epochs=150`, ...) |
| Section 0b | Setup and run-directory creation |
| Section 1 | Auto-discover RP02 fold indices (A2 alignment) |
| Section 2 | Python featurization (`featurize_bbbp` -> `featurized.json`) |
| Section 3 | MATLAB GCN 5-fold CV training |
| Section 4 | Python vs MATLAB comparison report |
| Section 5 | Visualization (per-fold AUC with Python reference line) |
| Section 6 | Save results (`metrics_matlab.json` / `fold_predictions.csv`) |

---

## Results

### Python GCN (PyG / torch_geometric)

| Metric | Value | Status |
|---|---|---|
| **ROC-AUC CV (5-fold)** | **0.9038 +/- 0.0203** | **PASS (>= 0.85)** |
| vs RP02 sklearn LR (delta) | -0.0080 | gap < 1 std (practically comparable) |

**Reference note**: This is the canonical Python reference used for the validated MATLAB comparison
(`R2_SPEC.md` Run-12 Ref; `seed=42`, `batch_size=64`, `dropout=0.3`, `StepLR`, `n_epochs=150`).
Later direct reruns can differ slightly and should not replace this reference in cross-implementation
comparisons without being labeled explicitly.
The corresponding per-fold AUC values for this pinned reference are not restated here because the
README should not mix them with later rerun fold values.

**Example later direct-rerun environment (not the canonical reference above)**:

| Item | Value |
|---|---|
| MATLAB | not recorded in the direct `python_env` rerun |
| Python | 3.10.11 |
| RDKit | 2024.03.6 |
| PyTorch | 2.11.0+cpu |
| torch_geometric | 2.7.0 |
| Commit | 4033717 |

### MATLAB DLT GCN

| Metric | Value | Status |
|---|---|---|
| **ROC-AUC CV (5-fold)** | **0.8872 +/- 0.0151** | **PASS (>= 0.85)** |
| vs Python GCN (delta) | -0.0166 | gap < 1 std (practically equivalent) |
| Runtime (all 5 folds) | ~9 min (GPU) | - |

Per-fold test AUC: 0.8994 / 0.8641 / 0.8796 / 0.8962 / 0.8967

---

## Evaluation Design

### Early Stopping and Inner Validation Split

Each outer fold (shared with RP02) is further split 80/20 into sub-train and inner validation.

- Early stopping monitors inner validation AUC only; the test fold is never accessed during training
- The best model by validation AUC is saved with `copy.deepcopy` and restored after early stopping
- Each test fold is evaluated exactly once after training completes
- Validation split uses a distinct seed per outer fold: `train_test_split(random_state=seed+fold_idx)`

> **Training-size note**: RP03 sub-train (about 80% of the RP02 outer-train set) is smaller
> than the RP02 training set. The paired comparison remains valid on the shared outer test
> folds, but the models are trained on different effective sample sizes. See
> `metrics.json -> fold_curves[].train_size` and `val_size` for exact counts.

### Reproducibility

- Model initialization RNG: `torch.manual_seed(seed + fold_idx)` provides fold-specific weight initialization
- DataLoader shuffling also uses a fold-specific generator (`fold_gen`)
- `pos_weight` is computed from sub-train only so validation label balance does not leak into the loss scale

---

## Verification (RF03 Numerical Criteria)

> RF03 category: **Cat A (absolute threshold)**. No Cat B or Cat C.  
> RF04: compliant (RF01 / RF02 / RF03 Cat A satisfied)

| Metric | Criterion | Rationale |
|---|---|---|
| ROC-AUC CV | >= 0.85 | Same data and split method as RP02. GCN is expected to be at least comparable to LR+ECFP4 (0.8826). The threshold is set conservatively, accounting for Yang et al.'s scaffold-split GNN gain (+0.042). |

**Result**: The pinned Python reference (0.9038) and MATLAB DLT GCN (0.8872) both **PASS**.

**Tolerance rationale**:
1. This repro uses GCNConv (undirected, no edge features), whereas the paper uses a more expressive directed MPNN
2. Random split tends to yield higher AUC than scaffold split
3. On small datasets (2039 molecules), GNNs can underperform LR+ECFP4 (Hu et al. 2020)

---

## Discussion

### Key Differences from Yang et al.

| Difference | Details |
|---|---|
| Architecture | Paper: directed MPNN (D-MPNN / Chemprop). This repro: 3-layer GCNConv (undirected, no edge features). |
| Split | Paper: Bemis-Murcko scaffold split (train 80% / val 10% / test 10%, 20 runs). This repro: random 5-fold CV. |
| Metric | Paper: mean test-set score over 20 runs. This repro: 5-fold CV mean. |
| Implementation | Paper: Chemprop (custom). This repro: torch_geometric GCNConv. |

### GNN vs Fingerprints (Comparison to RP02)

GCN processes molecules as atom graphs and preserves structural context that ECFP4 discards,
including local chemical environments and relative atom connectivity. A simple expectation would
be GCN >= LR+ECFP4, but several counterarguments apply:

- ECFP4 is unlikely to overfit at n=2039 because it is lower-dimensional and heavily regularized
- GCN has many more parameters and may generalize poorly on small datasets
- Prior MoleculeNet benchmarks already show that ECFP baselines remain competitive

**Outcome**: GCN reference AUC = 0.9038 vs LR-rev AUC = 0.9118, delta = -0.0080 (< 1 std = 0.0203).  
**"FP < GNN" is not supported on BBBP random 5-fold CV at n = 2039 (FP approx. GNN).**

### MATLAB DLT GCN vs Python/PyG

MATLAB DLT GCN (0.8872) is within 1 std of the pinned Python/PyG reference (0.9038), so the two
implementations are practically equivalent for this RP.

Key implementation differences:

- Mini-batch representation: 3D padded tensors (MATLAB) vs PyG DataLoader (Python)
- BatchNorm: not required for the MATLAB pass result (`BN=off` still reaches 0.8872)
- Runtime: MATLAB about 9 min (GPU) vs Python about 109 sec in a later direct `python_env` rerun
- The Python runtime note is not a like-for-like timing baseline for the pinned reference; `rp03_gnn.m`
  still documents an expected CPU runtime of about 8-15 min, so environment differences should be assumed

### Lessons Learned

- **"FP < GNN" is not supported**: the pinned GCN reference (0.9038) and LR-rev (0.9118) remain close on BBBP random CV
- **Early stopping matters**: the shortest fold converged at epoch 30 out of 150
- **`BCEWithLogitsLoss` with `pos_weight`** handles the 76.5% / 23.5% class imbalance without collapsing AUC
- **MATLAB `table()`** should receive `VariableNames` as a named argument, not as a string literal
- **PyTorch scalar conversion** should use `loss.item()` instead of `float(loss)` to avoid `requires_grad` warnings

---

## MATLAB DLT Implementation Notes

### `pagemtimes` vs `reshape`

This was not a minor implementation detail. In Run-1 through Run-11, the MATLAB path used
`reshape(AX, N*B, Fin) * W` inside `gcnConvBatch`, which is BUG-5 in `R2_SPEC.md`.
Because MATLAB is column-major, that reshape broke the batch/node-to-row correspondence and
invalidated the earlier MATLAB runs; after switching to `pagemtimes`, AUC improved from
Run-11 `0.5995` to validated Run-12 `0.8872`.

The MATLAB DLT GCN implementation uses `pagemtimes` for the linear transform after
adjacency-based message aggregation:

```matlab
AX = pagemtimes(A, X);       % N x Fin x B
H  = pagemtimes(AX, W) + b;  % N x Fout x B
```

`reshape(AX, N*B, Fin) * W` is not valid here. Because MATLAB is column-major,
`reshape(N x Fin x B, N*B, Fin)` interleaves pages and breaks the node/batch-to-row mapping.
`pagemtimes` correctly broadcasts `W` (`Fin x Fout`) across all pages.

### `dlarray` Constraints

| Constraint | Details | Workaround |
|---|---|---|
| No sparse support | `dlarray` does not accept sparse matrices | Convert with `full()` during preprocessing (acceptable at BBBP scale) |
| No built-in graph BN | `batchNormalizationLayer` targets image/sequence pipelines | Custom `applyBNLayer`; `BN=off` still achieves PASS |

---

## Related Files

| File | Content |
|---|---|
| `README.md` | English canonical README |
| `README.jp.md` | Japanese companion README |
| `rp03_gnn.m` | Python GCN orchestration |
| `rp03_gnn_core.py` | Python core (validation split + 5-fold CV) |
| `rp03_gnn_matlab.m` | MATLAB DLT GCN orchestration |
| `rp03_gnn_matlab_core.m` | MATLAB DLT GCN core (runMatlabGcn + local functions) |
| `diag1_matlab_forward.m` | MATLAB single-molecule forward diagnostic used in Section 9 Step 1 of `R2_SPEC.md` |
| `diag1_python_forward.py` | Python counterpart for the diag1 forward-equivalence check |
| `diag3_weight_transfer.py` | Python diagnostic used to verify MATLAB best-model weights and expose BUG-5 |
| `lock_template.json` | RF02 version-lock schema |
| `R2_SPEC.md` | Detailed MATLAB DLT implementation notes and diagnostic history |
| `result/runs/<ts>/lock_snapshot.json` | Runtime version info |
| `result/runs/<ts>/metrics.json` | Python GCN metrics (AUC CV / fold AUCs) |
| `result/runs/<ts>/metrics_matlab.json` | MATLAB GCN metrics |
| `result/runs/<ts>/fold_predictions.csv` | Test predictions for all folds (`fold / mol_idx / score / label`) |
| `result/runs/<ts>/featurized.json` | Python featurization cache for MATLAB |
| `result/runs/<ts>/learning_curves.csv` | Per-epoch train loss / validation AUC |
| `result/runs/<ts>/fold_auc_matlab.png` | Per-fold AUC bar chart with Python reference line |
| `repro/rp02_bbbp/` | RP02 ECFP4+LR comparison target |
