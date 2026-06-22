# RP03: Graph Neural Network — GCN on BBBP Classification

> **Goal**: Reproduce Blood-Brain Barrier (BBB) classification using a Graph
> Convolutional Network (GCN) following Yang et al. (2019) Chemprop, and
> compare against the ECFP4+LR baseline from RP02 (AUC CV = 0.883).
> Validates what graph learning gains over fixed-length fingerprints.

---

## Overview

| Item | Content |
|---|---|
| Paper | Yang, K. et al. (2019). Analyzing Learned Molecular Representations for Property Prediction. *J. Chem. Inf. Model.* 59(8):3370–3388. |
| DOI | [10.1021/acs.jcim.9b00237](https://doi.org/10.1021/acs.jcim.9b00237) |
| Task | Binary classification — BBB permeability (BBB+ = permeable, BBB- = non-permeable) |
| Model | 3-layer GCNConv + BatchNorm + GlobalMeanPool + 2FC (torch_geometric) |
| Data | BBBP / MoleculeNet (Wu et al. 2018, 2039 molecules) |
| Published metric | Yang et al. scaffold split: Chemprop 0.919, ECFP 0.877 |
| RP02 baseline | ECFP4+LR random 5-fold CV: AUC = 0.8826 |

---

## Environment (RF02 Version Lock)

Actual versions are recorded in `result/runs/<timestamp>/lock_snapshot.json` after each run.

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later |
| Python | 3.10 (EasyMolKit Embedded Python) |
| RDKit | 2022.03 or later |
| PyTorch | 2.0 or later (CPU-only) |
| torch_geometric | 2.0 or later |
| Toolbox | None (all computation runs in Python via pyrun) |

### Feature Definitions (RF01 required)

| Feature | Tool | Version | Definition |
|---|---|---|---|
| Atom features (25-dim) | RDKit | logged | Atom type onehot ×12 + degree onehot ×6 + formal charge onehot ×5 + is_aromatic ×1 + H count (0-4 normalized) ×1 |

> **Note**: GCNConv does not use edge features (bond information is implicit in graph topology).
> Extensions to MPNN / GATConv with edge features are future RP candidates.

---

## Data

- **Source**: MoleculeNet (Wu et al. 2018) / DeepChem distribution (same file as RP02)
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv`
- **License**: MoleculeNet benchmark (research and educational use)
- **Cache**: `data/benchmark/bbbp.csv` (shared with RP02, auto-downloaded by `emk.dataset.bbbp()`)
- **Count**: 2039 molecules (valid SMILES) / 2050 total (matches RP02)
- **Data hash**: see `result/runs/<ts>/lock_snapshot.json`

---

## Script

```
repro/rp03_gnn/rp03_gnn.m          MATLAB orchestration + visualization
repro/rp03_gnn/rp03_gnn_core.py    Python core (graph construction + GCN + 5-fold CV)
```

**How to run**: Open MATLAB with the project root as CWD and press Ctrl+Enter section by section.
**Note**: GCN training takes 5-10 min on CPU (5 folds × up to 150 epochs with early stopping).

| Section | Content |
|---|---|
| Section 0 | Setup (`emk.setup.snapshot()` for RF02 version capture) |
| Section 1 | Resolve BBBP CSV path + auto-discover A2 fold indices |
| Section 2 | Python: SMILES→graph + GCN training + 5-fold CV AUC |
| Section 3 | Comparison report (RP02 orig / RP02-rev / RP03 orig / RP03-rev audit numbers) |
| Section 4 | RF03 verification (`emk.repro.verify()`) |
| Section 5 | Visualization (learning curves, fold AUC comparison with RP02 reference line) |
| Section 6 | Save results (`makeRunDir()` → metrics.json / learning_curves.csv / lockfile) |

---

## Result (first run 2026-06-21) ⚠️ test-leak present

> **Note**: These results are from the original `rp03_gnn.m` before the A3 audit fix.
> Test-leak and fitclinear bias inflate the apparent GCN advantage.
> See **M-REPRO-AUDIT A3 Fix** below for the fair evaluation.

| Metric | Value | Status |
|---|---|---|
| **ROC-AUC CV (5-fold)** | **0.9151 ± 0.0190** | **✅ PASS (≥ 0.85)** |
| RP02 baseline (ECFP4+LR, biased) | 0.8826 | — |
| GCN vs baseline (Δ) | **+0.0325** ⚠️ | biased |

Per-fold AUCs: 0.9182 / 0.8816 / 0.9121 / 0.9249 / 0.9388
Shortest fold epochs: 39 / 150 (early stopping patience=20 effective)

**Environment (from lock_snapshot.json):**

| Item | Value |
|---|---|
| MATLAB | R2026a |
| Python | 3.10 |
| RDKit | 2024.03.6 |
| PyTorch | 2.11.0+cpu |
| torch_geometric | 2.7.0 |
| Commit | 9733c81 |

---

## M-REPRO-AUDIT A3 Fix (2026-06-21)

**Problem**: The original `rp03_gnn.m / rp03_gnn_core.py` had test-leak and inner-split issues:
- Early stopping used `te_loader` (test fold) each epoch → `fold_aucs` held best test AUC over all epochs.
- `StratifiedKFold.next()` for inner split always took the same first fold (fixed seed), so all outer folds used the same 20% partition pattern as the val set.

**Fix (A3)**:
1. Outer train fold split 80/20 into sub-train / val using `train_test_split(test_size=0.2, stratify=y, random_state=seed+fold_idx)`.
   Each outer fold uses a distinct seed so val sets differ across folds.
2. Early stopping monitors `val_loader` (inner val) AUC only.
3. `copy.deepcopy(model.state_dict())` saves best model; restored after early stopping.
4. Test fold evaluated exactly once after training completes.
5. Also integrates A2: `fold_indices_path` aligns outer folds with RP02-rev.

> **Training size note (A2 asymmetry)**: RP03 sub-train (~80% of RP02 outer train) differs
> in size from RP02 outer train. Paired comparison is valid on shared outer test folds, but
> models are trained on different-sized datasets. See `metrics.json` → `train_size_note`
> and `fold_curves[].train_size` / `val_size` for exact sizes per fold.

| Script | Model | AUC CV | Baseline | Δ | Status |
|---|---|---|---|---|---|
| original | GCN (test-leaked) | 0.9151 ± 0.0190 | fitclinear 0.8826 | +0.0325 | biased |
| current | GCN (leak-fixed) | **0.9038 ± 0.0203** | sklearn LR 0.9118 | **−0.0080** | **fair** |

Per-fold test AUC (A3 fixed): 0.9047 / 0.8670 / 0.9177 / 0.9034 / 0.9263

**Key finding**: The apparent GCN advantage of +0.0325 decomposed as:
- fitclinear solver bias: +0.0292 (resolved in A1)
- test-leak: +0.0113 (0.9151 → 0.9038)
- Fair difference after both fixes: **−0.0080** (gap < 1σ=0.0203, p=0.35, not significant)

> **Implication**: "FP < GNN" is **not supported** by BBBP random 5-fold CV at n=2039.
> GCN and LR+ECFP4 are statistically equivalent on this dataset.

---

## Verification (RF03 Numerical Criteria)

| Metric | Criterion | Rationale |
|---|---|---|
| ROC-AUC CV | ≥ 0.85 | Same data and split method as RP02. GCN expected to match or exceed LR+ECFP4 (0.8826). Conservative threshold accounting for Yang et al. scaffold-split GNN gain (+0.042). |

**Tolerance rationale**:
1. GCN (this repro) vs Chemprop D-MPNN (paper): directed MPNN has higher expressive power, so this repro may score slightly lower
2. Random split vs scaffold split: random split inflates AUC (see RP02 Discussion)
3. Small dataset (2039 molecules): GNN can underperform LR+ECFP4 (Hu et al. 2020)

---

## Discussion

### Key differences from Yang et al.

| Difference | Details |
|---|---|
| Architecture | Paper uses directed MPNN (D-MPNN / Chemprop). This repro uses 3-layer GCNConv (undirected, no edge features). |
| Split | Paper: Bemis-Murcko scaffold split (fixed train 80%/val 10%/test 10%, 20 runs). This repro: random 5-fold stratified CV. |
| Metric | Paper: fixed test set mean over 20 runs. This repro: 5-fold CV mean. |
| Implementation | Paper: Chemprop (custom). This repro: torch_geometric GCNConv. |

### Comparison with RP02 (GNN positioning)

GCN processes molecules as atom graphs, preserving structural context that ECFP4 discards
(relative atom positions, local chemical environment). GCN ≥ LR+ECFP4 is theoretically expected,
but counterarguments exist: ECFP4 is unlikely to overfit at n=2039, and GCN with many parameters
may generalize poorly on small datasets. Prior MoleculeNet benchmarks show ECFP baselines are competitive.

### Lessons learned

- [x] ~~GCN AUC=0.9151 exceeded the ECFP4+LR baseline (0.8826) by +0.0325~~ (pre-A3, biased)
  **A3 fix**: GCN-rev AUC=0.9038 vs LR-rev AUC=0.9118. Δ=−0.0080 (<1σ=0.0203), paired t(4)=−1.07, p=0.35. **"FP < GNN" is not supported on BBBP 5-fold CV (FP ≈ GNN).**
- [x] Early stopping worked: shortest fold converged at epoch 39 (patience=20). Only 39/150 epochs used → fast training (~2 min for 5 folds on CPU).
- [x] BCEWithLogitsLoss `pos_weight = n_neg/n_pos < 1` when BBB+ is majority: down-weights majority class (BBB+) loss so the minority class (BBB-) gets relatively higher gradient contribution. Achieves AUC ~0.90 despite 76.5%/23.5% imbalance.
- [x] MATLAB `table()`: do NOT pass `"VariableNames"` as a string literal — use `VariableNames=value` (named argument syntax) to avoid row-count mismatch errors.
- [x] PyTorch: use `loss.item()` instead of `float(loss)` to avoid UserWarning about requires_grad tensor scalar conversion.
- [x] Bridge to RP04 (ChemBERTa): GCN-rev AUC=0.9038 (A3 fair value) becomes the graph-learning baseline for language model comparison.

---

## Related files

| File | Content |
|---|---|
| `rp03_gnn.m` | MATLAB orchestration script (A3 fixed: leak-free CV + fold alignment). See git log for pre-A3 version. |
| `rp03_gnn_core.py` | Python core (train_test_split inner split + GCN + 5-fold CV) |
| `lock_template.json` | RF02 version lock schema |
| `result/runs/<ts>/lock_snapshot.json` | Runtime version info |
| `result/runs/<ts>/metrics.json` | Metrics (AUC CV / fold AUCs / audit comparison / train_size_note) |
| `result/runs/<ts>/learning_curves.csv` | Per-epoch avg Train Loss / Val AUC (truncated to shortest fold; see `learning_curve_note` in metrics.json) |
| `result/runs/<ts>/learning_curves.png` | Learning curves (dual-axis Loss + AUC) |
| `result/runs/<ts>/fold_auc_comparison.png` | Per-fold AUC bar chart with RP02 reference line |
| `repro/rp02_bbbp/` | RP02 ECFP4+LR baseline (comparison target) |
