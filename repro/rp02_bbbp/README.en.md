# RP02: MoleculeNet BBBP Classification Baseline — Standard Reproduction (RF01)

> **Purpose**: Reproduce the MoleculeNet (Wu et al. 2018) BBBP (Blood-Brain Barrier Permeability)
> classification task under the RF01-RF04 framework, establishing a Morgan ECFP4 +
> Logistic Regression baseline. Differentiates from RP01 (regression): classification task,
> ROC-AUC metric, fingerprint features.

---

## Overview

| Field | Value |
|---|---|
| Paper | Wu, Z. et al. (2018). MoleculeNet: A Benchmark for Molecular Machine Learning. *Chem. Sci.* 9:513-530. |
| DOI | [10.1039/C7SC02664A](https://doi.org/10.1039/C7SC02664A) |
| Task | Binary classification — Blood-Brain Barrier permeability (BBB+ = permeable, BBB- = non-permeable) |
| Model | Logistic Regression (L2 ridge) + Morgan ECFP4 (Wu et al. Table 4 "Circular FP" baseline) |
| Data | BBBP / MoleculeNet (2,039 valid molecules) |
| Reported metric | Wu et al. scaffold split: ROC-AUC = 0.690 (Logreg + Circular FP) |

---

## Environment (RF02 Version Lock)

Actual versions are recorded in `result/runs/<timestamp>/lock_snapshot.json` after execution.

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later |
| Python | 3.10 (EasyMolKit Embedded Python) |
| RDKit | 2022.03 or later |
| Toolbox | Statistics and Machine Learning Toolbox (`perfcurve`) |

### Fingerprint Definition (RF01 Required)

| Descriptor | Tool | Version | Definition |
|---|---|---|---|
| Morgan ECFP4 | RDKit | To be recorded | `rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)` — Circular FP; radius=2 corresponds to ECFP4 |

> **Note**: Uses the `rdFingerprintGenerator` API (new API, RDKit 2022.03+).
> Compatibility with the legacy `AllChem.GetMorganFingerprintAsBitVect` is confirmed in `emk.fingerprint.morgan()`.
> Batch FP computation processes all molecules in a single IPC round-trip via `pyrun`
> (same approach as the QED batch in RP01).

---

## Data

- **Source**: Martins et al. (2012) BBBP dataset (hosted by DeepChem/MoleculeNet)
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv`
- **License**: MIT / CC-BY (DeepChem distribution; original paper: J. Chem. Inf. Model. 2012)
- **Cache path**: `data/benchmark/bbbp.csv` (auto-downloaded by `emk.dataset.bbbp()` on first run)
- **Count**: 2,039 valid molecules (BBB+: 1,560 [76.5%], BBB-: 479 [23.5%]; 11 SMILES parse failures)
- **Data hash**: See `dataset_sha256` field in `result/runs/<ts>/lock_snapshot.json`

---

## Script

```
repro/rp02_bbbp/rp02_bbbp.m
```

**How to run**: Open MATLAB with the project root as CWD and run with `Ctrl+Enter` (section by section).

| Section | Content |
|---|---|
| Section 0 | Setup (Python engine warmup, `emk.setup.snapshot()` → RF02 version recording) |
| Section 1 | BBBP CSV cache check (`emk.dataset.bbbp()`) |
| Section 2 | Python nested CV execution (`rp02_sklearn_core.py`: outer 5-fold × inner 3-fold C selection) |
| Section 3 | RF03 verification (`emk.repro.verify()`, using `aucCV` only) |
| Section 4 | Pseudo-ROC curve generation (concatenated test-fold predictions, `makeRunDir()`) |
| Section 5 | Save results (metrics.json / outer_fold_indices.json / lock_snapshot.json) |

---

## Result (First run: 2026-06-20)

### Valid Molecules

| Item | Value |
|---|---|
| Total molecules | 2,050 |
| SMILES parse errors | 11 (BBB+: 7, BBB-: 4) |
| Valid molecules | 2,039 (BBB+: 1,560 [76.5%], BBB-: 479 [23.5%]) |

### Model Evaluation Metrics (sklearn LR nested CV, 2026-06-21 onward)

| Metric | Value | Pass/Fail |
|---|---|---|
| **ROC-AUC CV (outer 5-fold, nested CV)** | **0.9118 ± 0.0075** | **✅ PASS (≥ 0.85)** |
| Pseudo-ROC AUC (concatenated fold predictions) | See `auc_pseudo_roc` in metrics.json | — |

> **Legacy reference (fitclinear)**: ROC-AUC CV = 0.8826 ± 0.022 (deprecated by M-REPRO-AUDIT A1).
> Retained in code as `FITCLINEAR_AUC_HISTORICAL = 0.8826` for historical comparison only.

**Environment (from lock_snapshot.json):**

| Item | Value |
|---|---|
| MATLAB | R2026a |
| Python | 3.10 |
| RDKit | 2024.03.6 |
| Commit | 31383c0 |

---

## Verification (RF03 Acceptance Criteria)

| Metric | Criterion | Rationale |
|---|---|---|
| ROC-AUC CV | ≥ 0.85 | Higher than Wu et al. scaffold split (0.690) because random split has less distribution shift, making the task easier. Typical range for LR+ECFP4 with random split on BBBP: 0.88–0.92 |

**Tolerance rationale**:

1. **Split method difference (primary factor)**: Wu et al. use Bemis-Murcko scaffold split, where test scaffolds differ from training scaffolds, making generalization harder. This reproduction uses random 5-fold stratified CV, which yields higher AUC.
2. **Class imbalance**: BBB+ accounts for 76.5% of the dataset. ROC-AUC is robust to imbalance, but Accuracy can appear high simply by predicting the majority class.
3. **Regularization hyperparameter**: Wu et al.'s lambda setting is not published. This reproduction selects C by inner 3-fold CV from {0.01, …, 5.0}; the historical `fitclinear` default (Lambda = 1/n, equivalent to C = n) is a reference note only and is no longer used.

> **Note**: If the threshold is adjusted based on first-run results, record the reasoning here.

---

## Discussion

### Key Differences from Wu et al.

| Difference | Details |
|---|---|
| Split method | Wu et al.: Bemis-Murcko scaffold split (train/valid/test). This reproduction: random 5-fold stratified CV. This is the primary cause of the AUC difference |
| Evaluation | Wu et al.: mean of 3 runs on a fixed test set. This reproduction: 5-fold CV mean |
| Regularization | Wu et al.: sklearn LogisticRegression default (C=1.0). This reproduction: sklearn LR L2 ridge; C selected by inner 3-fold CV from {0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0}. See Lambda equivalence below |
| FP implementation | Wu et al.: DeepChem circular FP (RDKit-based). This reproduction: RDKit `rdFingerprintGenerator` (new API) |

### Comparison with RP01 (RP02 differentiation)

| Item | RP01 (regression) | RP02 (classification) |
|---|---|---|
| Task | logS prediction (continuous) | BBB permeability (binary) |
| Features | Physicochemical descriptors (4 or 9 dims) | Morgan ECFP4 (2,048 dims) |
| Primary metric | RMSE / R² | ROC-AUC |
| Data | ESOL 1,128 molecules | BBBP 2,039 molecules (class imbalance) |
| Model | `fitlm` (linear regression) | sklearn `LogisticRegression` (high-dimensional LR) |

### Lambda=1/n rationale and difference from sklearn C=1.0 (M-REPRO-AUDIT B3)

`fitclinear` objective (L2 ridge):
```
minimize:  (1/n) Σ loss(i)  +  Lambda * ||w||² / 2
```

sklearn `LogisticRegression` objective (parameter C):
```
minimize:  Σ loss(i)  +  (n / (2C)) * ||w||²
```

Comparing the two: `Lambda = 1/C` (the `n` terms cancel).
`fitclinear` default `Lambda = 1/n` corresponds to sklearn `C = n = 2039` — extremely weak regularization.
In contrast, Wu et al. / RP05 use sklearn default `C = 1.0` (`Lambda = 1.0` — strong regularization).

For ECFP4 2048-bit sparse high-dimensional features, this difference produces a measurable AUC gap
(A1 diagnosis: ~0.019). **RP02-rev uses nested CV to select C ∈ [0.01, 0.5]**, finding that
much stronger regularization is optimal (best_C per fold: [0.1, 0.5, 0.5, 0.1, 0.5];
AUC improved from 0.8826 → 0.9118).

> The original RP02 used `Lambda=1/n` because Wu et al.'s regularization setting was not published
> and the `fitclinear` default was the most neutral choice available. Nested CV is now the default for all subsequent RPs.

### Lessons Learned (handoff to subsequent RP / M-REPRO-SCALE)

- [x] sklearn LR nested CV AUC = 0.9118 vs Wu et al. scaffold split 0.690 → difference +0.22. fitclinear legacy value (0.8826) underestimated due to solver bias (M-REPRO-AUDIT A1)
- [x] `aucCV` (nested CV mean) and `pseudoAuc` (concatenated fold predictions AUC) have different meanings. RF03 pass/fail uses `aucCV` only (M2 fix)
- [x] Passing a logical matrix to `fitclinear` raises "X must be a numeric matrix" → `double()` conversion required (see Gotcha G41)
- [x] `rf03crit` struct polluted across sessions → resolved by explicit initialization with `struct("field", value)` syntax (see Gotcha G42)
- [x] Batch FP IPC: `pyrun` retrieves ECFP4 2048-bit FPs for n=2,039 molecules in 1 round-trip (practically faster than RP01-style per-molecule IPC)
- [ ] RP03 (Graph Learning) bridge: sklearn LR AUC=0.9118 is the new baseline for evaluating GNN improvement (legacy 0.8826 deprecated)

---

## Related Files

| File | Content |
|---|---|
| `rp02_bbbp.m` | Reproduction script (ECFP4 + LR classification) |
| `lock_template.json` | RF02 version lock schema template |
| `result/runs/<ts>/lock_snapshot.json` | Actual version information recorded at runtime |
| `result/runs/<ts>/metrics.json` | Evaluation metrics (auc_cv, best C per fold, RF03 pass/fail, comparison) |
| `result/runs/<ts>/roc_curve.png` | Pseudo-ROC curve (concatenated test-fold predictions) |
| `repro/rp01_esol/` | RP01 regression reproduction (baseline for RP02 differentiation) |
