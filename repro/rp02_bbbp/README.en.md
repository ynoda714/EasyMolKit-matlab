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
| Toolbox | Statistics and Machine Learning Toolbox (`fitclinear`, `perfcurve`, `cvpartition`) |

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
| Section 0 | Setup (`emk.setup.snapshot()` → RF02 version recording) |
| Section 1 | BBBP dataset loading (`emk.dataset.bbbp()`) |
| Section 2 | SMILES parsing + Morgan ECFP4 computation (`batchMorganFP_` → `pyrun` batch, 1 IPC round-trip) |
| Section 3 | Full-dataset logistic regression (`fitclinear` + L2 ridge) |
| Section 4 | 5-fold stratified CV (seed=42, `cvpartition` stratified + `perfcurve` ROC-AUC) |
| Section 5 | RF03 verification (`emk.repro.verify()`) |
| Section 6 | Save results (`makeRunDir()` → predictions.csv / roc_curve.png / metrics.json / `emk.setup.lockfile()`) |

---

## Result (First run: 2026-06-20)

### Valid Molecules

| Item | Value |
|---|---|
| Total molecules | 2,050 |
| SMILES parse errors | 11 (BBB+: 7, BBB-: 4) |
| Valid molecules | 2,039 (BBB+: 1,560 [76.5%], BBB-: 479 [23.5%]) |

### Model Evaluation Metrics

| Metric | Value | Pass/Fail |
|---|---|---|
| ROC-AUC (full dataset, train) | 0.9369 | — |
| Accuracy (full dataset, train) | 0.8960 | — |
| Balanced Accuracy (full dataset, train) | 0.8062 | — |
| **ROC-AUC CV (5-fold)** | **0.8826 ± 0.022** | **✅ PASS (≥ 0.85)** |
| Accuracy CV (5-fold) | 0.8725 ± 0.021 | — |
| Balanced Accuracy CV | 0.7502 ± 0.043 | — |

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
3. **Regularization hyperparameter**: Wu et al.'s lambda setting is not published. This reproduction uses `fitclinear` default (Lambda = 1/n).

> **Note**: If the threshold is adjusted based on first-run results, record the reasoning here.

---

## Discussion

### Key Differences from Wu et al.

| Difference | Details |
|---|---|
| Split method | Wu et al.: Bemis-Murcko scaffold split (train/valid/test). This reproduction: random 5-fold stratified CV. This is the primary cause of the AUC difference |
| Evaluation | Wu et al.: mean of 3 runs on a fixed test set. This reproduction: 5-fold CV mean |
| Regularization | Wu et al.: sklearn LogisticRegression default (C=1.0). This reproduction: `fitclinear` L2 ridge (Lambda = 1/n) |
| FP implementation | Wu et al.: DeepChem circular FP (RDKit-based). This reproduction: RDKit `rdFingerprintGenerator` (new API) |

### Comparison with RP01 (RP02 differentiation)

| Item | RP01 (regression) | RP02 (classification) |
|---|---|---|
| Task | logS prediction (continuous) | BBB permeability (binary) |
| Features | Physicochemical descriptors (4 or 9 dims) | Morgan ECFP4 (2,048 dims) |
| Primary metric | RMSE / R² | ROC-AUC |
| Data | ESOL 1,128 molecules | BBBP 2,039 molecules (class imbalance) |
| Model | `fitlm` (linear regression) | `fitclinear` (high-dimensional LR) |

### Lessons Learned (handoff to subsequent RP / M-REPRO-SCALE)

- [x] Random split AUC CV = 0.8826 vs Wu et al. scaffold split 0.690 → difference +0.19. Numerically confirmed that split method has a large impact on AUC
- [x] Balanced Accuracy CV = 0.7502 (notably lower than Accuracy 0.8725) → effect of 76.5% class imbalance; predicting the majority class BBB+ inflates Accuracy
- [x] Batch FP IPC: `pyrun` retrieves ECFP4 2048-bit FPs for n=2,039 molecules in 1 round-trip (~23 s). Practically faster than RP01-style per-molecule IPC
- [x] Passing a logical matrix to `fitclinear` raises "X must be a numeric matrix" → `double()` conversion required (see Gotcha G41)
- [x] `rf03crit` struct polluted across sessions → resolved by explicit initialization with `struct("field", value)` syntax (see Gotcha G42)
- [ ] RP03 (Graph Learning) bridge: this reproduction's ECFP4 LR AUC=0.8826 serves as the baseline for evaluating GNN improvement
- [ ] `fitclinear` Lambda = 1/n varies slightly across folds (acceptable)

---

## Related Files

| File | Content |
|---|---|
| `rp02_bbbp.m` | Reproduction script (ECFP4 + LR classification) |
| `lock_template.json` | RF02 version lock schema template |
| `result/runs/<ts>/lock_snapshot.json` | Actual version information recorded at runtime |
| `result/runs/<ts>/metrics.json` | Evaluation metrics (AUC / Acc / BalAcc) |
| `result/runs/<ts>/predictions.csv` | True labels, predicted labels, and scores for all molecules |
| `result/runs/<ts>/roc_curve.png` | ROC curve (full-data model, with CV AUC in title) |
| `repro/rp01_esol/` | RP01 regression reproduction (baseline for RP02 differentiation) |
