# RP02: MoleculeNet BBBP Classification Baseline

[日本語 README](./README.jp.md)

This reproduction builds a public EasyMolKit baseline for the MoleculeNet BBBP task from Wu et al. (2018).
The canonical path is a Morgan ECFP4 plus logistic-regression classifier evaluated with nested cross-validation.
It also keeps the diagnostic runs needed to explain why the current sklearn baseline replaced the earlier MATLAB `fitclinear` baseline.

---

## Overview

| Field | Value |
|---|---|
| Paper | Wu, Z. et al. (2018). MoleculeNet: A Benchmark for Molecular Machine Learning. *Chem. Sci.* 9:513-530. |
| DOI | [10.1039/C7SC02664A](https://doi.org/10.1039/C7SC02664A) |
| Task | Binary blood-brain barrier permeability classification (`BBB+` vs `BBB-`) |
| Model | Logistic regression (L2) with Morgan ECFP4 fingerprints |
| Data | MoleculeNet BBBP (`2,039` valid molecules) |
| Main takeaway | Representative validated run `20260701_210031_rp02_bbbp` gives nested-CV ROC-AUC `0.9143 +/- 0.0089` from `rp02_bbbp.m` Section 2 / `metrics.json`, using `outer_seed=42` and `inner_seed=7` |

---

## What This Reproduction Covers

The main script establishes one public baseline and three supporting diagnostics.

- Canonical baseline: sklearn `LogisticRegression` with outer 5-fold stratified CV and inner 3-fold `C` selection
- Historical baseline audit: comparison against the older MATLAB `fitclinear` path
- MATLAB parity study: diagnostic runs showing MATLAB `lbfgs` can approach the sklearn result when it uses its own optimal regularization
- Scaffold-split study: a separate estimate showing how much of the Wu et al. gap comes from split strategy

For GitHub-facing use, the primary result of RP02 is the sklearn nested-CV baseline.
The diagnostic scripts remain in the directory for reproducibility and auditability.

---

## Environment

Actual runtime versions are written to `result/runs/<timestamp>/lock_snapshot.json` after execution.
The representative README values below use `result/runs/20260701_210031_rp02_bbbp/lock_snapshot.json`.

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later |
| Python | 3.10 (EasyMolKit embedded Python) |
| RDKit | 2024.03.6 or later in the representative run |
| Toolbox | Statistics and Machine Learning Toolbox |

### Fingerprint Definition

| Feature | Tool | Definition |
|---|---|---|
| Morgan ECFP4 | RDKit `rdFingerprintGenerator.GetMorganGenerator` | Radius `2`, `2048` bits; radius 2 corresponds to ECFP4 |

---

## Data

- **Source**: DeepChem / MoleculeNet BBBP (`BBBP.csv`)
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv`
- **License**: Distributed by DeepChem; see the upstream dataset and Martins et al. (2012) for original provenance
- **Cache path**: `data/benchmark/bbbp.csv`
- **Count used here**: `2,039` valid molecules (`metrics.json -> n_valid`)
- **Excluded rows**: `11` SMILES parse failures; the total exclusion count is reproducible from the loader, but the per-label breakdown is not currently saved as a dedicated artifact

---

## How to Run

Use MATLAB from the project root and run:

```matlab
cd repro/rp02_bbbp
edit rp02_bbbp.m
```

Then execute the script section by section.

Prerequisites:

- EasyMolKit project dependencies are installed
- `emk.setup.initPython()` works in your local environment
- RDKit is available in the configured Python environment

Main scripts:

```text
repro/rp02_bbbp/rp02_bbbp.m
repro/rp02_bbbp/rp02_sklearn_core.py
repro/rp02_bbbp/a1_diagnosis_run.m
repro/rp02_bbbp/a1_diagnosis.py
repro/rp02_bbbp/r1b_solver.m
repro/rp02_bbbp/r1c_matlab_nested_cv.m
```

High-level flow:

| Section | Content |
|---|---|
| Section 0 | Setup and environment capture |
| Section 1 | Cache and validate the BBBP dataset |
| Section 2 | Run sklearn nested CV from `rp02_sklearn_core.py` |
| Section 2b | Re-run MATLAB `fitclinear` on the same folds for diagnosis |
| Section 2c | Run scaffold-based nested CV for comparison with Wu et al. |
| Section 3 | RF03 verification and comparison logging |
| Section 4 | Save pseudo-ROC artifacts |
| Section 5 | Save metrics, fold indices, and version lock data |

---

## Result

Representative validated baseline run: `result/runs/20260701_210031_rp02_bbbp` (`2026-07-01`)

| Metric | Value | Provenance |
|---|---|---|
| ROC-AUC (nested CV, outer 5-fold mean) | `0.9143 +/- 0.0089` | `rp02_bbbp.m` Section 2, `metrics.json -> auc_cv / auc_cv_std`, `rp02_sklearn_core.py` with `outer_seed=42`, `inner_seed=7` |
| Pseudo-ROC AUC | `0.9130` | `rp02_bbbp.m` Section 4, `metrics.json -> auc_pseudo_roc`; descriptive only |
| Scaffold 5-fold ROC-AUC | `0.8832 +/- 0.0173` | `rp02_bbbp.m` Section 2c, `metrics.json -> scaffold_cv.auc_cv / auc_cv_std` |
| Historical MATLAB `fitclinear` baseline | `0.8826 +/- 0.0220` | `a1_diagnosis_run.m`, artifact `result/runs/20260621_204831_a1_auc_gap/a1_diagnosis.json -> auc_rp02_matlab` and fold spread |
| MATLAB `lbfgs` diagnostic best case | `0.9084 +/- 0.0095` | `r1c_matlab_nested_cv.m` Section 3b fair nested CV, artifact `result/runs/20260626_103836_r1c_nested/r1c_results.json -> section3b` |

Representative runtime environment from `lock_snapshot.json`:

| Item | Value |
|---|---|
| MATLAB | `R2026a` |
| Python | `3.10` |
| RDKit | `2024.03.6` |
| Commit | `54ca766` |

Notes:

- The current public baseline is the `inner_seed=7` path used by `rp02_bbbp.m` and saved in `20260701_210031_rp02_bbbp`.
- `rp02_sklearn_core.py` also documents a historical certified run with `inner_seed=42`, which produced `AUC=0.9118`. Do not mix that historical run with the current README headline value.
- `GCN_AUC_RP03 = 0.9151` in `rp02_bbbp.m` Section 3 is a manually copied comparison constant from RP03. If RP03 is re-run, that value must be updated there before reusing the comparison text.
- `r1c_matlab_nested_cv.m` Section 3a fixed-Lambda sweep is a reference experiment that selects Lambda by directly looking at outer test-fold AUC. Even when it lands on the same numeric value as Section 3b in a particular run, treat Section 3a as optimistic support only, not as an unbiased generalization estimate.

---

## Verification

This RP uses a single RF03 Cat A acceptance rule for the public baseline:

| Metric | Practical target |
|---|---|
| Nested-CV ROC-AUC | `>= 0.85` |

Why this threshold is reasonable:

- The public RP02 baseline uses random stratified nested CV, which is easier than the fixed scaffold-split setting reported by Wu et al.
- The representative sklearn baseline is stable around the low `0.91` range, so `0.85` leaves room for ordinary environment drift without weakening the check too far
- Pseudo-ROC AUC is not used for pass/fail; only the nested-CV mean is authoritative

---

## Discussion

### Main differences from Wu et al.

| Difference | Details |
|---|---|
| Split strategy | Wu et al. report a scaffold-split test result; RP02 uses random nested CV as the canonical baseline and scaffold CV as a diagnostic |
| Evaluation style | Wu et al. report a fixed test-set result; RP02 emphasizes cross-validation stability |
| Hyperparameter handling | RP02 selects sklearn `C` by inner CV instead of relying on a fixed default |
| Fingerprint implementation | RP02 uses RDKit's current Morgan fingerprint generator API |

### Why the old MATLAB baseline was replaced

- The historical A1 audit isolates the original `0.8826` gap at fixed `C=1.0`: `solver_contribution = +0.0191` and `fold_split_contribution = +0.0079`, both from `result/runs/20260621_204831_a1_auc_gap/a1_diagnosis.json`.
- A positive `fold_split_contribution` means the sklearn stratified folds were more favorable than the MATLAB `cvpartition` folds in that audit setting, raising mean AUC by about `0.0079`.
- The fair MATLAB parity check is `r1c_matlab_nested_cv.m` Section 3b, where MATLAB `lbfgs` selects its own `Lambda` by inner CV and reaches `0.9084 +/- 0.0095`, leaving a residual sklearn gap of `+0.0058`.
- `r1c_matlab_nested_cv.m` Section 3a fixed-Lambda sweep is methodologically optimistic because it chooses the best Lambda while looking at outer test-fold performance. Use Section 3b nested CV as the primary comparison result.
- The terms are intentionally not interchangeable: A1 uses `solver_contribution` for a fixed-`C` historical audit, while `rp02_bbbp.m` Section 2b uses `solver_gap` and `regularization_gap` for a different matched-fold diagnostic.

### Practical takeaway

If you need the public RP02 baseline for later comparisons, use the sklearn nested-CV ROC-AUC of `0.9143 +/- 0.0089` from `20260701_210031_rp02_bbbp`.
If you need to understand solver sensitivity, fold sensitivity, or the scaffold-split gap, use the diagnostic artifacts in this directory rather than treating those studies as the primary benchmark.

---

## Files

| File | Content |
|---|---|
| `README.md` | English canonical README |
| `README.jp.md` | Japanese companion README |
| `rp02_bbbp.m` | Main reproduction script |
| `rp02_sklearn_core.py` | Python nested-CV core |
| `r1b_solver.m` | Solver comparison diagnostic |
| `r1c_matlab_nested_cv.m` | MATLAB nested-CV diagnostic |
| `a1_diagnosis_run.m` / `a1_diagnosis.py` | Historical A1 audit scripts |
| `lock_template.json` | Runtime metadata template |

Running the scripts also creates local outputs under `result/runs/<ts>/`, including metrics, pseudo-ROC artifacts, fold indices, diagnostic JSON files, and environment snapshots.
