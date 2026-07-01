# RP05: SHAP Reproduction for BBBP ECFP4

[日本語 README](./README.jp.md)

This reproduction evaluates how far EasyMolKit and MATLAB can reproduce a SHAP-based interpretation workflow for BBBP classification on ECFP4 fingerprints.
The primary validated result is the F2 analytical linear-SHAP reproduction. F3 is retained as a separate runtime-feasibility study for MATLAB `shapley()`.

---

## Overview

| Field | Value |
|---|---|
| Paper | Rodriguez-Perez, R. & Bajorath, J. (2020). Interpretation of Machine Learning Models Using Shapley Values: Application to Compound Potency and Multi-target Activity Predictions. *J. Comput.-Aided Mol. Des.* 34:1013-1026. |
| DOI | [10.1007/s10822-020-00314-0](https://doi.org/10.1007/s10822-020-00314-0) |
| Task | SHAP-based interpretation for BBBP classification |
| Model | Logistic Regression, Random Forest, SHAP / TreeSHAP |
| Data | MoleculeNet / DeepChem BBBP (`2,039` valid molecules) |
| Main takeaway | MATLAB reproduces the linear-SHAP path well when the regularization setting is matched; tree-SHAP remains runtime-limited in the current MATLAB workflow |

---

## What This Reproduction Covers

RP05 is a composite reproduction with three layers:

- Python logistic-regression SHAP as the reference path
- F2 analytical MATLAB linear-SHAP reproduction as the primary validated result
- F3 MATLAB `shapley()` vs Python TreeSHAP as a bounded feasibility study

This should be read as a public-dataset reproduction and comparison study, not as a claim that MATLAB fully matches every Python SHAP implementation.

Detailed F3 benchmark design, aggregation policy, and Zone C / D interpretation are documented in [F3_spec.md](./F3_spec.md).
The current implementation, however, runs F3 only in `exploratory` mode with `n_eval=[1, 2, 4]` and `n_repeats=1`; it does not execute the larger recommended grid described in the spec by default.

---

## Environment

Actual runtime versions are written to `result/runs/<timestamp>/lock_snapshot.json` after execution.

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later |
| Python | EasyMolKit Embedded Python |
| Python packages | `shap>=0.49.1`, `scikit-learn>=1.7` |
| Libraries | RDKit, NumPy, SciPy, pandas |
| Toolbox | Statistics and Machine Learning Toolbox |

### Descriptor Definitions

| Descriptor | Tool | Definition |
|---|---|---|
| Morgan ECFP4 | RDKit | Morgan fingerprint with `radius=2` and `nBits=2048` |

---

## Data

- **Source**: MoleculeNet / DeepChem BBBP
- **URL**: Acquired through the EasyMolKit dataset flow on first use
- **License**: Subject to the original dataset provider terms
- **Cache path**: `data/benchmark/bbbp.csv`
- **Count used here**: `2,039` valid molecules after SMILES validation
- **Split used for SHAP comparison**: stratified 80/20 split (`train=1,631`, `test=408`, `seed=42`)

---

## How to Run

Use MATLAB from the project root and run:

```matlab
cd repro/rp05_shap
edit rp05_shap.m
```

Then execute the script section by section.

Prerequisites:

- EasyMolKit project dependencies are installed
- `emk.setup.initPython()` works in your local setup
- The configured Python environment has RDKit and `shap`

Main files:

```text
repro/rp05_shap/rp05_shap.m
repro/rp05_shap/rp05_shap_core.py
repro/rp05_shap/F3_spec.md
repro/rp05_shap/lock_template.json
```

High-level flow:

| File / Section | Content |
|---|---|
| `rp05_shap.m` Section 0 | Setup, Python initialization, environment capture |
| `rp05_shap.m` Section 1 | BBBP dataset resolution and helper checks |
| `rp05_shap.m` Section 2 | Python logistic-regression SHAP reference run |
| `rp05_shap.m` Section 2b | Revised (`task_b`) LR + SHAP comparison |
| `rp05_shap.m` Section 2c | F2 analytical MATLAB linear-SHAP reproduction |
| `rp05_shap.m` Section 2d | F3 MATLAB `shapley()` benchmark and early-stop guardrail |
| `rp05_shap_core.py` | BBBP loading, ECFP4 generation, Python SHAP export |
| `F3_spec.md` | Detailed F3 benchmark and interpretation rules |

---

## Result

Representative verified run: `result/runs/20260629_160127_rp05_shap` (`2026-06-29`)

| Metric | Value | Note |
|---|---|---|
| `auc_cv` | `0.9096 +/- 0.0083` | Baseline LR cross-validation |
| `shap_lr_spearman` | `0.8538` | Python `LinearExplainer` vs LR-weight-derived ranking |
| `task_b.auc_cv_rev` | `0.9143 +/- 0.0079` | Revised LR comparison path |
| `task_b.ranking_spearman_rho` | `0.9149` | Baseline vs revised ranking comparison |
| `task_b.top20_overlap` | `14 / 20` | Overlap of the highest-ranked fingerprint bits |

### F2 Validated Result

| Comparison | Spearman rho |
|---|---:|
| MATLAB analytical SHAP vs Python LR SHAP (`C=1.0`) | `0.9274` |
| MATLAB analytical SHAP vs Python LR SHAP (`C=0.10`) | `0.8355` |

The matched-regularization case shows that the MATLAB analytical formulation reproduces the Python linear-SHAP ranking at a useful level.
This status is recorded in `metrics.json` as `f2_matlab_shap.zone_b_confirmed = true`.

### F3 Runtime-Limited Study

| Item | Value |
|---|---:|
| Configured F3 mode | `exploratory` |
| Configured `n_eval` grid | `[1, 2, 4]` |
| Configured `n_repeats` | `1` |
| Benchmark `n_eval` | `1` |
| Benchmark runtime | `469.7 s` |
| Projected total runtime | `54.8 min` |
| Benchmark rho | `0.5823` |
| MATLAB method | `interventional-tree` |
| Status | `SKIPPED` |

The `n_eval=1` row above is the guardrail benchmark used to estimate the configured exploratory run, not evidence that the larger spec-recommended grid was attempted.
In the current code path, F3 is configured for `n_eval=[1, 2, 4]` with a single repeat, and the script runs a benchmark first before deciding whether to continue with that exploratory grid.
F3 runs a benchmark first and exits early when the projected full runtime exceeds the configured guardrail.
In the representative run, the script stopped intentionally because the projected runtime exceeded `20.0` minutes.

---

## Verification

This RP uses Cat A and Cat C as the primary RF03 gates.
Cat B is retained as reference information for the baseline vs revised comparison.

| Metric | Acceptance rule |
|---|---|
| `auc_cv` | `>= 0.85` |
| `task_b.auc_cv_rev` | `>= 0.85` |
| `shap_lr_spearman` | `>= 0.85` |
| `f2_matlab_shap.zone_b_confirmed` | `true` |

Interpretation:

- The LR baseline remains comfortably above the minimum AUC threshold
- The F2 path is the primary validated MATLAB-side reproduction result
- F3 is not yet a completed zone-evaluable study under the current runtime guardrail
- In addition to the current runtime guardrail, the current exploratory grid `n_eval=[1, 2, 4]` does not meet the spec-side minimum formal Zone-evaluation requirement of `n_eval>=64`

---

## Discussion

### Main differences from the paper

| Difference | Details |
|---|---|
| Reference basis | This RP uses Python SHAP as the operational reference and asks how closely MATLAB can reproduce it |
| F2 scope | Linear-SHAP is reproduced analytically rather than by calling a Python explainer from MATLAB |
| F3 scope | Tree-based SHAP is treated as a feasibility-bounded comparison experiment, not as a claim of full MATLAB equivalence |

### Practical takeaway

- RP05 is already useful if the goal is to validate a MATLAB-side linear-SHAP workflow against a Python reference
- RP05 is still runtime-constrained if the goal is a full tree-SHAP reproduction in MATLAB on 2048-bit BBBP features
- F3 redesign ideas and benchmark-policy details belong in [F3_spec.md](./F3_spec.md), not in the top-level README

---

## Files

| File | Content |
|---|---|
| `README.md` | English canonical README |
| `README.jp.md` | Japanese companion README |
| `rp05_shap.m` | Main MATLAB reproduction script |
| `rp05_shap_core.py` | Python reference workflow and MAT export |
| `F3_spec.md` | Detailed F3 specification |
| `lock_template.json` | Runtime metadata template |

Running the script also creates local outputs under `result/runs/<ts>/`, including metrics, status records, and environment snapshots.
