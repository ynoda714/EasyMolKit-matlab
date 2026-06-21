# RP05: Explainable AI — SHAP for BBBP ECFP4 Logistic Regression

> **Purpose**: Apply `shap.LinearExplainer` to a BBBP classification LR model
> following Rodriguez-Perez & Bajorath (2020) to interpret ECFP4 bit-level
> feature attributions.

---

## Overview

| Item | Details |
|---|---|
| Paper | Rodriguez-Perez, R. & Bajorath, J. (2020). Interpretation of Machine Learning Models Using Shapley Values: Application to Compound Potency and Multi-target Activity Predictions. *J. Comput.-Aided Mol. Des.* 34:1013–1026. |
| DOI | [10.1007/s10822-020-00314-0](https://doi.org/10.1007/s10822-020-00314-0) |
| Task | Feature attribution / explainable AI (interpreting BBBP classification model) |
| Model | sklearn LogisticRegression (L2 ridge, C=1.0) + Morgan ECFP4 (radius=2, 2048 bits) |
| SHAP | `shap.LinearExplainer` (exact solution for linear models) |
| Data | BBBP (MoleculeNet / DeepChem, 2039 molecules) |
| Published metric | Paper uses ChEMBL data + SHAP analysis across RF/LR/NB/SVM. No directly comparable numeric target (methodology reproduction) |

---

## Environment (RF02 Version Lock)

Actual versions are recorded in `result/runs/<timestamp>/lock_snapshot.json` after execution.

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later |
| Python | 3.10 (EasyMolKit Embedded Python) |
| RDKit | 2022.03 or later |
| shap | 0.49.1 or later (install: `python_env/python.exe -m pip install shap`) |
| scikit-learn | 1.7.0 or later (auto-installed as shap dependency) |
| Toolbox | None required (LR runs in sklearn, no Statistics Toolbox needed) |

### Descriptor Definitions (RF01 required)

| Descriptor | Tool | Version | Definition |
|---|---|---|---|
| Morgan ECFP4 | RDKit | logged | `rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)` |

---

## Data

- **Source**: MoleculeNet (Wu et al. 2018) / DeepChem distribution
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv`
- **License**: MoleculeNet benchmark (research/education use)
- **Cache**: `data/benchmark/bbbp.csv` (shared with RP02)
- **Count**: 2039 valid SMILES / 2050 total (consistent with RP02)
- **Data hash**: see `result/runs/<ts>/lock_snapshot.json`

---

## Script

```
repro/rp05_shap/rp05_shap.m          MATLAB orchestration + visualization
repro/rp05_shap/rp05_shap_core.py    Python core (sklearn LR fit + SHAP)
```

**How to run**: Open MATLAB with project root as CWD and execute sections with Ctrl+Enter.

| Section | Contents |
|---|---|
| Section 0 | Setup, Python init, `emk.setup.snapshot()` |
| Section 1 | BBBP CSV check, path resolution |
| Section 2 | Python: ECFP4 compute → sklearn LR fit → 5-fold CV AUC |
| Section 3 | Extract Python results (global importance, example molecule SHAP) |
| Section 4 | MATLAB visualization (global importance bar, local waterfall) |
| Section 5 | RF03 verification (`emk.repro.verify()`) |
| Section 6 | Save results (`makeRunDir()` → lockfile) |

---

## Result (first run 2026-06-20)

| Metric | Value | Pass/Fail |
|---|---|---|
| ROC-AUC CV (5-fold) | 0.9096 ± 0.0083 | ✅ PASS (>= 0.85) |
| SHAP-LR Spearman ρ | 0.9015 | ✅ PASS (>= 0.90) |

**Environment (from lock_snapshot.json):**

| Item | Value |
|---|---|
| MATLAB | R2026a |
| Python | 3.10 |
| RDKit | 2024.03.6 |
| shap | 0.49.1 |
| scikit-learn | 1.7.2 |
| Commit | 086a1cf |

---

## Verification (RF03 Acceptance Criteria)

| Metric | Criterion | Rationale |
|---|---|---|
| `auc_cv` | >= 0.85 | Same data/method as RP02; standard LR+ECFP4 achievement on BBBP |
| `shap_lr_spearman` | >= 0.90 | LinearExplainer is exact for linear models; rank(mean\|SHAP\|) and rank(\|coef\|×std(X)) are theoretically strongly correlated |

---

## Discussion

### Differences from the Paper

| Difference | Details |
|---|---|
| Dataset | Paper uses ChEMBL proprietary data; this reproduction uses MoleculeNet BBBP |
| Model scope | Paper compares RF/NB/SVM/LR; this reproduction covers LR only (consistent with RP02) |
| Evaluation | Paper uses visual comparison of global importance; this adds Spearman ρ for quantitative verification |

### Lessons Learned

- `shap.LinearExplainer` triggers numba JIT compilation on first run (~30 s); fast thereafter
- `shap_values()` returns a list for multiclass models — requires `isinstance` check for binary case
- MATLAB `jsondecode` automatically converts a JSON 2D array to a double matrix (e.g., 3×20)
- Faithful methodology reproduction on chemogenomics tasks requires ChEMBL access — future candidate

---

## Related Files

| File | Contents |
|---|---|
| `rp05_shap.m` | MATLAB reproduction script |
| `rp05_shap_core.py` | Python: ECFP4 + sklearn LR + shap.LinearExplainer |
| `lock_template.json` | RF02 version lock schema |
| `result/runs/<ts>/lock_snapshot.json` | Runtime version snapshot |
| `result/runs/<ts>/metrics.json` | Evaluation metrics (AUC, Spearman ρ, RF03 pass/fail) |
| `result/runs/<ts>/top_bits_shap.csv` | Top-20 ECFP4 bits: mean\|SHAP\| and LR coefficients |
| `result/runs/<ts>/shap_global_importance.png` | Global importance bar chart |
| `result/runs/<ts>/shap_local_waterfall.png` | Local waterfall (TP / TN / misclassified) |
