# RP01: ESOL Classical QSAR + L05 Extended Descriptors

[日本語版 README](README.jp.md)

> **Purpose**: Reproduce the Delaney (2004) ESOL linear QSAR on EasyMolKit and
> evaluate whether the `L05` descriptor extension (`TPSA`, `HBD`, `HBA`,
> `FractionCSP3`, `QED`) improves `logS` prediction.

---

## Overview

| Item | Details |
|---|---|
| Paper | Delaney, J.S. (2004). *ESOL: Estimating Aqueous Solubility Directly from Molecular Structure.* *J. Chem. Inf. Comput. Sci.* 44(3):1000-1005 |
| DOI | [10.1021/ci034243x](https://doi.org/10.1021/ci034243x) |
| Task | Aqueous solubility prediction (regression, `logS`) |
| Model A | 4-descriptor linear regression: `LogP + MolWt + RotBonds + AroProp` |
| Model B | 9-descriptor linear regression: Model A + `TPSA + HBD + HBA + FractionCSP3 + QED` |
| Data | ESOL / MoleculeNet (1,128 molecules) |
| Published metric | Delaney Table 2: training RMSE = 0.996 |

---

## Environment

Actual runtime versions are recorded in `result/runs/<timestamp>/lock_snapshot.json`.

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later (validated on R2026a) |
| Python | 3.10 (EasyMolKit Embedded Python) |
| RDKit | 2022.03 or later (validated on 2024.03.6) |
| Toolbox | Statistics and Machine Learning Toolbox |

### Reproducibility Controls

- MATLAB RNG is initialized with `rng(42, "twister")` in Section 0.
- Before `cvpartition` in Section 5, the script re-applies `rng(42, "twister")` so the 5-fold split stays reproducible even if earlier sections consume RNG state.

### Descriptor Definitions

| Descriptor | Tool | Version | Definition |
|---|---|---|---|
| LogP | RDKit | recorded at runtime | Wildman-Crippen `MolLogP` |
| MolWt | RDKit | recorded at runtime | `Descriptors.MolWt` |
| NumRotatableBonds | RDKit | recorded at runtime | `CalcNumRotatableBonds` strict SMARTS |
| HeavyAtomCount | RDKit | recorded at runtime | number of non-hydrogen atoms |
| AromaticProportion | RDKit | recorded at runtime | `NumAromaticAtoms / HeavyAtomCount` |
| TPSA | RDKit | recorded at runtime | `CalcTPSA` |
| NumHDonors | RDKit | recorded at runtime | `CalcNumHBD` |
| NumHAcceptors | RDKit | recorded at runtime | `CalcNumHBA` |
| FractionCSP3 | RDKit | recorded at runtime | `CalcFractionCSP3` |
| QED | RDKit | recorded at runtime | `rdkit.Chem.QED.qed` |

> **Note**: `AromaticProportion` is computed from RDKit `GetIsAromatic()`.
> Delaney's original `clogP` and this RP's RDKit `MolLogP` are not identical.

---

## Data

- **Source**: Delaney (2004) ESOL dataset, hosted by DeepChem / MoleculeNet
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv`
- **License**: no explicit license statement, but widely used as a public benchmark
- **Cache path**: `data/benchmark/esol.csv`
- **Count**: 1,128 molecules (16 fewer than the paper's 1,144)
- **Data hash**: `dataset_sha256` in `result/runs/<ts>/lock_snapshot.json`

---

## Script

```text
repro/rp01_esol/rp01_esol.m
```

**How to run**: Open MATLAB with the project root as the current working directory and run section by section.

| Section | Content |
|---|---|
| Section 0 | setup, environment capture, RF03 criteria load from `lock_template.json`, and RNG initialization with `rng(42, "twister")` |
| Section 1 | load ESOL data and compute dataset SHA-256 |
| Section 2 | parse SMILES and compute descriptors |
| Section 3 | train Model A / Model B, compute full-dataset fit metrics, extract coefficients, and run TPSA sign verification |
| Section 4 | inspect multicollinearity via VIF |
| Section 5 | re-seed with `rng(42, "twister")` and run 5-fold CV |
| Section 6 | paired t-test for Model A vs Model B |
| Section 7 | RF03 verification |
| Section 8 | save outputs (`metrics.json`, `predictions.csv`, `lock_snapshot.json`) |

### Section 3 diagnostic note

Section 3 includes a `TPSA` sign verification (`B2`) diagnostic. The script compares the marginal Pearson correlation `corr(TPSA, logS)` with the partial regression coefficient sign of `TPSA` in Model B. If the signs reverse, the script logs the correlations of `TPSA` with `LogP`, `HBD`, and `HBA` to support a multicollinearity-based interpretation. The same summary is archived in `metrics.json` under `tpsa_b2`.

---

## Result

Initial run date: 2026-06-19. Latest verified run: 2026-06-30.
Validated artifact: `result/runs/20260630_031306_rp01_esol/`.

### Current validated status

| Metric | Value | Status |
|---|---|---|
| Model A RMSE (full-dataset fit) | 1.0094 | reference to Delaney training RMSE |
| Model B RMSE (full-dataset fit) | 0.9655 | reference |
| Model A RMSE CV | 1.0166 +/- 0.024 | PASS (`<= 1.20`) |
| Model A R^2 CV | 0.7638 +/- 0.022 | PASS (`>= 0.75`) |
| Model B RMSE CV | 0.9798 +/- 0.027 | PASS (`<= 1.20`) |
| Model B R^2 CV | 0.7804 +/- 0.023 | PASS (`>= 0.75`) |
| L05 delta RMSE | -0.0368 | reference |
| L05 delta R^2 | +0.0166 | reference |

### Direct comparison to the paper

| Metric | Paper | This RP | Note |
|---|---|---|---|
| Training RMSE | 0.996 | Model A: 1.0094 | same 4-descriptor structure, different descriptor implementation and dataset variant |
| Training RMSE | not reported for L05 | Model B: 0.9655 | extension model used only in this RP |

### Environment from validated run

| Item | Value |
|---|---|
| MATLAB | R2026a |
| Python | 3.10 |
| RDKit | 2024.03.6 |
| Commit | `a4e3305` |

---

## Verification

Under RF03, Cat A is the PASS gate. The Model A vs Model B comparison is recorded as Cat B reference information.

### Cat A: Absolute thresholds

| Metric | Criterion | Rationale |
|---|---|---|
| RMSE CV | `<= 1.20` | published training RMSE 0.996 plus tolerance for implementation and dataset differences |
| R^2 CV | `>= 0.75` | fixed lower bound based on the earlier RP00 pilot verification |

### Cat B: Relative comparison

| Metric | Model A | Model B | Delta | t(4) | p(one-sided) |
|---|---|---|---|---|---|
| CV RMSE | 1.0166 | 0.9798 | -0.0368 | 6.591 | 0.001 |
| CV R^2 | 0.7638 | 0.7804 | +0.0166 | -7.031 | 0.001 |

**Interpretation**: The L05 extension reduces CV RMSE and increases CV R^2. Both effects are significant in the fold-wise paired t-test.

### Cat B2: TPSA sign verification

This diagnostic is not an RF03 PASS/FAIL gate, but it is part of the interpretation record for Model B. It checks whether the partial `TPSA` coefficient in the multivariate fit keeps the same sign as the marginal `TPSA`-`logS` correlation. When the sign reverses, the script logs `TPSA` correlations with `LogP`, `HBD`, and `HBA` because the most likely explanation is multicollinearity rather than a chemically meaningful isolated effect. See `metrics.json` field `tpsa_b2` and `model_b_coefficients.csv`.

### Tolerance rationale

1. Delaney used `clogP`, while this RP uses RDKit `MolLogP`.
2. The MoleculeNet ESOL table contains 1,128 molecules instead of the paper's 1,144.
3. Rotatable-bond definitions may differ and can affect both coefficients and explained variance.

---

## Discussion

### Differences from the paper

| Difference | Details |
|---|---|
| Evaluation setup | the paper reports training RMSE; this RP uses 5-fold CV as the primary metric, but also records full-dataset fit RMSE for direct comparison |
| LogP implementation | the paper used `clogP`; this RP uses RDKit `MolLogP` |
| Dataset size | paper: 1,144 molecules; this RP: MoleculeNet version with 1,128 |
| Diagnostic scope | this RP adds VIF and TPSA sign-verification diagnostics that were not part of the original paper |

### Main takeaways

- The original 4-descriptor Delaney model is reproducible on EasyMolKit.
- The L05 extension improves predictive accuracy under this implementation and dataset configuration.
- Model B shows substantial multicollinearity around `TPSA`; the `tpsa_b2` diagnostic is included so coefficient-sign interpretation is not detached from that context.

---

## Related Files

| File | Content |
|---|---|
| `README.md` | canonical English README |
| `README.jp.md` | Japanese companion README |
| `rp01_esol.m` | main reproduction script |
| `lock_template.json` | RF02 lock schema and RF03 criteria source |
| `result/runs/20260630_031306_rp01_esol/metrics.json` | validated evaluation metrics, including `tpsa_b2` |
| `result/runs/<ts>/metrics.json` | evaluation metrics for other runs |
| `result/runs/<ts>/predictions.csv` | measured values, predictions, and descriptors |
| `result/runs/<ts>/model_b_coefficients.csv` | Model B coefficient table |
| `result/runs/<ts>/predicted_vs_actual.png` | predicted vs actual scatter plot |
| `repro/rp00_esol/` | earlier pilot RP |
