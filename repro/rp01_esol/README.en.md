# RP01: ESOL Classical QSAR + L05 Extended Descriptors — Standard Reproduction (RF01)

> **Purpose**: Formally reproduce the Delaney (2004) ESOL model under the RF01-RF04 framework,
> and verify whether L05 extended descriptors (TPSA / HBD / HBA / FractionCSP3 / QED)
> improve logS prediction accuracy.

---

## Overview

| Field | Value |
|---|---|
| Paper | Delaney, J.S. (2004). ESOL: Estimating Aqueous Solubility Directly from Molecular Structure. *J. Chem. Inf. Comput. Sci.* 44(3):1000-1005. |
| DOI | [10.1021/ci034243x](https://doi.org/10.1021/ci034243x) |
| Task | Aqueous solubility prediction (regression: logS) |
| Model A | Linear regression, 4 features (Delaney original): LogP + MolWt + RotBonds + AroProp |
| Model B | Linear regression, 9 features (L05 extended): Model A + TPSA + HBD + HBA + FractionCSP3 + QED |
| Data | ESOL / MoleculeNet (1,128 molecules) |
| Reported metric | Delaney (2004) Table 2: training RMSE = 0.996 (4-descriptor linear regression) |

---

## Environment (RF02 Version Lock)

Actual versions are recorded in `result/runs/<timestamp>/lock_snapshot.json` after execution.

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later (tested on R2026a) |
| Python | 3.10 (EasyMolKit Embedded Python) |
| RDKit | 2022.03 or later (tested on 2024.03.6) |
| Toolbox | Statistics and Machine Learning Toolbox (`fitlm`, `cvpartition`) |

### Descriptor Definitions (RF01 Required)

| Descriptor | Tool | Version | Definition |
|---|---|---|---|
| LogP | RDKit | To be recorded | Wildman-Crippen `MolLogP`. Different from Delaney's clogP |
| MolWt | RDKit | To be recorded | `Descriptors.MolWt` (IUPAC average atomic weight) |
| NumRotatableBonds | RDKit | To be recorded | `CalcNumRotatableBonds` strict SMARTS definition |
| HeavyAtomCount | RDKit | To be recorded | Number of non-hydrogen atoms (denominator of AromaticProportion) |
| AromaticProportion | RDKit | To be recorded | NumAromaticAtoms / HeavyAtomCount (computed via `pyrun` batch) |
| TPSA | RDKit | To be recorded | `CalcTPSA` (Ertl SMARTS-based) |
| NumHDonors | RDKit | To be recorded | `CalcNumHBD` (Ertl definition) |
| NumHAcceptors | RDKit | To be recorded | `CalcNumHBA` (Ertl definition) |
| FractionCSP3 | RDKit | To be recorded | `CalcFractionCSP3` (Lovering Fsp3) |
| QED | RDKit | To be recorded | `rdkit.Chem.QED.qed` (Bickerton 2012; composite of 8 physicochemical properties) |

> **Note**: `CalcNumAromaticAtoms` is absent in RDKit 2024.03.6.
> AromaticProportion is computed using `GetIsAromatic()` via a `pyrun` batch call (confirmed in RP00).

---

## Data

- **Source**: Delaney (2004) ESOL dataset (hosted by DeepChem/MoleculeNet)
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv`
- **License**: Public domain (no explicit license; widely used benchmark)
- **Cache path**: `data/benchmark/esol.csv` (auto-downloaded by `emk.dataset.esol()` on first run)
- **Count**: 1,128 molecules (MoleculeNet version has 16 fewer than the paper's 1,144)
- **Data hash**: See `dataset_sha256` field in `result/runs/<ts>/lock_snapshot.json`

---

## Script

```
repro/rp01_esol/rp01_esol.m
```

**How to run**: Open MATLAB with the project root as CWD and run with `Ctrl+Enter` (section by section),
or run the full script via the MATLAB MCP tool.

| Section | Content |
|---|---|
| Section 0 | Setup (`emk.setup.snapshot()` → RF02 version recording) |
| Section 1 | ESOL dataset loading (`emk.dataset.esol()`) |
| Section 2 | SMILES parsing + descriptor computation (`batchCalculate` + `pyrun` batch for QED/AroProp) |
| Section 3 | Full-dataset model training (Model A / Model B) |
| Section 4 | VIF analysis (multicollinearity diagnostics, immediately after model fitting) |
| Section 5 | 5-fold CV (seed=42, `cvpartition`; rng re-seeded immediately before) |
| Section 6 | Paired t-test (A vs B; RMSE tests A−B > 0, R² tests A−B < 0) |
| Section 7 | RF03 verification (`emk.repro.verify()`) |
| Section 8 | Save results (`makeRunDir()` → predictions.csv / metrics.json / `emk.setup.lockfile()`) |

---

## Result (First run: 2026-06-19)

### Model A — Delaney Original (4 features)

| Metric | Value | Pass/Fail |
|---|---|---|
| RMSE (full dataset) | 1.0116 | — |
| R² (full dataset) | 0.7680 | — |
| RMSE CV (5-fold) | 1.0166 ± 0.024 | ✅ PASS (≤ 1.20) |
| R² CV (5-fold) | 0.7622 ± 0.022 | ✅ PASS (≥ 0.75) |

### Model B — L05 Extended (9 features)

| Metric | Value | Pass/Fail |
|---|---|---|
| RMSE (full dataset) | 0.9698 | — |
| R² (full dataset) | 0.7877 | — |
| RMSE CV (5-fold) | 0.9798 ± 0.027 | ✅ PASS (≤ 1.20) |
| R² CV (5-fold) | 0.7790 ± 0.023 | ✅ PASS (≥ 0.75) |

### L05 Extension Effect

| Metric | Model A | Model B | Improvement | t(4) | p (one-sided) |
|---|---|---|---|---|---|
| CV RMSE | 1.0166 | 0.9798 | **−0.037 (−3.6%)** | 6.591 | **0.001** ✅ significant |
| CV R² | 0.7622 | 0.7790 | **+0.017 (+1.7pp)** | −7.116 | **0.001** ✅ significant |

> paired t-test (per-fold, n=5, one-sided: RMSE tests A−B > 0, R² tests A−B < 0): L05 extension shows statistically significant improvement on both metrics (M-REPRO-AUDIT B2)

**Environment (from lock_snapshot.json):**

| Item | Value |
|---|---|
| MATLAB | R2026a |
| Python | 3.10 |
| RDKit | 2024.03.6 |
| Commit | 1aac63d |

---

## Verification (RF03 Acceptance Criteria)

> **Criterion origin (M-REPRO-AUDIT B1)**: RMSE ≤ 1.20 and R² ≥ 0.75 were set by
> **post-hoc calibration** from the RP00 pilot result (RMSE=1.017, R²=0.762).
> RP01 is not an evaluation against an independently pre-set criterion; it is the first
> verification run against the RP00-calibrated threshold. This is intentional design
> (pilot → production). From RP02 onward, these thresholds are fixed.

| Metric | Criterion | Rationale |
|---|---|---|
| RMSE CV | ≤ 1.20 | Delaney ref 0.996 (training) + tolerance 0.20 (RDKit LogP ≠ clogP; 16-molecule dataset difference) |
| R² CV | ≥ 0.75 | Realistic ceiling for 5-fold CV with RDKit LogP (set from RP00 result of 0.762) |

**Tolerance rationale (informed by RP00 experience):**
1. **Descriptor implementation difference**: Delaney used clogP; this reproduction uses RDKit Crippen-Wildman LogP. Differences up to 0.5 for polar and large molecules.
2. **Dataset version difference**: Paper 1,144 molecules vs MoleculeNet 1,128 (−16 molecules).
3. **RotBonds definition difference**: RDKit strict SMARTS excludes terminal bonds (may differ from Delaney's definition).

---

## Discussion

### Differences from the Paper

| Difference | Details |
|---|---|
| Reported logS metric | Delaney training RMSE = 0.996 (this reproduction: Model A CV = 1.017; difference +2.1%, within tolerance) |
| RotBonds non-significant | Model A: p=0.85; Model B: p=0.71. Was significant in Delaney. Likely due to definition difference (same finding as RP00) |
| AroProp coefficient sign | Negative in both models (higher aromatic fraction → lower logS). Chemically plausible: aromatic compounds are generally less soluble |

### L05 Descriptor Effects

| Descriptor | Coefficient | p-value | Interpretation |
|---|---|---|---|
| TPSA | −0.016 | 1.3e-08 ✅ significant | Partial-regression sign is negative (sign reversal). Marginal r=+0.123 is positive, but high collinearity with LogP (r=−0.445, VIF=12.12) reverses the sign in the multivariate context. LogP absorbs the polarity signal; TPSA captures the residual. |
| HBA | +0.204 | 6.9e-07 ✅ significant | More H-bond acceptors → stronger water interaction |
| QED | +1.355 | 8.9e-11 ✅ significant | Higher drug-likeness reflects solubility-conscious design |
| FractionCSP3 | −0.439 | 0.009 ✅ significant | Higher sp3 fraction → increased steric bulk and crystallinity → lower solubility |
| HBD | −0.057 | 0.239 ❌ non-significant | Multicollinearity with TPSA and HBA |

### VIF Analysis (multicollinearity, M-REPRO-AUDIT B2)

**Model A (4 features)**: All VIF 1.25–1.49. No multicollinearity.

**Model B (9 features)**:

| Descriptor | VIF | Assessment |
|---|---|---|
| TPSA | 12.12 | ⚠️ HIGH (> 10) |
| LogP | 8.45 | Caution (5–10) |
| HBA | 9.30 | Caution (5–10) |
| MolWt | 7.48 | Caution (5–10) |
| AroProp | 4.61 | Acceptable |
| FractionCSP3 | 4.74 | Acceptable |
| HBD | 3.33 | Acceptable |
| RotBonds | 2.07 | Acceptable |
| QED | 1.14 | No issue |

TPSA's high VIF is the primary cause of the sign reversal noted above.
The L05 extension's statistical significance (paired t-test p=0.001) holds despite
multicollinearity, but directional interpretation of individual coefficients requires caution.

### Lessons Learned (handoff to subsequent RP / M-REPRO-SCALE)

- [x] RF01-RF04 framework confirmed to work for formal reproduction (first production run)
- [x] L05 descriptors (TPSA, HBA, QED) contribute significantly to logS prediction (p < 0.01)
- [x] QED `pyrun` batch call (same round-trip as AroProp) shows no performance issues
- [ ] Large datasets (RP02 MoleculeNet: thousands of molecules) may require batch size limits
- [ ] Non-significance of RotBonds should be verified in RP02 and beyond

---

## Related Files

| File | Content |
|---|---|
| `rp01_esol.m` | Reproduction script (Model A / B comparison) |
| `lock_template.json` | RF02 version lock schema template |
| `result/runs/<ts>/lock_snapshot.json` | Actual version information recorded at runtime |
| `result/runs/<ts>/metrics.json` | Evaluation metrics (all scalars) |
| `result/runs/<ts>/predictions.csv` | Measured values, predictions, and descriptors for all molecules |
| `result/runs/<ts>/predicted_vs_actual.png` | Scatter plot comparison: Model A vs B |
| `result/runs/<ts>/model_b_coefficients.csv` | Model B regression coefficients (TPSA sign verification archive) |
| `repro/rp00_esol/` | RP00 pilot (prototype for RF01-RF04 design) |
