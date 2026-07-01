# RP00: ESOL Aqueous Solubility Prediction

[日本語版 README](./README.jp.md)

This reproduction revisits the ESOL linear model from Delaney (2004) using EasyMolKit, MATLAB, and RDKit-derived descriptors.
It is the baseline reproduction for the ESOL series: the main goal is to show what can be reproduced with a simple four-descriptor linear model and to document the expected accuracy range on the public MoleculeNet ESOL dataset.

---

## Overview

| Field | Value |
|---|---|
| Paper | Delaney, J.S. (2004). ESOL: Estimating Aqueous Solubility Directly from Molecular Structure. *J. Chem. Inf. Comput. Sci.* 44(3):1000-1005. |
| DOI | [10.1021/ci034243x](https://doi.org/10.1021/ci034243x) |
| Task | Aqueous solubility prediction (`logS`, mol/L) |
| Model | Linear regression with 4 molecular descriptors |
| Data | MoleculeNet ESOL (`1,128` molecules total; regression count may be smaller after filtering) |
| Main takeaway | The reproduced 5-fold CV performance is close to the historical ESOL baseline, but descriptor-definition differences prevent an exact coefficient match |

---

## What This Reproduction Covers

The script reproduces the standard ESOL workflow with these descriptors:

- `LogP`
- `MolWt`
- `NumRotatableBonds`
- `AromaticProportion`

It also documents two practical differences from the original paper:

- Delaney's original `clogP` source differs from RDKit `MolLogP`
- The public MoleculeNet/DeepChem ESOL table contains `1,128` molecules, not the `1,144` reported in the paper

Because of those differences, this RP should be read as a faithful public-dataset reproduction of the ESOL idea, not as a byte-for-byte reconstruction of the original experimental environment.

---

## Environment

Actual runtime versions are written to `result/runs/<timestamp>/lock_snapshot.json` after execution.

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later |
| Python | 3.10 (EasyMolKit Embedded Python) |
| RDKit | 2024.03.6 or later |
| Toolbox | Statistics and Machine Learning Toolbox |

### Descriptor Definitions

| Descriptor | Tool | Definition |
|---|---|---|
| LogP | RDKit `Descriptors.MolLogP` | Crippen-Wildman `MolLogP`; this is not the same as Delaney's original `clogP` source |
| MolWt | RDKit `Descriptors.MolWt` | Average molecular weight including implicit hydrogens (IUPAC atomic weights) |
| NumRotatableBonds | RDKit `rdMolDescriptors.CalcNumRotatableBonds` | Strict SMARTS-based rotatable-bond count |
| AromaticProportion | `pyrun` batch | Aromatic atom count divided by heavy-atom count |

---

## Data

- **Source**: DeepChem / MoleculeNet `delaney-processed.csv`
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv`
- **License**: Public-domain source data from Delaney (2004)
- **Cache path**: `data/benchmark/esol.csv`
- **Total source count**: `1,128` molecules
- **Actual regression count**: depends on the run after SMILES-parse and descriptor-NaN filtering; see `metrics.json` field `n_molecules` and `excluded.csv`

---

## How to Run

Use MATLAB from the project root and run:

```matlab
cd repro/rp00_esol
edit rp00_esol_pilot.m
```

Then execute the script section by section, or run the full script in MATLAB if your local setup is already initialized.

Prerequisites:

- EasyMolKit project dependencies are installed
- Python initialization through `emk.setup.initPython()` works on your machine
- RDKit is available in the configured Python environment

Main script:

```text
repro/rp00_esol/rp00_esol_pilot.m
```

High-level flow:

| Section | Content |
|---|---|
| Section 0 | Setup and environment capture |
| Section 1 | Load the ESOL dataset and compute the dataset hash |
| Section 2 | Parse SMILES and compute descriptors |
| Section 3 | Fit the linear regression model on all molecules |
| Section 4 | Run 5-fold cross-validation |
| Section 5 | Summarize acceptance metrics |
| Section 6 | Save runtime outputs |

---

## Result

Representative first-run result (`2026-06-19`):

| Metric | Value | Note |
|---|---|---|
| RMSE (full dataset fit) | `1.0116` | Descriptive only |
| R^2 (full dataset fit) | `0.7680` | Descriptive only |
| RMSE (5-fold CV, pooled OOF) | about `1.02` | Primary RF03 RMSE criterion |
| RMSE (5-fold CV, fold average) | `1.0166 +/- 0.0243` | Reference summary across folds |
| R^2 (5-fold CV, pooled OOF) | about `0.76` | Main cross-validation signal |

Interpretation:

- The reproduced error level is close to the classic ESOL baseline
- The reproduced `R^2` is lower than the value implied by the original paper because the descriptor source differs, especially for `LogP`
- Coefficient direction is broadly consistent with the paper even when magnitudes differ

Example coefficients from the representative run:

```text
(Intercept)  +0.255
LogP         -0.745
MolWt        -0.0065
RotBonds     +0.0026
AroProp      -0.422
```

---

## Verification

This RP is a baseline calibration run for the ESOL workflow, not a pass/fail verification run by itself.
The practical acceptance range established here becomes the binding RF03 criterion for RP01 onward:

| Metric | Practical target |
|---|---|
| 5-fold CV RMSE | `<= 1.20` |
| 5-fold CV R^2 | `>= 0.75` |

The RMSE target is anchored to Delaney (2004) training RMSE `0.996` with an added tolerance of `0.20` for the RDKit `MolLogP` proxy.
The `R^2` target reflects the observed RDKit-based reproduction level of about `0.76` relative to the paper's roughly `0.84`.
These targets reflect the public dataset and RDKit-based descriptor choices used here, and they should not be interpreted as the exact original-paper thresholds.

---

## Discussion

### Main differences from the original paper

| Difference | Details |
|---|---|
| LogP implementation | Delaney used a different `clogP` source; this RP uses RDKit Crippen-Wildman `MolLogP` |
| Molecule count | Paper: `1,144`; MoleculeNet version used here: `1,128` |
| Evaluation style | The paper emphasizes training-set fit; this RP emphasizes 5-fold CV behavior |
| Rotatable-bond definition | RDKit strict rotatable-bond counting may differ from the paper's implicit definition |

### Why the coefficient match is imperfect

- `LogP` is the largest source of systematic mismatch
- `NumRotatableBonds` is sensitive to definition details and is not statistically important in this run
- Dataset curation differences between the public table and the paper also shift the fit

### Practical takeaway

If you want a compact public ESOL baseline in EasyMolKit, this RP is suitable.
If you want stronger ESOL performance, use later RP variants rather than expecting this exact four-descriptor model to close the remaining gap.

---

## Files

| File | Content |
|---|---|
| `README.md` | English canonical README |
| `README.jp.md` | Japanese companion README |
| `rp00_esol_pilot.m` | Reproduction script |
| `lock_template.json` | Runtime metadata template |

Running the script also creates local outputs under `result/runs/<ts>/`, including metrics, predictions, and environment snapshots.
