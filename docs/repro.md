# Reproducible Research: RP00–RP05

> Full listing of MATLAB reproductions of published Chemoinformatics papers.
> Each entry has a locked environment snapshot (RF02) and defined success criteria (RF03).
> For methods, discussion, and reproduction steps, see the individual `repro/<id>/README.md`.

## Reproduction Index

| ID | Paper (Year) | Task | Method | Dataset | Result | Zone |
|---|---|---|---|---|---|---|
| [RP00](../repro/rp00_esol/) | Delaney (2004) | Regression (solubility) | Linear regression on physicochemical descriptors | ESOL 1,128 cpds | CV RMSE=1.017, R²=0.762 | A |
| [RP01](../repro/rp01_esol/) | Delaney (2004) extended | Regression (solubility) | Linear regression + TPSA / QED / SA Score | ESOL 1,128 cpds | CV RMSE=0.980, R²=0.780 | A |
| [RP02](../repro/rp02_bbbp/) | Wu et al. (2018) MoleculeNet | Classification (BBBP) | ECFP4 (2048-bit) + Logistic Regression (Scaffold 5-fold) | BBBP 2,039 cpds | ROC-AUC CV=0.883 | B/C |
| [RP03](../repro/rp03_gnn/) | Yang et al. (2019) | Classification (BBBP) | Graph Convolutional Network (3-layer) | BBBP 2,039 cpds | ROC-AUC CV=0.915 | C |
| [RP04](../repro/rp04_chemberta/) | Chithrananda et al. (2020) | Classification (BBBP) | ChemBERTa CLS embedding + Logistic Regression | BBBP 2,039 cpds | ROC-AUC CV=0.927 | C |
| [RP05](../repro/rp05_shap/) | SHAP analysis on BBBP | Explainability | shap.LinearExplainer + LR model | BBBP 2,039 cpds | ROC-AUC CV=0.909, Spearman ρ=0.902 | B/D |

## Zone Legend

Zone classification for "What MATLAB Can Cover" — see [README.md](../README.md#what-matlab-can-cover) for the full framework and Zone B conditions.

| Zone | Meaning |
|---|---|
| A | MATLAB native — no Python-side ML needed |
| B | Conditionally equivalent — solver and regularization configured explicitly |
| C | Division of labor — Python for featurization, MATLAB for model training |
| D | Python advantage — not feasible or impractical in MATLAB |

> Some RP entries span multiple zones — e.g., RP02 evaluates both Zone B (MATLAB LR) and Zone C (pipeline).
