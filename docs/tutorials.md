# Tutorial Index — EasyMolKit

Full listing of all tutorial modules across all four layers.
For a layer overview, see [README.md](../../README.md).

---

## Layer 1: Foundation (Base MATLAB only)

| # | Title | What You Learn | Desktop | Online |
|---|---|---|:---:|:---:|
| F01 | Drawing Molecules with SMILES | Molecular representation, SMILES syntax | ✔ | ✔ |
| F02 | Calculating Molecular Properties | MW / LogP / TPSA — meaning and calculation | ✔ | ✔ |
| F03 | Introduction to Fingerprints | Bit vector representation, Morgan vs MACCS | ✔ | ✔ |
| F04 | Comparing Molecules by Similarity | Quantifying Tanimoto / Dice similarity | ✔ | ✔ |
| F05 | Substructure Search | SMARTS pattern matching | ✔ | ✔ |
| F06 | Reading Molecules from Files | SDF / SMILES file operations | ✔ | ✔ |

No additional Toolboxes required. All modules run on MATLAB Online Basic (free tier).

---

## Layer 2: Application Stories (Base MATLAB only)

| # | Title | Domain | Desktop | Online |
|---|---|---|:---:|:---:|
| S01 | Find Relatives of Caffeine | Everyday chemistry | ✔ | ✔ |
| S02 | Drug Filters: Lipinski's Rule of Five | Pharmacology | ✔ | ✔ |
| S03 | Structural Alerts for Hazardous Compounds | Safety | ✔ | ✔ |
| S04 | Introduction to Virtual Screening | Drug discovery | ✔ | ✔ |
| S05 | Unknown Compound Identification Challenge | Forensics | ✔ | ✔ |
| S06 | Search Compounds on PubChem | Databases | ✔ | ✔ |
| S07 | Analyze ChEMBL Activity Data | Drug discovery | ✔ | ✔ |

No additional Toolboxes required. All modules run on MATLAB Online Basic (free tier).

---

## Layer 3: Analytics

| # | Title | Required Toolbox | Topic | Desktop | Online |
|---|---|---|---|:---:|:---:|
| A01 | Chemical Space Mapping by PCA | Statistics and ML Toolbox | Dimensionality reduction | ✔ | ✔ |
| A02 | Molecular Clustering | Statistics and ML Toolbox | Structural similarity | ✔ | ✔ |
| A03 | QSAR Regression | Statistics and ML Toolbox | Property prediction | ✔ | ✔ |
| A04 | Drug Classification | Statistics and ML Toolbox | ML classification | ✔ | ✔ |
| A05 | Neural Network Property Prediction | Deep Learning Toolbox | Neural networks | ✔ | ✔ |
| A06 | Dose-Response Curve Fitting | Curve Fitting Toolbox | Pharmacology | ✔ | ✔ |
| A07 | Scaffold Analysis & R-Group Decomposition | Statistics and ML Toolbox | Medicinal chemistry | ✔ | ✔ |
| A08 | Mass Spectrometry × Cheminformatics | Signal Processing + Statistics and ML Toolbox | MS analysis | ✔ | ✔ |
| A09 | PFAS & Environmental Screening | Optimization Toolbox + Statistics and ML Toolbox (optional) | Environmental | ✔ | ✔ |
| A10 | Lead Optimization | Optimization Toolbox (optional) | Drug discovery | ✔ | ✔ |

All modules run on MATLAB Online Basic (free tier).

---

## Layer 4: Research

| # | Title | Required Toolbox | Desktop | Online |
|---|---|---|:---:|:---:|
| R01 | Large-scale Similarity Screening (GPU) | Parallel Computing Toolbox (GPU) | ✔ | △ (CPU only) |
| R02 | PK/PD Simulation | SimBiology | ✔ | ✔ |
| R03 | Forensic Chemometrics | Statistics and ML Toolbox + Parallel Computing Toolbox | ✔ | ✔ |
| R04 | Protein-Ligand Analysis † | Bioinformatics Toolbox | ✔ | ✔ |
| R05 | Molecular Language Model: SMILES Generation | Deep Learning Toolbox | ✔ | ✔ |
| R06 | REINFORCE Molecular Design | Deep Learning Toolbox + Reinforcement Learning Toolbox | ✔ | ✔ |
| R07 | Metabolomics † | Bioinformatics Toolbox + SimBiology | ✔ | ✔ |
| R08 | Protein-Ligand Docking Simulation ‡ | None (Track 1: meeko + vina + pdbfixer) | ✕ | ✔ |
| R09 | GNN Molecular Property Prediction § | Deep Learning Toolbox | ✔ | ✔ |
| R10 | ChemBERTa Transfer Learning § | Deep Learning Toolbox | ✔ | ✔ |

> **†** Requires `emk.setup.installExtra("biopython")` before first run (Track 1 add-on; independent of MATLAB licensing).
>
> **‡ MATLAB Online only** (not supported on Windows Desktop: vina has no Windows PyPI wheel; pdbfixer's openmm dependency is blocked by Smart App Control).
> Setup: in `main_emk.m`, set `cfg.optionalLibraries.meeko/vina/pdbfixer = true`, then run `installOnline(Config=cfg)`.
>
> **§** Requires the PyTorch + HuggingFace stack. Install in order:
> `emk.setup.installExtra("torch")` → `emk.setup.installExtra("torch_geometric")` → `emk.setup.installExtra("transformers")` → `emk.setup.installExtra("datasets")`.
> R10 requires R09's torch environment as a prerequisite.

---

## Reproducible Research

| ID | Paper | Method | Result |
|---|---|---|---|
| RP00 | Delaney (2004) ESOL — aqueous solubility | Linear regression on physicochemical descriptors | CV RMSE=1.017, R²=0.762 |
| RP01 | Delaney (2004) ESOL — extended | Linear regression + TPSA / QED / SA Score | CV RMSE=0.584, R²=0.906 |
| RP02 | Wu et al. (2018) MoleculeNet BBBP baseline | Morgan FP (ECFP4) + Random Forest | ROC-AUC CV=0.883 |
| RP03 | Yang et al. (2019) GNN on BBBP | Graph Convolutional Network | ROC-AUC CV=0.915 |
| RP04 | Chithrananda et al. (2020) ChemBERTa | Frozen CLS embedding + Logistic Regression | ROC-AUC CV=0.927 |
| RP05 | SHAP explainability on BBBP | shap.LinearExplainer + LR model | ROC-AUC CV=0.909, Spearman ρ=0.902 |

Each entry lives under `repro/rp*/` and includes a MATLAB script, environment lock (`lock_template.json`), and `README.en.md` with full context and results.
