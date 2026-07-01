# RP04: ChemBERTa BBBP Linear Probe

[日本語 README](./README.jp.md)

This reproduction evaluates whether a frozen ChemBERTa molecular language model can provide a strong BBBP classification baseline in EasyMolKit.
The main validated path is a Python ChemBERTa embedding pipeline followed by logistic regression, and the directory also includes Zone C follow-up variants that move more of the workflow into MATLAB.

---

## Overview

| Field | Value |
|---|---|
| Paper | Chithrananda, S. et al. (2020). ChemBERTa: Large-Scale Self-Supervised Pretraining for Molecular Property Prediction. *arXiv*:2010.09885. |
| DOI | [10.48550/arXiv.2010.09885](https://doi.org/10.48550/arXiv.2010.09885) |
| Task | BBB permeability classification (`BBB+` vs `BBB-`) |
| Model | `seyonec/ChemBERTa-zinc-base-v1` with logistic-regression linear probe |
| Data | MoleculeNet / DeepChem BBBP (`2,039` valid molecules) |
| Main takeaway | Frozen ChemBERTa CLS embeddings give a strong BBBP baseline; the representative AUC is about `0.927`, above the fair RP02 and RP03 baselines |

---

## What This Reproduction Covers

RP04 is a composite reproduction with three practical layers:

- Main path: Python ChemBERTa embedding extraction plus logistic regression in the standard RP04 workflow
- F1-b: Python embeddings with MATLAB `fitclinear` as the classifier
- F1-a: Python tokenization only, then MATLAB ONNX inference plus MATLAB `fitclinear`

This should be read as a public-dataset reproduction and transfer-learning comparison study.
It is not a byte-for-byte recreation of the exact ChemBERTa paper setup because the split protocol, model variant, and downstream training method differ from the original publication.

---

## Environment

Actual runtime versions are written to `result/runs/<timestamp>/lock_snapshot.json` after execution.

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later |
| Python | 3.10 (EasyMolKit Embedded Python) |
| RDKit | 2022.03 or later |
| transformers | 5.0 or later |
| scikit-learn | 1.0 or later |
| onnx | Required for F1-a ONNX export / patch helpers in Python |
| Toolbox | Main RP04 path: none required beyond EasyMolKit Python interop. F1-b requires Statistics and Machine Learning Toolbox (`fitclinear`, `perfcurve`, `cvpartition`). F1-a additionally requires Deep Learning Toolbox (`importNetworkFromONNX`) |

### Feature Definition

| Feature | Tool | Definition |
|---|---|---|
| ChemBERTa CLS embedding | `seyonec/ChemBERTa-zinc-base-v1` | Final-layer `[CLS]` representation from a frozen RoBERTa-style encoder after SMILES tokenization |

---

## Data

- **Source**: MoleculeNet / DeepChem BBBP
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv`
- **License**: MoleculeNet benchmark data for research and educational use
- **Cache path**: `data/benchmark/bbbp.csv`
- **Count used here**: `2,039` valid molecules from `2,050` rows
- **Data hash**: See `dataset_sha256` in `result/runs/<ts>/lock_snapshot.json`

---

## How to Run

Use MATLAB from the project root and open the main script:

```matlab
cd repro/rp04_chemberta
edit rp04_chemberta.m
```

Then execute the script section by section.

Prerequisites:

- EasyMolKit project dependencies are installed
- `emk.setup.initPython()` works on your machine
- The configured Python environment has RDKit, `transformers`, and `scikit-learn`
- F1-a additionally needs the Python `onnx` package
- The ChemBERTa model is available locally or can be downloaded through Hugging Face

Main files:

```text
repro/rp04_chemberta/rp04_chemberta.m
repro/rp04_chemberta/rp04_chemberta_core.py
repro/rp04_chemberta/rp04_chemberta_f1.m
repro/rp04_chemberta/rp04_chemberta_f1a.m
repro/rp04_chemberta/lock_template.json
```

High-level flow for `rp04_chemberta.m`:

| Section | Content |
|---|---|
| Section 0 | Setup and environment capture |
| Section 1 | Load BBBP dataset and resolve paths |
| Section 2 | Load fair RP02 / RP03 baselines |
| Section 3 | Extract ChemBERTa CLS embeddings and run 5-fold CV |
| Section 4 | Validate token length and compare `max_length=128` vs `512` |
| Section 5 | Summarize fair baseline comparison |
| Section 6 | Save visualization figures |
| Section 7 | Run RF03 verification |
| Section 8 | Save metrics and lock snapshot |

---

## Result

Representative validated main-path run: `2026-06-26` (`max_length=512`)
Comparison provenance for this table: `auc_rp02_rev` and `auc_rp03_rev` are the values loaded by that representative run from the latest available fair upstream artifacts at execution time. They are not the hardcoded fallback constants used only when no upstream run is found.

| Metric | Value | Note |
|---|---|---|
| `auc_cv` | `0.9270 +- 0.0107` | Main RP04 ChemBERTa result |
| `auc_rp02_rev` | `0.9143 +- 0.0079` | Fair LR+ECFP4 baseline loaded in the representative run |
| `auc_rp03_rev` | `0.9038 +- 0.0203` | Fair leak-fixed GCN baseline |
| `delta_vs_rp02_rev` | `+0.0128` | Main comparison signal |
| `delta_vs_rp03_rev` | `+0.0232` | Main comparison signal |

Token-length validation summary:

| Setting | Truncated molecules | AUC CV |
|---|---:|---:|
| `max_length=128` | `26 / 2039` (`1.3%`) | `0.9271` |
| `max_length=512` | `0 / 2039` (`0.0%`) | `0.9270` |

The difference between `max_length=128` and `512` is negligible (`delta_auc = -0.0001`), so the original truncation issue was disclosed and then quantitatively bounded.

### F1 Follow-up Status

| Variant | Python role | MATLAB role | AUC CV | Status |
|---|---|---|---:|---|
| F1-b (`2026-06-28`) | Tokenize + embed | `fitclinear` LR | `0.9138 +- 0.0088` | RF03 PASS |
| F1-a (`2026-06-28`, header-only record) | Tokenize only | ONNX inference + `fitclinear` LR | `0.9138 +- 0.0087` | RF03 PASS; current workspace does not retain the corresponding saved artifact |

The current workspace preserves a recorded `metrics.json` artifact for F1-b, but not for F1-a.
For F1-a, the script header records the last observed successful run result (`2026-06-28`), while the saved `metrics.json` / `lock_snapshot.json` for that run are not currently present in the workspace.
The F1-a script is already set up to emit both `metrics.json` and `lock_snapshot.json` on future successful runs.

The RP02 and RP03 comparison baselines are loaded dynamically from the latest available fair runs, so those comparison rows can shift if the upstream reference runs are updated.
If no upstream fair run is found, `rp04_chemberta.m` falls back to hardcoded reference constants (`RP02=0.9118`, `RP03=0.9038`) only as a contingency path.

---

## Verification

This RP is validated primarily through Cat A absolute-threshold checking.
The practical acceptance rule used here is:

| Metric | Acceptance rule |
|---|---|
| `ROC-AUC CV` | `>= 0.85` |

Rationale:

- The original paper reports scaffold-split results, which are harder than this RP's random 5-fold CV setup
- This RP uses a frozen-embedding linear probe rather than full fine-tuning
- The model checkpoint used here is a public ChemBERTa variant, not necessarily the exact paper variant

The representative main-path run passes RF03.
F1-b also has a retained artifact showing RF03 PASS.
F1-a has a script-header record of a `2026-06-28` successful run above the same threshold, but that specific run's saved artifact is not currently retained in the workspace.

---

## Discussion

### Main differences from the paper

| Difference | Details |
|---|---|
| Model variant | The paper discusses ChemBERTa variants such as 5M / 10M pretraining setups, while this RP uses `seyonec/ChemBERTa-zinc-base-v1` |
| Downstream method | The paper focuses on fine-tuning; this RP uses frozen embeddings plus logistic regression |
| Split protocol | The paper reports scaffold split; this RP uses random 5-fold CV for comparability with RP02 and RP03 |
| Evaluation basis | This RP emphasizes fair comparison against revised EasyMolKit baselines rather than headline replication of the original table |

### Practical takeaway

- RP04 shows that a frozen ChemBERTa representation is already competitive on BBBP in this EasyMolKit setting
- On the representative fair comparison, ChemBERTa outperforms the revised LR+ECFP4 and revised GCN baselines in mean AUC
- The F1 variants show that more of the downstream workflow can be moved into MATLAB, though the main public result still relies on Python for the standard path

---

## Files

| File | Content |
|---|---|
| `README.md` | English canonical README |
| `README.jp.md` | Japanese companion README |
| `rp04_chemberta.m` | Main ChemBERTa reproduction script |
| `rp04_chemberta_core.py` | Python embedding, tokenization, and ONNX helper code |
| `rp04_chemberta_f1.m` | F1-b Zone C variant |
| `rp04_chemberta_f1a.m` | F1-a ONNX-based Zone C variant |
| `lock_template.json` | Runtime metadata template |
| `metrics_biased_historical.json` | Historical biased reference values saved in `result/runs/<ts>/`; kept only to prevent loss of audit context, and not for analysis or comparison |

Running the scripts also creates local outputs under `result/runs/<ts>/`, including metrics, figures, and environment snapshots.
