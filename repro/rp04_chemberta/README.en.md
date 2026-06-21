# RP04: ChemBERTa Molecular Language Model — BBBP Linear Probe

> **Goal**: Reproduce BBBP BBB permeability classification using ChemBERTa
> (Chithrananda et al. 2020). Extract frozen CLS embeddings from pre-trained
> ChemBERTa (seyonec/ChemBERTa-zinc-base-v1) and train a logistic regression
> classifier (linear probe), comparing against RP02 (ECFP4+LR) and RP03 (GCN).

---

## Overview

| Item | Content |
|---|---|
| Paper | Chithrananda, S. et al. (2020). ChemBERTa: Large-Scale Self-Supervised Pretraining for Molecular Property Prediction. *arXiv*:2010.09885. |
| DOI | [10.48550/arXiv.2010.09885](https://doi.org/10.48550/arXiv.2010.09885) |
| Task | Binary classification — BBB permeability (BBB+ = permeable, BBB- = non-permeable) |
| Model | `seyonec/ChemBERTa-zinc-base-v1` (44M params, RoBERTa-base, pre-trained on ZINC SMILES) |
| Approach | Frozen CLS embedding (768-dim) → StandardScaler + LogisticRegression (linear probe) |
| Data | BBBP / MoleculeNet (Wu et al. 2018, 2039 molecules) |
| Published metric | ChemBERTa BBBP (scaffold split): ROC-AUC = 0.644–0.746 (model-size dependent) |
| Comparison | RP02 ECFP4+LR: 0.8826 / RP03 GCN: 0.9151 (random 5-fold CV) |

---

## Environment (RF02 Version Lock)

Actual versions are recorded in `result/runs/<timestamp>/lock_snapshot.json` after each run.

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later |
| Python | 3.10 (EasyMolKit Embedded Python) |
| RDKit | 2022.03 or later |
| transformers | 5.0 or later |
| scikit-learn | 1.0 or later (installed as shap dependency) |
| Toolbox | None (all computation in Python via pyrun) |

### Feature Definition (RF01 required)

| Feature | Tool | Version | Definition |
|---|---|---|---|
| ChemBERTa CLS embedding (768-dim) | `seyonec/ChemBERTa-zinc-base-v1` | logged | Last hidden state of the `[CLS]` token from RoBERTa encoder. SMILES tokenized with RoBERTa tokenizer; CLS representation from frozen ChemBERTa final layer. |

---

## Data

- **Source**: MoleculeNet (Wu et al. 2018) / DeepChem distribution (same file as RP02/RP03)
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv`
- **License**: MoleculeNet benchmark (research and educational use)
- **Cache**: `data/benchmark/bbbp.csv` (auto-downloaded by `emk.dataset.bbbp()`)
- **Count**: 2039 molecules (valid SMILES) / 2050 total (matches RP02/RP03)
- **Data hash**: see `result/runs/<ts>/lock_snapshot.json`

---

## Script

```
repro/rp04_chemberta/rp04_chemberta.m          MATLAB orchestration + visualization
repro/rp04_chemberta/rp04_chemberta_core.py    Python core (embedding extraction + LR + 5-fold CV)
```

**How to run**: Open MATLAB with project root as CWD and press Ctrl+Enter section by section.
**Note**: ChemBERTa model is cached locally. If downloading for the first time, ~280MB is needed.
CLS embedding extraction takes ~2-4 min on CPU (forward passes only, no training).

| Section | Content |
|---|---|
| Section 0 | Setup (`emk.setup.snapshot()` for RF02 version capture) |
| Section 1 | Resolve BBBP CSV path |
| Section 2 | Python: load ChemBERTa → extract CLS embeddings → 5-fold CV AUC |
| Section 3 | Three-method comparison display (RP02/03/04) |
| Section 4 | Visualization (3-method AUC bar chart + per-fold AUC) |
| Section 5 | RF03 verification (`emk.repro.verify()`) |
| Section 6 | Save results (`makeRunDir()` → metrics.json / lockfile) |

---

## Result (first run 2026-06-21)

| Metric | Value | Status |
|---|---|---|
| **ROC-AUC CV (5-fold)** | **0.9271 ± 0.0107** | **✅ PASS (≥ 0.85)** |
| vs RP02 ECFP4+LR (0.8826) | **+0.0445** | — |
| vs RP03 GCN (0.9151) | **+0.0120** | — |

Per-fold AUCs: 0.9213 / 0.9252 / 0.9314 / 0.9447 / 0.9128
Runtime: ~38 seconds (model load + embedding extraction + 5-fold LR)

**Environment (from lock_snapshot.json):**

| Item | Value |
|---|---|
| MATLAB | R2026a |
| Python | 3.10 |
| transformers | 5.7.0 |
| scikit-learn | 1.7.2 |
| Commit | 8a5515a |

---

## Verification (RF03 Numerical Criteria)

| Metric | Criterion | Rationale |
|---|---|---|
| ROC-AUC CV | ≥ 0.85 | Same data and split as RP02/RP03. ChemBERTa paper scaffold split results (0.64–0.75) are harder than random split, so ≥ 0.85 is expected with random 5-fold CV. |

**Tolerance rationale**:
1. **Split method**: Paper uses scaffold split (harder); this repro uses random 5-fold CV → AUC is higher
2. **Approach**: Paper uses fine-tuning; this repro uses frozen CLS + LR (linear probe) → quality depends on pre-trained representation
3. **Model**: Paper's ChemBERTa-5M/10M vs `seyonec/ChemBERTa-zinc-base-v1` (44M params) differ in size

---

## Discussion

### Key differences from Chithrananda et al.

| Difference | Details |
|---|---|
| Model | Paper uses ChemBERTa-5M / 10M (PubChem SMILES pre-training). This repro uses `seyonec/ChemBERTa-zinc-base-v1` (44M params, ZINC pre-training). |
| Approach | Paper uses fine-tuning (all parameters updated). This repro uses frozen CLS linear probe. |
| Split | Paper uses scaffold split. This repro uses random 5-fold CV (for consistency with RP02/RP03). |
| Metric | Paper uses fixed test set. This repro uses 5-fold CV mean. |

### Three-method comparison (RP02/03/04)

| Method | Features | Dim | Learning | AUC CV |
|---|---|---|---|---|
| ECFP4+LR (RP02) | Morgan fingerprint | 2048 | structural descriptor + linear | 0.8826 |
| GCN (RP03) | atom graph | — | end-to-end graph learning | 0.9151 |
| ChemBERTa+LR (RP04) | pre-trained SMILES embedding | 768 | frozen features + linear probe | — |

### Lessons learned

- [x] Frozen ChemBERTa CLS embedding + LR achieved AUC=0.9271, outperforming GCN (0.9151) by +0.012. Pre-trained ZINC representations transfer effectively to BBBP.
- [x] Runtime ~38 seconds (vs GCN ~2 min). Forward passes only, no training needed → extremely fast.
- [x] std=0.011 (vs GCN 0.019, ECFP4+LR 0.022). Pre-trained features show more stable fold-to-fold variance.
- [x] `lm_head.*` keys reported as UNEXPECTED — ignore safely. `AutoModel` loads only the RoBERTa encoder; MLM head weights are intentionally excluded.
- [x] Three-method summary: ECFP4+LR (0.883) < GCN (0.915) < ChemBERTa+LR (0.927). Pre-trained LM achieves best performance. ZINC 77M SMILES pre-training transfers effectively even to the small BBBP dataset (2039 molecules).
- [x] M-REPRO-SCALE all 5 RPs (RP01–RP05) complete.

---

## Related files

| File | Content |
|---|---|
| `rp04_chemberta.m` | MATLAB reproduction script |
| `rp04_chemberta_core.py` | Python: ChemBERTa loading + CLS extraction + LR |
| `lock_template.json` | RF02 version lock schema |
| `result/runs/<ts>/lock_snapshot.json` | Runtime version info |
| `result/runs/<ts>/metrics.json` | Metrics (AUC CV / fold AUCs / 3-method comparison) |
| `result/runs/<ts>/method_comparison.png` | 3-method AUC bar chart |
| `result/runs/<ts>/fold_auc_chemberta.png` | Per-fold AUC (ChemBERTa + RP02/RP03 reference) |
| `repro/rp02_bbbp/` | RP02 ECFP4+LR (comparison target) |
| `repro/rp03_gnn/` | RP03 GCN (comparison target) |
