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
| Comparison (biased) | RP02 ECFP4+LR: 0.8826 / RP03 GCN: 0.9151 ⚠️ biased (pre-B4) |
| Comparison (B4 fair) | RP02-rev: 0.9118 / RP03-rev: 0.9038 (fair baselines, loaded from metrics.json) |

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
**Note**: ChemBERTa model is cached locally (`~/.cache/huggingface/hub/`).
First download requires ~280 MB. CLS embedding extraction takes ~2–4 min on CPU (forward passes only, no training).

| Section | Content |
|---|---|
| Section 0 | Setup (`emk.setup.snapshot()` → RF02 version capture) |
| Section 1 | Load BBBP dataset and resolve CSV path |
| Section 2 | Load fair baselines (RP02-rev / RP03-rev from metrics.json, B4) |
| Section 3 | ChemBERTa embedding extraction + 5-fold CV (Python) |
| Section 4 | Token length validation (impact of max_length truncation, B4) |
| Section 5 | Three-method comparison report (fair baselines) |
| Section 6 | Visualization (3-method AUC bar chart + per-fold AUC) |
| Section 7 | RF03 verification (`emk.repro.verify()`) |
| Section 8 | Save results (`makeRunDir()` → metrics.json / lockfile) |

---

## M-REPRO-AUDIT B4 Fixes (2026-06-22)

**Issue 1 (hardcoded comparison values)**: Original `rp04_chemberta.m` used biased baselines.
- 0.8826: fitclinear solver bias (A1 confirmed: fair value 0.9118)
- 0.9151: GCN test leakage (A3 confirmed: fair value 0.9038)

**Fix**: `rp04_chemberta.m` dynamically loads `result/runs/*rp02_bbbp*/metrics.json` and
`result/runs/*rp03_gnn*/metrics.json`. Falls back to hardcoded values only when no run exists.

**Issue 2 (max_length unvalidated)**: Whether `max_length=128` is adequate for BBBP SMILES was unknown.
SMILES exceeding 128 tokens lose trailing structure, distorting their embeddings.

**Fix**: `rp04_chemberta_core.py` collects raw token lengths in Phase 1 of `_extract_cls_embeddings()`
and returns them as `token_length_stats` (keys: `n_truncated`, `frac_truncated`).

### B4 Results (2026-06-22)

**Token length validation (max_length=128):**

| Statistic | Value |
|---|---|
| min / max | 3 / 326 |
| mean / p50 | 40.9 / 37 |
| p90 / p95 / p99 | 69 / 80 / 134 |
| Truncated | **26 / 2039 (1.3%) ⚠️** |

1.3% of molecules (26) are truncated at max_length=128. p99=134 shows truncation is confined
to the tail of large molecules. Impact on AUC is expected to be small, but **disclosed as a known limitation**.
To fully avoid truncation, use max_length=512 (requires re-run).

**Fair comparison (dynamic load, confirmed values):**

| Script | Baseline | AUC CV | ChemBERTa advantage |
|---|---|---|---|
| old (deleted) | RP02=0.8826 ⚠️, RP03=0.9151 ⚠️ | 0.9271 | +0.0445 / +0.0120 ⚠️ overestimated |
| `rp04_chemberta.m` (current) | RP02=0.9118, RP03=0.9038 | 0.9271 | **+0.0153 (1.43σ) / +0.0233 (2.18σ)** |

---

## Verification (RF03 Numerical Criteria)

> RF03 category: **Cat A (absolute threshold)**. No Cat B / Cat C.
> RF04: ✅ compliant (RF01 / RF02 / RF03 Cat A satisfied)

| Metric | Criterion | Rationale |
|---|---|---|
| ROC-AUC CV | ≥ 0.85 | Same data and split as RP02/RP03. ChemBERTa paper scaffold split results (0.64–0.75) are harder than random split, so ≥ 0.85 is expected with random 5-fold CV. |

**Tolerance rationale**:
1. **Split method**: Paper uses scaffold split (harder); this repro uses random 5-fold CV → AUC is higher
2. **Approach**: Paper uses fine-tuning; this repro uses frozen CLS + LR (linear probe)
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

### Three-method comparison (fair baselines, M-REPRO-AUDIT B4)

| Method | Features | Dim | Learning | AUC CV (fair) | AUC CV (biased) |
|---|---|---|---|---|---|
| ECFP4+LR (RP02-rev) | Morgan fingerprint | 2048 | structural descriptor + nested CV LR | 0.9118 ± 0.0075 | 0.8826 ⚠️ |
| GCN (RP03-rev) | atom graph | — | end-to-end graph learning (leak-fixed) | 0.9038 ± 0.0203 | 0.9151 ⚠️ |
| ChemBERTa+LR (RP04) | pre-trained SMILES embedding | 768 | frozen features + linear probe | 0.9271 ± 0.0107 | same |

**C3 official verdict (M-REPRO-AUDIT C3, 2026-06-22)**:

| Claim | Paired t-test (df=4, one-sided) | Verdict |
|---|---|---|
| FP < GNN: GCN > LR | t=−1.07, p=0.35 | ❌ **Rejected** |
| LM > FP: ChemBERTa > LR | t=+2.99, p=**0.020** | ✅ **Significantly supported** |
| LM > GNN: ChemBERTa > GCN | t=+1.89, p=0.066 | ⚠️ **Borderline (insufficient support)** |

**Revised conclusion**: "FP < GNN < LM" does not hold.
Correct statement: **"FP ≈ GNN; LM > FP (significant at α=0.05)"**.
"LM > GNN" has a mean difference of +0.0233, but GCN's high fold variance (fold 2 = 0.867) leaves n=5 underpowered.

### Lessons learned

- [x] Frozen ChemBERTa CLS + LR achieved AUC=0.9271, exceeding RP03-rev GCN (0.9038) by +0.0233 and RP02-rev LR+ECFP4 (0.9118) by +0.0153 (fair baselines). ZINC pre-trained representations transfer effectively to BBBP.
- [x] Runtime ~38 seconds (vs GCN ~2 min). Forward passes only, no training → extremely fast.
- [x] std=0.011 (vs RP03-rev GCN 0.019). Pre-trained features show more stable fold-to-fold variance.
- [x] `lm_head.*` keys reported as UNEXPECTED — ignore safely. `AutoModel` loads only the RoBERTa encoder; MLM head weights are intentionally excluded.
- [x] Three-method summary (fair baselines): RP03-rev GCN (0.9038) ≈ RP02-rev LR+ECFP4 (0.9118) < ChemBERTa+LR (0.9271). "FP < GNN < LM" ranking does not hold; the accurate statement is "FP ≈ GNN; LM > FP (α=0.05 significant)".
- [x] M-REPRO-SCALE: all 4 RPs (RP01–RP04) are formal SCALE targets. RP00 is M-REPRO-PILOT; RP05 completed separately outside M-REPRO-SCALE.

---

## Related Files

| File | Content |
|---|---|
| `rp04_chemberta.m` | MATLAB script (dynamic metrics.json load + token length validation) |
| `rp04_chemberta_core.py` | Python: ChemBERTa + CLS extraction + LR + `token_length_stats` (B4) |
| `lock_template.json` | RF02 version lock schema |
| `result/runs/<ts>/lock_snapshot.json` | Runtime version info |
| `result/runs/<ts>/metrics.json` | Metrics (AUC CV / fold AUCs / fair comparison / token_length_stats) |
| `result/runs/<ts>/method_comparison_fair.png` | 3-method AUC bar chart (fair baselines) |
| `result/runs/<ts>/fold_auc_fair.png` | Per-fold AUC (ChemBERTa vs rev reference lines) |
| `repro/rp02_bbbp/rp02_bbbp.m` | RP02 (sklearn LR fair baseline, A1 fix applied) |
| `repro/rp03_gnn/rp03_gnn.m` | RP03 (GCN leak-fixed, A3 fix applied) |
