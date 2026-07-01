# RP<XX>: <Title> — Standard Reproduction Template (RF01)

> **Purpose**: RF01 standard reproduction template. Every reproduction study must follow this structure.
> Established in M-REPRO-FOUND (2026-06-19) based on lessons from RP00. RF03/RF04 were formally defined in M-REPRO-AUDIT C1/C2 (2026-06-22).
> Use `README.md` as the English canonical README and `README.jp.md` as the Japanese companion, with mutual links between them.
> The two README files should, in principle, remain symmetric in section order, table granularity, and covered content.

> **How to find candidate papers**: AnyResearch (with OpenAlex integration) can be used
> to systematically search for papers suitable for chemoinformatics reproduction work.
> See AnyResearch `docs/workflows/repro_discovery.md` for the workflow.

---

## Reproduction Framework Summary (RF00–RF04)

| # | Name | Definition | Main artifact |
|---|---|---|---|
| RF00 | Paper selection criteria | Scoring policy for choosing target papers | `docs/algorithms/repro_selection.md` |
| RF01 | Reproduction template | This document (README structure and descriptor definition table) | `repro/Template.md` / `repro/Template.jp.md` |
| RF02 | Version lock | Record runtime environment in `lock_snapshot.json` | `emk.setup.snapshot()` / `emk.setup.lockfile()` |
| RF03 | Verification criteria | Three-category validation of evaluation metrics (see below) | `emk.repro.verify()` |
| RF04 | External reproduction protocol | Procedure for external contributors to submit a reproduction (see below) | `CONTRIBUTING.md` |

### RF03 Verification Categories (M-REPRO-AUDIT C1)

| Category | Name | Content | Impact on RP PASS |
|---|---|---|---|
| **Cat A** | Absolute thresholds | Check upper/lower bounds with `emk.repro.verify()` (RMSE ≤ X, R² ≥ Y, AUC ≥ Z, etc.) | **Required**. Every metric must satisfy the criterion |
| **Cat B** | Relative comparison | Paired t-test and p-value based significance test between two conditions (e.g., Model A vs B) | **Reference only**. Not a direct PASS/FAIL gate |
| **Cat C** | Implementation sanity | Mathematical correctness sanity checks (e.g., SHAP Spearman ρ) | **Required when applicable** |

> **RP PASS rule**: All Cat A checks must pass, and all Cat C checks must pass when applicable. Cat B is recorded as contextual information for comparative studies.
>
> **Cat B policy for `_rev` scripts (design fixed)**: Before/after comparisons reported by revised scripts (`_rev`) via `logInfo` are treated as Cat B reference information. They do not become PASS/FAIL gates. The purpose of a revision is to restore fairness (remove bias, fix leakage), not necessarily to improve headline accuracy, so PASS is retained even if revised AUC is lower than the original.

### RF04 External Reproduction Protocol (M-REPRO-AUDIT C2)

**Definition of RF04 compliance**: An RP is recognized as externally reproducible only when all three requirements below are met.

| Requirement | Content | How to confirm |
|---|---|---|
| RF01 ✓ | `README.md` and `README.jp.md` exist in this template format | Overview / Environment / Data / Script / Result / Verification / Discussion are all filled in for both languages |
| RF02 ✓ | `lock_snapshot.json` exists under `result/runs/<ts>/` | `emk.setup.lockfile()` has been called |
| RF03 ✓ | Cat A (and Cat C when applicable) passed | `metrics.json` contains `rf03_pass: true` |

> For the procedure external contributors should follow when submitting a new reproduction, see `CONTRIBUTING.md` under "Submitting a Reproduction (RF04)".

---

## Script Authoring Checklist (developer-only; do not include in per-RP README)

> **Basis**: A cross-review of fixes across RP00–RP07 showed that the following 10 categories of issues recurred independently in multiple RP scripts (extracted on 2026-06-25).
> Check this list immediately after scaffolding a new RP script, and satisfy it from the initial implementation.

### A — Reproducibility and RNG management

- [ ] Place `rng(seed)` immediately before `cvpartition` (declaring it only in Section 0 is insufficient)
- [ ] At the start of each Python fold loop: `torch.manual_seed(seed+fold_idx)` and `np.random.seed(seed+fold_idx)`
- [ ] For DataLoader: specify `generator=torch.Generator().manual_seed(seed+fold_idx)`
- [ ] Immediately before model loading (`from_pretrained`, etc.): place `random.seed(seed)`

### B — NaN and data-quality guards

- [ ] Exclude rows with any missing descriptor in one pass: `any(ismissing(descTbl), 2)`
- [ ] On the Python side: check `notna()` before `MolFromSmiles`
- [ ] For per-molecule conversion failures: initialize with NaN + `try/except` + `logWarn` fallback

### C — File handles and figure management

- [ ] Wrap `fopen/fprintf/fclose` in `try/finally` (or replace with `writelines`)
- [ ] Immediately after `fopen`, guard with `if fid == -1; error(...); end`
- [ ] After `saveas`, call `close(fig)` to avoid accumulating windows during repeated section runs
- [ ] Guarantee figure close with `try/finally`

### D — Metric definition precision

- [ ] For CV R², baseline = `mean(y_all)` (global). Per-fold local mean is not allowed
- [ ] Training RMSE = `sqrt(RSS/n)` (biased). `sqrt(mdl.MSE)` is not allowed because it uses RSS/(n−p−1)
- [ ] RF03 pass condition is the logical AND across relevant models: `resA.pass && resB.pass`
- [ ] Compute OOF RMSE after concatenating predictions from all folds (not as a mean of fold metrics)

### E — Path and `runDir` consistency

- [ ] Use the return value of `makeRunDir` as `runDir`; do not introduce `absRunDir`
- [ ] At the start of any section that references `runDir`, guard with `if ~exist('runDir','var'); error(...); end`
- [ ] Call `makeRunDir` at the start of the result-saving section (not in a visualization section)

### F — `batchMorganFP_` / IPC validation

- [ ] Pass IPC input as a list object directly (concatenated-string IPC is forbidden)
- [ ] On the Python side: `assert len(bits) == n_molecules * n_bits` before reshape
- [ ] On the MATLAB side: add a per-molecule length check before reshape

### G — Data leakage prevention (XAI / CV workflows)

- [ ] For SHAP `LinearExplainer`, use `X_tr` only as background (never mix in test data)
- [ ] Pass an unfitted estimator to `cross_val_score`
- [ ] Explicitly separate CV scope and evaluation scope (for example, 80/20 split -> CV on train side only)

### H — Hyperparameter explicitness

- [ ] For `fitrlinear`, explicitly set `Lambda=1.0` (avoid MATLAB-version-dependent defaults)
- [ ] For `TreeBagger`, explicitly set `'OOBPrediction','off'` on all models
- [ ] Pass all model hyperparameters explicitly; do not rely on defaults

### I — Python core robustness

- [ ] Catch exceptions with `try/except` and always return JSON like `{"success": false, "error": "..."}`
- [ ] For PyTorch, standardize on `.detach().numpy()` (PyTorch >= 2.0 compatible)
- [ ] If `max_len` is variable, auto-scale with `batch_size = max(1, 32*128//max_len)`

### J — README synchronization

- [ ] Create the English README (`README.md`) and the Japanese README (`README.jp.md`) at the same time
- [ ] Keep heading order, table columns, and covered content symmetric between languages
- [ ] Reflect post-run Result values and record the revision date
- [ ] Keep the section table aligned with the actual code structure (section numbers and descriptions must match)
- [ ] Declare all RF03/RF04 metric fields in `lock_template.json` (single source of truth)

---

## Overview

| Field | Value |
|---|---|
| Paper | Author(s) (Year). Title. *Journal* Vol(Issue):Pages. |
| DOI | [10.xxxx/yyyy](https://doi.org/10.xxxx/yyyy) |
| Task | Regression / Classification / Descriptor calculation / Other |
| Model | Model name / method |
| Data | Dataset name (number of molecules) |
| Reported metric | Evaluation metric and value as reported in the paper |

---

## Environment (RF02 Version Lock)

Actual versions are recorded in `result/runs/<timestamp>/lock_snapshot.json` after execution.

| Item | Requirement |
|---|---|
| MATLAB | R2025a or later |
| Python | 3.10 (EasyMolKit Embedded Python) |
| RDKit | Specify version (e.g., 2022.03 or later) |
| Toolbox | List required toolboxes |

### Descriptor Definitions (RF01 Required)

For each descriptor used in the reproduction, specify the following three items.

| Descriptor | Tool | Version | Definition |
|---|---|---|---|
| LogP | RDKit | To be recorded | Crippen-Wildman `MolLogP`. Different from clogP / ALogP / XLogP |
| MolWt | RDKit | To be recorded | `Descriptors.MolWt` (IUPAC average atomic weight of all heavy atoms + implicit H) |
| NumRotatableBonds | RDKit | To be recorded | `CalcNumRotatableBonds` strict SMARTS definition (excludes terminal bonds) |
| (additional descriptor) | (tool) | — | (definition) |

> **Note**: Some APIs may not exist in older RDKit versions.
> If an alternative was used, document it here (e.g., `CalcNumAromaticAtoms` absent → replaced with `pyrun` batch).

---

## Data

- **Source**: Dataset name and provider
- **URL**: Download URL (or note if unavailable due to license restrictions)
- **License**: License name and conditions
- **Cache path**: `data/benchmark/<name>.csv` (auto-downloaded on first run)
- **Count**: N molecules (explain any discrepancy with the paper's count)
- **Data hash**: See `dataset_sha256` field in `result/runs/<ts>/lock_snapshot.json`

---

## Script

```
repro/rp<XX>_<name>/rp<XX>_<name>.m
```

**How to run**: Open MATLAB with the project root as CWD and run with Ctrl+Enter (section by section).

| Section | Content |
|---|---|
| Section 0 | Setup and environment capture (`emk.setup.snapshot()`) |
| Section 1 | Dataset loading |
| Section 2 | Preprocessing and descriptor/feature computation |
| Section 3 | Model training |
| Section 4 | Evaluation (CV / test set) |
| Section 5 | RF03 verification (`emk.repro.verify()`) |
| Section 6 | Save results (`makeRunDir()` → `emk.setup.lockfile()`) |

---

## Result (First run: YYYY-MM-DD)

| Metric | Value | Pass/Fail |
|---|---|---|
| (metric name) | (value) | ✅ PASS / ❌ FAIL |

**Environment (from lock_snapshot.json):**

| Item | Value |
|---|---|
| MATLAB | (version at runtime) |
| Python | (version at runtime) |
| RDKit | (version at runtime) |
| Commit | (commit hash at runtime) |

---

## Verification (RF03 Numerical Acceptance Criteria)

> RF03 category: **Cat A (absolute threshold)**. Add Cat B / Cat C subsections when applicable.
> RF04: mark as compliant only when RF01 / RF02 / RF03 Cat A are satisfied.

### Cat A — Absolute Thresholds (required)

| Metric | Criterion | Rationale |
|---|---|---|
| (metric name) | (upper or lower bound) | (rationale: paper value and tolerance justification) |

**Tolerance rationale**:
1. (Descriptor implementation differences, etc.)
2. (Dataset version differences, etc.)

> **Note**: If the threshold was adjusted based on the first-run results, record the reasoning here.

<!-- Add this only when a two-model comparison is needed -->
<!-- ### Cat B — Relative Comparison (reference only) -->
<!-- | Metric | Model A | Model B | delta | t(df) | p (one-sided) | -->

<!-- Add this when an implementation sanity check is required -->
<!-- ### Cat C — Implementation Sanity (required when applicable) -->

---

## Discussion

### Differences from the Paper

| Difference | Details |
|---|---|
| (type of difference) | (details and impact) |

### Lessons Learned (handoff to subsequent RP / M-REPRO-FOUND)

- [ ] (lesson / open issue)

---

## Related Files

| File | Content |
|---|---|
| `rp<XX>_<name>.m` | Reproduction script |
| `lock_template.json` | RF02 version lock schema template |
| `result/runs/<ts>/lock_snapshot.json` | Actual version information recorded at runtime |
| `result/runs/<ts>/metrics.json` | Evaluation metrics |
