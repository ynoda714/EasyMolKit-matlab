# RP<XX>: <Title> — Standard Reproduction Template (RF01)

> **Purpose**: RF01 standard reproduction template. Each reproduction study must follow this structure.
> Established in M-REPRO-FOUND (2026-06-19) based on experience from RP00.

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

## Verification (RF03 Acceptance Criteria)

| Metric | Criterion | Rationale |
|---|---|---|
| (metric name) | (upper or lower bound) | (rationale: paper value and tolerance justification) |

**Tolerance rationale**:
1. (Descriptor implementation differences, etc.)
2. (Dataset version differences, etc.)

> **Note**: If the threshold was adjusted based on the first-run results, record the reasoning here.

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
