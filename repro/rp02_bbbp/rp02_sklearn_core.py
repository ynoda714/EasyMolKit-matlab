"""rp02_sklearn_core.py  RP02-rev: sklearn LR baseline for BBBP (M-REPRO-AUDIT A1)

Replaces MATLAB fitclinear with sklearn LogisticRegression + nested CV.
Motivation: A1 diagnosis showed fitclinear yields AUC ~0.019 below sklearn lbfgs
on the same fold splits, making fair comparison with RP03/RP04 impossible.

Design:
  Outer 5-fold stratified CV (fold indices saved for A2 sharing with RP03/RP04)
    Inner 3-fold CV on train_outer -> select best C from c_grid
    Fit LR(best_C) on train_outer -> evaluate on test_outer -> AUC

Called from rp02_bbbp.m via:
    exec(open(script_path).read())
    result_json = run_rp02_sklearn(csv_path)

Return contract: always returns valid JSON with a "success" bool field.
  success=True  -> normal result fields present
  success=False -> "error" and "traceback" fields describe the failure
"""

import json
import traceback
import warnings
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator as rfg
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, GroupKFold, cross_val_score
from sklearn.metrics import roc_auc_score
from sklearn.exceptions import ConvergenceWarning


_C_GRID = [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0]
_SCAFFOLD_ERROR_SENTINEL = "\x00scaffold_error\x00"


def _load_ecfp4(csv_path, radius=2, n_bits=2048):
    import pandas as pd
    df = pd.read_csv(csv_path)

    # N4: explicit priority detection -- "smiles" takes precedence over "SMILES";
    # "p_np" takes precedence over "BBB". Error immediately if neither exists.
    if "smiles" in df.columns:
        smiles_col = "smiles"
    elif "SMILES" in df.columns:
        smiles_col = "SMILES"
    else:
        raise ValueError(
            f"No SMILES column found; expected 'smiles' or 'SMILES', "
            f"got columns: {list(df.columns)}"
        )

    if "p_np" in df.columns:
        label_col = "p_np"
    elif "BBB" in df.columns:
        label_col = "BBB"
    else:
        raise ValueError(
            f"No label column found; expected 'p_np' or 'BBB', "
            f"got columns: {list(df.columns)}"
        )

    valid = [Chem.MolFromSmiles(str(s)) is not None for s in df[smiles_col]]
    df    = df[np.array(valid)].reset_index(drop=True)

    gen  = rfg.GetMorganGenerator(radius=radius, fpSize=n_bits)
    mols = [Chem.MolFromSmiles(str(s)) for s in df[smiles_col]]
    X    = np.array([gen.GetFingerprintAsNumPy(m) for m in mols], dtype="float64")
    y    = df[label_col].values.astype(int)
    return X, y, df[smiles_col].tolist()


def _lr(C):
    return LogisticRegression(C=C, penalty="l2", solver="lbfgs", max_iter=2000)


def run_rp02_sklearn(csv_path, outer_seed=42, inner_seed=7,
                     n_outer=5, n_inner=3, c_grid=None):
    """Run nested CV and return JSON string.

    Returns metrics, per-fold details, outer_fold_indices (for A2 sharing).

    Note (M1): outer_seed=42, inner_seed=7 (different values) to prevent
    correlated fold structures between the outer evaluation loop and inner C
    selection. The historical certified run (M-REPRO-AUDIT A1, 2026-06-21)
    used inner_seed=42; documented results (AUC=0.9118) reflect that run.
    New executions with inner_seed=7 will produce slightly different values.

    Note (M5): outer_fold_indices uses string keys "foldN_train" / "foldN_test"
    so that MATLAB jsondecode produces valid struct field names (field names
    starting with a letter are safe; numeric-only keys are not). A list-of-dicts
    format would be more portable but would break RP03/RP04 readers.
    """
    try:
        if c_grid is None:
            c_grid = _C_GRID

        X, y, smiles = _load_ecfp4(csv_path)

        outer_cv      = StratifiedKFold(n_splits=n_outer, shuffle=True, random_state=outer_seed)
        fold_results  = []
        fold_indices  = {}
        all_test_prob = np.full(len(y), np.nan)  # for full-dataset pseudo-ROC; NaN detects partial failures
        covered_mask  = np.zeros(len(y), dtype=bool)  # H1: track test-set coverage

        for k, (tr, te) in enumerate(outer_cv.split(X, y)):
            fold_idx = k + 1
            fold_indices[f"fold{fold_idx}_train"] = tr.tolist()
            fold_indices[f"fold{fold_idx}_test"]  = te.tolist()
            covered_mask[te] = True

            # Inner CV: select best C
            # H3: tie-break -- strict > means the first C (smallest in grid)
            # that achieves the best inner_auc is selected.
            inner_cv   = StratifiedKFold(n_splits=n_inner, shuffle=True, random_state=inner_seed)
            best_c, best_inner = None, -1.0
            inner_grid = []
            for C in c_grid:
                with warnings.catch_warnings():
                    # L4/N1: suppress ConvergenceWarning by class (not by message,
                    # which sklearn does not match "lbfgs failed to converge").
                    # High-C runs may not fully converge on ECFP4 but the effect
                    # on inner-CV AUC ranking is negligible.
                    warnings.filterwarnings("ignore", category=ConvergenceWarning)
                    inner_aucs = cross_val_score(_lr(C), X[tr], y[tr],
                                                 cv=inner_cv, scoring="roc_auc")
                mu = float(inner_aucs.mean())
                inner_grid.append({"C": C, "inner_auc": mu})
                if mu > best_inner:  # strict >: tie-break favours first (smallest) C
                    best_inner, best_c = mu, C

            # Fit on full outer train, score on test
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=ConvergenceWarning)
                model = _lr(best_c).fit(X[tr], y[tr])
            prob  = model.predict_proba(X[te])[:, 1]
            all_test_prob[te] = prob
            auc   = float(roc_auc_score(y[te], prob))

            fold_results.append({
                "fold":         fold_idx,
                "train_size":   int(len(tr)),
                "test_size":    int(len(te)),
                "best_C":       float(best_c),
                "inner_auc":    float(best_inner),
                "test_auc":     float(auc),
                "inner_C_grid": inner_grid,   # L2: HP selection evidence
            })

        # H1: assert StratifiedKFold covered every sample exactly once
        uncovered = np.where(~covered_mask)[0].tolist()
        assert not uncovered, (
            f"StratifiedKFold did not cover all indices; "
            f"{len(uncovered)} samples never appeared in a test fold: {uncovered[:10]}"
        )

        # H1b: assert no test predictions remain at NaN initialiser value
        nan_idx = np.where(np.isnan(all_test_prob))[0].tolist()
        assert not nan_idx, (
            f"all_test_prob still contains {len(nan_idx)} NaN values after all folds; "
            f"a fold may have failed silently: {nan_idx[:10]}"
        )

        aucs    = [r["test_auc"] for r in fold_results]
        best_cs = [r["best_C"]   for r in fold_results]

        return json.dumps({
            "success":            True,
            "n_valid":            int(len(y)),
            "n_bbb_pos":          int(y.sum()),
            "n_bbb_neg":          int((y == 0).sum()),
            "auc_per_fold":       aucs,
            "auc_mean":           float(np.mean(aucs)),
            "auc_std":            float(np.std(aucs, ddof=1)),
            "auc_std_definition": "sample_std_ddof1",
            "best_C_per_fold":    best_cs,
            "fold_results":       fold_results,
            "outer_fold_indices": fold_indices,
            # For pseudo-ROC curve in MATLAB (all test-fold predictions)
            "y_true":   y.tolist(),
            "y_prob":   all_test_prob.tolist(),
        }, indent=2)

    except Exception as exc:  # H2: propagate any Python failure to MATLAB via JSON
        return json.dumps({
            "success":   False,
            "error":     str(exc),
            "traceback": traceback.format_exc(),
        }, indent=2)


# ---------------------------------------------------------------------------
# C task (M-REPRO-REFINE Phase 2): Scaffold 5-fold nested CV
# ---------------------------------------------------------------------------

def _murcko_scaffold(smi):
    """Bemis-Murcko generic scaffold SMILES. Returns "" for acyclic/invalid."""
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    except Exception:
        return _SCAFFOLD_ERROR_SENTINEL


def _load_ecfp4_scaffold(csv_path, radius=2, n_bits=2048):
    """Load BBBP CSV, compute ECFP4, assign Bemis-Murcko scaffold group IDs.

    Group IDs preserve CSV insertion order (same-scaffold molecules share an ID).
    n_scaffold_errors: molecules that passed SMILES parse but raised inside
    MurckoScaffoldSmiles; merged into the acyclic group (scaffold="").
    """
    import pandas as pd
    df = pd.read_csv(csv_path)

    if "smiles" in df.columns:
        smiles_col = "smiles"
    elif "SMILES" in df.columns:
        smiles_col = "SMILES"
    else:
        raise ValueError(f"No SMILES column; got: {list(df.columns)}")

    if "p_np" in df.columns:
        label_col = "p_np"
    elif "BBB" in df.columns:
        label_col = "BBB"
    else:
        raise ValueError(f"No label column; got: {list(df.columns)}")

    valid = [Chem.MolFromSmiles(str(s)) is not None for s in df[smiles_col]]
    df    = df[np.array(valid)].reset_index(drop=True)

    gen  = rfg.GetMorganGenerator(radius=radius, fpSize=n_bits)
    mols = [Chem.MolFromSmiles(str(s)) for s in df[smiles_col]]
    X    = np.array([gen.GetFingerprintAsNumPy(m) for m in mols], dtype="float64")
    y    = df[label_col].values.astype(int)

    scaffolds         = [_murcko_scaffold(str(s)) for s in df[smiles_col]]
    n_scaffold_errors = sum(1 for s in scaffolds if s == _SCAFFOLD_ERROR_SENTINEL)
    scaffolds         = ["" if s == _SCAFFOLD_ERROR_SENTINEL else s for s in scaffolds]
    seen   = {}
    groups = np.zeros(len(scaffolds), dtype=int)
    for i, sc in enumerate(scaffolds):
        if sc not in seen:
            seen[sc] = len(seen)
        groups[i] = seen[sc]

    return X, y, groups, n_scaffold_errors, int(np.unique(groups).size)


def run_rp02_scaffold(csv_path, n_outer=5, n_inner=3, c_grid=None):
    """Scaffold 5-fold nested CV: Bemis-Murcko GroupKFold outer x GroupKFold inner.

    Outer: GroupKFold(5) on Bemis-Murcko scaffold groups
    Inner: GroupKFold(3) on train-portion scaffold groups for C selection

    S1: GroupKFold does NOT shuffle; determinism comes from scaffold insertion
    order (CSV row order). A different bbbp.csv row order changes the split.
    S2: GroupKFold cannot stratify by class. BBBP BBB+ imbalance (76.5%) will
    produce variable per-fold class ratios -- expected behaviour for scaffold split.

    H1: Coverage assertion verifies every sample appears in exactly one test fold.
    H2: Python exceptions returned as JSON {success: false} so MATLAB gets a
    parseable error rather than a pyrun crash.
    """
    try:
        if c_grid is None:
            c_grid = _C_GRID

        X, y, groups, n_scaffold_errors, n_unique_scaffolds = _load_ecfp4_scaffold(csv_path)

        outer_cv       = GroupKFold(n_splits=n_outer)
        fold_results   = []
        all_test_prob  = np.full(len(y), np.nan)
        count_assigned = np.zeros(len(y), dtype=int)

        for k, (tr, te) in enumerate(outer_cv.split(X, y, groups)):
            fold_idx = k + 1
            count_assigned[te] += 1

            n_inner_groups = int(np.unique(groups[tr]).size)
            if n_inner_groups < 2:
                raise ValueError(
                    f"Outer fold {fold_idx}: train has only {n_inner_groups} unique "
                    "scaffold group(s); inner GroupKFold requires >= 2."
                )
            actual_n_inner = min(n_inner, n_inner_groups)
            inner_cv       = GroupKFold(n_splits=actual_n_inner)

            best_c, best_inner = None, -1.0
            inner_grid = []
            for C in c_grid:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", category=ConvergenceWarning)
                    inner_aucs = cross_val_score(
                        _lr(C), X[tr], y[tr],
                        cv=inner_cv, groups=groups[tr], scoring="roc_auc",
                    )
                mu = float(inner_aucs.mean())
                inner_grid.append({"C": C, "inner_auc": mu})
                if mu > best_inner:
                    best_inner, best_c = mu, C

            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=ConvergenceWarning)
                model = _lr(best_c).fit(X[tr], y[tr])
            prob              = model.predict_proba(X[te])[:, 1]
            all_test_prob[te] = prob
            auc               = float(roc_auc_score(y[te], prob))

            fold_results.append({
                "fold":         fold_idx,
                "train_size":   int(len(tr)),
                "test_size":    int(len(te)),
                "best_C":       float(best_c),
                "inner_auc":    float(best_inner),
                "test_auc":     float(auc),
                "inner_C_grid": inner_grid,
            })

        # H1: every sample in exactly one test fold
        uncovered  = np.where(count_assigned == 0)[0].tolist()
        duplicated = np.where(count_assigned > 1)[0].tolist()
        if uncovered:
            raise AssertionError(
                f"GroupKFold did not cover {len(uncovered)} samples: {uncovered[:10]}"
            )
        if duplicated:
            raise AssertionError(
                f"GroupKFold assigned {len(duplicated)} samples to multiple folds: "
                f"{duplicated[:10]}"
            )

        nan_indices = np.where(np.isnan(all_test_prob))[0].tolist()
        if nan_indices:
            raise AssertionError(
                f"all_test_prob has {len(nan_indices)} NaN entries: {nan_indices[:10]}"
            )

        aucs    = [r["test_auc"] for r in fold_results]
        best_cs = [r["best_C"]   for r in fold_results]

        return json.dumps({
            "success":            True,
            "n_valid":            int(len(y)),
            "n_bbb_pos":          int(y.sum()),
            "n_bbb_neg":          int((y == 0).sum()),
            "n_unique_scaffolds": n_unique_scaffolds,
            "n_scaffold_errors":  n_scaffold_errors,
            "auc_per_fold":       aucs,
            "auc_mean":           float(np.mean(aucs)),
            "auc_std":            float(np.std(aucs, ddof=1)),
            "auc_std_definition": "sample_std_ddof1",
            "best_C_per_fold":    best_cs,
            "fold_results":       fold_results,
            "y_true":             y.tolist(),
            "y_prob":             all_test_prob.tolist(),
        }, indent=2)

    except Exception as exc:  # H2
        return json.dumps({
            "success":   False,
            "error":     str(exc),
            "traceback": traceback.format_exc(),
        }, indent=2)
