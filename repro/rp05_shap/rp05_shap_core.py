"""rp05_shap_core.py  RP05 SHAP Computation Core

Fit sklearn LogisticRegression on BBBP ECFP4 fingerprints and compute
feature-level SHAP values via shap.LinearExplainer.

Called from rp05_shap.m via pyrun:
    exec(open(helper_path).read())
    result_json = run_rp05(csv_path)

Returns a JSON string with metrics, global importance, and per-molecule
SHAP for 3 representative examples (TP / TN / MIS or TN2).
"""

import numpy as np
import pandas as pd
import json
import warnings
import shap
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator as rfg
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.exceptions import ConvergenceWarning
from scipy.stats import spearmanr


def run_rp05(csv_path, radius=2, n_bits=2048, C=1.0,
             n_fold=5, seed=42, n_top=20, test_size=0.2):
    """Fit LR on BBBP ECFP4, compute SHAP values, return JSON string."""

    # --- Load and validate data ---
    df = pd.read_csv(csv_path)
    # Cast smiles to str uniformly (NaN → "nan", stray float → "1.0", etc.)
    # so MolFromSmiles is the sole validity gate — it returns None for any
    # non-parseable input including "nan", making a separate notna() guard redundant.
    df["smiles"] = df["smiles"].astype(str)
    all_mols = [Chem.MolFromSmiles(s) for s in df["smiles"]]
    valid_mask = np.array([m is not None for m in all_mols])
    df_v = df[valid_mask].reset_index(drop=True)
    n_valid = len(df_v)

    # --- Compute Morgan ECFP4 fingerprints ---
    gen = rfg.GetMorganGenerator(radius=radius, fpSize=n_bits)
    mols = [m for m in all_mols if m is not None]   # reuse from validation pass
    X = np.array([gen.GetFingerprintAsNumPy(m) for m in mols], dtype="float64")
    y = df_v["p_np"].values.astype(int)

    # --- 5-fold stratified CV ROC-AUC ---
    # Pass an unfitted estimator so intent matches behavior (cross_val_score
    # clones and refits internally; passing a fitted model is misleading).
    estimator = LogisticRegression(C=C, penalty="l2", solver="lbfgs", max_iter=1000)
    skf = StratifiedKFold(n_splits=n_fold, shuffle=True, random_state=seed)
    cv_aucs = cross_val_score(estimator, X, y, cv=skf, scoring="roc_auc")
    auc_cv  = float(cv_aucs.mean())
    auc_std = float(cv_aucs.std())

    # --- Train/test split for SHAP (avoids background data leakage) ---
    # LinearExplainer baseline E[f(X)] must be computed from train data only;
    # using the full dataset leaks test-set distribution into the SHAP baseline.
    idx_all = np.arange(n_valid)
    X_tr, X_te, y_tr, y_te, idx_tr, idx_te = train_test_split(
        X, y, idx_all, test_size=test_size, stratify=y, random_state=seed)
    n_train = len(y_tr)
    n_test  = len(y_te)

    # --- Fit final model on train set only ---
    model = LogisticRegression(C=C, penalty="l2", solver="lbfgs", max_iter=1000)
    model.fit(X_tr, y_tr)

    # --- SHAP: LinearExplainer (exact for linear models) ---
    # Background = X_tr; shap_values evaluated on X_te.
    # NOTE: CV uses full X (n_valid) while SHAP background is X_tr (80% split).
    # These are intentionally asymmetric: CV provides an unbiased AUC estimate
    # while SHAP avoids leaking test-set distribution into the LinearExplainer
    # baseline E[f(X)]. They measure different things and are not directly comparable.
    explainer = shap.LinearExplainer(model, X_tr)
    sv = explainer.shap_values(X_te)
    if isinstance(sv, shap.Explanation):
        sv = sv.values
        if sv.ndim == 3:   # (n_samples, n_features, n_classes)
            sv = sv[:, :, 1]
    elif isinstance(sv, list):
        sv = sv[1]   # positive class for binary output

    global_imp = np.mean(np.abs(sv), axis=0)   # (n_bits,)
    coefs      = model.coef_[0]                 # (n_bits,)
    std_x      = X_tr.std(axis=0)              # train-set std (n_bits,)

    # Spearman rho: rank(global_imp) vs rank(|coef| * std_train(X)).
    # For binary X_j in {0,1} with feature frequency p_j = mean(X[:,j]):
    #   E[|X_j - p_j|] = 2*p_j*(1-p_j),  std(X_j) = sqrt(p_j*(1-p_j))
    # Both are monotone functions of p_j, so ranks agree even though the
    # proportionality constant is not exactly std(X_j).
    lr_imp = np.abs(coefs) * std_x
    # Exclude bits with std_x == 0 (constant across all train molecules) to
    # avoid Spearman ρ distortion when lr_imp == 0 but global_imp != 0.
    # These bits are also absent from top_n_idx: LinearExplainer yields
    # SHAP[i,j] = coef[j] * (X[i,j] - E[X[:,j]]) = 0 when std_x[j] == 0,
    # so global_imp[j] == 0 and they cannot rank in the top-N.
    nonzero_std = std_x != 0
    if nonzero_std.sum() >= 2:
        rho, _ = spearmanr(global_imp[nonzero_std], lr_imp[nonzero_std])
    else:
        rho = float("nan")
    shap_lr_spearman = float(rho)

    # --- Select 3 representative molecules from test set ---
    y_pred = model.predict(X_te)
    y_prob  = model.predict_proba(X_te)[:, 1]

    # Best TP (high confidence BBB+)
    tp_scores = np.where((y_te == 1) & (y_pred == 1), y_prob, -np.inf)
    tp_idx    = int(np.argmax(tp_scores))

    # Best TN (high confidence BBB-)
    tn_scores = np.where((y_te == 0) & (y_pred == 0), 1.0 - y_prob, -np.inf)
    tn_idx    = int(np.argmax(tn_scores))

    # Most confidently misclassified; fall back to 2nd-best TN when test set
    # is fully correct to prevent TP appearing under the "MIS" label.
    mis_mask = y_te != y_pred
    if mis_mask.any():
        # |prob - 0.5| uniformly measures prediction confidence regardless of
        # direction, so the most confidently wrong sample is always selected.
        mis_conf = np.abs(y_prob - 0.5)
        mis_conf_masked = np.where(mis_mask, mis_conf, -np.inf)
        mis_idx  = int(np.argmax(mis_conf_masked))
        ex_types = ["TP", "TN", "MIS"]
    else:
        tn_scores2 = tn_scores.copy()
        tn_scores2[tn_idx] = -np.inf          # exclude 1st TN
        mis_idx  = int(np.argmax(tn_scores2))
        ex_types = ["TP", "TN", "TN2"]
        warnings.warn("No misclassified samples in test set; MIS replaced by TN2.")

    # Local indices into X_te / sv / y_pred / y_prob
    local_ex_idx = [tp_idx, tn_idx, mis_idx]
    # Global indices into df_v for name / smiles lookup
    global_ex_idx = [int(idx_te[i]) for i in local_ex_idx]

    # Top-n bits by mean |SHAP| for bar chart axes
    top_n_idx = np.argsort(global_imp)[::-1][:n_top].tolist()

    result = {
        # Metrics
        "auc_cv":            auc_cv,
        "auc_cv_std":        auc_std,
        "n_valid":           n_valid,
        "n_train":           n_train,
        "n_test":            n_test,
        "n_bbb_pos":         int(y.sum()),
        "n_bbb_pos_train":   int(y_tr.sum()),
        "n_bbb_pos_test":    int(y_te.sum()),
        "shap_lr_spearman":  shap_lr_spearman,
        # Global importance (top-n bits)
        "top_n_idx":   top_n_idx,               # 0-based bit indices
        "top_n_imp":   [float(global_imp[i]) for i in top_n_idx],
        "top_n_coef":  [float(coefs[i])      for i in top_n_idx],
        # Full 2048-bit importance array (for ranking comparison in task B)
        "global_imp_all": global_imp.tolist(),
        # 3 example molecules
        "ex_idx":    global_ex_idx,
        "ex_types":  ex_types,
        "ex_labels": [int(y_te[i])      for i in local_ex_idx],
        "ex_preds":  [int(y_pred[i]) for i in local_ex_idx],
        "ex_probs":  [float(y_prob[i]) for i in local_ex_idx],
        "ex_smiles": [str(df_v["smiles"].iloc[global_ex_idx[k]]) for k in range(3)],
        "ex_names":  [str(df_v["name"].iloc[global_ex_idx[k]])   for k in range(3)],
        # Per-example SHAP for top-n bits (n_ex x n_top)
        "ex_top_shap": [[float(sv[i, j]) for j in top_n_idx] for i in local_ex_idx],
        "ex_top_xval": [[float(X_te[i,  j]) for j in top_n_idx] for i in local_ex_idx],
    }

    return json.dumps(result)


# ---------------------------------------------------------------------------
# Task B (M-REPRO-REFINE Phase 2): SHAP on RP02-rev model
# ---------------------------------------------------------------------------

_C_GRID_REV = [0.01, 0.05, 0.1, 0.5, 1.0, 2.0, 5.0]


def run_rp05_rev(csv_path, radius=2, n_bits=2048,
                 n_outer=5, n_inner=3,
                 outer_seed=42, inner_seed=7,
                 n_top=20, test_size=0.2, c_grid=None):
    """SHAP on RP02-rev model: C optimized via nested CV.

    AUC: outer 5-fold stratified CV x inner 3-fold C selection (RP02-rev design,
         outer_seed=42, inner_seed=7 matching current rp02_sklearn_core.py).

    SHAP model: inner 3-fold CV on X_tr (the same 80% split as run_rp05,
         random_state=outer_seed=42) selects best C; final model fit on X_tr.

    Called from rp05_shap.m via:
        exec(open(hp).read(), globals()); rev_json = run_rp05_rev(cp)
    """
    if c_grid is None:
        c_grid = _C_GRID_REV

    # --- Load and validate data ---
    df = pd.read_csv(csv_path)
    df["smiles"] = df["smiles"].astype(str)
    all_mols = [Chem.MolFromSmiles(s) for s in df["smiles"]]
    valid_mask = np.array([m is not None for m in all_mols])
    df_v = df[valid_mask].reset_index(drop=True)
    n_valid = len(df_v)

    # --- Compute ECFP4 fingerprints ---
    gen = rfg.GetMorganGenerator(radius=radius, fpSize=n_bits)
    mols = [m for m in all_mols if m is not None]
    X = np.array([gen.GetFingerprintAsNumPy(m) for m in mols], dtype="float64")
    y = df_v["p_np"].values.astype(int)

    # --- AUC: nested CV (RP02-rev design) ---
    outer_cv = StratifiedKFold(n_splits=n_outer, shuffle=True, random_state=outer_seed)
    inner_cv = StratifiedKFold(n_splits=n_inner, shuffle=True, random_state=inner_seed)

    cv_aucs = []
    best_cs_cv = []
    for tr, te in outer_cv.split(X, y):
        best_c, best_inner = None, -1.0
        for C in c_grid:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=ConvergenceWarning)
                aucs_c = cross_val_score(
                    LogisticRegression(C=C, penalty="l2", solver="lbfgs", max_iter=2000),
                    X[tr], y[tr], cv=inner_cv, scoring="roc_auc")
            mu = float(aucs_c.mean())
            if mu > best_inner:
                best_inner, best_c = mu, C
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ConvergenceWarning)
            m_fold = LogisticRegression(C=best_c, penalty="l2", solver="lbfgs",
                                        max_iter=2000)
            m_fold.fit(X[tr], y[tr])
        auc_fold = float(roc_auc_score(y[te], m_fold.predict_proba(X[te])[:, 1]))
        cv_aucs.append(auc_fold)
        best_cs_cv.append(float(best_c))

    auc_cv  = float(np.mean(cv_aucs))
    auc_std = float(np.std(cv_aucs))

    # --- 80/20 split for SHAP (same seed as run_rp05 for direct comparison) ---
    idx_all = np.arange(n_valid)
    X_tr, X_te, y_tr, y_te, idx_tr, idx_te = train_test_split(
        X, y, idx_all, test_size=test_size, stratify=y, random_state=outer_seed)

    # Inner CV on X_tr to select best C for the SHAP model
    best_c_shap, best_inner_shap = None, -1.0
    for C in c_grid:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ConvergenceWarning)
            aucs_c = cross_val_score(
                LogisticRegression(C=C, penalty="l2", solver="lbfgs", max_iter=2000),
                X_tr, y_tr, cv=inner_cv, scoring="roc_auc")
        mu = float(aucs_c.mean())
        if mu > best_inner_shap:
            best_inner_shap, best_c_shap = mu, C

    # Fit final SHAP model on X_tr with best C
    model = LogisticRegression(C=best_c_shap, penalty="l2", solver="lbfgs",
                                max_iter=2000)
    model.fit(X_tr, y_tr)

    # --- SHAP: LinearExplainer (exact for linear models) ---
    explainer = shap.LinearExplainer(model, X_tr)
    sv = explainer.shap_values(X_te)
    if isinstance(sv, shap.Explanation):
        sv = sv.values
        if sv.ndim == 3:
            sv = sv[:, :, 1]
    elif isinstance(sv, list):
        sv = sv[1]

    global_imp = np.mean(np.abs(sv), axis=0)
    coefs      = model.coef_[0]
    std_x      = X_tr.std(axis=0)

    lr_imp = np.abs(coefs) * std_x
    nonzero_std = std_x != 0
    if nonzero_std.sum() >= 2:
        rho, _ = spearmanr(global_imp[nonzero_std], lr_imp[nonzero_std])
    else:
        rho = float("nan")

    top_n_idx = np.argsort(global_imp)[::-1][:n_top].tolist()

    return json.dumps({
        "auc_cv":           auc_cv,
        "auc_cv_std":       auc_std,
        "auc_per_fold":     cv_aucs,
        "best_C_per_fold":  best_cs_cv,
        "best_C_shap":      float(best_c_shap),
        "n_valid":          n_valid,
        "n_train":          int(len(y_tr)),
        "n_test":           int(len(y_te)),
        "n_bbb_pos":        int(y.sum()),
        "n_bbb_pos_train":  int(y_tr.sum()),
        "n_bbb_pos_test":   int(y_te.sum()),
        "shap_lr_spearman": float(rho),
        "top_n_idx":        top_n_idx,
        "top_n_imp":        [float(global_imp[i]) for i in top_n_idx],
        "top_n_coef":       [float(coefs[i])      for i in top_n_idx],
        "global_imp_all":   global_imp.tolist(),
    })


# ---------------------------------------------------------------------------
# F2 (M-REPRO-REFINE Phase 5): Data for MATLAB Analytical SHAP Comparison
# ---------------------------------------------------------------------------

def run_rp05_f2(csv_path, mat_path, radius=2, n_bits=2048, C_shap=0.10,
                test_size=0.2, outer_seed=42):
    """F2: Compute Python SHAP and save data for MATLAB analytical SHAP comparison.

    Fits LR with C=C_shap (=0.10, matching best_c_shap from run_rp05_rev inner CV).
    Saves X_tr, X_te, y_tr, y_te, and Python LinearExplainer global importance
    (evaluated on X_te) to mat_path via scipy.io.savemat.

    MATLAB loads mat_path, fits fitclinear (lbfgs, Lambda=1/n), computes analytical
    SHAP on X_te (same formula as shapley(..., Method='interventional-linear')), and
    compares feature importance rankings against Python output (Spearman rho >= 0.85
    confirms Zone B).

    Background for both Python and MATLAB SHAP: X_tr (80% train split).
    Query set for global importance comparison: X_te (20% held-out test set).
    Full 2039-mol SHAP (background=X_tr, query=X_all) computed on MATLAB side only,
    as an efficiency check (analytical formula is O(n*p) = fast).

    Same 80/20 split as run_rp05 and run_rp05_rev (random_state=outer_seed=42).
    """
    try:
        import scipy.io

        # Load + ECFP4 (identical pipeline to run_rp05 / run_rp05_rev)
        df = pd.read_csv(csv_path)
        df["smiles"] = df["smiles"].astype(str)
        all_mols = [Chem.MolFromSmiles(s) for s in df["smiles"]]
        valid_mask = np.array([m is not None for m in all_mols])
        df_v = df[valid_mask].reset_index(drop=True)
        gen = rfg.GetMorganGenerator(radius=radius, fpSize=n_bits)
        mols = [m for m in all_mols if m is not None]
        X = np.array([gen.GetFingerprintAsNumPy(m) for m in mols], dtype="float64")
        y = df_v["p_np"].values.astype(int)

        # Same 80/20 split as run_rp05 and run_rp05_rev
        X_tr, X_te, y_tr, y_te = train_test_split(
            X, y, test_size=test_size, stratify=y, random_state=outer_seed)

        # Fit Python LR (C=C_shap, matches best_c_shap from run_rp05_rev)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=ConvergenceWarning)
            model = LogisticRegression(C=C_shap, penalty="l2", solver="lbfgs",
                                       max_iter=2000, random_state=outer_seed)
            model.fit(X_tr, y_tr)

        # Python LinearExplainer SHAP on X_te (background = X_tr)
        explainer = shap.LinearExplainer(model, X_tr)
        sv_te = explainer.shap_values(X_te)
        if isinstance(sv_te, shap.Explanation):
            sv_te = sv_te.values
            if sv_te.ndim == 3:
                sv_te = sv_te[:, :, 1]
        elif isinstance(sv_te, list):
            sv_te = sv_te[1]
        global_imp_te = np.mean(np.abs(sv_te), axis=0)   # (n_bits,)

        # Save to .mat file for MATLAB load()
        scipy.io.savemat(str(mat_path), {
            "X_tr":                 X_tr,
            "X_te":                 X_te,
            "y_tr":                 y_tr.reshape(-1, 1).astype(float),
            "y_te":                 y_te.reshape(-1, 1).astype(float),
            "global_imp_python_te": global_imp_te.reshape(-1, 1),
            "coef_python":          model.coef_[0].reshape(-1, 1),
            "n_train":              float(len(y_tr)),
            "n_test":               float(len(y_te)),
            "C_shap":               float(C_shap),
        })

        return json.dumps({
            "success":  True,
            "n_train":  int(len(y_tr)),
            "n_test":   int(len(y_te)),
            "C_shap":   float(C_shap),
            "mat_path": str(mat_path),
        })

    except Exception as exc:
        import traceback as tb
        return json.dumps({
            "success":   False,
            "error":     str(exc),
            "traceback": tb.format_exc(),
        })


# ---------------------------------------------------------------------------
# F3 (M-REPRO-REFINE Phase 5): RF TreeSHAP for MATLAB shapley() comparison
# ---------------------------------------------------------------------------

def run_rp05_f3(csv_path, mat_path, radius=2, n_bits=2048,
                n_estimators=100, test_size=0.2, seed=42, n_fold=5,
                n_background=200):
    """F3: sklearn RF + TreeSHAP reference data for MATLAB shapley() sampling study.

    Trains RandomForestClassifier (n_estimators=100, random_state=42) on BBBP ECFP4.
    Computes Python reference SHAP values via shap.TreeExplainer using a 200-row
    background subset of X_tr. The reference global importance is evaluated on the
    full held-out test set. MATLAB then reuses the same X_tr / X_te / X_bg and
    compares subsampled shapley(...).fit(X_eval) rankings against this reference.

    Saved MAT payload intentionally contains raw matrices plus the Python full-test
    importance vector so MATLAB can run repeated n_eval-grid experiments without
    recomputing Python TreeSHAP.

    Same 80/20 split as run_rp05 / run_rp05_f2 (random_state=seed=42).

    Called from rp05_shap.m via:
        exec(open(hp).read(), globals()); f3_json = run_rp05_f3(cp, mp)
    """
    try:
        import scipy.io
        from sklearn.ensemble import RandomForestClassifier

        # Load + ECFP4 (identical pipeline to existing functions)
        df = pd.read_csv(csv_path)
        df["smiles"] = df["smiles"].astype(str)
        all_mols = [Chem.MolFromSmiles(s) for s in df["smiles"]]
        valid_mask = np.array([m is not None for m in all_mols])
        df_v = df[valid_mask].reset_index(drop=True)
        gen = rfg.GetMorganGenerator(radius=radius, fpSize=n_bits)
        mols = [m for m in all_mols if m is not None]
        X = np.array([gen.GetFingerprintAsNumPy(m) for m in mols], dtype="float64")
        y = df_v["p_np"].values.astype(int)

        # Same 80/20 split as run_rp05 and run_rp05_f2
        X_tr, X_te, y_tr, y_te = train_test_split(
            X, y, test_size=test_size, stratify=y, random_state=seed)

        # 5-fold CV AUC (informational; uses full dataset like run_rp05)
        rf_cv = RandomForestClassifier(n_estimators=n_estimators, random_state=seed)
        skf = StratifiedKFold(n_splits=n_fold, shuffle=True, random_state=seed)
        cv_aucs = cross_val_score(rf_cv, X, y, cv=skf, scoring="roc_auc")

        # Final model on training set
        rf = RandomForestClassifier(n_estimators=n_estimators, random_state=seed)
        rf.fit(X_tr, y_tr)

        # Subsample background to n_background (default 200) for speed.
        # Full X_tr (~1211 samples) gives ~6x slower computation with no
        # meaningful ranking change for binary ECFP4 features.
        # Same seed ensures Python and MATLAB use identical background rows.
        bg_rng = np.random.RandomState(seed)
        bg_idx = bg_rng.choice(len(X_tr), size=n_background, replace=False)
        X_bg = X_tr[bg_idx]

        # TreeSHAP: exact for tree ensembles; interventional with background=X_bg
        # avoids conditioning on training distribution leaking into test SHAP values.
        explainer = shap.TreeExplainer(rf, data=X_bg,
                                       feature_perturbation="interventional")
        sv = explainer.shap_values(X_te)
        # sklearn RF returns list of 2 arrays (class 0, class 1) for binary
        if isinstance(sv, list):
            sv_pos = sv[1]          # class 1 = BBB+
        elif hasattr(sv, "ndim") and sv.ndim == 3:
            sv_pos = sv[:, :, 1]   # (n_samples, n_features, n_classes)[:, :, 1]
        else:
            sv_pos = sv

        global_imp = np.mean(np.abs(sv_pos), axis=0)   # (n_bits,)

        scipy.io.savemat(str(mat_path), {
            "X_tr":                     X_tr,
            "X_te":                     X_te,
            "X_bg":                     X_bg,
            "y_tr":                     y_tr.reshape(-1, 1).astype(float),
            "y_te":                     y_te.reshape(-1, 1).astype(float),
            "global_imp_treeshap_full": global_imp.reshape(-1, 1),
            "bg_idx":                   bg_idx.reshape(-1, 1).astype(float),
            "rf_auc_cv":                float(cv_aucs.mean()),
            "rf_auc_std":               float(cv_aucs.std()),
            "n_background":             float(n_background),
            "seed":                     float(seed),
            "n_estimators":             float(n_estimators),
        })

        return json.dumps({
            "success":      True,
            "rf_auc_cv":    float(cv_aucs.mean()),
            "rf_auc_std":   float(cv_aucs.std()),
            "n_train":      int(len(y_tr)),
            "n_test":       int(len(y_te)),
            "n_background": int(n_background),
            "n_estimators": int(n_estimators),
        })

    except Exception as exc:
        import traceback as tb
        return json.dumps({
            "success":   False,
            "error":     str(exc),
            "traceback": tb.format_exc(),
        })
