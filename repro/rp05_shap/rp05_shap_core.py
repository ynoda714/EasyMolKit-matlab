"""rp05_shap_core.py  RP05 SHAP Computation Core

Fit sklearn LogisticRegression on BBBP ECFP4 fingerprints and compute
feature-level SHAP values via shap.LinearExplainer.

Called from rp05_shap.m via pyrun:
    exec(open(helper_path).read())
    result_json = run_rp05(csv_path)

Returns a JSON string with metrics, global importance, and per-molecule
SHAP for 3 representative examples (TP / TN / misclassified).
"""

import numpy as np
import pandas as pd
import json
import shap
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator as rfg
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from scipy.stats import spearmanr


def run_rp05(csv_path, radius=2, n_bits=2048, C=1.0,
             n_fold=5, seed=42, n_top=20):
    """Fit LR on BBBP ECFP4, compute SHAP values, return JSON string."""

    # --- Load and validate data ---
    df = pd.read_csv(csv_path)
    valid_mask = np.array(
        [Chem.MolFromSmiles(str(s)) is not None for s in df["smiles"]]
    )
    df_v = df[valid_mask].reset_index(drop=True)
    n_valid = len(df_v)

    # --- Compute Morgan ECFP4 fingerprints ---
    gen = rfg.GetMorganGenerator(radius=radius, fpSize=n_bits)
    mols = [Chem.MolFromSmiles(str(s)) for s in df_v["smiles"]]
    X = np.array([gen.GetFingerprintAsNumPy(m) for m in mols], dtype="float64")
    y = df_v["p_np"].values.astype(int)

    # --- Fit LogisticRegression (L2 ridge, same regularization as RP02) ---
    model = LogisticRegression(C=C, penalty="l2", solver="lbfgs", max_iter=1000)
    model.fit(X, y)

    # --- 5-fold stratified CV ROC-AUC ---
    skf = StratifiedKFold(n_splits=n_fold, shuffle=True, random_state=seed)
    cv_aucs = cross_val_score(model, X, y, cv=skf, scoring="roc_auc")
    auc_cv  = float(cv_aucs.mean())
    auc_std = float(cv_aucs.std())

    # --- SHAP: LinearExplainer (exact for linear models) ---
    explainer = shap.LinearExplainer(model, X)
    sv = explainer.shap_values(X)   # (n, n_bits) or list for multiclass
    if isinstance(sv, list):
        sv = sv[1]   # take positive class for binary output

    global_imp = np.mean(np.abs(sv), axis=0)   # (n_bits,)
    coefs      = model.coef_[0]                 # (n_bits,)
    std_x      = X.std(axis=0)                  # (n_bits,)

    # Spearman rho: rank(global_imp) vs rank(|coef| * std(X))
    # For LinearExplainer: global_imp_j ~ |coef_j| * E[|X_j - mean_j|]
    #   which is monotone in |coef_j| * std_j for binary features.
    lr_imp = np.abs(coefs) * std_x
    rho, _  = spearmanr(global_imp, lr_imp)
    shap_lr_spearman = float(rho)

    # --- Select 3 representative molecules ---
    y_pred = model.predict(X)
    y_prob  = model.predict_proba(X)[:, 1]

    # Best TP (high confidence BBB+)
    tp_scores = np.where((y == 1) & (y_pred == 1), y_prob, -np.inf)
    tp_idx    = int(np.argmax(tp_scores))

    # Best TN (high confidence BBB-)
    tn_scores = np.where((y == 0) & (y_pred == 0), 1.0 - y_prob, -np.inf)
    tn_idx    = int(np.argmax(tn_scores))

    # First misclassified molecule
    mis_mask = y != y_pred
    mis_idx  = int(np.where(mis_mask)[0][0]) if mis_mask.any() else tp_idx

    ex_idx   = [tp_idx, tn_idx, mis_idx]
    ex_types = ["TP", "TN", "MIS"]

    # Top-n bits by mean |SHAP| for bar chart axes
    top_n_idx = np.argsort(global_imp)[::-1][:n_top].tolist()

    result = {
        # Metrics
        "auc_cv":            auc_cv,
        "auc_cv_std":        auc_std,
        "n_valid":           n_valid,
        "n_bbb_pos":         int(y.sum()),
        "shap_lr_spearman":  shap_lr_spearman,
        "shap_coverage":     1.0,
        # Global importance (top-n bits)
        "top_n_idx":   top_n_idx,               # 0-based bit indices
        "top_n_imp":   [float(global_imp[i]) for i in top_n_idx],
        "top_n_coef":  [float(coefs[i])      for i in top_n_idx],
        # 3 example molecules
        "ex_idx":    ex_idx,
        "ex_types":  ex_types,
        "ex_labels": [int(y[i])      for i in ex_idx],
        "ex_preds":  [int(y_pred[i]) for i in ex_idx],
        "ex_probs":  [float(y_prob[i]) for i in ex_idx],
        "ex_smiles": [str(df_v["smiles"].iloc[i]) for i in ex_idx],
        "ex_names":  [str(df_v["name"].iloc[i])   for i in ex_idx],
        # Per-example SHAP for top-n bits (n_ex x n_top)
        "ex_top_shap": [[float(sv[i, j]) for j in top_n_idx] for i in ex_idx],
        "ex_top_xval": [[float(X[i,  j]) for j in top_n_idx] for i in ex_idx],
    }

    return json.dumps(result)
