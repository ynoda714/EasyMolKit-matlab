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
from scipy.stats import spearmanr


def run_rp05(csv_path, radius=2, n_bits=2048, C=1.0,
             n_fold=5, seed=42, n_top=20, test_size=0.2):
    """Fit LR on BBBP ECFP4, compute SHAP values, return JSON string."""

    # --- Load and validate data ---
    df = pd.read_csv(csv_path)
    # Explicit NaN check first; after notna() filter str() cast is safe but the
    # MolFromSmiles call itself is the definitive validity gate.
    valid_mask = df["smiles"].notna() & np.array(
        [Chem.MolFromSmiles(s) is not None for s in df["smiles"]]
    )
    df_v = df[valid_mask].reset_index(drop=True)
    n_valid = len(df_v)

    # --- Compute Morgan ECFP4 fingerprints ---
    gen = rfg.GetMorganGenerator(radius=radius, fpSize=n_bits)
    mols = [Chem.MolFromSmiles(str(s)) for s in df_v["smiles"]]
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
    explainer = shap.LinearExplainer(model, X_tr)
    sv = explainer.shap_values(X_te)
    if isinstance(sv, list):
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
    rho, _  = spearmanr(global_imp, lr_imp)
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
        mis_conf = np.where(y_pred == 1, y_prob, 1.0 - y_prob)
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
        "shap_lr_spearman":  shap_lr_spearman,
        # Global importance (top-n bits)
        "top_n_idx":   top_n_idx,               # 0-based bit indices
        "top_n_imp":   [float(global_imp[i]) for i in top_n_idx],
        "top_n_coef":  [float(coefs[i])      for i in top_n_idx],
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
