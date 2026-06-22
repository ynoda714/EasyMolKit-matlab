"""rp04_chemberta_core.py  RP04 ChemBERTa Molecular Language Model Core

Extract CLS token embeddings from frozen pre-trained ChemBERTa
(seyonec/ChemBERTa-zinc-base-v1) and train a logistic regression
classifier for BBBP BBB permeability prediction.

Approach: linear probe on frozen pre-trained representations.
  1. Tokenize SMILES -> ChemBERTa (no gradients) -> CLS embedding (768-dim)
  2. StandardScaler + LogisticRegression (L2, C=1.0) on CLS embeddings
  3. 5-fold stratified random CV, return ROC-AUC

Reference: Chithrananda, S. et al. (2020). ChemBERTa: Large-Scale
  Self-Supervised Pretraining for Molecular Property Prediction.
  arXiv:2010.09885. -- same author's model: seyonec/ChemBERTa-zinc-base-v1

Called from rp04_chemberta.m via pyrun:
    exec(open(helper_path).read())
    result_json = run_rp04(csv_path)

Returns JSON string with CV metrics, fold AUCs, and model info.
"""

import json
import numpy as np
import pandas as pd
import torch
from transformers import AutoTokenizer, AutoModel
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from rdkit import Chem

MODEL_NAME = 'seyonec/ChemBERTa-zinc-base-v1'


def _filter_valid(df):
    """Return (smiles_list, y_array) for rows where SMILES parses successfully."""
    smiles, labels = [], []
    for _, row in df.iterrows():
        if Chem.MolFromSmiles(str(row['smiles'])) is not None:
            smiles.append(str(row['smiles']))
            labels.append(int(row['p_np']))
    return smiles, np.array(labels, dtype=np.int32)


def _extract_cls_embeddings(smiles_list, tokenizer, model,
                             batch_size=32, max_len=128):
    """Forward-pass all SMILES through frozen ChemBERTa.

    Two explicit phases:
      Phase 1 — tokenize all SMILES without truncation to collect lengths
                 (pure tokenizer, no model, outside torch.no_grad()).
      Phase 2 — batch forward pass with truncation to extract CLS embeddings
                 (model-only work, inside torch.no_grad()).

    Returns:
        emb: ndarray (n, hidden) -- CLS token embeddings
        raw_lengths: ndarray (n,) int -- token counts without truncation
    """
    # Phase 1: raw token lengths (no model involvement).
    # Items are tokenized individually rather than in batches because
    # batch tokenization with padding=False requires return_tensors=None
    # (list-of-lists output) and adds complexity with no measurable gain
    # at ~2000 molecules. The one-at-a-time approach keeps intent clear.
    raw_lengths = []
    for s in smiles_list:
        ids = tokenizer(s, return_tensors='pt', padding=False, truncation=False)
        raw_lengths.append(int(ids['input_ids'].shape[1]))

    # Phase 2: batch forward pass (frozen model, truncation applied)
    model.eval()
    all_emb = []
    with torch.no_grad():
        for start in range(0, len(smiles_list), batch_size):
            batch = smiles_list[start:start + batch_size]
            enc = tokenizer(batch, return_tensors='pt', padding=True,
                            truncation=True, max_length=max_len)
            out = model(**enc)
            # .detach() required before .numpy() in PyTorch >= 2.0
            cls = out.last_hidden_state[:, 0, :].detach().numpy()
            all_emb.append(cls)
    return np.concatenate(all_emb, axis=0), np.array(raw_lengths)


def _summarise_token_lengths(raw_lengths, max_len=128):
    """Compute distribution statistics for raw (untruncated) token lengths.

    Keys 'n_truncated' / 'frac_truncated' are independent of the max_len value,
    so the JSON schema stays stable even when max_len is varied (e.g. 512).
    """
    arr = raw_lengths
    return {
        'min':           int(arr.min()),
        'max':           int(arr.max()),
        'mean':          round(float(arr.mean()), 2),
        'p50':           int(np.percentile(arr, 50)),
        'p90':           int(np.percentile(arr, 90)),
        'p95':           int(np.percentile(arr, 95)),
        'p99':           int(np.percentile(arr, 99)),
        'n_truncated':   int((arr > max_len).sum()),
        'frac_truncated': round(float((arr > max_len).mean()), 4),
    }


def run_rp04(csv_path, model_name=MODEL_NAME,
             batch_size=None, max_len=128,
             C=1.0, n_fold=5, seed=42):
    """Extract ChemBERTa CLS embeddings + 5-fold CV LR. Return JSON string.

    Args:
        batch_size: SMILES per forward pass. Defaults to max(1, 32 * 128 // max_len)
                    so the token budget stays roughly constant when max_len changes
                    (e.g. max_len=512 -> batch_size=8 to avoid OOM).
                    Pass an explicit integer to override the auto-scaling.
    """

    np.random.seed(seed)
    torch.manual_seed(seed)  # Fix PyTorch RNG for forward-pass consistency

    # Auto-scale batch_size to keep token budget ~constant across max_len values.
    # Baseline: batch_size=32, max_len=128 -> 4096 tokens per batch.
    if batch_size is None:
        batch_size = max(1, 32 * 128 // max_len)

    # Load data & filter valid SMILES
    df = pd.read_csv(csv_path)
    smiles, y = _filter_valid(df)
    n_valid = len(smiles)
    n_pos   = int(y.sum())

    # Load frozen ChemBERTa (model cached in ~/.cache/huggingface/hub/)
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    mdl       = AutoModel.from_pretrained(model_name)
    mdl.eval()

    hidden_size = mdl.config.hidden_size
    n_params    = int(sum(p.numel() for p in mdl.parameters()))

    # Phase 1: token lengths; Phase 2: CLS embeddings (see _extract_cls_embeddings)
    X, raw_lengths = _extract_cls_embeddings(smiles, tokenizer, mdl,
                                              batch_size=batch_size,
                                              max_len=max_len)
    token_stats = _summarise_token_lengths(raw_lengths, max_len=max_len)

    # 5-fold stratified CV: StandardScaler -> LogisticRegression (L2, lbfgs)
    pipe = Pipeline([
        ('scaler', StandardScaler()),
        ('lr', LogisticRegression(C=C, penalty='l2', solver='lbfgs',
                                  max_iter=2000, random_state=seed)),
    ])
    skf     = StratifiedKFold(n_splits=n_fold, shuffle=True, random_state=seed)
    cv_aucs = cross_val_score(pipe, X, y, cv=skf, scoring='roc_auc')

    # Per-fold AUC for comparison with RP02 / RP03
    auc_cv  = float(cv_aucs.mean())
    auc_std = float(cv_aucs.std())

    result = {
        'n_valid':           n_valid,
        'n_bbb_pos':         n_pos,
        'n_bbb_neg':         n_valid - n_pos,
        'auc_cv':            auc_cv,
        'auc_cv_std':        auc_std,
        'fold_aucs':         cv_aucs.tolist(),
        'model_name':        model_name,
        'hidden_size':       hidden_size,
        'n_params_M':        round(n_params / 1e6, 2),
        'embed_dim':         X.shape[1],
        'token_length_stats': token_stats,
        'hyperparams': {
            'batch_size': batch_size,
            'max_len':    max_len,
            'C':          C,
            'seed':       seed,
            'approach':   'frozen_CLS_embedding + LR',
        },
    }
    return json.dumps(result)
