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
import os
import random
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
MAX_LEN = 128  # default token limit; raise to 512 to capture longer SMILES (see README B4)


def _hf_local_only():
    """Return True when the current process is configured for offline HF use.

    This keeps the repro compatible with ordinary online runs while avoiding
    background metadata requests in restricted environments where the model is
    already cached locally.
    """
    truthy = {'1', 'true', 'yes', 'on'}
    return (
        str(os.environ.get('HF_HUB_OFFLINE', '')).strip().lower() in truthy or
        str(os.environ.get('TRANSFORMERS_OFFLINE', '')).strip().lower() in truthy or
        str(os.environ.get('EMK_HF_LOCAL_ONLY', '')).strip().lower() in truthy
    )


def _filter_valid(df):
    """Return (smiles_list, y_array) for rows where SMILES parses successfully."""
    smiles, labels = [], []
    for _, row in df.iterrows():
        if Chem.MolFromSmiles(str(row['smiles'])) is not None:
            smiles.append(str(row['smiles']))
            labels.append(int(row['p_np']))
    return smiles, np.array(labels, dtype=np.int32)


def _extract_cls_embeddings(smiles_list, tokenizer, model,
                             batch_size=32, max_len=MAX_LEN):
    """Forward-pass all SMILES through frozen ChemBERTa.

    Two explicit phases:
      Phase 1 — tokenize all SMILES without truncation to collect lengths
                 (pure tokenizer, no model, outside torch.no_grad()).
      Phase 2 — batch forward pass with truncation to extract CLS embeddings
                 (model-only work, inside torch.no_grad()).

    Caller must have set model.eval() before calling.

    Returns:
        emb: ndarray (n, hidden) -- CLS token embeddings
        raw_lengths: ndarray (n,) int -- token counts without truncation
    """
    # Phase 1: batch tokenization to collect raw lengths (no model involvement).
    # Use tokenizer() callable (works across transformers versions; batch_encode_plus deprecated).
    encoded = tokenizer(smiles_list, padding=False, truncation=False)
    raw_lengths = [len(ids) for ids in encoded['input_ids']]

    # Phase 2: batch forward pass (frozen model, truncation applied).
    # model.eval() is set by the caller (run_rp04); no redundant call here.
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


def _summarise_token_lengths(raw_lengths, max_len=MAX_LEN):
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
             batch_size=None, max_len=MAX_LEN,
             C=1.0, n_fold=5, seed=42):
    """Extract ChemBERTa CLS embeddings + 5-fold CV LR. Return JSON string.

    Args:
        batch_size: SMILES per forward pass. Defaults to max(1, 32 * 128 // max_len)
                    as a heuristic to keep the token budget ~constant across max_len
                    values (e.g. max_len=512 -> batch_size=8). Actual memory scales
                    with sequence_len * batch * hidden_size, which the formula ignores;
                    tune explicitly on memory-constrained hardware.
    """
    # Fix all Python/NumPy/PyTorch RNGs before any library call.
    # transformers uses random, numpy, and torch internally.
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    # Auto-scale batch_size to keep token budget ~constant across max_len values.
    # Baseline: batch_size=32, max_len=128 -> 4096 tokens per batch.
    # KNOWN LIMITATION: hidden_size (768) is not factored in; OOM is possible at
    # max_len=512 on low-memory hardware. Pass batch_size explicitly if needed.
    if batch_size is None:
        batch_size = max(1, 32 * 128 // max_len)

    # Load data & filter valid SMILES
    df = pd.read_csv(csv_path)
    smiles, y = _filter_valid(df)
    n_valid = len(smiles)
    n_pos   = int(y.sum())

    # Load frozen ChemBERTa (model cached in ~/.cache/huggingface/hub/)
    local_only = _hf_local_only()
    tokenizer = AutoTokenizer.from_pretrained(model_name, local_files_only=local_only)
    mdl       = AutoModel.from_pretrained(model_name, local_files_only=local_only)
    mdl.eval()  # set once here; _extract_cls_embeddings expects eval mode

    hidden_size = mdl.config.hidden_size
    # Counts all parameters regardless of requires_grad; equals total frozen count
    # since no optimizer is attached and no .backward() is called.
    n_params    = int(sum(p.numel() for p in mdl.parameters()))

    # Phase 1: token lengths; Phase 2: CLS embeddings (see _extract_cls_embeddings)
    X, raw_lengths = _extract_cls_embeddings(smiles, tokenizer, mdl,
                                              batch_size=batch_size,
                                              max_len=max_len)
    token_stats = _summarise_token_lengths(raw_lengths, max_len=max_len)

    # 5-fold stratified CV: StandardScaler -> LogisticRegression (L2, lbfgs)
    pipe = Pipeline([
        ('scaler', StandardScaler()),
        # lbfgs is deterministic; random_state has no effect but is kept for
        # forward-compatibility if the solver is ever changed to a stochastic one.
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


class _ChemBERTaExportWrapper(torch.nn.Module):
    """Thin wrapper to present a clean (input_ids, attention_mask) -> last_hidden_state
    signature for torch.onnx.export tracing.

    RobertaModel.forward() has many optional kwargs that confuse TorchScript tracing
    (e.g. use_cache gets passed both positionally and as a keyword).  The wrapper
    fixes the interface to exactly two inputs and one output.
    """

    def __init__(self, encoder):
        super().__init__()
        self.encoder = encoder

    def forward(self, input_ids, attention_mask):
        out = self.encoder(input_ids=input_ids, attention_mask=attention_mask)
        return out.last_hidden_state


def export_chemberta_onnx(save_path, model_name=MODEL_NAME, max_len=MAX_LEN):
    """Export frozen ChemBERTa to ONNX format for MATLAB importNetworkFromONNX.

    Uses fixed seq_len=max_len (dynamic batch_size) to avoid MATLAB ONNX
    importer issues with fully-dynamic shapes.  Opset 14 is chosen for broad
    MATLAB R2025a/R2026a compatibility.

    Returns JSON string with export metadata.
    """
    import os

    local_only = _hf_local_only()
    tokenizer = AutoTokenizer.from_pretrained(model_name, local_files_only=local_only)
    mdl       = AutoModel.from_pretrained(model_name, local_files_only=local_only)
    mdl.eval()

    hidden_size = mdl.config.hidden_size
    n_params    = int(sum(p.numel() for p in mdl.parameters()))

    wrapper = _ChemBERTaExportWrapper(mdl)
    wrapper.eval()

    # Dummy input: batch=1, fixed seq_len=max_len for tracing
    dummy = ['CCO']
    enc   = tokenizer(dummy, return_tensors='pt',
                      padding='max_length', truncation=True, max_length=max_len)
    input_ids      = enc['input_ids']       # (1, max_len) int64
    attention_mask = enc['attention_mask']  # (1, max_len) int64

    # dynamo=True (default in PyTorch >= 2.9) uses torch.export + onnxscript backend.
    # opset_version is advisory for the dynamo exporter; actual opset is >= 18.
    with torch.no_grad():
        torch.onnx.export(
            wrapper,
            (input_ids, attention_mask),
            save_path,
            input_names=['input_ids', 'attention_mask'],
            output_names=['last_hidden_state'],
            dynamic_axes={
                'input_ids':         {0: 'batch_size'},
                'attention_mask':    {0: 'batch_size'},
                'last_hidden_state': {0: 'batch_size'},
            },
            opset_version=14,
            verbose=False,
        )

    file_size_mb = round(os.path.getsize(save_path) / 1e6, 1)
    result = {
        'save_path':    save_path,
        'model_name':   model_name,
        'max_len':      max_len,
        'hidden_size':  hidden_size,
        'n_params_M':   round(n_params / 1e6, 2),
        'opset':        14,
        'file_size_mb': file_size_mb,
        'input_names':  ['input_ids', 'attention_mask'],
        'output_names': ['last_hidden_state'],
        'input_shape':  [1, max_len],
        'output_shape': [1, max_len, hidden_size],
    }
    return json.dumps(result)


def patch_onnx_for_matlab(input_path, output_path):
    """Remove IsNaN -> Where patterns from ONNX graph for MATLAB compatibility.

    MATLAB R2026a's ONNX importer does not support the 'IsNaN' operator.
    In ChemBERTa the pattern arises from nan_to_num() after Softmax in each
    attention layer (guards against all-masked padding rows).

    Safety: CLS token (position 0) is always a real token, so softmax over
    at least the CLS column never produces NaN for valid SMILES input.
    Bypassing IsNaN -> Where is therefore a no-op for our use case.

    Procedure:
      For each IsNaN node whose output feeds a Where node:
        1. Record that Where's output should be replaced by IsNaN's input
           (i.e., the raw softmax output).
        2. Remove IsNaN and Where nodes from the graph.
        3. Rewire downstream consumers to use the softmax output directly.

    Returns JSON string with patch metadata.

    Implementation note: protobuf repeated-field elements return new Python
    wrapper objects on each access, so id(node) is not stable across two
    iterations of graph.node.  Index-based tracking avoids this pitfall.
    Similarly, copy.deepcopy() is required before del graph.node[:] because
    clearing the repeated field may invalidate existing Python references to
    the underlying C++ protobuf objects.
    """
    import onnx
    import copy

    model = onnx.load(input_path)
    graph = model.graph

    # Snapshot the node list once to get stable indices and references.
    node_list = list(graph.node)

    # Map: isnan_output_name -> softmax_output_name (IsNaN's input)
    isnan_to_src = {}
    for node in node_list:
        if node.op_type == 'IsNaN':
            isnan_to_src[node.output[0]] = node.input[0]

    # Map: where_output_name -> softmax_output_name (bypass NaN handling)
    # Track removal by index (id() is unreliable for protobuf repeated elements).
    output_remap = {}
    skip_indices = set()
    for i, node in enumerate(node_list):
        if node.op_type == 'Where' and node.input[0] in isnan_to_src:
            softmax_out = isnan_to_src[node.input[0]]
            output_remap[node.output[0]] = softmax_out
            skip_indices.add(i)
    for i, node in enumerate(node_list):
        if node.op_type == 'IsNaN' and node.output[0] in isnan_to_src:
            skip_indices.add(i)

    # Rebuild node list: deep-copy every kept node (protobuf refs are
    # invalidated by del graph.node[:]) then remap inputs.
    new_nodes = []
    for i, node in enumerate(node_list):
        if i in skip_indices:
            continue
        node_copy = copy.deepcopy(node)
        for j, inp in enumerate(node_copy.input):
            if inp in output_remap:
                node_copy.input[j] = output_remap[inp]
        new_nodes.append(node_copy)

    del graph.node[:]
    graph.node.extend(new_nodes)
    onnx.save(model, output_path)

    result = {
        'input_path':          input_path,
        'output_path':         output_path,
        'isnan_pairs_removed': len(output_remap),
    }
    return json.dumps(result)


def prepare_onnx_for_matlab(input_path, output_path):
    """Patch ChemBERTa ONNX for MATLAB importNetworkFromONNX compatibility.

    Two-step pipeline applied to the single-file ONNX from export_chemberta_onnx:
      1. patch_onnx_for_matlab  -- removes 6 IsNaN->Where nan_to_num patterns
      2. IR version downgrade   -- sets ir_version=9 (MATLAB max fully-supported)

    Wrapping both steps here avoids an intermediate file and ensures the
    combination is always applied together.

    Returns JSON string with patch metadata.
    """
    import json, os, tempfile, onnx

    with tempfile.NamedTemporaryFile(suffix='.onnx', delete=False) as tf:
        tmp_path = tf.name

    try:
        patch_result = json.loads(patch_onnx_for_matlab(input_path, tmp_path))
        m = onnx.load(tmp_path)
        ir_before = m.ir_version
        m.ir_version = 9
        onnx.checker.check_model(m)
        onnx.save(m, output_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    result = {
        'input_path':          input_path,
        'output_path':         output_path,
        'isnan_pairs_removed': patch_result['isnan_pairs_removed'],
        'ir_version_before':   ir_before,
        'ir_version_after':    9,
        'file_size_mb':        round(os.path.getsize(output_path) / 1e6, 1),
    }
    return json.dumps(result)


def tokenize_for_matlab(csv_path, model_name=MODEL_NAME, max_len=MAX_LEN):
    """Tokenise BBBP SMILES only (no model inference). Return JSON.

    For F1-a: Python handles tokenisation only; MATLAB handles ChemBERTa
    inference (via imported ONNX network) and LR classification.

    Returns JSON with:
        input_ids      -- n_valid x max_len int matrix (padded to max_len)
        attention_mask -- n_valid x max_len int matrix
        labels         -- n_valid binary int list (0=BBB-, 1=BBB+)
    """
    df = pd.read_csv(csv_path)
    smiles, y = _filter_valid(df)
    n_valid = len(smiles)
    n_pos   = int(y.sum())

    local_only = _hf_local_only()
    tokenizer = AutoTokenizer.from_pretrained(model_name, local_files_only=local_only)
    enc = tokenizer(smiles, padding='max_length', truncation=True,
                    max_length=max_len, return_tensors='pt')

    result = {
        'n_valid':        n_valid,
        'n_bbb_pos':      n_pos,
        'n_bbb_neg':      n_valid - n_pos,
        'input_ids':      enc['input_ids'].numpy().astype(int).tolist(),
        'attention_mask': enc['attention_mask'].numpy().astype(int).tolist(),
        'labels':         y.tolist(),
        'max_len':        max_len,
        'model_name':     model_name,
    }
    return json.dumps(result)


def extract_embeddings_only(csv_path, model_name=MODEL_NAME,
                             batch_size=None, max_len=MAX_LEN, seed=42):
    """Extract ChemBERTa CLS embeddings only (no LR). Return JSON string.

    Zone C (F1): Python handles tokenisation + CLS inference.
    MATLAB uses the returned embeddings for fitclinear logistic regression.

    Returns a JSON string containing:
        embeddings  -- n_valid x hidden_size float matrix (as nested list)
        labels      -- n_valid binary int list (0=BBB-, 1=BBB+)
        metadata    -- model info and token statistics
    """
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)

    if batch_size is None:
        batch_size = max(1, 32 * 128 // max_len)

    df = pd.read_csv(csv_path)
    smiles, y = _filter_valid(df)
    n_valid = len(smiles)
    n_pos   = int(y.sum())

    local_only = _hf_local_only()
    tokenizer = AutoTokenizer.from_pretrained(model_name, local_files_only=local_only)
    mdl       = AutoModel.from_pretrained(model_name, local_files_only=local_only)
    mdl.eval()

    hidden_size = mdl.config.hidden_size
    n_params    = int(sum(p.numel() for p in mdl.parameters()))

    X, raw_lengths = _extract_cls_embeddings(smiles, tokenizer, mdl,
                                              batch_size=batch_size,
                                              max_len=max_len)
    token_stats = _summarise_token_lengths(raw_lengths, max_len=max_len)

    result = {
        'n_valid':            n_valid,
        'n_bbb_pos':          n_pos,
        'n_bbb_neg':          n_valid - n_pos,
        'embeddings':         X.tolist(),   # n_valid x hidden_size (JSON matrix)
        'labels':             y.tolist(),   # n_valid binary ints
        'model_name':         model_name,
        'hidden_size':        hidden_size,
        'n_params_M':         round(n_params / 1e6, 2),
        'embed_dim':          int(X.shape[1]),
        'token_length_stats': token_stats,
        'hyperparams': {
            'batch_size': batch_size,
            'max_len':    max_len,
            'seed':       seed,
            'approach':   'frozen_CLS_embedding_only',
        },
    }
    return json.dumps(result)
