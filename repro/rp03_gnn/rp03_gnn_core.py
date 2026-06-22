"""rp03_gnn_core.py  RP03: GCN BBBP Classification

Implements leak-free CV for GCN on BBBP:
  - Inner 80/20 val split (train_test_split, fold-specific seed) for early stopping.
    Each outer fold uses random_state=seed+fold_idx so val sets differ across folds.
  - Test fold evaluated once only after best model is restored.
  FIXED:    outer train fold is split 80/20 into sub-train / val;
            early stopping uses val AUC only;
            best model state is saved and restored;
            test fold is evaluated EXACTLY ONCE after early stopping.

Also integrates A2: accepts outer fold indices from rp02_bbbp
(outer_fold_indices.json) so that RP02-rev and RP03-rev use identical splits.
Note: RP03 sub-train (~80% of outer train) differs in size from RP02 outer train;
paired comparison is valid only on the shared outer test folds.

Architecture is unchanged from original:
  3-layer GCNConv + BatchNorm + ReLU + Dropout(0.3)
  -> GlobalMeanPool -> FC(64) -> FC(1)
  BCEWithLogitsLoss with pos_weight, Adam lr=1e-3, StepLR(50, 0.5)
"""

import copy
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from rdkit import Chem
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, train_test_split
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader as PyGLoader
from torch_geometric.nn import GCNConv, global_mean_pool


# ---------------------------------------------------------------------------
# Featurization (unchanged from original)
# ---------------------------------------------------------------------------

_ATOM_TYPES = ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'P', 'Si', 'B']


def _atom_features(atom):
    sym    = atom.GetSymbol()
    at_enc = [int(sym == t) for t in _ATOM_TYPES] + [int(sym not in _ATOM_TYPES)]
    deg    = min(atom.GetDegree(), 5)
    deg_enc = [int(deg == d) for d in range(6)]
    fc     = atom.GetFormalCharge()
    fc_enc = [int(fc == c) for c in [-2, -1, 0, 1, 2]]
    aro    = [int(atom.GetIsAromatic())]
    hcnt   = [min(atom.GetTotalNumHs(), 4) / 4.0]
    return at_enc + deg_enc + fc_enc + aro + hcnt


def _smiles_to_data(smiles, label):
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    x    = torch.tensor([_atom_features(a) for a in mol.GetAtoms()], dtype=torch.float)
    srcs, dsts = [], []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        srcs += [i, j]; dsts += [j, i]
    edge_index = (torch.tensor([srcs, dsts], dtype=torch.long) if srcs
                  else torch.zeros((2, 0), dtype=torch.long))
    return Data(x=x, edge_index=edge_index, y=torch.tensor([label], dtype=torch.float))


# ---------------------------------------------------------------------------
# Model (unchanged from original)
# ---------------------------------------------------------------------------

class GCNClassifier(torch.nn.Module):
    def __init__(self, in_channels=25, hidden=128, n_layers=3, dropout=0.3):
        super().__init__()
        self.convs   = torch.nn.ModuleList()
        self.bns     = torch.nn.ModuleList()
        self.dropout = torch.nn.Dropout(dropout)
        ch = in_channels
        for _ in range(n_layers):
            self.convs.append(GCNConv(ch, hidden))
            self.bns.append(torch.nn.BatchNorm1d(hidden))
            ch = hidden
        self.fc1 = torch.nn.Linear(hidden, hidden // 2)
        self.fc2 = torch.nn.Linear(hidden // 2, 1)

    def forward(self, x, edge_index, batch):
        for conv, bn in zip(self.convs, self.bns):
            x = self.dropout(F.relu(bn(conv(x, edge_index))))
        x = global_mean_pool(x, batch)
        x = self.dropout(F.relu(self.fc1(x)))
        return self.fc2(x).squeeze(-1)


# ---------------------------------------------------------------------------
# Training helpers (unchanged from original)
# ---------------------------------------------------------------------------

def _train_epoch(model, loader, optimizer, loss_fn, device):
    model.train()
    total = 0.0
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        out  = model(batch.x, batch.edge_index, batch.batch)
        loss = loss_fn(out, batch.y.squeeze(-1))
        loss.backward()
        optimizer.step()
        total += loss.item() * batch.num_graphs
    return total / len(loader.dataset)


@torch.no_grad()
def _eval_auc(model, loader, device):
    model.eval()
    probs, labels = [], []
    for batch in loader:
        batch = batch.to(device)
        out   = torch.sigmoid(model(batch.x, batch.edge_index, batch.batch))
        probs.extend(out.cpu().tolist())
        labels.extend(batch.y.squeeze(-1).cpu().tolist())
    if len(set(labels)) < 2:
        return 0.5
    return float(roc_auc_score(labels, probs))


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def run_rp03(csv_path, fold_indices_path=None,
                  hidden=128, n_layers=3, dropout=0.3,
                  lr=1e-3, n_epochs=150, patience=20,
                  n_fold=5, batch_size=64, seed=42):
    """Train GCN with proper train/val/test separation. Return JSON string.

    Args:
        csv_path:          Path to data/benchmark/bbbp.csv.
        fold_indices_path: Path to outer_fold_indices.json from rp02_bbbp (A2).
                           If None, generates fresh StratifiedKFold(seed).

    Inner val split: train_test_split(test_size=0.2, stratify=y, random_state=seed+fold_idx).
    Each outer fold uses a distinct seed so val sets differ across folds.

    Reproducibility note: DataLoader uses num_workers=0 (default). The fold-specific
    Generator guarantees deterministic shuffle within a fold at num_workers=0.
    If num_workers > 0 is added in future, a worker_init_fn seeded per worker is
    also required for bit-exact reproducibility.
    """
    torch.manual_seed(seed)
    np.random.seed(seed)
    device = torch.device('cpu')

    # Load and featurize
    df = pd.read_csv(csv_path)
    data_list, y_list = [], []
    for _, row in df.iterrows():
        d = _smiles_to_data(row['smiles'], int(row['p_np']))
        if d is not None:
            data_list.append(d)
            y_list.append(int(row['p_np']))

    y_arr    = np.array(y_list, dtype=np.int32)
    n_valid  = len(data_list)
    n_pos    = int(y_arr.sum())
    atom_dim = data_list[0].x.shape[1]

    # Outer fold splits (A2: use shared indices if available)
    if fold_indices_path and Path(str(fold_indices_path)).exists():
        with open(str(fold_indices_path)) as f:
            fi = json.load(f)
        outer_splits = [
            (np.array(fi[f"fold{k}_train"], dtype=int),
             np.array(fi[f"fold{k}_test"],  dtype=int))
            for k in range(1, n_fold + 1)
        ]
        fold_source = "a2_shared_rp02rev"
    else:
        skf = StratifiedKFold(n_splits=n_fold, shuffle=True, random_state=seed)
        outer_splits = list(skf.split(np.zeros(n_valid), y_arr))
        fold_source  = "internal_seed42"

    fold_aucs  = []
    fold_curves = []

    for fold_idx, (tr_val_idx, te_idx) in enumerate(outer_splits):
        # Fold-specific seeds ensure bit-exact per-fold reproducibility.
        torch.manual_seed(seed + fold_idx)
        np.random.seed(seed + fold_idx)
        # Inner 80/20 split: fold-specific seed avoids identical val sets across folds.
        sub_tr_local, sub_val_local = train_test_split(
            np.arange(len(tr_val_idx)), test_size=0.2,
            stratify=y_arr[tr_val_idx], random_state=seed + fold_idx)
        train_idx = tr_val_idx[sub_tr_local]
        val_idx   = tr_val_idx[sub_val_local]

        # generator pins DataLoader shuffle order independently of global RNG state.
        fold_gen  = torch.Generator().manual_seed(seed + fold_idx)
        tr_loader  = PyGLoader([data_list[i] for i in train_idx],
                                batch_size=batch_size, shuffle=True,
                                generator=fold_gen)
        val_loader = PyGLoader([data_list[i] for i in val_idx],
                                batch_size=batch_size, shuffle=False)
        te_loader  = PyGLoader([data_list[i] for i in te_idx],
                                batch_size=batch_size, shuffle=False)

        # Class weight from sub-train only
        n_tr_pos = int(y_arr[train_idx].sum())
        n_tr_neg = len(train_idx) - n_tr_pos
        # BBB+ (label=1) is the majority class (~76.5%).  pos_weight = n_neg/n_pos < 1
        # down-weights the majority class loss, giving the minority class (BBB-)
        # relatively higher gradient contribution to counteract the imbalance.
        pw       = torch.tensor([n_tr_neg / max(n_tr_pos, 1)], dtype=torch.float)
        loss_fn  = torch.nn.BCEWithLogitsLoss(pos_weight=pw)

        model     = GCNClassifier(in_channels=atom_dim, hidden=hidden,
                                   n_layers=n_layers, dropout=dropout).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=50, gamma=0.5)

        best_val_auc = 0.0
        best_state   = None
        no_improve   = 0
        tr_losses    = []
        val_aucs     = []

        for epoch in range(n_epochs):
            loss    = _train_epoch(model, tr_loader, optimizer, loss_fn, device)
            val_auc = _eval_auc(model, val_loader, device)   # VAL only; test unseen
            scheduler.step()

            tr_losses.append(float(loss))
            val_aucs.append(float(val_auc))

            if val_auc > best_val_auc:
                best_val_auc = val_auc
                best_state   = copy.deepcopy(model.state_dict())
                no_improve   = 0
            else:
                no_improve += 1
            if no_improve >= patience:
                break

        # Restore best model; evaluate test fold ONCE
        if best_state is not None:
            model.load_state_dict(best_state)
        test_auc = _eval_auc(model, te_loader, device)

        fold_aucs.append(float(test_auc))
        fold_curves.append({
            'fold':         fold_idx + 1,
            'train_size':   int(len(train_idx)),
            'val_size':     int(len(val_idx)),
            'test_size':    int(len(te_idx)),
            'best_val_auc': float(best_val_auc),
            'test_auc':     float(test_auc),
            'n_epochs':     len(tr_losses),
            'train_loss':   tr_losses,
            'val_auc_curve': val_aucs,
        })

    fold_epoch_counts = [fc['n_epochs'] for fc in fold_curves]
    min_len  = min(len(fc['train_loss']) for fc in fold_curves)
    avg_loss = np.mean([fc['train_loss'][:min_len] for fc in fold_curves], axis=0).tolist()
    avg_vauc = np.mean([fc['val_auc_curve'][:min_len] for fc in fold_curves], axis=0).tolist()

    lc_note = (
        f"avg_train_loss / avg_val_auc truncated to {min_len} epochs "
        f"(shortest early-stopped fold); fold epoch counts: {fold_epoch_counts}."
    )

    return json.dumps({
        'n_valid':              n_valid,
        'n_bbb_pos':            n_pos,
        'n_bbb_neg':            n_valid - n_pos,
        'fold_source':          fold_source,
        'auc_cv':               float(np.mean(fold_aucs)),
        'auc_cv_std':           float(np.std(fold_aucs)),
        'fold_aucs':            fold_aucs,
        'fold_curves':          fold_curves,
        'atom_feat_dim':        int(atom_dim),
        'avg_train_loss':       avg_loss,
        'avg_val_auc':          avg_vauc,
        'n_epochs_run':         min_len,
        'learning_curve_note':  lc_note,
        'hyperparams': {
            'hidden': hidden, 'n_layers': n_layers, 'dropout': dropout,
            'lr': lr, 'n_epochs': n_epochs, 'patience': patience,
            'batch_size': batch_size, 'seed': seed,
            'inner_val_split': 'train_test_split(test_size=0.2, random_state=seed+fold_idx)',
        },
    })
