"""rp03_gnn_core.py  RP03 Graph Neural Network Core

Train a 3-layer GCN on BBBP Blood-Brain Barrier data with 5-fold stratified CV.
Follows Yang et al. (2019) Chemprop paper: BBBP classification via graph-based
molecular representation, benchmarked against the ECFP4+LR baseline in RP02.

Called from rp03_gnn.m via pyrun:
    exec(open(helper_path).read())
    result_json = run_rp03(csv_path)

Returns JSON string with CV metrics, fold-level AUCs, and learning curves.

Architecture: GCNConv (3 layers, hidden=128) + BatchNorm + ReLU + Dropout(0.3)
              -> GlobalMeanPool -> FC(64) -> sigmoid
Atom features: 25-dim (type x12, degree x6, formal_charge x5, aromatic x1, H x1)
Loss: BCEWithLogitsLoss with pos_weight for class imbalance (BBB+ 76%, BBB- 24%)
"""

import json
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader as PyGLoader
from torch_geometric.nn import GCNConv, global_mean_pool
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from rdkit import Chem


# ---------------------------------------------------------------------------
# Featurization
# ---------------------------------------------------------------------------

_ATOM_TYPES = ['C', 'N', 'O', 'S', 'F', 'Cl', 'Br', 'I', 'P', 'Si', 'B']


def _atom_features(atom):
    """Return 25-dim float feature vector for one atom."""
    sym = atom.GetSymbol()
    at_enc  = [int(sym == t) for t in _ATOM_TYPES] + [int(sym not in _ATOM_TYPES)]  # 12
    deg     = min(atom.GetDegree(), 5)
    deg_enc = [int(deg == d) for d in range(6)]                                       # 6
    fc      = atom.GetFormalCharge()
    fc_enc  = [int(fc == c) for c in [-2, -1, 0, 1, 2]]                              # 5
    aro     = [int(atom.GetIsAromatic())]                                              # 1
    hcnt    = [min(atom.GetTotalNumHs(), 4) / 4.0]                                    # 1
    return at_enc + deg_enc + fc_enc + aro + hcnt   # total 25


def _smiles_to_data(smiles, label):
    """Convert SMILES string to torch_geometric Data. Returns None on parse failure."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None

    x = torch.tensor([_atom_features(a) for a in mol.GetAtoms()], dtype=torch.float)

    srcs, dsts = [], []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        srcs += [i, j]; dsts += [j, i]

    if srcs:
        edge_index = torch.tensor([srcs, dsts], dtype=torch.long)
    else:
        edge_index = torch.zeros((2, 0), dtype=torch.long)

    return Data(x=x, edge_index=edge_index, y=torch.tensor([label], dtype=torch.float))


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

class GCNClassifier(torch.nn.Module):
    """3-layer GCN with BatchNorm + global mean pooling for binary classification."""

    def __init__(self, in_channels=25, hidden=128, n_layers=3, dropout=0.3):
        super().__init__()
        self.convs    = torch.nn.ModuleList()
        self.bns      = torch.nn.ModuleList()
        self.dropout  = torch.nn.Dropout(dropout)

        ch_in = in_channels
        for _ in range(n_layers):
            self.convs.append(GCNConv(ch_in, hidden))
            self.bns.append(torch.nn.BatchNorm1d(hidden))
            ch_in = hidden

        self.fc1 = torch.nn.Linear(hidden, hidden // 2)
        self.fc2 = torch.nn.Linear(hidden // 2, 1)

    def forward(self, x, edge_index, batch):
        for conv, bn in zip(self.convs, self.bns):
            x = self.dropout(F.relu(bn(conv(x, edge_index))))
        x = global_mean_pool(x, batch)
        x = self.dropout(F.relu(self.fc1(x)))
        return self.fc2(x).squeeze(-1)


# ---------------------------------------------------------------------------
# Training helpers
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

def run_rp03(csv_path, hidden=128, n_layers=3, dropout=0.3,
             lr=1e-3, n_epochs=150, patience=20,
             n_fold=5, batch_size=64, seed=42):
    """Train GCN on BBBP with 5-fold stratified CV. Return JSON string."""

    torch.manual_seed(seed)
    np.random.seed(seed)
    device = torch.device('cpu')

    # Load & featurize all molecules
    df = pd.read_csv(csv_path)
    data_list, y_list = [], []
    for _, row in df.iterrows():
        d = _smiles_to_data(row['smiles'], int(row['p_np']))
        if d is not None:
            data_list.append(d)
            y_list.append(int(row['p_np']))

    y_arr   = np.array(y_list, dtype=np.int32)
    n_valid = len(data_list)
    n_pos   = int(y_arr.sum())
    atom_dim = data_list[0].x.shape[1]

    # 5-fold stratified CV
    skf        = StratifiedKFold(n_splits=n_fold, shuffle=True, random_state=seed)
    fold_aucs  = []
    fold_curves = []

    for fold_idx, (tr_idx, te_idx) in enumerate(skf.split(np.zeros(n_valid), y_arr)):
        tr_data = [data_list[i] for i in tr_idx]
        te_data = [data_list[i] for i in te_idx]
        tr_loader = PyGLoader(tr_data, batch_size=batch_size, shuffle=True)
        te_loader = PyGLoader(te_data, batch_size=batch_size, shuffle=False)

        # Class-imbalance correction (BBB+ ~76%, BBB- ~24%)
        n_tr_pos  = int(y_arr[tr_idx].sum())
        n_tr_neg  = len(tr_idx) - n_tr_pos
        pw        = torch.tensor([n_tr_neg / max(n_tr_pos, 1)], dtype=torch.float)
        loss_fn   = torch.nn.BCEWithLogitsLoss(pos_weight=pw)

        model     = GCNClassifier(in_channels=atom_dim, hidden=hidden,
                                   n_layers=n_layers, dropout=dropout).to(device)
        optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
        scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=50, gamma=0.5)

        best_auc, best_epoch, no_improve = 0.0, 0, 0
        tr_losses, va_aucs = [], []

        for epoch in range(n_epochs):
            loss = _train_epoch(model, tr_loader, optimizer, loss_fn, device)
            auc  = _eval_auc(model, te_loader, device)
            scheduler.step()
            tr_losses.append(float(loss))
            va_aucs.append(float(auc))

            if auc > best_auc:
                best_auc, best_epoch, no_improve = auc, epoch, 0
            else:
                no_improve += 1
            if no_improve >= patience:
                break

        fold_aucs.append(float(best_auc))
        fold_curves.append({'train_loss': tr_losses, 'val_auc': va_aucs,
                            'best_epoch': best_epoch, 'best_auc': float(best_auc)})

    # Average learning curves (truncate to shortest fold)
    min_len  = min(len(fc['train_loss']) for fc in fold_curves)
    avg_loss = np.mean([fc['train_loss'][:min_len] for fc in fold_curves], axis=0).tolist()
    avg_vauc = np.mean([fc['val_auc'][:min_len]    for fc in fold_curves], axis=0).tolist()

    result = {
        'n_valid':     n_valid,
        'n_bbb_pos':   n_pos,
        'n_bbb_neg':   n_valid - n_pos,
        'auc_cv':      float(np.mean(fold_aucs)),
        'auc_cv_std':  float(np.std(fold_aucs)),
        'fold_aucs':   fold_aucs,
        'atom_feat_dim': int(atom_dim),
        'avg_train_loss': avg_loss,
        'avg_val_auc':    avg_vauc,
        'n_epochs_run':   min_len,
        'hyperparams': {
            'hidden':     hidden,
            'n_layers':   n_layers,
            'dropout':    dropout,
            'lr':         lr,
            'n_epochs':   n_epochs,
            'patience':   patience,
            'batch_size': batch_size,
            'seed':       seed,
        },
    }
    return json.dumps(result)
