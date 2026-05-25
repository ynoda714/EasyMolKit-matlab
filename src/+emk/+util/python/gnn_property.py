"""
gnn_property.py  --  GCN molecular property prediction for EasyMolKit R09.

Called from MATLAB (dev_sandbox/research/r09_gnn_property.m) via system():

    python gnn_property.py --input <smiles_csv> --output <results_json> [options]

Input CSV must have columns: SMILES, ALogP  (header required).
The script builds PyTorch Geometric Data objects from RDKit, trains a
three-layer GCN with global mean pooling, and writes results to JSON.

Output JSON keys:
    num_train, num_val, num_test       -- dataset sizes after filtering
    train_rmse, val_rmse, test_rmse    -- root mean squared error
    train_mae, val_mae, test_mae       -- mean absolute error
    test_r2                            -- R-squared on test set
    num_atom_features                  -- size of atom feature vector
    hidden_channels, num_layers        -- model hyperparameters
    num_epochs                         -- actual epochs trained
    learning_curves                    -- list of {epoch, train_loss, val_loss}
    predictions                        -- list of {smiles, y_true, y_pred, split}

Requirements:
    torch, torch_geometric (+ torch_scatter, torch_sparse, torch_cluster)
    rdkit
    Install via emk.setup.installExtra("torch") then installExtra("torch_geometric").
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from typing import Optional

# torch must be imported at module level for the @torch.no_grad() decorator
import torch

# ---------------------------------------------------------------------------
# Atom and bond featurisation
# ---------------------------------------------------------------------------

_ELEM_LIST = ["C", "N", "O", "S", "F", "Cl", "Br", "I", "P"]   # 9 common + "other"

_HYBRIDIZATION_MAP = {
    "S":    0,
    "SP":   1,
    "SP2":  2,
    "SP3":  3,
    "SP3D": 4,
    "SP3D2": 5,
}


def _one_hot(value, choices):
    """Return a one-hot list of length len(choices)+1 (last = 'other')."""
    enc = [0] * (len(choices) + 1)
    idx = choices.index(value) if value in choices else len(choices)
    enc[idx] = 1
    return enc


def atom_features(atom):
    """Build a float feature vector for a single RDKit atom.

    Feature layout (33 dimensions total):
      [0:10]  element one-hot (_ELEM_LIST + other)
      [10:17] degree one-hot (0-5 + other)
      [17:23] formal charge one-hot (-2,-1,0,+1,+2 + other)
      [23:30] hybridization one-hot (S,SP,SP2,SP3,SP3D,SP3D2,other)
      [30]    is_aromatic
      [31]    is_in_ring
      [32]    num_explicit_Hs / 4  (normalised to ~[0,1])
    """
    from rdkit.Chem import rdchem

    hyb_name = str(atom.GetHybridization()).split(".")[-1]

    feats = (
        _one_hot(atom.GetSymbol(), _ELEM_LIST)                          # 10
        + _one_hot(atom.GetDegree(), [0, 1, 2, 3, 4, 5])               # 7
        + _one_hot(atom.GetFormalCharge(), [-2, -1, 0, 1, 2])          # 6
        + _one_hot(hyb_name, ["S", "SP", "SP2", "SP3", "SP3D", "SP3D2"])  # 7
        + [int(atom.GetIsAromatic())]                                   # 1
        + [int(atom.IsInRing())]                                        # 1
        + [atom.GetTotalNumHs() / 4.0]                                  # 1
    )
    return feats


NUM_ATOM_FEATURES = 33   # matches the layout above


def smiles_to_data(smiles: str, y_value: float, device=None):
    """Convert a SMILES string to a PyG Data object.

    Returns None if the SMILES is invalid or the molecule has no atoms.
    """
    import torch
    from torch_geometric.data import Data
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    mol = Chem.AddHs(mol)
    mol = Chem.RemoveHs(mol)   # keep implicit Hs only (lighter graph)

    n_atoms = mol.GetNumAtoms()
    if n_atoms == 0:
        return None

    # Node feature matrix [n_atoms, NUM_ATOM_FEATURES]
    x = torch.tensor([atom_features(a) for a in mol.GetAtoms()],
                     dtype=torch.float)

    # Edge index (undirected: each bond adds two directed edges)
    src_list, dst_list = [], []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        src_list += [i, j]
        dst_list += [j, i]

    if len(src_list) == 0:
        # Single-atom molecule -- add self-loop to avoid empty edge_index
        src_list = [0]
        dst_list = [0]

    edge_index = torch.tensor([src_list, dst_list], dtype=torch.long)

    data = Data(x=x, edge_index=edge_index,
                y=torch.tensor([y_value], dtype=torch.float),
                smiles=smiles)
    if device is not None:
        data = data.to(device)
    return data


# ---------------------------------------------------------------------------
# GCN model
# ---------------------------------------------------------------------------

def build_gcn(in_channels: int, hidden_channels: int, num_layers: int):
    """Build a GCN regression model for graph-level property prediction.

    Architecture:
        GCNConv -> BatchNorm -> ReLU   (repeated num_layers times)
        global_mean_pool
        Linear -> scalar output (ALogP)
    """
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.nn import GCNConv, BatchNorm, global_mean_pool

    class GCNModel(nn.Module):
        def __init__(self):
            super().__init__()
            self.convs = nn.ModuleList()
            self.bns   = nn.ModuleList()
            for i in range(num_layers):
                in_ch = in_channels if i == 0 else hidden_channels
                self.convs.append(GCNConv(in_ch, hidden_channels))
                self.bns.append(BatchNorm(hidden_channels))
            self.lin = nn.Linear(hidden_channels, 1)

        def forward(self, x, edge_index, batch):
            for conv, bn in zip(self.convs, self.bns):
                x = conv(x, edge_index)
                x = bn(x)
                x = F.relu(x)
            x = global_mean_pool(x, batch)
            return self.lin(x).squeeze(-1)

    return GCNModel()


# ---------------------------------------------------------------------------
# Training loop helpers
# ---------------------------------------------------------------------------

def train_epoch(model, loader, optimizer, device):
    """Run one training epoch. Returns mean MSE loss over all graphs."""
    import torch
    import torch.nn.functional as F

    model.train()
    total_loss = 0.0
    n_graphs   = 0
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad()
        pred = model(batch.x, batch.edge_index, batch.batch)
        loss = F.mse_loss(pred, batch.y)
        loss.backward()
        optimizer.step()
        total_loss += loss.item() * batch.num_graphs
        n_graphs   += batch.num_graphs
    return total_loss / n_graphs if n_graphs > 0 else float("nan")


@torch.no_grad()
def evaluate(model, loader, device):
    """Evaluate model on a DataLoader. Returns (rmse, mae, y_true, y_pred)."""
    import torch

    model.eval()
    all_pred, all_true = [], []
    for batch in loader:
        batch = batch.to(device)
        pred = model(batch.x, batch.edge_index, batch.batch)
        all_pred.append(pred.cpu())
        all_true.append(batch.y.cpu())

    if not all_pred:
        return float("nan"), float("nan"), [], []

    y_pred = torch.cat(all_pred).numpy()
    y_true = torch.cat(all_true).numpy()
    err    = y_pred - y_true
    rmse   = float(math.sqrt((err ** 2).mean()))
    mae    = float(abs(err).mean())
    return rmse, mae, y_true.tolist(), y_pred.tolist()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main(args):
    import torch
    from torch_geometric.loader import DataLoader

    device = torch.device("cpu")   # R09 targets CPU-only Embedded Python

    # ------------------------------------------------------------------
    # 1. Read input CSV
    # ------------------------------------------------------------------
    smiles_list, y_list = [], []
    with open(args.input, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            s = row.get("SMILES", "").strip()
            try:
                y = float(row.get("ALogP", "nan"))
            except ValueError:
                continue
            if s and not math.isnan(y):
                smiles_list.append(s)
                y_list.append(y)

    print(f"[gnn_property] Read {len(smiles_list)} SMILES from {args.input}",
          flush=True)

    # ------------------------------------------------------------------
    # 2. Build PyG dataset
    # ------------------------------------------------------------------
    dataset, valid_smiles, valid_y = [], [], []
    for s, y in zip(smiles_list, y_list):
        d = smiles_to_data(s, y)
        if d is not None:
            dataset.append(d)
            valid_smiles.append(s)
            valid_y.append(y)

    n_total = len(dataset)
    print(f"[gnn_property] Valid PyG graphs: {n_total} / {len(smiles_list)}",
          flush=True)

    if n_total < 10:
        print("[gnn_property] ERROR: too few valid molecules (< 10). Aborting.",
              flush=True)
        sys.exit(1)

    # ------------------------------------------------------------------
    # 3. Train / val / test split  (70 / 15 / 15)
    # ------------------------------------------------------------------
    import random
    random.seed(42)
    indices = list(range(n_total))
    random.shuffle(indices)

    n_train = int(round(0.70 * n_total))
    n_val   = int(round(0.15 * n_total))
    n_test  = n_total - n_train - n_val

    train_idx = indices[:n_train]
    val_idx   = indices[n_train:n_train + n_val]
    test_idx  = indices[n_train + n_val:]

    split_label = [""] * n_total
    for i in train_idx: split_label[i] = "train"
    for i in val_idx:   split_label[i] = "val"
    for i in test_idx:  split_label[i] = "test"

    train_set = [dataset[i] for i in train_idx]
    val_set   = [dataset[i] for i in val_idx]
    test_set  = [dataset[i] for i in test_idx]

    train_loader = DataLoader(train_set, batch_size=args.batch_size, shuffle=True)
    val_loader   = DataLoader(val_set,   batch_size=args.batch_size, shuffle=False)
    test_loader  = DataLoader(test_set,  batch_size=args.batch_size, shuffle=False)

    print(f"[gnn_property] Split: train={n_train}  val={n_val}  test={n_test}",
          flush=True)

    # ------------------------------------------------------------------
    # 4. Build model and optimiser
    # ------------------------------------------------------------------
    model     = build_gcn(NUM_ATOM_FEATURES, args.hidden_channels, args.num_layers)
    model     = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr,
                                 weight_decay=args.weight_decay)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer,
                                                 step_size=args.lr_step,
                                                 gamma=args.lr_gamma)

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"[gnn_property] GCN: hidden={args.hidden_channels}  "
          f"layers={args.num_layers}  params={n_params}", flush=True)

    # ------------------------------------------------------------------
    # 5. Training loop
    # ------------------------------------------------------------------
    learning_curves = []
    best_val_rmse   = float("inf")
    best_state      = None

    for epoch in range(1, args.epochs + 1):
        train_loss = train_epoch(model, train_loader, optimizer, device)
        val_rmse, val_mae, _, _ = evaluate(model, val_loader, device)
        scheduler.step()

        learning_curves.append({
            "epoch":      epoch,
            "train_loss": round(train_loss, 6),
            "val_rmse":   round(val_rmse,   6),
        })

        if val_rmse < best_val_rmse:
            best_val_rmse = val_rmse
            import copy
            best_state = copy.deepcopy(model.state_dict())

        if epoch % max(1, args.epochs // 10) == 0 or epoch == 1:
            print(f"[gnn_property] Epoch {epoch:4d}/{args.epochs}  "
                  f"train_loss={train_loss:.4f}  val_rmse={val_rmse:.4f}",
                  flush=True)

    # Restore best checkpoint
    if best_state is not None:
        model.load_state_dict(best_state)

    # ------------------------------------------------------------------
    # 6. Final evaluation
    # ------------------------------------------------------------------
    tr_rmse, tr_mae, tr_yt, tr_yp = evaluate(model, train_loader, device)
    va_rmse, va_mae, va_yt, va_yp = evaluate(model, val_loader,   device)
    te_rmse, te_mae, te_yt, te_yp = evaluate(model, test_loader,  device)

    # R-squared on test set
    if te_yt:
        yt_arr = te_yt
        mean_y = sum(yt_arr) / len(yt_arr)
        ss_tot = sum((y - mean_y) ** 2 for y in yt_arr)
        ss_res = sum((yt - yp) ** 2 for yt, yp in zip(te_yt, te_yp))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    else:
        r2 = float("nan")

    print(f"[gnn_property] FINAL  train_rmse={tr_rmse:.4f}  "
          f"val_rmse={va_rmse:.4f}  test_rmse={te_rmse:.4f}  "
          f"test_r2={r2:.4f}", flush=True)

    # ------------------------------------------------------------------
    # 7. Build predictions list
    # ------------------------------------------------------------------
    predictions = []
    # train
    for smiles, yt, yp in zip(
            [valid_smiles[i] for i in train_idx], tr_yt, tr_yp):
        predictions.append({"smiles": smiles, "y_true": round(yt, 4),
                             "y_pred": round(yp, 4), "split": "train"})
    # val
    for smiles, yt, yp in zip(
            [valid_smiles[i] for i in val_idx], va_yt, va_yp):
        predictions.append({"smiles": smiles, "y_true": round(yt, 4),
                             "y_pred": round(yp, 4), "split": "val"})
    # test
    for smiles, yt, yp in zip(
            [valid_smiles[i] for i in test_idx], te_yt, te_yp):
        predictions.append({"smiles": smiles, "y_true": round(yt, 4),
                             "y_pred": round(yp, 4), "split": "test"})

    # ------------------------------------------------------------------
    # 8. Write JSON output
    # ------------------------------------------------------------------
    results = {
        "num_train":          n_train,
        "num_val":            n_val,
        "num_test":           n_test,
        "num_atom_features":  NUM_ATOM_FEATURES,
        "hidden_channels":    args.hidden_channels,
        "num_layers":         args.num_layers,
        "num_epochs":         args.epochs,
        "num_params":         n_params,
        "train_rmse":         round(tr_rmse, 6),
        "val_rmse":           round(va_rmse, 6),
        "test_rmse":          round(te_rmse, 6),
        "train_mae":          round(tr_mae,  6),
        "val_mae":            round(va_mae,  6),
        "test_mae":           round(te_mae,  6),
        "test_r2":            round(r2, 6) if not math.isnan(r2) else None,
        "learning_curves":    learning_curves,
        "predictions":        predictions,
    }

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print(f"[gnn_property] Results written to {args.output}", flush=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Train a GCN to predict molecular properties (EasyMolKit R09)."
    )
    parser.add_argument("--input",   required=True,
                        help="Path to input CSV (columns: SMILES, ALogP)")
    parser.add_argument("--output",  required=True,
                        help="Path to output JSON with results")
    parser.add_argument("--hidden-channels", dest="hidden_channels",
                        type=int, default=64)
    parser.add_argument("--num-layers",      dest="num_layers",
                        type=int, default=3)
    parser.add_argument("--epochs",          type=int, default=150)
    parser.add_argument("--lr",              type=float, default=1e-3)
    parser.add_argument("--weight-decay",    dest="weight_decay",
                        type=float, default=1e-4)
    parser.add_argument("--batch-size",      dest="batch_size",
                        type=int, default=32)
    parser.add_argument("--lr-step",         dest="lr_step",
                        type=int, default=50)
    parser.add_argument("--lr-gamma",        dest="lr_gamma",
                        type=float, default=0.5)
    args = parser.parse_args()

    try:
        import torch
    except ImportError:
        print("ERROR: torch is not installed. Run: emk.setup.installExtra('torch')",
              flush=True)
        sys.exit(2)

    try:
        import torch_geometric
    except ImportError:
        print("ERROR: torch_geometric is not installed. "
              "Run: emk.setup.installExtra('torch_geometric')", flush=True)
        sys.exit(2)

    try:
        from rdkit import Chem
    except ImportError:
        print("ERROR: rdkit is not installed.", flush=True)
        sys.exit(2)

    main(args)
