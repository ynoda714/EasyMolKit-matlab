"""
chemberta.py  --  ChemBERTa embedding and fine-tuning for EasyMolKit R10.

Called from MATLAB (dev_sandbox/research/r10_chemberta.m) via system():

Modes
-----
1. Embed mode (--mode embed):
       python chemberta.py --mode embed --input <smiles_csv> --output <results_json>
                           [--model-name <hf_model_name>] [--batch-size <N>]

   Input CSV: must have a SMILES column (header required).
   Output JSON keys:
       model_name          -- HuggingFace model name used
       num_smiles          -- number of input SMILES
       num_valid           -- SMILES successfully tokenised
       embedding_dim       -- hidden dimension of the [CLS] embedding
       embeddings          -- list of {smiles, embedding (float list), valid (bool)}

2. Fine-tune mode (--mode finetune):
       python chemberta.py --mode finetune --input <smiles_csv> --output <results_json>
                           [--target-col <col>] [--model-name <hf_model_name>]
                           [--epochs <N>] [--lr <float>] [--batch-size <N>]
                           [--freeze-encoder] [--seed <int>]

   Input CSV: must have SMILES column and a numeric target column (default: ALogP).
   Output JSON keys:
       model_name          -- HuggingFace model name used
       target_col          -- regression target column
       num_train, num_val, num_test  -- dataset sizes
       embedding_dim       -- [CLS] embedding dimension
       freeze_encoder      -- whether encoder weights were frozen during training
       num_epochs          -- epochs trained
       num_params          -- total trainable parameters
       train_rmse, val_rmse, test_rmse
       train_mae,  val_mae,  test_mae
       test_r2
       learning_curves     -- list of {epoch, train_loss, val_rmse}
       predictions         -- list of {smiles, y_true, y_pred, split}

Requirements:
    torch, transformers
    rdkit  (for SMILES validation only)
    Install via:
        emk.setup.installExtra("torch")
        emk.setup.installExtra("transformers")
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from typing import Optional

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
from transformers import AutoTokenizer, AutoModel

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_MODEL = "seyonec/ChemBERTa-zinc-base-v1"
MAX_LENGTH = 128          # max tokeniser length for SMILES
SMILES_COLUMN = "SMILES"  # expected column name in input CSV


# ---------------------------------------------------------------------------
# SMILES tokenisation helpers
# ---------------------------------------------------------------------------

def load_tokenizer(model_name: str):
    """Load the HuggingFace tokeniser for model_name."""
    return AutoTokenizer.from_pretrained(model_name)


def load_encoder(model_name: str):
    """Load the pre-trained encoder (transformer backbone)."""
    model = AutoModel.from_pretrained(model_name)
    return model


# ---------------------------------------------------------------------------
# Dataset classes
# ---------------------------------------------------------------------------

class SMILESDataset(Dataset):
    """Dataset that tokenises SMILES strings."""

    def __init__(self, smiles_list: list[str], tokenizer, max_length: int = MAX_LENGTH):
        self.smiles = smiles_list
        self.tokenizer = tokenizer
        self.max_length = max_length

    def __len__(self) -> int:
        return len(self.smiles)

    def __getitem__(self, idx: int):
        enc = self.tokenizer(
            self.smiles[idx],
            max_length=self.max_length,
            padding="max_length",
            truncation=True,
            return_tensors="pt",
        )
        return {
            "input_ids":      enc["input_ids"].squeeze(0),
            "attention_mask": enc["attention_mask"].squeeze(0),
        }


class SMILESRegressionDataset(Dataset):
    """Dataset for regression fine-tuning (SMILES + scalar target)."""

    def __init__(self, smiles_list: list[str], targets: list[float],
                 tokenizer, max_length: int = MAX_LENGTH):
        self.smiles  = smiles_list
        self.targets = targets
        self.tokenizer = tokenizer
        self.max_length = max_length

    def __len__(self) -> int:
        return len(self.smiles)

    def __getitem__(self, idx: int):
        enc = self.tokenizer(
            self.smiles[idx],
            max_length=self.max_length,
            padding="max_length",
            truncation=True,
            return_tensors="pt",
        )
        return {
            "input_ids":      enc["input_ids"].squeeze(0),
            "attention_mask": enc["attention_mask"].squeeze(0),
            "target":         torch.tensor(self.targets[idx], dtype=torch.float32),
        }


# ---------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------

class ChemBERTaRegressor(nn.Module):
    """ChemBERTa backbone + linear regression head.

    Architecture:
        SMILES -> ChemBERTa encoder -> [CLS] embedding -> Dropout -> Linear -> scalar

    When freeze_encoder=True, only the linear head is trained (few-shot regime).
    When freeze_encoder=False, all layers are fine-tuned end-to-end.
    """

    def __init__(self, encoder: nn.Module, hidden_dim: int, freeze_encoder: bool = False):
        super().__init__()
        self.encoder = encoder
        if freeze_encoder:
            for param in self.encoder.parameters():
                param.requires_grad = False
        self.dropout    = nn.Dropout(p=0.1)
        self.regression = nn.Linear(hidden_dim, 1)

    def forward(self, input_ids: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
        out = self.encoder(input_ids=input_ids, attention_mask=attention_mask)
        # Use [CLS] token embedding (first token)
        cls_embedding = out.last_hidden_state[:, 0, :]
        cls_embedding = self.dropout(cls_embedding)
        return self.regression(cls_embedding).squeeze(-1)

    def encode(self, input_ids: torch.Tensor, attention_mask: torch.Tensor) -> torch.Tensor:
        """Return [CLS] embeddings without the regression head."""
        with torch.no_grad():
            out = self.encoder(input_ids=input_ids, attention_mask=attention_mask)
        return out.last_hidden_state[:, 0, :]


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------

def rmse(y_true, y_pred) -> float:
    err = [a - b for a, b in zip(y_true, y_pred)]
    return math.sqrt(sum(e * e for e in err) / len(err))


def mae(y_true, y_pred) -> float:
    return sum(abs(a - b) for a, b in zip(y_true, y_pred)) / len(y_true)


def r2_score(y_true, y_pred) -> float:
    mean_y = sum(y_true) / len(y_true)
    ss_tot = sum((y - mean_y) ** 2 for y in y_true)
    ss_res = sum((y - yhat) ** 2 for y, yhat in zip(y_true, y_pred))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0


# ---------------------------------------------------------------------------
# CSV helpers
# ---------------------------------------------------------------------------

def read_csv(path: str) -> list[dict]:
    """Read CSV into a list of row dicts."""
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Embed mode
# ---------------------------------------------------------------------------

def run_embed(args) -> dict:
    """Extract ChemBERTa [CLS] embeddings for all SMILES in the input CSV."""
    rows = read_csv(args.input)
    all_smiles = [r[SMILES_COLUMN].strip() for r in rows]

    print(f"[chemberta] Loading tokeniser: {args.model_name}", flush=True)
    tokenizer = load_tokenizer(args.model_name)
    print(f"[chemberta] Loading encoder: {args.model_name}", flush=True)
    encoder = load_encoder(args.model_name)
    encoder.eval()

    hidden_dim = encoder.config.hidden_size
    print(f"[chemberta] Encoder loaded. hidden_size={hidden_dim}", flush=True)

    dataset = SMILESDataset(all_smiles, tokenizer, max_length=MAX_LENGTH)
    loader  = DataLoader(dataset, batch_size=args.batch_size, shuffle=False)

    all_embeddings = []
    with torch.no_grad():
        for batch in loader:
            out = encoder(
                input_ids=batch["input_ids"],
                attention_mask=batch["attention_mask"],
            )
            cls_emb = out.last_hidden_state[:, 0, :].cpu().numpy()
            for emb in cls_emb:
                all_embeddings.append(emb.tolist())

    print(f"[chemberta] Embedded {len(all_embeddings)} SMILES.", flush=True)

    predictions = []
    for smi, emb in zip(all_smiles, all_embeddings):
        predictions.append({
            "smiles":    smi,
            "embedding": emb,
            "valid":     True,
        })

    result = {
        "model_name":    args.model_name,
        "num_smiles":    len(all_smiles),
        "num_valid":     len(predictions),
        "embedding_dim": hidden_dim,
        "embeddings":    predictions,
    }
    return result


# ---------------------------------------------------------------------------
# Fine-tune mode
# ---------------------------------------------------------------------------

def run_finetune(args) -> dict:
    """Fine-tune ChemBERTa for regression on the input SMILES + target CSV."""
    rows = read_csv(args.input)

    smiles_list = []
    targets     = []
    for r in rows:
        smi = r.get(SMILES_COLUMN, "").strip()
        tgt = r.get(args.target_col, "").strip()
        if not smi or not tgt:
            continue
        try:
            t = float(tgt)
        except ValueError:
            continue
        if not math.isfinite(t):
            continue
        smiles_list.append(smi)
        targets.append(t)

    n = len(smiles_list)
    print(f"[chemberta] {n} valid molecules for fine-tuning.", flush=True)
    if n < 10:
        raise ValueError(f"Too few valid molecules: {n}")

    # Reproducible train/val/test split (70/15/15)
    torch.manual_seed(args.seed)
    perm     = torch.randperm(n).tolist()
    n_train  = int(0.70 * n)
    n_val    = int(0.15 * n)
    idx_train = perm[:n_train]
    idx_val   = perm[n_train:n_train + n_val]
    idx_test  = perm[n_train + n_val:]

    print(f"[chemberta] Split: train={len(idx_train)}, val={len(idx_val)}, "
          f"test={len(idx_test)}", flush=True)

    print(f"[chemberta] Loading tokeniser: {args.model_name}", flush=True)
    tokenizer = load_tokenizer(args.model_name)
    print(f"[chemberta] Loading encoder: {args.model_name}", flush=True)
    encoder   = load_encoder(args.model_name)
    hidden_dim = encoder.config.hidden_size

    model = ChemBERTaRegressor(encoder, hidden_dim,
                                freeze_encoder=args.freeze_encoder)
    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"[chemberta] Model ready. trainable params={n_params:,}, "
          f"freeze_encoder={args.freeze_encoder}", flush=True)

    def make_subset(indices):
        return SMILESRegressionDataset(
            [smiles_list[i] for i in indices],
            [targets[i] for i in indices],
            tokenizer,
        )

    train_ds = make_subset(idx_train)
    val_ds   = make_subset(idx_val)
    test_ds  = make_subset(idx_test)

    train_loader = DataLoader(train_ds, batch_size=args.batch_size, shuffle=True)
    val_loader   = DataLoader(val_ds,   batch_size=args.batch_size)
    test_loader  = DataLoader(test_ds,  batch_size=args.batch_size)

    optimizer = torch.optim.Adam(
        filter(lambda p: p.requires_grad, model.parameters()),
        lr=args.lr,
        weight_decay=1e-4,
    )
    criterion = nn.MSELoss()

    # Learning rate scheduler: halve LR every lr_step_size epochs
    scheduler = torch.optim.lr_scheduler.StepLR(
        optimizer, step_size=max(1, args.epochs // 3), gamma=0.5
    )

    learning_curves = []

    for epoch in range(1, args.epochs + 1):
        # --- train ---
        model.train()
        train_loss_sum = 0.0
        n_train_batches = 0
        for batch in train_loader:
            optimizer.zero_grad()
            preds = model(batch["input_ids"], batch["attention_mask"])
            loss  = criterion(preds, batch["target"])
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()
            train_loss_sum += loss.item()
            n_train_batches += 1
        train_loss_avg = train_loss_sum / max(1, n_train_batches)

        # --- val ---
        model.eval()
        val_preds, val_true = [], []
        with torch.no_grad():
            for batch in val_loader:
                p = model(batch["input_ids"], batch["attention_mask"])
                val_preds.extend(p.tolist())
                val_true.extend(batch["target"].tolist())
        val_rmse_val = rmse(val_true, val_preds)

        scheduler.step()

        learning_curves.append({
            "epoch":      epoch,
            "train_loss": round(train_loss_avg, 6),
            "val_rmse":   round(val_rmse_val,   6),
        })

        if epoch % max(1, args.epochs // 10) == 0 or epoch == 1:
            print(f"[chemberta] Epoch {epoch:3d}/{args.epochs} "
                  f"train_loss={train_loss_avg:.4f}  "
                  f"val_rmse={val_rmse_val:.4f}", flush=True)

    # --- final evaluation ---
    model.eval()
    all_predictions = []

    def eval_split(loader, split_name, indices_src):
        src_smiles = [smiles_list[i] for i in indices_src]
        src_targets = [targets[i] for i in indices_src]
        preds_out, true_out = [], []
        with torch.no_grad():
            for batch in loader:
                p = model(batch["input_ids"], batch["attention_mask"])
                preds_out.extend(p.tolist())
                true_out.extend(batch["target"].tolist())
        for smi, yt, yp in zip(src_smiles, true_out, preds_out):
            all_predictions.append({
                "smiles": smi,
                "y_true": round(yt, 6),
                "y_pred": round(yp, 6),
                "split":  split_name,
            })
        return true_out, preds_out

    tr_true, tr_pred = eval_split(train_loader, "train", idx_train)
    va_true, va_pred = eval_split(val_loader,   "val",   idx_val)
    te_true, te_pred = eval_split(test_loader,  "test",  idx_test)

    result = {
        "model_name":    args.model_name,
        "target_col":    args.target_col,
        "num_train":     len(idx_train),
        "num_val":       len(idx_val),
        "num_test":      len(idx_test),
        "embedding_dim": hidden_dim,
        "freeze_encoder": args.freeze_encoder,
        "num_epochs":    args.epochs,
        "num_params":    n_params,
        "train_rmse":    round(rmse(tr_true, tr_pred), 6),
        "val_rmse":      round(rmse(va_true, va_pred), 6),
        "test_rmse":     round(rmse(te_true, te_pred), 6),
        "train_mae":     round(mae(tr_true, tr_pred),  6),
        "val_mae":       round(mae(va_true, va_pred),  6),
        "test_mae":      round(mae(te_true, te_pred),  6),
        "test_r2":       round(r2_score(te_true, te_pred), 6),
        "learning_curves": learning_curves,
        "predictions":     all_predictions,
    }
    return result


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="ChemBERTa embedding / fine-tuning")
    parser.add_argument("--mode",         choices=["embed", "finetune"], required=True)
    parser.add_argument("--input",        required=True, help="Input CSV path")
    parser.add_argument("--output",       required=True, help="Output JSON path")
    parser.add_argument("--model-name",   default=DEFAULT_MODEL)
    parser.add_argument("--target-col",   default="ALogP")
    parser.add_argument("--epochs",       type=int,   default=30)
    parser.add_argument("--lr",           type=float, default=2e-5)
    parser.add_argument("--batch-size",   type=int,   default=16)
    parser.add_argument("--freeze-encoder", action="store_true",
                        help="Freeze encoder; train only the regression head")
    parser.add_argument("--seed",         type=int,   default=42)
    args = parser.parse_args()

    if args.mode == "embed":
        result = run_embed(args)
    else:
        result = run_finetune(args)

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.output, "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    print(f"[chemberta] Results written: {args.output}", flush=True)


if __name__ == "__main__":
    main()
