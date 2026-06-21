# Algorithm Guide — emk.dataset

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## 概要

`emk.dataset` モジュールは機械学習ベンチマークデータセット（ESOL, FreeSolv, BBBP, Tox21）を
DeepChem S3 バケットからダウンロードし、`data/benchmark/` にキャッシュして table として返す。
すべての関数は `CacheDir` と `ForceDownload` オプションを持つ。

**共通アルゴリズム**:
1. キャッシュファイルの存在確認 (`exist(cacheFile, "file")`)
2. キャッシュなし or `ForceDownload=true` の場合 → `websave(tmpFile, url)` でダウンロード
3. 成功後 `tmpFile` → `cacheFile` に rename（`movefile`)
4. `readtable(cacheFile, ...)` で CSV を読み込み、列リネーム・型変換

**キャッシュパスの解決**: `which("logInfo")` から `fileparts×3` でプロジェクトルートを特定し `data/benchmark/` を結合する。

---

## `esol`

**データセット**: Delaney ESOL (Estimated SOLubility) — 有機化合物 ~1128 件の水溶性実測値。

| 列名 | 型 | 説明 |
|---|---|---|
| SMILES | string | 化合物 SMILES |
| Name | string | 化合物名 |
| logS | double | 測定 log S（mol/L） |
| logS_Delaney | double | Delaney 式による予測値 |
| MolWt | double | 分子量 |

**引用文献**:
- Delaney, J.S. (2004). ESOL: Estimating Aqueous Solubility Directly from Molecular Structure. *J. Chem. Inf. Comput. Sci.* 44(3):1000-1005. DOI: 10.1021/ci034243x

---

## `freesolv`

**データセット**: FreeSolv — 有機化合物 ~643 件の水和自由エネルギー（実験値・計算値）。

| 列名 | 型 | 説明 |
|---|---|---|
| SMILES | string | 化合物 SMILES |
| Name | string | 化合物名 |
| DeltaG_exp | double | 実験値（kcal/mol） |
| DeltaG_calc | double | GAFF/AM1-BCC 計算値 |
| DeltaG_exp_sem | double | 実験値の標準誤差 |

**引用文献**:
- Mobley, D.L. & Guthrie, J.P. (2014). FreeSolv: A Database of Experimental and Calculated Hydration Free Energies, with Input Files. *J. Comput. Aided Mol. Des.* 28(7):711-720. DOI: 10.1007/s10822-014-9747-x

---

## `bbbp`

**データセット**: BBBP (Blood-Brain Barrier Permeability) — 化合物 ~2050 件の BBB 透過性ラベル。

| 列名 | 型 | 説明 |
|---|---|---|
| SMILES | string | 化合物 SMILES |
| Name | string | 化合物名 |
| BBB | logical | true = BBB 透過 (p), false = 非透過 (n) |

**引用文献**:
- Martins, I.F. et al. (2012). A Bayesian Approach to in Silico Blood-Brain Barrier Penetration Modeling. *J. Chem. Inf. Model.* 52(6):1686-1697. DOI: 10.1021/ci300124c

---

## `tox21`

**データセット**: Tox21 — 化合物 ~7831 件の 12 毒性エンドポイント（マルチラベル）。欠損値は NaN。

**ダウンロードフロー**:
1. `websave(gzFile, url)` — `.csv.gz` ファイルをダウンロード
2. `gunzip(gzFile, cacheDir)` — `.csv` に展開
3. `readtable(..., VariableNamingRule="modify")` — ハイフン含む列名を自動修正
4. 列名マッピング: `"NR-AR"` → `"NR_AR"` 等

**エンドポイント**:

| 列名 | ターゲット |
|---|---|
| NR_AR | Androgen Receptor |
| NR_AR_LBD | Androgen Receptor LBD |
| NR_AhR | Aryl Hydrocarbon Receptor |
| NR_Aromatase | Aromatase |
| NR_ER | Estrogen Receptor Alpha |
| NR_ER_LBD | Estrogen Receptor LBD |
| NR_PPAR_gamma | PPAR-gamma |
| SR_ARE | Antioxidant Response Element |
| SR_ATAD5 | Genotoxicity (ATAD5) |
| SR_HSE | Heat Shock Response |
| SR_MMP | Mitochondrial Membrane Potential |
| SR_p53 | p53 Tumor Suppressor |

**引用文献**:
- Tox21 Challenge: https://tripod.nih.gov/tox21/challenge/
- Wu, Z. et al. (2018). MoleculeNet: A Benchmark for Molecular Machine Learning. *Chem. Sci.* 9(2):513-530. DOI: 10.1039/C7SC02664A

**テスト戦略** (`tests/unit/TestDataset.m`):

| カテゴリ | 内容 |
|---|---|
| TC-1: バリデーション | 無効 CacheDir 型・不明オプション（ネットワーク不要） |
| TC-2: ネットワーク | スキーマ確認・行数 bounds・キャッシュ再使用・Tox21 エンドポイント値 |
