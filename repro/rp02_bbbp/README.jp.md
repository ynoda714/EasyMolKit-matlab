# RP02: MoleculeNet BBBP 分類ベースライン

[English README](./README.md)

この repro は、Wu et al. (2018) の MoleculeNet BBBP タスクに対する EasyMolKit の公開ベースラインを構築するものです。
正本となる経路は、Morgan ECFP4 とロジスティック回帰を nested cross-validation で評価する sklearn ベースラインです。
あわせて、なぜ現行 sklearn ベースラインが旧 MATLAB `fitclinear` ベースラインを置き換えたのかを説明する診断実行も保持しています。

---

## Overview

| 項目 | 内容 |
|---|---|
| 論文 | Wu, Z. et al. (2018). MoleculeNet: A Benchmark for Molecular Machine Learning. *Chem. Sci.* 9:513-530. |
| DOI | [10.1039/C7SC02664A](https://doi.org/10.1039/C7SC02664A) |
| タスク | 血液脳関門透過性の二値分類 (`BBB+` vs `BBB-`) |
| モデル | Morgan ECFP4 フィンガープリント + ロジスティック回帰 (L2) |
| データ | MoleculeNet BBBP (`2,039` 有効分子) |
| 要点 | 代表検証 run `20260701_210031_rp02_bbbp` では、`rp02_bbbp.m` Section 2 / `metrics.json` に基づく nested-CV ROC-AUC が `0.9143 +/- 0.0089`。条件は `outer_seed=42`, `inner_seed=7` |

---

## この Reproduction が扱う範囲

主スクリプトは 1 つの公開ベースラインと 3 つの補助診断を提供します。

- 正本ベースライン: sklearn `LogisticRegression` による outer 5-fold stratified CV + inner 3-fold `C` 選択
- 履歴ベースライン監査: 旧 MATLAB `fitclinear` 経路との比較
- MATLAB parity study: MATLAB `lbfgs` が独自の最適正則化を使うと sklearn にどこまで近づくかの診断
- scaffold split 診断: Wu et al. との差のうち split 戦略がどの程度効いているかの確認

GitHub 向けの主結果は sklearn nested-CV ベースラインです。
診断用スクリプトは再現性と監査可能性のために同じディレクトリに保持しています。

---

## Environment

実際の実行環境は `result/runs/<timestamp>/lock_snapshot.json` に保存されます。
以下の代表値は `result/runs/20260701_210031_rp02_bbbp/lock_snapshot.json` を参照しています。

| 項目 | 要件 |
|---|---|
| MATLAB | R2025a 以降 |
| Python | 3.10 (EasyMolKit 組み込み Python) |
| RDKit | 代表 run では 2024.03.6 以降 |
| Toolbox | Statistics and Machine Learning Toolbox |

### Fingerprint Definition

| 特徴量 | ツール | 定義 |
|---|---|---|
| Morgan ECFP4 | RDKit `rdFingerprintGenerator.GetMorganGenerator` | radius `2`, `2048` bit。radius 2 は ECFP4 に対応 |

---

## Data

- **Source**: DeepChem / MoleculeNet BBBP (`BBBP.csv`)
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv`
- **License**: DeepChem 配布物。元データの由来は upstream dataset と Martins et al. (2012) を参照
- **Cache path**: `data/benchmark/bbbp.csv`
- **本 RP での使用件数**: `2,039` 有効分子 (`metrics.json -> n_valid`)
- **除外行**: SMILES parse failure `11` 件。総数はローダーから再計算できますが、ラベル別内訳は現在の保存 artifact には出力していません

---

## How to Run

プロジェクトルートで MATLAB を開き、次を実行します。

```matlab
cd repro/rp02_bbbp
edit rp02_bbbp.m
```

その後、セクションごとに実行してください。

前提条件:

- EasyMolKit の依存関係がインストール済みであること
- `emk.setup.initPython()` がローカル環境で動作すること
- 設定済み Python 環境で RDKit が利用できること

主なスクリプト:

```text
repro/rp02_bbbp/rp02_bbbp.m
repro/rp02_bbbp/rp02_sklearn_core.py
repro/rp02_bbbp/a1_diagnosis_run.m
repro/rp02_bbbp/a1_diagnosis.py
repro/rp02_bbbp/r1b_solver.m
repro/rp02_bbbp/r1c_matlab_nested_cv.m
```

高水準フロー:

| Section | 内容 |
|---|---|
| Section 0 | セットアップと環境情報取得 |
| Section 1 | BBBP データセットのキャッシュと検証 |
| Section 2 | `rp02_sklearn_core.py` による sklearn nested CV |
| Section 2b | 同一 fold 上での MATLAB `fitclinear` 再実行による診断 |
| Section 2c | Wu et al. 比較用の scaffold-based nested CV |
| Section 3 | RF03 検証と比較ログ |
| Section 4 | pseudo-ROC artifact 保存 |
| Section 5 | metrics、fold index、version lock の保存 |

---

## Result

代表的な検証済みベースライン run: `result/runs/20260701_210031_rp02_bbbp` (`2026-07-01`)

| 指標 | 値 | 出典 |
|---|---|---|
| ROC-AUC（nested CV, outer 5-fold 平均） | `0.9143 +/- 0.0089` | `rp02_bbbp.m` Section 2、`metrics.json -> auc_cv / auc_cv_std`、`rp02_sklearn_core.py` の `outer_seed=42`, `inner_seed=7` |
| Pseudo-ROC AUC | `0.9130` | `rp02_bbbp.m` Section 4、`metrics.json -> auc_pseudo_roc`。記述用のみ |
| Scaffold 5-fold ROC-AUC | `0.8832 +/- 0.0173` | `rp02_bbbp.m` Section 2c、`metrics.json -> scaffold_cv.auc_cv / auc_cv_std` |
| 旧 MATLAB `fitclinear` ベースライン | `0.8826 +/- 0.0220` | `a1_diagnosis_run.m`、artifact `result/runs/20260621_204831_a1_auc_gap/a1_diagnosis.json -> auc_rp02_matlab` と fold spread |
| MATLAB `lbfgs` 診断の最良値 | `0.9084 +/- 0.0095` | `r1c_matlab_nested_cv.m` Section 3b の fair nested CV、artifact `result/runs/20260626_103836_r1c_nested/r1c_results.json -> section3b` |

代表的な実行環境（`lock_snapshot.json` より）:

| 項目 | 値 |
|---|---|
| MATLAB | `R2026a` |
| Python | `3.10` |
| RDKit | `2024.03.6` |
| Commit | `54ca766` |

注記:

- 現在の公開ベースラインは `rp02_bbbp.m` が使う `inner_seed=7` の経路であり、README の主数値は `20260701_210031_rp02_bbbp` に対応します。
- `rp02_sklearn_core.py` には、`inner_seed=42` を使った履歴上の certified run (`AUC=0.9118`) も記録されています。現行 README の主値と混同しないでください。
- `rp02_bbbp.m` Section 3 の `GCN_AUC_RP03 = 0.9151` は RP03 から手動転記した比較定数です。RP03 を再実行した場合は、この値を更新してから比較文を再利用する必要があります。
- `r1c_matlab_nested_cv.m` Section 3a の fixed-Lambda sweep は outer test fold の AUC を直接見て最良 Lambda を選ぶ参考実験です。Section 3b と同値になった run であっても、3a 単独を不偏な汎化性能推定として扱うべきではありません。

---

## Verification

この RP では、公開ベースラインに対して RF03 Cat A の単一基準を使います。

| 指標 | 実務上の目標 |
|---|---|
| Nested-CV ROC-AUC | `>= 0.85` |

この閾値が妥当な理由:

- 公開 RP02 ベースラインは random stratified nested CV を使っており、Wu et al. の固定 scaffold split より容易です
- 代表的な sklearn ベースラインは `0.91` 台前半で安定しているため、`0.85` は通常の環境差分を吸収しつつ緩すぎない基準です
- PASS/FAIL には pseudo-ROC AUC を使わず、nested-CV 平均のみを正とします

---

## Discussion

### Wu et al. との主な差分

| 差分 | 内容 |
|---|---|
| Split 戦略 | Wu et al. は scaffold split の test 結果を報告。RP02 は random nested CV を正本とし、scaffold CV は診断扱い |
| 評価形式 | Wu et al. は固定 test-set 結果、RP02 は cross-validation の安定性を重視 |
| ハイパーパラメータ | RP02 は sklearn `C` を inner CV で選択し、固定デフォルトに依存しない |
| フィンガープリント実装 | RP02 は RDKit の現行 Morgan fingerprint generator API を利用 |

### なぜ旧 MATLAB ベースラインを置き換えたか

- 履歴 A1 監査では、元の `0.8826` の差を固定 `C=1.0` 条件で分解しており、`solver_contribution = +0.0191`、`fold_split_contribution = +0.0079` でした。出典は `result/runs/20260621_204831_a1_auc_gap/a1_diagnosis.json` です。
- `fold_split_contribution` が正ということは、この監査条件では MATLAB `cvpartition` fold よりも sklearn 側の stratified fold の方が平均 AUC を約 `0.0079` 押し上げたことを意味します。
- MATLAB 側のフェアな parity check は `r1c_matlab_nested_cv.m` Section 3b で、MATLAB `lbfgs` が inner CV で独自に `Lambda` を選ぶと `0.9084 +/- 0.0095` に到達し、sklearn との差は `+0.0058` まで縮みます。
- `r1c_matlab_nested_cv.m` Section 3a の fixed-Lambda sweep は outer test fold を見ながら最良 Lambda を選ぶため、方法論的には optimistic な参考値です。Section 3b の nested CV 結果を主たる比較値として扱ってください。
- 用語は意図的に分けています。A1 の `solver_contribution` は固定 `C` の履歴監査、`rp02_bbbp.m` Section 2b の `solver_gap` / `regularization_gap` は別定義の matched-fold 診断です。同じ意味ではありません。

### 実務上の使い分け

後続 RP の比較用に RP02 ベースラインが必要なら、`20260701_210031_rp02_bbbp` の sklearn nested-CV ROC-AUC `0.9143 +/- 0.0089` を使ってください。
ソルバー感度、fold 感度、scaffold split 差分を見たい場合は、主ベンチマークと混同せず、このディレクトリの診断 artifact を参照してください。

---

## Files

| ファイル | 内容 |
|---|---|
| `README.md` | 英語版正本 README |
| `README.jp.md` | 日本語版 README |
| `rp02_bbbp.m` | 主 reproduction スクリプト |
| `rp02_sklearn_core.py` | Python nested-CV コア |
| `r1b_solver.m` | ソルバー比較診断 |
| `r1c_matlab_nested_cv.m` | MATLAB nested-CV 診断 |
| `a1_diagnosis_run.m` / `a1_diagnosis.py` | 履歴 A1 監査スクリプト |
| `lock_template.json` | 実行メタデータテンプレート |

スクリプト実行時には `result/runs/<ts>/` 以下に、metrics、pseudo-ROC artifact、fold index、診断 JSON、環境 snapshot などのローカル出力も作成されます。
