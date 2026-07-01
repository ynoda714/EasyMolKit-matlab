# RP05: BBBP ECFP4 に対する SHAP 再現

[English README](./README.md)

この repro では、BBBP 分類を対象とした ECFP4 指紋の SHAP 解釈ワークフローを、EasyMolKit と MATLAB でどこまで再現できるかを検証します。
主たる検証済み成果は F2 の解析的 linear-SHAP 再現であり、F3 は MATLAB `shapley()` の実行時間上の実現可能性検証として分離しています。

---

## Overview

| 項目 | 内容 |
|---|---|
| 論文 | Rodriguez-Perez, R. & Bajorath, J. (2020). Interpretation of Machine Learning Models Using Shapley Values: Application to Compound Potency and Multi-target Activity Predictions. *J. Comput.-Aided Mol. Des.* 34:1013-1026. |
| DOI | [10.1007/s10822-020-00314-0](https://doi.org/10.1007/s10822-020-00314-0) |
| タスク | BBBP 分類に対する SHAP ベース解釈 |
| モデル | Logistic Regression, Random Forest, SHAP / TreeSHAP |
| データ | MoleculeNet / DeepChem BBBP（有効分子 `2,039` 件） |
| 要点 | 正則化設定を揃えれば、MATLAB は linear-SHAP 経路を良好に再現できる。一方で tree-SHAP は現行 MATLAB ワークフローでは実行時間制約が大きい |

---

## この Reproduction が扱う範囲

RP05 は次の 3 層からなる複合 repro です。

- Python による logistic-regression SHAP を基準経路として扱う
- F2 の解析的 MATLAB linear-SHAP 再現を主たる検証対象とする
- F3 の MATLAB `shapley()` と Python TreeSHAP の比較を、制約付きの実現可能性検証として扱う

したがって本 RP は、公開データセット上での再現・比較研究として読むべきものであり、MATLAB が Python 側のあらゆる SHAP 実装を完全再現したと主張するものではありません。

F3 の詳細なベンチマーク設計、集計方針、Zone C / D の解釈は [F3_spec.md](./F3_spec.md) に分離しています。
ただし現行実装では、F3 は `exploratory` モードで `n_eval=[1, 2, 4]`、`n_repeats=1` のみを対象としており、spec にあるより大きな推奨グリッドを既定では実行しません。

---

## Environment

実行時の実バージョンは、実行後に `result/runs/<timestamp>/lock_snapshot.json` へ保存されます。

| 項目 | 要件 |
|---|---|
| MATLAB | R2025a 以降 |
| Python | EasyMolKit Embedded Python |
| Python パッケージ | `shap>=0.49.1`, `scikit-learn>=1.7` |
| ライブラリ | RDKit, NumPy, SciPy, pandas |
| Toolbox | Statistics and Machine Learning Toolbox |

### Descriptor Definitions

| 記述子 | ツール | 定義 |
|---|---|---|
| Morgan ECFP4 | RDKit | `radius=2`, `nBits=2048` の Morgan fingerprint |

---

## Data

- **Source**: MoleculeNet / DeepChem BBBP
- **URL**: 初回利用時に EasyMolKit のデータ取得フローを通じて取得
- **License**: 元データ提供元の条件に従う
- **Cache path**: `data/benchmark/bbbp.csv`
- **本 RP での使用件数**: SMILES 妥当性確認後の有効分子 `2,039` 件
- **SHAP 比較に使う split**: stratified 80/20 split（`train=1,631`, `test=408`, `seed=42`）

---

## How to Run

プロジェクトルートで MATLAB を開き、次を実行します。

```matlab
cd repro/rp05_shap
edit rp05_shap.m
```

その後、スクリプトをセクションごとに実行します。

前提条件:

- EasyMolKit の依存関係が導入済みである
- `emk.setup.initPython()` が手元の環境で動作する
- 設定済み Python 環境に RDKit と `shap` が入っている

主要ファイル:

```text
repro/rp05_shap/rp05_shap.m
repro/rp05_shap/rp05_shap_core.py
repro/rp05_shap/F3_spec.md
repro/rp05_shap/lock_template.json
```

大まかな流れ:

| ファイル / Section | 内容 |
|---|---|
| `rp05_shap.m` Section 0 | セットアップ、Python 初期化、環境記録 |
| `rp05_shap.m` Section 1 | BBBP データセット解決と補助チェック |
| `rp05_shap.m` Section 2 | Python logistic-regression SHAP 基準実行 |
| `rp05_shap.m` Section 2b | 改訂版 (`task_b`) LR + SHAP 比較 |
| `rp05_shap.m` Section 2c | F2 の解析的 MATLAB linear-SHAP 再現 |
| `rp05_shap.m` Section 2d | F3 の MATLAB `shapley()` ベンチマークと早期停止ガードレール |
| `rp05_shap_core.py` | BBBP 読み込み、ECFP4 生成、Python SHAP 出力 |
| `F3_spec.md` | F3 の詳細ベンチマーク条件と解釈ルール |

---

## Result

代表的な検証済み実行: `result/runs/20260629_160127_rp05_shap`（`2026-06-29`）

| 指標 | 値 | 注記 |
|---|---|---|
| `auc_cv` | `0.9096 +/- 0.0083` | ベースライン LR の交差検証 |
| `shap_lr_spearman` | `0.8538` | Python `LinearExplainer` と LR 重み由来順位の比較 |
| `task_b.auc_cv_rev` | `0.9143 +/- 0.0079` | 改訂版 LR 比較経路 |
| `task_b.ranking_spearman_rho` | `0.9149` | baseline と改訂版の順位比較 |
| `task_b.top20_overlap` | `14 / 20` | 上位 fingerprint bit の重なり |

### F2: 検証済み成果

| 比較 | Spearman rho |
|---|---:|
| MATLAB 解析的 SHAP vs Python LR SHAP (`C=1.0`) | `0.9274` |
| MATLAB 解析的 SHAP vs Python LR SHAP (`C=0.10`) | `0.8355` |

正則化を一致させたケースでは、MATLAB の解析的定式化が Python 側 linear-SHAP の順位を実用上十分な水準で再現できることを示しています。
この状態は `metrics.json` に `f2_matlab_shap.zone_b_confirmed = true` として記録されます。

### F3: 実行時間制約付き検証

| 項目 | 値 |
|---|---:|
| F3 設定モード | `exploratory` |
| 設定 `n_eval` グリッド | `[1, 2, 4]` |
| 設定 `n_repeats` | `1` |
| ベンチマーク `n_eval` | `1` |
| ベンチマーク実行時間 | `469.7 s` |
| 全体実行時間の予測値 | `54.8 min` |
| ベンチマーク rho | `0.5823` |
| MATLAB 手法 | `interventional-tree` |
| 状態 | `SKIPPED` |

上表の `n_eval=1` は、設定済み exploratory 実行の所要時間を見積もるためのガードレール用ベンチマークであり、spec 推奨のより大きなグリッドを実行したことを意味しません。
現行コード経路では、F3 は `n_eval=[1, 2, 4]`・単回反復の exploratory 設定として扱われ、その exploratory グリッドへ進む前にまずベンチマークを実行します。
F3 はまずベンチマークを実行し、全量実行時間の予測値が設定ガードレールを超える場合は早期終了します。
代表ランでは、予測実行時間が `20.0` 分を超えたため、意図通り停止しました。

---

## Verification

この RP では、RF03 の主たる判定として Cat A と Cat C を用います。
Cat B は baseline と改訂版の比較に関する参考情報として保持します。

| 指標 | 受け入れ条件 |
|---|---|
| `auc_cv` | `>= 0.85` |
| `task_b.auc_cv_rev` | `>= 0.85` |
| `shap_lr_spearman` | `>= 0.85` |
| `f2_matlab_shap.zone_b_confirmed` | `true` |

解釈:

- LR ベースラインは最小 AUC 閾値を十分に上回っている
- F2 経路が MATLAB 側の主たる検証済み再現結果である
- F3 は現行ガードレール下では、まだ zone 判定可能な完了研究には至っていない
- 現行の exploratory グリッド `n_eval=[1, 2, 4]` 自体も、spec 側の formal な Zone 判定最小要件である `n_eval>=64` を満たしていない

---

## Discussion

### 論文との差分

| 差分 | 詳細 |
|---|---|
| 基準の置き方 | 本 RP では Python SHAP を運用上の基準とし、MATLAB がどこまで近づけるかを評価している |
| F2 の範囲 | linear-SHAP は MATLAB から Python explainer を直接呼ばず、解析的に再現している |
| F3 の範囲 | tree 系 SHAP は完全一致主張ではなく、実現可能性に制約のある比較実験として扱っている |

### 実務上の読み取り

- Python 参照に対する MATLAB 側 linear-SHAP ワークフローの妥当性確認が目的なら、RP05 はすでに有用です
- 2048 bit BBBP 特徴量に対する完全な tree-SHAP 再現を MATLAB 単独で行う目的では、現状ワークフローはまだ実行時間制約が強いです
- F3 の再設計案やベンチマーク方針の詳細は、トップ README ではなく [F3_spec.md](./F3_spec.md) に置いています

---

## Files

| ファイル | 内容 |
|---|---|
| `README.md` | 英語版の正本 README |
| `README.jp.md` | 日本語版 README |
| `rp05_shap.m` | 主たる MATLAB 再現スクリプト |
| `rp05_shap_core.py` | Python 基準ワークフローと MAT 出力 |
| `F3_spec.md` | F3 詳細仕様 |
| `lock_template.json` | 実行時メタデータ用テンプレート |

スクリプト実行後には `result/runs/<ts>/` 以下に、metrics、status records、environment snapshots などのローカル出力も生成されます。
