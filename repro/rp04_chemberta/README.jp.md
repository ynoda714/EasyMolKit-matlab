# RP04: ChemBERTa BBBP 線形プローブ

[English README](./README.md)

この repro は、凍結した ChemBERTa 分子言語モデルが BBBP 分類でどの程度強いベースラインになるかを EasyMolKit 上で評価するものです。
主たる検証対象は、Python で ChemBERTa 埋め込みを抽出してロジスティック回帰をかける標準 RP04 経路であり、加えて MATLAB 側への移管を進めた Zone C 派生経路も収録しています。

---

## Overview

| 項目 | 内容 |
|---|---|
| 論文 | Chithrananda, S. et al. (2020). ChemBERTa: Large-Scale Self-Supervised Pretraining for Molecular Property Prediction. *arXiv*:2010.09885. |
| DOI | [10.48550/arXiv.2010.09885](https://doi.org/10.48550/arXiv.2010.09885) |
| タスク | BBB 透過性分類（`BBB+` / `BBB-`） |
| モデル | `seyonec/ChemBERTa-zinc-base-v1` に対するロジスティック回帰線形プローブ |
| データ | MoleculeNet / DeepChem BBBP（有効分子 `2,039` 件） |
| 要点 | 凍結 ChemBERTa の CLS 埋め込みは BBBP で強いベースラインとなり、代表ランでは公正基準の RP02 / RP03 を上回る `AUC ≈ 0.927` を示す |

---

## この Reproduction が扱う範囲

RP04 は 3 層構成の composite repro です。

- 主経路: Python で ChemBERTa 埋め込みを抽出し、標準 RP04 ワークフローでロジスティック回帰を行う
- F1-b: Python で埋め込みまで作り、MATLAB `fitclinear` で分類する
- F1-a: Python は tokenization のみ、MATLAB 側で ONNX 推論と `fitclinear` を行う

したがってこれは、公開データセット上での再現および転移学習比較研究として読むべきものであり、
元論文の split、モデル変種、下流学習法をそのまま完全復元したものではありません。

---

## Environment

実行時の実バージョンは `result/runs/<timestamp>/lock_snapshot.json` に保存されます。

| 項目 | 要件 |
|---|---|
| MATLAB | R2025a 以降 |
| Python | 3.10（EasyMolKit Embedded Python） |
| RDKit | 2022.03 以降 |
| transformers | 5.0 以降 |
| scikit-learn | 1.0 以降 |
| onnx | F1-a の ONNX export / patch helper に Python 側で必要 |
| Toolbox | 主経路の RP04 は EasyMolKit の Python 連携以外は追加不要。F1-b では Statistics and Machine Learning Toolbox（`fitclinear`, `perfcurve`, `cvpartition`）が必要。F1-a ではこれに加えて Deep Learning Toolbox（`importNetworkFromONNX`）が必要 |

### Feature Definition

| 特徴量 | ツール | 定義 |
|---|---|---|
| ChemBERTa CLS 埋め込み | `seyonec/ChemBERTa-zinc-base-v1` | SMILES を tokenization した後、凍結 RoBERTa 系エンコーダ最終層の `[CLS]` 表現を使用 |

---

## Data

- **Source**: MoleculeNet / DeepChem BBBP
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/BBBP.csv`
- **License**: 研究・教育利用向け MoleculeNet ベンチマークデータ
- **Cache path**: `data/benchmark/bbbp.csv`
- **本 RP での使用件数**: `2,039` 件（元の `2,050` 行から有効 SMILES を使用）
- **Data hash**: `result/runs/<ts>/lock_snapshot.json` の `dataset_sha256` を参照

---

## How to Run

MATLAB をプロジェクトルートで開き、主スクリプトを開きます。

```matlab
cd repro/rp04_chemberta
edit rp04_chemberta.m
```

その後、スクリプトをセクションごとに実行してください。

前提条件:

- EasyMolKit の依存関係が導入済みである
- `emk.setup.initPython()` がローカル環境で動作する
- 設定済み Python 環境に RDKit、`transformers`、`scikit-learn` が入っている
- F1-a を実行する場合は、Python の `onnx` パッケージも追加で必要
- ChemBERTa モデルがローカルにあるか、Hugging Face から取得可能である

主なファイル:

```text
repro/rp04_chemberta/rp04_chemberta.m
repro/rp04_chemberta/rp04_chemberta_core.py
repro/rp04_chemberta/rp04_chemberta_f1.m
repro/rp04_chemberta/rp04_chemberta_f1a.m
repro/rp04_chemberta/lock_template.json
```

`rp04_chemberta.m` の大まかな流れ:

| Section | 内容 |
|---|---|
| Section 0 | セットアップと環境記録 |
| Section 1 | BBBP データ読み込みとパス解決 |
| Section 2 | 公正基準の RP02 / RP03 を読み込む |
| Section 3 | ChemBERTa CLS 埋め込み抽出と 5-fold CV |
| Section 4 | トークン長検証と `max_length=128` / `512` 比較 |
| Section 5 | 公正基準との比較要約 |
| Section 6 | 可視化図を保存 |
| Section 7 | RF03 検証 |
| Section 8 | metrics と lock snapshot を保存 |

---

## Result

代表的な主経路の検証済みラン: `2026-06-26`（`max_length=512`）
この表の比較値の出典: `auc_rp02_rev` と `auc_rp03_rev` は、この代表ラン実行時に利用可能だった最新の公正 upstream artifact から動的に読み込まれた値です。upstream run が見つからない場合だけ使うハードコード fallback 定数とは別です。

| 指標 | 値 | 注記 |
|---|---|---|
| `auc_cv` | `0.9270 +- 0.0107` | 主たる RP04 ChemBERTa 結果 |
| `auc_rp02_rev` | `0.9143 +- 0.0079` | 代表ランで読み込まれた公正基準 LR+ECFP4 |
| `auc_rp03_rev` | `0.9038 +- 0.0203` | 公正基準 leak 修正済み GCN |
| `delta_vs_rp02_rev` | `+0.0128` | 主比較指標 |
| `delta_vs_rp03_rev` | `+0.0232` | 主比較指標 |

トークン長検証の要約:

| 設定 | 切り捨て分子数 | AUC CV |
|---|---:|---:|
| `max_length=128` | `26 / 2039`（`1.3%`） | `0.9271` |
| `max_length=512` | `0 / 2039`（`0.0%`） | `0.9270` |

`max_length=128` と `512` の差は `delta_auc = -0.0001` と無視できる範囲であり、
当初の切り捨て問題は開示したうえで定量的に上限を確認済みです。

### F1 派生の現状

| 派生 | Python の役割 | MATLAB の役割 | AUC CV | 状態 |
|---|---|---|---:|---|
| F1-b（`2026-06-28`） | tokenization + embedding | `fitclinear` LR | `0.9138 +- 0.0088` | RF03 PASS |
| F1-a（`2026-06-28`, ヘッダ記録のみ） | tokenization のみ | ONNX 推論 + `fitclinear` LR | `0.9138 +- 0.0087` | RF03 PASS。ただし対応する保存 artifact は現ワークスペースに残っていない |

現ワークスペースでは F1-b の `metrics.json` は保持されていますが、F1-a の保持済み artifact は見当たりません。
F1-a については、script ヘッダに直近の成功 run（`2026-06-28`）の結果が記録されていますが、その run の `metrics.json` / `lock_snapshot.json` は現ワークスペースに残っていません。
F1-a script 自体は、今後の成功 run で `metrics.json` と `lock_snapshot.json` を出力する状態になっています。

RP02 / RP03 の比較基準は最新の公正 run から動的に読み込むため、上流の基準 run が更新されると比較欄の値も変わり得ます。
上流の公正 run が見つからない場合に限り、`rp04_chemberta.m` は contingency としてハードコード定数（`RP02=0.9118`, `RP03=0.9038`）を使います。

---

## Verification

この RP は主として Cat A の絶対閾値で検証します。
実務上の受け入れ条件は次の 1 点です。

| 指標 | 受け入れ条件 |
|---|---|
| `ROC-AUC CV` | `>= 0.85` |

根拠:

- 元論文は scaffold split を報告しており、本 RP の random 5-fold CV より厳しい
- 本 RP は full fine-tuning ではなく、凍結埋め込みに対する線形プローブである
- 使用チェックポイントは公開 ChemBERTa 変種であり、論文内の全設定と完全一致ではない

代表主経路ランは RF03 を満たしています。
F1-b も保持済み artifact で RF03 PASS を確認できます。
F1-a も script ヘッダ上は `2026-06-28` の成功 run が同じ閾値を上回っていますが、その run の保存 artifact は現ワークスペースに残っていません。

---

## Discussion

### 論文との主な差分

| 差分 | 詳細 |
|---|---|
| モデル変種 | 論文では 5M / 10M 系の ChemBERTa が議論される一方、本 RP は `seyonec/ChemBERTa-zinc-base-v1` を使う |
| 下流学習法 | 論文は fine-tuning を主眼とするが、本 RP は凍結埋め込み + ロジスティック回帰を使う |
| split | 論文は scaffold split、本 RP は RP02 / RP03 との比較を優先して random 5-fold CV を使う |
| 評価の主眼 | 元論文表の再掲よりも、EasyMolKit 内の改訂済み公正基準との比較を重視する |

### 実務上の読み取り

- RP04 は、凍結 ChemBERTa 表現だけでもこの EasyMolKit 条件下で BBBP に十分競争力があることを示しています
- 代表的な公正比較では、ChemBERTa は改訂済み LR+ECFP4 と改訂済み GCN を平均 AUC で上回ります
- F1 派生により、下流ワークフローのより大きな部分を MATLAB 側に移せることも確認していますが、公開向け主結果は標準 RP04 経路を中心に読むのが適切です

---

## Files

| ファイル | 内容 |
|---|---|
| `README.md` | 英語 canonical README |
| `README.jp.md` | 日本語 companion README |
| `rp04_chemberta.m` | 主たる ChemBERTa 再現スクリプト |
| `rp04_chemberta_core.py` | Python 側の埋め込み、tokenization、ONNX 補助コード |
| `rp04_chemberta_f1.m` | F1-b Zone C 派生 |
| `rp04_chemberta_f1a.m` | F1-a ONNX ベース Zone C 派生 |
| `lock_template.json` | 実行時メタデータのテンプレート |
| `metrics_biased_historical.json` | `result/runs/<ts>/` に保存される historical biased 値。監査経緯を失わないためのもので、解析や比較には使わない |

スクリプト実行により `result/runs/<ts>/` 配下に metrics、図、environment snapshot などのローカル成果物も生成されます。
