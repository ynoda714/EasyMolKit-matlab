# RP00: ESOL 水溶性予測

[English README](./README.en.md)

この repro は、Delaney (2004) の ESOL 線形モデルを EasyMolKit と MATLAB で再現し、
公開データセット上でどの程度まで同等の挙動を確認できるかを示すベースラインです。
ESOL 系 RP の出発点として、4 つの分子記述子だけでどこまで再現できるかを確認する位置づけです。

---

## Overview

| 項目 | 内容 |
|---|---|
| 論文 | Delaney, J.S. (2004). ESOL: Estimating Aqueous Solubility Directly from Molecular Structure. *J. Chem. Inf. Comput. Sci.* 44(3):1000-1005. |
| DOI | [10.1021/ci034243x](https://doi.org/10.1021/ci034243x) |
| タスク | 水溶性予測（`logS`, mol/L） |
| モデル | 4 分子記述子による線形回帰 |
| データ | MoleculeNet ESOL（総数 `1,128` 分子。filter 後の回帰件数はこれより少なくなり得る） |
| 要点 | 古典的 ESOL ベースラインに近い 5-fold CV 精度は再現できるが、記述子定義の違いにより係数は完全一致しない |

---

## この Reproduction が扱う範囲

このスクリプトは、次の 4 記述子を使った標準的な ESOL ワークフローを再現します。

- `LogP`
- `MolWt`
- `NumRotatableBonds`
- `AromaticProportion`

また、元論文との実務上の差分として次の 2 点を明示します。

- Delaney の `clogP` と RDKit の `MolLogP` は同一ではない
- 公開されている MoleculeNet / DeepChem ESOL 表は、論文記載の `1,144` 件ではなく `1,128` 件である

そのため、この RP は「ESOL の考え方を公開データ上で再現したもの」であり、
元論文の実験環境を完全にそのまま復元したものではありません。

---

## Environment

実行時の実バージョンは `result/runs/<timestamp>/lock_snapshot.json` に保存されます。

| 項目 | 要件 |
|---|---|
| MATLAB | R2025a 以降 |
| Python | 3.10（EasyMolKit Embedded Python） |
| RDKit | 2024.03.6 以降 |
| Toolbox | Statistics and Machine Learning Toolbox |

### Descriptor Definitions

| 記述子 | ツール | 定義 |
|---|---|---|
| LogP | RDKit `Descriptors.MolLogP` | Crippen-Wildman `MolLogP`。Delaney の元論文で使われた `clogP` と同一ではない |
| MolWt | RDKit `Descriptors.MolWt` | 暗黙水素を含む平均分子量（IUPAC atomic weights） |
| NumRotatableBonds | RDKit `rdMolDescriptors.CalcNumRotatableBonds` | strict SMARTS に基づく回転可能結合数 |
| AromaticProportion | `pyrun` batch | 芳香族原子数を heavy-atom 数で割った比率 |

---

## Data

- **Source**: DeepChem / MoleculeNet `delaney-processed.csv`
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv`
- **License**: Delaney (2004) 由来の公開データ
- **Cache path**: `data/benchmark/esol.csv`
- **元データ総数**: `1,128` 分子
- **実際の回帰対象件数**: SMILES パース失敗や descriptor NaN 除外後の run 依存値。`metrics.json` の `n_molecules` と `excluded.csv` を参照

---

## How to Run

MATLAB でプロジェクトルートから次を開きます。

```matlab
cd repro/rp00_esol
edit rp00_esol_pilot.m
```

その後、セクションごとに実行するか、環境初期化が済んでいればスクリプト全体を実行してください。

前提条件:

- EasyMolKit の依存関係が導入済みである
- `emk.setup.initPython()` がローカル環境で正常に動く
- 設定済み Python 環境に RDKit が入っている

対象スクリプト:

```text
repro/rp00_esol/rp00_esol_pilot.m
```

処理の流れ:

| Section | 内容 |
|---|---|
| Section 0 | セットアップと環境情報取得 |
| Section 1 | ESOL データセット読込とデータハッシュ計算 |
| Section 2 | SMILES パースと記述子計算 |
| Section 3 | 全分子で線形回帰を学習 |
| Section 4 | 5-fold 交差検証 |
| Section 5 | 受入指標の要約 |
| Section 6 | 実行成果物を保存 |

---

## Result

代表実行結果（`2026-06-19`）:

| 指標 | 値 | 注記 |
|---|---|---|
| RMSE（全データ fit） | `1.0116` | 記述用 |
| R^2（全データ fit） | `0.7680` | 記述用 |
| RMSE（5-fold CV, pooled OOF） | 約 `1.02` | RF03 の主 RMSE 判定指標 |
| RMSE（5-fold CV, fold 平均） | `1.0166 +/- 0.0243` | fold 間の参考要約 |
| R^2（5-fold CV, pooled OOF） | 約 `0.76` | 主たる CV 指標 |

解釈:

- 再現誤差は古典的な ESOL ベースラインに近い
- `LogP` の定義差があるため、元論文ほど高い `R^2` にはならない
- 係数の大きさはずれても、主要な傾向の向きは概ね一致する

代表実行での係数例:

```text
(Intercept)  +0.255
LogP         -0.745
MolWt        -0.0065
RotBonds     +0.0026
AroProp      -0.422
```

---

## Verification

この RP は ESOL ワークフローに対するベースライン calibration run であり、RP00 自体は pass/fail 判定対象ではありません。
ここで較正した実務上の目安は、RP01 以降では拘束力のある RF03 条件として使います。

| 指標 | 実務上の目安 |
|---|---|
| 5-fold CV RMSE | `<= 1.20` |
| 5-fold CV R^2 | `>= 0.75` |

RMSE 閾値は Delaney (2004) の training RMSE `0.996` に、RDKit `MolLogP` 代替による許容差 `0.20` を足した目安です。
`R^2` 閾値は、論文のおよそ `0.84` に対し、RDKit ベース再現ではおよそ `0.76` を見込むという整理に基づきます。
これらは、この RP で使う公開データセットと RDKit ベース記述子の組み合わせに対する目安であり、元論文の閾値そのものではありません。

---

## Discussion

### 元論文との主な差分

| 差分 | 詳細 |
|---|---|
| LogP 実装 | Delaney は別系統の `clogP`、本 RP は RDKit Crippen-Wildman `MolLogP` を使用 |
| 分子数 | 論文は `1,144` 件、本 RP は MoleculeNet 版 `1,128` 件 |
| 評価方法 | 論文は学習データ fit の記述が中心、本 RP は 5-fold CV を重視 |
| Rotatable bond 定義 | RDKit strict 定義は論文中の暗黙定義と一致しない可能性がある |

### 係数が完全一致しない主因

- 最大の要因は `LogP` 実装差
- `NumRotatableBonds` は定義差に敏感で、この run では統計的にも重要ではない
- 公開データセットの curated 版と論文データとの差も fit をずらす

### 実務上の見方

EasyMolKit 上で簡潔な公開 ESOL ベースラインが欲しいなら、この RP は適しています。
より高い ESOL 精度を狙う場合は、この 4 記述子モデルに厳密一致を求めるより、後続 RP を参照する方が適切です。

---

## Files

| ファイル | 内容 |
|---|---|
| `README.md` | 英語正本 README |
| `README.jp.md` | 日本語副本 README |
| `rp00_esol_pilot.m` | 再現スクリプト |
| `lock_template.json` | 実行メタデータ用テンプレート |

スクリプト実行後には、`result/runs/<ts>/` 以下に metrics、predictions、environment snapshot などのローカル成果物が生成されます。

