# 再現可能研究：RP00–RP05

> 公開論文の MATLAB 再現実験の全一覧です。
> 各エントリに環境ロックスナップショット（RF02）と成功基準（RF03）を定義しています。
> 手法・考察・再現手順の詳細は各 `repro/<id>/README.md` を参照してください。

## 再現実験一覧

| ID | 論文（年） | タスク | 手法 | データセット | 結果 | Zone |
|---|---|---|---|---|---|---|
| [RP00](../../repro/rp00_esol/) | Delaney (2004) | 回帰（水溶性） | 線形回帰（物理化学記述子） | ESOL 1,128 件 | CV RMSE=1.017, R²=0.762 | A |
| [RP01](../../repro/rp01_esol/) | Delaney (2004) 拡張 | 回帰（水溶性） | 線形回帰 + TPSA/QED/SA Score | ESOL 1,128 件 | CV RMSE=0.980, R²=0.780 | A |
| [RP02](../../repro/rp02_bbbp/) | Wu et al. (2018) MoleculeNet | 分類（BBBP） | ECFP4 + Logistic Regression（Scaffold 5-fold） | BBBP 2,039 件 | ROC-AUC CV=0.883 | B/C |
| [RP03](../../repro/rp03_gnn/) | Yang et al. (2019) | 分類（BBBP） | Graph Convolutional Network（3 層） | BBBP 2,039 件 | ROC-AUC CV=0.915 | C |
| [RP04](../../repro/rp04_chemberta/) | Chithrananda et al. (2020) | 分類（BBBP） | ChemBERTa CLS 埋め込み + 線形回帰 | BBBP 2,039 件 | ROC-AUC CV=0.927 | C |
| [RP05](../../repro/rp05_shap/) | SHAP 解析 on BBBP | 説明可能 AI | shap.LinearExplainer + LR | BBBP 2,039 件 | ROC-AUC CV=0.909, Spearman ρ=0.902 | B/D |

## Zone 凡例

Zone 分類の詳細は [README.md](../../README.md#matlab-がカバーできる範囲) の「MATLAB がカバーできる範囲」を参照してください。

| Zone | 意味 |
|---|---|
| A | MATLAB 主体（Python 側 ML 不要） |
| B | 条件付き同等（ソルバー・正則化を明示設定） |
| C | 分業（Python で特徴量化、MATLAB でモデル学習） |
| D | Python 専有（MATLAB 実装が困難または不可） |

> 一部の RP は複数 Zone にまたがります（例: RP02 は Zone B・C の両方を評価）。
