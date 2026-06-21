# チュートリアル一覧 — EasyMolKit

全レイヤーのチュートリアルモジュール一覧です。
層の概要は [README.ja.md](README.ja.md) を参照してください。

---

## Layer 1: Foundation（Base MATLAB のみ）

| # | タイトル | 学習内容 | Desktop | Online |
|---|---|---|:---:|:---:|
| F01 | SMILES で分子を描く | 分子表現、SMILES 構文 | ✔ | ✔ |
| F02 | 分子物性を計算する | MW / LogP / TPSA の意味と計算 | ✔ | ✔ |
| F03 | フィンガープリント入門 | ビットベクトル表現、Morgan vs MACCS | ✔ | ✔ |
| F04 | 類似度で分子を比較する | Tanimoto / Dice 類似度の定量化 | ✔ | ✔ |
| F05 | 部分構造検索 | SMARTS パターンマッチング | ✔ | ✔ |
| F06 | ファイルから分子を読み込む | SDF / SMILES ファイル操作 | ✔ | ✔ |

追加ツールボックス不要。MATLAB Online Basic（無料枠）で完全動作します。

---

## Layer 2: Application Stories（Base MATLAB のみ）

| # | タイトル | ドメイン | Desktop | Online |
|---|---|---|:---:|:---:|
| S01 | カフェインの仲間を探す | 身近な化学 | ✔ | ✔ |
| S02 | 創薬フィルター：リピンスキーの Rule of Five | 薬理学 | ✔ | ✔ |
| S03 | 危険化合物の構造アラート | 安全性 | ✔ | ✔ |
| S04 | バーチャルスクリーニング入門 | 創薬 | ✔ | ✔ |
| S05 | 未知化合物同定チャレンジ | 法科学 | ✔ | ✔ |
| S06 | PubChem で化合物を検索する | データベース | ✔ | ✔ |
| S07 | ChEMBL の活性データを解析する | 創薬 | ✔ | ✔ |

追加ツールボックス不要。MATLAB Online Basic（無料枠）で完全動作します。

---

## Layer 3: Analytics

| # | タイトル | 必要 Toolbox | 内容 | Desktop | Online |
|---|---|---|---|:---:|:---:|
| A01 | PCA による化学空間マッピング | Statistics and Machine Learning Toolbox | 次元削減・化学空間可視化 | ✔ | ✔ |
| A02 | 分子クラスタリング | Statistics and Machine Learning Toolbox | 階層クラスタリング・構造類似性 | ✔ | ✔ |
| A03 | QSAR 回帰 | Statistics and Machine Learning Toolbox | LogP 予測・回帰モデル評価 | ✔ | ✔ |
| A04 | 薬物分類 | Statistics and Machine Learning Toolbox | SVM / ランダムフォレスト・ROC 曲線 | ✔ | ✔ |
| A05 | ニューラルネットワーク物性予測 | Deep Learning Toolbox | フィードフォワード NN・物性回帰 | ✔ | ✔ |
| A06 | 用量反応カーブフィット | Curve Fitting Toolbox | Hill 式・EC50 推定 | ✔ | ✔ |
| A07 | Scaffold 分析と R グループ分解 | Statistics and Machine Learning Toolbox | 創薬化学分析 | ✔ | ✔ |
| A08 | 質量分析×ケモインフォマティクス | Signal Processing Toolbox + Statistics and Machine Learning Toolbox | 同位体パターンマッチング・MS アノテーション | ✔ | ✔ |
| A09 | PFAS・環境スクリーニング | Optimization Toolbox + Statistics and Machine Learning Toolbox（任意） | SMARTS スクリーニング・Pareto 最適化 | ✔ | ✔ |
| A10 | リード最適化 | Optimization Toolbox（任意） | 多目的最適化・Derringer-Suich 法 | ✔ | ✔ |

MATLAB Online Basic（無料枠）で完全動作します。

---

## Layer 4: Research

| # | タイトル | 必要 Toolbox | Desktop | Online |
|---|---|---|:---:|:---:|
| R01 | 大規模類似度スクリーニング（GPU） | Parallel Computing Toolbox (GPU) | ✔ | △（CPU Only） |
| R02 | PK/PD シミュレーション | SimBiology | ✔ | ✔ |
| R03 | 法科学ケモメトリクス | Statistics and Machine Learning Toolbox + Parallel Computing Toolbox | ✔ | ✔ |
| R04 | タンパク質-リガンド解析 † | Bioinformatics Toolbox | ✔ | ✔ |
| R05 | 分子言語モデル：SMILES 生成 | Deep Learning Toolbox | ✔ | ✔ |
| R06 | REINFORCE 分子設計 | Deep Learning Toolbox + Reinforcement Learning Toolbox | ✔ | ✔ |
| R07 | メタボロミクス † | Bioinformatics Toolbox + SimBiology | ✔ | ✔ |
| R08 | タンパク質-リガンド ドッキングシミュレーション ‡ | なし（Track 1: meeko + vina + pdbfixer） | ✕ | ✔ |
| R09 | GNN 分子性質予測 § | Deep Learning Toolbox | ✔ | ✔ |
| R10 | ChemBERTa 転移学習 § | Deep Learning Toolbox | ✔ | ✔ |

> **†** 初回実行前に `emk.setup.installExtra("biopython")` が必要です（Track 1 追加ライブラリ。MATLAB ライセンスとは別要件）。
>
> **‡ MATLAB Online 限定**（Windows Desktop 非対応: vina は Windows PyPI ホイールなし・pdbfixer の openmm は Smart App Control でブロック）。
> セットアップ: `main_emk.m` で `cfg.optionalLibraries.meeko/vina/pdbfixer = true` に設定 → `installOnline(Config=cfg)` で一括導入。
>
> **§** PyTorch + HuggingFace スタックが必要です。以下の順でインストールしてください:
> `emk.setup.installExtra("torch")` → `emk.setup.installExtra("torch_geometric")` → `emk.setup.installExtra("transformers")` → `emk.setup.installExtra("datasets")`。
> R10 は R09 の torch 環境が前提です。

---

## 再現可能研究

| ID | 論文 | 手法 | 結果 |
|---|---|---|---|
| RP00 | Delaney (2004) ESOL — 水溶性 | 物理化学記述子による線形回帰 | CV RMSE=1.017, R²=0.762 |
| RP01 | Delaney (2004) ESOL — 拡張版 | 線形回帰 + TPSA / QED / SA Score | CV RMSE=0.584, R²=0.906 |
| RP02 | Wu et al. (2018) MoleculeNet BBBP ベースライン | Morgan FP (ECFP4) + Random Forest | ROC-AUC CV=0.883 |
| RP03 | Yang et al. (2019) GNN on BBBP | Graph Convolutional Network | ROC-AUC CV=0.915 |
| RP04 | Chithrananda et al. (2020) ChemBERTa | 凍結 CLS 埋め込み + ロジスティック回帰 | ROC-AUC CV=0.927 |
| RP05 | SHAP 説明可能 AI on BBBP | shap.LinearExplainer + LR モデル | ROC-AUC CV=0.909, Spearman ρ=0.902 |

各エントリは `repro/rp*/` 以下に MATLAB スクリプト・環境ロック（`lock_template.json`）・結果詳細（`README.en.md`）を含みます。
