# EasyMolKit

🇬🇧 [English](../en/README.md)

**MATLAB から RDKit を簡単に使えるケモインフォマティクス統合環境**

> リポジトリをクローンしてセットアップコマンドを 1 本実行するだけで、RDKit が MATLAB 上で使えるようになります

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=ynoda714/EasyMolKit-matlab)

> *[MATLAB Home License](https://www.mathworks.com/products/matlab-home.html) を使った個人の余暇プロジェクトです。学習・趣味の目的で公開しており、MathWorks や商業組織とは無関係です。*

## なぜ EasyMolKit？

ケモインフォマティクス分野では Python + RDKit が事実上の標準ツールチェーンです。
しかし、その環境構築には多くのハードルがあります。

- Python バージョンや仮想環境の管理
- conda / pip による複雑な RDKit インストール
- 商用ツール（PyMOL など）との環境競合
- MATLAB ユーザーには馴染みのない Python エコシステム

**EasyMolKit はこれらのハードルをすべて解消します。** MATLAB の `pyenv` ベースの Python 連携を活用し、
RDKit の機能を通常の MATLAB 関数として提供します — Python の知識は一切不要です。

## 特徴

- **ゼロ設定**: `emk.setup.install()` を 1 回呼ぶだけで Python + RDKit を自動配備
- **MATLAB ネイティブ**: 結果は `table` / `struct` / `double` で返却 — ワークスペースで即利用可能
- **Desktop & Online 対応**: Windows Desktop および MATLAB Online をサポート（macOS / Linux Desktop は未検証）
- **豊富な API**: 15 モジュール 76 関数 — 記述子・フィンガープリント・スキャフォルド・フィルター・クラスタリング・3D コンフォーマーなど
- **再現可能研究**: `repro/` に環境ロック付きで 6 本の論文再現実験（RP00〜RP05）を収録

## 対象ユーザー

- **化学・薬学・医学分野の研究者**で MATLAB をメインの研究環境として使用している方
- **ケモインフォマティクスを学ぶ学生**（MATLAB Online Basic 無料枠で Layer 1〜3 が動作）
- **Python 環境構築に時間をかけたくない MATLAB ユーザー**

## 動作要件

| 項目 | Desktop | MATLAB Online |
|---|---|---|
| MATLAB | R2025b 以降 | R2025b 以降 |
| Python | 自動配備（Embedded Python） | システム Python が使用される |
| RDKit | 自動配備 | `emk.setup.installOnline()` でインストール |
| OS | Windows | — |

## クイックスタート

```matlab
% 1. リポジトリをクローン
%    git clone https://github.com/ynoda714/EasyMolKit-matlab.git
%    cd EasyMolKit-matlab

% 2. main_emk.m を MATLAB で開き、各セクションを Ctrl+Enter で実行:
%
%   Section 0a  — パスセットアップ & 設定      (必要に応じて cfg.useCase.* を編集)
%   Section 0b  — Python + RDKit セットアップ  (初回のみ; 約 2〜5 分)
%   Section 1   — 分子操作の基礎
%   Section 2   — 記述子計算
%   Section 3   — フィンガープリント & 類似度
```

> ⚠️ **Ctrl+Enter（セクション実行）** を使用してください — **F5（ファイル全体実行）** は初回セットアップ時に失敗します。

詳細は [docs/quickstart.md](../quickstart.md) を参照してください。

### 企業 PC・制限付きネットワーク環境

**企業 PC** をご利用の場合、IT セキュリティポリシーによってローカル Python の配備がブロックされる可能性があります:

| 問題 | 症状 | 対処法 |
|---|---|---|
| プロキシサーバー | `pip install` タイムアウト / SSL エラー | `main_emk.m` Section 0a で `cfg.python.proxy = "http://proxy.example.com:8080"` を設定 |
| Windows Defender / Smart App Control | Embedded Python 展開時に隔離・ブロック | `python_env/` をセキュリティソフトの除外対象に追加、または MATLAB Online を使用 |
| IT ポリシー（実行ファイルダウンロード禁止） | セットアップのダウンロード手順で失敗 | **MATLAB Online** を使用（ローカル Python 配備が不要） |
| ウイルス対策ソフトによる隔離 | 展開後に Python バイナリが消える | `python_env/` を除外対象に追加、または MATLAB Online を使用 |

> 💡 **企業環境での推奨**: [MATLAB Online](https://matlab.mathworks.com/open/github/v1?repo=ynoda714/EasyMolKit-matlab) をご利用ください — ローカル Python インストールが不要で、L1〜L3 チュートリアルはすべて無料の Basic プランで動作します。

## 追加ライブラリ（Track 1 & Track 2）

EasyMolKit はアドオンライブラリを 2 つのトラックで管理します。

| トラック | ライブラリ | インストール方法 | ライセンス |
|---|---|---|---|
| **Track 1** | pubchempy, mordred, biopython, torch, torch_geometric, transformers, datasets 等 | `emk.setup.installExtra()` — Embedded Python に直接追加 | MIT / BSD-3 / Apache-2.0 |
| **Track 2** | Open Babel, MDAnalysis, PyMOL OSS | 別途 CPython 環境が必要; `emk.setup.useExternal()` で接続 | GPLv2 / GPLv2+ / BSD |

### Track 1: Embedded Python への追加インストール

```matlab
% インストール前にレシピとライセンスを確認
emk.setup.recipe("pubchempy")          % インストール手順とライセンスを表示
emk.setup.installExtra("pubchempy")    % Embedded Python にインストール
emk.setup.installExtra("mordred")      % 1800 種以上の記述子ライブラリ
emk.setup.installExtra("biopython")    % PDB / 配列解析

% PyTorch + HuggingFace スタック（R09 / R10 で必要）
emk.setup.installExtra("torch")           % CPU のみ、約 800 MB（最初にインストール）
emk.setup.installExtra("torch_geometric") % GNN ライブラリ（torch が必要）
emk.setup.installExtra("transformers")    % HuggingFace Transformers
emk.setup.installExtra("datasets")        % HuggingFace Datasets

% インストール状況を確認
T = emk.setup.validate()
```

> 全ライブラリの一括インストールは `main_setup_extra.m` を参照。

### Track 2: 外部 CPython 環境が必要なライブラリ

GPL ライセンスや技術的制約により、Open Babel・MDAnalysis・PyMOL は Embedded Python に
バンドルできません。**別途 CPython 3.10+ 環境**にインストールし、`emk.setup.useExternal()` で接続します。
セットアップ手順の詳細は [docs/quickstart.md — Track 2](../quickstart.md) を参照してください。

## チュートリアル & サンプル

EasyMolKit は `examples/` 配下に段階的な学習コンテンツを提供します。

| 層 | 対象者 | コンテンツ | リリース |
|---|---|---|---|
| **L1 Foundation** | 全ユーザー | API を 1 概念ずつ学ぶ（6 モジュール、各 5〜15 分） | ✅ v1.0.0 |
| **L2 Application Stories** | Foundation 修了後 | 複数機能を組み合わせた実践ワークフロー（7 モジュール、各 20〜40 分） | ✅ v1.1.0 |
| **L3 Analytics** | 全ユーザー | QSAR・クラスタリング・MS 解析・最適化など（A01〜A10、各 30〜60 分） | ✅ v1.2.0 |
| **L4 Research** | 全ユーザー | 研究レベルの応用（R01〜R10、各 30〜90 分） | ✅ v1.3.0 |

*L1〜L3 は MATLAB Online Basic（無料枠）で完全動作します。*

モジュール別のツールボックス要件・プラットフォーム対応は [docs/ja/tutorials.ja.md](tutorials.ja.md) を参照してください。

## 再現可能研究

`repro/` に公開論文の MATLAB 再現実験（6 本）を収録しています。
各エントリに環境ロックスナップショット（RF02）と成功基準（RF03）を定義しています。

手法・データセット・結果の詳細は [docs/ja/repro.ja.md](repro.ja.md) を参照してください。

## MATLAB がカバーできる範囲

`repro/` の実験から、MATLAB のケモインフォマティクスにおける適用範囲を実証的に整理しました。
優劣の判定ではなく条件の記述です——「MATLAB は X、Python は Y」。

| ゾーン | タスク種別 | 主な条件 | 結果 |
|---|---|---|---|
| **A — MATLAB 主体** | 統計解析・フィルター・可視化 | デフォルト設定 | Python 代替と同等以上 |
| **B — 条件付き同等** | 古典的 ML（LR・Ridge・RF）・線形 SHAP | ソルバーと正則化を明示設定 | gap < 1σ（practical tie）；SHAP Spearman ρ = 0.915–0.927 |
| **C — 分業** | 全 ML / DL / LLM パイプライン | Python（RDKit）で特徴量化、MATLAB でモデル学習 | 完全動作：GCN Δ = −0.017 < 1σ、ChemBERTa AUC = 0.914 |
| **D — Python 専有** | 非線形 SHAP（TreeSHAP / KernelSHAP）| `shap` ライブラリが必要 | `TreeExplainer` / `KernelExplainer` は MATLAB に存在しない |

**Zone B — practical parity に必要な 3 条件:**
1. ソルバーを `lbfgs` に明示指定——デフォルトの SGD は高次元スパース特徴量で性能劣化
2. 正則化スケールを意識——`Lambda = 1/n`（MATLAB）≡ `C ≈ n`（sklearn）；規約が逆だが性能は同等
3. 前処理（SMILES → フィンガープリント）は RDKit 必須；ML 本体のみ MATLAB が担う分業が前提

**Zone C — 実証済みパイプライン（RP03・RP04・全 RP）:**
- **記述子パイプライン**（SMILES → 特徴量 → ML/統計）: 全 RP で完全動作——機能ギャップなし
- **GCN / ディープラーニング**（RP03）: MATLAB DLT 3 層 GCN が AUC = 0.887 ± 0.015 を達成（Python 0.904 ± 0.020、Δ = −0.017 < 1σ）；グラフ特徴量化は Python 必須
- **LLM 埋め込み**（RP04）: Python でトークン化 → MATLAB が ONNX 推論 + ロジスティック回帰、AUC = 0.914 ± 0.009（RF03 PASS）；ONNX 忠実性確認済み（F1-a = F1-b の結果が完全一致）

## API 概要

| モジュール | 主な関数 | 説明 |
|---|---|---|
| `emk.setup` | `install()`, `verify()`, `snapshot()`, `verifyLock()` | Python 環境配備・初期化・RF02 バージョンロック |
| `emk.mol` | `fromSmiles()`, `toSmiles()`, `isValid()`, `hasSubstruct()` | 分子オブジェクトの作成 & 変換 |
| `emk.descriptor` | `molWeight()`, `calculate()`, `qed()`, `saScore()`, `bcut()` | 分子記述子の計算 |
| `emk.fingerprint` | `morgan()`, `maccs()`, `toArray()` | フィンガープリント生成 |
| `emk.similarity` | `tanimoto()`, `dice()`, `rankBy()`, `matrix()` | 分子類似度計算 |
| `emk.scaffold` | `genericMurcko()`, `brics()`, `rgroup()` | Scaffold 分析・フラグメント分解 |
| `emk.dataset` | `esol()`, `freesolv()`, `bbbp()`, `tox21()` | ベンチマークデータセット（ローカルキャッシュ付き） |
| `emk.filter` | `lipinski()`, `veber()`, `pains()`, `reos()` | 創薬化学フィルター |
| `emk.cluster` | `butina()` | Butina 球面排除クラスタリング |
| `emk.diversity` | `pick()` | MaxMin 多様性サブセット選択 |
| `emk.conformer` | `embed()`, `optimize()` | 3D コンフォーマー生成・力場最適化 |
| `emk.shape` | `compare()` | 3D 形状類似度比較 |
| `emk.repro` | `verify()` | RF03 再現成功条件の検証 |
| `emk.io` | `readSdf()`, `writeSdf()`, `readSmilesList()` | SDF / SMILES ファイル入出力 |
| `emk.viz` | `draw2d()` | 2D 構造描画（※） |

> **※ 描画に関する注意**: `emk.viz.draw2d()` は RDKit（Python）で PNG を生成し MATLAB に転送します。
> 描画には **1 分子あたり 0.5〜2 秒**かかります。MATLAB Online では
> プロセス間通信のオーバーヘッドがさらに加わります（構造上の制約であり改善不可）。

全 API の詳細は [docs/function_reference.md](../function_reference.md) を参照してください。

## ディレクトリ構成

```
EasyMolKit/
├─ main_emk.m               # RDKit セットアップ & 基本操作（セクション実行）
├─ config/
│   └─ settings.example.json  # 設定テンプレート
├─ examples/
│   ├─ japanese/              # 配布教材（日本語版、plain-text Live Code）
│   └─ english/               # 配布教材（英語版、コメントのみ異なる）
├─ repro/                     # 再現可能研究（RP00〜RP05）
├─ src/
│   └─ +emk/                  # メインパッケージ（15 モジュール、76 関数）
├─ tests/
│   ├─ unit/                  # matlab.unittest クラスベーステスト
│   └─ smoke/                 # スモークテスト
├─ data/                      # キュレーション済みサンプルデータ
└─ docs/                      # ドキュメント
```

## ライセンス

EasyMolKit: [MIT License](../../LICENSE)

### サードパーティライセンス

| ライブラリ | ライセンス | 用途 |
|---|---|---|
| RDKit | BSD-3-Clause | ケモインフォマティクスコア |
| Python (CPython) | PSF License | ランタイム環境 |

詳細は [THIRD_PARTY_NOTICES.md](../../THIRD_PARTY_NOTICES.md) および [docs/compliance.md](../compliance.md) を参照してください。

## コントリビューション

バグ報告・機能リクエスト・プルリクエストを歓迎します。
[CONTRIBUTING.md](../../CONTRIBUTING.md) のガイドラインを参照してください。

## ドキュメント

| ファイル | 説明 |
|---|---|
| [docs/ja/quickstart.ja.md](quickstart.ja.md) | セットアップ手順・Track 2 設定・FAQ（日本語版） |
| [docs/ja/tutorials.ja.md](tutorials.ja.md) | チュートリアル全覧（F01〜R10、RP00〜RP05） |
| [docs/ja/repro.ja.md](repro.ja.md) | 再現可能研究一覧（RP00〜RP05 手法・結果付き） |
| [docs/ja/function_catalog.ja.md](function_catalog.ja.md) | コンパクト関数カタログ（76 関数） |
| [docs/ja/function_reference.ja.md](function_reference.ja.md) | 全関数シグネチャリファレンス |
| [docs/ja/test_catalog.ja.md](test_catalog.ja.md) | テストクラスカタログ |
| [docs/python_integration.md](../python_integration.md) | Python 連携アーキテクチャ |
| [docs/platform_support.md](../platform_support.md) | Desktop / Online プラットフォーム対応 |
| [docs/compliance.md](../compliance.md) | ライセンス & コンプライアンス |

### 英語版ドキュメント

| ファイル | 説明 |
|---|---|
| [docs/en/README.md](../en/README.md) | Repository overview (English) |
| [docs/en/tutorials.md](../en/tutorials.md) | Full tutorial listing (English) |
| [docs/quickstart.md](../quickstart.md) | Setup guide (English) |
| [docs/function_catalog.md](../function_catalog.md) | Compact function catalog (English) |
| [docs/function_reference.md](../function_reference.md) | Full function reference (English) |
