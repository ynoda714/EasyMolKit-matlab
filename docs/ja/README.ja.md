# EasyMolKit

🇬🇧 [English](../../README.md)

**MATLAB から RDKit を簡単に使える、統合ケモインフォマティクス環境**

> リポジトリをクローン → セットアップコマンド 1 本 → RDKit が MATLAB 上で即利用可能（Python は自動配備、手動設定不要）

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=ynoda714/EasyMolKit-matlab&file=main_rdkit.m)

---

## なぜ EasyMolKit?

ケモインフォマティクス分野では Python + RDKit が事実上の標準ツールチェーンです。
しかし、まず Python のインストールと設定が必要で、さらに複数の障壁があります。

- Python の事前インストールが必要 — ケモインフォマティクスのコードを動かす前に必須
- Python のバージョン管理・仮想環境の設定
- conda / pip を使った RDKit インストール
- 商用ツール（PyMOL 等）との環境競合
- Python のエコシステムに不慣れな ユーザー

**EasyMolKit はこれらの障壁を取り除きます。** MATLAB の `pyenv` ベースの Python 連携を活用し、
RDKit の機能を標準的な MATLAB 関数のように利用できます。Python の手動インストールも知識も原則不要です。

---

## 主な機能

- **ゼロ設定**: Python の事前インストール不要 — `emk.setup.install()` の 1 回呼び出しで Embedded Python + RDKit を自動ダウンロード・配備
- **MATLAB ネイティブ**: 結果は MATLAB の `table` / `struct` / `double` として返却。ワークスペースですぐに利用可能
- **Desktop & Online 対応**: Windows Desktop および MATLAB Online をサポート（macOS / Linux Desktop は未テスト）
- **RDKit ラッパー**: 分子解析・記述子計算・フィンガープリント・類似度検索をすべて MATLAB 関数として提供
- **拡張性**: 将来バージョンで PyMOL Open-Source による 3D 構造可視化を予定

---

## 対象ユーザー

- **化学・薬学・医学の研究者** で MATLAB を主要な研究環境として使用している方
- **ケモインフォマティクスを学ぶ学生**
- **Python 環境構築に時間を取られたくない MATLAB ユーザー**

---

## 動作環境

| 項目 | Desktop | MATLAB Online |
|---|---|---|
| MATLAB | R2025b 以降 | 最新版（MathWorks 管理） |
| Python | 手動設定不要（Embedded Python を自動配備） | 事前インストール済み |
| RDKit | 自動配備 | `emk.setup.installOnline()` でインストール |
| ネット接続 | 初回セットアップ時に必要（Python・RDKit を自動ダウンロード） | `emk.setup.installOnline()` 実行時に必要 |
| OS | Windows | — |

---

## クイックスタート

```matlab
% 1. リポジトリをクローン
%    git clone https://github.com/ynoda714/EasyMolKit-matlab.git
%    cd EasyMolKit-matlab

% 2. MATLAB で開き、初回セットアップを実行
addpath(genpath("src"));
emk.setup.install();          % Desktop
% emk.setup.installOnline();  % MATLAB Online

% 3. 試してみよう
mol = emk.mol.fromSmiles("CCO");          % エタノール
mw  = emk.descriptor.molWeight(mol);      % 分子量
fp  = emk.fingerprint.morgan(mol);        % Morgan フィンガープリント

% 4. 分子間の類似度を比較
mol2  = emk.mol.fromSmiles("CCCO");       % プロパノール
fp2   = emk.fingerprint.morgan(mol2);
score = emk.similarity.tanimoto(fp, fp2); % Tanimoto 係数
```

詳細は [クイックスタートガイド（日本語）](quickstart_ja.md) を参照してください。

---

## 追加ライブラリ（Track 1 & Track 2）

EasyMolKit はアドオンライブラリを 2 つのトラックで管理します。

| トラック | ライブラリ | インストール方法 | ライセンス |
|---|---|---|---|
| **Track 1** | pubchempy, mordred, biopython, torch, torch_geometric, transformers, datasets 等 | `emk.setup.installExtra()` で Embedded Python に追加 | MIT / BSD-3 / Apache-2.0 |
| **Track 2** | Open Babel, MDAnalysis, PyMOL OSS | 別途 CPython 環境が必要。`emk.setup.useExternal()` で接続 | GPLv2 / GPLv2+ / BSD |

---

## チュートリアル＆サンプル（4 層構成）

`examples/` 以下に段階的な学習コンテンツを用意しています。

| Layer | 対象 | 必要 Toolbox | 内容 | リリース |
|---|---|---|---|---|
| **L1 Foundation** | 全ユーザー | なし | API を 1 つずつ体験（6 モジュール、各 5〜15 分） | ✅ v1.0.0 |
| **L2 Application Stories** | Foundation 修了後 | なし | 複数機能を組み合わせた実践ワークフロー（7 モジュール、各 20〜40 分） | 🔜 v1.1.0 |
| **L3 Analytics** | 全ユーザー | 統計・ML 等（種別による） | QSAR・クラスタリング・MS 解析・最適化 等（A01〜A10、各 30〜60 分） | 🔜 v1.2.0 |
| **L4 Research** | 全ユーザー | 種別による | 研究レベルの応用（R01〜R10、各 30〜90 分） | 🔜 v1.3.0 |



### Layer 1: Foundation（ベース MATLAB のみ）

| # | タイトル | 必要 Toolbox | 学べること | Desktop | Online |
|---|---|---|---|:---:|:---:|
| F01 | SMILES で分子を描く | なし | 分子表現、SMILES 構文 | ✔ | ✔ |
| F02 | 分子特性を計算する | なし | MW / LogP / TPSA の意味と計算方法 | ✔ | ✔ |
| F03 | フィンガープリント入門 | なし | ビットベクトル表現、Morgan vs MACCS | ✔ | ✔ |
| F04 | 類似度で分子を比較する | なし | Tanimoto / Dice 類似度の定量化 | ✔ | ✔ |
| F05 | 部分構造検索 | なし | SMARTS パターンマッチング | ✔ | ✔ |
| F06 | ファイルから分子を読み込む | なし | SDF / SMILES ファイル操作 | ✔ | ✔ |

---

## ディレクトリ構成

```
EasyMolKit-matlab/
├─ main_rdkit.m               # RDKit セットアップ＆基本操作（セクション実行）
├─ config/
│   └─ settings.example.json  # 設定テンプレート
├─ examples/
│   ├─ japanese/              # 配布教材 — 日本語版（plain-text Live Code）
│   │   ├─ foundation/        #   L1: API 基礎（ベース MATLAB のみ）
│   │   ├─ stories/           #   L2: 応用ストーリー
│   │   ├─ analytics/         #   L3: 統計・ML 連携
│   │   └─ research/          #   L4: 研究レベル
│   └─ english/               # 配布教材 — 英語版（コメントのみ異なる）
├─ src/
│   ├─ +emk/                  # メインパッケージ
│   │   ├─ +setup/            # Python 環境セットアップ
│   │   ├─ +mol/              # 分子オブジェクト操作
│   │   ├─ +descriptor/       # 記述子計算
│   │   ├─ +fingerprint/      # フィンガープリント生成
│   │   ├─ +similarity/       # 類似度計算
│   │   ├─ +filter/           # 分子フィルタリング（Lipinski 等）
│   │   ├─ +io/               # ファイル入出力
│   │   ├─ +viz/              # 可視化
│   │   └─ +util/             # パッケージ内ユーティリティ
│   ├─ config/                # 設定ロード
│   └─ util/                  # ログヘルパー・共通ユーティリティ
├─ result/                    # 実行成果物（Git 非追跡）
├─ tests/
│   ├─ unit/                  # matlab.unittest クラスベーステスト
│   └─ smoke/                 # スモークテスト
├─ data/                      # キュレーション済みサンプルデータ
└─ docs/                      # ドキュメント
    └─ ja/                    # 日本語ドキュメント
```

---

## API 一覧

| モジュール | 主な関数 | 説明 |
|---|---|---|
| `emk.setup` | `install()`, `verify()`, `initPython()` | Python 環境の自動配備・初期化 |
| `emk.mol` | `fromSmiles()`, `toSmiles()`, `isValid()`, `hasSubstruct()` | 分子オブジェクトの生成・変換 |
| `emk.descriptor` | `molWeight()`, `calculate()`, `batchCalculate()` | 分子記述子の計算 |
| `emk.fingerprint` | `morgan()`, `maccs()`, `toArray()` | フィンガープリント生成 |
| `emk.similarity` | `tanimoto()`, `dice()` | 分子類似度計算 |
| `emk.io` | `readSdf()`, `writeSdf()`, `readSmilesList()` | SDF / SMILES ファイル入出力 |
| `emk.viz` | `draw2d()` | 2D 構造描画 |

完全な API 仕様は [function_reference.md (English)](../function_reference.md) を参照してください。
関数リファレンス（コンパクト版・日本語）は [docs/ja/function_reference_ja.md](function_reference_ja.md) を参照してください。

---

## ライセンス

EasyMolKit: [MIT License](../../LICENSE)

### サードパーティライセンス

| ライブラリ | ライセンス | 用途 |
|---|---|---|
| RDKit | BSD-3-Clause | ケモインフォマティクスコア |
| Python (CPython) | PSF License | ランタイム環境 |
| PyMOL Open-Source | Python/BSD 系 | 3D 可視化（将来リリース） |

詳細は [THIRD_PARTY_NOTICES.md](../../THIRD_PARTY_NOTICES.md) および [docs/compliance.md](../compliance.md) を参照してください。

---

## コントリビューション

バグ報告・機能要望・Pull Request を歓迎します。
詳細は [CONTRIBUTING.md](../../CONTRIBUTING.md) をご覧ください。

---

## 免責事項

- 本プロジェクトは個人の余暇活動として開発・維持されています。
- 結果の正確性については利用者の責任で検証してください
- サポートは best-effort です。応答に時間がかかる場合があります
- 本ソフトウェアは「現状のまま」提供されます（MIT License の免責条項を参照）
- 明示・黙示を問わず、いかなる保証も提供しません
- 利用によって生じた損害について開発者は責任を負いません
- 外部データソース（PubChem, ChEMBL 等）の利用はそれぞれの利用規約に従ってください
- 医療・安全上の判断への直接利用には専門家の確認が必要です

---

## ドキュメント

| ファイル | 内容 |
|---|---|
| [クイックスタートガイド（日本語）](quickstart_ja.md) | セットアップ手順・FAQ |
| [docs/function_reference.md](../function_reference.md) | 関数シグネチャ一覧 |
| [docs/ja/function_reference_ja.md](function_reference_ja.md) | 関数リファレンス（日本語コンパクト版） |
| [docs/ja/test_catalog_ja.md](test_catalog_ja.md) | テストカタログ（日本語） |
| [docs/python_integration.md](../python_integration.md) | Python 連携アーキテクチャ |
| [docs/platform_support.md](../platform_support.md) | Desktop / Online プラットフォーム詳細 |
| [docs/compliance.md](../compliance.md) | ライセンス・コンプライアンス |
