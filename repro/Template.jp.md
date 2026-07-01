# RP<XX>: <タイトル> — 標準再現テンプレート（RF01）

> **目的**: RF01 標準再現テンプレート。各再現研究はこの構造に従うこと。
> RP00 の経験を踏まえ M-REPRO-FOUND (2026-06-19) で確立。RF03/RF04 は M-REPRO-AUDIT C1/C2 (2026-06-22) で正式定義。
> `README.md` を英語正本、`README.jp.md` を日本語副本とし、相互リンクを置くこと。
> 両 README は原則として見出し順・表の粒度・記載対象を対称に保つこと。

> **候補論文の探し方**: AnyResearch（OpenAlex 連携ツール）を使うと、
> ケモインフォクエリで再現可能な論文を系統的に探索できる。
> 手順は AnyResearch の `docs/workflows/repro_discovery.md` を参照。

---

## 再現フレームワーク（RF00–RF04）概要

| # | 名称 | 定義 | 主な成果物 |
|---|---|---|---|
| RF00 | 論文選定基準 | 再現対象論文のスコアリングポリシー | `docs/algorithms/repro_selection.md` |
| RF01 | 再現テンプレート | 本ドキュメント（README 構造・記述子定義表） | `repro/Template.md` / `repro/Template.jp.md` |
| RF02 | バージョンロック | 実行環境を `lock_snapshot.json` に記録 | `emk.setup.snapshot()` / `emk.setup.lockfile()` |
| RF03 | 検証基準 | 評価指標の 3 カテゴリ検証（下記参照） | `emk.repro.verify()` |
| RF04 | 外部再現プロトコル | 外部貢献者が再現を提出する手順（下記参照） | `CONTRIBUTING.md` |

### RF03 検証カテゴリ（M-REPRO-AUDIT C1）

| カテゴリ | 名称 | 内容 | RP PASS への影響 |
|---|---|---|---|
| **Cat A** | 絶対閾値 | `emk.repro.verify()` による upper/lower bounds 照合（RMSE ≤ X、R² ≥ Y、AUC ≥ Z など）| **必須**。全指標が基準内 |
| **Cat B** | 相対比較 | paired t-test + p 値による 2 条件間の有意差検定（例: Model A vs B）| **参考情報**（PASS/FAIL に非直結）|
| **Cat C** | 実装健全性 | 数学的正確さのサニティチェック（例: SHAP Spearman ρ）| 該当実装がある場合は**必須** |

> **RP PASS 判定**: Cat A（必須）+ Cat C（該当時）が全て PASS であること。Cat B は比較研究の文脈情報として記録する。
>
> **`_rev` スクリプトでの Cat B 扱い（設計確定）**: 修正版スクリプト（`_rev`）が元スクリプトとの差分を `logInfo` で報告する「前後比較」は Cat B（参考情報）扱いとする。PASS/FAIL ゲートには昇格しない。修正の目的は精度向上でなく公正性回復（バイアス除去・リーク修正）であるため、修正後 AUC が元より低くなっても PASS は維持される。

### RF04 外部再現プロトコル（M-REPRO-AUDIT C2）

**RF04 準拠の定義**: 以下 3 要件を全て満たすとき、当該 RP は「外部再現可能」と認定される。

| 要件 | 内容 | 確認方法 |
|---|---|---|
| RF01 ✓ | 本テンプレート形式の `README.md` / `README.jp.md` が存在する | Overview / Environment / Data / Script / Result / Verification / Discussion の全セクションが両言語で記入済み |
| RF02 ✓ | `lock_snapshot.json` が `result/runs/<ts>/` に存在する | `emk.setup.lockfile()` 呼び出し済み |
| RF03 ✓ | Cat A（+ 該当時 Cat C）が PASS | `metrics.json` の `rf03_pass: true` |

> 外部貢献者が新たな再現を提出する手順は `CONTRIBUTING.md` の「Submitting a Reproduction (RF04)」を参照。

---

## スクリプト作成時チェックリスト（開発者用・README に含めない）

> **根拠**: RP00〜RP07 のコードレビュー修正を横断分析した結果、以下 10 カテゴリの問題が複数 RP で独立して再発した（2026-06-25 抽出）。
> 新規 RP のスクリプトを scaffold した直後にこのリストを確認し、**初回作成時から満たす**こと。

### A — 再現性・RNG 管理

- [ ] `rng(seed)` を `cvpartition` 直前に配置（Section 0 の宣言だけでは不足）
- [ ] Python fold ループ先頭: `torch.manual_seed(seed+fold_idx)` + `np.random.seed(seed+fold_idx)`
- [ ] DataLoader: `generator=torch.Generator().manual_seed(seed+fold_idx)` を指定
- [ ] モデルロード（`from_pretrained` 等）直前: `random.seed(seed)` を配置

### B — NaN / データ品質ガード

- [ ] 全 descriptor 列をまとめて除外: `any(ismissing(descTbl), 2)` で一括マスク
- [ ] Python 側: `MolFromSmiles` 前に `notna()` チェックを置く
- [ ] per-molecule 変換失敗: NaN init + `try/except` + `logWarn` フォールバック

### C — ファイルハンドル / Figure 管理

- [ ] `fopen/fprintf/fclose` は `try/finally` で囲む（または `writelines` に置換）
- [ ] `fopen` 直後に `if fid == -1; error(...); end` ガード
- [ ] `saveas` 直後に `close(fig)`（繰り返しセクション実行でウィンドウ蓄積防止）
- [ ] Figure close も `try/finally` で保証

### D — メトリクス定義の精確化

- [ ] CV R² の baseline = `mean(y_all)`（global）。per-fold 局所平均は NG
- [ ] Training RMSE = `sqrt(RSS/n)`（バイアスあり）。`sqrt(mdl.MSE)` は RSS/(n−p−1) なので NG
- [ ] RF03 pass 条件は複数モデルの論理積: `resA.pass && resB.pass`
- [ ] OOF RMSE は全 fold の予測を結合してから計算（fold 平均でない）

### E — パス・runDir 統一

- [ ] `makeRunDir` 戻り値を `runDir` として使用。`absRunDir` は作らない
- [ ] `runDir` を参照する Section の先頭に `if ~exist('runDir','var'); error(...); end` ガード
- [ ] `makeRunDir` の呼び出しは結果保存 Section の先頭（可視化 Section ではない）

### F — batchMorganFP\_ / IPC バリデーション

- [ ] IPC は list オブジェクトをそのまま渡す（concatenated-string IPC 禁止）
- [ ] Python 側: reshape 前に `assert len(bits) == n_molecules * n_bits`
- [ ] MATLAB 側: reshape 前に per-molecule 長チェックを追加

### G — データリーク防止（説明可能 AI / CV 系）

- [ ] SHAP LinearExplainer の background = `X_tr` のみ（test 混入禁止）
- [ ] `cross_val_score` には unfitted estimator を渡す
- [ ] CV スコープと eval スコープを明示的に分離（80/20 split → CV は train 側のみ）

### H — ハイパーパラメータの明示

- [ ] `fitrlinear`: `Lambda=1.0` を明示（MATLAB バージョン間のデフォルト変動を排除）
- [ ] `TreeBagger`: `'OOBPrediction','off'` を全モデルに明示
- [ ] モデルの全ハイパーパラメータを明示引数で渡す（デフォルト依存禁止）

### I — Python コアの堅牢性

- [ ] `try/except` で捕捉 → `{"success": false, "error": "..."}` JSON を常に return
- [ ] PyTorch: `.detach().numpy()` に統一（PyTorch >= 2.0 対応）
- [ ] `max_len` 可変設計なら `batch_size = max(1, 32*128//max_len)` で自動スケール

### J — README 同期

- [ ] 英語 README (`README.md`) と日本語 README (`README.jp.md`) を同時に作成
- [ ] 見出し順・表の列構成・記載対象を日英で対称に揃える
- [ ] 実行後の Result 値を反映し「改訂日」を明記
- [ ] Section テーブルをコード構成に合わせる（Section 番号・説明の一致）
- [ ] `lock_template.json` に全 RF03/RF04 メトリクスフィールドを宣言（single source of truth）

---

## Overview

| 項目 | 内容 |
|---|---|
| 論文 | 著者 (年). タイトル. *ジャーナル* 巻(号):ページ. |
| DOI | [10.xxxx/yyyy](https://doi.org/10.xxxx/yyyy) |
| タスク | 回帰 / 分類 / 記述子計算 / その他 |
| モデル | モデル名・手法 |
| データ | データセット名（件数）|
| 公開指標 | 論文記載の評価指標と値 |

---

## Environment（RF02 バージョンロック）

実行後 `result/runs/<timestamp>/lock_snapshot.json` に実際のバージョンが記録される。

| 項目 | 要件 |
|---|---|
| MATLAB | R2025a 以降 |
| Python | 3.10（EasyMolKit Embedded Python）|
| RDKit | バージョンを明記（例: 2022.03 以降）|
| Toolbox | 必要 Toolbox を列挙 |

### 使用記述子の定義（RF01 必須項目）

再現研究で用いた各記述子について、以下の三点セットを明記する。

| 記述子 | 計算ツール | バージョン | 定義・算出定義 |
|---|---|---|---|
| LogP | RDKit | 記録要 | Crippen-Wildman `MolLogP`。clogP / ALogP / XLogP とは異なる |
| MolWt | RDKit | 記録要 | `Descriptors.MolWt`（全重原子 + 暗黙的H の IUPAC 平均原子量）|
| NumRotatableBonds | RDKit | 記録要 | `CalcNumRotatableBonds` strict SMARTS 定義（末端結合を含まない）|
| （追加記述子） | （ツール）| — | （定義）|

> **ノート**: RDKit のバージョンによっては一部 API が存在しない場合がある。
> 代替手段を用いた場合は本セクションに注記すること（例: `CalcNumAromaticAtoms` 不在 → `pyrun` バッチで代替）。

---

## Data

- **ソース**: データセット名・提供元
- **URL**: ダウンロードURL（またはライセンスにより非公開の場合はその旨）
- **ライセンス**: ライセンス名・条件
- **キャッシュ先**: `data/benchmark/<name>.csv`（初回実行時に自動ダウンロード）
- **件数**: 〇〇 件（論文の〇〇 件と一致/不一致の場合は差異を説明）
- **データハッシュ**: `result/runs/<ts>/lock_snapshot.json` の `dataset_sha256` フィールドを参照

---

## Script

```
repro/rp<XX>_<name>/rp<XX>_<name>.m
```

**実行方法**: MATLAB でプロジェクトルートを CWD にして Ctrl+Enter でセクション実行。

| セクション | 内容 |
|---|---|
| Section 0 | セットアップ・環境情報取得（`emk.setup.snapshot()` 呼び出し）|
| Section 1 | データセット読み込み |
| Section 2 | 前処理・記述子計算 |
| Section 3 | モデル学習 |
| Section 4 | 評価（CV / test set）|
| Section 5 | RF03 検証（`emk.repro.verify()` 呼び出し）|
| Section 6 | 結果保存（`makeRunDir()` → `emk.setup.lockfile()`）|

---

## Result（初回実行 YYYY-MM-DD）

| 指標 | 値 | 合否 |
|---|---|---|
| （指標名） | （値） | ✅ PASS / ❌ FAIL |

**環境 (lock_snapshot.json より):**

| 項目 | 値 |
|---|---|
| MATLAB | （実行時のバージョン）|
| Python | （実行時のバージョン）|
| RDKit | （実行時のバージョン）|
| Commit | （実行時のハッシュ）|

---

## Verification（RF03 数値検証基準）

> RF03 カテゴリ: **Cat A（絶対閾値）** — Cat B / Cat C がある場合は対応サブセクションを追加すること。
> RF04: ✅ 準拠（RF01 / RF02 / RF03 Cat A 充足時に記入）

### Cat A — 絶対閾値（必須）

| 指標 | 基準 | 根拠 |
|---|---|---|
| （指標名） | （上限 or 下限）| （根拠：論文値・許容差の理由）|

**許容差の根拠**:
1. （記述子実装の差異など）
2. （データセット版数の差異など）

> **注意**: 初回実行の結果を踏まえて閾値を再調整した場合は、ここにその経緯を記録すること。

<!-- Cat B（相対比較）: 2 モデル比較が必要な場合に追加 -->
<!-- ### Cat B — 相対比較（参考情報） -->
<!-- | 指標 | Model A | Model B | delta | t(df) | p（片側） | -->
<!-- Cat C（実装健全性）: サニティチェックが必要な場合に追加 -->
<!-- ### Cat C — 実装健全性（該当時必須） -->

---

## Discussion

### 論文との差異

| 差異 | 詳細 |
|---|---|
| （差異の種類）| （詳細・影響の説明）|

### 学習事項（後続 RP / M-REPRO-FOUND への引き継ぎ）

- [ ] （学習点・課題）

---

## 関連ファイル

| ファイル | 内容 |
|---|---|
| `rp<XX>_<name>.m` | 再現スクリプト本体 |
| `lock_template.json` | RF02 バージョンロックのスキーマテンプレート |
| `result/runs/<ts>/lock_snapshot.json` | 実行時の実際のバージョン情報 |
| `result/runs/<ts>/metrics.json` | 評価指標 |
