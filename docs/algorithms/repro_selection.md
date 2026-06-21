# RF00: Reproduction Selection Policy

> **権威**: `docs/algorithms/repro_selection.md`（本ファイル）  
> **概要**: EasyMolKit で再現する論文を選定するためのスコアリング基準（RF00）。  
> DS03 はこの定義を OpenAlex 上で自動計算するのみ（独自定義禁止 / single source of truth）。

---

## スコアリング式

$$\text{Reproduction Score} = 0.30 \times S_{\text{Citation}} + 0.20 \times S_{\text{Data}} + 0.20 \times S_{\text{Code}} + 0.20 \times S_{\text{MATLAB}} + 0.10 \times S_{\text{Edu}}$$

各スコア $S_i \in [0, 5]$（整数）。合計は最大 5.0 点。

---

## 評価項目定義

### Citation Impact（$S_{\text{Citation}}$）

| 点数 | 基準 |
|---|---|
| 5 | 引用数 ≥ 1000 または分野での landmark 論文 |
| 4 | 引用数 500〜999 |
| 3 | 引用数 100〜499 |
| 2 | 引用数 20〜99 |
| 1 | 引用数 5〜19 |
| 0 | 引用数 < 5 または査読なし |

データソース: OpenAlex `cited_by_count`（DS03 で自動取得）。

### Dataset Availability（$S_{\text{Data}}$）

| 点数 | 基準 |
|---|---|
| 5 | `emk.dataset.*` または `data/benchmark/` に既収録 |
| 4 | パブリックドメイン / CC0 で直接ダウンロード可能 |
| 3 | 登録不要の無料ダウンロード（例: MoleculeNet, ChEMBL）|
| 2 | 要登録・無料（例: PubChem, ZINC）|
| 1 | 要申請・限定公開（例: 要論文著者への問い合わせ）|
| 0 | 非公開・商用データのみ |

### Code Availability（$S_{\text{Code}}$）

| 点数 | 基準 |
|---|---|
| 5 | GitHub 等に実装コード公開（テスト付き）|
| 4 | GitHub 等に実装コード公開（テストなし）|
| 3 | 論文補足資料にコードあり |
| 2 | 疑似コード / アルゴリズム詳細の記載あり |
| 1 | 手法説明は十分だが実装詳細なし |
| 0 | ブラックボックス / 独自ツールのみ |

### MATLAB Advantage（$S_{\text{MATLAB}}$）

| 点数 | 基準 |
|---|---|
| 5 | MATLAB ネイティブ関数（`fitlm`, `cvpartition` 等）で完結 |
| 4 | RDKit + MATLAB 標準 Toolbox のみで実装可能 |
| 3 | `installExtra()` で追加ライブラリを入れれば実装可能 |
| 2 | Track 2（venv 隔離）ライブラリで実装可能 |
| 1 | MATLAB より Python / R で実装した方が自然（移植コスト高い）|
| 0 | MATLAB での実装が現実的でない（GPU 特化等）|

### Educational Value（$S_{\text{Edu}}$）

| 点数 | 基準 |
|---|---|
| 5 | Chemoinformatics の基礎概念（QSAR, Fingerprint 等）を例示する landmark |
| 4 | 教育教材として十分な再利用性あり |
| 3 | 研究レベルだが説明しやすい手法 |
| 2 | 実用的だが教育的文脈が薄い |
| 1 | 高度すぎて初学者への説明困難 |
| 0 | 教育的意義が低い（ニッチ応用等）|

---

## 選定優先度区分

| Tier | Reproduction Score | 方針 |
|---|---|---|
| A | ≥ 4.0 | 優先着手（M-REPRO-SCALE での主要対象）|
| B | 3.0〜3.9 | 次期候補（Tier A 完了後に着手）|
| C | 2.0〜2.9 | バックログ（コミュニティ貢献候補）|
| D | < 2.0 | 対象外 |

---

## 初期重みの根拠

初期重み（0.30, 0.20, 0.20, 0.20, 0.10）は RP00-RP02 経験前の暫定値。

| 重み配分 | 根拠 |
|---|---|
| Citation × 0.30 | 再現価値の主要指標。引用数 = コミュニティの検証需要 |
| Data × 0.20 | データなしでは再現不可能。必要条件に近い |
| Code × 0.20 | 元実装参照でアルゴリズム定義が確定する |
| MATLAB × 0.20 | EasyMolKit の価値提案（MATLAB で再現できること）|
| Education × 0.10 | 副次的効果。教材量産はすでに Phase 2 で達成 |

RP00-RP02 完了後に重みを再調整する（初期値はゼロから手動で設定した暫定値）。

---

## DS03 との連携

DS03（M-DISCOVER フェーズ以降）は RF00 の定義を OpenAlex API 上で自動計算する。

- $S_{\text{Citation}}$: `cited_by_count` から自動マッピング
- $S_{\text{Data}}$, $S_{\text{Code}}$: OpenAlex の `has_fulltext`, `open_access` 等を参照
- $S_{\text{MATLAB}}$, $S_{\text{Edu}}$: 人手評価（DS03 の自動化対象外）

RF00 の定義変更は必ず本ファイル（権威ファイル）を先に改訂してから DS03 の実装を変更すること。

---

## 改訂履歴

| 日付 | 変更内容 |
|---|---|
| 2026-06-19 | M-REPRO-FOUND で初版確立（RP00 経験後） |
