# Log Format

## 1. 目的
パイプラインの実行ログを run 単位で統一し、再現性・監査性を担保する。

---

## 2. コンソールログ書式

```
[HH:MM:SS][INFO]  通常情報メッセージ
[HH:MM:SS][WARN]  警告メッセージ
[HH:MM:SS][ERROR] エラーメッセージ
[HH:MM:SS][DEBUG] デバッグ詳細（VERBOSE時のみ）
[####------]  40% ( 4/10) ループ進捗ラベル
```

### ヘルパー関数（`src/util/` に配置）

| 関数 | 用途 |
|---|---|
| `logInfo(msg, ...)` | 通常情報 |
| `logWarn(msg, ...)` | 警告 |
| `logError(msg, ...)` | エラー |
| `logDebug(msg, ...)` | デバッグ詳細 (`EMK_LOG_VERBOSE=1` 時のみ出力) |
| `logProgress(i, n, label)` | ループ進捗バー（`\r` 上書き） |
| `logSection(scriptId, label, layer)` | チュートリアルセクション開始バナー |

`logSection` 出力例:
```
[11:36:40][INFO]  --- R01 | Section 0: Setup  [Research L4] ---
```
チュートリアルスクリプト（`examples/`）の各 `%%` セクション先頭で呼び出す。引数: `scriptId`="R01"、`label`="Section 0: Setup"、`layer`="Research L4"。

- `msg` は `sprintf` 書式対応（`logInfo("rows=%d", n)` のように可変引数可）
- `fprintf` 直書きは原則禁止
- `.m` ファイル内のログメッセージは英語のみ

---

## 3. ファイル構成（実行単位）

`result/runs/<YYYYMMDD_HHMMSS>/` に以下を保存する。

- **主成果物 CSV** — run ごとの出力
- `run_meta.json` — 実行メタデータ（JSONログ）

注記:
- `result/` はGit追跡対象外（`.gitkeep` のみ追跡）
- `*.csv` / `*.xlsx` もGit追跡対象外

---

## 4. run_meta.json 仕様

### 4.1 必須キー

| key | type | 説明 |
|---|---|---|
| `run_id` | string | 実行ID（`yyyyMMdd_HHmmss` 推奨） |
| `run_timestamp` | string | 実行時刻（ISO 8601） |
| `eval_mode` | string | `prod` / `ref` |
| `input_file` | string | 入力ファイルパス |
| `output_file` | string | 主成果物パス |

### 4.2 推奨キー

| key | type | 説明 |
|---|---|---|
| `row_count_input` | int | 入力行数 |
| `row_count_output` | int | 出力行数 |
| `elapsed_sec` | float | 実行時間（秒） |
| `host` | string | 実行ホスト識別 |
| `matlab_version` | string | MATLAB バージョン |

---

## 5. 監査カウンタ（任意）

品質ゲートが必要な場合、`run_meta.json` に以下を追加する。

| key | 説明 |
|---|---|
| `invalid_row_count` | 無効行数（目標: 0） |
| `error_count` | エラー発生件数（目標: 0） |
| `warn_count` | 警告発生件数（参考値） |
