# Algorithm Guide — emk.util / src/util / src/config

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## 9.1 `pyToMatlab`

**設計意図**: Python → MATLAB 型変換の統一レイヤ（ADR-002 rev.3）。

**型変換マッピング**:

| Python 型 | MATLAB 型 |
|---|---|
| `py.int` | `double` |
| `py.float` | `double` |
| `py.str` | `string` |
| `py.list` | `cell` または `table` |
| `py.dict` | `struct` |
| `py.None` | `missing` |
| `py.numpy.ndarray` | `double` array |

**テスト戦略**:
- 各型変換の正確性テスト（入力→出力の型・値一致）
- `py.None` → `missing` の特殊ケース

---

## 9.2 `emk.util.isOnline`

**設計意図**: Desktop / MATLAB Online を自動判別し、`emk.setup.*` が適切なインストールパスを選択できるようにする（`platform_support.md` §3 参照）。

**検出アルゴリズム（3 段階）**:

| 優先順位 | 手段 | 根拠 |
|---|---|---|
| 1 | `ismatlabonline()` 組み込み関数 | R2023b+ で公式サポート。最も確実 |
| 2 | 環境変数 `MATLAB_ONLINE` | MATLAB Online インフラが自動設定。テスト時の注入にも使用 |
| 3 | `~ispc() && computer("arch")=="glnxa64"` heuristic | `ismatlabonline` 未搭載の旧 MATLAB 向け。Linux Desktop との区別は不能 |

**テスト戦略**:
- env var `MATLAB_ONLINE=true/1/false` による分岐を Unit テストでカバー
- Windows Desktop では `MATLAB_ONLINE` 未設定時に `false` を返すことを確認
- 実際の MATLAB Online での動作確認は M1 スモークテスト時に実施

---

## 9.3 `emkLoadConfig`

**設計意図**: 設定の優先順位（環境変数 > JSON > デフォルト）を一元管理し、全モジュールが同一の設定構造体を参照できるようにする。

**アルゴリズム概要**:
1. `buildDefaults()` でハードコードされたデフォルト値を生成
2. `config/settings.json` が存在すればパース → `mergeStruct` でデフォルトに上書き（未知キーは無視）
3. `EMK_<SECTION>_<KEY>` パターンの環境変数をスキャン → 型推論（logical / numeric / string）で上書き

**boolean 判定ルール**: `"true"` または `"1"` で `true`。`strcmp` を使用（`getenv` は常に char を返すため `isequal` との型不一致を回避）。

**テスト戦略**:
- デフォルト値の構造確認
- env var による文字列・boolean・数値の上書き
- settings.json 読み込みと未知キー無視の確認

---

## 9.4 ログヘルパー群 (`logInfo`, `logWarn`, `logError`, `logDebug`, `logProgress`)

**設計意図**: `fprintf` 直書きを排除し、ログレベル・タイムスタンプ・verbose 制御を統一する（`docs/log_format.md` 参照）。

**出力形式**: `[HH:MM:SS][LEVEL]  message`

**verbose 制御**: `logDebug` は環境変数 `EMK_LOG_VERBOSE=1` 時のみ出力。`strcmp` でチェック（char との型一致保証）。

**テスト戦略**:
- 各レベルのタイムスタンプ形式・レベルラベルの出力確認
- `logDebug` の verbose ON/OFF 切り替え
- `logProgress` の割合計算（25%/100%）と `n=1` エッジケース

---

## 9.5 `makeRunDir`

**設計意図**: 成果物ディレクトリのパス生成と作成を一元化し、固定パスや `mkdir` 直書きを排除する。

**アルゴリズム概要**:
- `datetime` の `yyyyMMdd_HHmmss` フォーマットでタイムスタンプを生成
- `Prefix` が指定された場合、`<timestamp>_<prefix>` をディレクトリ名とする
- `BaseDir` 配下に `mkdir` で作成（中間ディレクトリも含めて作成）

**テスト戦略**:
- Prefix の有無でディレクトリ名パターンを確認
- 連続 2 回の呼び出しで異なるパスが返ること（`Prefix` を変えて一意性を担保）
