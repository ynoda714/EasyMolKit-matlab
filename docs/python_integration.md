# Python 連携アーキテクチャ

> EasyMolKit における MATLAB – Python 連携の技術設計と実装方針。

---

## 1. 概要

EasyMolKit は MATLAB の `pyenv` 機能を使って Python (RDKit) を呼び出す。
ユーザーは Python の存在を意識せず、`emk.*` MATLAB 関数として利用する。

```
ユーザー → emk.mol.fromSmiles("CCO")
              → pyenv (OutOfProcess)
                  → py.rdkit.Chem.MolFromSmiles("CCO")
              ← py.rdkit.Chem.Mol (MATLAB 参照)
         ← mol (ユーザーは Python 型を意識しない)
```

---

## 2. pyenv の動作モード

| モード | 説明 | 採用 |
|---|---|---|
| **OutOfProcess** | Python を別プロセスで実行。クラッシュ耐性あり | ✅ 採用 |
| InProcess | Python を MATLAB プロセス内で実行。高速だがクラッシュで MATLAB ごと落ちる | ❌ |

**理由**: RDKit は C++ ベースで、まれにセグフォルトが発生しうる。OutOfProcess なら Python プロセスのみ落ちて MATLAB は安全。

```matlab
% pyenv 設定の基本パターン
pyenv(Version="python_env/python.exe", ExecutionMode="OutOfProcess");
```

> **制約**: MATLAB セッションで `pyenv(Version=...)` は **一度しか設定できない**。
> 変更するには MATLAB の再起動が必要。`emk.setup.initPython()` はこの制約を考慮して設計する。
> 起動直後に一度だけ呼ぶ運用を徹底する。
>
> **二重呼び出し防御**: `initPython()` 内で `pe = pyenv; if pe.Status ~= "NotLoaded"` の場合は
> 再設定をスキップし、現在の環境で RDKit が import 可能かだけ検証する。
> ユーザーが `install()` を 2 回呼んだ場合や、既に `py.xxx` を使った後の呼び出しでエラーにならないようにする。

---

## 3. Python 配備戦略

### 3.1 Desktop: Embedded Python (Windows)

python.org が公開する **Embeddable Package** を使用する。
**Python バージョン: 3.10 固定**（MATLAB Online のデフォルトに合わせる）。
venv は不使用（Embedded Python 自体が隔離された環境であるため冗長）。

```
EasyMolKit/
└─ python_env/                    ← Git 非追跡
   ├─ python.exe               ← Embedded Python 本体
   ├─ python310.dll
   ├─ python310.zip            ← 標準ライブラリ
   ├─ Lib/
   │   └─ site-packages/
   │       └─ rdkit/           ← pip install で追加
   └─ Scripts/
       └─ pip.exe
```

**特徴**:
- zip 展開のみでレジストリ・PATH・環境変数を一切変更しない。他の Python 環境への影響ゼロ
- クリーンアップはフォルダ削除のみ（アンインストーラ不要）

**配備フロー**:
1. `ispc()` で Windows を確認
2. python.org から Embeddable Package (`python-3.10.x-embed-amd64.zip`) をダウンロード
3. `python_env/` に展開
4. **`python310._pth` を編集して `import site` 行のコメントを解除**（これをしないと `site-packages` が認識されず pip install したパッケージが import できない）
5. `get-pip.py` をダウンロードして pip をブートストラップ
6. `pip install rdkit-pypi==<version>` を実行（バージョンは `settings.json` の `rdkitVersion` で固定。proxy 設定があれば `--proxy` 引数に反映）
7. `pip show rdkit-pypi` でインストール済みバージョンを照合。不一致時は警告を出力
8. `pyenv(Version="python_env/python.exe")` を設定
9. `py.importlib.import_module("rdkit")` で検証

> **パス長チェック（ADR-001 rev.3）**: ステップ 1 の前にプロジェクトルートのフルパス長を検査。
> `fullfile(pwd, "python_env", "Lib", "site-packages", "rdkit")` が 200 文字超で警告、240 文字超でエラー。
> Windows の MAX_PATH (260) 制限に起因する問題を未然防止する。

> **重要**: ステップ 4 の `_pth` ファイル編集は必須。Embeddable Package はデフォルトで `site-packages` を無効化しており、
> `python310._pth` 内の `#import site` の `#` を削除して `import site` にする必要がある。

**サイズ見積もり**: Python 本体 ~30MB + RDKit ~200MB = **合計 ~230MB**

> **macOS / Linux Desktop**: Embedded Python は Windows 専用形式のため、macOS / Linux Desktop は **Deferred**。
> 将来対応時は conda-forge / miniforge 等の代替手段を ADR で決定する。

### 3.2 Online: System Python + pip

MATLAB Online には Python 3.10 がプリインストールされている。
pip はデフォルトでは存在しないため、`get-pip.py` でブートストラップが必要。

**配備フロー**:
1. `pyenv` で Python の存在を確認
2. `websave("get-pip.py", "https://bootstrap.pypa.io/get-pip.py")` で get-pip.py を取得
3. `!python get-pip.py --user` で pip をブートストラップ
4. `!~/.local/bin/pip install rdkit-pypi==<version>` を毎回実行（バージョンは `settings.json` の `rdkitVersion` で固定。既インストールなら高速スキップ）
5. `insert(py.sys.path, int32(0), '/home/matlab/.local/lib/python3.10/site-packages')` でパスを追加
6. `py.importlib.import_module("rdkit")` で検証

> **セッション揮発性**: MATLAB Online のセッションは揮発性があるため、起動のたびに pip ブートストラップ + install を実行する。
> 既にインストール済みであれば pip は「Requirement already satisfied」で数秒で完了する。
> **重要**: `--user` でインストールしたパッケージは MATLAB の `pyenv` から自動的には見えないため、
> `py.sys.path` への明示的なパス挿入が必須（MLChem.m で実証済み）。
> **注意**: MATLAB Online の Python バージョンは MathWorks が管理。
> RDKit との互換性問題が発生する可能性があり、バージョン確認が必要。
> **pip の呼び出し規約**: `.m` コード内では `!` コマンドを使用（MATLAB Online で実証済み）。
> `system()` は戻り値が必要な場合（エラーコード取得等）にのみ使用。

### 3.3 プラットフォーム検出

```matlab
% emk.setup 内部のプラットフォーム検出ロジック（概念）
%
% 検出優先順位:
% 1. ismatlabonline() が存在すれば最も確実
% 2. 環境変数 (MATLAB_ONLINE 等) の存在チェック
% 3. ispc() + computer("arch") の組み合わせ

if exist('ismatlabonline','file') && ismatlabonline()
    % → installOnline() フロー
elseif ~isempty(getenv('MATLAB_ONLINE'))
    % → installOnline() フロー (フォールバック)
elseif ispc()
    % → Windows Embedded Python フロー
else
    % → macOS/Linux Desktop: 未サポート (Deferred)
    error('emk:setup:unsupportedPlatform', ...
        'macOS/Linux Desktop is not yet supported. Use MATLAB Online instead.');
end
```

> `ismatlabonline()` の R2025b での利用可能性は M0-9 タスクで実機検証する。
> 利用不可の場合は環境変数フォールバックを使用。

---

## 4. MATLAB ↔ Python 型変換

### 4.1 自動変換（MATLAB pyenv 組み込み）

| MATLAB 型 | Python 型 | 方向 |
|---|---|---|
| `double` | `float` / `int` | 双方向 |
| `string` | `str` | 双方向 |
| `logical` | `bool` | 双方向 |
| `cell` | `list` | MATLAB → Python |
| `struct` | — | 手動変換が必要 |

### 4.2 RDKit 固有の型変換

| RDKit 型 | EasyMolKit での扱い | 変換タイミング |
|---|---|---|
| `rdkit.Chem.Mol` | `py.rdkit.Chem.Mol` として参照保持 | 中間状態（変換しない） |
| `DataStructs.ExplicitBitVect` | 参照保持 or `logical` 配列 | 類似度計算用は参照、表示用は変換 |
| 記述子値（float） | `double` | 計算結果で即変換 |
| SMILES（str） | `string` | 即変換 |
| PIL Image（2D 構造描画） | `uint8` 配列 | Python側で PIL → `BytesIO` → PNGバイト列 → MATLAB `uint8` → 一時ファイル → `imshow()` |

**方針（ADR-002 rev.3）**: 中間オブジェクト（Mol, Fingerprint）は Python 参照のまま保持。
最終出力（数値・テーブル）段階で **`emk.util.pyToMatlab()` 統一変換レイヤ**を経由して MATLAB ネイティブ型に変換。
型変換マッピング: `py.int`→`double`, `py.float`→`double`, `py.str`→`string`, `py.list`→`cell`/`table`, `py.dict`→`struct`, `py.None`→`missing`, `py.numpy.ndarray`→`double` array。
Python の `list` / `dict` は直接使わず、`pyToMatlab()` で `table` / `struct` へ変換する。
画像（PIL Image）は Python 側で `BytesIO` 経由で PNG バイト列に変換し、MATLAB で `uint8` 配列として受け取り、
一時ファイル経由で `imshow()` で表示する（MLChem.m で実証済みパターン）。

### 4.3 Python 例外の MATLAB 変換

```matlab
% パターン: try-catch でラップ
try
    pyMol = py.rdkit.Chem.MolFromSmiles(smiles);
catch pyErr
    throwAsCaller(MException( ...
        "emk:mol:fromSmiles:pythonError", ...
        "Failed to parse SMILES: %s", string(pyErr.message)));
end
```

---

## 5. パフォーマンス考慮事項

### 5.1 OutOfProcess のオーバーヘッド

- 各 `py.xxx` 呼び出しでプロセス間通信が発生（数 ms/call）
- 単一分子操作は問題なし。大量バッチ（数千分子）では影響あり

### 5.2 バッチ処理の最適化（ADR-002 rev.3: IPC 最小化原則）

- **Python 側でバッチ実行**: 1000 分子を 1 つの Python 関数で一括処理し、結果配列を一括で返す
- **MATLAB 側でループ**: 個別に `emk.*` を呼ぶ（シンプルだが IPC コストが累積）
- **方針**: バッチ API（`batchCalculate`, `batchMorgan` 等）は Python 側でループを集約し、IPC 往復を 1 回に抑える設計を標準とする。Python ヘルパースクリプトは `src/+emk/+util/python/` に配置

### 5.3 メモリ管理

- OutOfProcess の Python オブジェクトは MATLAB 変数のスコープに紐づく
- `clear mol` で Python 側の参照も解放される
- 大量の Mol オブジェクトを保持する場合はメモリ使用量に注意

---

## 6. セキュリティ考慮事項

- Embedded Python のダウンロード元は **python.org 公式のみ**
- pip install は **PyPI 公式のみ**。カスタムインデックスは設定可能だがデフォルト無効
- `emk.setup.install()` はネットワーク接続が必要。オフラインインストールは将来対応
- ユーザー入力（SMILES 等）は RDKit に渡す前にバリデーション不要（RDKit 自体が安全にパース）
- ファイルパスは ProjectRoot 外へのアクセスを禁止

---

## 7. 未解決事項・要確認

> 実装前に回答が必要な項目。

| # | 項目 | 説明 | 状態 |
|---|---|---|---|
| PI-1 | ~~Embedded Python のバージョン方針~~ | **3.10 固定** (MATLAB Online のデフォルトに合わせる) | ✅ 決定済 |
| PI-2 | ~~macOS / Linux 対応の優先度~~ | **Deferred**。Windows Desktop + MATLAB Online を P0 とする | ✅ 決定済 |
| PI-3 | オフラインインストール | ADR-001 rev.3 で選択肢 E（vendor 化）を将来オプションとして記載。現時点では不採用 | ⬜ 将来検討 |
| PI-4 | ~~プロキシ対応~~ | `settings.json` の `proxy` 項目（単一文字列）で対応。pip の `--proxy` 引数に反映 | ✅ 決定済 |
| PI-5 | `ismatlabonline()` 互換性 | R2025b で利用可能か。代替検出方法。M0-9 で実機検証予定 | ⬜ 要検証 |
| PI-6 | Python アップデート戦略 | 配備済み Python のバージョンアップ方法 | ⬜ 要設計 |
| PI-7 | ~~MATLAB Online の Python バージョン~~ | **Python 3.10** であることを確認済み。これに合わせて Desktop も 3.10 に固定 | ✅ 確認済 |
