# Algorithm Guide — emk.setup

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## 2.1 `install()`

**設計意図**: ユーザーが Python の知識なしに RDKit 環境を構築できるようにする。

**アルゴリズム概要**:
1. Online ガード（Desktop 専用チェック）
2. パス長チェック（MAX_PATH 防御: warn > 200, error > 240 chars — ADR-001 rev.3）
3. 冪等チェック（`python.exe` 存在 + `pip show` バージョン照合 → 一致なら再インストールをスキップ）
4. Embedded Python zip ダウンロード → `python_env/` に展開
5. `python310._pth` 編集（`#import site` → `import site` にコメント解除。Embedded Python はデフォルトで site-packages が無効なため必須）
6. `get-pip.py` で pip ブートストラップ
7. `pip install rdkit-pypi==<version>` でバージョン固定インストール（proxy 設定対応）
8. バージョン照合（`pip show rdkit-pypi` の出力と `settings.json` の `rdkitVersion` を比較、不一致は警告）
9. `emk.setup.initPython()` で pyenv を設定
10. `emk.setup.verify()` で最終確認

**引用文献**:
- [Python Embeddable Package](https://docs.python.org/3/using/windows.html#the-embeddable-package) — python.org 公式ドキュメント
- ADR-001 rev.3

**テスト戦略** (`tests/unit/TestInstall.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| `test_onlineEnv_throwsNotDesktop` | `MATLAB_ONLINE=true` で `notDesktop` | なし |
| `test_onlineEnv_errorMessage_mentionsInstallOnline` | エラーメッセージに `installOnline` を含む | なし |
| `test_pathTooLong_throwsPathTooLong` | pwd + サフィックス > 240 で `pathTooLong` | なし |
| `test_pathTooLong_errorMessage_containsLength` | エラーメッセージに 3 桁の数値を含む | なし |
| `test_shortPath_doesNotThrowPathTooLong` | 短い pwd で `pathTooLong` が発生しない | なし |
| `test_unsupportedVersion_throwsDownloadFailed` | `"3.11"` 指定で `downloadFailed` | なし |
| `test_unsupportedVersion_errorMessage_containsVersion` | エラーメッセージに `"3.11"` を含む | なし |

TC1/TC2: ネットワーク不要。Online ガードとパス長チェックは `downloadPython_()` 呼び出し前に発火する。  
TC3: `resolvePatchVersion_()` が `downloadFailed` を投げるのは、ネットワーク合理が発生する前のローカル処理であることを間接的に検証。  
実際の DL/展開/pip 処理の検証は Smoke テストスイートで実施。

---

## 2.2 `installOnline(Config=struct())`

**設計意図**: MATLAB Online ユーザーが Python の知識なしに毎セッション安全に RDKit を利用できるようにする。
Online セッションは揮発性があるため、毎回呼び出しが想定される。すでに正しいバージョンがインストール済みの場合は pip 処理をスキップして高速に完了する。
`Config` パラメータを受け取ることで、`main_rdkit.m` から一貫した設定を渡し、オプションライブラリも一括インストールできる。

**アルゴリズム概要**:
1. Online ガード（Desktop 専用チェック。`emk.util.isOnline()` が false → 即 error）
2. Config 解決（`options.Config` が空なら `emkLoadConfig()` を呼ぶ。空でなければそのまま使用）
3. バージョンチェック（`pip show rdkit-pypi` の Version: フィールドを `cfg.rdkit.version` と比較。一致なら Step 4-5 をスキップ）
4. `get-pip.py` で pip をブートストラップ（`python get-pip.py --user`。pip が未インストールのセッションでも動作）
5. `pip install rdkit-pypi==<version> --user`（バージョン固定 + proxy 設定対応）
6. `emk.setup.initPython()` で pyenv を設定（`py.sys.path` アクセスの前に必須）
7. `py.sys.path().insert(int32(0), sitePackages)` でパスを挿入（MATLAB の pyenv は `--user` インストール先を自動認識しないため必須）
8. `emk.setup.verify()` で最終確認
9. オプションライブラリ（STEP 8）— `cfg.optionalLibraries.<name>=true` のエントリについて `pip install --user <pipspec>` を実行。`resolveOnlinePipName_()` で論理名 → pip 指定子を解決。

**py.sys.path 挿入が必須な根拠**:
MATLAB Online の `pyenv` は `/home/matlab/.local/lib/python3.10/site-packages`（`--user` インストール先）を
デフォルトの検索パスに含めない。明示的な挿入なしには `import rdkit` が失敗する（MLChem.m 先行事例で実証済み）。

**バージョンチェックの根拠**:
MATLAB Online のセッションは揮発性があり毎回 `installOnline()` が呼ばれる想定だが、
`get-pip.py` ダウンロード + `pip install`（既存でも数秒）のネットワーク往復を削減するため
`pip show` によるバージョン照合でスキップ判定を導入する（ADR-001 rev.3）。

**引用文献**:
- [pip install --user](https://pip.pypa.io/en/stable/cli/pip_install/#install-options) — pip 公式ドキュメント
- [get-pip.py bootstrap](https://bootstrap.pypa.io/) — PyPA 公式
- ADR-001 rev.3
- `20250513_Python_on_ML/MLChem.m` — MATLAB Online + RDKit の先行 PoC

**テスト戦略** (`tests/unit/TestInstallOnline.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| `test_desktopEnv_throwsNotOnline` | `MATLAB_ONLINE` 未設定で `notOnline` エラー | なし |
| `test_desktopEnv_errorMessage_mentionsEmkInstall` | エラーメッセージに `"emk.setup.install"` を含む | なし |
| `test_desktopEnv_errorMessage_mentionsOnline` | エラーメッセージに `"online"` を含む | なし |
| `test_onlineEnv_guardPasses_notOnlineNotThrown` | `MATLAB_ONLINE=true` のとき `notOnline` が**投げられない** | なし |
| `test_emptyConfigStruct_desktopGuardStillFires` | `Config=struct()` でも `notOnline` が投げられる（Config は Desktop ガードをバイパスしない）(TC3) | Desktop |
| `test_configWithOptionalLibs_desktopGuardStillFires` | `optionalLibraries` 付き Config でも `notOnline`（TC4）| Desktop |
| `test_configKeyword_acceptedByArgumentsBlock` | `Config=struct()` が arguments ブロックで受理される（TC5）| Desktop |

---

## 2.3 `verify()`

**設計意図**: Python 環境と RDKit の可用性をユーザーに可視化する。
非スロー設計（エラーを投げず struct を返す）とすることで、`install()` の末尾や
診断スクリプトから安全に呼び出せる。

**アルゴリズム概要**:
1. `pyenv()` で現在の Status を確認
2. Status が `NotLoaded` なら `emk.setup.initPython()` を試みる
3. Python が起動したら `result.python = true` に設定
4. `py.importlib.import_module("rdkit.Chem")` で RDKit の import を試みる
   - 成功: `result.rdkit = true`
   - 失敗（ImportError 等）: catch して `logWarn` を出力し `result.rdkit = false` を維持

**非スロー設計の根拠**:
`install()` の step 11 として呼ばれるため、verify() が例外を投げると
インストール成功後でも install() が失敗扱いになってしまう。
また、診断用途（スタンドアロン呼び出し）では呼び出し元が結果を元に判断できる方が柔軟性が高い。

**テスト戦略** (`tests/unit/TestVerify.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| `test_returnStruct_hasField_python/rdkit/version` | struct に 3 フィールドがある | なし |
| `test_pythonField_isLogical` / `test_rdkitField_isLogical` | `.python`/`.rdkit` が scalar logical | なし |
| `test_versionField_isString` | `.version` が scalar string | なし |
| `test_neverThrows` | どの Python 状態でも例外を投げない（全最重要） | なし |
| `test_initPythonFails_pythonField_isFalse` | `initPython` 失敗時に `.python=false` | NotLoaded |
| `test_initPythonFails_rdkitField_isFalse` | 同上、`.rdkit=false` | NotLoaded |
| `test_initPythonFails_versionField_isEmpty` | 同上、`.version=""` | NotLoaded |
| `test_pythonLoaded_pythonField_isTrue` | Loaded 状態で `.python=true` | Loaded |
| `test_pythonLoaded_versionField_nonEmpty` | Loaded 状態で `.version` 非空 | Loaded |

`test_neverThrows` は前提条件なしで常に実行される。非スロー設計の最重要検証。

---

## 2.4 `initPython()`

**設計意図**: `pyenv` の設定をプラットフォームに応じて自動化し、呼び出し元がプラットフォーム差異を意識しなくて済むようにする。

**アルゴリズム概要**:
1. `pyenv().Status` で二重呼び出しガード（最初に実行。Status ≠ "NotLoaded" → 冪等スキップ）
2. `emkLoadConfig()` で設定をロード
3. **Track 2 検出**（ADR-007）: `cfg.python.external_path` が非空かつ非空白 → `emk.setup.useExternal(extPath)` を呼び出してリターン
   - `EMK_PYTHON_EXTERNAL_PATH` 環境変数（または `settings.json` の `python.external_path`）で有効化
   - 相対パスの場合はプロジェクトルートを基準に絶対パスへ変換
4. **Track 1**（Embedded Python）: `python_env/python.exe` の存在確認 → 不在なら `notInstalled`
5. MATLAB Online: ファイル確認をスキップし、システム Python で `pyenv` を設定
6. `pyenv(Version=pyPath, ExecutionMode=execMode)` でセッションを初期化

**二重呼び出し防御の根拠**:
MATLAB の `pyenv` はセッション中に一度しか `Version`/`ExecutionMode` を変更できない（Python が "Loaded" 状態になると設定変更不可）。`install()` の末尾や `main_rdkit.m` の Section 0a から繰り返し呼ばれる可能性があるため、冪等性を保証するガードが必要（ADR-001 rev.3）。

**実装上の制約（ガード優先の根拠）**:
ガードは関数の**最初**に置く必要がある。理由: Python が `Loaded` 状態でも `python_env/python.exe` が存在しない
デスクトップ環境（例: 別 PC で再実行、`install()` 後に `python_env/` を削除した場合）では
ファイルチェックが先に来ると `notInstalled` が誤発生する。冪等性を保証するにはガード最優先が必須。

**テスト戦略** (`tests/unit/TestInitPython.m`):

| テスト名 | 状態条件 | 検証内容 |
|---|---|---|
| `test_desktop_noPythonEnv_throwsNotInstalled` | NotLoaded + Desktop | エラーID = `notInstalled` |
| `test_desktop_noPythonEnv_errorMessage_containsPath` | 同上 | メッセージに `python_env` パスを含む |
| `test_desktop_noPythonEnv_errorMessage_mentionsInstall` | 同上 | メッセージに `install` ヒントを含む |
| `test_guardFiresBeforeFileCheck_noNotInstalledWhenLoaded` | **Loaded + Desktop + 不在** | ガードが先に発動し `notInstalled` を投げない（最重要）|
| `test_alreadyLoaded_desktop_silentReturn` | Loaded + Desktop | エラーなし・pyenv Status 変化なし |
| `test_alreadyLoaded_online_silentReturn` | Loaded + Online | エラーなし |
| `test_online_noPythonEnv_notInstalledNeverThrown` | any + Online + 不在 | `notInstalled` を投げない（`pyenvFailed` は許容）|
| `test_externalPathSet_nonEmpty_desktop_callsUseExternal` (**EXT1**) | NotLoaded + Desktop + `external_path` 非空 | `useExternal` エラーに到達、`notInstalled` が出ない（Track 2 パスを確認）|
| `test_externalPathEmpty_desktop_fallsBackToEmbedded` (**EXT2**) | NotLoaded + Desktop + `external_path` 空 | `notInstalled` が出る（Track 1 フォールバック確認）|
| `test_externalPathWhitespace_desktop_treatedAsEmpty` (**EXT3**) | NotLoaded + Desktop + `external_path` 空白のみ | `notInstalled` が出る（`strtrim` による空判定を確認）|

EXT1–EXT3 の前提: `EMK_PYTHON_EMBEDDED_DIR` に不在ディレクトリを設定して Track 1 の影響を排除。`EMK_PYTHON_EXTERNAL_PATH` で `external_path` を注入（`emkLoadConfig` の `buildDefaults()` に `external_path=""` を追加することで env var オーバーライドが有効化される）。

---

## 2.5 `installExtra(name)`

**設計意図**: Track 1 オプションライブラリ（PubChemPy / mordred / Biopython / torch 系）を Embedded Python に追加インストールし、ユーザーが pip を意識せずに使えるようにする。

**アルゴリズム概要**:
1. Online ガード（Desktop 専用チェック）
2. **Windows 早期 return**（`ispc()=true` のとき、`vina` と `pdbfixer` は pip より前に `logWarn` + `recipe("docking")` を呼んで `return`）
   - `vina`: Windows 向け PyPI ホイールが存在しない（Boost C++ ビルドが必要）
   - `pdbfixer`: openmm の `.pyd` が Windows Smart App Control (SAC) にブロックされる
   - `meeko` は Windows でも PyPI ホイールが存在するため早期 return なし
3. `resolveLibInfo_()` でライブラリ名 → pip パッケージ名・import 名を解決（不明なら `unknownLibrary` エラー）
4. `mfilename("fullpath")` から 4 階層上のプロジェクトルートを計算し `python_env/python.exe` を特定（ADR-005）
5. `cfg.extraLibraries.<name>` からバージョンを取得。未定義ならバージョン未固定でインストール
6. `python.exe -m pip install <spec> --no-warn-script-location` を実行
   - `torch_geometric` は特殊フロー（下記 §2.5.1 参照）でこのステップを置き換える
7. `python.exe -c "import <importName>"` でインポート検証

**ライブラリ名マッピング**:

| 名前 (name) | pip パッケージ | import 名 | ライセンス |
|---|---|---|---|
| `pubchempy` | `pubchempy` | `pubchempy` | MIT |
| `mordred` | `mordredcommunity` | `mordred` | BSD-3 |
| `biopython` | `biopython` | `Bio` | Biopython License |
| `torch` | `torch` | `torch` | BSD-3 |
| `torch_geometric` | `torch_geometric` + extras | `torch_geometric` | MIT |
| `transformers` | `transformers` | `transformers` | Apache-2.0 |
| `meeko` | `meeko` | `meeko` | **LGPL-2.1** (動的インポートは MIT に安全, CL-7 参照) |
| `vina` | `vina` | `vina` | Apache-2.0（Windows は PyPI ホイールなし → 早期 return）|
| `pdbfixer` | `pdbfixer` | `pdbfixer` | MIT（openmm >= 8.2 を自動依存として引き込む；Windows SAC にブロック → 早期 return）|

### §2.5.1 `torch_geometric` 動的 URL 生成

**課題**: PyG (PyTorch Geometric) のホイールは PyTorch バージョンと ABI に紐付いており、
PyPI に公開されていない。`pip install torch_geometric` だけでは依存パッケージが解決できず失敗する。

**アルゴリズム**:
1. `python.exe -c "import torch; print(torch.__version__)"` で torch バージョン文字列を取得
2. `+cpu` などのビルドサフィックスを除去（例: `"2.5.0+cpu"` → `"2.5.0"`）
3. PyG の CPU 版 find-links URL を構築:
   ```
   https://data.pyg.org/whl/torch-<X.Y.Z>+cpu.html
   ```
4. `torch_geometric`, `torch_scatter`, `torch_sparse`, `torch_cluster` を一括インストール:
   ```
   python.exe -m pip install torch_geometric torch_scatter torch_sparse torch_cluster
       --no-warn-script-location -f <pygUrl>
   ```
5. `import torch_geometric` でインポート検証

**companion パッケージが必要な理由**:
`torch_scatter`, `torch_sparse`, `torch_cluster` は PyG のコアアルゴリズム（メッセージパッシング・スパース行列演算・グラフクラスタリング）に必要。PyPI ではソースのみ公開されており、コンパイル済みホイールは PyG の find-links URL から取得する必要がある。R08（GNN 分子性質予測）の完全動作に必須。

**バージョン固定方針**:
`config/settings.example.json` の `extraLibraries.torch` に `"2.5.0"` を設定することで、
torch バージョンを決定論的に管理し、PyG URL を再現可能にする。`torch_geometric` 自体のバージョン固定も同フィールドで対応可能。

**引用文献**:
- [PyTorch Geometric Installation Guide](https://pytorch-geometric.readthedocs.io/en/latest/install/installation.html) — PyTorch Geometric 公式（動的 URL パターンの根拠）
- ADR-006 rev.2（torch 系 Track 1 追加・PyG 動的 URL 構築の設計決定）

**テスト戦略** (`tests/unit/TestInstallExtra.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| `test_onlineEnv_throwsNotDesktop` | Online で `notDesktop` エラー | なし |
| `test_onlineEnv_errorMessage_guidesPipInstall` | `notDesktop` メッセージに pip install ガイダンスを含む | なし |
| `test_unknownLibrary_throwsUnknownLibrary` | 不明名で `unknownLibrary` エラー | なし |
| `test_unknownLibrary_errorMessage_containsName` | エラーに渡した名前を含む | なし |
| `test_unknownLibrary_errorMessage_mentionesSupportedNames` | Track 1 の名前が 1 つ以上含まれる | なし |
| `test_unknownLibrary_errorMessage_mentionsRecipe` | Track 2 への誘導として `recipe` が含まれる | なし |
| `test_noPythonEnv_throwsInstallFailed` | `python.exe` 不在で `installFailed`（`EMK_PYTHON_EMBEDDED_DIR` で模擬）| Desktop |
| `test_noPythonEnv_errorMessage_containsPath` | `installFailed` メッセージにパスを含む | Desktop |
| `test_noPythonEnv_errorMessage_mentionsInstall` | `installFailed` メッセージに `install` ヒントを含む | Desktop |
| `test_charInput_acceptedByArguments` | char 型入力が `arguments` ブロックで string に変換される | なし |
| `test_newLibNames_areKnown_notUnknownLibraryError` | torch/torch_geometric/transformers が `unknownLibrary` を投げない（installFailed を投げる）| Desktop |
| `test_torchGeometric_onlineEnv_throwsNotDesktop` | Online で `torch_geometric` が `notDesktop` を投げる（`installTorchGeometric_()` 到達前に止まる）| なし |
| `test_settingsExample_torch_versionFieldExists` | `settings.example.json` に `extraLibraries.torch` フィールドが存在する | なし |
| `test_settingsExample_torch_notDateBased` | torch バージョンが X.Y.Z semver（日付ベースでない）| なし |
| `test_unknownLibrary_errorMessage_mentionsTorchName` | `unknownLibrary` メッセージに `"torch"` が含まれる（A-1 追加名の案内）| なし |
| `test_meeko_isKnown_throwsInstallFailed` (TC10) | meeko が `installFailed` を投げる（`unknownLibrary` でない。A-3 登録確認）| Desktop |
| `test_dockingLibNames_onlineEnv_throwNotDesktop` (TC10b) | meeko/vina/pdbfixer 各々が Online で `notDesktop` を投げる | なし |
| `test_vina_windowsDesktop_returnsNormally` (TC10c) | vina が Windows Desktop で正常 return（エラーなし。PyPI ホイール不在 → 早期 return）| Desktop, Windows |
| `test_pdbfixer_windowsDesktop_returnsNormally` (TC10d) | pdbfixer が Windows Desktop で正常 return（エラーなし。openmm SAC ブロック → 早期 return）| Desktop, Windows |
| `test_vina_nonWindowsDesktop_isKnown_throwsInstallFailed` (TC10e) | vina が非 Windows Desktop で `installFailed`（Python 不在）| Desktop, 非 Windows |
| `test_pdbfixer_nonWindowsDesktop_isKnown_throwsInstallFailed` (TC10f) | pdbfixer が非 Windows Desktop で `installFailed`（Python 不在）| Desktop, 非 Windows |
| `test_unknownLibrary_errorMessage_mentionsDockingName` | `unknownLibrary` メッセージに meeko・vina・pdbfixer の**全て**が含まれる（個別 `verifyTrue`）| なし |
| `test_meeko_lgplWarning_loggedBeforeInstall` | `installExtra("meeko")` が `installFailed` に到達する（LGPL 警告が pip ステップ前にあることの間接確認）| Desktop |

---

## 2.6 `useExternal(pythonPath)`

**設計意図**: Track 2 ライブラリ（Open Babel, MDAnalysis, PyMOL-OSS 等）のために外部 CPython 環境に `pyenv` を切り替え、EasyMolKit の統一 API でアクセスできるようにする。

**アルゴリズム概要**:
1. 入力検証（string/char 型チェック・空文字ガード）
2. `pyenv().Status` が `NotLoaded` でなければ冪等警告のみで返す（pyenv 制約: 一度 Loaded になると変更不可）
3. `isfile(pythonPath)` で実行ファイルの存在確認
4. `pyenv(Version=pythonPath, ExecutionMode=cfg.python.execution_mode)` で外部 Python を設定

**冪等ガードの根拠**: `pyenv(Version=...)` はセッション中に一度しか有効でないため、すでに Loaded な場合に再実行してもエラーになるだけ。警告のみで処理を続行する方がユーザー体験が良い（ADR-001 rev.3 と同様の設計）。

**テスト戦略** (`tests/unit/TestUseExternal.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| `test_numericInput_throwsInvalidInput` | 数値入力で `invalidInput` | なし |
| `test_cellInput_throwsInvalidInput` | cell 配列入力で `invalidInput` | なし |
| `test_emptyString_throwsInvalidInput` | 空文字 `""` で `invalidInput` | なし |
| `test_whitespaceOnlyString_throwsInvalidInput` | スペースのみ文字列で `invalidInput`（strtrim ガード） | なし |
| `test_emptyChar_throwsInvalidInput` | 空 char `''` で `invalidInput` | なし |
| `test_nonexistentPath_throwsFileNotFound` | 不在パスで `fileNotFound`（NotLoaded 時）| NotLoaded（assumeTrue）|
| `test_nonexistentPath_errorMessage_containsPath` | `fileNotFound` メッセージにパスを含む | NotLoaded |
| `test_nonexistentPath_errorMessage_guidesUser` | メッセージに Python/install 文言を含む | NotLoaded |
| `test_charPath_acceptedAsString` | char パスが受け入れられ `invalidInput` が出ない | なし |
| `test_alreadyLoaded_noThrow` | Python Loaded 時に例外なし・冪等スキップ | Loaded |

---

## 2.7 `installTrack2(name)`

**設計意図**: Track 2 ライブラリ（MDAnalysis / PyMOL-OSS）を 1 コマンドで導入できるようにする。
`python -m venv` + `pip install` の手動作業をラッパーが担い、`settings.json` に接続先を永続化することで
次回セッションから `initPython()` が自動検出する（ADR-007）。

**アルゴリズム概要**:
1. ライブラリ名バリデーション（`unknownLibrary` を最初に発火させ精度の高いエラーを返す）
2. Desktop-only ガード（`emk.util.isOnline()` が true → `notDesktop`）
3. Base Python 解決 — `BasePython` 引数 > `py -3` ランチャー > `python` in PATH の順で検索
4. venv ディレクトリ決定: `<projectRoot>/python_env_t2/<name>/`（`mfilename` ベース。ADR-005 準拠）
5. GPL ライセンス警告をログに出力（`logWarn`）
6. venv 作成（冪等: `Scripts/python.exe` が既存なら SKIP）: `basePy -m venv venvDir`
7. pip install（venv の pip を使用）
8. import 検証（`venvPy -c "import <name>"`）
9. `config/settings.json` に `python.external_path = venvPy` を書込（存在しなければ新規作成）
10. `emk.setup.useExternal(venvPy)` で接続

**次回セッションの自動検出フロー** (`initPython()` の拡張):
- `cfg.python.external_path` が非空の場合 → `useExternal(extPath)` を呼び出してリターン
- Track 1（Embedded Python）は `external_path` が空の場合にのみ使用される

**venv 冪等化の根拠**: `installTrack2()` は繰り返し呼び出されても安全でなければならない。
`Scripts/python.exe` が存在する場合は venv 作成をスキップし pip install のみ実行する。
これによりバージョン更新（`pip install --upgrade <pkg>`）にも対応できる。

**GPL 汚染回避の根拠**: EasyMolKit のソースコードには GPL コードを一切含まない。
Track 2 ライブラリはユーザーの外部 venv にインストールされるため、EasyMolKit MIT ライセンスへの汚染は生じない。
ただし `logWarn` でユーザーに GPL ライセンスを通知し、利用判断を促す（ADR-006 / ADR-007）。

**引用文献**:
- [Python venv](https://docs.python.org/3/library/venv.html) — python.org 公式ドキュメント
- [MDAnalysis Installation](https://docs.mdanalysis.org/stable/documentation_pages/installation.html) — MDAnalysis 公式
- [PyMOL Open-Source](https://github.com/schrodinger/pymol-open-source) — Schrodinger GitHub
- ADR-006 / ADR-007

**テスト戦略** (`tests/unit/TestInstallTrack2.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| `test_onlineEnv_throwsNotDesktop` | Online で `notDesktop` | なし |
| `test_onlineEnv_errorMessage_guidesPipInstall` | エラーメッセージに `pip` を含む | なし |
| `test_onlineEnv_notDesktopMessage_containsLibPipPackage_mdanalysis` (**TC1c**) | `notDesktop` メッセージに `"MDAnalysis"` を含む（pip パッケージ名） | なし |
| `test_onlineEnv_notDesktopMessage_containsLibPipPackage_pymol` (**TC1d**) | `notDesktop` メッセージに `"pymol-open-source"` を含む | なし |
| `test_unknownLibrary_throwsUnknownLibrary` | 未対応名で `unknownLibrary` | なし |
| `test_unknownLibrary_errorMessage_containsName` | エラーメッセージに不正名を含む | なし |
| `test_unknownLibrary_errorMessage_listsSupportedNames` | サポート名（mdanalysis/pymol）を列挙 | なし |
| `test_unknownLibrary_errorMessage_mentionsOpenBabelRecipe` | `recipe()` への案内を含む | なし |
| `test_unknownLibrary_firesBeforeNotDesktop_onOnline` (**TC2e**) | MATLAB Online + 不正名 → `unknownLibrary`（`notDesktop` より先に発火） | なし |
| `test_mdanalysis_isRecognized_notUnknownLibrary` | `"mdanalysis"` が `notDesktop` まで進む | なし |
| `test_pymol_isRecognized_notUnknownLibrary` | `"pymol"` が `notDesktop` まで進む | なし |
| `test_charInput_isCoerced_reachesUnknownLibrary` | `'char'` リテラル（単引用符）が coerce されて `unknownLibrary` へ | なし |
| `test_nonExistentBasePython_throwsBasePythonNotFound` | 不在パスで `basePythonNotFound` | Desktop |
| `test_nonExistentBasePython_errorMessage_containsPath` | エラーメッセージにパスを含む | Desktop |
| `test_settingsExample_hasPythonExternalPath` | `settings.example.json` に `external_path` フィールドが存在 | なし |
| `test_pymolOnlineEnv_throwsNotDesktop` | `"pymol"` で Online 時に `notDesktop` | なし |

バリデーション順: `unknownLibrary`（RDKit 不要）→ `notDesktop`（RDKit 不要）→ `basePythonNotFound` の順で発火する。
venv 作成・pip install は実際のネットワーク・Python 環境が必要なため Smoke テストスイートで実施する。

---

## 2.7 `validate(Libraries=[])`

**設計意図**: 現在の Python 環境に必要なライブラリが揃っているかを一覧形式で診断し、Track 1 / Track 2 の区分と共に提示する。非スロー設計により `install()` の後処理や診断スクリプトから安全に呼べる。

**アルゴリズム概要**:
1. `Libraries` が指定されない場合、既知の全ライブラリリスト（rdkit / pubchempy / mordred / biopython / openbabel / mdanalysis / pymol）を使用
2. `resolvePythonExe_()` で現在の Python 実行ファイルを取得（非スロー; initPython 試行 → embedded fallback → "" の順）
3. 各ライブラリについて `python -m pip show <pip-show-name>` を実行し、Version フィールドをパース
4. `table(Library, Installed, Version, Track)` を返す

**pip show 名マッピング**:

| name | pip show 名 | Track |
|---|---|---|
| `rdkit` | `rdkit` | 1 |
| `pubchempy` | `pubchempy` | 1 |
| `mordred` | `mordredcommunity` | 1 |
| `biopython` | `biopython` | 1 |
| `openbabel` | `openbabel` | 2 |
| `mdanalysis` | `MDAnalysis` | 2 |
| `pymol` | `pymol` | 2 |

**非スロー設計の根拠**: Python が未設定の場合でも `Installed=false` の行を返すだけで済み、診断ツールとしての利用に支障がない。

**テスト戦略** (`tests/unit/TestValidate.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| `test_returnType_isTable` | 戻り値が table | なし |
| `test_returnTable_hasColumn_Library` / `_Installed` / `_Version` / `_Track` | 各列が存在 | なし |
| `test_returnTable_hasExactly4Columns` | 列数が正確に 4 | なし |
| `test_installedColumn_isLogical` | Installed が logical（数値ではない）| なし |
| `test_libraryColumn_isString` / `test_versionColumn_isString` / `test_trackColumn_isString` | 各文字列列の型 | なし |
| `test_trackColumn_valuesAreKnownTokens` | Track 値が `"1"` / `"2"` / `"?"` のみ | なし |
| `test_invalidLibraries_nonString_throwsInvalidInput` | 数値 Libraries で `invalidInput` | なし |
| `test_invalidLibraries_cellArray_throwsInvalidInput` | cell 配列で `invalidInput` | なし |
| `test_defaultList_containsRdkitRow` | デフォルトで rdkit 行を含む | なし |
| `test_defaultList_returns7Rows` | デフォルトで 7 行 | なし |
| `test_defaultList_containsAllKnownLibraries` | 既知 7 名が全て含まれる | なし |
| `test_customLibraries_rowCountMatchesInput` | 指定数と返り行数が一致 | なし |
| `test_customLibraries_rowOrderMatchesInput` | 返り行の順序が入力順と一致 | なし |
| `test_customLibraries_onlyRequestedLibrariesReturned` | 余分な行が追加されない | なし |
| `test_neverThrows` / `test_customLibraries_neverThrows` | 例外を投げない（非スロー設計）| なし |
| `test_emptyStringArray_usesDefaultList` | `string.empty` はデフォルトリスト（7 行）と同等 | なし |

---

## 2.8 `recipe(name)`

**設計意図**: Track 2 ライブラリの手動インストール手順（コマンド・ライセンス・注意事項）を MATLAB コマンドウィンドウに表示し、ドキュメントを開かなくても作業を進められるようにする。Track 1 ライブラリは `installExtra()` の呼び出し方を案内する。

**アルゴリズム概要**:
1. `name` を switch-case で照合し、Track 1 / Track 2 用のレシピ情報を選択
2. `printTrack1_()` または `printTrack2_()` で `logInfo()` を使い整形表示
3. 不明な名前は `unknownLibrary` エラー

**GPL ライセンス表示ポリシー**: Track 2 の GPL ライブラリ（openbabel / mdanalysis）については、EasyMolKit 本体（MIT）とは独立している旨を明示的に Notes に含める（ADR-006 rev.1）。

**テスト戦略** (`tests/unit/TestRecipe.m`):

| テスト名 | 検証内容 | 前提条件 |
|---|---|---|
| `test_unknownLibrary_throwsUnknownLibrary` | 不明名で `unknownLibrary` エラー | なし |
| `test_unknownLibrary_errorMessage_containsName` | エラーメッセージに渡した名前を含む | なし |
| `test_unknownLibrary_errorMessage_listsSupportedNames` | エラーに Track 1 の名前が 1 つ以上含まれる | なし |
| `test_pubchempy_noThrow` / `test_mordred_noThrow` / `test_biopython_noThrow` | Track 1 各名で例外なし | なし |
| `test_openbabel_noThrow` / `test_mdanalysis_noThrow` / `test_pymol_noThrow` | Track 2 各名で例外なし | なし |
| `test_allCanonicalNames_noThrow` | 全サポート名をループで確認（リグレッションガード）| なし |
| `test_charInput_unknownLibrary_throwsUnknownLibrary` | char 入力が `arguments` で string に変換される | なし |
| `test_meeko_noThrow` | meeko (A-3 Track 1) で例外なし | なし |
| `test_meeko_recipe_mentionsLGPL` | meeko case ブロックに `"LGPL"` の文字を含む（ソース静的検証）| なし |
| `test_meeko_recipe_mentionsComplianceDoc` | meeko case ブロックに `"compliance"` の文字を含む（CL-7 参照確認）| なし |
| `test_meeko_recipe_mentionsTrack1` | meeko case が `printTrack1_()` を呼ぶ（Track 1 分類確認）| なし |
| `test_vina_noThrow` | vina (A-3 Track 1) で例外なし | なし |
| `test_pdbfixer_noThrow` | pdbfixer (A-3 Track 1) で例外なし | なし |
| `test_pdbfixer_recipe_mentionsOpenmm` | pdbfixer case ブロックに `"openmm"` を含む（自動依存の通知確認）| なし |
| `test_docking_noThrow` | docking (combined A-3 recipe) で例外なし | なし |
| `test_docking_recipe_mentionsAllThreeComponents` | `printDocking_()` 関数本体に meeko・vina・pdbfixer の**全て**を含む（個別 `verifyTrue`）| なし |

