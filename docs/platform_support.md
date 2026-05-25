# プラットフォームサポート

> EasyMolKit の Desktop (Windows) および MATLAB Online 対応方針。
> macOS / Linux Desktop は Deferred。

---

## 1. サポートマトリクス

| プラットフォーム | MATLAB | Python | RDKit 配備 | 優先度 | 状態 |
|---|---|---|---|---|---|
| Windows x64 (Desktop) | R2025b+ | Embedded Python 3.10 (python.org) | `emk.setup.install()` | **P0** | ⬜ |
| MATLAB Online (Linux) | R2025b+ | Built-in Python 3.10 | `emk.setup.installOnline()` | **P0** | ⬜ |
| macOS arm64 (Desktop) | R2025b+ | 未定 | 未定 | Deferred | 保留 |
| macOS x64 (Desktop) | R2025b+ | 未定 | 未定 | Deferred | 保留 |
| Linux x64 (Desktop) | R2025b+ | 未定 | 未定 | Deferred | 保留 |

> **P0** = MVP で対応必須。 **Deferred** = 将来要望に応じて対応。

---

## 2. プラットフォーム別アーキテクチャ

> 配備フローの技術詳細 → [`docs/python_integration.md`](python_integration.md)

### 2.1 Windows x64 (Desktop)

Embedded Python 3.10 を `python_env/` に zip 展開して配備。  
`python310._pth` 編集 → `get-pip.py` → `pip install rdkit-pypi`。  
詳細手順 → [`python_integration.md` §3.1](python_integration.md#31-desktop-embedded-python-windows)

**特記事項**: zip 展開のみでレジストリ・PATH 不変 / Windows Defender のブロック可能性 / 長いパス(260文字超)注意 / プロキシ環境は `settings.json` で対応

### 2.2 macOS / Linux (Desktop) — Deferred

対応時に ADR 追加が必要。現時点では MATLAB Online 経由での利用を推奨。

### 2.3 MATLAB Online

プリインストール Python 3.10 に `get-pip.py` ブートストラップ → `pip install rdkit-pypi`。  
詳細手順 → [`python_integration.md` §3.2](python_integration.md#32-online-system-python--pip)

**特記事項**: セッション揮発性で毎回 pip install 実行 / `py.sys.path` 挿入必須 / ストレージ容量制限あり

---

## 3. プラットフォーム検出フロー

> 検出ロジックの詳細 → [`python_integration.md` §3.3](python_integration.md#33-プラットフォーム検出)

検出優先順位: `ismatlabonline()` → 環境変数 `MATLAB_ONLINE` → `ispc()` + `computer("arch")`

---

## 4. プラットフォーム間の API 差異

**目標**: ユーザー向け API はプラットフォーム間で **完全に同一** とする。

| API | Desktop | Online | 差異 |
|---|---|---|---|
| `emk.setup.install()` | Embedded Python 配備 | pip install | 内部実装のみ異なる |
| `emk.mol.fromSmiles()` | 同一 | 同一 | なし |
| `emk.descriptor.*()` | 同一 | 同一 | なし |
| `emk.viz.draw2d()` | MATLAB figure | MATLAB figure | 表示方法に差異あり（※） |
| `emk.io.readSdf()` | ローカルファイル | Online ストレージ | ファイルパスの制約 |

> ※ MATLAB Online では figure の interactivity に制限がある場合あり。

---

## 5. テスト戦略

### 5.1 プラットフォーム固有テスト

| テスト | 目的 | 実行環境 |
|---|---|---|
| `TestSetup.m` | Python 配備・初期化 | 全プラットフォーム |
| `TestSetupWindows.m` | Embedded Python (Windows) | Windows |
| `TestSetupOnline.m` | pip install (Online) | MATLAB Online |
| `TestPlatformDetect.m` | プラットフォーム検出 | 全プラットフォーム |

### 5.2 CI/CD

- **GitHub Actions**: Windows / Ubuntu / macOS のマトリクスビルド（将来）
- **MATLAB Online**: 手動テスト（CI 自動化は困難）

---

## 6. 既知の制約

| # | 制約 | 影響 | 対策 |
|---|---|---|---|
| PS-1 | `pyenv(Version=...)` はセッション中 1 回のみ | 別の Python を使用中だと切替不可 | MATLAB 再起動を案内。起動直後に一度だけ呼ぶ運用を徹底 |
| PS-2 | Embedded Python は Windows 専用形式 | macOS/Linux Desktop は別方式が必要 | macOS/Linux は Deferred。対応時に ADR で代替方式を決定 |
| PS-3 | MATLAB Online の Python バージョン固定 | RDKit 互換性問題の可能性 | バージョン確認を setup に組込。Desktop も 3.10 に固定して差異を最小化 |
| PS-4 | MATLAB Online のストレージ容量 | RDKit (~200MB) で圧迫 | 最小構成での install を検討 |
| PS-5 | プロキシ環境での pip install | ダウンロード失敗 | `settings.json` の proxy 設定（単一文字列）で対応 |
| PS-6 | MATLAB Online のセッション揮発性 | インストール済み RDKit が消える可能性 | 毎回 `get-pip.py` ブートストラップ + `!~/.local/bin/pip install rdkit-pypi==<version>` を実行（既存なら高速スキップ）。`installOnline()` にバージョンチェック＋スキップロジックを実装（M1 タスク） |
| PS-7 | Windows パス長制限 (MAX_PATH 260) | `python_env/` の深い階層でパス長超過 | `install()` 冒頭でパス長チェック。200文字超で警告、240文字超でエラー（ADR-001 rev.3） |

---

## 7. 未解決事項

| # | 項目 | 説明 | 状態 |
|---|---|---|---|
| PS-Q1 | ~~macOS 配備方式の確定~~ | **Deferred**。対応時に ADR で決定 | 保留 |
| PS-Q2 | ~~Linux 配備方式の確定~~ | **Deferred**。対応時に ADR で決定 | 保留 |
| PS-Q3 | ~~MATLAB Online の Python バージョン~~ | **Python 3.10** 確認済み。Desktop も同バージョンに固定 | ✅ 確認済 |
| PS-Q4 | MATLAB Online セッション永続性 | `--user` install した RDKit は次回セッションで使えるか | ⬜ 要検証（毎回 pip install で対応済み） |
| PS-Q5 | `ismatlabonline()` 関数の有無 | R2025b で使用可能か。代替検出手段。M0-9 で実機検証予定 | ⬜ 要検証 |
| PS-Q6 | MATLAB バージョン下限 | R2025b で検証済み。下位互換は未確認 | ✅ R2025b 確定 |
| PS-Q7 | ~~Apple Silicon 対応~~ | **Deferred** (macOS 全体が保留) | 保留 |

---

## 8. MATLABライセンス別 Toolbox 収録一覧

> チュートリアル **Layer 分類の基準**。Source: MathWorks 公式 (2026-04 確認)。

| Toolbox | Online Basic（無料） | Student Suite | Campus / Individual |
|---|---|---|---|
| Statistics and Machine Learning | ✓ | ✓ | ✓ |
| Deep Learning | ✓ | ✓ | ✓ |
| Curve Fitting | ✓ | ✓ | ✓ |
| **Optimization** | **✓** | ✓ | ✓ |
| **Signal Processing** | **✓** | ✓ | ✓ |
| Image Processing | ✓ | ✓ | ✓ |
| Symbolic Math | ✓ | ✓ | ✓ |
| Control System | ✓ | ✓ | ✓ |
| Text Analytics | ✓ | — | ✓ |
| Simulink | Online 版 | Desktop 版 | ✓ |
| Parallel Computing | — | ✓ | ✓ |
| DSP System | — | ✓ | ✓ |
| Simscape / Electrical | — | ✓ | ✓ |
| SimBiology | — | — | ✓ |
| Bioinformatics | — | — | ✓ |
| Reinforcement Learning | — | — | ✓ |

### Layer 判定早見表

| Layer | 必要ライセンス | 使用 Toolbox |
|---|---|---|
| L1–L2 | Base MATLAB | なし（RDKit のみ） |
| L3 | **Online Basic（無料）** | Statistics and ML, Deep Learning, Curve Fitting, Optimization, Signal Processing |
| L4-Student | Student Suite | + Parallel Computing |
| L4-Campus | Campus / Individual | + SimBiology, Bioinformatics, Reinforcement Learning |

> ⚠️ **旧版誤記**: Optimization / Signal Processing は Campus-Wide 相当として扱われていたが、
> いずれも **Online Basic（無料）に含まれる**。L3 が正しい配置。
