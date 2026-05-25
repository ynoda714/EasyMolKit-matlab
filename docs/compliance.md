# ライセンス・コンプライアンス

> EasyMolKit が依存する OSS のライセンス要件と再配布条件の整理。
> 実装・配布前に本文書を確認し、要件を満たすこと。

---

## 1. EasyMolKit 本体のライセンス

| 項目 | 内容 |
|---|---|
| ライセンス | MIT License（予定） |
| 理由 | 学術 MATLAB ツールの標準。制約が最小で広く採用されている |

> **TODO**: `LICENSE` ファイルをルートに配置する（M-PUB で対応）

---

## 2. サードパーティライセンス一覧

### 2.1 Python (CPython)

| 項目 | 内容 |
|---|---|
| ライセンス | PSF License (Python Software Foundation License) |
| 種別 | Permissive (BSD 類似) |
| ソース | https://www.python.org/ |
| 再配布条件 | ライセンス表示・著作権表示の保持 |
| EasyMolKit での利用 | Embedded Python として `python_env/` に配備（Desktop） |

**再配布時の義務**:
- Python の LICENSE ファイルを含めること
- 著作権表示を `THIRD_PARTY_NOTICES.md` に記載すること
- Embedded Python は "redistribution" に該当するため、PSF License の表示が必須

> **参照**: https://docs.python.org/3/license.html

### 2.2 RDKit

| 項目 | 内容 |
|---|---|
| ライセンス | BSD-3-Clause |
| ソース | https://github.com/rdkit/rdkit |
| 再配布条件 | ライセンス表示・著作権表示・免責事項の保持 |
| EasyMolKit での利用 | pip install で配備。MATLAB からラッパー経由で呼び出し |

**再配布時の義務**:
- BSD-3-Clause ライセンス文を `THIRD_PARTY_NOTICES.md` に記載
- RDKit の著作権者名を明記
- **制限**: RDKit の名称を EasyMolKit の宣伝目的に使用しない（BSD-3 の第三条項）

> **参照**: https://github.com/rdkit/rdkit/blob/master/license.txt

### 2.3 PyMOL Open-Source（将来対応）

| 項目 | 内容 |
|---|---|
| ライセンス | Python-2.0 (PSF-like) / BSD-like |
| ソース | https://github.com/schrodinger/pymol-open-source |
| 再配布条件 | ライセンス表示の保持 |
| EasyMolKit での利用 | Phase 2 で 3D 構造可視化に利用予定 |

> **注意**: PyMOL には **商用版**（Schrödinger 社が販売）と **OSS 版** がある。
> EasyMolKit は OSS 版のみを対象とする。商用版のライセンスとの混同に注意。

### 2.4 ChEMBL データ（サンプルデータセット）

| 項目 | 内容 |
|---|---|
| ライセンス | **CC-BY-SA 3.0** (Creative Commons Attribution-ShareAlike 3.0 Unported) |
| ソース | https://www.ebi.ac.uk/chembl/ |
| 対象ファイル | `data/list/fda_drugs.csv`, `data/list/forensic_challenge.csv` |
| 再配布条件 | 帰属表示（ChEMBL / EMBL-EBI）必須 + **ShareAlike**: 派生データを同じ CC-BY-SA 3.0 で配布 |

**重要: ShareAlike 条件**:
- `fda_drugs.csv` / `forensic_challenge.csv` の派生物（加工・フィルタリングしたデータ等）を配布する場合、CC-BY-SA 3.0 ライセンス下で配布する必要がある
- EasyMolKit 本体（MIT）のコードとは **別ライセンス** として README / THIRD_PARTY_NOTICES.md に明示すること
- 参照: https://chembl.gitbook.io/chembl-interface-documentation/about#data-licensing

### 2.5 PubChem データ（サンプルデータセット）

| 項目 | 内容 |
|---|---|
| ライセンス | CC0（パブリックドメイン） |
| ソース | https://pubchem.ncbi.nlm.nih.gov/ |
| 対象ファイル | `data/list/everyday_chemicals.csv` |
| 再配布条件 | 制限なし |

### 2.6 pip / setuptools

| 項目 | 内容 |
|---|---|
| ライセンス | MIT License |
| 利用 | Embedded Python への RDKit インストール手段 |

### 2.7 meeko（A-3 Track 1 ライブラリ）

| 項目 | 内容 |
|---|---|
| ライセンス | **LGPL-2.1** (Lesser GNU General Public License v2.1) |
| ソース | https://github.com/forlilab/Meeko |
| 利用 | AutoDock Vina 向けリガンド PDBQT 変換。`installExtra("meeko")` で Embedded Python に配備 |
| 再配布条件 | ライセンス表示。LGPL 成果物は EasyMolKit コードとは独立して配布 |

**LGPL-2.1 と EasyMolKit (MIT) の互換性**:
- Python の `import meeko` は **動的リンク** に相当する。LGPL-2.1 Section 6 により、  
  動的リンク利用者は LGPL ライブラリを改変して再コンパイルする権利があれば copyleft を回避できる
- EasyMolKit はバイナリを改変・静的リンクしない。`pip install` によりユーザーが meeko を取得し、  
  Python 実行時に動的にロードされるため、EasyMolKit (MIT) への copyleft 汚染は**発生しない**
- 参照: FSF「Can I use GPL-covered tools to build non-GPL works?」および LGPL Section 6

> **結論 (CL-7)**: `installExtra("meeko")` による動的リンク利用は EasyMolKit の MIT ライセンスと互換。

### 2.8 vina（A-3 Track 1 ライブラリ）

| 項目 | 内容 |
|---|---|
| ライセンス | Apache-2.0 |
| ソース | https://github.com/ccsb-scripps/AutoDock-Vina |
| 利用 | 分子ドッキングエンジン。`installExtra("vina")` で Embedded Python に配備 |
| 再配布条件 | ライセンス表示・変更表示 |

**Apache-2.0 と MIT の互換性**: Permissive ライセンス同士であり互換性あり。  
特許免除条項（Section 3）を含むが、再配布制限は課さない。

### 2.9 pdbfixer（A-3 Track 1 ライブラリ）

| 項目 | 内容 |
|---|---|
| ライセンス | MIT License |
| ソース | https://github.com/openmm/pdbfixer |
| 利用 | PDB 構造前処理（欠損残基補完・水素付加）。`installExtra("pdbfixer")` で自動インストール |
| 依存 | openmm >= 8.2 を自動インストール（§2.10 参照） |
| 再配布条件 | ライセンス表示 |

### 2.10 openmm（pdbfixer の自動依存）

| 項目 | 内容 |
|---|---|
| ライセンス | MIT + LGPL-2.1+ (コアは MIT、一部 LGPL コンポーネントを含む) |
| ソース | https://github.com/openmm/openmm |
| 利用 | pdbfixer の依存パッケージとして `pip install pdbfixer` 時に自動取得 |
| 再配布条件 | ライセンス表示。LGPL 部分は meeko と同様の動的リンク判断（§2.7 参照） |

> **結論**: `pdbfixer` + `openmm` の `import` による動的リンクは EasyMolKit (MIT) に影響しない。

---

## 3. ライセンス互換性マトリクス

| ライブラリ / データ | ライセンス | MIT 互換 | 再配布義務 |
|---|---|---|---|
| Python (CPython) | PSF | ✅ | ライセンス表示 |
| RDKit | BSD-3 | ✅ | ライセンス表示・名称使用制限 |
| ChEMBL データ | CC-BY-SA 3.0 | ⚠️ SA 条件あり | 帰属表示 + ShareAlike（データのみに適用）|
| PubChem データ | CC0 | ✅ | なし |
| PyMOL-OSS | PSF/BSD | ✅ | ライセンス表示 |
| pip | MIT | ✅ | ライセンス表示 |
| mordredcommunity | BSD-3-Clause | ✅ | ライセンス表示・名称使用制限 |
| Biopython | Biopython License (BSD-like) | ✅ | ライセンス表示・名称使用制限 |
| PubChemPy | MIT | ✅ | ライセンス表示 |
| PyTorch (CPU-only) | BSD-3-Clause | ✅ | ライセンス表示・名称使用制限 |
| PyTorch Geometric + scatter/sparse/cluster | MIT | ✅ | ライセンス表示 |
| HuggingFace Transformers | Apache-2.0 | ✅ | ライセンス表示・変更表示 |
| Open Babel | **GPLv2** | ⚠️ Track 2 のみ | Track 2（ユーザー導入）。EasyMolKit 本体に含めない |
| MDAnalysis | **GPLv2+** | ⚠️ Track 2 のみ | Track 2（ユーザー導入）。EasyMolKit 本体に含めない |
| meeko | **LGPL-2.1** | ✅ 動的リンク安全 | ライセンス表示。動的 import は copyleft 汚染なし（§2.7 CL-7）|
| vina | Apache-2.0 | ✅ | ライセンス表示・変更表示 |
| pdbfixer | MIT | ✅ | ライセンス表示 |
| openmm | MIT + LGPL-2.1+ | ✅ 動的リンク安全 | pdbfixer の自動依存。動的 import は安全（§2.10）|

**結論**: 全 Track 1 依存ライブラリが Permissive License（MIT / BSD-3 / Apache-2.0 / PSF）であり、EasyMolKit を MIT License で配布することに問題はない。Apache-2.0 は MIT と互換性がある（特許免除条項を含むが再配布制限なし）。

---

## 4. 再配布形態と義務

### 4.1 ソースコード配布（GitHub）

EasyMolKit のソースコードのみを配布する場合（`python_env/` は `.gitignore` で除外）:

- ✅ `LICENSE` (MIT) をルートに配置
- ✅ `THIRD_PARTY_NOTICES.md` に依存ライブラリのライセンスを記載
- ✅ `emk.setup.install()` がユーザー環境で Python + RDKit をダウンロード
- → Python / RDKit のバイナリ再配布には **該当しない**（ユーザーが自分でダウンロード）

### 4.2 バイナリ同梱配布（将来検討）

Embedded Python + RDKit を含むパッケージを配布する場合:

- ⚠️ Python の PSF License 全文を同梱
- ⚠️ RDKit の BSD-3 ライセンス全文を同梱
- ⚠️ `THIRD_PARTY_NOTICES.md` に全著作権表示を記載
- ⚠️ python.org / RDKit の配布条件を再確認

> **推奨**: 当面は **ソースコード配布のみ** とし、`emk.setup.install()` でユーザー環境に
> ダウンロードする方式を採用。バイナリ同梱配布は追加的な法的レビューが必要。

---

## 5. 誤解しやすいポイント

### 5.1 PyMOL の商用版と OSS 版

- **商用版**: Schrödinger 社が販売。有償ライセンス。EasyMolKit では使用しない
- **OSS 版**: GitHub で公開。ビルドが必要だが無償利用可能
- `pymol-open-source` パッケージを明示的に指定すること

### 5.2 RDKit のライセンスと学術利用

- RDKit は BSD-3 であり、学術・商用を問わず自由に利用可能
- ただし、RDKit の名前を使って EasyMolKit を宣伝することは禁止（BSD-3 第三条項）
- 「RDKit を利用しています」という事実の記述は問題ない

### 5.3 MATLAB ライセンス

- EasyMolKit は MATLAB のライセンスを同椅・再配布しない
- ユーザーが MATLAB ライセンスを保有していることが前提
- MATLAB Online は MathWorks のサービスとして利用

---

## 6. THIRD_PARTY_NOTICES.md テンプレート

> M-PUB フェーズで以下のテンプレートを使用して作成する。

```markdown
# Third-Party Notices

EasyMolKit uses the following third-party open-source software:

## Python (CPython)
- Copyright (c) 2001-2026 Python Software Foundation
- License: PSF License
- https://www.python.org/

## RDKit
- Copyright (c) 2006-2026 Greg Landrum and RDKit contributors
- License: BSD-3-Clause
- https://github.com/rdkit/rdkit

## ChEMBL Data (data/list/fda_drugs.csv, data/list/forensic_challenge.csv)
- Source: ChEMBL database, EMBL-EBI (https://www.ebi.ac.uk/chembl/)
- License: Creative Commons Attribution-ShareAlike 3.0 Unported (CC-BY-SA 3.0)
- These data files are distributed under CC-BY-SA 3.0, separate from the MIT-licensed code.

## PyMOL Open-Source (future)
- Copyright (c) Schrödinger, LLC
- License: Python-2.0
- https://github.com/schrodinger/pymol-open-source
```

---

## 7. チェックリスト

開発・リリース時に確認すること:

- [ ] `LICENSE` ファイルがルートに存在する
- [ ] `THIRD_PARTY_NOTICES.md` が全依存ライブラリを網羅している
- [ ] `emk.setup.install()` がダウンロード元として公式サイトのみを使用している
- [ ] README に Disclaimer / License セクションがある
- [ ] RDKit の名前を宣伝目的に使用していない
- [ ] PyMOL は OSS 版のみを対象としていることが明記されている
- [ ] Python のバイナリを Git リポジトリに含めていない（`.gitignore` 確認）

---

## 8. 未解決事項

| # | 項目 | 説明 | 状態 |
|---|---|---|---|
| CL-1 | EasyMolKit 本体のライセンス最終決定 | MIT 予定だが正式決定は公開前 | ⬜ 要決定 |
| CL-2 | 大学・研究機関での配布制約 | 所属機関の知財ポリシーとの整合確認 | ⬜ 要確認 |
| CL-3 | PyMOL-OSS のビルド配布 | ビルド済みバイナリの再配布可否 | ⬜ 要調査 |
| CL-4 | 法的レビュー | バイナリ同梱配布時の法的確認 | ⬜ 要対応（将来）|
| CL-5 | mordredcommunity ライセンス確認 | **✅ 確認済み**: BSD-3-Clause（JacksonBurns/mordred-community）。fork 時にライセンス所有者名を更新。MIT 互換あり | ✅ 確認済み |
| CL-6 | Biopython ライセンス確認 | **✅ 確認済み**: 「Biopython License Agreement」（BSD-like Permissive）+ 一部ファイルは BSD-3 との二重ライセンス。MIT 互換あり。名称使用制限あり（BSD-3 同等） | ✅ 確認済み |
| CL-7 | PyTorch ライセンス確認 | **✅ 確認済み**: BSD-3-Clause（Meta Platforms, Inc.）。`pip install torch` で CPU-only wheel を取得。MIT 互換あり。名称使用制限あり（BSD-3 同等）。参照: https://github.com/pytorch/pytorch/blob/main/LICENSE | ✅ 確認済み |
| CL-8 | PyTorch Geometric ライセンス確認 | **✅ 確認済み**: MIT（Matthias Fey）。companion パッケージ（torch_scatter/sparse/cluster）も MIT。MIT 互換あり。参照: https://github.com/pyg-team/pytorch_geometric/blob/master/LICENSE | ✅ 確認済み |
| CL-9 | HuggingFace Transformers ライセンス確認 | **✅ 確認済み**: Apache-2.0（Hugging Face, Inc.）。MIT との互換性あり（特許免除条項を含むが再配布制限なし）。変更ファイルへの変更表示義務あり。参照: https://github.com/huggingface/transformers/blob/main/LICENSE | ✅ 確認済み |
