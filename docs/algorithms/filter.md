# Algorithm Guide — emk.filter

> アルゴリズム根拠インデックス → [algorithm_guide.md](../algorithm_guide.md)  
> API シグネチャ → [function_reference.md](../function_reference.md)

---

## 10.1 `lipinski`

**設計意図**: Lipinski の「Rule of Five」(Ro5) を MATLAB table に適用し、経口吸収性の初期スクリーニングを行う。医薬品候補化合物の ADMET 評価の第一段階として使用する。

**アルゴリズム概要**:
1. `istable(tbl)` を確認（RDKit 不要）
2. `MaxViolations` が整数 ∈ [0, 4] であることを確認（RDKit 不要）
3. 必須カラムの存在確認: `MolWt`, `LogP`, `NumHDonors`, `NumHAcceptors`（RDKit 不要）
4. 各行に対して Ro5 違反数を計算（NaN の扱いに注意）
5. `Violations <= MaxViolations` で Pass/Fail を決定
6. 元の table に `Violations_Ro5` (double) と `Pass_Ro5` (logical) カラムを追加して返す

**Ro5 閾値と根拠**:

| 記述子 | 閾値 | 違反条件 |
|---|---|---|
| MolWeight (MW) | 500 | MW > 500 |
| LogP | 5 | LogP > 5 |
| HBD (H-bond donors) | 5 | HBD > 5 |
| HBA (H-bond acceptors) | 10 | HBA > 10 |

最大違反数 `MaxViolations = 0` がデフォルト（厳格モード）。`MaxViolations = 1` は「Ro5 の 1 違反許容」のソフトモードに対応（ペプチド模倣化合物の評価等で使用）。

**NaN の扱いに関する重要な注意事項**:
MATLAB では `NaN > 500` は `false` を返す（MATLAB の比較演算子の仕様）。
これは NaN を「違反なし」として扱うことを意味する。`lipinski` 関数はこの動作を許容し、
NaN を含む行が `Pass_Ro5 = true` になる可能性についてユーザーは `logWarn` で警告を受ける。

**バリデーション順序の根拠**:
「RDKit 不要で発火するバリデーション → RDKit 必要なもの」の原則に従い、
`istable` → `MaxViolations` → カラム存在確認 → 実際の計算 の順で検証する。
これにより RDKit 未起動セッションでも入力エラーが確実に捕捉される。

**引用文献**:
- Lipinski, C.A., Lombardo, F., Dominy, B.W. & Feeney, P.J. (1997). Experimental and computational approaches to estimate solubility and permeability in drug discovery and development settings. *Adv. Drug Deliv. Rev.* 23(1-3):3–25. DOI: 10.1016/S0169-409X(96)00423-1
- Veber, D.F. et al. (2002). Molecular properties that influence the oral bioavailability of drug candidates. *J. Med. Chem.* 45(12):2615–2623. DOI: 10.1021/jm020017n

**テスト戦略** (`tests/unit/TestFilter.m`):

テストカテゴリ:

| カテゴリ | 内容 |
|---|---|
| TC-1: 入力バリデーション | `istable`・`MaxViolations`・カラム存在確認 |
| TC-2: 出力スキーマ | `Violations_Ro5` (double) + `Pass_Ro5` (logical) カラム追加 |
| TC-3: 全合格ケース | 全 Ro5 準拠 → 全て `Pass_Ro5 = true` |
| TC-4: 全不合格ケース | 全 Ro5 違反 → 全て `Pass_Ro5 = false` |
| TC-5: 混在ケース | 一部合格・一部不合格 |
| TC-6: 境界値 | MW=500, LogP=5, HBD=5, HBA=10（境界値は違反ではない） |
| TC-7: MaxViolations | MaxViolations=1 で 1 違反を許容 |
| TC-8: NaN 境界 | NaN 記述子が `Pass_Ro5 = true` になる（false negative 警告） |
| TC-9: RDKit 統合 | `toTable` 出力を直接 `lipinski` に渡す統合テスト |

具体的なテスト一覧:

| テスト名 | 検証内容 | RDKit 要否 |
|---|---|---|
| `test_lipinski_nonTableInput_throwsInvalidInput` | array 入力 → `invalidInput` | 不要 |
| `test_lipinski_invalidMaxViolations_throwsInvalidMaxViol` | `MaxViolations = 5` → `invalidMaxViol` | 不要 |
| `test_lipinski_negativeMaxViolations_throwsInvalidMaxViol` | `MaxViolations = -1` → `invalidMaxViol` | 不要 |
| `test_lipinski_nonIntMaxViolations_throwsInvalidMaxViol` | `MaxViolations = 1.5` → `invalidMaxViol` | 不要 |
| `test_lipinski_missingMolWeightCol_throwsMissingColumns` | MW 列欠如 → `missingColumns` | 不要 |
| `test_lipinski_missingLogPCol_throwsMissingColumns` | LogP 列欠如 → `missingColumns` | 不要 |
| `test_lipinski_missingHBDCol_throwsMissingColumns` | HBD 列欠如 → `missingColumns` | 不要 |
| `test_lipinski_missingHBACol_throwsMissingColumns` | HBA 列欠如 → `missingColumns` | 不要 |
| `test_lipinski_returnsTable` | 戻り値が `table` | 不要 |
| `test_lipinski_outputHasViolationsRo5Column` | `Violations_Ro5` カラム存在 | 不要 |
| `test_lipinski_outputHasPass_Ro5Column` | `Pass_Ro5` カラム存在 | 不要 |
| `test_violationsRo5Column_isDouble` | `Violations_Ro5` が double ★ | 不要 |
| `test_lipinski_ro5PassIsLogical` | `Pass_Ro5` が logical ★ | 不要 |
| `test_lipinski_allPass_allTrue` | 全 Ro5 準拠 → `Pass_Ro5 = [true, true, true]` | 不要 |
| `test_lipinski_allFail_allFalse` | 全 Ro5 違反 → `Pass_Ro5 = [false, false, false]` | 不要 |
| `test_lipinski_mixedResults_correctPassFail` | 混在 → 正確な Pass/Fail ベクトル | 不要 |
| `test_lipinski_boundary_MW500_notViolation` | MW=500 → 違反なし（> ではなく >= でない境界確認）★ | 不要 |
| `test_lipinski_boundary_LogP5_notViolation` | LogP=5 → 違反なし | 不要 |
| `test_lipinski_boundary_HBD5_notViolation` | HBD=5 → 違反なし | 不要 |
| `test_lipinski_boundary_HBA10_notViolation` | HBA=10 → 違反なし | 不要 |
| `test_lipinski_maxViol1_allows1Violation` | 1 違反分子が `MaxViolations=1` で Pass | 不要 |
| `test_lipinski_maxViol1_rejects2Violations` | 2 違反分子が `MaxViolations=1` で Fail | 不要 |
| `test_lipinski_nanMW_doesNotThrow` | NaN 記述子が例外なし ★★ | 不要 |
| `test_lipinski_nanMW_treatedAsNoViolation` | NaN MW → `Violations_Ro5 = 0`（false negative 警告）★★ | 不要 |
| `test_lipinski_rdkit_toTable_integration` | `toTable` 出力 → `lipinski` → `Pass_Ro5` が全 true | 要 |

**★ 型検証テスト (TC-2) の重要性**:
`test_violationsRo5Column_isDouble` は `double` 型を明示検証する。`int32` や `uint8` が
返った場合（算術演算でのデフォルト変換等）を検出する。

**★★ NaN テスト (TC-8) の重要性**:
`test_lipinski_nanMW_treatedAsNoViolation` は「NaN が false negative になる既知の動作」を
テストで文書化する役割を持つ。ユーザーへの `logWarn` 出力を確認することも含める。
