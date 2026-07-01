# RP01: ESOL Classical QSAR + L05 Extended Descriptors

[English README](README.en.md)

> **逶ｮ逧・*: Delaney (2004) 縺ｮ ESOL 邱壼ｽ｢ QSAR 繧・EasyMolKit 荳翫〒蜀咲樟縺励・> `L05` 險倩ｿｰ蟄先僑蠑ｵ (`TPSA`, `HBD`, `HBA`, `FractionCSP3`, `QED`) 縺・> `logS` 莠域ｸｬ繧呈隼蝟・☆繧九°繧定ｩ穂ｾ｡縺吶ｋ縲・
---

## Overview

| 鬆・岼 | 蜀・ｮｹ |
|---|---|
| 隲匁枚 | Delaney, J.S. (2004). *ESOL: Estimating Aqueous Solubility Directly from Molecular Structure.* *J. Chem. Inf. Comput. Sci.* 44(3):1000-1005 |
| DOI | [10.1021/ci034243x](https://doi.org/10.1021/ci034243x) |
| 繧ｿ繧ｹ繧ｯ | 豌ｴ貅ｶ隗｣蠎ｦ莠域ｸｬ・亥屓蟶ｰ, `logS`・・|
| Model A | 4 險倩ｿｰ蟄千ｷ壼ｽ｢蝗槫ｸｰ: `LogP + MolWt + RotBonds + AroProp` |
| Model B | 9 險倩ｿｰ蟄千ｷ壼ｽ｢蝗槫ｸｰ: Model A + `TPSA + HBD + HBA + FractionCSP3 + QED` |
| 繝・・繧ｿ | ESOL / MoleculeNet・・,128 蛻・ｭ撰ｼ・|
| 隲匁枚蜈ｬ陦ｨ蛟､ | Delaney Table 2: training RMSE = 0.996 |

---

## Environment

螳溯｡梧凾縺ｮ螳溘ヰ繝ｼ繧ｸ繝ｧ繝ｳ縺ｯ `result/runs/<timestamp>/lock_snapshot.json` 縺ｫ險倬鹸縺輔ｌ繧九・
| 鬆・岼 | 隕∽ｻｶ |
|---|---|
| MATLAB | R2025a 莉･髯搾ｼ域､懆ｨｼ螳溯｡後・ R2026a・・|
| Python | 3.10・・asyMolKit Embedded Python・・|
| RDKit | 2022.03 莉･髯搾ｼ域､懆ｨｼ螳溯｡後・ 2024.03.6・・|
| Toolbox | Statistics and Machine Learning Toolbox |

### Reproducibility Controls

- Section 0 縺ｧ MATLAB RNG 繧・`rng(42, "twister")` 縺ｫ蛻晄悄蛹悶☆繧九・- Section 5 縺ｧ `cvpartition` 縺ｮ逶ｴ蜑阪↓蜀榊ｺｦ `rng(42, "twister")` 繧帝←逕ｨ縺励∝燕谿ｵ縺ｧ RNG 迥ｶ諷九′豸郁ｲｻ縺輔ｌ縺ｦ繧・5-fold 蛻・牡縺悟・迴ｾ縺輔ｌ繧九ｈ縺・↓縺励※縺・ｋ縲・
### Descriptor Definitions

| 險倩ｿｰ蟄・| 繝・・繝ｫ | 繝舌・繧ｸ繝ｧ繝ｳ | 螳夂ｾｩ |
|---|---|---|---|
| LogP | RDKit | 螳溯｡梧凾險倬鹸 | Wildman-Crippen `MolLogP` |
| MolWt | RDKit | 螳溯｡梧凾險倬鹸 | `Descriptors.MolWt` |
| NumRotatableBonds | RDKit | 螳溯｡梧凾險倬鹸 | `CalcNumRotatableBonds` strict SMARTS |
| HeavyAtomCount | RDKit | 螳溯｡梧凾險倬鹸 | 髱樊ｰｴ邏蜴溷ｭ先焚 |
| AromaticProportion | RDKit | 螳溯｡梧凾險倬鹸 | `NumAromaticAtoms / HeavyAtomCount` |
| TPSA | RDKit | 螳溯｡梧凾險倬鹸 | `CalcTPSA` |
| NumHDonors | RDKit | 螳溯｡梧凾險倬鹸 | `CalcNumHBD` |
| NumHAcceptors | RDKit | 螳溯｡梧凾險倬鹸 | `CalcNumHBA` |
| FractionCSP3 | RDKit | 螳溯｡梧凾險倬鹸 | `CalcFractionCSP3` |
| QED | RDKit | 螳溯｡梧凾險倬鹸 | `rdkit.Chem.QED.qed` |

> **豕ｨ險・*: `AromaticProportion` 縺ｯ RDKit `GetIsAromatic()` 縺九ｉ險育ｮ励☆繧九・> Delaney 隲匁枚縺ｮ `clogP` 縺ｨ譛ｬ RP 縺ｮ RDKit `MolLogP` 縺ｯ蜷御ｸ縺ｧ縺ｯ縺ｪ縺・・
---

## Data

- **Source**: Delaney (2004) ESOL dataset, hosted by DeepChem / MoleculeNet
- **URL**: `https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/delaney-processed.csv`
- **License**: 譏守､ｺ逧・↑繝ｩ繧､繧ｻ繝ｳ繧ｹ陦ｨ險倥・縺ｪ縺・′縲∝・髢九・繝ｳ繝√・繝ｼ繧ｯ縺ｨ縺励※蠎・￥蛻ｩ逕ｨ縺輔ｌ縺ｦ縺・ｋ
- **Cache path**: `data/benchmark/esol.csv`
- **Count**: 1,128 蛻・ｭ撰ｼ郁ｫ匁枚縺ｮ 1,144 繧医ｊ 16 蟆代↑縺・ｼ・- **Data hash**: `result/runs/<ts>/lock_snapshot.json` 縺ｮ `dataset_sha256`

---

## Script

```text
repro/rp01_esol/rp01_esol.m
```

**螳溯｡梧婿豕・*: 繝励Ο繧ｸ繧ｧ繧ｯ繝医Ν繝ｼ繝医ｒ繧ｫ繝ｬ繝ｳ繝医ョ繧｣繝ｬ繧ｯ繝医Μ縺ｫ縺励※ MATLAB 縺ｧ section 縺斐→縺ｫ螳溯｡後☆繧九・
| Section | 蜀・ｮｹ |
|---|---|
| Section 0 | 繧ｻ繝・ヨ繧｢繝・・縲∫腸蠅・叙蠕励～lock_template.json` 縺九ｉ縺ｮ RF03 蝓ｺ貅冶ｪｭ霎ｼ縲～rng(42, "twister")` 縺ｫ繧医ｋ RNG 蛻晄悄蛹・|
| Section 1 | ESOL 繝・・繧ｿ隱ｭ霎ｼ縺ｨ dataset SHA-256 險育ｮ・|
| Section 2 | SMILES 隗｣譫舌→險倩ｿｰ蟄占ｨ育ｮ・|
| Section 3 | Model A / Model B 縺ｮ蟄ｦ鄙偵∝・繝・・繧ｿ fit 謖・ｨ吶・邂怜・縲∽ｿよ焚謚ｽ蜃ｺ縲ゝPSA 隨ｦ蜿ｷ險ｺ譁ｭ |
| Section 4 | VIF 縺ｫ繧医ｋ螟夐㍾蜈ｱ邱壽ｧ遒ｺ隱・|
| Section 5 | `rng(42, "twister")` 縺ｧ蜀阪す繝ｼ繝峨＠縺ｦ 5-fold CV 螳溯｡・|
| Section 6 | Model A vs Model B 縺ｮ paired t-test |
| Section 7 | RF03 讀懆ｨｼ |
| Section 8 | 蜃ｺ蜉帑ｿ晏ｭ假ｼ・metrics.json`, `predictions.csv`, `lock_snapshot.json`・・|

### Section 3 diagnostic note

Section 3 縺ｫ縺ｯ `TPSA` 縺ｮ隨ｦ蜿ｷ險ｺ譁ｭ (`B2`) 繧貞性繧縲ゅ％繧後・ `corr(TPSA, logS)` 縺ｮ蜻ｨ霎ｺ逶ｸ髢｢縺ｮ隨ｦ蜿ｷ縺ｨ縲｀odel B 縺ｫ縺翫￠繧・`TPSA` 縺ｮ蛛丞屓蟶ｰ菫よ焚縺ｮ隨ｦ蜿ｷ繧呈ｯ碑ｼ・☆繧九ｂ縺ｮ縺ｧ縲∫ｬｦ蜿ｷ縺碁・ｻ｢縺励◆蝣ｴ蜷医・ `TPSA` 縺ｨ `LogP`縲～HBD`縲～HBA` 縺ｮ逶ｸ髢｢繧定ｿｽ蜉縺ｧ繝ｭ繧ｰ蜃ｺ蜉帙＠縲∝､夐㍾蜈ｱ邱壽ｧ縺ｫ蝓ｺ縺･縺剰ｧ｣驥医ｒ陬懷勧縺吶ｋ縲ょ酔縺倩ｦ∫ｴ・・ `metrics.json` 縺ｮ `tpsa_b2` 縺ｫ菫晏ｭ倥＆繧後ｋ縲・
---

## Result

蛻晏屓螳溯｡梧律縺ｯ 2026-06-19縲∵怙譁ｰ縺ｮ讀懆ｨｼ螳溯｡梧律縺ｯ 2026-06-30縲・讀懆ｨｼ貂医∩ artifact 縺ｯ `result/runs/20260630_031306_rp01_esol/`縲・
### Current validated status

| 謖・ｨ・| 蛟､ | 繧ｹ繝・・繧ｿ繧ｹ |
|---|---|---|
| Model A RMSE・亥・繝・・繧ｿ fit・・| 1.0094 | Delaney training RMSE 縺ｨ縺ｮ豈碑ｼ・畑 |
| Model B RMSE・亥・繝・・繧ｿ fit・・| 0.9655 | 蜿り・､ |
| Model A RMSE CV | 1.0166 +/- 0.024 | PASS (`<= 1.20`) |
| Model A R^2 CV | 0.7638 +/- 0.022 | PASS (`>= 0.75`) |
| Model B RMSE CV | 0.9798 +/- 0.027 | PASS (`<= 1.20`) |
| Model B R^2 CV | 0.7804 +/- 0.023 | PASS (`>= 0.75`) |
| L05 delta RMSE | -0.0368 | 蜿り・ュ蝣ｱ |
| L05 delta R^2 | +0.0166 | 蜿り・ュ蝣ｱ |

### 隲匁枚縺ｨ縺ｮ逶ｴ謗･豈碑ｼ・
| 謖・ｨ・| 隲匁枚 | 譛ｬ RP | 豕ｨ險・|
|---|---|---|---|
| Training RMSE | 0.996 | Model A: 1.0094 | 4 險倩ｿｰ蟄先ｧ区・縺ｯ蜷後§縺縺後∬ｨ倩ｿｰ蟄仙ｮ溯｣・→繝・・繧ｿ迚医′逡ｰ縺ｪ繧・|
| Training RMSE | L05 蛛ｴ縺ｯ隲匁枚險倩ｼ峨↑縺・| Model B: 0.9655 | 諡｡蠑ｵ繝｢繝・Ν縺ｯ譛ｬ RP 迢ｬ閾ｪ |

### Environment from validated run

| 鬆・岼 | 蛟､ |
|---|---|
| MATLAB | R2026a |
| Python | 3.10 |
| RDKit | 2024.03.6 |
| Commit | `a4e3305` |

---

## Verification

RF03 縺ｧ縺ｯ Cat A 縺・PASS 蛻､螳壽擅莉ｶ縺ｧ縺ゅｊ縲｀odel A vs Model B 縺ｮ豈碑ｼ・・ Cat B 縺ｮ蜿り・ュ蝣ｱ縺ｨ縺励※謇ｱ縺・・
### Cat A: Absolute thresholds

| 謖・ｨ・| 蝓ｺ貅・| 譬ｹ諡 |
|---|---|---|
| RMSE CV | `<= 1.20` | 隲匁枚縺ｮ training RMSE 0.996 縺ｫ蟇ｾ縺励∝ｮ溯｣・ｷｮ縺ｨ繝・・繧ｿ蟾ｮ縺ｮ險ｱ螳ｹ蟷・ｒ蜉縺医◆ |
| R^2 CV | `>= 0.75` | 蜈郁｡後☆繧・RP00 繝代う繝ｭ繝・ヨ讀懆ｨｼ縺ｫ蝓ｺ縺･縺丞崋螳壻ｸ矩剞 |

### Cat B: Relative comparison

| 謖・ｨ・| Model A | Model B | 蟾ｮ蛻・| t(4) | p(one-sided) |
|---|---|---|---|---|---|
| CV RMSE | 1.0166 | 0.9798 | -0.0368 | 6.591 | 0.001 |
| CV R^2 | 0.7638 | 0.7804 | +0.0166 | -7.031 | 0.001 |

**隗｣驥・*: L05 諡｡蠑ｵ縺ｯ CV RMSE 繧剃ｽ惹ｸ九＆縺帙，V R^2 繧剃ｸ頑・縺輔○縺溘ゅ←縺｡繧峨ｂ fold 蜊倅ｽ阪・ paired t-test 縺ｧ縺ｯ譛画э縺ｧ縺ゅｋ縲・
### Cat B2: TPSA sign verification

縺薙・險ｺ譁ｭ縺ｯ RF03 縺ｮ PASS/FAIL 譚｡莉ｶ縺ｧ縺ｯ縺ｪ縺・′縲｀odel B 縺ｮ隗｣驥郁ｨ倬鹸縺ｨ縺励※菫晄戟縺吶ｋ縲ょ､壼､蛾㍼蝗槫ｸｰ荳ｭ縺ｮ `TPSA` 蛛丞屓蟶ｰ菫よ焚縺ｮ隨ｦ蜿ｷ縺後～TPSA` 縺ｨ `logS` 縺ｮ蜻ｨ霎ｺ逶ｸ髢｢縺ｮ隨ｦ蜿ｷ繧堤ｶｭ謖√＠縺ｦ縺・ｋ縺九ｒ遒ｺ隱阪＠縲・・ｻ｢縺励◆蝣ｴ蜷医・ `LogP`縲～HBD`縲～HBA` 縺ｨ縺ｮ逶ｸ髢｢繧剃ｽｵ險倥☆繧九ゆｸｻ縺溘ｋ諢丞峙縺ｯ縲∝喧蟄ｦ逧・↓蜊倡ｴ斐↑蝗譫懆ｧ｣驥医〒縺ｯ縺ｪ縺上∝､夐㍾蜈ｱ邱壽ｧ縺ｮ蠖ｱ髻ｿ繧貞・繧雁・縺代ｋ縺薙→縺ｫ縺ゅｋ縲ょ盾辣ｧ蜈医・ `metrics.json` 縺ｮ `tpsa_b2` 縺ｨ `model_b_coefficients.csv`縲・
### Tolerance rationale

1. Delaney 縺ｯ `clogP` 繧剃ｽｿ逕ｨ縺励∵悽 RP 縺ｯ RDKit `MolLogP` 繧剃ｽｿ逕ｨ縺吶ｋ縲・2. MoleculeNet 迚・ESOL 縺ｯ隲匁枚縺ｮ 1,144 蛻・ｭ舌〒縺ｯ縺ｪ縺・1,128 蛻・ｭ舌〒縺ゅｋ縲・3. Rotatable bond 縺ｮ螳夂ｾｩ蟾ｮ縺御ｿよ焚繧・ｪｬ譏守紫縺ｫ蠖ｱ髻ｿ縺励≧繧九・
---

## Discussion

### Differences from the paper

| 蟾ｮ蛻・| 蜀・ｮｹ |
|---|---|
| 隧穂ｾ｡險ｭ螳・| 隲匁枚縺ｯ training RMSE 繧貞ｱ蜻翫＠縲∵悽 RP 縺ｯ 5-fold CV 繧剃ｸｻ謖・ｨ吶→縺吶ｋ縲ゅ◆縺縺礼峩謗･豈碑ｼ・・縺溘ａ縺ｫ蜈ｨ繝・・繧ｿ fit RMSE 繧ゆｽｵ險倥☆繧・|
| LogP 螳溯｣・| 隲匁枚縺ｯ `clogP`縲∵悽 RP 縺ｯ RDKit `MolLogP` |
| 繝・・繧ｿ莉ｶ謨ｰ | 隲匁枚縺ｯ 1,144縲∵悽 RP 縺ｯ MoleculeNet 迚・1,128 |
| 險ｺ譁ｭ遽・峇 | 譛ｬ RP 縺ｧ縺ｯ蜴溯送縺ｫ縺ｪ縺・VIF 縺ｨ TPSA 隨ｦ蜿ｷ險ｺ譁ｭ繧定ｿｽ蜉縺励※縺・ｋ |

### Main takeaways

- Delaney 縺ｮ 4 險倩ｿｰ蟄舌Δ繝・Ν縺ｯ EasyMolKit 荳翫〒蜀咲樟蜿ｯ閭ｽ縺ｧ縺ゅｋ縲・- L05 諡｡蠑ｵ縺ｯ縺薙・螳溯｣・→繝・・繧ｿ譚｡莉ｶ縺ｧ縺ｯ莠域ｸｬ諤ｧ閭ｽ繧呈隼蝟・☆繧九・- Model B 縺ｧ縺ｯ `TPSA` 蜻ｨ霎ｺ縺ｮ螟夐㍾蜈ｱ邱壽ｧ縺悟ｼｷ縺上∽ｿよ焚隨ｦ蜿ｷ縺ｮ隗｣驥医ｒ `tpsa_b2` 險ｺ譁ｭ縺九ｉ蛻・ｊ髮｢縺励※縺ｯ縺ｪ繧峨↑縺・・
---

## Related Files

| 繝輔ぃ繧､繝ｫ | 蜀・ｮｹ |
|---|---|
| `README.md` | 闍ｱ隱樒沿 README |
| `README.jp.md` | 譌･譛ｬ隱樒沿 README |
| `rp01_esol.m` | 蜀咲樟繧ｹ繧ｯ繝ｪ繝励ヨ譛ｬ菴・|
| `lock_template.json` | RF02 繝ｭ繝・け繧ｹ繧ｭ繝ｼ繝槫・ RF03 蝓ｺ貅悶た繝ｼ繧ｹ |
| `result/runs/20260630_031306_rp01_esol/metrics.json` | 讀懆ｨｼ貂医∩縺ｮ隧穂ｾ｡謖・ｨ吩ｸ蠑上Ａtpsa_b2` 繧貞性繧 |
| `result/runs/<ts>/metrics.json` | 縺昴・莉悶・螳溯｡後↓蟇ｾ縺吶ｋ隧穂ｾ｡謖・ｨ・|
| `result/runs/<ts>/predictions.csv` | 螳滓ｸｬ蛟､縲∽ｺ域ｸｬ蛟､縲∬ｨ倩ｿｰ蟄・|
| `result/runs/<ts>/model_b_coefficients.csv` | Model B 菫よ焚陦ｨ |
| `result/runs/<ts>/predicted_vs_actual.png` | 螳滓ｸｬ vs 莠域ｸｬ謨｣蟶・峙 |
| `repro/rp00_esol/` | 蜈郁｡後ヱ繧､繝ｭ繝・ヨ RP |

