# Function Catalog — EasyMolKit

> Compact reference for all 51 public functions across all modules.
> Full signatures and error IDs → [function_reference.md](function_reference.md)
> Japanese version → [function_catalog.ja.md](function_catalog.ja.md)

---

## emk.setup (9 functions)

| Function | Summary |
|---|---|
| `install(PythonVersion="3.10")` | Deploy Embedded Python + RDKit to `python_env/` (Desktop only) |
| `installOnline(Config=struct())` | Bootstrap pip and install RDKit on MATLAB Online |
| `initPython()` | Detect platform and configure `pyenv` (idempotent) |
| `verify()` | Returns struct with `.python`, `.rdkit`, `.version`; never throws |
| `installExtra(name)` | pip-install a Track 1 add-on into Embedded Python |
| `useExternal(pythonPath)` | Switch `pyenv` to an external CPython (Track 2; must be called before Python loads) |
| `validate(Libraries=string.empty)` | Returns install-status table for all or specified libraries; never throws |
| `recipe(name)` | Print installation recipe and license info to Command Window |
| `installTrack2(name, BasePython="")` | Create venv for a Track 2 library and write `settings.json` |

---

## emk.mol (8 functions)

| Function | Summary |
|---|---|
| `fromSmiles(smiles)` | Parse SMILES string → `py.rdkit.Chem.Mol` |
| `toSmiles(mol)` | Canonical SMILES string from a Mol object |
| `isValid(smiles)` | Returns `logical`; never throws (except `rdkitError`) |
| `hasSubstruct(mol, query)` | SMARTS/Mol substructure match; cell input → `logical(1,N)` |
| `toStruct(mol, Format="molblock")` | Mol → MATLAB struct (molblock or pickle format) |
| `fromStruct(s)` | Reverse of `toStruct`; auto-detects molblock or pickle |
| `toTable(mols, Properties=…)` | Cell of Mols → table (SMILES + selected descriptors) |
| `scaffold(mol)` | Returns Bemis-Murcko scaffold as Mol |

---

## emk.descriptor (6 functions)

| Function | Summary |
|---|---|
| `molWeight(mol)` | Average molecular weight [g/mol] |
| `calculate(mol, descriptorNames=<all10>)` | Compute up to 11 standard descriptors → struct |
| `batchCalculate(mols, descriptorNames=<all10>)` | Multiple Mols → table (NaN for invalid rows) |
| `mordred(mol, descriptorNames=[])` | ~1800 Mordred 2D descriptors → struct; requires `installExtra("mordred")` |
| `mordredBatch(mols, descriptorNames=[])` | Batch Mordred computation → table; requires `installExtra("mordred")` |
| `mordredNames()` | Returns sorted `string(1×N)` list of ~1800 available Mordred descriptor names |

**Standard descriptors** (`all10`): `MolWt`, `ExactMolWt`, `LogP`, `TPSA`, `NumHAcceptors`, `NumHDonors`,
`NumRotatableBonds`, `RingCount`, `FractionCSP3`, `HeavyAtomCount`. `MolFormula` available on explicit request.

---

## emk.fingerprint (3 functions)

| Function | Summary |
|---|---|
| `morgan(mol, Radius=2, NBits=2048)` | Morgan (ECFP4) fingerprint → `py.rdkit.DataStructs.ExplicitBitVect` |
| `maccs(mol)` | 167-bit MACCS keys fingerprint |
| `toArray(fp)` | Convert fingerprint → `logical(1,N)` MATLAB array |

---

## emk.similarity (4 functions)

| Function | Summary |
|---|---|
| `tanimoto(fp1, fp2)` | Tanimoto coefficient ∈ [0, 1] |
| `dice(fp1, fp2)` | Dice coefficient ∈ [0, 1] |
| `rankBy(queryFp, dbFps, N=Inf, Metric="tanimoto")` | Rank database FPs by similarity; returns `.Indices`, `.Scores`, `.Metric` |
| `matrix(fps, Metric="tanimoto")` | N×N symmetric pairwise similarity matrix |

---

## emk.filter (1 function)

| Function | Summary |
|---|---|
| `lipinski(tbl, MaxViolations=0)` | Add `Violations_Ro5` (double) and `Pass_Ro5` (logical) columns to an input table |

Requires columns: `MolWt`, `LogP`, `NumHDonors`, `NumHAcceptors`.

---

## emk.io (3 functions)

| Function | Summary |
|---|---|
| `readSdf(filePath)` | SDF file → `1×N cell` of Mol objects (failed molecules skipped with warning) |
| `writeSdf(mols, filePath)` | Write Mol cell array to SDF file |
| `readSmilesList(filePath)` | One-SMILES-per-line file → `1×N cell` of Mol objects |

---

## emk.viz (1 function)

| Function | Summary |
|---|---|
| `draw2d(mol, Title="", Width=300, Height=300)` | Render 2D structure to a MATLAB figure window |

---

## emk.db (5 functions)

| Function | Summary |
|---|---|
| `searchPubchem(query, Type="name")` | PubChem search → 5-column table; no Python required |
| `pubchemFetch(identifier, NameSpace="name", MaxSynonyms=10)` | Extended PubChem fetch via PubChemPy → struct |
| `searchChembl(query, Type="name")` | ChEMBL search → 8-column table; no Python required |
| `searchChemblTarget(query, TargetType="SINGLE PROTEIN", MaxRows=10)` | ChEMBL target search → 4-column table |
| `getChemblActivity(targetId, ActivityType="IC50", MaxRows=50, MinActivity_nM=Inf)` | ChEMBL bioactivity data → 5-column table |

`pubchemFetch` requires `installExtra("pubchempy")`.

---

## emk.util (2 functions)

| Function | Summary |
|---|---|
| `isOnline()` | Returns `logical` — `true` when running on MATLAB Online |
| `benchmarkBatch(smilesList, descriptorNames=<all10>)` | Throughput benchmark → struct with timing fields |

---

## src/config (1 function)

| Function | Summary |
|---|---|
| `emkLoadConfig()` | Load settings (env vars > `settings.json` > defaults) → struct |

---

## src/util (8 functions)

| Function | Summary |
|---|---|
| `logInfo(msg, ...)` | `[HH:MM:SS][INFO]  message` to console |
| `logWarn(msg, ...)` | `[HH:MM:SS][WARN]  message` to console |
| `logError(msg, ...)` | `[HH:MM:SS][ERROR]  message` to console |
| `logDebug(msg, ...)` | Output only when `EMK_LOG_VERBOSE=1` |
| `logProgress(i, n, label)` | Progress bar for loop iterations |
| `logSection(scriptId, label, layer)` | Section banner for tutorial scripts |
| `makeRunDir(Prefix="", BaseDir="result/runs")` | Create timestamped run directory; never use `mkdir` directly |
| `resolveProjectRoot()` | 3-stage fallback to locate and `cd` to project root |
