# Function Reference — EasyMolKit

> Complete list of function signatures and options. Primary reference during development.
> Algorithm rationale → [algorithm_guide.md](algorithm_guide.md)

```
src/+emk/+setup/       src/+emk/+mol/        src/+emk/+descriptor/
src/+emk/+fingerprint/ src/+emk/+similarity/  src/+emk/+filter/
src/+emk/+scaffold/    src/+emk/+dataset/     src/+emk/+cluster/
src/+emk/+diversity/   src/+emk/+conformer/   src/+emk/+shape/
src/+emk/+io/          src/+emk/+viz/         src/+emk/+util/
src/+emk/+repro/       src/config/            src/util/
```

---

## emk.setup

### `emk.setup.install(PythonVersion="3.10")`
Automatically deploys Embedded Python + RDKit into `python_env/` for Desktop use.
Runs path-length check → zip extraction → pip install → verify in one call.  
Errors: `notDesktop` (Online environment), `pathTooLong` (>240 chars), `downloadFailed`, `installFailed`.

### `emk.setup.installOnline(Config=struct())`
For MATLAB Online. Runs `get-pip.py` → `pip install rdkit-pypi==<ver> --user` → `py.sys.path` insertion → verify.
Skips if an already-matching version is detected.  
When `Config` is passed, libraries with `cfg.optionalLibraries.<name>=true` are additionally installed via
`pip install --user` (symmetric API to the Desktop `install(Config=cfg)`).  
Errors: `notOnline` (Desktop environment), `bootstrapFailed`, `installFailed`.

### `emk.setup.initPython()`
Detects the platform → calls `pyenv(Version=..., ExecutionMode="OutOfProcess")`.
Idempotent: skips silently if Python is already Loaded.  
Errors: `notInstalled` (Desktop + `python_env/python.exe` missing), `pyenvFailed`.

### `emk.setup.verify()`
Returns: `struct` — `.python` (logical), `.rdkit` (logical), `.version` (string).
Never throws (non-throwing design).

### `emk.setup.installExtra(name)`
Installs a Track 1 optional library into the Embedded Python (`python_env/`) via pip.  
`name`: `"pubchempy"` / `"mordred"` / `"biopython"` / `"meeko"` / `"vina"` / `"pdbfixer"` / `"torch"` / `"torch_geometric"` / `"transformers"`.
Versions are managed via `config/settings.json` under `extraLibraries.<name>`. Verifies via `import` after installation.  
`"torch_geometric"` uses a special flow: auto-detects the installed torch version, dynamically constructs the
PyG wheel URL (`https://data.pyg.org/whl/torch-<X.Y.Z>+cpu.html`), and also installs
`torch_scatter` / `torch_sparse` / `torch_cluster`.  
`"pdbfixer"` automatically fetches openmm >= 8.2 as a dependency (~70 MB).  
Errors: `notDesktop` (Online environment), `unknownLibrary`, `installFailed`, `importVerifyFailed`.

### `emk.setup.useExternal(pythonPath)`
For Track 2. Switches `pyenv` to the external CPython executable at `pythonPath`.
Idempotent within a session (logs a warning if Python is already Loaded).  
`pythonPath`: string or char (absolute path). Calling `validate()` afterwards is recommended.  
Errors: `invalidInput` (non-string/char or empty), `fileNotFound`, `pyenvFailed`.

### `emk.setup.validate(Libraries=string.empty)`
Diagnoses the presence and version of specified libraries using Python `importlib.metadata`.
When `Libraries` is omitted, checks RDKit + all currently installed Track 1 libraries.  
Returns: `table` — columns: `Library`(string), `Installed`(logical), `Version`(string), `Track`(string).
Never throws (non-throwing design).  
Errors: `invalidInput` (non-string array).

### `emk.setup.recipe(name)`
Displays the installation recipe (track classification, commands, license, notes) for the specified
library in the MATLAB Command Window.  
`name`: `"pubchempy"` / `"mordred"` / `"biopython"` / `"meeko"` / `"vina"` / `"pdbfixer"` /
`"docking"` (meeko+vina+pdbfixer combined) / `"openbabel"` / `"mdanalysis"` / `"pymol"`.  
Returns: nothing (display only).  
Errors: `unknownLibrary`.

### `emk.setup.installTrack2(name, BasePython="")`
Creates a venv at `python_env_t2/<name>/` for Track 2 libraries, runs pip install, writes
`settings.json` with `python.external_path`, and connects via `useExternal()` (ADR-007).  
`name`: `"mdanalysis"` (GPLv2+) / `"pymol"` (PSF/BSD). Open Babel is excluded (requires MSI — see `recipe("openbabel")`).  
`BasePython`: if omitted, auto-detects `py` / `python` from PATH.  
On subsequent sessions, `initPython()` auto-detects `external_path` and calls `useExternal()`.  
Errors: `notDesktop`, `unknownLibrary`, `basePythonNotFound`, `venvFailed`, `installFailed`, `importVerifyFailed`, `settingsWriteFailed`.

### `emk.setup.snapshot()`
Returns: `struct` — `.matlab` `.python` `.rdkit` `.commit` `.toolboxes` `.timestamp`.
Captures current environment for RF02 version lock. Never throws (non-throwing design).  
See also: `emk.setup.lockfile`, `emk.setup.verifyLock`.

### `emk.setup.lockfile(snap, filePath)` / `emk.setup.lockfile(filePath)`
Save mode (2 args): serialises `snap` struct as pretty-printed JSON to `filePath`.  
Load mode (1 string arg): reads a lock JSON file and returns the struct.  
Errors: `invalidInput`, `writeError`, `fileNotFound`.

### `emk.setup.verifyLock(lockRef)`
Compares the current environment against a saved RF02 lock.
`lockRef`: a JSON file path (string/char) or lock struct.  
Returns: `struct` — `.pass` (logical), `.details` (per-field: `.expected .actual .match`), `.warnings` (string).  
Critical fields: `matlab`, `python`, `rdkit`. Non-critical (warning only): `commit`.  
Errors: `invalidInput`.

---

## emk.mol

### `emk.mol.fromSmiles(smiles)`
Returns: `py.rdkit.Chem.rdchem.Mol` — RDKit Mol object (Python reference, ADR-002). Input: `string|char`.  
Errors: `invalidInput` (wrong type), `invalidSmiles` (empty or RDKit returns None).

### `emk.mol.toSmiles(mol)`
Returns: `string` — Canonical SMILES. Canonicalized via the Morgan algorithm; may differ from the input SMILES.  
Errors: `invalidInput` (non-Mol object), `rdkitError`.

### `emk.mol.isValid(smiles)`
Returns: `logical` — `true` if the SMILES is valid, `false` otherwise. Never throws (only rethrows `rdkitError`).  
Errors: `invalidInput` (non-string/char), `rdkitError` (only when RDKit is not running).

### `emk.mol.hasSubstruct(mol, query)`
Returns: `logical` or `logical(1,N)` — Substructure match result. `query` is a SMARTS string or Mol.
When `mol` is a cell array, returns a 1×N logical vector.  
Errors: `invalidMol`, `invalidQuery` (invalid SMARTS), `rdkitError`.

### `emk.mol.toStruct(mol, Format="molblock")`
Returns: `struct` — `.smiles` `.formula` `.numAtoms` `.numBonds` `.molblock` (or `.pickle`).
Use `Format="pickle"` for high-fidelity storage.  
Errors: `invalidInput`, `rdkitError`.

### `emk.mol.fromStruct(s)`
Returns: `py.rdkit.Chem.Mol` — Inverse of `toStruct`. Auto-detects and restores from molblock/pickle fields.  
Errors: `invalidInput` (not a struct or missing required fields), `rdkitError`.

### `emk.mol.toTable(mols, Properties=["SMILES","MolWt",...all10])`
Returns: `table` (N×M) — SMILES (string) + descriptors (double).
Invalid Mol rows have SMILES=`"<invalid>"` and descriptor=NaN.  
Errors: `invalidInput`, `unknownProperty`, `emptyProperties`.

### `emk.mol.scaffold(mol)`
Returns: `py.rdkit.Chem.rdchem.Mol` — Bemis-Murcko scaffold. Acyclic molecules return a 0-atom Mol (no exception).  
Errors: `invalidInput` (non-Mol), `rdkitError`.

---

## emk.descriptor

### `emk.descriptor.molWeight(mol)`
Returns: `double` — Average molecular weight [g/mol]. Uses `Descriptors.MolWt`
(sum of IUPAC average atomic weights for all atoms including implicit H).  
Errors: `invalidInput` (non-Mol).

### `emk.descriptor.calculate(mol, descriptorNames=<all10>)`
Returns: `struct` — field names = descriptor names, values = double (`MolFormula` only is string).  
Default (no argument) computes the following 10 numeric descriptors. `MolFormula` is returned only when explicitly requested.

| Name | Description | RDKit API | Type |
|---|---|---|---|
| `MolWt` | Average molecular weight (g/mol) | `Descriptors.MolWt` | double |
| `ExactMolWt` | Monoisotopic molecular weight | `Descriptors.ExactMolWt` | double |
| `LogP` | Wildman-Crippen LogP | `Descriptors.MolLogP` | double |
| `TPSA` | Topological polar surface area (Å²) | `Descriptors.TPSA` | double |
| `NumHAcceptors` | H-bond acceptor count | `rdMolDescriptors.CalcNumHBA` | double |
| `NumHDonors` | H-bond donor count | `rdMolDescriptors.CalcNumHBD` | double |
| `NumRotatableBonds` | Rotatable bond count | `rdMolDescriptors.CalcNumRotatableBonds` | double |
| `RingCount` | Total ring count (SSSR) | `rdMolDescriptors.CalcNumRings` | double |
| `FractionCSP3` | Fraction of sp3 carbons (Fsp3) | `Descriptors.FractionCSP3` | double |
| `HeavyAtomCount` | Heavy atom count | `mol.GetNumHeavyAtoms()` | double |
| `MolFormula` | Molecular formula string (e.g., "C2H6O") | `rdMolDescriptors.CalcMolFormula` | **string** |

Errors: `invalidInput`, `unknownDescriptor`.

### `emk.descriptor.batchCalculate(mols, descriptorNames=<all10>)`
Returns: `table` (N×M) — Batch descriptor calculation for multiple molecules.
Column names = descriptor names, values = double. Invalid mol rows are NaN.  
Errors: `invalidInput` (non-cell), `unknownDescriptor`, `allMolsFailed`.

### `emk.descriptor.mordred(mol, descriptorNames=[])`
Calculates Mordred 2D descriptors and returns a struct. Computes ~1800 descriptors when
`descriptorNames` is omitted. Failed descriptors are NaN.  
Requires: `emk.setup.installExtra("mordred")` must be installed first.  
Errors: `invalidInput` (non-Mol), `libraryNotFound`, `pythonError`.

### `emk.descriptor.mordredBatch(mols, descriptorNames=[])`
Returns: `table` (N×M) — Batch Mordred descriptor calculation. Uses a `run_mordred.py` helper
for a single IPC round trip (ADR-002 rev.3). Invalid mol rows are NaN.  
Requires: `emk.setup.installExtra("mordred")` must be installed.  
Errors: `invalidInput` (non-cell), `libraryNotFound`, `allMolsFailed`, `pythonError`.

### `emk.descriptor.mordredNames()`
Returns: `string(1×N)` — Sorted list of available 2D Mordred descriptor names (~1800 entries).
Use as the `descriptorNames` argument to `mordred()` / `mordredBatch()`.  
Requires: `emk.setup.installExtra("mordred")` must be installed.  
Errors: `libraryNotFound`, `pythonError`.

### `emk.descriptor.qed(mol)`
Returns: `double` ∈ [0,1] — QED (Quantitative Estimate of Drug-likeness) score. Combines 8
molecular properties via desirability functions. Score ≥ 0.67 is a common drug-likeness threshold.
Uses `rdkit.Chem.QED.qed`.  
Errors: `invalidInput` (non-Mol), `rdkitError`.

### `emk.descriptor.saScore(mol)`
Returns: `double` ∈ [1,10] — SA (Synthetic Accessibility) Score. Lower is easier to synthesise
(1-3: easy, 3-6: moderate, 6-10: difficult/impossible). Combines fragment contributions with
ring/stereo complexity penalties.  
Requires: `rdkit.Contrib.SA_Score.sascorer` on the Python path (included in standard RDKit).  
Errors: `invalidInput`, `rdkitError` (including module-not-found).

### `emk.descriptor.bcut(mol)`
Returns: `double(1,8)` — BCUT2D descriptors. Eigenvalues of the Burden matrix weighted by MW,
partial charge, LogP, and molar refractivity. Order: `BCUT2D_MWHI`, `BCUT2D_MWLOW`, `BCUT2D_CHGHI`, `BCUT2D_CHGLO`,
`BCUT2D_LOGPHI`, `BCUT2D_LOGPLOW`, `BCUT2D_MRHI`, `BCUT2D_MRLOW`. Uses `rdMolDescriptors.BCUT2D`.  
Errors: `invalidInput`, `rdkitError`.

### `emk.descriptor.fragmentCount(mol)`
Returns: `struct` — Ring system and functional group fragment counts. Fields: `NumRings`,
`NumAromaticRings`, `NumAliphaticRings`, `NumHeteroRings`, `NumCarbonyl`, `NumAmine`,
`NumHydroxyl`, `NumHalogen`, `NumNitrile`, `NumSulfonamide`, `NumAmide`.
Uses `rdMolDescriptors` + `Fragments` (fr_* SMARTS).  
Errors: `invalidInput`, `rdkitError`.

---

## emk.scaffold

### `emk.scaffold.genericMurcko(mol)`
Returns: `py.rdkit.Chem.rdchem.Mol` — Generic Murcko scaffold (all atoms set to C, all bonds to
single). Two-step: `GetScaffoldForMol` then `MakeScaffoldGeneric`. Acyclic molecules return a
0-atom Mol without error.  
Note: differs from `emk.mol.scaffold` which preserves atom/bond types.  
Errors: `invalidInput`, `rdkitError`.

### `emk.scaffold.brics(mol)`
Returns: `string(1,N)` — BRICS fragment SMILES array (non-deterministic order; frozenset origin).
Uses `rdkit.Chem.BRICS.BRICSDecompose`.  
Errors: `invalidInput`, `rdkitError`.

### `emk.scaffold.rgroup(mols, coreSmiles)`
Returns: `[tbl, unmatchedIdx]` — R-group decomposition result.  
`tbl`: `table` with columns `Core`, `R1`, `R2`, ... (string SMILES for each R-group).  
`unmatchedIdx`: `double(1,K)` — 1-based indices of molecules that did not match the core.  
Recommend using `[*:1]`-style R-mapped cores in `coreSmiles`.  
Errors: `invalidInput`, `invalidCore`, `noMatch` (all molecules unmatched), `rdkitError`.

---

## emk.dataset

### `emk.dataset.esol(CacheDir="", ForceDownload=false)`
Returns: `table` — ESOL (Delaney) aqueous solubility dataset (~1128 compounds). Columns: `SMILES`,
`Name` (string), `logS`, `logS_Delaney`, `MolWt` (double). Downloads `delaney-processed.csv` from
DeepChem S3; caches to `data/benchmark/esol.csv`.  
Errors: `invalidInput`, `downloadFailed`, `parseFailed`.

### `emk.dataset.freesolv(CacheDir="", ForceDownload=false)`
Returns: `table` — FreeSolv hydration free energy dataset (~643 compounds). Columns: `SMILES`,
`Name` (string), `DeltaG_exp`, `DeltaG_calc`, `DeltaG_exp_sem` (double).  
Errors: `invalidInput`, `downloadFailed`, `parseFailed`.

### `emk.dataset.bbbp(CacheDir="", ForceDownload=false)`
Returns: `table` — BBBP blood-brain barrier permeability dataset (~2050 compounds). Columns:
`SMILES`, `Name` (string), `BBB` (logical, true = permeable).  
Errors: `invalidInput`, `downloadFailed`, `parseFailed`.

### `emk.dataset.tox21(CacheDir="", ForceDownload=false)`
Returns: `table` — Tox21 multi-target toxicity dataset (~7831 compounds). Columns: `SMILES`,
`MolID` (string) + 12 toxicity endpoints (double, 0/1/NaN). Downloads `.csv.gz`, decompresses
with `gunzip`.  
Endpoints: `NR_AR`, `NR_AR_LBD`, `NR_AhR`, `NR_Aromatase`, `NR_ER`, `NR_ER_LBD`,
`NR_PPAR_gamma`, `SR_ARE`, `SR_ATAD5`, `SR_HSE`, `SR_MMP`, `SR_p53`.  
Errors: `invalidInput`, `downloadFailed`, `parseFailed`.

---

## emk.filter

### `emk.filter.lipinski(tbl, MaxViolations=0)`
Returns: `table` — Input table with `Violations_Ro5`(double) and `Pass_Ro5`(logical) columns appended.
Required columns: `MolWt`, `LogP`, `NumHDonors`, `NumHAcceptors`.

Ro5 thresholds: one violation each for MW>500, LogP>5, HBD>5, HBA>10.
Valid range for `MaxViolations` is [0,4].  
> NaN descriptors evaluate as `NaN > 500 = false`, so they count as "no violation" (false negative risk).
> Using `rmmissing(tbl)` as pre-processing is recommended.  
Errors: `invalidInput` (non-table), `invalidMaxViol` (outside [0,4] / NaN / Inf), `missingColumns`.

### `emk.filter.veber(tbl)`
Returns: `table` — Input table with `Violations_Veber`(double, 0-2) and `Pass_Veber`(logical) appended.
Required columns: `NumRotatableBonds`, `TPSA`. Pure MATLAB; no RDKit required.  
Veber criteria: `NumRotatableBonds > 10` or `TPSA > 140 Å²` counts as one violation each.  
Errors: `invalidInput`, `missingColumns`.

### `emk.filter.pains(tbl)`
Returns: `table` — Input table with `NumPainsAlerts`(double), `PainsAlerts`(string, comma-separated),
and `HasPains`(logical) columns appended. Required column: `SMILES`.
Uses RDKit `FilterCatalog(PAINS)`.  
Errors: `invalidInput`, `missingColumns`, `rdkitError`.

### `emk.filter.reos(tbl)`
Returns: `table` — Input table with `Violations_REOS`(double, 0-6) and `Pass_REOS`(logical) appended.
Required columns: `MolWt`, `LogP`, `NumHDonors`, `NumHAcceptors`, `NumRotatableBonds`, `HeavyAtomCount`.
Pure MATLAB; no RDKit required.  
REOS ranges (6/7 criteria; FormalCharge excluded): MW [200,500], LogP [-5,5], HBD [0,5],
HBA [0,10], RotBonds [0,8], HeavyAtoms [15,50].  
Errors: `invalidInput`, `missingColumns`.

---

## emk.cluster

### `emk.cluster.butina(fps, Threshold=0.2, Metric="tanimoto")`
Returns: `cell(1,C)` — C clusters. Each element is a `double(1,K)` array of 1-based indices
into `fps`. `clusters{1}` is the largest cluster; the first index in each cluster is the centroid.  
Butina sphere-exclusion clustering. `Threshold` is a Tanimoto **distance** threshold
(80% similarity cutoff at default 0.2).  
`Metric`: `"tanimoto"` only.  
Errors: `invalidInput`, `invalidThreshold` (outside (0,1]), `invalidMetric`, `rdkitError`.

---

## emk.diversity

### `emk.diversity.pick(fps, N, Metric="tanimoto", Seed=-1)`
Returns: `double(1,N)` — 1-based indices of N maximally diverse molecules selected from `fps`.  
MaxMin (Kennard-Stone variant) algorithm via `MaxMinPicker.LazyBitVectorPick`.  
`Seed=-1`: random start; `Seed≥0`: reproducible.  
`Metric`: `"tanimoto"` only.  
Errors: `invalidInput`, `invalidN` (outside [1,M]), `invalidMetric`, `rdkitError`.

---

## emk.conformer

### `emk.conformer.embed(mol, Method="ETKDGv3", RandomSeed=-1)`
Returns: `py.rdkit.Chem.rdchem.Mol` — Molecule with one 3D conformer (hydrogens removed).
Pipeline: `AddHs` → `EmbedMolecule(params)` → `RemoveHs`.  
`Method`: `"ETKDGv3"` (recommended), `"ETKDGv2"`, `"ETKDG"`, `"KDG"`.
`RandomSeed=-1`: random; `RandomSeed≥0`: reproducible.  
Errors: `invalidInput`, `invalidMethod`, `embeddingFailed` (RDKit returned -1), `rdkitError`.

### `emk.conformer.optimize(mol, ForceField="MMFF94", MaxIter=2000)`
Returns: `py.rdkit.Chem.rdchem.Mol` — Force-field minimized molecule.
Pipeline: `AddHs(addCoords=true)` → minimize → `RemoveHs`.  
`ForceField`: `"MMFF94"` (recommended for drug-like molecules) or `"UFF"` (fallback for unusual elements).
Non-convergence (status=1) logs a warning only; FF setup failure (status=-1) raises `optimizeFailed`.  
Errors: `invalidInput` (non-Mol or no conformer), `invalidForceField`, `optimizeFailed`, `rdkitError`.

---

## emk.shape

### `emk.shape.compare(mol1, mol2, Method="protrude")`
Returns: `double` ∈ [0,1] — 3D shape similarity score (higher = more similar). Both molecules
must have at least one 3D conformer (from `emk.conformer.embed`).  
`Method="protrude"`: `1 - ShapeProtrudeDist(mol1, mol2, allowReordering=false)`.  
`Method="tanimoto"`: `ShapeTverskyIndex(mol1, mol2, 1.0, 1.0)` (shape Tanimoto).  
Score is clamped to [0,1] for floating-point safety.  
Errors: `invalidInput` (non-Mol or no conformer), `invalidMethod`, `rdkitError`.

---

## emk.fingerprint

### `emk.fingerprint.morgan(mol, Radius=2, NBits=2048)`
Returns: `py.rdkit.DataStructs.ExplicitBitVect` — Morgan (ECFP) fingerprint. Radius=2 ≈ ECFP4.  
Errors: `invalidInput`, `rdkitError`.

### `emk.fingerprint.maccs(mol)`
Returns: `py.rdkit.DataStructs.ExplicitBitVect` — 167-bit MACCS keys fingerprint (bit 0 is unused).  
Errors: `invalidInput`, `rdkitError`.

### `emk.fingerprint.toArray(fp)`
Returns: `logical(1,N)` — Converts a fingerprint to a MATLAB logical array. N = `fp.GetNumBits()`.  
Errors: `invalidInput` (not a Python object or does not have `ToBitString()`).

---

## emk.similarity

### `emk.similarity.tanimoto(fp1, fp2)`
Returns: `double` ∈ [0,1] — Tanimoto coefficient. Returns 1.0 for identical fingerprints.
Mismatched bit lengths raise `rdkitError`.  
Errors: `invalidInput` (non-Python object), `rdkitError`.

### `emk.similarity.dice(fp1, fp2)`
Returns: `double` ∈ [0,1] — Dice coefficient. For binary vectors, Dice ≥ Tanimoto always holds.  
Errors: `invalidInput`, `rdkitError`.

### `emk.similarity.rankBy(queryFp, dbFps, N=Inf, Metric="tanimoto")`
Returns: `struct` — `.Indices(1×K)` `.Scores(1×K)` `.Metric(string)`. Sorted in descending score order.
Uses `BulkTanimotoSimilarity` (single IPC round trip).  
Errors: `invalidQueryFp`, `invalidDbFps`, `invalidN` (0/negative/non-integer), `invalidMetric`, `rdkitError`.

### `emk.similarity.matrix(fps, Metric="tanimoto")`
Returns: `double(N×N)` — Symmetric similarity matrix. Diagonal = 1.0. Floating-point errors removed via `(S+S')/2`.  
Errors: `invalidInput` (non-empty cell), `invalidMetric`, `rdkitError`.

---

## emk.filter

### `emk.filter.lipinski(tbl, MaxViolations=0)`
Returns: `table` — Input table with `Violations_Ro5`(double) and `Pass_Ro5`(logical) columns appended.
Required columns: `MolWt`, `LogP`, `NumHDonors`, `NumHAcceptors`.

Ro5 thresholds: one violation each for MW>500, LogP>5, HBD>5, HBA>10.
Valid range for `MaxViolations` is [0,4].  
> NaN descriptors evaluate as `NaN > 500 = false`, so they count as "no violation" (false negative risk).
> Using `rmmissing(tbl)` as pre-processing is recommended.  
Errors: `invalidInput` (non-table), `invalidMaxViol` (outside [0,4] / NaN / Inf), `missingColumns`.

### `emk.filter.veber(tbl)`
Returns: `table` — Appends `Violations_Veber`(double, 0-2) and `Pass_Veber`(logical). Pure MATLAB.
Required: `NumRotatableBonds`, `TPSA`. Violations: `>10` or `>140 Å²`.  
Errors: `invalidInput`, `missingColumns`.

### `emk.filter.pains(tbl)`
Returns: `table` — Appends `NumPainsAlerts`, `PainsAlerts` (string), `HasPains` (logical).
Required: `SMILES`. Uses RDKit FilterCatalog(PAINS).  
Errors: `invalidInput`, `missingColumns`, `rdkitError`.

### `emk.filter.reos(tbl)`
Returns: `table` — Appends `Violations_REOS`(double, 0-6) and `Pass_REOS`(logical). Pure MATLAB.
Required: `MolWt`, `LogP`, `NumHDonors`, `NumHAcceptors`, `NumRotatableBonds`, `HeavyAtomCount`.
REOS ranges (6/7; FormalCharge excluded): MW[200,500], LogP[-5,5], HBD[0,5], HBA[0,10],
RotBonds[0,8], HeavyAtoms[15,50].  
Errors: `invalidInput`, `missingColumns`.

---

## emk.io

### `emk.io.readSdf(filePath)`
Returns: `1×N cell` — Mol cell array loaded from an SDF file.
Uses `SDMolSupplier(removeHs=True, sanitize=True)`. Failed molecules are skipped with `logWarn`.  
Errors: `invalidInput`, `fileNotFound`, `rdkitError`.

### `emk.io.writeSdf(mols, filePath)`
Writes molecules to an SDF file. The parent directory must exist. Overwrites existing files.
`writer.close()` is always called, even on error.  
Errors: `invalidInput`, `invalidMol` (non-Python-object element), `dirNotFound`, `rdkitError`.

### `emk.io.readSmilesList(filePath)`
Returns: `1×N cell` — Reads a one-SMILES-per-line file. Skips `#` comments and blank lines.
The name column after tab/space is ignored. Failed lines are skipped with `logWarn`.  
Errors: `invalidInput`, `fileNotFound`, `allLinesFailed`.

---

## emk.viz

### `emk.viz.draw2d(mol, Title="", Width=300, Height=300)`
Returns: `matlab.ui.Figure` — Renders a 2D structure in a MATLAB figure.
Pipeline: `Compute2DCoords` → `MolToFile` → `imread` → `imshow` (via temp file).  
Errors: `invalidInput` (non-Mol or Width/Height < 1), `rdkitError`.

---

## emk.db

### `emk.db.searchPubchem(query, Type="name")`
Returns: `table` (5 columns: `CID` double, `IUPACName` `MolecularFormula` `IsomericSMILES` string, `MolecularWeight` double).  
`Type`: `"name"` / `"smiles"` / `"cid"` / `"inchikey"`. Uses `webread` only (no Python required).  
Errors: `invalidInput`, `invalidType`, `notFound` (HTTP 404), `networkError`.

### `emk.db.pubchemFetch(identifier, NameSpace="name", MaxSynonyms=10)`
Uses PubChemPy to fetch an extended property set from PubChem (synonyms, InChI/InChIKey, XLogP, TPSA, HBD/HBA, etc.) and returns a struct.
An enhanced version of `searchPubchem` (requires Python).  
`NameSpace`: `"name"` / `"smiles"` / `"cid"` / `"inchi"` / `"inchikey"` / `"formula"`. Also accepts a numeric scalar CID.  
Main struct fields: `CID`(double), `IUPACName` `MolecularFormula` `IsomericSMILES` `InChI` `InChIKey`(string),
`MolecularWeight` `XLogP` `TPSA` `HBondDonors` `HBondAcceptors` `RotatableBonds` `HeavyAtomCount` `Charge` `Complexity`(double),
`Synonyms`(string array).  
Requires: `emk.setup.installExtra("pubchempy")` must be installed first.  
Errors: `invalidInput`, `invalidNamespace`, `libraryNotFound`, `notFound`, `pythonError`.

### `emk.db.searchChembl(query, Type="name")`
Returns: `table` (8 columns: `ChEMBLID` `Name` `MolecularFormula` `SMILES` `InChIKey` string,
`MolecularWeight` `ALogP` `HBondDonors` `HBondAcceptors` double).  
`Type`: `"name"` / `"chemblid"` / `"smiles"` / `"inchikey"`. Uses `webread` only.  
Errors: `invalidInput`, `invalidType`, `notFound`, `networkError`.

### `emk.db.searchChemblTarget(query, TargetType="SINGLE PROTEIN", MaxRows=10)`
Returns: `table` (4 columns: `TargetChEMBLID` `PreferredName` `Organism` `TargetType` string).
Partial-match search via `pref_name__icontains`. Uses `webread` only (no Python required).  
Errors: `invalidInput`, `invalidOptions`, `notFound`, `networkError`.

### `emk.db.getChemblActivity(targetId, ActivityType="IC50", MaxRows=50, MinActivity_nM=Inf)`
Returns: `table` (5 columns: `MoleculeChEMBLID` `Name` `SMILES` `ActivityType` string, `Value_nM` double).
Returns only rows where `standard_relation="="` and units are nM. Retains only rows with `Value_nM <= MinActivity_nM` (default: all rows).
Uses `webread` only.  
Errors: `invalidInput`, `invalidOptions`, `notFound`, `networkError`.

---

## emk.repro

### `emk.repro.verify(metrics, criteria)`
Validates reproduction metrics against RF03 acceptance criteria.  
`metrics`: struct with numeric metric fields (e.g. `metrics.rmse_cv = 1.017`).  
`criteria`: struct of acceptance bounds; each field is a sub-struct with `"upper"` and/or `"lower"` keys.  
Returns: `struct` — `.pass` (logical), `.details` (per-metric: `.value .pass .criteria`), `.report` (char, formatted summary).  
Errors: `invalidInput` (non-struct), `missingMetric` (criteria field absent from metrics).

---

## emk.util

### `emk.util.isOnline()`
Returns: `logical` — `true` on MATLAB Online. Detection order: `ismatlabonline()` → `MATLAB_ONLINE` env var → Linux x64 heuristic.

### `emk.util.benchmarkBatch(smilesList, descriptorNames=<all10>)`
Returns: `struct` — `.nMols` `.parseSec` `.batchSec` `.totalSec` `.molsPerSec` `.tbl`. For throughput measurement.  
Errors: `invalidInput`, `emptyInput`.

---

## src/config

### `emkLoadConfig()`
Returns: `struct` — `.python` `.rdkit` `.runtime` `.output` `.run`.
Priority: `EMK_<SECTION>_<KEY>` environment variables > `config/settings.json` > defaults.

---

## src/util

### `logInfo(msg, ...)` / `logWarn(msg, ...)` / `logError(msg, ...)` / `logDebug(msg, ...)`
Outputs in `[HH:MM:SS][LEVEL]  message` format. `logDebug` only outputs when `EMK_LOG_VERBOSE=1`.
Direct use of `fprintf` is prohibited.

### `logProgress(i, n, label)`
Displays a loop progress bar as a percentage of `i/n`.

### `logSection(scriptId, label, layer)`
Banner log called at the beginning of each `%%` section in tutorial scripts.  
Example output: `[11:36:40][INFO]  --- R01 | Section 0: Setup  [Research L4] ---`  
Arguments: `scriptId` = "F01"/"S01"/"R01"/"A01", `label` = section header text,
`layer` = "Foundation L1" / "Stories L2" / "Analytics L3" / "Research L4".

### `makeRunDir(Prefix="", BaseDir="result/runs")`
Returns: `char` — Creates and returns the path `result/runs/yyyyMMdd_HHmmss[_Prefix]`.
Direct use of `mkdir` is prohibited.
