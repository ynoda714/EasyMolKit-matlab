# Test Catalog — EasyMolKit

> 30 unit test classes + 2 smoke scripts.
> Japanese version → [../test_catalog.md](../test_catalog.md)

## Running Tests

```matlab
addpath(genpath("src"));
addpath(genpath("tests"));   % required for class-name lookup

% Run all unit tests
suite   = testsuite("tests/unit");
runner  = matlab.unittest.TestRunner.withNoPlugins;
results = runner.run(suite);
fprintf("RESULT: %d PASS / %d FAIL / %d Total\n", ...
    sum([results.Passed]), sum([results.Failed]), numel(results));

% Run smoke tests
run("tests/smoke/test_mvp_smoke.m");
run("tests/smoke/test_m2_smoke.m");
```

---

## Unit Tests (`tests/unit/`) — 30 classes

| Class | Target Module | Scope |
|---|---|---|
| `TestCluster` | `emk.cluster` | Butina sphere-exclusion clustering, threshold, validation |
| `TestConformer` | `emk.conformer` / `emk.shape` | embed/optimize/compare: validation, integration, seed reproducibility |
| `TestDataset` | `emk.dataset` | ESOL / FreeSolv / BBBP / Tox21 validation, download guards |
| `TestDb` | `emk.db` | PubChem / ChEMBL search, fetch, error handling |
| `TestDescriptor` | `emk.descriptor` | 10 standard descriptors + qed / saScore / bcut(1×8) / fragmentCount; batch, NaN for invalid mol |
| `TestDiversity` | `emk.diversity` | MaxMin diversity selection, Seed reproducibility, validation |
| `TestFilter` | `emk.filter` | Lipinski / Veber / REOS / PAINS: violation count, pass/fail, edge cases |
| `TestFingerprint` | `emk.fingerprint` | Morgan (radius, nBits), MACCS, `toArray` type/shape |
| `TestInitPython` | `emk.setup.initPython` | Platform detection, idempotent `pyenv` configuration |
| `TestInstall` | `emk.setup.install` | Embedded Python deployment, path length guard |
| `TestInstallExtra` | `emk.setup.installExtra` | Track 1 library install and import verification |
| `TestInstallOnline` | `emk.setup.installOnline` | MATLAB Online bootstrap, skip-if-current logic |
| `TestInstallTrack2` | `emk.setup.installTrack2` | venv creation, `settings.json` write, useExternal call |
| `TestIo` | `emk.io` | SDF round-trip, SMILES list read, invalid mol skip |
| `TestIsOnline` | `emk.util.isOnline` | Desktop / Online detection heuristics |
| `TestLoadConfig` | `emkLoadConfig` | env var override, JSON loading, defaults |
| `TestLogHelpers` | `logInfo / logWarn / logError / logDebug` | Output format, verbose flag |
| `TestMakeRunDir` | `makeRunDir` | Directory creation, timestamp format, prefix |
| `TestMol` | `emk.mol` | fromSmiles, toSmiles, isValid, hasSubstruct, scaffold |
| `TestMordred` | `emk.descriptor.mordred*` | ~1800 descriptor names, batch, NaN handling |
| `TestPubchemFetch` | `emk.db.pubchemFetch` | Extended fetch, all NameSpace types, not-found error |
| `TestRdkitModule` | `emk` (general) | RDKit availability check, version string |
| `TestRecipe` | `emk.setup.recipe` | Output for all known library names, unknownLibrary error |
| `TestScaffold` | `emk.scaffold` | genericMurcko / brics / rgroup: validation, SMARTS matching |
| `TestSimilarity` | `emk.similarity` | tanimoto, dice, rankBy (N, Metric), matrix symmetry |
| `TestToTable` | `emk.mol.toTable` | Column names, types, NaN for invalid mol, empty input |
| `TestUseExternal` | `emk.setup.useExternal` | Path validation, fileNotFound, idempotent skip |
| `TestValidate` | `emk.setup.validate` | Table structure, all-libraries scan, specific names |
| `TestVerify` | `emk.setup.verify` | Struct fields, logical values, no-throw guarantee |
| `TestViz` | `emk.viz` | draw2d returns Figure, invalid input errors |

### Testing Conventions

- Tests that do not require RDKit are ordered **before** RDKit-dependent tests within each class.
- Every RDKit-dependent test method begins with:
  ```matlab
  tc.assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available");
  ```
- All tests use `verifyClass` / `verifySize` in addition to value checks (type/shape bugs are caught explicitly).
- No direct `py.rdkit.*` calls in test code — always use `emk.*` wrappers.

---

## Smoke Tests (`tests/smoke/`) — 2 scripts

| Script | Purpose | Key Steps |
|---|---|---|
| `test_mvp_smoke.m` | Core end-to-end workflow | `install` → `fromSmiles` → `calculate` → `tanimoto` → `draw2d` |
| `test_m2_smoke.m` | Milestone 2 workflow | `readSdf` / `readSmilesList` → `batchCalculate` → `lipinski` → `searchPubchem` / `searchChembl` |
