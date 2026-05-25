# Test Catalog — EasyMolKit

🇯🇵 [日本語](ja/test_catalog_ja.md)

> Complete list of test classes (25 classes, 869 tests). For how to run the test suite, see [quickstart.md](quickstart.md#running-tests).

```matlab
% Run tests
addpath(genpath("src"));
addpath(genpath("tests"));
suite   = testsuite("tests/unit");
runner  = matlab.unittest.TestRunner.withNoPlugins;
results = runner.run(suite);
fprintf("RESULT: %d PASS / %d FAIL / %d Total\n", ...
    sum([results.Passed]), sum([results.Failed]), numel(results));
```

---

## Test Classes

| Class | Tests | Module | Key Coverage |
|---|:---:|---|---|
| `TestMol` | 97 | `emk.mol` | SMILES parsing, canonicalization, substructure search, Mol serialization, scaffold extraction |
| `TestDb` | 95 | `emk.db` | PubChem / ChEMBL search & activity retrieval (mock + live server) |
| `TestSimilarity` | 93 | `emk.similarity` | Tanimoto / Dice coefficients, ranking search, N×N similarity matrix |
| `TestDescriptor` | 83 | `emk.descriptor` | 10 standard descriptors, batch calculation, invalid Mol handling |
| `TestFilter` | 72 | `emk.filter` | Lipinski Ro5 filter, violation count, NaN handling |
| `TestFingerprint` | 64 | `emk.fingerprint` | Morgan / MACCS generation, `toArray()` conversion, bit-length validation |
| `TestIo` | 49 | `emk.io` | SDF read/write, SMILES file loading, invalid Mol skipping |
| `TestPubchemFetch` | 36 | `emk.db` | PubChemPy detailed fetch (CID / name / SMILES / InChIKey) |
| `TestToTable` | 34 | `emk.mol` | `toTable()` column structure & type validation, NaN substitution for invalid Mol |
| `TestMordred` | 27 | `emk.descriptor` | Mordred 2D descriptor calculation (~1800 descriptors), batch processing |
| `TestInstallExtra` | 24 | `emk.setup` | Track 1 library install, verification, version management |
| `TestRdkitModule` | 25 | `emk.setup` | RDKit submodule importlib check, version retrieval |
| `TestRecipe` | 23 | `emk.setup` | `recipe()` display, unknown library error |
| `TestValidate` | 22 | `emk.setup` | `validate()` table format & installation state verification |
| `TestLoadConfig` | 15 | `src/config` | Config load priority (env var > JSON > defaults) |
| `TestInstallTrack2` | 16 | `emk.setup` | Track 2 venv creation, connection, settings write |
| `TestLogHelpers` | 12 | `src/util` | Log output format (INFO / WARN / ERROR / DEBUG) |
| `TestVerify` | 12 | `emk.setup` | Non-throwing behavior of `verify()`, Python / RDKit state detection |
| `TestInitPython` | 10 | `emk.setup` | `initPython()` idempotency, pyenv configuration, platform detection |
| `TestUseExternal` | 10 | `emk.setup` | External CPython path specification, invalid path errors |
| `TestIsOnline` | 6 | `emk.util` | MATLAB Online / Desktop environment detection |
| `TestInstall` | 7 | `emk.setup` | `install()` Desktop deployment flow |
| `TestInstallOnline` | 7 | `emk.setup` | `installOnline()` Online deployment flow |
| `TestViz` | 25 | `emk.viz` | `draw2d()` rendering, Figure return, input validation |
| `TestMakeRunDir` | 5 | `src/util` | `makeRunDir()` directory creation, timestamp format |

**Total: 25 classes / 869 tests**

---

## Test Categories

Each test uses `assumeTrue` so that **RDKit-independent tests** run first.

| Category | Examples |
|---|---|
| No RDKit required | Input validation, config loading, logging, utilities |
| RDKit required | Molecular operations, descriptor calculation, fingerprints, similarity |
| Network required | PubChem / ChEMBL live server tests (subset of `emk.db`) |

> In environments without RDKit, `assumeTrue(tc.rdkitAvailable(), "Skipped: RDKit not available")`
> automatically skips affected tests.

---

## Related Documents

- [function_reference.md](function_reference.md) — Full function signatures and error specs
- [algorithm_guide.md](algorithm_guide.md) — Algorithm rationale and test validation strategy
