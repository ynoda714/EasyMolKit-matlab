# Contributing

Thank you for your interest in contributing to EasyMolKit!

## How to Contribute

### Reporting Bugs

Please use the [Bug Report](../../issues/new?template=bug_report.yml) issue template. Include:
- Steps to reproduce
- Expected vs actual behavior
- MATLAB version and OS
- RDKit version (`emk.setup.verify()` output)

### Feature Requests

Please use the [Feature Request](../../issues/new?template=feature_request.yml) issue template.

### Pull Requests

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature`)
3. Follow the coding conventions below
4. Add tests for new functionality
5. Ensure all tests pass
6. Submit a Pull Request

## Coding Conventions

- **Language**: All `.m` file content (comments, logs, errors) must be in English. Exception: `examples/japanese/` allows Japanese comments (educational content only; variable names remain English)
- **Naming**: `camelCase` for functions, `Test<Feature>.m` for test classes
- **Logging**: Use log helpers (`logInfo`, `logWarn`, `logError`). Never use `fprintf` directly
- **File placement**: All logic under `src/`. Domain logic in `src/+emk/+<module>/`
- **Artifacts**: Output via `makeRunDir()` to `result/runs/`. Never hardcode paths
- **Python**: Never call `py.rdkit.*` directly — use `emk.*` wrappers

## Development Setup

```matlab
addpath(genpath("src"));
emk.setup.install();  % First time only

% Run unit tests
suite = testsuite("tests/unit");
runner = matlab.unittest.TestRunner.withNoPlugins;
results = runner.run(suite);
```

## Submitting a Reproduction (RF04)

Reproducing a published cheminformatics result is one of the most valuable contributions.

### Required Components

Each reproduction must include:

1. **RF01 README** — follow the standard template at `repro/TEMPLATE.en.md`.
   - Overview, Environment, Data, Script, Result, Verification, Discussion sections
   - Descriptor definitions table (tool + version + definition for every descriptor used)

2. **RF02 Lock** — capture and save a version lock with `emk.setup.snapshot()`:
   ```matlab
   snap = emk.setup.snapshot();
   emk.setup.lockfile(snap, fullfile(runDir, "lock_snapshot.json"));
   ```
   The lock records MATLAB, Python, RDKit, and Toolbox versions.

3. **RF03 Verification** — report pass/fail against numerical criteria:
   ```matlab
   crit.rmse_cv = struct("upper", 1.20);
   crit.r2_cv   = struct("lower", 0.75);
   met.rmse_cv  = rmseCV;
   met.r2_cv    = r2CV;
   result = emk.repro.verify(met, crit);
   ```
   Save the result in `metrics.json` with `rf03_pass` and `rf03_criteria` fields.

4. **Script** — place the reproduction script in `repro/rp<XX>_<name>/`.
   The script must be a plain `.m` file with `%% Section` markers (no interactive input).

### Directory Structure

```
repro/
  rp<XX>_<name>/
    README.md           (RF01 -- based on repro/TEMPLATE.en.md)
    rp<XX>_<name>.m     (reproduction script)
    lock_template.json  (RF02 schema)
    result/             (.gitignore managed -- not committed)
```

### Automated Checks

Before submitting a pull request, verify:
- [ ] `repro/<folder>/README.md` contains all RF01 sections
- [ ] `lock_snapshot.json` exists in at least one `result/runs/<ts>/` folder
- [ ] `metrics.json` contains `rf03_pass` (true/false)
- [ ] All unit tests pass: `testsuite("tests/unit")`

### Community Labels

- **`Good First Reproduction`**: lightweight Tier A papers (e.g., linear regression, descriptor calculation)
- **`Good First Issue`**: doc fixes, test additions, descriptor additions, dataset additions

## Code of Conduct

Please be respectful and constructive in all interactions.
