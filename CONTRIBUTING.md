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

## Code of Conduct

Please be respectful and constructive in all interactions.
