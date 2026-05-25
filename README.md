# EasyMolKit

🇯🇵 [日本語](docs/ja/README.ja.md)

**An integrated Chemoinformatics environment that makes RDKit easy to use from MATLAB**

> Clone the repo → run one setup command → RDKit is ready to use in MATLAB — Python auto-deployed, no manual setup needed

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=ynoda714/EasyMolKit-matlab&file=main_rdkit.m)

## Why EasyMolKit?

In the Chemoinformatics field, Python + RDKit is the de-facto standard toolchain.
However, Python must be installed and configured first, and the setup process involves several barriers:

- Python pre-installation required — even before the first line of chemistry code
- Managing Python versions and virtual environments
- Complex RDKit installation via conda/pip
- Environment conflicts with commercial tools (e.g., PyMOL)
- The Python ecosystem is unfamiliar to MATLAB users

**EasyMolKit removes these barriers.** Leveraging MATLAB's `pyenv`-based Python integration,
users can access RDKit functionality as standard MATLAB functions — no manual Python setup or knowledge required.

## Features

- **Zero configuration**: No Python pre-installation needed — one call to `emk.setup.install()` automatically downloads and deploys an isolated Embedded Python + RDKit
- **MATLAB native**: Results are returned as MATLAB `table` / `struct` / `double` — immediately usable in your workspace
- **Desktop & Online**: Supports Windows Desktop and MATLAB Online (macOS / Linux Desktop untested)
- **RDKit wrappers**: Molecular analysis, descriptor calculation, fingerprints, and similarity search — all as MATLAB functions
- **Future-ready**: 3D structure visualization via PyMOL Open-Source planned for a future release

## Target Users

- **Chemistry, pharmacy, and medical researchers** who use MATLAB as their primary research environment
- **Students** learning Chemoinformatics
- **MATLAB users** who want to avoid spending time on Python environment setup

## Requirements

| Item | Desktop | MATLAB Online |
|---|---|---|
| MATLAB | R2025b or later | Latest release (managed by MathWorks) |
| Python | No setup needed — auto-deployed (Embedded Python) | Pre-installed |
| RDKit | Auto-deployed | Installed via `emk.setup.installOnline()` |
| Internet | Required for initial setup (auto-downloads Python + RDKit) | Required for `emk.setup.installOnline()` |
| OS | Windows | — |

## Quick Start

```matlab
% 1. Clone the repository
%    git clone https://github.com/ynoda714/EasyMolKit-matlab.git
%    cd EasyMolKit-matlab

% 2. Open in MATLAB and run the one-time setup
addpath(genpath("src"));
emk.setup.install();          % Desktop
% emk.setup.installOnline();  % MATLAB Online

% 3. Try it out
mol = emk.mol.fromSmiles("CCO");          % Ethanol
mw  = emk.descriptor.molWeight(mol);      % Molecular weight
fp  = emk.fingerprint.morgan(mol);        % Morgan fingerprint

% 4. Compare similarity between molecules
mol2  = emk.mol.fromSmiles("CCCO");       % Propanol
fp2   = emk.fingerprint.morgan(mol2);
score = emk.similarity.tanimoto(fp, fp2); % Tanimoto coefficient
```

For more details, see [docs/quickstart.md](docs/quickstart.md).

## Additional Libraries (Track 1 & Track 2)

EasyMolKit manages add-on libraries via two tracks.

| Track | Libraries | Installation | License |
|---|---|---|---|
| **Track 1** | pubchempy, mordred, biopython, torch, torch_geometric, transformers, datasets, etc. | `emk.setup.installExtra()` — added directly to Embedded Python | MIT / BSD-3 / Apache-2.0 |
| **Track 2** | Open Babel, MDAnalysis, PyMOL OSS | Requires a separate CPython environment; connect with `emk.setup.useExternal()` | GPLv2 / GPLv2+ / BSD |

## Tutorials & Examples (4-Layer Structure)

EasyMolKit provides progressive learning content under `examples/`.

| Layer | Audience | Required Toolbox | Content | Release |
|---|---|---|---|---|
| **L1 Foundation** | All users | None | Learn one API concept at a time (6 modules, 5–15 min each) | ✅ v1.0.0 |
| **L2 Application Stories** | After Foundation | None | Practical workflows combining multiple features (7 modules, 20–40 min each) | 🔜 v1.1.0 |
| **L3 Analytics** | All users | Varies (Statistics and ML, etc.) | QSAR, clustering, MS analysis, optimization, and more (A01–A10, 10 modules, 30–60 min each) | 🔜 v1.2.0 |
| **L4 Research** | All users | Varies (Parallel Computing, etc.) | Research-level applications (R01–R10, 30–90 min each) | 🔜 v1.3.0 |

### Layer 1: Foundation (Base MATLAB only)

| # | Title | Required Toolbox | What You Learn | Desktop | Online |
|---|---|---|---|:---:|:---:|
| F01 | Drawing Molecules with SMILES | None | Molecular representation, SMILES syntax | ✔ | ✔ |
| F02 | Calculating Molecular Properties | None | MW / LogP / TPSA — meaning and calculation | ✔ | ✔ |
| F03 | Introduction to Fingerprints | None | Bit vector representation, Morgan vs MACCS | ✔ | ✔ |
| F04 | Comparing Molecules by Similarity | None | Quantifying Tanimoto / Dice similarity | ✔ | ✔ |
| F05 | Substructure Search | None | SMARTS pattern matching | ✔ | ✔ |
| F06 | Reading Molecules from Files | None | SDF / SMILES file operations | ✔ | ✔ |

## Directory Structure

```
EasyMolKit-matlab/
├─ main_rdkit.m               # RDKit setup & basic operations (run section by section)
├─ config/
│   └─ settings.example.json  # Configuration template
├─ examples/
│   ├─ japanese/              # Distribution materials — Japanese (plain-text Live Code)
│   │   ├─ foundation/        #   L1: API basics (Base MATLAB only)
│   │   ├─ stories/           #   L2: Application stories
│   │   ├─ analytics/         #   L3: Statistics & ML integration
│   │   └─ research/          #   L4: Research level
│   └─ english/               # Distribution materials — English (comments differ only)
├─ src/
│   ├─ +emk/                  # Main package
│   │   ├─ +setup/            # Python environment setup
│   │   ├─ +mol/              # Molecular object operations
│   │   ├─ +descriptor/       # Descriptor calculation
│   │   ├─ +fingerprint/      # Fingerprint generation
│   │   ├─ +similarity/       # Similarity calculation
│   │   ├─ +filter/           # Molecular filtering (Lipinski, etc.)
│   │   ├─ +io/               # File I/O
│   │   ├─ +viz/              # Visualization
│   │   └─ +util/             # Package-level utilities
│   ├─ config/                # Configuration loader
│   └─ util/                  # Log helpers & common utilities
├─ result/                    # Run artifacts (not tracked by Git)
├─ tests/
│   ├─ unit/                  # matlab.unittest class-based tests
│   └─ smoke/                 # Smoke tests
├─ data/                      # Curated sample data
└─ docs/                      # Documentation
    └─ ja/                    # Japanese documentation
```

## API Overview

| Module | Example functions | Description |
|---|---|---|
| `emk.setup` | `install()`, `verify()`, `initPython()` | Automatic Python environment deployment & initialization |
| `emk.mol` | `fromSmiles()`, `toSmiles()`, `isValid()`, `hasSubstruct()` | Molecular object creation & conversion |
| `emk.descriptor` | `molWeight()`, `calculate()`, `batchCalculate()` | Molecular descriptor calculation |
| `emk.fingerprint` | `morgan()`, `maccs()`, `toArray()` | Fingerprint generation |
| `emk.similarity` | `tanimoto()`, `dice()` | Molecular similarity calculation |
| `emk.io` | `readSdf()`, `writeSdf()`, `readSmilesList()` | SDF / SMILES file I/O |
| `emk.viz` | `draw2d()` | 2D structure rendering |

For full API details, see [docs/function_reference.md](docs/function_reference.md).

## License

EasyMolKit: [MIT License](LICENSE)

### Third-party licenses

| Library | License | Purpose |
|---|---|---|
| RDKit | BSD-3-Clause | Chemoinformatics core |
| Python (CPython) | PSF License | Runtime environment |
| PyMOL Open-Source | Python/BSD-like | 3D visualization (future release) |

For details, see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and [docs/compliance.md](docs/compliance.md).

## Contributing

Bug reports, feature requests, and pull requests are welcome.
See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

## Disclaimer

- This project is developed and maintained as a personal leisure activity
- Please verify the accuracy of results at your own responsibility
- Support is best-effort; responses may be delayed
- This software is provided "AS IS" (see MIT License disclaimer)
- No warranty of any kind, express or implied
- The developers are not liable for any damages arising from the use of this software
- Use of external data sources (PubChem, ChEMBL, etc.) is subject to their respective terms of service
- Direct use for medical or safety decisions requires expert review

## Documents

| File | Description |
|---|---|
| [docs/quickstart.md](docs/quickstart.md) | Setup steps & FAQ |
| [docs/function_reference.md](docs/function_reference.md) | Function signature reference |
| [docs/python_integration.md](docs/python_integration.md) | Python integration architecture |
| [docs/platform_support.md](docs/platform_support.md) | Desktop / Online platform support |
| [docs/compliance.md](docs/compliance.md) | License & compliance |
