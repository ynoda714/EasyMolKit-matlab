# EasyMolKit

🇯🇵 [日本語はこちら](docs/ja/README.ja.md)

**An integrated Chemoinformatics environment that makes RDKit easy to use from MATLAB**

> Clone the repo → run one setup command → RDKit is ready to use in MATLAB

## Why EasyMolKit?

In the Chemoinformatics field, Python + RDKit is the de-facto standard toolchain.
However, getting started involves several barriers:

- Managing Python versions and virtual environments
- Complex RDKit installation via conda/pip
- Environment conflicts with commercial tools (e.g., PyMOL)
- The Python ecosystem is unfamiliar to MATLAB users

**EasyMolKit removes these barriers.** Leveraging MATLAB's `pyenv`-based Python integration,
users can access RDKit functionality as standard MATLAB functions — no Python knowledge required.

## Features

- **Zero configuration**: One call to `emk.setup.install()` automatically deploys Python + RDKit
- **MATLAB native**: Results are returned as MATLAB `table` / `struct` / `double` — immediately usable in your workspace
- **Desktop & Online**: Supports Windows Desktop and MATLAB Online (macOS / Linux Desktop untested)
- **RDKit wrappers**: Molecular analysis, descriptor calculation, fingerprints, and similarity search — all as MATLAB functions
- **Future-ready**: 3D structure visualization via PyMOL Open-Source planned for a future release

## Target Users

- **Chemistry, pharmacy, and medical researchers** who use MATLAB as their primary research environment
- **Students** learning Chemoinformatics (MATLAB Online Basic free tier covers Layers 1–3)
- **MATLAB users** who want to avoid spending time on Python environment setup

## Requirements

| Item | Desktop | MATLAB Online |
|---|---|---|
| MATLAB | R2025b or later | R2025b or later |
| Python | Auto-deployed (Embedded Python) | Pre-installed |
| RDKit | Auto-deployed | Installed via `emk.setup.installOnline()` |
| OS | Windows | — |

## Quick Start

```matlab
% 1. Clone the repository
%    git clone https://github.com/ynoda/EasyMolKit.git
%    cd EasyMolKit

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

### Track 1: Install additional packages into Embedded Python

```matlab
% Review installation steps and license before installing
emk.setup.recipe("pubchempy")          % Show installation recipe and license
emk.setup.installExtra("pubchempy")    % Install into Embedded Python
emk.setup.installExtra("mordred")      % 1800+ descriptor library
emk.setup.installExtra("biopython")    % PDB / sequence analysis

% PyTorch + HuggingFace stack (required for R09 / R10)
emk.setup.installExtra("torch")           % CPU-only, ~800 MB (must be installed first)
emk.setup.installExtra("torch_geometric") % GNN library (requires torch)
emk.setup.installExtra("transformers")    % HuggingFace Transformers
emk.setup.installExtra("datasets")        % HuggingFace Datasets (used with transformers)

% Verify installation
T = emk.setup.validate()
```

> For bulk installation of all libraries, see `main_setup_extra.m`.

### Track 2: Libraries that require an external CPython environment

Open Babel, MDAnalysis, and PyMOL cannot be bundled with Embedded Python due to GPL licensing
or technical constraints. Install them in a **separate CPython 3.10+ environment** and connect
via `emk.setup.useExternal()`.

> ⚠️ **Important constraint**: `emk.setup.useExternal()` must be called **before Python is loaded**
> in the MATLAB session. Once loaded, it cannot be changed within the same session (MATLAB `pyenv`
> limitation). When using Track 2, call it immediately after `addpath` at the top of your script.

#### Step 1: Set up an external Python environment

We recommend separate environments for each library. Examples:

**MDAnalysis (MD trajectory analysis)**
```powershell
# Run in PowerShell / Command Prompt
python -m venv C:\envs\mdenv
C:\envs\mdenv\Scripts\pip install MDAnalysis
```

**Open Babel (chemical file format conversion & 3D coordinate generation)**
```
1. Download the Windows installer (with Python bindings) from
   https://github.com/openbabel/openbabel/releases
2. Install into a CPython 3.10+ environment (python_env/ cannot be used)
3. Follow the same steps below to call useExternal()
```

**PyMOL OSS (3D visualization)**
```powershell
python -m venv C:\envs\pymolenv
C:\envs\pymolenv\Scripts\pip install pymol-open-source
```

#### Step 2: Connect at the start of a MATLAB session

```matlab
addpath(genpath("src"));

% Call useExternal() here — before Python is loaded
emk.setup.useExternal("C:\envs\mdenv\python.exe")

% Use emk.* functions normally from here
mol = emk.mol.fromSmiles("CCO");
```

> To avoid calling this every session, add it to `config/settings.json` for automatic loading:
>
> ```json
> {
>   "python": {
>     "external_path": "C:\\envs\\mdenv\\python.exe"
>   }
> }
> ```

#### Step 3: Verify the connection

```matlab
% Confirm the library is accessible from the connected environment
T = emk.setup.validate()

% Review detailed instructions for each library
emk.setup.recipe("openbabel")
emk.setup.recipe("mdanalysis")
emk.setup.recipe("pymol")
```

#### Track 2 available libraries

| Library | Purpose | pip command | License |
|---|---|---|---|
| **openbabel** | Chemical file conversion (110+ formats) & 3D coordinate generation | Windows installer required (see above) | GPLv2 |
| **mdanalysis** | MD trajectory analysis (GROMACS / AMBER / NAMD, etc.) | `pip install MDAnalysis` | GPLv2+ |
| **pymol** | 3D molecular visualization (PyMOL Open-Source) | `pip install pymol-open-source` | Python-2.0/BSD |

> ⚠️ **GPL license notice**: Open Babel and MDAnalysis are GPLv2 / GPLv2+ libraries.
> EasyMolKit itself (MIT) is not affected, but scripts that use these libraries may be subject
> to GPL conditions. Review [docs/compliance.md](docs/compliance.md) before commercial use.

## Tutorials & Examples (4-Layer Structure)

EasyMolKit provides progressive learning content under `examples/`.

| Layer | Audience | Required Toolbox | Content | Release |
|---|---|---|---|---|
| **L1 Foundation** | All users | None | Learn one API concept at a time (6 modules, 5–15 min each) | ✅ v1.0.0 |
| **L2 Application Stories** | After Foundation | None | Practical workflows combining multiple features (7 modules, 20–40 min each) | ✅ v1.1.0 |
| **L3 Analytics** | All users | Varies (Statistics and ML, etc.) | QSAR, clustering, MS analysis, optimization, and more (A01–A10, 10 modules, 30–60 min each) | 🔜 v1.2.0 |
| **L4 Research** | All users | Varies (Parallel Computing, etc.) | Research-level applications (R01–R10, 30–90 min each) | 🔜 v1.3.0 |

*L1–L3 run entirely on MATLAB Online Basic (free tier).*

### Layer 1: Foundation (Base MATLAB only)

| # | Title | Required Toolbox | What You Learn | Desktop | Online |
|---|---|---|---|:---:|:---:|
| F01 | Drawing Molecules with SMILES | None | Molecular representation, SMILES syntax | ✔ | ✔ |
| F02 | Calculating Molecular Properties | None | MW / LogP / TPSA — meaning and calculation | ✔ | ✔ |
| F03 | Introduction to Fingerprints | None | Bit vector representation, Morgan vs MACCS | ✔ | ✔ |
| F04 | Comparing Molecules by Similarity | None | Quantifying Tanimoto / Dice similarity | ✔ | ✔ |
| F05 | Substructure Search | None | SMARTS pattern matching | ✔ | ✔ |
| F06 | Reading Molecules from Files | None | SDF / SMILES file operations | ✔ | ✔ |

### Layer 2: Application Stories (Base MATLAB only)

| # | Title | Required Toolbox | Domain | Desktop | Online |
|---|---|---|---|:---:|:---:|
| S01 | Find Relatives of Caffeine | None | Everyday chemistry | ✔ | ✔ |
| S02 | Drug Filters: Lipinski's Rule of Five | None | Pharmacology | ✔ | ✔ |
| S03 | Structural Alerts for Hazardous Compounds | None | Safety | ✔ | ✔ |
| S04 | Introduction to Virtual Screening | None | Drug discovery | ✔ | ✔ |
| S05 | Unknown Compound Identification Challenge | None | Forensics | ✔ | ✔ |
| S06 | Search Compounds on PubChem | None | Databases | ✔ | ✔ |
| S07 | Analyze ChEMBL Activity Data | None | Drug discovery | ✔ | ✔ |

### Layer 3 & 4: Analytics and Research

> 🔜 **Not yet released.** Detailed module lists will be published when available.
>
> L3 Analytics (A01–A10) is planned for **v1.2.0**. L4 Research (R01–R10) is planned for **v1.3.0**.

## Directory Structure

```
EasyMolKit/
├─ main_rdkit.m               # RDKit setup & basic operations (run section by section)
├─ config/
│   └─ settings.example.json  # Configuration template
├─ examples/
│   ├─ japanese/              # Distribution materials — Japanese (plain-text Live Code)
│   │   ├─ foundation/        #   L1: API basics (Base MATLAB only)
│   │   ├─ stories/           #   L2: Application stories (Base MATLAB only)
│   │   ├─ analytics/         #   L3: Statistics & ML integration (free toolboxes)
│   │   └─ research/          #   L4: Research level (Campus License)
│   └─ english/               # Distribution materials — English (comments differ only)
│       ├─ foundation/
│       ├─ stories/
│       ├─ analytics/
│       └─ research/
├─ src/
│   ├─ +emk/                  # Main package
│   │   ├─ +setup/            # Python environment setup
│   │   ├─ +mol/              # Molecular object operations
│   │   ├─ +descriptor/       # Descriptor calculation
│   │   ├─ +fingerprint/      # Fingerprint generation
│   │   ├─ +similarity/       # Similarity calculation
│   │   ├─ +filter/           # Molecular filtering (Lipinski, etc.)
│   │   ├─ +io/               # File I/O
│   │   ├─ +viz/              # Visualization (future PyMOL integration)
│   │   └─ +util/             # Package-level utilities
│   ├─ config/                # Configuration loader
│   └─ util/                  # Log helpers & common utilities
├─ result/                    # Run artifacts (not tracked by Git)
├─ tests/
│   ├─ unit/                  # matlab.unittest class-based tests
│   └─ smoke/                 # Smoke tests
├─ data/                      # Curated sample data
└─ docs/                      # Documentation
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
| `emk.viz` | `draw2d()` | 2D structure rendering (※) |

> **※ Rendering note**: `emk.viz.draw2d()` generates a PNG via RDKit (Python) and transfers it to MATLAB.
> Rendering takes **0.5–2 seconds per molecule**. On MATLAB Online, inter-process communication overhead
> adds further latency when rendering many molecules in sequence (this is a structural constraint and
> cannot be improved).

For full API details, see [docs/function_reference.md](docs/function_reference.md).

## Key Conventions

- All logic lives under `src/` — no `.m` files in the project root (except `main_<feature>.m` entry points)
- `.m` files use English only (comments, logs, error messages)
- Prefer `string` type (`"..."` literals); requires R2025b+
- Never call `py.rdkit.*` directly — always go through `emk.*` wrappers
- Use `logInfo` / `logWarn` / `logError` for all output (`fprintf` is prohibited)
- Run artifacts are saved to `result/runs/<YYYYMMDD_HHMMSS>/` (not tracked by Git)

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

This software is provided for research and educational purposes.

- This software is provided "AS IS" without warranty of any kind, express or implied
- The developers are not liable for any damages arising from the use of this software
- Use of external data sources (PubChem, ChEMBL, etc.) is subject to their respective terms of service
- Predictions and calculation results are for research purposes only; direct application to medical or safety decisions requires expert review

## Documents

| File | Description |
|---|---|
| [docs/quickstart.md](docs/quickstart.md) | Setup steps & FAQ |
| [docs/function_reference.md](docs/function_reference.md) | Full function signature reference |
| [docs/function_catalog.md](docs/function_catalog.md) | Compact function catalog (51 functions, all modules) |
| [docs/test_catalog.md](docs/test_catalog.md) | Test class catalog (25 classes, unit + smoke) |
| [docs/python_integration.md](docs/python_integration.md) | Python integration architecture |
| [docs/platform_support.md](docs/platform_support.md) | Desktop / Online platform support |
| [docs/compliance.md](docs/compliance.md) | License & compliance |

### Japanese Documents

| File | Description |
|---|---|
| [docs/ja/README.ja.md](docs/ja/README.ja.md) | このリポジトリの概要（日本語版） |
| [docs/function_catalog.ja.md](docs/function_catalog.ja.md) | コンパクト関数カタログ（日本語版） |
| [docs/test_catalog.ja.md](docs/test_catalog.ja.md) | テストクラスカタログ（日本語版） |
