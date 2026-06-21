# EasyMolKit

🇯🇵 [日本語はこちら](docs/ja/README.ja.md)

**An integrated Chemoinformatics environment that makes RDKit easy to use from MATLAB**

> Clone the repo → run one setup command → RDKit is ready to use in MATLAB

[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=ynoda714/EasyMolKit-matlab)

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
- **Rich API**: 76 functions across 15 modules — descriptors, fingerprints, scaffolds, filters, clustering, 3D conformers, and more
- **Reproducible research**: 6 published papers reproduced with locked environments under `repro/`

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

% 2. Open main_rdkit.m in MATLAB, then run each section with Ctrl+Enter:
%
%   Section 0a  — Path setup & config         (edit cfg.useCase.* here if needed)
%   Section 0b  — Python + RDKit setup        (first time only; ~2-5 min)
%   Section 1   — Basic molecule operations
%   Section 2   — Descriptor calculation
%   Section 3   — Fingerprints & similarity
```

> ⚠️ Use **Ctrl+Enter (Run Section)**, not **F5 (Run File)**.
> Running all sections at once will fail on first setup.

For more details, see [docs/quickstart.md](docs/quickstart.md).

### Corporate / Restricted Network Environments

If you are on a **corporate PC**, local Python deployment may be blocked by IT policy:

| Issue | Symptom | Solution |
|---|---|---|
| Proxy server | `pip install` times out / SSL error | Set `cfg.python.proxy = "http://proxy.example.com:8080"` in Section 0a of `main_rdkit.m` |
| Windows Defender / Smart App Control | Embedded Python extraction is quarantined | Whitelist the `python_env/` directory, or use MATLAB Online |
| IT policy (executable downloads blocked) | Setup fails at the download step | Use **MATLAB Online** — no local Python deployment needed |
| Antivirus quarantine | Python binaries disappear after extraction | Whitelist `python_env/`, or use MATLAB Online |

> 💡 **Recommended for corporate environments**: Use [MATLAB Online](https://matlab.mathworks.com/open/github/v1?repo=ynoda714/EasyMolKit-matlab) — no local Python installation is needed, and all L1–L3 tutorials run on the free Basic tier.

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

Open Babel, MDAnalysis, and PyMOL require a **separate CPython 3.10+ environment** due to GPL
licensing or technical constraints. Connect via `emk.setup.useExternal()` — which must be called
**before Python is loaded** in the MATLAB session.

For step-by-step setup instructions, see [docs/quickstart.md — Track 2](docs/quickstart.md).

## Tutorials & Examples

EasyMolKit provides progressive learning content under `examples/`.

| Layer | Audience | Content | Release |
|---|---|---|---|
| **L1 Foundation** | All users | One API concept at a time (6 modules, 5–15 min each) | ✅ v1.0.0 |
| **L2 Application Stories** | After Foundation | Practical workflows combining multiple features (7 modules, 20–40 min each) | ✅ v1.1.0 |
| **L3 Analytics** | All users | QSAR, clustering, MS analysis, optimization (A01–A10, 30–60 min each) | ✅ v1.2.0 |
| **L4 Research** | All users | Research-level applications (R01–R10, 30–90 min each) | ✅ v1.3.0 |

*L1–L3 run entirely on MATLAB Online Basic (free tier).*

For the full per-module listing with Toolbox requirements and platform support,
see [docs/tutorials.md](docs/tutorials.md).

## Reproducible Research

`repro/` contains MATLAB reproductions of published Chemoinformatics papers, each with a
locked environment snapshot and defined success criteria.

| ID | Paper | Method | Result |
|---|---|---|---|
| RP00 | Delaney (2004) ESOL | Linear regression on physicochemical descriptors | CV RMSE=1.017, R²=0.762 |
| RP01 | Delaney (2004) ESOL extended | Linear regression + TPSA / QED / SA Score | CV RMSE=0.584, R²=0.906 |
| RP02 | Wu et al. (2018) MoleculeNet BBBP | Morgan FP (ECFP4) + Random Forest | ROC-AUC CV=0.883 |
| RP03 | Yang et al. (2019) GNN on BBBP | Graph Convolutional Network | ROC-AUC CV=0.915 |
| RP04 | Chithrananda et al. (2020) ChemBERTa | Frozen CLS embedding + Logistic Regression | ROC-AUC CV=0.927 |
| RP05 | SHAP explainability on BBBP | shap.LinearExplainer + LR model | ROC-AUC CV=0.909, Spearman ρ=0.902 |

## API Overview

| Module | Example functions | Description |
|---|---|---|
| `emk.setup` | `install()`, `verify()`, `snapshot()`, `verifyLock()` | Python environment deployment, initialization & RF02 version lock |
| `emk.mol` | `fromSmiles()`, `toSmiles()`, `isValid()`, `hasSubstruct()` | Molecular object creation & conversion |
| `emk.descriptor` | `molWeight()`, `calculate()`, `qed()`, `saScore()`, `bcut()` | Molecular descriptor calculation |
| `emk.fingerprint` | `morgan()`, `maccs()`, `toArray()` | Fingerprint generation |
| `emk.similarity` | `tanimoto()`, `dice()`, `rankBy()`, `matrix()` | Molecular similarity calculation |
| `emk.scaffold` | `genericMurcko()`, `brics()`, `rgroup()` | Scaffold analysis & fragment decomposition |
| `emk.dataset` | `esol()`, `freesolv()`, `bbbp()`, `tox21()` | Benchmark dataset loaders with local cache |
| `emk.filter` | `lipinski()`, `veber()`, `pains()`, `reos()` | Medicinal chemistry filters |
| `emk.cluster` | `butina()` | Butina sphere-exclusion clustering |
| `emk.diversity` | `pick()` | MaxMin diverse subset selection |
| `emk.conformer` | `embed()`, `optimize()` | 3D conformer generation & force-field optimization |
| `emk.shape` | `compare()` | 3D shape similarity comparison |
| `emk.repro` | `verify()` | RF03 reproduction success verification |
| `emk.io` | `readSdf()`, `writeSdf()`, `readSmilesList()` | SDF / SMILES file I/O |
| `emk.viz` | `draw2d()` | 2D structure rendering (※) |

> **※ Rendering note**: `emk.viz.draw2d()` generates a PNG via RDKit (Python) and transfers it to MATLAB.
> Rendering takes **0.5–2 seconds per molecule**. On MATLAB Online, inter-process communication overhead
> adds further latency when rendering many molecules in sequence (this is a structural constraint and
> cannot be improved).

For full API details, see [docs/function_reference.md](docs/function_reference.md).

## Directory Structure

```
EasyMolKit/
├─ main_rdkit.m               # RDKit setup & basic operations (run section by section)
├─ config/
│   └─ settings.example.json  # Configuration template
├─ examples/
│   ├─ japanese/              # Distribution materials — Japanese (plain-text Live Code)
│   └─ english/               # Distribution materials — English (comments differ only)
├─ repro/                     # Reproducible research (RP00–RP05)
├─ src/
│   └─ +emk/                  # Main package (15 modules, 76 functions)
├─ tests/
│   ├─ unit/                  # matlab.unittest class-based tests
│   └─ smoke/                 # Smoke tests
├─ data/                      # Curated sample data
└─ docs/                      # Documentation
```

## License

EasyMolKit: [MIT License](LICENSE)

### Third-party licenses

| Library | License | Purpose |
|---|---|---|
| RDKit | BSD-3-Clause | Chemoinformatics core |
| Python (CPython) | PSF License | Runtime environment |

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
| [docs/quickstart.md](docs/quickstart.md) | Setup steps, Track 2 setup & FAQ |
| [docs/tutorials.md](docs/tutorials.md) | Full tutorial listing (F01–R10, RP00–RP05) |
| [docs/function_reference.md](docs/function_reference.md) | Full function signature reference |
| [docs/function_catalog.md](docs/function_catalog.md) | Compact function catalog (76 functions) |
| [docs/python_integration.md](docs/python_integration.md) | Python integration architecture |
| [docs/platform_support.md](docs/platform_support.md) | Desktop / Online platform support |
| [docs/compliance.md](docs/compliance.md) | License & compliance |

### Japanese Documents

| File | Description |
|---|---|
| [docs/ja/README.ja.md](docs/ja/README.ja.md) | このリポジトリの概要（日本語版） |
| [docs/ja/tutorials.ja.md](docs/ja/tutorials.ja.md) | チュートリアル一覧（日本語版） |
| [docs/function_catalog.ja.md](docs/function_catalog.ja.md) | コンパクト関数カタログ（日本語版） |
