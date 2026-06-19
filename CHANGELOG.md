# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [1.3.1] - 2026-06-19

### Changed

- **README / README.ja.md**: Rewrote Quick Start to center on `main_rdkit.m` section-by-section execution
  - Now guides users to open `main_rdkit.m` and run Section 0a (path & config) then Section 0b (Python setup) with Ctrl+Enter
  - Removed the standalone `emk.setup.install()` code snippet that bypassed the `cfg.useCase.*` configuration in Section 0a
  - Added explicit warning to use Ctrl+Enter (Run Section), not F5 (Run File)
- **README / README.ja.md**: Added "Corporate / Restricted Network Environments" guidance under Quick Start
  - Proxy configuration via `cfg.python.proxy` in Section 0a of `main_rdkit.m`
  - Windows Defender / Smart App Control may quarantine Embedded Python extraction; whitelist `python_env/` as workaround
  - IT policy blocking executable downloads — MATLAB Online recommended as an alternative
  - Antivirus quarantine of Python binaries — whitelist `python_env/` or use MATLAB Online

---

## [1.3.0] - 2026-06-16

### Added

**Tutorial content**
- **Layer 4: Research (R01–R10)** — 10 research-level tutorials (Japanese & English, 30–90 min each).
  - **R01** Large-scale Similarity Screening (GPU) — GPU-accelerated batch Tanimoto; requires Parallel Computing Toolbox (GPU); CPU fallback on MATLAB Online
  - **R02** PK/PD Simulation — pharmacokinetic/pharmacodynamic modeling; requires SimBiology
  - **R03** Forensic Chemometrics — mixture deconvolution and classification; requires Statistics and ML Toolbox + Parallel Computing Toolbox
  - **R04** Protein-Ligand Analysis † — PDB parsing and binding-site analysis; requires Bioinformatics Toolbox
  - **R05** Molecular Language Model: SMILES Generation — character-level RNN generative model; requires Deep Learning Toolbox
  - **R06** REINFORCE Molecular Design — policy-gradient de novo design; requires Deep Learning Toolbox + Reinforcement Learning Toolbox
  - **R07** Metabolomics † — untargeted metabolomics workflow; requires Bioinformatics Toolbox + SimBiology
  - **R08** Protein-Ligand Docking Simulation ‡ — automated docking via AutoDock Vina; MATLAB Online only (Track 1: meeko + vina + pdbfixer)
  - **R09** GNN Molecular Property Prediction § — message-passing GNN for property regression; requires Deep Learning Toolbox + PyTorch stack
  - **R10** ChemBERTa Transfer Learning § — HuggingFace transformer fine-tuning for SMILES; requires Deep Learning Toolbox + PyTorch stack

  > **†** Requires `emk.setup.installExtra("biopython")` before first run (Track 1 add-on).
  > **‡** MATLAB Online only (Windows Desktop unsupported: no Windows PyPI wheel for vina; pdbfixer's openmm blocked by Smart App Control). Setup: set `cfg.optionalLibraries.meeko/vina/pdbfixer = true` in `main_rdkit.m`, then run `installOnline(Config=cfg)`.
  > **§** Requires PyTorch + HuggingFace stack. Install in order: `torch` → `torch_geometric` → `transformers` → `datasets`. R10 requires R09's torch environment.

**README**
- Added "Open in MATLAB Online" badge to both English and Japanese READMEs — opens the full repository in MATLAB Online with one click

### Changed

- **Analytics (A01–A10) English** — Applied Step 2 English cultural adaptation (idiomatic phrasing, contextual rewrites throughout); fixed `logInfo` argument order

### Removed

- Answer scripts for Analytics modules (`examples/*/analytics/answers/`) — removed from distribution

---

## [1.2.0] - 2026-06-09

### Added

**Tutorial content**
- **Layer 3: Analytics (A01–A10)** — 10 analytical tutorials (Japanese & English, 30–60 min each).
  Requires Statistics and ML Toolbox (A01–A04, A07–A09), Deep Learning Toolbox (A05),
  Curve Fitting Toolbox (A06), and Optimization Toolbox (A09–A10, optional).
  All 10 modules run on MATLAB Online Basic (free tier).
  - **A01** Chemical Space Mapping by PCA — dimensionality reduction, chemical space visualization
  - **A02** Molecular Clustering — hierarchical clustering, structural similarity grouping
  - **A03** QSAR Regression — LogP prediction, regression model evaluation
  - **A04** Drug Classification — SVM / Random Forest, ROC curve analysis
  - **A05** Neural Network Property Prediction — feedforward NN, property regression
  - **A06** Dose-Response Curve Fitting — Hill equation, EC50 estimation
  - **A07** Scaffold Analysis & R-Group Decomposition — medicinal chemistry analytics
  - **A08** Mass Spectrometry × Cheminformatics — isotope pattern matching, MS annotation
  - **A09** PFAS & Environmental Screening — SMARTS-based screening, Pareto optimization
  - **A10** Lead Optimization — multi-objective optimization, Derringer-Suich desirability
  - Answer scripts for all Analytics modules (`examples/*/analytics/answers/`)

---

## [1.1.0] - 2026-05-31

### Added

**Tutorial content**
- **Layer 2: Application Stories (S01–S07)** — 7 story-based tutorials (Japanese & English, 20–40 min each).
  All modules require Base MATLAB only (no additional Toolboxes).
  - **S01** Finding Caffeine's Molecular Cousins — similarity search, top-N ranking
  - **S02** Drug-likeness Filter — Lipinski's Rule of Five
  - **S03** Structural Alert Detection for Safety Screening — SMARTS batch matching
  - **S04** Virtual Screening Workflow — fingerprint similarity, Top-N filtering
  - **S05** Unknown Compound Identification Challenge — forensic chemistry scenario
  - **S06** Search Compounds on PubChem — PUG REST API integration
  - **S07** Analyze ChEMBL Activity Data — bioactivity retrieval and analysis
  - Answer scripts for all Stories modules (`examples/*/stories/answers/`)

**Core API**
- `resolveProjectRoot()` (`src/util/resolveProjectRoot.m`) — Three-stage fallback for stable
  project root detection across Windows Desktop, MATLAB MCP, and MATLAB Online environments.

### Changed

- **Foundation (F01–F06) Section 0** — Replaced per-file inline root-detection block (~25 lines)
  with a concise 9-line template using `resolveProjectRoot()`. MATLAB Online path resolution is now stable.

---

## [1.0.0] - 2026-05-25

Initial public release. Includes the **Layer 1 Foundation** tutorial series (F01–F06)
and the complete `emk.*` API for RDKit-based Chemoinformatics from MATLAB.

### Added

**Core API (`src/+emk/`)**
- `emk.setup.install()` — One-command Embedded Python 3.10 + RDKit deployment for Windows Desktop
- `emk.setup.installOnline()` — Automated RDKit setup for MATLAB Online via `get-pip.py` bootstrap
- `emk.setup.initPython()` — Platform-aware `pyenv` configuration (OutOfProcess mode)
- `emk.setup.verify()` — Non-throwing Python/RDKit status check
- `emk.setup.installExtra()` — Track 1 optional library installer (pubchempy, mordred, biopython, torch, etc.)
- `emk.setup.useExternal()` — Track 2 external CPython connector for GPL libraries
- `emk.setup.validate()` — Installed library diagnostics table
- `emk.setup.recipe()` — Per-library installation guide display
- `emk.setup.installTrack2()` — Automated venv creation for Track 2 libraries
- `emk.mol.fromSmiles()` — SMILES → RDKit Mol object
- `emk.mol.toSmiles()` — Canonical SMILES from Mol
- `emk.mol.isValid()` — SMILES validity check (non-throwing)
- `emk.mol.hasSubstruct()` — SMARTS substructure matching
- `emk.mol.toStruct()` / `fromStruct()` — Mol serialization (molblock / pickle)
- `emk.mol.toTable()` — Mol array → MATLAB table with descriptors
- `emk.mol.scaffold()` — Bemis-Murcko scaffold extraction
- `emk.descriptor.molWeight()` — Average molecular weight
- `emk.descriptor.calculate()` — 10 standard descriptors (MolWt, LogP, TPSA, HBD/HBA, etc.)
- `emk.descriptor.batchCalculate()` — Batch descriptor calculation → table
- `emk.descriptor.mordred()` / `mordredBatch()` / `mordredNames()` — Mordred 2D descriptors (~1800)
- `emk.fingerprint.morgan()` — Morgan (ECFP) fingerprint (Radius=2, NBits=2048 default)
- `emk.fingerprint.maccs()` — 167-bit MACCS keys fingerprint
- `emk.fingerprint.toArray()` — Fingerprint → MATLAB logical array
- `emk.similarity.tanimoto()` / `dice()` — Pairwise similarity coefficients
- `emk.similarity.rankBy()` — Top-N similarity ranking against a database
- `emk.similarity.matrix()` — N×N symmetric similarity matrix
- `emk.filter.lipinski()` — Lipinski Rule of Five filter with violation count
- `emk.io.readSdf()` / `writeSdf()` — SDF file read/write
- `emk.io.readSmilesList()` — One-SMILES-per-line file reader
- `emk.viz.draw2d()` — 2D structure rendering via RDKit PNG → MATLAB figure
- `emk.db.searchPubchem()` / `pubchemFetch()` — PubChem compound search (REST + PubChemPy)
- `emk.db.searchChembl()` / `searchChemblTarget()` / `getChemblActivity()` — ChEMBL REST search
- `emk.util.isOnline()` — MATLAB Online detection
- `emk.util.benchmarkBatch()` — Throughput measurement utility

**Tutorial content (`examples/`)**
- **F01** Drawing Molecules with SMILES — molecular representation, SMILES syntax
- **F02** Calculating Molecular Properties — MW / LogP / TPSA
- **F03** Introduction to Fingerprints — bit vectors, Morgan vs MACCS
- **F04** Comparing Molecules by Similarity — Tanimoto / Dice
- **F05** Substructure Search — SMARTS pattern matching
- **F06** Reading Molecules from Files — SDF / SMILES file I/O
- Answer scripts for all Foundation modules (`examples/*/foundation/answers/`)
- Both Japanese (`examples/japanese/`) and English (`examples/english/`) versions provided

**Documentation**
- `docs/quickstart.md` — Setup steps & FAQ
- `docs/function_reference.md` — Full API signature reference
- `docs/algorithm_guide.md` — Algorithm rationale & test strategy
- `docs/python_integration.md` — Python integration architecture
- `docs/platform_support.md` — Desktop / Online platform details
- `docs/compliance.md` — License & compliance notes

**Infrastructure**
- `main_rdkit.m` — Section-executable entry point (Ctrl+Enter workflow)
- `config/settings.example.json` — Configuration template
- `data/list/` — Curated sample datasets (CC0 / BSD-3 / CC-BY-SA 3.0)
- `tests/unit/` — `matlab.unittest` class-based test suite
- `tests/smoke/` — Smoke test suite

### Notes

- Requires MATLAB R2025b or later
- Windows Desktop and MATLAB Online supported; macOS / Linux Desktop untested
- Layer 2 (Stories) was released in v1.1.0; Layer 3 (Analytics) in v1.2.0
- Layer 4 (Research) tutorials were released in v1.3.0
