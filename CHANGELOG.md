# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

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
- Layer 2 (Stories), Layer 3 (Analytics), and Layer 4 (Research) tutorials are planned for
  v2.0.0, v3.0.0, and v4.0.0 respectively
