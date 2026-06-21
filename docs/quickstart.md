# Quick Start Guide — EasyMolKit

> Version: v0.2.0 / Updated: 2026-04-19

---

## Prerequisites

| Item | Desktop | MATLAB Online |
|---|---|---|
| MATLAB | R2025b or later | R2025b or later |
| Network | Required for initial setup | Required |
| Python | Not required (auto-deployed) | Pre-installed |
| `config/settings.json` | Copy from `settings.example.json` | Same |

---

## Initial Setup

### MATLAB Desktop

```matlab
% 1. Clone the repository
%    git clone https://github.com/<owner>/EasyMolKit.git

% 2. Open the EasyMolKit folder in MATLAB

% 3. Run setup (one time only — automatically deploys Python + RDKit)
addpath(genpath("src"));
emk.setup.install();
```

> `emk.setup.install()` performs the following automatically:
> 1. Platform detection (Windows Desktop / MATLAB Online)
> 2. Download and extract Embedded Python 3.10 into `python_env/`
> 3. Install RDKit via pip
> 4. Configure `pyenv` (OutOfProcess mode)
> 5. Verify the RDKit import

### MATLAB Online

```matlab
addpath(genpath("src"));
emk.setup.installOnline();
```

> On MATLAB Online, `get-pip.py` bootstrap + `!~/.local/bin/pip install rdkit-pypi` is run
> against the pre-installed Python.

### MCP Server Configuration (developers only)

Open `.vscode/mcp.json` and update the two path entries to the absolute path of your project:

```json
"--initial-working-folder", "D:\\workspace\\EasyMolKit",
"--log-folder", "D:\\workspace\\EasyMolKit\\logs\\matlab-mcp",
```

### Prepare the configuration file

```powershell
copy config\settings.example.json config\settings.json
# Edit settings.json to configure options
```

---

## Basic Usage

### Step 1. Experience it with section execution

```matlab
addpath(genpath("src"));
main_emk    % Use Run Section (Ctrl+Enter) to execute sections one at a time
```

> ⚠️ Use **Run Section (Ctrl+Enter)**, not "Run" (F5).
> Pressing F5 evaluates all sections at once.

### Step 2. Use the API directly

```matlab
% Create a molecule
mol = emk.mol.fromSmiles("CCO");          % Ethanol

% Calculate descriptors
mw   = emk.descriptor.molWeight(mol);     % Molecular weight
logP = emk.descriptor.calculate(mol, "LogP");

% Fingerprints
fp = emk.fingerprint.morgan(mol);

% Similarity calculation
mol2  = emk.mol.fromSmiles("CCCO");
fp2   = emk.fingerprint.morgan(mol2);
score = emk.similarity.tanimoto(fp, fp2);

% SDF file I/O
mols = emk.io.readSdf("data/sample.sdf");
```

### Step 3. Check the output artifacts

```
result/runs/<YYYYMMDD_HHMMSS>/
  ├─ run_meta.json
  └─ *.csv / *.mat
```

---

## Running Tests

```matlab
% Unit tests
addpath(genpath("src"));
suite   = testsuite("tests/unit");
runner  = matlab.unittest.TestRunner.withNoPlugins;
results = runner.run(suite);
fprintf("RESULT: %d PASS / %d FAIL / %d Total\n", ...
    sum([results.Passed]), sum([results.Failed]), numel(results));

% Smoke tests
addpath(genpath("src")); addpath("tests/smoke");
test_mvp_smoke();
```

---

## Three-tier UX

| Tier | Entry point | Audience | Characteristics |
|---|---|---|---|
| **Tier 1: Try it now** | `main_emk.m` | First-time users | Run sections in order with Ctrl+Enter; preset SMILES for an instant experience |
| **Tier 2: Customize** | `main_<feature>.m` | Regular users | Parameter control, section execution |
| **Tier 3: Direct API** | `emk.*` functions | Developers & advanced users | Batch processing, custom scripts |

---

## FAQ

**Q: `emk.setup.install()` fails**
A: Check your network connection. In a proxy environment, set `python.proxy` in `config/settings.json`.

**Q: Path is not recognized**
A: Run `addpath(genpath("src"))`. This is automatically executed in Section 0a of `main_emk.m`.

**Q: Configuration file not found**
A: Copy `config/settings.example.json` to `config/settings.json`.

**Q: Error from `pyenv`**
A: If another Python is already configured in `pyenv`, restart MATLAB and run `emk.setup.initPython()`.
The `pyenv` Version can only be set once per MATLAB session.

**Q: RDKit is not available on MATLAB Online**
A: Run `emk.setup.installOnline()`. It internally runs `get-pip.py` bootstrap + `!~/.local/bin/pip install rdkit-pypi`.

**Q: Cannot connect to the MCP server**
A: Verify that the paths in `.vscode/mcp.json` are correct. MATLAB must be running.

**Q: Output artifacts not found**
A: Check the latest timestamped folder under `result/runs/`.

---

## Track 2: External CPython Environment

Open Babel, MDAnalysis, and PyMOL cannot be bundled with Embedded Python due to GPL licensing
or technical constraints. Install them in a **separate CPython 3.10+ environment** and connect
via `emk.setup.useExternal()`.

> ⚠️ **Important constraint**: `emk.setup.useExternal()` must be called **before Python is loaded**
> in the MATLAB session. Once loaded, it cannot be changed within the same session (MATLAB `pyenv`
> limitation). Call it immediately after `addpath` at the top of your script.

### Step 1: Set up an external Python environment

We recommend separate environments for each library.

**MDAnalysis (MD trajectory analysis)**
```powershell
python -m venv C:\envs\mdenv
C:\envs\mdenv\Scripts\pip install MDAnalysis
```

**Open Babel (chemical file format conversion & 3D coordinate generation)**
```
1. Download the Windows installer (with Python bindings) from
   https://github.com/openbabel/openbabel/releases
2. Install into a CPython 3.10+ environment (python_env/ cannot be used)
3. Follow Step 2 below to connect
```

**PyMOL OSS (3D visualization)**
```powershell
python -m venv C:\envs\pymolenv
C:\envs\pymolenv\Scripts\pip install pymol-open-source
```

### Step 2: Connect at the start of a MATLAB session

```matlab
addpath(genpath("src"));

% Call useExternal() before Python is loaded
emk.setup.useExternal("C:\envs\mdenv\python.exe")

% Use emk.* functions normally from here
mol = emk.mol.fromSmiles("CCO");
```

> To persist the setting across sessions, add it to `config/settings.json`:
> ```json
> { "python": { "external_path": "C:\\envs\\mdenv\\python.exe" } }
> ```

### Step 3: Verify the connection

```matlab
T = emk.setup.validate()
emk.setup.recipe("mdanalysis")   % detailed instructions per library
```

### Available Track 2 libraries

| Library | Purpose | License |
|---|---|---|
| **openbabel** | Chemical file conversion (110+ formats) & 3D coordinate generation | GPLv2 |
| **mdanalysis** | MD trajectory analysis (GROMACS / AMBER / NAMD, etc.) | GPLv2+ |
| **pymol** | 3D molecular visualization (PyMOL Open-Source) | Python-2.0/BSD |

> ⚠️ Open Babel and MDAnalysis are GPLv2 / GPLv2+ libraries. Scripts that use them may be subject
> to GPL conditions. See [docs/compliance.md](compliance.md) before commercial use.

---

## Log Output Format

```
[HH:MM:SS][INFO]  Processing completed (42 molecules)
[HH:MM:SS][WARN]  Invalid SMILES at row 5 — skipped
[HH:MM:SS][ERROR] RDKit import failed: module not found
[####------]  40% ( 4/10) molecules
```
