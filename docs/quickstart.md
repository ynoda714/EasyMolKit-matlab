# Quick Start Guide — EasyMolKit

🇯🇵 [日本語](ja/quickstart_ja.md)

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
main_rdkit    % Use Run Section (Ctrl+Enter) to execute sections one at a time
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

> Full test class list → [test_catalog.md](test_catalog.md)

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
| **Tier 1: Try it now** | `main_rdkit.m` | First-time users | Run sections in order with Ctrl+Enter; preset SMILES for an instant experience |
| **Tier 2: Customize** | `main_<feature>.m` | Regular users | Parameter control, section execution |
| **Tier 3: Direct API** | `emk.*` functions | Developers & advanced users | Batch processing, custom scripts |

---

## FAQ

**Q: `emk.setup.install()` fails**
A: Check your network connection. In a proxy environment, set `python.proxy` in `config/settings.json`.

**Q: Path is not recognized**
A: Run `addpath(genpath("src"))`. This is automatically executed in Section 0a of `main_rdkit.m`.

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

## Log Output Format

```
[HH:MM:SS][INFO]  Processing completed (42 molecules)
[HH:MM:SS][WARN]  Invalid SMILES at row 5 — skipped
[HH:MM:SS][ERROR] RDKit import failed: module not found
[####------]  40% ( 4/10) molecules
```
