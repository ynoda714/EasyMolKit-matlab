%% EasyMolKit -- Optional Library Setup
% Use this file to install and manage Track 1 / Track 2 optional libraries.
% Run each section with Ctrl+Enter (MATLAB Run Section).
% Always run Section 0 first to initialize paths.
%
% Track 1 (Embedded Python, Desktop only):
%   Install via emk.setup.installExtra(name).
%   These packages go into python_env/ alongside RDKit.
%   Library name -> license -> approximate download size:
%     pubchempy       MIT              ~1 MB
%     mordred         BSD-3            ~5 MB
%     biopython       Biopython Lic.  ~50 MB
%     torch           BSD-3          ~800 MB  (CPU-only)
%     torch_geometric MIT            ~200 MB  (requires torch first)
%     transformers    Apache-2.0     ~200 MB
%     datasets        Apache-2.0      ~50 MB  (pairs with transformers)
%     meeko           LGPL-2.1        ~5 MB   (see docs/compliance.md CL-7)
%     vina            Apache-2.0      ~30 MB
%     pdbfixer        MIT             ~70 MB  (auto-installs openmm)
%
% Track 2 (separate venv, Desktop only):
%   Install via emk.setup.installTrack2(name).
%   Creates python_env_t2/<name>/ and updates settings.json.
%   Library name -> license:
%     mdanalysis      GPL-2.0
%     pymol           PSF (open-source build)
%
% Usage:
%   1. Run Section 0 to set up paths (required before any other section)
%   2. Run the section for the library you want to install
%   3. Each section is idempotent -- safe to re-run
%
% Notes:
%   - Desktop only. MATLAB Online does not support installExtra or installTrack2.
%   - For install recipes (manual steps, notes), run Section R.
%   - torch must be installed BEFORE torch_geometric or transformers.

%% Section 0: Path Init  [RUN THIS FIRST]
% ---------------------------------------------------------
addpath(genpath("src"));
emk.setup.initPython();
logInfo("main_setup_extra: paths initialized, Python ready");

%% Section R: Show Install Recipe (reference only, no install)
% ---------------------------------------------------------
% Displays install commands, license notes, and usage hints.
% Uncomment the library you want to inspect:

emk.setup.recipe("pubchempy")
emk.setup.recipe("mordred")
emk.setup.recipe("biopython")
emk.setup.recipe("torch")
emk.setup.recipe("torch_geometric")
emk.setup.recipe("transformers")
emk.setup.recipe("datasets")
emk.setup.recipe("meeko")      % LGPL-2.1 -- see docs/compliance.md CL-7
emk.setup.recipe("vina")
emk.setup.recipe("pdbfixer")
emk.setup.recipe("docking")    % combined: meeko + vina + pdbfixer
emk.setup.recipe("openbabel")
emk.setup.recipe("mdanalysis")
emk.setup.recipe("pymol")

%% Section 1a: Track 1 -- PubChem / Descriptors / Bio
% ---------------------------------------------------------
emk.setup.installExtra("pubchempy")
emk.setup.installExtra("mordred")
emk.setup.installExtra("biopython")

%% Section 1b: Track 1 -- PyTorch (CPU-only, ~800 MB)
% ---------------------------------------------------------
% Install torch before torch_geometric or transformers.
emk.setup.installExtra("torch")

%% Section 1c: Track 1 -- PyTorch Geometric / Transformers
% ---------------------------------------------------------
% Requires torch to be installed first (Section 1b).
emk.setup.installExtra("torch_geometric")
emk.setup.installExtra("transformers")
emk.setup.installExtra("datasets")

%% Section 1d: Track 1 -- Docking Pipeline (MATLAB Online only)
% ---------------------------------------------------------
% meeko / vina / pdbfixer are NOT installed from this Desktop-only file.
%
%   Reason:
%     vina   -- no Windows PyPI wheel (Boost C++ build required)
%     pdbfixer -- requires openmm; blocked by Windows Smart App Control
%     meeko  -- ligand prep tool; requires vina to be useful
%
%   To install on MATLAB Online:
%     In main_rdkit.m Section 0a, set the flag below to true,
%     then run Section 0b:
%       cfg.useCase.docking = true
%
%   For install notes and standalone binary links, see Section R:
%     emk.setup.recipe("docking")

%% Section 2a: Track 2 -- MDAnalysis (GPL-2.0, separate venv)
% ---------------------------------------------------------
% Creates python_env_t2/mdanalysis/ and updates settings.json.
% emk.setup.installTrack2("mdanalysis")

%% Section 2b: Track 2 -- PyMOL open-source (separate venv)
% ---------------------------------------------------------
% Creates python_env_t2/pymol/ and updates settings.json.
% emk.setup.installTrack2("pymol")

%% Section 3: Verify Installed Libraries
% ---------------------------------------------------------
% Calls verify() and shows which optional libraries are importable.
result = emk.setup.verify();
logInfo("Python loaded : %d", result.python);
logInfo("RDKit loaded  : %d", result.rdkit);
logInfo("Python version: %s", result.version);
