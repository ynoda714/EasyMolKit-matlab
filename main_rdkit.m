%% EasyMolKit -- Main Entry Point
% Run each section with Ctrl+Enter (MATLAB Run Section).
% Always run Section 0a first to initialize paths and configuration.
%
% Usage:
%   1. Open this file in MATLAB
%   2. Edit Section 0a: toggle optional libraries on/off
%   3. Run Section 0a (Ctrl+Enter) -- path setup & config
%   4. Run Section 0b (Ctrl+Enter) -- Python setup (first time only)
%   5. Run subsequent sections as needed

%% Section 0a: Path Setup & User Configuration  [RUN THIS FIRST]
% ---------------------------------------------------------
% Edit the "User Settings" block below to choose optional libraries
% and adjust settings.  No settings.json file is required.

addpath(genpath("src"));
cfg = loadConfig();   % loads defaults (+ settings.json if present)

% ---- User Settings: edit this block --------------------------------
%
%   Use-case groups -- set true to install all libraries for that workflow.
%   Install order: torch is resolved before torch_geometric automatically.
%
%   qsar    : pubchempy + mordred               (~6 MB)
%   bio     : biopython                          (~50 MB)
%   ml      : torch (CPU-only) + torch_geometric (~800 MB)  -- R08/R09
%   nlp     : transformers + datasets            (~2 GB)    -- R09
%   docking : scipy + meeko + prody + vina + pdbfixer       -- Online only
%             (Desktop: vina/pdbfixer skipped; see emk.setup.recipe("docking"))
%
cfg.useCase.qsar    = true;   % pubchempy + mordred
cfg.useCase.bio     = true;   % biopython
cfg.useCase.ml      = true;   % torch (CPU-only) + torch_geometric  -- R08/R09
cfg.useCase.nlp     = true;  % transformers + datasets              -- R09
cfg.useCase.docking = true;   % scipy + meeko + prody + vina + pdbfixer (Online only)
%
%   Individual overrides (optional -- prefer useCase above for group installs):
% cfg.optionalLibraries.scipy = true;   % standalone scipy without full docking stack
%
%   Proxy (uncomment and fill in if behind a corporate/university proxy):
% cfg.python.proxy = "http://proxy.example.com:8080";
%
%   RDKit version pin (change only if you need a specific version):
% cfg.rdkit.version = "2024.03.6";
%
% --------------------------------------------------------------------

logInfo("EasyMolKit initialized");
logInfo("  Platform : %s (Online=%d)", computer("arch"), emk.util.isOnline());
logInfo("  Python   : %s (%s)", cfg.python.version, cfg.python.embedded_dir);
logInfo("  RDKit    : %s", cfg.rdkit.version);
ucNames_ = fieldnames(cfg.useCase);
ucParts_ = cellfun(@(f) sprintf('%s=%d', f, cfg.useCase.(f)), ucNames_, "UniformOutput", false);
logInfo("  UseCases : %s", strjoin(ucParts_, "  "));
clear ucNames_ ucParts_;
logInfo("  Mode     : %s", cfg.runtime.eval_mode);

%% Section 0b: Python Environment Setup  [FIRST TIME / AFTER CLEAN]
% ---------------------------------------------------------
% Desktop: downloads Embedded Python 3.10 + RDKit + selected optional libs
% Online : pip-installs rdkit-pypi into user site-packages

if emk.util.isOnline()
    emk.setup.installOnline(Config=cfg);
else
    emk.setup.install(Config=cfg);   % <-- passes your settings above
end

%% Section 1: Basic Molecule Operations
% ---------------------------------------------------------

mol    = emk.mol.fromSmiles("CCO");
smiles = emk.mol.toSmiles(mol);
logInfo("Canonical SMILES: %s", smiles);

%% Section 2: Descriptor Calculation
% ---------------------------------------------------------

mw  = emk.descriptor.molWeight(mol);
logInfo("Molecular Weight: %.3f", mw);

%% Section 3: Fingerprint & Similarity
% ---------------------------------------------------------

mol2 = emk.mol.fromSmiles("CCCO");
fp1  = emk.fingerprint.morgan(mol);
fp2  = emk.fingerprint.morgan(mol2);
sim  = emk.similarity.tanimoto(fp1, fp2);
logInfo("Tanimoto similarity: %.4f", sim);

