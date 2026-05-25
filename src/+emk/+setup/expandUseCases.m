function cfg = expandUseCases(cfg)
% expandUseCases  Expand cfg.useCase.* group flags into cfg.optionalLibraries.*.
%
%   cfg = emk.setup.expandUseCases(cfg)
%
%   Translates high-level use-case toggles set in main_rdkit.m into the
%   individual cfg.optionalLibraries.* flags consumed by install() and
%   installOnline().  Called automatically by both install functions before
%   the optional library install loop.
%
%   Use-case groups and the libraries they activate:
%
%     useCase.qsar    -> pubchempy, mordred
%     useCase.bio     -> biopython
%     useCase.ml      -> torch (CPU-only), torch_geometric
%     useCase.nlp     -> transformers, datasets
%     useCase.docking -> scipy, meeko, prody, vina, pdbfixer
%                        (Online only; Desktop skips vina/pdbfixer gracefully)
%
%   Arguments:
%     cfg (struct) - config struct from loadConfig(). Must contain cfg.useCase.
%                   When cfg.useCase is absent the struct is returned unchanged.
%
%   Returns:
%     cfg (struct) - same struct with cfg.optionalLibraries.* merged in.
%
%   See also: emk.setup.install, emk.setup.installOnline

    if ~isfield(cfg, "useCase")
        return;
    end

    uc = cfg.useCase;

    if isfield(uc, "qsar") && uc.qsar
        cfg.optionalLibraries.pubchempy = true;
        cfg.optionalLibraries.mordred   = true;
    end

    if isfield(uc, "bio") && uc.bio
        cfg.optionalLibraries.biopython = true;
    end

    if isfield(uc, "ml") && uc.ml
        cfg.optionalLibraries.torch           = true;
        cfg.optionalLibraries.torch_geometric = true;
    end

    if isfield(uc, "nlp") && uc.nlp
        cfg.optionalLibraries.transformers = true;
        cfg.optionalLibraries.datasets     = true;
    end

    if isfield(uc, "docking") && uc.docking
        cfg.optionalLibraries.scipy    = true;
        cfg.optionalLibraries.meeko    = true;
        % prody: optional enhancement for mk_prepare_receptor.
        % meeko 0.7.x uses --read_pdb (no prody) by default; prody enables
        % --read_with_prody for advanced use cases (altloc, mmCIF, etc.).
        % Best-effort only: install skipped gracefully if no pre-built wheel.
        cfg.optionalLibraries.prody    = true;
        cfg.optionalLibraries.vina     = true;
        cfg.optionalLibraries.pdbfixer = true;
    end
end
