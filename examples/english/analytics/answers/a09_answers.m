%[text] # A09 Answers: PFAS and Environmental Chemical Screening
%[text] Reference answers for the "Try It Yourself" exercise in a09_pfas_screening.m.
%[text] First, run a09_pfas_screening.m (at least sections 0-5) to
%[text] build the necessary workspace variables. Then check with this file.
addpath(genpath("src"));
emk.setup.initPython();

mol_warmup = emk.mol.fromSmiles("C"); clear mol_warmup;
hasOptTbx   = license("test", "optimization_toolbox");
hasStatsTbx = license("test", "statistics_toolbox");
logInfo("A09 Answers: Setup complete  (OptTbx=%d  StatsTbx=%d)", hasOptTbx, hasStatsTbx);
%%
%[text] ## Reconstructing Prerequisite Variables (Reproduction of a09 Sections 0-5)

CHEMICALS = { ...
    "PFBA",       "OC(=O)C(F)(F)C(F)(F)C(F)(F)F",                                   "PFCA"; ...
    "PFHxA",      "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                    "PFCA"; ...
    "PFOA",       "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",      "PFCA"; ...
    "PFNA",       "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "PFCA"; ...
    "PFDA",       "OC(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "PFCA"; ...
    "PFBS",       "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                        "PFSA"; ...
    "PFHxS",      "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",         "PFSA"; ...
    "PFOS",       "OS(=O)(=O)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "PFSA"; ...
    "6:2FTS",     "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",               "FTS"; ...
    "8:2FTS",     "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F", "FTS"; ...
    "4:2FTS",     "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)F",                              "FTS"; ...
    "Fluoxetine", "CNCCC(c1ccccc1)Oc1ccc(cc1)C(F)(F)F",                            "NonPFAS"; ...
    "Ciprofloxacin","OC(=O)c1cn(C2CCNCC2)c2cc(F)c(nc12)N1CCCC1",                   "NonPFAS"; ...
    "Flurbiprofen","OC(=O)C(C)c1ccc(cc1)-c1cccc(F)c1",                             "NonPFAS"; ...
    "Diflunisal",  "OC(=O)c1ccc(cc1O)-c1ccc(F)cc1F",                               "NonPFAS"; ...
    "Trifluridine","OC1C(CO)OC(n2cc(c(=O)[nH]2)C(F)(F)F)C1F",                      "NonPFAS"; ...
    "Halothane",   "FC(F)(F)C(Cl)Br",                                               "NonPFAS"; ...
    "Sevoflurane", "FC(F)(F)C(F)OCC(F)(F)F",                                        "NonPFAS"; ...
    "Flecainide",  "OCC(NC(=O)c1cc(OCC(F)(F)F)c(cc1OCC(F)(F)F)OCC(F)(F)F)CC",     "NonPFAS"; ...
    "Ethanol",     "CCO",                                                            "NonPFAS"; ...
};

nChem   = size(CHEMICALS, 1);
names   = string(CHEMICALS(:,1));
smiles  = string(CHEMICALS(:,2));
trueCls = string(CHEMICALS(:,3));

mols  = cell(1, nChem);
valid = false(1, nChem);
for k = 1:nChem
    try; mols{k} = emk.mol.fromSmiles(smiles(k)); valid(k) = true; catch; end
end

validMols  = mols(valid);
validNames = names(valid);
validSmiles = smiles(valid);
validCls   = trueCls(valid);
nValid     = sum(valid);

SMARTS_PFCHAIN    = "[#6](F)(F)-[#6](F)(F)";
SMARTS_CF3        = "[#6](F)(F)F";
SMARTS_SULFONYL   = "[#16](=O)(=O)";
SMARTS_POLYFLUORO = "[#6](F)(F)[#6H2]";

hasChain     = emk.mol.hasSubstruct(validMols, SMARTS_PFCHAIN);
hasCF3       = emk.mol.hasSubstruct(validMols, SMARTS_CF3);
hasSulfonyl  = emk.mol.hasSubstruct(validMols, SMARTS_SULFONYL);
hasPolyFluoro = emk.mol.hasSubstruct(validMols, SMARTS_POLYFLUORO);

isPFAS     = hasChain | (hasCF3 & hasChain);
classGuess = repmat("NonPFAS", nValid, 1);
for k = 1:nValid
    if isPFAS(k)
        if hasSulfonyl(k),      classGuess(k) = "PFSA";
        elseif hasPolyFluoro(k), classGuess(k) = "FTS";
        else,                    classGuess(k) = "PFCA";
        end
    end
end

DESCS       = ["LogP", "TPSA", "HeavyAtomCount"];
descTbl     = emk.descriptor.batchCalculate(validMols, DESCS);
logp_vec    = descTbl.LogP;
tpsa_vec    = descTbl.TPSA;
nF_vec      = zeros(nValid, 1);
for k = 1:nValid; nF_vec(k) = count(validSmiles(k), "F"); end

score_logP  = min(max(logp_vec / 10, 0), 1);
score_F     = min(max(nF_vec / 17,   0), 1);
score_TPSA  = exp(-((tpsa_vec - 40).^2) / (2 * 30^2));
scoresMat   = [score_logP, score_F, score_TPSA];

targetScores = zeros(nValid, 1);
for k = 1:nValid
    cls = validCls(k); nm = validNames(k);
    if cls == "NonPFAS",                              targetScores(k) = 0.10;
    elseif any(nm == ["PFOA","PFOS","PFNA","PFDA"]),  targetScores(k) = 0.80;
    elseif any(nm == ["PFHxA","PFHxS","8:2FTS"]),     targetScores(k) = 0.60;
    elseif any(nm == ["PFBA","PFBS","6:2FTS","4:2FTS"]),targetScores(k) = 0.40;
    else,                                             targetScores(k) = 0.50;
    end
end

if hasOptTbx
    objFun = @(w) sum((scoresMat * w - targetScores).^2);
    Aeq = [1,1,1]; beq = 1;
    lb = [0.05;0.05;0.05]; ub = [0.80;0.80;0.80]; w0 = [1/3;1/3;1/3];
    opts = optimoptions("fmincon", Display="off", Algorithm="sqp");
    [w_opt, ~] = fmincon(objFun, w0, [], [], Aeq, beq, lb, ub, [], opts);
else
    w_opt = [1/3;1/3;1/3];
end
concernScore = scoresMat * w_opt;

pfasHitIdx = find(isPFAS);
nHits      = numel(pfasHitIdx);

logInfo("Reconstructing prerequisites: %d chemicals, %d PFAS hits.", nValid, nHits);
%%
%[text] ## Let's Try 1: Carbon count in PFOS vs PFHxS; Fluoxetine non-PFAS; Halothane OECD

logInfo("Let's Try 1 -- Carbon chain length in PFOS and PFHxS:");
for nm = ["PFOS", "PFHxS"]
    idx  = find(validNames == nm);
    fc   = parseSMILES_FCount(validSmiles(idx));
    mol  = validMols{idx};
    d    = emk.descriptor.calculate(mol, ["HeavyAtomCount","MolFormula"]);
    nCF2 = count(validSmiles(idx), "C(F)(F)");
    logInfo("  %s: formula=%s  heavy=%d  CF2 units in SMILES=%d", ...
        nm, d.MolFormula, d.HeavyAtomCount, nCF2);
end
%[text] Answer: PFOS (Perfluorooctanesulfonic acid) has 8 CF2 units in addition to the sulfonyl head group,
%[text]    making a total of 8 fluorinated carbons. PFHxS has 6.
%[text]    Name prefixes: PF = Perfluoro, Hx = Hexa (6), S = Sulfonate,
%[text]    O = Octa (8). These numbers indicate the carbon chain length.
logInfo("Let's Try 1 -- Fluoxetine: Why is it non-PFAS?");
idx_fluox = find(validNames == "Fluoxetine");
mol_fluox = validMols{idx_fluox};
d_fluox   = emk.descriptor.calculate(mol_fluox, "MolFormula");
hasCF3_fluox  = emk.mol.hasSubstruct({mol_fluox}, SMARTS_CF3);
hasChain_fluox = emk.mol.hasSubstruct({mol_fluox}, SMARTS_PFCHAIN);
logInfo("  Fluoxetine formula: %s  hasCF3=%d  hasPFchain=%d", ...
    d_fluox.MolFormula, hasCF3_fluox, hasChain_fluox);
%[text] Answer: Fluoxetine has a single isolated -CF3 group (trifluoromethyl), but
%[text]    the OECD PFAS definition requires at least one -CF2-CF2- unit (two consecutive fully fluorinated carbons).
%[text]    A lone CF3 does not form a perfluoroalkyl chain. hasCF3 is true but hasPFchain is false
%[text]    => Classified as NonPFAS.
%[text]    Important distinction: Many pharmaceuticals contain isolated CF3 or CF2 groups,
%[text]    but are not of concern for environmental persistence like PFAS.
logInfo("Let's Try 1 -- Halothane (FC(F)(F)C(Cl)Br): Is it PFAS by OECD?");
idx_halo  = find(validNames == "Halothane");
mol_halo  = validMols{idx_halo};
hasCF3_halo   = emk.mol.hasSubstruct({mol_halo}, SMARTS_CF3);
hasChain_halo = emk.mol.hasSubstruct({mol_halo}, SMARTS_PFCHAIN);
logInfo("  Halothane: hasCF3=%d  hasPFchain=%d  -> isPFAS=%d", ...
    hasCF3_halo, hasChain_halo, isPFAS(idx_halo));
%[text] Answer: Halothane has a CF3 on one carbon, but the adjacent carbons have Cl and Br, not F.
%[text]    There is no -CF2-CF2- chain. According to the OECD 2021 definition, Halothane is not PFAS.
%[text]    OECD excludes substances where all C-F bonds are directly attached to heteroatoms,
%[text]    but Halothane has a C-C bond (to CHClBr), so this exclusion does not apply.
%[text]    Correct judgment: No perfluoroalkyl chain =>
%[text]    Presence of -CF3 does not make it PFAS.
%%
%[text] ## Let's Try 2: Control hasCF3; Inadequacy of SMARTS_CF3; Improved FTS SMARTS

logInfo("Let's Try 2 -- Non-PFAS control hasCF3:");
nonPFASidx = find(validCls == "NonPFAS");
for k = nonPFASidx'
    if hasCF3(k)
        logInfo("  %-18s: hasCF3=1  isPFAS=%d", validNames(k), isPFAS(k));
    end
end
%[text] Answer: Controls such as Fluoxetine, Trifluridine, Halothane, Sevoflurane, and Flecainide
%[text]    return hasCF3 = true. However, they do not have a perfluoroalkyl chain (SMARTS_PFCHAIN),
%[text]    so they are all classified as NonPFAS. Using hasCF3 alone for PFAS detection
%[text]    causes a very high false positive rate.
logInfo("Let's Try 2 -- Improved FTS SMARTS  [#6](F)(F)-[#6H2]:");
hasFTSjunction = emk.mol.hasSubstruct(validMols, "[#6](F)(F)-[#6H2]");
ftsIdx = find(validCls == "FTS");
logInfo("  True FTS compounds:");
for k = ftsIdx'
    logInfo("    %-10s: hasFTSjunction=%d  classGuess=%s", ...
        validNames(k), hasFTSjunction(k), classGuess(k));
end
logInfo("  Non-FTS compounds that hit the FTS junction SMARTS:");
for k = 1:nValid
    if hasFTSjunction(k) && validCls(k) ~= "FTS"
        logInfo("    %-18s (true=%s)", validNames(k), validCls(k));
    end
end
%[text] Answer: The junction pattern [#6](F)(F)-[#6H2] (CF2 bonded to CH2) is
%[text]    a structural signature of fluorotelomer substances.
%[text]    It correctly identifies 6:2FTS, 8:2FTS, and 4:2FTS, and does not hit PFCA
%[text]    (no CH2 in the fluorinated chain, entirely -CF2-CF2-) or PFSA.
%[text]    This more specific pattern reduces misclassification of FTS as PFCA.
%%
%[text] ## Let's Try 3: Maximum F Score; Ethanol Proxy Score; TPSA Proxy Redesign

[maxFscore, maxFidx] = max(score_F);
logInfo("Let's Try 3 -- Maximum F Score:");
logInfo("  Compound: %s  (nF=%d  score_F=%.3f)", ...
    validNames(maxFidx), nF_vec(maxFidx), maxFscore);
%[text] Answer: score_F is capped at 1.0 for nF >= 17 (normalization constant).
%[text]    Both PFNA (nF=17) and PFDA (nF=19) reach score_F = 1.000. Since max() returns the first
%[text]    occurrence, PFNA is displayed as the winner.
%[text]    In reality, PFDA has the most fluorine atoms (19F: 9 CF2 units + 1 CF3 = 19F).
%[text]    The longest chain compound is the most resistant to degradation -- consistent with Buck et al. (2011).
%[text]    To avoid saturation, it is recommended to normalize by the dataset's nF_max
%[text]    (e.g., nF/19).
ethanolIdx = find(validNames == "Ethanol");
if ~isempty(ethanolIdx)
    logInfo("Let's Try 3 -- Ethanol Proxy Score:");
    logInfo("  nF=%d  score_logP=%.3f  score_F=%.3f  score_TPSA=%.3f", ...
        nF_vec(ethanolIdx), score_logP(ethanolIdx), ...
        score_F(ethanolIdx), score_TPSA(ethanolIdx));
end
%[text] Answer: Ethanol (CCO): No fluorine atoms => score_F = 0; LogP ~ -0.31 =>
%[text]    score_logP = 0 (lower bound fixed); TPSA ~ 20 A^2 => ~0.64 with Gaussian centered at 40
%[text]    -- this is the only non-zero score. Ethanol functions as a negative control with no persistence concerns.
%[text]    All three proxies should be close to 0, but TPSA = 20 is near the Gaussian center (40 A^2),
%[text]    so the TPSA proxy returns a small score.
%[text]    This is a known weakness of this proxy design.
logInfo("Let's Try 3 -- Discussion on TPSA Proxy Redesign:");
pfosIdx = find(validNames == "PFOS");
if ~isempty(pfosIdx)
    logInfo("  PFOS: TPSA=%.1f  score_TPSA=%.3f", tpsa_vec(pfosIdx), score_TPSA(pfosIdx));
end
%[text] Answer: PFOS has TPSA ~ 115 A^2 (sulfonate head group) but high environmental persistence.
%[text]    The Gaussian centered at 40 A^2 assigns a lower TPSA score to PFOS than to short-chain PFCA (TPSA ~ 37 A^2).
%[text]    This is counterintuitive for concern assessment.
%[text]    Redesign proposals:
%[text]    (a) Use a sigmoid: score_TPSA = 1 / (1 + exp(-(TPSA - 60) / 10))
%[text]        -> High concern at both extremes.
%[text]    (b) Use solubility prediction (`emk.descriptor`: MolLogS) instead of TPSA.
%[text]    (c) Replace TPSA with a binary flag for ionic head groups (-COOH / -SO3H).
%%
%[text] ## Let's Try 4: Maximum Weight Proxy; Impact of Changing PFOS Target; Fourth Proxy

logInfo("Let's Try 4 -- Optimizing Proxy Weights:");
[maxW, maxWidx] = max(w_opt);
proxyNames = ["LogP", "F-count", "TPSA"];
logInfo("  w = [%.3f, %.3f, %.3f]  ->  Maximum: %s (%.3f)", ...
    w_opt(1), w_opt(2), w_opt(3), proxyNames(maxWidx), maxW);
%[text] Answer: The reason why w_F (F-count) usually receives the maximum weight:
%[text]    (1) Long-chain PFAS have high nF and high target scores (0.7-0.8).
%[text]    (2) Short-chain PFAS and FTS have moderate nF and moderate target scores (0.4-0.6).
%[text]    (3) Non-PFAS controls have nF = 0 and low target scores (0.1).
%[text]    The F-count proxy completely separates the three concern layers, so fmincon assigns the maximum weight here.
%[text]    This is consistent with Buck et al. (2011) -- C-F bond density is identified as a major factor in environmental persistence.
if hasOptTbx
    logInfo("Let's Try 4 -- Impact of Raising PFOS Target from 0.80 to 0.95:");
    targetScores_mod          = targetScores;
    pfos_idx_local            = find(validNames == "PFOS");
    targetScores_mod(pfos_idx_local) = 0.95;
    objFun_mod = @(w) sum((scoresMat * w - targetScores_mod).^2);
    [w_opt_mod, ~] = fmincon(objFun_mod, w0, [], [], Aeq, beq, lb, ub, [], opts);
    logInfo("  Original Weights: [%.3f, %.3f, %.3f]", w_opt(1), w_opt(2), w_opt(3));
    logInfo("  Modified Weights: [%.3f, %.3f, %.3f]", ...
        w_opt_mod(1), w_opt_mod(2), w_opt_mod(3));
    % Answer: Raising the PFOS target causes the optimizer to assign more weight to the F-count proxy
    %    (PFOS has 17 F atoms). The weight for LogP usually remains unchanged or decreases
    %    (since PFOS already has a high LogP).
    %    This sensitivity analysis shows that the regulatory concern rubric directly determines proxy weights
    %    —— a form of encoding expertise.

    logInfo("Let's Try 4 -- Fourth Proxy: CF2 Chain Length:");
    nCF2_vec = zeros(nValid, 1);
    for k = 1:nValid; nCF2_vec(k) = count(validSmiles(k), "C(F)(F)"); end
    score_CF2 = min(max(nCF2_vec / 10, 0), 1);  % normalize to ~10 CF2 units
    scoresMat4 = [score_logP, score_F, score_TPSA, score_CF2];
    Aeq4 = [1,1,1,1]; beq4 = 1;
    lb4 = [0.05;0.05;0.05;0.05]; ub4 = [0.80;0.80;0.80;0.80];
    w0_4 = [0.25;0.25;0.25;0.25];
    objFun4 = @(w) sum((scoresMat4 * w - targetScores).^2);
    [w_opt4, fval4] = fmincon(objFun4, w0_4, [], [], Aeq4, beq4, lb4, ub4, [], opts);
    logInfo("  4 Proxy Weights: LogP=%.3f  F=%.3f  TPSA=%.3f  CF2=%.3f  (RSS=%.4f)", ...
        w_opt4(1), w_opt4(2), w_opt4(3), w_opt4(4), fval4);
    % Answer: Adding a dedicated CF2 chain length proxy usually reduces the residual (RSS).
    %    Since the chain length proxy and F-count proxy are correlated, the optimizer distributes
    %    weights between them. This indicates that optimal weighting is not unique when multicollinearity is present.
end
%%
%[text] ## Let's Try 5: Highest Tanimoto Pair; PFOA vs PFOS Similarity; Effect of Radius

if nHits > 1
    fps = cell(1, nHits);
    for k = 1:nHits
        fps{k} = emk.fingerprint.morgan(validMols{pfasHitIdx(k)});
    end
    simMat = zeros(nHits, nHits);
    for i = 1:nHits
        for j = 1:nHits
            simMat(i,j) = emk.similarity.tanimoto(fps{i}, fps{j});
        end
    end
    hitNames = validNames(pfasHitIdx);

    logInfo("Let's Try 5 -- Most Similar PFAS Pair:");
    triMask  = triu(true(nHits), 1);
    triVals  = simMat(triMask);
    [maxSim, maxLinIdx] = max(triVals);
    [ii, jj] = ind2sub([nHits, nHits], find(triMask));
    logInfo("  Maximum Tanimoto = %.4f between %s and %s", ...
        maxSim, hitNames(ii(maxLinIdx)), hitNames(jj(maxLinIdx)));
    % Answer: The most similar pair is PFHxA (C6 PFCA) and PFOA (C8 PFCA) with T = 1.0.
    %    This is not an error but reflects how Morgan/ECFP4 functions with repeating chains.
    %    Both molecules have a structure of -CF2- repeats with one -COOH added.
    %    After bit enumeration with radius 2, the unique cyclic substructure set of PFHxA is a strict subset of PFOA
    %    (the extra -CF2- units in PFOA do not add new bit patterns once the chain is long enough).
    %    Tanimoto = |A AND B| / |A OR B| = |PFHxA bits| / |PFOA bits| = 1.0
    %    because all bits of PFHxA appear in PFOA, and PFOA adds no new bits.
    %    This fingerprint saturation is a known limitation of homologous series (Maggiora et al. 2014).
    %    It can be partially mitigated with higher radii or count-based fingerprints (e.g., ECFP6).

    pfoa_idx  = find(hitNames == "PFOA");
    pfos_idx  = find(hitNames == "PFOS");
    if ~isempty(pfoa_idx) && ~isempty(pfos_idx)
        t_pfoa_pfos = simMat(pfoa_idx, pfos_idx);
        logInfo("Let's Try 5 -- PFOA vs PFOS Tanimoto = %.4f", t_pfoa_pfos);
    end
    % Answer: PFOA (C8 PFCA) and PFOS (C8 PFSA) have the same carbon chain length but
    %    different head groups (-COOH vs -SO3H). Tanimoto is usually 0.3–0.5:
    %    the chain substructures overlap, but the head group substructures (carboxylic acid vs sulfonate)
    %    differ, reducing similarity.
    %    Radius 2 captures neighborhoods of 2 bonds. With radius 1, the head group difference has less impact,
    %    and T rises to 0.5–0.7.

    logInfo("Let's Try 5 -- Effect of Reduced Radius (radius=1):");
    fps_r1 = cell(1, nHits);
    for k = 1:nHits
        fps_r1{k} = emk.fingerprint.morgan(validMols{pfasHitIdx(k)}, Radius=1);
    end
    if ~isempty(pfoa_idx) && ~isempty(pfos_idx)
        t_r1 = emk.similarity.tanimoto(fps_r1{pfoa_idx}, fps_r1{pfos_idx});
        logInfo("  PFOA vs PFOS T(radius=1) = %.4f  vs T(radius=2) = %.4f", ...
            t_r1, t_pfoa_pfos);
    end
    % Answer: With radius 1, the algorithm considers only neighborhoods of 1 bond.
    %    The shared -CF2-CF2- chain contributes strongly, and the head group differences
    %    are confined to 1–2 atoms, reducing discriminative power.
    %    Reducing the radius generally increases similarity. This illustrates why Morgan radius 2 (ECFP4) is
    %    a standard choice: it captures a balanced amount of local (bond) and mid-range (substructure)
    %    structural information.
end
%%
%[text] ## Let's Try 6: High Concern Clusters; Anesthetic Exemption; CSV Export

logInfo("Let's Try 6 -- High Concern Clusters (LogP vs F count):");
highConcernMask = concernScore > 0.5 & isPFAS';
logInfo("  Compounds with concern > 0.5 and isPFAS:");
for k = find(highConcernMask)'
    logInfo("    %-12s: LogP=%.2f  nF=%d  score=%.3f  [%s]", ...
        validNames(k), logp_vec(k), nF_vec(k), concernScore(k), validCls(k));
end
%[text] Answer: Long-chain PFCA and PFSA (PFOA, PFOS, PFNA, PFDA, PFHxS) cluster in the
%[text]    upper right region of the LogP vs F count scatter plot (high LogP and high nF).
%[text]    A simple binary screen with a threshold of nF >= 10 captures all long-chain PFAS
%[text]    with zero false negatives in this dataset. However, F count alone misses FTS
%[text]    (few Fs but a PFCA precursor).
logInfo("Let's Try 6 -- Anesthetic Exemption Design:");
volatileMask  = zeros(nValid, 1);
for k = 1:nValid
    d_mw = emk.descriptor.calculate(validMols{k}, "MolWt");
    if d_mw.MolWt < 200
        volatileMask(k) = 1;
    end
end
exemptIdx = find(volatileMask & ~isPFAS');
logInfo("  Volatile (MW<200) non-PFAS compounds in the dataset:");
for k = exemptIdx'
    logInfo("    %-14s: nF=%d  isPFAS=%d", validNames(k), nF_vec(k), isPFAS(k));
end
%[text] Answer: Halothane (MW=197) and Sevoflurane (MW=200) are volatile anesthetic gases.
%[text]    The "volatile exemption" flag (MW < 200) excludes them from high concern layers
%[text]    even with moderate F count scores. This mimics real regulatory practice:
%[text]    Anesthetic gases are regulated by occupational exposure limits, not environmental persistence frameworks.
%[text]    Exemption rubric: volatilityScore = (MW < 200) ? 0 : 1
%[text]    Application: finalScore = concernScore .* volatilityScore
%[text]
%[text] CSV Export
[~, rankOrder] = sort(concernScore, "descend");
rankVec = zeros(nValid, 1);
rankVec(rankOrder) = (1:nValid)';
summaryTbl = table( ...
    rankVec, validNames, validCls, classGuess, ...
    round(logp_vec,2), nF_vec, round(tpsa_vec,1), round(concernScore,3), ...
    VariableNames=["Rank","Name","TrueClass","FlaggedClass","LogP","nF","TPSA_A2","ConcernScore"]);
summaryTbl = sortrows(summaryTbl, "Rank");

runDir = makeRunDir("Prefix", "a09_answers");
csvPath = fullfile(runDir, "pfas_screening_report.csv");
writetable(summaryTbl, csvPath);
logInfo("Let's Try 6 -- Exported summary table to: %s", csvPath);
disp(summaryTbl);
%[text] Answer: The exported CSV can be loaded into Excel or any analysis tool.
%[text]    Recommended columns to retain for non-expert stakeholders:
%[text]    Rank, Name, FlaggedClass, ConcernScore.
%[text]    Columns needing explanation in footnotes: LogP (lipophilicity), nF (fluorine count).
%[text]    TrueClass (for internal validation) and TPSA_A2 can be removed if unnecessary.
%[text]    Example introduction for the report: "Five chemicals with scores > 0.5 have been
%[text]    flagged for priority review as long-chain PFAS (Table 1)."
logInfo("A09 Answer: Completed.");
%%
%[text] ## Local Helper Function

function nF = parseSMILES_FCount(smi)
%[text] Returns the number of occurrences of "F" in a SMILES string.
    nF = count(smi, "F");
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]%   data: {"layout":"inline","rightPanelPercent":40}
%---
