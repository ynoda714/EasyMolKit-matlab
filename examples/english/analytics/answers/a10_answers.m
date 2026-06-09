%[text] # A10 Answers: Lead Optimization -- Multi-objective Property Optimization
%[text] Reference answers for the "Try It Yourself" exercise in a10_lead_optimization.m.
%[text] First, execute a10_lead_optimization.m (at least sections 0-3) to
%[text] build the necessary workspace variables. Then check with this file.
%[text] addpath(genpath("src"));
%[text] emk.setup.initPython();

%[text] mol_warmup = emk.mol.fromSmiles("C"); clear mol_warmup;
%[text] hasOptTbx  = license("test", "optimization_toolbox");
%[text] logInfo("A10 Answers: Setup complete (OptTbx=%d)", hasOptTbx);
%%
%[text] ## Reconstructing Prerequisites (Reproduction of a10 Sections 0-3)

LEAD_SMILES = "CC(=O)Oc1ccccc1C(=O)O";
LEAD_NAME   = "Aspirin";
TARGET_LOGP = 2.0;
TARGET_TPSA = 70.0;
TARGET_MW   = 350.0;

mol_lead = emk.mol.fromSmiles(LEAD_SMILES);
d_lead   = emk.descriptor.calculate(mol_lead);

DATA_FILE = "data/list/fda_drugs.csv";
rawTbl    = readtable(DATA_FILE, TextType="string");
nLib      = height(rawTbl);
mols      = cell(1, nLib);
valid     = false(1, nLib);
for k = 1:nLib
    try; mols{k} = emk.mol.fromSmiles(rawTbl.SMILES(k)); valid(k) = true; catch; end
end
validIdx  = find(valid);
mols      = mols(validIdx);
libNames  = rawTbl.Name(validIdx);
nValid    = numel(mols);

DESCS   = ["LogP", "TPSA", "MolWt"];
descTbl = emk.descriptor.batchCalculate(mols, DESCS);
logp_vec = descTbl.LogP;
tpsa_vec = descTbl.TPSA;
mw_vec   = descTbl.MolWt;
propMat  = [logp_vec, tpsa_vec, mw_vec];

d_logp = @(p) ...
    ((p >= 0.5) & (p <  2.0)) .* ((p - 0.5) ./ 1.5) + ...
    ((p >= 2.0) & (p <= 4.0)) .* ((4.0 - p) ./ 2.0) + ...
    (p == 2.0);
d_tpsa = @(p) (p <= 60) .* 1.0 + (p > 60 & p <= 90) .* ((90 - p) ./ 30);
d_mw   = @(p) ...
    ((p >= 150) & (p <  350)) .* ((p - 150) ./ 200) + ...
    ((p >= 350) & (p <= 500)) .* ((500 - p) ./ 150) + ...
    (p == 350);

d1 = d_logp(logp_vec);
d2 = d_tpsa(tpsa_vec);
d3 = d_mw(mw_vec);
D  = (d1 .* d2 .* d3) .^ (1/3);

logInfo("A10 Prerequisite reconstruction complete (%d valid molecules).", nValid);
%%
%[text] ## Let's Try 1: Aspirin MW "Budget", Ibuprofen Comparison, CNS vs Peripheral Targets
%[text] Q: Aspirin MW = 180 g/mol is significantly below the 350 g/mol target.
%[text] Answer: A low MW gives medicinal chemists a "molecular weight budget" (~170 Da),
%[text]    allowing the addition of substituents to improve potency, selectivity, and PK
%[text]    without violating Lipinski's Ro5 (MW <= 500). Starting with a small lead compound
%[text]    is strongly recommended in modern drug discovery.
mol_ibu = emk.mol.fromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
d_ibu   = emk.descriptor.calculate(mol_ibu);

logInfo("Let's Try 1 -- Aspirin vs Ibuprofen vs Target:");
logInfo("  %12s: LogP=%+5.2f TPSA=%5.1f MW=%5.1f", ...
    LEAD_NAME, d_lead.LogP, d_lead.TPSA, d_lead.MolWt);
logInfo("  %12s: LogP=%+5.2f TPSA=%5.1f MW=%5.1f", ...
    "Ibuprofen", d_ibu.LogP, d_ibu.TPSA, d_ibu.MolWt);
logInfo("  %12s: LogP=%+5.2f TPSA=%5.1f MW=%5.1f", ...
    "Target", TARGET_LOGP, TARGET_TPSA, TARGET_MW);

%[text] Answer: Ibuprofen (MW=206, LogP~3.1, TPSA~37) is closer to the LogP target than
%[text]    Aspirin (MW=180, LogP=1.3, TPSA=63). TPSA is significantly below the 70 A^2 target.
%[text]    However, neither reaches MW=350 ── both are "lean" NSAIDs.
%[text]
%[text] Peripheral anti-inflammatory drug TARGET_TPSA: < 90 A^2 (Veber rule).
%[text] CNS drugs: < 60 A^2 (blood-brain barrier criteria).
%[text] Both Aspirin and Ibuprofen meet the peripheral criteria. For CNS criteria,
%[text] Aspirin is borderline (TPSA=63, close to the 60 A^2 limit).
%%
%[text] ## ✏️ Try It 2: LogP-TPSA Correlation; Number of Molecules Meeting All 3 Criteria; Ro5 Box

r_logp_tpsa = corr(logp_vec, tpsa_vec);
logInfo("Try It 2 -- Pearson r(LogP, TPSA) = %+.3f", r_logp_tpsa);

%[text] **Answer:** r ~ -0.3 to -0.5 (weak negative correlation). Lipophilic compounds (high LogP) tend
%[text]    to have fewer polar groups (low TPSA). However, the correlation is weak with many exceptions
%[text]    (e.g., chlorinated aromatics have high LogP, but chlorine substitution keeps TPSA moderate).
n_all_three = sum((logp_vec >= 1.0 & logp_vec <= 3.5) & ...
                   tpsa_vec <= 90 & mw_vec <= 400);
logInfo("Try It 2 -- Molecules meeting all 3 strict criteria: %d / %d", ...
    n_all_three, nValid);

%[text] **Answer:** In this FDA library, 75 / 200 (37.5%%) molecules satisfy all 3 criteria simultaneously.
%[text]    The "oral drug sweet spot" is strongly reflected in approved drugs —
%[text]    molecules that fail ADMET rarely make it to FDA approval.
%[text]
%[text] Lipinski vs Veber comparison:
n_lipinski = sum(logp_vec <= 5 & mw_vec <= 500);
n_veber    = sum(tpsa_vec <= 90 & mw_vec <= 400);
logInfo("Try It 2 -- Lipinski (MW<=500, LogP<=5): %d molecules", n_lipinski);
logInfo("Try It 2 -- Veber (TPSA<=90, MW<=400):   %d molecules", n_veber);

%[text] **Answer:** Most FDA drugs pass Lipinski (~90%%). The stricter
%[text]    Veber / 400 Da guideline passes fewer (~60-70%%).
%[text]    The FDA library is biased toward approved drugs, which by definition have acceptable human absorption.
%%
%[text] ## Let's Try 3: Arithmetic Mean vs Geometric Mean; CNS TPSA Target; Max d_LogP & Min d_TPSA
%[text] Desirability of default (geometric) vs arithmetic mean
D_arith = (d1 + d2 + d3) / 3;
[D_geo_top10,  idx_geo]   = maxk(D,       10);
[D_arith_top10,idx_arith] = maxk(D_arith, 10);

overlap_top10 = numel(intersect(idx_geo, idx_arith));
logInfo("Let's Try 3 -- Top 10 Overlap (Geometric vs Arithmetic Mean): %d / 10", overlap_top10);

%[text] Answer: Arithmetic mean allows d_i = 0 if the other two are high.
%[text]    Geometric mean collapses to 0 if any d_i = 0.
%[text]    Thus, arithmetic mean promotes "molecules with one failing property but two excellent ones"
%[text]    — behavior usually undesirable in drug discovery.
%[text]    In this FDA library, top 10 overlap is 9/10 —
%[text]    top candidates have all three properties within range (no d_i = 0),
%[text]    so ranking remains mostly unchanged.
%[text]    In libraries with many "outlier properties" (e.g., very high MW or extreme LogP),
%[text]    overlap significantly decreases (5-7).
%[text]
%[text] CNS TPSA Target: Narrow down to [0, 60]
d_tpsa_cns = @(p) ...
    (p <= 45) .* 1.0 + ...
    (p > 45 & p <= 60) .* ((60 - p) ./ 15) + ...
    (p > 60) .* 0;

D_cns = (d_logp(logp_vec) .* d_tpsa_cns(tpsa_vec) .* d_mw(mw_vec)) .^ (1/3);
n_cns_survive = sum(D_cns > 0);
logInfo("Let's Try 3 -- Molecules with D_CNS > 0 (TPSA < 60): %d / %d", ...
    n_cns_survive, nValid);

%[text] Answer: In this library, 50 / 200 (25%) pass the CNS criteria (TPSA < 60).
%[text]    Many FDA drugs target peripheral tissues,
%[text]    not optimized for CNS permeability.
%[text]
%[text] Molecule with highest d_LogP but lowest d_TPSA
[~, i_top_logp] = max(d1);
[~, i_low_tpsa] = min(d2);
logInfo("Let's Try 3 -- Highest d_LogP: %s (LogP=%.2f, TPSA=%.1f)", ...
    libNames(i_top_logp), logp_vec(i_top_logp), tpsa_vec(i_top_logp));
logInfo("Let's Try 3 -- Lowest d_TPSA:  %s (LogP=%.2f, TPSA=%.1f)", ...
    libNames(i_low_tpsa), logp_vec(i_low_tpsa), tpsa_vec(i_low_tpsa));

%[text] Answer: The molecule with lowest d_TPSA is PRAZOSIN (LogP=1.78, TPSA=107, d_TPSA=0).
%[text]    LogP is moderate, but it has a piperazine ring, furan ring, and two amide groups,
%[text]    each contributing to TPSA. This shows TPSA is determined by the number of polar groups,
%[text]    not lipophilicity. The molecule with highest d_LogP is METOCLOPRAMIDE (LogP=2.00, d_LogP=1.0).
%[text]    LogP perfectly matches the target peak value. They are different molecules,
%[text]    confirming d_LogP and d_TPSA as independent axes of drug-likeness.
%%
%[text] ## Let's Try 4: Interpretation of fgoalattain / gamma; Change MW Target; Change Weights

if ~hasOptTbx
    logWarn("Let's Try 4 -- Optimization Toolbox unavailable. Displaying concept only.");
    logInfo("Let's Try 4 -- gamma <= 0 means all goals are within the convex hull.");
    logInfo("           gamma > 0 means the target is outside.");
    logInfo("           Verification method: Execute section 4 of a10 with Optimization Toolbox.");
else
    % --- Impact of Tightening MW Target ---
    propGoal_tight = [TARGET_LOGP; TARGET_TPSA; 250.0];  % MW = 250 instead of 350
    propWeight_def = [1.0; 1.0; 1.0];

    fun_blend = @(x) propMat' * x(:);
    n = nValid;
    x0  = ones(n, 1) / n;
    opt = optimoptions("fgoalattain", Display="off");

    [~, ~, gamma_350, ~] = fgoalattain(fun_blend, x0, ...
        [TARGET_LOGP; TARGET_TPSA; 350], propWeight_def, ...
        [], [], ones(1,n), 1, zeros(n,1), ones(n,1), [], opt);
    [~, ~, gamma_250, ~] = fgoalattain(fun_blend, x0, ...
        propGoal_tight, propWeight_def, ...
        [], [], ones(1,n), 1, zeros(n,1), ones(n,1), [], opt);

    logInfo("Let's Try 4 -- gamma (MW Target=350): %.4f", gamma_350);
    logInfo("Let's Try 4 -- gamma (MW Target=250): %.4f", gamma_250);

    % Answer: Even if TARGET_MW is narrowed from 350 to 250, gamma does not change.
    %    This is because the library's convex hull spans the range MW = [46, 924] g/mol.
    %    Since any MW within that range can be trivially achieved by convex combination,
    %    both targets are deep inside the convex hull (gamma << 0).
    %    To observe gamma approaching zero, set TARGET_MW below the library's minimum
    %    (e.g., 40 g/mol) or above the maximum (e.g., 1000 g/mol).

    % --- Impact of Doubling TPSA Weight ---
    propWeight_tpsa = [0.5; 2.0; 1.0];
    [~, fval_wtpsa, gamma_wtpsa, ~] = fgoalattain(fun_blend, x0, ...
        [TARGET_LOGP; TARGET_TPSA; TARGET_MW], propWeight_tpsa, ...
        [], [], ones(1,n), 1, zeros(n,1), ones(n,1), [], opt);

    propSigma = std(propMat);
    propMatN  = propMat ./ propSigma;
    fvalN     = fval_wtpsa(:)' ./ propSigma;
    dist_tw   = sqrt(sum((propMatN - fvalN).^2, 2));
    [~, top_tw] = min(dist_tw);
    logInfo("Let's Try 4 -- Top 1 (TPSA Weight Doubled): %s", libNames(top_tw));

    % Answer: With double TPSA weight, fgoalattain penalizes deviations from
    %    TARGET_TPSA more strongly. The optimal blend shifts toward
    %    molecules with TPSA close to 70. The top-1 real molecule may
    %    change if there is a drug with good TPSA but only average LogP/MW.
end
%%
%[text] ## Let's Try 5: Interpretation of the Pareto Front (Conceptual Questions -- No Code Required)

logInfo("Let's Try 5 -- Pareto Front Conceptual Answers:");
logInfo("  Q1: TPSA decreases with an increase in LogP (negative trade-off).");
logInfo("      In a well-organized distribution, it's monotonic; exceptions exist with diverse scaffolds.");
logInfo("      Flat regions = many drugs with similar TPSA in the same LogP range.");
logInfo("      ");
logInfo("  Q2: Ideal [LogP=2, TPSA=70] is likely on the Pareto front (difficult to achieve).");
logInfo("      Libraries usually do not contain drugs with these combined properties.");
logInfo("      ");
logInfo("  Q3: LogP-MW correlation is +0.337 in this library (weak positive correlation).");
logInfo("      Similar magnitude to LogP-TPSA (-0.333) but opposite sign.");
logInfo("      Larger, heavier scaffolds tend to be more lipophilic (alkyl chains/ring additions).");
logInfo("      Polar groups lower LogP but have less impact on MW than TPSA increase.");
logInfo("      ");

%[text] Verification of Pearson Correlation
r_logp_mw = corr(logp_vec, mw_vec);
logInfo("  Pearson r(LogP,TPSA)=%.3f  r(LogP,MW)=%.3f", ...
    corr(logp_vec,tpsa_vec), r_logp_mw);

%[text] Q4: Knee Point -- Maximum |Δ TPSA| / |Δ LogP|
%[text] Requires pareto_logp / pareto_tpsa vectors from Section 5.
%[text] Execute a10 Section 5 first, then run the following:
%[text]   dTPSA = diff(pareto_tpsa);
%[text]   dLogP = diff(pareto_logp);
%[text]   [~, knee] = max(abs(dTPSA ./ dLogP));
%[text]   logInfo("Knee at LogP = %.2f, TPSA = %.1f", pareto_logp(knee), pareto_tpsa(knee));
%%
%[text] ## Let's Try 6: d_TPSA vs d_LogP Champion; PAINS Filter; Synthetic Complexity
%[text] Molecules with highest d_TPSA vs highest d_LogP
[~, i_best_tpsa] = max(d2);
[~, i_best_logp] = max(d1);
same_molecule = isequal(i_best_tpsa, i_best_logp);

logInfo("Let's Try 6 -- Highest d_TPSA: %s (TPSA=%.1f, LogP=%.2f, D=%.3f)", ...
    libNames(i_best_tpsa), tpsa_vec(i_best_tpsa), logp_vec(i_best_tpsa), D(i_best_tpsa));
logInfo("Let's Try 6 -- Highest d_LogP: %s (LogP=%.2f, TPSA=%.1f, D=%.3f)", ...
    libNames(i_best_logp), logp_vec(i_best_logp), tpsa_vec(i_best_logp), D(i_best_logp));
logInfo("Let's Try 6 -- Same molecule? %d  (Indicates multi-objective tension)", same_molecule);

%[text] Answer: d_TPSA = 1.0 applies to all molecules with TPSA <= 60, so max(d2) returns
%[text]    the first matching molecule in array order (NICOTINE, TPSA=16.1, LogP=1.85).
%[text]    "Highest d_TPSA" is not unique ── 50 molecules tie with d_TPSA=1.
%[text]    LogP champion is METOCLOPRAMIDE (LogP=2.00, d_LogP=1.0).
%[text]    LogP perfectly matches the target value 2.0, achieving the highest d_LogP score.
%[text]    They are different molecules (same_molecule=0),
%[text]    illustrating the essence of multi-objective tension in lead optimization.
%[text]
%[text] PAINS Filter Integration (execute if emk.filter.loadPainsSmarts is available)
try
    pains_smarts = emk.filter.loadPainsSmarts();
    [D_sorted_f, rank_by_D] = sort(D, "descend");
    hasPains_top10 = false(1, 10);
    for k = 1:10
        idx_k = rank_by_D(k);
        hasPains_top10(k) = emk.mol.hasSubstruct(mols{idx_k}, pains_smarts);
    end
    n_pains_free = sum(~hasPains_top10);
    logInfo("Let's Try 6 -- PAINS-free in top 10 candidates: %d / 10", n_pains_free);
catch
    logWarn("Let's Try 6 -- PAINS filter unavailable (emk.filter.loadPainsSmarts).");
    logInfo("            Install PAINS SMARTS data or use S03 workflow.");
end

%[text] Q: Your Lead Compound ── Template
%[text] Caffeine as a lead for CNS target example (LogP 1-3, TPSA < 60, MW 200-400):
%[text]   LEAD_SMILES = "Cn1cnc2c1c(=O)n(C)c(=O)n2C";
%[text]   TARGET_LOGP = 1.5; TARGET_TPSA = 45.0; TARGET_MW = 280.0;
%[text]   Re-execute from Section 1.
logInfo("A10 Answer: Completed.");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]%   data: {"layout":"inline","rightPanelPercent":40}
%---
