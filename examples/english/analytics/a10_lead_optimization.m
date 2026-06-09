%[text] # A10: Lead Optimization — Multi-objective Property Optimization
%[text] EasyMolKit Analytics — Layer 3
%[text]
%[text] **Story**
%[text] A drug discovery chemist has selected Aspirin as a structural lead for developing a new anti-inflammatory analgesic. Aspirin is effective, but its molecular weight (180 g/mol) is low for a modern drug candidate, leaving room to add substituents to enhance specificity.
%[text] Challenge: Adding functional groups to improve target binding can simultaneously alter multiple physicochemical properties.
%[text] - Adding lipophilic groups increases LogP (the logarithm of the octanol/water partition coefficient), improving membrane permeability, but also increases molecular weight (MW) and decreases water solubility.
%[text] - Adding polar groups to improve solubility increases the topological polar surface area (TPSA), potentially reducing oral absorption.
%[text] These trade-offs define the "multi-objective" nature of lead optimization.
%[text] Since it is difficult for a single molecule to simultaneously optimize all desirable properties, chemists need to find the "least bad" compromise.
%[text] In this exercise, you will:
%[text]
%[text] 1. Evaluate Aspirin against three oral drug design targets:
%[text] LogP = 2.0, TPSA = 70 Å², MW = 350 g/mol.
%[text] 2. Load a library of 200 FDA-approved drugs as an analog reference and explore their 3D property space.
%[text] 3. Convert raw properties into dimensionless desirability scores using the Derringer-Suich method.
%[text] 4. Formulate and solve a goal attainment problem with fgoalattain (Optimization Toolbox) to find the best blend of property space that meets all targets simultaneously.
%[text] 5. Track the Pareto front of LogP and TPSA trade-offs using the epsilon constraint method (linprog).
%[text] 6. Select and visualize top lead optimization candidates.
%[text]
%[text] **Learning Objectives**
%[text] - Understand why lead optimization is inherently multi-objective
%[text] - Encode pharmacological intuition as Derringer-Suich desirability
%[text] - Minimize maximum slack across all property goals with fgoalattain (Optimization Toolbox)
%[text] - Interpret the achievement factor gamma as a progress indicator
%[text] - Apply the epsilon constraint method and linprog to track the Pareto front
%[text] - Distinguish between Pareto-optimal and suboptimal candidates
%[text]
%[text] **Prerequisites**
%[text] - Completion of S02 (Drug Filter) — Context of Lipinski's Ro5
%[text] - Recommended: A03 (QSAR Regression) — Context of QSAR modeling
%[text] - Optimization Toolbox (fgoalattain, linprog) — If unavailable, Sections 4-5 use manual fallback, so all concepts are learnable.
%[text]
%[text] **Operating Environment**
%[text] Compatibility Summary:
%[text] Student License / Campus-wide — Full Optimization Toolbox
%[text] MATLAB Online Free Version — Optimization Toolbox not supported
%[text] MATLAB Online (Individual / Campus) — Optimization Toolbox available
%[text] Without the toolbox, Sections 4-5 switch to manual alternatives via exhaustive search.
%[text]
%[text] Estimated Time Required: 45-90 minutes
%[text]
%[text] **Data:**
%[text] data/list/fda_drugs.csv — 200 FDA-approved drugs (ChEMBL, CC-BY-SA 3.0)
%[text]
%[text] **References**
%[text] Derringer G (1980) Simultaneous optimization of several response variables. J Quality Technology 12:214-219. doi:10.1080/00224065.1980.11980968
%[text] Veber DF et al. (2002) Molecular properties that influence the oral bioavailability of drug candidates. J Med Chem 45:2615-2623. doi:10.1021/jm020017n
%[text] Lipinski CA et al. (2001) Experimental and computational approaches to estimate solubility and permeability in drug discovery and development. Adv Drug Deliv Rev 46:3-26. doi:10.1016/S0169-409X(00)00129-0
%[text] Cohon JL & Marks DH (1975) A review and evaluation of multiobjective programming techniques. Water Resour Res 11:208-220. doi:10.1029/WR011i002p00208 (epsilon constraint method)
%[text] Charnes A & Cooper WW (1977) Goal programming and multiple objective optimization. Eur J Oper Res 1:39-54. doi:10.1016/S0377-2217(77)81007-2 (foundation of goal attainment)
%[text]
%[text] Execution: Run sections one by one with Ctrl+Enter
%%
%[text] ## Section 0: Setup

% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
emk.setup.initPython();

%[text] Warm up the Python/RDKit process before main execution
mol_warmup = emk.mol.fromSmiles("C");   % Methane -- lightweight
clear mol_warmup;

%[text] Check if the Optimization Toolbox is available
hasOptTbx = license("test", "optimization_toolbox");
if hasOptTbx
    logInfo("A10: Optimization Toolbox detected -- full optimization enabled.");
else
    logWarn("A10: Optimization Toolbox not detected.");
    logWarn("     Sections 4-5 will use manual grid search fallback.");
    logWarn("     Install Optimization Toolbox for fgoalattain / linprog.");
end

logSection("A10", "Section 0: Setup", "Analytics L3");
logInfo("A10: Setup complete.");
%%
%[text] ## Section 1: Lead Compounds and the Multi-Objective Dilemma
%[text]
%[text] Setup is complete. First, let's calculate the properties of the lead compound, Aspirin, and see how far it deviates from the targets.
%[text] Understand the structure of the multi-objective trade-off on three axes: LogP (logarithm of lipophilicity), TPSA (topological polar surface area), and MW (molecular weight).
%[text]
%[text] ### Concept: The Three Pillars of Oral Drug-Likeness
%[text]
%[text] In modern oral drug discovery, molecules are evaluated based on three major physicochemical thresholds (based on Veber 2002 and Lipinski 2001).
%[text]
%[text] (1)  LogP (logarithm of lipophilicity)
%[text] Range: 1.0 <= LogP <= 3.5 (optimal balance of membrane permeability and solubility)
%[text] Too low: Low membrane permeability
%[text] Too high: Low solubility, risk of CYP metabolism
%[text]
%[text] (2)  TPSA (topological polar surface area, Å²)
%[text] Target: <= 90 Å² (Veber rule for good oral absorption)
%[text] Target: <= 60 Å² (CNS penetration -- blood-brain barrier)
%[text] Too high: TPSA > 130 Å² -> generally low passive absorption
%[text]
%[text] (3)  MW (molecular weight, g/mol)
%[text] Target: <= 400 g/mol (stricter "beyond-Ro5" guideline)
%[text] Lipinski: <= 500 g/mol
%[text] Too high: Reduced absorption, may be recognized by efflux transporters
%[text]
%[text] Target achievement goals (used in Sections 3-5):
%[text] LogP* = 2.0, TPSA* = 70.0 Å², MW* = 350.0 g/mol
%[text]
%[text] Lead Compound
logSection("A10", "Section 1: Lead Compounds and the Multi-Objective Dilemma", "Analytics L3");
LEAD_SMILES = "CC(=O)Oc1ccccc1C(=O)O";   % Aspirin
LEAD_NAME   = "Aspirin";

%[text] Property Targets
TARGET_LOGP = 2.0;    % Optimal lipophilicity
TARGET_TPSA = 70.0;   % A^2, balanced polar surface area
TARGET_MW   = 350.0;  % g/mol, strict drug-likeness

mol_lead  = emk.mol.fromSmiles(LEAD_SMILES);
d_lead    = emk.descriptor.calculate(mol_lead);

logInfo("Lead Compound: %s", LEAD_NAME);
logInfo("  LogP = %+.2f   Target = %.1f  Deviation = %+.2f", ...
    d_lead.LogP, TARGET_LOGP, d_lead.LogP - TARGET_LOGP);
logInfo("  TPSA = %5.1f A^2   Target = %.1f  Deviation = %+.1f", ...
    d_lead.TPSA, TARGET_TPSA, d_lead.TPSA - TARGET_TPSA);
logInfo("  MW   = %5.1f g/mol  Target = %.1f  Deviation = %+.1f", ...
    d_lead.MolWt, TARGET_MW, d_lead.MolWt - TARGET_MW);

%[text] **💡 Observation Point 1**
%[text] Aspirin's MW (180 g/mol) is significantly below the target of 350 g/mol.
%[text] What does this imply for medicinal chemists developing analogs?
%[text] (Hint: A lead with low molecular weight has more "room" for adding substituents.)
%[text] Let's calculate the properties of Ibuprofen and compare it with Aspirin:
%[text] mol_ibu = emk.mol.fromSmiles("CC(C)Cc1ccc(cc1)C(C)C(=O)O");
%[text] d_ibu   = emk.descriptor.calculate(mol_ibu);
%[text] Determine which molecule is closer to all three targets simultaneously.
%[text] The targets here are for a hypothetical CNS drug.
%[text] How should TARGET_TPSA be adjusted for a peripheral anti-inflammatory drug?
%[text] (Hint: CNS drugs require TPSA < 60, but peripheral drugs can allow < 90.)
% ... (Try writing code here)
%%
%[text] ## Section 2: Building Reference Library and Calculating Properties
%[text]
%[text] We have checked the properties of the lead compound. Next, we will load 200 FDA-approved drugs as a reference library and explore the chemical space.
%[text] Let's see how actual drugs are distributed in the property space.
%[text]
%[text] ### Concept: Using FDA Drugs as "Reference Chemical Space"
%[text]
%[text] In lead optimization, the reference library serves two purposes:
%[text] (a) To define a realistically accessible chemical space --
%[text] If similar drugs already exist and have passed clinical trials, nearby analogs are more likely to succeed.
%[text] (b) To provide a "target distribution" for property optimization:
%[text] The library shows the LogP (logarithm of partition coefficient) / TPSA (topological polar surface area) / MW (molecular weight) values achieved by actual drugs.
%[text]
%[text] Here, we use 200 FDA-approved drugs from ChEMBL as a reference.
%[text] In actual projects, commercially available building block libraries or computationally generated sets of analogs are used.
logSection("A10", "Section 2: Building Reference Library and Calculating Properties", "Analytics L3");
DATA_FILE = "data/list/fda_drugs.csv";
logInfo("Loaded %d entries from %s", height(rawTbl), DATA_FILE);

%[text] Parsing molecules.
nLib  = height(rawTbl);
mols  = cell(1, nLib);
valid = false(1, nLib);

for k = 1:nLib
    try
        mols{k} = emk.mol.fromSmiles(rawTbl.SMILES(k));
        valid(k) = true;
    catch
        % Skip unparsable entries
    end
end

validIdx = find(valid);
mols     = mols(validIdx);
libNames = rawTbl.Name(validIdx);
libSmiles = rawTbl.SMILES(validIdx);
nValid   = numel(mols);
logInfo("Parsed %d / %d molecules.", nValid, nLib);

%[text] Calculating descriptors: LogP, TPSA, MW
DESCS   = ["LogP", "TPSA", "MolWt"];
descTbl = emk.descriptor.batchCalculate(mols, DESCS);

logp_vec = descTbl.LogP;    % N x 1 double
tpsa_vec = descTbl.TPSA;    % N x 1 double
mw_vec   = descTbl.MolWt;   % N x 1 double

%[text] Property matrix: N x 3  [LogP, TPSA, MW]
propMat = [logp_vec, tpsa_vec, mw_vec];

logInfo("Library property range:");
logInfo("  LogP : [%.2f, %.2f]  Mean = %.2f", ...
    min(logp_vec), max(logp_vec), mean(logp_vec));
logInfo("  TPSA : [%.1f, %.1f] A^2  Mean = %.1f", ...
    min(tpsa_vec), max(tpsa_vec), mean(tpsa_vec));
logInfo("  MW   : [%.1f, %.1f] g/mol  Mean = %.1f", ...
    min(mw_vec), max(mw_vec), mean(mw_vec));

%[text] --- 3D Property Space Plot ---
figure("Name", "A10 Sec2: FDA Drug Property Space (LogP vs TPSA vs MW)");
scatter3(logp_vec, tpsa_vec, mw_vec, 30, mw_vec, "filled", ...
    MarkerFaceAlpha=0.6);
hold on;
%[text] Mark the lead compound.
scatter3(d_lead.LogP, d_lead.TPSA, d_lead.MolWt, 120, "r", ...
    "^", "filled", DisplayName=LEAD_NAME);
%[text] Mark the target point.
scatter3(TARGET_LOGP, TARGET_TPSA, TARGET_MW, 150, "k", ...
    "pentagram", LineWidth=2, DisplayName="Target");
xlabel("LogP");
ylabel("TPSA (A^2)");
zlabel("MW (g/mol)");
title("FDA Drug Property Space");
cb = colorbar();
cb.Label.String = "MW (g/mol)";
legend(Location="best");
grid on;

%[text] **💡 Observation Point 2**
%[text] Rotate the 3D scatter plot to see if there is a visible correlation between LogP and TPSA.
%[text] Is there a visible correlation between LogP and TPSA?
%[text] (Hint: Lipophilic compounds tend to have fewer polar groups.)
%[text] Let's calculate the Pearson correlation: corr(logp_vec, tpsa_vec)
%[text] Check how many library molecules already meet all three targets:
%[text] Count the number of molecules satisfying 1.0 <= LogP <= 3.5 AND TPSA <= 90 AND MW <= 400.
%[text] Count: sum((logp_vec >= 1.0 & logp_vec <= 3.5) & ...
%[text] tpsa_vec <= 90 & mw_vec <= 400)
%[text] Try adding Lipinski's Ro5 boundaries as a transparent box in the 3D plot:
%[text] Consider the conditions [LogP <= 5, TPSA (not an Ro5 criterion), MW <= 500].
%[text] See how it compares with the Veber/strict "350 g/mol" guideline.
% ... (Try writing the code here)
%%
%[text] ## Section 3: Desirability Function -- From Properties to Scores
%[text]
%[text] The property space of the reference library has been revealed. Next, we will convert each property to a score in [0,1] using the Derringer-Suich desirability function.
%[text] Learn how to optimize multiple properties simultaneously and express them with a single metric.
%[text]
%[text] ### Concept: Derringer-Suich Desirability (1980)
%[text]
%[text] The desirability function maps raw property values to a dimensionless score d_i in [0, 1].
%[text] d_i = 1.0  --  when the property perfectly matches the ideal target
%[text] d_i = 0.0  --  when the property is outside the acceptable range (below lower limit or above upper limit)
%[text]
%[text] The standard form for a "two-sided" (target value) goal is as follows:
%[text]
%[text] d_i(p) = 0                                when p < L_i
%[text] = ((p - L_i) / (T_i - L_i))^s     when L_i <= p <= T_i
%[text] = ((U_i - p) / (U_i - T_i))^t     when T_i < p <= U_i
%[text] = 0                                when p > U_i
%[text]
%[text] Here, L_i is the lower limit, T_i is the ideal target, U_i is the upper limit, and s and t control the curvature.
%[text]
%[text] The composite desirability D is calculated by combining all individual d_i.
%[text]
%[text] D = (d_1 * d_2 * ... * d_k)^(1/k)    (geometric mean)
%[text]
%[text] When D = 1, all properties are ideal simultaneously.
%[text] When D = 0, at least one property is unacceptable.
%[text] The geometric mean imposes a stronger penalty than the arithmetic mean for a single very low score.
%[text]
%[text] Definition of property desirability (this exercise):
%[text]
%[text] LogP:   L = 0.5, T = 2.0, U = 4.0   (s=t=1, linear ramp)
%[text] TPSA:   One-sided upper limit (lower limit always acceptable):
%[text] d_TPSA = 1              when TPSA <= 60
%[text] = (90 - TPSA)/30 when 60 < TPSA <= 90
%[text] = 0              when TPSA > 90
%[text] MW:     L = 150, T = 350, U = 500   (s=1, t=1)
%[text]
logSection("A10", "Section 3: Desirability Function — From Properties to Scores", "Analytics L3");
%[text] --- Individual Desirability Functions (Vectorized) ---
d_logp = @(p) max(0, min(1, ...
    (p >= 0.5 & p <= 2.0) .* ((p - 0.5) ./ 1.5) + ...
    (p >  2.0 & p <= 4.0) .* ((4.0 - p) ./ 2.0) + ...
    (p == 2.0)));
%[text] Note: The above expression may overlap at p=2. Use piecewise.
d_logp = @(p) ...
    ((p >= 0.5) & (p <  2.0)) .* ((p - 0.5) ./ 1.5) + ...
    ((p >= 2.0) & (p <= 4.0)) .* ((4.0 - p) ./ 2.0) + ...
    (p == 2.0);

d_tpsa = @(p) ...
    (p <= 60) .* 1.0 + ...
    (p >  60 & p <= 90) .* ((90 - p) ./ 30) + ...
    (p >  90) .* 0;

d_mw   = @(p) ...
    ((p >= 150) & (p <  350)) .* ((p - 150) ./ 200) + ...
    ((p >= 350) & (p <= 500)) .* ((500 - p) ./ 150) + ...
    (p == 350);

%[text] Composite Desirability (Geometric Mean of 3 Components)
D_composite = @(logp, tpsa, mw) ...
    (d_logp(logp) .* d_tpsa(tpsa) .* d_mw(mw)) .^ (1/3);

%[text] Calculate the library scores.
d1 = d_logp(logp_vec);
d2 = d_tpsa(tpsa_vec);
d3 = d_mw(mw_vec);
D  = (d1 .* d2 .* d3) .^ (1/3);

%[text] Calculate the lead compound score.
D_lead = D_composite(d_lead.LogP, d_lead.TPSA, d_lead.MolWt);
logInfo("Lead (%s) Desirability: D=%.3f (d_LogP=%.2f, d_TPSA=%.2f, d_MW=%.2f)", ...
    LEAD_NAME, D_lead, d_logp(d_lead.LogP), d_tpsa(d_lead.TPSA), d_mw(d_lead.MolWt));

%[text] Display library statistics.
logInfo("Library Desirability: Mean=%.3f, Max=%.3f, Median=%.3f", ...
    mean(D), max(D), median(D));

%[text] Display the top 10 library molecules by desirability.
[D_sorted, D_rank] = sort(D, "descend");
topN_des = min(10, nValid);
logInfo("Top %d Library Molecules by Composite Desirability:", topN_des);
for k = 1:topN_des
    idx = D_rank(k);
    logInfo("  %2d. D=%.3f  LogP=%5.2f  TPSA=%5.1f  MW=%5.1f  %s", ...
        k, D_sorted(k), logp_vec(idx), tpsa_vec(idx), mw_vec(idx), libNames(idx));
end

%[text] --- Desirability Function Plot ---
p_range_logp = linspace(-1, 6, 300);
p_range_tpsa = linspace(0, 150, 300);
p_range_mw   = linspace(100, 700, 300);

figure("Name", "A10 Sec3: Desirability Function");
subplot(1, 3, 1);
plot(p_range_logp, d_logp(p_range_logp), "b-", LineWidth=2);
xline(TARGET_LOGP, "k--", Label="Target", LabelVerticalAlignment="bottom");
xline(d_lead.LogP, "r:", LineWidth=1.5, Label=LEAD_NAME);
xlabel("LogP");  ylabel("d_{LogP}");
title("LogP Desirability");  ylim([0, 1.1]);  grid on;

subplot(1, 3, 2);
plot(p_range_tpsa, d_tpsa(p_range_tpsa), "g-", LineWidth=2);
xline(TARGET_TPSA, "k--", Label="Target", LabelVerticalAlignment="bottom");
xline(d_lead.TPSA, "r:", LineWidth=1.5, Label=LEAD_NAME);
xlabel("TPSA (A^2)");  ylabel("d_{TPSA}");
title("TPSA Desirability");  ylim([0, 1.1]);  grid on;

subplot(1, 3, 3);
plot(p_range_mw, d_mw(p_range_mw), "m-", LineWidth=2);
xline(TARGET_MW, "k--", Label="Target", LabelVerticalAlignment="bottom");
xline(d_lead.MolWt, "r:", LineWidth=1.5, Label=LEAD_NAME);
xlabel("MW (g/mol)");  ylabel("d_{MW}");
title("MW Desirability");  ylim([0, 1.1]);  grid on;

sgtitle("Derringer-Suich Desirability Function");

%[text] **💡 Observation Point 3**
%[text] The geometric mean D = (d1*d2*d3)^(1/3) penalizes a single poor score.
%[text] Replace with the arithmetic mean and observe the results.
%[text] D_arith = (d1 + d2 + d3) / 3;
%[text] Check how the top 10 ranking changes and consider which formulation is more conservative for compounds with low properties.
%[text] If you want to favor CNS drugs by changing LogP desirability (TPSA < 60, LogP 1-3):
%[text] Adjust the target of d_tpsa to 45 and the upper limit to 60.
%[text] Check how many library molecules remain under this condition.
%[text] Find the molecule with the highest d_LogP score and the lowest d_TPSA score.
%[text] Consider the chemical features explaining the low TPSA score.
% ... (Try writing code here)
%%
%[text] ## Section 4: Goal Attainment Optimization (fgoalattain)
%[text]
%[text] Candidates were ranked by desirability score. Next, we use fgoalattain to find the optimal blend that simultaneously meets three property goals.
%[text] Let's understand the role of the attainment factor gamma.
%[text]
%[text] ### Concept: Goal Attainment and Attainment Factor
%[text]
%[text] The goal attainment method (Charnes & Cooper 1977) solves the following problem:
%[text]
%[text] min  gamma
%[text] s.t. f_j(x) - gamma * w_j <= goal_j    for all j = 1..m
%[text] x in feasible set X
%[text]
%[text] Where:
%[text] x = decision variables
%[text] f_j(x) = j-th objective function
%[text] goal_j = desired target for the j-th objective
%[text] w_j = weight (positive, scales slack)
%[text] gamma = attainment factor (optimization variable)
%[text]
%[text] Interpretation of gamma:
%[text] gamma <= 0  =>  All goals are simultaneously achieved
%[text] gamma >  0  =>  Some goals are overachieved, the "least" goal is gamma * w units unmet
%[text]
%[text] Problem formulation (this exercise):
%[text]
%[text] Represent a "virtual compound" as a convex mixture of library molecules.
%[text] Decision variables x in R^N (N = library size) satisfy:
%[text]
%[text] x >= 0, sum(x) = 1    (mixture constraint)
%[text]
%[text] The blended properties are:
%[text]
%[text] f(x) = [sum_i x_i * LogP_i, sum_i x_i * TPSA_i, sum_i x_i * MW_i]
%[text] = PropMat' * x     (3 x 1 vector)
%[text]
%[text] Where PropMat is an N x 3 property matrix.
%[text]
%[text] fgoalattain finds x* such that PropMat'*x* is closest to the goal vector [TARGET_LOGP, TARGET_TPSA, TARGET_MW] in a minimax sense.
%[text]
%[text] Real molecule closest to the optimal blend point:
%[text] After computing the optimal blend point p* = PropMat' * x*, find the library molecule with the smallest Euclidean distance to p* in the normalized property space (zero mean, unit variance).
%[text]
%[text] Note: The convex hull of the library may not contain the target point.
%[text] The attainment factor gamma quantifies how far the target is outside the convex hull.
logSection("A10", "Section 4: Goal Attainment Optimization (fgoalattain)", "Analytics L3");
n = nValid;
propGoal   = [TARGET_LOGP; TARGET_TPSA; TARGET_MW];
propWeight = [1.0; 1.0; 1.0];   % Equal weight for all objectives

%[text] Mixture constraint: sum(x) = 1, x >= 0
Aeq_mix = ones(1, n);
beq_mix = 1;
lb_mix  = zeros(n, 1);
ub_mix  = ones(n, 1);

%[text] Objective function: Blended property vector (linear in x)
fun_blend = @(x) propMat' * x(:);   % returns 3 x 1

if hasOptTbx
    % --- OPTIMIZATION TOOLBOX path ---
    x0_mix   = ones(n, 1) / n;   % Start from uniform mixture
    opts_ga  = optimoptions("fgoalattain", ...
        Display="off", MaxFunctionEvaluations=n*50, ...
        OptimalityTolerance=1e-6, ConstraintTolerance=1e-6);

    logInfo("Running fgoalattain (%d variables, 3 objectives)...", n);
    tic;
    [x_opt, fval_blend, attainfactor, exitflag] = fgoalattain( ...
        fun_blend, x0_mix, propGoal, propWeight, ...
        [], [], Aeq_mix, beq_mix, lb_mix, ub_mix, [], opts_ga);
    t_opt = toc;

    logInfo("fgoalattain completed in %.2f seconds (exitflag=%d).", t_opt, exitflag);
    logInfo("Attainment factor gamma = %.4f", attainfactor);
    if attainfactor <= 0
        logInfo("  => All property goals achieved within the convex hull.");
    else
        logInfo("  => Goals partially unmet -- target is outside the convex hull.");
        logInfo("     Nearest blend point: LogP=%.2f, TPSA=%.1f, MW=%.1f", ...
            fval_blend(1), fval_blend(2), fval_blend(3));
    end

else
    % --- Manual fallback: Nearest molecule by weighted Euclidean distance ---
    logWarn("fgoalattain unavailable -- using weighted nearest neighbor fallback.");
    propSigma = std(propMat);           % normalisation scale
    propGoalN = propGoal(:)' ./ propSigma;
    propMatN  = propMat ./ propSigma;
    dist_all  = sqrt(sum((propMatN - propGoalN).^2, 2));
    [~, nn_idx] = min(dist_all);
    fval_blend = propMat(nn_idx, :)';  % 3 x 1
    attainfactor = max((fval_blend - propGoal) ./ propWeight);
    logInfo("Nearest: %s  (D_composite=%.3f)", libNames(nn_idx), D(nn_idx));
end

%[text] --- Find the real molecule closest to the optimal blend point ---
%[text] Normalize by library standard deviation (unit variance distance)
propSigma = std(propMat);
propMatN  = propMat ./ propSigma;       % N x 3 normalised
fvalN     = fval_blend(:)' ./ propSigma;  % 1 x 3

dist_to_blend = sqrt(sum((propMatN - fvalN).^2, 2));
[~, nearestRank] = sort(dist_to_blend);

%[text] Display top 5 nearest real molecules
TOP_GA = 5;
logInfo("Top %d real molecules nearest to the optimal blend point:", TOP_GA);
logInfo("  (Optimal blend: LogP=%.2f, TPSA=%.1f, MW=%.1f)", ...
    fval_blend(1), fval_blend(2), fval_blend(3));
for k = 1:min(TOP_GA, nValid)
    idx = nearestRank(k);
    logInfo("  %d. Distance=%.3f  D=%.3f  LogP=%5.2f  TPSA=%5.1f  MW=%5.1f  %s", ...
        k, dist_to_blend(idx), D(idx), logp_vec(idx), tpsa_vec(idx), ...
        mw_vec(idx), libNames(idx));
end

%[text] --- Radar (Spider) Chart: Lead vs Top Desirability Candidates ---
%[text] Note: When gamma << 0, the fgoalattain blend nearest neighbor (nearestRank(1)) reaches LogP~-3.5 (as fgoalattain minimizes the maximum gap of unmet goals).
%[text] For meaningful comparison, use the top desirability molecule instead.
[~, topGaIdx]  = max(D);
radarProps = ["LogP", "TPSA", "MW"];
radarLead  = [d_lead.LogP, d_lead.TPSA, d_lead.MolWt];
radarCand  = [logp_vec(topGaIdx), tpsa_vec(topGaIdx), mw_vec(topGaIdx)];
radarGoal  = [TARGET_LOGP, TARGET_TPSA, TARGET_MW];

%[text] Normalize to [0,1] using fixed drug-like space boundaries.
%[text] Using library min/max stretches the chart due to extreme outliers (MW=924, TPSA=319). Fixed boundaries keep the chart in a human-readable region for comparing 3 compounds.
radarMin   = [-1,   0, 100];   % [LogP_lo, TPSA_lo, MW_lo]
radarMax   = [ 5, 120, 500];   % [LogP_hi, TPSA_hi, MW_hi]
radarNorm  = @(v) (v - radarMin) ./ (radarMax - radarMin);
radarNorm  = @(v) max(0, min(1, (v - radarMin) ./ (radarMax - radarMin)));  % clip to [0,1]

figure("Name", "A10 Sec4: Property Profile -- Lead vs Candidate vs Goal");
theta    = linspace(0, 2*pi, numel(radarProps) + 1);
nLead    = radarNorm(radarLead);   % 1x3
nCand    = radarNorm(radarCand);   % 1x3
nGoal    = radarNorm(radarGoal);   % 1x3
rLead    = [nLead, nLead(1)];      % 1x4  (close the polygon)
rCand    = [nCand, nCand(1)];      % 1x4
rGoal    = [nGoal, nGoal(1)];      % 1x4

polarplot(theta, rLead, "r-o", LineWidth=2, DisplayName=LEAD_NAME);
hold on;
polarplot(theta, rCand, "b-s", LineWidth=2, DisplayName=libNames(topGaIdx));
polarplot(theta, rGoal, "k--^", LineWidth=1.5, DisplayName="Goal");
legend(Location="southoutside");
title(sprintf("Property Profile (Normalized)\n%s vs %s vs Goal", ...
    LEAD_NAME, libNames(topGaIdx)));
pax = gca;
pax.ThetaTick      = [0, 120, 240];
pax.ThetaTickLabel = {"LogP", "TPSA (A^2)", "MW (g/mol)"};

%[text] **💡 Observation Point 4**
%[text] Let's confirm what the attainment factor gamma is.
%[text] If gamma <= 0, all goals are achievable; if gamma > 0, some are unmet.
%[text] Change TARGET_MW from 350 to 250 and rerun to observe changes in gamma.
%[text] Observe how gamma changes and which goal becomes the hardest to achieve.
%[text] Note: When gamma << 0 (e.g., -5.5), fgoalattain drives the blend to fall below all goals simultaneously. This means the "optimal blend" point may be far from the goals (e.g., LogP = 2.0 instead of -3.5). The top 5 nearest molecules to that blend point often have D = 0 (outside the desirability zone).
%[text] To find the nearest feasible compound, use the top D ranking (BETAXOLOL) instead of the fgoalattain blend nearest neighbor.
%[text] Check what happens when you prioritize LogP over TPSA by changing weights:
%[text] propWeight = [0.5; 2.0; 1.0];  % Double weight for TPSA
%[text] This tells fgoalattain that "missing the TPSA goal costs twice as much."
%[text] See which molecule moves to the top rank.
%[text] Consider what the interpretation of "mixture" physically means.
%[text] x_opt assigns weights to each library molecule. If x_opt has three non-zero entries with weights [0.3, 0.4, 0.3], it means the "ideal analog" combines structural features of all three molecules.
%[text] In actual lead optimization, it guides which substituents to merge.
% ... (Try writing code here)
%%
%[text] ## Section 5: Pareto Front -- LogP vs TPSA Trade-off
%[text]
%[text] The results of the goal achievement optimization have been obtained. Next, let's plot the Pareto front of LogP (logarithm of the partition coefficient) and TPSA (topological polar surface area) to understand the overall trade-off.
%[text] Determine which compounds are Pareto optimal and which are dominated.
%[text]
%[text] ### Concept: Pareto Optimality and Epsilon Constraint Method
%[text]
%[text] A solution x is Pareto optimal (non-dominated) if there is no other feasible solution that is better in all objectives simultaneously.
%[text]
%[text] In a 2-objective context (minimizing |LogP - 2.0| and TPSA):
%[text] Molecule A dominates molecule B if:
%[text] A is lower than B in both |LogP - target| and TPSA.
%[text]
%[text] The Pareto front is the set of all non-dominated solutions.
%[text] It represents the true trade-off between the two objectives.
%[text]
%[text] EPSILON Constraint Method (Cohon & Marks 1975):
%[text] Instead of combining objectives into one (weighted sum), optimize one objective while constraining the other:
%[text]
%[text] min  TPSA(x) = tpsa' * x
%[text] s.t. LogP(x) = logp' * x >= epsilon      (LogP lower bound)
%[text] sum(x) = 1, x >= 0
%[text]
%[text] Sweep epsilon from LogP_min to LogP_max of the library to trace the complete Pareto front. Solve the linear program (LP) at each epsilon value using linprog (Optimization Toolbox).
%[text]
%[text] Fallback (without Optimization Toolbox):
%[text] Verify non-dominance and extract the empirical Pareto front directly from the library:
%[text] Molecule i is non-dominated if there is no molecule j such that TPSA_j < TPSA_i and |LogP_j - 2| < |LogP_i - 2|.
%[text]
%[text] Definition of objective functions concerning mixed variable x
%[text] obj1 = |LogP(x) - TARGET_LOGP|  (one-sided use: LogP(x) >= eps)
%[text] obj2 = TPSA(x)                  (minimize TPSA)
%[text]
%[text] Note on linprog mixed approach:
%[text] The epsilon constraint method of linprog for 200 molecule mixtures is mathematically valid but may misleadingly present a Pareto front in practice.
%[text] Convex combinations always assign ~100% weight to the lowest TPSA molecule (TPSA=3.2), making the optimal TPSA trivially near zero for any LogP. The front shrinks to a flat line at the plot's bottom, diverging from actual drugs.
%[text] Instead, use the empirical Pareto front (non-dominated real molecules).
%[text]
%[text] --- Empirical Pareto Front: Non-dominated Set from Library ---
logSection("A10", "Section 5: Pareto Front — LogP vs TPSA Trade-off", "Analytics L3");
obj1 = abs(logp_vec - TARGET_LOGP);   % minimise deviation from target LogP
obj2 = tpsa_vec;                        % minimise TPSA

nondominatedMask = true(nValid, 1);
for ii = 1:nValid
    for jj = 1:nValid
        if jj ~= ii && obj1(jj) <= obj1(ii) && obj2(jj) <= obj2(ii) && ...
           (obj1(jj) < obj1(ii) || obj2(jj) < obj2(ii))
            nondominatedMask(ii) = false;
            break;
        end
    end
end
pIdx        = find(nondominatedMask);
[~, sortP]  = sort(logp_vec(pIdx));
pIdx        = pIdx(sortP);
pareto_logp = logp_vec(pIdx);
pareto_tpsa = tpsa_vec(pIdx);
logInfo("Empirical Pareto Front: %d non-dominated molecules.", numel(pIdx));

if hasOptTbx
    % linprog concept demo (NOT used for the plot -- see NOTE above).
    % Uncomment to observe the mixture-front collapse:
    %   N_LP = 20;
    %   eps_sw = linspace(min(logp_vec)+0.1, max(logp_vec)-0.1, N_LP);
    %   opt_lp = optimoptions("linprog", Display="off");
    %   for k = 1:N_LP
    %       [~, tv, fk] = linprog(tpsa_vec, -logp_vec', -eps_sw(k), ...
    %           ones(1,n), 1, zeros(n,1), ones(n,1), opt_lp);
    %       if fk==1, fprintf("eps=%.2f TPSA_mix=%.1f\n",eps_sw(k),tv); end
    %   end
end

%[text] --- Pareto Front Plot ---
figure("Name", "A10 Sec5: LogP vs TPSA Pareto Front");

%[text] Library Background
scatter(logp_vec, tpsa_vec, 25, D, "filled", ...
    MarkerFaceAlpha=0.5, DisplayName="Library (colour = D)");
colormap("cool");
cb2 = colorbar();
cb2.Label.String = "Composite Desirability D";
hold on;

%[text] Pareto Front Curve
plot(pareto_logp, pareto_tpsa, "k-o", LineWidth=2, ...
    MarkerFaceColor="k", MarkerSize=5, DisplayName="Pareto front");

%[text] Ideal (Target) Point
plot(TARGET_LOGP, TARGET_TPSA, "rp", MarkerSize=18, ...
    LineWidth=2, DisplayName="Ideal target");

%[text] Lead Compound
plot(d_lead.LogP, d_lead.TPSA, "rv", MarkerSize=12, ...
    MarkerFaceColor="r", DisplayName=LEAD_NAME);

%[text] Veber TPSA Threshold
yline(90, "g--", LineWidth=1.5, Label="Veber TPSA=90", DisplayName="Veber TPSA=90");
yline(60, "b--", LineWidth=1.5, Label="CNS TPSA=60",  DisplayName="CNS TPSA=60");

xlabel("LogP");
ylabel("TPSA (A^2)");
title("LogP vs TPSA: Pareto Front (epsilon constraint)");
legend(Location="northeast");
grid on;

%[text] **💡 Observation Point 5**
%[text] Observe the Pareto front curve and see how TPSA changes as you move from left to right (increasing LogP). Is the trade-off monotonic?
%[text] What does this indicate about the diversity of the library?
%[text] The ideal point [LogP=2, TPSA=70] is marked with a star.
%[text] Check if this point is on the Pareto front, and if not, whether it is above or below.
%[text] (Points below the front are infeasible, indicating that the library molecules or mixtures cannot achieve that combination.)
%[text] Try changing the Pareto sweep to use MW (molecular weight) as the second objective.
%[text] Replace tpsa_vec with mw_vec and TARGET_TPSA with TARGET_MW.
%[text] Rerun linprog and see if the LogP-MW trade-off is stronger or weaker than the LogP-TPSA trade-off. (Hint: Check the Pearson correlation first.)
%[text] The epsilon constraint method samples the Pareto front uniformly in LogP space.
%[text] Points in the "flat" region of the front indicate little gain from accepting higher LogP.
%[text] See if you can identify the "knee point" -- the Pareto point where increasing LogP results in the maximum TPSA reduction.
% ... (Try writing code here)
%%
%[text] ## Section 6: Lead Optimization Summary and Candidate Visualization
%[text]
%[text] The Pareto front has been obtained. We select top candidates by combining desirability, goal distance, and Pareto optimality.
%[text] Display the final candidate table and structures to complete the exercise.
%[text]
%[text] ### Concept: Candidate Selection — Combining Desirability and Proximity
%[text]
%[text] The final selection combines the following three criteria:
%[text] (a)  Composite desirability D (see Section 3)
%[text] (b)  Proximity to optimal goal values (see Section 4)
%[text] (c)  Pareto non-dominance of LogP vs TPSA (see Section 5)
%[text]
%[text] In actual projects, the following additional filters are applied:
%[text] - Synthetic Accessibility Score (SAS)
%[text] - hERG cardiotoxicity prediction
%[text] - PAINS / structural alert filters (see S03)
%[text] - ADMET modeling (Absorption, Distribution, Metabolism, Excretion, Toxicity)
%[text] Example: pkCSM or SwissADME web servers
%[text]
%[text] --- Construct Final Candidate Table ---
logSection("A10", "Section 6: Lead Optimization Summary and Candidate Visualization", "Analytics L3");
TOP_FINAL = 10;

%[text] Ranking based on composite desirability
[D_sorted_f, rank_by_D] = sort(D, "descend");

candidateTbl = table( ...
    (1:TOP_FINAL)', ...
    libNames(rank_by_D(1:TOP_FINAL)), ...
    round(logp_vec(rank_by_D(1:TOP_FINAL)), 2), ...
    round(tpsa_vec(rank_by_D(1:TOP_FINAL)), 1), ...
    round(mw_vec(rank_by_D(1:TOP_FINAL)), 1), ...
    round(D_sorted_f(1:TOP_FINAL), 3), ...
    round(dist_to_blend(rank_by_D(1:TOP_FINAL)), 3), ...
    VariableNames=["Rank", "Name", "LogP", "TPSA_A2", "MW_gmol", ...
                   "Desirability", "Dist2Goal"]);

logInfo("--- Top %d Candidates by Composite Desirability ---", TOP_FINAL);
disp(candidateTbl);

%[text] --- Comparison of Lead vs Top Candidate ---
bestIdx = rank_by_D(1);

logInfo("--- Lead vs Best Candidate ---");
logInfo("               %-25s  %-25s  Target", LEAD_NAME, libNames(bestIdx));
logInfo("  LogP       : %8.2f                      %8.2f          %.1f", ...
    d_lead.LogP, logp_vec(bestIdx), TARGET_LOGP);
logInfo("  TPSA (A^2) : %8.1f                      %8.1f          %.1f", ...
    d_lead.TPSA, tpsa_vec(bestIdx), TARGET_TPSA);
logInfo("  MW (g/mol) : %8.1f                      %8.1f          %.1f", ...
    d_lead.MolWt, mw_vec(bestIdx), TARGET_MW);
logInfo("  Desirability: %7.3f                      %7.3f          1.000", ...
    D_lead, D(bestIdx));

%[text] --- Visualize Top 5 Candidates ---
TOP_VIZ = 5;
vizMols  = cell(1, TOP_VIZ);
vizNames = strings(1, TOP_VIZ);
for k = 1:TOP_VIZ
    idx          = rank_by_D(k);
    vizMols{k}   = mols{idx};
    vizNames(k)  = sprintf("%s\nD=%.3f LogP=%.2f MW=%.0f", ...
        libNames(idx), D(idx), logp_vec(idx), mw_vec(idx));
end

logInfo("Rendering Top %d Candidates...", TOP_VIZ);
for k = 1:TOP_VIZ
    try
        figure();  % create a fresh figure so draw2d does not reuse the Pareto axes
        fig = emk.viz.draw2d(vizMols{k}, Title=vizNames(k));
        set(fig, "Name", sprintf("A10 Sec6: Candidate %d -- %s", ...
            k, libNames(rank_by_D(k))));
    catch ME
        logWarn("draw2d failed for Candidate %d: %s", k, ME.message);
    end
end

%[text] --- Summary Bar Graph: Desirability Breakdown ---
barData = [d_logp(d_lead.LogP),  d_tpsa(d_lead.TPSA),  d_mw(d_lead.MolWt);   ...
           d_logp(logp_vec(rank_by_D(1:TOP_VIZ-1))), ...
           d_tpsa(tpsa_vec(rank_by_D(1:TOP_VIZ-1))), ...
           d_mw(mw_vec(rank_by_D(1:TOP_VIZ-1)))];

barLabels = [LEAD_NAME; libNames(rank_by_D(1:TOP_VIZ-1))];
%[text] Shorten long names for display
for k = 1:numel(barLabels)
    if strlength(barLabels(k)) > 18
        barLabels(k) = extractBefore(barLabels(k), 18) + "...";
    end
end

figure("Name", "A10 Sec6: Desirability Breakdown -- Lead vs Top Candidates");
b = bar(barData, "grouped");
b(1).FaceColor = [0.2 0.6 0.9];   % d_LogP
b(2).FaceColor = [0.3 0.8 0.4];   % d_TPSA
b(3).FaceColor = [0.9 0.5 0.2];   % d_MW
legend(["d_{LogP}", "d_{TPSA}", "d_{MW}"], Location="northeast");
xticklabels(barLabels);
xtickangle(20);
ylabel("Individual Desirability Score");
title("Desirability Breakdown: Lead vs Top Candidates");
ylim([0, 1.05]);
yline(1.0, "k--", HandleVisibility="off");
grid on;

%[text] **💡 Observation Point 6**
%[text] Check the molecules with the highest d_TPSA score and the highest d_LogP score. Are they the same molecule?
%[text] What does this indicate about the challenges of multi-objective selection?
%[text] Try adding the PAINS filter from S03 to the candidate list:
%[text] pains_smarts = emk.filter.loadPainsSmarts();  % (if available)
%[text] hasPains = cellfun(@(m) emk.mol.hasSubstruct(m, pains_smarts), mols);
%[text] Remove PAINS-positive candidates from the top 10 list.
%[text] How many remain?
%[text] In actual lead optimization, synthetic accessibility is also considered.
%[text] Try manually ranking top candidates based on synthetic complexity.
%[text] Does high desirability correlate with easy synthesis?
%[text] Choose your lead compound and property target set.
%[text] Re-run the full pipeline from Section 1.
%[text] Recommended Leads:
%[text] Caffeine: "Cn1cnc2c1c(=O)n(C)c(=O)n2C"   — CNS Stimulant
%[text] Ibuprofen: "CC(C)Cc1ccc(cc1)C(C)C(=O)O"   — NSAID
%[text] Metformin: "CN(C)C(=N)NC(=N)N"              — Antidiabetic
%[text]
%[text] **Summary**
%[text] We applied the Derringer-Suich desirability function to convert pharmacological intuition into dimensionless scores and used fgoalattain for goal attainment optimization to simultaneously minimize all property goals.
%[text] Tracked the Pareto front of LogP vs TPSA using the epsilon constraint method, excluding inferior candidates.
logInfo("A10: Completed.");

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
