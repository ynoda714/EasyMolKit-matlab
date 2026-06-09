%[text] # A06 Answer: Dose-Response Curve Fitting
%[text] Reference answer for the "Try It" exercise in a06_dose_response.m.
%[text]
%[text] Note: This file uses only the Curve Fitting Toolbox -- no RDKit / Python
%[text]       dependency. `emk.setup.initPython()` is not required.
%[text]
%[text] The expression and parameter name arguments for `fittype()` and `fitoptions()`
%[text] must all be char literals (not double-quoted strings).
addpath(genpath("src"));
logInfo("A06 Answer: Setup complete (Python not required)");

%[text] Redefine shared constants and helper functions (for standalone use)
TRUE_EMAX  = 100;
TRUE_EC50  = 1.0;
TRUE_N     = 1.5;
TRUE_EMIN  = 5;
NOISE_STD  = 8.0;
N_CONC     = 8;
rng(42);

hillModel = fittype( ...
    'Emin + (Emax - Emin) ./ (1 + (EC50 ./ x).^n)', ...
    'independent', 'x', ...
    'coefficients', {'Emax', 'n', 'EC50', 'Emin'});

hillFit_ = @(conc, resp) fit(conc(:), resp(:), hillModel, ...
    fitoptions(hillModel, ...
        'Lower',       [50,   0.1,  0.001, -10], ...
        'Upper',       [150,  10,   100,    50], ...
        'StartPoint',  [90,   1.0,  1.0,    10], ...
        'TolFun',      1e-8, ...
        'MaxIter',     2000, ...
        'Robust',      'off'));

conc = logspace(-3, 2, N_CONC)';   % 8 concentrations 0.001 to 100 uM
trueResp = @(c) TRUE_EMIN + (TRUE_EMAX - TRUE_EMIN) ./ ...
    (1 + (TRUE_EC50 ./ c).^TRUE_N);
resp = trueResp(conc) + NOISE_STD * randn(N_CONC, 1);
%%
%[text] ## Let's Try 1: Biology of Hill Slope; y at x=10*EC50 and x=0.1*EC50; Proof of Midpoint

logInfo("=== Let's Try 1: Hill Slope n Value ===");
logInfo("n = 1  : simple Michaelis-Menten (enzyme kinetics)");
logInfo("n = 2.8: Hemoglobin O2 binding (cooperative)");
logInfo("n = 1.5: Weak cooperativity (our true model)");

%[text] Plot the comparison of sigmoid shapes by Hill slope (n)
concFine = logspace(-2, 2, 200)';
figure("Name", "A06 Hill Slope Effect");
hold on;
for nVal = [0.5, 1, 1.5, 2, 4]
    yFine = TRUE_EMIN + (TRUE_EMAX - TRUE_EMIN) ./ (1 + (TRUE_EC50./concFine).^nVal);
    semilogx(concFine, yFine, LineWidth=1.5, DisplayName=sprintf("n=%.1f", nVal));
end
hold off;
xlabel("Concentration (uM)"); ylabel("Response (%)");
title("Shape of Hill Equation vs Hill Slope n");
legend(Location="southeast"); grid on;

%[text] Calculate the response at x=10*EC50 and x=0.1*EC50
y10  = TRUE_EMIN + (TRUE_EMAX-TRUE_EMIN) / (1 + (TRUE_EC50/(10*TRUE_EC50))^TRUE_N);
y01  = TRUE_EMIN + (TRUE_EMAX-TRUE_EMIN) / (1 + (TRUE_EC50/(0.1*TRUE_EC50))^TRUE_N);
logInfo("y at x=10*EC50 : %.1f%%  (expected ~%d%% for Emax=%d)", ...
    y10, 100*(TRUE_EMAX-TRUE_EMIN)/(TRUE_EMAX-TRUE_EMIN+1)*TRUE_EMIN, TRUE_EMAX);
logInfo("y at x=0.1*EC50: %.1f%%  (close to Emin=%.1f)", y01, TRUE_EMIN);

%[text] Algebraic proof of midpoint at x=EC50
%[text] Hill Equation: y = Emin + (Emax-Emin)/(1+(EC50/x)^n)
%[text] When x=EC50: Denominator = 1+(EC50/EC50)^n = 1+1 = 2
%[text] => y = Emin + (Emax-Emin)/2 = (Emax+Emin)/2
midpoint_analytical = (TRUE_EMAX + TRUE_EMIN) / 2;
midpoint_formula    = TRUE_EMIN + (TRUE_EMAX - TRUE_EMIN) / 2;
logInfo("Midpoint proof: Emin+(Emax-Emin)/2 = %.1f  ==  (Emax+Emin)/2 = %.1f", ...
    midpoint_formula, midpoint_analytical);

%[text] Answer: A large n (steep sigmoid) means the response switches sharply over a narrow concentration range.
%[text]    Hemoglobin (n~2.8) allows efficient O2 binding in the lungs (high pO2) and
%[text]    release in tissues (low pO2).
%[text]    At x=10*EC50 (very high dose), y approaches Emax but never fully reaches it.
%[text]    At x=0.1*EC50 (below effective dose), y ≈ Emin.
%[text]    Algebraically: When x=EC50, (EC50/x)^n = 1, denominator = 2, and
%[text]    y = Emin + (Emax-Emin)/2 = (Emax+Emin)/2 -- this is exactly the midpoint.
%%
%[text] ## Let's Try 2: Effect of Noise; Triplicate Design Using SEM

logInfo("=== Let's Try 2: Noise and Measurement Design ===");
for noiseLevel = [0, 4, 8, 15, 25]
    rng(42);
    respN = trueResp(conc) + noiseLevel * randn(N_CONC, 1);
    try
        fN = hillFit_(conc, respN);
        logInfo("noiseStd=%-3d  -> fitted EC50=%.3f (true=%.3f)  error=%+.3f", ...
            noiseLevel, fN.EC50, TRUE_EC50, fN.EC50 - TRUE_EC50);
    catch ME
        logWarn("noiseStd=%-3d  -> fit failed: %s", noiseLevel, ME.message);
    end
end

%[text] Triplicate design (3 measurements per concentration)
rng(42);
nReps = 3;
respTri = trueResp(conc) * ones(1, nReps) + NOISE_STD * randn(N_CONC, nReps);
semPerConc = std(respTri, 0, 2) / sqrt(nReps);
logInfo("Triplicate SEM (mean across concentrations): %.2f", mean(semPerConc));
logInfo("Single point noise standard deviation: %.2f", NOISE_STD);

%[text] Fit using the triplicate means
respMean = mean(respTri, 2);
fTri = hillFit_(conc, respMean);
logInfo("Fit on triplicate means: EC50=%.3f (true=%.3f)", fTri.EC50, TRUE_EC50);

%[text] Answer: As noise increases, the confidence interval of EC50 widens. At noiseStd=25,
%[text]    the fit may not converge, or biologically unreasonable parameters
%[text]    (EC50 hitting boundary values) may be returned.
%[text]    Triplicate measurements reduce effective noise by sqrt(3) ≈ 1.73 times.
%[text]    SEM = std/sqrt(n) quantifies the uncertainty at each concentration point.
%[text]    Increasing the number of measurements improves fit quality but triples reagent costs.
%%
%[text] ## Let's Try 3: Alternative Starting Points; Risks of Unconstrained Fit

logInfo("=== Let's Try 3: Convergence Sensitivity ===");
startPoints = {[90, 1.0, 1.0, 10], [60, 2.0, 5.0, 0], [99, 0.5, 0.1, 1]};
for k = 1:numel(startPoints)
    sp = startPoints{k};
    try
        fk = fit(conc, resp, hillModel, ...
            fitoptions(hillModel, ...
                'Lower',      [50,  0.1, 0.001, -10], ...
                'Upper',      [150, 10,  100,    50], ...
                'StartPoint', sp, ...
                'TolFun',     1e-8, 'MaxIter', 2000));
        logInfo("StartPoint [Emax=%.0f n=%.1f EC50=%.1f Emin=%.0f] -> EC50=%.3f", ...
            sp(1), sp(2), sp(3), sp(4), fk.EC50);
    catch ME
        logWarn("StartPoint failed: %s", ME.message);
    end
end

%[text] Check the risks of unconstrained fit
try
    fUC = fit(conc, resp, hillModel, ...
        fitoptions(hillModel, 'MaxIter', 1000));
    logInfo("Unconstrained fit: EC50=%.3f  n=%.3f", fUC.EC50, fUC.n);
catch ME
    logWarn("Unconstrained fit failed: %s", ME.message);
end

%[text] Answer: Appropriate starting points close to the true value will ensure convergence.
%[text]    Inappropriate starting points (e.g., EC50=0.1 when the true value is 1.0) may converge
%[text]    to a local minimum or may not converge within MaxIter.
%[text]    Unconstrained fit may return negative EC50 (physically meaningless) or n > 10 (no known
%[text]    biological mechanism produces such a steep sigmoid). Always provide physiologically
%[text]    meaningful bounds to fittype/fit.
%[text]    Derivation of Hill equation: Assuming rapid equilibrium of ligand-receptor
%[text]    EC50^n = [L]^n * (1 - y)/y (y is the fractional occupancy)
%[text]    => y = [L]^n / ([L]^n + EC50^n)
%%
%[text] ## Let's Try 4: CI Width with N_CONC=5; Zero Noise CI; Meaning of predint

logInfo("=== Let's Try 4: Confidence Interval Width ===");
for nPoints = [5, 8, 12, 20]
    rng(42);
    cTmp = logspace(-3, 2, nPoints)';
    rTmp = trueResp(cTmp) + NOISE_STD * randn(nPoints, 1);
    try
        fTmp = hillFit_(cTmp, rTmp);
        ciTmp = confint(fTmp, 0.95);
        width = ciTmp(2,3) - ciTmp(1,3);   % column 3 = EC50
        logInfo("N_CONC=%2d -> EC50 CI Width = %.3f", nPoints, width);
    catch
        logWarn("N_CONC=%2d -> fit failed", nPoints);
    end
end

%[text] Check the case of zero noise
rng(42);
respZero = trueResp(conc);   % no noise
fZero  = hillFit_(conc, respZero);
ciZero = confint(fZero, 0.95);
logInfo("Zero noise CI for EC50: [%.4f, %.4f]", ciZero(1,3), ciZero(2,3));
logInfo("(Note: CI may be degenerate when residuals=0 -- tool-dependent)");

%[text] Comparison of predint: observation (Observed PI) vs functional (Mean Curve PI)
fMain = hillFit_(conc, resp);
concPred = linspace(1e-3, 100, 50)';
predObs  = predint(fMain, concPred, 0.95, "observation");
predFun  = predint(fMain, concPred, 0.95, "functional");
logInfo("Prediction Interval Width at EC50 (Observation vs Functional):");
[~, idxEC50] = min(abs(concPred - TRUE_EC50));
logInfo("  Observation: %.3f  |  Functional: %.3f", ...
    predObs(idxEC50,2) - predObs(idxEC50,1), ...
    predFun(idxEC50,2) - predFun(idxEC50,1));

%[text] Answer: The more concentration points, the narrower the EC50 CI (as the data more strongly constrains the sigmoid shape).
%[text]    Doubling from 8 to 16 points approximately halves the CI width.
%[text]    With zero noise, residuals become zero, and the algorithm cannot estimate noise variance,
%[text]    resulting in a degenerate CI.
%[text]    "Functional" predint is narrower (provides CI for the true mean curve).
%[text]    "Observation" predint is wider as it adds measurement noise variance,
%[text]    accurately representing where a new single measurement might fall.
%%
%[text] ## Let's Try 5: Batch EC50 Profiling; Worst Estimate; Partial Agonist Boundary

logInfo("=== Let's Try 5: Batch EC50 Profiling ===");

%[text] Four virtual compounds (Note: Data unique to this answer file. EC50 values differ from compounds A-D in the materials)
cpdNames    = ["Compound A", "Compound B", "Compound C", "Compound D"];
cpdEC50True = [0.5,  2.0, 0.1,  8.0];
cpdEmax     = [100,  100, 100,  60];   % Compound D is a partial agonist
cpdEmin     = [5,    5,   5,    5];
nCpd        = numel(cpdNames);
CPD_EMIN_CONST = 5;

rng(42);
ec50_fitted = nan(1, nCpd);
for c = 1:nCpd
    cpdConc = logspace(-3, 2, N_CONC)';
    cpdResp = CPD_EMIN_CONST + ...
              (cpdEmax(c) - CPD_EMIN_CONST) ./ ...
              (1 + (cpdEC50True(c) ./ cpdConc).^TRUE_N) + ...
              NOISE_STD * randn(N_CONC, 1);
    % For partial agonist, allow higher Upper bound on EC50 (less constrained)
    upperEC50 = ternary_(cpdEmax(c) < 100, 200, 100);
    try
        fc = fit(cpdConc, cpdResp, hillModel, ...
            fitoptions(hillModel, ...
                'Lower',      [max(cpdEmax(c)*0.5, 30),  0.1,  0.001, -10], ...
                'Upper',      [min(cpdEmax(c)*1.5, 150), 10,   upperEC50, 50], ...
                'StartPoint', [cpdEmax(c)*0.9, 1.0, 1.0, CPD_EMIN_CONST], ...
                'TolFun',     1e-8, 'MaxIter', 2000));
        ec50_fitted(c) = fc.EC50;
        logInfo("%-12s  true EC50=%.2f  fitted EC50=%.3f  error=%+.3f", ...
            cpdNames(c), cpdEC50True(c), fc.EC50, fc.EC50 - cpdEC50True(c));
    catch ME
        logWarn("%-12s  fit failed: %s", cpdNames(c), ME.message);
    end
end

%[text] Identify the EC50 estimate furthest from the true value
validFit = ~isnan(ec50_fitted);
[~, worstRel] = max(abs(ec50_fitted(validFit) - cpdEC50True(validFit)));
idxValid = find(validFit);
logInfo("EC50 furthest from true value: %s", cpdNames(idxValid(worstRel)));

%[text] Calculate pEC50 in bulk (EC50: uM → M conversion then pEC50 = -log10)
pec50 = -log10(ec50_fitted * 1e-6);
logInfo("pEC50 Ranking (higher is more potent):");
[~, potOrd] = sort(pec50, "descend", "MissingPlacement", "last");
for k = 1:nCpd
    c = potOrd(k);
    if ~isnan(pec50(c))
        logInfo("  %d. %-12s  pEC50=%.2f", k, cpdNames(c), pec50(c));
    end
end

%[text] Answer: Compound D (partial agonist, Emax=60%) may have larger EC50 fit errors.
%[text]    The upper plateau of the sigmoid is unclear, resulting in fewer data points constraining the upper asymptote.
%[text]    Partial agonists require a wider Upper boundary for Emax
%[text]    (to prevent the optimizer from clipping the fit early).
%[text]    pEC50 = -log10(EC50 in mol/L): Higher pEC50 means lower EC50 and higher potency.
%[text]    Conversion example: EC50=1.0 uM = 1e-6 M -> pEC50=6.0.
%[text]    Compound C (EC50=0.1 uM = 1e-7 M) -> pEC50=7.0 (most potent in this set).
logInfo("A06 Answer Completed.");
%%
%[text] ## Let's Try 6: Visualization and Unit Conversion (Conceptual Questions)
%[text]
%[text] Let's Try 6 involves only conceptual and illustrative questions (no new code needed).
%[text]
%[text] - **Q1**: To convert fitted EC50 (uM) to IC50 (nM), multiplying by `* 1000` is sufficient.
%[text]   The CI in Section 4 can also be converted with the same factor (only the scale changes, the ratio remains constant).
%[text]
%[text] - **Q2**: Refer to the compound curve plots and pEC50 bar chart in Section 6.
%[text]   Compounds with larger EC50 CI error bars have lower ranking reliability.
%[text]
%[text] - **Q3**: Cooperative binding is quantified by a Hill slope n > 1.
%[text]   Allosteric modulation is a mechanism where ligand binding at a site distant from the binding site changes the shape of the target protein.
%[text]   Hemoglobin's O2 binding (n≈2.8) is known as a representative example.

%[appendix]{"version":"1.0"}
%[metadata:view]%---

%   data: {"layout":"inline","rightPanelPercent":40}
%---
