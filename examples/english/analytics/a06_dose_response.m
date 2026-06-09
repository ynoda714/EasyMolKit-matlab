%[text] # A06: Dose-Response Curve Fitting
%[text] EasyMolKit Analytics — Layer 3
%[text] 
%[text] "How does the effect change if the concentration of this drug is increased tenfold?" — This is a common question in drug discovery. The answer is quantified by **dose-response curves** and **EC50** (half-maximal effective concentration).
%[text] By using the Hill equation (4PL model), EC50 can be statistically estimated from experimental data with noise.
%[text] In this script, you will generate assay data for 4 synthetic compounds, fit them, and use `fit()` to determine EC50 with a 95% confidence interval.
%[text] 
%[text] **Story**
%[text] A pharmacologist is testing the enzyme inhibitory activity of 4 drug candidates. Each compound is applied at various concentrations, and the inhibition rate of the target enzyme is measured. The goal is to determine the EC50 (half-maximal effective concentration) for each compound, which is one of the most important efficacy indicators in drug discovery.
%[text] The Hill equation (also known as the 4-parameter logistic, 4PL, or sigmoidal dose-response model) describes changes in biological response with drug concentration.
%[text] It captures three main characteristics:
%[text] - Baseline response (no drug)
%[text] - Maximum response (saturating drug amount)
%[text] - Steepness of transition (Hill slope / cooperativity) \
%[text] In this exercise, you will:
%[text] 1. Simulate realistic dose-response datasets for 4 compounds with different EC50 values (spanning 1000-fold) and add measurement noise to mimic actual assay variability.
%[text] 2. Fit a single compound using a custom `fittype` with the 4PL model.
%[text] 3. Extract EC50 and Hill slope with 95% confidence intervals.
%[text] 4. Perform batch fitting across all compounds and summarize the results in a table.
%[text] 5. Visualize dose-response curves on a log concentration axis and rank compounds by potency. \
%[text] 
%[text] **Learning Objectives**
%[text] - Understand the Hill equation and its pharmacological interpretation
%[text] - Use MATLAB's `fit()` with custom nonlinear models (`fittype`)
%[text] - Set physically meaningful lower/upper bounds in curve fitting
%[text] - Extract and interpret EC50 confidence intervals (`confint`)
%[text] - Read and create standard dose-response plots (log scale x-axis)
%[text] - Rank compounds by potency and discuss uncertainties \
%[text] 
%[text] **Note**: This exercise does not require Python / RDKit.
%[text] It is purely MATLAB + Curve Fitting Toolbox.
%[text] 
%[text] **Prerequisites**
%[text] - Basic statistical knowledge (fitting, residuals)
%[text] - Curve Fitting Toolbox (`fit`, `fittype`, `fitoptions`, `confint`, `predint`)
%[text] - No internet connection required \
%[text] 
%[text] Estimated time: 30–45 minutes
%[text] 
%[text] **References**
%[text] - Hill AV (1910) The possible effects of the aggregation of the molecules of haemoglobin on its dissociation curves.*J Physiol* 40:iv-vii. (original Hill equation)
%[text] - Motulsky H & Christopoulos A (2004) Fitting Models to Biological Data Using Linear and Nonlinear Regression. Oxford University Press.(standard reference for dose-response fitting)
%[text] - Sebaugh JL (2011) Guidelines for accurate EC50/IC50 estimation.*Pharm Stat* 10:128-134. doi:10.1002/pst.426
%[text] - Ritz C, Baty F, Streibig JC, Gerhard D (2015) Dose-response analysis using R. *PLoS ONE* 10:e0146021. doi:10.1371/journal.pone.0146021 \
%[text] 
%[text] How to run: Execute sections one by one with Ctrl+Enter
%%
%[text] ## Section 0: Setup
%[text] **Note**: In this section, `emk.setup.initPython()` is not used. This exercise is completed using only MATLAB, and RDKit is not used.
% Resolve project root (works for Desktop, MCP, and MATLAB Online)
sDir = fileparts(mfilename('fullpath'));
if strlength(sDir) > 0
    addpath(genpath(fullfile(sDir, '..', '..', '..', 'src')));
elseif ~isempty(which("logInfo"))
    addpath(genpath(fileparts(fileparts(which("logInfo")))));
end
projectRoot = resolveProjectRoot();
addpath(genpath(fullfile(projectRoot, 'src')));
logSection("A06", "Section 0: Setup", "Analytics L3");
logInfo("A06: Setup complete (Curve Fitting Toolbox required, Python not needed)");
%%
%[text] ## Section 1: Hill Equation and EC50
%[text] Setup is complete. First, let's understand the meaning of the Hill equation and EC50 (half-maximal effective concentration), which are central to this exercise. Observe the differences in shape on a graph to develop intuition.
%[text] ### Concept: 4-Parameter Logistic (4PL) Hill Equation
%[text] The Hill equation describes a sigmoidal dose-response curve.
%[text]{"align":"center"}   y = Emin + (Emax - Emin) / (1 + (EC50 / x)^n)
%[text] 
%[text] ### Parameters:
%[text] - Emax  -- Maximum response (upper plateau, e.g., 100% inhibition)
%[text] - Emin  -- Minimum response (lower plateau, e.g., 0% at zero dose)
%[text] - EC50  -- Half-maximal effective concentration \[same units as x\]
%[text] - x = EC50: y = (Emax + Emin) / 2 (midpoint)  n     -- Hill slope (cooperativity coefficient)
%[text] -  n = 1: Simple bimolecular binding (no cooperativity)
%[text] - n \> 1: Positive cooperativity (steep transition, e.g., ion channels)
%[text] -  n \< 1: Negative cooperativity / mixed mechanism (gradual transition) \
%[text] ### Interpretation of EC50:
%[text] Lower EC50 = higher potency (less drug needed for half-maximal effect). EC50 is estimated from noisy data, so always report with a 95% confidence interval. A compound with EC50 = 1 uM ± 0.5 uM is semantically different from EC50 = 1 uM ± 10 uM.
%[text] ### IC50 vs EC50:
%[text] -   IC50: Concentration for 50% inhibition (antagonists, enzyme inhibitors)
%[text] -   EC50: Concentration for 50% maximal effect (agonists, activators)
%[text] -   Mathematically identical, but labels depend on assay context.
%[text] -   In this exercise, using enzyme inhibition, EC50 = IC50. \
%[text] 
%[text] Concentration axis: Always plot on a log10 scale. The sigmoidal shape (S-curve) only clearly appears on a logarithmic scale. On a linear scale, the curve appears as a sharp step at a single point.
%[text] To develop intuition, compare the shape of the Hill equation with different Hill slopes (n).
xDemo  = logspace(-2, 2, 200);   % concentration / EC50 (dimensionless)
nVals  = [0.5, 1.0, 2.0, 4.0];
colors = [0.8 0.2 0.2; 0.2 0.6 0.2; 0.2 0.4 0.9; 0.8 0.5 0.0];

figure("Name", "A06 Hill Equation Shapes");
for k = 1:numel(nVals)
    y_demo = 100 ./ (1 + (1 ./ xDemo).^nVals(k));
    semilogx(xDemo, y_demo, "-", LineWidth=2.0, Color=colors(k,:), ...
        DisplayName=sprintf("n = %.1f", nVals(k)));
    hold on;
end
xline(1.0, "k--", "EC50", LabelHorizontalAlignment="right", LineWidth=1.5, ...
    HandleVisibility="off");
yline(50,  "k:",  "50%",  LabelVerticalAlignment="bottom",  LineWidth=1.0, ...
    HandleVisibility="off");
xlabel("Concentration / EC50 (log scale)"); ylabel("Response (%)");
title("Hill Equation Shapes with Different Hill Slopes (n)");
legend(Location="northwest"); grid on; ylim([0 110]);

logInfo("A06: Plotted Hill equation demo. Observe the S-curve on a log scale.");
%[text] **💡 Observation Point 1**
%[text] Note that with n = 4, the transition is very steep (close to all-or-nothing).
%[text] Consider biological mechanisms that result in Hill slopes \> 3.
%[text] (Hint: Hemoglobin O2 binding is n ~ 2.8; voltage-dependent Na+ channels can have effective n \> 4)
%[text] For n = 1, observe the response at x = 10*EC50 and x = 0.1*EC50.
%[text] Verify analytically: y = 100 / (1 + (EC50/x)^1)
%[text] Consider why EC50 is called "half-maximal" and algebraically show that
%[text] for any value of n, y = (Emin + Emax)/2 when x = EC50.
% ... (Try writing code here)
%%
%[text] ## Section 2: Generating Synthetic Dose-Response Data
%[text] 
%[text] We have deepened our understanding of the Hill equation. Next, we will create synthetic data that mimics actual assays.
%[text] By adding noise, we replicate the variability found in experiments.
%[text] 
%[text] ### Concept: Simulation of Realistic Assay Measurements
logSection("A06", "Section 2: Generating Synthetic Dose-Response Data", "Analytics L3");
%[text] In actual dose-response assays, bioactivity is measured at 8-12 concentrations per compound, typically sampled at logarithmic intervals (semi-log or full-log steps) to cover 3-4 orders of magnitude around the EC50.
%[text] 
%[text] Measurement noise in cell-based assays is usually expressed as additive Gaussian noise in response units (~3-8% standard deviation).
%[text] In enzyme assays, a 5% standard deviation is common.
%[text] 
%[text] Four synthetic compounds represent a potency range over 3 log units.
%[text] 
%[text] Compound A -- High potency    (EC50 = 0.5  uM, n = 1.2)
%[text] Compound B -- Moderate potency (EC50 = 5.0  uM, n = 1.0)
%[text] Compound C -- Low potency    (EC50 = 50.0 uM, n = 2.0)
%[text] Compound D -- Very low (EC50 = 500  uM, n = 0.8, partial agonist)
%[text] 
%[text] Partial agonist (Compound D): Emax = 80% means it cannot achieve full inhibition even at saturating concentrations. This is common for compounds that only partially block the active site.
rng(2026);  % fix seed for reproducibility
%[text] Set the true parameters for each compound.
compNames = ["Compound A", "Compound B", "Compound C", "Compound D"];
EC50_true = [0.5, 5.0, 50.0, 500.0];   % uM
n_true    = [1.2, 1.0,  2.0,   0.8];   % Hill slope
Emax_true = [100, 100, 100,   80.0];   % % response (top)
Emin_true = [  0,   0,   0,    5.0];   % % response (bottom)
noiseStd  = 5.0;                        % assay noise (% response units)

N_CONC    = 10;   % concentrations per compound
%[text] Store all data in a structure for batch fitting.
data = struct();
for c = 1:numel(compNames)
    % Log-spaced concentrations: 2 log units below to 2 log units above EC50
    logMin = log10(EC50_true(c)) - 2;
    logMax = log10(EC50_true(c)) + 2;
    conc   = logspace(logMin, logMax, N_CONC)';  % uM, column vector

    % True Hill response + noise
    yTrue  = Emin_true(c) + (Emax_true(c) - Emin_true(c)) ./ ...
             (1 + (EC50_true(c) ./ conc).^n_true(c));
    yNoisy = yTrue + noiseStd * randn(N_CONC, 1);
    yNoisy = max(0, min(100, yNoisy));  % clamp to [0, 100]%

    data(c).name     = compNames(c);
    data(c).conc     = conc;     % observed concentrations (uM)
    data(c).response = yNoisy;  % observed responses (% inhibition)
    data(c).EC50     = EC50_true(c);
    data(c).n        = n_true(c);
    data(c).Emax     = Emax_true(c);
    data(c).Emin     = Emin_true(c);
end

logInfo("Generated dose-response data for %d compounds (each %d points, noise=%.0f%%)", ...
    numel(compNames), N_CONC, noiseStd);
%[text] Let's take a quick look at the data (display all 4 compounds on a single log-scale graph).
figure("Name", "A06 Raw Data");
set(gcf, "Position", [100 100 900 520]);
tiledlayout(2, 2, "TileSpacing", "loose");
for c = 1:numel(compNames)
    nexttile;
    semilogx(data(c).conc, data(c).response, "o", MarkerSize=8, ...
        MarkerFaceColor=[0.3 0.5 0.9], MarkerEdgeColor="k");
    xlabel("Concentration (uM, log scale)"); ylabel("Inhibition Rate (%)");
    title(sprintf("%s\nEC50 = %.1f uM (True Value)", data(c).name, data(c).EC50));
    ylim([-5 110]); grid on;
end
sgtitle("A06: Raw Dose-Response Data (Synthetic, Noise = 5%)");
%[text] **💡 Observation Point 2**
%[text] Try noiseStd = 15.0 (high-noise assay) and replot to observe changes in scatter.
%[text] Determine which compound's EC50 is most difficult to estimate accurately and why.
%[text] In actual assays, measurements are often taken three times at each concentration.
%[text] Modify the code to have three replicate measurements per concentration (N\_CONC x 3).
%[text] Consider how to calculate and display the mean ± SEM for each concentration.
% ... (Try writing the code here)
%%
%[text] ## Section 3: Fitting a Single Compound with the 4PL Hill Equation
%[text] 
%[text] The data is ready. Let's learn the fitting process using compound A.
%[text] Check how to use `fit()`, set initial values (StartPoint), and constraints.
%[text] 
%[text] ### Concept: Nonlinear Least Squares Fitting with `fit()`
%[text] MATLAB's `fit()` (Curve Fitting Toolbox) minimizes the sum of squared residuals between the model and observations.
%[text] Minimize for {Emax, n, EC50, Emin}: sum\_i (y\_i - y\_model(x\_i))^2
%[text] 
%[text] The Hill equation is nonlinear with respect to its parameters, so an iterative algorithm (default is Levenberg-Marquardt) is used for minimization.
%[text] This algorithm is sensitive to initial values.
%[text] 
%[text] StartPoint: Initial estimates for each parameter
%[text] Lower/Upper: Box constraints to enforce physical limits
%[text] 
%[text] Physical boundary conditions for inhibition assays:
%[text] Emax  \[50, 120\]  -- Upper plateau (close to ~100%; slight excess allowed)
%[text] n     \[0.1, 5.0\] -- Hill slope (rarely > 4 in practice)
%[text] EC50  \[1e-4, 1e6\] -- Assay concentration range (uM)
%[text] Emin  \[-10, 30\]  -- Lower plateau (close to ~0%; baseline noise allowed)
%[text] 
%[text] Fit compound A (index 1) as a working sample.
c = 1;
xFit = data(c).conc;
yFit = data(c).response;

%[text] Define the 4-parameter Hill model as a fittype.
%[text] Order of coefficients: {Emax, n, EC50, Emin}
%[text] fittype requires char literals (not strings) for expressions and parameter names.
hillModel = fittype( ...
    'Emin + (Emax - Emin) ./ (1 + (EC50 ./ x).^n)', ...
    'independent', 'x', ...
    'coefficients', {'Emax', 'n', 'EC50', 'Emin'});

fitOpts = fitoptions(hillModel);
fitOpts.Lower      = [50,  0.1, 1e-4, -10];    % [Emax, n, EC50, Emin]
fitOpts.Upper      = [120, 5.0, 1e6,   30];
fitOpts.StartPoint = [100, 1.0, median(xFit), 0];

[fitResult, gof] = fit(xFit, yFit, hillModel, fitOpts);

%[text] Extract the fitted parameters.
ec50Fit  = fitResult.EC50;
nFit     = fitResult.n;
EmaxFit  = fitResult.Emax;
EminFit  = fitResult.Emin;

logInfo("Compound A fit: EC50=%.3f uM (true value=%.1f)  n=%.2f (true value=%.1f)", ...
    ec50Fit, data(c).EC50, nFit, data(c).n);
logInfo("               Emax=%.1f%%  Emin=%.1f%%  R2=%.4f  RMSE=%.2f%%", ...
    EmaxFit, EminFit, gof.rsquare, gof.rmse);

%[text] Overlay the fitted curve on the data and plot.
xPlot  = logspace(log10(min(xFit)) - 0.5, log10(max(xFit)) + 0.5, 500);
yPlot  = EminFit + (EmaxFit - EminFit) ./ (1 + (ec50Fit ./ xPlot).^nFit);

figure("Name", "A06 Single Fit -- Compound A");
semilogx(xFit,  yFit,  "o", MarkerSize=10, MarkerFaceColor=[0.2 0.5 0.9], ...
    MarkerEdgeColor="k", DisplayName="Observed"); hold on;
semilogx(xPlot, yPlot, "b-", LineWidth=2.5, DisplayName="4PL Hill fit");
xline(ec50Fit, "--r", sprintf("EC50 = %.2f uM", ec50Fit), ...
    LineWidth=1.5, LabelHorizontalAlignment="right", HandleVisibility="off");
yline(50, ":k", LineWidth=1.0, HandleVisibility="off");
xlabel("Concentration (uM, log scale)"); ylabel("Inhibition Rate (%)");
title(sprintf("A06: %s -- 4PL Hill Fit  (R^2=%.3f)", data(c).name, gof.rsquare));
legend(Location="northwest"); ylim([-5 110]); grid on;

%[text] **💡 Observation Point 3**
%[text] Check the result of fitting with StartPoint = \[100, 3.0, median(xFit), 0\].
%[text] Did the algorithm converge to the same EC50? If not, consider why.
%[text] (Hint: The Levenberg-Marquardt algorithm can get trapped in local minima)
%[text] Check the result of refitting without Lower/Upper constraints. How does the result change without constraints?
%[text] Is it possible for Emax to become unrealistically large (e.g., 500%)?
%[text] Consider why constraints are important in biological dose-response fitting.
%[text] The Hill equation was originally derived for oxygen binding to hemoglobin.
%[text] Investigate its derivation and consider the assumption that the Hill slope equals the number of ligand binding sites.
% ... (Try writing code here)
%%
%[text] ## Section 4: Confidence Interval of EC50
%[text] 
%[text] The fitting for a single compound is complete. It is essential to report not only the fit values but also their uncertainty (95% confidence interval).
%[text] Visualize the confidence interval of EC50 and evaluate the reliability of the estimation.
%[text] 
%[text] ### Concept: Quantification of Uncertainty in Nonlinear Fitting
logSection("A06", "Section 4: Confidence Interval of EC50", "Analytics L3");
%[text] `confint()` calculates the 95% confidence interval for each fit parameter based on the curvature of the sum of squares surface at the minimum.
%[text] A narrow CI indicates well-determined parameters, while a wide CI means the data cannot precisely fix the parameters.
%[text] 
%[text] The width of the CI depends on the following factors:
%[text] - Number of data points (more points result in a narrower CI)
%[text] - Noise level (lower noise results in a narrower CI)
%[text] - Data coverage (more measurement points near EC50 result in a narrower CI for EC50)
%[text] - Model identifiability (whether all parameters can be independently estimated) \
%[text] 
%[text] Guideline: If the 95% CI of log10(EC50) exceeds 1 unit (10-fold uncertainty range), the EC50 estimate is unreliable.
%[text] Consider repeating the assay or extending the concentration range.
ci = confint(fitResult, 0.95);   % 2 x 4 matrix: [lower; upper] for [Emax, n, EC50, Emin]
ec50Lo = ci(1, 3);
ec50Hi = ci(2, 3);
nLo    = ci(1, 2);
nHi    = ci(2, 2);

logInfo("Compound A -- 95%% CI:");
logInfo("  EC50 = %.3f uM  [%.3f, %.3f]  (Uncertainty %.1f-fold)", ...
    ec50Fit, ec50Lo, ec50Hi, ec50Hi / max(ec50Lo, 1e-9));
logInfo("  n    = %.2f     [%.2f, %.2f]", nFit, nLo, nHi);

%[text] Plot with 95% prediction interval band.
yPI = predint(fitResult, xPlot, 0.95, "functional", "off");  % 95% PI of the curve

figure("Name", "A06 Confidence Interval Band -- Compound A");
fill([xPlot, fliplr(xPlot)], [yPI(:,1)', fliplr(yPI(:,2)')], ...
    [0.7 0.8 1.0], FaceAlpha=0.4, EdgeColor="none", DisplayName="95% PI band");
hold on;
semilogx(xFit,  yFit,  "o", MarkerSize=10, MarkerFaceColor=[0.2 0.5 0.9], ...
    MarkerEdgeColor="k", DisplayName="Observed values");
semilogx(xPlot, yPlot, "b-", LineWidth=2.5, DisplayName="4PL Fit");
xline(ec50Fit, "--r", sprintf("EC50 = %.2f [%.2f, %.2f] uM", ec50Fit, ec50Lo, ec50Hi), ...
    LineWidth=1.5, LabelHorizontalAlignment="right", ...
    LabelVerticalAlignment="bottom", HandleVisibility="off");
xlabel("Concentration (uM, log scale)"); ylabel("Inhibition rate (%)");
title(sprintf("A06: %s -- 95%% Confidence Interval Band", data(c).name));
legend(Location="northwest"); ylim([-10 115]); grid on;
set(gca, "XScale", "log");

logInfo("Displaying 95%% prediction interval band.");

%[text] **💡 Observation Point 4**
%[text] When the simulation is rerun with N\_CONC = 5 (half the number of data points),
%[text] observe how the width of the 95% CI of EC50 changes.
%[text] Determine the typically recommended minimum number of concentrations.
%[text] (Reference: Sebaugh 2011 recommends \>= 10 for reliable 4PL fitting)
%[text] When noiseStd = 0 is set and refitted, check if the CI shrinks to nearly zero.
%[text] Consider why the CI does not become exactly zero in practice.
%[text] "functional" in predint provides the CI of the mean curve.
%[text] Change to "observation" to obtain the CI for a single new observation and note the difference.
%[text] Explain why the observation PI is always wider than the functional PI.
% ... (Try writing code here)
%%
%[text] ## Section 5: Batch Fitting for All Compounds
%[text] 
%[text] You have learned the fitting procedure for a single compound. We will now automatically apply this to all 4 compounds and compare their potencies. Let's check the results in a summary table.
%[text] 
%[text] ### Concept: Systematic Potency Profiling
%[text] In actual drug discovery campaigns, the same fitting procedure is applied to dozens to hundreds of compounds.
%[text] Batch fitting automatically extracts the EC50 values and their confidence intervals for each compound.
%[text] 
%[text] The results are reported in a summary table ranking compounds by potency (lowest EC50 = highest potency).
%[text] The width of the confidence interval (CI) indicates reliability.
%[text] 
%[text] CI width = EC50\_upper / EC50\_lower
%[text] < 2 times: Excellent precision
%[text] 2-5 times: Acceptable range (typical for noisy cell assays)
%\[text]   > 10 times: Unreliable. Consider extending the concentration range or repeating the assay.
fitSummary = struct();

for c = 1:numel(compNames)
    xc = data(c).conc;
    yc = data(c).response;

    % Adaptive StartPoint: Start EC50 from the geometric mean of concentrations
    startEC50 = exp(mean(log(xc)));

    fopts = fitoptions(hillModel);
    fopts.Lower      = [50,  0.1, 1e-4, -10];
    fopts.Upper      = [120, 5.0, 1e6,   30];
    fopts.StartPoint = [100, 1.0, startEC50, 0];

    try
        [fr, gf] = fit(xc, yc, hillModel, fopts);
        ci_c     = confint(fr, 0.95);

        fitSummary(c).Name      = data(c).name;
        fitSummary(c).EC50      = fr.EC50;
        fitSummary(c).EC50_lo   = ci_c(1, 3);
        fitSummary(c).EC50_hi   = ci_c(2, 3);
        fitSummary(c).HillSlope = fr.n;
        fitSummary(c).Emax      = fr.Emax;
        fitSummary(c).R2        = gf.rsquare;
        fitSummary(c).RMSE      = gf.rmse;
        fitSummary(c).fitObj    = fr;

        logInfo("%-12s  EC50=%7.2f uM [%6.2f, %7.1f]  n=%.2f  R2=%.3f", ...
            data(c).name, fr.EC50, ci_c(1,3), ci_c(2,3), fr.n, gf.rsquare);
    catch ME
        logWarn("Fit failed for %s: %s", data(c).name, ME.message);
        fitSummary(c).Name = data(c).name;
        fitSummary(c).EC50 = NaN;
    end
end

%[text] Construct the summary table.
ec50Vec = [fitSummary.EC50];
loVec   = [fitSummary.EC50_lo];
hiVec   = [fitSummary.EC50_hi];
nVec    = [fitSummary.HillSlope];
r2Vec   = [fitSummary.R2];

summTbl = table(compNames(:), EC50_true(:), ec50Vec(:), loVec(:), hiVec(:), ...
    n_true(:), nVec(:), r2Vec(:), ...
    VariableNames=["Compound","EC50_true_uM","EC50_fit_uM", ...
                   "CI_lower","CI_upper","n_true","n_fit","R2"]);
disp(summTbl);

%[text] **💡 Observation Point 5**
%[text] Check which compound has the fitted EC50 furthest from the true value.
%[text] Verify if the true value is always within the 95% confidence interval.
%[text] If not, consider whether this indicates a bug. (Hint: 95% CI does not guarantee inclusion in all experiments -- only 95% of the time on average)
%[text] Compound D has Emax = 80% (partial inhibition). Consider how to adjust the lower bound to properly allow fitting of partial agonists.
%[text] (The current constraint forces Emax \>= 50%, which may bias EC50 estimation)
%[text] Extend the batch fitting loop to also calculate pEC50 = -log10(EC50\_uM) and confirm that pEC50 is commonly used in medicinal chemistry to convert potency to a linear scale.
%[text] (pEC50 = 6 means EC50 = 1 uM).
% ... (Try writing the code here)


%%
%[text] ## Section 6: Potency Comparison and Multi-Panel Visualization
%[text] 
%[text] We have obtained the EC50 (half-maximal effective concentration) for all compounds. Finally, we will overlay the dose-response curves of all compounds for visual comparison and rank them by potency.
%[text] We will also create a ranking bar chart on the pEC50 scale.
%[text] 
%[text] ### Concept: Visualization and Ranking of Compound Potency
%[text] Dose-response curves are plotted on a logarithmic concentration axis. On this scale, a well-fitted 4PL Hill curve appears as a symmetric S-shape (sigmoid), with the EC50 corresponding to the inflection point.
%[text] 
%[text] **Ranking by Potency:** EC50 is reported in ascending order. Since the EC50 distribution is log-normal, it is not reported on a linear scale. Example report: `Compound A: EC50 = 0.5 uM (95% CI: 0.3 — 0.9 uM)`
%[text] 
logSection("A06", "Section 6: Potency Comparison and Multi-Panel Visualization", "Analytics L3");
%[text] Overlay and plot the dose-response curves of all compounds on a single graph.
cmap = lines(numel(compNames));

figure("Name", "A06 All Dose-Response Curves");
xPlotGlobal = logspace(-3, 4, 500);  % cover full concentration range (uM)

for c = 1:numel(compNames)
    if isnan(fitSummary(c).EC50), continue; end
    fr = fitSummary(c).fitObj;
    yc = fr.Emin + (fr.Emax - fr.Emin) ./ ...
         (1 + (fr.EC50 ./ xPlotGlobal).^fr.n);
    semilogx(xPlotGlobal, yc, "-", LineWidth=2.5, Color=cmap(c, :), ...
        DisplayName=sprintf("%s  EC50=%.1f uM", compNames(c), fr.EC50));
    hold on;
    % Mark EC50 on the curve (midpoint = (Emin+Emax)/2 is the exact y-coordinate)
    semilogx(fr.EC50, (fr.Emax + fr.Emin) / 2, "v", MarkerSize=12, MarkerFaceColor=cmap(c, :), ...
        MarkerEdgeColor="k", HandleVisibility="off");
    % Raw data points (MarkerFaceAlpha is R2020b+; use set() for compatibility)
    h_pts = semilogx(data(c).conc, data(c).response, "o", MarkerSize=6, ...
        MarkerFaceColor=cmap(c, :), MarkerEdgeColor="k", ...
        HandleVisibility="off");
    try; set(h_pts, "MarkerFaceAlpha", 0.6); catch; end
end
yline(50, "--k", "50% (EC50 threshold)", LabelHorizontalAlignment="left", LineWidth=1.0, HandleVisibility="off");
xlabel("Concentration (uM, log scale)"); ylabel("Inhibition Rate (%)");
title("A06: All Compounds -- Fitted Dose-Response Curves");
legend(Location="east"); ylim([-5 110]); grid on;

%[text] Create a potency ranking bar graph (pEC50 scale).
pEC50_true = -log10(EC50_true * 1e-6);  % convert uM to M first
pEC50_fit  = -log10(ec50Vec  * 1e-6);

[~, sortIdx] = sort(pEC50_fit, "descend");  % rank by potency (highest pEC50 first)

figure("Name", "A06 Potency Ranking");
colors_bar = cmap(sortIdx, :);
bh = bar(pEC50_fit(sortIdx), FaceColor="flat", DisplayName="Fitted pEC50");
bh.CData = colors_bar;
hold on;
%[text] Add CI error bars (asymmetric on log scale → symmetric in log space).
pLo = -log10(hiVec(sortIdx) * 1e-6);  % note: EC50_hi -> lower pEC50
pHi = -log10(loVec(sortIdx) * 1e-6);  % note: EC50_lo -> higher pEC50
errorbar(1:numel(compNames), pEC50_fit(sortIdx), ...
    pEC50_fit(sortIdx) - pLo, pHi - pEC50_fit(sortIdx), ...
    ".k", LineWidth=1.5, DisplayName="95% CI");
plot(1:numel(compNames), pEC50_true(sortIdx), "r^", MarkerSize=10, ...
    MarkerFaceColor="r", DisplayName="True pEC50");
set(gca, "XTickLabel", compNames(sortIdx), "XTickLabelRotation", 20);
ylabel("pEC50 = -log_{10}[EC50 (M)]"); title("A06: Compound Potency Ranking");
legend(Location="northeast"); grid on;
logInfo("Potency ranking complete. Higher pEC50 indicates greater potency.");

%[text] **💡 Observation Point 6**
%[text] Let's convert the fitted EC50 values (uM) to IC50 (nM) and verify.
%[text] (1 uM = 1000 nM) Identify which compounds are referred to as having "nanomolar potency."
%[text] Drug candidates typically advance if IC50 < 1 uM (pEC50 > 6).
%[text] Compound D has a small Emax (partial inhibition). Consider how this affects pEC50 estimation.
%[text] Consider if it is fair to rank compound D by pEC50 alongside full inhibitors.
%[text] Consider what additional metrics should be reported.
%[text] In actual campaigns, compounds are often retested in repeated experiments.
%[text] If Compound A gives EC50 = 0.3 uM in one run and EC50 = 0.9 uM in another, check if these results are consistent.
%[text] Compare with the 95% CI calculated in Section 4.
%[text] The Hill equation assumes a single binding site when n = 1.
%[text] Consider what structural information can explain why Compound C has n ~ 2.
%[text] Investigate "cooperative binding" and "allosteric regulation."
% ... (Try writing code here)

%[appendix]{"version":"1.0"}
%[metadata:view]%---

%   data: {"layout":"inline","rightPanelPercent":40}
%---
