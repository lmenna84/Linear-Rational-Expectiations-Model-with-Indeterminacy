%% EXAMPLES FOR lubikschorfheide.m
% This script demonstrates how to use the lubikschorfheide function
% to solve DSGE models using the Lubik-Schorfheide (2003) method.
% Each example uses a basic New Keynesian model with different options.

clear
addpath('../functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 1: Basic 3-Equation NK Model (Linear, Determinacy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define symbolic variables using the eval-loop pattern
% (must match the pattern used internally by lubikschorfheide.m)
% _f denotes E_t[x_{t+1}] (forward), _l denotes x_{t-1} (lag)
vars = {'y' 'p' 'ii' 'rn'};

var   = zeros(0,0);
var_l = zeros(0,0);
var_f = zeros(0,0);
for jj = 1:size(vars,2)
    eval(['syms ' vars{jj}]);
    eval(['var   = [var '   vars{jj} '];']);
    eval(['syms ' vars{jj} '_l']);
    eval(['var_l = [var_l ' vars{jj} '_l];']);
    eval(['syms ' vars{jj} '_f']);
    eval(['var_f = [var_f ' vars{jj} '_f];']);
end

% Define symbolic shocks
syms eps_r eps_y;  % monetary policy shock, natural rate shock

% Model parameters (canonical NK calibration)
bet   = 0.99;    % discount factor
sig   = 1;       % inverse EIS
kap   = 0.1;     % slope of Phillips curve
phi_pi = 1.5;    % Taylor rule response to inflation
phi_y  = 0.5;    % Taylor rule response to output gap
rho_rn = 0.8;    % persistence of natural rate

% Model equations (already linear, in deviation from steady state):
% 1. IS curve:      y_t = E_t[y_{t+1}] - sig*(ii_t - E_t[p_{t+1}] - rn_t)
% 2. Phillips curve: p_t = bet*E_t[p_{t+1}] + kap*y_t
% 3. Taylor rule:   ii_t = phi_pi*p_t + phi_y*y_t + eps_r_t
% 4. Natural rate:   rn_t = rho_rn*rn_{t-1} + eps_y_t

EQ1 = y - y_f + sig*(ii - p_f - rn);
EQ2 = p - bet*p_f - kap*y;
EQ3 = ii - phi_pi*p - phi_y*y - eps_r;
EQ4 = rn - rho_rn*rn_l - eps_y;

EQ = [EQ1; EQ2; EQ3; EQ4];
% var is already built by the eval-loop above
shocks = [eps_r eps_y];
std_shocks = [0.0025 0.005];  % 25bp monetary shock, 50bp demand shock

% Solve (linear model, no steady state needed)
Q1 = lubikschorfheide(EQ, var, shocks, std_shocks, [], 1);
close all

disp('Example 1: Basic NK model solved successfully');
disp(['Degree of indeterminacy: ' num2str(Q1.indetdegree)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 2: Nonlinear NK Model with Automatic Linearization (linear=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A simple New Keynesian model written in levels with nonlinear terms.
% lubikschorfheide automatically linearizes around the provided steady state.
%
% When using linear=0, two requirements must hold:
%   1. Equations are in LEVELS (not deviations from ss). The function
%      internally redefines variables as deviations and Taylor-expands.
%   2. The steady state ss must satisfy EQ(ss, ss, ss, 0) = 0.
%
% Model equations (nonlinear):
%   1. Euler:    y^(-sig)=bet*y_f^(-sig)(i/p_f)
%   2. Phillips: y^(eta-sig)=1+K/(mu-1)*(p-1)*p-bet*K/(mu-1)*(y_f/y)^(-sig)*(p_f-1)*p_f*y_f/y
%   3. Taylor:   ii = ii_bar*p^phi_pi*y^phi_y*exp(eps_r)  
%
% The log-linearized model matches Example 1 (without the natural rate process),
% so the results should be comparable, once the irfs are reported in percentage change.
K=77;
eta=2;
% The log-linearized Rotemberg NKPC slope is (eta-sig)*(mu-1)/K.
% To match the linear model's kap, set mu so that (eta-sig)*(mu-1)/K = kap:
mu=kap*K/(eta-sig)+1;

% Steady state (analytical):
y_ss=1;
p_ss=1;
ii_ss=1/bet;
ss_nl=[y_ss,p_ss,ii_ss];

EQ1_nl = y^(-sig)-bet*y_f^(-sig)*(ii/p_f);
EQ2_nl = -y^(eta-sig)+1+K/(mu-1)*(p-1)*p-bet*K/(mu-1)*(y_f/y)^(-sig)*(p_f-1)*p_f*y_f/y;
EQ3_nl = -ii + ii_ss*p^phi_pi*y^phi_y*exp(eps_r);

EQ_nl = [EQ1_nl; EQ2_nl; EQ3_nl];
% Build var_nl using the same eval-loop pattern
var_nl = zeros(0,0);
vars_nl = {'y' 'p' 'ii'};
for jj = 1:size(vars_nl,2)
    eval(['var_nl = [var_nl ' vars_nl{jj} '];']);
end
shocks_nl = [eps_r];
std_nl = [0.0025];  % 25bp monetary shock


fprintf('\nExample 2: Steady state values:\n');
fprintf('  y_ss = %.6f, p_ss = %.6f, ii_ss = %.6f\n', y_ss, p_ss, ii_ss);


% Solve with automatic linearization (linear=0)
Q2 = lubikschorfheide(EQ_nl, var_nl, shocks_nl, std_nl, ss_nl, 0, 40, 1);
close all

disp('Example 2: Nonlinear NK model with linearization solved');
disp(['Degree of indeterminacy: ' num2str(Q2.indetdegree)]);

% --- Comparison plot: nonlinear (% deviations) vs linear (log-deviations) ---
% Example 1 variables are log-deviations from ss (standard log-linearized NK).
% Example 2 variables are absolute deviations from ss (lubikschorfheide
% linearizes around ss). To convert to percentage deviations, divide by ss.
% Since y_ss = 1 and p_ss = 1, output and inflation IRFs coincide directly.
% For the interest rate, divide by ii_ss = 1/bet.

T_irf = 40;
t_plot = 1:T_irf;

% Example 1: log-deviation IRFs to eps_r (first 3 variables: y, p, ii)
irf1_y  = Q1.irf_lubik.eps_r(1, 2:end);
irf1_p  = Q1.irf_lubik.eps_r(2, 2:end);
irf1_ii = Q1.irf_lubik.eps_r(3, 2:end);

% Example 2: absolute deviation IRFs to eps_r, converted to % deviations
irf2_y  = log((Q2.irf_lubik.eps_r(1, 2:end)+y_ss) / y_ss);
irf2_p  = log((Q2.irf_lubik.eps_r(2, 2:end) +p_ss)/ p_ss);
irf2_ii = log((Q2.irf_lubik.eps_r(3, 2:end)+ii_ss) / ii_ss);

figure('Name', 'Example 2: Nonlinear vs Linear IRFs', 'NumberTitle', 'off');

subplot(1,3,1);
plot(t_plot, irf1_y, 'b-', 'LineWidth', 2); hold on;
plot(t_plot, irf2_y, 'r--', 'LineWidth', 2);
title('Output (y)');
xlabel('Periods'); ylabel('% deviation from SS');
legend('Linear (Ex.1)', 'Nonlinear (Ex.2)', 'Location', 'best');
grid on;

subplot(1,3,2);
plot(t_plot, irf1_p, 'b-', 'LineWidth', 2); hold on;
plot(t_plot, irf2_p, 'r--', 'LineWidth', 2);
title('Inflation (p)');
xlabel('Periods'); ylabel('% deviation from SS');
legend('Linear (Ex.1)', 'Nonlinear (Ex.2)', 'Location', 'best');
grid on;

subplot(1,3,3);
plot(t_plot, irf1_ii, 'b-', 'LineWidth', 2); hold on;
plot(t_plot, irf2_ii, 'r--', 'LineWidth', 2);
title('Interest rate (ii)');
xlabel('Periods'); ylabel('% deviation from SS');
legend('Linear (Ex.1)', 'Nonlinear (Ex.2)', 'Location', 'best');
grid on;

sgtitle('Monetary shock IRFs: Linear (Ex.1) vs Nonlinear in % dev. (Ex.2)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 3: Custom IRF Length and No Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve with 60-period IRFs instead of default 40, no plots
Q3 = lubikschorfheide(EQ, var, shocks, std_shocks, [], 1, 60, 0);

disp('Example 3: Solved with 60-period IRFs, no plots');
disp(['IRF size: ' num2str(size(Q3.irf_lubik.eps_r, 2)) ' periods']);

% Manually plot selected IRFs
figure('Name', 'Example 3: Custom IRF Plot');
subplot(2,2,1);
plot(1:60, Q3.irf_lubik.eps_r(1, 2:end), 'LineWidth', 2);
title('Output response to monetary shock');
xlabel('Periods'); ylabel('Deviation from SS');
grid on;

subplot(2,2,2);
plot(1:60, Q3.irf_lubik.eps_r(2, 2:end), 'LineWidth', 2);
title('Inflation response to monetary shock');
xlabel('Periods'); ylabel('Deviation from SS');
grid on;

subplot(2,2,3);
plot(1:60, Q3.irf_lubik.eps_y(1, 2:end), 'LineWidth', 2);
title('Output response to demand shock');
xlabel('Periods'); ylabel('Deviation from SS');
grid on;

subplot(2,2,4);
plot(1:60, Q3.irf_lubik.eps_y(2, 2:end), 'LineWidth', 2);
title('Inflation response to demand shock');
xlabel('Periods'); ylabel('Deviation from SS');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 4: Indeterminacy Case (Passive Monetary Policy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Modify Taylor rule to violate Taylor principle: phi_pi < 1
phi_pi_passive = 0.8;  % Passive policy

EQ1_indet = y - y_f + sig*(ii - p_f - rn);
EQ2_indet = p - bet*p_f - kap*y;
EQ3_indet = ii - phi_pi_passive*p - phi_y*y - eps_r;  % Passive rule
EQ4_indet = rn - rho_rn*rn_l - eps_y;

EQ_indet = [EQ1_indet; EQ2_indet; EQ3_indet; EQ4_indet];

Q4 = lubikschorfheide(EQ_indet, var, shocks, std_shocks, [], 1, 40, 1);

disp('Example 4: Indeterminate model (passive monetary policy)');
disp(['Degree of indeterminacy: ' num2str(Q4.indetdegree)]);

if Q4.indetdegree > 0
    disp('Model exhibits sunspot equilibria!');
    disp('IRFs include belief shock responses:');
    disp(fieldnames(Q4.irf_lubik));
end

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 5: Instability Case with Forced Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To illustrate instability, we use a simple model with a predetermined
% variable whose dynamics are explosive. Consider:
%   (1) d_t = a * d_{t-1} + b * c_t + eps_t   (debt accumulation)
%   (2) c_t = phi * E_t[c_{t+1}]               (Euler equation)
%
% Variable d is predetermined (depends on d_{t-1}), c is forward-looking.
% With a > 1, debt dynamics are explosive. This makes the predetermined
% eigenvalue too large, producing fewer stable roots than needed.

vars_ex5 = {'d' 'c'};
var_ex5 = zeros(0,0); var_l_ex5 = zeros(0,0); var_f_ex5 = zeros(0,0);
for jj = 1:size(vars_ex5,2)
    eval(['syms ' vars_ex5{jj}]);
    eval(['var_ex5 = [var_ex5 ' vars_ex5{jj} '];']);
    eval(['syms ' vars_ex5{jj} '_l']);
    eval(['var_l_ex5 = [var_l_ex5 ' vars_ex5{jj} '_l];']);
    eval(['syms ' vars_ex5{jj} '_f']);
    eval(['var_f_ex5 = [var_f_ex5 ' vars_ex5{jj} '_f];']);
end
syms eps_d;

a_debt = 1.2;   % explosive debt dynamics (> 1)
b_debt = 0.3;
phi_c  = 0.5;

EQ_unstable = [d - a_debt*d_l - b_debt*c - eps_d;
               c - phi_c*c_f];

shocks_ex5 = [eps_d];
std_ex5    = [0.01];

% 5a: Without forcing a solution (instsol=0). lubikschorfheide returns
%     early under instability without assigning Q, so we use try-catch.
disp('Example 5a: Attempting to solve unstable model...');
try
    Q5a = lubikschorfheide(EQ_unstable, var_ex5, shocks_ex5, std_ex5, [], 1, 40, 0, 0);
catch
    disp('  As expected, no solution exists (instability).');
end

% 5b: Force a solution by treating the least explosive roots as stable
%     (instsol=1). This is NOT a valid REE but can be useful for inspection.
disp('Example 5b: Forcing solution for unstable model (instsol=1)...');
Q5b = lubikschorfheide(EQ_unstable, var_ex5, shocks_ex5, std_ex5, [], 1, 40, 1, 1);

disp('Forced solution obtained.');
disp('Warning: This is not a valid REE solution.');

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 6: Analyzing the Solution Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The solution takes the form (in the stable block W1 of the QZ decomposition):
%   W1_t = V1lubik * W1_{t-1} + V2lubik_shocks * shocks_t
%                               + V2lubik_belief * belief_t   (indeterminacy only)
% Original variables are recovered via: X_t = Z * [W1_t; W2_t]

disp(' ');
disp('=================================================================');
disp('Example 6: Understanding the solution structure');
disp('=================================================================');

% --- 6a: Determinacy (Q1 from Example 1) ---
disp(' ');
disp('--- 6a: DETERMINACY (Example 1, phi_pi = 1.5) ---');
disp(' ');
disp('Autoregressive matrix V1lubik:');
disp(Q1.V1lubik);
fprintf('  Eigenvalues of V1lubik: ');
eig_V1 = eig(Q1.V1lubik);
fprintf('%.4f  ', eig_V1);
fprintf('\n  All inside unit circle => stable dynamics.\n');

disp(' ');
disp('Shock impact matrix V2lubik_shocks (numeric, no free parameters):');
disp(double(Q1.V2lubik_shocks));

disp('Under determinacy, the solution is UNIQUE:');
disp('  - V2lubik_shocks is fully determined (no symbolic M1).');
disp('  - V2lubik_belief is empty (no belief shocks).');
disp(['  - Degree of indeterminacy: ' num2str(Q1.indetdegree)]);

% --- 6b: Indeterminacy (Q4 from Example 4) ---
disp(' ');
disp('--- 6b: INDETERMINACY (Example 4, phi_pi = 0.8) ---');
disp(' ');
disp('Autoregressive matrix V1lubik:');
disp(Q4.V1lubik);
fprintf('  Eigenvalues of V1lubik: ');
eig_V4 = eig(Q4.V1lubik);
fprintf('%.4f  ', eig_V4);
fprintf('\n');

disp(' ');
disp('Shock impact matrix V2lubik_shocks (SYMBOLIC, contains free M1):');
disp(Q4.V2lubik_shocks);
disp('  The M1 entries are FREE parameters: each choice gives a different REE.');
disp('  Setting M1 = 0 gives the orthogonal solution (beliefs ignore shocks).');

% Evaluate V2 at M1 = 0 for comparison
M1_syms_ex6 = symvar(Q4.V2lubik_shocks);
V2_orth_ex6 = double(subs(Q4.V2lubik_shocks, M1_syms_ex6, zeros(1, length(M1_syms_ex6))));
disp(' ');
disp('V2lubik_shocks evaluated at M1 = 0 (orthogonal solution):');
disp(V2_orth_ex6);

disp(' ');
disp('Belief shock impact matrix V2lubik_belief:');
disp(double(Q4.V2lubik_belief));
disp('  This matrix maps sunspot/belief shocks into the state variables.');
disp(['  Degree of indeterminacy: ' num2str(Q4.indetdegree)]);
disp(['  Number of belief shocks: ' num2str(Q4.indetdegree)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 7: Comparing Different Shock Sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Small shock
std_small = [0.001 0.001];
Q7a = lubikschorfheide(EQ, var, shocks, std_small, [], 1, 40, 0);

% Large shock
std_large = [0.01 0.01];
Q7b = lubikschorfheide(EQ, var, shocks, std_large, [], 1, 40, 0);

figure('Name', 'Example 7: Shock Size Comparison');
subplot(1,2,1);
plot(1:40, Q7a.irf_lubik.eps_r(2, 2:end), 'b-', 'LineWidth', 2); hold on;
plot(1:40, Q7b.irf_lubik.eps_r(2, 2:end), 'r--', 'LineWidth', 2);
title('Inflation response to monetary shock');
legend('0.1% shock', '1% shock', 'Location', 'best');
xlabel('Periods'); grid on;

subplot(1,2,2);
plot(1:40, Q7a.irf_lubik.eps_r(2, 2:end)./std_small(1), 'b-', 'LineWidth', 2); hold on;
plot(1:40, Q7b.irf_lubik.eps_r(2, 2:end)./std_large(1), 'r--', 'LineWidth', 2);
title('Inflation response (normalized by shock size)');
legend('0.1% shock (normalized)', '1% shock (normalized)', 'Location', 'best');
xlabel('Periods'); grid on;

disp('Example 7: IRFs are linear in shock size (responses scale proportionally)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 8: Multiple Solutions Under Indeterminacy (Non-Orthogonal M1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
disp('=================================================================');
disp('Example 8: Exploring Alternative Solutions Under Indeterminacy');
disp('=================================================================');

% Use indeterminate model from Example 4
% Q4 has indeterminacy degree > 0, so V2lubik_shocks contains symbolic M1

if Q4.indetdegree == 0
    disp('Model is determinate. Skipping Example 9.');
else
    disp(['Degree of indeterminacy: ' num2str(Q4.indetdegree)]);
    disp(' ');

    % Extract symbolic M1 variables from V2lubik_shocks
    M1_syms = symvar(Q4.V2lubik_shocks);
    disp('Symbolic M1 matrix has dimensions:');
    disp(['  (' num2str(Q4.indetdegree) ' x ' num2str(length(shocks)) ')']);

    % Define different M1 matrices for different solutions:

    % Solution 1: Orthogonal solution (M1 = 0) - already computed by function
    % Interpretation: beliefs do NOT respond to fundamental shocks.
    M1_orthogonal = zeros(Q4.indetdegree, length(shocks));

    % Solution 2: Alternative arbitrary solution
    M1_alternative1 = 0.3 * ones(Q4.indetdegree, length(shocks));

    % Solution 3: Custom solution - beliefs respond asymmetrically
    M1_alternative2 = zeros(Q4.indetdegree, length(shocks));
    M1_alternative2(:, 1) = 0.8;   % Response to eps_r
    M1_alternative2(:, 2) = -0.2;  % Response to eps_y (can be negative!)

    % Get policy matrices
    V1 = Q4.V1lubik;
    V2_sym = Q4.V2lubik_shocks;  % Symbolic (contains M1)

    % Evaluate V2 for each solution
    V2_orth = double(subs(V2_sym, M1_syms, M1_orthogonal(:)'));
    V2_alt1 = double(subs(V2_sym, M1_syms, M1_alternative1(:)'));
    V2_alt2 = double(subs(V2_sym, M1_syms, M1_alternative2(:)'));

    % Simulate IRFs for monetary shock under each solution
    T = 40;  % IRF horizon
    nvar = length(var);
    shock_vec = zeros(length(shocks), 1);
    shock_vec(1) = std_shocks(1);  % Monetary shock

    % Initialize IRFs (in W1 coordinates)
    W1_orth = zeros(size(V1, 1), T+1);
    W1_alt1 = zeros(size(V1, 1), T+1);
    W1_alt2 = zeros(size(V1, 1), T+1);

    % Simulate
    for t = 1:T
        W1_orth(:, t+1) = V1 * W1_orth(:, t) + V2_orth * shock_vec;
        W1_alt1(:, t+1) = V1 * W1_alt1(:, t) + V2_alt1 * shock_vec;
        W1_alt2(:, t+1) = V1 * W1_alt2(:, t) + V2_alt2 * shock_vec;

        % Shock only in period 1
        if t == 1
            shock_vec = zeros(length(shocks), 1);
        end
    end

    % Transform back to original variables using Z matrix
    X_orth = Q4.Z * [W1_orth; zeros(nvar*2 - size(V1,1), T+1)];
    X_alt1 = Q4.Z * [W1_alt1; zeros(nvar*2 - size(V1,1), T+1)];
    X_alt2 = Q4.Z * [W1_alt2; zeros(nvar*2 - size(V1,1), T+1)];

    % Plot comparison: output and inflation responses
    figure('Name', 'Example 9: Multiple Solutions Under Indeterminacy');

    subplot(2,2,1);
    plot(0:T, X_orth(1, :), 'b-', 'LineWidth', 2); hold on;
    plot(0:T, X_alt1(1, :), 'r--', 'LineWidth', 2);
    plot(0:T, X_alt2(1, :), 'g-.', 'LineWidth', 2);
    title('Output response to monetary shock');
    xlabel('Periods'); ylabel('Deviation from SS');
    legend('Orthogonal (M1=0)', 'Alternative 1', 'Alternative 2', 'Location', 'best');
    grid on;

    subplot(2,2,2);
    plot(0:T, X_orth(2, :), 'b-', 'LineWidth', 2); hold on;
    plot(0:T, X_alt1(2, :), 'r--', 'LineWidth', 2);
    plot(0:T, X_alt2(2, :), 'g-.', 'LineWidth', 2);
    title('Inflation response to monetary shock');
    xlabel('Periods'); ylabel('Deviation from SS');
    legend('Orthogonal (M1=0)', 'Alternative 1', 'Alternative 2', 'Location', 'best');
    grid on;

    subplot(2,2,3);
    plot(0:T, X_orth(3, :), 'b-', 'LineWidth', 2); hold on;
    plot(0:T, X_alt1(3, :), 'r--', 'LineWidth', 2);
    plot(0:T, X_alt2(3, :), 'g-.', 'LineWidth', 2);
    title('Interest rate response to monetary shock');
    xlabel('Periods'); ylabel('Deviation from SS');
    legend('Orthogonal (M1=0)', 'Alternative 1', 'Alternative 2', 'Location', 'best');
    grid on;

    subplot(2,2,4);
    plot(0:T, X_orth(4, :), 'b-', 'LineWidth', 2); hold on;
    plot(0:T, X_alt1(4, :), 'r--', 'LineWidth', 2);
    plot(0:T, X_alt2(4, :), 'g-.', 'LineWidth', 2);
    title('Natural rate response to monetary shock');
    xlabel('Periods'); ylabel('Deviation from SS');
    legend('Orthogonal (M1=0)', 'Alternative 1', 'Alternative 2', 'Location', 'best');
    grid on;

    % Display interpretation
    disp(' ');
    disp('Interpretation of M1 and multiple solutions:');
    disp('-------------------------------------------');
    disp('Under indeterminacy, M1 is a FREE MATRIX that indexes different');
    disp('rational expectations equilibria. Each choice of M1 gives a valid REE.');
    disp(' ');
    disp('1. ORTHOGONAL SOLUTION (M1 = 0):');
    disp('   Beliefs do not respond to fundamental shocks.');
    disp('   Default in this code. Used in Lubik-Schorfheide (2003) Fig. 1.');
    disp(' ');
    disp('2. ALTERNATIVE SOLUTIONS (M1 ~= 0):');
    disp('   Beliefs respond to fundamental shocks according to M1.');
    disp('   Can amplify or dampen the effects of fundamental shocks.');
    disp('   M1 elements can be positive or negative.');
    disp(' ');
end
