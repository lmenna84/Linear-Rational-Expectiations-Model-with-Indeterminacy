%% EXAMPLES FOR chomoreno.m
% This script demonstrates how to use the chomoreno function
% to solve DSGE models using the forward method of Cho and Moreno (2011).
% Each example uses a basic New Keynesian model with different options.
% The examples parallel those in examples_lubikschorfheide.m so that
% users can compare the two solution methods.

clear
addpath('../functions')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 1: Basic 3-Equation NK Model (Linear, Determinacy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define symbolic variables using the eval-loop pattern
% (must match the pattern used internally by chomoreno.m)
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
shocks = [eps_r eps_y];
std_shocks = [0.0025 0.005];  % 25bp monetary shock, 50bp demand shock

% Solve (linear model, no steady state needed)
Q1 = chomoreno(EQ, var, shocks, std_shocks, [], 1);
close all

disp('Example 1: Basic NK model solved successfully');

% Also solve with lubikschorfheide for comparison
Q1_lubik = lubikschorfheide(EQ, var, shocks, std_shocks, [], 1, 40, 0);
close all

% Compare IRFs between the two methods
% (V1lubik operates in the transformed QZ space, so direct matrix
%  comparison is not meaningful. IRFs are in the original variable space.)
fprintf('\nComparison with lubikschorfheide (determinacy):\n');
irf_cm_r  = Q1.irf.eps_r(1:4, 2:end);
irf_ls_r  = Q1_lubik.irf_lubik.eps_r(1:4, 2:end);
irf_cm_y  = Q1.irf.eps_y(1:4, 2:end);
irf_ls_y  = Q1_lubik.irf_lubik.eps_y(1:4, 2:end);
fprintf('  Max |IRF_chomoreno - IRF_lubik| (eps_r): %.2e\n', ...
    max(abs(irf_cm_r - irf_ls_r), [], 'all'));
fprintf('  Max |IRF_chomoreno - IRF_lubik| (eps_y): %.2e\n', ...
    max(abs(irf_cm_y - irf_ls_y), [], 'all'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 2: Nonlinear NK Model with Automatic Linearization (linear=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A simple New Keynesian model written in levels with nonlinear terms.
% chomoreno automatically linearizes around the provided steady state.
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
Q2 = chomoreno(EQ_nl, var_nl, shocks_nl, std_nl, ss_nl, 0, 40, 1);
close all

disp('Example 2: Nonlinear NK model with linearization solved');

% --- Comparison plot: nonlinear (% deviations) vs linear (log-deviations) ---
% Example 1 variables are log-deviations from ss (standard log-linearized NK).
% Example 2 variables are absolute deviations from ss (chomoreno
% linearizes around ss). To convert to percentage deviations, divide by ss.
% Since y_ss = 1 and p_ss = 1, output and inflation IRFs coincide directly.
% For the interest rate, divide by ii_ss = 1/bet.

T_irf = 40;
t_plot = 1:T_irf;

% Example 1: log-deviation IRFs to eps_r (first 3 variables: y, p, ii)
irf1_y  = Q1.irf.eps_r(1, 2:end);
irf1_p  = Q1.irf.eps_r(2, 2:end);
irf1_ii = Q1.irf.eps_r(3, 2:end);

% Example 2: absolute deviation IRFs to eps_r, converted to % deviations
irf2_y  = log((Q2.irf.eps_r(1, 2:end)+y_ss) / y_ss);
irf2_p  = log((Q2.irf.eps_r(2, 2:end) +p_ss)/ p_ss);
irf2_ii = log((Q2.irf.eps_r(3, 2:end)+ii_ss) / ii_ss);

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

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 3: Custom IRF Length and No Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve with 60-period IRFs instead of default 40, no plots
Q3 = chomoreno(EQ, var, shocks, std_shocks, [], 1, 60, 0);

disp('Example 3: Solved with 60-period IRFs, no plots');
disp(['IRF size: ' num2str(size(Q3.irf.eps_r, 2)) ' periods']);

% Manually plot selected IRFs
figure('Name', 'Example 3: Custom IRF Plot');
subplot(2,2,1);
plot(1:60, Q3.irf.eps_r(1, 2:end), 'LineWidth', 2);
title('Output response to monetary shock');
xlabel('Periods'); ylabel('Deviation from SS');
grid on;

subplot(2,2,2);
plot(1:60, Q3.irf.eps_r(2, 2:end), 'LineWidth', 2);
title('Inflation response to monetary shock');
xlabel('Periods'); ylabel('Deviation from SS');
grid on;

subplot(2,2,3);
plot(1:60, Q3.irf.eps_y(1, 2:end), 'LineWidth', 2);
title('Output response to demand shock');
xlabel('Periods'); ylabel('Deviation from SS');
grid on;

subplot(2,2,4);
plot(1:60, Q3.irf.eps_y(2, 2:end), 'LineWidth', 2);
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

Q4 = chomoreno(EQ_indet, var, shocks, std_shocks, [], 1, 40, 1);
close all

disp('Example 4: Indeterminate model (passive monetary policy)');
disp('Under indeterminacy, chomoreno returns the unique forward (no-bubble) solution.');

% Compare with lubikschorfheide: find the M1 that reproduces the forward solution
Q4_lubik = lubikschorfheide(EQ_indet, var, shocks, std_shocks, [], 1, 40, 0);
close all
M1FS = find_M1_forward(Q4, Q4_lubik);

fprintf('\nM1 values that reproduce the forward solution:\n');
disp(M1FS);
fprintf('These are the belief parameters in the Lubik-Schorfheide framework\n');
fprintf('that make its solution coincide with the Cho-Moreno forward solution.\n');

% Compare IRFs
fprintf('\nIRF comparison (chomoreno vs lubik continuity solution):\n');
% Build continuity IRFs from lubik using M1FS
V1l = Q4_lubik.V1lubik;
V2l_sym = Q4_lubik.V2lubik_shocks;
M1_syms = symvar(V2l_sym);
V2l_cont = double(subs(V2l_sym, M1_syms, M1FS(:)'));
ZS = Q4_lubik.Z;
count_stable = size(V1l,1);
nvar = size(var,2);

% Simulate lubik continuity IRFs for eps_r
T = 40;
W1 = zeros(count_stable, T+1);
shock_vec = zeros(length(shocks),1);
shock_vec(1) = std_shocks(1);
for tt = 1:T
    W1(:, tt+1) = V1l * W1(:, tt) + V2l_cont * shock_vec;
    if tt == 1, shock_vec = zeros(length(shocks),1); end
end
X_lubik_cont = ZS * [W1; zeros(nvar*2-count_stable, T+1)];

fprintf('  Max |IRF_chomoreno - IRF_lubik_continuity| (eps_r, y): %.2e\n', ...
    max(abs(Q4.irf.eps_r(1,2:end) - X_lubik_cont(1,2:end))));
fprintf('  Max |IRF_chomoreno - IRF_lubik_continuity| (eps_r, p): %.2e\n', ...
    max(abs(Q4.irf.eps_r(2,2:end) - X_lubik_cont(2,2:end))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 5: Comparing Different Shock Sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Small shock
std_small = [0.001 0.001];
Q5a = chomoreno(EQ, var, shocks, std_small, [], 1, 40, 0);

% Large shock
std_large = [0.01 0.01];
Q5b = chomoreno(EQ, var, shocks, std_large, [], 1, 40, 0);

figure('Name', 'Example 5: Shock Size Comparison');
subplot(1,2,1);
plot(1:40, Q5a.irf.eps_r(2, 2:end), 'b-', 'LineWidth', 2); hold on;
plot(1:40, Q5b.irf.eps_r(2, 2:end), 'r--', 'LineWidth', 2);
title('Inflation response to monetary shock');
legend('0.1% shock', '1% shock', 'Location', 'best');
xlabel('Periods'); grid on;

subplot(1,2,2);
plot(1:40, Q5a.irf.eps_r(2, 2:end)./std_small(1), 'b-', 'LineWidth', 2); hold on;
plot(1:40, Q5b.irf.eps_r(2, 2:end)./std_large(1), 'r--', 'LineWidth', 2);
title('Inflation response (normalized by shock size)');
legend('0.1% shock (normalized)', '1% shock (normalized)', 'Location', 'best');
xlabel('Periods'); grid on;

disp('Example 5: IRFs are linear in shock size (responses scale proportionally)');

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 6: Unit Root in a Simple Small Open Economy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This example shows how eigenvalue methods and the forward method handle
% unit roots differently, using a simplified version of the Schmitt-Grohe
% and Uribe (2003) small open economy model.
%
% The representative household maximizes E_0 sum beta^t log(c_t)
% subject to c_t + B_t = R_{t-1}*B_{t-1} + y_bar
% where y_bar is a constant endowment, B_t is international bond holdings,
% and R_t is the (gross) world interest rate.
%
% Steady state: R_ss = 1/beta, B_ss = 0, c_ss = y_bar.
%
% Log-linearized Euler equation: c_f - c - R = 0
% Exogenous interest rate:       R - rho*R_l - eps = 0
%
% The model is subject to a unit root as discussed by Schmitt-Grohe and
% Uribe (2003).
%
% Key result:
%   - Cho-Moreno: DETERMINACY. The forward method pins down consumption
%     via the no-bubble condition: c_t = -R_t/(1-rho). V1 eigenvalues
%     are 0 and rho (both < 1). The unit root does not appear in V1; it
%     shows up in the F matrix (eigenvalue 1.0 <= threshold).
%   - Lubik-Schorfheide (and Dynare): INDETERMINACY. The QZ eigenvalues
%     include 1.0 from the unit root. With threshold 1+1e-6, this counts
%     as "stable," giving 3 stable eigenvalues vs 2 needed. The extra
%     degree of freedom allows permanent shifts in c (sunspot/belief
%     shocks).
%
% See also Cho and McCallum (2015, footnote 9): "There seems to be no
% general agreement as to whether this case [r(F)=1] is considered as
% indeterminate or determinate."

fprintf('\n=== Example 6: Unit Root in Simple SOE ===\n');

% --- Parameters ---
bet_ex6 = 0.99;
rho_ex6 = 0.8;

% --- Define variables ---
vars_ex6 = {'c6' 'R6'};
var_ex6   = zeros(0,0);
for jj = 1:size(vars_ex6,2)
    eval(['syms ' vars_ex6{jj}]);
    eval(['var_ex6 = [var_ex6 ' vars_ex6{jj} '];']);
    eval(['syms ' vars_ex6{jj} '_l']);
    eval(['syms ' vars_ex6{jj} '_f']);
end

syms eps6;
shocks_ex6 = [eps6];
std_ex6 = [0.01];

% Euler: c_f - c - R = 0
% R:     R - rho*R_l - eps = 0
EQ_ex6 = [c6_f - c6 - R6; ...
          R6 - rho_ex6*R6_l - eps6];

% --- Cho-Moreno ---
fprintf('\nCho-Moreno:\n');
Q6_cm = chomoreno(EQ_ex6, var_ex6, shocks_ex6, std_ex6, [], 1, 40, 1);
fprintf('V1 eigenvalues: ');
fprintf('%.6f  ', abs(eig(Q6_cm.V1)));
fprintf('\n');
fprintf('V1 =\n');
disp(Q6_cm.V1);
fprintf('V2 =\n');
disp(Q6_cm.V2);
fprintf('=> c_t = %.1f * R_t = -R_t/(1-rho). No unit root in the solution.\n', ...
    Q6_cm.V2(1)/std_ex6(1));
close all

% --- Lubik-Schorfheide ---
fprintf('\nLubik-Schorfheide:\n');
Q6_ls = lubikschorfheide(EQ_ex6, var_ex6, shocks_ex6, std_ex6, [], 1, 40, 1);
fprintf('Generalized eigenvalues (sorted): ');
fprintf('%.6f  ', Q6_ls.gen_eigs);
fprintf('\n');
fprintf('Indeterminacy degree: %d\n', Q6_ls.indetdegree);
fprintf('Need %d stable eigenvalues for determinacy; have %d.\n', ...
    size(var_ex6,2), sum(Q6_ls.gen_eigs <= 1+1e-6));
close all

% Compare IRFs: CM forward vs Lubik orthogonal (M1=0)
fprintf('\nIRF comparison (CM forward vs Lubik orthogonal):\n');
irf_cm6  = Q6_cm.irf.eps6(:, 2:end);
irf_ls6  = Q6_ls.irf_lubik.eps6(1:2, 2:end);
fprintf('  Max |diff|: %.2e\n', max(abs(irf_cm6 - irf_ls6), [], 'all'));

% Find M1 that matches forward solution
M1_ex6 = find_M1_forward(Q6_cm, Q6_ls);
fprintf('M1 (forward solution): ');
fprintf('%.6f  ', M1_ex6);
fprintf('\n');

% Plot: CM forward vs Lubik orthogonal + Lubik belief shock
t_ex6 = 1:40;
var_labels_ex6 = {'c (consumption)', 'R (interest rate)'};
figure('Name', 'Example 6: Unit Root in Simple SOE', ...
       'NumberTitle', 'off', 'Position', [100 100 900 400]);

% Fundamental shock IRFs
for jj = 1:2
    subplot(1,3,jj);
    hold on;
    plot(t_ex6, irf_cm6(jj,:), 'b-', 'LineWidth', 2);
    plot(t_ex6, irf_ls6(jj,:), 'r--', 'LineWidth', 1.5);
    title([var_labels_ex6{jj} ' (eps shock)']);
    xlabel('Periods'); ylabel('Deviation');
    if jj == 1, legend('CM forward', 'Lubik orth.', 'Location', 'best'); end
    grid on; hold off;
end

% Belief shock (only exists in Lubik, not in CM)
subplot(1,3,3);
belief_irf = Q6_ls.irf_lubik.belief1(1:2, 2:end);
hold on;
plot(t_ex6, belief_irf(1,:), 'b-', 'LineWidth', 2);
plot(t_ex6, belief_irf(2,:), 'r-', 'LineWidth', 2);
title('Belief shock (Lubik only)');
xlabel('Periods'); ylabel('Deviation');
legend('c', 'R', 'Location', 'best');
grid on; hold off;

sgtitle({'Example 6: CM finds Determinacy, Lubik finds Indeterminacy (degree 1)', ...
         'Unit root => Lubik allows permanent c shifts via belief shocks'});

close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example 7: Block-Recursive Structure with Budget Constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This example adds the linearized budget constraint to the model of
% Example 6:
%   (3) Budget: c + b - (1/beta)*b_l = 0
%
% where b is the bond position (linearized around B_ss = 0, notice b is deviation from steady state, not log-deviation).
%
% Crucially, b does NOT appear in equations (1)-(2): the core block
% (c, R) is self-contained. b is a "non-interacting" appendix variable.
% Including or excluding b should not affect the solution for c and R.
%
% The coefficient 1/beta > 1 on b_l means b has an explosive eigenvalue.
%
% Key result:
%   - Cho-Moreno WITHOUT b: Determinacy (as in Example 6)
%   - Cho-Moreno WITH b:    Instability (V1 eigenvalue 1/beta > 1)
%     But the core variables c, R have IDENTICAL IRFs either way.
%     The forward method is invariant to the non-interacting block.
%
%   - Lubik/Dynare WITHOUT b: Indeterminacy (as in Example 6)
%   - Lubik/Dynare WITH b:    DETERMINACY! The explosive 1/beta eigenvalue
%     provides the "missing" unstable root. Adding a non-interacting
%     equation flips the classification from indeterminacy to determinacy.
%
% This is the block-recursive structure phenomenon of Cho and McCallum
% (2015, Section 5, p.165): "the determinate solution may differ from
% the forward solution only when a model has a block-recursive structure."
%

fprintf('\n=== Example 7: Block-Recursive Structure (Budget Constraint) ===\n');

% --- Define variables (c, R, b) ---
vars_ex7 = {'c6' 'R6' 'b6'};
var_ex7   = zeros(0,0);
for jj = 1:size(vars_ex7,2)
    eval(['syms ' vars_ex7{jj}]);
    eval(['var_ex7 = [var_ex7 ' vars_ex7{jj} '];']);
    eval(['syms ' vars_ex7{jj} '_l']);
    eval(['syms ' vars_ex7{jj} '_f']);
end

% Euler:  c_f - c - R = 0
% R:      R - rho*R_l - eps = 0
% Budget: c + b - (1/beta)*b_l = 0
EQ_ex7 = [c6_f - c6 - R6; ...
          R6 - rho_ex6*R6_l - eps6; ...
          c6 + b6 - 1/bet_ex6*b6_l];

% --- Cho-Moreno WITH b ---
fprintf('\nCho-Moreno WITH b:\n');
Q7_cm = chomoreno(EQ_ex7, var_ex7, shocks_ex6, std_ex6, [], 1, 40, 1);
fprintf('V1 eigenvalues: ');
fprintf('%.6f  ', abs(eig(Q7_cm.V1)));
fprintf('\n');
close all

% Invariance check
fprintf('\nCho-Moreno invariance (c, R with vs without b):\n');
fprintf('  Max |V1(1:2,1:2) - V1_ex6|: %.2e\n', ...
    max(abs(Q7_cm.V1(1:2,1:2) - Q6_cm.V1), [], 'all'));
fprintf('  Max |V2(1:2,:) - V2_ex6|:   %.2e\n', ...
    max(abs(Q7_cm.V2(1:2,:) - Q6_cm.V2), [], 'all'));
fprintf('  Max |IRF(1:2,:) - IRF_ex6|: %.2e\n', ...
    max(abs(Q7_cm.irf.eps6(1:2,2:end) - Q6_cm.irf.eps6(:,2:end)), [], 'all'));
fprintf('  => Forward solution INVARIANT to including non-interacting b.\n');

% --- Lubik-Schorfheide WITH b ---
fprintf('\nLubik-Schorfheide WITH b:\n');
Q7_ls = lubikschorfheide(EQ_ex7, var_ex7, shocks_ex6, std_ex6, [], 1, 40, 1);
fprintf('Generalized eigenvalues (sorted): ');
fprintf('%.6f  ', Q7_ls.gen_eigs);
fprintf('\n');
fprintf('Indeterminacy degree: %d\n', Q7_ls.indetdegree);
fprintf('Need %d stable for determinacy; 1/beta = %.4f provides the missing unstable root.\n', ...
    size(var_ex7,2), 1/bet_ex6);
close all

% Classification comparison
fprintf('\n--- Classification summary ---\n');
fprintf('                   Without b       With b\n');
fprintf('  Cho-Moreno:      Determinacy     Instability (but c,R identical)\n');
fprintf('  Lubik/Dynare:    Indeterminacy   Determinacy\n');
fprintf('  => Adding a non-interacting equation changes Lubik classification!\n');

% Compare Lubik solutions: without b (indeterminacy) vs with b (determinacy)
fprintf('\nLubik IRF comparison (c and R: without b vs with b):\n');
irf_ls6_core = Q6_ls.irf_lubik.eps6(1:2, 2:end);  % indeterminacy, orthogonal
irf_ls7_core = Q7_ls.irf_lubik.eps6(1:2, 2:end);   % determinacy
for jj = 1:2
    fprintf('  %s: max|diff| = %.2e\n', var_labels_ex6{jj}, ...
        max(abs(irf_ls6_core(jj,:) - irf_ls7_core(jj,:))));
end

% Compare b IRF: CM (explosive) vs Lubik (stable)
fprintf('\nb IRF comparison (Cho-Moreno vs Lubik):\n');
b_irf_cm    = Q7_cm.irf.eps6(3, 2:end);
b_irf_lubik = Q7_ls.irf_lubik.eps6(3, 2:end);
fprintf('  CM b at t=40:    %.6f (explosive: eigenvalue 1/beta = %.4f)\n', ...
    b_irf_cm(end), 1/bet_ex6);
fprintf('  Lubik b at t=40: %.6f (stable: Lubik forces all variables stable)\n', ...
    b_irf_lubik(end));

% --- Plot ---
t_ex7 = 1:40;
figure('Name', 'Example 7: Block-Recursive Structure', ...
       'NumberTitle', 'off', 'Position', [100 100 1000 500]);

% c: CM with/without b + Lubik with/without b
subplot(1,3,1);
hold on;
plot(t_ex7, Q6_cm.irf.eps6(1, 2:end), 'b-', 'LineWidth', 2);
plot(t_ex7, Q7_cm.irf.eps6(1, 2:end), 'b--', 'LineWidth', 1.5);
plot(t_ex7, irf_ls6_core(1,:), 'r-', 'LineWidth', 2);
plot(t_ex7, irf_ls7_core(1,:), 'r--', 'LineWidth', 1.5);
title('c (consumption)');
xlabel('Periods'); ylabel('Deviation');
legend('CM no b', 'CM with b', 'Lubik no b (orth)', 'Lubik with b', 'Location', 'best');
grid on; hold off;

% R: similar
subplot(1,3,2);
hold on;
plot(t_ex7, Q6_cm.irf.eps6(2, 2:end), 'b-', 'LineWidth', 2);
plot(t_ex7, Q7_cm.irf.eps6(2, 2:end), 'b--', 'LineWidth', 1.5);
plot(t_ex7, irf_ls6_core(2,:), 'r-', 'LineWidth', 2);
plot(t_ex7, irf_ls7_core(2,:), 'r--', 'LineWidth', 1.5);
title('R (interest rate)');
xlabel('Periods'); ylabel('Deviation');
grid on; hold off;

% b: CM (explosive) vs Lubik (stable)
subplot(1,3,3);
hold on;
plot(t_ex7, b_irf_cm, 'b-', 'LineWidth', 2);
plot(t_ex7, b_irf_lubik, 'r-', 'LineWidth', 2);
title('b (bond position)');
xlabel('Periods'); ylabel('Deviation');
legend('CM (explosive)', 'Lubik (forced stable)', 'Location', 'best');
grid on; hold off;

sgtitle({'Example 7: Adding budget constraint (non-interacting b)', ...
         'CM invariant; Lubik flips from Indeterminacy to Determinacy'});

close all
