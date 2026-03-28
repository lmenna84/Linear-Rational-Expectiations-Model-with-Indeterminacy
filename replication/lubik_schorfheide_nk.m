%% lubik_schorfheide_nk.m
% =========================================================================
% Replicates impulse responses from:
%   Lubik, T.A. and Schorfheide, F. (2003), "Computing sunspot equilibria
%   in linear rational expectations models," Journal of Economic Dynamics
%   & Control, 28, 273-285.
%
% Model (Section 2 of the paper):
%   (1) E_t[y_{t+1}] + sigma*E_t[pi_{t+1}] = y_t + sigma*R_t   [IS curve]
%   (2) beta*E_t[pi_{t+1}] = pi_t - kappa*y_t                   [Phillips curve]
%   (3) R_t = psi*pi_t + eps_t                                   [Taylor rule]
%
% After substituting (3) into (1), the system reduces to 2 equations
% in 2 forward-looking variables (y, pi):
%   y_f + sigma*pi_f - y - sigma*psi*pi = sigma*eps
%   beta*pi_f - pi + kappa*y = 0
%
% Under a passive monetary policy rule (psi < 1), the model is
% indeterminate and sunspot shocks can affect equilibrium dynamics.
% Under an active rule (psi > 1), the equilibrium is determinate.
%
% This script replicates Figure 1 from the paper:
%   - Orthogonal solution (M1 = 0): dashed lines
%   - Continuity/forward solution (via Cho-Moreno): solid lines
%   - Sunspot shock (scaled by 0.5): dotted lines
%
% Parameters (Figure 1): beta=0.99, kappa=0.5, psi=0.95, sigma=1
% Shock: unanticipated 25-basis-point interest rate cut (eps = -0.25)
% =========================================================================

clear
addpath('../functions')

%% Parameters (from Figure 1 of the paper)
bet  = 0.99;    % discount factor (beta)
kap  = 0.5;     % Phillips curve slope (kappa)
psi0 = 0.95;    % Taylor rule coefficient on inflation (psi < 1 => indeterminacy)
sig  = 1;       % intertemporal elasticity of substitution (sigma)

%% Analytical eigenvalues (eq. 27 of the paper)
% lambda_1, lambda_2 = (1/2)(1 + (kappa*sigma+1)/beta)
%                       -/+ (1/2)*sqrt( ((kappa*sigma+1)/beta - 1)^2
%                                       + 4*kappa*sigma*(1-psi)/beta )
temp1 = (kap*sig + 1) / bet;
discriminant = (temp1 - 1)^2 + 4*kap*sig*(1 - psi0)/bet;
lambda1 = 0.5*(1 + temp1) - 0.5*sqrt(discriminant);
lambda2 = 0.5*(1 + temp1) + 0.5*sqrt(discriminant);

fprintf('\n--- Analytical eigenvalues (eq. 27) ---\n');
fprintf('lambda_1 = %.6f  (stable if |lambda_1| < 1)\n', lambda1);
fprintf('lambda_2 = %.6f  (stable if |lambda_2| < 1)\n', lambda2);
if psi0 < 1
    fprintf('psi = %.2f < 1: passive rule => INDETERMINACY expected\n', psi0);
    fprintf('  (lambda_1 < 1, lambda_2 > 1 => one explosive root, but k=2 > m=1)\n');
else
    fprintf('psi = %.2f >= 1: active rule => DETERMINACY expected\n', psi0);
end

%% Define symbolic variables
% y = output gap (log-deviation from steady state)
% p = inflation  (log-deviation from steady state)
vars = {'y' 'p'};

var   = zeros(0,0);
var_l = zeros(0,0);
var_f = zeros(0,0);
for ii = 1:size(vars,2)
    eval(['syms ' vars{ii}]);
    eval(['var   = [var '   vars{ii} '];']);
    eval(['syms ' vars{ii} '_l']);
    eval(['var_l = [var_l ' vars{ii} '_l];']);
    eval(['syms ' vars{ii} '_f']);
    eval(['var_f = [var_f ' vars{ii} '_f];']);
end

syms eps
shocks     = [eps];
std_shocks = [-0.25];  % 25-basis-point interest rate CUT (negative = expansionary)

%% Model equations (linear form, Taylor rule substituted out)
% After substituting R = psi*pi + eps into the IS equation:
%   Eq1 (IS):       y_f + sig*p_f - y - sig*psi*p - sig*eps = 0
%   Eq2 (Phillips): bet*p_f - p + kap*y = 0

EQ = [y_f + sig*p_f - y - sig*psi0*p - sig*eps; ...
      bet*p_f - p + kap*y];

%% Solve using Lubik-Schorfheide method
lenght_irf = 40;    % IRF horizon (matches Figure 1)
linear     = 1;     % model is already linear
plotirf    = 0;     % suppress default plots; we produce our own
instsol    = 0;

Q = lubikschorfheide(EQ, var, shocks, std_shocks, [], linear, lenght_irf, plotirf, instsol);

%% Solve using Cho-Moreno forward method
Q_cm = chomoreno(EQ, var, shocks, std_shocks, [], linear, lenght_irf, 0);

%% Find M1 values that reproduce the forward solution
% In this purely forward-looking model, the forward solution coincides
% with the continuity solution. This need not hold in general.
M1FS = find_M1_forward(Q_cm, Q);
fprintf('\nM1 (forward solution):\n');
disp(M1FS);

%% The forward solution IRFs from chomoreno
% In this model they also correspond to the continuity solution of the paper.

%% Extract IRFs
% The augmented system has 2*nvar = 4 variables: [y, p, xi_y, xi_p]
% Rows 1-2 are the structural variables (y, p)
% Column 1 = initial condition (zeros), columns 2:end = periods 1:lenght_irf

nvar = size(var, 2);
t    = 1:lenght_irf;

% --- Policy shock IRFs (orthogonal solution, M1 = 0) ---
irf_policy = Q.irf_lubik.eps;
y_orth = irf_policy(1, 2:end);   % output
p_orth = irf_policy(2, 2:end);   % inflation

% Reconstruct interest rate: R_t = psi*pi_t + eps_t (eps only at t=1)
R_orth    = psi0 * p_orth;
R_orth(1) = R_orth(1) + std_shocks(1);

% --- Sunspot/belief shock IRFs ---
if Q.indetdegree > 0
    % Identify the belief shock field (named 'belief1', 'belief2', etc.)
    fnames = fieldnames(Q.irf_lubik);
    belief_field = '';
    for ii = 1:length(fnames)
        if ~strcmp(fnames{ii}, char(shocks(1)))
            belief_field = fnames{ii};
            break;
        end
    end
    irf_sunspot_raw = Q.irf_lubik.(belief_field);

    % Rescale: the code uses a belief shock of 0.01 by default;
    % the paper plots the response to a reduced-form sunspot scaled by 0.5
    sunspot_scale = 0.5 / 0.01;

    y_sun = irf_sunspot_raw(1, 2:end) * sunspot_scale;
    p_sun = irf_sunspot_raw(2, 2:end) * sunspot_scale;
    R_sun = psi0 * p_sun;  % no fundamental shock in sunspot response

    fprintf('\nDegree of indeterminacy: %d\n', Q.indetdegree);
else
    fprintf('\nModel is determinate: no sunspot shocks.\n');
end

% --- Continuity/forward solution IRFs (from Cho-Moreno) ---
irf_cm = Q_cm.irf.eps;
y_cont = irf_cm(1, 2:end);   % output
p_cont = irf_cm(2, 2:end);   % inflation
R_cont    = psi0 * p_cont;
R_cont(1) = R_cont(1) + std_shocks(1);

%% Display numerical eigenvalues from the solution
fprintf('\n--- Eigenvalues of V1lubik (autoregressive matrix) ---\n');
eig_V1 = eig(double(Q.V1lubik));
for ii = 1:length(eig_V1)
    fprintf('  eig_%d = %.6f\n', ii, eig_V1(ii));
end

%% Plot Figure 1 (replication: orthogonality + continuity + sunspot)
figure('Name', 'Lubik-Schorfheide (2003) Figure 1', 'NumberTitle', 'off', ...
       'Position', [100 100 600 700])

% --- Response of Output ---
subplot(3,1,1)
hold on;
yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(t, y_cont, 'k-', 'LineWidth', 1.5);
plot(t, y_orth, 'k--', 'LineWidth', 1.5);
if Q.indetdegree > 0
    plot(t, y_sun, 'k:', 'LineWidth', 1.5);
    legend('Policy [Continuity]', 'Policy [Orthogonality]', 'Sunspot', 'Location', 'best');
end
title('Response of Output');
ylabel('y');
hold off;

% --- Response of Inflation ---
subplot(3,1,2)
hold on;
yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(t, p_cont, 'k-', 'LineWidth', 1.5);
plot(t, p_orth, 'k--', 'LineWidth', 1.5);
if Q.indetdegree > 0
    plot(t, p_sun, 'k:', 'LineWidth', 1.5);
    legend('Policy [Continuity]', 'Policy [Orthogonality]', 'Sunspot', 'Location', 'best');
end
title('Response of Inflation');
ylabel('pi');
hold off;

% --- Response of Interest Rate ---
subplot(3,1,3)
hold on;
yline(0, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
plot(t, R_cont, 'k-', 'LineWidth', 1.5);
plot(t, R_orth, 'k--', 'LineWidth', 1.5);
if Q.indetdegree > 0
    plot(t, R_sun, 'k:', 'LineWidth', 1.5);
    legend('Policy [Continuity]', 'Policy [Orthogonality]', 'Sunspot', 'Location', 'best');
end
title('Response of Interest Rate');
ylabel('R');
xlabel('Periods');
hold off;

sgtitle({'Lubik-Schorfheide (2003) Figure 1', ...
         ['beta = ' num2str(bet) ', kappa = ' num2str(kap) ...
          ', psi = ' num2str(psi0) ', sigma = ' num2str(sig)]});
