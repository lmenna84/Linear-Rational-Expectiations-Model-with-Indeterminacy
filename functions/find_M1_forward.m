function M1FS = find_M1_forward(Q_chomoreno, Q_lubik)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Lorenzo Menna %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Purpose:
% Given the output of chomoreno.m (forward solution) and lubikschorfheide.m
% (Lubik-Schorfheide solution under indeterminacy), this function finds the
% values of the belief matrix M1 that make the Lubik-Schorfheide solution
% coincide with the Cho-Moreno forward (continuity) solution.
% The computation is done via OLS regression.
% -----------------------------------
% Inputs:
% Q_chomoreno = output structure from chomoreno.m. Must contain:
%   Q_chomoreno.V2 = NxM forward solution shock impact matrix
% Q_lubik = output structure from lubikschorfheide.m (indeterminacy case).
%   Must contain:
%   Q_lubik.V2lubik_shocks = symbolic matrix (function of M1) of shock
%       responses in the Lubik-Schorfheide transformed system
%   Q_lubik.Z = transformation matrix Z from the ordered QZ decomposition
%   Q_lubik.indetdegree = degree of indeterminacy (kk-rr > 0 required)
% -----------------------------------
% Returns:
% M1FS = (indetdegree x M) matrix of belief parameters that reproduce
%   the forward solution within the Lubik-Schorfheide framework.
%   Each column corresponds to a shock.
% -----------------------------------

% Extract needed objects
V2n = Q_chomoreno.V2;
V2lubik_shocks = Q_lubik.V2lubik_shocks;
ZS = Q_lubik.Z;
indetdegree = Q_lubik.indetdegree;

nvar = size(V2n, 1);
nshocks = size(V2n, 2);

if indetdegree == 0
    error('find_M1_forward:determinacy', ...
        'The model is determinate (indetdegree=0). M1 does not exist.');
end

% Get symbolic M1 variables from V2lubik_shocks
M1 = sym(zeros(indetdegree, nshocks));
for yy = 1:nshocks
    for xx = 1:indetdegree
        eval(['syms M_' int2str(xx) '_' int2str(yy) ';']);
        eval(['M1(' int2str(xx) ',' int2str(yy) ')=M_' int2str(xx) '_' int2str(yy) ';']);
    end
end

% Number of stable eigenvalues (size of the W1 block)
count = size(V2lubik_shocks, 1);
ZS11 = ZS(1:count, 1:count);

% Transform Lubik-Schorfheide shock responses to original variable space
adj_V2lubik_shocks = ZS11 * V2lubik_shocks;

% Set up the regression: V2n = adj_V2lubik_shocks(1:nvar,:)
% This is linear in M1, so we use OLS shock by shock
for xx = 1:nshocks
    eqns(:,:,xx) = V2n(:,xx) - adj_V2lubik_shocks(1:nvar,xx);
end

M1FS = zeros(indetdegree, nshocks);
for xx = 1:nshocks
    regressori = jacobian(eqns(:,:,xx), M1(:,xx).');
    term_noto = subs(eqns(:,:,xx), M1(:,xx).', zeros(1, indetdegree));
    regressori = double(regressori);
    term_noto = double(term_noto);
    regressori = -regressori;
    M1FS(:,xx) = ((regressori'*regressori)\(regressori'*term_noto))';
end

% Verification: residual should be close to zero
prova3 = double(subs(adj_V2lubik_shocks(1:nvar,:), M1, M1FS) - V2n);
if max(abs(prova3(:))) > 1e-6
    warning('find_M1_forward:residual', ...
        'Max residual = %.2e. The forward solution may not be exactly reproducible.', ...
        max(abs(prova3(:))));
end

fprintf('M1 forward solution found. Max residual: %.2e\n', max(abs(prova3(:))));

end
