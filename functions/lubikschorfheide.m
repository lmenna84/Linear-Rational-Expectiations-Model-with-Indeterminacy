
function Q=lubikschorfheide(EQ,var,shocks,std_shocks,ss,linear,lenght_irf,plotirf,instsol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Lorenzo Menna %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Purpose:
% Solve DSGE models to first order using Lubik and Schorfheide (2003),
% "Computing sunspot equilibria in linear rational expectations models,"
% Journal of Economic Dynamics & Control, 28, 273-285.
% -----------------------------------
% Inputs:
% EQ = Nx1 vector containing the system of symbolic equations of the DSGE
% model, set equal to zero. Variables must be symbolic, and must add _f
% when forwarded and _l when lagged. Shocks must also appear as symbolic.
% Parameters must already have been given a value, so they should not
% appear as symbolic in the system.
% var = 1xN symbolic vector containing the N variables
% shocks = 1xM symbolic vector containing shocks
% std_shocks = 1xM vector containing the size of the impulse for each
% shock in the IRFs (e.g., one standard deviation; can be negative)
% Optional:
% ss = 1xN vector containing steady state values. Variables must be ordered
% as in var. Not optional if the model is not linear (linear=0). Can be
% left empty ([]) if the model is already linear (linear=1).
% linear = 1 if the model is already in linear form, 0 otherwise. Default
% is 0 (or 1 if ss is not provided).
% lenght_irf = scalar containing the desired length of the IRFs. Default is
% 40.
% plotirf = 1 if one wants to have the IRFs plotted, 0 otherwise. Default
% is 1.
% instsol = 1 if you want to force a solution under instability (by
% treating the least explosive roots as stable). Default is 0.
% -----------------------------------
% Returns:
% Q = structure containing:
%
% Q.V1lubik = count x count autoregressive matrix of the policy function
% in the stable block of the transformed system:
%   W1_t = V1lubik * W1_{t-1} + V2lubik_shocks * shocks + V2lubik_belief * beliefs
% where W = [W1; W2=0] = Z' * X and X = [Y; LAMBDA] is the augmented state
% vector (Y = original variables, LAMBDA = conditional expectations E_{t-1}[Y_t]).
% count = number of stable eigenvalues.
%
% Q.V2lubik_shocks = count x M matrix mapping fundamental shocks to the
% stable block W1. This matrix is a symbolic function of the free matrix M1,
% whose dimension is (k-r) x M, where k-r is the degree of indeterminacy.
% M1 indexes different solutions: M1=0 gives the orthogonal solution
% (beliefs do not respond to fundamental shocks).
%
% Q.V2lubik_belief = count x (k-r) matrix mapping sunspot/belief shocks
% to the stable block W1. Under determinacy (k-r=0) this matrix is empty.
%
% Q.Z = (2N) x (2N) matrix from the QZ decomposition, used to transform
% between the original variables X = [Y; LAMBDA] and the canonical W = Z' * X.
%
% Q.ETAstar = N x 1 symbolic vector containing the solution for the
% forecast errors eta as a function of fundamental shocks, the free
% matrix M1, and sunspot shocks (only returned if k-r > 0).
%
% Q.indetdegree = degree of indeterminacy (k - r). Zero under determinacy.
%
% Q.irf_lubik = structure containing the impulse responses. Has one field
% for each fundamental shock (named after the shock symbol) and one field
% for each belief shock (named belief1, belief2, ...). Each field is a
% (2N) x (lenght_irf+1) matrix; rows 1:N are the original variables,
% rows (N+1):2N are the conditional expectations. Column 1 is the initial
% state (zeros); columns 2:(lenght_irf+1) are periods 1:lenght_irf.
% The fundamental shock IRFs use the orthogonal solution (M1 = 0).
%
% Note: the IRFs plotted (when plotirf=1) use the orthogonal solution (OS),
% obtained by setting M1 = 0. Under orthogonality, beliefs do not respond
% to fundamental shocks. Other solutions (e.g., the continuity/forward
% solution) can be obtained by substituting the appropriate M1 values into
% the symbolic Q.V2lubik_shocks matrix.
% ------------------------------------


if nargin<9
    instsol=0;
end
if nargin<8
    instsol=0;
    plotirf=1;
end
if nargin<7
    instsol=0;
    plotirf=1;
    lenght_irf=40;
end
if nargin<6
    instsol=0;
    plotirf=1;
    linear=0;
    lenght_irf=40;
end
if nargin<5
    instsol=0;
    plotirf=1;
    linear=1;
    lenght_irf=40;
end

if isempty(instsol)==1
    instsol=0;
end

if isempty(lenght_irf)==1
    lenght_irf=40;
end
if isempty(linear)==1
    linear=0;
end
if isempty(plotirf)==1
    plotirf=1;
end

Y=var;
nvar=size(var,2);

for xx=1:nvar
    eval(['syms ' char(var(xx)) '_f;']);
    eval(['syms ' char(var(xx)) '_l;']);
    eval(['Y_f(1,' int2str(xx) ')=' char(var(xx)) '_f;']);
    eval(['Y_l(1,' int2str(xx) ')=' char(var(xx)) '_l;']);
end
for xx=1:nvar
    eval(['syms ' char(var(xx)) '_ss;']);
    eval(['Y_ss(1,' int2str(xx) ')=' char(var(xx)) '_ss;']);
end

variables=[Y Y_f Y_l shocks];

% --- Steady-state check (nonlinear models only) ---
% Verify that EQ(ss, ss, ss, 0) = 0 before linearizing.
if linear==0
    EQ_check = subs(EQ, Y, ss);
    EQ_check = subs(EQ_check, Y_f, ss);
    EQ_check = subs(EQ_check, Y_l, ss);
    EQ_check = subs(EQ_check, shocks, zeros(1, size(shocks,2)));
    ss_residual = double(EQ_check);
    if max(abs(ss_residual)) > 1e-6
        fprintf('\nSteady-state residuals:\n');
        for xx=1:length(ss_residual)
            fprintf('  Equation %d: %.2e\n', xx, ss_residual(xx));
        end
        error('lubikschorfheide:badSS', ...
            ['The provided steady state does not satisfy the model equations.\n' ...
             'Max absolute residual: %.2e.\n' ...
             'Ensure that EQ(ss, ss, ss, 0) = 0.'], max(abs(ss_residual)));
    end
end

if linear==0
    % Redefine all variables as deviations from the steady state
    for xx=1:nvar
        EQ=subs(EQ,Y(1,xx),Y_ss(1,xx)+Y(1,xx));
        EQ=subs(EQ,Y_l(1,xx),Y_ss(1,xx)+Y_l(1,xx));
        EQ=subs(EQ,Y_f(1,xx),Y_ss(1,xx)+Y_f(1,xx));
    end
end

if linear==0
    EQLin=taylor(EQ,variables,[zeros(1,nvar*3) zeros(1,size(shocks,2))],'Order',2);
else EQLin=EQ;
end

if linear==0
    EQLin=subs(EQLin,Y_ss,ss);
end

for xx=1:size(var,2)
    name_vars{xx}=char(var(xx));
end
for xx=1:size(shocks,2)
    name_shocks{xx}=char(shocks(xx));
end

for xx=1:lenght_irf
    for yy=1:size(shocks,2)
        shock_irf(yy,xx)=0;
    end
end

% --- Step 1: Augment the system with expectation errors ---
% Introduce conditional expectations LAMBDA_t = E_{t-1}[Y_t] and forecast
% errors ETA_t = Y_t - E_{t-1}[Y_t], so that Y_t = LAMBDA_{t-1} + ETA_t
% and E_t[Y_{t+1}] = LAMBDA_t. The augmented variable is [Y; LAMBDA] and the
% system is recast into the Sims (2000) canonical form:
%   LAMBDA0 * [Y; LAMBDA]_t = LAMBDA1 * [Y; LAMBDA]_{t-1} + PSI * shocks + PAI * ETA
% This representation is used by Lubik and Schorfheide (2003) to derive the
% full set of solutions under indeterminacy.
LAMBDA=zeros(0,0);
LAMBDA_l=zeros(0,0);
ETA=zeros(0,0);
for xx=1:size(Y,2)
    eval(['syms epsilon' char(Y(xx)) ';']);
    eval(['LAMBDA=[LAMBDA epsilon' char(Y(xx)) '];']);
    eval(['syms epsilon' char(Y(xx)) '_l;']);
    eval(['LAMBDA_l=[LAMBDA_l epsilon' char(Y(xx)) '_l];']);
    eval(['syms eta' char(Y(xx)) ';']);
    eval(['ETA=[ETA eta' char(Y(xx)) '];']);
end
EQ_indet=subs(EQLin,Y_f,LAMBDA);
EQ_indet=subs(EQ_indet,Y,LAMBDA_l+ETA);
EQ_add=zeros(0,0);
for xx=1:size(Y,2)
    EQ_add=[EQ_add;-ETA(xx)+Y(xx)-LAMBDA_l(xx)];
end
EQ_indet=[EQ_add;EQ_indet];
Y_new=[Y,LAMBDA];
Y_l_new=[Y_l,LAMBDA_l];
LAMBDA0=double(jacobian(EQ_indet,Y_new));
LAMBDA1=-double(jacobian(EQ_indet,Y_l_new));
PSI=-double(jacobian(EQ_indet,shocks));
PAI=-double(jacobian(EQ_indet,ETA));
% --- Step 2: QZ (generalized Schur) decomposition ---
% Decompose the pencil (LAMBDA0, LAMBDA1) via QZ: Q'*LAMBDA0*Z = LAMBDONA,
% Q'*LAMBDA1*Z = OMEGA, with LAMBDONA and OMEGA upper triangular.
% Generalized eigenvalues are diag(OMEGA)./diag(LAMBDONA).
% Reorder so that stable roots (|eigenvalue| <= 1) come first.
[LAMBDONA,OMEGA,QQ,Z]=qz(LAMBDA0,LAMBDA1);
eigs=ordeig(OMEGA,LAMBDONA);
[ordeigs,order1]=sort(abs(eigs)); 
for xx=1:size(eigs,1)
    if abs(eigs(xx))<=1+1e-6
        order(xx)=10;
    else order(xx)=1;
    end
end
[LAMBDONAS,OMEGAS,QS,ZS]=ordqz(LAMBDONA,OMEGA,QQ,Z,order');
count=0;
for xx=1:size(ordeigs,1)
    if abs(ordeigs(xx,1))<=1+1e-6
        count=count+1;
    end
end
% --- Step 3: Determinacy classification ---
% count = number of stable eigenvalues, n = N (number of original variables)
% count == n: Determinacy (unique stable solution)
% count <  n: Instability (no stable solution; too many explosive roots)
% count >  n: Indeterminacy (multiple stable solutions; sunspots matter)
if count==size(eigs,1)/2
    display('Determinacy');
elseif count<size(eigs,1)/2
    display('Instability');
    if instsol==0
        return
    elseif instsol==1
        threshold=ordeigs(size(eigs,1)/2)+0.000000000001;
        count=0;
        for xx=1:size(ordeigs,1)
            if abs(ordeigs(xx,1))<=threshold
                count=count+1;
            end
        end
    end
else display('Indeterminacy');
end
LAMBDONAS11=LAMBDONAS(1:count,1:count);
LAMBDONAS12=LAMBDONAS(1:count,count+1:size(LAMBDONAS,2));
LAMBDONAS22=LAMBDONAS(count+1:size(LAMBDONAS,1),count+1:size(LAMBDONAS,2));
OMEGAS11=OMEGAS(1:count,1:count);
OMEGAS12=OMEGAS(1:count,count+1:size(OMEGAS,2));
OMEGAS22=OMEGAS(count+1:size(OMEGAS,1),count+1:size(OMEGAS,2));
QS1=QS(1:count,:);
QS2=QS(count+1:size(QQ,1),:);
% --- Step 4: SVD of the impact matrix Q2*PAI ---
% The stability condition requires Q2*(PSI*eps + PAI*eta) = 0.
% The SVD of Q2*PAI determines how many independent restrictions the
% explosive block places on the forecast errors eta. The rank r of Q2*PAI
% gives the number of pinned-down directions; the degree of indeterminacy
% is k - r, where k = dim(eta). The columns of V2 (from the SVD) span the
% free directions along which sunspot/belief shocks can perturb eta.
mm=size(Y_new,2)-count;
kk=size(ETA,2);
[U,S,V]=svd(QS2*PAI);
minn=min(mm,kk);
Stemp=S(1:minn,1:minn);
count2=0;
for xx=1:minn
    if abs(Stemp(minn-xx+1,minn-xx+1))<1e-10
        count2=count2+1;
    end
end
rr=minn-count2;
S11=Stemp(1:rr,1:rr);
U1=U(1:mm,1:rr);
U2=U(1:mm,rr+1:size(U,2));
V1=V(1:kk,1:rr);
V2=V(1:kk,rr+1:size(V,2));
M1=zeros(0,0);
for yy=1:size(shocks,2)
    for xx=1:kk-rr
        eval(['syms M_' int2str(xx) '_' int2str(yy) ';']);
        eval(['M1=[M1 M_' int2str(xx) '_' int2str(yy) '];']);
    end
end
M1=reshape(M1,[kk-rr,size(shocks,2)]);
ZETAstar=zeros(0,0);
for xx=1:kk-rr
    eval(['syms ZETA_' int2str(xx) ';']);
    eval(['ZETAstar=[ZETAstar ZETA_' int2str(xx) '];']);
end
if kk-rr>0
    ETAstar=(-V1*inv(S11)*U1'*QS2*PSI+V2*M1)*shocks.'+V2*ZETAstar.';
end
% --- Step 5: Construct the full set of stable solutions (Proposition 1) ---
% The forecast errors are: eta = (-V1*D11^{-1}*U1'*Q2*PSI + V2*M1)*eps + V2*zeta
% where M1 is a (k-r) x l matrix that indexes each solution (free under
% indeterminacy), and zeta is a (k-r) x 1 sunspot/belief shock.
% The response of the canonical system to shocks and beliefs:
%   resp_shocks (depends on M1): total shock impact including forecast error adjustment
%   resp_belief: impact of sunspot shocks through the free directions V2
% The policy function in the stable block (transformed coordinates W1) is:
%   W1_t = V1lubik * W1_{t-1} + V2lubik_shocks * eps + V2lubik_belief * zeta
resp_shocks=PSI-PAI*V1*inv(S11)*U1'*QS2*PSI+PAI*V2*M1;
resp_belief=PAI*V2;
ZS11=ZS(1:count,1:count);
ZS21=ZS(count+1:size(ZS,1),1:count);
ZS22=ZS(count+1:size(ZS,1),count+1:size(ZS,2));
%V1lubik=[ZS11*inv(LAMBDONAS11)*OMEGAS11*inv(ZS11);ZS21*inv(ZS11)*ZS11*inv(LAMBDONAS11)*OMEGAS11*inv(ZS11)];
%V2lubik_shocks=[ZS11*inv(LAMBDONAS11)*QS1*resp_shocks;ZS21*inv(ZS11)*ZS11*inv(LAMBDONAS11)*QS1*resp_shocks];
%V2lubik_belief=[ZS11*inv(LAMBDONAS11)*QS1*resp_belief;ZS21*inv(ZS11)*ZS11*inv(LAMBDONAS11)*QS1*resp_belief];
V1lubik=inv(LAMBDONAS11)*OMEGAS11;
V2lubik_shocks=inv(LAMBDONAS11)*QS1*resp_shocks;
V2lubik_belief=inv(LAMBDONAS11)*QS1*resp_belief;
belief=zeros(0,0);
for xx=1:kk-rr
    eval(['syms belief' int2str(xx) ';']);
    eval(['belief=[belief belief' int2str(xx) '];']);
end
for xx=1:kk-rr
    name_belief{xx}=char(belief(xx));
end
for xx=1:size(var,2)
    name_vars{nvar+xx}=char(Y_new(nvar+xx));
end
for xx=1:size(shocks,2)
    eval(['Y_irf_indet.' char(shocks(1,xx)) '=zeros(nvar*2,lenght_irf);']);
    eval(['Y_irf_indet_temp.' char(shocks(1,xx)) '=zeros(nvar*2,lenght_irf);']);
end
for xx=1:kk-rr
    eval(['Y_irf_indet.' char(belief(1,xx)) '=zeros(nvar*2,lenght_irf);']);
    eval(['Y_irf_indet_temp.' char(belief(1,xx)) '=zeros(nvar*2,lenght_irf);']);
end
for xx=1:lenght_irf
    for yy=1:kk-rr
        belief_irf(yy,xx)=0;
    end
end
% --- Step 6: Compute impulse responses ---
% For fundamental shocks: use the orthogonal solution (M = 0), meaning
%   beliefs do not respond to fundamental shocks.
% For belief/sunspot shocks: apply a unit belief shock through V2lubik_belief.
% All IRFs are transformed back to original variables via the Z matrix.
for yy=1:size(shocks,2)
    shock_irf(:,1)=zeros(size(shocks,2),1);
    shock_irf(yy,1)=std_shocks(1,yy);
    for xx=1:lenght_irf
        eval(['Y_irf_indet.' char(shocks(1,yy)) '(1:count,xx+1)=V1lubik(1:count,:)*(Y_irf_indet.' char(shocks(1,yy)) '(1:count,xx))+subs(V2lubik_shocks(1:count,:),M1,zeros(size(M1,1),size(M1,2)))*(shock_irf(:,xx));']);
        %eval(['Y_irf_indet.' char(shocks(1,yy)) '(count+1:nvar*2,xx+1)=V1lubik(count+1:nvar*2,:)*(Y_irf_indet.' char(shocks(1,yy)) '(1:count,xx))+subs(V2lubik_shocks(count+1:nvar*2,:),M1,zeros(size(M1,1),size(M1,2)))*(shock_irf(:,xx));']);
        eval(['Y_irf_indet.' char(shocks(1,yy)) '(count+1:nvar*2,xx+1)=zeros(nvar*2-count,1);']);
        eval(['Y_irf_indet_temp.' char(shocks(1,yy)) '(:,xx+1)=ZS*Y_irf_indet.' char(shocks(1,yy)) '(:,xx+1);']);
    end
end
for yy=1:kk-rr
    belief_irf(:,1)=zeros(kk-rr,1);
    belief_irf(yy,1)=0.01;
    for xx=1:lenght_irf
        eval(['Y_irf_indet.' char(belief(1,yy)) '(1:count,xx+1)=V1lubik(1:count,:)*(Y_irf_indet.' char(belief(1,yy)) '(1:count,xx))+V2lubik_belief(1:count,:)*(belief_irf(:,xx));']);
        %eval(['Y_irf_indet.' char(belief(1,yy)) '(count+1:nvar*2,xx+1)=V1lubik(count+1:nvar*2,:)*(Y_irf_indet.' char(belief(1,yy)) '(1:count,xx))+V2lubik_belief(count+1:nvar*2,:)*(belief_irf(:,xx));']);
        eval(['Y_irf_indet.' char(belief(1,yy)) '(count+1:nvar*2,xx+1)=zeros(nvar*2-count,1);']);
        eval(['Y_irf_indet_temp.' char(belief(1,yy)) '(:,xx+1)=ZS*Y_irf_indet.' char(belief(1,yy)) '(:,xx+1);']);
    end
end
for xx=1:size(shocks,2)
    eval(['Y_irf_indet.' char(shocks(1,xx)) '=Y_irf_indet_temp.' char(shocks(1,xx)) ';']);
end
for xx=1:kk-rr
    eval(['Y_irf_indet.' char(belief(1,xx)) '=Y_irf_indet_temp.' char(belief(1,xx)) ';']);
end
if plotirf==1
    for xx=1:size(shocks,2)
        count=0;
       % figure(size(shocks,2)+xx)
       figure(xx)
        for yy=1:nvar*2
            count=count+1;
            subplot(ceil((2*nvar)^(1/2)),ceil((2*nvar)^(1/2)),count)
            eval(['plot(1:lenght_irf,Y_irf_indet.' char(shocks(1,xx)) '(' int2str(yy) ',2:lenght_irf+1));']);
          %  eval(['irf_lubik_' char(shocks(1,xx)) '=Y_irf_indet.' char(shocks(1,xx)) '(' int2str(yy) ',2:lenght_irf+1));']);
            eval(['title(''' name_vars{yy} '-' name_shocks{xx} '-OS'');']);
        end
    end
    for xx=1:kk-rr
        count=0;
       % figure(size(shocks,2)*2+xx)
       figure(size(shocks,2)+xx)
        for yy=1:nvar*2
            count=count+1;
            subplot(ceil((2*nvar)^(1/2)),ceil((2*nvar)^(1/2)),count)
            eval(['plot(1:lenght_irf,Y_irf_indet.' char(belief(1,xx)) '(' int2str(yy) ',2:lenght_irf+1));']);
       %     eval(['irf_lubik' char(belief(1,xx)) '=Y_irf_indet.' char(belief(1,xx)) '(' int2str(yy) ',2:lenght_irf+1));']);
            eval(['title(''' name_vars{yy} '-' name_belief{xx} '-OS'');']);
        end
    end
end
Q.irf_lubik=Y_irf_indet;
Q.V1lubik=V1lubik;
Q.V2lubik_shocks=V2lubik_shocks;
Q.V2lubik_belief=V2lubik_belief;
Q.Z=ZS;
if kk-rr>0
    Q.ETAstar=ETAstar;
end
Q.indetdegree=kk-rr;
Q.gen_eigs=ordeigs;


end
