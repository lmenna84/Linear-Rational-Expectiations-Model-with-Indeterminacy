
function Q=chomoreno(EQ,var,shocks,std_shocks,ss,linear,lenght_irf,plotirf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Author: Lorenzo Menna %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Purpose: 
% Solve DSGE models to first order using the forward method of Cho and
% Mccallum
% -----------------------------------
% Inputs:
% EQ = Nx1 vector containing the system of symbolic equations of the DSGE
% model. Variables must be symbolic, and must add _f when forwarded and _l
% when lagged. Shocks must also appear as symbolic. Parameters must already have been given a value, so they
% should not appear as symbolic in the system
% var = 1xN symbolic vector containing the N variables
% shocks = 1xM symbolic vector containing shocks
% std_shocks = 1xM vector containing shocks standard deviations
% Optional:
% ss = 1xN vector containing steady state values. Variables must be ordered
% as in var. It is not optional if the model is not linear.
% lenght_irf = scalar containing the desired lenght of the irfs. Default is
% 40.
% linear = 1 if the model is already in linear form, 0 otherwise. Default
% is 0.
% plotirf = 1 if one wants to have the irfs plotted, 0 otherwise. Default
% is 1.
% -----------------------------------
% Returns:
% Q = structure
% Q.V1 = NxN matrix of the policy function Y=V1*Y_l+V2*shocks
% Q.V2 = NxM matrix of the policy function Y=V1*Y_l+V2*shocks
% Q.var = NxN variance covariance matrix of the variables
% Q.corr = NxN correlation matrix of the variables
% Q.irf = structure containing as many fields as shocks. Each field
% contains a Nx(lenght_irf+1) matrix where each row contains the impulse
% response of variables ordered as in var to the shock
% Q.prova1 = NxN matrix containing residuals of the alghorithm which
% computes V1. Values must all be close to zero if the computation is
% correct.
% Q.prova2 = NxM matrix containing residuals of the alghorithm which
% computes V2. Values must all be close to zero if the computation is
% correct.
% Q.ss = 1xN vector containing steady state values. Present only if the
% model is in levels.
% Q.dev_ss = Nx1 vector containing equation residuals. All residuals must
% be close to zero if the steady state provided by the user is correct.
% ------------------------------------

threshold=1+1e-6;

if nargin<8
    plotirf=1;
end
if nargin<7
    plotirf=1;
    lenght_irf=40;
end
if nargin<6
    plotirf=1;
    linear=0;
    lenght_irf=40;
end
if nargin<5
    plotirf=1;
    linear=1;
    lenght_irf=40;
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
        error('chomoreno:badSS', ...
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

B=double(jacobian(EQLin,Y));
C=double(jacobian(EQLin,Y_l));
D=double(jacobian(EQLin,shocks));
E=double(subs(EQLin,variables,zeros(1,size(variables,2))));
A=double(jacobian(EQLin,Y_f));

% Now use the cho's method to obtain V1
% Rewrite the structural model as Y=Anew*Y_f+Bnew*Y_l+Cnew*Z
Anew=inv(B)*(-A);
Bnew=inv(B)*(-C);
Cnew=inv(B)*(-D);
Dnew=inv(B)*(-E);

beginVn(:,:,1)=Bnew;

distance2=1000;
countxx=1;

while distance2>0.00000001
LHS=[inv(eye(nvar,nvar)-Anew*beginVn(:,:,countxx))*Bnew];
beginVn(:,:,countxx+1)=double(LHS(1:nvar,:));
distance2=[abs(beginVn(:,:,countxx+1)-beginVn(:,:,countxx))];
distance2=double(distance2);
distance2=max(distance2);
distance2=max(distance2);
countxx=countxx+1;
if countxx>100000
    display('Model does not satisfy the FCC')
    return
end
end

V1n(:,:)=beginVn(:,:,countxx);
clear beginVn LHS

prova1=double([A*V1n(:,:)*V1n(:,:)+B*V1n(:,:)+C]);

% Check determinacy using Cho Mccallum
radius_sigma=max(abs(eig(V1n)));
F=inv(eye(nvar)-Anew*V1n)*Anew;
radius_F=max(abs(eig(F)));

if radius_sigma<1 & radius_F<=threshold;%1.000001
    display('Determinacy');
elseif radius_sigma<1 & radius_F>threshold;%1.000001
    display('Indeterminacy');
elseif radius_sigma>1
    display('Instability');
end

%rho=0;
beginV2n(:,:,1)=Cnew;

distance2=1000;
countxx=1;
while distance2>0.00000000001
LHS=[inv(eye(nvar,nvar)-Anew*V1n(:,:))*Cnew];%+inv(eye(nvar,nvar)-Anew*V1n(:,:))*(Anew*beginV2n(:,:,countxx))*rho];

beginV2n(:,:,countxx+1)=double(LHS(1:nvar,:));
distance2=[abs(beginV2n(:,:,countxx+1)-beginV2n(:,:,countxx))];
distance2=double(distance2);
distance2=max(distance2);
distance2=max(distance2);
countxx=countxx+1;
if countxx>100000
    break
end
end

V2n(:,:)=beginV2n(:,:,countxx);
clear beginV2n 

prova2=double([A*(V1n(:,:)*V2n(:,:))+B*V2n(:,:)+D]);

% IRF
for xx=1:size(shocks,2)
    eval(['Y_irf.' char(shocks(1,xx)) '=zeros(nvar,lenght_irf);']);
end

for xx=1:lenght_irf
    for yy=1:size(shocks,2)
        shock_irf(yy,xx)=0;
    end
end

for yy=1:size(shocks,2)
    shock_irf(:,1)=zeros(size(shocks,2),1);
    shock_irf(yy,1)=std_shocks(1,yy);
    for xx=1:lenght_irf
        eval(['Y_irf.' char(shocks(1,yy)) '(:,xx+1)=V1n(:,:)*(Y_irf.' char(shocks(1,yy)) '(:,xx))+V2n(:,:)*(shock_irf(:,xx));']);
    end
end

for xx=1:size(var,2)
    name_vars{xx}=char(var(xx));
end
for xx=1:size(shocks,2)
    name_shocks{xx}=char(shocks(xx));
end

if plotirf==1
    for xx=1:size(shocks,2)
        count=0;
        figure(xx)
        for yy=1:size(var,2)
            count=count+1;
            subplot(ceil(nvar^(1/2)),ceil(nvar^(1/2)),count)
            eval(['plot(1:lenght_irf,Y_irf.' char(shocks(1,xx)) '(' int2str(yy) ',2:lenght_irf+1));']);
            eval(['title(''' name_vars{yy} '-' name_shocks{xx} '-FS'');']);
        end
    end
end

% The theoretical variance matrix can be built as follows
% First build the variance covariance matrix of the shock multiplied by
% V2n.
varcov_shock=V2n*diag(std_shocks).^2*V2n';
% Using vectorization and kronecker products:
vec_varcov=inv(eye(nvar*nvar)-kron(V1n,V1n))*reshape(varcov_shock,[nvar*nvar,1]);
varcov=reshape(vec_varcov,[nvar,nvar]);

DD = sqrt(diag(varcov));
DD=diag(DD);
R=inv(DD)*varcov*inv(DD); 


Q.V1=V1n;
Q.V2=V2n;
Q.var=varcov;
Q.corr=R;
Q.irf=Y_irf;
Q.prova1=prova1;
Q.prova2=prova2;
if isempty(ss)==0
Q.ss=ss;
end
Q.dev_ss=E;

end
