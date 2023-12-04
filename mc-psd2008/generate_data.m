%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate data for Pesendorfer and Schmidt-Dengler (2008) model    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%allow setting nsims, nobs, and eqm_dgp in scripts
keepvars = {'nsims', 'nobs', 'eqm_dgp'};
clearvars('-except', keepvars{:});

if ~exist('nsims')
    nsims=1000;  %number of simulations
end
if ~exist('nobs')
    nobs=1000;
end

nstates=4;

%parameters
theta_true=[1.2;-2.4;-0.2];
beta=0.9;


%states are arranged as follows (defined by LAST period's actions):
% 1: exit/exit (xx)
% 2: exit/enter (xe)
% 4: enter/exit (ex)
% 4: enter/enter (ee)


%action 1 is ENTER


%exit/entry probability by state

if ~exist('eqm_dgp')
    eqm_dgp=3; %in {1,2,3}
end

%Equilibrium 1
if eqm_dgp==1
    p1=1-[0.27;0.39;0.20;0.25];
    p1=[1-p1,p1];
    % p2=1-[0.72;0.78;0.58;0.71];
    p2=1-[0.72; 0.58; 0.78; 0.71];
    p2=[1-p2,p2];
end

%Equilibrium 2
if eqm_dgp==2
    p1=1-[0.38; 0.69; 0.17; 0.39];
    p1=[1-p1,p1];
    % p2=1-[0.47; 0.70; 0.16; 0.42];
    p2=1-[0.47; 0.16; 0.70; 0.42];
    p2=[1-p2,p2];
end

%Equilibrium 3
if eqm_dgp==3
    p1=1-[0.42; 0.70; 0.16; 0.41];
    p1=[1-p1,p1];
    p2=1-[0.42; 0.16; 0.70; 0.41];
    p2=[1-p2,p2];
end


%%

%%%%%%%%%%%%%%%%%%%%%
% Check Equilibrium %
%%%%%%%%%%%%%%%%%%%%%

%flow utility from exit
u1_0=[0;0;0.1;0.1]; %u(a=0) for agent 1
u2_0=[0;0.1;0;0.1]; %u(a=0) for agent 2

%u_i=h_i*theta
h1=[ones(nstates,1), p2(:,2), [1;1;0;0]];
h2=[ones(nstates,1), p1(:,2), [1;0;1;0]];

%data
Data.beta=beta;
Data.nstates=nstates;
Data.nobs=nobs;
Data.nsims=nsims;
Data.u1_0=u1_0;
Data.u2_0=u2_0;
Data.eqm_dgp=eqm_dgp;

%conditional transition matrices
F1_0=[p2(:,1), p2(:,2), zeros(nstates,2)]; %agent 1, a=0
F1_1=[zeros(nstates,2), p2(:,1), p2(:,2)]; %agent 1, a=1
F2_0=[p1(:,1), zeros(nstates,1), p1(:,2), zeros(nstates,1)]; %agent 2, a=0
F2_1=[zeros(nstates,1), p1(:,1), zeros(nstates,1), p1(:,2)]; %agent 2, a=1

%unconditional state transitions
FU=[p1(:,1).*p2(:,1), p1(:,1).*p2(:,2), p1(:,2).*p2(:,1), p1(:,2).*p2(:,2)];

% psi1=-log(p1);
% psi2=-log(p2);

%note: not actually psi . . . this is e(p)
psi1=normpdf(norminv(p1))./(2.*p1);
psi2=normpdf(norminv(p2))./(2.*p2);

%build H matrices
DG=eye(nstates)-beta.*FU;
DGinv=inv(DG);
H1=h1+beta.*(F1_1-F1_0)*(DG\(p1(:,2).*h1));
z1=-u1_0+beta.*(F1_1-F1_0)*(DG\(p1(:,1).*u1_0+sum(p1.*psi1,2)));
H2=h2+beta.*(F2_1-F2_0)*(DG\(p2(:,2).*h2));
z2=-u2_0+beta.*(F2_1-F2_0)*(DG\(p2(:,1).*u2_0+sum(p2.*psi2,2)));

%new choice probs
v1diff=H1*theta_true+z1;
v2diff=H2*theta_true+z2;
p1_check=[1-normcdf(v1diff), normcdf(v1diff)];
p2_check=[1-normcdf(v2diff), normcdf(v2diff)];


%%

%%%%%%%%%%%%%%%%%%%%
% Find Equilibrium %
%%%%%%%%%%%%%%%%%%%%

%finds eq'm entry probabilities

p_init=[p1(:,2);p2(:,2)];

p_eqm=fsolve(@(p) EqmDiff(p,theta_true,Data),p_init);

% p_eqm=fsolve(@(p) EqmDiff(p,theta_true,Data),0.3.*ones(size(p_init))); %finds equilibrium 3

% p_eqm=fminsearch(@(p) norm(EqmDiff(p,theta_true,Data),2),p_init);


%assign to choice probs
p1=[1-p_eqm(1:4),p_eqm(1:4)];
p2=[1-p_eqm(5:8),p_eqm(5:8)];


%u_i=h_i*theta
h1=[ones(nstates,1), p2(:,2), [1;1;0;0]];
h2=[ones(nstates,1), p1(:,2), [1;0;1;0]];

%conditional transition matrices
F1_0=[p2(:,1), p2(:,2), zeros(nstates,2)]; %agent 1, a=0
F1_1=[zeros(nstates,2), p2(:,1), p2(:,2)]; %agent 1, a=1
F2_0=[p1(:,1), zeros(nstates,1), p1(:,2), zeros(nstates,1)]; %agent 2, a=0
F2_1=[zeros(nstates,1), p1(:,1), zeros(nstates,1), p1(:,2)]; %agent 2, a=1

%unconditional state transitions
FU=[p1(:,1).*p2(:,1), p1(:,1).*p2(:,2), p1(:,2).*p2(:,1), p1(:,2).*p2(:,2)];

% psi1=-log(p1);
% psi2=-log(p2);
% psi1=p1.*norminv(p1)+normpdf(norminv(p1(:,1)));
% psi2=p2.*norminv(p2)+normpdf(norminv(p2(:,1)));

%note: not actually psi . . . this is e(p)
psi1=normpdf(norminv(p1))./(2.*p1);
psi2=normpdf(norminv(p2))./(2.*p2);

%build H matrices
DG=eye(nstates)-beta.*FU;
DGinv=inv(DG);
H1=h1+beta.*(F1_1-F1_0)*(DG\(p1(:,2).*h1));
z1=-u1_0+beta.*(F1_1-F1_0)*(DG\(p1(:,1).*u1_0+sum(p1.*psi1,2)));
H2=h2+beta.*(F2_1-F2_0)*(DG\(p2(:,2).*h2));
z2=-u2_0+beta.*(F2_1-F2_0)*(DG\(p2(:,1).*u2_0+sum(p2.*psi2,2)));


%useful stuff
%(integrated) value functions
V1=DG\(p1(:,1).*u1_0+p1(:,2).*(h1*theta_true)+sum(p1.*psi1,2));
V2=DG\(p2(:,1).*u2_0+p2(:,2).*(h2*theta_true)+sum(p2.*psi2,2));

%v's
v1_0=u1_0+beta.*(F1_0*V1);
v1_1=h1*theta_true+beta.*(F1_1*V1);
v1=[v1_0, v1_1];
v2_0=u2_0+beta.*(F2_0*V2);
v2_1=h2*theta_true+beta.*(F2_1*V2);
v2=[v2_0, v2_1];
v_true=[v1_0;v1_1;v2_0;v2_1];

%%

%generate data

%unconditional state transitions
FU=[p1(:,1).*p2(:,1), p1(:,1).*p2(:,2), p1(:,2).*p2(:,1), p1(:,2).*p2(:,2)];

%stationary distribution
[temp1,temp2]=eig(FU');
temp2=real(diag(temp2));
[B,I]=sort(temp2);
statdist=temp1(:,I(end))./sum(temp1(:,I(end)));
statdist_CDF=cumsum(statdist);
statedges=[0,statdist_CDF'];



%generate data
rng(0,'twister');

%for sampling
xdigi=(1:4);


%each observation we: 
% 1: draw a state from the stationary distribution
% 2: draw actions for players 1 and 2

%this makes the data i.i.d.


%draw states
x=randsample(xdigi,nobs*nsims,true,statdist');
X_sims=reshape(x,nobs,nsims);  %realized states

%draw actions
a1draws=rand(nobs,nsims);
a2draws=rand(nobs,nsims);

p11=p1(:,2);
p21=p2(:,2);

A1_sims=1.*((a1draws-p11(X_sims))<=0);  %agent 1 actions
A2_sims=1.*((a2draws-p21(X_sims))<=0);  %agent 2 actions

clear x a1draws a2draws B I temp1 temp2 statedges xdigi



%%
%save the data
save('psd2008_data');
