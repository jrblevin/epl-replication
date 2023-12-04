%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pesendorfer and Schmidt-Denger 2010 Model Estimation  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

%parameter
theta=-2;  %in [-10,-1]

%eq'm choice probs (2 players)
p=1./(1-theta);  %symmetric eq'm
p1=p;
p2=p;

%outer tolerance for iterative estimation
tolouter=1e-6;

%%

%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate Observations %
%%%%%%%%%%%%%%%%%%%%%%%%%

%reset rng
rng(0,'twister');

%number of obs. and sims
nobs=10000; %divisible by 2
nsims=500;

%sample observations
temp=rand(nobs,nsims);
A=1.*(temp<=p);   %observations;
A1=A(1:nobs/2,:);  %player 1
A2=A(nobs/2+1:end,:);  %player 2

%%
%%%%%%%%%%%%
% Full MLE %
%%%%%%%%%%%%

%setup
opts1=optimoptions('fmincon','Display','off');
opts2=optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',true);

theta_MLE=zeros(nsims,1);
p_MLE=zeros(nsims,1);

%loop over simulations
tic
for ii=1:nsims

    %data for estimation
    a=A(:,ii);
    a1=A1(:,ii);
    a2=A2(:,ii);

    %starting guess
%     theta0=theta;
    p10=mean(a1);
    p20=mean(a2);
    theta0=((p10-1)./p20+(p20-1)./p10)./2;  %mean of implied thetas
    theta0(theta0>=-1)=-1.001;
    theta0(theta<=-10)=-9.999;

    %data to pass to fcn
    Data.a=a;

    %anonymous objective function
    an_ll = @(theta) LL_MLE(theta,Data);

    %Estimation
    %theta1=fmincon(an_ll,theta0,[],[],[],[],-10,-1,[],opts1);
    theta1=fmincon(an_ll,theta0,[],[],[],[],-10,-1,[],opts2);  %analytic gradient

    %store estimates
    theta_MLE(ii)=theta1;
    p_MLE(ii)=1./(1-theta1);

end
time_MLE=toc

mean_MLE=mean(theta_MLE)
MSE_MLE=mean((theta_MLE-theta).^2)

%%

%%%%%%%%%%%%%%%%%%%%
% EPL Estimation   %
%%%%%%%%%%%%%%%%%%%%

%options
opts1=optimoptions('fmincon','Display','off');
opts2=optimoptions('fmincon','Display','off','SpecifyObjectiveGradient',true);

%storing estimates
theta_EPL=zeros(nsims,1); %converged estimates
theta_EPL_1step=zeros(nsims,1);  %1-step estimates
p1_EPL=zeros(nsims,1);   %converged estimates
p2_EPL=zeros(nsims,1);
theta0_EPL=zeros(nsims,1);  %initial estimates
p10_EPL=zeros(nsims,1);
p20_EPL=zeros(nsims,1);

%matrices for building jacobian estimate
idmat=eye(2);      %identity matrix
mymat=[0,1,;... %placeholder matrix
    1,0];

%maximum EPL iterations
maxiter=20;
numiter_EPL=zeros(nsims,1);  %number of iterations
nonconv_EPL=zeros(nsims,1);  %will equal one if max iterations hit

%loop over simulations
tic
for ii=1:nsims

    %get observations
    a1=A1(:,ii);
    a2=A2(:,ii);

    %initial estimates of choice probs
    p10=mean(a1);
    p20=mean(a2);
    p10_EPL(ii)=p10;
    p20_EPL(ii)=p20;

    %initial theta and v's
    theta0=((p10-1)./p20+(p20-1)./p10)./2;  %mean of implied thetas
    theta0(theta0>=-1)=-0.999;
    theta0(theta<=-10)=-9.999;
    theta0_EPL(ii)=theta0;
    v10=theta0.*p20;
    v20=theta0.*p10;

    %collecting them
    Y0=[v10;v20];

    %diff and iter
    diff=1;
    iter=0;

    %Estimation loop
    while diff>tolouter && iter<maxiter

        %jacobian estimate
        GY=idmat-theta0.*mymat;

        %Y1=C+B.*theta;
        C=Y0-GY\Y0;
        B=GY\([v20;v10]+[1;1]);

        %stuff to pass thru to function
        Data.nobs=nobs;
        Data.a1=a1;
        Data.a2=a2;
        Data.C=C;
        Data.B=B;

        %define objective function
        an_ll=@(theta) LL_EPL(theta,Data);

        %Estimation
        %theta1=fmincon(an_ll,theta0,[],[],[],[],-10,-1,[],opts1);
        theta1=fmincon(an_ll,theta0,[],[],[],[],-10,-1,[],opts2); %analytic gradient

        %update everything
        diff=norm(theta1-theta0,Inf);
        iter=iter+1;
        Y0=C+B.*theta1;
        v10=Y0(1);
        v20=Y0(2);
        theta0=theta1;  %mean of implied thetas
        theta0(theta0>=-1)=-0.999;
        theta0(theta0<=-10)=-9.999;

        if iter==1
            theta_EPL_1step(ii)=theta1;
            diff=1; %must have at least 2 iterations
        end

    end

    %store estimates and flag
    theta_EPL(ii)=theta1;
    p1_EPL(ii)=1+v10;
    p2_EPL(ii)=1+v20;
    numiter_EPL(ii)=iter;
    nonconv_EPL(ii)=(iter==maxiter);

end
time_EPL=toc

mean_EPL=mean(theta_EPL)
MSE_EPL=mean((theta_EPL-theta).^2)
convfails_EPL=sum(nonconv_EPL)
maxits_EPL=max(numiter_EPL,[],1)

%%

%%%%%%%%%%%%%%%%%%%%
% NPL Estimation   %
%%%%%%%%%%%%%%%%%%%%

%options
opts1=optimoptions('fmincon','Display','off');

%storing estimates
theta_NPL=zeros(nsims,1); %converged estimates
theta_NPL_1step=zeros(nsims,1);  %1-step estimates
p1_NPL=zeros(nsims,1);   %converged estimates
p2_NPL=zeros(nsims,1);
theta0_NPL=zeros(nsims,1);  %initial estimates
p10_NPL=zeros(nsims,1);
p20_NPL=zeros(nsims,1);

%maximum EPL iterations
maxiter=40;
numiter_NPL=zeros(nsims,1);  %number of iterations
nonconv_NPL=zeros(nsims,1);  %will equal one if max iterations hit

%loop over simulations
tic;
for ii=1:nsims

    %get observations
    a1=A1(:,ii);
    a2=A2(:,ii);

    %initial estimates of choice probs
    p10=mean(a1);
    p20=mean(a2);
    p10_NPL(ii)=p10;
    p20_NPL(ii)=p20;

    %initial theta and v's
    theta0=((p10-1)./p20+(p20-1)./p10)./2;  %mean of implied thetas
    theta0(theta0>=-1)=-0.999;
    theta0(theta<=-10)=-9.999;
    theta0_NPL(ii)=theta0;

    %diff and iter
    diff=1;
    iter=0;

    %Estimation loop
    while diff>tolouter && iter<maxiter

        %stuff to pass thru to function
        Data.nobs=nobs;
        Data.a1=a1;
        Data.a2=a2;
        Data.p10=p10;
        Data.p20=p20;

        %define objective function
        an_ll=@(theta) LL_NPL(theta,Data);

        %Estimation
        theta1=fmincon(an_ll,theta0,[],[],[],[],-10,-1,[],opts1);

        %update everything
        diff=norm(theta1-theta0,Inf);
        iter=iter+1;
        p11=1-myCDF(-theta1.*p20);
        p21=1-myCDF(-theta1.*p10);
        p10=p11;
        p20=p21;
        theta0=theta1;
        if iter==1
            theta_NPL_1step(ii)=theta1;
            diff=1; %must have at least 2 iterations
        end

    end

    %store estimates and flag
    theta_NPL(ii)=theta1;
    p1_NPL(ii)=p10;
    p2_NPL(ii)=p20;
    numiter_NPL(ii)=iter;
    nonconv_NPL(ii)=(iter==maxiter);

end
time_NPL=toc

mean_NPL=mean(theta_NPL)
MSE_NPL=mean((theta_NPL-theta).^2)
convfails_NPL=sum(nonconv_NPL)
maxits_NPL=max(numiter_NPL,[],1)


% Report results
title = {'Estimator','Mean','MSE'};
statlist = {'MLE','EPL','NPL'};
results = [
        mean_MLE, MSE_MLE;
        mean_EPL, MSE_EPL;
        mean_NPL, MSE_NPL;
];
disp(' ');
disp([title; statlist', num2cell(results)]);
