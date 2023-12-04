%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimates the Pesendorfer and Schmidt-Dengler (2008) Model %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Allow setting c in scripts
keepvars = {'c', 'nstart'};
clearvars('-except', keepvars{:});

load('psd2008_data');

disp(' ')
disp('*********************************************')

disp(sprintf('Equilibrium: eqm_dgp = %d', Data.eqm_dgp))
disp(sprintf('Replications: nsims = %d', Data.nsims))
disp(sprintf('Observations: nobs = %d', Data.nobs))
disp(sprintf('Weight on noise in initial estimates: c = %4.2f', c))

rng(0, 'twister');

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Preliminaries        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Tolerances
tolouter = 1e-6;  % Tolerance for outer loops

% Solver options
tolsolver = 1e-8;
opts1 = optimoptions('fminunc','Display','off',...
    'OptimalityTolerance',tolsolver,'StepTolerance',tolsolver,...
    'FunctionTolerance', tolsolver);
opts2 = optimoptions('fminunc','Display','notify','SpecifyObjectiveGradient',true,...
    'OptimalityTolerance',tolsolver,'StepTolerance',tolsolver,...
    'FunctionTolerance', tolsolver);
opts3 = optimoptions('fminunc','Display','notify','SpecifyObjectiveGradient',true);

% Flow utility specification
u1_0 = [0; 0; 0.1; 0.1]; % u(a = 0) for agent 1
u2_0 = [0; 0.1; 0; 0.1]; % u(a = 0) for agent 2

% Number of starting values (first is estimated, remaining random)
if ~exist('nstart')
    nstart = 5; % Aguirregabiria and Marcoux (2019) use 5 starting values.
end
disp(sprintf('Number of starting values: %d', nstart));

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate Choice Probablities %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Weight on noise (c = 0: no additional noise)
if ~exist('c')
    c = 0;
end
disp(sprintf('Starting value weights: %4.2f * initial estimates + %4.2f * noise ', [1-c, c]));

% Storing initial probabilities
phat1_init = zeros(nstates, 2, nstart, nsims);
phat2_init = zeros(nstates, 2, nstart, nsims);

disp(' ')
disp('CCP Estimation...')

for j = 1:nsims
    for r = 1:nstart
        if r == 1
            % Get correct observations
            X = X_sims(:, j);
            A1 = A1_sims(:, j);
            A2 = A2_sims(:, j);

            % Frequency estimators
            for ii = 1:4
                instate = (X == ii);
                phat1_init(ii, 2, 1, j) = mean(A1(instate), 1);
                phat2_init(ii, 2, 1, j) = mean(A2(instate), 1);
            end

            % Add noise to estimated CCPs (avg of CCP's and uniform[0,1])
            phat1_init(:, 2, 1, j) = (1 - c) .* phat1_init(:, 2, 1, j) + c .* rand(nstates, 1);
            phat2_init(:, 2, 1, j) = (1 - c) .* phat2_init(:, 2, 1, j) + c .* rand(nstates, 1);
        else
            % Random starting values
            phat1_init(:, 2, r, j) = rand(nstates, 1);
            phat2_init(:, 2, r, j) = rand(nstates, 1);
        end
    end
end
phat1_init(:, 1, :, :) = 1 - phat1_init(:, 2, :, :);
phat2_init(:, 1, :, :) = 1 - phat2_init(:, 2, :, :);

% Correct for 0 or 1 values
myzero = 1e-9;
phat1_init(phat1_init < myzero) = myzero;
phat2_init(phat2_init < myzero) = myzero;
phat1_init(phat1_init > 1 - myzero) = 1 - myzero;
phat2_init(phat2_init > 1 - myzero) = 1 - myzero;

%% NPL Estimation %%

maxits = 100; %maximum number of iterations

theta_NPL = zeros(3,nsims);
theta_NPL_1step = zeros(3,nsims);  %store estimates from first iteration
theta_NPL_1step_all = zeros(3,nstart,nsims);  %store estimates from first iteration

phat1_NPL = zeros(nstates,2,nsims);
phat2_NPL = zeros(nstates,2,nsims);
phat1_NPL_1step = zeros(nstates,2,nsims);
phat2_NPL_1step = zeros(nstates,2,nsims);
phat1_NPL_1step_all = zeros(nstates,2,nstart,nsims);
phat2_NPL_1step_all = zeros(nstates,2,nstart,nsims);

its_NPL = zeros(nstart, nsims);

disp('NPL Estimation...')

tic;
parfor j = 1:nsims

    % Get correct observations
    X = X_sims(:,j);
    A1 = A1_sims(:,j);
    A2 = A2_sims(:,j);

    % Plug into Data 
    Data_sim = Data;
    Data_sim.beta = beta;
    Data_sim.u1_0 = u1_0;
    Data_sim.u2_0 = u2_0;
    Data_sim.X = X;
    Data_sim.A1 = A1;
    Data_sim.A2 = A2;

    % Store temporary estimates for each starting value
    fval_temp = zeros(nstart,1);
    theta_temp = zeros(3, nstart);
    phat1_temp = zeros(nstates, 2, nstart);
    phat2_temp = zeros(nstates, 2, nstart);

    % For each starting value
    for r = 1:nstart

        % Initialize diff and tol
        diff = 1;
        its = 0;

        % Initial guesses
        theta0 = theta_true;
        thetahat = zeros(size(theta0)); % Prevent Matlab warning about uninitialized temporary
        fval = Inf; % Prevent Matlab warning about uninitialized temporary
        phat1 = phat1_init(:,:,r,j);
        phat2 = phat2_init(:,:,r,j);

        % Estimation
        while diff > tolouter && its < maxits
            % Get choice-specific value functions making social surplus zero
            % psi1 = -log(phat1);
            % psi2 = -log(phat2);

            % Note: not actually psi . . . this is e(p)
            psi1 = normpdf(norminv(phat1))./(2.*phat1);
            psi2 = normpdf(norminv(phat2))./(2.*phat2);

            %u_i = h_i*theta
            h1 = [ones(nstates,1), phat2(:,2), [1;1;0;0]];
            h2 = [ones(nstates,1), phat1(:,2), [1;0;1;0]];
            
            %conditional transition matrices
            F1_0 = [phat2(:,1), phat2(:,2), zeros(nstates,2)]; %agent 1, a = 0
            F1_1 = [zeros(nstates,2), phat2(:,1), phat2(:,2)]; %agent 1, a = 1
            F2_0 = [phat1(:,1), zeros(nstates,1), phat1(:,2), zeros(nstates,1)]; %agent 2, a = 0
            F2_1 = [zeros(nstates,1), phat1(:,1), zeros(nstates,1), phat1(:,2)]; %agent 2, a = 1
            
            %unconditional transition matrix
            F = [phat1(:,1).*phat2(:,1), phat1(:,1).*phat2(:,2), ...
               phat1(:,2).*phat2(:,1), phat1(:,2).*phat2(:,2)];
            
            %objects for NPL approximation to V
            %u[i]_1-u[i]_0=H[i]*theta+z[i]
            DG = eye(nstates)-beta.*F;
            DGinv = inv(DG);
            H1 = h1+beta.*(F1_1-F1_0)*(DG\(phat1(:,2).*h1));
            z1 = -u1_0+beta.*(F1_1-F1_0)*(DG\(phat1(:,1).*u1_0+sum(phat1.*psi1,2)));
            H2 = h2+beta.*(F2_1-F2_0)*(DG\(phat2(:,2).*h2));
            z2 = -u2_0+beta.*(F2_1-F2_0)*(DG\(phat2(:,1).*u2_0+sum(phat2.*psi2,2)));
            Data_sim.DG = DG;
            Data_sim.DGinv = DGinv;
            Data_sim.H1 = H1;
            Data_sim.H2 = H2;
            Data_sim.z1 = z1;
            Data_sim.z2 = z2;
            Data_sim.phat1 = phat1;
            Data_sim.phat2 = phat2;
            Data_sim.psi1 = psi1;
            Data_sim.psi2 = psi2;
            Data_sim.u1_0 = u1_0;
            Data_sim.u2_0 = u2_0;
            Data_sim.nstates = nstates;

            % Estimation
            anon_ll = @(theta) LL_NPL(theta,Data_sim); %anonymous log-likelihood

            % Avoid any problems with random starting values
            try
                % thetahat = fminunc(anon_ll,theta0,opts1);  % no gradient
                [ thetahat, fval, exflag] = fminunc(anon_ll, theta0, opts2);  % with gradient
            catch
                warning(sprintf('Starting value %d resulted in an error.', r));
                thetahat = theta0;
                fval = Inf;
            end

            % Update diff and its
            diff = norm(thetahat-theta0, Inf);
            its = its + 1;
            
            % Update choice probs and theta0
            phat1(:,2) = normcdf(H1*thetahat+z1);
            phat1(:,1) = 1 - phat1(:,2);
            phat2(:,2) = normcdf(H2*thetahat+z2);
            phat2(:,1) = 1 - phat2(:,2);
            theta0 = thetahat;

            % Always take 1-NPL from initial consistent CCPs
            if its == 1
                if r == 1
                    theta_NPL_1step(:, j) = thetahat;
                    phat1_NPL_1step(:, :, j) = phat1;   %uses updated value of its
                    phat2_NPL_1step(:, :, j) = phat2;
                end
                theta_NPL_1step_all(:, r, j) = thetahat;
                phat1_NPL_1step_all(:, :, r, j) = phat1;
                phat2_NPL_1step_all(:, :, r, j) = phat2;
            end

        end

        % Store results for each starting value
        fval_temp(r) = fval;
        theta_temp(:,r) = thetahat;
        phat1_temp(:,:,r) = phat1;
        phat2_temp(:,:,r) = phat2;
        its_NPL(r,j) = its;

    end % End loop over multiple starting values

    % Calculate best overall
    [fval_best, r_best] = min(fval_temp);

    % Store estimates
    theta_NPL(:,j) = theta_temp(:,r_best);
    phat1_NPL(:,:,j) = phat1_temp(:,:,r_best);
    phat2_NPL(:,:,j) = phat2_temp(:,:,r_best);
end
time_NPL = toc / 60;

bias_NPL_1step = mean(theta_NPL_1step-theta_true,2);
MSE_NPL_1step = mean((theta_NPL_1step-theta_true).^2,2);
MAD_NPL_1step = mean(abs(theta_NPL_1step-theta_true),2);

MSEtot_NPL_1step = sum(MSE_NPL_1step);

bias_NPL = mean(theta_NPL-theta_true,2);
MSE_NPL = mean((theta_NPL-theta_true).^2,2);
MAD_NPL = mean(abs(theta_NPL-theta_true),2);

MSEtot_NPL = sum(MSE_NPL);

convfail_NPL = sum(sum(its_NPL == maxits)) / (nstart * nsims);
its_NPL = sum(its_NPL, 1); % Sum across starting values

%% EPL Estimation %%

maxits = 100; %maximum number of iterations

theta_EPL = zeros(3,nsims);
theta_EPL_1step = zeros(3,nsims);
phat1_EPL = zeros(nstates,2,nsims);
phat2_EPL = zeros(nstates,2,nsims);

its_EPL = zeros(nstart, nsims);

%for use in taking value function differences
diffmat = [-eye(nstates),eye(nstates)];

disp('EPL Estimation...')

tic;
parfor j = 1:nsims
    
    % Get correct observations
    X = X_sims(:,j);
    A1 = A1_sims(:,j);
    A2 = A2_sims(:,j);

    %plug into Data 
    Data_sim = Data
    Data_sim.beta = beta;
    Data_sim.u1_0 = u1_0;
    Data_sim.u2_0 = u2_0;
    Data_sim.X = X;
    Data_sim.A1 = A1;
    Data_sim.A2 = A2;

    % Store temporary estimates for each starting value
    fval_temp = zeros(nstart,1);
    theta_temp = zeros(3, nstart);
    phat1_temp = zeros(nstates, 2, nstart);
    phat2_temp = zeros(nstates, 2, nstart);

    % For each starting value
    for r = 1:nstart

        %initialize diff and tol
        diff = 1;
        its = 0;
        
        % Initial guesses
        theta0 = theta_NPL_1step_all(:,r,j);  % Single-iteration NPL estimate for starting value r
        thetahat = zeros(size(theta0)); % Prevent Matlab warning about uninitialized temporary
        fval = Inf; % Prevent Matlab warning about uninitialized temporary
        phat1 = phat1_init(:,:,r,j);
        phat2 = phat2_init(:,:,r,j);
        %phat1 = phat1_NPL(:,:,j);
        %phat2 = phat2_NPL(:,:,j);

        %note: not actually psi . . . this is e(p)
        psi1 = normpdf(norminv(phat1))./(2.*phat1);
        psi2 = normpdf(norminv(phat2))./(2.*phat2);
        
        %u_i = h_i*theta
        h1 = [ones(nstates,1), phat2(:,2), [1;1;0;0]];
        h2 = [ones(nstates,1), phat1(:,2), [1;0;1;0]];
        
        %conditional transition matrices
        F1_0 = [phat2(:,1), phat2(:,2), zeros(nstates,2)]; %agent 1, a = 0
        F1_1 = [zeros(nstates,2), phat2(:,1), phat2(:,2)]; %agent 1, a = 1
        F2_0 = [phat1(:,1), zeros(nstates,1), phat1(:,2), zeros(nstates,1)]; %agent 2, a = 0
        F2_1 = [zeros(nstates,1), phat1(:,1), zeros(nstates,1), phat1(:,2)]; %agent 2, a = 1
        
        %unconditional transition matrix
        F = [phat1(:,1).*phat2(:,1), phat1(:,1).*phat2(:,2), ...
           phat1(:,2).*phat2(:,1), phat1(:,2).*phat2(:,2)];
        DG = eye(nstates)-beta.*F;
        
        
        %     H1=h1+beta.*(F1_1-F1_0)*(DG\(phat1(:,2).*h1));
        %         z1=-u1_0+beta.*(F1_1-F1_0)*(DG\(phat1(:,1).*u1_0+sum(phat1.*psi1,2)));
        %(integrated) value functions
        V1 = DG\(phat1(:,1).*u1_0+phat1(:,2).*(h1*theta0)+sum(phat1.*psi1,2));
        V2 = DG\(phat2(:,1).*u2_0+phat2(:,2).*(h2*theta0)+sum(phat2.*psi2,2));
        %v's
        v1_0 = u1_0+beta.*(F1_0*V1);
        v1_1 = h1*theta0+beta.*(F1_1*V1);
        v1 = [v1_0, v1_1];
        v2_0 = u2_0+beta.*(F2_0*V2);
        v2_1 = h2*theta0+beta.*(F2_1*V2);
        v2 = [v2_0, v2_1];
        v = [v1_0;v1_1;v2_0;v2_1];

        %estimation
        while diff>tolouter && its<maxits
            
            %matrix to invert (dG/dv)
            DG = eye(nstates*4)-Gderiv(theta0,v,Data_sim);
            
            %u_i=h_i*theta
            h1 = [ones(nstates,1), phat2(:,2), [1;1;0;0]];
            h2 = [ones(nstates,1), phat1(:,2), [1;0;1;0]];

            %conditional transition matrices
            F1_0 = [phat2(:,1), phat2(:,2), zeros(nstates,2)]; %agent 1, a = 0
            F1_1 = [zeros(nstates,2), phat2(:,1), phat2(:,2)]; %agent 1, a = 1
            F2_0 = [phat1(:,1), zeros(nstates,1), phat1(:,2), zeros(nstates,1)]; %agent 2, a = 0
            F2_1 = [zeros(nstates,1), phat1(:,1), zeros(nstates,1), phat1(:,2)]; %agent 2, a = 1

            %terms for choice-specific value functions 
            %G(theta,v1,v2)=A*theta+b
            A = [zeros(nstates,3); h1; zeros(nstates,3); h2];
            b = [u1_0+beta.*(F1_0*normsurplus(v1)); beta.*(F1_1*normsurplus(v1)); ...
               u2_0+beta.*(F2_0*normsurplus(v2)); beta.*(F2_1*normsurplus(v2))];

            %v[i]_1-v[i]_0=H[i]*theta+z[i]
            Hfull = DG\A;
            H1 = diffmat*Hfull(1:nstates*2,:);
            H2 = diffmat*Hfull(nstates*2+1:end,:);

            zfull = DG\(b-v)+v;
            z1 = diffmat*zfull(1:nstates*2,:);
            z2 = diffmat*zfull(nstates*2+1:end,:);
            
            
            %objects for NPL approximation to V
            %v[i]_1-v[i]_0=H[i]*theta+z[i]
            %         DG=eye(nstates)-beta.*F;
            %         DGinv=inv(DG);
            %         H1=h1+beta.*(F1_1-F1_0)*(DG\(phat1(:,2).*h1));
            %         z1=-u1_0+beta.*(F1_1-F1_0)*(DG\(phat1(:,1).*u1_0+sum(phat1.*psi1,2)));
            %         H2=h2+beta.*(F2_1-F2_0)*(DG\(phat2(:,2).*h2));
            %         z2=-u2_0+beta.*(F2_1-F2_0)*(DG\(phat2(:,1).*u2_0+sum(phat2.*psi2,2)));
            Data_sim.DG = DG;
            Data_sim.DGinv = DGinv;
            Data_sim.Hfull = Hfull;
            Data_sim.H1 = H1;
            Data_sim.H2 = H2;
            Data_sim.zfull = zfull;
            Data_sim.z1 = z1;
            Data_sim.z2 = z2;
            Data_sim.phat1 = phat1;
            Data_sim.phat2 = phat2;
            Data_sim.psi1 = psi1;
            Data_sim.psi2 = psi2;
            Data_sim.u1_0 = u1_0;
            Data_sim.u2_0 = u2_0;
            Data_sim.nstates = nstates;
            
            %initial guess of theta used by optimizer
            %helps when obj. fcn is undefined at theta0
            tempH = [H1;H2];
            tempz = [z1;z2];
            vdiff0 = tempH*theta0+tempz;
            if sum(sum(([1-normcdf(vdiff0),normcdf(vdiff0)])>0.9999))==0
                theta00 = theta0;
            else
                theta00 = (tempH'*tempH)\(-tempH'*tempz);
            end

            % Estimation
            anon_ll = @(theta) LL_EPL(theta,Data_sim); %anonymous log-likelihood

            % Avoid any problems with random starting values
            try
                [thetahat, fval, exflag]=fminunc(anon_ll,theta00,opts2);  %with gradient
            catch
                warning(sprintf('Starting value %d resulted in an error.', r));
                thetahat = theta00;
                fval = Inf;
            end
         
            % Update diff and its
            diff = norm(thetahat-theta0,Inf);
            its = its+1;
            
            % Update everything
            v = Hfull*thetahat+zfull;
            v1_0 = v(1:4);
            v1_1 = v(5:8);
            v2_0 = v(9:12);
            v2_1 = v(13:16);
            v1 = [v1_0,v1_1];
            v2 = [v2_0,v2_1];
            phat1 = [1-normcdf(v1_1-v1_0),normcdf(v1_1-v1_0)];
            phat2 = [1-normcdf(v2_1-v2_0),normcdf(v2_1-v2_0)];

            %update theta
            theta0 = thetahat;

            if its == 1 & r == 1
                theta_EPL_1step(:,j) = thetahat;
            end

        end

        % Store results for each starting value
        fval_temp(r) = fval;
        theta_temp(:,r) = thetahat;
        phat1_temp(:,:,r) = phat1;
        phat2_temp(:,:,r) = phat2;
        its_EPL(r,j) = its;

    end % End loop over multiple starting values

    % Calculate best overall
    [fval_best, r_best] = min(fval_temp);

    % Store estimates
    theta_EPL(:,j) = theta_temp(:,r_best);
    phat1_EPL(:,:,j) = phat1_temp(:,:,r_best);
    phat2_EPL(:,:,j) = phat2_temp(:,:,r_best);
end
time_EPL = toc/60;

bias_EPL_1step = mean(theta_EPL_1step-theta_true,2);
MSE_EPL_1step = mean((theta_EPL_1step-theta_true).^2,2);
MAD_EPL_1step = mean(abs(theta_EPL_1step-theta_true),2);

MSEtot_EPL_1step = sum(MSE_EPL_1step);

bias_EPL = mean(theta_EPL-theta_true,2);
MSE_EPL = mean((theta_EPL-theta_true).^2,2);
MAD_EPL = mean(abs(theta_EPL-theta_true),2);

MSEtot_EPL = sum(MSE_EPL);

convfail_EPL = sum(sum(its_EPL == maxits)) / (nstart * nsims);
its_EPL = sum(its_EPL, 1); % Sum across starting values


%% Report Results %%

disp(sprintf('Equilibrium: eqm_dgp = %d', Data.eqm_dgp))
disp(sprintf('Replications: nsims = %d', Data.nsims))
disp(sprintf('Observations: nobs = %d', Data.nobs))
disp(sprintf('Weight on noise in initial estimates: c = %4.2f', c))

%Bias
title = {'Param','1-NPL','1-EPL','NPL','EPL'};
paramlist = {'theta_M','theta_C','theta_EC'};
results = [bias_NPL_1step,bias_EPL_1step,bias_NPL,bias_EPL];
disp(' ');
disp('Bias of Estimates');
disp([title;paramlist',num2cell(results)]);

%MSE
title = {'Param','1-NPL','1-EPL','NPL','EPL'};
paramlist = {'theta_M','theta_C','theta_EC'};
results = [MSE_NPL_1step,MSE_EPL_1step,MSE_NPL,MSE_EPL];
disp(' ');
disp('MSE of Estimates');
disp([title;paramlist',num2cell(results)]);

%MAD
title = {'Param','1-NPL','1-EPL','NPL','EPL'};
paramlist = {'theta_M','theta_C','theta_EC'};
results = [MAD_NPL_1step,MAD_EPL_1step,MAD_NPL,MAD_EPL];
disp(' ');
disp('MAD of Estimates');
disp([title;paramlist',num2cell(results)]);

%Iterations
title = {'Param','1-NPL','1-EPL','NPL','EPL'};
statlist = {'Mean','Median','Std. Dev.','IQR'};
results = [1,1,mean(its_NPL),mean(its_EPL);...
    1,1,median(its_NPL),median(its_EPL);...
    0,0,std(its_NPL),std(its_EPL);...
    0,0,iqr(its_NPL),iqr(its_EPL)];
disp(' ');
disp('Iterations');
disp([title;statlist',num2cell(results)]);

%Time and Convergence
title = {'Param','1-NPL','1-EPL','NPL','EPL'};
statlist = {'Time','Conv %'};
results = [0,0,time_NPL,time_EPL;...
    1,1,1-convfail_NPL,1-convfail_EPL];
disp(' ');
disp('Time and Non-Convergence');
disp([title;statlist',num2cell(results)]);

% Save workspace
save(sprintf('psd2008_estimate_eqm%d_N%d_nstart%d_c%4.2f.mat',Data.eqm_dgp,Data.nobs,nstart,c))
