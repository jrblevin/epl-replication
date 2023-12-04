% main.m
%
% Adam Dearing and Jason Blevins
% Columbus, Ohio
% May 19, 2023

% ------------------------------------------------------------------------------
%                              MODEL
% ------------------------------------------------------------------------------
%
%  MAIN FEATURES OF THE MODEL
%
%  - Dynamic game of firm entry-exit in a market
%
%  -   i = firm index in { 1, 2, 3, 4, 5 }
%      t = time index
%
%  -   Decision variable:
%      a[it] = Indicator of the event 'firm i operates in the market at period t'
%
%          a[it] = 0 ----> Not active in market
%          a[it] = 1 ----> Active in market
%
%  -   The game is dynamic because there is a sunk cost of entry
%
%  -   The state variables of this game are the indicators of whether
%      the firms were active or not in the market at previous period,
%      and market size. These are payoff relevant state variables because
%      they determine whether a firm has to pay an entry cost to operate
%      in the market.
%
%  -   State variables
%      x[t] = ( s[t], a[1,t-1], a[2,t-1], a[3,t-1], a[4,t-1], a[5,t-1] )
%
%      e0[it] and e1[it] = Firm i's private information.
%
%  - Profit function
%
%      If a[it] = 0 (firm i is not active at period t)
%
%          Profit(0) = e0[it]
%
%      If a[it] = 1 (firm i is not active at period t)
%
%          Profit(1) = theta_fc_i - theta_ec * (1-a[i,t-1])
%                    + theta_rs * s[t]
%                    - theta_rn * ln(N[-it]) + 1) + e1[it]
%  where:
%          theta_fc_i, theta_ec, theta_rs, theta_rn are parameters
%
%          theta_fc_i  =   Firm i's fixed effect
%          theta_ec    =   Entry cost
%          N[-it]      =   Number of firms, other than i, operating in the market at period t.
%
%  eps_i(0) and eps_i(1) are private information shocks that are
%  i.i.d. over time, markets, and players, and independent of each other,
%  with Extreme Value Type 1 distribution.
%
%  - The following values are held fixed
%
%      Number of local markets (M) =   1,229 zip codes
%      Number of time periods  (T) =   24 years (1998-2021)
%      Number of players (N)       =   3
%      Discount factor             =   0.95
%      Support of s[t]             =   { 1, 2, 3, 4, 5 }
%      Transition matrix for s[t]  =   (defined below)
%
%  - We estimate the following parameters:
%
%      theta_fc_1
%      theta_fc_2
%      theta_fc_3
%      theta_ec
%      theta_rn
%      theta_rs

% -----------------------------------------------
%  VALUES OF PARAMETERS AND OTHER CONSTANTS
% -----------------------------------------------
nplayer = 3;    % Number of players
naction = 2;    % Number of actions
nstart = 5;     % Number of starting values
miniter = 3;    % Minimum number of NPL or EPL iterations
maxiter = 100;  % Maximum number of NPL or EPL iterations
disfact = 0.95; % Discount factor
sigmaeps = 1;   % Standard deviation of epsilon
kparam = 6;     % Number of parameters to estimate
geo = "county"; % County-level data
nrepli = 251;   % Bootstrap replications

% Points of support and transition probability of state variable s[t], market size
sval = [ 1:1:5 ]';  % Support of market size
numsval = size(sval, 1);  % Number of possible market sizes
nstate = numsval * (2^nplayer);  % Number of points in the state space

% Read transition matrix from file ptrans.txt
dat = readmatrix("ptrans.txt");
ptrans = dat(2:end,2:end);              % Drop column and row labels
ptrans = ptrans ./ sum(ptrans,2);       % Row normalize

% Vector with names of parameters of profit function
namesb = [ 'FC_SC'; 'FC_CC'; 'FC_BJ'; '   RS'; '   RN'; '   EC' ];

% Structure for storing parameters and settings
param.disfact = disfact;
param.sigmaeps = sigmaeps;
param.sval = sval;
param.ptrans = ptrans;
param.verbose = 0;

% State space
numa = 2^nplayer;
aval = zeros(numa, nplayer);
for i = 1:nplayer
  aval(:,i) = kron(kron(ones(2^(i-1),1), [ 0; 1 ]), ones(2^(nplayer - i), 1));
end
vstate = zeros(nstate, nplayer + 1);
vstate(:, 1) = kron(param.sval, ones(numa, 1));
vstate(:, 2:nplayer+1) = kron(ones(numsval, 1), aval);

% -----------------------------------------------------------
%  Load data
% -----------------------------------------------------------

dat = readmatrix(sprintf('clubstore_%s.csv', geo));
nobs = size(dat,1);
nmarket = size(unique(dat(:,1)),1); % Unique market indices
nyear = nobs / nmarket;

% --------------------------------------------------------------------------
%  Estimation
% --------------------------------------------------------------------------

% Cell arrays of starting values
start_name = {}; % Name
start_name_short = {}; % Short name for starting value
start_p = {}; % CCPs
start_theta = zeros(kparam, nstart); % Parameters
start_v = {}; % Value function
start_1npl = zeros(nstart, 1); % Whether or not to use 1-NPL values to initialize EPL
start_err = zeros(nstart, 1); % Failure
start_iter = zeros(nstart, 1); % Iteration count
start_fail = zeros(nstart, 1); % Failure
start_ll = zeros(nstart, 1);
start_theta1 = zeros(kparam, nstart); % 1-step estimates for each starting value
start_thetac = zeros(kparam, nstart); % Converged estimates for each starting value
start_v1npl = {}; % 1-NPL value function for each starting value
theta_seq = {}; % Sequences of estimates for each starting value
varb = {}; % Variance-covariance matrices for each starting value
start_times = {}; % Sequence of computational times per estimate

% Matrices of estimates for each iteration
bmat_1npl = zeros(nrepli,kparam); % 1-NPL estimates
bmat_2npl = zeros(nrepli,kparam); % 2-NPL estimates
bmat_3npl = zeros(nrepli,kparam); % 3-NPL estimates
bmat_cnpl = zeros(nrepli,kparam); % Converged NPL estimates
bmat_1epl = zeros(nrepli,kparam); % 1-EPL estimates
bmat_2epl = zeros(nrepli,kparam); % 2-EPL estimates
bmat_3epl = zeros(nrepli,kparam); % 3-EPL estimates
bmat_cepl = zeros(nrepli,kparam); % Converged EPL estimates

% Total computational times for each method
time_cnpl = zeros(nrepli, 1);
time_cepl = zeros(nrepli, 1);

% Per iteration times for each method
time_cnpl_iter = zeros(nrepli, 1);
time_cepl_iter = zeros(nrepli, 1);

% Counterfactuals
mean_nf_est_epl = zeros(nrepli, 1);
mean_nf_est_npl = zeros(nrepli, 1);
mean_nf_cf_epl = zeros(nrepli, 1);
mean_nf_cf_npl = zeros(nrepli, 1);
mean_entries_est_epl = zeros(nrepli, 1);
mean_entries_est_npl = zeros(nrepli, 1);
mean_entries_cf_epl = zeros(nrepli, 1);
mean_entries_cf_npl = zeros(nrepli, 1);
mean_exits_est_epl = zeros(nrepli, 1);
mean_exits_est_npl = zeros(nrepli, 1);
mean_exits_cf_epl = zeros(nrepli, 1);
mean_exits_cf_npl = zeros(nrepli, 1);
nf_s_est_epl = zeros(nplayer, numsval, nrepli);
nf_s_est_npl = zeros(nplayer, numsval, nrepli);
nf_s_cf_epl = zeros(nplayer, numsval, nrepli);
nf_s_cf_npl = zeros(nplayer, numsval, nrepli);
freq_active_s_est_epl = zeros(nplayer, numsval, nrepli);
freq_active_s_est_npl = zeros(nplayer, numsval, nrepli);
freq_active_s_cf_epl = zeros(nplayer, numsval, nrepli);
freq_active_s_cf_npl = zeros(nplayer, numsval, nrepli);

std_nf_est_epl = zeros(nrepli, 1);
std_nf_est_npl = zeros(nrepli, 1);
std_nf_cf_epl = zeros(nrepli, 1);
std_nf_cf_npl = zeros(nrepli, 1);
bareg_nf_est_epl = zeros(nrepli, 1);
bareg_nf_est_npl = zeros(nrepli, 1);
bareg_nf_cf_epl = zeros(nrepli, 1);
bareg_nf_cf_npl = zeros(nrepli, 1);
excess_est_epl = zeros(nrepli, 1);
excess_est_npl = zeros(nrepli, 1);
excess_cf_epl = zeros(nrepli, 1);
excess_cf_npl = zeros(nrepli, 1);
corr_ent_exit_est_epl = zeros(nrepli, 1);
corr_ent_exit_est_npl = zeros(nrepli, 1);
corr_ent_exit_cf_epl = zeros(nrepli, 1);
corr_ent_exit_cf_npl = zeros(nrepli, 1);

nf_hist_est_epl = zeros(nplayer+1, nrepli);
nf_hist_est_npl = zeros(nplayer+1, nrepli);
nf_hist_cf_epl = zeros(nplayer+1, nrepli);
nf_hist_cf_npl = zeros(nplayer+1, nrepli);

avgv_est_epl = zeros(nplayer, nrepli);
avgv_est_npl = zeros(nplayer, nrepli);
avgv_cf_epl = zeros(nplayer, nrepli);
avgv_cf_npl = zeros(nplayer, nrepli);

% Print a header
fprintf('=============================================================================================================\n');
fprintf('Club Store Entry/Exit Model of Dearing and Blevins (2024)\n');
fprintf('=============================================================================================================\n');
fprintf('%25s%10s\n', "Geography", geo);
fprintf('%25s%10d\n', "Markets", nmarket);
fprintf('%25s%10d\n', "Years", nyear);
fprintf('%25s%10d\n', "Observations", nobs);
fprintf('%25s%10d\n', "Replications", nrepli);
fprintf('\n');
fprintf('          Profit(1) = theta_fc_i - theta_ec * (1-a[i,t-1])\n');
fprintf('                    + theta_rs * s[t]\n');
fprintf('                    - theta_rn * ln(N[-it]) + 1) + e1[it]\n');

draw = 1;
while draw < nrepli + 1

    fprintf('=============================================================================================================\n');
    % Resample data after first replication
    if draw > 1
        fprintf("Bootstrap Replication %0d\n", draw-1);
        replacement = 1;
        markets = randsample(nmarket, nmarket, replacement);
        indices = [];
        % Calculate the indices for the observations for each market.
        for m = 1:nmarket
            range = 1:nyear;
            range = range + nyear*(markets(m)-1);
            indices = [ indices, range ];
        end
        dat_resamp = dat(indices,:);
    else
        fprintf('Estimation with Observed Sample\n');
        dat_resamp = dat;
    end
    fprintf('=============================================================================================================\n\n');

    market = dat_resamp(:,1);
    year = dat_resamp(:,2);
    aobs = dat_resamp(:,3:5);
    aobs_1 = dat_resamp(:,6:8);
    sobs = dat_resamp(:,9);

    % Calculate summary statistics of sample
    if (draw == 1)
        [ obs_nf, obs_entries, obs_exits, obs_nf_s, obs_freq_active_s, obs_std_nf, obs_bareg_nf, obs_excess, obs_corr_ent_exit, obs_nf_hist, obs_avgv ] = mpestat(aobs, aobs_1, sobs, 1, sval, nmarket);
    end

    if param.verbose > 0
        fprintf('-----------------------------------------------------------------------------------------\n');
        fprintf('         (a.1)   Estimation of initial CCPs (Semi-Parametric: Logit)\n');
    end
    % Construct dependent (aobsSP) and explanatory variables (xobsSP)
    aobsSP = reshape(aobs', nobs*nplayer, 1);
    alphai = kron(ones(nobs,1), eye(nplayer));
    xobsSP = kron(sobs, ones(nplayer,1));
    nfirms_1 = kron(sum(aobs_1, 2), ones(nplayer,1));
    aobs_1SP = reshape(aobs_1', nobs*nplayer, 1);
    xobsSP = [ alphai, xobsSP, aobs_1SP, nfirms_1 ];
    % Logit estimation
    [ best_logit, varest ] = milogit(aobsSP, xobsSP);
    % Construct probabilities
    vstateSP = [ ones(size(vstate,1), nplayer), vstate, ...
                 sum(vstate(:,2:nplayer+1)')' ];
    best = [ diag(best_logit(1:nplayer)) ; ...
             ones(1,nplayer) * best_logit(nplayer+1); ...
             eye(nplayer) * best_logit(nplayer+2); ...
             ones(1,nplayer) * best_logit(nplayer+3) ];
    v0SP = zeros(nstate, nplayer, naction);
    v0SP(:,:,2) = vstateSP*best;
    prob0SP = 1 ./ (1+exp(-vstateSP*best));

    start_name{1} = 'Semi-Parametric: Logit';
    start_name_short{1} = 'SP';
    start_p{1} = prob0SP;
    start_1npl(1) = 1;
    start_theta(:,1) = best_logit;
    start_v{1} = v0SP;

    if param.verbose > 0
        fprintf('-----------------------------------------------------------------------------------------\n');
        fprintf('         (a.2)   Perturbations of Semi-Parametric Logit CCPs (10)\n');
    end

    for r = 2:min(4, nstart)
        logit_perturb = best_logit + 0.2*randn(kparam, 1);
        best_perturb = [ diag(logit_perturb(1:nplayer)) ; ...
                         ones(1,nplayer) * logit_perturb(nplayer+1); ...
                         eye(nplayer) * logit_perturb(nplayer+2); ...
                         ones(1,nplayer) * logit_perturb(nplayer+3) ];
        prob0SP_perturb = 1 ./ (1+exp(-vstateSP*best_perturb));
        v0SP_perturb = zeros(nstate, nplayer, naction);
        v0SP_perturb(:,:,2) = vstateSP*best_perturb;

        start_name{r} = sprintf('Logit Perturbation %d', r - 2);
        start_name_short{r} = 'SP+Rnd';
        start_p{r} = prob0SP_perturb;
        start_1npl(r) = 0;
        start_theta(:,r) = logit_perturb;
        start_v{r} = v0SP_perturb;
    end

    if param.verbose > 0
        fprintf('-----------------------------------------------------------------------------------------\n');
        fprintf('         (a.3)   Estimation of initial CCPs (Non-Parametric)\n');
    end
    prob0NP = freqprob(aobs, [ sobs, aobs_1 ], vstate);

    start_name{5} = 'Non-Parametric';
    start_name_short{5} = 'NP';
    start_p{5} = prob0NP;
    start_1npl(5) = 1;

    if param.verbose > 0
        fprintf('-----------------------------------------------------------------------------------------\n');
        fprintf('         (a.4)   Completely random starting values...\n');
    end

    for r = 6:nstart
        random_theta = 0.25*randn(kparam, 1);
        random_mat = [ diag(random_theta(1:nplayer)) ; ...
                       ones(1,nplayer) * random_theta(nplayer+1); ...
                       eye(nplayer) * random_theta(nplayer+2); ...
                       ones(1,nplayer) * random_theta(nplayer+3) ];
        prob0R = 1 ./ (1+exp(-vstateSP*random_mat));
        v0R = zeros(nstate, nplayer, naction);
        v0R(:,:,2) = vstateSP*random_mat;

        start_name{r} = sprintf('Random #%d', r-4);
        start_name_short{r} = 'Rand';
        start_p{r} = prob0R;
        start_1npl(r) = 0;
        start_theta(:,r) = random_theta;
        start_v{r} = v0R;
    end

    %
    % NPL
    %

    for r = 1:nstart
        if param.verbose > 0
            fprintf('-----------------------------------------------------------------------------------------\n');
            fprintf('         (b.%d)   NPL Algorithm: %s\n', r, start_name{r});
        end
        [ theta_seq{r}, varb{r}, pcnpl{r}, start_v1npl{r}, vcnpl{r}, start_ll(r), start_iter(r), start_err(r), start_times{r} ] = npldygam(aobs, sobs, aobs_1, sval, ptrans, start_p{r}, disfact, miniter, maxiter, param.verbose);
    end

    max_llike = max(start_ll);
    best_r = find(start_ll == max_llike, 1, 'first');
    if best_r
        bmat_1npl(draw,:) = theta_seq{best_r}(1,:);
        bmat_2npl(draw,:) = theta_seq{best_r}(2,:);
        bmat_3npl(draw,:) = theta_seq{best_r}(3,:);
        bmat_cnpl(draw,:) = theta_seq{best_r}(start_iter(best_r),:);
        time_cnpl(draw,1) = sum(start_times{best_r});
        iter_cnpl(draw,1) = start_iter(best_r);
        time_cnpl_iter(draw,1) = time_cnpl(draw,1) / iter_cnpl(draw,1);
        fail_cnpl(draw,1) = 0;
        pequil_npl = pcnpl{best_r};
        vequil_npl = vcnpl{best_r};
    else
        fail_cnpl(draw,1) = 1;
    end

    fprintf('-----------------------------------------------------------------\n');
    fprintf('%10s', 'NPL');
    for r = 1:nstart
        fprintf(' %10s', start_name_short{r});
    end
    fprintf('\n');
    fprintf('-----------------------------------------------------------------\n');
    for k = 1:kparam
        fprintf('%10s', namesb(k,:))
        for r = 1:nstart
            if (start_iter(r) > 0)
                fprintf(' %10.4f', theta_seq{r}(start_iter(r),k));
            else
                fprintf(' %10s', '-');
            end
        end
        fprintf('\n');
    end
    fprintf('-----------------------------------------------------------------\n');
    fprintf('%10s', 'Iter');
    for r = 1:nstart
        fprintf(' %10d', start_iter(r));
    end
    fprintf('\n');

    fprintf('%10s', 'LL/n');
    for r = 1:nstart
        fprintf(' %10.3f', start_ll(r)/nobs);
    end
    fprintf('\n');

    fprintf('%10s', 'Best');
    for r = 1:nstart
        fprintf(' %10d', best_r == r);
    end
    fprintf('\n');

    fprintf('%10s', 'MaxIter');
    for r = 1:nstart
        fprintf(' %10d', start_iter(r) == maxiter);
    end
    fprintf('\n');

    fprintf('%10s', 'Error');
    for r = 1:nstart
        fprintf(' %10d', start_err(r));
    end
    fprintf('\n');

    %
    % EPL
    %

    for r = 1:nstart
        if param.verbose > 0
            fprintf('-----------------------------------------------------------------------------------------\n');
            fprintf('         (b.%d)   EPL Algorithm: %s\n', r, start_name{r});
        end
        if (start_1npl(r) > 0)
            v0 = start_v1npl{r};
            theta0 = start_theta1(:,r);
        else
            v0 = start_v{r};
            theta0 = start_theta(:,r);
        end
        [ theta_seq{r}, varb{r}, start_ll(r), start_iter(r), start_err(r), start_times{r}, pchoice{r}, v{r} ] = epldygam(aobs, sobs, aobs_1, sval, ptrans, theta0, v0, disfact, miniter, maxiter, param.verbose);
    end

    max_llike = max(start_ll);
    best_r = find(start_ll == max_llike, 1, 'first');
    if best_r
        bmat_1epl(draw,:) = theta_seq{best_r}(1,:);
        bmat_2epl(draw,:) = theta_seq{best_r}(2,:);
        bmat_3epl(draw,:) = theta_seq{best_r}(3,:);
        bmat_cepl(draw,:) = theta_seq{best_r}(start_iter(best_r),:);
        time_cepl(draw,1) = sum(start_times{best_r});
        iter_cepl(draw,1) = start_iter(best_r);
        time_cepl_iter(draw,1) = time_cepl(draw,1) / iter_cepl(draw,1);
        fail_cepl(draw,1) = 0;
        pequil_epl = pchoice{best_r};
        vequil_epl = v{best_r};
    else
        fail_cepl(draw,1) = 1;
    end

    fprintf('-----------------------------------------------------------------\n');
    fprintf('%10s', 'EPL');
    for r = 1:nstart
        fprintf(' %10s', start_name_short{r});
    end
    fprintf('\n');
    fprintf('-----------------------------------------------------------------\n');
    for k = 1:kparam
        fprintf('%10s', namesb(k,:))
        for r = 1:nstart
            if (start_iter(r) > 0)
                fprintf(' %10.4f', theta_seq{r}(start_iter(r),k));
            else
                fprintf(' %10s', '-');
            end
        end
        fprintf('\n');
    end
    fprintf('-----------------------------------------------------------------\n');
    fprintf('%10s', 'Iter');
    for r = 1:nstart
        fprintf(' %10d', start_iter(r));
    end
    fprintf('\n');

    fprintf('%10s', 'LL/n');
    for r = 1:nstart
        fprintf(' %10.3f', start_ll(r)/nobs);
    end
    fprintf('\n');

    fprintf('%10s', 'Best');
    for r = 1:nstart
        fprintf(' %10d', best_r == r);
    end
    fprintf('\n');

    fprintf('%10s', 'MaxIter');
    for r = 1:nstart
        fprintf(' %10d', start_iter(r) == maxiter);
    end
    fprintf('\n');

    fprintf('%10s', 'Error');
    for r = 1:nstart
        fprintf(' %10d', start_err(r));
    end
    fprintf('\n');
    fprintf('-----------------------------------------------------------------\n');

    % --------------------------------------------------------------------------
    % EPL Estimates: Simulate forward from initial conditions
    % --------------------------------------------------------------------------
    indobs = find(year==min(year));
    aobs_1_tmp = aobs_1(indobs,:);
    sobs_tmp = sobs(indobs,:);
    sim_aobs_1 = [];
    sim_sobs = [];
    sim_aobs = [];
    for t = 1:nyear
        [ aobs_out, sobs_out ] = forwardsim(aobs_1_tmp, sobs_tmp, pequil_epl, ptrans, sval, nobs);
        sim_aobs_1 = [ sim_aobs_1; aobs_1_tmp ];
        sim_sobs = [ sim_sobs; sobs_tmp ];
        sim_aobs = [ sim_aobs; aobs_out ];
        aobs_1_tmp = aobs_out;
        sobs_tmp = sobs_out;
    end

    vmat_epl = reshape(vequil_epl, [nstate, nplayer, naction]);
    [ mean_nf_est_epl(draw), mean_entries_est_epl(draw), mean_exits_est_epl(draw), nf_s_est_epl(:,:,draw), freq_active_s_est_epl(:,:,draw), std_nf_est_epl(draw), bareg_nf_est_epl(draw), excess_est_epl(draw), corr_ent_exit_est_epl(draw), nf_hist_est_epl(:,draw), avgv_est_epl(:,draw) ] =  mpestat(sim_aobs, sim_aobs_1, sim_sobs, param.verbose, sval, nmarket, vmat_epl, pequil_epl);
    nf_s_est_epl(:,:,draw) = nf_s_est_epl(:,:,draw) / nyear; % Divide observations by years to get #markets

    % -----------------------------------------------------------------------
    %  EPL Counterfactual: No competitive effect
    % -----------------------------------------------------------------------
    cf_param = bmat_cepl(draw,:)';
    cf_param(5) = 0;
    eqcond = @(v) Gfunc(sval, ptrans, cf_param, v, disfact, nplayer, naction);
    vin = reshape(vequil_epl, [nstate*nplayer*naction,1]);
    vequil = fsolve(eqcond, vin);
    vmat_cf = reshape(vequil, [nstate, nplayer, naction]);
    pequil_cf = exp(vmat_cf(:,:,2))./(exp(vmat_cf(:,:,1))+exp(vmat_cf(:,:,2)));

    % --------------------------------------------------------------------------
    % EPL Counterfactual: Simulate forward from initial conditions
    % --------------------------------------------------------------------------
    indobs = find(year==min(year));
    aobs_1_tmp = aobs_1(indobs,:);
    sobs_tmp = sobs(indobs,:);
    sim_aobs_1 = [];
    sim_sobs = [];
    sim_aobs = [];
    for t = 1:nyear
        [ aobs_out, sobs_out ] = forwardsim(aobs_1_tmp, sobs_tmp, pequil_cf, ptrans, sval, nobs);
        sim_aobs_1 = [ sim_aobs_1; aobs_1_tmp ];
        sim_sobs = [ sim_sobs; sobs_tmp ];
        sim_aobs = [ sim_aobs; aobs_out ];
        aobs_1_tmp = aobs_out;
        sobs_tmp = sobs_out;
    end

    [ mean_nf_cf_epl(draw), mean_entries_cf_epl(draw), mean_exits_cf_epl(draw), nf_s_cf_epl(:,:,draw), freq_active_s_cf_epl(:,:,draw), std_nf_cf_epl(draw), bareg_nf_cf_epl(draw), excess_cf_epl(draw), corr_ent_exit_cf_epl(draw), nf_hist_cf_epl(:,draw), avgv_cf_epl(:,draw) ] =  mpestat(sim_aobs, sim_aobs_1, sim_sobs, param.verbose, sval, nmarket, vmat_cf, pequil_cf);
    nf_s_cf_epl(:,:,draw) = nf_s_cf_epl(:,:,draw) / nyear; % Divide observations by years to get #markets

    % --------------------------------------------------------------------------
    % NPL Estimates: Simulate forward from initial conditions
    % --------------------------------------------------------------------------
    indobs = find(year==min(year));
    aobs_1_tmp = aobs_1(indobs,:);
    sobs_tmp = sobs(indobs,:);
    sim_aobs_1 = [];
    sim_sobs = [];
    sim_aobs = [];
    for t = 1:nyear
        [ aobs_out, sobs_out ] = forwardsim(aobs_1_tmp, sobs_tmp, pequil_npl, ptrans, sval, nobs);
        sim_aobs_1 = [ sim_aobs_1; aobs_1_tmp ];
        sim_sobs = [ sim_sobs; sobs_tmp ];
        sim_aobs = [ sim_aobs; aobs_out ];
        aobs_1_tmp = aobs_out;
        sobs_tmp = sobs_out;
    end

    vmat_npl = reshape(vequil_npl, [nstate, nplayer, naction]);
    [ mean_nf_est_npl(draw), mean_entries_est_npl(draw), mean_exits_est_npl(draw), nf_s_est_npl(:,:,draw), freq_active_s_est_npl(:,:,draw), std_nf_est_npl(draw), bareg_nf_est_npl(draw), excess_est_npl(draw), corr_ent_exit_est_npl(draw), nf_hist_est_npl(:,draw), avgv_est_npl(:,draw) ] =  mpestat(sim_aobs, sim_aobs_1, sim_sobs, param.verbose, sval, nmarket, vmat_npl, pequil_npl);
    nf_s_est_npl(:,:,draw) = nf_s_est_npl(:,:,draw) / nyear; % Divide observations by years to get #markets

    % -----------------------------------------------------------------------
    %  NPL Counterfactual: No competitive effect
    % -----------------------------------------------------------------------

    cf_param = bmat_cnpl(draw,:)';
    cf_param(5) = 0;
    eqcond = @(v) Gfunc(sval, ptrans, cf_param, v, disfact, nplayer, naction);
    vin = reshape(vequil_npl, [nstate*nplayer*naction,1]);
    vequil = fsolve(eqcond, vin);
    vmat_cf = reshape(vequil, [nstate, nplayer, naction]);
    pequil_cf = exp(vmat_cf(:,:,2))./(exp(vmat_cf(:,:,1))+exp(vmat_cf(:,:,2)));

    % --------------------------------------------------------------------------
    % NPL Counterfactual: Simulate forward from initial conditions
    % --------------------------------------------------------------------------
    indobs = find(year==min(year));
    aobs_1_tmp = aobs_1(indobs,:);
    sobs_tmp = sobs(indobs,:);
    sim_aobs_1 = [];
    sim_sobs = [];
    sim_aobs = [];
    for t = 1:nyear
        [ aobs_out, sobs_out ] = forwardsim(aobs_1_tmp, sobs_tmp, pequil_cf, ptrans, sval, nobs);
        sim_aobs_1 = [ sim_aobs_1; aobs_1_tmp ];
        sim_sobs = [ sim_sobs; sobs_tmp ];
        sim_aobs = [ sim_aobs; aobs_out ];
        aobs_1_tmp = aobs_out;
        sobs_tmp = sobs_out;
    end

    [ mean_nf_cf_npl(draw), mean_entries_cf_npl(draw), mean_exits_cf_npl(draw), nf_s_cf_npl(:,:,draw), freq_active_s_cf_npl(:,:,draw), std_nf_cf_npl(draw), bareg_nf_cf_npl(draw), excess_cf_npl(draw), corr_ent_exit_cf_npl(draw), nf_hist_cf_npl(:,draw), avgv_cf_npl(:,draw) ] =  mpestat(sim_aobs, sim_aobs_1, sim_sobs, param.verbose, sval, nmarket, vmat_cf, pequil_cf);
    nf_s_cf_npl(:,:,draw) = nf_s_cf_npl(:,:,draw) / nyear; % Divide observations by years to get #markets

    fprintf('\n');
    draw = draw + 1;

end

fprintf('=============================================================================================================\n');
fprintf('Final Results: 1-NPL and 1-EPL\n');
fprintf('=============================================================================================================\n');
fprintf('%25s%10s\n', "Geography", geo);
fprintf('%25s%10d\n', "Markets", nmarket);
fprintf('%25s%10d\n', "Years", nobs/nmarket);
fprintf('%25s%10d\n', "Observations", nobs);
fprintf('%25s%10d\n', "Bootstrap replications", nrepli);
fprintf('\n');
fprintf('          Profit(1) = theta_fc_i - theta_ec * (1-a[i,t-1])\n');
fprintf('                    + theta_rs * s[t]\n');
fprintf('                    - theta_rn * ln(N[-it]) + 1) + e1[it]\n');
fprintf('=============================================================================================================\n');
fprintf('%10s %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n', 'Param', 'NPL', 'Mean', 'Med', 'S.E.', 'EPL', 'Mean', 'Med', 'S.E.');
fprintf('-------------------------------------------------------------------------------------------------------------\n');
for k = 1:kparam
    fprintf('%10s %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            namesb(k,:), ...
            bmat_1npl(1,k), mean(bmat_1npl(2:nrepli,k)), median(bmat_1npl(2:nrepli,k)), std(bmat_1npl(2:end,k)), ...
            bmat_1epl(1,k), mean(bmat_1epl(2:nrepli,k)), median(bmat_1epl(2:nrepli,k)), std(bmat_1epl(2:end,k)));
end
fprintf('=============================================================================================================\n');
fprintf('Final Results: 2-NPL and 1-EPL\n');
fprintf('=============================================================================================================\n');
fprintf('%25s%10s\n', "Geography", geo);
fprintf('%25s%10d\n', "Markets", nmarket);
fprintf('%25s%10d\n', "Years", nobs/nmarket);
fprintf('%25s%10d\n', "Observations", nobs);
fprintf('%25s%10d\n', "Bootstrap replications", nrepli);
fprintf('\n');
fprintf('          Profit(1) = theta_fc_i - theta_ec * (1-a[i,t-1])\n');
fprintf('                    + theta_rs * s[t]\n');
fprintf('                    - theta_rn * ln(N[-it]) + 1) + e1[it]\n');
fprintf('=============================================================================================================\n');
fprintf('%10s %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n', 'Param', 'NPL', 'Mean', 'Med', 'S.E.', 'EPL', 'Mean', 'Med', 'S.E.');
fprintf('-------------------------------------------------------------------------------------------------------------\n');
for k = 1:kparam
    fprintf('%10s %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            namesb(k,:), ...
            bmat_2npl(1,k), mean(bmat_2npl(2:nrepli,k)), median(bmat_2npl(2:nrepli,k)), std(bmat_2npl(2:nrepli,k)), ...
            bmat_1epl(1,k), mean(bmat_1epl(2:nrepli,k)), median(bmat_1epl(2:nrepli,k)), std(bmat_1epl(2:nrepli,k)));
end
fprintf('=============================================================================================================\n');
fprintf('Final Results: Converged\n');
fprintf('=============================================================================================================\n');
fprintf('%25s%10s\n', "Geography", geo);
fprintf('%25s%10d\n', "Markets", nmarket);
fprintf('%25s%10d\n', "Years", nobs/nmarket);
fprintf('%25s%10d\n', "Observations", nobs);
fprintf('%25s%10d\n', "Bootstrap replications", nrepli);
fprintf('\n');
fprintf('          Profit(1) = theta_fc_i - theta_ec * (1-a[i,t-1])\n');
fprintf('                    + theta_rs * s[t]\n');
fprintf('                    - theta_rn * ln(N[-it]) + 1) + e1[it]\n');
fprintf('=============================================================================================================\n');
fprintf('%10s %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n', 'Param', 'NPL', 'S.E.', 'CI_LB', 'CI_UB', 'EPL', 'S.E.', 'CI_LB', 'CI_UB');
fprintf('-------------------------------------------------------------------------------------------------------------\n');
for k = 1:kparam
    fprintf('%10s %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            namesb(k,:), ...
            bmat_cnpl(1,k), std(bmat_cnpl(2:nrepli,k)), prctile(bmat_cnpl(2:nrepli,k),[2.5]), prctile(bmat_cnpl(2:nrepli,k),[97.5]), ...
            bmat_cepl(1,k), std(bmat_cepl(2:nrepli,k)), prctile(bmat_cepl(2:nrepli,k),[2.5]), prctile(bmat_cepl(2:nrepli,k),[97.5]));
end
fprintf('=====================================================================================================================================================\n');
fprintf('Estimated Simulations: EPL vs NPL\n');
fprintf('=====================================================================================================================================================\n');
fprintf('%25s%10s\n', "Geography", geo);
fprintf('%25s%10d\n', "Markets", nmarket);
fprintf('%25s%10d\n', "Years", nobs/nmarket);
fprintf('%25s%10d\n', "Observations", nobs);
fprintf('%25s%10d\n', "Bootstrap replications", nrepli);
fprintf('=====================================================================================================================================================\n');
fprintf('%30s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n', '', 'S=1', '', 'S=2', '', 'S=3', '', 'S=4', '', 'S=5', '');
fprintf('%30s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n', '', 'Mean', 'S.E.', 'Mean', 'S.E.', 'Mean', 'S.E.', 'Mean', 'S.E.', 'Mean', 'S.E.');
fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%30s %10.4f  %10.4f\n', "Mean Active Firms (EPL)", mean(mean_nf_est_epl), std(mean_nf_est_epl));
fprintf('%30s %10.4f  %10.4f\n', "Mean Active Firms (NPL)", mean(mean_nf_est_npl), std(mean_nf_est_npl));
fprintf('%30s %10.4f  %10.4f\n', "Mean Entries (EPL)", mean(mean_entries_est_epl), std(mean_entries_est_epl));
fprintf('%30s %10.4f  %10.4f\n', "Mean Entries (NPL)", mean(mean_entries_est_npl), std(mean_entries_est_npl));
fprintf('%30s %10.4f  %10.4f\n', "Mean Exits (EPL)", mean(mean_exits_est_epl), std(mean_exits_est_epl));
fprintf('%30s %10.4f  %10.4f\n', "Mean Exits (NPL)", mean(mean_exits_est_npl), std(mean_exits_est_npl));

fprintf('%30s %10.4f  %10.4f\n', "Mean Empty Markets (EPL)", mean(nf_hist_est_epl(1,:)), std(nf_hist_est_epl(1,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Empty Markets (NPL)", mean(nf_hist_est_npl(1,:)), std(nf_hist_est_npl(1,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Monopoly Markets (EPL)", mean(nf_hist_est_epl(2,:)), std(nf_hist_est_epl(2,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Monopoly Markets (NPL)", mean(nf_hist_est_npl(2,:)), std(nf_hist_est_npl(2,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Duopoly Markets (EPL)", mean(nf_hist_est_epl(3,:)), std(nf_hist_est_epl(3,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Duopoly Markets (NPL)", mean(nf_hist_est_npl(3,:)), std(nf_hist_est_npl(3,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Triopoly Markets (EPL)", mean(nf_hist_est_epl(4,:)), std(nf_hist_est_epl(4,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Triopoly Markets (NPL)", mean(nf_hist_est_npl(4,:)), std(nf_hist_est_npl(4,:)));

fprintf('%30s %10.4f  %10.4f\n', "Reg N[t] ~ N[t-1] (EPL)", mean(bareg_nf_est_epl), std(bareg_nf_est_epl));
fprintf('%30s %10.4f  %10.4f\n', "Reg N[t] ~ N[t-1] (NPL)", mean(bareg_nf_est_npl), std(bareg_nf_est_npl));
fprintf('%30s %10.4f  %10.4f\n', "Excess Turnover (EPL)", mean(excess_est_epl), std(excess_est_epl));
fprintf('%30s %10.4f  %10.4f\n', "Excess Turnover (NPL)", mean(excess_est_npl), std(excess_est_npl));
fprintf('%30s %10.4f  %10.4f\n', "Corr Ent Exit (EPL)", mean(corr_ent_exit_est_epl), std(corr_ent_exit_est_epl));
fprintf('%30s %10.4f  %10.4f\n', "Corr Ent Exit (NPL)", mean(corr_ent_exit_est_npl), std(corr_ent_exit_est_npl));

for i = 1:nplayer
    nf_s_est_epl_i = reshape(nf_s_est_epl(i,:,:), [numsval, nrepli]);
    nf_s_est_npl_i = reshape(nf_s_est_npl(i,:,:), [numsval, nrepli]);
    freq_active_s_est_epl_i = reshape(freq_active_s_est_epl(i,:,:), [numsval, nrepli]);
    freq_active_s_est_npl_i = reshape(freq_active_s_est_npl(i,:,:), [numsval, nrepli]);
    fprintf('%30s %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            sprintf("Markets Active %d (EPL)", i), ...
            mean(nf_s_est_epl_i(1,:)), std(nf_s_est_epl_i(1,:)), ...
            mean(nf_s_est_epl_i(2,:)), std(nf_s_est_epl_i(2,:)), ...
            mean(nf_s_est_epl_i(3,:)), std(nf_s_est_epl_i(3,:)), ...
            mean(nf_s_est_epl_i(4,:)), std(nf_s_est_epl_i(4,:)), ...
            mean(nf_s_est_epl_i(5,:)), std(nf_s_est_epl_i(5,:)));
    fprintf('%30s %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            sprintf("Markets Active %d (NPL)", i), ...
            mean(nf_s_est_npl_i(1,:)), std(nf_s_est_npl_i(1,:)), ...
            mean(nf_s_est_npl_i(2,:)), std(nf_s_est_npl_i(2,:)), ...
            mean(nf_s_est_npl_i(3,:)), std(nf_s_est_npl_i(3,:)), ...
            mean(nf_s_est_npl_i(4,:)), std(nf_s_est_npl_i(4,:)), ...
            mean(nf_s_est_npl_i(5,:)), std(nf_s_est_npl_i(5,:)));
    fprintf('%30s %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            sprintf("Frequency Active %d (EPL)", i), ...
            mean(freq_active_s_est_epl_i(1,:)), std(freq_active_s_est_epl_i(1,:)), ...
            mean(freq_active_s_est_epl_i(2,:)), std(freq_active_s_est_epl_i(2,:)), ...
            mean(freq_active_s_est_epl_i(3,:)), std(freq_active_s_est_epl_i(3,:)), ...
            mean(freq_active_s_est_epl_i(4,:)), std(freq_active_s_est_epl_i(4,:)), ...
            mean(freq_active_s_est_epl_i(5,:),"omitnan"), std(freq_active_s_est_epl_i(5,:),"omitnan"));
    fprintf('%30s %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            sprintf("Frequency Active %d (NPL)", i), ...
            mean(freq_active_s_est_npl_i(1,:)), std(freq_active_s_est_npl_i(1,:)), ...
            mean(freq_active_s_est_npl_i(2,:)), std(freq_active_s_est_npl_i(2,:)), ...
            mean(freq_active_s_est_npl_i(3,:)), std(freq_active_s_est_npl_i(3,:)), ...
            mean(freq_active_s_est_npl_i(4,:)), std(freq_active_s_est_npl_i(4,:)), ...
            mean(freq_active_s_est_npl_i(5,:),"omitnan"), std(freq_active_s_est_npl_i(5,:),"omitnan"));

    fprintf('%30s %10.4f  %10.4f\n', ...
            sprintf("Avg. Profit %d (EPL)", i), ...
            mean(avgv_est_epl(i,:)), std(avgv_est_epl(1,:)));
    fprintf('%30s %10.4f  %10.4f\n', ...
            sprintf("Avg. Profit %d (NPL)", i), ...
            mean(avgv_est_npl(i,:)), std(avgv_est_npl(1,:)));

end
fprintf('=====================================================================================================================================================\n');
fprintf('Counterfactual Simulations: EPL vs NPL\n');
fprintf('=====================================================================================================================================================\n');
fprintf('%25s%10s\n', "Geography", geo);
fprintf('%25s%10d\n', "Markets", nmarket);
fprintf('%25s%10d\n', "Years", nobs/nmarket);
fprintf('%25s%10d\n', "Observations", nobs);
fprintf('%25s%10d\n', "Bootstrap replications", nrepli);
fprintf('%25s%25s\n', "Counterfactual", "No competitive effect");
fprintf('=====================================================================================================================================================\n');
fprintf('%30s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n', '', 'S=1', '', 'S=2', '', 'S=3', '', 'S=4', '', 'S=5', '');
fprintf('%30s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s  %10s\n', '', 'Mean', 'S.E.', 'Mean', 'S.E.', 'Mean', 'S.E.', 'Mean', 'S.E.', 'Mean', 'S.E.');
fprintf('-----------------------------------------------------------------------------------------------------------------------------------------------------\n');
fprintf('%30s %10.4f  %10.4f\n', "Mean Active Firms (EPL)", mean(mean_nf_cf_epl), std(mean_nf_cf_epl));
fprintf('%30s %10.4f  %10.4f\n', "Mean Active Firms (NPL)", mean(mean_nf_cf_npl), std(mean_nf_cf_npl));
fprintf('%30s %10.4f  %10.4f\n', "Mean Entries (EPL)", mean(mean_entries_cf_epl), std(mean_entries_cf_epl));
fprintf('%30s %10.4f  %10.4f\n', "Mean Entries (NPL)", mean(mean_entries_cf_npl), std(mean_entries_cf_npl));
fprintf('%30s %10.4f  %10.4f\n', "Mean Exits (EPL)", mean(mean_exits_cf_epl), std(mean_exits_cf_epl));
fprintf('%30s %10.4f  %10.4f\n', "Mean Exits (NPL)", mean(mean_exits_cf_npl), std(mean_exits_cf_npl));

fprintf('%30s %10.4f  %10.4f\n', "Mean Empty Markets (EPL)", mean(nf_hist_cf_epl(1,:)), std(nf_hist_cf_epl(1,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Empty Markets (NPL)", mean(nf_hist_cf_npl(1,:)), std(nf_hist_cf_npl(1,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Monopoly Markets (EPL)", mean(nf_hist_cf_epl(2,:)), std(nf_hist_cf_epl(2,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Monopoly Markets (NPL)", mean(nf_hist_cf_npl(2,:)), std(nf_hist_cf_npl(2,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Duopoly Markets (EPL)", mean(nf_hist_cf_epl(3,:)), std(nf_hist_cf_epl(3,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Duopoly Markets (NPL)", mean(nf_hist_cf_npl(3,:)), std(nf_hist_cf_npl(3,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Triopoly Markets (EPL)", mean(nf_hist_cf_epl(4,:)), std(nf_hist_cf_epl(4,:)));
fprintf('%30s %10.4f  %10.4f\n', "Mean Triopoly Markets (NPL)", mean(nf_hist_cf_npl(4,:)), std(nf_hist_cf_npl(4,:)));

fprintf('%30s %10.4f  %10.4f\n', "Reg N[t] ~ N[t-1] (EPL)", mean(bareg_nf_cf_epl), std(bareg_nf_cf_epl));
fprintf('%30s %10.4f  %10.4f\n', "Reg N[t] ~ N[t-1] (NPL)", mean(bareg_nf_cf_npl), std(bareg_nf_cf_npl));
fprintf('%30s %10.4f  %10.4f\n', "Excess Turnover (EPL)", mean(excess_cf_epl), std(excess_cf_epl));
fprintf('%30s %10.4f  %10.4f\n', "Excess Turnover (NPL)", mean(excess_cf_npl), std(excess_cf_npl));
fprintf('%30s %10.4f  %10.4f\n', "Corr Ent Exit (EPL)", mean(corr_ent_exit_cf_epl), std(corr_ent_exit_cf_epl));
fprintf('%30s %10.4f  %10.4f\n', "Corr Ent Exit (NPL)", mean(corr_ent_exit_cf_npl), std(corr_ent_exit_cf_npl));

for i = 1:nplayer
    nf_s_cf_epl_i = reshape(nf_s_cf_epl(i,:,:), [numsval, nrepli]);
    nf_s_cf_npl_i = reshape(nf_s_cf_npl(i,:,:), [numsval, nrepli]);
    freq_active_s_cf_epl_i = reshape(freq_active_s_cf_epl(i,:,:), [numsval, nrepli]);
    freq_active_s_cf_npl_i = reshape(freq_active_s_cf_npl(i,:,:), [numsval, nrepli]);
    fprintf('%30s %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            sprintf("Markets Active %d (EPL)", i), ...
            mean(nf_s_cf_epl_i(1,:)), std(nf_s_cf_epl_i(1,:)), ...
            mean(nf_s_cf_epl_i(2,:)), std(nf_s_cf_epl_i(2,:)), ...
            mean(nf_s_cf_epl_i(3,:)), std(nf_s_cf_epl_i(3,:)), ...
            mean(nf_s_cf_epl_i(4,:)), std(nf_s_cf_epl_i(4,:)), ...
            mean(nf_s_cf_epl_i(5,:)), std(nf_s_cf_epl_i(5,:)));
    fprintf('%30s %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            sprintf("Markets Active %d (NPL)", i), ...
            mean(nf_s_cf_npl_i(1,:)), std(nf_s_cf_npl_i(1,:)), ...
            mean(nf_s_cf_npl_i(2,:)), std(nf_s_cf_npl_i(2,:)), ...
            mean(nf_s_cf_npl_i(3,:)), std(nf_s_cf_npl_i(3,:)), ...
            mean(nf_s_cf_npl_i(4,:)), std(nf_s_cf_npl_i(4,:)), ...
            mean(nf_s_cf_npl_i(5,:)), std(nf_s_cf_npl_i(5,:)));
    fprintf('%30s %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            sprintf("Frequency Active %d (EPL)", i), ...
            mean(freq_active_s_cf_epl_i(1,:)), std(freq_active_s_cf_epl_i(1,:)), ...
            mean(freq_active_s_cf_epl_i(2,:)), std(freq_active_s_cf_epl_i(2,:)), ...
            mean(freq_active_s_cf_epl_i(3,:)), std(freq_active_s_cf_epl_i(3,:)), ...
            mean(freq_active_s_cf_epl_i(4,:)), std(freq_active_s_cf_epl_i(4,:)), ...
            mean(freq_active_s_cf_epl_i(5,:),"omitnan"), std(freq_active_s_cf_epl_i(5,:),"omitnan"));
    fprintf('%30s %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
            sprintf("Frequency Active %d (NPL)", i), ...
            mean(freq_active_s_cf_npl_i(1,:)), std(freq_active_s_cf_npl_i(1,:)), ...
            mean(freq_active_s_cf_npl_i(2,:)), std(freq_active_s_cf_npl_i(2,:)), ...
            mean(freq_active_s_cf_npl_i(3,:)), std(freq_active_s_cf_npl_i(3,:)), ...
            mean(freq_active_s_cf_npl_i(4,:)), std(freq_active_s_cf_npl_i(4,:)), ...
            mean(freq_active_s_cf_npl_i(5,:),"omitnan"), std(freq_active_s_cf_npl_i(5,:),"omitnan"));

    fprintf('%30s %10.4f  %10.4f\n', ...
            sprintf("Avg. Profit %d (EPL)", i), ...
            mean(avgv_cf_epl(i,:)), std(avgv_cf_epl(i,:)));
    fprintf('%30s %10.4f  %10.4f\n', ...
            sprintf("Avg. Profit %d (NPL)", i), ...
            mean(avgv_cf_npl(i,:)), std(avgv_cf_npl(i,:)));
end
fprintf('=====================================================================================================================================================\n');

%% Estimates Table

namesb_latex = [
    "\theta_{\text{FC},\text{SC}}",
    "\theta_{\text{FC},\text{CC}}",
    "\theta_{\text{FC},\text{BJ}}",
    "\theta_{\text{RS}}",
    "\theta_{\text{RN}}",
    "\theta_{\text{EC}}" ];

fprintf("\\begin{table}[h]\n");
fprintf("\\centering\n");
fprintf("\\caption{Warehouse Clubs: Parameter Estimates}\n");
fprintf("\\label{tab:app:epl}\n");
fprintf("\\begin{tabular}{lrrr}\n");
fprintf("\\toprule\n");
fprintf("Parameter                    & Estimate    & S.E. & 95\\%% CI \\\\\n");
fprintf("\\midrule\n");
for k = 1:kparam
    fprintf("$%0s$ & %7.3f & (%0.3f) & [%0.3f, %0.3f] \\\\\n", ...
            namesb_latex(k,:), ...
            bmat_cepl(1,k), std(bmat_cepl(2:nrepli,k)), ...
            prctile(bmat_cepl(2:nrepli,k),[2.5]), prctile(bmat_cepl(2:nrepli,k),[97.5]));
end
fprintf("\\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\end{table}\n");

fprintf('%%==============================%%\n');

%% Counterfactual Table 1

fprintf("\\begin{table}[h]\n");
fprintf("\\centering\n");
fprintf("\\caption{Warehouse Clubs: Counterfactual (Aggregate)}\n");
fprintf("\\label{tab:counterfactual_1}\n");
fprintf("\\begin{tabular}{lrrrrr}\n");
fprintf("\\toprule\n");
fprintf("              & & \\multicolumn{2}{c}{Estimated} &  \\multicolumn{2}{c}{Counterfactual} \\\\\n");
fprintf("              & Observed & Mean & S.E. & Mean & S.E. \\\\\n");
fprintf("\\midrule\n");
fprintf("%20s & %d & %7.3f & (%0.3f) & %7.3f & (%0.3f) \\\\\n", "Active Firms", ...
        obs_nf, ...
        mean(mean_nf_est_epl), std(mean_nf_est_epl), ...
        mean(mean_nf_cf_epl), std(mean_nf_cf_epl));
fprintf("%20s & %d & %7.3f & (%0.3f) & %7.3f & (%0.3f) \\\\\n", "Entries", ...
        obs_entries, ...
        mean(mean_entries_est_epl), std(mean_entries_est_epl), ...
        mean(mean_entries_cf_epl), std(mean_entries_cf_epl));
fprintf("%20s & %d & %7.3f & (%0.3f) & %7.3f & (%0.3f) \\\\\n", "Exits", ...
        obs_exits, ...
        mean(mean_exits_est_epl), std(mean_exits_est_epl), ...
        mean(mean_exits_cf_epl), std(mean_exits_cf_epl));
fprintf("\\midrule\n");
fprintf("Markets with  & &      &      &      &      \\\\\n");
for j = 0:3
    fprintf("\\quad %d Firms & %d & %7.3f & (%0.3f) & %7.3f & (%0.3f) \\\\\n", j, ...
            obs_nf_hist(j+1), ...
            mean(nf_hist_est_epl(j+1,:)), std(nf_hist_est_epl(j+1,:)), ...
            mean(nf_hist_cf_epl(j+1,:)), std(nf_hist_cf_epl(j+1,:)));
end
fprintf("\\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\end{table}\n");

fprintf('%%==============================%%\n');

%% Counterfactual Table 2

fprintf("\\begin{sidewaystable}[p]\n");
fprintf("\\centering\n");
fprintf("\\caption{Wholesale Clubs: Counterfactual (By Market Size)}\n");
fprintf("\\label{tab:counterfactual_2}\n");
fprintf("\\begin{tabular}{@{}lrrrrrrrrrr@{}}\n");
fprintf("\\toprule\n");
fprintf("& \\multicolumn{2}{c}{s=1} & \\multicolumn{2}{c}{s=2} & \\multicolumn{2}{c}{s=3} & \\multicolumn{2}{c}{s=4} & \\multicolumn{2}{c}{s=5} \\\\\n");
fprintf("\\cmidrule(lr){2-3} \\cmidrule(lr){4-5} \\cmidrule(lr){6-7} \\cmidrule(lr){8-9} \\cmidrule(lr){10-11}\n");
fprintf("                 & Mean & S.E. & Mean & S.E. & Mean & S.E. & Mean & S.E. & Mean & S.E. \\\\\n");

for j = 1:2
    if j == 1
        type = "Estimated";
        nf_data = nf_s_est_epl;
    else
        type = "Counterfactual";
        nf_data = nf_s_cf_epl;
    end

    fprintf("\\midrule\n");
    fprintf("\\multicolumn{11}{l}{%s} \\\\\n", type);

    for i = 1:nplayer
        if i == 1
            firm = "Sam's Club";
        elseif i == 2
            firm = "Costco";
        else
            firm = "BJ's";
        end

        nf_data_i = reshape(nf_data(i,:,:), [numsval, nrepli]);

        fprintf("\\quad %20s & %7.3f & (%0.3f) & %7.3f & (%0.3f) & %7.3f & (%0.3f) & %7.3f & (%0.3f) & %7.3f & (%0.3f) \\\\\n", firm, ...
                mean(nf_data_i(1,:)), std(nf_data_i(1,:)), ...
                mean(nf_data_i(2,:)), std(nf_data_i(2,:)), ...
                mean(nf_data_i(3,:)), std(nf_data_i(3,:)), ...
                mean(nf_data_i(4,:)), std(nf_data_i(4,:)), ...
                mean(nf_data_i(5,:)), std(nf_data_i(5,:)));
    end
end
fprintf("\\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\end{sidewaystable}\n");

fprintf('%%==============================%%\n');

fprintf("\\begin{table}[h]\n");
fprintf("\\centering\n");
fprintf("\\caption{Warehouse Clubs: Counterfactual (Profits)}\n");
fprintf("\\label{tab:counterfactual_3}\n");
fprintf("\\begin{tabular}{lrrrr}\n");
fprintf("  \\toprule\n");
fprintf("             & \\multicolumn{2}{c}{Simulated} & \\multicolumn{2}{c}{Counterfactual} \\\\\n");
fprintf("             & Mean  & S.E.  & Mean    & S.E.            \\\\\n");
fprintf("  \\midrule\n");
fprintf("  Sam's Club & %7.3f & (%0.3f) & %7.3f & (%0.3f) \\\\\n", ...
        mean(avgv_est_epl(1,:)), std(avgv_est_epl(1,:)), ...
        mean(avgv_cf_epl(1,:)), std(avgv_cf_epl(1,:)));
fprintf("  Costco     & %7.3f & (%0.3f) & %7.3f & (%0.3f) \\\\\n", ...
        mean(avgv_est_epl(2,:)), std(avgv_est_epl(2,:)), ...
        mean(avgv_cf_epl(2,:)), std(avgv_cf_epl(2,:)));
fprintf("  BJ's       & %7.3f & (%0.3f) & %7.3f & (%0.3f) \\\\\n", ...
        mean(avgv_est_epl(3,:)), std(avgv_est_epl(3,:)), ...
        mean(avgv_cf_epl(3,:)), std(avgv_cf_epl(3,:)));
fprintf("  \\bottomrule\n");
fprintf("\\end{tabular}\n");
fprintf("\\end{table}\n");
