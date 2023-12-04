% mc_epl.m -- NPL and EPL Monte Carlo
%
% This is based on a Matlab translation of the Gauss program of
% Aguirregabiria and Mira (2007) originally named
% am_econometrica_2007_montecarlo.prg, with the EPL estimator added.
% Comments from the original are preserved below.
%
% Adam Dearing and Jason Blevins
% Columbus, Ohio
% June 30, 2020

% ****************************************************************************
%
%  am_econometrica_2007_montecarlo.prg
%
%      GAUSS program that implements the Monte Carlo experiments reported in
%      Section 4 of the paper: Aguirregabiria, V., and P. Mira (2007):
%      'Sequential Estimation of Dynamic Discrete Games,' Econometrica,
%      Vol. 75, No. 1, pages 1-53.
%
%  by Victor Aguirregabiria
%
%  - Last version: May 2005
%  - Comments and explanations added on August 2011
%
% ****************************************************************************


% ------------------------------------------------------------------------------
%                              MODEL
% ------------------------------------------------------------------------------
%
%  MAIN FEATURES OF THE MODEL
%
%  - Dynamic game of firm entry-exit in a market
%
%  -   i = Firm index, i belongs to (1,2,3,4,5)
%      t = Time index
%
%  -   Decision variable:
%      a[it] = Indicator of the event 'firm i operates in the market at period t'
%
%          a[it] = 0 ----> No active in market
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
%
% ------------------------------------------------------------------------------
%              MONTE CARLO EXPERIMENTS
% ------------------------------------------------------------------------------
%
%  - We implement 6 experiments.
%
%  - The following parameters are the same in the 6 experiments:
%
%      Number of local markets (M) =   400
%      Number of time periods  (T) =   1
%      Number of players (N)       =   5
%      Number of Monte Carlo simulations   =   1000
%
%      Discount factor     =   0.95
%      theta_fc_1          =   -1.9
%      theta_fc_2          =   -1.8
%      theta_fc_3          =   -1.7
%      theta_fc_4          =   -1.6
%      theta_fc_5          =   -1.5
%      theta_rs            =   1.0
%
%      Values for s[t]     =   (1,2,3,4,5)
%
%      Transition probability of s[t]  =   ( 0.8 ~ 0.2 ~ 0.0 ~ 0.0 ~ 0.0 )
%                                        | ( 0.2 ~ 0.6 ~ 0.2 ~ 0.0 ~ 0.0 )
%                                        | ( 0.0 ~ 0.2 ~ 0.6 ~ 0.2 ~ 0.0 )
%                                        | ( 0.0 ~ 0.0 ~ 0.2 ~ 0.6 ~ 0.2 )
%                                        | ( 0.0 ~ 0.0 ~ 0.0 ~ 0.2 ~ 0.8 )
%
%  - The following values of the parameters theta_rn and theta_ec
%    defines experiments 1 to 6:
%
%          Experiment 1:   theta_ec = 1.0      and     theta_rn = 1.0
%          Experiment 2:   theta_ec = 1.0      and     theta_rn = 3.0
%          Experiment 3:   theta_ec = 1.0      and     theta_rn = 6.0

% -----------------------------------------------------------
%  SELECTION OF THE MONTE CARLO EXPERIMENT TO IMPLEMENT
% -----------------------------------------------------------

if ~exist('selexper')
    selexper = 1; % Value from 1 to 6 that represents the index of the
                  % Monte Carlo experiment to implement. A run of this
                  % program implements one experiment.
end

if ~exist('nobs')
    nobs = 1600;  % Number of markets (observations)
end

% -----------------------------------------------
%  VALUES OF PARAMETERS AND OTHER CONSTANTS
% -----------------------------------------------
nrepli = 1000;  % Number of Monte carlo replications
nplayer = 5;    %  Number of players
naction = 2;    % Number of actions
numexp = 3;     %  Total number of Monte Carlo experiments
nstart = 1;     % Number of starting values
miniter = 3;    %  Minimum number of NPL or EPL iterations
maxiter = 100;  %  Maximum number of NPL or EPL iterations

theta_fc = zeros(numexp, nplayer);
theta_fc(:, 1) = -1.9 * ones(numexp, 1); % Vector with values of theta_fc_1 for each experiment
theta_fc(:, 2) = -1.8 * ones(numexp, 1); % Vector with values of theta_fc_2 for each experiment
theta_fc(:, 3) = -1.7 * ones(numexp, 1); % Vector with values of theta_fc_3 for each experiment
theta_fc(:, 4) = -1.6 * ones(numexp, 1); % Vector with values of theta_fc_4 for each experiment
theta_fc(:, 5) = -1.5 * ones(numexp, 1); % Vector with values of theta_fc_1 for each experiment

theta_rs = 1.0 * ones(numexp, 1); % Vector with values of theta_rs for each experiment
disfact = 0.95 * ones(numexp, 1); % Vector with values of discount factor for each experiment
sigmaeps = 1 * ones(numexp, 1);   % Vector with values of std. dev. epsilon for each experiment

theta_rn = [ 1.0; 2.5; 4.0 ]; % Vector with values of theta_rn for each experiment
theta_ec = [ 1.0; 1.0; 1.0 ]; % Vector with values of theta_ec for each experiment

% Points of support and transition probability of state variable s[t], market size
sval = [ 1:1:5 ]';  % Support of market size
numsval = size(sval, 1);  % Number of possible market sizes
nstate = numsval * (2^nplayer);  % Number of points in the state space
ptrans = [ 0.8, 0.2, 0.0, 0.0, 0.0 ;
           0.2, 0.6, 0.2, 0.0, 0.0 ;
           0.0, 0.2, 0.6, 0.2, 0.0 ;
           0.0, 0.0, 0.2, 0.6, 0.2 ;
           0.0, 0.0, 0.0, 0.2, 0.8 ];

% Selecting the parameters for the experiment
theta_fc = theta_fc(selexper, :)';
theta_rs = theta_rs(selexper);
disfact = disfact(selexper);
sigmaeps = sigmaeps(selexper);
theta_ec = theta_ec(selexper);
theta_rn = theta_rn(selexper);

% Vector with true values of parameters
trueparam = [ theta_fc; theta_rs; theta_rn; theta_ec; disfact; sigmaeps ];

% Vector with names of parameters of profit function
namesb = [ 'FC_1'; 'FC_2'; 'FC_3'; 'FC_4'; 'FC_5'; '  RS'; '  RN'; '  EC' ];

% Number of parameters to estimate
kparam = size(trueparam, 1) - 2;

% Structure for storing true parameters and settings
param.theta_fc = theta_fc;
param.theta_rs = theta_rs;
param.theta_rn = theta_rn;
param.theta_ec = theta_ec;
param.disfact = disfact;
param.sigmaeps = sigmaeps;
param.sval = sval;
param.ptrans = ptrans;
param.verbose = 1;

% Seed for (pseudo) random number generation
rand('seed', 20150403);

% -----------------------------------------------------------------------
%  COMPUTING A MARKOV PERFECT EQUILIBRIUM OF THIS DYNAMIC GAME
% -----------------------------------------------------------------------
eqfile = sprintf('mc_epl_eq%d.mat', selexper);
if isfile(eqfile)
    disp('Loading pre-calculated equilibrium from file.')
    load(eqfile);
else
    disp('Solving for an equilibrium of the game.')
    v0 = zeros(nstate*nplayer*naction,1);
    eqcond = @(v) Gfunc(sval, ptrans, trueparam(1:kparam), v, disfact, nplayer, naction);
    vequil = fsolve(eqcond, v0, optimoptions('fsolve', 'Display', 'iter'));
    vmat = reshape(vequil, [nstate, nplayer, naction]);
    pequil = exp(vmat(:,:,2))./(exp(vmat(:,:,1))+exp(vmat(:,:,2)));
    save(eqfile, 'vequil', 'vmat', 'pequil');
end

% -----------------------------------------------------------------------
%  COMPUTING THE STEADY STATE DISTRIBUTION OF THIS DYNAMIC GAME
% -----------------------------------------------------------------------
numa = 2^nplayer;
aval = zeros(numa, nplayer);
for i = 1:nplayer
  aval(:,i) = kron(kron(ones(2^(i-1),1), [ 0; 1 ]), ones(2^(nplayer - i), 1));
end
ptrana = ones(nstate, numa);
for i = 1:nplayer;
  mi = aval(:,i)';
  ppi = pequil(:,i);
  ppi1 = repmat(ppi, 1, numa);
  ppi0 = 1 - ppi1;
  mi1 = repmat(mi, nstate, 1);
  mi0 = 1 - mi1;
  ptrana = ptrana .* (ppi1 .^ mi1) .* (ppi0 .^ mi0);
end
ptrana = (kron(param.ptrans, ones(numa, numa))) .* (kron(ones(1,numsval), ptrana));
critconv = (1e-3)*(1/nstate);
criter = 1000;
psteady0 = (1/nstate) * ones(nstate, 1);
while (criter > critconv)
  psteady1 = ptrana' * psteady0;
  criter = max(abs(psteady1 - psteady0), [], 1);
  psteady0 = psteady1;
end
psteady = psteady0;

vstate = zeros(nstate, nplayer + 1);
vstate(:, 1) = kron(param.sval, ones(numa, 1));
vstate(:, 2:nplayer+1) = kron(ones(numsval, 1), aval);

% --------------------------------------------------------------------------
%  SIMULATING DATA (50,000 OBSERVATIONS) FROM THE EQUILIBRIUM
%  TO OBTAIN DESCRIPTIVE STATISTICS ON THE DYNAMICS OF MARKET STRUCTURE
% --------------------------------------------------------------------------
nobsfordes = 50000;
[ aobs, aobs_1, sobs ] = simdygam(nobsfordes, pequil, psteady, vstate);
mpestat(aobs, aobs_1);

% --------------------------------------------------------------------------
%  MONTE CARLO EXPERIMENT
% --------------------------------------------------------------------------

fprintf('\n');
fprintf('*****************************************************************************************\n');
fprintf('       MONTE CARLO EXPERIMENT #%d\n', selexper);
fprintf('*****************************************************************************************\n');
fprintf('\n');
param.verbose = 1;

% Cell arrays of starting values for each iteration
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
start_times = {}; % Sequence of computational times per estimate

% Matrices of estimates for each iteration
bmat_1npl = zeros(nrepli, kparam); % 1-NPL estimates
bmat_2npl = zeros(nrepli, kparam); % 2-NPL estimates
bmat_3npl = zeros(nrepli, kparam); % 3-NPL estimates
bmat_cnpl = zeros(nrepli, kparam); % Converged NPL estimates
bmat_1epl = zeros(nrepli, kparam); % 1-EPL estimates
bmat_2epl = zeros(nrepli, kparam); % 2-EPL estimates
bmat_3epl = zeros(nrepli, kparam); % 3-EPL estimates
bmat_cepl = zeros(nrepli, kparam); % Converged EPL estimates
btrue_npl = zeros(nrepli, kparam); % 1-NPL estimates using true CCPs
btrue_epl = zeros(nrepli, kparam); % 1-EPL estimates using true value function and parameters

% Iteration counts for each method
iter_cnpl = zeros(nrepli, 1);
iter_cepl = zeros(nrepli, 1);

% Total computational times for each method
time_1npl = zeros(nrepli, 1);
time_1epl = zeros(nrepli, 1);
time_2npl = zeros(nrepli, 1);
time_2epl = zeros(nrepli, 1);
time_3npl = zeros(nrepli, 1);
time_3epl = zeros(nrepli, 1);
time_cnpl = zeros(nrepli, 1);
time_cepl = zeros(nrepli, 1);

% Per iteration times for each method
time_2npl_iter = zeros(nrepli, 1);
time_2epl_iter = zeros(nrepli, 1);
time_3npl_iter = zeros(nrepli, 1);
time_3epl_iter = zeros(nrepli, 1);
time_cnpl_iter = zeros(nrepli, 1);
time_cepl_iter = zeros(nrepli, 1);

% Failure counts for each method
fail_cnpl = zeros(nrepli, 1);
fail_cepl = zeros(nrepli, 1);

% Note: When there are multicollinearity problems in a Monte Carlo
% sample we ignore that sample and take a new one. We want to check
% for the number of times we have to make these redraws.
redraws = 0;

for draw = 1:nrepli
  fprintf('=========================================================================================\n');
  fprintf('     Replication = %d\n', draw);
  fprintf('-----------------------------------------------------------------------------------------\n');
  fprintf('         (a.0)   Simulations of x''s and a''s\n');
  flag = 1;
  while (flag==1)
    [ aobs, aobs_1, sobs] = simdygam(nobs, pequil, psteady, vstate);
    check = sum(sum([ aobs, aobs_1 ]) == zeros(1, 2*nplayer));
    check = check + sum(sum([ aobs, aobs_1 ]) == (nobs .* ones(1, 2*nplayer)));
    if (check > 0)
        flag = 1;
    elseif (check == 0)
        flag = 0;
    end
    redraws = redraws + flag; % Counts the number re-drawings
  end

  fprintf('-----------------------------------------------------------------------------------------\n');
  fprintf('         (a.1)   Estimation of initial CCPs (Semi-Parametric: Logit)\n');
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

  fprintf('-----------------------------------------------------------------------------------------\n');
  fprintf('         (a.2)   Perturbations of Semi-Parametric Logit CCPs (10)\n');

  for r = 2:min(3, nstart)
      logit_perturb = best_logit + randn(kparam, 1);
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

  fprintf('-----------------------------------------------------------------------------------------\n');
  fprintf('         (a.3)   Estimation of initial CCPs (Non-Parametric)\n');
  prob0NP = freqprob(aobs, [ sobs, aobs_1 ], vstate);

  start_name{4} = 'Non-Parametric';
  start_name_short{4} = 'NP';
  start_p{4} = prob0NP;
  start_1npl(4) = 1;

  fprintf('-----------------------------------------------------------------------------------------\n');
  fprintf('         (a.4)   Completely random starting values...\n');

  for r = 5:nstart
      random_theta = randn(kparam, 1);
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

  %--------------------------------------------------
  % NPL
  %--------------------------------------------------

  for r = 1:nstart
      fprintf('-----------------------------------------------------------------------------------------\n');
      fprintf('         (b.%d)   NPL Algorithm: %s\n', r, start_name{r});
      [ theta_seq{r}, varb, pnpl, start_v1npl{r}, vcnpl, start_ll(r), start_iter(r), start_err(r), start_times{r} ] = npldygam(aobs, sobs, aobs_1, sval, ptrans, start_p{r}, disfact, miniter, maxiter);
  end

  max_llike = max(start_ll);
  best_r = find(start_ll == max_llike, 1, 'first');
  if best_r && sum(start_err) < nstart
      bmat_1npl(draw,:) = theta_seq{best_r}(1,:);
      bmat_2npl(draw,:) = theta_seq{best_r}(2,:);
      bmat_3npl(draw,:) = theta_seq{best_r}(3,:);
      bmat_cnpl(draw,:) = theta_seq{best_r}(start_iter(best_r),:);
      time_1npl(draw,1) = start_times{best_r}(1);
      time_2npl(draw,1) = sum(start_times{best_r}(1:2));
      time_3npl(draw,1) = sum(start_times{best_r}(1:3));
      time_cnpl(draw,1) = sum(start_times{best_r});
      iter_cnpl(draw,1) = start_iter(best_r);
      time_2npl_iter(draw,1) = time_2npl(draw,1) / 2;
      time_3npl_iter(draw,1) = time_3npl(draw,1) / 3;
      time_cnpl_iter(draw,1) = time_cnpl(draw,1) / iter_cnpl(draw,1);
      fail_cnpl(draw,1) = 0;
  else
      fail_cnpl(draw,1) = 1;
  end

  fprintf('-----------------------------------------------------\n');
  fprintf('NPL Summary (Replication %d)\n', draw);
  fprintf('-----------------------------------------------------\n');
  fprintf('%10s %10s', 'Start Val.', 'True');
  for r = 1:nstart
      fprintf(' %10s', start_name_short{r});
  end
  fprintf('\n');
  fprintf('-----------------------------------------------------\n');
  for k = 1:kparam
      fprintf('%10s %10.4f', namesb(k,:), trueparam(k))
      for r = 1:nstart
          fprintf(' %10.4f', theta_seq{r}(start_iter(r),k));
      end
      fprintf('\n');
  end

  fprintf('%10s %10s', 'Iter', '');
  for r = 1:nstart
      fprintf(' %10d', start_iter(r));
  end
  fprintf('\n');

  fprintf('%10s %10s', 'LL', '');
  for r = 1:nstart
      fprintf(' %10.3f', start_ll(r));
  end
  fprintf('\n');

  fprintf('%10s %10s', 'Best', '');
  for r = 1:nstart
      fprintf(' %10d', best_r == r);
  end
  fprintf('\n');

  fprintf('%10s %10s', 'MaxIter', '');
  for r = 1:nstart
      fprintf(' %10d', start_iter(r) == maxiter);
  end
  fprintf('\n');

  fprintf('%10s %10s', 'Error', '');
  for r = 1:nstart
      fprintf(' %10d', start_err(r));
  end
  fprintf('\n');

  %--------------------------------------------------
  % EPL
  %--------------------------------------------------

  for r = 1:nstart
      fprintf('-----------------------------------------------------------------------------------------\n');
      fprintf('         (b.%d)   EPL Algorithm: %s\n', r, start_name{r});
      if (start_1npl(r) > 0)
          v0 = start_v1npl{r};
          theta0 = start_theta1(:,r);
      else
          v0 = start_v{r};
          theta0 = start_theta(:,r);
      end
      [ theta_seq{r}, varb, start_ll(r), start_iter(r), start_err(r), start_times{r} ] = epldygam(aobs, sobs, aobs_1, sval, ptrans, theta0, v0, disfact, miniter, maxiter);
  end

  max_llike = max(start_ll);
  best_r = find(start_ll == max_llike, 1, 'first');
  if best_r && sum(start_err) < nstart
      bmat_1epl(draw,:) = theta_seq{best_r}(1,:);
      bmat_2epl(draw,:) = theta_seq{best_r}(2,:);
      bmat_3epl(draw,:) = theta_seq{best_r}(3,:);
      bmat_cepl(draw,:) = theta_seq{best_r}(start_iter(best_r),:);
      time_1epl(draw,1) = start_times{best_r}(1);
      time_2epl(draw,1) = sum(start_times{best_r}(1:2));
      time_3epl(draw,1) = sum(start_times{best_r}(1:3));
      time_cepl(draw,1) = sum(start_times{best_r});
      iter_cepl(draw,1) = start_iter(best_r);
      time_2epl_iter(draw,1) = time_2epl(draw,1) / 2;
      time_3epl_iter(draw,1) = time_3epl(draw,1) / 3;
      time_cepl_iter(draw,1) = time_cepl(draw,1) / iter_cepl(draw,1);
      fail_cepl(draw,1) = 0;
  else
      fail_cepl(draw,1) = 1;
  end

  fprintf('-----------------------------------------------------\n');
  fprintf('EPL Summary (Replication %d)\n', draw);
  fprintf('-----------------------------------------------------\n');
  fprintf('%10s %10s', 'Start Val.', 'True');
  for r = 1:nstart
      fprintf(' %10s', start_name_short{r});
  end
  fprintf('\n');
  fprintf('-----------------------------------------------------\n');
  for k = 1:kparam
      fprintf('%10s %10.4f', namesb(k,:), trueparam(k))
      for r = 1:nstart
          fprintf(' %10.4f', theta_seq{r}(start_iter(r),k));
      end
      fprintf('\n');
  end

  fprintf('%10s %10s', 'Iter', '');
  for r = 1:nstart
      fprintf(' %10d', start_iter(r));
  end
  fprintf('\n');

  fprintf('%10s %10s', 'LL', '');
  for r = 1:nstart
      fprintf(' %10.3f', start_ll(r));
  end
  fprintf('\n');

  fprintf('%10s %10s', 'Best', '');
  for r = 1:nstart
      fprintf(' %10d', best_r == r);
  end
  fprintf('\n');

  fprintf('%10s %10s', 'MaxIter', '');
  for r = 1:nstart
      fprintf(' %10d', start_iter(r) == maxiter);
  end
  fprintf('\n');

  fprintf('%10s %10s', 'Error', '');
  for r = 1:nstart
      fprintf(' %10d', start_err(r));
  end
  fprintf('\n');

  fprintf('-----------------------------------------------------------------------------------------\n');
  fprintf('         (d.1)   NPL algorithm using true values as initial CCPs\n');
  [ b1npl, varb, pnpl, v1npl, vcnpl, iter ] = npldygam(aobs, sobs, aobs_1, sval, ptrans, pequil, disfact, 1, 1);
  btrue_npl(draw,:) = b1npl;

  fprintf('-----------------------------------------------------------------------------------------\n');
  fprintf('         (d.2)   EPL algorithm using theta from 1-NPL w/true CCPs\n');
  [ b1epl, varb, iter, err ] = epldygam(aobs, sobs, aobs_1, sval, ptrans, trueparam(1:kparam,1), vmat, disfact, 1, 1);
  btrue_epl(draw,:) = b1epl;
end

fprintf('--------------------------------------------------------------------------------------------------------------\n');
fprintf('           Number of Re-drawings due to Multicollinearity = %d\n',  redraws);

% True Parameters
theta_true = trueparam(1:kparam,1)';

% Replication selection
keep_cnpl = fail_cnpl == 0;
keep_cepl = fail_cepl == 0;

% Empirical Means
mean_1npl = mean(bmat_1npl);
mean_2npl = mean(bmat_2npl);
mean_3npl = mean(bmat_3npl);
mean_cnpl = mean(bmat_cnpl(keep_cnpl,:));
mean_1epl = mean(bmat_1epl);
mean_2epl = mean(bmat_2epl);
mean_3epl = mean(bmat_3epl);
mean_cepl = mean(bmat_cepl(keep_cepl,:));
mean_1npl_true = mean(btrue_npl);
mean_1epl_true = mean(btrue_epl);

% Mean Bias
meanbias_1npl = mean_1npl - theta_true;
meanbias_2npl = mean_2npl - theta_true;
meanbias_3npl = mean_3npl - theta_true;
meanbias_cnpl = mean_cnpl - theta_true;
meanbias_1epl = mean_1epl - theta_true;
meanbias_2epl = mean_2epl - theta_true;
meanbias_3epl = mean_3epl - theta_true;
meanbias_cepl = mean_cepl - theta_true;
meanbias_1npl_true = mean_1npl_true - theta_true;
meanbias_1epl_true = mean_1epl_true - theta_true;

% Empirical Variances
var_1npl = var(bmat_1npl);
var_2npl = var(bmat_2npl);
var_3npl = var(bmat_3npl);
var_cnpl = var(bmat_cnpl(keep_cnpl,:));
var_1epl = var(bmat_1epl);
var_2epl = var(bmat_2epl);
var_3epl = var(bmat_3epl);
var_cepl = var(bmat_cepl(keep_cepl,:));
var_1npl_true = var(btrue_npl);
var_1epl_true = var(btrue_epl);

% MSE
mse_1npl = meanbias_1npl.^2 + var_1npl;
mse_2npl = meanbias_2npl.^2 + var_2npl;
mse_3npl = meanbias_3npl.^2 + var_3npl;
mse_cnpl = meanbias_cnpl.^2 + var_cnpl;
mse_1epl = meanbias_1epl.^2 + var_1epl;
mse_2epl = meanbias_2epl.^2 + var_2epl;
mse_3epl = meanbias_3epl.^2 + var_3epl;
mse_cepl = meanbias_cepl.^2 + var_cepl;
mse_1npl_true = meanbias_1npl_true.^2 + var_1npl_true;
mse_1epl_true = meanbias_1epl_true.^2 + var_1epl_true;

% Iterations
median_iter_cnpl = median(iter_cnpl(keep_cnpl));
median_iter_cepl = median(iter_cepl(keep_cepl));
max_iter_cnpl = max(iter_cnpl(keep_cnpl));
max_iter_cepl = max(iter_cepl(keep_cepl));
iqr_iter_cnpl = iqr(iter_cnpl(keep_cnpl));
iqr_iter_cepl = iqr(iter_cepl(keep_cepl));
iter_fail_cnpl = sum(iter_cnpl(keep_cnpl) == maxiter);
iter_fail_cepl = sum(iter_cepl(keep_cepl) == maxiter);

% Total time and time per iteration
median_time_1npl = median(time_1npl(keep_cnpl));
median_time_1epl = median(time_1epl(keep_cepl));
median_time_2npl = median(time_2npl(keep_cnpl));
median_time_2epl = median(time_2epl(keep_cepl));
median_time_3npl = median(time_3npl(keep_cnpl));
median_time_3epl = median(time_3epl(keep_cepl));
median_time_cnpl = median(time_cnpl(keep_cnpl));
median_time_cepl = median(time_cepl(keep_cepl));

median_time_1npl_iter = median_time_1npl;
median_time_1epl_iter = median_time_1epl;
median_time_2npl_iter = median(time_2npl_iter(keep_cnpl));
median_time_2epl_iter = median(time_2epl_iter(keep_cepl));
median_time_3npl_iter = median(time_3npl_iter(keep_cnpl));
median_time_3epl_iter = median(time_3epl_iter(keep_cepl));
median_time_cnpl_iter = median(time_cnpl_iter(keep_cnpl));
median_time_cepl_iter = median(time_cepl_iter(keep_cepl));

mean_time_1npl = mean(time_1npl(keep_cnpl));
mean_time_1epl = mean(time_1epl(keep_cepl));
mean_time_2npl = mean(time_2npl(keep_cnpl));
mean_time_2epl = mean(time_2epl(keep_cepl));
mean_time_3npl = mean(time_3npl(keep_cnpl));
mean_time_3epl = mean(time_3epl(keep_cepl));
mean_time_cnpl = mean(time_cnpl(keep_cnpl));
mean_time_cepl = mean(time_cepl(keep_cepl));

mean_time_1npl_iter = mean_time_1npl;
mean_time_1epl_iter = mean_time_1epl;
mean_time_2npl_iter = mean(time_2npl_iter(keep_cnpl));
mean_time_2epl_iter = mean(time_2epl_iter(keep_cepl));
mean_time_3npl_iter = mean(time_3npl_iter(keep_cnpl));
mean_time_3epl_iter = mean(time_3epl_iter(keep_cepl));
mean_time_cnpl_iter = mean(time_cnpl_iter(keep_cnpl));
mean_time_cepl_iter = mean(time_cepl_iter(keep_cepl));

total_time_1npl = sum(time_1npl(keep_cnpl));
total_time_1epl = sum(time_1epl(keep_cepl));
total_time_2npl = sum(time_2npl(keep_cnpl));
total_time_2epl = sum(time_2epl(keep_cepl));
total_time_3npl = sum(time_3npl(keep_cnpl));
total_time_3epl = sum(time_3epl(keep_cepl));
total_time_cnpl = sum(time_cnpl(keep_cnpl));
total_time_cepl = sum(time_cepl(keep_cepl));

% Report
fprintf('=============================================================================================================\n');
fprintf('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', 'Bias', 'True', '1-NPL', '1-EPL', '2-NPL', '2-EPL', '3-NPL', '3-EPL', 'C-NPL', 'C-EPL');
fprintf('-------------------------------------------------------------------------------------------------------------\n');
for k = 1:kparam
    fprintf('%10s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n', ...
            namesb(k,:), theta_true(k), ...
            meanbias_1npl(k), meanbias_1epl(k), meanbias_2npl(k), meanbias_2epl(k), ...
            meanbias_3npl(k), meanbias_3epl(k), meanbias_cnpl(k), meanbias_cepl(k));
end
fprintf('=============================================================================================================\n');
fprintf('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', 'MSE', 'True', '1-NPL', '1-EPL', '2-NPL', '2-EPL', '3-NPL', '3-EPL', 'C-NPL', 'C-EPL');
fprintf('-------------------------------------------------------------------------------------------------------------\n');
for k = 1:kparam
    fprintf('%10s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n', ...
            namesb(k,:), theta_true(k), ...
            mse_1npl(k), mse_1epl(k), mse_2npl(k), mse_2epl(k), ...
            mse_3npl(k), mse_3epl(k), mse_cnpl(k), mse_cepl(k));
end
fprintf('=============================================================================================================\n');
fprintf('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', 'Iterations', 'True', '1-NPL', '1-EPL', '2-NPL', '2-EPL', '3-NPL', '3-EPL', 'C-NPL', 'C-EPL');
fprintf('-------------------------------------------------------------------------------------------------------------\n');
fprintf('%10s %10s %10d %10d %10d %10d %10d %10d %10.1f %10.1f\n', ...
            'Median', '', 1, 1, 2, 2, 3, 3, median_iter_cnpl, median_iter_cepl);
fprintf('%10s %10s %10d %10d %10d %10d %10d %10d %10d %10d\n', ...
            'Max', '', 1, 1, 2, 2, 3, 3, max_iter_cnpl, max_iter_cepl);
fprintf('%10s %10s %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f %10.1f\n', ...
            'IQR', '', 0, 0, 0, 0, 0, 0, iqr_iter_cnpl, iqr_iter_cepl);
fprintf('%10s %10s %10s %10s %10s %10s %10s %10s %10d %10d\n', ...
            'Non-Conv.', '', '', '', '', '', '', '', iter_fail_cnpl, iter_fail_cepl);
fprintf('=============================================================================================================\n');
fprintf('%10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n', 'Time (s)', 'True', '1-NPL', '1-EPL', '2-NPL', '2-EPL', '3-NPL', '3-EPL', 'C-NPL', 'C-EPL');
fprintf('-------------------------------------------------------------------------------------------------------------\n');
fprintf('%10s %10s %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n', ...
            'Total', '', total_time_1npl, total_time_1epl, total_time_2npl, total_time_2epl, ...
                         total_time_3npl, total_time_3epl, total_time_cnpl, total_time_cepl);
fprintf('%10s %10s %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n', ...
            'Mean', '', mean_time_1npl, mean_time_1epl, mean_time_2npl, mean_time_2epl, ...
                        mean_time_3npl, mean_time_3epl, mean_time_cnpl, mean_time_cepl);
fprintf('%10s %10s %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f\n', ...
            'Median', '', median_time_1npl, median_time_1epl, median_time_2npl, median_time_2epl, ...
                          median_time_3npl, median_time_3epl, median_time_cnpl, median_time_cepl);
fprintf('%10s %10s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n', ...
            'Med/Iter', '', median_time_1npl_iter, median_time_1epl_iter, median_time_2npl_iter, median_time_2epl_iter, ...
                            median_time_3npl_iter, median_time_3epl_iter, median_time_cnpl_iter, median_time_cepl_iter);
fprintf('=============================================================================================================\n');


fprintf('==================================================================================================\n');
fprintf('Tests for Normality\n');
fprintf('%10s %43s %43s\n', '', 'C-NPL', 'C-EPL');
fprintf('%10s %10s %10s %10s %10s %10s %10s %10s %10s\n', 'Param', 'KS', 'CV', 'Reject', 'p', 'KS', 'CV', 'Reject', 'p');
fprintf('--------------------------------------------------------------------------------------------------\n');
for k = 1:kparam
    mean_npl_k = mean(bmat_cnpl(keep_cnpl, k));
    std_npl_k = std(bmat_cnpl(keep_cnpl, k));
    npl_k_std = (bmat_cnpl(keep_cnpl, k) - mean_npl_k) / std_npl_k;
    [h_npl, p_npl, k_npl, cv_npl] = kstest(npl_k_std);

    mean_epl_k = mean(bmat_cepl(keep_cepl, k));
    std_epl_k = std(bmat_cepl(keep_cepl, k));
    epl_k_std = (bmat_cepl(keep_cepl, k) - mean_epl_k) / std_epl_k;
    [h_epl, p_epl, k_epl, cv_epl] = kstest(epl_k_std);

    fprintf('%10s %10.3f %10.3f %10i %10.3f %10.3f %10.3f %10i %10.3f\n', namesb(k,:), k_npl, cv_npl, h_npl, p_npl, k_epl, cv_epl, h_epl, p_epl);
end
fprintf('==================================================================================================\n');
