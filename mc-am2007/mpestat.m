% --------------------------------------------------------------------------
% MPESTAT        Reports summary statistics from simulated dynamic game data
%
% by Victor Aguirregabiria
%
% Converted from Gauss to Matlab by Jason Blevins.
%
% --------------------------------------------------------------------------
%
% FORMAT:
%         mpestat(aobs, aobs_1)
%
% INPUTS:
%
%     aobs    - (nobs x nplayer) matrix with observations of
%               firms' activity decisions (1=active ; 0=no active)
%
%     aobs_1  - (nobs x nplayer) matrix with observations of
%               firms' initial states (1=incumbent; 0=potential entrant)
%
% --------------------------------------------------------------------------

function mpestat(aobs, aobs_1)

  global sval theta_fc theta_rs theta_rn theta_ec disfact sigmaeps ptrans;

  nobsfordes = size(aobs, 1);

  fprintf('\n');
  fprintf('*****************************************************************************************\n');
  fprintf('   DESCRIPTIVE STATISTICS FROM THE EQUILIBRIUM\n');
  fprintf('   BASED ON %d OBSERVATIONS\n', nobsfordes);
  fprintf('\n');
  fprintf('   TABLE 2 OF THE PAPER AGUIRREGABIRIA AND MIRA (2007)\n');
  fprintf('*****************************************************************************************\n');
  fprintf('\n');
  nf = sum(aobs')';      % Number of active firms in the market at t
  nf_1 = sum(aobs_1')';  % Number of active firms in the market at t-1

  %  Regression of (number of firms t) on (number of firms t-1)
  [ b, sigma, resid ] = mvregress([ ones(nobsfordes, 1), nf_1 ], nf);
  bareg_nf = b(2);  % Estimate of autorregressive parameter
  entries = sum((aobs.*(1-aobs_1))')';   % Number of new entrants at t
  exits = sum(((1-aobs).*aobs_1)')';     % Number of firm exits at t
  excess = mean(entries+exits-abs(entries-exits))'; % Excess turnover
  buff = corr([ entries, exits ]);
  corr_ent_exit = buff(1,2); % Correlation entries and exits
  freq_active = mean(aobs)'; % Frequencies of being active

  fprintf('\n');
  fprintf('----------------------------------------------------------------------------------------\n');
  fprintf('       (1)    Average number of active firms   = %12.4f\n', mean(nf));
  fprintf('----------------------------------------------------------------------------------------\n');
  fprintf('       (2)    Std. Dev. number of firms        = %12.4f\n', std(nf));
  fprintf('----------------------------------------------------------------------------------------\n');
  fprintf('       (3)    Regression N[t] on N[t-1]        = %12.4f\n', bareg_nf);
  fprintf('----------------------------------------------------------------------------------------\n');
  fprintf('       (4)    Average number of entrants       = %12.4f\n', mean(entries)');
  fprintf('----------------------------------------------------------------------------------------\n');
  fprintf('       (5)    Average number of exits          = %12.4f\n', mean(exits)');
  fprintf('----------------------------------------------------------------------------------------\n');
  fprintf('       (6)    Excess turnover (in # of firms)  = %12.4f\n', excess);
  fprintf('----------------------------------------------------------------------------------------\n');
  fprintf('       (7)    Correlation entries and exits    = %12.4f\n', corr_ent_exit);
  fprintf('----------------------------------------------------------------------------------------\n');
  fprintf('       (8)    Probability of being active      =\n');
  disp(freq_active)
  fprintf('----------------------------------------------------------------------------------------\n');
  fprintf('\n');

end
