% --------------------------------------------------------------------------
% MPESTAT        Reports summary statistics from data
%
% Originally by Victor Aguirregabiria.
% Converted from Gauss to Matlab by Jason Blevins.
% Updated to include statistics about market size.
% --------------------------------------------------------------------------

function [ mean_nf, mean_entries, mean_exits, nf_s, freq_active_s, std_nf, bareg_nf, excess, corr_ent_exit, hist, avgv ] =  mpestat(aobs, aobs_1, sobs, verbose, sval, nm, vequil, pequil)

  eulerc = 0.5772;
  nobs = size(aobs, 1); % total observations
  nums = size(sval, 1);
  nt = nobs / nm; % time periods per market

  if verbose > 0
      fprintf('\n');
      fprintf('-----------------------------------------------------------------------------------------\n');
      fprintf('   Descriptive Statistics (%d obs.)\n', nobs);
      fprintf('-----------------------------------------------------------------------------------------\n');
      fprintf('\n');
  end
  nplayer = size(aobs, 2);
  nf = sum(aobs')';      % Number of active firms in the market at t
  nf_1 = sum(aobs_1')';  % Number of active firms in the market at t-1

  % Regression of (number of firms t) on (number of firms t-1)
  [ b, sigma, resid ] = mvregress([ ones(nobs, 1), nf_1 ], nf);
  bareg_nf = b(2);  % Estimate of autorregressive parameter
  entries = sum((aobs.*(1-aobs_1))')';   % Number of new entrants at t
  exits = sum(((1-aobs).*aobs_1)')';     % Number of firm exits at t
  excess = mean(entries+exits-abs(entries-exits))'; % Excess turnover
  buff = corr([ entries, exits ]);
  corr_ent_exit = buff(1,2); % Correlation entries and exits
  freq_active = mean(aobs); % Unconditional frequencies of being active
  % Summary statistics
  mean_nf = mean(nf);
  std_nf = std(nf);
  mean_entries = mean(entries)';
  mean_exits = mean(exits)';
  % Distribution of market size
  freq_s = zeros(1, nums);
  for i = 1:nums
      freq_s(i) = mean(sobs == i);
  end
  % Number of firms and activity frequency by market size
  nf_s = zeros(nplayer, nums);
  freq_active_s = zeros(nplayer, nums);
  for i = 1:nums
      inds = (sobs == i);
      counts = sum(inds);
      nf_s(:,i) = sum(repmat((sobs == i), 1, nplayer) .* aobs);
      freq_active_s(:,i) = nf_s(:,i) / counts;
  end
  % Number of markets with 0, 1, 2, or 3 firms at end of sample
  hist = zeros(nplayer+1,1);
  ind_final = nt:nt:nm*nt;
  nf_final = nf(ind_final); % final NM observations (final period, all markets)
  for i = 0:nplayer
      hist(i+1) = sum(nf_final == i);
  end

  % Calculate average profits by firm
  avgv = zeros(nplayer, 1);
  if exist('vequil', 'var') & exist('pequil', 'var')
      numa = 2^nplayer;
      numx = nums*numa;

      % ------------------------------------------------
      % Matrix with observed indexes of state variables
      % ------------------------------------------------
      indsobs = (repmat(sobs, 1, nums)==repmat(sval', nobs, 1))*[ 1:nums ]';
      twop = kron(ones(nobs, 1), 2.^(nplayer-[1:nplayer]));
      indobs = sum(aobs_1.*twop, 2);
      indobs = (indsobs-1).*(2^nplayer) + indobs + 1;

      % ------------------------------------------------
      % Average profits for each player
      % ------------------------------------------------
      for i = 1:nplayer
          for t = 1:nobs
              v = vequil(indobs(t), i, (aobs(t)+1));
              p = pequil(indobs(t), i);
              ep = eulerc - ((1-aobs(t)) * log(1-p) + aobs(t) * log(p));
              avgv(i) = avgv(i) + v + ep;
          end
      end
      avgv = avgv / (1000*nt); % Average over time periods and rescale
  end

  if verbose > 0
      fprintf('\n');
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('       (1)    Average number of active firms   = %12.4f\n', mean_nf);
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('       (2)    Std. Dev. number of firms        = %12.4f\n', std_nf);
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('       (3)    Regression N[t] on N[t-1]        = %12.4f\n', bareg_nf);
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('       (4)    Average number of entrants       = %12.4f\n', mean_entries);
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('       (5)    Average number of exits          = %12.4f\n', mean_exits);
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('       (6)    Excess turnover (in # of firms)  = %12.4f\n', excess);
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('       (7)    Correlation entries and exits    = %12.4f\n', corr_ent_exit);
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('       (8)    Probability of being active      =\n');
      disp(freq_active)
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('       (9)    Distribution of market size      =\n');
      disp(freq_s)
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('      (10)    Obs. active by market size       =\n');
      disp(nf_s)
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('      (11)    Active probability by size       =\n');
      disp(freq_active_s)
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('      (12)    Histogram of number of firms     =\n');
      disp(hist')
      fprintf('----------------------------------------------------------------------------------------\n');
      fprintf('\n');

      %% Summary Statistics Table

      fprintf("\\begin{table}[htbp]\n");
      fprintf("\\centering\n");
      fprintf("\\hspace*{-0.7cm}\n");
      fprintf("\\begin{threeparttable}\n");
      fprintf("\\caption{Warehouse Clubs: Summary Statistics}\n");
      fprintf("\\label{tab:summary_staistics}\n");
      fprintf("\\begin{tabular}{lr}\n");
      fprintf("\\toprule\n");
      fprintf("Statistic                             & Value   \\\\\n");
      fprintf("\\midrule\n");
      fprintf("Average active firms                & %7.3f \\\\\n", mean_nf);
      fprintf("S.D. active firms                   & %7.3f \\\\\n", std_nf);
      fprintf("AR(1) for active firms\\tnote{1}     & %7.3f \\\\\n", bareg_nf);
      fprintf("Average entrants                    & %7.3f \\\\\n", mean_entries);
      fprintf("Average exits                       & %7.3f \\\\\n", mean_exits);
      fprintf("Excess turnover\\tnote{2}              & %7.3f \\\\\n", excess);
      fprintf("Correlation between entries and exits & %7.3f \\\\\n", corr_ent_exit);
      fprintf("Probability of being active           & \\\\\n");
      fprintf("\\quad Sam's Club                      & %7.3f \\\\\n", freq_active(1));
      fprintf("\\quad Costco                          & %7.3f \\\\\n", freq_active(2));
      fprintf("\\quad BJ's                            & %7.3f \\\\\n", freq_active(3));
      fprintf("Distribution of market size           &  \\\\\n");
      fprintf("\\quad $s = 1$                         & %7.3f \\\\\n", freq_s(1));
      fprintf("\\quad $s = 2$                         & %7.3f \\\\\n", freq_s(2));
      fprintf("\\quad $s = 3$                         & %7.3f \\\\\n", freq_s(3));
      fprintf("\\quad $s = 4$                         & %7.3f \\\\\n", freq_s(4));
      fprintf("\\quad $s = 5$                         & %7.3f \\\\\n", freq_s(5));
      fprintf("Markets                                & %d \\\\\n", nm);
      fprintf("Years                                  & %d \\\\\n", nt);
      fprintf("Observations (Markets $\\times$ Years)  & %d \\\\\n", nobs);
      fprintf("\\bottomrule\n");
      fprintf("\\end{tabular}\n");
      fprintf("    \\begin{tablenotes}\n");
      fprintf("\\footnotesize\n");
      fprintf("      \\item[1] AR(1) for \\#active is the autoregressive coefficient regressing the number of current active firms on the number of active firms in previous period. \n");
      fprintf("      \\item[2] Excess turnover is defined as (\\#entrants + \\#exits) - abs(\\#entrants - \\#exits). \n");
      fprintf("    \\end{tablenotes}\n");
      fprintf("\\end{threeparttable}\n");
      fprintf("\\end{table}\n");
  end
end
