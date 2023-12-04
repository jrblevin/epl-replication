% -----------------------------------------------------------------------
%
% FORWARDSIM      Simulates data on state transition and decisions.
%
% --------------------------------------------------------------------------
%
% FORMAT:
% [ aobs, aobs_1, sobs ] = forwardsim(aobs_1, sobs, pchoice, ptrans, mstate, sval, nobs)
%
% INPUTS
%     aobs_1 - (nobs x nplayer) matrix with players' initial states.
%
%     sobs    - (nobs x 1) vector of initial market states
%
%     pchoice - (nstate x nplayer) matrix with MPE probs of entry
%
%     ptrans  - (nums x nums) matrix of transition probabilities
%               of market characteristics.
%
%     mstate  - (nstate x (nplayer+1)) matrix with values of state variables {s[t],a[t-1]}.
%              The states in the rows of prob are ordered as follows.
%              Example: sval=(1|2) and 3 players:
%                       s[t]    a[1,t-1]    a[2,t-1]    a[3,t-1]
%              Row 1:     1           0           0           0
%              Row 2:     1           0           0           1
%              Row 3:     1           0           1           0
%              Row 4:     1           0           1           1
%              Row 5:     1           1           0           0
%              Row 6:     1           1           0           1
%              Row 7:     1           1           1           0
%              Row 8:     1           1           1           1
%              Row 9:     2           0           0           0
%              Row 10:    2           0           0           1
%              Row 11:    2           0           1           0
%              Row 12:    2           0           1           1
%              Row 13:    2           1           0           0
%              Row 14:    2           1           0           1
%              Row 15:    2           1           1           0
%              Row 16:    2           1           1           1
%
%     sval    - (nums x 1) vector with values of market characteristics
%
%     nobs    - Number of simulations (markets)
%
% OUTPUTS:
%     aobs   - (nobs x nplayer) matrix with players' choices.
%
%     aobs_1 - (nobs x nplayer) matrix with players' initial states.
%
%     sobs   - (nobs x 1) vector with simulated values of s[t]
%
% ----------------------------------------------------------------------------

function [ aobs, sobs ] = forwardsim(aobs_1, sobs, pchoice, ptrans, sval, nobs)

  nplayer = size(pchoice, 2);
  numx = size(pchoice, 1);
  numa = 2^nplayer;
  nums = numx / numa;
  nobs = size(sobs, 1);

  % ------------------------------------------------
  % Matrix with observed indexes of state variables
  % ------------------------------------------------
  indsobs = (repmat(sobs, 1, nums)==repmat(sval', nobs, 1))*[ 1:nums ]';
  twop = kron(ones(nobs, 1), 2.^(nplayer-[1:nplayer]));
  indobs = sum(aobs_1.*twop, 2);
  indobs = (indsobs-1).*(2^nplayer) + indobs + 1;

  % --------------------------------------------------------
  % Generate random draws for a[t] (given s[t],a[t-1])
  % --------------------------------------------------------
  pchoice = pchoice(indobs,:);
  uobs = rand(nobs, nplayer);
  aobs = (uobs <= pchoice);

  % --------------------------------------------------------
  % Simulate state transitions for s[t] given s[t-1]
  % --------------------------------------------------------
  pmat1 = cumsum(ptrans, 2);
  pmat0 = cumsum([ zeros(nums,1), ptrans(:,1:nums-1) ], 2);
  pbuff1 = pmat1(sobs,:);
  pbuff0 = pmat0(sobs,:);
  uobs = rand(nobs, 1);
  uobs = kron(uobs, ones(1, nums));
  uobs = (uobs>=pbuff0) .* (uobs<=pbuff1);
  sobs = uobs * [ 1:nums ]';
end
