% -----------------------------------------------------------------------
%
% SIMDYGAM        Simulates data of state and decision variables
%                 from the steady-state distribution of
%                 a Markov Perfect Equilibrium in a dynamic
%                 game of firms' market entry/exit with
%                 incomplete information
%
% by Victor Aguirregabiria
%
% Converted from Gauss to Matlab by Jason Blevins.
%
% --------------------------------------------------------------------------
%
% FORMAT:
% { aobs, aobs_1, indsobs } = simdygam(nobs,pequil,psteady,mstate)
%
% INPUTS
%     nobs    - Number of simulations (markets)
%
%     pchoice - (nstate x nplayer) matrix with MPE probs of entry
%
%     psteady - (nstate x 1) vector with steady-state distribution of {s[t],a[t-1]}
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
% OUTPUTS:
%     aobs   - (nobs x nplayer) matrix with players' choices.
%
%     aobs_1 - (nobs x nplayer) matrix with players' initial states.
%
%     indsobs - (nobs x 1) vector with simulated values of s[t]
%
% ----------------------------------------------------------------------------

function [ aobs, aobs_1, sobs, xobs ] = simdygam(nobs, pchoice, psteady, mstate)

  nplay = size(pchoice, 2);
  nums = size(pchoice, 1);
  numa = 2^nplay;
  numx = nums / numa;

  % ----------------------------------------------------------------------
  % a. Generating random draws from ergodic distribution of (s[t],a[t-1])
  % ----------------------------------------------------------------------
  pbuff1 = cumsum(psteady);
  pbuff0 = cumsum([ 0; psteady(1:nums-1) ]);
  uobs = rand(nobs, 1);
  pbuff1 = kron(pbuff1, ones(1, nobs));
  pbuff0 = kron(pbuff0, ones(1, nobs));
  uobs = kron(uobs, ones(1, nums));
  uobs = (uobs>=(pbuff0')) .* (uobs<=(pbuff1'));
  xobs = uobs * [ 1:nums ]';
  sobs = mstate(xobs, 1);
  aobs_1 = mstate(xobs, 2:nplay+1);
  clear pbuff0 pbuff1;

  % --------------------------------------------------------
  % b. Generating random draws for a[t] (given s[t],a[t-1])
  % --------------------------------------------------------
  pchoice = pchoice(xobs,:);
  uobs = rand(nobs, nplay);
  aobs = (uobs <= pchoice);

end
