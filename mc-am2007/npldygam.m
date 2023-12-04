% ----------------------------------------------------------------------------
% NPLDYGAM       Procedure that estimates the structural parameters
%                of dynamic game of firms' entry/exit
%                using a Nested Pseudo-Likelihood (NPL) algorithm
%
% Based on the NPL Gauss code by Victor Aguirregabiria.
% Converted from Gauss to Matlab by Jason Blevins.
% Modified to include CCPs and (un-normalized) choice-specific value
% functions from 1-step NPL along with return values.
%
% ----------------------------------------------------------------------------
%
% FORMAT:
% [bnpl,varb,pchoice,v1npl,vnpl,iter,err,times] = npldygam(aobs,zobs,aobs_1,zval,ptranz,pchoice,bdisc,miniter,maxiter)
%
% INPUTS:
%     aobs    - (nobs x nplayer) matrix with observations of
%               firms' activity decisions (1=active ; 0=no active)
%     zobs    - (nobs x 1) vector with observations of market
%               exogenous characteristics.
%     aobs_1  - (nobs x nplayer) matrix with observations of
%               firms' initial states (1=incumbent; 0=potential entrant)
%     zval    - (numz x 1) vector with values of market characteristics
%     ptranz  - (numz x numz) matrix of transition probabilities
%               of market characteristics.
%     pchoice - ((numz*2^nplayer) x nplayer) matrix of players'
%               choice probabilities used to initialize the procedure.
%     bdisc   - Discount factor
%     miniter - Minimum number of NPL iterations
%     maxiter - Maximum number of NPL iterations
%
% OUTPUTS:
%     bnpl    - (kparam x niter) matrix with NPL estimates as rows:
%               (theta_fc_1; theta_fc_2; ... ; theta_fc_n; theta_rs; theta_rn; theta_ec).
%     varb    - (kparam x kparam) variance matrix from clogit estimation.
%     pchoice - ((numz*2^nplayer) x nplayer) matrix of players'
%               choice probabilities consistent with estimates.
%     v1npl   - (nplayer x 2 x (numz*2^nplayer)) matrix of players'
%               choice-specific value functions consistent with estimates
%               from the first iteration of NPL.
%     vnpl    - (nplayer x 2 x (numz*2^nplayer)) matrix of players'
%               choice-specific value functions consistent with the
%               converged NPL estimates.
%     iter    - Number of NPL iterations completed.
%     err     - Error flag indicator
%     times   - Vector of computational time (sec.) per iteration
%
% ----------------------------------------------------------------------------

function [ bnpl, varb, pchoice, v1npl, vnpl, llike, iter, err, times ] = npldygam(aobs, zobs, aobs_1, zval, ptranz, pchoice, bdisc, miniter, maxiter)
  % ---------------
  % Some constants
  % ---------------
  eulerc = 0.5772;
  myzero = 1e-16;
  nobs = size(aobs, 1);
  nplayer = size(aobs, 2);
  naction = 2;
  numa = 2^nplayer;
  numz = size(zval, 1);
  numx = numz*numa;
  kparam = nplayer + 3;
  critconv = (1e-2)*(1/kparam);
  bcnpl = zeros(kparam, 1);
  bnpl = [];
  times = [];
  varb = zeros(kparam, kparam);
  v1npl = zeros(numx, nplayer, naction);
  vnpl = zeros(numx, nplayer, naction);
  llike = -1e300;

  % ----------------------------------------------
  % Matrix with values of states of (s[t],a[t-1])
  % ----------------------------------------------
  aval = zeros(numa, nplayer);
  for i = 1:nplayer
    aval(:,i) = kron(ones(2^(i-1),1), kron([ 0; 1 ], ones(2^(nplayer-i),1)));
  end
  mstate = zeros(numx, nplayer+1);
  mstate(:,1) = kron(zval, ones(numa, 1));
  mstate(:,2:nplayer+1) = kron(ones(numz, 1), aval);

  % ------------------------------------------------
  % Matrix with observed indexes of state variables
  % ------------------------------------------------
  indzobs = (repmat(zobs, 1,numz)==repmat(zval', nobs, 1))*[ 1:numz ]';
  twop = kron(ones(nobs, 1), 2.^(nplayer-[1:nplayer]));
  indobs = sum(aobs_1.*twop, 2);
  indobs = (indzobs-1).*(2^nplayer) + indobs + 1;

  % -------------
  % NPL algorithm
  % -------------
  aobs = 1 + reshape(aobs, [nobs*nplayer,1]);
  u0 = zeros(numx*nplayer,kparam);
  u1 = zeros(numx*nplayer,kparam);
  e0 = zeros(numx*nplayer,1);
  e1 = zeros(numx*nplayer,1);
  v = zeros(numx,nplayer,naction);

  % Iterate until convergence criterion met
  criter = 1000;
  iter = 1;
  while ((criter > critconv) && (iter <= maxiter))
    %fprintf('\n');
    fprintf('-----------------------------------------------------\n');
    fprintf('NPL ESTIMATOR: K = %d\n', iter);
    fprintf('-----------------------------------------------------\n');
    fprintf('\n');
    tic;

    % -----------------------------------------------------------
    % (a) Matrix of transition probabilities Pr(a[t]|s[t],a[t-1])
    % -----------------------------------------------------------
    ptrana = ones(numx, numa);
    for i = 1:nplayer;
      mi = aval(:,i)';
      ppi = pchoice(:,i);
      ppi1 = repmat(ppi, 1, numa);
      ppi0 = 1 - ppi1;
      mi1 = repmat(mi, numx, 1);
      mi0 = 1 - mi1;
      ptrana = ptrana .* (ppi1 .^ mi1) .* (ppi0 .^ mi0);
    end

    % -----------------------
    %  (b) Storing I-b*F
    % -----------------------
    i_bf = kron(ptranz, ones(numa, numa)).*kron(ones(1,numz), ptrana);
    i_bf = eye(numx) - bdisc*i_bf;

    % -----------------------------------------
    %  (c) Construction of explanatory variables
    % -----------------------------------------
    uobs0 = zeros(nobs*nplayer,kparam);
    uobs1 = zeros(nobs*nplayer,kparam);
    eobs0 = zeros(nobs*nplayer,1);
    eobs1 = zeros(nobs*nplayer,1);
    for i = 1:nplayer
      % --------------------------------------------
      %  (c.1) Matrices Pr(a[t] | s[t],a[t-1], ai[t])
      % --------------------------------------------
      mi = aval(:,i)';
      ppi = pchoice(:,i);
      ppi= (ppi>=myzero).*(ppi<=(1-myzero)).*ppi ...
         + (ppi<myzero).*myzero ...
         + (ppi>(1-myzero)).*(1-myzero);
      ppi1 = repmat(ppi, 1, numa);
      ppi0 = 1 - ppi1;
      mi1 = repmat(mi, numx, 1);
      mi0 = 1 - mi1;
      ptrani = ((ppi1 .^ mi1) .* (ppi0 .^ mi0));
      ptranai0 = ptrana .* (mi0 ./ ptrani);
      ptranai1 = ptrana .* (mi1 ./ ptrani);
      clear mi;

      % ------------------------------------
      %  (c.2) Computing hi = E(ln(Sum(aj)+1))
      % ------------------------------------
      hi = aval;
      hi(:,i) = ones(numa, 1);
      hi = ptranai1 * log(sum(hi, 2));

      % ---------------------------
      %  (c.3) Creating U0 and U1
      % ---------------------------
      umat0 = zeros(numx, nplayer+3);
      umat1 = eye(nplayer);
      umat1 = umat1(i,:);
      umat1 = [ repmat(umat1, numx, 1), mstate(:,1), (-hi), (mstate(:,i+1)-1) ];
      clear hi;

      % ------------------------------
      %  (c.4) Creating sumu and sume
      % ------------------------------
      ppi1 = kron(ppi, ones(1, nplayer+3));
      ppi0 = 1 - ppi1;
      sumu = ppi0.*umat0 + ppi1.*umat1;
      sume = (1-ppi).*(eulerc-log(1-ppi)) + ppi.*(eulerc-log(ppi));
      clear ppi ppi0 ppi1;

      % -------------------
      %  (c.5) Creating ww
      % -------------------
      rhs = [ sumu, sume ];
      ww = i_bf \ rhs;
      clear sumu sume;

      % ----------------------------------
      %  (c.6) Creating utilda and etilda
      % ----------------------------------
      ptranai0 = kron(ptranz, ones(numa, numa)) .* kron(ones(1,numz), ptranai0);
      utilda0 = umat0 + bdisc*(ptranai0*ww(:,1:kparam));
      etilda0 = bdisc*(ptranai0*ww(:,kparam+1));
      clear umat0 ptranai0;

      ptranai1 = kron(ptranz, ones(numa, numa)) .* kron(ones(1,numz), ptranai1);
      utilda1 = umat1 + bdisc*(ptranai1*ww(:,1:kparam));
      etilda1 = bdisc*(ptranai1*ww(:,kparam+1));
      clear umat1 ptranai1;

      % -------------------------------------------
      %  (c.7) Creating observations uobs and eobs
      % -------------------------------------------
      uobs0((i-1)*nobs+1:i*nobs,:) = utilda0(indobs,:);
      uobs1((i-1)*nobs+1:i*nobs,:) = utilda1(indobs,:);
      eobs0((i-1)*nobs+1:i*nobs,:) = etilda0(indobs,:);
      eobs1((i-1)*nobs+1:i*nobs,:) = etilda1(indobs,:);
      u0((i-1)*numx+1:i*numx,:) = utilda0;
      u1((i-1)*numx+1:i*numx,:) = utilda1;
      e0((i-1)*numx+1:i*numx,:) = etilda0;
      e1((i-1)*numx+1:i*numx,:) = etilda1;
      clear utilda0 utilda1 etilda0 etilda1;
    end

    % ------------------------------------------
    %  (d) Pseudo Maximum Likelihood Estimation
    % ------------------------------------------
    [ thetaest, varb, llike, cl_iter, err ] = clogit(aobs, [uobs0, uobs1], [eobs0, eobs1], 100, bcnpl);
    if err > 0
        return;
    end

    % ----------------------------
    %  (e) Updating probabilities
    % ----------------------------
    for i = 1:nplayer
      v(:,i,2) = u1((i-1)*numx+1:i*numx,:)*thetaest + e1((i-1)*numx+1:i*numx,:);
      v(:,i,1) = u0((i-1)*numx+1:i*numx,:)*thetaest + e0((i-1)*numx+1:i*numx,:);
      buff = v(:,i,2) - v(:,i,1);
      pchoice_est(:,i) = exp(buff)./(1+exp(buff));
    end

    % Check for convergence after minimum number of iterations
    if (iter > miniter)
        criter1 = max(abs(thetaest - bcnpl));
        criter2 = max(max(abs(pchoice_est - pchoice)));
        criter = max(criter1, criter2);
    end

    % Save values for next iteration
    bcnpl = thetaest;
    vnpl = v;
    pchoice = pchoice_est;

    % Collect estimates
    bnpl = [ bnpl; thetaest' ];

    % Collect times
    times = [ times; toc ];

    % Separately save value function from first iteration
    if iter == 1
        v1npl = v;
    end

    disp('theta = ')
    disp(bcnpl')
    disp('llike = ')
    disp(llike)

    % Proceed to the next iteration
    iter = iter + 1;
  end

  % Undo final iteration increment
  iter = iter - 1;
end
