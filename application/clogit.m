% -------------------------------------------------------------------------------------
% CLOGIT -    Maximum Likelihood estimation of McFadden's Conditional Logit
%
%                  Optimization algorithm: Newton's method with analytical
%                  gradient and hessian
%
% by Victor Aguirregabiria
%
% Converted from Gauss to Matlab by Jason Blevins.
%
% -------------------------------------------------------------------------------------
%
% Format      [ best, varest, iter, err ] = clogit(ydum, x, restx, maxiter, xinit)
%
% Input        ydum    - (nobs x 1) vector of observations of dependet variable
%                        Categorical variable with values: {1, 2, ..., nalt}
%
%              x       - (nobs x (k * nalt)) matrix of explanatory variables
%                        associated with unrestricted parameters.
%                        First k columns correspond to alternative 1, and so on
%
%              restx   - (nobs x nalt) vector of the sum of the explanatory
%                        variables whose parameters are restricted to be
%                        equal to 1.
%
%  Output      best    - (k x 1) vector with ML estimates.
%
%              varest  - (k x k) matrix with estimate of covariance matrix
%
% -------------------------------------------------------------------------------------

function [ b0, Avarb, llike, iter, err ] = clogit(ydum, x, restx, maxiter, xinit)
  cconvb = 1e-6;
  myzero = 1e-16;
  nobs = size(ydum, 1);
  nalt = max(ydum);
  npar = size(x, 2) / nalt;
  err = 0;
  lsopts.SYM = true;
  lsopts.POSDEF = true;

  xysum = 0;
  for j = 1:nalt
    xysum = xysum + sum( repmat(ydum==j, 1, npar) .* x(:,npar*(j-1)+1:npar*j) );
  end

  iter = 1;
  criter = 1000;
  llike = 0;
  b0 = xinit;
  while (criter > cconvb && iter <= maxiter)
    % fprintf('\n');
    % fprintf('Iteration                = %d\n', iter);
    % fprintf('Log-Likelihood function  = %12.4f\n', llike);
    % fprintf('Norm of b(k)-b(k-1)      = %12.4f\n', criter);
    % fprintf('\n');

    % Computing probabilities
    phat = zeros(nobs, nalt);
    for j = 1:nalt
      phat(:,j) = x(:,npar*(j-1)+1:npar*j)*b0 + restx(:,j);
    end
    phatmax = repmat(max(phat, [], 2), 1, nalt);
    phat = phat - phatmax;
    phat = exp(phat) ./ repmat(sum(exp(phat), 2), 1, nalt);

    % Computing xmean
    sumpx = zeros(nobs, npar);
    xxm = 0;
    llike = -nobs;
    for j = 1:nalt
      xbuff = x(:,npar*(j-1)+1:npar*j);
      sumpx = sumpx + repmat(phat(:,j), 1, npar) .* xbuff;
      xxm = xxm + (repmat(phat(:,j), 1, npar).*xbuff)'*xbuff;
      llike = llike ...
              + sum( (ydum==j) ...
                     .* log( (phat(:,j) > myzero).*phat(:,j) ...
                             + (phat(:,j) <= myzero).*myzero) );
    end

    % Computing gradient
    d1llike = xysum - sum(sumpx);

    % Computing hessian
    d2llike = - (xxm - sumpx'*sumpx);

    % Gauss iteration
    try
        b1 = b0 + linsolve(-d2llike, d1llike', lsopts);
    catch ME
        err = 1;
        break;
    end
    criter = sqrt( (b1-b0)'*(b1-b0) );
    b0 = b1;
    iter = iter + 1;

    % fprintf('Coefficients  = ');
    % b1'
    % fprintf('\n');
  end

  Avarb  = inv(-d2llike);
  % sdb    = sqrt(diag(Avarb));
  % tstat  = b0./sdb;
  % numyj  = sum(kron(ones(1,nalt), ydum)==kron(ones(nobs,1),(1:nalt)));
  % logL0  = sum(numyj.*log(numyj./nobs));
  % lrindex = 1 - llike/logL0;
end
