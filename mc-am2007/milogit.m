% -----------------------------------------------------------------
% MILOGIT - Estimation of a Logit Model by Maximum Likelihood
%           The optimization algorithm is a Newton's method.
%
% by Victor Aguirregabiria
%
% Converted from Gauss to Matlab by Jason Blevins.
%
% -----------------------------------------------------------------
%
% Format       [ best, varest ] = milogit(ydum, x)
%
% Input        ydum    - vector of observations of the dependent variable
%              x       - matrix of explanatory variables
%
% Output       best    - ML estimates
%              varest  - estimate of the covariance matrix
% -----------------------------------------------------------------

function [ b0, Avarb ] = milogit(ydum, x)
  nobs = size(ydum, 1);
  nparam = size(x, 2);
  eps1 = 1e-4;
  eps2 = 1e-2;
  b0 = zeros(nparam, 1);
  iter = 1;
  criter1 = 1000;
  criter2 = 1000;
  lsopts.SYM = true; lsopts.POSDEF = true;

  while ((criter1 > eps1) || (criter2 > eps2))
    % printf("\n");
    % printf("Iteration                = %d\n", iter);
    % printf("Log-Likelihood function  = %12.4f\n", loglogit(ydum,x,b0));
    % printf("Norm of b(k)-b(k-1)      = %12.4f\n", criter1);
    % printf("Norm of Gradient         = %12.4f\n", criter2);
    % printf("\n");
    expxb0 = exp(-x*b0);
    Fxb0 = 1./(1+expxb0);
    dlogLb0 = x'*(ydum - Fxb0);
    d2logLb0 = ( repmat(Fxb0 .* (1-Fxb0), 1, nparam) .* x )'*x;
    b1 = b0 + linsolve(d2logLb0, dlogLb0, lsopts);
    criter1 = sqrt( (b1-b0)'*(b1-b0) );
    criter2 = sqrt( dlogLb0'*dlogLb0 );
    b0 = b1;
    iter = iter + 1;
  end

  expxb0 = exp(-x*b0);
  Fxb0 = 1./(1+expxb0);
  Avarb = -d2logLb0;
  Avarb = inv(-Avarb);
  % sdb = sqrt(diag(Avarb));
  % tstat = b0./sdb;
  % llike = loglogit(ydum,x,b0);
  % numy1 = sum(ydum);
  % numy0 = nobs - numy1;
  % logL0 = numy1*log(numy1) + numy0*log(numy0) - nobs*log(nobs);
  % LRI = 1 - llike/logL0;
  % pseudoR2 = 1 - ( (ydum - Fxb0)'*(ydum - Fxb0) )/numy1;
end
