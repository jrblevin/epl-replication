% --------------------------------------------------------------------------
% FREQPROB       Procedure that obtains a frequency estimation
%                of Prob(Y|X) where Y is a vector of binary
%                variables and X is a vector of discrete variables
%
% by Victor Aguirregabiria
%
% Converted from Gauss to Matlab by Jason Blevins.
%
% --------------------------------------------------------------------------
%
% FORMAT:
%
%         freqp = freqprob(yobs,xobs,xval)
%
% INPUTS:
%
%     yobs    - (nobs x q) vector with sample observations
%               of Y = Y1 ~ Y2 ~ ... ~ Yq
%
%     xobs    - (nobs x k) matrix with sample observations of X
%
%     xval    - (numx x k) matrix with the values of X for which
%               we want to estimate Prob(Y|X).
%
% OUTPUTS:
%
%     freqp   - (numx x q) vector with frequency estimates of
%               Pr(Y|X) for each value in xval.
%               Pr(Y1=1|X) ~ Pr(Y2=1|X) ~ ... ~ Pr(Yq=1|X)
%
% --------------------------------------------------------------------------

function prob1 = freqprob(yobs, xobs, xval)
  numx = size(xval, 1);
  numq = size(yobs, 2);
  numobs = size(xobs, 1);
  prob1 = zeros(numx, numq);
  for t = 1:numx
    xvalt = kron(ones(numobs,1), xval(t,:));
    selx = prod(xobs==xvalt, 2);
    denom = sum(selx);
    if (denom == 0)
      prob1(t,:) = zeros(1, numq);
    else
      numer = sum(kron(selx, ones(1,numq)) .* yobs);
      prob1(t,:) = (numer') ./ denom;
    end
  end
end
