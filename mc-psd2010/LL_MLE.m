function [f,g]=LL_MLE(theta,data)

%LL_MLE: log likelihood for MLE

%uses analytic equilibrium solution

v=theta./(1-theta);
Dv=1./((1-theta).^2);

temp=data.a.*log(1+v)+(1-data.a).*log(-v);

f=-mean(temp);

if nargin>1

    temp=(data.a./(1+v)+(1-data.a)./(-v)).*Dv;  %n x k
    g=-mean(temp,1)';

end

end
