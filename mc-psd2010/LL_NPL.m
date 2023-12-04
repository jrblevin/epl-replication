function f = LL_NPL(theta,data)
%LL_NPL: NPL log likelihood function

% v's
v1=theta.*data.p20;
v2=theta.*data.p10;

% Log Likelihood
% temp1=data.a1.*log(1-myCDF(-v1))+(1-data.a1).*log(myCDF(-v1));
% temp2=data.a2.*log(1-myCDF(-v2))+(1-data.a2).*log(myCDF(-v2));

% Log Likelihood with bounded p's
myzero = 1.0e-6;
p1 = min(max(myCDF(-v1), myzero), 1-myzero);
p2 = min(max(myCDF(-v2), myzero), 1-myzero);
p1m = 1 - p1;
p2m = 1 - p2;
temp1 = data.a1.*log(p1m) + (1 - data.a1).*log(p1);
temp2 = data.a2.*log(p2m) + (1 - data.a2).*log(p2);

f = -((sum(temp1)+sum(temp2))./data.nobs);

end
