function [f,g] = LL_NPL(theta,data)
%LL_NPL: log likelihood for NPL estimation

%differenced choice-specific value functions
v1diff=data.H1*theta+data.z1;
v1diff=v1diff(data.X); %pull observations
v2diff=data.H2*theta+data.z2;
v2diff=v2diff(data.X); %pull observations


%log likelihood (normal)
p1=normcdf(v1diff);
p2=normcdf(v2diff);
myzero = 1.0e-9;
p1 = min(max(p1, myzero), 1 - myzero);
p2 = min(max(p2, myzero), 1 - myzero);
temp1=data.A1.*log(p1)+(1-data.A1).*log(1-p1);
temp2=data.A2.*log(p2)+(1-data.A2).*log(1-p2);
f=-mean([temp1;temp2]);

%build analytic gradient
if nargin>1   
    
    temp1=(data.A1.*normpdf(v1diff)./p1-(1-data.A1).*normpdf(v1diff)./(1-p1)).*data.H1(data.X,:);
    temp2=(data.A2.*normpdf(v2diff)./p2-(1-data.A2).*normpdf(v2diff)./(1-p2)).*data.H2(data.X,:);
    
    g=-mean([temp1;temp2]);
    
end
   
end

