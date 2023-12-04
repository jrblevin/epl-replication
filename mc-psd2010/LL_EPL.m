function [f,g] = LL_EPL(theta,data)
%LL_EPL: EPL log likelihood function

%Y=C+B.*theta;
Y=data.C+data.B*theta;

%v's
v1=Y(1);
v2=Y(2);

%derivatives of v's
Dv1=data.B(1,:);
Dv2=data.B(2,:);

%log likelihood
temp1=data.a1.*log(1+v1)+(1-data.a1).*log(-v1);
temp2=data.a2.*log(1+v2)+(1-data.a2).*log(-v2);

f=-((sum(temp1)+sum(temp2))./data.nobs);


%gradient
if nargin>1
    temp1=(data.a1./(1+v1)+(1-data.a1)./(-v1)).*Dv1;  %n1 x k
    temp2=(data.a2./(1+v2)+(1-data.a2)./(-v2)).*Dv2;  %n2 x k

    g=-(sum(temp1,1)+sum(temp2,1))'./data.nobs;
end

end
