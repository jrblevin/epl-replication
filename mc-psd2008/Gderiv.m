function f=Gderiv(theta,v,data)
%Gderiv: derivative of Gfun w.r.t. v

%note: does not match DG in the paper
%In paper, DG = eye - Gderiv

zeromat=zeros(data.nstates,1);

v1=[v(1:4),v(5:8)];
v2=[v(9:12),v(13:16)];

v1diff=v1(:,2)-v1(:,1);
v2diff=v2(:,2)-v2(:,1);

p1=[1-normcdf(v1diff), normcdf(v1diff)];
p2=[1-normcdf(v2diff), normcdf(v2diff)];

%derivatives of choice probabilities
dp1dv1_0=[normpdf(v1diff),-normpdf(v1diff)];
dp1dv1_1=-dp1dv1_0;
dp2dv2_0=[normpdf(v2diff),-normpdf(v2diff)];
dp2dv2_1=-dp2dv2_0;

%for reference
% %u[i]_1=h[i]*theta
% h1=[ones(data.nstates,1), p2(:,2), [1;1;0;0]];
% h2=[ones(data.nstates,1), p1(:,2), [1;0;1;0]];

%conditional transition matrices
F1_0=[p2, zeromat, zeromat]; %agent 1, a=0
F1_1=[zeromat, zeromat, p2]; %agent 1, a=1
F2_0=[p1(:,1), zeromat, p1(:,2), zeromat]; %agent 2, a=0
F2_1=[zeromat, p1(:,1), zeromat, p1(:,2)]; %agent 2, a=1

% %for reference
% %G1: player 1
% G1_0=data.u1_0+data.beta.*(F1_0*normsurplus(v1));
% G1_1=h1*theta+data.beta.*(F1_1*normsurplus(v1));
% G2_0=data.u2_0+data.beta.*(F2_0*normsurplus(v2));
% G2_1=h2*theta+data.beta.*(F2_1*normsurplus(v2));

%dG1_0dv
temp1=data.beta.*(F1_0*diag(p1(:,1)));
temp2=data.beta.*(F1_0*diag(p1(:,2)));
temp3=data.beta.*([dp2dv2_0,zeromat,zeromat]*normsurplus(v1));
temp4=data.beta.*([dp2dv2_1,zeromat,zeromat]*normsurplus(v1));
dG1_0dv=[temp1,temp2,diag(temp3),diag(temp4)];

%dG1_1dv
temp1=data.beta.*(F1_1*diag(p1(:,1)));
temp2=data.beta.*(F1_1*diag(p1(:,2)));
temp3=[zeromat,dp2dv2_0(:,2),zeromat]*theta...
    +data.beta.*([zeromat,zeromat,dp2dv2_0]*normsurplus(v1));
temp4=[zeromat,dp2dv2_1(:,2),zeromat]*theta...
    +data.beta.*([zeromat,zeromat,dp2dv2_1]*normsurplus(v1));
dG1_1dv=[temp1,temp2,diag(temp3),diag(temp4)];

%dG2_0dv
temp1=data.beta.*([dp1dv1_0(:,1),zeromat,dp1dv1_0(:,2),zeromat]*normsurplus(v2));
temp2=data.beta.*([dp1dv1_1(:,1),zeromat,dp1dv1_1(:,2),zeromat]*normsurplus(v2));
temp3=data.beta.*(F2_0*diag(p2(:,1)));
temp4=data.beta.*(F2_0*diag(p2(:,2)));
dG2_0dv=[diag(temp1),diag(temp2),temp3,temp4];

%dG2_1dv
temp1=[zeromat,dp1dv1_0(:,2),zeromat]*theta...
    +data.beta.*([zeromat,dp1dv1_0(:,1),zeromat,dp1dv1_0(:,2)]*normsurplus(v2));
temp2=[zeromat,dp1dv1_1(:,2),zeromat]*theta...
    +data.beta.*([zeromat,dp1dv1_1(:,1),zeromat,dp1dv1_1(:,2)]*normsurplus(v2));
temp3=data.beta.*(F2_1*diag(p2(:,1)));
temp4=data.beta.*(F2_1*diag(p2(:,2)));
dG2_1dv=[diag(temp1),diag(temp2),temp3,temp4];


%dGdv
f=[dG1_0dv;dG1_1dv;dG2_0dv;dG2_1dv];


end

