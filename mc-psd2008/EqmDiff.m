function f = EqmDiff(p,theta,data)
%EqmDiff: we have an equilibrium when f=0

p1=[1-p(1:4),p(1:4)];
p2=[1-p(5:8),p(5:8)];
    
%u_i=h_i*theta
h1=[ones(data.nstates,1), p2(:,2), [1;1;0;0]];
h2=[ones(data.nstates,1), p1(:,2), [1;0;1;0]];

%conditional transition matrices
F1_0=[p2(:,1), p2(:,2), zeros(data.nstates,2)]; %agent 1, a=0
F1_1=[zeros(data.nstates,2), p2(:,1), p2(:,2)]; %agent 1, a=1
F2_0=[p1(:,1), zeros(data.nstates,1), p1(:,2), zeros(data.nstates,1)]; %agent 2, a=0
F2_1=[zeros(data.nstates,1), p1(:,1), zeros(data.nstates,1), p1(:,2)]; %agent 2, a=1

%unconditional state transitions
FU=[p1(:,1).*p2(:,1), p1(:,1).*p2(:,2), p1(:,2).*p2(:,1), p1(:,2).*p2(:,2)];

%note: not actually psi . . . this is e(p)
psi1=normpdf(norminv(p1))./(2.*p1);
psi2=normpdf(norminv(p2))./(2.*p2);

%build H matrices
DG=eye(data.nstates)-data.beta.*FU;
H1=h1+data.beta.*(F1_1-F1_0)*(DG\(p1(:,2).*h1));
z1=-data.u1_0+data.beta.*(F1_1-F1_0)*(DG\(p1(:,1).*data.u1_0+sum(p1.*psi1,2)));
H2=h2+data.beta.*(F2_1-F2_0)*(DG\(p2(:,2).*h2));
z2=-data.u2_0+data.beta.*(F2_1-F2_0)*(DG\(p2(:,1).*data.u2_0+sum(p2.*psi2,2)));

%new choice probs
v1diff=H1*theta+z1;
v2diff=H2*theta+z2;

f=p-[normcdf(v1diff);normcdf(v2diff)];

end

