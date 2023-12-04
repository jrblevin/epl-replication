function f = normsurplus(v)
%normsurplus: social surplus function with normal distribution

vdiff=v(:,2)-v(:,1);

p=[1-normcdf(vdiff),normcdf(vdiff)];  %Pr(a=1)

% %social surplus
% f=(1-p).*v(:,1)+p.*v(:,2)+normpdf(vdiff);

%e(p) (not actually psi)
psi=normpdf(norminv(p))./(2.*p);

f=sum(p.*(v+psi),2);



end

