function f = myCDF(x)

%myCDF: CDF from Pesendorfer and Schmidt-Dengler (2010)

%simga=1

alpha=1e-10;

if x>=alpha && x<1-alpha

    f=x;

else

    if x<alpha

        f=2*alpha*normcdf(x-alpha);

    else

        f=1-alpha+2*alpha*(normcdf(x-1+alpha)-0.5);

    end

end

end
