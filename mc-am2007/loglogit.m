function llik = loglogit(ydum, x, b)
  myzero = 1e-12;
  expxb = exp(-x*b);
  Fxb = 1./(1+expxb);
  Fxb = Fxb + (myzero - Fxb).*(Fxb<myzero) ...
            + (1-myzero - Fxb).*(Fxb>1-myzero);
  llik = ydum'*ln(Fxb) + (1-ydum)'*ln(1-Fxb);
end
