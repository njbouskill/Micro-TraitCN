function p=GetGdtScalar(aj,nJ,pmax)

Extra.nJ=nJ;
Extra.aj=aj;
Extra.iJ=1./nJ;
p = brent_zero ( 0.0, pmax, eps, eps, @gdtfun ,Extra);

end