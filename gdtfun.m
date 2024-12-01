function value=gdtfun(p,Extra)
%gradient function for mBBKS1
value=1.0;
for jj = 1 : Extra.nJ
    value=value*(1.0+Extra.aj(jj).*p.^(Extra.iJ));
end
value=value-p;
end