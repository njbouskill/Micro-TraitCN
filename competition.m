function scal=competition(inflx,ss,opt)
%do substrate competition
%inflx is totoal request, ss is total resources available
if(nargin==2)
    opt=1;
end
%based on incoming request inflx and available substrate ss, the function
%returns scaling factors in scal, such that sum(inflx.*scal)<=ss
scal=ones(size(inflx));
if(opt==1)
    %do equal competition
    inss=sum(inflx);%total request for resources
    if(ss>=inss)
        return;
    else
        scal=scal.*ss./(inss+eps);
    end
else
    
end