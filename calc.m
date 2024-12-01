function c=calc(tmin,tmax,topt)
c=0.3;

while(1)
    f=1-exp(c*(topt-tmax))*(1+c*(topt-tmin));
    fp=exp(c*(topt-tmax))*(-(topt-tmax)*(1+c*(topt-tmin))-(topt-tmin));
    dc=-f/fp;    
    c=c+dc;
    if(abs(dc/c)<1d-6)
        break;
    end
end