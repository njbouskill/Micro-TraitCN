function tscal=ratkowsky(tsoi,tmin,tmax,c)

tscal=zeros(size(tsoi));
for j = 1 : length(tsoi)
    if(tsoi(j)>tmin && tsoi(j) <tmax)    
        tscal(j)=(tsoi(j)-tmin).*(1-exp(c.*(tsoi(j)-tmax)));    
        tscal(j)=tscal(j).^2;    
    end
end
tscal=tscal./max(tscal);