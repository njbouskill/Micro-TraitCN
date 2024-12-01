function tscal=ratkowsky_scal(tsoi,tmin,tmax,topt,c)



if(tsoi>tmin && tsoi <tmax)    
    tscal=(tsoi-tmin).*(1-exp(c.*(tsoi-tmax)));    
    tscal_opt=(topt-tmin).*(1-exp(c.*(topt-tmax)));
    tscal=(tscal/tscal_opt).^2;  
else
    tscal=0.0;
end
