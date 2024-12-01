function tscal=ctmi(temp,tmin,tmax,topt)

if(temp <= tmin || temp >= tmax)
   tscal = 0d0;
else
   tscal = (temp-tmax).*(temp-tmin).^2d0./((topt-tmin).*((topt-tmin).*(temp-topt)-(topt-tmax).*(topt+tmin-2d0*temp)));
end
   
   