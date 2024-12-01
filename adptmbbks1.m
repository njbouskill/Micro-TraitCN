function y=adptmbbks1(odefun,y0,t,dt)
%first order adaptive mbbks ode integrator
y=y0;

dt2=dt;
crpdif=1.d-2;
dtmin=dt/2.^0;
dtr=dt;
neq=length(y);
while(1)
    dydt=feval(odefun,t,y);
    %identify pmax
    dt1=dt2;
    nJ=0;
    pmax=0.0;
    aj=zeros(size(dydt));
    for jj = 1 : neq
        if(dydt(jj)<0.)
            nJ=nJ+1;
            pm=-y(jj)/(dydt(jj)*dt1);
            aj(nJ)=-1./pm;
            if(nJ==1)
                pmax=pm;
            else
                pmax=min(pm,pmax);
            end
        end
    end
    if(nJ>0)
        pmax=min(1.0,pmax.^nJ);
        %solve
        p=GetGdtScalar(aj,nJ,pmax);
        if(abs(p-1.0)<crpdif || dt1<=dtmin)
            y=y+dydt.*dt1.*p.^(1./nJ);
            dtr=dtr-dt1;
            dt2=dtr;
        else
            dt2=dt1*0.5;
        end
    else
        p=1.0;
        y=y+dydt.*dt1;
        dtr=dtr-dt1;
    end

    if(dtr==0.0)
        break
    end
end
end