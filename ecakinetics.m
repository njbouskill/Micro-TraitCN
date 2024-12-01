function [siej]=ecakinetics(kd,ss,ee)
%the first order approximate new enzyme substrate kinetics

[I,J]=size(kd);
siej=zeros(I,J);
dnrm=zeros(I,J);


for i = 1 : I
    dnm1=0.0;
    for k = 1 : J
        if(kd(i,k)>0)
            dnm1=dnm1+ee(k)/kd(i,k);
        end
    end
    
    for j = 1 : J
        dnm2=0.0;
        if(kd(i,j)>0)
            for k = 1 : I
                if(kd(k,j)>0.0)                
                    dnm2=dnm2+ss(k)/kd(k,j);
                end
            end    
            dnrm(i,j)=kd(i,j)*(1+dnm1+dnm2);
            siej(i,j)=ss(i)*ee(j)/dnrm(i,j);
        end
    end
end