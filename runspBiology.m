function [y, y_rst, time] = runspBiology(x0, dt, kend1, kend, ns)
%runBiology
%%
% Model to frame the different autotrophic ammonium oxidizer guilds
% Simple chemostat model with initially well defined chemistry
% The main factor examined is the behaviour of guilds under different NH4/
% DON concentrations
% I consider six guilds that make up the physiology of different
% nitrifiers. These guilds are made up by actual AOB species and one AOA 
% guild  
% The involved biogeochemistries are:
% biological reactions
% NH4(+) + 3/2O2(aq) -> NO2(-) + H2O + 2H(+)
% a fraction of the NH4 uptaken will be assimilated as biomass, which
% supports the uptake of CO2 (via rbcL).
% The rate of CO2 fixation (and biomass development) is dependent on
% achieving a minimal CN quota of 6.6. Energy generation for fixing CO2 is
% from oxidizing ammonia therefore biomass development is directly tied to
% oxidation rates.
global par;
global env;
global par0;
global diag;
par0=par(ns);
xsiz = length(x0);
y   = zeros((kend-kend1)/par0.hist_frq,xsiz);
time= zeros((kend-kend1)/par0.hist_frq,1);
t = kend1*dt;
%set initial conditions
k = kend1+1;
k2 =kend1+1;
k3 =kend1+1;
k4 =kend1+1;
k5 =kend1+1;
k6 =kend1+1;
k7 =kend1+1;
k8 =kend1+1;
x(1,:) = x0;
y(1,:) = x0;
kk = -par0.spinup+1;

while(1)
    if(kk>0)
        %do pulse fertilization
        if(k2 <= length(par0.nh4_add_time))        
            if(kk == par0.nh4_add_time(k2))            
                x(1,par0.id_nh3x)=x(1,par0.id_nh3x)+par0.nh4_add_amount(k2);            
                k2=k2+1;        
            end            
        end        
        if(k3 <= length(par0.o2_add_time))        
            if(kk == par0.o2_add_time(k3))            
                x(1,par0.id_o2x)=x(1,par0.id_o2x)+par0.o2_add_amount(k3);            
                k3=k3+1;        
            end            
        end
        if(k4 <= length(par0.no2_add_time))
            if(kk == par0.no2_add_time(k4))
                x(1,par0.id_no2x)=x(1,par0.id_no2x)+par0.no2_add_amount(k4);
                k4=k4+1;
            end
        end
        if(k5 <= length(par0.DON_add_time))
           if(kk == par0.DON_add_time(k5))
              x(1,par0.id_DON)=x(1,par0.id_DON)+par0.DON_add_amount(k5);
             k5=k5+1;
           end
        end
        if(k6 <= length(par0.glutamate_add_time))
           if(kk == par0.glutamate_add_time(k6))
              x(1,par0.id_glutamate)=x(1,par0.id_glutamate)+par0.glutamate_add_amount(k6);
             k6=k6+1;
           end
        end
        if(k7 <= length(par0.glucose_add_time))
           if(kk == par0.glucose_add_time(k7))
              x(1,par0.id_glucose)=x(1,par0.id_glucose)+par0.glucose_add_amount(k7);
             k7=k7+1;
           end
        end
        if(k8 <= length(par0.acetate_add_time))
           if(kk == par0.acetate_add_time(k8))
              x(1,par0.id_acetate)=x(1,par0.id_acetate)+par0.acetate_add_amount(k8);
             k8=k8+1;
           end
        end
    end
    if(kk<=1)
        kk1=1;
    else
        
        kk1=k-kend1;
    end
    
    %add solubility control from temperature
    bsc_nh3=5.6e1*exp(-4100.*(1./env.tsoi(kk1)-1./298.15))*env.tsoi(kk1)/12.2;  %busen solubility, not used right now
    bsc_o2 =1.3e-3*exp(-1500.*(1./env.tsoi(kk1)-1./298.15))*env.tsoi(kk1)/12.2; %bunsen solubility of oxygen
    bsc_no =1.9e-3*exp(-1500.*(1./env.tsoi(kk1)-1./298.15))*env.tsoi(kk1)/12.2; %bunsen solubility of no
    bsc_don =1.9e-3*exp(-1500.*(1./env.tsoi(kk1)-1./298.15))*env.tsoi(kk1)/12.2;    
    par0.Tsoi = env.tsoi(kk1);
    par0.pH   = env.pH(kk1);
    %par0.bgflx_o2=env.bgflx_o2(kk1);
    %fraction of aqueous tracer in their bulk concentrations
    o2_aq_frac    = (env.h2ovol(kk1)*bsc_o2)/(env.h2ovol(kk1)*bsc_o2...
        +env.airvol(kk1));                                                   %aqueous o2
    no_aq_frac    = (env.h2ovol(kk1)*bsc_no)/(env.h2ovol(kk1)*bsc_no...
        +env.airvol(kk1));                                                   %aqueous no
    %n2o_aq_frac   = (env.h2ovol(kk1)*bsc_n2o)/(env.h2ovol(kk1)*bsc_n2o+env.airvol(kk1));      %aqueous n2o
    nh4_2_nh3aq   = par0.ke_nh4_nh3*10.^(-env.pH(kk1));
    
    par0.nox_2_aqu = no_aq_frac;                                            %no<->no(aq), bunsen solubility, moisture, temperature control
    par0.o2x_2_aqu = o2_aq_frac;                                            %o2<->o2(aq), bunsen solubility, moisture, temperature control 
    bsc_nh3= (1+nh4_2_nh3aq);
    par0.nh4_fr_frac= nh4_2_nh3aq*env.h2ovol(kk1)*bsc_nh3...
        /(env.h2ovol(kk1)*bsc_nh3+env.airvol(kk1))/(nh4_2_nh3aq+1);            %nh3x:= nh3(g)+nh3(aq)+nh4(+)
    par0.nh3_fr_frac=par0.nh4_fr_frac./nh4_2_nh3aq;
    %do biogeochemistry
    y1(1,:)=adptmbbks1(@biology,x(1,:),t,dt);
    %update time
    t= t+dt;
    y_rst(1,:)=y1(1,:);  %store the restart array
    if(k>kend+par0.spinup)
        break;
    end    
    k=k+1;
    if(k>par0.spinup)    
        kk=k-par0.spinup-kend1;
        if(par0.is_stamp)                    
            %time stamp sampling                
            if(mod(kk,par0.hist_frq)==0)            
                ii=ceil(kk/par0.hist_frq);            
                y(ii,:)=y1(1,:);            
                diag(ns).nh3aq(ii)=y1(1,par0.id_nh3x).*par0.nh3_fr_frac;            
                diag(ns).nh4(ii)  =y1(1,par0.id_nh3x).*par0.nh4_fr_frac;     
                time(ii)=t;
            end            
        else            
            %average                
            ii = ceil(kk/par0.hist_frq);        
            y(ii,:)=y(ii,:)+y1(1,:)./par0.hist_frq;        
            diag(ns).nh3aq(ii)=diag(ns).nh3aq(ii)+y1(1,par0.id_nh3x).*par0.nh3_fr_frac./par0.hist_frq;        
            diag(ns).nh4(ii)  =diag(ns).nh4(ii)+y1(1,par0.id_nh3x).*par0.nh4_fr_frac./par0.hist_frq;       
            time(ii)=t;
        end
    end
    
    x=y1;    
end




