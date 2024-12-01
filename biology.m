function dxdt=biology(t,x,dt)
%In the code, I have not tried to make the assumption that the cell
%division is in steady state, i.e. nutrient/substrate uptake for biomass
%equals the growth. Rather, the micobes are allowed to hold some extra
%nutrient and substrate before division actually occur. This could in
%theory model the delay between extracellular enzyme production and
%nutrient/substrate uptake.
%input:
%t    : time
%x    : vector of state variables
%dt   : time step
%return:
%dxdt  : changing rate of the state variables
%Note:
%%
global par0;                                                               % used to pass model parameters


dxdt=zeros(size(x));                                                       % initialize return variable

%o2_input = 1d-4;
%no3_input = 1d-4;
%glutamate_input = 1d-8;
%glucose_input = 1d-8;

%get the substrates
nh3_fr  =x(par0.id_nh3x);                                                  % bioavailable NH3
co2_fr  =x(par0.id_co2x);                                                  % bioavailable CO2 for biomass synthesis, assuming they are autotrophs
no3_fr  =x(par0.id_no3x);                                                  % free NO3(-) pool
no2_fr  =x(par0.id_no2x);                                                  % free aqueous NO2 pool
no_fr   =x(par0.id_nox)*par0.nox_2_aqu;                                    % free aqueous NO pool
o2_fr   =x(par0.id_o2x)*par0.o2x_2_aqu;                                    % free aqueous O2 pool
n2o_fr  =x(par0.id_n2ox);                                                  % free N2O(-) pool           


if(par0.id_DOC>0)
    doc_fr  =x(par0.id_DOC);                                               % DOC pool
end
if(par0.id_UREA>0)
    UREA_fr =x(par0.id_UREA);                                              % Urea pool
end
if(par0.id_AMA>0)
    AMA_fr  =x(par0.id_AMA);                                               % Amino acid pool
end
if(par0.id_DON>0)
    DON_fr  =x(par0.id_DON);                                               % DON pool
end
if(par0.id_glucose>0)
    glucose_fr = x(par0.id_glucose);
end
if(par0.id_acetate>0)
    acetate_fr = x(par0.id_acetate);
end
if(par0.id_glutamate>0)
    glutamate_fr = x(par0.id_glutamate);
end
%Calculate activities coefficient of the chemical species using BDOT
%equations. Look at activity-coefficients_BDOT.xls in folder for more
%information

capA = 0.5092;
capB = 0.3283;
bdot = 0.0410;

   %calculate ionic strength of solution
   ionicstr = 0.0;

   for(i=1:par0.id_glutamate)
           ionicstr = ionicstr + par0.charge(i)*par0.charge(i)*x(i);
   end
   
     ionicstr = 0.5*ionicstr;
     capI = sqrt(ionicstr);
     %calculate activity coefficients for each chemical species
     for(i=1:par0.id_glutamate)
            par0.activitycoeff(i) = 1.0;
           par0.activitycoeff(i) = 10.0^(-capA*par0.charge(i)*par0.charge(i)*(capI/(1+capB*par0.ionradius(i)*capI)) + bdot*capI*capI);
     end

%Thermodynamics parameter, following LaRowe et al, 2012
     
    sai = 0.12;            %membrane potential? Volts, LaRowe et al, 2012 suggested a universal value of 120mV
    gasconstant = 8.314;    %!J K-1 mol-1
    faradays = 96485;       %Faraday's constant joules per volt gram equivalent    
         
%If assuming nh3 is the substrate for aob to generate energy (through  
%passive uptake), then the pH would show its effect through regulating 
%the reaction NH3+H(+)->NH4(+),
%which thus avoids the specification of an pH function to control the AOB
%activity. This can be implemented with a substrate preference function
%which 
nh3_fr = nh3_fr .* par0.nh3_fr_frac;                                       %free nh3h2o for nitrification
ratio_no2_o2=no2_fr./max(o2_fr,1.e-20);                                    %ratio of [NO2(-)]/[O2} used to indicate the level toxicity 
ratio_no_o2 =no_fr./max(o2_fr,1.e-20);                                     %ratio of [NO]/[O2] used to indicate the level of toxicity

par1.kinetics=par0.kinetics_nob_no2;
par2.kinetics=par0.kinetics_co2;

%
%%
%NOB energy genesis and substrate assimilation
dNO2_energy_nob = zeros(par0.nb_nob,1);
dO2_energy_nob  = zeros(par0.nb_nob,1);
dNH3_biomass_nob= zeros(par0.nb_nob,1);
dCO2_biomass_nob= zeros(par0.nb_nob,1);
inhib_nob=0.0;
if(par0.id_DON>0)
    dDON_energy_nob        = zeros(par0.nb_nob,1);                         %uptake of general don (other than urea and amino acid) by aob for energy
    dDON_biomass_nob       = zeros(par0.nb_nob,1);                         %uptake of don by aob for biomass
    dO2_DON_energy_nob     = zeros(par0.nb_nob,1);                         %O2 demand for DON oxidation
    dNH3_DON_biomass_nob   = zeros(par0.nb_nob,1);                         %CO2 yield from DON oxidation
    dDON_N_biomass_nob     = zeros(par0.nb_nob,1);                         %nitrogen from DON used for keeping biomass stoichiometry
end
if(par0.do_nob)
    for j = 1 : par0.nb_nob
        %do temperature modification on the transporters
        tscal=ctmi(par0.Tsoi, ...
            par0.no2_nob_tmin(j),par0.no2_nob_tmax(j), ...
            par0.no2_nob_topt(j));
        vmax_no2_nob=par0.vmax_no2_nob(j)*tscal;
        
        %NO2(-)+0.5O2(aq)         -> NO3(-)
        par1.km=mm_km(par0.kme_no2_nob(j),vmax_no2_nob,par0.no2_nob_porter_dd(j));        
        par1.vmax=vmax_no2_nob;
        par1.km1 =par0.kmscal_no2_nob(j)+x(par0.id_nob_cell(j));           %apparent km
        if(strcmp(par1.kinetics,'haldane'))
            par1.ki=par0.ki_no2_nob(j);
        end
        %consumption NO2 for energy genesis
        %here I assume the competition from detoxification is not significant
        scal_monod=vmax_no2_nob*x(par0.id_nob_cell(j))...
            *monod(o2_fr,par0.kme_o2_nob(j));
        dNO2_energy_nob(j)=substrate_kinetics(no2_fr,par1)...
            *scal_monod;
        inhib_nob=inhib_nob+scal_monod;
        dO2_energy_nob(j) = dNO2_energy_nob(j).*0.5;
        
        %NOB autotrophic pathway, the uptake of CO2 is
        Qn =x(par0.id_nob_n(j))./x(par0.id_nob_cell(j));                   %cell based nitrogen quota    
        
        rcn=x(par0.id_nob_c(j))./x(par0.id_nob_n(j));                      %cn ratio at present time step
    
        scal_co2 = 1.0-(rcn-par0.rcn_nob_min(j))./(par0.rcn_nob_max(j)...
            -par0.rcn_nob_min(j));                                          %stoichiometry constraint on C uptake
        scal_co2 = max(scal_co2,0.0);
        scal_nh3=(1./par0.rcn_nob_min(j)-1./rcn)...                         %stoichiometry constraint on N assimilation
            /(1/par0.rcn_nob_min(j)-1./par0.rcn_nob_max(j));
        scal_nh3=max(scal_nh3,0.0);
        
        vmax_NH3_nob = 0.0;
        if(par0.id_DON>0)
            %heterotrophic pathway
            %C+O2->CO2
            km_DON_nob = par0.km_DON_nob(j);    
            dDON_energy_nob(j) = par0.vmax_DON_nob(j)*tscal...
               .*monod(DON_fr, km_DON_nob)*x(par0.id_nob_cell(j));         %for DON
            dDON_biomass_nob(j)= par0.yld_DON_nob(j).* dDON_energy_nob(j); %carbon yield for biomass
            
            dNH3_DON_biomass_nob(j)=dDON_biomass_nob(j)/rcn;               %maximum nitrogen demand from DON uptake
            % 
            dDON_N_biomass_nob(j) = dDON_energy_nob(j) *par0.NC_DON ...
                * par0.alpahN_DON_nob(j); 
            dDON_N_biomass_nob(j) = min(dNH3_DON_biomass_nob(j),dDON_N_biomass_nob(j));
            dNH3_DON_biomass_nob(j) = dNH3_DON_biomass_nob(j) - dDON_N_biomass_nob(j); %maximum nh3 demand from DON processing
            dDON_biomass_nob(j) = dDON_biomass_nob(j) * scal_co2;          %stoichiometry adjustment
            dDON_energy_nob(j) = dDON_energy_nob(j) - dDON_biomass_nob(j); 
            
            vmax_NH3_nob = dNH3_DON_biomass_nob(j);
            dO2_DON_energy_nob(j) = dDON_energy_nob(j);
        end        
        %uptake ammonia? /NO3?, NO3 is too expensive to assimilate, so I assume
        %NH4 here
        dCO2_biomass_nob(j)=dNO2_energy_nob(j)*par0.yld_co2_nob(j);        %nh4 needed for biomass synthesis
                
        vmax_NH3_nob = vmax_NH3_nob + dCO2_biomass_nob(j)./rcn;            %maximum NH3 uptake rate
        par2.km=par0.kmb_NH3_nob(j);
        par2.vmax=vmax_NH3_nob;
        par2.km1 =par0.kmscal_NH3_nob(j);
    
        dNH3_biomass_nob(j)=vmax_NH3_nob*substrate_kinetics(nh3_fr,par2)...
            .*scal_nh3;                                                    %carbon needed for growth with the energy generated

        SS    =[Qn,x(par0.id_nob_c(j))/x(par0.id_nob_cell(j))];            %state variable to compute potential rate of cell division
        SS_min=[1./par0.rcn_nob_max(j),1.0];                               %minimum quota for cell division
    
        %Cell divistion rate    
        dxdt(par0.id_nob_cell(j))=par0.gmax_cell_nob(j)*growth_droop(SS,SS_min)...
        *x(par0.id_nob_cell(j));
        if(abs(dxdt(par0.id_nob_cell(j)))>1d3)
            fprintf('growth=%f,SS=%f,%f,SS_min=%f,%f,cell=%f\n',growth_droop(SS,SS_min),SS,SS_min,x(par0.id_nob_cell(j)));
            fprintf('rcn=%f,Qn=%f,c=%f,cell=%f\n',1./par0.rcn_nob_max(j),Qn,x(par0.id_nob_c(j)),x(par0.id_nob_cell(j)));
            error('xx');
        end
        

        %fprintf('droop=%e,%e,%e,%e,%e\n',growth_droop(SS,SS_min),SS,SS_min)
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Heterotrophs energy genesis and substrate assimilation
dNO3_energy_hetero = zeros(par0.nb_hetero,1);
dO2_energy_hetero  = zeros(par0.nb_hetero,1);
dNH3_biomass_hetero= zeros(par0.nb_hetero,1);
if(par0.id_glucose>0)
    dglucose_NO3_energy_hetero        = zeros(par0.nb_hetero,1);                   %uptake of glucose for energy coupled with NO3
    dglucose_N2O_energy_hetero        = zeros(par0.nb_hetero,1);                   %uptake of glucose for energy coupled with N2O
    dglucose_O2_energy_hetero        = zeros(par0.nb_hetero,1);                   %uptake of glucose for energy coupled with O2
    dglucose_NO3_biomass_hetero       = zeros(par0.nb_hetero,1);                       %uptake of glucose for biomass coupled with NO3
    dglucose_N2O_biomass_hetero       = zeros(par0.nb_hetero,1);                       %uptake of glucose for biomass coupled with N2O
    dglucose_O2_biomass_hetero       = zeros(par0.nb_hetero,1);                       %uptake of glucose for biomass coupled with O2
    dO2_glucose_energy_hetero     = zeros(par0.nb_hetero,1);                       %O2 demand for glucose oxidation
    dNO3_glucose_energy_hetero     = zeros(par0.nb_hetero,1);                      %NO3 demand for glucose oxidation
    dN2O_glucose_energy_hetero     = zeros(par0.nb_hetero,1);                      %N2O demand for glucose oxidation
    dNH3_glucose_biomass_hetero   = zeros(par0.nb_hetero,1);                       %N demand generated from c acquition through glucose oxidation
    dglucose_N_biomass_hetero = zeros(par0.nb_hetero,1);
end
if(par0.id_acetate>0)
    dacetate_NO3_energy_hetero        = zeros(par0.nb_hetero,1);                       %uptake of acetate for energy coupled with NO3
    dacetatec_NO3_energy_hetero        = zeros(par0.nb_hetero,1);          %uptake of acetate for energy coupled with NO3,complete pathway
    dacetate_N2O_energy_hetero        = zeros(par0.nb_hetero,1);                       %uptake of acetate for energy coupled with N2O
    dacetate_O2_energy_hetero        = zeros(par0.nb_hetero,1);                       %uptake of acetate for energy with O2
    dacetate_NO3_biomass_hetero       = zeros(par0.nb_hetero,1);                       %uptake of acetate for biomass coupled with NO3
    dacetatec_NO3_biomass_hetero       = zeros(par0.nb_hetero,1);          %uptake of acetate for biomass coupled with NO3,complete pathway
    dacetate_N2O_biomass_hetero       = zeros(par0.nb_hetero,1);                       %uptake of acetate for biomass coupled with N2O
    dacetate_O2_biomass_hetero       = zeros(par0.nb_hetero,1);                       %uptake of acetate for biomass with O2
    dO2_acetate_energy_hetero     = zeros(par0.nb_hetero,1);                       %O2 demand for acetate oxidation
    dNO3_acetate_energy_hetero     = zeros(par0.nb_hetero,1);                      %NO3 demand for acetate oxidation
    dNO3_acetatec_energy_hetero     = zeros(par0.nb_hetero,1);             %NO3 demand for acetate oxidation, complete pathway
    dN2O_acetate_energy_hetero     = zeros(par0.nb_hetero,1);                      %N2O demand for acetate oxidation
    dNH3_acetate_biomass_hetero   = zeros(par0.nb_hetero,1);                       %N demand generated from c acquition through acetate oxidation
    dNH3_acetatec_biomass_hetero   = zeros(par0.nb_hetero,1);              %N demand generated from c acquition through acetate oxidation,complete pathway
    dacetate_N_biomass_hetero = zeros(par0.nb_hetero,1);
end
if(par0.id_glutamate>0)
    dglutamate_NO3_energy_hetero        = zeros(par0.nb_hetero,1);                       %uptake of glutamate for energy coupled with NO3
    dglutamatec_NO3_energy_hetero        = zeros(par0.nb_hetero,1);            %uptake of glutamate for energy coupled with NO3, complete pathway
    dglutamate_N2O_energy_hetero        = zeros(par0.nb_hetero,1);                       %uptake of glutamate for energy coupled with N2O
    dglutamate_O2_energy_hetero        = zeros(par0.nb_hetero,1);                       %uptake of glutamate for energy with O2
    dglutamate_NO3_biomass_hetero       = zeros(par0.nb_hetero,1);                       %uptake of glutamate for biomass coupled with NO3
    dglutamatec_NO3_biomass_hetero       = zeros(par0.nb_hetero,1);         %uptake of glutamate for biomass coupled with NO3, , complete pathway
    dglutamate_N2O_biomass_hetero       = zeros(par0.nb_hetero,1);                       %uptake of glutamate for biomass coupled with N2O
    dglutamate_O2_biomass_hetero       = zeros(par0.nb_hetero,1);                       %uptake of glutamate for biomass with O2
    dO2_glutamate_energy_hetero     = zeros(par0.nb_hetero,1);                       %O2 demand for glutamate oxidation
    dNO3_glutamate_energy_hetero     = zeros(par0.nb_hetero,1);                      %NO3 demand for glutamate oxidation
    dNO3_glutamatec_energy_hetero     = zeros(par0.nb_hetero,1);            %NO3 demand for glutamate oxidation, , complete pathway
    dN2O_glutamate_energy_hetero     = zeros(par0.nb_hetero,1);                      %N2O demand for glutamate oxidation
    dNH3_glutamate_biomass_hetero   = zeros(par0.nb_hetero,1);                       %N demand generated from c acquition through glutamate oxidation
    dNH3_glutamatec_biomass_hetero   = zeros(par0.nb_hetero,1);            %N demand generated from c acquition through glutamate oxidation, complete pathway
    dglutamate_N_biomass_hetero = zeros(par0.nb_hetero,1);
    dO2_glutamate_N_biomass_hetero = zeros(par0.nb_hetero,1);
    dNO3_glutamate_N_biomass_hetero = zeros(par0.nb_hetero,1);
    dN2O_glutamate_N_biomass_hetero = zeros(par0.nb_hetero,1);
    dNO3_glutamatec_N_biomass_hetero = zeros(par0.nb_hetero,1);
end
if(par0.do_hetero)
    for j = 1 : par0.nb_hetero
        
        km_no3_hetero = par0.km_no3_hetero(j);
        km_o2_hetero = par0.km_o2_hetero(j);
        km_nh3_hetero = par0.km_nh3_hetero(j);
        km_n2o_hetero = par0.km_n2o_hetero(j);
        km_glutamate_hetero = par0.km_glutamate_hetero(j);
        
        Qn =x(par0.id_hetero_n(j))./x(par0.id_hetero_cell(j));                   %cell based nitrogen demand    
        
        rcn=x(par0.id_hetero_c(j))./x(par0.id_hetero_n(j));                      %cn ratio at present time step
    
        scal_c = 1.0-(rcn-par0.rcn_hetero_min(j))./(par0.rcn_hetero_max(j)...
            -par0.rcn_hetero_min(j));                                           %stoichiometry constraint on C uptake
        scal_c = max(scal_c,0.0);
        scal_nh3=(1./par0.rcn_hetero_min(j)-1./rcn)...                         %stoichiometry constraint on N assimilation
            /(1/par0.rcn_hetero_min(j)-1./par0.rcn_hetero_max(j));
        scal_nh3=max(scal_nh3,0.0);
        
        vmax_NH3_hetero = 1d-6;
      if(par0.id_glucose>0)    %do temperature modification on the transporters for glucose uptake
            tscal=ctmi(par0.Tsoi, ...
            par0.glucose_hetero_tmin(j),par0.glucose_hetero_tmax(j), ...
            par0.glucose_hetero_topt(j));
            vmax_glucose_hetero=par0.vmax_glucose_hetero(j)*tscal;
            km_glucose_hetero = par0.km_glucose_hetero(j);  
            ki_o2_hetero = par0.ki_o2_hetero(j);
            
            %heterotrophic pathway 1a
            %6 NO3- + glucose + 6 H+ <--> 6 CO2 + 9 H2O + 3 N2O 
            keq = 7.709904E+16;   %use excel to calculate Keq first
            qql =  par0.activitycoeff(par0.id_no3x)^(-6) * par0.activitycoeff(par0.id_glucose)^(-1)* ...
                   par0.activitycoeff(par0.id_co2x)^6 * par0.activitycoeff(par0.id_n2ox)^3;
            delg = -gasconstant*par0.Tsoi * log( keq / qql );
            thermodrv = 1.0 / ( exp( (delg+faradays*sai)/(gasconstant*par0.Tsoi)) + 1.0);
            
            
            dglucose_NO3_energy_hetero(j) = vmax_glucose_hetero...
               .*monod(glucose_fr, km_glucose_hetero).*monod(no3_fr, km_no3_hetero).* ...
                x(par0.id_hetero_cell(j)).*inhibition(o2_fr, ki_o2_hetero)*thermodrv;  %for glucose-NO3, inhibition by O2
           
            dglucose_NO3_biomass_hetero(j)= par0.yld_cglucose_hetero(j).* dglucose_NO3_energy_hetero(j); %carbon yield for biomass
            
            dNH3_glucose_biomass_hetero(j)=dglucose_NO3_biomass_hetero(j)/rcn;               %maximum nitrogen demand from DON uptake
          
            dglucose_N_biomass_hetero(j) = dglucose_NO3_energy_hetero(j) *par0.NC_glucose ...
                * par0.alpahN_glucose_hetero(j); 
            dglucose_N_biomass_hetero(j) = min(dNH3_glucose_biomass_hetero(j),dglucose_N_biomass_hetero(j));
            dNH3_glucose_biomass_hetero(j) = dNH3_glucose_biomass_hetero(j) - dglucose_N_biomass_hetero(j); %maximum nh3 demand from DON processing
             
            dglucose_NO3_biomass_hetero(j) = dglucose_NO3_biomass_hetero(j) * scal_c;          %stoichiometry adjustment
            dglucose_NO3_energy_hetero(j) = dglucose_NO3_energy_hetero(j) - dglucose_NO3_biomass_hetero(j); 
            
            vmax_NH3_hetero = dNH3_glucose_biomass_hetero(j);
            dNO3_glucose_energy_hetero(j) = dglucose_NO3_energy_hetero(j).*6;     %6 moles of NO3 per mole of glucose
            
            %heterotrophic pathway 1b
            %12 N2O + glucose <--> 6 CO2 + 6 H2O + 12 N2  
            
            keq = 3.58201E+25;   %use excel to calculate Keq first
            qql =  par0.activitycoeff(par0.id_n2ox)^(-12) * par0.activitycoeff(par0.id_glucose)^(-1)* ...
                   par0.activitycoeff(par0.id_co2x)^6;
            delg = -gasconstant*par0.Tsoi * log( keq / qql );
            thermodrv = 1.0 / ( exp( (delg+faradays*sai)/(gasconstant*par0.Tsoi)) + 1.0);
            
            dglucose_N2O_energy_hetero(j) = vmax_glucose_hetero...
               .*monod(glucose_fr, km_glucose_hetero).*monod(n2o_fr, km_n2o_hetero).* ...
                x(par0.id_hetero_cell(j))*thermodrv;  %for glucose-N2O
           
            dglucose_N2O_biomass_hetero(j)= par0.yld_cglucose_hetero(j).* dglucose_N2O_energy_hetero(j); %carbon yield for biomass
            
            dNH3_glucose_biomass_hetero(j)=dglucose_N2O_biomass_hetero(j)/rcn;               %maximum nitrogen demand from DON uptake
            
            dglucose_N_biomass_hetero(j) = dglucose_N2O_energy_hetero(j) *par0.NC_glucose ...
                * par0.alpahN_glucose_hetero(j); 
            dglucose_N_biomass_hetero(j) = min(dNH3_glucose_biomass_hetero(j),dglucose_N_biomass_hetero(j));
            dNH3_glucose_biomass_hetero(j) = dNH3_glucose_biomass_hetero(j) - dglucose_N_biomass_hetero(j); %maximum nh3 demand from DON processing
             
            dglucose_N2O_biomass_hetero(j) = dglucose_N2O_biomass_hetero(j) * scal_c;          %stoichiometry adjustment
            dglucose_N2O_energy_hetero(j) = dglucose_N2O_energy_hetero(j) - dglucose_N2O_biomass_hetero(j); 
            
            vmax_NH3_hetero = vmax_NH3_hetero + dNH3_glucose_biomass_hetero(j);
            dN2O_glucose_energy_hetero(j) = dglucose_N2O_energy_hetero(j).*12;     %6 moles of N2O per mole of glucose
            
            %heterotrophic pathway 2
            %6 O2 + glucose <--> 6 CO2 + 6 H2O
            keq = 2.32519E+18;   %use excel to calculate Keq first
            qql =  par0.activitycoeff(par0.id_o2x)^(-6) * par0.activitycoeff(par0.id_glucose)^(-1)* ...
                   par0.activitycoeff(par0.id_co2x)^6;
            delg = -gasconstant*par0.Tsoi * log( keq / qql );
            thermodrv = 1.0 / ( exp( (delg+faradays*sai)/(gasconstant*par0.Tsoi)) + 1.0);
            
            dglucose_O2_energy_hetero(j) = vmax_glucose_hetero...
                  .*monod(glucose_fr, km_glucose_hetero).*monod(o2_fr, km_o2_hetero)*x(par0.id_hetero_cell(j)) * thermodrv;  %for glucose-O2
              
            dglucose_O2_biomass_hetero(j)= par0.yld_cglucose_hetero(j).* dglucose_O2_energy_hetero(j); %carbon yield for biomass assume same as glucose-NO3 for now 
            
            dNH3_glucose_biomass_hetero(j)=dglucose_O2_biomass_hetero(j)/rcn;               %maximum nitrogen demand from DON uptake
            
            dglucose_N_biomass_hetero(j) = dglucose_O2_energy_hetero(j) *par0.NC_glucose ...
                * par0.alpahN_glucose_hetero(j); 
            dglucose_N_biomass_hetero(j) = min(dNH3_glucose_biomass_hetero(j),dglucose_N_biomass_hetero(j));
            dNH3_glucose_biomass_hetero(j) = dNH3_glucose_biomass_hetero(j) - dglucose_N_biomass_hetero(j); %maximum nh3 demand from DON processing
              
            dglucose_O2_biomass_hetero(j) = dglucose_O2_biomass_hetero(j) * scal_c;          %stoichiometry adjustment
            dglucose_O2_energy_hetero(j) = dglucose_O2_energy_hetero(j) - dglucose_O2_biomass_hetero(j);
            vmax_NH3_hetero = vmax_NH3_hetero + dNH3_glucose_biomass_hetero(j);
            
            dO2_glucose_energy_hetero(j) = dglucose_O2_energy_hetero(j).*6;     %6 moles of O2 per mole of glucose
      end

      
           %vmax_NH3_hetero = vmax_NH3_hetero;
          if(par0.id_acetate>0)    %do temperature modification on the transporters for acetate uptake
            tscal=ctmi(par0.Tsoi, ...
            par0.acetate_hetero_tmin(j),par0.acetate_hetero_tmax(j), ...
            par0.acetate_hetero_topt(j));
            vmax_acetate_hetero=par0.vmax_acetate_hetero(j)*tscal;
            vmax_acetatec_hetero=par0.vmax_acetatec_hetero(j)*tscal;
            km_acetate_hetero = par0.km_acetate_hetero(j);
            
            %heterotrophic pathway 1a
            %2 NO3- + 3 H+ + acetate <--> 2 CO2 + N2O + 3 H2O 
            keq = 7.74451E+19;   %use excel to calculate Keq first
            qql =  par0.activitycoeff(par0.id_no3x)^(-2) * par0.activitycoeff(par0.id_acetate)^(-1)* ...
                   par0.activitycoeff(par0.id_co2x)^2 * par0.activitycoeff(par0.id_n2ox)^1;
            delg = -gasconstant*par0.Tsoi * log( keq / qql );
            thermodrv = 1.0 / ( exp( (delg+faradays*sai)/(gasconstant*par0.Tsoi)) + 1.0);
            dacetate_NO3_energy_hetero(j) = vmax_acetate_hetero...
               .*monod(acetate_fr, km_acetate_hetero).*monod(no3_fr, km_no3_hetero).*inhibition(o2_fr, ki_o2_hetero)*x(par0.id_hetero_cell(j))*thermodrv;  %for acetate-NO3
           
            dacetate_NO3_biomass_hetero(j)= par0.yld_cacetate_hetero(j).* dacetate_NO3_energy_hetero(j); %carbon yield for biomass
            
            dNH3_acetate_biomass_hetero(j)= dacetate_NO3_biomass_hetero(j)/rcn;               %maximum nitrogen demand from DON uptake
          
            dacetate_N_biomass_hetero(j) = dacetate_NO3_energy_hetero(j) *par0.NC_acetate ...
                * par0.alpahN_acetate_hetero(j); 
            dacetate_N_biomass_hetero(j) = min(dNH3_acetate_biomass_hetero(j),dacetate_N_biomass_hetero(j));
            dNH3_acetate_biomass_hetero(j) = dNH3_acetate_biomass_hetero(j) - dacetate_N_biomass_hetero(j); %maximum nh3 demand from DON processing

            dacetate_NO3_biomass_hetero(j) = dacetate_NO3_biomass_hetero(j) * scal_c;          %stoichiometry adjustment
            dacetate_NO3_energy_hetero(j) = dacetate_NO3_energy_hetero(j) - dacetate_NO3_biomass_hetero(j); 
            
            vmax_NH3_hetero = vmax_NH3_hetero + dNH3_acetate_biomass_hetero(j);
            dNO3_acetate_energy_hetero(j) = dacetate_NO3_energy_hetero(j).*2;     %2 moles of NO3 per mole of acetate
            
            %heterotrophic pathway 1b
            %4 N2O + H+ + acetate <--> 2 CO2 + 4 N2 + 2 H2O 
            keq = 5.74319E+26;   %use excel to calculate Keq first
            qql =  par0.activitycoeff(par0.id_n2ox)^(-4) * par0.activitycoeff(par0.id_acetate)^(-1)* ...
                   par0.activitycoeff(par0.id_co2x)^2;
            delg = -gasconstant*par0.Tsoi * log( keq / qql );
            thermodrv = 1.0 / ( exp( (delg+faradays*sai)/(gasconstant*par0.Tsoi)) + 1.0);
            dacetate_N2O_energy_hetero(j) = vmax_acetate_hetero...
               .*monod(acetate_fr, km_acetate_hetero).*monod(n2o_fr, km_n2o_hetero)*x(par0.id_hetero_cell(j))*thermodrv;  %for acetate-NO3
           
            dacetate_N2O_biomass_hetero(j)= par0.yld_cacetate_hetero(j).* dacetate_N2O_energy_hetero(j); %carbon yield for biomass
            
            dNH3_acetate_biomass_hetero(j)= dacetate_N2O_biomass_hetero(j)/rcn;               %maximum nitrogen demand from DON uptake

            dacetate_N_biomass_hetero(j) = dacetate_N2O_energy_hetero(j) *par0.NC_acetate ...
                * par0.alpahN_acetate_hetero(j); 
            dacetate_N_biomass_hetero(j) = min(dNH3_acetate_biomass_hetero(j),dacetate_N_biomass_hetero(j));
            dNH3_acetate_biomass_hetero(j) = dNH3_acetate_biomass_hetero(j) - dacetate_N_biomass_hetero(j); %maximum nh3 demand from DON processing

            dacetate_N2O_biomass_hetero(j) = dacetate_N2O_biomass_hetero(j) * scal_c;          %stoichiometry adjustment
            dacetate_N2O_energy_hetero(j) = dacetate_N2O_energy_hetero(j) - dacetate_N2O_biomass_hetero(j); 
            
            vmax_NH3_hetero = vmax_NH3_hetero + dNH3_acetate_biomass_hetero(j);
            dN2O_acetate_energy_hetero(j) = dacetate_N2O_energy_hetero(j).*4;     %4 moles of N2O per mole of acetate
            
            %heterotrophic pathway 1c
            %1.6NO3(-) + 2.6H(+) + CH3COO(-) <> 2CO2 + 0.8N2 + 2.8H20 
            keq = 3.59142E+17;   %use excel to calculate Keq first
            qql =  par0.activitycoeff(par0.id_no3x)^(-1.6) * par0.activitycoeff(par0.id_acetate)^(-1)* ...
                   par0.activitycoeff(par0.id_co2x)^2;
            delg = -gasconstant*par0.Tsoi * log( keq / qql );
            thermodrv = 1.0 / ( exp( (delg+faradays*sai)/(gasconstant*par0.Tsoi)) + 1.0);
            dacetatec_NO3_energy_hetero(j) = vmax_acetatec_hetero...
               .*monod(acetate_fr, km_acetate_hetero).*monod(no3_fr, km_no3_hetero).*inhibition(o2_fr, ki_o2_hetero)*x(par0.id_hetero_cell(j))*thermodrv;  %for acetate-NO3
           
            dacetatec_NO3_biomass_hetero(j)= par0.yld_cacetatec_hetero(j).* dacetatec_NO3_energy_hetero(j); %carbon yield for biomass
            
            dNH3_acetatec_biomass_hetero(j)= dacetatec_NO3_biomass_hetero(j)/rcn;               %maximum nitrogen demand from DON uptake
     
            dacetatec_N_biomass_hetero(j) = dacetatec_NO3_energy_hetero(j) *par0.NC_acetate ...
                * par0.alpahN_acetate_hetero(j); 
            dacetatec_N_biomass_hetero(j) = min(dNH3_acetatec_biomass_hetero(j),dacetatec_N_biomass_hetero(j));
            dNH3_acetatec_biomass_hetero(j) = dNH3_acetatec_biomass_hetero(j) - dacetatec_N_biomass_hetero(j); %maximum nh3 demand from DON processing
            
            dacetatec_NO3_biomass_hetero(j) = dacetatec_NO3_biomass_hetero(j) * scal_c;          %stoichiometry adjustment
            dacetatec_NO3_energy_hetero(j) = dacetatec_NO3_energy_hetero(j) - dacetatec_NO3_biomass_hetero(j); 
            
            vmax_NH3_hetero = vmax_NH3_hetero + dNH3_acetatec_biomass_hetero(j);
            dNO3_acetatec_energy_hetero(j) = dacetatec_NO3_energy_hetero(j).*1.6;     %2 moles of NO3 per mole of acetate
            
            %heterotrophic pathway 2
            %2 O2 + acetate + H+ <--> 2 CO2 + 2 H2O
            keq = 3.72807E+19;   %use excel to calculate Keq first
            qql =  par0.activitycoeff(par0.id_o2x)^(-2) * par0.activitycoeff(par0.id_acetate)^(-1)* ...
                   par0.activitycoeff(par0.id_co2x)^2;
            delg = -gasconstant*par0.Tsoi * log( keq / qql );
            thermodrv = 1.0 / ( exp( (delg+faradays*sai)/(gasconstant*par0.Tsoi)) + 1.0);
            dacetate_O2_energy_hetero(j) = vmax_acetate_hetero...
                  .*monod(acetate_fr, km_acetate_hetero).*monod(o2_fr, km_o2_hetero)*x(par0.id_hetero_cell(j))*thermodrv;  %for glucose-O2
              
            dacetate_O2_biomass_hetero(j)= par0.yld_cacetate_hetero(j).* dacetate_O2_energy_hetero(j); %carbon yield for biomass assume same as glucose-NO3 for now 
            
            dNH3_acetate_biomass_hetero(j)=dacetate_O2_biomass_hetero(j)/rcn;               %maximum nitrogen demand from DON uptake
            
            dacetate_N_biomass_hetero(j) = dacetate_O2_energy_hetero(j) *par0.NC_acetate ...
                * par0.alpahN_acetate_hetero(j); 
            dacetate_N_biomass_hetero(j) = min(dNH3_acetate_biomass_hetero(j),dacetate_N_biomass_hetero(j));
            dNH3_acetate_biomass_hetero(j) = dNH3_acetate_biomass_hetero(j) - dacetate_N_biomass_hetero(j); %maximum nh3 demand from DON processing
            
            dacetate_O2_biomass_hetero(j) = dacetate_O2_biomass_hetero(j) * scal_c;          %stoichiometry adjustment
            dacetate_O2_energy_hetero(j) = dacetate_O2_energy_hetero(j) - dacetate_O2_biomass_hetero(j);
            vmax_NH3_hetero = vmax_NH3_hetero + dNH3_acetate_biomass_hetero(j);
            
            dO2_acetate_energy_hetero(j) = dacetate_O2_energy_hetero(j).*2;     %2 moles of O2 per mole of acetate
            
          end
          
     if(par0.id_glutamate>0)    %do temperature modification on the transporters for glutamate uptake
            tscal=ctmi(par0.Tsoi, ...
            par0.glutamate_hetero_tmin(j),par0.glutamate_hetero_tmax(j), ...
            par0.glutamate_hetero_topt(j));
            %tscal = 0.5;
            vmax_glutamate_hetero=par0.vmax_glutamate_hetero(j)*tscal;
            vmax_glutamatec_hetero=par0.vmax_glutamatec_hetero(j)*tscal;
            km_glutamate_hetero = par0.km_glutamate_hetero(j);
            
            dglutamate_N_biomass_hetero(j) = 0.0;
            
            %heterotrophic pathway 1a -- Partial denitrification
            %4.5NO3(-) + 4.5H(+) + glutamate <> NH3 + 5CO2 + 2.25 N2O+ 5.256 H2O 
            keq = 2.24582E+16;   %use excel to calculate Keq first
            qql =  par0.activitycoeff(par0.id_no3x)^(-4.5) * par0.activitycoeff(par0.id_glutamate)^(-1)* ...
                   par0.activitycoeff(par0.id_co2x)^5 * par0.activitycoeff(par0.id_n2ox)^2.25 * par0.activitycoeff(par0.id_nh3x)^1;
            delg = -gasconstant*par0.Tsoi * log( keq / qql );
            thermodrv = 1.0 / ( exp( (delg+faradays*sai)/(gasconstant*par0.Tsoi)) + 1.0);
            dglutamate_NO3_energy_hetero(j) = vmax_glutamate_hetero...
               .*monod(glutamate_fr, km_glutamate_hetero).*monod(no3_fr, km_no3_hetero).*inhibition(o2_fr, ki_o2_hetero)*x(par0.id_hetero_cell(j))*thermodrv;  %for acetate-NO3
           
            dglutamate_NO3_biomass_hetero(j)= par0.yld_cglutamate_hetero(j).* dglutamate_NO3_energy_hetero(j); %carbon yield for biomass
            
            dNH3_glutamate_biomass_hetero(j)= dglutamate_NO3_biomass_hetero(j)/rcn;               %maximum nitrogen demand from DON uptake
            
             dNO3_glutamate_N_biomass_hetero(j) = dglutamate_NO3_energy_hetero(j) *par0.NC_glutamate ...
                * par0.alpahN_glutamate_hetero(j); 
            dNO3_glutamate_N_biomass_hetero(j) = min(dNH3_glutamate_biomass_hetero(j),dNO3_glutamate_N_biomass_hetero(j));
            dNH3_glutamate_biomass_hetero(j) = dNH3_glutamate_biomass_hetero(j) - dNO3_glutamate_N_biomass_hetero(j); %maximum nh3 demand from DON processing
            
            dglutamate_NO3_biomass_hetero(j) = dglutamate_NO3_biomass_hetero(j) * scal_c;          %stoichiometry adjustment
            dglutamate_NO3_energy_hetero(j) = dglutamate_NO3_energy_hetero(j) - dglutamate_NO3_biomass_hetero(j); 
            
            vmax_NH3_hetero = vmax_NH3_hetero + dNH3_glutamate_biomass_hetero(j);
            dNO3_glutamate_energy_hetero(j) = dglutamate_NO3_energy_hetero(j).*4.5;     %4.5 moles of NO3 per mole of acetate
            
            dglutamate_N_biomass_hetero(j) = dglutamate_N_biomass_hetero(j) + dNO3_glutamate_N_biomass_hetero(j);
            
            %heterotrophic pathway 1b -- N2O reduction
            %9 N2O + glutamate + 2 H(+) <> 5 CO2+ 3 H2O + NH4(+) + 9 N2 
            keq = 3.08112E+26;   %use excel to calculate Keq first
            qql =  par0.activitycoeff(par0.id_n2ox)^(-9) * par0.activitycoeff(par0.id_glutamate)^(-1)* ...
                   par0.activitycoeff(par0.id_co2x)^5 * par0.activitycoeff(par0.id_nh3x)^1;
            delg = -gasconstant*par0.Tsoi * log( keq / qql );
            thermodrv = 1.0 / ( exp( (delg+faradays*sai)/(gasconstant*par0.Tsoi)) + 1.0);

            dglutamate_N2O_energy_hetero(j) = vmax_glutamate_hetero...
               .*monod(glutamate_fr, km_glutamate_hetero).*monod(n2o_fr, km_n2o_hetero)*x(par0.id_hetero_cell(j))*thermodrv;  %for acetate-NO3
           
            dglutamate_N2O_biomass_hetero(j)= par0.yld_cglutamate_hetero(j).* dglutamate_N2O_energy_hetero(j); %carbon yield for biomass
            
            dNH3_glutamate_biomass_hetero(j)= dglutamate_N2O_biomass_hetero(j)/rcn;               %maximum nitrogen demand from DON uptake
            
             dN2O_glutamate_N_biomass_hetero(j) = dglutamate_N2O_energy_hetero(j) *par0.NC_glutamate ...
                * par0.alpahN_glutamate_hetero(j); 
            dN2O_glutamate_N_biomass_hetero(j) = min(dNH3_glutamate_biomass_hetero(j),dN2O_glutamate_N_biomass_hetero(j));
            dNH3_glutamate_biomass_hetero(j) = dNH3_glutamate_biomass_hetero(j) - dN2O_glutamate_N_biomass_hetero(j); %maximum nh3 demand from DON processing
            
            dglutamate_N2O_biomass_hetero(j) = dglutamate_N2O_biomass_hetero(j) * scal_c;          %stoichiometry adjustment
            dglutamate_N2O_energy_hetero(j) = dglutamate_N2O_energy_hetero(j) - dglutamate_N2O_biomass_hetero(j); 
            
            vmax_NH3_hetero = vmax_NH3_hetero + dNH3_glutamate_biomass_hetero(j);
            dN2O_glutamate_energy_hetero(j) = dglutamate_N2O_energy_hetero(j).*9;     %9 moles of N2O per mole of acetate
            
            dglutamate_N_biomass_hetero(j) = dglutamate_N_biomass_hetero(j) + dN2O_glutamate_N_biomass_hetero(j);

            %heterotrophic pathway 1c  - Canonical denitrification
            %3.6 NO3(-)+ 5.6 H(+) glutamate <> 5 CO2+ NH4(+)+ 1.8 N2 + 4.8 H2O 
            keq = 1.43875E+18;   %use excel to calculate Keq first
            qql =  par0.activitycoeff(par0.id_no3x)^(-3.6) * par0.activitycoeff(par0.id_glutamate)^(-1)* ...
                   par0.activitycoeff(par0.id_co2x)^5 * par0.activitycoeff(par0.id_nh3x)^1;
            delg = -gasconstant*par0.Tsoi * log( keq / qql );
            thermodrv = 1.0 / ( exp( (delg+faradays*sai)/(gasconstant*par0.Tsoi)) + 1.0);

            dglutamatec_NO3_energy_hetero(j) = vmax_glutamatec_hetero...
               .*monod(glutamate_fr, km_glutamate_hetero).*monod(no3_fr, km_no3_hetero).*inhibition(o2_fr, ki_o2_hetero)*x(par0.id_hetero_cell(j))*thermodrv;  %for acetate-NO3
           
            dglutamatec_NO3_biomass_hetero(j)= par0.yld_cglutamatec_hetero(j).* dglutamatec_NO3_energy_hetero(j); %carbon yield for biomass
            
            dNH3_glutamatec_biomass_hetero(j)= dglutamatec_NO3_biomass_hetero(j)/rcn;               %maximum nitrogen demand from DON uptake
            
            dNO3_glutamatec_N_biomass_hetero(j) = dglutamatec_NO3_energy_hetero(j) *par0.NC_glutamate ...
                * par0.alpahN_glutamate_hetero(j); 
            dNO3_glutamatec_N_biomass_hetero(j) = min(dNH3_glutamatec_biomass_hetero(j),dNO3_glutamatec_N_biomass_hetero(j));
            dNH3_glutamatec_biomass_hetero(j) = dNH3_glutamatec_biomass_hetero(j) - dNO3_glutamatec_N_biomass_hetero(j); %maximum nh3 demand from DON processing
            
            dglutamatec_NO3_biomass_hetero(j) = dglutamatec_NO3_biomass_hetero(j) * scal_c;          %stoichiometry adjustment
            dglutamatec_NO3_energy_hetero(j) = dglutamatec_NO3_energy_hetero(j) - dglutamatec_NO3_biomass_hetero(j); 
            
            vmax_NH3_hetero = vmax_NH3_hetero + dNH3_glutamatec_biomass_hetero(j);
            dNO3_glutamatec_energy_hetero(j) = dglutamatec_NO3_energy_hetero(j).*3.6;     %3.6 moles of NO3 per mole of glutamate
                        
            dglutamate_N_biomass_hetero(j) = dglutamate_N_biomass_hetero(j) + dNO3_glutamatec_N_biomass_hetero(j);

            %heterotrophic pathway 2 -- aerobic heterotrophy
            %4.5 O2+ glutamate + 2 H(+) <> 5 CO2 + NH4(+) + 2 H2O
            keq = 9.82831E+16;   %use excel to calculate Keq first
            qql =  par0.activitycoeff(par0.id_o2x)^(-4.5) * par0.activitycoeff(par0.id_glutamate)^(-1)* ...
                   par0.activitycoeff(par0.id_co2x)^5 * par0.activitycoeff(par0.id_nh3x)^1;
            delg = -gasconstant*par0.Tsoi * log( keq / qql );
            thermodrv = 1.0 / ( exp( (delg+faradays*sai)/(gasconstant*par0.Tsoi)) + 1.0);

            dglutamate_O2_energy_hetero(j) = vmax_glutamate_hetero...
                  .*monod(glutamate_fr, km_glutamate_hetero).*monod(o2_fr, km_o2_hetero)*x(par0.id_hetero_cell(j))*thermodrv;  %for glucose-O2
              
            dglutamate_O2_biomass_hetero(j)= par0.yld_cglutamate_hetero(j).* dglutamate_O2_energy_hetero(j); %carbon yield for biomass assume same as glucose-NO3 for now 
            
            dNH3_glutamate_biomass_hetero(j)=dglutamate_O2_biomass_hetero(j)/rcn;               %maximum nitrogen demand from DON uptake
            
            dO2_glutamate_N_biomass_hetero(j) = dglutamate_O2_energy_hetero(j) *par0.NC_glutamate ...
                * par0.alpahN_glutamate_hetero(j); 
            dO2_glutamate_N_biomass_hetero(j) = min(dNH3_glutamate_biomass_hetero(j),dO2_glutamate_N_biomass_hetero(j));
            dNH3_glutamate_biomass_hetero(j) = dNH3_glutamate_biomass_hetero(j) - dO2_glutamate_N_biomass_hetero(j); %maximum nh3 demand from DON processing
            
            dglutamate_O2_biomass_hetero(j) = dglutamate_O2_biomass_hetero(j) * scal_c;          %stoichiometry adjustment
            dglutamate_O2_energy_hetero(j) = dglutamate_O2_energy_hetero(j) - dglutamate_O2_biomass_hetero(j);
            vmax_NH3_hetero = vmax_NH3_hetero + dNH3_glutamate_biomass_hetero(j);
            
            dO2_glutamate_energy_hetero(j) = dglutamate_O2_energy_hetero(j).*4.5;     %4.5 moles of O2 per mole of acetate

            dglutamate_N_biomass_hetero(j) = dglutamate_N_biomass_hetero(j) + dO2_glutamate_N_biomass_hetero(j);

          end
      
       dNH3_biomass_hetero(j)=vmax_NH3_hetero.*monod(nh3_fr, km_nh3_hetero)*x(par0.id_hetero_cell(j))...
        .*scal_nh3;  
        
        SS    =[Qn,x(par0.id_hetero_c(j))/x(par0.id_hetero_cell(j))];            %state variable to compute potential rate of cell division
        SS_min=[1./par0.rcn_hetero_max(j),1.0];                               %minimum quota for cell division
    
        %Cell divistion rate    
        dxdt(par0.id_hetero_cell(j))=par0.gmax_cell_hetero(j)*growth_droop(SS,SS_min)...
        *x(par0.id_hetero_cell(j));
        if(abs(dxdt(par0.id_hetero_cell(j)))>1d3)
            fprintf('growth=%f,SS=%f,%f,SS_min=%f,%f,cell=%f\n',growth_droop(SS,SS_min),SS,SS_min,x(par0.id_hetero_cell(j)));
            fprintf('rcn=%f,Qn=%f,c=%f,cell=%f\n',1./par0.rcn_hetero_max(j),Qn,x(par0.id_hetero_c(j)),x(par0.id_hetero_cell(j)));
            error('xx');
        end
        
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AOB energy genesis and substrate assimilation
dNH3_energy_aob        = zeros(par0.nb_aob,1);                              %nh4 demand for energy genesis
dNH3_biomass_aob       = zeros(par0.nb_aob,1);                              %nh4 demand for biomass synthesis
dNH3_hydroxyl_aob      = zeros(par0.nb_aob,1);                              %NH4 consumption through hydroxylamine decomposition
dNO2_detox_aob         = zeros(par0.nb_aob,1);                              %NO2 consumption through detoxification
dNO_detox_aob          = zeros(par0.nb_aob,1);                              %NO consumption through detoxification
if(par0.id_UREA>0)
    dUREA_energy_aob       = zeros(par0.nb_aob,1);                              %uptake of urea by aob for energy
    dUREA_biomass_aob      = zeros(par0.nb_aob,1);                              %uptake of urea by aob for biomass
    dO2_UREA_energy_aob    = zeros(par0.nb_aob,1);                              %O2 demand for UREA oxidation
    dCO2_UREA_biomass_aob  = zeros(par0.nb_aob,1);                              %CO2 yield from UREA oxidation
end
if(par0.id_AMA>0)
    dAMA_energy_aob        = zeros(par0.nb_aob,1);                              %uptake of amino acid by aob for energy
    dAMA_biomass_aob       = zeros(par0.nb_aob,1);                              %uptake of amino acid by aob for biomass
    dO2_AMA_energy_aob     = zeros(par0.nb_aob,1);                              %O2 demand for AMA oxidation
    dCO2_AMA_biomass_aob   = zeros(par0.nb_aob,1);                              %CO2 yield from AMA oxidation
end
if(par0.id_DON>0)
    dDON_energy_aob        = zeros(par0.nb_aob,1);                              %uptake of general don (other than urea and amino acid) by aob for energy
    dDON_biomass_aob       = zeros(par0.nb_aob,1);                              %uptake of don by aob for biomass
    dO2_DON_energy_aob     = zeros(par0.nb_aob,1);                              %O2 demand for DON oxidation
    dCO2_DON_biomass_aob   = zeros(par0.nb_aob,1);                              %CO2 yield from DON oxidation
end
dO2_NH3_hydroxyl_aob   = zeros(par0.nb_aob,1);                              %O2 demand for N2O production from NH3 thru hydroxylamine decomposition
dO2_NH3_energy_aob     = zeros(par0.nb_aob,1);                              %O2 demand for NH4 oxidation


dCO2_NH3_biomass_aob   = zeros(par0.nb_aob,1);                              %CO2 yield from NH4 oxidation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%note about substrate uptake by AOB
%for NH3 and NH4(+)
%because AOB produces H(+) during the energy production step, therefore, by
%assuming passive uptake of substrate, uptake NH4(+) is less efficient than
%uptake NH3H2O. This makes the affinity to NH4(+) lower at low pH. However,
%at low pH, NH4(+) is more abundant than NH3H2O. At this stage the uptake
%of DON will be turned on. In principle, under low pH, uptake DON will be
%easier than uptake NH4(+). Therefore, I would assume the NH4 used for
%energy production and biomass synthesis are all from DON. This will also
%make the competition strategy linear, which can be easily solved. 
%Also, for simplicity, it assumed each energy production process is related
%to a request for CO2 uptake, which are scaled by the same factor that
%determined on the CN stoichiometry of AOB. This also makes it easier to 
%implement the nitrogen downregulation on AOB.

%inhib_nob,x(par0.id_aob_cell)

par1.kinetics=par0.kinetics_aob_nh3;
%assumming the microbes are not C and O2 limited
for j = 1 : par0.nb_aob
    %compute the nh4 uptake rate for each guild
    %NH4(+)   + 3/2O2(aq)     -> NO2(-)  + H2O + 2H(+)      , aob,
    tscal=ctmi(par0.Tsoi, ...
        par0.nh3_aob_tmin(j),par0.nh3_aob_tmax(j),par0.nh3_aob_topt(j));

    vmax_nh3_aob=par0.vmax_nh3_aob(j)*tscal;
    par1.km=mm_km(par0.km_nh3_aob(j),vmax_nh3_aob,par0.nh3_aob_porter_dd(j));
    
   %modify km based on ratio [NO2(-)]/[O2], 
    %par1.km =par1.km*(1+par0.no2_tox_aob(j)*max(ratio_no2_o2,ratio_no_o2));
    
    par1.vmax=vmax_nh3_aob;
    par1.km1 =par0.kmscal_nh3_aob(j);
   if(strcmp(par1.kinetics,'haldane'))        
        par1.ki=par0.ki_nh3_aob(j);
  end
    
    dNH3_energy_aob(j) =vmax_nh3_aob*substrate_kinetics(nh3_fr,par1)...
        *x(par0.id_aob_cell(j))*monod(o2_fr,par0.km_o2_aob(j));   
   
    %hydroxylamine decomposition, assume it is proporational to NH3
    %oxidation to NO2
    dNH3_hydroxyl_aob(j)=dNH3_energy_aob(j)*par0.hydroxyl_dcmp_aob(j);
    dNH3_energy_aob(j)  =dNH3_energy_aob(j)-dNH3_hydroxyl_aob(j);
    dO2_NH3_energy_aob(j)  =1.5*dNH3_energy_aob(j);
    dO2_NH3_hydroxyl_aob(j)=dNH3_hydroxyl_aob(j);
    %compute traits of the microbe
    Qn =x(par0.id_aob_n(j))./x(par0.id_aob_cell(j));                       %nitrogen quota
    rcn=x(par0.id_aob_c(j))./x(par0.id_aob_n(j));                          %cn ratio 
    scal_co2 = max([1.0-(rcn-par0.rcn_aob_min(j))./(par0.rcn_aob_max(j)...   
        -par0.rcn_aob_min(j)),0.0]);                                        %stoichiometry constraint on C uptake
    scal_nbm=(1./par0.rcn_aob_min(j)-1./rcn)...                             %stoichiometry constraint on N assimilation
        /(1/par0.rcn_aob_min(j)-1./par0.rcn_aob_max(j));
    scal_nbm=max(scal_nbm,0.0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Compute biomass assimilation. In conditions of abundant substrate, the
    %uptake rate should be determined by the minimum stoichiometry.
    %I assumed the uptake of NH4 for biomass
    %synthesis is proportional to energy yield, C uptake follows the M-M kinetics
    %This makes the CN stoichiometry fluctuatable. 
    %Remark: when CO2 is not limited, dC_biomass(j)/dN_biomass =
    %aob_Cmin(j)/aob_Nmin(j);
    
    dNH3_biomass_aob(j)=dNH3_energy_aob(j).*par0.yld_nh3_aob(j);           %nh4 needed for biomass synthesis

    
    vmax_CO2_NH3_aob =dNH3_biomass_aob(j)./rcn;                            %maximum CO2 uptake rate, the uptake of CO2 is assumed to maintain the CN stoichiometry
    

    par2.km=par0.km_co2_aob(j);
    par2.vmax=vmax_CO2_NH3_aob;
    par2.km1 =par0.kmscal_co2_aob(j);
    if(strcmp(par2.kinetics,'haldane'))        
        par2.ki=par0.ki_co2_aob(j);
    end    
    
    dCO2_NH3_biomass_aob(j)=vmax_CO2_NH3_aob*substrate_kinetics(co2_fr,par2)...
        .*scal_co2;                                                        %carbon needed for growth with the energy generated
    
    
    dNH3_biomass_aob(j)=dNH3_biomass_aob(j).*scal_nbm;
    
    %do DON uptake
    %I below consider only three groups of DON, urea, amino acid, and others
    %(longer chains)
    %urea > peptides > amino acids, peptides is a special type of amino
    %the basic assumption is 
    %acids
    if(par0.id_UREA>0)
    %Urea: (NH2)2CO + H2O->2NH3+CO2    
        km_UREA_aob = par0.km_UREA_aob(j,1)*(1.0+par0.km_UREA_aob(j,2)....
          *10.^(-par0.km_UREA_aob(j,3)/par0.pH));                           %affinity parameter for UREA, a function of pH, high pH low affinity
        dUREA_energy_aob(j) = par0.vmax_UREA_aob(j).*tscal...
          *monod(UREA_fr,km_UREA_aob)...
          *x(par0.id_aob_cell(j));                                         %amount of UREA being captured to process nitrogen for energy

        dUREA_biomass_aob(j) = par0.yld_UREA_aob(j).* dUREA_energy_aob(j); %amount of urea-N assimilated into biomass     
    
        vmax_CO2_UREA_aob = dUREA_biomass_aob(j)*2./rcn;                   %because one Urea has two -NH2
        par2.vmax=vmax_CO2_NH3_aob;
        dCO2_UREA_biomass_aob(j)=vmax_CO2_UREA_aob...
           *substrate_kinetics(co2_fr,par2).*scal_co2;                     %carbon needed for growth with the energy generated
    
        dUREA_biomass_aob(j)=dUREA_biomass_aob(j).*scal_nbm;               %because one Urea molecule has two NH2- bond
        dO2_UREA_energy_aob(j) = dUREA_energy_aob(j).*3.;                   %so the oxygen consumption is 3 mole per mole urea consumed
    end
    if(par0.id_AMA>0)
    %amino acid
        km_AMA_aob = par0.km_AMA_aob(j,1).*(1.0+par0.km_AMA_aob(j,2)....
           *10.^(-par0.km_AMA_aob(j,3)/par0.pH));                          %affinity parameter for amino acid, a function of pH
        dAMA_energy_aob(j)  = par0.vmax_AMA_aob(j)*tscal...
            *monod(AMA_fr, km_AMA_aob)...
           *x(par0.id_aob_cell(j));                                        %for amino acid,   
        dAMA_biomass_aob(j) = par0.yld_AMA_aob(j)*dAMA_energy_aob(j);
        vmax_CO2_AMA_aob = dAMA_biomass_aob(j)*par0.NC_AMA/rcn;            %
        par2.vmax=vmax_CO2_AMA_aob;
    
        dCO2_AMA_biomass_aob(j)=vmax_CO2_AMA_aob*substrate_kinetics(co2_fr,par2)...
          .*scal_co2;                                                        %carbon needed for growth with the energy generated
        dAMA_biomass_aob(j) = dAMA_biomass_aob(j).*scal_nbm;
    
        dO2_AMA_energy_aob(j)=dAMA_energy_aob(j)*1.5*par0.NC_AMA;
    end
    if(par0.id_DON>0)
    %DON other than amino acid and urea
        km_DON_aob = par0.km_DON_aob(j,1).*(1+par0.km_DON_aob(j,2)...
           .*10.^(-par0.km_DON_aob(j,3)/par0.pH));
    
        dDON_energy_aob(j) = par0.vmax_DON_aob(j)*tscal...
            .*monod(DON_fr, km_DON_aob)...
          *x(par0.id_aob_cell(j));                                         %for DON
        dDON_biomass_aob(j)= par0.yld_DON_aob(j).* dDON_energy_aob(j);
        vmax_CO2_DON_aob = dDON_biomass_aob(j)*par0.NC_DON/rcn;
    
        par0.vmax=vmax_CO2_DON_aob;
        dCO2_DON_biomass_aob(j)=vmax_CO2_DON_aob...
           *substrate_kinetics(co2_fr,par2).*scal_co2;                     %carbon needed for growth with the energy generated
    
        dDON_biomass_aob(j)= dDON_biomass_aob(j).*scal_nbm;                %DON assimilated for nitrogen biomass
        dO2_DON_energy_aob(j) = dDON_energy_aob(j)*1.5*par0.NC_DON;        %oxygen consumption for DON oxidation
    %release extra carbon associated with DON, if it is needed
    end    
    
    %cell division
    SS    =[Qn,x(par0.id_aob_c(j))/x(par0.id_aob_cell(j))];                  %state variable to compute potential rate of cell division
    SS_min=[1./par0.rcn_aob_max(j),1.0];                                    %minimum quota for cell division
    %Compute the maximum growth rate
    %dcell_max=growth_max(SS,SS_min,x(par0.id_aob_cell(j)),dt);             %maximum allowed cell division rate with available quota
    
    %Cell divistion rate    
    dxdt(par0.id_aob_cell(j))=par0.gmax_cell_aob(j)*growth_droop(SS,SS_min)...
        *x(par0.id_aob_cell(j));
    
    %do nitrifier denitrification
    %Assumption 1: the carbon to donate electron comes
    %from biomass. As a result, the payment for detoxify NO2(-) is a
    %loss/consumption of assimilated biomass (kiil themselves?).
    %Assumption 2: the NO generated from detoxficiation is negligible? This
    %is not right under very low O2 environment, so a two-stage
    %denitrification is used
    %Assumption 3: if there's sufficient number of nob working with proper
    %level of oxygen provided, then the first step of detoxification will
    %be surpressed, this is represented as an inhibition on affinity,
    %higher no2/o2 ratio, higher affinity, greater value of 
    % sum(x(par0.id_nob_cell).*monod(o2_fr,par0.kme_o2_nob)) lower affinity
    %step1: NO2(-) + 0.25CH2O + H(+) -> NO + 0.25CO2 + 0.75H2O
    km_NO2_aob = par0.km_NO2_aob(j)*...
        (1+par0.no2_detox_aob(j,1)./max(ratio_no2_o2,1d-12)+...
        par0.no2_detox_aob(j,2).*inhib_nob/x(par0.id_aob_cell(j)));
    km_NO_aob  = par0.km_NO_aob(j)*(1+par0.no_detox_aob(j)./max(ratio_no_o2,1d-12));
    dNO2_detox_aob(j)=par0.NO2_vmax_aob(j)*monod(no2_fr,km_NO2_aob)*x(par0.id_aob_cell(j));  %because the C is provided from inside the cell, it will not constrain the reaction rate
    %step2: NO + 0.25CH2O  -> 0.5N2O + 0.25CO2 + 0.25H2O
    dNO_detox_aob(j) =par0.NO_vmax_aob(j)*monod(no_fr,km_NO_aob)*x(par0.id_aob_cell(j));
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%par1.kinetics=par0.kinetics_anam_nh3;
par1.kinetics=par0.kinetics_anam_o2;
%anam energy genesis and substrate assimilation
dNO2_energy_anam = zeros(par0.nb_anam,1);
dNH3_energy_anam  = zeros(par0.nb_anam,1);
dNH3_biomass_anam= zeros(par0.nb_anam,1);
dCO2_biomass_anam= zeros(par0.nb_anam,1);

if(par0.do_anam)
    for j = 1 : par0.nb_anam
        %do temperature modification on the transporters
           % Ammonium uptake 
        tscal=ctmi(par0.Tsoi, ...
            par0.nh3_anam_tmin(j),par0.nh3_anam_tmax(j), ...
            par0.nh3_anam_topt(j));
        vmax_nh3_anam=par0.vmax_nh3_anam(j)*tscal;
        ki_no2_anam = par0.ki_no2_anam(j);
        ki_o2_anam = par0.ki_o2_anam(j);
            
            % Nitrite uptake 
        tscal=ctmi(par0.Tsoi, ...
            par0.no2_anam_tmin(j),par0.no2_anam_tmax(j), ...
            par0.no2_anam_topt(j));
        vmax_no2_anam=par0.vmax_no2_anam(j)*tscal;
        
        % The anammox reaction can be summed up as 1 mole of ammonium
        % oxidized for 1.3 moles of NO2 reduced. 
        % 1NH4 + 1.3 NO2 -> 1 N2 + 0.3 NO3 + 2H2O.
        
        par1.km=mm_km(par0.km_no2_anam(j),vmax_no2_anam,par0.no2_anam_porter_dd(j));        
        par1.vmax=vmax_no2_anam;
        par1.km1 =par0.kmscal_no2_anam(j)+x(par0.id_anam_cell(j));           %apparent km
 %       if(strcmp(par1.kinetics,'haldane'))
 %           par1.ki=par0.ki_no2_anam(j);
 %       end
 %       if(strcmp(par1.kinetics,'haldane'))
 %           par1.ki=par0.ki_o2_anam(j);
 %       end
        
        %oxidation of NH3 for energy genesis
        
        scal_monod=vmax_nh3_anam*x(par0.id_anam_cell(j));
        dNH3_energy_anam(j)=monod(no2_fr,par0.km_no2_anam(j))*monod(nh3_fr,par0.km_nh3_anam(j))...
            *scal_monod.*inhibition(o2_fr, ki_o2_anam);
        
        dNO2_energy_anam(j) = dNH3_energy_anam(j).*1.3;    %based on the stoichiometric relationship
        
        %anam autotrophic pathway, the uptake of CO2 is
        Qn =x(par0.id_anam_n(j))./x(par0.id_anam_cell(j));                   %cell based nitrogen quota    
        
        rcn=x(par0.id_anam_c(j))./x(par0.id_anam_n(j));                      %cn ratio at present time step
    
        scal_co2 = 1.0-(rcn-par0.rcn_anam_min(j))./(par0.rcn_anam_max(j)...
            -par0.rcn_anam_min(j));                                          %stoichiometry constraint on C uptake
        scal_co2 = max(scal_co2,0.0);
        scal_nh3=(1./par0.rcn_anam_min(j)-1./rcn)...                         %stoichiometry constraint on N assimilation
            /(1/par0.rcn_anam_min(j)-1./par0.rcn_anam_max(j));
        scal_nh3=max(scal_nh3,0.0);
        
        vmax_NH3_anam = 0.0;
             
        % Uptake ammonia for biomass 
        dCO2_biomass_anam(j)=dNO2_energy_anam(j)*par0.yld_co2_anam(j);        %nh4 needed for biomass synthesis
                
        vmax_NH3_anam = vmax_NH3_anam + dCO2_biomass_anam(j)./rcn;            %maximum NH3 uptake rate
        par2.km=par0.kmb_NH3_anam(j);
        par2.vmax=vmax_NH3_anam;
        par2.km1 =par0.kmscal_NH3_anam(j);
    
        dNH3_biomass_anam(j)=vmax_NH3_anam*substrate_kinetics(nh3_fr,par2)...
            .*scal_nh3;                                                    %carbon needed for growth with the energy generated

        SS    =[Qn,x(par0.id_anam_c(j))/x(par0.id_anam_cell(j))];            %state variable to compute potential rate of cell division
        SS_min=[1./par0.rcn_anam_max(j),1.0];                               %minimum quota for cell division
    
        %Cell divistion rate    
        dxdt(par0.id_anam_cell(j))=par0.gmax_cell_anam(j)*growth_droop(SS,SS_min)...
        *x(par0.id_anam_cell(j));
        if(abs(dxdt(par0.id_anam_cell(j)))>1d3)
            fprintf('growth=%f,SS=%f,%f,SS_min=%f,%f,cell=%f\n',growth_droop(SS,SS_min),SS,SS_min,x(par0.id_anam_cell(j)));
            fprintf('rcn=%f,Qn=%f,c=%f,cell=%f\n',1./par0.rcn_anam_max(j),Qn,x(par0.id_anam_c(j)),x(par0.id_anam_cell(j)));
            error('xx');
        end
        

        %fprintf('droop=%e,%e,%e,%e,%e\n',growth_droop(SS,SS_min),SS,SS_min)
    end
end
%%
%substrate competition
%dNH3_demand = [dNH3_energy_aob+dNH3_biomass_aob+dNH3_hydroxyl_aob;...
%    dNH3_biomass_nob].*dt;
%dCO2_demand = [dCO2_NH3_biomass_aob];
%dO2_demand  = [dO2_NH3_energy_aob+dO2_NH3_hydroxyl_aob];
%do2dsp=1;
%if(par0.id_UREA>0)
%    dCO2_demand = [dCO2_demand;dCO2_UREA_biomass_aob];
%    dO2_demand = [dO2_demand;dO2_UREA_energy_aob];
%    dUREA_demand= (dUREA_energy_aob+dUREA_biomass_aob).*dt;
%    UREA_scal=competition(dUREA_demand,UREA_fr);
%    do2dsp=do2dsp+1;
%end
%if(par0.id_AMA>0)
%    dCO2_demand = [dCO2_demand;dCO2_AMA_biomass_aob];
%    dO2_demand = [dO2_demand;dO2_AMA_energy_aob];
%    dAMA_demand = (dAMA_energy_aob+dAMA_biomass_aob).*dt;
%    AMA_scal =competition(dAMA_demand,AMA_fr);
%    do2dsp=do2dsp+1;    
%end
%if(par0.id_DON>0)
%    dCO2_demand = [dCO2_demand;dCO2_DON_biomass_aob];
%    dO2_demand = [dO2_demand;dO2_DON_energy_aob];
%    dDON_demand = (dDON_energy_aob+dDON_biomass_aob).*dt;
%    DON_scal=competition(dDON_demand,DON_fr);
%    do2dsp=do2dsp+1;    
%end
%dCO2_demand=[dCO2_demand;dCO2_biomass_nob].*dt;
%dO2_demand =[dO2_demand;dO2_energy_nob].*dt;
%dNO2_demand = [dNO2_detox_aob; dNO2_energy_nob].*dt;

%the competation is assumed to happen among inter and intra microbes
%O2_scal  =competition(dO2_demand, o2_fr);                                  %scaling factor for O2
%NH3_scal =competition(dNH3_demand,nh3_fr);                                 %scaling factor from NH4 competition, it should hold O2_scal<=NH4_scal

%CO2_scal =competition(dCO2_demand,co2_fr);                                 %note it should hold CO2_scal <=NH4_scal, 
%NO2_scal =competition(dNO2_demand,no2_fr);        
%dealing with multiple constraints

%NH3_scal =min(NH3_scal,O2_scal([(1:par0.nb_aob),par0.nb_aob*do2dsp+(1:par0.nb_nob)]));
%eNO2=0.0;  %residual NO2 from NOB
%for j = 1 : par0.nb_nob
%    eNO2=eNO2+ max(NO2_scal(par0.nb_aob+j)- O2_scal(par0.nb_aob*do2dsp+j),0.0);
%    NO2_scal(par0.nb_aob+j) = min(NO2_scal(par0.nb_aob+j),O2_scal(par0.nb_aob*do2dsp+j)); 
%end
%redistribute the excessive NO2 for detoxification
%NO2_scal(1:par0.nb_aob)=NO2_scal(1:par0.nb_aob).*(1+eNO2./sum(NO2_scal(1:par0.nb_aob)));


%dscal=NH3_scal;
%if(par0.id_UREA>0)
%    UREA_scal=min(UREA_scal,O2_scal(1:par0.nb_aob));
%    dscal =[dscal;UREA_scal];          %the CO2 uptake is restricted by energy production and CO2 availability
%end
%if(par0.id_AMA>0)
%    AMA_scal =min(AMA_scal,O2_scal(1:par0.nb_aob));
%    dscal =[dscal;AMA_scal];          %the CO2 uptake is restricted by energy production and CO2 availability
%end
%if(par0.id_DON>0)
%    DON_scal =min(DON_scal,O2_scal(1:par0.nb_aob));
%    dscal =[dscal;DON_scal];          %the CO2 uptake is restricted by energy production and CO2 availability
%end

%CO2_scal =min(CO2_scal,dscal);

%bring down the rate
%dNH3_energy_aob      = dNH3_energy_aob .*NH3_scal(1:par0.nb_aob);
%dNH3_hydroxyl_aob    = dNH3_hydroxyl_aob.*NH3_scal(1:par0.nb_aob);
%do2dsp=1;
%if(par0.id_UREA>0)
%    dUREA_energy_aob     = dUREA_energy_aob.*UREA_scal(1:par0.nb_aob);
%    dUREA_biomass_aob    = dUREA_biomass_aob.*UREA_scal(1:par0.nb_aob);
%    dCO2_UREA_biomass_aob= dCO2_UREA_biomass_aob.*CO2_scal((1:par0.nb_aob)+par0.nb_aob*do2dsp);
%    do2dsp=do2dsp+1;    
%end
%if(par0.id_AMA>0)
%    dAMA_energy_aob      = dAMA_energy_aob.*AMA_scal(1:par0.nb_aob);
%    dAMA_biomass_aob     = dAMA_biomass_aob.*AMA_scal(1:par0.nb_aob);
%    dCO2_AMA_biomass_aob = dCO2_AMA_biomass_aob.*CO2_scal((1:par0.nb_aob)+par0.nb_aob*do2dsp);
%    do2dsp=do2dsp+1;    
%end
%if(par0.id_DON>0)
%    dDON_energy_aob      = dDON_energy_aob.*DON_scal(1:par0.nb_aob);
%    dDON_biomass_aob     = dDON_biomass_aob.*DON_scal(1:par0.nb_aob);    
%    dCO2_DON_biomass_aob = dCO2_DON_biomass_aob.*CO2_scal((1:par0.nb_aob)+par0.nb_aob*do2dsp);
%    do2dsp=do2dsp+1;    
%end

%dCO2_NH3_biomass_aob = dCO2_NH3_biomass_aob.*CO2_scal(1:par0.nb_aob);
%dCO2_biomass_nob     = dCO2_biomass_nob.*CO2_scal((1:par0.nb_nob)+par0.nb_aob*do2dsp);


%dNH3_biomass_aob = dNH3_biomass_aob.*NH3_scal(1:par0.nb_aob);
%dNH3_biomass_nob = dNH3_biomass_nob.*NH3_scal(par0.nb_aob+1:end);
%dNO2_energy_nob  = dNO2_energy_nob .*NO2_scal(par0.nb_aob+1:end);
%dNO2_detox_aob   = dNO2_detox_aob  .*NO2_scal(1:par0.nb_aob);

%%
%update state variable
%Compute the total demand and the temporal trends for the state variables

dxdt(par0.id_nh3x)=dxdt(par0.id_nh3x)-sum(dNH3_energy_aob)-sum(dNH3_hydroxyl_aob);

dxdt(par0.id_n2ox)=dxdt(par0.id_n2ox)+sum(dNH3_hydroxyl_aob).*0.5...
    +sum(dNO_detox_aob).*0.5;

dxdt(par0.id_nox) =dxdt(par0.id_nox)+sum(dNO2_detox_aob)-sum(dNO_detox_aob);
dxdt(par0.id_nitrificationflux)=sum(dNH3_hydroxyl_aob).*0.5...
    +sum(dNO_detox_aob).*0.5;
dxdt(par0.id_no2x)=dxdt(par0.id_no2x)+sum(dNH3_energy_aob);
if(par0.id_UREA>0)
    dxdt(par0.id_no2x)=dxdt(par0.id_no2x)+sum(dUREA_energy_aob).*2;
end
if(par0.id_AMA>0)
    dxdt(par0.id_no2x)=dxdt(par0.id_no2x)+sum(dAMA_energy_aob)*par0.NC_AMA;
end
if(par0.id_DON>0)
    dxdt(par0.id_no2x)=dxdt(par0.id_no2x)+sum(dDON_energy_aob)*par0.NC_DON;
end
dxdt(par0.id_no2x)=dxdt(par0.id_no2x)-sum(dNO2_energy_nob)-sum(dNO2_detox_aob);

dxdt(par0.id_no3x)=dxdt(par0.id_no3x)+sum(dNO2_energy_nob);

dxdt(par0.id_o2x) =-sum(dO2_NH3_energy_aob);
if(par0.id_UREA>0)
    dxdt(par0.id_o2x)=dxdt(par0.id_o2x)-sum(dO2_UREA_energy_aob);
    dxdt(par0.id_UREA)=-sum(dUREA_energy_aob);
end
if(par0.id_AMA>0)
    dxdt(par0.id_o2x) = dxdt(par0.id_o2x)-sum(dO2_AMA_energy_aob);
    dxdt(par0.id_AMA) =-sum(dAMA_energy_aob);
end
if(par0.id_DON>0)
    dxdt(par0.id_o2x)=dxdt(par0.id_o2x)-sum(dO2_DON_energy_aob);
    dxdt(par0.id_DON) =-sum(dDON_energy_aob);
    if(par0.do_nob)    
        dxdt(par0.id_DON) = dxdt(par0.id_DON) - sum(dDON_energy_nob)...
            - sum(dDON_biomass_nob);
        dxdt(par0.id_o2x)=dxdt(par0.id_o2x) -sum(dO2_DON_energy_nob);
    end
end
if(par0.do_nob)    
    dxdt(par0.id_o2x)=dxdt(par0.id_o2x)-sum(dO2_energy_nob);
end


%fprintf('df=%e,%e,%e,%e\n',dxdt(par0.id_nh3x),dxdt(par0.id_no2x),sum(dN_biomass),sum(dC_biomass));
%update the substrate storage
for j = 1 : par0.nb_aob
    dxdt(par0.id_aob_c(j))=dxdt(par0.id_aob_c(j))+dCO2_NH3_biomass_aob(j);
    dxdt(par0.id_aob_n(j))=dxdt(par0.id_aob_n(j))+dNH3_biomass_aob(j);
    dxdt(par0.id_co2x) = dxdt(par0.id_co2x) - dCO2_NH3_biomass_aob(j);
    if(par0.id_UREA>0)
        dxdt(par0.id_aob_c(j))=dxdt(par0.id_aob_c(j))+dCO2_UREA_biomass_aob(j);
        dxdt(par0.id_aob_n(j))=dxdt(par0.id_aob_n(j))+dUREA_biomass_aob(j).*2.0;        
        dxdt(par0.id_UREA) = dxdt(par0.id_UREA) - dUREA_biomass_aob(j);
        dxdt(par0.id_co2x) = dxdt(par0.id_co2x) - dCO2_UREA_biomass_aob(j)+...
            dUREA_biomass_aob(j)+dUREA_energy_aob(j);
    end
    
    if(par0.id_AMA>0)
        dxdt(par0.id_aob_c(j))=dxdt(par0.id_aob_c(j))+dCO2_AMA_biomass_aob(j);
        dxdt(par0.id_aob_n(j))=dxdt(par0.id_aob_n(j))+dAMA_biomass_aob(j)*par0.NC_AMA;    
        dxdt(par0.id_AMA)  = dxdt(par0.id_AMA) - dAMA_biomass_aob(j);
        dxdt(par0.id_co2x) = dxdt(par0.id_co2x) - dCO2_AMA_biomass_aob(j);
    end
    if(par0.id_DON>0)
        dxdt(par0.id_aob_c(j))=dxdt(par0.id_aob_c(j))+dCO2_DON_biomass_aob(j);
        dxdt(par0.id_aob_n(j))=dxdt(par0.id_aob_n(j))+dDON_biomass_aob(j)*par0.NC_DON;
        dxdt(par0.id_DON)  = dxdt(par0.id_DON) - dDON_biomass_aob(j);
        dxdt(par0.id_co2x) = dxdt(par0.id_co2x) - dCO2_DON_biomass_aob(j);
    end
    
    dxdt(par0.id_nh3x) = dxdt(par0.id_nh3x) - dNH3_biomass_aob(j);    
    
    
    if(par0.id_DOC>0)    
        dxdt(par0.id_DOC) = dxdt(par0.id_DOC);
        if(par0.id_AMA>0)
            dxdt(par0.id_DOC) = dxdt(par0.id_DOC)+ (dAMA_energy_aob(j)+...
                dAMA_biomass_aob(j));
        end
        if(par0.id_DON>0)
            dxdt(par0.id_DOC) = dxdt(par0.id_DOC)+(dDON_energy_aob(j)+...
                dDON_biomass_aob(j));

        end
    end
end

% ANAMMOX

if(par0.do_anam)
    for j = 1 : par0.nb_anam
        dxdt(par0.id_anam_c(j))=dxdt(par0.id_anam_c(j))+dCO2_biomass_anam(j);
        dxdt(par0.id_anam_n(j))=dxdt(par0.id_anam_n(j))+dNH3_biomass_anam(j);
        dxdt(par0.id_nh3x) = dxdt(par0.id_nh3x) - dNH3_energy_anam(j) - dNH3_biomass_anam(j);  
        dxdt(par0.id_co2x) = dxdt(par0.id_co2x) - dCO2_biomass_anam(j);
        dxdt(par0.id_no2x) = dxdt(par0.id_no2x) - dNO2_energy_anam(j);
    end
end

% NOB

if(par0.do_nob)
    for j = 1 : par0.nb_nob
        dxdt(par0.id_nob_c(j))=dxdt(par0.id_nob_c(j))+dCO2_biomass_nob(j);
        if(par0.id_DON>0)
            dxdt(par0.id_nob_c(j))=dxdt(par0.id_nob_c(j))+dDON_biomass_nob(j);
            dxdt(par0.id_nob_n(j))=dxdt(par0.id_nob_n(j))+dDON_N_biomass_nob(j);
            dxdt(par0.id_co2x) = dxdt(par0.id_co2x) + dDON_energy_nob(j);
            dxdt(par0.id_nh3x) = dxdt(par0.id_nh3x) + ...
                (dDON_energy_nob(j) + dDON_biomass_nob(j))*par0.NC_DON-...
                dDON_N_biomass_nob(j);
        end
        dxdt(par0.id_nob_n(j))=dxdt(par0.id_nob_n(j))+dNH3_biomass_nob(j);
    
        dxdt(par0.id_nh3x) = dxdt(par0.id_nh3x) - dNH3_biomass_nob(j);  
        dxdt(par0.id_co2x) = dxdt(par0.id_co2x) - dCO2_biomass_nob(j);    
    end
end

% Heterotrophs

if(par0.do_hetero)
    for j = 1 : par0.nb_hetero
        
        if(par0.id_glucose>0)
            dxdt(par0.id_hetero_c(j))=dxdt(par0.id_hetero_c(j))+ dglucose_NO3_biomass_hetero(j) + dglucose_O2_biomass_hetero(j) + dglucose_N2O_biomass_hetero(j);
            dxdt(par0.id_no3x) = dxdt(par0.id_no3x) - dNO3_glucose_energy_hetero(j);
            dxdt(par0.id_n2ox) = dxdt(par0.id_n2ox) - dN2O_glucose_energy_hetero(j) + dglucose_NO3_energy_hetero(j).*3;
            dxdt(par0.id_co2x) = dxdt(par0.id_co2x) + dglucose_NO3_energy_hetero(j).*6 + dglucose_N2O_energy_hetero(j).*6 + dglucose_O2_energy_hetero(j).*6;
            dxdt(par0.id_o2x) = dxdt(par0.id_o2x) - dO2_glucose_energy_hetero(j);
            dxdt(par0.id_glucose) = dxdt(par0.id_glucose) - dglucose_O2_energy_hetero(j)- dglucose_O2_biomass_hetero(j) - ...
                dglucose_NO3_energy_hetero(j)- dglucose_NO3_biomass_hetero(j) - dglucose_N2O_energy_hetero(j)- dglucose_N2O_biomass_hetero(j);
        end
        if(par0.id_acetate>0)
            dxdt(par0.id_hetero_c(j))=dxdt(par0.id_hetero_c(j))+ dacetate_NO3_biomass_hetero(j) + dacetatec_NO3_biomass_hetero(j) + dacetate_O2_biomass_hetero(j) + dacetate_N2O_biomass_hetero(j);
            dxdt(par0.id_no3x) = dxdt(par0.id_no3x) - dNO3_acetate_energy_hetero(j)- dNO3_acetatec_energy_hetero(j);
            dxdt(par0.id_n2ox) = dxdt(par0.id_n2ox) - dN2O_acetate_energy_hetero(j) + dacetate_NO3_energy_hetero(j);
            dxdt(par0.id_co2x) = dxdt(par0.id_co2x) + dacetate_NO3_energy_hetero(j).*2 + dacetatec_NO3_energy_hetero(j).*2 + dacetate_N2O_energy_hetero(j).*2 + dacetate_O2_energy_hetero(j).*2;
            dxdt(par0.id_o2x) = dxdt(par0.id_o2x) - dO2_acetate_energy_hetero(j);
            dxdt(par0.id_acetate) = dxdt(par0.id_acetate) - dacetate_O2_energy_hetero(j)- dacetate_O2_biomass_hetero(j) - ...
                dacetate_NO3_energy_hetero(j)- dacetate_NO3_biomass_hetero(j) - dacetatec_NO3_energy_hetero(j)- dacetatec_NO3_biomass_hetero(j) - dacetate_N2O_energy_hetero(j)- dacetate_N2O_biomass_hetero(j);
        end
        
        if(par0.id_glutamate>0)
            dxdt(par0.id_hetero_c(j))=dxdt(par0.id_hetero_c(j))+ dglutamate_NO3_biomass_hetero(j) + dglutamatec_NO3_biomass_hetero(j) + dglutamate_O2_biomass_hetero(j) + dglutamate_N2O_biomass_hetero(j);
            dxdt(par0.id_no3x) = dxdt(par0.id_no3x) - dNO3_glutamate_energy_hetero(j)- dNO3_glutamatec_energy_hetero(j);
            dxdt(par0.id_n2ox) = dxdt(par0.id_n2ox) - dN2O_glutamate_energy_hetero(j) + dglutamate_NO3_energy_hetero(j);
            dxdt(par0.id_co2x) = dxdt(par0.id_co2x) + dglutamate_NO3_energy_hetero(j).*5 + dglutamatec_NO3_energy_hetero(j).*5 + dglutamate_N2O_energy_hetero(j).*5 + dglutamate_O2_energy_hetero(j).*5;
            dxdt(par0.id_o2x) = dxdt(par0.id_o2x) - dO2_glutamate_energy_hetero(j);
            dxdt(par0.id_glutamate) = dxdt(par0.id_glutamate) - dglutamate_O2_energy_hetero(j)- dglutamate_O2_biomass_hetero(j) - ...
                dglutamate_NO3_energy_hetero(j)- dglutamate_NO3_biomass_hetero(j)-dglutamatec_NO3_energy_hetero(j)- dglutamatec_NO3_biomass_hetero(j) - dglutamate_N2O_energy_hetero(j)- dglutamate_N2O_biomass_hetero(j);
            dxdt(par0.id_denitrificationflux) = dxdt(par0.id_denitrificationflux) - dN2O_glutamate_energy_hetero(j) + dglutamate_NO3_energy_hetero(j);
            dxdt(par0.id_hetero_n(j)) = dxdt(par0.id_hetero_n(j))+ dglutamate_N_biomass_hetero(j);
        end     
         
        dxdt(par0.id_hetero_n(j))=dxdt(par0.id_hetero_n(j))+dNH3_biomass_hetero(j);
        dxdt(par0.id_nh3x) = dxdt(par0.id_nh3x) - dNH3_biomass_hetero(j)+dglutamate_O2_energy_hetero(j)+ dglutamate_O2_biomass_hetero(j) + ...
                dglutamate_NO3_energy_hetero(j)+ dglutamate_NO3_biomass_hetero(j)+dglutamatec_NO3_energy_hetero(j)+ dglutamatec_NO3_biomass_hetero(j) + dglutamate_N2O_energy_hetero(j)+ dglutamate_N2O_biomass_hetero(j);  
            
    end
end

%dxdt(par0.id_glutamate) = dxdt(par0.id_glutamate) + glutamate_input;
%dxdt(par0.id_glutamate) = 0.0;
%dxdt(par0.id_glucose) = 0.0;
%dxdt(par0.id_acetate) = 0.0;
%dxdt(par0.id_o2x) = 0.0;
%dxdt(par0.id_no3x) = 0.0;
%dxdt(par0.id_no2x) = 0.0;

%let the micobial die
%here I assume in order to do detoxification, the AOB needs to commit
%suicide in a certain sense to release C to accept the electrons. Another
%hypothesis might be only to release C, Jinyun Tang. Feb. 10, 2012.
for j = 1 : par0.nb_aob
    c_death_detox=(dNO2_detox_aob(j)+dNO_detox_aob(j))*0.25;
    
    c_death_mort=par0.kd_aob(j)*x(par0.id_aob_c(j));
    c_death = c_death_detox+c_death_mort;
    rnc = x(par0.id_aob_n(j)) / x(par0.id_aob_c(j));
    
    dxdt(par0.id_aob_c(j)) = dxdt(par0.id_aob_c(j)) - c_death; 
    dxdt(par0.id_aob_n(j)) = dxdt(par0.id_aob_n(j)) - c_death*rnc;
    dxdt(par0.id_aob_cell(j))=dxdt(par0.id_aob_cell(j))-c_death...
        *x(par0.id_aob_cell(j)) / x(par0.id_aob_c(j));

    %where does the n go from death of the bacteria?
    %for N from detoxification, I assume the nitrogen return to the NH4
    %pool, but for other natural mortality, it goes to microbial residual
    %pool, which currently is not tracked.
    
    dxdt(par0.id_nh3x) = dxdt(par0.id_nh3x) + c_death_detox.*rnc;
    dxdt(par0.id_co2x) = dxdt(par0.id_co2x) + c_death_detox;

    dxdt(par0.id_mcb_rsc)=dxdt(par0.id_mcb_rsc)+c_death_mort;
    dxdt(par0.id_mcb_rsn)=dxdt(par0.id_mcb_rsn)+c_death_mort*rnc;
    
    %recycle the nh4?
    
end
if(par0.do_nob)
    for j = 1 : par0.nb_nob
        c_death=par0.kd_nob(j)*x(par0.id_nob_c(j));
        rnc = x(par0.id_nob_n(j))/x(par0.id_nob_c(j));
        dxdt(par0.id_nob_cell(j))=dxdt(par0.id_nob_cell(j))-c_death...
            *x(par0.id_nob_cell(j))/x(par0.id_nob_c(j));

        dxdt(par0.id_nob_c(j)) = dxdt(par0.id_nob_c(j)) - c_death; 
        dxdt(par0.id_nob_n(j)) = dxdt(par0.id_nob_n(j)) - c_death*rnc; 
    
        dxdt(par0.id_mcb_rsc)=dxdt(par0.id_mcb_rsc)+c_death;
        dxdt(par0.id_mcb_rsn)=dxdt(par0.id_mcb_rsn)+c_death*rnc;    
        %recycle the nh4?
        
    end
end

if(par0.do_anam)
    for j = 1 : par0.nb_anam
        c_death=par0.kd_anam(j)*x(par0.id_anam_c(j));
        rnc = x(par0.id_anam_n(j))/x(par0.id_anam_c(j));
        dxdt(par0.id_anam_cell(j))=dxdt(par0.id_anam_cell(j))-c_death...
            *x(par0.id_anam_cell(j))/x(par0.id_anam_c(j));

        dxdt(par0.id_anam_c(j)) = dxdt(par0.id_anam_c(j)) - c_death; 
        dxdt(par0.id_anam_n(j)) = dxdt(par0.id_anam_n(j)) - c_death*rnc; 
    
        dxdt(par0.id_mcb_rsc)=dxdt(par0.id_mcb_rsc)+c_death;
        dxdt(par0.id_mcb_rsn)=dxdt(par0.id_mcb_rsn)+c_death*rnc;    
    end
end


if(par0.do_hetero)
    for j = 1 : par0.nb_hetero
        c_death=par0.kd_hetero(j)*x(par0.id_hetero_c(j));
        rnc = x(par0.id_hetero_n(j))/x(par0.id_hetero_c(j));
        dxdt(par0.id_hetero_cell(j))=dxdt(par0.id_hetero_cell(j))-c_death...
            *x(par0.id_hetero_cell(j))/x(par0.id_hetero_c(j));

        dxdt(par0.id_hetero_c(j)) = dxdt(par0.id_hetero_c(j)) - c_death; 
        dxdt(par0.id_hetero_n(j)) = dxdt(par0.id_hetero_n(j)) - c_death*rnc; 
    
        dxdt(par0.id_mcb_rsc)=dxdt(par0.id_mcb_rsc)+c_death;
        dxdt(par0.id_mcb_rsn)=dxdt(par0.id_mcb_rsn)+c_death*rnc;    
        %recycle the nh4?
        
    end
end
%fprintf('t1=%e,%e,%e\n',dxdt(par0.id_DON),sum(dDON_energy_nob)+sum(dDON_biomass_nob),...
%    sum(dDON_energy_aob)+sum(dDON_biomass_aob));

%codes below is for mass balance check
if(par0.cnbalance_check)
    dC=dxdt(par0.id_mcb_rsc)+sum(dxdt(par0.id_aob_c))+sum(dxdt(par0.id_nob_c))...
        +dxdt(par0.id_co2x);
    if(par0.id_UREA>0)
        dC=dC+dxdt(par0.id_UREA);
    end
    if(par0.id_AMA>0)
        dC=dC+dxdt(par0.id_AMA);
    end
    if(par0.id_DON>0)
        dC=dC+dxdt(par0.id_DON);        
    end
    if(par0.id_DOC>0)
        dC=dC+dxdt(par0.id_DOC);
    end    
    dN=dxdt(par0.id_mcb_rsn)+sum(dxdt(par0.id_aob_n))+sum(dxdt(par0.id_nob_n))...
        +dxdt(par0.id_nh3x)+dxdt(par0.id_no3x)+dxdt(par0.id_no2x)...
        +dxdt(par0.id_nox)+dxdt(par0.id_n2ox)*2;
    if(par0.id_UREA>0)
        dN=dN+dxdt(par0.id_UREA)*2.;
    end
    if(par0.id_AMA>0)    
        dN=dN+dxdt(par0.id_AMA)*par0.NC_AMA;
    end
    if(par0.id_DON>0)
        dN=dN+dxdt(par0.id_DON)*par0.NC_DON;        
    end
    fprintf('dC=%e,dN=%e,maxct=%e,maxnt=%e\n',dC,dN,...
    max([dxdt(par0.id_mcb_rsc),sum(dxdt(par0.id_aob_c)),sum(dxdt(par0.id_nob_c)),dxdt(par0.id_co2x)]),...
    max([dxdt(par0.id_mcb_rsn),sum(dxdt(par0.id_aob_n)),sum(dxdt(par0.id_nob_n)),dxdt(par0.id_nh3x),dxdt(par0.id_no3x),dxdt(par0.id_no2x),dxdt(par0.id_nox),dxdt(par0.id_n2ox)]));
end

%dxdt(par0.id_o2x) = dxdt(par0.id_o2x) + par0.bgflx_o2;

end