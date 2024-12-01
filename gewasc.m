%function gewasc(targetvar,filenum)
close all;
clear all;
clc;
tstart = cputime;
%%
%01/12/15
global par;
global env;
global diag;
addpath('../sampling/');

%%
%simulator setup
nsamp=1;                                                                   % number of parameter samples, or size of the ensembles
usenormal=1;
par.do_nob = 1;                                                            % set it to zero to turn off nob
par.do_anam = 1;
par.do_hetero = 1;                                                         % % set it to zero to turn off heterotrophs
dt  = 3600;                                                                % time step (s) 
kend= 720;                                                                 % total run time in hrs (720 hrs = 1 month), 
%t   = 0.0;
par.hist_frq=2;                                                            % Sampling frequency what units?
par.is_stamp=1;                                                            
nhist=kend/par.hist_frq;
par.spinup  = 0;
par.cnbalance_check = 0;
par.restart=0;
kend1=0;
if(par.restart)
    par.restart_inputfile='restart_new2.mat';                              %if it is a restart run, put the correct name of the restart file here
    par.spinup=0;                                                          %no spinup is needed for a restart run
end
par.restart_outputfile ='';                                                %name the output restart file, set it to empty when it is not requested
%%
% Instructions for substrate uptake
% monod1     : classic M-M kinetics
% monod2     : M-M kinetics with diffusion limiting
% haldane    : haldane kinetics
% haldanedifl: haldane kinetics with diffusion limiting
% ecakinetics: See Tang and Riley Biogeosciences 2014
par.kinetics_aob_nh3='haldane';
par.kinetics_anam_nh3='haldane';
par.kinetics_anam_o2='haldane';
par.kinetics_nob_no2='monod2';
par.kinetics_co2='monod2';
par.kinetics_hetero='monod1';
%par.kinetics_aob_nh3='ecakinetics';
%par.kinetics_anam_nh3='ecakinetics';
%par.kinetics_anam_o2='ecakinetics';
%par.kinetics_nob_no2='ecakinetics';
%par.kinetics_co2='ecakinetics';
%par.kinetics_hetero='ecakinetics';

%%
% Stoichiometric ratios of different organic substrates
% NOTE: Names are for convenience, & have little to do with actual
% substrate.
par.NC_DON = 1./5;                                                         
par.NC_AMA = 1./5;
par.NC_UREA = 1./5;
par.NC_glutamate = 1./3;                                               
par.NC_acetate = 1./15;                                               
par.NC_glucose = 1./10;                                               

%%
% Trade-offs implicit in the code. 
% 1. Specific affinity trades-off: A positive relationship is noted between
% KM and VMAX. 
% 2. Mixotrophy trade-off against growth (AOB/ NOB), as does the capacity to take-up
% and utilize multiple substrates (heterotrophs). 
% 3. Specific affinity Vs. growth rate: 
% 4. Growth rate Vs. efficiency. Rate yield trade off (RYTO). Negative relationship. 
%%
% AMMONIA-OXIDIZING BACTERIA/ ARCHAEA
% 4 different guilds are parameterized and represent two nitrosomonas, one
% nitrosospira and an ammonia-oxidizing archaea. See references from
% Martens-Habbena, Prosser, and also  Koops and Pommerening-Roser.
par.km_nh3_aob      = [200, 100, 40, 0.05]*1d-6;                           % NH3 affinity constants (M) Most are in the mid uM range. 'AOA' is in the mid nM range
par.vmax_nh3_aob    = [10, 7, 2, 1.7]*1d-5;                                % max uptake rates (M/s)
par.kmscal_nh3_aob  = [0, 0, 0, 0];                                        % scaling parameter accounting for diffusion barrier
par.ki_nh3_aob      = [300, 300, 80, 1d-3];                                % NH3 inhition factor haldane kinetics.
par.nb_aob          = length(par.km_nh3_aob);                              % number of aob guilds
%Qmin = 6;                                                                 % Cross guild minimum stoichiometric quota for cell growth
%Qmax = 20;                                                                % Cross guild maximum stoichiometric quota for cell growth
par.rcn_aob_min     = [6.6, 6.6, 6.6, 6.6];                                % Guild specific minimal stoichiometry     
par.rcn_aob_max     = [7.5, 7.5, 7.8, 7];                                  % Guild specific maximal stoichiometry
par.kmscal_co2_aob  = [0.0, 0.0, 0.0, 0.0];                                % diffusion limiting parameter
par.km_co2_aob      = [1,1,1,1].*1d-6;                                     % CO2 affinity
par.ki_co2_aob      = par.km_co2_aob.*1d1000;                              % inhibition factor for co2 if haldane kinetics is used, set extremely high to avoid this.
par.km_o2_aob       = [3, 3, 3, 0.3].*1d-6;                                % O2 affinity
par.gmax_cell_aob=[0.9, 0.8, 0.3, 0.25]./3600;                             % Maximum cell division rate - taken from genome informed generation times, divided by hrs.
par.yld_nh3_aob = [0.04, 0.045, 0.05, 0.06]...                             % [0.05, 0.04, 0.075, 0.1]
    ./par.rcn_aob_min;                                                     % yield for biomass synthesis, inside the model, the ratio is flexible 
par.kd_aob      = [0, 0, 0, 0];                                            % First order mortality decay rate [1, 1, 1, 1].*1d-10; 
par.Eact_aob    = [50, 50, 50, 50].*1d4;                                   % Arrhenius equation writes logK=-Ea/(RT)+log(A), the range of Ea is 38.6-87.1 kJ mol ^-1 
par.nh3_aob_porter_dd =[0,0,0,0];                                          % Cell's substarte uptake capacity through diffusion, Bonachela et al., 2011
par.nh3_aob_tmin =[274,274,274,274];                                       % Minimum temperature for enzyme activity, see ref. Ratkowsky 1982, 1983
par.nh3_aob_tmax =[305,305,305,305];                                       % Maximum temperature for enzyme activity
par.nh3_aob_topt=[299,299,299,299];                                        % Optimal temperature for enzyme activity
%%
% N2O production pathways: Ammonia-oxidation detoxification of nitrite through nitrifier
% denitrification & basel rate of NH2OH decomposition 
% There are 3 main assumptions of N2O from nitrifier denitrification
    % Assumption 1: The cost of detoxifying NO2 is borne out from the
    % biomass, resulting in a loss of assimilated biomass.
    % Assumption 2: The NO generated from detoxficiation is not negligible, 
    % particularly in a very low O2 environment, which is two-stage
    % detoxification appraoch  is used
    % Assumption 3: If the available O2 also support a sufficent NOB 
    % population, then the first step of detoxification will
    % be surpressed. This is represented as an inhibition on affinity,
    % with higher no2/o2 ratio, higher affinity, greater value of 
    % sum(x(par0.id_nob_cell).*monod(o2_fr,par0.kme_o2_nob)) lower affinity
    % step1: NO2(-) + 0.25CH2O + H(+) -> NO + 0.25CO2 + 0.75H2O
    % Few studies specifically measure these parameters, therefore, they
    % need to be calibrated before using them. 
% Parameters for NO2(-) toxicity
% par.no2_detox_aob=[10,10,10,10;                                            % Inhibition factor from toxicity, high toxicity, high affinity
%                    1d9,1d9,1d9,1d9]';                                      % Inhibition factor from nob competition for no2, high nob activity, low affinity
% par.no_detox_aob =[10,10,10,10];                                           % Inhibition factor from toxicity for no
% par.km_NO2_aob   =[7, 7, 7, 7].*1d-6;                                        % NO2 affinity during detoxification
% par.km_NO_aob    =[7, 7, 7, 7].*1d-6;                                        % NO affinity during detoxification
% par.no2_tox_aob  =[10,10,10,10];                                           % Toxicity factor - NO2 inhibits NH4 oxidation 
% par.NO_vmax_aob  =[1, 1, 1, 1].*1d-7;                                         % Maximum rate for NO detoxification, set to zero to turn off detoxification
% par.NO2_vmax_aob =[1, 1, 1, 1].*1d-7;                                         % Maximum rate for NO2 detoxification,set to zero to turn off detoxification
% par.hydroxyl_dcmp_aob=[1,1,1,1].*1d-3;                                     % Fraction of hydroxylamine decomposition (No data to constrain this - it might have to be approximated from N2O data)
%
par.no2_detox_aob=[10,10,10,10;                                            % Inhibition factor from toxicity, high toxicity, high affinity
                   1d9,1d9,1d9,1d9]';                                      % Inhibition factor from nob competition for no2, high nob activity, low affinity
par.no_detox_aob =[10,10,10,10];                                           % Inhibition factor from toxicity for no
par.km_NO2_aob   =[7, 7, 7, 7].*1d-60;                                        % NO2 affinity during detoxification
par.km_NO_aob    =[7, 7, 7, 7].*1d-60;                                        % NO affinity during detoxification
par.no2_tox_aob  =[10,10,10,10];                                           % Toxicity factor - NO2 inhibits NH4 oxidation 
par.NO_vmax_aob  =[1, 1, 1, 1].*1d-70;                                         % Maximum rate for NO detoxification, set to zero to turn off detoxification
par.NO2_vmax_aob =[1, 1, 1, 1].*1d-70;                                         % Maximum rate for NO2 detoxification,set to zero to turn off detoxification
par.hydroxyl_dcmp_aob=[1,1,1,1].*1d-30;                                     % Fraction of hydroxylamine decomposition (No data to constrain this - it might have to be approximated from N2O data)
%
%%
% Parameters for AOB/ AOA growing nitrifying organic nitrogen (hydrolysis
%and oxidation)
par.km_UREA_aob   = [1d-6,1d-6,1d-6,1d-6;
    1d-6,1d-6,1d-6,1d-6;
    1d-6,1d-6,1d-6,1d-6]';                                                 % Urea affinty, kmf=km(1)*(1+km(2)*10.^(-km(3)./pH)), increases with pH
par.vmax_UREA_aob = [10, 7, 2, 1.7]*1d-6;                                  % Urea uptake rate.
par.yld_UREA_aob  = [0.05, 0.04, 0.075, 0.1]...
    ./par.rcn_aob_min;                                                     % Yield rate for AOB growing on urea, 1 M UREA leads to N mol NH3 assimilation

par.km_AMA_aob    = [1d-7,1d-7,1d-7,1d-7;
    10,10,10,10;
    1,1,1,1]';                                                             % Amino acid affinty.
par.vmax_AMA_aob  = [0,0,0,0];                                             % Amino acid uptake rate
par.yld_AMA_aob   = [0.02,0.05,0.073,0.1]...
    ./par.rcn_aob_min;                                                     % Yield rate for AOB growing on Amino acids.

par.km_DON_aob    = [1d-7,1d-7,1d-7,1d-7;
    20,20,20,20;
    1,1,1,1]';                                                             % DON (of predetermined CN) affinty.
par.vmax_DON_aob  = [0,0,0,0].*1;                                          % DON uptake rate
par.yld_DON_aob   = [0.02,0.05,0.073,0.1]...
    .*0./par.rcn_aob_min;                                                  % Yield rate for AOB growing on DON.
                                
%%
% NITRITE-OXIDIZING BACTERIA
% NOB - parameter values and guild set up by Xavier Le Roux @ Lyon (Le Roux et al., Frontiers).
% We parameterize three different guilds, an autotrophic and a mixotrophic copiotroph
% and an autotrophic oligotroph.
% NOTE: The fourth guild here was added prior to the heterotrophs below. This 
% served as a heterotrophic analog mineralizing organic matter to release NH3. 
% It is currently turned off
% par.kme_no2_nob=[349, 1035, 19, 1d9].*1d-6;                                % copiotrophic autotroph, copiotrophic mixotroph, oligotrophic autotroph, heterotroph
% par.vmax_no2_nob=[6, 8, 2, 0].*1d-5;                                       % NO2 uptake rate in M/s
% par.kme_o2_nob=[8, 10, 6, 100].*1d-6;                                      % NO2 affinity in M
% par.nb_nob  = length(par.vmax_no2_nob);                                    % Number of NOB guilds
% par.rcn_nob_max=[8, 7, 7.5, 9];                                            % Guild specific maximal stoichiometry
% par.rcn_nob_min=[6.6, 6.6, 6.6, 6.6];                                      % Guild specific maximal stoichiometry
% par.gmax_cell_nob=[0.7, 0.5, 0.3, 0.3]/3600;                               % Guild specific growth rate 
% par.kd_nob =[1, 1, 1, 1].*1d-100;                                          % Mortality rate? 
% par.yld_co2_nob=[0.04, 0.04, 0.05, 0.02]./par.rcn_nob_min;                 % CO2 based biomass yield rate for biomass assimilation = 50, 100, 17 moles of nitrite/ mole of C [0.02,0.01,0.06,0.06], [0.011, 0.015, 0.02, 0.04
% par.kmscal_no2_nob  = [0.,0.,0.,0.];
% par.kmscal_NH3_nob  = [0.,0.,0.,0.];
% par.kmb_NH3_nob  = [2,2,2,1].*1d-6;                                        % NH3 affinity in M. NH3 required to balance cellular stoichio.
% par.km_DON_nob = [1d100,5,1d100,10].*1d-6;                                 % Organic nitrogen affinity in M.
% par.yld_DON_nob=[0.05,0.05,0.05,0.5];                                      % Biomass yield growing on DON.
% par.vmax_DON_nob = [0.0,3.0,0.0, 6.0].*1d-6;                               % Organic nitrogen uptake rate in M/s.
% par.alpahN_DON_nob =[0.0,1.0,0.0,1.0];                                     % Inorganic N assimilation rate from organic nitrogen to maintain cellular CN stoichiometry
% %par.Eact_nob    = [34.2, 34.2]*1d3;                        
% %par.Eact_nob = par.Eact_nob./(8.3144621*273.15);
% % The activation energy for enzyme activity is not used here because it's 
% % difficult to incorporate the asymmetric temperature response curve
% par.no2_nob_porter_dd=[0.,0.,0., 0.];                                      % Cell's substarte uptake capacity through diffusion, idea from Bonachela et al., 2011, 1/(4*pi*D*rc)
% par.no2_nob_tmin =[274,274,274,274];                                       % Minimum temperature for enzyme activity, see ref. Ratkowsky 1982, 1983 
% par.no2_nob_tmax =[305,305,305,305];                                       % Maximum temperature for enzyme activity
% par.no2_nob_topt= [299,299,299,299];                                       % Optimal temperature for enzyme activity
%
% Create Commamox
par.kme_no2_nob=[1, 1, 1, 1].*1d100;                                        % copiotrophic autotroph, copiotrophic mixotroph, oligotrophic autotroph, heterotroph
par.vmax_no2_nob=[1, 1, 1, 1].*1d100;                                       % NO2 uptake rate in M/s
par.kme_o2_nob=[1, 1, 1, 1].*1d100;                                          % NO2 affinity in M
par.nb_nob  = length(par.vmax_no2_nob);                                    % Number of NOB guilds
par.rcn_nob_max=[8, 7, 7.5, 9];                                            % Guild specific maximal stoichiometry
par.rcn_nob_min=[6.6, 6.6, 6.6, 6.6];                                      % Guild specific maximal stoichiometry
par.gmax_cell_nob=[0.7, 0.5, 0.3, 0.3]/3600;                               % Guild specific growth rate 
par.kd_nob =[1, 1, 1, 1].*1d-100;                                          % Mortality rate? 
par.yld_co2_nob=[0.04, 0.04, 0.05, 0.02]./par.rcn_nob_min;                 % CO2 based biomass yield rate for biomass assimilation = 50, 100, 17 moles of nitrite/ mole of C [0.02,0.01,0.06,0.06], [0.011, 0.015, 0.02, 0.04
par.kmscal_no2_nob  = [0.,0.,0.,0.];
par.kmscal_NH3_nob  = [0.,0.,0.,0.];
par.kmb_NH3_nob  = [300,2,2,1].*1d-6;                                        % NH3 affinity in M. NH3 required to balance cellular stoichio.
par.km_DON_nob = [1d100,5,1d100,10].*1d-6;                                 % Organic nitrogen affinity in M.
par.yld_DON_nob=[0.05,0.05,0.05,0.5];                                      % Biomass yield growing on DON.
par.vmax_DON_nob = [0.0,3.0,0.0, 6.0].*1d-6;                               % Organic nitrogen uptake rate in M/s.
par.alpahN_DON_nob =[0.0,1.0,0.0,1.0];                                     % Inorganic N assimilation rate from organic nitrogen to maintain cellular CN stoichiometry
%par.Eact_nob    = [34.2, 34.2]*1d3;                        
%par.Eact_nob = par.Eact_nob./(8.3144621*273.15);
% The activation energy for enzyme activity is not used here because it's 
% difficult to incorporate the asymmetric temperature response curve
par.no2_nob_porter_dd=[0.,0.,0., 0.];                                      % Cell's substarte uptake capacity through diffusion, idea from Bonachela et al., 2011, 1/(4*pi*D*rc)
par.no2_nob_tmin =[274,274,274,274];                                       % Minimum temperature for enzyme activity, see ref. Ratkowsky 1982, 1983 
par.no2_nob_tmax =[305,305,305,305];                                       % Maximum temperature for enzyme activity
par.no2_nob_topt= [299,299,299,299];                                       % Optimal temperature for enzyme activity

%%
% ANAEROBIC AMMONIUM OXIDATION
% Returns the biomass for anammox planctomycetes growing on ammonia and
% nitrite yielding N2. Trait parameters are taken from few sources
% measuring it at the present time. We resolve two different guilds of
% anammox representing Brocadia anammoxidans and Kuenenia stuttgartiensis.
% Kartal et al., 2012: B. sinica has a growth rate ~ 1.5 x higher than the 
% other two, but with lower KM values for substrate. K. suttgartiensis has 
% an affinity approach with a lower growth rate. 
% Represent as clear r and k strat. for anammox. 
% This approach assumes that NH3 goes toward the microbial biomass (rather 
% than equal contribution from NO2). The process is inhibited by nitrite
% and oxygen.
% Guild #1 below: B. sinica; Guild #2: K. stuttgartiensis.
par.km_no2_anam = [80, 0.5].*1d-6;                                         % NO2 affinity constants (M) Range: 86+/- 4; 0.3 - 3.                               
par.vmax_no2_anam = [4, 1].*1d-5;                                         % Max uptake rate of NO2 (M/s) Value: unknown; Range: Unknown
par.kmscal_no2_anam  = [0., 0];                                            % Scaling parameter accounts for diffusion barrier
par.km_nh3_anam = [20, 0.5].*1d-6;                                         % NH3 affinity constants (M) No data for guild #2, based on NO2
par.vmax_nh3_anam = [4, 1].*1d-5;                                          % Max uptake rate of NH3 (M/s)
par.kmscal_NH3_anam  = [0., 0]; 
par.ki_no2_anam = [5, 20].*1d6;                                            % Inhibition constant for NO2 (M) 
par.ki_o2_anam = [1, 1].*1d-9;                                             % O2 inhibition constant (M) Measured at 30 and 200 d-6M -- but very difficult to implement. 
par.km_co2_anam = [1,1].*1d-6;                                             % CO2 affinity constant (M) 
par.nb_anam = length(par.km_no2_anam); 
par.rcn_anam_min = [6.6, 6.6];                                             % Stoichiometry - No data assume similar to nitrifiers                               
par.rcn_anam_max = [7, 7];  
par.kmscal_co2_anam  = [0., 0]; 
par.gmax_cell_anam = [0.008, 0.004]/3600;                                  % Growth rates (d-1) basically ~ 2d-6 & 1 d-6
par.kd_anam = [1,1].*1d-12;                            
par.yld_co2_anam = [0.01,0.01]./par.rcn_anam_min;                          % Yield coefficient - moles of NO2/ mole CO2 fixed, measured for guild #1 not #2. 
par.kmb_NH3_anam  = [2,2].*1d-6;
par.no2_anam_porter_dd = [0,0];
par.no2_anam_tmin = [275,277];                                        
par.no2_anam_tmax = [307,308];                                       
par.no2_anam_topt = [293,296];  
par.nh3_anam_porter_dd = [0,0]; 
par.nh3_anam_tmin = [274,274];                                        
par.nh3_anam_tmax = [305,305];                                       
par.nh3_anam_topt = [299,299];  
%%
% AEROBIC AND ANAEROBIC HETEROTROPHS
% Heterotrophs -- 12 different groups -- broken down into 4 groups of three ecological strategies.
% Four groups: (1) Obligate aerobes (use for now, but these probably don't exist?), (2) Facultative canonical denitrifiers (NO3 -> N2), 
% (3) Facultative partial denitrifiers, and (4) Aerobic N2O consumers (see Sanford et al., PNAS, 2012 or Jones et al., Nat. Climate Change, 2013).
% These groups account for competition, N2O production and consumption.
% Ecological strategies - Specialist (one carbon compound, i.e., just glutamate), Intermediate (2 carbon compounds, glutamate and glucose), 
% and generalist (three carbon compounds, glutamate). Increased metabolic diversity trades-off against growth rate  
% See note on carbon compund names above.
% Organized below by group [1,1,1,2,2,2,3,3,3,4,4,4] == [specialist, intermediate, generalist] 
% These reactions are stoichiometrically balanced (see biology.m ln 206
% onwards)
%
%par.vmax_no3_hetero=[3, 5, 1, 0, 0].*1d-6;
% par.km_no3_hetero=[100, 100, 100, 10, 10, 10,...
%     10, 10, 10, 10, 10, 10].*1d-6;                                      % NO3 affinity constant for denitrifiers (M)
par.km_no3_hetero=[10, 10, 10, 1, 1, 1,...
    1, 1, 1, 10, 10, 10].*1d-6;
%par.vmax_acetate_hetero=[1d-100, 1d-100, 1, 1d-100, 1d-100, 1d-100,...
%    1d-100, 1d-100, 1, 1d-100, 1d-100, 1].*1d-6;
% Acetate uptake
par.km_acetate_hetero=[1d9, 1d9, 100, 1d9, 1d9, 100,... 
    1d9, 1d9, 100, 1d9, 1d9, 100].*1d-6;                                   % OM compound #1 affinity (M)                         
par.vmax_acetate_hetero=[0,0,10,0,0,10,...
    0,0,10,1d-100,1d-100,1d-100].*1d-6;                                                  % OM compound #1 uptake rate (M/s) all guilds aerobic/ partial denitrification pathway                   
par.vmax_acetatec_hetero=[1d-100, 1d-100, 1d-100, 1d-100, 1d-100, 10,...
    1d-100, 1d-100, 1d-100, 1d-100, 1d-100, 1d-100].*1d-6;                 % OM compound #1 uptake rate (M/s) complete denitrification pathway                       
% Glucose uptake
par.km_glucose_hetero=[1d9, 50, 100, 1d9, 50, 100,... 
    1d9, 50, 100, 1d9, 50, 100].*1d-6;                                     % OM compound #2 affinity (M)
par.vmax_glucose_hetero=[1d-100, 10, 5, 1d-100, 10, 5,...
    1d-100, 10, 5, 1d-100, 1d-100, 1d-100].*1d-6;                                    % OM compound #2 uptake rate (M/s) all guilds aerobic/ partial denitrification pathway 
par.vmax_glucosec_hetero=[1d-100, 1d-100, 1d-100, 1d-100, 50, 1,...
    1d-100, 1d-100, 1d-100, 1d-100, 1d-100, 1d-100].*1d-6;                 % OM compound #2 uptake rate (M/s) complete denitrification pathway 
% Glutamate uptakte
par.km_glutamate_hetero=[50, 250, 500, 50, 250, 500,...
    50, 250, 500, 50, 500, 1000].*1d-6;                                    % OM compound #3 affinity (M)
par.vmax_glutamate_hetero=[28, 15, 8, 20, 10, 5,...
    21, 12, 7, 1d-100, 1d-100, 1d-100].*1d-6;                                             % OM compound #3 uptake rate (M/s) all guilds aerobic/ partial denitrification pathway
par.vmax_glutamatec_hetero=[1d-100, 1d-100, 1d-100, 10, 5, 1,...
    1d-100, 1d-100, 1d-100, 1d-100, 1d-100, 1d-100].*1d-6;                 % OM compound #3 uptake rate (M/s) complete denitrification pathway 
%par.vmax_glutamate_hetero=[30, 15, 10, 10, 5, 1,...
%    10, 5, 1, 10, 5, 2].*1d-6;
%par.vmax_glutamatec_hetero=[1d-100, 1d-100, 1d-100, 15, 8, 4,...
%    1d-100, 1d-100, 1d-100, 1d-100, 1d-100, 1d-100].*1d-6;
par.km_n2o_hetero=[1d9, 1d9, 1d9, 1d9, 1d9, 1d9,...
    1d9, 1d9, 1d9, 10, 50, 100].*1d-6;                                     % N2O affinity constant (M)                      

%par.vmax_o2_hetero=[3, 5, 1, 0, 0].*1d-6;
par.km_o2_hetero=[0.5, 0.5, 0.5, 5, 5, 5, 5, 5, 5, 1d9, 1d9, 1d9].*1d-6;   % O2 affinity constant (M)
par.nb_hetero = length(par.vmax_glucose_hetero);                           % Number of heterotrophic guilds
par.rcn_hetero_max=[9,9,9,9,9,9,9,9,9,9,9,9];                              % Guild specific maximal stoichiometry                            
par.rcn_hetero_min=[5,5,5,5,5,5,5,5,5,5,5,5];                              % Guild specific maximal stoichiometry      
par.gmax_cell_hetero=[8, 5.5, 2.5, 6, 4, 1,... 
    7, 5, 2, 4, 0.8, 0.3]/3600;                                             % Guild specific growth rate
par.kd_hetero =[4d-9, 4d-9, 4d-9, 4d-9, 4d-9, 4d-9,...
    4d-9, 4d-9, 4d-9, 4d-9, 4d-9, 4d-9];  
%par.kd_hetero =[9.9d-8, 9.9d-8, 9.9d-8, 9.9d-8, 9.9d-8, 9.9d-8,...
%    9.9d-8, 9.9d-8, 9.9d-8, 9.9d-8, 9.9d-8, 9.9d-8];                      % Guild specific first order mortality rate.

% Below are the biomass yield factors for organisms growing with the
% different carbon sources and electron acceptor pathways (O2, NO3->N2O or (seperately) NO3->N2). 
par.yld_cglucose_hetero=[0, 0.025, 0.035, 0, 0.035, 0.045,...
   0, 0.03, 0.04, 0, 0.03, 0.04]./par.rcn_hetero_min;                      
%par.yld_cglucosec_hetero=[0, 0, 0, 0, 0.25, 0.25,...
%    0, 0, 0, 0, 0, 0]./par.rcn_hetero_min;
par.yld_cacetate_hetero=[0, 0, 0.035, 0, 0, 0.05,...
   0, 0, 0.04, 0, 0, 0.05]./par.rcn_hetero_min;                             
par.yld_cacetatec_hetero=[0, 0, 0, 0, 0, 0.065,...
    0, 0, 0, 0, 0, 0]./par.rcn_hetero_min;                                  
par.yld_cglutamate_hetero=[0.035, 0.045, 0.055, 0.045, 0.055, 0.065,...
   0.042, 0.052, 0.062, 0.04, 0.05, 0.06]./par.rcn_hetero_min;
par.yld_cglutamatec_hetero=[0, 0, 0, 0.045, 0.055, 0.065,...
    0, 0, 0, 0, 0, 0]./par.rcn_hetero_min;    

par.km_nh3_hetero = [100, 100, 100, 100, 100, 100,...
    100, 100, 100, 100, 100, 100]*1d-6;                                    % NH3 affinity (M) to maintain cellular stoichiometry 
%par.vmax_nh3_hetero = [1, 1, 1, 1, 1, 1, 1, 1, 1,...
%    1, 1, 1]*1d-6;
par.ki_o2_hetero = [1, 1, 1, 0.1, 0.1, 0.1,...
    0.1, 0.1, 0.1, 1, 1, 1]*1d-6;                                          % Inhibition of NO3-reduction by oxygen concentrations - This can be in the nM range (Stolper et al., PNAS).

par.yld_N_hetero=[0.03, 0.03, 0.03, 0.03, 0.03, 0.03,...
    0.03, 0.03, 0.03, 0.03, 0.03, 0.03];
par.alpahN_N_hetero =[0.5, 0.5, 0.5, 0.5, 0.5, 0.5,...
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5];                                         % Nitrogen assimilation rate from a compounds N to maintain CN stoichiometry
par.alpahN_glutamate_hetero =[0.5, 0.5, 0.5, 0.5, 0.5, 0.5,...
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5];                                         % Nitrogen assimilation rate from a compounds N to maintain CN stoichiometry
par.alpahN_glucose_hetero =[0.5, 0.5, 0.5, 0.5, 0.5, 0.5,...
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5];                                         % Nitrogen assimilation rate from a compounds N to maintain CN stoichiometry
par.alpahN_acetate_hetero =[0.5, 0.5, 0.5, 0.5, 0.5, 0.5,...
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5];                                         % Nitrogen assimilation rate from a compounds N to maintain CN stoichiometry


% Below is information about transporter density and enzyme temperature
% responses. Similar to that described above, except this takes into
% account the membrane enzymes involved in carbon compound and nutrient uptake. 

par.nh3_hetero_porter_dd=[0.,0.,0., 0., 0., 0.,0.,0., 0., 0., 0., 0.];     % % Cell's substarte uptake capacity through diffusion, idea from Bonachela et al., 2011, 1/(4*pi*D*rc)   
par.nh3_hetero_tmin =[274,274,274,274,274,274,...
    274,274,274,274,274,274];                                              % Minimum temperature for enzyme activity, see ref. Ratkowsky 1982, 1983        
par.nh3_hetero_tmax =[305,305,305,305,305,305,...
    305,305,305,305,305,305];                                              % Maximum temperature for enzyme activity.               
par.nh3_hetero_topt= [299,299,299,299,299,299,...
    299,299,299,299,299,299];                                              % Optimal temperature for enzyme activity

par.no3_hetero_porter_dd=[0.,0.,0., 0., 0., 0.,0.,0., 0., 0., 0., 0.];                                             
par.no3_hetero_tmin =[274,274,274,274,274,274,...
    274,274,274,274,274,274];                                               
par.no3_hetero_tmax =[305,305,305,305,305,305,...
    305,305,305,305,305,305];                                              
par.no3_hetero_topt= [299,299,299,299,299,299,...
    299,299,299,299,299,299];

par.o2_hetero_porter_dd=[0.,0.,0., 0., 0., 0.,0.,0., 0., 0., 0., 0.];                                             
par.o2_hetero_tmin =[274,274,274,274,274,274,...
    274,274,274,274,274,274];      
par.o2_hetero_tmax =[305,305,305,305,305,305,...
    305,305,305,305,305,305];                                               
par.o2_hetero_topt= [299,299,299,299,299,299,...
    299,299,299,299,299,299];

par.glucose_hetero_porter_dd=[0.,0.,0., 0., 0., 0.,0.,0., 0., 0., 0., 0.];                                             
par.glucose_hetero_tmin =[274,274,274,274,274,274,...
    274,274,274,274,274,274];                                               
par.glucose_hetero_tmax =[305,305,305,305,305,305,...
    305,305,305,305,305,305];                                               
par.glucose_hetero_topt= [299,299,299,299,299,299,...
    299,299,299,299,299,299];

par.acetate_hetero_porter_dd=[0.,0.,0., 0., 0., 0.,0.,0., 0., 0., 0., 0.];                                             
par.acetate_hetero_tmin =[274,274,274,274,274,274,...
    274,274,274,274,274,274];                                               
par.acetate_hetero_tmax =[305,305,305,305,305,305,...
    305,305,305,305,305,305];                                               
par.acetate_hetero_topt= [299,299,299,299,299,299,...
    299,299,299,299,299,299];

par.glutamate_hetero_porter_dd=[0.,0.,0., 0., 0.];                          
par.glutamate_hetero_tmin =[274,274,274,274,274,274,...
    274,274,274,274,274,274];                                               
par.glutamate_hetero_tmax =[305,305,305,305,305,305,...
    305,305,305,305,305,305];                                             
par.glutamate_hetero_topt= [299,299,299,299,299,299,...
    299,299,299,299,299,299];

%%
%set up id of state variables
par.id_nh3x=1;                                                             %nh4(+), nh3(g), nh4(s), nh3(aq)
par.id_no2x=2;                                                             %no2(-), hno3(aq)
par.id_nox =3;                                                             %no(g), no(aq)
par.id_n2ox=4;
par.id_o2x =5;                                                             %o2(g), o2(aq)
par.id_no3x=6;
par.id_mcb_rsc=7;                                                          %c pool of microbial residues
par.id_mcb_rsn=8;                                                          %n pool of microbial residues
par.id_don_c = 9;                                                          %c pool of don
par.id_don_n = 10;                                                         %n pool of don
par.id_UREA = 11;
par.id_AMA  = 12;
par.id_DON  = 13;
par.id_DOC  = 14;
par.id_co2x  =15;                                                          %co2(aq),hco3(-),co3(2-),h2co3
par.id_glucose = 16;
par.id_acetate = 17;
par.id_glutamate = 18;

%stand free energy of formation
par.freeEform = zeros(par.id_glutamate);

%charge of the ion
par.charge = zeros(par.id_glutamate);
par.charge(par.id_no2x) = -1;
par.charge(par.id_no3x) = -1;
par.charge(par.id_glutamate) = -1;

%alpha - hydrated ion radius (nm)
par.ionradius = zeros(par.id_glutamate);
par.ionradius(par.id_no2x) = 0.3;
par.ionradius(par.id_no3x) = 0.3;
par.ionradius(par.id_acetate) = 0.45;
par.ionradius(par.id_glutamate) = 0.45;                                     % NOTE: Assumed. I'm not sure what  the ionic radius of this for now

par.activitycoeff = ones(par.id_glutamate);


%%
for j = 1 : par.nb_aob
    par.id_aob_c(j)=par.id_glutamate+3*j-2;                                %mol C
    par.id_aob_n(j)=par.id_aob_c(j)+1;                                     %mol N
    par.id_aob_cell(j)=par.id_aob_n(j)+1;                                  %mol C equivalent
end
xsiz=par.id_aob_cell(par.nb_aob);

if(par.do_nob)
    for j = 1 : par.nb_nob    
        par.id_nob_c(j)=xsiz+3*j-2;    
        par.id_nob_n(j)=par.id_nob_c(j)+1;    
        par.id_nob_cell(j)=par.id_nob_n(j)+1;
    end    
    xsiz=par.id_nob_cell(j);
end

if(par.do_anam)
for j = 1 : par.nb_anam
    par.id_anam_c(j)=xsiz+3*j-2;                                %mol C
    par.id_anam_n(j)=par.id_anam_c(j)+1;                                     %mol N
    par.id_anam_cell(j)=par.id_anam_n(j)+1;                                  %mol C equivalent
end
xsiz=par.id_anam_cell(j);
end

if(par.do_hetero)
    for j = 1 : par.nb_hetero    
        par.id_hetero_c(j)=xsiz+3*j-2;    
        par.id_hetero_n(j)=par.id_hetero_c(j)+1;    
        par.id_hetero_cell(j)=par.id_hetero_n(j)+1;
    end    
    xsiz=par.id_hetero_cell(j);
end
%Initialize N2O from either nitrification or denitrification
par.id_nitrificationflux = xsiz + 1;
par.id_denitrificationflux = par.id_nitrificationflux+1;
xsiz = par.id_denitrificationflux;

%solubility control
par.ke_nh4_nh3=10^(9.24);                                                  %nh4(+) -> nh3(aq)+H(+),[nh4(+)]=ke*[nh3(aq)]*[h(+)]

%%
%Magnitude and timestep for pulsed additions
% See spbiology.m for structure. 

% NH4 pulses
par.nh4_add_amount =  []; 
par.nh4_add_time   = []; 

% NO2 pulse
par.no2_add_amount =  []; 
par.no2_add_time   = [];  

% NO3 pulse
par.no3_add_amount =  []; 
par.no3_add_time   = [];

% O2 pulse
par.o2_add_amount = [1d-2];                                                    % 1d-3, 1d-3, 5d-3, 5d-3, 1d-2 2d-3, 1d-3, 1d-3, 1d-3, 1d-3  5d-3, 5d-3, 5d-3, 5d-3
par.o2_add_time  = [2000];                                                     % 1650, 1700, 1750, 1800, 1850 650, 1600, 2000, 2200, 2400, 2600, 2800, 3180, 3540
        
% DON pulse
par.DON_add_amount = [];
par.DON_add_time  = [];

% ED_1: Glutamate pulse
par.glutamate_add_amount = [];                                             % [2d-5, 2d-5, 2d-5, 2d-5, 2d-5, 2d-5,...
                                                                           %  2d-5, 2d-5, 2d-5, 2d-5, 2d-5, 2d-5];
                                                                           % [5d-5, 5d-5, 5d-5, 5d-5, 5d-5, 5d-5,...
                                                                           % 5d-5, 5d-5, 5d-5, 5d-5, 5d-5, 5d-5]; 
                                                                           % [1d-4, 1d-4, 1d-4, 1d-4,...
                                                                           % 1d-4, 1d-4];
par.glutamate_add_time = [];                                               % [360, 720, 1080, 1440, 1800, 2160. 2520];                                                
                                                                           % [720, 1440, 2160, 2880, 3600, 4320]; 
                                                                           % [360, 720, 1080, 1440, 1800, 2160,...
                                                                           % 2520, 2880, 3240, 3600, 3960, 4320]; 
                                                                           % [360,720,1020,1440,...
                                                                           % 1800, 2340];                                     
                                                                           % 180, 360, 540, 720, 900, 1020, 1260, 1440,...
                                                                           % 1620, 1800, 1980, 2340,2880, 4500, 4680, 5040

% ED_2: Glucose pulse
par.glucose_add_amount = [2d-5, 2d-5, 2d-5, 2d-5, 2d-5, 2d-5,...
                           2d-5, 2d-5, 2d-5, 2d-5, 2d-5, 2d-5];
                                                                           % [5d-5, 5d-5, 5d-5, 5d-5, 5d-5, 5d-5,...
                                                                           % 5d-5, 5d-5, 5d-5, 5d-5, 5d-5, 5d-5];                                     
                                                                           % [5d-6, 5d-6, 5d-6, 5d-6, 5d-6, 5d-6, 5d-6, 5d-6,...
                                                                           % 5d-6, 5d-6, 5d-6, 5d-6, 5d-6, 5d-6, 5d-6, 5d-6]; Pulsing of organic matter into the aquifer [5d-4, 5d-4, 5d-4, 5d-4, 5d-4, 5d-4]
% Added every two weeks.
par.glucose_add_time = [360, 720, 1080, 1440, 1800, 2160,...
                        2520, 2880, 3240, 3600, 3960, 4320];                % [360, 720, 1080, 1440, 1800, 2160. 2520]
                                                                           % [360, 720, 1080, 1440, 1800, 2160,...
                                                                           % 2520, 2880, 3240, 3600, 3960, 4320];   
                                                                           % 2880, 4500, 4680, 5040  10 uMs every 2 weeks pulsed (wk 2, 4, 6, 8, 10, 12)

% ED_3: Acetate pulse
par.acetate_add_amount = [2d-5, 2d-5, 2d-5, 2d-5, 2d-5, 2d-5,...
                           2d-5, 2d-5, 2d-5, 2d-5, 2d-5, 2d-5];
                                                                           % [1d-5, 1d-5, 1d-5, 1d-5, 1d-5, 1d-5,...
                                                                           % 1d-5, 1d-5, 1d-5, 1d-5, 1d-5, 1d-5];
par.acetate_add_time  = [360, 720, 1080, 1440, 1800, 2160,...
                        2520, 2880, 3240, 3600, 3960, 4320];                                               % [360, 720, 1080, 1440, 1800, 2160,...
                                                                           % 2520, 2880, 3240, 3600, 3960, 4320];  
                                                                           % 2880, 4500, 4680, 5040
%%
%physical forcing setup, pH, moisture and temperature
env.h2ovol=0.03+zeros(kend+1);                                             %water filled pores
env.airvol=0.1+zeros(kend+1);                                              %air filled pores
env.pH = 7.8+zeros(kend+1);                                                %pH values Don't use above ~ 8 the ammonimu: ammonia ratio is not realistic.

% Specific for the 2d code. NOTE: Flow and transport between grids is
% waiting to be incorporated.
x0=zeros(nsamp,xsiz);

if(par.restart)
    load(par.restart_inputfile);
    if(size(y_rst,1) ~= size(x0,1) || size(y_rst,2)~=size(x0,2)|| size(y_rst,3)~=size(x0,3)|| size(y_rst,4)~=size(x0,4))
        error('bad restart input file specified');
    end
    kend=kend+kend1;
    env.tsoi=297+temperature_generator(dt.*(kend1+1:kend+1),10);           %temperature
    x0=y_rst;
else
end
    env.tsoi=297+temperature_generator(dt.*(1:kend+1),10);                 %temperature    
    %env.bgflx_nh3=nh3_mineralization(dt.*(1:kend+1),1.d-16);              %background mineralization, put a time series to replace it if needed 
    %env.bgflx_no2=no2_mineralization(dt.*(1:kend+1),1.d-16);
    %env.bgflx_o2 = 1+o2_mineralization(dt.*(1:kend+1),1d-1);
    %env.bgflx_don=don_mineralization(dt.*(1:kend+1),1.d-16);
    
    % Set-up the initial conditions
    x0(1,par.id_nh3x)=2d-6;                                                % M NH4 targetvar+zeros(kend+1,1)
    x0(1,par.id_co2x)=1;                                                   % M CO2 - assumes non-limited by CO2
    x0(1,par.id_o2x) =5d-5;                                                % M O2 - assumes non-limited by O2
    x0(1,par.id_no2x)=1d-3;                                                % NO2 - too much with inhibit AOB growth.
    x0(1,par.id_no3x)=5d-3;                                                %
    x0(1,par.id_nox) =0;
    x0(1,par.id_n2ox)=0;
    x0(1,par.id_glucose)= 1.5d-4;
    x0(1,par.id_acetate)= 1.5d-4;
    x0(1,par.id_glutamate)= 0;    

    if(par.id_UREA>0)
        x0(1,par.id_UREA)=0;
    end
    if(par.id_AMA>0)    
        x0(1,par.id_AMA) =0;
    end
    if(par.id_DON>0)    
        x0(1,par.id_DON) =0;
    end
    
    %% 
    % Set individual biomasses based on the spin-up data.
  % AOB biomass
    for j = 1 : par.nb_aob
        x0(1,par.id_aob_c(j))=1d-5;      %mol C
        x0(1,par.id_aob_n(j))=1d-5/6.6;  %mol N
        x0(1,par.id_aob_cell(j))=1d-5;   %mol C
    end
    
    
   %ANAMMOX Biomass
    if(par.do_anam)
        for j = 1 : par.nb_anam
          x0(1,par.id_anam_c(j))=1d-6;      %mol C
          x0(1,par.id_anam_n(j))=1d-6/6.6;  %mol N
          x0(1,par.id_anam_cell(j))=1d-6;   %mol C
        end
    end
    
    % NOB Biomass
    if(par.do_nob)
        for j = 1 : par.nb_nob    
             x0(1,par.id_nob_c(j))=1d-5;      %mol C    
             x0(1,par.id_nob_n(j))=1d-5/par.rcn_nob_min(j);  %mol N    
             x0(1,par.id_nob_cell(j))=1d-5;   %mol C
         end
    end
%     
    % Heterotrophic biomasses
     if(par.do_hetero)
%         for j = 1 : par.nb_hetero    
%             x0(1,par.id_hetero_c(j))=1d-6;      %mol C    
%             x0(1,par.id_hetero_n(j))=1d-6/par.rcn_hetero_min(j);  %mol N    
%             x0(1,par.id_hetero_cell(j))=1d-6;   %mol C            
%          end 
    % Hetero 1    
            x0(1,par.id_hetero_c(1))=1d-6;      %mol C    
            x0(1,par.id_hetero_n(1))=1d-6/par.rcn_hetero_min(1);  %mol N    
            x0(1,par.id_hetero_cell(1))=1d-6;   %mol C            
    % Hetero 2    
            x0(1,par.id_hetero_c(2))=1.15d-6;      %mol C    
            x0(1,par.id_hetero_n(2))=1.15d-6/par.rcn_hetero_min(2);  %mol N    
            x0(1,par.id_hetero_cell(2))=1.15d-6;   %mol C 
    % Hetero 3    
            x0(1,par.id_hetero_c(3))=1.42d-6;      %mol C    
            x0(1,par.id_hetero_n(3))=1.42d-6/par.rcn_hetero_min(3);  %mol N    
            x0(1,par.id_hetero_cell(3))=1.42d-6;   %mol C  
    % Hetero 4    
            x0(1,par.id_hetero_c(4))=1d-6;      %mol C    
            x0(1,par.id_hetero_n(4))=1d-6/par.rcn_hetero_min(4);  %mol N    
            x0(1,par.id_hetero_cell(4))=1d-6;   %mol C 
    % Hetero 5    
            x0(1,par.id_hetero_c(5))=1.2d-6;      %mol C    
            x0(1,par.id_hetero_n(5))=1.2d-6/par.rcn_hetero_min(5);  %mol N    
            x0(1,par.id_hetero_cell(5))=1.2d-6;   %mol C 
    % Hetero 6    
            x0(1,par.id_hetero_c(6))=1.55d-6;      %mol C    
            x0(1,par.id_hetero_n(6))=1.55d-6/par.rcn_hetero_min(6);  %mol N    
            x0(1,par.id_hetero_cell(6))=1.55d-6;   %mol C  
    % Hetero 7    
            x0(1,par.id_hetero_c(7))=1d-6;      %mol C    
            x0(1,par.id_hetero_n(7))=1d-6/par.rcn_hetero_min(7);  %mol N    
            x0(1,par.id_hetero_cell(7))=1d-6;   %mol C  
    % Hetero 8   
            x0(1,par.id_hetero_c(8))=1.13d-6;      %mol C    
            x0(1,par.id_hetero_n(8))=1.13d-6/par.rcn_hetero_min(8);  %mol N    
            x0(1,par.id_hetero_cell(8))=1.13d-6;   %mol C
    % Hetero 9    
            x0(1,par.id_hetero_c(9))=1.46d-6;      %mol C    
            x0(1,par.id_hetero_n(9))=1.46d-6/par.rcn_hetero_min(9);  %mol N    
            x0(1,par.id_hetero_cell(9))=1.46d-6;   %mol C  
    % Hetero 10    
            x0(1,par.id_hetero_c(10))=1.1d-6;      %mol C    
            x0(1,par.id_hetero_n(10))=1.1d-6/par.rcn_hetero_min(10);  %mol N    
            x0(1,par.id_hetero_cell(10))=1.1d-6;   %mol C  
    % Hetero 11    
            x0(1,par.id_hetero_c(11))=1.1d-6;      %mol C    
            x0(1,par.id_hetero_n(11))=1.1d-6/par.rcn_hetero_min(11);  %mol N    
            x0(1,par.id_hetero_cell(11))=1.1d-6;   %mol C  
    % Hetero 12    
            x0(1,par.id_hetero_c(12))=1.1d-6;      %mol C    
            x0(1,par.id_hetero_n(12))=1.1d-6/par.rcn_hetero_min(12);  %mol N    
            x0(1,par.id_hetero_cell(12))=1.1d-6;   %mol C  
     end
    
    for ns = 2 : nsamp
        x0(ns,:)=x0(1,:);
    end


x0(1,par.id_nitrificationflux)=0.0;
x0(1,par.id_denitrificationflux)=0.0;

y=zeros(nhist,xsiz,nsamp);
y_rst=zeros(nsamp,xsiz);
diag.nh4  =zeros(nhist,1);
diag.nh3aq=zeros(nhist,1);
if(nsamp>1)
    par(2:nsamp)=par;
    diag(2:nsamp)=diag;
end

%%
% Below we consider the liklihood of genotypic diversity within guilds.
% nsamp (ln 16) represents the number of individuals/ guild. If that is set
% greater than 1 the trait values will be taken from below and selected by
% MCMC. 
% For the current code we consider the variance in trait space to be equal
% to 2 std dev. of the mean trait value via a unimodel distribution.
% This is for AOB and NOB.

%% AOB

if(nsamp>1)    
    for ns = 2 : nsamp    
        par(ns)=par(1);
    end

    % Vmax
    lb=[9, 6.3, 1.8, 1.6];
    ub=[11, 7.7, 2.2, 1.9];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
         par(ns).vmax_nh3_aob([1,2,3,4])=s(ns,:).*1d-6;
    end

    % Km(NH3)
    lb=[160, 90, 36, 0.04];
    ub=[240, 110, 44, 0.06];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).km_nh3_aob([1,2,3,4])=s(ns,:).*1d-6;
    end

    % Growth rate
    lb=[0.8, 0.72, 0.27, 0.22];
    ub=[0.99, 0.89, 0.33, 0.28];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));        
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).gmax_cell_aob([1,2,3,4]) = s(ns,:)./3600;
    end

    % Km(O2)
    lb=[2.7, 2.7, 2.7, 0.27];
    ub=[3.3, 3.3, 3.3, 0.33];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));        
    else        
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).km_o2_aob([1,2,3,4]) =s(ns,:).*1d-6;
    end

    % SUE
    lb=[0.036, 0.04, 0.045, 0.054];
    ub=[0.044, 0.05, 0.055, 0.066];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));        
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).yld_nh3_aob([1,2,3,4]) = s(ns,:)./par(ns).rcn_aob_min([1,2,3,4]);
    end
    
    
    %NOB
    
    %Vmax(NO2)
    lb=[5.4, 7.2, 1.8, 0];
    ub=[6.6, 8.8, 2.2, 0];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));        
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).vmax_no2_nob([1,2,3,4]) = s(ns,:).*1d-6;
    end

    %Km(NO2)
    lb=[314, 932, 1.8, 1d9];
    ub=[384, 1138, 2.2, 1d9];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));        
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).kme_no2_nob([1,2,3,4]) = s(ns,:).*1d-6;
    end

    %Growth rate
    lb=[0.63, 0.45, 0.27, 0];
    ub=[0.77, 0.55, 0.33, 0];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).gmax_cell_nob([1,2,3,4]) = s(ns,:).*1d-6;
    end

    %O2_affinity
    lb=[7.2,  9, 5.4, 0];
    ub=[8.8, 11, 6.6, 0];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).vmax_DON_nob([1,2,3,4]) = s(ns,:).*1d-6;
    end

    % SUE
    lb=[0.036, 0.036, 0.045, 0];
    ub=[0.044, 0.044, 0.055, 0];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));        
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).yld_co2_nob([1,2,3,4]) = s(ns,:)./par(ns).rcn_nob_min([1,2,3,4]);
    end
    
    
    % Heterotroph

    % KM(NH3)
    lb=[90, 90, 90, 9, 9, 9, 9, 9, 9, 9, 9, 9];
    ub=[110, 110, 110, 11, 11, 11, 11, 11, 11, 11, 11, 11];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).km_nh3_hetero([1,2,3,4,5,6,7,8,9,10,11,12])=s(ns,:).*1d-6;
    end

    % KM(O2)
    lb=[0.45, 0.45, 0.45, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 1d-9, 1d-9, 1d-9];
    ub=[0.55, 0.55, 0.55, 5.5, 5.5, 5.5, 5.5, 5.5, 5.5, 1d-9, 1d-9, 1d-9];
        if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).km_o2_hetero([1,2,3,4,5,6,7,8,9,10,11,12])=s(ns,:).*1d-6;
    end

    % SUE(Glutamate)
    lb=[0.0315, 0.0405, 0.0495, 0.0405, 0.0495, 0.0585, 0.0378, 0.0468, 0.0558, 0.0378, 0.0468, 0.0558];
    ub=[0.0385, 0.0495, 0.0605, 0.0495, 0.0605, 0.0715, 0.0462, 0.0572, 0.0682, 0.0462, 0.0572, 0.0682];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).yld_cglutamate_hetero([1,2,3,4,5,6,7,8,9,10,11,12])=s(ns,:)./par(ns).rcn_hetero_min([1,2,3,4,5,6,7,8,9,10,11,12]);
    end

    % Growth rate
    lb=[7.2, 7.95, 2.25, 5.4, 3.6, 0.9, 6.3, 4.5, 1.8, 6.3, 4.5, 1.8];
    ub=[8.8, 6.05, 2.75, 6.6, 4.4, 1.1, 7.7, 5.5, 2.2, 7.7, 5.5, 2.2];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).gmax_cell_hetero([1,2,3,4,5,6,7,8,9,10,11,12])=s(ns,:)./3600;
    end

    % Ki(O2)
    lb=[0.9, 0.9, 0.9, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09, 0.09];
    ub=[1.1, 1.1, 1.1, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011, 0.011];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).ki_o2_hetero([1,2,3,4,5,6,7,8,9,10,11,12])=s(ns,:).*1d-8;
    end
    
    
    % Anammox

    %KM(NH3)
    lb=[18, 0.45];
    ub=[22, 0.55];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).km_nh3_anam([1,2]) = s(ns,:).*1d-6;
    end

    %Vmax(NH3)
    lb=[3.6, 0.9];
    ub=[4.4, 1.1];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).vmax_nh3_anam([1,2]) = s(ns,:).*1d-5;
    end

    %KM(NO2)
    lb=[72, 0.45];
    ub=[88, 0.55];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).km_no2_anam([1,2]) = s(ns,:).*1d-6;
    end

    %Vmax(NH3)
    lb=[3.6, 0.9];
    ub=[4.4, 1.1];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).vmax_no2_anam([1,2]) = s(ns,:).*1d-5;
    end

    %Ki(O2)
    lb=[1d-7, 1d-7];
    ub=[1d-9, 1d-9];
    if(usenormal)
        mm=(lb+ub)./2;
        sgm=(ub-lb)./6;
        s=latin_hs(mm,sgm,nsamp,length(lb));
    else
        s=lhsu(lb,ub,nsamp);
    end
    for ns = 1 : nsamp
        par(ns).ki_o2_anam([1,2]) = s(ns,:);
    end
end

%compute the shape parameter for temperature scaling
%for ns = 1 : nsamp
%   for j = 1 : par(ns).nb_aob
%       par(ns).nh3_aob_cc(j)=calc(par(ns).nh3_aob_tmin(j),...
%          par(ns).nh3_aob_tmax(j),par(ns).nh3_aob_topt(j));
%   end
%   for j = 1 : par(ns).nb_nob
%       par(ns).no2_nob_cc(j)=calc(par(ns).no2_nob_tmin(j),...
%          par(ns).no2_nob_tmax(j),par(ns).no2_nob_topt(j));
%   end        
%end

        for ns = 1 : nsamp
            [y(:,:,ns),y_rst(ns,:),time]=runspBiology(x0(ns,:), dt, kend1,kend,ns);        
        end


%write restart file if requested
if(~isempty(strfind(par(1).restart_outputfile,'.mat')))
    kend1=kend;
    save(par(1).restart_outputfile,'y_rst','kend1','s');
    
end
telap = cputime-tstart;

time=time./86400;

%%
% Output - Figures and csv files.


figure; 
subplot(3,1,1);
    plot(time,y(:,par(1).id_co2x,ns),'r'); title(['CO_2'])
subplot(3,1,2);
    plot(time,y(:,par(1).id_nh3x,ns),'g'); title(['NH_3'])
subplot(3,1,3);
    plot(time,y(:,par(1).id_o2x,ns),'b'); title(['O_2'])


% figure;
% subplot(4,1,1);
%     plot(time,y(:,par(1).id_no2x,ns),'m'); title(['NO_2'])
% subplot(4,1,2);
%     plot(time,y(:,par(1).id_n2ox,ns),'k'); title(['N2O'])
% subplot(4,1,3);
%     plot(time,y(:,par(1).id_no3x,ns),'y'); title(['NO_3'])
% subplot(4,1,4);
%     plot(time,y(:,par(1).id_DON,ns),'b'); title(['OC'])

figure;
subplot(3,1,1);
    plot(time,y(:,par(1).id_no2x,ns),'m'); title(['NO_2'])
subplot(3,1,2);
    plot(time,y(:,par(1).id_n2ox,ns),'k'); title(['N2O'])
subplot(3,1,3);
    plot(time,y(:,par(1).id_no3x,ns),'y'); title(['NO_3'])


figure;
subplot(3,1,1);
plot(time,y(:,par(1).id_glucose,ns),'r'); title(['Glucose'])
subplot(3,1,2);
plot(time,y(:,par(1).id_acetate,ns),'b'); title(['Acetate'])
subplot(3,1,3);
plot(time,y(:,par(1).id_glutamate,ns),'b'); title(['Glutamate'])

figure
hold on
plot(time(2:end),diff(y(:,par(1).id_nitrificationflux,ns)),'-r');
plot(time(2:end),diff(y(:,par(1).id_denitrificationflux,ns)),'--b');
legend ('N2O-NTR', 'N2O-DNTR')


figure;
for ns=1:nsamp;
    hold on
    plot(time,y(:,par(ns).id_aob_cell(1),ns),'-r');
    plot(time,y(:,par(ns).id_aob_cell(2),ns),'-b');
    plot(time,y(:,par(ns).id_aob_cell(3),ns),'-g');
    plot(time,y(:,par(ns).id_aob_cell(4),ns),'-m');
    title('AOB Biomass');
end


figure
for ns=1:nsamp;
    hold on
    plot(time,y(:,par(ns).id_nob_cell(1),ns),'-r');
    plot(time,y(:,par(ns).id_nob_cell(2),ns),'-b');
    plot(time,y(:,par(ns).id_nob_cell(3),ns),'-g');
    title('NOB Biomass');
end


figure
for ns=1:nsamp;
    hold on
    plot(time,y(:,par(ns).id_anam_cell(1),ns),'-r');
    plot(time,y(:,par(ns).id_anam_cell(2),ns),'-b');
    title('Anammox Biomass');
end

figure
for ns=1:nsamp;
    hold on
    plot(time,y(:,par(ns).id_hetero_cell(1),ns),'-r');
    plot(time,y(:,par(ns).id_hetero_cell(2),ns),'-b');
    plot(time,y(:,par(ns).id_hetero_cell(3),ns),'-g');
    plot(time,y(:,par(ns).id_hetero_cell(4),ns),'--r');
    plot(time,y(:,par(ns).id_hetero_cell(5),ns),'--b');
    plot(time,y(:,par(ns).id_hetero_cell(6),ns),'--g');
    plot(time,y(:,par(ns).id_hetero_cell(7),ns),':r');
    plot(time,y(:,par(ns).id_hetero_cell(8),ns),':b');
    plot(time,y(:,par(ns).id_hetero_cell(9),ns),':g');
    title('Heterotrophic Biomass');
end

% % Export the data as csv files.
% for ns =1:nsamp
% 	 a1(:,:,ns)= y(:,par(ns).id_aob_cell(1),ns);
%      a2(:,:,ns)= y(:,par(ns).id_aob_cell(2),ns);
%      a3(:,:,ns)= y(:,par(ns).id_aob_cell(3),ns);
%      a4(:,:,ns)= y(:,par(ns).id_aob_cell(4),ns);
% end
% % a1 = a(:,:,1);
% % a2 = a(:,:,2);
% % a3 = a(:,:,3);
% %save ('aob.mat', 'a');
% for ns =1:nsamp
%     n1(:,:,ns) = y(:,par(ns).id_nob_cell(1),ns);
%     n2(:,:,ns) = y(:,par(ns).id_nob_cell(2),ns);
%     n3(:,:,ns) = y(:,par(ns).id_nob_cell(3),ns);
% end
% % n1 = n(:,:,1);
% % n2 = n(:,:,2);
% % n3 = n(:,:,3);
% %save ('nob.mat', 'b');
% % for j = nsamp
% %     h(:,:,ns) = y(:,par(ns).id_hetero_cell(1:12),ns);   
% % end
% % h1 = h(:,:,1);
% % h2 = h(:,:,2);
% % h3 = h(:,:,3);
% for ns=1:nsamp;
%     h1(:,:,ns) = y(:,par(ns).id_hetero_cell(1),ns);
%     h2(:,:,ns) = y(:,par(ns).id_hetero_cell(2),ns);
%     h3(:,:,ns) = y(:,par(ns).id_hetero_cell(3),ns);
%     h4(:,:,ns) = y(:,par(ns).id_hetero_cell(4),ns);
%     h5(:,:,ns) = y(:,par(ns).id_hetero_cell(5),ns);
%     h6(:,:,ns) = y(:,par(ns).id_hetero_cell(6),ns);
%     h7(:,:,ns) = y(:,par(ns).id_hetero_cell(7),ns);
%     h8(:,:,ns) = y(:,par(ns).id_hetero_cell(8),ns);
%     h9(:,:,ns) = y(:,par(ns).id_hetero_cell(9),ns);
% end
% for ns = 1: nsamp
%     am1(:,:,ns) = y(:,par(ns).id_anam_cell(1),ns);
%     am2(:,:,ns) = y(:,par(ns).id_anam_cell(2),ns);
% end
% % am1 = am(:,:,1);
% % am2 = am(:,:,2);
% % am3 = am(:,:,3);
% c = y(:,par(1).id_no2x,ns);
% d = y(:,par(1).id_nh3x,ns);
% e = y(:,par(1).id_n2ox,ns);
% f = (time);
% g = y(:,par(1).id_no3x,ns);
% t = (env.tsoi);
% p = (env.pH);
% u = y(:,par(1).id_UREA,ns);
% x = y(:,par(1).id_DON,ns);
% co = y(:,par(1).id_co2x,ns);
% na = y(:,par(1).id_nitrificationflux,ns);
% nd = y(:,par(1).id_denitrificationflux,ns);
% % csvwrite('AOB.csv', a);
% % csvwrite('NOB.csv', n);
% %csvwrite('Hetero.csv', h);
% csvwrite('NO2.csv', c);
% csvwrite('NH3.csv', d);
% csvwrite('N2O.csv', e);
% csvwrite('time.csv', f);
% csvwrite('temp.csv', t);
% csvwrite('pH.csv', p);
% csvwrite('NO3.csv',g);
% csvwrite('CO2.csv', co);
% csvwrite('UREA.csv',u);
% csvwrite('DON.csv', x);
% % csvwrite('Amx.csv', am);
% csvwrite('N2O_NTR.csv', na);
% csvwrite('N2O_DNTR.csv', nd);
% % csvwrite('hetero1.csv', h1);
% % csvwrite('hetero2.csv', h2);
% % csvwrite('hetero3.csv', h3);
% %out2file(t2,y2,filenum);
% csvwrite('Hetero1.csv', h1);
% csvwrite('Hetero2.csv', h2);
% csvwrite('Hetero3.csv', h3);
% csvwrite('Hetero4.csv', h4);
% csvwrite('Hetero5.csv', h5);
% csvwrite('Hetero6.csv', h6);
% csvwrite('Hetero7.csv', h7);
% csvwrite('Hetero8.csv', h8);
% csvwrite('Hetero9.csv', h9);
% csvwrite('aob1.csv', a1);
% csvwrite('aob2.csv', a2);
% csvwrite('aob3.csv', a3);
% csvwrite('nob1.csv', n1);
% csvwrite('nob2.csv', n2);
% csvwrite('nob3.csv', n3);
% csvwrite('am1.csv', am1);
% csvwrite('am2.csv', am2);



