%% Alcott, Mills and Poulton 2019, Science
% Model based on Slomp and VanCappellen, 2007, Biogeosciences; Tsandev et al., 2008, GBC; Tsandev and Slomp, 2009, EPSL.
% Model equations (do not run this script)

function dy = Alcott_et_al_2019_Science(t,y)
% Set up dy array
dy = zeros(22,1) ;

%%% Set up Global parameters
global stepnumber
global pars
global workingstate
global starting
global per
global present
global sensparams

%% Reservoirs
Water_P = y(1) ;
Water_D = y(2) ;
Water_S = y(3) ;
Water_DP = y(4) ;
POC_P = y(5) ;
POC_D = y(6) ;
POC_S = y(7) ;
POC_DP = y(8) ;
O2_DP = y(12) ;
SRP_P = y(13) ;
OP_P = y(14) ;
SRP_D = y(15) ;
OP_D = y(16) ;
SRP_S = y(17) ;
OP_S = y(18) ;
SRP_DP = y(19) ;
OP_DP = y(20) ;
O2_A = y(21) ;
A = y(22) ;
Aiso = y(23) ;
% Linear relationship for oxygen content in boxes in contact with atmosphere 
O2_P = ( present.O2_P * ( O2_A / present.O2_A ) ) ;
O2_S = ( present.O2_S * ( O2_A / present.O2_A ) ) ;
O2_D = ( present.O2_D * ( O2_A / present.O2_A ) );
Norm_O2_D = O2_D / present.O2_D ;


%% Forcings
tgeol = t/1e6 ; % For now but set to 570 for neoprot

Pforce = per.P ;

% pars.fbiota = interp1([-4.2e9 -450e6 -350e6 0], [0.75 0.75 1 1],t) ; %Need higher degassing rates to accomodate higher pre-land plant weathering
% locb = interp1([-4.2e9 -450e6 -350e6 0], [0 0 1 1],t) ;
carbconst = 0.9 ;
silconst = 0.33 ;

% EXPOSED = interp1([-4.2e9 -3e9 -2.8e9 -1.8e9 -1.79e9 0],[sensparams.EXP sensparams.EXP sensparams.EXP2 sensparams.EXP2 1 1],t) ; %Rough estimate of rapid emergence prior to GOE
EXPOSED = interp1([-4.2e9 -3e9 (-1e9*sensparams.EXPtiming) -1.7e9 -1.6e9 0],[sensparams.EXP sensparams.EXP sensparams.EXP2 1 1 1],t) ; %Rough estimate of rapid emergence prior to GOE

% EXPOSED = interp1([-4.2e9 -3e9 -0.82e9 -0.81e9 -0.8e9 0],[0.1 0.1 0.1 1 1 1],t) ; %Rough estimate of rapid emergence prior to GOE
% EXPOSED = interp1([-4.2e9 -2.5e9 -2.3e9 -1.5e9 0],[0.05 0.2 0.6 1 1],t); %Jacobson, 1988 Use 25%redoxdepend
% EXPOSED = interp1([-4.2e9 -3.1e9 -2.9e9 -1.9e9 -1.7e9 0],[0.05 0.1 0.4 0.5 0.8 1],t); %Ying et al 2011 Use 25%redoxdepend
% EXPOSED = interp1([-4.2e9 -2.7e9 -2.4e9 0],[0.1 0.35 0.75 1],t);% Taylorand McLennan 1985 Use 25%redoxdepend
% EXPOSED = interp1([-4.2e9 -3.5e9 -3e9 -2e9 0],[0.3 0.5 0.7 0.8 1],t); %Campbell 2003 Use 50%redoxdepend
% % EXPOSED = interp1([-4.2e9 -2e9 -1e9 0], [0.01 0.5 0.5 1],t) ;

% D = interp1([-4.2e9 -3e9 -2e9 0],[(12* sensparams.D) (5* sensparams.D) (2* sensparams.D) 1],t) * sensparams.D;
D = interp1([-4.2e9 -3e9 -2e9 0],[12 5 2 1],t) * sensparams.D;

% SET UP FOR NEOPROTEROZOIC LIMIT CYCLES 0.5 0.5
% pars.CPanoxic_prox = 1100 * sensparams.CP ;
% pars.CPanoxic_dist = 1100 * sensparams.CP ;
% pars.CPanoxic_deep = 1100 * sensparams.CP ;
pars.CPanoxic_prox = sensparams.CP ;
pars.CPanoxic_dist = sensparams.CP ;
pars.CPanoxic_deep = sensparams.CP ;

% D = interp1([-4.2e9 -3e9 -2e9 0],[24 10 4 1],t) * sensparams.D;
U = 1 ;
PG = 1 ;
% EXPOSED = 1 ;
% C = interp1([-4.2e9 -4e9   0],[0.001 0.001 1],t) ^ sensparams.C;

C = interp1([per.C_HW_2006_time],[per.C_HW_2006_data],t)^ sensparams.C ;
% C = interp1([per.C_HW_2006_time],[per.C_HW_2006_data],t)^ 4 ;

% pars.fbiota = interp1([-4.2e9 -450e6 -350e6 0], [0.15 0.15 1 1],t) * sensparams.fbiota;
pars.fbiota = interp1([-4.2e9 -450e6 -350e6 0],[sensparams.fbiota sensparams.fbiota 1 1],t) ;
locb = interp1([-4.2e9 -450e6 -350e6 0], [0 0 1 1],t) ;
% locb = interp1([-4.2e9 -450e6 -350e6 0], [sensparams.locb sensparams.locb 1 1],t) ;
% 
% % SET UP FOR GOE LIMIT CYCLES 0.5Porg 0.5Pauth
% pars.CPanoxic_prox = 4000 ;
% pars.CPanoxic_dist = 4000 ;
% pars.CPanoxic_deep = 4000  ;
% D = interp1([-4.2e9 -3e9 -2e9 0],[24 10 4 1],t) ;
% U = 1 ;
% PG = 1 ;
% EXPOSED = 1 ;
% C = interp1([-4.2e9  0],[0.01 1],t);
% pars.fbiota = interp1([-4.2e9 -450e6 -350e6 0], [0.15 0.15 1 1],t) ;
% locb = interp1([-4.2e9 -450e6 -350e6 0], [0 0 1 1],t) ;

% % %0.5SCAV on with limiter 0.9Porg 0.9Pauth
% pars.CPanoxic_prox = 1100 ;
% pars.CPanoxic_dist = 1100 ;
% pars.CPanoxic_deep = 1100  ;
% D = interp1([-4.2e9 -3e9 -2e9 0],[16 8 5 1],t) ;
% U = 1 ;
% PG = 1 ;
% EXPOSED = 1 ;
% C = interp1([-4.2e9  0],[0.01 1],t);
% pars.fbiota = interp1([-4.2e9 -450e6 -350e6 0], [0.75 0.75 1 1],t) ;
% locb = interp1([-4.2e9 -450e6 -350e6 0], [0 0 1 1],t) ;

% % %EDIT
% pars.CPanoxic_prox = 250 ;
% pars.CPanoxic_dist = 250 ;
% pars.CPanoxic_deep = 250  ;
% D = interp1([-4.2e9 -3e9 -2e9 0],[12 5 2 1],t) ;
% U = 1 ;
% PG = 1 ;
% EXPOSED = 1 ;
% C = interp1([-4.2e9  0],[0.01 1],t);
% pars.fbiota = interp1([-4.2e9 -450e6 -350e6 0], [0.75 0.75 1 1],t) ;
% locb = interp1([-4.2e9 -450e6 -350e6 0], [0 0 1 1],t) ;

% %Present day
% pars.CPanoxic_prox = 106 ;
% pars.CPanoxic_dist = 106 ;
% pars.CPanoxic_deep = 106  ;
% pars.CPanoxic_prox = 1100 ;
% pars.CPanoxic_dist = 1100 ;
% pars.CPanoxic_deep = 1100  ;
% D = 1 ;
% U = 1 ;
% PG = 1 ;
% EXPOSED = 1 ;
% C = 1;
% pars.fbiota = 1;
% locb = 1;
% tgeol = 0 ;

%% Concentrations 
present.Conc_O2_deep = present.O2_DP / starting.Water_DP ;
O2_Pconc = O2_P / Water_P ;
O2_Dconc = O2_D / Water_D ;
O2_Sconc = O2_S / Water_S ;
O2_DPconc = O2_DP/Water_DP;                        
OP_Dconc  = y(16)/y(2);                           
OP_Pconc  = y(14)/y(1);                           
SRP_DPconc = y(19)/y(4);                           
SRP_Dconc = y(15)/y(2);                            
SRP_Pconc = y(13)/y(1);                            
SRP_Sconc = y(17)/y(3);                            
    

%% Oceanic Water Cycle
%%% As in Slomp and Van Cappellen, 2007

%River flow entering Water_P
River_Water = 37e12 ; %Berner and Berner1996

%flow from Water_P to Water_D
Water_P_D = River_Water ; 

%Coastal Upwelling
Water_DP_D = pars.kWF6 * Water_DP ; %12Sv -  Brink et al1995

%flow from Water_D to Water_S
Water_D_S = Water_P_D + Water_DP_D ;

%Oceanic Upwelling
Water_DP_S = Water_DP* pars.kWF5 ; %120Sv - Brink et al1995

%downwelling (Water_S to Water_DP)
Water_S_DP = pars.kWF4 * Water_S ;

%Evaporation from Water_S
Evaporation_Water = River_Water ;

%% Hydrological Cycle

%Proximal Zone
dy(1) = River_Water - Water_P_D ;

%Distal Zone
dy(2) = Water_P_D + Water_DP_D - Water_D_S ;

%Surface Ocean
dy(3) = Water_DP_S + Water_D_S - Water_S_DP - Evaporation_Water ;

%Deep Ocean
dy(4) = Water_S_DP - Water_DP_S - Water_DP_D ;


%% Normalised to present day reservoirs
Norm_SRP_P = SRP_P / present.SRP_P ;
Norm_SRP_D = SRP_D / present.SRP_D ;
Norm_SRP_S = SRP_S / present.SRP_S ;
Norm_OP_P = OP_P / present.OP_P ;
Norm_OP_D = OP_D / present.OP_D ;
Norm_OP_S = OP_S / present.OP_S ;
Norm_POC_P = POC_P / present.POC_P ;
Norm_POC_D = POC_D / present.POC_D ;
Norm_POC_S = POC_S / present.POC_S ;
Norm_POC_DP = POC_DP / present.POC_DP ;
Norm_O2_A = O2_A / present.O2_A ;


%% Marine Carbon Cycle

% Primary production in Proximal
PP_P = pars.kPhotoprox * Norm_SRP_P * pars.Redfield_CP ; 

% POC mineralisation in Proximal
POC_Min_P = pars.kminprox * Norm_POC_P ;

% POC export from Proximal to Distal
OP_P_D = Water_P_D * OP_Pconc ;
XP_P_D = OP_P_D * pars.Redfield_CP ; 

% Proximal sediment POC burial
POC_P_Burial = pars.Prox_C_Bur * PP_P ; 

% Primary Production in Distal
PP_D = pars.kPhotodist * Norm_SRP_D * pars.Redfield_CP ; 

% POC mineralisation in Distal
POC_Min_D = pars.kmindist * Norm_POC_D ; 

% POC Export from Distal to Surface
OP_D_S = Water_D_S * OP_Dconc ;
XP_D_S = OP_D_S * pars.Redfield_CP ; 

% Distal sediment POC burial
POC_D_Burial = pars.Dist_C_Bur * ( XP_P_D + PP_D ) ; 

% Primary Production in Surface
PP_S = pars.kPhotosurf * Norm_SRP_S * pars.Redfield_CP ; 

% POC mineralisation in Surface
POC_Min_S = pars.kminsurf * Norm_POC_S ; 

% POC export from Surface to Deep
XP_S_DP = pars.Surf_Deep_XP * ( XP_D_S + PP_S ) ; 

% POC respiration in Water_DP
POC_DP_Resp = pars.kCF12 * Norm_POC_DP ; 

%% Deep Carbon Burial

% C:P ratios for oxic and anoxic conditions

% invlogO2 = -(log10(Norm_O2_A));
% 
% if t <= -3e9 
%     pars.CPanoxic_prox = 1000  ;
%     pars.CPanoxic_dist = 1000  ;
%     pars.CPanoxic_deep = 1000 ;
% else
%     CPanoxic = (4 / (1 + exp(-2*(invlogO2 - 4))) ) *1000;  
%     pars.CPanoxic_prox = 1000 + CPanoxic ;
%     pars.CPanoxic_dist = 1000 + CPanoxic ;
%     pars.CPanoxic_deep = 1000 ; %assume deep ocean is always ferruginous
% end





if O2_DPconc < present.Conc_O2_deep
    
    % Deep sediment POP burial
    % Redox dependancy is set up according to Slomp and VC, 2007 and Tsandev et al., 2009. 
    OP_DP_Burial = ( pars.kPOP_Bur_Deep * XP_S_DP / pars.CPoxic ) * ( (1-per.POP_deep_feedback) + (per.POP_deep_feedback * O2_DPconc / present.Conc_O2_deep ));
    
    % Deep sediment FeP burial
    P_FeP_DP = pars.kFeP_Deep * ( O2_DPconc / present.Conc_O2_deep );
    
    % Deep sediment POC burial
    POC_DP_Burial = OP_DP_Burial * ( ( pars.CPanoxic_deep * pars.CPoxic ) / ( ( ( O2_DPconc / present.Conc_O2_deep ) * pars.CPanoxic_deep ) + ( ( 1 - ( O2_DPconc / present.Conc_O2_deep ) ) * pars.CPoxic ) ) ) ;
        
else
     OP_DP_Burial= pars.kPOP_Bur_Deep * XP_S_DP / pars.CPoxic;  
     P_FeP_DP = pars.kFeP_Deep ;
     POC_DP_Burial = pars.CPoxic * OP_DP_Burial ;    
end

%Carbon Proximal Zone
dy(5) = PP_P - POC_Min_P - POC_P_Burial  - XP_P_D ;

%Carbon Distal Zone
dy(6) = XP_P_D + PP_D - POC_D_Burial - POC_Min_D - XP_D_S ; 

%Carbon Surface Ocean
dy(7) = XP_D_S + PP_S - POC_Min_S - XP_S_DP ;

%Carbon Deep Ocean
dy(8) = XP_S_DP - POC_DP_Resp - POC_DP_Burial ;

POCTotal = PP_P - POC_Min_P - POC_P_Burial  - XP_P_D + XP_P_D + PP_D - POC_D_Burial - POC_Min_D - XP_D_S +XP_D_S + PP_S - POC_Min_S - XP_S_DP + XP_S_DP - POC_DP_Resp - POC_DP_Burial;


%% Oxygen Cycle



% fanoxic parameters (From Watson et al., 2017)
kanox = 10 ; 
O2O20 = Norm_O2_A ;
kU = 0.4 ;

% fanoxic calculation from Watson et al., 2017
fanoxicdist = 1 / ( 1 + exp(-kanox * ( kU * Norm_SRP_D - O2O20 ) ) ) ; 
fanoxicprox = 1 / ( 1 + exp(-kanox * ( kU * Norm_SRP_P - O2O20 ) ) ) ; 

present.fanoxicprox = 0.0025 ;
present.fanoxicdist = 0.0025 ;

% if Norm_O2_A <=1e-6
%     pars.CPanoxic_prox = 1000  ;
%     pars.CPanoxic_dist = 1000  ;
%    
% else
%     pars.CPanoxic_prox = (4 / (1 + exp(-100*(fanoxicprox - 0.1))) ) *1000;
%     pars.CPanoxic_dist = (4 / (1 + exp(-100*(fanoxicdist - 0.1))) ) *1000;
% end


%Concentration of oxygen 
O2_Sconc = O2_S/Water_S ;

%O2 Downwelling
O2_S_DP = Water_S_DP * O2_Sconc ;  

%O2 coastal upwelling
O2_DP_D =  Water_DP_D * O2_DPconc; 

%O2 oceanic upwelling
O2_DP_S = Water_DP_S * O2_DPconc; 

%Aerobic O2 respiration
Respiration_O21 = pars.kCF12 * Norm_POC_DP / pars.Redfield_CO2; % Preliminary Respiration quanitity before Monod inclusion
KmO2 = 0.0001 ; % monod constant for oxic respiration in mol/m3 
Mon_O2_deep = O2_DPconc / ( KmO2 + O2_DPconc ) ;
Respiration_O2 = Respiration_O21 * Mon_O2_deep ;

FeO = pars.FeO * D * ( O2_DPconc / present.Conc_O2_deep ) ;

%% Oxygen dys

%Oxygen Deep Ocean 
dy(12) = - Respiration_O2 - O2_DP_S - O2_DP_D + O2_S_DP - FeO ;

Atmos_Weather = pars.O2_A_Weathering * (sqrt(O2_A/ present.O2_A)) * EXPOSED ;

% generic_reductant_flux = generic_reductant_flux_i * sigmf((log10(Norm_O2_A)),[3,-5]) ; 

Total_POC_Burial = POC_P_Burial + POC_D_Burial + POC_DP_Burial ; 


logNormAtmosO2 = log10(Norm_O2_A) ;
if logNormAtmosO2 <=0
    logNormAtmosO2 = -10 ;
else
end
    save = sigmf(logNormAtmosO2,[3,-5]) ;


FrgfO2 = pars.rgf * D * save ;
Frgf = pars.rgf * D ;

% Land Corg burial 
Flocb = pars.Flocb_0 * locb ;

%% Oxygen atmosphere
dy(21) = Total_POC_Burial - Atmos_Weather - FrgfO2 + Flocb ;

%% Carbon isotope

%Oxidative weathering
Foxidw = Atmos_Weather ; %Used to match oxidative weathering for the oxygen model
% Fmocb = Total_POC_Burial ;
%CO2 ppm calculation
CO2 = (A / starting.A_0)^2 ;
CO2ppm = CO2 * 280 ;

%%%%%%% atmospheric fraction of total CO2, atfrac(A)
atfrac0 = 0.01614 ;
%%%%%%% variable
atfrac = atfrac0 * (A/starting.A_0) ;

%%%%%%%% calculations for pCO2, pO2
RCO2 = (A/starting.A_0)*(atfrac/atfrac0) ;
CO2atm = RCO2*(280e-6) ;

%Climate temp adjustment from COPSE
kclim = 5 ; %Kelvin

%COPSE method of temperature.

constant = 7.4 ;
% constant = 3 ;
% GAST = 288 + (kclim * (  ( log10 (CO2ppm / 280) ) / log10(2) )) - (( constant * ( tgeol / -570 ) )) ;

GAST = 288 + (kclim * (  ( log(RCO2) ) / log(2) )) - (( constant * ( tgeol / -570 ) )) ;
%%%% relcalculate low lat temp for surface processes 
tgrad = 0.66 ;
tgrad = 1 ;
tc = 288*(1-tgrad) ;
Tsurf = GAST*tgrad + tc + 10 ; %%% low lat temp 25C at present
DT = GAST - 288 ;
DTsurf = Tsurf - 298 ;

%Basalt weathering temp dependance
fTbas = exp(pars.kTbas * (DT))  ;

%Granite weathering temp dependance
fTgran = exp(pars.kTgran * (DT))  ;

%Additional weathering flux dependances
% k2 = 0.038 ;
% frunoff = (1+  exp( k2*DT ) )^0.65 ; %From Bens Snowball paper
% % frunoff = ( exp ( 0.038*( DTsurf ) ) )^0.4 ;
% % grunoff = (1+ exp( 0.1*DT ) ) ^0.5 ; %Best guess quick fit to previous linear grunoff
% g_T = (1+  exp( k2*DT ) )^0.65  ;


f_T_bas =  exp(0.0608*(DTsurf)) * ( (1 + 0.038*(DTsurf))^0.65 ) ; %%% 42KJ/mol
f_T_gran =  exp(0.0724*(DTsurf)) * ( (1 + 0.038*(DTsurf))^0.65 ) ; %%% 50 KJ/mol
% g_T = (1+  exp( 0.038*DTsurf ) )^0.65 ;
% g_T = 1 + 0.087*(DT) ;
g_T = exp( 0.05*DTsurf ) ; %close fit to COPSE but doesnt go below zero


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%QUICK FIX FOR NOW
if f_T_bas <=0
    f_T_bas = 0 ;
else
    f_T_bas = f_T_bas ;
end
if f_T_gran <=0
    f_T_gran = 0 ;
else
    f_T_gran = f_T_gran ;
end

%Carbonate weathering
Fcarbw = pars.Fcarbw_0 * ( C ) * (U^carbconst) * PG * pars.fbiota * g_T * EXPOSED ; 

% %Granite weathering
% Fgranw = pars.k_granw * (U^silconst) * PG * fTgran * pars.fbiota * frunoff * EXPOSED; % 70% is granite weathering of the total 4e12 FREE PARAM
% 
% %Basalt weathering
% Fbasw = pars.k_basw * PG * fTbas * pars.fbiota * frunoff * EXPOSED; %30% is basalt weathering FREE PARAM


Fbasw = pars.k_basw * PG * pars.fbiota * f_T_bas * EXPOSED ;
Fgranw = pars.k_granw * (U^silconst) * PG * pars.fbiota * f_T_gran * EXPOSED ;
%Silicate weathering
Fsilw = Fgranw + Fbasw ;

%Total weathering
Fw = Fcarbw + Fsilw ; 

%Carbonate carbon degassing
Fccdeg = pars.Fccdeg_0 * D * C ; %6e12 used by Shields et al for NeoProt.

Focdeg = pars.Focdeg_0 * D  ;

%Total organic carbon burial
Fmocb = (POC_P_Burial + POC_D_Burial + POC_DP_Burial) ; %Used to match organic carbon burial in the current model

%Carbonate burial
Fmccb = Fw ; 

%Sea floor weathering
Fsfw = pars.Fsfw_0 * D * exp(pars.kTsfw * (DT));

%Carbon isotope fractionations for 1 box COPSE
dG = -25 ;
dC = 1 ;
dA = Aiso ;
DB = 25 ;
dM = -5 ; %Mantle value

%%%% carbonate fractionation
delta_o = atfrac*( (9483/GAST)  - 23.89 ) ; 
d_mccb = delta_o + 15.1 - (4232/GAST) ; %%%% calcite burial
% delta_mccb = dA + d_mccb ;
delta_mccb = dA ;

%%%% marine organic
capdelB_0 = 33 ;
Jparam =  5 ;
e_B_co2 = -9 / sqrt(RCO2) ;
e_o2 = Jparam * ( (O2_A/starting.O2_A) - 1) ;
DB = capdelB_0 + e_B_co2 + e_o2 ;
%%%% final calc
d_mocb = d_mccb - DB ;
delta_mocb = dA + d_mocb ;
% d_mocb = dA ;

%%%% land plant
capdelP_0 = 19 ;
capdelP = capdelP_0 + e_o2;
%%%% atmospheric
delta_a = (atfrac-1)*( (9483/GAST)  - 23.89 ) ; 
%%%% final calc
d_locb = delta_a - capdelP ;
delta_locb = dA + d_locb ;


dy(22) = Foxidw  + Focdeg + Fcarbw + Fccdeg - Fmocb - Fmccb - Fsfw + Frgf - Flocb ;
% dy(23) = ((Foxidw*dG) + (Focdeg*dG) + (Fcarbw*dC) + (Fccdeg*dC) + (Frgf*dG) - (Fmocb*(dA-DB)) - (Fmccb*dA) - (Fsfw*dA) - (Flocb*(dA-DB)))/A ;
dy(23) = ((Foxidw*dG) + (Focdeg*dG) + (Fcarbw*dC) + (Fccdeg*dC) + (Frgf*dM) - (Fmocb*delta_mocb) - (Fmccb*delta_mccb) - (Fsfw*delta_mccb) - (Flocb*(delta_locb)) )/A ;
% dy(23) = ((Foxidw*dG) + (Focdeg*dG) + (Fcarbw*dC) + (Fccdeg*dC) + (Frgf*dM) - (0*delta_mocb) - (Fmccb*delta_mccb) - (Fsfw*delta_mccb) - (0*(delta_locb)) )/A ;


%% P Cycle

% SRP upwelling Water_DP to Water_S
SRP_DP_S = SRP_DPconc * Water_DP_S ; 

% SRP upwelling Water_DP to Water_D
SRP_DP_D = Water_DP_D * SRP_DPconc ; 


%% Proximal Coastal
%O2 limit for scavenging
scavlim = 1e-3 ;
Pfrac_silw = 2/12 ;
Pfrac_carbw = 5/12 ;
Pfrac_oxidw = 5/12 ;

% Forcing to vary Riverine P input.
% River_SRP = pars.River_SRP_0 * Pforce ;
River_SRP = pars.River_SRP_0 * Pforce * ( (Pfrac_silw * (Fsilw/pars.Fsilw_0) ) + ( Pfrac_carbw * (Fcarbw/pars.Fcarbw_0) ) + ( Pfrac_oxidw * (Foxidw/pars.O2_A_Weathering) ) );

% Primary Production proximal
P_PP_P =  PP_P/pars.Redfield_CP ; 

% POP mineralisation
OP_P_Min = Norm_OP_P * pars.kPrel_prox  ; 

% SRP transport from Proximal to Distal
SRP_P_D = SRP_Pconc * Water_P_D ;

% Proximal sediment POP burial 
OP_P_Burial = pars.Prox_C_Bur * PP_P * ( ( ( 1-fanoxicprox ) / pars.CPoxic ) + ( fanoxicprox / pars.CPanoxic_prox ) ) ;
% OP_P_Burial = pars.kOP_P_0 * PP_P * ( ( ( 1-fanoxicprox ) / pars.CPoxic ) + ( fanoxicprox / pars.CPanoxic_prox ) ) ;

% Proximal FeP burial
P_FeP_P = pars.kFePprox * SRP_P * ( 1- fanoxicprox );

% Proximal CaP burial
% P_AuthP_P = pars.kPrel_prox * OP_P * pars.kCaP_prox ;  %OLD
P_AuthP_P = pars.kCaP_P * OP_P_Min * ( 1 - fanoxicprox ) ; %NEW

%% Distal Coastal

% Primary Production in Water_D
P_PP_D = PP_D / pars.Redfield_CP ; 

% POP mineralisation in Water_D
OP_D_Min = Norm_OP_D * pars.kPrel_dist ; 

% SRP transport Water_D to Water_S
SRP_D_S = SRP_Dconc * Water_D_S ;

% Distal sediment POP burial
OP_D_Burial = pars.kPOPDOADist * (PP_D + XP_P_D) * ( ( ( 1-fanoxicdist ) / pars.CPoxic ) + ( fanoxicdist / pars.CPanoxic_dist ) ) ;

% Distal sediment FeP burial
P_FeP_D = pars.kFePDOADist * SRP_D * ( 1- fanoxicdist ) ; 

% Distal Sediment CaP burial
P_AuthP_D = pars.kCaPDOADist * OP_D_Min * ( 1 - fanoxicdist ) ;

%% Surface Ocean

% Primary Production 
P_PP_S = PP_S / pars.Redfield_CP; % Links P to C cycle

% POP mineralisation 
OP_S_Min = Norm_OP_S * pars.kPrel_surf ; 

% SRP downwelling from Water_S to Water_DP
SRP_S_DP = SRP_Sconc * Water_S_DP ;

% POP export from Water_S to Water_DP
OP_S_DP = pars.kCF11 * ( PP_S + XP_D_S ) / pars.Redfield_CP ; % Links P to C cycle

% Scavenging flux of FeP from surface ocean.
eSCAV_SURF = 1 / ( 1 + exp(5000*(O2_DPconc - scavlim ))) ; % Simple logistic curve, similar to Reinhard et al., 2017
Fe_SCAV_SURF = min( eSCAV_SURF * SRP_DP_S * per.sig_SCAV , 4.5e10 ); % upper limit variable
% Fe_SCAV_SURF = min( eSCAV_SURF * SRP_DP_S * per.sig_SCAV ); % upper limit variable

%% Deep Ocean

% POP mineralisation in Water_DP
OP_DP_Min = pars.kPrel_deep * OP_DP ;

% Deep sediment CaP burial
% Redox dependancy is set up according to Slomp and VC, 2007 and Tsandev et al., 2009. 
P_AuthP_DP =  pars.fPF34 * OP_DP_Min * ( (1-per.CaP_deep_feedback) + ( per.CaP_deep_feedback * ( O2_DPconc/present.Conc_O2_deep ) ) ) ; 



%% Phosphorus Differentials

    %SRP Proximal
    dy(13) = River_SRP - P_PP_P + OP_P_Min - P_FeP_P - P_AuthP_P - SRP_P_D ; 

    %POP Proximal
    dy(14) = P_PP_P - OP_P_Min - OP_P_Burial - OP_P_D  ;

    %SRP Distal
    dy(15) = SRP_P_D - P_PP_D + OP_D_Min - P_FeP_D - P_AuthP_D - SRP_D_S + SRP_DP_D ;

    %POP Distal
    dy(16) = OP_P_D + P_PP_D - OP_D_Min - OP_D_Burial - OP_D_S ;

    %SRP Surface Ocean
    dy(17) = SRP_D_S - P_PP_S + OP_S_Min - SRP_S_DP + SRP_DP_S ;

    %POP Surface Ocean
    dy(18) = OP_D_S + P_PP_S - OP_S_Min - OP_S_DP;

    %SRP Deep Ocean
    dy(19) = SRP_S_DP + OP_DP_Min - P_FeP_DP - P_AuthP_DP - SRP_DP_S - SRP_DP_D - Fe_SCAV_SURF;

    %POP Deep Ocean
    dy(20) = OP_S_DP - OP_DP_Min - OP_DP_Burial ;


    
%     
% dy(1) = 0 ;
% dy(2) = 0 ;
% dy(3) = 0 ;
% dy(4) = 0 ;
% dy(5) = 0 ;
% dy(6) = 0 ;
% dy(7) = 0 ;
% dy(8) = 0 ;
% dy(9) = 0 ;
% dy(10) = 0 ;
% dy(11) = 0 ;
% dy(12) = 0 ;
% dy(13) = 0 ;
% dy(14) = 0 ;
% dy(15) = 0 ;
% dy(16) = 0 ;
% dy(17) = 0 ;
% dy(18) = 0 ;
% dy(19) = 0 ;
% dy(20) = 0 ;
% dy(21) = 0 ;
% dy(22) = 0 ;
% dy(23) = 0 ;
%% Saving data
workingstate.FrgfO2(stepnumber,1) = FrgfO2 ;
% 
% workingstate.POC_P(stepnumber,1) = POC_P ;
% workingstate.POC_D(stepnumber,1) = POC_D ;
% workingstate.POC_S(stepnumber,1) = POC_S ;
% workingstate.POC_DP(stepnumber,1) = POC_DP ;
% workingstate.PP_P(stepnumber,1) = PP_P ;
% 
% workingstate.PP_D(stepnumber,1) = PP_D ;
% 
% workingstate.PP_S(stepnumber,1) = PP_S ;

% workingstate.O2_P(stepnumber,1) = O2_P ;
% workingstate.O2_D(stepnumber,1) = O2_D ;
% workingstate.O2_S(stepnumber,1) = O2_S ;
workingstate.O2_DP(stepnumber,1) = O2_DP ;
workingstate.O2_A(stepnumber,1) = O2_A ;
% workingstate.Respiration_O2(stepnumber,1) = Respiration_O2 ;
% workingstate.SRP_P(stepnumber,1) = SRP_P ;
% workingstate.OP_P(stepnumber,1) = OP_P ;
% workingstate.SRP_D(stepnumber,1) = SRP_D ;
% workingstate.OP_D(stepnumber,1) = OP_D ;
% workingstate.SRP_S(stepnumber,1) = SRP_S ;
% workingstate.OP_S(stepnumber,1) = OP_S ;
workingstate.SRP_DP(stepnumber,1) = SRP_DP ;
% workingstate.OP_DP(stepnumber,1) = OP_DP ;
workingstate.River_SRP(stepnumber,1) = River_SRP ;

% workingstate.CP_Dist(stepnumber,1) = ( ( ( 1-fanoxicdist ) * pars.CPoxic ) + ( fanoxicdist * pars.CPanoxic_dist ) ) ;
% workingstate.CP_Deep(stepnumber,1) = POC_DP_Burial / OP_DP_Burial ;
% workingstate.CP_Prox(stepnumber,1) = ( ( ( 1-fanoxicprox ) * pars.CPoxic ) + ( fanoxicprox * pars.CPanoxic_prox ) ) ;
workingstate.Deep_Preac_Burial(stepnumber,1) = P_AuthP_DP + P_FeP_DP + OP_DP_Burial ;
workingstate.Dist_Preac_Burial(stepnumber,1) = P_AuthP_D + P_FeP_D + OP_D_Burial ;
workingstate.Prox_Preac_Burial(stepnumber,1) = P_AuthP_P + P_FeP_P + OP_P_Burial ;
workingstate.fanoxicdist(stepnumber,1) = fanoxicdist ;
workingstate.fanoxicprox(stepnumber,1) = fanoxicprox ;
% workingstate.Atmos_Weather(stepnumber,1) = Atmos_Weather ;
% workingstate.Fe_SCAV_SURF(stepnumber,1)= Fe_SCAV_SURF ;

workingstate.CO2atm(stepnumber,1) = CO2atm ;
workingstate.GAST(stepnumber,1) = GAST ;
workingstate.A(stepnumber,1) = A ;
workingstate.Aiso(stepnumber,1) = Aiso ;
workingstate.Fmocb(stepnumber,1) = Fmocb ;
% workingstate.Flocb(stepnumber,1) = Flocb ;
% workingstate.Fcarbw(stepnumber,1) = Fcarbw ; 
% workingstate.Fgranw(stepnumber,1) = Fgranw ;
% workingstate.Fbasw(stepnumber,1) = Fbasw ;
% workingstate.Fsilw(stepnumber,1) = Fsilw ;
% workingstate.Foxidw(stepnumber,1) = Foxidw ;
% workingstate.Fsfw(stepnumber,1) = Fsfw ;
% workingstate.Fmccb(stepnumber,1) = Fmccb ;
% workingstate.Fccdeg(stepnumber,1) = Fccdeg ;
% workingstate.Frgf(stepnumber,1) = Frgf ;

% workingstate.fbiota(stepnumber,1) = pars.fbiota ;
% workingstate.D(stepnumber,1) = D ;
% workingstate.C(stepnumber,1) = C ;
% workingstate.U(stepnumber,1) = U ;
% workingstate.EXPOSED(stepnumber,1) = EXPOSED ;
% workingstate.PG(stepnumber,1) = PG ;
workingstate.CPanoxic(stepnumber,1) = pars.CPanoxic_deep ;

% workingstate.doxidw(stepnumber,1) = (Foxidw*dG) ;
% workingstate.docdeg(stepnumber,1) = (Focdeg*dG) ;
% workingstate.dcarbw(stepnumber,1) = (Fcarbw*dC) ;
% workingstate.dccdeg(stepnumber,1) = (Fccdeg*dC) ;
% workingstate.drgf(stepnumber,1) = (Frgf*dG) ;
% workingstate.dmocb(stepnumber,1) = (Fmocb*(delta_mocb)) ;
% workingstate.dmccb(stepnumber,1) = (Fmccb*delta_mccb) ;
% workingstate.dsfw(stepnumber,1) = (Fsfw*dA);
% workingstate.dlocb(stepnumber,1) = (Flocb*(delta_locb));
% workingstate.delta_mccb(stepnumber,1) = delta_mccb ;
% workingstate.delta_mocb(stepnumber,1) = delta_mocb ;
% workingstate.delta_locb(stepnumber,1) = delta_locb ;
workingstate.forg(stepnumber,1) = Fmocb / (Fmccb+Fmocb+Fsfw);


% workingstate.P_PP_P(stepnumber,1) = P_PP_P ;
% workingstate.OP_P_Min(stepnumber,1) = OP_P_Min ;
% workingstate.OP_P_Burial(stepnumber,1) = OP_P_Burial ;
% workingstate.P_FeP_P(stepnumber,1) = P_FeP_P ;
% workingstate.P_AuthP_P(stepnumber,1) = P_AuthP_P ;
% workingstate.SRP_P_D(stepnumber,1) = SRP_P_D ;
% workingstate.OP_P_D(stepnumber,1) = OP_P_D ;
% workingstate.P_PP_D(stepnumber,1) = P_PP_D ;
% workingstate.OP_D_Min(stepnumber,1) = OP_D_Min ;
% workingstate.OP_D_Burial(stepnumber,1) = OP_D_Burial ;
% workingstate.P_FeP_D(stepnumber,1) = P_FeP_D ;
% workingstate.P_AuthP_D(stepnumber,1) = P_AuthP_D ;
% workingstate.SRP_D_S(stepnumber,1) = SRP_D_S ;
% workingstate.OP_D_S(stepnumber,1) = OP_D_S ;
% workingstate.P_PP_S(stepnumber,1) = P_PP_S ;
% workingstate.OP_S_Min(stepnumber,1) = OP_S_Min ;
% workingstate.SRP_S_DP(stepnumber,1) = SRP_S_DP ;
% workingstate.OP_S_DP(stepnumber,1) = OP_S_DP ;
% workingstate.OP_DP_Min(stepnumber,1) = OP_DP_Min ;
% workingstate.OP_DP_Burial(stepnumber,1) = OP_DP_Burial ;
% workingstate.P_FeP_DP(stepnumber,1) = P_FeP_DP ;
% workingstate.P_AuthP_DP(stepnumber,1) = P_AuthP_DP ;
% workingstate.SRP_DP_S(stepnumber,1) = SRP_DP_S ;
% workingstate.SRP_DP_D(stepnumber,1) = SRP_DP_D ;

% workingstate.O2_S_DP(stepnumber,1) = O2_S_DP ;
% workingstate.O2_DP_D(stepnumber,1) = O2_DP_D ;
% workingstate.O2_DP_S(stepnumber,1) = O2_DP_S ;

% workingstate.POC_Min_S(stepnumber,1) = POC_Min_S ;
% workingstate.XP_S_DP(stepnumber,1) = XP_S_DP ;
% workingstate.POC_DP_Resp(stepnumber,1) = POC_DP_Resp ;
% workingstate.POC_DP_Burial(stepnumber,1) = POC_DP_Burial ;
% workingstate.POC_Min_D(stepnumber,1) = POC_Min_D ;
% workingstate.POC_D_Burial(stepnumber,1) = POC_D_Burial ;
% workingstate.XP_D_S(stepnumber,1) = XP_D_S ;
% workingstate.POC_Min_P(stepnumber,1) = POC_Min_P ;
% workingstate.POC_P_Burial(stepnumber,1) = POC_P_Burial ;
% workingstate.XP_P_D(stepnumber,1) = XP_P_D ;

% workingstate.Water_P(stepnumber,1) = Water_P ;
% workingstate.Water_D(stepnumber,1) = Water_D ;
% workingstate.Water_S(stepnumber,1) = Water_S ;
% workingstate.Water_DP(stepnumber,1) = Water_DP;
% workingstate.River_Water(stepnumber,1) = River_Water ;
% workingstate.Prox_Water_D(stepnumber,1) = Water_P_D ;
% workingstate.Dist_Water_S(stepnumber,1) = Water_D_S ;
% workingstate.Surf_Water_DP(stepnumber,1) = Water_S_DP ;
% workingstate.Water_DP_S(stepnumber,1) = Water_DP_S ;
% workingstate.Water_DP_D(stepnumber,1) = Water_DP_D ;
% workingstate.Evaporation_Water(stepnumber,1) = Evaporation_Water ;
workingstate.Total_POC_Burial(stepnumber,1) = Total_POC_Burial ;

%%%%%%% record time
workingstate.time(stepnumber,1) = t ;

%%% final action: record current model step
stepnumber = stepnumber + 1 ;



end
