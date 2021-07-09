%% LCOW
%CAPEX
BatchelorsBudget = 95e3; % dollars - coupling and desal
OSWEC_Cost = 3877896; % dollars - WEC
CAPEX = BatchelorsBudget + OSWEC_Cost;
FCR = 0.108; % 10.8% fixed charge rate

%OPEX
OPEX_Car_D = 1.45*75*365; %OPEX Caribbean desal annual 100 m^3/day
OPEX_Greece_D = 1.19*200*365; %OPEX Greece desal annual 100 m^3/day
Prop_Car = 68107/477843; %Proportion between OPEX desal and OPEX WEC with NREL data
OPEX_Car_WEC = OPEX_Car_D*Prop;  %OPEX Caribbean WEC
OPEX_Greece_WEC = OPEX_Greece_D*Prop; %OPEX Greece WEC
OPEX_Car_tot = OPEX_Car_D + OPEX_Car_WEC; %OPEX Caribbean total 
OPEX_Greece_tot = OPEX_Greece_D + OPEX_Greece_WEC; %OPEX Greece total

LCOW_num_Car = FCR*CAPEX + OPEX_Car_tot; %LCOW numerator Caribbean
LCOW_num_Greece = FCR*CAPEX + OPEX_Greece_tot; %LCOW numerator Greece

V_cycle = V_p_tot.data(end); %Volume produced per cycle
t_cycle = V_p_tot.time(end); %Total cycle length
annual_cycles = (t_cycle/86400/365)^(-1); %cycles per year 
LCOW_den = V_cycle*annual_cycles; %LCOW denominator 

LCOW_Car = LCOW_num_Car/LCOW_den; %LCOW Caribbean
LCOW_Greece = LCOW_num_Greece/LCOW_den; %LCOW Greece






