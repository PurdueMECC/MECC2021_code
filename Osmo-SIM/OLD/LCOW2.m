%CapEx calculations RO
dailyvol = 3100*100; %m^3/day
capFac = 0.49; 
CapRO = dailyvol; % capacity factor * m^3/day
capexRO = 3684700*100; % CapEx budget for 100 m^3/day BRO plant
AWP = dailyvol*capFac*365

% OpEx calculations RO
Nl = 16;%(CapRO*264/6e+6)^0.4 + 18/1.4; % Number of laborers
Nm = 31;%(5+ CapRO/5500)/2; % Number of managers
labor = 29700*Nl; % Direct labor costs
manage = 66000*Nm; % Management labor costs
spare = 0.04*AWP; % Spare parts
pre = 0.03*AWP; % Pretreatment
post = 0.01*AWP; % Posttreatment
memb = 0.07*AWP; % Membranes
ins = 0.005*capexRO; % Insurance
opexRO = labor+ manage + pre + post + memb + spare + ins

%WEC
capexWEC = 3877896*100;
opexWEC = 68107*100;

%Totals
capex = capexWEC + capexRO;
opex = opexWEC + opexRO;

FCR = 0.108; % 10.8% fixed charge rate

%LCOW
LCOW = (FCR*capex + opex)/AWP