%% PTO-Sim Input File 
% Edited by: Emily Bywater, Alec Lanter, Hayden Schennum
% 4/8/21

%% Non-Compressible Fluid Hydraulic PTO-Sim
ptosim = ptoSimClass('Non-Compressible Fluid Hydraulic');

%% Rotary to Linear Crank
ptosim.motionMechanism.crank = 3;
ptosim.motionMechanism.offset = 1.3;
ptosim.motionMechanism.rodLength = 5;

%% Piston 
ptosim.pistonNCF.topA = 0.0378;                                                      % Top piston area [m^2]
ptosim.pistonNCF.botA = 0.0378;                                                      % Bottom piston area [m^2]

%% Accumulator (from low pressure default parameters)
V_0 = 6;                                                           % Initial volume [m^3]
p_rated = 16e6;                                                    % Rated working pressure (Pa)
p_upper_limit = (4/3)*p_rated;               % Upper working pressure
p_lower_limit = (0.5)*p_upper_limit;         % Lower working pressure
p_precharge = 0.9*p_lower_limit;             % Precharge pressure
V_max = V_0*(1-(p_precharge/p_upper_limit)^(1/1.4)); % Volume corresponding to upper working pressure
V_min = V_0*(1-(p_precharge/p_lower_limit)^(1/1.4)); % Volume corresponding to lower working pressure
V_eq = V_max;                           % Start at the volume corresponding to upper working pressure
p_eq = p_precharge/(1-V_eq/V_0)^(1.4);  % Start at the upper working pressure
V_rated = V_0*(1-(p_precharge/p_rated)^(1/1.4)); % Volume corresponding to rated pressure


%% Desired Shaft Speed
%n_shaft_ref_m = 130; % rev/s 
n_shaft_ref_m = 50;%J_w_des*A_mem*n_ser*n_par/(V_D_hpp); % rev/s (near 50 rev/s)

%% Desired Recovery Ratio and Flux
RR_tot = 0.5;                      %Final desired recovery ratio (m^3 perm / m^3 feed)
%J_w_des = 4.17e-6;     %15 lmh        %Desired flux (m/s)   % LMH /3.6e6 = m/s
J_w_des = 8.33e-06;

%% Pump Efficiencies
eta_hp = 0.85;                     %mechanical efficiency of HPP
eta_cp = 0.65;                     %mechanical efficiency of CP

%% Membrane Sizing Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Q_in_avg = .014;                % Average input flowrate (m^3).  Measured after simulation, input here.
n_par = 320;                        %Number of membranes in parallel
n_ser = 1;                         %Number of membranes in series

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thick = 7.112e-4;                         %Spacer thickness [m] % thick = %0.028; 
L_mem = .96;                      %Length of membrane [m] % L_mem = 1.016; 
A_mem = 7.4;                        %Area of membrane [m^2] % A_mem = 40.88;
V_mem = A_mem * thick / 2;                         %Volume of membrane [m^3]

%% Desired Flow Rates
Q_p_des = J_w_des * A_mem * n_par * n_ser;                  %Desired permeate flow rate [m^3/s]
RR_inst_des = 0.1;                          %Desired instantaneous recovery ratio (m^3/s perm / m^3/s feed)
Q_cp = Q_p_des / RR_inst_des - Q_p_des;     %Circulation pump flow rate (fixed) [m^3/s]


%% Hi-P Pump Displacement Volume
%V_D_hpp = 2.27e-6 * 10 * 10 / 10;                           %volumetric displacement of HPP [m^3/rev]
%V_D_hpp = 5.98e-4; %m^3/rev
V_D_hpp = Q_p_des/n_shaft_ref_m ;      %Volumetric displacement of HPP [m^3/rev]

%% Motor Displacement Volume 
desPropThruMain = .8;  %Desired proportion of flowrate thru main loop (0.8 provides safety against a 0 flowrate thru kidney loop)
Q_2_des = Q_in_avg*desPropThruMain; %Desired flowrate thru main loop (m^3/s)
%V_D_m = 3.25e-6 / 2.5 * 10 * 10 / 10; %m^3/rev
%V_D_m = 4.42e-4; %m^3/rev
V_D_m = Q_2_des/n_shaft_ref_m;  %Volumetric displacement of hydraulic motor [m^3/rev]

%% Variable Displacement Cylinder Dimensions
L_tank = 1;  %m                                        %Length (stroke) of VDC [m]
V_tank_i = (V_mem * n_ser * n_par * RR_tot)/(1 - RR_tot);   %Piston tank volume [m^3]
A_p = V_tank_i / L_tank; %Area of VDC piston (m^2)A_p = 0.0081; %L_tank = V_tank_i / A_p;

%% Constants: Salinity Dynamics 1
A_w = 5.56e-12;                              %Membrane water permeability coefficient [m/s-Pa] 
D_h = 2*thick;                               %Hydraulic diameter for a narrow channel [m]
nu = 8.56e-7;                                %Kinematic viscosity of water [m^2/s]
D = 1.47e-9;                                 %Diffusion coefficient of salt in water https://doi.org/10.1021/ja01589a011 m/s
Sc = nu/D;                                   %Schmidt number for water
C_f_i = 35;                                  %Initial feed salinity [g/kg] %C_f_i = 5;  

M_NaCl = 58.55;                              %Molar mass of NaCl [kg/kmol]
rho_f = 1025;                                %feed density [kg/m^3]
R = 8.314;                                   %Universal gas constant [kJ/kmol-K]
T = 300;                                     %Temperature [K]
i = 2;                                       %van't Hoff factor (2 for a monovalent binary salt)
Phi = 0.93;                                  %Osmotic coefficient

b_prime = i * R * T / M_NaCl * Phi;          %kJ/kg salt
b = b_prime / 1e3 * 1e3 * rho_f;       %factor for calculating osmotic pressure (Pa * kg water / kg salt) 

%% Hydraulic Motor
J = 20*100;                               % Shaft (motor & generator) rotational moment of inertia  [kg-m^2]
angVelInit = n_shaft_ref_m;                         % Initial shaft speed (rev/s)
eta_motor = .95;  % Motor efficiency (from PTO-Sim paper https://energy.sandia.gov/wp-content/uploads/2014/06/SAND2015-2069C.pdf)
% ditto in https://www.nrel.gov/docs/fy19osti/71078.pdf

p_c_init = p_rated/4;  % Pa, based on viewed results
p_h_init = p_rated;

%% Control
Q_2_init = V_D_m*n_shaft_ref_m;
Omega_2_init = Q_2_init/(0.7*sqrt(2*p_c_init/rho_f));
Q_1_init = Q_in_avg - Q_2_init;
Omega_1_init = Q_1_init/(0.7*sqrt(2*p_h_init/rho_f)); 

% MOTOR VALVE (V2)
Kp = 1e-5;
% Kp = 1e-8*n_shaft_ref; % (m^2)/(rev/s)
Ki = 0;
Kd = Kp*100;
% Kd = 0;

% KIDNEY VALVE (V1)
% Omega_1_init = 3.5e-6; %% UNKNOWN how to get; changes based on load. (m^2)
% Kp_k = 0;
Kp_k = (-1e-13)/1.5; % (m^2)/(Pa)
% Ki_k = -1e-11;
Ki_k = 0;
% Kd_k = 0;
Kd_k = Kp_k*100;
% Kd_k = -1e-11;

tFactor = 1;
num_batches = 10*tFactor; % BRO cycles to simulate