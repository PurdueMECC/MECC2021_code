
%% Pre-processing Directives
close all;
clear;
clc;

%% Desired Recovery Ratio and Flux
RR_tot = 0.5;                      %Final desired recovery ratio
J_w_des = 4.17e-6;                  %Desired flux (m/s)

%% Pump Flow Rate Conversion
V_D_hpp = 2.27e-6;                           %volumetric displacement of HPP [m^3/rev]

%% Pump Efficiencies
% eta_volu_hp = 1;                   %volumetric efficiency of HPP (1)
eta_hp = 0.85;                     %mechanical efficiency of HPP
% eta_volu_cp = 1;                   %volumetric efficiency of CP (1)
eta_cp = 0.65;                     %mechanical efficiency of CP

%% Membrane Sizing Parameters
n_par = 2;                         %Number of membranes in parallel
n_ser = 4;                         %Number of membranes in series
thick = 0.028;                         %Spacer thickness [m]
L_mem = 1.016;                         %Length of membrane [m]
A_mem = 40.88;                         %Area of membrane [m^2]
V_mem = A_mem * thick / 2;                         %Volume of membrane [m^3]

%% Desired Flow Rates
Q_p_des = J_w_des * A_mem;                  %Desired permeate flow rate [m^3/s]
RR_inst_des = 0.1;                          %Desired instantaneous recovery ratio
Q_cp = Q_p_des / RR_inst_des - Q_p_des;     %Circulation pump flow rate (fixed) [m^3/s]

%% Variable Displacement Cylinder Dimensions
A_p = 1;                                               %Piston area [m^2]
V_tank_i = (V_mem * n_ser * n_par * RR_tot)/(1 - RR_tot);   %Piston tank volume [m^3]
L_tank = V_tank_i / A_p;                                    %Length (stroke) of VDC [m]

%% Constants: Salinity Dynamics 1
A_w = 5.56e-7;                               %Membrane water permeability coefficient [m/s-bar]  %conv to bar: divide by 1e5
D_h = 2*thick;                               %Hydraulic diameter for a narrow channel [m]
nu = 8.56e-7;                                %Kinematic viscosity of water [m^2/s]
D = 1.47e-9;                                 %Diffusion coefficient of salt in water https://doi.org/10.1021/ja01589a011 m/s
Sc = nu/D;                                   %Schmidt number for water
C_f_i = 35;                                  %Initial feed salinity [g/kg] roughly g/L
M_NaCl = 58.55;                              %Molar mass of NaCl [kg/kmol]
rho_f = 1025;                                %feed density [kg/m^3]
R = 8.314;                                   %Universal gas constant [kJ/kmol-K]
T = 300;                                     %Temperature [K]
i = 2;                                       %van't Hoff factor (2 for a monovalent binary salt)
Phi = 0.93;                                  %Osmotic coefficient

b_prime = i * R * T / M_NaCl * Phi;          %kJ/kg salt
b = b_prime / 1e3 * 1e3 * rho_f / 1e5;       %factor for calculating osmotic pressure (bar * kg water / kg salt) %conv to Pa: multiply by 1e5
 
%% Simulink Running Commands
time_max = 500000;
BROsim = sim('BRO_piston_MOD_test.slx',time_max);

% SEC = BROsim.Energy_kJ.data(end)/BROsim.Permeate_m3.data(end)/3600; %kW*hr/m^3
fprintf('Run complete!\n')

%% Plots
% set(0,'DefaultFigureWindowStyle','docked')
% 
% figure();
% plot(BROsim.power)
% xlim([0,time_max])
% set(findall(gcf,'type','axes'),'fontsize',16)
% xlabel('Time (s)')
% ylabel('Power (kW)')
% title('Hi-P Pump Power Required')
% grid on
% 
% figure();
% plot(BROsim.pressure)
% xlim([0,time_max])
% set(findall(gcf,'type','axes'),'fontsize',16)
% xlabel('Time (s)')
% ylabel('Pressure (bar)')
% title('Absolute Pressure Across RO Membrane')
% grid on
% 
% figure();
% plot(BROsim.pos)
% xlim([0,time_max])
% set(findall(gcf,'type','axes'),'fontsize',16)
% xlabel('Time (s)')
% ylabel('Position (m)')
% title('Position of Piston in Variable Displacement Cylinder')
% grid on
% 
% figure();
% plot(BROsim.flowrate)
% xlim([0,time_max])
% set(findall(gcf,'type','axes'),'fontsize',16)
% xlabel('Time (s)')
% ylabel('Flow Rate (m^3/s)')
% title('Flow Rate Thru Hi P Pump')
% grid on