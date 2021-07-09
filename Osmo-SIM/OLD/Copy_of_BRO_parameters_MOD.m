close all;
clear;
clc;
%commented HS 3/17/21

%only used for sec
eta_hp = 0.85;                  %high-pressure pump efficiency
eta_cp = 0.65;                  %circulation pump efficiency

RR_pass = 0.1;                     %Recovery ratio per pass [-] %NOT set
RR = 0.5;                           %Final recovery ratio[-] %NOT set
J_w = 20;                            %Membrane flux [LMH] %NOT set
%A_mem = 7.4;                        %RO membrane area [m^2]
A_mem = 30;                         %large-scale RO area (m^2)

n_series = 1;
n_parallel = 90;

V_dot_hpp = J_w*A_mem*1e-3/3600;    %High pressure pump flow rate [m^3/s] %change eqn to calc J_w    %Do we not want to add a multiply by n_para in front?
% assumes that there is one long mem. module; each ind. piece has area A_mem
% flux is defined as avg flux over all pieces

V_dot_cp = V_dot_hpp*(1-RR_pass)/RR_pass; %this WILL BE SET, roughly 10x desired V_dot_perm

%we will change the input parameters for scaling
t = 28*0.0254/1000;             %Spacer thickness [m]
V_mem = A_mem*t/2;              %Feed side RO membrane module volume [m^3]
L_mem = 0.96;                   %RO membrane module length [m]


d_p = 4*0.0254;                 %Piston diameter [m] %Maybe we should scale up diameter too?
A_p = 0.25*pi*d_p^2;            %Piston area [m^2]
L_max = 24*0.0254;              %Piston stroke [m]
V_tank_i = A_p*L_max;           %Piston tank volume [m^3]

L = RR*(L_max + V_mem/A_p);     %Recovery displacement [m] %change to solve for RR

A_w = 2/3600*1e-3;              %Membrane water permeability coefficient [m/s-bar]
D_h = 2*t;                      %Hydraulic diameter for a narrow channel [m]

nu = 8.56e-7;                   %Kinematic viscosity of water [m^2/s]
D = 1.47e-9;                    %Diffusion coefficient of salt in water https://doi.org/10.1021/ja01589a011 m/s
Sc = nu/D;                      %Schmidt number for water

C_f_i = 35;                     %Initial feed salinity [g/kg] roughly g/L
M_NaCl = 58.55;                 %Molar mass of NaCl [kg/kmol]
rho_f = 1025;                   %feed density [kg/m^3]
R = 8.314;                      %Universal gas constant [kJ/kmol-K]
T = 300;                        %Temperature [K]
i = 2;                          %van't Hoff factor (2 for a monovalent binary salt)
b = i*R*T*1e-2/M_NaCl*.93;          %factor for calculating osmotic pressure (bar * g salt / kg water)
%.93 is the osmotic coefficient (phi - see link)

% will need sim parameters to match WEC-Sim
time_max = 0.97*A_p*L_max/V_dot_hpp;
BROsim = sim('BRO_piston_MOD.slx',time_max);

SEC = BROsim.Energy_kJ.data(end)/BROsim.Permeate_m3.data(end)/3600; %kW*hr/m^3
fprintf('run complete\n')

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