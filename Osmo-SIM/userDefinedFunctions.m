%% Housekeeping / Diagnostic
close all
clear table  
%set(0,'DefaultFigureWindowStyle','docked')

Q_in_avg_meas = mean(Q_in)
Q_1_avg_meas = mean(Q_1);
Q_2_avg_meas = mean(Q_2);
Q2_Qin = Q_2_avg_meas/Q_in_avg_meas
Q1_min = min(Q_1);

n_mem_max = Q_in_avg_meas*V_D_hpp/(J_w_des*A_mem*V_D_m);
n_mem_est = n_mem_max * .7;
P_kvalve_min = min(P_kvalve_flow)/1000 %kW
P_mvalve_min = min(P_mvalve_flow)/1000 %kW


%Q_in_avg_m3_day = mean(Q_in)*86400 %m^3/day
Q_p_avg_m3_day = mean(Q_p)*86400 %m^3/day

t_simu = V_p_tot.time(end);
%% Time tracker
t_active = V_p_tot.time(end); %Total cycle length (all cycles)
V_p_active = V_p_tot.data(end); %Total volume produced (all cycles)

Q_cp_flushing = 10*Q_cp; %flow rate flushing (m^3/s)
t_flushing = V_tank_i/Q_cp_flushing * num_batches; %(s) (all cycles, note num_batches)
t_total = t_active + t_flushing; % total time across all cycles
Q_p_avg = V_p_active/t_total; 
%Q_p_avg_m3_day = Q_p_avg * 86400
daily_cycles = 86400/t_total * num_batches; % cycles per day ("cycle" = "batch")

E_circ_active = E_cp.data(end);
u_avg_flushing = Q_cp_flushing*2*L_mem/(thick*A_mem);                      %Average velocity on feed side [m/s]
Re_flushing = u_avg_flushing * thick / nu;                                      %Reynolds number                     
f_flushing = 3.23*Re_flushing^-0.3;                                        %friction factor
p_drop_flushing = n_ser*f_flushing*L_mem*rho_f*u_avg_flushing^2/(2*D_h);                 %Pressure drop across all membrane modules [Pa]
P_flushing = Q_cp_flushing * p_drop_flushing / eta_cp;
E_circ_flushing = P_flushing * t_flushing; 
E_circ_tot = E_circ_active + E_circ_flushing; 

%% SEC

SEC_hpp = mean(P_hpp_shaft)/(mean(Q_p))/3600/1000 ;%kWh/m^3

SEC = mean(P_WEC_piston)/Q_p_avg/3600/1000;
WECPower = mean(P_WEC_piston.data);
MotorPower = mean(P_motor_flow.data);
Recovery1 = mean(P_mvalve_flow.data);
Recovery2 = mean(P_kvalve_flow.data);
P_Back = (Recovery1+Recovery2)*eta_gen;
SEC_final = (WECPower-P_Back)/Q_p_avg/3600/1000 %net power consumed/permeate produced
NREL = 2.8; 
Interview = 3.18;
x = categorical({'Osmocean-Gen', 'WEC-RO', 'Elec-RO', 'Osmocean'});
x = reordercats(x, {'Osmocean-Gen', 'WEC-RO', 'Elec-RO', 'Osmocean'});
y = [SEC_final NREL Interview SEC];
bar(x,y, 'FaceColor', [0.9290, 0.6940, 0.1250])
ylabel('SEC')
title('Reported vs. Osmocean SEC')


%% LCOW
%CapEx calculations RO
WEC_CapEx = 3877896; 
WEC_OpEx = 68107; 
WEC_CapEx_array = 100*WEC_CapEx;
WEC_OpEx_array = 100*WEC_OpEx;
CapRO = Q_p_avg_m3_day;
CapRO_array = 100*CapRO;
BatchelorsBudget_array = 107238.01*Q_p_avg_m3_day;
capFac = 0.49; 
capexRO = BatchelorsBudget_array; % CapEx budget for 100 m^3/day BRO plant
AWP = CapRO*capFac*365;
AWP_array = AWP*100;

% OpEx calculations RO
Nl = round((CapRO_array*264/6e+6)^0.4 + 18/1.4); % Number of laborers
Nm = round((5+ CapRO_array/5500)/2); % Number of managers
labor = 29700*Nl; % Direct labor costs
manage = 66000*Nm; % Management labor costs
spare = 0.04*AWP_array; % Spare parts
pre = 0.03*AWP_array; % Pretreatment
post = 0.01*AWP_array; % Posttreatment
memb = 1.2*0.07*AWP_array; % Membranes
ins = 0.005*capexRO; % Insurance

CostBack = P_Back*100*8e-5*8760; %W * $/Whr * hr/year

opexRO = labor+ manage + pre + post + memb + spare + ins - CostBack;


LCOW = (0.108*(capexRO+WEC_CapEx_array)+opexRO+WEC_OpEx_array)/AWP_array %$/m^3/yr


%% Plots 

set(0,'DefaultFigureWindowStyle','docked')
set(0, 'DefaultFigureColor', 'w')
time = F_WEC.time;

% %Piston Displacement (redundant)
% figure()
% plot(time, x_tank.data, 'r-', 'LineWidth', 1.04)
% title('Piston Displacement vs. Time')
% xlabel('Time (sec)')
% ylabel('Displacement (m)')
% xlim([0 t_simu])


% %PTO Force (redundant)
% figure()
% plot(time(2000:3000), F_WEC.data(2000:3000)/1e6, 'b-', 'LineWidth', 1.04)
% title('PTO Force vs. Time')
% xlabel('Time (sec)')
% ylabel('Force (MN)')
% xlim([200 300])

%Shaft Speed (redundant)
figure()
plot(time, n.data, 'b-', 'LineWidth', 1.04)
hold on
yline(n_shaft_ref_m, 'r--', 'LineWidth', 1.04)
title('Shaft Speed vs. Time')
xlabel('Time (sec)')
ylabel('Shaft Speed (rps)')
ylim([n_shaft_ref_m-1 n_shaft_ref_m+1])
legend('Actual', 'Desired', 'Location', 'northeast', 'Orientation', 'horizontal', 'Box', 'off');
xlim([0, 400*tFactor])
set(0, 'DefaultFigureColor', 'w')

%piston Pressure
figure()
plot(time, p_h.data/1e6, 'b-', 'LineWidth', 1.04)
hold on
plot(time, p_c.data/1e6, 'r-', 'LineWidth', 1.04)
title('WEC Piston Pressure = Accumulator Pressure vs. Time')
xlabel('Time (sec)')
ylabel('Pressure (MPa)')
ylim([0 25])
xlim([0 400*tFactor])
hold on
yline(21, 'k--')
hold on
yline(11, 'k--')
legend('p_h', 'p_c', 'Location', 'northeast', 'Orientation', 'horizontal', 'Box', 'off');
set(0, 'DefaultFigureColor', 'w')

% Valve Areas
figure()
plot(time, Omega_1.data, 'b-', 'LineWidth', 1.04)
hold on
plot(time, Omega_2.data, 'r-', 'LineWidth', 1.04)
title('Valve Area vs. Time')
xlabel('Time (sec)')
ylabel('Valve Area (m^2)')
ylim([0 inf])
xlim([0 400*tFactor])
legend('Kidney Valve Area','Motor Valve Area', 'Location', 'northeast', 'Orientation', 'vertical', 'Box', 'off');

%BRO Pressures
figure()
plot(time, p_f.data/1e5, 'b-', 'LineWidth', 1.04)
% hold on
% plot(time, p_b.data/1e5, 'r-', 'LineWidth', 1.04)
hold on
plot(time, pi.data/1e5, 'r-', 'LineWidth', 1.04)
title('BRO pressures vs. Time')
xlabel('Time (sec)')
ylabel('Pressure (bar)')
legend('Feed pressure and HPP outlet pressure', 'Osmotic pressure', 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
xlim([0 400*tFactor])
set(0, 'DefaultFigureColor', 'w')
%Flow Rates
figure()
plot(time, Q_in.data, 'r-', 'LineWidth', 1.04)
hold on
plot(time, Q_out.data, 'b-', 'LineWidth', 1.04)
hold on
plot(time, Q_1.data, 'm-', 'LineWidth', 1.04)
hold on
plot(time, Q_2.data, 'g-', 'LineWidth', 1.04)
title('Flowrate vs. Time')
xlabel('Time (sec)')
ylabel('Flowrate (m^3/sec)')
legend('Accumulator inlet flow', 'Accumulator outlet flow', 'Kidney loop flow', 'Motor loop flow', 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');
xlim([0 400*tFactor])

% %Power (SHAFT)
% figure()
% plot(time, P_WEC_shaft.data, 'r-', 'LineWidth', 1.04)
% hold on
% plot(time, P_motor_shaft.data, 'b-', 'LineWidth', 1.04)
% hold on
% plot(time, P_hpp_shaft.data, 'g-', 'LineWidth', 1.04)
% title('Shaft Power vs. Time')
% xlabel('Time (sec)')
% ylabel('Shaft Power (W)')
% legend('WEC Input Shaft Power','Motor Shaft Power','BRO Shaft Power', 'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off');

%Power (FLOW)
figure()
plot(time, P_WEC_flow.data/1e6, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1.04)
hold on
plot(time, P_motor_flow.data/1e6, 'r-', 'LineWidth', 1.04)
hold on
plot(time, P_hpp_flow.data/1e6, 'b-', 'LineWidth', 1.04)
hold on
plot(time, P_kvalve_flow.data/1e6, 'm-', 'LineWidth', 1.04)
hold on
plot(time, P_mvalve_flow.data/1e6, 'c-', 'LineWidth', 1.04)
title('Flow Power vs. Time')
xlabel('Time (sec)')
ylabel('Flow Power (MW)')
legend('WEC','Motor','High Pressure Pump','Kidney Valve','Motor Valve', 'Location', 'northoutside', 'Orientation', 'horizontal', 'Box', 'off');
xlim([0 400*tFactor])



set(findall(groot,'type','text'),'fontSize',12,'fontWeight','bold')

% test1 = mean(p_h.data)
% test2 = mean(Q_out.data)


% [Start_pressure_loss, PAT_pressure_loss, Kidney_Pressure_Loss, PAT_pressure, Kidney_pressure] = Mass_Flow_in_a_bow(4, mean(Q_out.data), mean(p_h.data));
% function [Start, PAT, Kidney, PATPres, KidneyPres] = Mass_Flow_in_a_bow(pipe_diam, InletFlowrate, InletPressure)
% clear all
% clc
% clear
%% Packaged Pipe flow calculations for use in larger model - N. Kiefer
% Subfunction version of pipe flow calculations to integrate with model
% inputs:  - Pipe diameter                     [inches ]
%          - pressure at the inlet             [Pa]
%          - Flowrate into the WEC             [m^3/sec]
%
% outputs: - Pressure Loss in inlet            [Pa]
%          - Pressure Loss in Turbine Loop     [Pa]
%          - Pressure Loss in Kidney Loop      [Pa]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Code Starts Here
% pressure_psi = pa2psi(InletPressure);
% Q_start_gpm = m3psec2gpm(InletFlowrate);
% %pressure_psi = 100:1:200; %Psi at the inlet from the WEC
% %Q_start = 1:1:500; %GPM from the WEC
% [Start, PAT, Kidney, PATPres, KidneyPres] = MassFlow(pressure_psi, Q_start_gpm, pipe_diam); %Pressure losses in PSI
% [Start, PAT, Kidney, PATPres, KidneyPres] = MassFlow(pressure_psi, Q_start_gpm, pipe_diam); %Pressure losses in PSI
% conversionConst_psi = 1 / 6894.76; %psi / Pa
% conversionConst_gpm = 1 / 15850.32314; % m^3/sec // gpm
% 
% Start = Start ./ conversionConst_psi;
% PAT = PAT ./ conversionConst_psi;
% Kidney= Kidney ./ conversionConst_psi;
% PATPres= PATPres ./ conversionConst_psi;
% KidneyPres= KidneyPres ./ conversionConst_psi;
% 
% function [pressure_bar] = pa2psi(InletPressure)
%     conversionConst = 1 / 6894.76; %psi / Pa
% 	pressure_bar = InletPressure .* conversionConst;
% end
% function [Q_start_gpm] = m3psec2gpm(Flowrate)
%     conversionConst = 1 / 15850.32314; % m^3/sec // gpm
%     Q_start_gpm = Flowrate ./ conversionConst;
% end
% function [H_start_loss, H_turb_loss, H_generator_loss, p_turbine, p_generator] = MassFlow(p_in, Q_start, pipe_diam);
% %% Description
% % This script, written by N. Kiefer for the Batchlors senior design project
% % is designed to take input pressure, density, and mass flow rate through
% % the ERD turbine and calculate mass flow rates and losses in the system.
% %
% % NOTE: ALL FITTING LOSS COEFFICIENTS ARE FROM Fox pg 307 unless otherwise stated
% %
% % INPUTS : input pressure, ERD flow rate, flow distribution, losses
% % OUTPUTS: Pressure Loss calcuations, local flow rates
% %% Inputs
% % Givens
% %p_in = 151; %psi
% %Q_turb = 11273; %gpm
% g = 9.81; %m/s^2
% meter_to_psi = 1 / 1.42197; % m_head / psi
% psi2bar = 0.0689476; %bar to psi
% % Assumptions / Found Values
% SG_water = 1.025; %dimless
% rho_water_reference = 1000; %kg/m^3
% roughness_specific = 0.0197E-3; %ft - SS, Turned - https://www.engineeringtoolbox.com/surface-roughness-ventilation-ducts-d_209.html
% mu_metric = 1.01E-3; % N*s / m^2
% % System Dimensions
% nip_length = 2 / 12; % nipple length in ft
% Start_length = 4 * nip_length; %ft - feet of tubing under pressure in section A
% Gener_length = 4 * nip_length; %ft - feet of tubing under pressure in section B
% Turb_length = 4 * nip_length; %ft - feet of tubing under pressure in section C
% %pipe_diam = 4:.25:12; %inch
% %% Preperation Calculations
% rho_water = calculateDensity(SG_water, rho_water_reference);
% %% Calculating Mass and Volumetric Flow Rates
% m_dot_start = massflow(rho_water, Q_start, "mass"); %kg/m^3
% m_dot_turbine = m_dot_start .* 0.9; %kg/m^3
% m_dot_generator = m_dot_start .* .1; %kg/m^3
% Q_turb = massflow(rho_water, m_dot_turbine, "volume"); %gpm
% Q_generator = massflow(rho_water, m_dot_generator, "volume"); %gpm
% %% Calculate pipe surface roughness
% pipe_diam_ft = pipe_diam ./ 12;
% % fprintf("The relative surface roughness is \n")
% relative_roughess = roughness_specific ./ pipe_diam_ft; %dimless
% %% Calculating pipe area in m^2
% pipe_diam_meter = pipe_diam ./ 39.3701; %m
% pipe_area_metric = 3.14 .* ((pipe_diam_meter ./ 2) .^2);
% %% Calculating pipe velocities
% velo_turb      = veloCalc(rho_water, pipe_area_metric, m_dot_turbine); %m/s
% velo_start     = veloCalc(rho_water, pipe_area_metric, m_dot_start); %m/s
% velo_generator = veloCalc(rho_water, pipe_area_metric, m_dot_generator);  %m/s
% %% Calculating Reynods number
% Re_turb      = Reynolds(rho_water, pipe_diam_meter, velo_turb, mu_metric); %dimless
% Re_start     = Reynolds(rho_water, pipe_diam_meter, velo_start, mu_metric); %dimless
% Re_generator = Reynolds(rho_water, pipe_diam_meter, velo_generator, mu_metric); %dimless
% %% FROM MOODY DIAGRAM FIND Friction Factor = 0.07
% f = 0.019;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Calculation of head loss in each section
% % Total flow loss
% [h_start, H_start] = pipeLoss(f, Start_length, pipe_diam_meter, velo_start, g); %metric
% [h_turb, H_turb] = pipeLoss(f, Turb_length, pipe_diam_meter, velo_turb, g); %metric
% [h_generator, H_generator] = pipeLoss(f, Gener_length, pipe_diam_meter, velo_generator, g); %metric
% % Minor Losses
% tee_thread_branch = 2.0;
% check_valve_swing = 2.30;
% % Start Fittings Losses
% % Addtl accounts for M2M thread addaptor
% K_start = 4 * tee_thread_branch + 1.0;
% [h_start_minor, H_start_minor] = minorLoss(K_start, velo_start, g);
% % Generator Fittings Losses
% %
% K_generator = 2 * tee_thread_branch + 1 * check_valve_swing;
% [h_generator_minor, H_generator_minor] = minorLoss(K_generator, velo_generator, g);
% % Turbine Fittings Losses
% %
% K_turb = 2 * tee_thread_branch + 1 * check_valve_swing;
% [h_turb_minor, H_turb_minor] = minorLoss(K_turb, velo_turb, g);
% %% Summing up Pressure Losses
% H_start_total = H_start + H_turb_minor;
% H_generator_total = H_generator + H_generator_minor;
% H_turb_total = H_turb + H_turb_minor;
% h_start_total = h_start + h_turb_minor;
% h_generator_total = h_generator + h_generator_minor;
% h_turb_total = h_turb + h_turb_minor;
% %% Making Sense of Pressure Loss
% H_start_loss = m2psi(H_start_total, meter_to_psi);
% H_turb_loss = m2psi(H_turb_total, meter_to_psi);
% H_generator_loss = m2psi(H_generator_total, meter_to_psi);
% % Calculating Pressure At generator and tubine
% p_generator = p_in - H_start_loss - H_generator_loss;
% p_turbine = p_in - H_start_loss - H_turb_loss;
%% Outputting to command window
% fprintf("\nWith inlet conditions of %.2f psi (%.3f bar) and %d inch pipes\n", p_in,  psitobar(p_in, psi2bar), pipe_diam)
% fprintf("\nFlow at generator inlet is %.2f psi (%.3f bar) @ %.2f gpm \n", p_generator, psitobar(p_generator, psi2bar), Q_generator)
% fprintf("Flow at turbine inlet is %.2f psi (%.3f bar) @ %.2f gpm \n\n", p_turbine, psitobar(p_turbine, psi2bar), Q_turb)
% fprintf("Pressure losses in:\n     Start - %.2f psi (%.3f bar)\n     Generator - %.2f psi (%.3f bar)\n     Turbine - %.2f psi (%.3f bar)\n\n", H_start_loss, psitobar(H_start_loss, psi2bar),H_generator_loss, psitobar(H_generator_loss, psi2bar),H_turb_loss, psitobar(H_turb_loss, psi2bar))
%% Graphing
% figure(1)
% hold on
% plot(p_in, H_start_loss, 'k')
% plot(p_in, H_turb_loss, 'r')
% plot(p_in, H_generator_loss, 'm')
% %xscale = [.25, 1.5];
% legend('Inlet Loss', 'Hydraulic Generator Loop Loss', 'Kidney Loop Loss')
% xlabel('p_in [psi]');
% ylabel('Pressure Loss [psi]');
% title('Pressure Loss as a Function of p_in');
% %fprintf('Pressure Loss as a Function of Pipe Diameter at %d gpm\n', Q_turb);
%
% set(gca, 'TitleFontSizeMultiplier', 2.5)
%set(gca, 'XTick', [Q_start])
%set(gca, 'YTick', [0:50:max(H_start_loss)+3*50])
%set(gca, 'YScale', 'log')
% % ax.Color = rand(1,3);%[0,0,0]
% set(gca, 'Color', [1,1,1])
% % set(gca, 'figureColor', [0,0,0])
% set(gca, 'XColorMode', 'auto')
% set(gca, 'YColorMode', 'auto')
% %% Subfunctions
% %Subfunction that converts psi to bar
% function output = psitobar(psi, psi2bar)
% output = psi .* psi2bar;
% end
% %Subfunction that converts meter head loss to psi
% function output = m2psi(m, m2p)
%  output = m .* m2p;
% end
% %Subfunction that calculates minor losses due to friction
% function [h, H] = minorLoss(K, V, g)
% h = K .* (V.^2)./2;
% H = h ./ g;
% end
% % Subfunction that calculates pipe head loss due to friction
% function [h, H] = pipeLoss(f, L, D, V, g)
% h = f .* (L ./ D) .* (V.^2) ./ 2;
% H = h ./ g;
% end
% % Subfunction that calcualtes velocity from mass flow rate
% function velo = veloCalc(rho, area, m_dot)
%     velo = m_dot ./ (rho .* area);
% end
% % Subfunction to calculate reynolds number
% function Re = Reynolds(rho, D, V, mu)
%     Re = (rho .* D .* V) ./ mu;
% end
% % Subfunction to calculate values for mass flow converstions
% function output = massflow(alpha, bravo, charlie)
% gal_to_mcube = 1/264; %m^3 / gal
% sec_to_min = 60; %sec / min
% output = "ERROR IN FUCNTION MASSFLOW";
%     if charlie == "mass"
%         output = alpha .* bravo .* gal_to_mcube .* (1 ./ sec_to_min);
%     end
%     if charlie == "volume"
%         output = (bravo ./ alpha) .* sec_to_min ./ gal_to_mcube;
%     end
% end
% % Subfunction to calculate density of water in system
% function rho = calculateDensity(SG, rho_ref)
%     rho = SG .* rho_ref;
% end
% end
% end