clear all
clc
clear
tic
%% Test Inputs
pipe_diam = 5; %pipe Diameter in inches
InletFlowrate = 0.024; %m^3/sec out of accumulator(avg)
InletPressure = 7.9e6:1e2:8.2e6; %Pressure at inlet to hydraulic Motor
[Start, PAT, Kidney, PATPres, KidneyPres] = Mass_Flow_in_a_bow(pipe_diam, InletFlowrate, InletPressure);
%%Plotting
figure(1)
hold on
%plot(InletPressure, InletPressure)
plot(InletPressure*1e-5, KidneyPres*0.0689476, '--k')
plot(InletPressure*1e-5, PATPres*0.0689476, '-.g')
grid on
xlabel("Inlet Pressure [bar]")
ylabel("Local Pressure [bar]")
title("Local Pressure as a function of Inlet Pressure at const flowrate")
legend_labels = ["Kidney Pressure", "PAT Pressure"];
legend(legend_labels, 'Location', 'southeast')
toc
%%The packaged Subfunction
function [Start, PAT, Kidney, PATPres, KidneyPres] = Mass_Flow_in_a_bow(pipe_diam, InletFlowrate, InletPressure)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Code Starts Here
pressure_psi = pa2psi(InletPressure);
Q_start_gpm = m3psec2gpm(InletFlowrate);
%pressure_psi = 100:1:200; %Psi at the inlet from the WEC
%Q_start = 1:1:500; %GPM from the WEC
[Start, PAT, Kidney, PATPres, KidneyPres] = MassFlow(pressure_psi, Q_start_gpm, pipe_diam); %Pressure losses in PSI
function [pressure_bar] = pa2psi(InletPressure)
    conversionConst = 1 / 6894.76; %psi / Pa
	pressure_bar = InletPressure .* conversionConst;
end
function [Q_start_gpm] = m3psec2gpm(Flowrate)
    conversionConst = 1 / 15850.32314; % m^3/sec // gpm
    Q_start_gpm = Flowrate ./ conversionConst;
end
function [H_start_loss, H_turb_loss, H_generator_loss, p_turbine, p_generator] = MassFlow(p_in, Q_start, pipe_diam);
%% Description
% This script, written by N. Kiefer for the Batchlors senior design project
% is designed to take input pressure, density, and mass flow rate through
% the ERD turbine and calculate mass flow rates and losses in the system.
%
% NOTE: ALL FITTING LOSS COEFFICIENTS ARE FROM Fox pg 307 unless otherwise stated
%
% INPUTS : input pressure, ERD flow rate, flow distribution, losses
% OUTPUTS: Pressure Loss calcuations, local flow rates
%% Inputs
% Givens
%p_in = 151; %psi
%Q_turb = 11273; %gpm
g = 9.81; %m/s^2
meter_to_psi = 1 / 1.42197; % m_head / psi
psi2bar = 0.0689476; %bar to psi
% Assumptions / Found Values
SG_water = 1.025; %dimless
rho_water_reference = 1000; %kg/m^3
roughness_specific = 0.0197E-3; %ft - SS, Turned - https://www.engineeringtoolbox.com/surface-roughness-ventilation-ducts-d_209.html
mu_metric = 1.01E-3; % N*s / m^2
% System Dimensions
nip_length = 2 / 12; % nipple length in ft
Start_length = 4 * nip_length; %ft - feet of tubing under pressure in section A
Gener_length = 4 * nip_length; %ft - feet of tubing under pressure in section B
Turb_length = 4 * nip_length; %ft - feet of tubing under pressure in section C
%pipe_diam = 4:.25:12; %inch
%% Preperation Calculations
rho_water = calculateDensity(SG_water, rho_water_reference);
%% Calculating Mass and Volumetric Flow Rates
m_dot_start = massflow(rho_water, Q_start, "mass"); %kg/m^3
m_dot_turbine = m_dot_start .* 0.9; %kg/m^3
m_dot_generator = m_dot_start .* .1; %kg/m^3
Q_turb = massflow(rho_water, m_dot_turbine, "volume"); %gpm
Q_generator = massflow(rho_water, m_dot_generator, "volume"); %gpm
%% Calculate pipe surface roughness
pipe_diam_ft = pipe_diam ./ 12;
% fprintf("The relative surface roughness is \n")
relative_roughess = roughness_specific ./ pipe_diam_ft; %dimless
%% Calculating pipe area in m^2
pipe_diam_meter = pipe_diam ./ 39.3701; %m
pipe_area_metric = 3.14 .* ((pipe_diam_meter ./ 2) .^2);
%% Calculating pipe velocities
velo_turb      = veloCalc(rho_water, pipe_area_metric, m_dot_turbine); %m/s
velo_start     = veloCalc(rho_water, pipe_area_metric, m_dot_start); %m/s
velo_generator = veloCalc(rho_water, pipe_area_metric, m_dot_generator);  %m/s
%% Calculating Reynods number
Re_turb      = Reynolds(rho_water, pipe_diam_meter, velo_turb, mu_metric); %dimless
Re_start     = Reynolds(rho_water, pipe_diam_meter, velo_start, mu_metric); %dimless
Re_generator = Reynolds(rho_water, pipe_diam_meter, velo_generator, mu_metric); %dimless
%% FROM MOODY DIAGRAM FIND Friction Factor = 0.07
f = 0.019;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of head loss in each section
% Total flow loss
[h_start, H_start] = pipeLoss(f, Start_length, pipe_diam_meter, velo_start, g); %metric
[h_turb, H_turb] = pipeLoss(f, Turb_length, pipe_diam_meter, velo_turb, g); %metric
[h_generator, H_generator] = pipeLoss(f, Gener_length, pipe_diam_meter, velo_generator, g); %metric
% Minor Losses
tee_thread_branch = 2.0;
check_valve_swing = 2.30;
% Start Fittings Losses
% Addtl accounts for M2M thread addaptor
K_start = 4 * tee_thread_branch + 1.0;
[h_start_minor, H_start_minor] = minorLoss(K_start, velo_start, g);
% Generator Fittings Losses
%
K_generator = 2 * tee_thread_branch + 1 * check_valve_swing;
[h_generator_minor, H_generator_minor] = minorLoss(K_generator, velo_generator, g);
% Turbine Fittings Losses
%
K_turb = 2 * tee_thread_branch + 1 * check_valve_swing;
[h_turb_minor, H_turb_minor] = minorLoss(K_turb, velo_turb, g);
%% Summing up Pressure Losses
H_start_total = H_start + H_turb_minor;
H_generator_total = H_generator + H_generator_minor;
H_turb_total = H_turb + H_turb_minor;
h_start_total = h_start + h_turb_minor;
h_generator_total = h_generator + h_generator_minor;
h_turb_total = h_turb + h_turb_minor;
%% Making Sense of Pressure Loss
H_start_loss = m2psi(H_start_total, meter_to_psi);
H_turb_loss = m2psi(H_turb_total, meter_to_psi);
H_generator_loss = m2psi(H_generator_total, meter_to_psi);
% Calculating Pressure At generator and tubine
p_generator = p_in - H_start_loss - H_generator_loss;
p_turbine = p_in - H_start_loss - H_turb_loss;
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
%% Subfunctions
%Subfunction that converts psi to bar
function output = psitobar(psi, psi2bar)
output = psi .* psi2bar;
end
%Subfunction that converts meter head loss to psi
function output = m2psi(m, m2p)
 output = m .* m2p;
end
%Subfunction that calculates minor losses due to friction
function [h, H] = minorLoss(K, V, g)
h = K .* (V.^2)./2;
H = h ./ g;
end
% Subfunction that calculates pipe head loss due to friction
function [h, H] = pipeLoss(f, L, D, V, g)
h = f .* (L ./ D) .* (V.^2) ./ 2;
H = h ./ g;
end
% Subfunction that calcualtes velocity from mass flow rate
function velo = veloCalc(rho, area, m_dot)
    velo = m_dot ./ (rho .* area);
end
% Subfunction to calculate reynolds number
function Re = Reynolds(rho, D, V, mu)
    Re = (rho .* D .* V) ./ mu;
end
% Subfunction to calculate values for mass flow converstions
function output = massflow(alpha, bravo, charlie)
gal_to_mcube = 1/264; %m^3 / gal
sec_to_min = 60; %sec / min
output = "ERROR IN FUCNTION MASSFLOW";
    if charlie == "mass"
        output = alpha .* bravo .* gal_to_mcube .* (1 ./ sec_to_min);
    end
    if charlie == "volume"
        output = (bravo ./ alpha) .* sec_to_min ./ gal_to_mcube;
    end
end
% Subfunction to calculate density of water in system
function rho = calculateDensity(SG, rho_ref)
    rho = SG .* rho_ref;
end
end
end