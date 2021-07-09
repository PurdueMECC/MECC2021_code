%% Simulation Data
simu = simulationClass();               
simu.simMechanicsFile = 'Merged_Model.slx';  % Specify Simulink Model File with PTO-Sim                  
simu.startTime = 0;                     
simu.rampTime = 1; %seconds                      
simu.endTime=5000;                       
simu.dt = 0.1;                         
simu.CITime = 30;
simu.explorer = 'off';

%% Wave Information
%Irregular Waves using PM Spectrum
waves = waveClass('irregular'); %always do REGULAR FIRST to determine the number of modules
waves.H = 1.25; %m
waves.T = 7.25; %sec
waves.spectrumType = 'PM';
waves.phaseSeed=1;
eta_gen = 0.85;
% eta_gen = 0;

%% Body Data
% Flap
body(1) = bodyClass('../hydroData/oswec.h5');   
body(1).geometryFile = '../geometry/flap.stl';    
body(1).mass = 127000;                         
body(1).momOfInertia = [1.85e6 1.85e6 1.85e6]; 
body(1).linearDamping(5,5) = 1*10^7;    % Specify damping on body 1 in pich

% Base
body(2) = bodyClass('../hydroData/oswec.h5');   
body(2).geometryFile = '../geometry/base.stl';    
body(2).mass = 'fixed';                        

%% PTO and Constraint Parameters
% Fixed Constraint
constraint(1)= constraintClass('Constraint1');  
constraint(1).loc = [0 0 -10];                  

% Rotational PTO
pto(1) = ptoClass('PTO1');                      % Initialize ptoClass for PTO1
pto(1).k = 0;                                   % PTO Stiffness Coeff [Nm/rad]
pto(1).c = 0;                                   % PTO Damping Coeff [Nsm/rad]
pto(1).loc = [0 0 -8.9];                        % PTO Location [m]
