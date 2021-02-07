%-----------------------------
% ME EN 4650
% Lab 2: Spark Ignition Engine
% Alia Zaki
% u1072443
%-----------------------------

close all
clear all
clc

% load data
load('Lab2_data');
air = load('AirProperties.mat');

% Figure 1a
w = data(3,1:5); % crankshaft speed (rpm)
torque = data(4,1:5); % torque (N-m)

plot(w,torque,'ko');
title('1a')
ylabel('Crankshaft Torque (N-m)');
xlabel('Crankshaft Speed (RPM)');

% Figure 1b
w_d = data(12,1:5); % crankshaft speed of dry run test(rpm)
torque_d = data(13,1:5); % torque of dry run test (N-m)

rpm2rads = 2*pi/60; % convert rpm to rad/s

W_dot_f = w_d*rpm2rads.*torque_d; % mechanical losses (Watts)
W_dot_b = w*rpm2rads.*torque; % brake power (Watts)

% to calculate the net power of the system, assume gas is ideal
% calculating Q_dot_out
m_dot_air = data(6,1:5)/60000; % mass flow rate of air (m^3/s)
T_air = data(7,1:5); % inlet air temperature (oC) - T1
T_exhaust = data(9,1:5); % exhaust air temperature (oC) - T4
Cv = 1005; % specific heat at constant volume at room temperature (J/kg*K)

Q_dot_out = m_dot_air.*Cv.*(T_exhaust - T_air);

% calculating Q_dot_in
rho_fuel = 726; % density of 91 octane gasoline (kg/m^3)
LHV = 44e6; % Lower Heating Value of 91 octane gasoline
V_fuel = data(5,1:5)*1e-6; % fuel consumption (m^3)
V_dot_fuel = V_fuel/60; % volumetric flow rate (m^3/s)
m_dot_fuel = rho_fuel*V_dot_fuel; % mass flow rate (kg/s)

E_dot_in = m_dot_fuel*LHV; 

W_dot_net = abs(E_dot_in - Q_dot_out); % net work of engine (Watts)

figure;
plot(w,W_dot_b/1000,'bo');
hold on 
plot(w,W_dot_f/1000,'rs');
plot(w,W_dot_net/1000,'k-');
hold off
title('1b');
ylabel('Power (kW)');
xlabel('Crankshaft Speed (RPM)');
legend('Brake Power','Mechanical Losses','Net Engine Power');

% Figure 1c
nth = W_dot_b./E_dot_in; % thermal efficiency
nth = nth*100; % thermal efficiency in %

k = 1.4; % ratio of specific heats (Cp/Cv)
r = 7; % compression ratio (V_BDC/V_TDC)
n_otto = 1 - (1/(r^(k-1))); % ideal thermal efficiency of the otto cycle
n_otto = n_otto*100; % ideal thermal efficiency in %
n_otto = n_otto*ones(length(nth));

figure;
plot(w,nth,'ko');
hold on
plot(w,n_otto,'k');
hold off
title('1c');
ylabel('Efficiency (%)');
xlabel('Crankshaft Speed (RPM)');
legend('Thermal Efficiency','Ideal Thermal Efficiency (Otto)');

% Figure 1d
Q_dot = E_dot_in - (W_dot_b + W_dot_f); % heat lost to the surroundings (W)

figure;
plot(w,E_dot_in./1000,'bo');
hold on
plot(w,W_dot_b./1000,'gd');
plot(w,Q_dot./1000,'rs');
plot(w,W_dot_f./1000,'kx');
hold off
title('1d');
ylabel('Work Rate (kW)');
xlabel('Crankshaft Speed (RPM)');
legend('Input Energy Rate','Brake Power','Heat Loss Rate','Mechanical Loss Rate');


% Figure 1e
n_c = 2; % number of revolutions per power stroke
B = 65.1/1000; % bore
S = 44.4/1000; % stroke
Vd = (pi/4)*B^2*S; % displacement volume of piston
MEP = torque*2*pi*n_c/Vd; % mean effective pressure 

figure;
plot(w,MEP/1000,'ko');
hold on 
plot(w(4),MEP(4)/1000,'k+');
title('1e');
ylabel('MEP: Mean Effective Pressure (kPa)');
xlabel('Crankshaft Engine Speed (RPM)');
legend('MEP','MEP at 2600 RPM');

% Question 2a
i = mean(W_dot_b./E_dot_in)
ii = mean(W_dot_f./E_dot_in)
iii = mean(Q_dot./E_dot_in)

% Question 2b
n_mech = W_dot_b./(W_dot_b + W_dot_f);
n_mech_ave = mean(n_mech)

discrepancy = (mean(n_otto) - mean(nth))/mean(n_otto)

% Question 2c
MEP_2600 = MEP(4); % MEP at 2600 RPM
A = (pi/4)*B^2; % area of piston
F = MEP_2600*A  % average force


















