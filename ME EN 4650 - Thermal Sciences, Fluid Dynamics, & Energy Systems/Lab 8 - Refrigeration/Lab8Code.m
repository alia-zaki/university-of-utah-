%%%%%%%%%%%%%%%%%%%%%
% Jacob Gonzalez-Perez
% u0913317
% ME EN 4610
% Lab 8: Refrigeration Lab
% March 17, 2020
%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;


%% ========================================%%
%                                           %
%              Data Processing              %
%                                           %
% ========================================= %

% Data Read
data = xlsread('RefrigerationLab_SampleData-2.xlsx');

% Load the CoolProp package into Matlab 
[v,e] = pyversion; system([e,' -m pip install --user -U CoolProp']); % CoolProp install

% Experimental Constants
T_amb = data(1,6);                     % Atmospheric temperature [C]
P_amb = Convert('p',data(2,6),'mmHg','Pa');    % Atmospheric pressure from sample data [Pa]
W_fan = data(3,6);                  % Electrical power to fans only [W]
eff = 0.78;

% Data Allocation
volFlow = Convert('mf',data(8:11,1),'gpm','m3ps');   % Volumetric flow [m^3/s]
W_tot = data(8:11,2);

% Temperatures & Pressures
T_1 = data(8:11,3);     % Temperatures in Fahrenheit
P_1 = data(8:11,4);     % Pressures in psig
T_2 = data(8:11,5);
P_2 = data(8:11,6);
T_3 = data(8:11,7);
P_3 = data(8:11,8);
T_4 = data(8:11,9);
P_4 = data(8:11,10);
T_5 = data(8:11,11);
P_5 = data(8:11,12);
T_c1 = data(8:11,13);                           % Temperatures in deg C
T_c2 = data(8:11,14);
T_e1 = data(8:11,15);
T_e2 = data(8:11,16);

% Finding rho and h
T = Convert('t',[T_1,T_2,T_3,T_4,T_5],'F','K');        % [K]
P = Convert('p',[P_1,P_2,P_3,P_4,P_5],'psig','Pa');    % [Pa]
rho = zeros(4,5);
h = zeros(4,5);
for i = 1:4
    for j = 1:5
        rho(i,j) = py.CoolProp.CoolProp.PropsSI('D','T',T(i,j),'P',P(i,j),'R134A');  % [kg/m^3]
        h(i,j) = py.CoolProp.CoolProp.PropsSI('H','T',T(i,j),'P',P(i,j),'R134A');    % [J/kg]
    end
end
h(:,4) = h(:,3);

% Data Analysis
massFlow = volFlow.*rho(:,3);

T_cAvg = 0.5*(T_c1 + T_c2); % [deg C]
T_eAvg = 0.5*(T_e2 + T_e2);

W_comp = W_tot - W_fan;
W_in = eff*W_comp;

w_in = W_in./massFlow;

q_H = -(h(:,3)-h(:,2));
q_L = h(:,1)-h(:,4);

q_loss = w_in + h(:,1) - h(:,2);

COP_R = q_L./w_in;

s_1 = zeros(4,1);
for i = 1:4
    s_1(i) = py.CoolProp.CoolProp.PropsSI('S','T',T(i,1),'P',P(i,1),'R134A');  % [kJ/kg-K]
end

s_2 = s_1;

h_2s = zeros(4,1);
for i = 1:4
    h_2s = py.CoolProp.CoolProp.PropsSI('H','T',T(i,2),'S',s_2(i),'R134A');  % [kJ/kg]
end

eff_c = ((h_2s - h(:,1))./(h(:,2) - h(:,1) + q_loss))*100;



%% ========================================%%
%                                           %
%                Figure 1a                  %
%                                           %
% ========================================= %

% Plotting
yline(T_amb,'k--');
hold on
plot(massFlow,T_3,'bo')
hold on
plot(massFlow,T_4,'rs')
hold on
scatter(massFlow,T_eAvg,'rs','filled')
hold on
scatter(massFlow,T_cAvg,'bo','filled')
hold on

xlabel('Mass Flow [kg/s]')
ylabel('Temperature [deg C]')
legend('T_{amb}','T_3','T_4','Average T_e','Average T_c','Location','southeast')



%% ========================================%%
%                                           %
%                Figure 1b                  %
%                                           %
% ========================================= %

% Plotting
figure
plot(massFlow,q_L/1000,'bo')
hold on
plot(massFlow,q_H/1000,'rs')
hold on
plot(massFlow,q_loss/1000,'gx')
hold on
plot(massFlow,w_in/1000,'kd')

xlabel('Mass Flow [kg/s]')
ylabel('Specific Energy [kJ/kg]')
legend('q_L','q_H','q_{loss}','w_{in}','Location','east')





%% ========================================%%
%                                           %
%                Figure 1c                  %
%                                           %
% ========================================= %

% Plotting
figure
plot(massFlow,COP_R,'ro')

xlabel('Mass Flow [kg/s]')
ylabel('COP_R [-]')



%% ========================================%%
%                                           %
%                Figure 1d                  %
%                                           %
% ========================================= %

% Plotting
figure
yyaxis left
plot(massFlow,eff_c,'ro')
ylabel('\eta_C [%]')

yyaxis right
plot(massFlow,W_tot,'bo')
ylabel('Total Work [W]')

xlabel('Mass Flow Rate [kg/s]')



%% ========================================%%
%                                           %
%                Figure 1e                  %
%                                           %
% ========================================= %

openfig('Ph_Diagram_R134a')
hold on
P = P/1e6; % [MPa]
h = h_1/1e3; % [kJ/kg]
semilogy(h(:,1),P(:,1),'ro--','MarkerFaceColor','r');
