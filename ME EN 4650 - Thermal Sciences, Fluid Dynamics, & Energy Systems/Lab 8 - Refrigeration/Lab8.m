% -------------------------------------
% Alia Zaki
% u1072443
% ME EN 4650
% Lab 8: Vapor-Compression Refrigeration Cycle
% -------------------------------------

clc
clear all
close all

%======================================
%              DATA
%======================================
data = xlsread('RefrigerationLab_SampleData-2');

% CONVERSIONS---------------------------------------
mmHg2Pa = 133.32;
psig2Pa = 6894.76;

% AMBIENT-------------------------------------------
Tamb = data(1,6);              % [deg C]
Pamb = data(2,6)*mmHg2Pa;      % [Pa]
W_dot_fan = data(3,6);         % electrical power to fans [W]

V_dot = data(8:11,1)*6.3e-5; % volumetric flow rate [m^3/s]
W_dot_total = data(8:11,2);  % total power [W]

% REFRIGERATION STATE POINTS------------------------
% temperatures [deg F]
T1 = data(8:11,3);
T2 = data(8:11,5);
T3 = data(8:11,7);  % condensor outlet/expansion valve inlet
T4 = data(8:11,9);  % expansion valve outlet/evaporator inlet
T5 = data(8:11,11);

T = [T1 T2 T3 T4 T5];
T = (T - 32)* 5/9; % [deg C]

% pressure
P1 = data(8:11,4);
P2 = data(8:11,6);
P3 = data(8:11,8);
P4 = data(8:11,10);
P5 = data(8:11,12);

P = [P1 P2 P3 P4 P5]; % [psig]
P = P*psig2Pa + Pamb; % gauge pressure - > absolute pressure [Pa]

% AIR TEMPERATURES-----------------------------------
% exiting condenser [deg C]
Tc1 = data(8:11,13);
Tc2 = data(8:11,14);

Tc = 0.5*(Tc1 + Tc2);

% exiting evaporator [deg C]
Te1 = data(8:11,15);
Te2 = data(8:11,16);

Te = 0.5*(Te1 + Te2);

%======================================
%           DATA ANALYSIS
%======================================
[v,e] = pyversion; system([e,' -m pip install --user -U CoolProp']); % CoolProp install

 T = T + 273.15; % [K]
for i = 1:4
    for j = 1:5
        rho(i,j) = py.CoolProp.CoolProp.PropsSI('D','T',T(i,j),'P',P(i,j),'R134a'); % [kg/m^3]
        h(i,j) = py.CoolProp.CoolProp.PropsSI('H','T',T(i,j),'P',P(i,j),'R134a'); % [J/kg]
    end
end
T = T - 273.15; % [deg C]
h(:,4) = h(:,3);
 % h4 = h3, otherwise need another propery for h4 [J/kg]
m_dot = V_dot.*rho(:,3);           % mass flow rate [kg/s]
W_dot_c = W_dot_total - W_dot_fan; % power supplied to compressor motor

n = 0.78; % nM*nH (mechanical + hydraulic)
W_dot_in = n*W_dot_c; % power supplied to refrigerant by compressor [W] 

w_in = W_dot_in./m_dot;  % specific work to refrigerant [J/kg]

qH = -(h(:,3) - h(:,2)); % heat per unit mass rejected [J/kg]
qL = h(:,1) - h(:,4);   % heat per unit mass transferred [J/kg]
q_loss = w_in + h(:,1) - h(:,2); % heat loss [J/kg]

COP_R = qL./w_in; % coefficient of performance [-]

T = T + 273.15; % [K]
for i = 1:4
    s1_ideal(i) = py.CoolProp.CoolProp.PropsSI('S','T',T(i,1),'P',P(i,1),'R134a');
end
T = T - 273.15; % [deg C]

s2_ideal = s1_ideal; % because ideal cycle [J/K]

T = T + 273.15;
for i = i:4
    h2_ideal(i) = py.CoolProp.CoolProp.PropsSI('H','T',T(i,2),'S',s2_ideal(i),'R134a');
end
T = T - 273.15;

nC_num = h2_ideal(4) - h(:,1);
nC_den = h(:,2) - h(:,1) + q_loss;
nC = (nC_num./nC_den).*100; % isentropic efficiency of compressor

%======================================
%              FIGURES
%======================================
% FIGURE 1a --------------------------------------------
yline(Tamb,'k--');
hold on
plot(m_dot,T(:,3),'bo');
plot(m_dot,T(:,4),'rs');
plot(m_dot,Tc,'bo','MarkerFaceColor','b');
plot(m_dot,Te,'rs','MarkerFaceColor','r');
hold off
title('Figure 1a');
ylabel('Temperature [C]');
xlabel('Mass Flow Rate [kg/s]');
legend('Tamb','T3','T4','Tc','Te','Location','Southeast');

% FIGURE 1b ----------------------------------------------
qL = qL/1000;         % [kJ/kg]
qH = qH/1000;         % [kJ/kg]
q_loss = q_loss/1000; % [kJ/kg]
w_in = w_in/1000;     % [kJ/kg]

figure;
plot(m_dot,qL,'bo');
hold on 
plot(m_dot,qH,'rs');
plot(m_dot,q_loss,'gx');
plot(m_dot,w_in,'kd');
title('Figure 1b');
ylabel('Specific Energy [kJ/kg]');
xlabel('Mass Flow Rate [kg/s]');
legend('qL','qH','q_{loss}','w_{in}');

% FIGURE 1c ----------------------------------------------
figure;
plot(m_dot,COP_R,'md','MarkerFaceColor','m');
title('Figure 1c');
ylabel('COP_R [-]');
xlabel('Mass Flow Rate [kg/s]');

% FIGURE 1d ----------------------------------------------
figure;
yyaxis left
plot(m_dot,nC,'bo','MarkerFaceColor','b');
ylabel('nC [%]');
yyaxis right
plot(m_dot,W_dot_total,'ro','MarkerFaceColor','r');
title('Figure 1d');
ylabel('Total Power [W]');
xlabel('Mass Flow Rate [kg/s]');
legend('nC','Total Power','Location','Southeast');
hold off

% FIGURE 1e ----------------------------------------------
openfig('Ph_Diagram_R134a');
hold on
P = P/1e6; % [MPa]
h = h/1e3; % [kJ/kg]
semilogy(h(1,:),P(1,:),'ro--','MarkerFaceColor','r');
text(h(1,1),P(1,1),'1','Color','k','FontSize',20)
text(h(1,2),P(1,2),'2','Color','k','FontSize',20)
text(h(1,3),P(1,3),'3','Color','k','FontSize',20)
text(h(1,4),P(1,4),'4','Color','k','FontSize',20)
hold off
