% -------------------------------------
% Alia Zaki
% ME EN 4650
% Lab 11 - Blackbody Radiation
% -------------------------------------

clc
clear all
close all

%======================================
%               DATA
%======================================
% CONVERSIONS---------------------------------
mmHgToPa = 133.32;   
inH2OToPa = 248.84; 
mmTom = 1e-3;
gpmTom3_s = 6.3e-5;
inTom = 0.0254;

% DIMENSIONS---------------------------------
Do_tube = 0.25*inTom;    % [m] outside tube diameter 
t = 0.028*inTom;         % [m] tube wall thickness
L = 9*inTom;             % [m] shell/tube length
Do_shell = 2.12*inTom;   % [m] shell outside diameter
s = 1.125*inTom;         % [m] baffle spacing

Atube = pi*(Do_tube - t)*L; % [m^2] inside area of tube 

% EXPERIMENTAL DATA----------------------------
data = xlsread('HeatExchanger_DataSheet.xlsx');
ind = 7:10;
Patm = data(1,1)*mmHgToPa;   % atmospheric pressure [Pa]

V_dot_c = data(ind,1)'*gpmTom3_s;  % [m^3/s] volumetric flow rate (cold)
V_dot_h = data(ind,2)'*gpmTom3_s;   % [m^3/s] volumetric flow rate (hot) 

Th_i = (data(ind,3) - 32)' * (5/9) + 273.15; % [K] inlet temperature (hot)
Th_o = (data(ind,4) - 32)' * (5/9) + 273.15; % [K] oulet temperature (hot) 

Tc_i = (data(ind,5) - 32)' * (5/9) + 273.15; % [K] inlet temperature (cold)
Tc_o = (data(ind,6) - 32)' * (5/9) + 273.15; % [K] outlet temperature (cold) 

Tshell = (data(ind,7) - 32)' * (5/9) + 273.15; % [K] shell casing temperature
Tamb = (data(ind,8) - 32)' * (5/9) + 273.15;   % [K] ambient temperature

%======================================
%             ANALYSIS
%======================================
Tc_ave = 0.5*(Tc_i + Tc_o); % [K] average cold temperature
Th_ave = 0.5*(Th_i + Th_o); % [K] average hot temperature

% rho [kg/m^3], cp [J/kg-K]
[rho_c,cp_c] = WaterProperties(Tc_ave-273.15); 
[rho_h,cp_h] = WaterProperties(Th_ave-273.15);

% mass flow rate [kg/s]
m_dot_c = rho_c.*V_dot_c;
m_dot_h = rho_h.*V_dot_h;

% heat capacities 
Cc = m_dot_c.*cp_c;
Ch = m_dot_h.*cp_h;

C = [Cc; Ch];

Cmin = min(C);
Cmax = max(C);

Cr = Cmin./Cmax;

% heat transfer rate [J/m^2-kg]
qh = Ch.*(Th_i - Th_o); 
qc = Cc.*(Tc_o - Tc_i);

eff = qh./(Cmin.*(Th_i - Tc_i));

delta_T1 = Th_i - Tc_o;
delta_T2 = Th_o - Tc_i;
delta_Tlm = (delta_T1 - delta_T2)./log(delta_T1/delta_T2);
%delta_Tlm = abs(delta_Tlm);
UA = qh/delta_Tlm; 

NTU = UA./Cmin;

eff_th_num = 1 - exp(-NTU.*(1-Cr));
eff_th_den = 1 - Cr.*exp(-NTU.*(1-Cr));

eff_th = eff_th_num./eff_th_den;

% CONVECTION COOLING---------------------
Tf = 0.5*(Tshell + Tamb);           % [K] film temperature 

% air fluid properties
for i = 1:4
    [rho_air(i),miu_air(i),k_air(i),cp_air(i)] = AirProperties(Tf(i),Patm);
    %[rho_air,miu_air,k_air,cp_air] = AirProperties(Tf(i),Patm);
end

nu_air = miu_air/rho_air;           % [m^2/s] kinematic viscosity 
alpha_air = k_air./(rho_air.*cp_air); % [m^2/s] thermal diffusivity
beta = 1./Tf;                        % [1/K] volumetric thermal expansion coefficient
g = 9.81;                           % [m/s^2] gravitational acceleration

% Raleigh Number 
RaD_num = g.*beta.*(Tshell - Tamb).*Do_shell^3;
RaD_den = nu_air.*alpha_air;
RaD = RaD_num/RaD_den;

% Nusselt Number 
NuD = 0.48.*RaD^(1/4);

% average heat transfer coefficient
h = k_air*NuD/Do_shell;

qconv = h*pi*Do_shell*L.*(Tshell - Tamb)

% THERMAL RADIATION----------------------
sigma = 5.6703e-8; % [W/m^2-K^4] Stefan-Boltzmann constant
e = 0.95;          % [-] emissivity of shell

qrad = e*sigma.*(Tshell.^4 - Tamb.^4)*pi*Do_shell*L

%`PERCENT UNCERTAINTY---------------------
delta_Tc = (Tc_o - Tc_i);
delta_Th = (Th_i - Th_o);

sigma_deltaT = 0.1; 
sigma_Vdot = 0.2*gpmTom3_s;

uncertainty_c = 100*sqrt((sigma_Vdot./V_dot_c).^2 + (sigma_deltaT./delta_Tc).^2);
uncertainty_h = 100*sqrt((sigma_Vdot./V_dot_h).^2 + (sigma_deltaT./delta_Th).^2);

% PERCENT DIFFERENCE 
delta_q = 100 * abs(qh - qc)./(0.5*(qh + qc));
delta_eff = 100* abs(eff - eff_th)./eff_th;

%======================================
%             FIGURES
%======================================

% FIGURE 1c------------------------------
p1 = openfig('EffectivenessNTU.fig');
hold on
p2 = plot(NTU,eff,'bo');
p3 = plot(NTU,eff_th,'k+');
p = [p1 p2 p3];
hold off
legend([p2 p3],'Experimental','Theoretical','Location','Southeast');

figure;
plot(NTU,eff,'bo');
hold on
plot(NTU,eff_th,'k+');
legend('Experimental','Theoretical');
ylabel('Effectiveness');
xlabel('NTU');
