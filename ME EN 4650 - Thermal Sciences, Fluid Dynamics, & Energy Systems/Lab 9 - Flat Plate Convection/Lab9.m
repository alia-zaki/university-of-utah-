% -------------------------------------
% Alia Zaki
% ME EN 4650
% Lab 8: 
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

% FLAT PLATE DIMENSIONS-----------------------
Wtot = mmTom*144; % total width [m]
Ltot = mmTom*297; % total length [m]

w = mmTom*68;     % heated width [m]
Lh = mmTom*153;   % heated length [m]
lead = mmTom*77;  % distance from leading edge [m]
trail = mmTom*67; % distance from trailing edge [m]\

% EXPERIMENTAL----------------------------------------
data = load('ConvectionData.dat');

num = data(:,1);            % thermocouple number [-]
Ts = data(:,2) + 273.15;    % thermocouple temperature [K]

Pbaro = mmHgToPa*(776 - 120.6);    % barometric pressure [Pa]
Tinf = 22.2 + 273.15;              % ambient temperature [K]
rho = AirProperties(Tinf,Pbaro);   % density of air [kg/m^3]

Pdyn = inH2OToPa*0.1;   % dynamic pressure [Pa]
Uinf = sqrt((2*Pdyn)/rho);  % freestream velocity [m/s] 

V = 40.03;        % voltage [VAC]
R = 157.1;        % resistance [Ohms]

%======================================
%             ANALYSIS
%======================================

% Experimental
%---------------
% HEAT FLUX TOP SURFACE-----------------------------
qs_flux = V^2/(2*R*Lh*w); % heat flux [W/m^2]
qs = (w*Lh)*qs_flux;      % heat transfer rate [W] 

% HEAT TRANSFER COEFFICIENT-------------------------
h = qs_flux./(Ts' - Tinf); % local heat transfer coefficient [W/m^2-K]

delta_xi = Lh/14;   % discretized lengths [m]
Lt = Lh - delta_xi; % distance from 1st thermocouple to last [m]

x = [85 92 102 112 123 123 134 143 153 162 173 173 186 196 209 219]*mmTom; % thermocouple location array

h_ave = (1/Lt)*trapz(x,h);  % average heat transfer coefficient [W/m^2-K]

% NUSSELT NUMBER--------------------------------------
Tf = (Ts + Tinf)/2; % local film temperature [K]
Tf_ave = mean(Tf);  % average film temperature [K]

for i = 1:length(Tf)
    [~,~,kf(i),~] = AirProperties(Tf(i),Pbaro); % thermal conductivity of air [W/m^3-K]
end

[~,~,kf_ave,~] = AirProperties(Tf_ave,Pbaro); % thermal conductivity of air
                                              % based on average film temperature

Nu = h.*x./kf;               % local nusselt number [-]
Nu_ave = h_ave*Ltot/kf_ave;  % average nusselt number [-]

% Theoretical
%-------------
% LOCAL NUSSELT NUMBER--------------------------------
[rho,miu,~,Cp] = AirProperties(Tf_ave,Pbaro);
nu = miu/rho;
alpha = kf_ave/(Cp*rho);

Pr = nu/alpha;
Re = Uinf*x/nu;

Re_c = 5e5;
Re_lam = Re < Re_c; % laminar check (looks like all is laminar)

Nu_th_top = 0.453 * (Re.^(1/2)) * (Pr^(1/3));
Nu_th_bot = 1 - (lead./x).^(9/10);
Nu_th = Nu_th_top./(Nu_th_bot.^(1/9));

h_th = (kf_ave./x).*Nu_th;

% AVERAGE NUSSELT NUMBER-----------------------------
L = lead+Lh; % length of heated plate + leading edge [m]
ReL = Uinf*L/nu;
Re_lam = ReL < Re_c % laminar check (looks like all laminar also)

h_ave_th = 2*(kf_ave/(L-lead))*(0.453*ReL^(1/2)*Pr^(1/3))*(1-(lead/L)^(3/4))^(2/3);
Nu_ave_th = h_ave_th*L/kf_ave;

% PREDICTED SURFACE TEMPERATURE DISTRIBUTION---------
Ts_th = Tinf + qs_flux./h_th;

% PREDICTED HEAT TRANSFER----------------------------
qs_flux_th = h_th.*(Ts_th - Tinf);
qs_th = (w*Lh/Lt)*trapz(x,qs_flux_th);

% Radiation
%-------------
% ESTIMATED HEAT TRANSFER DUE TO RADIATION-----------
sigma = 5.6703e-8; % Stefan-Boltzmann constant [W/m^2-K^4]
e = 0.7; % emissivity of plate [-]

qrad_flux = e*sigma.*(Ts.^4 - Tinf^4);    % local radiation heat flux [W/m^2]
qrad_ave_flux = (1/Lt)*trapz(x,qrad_flux); % average radiation heat flux [W/m^2]
qrad = Lh*w*qrad_ave_flux;                % average radiation heat transfer [W]

q_flux_th = qs_flux_th - qrad_flux;
Ts_th_corrected = Tinf + q_flux_th/h_th;

%======================================
%             FIGURES
%======================================

% FIGURE 1a-------------------------------------------
x_ = (x' - lead)/Lh;
x_top = [x_(1:5); x_(7:11); x_(13:16)];
x_bot = [x_(6); x_(12)];

Ts_top = [Ts(1:5); Ts(7:11); Ts(13:16)];
Ts_bot = [Ts(6); Ts(12)];

plot(x_top,Ts_top,'bo');
hold on
plot(x_bot,Ts_bot,'rs');
plot(x_,Ts_th,'k-');
plot(x_,Ts_th_corrected,'k--');
title('Figure 1a');
ylabel('Surface Temperature [K]');
xlabel('Nondimensional Length [-]');
legend('Measured Top','Measured Bottom','Theoretical','Theoretical + Radiation','Location','Northwest');
hold off

% FIGURE 1b-------------------------------------------
figure;
plot(x_,h,'bo');
hold on
plot(x_,h_th,'b-');
title('Figure 1b');
ylabel('Local Heat Transfer Coefficient [W/m^2-K]');
xlabel('Nondimensional Length [-]');
legend('Experimental','Theoretical');
hold off

% FIGURE 1c-------------------------------------------
figure;
plot(x_,Nu,'bo');
hold on
plot(x_,Nu_th,'b-');
title('Figure 1c');
ylabel('Local Nusselt Number');
xlabel('Nondimensional Length [-]');
legend('Experimental','Theoretical');

%======================================
%             QUESTIONS
%======================================

% QUESTION 2a-----------------------------------------
Nu_percent_diff = ((Nu - Nu_th)./Nu_th).*100;
Ts_percent_diff = ((Ts' - Ts_th)./Ts_th).*100;
h_percent_diff = ((h - h_th)./h_th).*100;

figure;
plot(x_,Nu_percent_diff,'rs');
hold on
plot(x_,Ts_percent_diff,'bo');
plot(x_,h_percent_diff,'k*');
title('Question 2a');
ylabel('Percent Difference [%]');
xlabel('Nondimensional Length');
legend('Nu','Ts','h');

% QUESTION 2b-----------------------------------------
h_ave_percent_diff = ((h_ave - h_ave_th)/h_ave_th)*100
Nu_ave_percent_diff = ((Nu_ave - Nu_ave_th)/Nu_ave_th)*100

% QUESTION 2c-----------------------------------------
heat_lost = (qrad_ave_flux/qs_flux)*100

% QUESTION 2d-----------------------------------------
