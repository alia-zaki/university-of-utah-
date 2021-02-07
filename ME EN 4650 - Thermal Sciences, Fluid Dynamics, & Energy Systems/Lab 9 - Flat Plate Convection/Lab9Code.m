%%%%%%%%%%%%%%%%%%%%%
% Jacob Gonzalez-Perez
% u0913317
% ME EN 4610
% Lab 9: Flat Plate Convection
% April 3, 2020
%%%%%%%%%%%%%%%%%%%%%

close all;
clear;
clc;


%% ========================================%%
%                                           %
%               DATA ANALYSIS               %
%                                           %
% ========================================= %

% Read Data__________________________________
data = readmatrix('ConvectionData.xlsx');

num_top = data([1:5 7:11 13:16],1);     % Number of thermocouple
T_top = Convert('t',data([1:5 7:11 13:16],2),'C','K');  % Top thermocouple temperature [C2K]
x_top = Convert('l',data([1:5 7:11 13:16],3),'mm','m'); % Location [mm]

num_bot = data([6 12],1);
T_bot = Convert('t',data([6 12],2),'C','K');    % Bottom thermocouple temperature [C2K]
x_bot = Convert('l',data([6 12],3),'mm','m');

P_baro = Convert('p',data(1,5),'mmHg','Pa'); % Barometric pressure [mmHg]
T_inf = Convert('t',data(2,5),'C','K');  % Ambient temperature [C2K]
P_dyn = Convert('p',data(3,5),'H20','Pa');  % Dynamic pressure [H202Pa]
[rho,~,~,~] = AirProperties(T_inf,P_dyn);
U_inf = sqrt((P_dyn*2)/rho); % Freestream velocity
V = data(4,5);      % Voltage [VAC]
R = data(5,5);      % Resistance [Ohm]




% Given______________________________________
xi = Convert('l',77,'mm','m');    % [mm2m]
L = Convert('l',230,'mm','m');    % Length [mm]
W = Convert('l',144,'mm','m');    % Width [mm]
L_h = Convert('l',153,'mm','m');  % Heated length [mm]
L_T = x_top(end) - x_top(1);
w = Convert('l',68,'mm','m');     % Heated width [mm]
t = Convert('l',13.8,'mm','m');   % Thickness [mm]


% Known_______________________________________
emis = 0.7;
sigma = 5.6703E-8;  % Stefan-Botzmann constant [W/m^2-K]
Re_c = 5E5;         % Critical Reynolds number


% Calculations from Measurement Data_________
% Heat flux from top surface
qFlux_s = (V^2)/(2*R*L_h*w);

% Heat transfer coefficient
h_x = qFlux_s./(T_top - T_inf);
hAvg_L = (1/L_T).*trapz(h_x,x_top);

% Nusselt Number
T_f = (T_top +T_inf)/2;
[~,~,k_f,~]=AirProperties(T_f,P_baro);
Nu_x = h_x.*x_top/k_f;
TAvg = mean(T_top);
[rho,mu,kAvg_f,c_p] = AirProperties(TAvg,P_baro);
NuAvg_L = hAvg_L*L/kAvg_f;


% Theoretical Predictions_____________________
% Local Nusselt Number and Heat Transfer Coefficient
Re_x = rho*U_inf*x_top/mu;
Pr = (mu/rho) / (kAvg_f/(c_p*rho));

Nu_xth = 1:length(x_top);
for i = 1:length(x_top)
    Re_x - Re_c
    if Re_x(i) < Re_c
        numer = 0.453*(Re_x(i)^(1/2))*(Pr^(1/3));
        denom = (1 - ((xi/x_top(i))^(3/4)) )^(1/3);
        Nu_xth(i) = numer/denom;
    else
        numer = 0.031*(Re_x(i)^(4/5))*(Pr^(1/3));
        denom = (1 - ((xi/x_top(i))^(9/10)) )^(1/9);
        Nu_xth(i) = numer/denom;
    end
end

h_xth = (kAvg_f/x_top).*Nu_xth;

%

% Predicted Surface Temperature Distribution 
% T_sth = T_inf + (qFlux_s/h_xth);


% Estimated Heat Transfer due to Radiation____
% T_sRad = qFlux_radL * L_T




%% ========================================%%
%                                           %
%                FIGURE 1A                  %
%                                           %
% ========================================= %

% Addtional Calculations_____________________
xNon_top = (x_top - xi)/L;
xNon_bot = (x_bot - xi)/L;


% Plotting___________________________________
figure
hold on
plot(xNon_top,T_top,'bo')
hold on
plot(xNon_bot,T_bot,'rs')

xlabel('Nondimensional Length [-]')
ylabel('Surface Temperature [K]')
legend('Top Surface','Bottom Surface','Location','east')

xlim([0 1])







