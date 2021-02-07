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

% DATA----------------------------------------
data = xlsread('BlackbodyRadiation_SampleData.xlsx');

Tatm = data(1,2);            % [deg C] atmospheric temperature
Patm = data(2,2)*mmHgToPa;   % [Pa] atmospheric pressure

% variable h -> separation distance
Th = data(6:10,1)' + 273.15;    % [K] blackbody temperature
Dh = data(6:10,2)'*inTom;       % [m] aperture diameter of source 
hh = data(6:10,3)'*inTom;       % [m] separation distance 
qh = data(6:10,4)';             % [uW] heat transfer rate

% variable T -> blackbody temperature
TT = data(6:14,6)' + 273.15;    % [K] black body temperature
DT = data(6:14,7)'*inTom;        % [m] aperture diameter of source 
hT = data(6:14,8)'*inTom;        % [m] aperture radius of source
qT = data(6:14,9)';              % [uW] heat transfer rate

% variable D -> aperture diameter of source
TD = data(23:26,1)'+ 273.15;     % [deg C] black body temperature
DD = data(23:26,2)'*inTom;       % [m] aperture diameter of source 
hD = data(23:26,3)'*inTom;       % [m] aperture radius of source
qD = data(23:26,4)';             % [uW] heat transfer rate


% DIMENSIONS----------------------------------
Ad = 1*(mmTom)^2; % [m^3] surface area of sensing chip

As_h = (pi/4)*Dh.^2;
As_T = (pi/4)*DT.^2;
As_D = (pi/4)*DD.^2;

%======================================
%               FIGURES
%======================================

e = 1;           % [-] emissivity of surface (blackbody, e = 1)
sigma = 5.67e-8; % [W/m^2-K^4] Stefan-Boltzmann 

% [W/m^2] total emissive power 
Eb_h = sigma*Th;
Eb_T = sigma*TT;
Eb_D = sigma*TD;

% [W/m^2] total intensity of emitted radiation
Ib_h = Eb_h/pi;
Ib_T = Eb_T/pi;
Ib_D = Eb_D/pi;

% [uW] theoretical heat transfer rate
qh_th = 1e6*(sigma.*(Th.^4).*Ad.*(Dh.^2)./(Dh.^2 + 4*hh.^2));
qT_th = 1e6*(sigma.*(TT.^4).*Ad.*(DT.^2)./(DT.^2 + 4*hT.^2));
qD_th = 1e6*(sigma.*(TD.^4).*Ad.*(DD.^2)./(DD.^2 + 4*hD.^2));

% FIGURE 1a-------------------------------------
x = 1:30;
y = x;

plot(qh,qh_th,'bo');
hold on
plot(qT,qT_th,'rs');
plot(qD,qD_th,'gd');
plot(x,y,'k');
hold off
legend('Exp #1','Exp #2','Exp #3','Location','Southeast');
xlabel('Experimental Heat Transfer Rate [uW]');
ylabel('Theoretical Heat Transfer Rate [uW]');
title('Figure 1a');

% FIGURE 1b-------------------------------------
inv_hh = 1./(hh.^2);
p1 = polyfit(inv_hh,qh,1);
y = polyval(p1,inv_hh);

% calculated R^2
qh_mean = mean(qh);

SR = sum((qh - y).^2);
ST = sum((qh - qh_mean).^2);

R2 = 1 - SR/ST;

woop = ['Linear Fit (R^2 = ', num2str(R2),')'];

figure;
plot(inv_hh,qh,'bo');
hold on
plot(inv_hh,y,'k');
hold off
legend('Data',woop,'Location','Southeast');
ylabel('Experimental Heat Transfer Rate [uW]');
xlabel('Reciprocal of the Squared Separation Distance [m^{-2}]');
title('Figure 1b');

% FIGURE 1c-------------------------------------
T4 = TT.^4;

p2 = polyfit(T4,qT,1);
y2 = polyval(p2,T4);

% calculated R^2
qT_mean = mean(qT);

SR = sum((qT - y2).^2);
ST = sum((qT - qT_mean).^2);

R2 = 1 - SR/ST;
woop = ['Linear Fit (R^2 = ', num2str(R2,4),')'];

figure;
plot(T4,qT,'rs');
hold on
plot(T4,y2,'k');
hold off
legend('Data',woop,'Location','Southeast');
ylabel('Experimental Heat Transfer Rate [uW]');
xlabel('Measured Blackbody Temperatue to the 4th Power [K^4]');
title('Figure 1c');

% FIGURE 1d-------------------------------------
D2 = DD.^2;

p3 = polyfit(D2,qD,1);
y3 = polyval(p3,D2);

% calculated R^2
qD_mean = mean(qD);

SR = sum((qD - y3).^2);
ST = sum((qD - qD_mean).^2);

R2 = 1 - SR/ST;
woop = ['Linear Fit (R^2 = ', num2str(R2),')'];

figure;
plot(D2,qD,'gd');
hold on
plot(D2,y3,'k');
hold off
legend('Data',woop,'Location','Southeast');
ylabel('Experimental Heat Transfer Rate [uW]');
xlabel('Squared Diameter of Aperture [m^2]');
title('Figure 1d');

%======================================
%              QUESTIONS
%======================================

perc_diff_h = (1/length(qh))*sum((abs(qh - qh_th)/qh_th)*100)
perc_diff_T = (1/length(qT))*sum((abs(qT - qT_th)/qT_th)*100)
perc_diff_D = (1/length(qD))*sum((abs(qD - qD_th)/qD_th)*100)





