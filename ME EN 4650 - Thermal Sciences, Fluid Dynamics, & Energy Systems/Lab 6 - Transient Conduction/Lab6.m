% -------------------------------------
% Alia Zaki
% u1072443
% ME EN 4650
% Lab 6: Transient Conduction
% -------------------------------------

clc
clear all
close all

%======================================
%              DATA
%======================================
data_al = load('4650TransientConductionAluminum.dat');
data_br = load('4650TransientConductionBrass.dat');

r0 = 1*2.54/100;   % NO IMPERIAL SCUM UNITS IN THIS CODE (m)
As = 4*pi*r0^2;    % surface area of sphere (m^2)
V = (4/3)*pi*r0^3; % volume of sphere (m^3)

% aluminum
rho_al = 2760;  % density (kg/m^3)
k_al = 121.4;   % thermal conductivity (W/m-k)
c_al = 895.8;   % specific heat capacity (J/kg-K)

Tbath_al = data_al(26:end,1); % bath temperature (deg C)
Tball_al = data_al(26:end,2); % aluminum ball temperature (deg C) 
t_al = data_al(26:end,3);     % time (s)   

% brass
rho_br = 8500;  % density (kg/m^3)
k_br = 116.0;   % thermal conductivity (W/m-K)
c_br = 382.6;   % specific heat capacity (J/kg-K)

Tbath_br = data_br(30:end,1); % bath temperature (deg C)
Tball_br = data_br(30:end,2); % aluminum ball temperature (deg C) 
t_br = data_br(30:end,3);     % time (s)   

%======================================
%             FIGURES
%======================================
%-------------------------------------- FIGURE 1a
plot(t_al,Tball_al,'bo');
hold on
plot(t_al,Tbath_al,'b--');
hold on
plot(t_br,Tball_br,'ro');
hold on
plot(t_br,Tbath_br,'r--');
legend('T aluminum','Tbath aluminum','T brass','Tbath brass','Location','Southeast');
title('Figure 1a');
ylabel('Temperature (deg C)');
xlabel('Time (s)');
hold off

%-------------------------------------- FIGURE 1b
for i = 1:length(t_al)
    thetaStar_al(i) = (Tball_al(i) - mean(Tbath_al))/(Tball_al(1) - mean(Tbath_al));
end

alpha_al = k_al/(rho_al*c_al);
Fo_al = alpha_al*t_al/(r0^2);

ind_al = find(thetaStar_al>=0.05); % indices of data to use in curve fit
tData_al = t_al(ind_al)';           % time data to use in curve fit
thetaData_al = thetaStar_al(ind_al);  % theta data to use incurve fit
FoData_al = Fo_al(ind_al);

y_al = log(thetaData_al);       % transformed variable
p = polyfit(tData_al,y_al,1);   % straight line fit
m_al = p(1);                    % slope of fit
b_al = p(2);                    % intercept of fit
taoLC_al = -1/m_al;                % thermal time constant
N_al = length(y_al);            % number of data points in curve fit

thetaFit_al = exp(m_al*t_al + b_al);

% standard deviation taoLC
err_data = y_al-(m_al*tData_al+b_al);             % fit errors
s2 = sum(err_data.^2)/(N_al-2);             % mean squared error
sum_ti2 = sum(tData_al.^2);                 % sum of squared time data
sum2_ti = sum(tData_al)^2;                  % squared sum of time data
err_m = sqrt(N_al*s2/(N_al*sum_ti2-sum2_ti));  % error in m
err_taoLC_al = 1/m_al^2*err_m;                 % error in taoLC
err_taoLC_al = 2*err_taoLC_al;

figure;
semilogy(Fo_al',thetaStar_al,'ro');
title('Figure 1b');
ylabel('Theta*');
xlabel('Fourier Number');
hold on
plot(Fo_al,thetaFit_al,'k');
legend('Theta*','Linear Fit');

%-------------------------------------- FIGURE 1c
for i = 1:length(t_br)
    thetaStar_br(i) = (Tball_br(i) - mean(Tbath_br))/(Tball_br(1) - mean(Tbath_br));
end

alpha_br = k_br/(rho_br*c_br);
Fo_br = alpha_br*t_br/(r0^2);

ind_br = find(thetaStar_br>=0.05); % indices of data to use in curve fit
tData_br = t_br(ind_br)';           % time data to use in curve fit
thetaData_br = thetaStar_br(ind_br);  % theta data to use incurve fit
FoData_br = Fo_br(ind_br);

y_br = log(thetaData_br);       % transformed variable
p = polyfit(tData_br,y_br,1);   % straight line fit
m_br = p(1);                        % slope of fit
b_br = p(2);                 % intercept of fit
taoLC_br = -1/m_br;          % thermal time constant
N_br = length(y_br);            % number of data points in curve fit

thetaFit_br = exp(m_br*t_br + b_br);

% standard deviation taoLC
err_data = y_br - (m_br*tData_br + b_br);          % fit errors
s2 = sum(err_data.^2)/(N_br-2);             % mean squared error
sum_ti2 = sum(tData_br.^2);                 % sum of squared time data
sum2_ti = sum(tData_br)^2;                  % squared sum of time data
err_m = sqrt(N_br*s2/(N_br*sum_ti2-sum2_ti));  % error in m
err_taoLC_br = 1/m_br^2*err_m;              % error in taoLC
err_taoLC_br = 2*err_taoLC_br;

figure;
semilogy(Fo_br',thetaStar_br,'ro');
title('Figure 1c');
ylabel('Theta*');
xlabel('Fourier Number');
hold on
plot(Fo_br,thetaFit_br,'k');
legend('Theta*','Linear Fit');

%-------------------------------------- FIGURE 1d
%-------------------------
% LUMPED CAPACITANCE MODEL
%-------------------------
% heat transfer coefficients
hLC_al = rho_al*V*c_al/(taoLC_al*As);   
hLC_br = rho_br*V*c_br/(taoLC_br*As);

% hLC standard deviation (95% confidence interval)
err_hLC_al = 2*(rho_al*V*c_al/(As*(taoLC_al^2)))*err_taoLC_al;
err_hLC_br = 2*(rho_br*V*c_br/(As*(taoLC_br^2)))*err_taoLC_al;

% validate LC model
Bi_al = hLC_al*r0/k_al
Bi_br = hLC_br*r0/k_br

%-------------------------
% ONE-TERM APPROXIMATION
%-------------------------

% ALUMINUM---------------------------------

% Fo_al=(tData_al*alpha_al/r0^2)';            %Fourier number for data set
% ind=Fo_al>=0.5 & thetaStar_al>=0.05;  %indices of data to use in curve fit
% FoData_al=Fo_al(ind);            %Fo data to use in curve fit
% thetaData_al=thetaStar_al(ind);      %theta data to use in curve fit
% N_al=length(ind);              %number of data points in curve fit

zeta1 = linspace(0.01,2,150)';              % parameter values
s2 = zeros(size(zeta1));                    % initialize s2 array
for k = 1:length(zeta1)
    z = zeta1(k);                             % current value of zeta1
    theta_model = (4*(sin(z)-z*cos(z)))/...
        (2*z-sin(2*z))*exp(-z.^2*FoData_al);   % theta based on model
    err = thetaData_al' - theta_model;              % curve fit error
    s2(k) = sum(err.^2)/(N_al-1);                % sum of squared error
end
[s2_min,i_min] = min(s2);                     % find minimum s2 value

% zeta1
zeta1_al = zeta1(i_min);                     % best-fit zeta1 value
err_zeta1_al=Error_Zeta1(zeta1_al,FoData_al,thetaData_al);

% taoES
taoES_al = r0^2/(alpha_al*zeta1_al^2);
err_taoES_al = ((2*r0^2)/(alpha_al*zeta1_al^3))*err_zeta1_al;
err_taoES_al = 2*err_taoES_al;

% hES
hES_al = (k_al/r0)*(1-zeta1_al*cot(zeta1_al));
err_hES_al = (k_br/r0)*(cot(zeta1_al) - zeta1*((cot(zeta1_al))^2+1))*err_zeta1_al;
err_hES_al = mean(err_hES_al);
err_hES_al = 2*err_hES_al;

Q_al = (rho_al*V*c_al)*(1-exp(t_al(end)/taoES_al));

% BRASS-----------------------------------
% Fo_br =tData_br*alpha_br/r0^2;            %Fourier number for data set
% ind=Fo_br>=0.5 & thetaStar_al>=0.05;  %indices of data to use in curve fit
% FoData_br=Fo_br(ind);            %Fo data to use in curve fit
% thetaData_br=thetaStar_br(ind);      %theta data to use in curve fit
% N_br=length(ind);              %number of data points in curve fit

zeta1 = linspace(0.01,2,150)';                % parameter values
s2 = zeros(size(zeta1));                      % initialize s2 array
for k = 1:length(zeta1)
    z = zeta1(k);                             % current value of zeta1
    theta_model = (4*(sin(z)-z*cos(z)))/...
        (2*z-sin(2*z))*exp(-z.^2*FoData_br);   % theta based on model
    err = thetaData_br' - theta_model;              % curve fit error
    s2(k) = sum(err.^2)/(N_br-1);                % sum of squared error
end
[s2_min,i_min] = min(s2);                     % find minimum s2 value

zeta1_br = zeta1(i_min);                     % best-fit zeta1 value
err_zeta1_br = Error_Zeta1(zeta1_br,FoData_br,thetaData_br);

taoES_br = r0^2/(alpha_br*zeta1_br^2);
err_taoES_br = ((2*r0^2)/(alpha_br*zeta1_br^3))*err_zeta1_br;
err_taoES_br = 2*err_taoES_br;

hES_br = (k_br/r0)*(1-zeta1_br*cot(zeta1_br));
err_hES_br = (k_br/r0)*(cot(zeta1_br) - zeta1*((cot(zeta1_br))^2+1))*err_zeta1_br;
err_hES_br = mean(err_hES_br);
err_hES_br = 2*err_hES_br;

Q_br = (rho_br*V*c_br)*(1-exp(t_br(end)/taoES_br));

%======================================
%             QUESTIONS
%======================================
%-------------------------------------- QUESTION 2a
percent_diff = (abs(hES_br - hES_al)/(0.5*(hES_br + hES_al)))*100

%-------------------------------------- QUESTION 2b
yi_al = log(thetaFit_al);
St = sum((yi_al-y_al).^2);
Sr = sum((yi_al - m_al*tData_al - b_al).^2);
R2_al = (St - Sr)/St

yi_br = log(thetaFit_br);
St = sum((yi_br-y_br).^2);
Sr = sum((yi_br - m_br*tData_br - b_br).^2);
R2_br = (St - Sr)/St

%-------------------------------------- QUESTION 2c
r0_al = Bi_al*k_al/hES_al
r0_br = Bi_br*k_br/hES_br

%-------------------------------------- QUESTION 2d
Q_dot_al = Q_al/t_al(end)
Q_dot_br = Q_br/t_br(end)