%-----------------------------
% ME EN 4650
% Lab 1: Cooling Tower
% Alia Zaki
% u1072443
%-----------------------------

close all
clear all
clc

data = xlsread('Lab 1 Cooling Tower Data Sheet');

% inlet mass flow rate
m_dot1 = 16;              % g/s
m_dot2 = 32;              % g/s
m_dot3 = 43;              % g/s
mdot_in = [m_dot1 m_dot2 m_dot3];

T_amb = data(1,7);              % C, ambient temperature
P_atm = (data(2,7)/7.501)*1e3;  % Pa, atmospheric pressure
qdot_in = data(3,7)*1000;       % W, heat transfer in
D_makeup = data(4,7)/100;       % m, makeup tube diameter

% point of temperature measurement heights (m)
A = 0/100;      
F = 24.8/100;
G = 48.3/100;
H = 71.8/100;
B = 100/100;

% Plot 1a
% x -> positions
% y -> wet bulb temperature (Twb), water temperature (Tw), range, and
% approach

Tw2 = data(22:26,8)';   % water temperature (C)
Twb2 = data(8:12,9)';   % wet bulb temperature (C)
h = [A B H G F];        % height of cooling tower (m)
plot(h,Tw2,'rs')
hold on
plot(h,Twb2,'bo')
legend('Water','Wet Bulb')
xlabel('Positions (m)')
ylabel('Temperature (C)')
title('Plot 1a')

% Range/Approach
Tw_out = Tw2(2);
Tw_in = Tw2(1);
Twb_in = Twb2(1);

yline(Tw_out,'k','LineWidth',1);
yline(Tw_in,'k','LineWidth',1);
yline(Twb_in,'k','LineWidth',1);

% Plot 1b
Twb_in = data(3,1:3)';  % inlet wet bulb temperature (C)
Tw_in = data(6,1:3)';   % inlet water temperature (C)
Tw_out = data(7,1:3)';  % oultet water temperature

for i = 1:3
    R = Tw_out(i) - Tw_in(i);   % range (C)
    A = Tw_out(i) - Twb_in(i);  % approach (C)
    
    eff(i) = R/(R+A);           % efficiency
end

figure
plot(mdot_in,eff,'bo-')
xlabel('Mass Flow Rate (g/s)')
ylabel('Efficiency (%)')
title('Plot 1b')

% Plot 1c
A = 0/100;      
F = 24.8/100;
G = 48.3/100;
H = 71.8/100;
B = 100/100;

z = [A B H G F];    % cooling tower heights (m)

Tdb = data(8:12,7:9);   % dry bulb temperatures (C)
Twb = data(15:19,7:9);  % wet bulb temperatures (C)

for i = 1:5
    for j = 1:3
        [Tdb_,w_,phi,h,Tdp,v,Twb_]=Psychrometrics('tdb',Tdb(i,j),'twb',Twb(i,j),'p',P_atm/1000);
        w(i,j) = w_;    % specific humidity (kg/kg)
    end
end

figure;
plot(z,w(1:5,1),'gd');
hold on
plot(z,w(1:5,2),'bo');
plot(z,w(1:5,3),'rs');
hold off
legend('m-dot = 16 g/s', 'm-dot = 32 g/s', 'm-dot = 43 g/s','Location','southeast')
xlabel('Positions (m)');
ylabel('Specific Humidity (kg/kg)');
title('Plot 1c')

% Plot 1d
A = 0/100;      
F = 24.8/100;
G = 48.3/100;
H = 71.8/100;
B = 100/100;

z = [A B H G F];    % heights of cooling tower (m)

figure
plot(z,Tdb(1:5,1)','gd')
hold on
plot(z,Tdb(1:5,2)','bo')
plot(z,Tdb(1:5,3)','rs')
hold off
legend('m-dot = 16 g/s', 'm-dot = 32 g/s', 'm-dot = 43 g/s','Location','southwest')
xlabel('Positions (m)');
ylabel('Dry Bulb Temperature (C)');
title('Plot 1d')

% Plot 1e
Tdb = data(8:9,7:9);    % dry bulb temperatures (C)
Twb = data(15:16,7:9);  % wet bulb temperatures (C)
w_in = w(1,:);          % inlet specific humidity (kg/kg)
w_out = w(2,:);         % outlet specific humidity (kg/kg)
deltaP = 10;            % mm H2O, manometer reading 

for i = 1:2
    for j = 1:3
        [Tdb_,w,phi,h,Tdp,v_,Twb_]=Psychrometrics('tdb',Tdb(i,j),'twb',Twb(i,j),'p',P_atm/1000);
        v(i,j) = v_; % spefici volume of dry air m3/kg
                 
    end
end

for i = 1:3
    % ma_dot_out = ma_dot_in due to mass conservation 
    % mass flow rate of dry air (g/s)    
    ma_dot_in(i) = 0.0137*sqrt(deltaP/(1+w_in(i))*v(1,i));
    ma_dot_out(i) = 0.0137*sqrt(deltaP/(1+w_out(i))*v(2,i));
    % mass flow rate of vapor (g/s)
    mv_dot_in = w_in(i)*ma_dot_in; 
    mv_dot_out = w_out(i)*ma_dot_out;
    % mass flow rate of water out of cooling tower 
    mw_dot_out = mdot_in(i) + mv_dot_in - mv_dot_out;
    ratio = mw_dot_out/mdot_in(i);
end

Tw_in = data(6,1:3); % inlet water temperatures

figure;
plot(Tw_in,ratio,'bo-');
xlabel('Inlet Temperature (C)');
ylabel('Ratio of Outlet to Inlet Mass Flow Rate');
title('Plot 1e');

% Plot 1f
Tdb = data(8:9,7:9);
Twb = data(15:16,7:9);

Tdb_in = Tdb(1,:);
Twb_in = Twb(1,:);
Tdb_out = Tdb(2,:);
Twb_out = Twb(2,:);

for i = 1:3
    [Tdb_,w,phi,ha_in,Tdp,v,Twb_]=Psychrometrics('tdb',Tdb_in(i),'twb',Twb_in(i),'p',P_atm/1000);
    [Tdb_,w,phi,ha_out,Tdp,v_,Twb_]=Psychrometrics('tdb',Tdb_out(i),'twb',Twb_out(i),'p',P_atm/1000);
    
    qa_dot(i) = ma_dot_in(i)*(ha_out-ha_in);  % heat transfer rate into dry air
    qamb_dot(i) = qdot_in - qa_dot(i);        % heat transfer rate into ambient 
end

figure
plot(Tw_in,qa_dot/1000,'bo-')
hold on
plot(Tw_in,qamb_dot/1000,'ro-')
hold off
title('Plot 1f')
ylabel('Heat Transfer Rate (kW)')
xlabel('Inlet Temperature (C)')
legend('Air','Ambient');


    


