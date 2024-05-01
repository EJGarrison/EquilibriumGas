%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   MAE 250F Final Project    %
%     Gas Dynamics Code       %
%                             %
%  Authored By Evan Garrison  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clearvars; close all;
%warning('off')

%% Part 1: Equilibrium Gas

% Given Knowns
T0 = 300:100:15000;                             % [K]
p0 = 101325*[0.0001 0.001 0.01 0.1 1 10 100];   % [Pa]

% Initialize
Z = zeros(length(p0),length(T0));
EZ = zeros(length(p0),length(T0));
SZ = zeros(length(p0),length(T0));
rho = zeros(length(p0),length(T0));
h = zeros(length(p0),length(T0));
e = zeros(length(p0),length(T0));
s = zeros(length(p0),length(T0));
R = zeros(length(p0),length(T0));

for i = 1:length(p0)
    for j = 1:length(T0)

        result = equilibriumgas(p0(i),T0(j));
        Z(i,j) = result(1);
        EZ(i,j) = result(2);
        SZ(i,j) = result(3);
        rho(i,j) = result(4);
        h(i,j) = result(5);
        e(i,j) = result(6);
        s(i,j) = result(7);
        R(i,j) = result(8);

    end
end

save result1.mat Z EZ SZ rho h e s R p0 T0

%% Plotting vs Hansen (NACA TN 4150)

% Compressibility
z_hansen = readtable("Digitized Data\compressibility_hansen.csv","VariableNamingRule","preserve");% figure 1 of Hansen
f1 = figure(1);
hold on
colors = [0 0 153; 28 134 238; 153 204 255; 255 204 153; 255 153 51; 255 128 0; 204 102 0]/255;
for k1 = 1:length(p0)
    plot(z_hansen.(2*k1-1), z_hansen.(2*k1), Color=colors(k1,:))
    hold on
end
for k2 = 1:length(p0)
    scatter(T0, Z(k2,:), [], colors(k2,:), "+")
    hold on
end
ylim([0 4.1])
title("Compressibility vs Temperature")
xlabel("Temperature [K]")
ylabel("Compressibility (Z)")
legend("0.0001 atm (Hansen)","0.001 atm (Hansen)","0.01 atm (Hansen)", ...
    "0.1 atm (Hansen)","1 atm (Hansen)","10 atm (Hansen)","100 atm (Hansen)", ...
    "0.0001 atm","0.001 atm","0.01 atm", ...
    "0.1 atm","1 atm","10 atm","100 atm",'Location','northwest','NumColumns',2)
hold off

% Dimensionless Energy
e_hansen = readtable("Digitized Data\energy_hansen.csv","VariableNamingRule","preserve");% figure 2 of Hansen (ZE/RT)
f2 = figure(2);
hold on
for l1 = 1:length(p0)
    plot(e_hansen.(2*l1-1), e_hansen.(2*l1), Color=colors(l1,:))
    hold on
end
for l2 = 1:length(p0)
    scatter(T0, EZ(l2,:), [], colors(l2,:), "+")
    hold on
end
title("Dimensionless Energy vs Temperature")
ylim([0 50])
xlabel("Temperature [K]")
ylabel("Dimensionless Energy EZ/RT")
legend("0.0001 atm (Hansen)","0.001 atm (Hansen)","0.01 atm (Hansen)", ...
   "0.1 atm (Hansen)","1 atm (Hansen)","10 atm (Hansen)","100 atm (Hansen)", ...
   "0.0001 atm","0.001 atm","0.01 atm", ...
   "0.1 atm","1 atm","10 atm","100 atm",'Location','northwest','NumColumns',2)
hold off

% Dimensionless Entropy
s_hansen = readtable("Digitized Data\entropy_hansen.csv","VariableNamingRule","preserve");% Figure 3 of Hansen(ZS/R)
f3 = figure(3);
hold on
for m1 = 1:length(p0)
    plot(s_hansen.(2*m1-1), s_hansen.(2*m1), Color=colors(m1,:))
    hold on
end
for m2 = 1:length(p0)
    scatter(T0, SZ(m2,:), [], colors(m2,:), "+")
    hold on
end
title("Dimensionless Entropy vs Temperature")
ylim([0 140])
xlabel("Temperature [K]")
ylabel("Dimensionless Entropy SZ/R")
legend("0.0001 atm (Hansen)","0.001 atm (Hansen)","0.01 atm (Hansen)", ...
   "0.1 atm (Hansen)","1 atm (Hansen)","10 atm (Hansen)","100 atm (Hansen)", ...
   "0.0001 atm","0.001 atm","0.01 atm", ...
   "0.1 atm","1 atm","10 atm","100 atm",'Location','northwest','NumColumns',2)
hold off

% figure(9)
% for mm = 1:length(p0)
%     scatter(T0, h(mm,:), [], colors(mm,:), "+")
%     hold on
% end

% figure(4)
%hold on
%a_hansen = readtable("Digitized Data\speedofsound_hansen_datasets.csv");% Figure 5 (a^2rho/p)

%% Part 2: Topic 1 - Normal Shock
% Define Universal Constants %

k = 1.380649*10^-23;        % Boltzmann Constant [J/K]
h = 6.62607015*10^-34;      % Planck's Consant [J-s]
NA = 6.02214076*10^23;      % Avogadro's Number [mol^-1]
Rhat = NA*k;                % Universal Gas Constant [J/mol-K]
M0 = 0.02896968;            % Molar Mass of Dry Air [kg/mol]

% Define Free Stream Conditions %
% ft to m from https://www.nist.gov/pml/us-surveyfoot/revised-unit-conversion-factors

Hft = [35900 59800 82200 100000 120300 154800 173500 ...
    200100 230400 259700 294800 322900];
Hm = Hft*0.3048;
Tatm = [216.65 216.65 223.65 230.89 249.37 270.65 255.268 ...
    230.05 205.65 188.983 189.83 212.89];                               % K
patm = [2.0018e2 5.8313e1 1.8474e1 7.7067 3.0121 6.6938e-1 3.2782e-1 ...
     9.2140e-2 1.9034e-2 4.4568e-3 5.8439e-4 1.2454e-4]*10^2; % Pa
rhoatm = [0.32190 9.3767e-2 2.8777e-2 1.1628e-2 4.2079e-3 8.6160e-4 ...
    4.4738e-4 1.3953e-4 3.2245e-5 8.2196e-6 1.061e-6 1.954e-7];        % kg/m^3

x1 = struct(u = (6000:10000:46000)*0.3048, T = Tatm, p = patm,...
    gamma = 1.4, Hm = Hm);

for ii = 1:length(Hm)
    result2 = equilibriumgas(x1.p(ii), x1.T(ii));
    x1.h(ii) = result2(5);
    x1.rho(ii) = result2(4);
    x1.R(ii) = result2(8);
end

h2 = zeros(length(Hm), length(x1.u));
p2 = zeros(length(Hm), length(x1.u));
T2 = zeros(length(Hm), length(x1.u));
rho2 = zeros(length(Hm), length(x1.u));
r1r2 = 0.1*ones(length(Hm), length(x1.u));
tol = 1e-15;

for i = 1:length(Hm)
    for j = 1:length(x1.u)
        
        r1r2last = 0;
        M1 = x1.u(j)/sqrt((x1.gamma*x1.R(i))*x1.T(i));
        n = 1;

        while abs(r1r2(i,j) - r1r2last) > tol && n < 50

            r1r2last = r1r2(i,j);
            p2(i,j) = x1.p(i) + x1.rho(i)*x1.u(j)^2*(1 - r1r2(i,j));
            h2(i,j) = x1.h(i) + 0.5*x1.u(j)^2*(1 - r1r2(i,j)^2);

            T2high = ((2*x1.gamma*M1^2 - (x1.gamma))*...
                ((x1.gamma - 1)*M1^2 + 2))...
                /((x1.gamma+1)*M1^2)*x1.T(i); %90*x1.T(i);
            T2low = 1.75*x1.T(i);
            
            h2new = 0;
            result2high = equilibriumgas(p2(i,j), T2high);
            h2high = result2high(5);
            result2low = equilibriumgas(p2(i,j), T2low);
            h2low = result2low(5);
            m = 1;

            while abs(h2(i,j) - h2new) > 10 && m < 30

                T2new = ((T2high - T2low)/(h2high - h2low))*(h2(i,j) - h2low) + T2low;
                result2 = equilibriumgas(p2(i,j), T2new);
                h2new = result2(5);

                if h2new < h2(i,j)
                   T2low = T2new;
                   h2low = h2new;
                else
                   T2high = T2new;
                   h2high = T2new;
                end

                if T2low < x1.T(i)
                    T2low = x1.T(i);
                end
  
                m = m + 1;
            end

            rho2(i,j) = result2(4);
            r1r2(i,j) = x1.rho(i)/rho2(i,j);
            n = n + 1;

        end

        T2(i,j) = T2new;
        r1r2(i,j) = rho2(i,j)/x1.rho(i);

    end
end

save result2.mat r1r2 T2

%% Plotting vs Anderson
% plot and compare to Fig 14.4, 14.5 of Anderson
% u1: 3000 to 46,000 ft/s (x-axis); T2: 500 to 14000 K and rho2/rho1: 4 to 18 (y-axis)
% p1: altitude from 35,900 ft to 322,900 ft (each line)

T2_anderson = readtable("Digitized Data\T2_anderson.csv","VariableNamingRule","preserve");
r2r1_anderson = readtable("Digitized Data\rho2rho1_anderson2.csv","VariableNamingRule","preserve");
f4 = figure(4);
hold on
colors = [linspace(100,255,length(Hm)); linspace(0,102,length(Hm)); zeros(1,length(Hm))]'/255;
for k3 = 1:length(Hm)
    plot(sort(T2_anderson.(2*k3-1))*0.3048, sort(T2_anderson.(2*k3)), Color=colors(k3,:))
    hold on
end
for k4 = 1:length(Hm)
    scatter(x1.u, T2(k4,:), [], colors(k4,:), "+")
    hold on
end
ylim([0 14000])
title("Temperature Behind Shock vs Freestream Velocity")
xlabel("$u_1$ [m/s]",Interpreter='latex',FontSize=16)
ylabel("$T_2$ [K]",Interpreter='latex',FontSize=16)
legend("35900 ft (Anderson)","59800 ft (Anderson)","82200 ft (Anderson)", ...
    "100000 ft (Anderson)","120300 ft (Anderson)","154800 ft (Anderson)","173500 ft (Anderson)", ...
    "200100 ft (Anderson)","230400 ft (Anderson)","259700 ft (Anderson)","294800 ft (Anderson)",...
    "322900 ft (Anderson)",'Location','northwest','NumColumns',2)
hold off

f5 = figure(5);
hold on
colors = [linspace(100,255,length(Hm)); linspace(0,102,length(Hm)); zeros(1,length(Hm))]'/255;
for l3 = 1:length(Hm)
    plot(sort(r2r1_anderson.(2*l3-1))*0.3048, sort(r2r1_anderson.(2*l3)), Color=colors(l3,:))
    hold on
end
for l4 = 1:length(Hm)
    scatter(x1.u, 1./r1r2(l4,:), [], colors(l4,:), "+")
    hold on
end
ylim([0 24])
title("Ratio of Densities vs Freestream Velocity")
xlabel("$u_1$ [ft/s]",Interpreter='latex',FontSize=16)
ylabel("$\frac{\rho_2}{\rho_1}$",Interpreter='latex',FontSize=20)
legend("35900 ft (Anderson)","59800 ft (Anderson)","82200 ft (Anderson)", ...
    "100000 ft (Anderson)","120300 ft (Anderson)","154800 ft (Anderson)","173500 ft (Anderson)", ...
    "200100 ft (Anderson)","230400 ft (Anderson)","259700 ft (Anderson)","294800 ft (Anderson)",...
    "322900 ft (Anderson)",'Location','northwest','NumColumns',2)
hold off