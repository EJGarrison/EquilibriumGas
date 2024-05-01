function flowstate = equilibriumgas(p,T)
    % Debugging and Info
    % equilibriumgas(p,T) uses the results from statistical mechanics to
    % numerically solve for the properties of air at 
    %p = 101325*100;
    %T = 15000;
    % Define Constants %
    % Constants are from Hansen (1958, NACA Technical Note 4359)
    % or are universal and taken from physics.nist.gov
    
    % Given Knowns
    V = 1;                      % [m^3]
    N2per = 0.79;               % Mixture percent of N2 at low temp (air)
    O2per = 0.21;               % Mixture percent of O2 at low temp (air)
    
    % Universal Constants
    k = 1.380649*10^-23;        % Boltzmann Constant [J/K]
    h = 6.62607015*10^-34;      % Planck's Consant [J-s]
    NA = 6.02214076*10^23;      % Avogadro's Number [mol^-1]
    Rhat = NA*k;                % Universal Gas Constant [J/mol-K]
    
    % Electron
    e = struct("m", 9.1093837015*10^-31, "R", Rhat/(5.4857990888*10^-7), ...
        "Mhat", 5.4857990888*10^-7, "g", [2 0]);
    % Electron Mass [kg], Specific gas constant for eletrons,
    % Molecular weight of an electron, Electronic degeneracy of an electron
    
    % Monatomic Nitrogen
    N = struct("thetai", 168800, "Mhat", 14*10^-3, "m", 14*10^-3/NA, "R", Rhat/(14*10^-3), ...
        "g", [4 10 6]);           
    % thetai [K], Mhat [kg/mol], m [kg], R_N [J/kg-K], ...
    % Electric degeneracy (g0, g1, g2) []  
    
    % Diatomic Nitrogen
    N2 = struct("thetad", 113200, "thetav", 3390, "thetar", 2.9, "thetai", 181000, ...
        "Mhat", 2*N.Mhat, "m", 2*N.m, "R", N.R/2, "g", [1 3]);
    % thetad [K], thetav [K] = h*nu/k, thetar [K] = h^2/(8*pi*I*k), thetai [K], ...
    % Mhat [kg/mol], m [kg], R_N2 [J/kg-K], Electric degeneracy []
    
    % Ionized Nitrogen
    Nplus = struct(N);
    Nplus.g = [1 3 5 5 1 5];
    
    % Monatomic Oxygen
    O = struct("thetai", 158000, "Mhat", 15.9994*10^-3, "m", 15.9994*10^-3/NA, ...
        "R", Rhat/(15.9994*10^-3), "g", [5 3 1 5 1]);
    % thetai [K], Mhat [kg/mol], m [kg], R [J/kg-K], Electric degeneracy []
    
    % Diatomic Oxygen
    O2 = struct("thetad", 59500, "thetav", 2270, "thetar", 2.1, "thetai", 142000, ...
        "Mhat", 2*O.Mhat, "m", 2*O.m, "R", O.R/2, "g", [3 2 2]);
    % thetad [K], thetav [K], thetar [K], thetai [K], ...
    % Mhat [kg/mol], m [kg], R [J/kg-K], Electric degeneracy []
    
    % Ionized Oxygen
    Oplus = struct(O);
    Oplus.g = [4 10 6];
    
    % Nitric Oxide
    NO = struct("thetad", 75500, "thetav", 2740, "thetar", 2.5, "thetai", 108000, ...
        "m", N.m+O.m, "R", Rhat/(N.Mhat+O.Mhat), "Mhat", (N.Mhat+O.Mhat), "g", [2 2 8]);
    % [K], thetav [K], thetad [K], thetai [K], Electric degeneracy []
            
    % Calculate Partition Functions %

    % Monatomic Nitrogen
    N.Qtr = V*(2*pi*N.m*k*T/h^2)^1.5;
    N.Qel = N.g(1) + N.g(2)*exp(-27700/T) + N.g(3)*exp(-41500/T);
    N.Q = N.Qtr*N.Qel;
    
    % Diatomic Nitrogen
    N2.Qtr = V*(2*pi*N2.m*k*T/h^2)^1.5;
    N2.Qel = N2.g(1);
    N2.Qr = T/(2*N2.thetar);
    N2.Qv = 1/(1 - exp(-N2.thetav/T));
    N2.Q = N2.Qtr*N2.Qel*N2.Qr*N2.Qv;
    
    % Ionized Nitrogen
    Nplus.Qtr = V*(2*pi*(N.m - e.m)*k*T/h^2)^1.5;
    Nplus.Qel = Nplus.g(1) + Nplus.g(2)*exp(-70.6/T) + ...
        Nplus.g(3)*exp(-188.9/T) + Nplus.g(4)*exp(-22000/T) + ...
        Nplus.g(5)*exp(-47000/T) + Nplus.g(6)*exp(-67900/T);
    Nplus.Q = Nplus.Qel*Nplus.Qtr;
    
    % Nitric Oxide
    NO.Qtr = V*(2*pi*NO.m*k*T/h^2)^1.5;
    NO.Qel = NO.g(1) + NO.g(2)*exp(-174/T);
    NO.Qr = T/(NO.thetar);
    NO.Qv = 1/(1 - exp(-NO.thetav/T));
    NO.Q = NO.Qtr*NO.Qel*NO.Qr*NO.Qv;
   
    % Monatomic Oxygen
    O.Qtr = V*(2*pi*O.m*k*T/h^2)^1.5;
    O.Qel = O.g(1) + O.g(2)*exp(-228/T) + O.g(3)*exp(-326/T) + ...
        O.g(4)*exp(-22800/T) + O.g(5)*exp(-48600/T);
    O.Q = O.Qtr*O.Qel;
    
    % Diatomic Oxygen
    O2.Qtr = V*(2*pi*O2.m*k*T/h^2)^1.5;
    O2.Qel = O2.g(1) + O2.g(2)*exp(-11390/T) + O2.g(3)*exp(-18990/T);
    O2.Qr = T/(2*O2.thetar);
    O2.Qv = 1/(1 - exp(-O2.thetav/T));
    O2.Q = O2.Qtr*O2.Qel*O2.Qr*O2.Qv;
    
    % Ionized Oxygen
    Oplus.Qtr = V*(2*pi*(O.m - e.m)*k*T/h^2)^1.5;
    Oplus.Qel = Oplus.g(1) + Oplus.g(2)*exp(-38600/T) + Oplus.g(3)*exp(-58200/T);
    Oplus.Q = Oplus.Qtr*Oplus.Qel;
    
    % Electron
    e.Qtr = V*(2*pi*e.m*k*T/h^2)^1.5;
    e.Qel = 2;
    e.Q = e.Qtr*e.Qel;
    
    % Reaction Constants
    N2.Kp = (k*T/V)*(N.Q^2 / N2.Q)*exp(-N2.thetad / T);
    O2.Kp = (k*T/V)*(O.Q^2 / O2.Q)*exp(-O2.thetad / T);
    NO.Kp = (k*T/V)*((N.Q*O.Q) / NO.Q)*exp(-NO.thetad / T);
    O.Kp = (k*T/V)*((Oplus.Q*e.Q) / O.Q)*exp(-O.thetai / T);
    N.Kp = (k*T/V)*((Nplus.Q*e.Q) / N.Q)*exp(-N.thetai / T);

    % Partial Pressures %
    if T < 1000
        O2.p = O2per*p;
        O.p = 0;
        NO.p = 0;
        N2.p = N2per*p;
        N.p = 0;
        Oplus.p =  0;
        Nplus.p =  0;
        e.p = 0;
    else
        syms pO2 pO pOplus pNO pN2 pN pNplus pe
        f1 = (pO^2)/pO2 - O2.Kp == 0;
        f2 = (pN^2)/pN2 - N2.Kp == 0;
        f3 = ((2*pO2 + pO + pOplus + pNO)/(pNO + 2*pN2 + pN + pNplus)) - O2per/N2per == 0;
        f4 = pOplus*pe/pO - O.Kp == 0;
        f5 = pNplus*pe/pN - N.Kp == 0;
        f6 = pOplus + pNplus - pe == 0;
        f7 = pO*pN/pNO - NO.Kp == 0;
        f8 = pO2 + pO + pOplus + pNO + pN2 + pN + pNplus + pe - p == 0;
    
        psol = vpasolve([f1 f2 f3 f4 f5 f6 f7 f8], ...
            [pO2 pO pOplus pNO pN2 pN pNplus pe], ...
            [0,p+1; 0,p+1; 0,p+1; 0,p+1; 0,p+1; 0,p+1; 0,p+1; 0,p+1]);
        
            O2.p = double(psol.pO2);
            O.p = double(psol.pO);
            NO.p = double(psol.pNO);
            N2.p = double(psol.pN2);
            N.p = double(psol.pN);
            Oplus.p = double(psol.pOplus);
            Nplus.p = double(psol.pNplus);
            e.p = double(psol.pe);

        if isempty(psol.pO2) == 1
            f3 = ((2*pO2 + pO + pNO)/(pNO + 2*pN2 + pN)) - O2per/N2per == 0;
            f8 = pO2 + pO + pNO + pN2 + pN - p == 0;
            psol = vpasolve([f1 f2 f3 f7 f8], ...
                [pO2 pO pNO pN2 pN], ...
                [0,p+1; 0,p+1; 0,p+1; 0,p+1; 0,p+1]);

            O2.p = double(psol.pO2); %O2.p = O2.p(1);
            O.p = double(psol.pO); %O.p = O.p(1);
            NO.p = double(psol.pNO); %NO.p = NO.p(1);
            N2.p = double(psol.pN2); %N2.p = N2.p(1);
            N.p = double(psol.pN); %N.p = N.p(1);
            Oplus.p = 0;
            Nplus.p = 0;
            e.p = 0;
        end
        
        if length(psol.pO2) > 1
            O2.p = O2.p(1);
            O.p = O.p(1);
            NO.p = NO.p(1);
            N2.p = N2.p(1);
            N.p = N.p(1);
            Oplus.p = Oplus.p(1);
            Nplus.p = Nplus.p(1);
            e.p = e.p(1);
        end

    end
    %ps = [O2.p O.p Oplus.p NO.p N2.p N.p Nplus.p e.p];
    %sum(ps)

    % Density %
     Mhat = (O2.p/p)*O2.Mhat + (N2.p/p)*N2.Mhat + (N.p/p)*N.Mhat + ...
         (O.p/p)*O.Mhat + (NO.p/p)*NO.Mhat + (Nplus.p/p)*N.Mhat + ...
         (Oplus.p/p)*Oplus.Mhat + (e.p/p)*e.Mhat;
    R = Rhat/Mhat;
    rho = p/(R*T);
    
    % Z %
    Z = p/(rho*287.052874*T);
    
    % Mass Fractions %
    O2.cs = ((O2.p/p)*O2.Mhat)/Mhat;
    O.cs = ((O.p/p)*O.Mhat)/Mhat;
    Oplus.cs = ((Oplus.p/p)*Oplus.Mhat)/Mhat;
    NO.cs = ((NO.p/p)*NO.Mhat)/Mhat;
    N2.cs = ((N2.p/p)*N2.Mhat)/Mhat;
    N.cs = ((N.p/p)*N.Mhat)/Mhat;
    Nplus.cs = ((Nplus.p/p)*N.Mhat)/Mhat;
    e.cs = ((e.p/p)*e.Mhat)/Mhat;
    %cs = [O2.cs O.cs Oplus.cs NO.cs N2.cs N.cs Nplus.cs e.cs];
    %sum(cs)

    % Specific and Non-Dimensional Energy %
    O2.etr = 1.5*O2.R*T;
    O2.evib = (O2.R*O2.thetav)/(exp(O2.thetav/T) - 1);
    O2.erot = O2.R*T; % T >> O2.thetar
    O2.eel = (O2.R*11390*(O2.g(2)/O2.g(1))*exp(-11390/T))/(1 + (O2.g(2)/O2.g(1))*exp(-11390/T));
    O2.e = O2.etr + O2.evib + O2.erot + O2.eel;

    O.etr = 1.5*O.R*T;
    O.eel = (O.R*228*(O.g(2)/O.g(1))*exp(-228/T))/(1 + (O.g(2)/O.g(1))*exp(-228/T));
    O.e0 = O.R*0.5*O2.thetad;
    O.e = O.etr + O.eel + O.e0;

    Oplus.etr = 1.5*Oplus.R*T;
    Oplus.e0 = Oplus.R*(O.thetai + 0.5*O2.thetad);
    Oplus.e = Oplus.etr + Oplus.e0;

    NO.etr = 1.5*NO.R*T;
    NO.evib = (NO.R*NO.thetav)/(exp(NO.thetav/T) - 1);
    NO.erot = NO.R*T;
    NO.eel = (NO.R*174*exp(-174/T))/(1 + exp(-174/T));
    NO.e0 = NO.R*(-NO.thetad + 0.5*N2.thetad + 0.5*O2.thetad);
    NO.e  = NO.etr + NO.evib + NO.erot + NO.eel + NO.e0;

    N.etr = 1.5*N.R*T;
    N.e0 = N.R*0.5*N2.thetad;
    N.e = N.etr + N.e0;

    N2.etr = 1.5*N2.R*T;
    N2.evib = (N2.R*N2.thetav)/(exp(N2.thetav/T) - 1);
    N2.erot = N2.R*T;
    N2.e = N2.etr + N2.evib + N2.erot;

    Nplus.etr = 1.5*Nplus.R*T;
    Nplus.eel = (Nplus.R*70.6*(Nplus.g(2)/Nplus.g(1))*exp(-70.6/T))/(1+(Nplus.g(2)/Nplus.g(1))*exp(-70.6/T));
    Nplus.e0 = Nplus.R*(N.thetai + 0.5*N2.thetad);
    Nplus.e = Nplus.etr + Nplus.eel + Nplus.e0;

    e.e = e.R*1.5*T;

    etotal = O2.e*O2.cs + O.e*O.cs + Oplus.e*Oplus.cs + NO.e*NO.cs + ...
        N2.e*N2.cs + N.e*N.cs + Nplus.e*Nplus.cs + e.e*e.cs;
    EZ = Z*etotal/(R*T);

    % Specific Enthalpy %
    htotal = etotal + R*T;
    
    % Specific and Non-Dimensional Entropy %
    % Translational, vibrational, rotational, electronic

    O2.str = O2.R*(2.5*log(T) - log(O2.p) + log((2*pi*O2.m/h^2)^1.5*k^2.5) + 2.5);
    O2.svib = O2.R*(-log(1 - exp(-O2.thetav/T)) + (O2.thetav/T)/(exp(O2.thetav/T) - 1));
    O2.srot = O2.R*(log(T/(2*O2.thetar)) + 1);
    O2.sel = O2.R*(log(O2.g(1)) + log((1 + (O2.g(2)/O2.g(1))*exp(-11390/T))...
        + ((O2.g(2)/O2.g(1))*(11390/T)*exp(-11390/T))/(1 + ((O2.g(2)/O2.g(1))*exp(-11390/T))^2)));
    O2.s = O2.str + O2.svib + O2.srot + O2.sel;
    
    O.str = O.R*(2.5*log(T) - log(O.p) + (log((2*pi*O.m/h^2)^1.5*k^2.5) + 2.5));
    O.sel = O.R*(log(O.g(1)));% + log(1 + (O.g(2)/O.g(1))*exp(-228/T))...
        %+ ((O.g(2)/O.g(1))*(228/T)*exp(-228/T))/(1 + ((O.g(2)/O.g(1))*exp(-228/T))^2));
    O.s = O.str + O.sel;

    Oplus.str = Oplus.R*(2.5*log(T) - log(Oplus.p) + log((2*pi*Oplus.m/h^2)^1.5*k^2.5) + 2.5);
    Oplus.sel = Oplus.R*(log(Oplus.g(1)));% + log(Oplus.g(1) + Oplus.g(2)*exp(-38600/T))...
        %+ ((Oplus.g(2)/Oplus.g(1))*(38600/T)*exp(-38600/T))/(1 + ((Oplus.g(2)/Oplus.g(1))*exp(-38600/T))^2));
    Oplus.s = Oplus.str + Oplus.sel;

    NO.str = NO.R*(2.5*log(T) - log(NO.p) + log((2*pi*NO.m/h^2)^1.5*k^2.5) + 2.5);
    NO.svib = NO.R*(-log(1 - exp(-NO.thetav/T)) + (NO.thetav/T)/(exp(NO.thetav/T) - 1));
    NO.srot = NO.R*(log(T/(2*NO.thetar)) + 1);
    NO.sel = NO.R*(log(NO.g(1)) + log(NO.g(1) + NO.g(2)*exp(-174/T))...
        + ((NO.g(2)/NO.g(1))*(174/T)*exp(-174/T))/(1 + ((NO.g(2)/NO.g(1))*exp(-174/T))^2));
    NO.s = NO.str + NO.svib + NO.srot + NO.sel;

    N2.str = N2.R*(2.5*log(T) - log(N2.p) + log((2*pi*N2.m/h^2)^1.5*k^2.5) + 2.5);
    N2.svib = N2.R*(-log(1 - exp(-N2.thetav/T)) + (N2.thetav/T)/(exp(N2.thetav/T) - 1));
    N2.srot = N2.R*(log(T/(2*N2.thetar)) + 1);
    N2.s = N2.str + N2.svib + N2.srot;

    N.str = N.R*(2.5*log(T) - log(N.p) + log((2*pi*N.m/h^2)^1.5*k^2.5) + 2.5);
    N.sel = N.R*(log(N.g(1)));% + log(N.g(1) + N.g(2)*exp(-27700/T))...
        %+ ((N.g(2)/N.g(1))*(27700/T)*exp(-27700/T))/(1 + ((N.g(2)/N.g(1))*exp(-27700/T))^2));
    N.s = N.str + N.sel;

    Nplus.str = Nplus.R*(2.5*log(T) - log(Nplus.p) + log((2*pi*Nplus.m/h^2)^1.5*k^2.5) + 2.5);
    Nplus.sel = Nplus.R*(log(Nplus.g(1))); %+ log(Nplus.g(1) + Nplus.g(2)*exp(-70.6/T))...
        %+ ((Nplus.g(2)/Nplus.g(1))*(70.6/T)*exp(-70.6/T))/(1 + ((Nplus.g(2)/Nplus.g(1))*exp(-70.6/T))^2));        
    Nplus.s = Nplus.str + Nplus.sel;

    e.s = e.R*(2.5*log(T) - log(e.p) + log((2*pi*e.m/h^2)^1.5*k^2.5) + 2.5) ...
        + e.R*(2*log(e.g(1)));   

    stotal = O2.s*O2.cs + O.s*O.cs + Oplus.s*Oplus.cs + NO.s*NO.cs + ...
        N2.s*N2.cs + N.s*N.cs + Nplus.s*Nplus.cs + e.s*e.cs; 
    SZ = (Z*stotal)/R;

    % Speed of Sound %
    
    % Specific Heats %
    
    % Effective Gamma %

    flowstate = [Z EZ SZ rho htotal etotal stotal R];
end