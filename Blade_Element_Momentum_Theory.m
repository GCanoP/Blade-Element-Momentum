% Blade Element Momentum Theory - Code Implementation. 
% Version 1.1 - Gerardo Cano Perea - 10/10/2020

% NACA 638-xx Hydrofoil Data
Stations = 8; 
NACA_63812 = readtable('.\Hydrofoil_Data\NACA_63812_Extended.csv');
NACA_63813 = readtable('.\Hydrofoil_Data\NACA_63813_Extended.csv');
NACA_63814 = readtable('.\Hydrofoil_Data\NACA_63814_Extended.csv');
NACA_63815 = readtable('.\Hydrofoil_Data\NACA_63815_Extended.csv');
NACA_63816 = readtable('.\Hydrofoil_Data\NACA_63816_Extended.csv');
NACA_63817 = readtable('.\Hydrofoil_Data\NACA_63817_Extended.csv');
NACA_63818 = readtable('.\Hydrofoil_Data\NACA_63818_Extended.csv');
NACA_63821 = readtable('.\Hydrofoil_Data\NACA_63821_Extended.csv');

% NACA 638-xx Hydrofoil Polar Data Objects, 
Re = 50000;  % Reynolds Number [-]
NACA_63812_Polar = Polar(Re, NACA_63812.AoA, NACA_63812.CL, NACA_63812.CD);
NACA_63813_Polar = Polar(Re, NACA_63813.AoA, NACA_63813.CL, NACA_63813.CD);
NACA_63814_Polar = Polar(Re, NACA_63814.AoA, NACA_63814.CL, NACA_63814.CD);
NACA_63815_Polar = Polar(Re, NACA_63815.AoA, NACA_63815.CL, NACA_63815.CD);
NACA_63816_Polar = Polar(Re, NACA_63816.AoA, NACA_63816.CL, NACA_63816.CD);
NACA_63817_Polar = Polar(Re, NACA_63817.AoA, NACA_63817.CL, NACA_63817.CD);
NACA_63818_Polar = Polar(Re, NACA_63818.AoA, NACA_63818.CL, NACA_63818.CD);
NACA_63821_Polar = Polar(Re, NACA_63821.AoA, NACA_63821.CL, NACA_63821.CD);
Hydrofoils = [NACA_63821_Polar, NACA_63818_Polar, NACA_63817_Polar, NACA_63816_Polar,...
              NACA_63815_Polar, NACA_63814_Polar, NACA_63813_Polar, NACA_63812_Polar];

% Importing Geometric Turbine Data
% Airfoil Name - Mu [-] - Chord [m] - Twist [deg] 
Geometric = readtable('.\Geometric_Turbine\Geometric_Turbine.csv');

% Turbine Static and Dynamic Variables. 
Turbine = Turbine();
Turbine.Blades = 3; % Blades[-]
Turbine.Radio = 0.15; % Radio[m]
Turbine.Radio_Hub = 0.03; % Radio_Hub[m]
Turbine.U_free = 0.7; % U_free[m/s]
Turbine.TSR = 6;  % TSR[-]
Turbine.Omega = Turbine.TSR*Turbine.U_free...
    /Turbine.Radio; % Angular Speed [rad/seg]

% Dynamic Geometric Variables
Geometric.r = Geometric.Mu.*Turbine.Radio;
Geometric.TSR_local = Geometric.r.*Turbine.Omega/Turbine.U_free;
% Change the Value in Last Station to Avoid an Error in the Loss Equations.
Geometric.r(Stations) = 0.95*Turbine.Radio;

% Preallocating Variables. 
Geometric.Phi(:) = 0;
Geometric.Alpha(:)= 0;
Geometric.Beta = Geometric.Twist;
Geometric.Sigma(:) = 0;
Geometric.CL(:) = 0; 
Geometric.CD(:) = 0; 
Geometric.AeroEfi(:)= 0;
Geometric.Cx(:) = 0; 
Geometric.Cy(:) = 0; 
Geometric.F_tip(:) = 0;
Geometric.F_root(:) = 0;
Geometric.F_total(:) = 0; 
Geometric.a(:) = 0;
Geometric.b(:) = 0; 
Geometric.a_new(:) = 0;
Geometric.b_new(:) = 0; 
Geometric.error_a(:) = 5; % Initial Error [5%]
Geometric.error_b(:) = 5; % Initial Error [5%]
Geometric.W(:) = 0; % Relative Velocity 

% Blade Element Momentum Theory - Iterative Calculation.
for i=1:(Stations)
    while(Geometric.error_b(i)>1) % Relative Error 1%
        % Calculating the Inflow Angle at the i-th Station.
        Geometric.Phi(i) = atand((1-Geometric.a(i))/((1+Geometric.b(i))*Geometric.TSR_local(i)));
        Geometric.Alpha(i) = Geometric.Phi(i) - Geometric.Beta(i);
        % Calculating the Blade Solidity at the i-th Station.
        Geometric.Sigma(i) = ((Turbine.Blades*Geometric.Chord(i))/(2*pi*Geometric.r(i)));
        % Claculating the Polar Coefficients at the i-th Station.
        Geometric.CL(i) = Hydrofoils(i).Interp_CL(Geometric.Alpha(i));
        Geometric.CD(i) = Hydrofoils(i).Interp_CD(Geometric.Alpha(i));
        Geometric.AeroEfi(i) = Geometric.CL(i)/Geometric.CD(i);
        % Calculating the Normal and Tangential Components at the i-th Station.
        Geometric.Cx(i) = (Geometric.CL(i)*cosd(Geometric.Phi(i)))+(Geometric.CD(i)*sind(Geometric.Phi(i)));
        Geometric.Cy(i) = (Geometric.CL(i)*sind(Geometric.Phi(i)))-(Geometric.CD(i)*cosd(Geometric.Phi(i)));
        % Calculating the Tip - Root Loss Factor
        Geometric.F_tip(i) = (2/pi).*acos(exp(-(((Turbine.Blades/2)*(1-(Geometric.r(i)/Turbine.Radio)))/...
            ((Geometric.r(i)/Turbine.Radio).*sind(Geometric.Phi(i))))));
        Geometric.F_root(i) = ((2/pi)*acos(exp(-((Turbine.Blades/2)*((Geometric.r(i)-Turbine.Radio_Hub)/...
            (Geometric.r(i)*sind(Geometric.Phi(i))))))));
        Geometric.F_total(i) = Geometric.F_root(i)*Geometric.F_tip(i);
        % Calculating the Axial and Tangential Induction Factors. 
        [Geometric.a_new(i), Geometric.b_new(i)] = Get_Induction_Factors(Geometric.Sigma(i), Geometric.Phi(i),...
            Geometric.F_total(i), Geometric.Cx(i), Geometric.Cy(i));
        % Calculating the Relative Errors.
        Geometric.error_a(i) = ((abs(Geometric.a(i)-Geometric.a_new(i)))/Geometric.a_new(i))*100;
        Geometric.error_b(i) = ((abs(Geometric.b(i)-Geometric.b_new(i)))/Geometric.b_new(i))*100;
        % New Variables. 
        Geometric.a(i) = Geometric.a_new(i);
        Geometric.b(i) = Geometric.b_new(i);
        % Calculating the Relative Velocity 
        Geometric.W(i) = ((Turbine.U_free*(1-Geometric.a(i)))^2+(Turbine.Omega*Geometric.r(i)*...
            (1+Geometric.b(i)))^2)^(1/2);
    end 
end 

% Calculating the Dynamic Loads Along Blade.
Turbine.Rho = 998.2; % Fluid Density [kg/m^3]
Blade_Loads = readtable('.\Blade_Loads\Blade_Loads.csv'); 
Geometric.r(Stations) = 1*Turbine.Radio;

% Preallocating Variables. 
% New Variables - Extended Vectors.
Blade_Loads.a(:)=0;
Blade_Loads.b(:)=0;
Blade_Loads.F_total(:)=0;
Blade_Loads.Cx(:)=0;
Blade_Loads.Cy(:)=0;
Blade_Loads.W(:)=0;
Blade_Loads.Chord(:)=0;
Blade_Loads.Phi(:)=0;
for j=1:(length(Blade_Loads.r))
    Blade_Loads.a(j) = interp1(Geometric.r, Geometric.a, Blade_Loads.r(j));
    Blade_Loads.b(j) = interp1(Geometric.r, Geometric.b, Blade_Loads.r(j));
    Blade_Loads.F_total(j) = interp1(Geometric.r, Geometric.F_total, Blade_Loads.r(j));
    Blade_Loads.Cx(j) = interp1(Geometric.r, Geometric.Cx, Blade_Loads.r(j));
    Blade_Loads.Cy(j) = interp1(Geometric.r, Geometric.Cy, Blade_Loads.r(j));
    Blade_Loads.W(j) = interp1(Geometric.r, Geometric.W, Blade_Loads.r(j));
    Blade_Loads.Chord(j) = interp1(Geometric.r, Geometric.Chord, Blade_Loads.r(j));
    Blade_Loads.Phi(j) = interp1(Geometric.r, Geometric.Phi, Blade_Loads.r(j));
end

% Calculating Thrust.
Blade_Loads.Q(:) = 0;
Blade_Loads.T(:) = 0;
for k=1:(length(Blade_Loads.r)-1)
    lb = Blade_Loads.r(k);
    ub = Blade_Loads.r(k+1);
    fun_1 = @(x) 4.*pi.*x.*Turbine.Rho.*Turbine.U_free.^2.*Blade_Loads.a(k).*...
        (1.-Blade_Loads.a(k)).*Blade_Loads.F_total(k);
    fun_2 = @(x) 4.*pi.*x.^3.*Turbine.Rho.*Turbine.U_free.*Turbine.Omega.*...
        Blade_Loads.b(k).*(1.-Blade_Loads.a(k)).*Blade_Loads.F_total(k);
%     fun_2 = @(x) 0.5.*Turbine.Rho.*Blade_Loads.W(k)^2.*Turbine.Blades.*Blade_Loads.Chord(k).*x.*Blade_Loads.Cy(k);
    Blade_Loads.T(k) = integral(fun_1, lb, ub);
    Blade_Loads.Q(k) = integral(fun_2, lb, ub);
end

% Calculating the Thrust and Power Coefficients.
Turbine.Cp = ((sum(Blade_Loads.Q)*Turbine.Omega)/(0.5*Turbine.Rho*pi*Turbine.Radio^2*Turbine.U_free^3));
Turbine.Ct = ((sum(Blade_Loads.T))/(0.5*Turbine.Rho*pi*Turbine.Radio^2*Turbine.U_free^2));

disp(Turbine);

% Writting Specific Turbine File.
Turbine_File = 'Turbine_Parameters.txt';
fid = fopen(Turbine_File,'w+');
fprintf(fid, 'Turbine Parameters.');
fprintf(fid, '\nBlades = %d [-]', Turbine.Blades);
fprintf(fid, '\nRadio = %d [m]', Turbine.Radio);
fprintf(fid, '\nRadio Hub = %d [m]', Turbine.Radio_Hub);
fprintf(fid, '\nTip Speed Ratio = %d [-]', Turbine.TSR);
fprintf(fid, '\nFree-Stream Velocity = %d [m/s]', Turbine.U_free);
fprintf(fid, '\nAngular Velocity = %d [rad/seg]', Turbine.Omega);
fprintf(fid, '\nDensity  = %d [kg/m^3]', Turbine.Rho);
fprintf(fid, '\nPower Coefficient = %d [-]', Turbine.Cp);
fprintf(fid, '\nThrust Coefficient = %d [-]', Turbine.Ct);
fclose(fid);

% Writting Geometric and Loads Data Tables.
writetable(Geometric,'Geometric_Turbine.csv');
writetable(Blade_Loads,'Blade_Loads.csv');




