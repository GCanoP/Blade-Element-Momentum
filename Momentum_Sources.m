% Calculating Momentum Sources.

% Importing Blade Loads Table. 
Blade_Loads = readtable('Blade_Loads.csv');  
Blade_Loads.Sa(:) = 0; % Axial Source Term [N]
Blade_Loads.St(:) = 0; % Tangential Source Term [N]
Blade_Loads.dA(:) = 0; % Differential Area [m^2]

% Turbine Static and Dynamic Variables. 
Turbine = Turbine();
Turbine.Blades = 3; % Blades[-]
Turbine.Radio = 0.15; % Radio[m]
Turbine.Radio_Hub = 0.03; % Radio_Hub[m]
Turbine.U_free = 0.7; % U_free[m/s]
Turbine.TSR = 5;  % TSR[-]
Turbine.Rho = 998.2; % Fluid Density [kg/m^3]
Turbine.Omega = Turbine.TSR*Turbine.U_free...
    /Turbine.Radio; % Angular Speed [rad/seg]

% Calculating the Force Sources for CFD-BEM Model [N]. 
for i=1:(length(Blade_Loads.Mu)-1)
    % Calculating Axial Force Source - [N]
    Blade_Loads.Sa(i) = ((((0.5)*(Turbine.Rho)*(Blade_Loads.W(i)^2)*(Turbine.Blades)*...
        (Blade_Loads.Chord(i))*(Blade_Loads.Cx(i))*(Blade_Loads.r(i+1)))-((0.5)*...
        (Turbine.Rho)*(Blade_Loads.W(i)^2)*(Turbine.Blades)*(Blade_Loads.Chord(i))*...
        (Blade_Loads.Cx(i))*(Blade_Loads.r(i)))));
    % Calculating Tangential Force Source - [N]
    Blade_Loads.St(i) = (((0.5)*(Turbine.Rho)*(Blade_Loads.W(i)^2)*(Turbine.Blades)*...
        (Blade_Loads.Chord(i))*(Blade_Loads.Cy(i))*(Blade_Loads.r(i+1)))-((0.5)*...
        (Turbine.Rho)*(Blade_Loads.W(i)^2)*(Turbine.Blades)*(Blade_Loads.Chord(i))*...
        (Blade_Loads.Cy(i))*(Blade_Loads.r(i))));
    % Claculating the Differential Area. 
    Blade_Loads.dA(i) = pi*((Blade_Loads.r(i+1))^2-(Blade_Loads.r(i))^2);
end
% Momentum Sources - Parameters. 
delta_x = 0.0075; % Actuator Disck Thickness [m] 
No_discs = (length(Blade_Loads.Mu)-1); % Number of annulus discs [-] 

% Creating Momentum Source Table 
Annular_Momentum = readtable('.\Momentum_Source\Momentum_Sources.csv');
Annular_Momentum.Axial_Force(:) = 0; % Annular Axial Force [N]  
Annular_Momentum.Tangent_Force(:) = 0; % Annular Tangential Force [N]
Annular_Momentum.Area(:) = 0; % Annular Area [m^2] 
Annular_Momentum.Volumen(:) = 0; % Annular Volumen [m^3]
Annular_Momentum.Sx(:) = 0; % Annular Axial Momentum Source [N/m^3]   
Annular_Momentum.St(:) = 0; % Annular Tangential Momentum Source [N/m^3]
i=1; % Iterator i

% % Calculating the Momentum Sources for CFD-BEM Model [N]. 
for j=1:10:No_discs
     Annular_Momentum.Axial_Force(i) = sum(Blade_Loads.Sa(j:(j+9)));
     Annular_Momentum.Tangent_Force(i) = sum(Blade_Loads.St(j:(j+9)));
     Annular_Momentum.Area(i) = sum(Blade_Loads.dA(j:(j+9)));
     Annular_Momentum.Volumen(i) = ((Annular_Momentum.Area(i))*(delta_x));
     % Calculating the Axial Momentum Source [N/m^3] 
     Annular_Momentum.Sx(i) = ((Annular_Momentum.Axial_Force(i))/(Annular_Momentum.Volumen(i)));
     % Calculating the Tangential Momentum Source [N/m^3] 
     Annular_Momentum.St(i) = ((Annular_Momentum.Tangent_Force(i))/(Annular_Momentum.Volumen(i)));
     i=i+1;
end
% 
disp(sum(Blade_Loads.T));
disp(sum(Blade_Loads.Sa));
disp(sum(Annular_Momentum.Axial_Force));

disp(sum(Blade_Loads.Q));
disp(sum(Blade_Loads.St));
disp(sum(Annular_Momentum.Tangent_Force));

writetable(Annular_Momentum,'Annular_Momentum_Sources.csv');
