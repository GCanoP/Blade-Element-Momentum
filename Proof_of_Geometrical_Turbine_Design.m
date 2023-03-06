% Program to Probe the Geometrical Design of a Turbine. 
% ========================================
% Basic Data: NACA 63815
% ========================================
U_inf = 0.7; % Free-Stream Velocity [m/s]
No_Blades = 3; % Number of Blades [-]
Turbine_Radio = 0.15; % Turbine Radio [m] 
Hub_Radio = 0.03;  % Hub Radio [m]
TSR = 4; % Tip Speed Ratio [-] 
Density = 998.2; % Water Density [kg/m^3]
Omega = (TSR*U_inf)/Turbine_Radio; % Angular Speed [rad/seg]
Mu = 0.7; % Dimensionless Radius [-] 
r = Turbine_Radio * Mu; % Local Radius [-]
AoA = linspace(0,20,201);
Chord = 0.024;
Bs = (No_Blades*Chord)/(2*pi*r);
TSR_local = (r*Omega)/U_inf;

NACA_63815 = readtable('.\Hydrofoil_Data\NACA_63815.csv');
NACA_63815.Aero = NACA_63815.CL./NACA_63815.CD; 
Aero = max(NACA_63815.Aero);
%Angle = NACA_63815.AoA(find(NACA_63815.Aero == Aero));
% ========================================
% Initilizing: BEM method
% ========================================
a_ini = 0;
b_ini = 0;
