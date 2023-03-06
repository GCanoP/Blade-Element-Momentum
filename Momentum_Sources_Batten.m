% Batten Momentum Source Terms. 

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

% Importing Geometric Data for an Specific TSR          
Geometric_Data = readtable('Geometric_Turbine.csv');
Density = 998.2; % Water Density [kg/m^3]
delta_x = 0.0075; % Thickness.

% Calculating the Axial and Tangential Momentum Sources. 
Geometric_Data.Sx(:) = 0; 
Geometric_Data.St(:) = 0; 
for i=2:(Stations)
    Geometric_Data.Sx(i) = (1/(2*delta_x))*(Density)*((Geometric_Data.Sigma(i)+Geometric_Data.Sigma(i-1))/2)...
        *((((Geometric_Data.W(i)+Geometric_Data.W(i-1))/2))^2)*((Geometric_Data.Cx(i)+Geometric_Data.Cx(i-1))/2) ; 
    Geometric_Data.St(i) = (1/(2*delta_x))*(Density)*((Geometric_Data.Sigma(i)+Geometric_Data.Sigma(i-1))/2)...
        *((((Geometric_Data.W(i)+Geometric_Data.W(i-1))/2))^2)*((Geometric_Data.Cy(i)+Geometric_Data.Cy(i-1))/2) ; 
end

% Writting the Momentum Sources. 
writetable(Geometric_Data, 'Annular_Momentum_Sources_Batten.csv');
       