classdef Turbine
    % Turbine Static and Dynamic Variables
    properties
        Blades % Number of Blades [-]
        Radio % Rotor Radius [m]
        Radio_Hub % Hub Radius [m]
        TSR % Tip Speed Ratio [-] 
        U_free % Free-Stream Velocity [m/s] 
        Omega % Angular Velocity [rad/seg]  
        Rho % Fluid Density [kg/m^3] 
        Cp % Power Coefficient [-] 
        Ct % Thrust Coefficient [-] 
    end
end

