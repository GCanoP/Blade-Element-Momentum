classdef Polar
    % Airfoil Polar Object.
    properties (Access = public) 
        Re % Reynolds Number [-] 
        AoA % Angle of Attack [deg]
        CL % Coefficient of Lift [-]
        CD % Coefficient of Drag [-]
    end
    methods
        % Method 1. Defining the Constructor.
        function obj = Polar(Re, AoA, CL, CD)
            if nargin == 4 % Number of Arguments.
                obj.Re = Re;
                obj.AoA = AoA;
                obj.CL = CL; 
                obj.CD = CD;
            end 
        end
        % Method 2. Interpolation for Lift Coefficients. 
        function Interp_CL = Interp_CL(obj, alpha) 
        Interp_CL = interp1(obj.AoA, obj.CL, alpha);
        end 
        % Method 3. Interpolation for Drag Coefficients. 
        function Interp_CD = Interp_CD(obj, alpha)
            Interp_CD = interp1(obj.AoA, obj.CD, alpha);
        end  
    end  
end

