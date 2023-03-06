function [a,b] = Get_Induction_Factors(Sigma, Phi, F_total, Cx, Cy)
syms ax at % Symbolic Variables for Axial and Tangential Factors. 
% if (Ct_turbine > 0.4)
% a = real(((18*F_total)-(20)-(3*sqrt((Cx*((50)-(36*F_total)))+(12*F_total*((3*F_total)-(4))))))/((36*F_total)-(50)));
% else
% Solving Axial Induction Factor.
eq_ax =(((ax)/(1-ax))==((Sigma*Cx)/(4*F_total*(sind(Phi))^2)));
a = real(solve(eq_ax, ax));
%end
% Solving Tangential Induction Factor. 
eq_at=((at)/(1+at))==((Sigma*Cy)/(4*F_total*(sind(Phi))*(cosd(Phi))));
b = real(solve(eq_at, at));
end

