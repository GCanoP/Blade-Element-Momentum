% Graphing Turbine Performance. 

% Graph_Performance = readtable('Turbine_Performance.csv');
% fig = figure;
% color_1 = [(39/256) (150/256) (143/256)];
% color_2 = [(74/256) (110/256) (162/256)];
% set(fig,'defaultAxesColorOrder',[color_1; color_2])
% yyaxis left  
% plot(Graph_Performance.TSR, Graph_Performance.Cp,'Color',color_1);
% ylim([0 0.375]);
% yyaxis right
% plot(Graph_Performance.TSR, Graph_Performance.Ct,'Color',color_2);
% ylim([0.125 1]);
% xlim([1 6]);
% grid on  
% grid minor 
% title('Parametros de Rendimiento')
% xlabel('TSR - \bf\lambda');
% ylabel('Coeficiente de Empuje - C_{T}');
% yyaxis left
% ylabel('Coeficiente de Potencia - C_{P}');

% Graph_Loads = readtable('.\Turbine_for_CFD\TSR_4\Blade_Loads.csv');
% fig = figure;
% color_1 = [(39/256) (150/256) (143/256)];
% color_2 = [(74/256) (110/256) (162/256)];
% set(fig,'defaultAxesColorOrder',[color_1; color_2])
% yyaxis left  
% plot(Graph_Loads.Mu(1:70), Graph_Loads.T(1:70),'Color',color_1);
% ylim([0.0625 0.24]);
% yyaxis right
% plot(Graph_Loads.Mu(1:70), Graph_Loads.Q(1:70),'Color',color_2);
% ylim([0.001 0.01]);
% xlim([0.3 1]);
% grid on  
% grid minor 
% title('Cargas Hidrodinámicas sobre el Alabe (\bf\lambda = 4)')
% xlabel('Distancia Radial \bf\mu');
% ylabel('Torque [Nm]');
% yyaxis left
% ylabel('Empuje [N]');


