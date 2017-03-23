 close all; clear all;

fs = 14; % Font Size
%% Scatter Plots for enerygy
%---------------------------------------------------
energyData = dlmread('DensOStates.txt');
eps = energyData(:,1);
dos = energyData(:,2);

plot(eps, dos, 'LineWidth', 2);
xlabel('\epsilon (eV)','FontSize', fs); ylabel('DoS (1/eV)','FontSize', fs);
title('Density of States', 'FontSize', fs);