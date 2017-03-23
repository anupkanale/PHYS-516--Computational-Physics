 close all; clear all;

fs = 14; % Font Size
%% Scatter Plots for enerygy
%---------------------------------------------------
energyData = dlmread('DensOStates.txt');
eps = energyData(:,1);
dos = energyData(:,2);

figure()
plot(eps, dos, 'LineWidth', 1.5);
xlabel('\epsilon (eV)','FontSize', fs); ylabel('DoS (1/eV)','FontSize', fs);
title('Density of States', 'FontSize', fs);

energyData = dlmread('fermiDist.txt');
d = energyData(:,1);
ferm = energyData(:,2);

figure()
plot(d, ferm, 'ro', 'LineWidth', 1.5);
xlim([d(1) d(end)])
xlabel('\epsilon (eV)','FontSize', fs); ylabel('Ferm','FontSize', fs);
title('Fermi Distribution', 'FontSize', fs);