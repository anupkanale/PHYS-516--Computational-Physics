 close all; clear all;

%% Scatter Plots
vacData = dlmread('mdOutVac.txt');
time = vacData(:,1);
vac = vacData(:,2);

plot(time, vac, 'r-', 'LineWidth',2);
xlabel('time (LJ unit)'); ylabel('Vel Auto Corr.');