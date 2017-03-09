 close all; clear all;

fs = 16; % Font Size
%% Scatter Plots
energyData = dlmread('qd1Out.txt');
time = energyData(:,1);
epot = energyData(:,2);
ekin = energyData(:,3);
etot = energyData(:,4);

plot(time, etot, time, epot, time, ekin, 'LineWidth', 2);
xlabel('time (au)','FontSize', fs); ylabel('Energy (au)','FontSize', fs);
legend('Total', 'Potential', 'Kinetic', 'Location','best');
title('Energy Conservation for 1D Square Barrier', 'FontSize', fs);