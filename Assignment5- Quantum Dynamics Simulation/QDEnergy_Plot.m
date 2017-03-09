 close all; clear all;

fs = 16; % Font Size
%% Scatter Plots for enerygy
%---------------------------------------------------
energyData = dlmread('qd1Out.txt');
time = energyData(:,1);
epot = energyData(:,2);
ekin = energyData(:,3);
etot = energyData(:,4);

plot(time, etot, time, epot, time, ekin, 'LineWidth', 2);
xlabel('time (au)','FontSize', fs); ylabel('Energy (au)','FontSize', fs);
legend('Total', 'Potential', 'Kinetic', 'Location','best');
title('Energy Conservation for 1D Square Barrier', 'FontSize', fs);

%% Scatter plots for Reflection and Transmission
%---------------------------------------------------
% coeffData = dlmread('qd1Out.txt');
% time = coeffData(:,1);
% refl = coeffData(:,2);
% tran = coeffData(:,3);
% 
% plot(time, refl, time, tran, 'LineWidth', 2);
% grid on;
% xlabel('time (au)','FontSize', fs); ylabel('Probability','FontSize', fs);
% legend('Refletion Coefficient', 'Transmission Coefficient', 'Location','best');
% title('Refletion & Transmission Coefficients', 'FontSize', fs);