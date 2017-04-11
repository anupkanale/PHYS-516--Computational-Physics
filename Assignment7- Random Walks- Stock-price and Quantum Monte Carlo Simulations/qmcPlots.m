 clc;
close all; clear;
fs = 14; % Font Size

%% Quantum Monte Carlo PLots
%------------------------------------
energyData = dlmread('vRef.txt');
mcstep = energyData(:,1);
energy = energyData(:,2);
figure();
plot(mcstep, energy,'LineWidth', 1.5);

xlabel('MC step','FontSize', fs); ylabel('Energy','FontSize', fs);
title('Enerygy variation', 'FontSize', fs-4);