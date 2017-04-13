% Quantum Monte Carlo Plots
%----------------------------------
clc;
close all; clear;
fs = 14; % Font Size

% Energy Plot
%------------------------------------
energyData = dlmread('vRef.txt');
mcstep = energyData(:,1);
energy = energyData(:,2);
figure(1);
plot(mcstep, energy,'LineWidth', 1.5);

xlabel('MC step','FontSize', fs); ylabel('Energy','FontSize', fs);
title('Enerygy variation', 'FontSize', fs-4);

% Histogram
%-----------------------
histData = dlmread('qmcHistogram.txt');
x = histData(:,1);
psi = histData(:,2);
figure(2);
plot(x, psi,'LineWidth', 1.5);

xlabel('x','FontSize', fs); ylabel('\psi','FontSize', fs);
title('Histogram', 'FontSize', fs-4);
