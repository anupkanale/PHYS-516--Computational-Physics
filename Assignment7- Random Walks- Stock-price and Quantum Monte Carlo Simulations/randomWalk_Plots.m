clc;
close all; clear all;

fs = 14; % Font Size
%% Scatter Plots for enerygy
%---------------------------------------------------
energyData = dlmread('stockData.txt');
day = energyData(:,1);
price = energyData(:,2);

figure()
plot(day, price, 'LineWidth', 1.5);
xlabel('day','FontSize', fs); ylabel('Stock price','FontSize', fs);
title('Stock price Monte Carlo simulation', 'FontSize', fs-4);