clc;
close all; clear all;

N = 10000;
fs = 14; % Font Size
%% Scatter Plots for Stock prices
%---------------------------------------------------
for walker=1:N
    filename = strcat('stockData_walker', int2str(walker), '.txt');
    stockData = dlmread(filename);
    day = stockData(:,1);
    price = stockData(:,2);
    plot(day, price, 'LineWidth', 1.5);
    hold on;
    delete(filename);
end

xlabel('day','FontSize', fs); ylabel('Stock price','FontSize', fs);
title('Stock price Monte Carlo simulation', 'FontSize', fs-4);
hold off;

%% Histogram data
%---------------------------------------------------
histData = dlmread('histogramData.txt');
investor = histData(:,1);
count = histData(:,2);
plot(investor, count, 'LineWidth', 1.5);

xlabel('investor','FontSize', fs); ylabel('Count','FontSize', fs);
title('Histogram', 'FontSize', fs-4);
