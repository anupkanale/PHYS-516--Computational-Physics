 close all; clear all;

%% Scatter Plots
vacData = dlmread('mdOutVac.txt');
time = vacData(:,1);
vac = vacData(:,2);
f = fit(time,vac,'exp1');

plot(f, time, vac);% 'b.', time, f);
xlabel('time (LJ unit)'); ylabel('Vel Auto Corr.');