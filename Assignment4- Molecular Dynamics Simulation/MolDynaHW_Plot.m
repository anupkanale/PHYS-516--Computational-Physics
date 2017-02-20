 close all; clear all;

%% Scatter Plots
eulerData = dlmread('mdOutEuler.txt');
verletData = dlmread('mdOutVerlet.txt');
time = eulerData(:,1);
energyEuler = eulerData(:,4);
energyVerlet = verletData(:,4);

plot(time, energyEuler, 'r-', time, energyVerlet, 'b--','LineWidth', 2);
xlabel('time (LJ unit)'); ylabel('Energy (LJ unit)');
legend('Euler', 'Velocity-Verlet', 'Location','best');