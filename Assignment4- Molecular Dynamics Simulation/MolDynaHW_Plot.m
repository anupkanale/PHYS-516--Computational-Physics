close all; clear all;

%% Scatter Plots
fileID = fopen('mdOutEuler', 'r');
formatSpec = '%d %f';
sizeA = [2 Inf];
[time, energy] = [];
fclose(fileID);

plot(time, energyEuler, 'r-', time, energyVerlet, 'b--','LineWidth', 2);
xlabel('time (LJ unit)'); ylabel('Energy (LJ unit)');