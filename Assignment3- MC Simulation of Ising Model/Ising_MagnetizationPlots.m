close all; clear all;

%% Scatter Plots
JDivT = [0.2:0.1:0.8];
M = [-9.427350e-01,1.361261e+00 ,5.233627e+01,...
    3.662830e+02,3.894057e+02,3.959868e+02,3.984030e+02];
yyaxis left;
plot(JDivT, M, 'r-o','LineWidth', 2);
xlabel('J_b/K_B T'); ylabel('|Mean magnetization|');

sd = [3.403888e+01, 5.418066e+01,1.523472e+02,...
    1.449735e+01,6.180052e+00,3.360630e+00,1.993231e+00];
yyaxis right;
plot(JDivT, sd, 'b-o','LineWidth', 2);
xlabel('J_b/K_B T'); ylabel('Standard Deviation');

%% Histogram
figure()
A = importdata('Magnetization_data.txt', '\n', 0);
histogram(A);
ylabel('MC samples'); xlabel('Magnetization');
xlim([-150 150]);