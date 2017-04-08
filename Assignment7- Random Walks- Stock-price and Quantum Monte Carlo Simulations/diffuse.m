close all;
clear; clc;

N_walker = 1000;
maxStep = 100;
position = zeros(maxStep,1);
hist = zeros(1001,1);

for walker=1:N_walker
    for step=2:maxStep
        if (rand(1)>0.5), position(step) = position(step-1)+1;
        else, position(step) = position(step-1)-1;
        end
    end
    k(step) = position(step)+500;
    hist(k(step)) = hist(k(step))+1;
end

x = -500:500;
plot(x,hist);