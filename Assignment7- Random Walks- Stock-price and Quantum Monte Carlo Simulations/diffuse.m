close all;
clear; clc;

N_walker = 1;
maxStep = 500;
position = zeros(maxStep,1);

for walker=1:N_walker
    for step=2:maxStep
        if (rand(1)>0.5), position(step) = position(step-1)+1;
        else, position(step) = position(step-1)-1;
        end
    end
end