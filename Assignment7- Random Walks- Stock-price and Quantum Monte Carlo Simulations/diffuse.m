close all; clear; clc;
N_walker = 1;
maxStep = 500;

for walker=1:N_walker
    position = 0;
    for step=2:maxStep
        if (rand(1)>0.5), position = position+1;
        else, position = position-1;
        end
        plot(step, position,'k*');
        hold on;
    end
end

