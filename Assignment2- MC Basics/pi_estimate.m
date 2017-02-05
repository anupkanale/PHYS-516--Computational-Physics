close all; clear all;
npts = 10.^[1:7]; % number of points for integration

% number of trials
for ii=1:length(npts)
    
    sum = 0; sum2 = 0;
    for jj=1:npts(ii)
        x = rand;
        fx = 4/(1+x^2);
        sum = sum+fx;
        sum2 = sum2+fx^2;
    end;
    pi(ii) = sum/npts(ii);
    pi2(ii) = sum2/npts(ii);
    stdv(ii) = (pi2(ii)-pi(ii)^2)/(npts(ii)-1);
    stdv(ii) = sqrt(stdv(ii));
end;

%% Plots
figure, set(gcf,'color','w'), hold on, box on
errorbar(log10(npts),pi,stdv, 'LineWidth', 2.5);
grid on;
hold on;
plot(log10(npts),3.14*ones(length(npts),1), '--r', 'LineWidth', 2.5);
xlabel('log M'); ylabel('\pi Estimate');

figure()
loglog(npts, stdv, 'o')