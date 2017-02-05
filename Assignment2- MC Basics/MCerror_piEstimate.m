close all;
npts = 10.^[1:6]; % number of points for integration
Nseed = 100;

% number of trials
for ii=1:length(npts)
    %number of repetitions
    term1=0; term2=0;
    for kk=1:Nseed
        %MC integration
        sum = 0;
        for jj=1:npts(ii)
            x = rand;
            fx = 4/(1+x^2);
            sum = sum+fx;
        end;
        pi(kk,ii) = sum/npts(ii);
        term1 = term1 + pi(kk,ii)*pi(kk,ii);
        term2 = term2 + pi(kk,ii);
    end;
    stdvMeas(ii) = term1/Nseed - (term2/Nseed)^2;
    stdvMeas(ii) = sqrt(stdvMeas(ii));
end;

%% Plots
plot(log10(npts),log10(stdvMeas),'ob');
h = lsline
p2 = polyfit(get(h,'xdata'),get(h,'ydata'),1)
grid on;
xlabel('log M'); ylabel('log SD')
hold on;
plot(log10(npts),log10(stdv), 'or')
legend('measured \sigma', 'LS Fit', 'unbiased estimate')