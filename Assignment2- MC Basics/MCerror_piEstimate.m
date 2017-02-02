ntry = [10,1e2,1e3,1e4,1e5,1e6];
Nseed=100;
for ii=1:length(ntry)
    sum = 0; sum2 = 0;
    term1 = 0; term2 = 0;
    for kk=1:Nseed
        for jj=1:ntry(ii)
            x = rand;
            fx = 4/(1+x^2);
            sum = sum+fx;
            sum2 = sum2+fx^2;
        end;
        pi(kk) = sum/(ntry(ii)*Nseed);
        term1 = term1 + pi(kk);
        pi2(kk) = sum2/(ntry(ii)*Nseed);
        term2 = term2 + pi2(kk);
    end;
    stdvExp(ii) = term2/Nseed - (term1/Nseed)^2;
    stdvExp(ii) = sqrt(stdvExp(ii));
end;
plot(log10(ntry), log10(stdvExp));
