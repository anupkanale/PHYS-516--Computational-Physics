close all; clear all;
ntry = [10,1e2,1e3,1e4,1e5,1e6];
for ii=1:length(ntry)
    sum = 0; sum2 = 0;
    for jj=1:ntry(ii)
        x = rand;
        fx = 4/(1+x^2);
        sum = sum+fx;
        sum2 = sum2+fx^2;
    end;
    pi(ii) = sum/ntry(ii);
    pi2(ii) = sum2/ntry(ii);
    stdv(ii) = (pi2(ii)-pi(ii)^2)/(ntry(ii)-1);
end;
errorbar(log10(ntry),pi,stdv);
