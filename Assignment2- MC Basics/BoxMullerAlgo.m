%generate two independent sets of random numbers
close all; clear all;

for ii=1:1000
    r1(ii) = rand;
    r2(ii) = rand;
    
    %Transformation
    zeta1(ii) = sqrt(-2*log(r1(ii)))*cos(2*pi*r2(ii));
    zeta2(ii) = sqrt(-2*log(r1(ii)))*sin(2*pi*r2(ii));
end;
figure()
plot(r1,r2, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r');
xlabel('r_1'); ylabel('r_2');

figure();
plot(zeta1,zeta2, 'bo', 'MarkerSize', 5, 'MarkerFaceColor', 'b');
xlabel('\zeta_1'); ylabel('\zeta_2');
