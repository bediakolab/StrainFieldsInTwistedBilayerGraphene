% critical_amplitude_Tadmor_test.m

a = 1.42;
ampcrit = @(r) 2.*r.*sin(a./(4.*r));

figure;
r = 0:0.01:10;
ampcrits = ampcrit(r);
plot(r,ampcrits);
xlabel('Distance from AA center');
ylabel('Cutoff radius for AA domain');



