% testingStirling.m
%
% Because I am a bit worried about losing precision through log(factorial)
% calculations.

values = 1:200;
storage = zeros(200,2);
for i = values
    storage(i,1) = log(factorial(i));
    storage(i,2) = i*log(i) - i;
end

figure;
plot(values',storage(:,1),'k',values',storage(:,2),'r');
legend('Explicit form','Stirlings approximation');
title('Computation of ln(n!)');

stirling_residuals = storage(:,1) - storage(:,2);
figure
plot(values',stirling_residuals);
legend('Absolute residuals of Stirlings approximation');
title('Computation of ln(n!)');

relative_stirling_residuals = stirling_residuals./storage(:,1);
figure
semilogy(values',relative_stirling_residuals);
legend('Relative residuals of Stirlings approximation');
title('Computation of ln(n!)');


