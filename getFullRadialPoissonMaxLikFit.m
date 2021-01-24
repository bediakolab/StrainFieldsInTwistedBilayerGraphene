function [log_lik] = getFullRadialPoissonMaxLikFit(c,fitdata)
% Implementation of maximum likelihood estimation

% Needs to obtain the power law structure as before.
matsize = size(fitdata);
calc_poisson_means = getFullRadialPowerLawPred(c,matsize);
data_to_use = fitdata;
removal_indices = (fitdata < 0);
data_to_use(removal_indices) = [];
calc_poisson_means(removal_indices) = [];
% logliks = log(poisspdf(data_to_use,calc_poisson_means));
% % logliks(logliks < -1000) = -1000;  % To fix my problem with ill conditioning.
% % % super hacky, and the right way is to expand the poisson distribution to
% % % reduce numerical error. Get around to it at some point
k = data_to_use;
lambda = calc_poisson_means;
% explicitly compute the logarithm of the poisson pmf
logliks = k.*log(lambda) - lambda - log(factorial(k));
% If the calculation fails with the explicit formula, try again using
% Stirling's approximation.
% Note to self: do some digging to see if Stirling's formula can be
% computationally superior even when the explicit formula may be within the
% overflow bounds. I'm guessing not.
use_stirling = isinf(logliks);
k_us = k(use_stirling);
lambda_us = lambda(use_stirling);
logliks(use_stirling) = k_us.*log(lambda_us) - lambda_us + k_us - k_us.*log(k_us);
log_lik = -sum(sum(logliks));  % function maximizes by default; hence throw in the negative.

end

