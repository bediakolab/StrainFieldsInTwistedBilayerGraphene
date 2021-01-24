function [log_lik] = getFullRadialPoissonMaxLikFit_bootstrap(c,linear_data)
% Implementation of maximum likelihood estimation
% Uses bootstrap

% Needs to obtain the power law structure as before.
calc_poisson_means = getFullRadialPowerLawPred_bootstrap(c,linear_data);
data_to_use = linear_data;
% Because we will have applied the mask filtering outside of the bootstrp
% call, there will be no need to filter in here.
% removal_indices = (linear_data(:,1) < 0);  % because now column 1 gives the intensities.
% data_to_use(removal_indices) = [];
% calc_poisson_means(removal_indices) = [];
% logliks = log(poisspdf(data_to_use(:,1),calc_poisson_means));
% log_lik = -sum(sum(logliks));  % function maximizes by default; hence throw in the negative.

%% The improved version of the Poisson loglik code.
k = data_to_use(:,1);
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

