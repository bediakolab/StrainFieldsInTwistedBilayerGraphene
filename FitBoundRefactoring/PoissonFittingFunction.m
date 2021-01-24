function [ log_lik ] = PoissonFittingFunction( lambda,values )
% Uses code developed over the summer on the Bediako lab computers

% explicitly compute the logarithm of the poisson pmf
k = values;
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

