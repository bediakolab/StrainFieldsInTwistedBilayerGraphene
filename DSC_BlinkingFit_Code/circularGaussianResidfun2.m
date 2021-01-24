function [ residuals ] = circularGaussianResidfun2( c,space,fitdata )
% 04/28/2020
%
% fitdata is an nx1 column of the displacement amplitude values we are
% trying to fit.
%
% c(1) is xmean
% c(2) is ymean
% sigma1 is x marginal var
% sigma2 is y marginal var
% p_covar is covariance weight, to be bounded in the range [-1,1].
% A is vertical translation constant
% B is undetermined prefactor scaling

[ gausspred ] = circularGaussianPredfun2( space,c );
residuals = fitdata - gausspred;

end

