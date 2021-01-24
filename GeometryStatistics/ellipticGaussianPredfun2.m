function [ gausspred ] = ellipticGaussianPredfun2( space,c )
% Nathanael Kazmierczak, 04/28/2020
%
% Written after the initial attempt to do this with trig rotations, which
% didn't work out for whatever reason.
%
% A is a vertical translation constant, while B is an unspecified
% prefactor.

% Unpack variables
xmean = c(1);
ymean = c(2);
FWHM1 = c(3);
FWHM2 = c(4);
p_covar = c(5);
A = c(6);
B = c(7);

sigma1 = FWHM1/(2*sqrt(2*log(2)));
sigma2 = FWHM2/(2*sqrt(2*log(2)));

% Construct covariance matrix for multivariate normal
S = zeros(2);
S(1,1) = sigma1^2;
S(2,2) = sigma2^2;
S(1,2) = p_covar*sigma1*sigma2;
S(2,1) = S(1,2);
mus = [xmean,ymean];

% space, the coordinates, needs to be an nx2 array with x in first col and
% y in second (same as before, supports deletion of points on the basis of
% saddle point filter or whatever else we want.
% x = space(:,1);
% y = space(:,2);
gauss = mvnpdf(space,mus,S);
gausspred = A - B*gauss./max(gauss); % B is now the true scaling.

end

