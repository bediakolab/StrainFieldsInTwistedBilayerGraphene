function [ gausspred ] = circularGaussianPredfun2( space,c )
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
FWHM = c(3);
A = c(4);
B = c(5);

sigma = FWHM/(2*sqrt(2*log(2)));



% space, the coordinates, needs to be an nx2 array with x in first col and
% y in second (same as before, supports deletion of points on the basis of
% saddle point filter or whatever else we want.)
% Turn these into radii for the fit.
x = space(:,1);
y = space(:,2);
r = sqrt((x-xmean).^2 + (y-ymean).^2);
gauss = exp(-0.5*(r./sigma).^2);
gausspred = A - B*gauss; % B is now the true scaling.

end

