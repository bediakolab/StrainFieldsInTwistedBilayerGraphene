function [predmat] = getFullRadialPowerLawPred(c,matsize)
% Assumes a spacing of one in each pixel direction.
% c(1) is x0, c(2) is y0, c(3) is log10(A), and c(4) is B (the exponent).
A = 10^c(3);
B = c(4);
x0 = c(1);
y0 = c(2);
xbase = 1:matsize(2);
ybase = 1:matsize(1);
[xspace,yspace] = meshgrid(xbase,ybase);
xspace_centered = xspace - x0;
yspace_centered = yspace - y0;
radii = sqrt(xspace_centered.^2 + yspace_centered.^2);
predmat = A*radii.^B;

end
