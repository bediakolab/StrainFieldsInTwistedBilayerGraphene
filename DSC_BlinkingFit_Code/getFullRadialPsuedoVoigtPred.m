function [pvoigt_predmat] = getFullRadialPsuedoVoigtPred(c,matsize)
% Assumes a spacing of one in each pixel direction. c(1) is x0, c(2) is y0,
% c(3) is A (scaling), c(4) is fG, c(5) is fL (the full width at half max
% for Gaussian and Lorentzian, respectively).
A = c(3);
fG = c(4);
fL = c(5);
eta = c(6);
x0 = c(1);
y0 = c(2);
xbase = 1:matsize(2);
ybase = 1:matsize(1);
[xspace,yspace] = meshgrid(xbase,ybase);
xspace_centered = xspace - x0;
yspace_centered = yspace - y0;
radii = sqrt(xspace_centered.^2 + yspace_centered.^2);
sigma = fG/(2*sqrt(2*log(2)));  % Gaussian width parameter
gamma = fL/2;  % Lorentz width parameter
% f = (fG^5 + 2.69269*fG^4*fL + 2.42843*fG^3*fL^2 + 4.47163*fG^2*fL^3 + 0.07842*fG*fL^4 + fL^5)^1/5;
% eta = 1.36603*(fL/f) - 0.47719*(fL/f)^2 + 0.11116*(fL/f)^3;
lorentz_predmat = A*(1/pi)*(0.5*gamma)./(radii.^2 + (0.5*gamma)^2);
gaussian_predmat = A*(1/(sigma*sqrt(2*pi)))*exp(-0.5*(radii./sigma).^2);
pvoigt_predmat = eta*lorentz_predmat + (1-eta)*gaussian_predmat;

end

