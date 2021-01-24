function [predmat] = getFullRadialPowerLawPred_bootstrap(c,linear_data)
% c(1) is x0, c(2) is y0, c(3) is log10(A), and c(4) is B (the exponent).
%
% Difference here from getFullRadialPowerLawPred is that we are operating
% on a linear input of data. Column one will be the intensity. Column 2 is
% ycoord on the meshgrid (array dimension 1), while column 3 is xcoord.
%
% Radii are still computed in here, since they change on each iteration.
A = 10^c(3);
B = c(4);
x0 = c(1);
y0 = c(2);

x_centered = linear_data(:,3) - x0;
y_centered = linear_data(:,2) - y0;
radii = sqrt(x_centered.^2 + y_centered.^2);
predmat = A*radii.^B;

end

