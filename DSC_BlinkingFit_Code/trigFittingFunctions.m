function [ pred_vals ] = trigFittingFunctions( DSC_guess, scaling_constant, bound_handling_flag )
% Fitting functions built on the basis of electron diffraction simulations
% showing a cos^2 oscillation in the disk blinking.
%
% Default scaling constant is 0.9 when fitting data, but 1 when predicting
% blinking patterns.
%
% Nathanael Kazmierczak, Dec 2019
%
% Modification Mar 2020 -- set scaling_constant to pointwise multiplication
% so that we can have different heights for the different trig functions.

% global hexagon_lattice_constant

hexagon_lattice_constant = 2.461;
[ ~, mat ] = getLatticePlaneNormalVectors(0);
[ ~, mags ] = project2( DSC_guess,mat );

plane_sep = zeros(12,1);
plane_sep(1:6,1) = 3*hexagon_lattice_constant/2/sqrt(3);  % inner disk lattice plane separation
plane_sep(7:12,1) = hexagon_lattice_constant/2;  % outer disks

mul = (pi/2)./(plane_sep/2);

pred_vals = (scaling_constant'.*cos(mags.*mul).^2)';

% bound the acceptible region (equivalent of forbidding extrapolation in
% the scatterInterpolant class).
t = 60/180*pi;
rotmat = [cos(t) sin(t); -sin(t) cos(t)];
v1 = [0,hexagon_lattice_constant/sqrt(3)];
v2 = v1*rotmat';
vl1 = vectorLine(v2-v1,v1);
vl2 = vectorLine(v1,v2);
vl3 = vectorLine(v2,-v1);
vl4 = vectorLine(v2-v1,-v1);
vl5 = vectorLine(v1,-v2);
vl6 = vectorLine(-v2,v1);

if bound_handling_flag == 0  % the original method

c = zeros(1,6);
% if scaling_constant == 1
%     pad = 1;
% end
pad = 0;  % previously 0.1, changed 02/02/2020 to try to clean up fit.
c(1) = vl1 < DSC_guess - pad;
c(2) = vl2 < DSC_guess - pad;
c(3) = vl3 > DSC_guess + pad;
c(4) = vl4 > DSC_guess + pad;
c(5) = vl5 > DSC_guess + pad;
c(6) = vl6 < DSC_guess - pad;
if any(c)
%     pred_vals = 10*ones(1,12);
    pred_vals = nan(1,12);
end


elseif bound_handling_flag == 1
    % This is currently the method where all of the predictions will be
    % blinking on.
end


    

end

