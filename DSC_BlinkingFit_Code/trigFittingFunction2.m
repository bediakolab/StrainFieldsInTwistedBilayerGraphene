function [ pred_vals ] = trigFittingFunction2( scaling_constants, DSC_guess )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

persistent vl1 vl2 vl3 vl4 vl5 vl6 mat mul
if isempty(vl1)
    hexagon_lattice_constant = 2.461;
    [ ~, mat ] = getLatticePlaneNormalVectors(0);
    plane_sep = zeros(12,1);
    plane_sep(1:6,1) = 3*hexagon_lattice_constant/2/sqrt(3);  % inner disk lattice plane separation
    plane_sep(7:12,1) = hexagon_lattice_constant/2;  % outer disks
    mul = (pi/2)./(plane_sep/2);
    
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
end


pad = 0;
c = [];
c(:,1) = vl1 < DSC_guess - pad;
c(:,2) = vl2 < DSC_guess - pad;
c(:,3) = vl3 > DSC_guess + pad;
c(:,4) = vl4 > DSC_guess + pad;
c(:,5) = vl5 > DSC_guess + pad;
c(:,6) = vl6 < DSC_guess - pad;

%     pred_vals = 10*ones(1,12);
pred_vals(any(c,2),:) = 10000;

     
[ mags ] = project2b( DSC_guess,mat );
pred_vals = (repmat(scaling_constants,1,size(mags,2)).*cos(mags.*repmat(mul,1,size(mags,2))).^2)';




end

