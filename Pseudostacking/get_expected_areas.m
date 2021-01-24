% get_expected_areas.m
%
% This script uses math for the Wigner-Seitz modified partition of the
% half-hexagon for pseudostacking. It covers the unreconstructed case where
% the moire evenly samples the fitting region.
%
% Nathanael Kazmierczak, 07/11/2020

[ ~, ~, a ] = getDSCBasisVectors();
abond = a/sqrt(3);
area_aa = pi/8*(abond)^2;
area_sp = 6*(3/8*tan(deg2rad(15))*abond^2 - pi/84*abond^2);
area_ab = 6*(sqrt(3)/8*abond^2*(1-sqrt(3)*tan(deg2rad(15))) - pi/84*abond^2);
area_total = area_aa + area_sp + area_ab;
AApercent = area_aa/area_total*100
ABpercent = area_ab/area_total*100
SPpercent = area_sp/area_total*100
