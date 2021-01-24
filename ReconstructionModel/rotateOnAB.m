function [ xy_lat_rotated ] = rotateOnAB( xy_lat, xshift, yshift, rotation_angles )
% Equivalent function to gaussRotateOnAA.m
%
% Assumes that the rotation angles have already been built in the caller
% function. Caller function parses whether the given basis point is a
% sublattice 1 or sublattice 2.
% 
% Nathanael Kazmierczak, 06/08/2020

xlat_centered = xy_lat(:,1) - xshift;
ylat_centered = xy_lat(:,2) - yshift;
% need to get coordinates for mask, which is assumed for now to have unit
% Angstrom spacing.
[r,c] = size(rotation_angles);
xbase = -(c-1)/2:(c-1)/2;
ybase = -(r-1)/2:(r-1)/2;
[xspace,yspace] = meshgrid(xbase,ybase);
% Default extrapolation value is zero angle.
angles_for_each_atom = interp2(xspace,yspace,rotation_angles,xlat_centered,ylat_centered,'linear',0);
if nnz(angles_for_each_atom) == 0
    xy_lat_rotated = xy_lat;
    return
end
rotated_centered_lat = zeros(size(xy_lat));
for i = 1:numel(xlat_centered)
    thisx = xlat_centered(i);
    thisy = ylat_centered(i);
    this_lat_vals = [thisx;thisy];
    this_angle = deg2rad(angles_for_each_atom(i));
    rotmat = [cos(this_angle), -sin(this_angle); sin(this_angle), cos(this_angle)];
    this_rotated_lat_vals = rotmat*this_lat_vals;
    rotated_centered_lat(i,:) = this_rotated_lat_vals;
end
% Undo the shift so the lattice is back in the original location.
xy_lat_rotated = rotated_centered_lat + [xshift,yshift];

end

