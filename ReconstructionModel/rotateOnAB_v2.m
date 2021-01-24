function [ xy_lat_recon ] = rotateOnAB_v2( xy_lat, xshift, yshift, rotation_angles, gamma, is_right_facing, corner_angle_deg, side_length )
% Equivalent function to gaussRotateOnAA.m
%
% Assumes that the rotation angles have already been built in the caller
% function. Caller function parses whether the given basis point is a
% sublattice 1 or sublattice 2.
%
% Modified 07/27/2020 to become rotateOnAB_v2.m. Is_up_facing is now going
% to determine whether the given triangle points up or down relative to y,
% so that we can compute the direction parallel to the boundary for each
% atom. 
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
    xy_lat_recon = xy_lat;
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
xy_lat_rotated_central = rotated_centered_lat + [xshift,yshift];


% Obtain the displacement change magnitude for each atom.
diff_circular = xy_lat_rotated_central - xy_lat;
diffmag = sum(diff_circular.^2,2).^0.5;
tol = 1e-6;
has_changed = diffmag > tol;
% Compute the direction of displacement change relative to the triangle.
centroid_angles = mod(atan2(ylat_centered,xlat_centered),2*pi);
corner_angle_rad = deg2rad(corner_angle_deg);
% Assign the type of angle 
if is_right_facing
    side_NEE_corner = false(size(diffmag));
    side_SEE_corner = false(size(diffmag));
    side_NNE_corner = false(size(diffmag));
    side_SSE_corner = false(size(diffmag));
    side_NNW_corner = false(size(diffmag));
    side_SSW_corner = false(size(diffmag));
    side_NE = false(size(diffmag));
    side_SE = false(size(diffmag));
    side_W = false(size(diffmag));
    NEE_angle = 0 + corner_angle_rad;
    SEE_angle = 2*pi - corner_angle_rad;
    NNE_angle = 2*pi/3 - corner_angle_rad;
    NNW_angle = 2*pi/3 + corner_angle_rad;
    SSW_angle = 4*pi/3 - corner_angle_rad;
    SSE_angle = 4*pi/3 + corner_angle_rad;
    side_NE(has_changed) = centroid_angles(has_changed) >= NEE_angle & centroid_angles(has_changed) < NNE_angle;
    side_W(has_changed) = centroid_angles(has_changed) >= NNW_angle & centroid_angles(has_changed) < SSW_angle;
    side_SE(has_changed) = centroid_angles(has_changed) >= SSE_angle & centroid_angles(has_changed) < SEE_angle;
    side_NEE_corner(has_changed) = centroid_angles(has_changed) >= 0 & centroid_angles(has_changed) < NEE_angle;
    side_NNE_corner(has_changed) = centroid_angles(has_changed) >= NNE_angle & centroid_angles(has_changed) < 2*pi/3;
    side_NNW_corner(has_changed) = centroid_angles(has_changed) >= 2*pi/3 & centroid_angles(has_changed) < NNW_angle;
    side_SSW_corner(has_changed) = centroid_angles(has_changed) >= SSW_angle & centroid_angles(has_changed) < 4*pi/3;
    side_SSE_corner(has_changed) = centroid_angles(has_changed) >= 4*pi/3 & centroid_angles(has_changed) < SSE_angle;
    side_SEE_corner(has_changed) = centroid_angles(has_changed) >= SEE_angle & centroid_angles(has_changed) < 2*pi;
    
%     NE_direction = -[cos(5*pi/6),sin(5*pi/6)];
%     SE_direction = [cos(7*pi/6),sin(7*pi/6)];
%     W_direction = [0,1];
    NE_angle = 11*pi/6;
    SE_angle = 7*pi/6;
    W_angle = 3*pi/6;


%     side_directions = zeros(size(diffmag,1),2);
%     side_directions(side_NE,:) = NE_direction.*ones(nnz(side_NE),2); % This assignment may or may not work
%     side_directions(side_SE,:) = SE_direction.*ones(nnz(side_SE),2); % This assignment may or may not work
%     side_directions(side_W,:) = W_direction.*ones(nnz(side_W),2); % This assignment may or may not work
    side_angles = zeros(size(diffmag,1),1);
    side_angles(side_NE,:) = NE_angle;
    side_angles(side_SE,:) = SE_angle;
    side_angles(side_W,:) = W_angle;

    
%     SEdir = SE_direction.*ones(nnz(side_NEE_corner | side_SEE_corner),2);
%     NEdir = NE_direction.*ones(nnz(side_NEE_corner | side_SEE_corner),2);
    ca_modpi = centroid_angles;
    ca_modpi(ca_modpi > pi) = ca_modpi(ca_modpi > pi) - 2*pi;
    a1 = SEE_angle - 2*pi;
    a2 = NEE_angle;
    t = (ca_modpi(side_NEE_corner | side_SEE_corner)-a1)./(a2-a1);
%     side_directions(side_NEE_corner | side_SEE_corner,:) = t.*NEdir + (1-t).*SEdir;
    side_angles(side_NEE_corner | side_SEE_corner,:) = t.*(NE_angle) + (1-t).*SE_angle;
    
    comboind = side_NNE_corner | side_NNW_corner;
%     NEdir = NE_direction.*ones(nnz(comboind),2);
%     Wdir = W_direction.*ones(nnz(comboind),2);
    a1 = NNE_angle;
    a2 = NNW_angle;
    t = (centroid_angles(comboind)-a1)./(a2-a1);
%     side_angles(comboind,:) = t.*Wdir + (1-t).*NEdir;
    side_angles(comboind,:) = mod(t.*W_angle + (1-t).*(NE_angle-2*pi),2*pi);
    
    comboind = side_SSW_corner | side_SSE_corner;
%     SEdir = SE_direction.*ones(nnz(comboind),2);
%     Wdir = W_direction.*ones(nnz(comboind),2);
    a1 = SSW_angle;
    a2 = SSE_angle;
    t = (centroid_angles(comboind)-a1)./(a2-a1);
    side_angles(comboind,:) = t.*SE_angle + (1-t).*W_angle;
    
    % Flip side angles if this is a positive rotation
    if max(rotation_angles(:)) > min(rotation_angles(:))
        side_angles = mod(side_angles + pi,2*pi);
    end
    
    % Now, must determine how far away we are from the triangle edge. This
    % will determine the degree to which we interpolate the angle. 
    ca1 = centroid_angles - pi/3;
    ca2 = centroid_angles - pi;
    ca3 = centroid_angles - 5*pi/3;
    effective_angle = min(abs([ca1,ca2,ca3]),[],2);
    d = side_length./(2*sqrt(3)*cos(effective_angle));
    centroid_dists = sum(xlat_centered.^2 + ylat_centered.^2,2).^0.5;
    dfrac = zeros(size(diffmag));
    dfrac(has_changed) = centroid_dists(has_changed)./d(has_changed);
    dfrac(dfrac > 1) = 1;
    
    % Get the angles
    circular_diff_angle = atan2(diff_circular(:,2),diff_circular(:,1));
    circular_diff_angle_p2pi = circular_diff_angle + 2*pi;
    dist_orgphase = abs(circular_diff_angle - side_angles);
    dist_newphase = abs(circular_diff_angle_p2pi - side_angles);
    use_org = dist_orgphase <= dist_newphase;
    use_new = dist_orgphase > dist_newphase;
    % Obtain interpolation variable (can put other functions here);
    td = dfrac.^gamma;
    newangles = zeros(size(diffmag));
    newangles(use_org) = (1-td(use_org)).*circular_diff_angle(use_org) + td(use_org).*side_angles(use_org);
    newangles(use_new) = (1-td(use_new)).*circular_diff_angle_p2pi(use_new) + td(use_new).*side_angles(use_new);
    % Must find new displacement
    newdisp = zeros(size(diffmag,1),2);
    newdisp(has_changed,:) = diffmag(has_changed).*[cos(newangles(has_changed)),sin(newangles(has_changed))];
    
    xy_lat_recon = xy_lat + newdisp;
    return
end


% Form output
% temp
xy_lat_recon = xy_lat_rotated_central;
end

