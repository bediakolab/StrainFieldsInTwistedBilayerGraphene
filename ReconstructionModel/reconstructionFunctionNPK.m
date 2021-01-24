function [ xy_lat, rotfield ] = reconstructionFunctionNPK( AA_angle, AA_distance, AB_angle, AB_buffer_distance, AB_smooth_distance, xy_lat, cellDimXY, plotrotfield )
% Simplified reconstruction model function that takes into account not just
% AA rotation (as in Tadmor, see gaussRotateLattice03.m), but also the
% fixed-body rotation in AB domains.
%
% Nathanael Kazmierczak, 06/08/2020

TEST = false;

%% First, build the rotation field.
perm_val = 9;
basis = permn(-perm_val:perm_val,2);
% [ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
w1 = [cellDimXY(1),0];
w2 = cellDimXY./2;
W = [w1',w2'];
genpoints = (W*basis')';
AB_lat1 = genpoints + [cellDimXY(1)/3,0];
AB_lat2 = genpoints + [cellDimXY(1)/6,cellDimXY(2)/2];

% Trim anything more than n unit cells away
unit_cell_max = 1;
to_delete = genpoints(:,1) > cellDimXY(1)*(unit_cell_max+1) | genpoints(:,1) < -cellDimXY(1)*unit_cell_max | genpoints(:,2) > cellDimXY(2)*(unit_cell_max+1) | genpoints(:,2) < -cellDimXY(2)*unit_cell_max;
to_delete_ABlat1 = AB_lat1(:,1) > cellDimXY(1)*(unit_cell_max+1) | AB_lat1(:,1) < -cellDimXY(1)*unit_cell_max | AB_lat1(:,2) > cellDimXY(2)*(unit_cell_max+1) | AB_lat1(:,2) < -cellDimXY(2)*unit_cell_max;
to_delete_ABlat2 = AB_lat2(:,1) > cellDimXY(1)*(unit_cell_max+1) | AB_lat2(:,1) < -cellDimXY(1)*unit_cell_max | AB_lat2(:,2) > cellDimXY(2)*(unit_cell_max+1) | AB_lat2(:,2) < -cellDimXY(2)*unit_cell_max;
genpoints(to_delete,:) = [];
AB_lat1(to_delete_ABlat1,:) = [];
AB_lat2(to_delete_ABlat2,:) = [];
f1 = figure; scatter(genpoints(:,1),genpoints(:,2));

% These are the AA rotation points. Center a Gaussian rotation field around
% each of them.
n_AA_rotations = size(genpoints,1);
for i = 1:n_AA_rotations
    x_shift = genpoints(i,1);
    y_shift = genpoints(i,2);
    [ xy_lat ] = gaussRotateOnAA( AA_angle, AA_distance, xy_lat, x_shift, y_shift );
end

if TEST
    figure; scatter(xy_lat(:,1),xy_lat(:,2),'filled');
end

% Rotate around AB rotation points in much the same way.

AB_rotation_points = vertcat(AB_lat1,AB_lat2);
if TEST
    figure(f1); hold on;
    scatter(AB_rotation_points(:,1),AB_rotation_points(:,2)); legend('AA','AB'); axis equal; xlabel('x'); ylabel('y'); set(gca,'yDir','normal');
end
% Method: translate to the origin for each, and then apply appropriate
% rotation field. The field is different for the two sublattices, because of the basis.

xbase = -cellDimXY(1):cellDimXY(1);
ybase = -cellDimXY(2):cellDimXY(2);
[xspace,yspace] = meshgrid(xbase,ybase);
vertices = [AB_buffer_distance - 1/6*cellDimXY(1), -1/2*cellDimXY(2) + sqrt(3)*AB_buffer_distance;
    AB_buffer_distance - 1/6*cellDimXY(1), 1/2*cellDimXY(2) - sqrt(3)*AB_buffer_distance;
    -sqrt(3)*AB_buffer_distance + 1/3*cellDimXY(1), 0];
% for testing
nobuf_vertices = [-1/6*cellDimXY(1), -1/2*cellDimXY(2);
    -1/6*cellDimXY(1), 1/2*cellDimXY(2);
    1/3*cellDimXY(1), 0];
% This will make them line up with the sublattices correctly. 
[ mask1 ] = isInTriangle( xspace, yspace, -vertices );
[ mask2 ] = isInTriangle( xspace, yspace, vertices );
if TEST
    [ mask1_nobuf ] = isInTriangle( xspace, yspace, nobuf_vertices );
    [ mask2_nobuf ] = isInTriangle( xspace, yspace, -nobuf_vertices );
    bufmask1 = double(mask1) + double(mask1_nobuf);
    bufmask2 = double(mask2) + double(mask2_nobuf);
    figure;
    imagesc(bufmask1); axis equal; set(gca,'yDir','normal'); title('mask1'); colorbar;
    figure;
    imagesc(bufmask2); axis equal; set(gca,'yDir','normal'); title('mask2'); colorbar;
end


anglemat_1 = double(mask1)*AB_angle;
anglemat_2 = double(mask2)*AB_angle;
anglemat_1 = imgaussfilt(anglemat_1,AB_smooth_distance);
anglemat_2 = imgaussfilt(anglemat_2,AB_smooth_distance);

if TEST
    figure;
    imagesc(anglemat_1); axis equal; set(gca,'yDir','normal'); title('Angle matrix 1'); colorbar;
    figure;
    imagesc(anglemat_2); axis equal; set(gca,'yDir','normal'); title('Angle matrix 2'); colorbar;
end

% Perform AB rotations
nLat1 = size(AB_lat1,1);
nLat2 = size(AB_lat2,1);
% Need to make sure this is the correct rotation matrix here.
for i = 1:nLat1
    if mod(i,10) == 0
        fprintf('Rotating around AB domain %d of %d, sublattice 1\n',i,nLat1);
    end
    xshift = AB_lat1(i,1);
    yshift = AB_lat1(i,2);
    [ xy_lat ] = rotateOnAB( xy_lat, xshift, yshift, anglemat_1 );
end
for i = 1:nLat2
    if mod(i,10) == 0
        fprintf('Rotating around AB domain %d of %d, sublattice 2\n',i,nLat2);
    end
    xshift = AB_lat2(i,1);
    yshift = AB_lat2(i,2);
    [ xy_lat ] = rotateOnAB( xy_lat, xshift, yshift, anglemat_2 );
end

%% Visualize the rotation field that has been imparted on the data
if nargout > 2 || plotrotfield
rotxbase = -unit_cell_max*cellDimXY(1):(unit_cell_max+1)*cellDimXY(1);
rotybase = -unit_cell_max*cellDimXY(2):(unit_cell_max+1)*cellDimXY(2);
[rotxspace,rotyspace] = meshgrid(rotxbase,rotybase);
rotfield = zeros(size(rotxspace));
for i = 1:n_AA_rotations
    x_shift = genpoints(i,1);
    y_shift = genpoints(i,2);
    dists = sqrt((rotxspace-x_shift).^2 + (rotyspace-y_shift).^2);
    thisAArot = AA_angle.*exp(-(dists.^2)./2.0./(AA_distance.^2));  % compute this in degrees
    rotfield = rotfield + thisAArot;
end
for i = 1:nLat1
    xshift = AB_lat1(i,1);
    yshift = AB_lat1(i,2);
    these_vertices = -vertices + [xshift,yshift];
    [ mask1 ] = isInTriangle( rotxspace, rotyspace, these_vertices );
    this_anglemat_1 = double(mask1)*AB_angle;
    this_anglemat_1 = imgaussfilt(this_anglemat_1,AB_smooth_distance);
    rotfield = rotfield + this_anglemat_1;
end
for i = 1:nLat2
    xshift = AB_lat2(i,1);
    yshift = AB_lat2(i,2);
    these_vertices = vertices + [xshift,yshift];
    [ mask2 ] = isInTriangle( rotxspace, rotyspace, these_vertices );
    this_anglemat_2 = double(mask2)*AB_angle;
    this_anglemat_2 = imgaussfilt(this_anglemat_2,AB_smooth_distance);
    rotfield = rotfield + this_anglemat_2;
end

if TEST || plotrotfield
    figure;
    imagesc(rotxbase,rotybase,rotfield); set(gca,'yDir','normal'); axis equal; colorbar;
    hold on;
    scatter(genpoints(:,1),genpoints(:,2),'k');
    scatter(AB_rotation_points(:,1),AB_rotation_points(:,2),'r'); legend('AA','AB'); xlabel('x'); ylabel('y');
    xlim([min(rotxbase)-5,max(rotxbase)+5]);
    ylim([min(rotybase)-5,max(rotybase)+5]);
end
end

end

