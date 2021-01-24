function [ masked_image, mask, out_lattice_points, beamstop_mask, lattice_vectors ] = makeHexagonalLatticeMask( myimage, radius, init_lattice_points, init_beamstop_mask )
% Function for making a hexagonal lattice image mask.
% Not for precise fitting of the lattice, just for processing.
%
% Nathanael Kazmierczak, Nov 2019

% Allow the user to select points via ginput until satisfied.
% 3rd argument optional.

if isa(myimage,'cell')
    beamstop_flag = 0;
    myimage = myimage{1};
else
    beamstop_flag = 1;
end

imagesize = size(myimage);
exponent = 0.3;
f1 = plotDP(myimage,exponent);
axis square;
colormap jet
colorbar


%% Obtain the locations the user wants to fit to a hexagon
if nargin <= 1
    disp('Please zoom graph in preparation to set the disk radius.');
    input('Press any key to continue.');
    disp('Please click two points to set the disk radius for the mask.');
    [xc,yc] = ginput(2);
    radius = sqrt((xc(1)-xc(2))^2 + (yc(1)-yc(2))^2);
    f1 = plotDP(myimage,exponent);
end
if nargin <= 2
    disp('Please click on the centers of disks to be cut out of the mask (data not included).');
    disp('The first clicked point will be fixed as the origin, while the second and third points define trial basis vectors.');
    disp('Thus the first three points should not be in a line.');
    disp('Beginning disk selection: click the first disk.');
    stored_locations = zeros(0,2);
    while true
        figure(f1);
        this_location = ginput(1);
        disp('Please choose an action: ');
        disp('0: Delete the current point and select again');
        disp('1: Accept the current point and select again');
        disp('2: Accept the current point and terminate point selection');
        choice = input('');
        
        if choice ~= 0
            stored_locations(end+1,:) = this_location;
        end
        
        if choice == 2
            break
        end
    end
    
    stored_locations = stored_locations';
    
    %% Mask beamstop and central area so they are not inadvertently included
    if beamstop_flag
        figure(f1)
        disp('Please draw any regions that should not be used, such as beamstop and central region.');
        beamstop_mask = logical(roipoly());
    else
        beamstop_mask = false(size(myimage));
    end
    
    %% Fit a hexagonal lattice to the selected points
    % First, index using the first point as the center and the second and third
    % points as the tips of the initial guess lattice vectors.
    v1_init = stored_locations(:,2) - stored_locations(:,1);
    v2_init = stored_locations(:,3) - stored_locations(:,1);
    v_inits = [v1_init, v2_init];
    xy0 = stored_locations(:,1);
    allpoints_centered = stored_locations - repmat(xy0,1,size(stored_locations,2));
    indices = v_inits\allpoints_centered;
    indices = round(indices);
    
    % Next, having indexed all points, figure out which lattice vectors best
    % fit those points.
    lattice_vectors = allpoints_centered/indices;
    % here v1 is first column, v2 second column, etc.
    % Note that the way this is done, we are indeed pinning the origin of the
    % lattice on the initial click. Wasn't there a way to avoid this in
    % py4DSTEM? Look into it.
    
    % Go out five lattice vectors in each direction with full permutation to
    % get the desired coverage.
    % NPK updated 04/09/2020 to expand the permutation basis, as it was not
    % fully covering some images.
    basevec = -10:10;
    combos = nchoosek(basevec,2);
    perm_idxs = vertcat(combos,horzcat(combos(:,2),combos(:,1))); % Manually enumerate permutations b/c Matlab doesn't have a good function here.
    perm_idxs_wrep = vertcat(perm_idxs,[basevec',basevec'])';
    calc_points = lattice_vectors*perm_idxs_wrep;
    calc_points = calc_points + repmat(xy0,1,size(calc_points,2));
    
    % Test and delete any points that have fallen out of bounds of the image.
    % NPK 03/15/2020: add some padding here based on the radius, because
    % even though the center may be outside of the image, there could still
    % be a contribution.
    calc_points(:,any(calc_points < -radius,1)) = [];
    calc_points(:,any(calc_points > repmat(imagesize'+radius,1,size(calc_points,2)),1)) = [];
    hold on; scatter(calc_points(1,:)',calc_points(2,:)','r','filled')
    calc_points = calc_points';
else
    calc_points = init_lattice_points;
    beamstop_mask = init_beamstop_mask;
end



% Use the constructed lattice to make the mask.
mask = false(imagesize);
rowbase = 1:imagesize(1);
colbase = 1:imagesize(2);
[rowspace,colspace] = meshgrid(rowbase,colbase);
for i = 1:size(calc_points,1)
    row0 = calc_points(i,1);
    col0 = calc_points(i,2);
    circle_mask = isInCircle(rowspace,colspace,row0,col0,radius);
    mask = mask | circle_mask;
end

% Modified by NPK on 03/15/2020 to remove this functionality, which was
% stunting the hBN mask.
% % Use the beamstop mask to filter.
% mask = mask & ~beamstop_mask;


% Plot results
factor = 10;
filtered_image = double(myimage).^exponent/factor;
filtered_image(mask) = filtered_image(mask)*factor;
% f2 = figure;
% imagesc(filtered_image);
f2 = plotDP(filtered_image,exponent);
axis square;
colormap jet
colorbar

% TODO: add an elimination feature where you can click on any disk that
% didn't come out the way out wanted and remove it from the mask.


out_lattice_points = calc_points; 
masked_image = filtered_image;

if nargin > 2
    lattice_vectors = [];
end





end

