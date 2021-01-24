function [ SP_line ] = registerUniqueSPLines( displacement_field, AA_mask, SP_id, plot_flag, overlay_image, upsample_flag, optional_SP_tolerance, tear_mask )
% Nathanael Kazmierczak, 04/02/2020
%
% Added 06/01/2020: section of code for protecting the edges from dissolve.

trimsize = 0; % Gives a one pixel wide boundary mask
boundary_mask_1 = logical(trimArray(zeros(size(AA_mask)),trimsize,ones(size(AA_mask))));
trimsize = 1; % Gives a one pixel wide boundary mask
boundary_mask_2 = logical(trimArray(zeros(size(AA_mask)),trimsize,ones(size(AA_mask))));
trimsize = 4; % Gives a five pixel wide boundary mask
boundary_mask_5 = logical(trimArray(zeros(size(AA_mask)),trimsize,ones(size(AA_mask))));
useNew = true;

if ~isempty(tear_mask)
    tear_boundary_mask = boundarymask(tear_mask);
    tear_boundary_mask_5 = tear_boundary_mask;
    for i = 1:4
        tear_boundary_mask_5 = boundarymask(tear_boundary_mask_5) | tear_boundary_mask_5;
        tear_boundary_mask_5(tear_mask) = false;
    end
end

SP_hsv = [0,1,1;
          0.33,1,1;
          0.66,1,1];
AB_hsv = [0,0,1;
          0.33,0,1;
          0.66,0,1];
[ ~, HSV_color_stack ] = getCustomDisplacementColor( displacement_field, SP_hsv, AB_hsv, 1, 1 );
H = HSV_color_stack(:,:,1);
S = HSV_color_stack(:,:,2);
V = HSV_color_stack(:,:,3);
if isempty(optional_SP_tolerance)
    mask1 = S > 0.45;
else
    mask1 = S > optional_SP_tolerance;
end
% NPK testing 04/25/2020
% mask1 = true(size(mask1));
switch SP_id
    case 1
        SPmask = (H <= 1/6 | H > 5/6) & mask1;
    case 2
        SPmask = (H <= 3/6 & H > 1/6) & mask1;
    case 3
        SPmask = (H <= 5/6 & H > 3/6) & mask1;
end
if ~isempty(tear_mask)
    SPmask(tear_mask) = false;
end
SPmorphfilt = bwmorph(SPmask,'clean');
SPmorphfilt = ~bwareaopen(~SPmorphfilt,50,4);  % close in dark spots
SPmorphfilt = bwareaopen(SPmorphfilt,20);  % clear out light spots
SPmorphfilt(AA_mask) = SPmask(AA_mask);  % Protect AA regions from the filter
% 06/01/2020
if useNew
    SPmorphfilt(boundary_mask_5) = SPmorphfilt(boundary_mask_5) | SPmask(boundary_mask_5);
    if ~isempty(tear_mask)
        SPmorphfilt(tear_boundary_mask_5) = SPmorphfilt(tear_boundary_mask_5) | SPmask(tear_boundary_mask_5);
    end
end
inclusive_thresh = 5;
n = 1;

% Fatten within the AA mask
% niter = 10;
% for i = 1:niter
%     filteres = ~bwneighborerode( ~SPmorphfilt, 7, 1, false );
%     SPmorphfilt(AA_mask) = filteres(AA_mask);
% end

% NPK 05/23/2020: attempted to fix the saddle point registration problems
% but deterministically populating the mask. Maybe this is actually an
% asymmetry problem on the erosion.
SPmorphfilt(AA_mask) = true;


flag = false;
if ~useNew  % this seems to have been doing more harm than good.
    SPmorphfilt = bwmorph(SPmorphfilt,'bridge');
    SPmorphfilt = ~bwneighborerode( ~SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = bwneighborerode( SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = ~bwneighborerode( ~SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = bwneighborerode( SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = ~bwneighborerode( ~SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = bwneighborerode( SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = ~bwneighborerode( ~SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = bwneighborerode( SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = ~bwneighborerode( ~SPmorphfilt, inclusive_thresh, n, flag );
else  % use bridge to protect these points
%     testmat = bwmorph(SPmorphfilt,'bridge');
% %     testmat = bwneighborerode( SPmorphfilt, inclusive_thresh, 2, flag );
% % %     bridge_pixels = SPmorphfilt & ~testmat;
% %     bridge_pixels = SPmorphfilt & ~testmat;
% % %     SPmorphfilt = ~bwneighborerode( ~SPmorphfilt, inclusive_thresh, 7, flag );
% %     thickened_bridges = ~bwneighborerode( ~bridge_pixels, inclusive_thresh, 8, flag );

    SPmorphfilt = ~bwneighborerode( ~SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = bwneighborerode( SPmorphfilt, inclusive_thresh, n, flag );
%     SPmorphfilt(bridge_pixels) = true;
    SPmorphfilt = ~bwneighborerode( ~SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = bwneighborerode( SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = ~bwneighborerode( ~SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = bwneighborerode( SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = ~bwneighborerode( ~SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = bwneighborerode( SPmorphfilt, inclusive_thresh, n, flag );
    SPmorphfilt = ~bwneighborerode( ~SPmorphfilt, inclusive_thresh, n, flag );
end
% pre_erosion_SPmorphfilt = SPmorphfilt;

% protect
% actually this doesn't seem needed
% if useNew
%     SPmorphfilt(boundary_mask_5) = SPmorphfilt(boundary_mask_5) | pre_erosion_SPmorphfilt(boundary_mask_5);
%     SPmorphfilt(tear_boundary_mask_5) = SPmorphfilt(tear_boundary_mask_5) | pre_erosion_SPmorphfilt(tear_boundary_mask_5);
% end


if upsample_flag
    % Upsample for smoothness before final pass
    upsample_factor = 4;
    SPmorphfilt = imresize(SPmorphfilt,upsample_factor);
end

pre_open_SPmorphfilt = SPmorphfilt;
SPmorphfilt = ~bwareaopen(~SPmorphfilt,50);  % close in dark spots
SPmorphfilt = bwareaopen(SPmorphfilt,20);  % clear out light spots
if useNew
    SPmorphfilt(boundary_mask_5) = SPmorphfilt(boundary_mask_5) | pre_open_SPmorphfilt(boundary_mask_5);
    if ~isempty(tear_mask)
        SPmorphfilt(tear_boundary_mask_5) = SPmorphfilt(tear_boundary_mask_5) | pre_open_SPmorphfilt(tear_boundary_mask_5);
    end
    SPmorphfilt = ~bwareaopen(~SPmorphfilt,10);  % close in dark spots
    SPmorphfilt = bwareaopen(SPmorphfilt,10);  % clear out light spots
end

% Erode until connectivity starts to change. Evaluate connectivity changes
% off the boundary
boundary_trim_value = 3;
if useNew
    edge_array = false(size(SPmorphfilt));
    edge_array(boundary_mask_1) = SPmorphfilt(boundary_mask_1);
    edge_points = bwmorph(edge_array,'shrink',inf);
    SPmorphfilt(boundary_mask_1) = edge_points(boundary_mask_1);
    if boundary_trim_value >= 0 
        conres_init = bwconncomp(trimArray(SPmorphfilt,boundary_trim_value));
    else
        conres_init = bwconncomp(SPmorphfilt);
    end
    init_conn_comp = conres_init.NumObjects;
else
    conres_init = bwconncomp(SPmorphfilt);
    init_conn_comp = conres_init.NumObjects;
end
count = 0;
pre_erosion_SPmorphfilt = SPmorphfilt;

while true
    count = count + 1;
    % NPK fixed this on 05/23/2020 to not allow breakages in the main
    % links.
    SPmorphfilt_stor = SPmorphfilt;
    SPmorphfilt = bwneighborerode( SPmorphfilt, inclusive_thresh, 1, false );
    
    if useNew
        SPmorphfilt(boundary_mask_1) = edge_points(boundary_mask_1);
        %         SPmorphfilt(boundary_mask_1) = SPmorphfilt(boundary_mask_1) | pre_erosion_SPmorphfilt(boundary_mask_1);
        if ~isempty(tear_mask)
            SPmorphfilt(tear_boundary_mask_5) = SPmorphfilt(tear_boundary_mask_5) | pre_erosion_SPmorphfilt(tear_boundary_mask_5);
        end
    end
    
    if useNew && boundary_trim_value >= 0
        conres = bwconncomp(trimArray(SPmorphfilt,boundary_trim_value));
    else
        conres = bwconncomp(SPmorphfilt);
    end
    this_conn_comp = conres.NumObjects;
    if this_conn_comp ~= init_conn_comp
        SPmorphfilt = SPmorphfilt_stor;  % Undo the result of the last filter operation.
        break
    end
    
end
% SPmorphfilt = bwneighborerode( SPmorphfilt, inclusive_thresh, 10, true );
if useNew
%     edge_select = logical(trimArray(zeros(size(SPmorphfilt)),0,ones(size(SPmorphfilt))));
    edge_array = false(size(SPmorphfilt));
    edge_array(boundary_mask_1) = SPmorphfilt(boundary_mask_1);
    edge_points = bwmorph(edge_array,'shrink',inf);
    SPmorphfilt(boundary_mask_1) = edge_points(boundary_mask_1);
    
    counter = 1;
    while true
        lastfilt = SPmorphfilt;
        SPmorphfilt = bwmorph( SPmorphfilt, 'thin', 1 );
        SPmorphfilt(boundary_mask_1) = edge_points(boundary_mask_1);
        if nnz(~(lastfilt == SPmorphfilt)) == 0 || counter > 100
            break
        end
        counter = counter + 1;
    end
else
    SPmorphfilt = bwmorph( SPmorphfilt, 'thin', inf );
end


if plot_flag
    figure; imagesc(SPmorphfilt);
    figure(overlay_image);
    hold on
    h = imagesc(cat(3,SPmorphfilt,SPmorphfilt,zeros(size(SPmorphfilt))));
    h.AlphaData = 0.25;
end

SP_line = SPmorphfilt;

end

