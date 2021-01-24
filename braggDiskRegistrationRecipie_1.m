function [disk_locations, disk_intensities, com_coordinates] = braggDiskRegistrationRecipie_1(DP,hexagonal_window_masks,kernel_radius,sigma,min_disk_distance,num_to_find_per_window,make_figures,min_intens_rel_thres,max_intens_rel_thres,subpixel)
% Requires data to already be read in.
% A function-based version of the initial testing script,
% bragg_disk_cross_correlation.m
%
% Presumes that the kernel will be a simple disk for a circular moving
% average treatment.
%
% COM coordinates assumes the use of inner hBNs and outer graphenes, in
% which case only eight rather than all ten are used.
%
% Note that there is a lot of code which is shared between this function
% and the function braggDiskRegistrationRecipe_2.m, so much of the code has
% been refactored into braggDiskRegistrationRecipe_innards.m. (On further
% contemplation, this seemed premature, so just pass an ROI kernel as
% kernel_radius and the program will figure it out.)
%
% HAS BEEN REFACTORED; NOW USES CONV2 RATHER THAN XCORR2.

TEST_FIGURES = 0;

DP_dim = size(DP);
idx = round(DP_dim/2);
xbase = 1:DP_dim(2);
ybase = 1:DP_dim(1);
[xspace,yspace] = meshgrid(xbase,ybase);
if numel(kernel_radius) == 1
    probe_kernel = isInCircle(xspace,yspace,idx(1),idx(2),kernel_radius);
else
    probe_kernel = kernel_radius;
end

if TEST_FIGURES
    figure
    pcolor(probe_kernel); shading flat;
end

%%%%% Code after construction of the kernel (convolution and maxima
%%%%% finding) has now been refactored into the function
%%%%% braggDiskRegistrationRecipe_innards.m to accomodate registration of
%%%%% non-disk shapes.

[BW,conv_map,intensities_subpixel,coords_subpixel,intensities_pixel,...
    coords_pixel] = braggDiskRegistrationRecipe_innards(DP,probe_kernel,sigma,subpixel);

%% Now, find the two most intense peaks in each window.

% Since we are also going to apply this to hBN, we want to be able to get
% only the most intense peak as well. At some point, come up with a
% general, elegant way to handle this, but for now just hack it out.

number_of_window_masks = size(hexagonal_window_masks,3);
if num_to_find_per_window == 2
    maximastorage = zeros(2*number_of_window_masks,4);  % add a fourth column which will be the ID specifying which mask region the point came from.
elseif num_to_find_per_window == 1
    maximastorage = zeros(number_of_window_masks,4);
end
for i = 1:number_of_window_masks
    this_mask = hexagonal_window_masks(:,:,i);
    if strcmp(subpixel,'False')
        these_maxima = BW & this_mask;
        %%% Refactored 09032019: generalized maxima finding.
        [ intensities, xinds, yinds, ~ ] = convMax( these_maxima, conv_map );
    else  % yee haw
        % We need the pixel values because we need to know if the subpixel
        % falls within the mask range; do this by testing the pixel value and accepting
        % subpixel if pixel is in mask. Possibly rounding subpixel would be a
        % better/faster way of doing this so we only need one call to
        % python; do that in the future.
        lininds = sub2ind(size(this_mask),coords_pixel(:,1),coords_pixel(:,2));
        userows = this_mask(lininds);
        xinds = coords_subpixel(logical(userows),1);
        yinds = coords_subpixel(logical(userows),2);
        intensities = intensities_subpixel(logical(userows),1);
    end
    combodata = [intensities,xinds,yinds,repmat(i,size(intensities,1),1)];  % where i is now going to index which window this originates from, despite any sortrows processing.
    sortedcombodata = flipud(sortrows(combodata));
    % Before storing, need to guard against too many minima near the same
    % location:
    [sortedcombodata] = removeTooNear(sortedcombodata,min_disk_distance);
    num_found = size(sortedcombodata,1);
    
    if num_to_find_per_window == 2
        switch num_found
            case 0
                maximastorage(i*2-1:i*2,:) = [nan(2,3),[i;i]];
            case 1
                maximastorage(i*2-1:i*2,:) = vertcat(sortedcombodata(1,:),[nan(1,3),i]);
            otherwise
                maximastorage(i*2-1:i*2,:) = sortedcombodata(1:2,:);
        end
    elseif num_to_find_per_window == 1
        switch num_found
            case 0
                maximastorage(i,:) = [nan(1,3),i];
            otherwise
                maximastorage(i,:) = sortedcombodata(1,:);
        end
    end
    
end

%% Threshold peak intensities both above and below
%Added 08132019: code for ensuring that we do not get peaks that
% are too low in intensity (probably off-target) or too high
% (probably hBN). Compute relative to median on the assumption that
% generally we will get things that make sense.

med_intensity = median(maximastorage(:,1),'omitnan');
% Those thresholds should be specified as multipliers
lb = min_intens_rel_thres*med_intensity;
ub = max_intens_rel_thres*med_intensity;
for i = 1:(number_of_window_masks*num_to_find_per_window) % we always populate all of them, even if some are nans
    this_disk_intensity = maximastorage(i,1);
    if lb > this_disk_intensity
        maximastorage(i,:) = [nan(1,3),maximastorage(i,4)];
    end
    if ub < this_disk_intensity
        maximastorage(i,:) = [nan(1,3),maximastorage(i,4)];
    end
end

%% Plot on top of the maxima map from correlation and smoothing

if make_figures
    figure;
    plotmask = sum(hexagonal_window_masks,3);
    plotmask(plotmask == 0) = 0.1;
    %plotdata = plotmask .* log10(smootheddata);
    plotdata = plotmask .* smootheddata;
    gamma = 0;
    if gamma 
        plotdata_orgmax = max(max(plotdata));
        plotdata_gamma = plotdata.^gamma;
        plotdata_gamma_max = max(max(plotdata_gamma));
        plotdata_gamma_normalized = plotdata_gamma./plotdata_gamma_max.*plotdata_orgmax;
        plotdata = plotdata_gamma_normalized;
    end
    image(plotdata,'CDataMapping','scaled'); shading flat;
    ax = gca;
    hold on;
    h = scatter(ax,maximastorage(:,3),maximastorage(:,2),'filled'); colormap(ax,hsv);
end

%% Plot on top of the raw data

if make_figures
    figure;
    plotmask = sum(hexagonal_window_masks,3);
    plotmask(plotmask == 0) = 0.1;
    plotdata = plotmask .* double(uint16(DP));
    image(plotdata,'CDataMapping','scaled'); shading flat;
    ax = gca;
    hold on;
    h = scatter(ax,maximastorage(:,3),maximastorage(:,2),'filled'); colormap(ax,pink);
end
switch num_to_find_per_window
    case 1  % using inner hBNs; average all of them
        com = mean(maximastorage(1:6,2:3),1,'omitnan');
    case 2  % using outer graphenes; use only 1:4, 7:10
        com = mean(maximastorage([1:4,7:10],2:3),1,'omitnan');
end
% Testing code: plot the center of mass from spots 1,3,7,9, and see if those coincide with the center of the beam.
% Because we have pulled columns 2:3, this ordering in the scatter is
% correct.
if make_figures
    scatter(ax,com(2),com(1),'r','filled');
end
%% Return what has been found 

disk_locations = fliplr(maximastorage(:,2:4));  % in which the window it came from is the first coordinate, x is the second, y is the third
% that is, disk_locations column 3 has xcoords, while column 4 has y coords
disk_intensities = maximastorage(:,1);
com_coordinates = fliplr(com);




end

