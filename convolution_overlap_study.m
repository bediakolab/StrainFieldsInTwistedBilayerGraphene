% convolution_overlap_study.m
%
% Script for studying the convolution map of two Bragg disks as they
% approach and eventually overlap.
%
% Nathanael Kazmierczak, 08/10/2019
% Bediako Lab, UC Berkeley

%% Define a 601-301 space with two very large disks.

xbase = 1:601;
ybase = 1:301;
[xspace,yspace] = meshgrid(xbase,ybase);  % helpfully aligns x and y the way we think of them.
x0 = round(mean([xbase(1),xbase(end)]));
y0 = round(mean([ybase(1),ybase(end)]));
r = 100;  % large disk fully fitting on the image.
spacings = 400:-10:0;  % values at which we will compute the convolution map.
% Nice thing about this symmetry is that we will be able to take only the
% slice down the middle. Manually normalize the correlation map for each
% slice instead of using normxcorr because people have been fussing about
% that.
numspacings = numel(spacings); 
convolution_kernel = double(isInCircle(xspace,yspace,x0,y0,r));  % fits the shape exactly by design.

additive_convolution_matrices_storage = zeros(size(xspace,1),size(xspace,2),numspacings);
additive_convolution_slices_storage = zeros(numspacings,size(xspace,2));
additive_convolution_localmax_storage = zeros(0,4); % coordinates are spacing, x, y.
not_additive_convolution_matrices_storage = zeros(size(xspace,1),size(xspace,2),numspacings);
not_additive_convolution_slices_storage = zeros(numspacings,size(xspace,2));
not_additive_convolution_localmax_storage = zeros(0,4); % coordinates are spacing, x, y, value.

%set range for producing convolution overlayed matrix.
num_x = length(xbase);
num_y = length(ybase);
idx_x = round(num_x/2);
idx_y = round(num_y/2);
starty = idx_y + 1;
startx = idx_x + 1;
rangey = starty:((starty+num_y)-1);
rangex = startx:((startx+num_x)-1);

% verify the correct convolution kernel
figure;
pcolor(xbase,ybase,convolution_kernel); shading flat;
title 'Convolution kernel'
xlabel('x')
ylabel('y');

indices_of_interest = [1,25,36];

for i = 1:numspacings
    i
    % Generate circular masks
    this_spacing = spacings(i);
    xcenter1 = x0 + this_spacing/2;  % ycenters are fixed (first matrix dimension);
    xcenter2 = x0 - this_spacing/2;
    
    circlemask1 = double(isInCircle(xspace,yspace,xcenter1,y0,r));
    circlemask2 = double(isInCircle(xspace,yspace,xcenter2,y0,r));
    
    % Now explore two cases: intensities are additive, and intensities are not additive.
    % I think the former is probably the case.
    additive_mask = circlemask1 + circlemask2;
    not_additive_mask = double(circlemask1 | circlemask2);
    
    % verify at key points that the masks are being generated according to
    % what we think.
    if ismember(i, indices_of_interest)
        figure; subplot(2,1,1);
        pcolor(xbase,ybase,additive_mask); shading flat; colorbar;
        titlestring = sprintf('Additive mask: spacing = %d',this_spacing);
        title(titlestring);
        xlabel('x');
        ylabel('y');
        
        subplot(2,1,2);
        pcolor(xbase,ybase,not_additive_mask); shading flat; colorbar;
        titlestring = sprintf('Non-additive mask: spacing = %d',this_spacing);
        title(titlestring);
        xlabel('x');
        ylabel('y');
    end
    
    % convolve (xcorr just flips a matrix upside down as compared to conv2,
    % it seems, since our matrices are entirely real).
    additive_convolution = xcorr2(additive_mask,convolution_kernel);
    not_additive_convolution = xcorr2(not_additive_mask,convolution_kernel);
    
    % turn these into overlays with the correct dimension
    additive_convolution_overlay = additive_convolution(rangey,rangex);
    not_additive_convolution_overlay = not_additive_convolution(rangey,rangex);
    
    % ascertain all the local maxima falling within a factor of 0.8 of the max
    % intensity value (so that we don't get baseline noise).
    threshold = 0.8;
    for j = 1:2  % for both additive and non-additive
        switch j
            case 1
                thismat = additive_convolution_overlay;
            case 2
                thismat = not_additive_convolution_overlay;
        end
        allmax = find(imregionalmax(thismat));
        thresholded_max = allmax;
        thresholded_max(thismat(allmax) < threshold) = []; % check and make sure this actually works when we start introducing error.
        values = thismat(thresholded_max);
        [yinds, xinds] = ind2sub(size(thismat),thresholded_max);  % thresholded_max is full of linear indices
        for k = 1:length(thresholded_max)  % because we don't know how many there are, store one at a time.
            switch j
                case 1
                    additive_convolution_localmax_storage(end+1,:) = [this_spacing,xinds(k),yinds(k),values(k)];
                case 2
                    not_additive_convolution_localmax_storage(end+1,:) = [this_spacing,xinds(k),yinds(k),values(k)];
            end
        end
    end
    
    % store the matrices and slices
    additive_convolution_matrices_storage(:,:,i) = additive_convolution_overlay;
    not_additive_convolution_matrices_storage(:,:,i) = not_additive_convolution_overlay;
    aconvslice = additive_convolution_overlay(y0,:)./(max(additive_convolution_overlay(y0,:))); % this will be normalized.
    naconvslice = not_additive_convolution_overlay(y0,:)./(max(not_additive_convolution_overlay(y0,:))); % this will be normalized.
    additive_convolution_slices_storage(i,:) = aconvslice;
    not_additive_convolution_slices_storage(i,:) = naconvslice;
end

% graphs of the slices
figure; 
contourf_thresholded_slicemax(xbase,spacings,additive_convolution_slices_storage,threshold,25); 
shading flat;
title 'Additive bragg disks: convolution map'
xlabel 'X coordinate on image'
ylabel 'Spacing between bragg disk centers'
colormap(pink);
colorbar;

figure; 
contourf_thresholded_slicemax(xbase,spacings,not_additive_convolution_slices_storage,threshold,25); 
title 'Non-additive bragg disks: convolution map'
xlabel 'X coordinate on image'
ylabel 'Spacing between bragg disk centers'
colormap(pink);
colorbar;
hold on


% directly examine the convolution map for the three cases of interest.
for j = 1:numel(indices_of_interest)
    this_index = indices_of_interest(j);
    this_spacing = spacings(this_index);
    figure; subplot(2,1,1);
    contourf(xbase,ybase,additive_convolution_matrices_storage(:,:,this_index),25);
    shading flat;
    titlestring = sprintf('Additive convolution map: spacing = %d',this_spacing);
    title(titlestring);
    xlabel('x');
    ylabel('y');
    hold on 
    % find all thresholded local max belonging to this convolution map as
    % defined by its spacing
    these_max = additive_convolution_localmax_storage;
    these_max(~(these_max(:,1) == this_spacing),:) = [];
%     scatter3(these_max(:,2),these_max(:,3),these_max(:,4),'r','filled');
    scatter(these_max(:,2),these_max(:,3),'r','filled');    

    subplot(2,1,2);
    contourf(xbase,ybase,not_additive_convolution_matrices_storage(:,:,this_index),25);
    shading flat;
    titlestring = sprintf('Not additive convolution map: spacing = %d',this_spacing);
    title(titlestring);
    xlabel('x');
    ylabel('y');
    hold on
    % find all thresholded local max belonging to this convolution map as
    % defined by its spacing
    these_max = not_additive_convolution_localmax_storage;
    these_max(~(these_max(:,1) == this_spacing),:) = [];
    scatter(these_max(:,2),these_max(:,3),'r','filled');
end    
