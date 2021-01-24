function [ average_probe_convolution,average_probe_stacking ] = getAveragedDiskShape_innards( DPROIs,subpixel )
% Refactored code from getAveragedDiskShape.m for group use there and in 
% getAveragedDiskShape_noIndices.m

average_probe_stacking = mean(DPROIs,3,'omitnan');
% All of the DPROIs are now the same size, but they are not yet aligned. To
% do this, we will successively convolve them against each other. At some
% point, make a general correlation function for Matlab that can handle the
% phase correlation, the hybrid correlation, etc.
aligned_DPROIs = zeros(size(DPROIs));
aligned_DPROIs(:,:,1) = DPROIs(:,:,1);  % start aligning to the first DP
num_DPs = size(DPROIs,3);
for i = 2:num_DPs
    ref_DPROI = mean(aligned_DPROIs(:,:,1:(i-1)),3);
    this_DPROI = DPROIs(:,:,i);
    % correlate to figure out the best alignment. Don't worry about
    % subpixel for now, as the averaging should help out.
    convmap = conv2(double(this_DPROI),double(ref_DPROI),'same');
    switch subpixel
        case 'False'
            these_maxima = imregionalmax(convmap);
            lininds_maxima = find(these_maxima);
            [xinds,yinds] = ind2sub(size(these_maxima),lininds_maxima);
            intensities = convmap(lininds_maxima);
        case 'True'  % This is going to be quite hard.
            error('not yet supported');
    end
    [~,maxind] = max(intensities);
    xind = xinds(maxind);
    yind = yinds(maxind);
    % Need to worry about odd and even
    if mod(size(convmap,2),2)
        refx = (size(convmap,2)+1)/2 - 1;
    else
        refx = size(convmap,2)/2 - 1;
    end
    if mod(size(convmap,1),2)
        refy = (size(convmap,1)+1)/2 - 1;
    else
        refy = size(convmap,1)/2 - 1;
    end
    shiftx = xind - refx;
    shifty = yind - refy;
    % Function I write to do the index shifting.
    shiftedDPROI = shiftImage(this_DPROI,shiftx,shifty);
    aligned_DPROIs(:,:,i) = shiftedDPROI;
    disp(i);
end


% All DPROIs have been aligned -- we will now average through them to get
% averaged probe.
average_probe_convolution = mean(aligned_DPROIs,3,'omitnan');
figure; image(average_probe_convolution,'CDataMapping','Scaled'); colormap(pink); colorbar;
pbaspect([1 1 1]);
title('Averaged Probe from Convolution Alignment');
figure; image(average_probe_stacking,'CDataMapping','Scaled'); colormap(pink); colorbar;
pbaspect([1 1 1]);
title('Averaged Probe from Stacking Registration Locations');
% Plot four other examples for giggles.
fourrands = randi(num_DPs,4,1);
figure; 
for i = 1:4
    subplot(2,2,i);
    image(aligned_DPROIs(:,:,i),'CDataMapping','Scaled'); colormap(pink);
    pbaspect([1 1 1]);
    title(strcat(['Aligned Disk #',num2str(fourrands(i))]));
end



end

