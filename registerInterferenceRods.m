function [ sortedmaxcombo ] = registerInterferenceRods( roicut,kernel,sigma,subpixel, rel_intensity_threshold, min_disk_distance )
% Refactored out of track_intereference_motion.m

[subpixel_mask,conv_map,~,~,~,~] = braggDiskRegistrationRecipe_innards(roicut,kernel,sigma,subpixel);
% going into the thresholding, the typical format is [intensity, Icoord?,
% Jcoord?, window id]
[ intensities, xinds, yinds, maxima_storage ] = convMax( subpixel_mask, conv_map );

%% (4) Filter for reliable registrations
% New function built for thresholding 
[ sortedmaxcombo ] = relIntensityThreshold( maxima_storage, rel_intensity_threshold );
[sortedmaxcombo] = removeTooNear(sortedmaxcombo,min_disk_distance);


end


