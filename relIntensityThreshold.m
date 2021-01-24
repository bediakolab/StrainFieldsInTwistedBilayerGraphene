function [ sortedmax ] = relIntensityThreshold( maxima_storage, rel_intensity_threshold )
% Threshold by intensity
%
% maxima_storage has the format [intensities, xcoords, ycoords]. Not
% entirely sure if those are legit xcoords or maybe Icoords.

sortedmax = sortrows(maxima_storage);
max_intensity = sortedmax(end,1);
sortedmax(sortedmax(:,1) < max_intensity*rel_intensity_threshold,:) = [];

end

