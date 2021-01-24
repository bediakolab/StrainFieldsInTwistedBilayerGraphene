function [ output_mat ] = filterMedianOutlier( input_mat, neighborhood_sidelength, threshold )
% Wrapper function for medfilt2, but only allows the replacement if the
% matrix value deviates by more than threshold from the median value.
%
% Nathanael Kazmierczak, 05/19/2020
% Modified 07/18/2020 to correct the problem where pixels in the corners
% were getting set to zero. 

medfilt_result = medfilt2(input_mat,[neighborhood_sidelength,neighborhood_sidelength]);
diff_mag = abs(medfilt_result - input_mat);
output_mat = input_mat;
logicals = (diff_mag > threshold);
output_mat(logicals) = medfilt_result(logicals);
n_replacements = nnz(logicals);
% Because of the quirk on trimArray, need to set 0 to get a one pixel
% boundary. Don't want filtering to take place around the boundary of the
% image, because the median outlier filter tends to loose pixels in the
% corner of the image. 
edge_protection_mask = trimArray(false(size(input_mat)),0,true(size(input_mat)));
output_mat(edge_protection_mask) = input_mat(edge_protection_mask);
fprintf('filterMedianOutlier() made %d replacements (%.2f%%) at a threshold of %.2f and a sidelength of %d.\n',...
    n_replacements,n_replacements/numel(logicals)*100,threshold,neighborhood_sidelength);

end

