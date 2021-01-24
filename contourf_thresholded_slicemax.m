function contourf_thresholded_slicemax(X,Y,Z,threshold,linenumber)
% Some plotting code for Nathanael's convolution simulations.
% Also provides some nice thresholding boilerplate code.

% detect slicewise local minima in slices for fixed first dimension
num_y = size(Z,1); % where it is a matrix
local_max_storage = zeros(0,3); % for xyz storage
for i = 1:num_y
    this_y = Y(i);
    slice = Z(i,:);
    maxima_indices = find(islocalmax(slice));
    maxima_indices(slice(maxima_indices) < threshold) = [];
    intensities = slice(maxima_indices);
    for j = 1:numel(maxima_indices)
        local_max_storage(end+1,:) = [maxima_indices(j),this_y,intensities(j)];
    end
end

% figure, axes will be built outside of this function
contourf(X,Y,Z,linenumber);
hold on
scatter(local_max_storage(:,1),local_max_storage(:,2),'filled')

end

