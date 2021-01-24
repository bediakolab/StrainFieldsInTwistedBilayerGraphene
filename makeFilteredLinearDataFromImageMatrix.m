function [linear_filtered_mat] = makeFilteredLinearDataFromImageMatrix(mat)
% Assuming the filtering flag is -1.
% Unwrap the matrices in identical fashion after having applied meshgrid
% (ind2sub would work too).
xbase = 1:size(mat,2);
ybase = 1:size(mat,1);
[xspace,yspace] = meshgrid(xbase,ybase);
linear_matrix = [mat(:),yspace(:),xspace(:)];
linear_filtered_mat = linear_matrix;
linear_filtered_mat(linear_filtered_mat(:,1) < 0,:) = [];


end

