function [ A_merged ] = merge4thdim( A )
% A input should be a 4-dimensional datacube.
% A_merged will be a 3-dimensional datacube, where the 3rd dimension is
% preserved contiguously, and the fourth dimension changes in small
% increments.

[d1,d2,d3,d4] = size(A);
A_merged = zeros(d1,d2,d3*d4);
for i = 1:d4
    idxstart = (i-1)*d3+1;  % because we are extracting all d3s from a set value of d4.
    idxend = d3*i;
    A_merged(:,:,idxstart:idxend) = A(:,:,:,i);
end

