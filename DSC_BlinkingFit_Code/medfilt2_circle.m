function [ filtered_mat ] = medfilt2_circle( filtermat, filter_diameter )
% Function as a part of filterDisplacement.m
% Nothing for this but a for loop.
%
% Nathanael Kazmierczak, 05/20/2020

radius = filter_diameter/2;
[r,c] = size(filtermat);
rbase = 1:r;
cbase = 1:c;  % indices
[rspace,cspace] = meshgrid(rbase,cbase);
filtered_mat = zeros(size(filtermat));
for i = 1:r
    for j = 1:c
        r0 = i;
        c0 = j;
        [tf] = isInCircle(rspace,cspace,r0,c0,radius);
        medval = median(filtermat(tf));
        filtered_mat(i,j) = medval;
    end
end



end

