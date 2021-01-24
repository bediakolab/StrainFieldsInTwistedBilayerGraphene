function [ filtered_mat ] = movmeanfilt2_circle( filtermat, filter_diameter )
% Function as a part of filterDisplacement.m
% Companion to medfilt2_circle.m, but this one can be done with fft.
%
% Nathanael Kazmierczak, 05/20/2020

if mod(filter_diameter,2) == 0
    filter_diameter = filter_diameter + 1;  % needs to be odd
    warning('Inflating filter diameter for parity in movmeanfilt2_circle()');
end

radius = filter_diameter/2;
rbase = 1:filter_diameter;
cbase = 1:filter_diameter;  % indices
[rspace,cspace] = meshgrid(rbase,cbase);
[tf] = isInCircle(rspace,cspace,(filter_diameter+1)/2,(filter_diameter+1)/2,radius);
H = zeros(size(tf));
H(tf) = 1;
numpixelsmat = ones(size(filtermat));

imsum = filter2(H,filtermat);
imcounts = filter2(H,numpixelsmat);
averages = imsum./imcounts;
filtered_mat = averages;

end

