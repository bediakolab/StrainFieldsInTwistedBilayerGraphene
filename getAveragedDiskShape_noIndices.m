function [ average_probe_convolution, average_probe_stacking ] = getAveragedDiskShape_noIndices( DPROIs, subpixel, max_convolve )
% Function principally for use when probing the interference patterns
% present in the dataset 2.
%
% In the usage we want, the DPROIs have already been fully constructed of
% the proper shape in the driver script.
%
% Nathanael Kazmierczak, 09032019 

% We may need to thin these down for a reasonable convolution time. 
max_convolve


% Leaving this as a function just in case we need to do anything else here.
[ average_probe_convolution, average_probe_stacking ] = getAveragedDiskShape_innards( DPROIs, subpixel );

end

