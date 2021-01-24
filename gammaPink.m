function [ colormapmat ] = gammaPink( gamma )
% Function for applying a gamma filter to the pink colormap.
% Nathanael Kazmierczak, 06/16/2020

mat = pink;
n = size(mat,1);
current_scales = linspace(0,1,n)';
new_scales = current_scales.^gamma;
colormapmat = mat./current_scales.*new_scales;
colormapmat(isnan(colormapmat)) = 0;

end

