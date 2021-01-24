function [ legend_handles ] = makeConcentricRings( axh )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

center = [0,0];
radii = 0.25:0.25:1.5;
colors = linspace(0.9,0,numel(radii));
axes(axh);
for i = 1:numel(radii)
    legend_handles(i) = viscircles(center,radii(i),'Color',repmat(colors(i),1,3));
end


end

