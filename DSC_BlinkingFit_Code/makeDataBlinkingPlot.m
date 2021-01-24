function makeDataBlinkingPlot( disk_average )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

disk_average = disk_average./max(disk_average);

figure;
contourf(disk_average,20,'LineStyle','None');
colormap(gray);
colorbar
axis equal


end

