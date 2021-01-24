function plotDPandDisks( DP, disk_locations )
% Pass disk locations in IJ format.

figure; image(uint8(DP),'CDataMapping','Scaled'); colormap(hsv); colorbar; pbaspect([1 1 1]);
hold on; 
scatter(disk_locations(:,2),disk_locations(:,1),'filled'); colormap(pink);

end

