% custom_displacement_color_testing_script.m

grid_density = 0.01;
% [ raster_points ] = getDSCHexagonRaster(grid_density,2.461);
xbase = -1.5:grid_density:1.5;
ybase = 0:grid_density:1.5;
[xspace,yspace] = meshgrid(xbase,ybase);
catmat = cat(3,xspace,yspace);

[ RGB_color_stack ] = getCustomDisplacementColor( catmat );
figure;
imagesc(xbase,ybase,RGB_color_stack);
set(gca,'ydir','normal');
axis equal
title('2D Colormap');
xlabel('x displacement');
ylabel('y displacement');

load('objectdataRadialResiduals.mat');
[ RGB_color_stack ] = getCustomDisplacementColor( m4.DSC_fit_storage );
figure
subplot(1,2,1);
imagesc(RGB_color_stack);
axis square
xlabel('Realspace x');
ylabel('Realspace y');
title('Colorized AB, SP, and AA regions');

SP_hsv = [0.8,0.8,1];
AB_hsv = [0.5,0.3,1];
[ RGB_color_stack ] = getTriangleColorLegend(grid_density,AB_hsv,SP_hsv);
subplot(1,2,2);
imagesc(RGB_color_stack);
