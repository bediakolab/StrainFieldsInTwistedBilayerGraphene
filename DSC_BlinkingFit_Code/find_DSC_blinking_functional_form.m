% find_DSC_blinking_functional_form.m

load('12202019rasteraverages.mat');

figure;
scatter(raster_points(:,1),raster_points(:,2));
xlabel('x');
ylabel('y');
axis equal

averages = averages - min(averages);
averages = averages./repmat(max(averages),size(averages,1),1);
raster_points = horzcat(raster_points,averages);

fixed_x1 = raster_points(abs(raster_points(:,2) - -0.020859012475669)<1e-6,:);
fixed_x2 = raster_points(abs(raster_points(:,2) - 0.279140987524331)<1e-6,:);
fixed_x3 = raster_points(abs(raster_points(:,2) - 0.579140987524331)<1e-6,:);
fixed_x4 = raster_points(abs(raster_points(:,2) - -0.320859012475669)<1e-6,:);
fixed_x5 = raster_points(abs(raster_points(:,2) - -0.620859012475669)<1e-6,:);

figure;
hold on
scatter(fixed_x1(:,1),fixed_x1(:,3));
scatter(fixed_x2(:,1),fixed_x2(:,4));
scatter(fixed_x3(:,1),fixed_x3(:,5));
scatter(fixed_x4(:,1),fixed_x4(:,6));
scatter(fixed_x5(:,1),fixed_x5(:,7));

