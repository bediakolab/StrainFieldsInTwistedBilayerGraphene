% make_fitting_functions.m
%
% Nathanael Kazmierczak, Dec 2019
% Using the results of the pseudostacking simulations.

load('12202019rasteraverages.mat');  % populates averages and raster_points
averages = averages-repmat(min(averages),size(averages,1),1);
averages = averages./repmat(max(averages),size(averages,1),1);
func_cell = cell(1,12);
clear si
for i = 1:12
    si = scatteredInterpolant(raster_points(:,1),raster_points(:,2),averages(:,i),'linear','nearest');
    thisfun = @(dsc_vector) si(dsc_vector(:,1),dsc_vector(:,2));
    func_cell{i} = thisfun;
end

