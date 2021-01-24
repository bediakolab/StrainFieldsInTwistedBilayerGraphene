% trblg_raster_driver.m
%
% Nathanael Kazmierczak, Dec 2019

input('HAAAAALLLLLTTTTT! Make sure you are not about to overwrite data.');
savename = 'test.mat';

load('12202019multislice_disk_centers.mat');  % will populate disk_centers and radius for the disk averaging.
addpath(genpath(pwd));

moire_cell_to_match = 1.2;
% grid_density = 0.2;
grid_density = 0.3;
hexagon_lattice_constant = 2.461;
[ raster_points ] = getDSCHexagonRaster(grid_density,hexagon_lattice_constant);
clear trblg_raster;

% 12 when using '12202019multislice_disk_centers.mat' file
disks_used = size(disk_centers,1);
averages = zeros(size(raster_points,1),disks_used);
mytimer = tic;
for i = 1:size(raster_points,1)
    fprintf('Beginning multislice iteration %d out of %d\n',i,size(raster_points,1));
    fprintf('\tElapsed time is %f seconds\n',toc(mytimer));
    thisDSCvector = raster_points(i,:);
    trblg_raster(i) = TranslatedBilayerGraphene(thisDSCvector,moire_cell_to_match);
    trblg_raster(i).simulate();
    trblg_raster(i).setAveragingCircles(disk_centers,radius);
    thisaverages = trblg_raster(i).getAverageDiskIntensities();
    averages(i,:) = thisaverages';
    save(savename);
end

% plot through interpolation. Recall that raster_points give the DSC
% coordinates
makeBlinkingPlot( averages(:,1), raster_points );
title 'Inner disk 1'
makeBlinkingPlot( averages(:,2), raster_points );
title 'Inner disk 2'
makeBlinkingPlot( averages(:,3), raster_points );
title 'Inner disk 3'
makeBlinkingPlot( averages(:,4), raster_points );
title 'Inner disk 4'
makeBlinkingPlot( averages(:,5), raster_points );
title 'Inner disk 5'
makeBlinkingPlot( averages(:,6), raster_points );
title 'Inner disk 6'
makeBlinkingPlot( averages(:,7), raster_points );
title 'Outer disk 1'
makeBlinkingPlot( averages(:,8), raster_points );
title 'Outer disk 2'
makeBlinkingPlot( averages(:,9), raster_points );
title 'Outer disk 3'
makeBlinkingPlot( averages(:,10), raster_points );
title 'Outer disk 4'
makeBlinkingPlot( averages(:,11), raster_points );
title 'Outer disk 5'
makeBlinkingPlot( averages(:,12), raster_points );
title 'Outer disk 6'

save(savename);

