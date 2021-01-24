% validateBlinkingFittingFunction.m

grid_density = 0.3;
hexagon_lattice_constant = 2.461;
[ raster_points ] = getDSCHexagonRaster(grid_density,hexagon_lattice_constant);
results = blinkingFittingFunction(raster_points);

% These are matching exactly with what we got straight off of the
% simulation.
makeBlinkingPlot( results(:,1), raster_points );
makeBlinkingPlot( results(:,10), raster_points );
