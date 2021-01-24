% trblg_driver.m
%
% Basic demonstration of trblg class.
% Nathanael Kazmierczak, Dec 2019

DSC_vector_or_type = 'AA';
moire_cell_to_match = 1.2;
trblg = TranslatedBilayerGraphene(DSC_vector_or_type,moire_cell_to_match);
trblg.plot();

% trblg.computeDSCField();
% figh = trblg.plotDSCField();
% trblg.plotDSCLattice();
% figh = trblg.assignPsuedostacking(figh);

grid_density = 0.2;
hexagon_lattice_constant = 2.461;
[ raster_points ] = getDSCHexagonRaster(grid_density,hexagon_lattice_constant);
clear trblg_raster;
for i = 1:size(raster_points,1)
    thisDSCvector = raster_points(i,:);
    trblg_raster(i) = TranslatedBilayerGraphene(thisDSCvector,moire_cell_to_match);
end


stack4D = trblg_raster(17).simulate();


plotDP(stack4D);
xlim([2.263892128279883e+02 3.706107871720116e+02]);
ylim([2.272580174927113e+02 3.680043731778425e+02]);
if ~exist('disk_centers','var')
    disk_centers = zeros(6,2);
    input('Zoom graph to prepare to click the disks.');
    disp('Click once for the first disk (top vertical position) and a second time to set the circular radius.');
    setcirc = ginput(2);  
    radius = sum((setcirc(1,:) - setcirc(2,:)).^2).^0.5;
    disk_centers(1,:) = setcirc(1,:);
    
    for i = 1:5
        input('Zoom graph to prepare to click the next disk center.');
        fprintf('Click once to set the center of the %dth disk.\n',i+1);
        disk_centers(i+1,:) = ginput(1); 
    end
    viscircles(disk_centers,radius*ones(size(disk_centers,1),1));
end



