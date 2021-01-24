function [ blinking_matrix ] = makeBlinkingPlot( disk_average, DSC_coords )
% Given a set of data that represents a hexagonally rastered blinking
% pattern, visualize their progression across the DSC lattice. Need to pass
% in only a single column vector for disk average.
%
% Nathanael Kazmierczak, Dec 2019.

hexagon_lattice_constant = 2.461;

% Turn off extrapolation so that we will get nan outside of the hexagonal
% hull. Later, a hexagonal tiling will make this much smoother.
disk_average = disk_average./max(disk_average);
blinkingInterp = scatteredInterpolant(DSC_coords(:,1),DSC_coords(:,2),disk_average,'nearest','none');
% Evaluate at the same coordinates as in the cell.
grid_density = 0.1;
t = 60/180*pi;
rotmat = [cos(t) sin(t); -sin(t) cos(t)];
v1 = [0,hexagon_lattice_constant/sqrt(3)];
v2 = v1*rotmat';
PADDING = 0.1;
xbase = (-max([v1(1),v2(1)])-PADDING):grid_density:(max([v1(1),v2(1)])+PADDING);
ybase = (-max([v1(2),v2(2)])-PADDING):grid_density:(max([v1(2),v2(2)])+PADDING);
[xspace,yspace] = meshgrid(xbase,ybase);
blinking_matrix = blinkingInterp(xspace,yspace);

figure;
% set(gcf,'color',[0.277777777777778,0.173600000000000,0.110555555555556]);
contourf(xbase,ybase,blinking_matrix,20,'LineStyle','None');
colormap(gray);
colorbar
axis equal

end

