function [ raster_points ] = getDSCHexagonRaster(grid_density,hexagon_lattice_constant)
    % Returns a Cartesian sampling of the hexagonal DSC lattice
    % I.e., those points lying inside the boundaries of the hexagon.
    % 
    % Nathanael Kazmierczak, Dec 2019

    PLOT_FLAG = 0;
    
    t = 60/180*pi;
    rotmat = [cos(t) sin(t); -sin(t) cos(t)];
    v1 = [0,hexagon_lattice_constant/sqrt(3)];
    v2 = v1*rotmat';
    
    % h defines all of the extremal vertices on the DSC lattice.
    h = zeros(6,2);
    h(1,:) = v1;
    h(2,:) = v2;
    h(3,:) = v2-v1;
    h(4,:) = -v1;
    h(5,:) = -v2;
    h(6,:) = v1-v2;
    
    vl1 = vectorLine(v2-v1,v1);
    vl2 = vectorLine(v1,v2);
    vl3 = vectorLine(v2,-v1);
    vl4 = vectorLine(v2-v1,-v1);
    vl5 = vectorLine(v1,-v2);
    vl6 = vectorLine(-v2,v1);
    
    PADDING = 0.1;
    if numel(grid_density) == 1
    xbase = (-max([v1(1),v2(1)])-PADDING):grid_density:(max([v1(1),v2(1)])+PADDING);
    ybase = (-max([v1(2),v2(2)])-PADDING):grid_density:(max([v1(2),v2(2)])+PADDING);
    else
        xpad = 0.001;
        ypad = 0.001;
    xbase = linspace(-v2(1)+xpad,v2(1)-xpad,grid_density(1));
    ybase = linspace(0+ypad,v1(2)-ypad,grid_density(2));
    end
    
    [xspace,yspace] = meshgrid(xbase,ybase);
    raster_points = [xspace(:),yspace(:)];
    
    raster_points(vl1 <= raster_points,:) = [];
    raster_points(vl2 <= raster_points,:) = [];
    raster_points(vl3 >= raster_points,:) = [];
    raster_points(vl4 >= raster_points,:) = [];
    raster_points(vl5 >= raster_points,:) = [];
    raster_points(vl6 <= raster_points,:) = [];

    if PLOT_FLAG
        figure;
        scatter(raster_points(:,1),raster_points(:,2));
        xlabel('x');
        ylabel('y');
        axis equal
    end
    
end
