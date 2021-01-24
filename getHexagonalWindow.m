function [sectormasks,hexagon_vertices] = getHexagonalWindow(graphdata,arcangle,annular_increment)
% Function calls ginput
disp('This function is for developing windows around clusters of Bragg disks for detection.');
disp('Please choose a set of peaks to center the window around. Click on disk centers.');
disp('Please proceed in a single direction around the hexagon.');

figure
pcolor(log10(double(graphdata))); shading flat; colormap(bone);
manual_vertices = ginput(6);

xbase = 1:size(graphdata,2);
ybase = 1:size(graphdata,1);
[xspace,yspace] = meshgrid(xbase,ybase);

% (1) Determine the center by fitting a regular hexagon to the chosen
% vertices. Note this will not be the center for any actual strain mapping,
% but rather for defining the origin for drawing lines.
d1 = sqrt(sum((manual_vertices(1,:) - manual_vertices(4,:)).^2));
d2 = sqrt(sum((manual_vertices(2,:) - manual_vertices(5,:)).^2));
d3 = sqrt(sum((manual_vertices(3,:) - manual_vertices(6,:)).^2));
d_guess = mean([d1,d2,d3])/2;
x0_guess = mean(manual_vertices(:,1));
y0_guess = mean(manual_vertices(:,2));
theta_guess = 0;
c_guess = [x0_guess, y0_guess, d_guess,theta_guess];


% hexpred = @(c) sortrows(orientHexagonToMatch(getHexagon(c(1),c(2),c(3)),manual_vertices));
% A better guess may allow full rotational freedom
hexpred = @(c) rotateHexagon(orderHexagonVertices(getHexagon(c(1),c(2),c(3))),c(4));
hexdev = @(c) rms(rms(hexpred(c) - orderHexagonVertices(manual_vertices)));
options = optimset;
options.Display = 'iter';
optimparams = fminsearch(hexdev,c_guess,options)
hexagon_vertices = hexpred(optimparams);
initial_hexagon_vertices = hexpred(c_guess);

figure
plot(hexagon_vertices(:,1),hexagon_vertices(:,2),'ro-');
hold on
plot(manual_vertices(:,1),manual_vertices(:,2),'ko-');
plot(initial_hexagon_vertices(:,1),initial_hexagon_vertices(:,2),'go-');
legend('hexagon fit','manual','initial guess');

x0 = optimparams(1);
y0 = optimparams(2);
circle_radius = sqrt(sum((hexagon_vertices(1,:) - [x0,y0]).^2));
% annular_increment = 10;
[annulus_mask] = isInAnnulus(xspace,yspace,x0,y0,circle_radius-annular_increment,circle_radius+annular_increment);
% produce six lines oriented to the six vertices of the hexagon.
origin = [x0,y0];
thetas = getVertexAngle(hexagon_vertices,origin);
% raymasks = zeros(size(xspace,1),size(xspace,2),6);s
sectormasks = zeros(size(xspace,1),size(xspace,2),6);
% thickness = 40;
% arcangle = pi/12;
% for i = 1:6
%     raymasks(:,:,i) = isInAngledThickRay(xspace,yspace,origin,thetas(i),thickness) & annulus_mask;
% end
for i = 1:6
    sectormasks(:,:,i) = isInAngledSector(xspace,yspace,origin,thetas(i),arcangle) & annulus_mask;
end

figure; pcolor(sum(sectormasks,3)); shading flat;
disp('Exiting function');

end

