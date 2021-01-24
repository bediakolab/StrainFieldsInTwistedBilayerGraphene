function [ output_args ] = scatterSpherical( coords, scaling, color, figh )
% Function for plotting spheres in scatter plot
%
% coords should be a nx3 array giving the center coordinates of the sphere.
% color should be a 1x3 RGB array in range 0 to 1.
%
% Nathanael Kazmierczak, 10/30/2020
SPHERERES = 5;

figure(figh)
hold on
[x,y,z] = sphere;
x = x*scaling;
y = y*scaling;
z = z*scaling;

for i = 1:size(coords,1)
    x0 = coords(i,1);
    y0 = coords(i,2);
    z0 = coords(i,3);
    xf = x + x0;
    yf = y + y0;
    zf = z - z0;
    surf(xf,yf,zf,'EdgeColor','none','FaceColor',color);
    
    
end


end

