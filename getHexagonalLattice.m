function [points, centers] = getHexagonalLattice(x0,y0,d,xnum,ynum)
%xnum, ynum are the number of contiguous hexagons, because this is the easiest
%way to think of it. x0 y0 are central coordinates of the bottom left hand hexagon. 

xgen_centers = x0:4*d:(x0+xnum*2*d);
ygen_centers = y0:d*sqrt(3):(y0+ynum*d*sqrt(3));
xcenters = x0:2*d:(x0+xnum*2*d);
ycenters = ygen_centers;
realcenters = [xgen_centers', repmat(y0,numel(xgen_centers),1)];
for i = 2:numel(ygen_centers)
    if ~mod(i,2)
        xshift = d;
    else
        xshift = 0;
    end
    realcenters = [realcenters; [(xgen_centers+xshift)', repmat(ygen_centers(i),numel(xgen_centers),1)]];
end

points = [];
for i = 1:size(realcenters,1)
    points = [points; getHexagon(realcenters(i,1),realcenters(i,2),d)];
end

[a,b] = meshgrid(xcenters,ycenters);
centers = [];
for i = 1:size(a,1)
    if mod(i,2)
        centers = [centers; [a(i,:)',b(i,:)']];
    else
        centers = [centers; [a(i,:)'+d,b(i,:)']];
    end
end

end

