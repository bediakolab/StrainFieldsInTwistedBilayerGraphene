% testHexagonalLatticeMask.m
xbase = -2:0.01:10;
ybase = -2:0.01:10;
[xspace,yspace] = meshgrid(xbase,ybase);

x0_init = 0;
y0_init = 0;
d = 1;
xnum = 10;
ynum = 10;
[points, centers] = getHexagonalLattice(x0_init,y0_init,d,xnum,ynum);


[hexmask] = generateHexagonalLatticeMask(xspace,yspace,centers(:,1),centers(:,2),0.3);

circlemask = isInCircle(xspace,yspace,2.3,2*sqrt(3),0.5);
bothmask = hexmask + 2*circlemask;

figure;
pcolor(xspace,yspace,hexmask); shading flat;

figure;
pcolor(xspace,yspace,bothmask); shading flat;
title 'Beamspot radius= 0.5, Reconstruction "radius" = 0.3'
xlabel('x')
ylabel('y')
