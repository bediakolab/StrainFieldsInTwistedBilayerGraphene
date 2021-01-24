% testIsInCircle.m

xbase = -2:0.01:10;
ybase = -2:0.01:10;
[xspace,yspace] = meshgrid(xbase,ybase);

x0 = 0;
y0 = 0;
r = 1;

[tf] = isInCircle(xspace,yspace,x0,y0,r);
figure;
pcolor(xspace,yspace,double(tf)); shading flat;



