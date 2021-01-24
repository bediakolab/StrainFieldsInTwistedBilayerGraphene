function [ cmap ] = blackredcolormap( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

len1 = 50;
hval = 0.98;
c2hsv = [hval,1,1];
c1hsv = [hval,1,0];

cmap = [linspace(c1hsv(1),c2hsv(1),len1)',linspace(c1hsv(2),c2hsv(2),len1)',linspace(c1hsv(3),c2hsv(3),len1)'];
cmap = hsv2rgb(cmap);

end

