function [ cmap ] = whiteburgundycolormap( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

len1 = 50;
% hval = 0. 1;
c2hsv = ([128,30,0]+1)./256*1;
c1hsv = [1,1,1];

cmap = [linspace(c1hsv(1),c2hsv(1),len1)',linspace(c1hsv(2),c2hsv(2),len1)',linspace(c1hsv(3),c2hsv(3),len1)'];
% cmap = cmap;

end

