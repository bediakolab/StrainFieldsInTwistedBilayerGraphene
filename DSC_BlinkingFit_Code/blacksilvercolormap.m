function [ cmap ] = blacksilvercolormap( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

len1 = 50;
% hval = 0. 1;
c2hsv = ([168,169,173]+1)./256.*1.3;
c1hsv = [0,0,0];

cmap = [linspace(c1hsv(1),c2hsv(1),len1)',linspace(c1hsv(2),c2hsv(2),len1)',linspace(c1hsv(3),c2hsv(3),len1)'];
% cmap = cmap;

end

