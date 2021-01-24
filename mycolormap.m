function [ cmap ] = mycolormap( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Originals
% c1hsv = [0.8,0.8,1];
% c2hsv = [0.8,1,0.4];

c2hsv = [0.8,1,0.4];
c1hsv = [0.4,1,0.6];


c1rgb = hsv2rgb(c1hsv);
c2rgb = hsv2rgb(c2hsv);

greyconst = 1;

% mycmap = [0.84,0.2,1;1,1,1;0.7,1,1];
mycmap = [c1rgb;greyconst,greyconst,greyconst;c2rgb];

seg1 = [linspace(mycmap(1,1),mycmap(2,1),50)',linspace(mycmap(1,2),mycmap(2,2),50)',linspace(mycmap(1,3),mycmap(2,3),50)'];
seg2 = [linspace(mycmap(2,1),mycmap(3,1),50)',linspace(mycmap(2,2),mycmap(3,2),50)',linspace(mycmap(2,3),mycmap(3,3),50)'];
cmap = vertcat(seg1,seg2);

end

