function [ RGB_color_stack ] = getTriangleColorLegend(density,AB_hsv,SP_hsv)
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    SP_hsv = [0.8,0.8,1];
    AB_hsv = [0.5,0.3,1];
end

vec1 = [1;0];
t = pi/3;
rotmat = [cos(t) -sin(t); sin(t) cos(t)];
vec2 = rotmat*vec1;
vec3 = rotmat*vec2;

v1 = vectorLine(vec1',[0,0]);
v2 = vectorLine(vec2',[0,0]);
v3 = vectorLine(vec3',[1,0]);

% density = 0.01;
xbase = 0:density:1;
ybase = 0:density:1;

[xspace,yspace] = meshgrid(xbase,ybase);
mask = false(size(xspace));
for i = 1:size(xspace,2)
    thispointcol = [xspace(:,i),yspace(:,i)];
    c1 = v1 <= thispointcol;
    c2 = v2 >= thispointcol;
    c3 = v3 >= thispointcol;
    mask(c1 & c2 & c3,i) = 1;
end

angles = atan(yspace./xspace);
amplitudes = sqrt(xspace.^2 + yspace.^2);
max_amplitudes = sqrt(3)./(2*sin(2*pi/3 - angles));
Values = amplitudes./max_amplitudes;
Hues = angles/(pi/3)*AB_hsv(1) + (1-angles/(pi/3))*SP_hsv(1);
Saturations = angles/(pi/3)*AB_hsv(2) + (1-angles/(pi/3))*SP_hsv(2);

Values(~mask) = 1;
Hues(~mask) = 0;
Saturations(~mask) = 0;

HSV_color_stack = cat(3,Hues,Saturations,Values);
RGB_color_stack = hsv2rgb(HSV_color_stack);

end

