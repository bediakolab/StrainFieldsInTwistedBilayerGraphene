function [mask] = isInAngledThickRay(xspace,yspace,origin,theta,thickness)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here


x1 = cos(theta) + origin(1);
y1 = sin(theta) + origin(2);
m = (origin(2)-y1)/(origin(1) - x1);

mask = (yspace > m*(xspace-x1) + y1 - thickness) & (yspace < m*(xspace-x1) + y1 + thickness);
% Should refine at somepoint to condition on whichever is larger.
if cos(theta) > 0
    mask(xspace - origin(1) < 0) = 0;
else
    mask(xspace - origin(1) > 0) = 0;
end

if sin(theta) > 0
    mask(yspace - origin(2) < 0) = 0;
else
    mask(yspace - origin(2) > 0) = 0;
end

end

