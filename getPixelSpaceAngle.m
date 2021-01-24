function [thetas] = getPixelSpaceAngle(xspace,yspace,origin)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here


xspace_centered = xspace - origin(1);
yspace_centered = yspace - origin(2);
% the four quadrant angle.
thetas = atan2(yspace_centered,xspace_centered);
thetas(thetas < 0) = thetas(thetas < 0) + 2*pi;
% so that angles range from [0, 2pi]

end

