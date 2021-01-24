function [tf] = isInCircle(x,y,x0,y0,r)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

rcalc = sqrt((x-x0).^2 + (y-y0).^2);
tf = (rcalc <= r);

end

