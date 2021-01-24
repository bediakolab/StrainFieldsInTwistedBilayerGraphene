function [f] = hs024_obj(x)
% Objective function should be prepared to receive a set of points

f=((x(1,:)-3).^2-9).*x(2,:).^3/(27*sqrt(3));

return