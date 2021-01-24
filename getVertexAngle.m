function [theta] = getVertexAngle(vertices,origin)
% Assumes that the first coordinate is x, the second coordinate is y


centered = vertices - origin;
% the four quadrant angle.
theta = atan2(centered(:,2),centered(:,1));

end

