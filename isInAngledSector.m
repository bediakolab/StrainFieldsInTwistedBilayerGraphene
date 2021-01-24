function [mask] = isInAngledSector(xspace,yspace,origin,theta,arcangle)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

% need to work out the appropriate modular arithmetic here. atan2 is going
% to produce angles between [-pi, pi], which I have converted to [0,2pi]
% for ease of modding.

upperangle = mod(theta + arcangle/2,2*pi);
lowerangle = mod(theta - arcangle/2,2*pi);
[pixel_thetas] = getPixelSpaceAngle(xspace,yspace,origin);

% maybe the clearest way is to handle the cases separately:
if upperangle > lowerangle   % as expected if no mod related ambiguity has taken place
    mask = (pixel_thetas < upperangle) & (pixel_thetas > lowerangle);
else  % lowerangle has been increased by a factor of 2pi to avoid going negative.
    mask = (pixel_thetas < upperangle) | (pixel_thetas > lowerangle);
    % the or condition arises because anything smaller than upperangle is
    % bounded below by zero, which converts to 2pi, which is an upper bound
    % on stuff above lower angle.
end

end

