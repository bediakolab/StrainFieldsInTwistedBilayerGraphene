function [ roicut, Irangecut, Jrangecut ] = cutROI_fromPosition( DP, central_position, window_width )
% Refactoring of code originally found in cutROI, when we began to want to
% do ROI cuts on things other than registered Bragg disks (the interference
% rods).

% central_position is in I,J format

% We will use the convention that if window_width is even, we will add the
% extra pixel to the top and right of the disk center.
if mod(window_width,2) % easy case
    increment = (window_width-1)/2;
    % Forming square rois
    Irangecut = (central_position(1)-increment):(central_position(1)+increment);
    Jrangecut = (central_position(2)-increment):(central_position(2)+increment);
else
    increment = (window_width-2)/2;
    Irangecut = (central_position(1)-increment):(central_position(1)+increment+1);
    Jrangecut = (central_position(2)-increment):(central_position(2)+increment+1);
end

roicut = DP(Irangecut,Jrangecut);

% This still needs to be tested better.

end

