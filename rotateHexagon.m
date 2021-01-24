function [rotated_vertices] = rotateHexagon(initial_vertices,theta)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

TEST = 0;

rotmat = [cos(theta), -sin(theta); sin(theta), cos(theta)];
meanvals = mean(initial_vertices);
centered_init = initial_vertices - meanvals;
centered_rotated_vertices = centered_init*rotmat;
rotated_vertices = centered_rotated_vertices + meanvals;

if TEST
    figure
    plot(initial_vertices(:,1),initial_vertices(:,2),'ko');
    hold on
    plot(rotated_vertices(:,1),rotated_vertices(:,2),'bo');
end

end

