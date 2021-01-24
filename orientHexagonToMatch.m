function [rotated_vertices] = orientHexagonToMatch(hexagon_to_rotate,hexagon_to_match)
% assume that the unrotated hexagon is in the standard orientation produced
% by getHexagon.
TEST = 0;
% First, find the counterclockwise angle by which we should rotate. 
centered_htm = hexagon_to_match - mean(hexagon_to_match);
ordered_vertices = orderHexagonVertices(centered_htm);
% % bottom_right_vertex_centered = getHexagonVertex('bottom right',centered_htm);
bottom_right_vertex_centered = ordered_vertices(1,:); % by convention
theta = atan(bottom_right_vertex_centered(1)/bottom_right_vertex_centered(2));
[rotated_vertices] = rotateHexagon(hexagon_to_rotate,theta);

if TEST
    figure
    plot(hexagon_to_rotate(:,1),hexagon_to_rotate(:,2),'ko');
    hold on
    plot(hexagon_to_match(:,1),hexagon_to_match(:,2),'bo');
    plot(rotated_vertices(:,1),rotated_vertices(:,2),'ro');
end

end

