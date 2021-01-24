function [ordered_vertices] = orderHexagonVertices(vertices)
% index has the convention that "1" is the vertex with the angle closest to
% 

% get vertex angles
origin = mean(vertices,1);
[thetas] = getVertexAngle(vertices,origin);
combined = [thetas,vertices];
combined = sortrows(combined);
reference_angle = 3/2 * pi + 1/12 * pi;
thetadist = abs(combined(:,1) - reference_angle);
[~,min_idx] = min(thetadist);

ordered_vertices = zeros(6,2);
ordered_vertices(1,:) = combined(min_idx,2:3);


% iteratively find the remaining vertices
% because we applied sortrows to the angle first, we should now be able to
% simply pull off the next index, (idx-1) mod 6 + 1
idx = min_idx;
for i = 2:6
    % find the next largest angle after the one
    idx = idx + 1;
    ordered_vertices(i,:) = combined(mod(idx-1,6)+1,2:3);
end
% % 
% % switch location
% %     case 'bottom right'
% %         vertices = sortrows(vertices,2);
% %         vertices(4:6,:) = []; % delete the top ones.
% %         vertices = sortrows(vertices,1);
% %         vertex = vertices(3,:);
% %     case 
% % end        
% % 

end

