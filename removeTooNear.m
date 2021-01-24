function [unique_positions] = removeTooNear(combodata,cutoff_distance)
% operates on dim 2 and 3; This is because of the need to track
% intensities, which were col 1 for sortrows. Furthermore dim 4 is the
% window thta we are in

unique_positions = zeros(0,size(combodata,2));
for i = 1:size(combodata,1)
    candidate = combodata(i,:);
    candidate_location = candidate(:,2:3);
    candidate_intensity = candidate(:,1);
    if i == 1
        unique_positions(end+1,:) = candidate;
    else
        okay_to_append = 1;
        for j = 1:size(unique_positions,1)
            this_location = unique_positions(j,2:3);
            this_intensity = unique_positions(j,1);
            distance = norm(this_location - candidate_location);
            if distance < cutoff_distance
                if candidate_intensity <= this_intensity
                    okay_to_append = 0;
                    break;
                else   % the new one is worthy of replacing the old
                    % It would be possible to have it be able to replace
                    % multiple, but hopefully that won't happen esp given
                    % that we have presorted, so it shouldn't happen
                    % anyways.
                    unique_positions(j,:) = candidate;
                    okay_to_append = 0;  % because we are not appending but replacing.
                    break;
                end
            end    
        end
        
        if okay_to_append
            unique_positions(end+1,:) = candidate;
        end
    end
end

end

