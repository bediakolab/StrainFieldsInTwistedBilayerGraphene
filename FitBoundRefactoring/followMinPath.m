function [ BW_final ] = followMinPath( matrix, BW_init, starting_point )
% In preparation for the domain pinning to reconstruct the extended zone
% displacement field.

AA_thresh = 1;

[r,c] = size(matrix);
% convert to indices
this_point = fliplr(round(starting_point));
BW = BW_init;
count = 0;
while true
    count = count + 1;
    this_val = matrix(this_point(1),this_point(2));
    BW(this_point(1),this_point(2)) = 1;
    % Is this point on the boundary? 
    if this_point(2) == c || this_point(2) == 1 || this_point(1) == r || this_point(1) == 1
        break
    end
    
    % Get neighbor values
    neighbors = matrix(this_point(1)-1 : this_point(1)+1,this_point(2)-1 : this_point(2)+1);
    neighborsBW = BW(this_point(1)-1 : this_point(1)+1,this_point(2)-1 : this_point(2)+1);
    diffvalues = neighbors - this_val;
    
    % Have we reached a local minimum at an AA stacking region?
    is_local_min = all(all(diffvalues >= 0))
    is_AA = this_val < AA_thresh
    if is_local_min && is_AA
        break
    end
    
    TEST = 1;
    if TEST 
        hold on;
        spy(BW);
        set(gca,'ydir','normal');
    end
    
    % If termination conditions have not been met, move to the nearest low
    % point that has not already been visited.
    diffvalues(neighborsBW) = 1000;
    [~,idx] = min(diffvalues(:));
    [ridx,cidx] = ind2sub(size(diffvalues),idx);
    radd = ridx - 2;
    cadd = cidx - 2;
    this_point = this_point + [radd,cadd];    
end

BW_final = BW;

end

