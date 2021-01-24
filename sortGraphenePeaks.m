function [first_graphene_disks,second_graphene_disks] = sortGraphenePeaks(graphene_locations,graphene_intensities,origin)
% Given some bragg disks registered for twisted bilayer, use knowledge of
% the angle from the center to sort into first and second graphene sets.
%
% We are not going to worry about weighted least squares just yet.
%
% The present algorithm is going to get really quite boogered up if we fail
% a disk detection. In such a situation, we will need to omit the real
% space pixel from consideration.
%
% Recall that graphene_locations come in in row, column

testnan = graphene_locations(:,1);
numnan = sum(isnan(testnan));
if numnan > 2
    first_graphene_disks = [];
    second_graphene_disks = [];
    return
else
    graphene_locations(isnan(testnan),:) = [];
    graphene_intensities(isnan(testnan),:) = [];
end
    
icombo = [graphene_locations,graphene_intensities];  % not sure we actually need the intensities in the end.
isort = sortrows(icombo);
% Under the new way of doing this, we are now sorting on hexagonal window
% index.


% % % % First get rid of those two stragglers from behind the beamstop
% % % icombo = [graphene_intensities,graphene_locations];
% % % isort = sortrows(icombo);
% % % isort(1:2,:) = [];  % Because these are the bit behind the beamstop.
% 
% if any(isnan(graphene_locations))
%     first_graphene_disks = [];
%     second_graphene_disks = [];
% end


%%% NPK 09/01/2019: Something has gone wrong here on version control --
%%% this is still getting called, but it is set up to handle an old version
%%%. 
graphene_locations = isort(:,2:3);
% graphene_locations = isort(:,1:2);  
% We know from the way these are generated that 1:2 belong to a certain
% window, 3:4 belong to a certain window, etc.
%
%  Note that I think the original fliplr below was wrong -- graphene peaks
%  were always in xcoord, ycoord format by the time we got here. So we
%  might get things that are the other way around now.
[thetas] = getVertexAngle(graphene_locations,origin);
thetas(thetas < 0) = thetas(thetas < 0) + 2*pi;
% so that angles range from [0, 2pi]

num_windows = size(graphene_locations,1)/2;
bottom_counter = 1;
top_counter = 2;
num_disks = numel(thetas);
first_graphene_disks = zeros(num_disks/2,2);
second_graphene_disks = zeros(num_disks/2,2);

for i = 1:num_windows
    % graphene_locations are still sorted by window.
    this_sorted = isort(bottom_counter:top_counter,:);
    these_thetas = thetas(bottom_counter:top_counter);
    % now that we have the bit beloning to each window, we have to sort
    % them by angle again
    windowcombo = [these_thetas,this_sorted];
    windowcombo = sortrows(windowcombo);
    if (0 <= windowcombo(1,1) && windowcombo(1,1) <= pi/2) && (3*pi/2 <= windowcombo(2,1) && windowcombo(2,1) <= 2*pi)
        % in such a case sorted(2,:) actually belongs to the first graphene
        % group as we are traversing counterclockwise.
        first_graphene_disks(i,:) = windowcombo(2,3:4);  % the first column being which of the windows it came from
        second_graphene_disks(i,:) = windowcombo(1,3:4);
    else % standard ordering: more advanced angle is the second graphene
        first_graphene_disks(i,:) = windowcombo(1,3:4);
        second_graphene_disks(i,:) = windowcombo(2,3:4);
    end
    top_counter = top_counter + 2;
    bottom_counter = bottom_counter + 2;
end

%% Delete any window where a detection was lost
% If a window only came up with one disk when it was supposed to have two,
% then it will be very hard to say computationally which was supposed to go
% with which set of graphene disks.
% Therefore if any row has a nan, delete the entire row from both disk
% storages, as both are unreliable.

todelete = [];
for i = 1:size(first_graphene_disks,1)
    if any(isnan(first_graphene_disks(i,:))) || any(isnan(second_graphene_disks(i,:)))
        todelete = [todelete, i];
    end
end
first_graphene_disks(todelete,:) = [];
second_graphene_disks(todelete,:) = [];

% In principle all nans should die at this stage.

% % % % % % if the first and the last should go together, this will actually cause
% % % % % % problems, so fix that here.
% % % % % if abs((sorted(end,1)-2*pi) - sorted(1,1)) < abs(sorted(end,1) - sorted(end-1,1))
% % % % %     sorted(end+1,:) = sorted(1,:) + [2*pi,0,0];
% % % % %     sorted(1,:) = [];
% % % % % end

% % % % % % Since the discontinuity happens at 2pi, it should be enough to ask
% % % % % % whether one angle is in 0 < angle < pi/2 and another is in 3pi/2 < 2pi.
% % % % % % Then we need to flip.
% % % % % num_disks = numel(thetas);
% % % % % uppercounter = 2;
% % % % % first_graphene_disks = zeros(num_disks/2,2);
% % % % % second_graphene_disks = zeros(num_disks/2,2);
% % % % % for i = 1:num_disks/2
% % % % % %     disks = graphene_locations(i:uppercounter,:);
% % % % % %     these_thetas = thetas(i:uppercounter);
% % % % %     this_sorted = sorted(uppercounter-1:uppercounter,:);
% % % % %     if (0 <= this_sorted(1,1) && this_sorted(1,1) <= pi/2) && (3*pi/2 <= this_sorted(2,1) && this_sorted(2,1) <= 2*pi)
% % % % %         % in such a case sorted(2,:) actually belongs to the first graphene
% % % % %         % group as we are traversing counterclockwise.
% % % % %         first_graphene_disks(i,:) = this_sorted(2,2:3);
% % % % %         second_graphene_disks(i,:) = this_sorted(1,2:3);
% % % % %     else % standard ordering: more advanced angle is the second graphene
% % % % %         first_graphene_disks(i,:) = this_sorted(1,2:3);
% % % % %         second_graphene_disks(i,:) = this_sorted(2,2:3);
% % % % %     end
% % % % %     uppercounter = uppercounter + 2;
% % % % % end



end

