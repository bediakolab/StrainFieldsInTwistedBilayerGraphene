function [graphene1_latticevectors,graphene2_latticevectors] = computeFullLatticeVectors(disk_locations,disk_intensities,com_coordinates_graphene,...
    graphene1_lattice_for_indexing,graphene2_lattice_for_indexing)
% Assumption that initial guess has already been drawn for the peaks.
[first_graphene_disks,second_graphene_disks] = sortGraphenePeaks(disk_locations,disk_intensities,com_coordinates_graphene);
% Catch in case too many nans
% % % if isempty(first_graphene_disks) || isempty(second_graphene_disks)
% % %     graphene1_latticevectors = nan;
% % %     graphene2_latticevectors = nan;
% % %     return
% % % end
%% Index all lattice peaks
[hk_graphene1,hk_graphene2] = indexGraphenePeaks(graphene1_lattice_for_indexing,graphene2_lattice_for_indexing,...
    first_graphene_disks,second_graphene_disks,com_coordinates_graphene);

%% Calculate lattice vectors
[graphene1_latticevectors,graphene2_latticevectors] = calculateLatticeVectors(...
    first_graphene_disks,second_graphene_disks,hk_graphene1,hk_graphene2);

end

