function [hk_graphene1,hk_graphene2] = indexGraphenePeaks(graphene1_lattice_for_indexing,graphene2_lattice_for_indexing,...
    first_graphene_disks,second_graphene_disks,origin)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% Center the whole lot of disks (lattice vectors already centered).
first_graphene_disks_centered = first_graphene_disks - repmat(origin,size(first_graphene_disks,1),1);
second_graphene_disks_centered = second_graphene_disks - repmat(origin,size(second_graphene_disks,1),1);

hk_graphene1 = (graphene1_lattice_for_indexing'\first_graphene_disks_centered')'; 
hk_graphene2 = (graphene2_lattice_for_indexing'\second_graphene_disks_centered')'; 

hk_graphene1 = round(hk_graphene1);
hk_graphene2 = round(hk_graphene2);

end

