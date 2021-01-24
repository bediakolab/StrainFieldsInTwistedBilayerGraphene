function [graphene1_latticevectors,graphene2_latticevectors] = calculateLatticeVectors(first_graphene_disks,second_graphene_disks,hk_graphene1,hk_graphene2)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

% Add the ones so the origin can move
design_graphene1 = [ones(size(hk_graphene1,1),1), hk_graphene1];
design_graphene2 = [ones(size(hk_graphene2,1),1), hk_graphene2];

% Solve for coordinates of the first graphene 
graphene1_latticevectors = design_graphene1\first_graphene_disks;
% second graphene
graphene2_latticevectors = design_graphene2\second_graphene_disks;

end

