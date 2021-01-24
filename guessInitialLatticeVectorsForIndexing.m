function [graphene1_lattice_for_indexing,graphene2_lattice_for_indexing] = guessInitialLatticeVectorsForIndexing(plotdata,graphene_disks_1,graphene_disks_2,origin)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

figure; image(uint8(plotdata),'CDataMapping','scaled'); colormap pink; hold on;
scatter(graphene_disks_1(:,1),graphene_disks_1(:,2),'r','filled');
scatter(graphene_disks_2(:,1),graphene_disks_2(:,2),'c','filled');
scatter(origin(1),origin(2),'g','filled');
legend('Graphene disks 1','Graphene disks 2','Origin');

disp('For graphene 1, please guess the initial lattice vectors.');
disp('Click on two spots on the graph that will be used to index all active peaks.');
g1_clicks = ginput(2);

disp('For graphene 2, please guess the initial lattice vectors.');
disp('Click on two spots on the graph that will be used to index all active peaks.');
g2_clicks = ginput(2);

% vectors will be represented from the origin
graphene1_lattice_for_indexing = g1_clicks - [origin;origin];  % origin should be in xy coords
graphene2_lattice_for_indexing = g2_clicks - [origin;origin];  % origin should be in xy coords

end

