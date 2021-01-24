function [ roicut, Irangecut, Jrangecut ] = cutROI( DP, disk_locations, hk_disks, index_choice, window_width )
% Function called by getAveragedDiskShape.m
%
% Goal is to cut an roi around the specified hk graphene peaks. 
%
% Disk locations must be submitted in I,J format, not x,y.

pos1_match = hk_disks(:,1) == index_choice(1);
pos2_match = hk_disks(:,2) == index_choice(2);
userow = pos1_match & pos2_match;   % Logical indexing.
disk_location = disk_locations(userow,:);

[ roicut, Irangecut, Jrangecut ] = cutROI_fromPosition( DP, disk_location, window_width );

end

