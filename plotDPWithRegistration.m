function plotDPWithRegistration(filename,DP_index_3,DP_index_4,hexagonal_window_masks)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

kernel_radius = 5;
sigma = 0;
min_disk_distance = 8;
num_to_find_per_window = 2;
make_figures = 1;
lower_threshold = 0.5;
upper_threshold = 2;

start = [1,1,DP_index_3,DP_index_4];
count = [512,512,1,1];  % the full dimensions of this cropped piece of data.
DP = h5read(filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);


[disk_locations, disk_intensities, com_coordinates_graphene] = braggDiskRegistrationRecipie_1(...
      DP,hexagonal_window_masks,kernel_radius,sigma,min_disk_distance,num_to_find_per_window,make_figures,lower_threshold,upper_threshold);
title(strcat(['Coordinates: (',num2str(DP_index_3),',',num2str(DP_index_4),')']));          

end

