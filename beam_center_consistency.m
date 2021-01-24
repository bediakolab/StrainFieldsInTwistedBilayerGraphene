% beam_center_consistency.m
%
% Script for ascertaining how consistently our center detection methods are
% working over multiple diffraction patterns collected.
%
% Nathanael Kazmierczak, 08/11/2019
% UC Berkeley, Bediako Lab.

%% Load in the data
start = [1,1,1,1];
count = [512,512,70,1]; % 100 diffraction patterns from the first row.
data = h5read('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
firstDP = data(:,:,1);

%% Make relevant masks
% Central power law mask
use_teeny_mask = 0;
teeny_mask_radius_factor = 0;
[beamstop_masks,Irng,Jrng,first_maskeddata] = makeBeamstopCenterMasks(firstDP,use_teeny_mask,teeny_mask_radius_factor);
% Outer graphene masks
disp('Please specify hexagonal window masks for outer graphene.');
hexagonal_window_masks_graphene = getHexagonalWindow(firstDP);

%% Iterate over the 100 diffraction patterns and collect origins.
numDP = size(data,3);
%numDP = 4;
poisson_origin_storage = zeros(numDP,2);
graphene_origin_storage = zeros(numDP,2);
graphene_disk_storage = zeros(0,3);
kernel_radius = 5;
sigma = 0;
min_disk_distance = 8;
num_to_find_per_window = 2;
make_figures = 0;

for i = 1:numDP
    thisDP = data(:,:,i);
    
    % Obtain the global Poisson MLE coordinates
    % Bootfun versions are simply used here because the syntax is compact,
    % and they should be robust after having fixed the numerical stability
    % issues.
    this_maskeddata = applyBeamstopCenterMasks(thisDP,beamstop_masks);
    linear_filtered_data = makeFilteredLinearDataFromImageMatrix(this_maskeddata);
    PoissonMLEorigin_coords = PoissonPowerlawBootfun(linear_filtered_data);
    PoissonMLEorigin_global_coords = convertMaskedCoordsToOriginalCoords(PoissonMLEorigin_coords,Irng,Jrng);
    poisson_origin_storage(i,:) = PoissonMLEorigin_global_coords;
    
    % Obtain the (global) graphene COM coordinates
    [disk_locations, disk_intensities, com_coordinates_graphene] = braggDiskRegistrationRecipie_1(...
        thisDP,hexagonal_window_masks_graphene,kernel_radius,sigma,min_disk_distance,num_to_find_per_window,make_figures);
    graphene_origin_storage(i,:) = fliplr(com_coordinates_graphene);
    disk_locations_id = [repmat(i,size(disk_locations,1),1), disk_locations];
    graphene_disk_storage = vertcat(graphene_disk_storage,disk_locations_id);
end

% useful for seeing how many disks detected for each DP
% size(graphene_disk_storage(graphene_disk_storage(:,1) == 4, :),1)

%% Make graphs
figure
plot(poisson_origin_storage(:,1),poisson_origin_storage(:,2),'-ok',graphene_origin_storage(:,1),graphene_origin_storage(:,2),'-or');
legend('Poisson MLE power law fit origin','Graphene COM origin');
title 'Comparison of center finding methods on 70 DPs, dataset 3';
xlabel('x (pixels)');
ylabel('y (pixels)');
