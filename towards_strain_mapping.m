% towards_strain_mapping.m

%
% As of 08/27/2019, this file is runnable off of NPK laptop because the
% appropriate .h5 files have been downloaded to his thumb drive.
% 
% towards_strain_mapping_fulldata.m I believe does pretty much the same
% thing, but loads the data in piece by piece and then deletes so as not to
% overwhelm the Matlab RAM.
%

addpath('/Volumes/SANDISK/BediakoLabLexar');
addpath('/Volumes/SANDISK/BediakoLabLexar/4DSTEMdata_updating');

%% Load
start = [1,1,1,1];
count = [512,512,1,1];  % the full dimensions of this cropped piece of data.
data = h5read('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
% data = ReadDMFile('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um.dm4');
firstDP = data; %data(:,:,1);


%% Find Disks

arcangle = pi/12;
annular_increment = 20;
hexagonal_window_masks = getHexagonalWindow(firstDP,arcangle,annular_increment);
kernel_radius = 5;
sigma = 0;
min_disk_distance = 8;
num_to_find_per_window = 2;
make_figures = 1;
lower_threshold = 0.5;
upper_threshold = 2;
subpixel = 'False';
[disk_locations, disk_intensities, com_coordinates_graphene] = braggDiskRegistrationRecipie_1(...
    firstDP,hexagonal_window_masks,kernel_radius,sigma,min_disk_distance,num_to_find_per_window,make_figures,lower_threshold,upper_threshold,subpixel);
% Fixed because this was stupid
% com_coordinates_graphene = fliplr(com_coordinates_graphene);  % so that they are xy coordinates.
% disk_locations = fliplr(disk_locations);  % so that they are xy coordinates.

%% Sort into different sets for each graphene
[first_graphene_disks,second_graphene_disks] = sortGraphenePeaks(disk_locations,disk_intensities,com_coordinates_graphene);
%Verify
figure; image(uint8(firstDP),'CDataMapping','scaled'); colormap pink; hold on;
scatter(first_graphene_disks(:,1),first_graphene_disks(:,2),'r','filled');
scatter(second_graphene_disks(:,1),second_graphene_disks(:,2),'c','filled');

%% Draw initial guess for best lattice vectors
[graphene1_lattice_for_indexing,graphene2_lattice_for_indexing] = guessInitialLatticeVectorsForIndexing(...
    firstDP,first_graphene_disks,second_graphene_disks,com_coordinates_graphene);

%% Index all lattice peaks
[hk_graphene1,hk_graphene2] = indexGraphenePeaks(graphene1_lattice_for_indexing,graphene2_lattice_for_indexing,...
    first_graphene_disks,second_graphene_disks,com_coordinates_graphene);

%% Calculate lattice vectors
[graphene1_latticevectors,graphene2_latticevectors] = calculateLatticeVectors(...
    first_graphene_disks,second_graphene_disks,hk_graphene1,hk_graphene2);

% This has completed the major mechanics. We will now need to package up
% our functions into a larger one in order to compute strain, because we
% need reference vectors, which we don't have.

%% Lattice vectors for multiple DPs
% lets go for storing these in a 3-way array, with vertical slices
% corresponding to a new DP.
graphene1_lattice_storage = zeros(3,2,size(data,3),size(data,4));
graphene2_lattice_storage = zeros(3,2,size(data,3),size(data,4));
make_figures = 0;
counter = 1;
for i = 1:size(data,3)
    for j = 1:size(data,4)
        [disk_locations, disk_intensities, com_coordinates_graphene] = braggDiskRegistrationRecipie_1(...
            data(:,:,i,j),hexagonal_window_masks,kernel_radius,sigma,min_disk_distance,num_to_find_per_window,make_figures,lower_threshold,upper_threshold);
        com_coordinates_graphene = fliplr(com_coordinates_graphene);  % so that they are xy coordinates.
        disk_locations = fliplr(disk_locations);  % so that they are xy coordinates.
        
        [graphene1_latticevectors,graphene2_latticevectors] = computeFullLatticeVectors(...
            disk_locations,disk_intensities,com_coordinates_graphene,graphene1_lattice_for_indexing,graphene2_lattice_for_indexing);
        if any(isnan(graphene1_latticevectors)) | any(isnan(graphene2_latticevectors))
            graphene1_lattice_storage(:,:,i,j) = nan(3,2);
            graphene2_lattice_storage(:,:,i,j) = nan(3,2);
            counter = counter + 1;
            disp(counter);
            continue
        end
        
        graphene1_lattice_storage(:,:,i,j) = graphene1_latticevectors;
        graphene2_lattice_storage(:,:,i,j) = graphene2_latticevectors;
        counter = counter + 1;
        disp(counter);
    end
end

%% Average down to get some references
% If necessary, average over the fourth dimension too first
graphene1_unstrained_lattice = mean(mean(graphene1_lattice_storage,4),3,'omitnan');
graphene2_unstrained_lattice = mean(mean(graphene2_lattice_storage,4),3,'omitnan');
graphene1_unstrained_lattice = graphene1_unstrained_lattice(2:3,:);
graphene2_unstrained_lattice = graphene2_unstrained_lattice(2:3,:);

%% Compute strain
[graphene1strain] = computeStrain(graphene1_unstrained_lattice,graphene1_lattice_storage);
[graphene2strain] = computeStrain(graphene2_unstrained_lattice,graphene2_lattice_storage);

%% Make some plots
plotStrain(graphene1strain,'Graphene1BeforeSubpixel')
plotStrain(graphene2strain,'Graphene2BeforeSubpixel')




