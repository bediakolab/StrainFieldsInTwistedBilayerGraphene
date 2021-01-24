% towards_strain_mapping_fulldata.m
%
% 08142019 revised so that only five window masks are drawn and we don't
% throw a data point out even if one of the vectors is missing, but rather
% just do the regression without it present, throwing out its Miller index
% additionally. (In reality the Miller index will never even exist, so it
% will all be fine.)

addpath('C:\Users\Bediako Lab\Desktop\4DSTEM data');
addpath('C:\Users\Bediako Lab\Desktop\4DSTEM data\Maddiedata07292019');
addpath('C:\Users\Bediako Lab\Desktop\PythonWorkspace\py4DSTEM\py4DSTEM\process\utils');

addpath('/Volumes/SANDISK/BediakoLabLexar');
addpath('/Volumes/SANDISK/BediakoLabLexar/4DSTEMdata_updating');


%% Load
start = [1,1,1,1];
% start = [1,1,30,50];  % Trying to test one of the ones that was failing.
count = [512,512,1,1];  % the full dimensions of this cropped piece of data.
data = h5read('3__70x70_ss4nm_0p05s_spot9)alpha=0p55_bin=4_CL=130_80kV_C2=40um.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);

% data = ReadDMFile('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um.dm4');
firstDP = data; %data(:,:,1);


%% Find Disks
arcangle = pi/12;
annular_increment = 20;
% hexagonal_window_masks = getHexagonalWindow(firstDP,arcangle,annular_increment);
numberofsites = 5;
[hexagonal_window_masks] = getFreehandMasks(firstDP,numberofsites);
kernel_radius = 5;
sigma = 0;
min_disk_distance = 8;
num_to_find_per_window = 2;
make_figures = 1;
lower_threshold = 0.5;
upper_threshold = 2;
subpixel = 'True';
[disk_locations, disk_intensities, com_coordinates_graphene] = braggDiskRegistrationRecipie_1(...
    firstDP,hexagonal_window_masks,kernel_radius,sigma,min_disk_distance,num_to_find_per_window,make_figures,lower_threshold,upper_threshold,subpixel);
% Fixed this in the function call itself: these are now xycoordinates
% naturally.
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

%% 08142019: loop over the datacube, loading in 70x10 chunks at a time.
division_size = 10;
graphene1_lattice_storage_outer = zeros(3,2,70,70);
graphene2_lattice_storage_outer = zeros(3,2,70,70);
for k = 1:7
    start = [1,1,1,(k-1)*10+1];
    count = [512,512,70,10];  % chunks are in the fourth dimension
    filename = '3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um.h5';
    this_data = h5read(filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
    
    %% Lattice vectors for multiple DPs
    % lets go for storing these in a 3-way array, with vertical slices
    % corresponding to a new DP.
    graphene1_lattice_storage = zeros(3,2,size(this_data,3),size(this_data,4));
    graphene2_lattice_storage = zeros(3,2,size(this_data,3),size(this_data,4));
    make_figures = 0;
    counter = 1;
    for i = 1:size(this_data,3)
        for j = 1:size(this_data,4)
            [disk_locations, disk_intensities, com_coordinates_graphene] = braggDiskRegistrationRecipie_1(...
                this_data(:,:,i,j),hexagonal_window_masks,kernel_radius,sigma,...
                min_disk_distance,num_to_find_per_window,make_figures,lower_threshold,upper_threshold,subpixel);
%             com_coordinates_graphene = fliplr(com_coordinates_graphene);  % so that they are xy coordinates.
%             disk_locations = fliplr(disk_locations);  % so that they are xy coordinates.
%             
            [graphene1_latticevectors,graphene2_latticevectors] = computeFullLatticeVectors(...
                disk_locations,disk_intensities,com_coordinates_graphene,graphene1_lattice_for_indexing,graphene2_lattice_for_indexing);
% % %             if any(isnan(graphene1_latticevectors)) | any(isnan(graphene2_latticevectors))
% % %                 graphene1_lattice_storage(:,:,i,j) = nan(3,2);
% % %                 graphene2_lattice_storage(:,:,i,j) = nan(3,2);
% % %                 counter = counter + 1;
% % %                 disp(counter);
% % %                 continue
% % %             end
            
            graphene1_lattice_storage(:,:,i,j) = graphene1_latticevectors;
            graphene2_lattice_storage(:,:,i,j) = graphene2_latticevectors;
            counter = counter + 1;
            disp(counter);
        end
    end
    
    % Now add these into a larger storage and delete the given piece of
    % data.
    graphene1_lattice_storage_outer(:,:,:,((k-1)*10+1):k*10) = graphene1_lattice_storage;
    graphene2_lattice_storage_outer(:,:,:,((k-1)*10+1):k*10) = graphene2_lattice_storage;
    clear this_data
end
%% Average down to get some references
% If necessary, average over the fourth dimension too first
% graphene1_unstrained_lattice = mean(mean(graphene1_lattice_storage_outer,4,'omitnan'),3,'omitnan');
% graphene2_unstrained_lattice = mean(mean(graphene2_lattice_storage_outer,4,'omitnan'),3,'omitnan');
graphene1_unstrained_lattice = median(median(graphene1_lattice_storage_outer,4),3);
graphene2_unstrained_lattice = median(median(graphene2_lattice_storage_outer,4),3);
graphene1_unstrained_lattice = graphene1_unstrained_lattice(2:3,:);
graphene2_unstrained_lattice = graphene2_unstrained_lattice(2:3,:);

%% Compute strain
[graphene1strain] = computeStrain(graphene1_unstrained_lattice,graphene1_lattice_storage_outer);
[graphene2strain] = computeStrain(graphene2_unstrained_lattice,graphene2_lattice_storage_outer);

%% Make some plots
figure;
hold on;
plotStrain(graphene1strain,'Graphene1BeforeSubpixel_nf',1,0)
plotStrain(graphene2strain,'Graphene2BeforeSubpixel_nf',2,0)
figure;
hold on;
plotStrain(graphene1strain,'Graphene1BeforeSubpixel_f',1,1)
plotStrain(graphene2strain,'Graphene2BeforeSubpixel_f',2,1)

