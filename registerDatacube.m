function [ graphene1_lattice_storage, graphene2_lattice_storage, hk_graphene1_storage, hk_graphene2_storage, disk_all_storage ] = registerDatacube( data )
% Wrapper function for getting disk registration and index positions.
% Not much control, just to clean it out of the way when we are trying to
% test something else. In particular used in conjunction with
% average_disk_shape_testing.m, where we don't care about doing the strain
% mapping calculation. (But it's free, so let's spit it out anyway.)
%% Find Disks
sizevec = size(data);
numdim = numel(sizevec);

if numdim == 3
    firstDP = data(:,:,1);
    idx4num = 1;
    numDPs = sizevec(3);
elseif numdim == 4
    firstDP = data(:,:,1,1);
    idx4num = sizevec(4);
    numDPs = sizevec(3)*sizevec(4);
end 

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
disk_all_storage = cell(numDPs,1);
hk_graphene1_storage = cell(numDPs,1);
hk_graphene2_storage = cell(numDPs,1);
make_figures = 0;
counter = 1;
subpixel = 'False';
for i = 1:size(data,3)
    for j = 1:idx4num
        if numdim == 3
            thisDP = data(:,:,i);
        elseif numdim == 4
            thisDP = data(:,:,i,j);
        end
        
        [disk_locations, disk_intensities, com_coordinates_graphene] = braggDiskRegistrationRecipie_1(...
            thisDP,hexagonal_window_masks,kernel_radius,sigma,...
            min_disk_distance,num_to_find_per_window,make_figures,lower_threshold,upper_threshold,subpixel);
        %             com_coordinates_graphene = fliplr(com_coordinates_graphene);  % so that they are xy coordinates.
        %             disk_locations = fliplr(disk_locations);  % so that they are xy coordinates.
        %
        
        % This is redundant as it is also performed in
        % computeFullLatticeVectors. But we want to extract this
        % information for averaging the disks, so compute here so we can
        % save.
        [first_graphene_disks,second_graphene_disks] = sortGraphenePeaks(disk_locations,disk_intensities,com_coordinates_graphene);
        [hk_graphene1,hk_graphene2] = indexGraphenePeaks(graphene1_lattice_for_indexing,graphene2_lattice_for_indexing,...
            first_graphene_disks,second_graphene_disks,com_coordinates_graphene);
        
        [graphene1_latticevectors,graphene2_latticevectors] = computeFullLatticeVectors(...
            disk_locations,disk_intensities,com_coordinates_graphene,graphene1_lattice_for_indexing,graphene2_lattice_for_indexing);
        % % %             if any(isnan(graphene1_latticevectors)) | any(isnan(graphene2_latticevectors))
        % % %                 graphene1_lattice_storage(:,:,i,j) = nan(3,2);
        % % %                 graphene2_lattice_storage(:,:,i,j) = nan(3,2);
        % % %                 counter = counter + 1;
        % % %                 disp(counter);
        % % %                 continue
        % % %             end
        disk_all_storage{counter} = disk_locations;
        hk_graphene1_storage{counter} = hk_graphene1;
        hk_graphene2_storage{counter} = hk_graphene2;
        
        graphene1_lattice_storage(:,:,i,j) = graphene1_latticevectors;
        graphene2_lattice_storage(:,:,i,j) = graphene2_latticevectors;
        counter = counter + 1;
        disp(counter);
    end
end



end

